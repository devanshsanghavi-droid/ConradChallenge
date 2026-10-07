"""
simulate.py — ground truth + everything the surrogate is allowed to see.

Two models are built from the same EPANET .inp:

  * TRUTH  : hidden per-pipe wall-decay coefficients (old/rough pipes decay faster),
             perturbed demands, perturbed source dose.  The surrogate never sees these.
  * NOMINAL: the operator's own (imperfect) hydraulic model.  From it we take everything
             EPANET gives us for free without a single chlorine measurement:
               - water age (AGE quality mode)
               - hydraulics: pressure, pipe flow and velocity, flow direction
               - the pipe table: length, diameter, Hazen-Williams roughness (material/age proxy)
               - topology: hydraulic distances, paths from sources, tanks
             and a chlorine simulator whose decay parameters are UNKNOWN and get calibrated
             from grab samples (see simgp.py).

Everything runs on WNTR (EPA/Sandia), which bundles the EPANET 2.2 engine.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd
import wntr

from .chemistry import (ER_TRUTH_RANGE_K, THETA_W_TRUTH_RANGE, TREF_C, TRUTH_SEED_OFFSET,
                        Chemistry, apply_water_temperature, arrhenius, remove_epanet_files, toc_ratio,
                        use_mg_per_litre)

DAY = 86400
DURATION_DAYS = 7                  # every run is 7 days; the last day is scored
QUALITY_STEP_S = 300
LIB = os.path.join(os.path.dirname(wntr.__file__), "library", "networks")


def net_path(name: str = "Net3") -> str:
    """Path to a network bundled with WNTR (Net1/2/3/6, ky4, ky10) or a user .inp."""
    return name if os.path.exists(name) else os.path.join(LIB, f"{name}.inp")


def load(name: str, duration_days: int = DURATION_DAYS) -> wntr.network.WaterNetworkModel:
    wn = wntr.network.WaterNetworkModel(net_path(name))
    wn.options.time.duration = duration_days * DAY
    wn.options.time.hydraulic_timestep = 3600
    wn.options.time.report_timestep = 3600
    wn.options.time.quality_timestep = QUALITY_STEP_S
    return wn


def roughness_factor(c: float, gamma: float = 1.0) -> float:
    """Hypothesis encoded in both truth and calibrated simulator: rougher (older) pipe -> faster
    wall decay.  C=130 -> 1x, C=110 -> 1.6x, C=199 -> 0.2x when gamma=1; gamma=0 switches it off."""
    return 2.0 ** (-gamma * (c - 130.0) / 30.0)


def source_nodes(wn) -> list[str]:
    """Where water (and chlorine) enters: reservoirs; else junctions with a negative demand (a pumped
    inflow modelled as a negative demand, as in Net2); else the tanks.  A CONCEN source only acts on
    external inflow at a node, so on a tank-fed file without a reservoir the source must sit at the
    inflow junction — at the tank itself it does nothing."""
    if wn.reservoir_name_list:
        return list(wn.reservoir_name_list)
    inflow = [j for j, n in wn.junctions() if any(ts.base_value < 0 for ts in n.demand_timeseries_list)]
    return inflow or list(wn.tank_name_list)


def hydraulic_graph(wn) -> nx.Graph:
    g = nx.Graph()
    for lname, link in wn.links():
        length = getattr(link, "length", None) or 1.0  # pumps/valves: 1 m
        a, b = link.start_node_name, link.end_node_name
        if not g.has_edge(a, b) or g[a][b]["weight"] > length:
            g.add_edge(a, b, weight=length, link=lname)
    return g


@dataclass
class Scenario:
    wn_name: str
    seed: int
    sample_hour: int
    junctions: list[str]
    truth_snapshot: pd.Series | None        # mg/L at sampling hour, last day (None for a real .inp)
    truth_daily_min: pd.Series | None       # mg/L daily minimum, last day
    truth_by_hour: pd.DataFrame | None      # hour x junction, last day
    age_by_hour_h: pd.DataFrame      # hour x junction, NOMINAL model
    hyd_dist: pd.DataFrame           # length-weighted shortest path (m), junction x junction
    graph: nx.Graph
    pipes: pd.DataFrame              # per pipe: start, end, length, diameter, roughness, flow, velocity
    node_hyd: pd.DataFrame           # per junction: pressure stats from nominal hydraulics
    coords: dict
    wn: wntr.network.WaterNetworkModel  # nominal model
    source_dose: float
    structural: dict | None = None      # what structural_noise did to the truth (None if off)
    truth_age_by_hour_h: pd.DataFrame | None = None   # hour x junction, water age on the TRUTH (truth_age=True)
    chem: dict | None = None            # the chemistry of the truth and its hidden draws (chem=... only)

    @property
    def age_snapshot_h(self):
        return self.age_by_hour_h.loc[self.sample_hour]

    @property
    def age_daily_mean_h(self):
        return self.age_by_hour_h.mean()


def _last_day(df: pd.DataFrame, cols) -> pd.DataFrame:
    out = df.loc[df.index >= 6 * DAY, cols].copy()
    out.index = ((out.index - 6 * DAY) // 3600).astype(int)
    return out.loc[out.index < 24]


def simulate_nominal_chlorine(name: str, kb_per_day: float, kw_m_per_day: float, gamma: float,
                              source_dose: float = 1.2, demand_mult: float = 1.0,
                              rough_mult: float = 1.0, file_prefix: str = "temp",
                              kb_scale: float = 1.0, kw_scale: float = 1.0,
                              temp_C: float | None = None) -> pd.DataFrame:
    """Chlorine on the operator's NOMINAL model for a candidate (kb, kw, gamma).  hour x junction.

    demand_mult / rough_mult are the hydraulic-mismatch axes of the grid: the operator's demands and
    Hazen-Williams C are scaled globally.  The wall-decay hypothesis stays on the operator's C table
    (as in the truth: roughness noise is a hydraulic error, not a chemistry one).

    kb_scale / kw_scale multiply the bulk and wall rates (e.g. an Arrhenius factor for this month's water,
    chemistry.Chemistry.kb_scale); temp_C multiplies the file's own VISCOSITY and DIFFUSIVITY by their ratios
    at that temperature (chemistry.apply_water_temperature).  The defaults write an identical .inp."""
    wn = load(name)
    wn.options.quality.parameter = "CHEMICAL"
    wn.options.reaction.bulk_coeff = -kb_per_day * kb_scale / DAY
    wn.options.reaction.wall_coeff = -kw_m_per_day * kw_scale / DAY
    for _, pipe in wn.pipes():
        pipe.wall_coeff = -kw_m_per_day * kw_scale * roughness_factor(pipe.roughness, gamma) / DAY
        pipe.roughness = pipe.roughness * rough_mult
    apply_water_temperature(wn, temp_C)
    if demand_mult != 1.0:
        for _, j in wn.junctions():
            for ts in j.demand_timeseries_list:
                ts.base_value = ts.base_value * demand_mult
    for res in source_nodes(wn):
        wn.add_source(f"src_{res}", res, "CONCEN", source_dose)
    q = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=file_prefix).node["quality"]
    return _last_day(q, wn.junction_name_list).clip(lower=0.0)


def closable_pipes(wn) -> list[str]:
    """Pipes whose closure keeps the network connected: not a bridge of the link graph, or with a
    parallel link between the same two nodes."""
    g = hydraulic_graph(wn)
    bridges = {frozenset(e) for e in nx.bridges(g)}
    ends = {}
    for lname, link in wn.links():
        ends.setdefault(frozenset((link.start_node_name, link.end_node_name)), []).append(lname)
    out = []
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        key = frozenset((p.start_node_name, p.end_node_name))
        if key not in bridges or len(ends[key]) > 1:
            out.append(pn)
    return out


def apply_structural_noise(wn, rng, mode: str = "spec") -> dict:
    """Task 3: errors in the operator's file that no grid axis can represent.  With probability 0.5 close
    one random non-bridge pipe (connectivity re-checked with networkx); always perturb one random tank.
    Applied to the TRUTH only.

    mode="spec"       : the tank's initial level x0.7 (clipped to its minimum level).  NOTE: a 7-day warm-up
                        forgets an initial level — on Net3 this changes the scored day by < 0.05 mg/L.
    mode="persistent" : the tank's diameter x0.7 (it holds half the water the operator's file says), a
                        geometric error that persists: different turnover, different pump cycling."""
    info = {"mode": mode, "closed_pipe": None, "tank": None, "tank_level_m": None, "tank_diameter_m": None}
    if rng.uniform() < 0.5:
        cands = closable_pipes(wn)
        rng.shuffle(cands)
        for pn in cands:
            g = hydraulic_graph(wn)
            p = wn.get_link(pn)
            g2 = g.copy()
            # remove only if this is the sole link between the two nodes
            if sum(1 for _, l in wn.links() if {l.start_node_name, l.end_node_name} == {p.start_node_name, p.end_node_name}) == 1:
                g2.remove_edge(p.start_node_name, p.end_node_name)
            if nx.is_connected(g2):
                p.initial_status = wntr.network.LinkStatus.Closed
                info["closed_pipe"] = pn
                break
    if wn.tank_name_list:
        tn = str(rng.choice(wn.tank_name_list))
        t = wn.get_node(tn)
        info["tank"] = tn
        if mode == "persistent":
            info["tank_diameter_m"] = (round(float(t.diameter), 2), round(float(0.7 * t.diameter), 2))
            t.diameter = 0.7 * t.diameter
        else:
            new_level = max(t.min_level + 0.05 * (t.max_level - t.min_level), 0.7 * t.init_level)
            info["tank_level_m"] = (round(float(t.init_level), 2), round(float(new_level), 2))
            t.init_level = new_level
    return info


def build_scenario(name: str = "Net3", seed: int = 0, sample_hour: int = 14,
                   source_dose: float = 1.2, kb_per_day: float = 0.40,
                   kw_m_per_day: float = 0.70, structural_noise: bool | str = False,
                   month_seed: int | None = None, chem: Chemistry | None = None,
                   truth_age: bool = False) -> Scenario:
    """structural_noise: False, True (= "spec") or "persistent"; see apply_structural_noise.
    month_seed: if given, the pipe-level truth (per-pipe wall decay, roughness) comes from `seed` and the
    operating truth (bulk decay, demand, dose) from `month_seed`: the same network in a different month.
    chem: None is the committed truth, draw for draw.  A chemistry.Chemistry gives the truth that chemistry
    (see _chem_truth); it consumes the same rng / rng_m draws in the same order, and every new draw comes from
    its own generator, default_rng(20_000 + seed) for free chlorine or default_rng(40_000 + seed) for chloramine
    (a chloramine truth is refused until task 12 builds its decay physics).
    truth_age: also run EPANET AGE on the truth's perturbed network after all draws (no new draws)."""
    rng = np.random.default_rng(seed)
    rng_m = np.random.default_rng(month_seed) if month_seed is not None else rng

    # ---------------- TRUTH (hidden) ----------------
    wn = load(name)
    wn.options.quality.parameter = "CHEMICAL"
    if truth_age:   # what the AGE run needs back after the chemistry recipe has changed them
        age_state = (wn.options.quality.tolerance, {n: node.initial_quality for n, node in wn.nodes()})
    structural = None
    if structural_noise:
        mode = structural_noise if isinstance(structural_noise, str) else "spec"
        structural = apply_structural_noise(wn, np.random.default_rng(10_000 + seed), mode)
    chem_info = None
    if chem is None:
        wn.options.reaction.bulk_coeff = -kb_per_day * rng_m.uniform(0.8, 1.2) / DAY
        wn.options.reaction.wall_coeff = -kw_m_per_day / DAY
        for _, pipe in wn.pipes():
            pipe.wall_coeff = -kw_m_per_day * roughness_factor(pipe.roughness, 1.0) * np.exp(rng.normal(0, 0.4)) / DAY
            pipe.roughness = pipe.roughness * np.exp(rng.normal(0, 0.10))   # hydraulic model mismatch
        global_mult = rng_m.uniform(0.85, 1.15)
        for _, j in wn.junctions():
            for ts in j.demand_timeseries_list:
                ts.base_value = ts.base_value * global_mult * np.exp(rng_m.normal(0.0, 0.15))
        for res in source_nodes(wn):
            wn.add_source(f"src_{res}", res, "CONCEN", source_dose * rng_m.uniform(0.9, 1.1))
        q = wntr.sim.EpanetSimulator(wn).run_sim().node["quality"]
    else:
        q, chem_info = _chem_truth(wn, seed, rng, rng_m, chem, source_dose, kb_per_day, kw_m_per_day)
    junctions = wn.junction_name_list
    truth_by_hour = _last_day(q, junctions).clip(lower=0.0)
    truth_age_h = _truth_age(wn, *age_state) if truth_age else None

    sc = nominal_scenario(name, sample_hour, source_dose)
    sc.seed, sc.structural = seed, structural
    sc.truth_snapshot, sc.truth_daily_min, sc.truth_by_hour = truth_by_hour.loc[sample_hour], truth_by_hour.min(), truth_by_hour
    sc.truth_age_by_hour_h, sc.chem = truth_age_h, chem_info
    return sc


def _run_quality(wn) -> pd.DataFrame:
    """One EPANET run in its own temp directory (no files left in the cwd); node quality as WNTR returns it."""
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as tmp:
        prefix = os.path.join(tmp, "run")
        try:
            return wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=prefix).node["quality"]
        finally:
            remove_epanet_files(prefix)


def hidden_chem_draws(disinfectant: str, seed: int) -> dict:
    """The truth's own chemistry draws for one seed, from default_rng(TRUTH_SEED_OFFSET[disinfectant] + seed),
    always in this fixed order whether or not a condition uses them, so later additions never shift them:
      E_true  ~ U(4660, 12104) K   bulk Arrhenius E/R of this network's water (the compiled span quoted by
                                   Cejas, Diaz & Gonzalez 2026)
      theta_w ~ U(1.00, 1.07)      wall temperature factor theta_w^(T - 20) (ASSUMPTION: Lee et al. 2014 give
                                   only the direction, wall decay rising with temperature)"""
    rng_seed = TRUTH_SEED_OFFSET[disinfectant] + seed
    rc = np.random.default_rng(rng_seed)
    er_true = float(rc.uniform(*ER_TRUTH_RANGE_K))
    theta_w = float(rc.uniform(*THETA_W_TRUTH_RANGE))
    return {"rng_seed": rng_seed, "E_true_K": er_true, "theta_w": theta_w}


def _chem_truth(wn, seed, rng, rng_m, chem: Chemistry, source_dose, kb_per_day, kw_m_per_day):
    """The hidden truth under a chemistry condition.  The committed draws are consumed exactly as in the
    default branch of build_scenario (bulk factor, per-pipe wall and roughness, global and per-node demand,
    dose per source); the new ones come from hidden_chem_draws (their own generator).
    At T = 20 C (or None), TOC = TOC_ref (or None) and first order, every factor is exactly 1 and the truth
    equals build_scenario(chem=None) bit for bit (a saved check).  kinetics 'first_si' runs the same first
    order with the corrected unit recipe (chemistry.use_mg_per_litre); 'order2' and 'clark' arrive with the
    tasks that use them, and a chloramine truth is refused until task 12 builds its decay physics."""
    chem.require_built()
    if chem.kinetics not in ("first", "first_si"):
        raise NotImplementedError(f"truth kinetics {chem.kinetics!r} is not built yet (tasks 11 and 12)")
    hd = hidden_chem_draws(chem.disinfectant, seed)
    er_true, theta_w = hd["E_true_K"], hd["theta_w"]
    T = chem.temp_C
    fb = 1.0 if T is None else arrhenius(T, er_true)
    fw = 1.0 if T is None else theta_w ** (T - TREF_C)
    tr = toc_ratio(chem.toc_mgL)

    u = rng_m.uniform(0.8, 1.2)
    wn.options.reaction.bulk_coeff = -kb_per_day * u * fb * tr / DAY
    wn.options.reaction.wall_coeff = -kw_m_per_day * fw / DAY
    for _, pipe in wn.pipes():
        pipe.wall_coeff = -kw_m_per_day * fw * roughness_factor(pipe.roughness, 1.0) * np.exp(rng.normal(0, 0.4)) / DAY
        pipe.roughness = pipe.roughness * np.exp(rng.normal(0, 0.10))
    global_mult = rng_m.uniform(0.85, 1.15)
    for _, j in wn.junctions():
        for ts in j.demand_timeseries_list:
            ts.base_value = ts.base_value * global_mult * np.exp(rng_m.normal(0.0, 0.15))
    doses = {res: source_dose * rng_m.uniform(0.9, 1.1) for res in source_nodes(wn)}
    v_ratio, d_ratio = apply_water_temperature(wn, T)
    if chem.kinetics == "first":
        for res, dose in doses.items():
            wn.add_source(f"src_{res}", res, "CONCEN", dose)
        scale = 1.0
    else:
        scale = use_mg_per_litre(wn, doses)
    q = _run_quality(wn) * scale
    info = {"disinfectant": chem.disinfectant, "species": chem.species, "kinetics": chem.kinetics,
            "temp_C": T, "toc_mgL": chem.toc_mgL, "rng_seed": hd["rng_seed"],
            "E_true_K": er_true, "theta_w": theta_w, "bulk_month_factor": float(u),
            "bulk_temp_factor": fb, "bulk_toc_factor": tr, "wall_temp_factor": fw,
            "viscosity_ratio": v_ratio, "diffusivity_ratio": d_ratio,
            "kb_per_day_effective": kb_per_day * u * fb * tr,
            "source_doses_mgL": {k: float(v) for k, v in doses.items()}}
    return q, info


def _truth_age(wn, tolerance, initial_quality) -> pd.DataFrame:
    """Water age (hours, hour x junction, last day) on the truth's own perturbed network.  The file's
    tolerance and initial qualities (initial ages) are restored first, as the nominal AGE run uses them."""
    wn.options.quality.parameter = "AGE"
    wn.options.quality.tolerance = tolerance
    for n, q0 in initial_quality.items():
        wn.get_node(n).initial_quality = q0
    return _last_day(_run_quality(wn), wn.junction_name_list) / 3600.0


def nominal_scenario(name: str, sample_hour: int = 14, source_dose: float = 1.2) -> Scenario:
    """The operator's side only: what a real .inp gives us with no chlorine measurement.  Used by the
    app (no hidden truth) and by build_scenario (which adds one)."""
    wn_nom = load(name)
    junctions = wn_nom.junction_name_list
    wn_nom.options.quality.parameter = "AGE"
    res = wntr.sim.EpanetSimulator(wn_nom).run_sim()
    age_by_hour = _last_day(res.node["quality"], junctions) / 3600.0

    pressure = _last_day(res.node["pressure"], junctions)
    node_hyd = pd.DataFrame({"pressure_mean_m": pressure.mean(), "pressure_min_m": pressure.min()})

    pipe_names = wn_nom.pipe_name_list
    flow = _last_day(res.link["flowrate"], pipe_names)
    vel = _last_day(res.link["velocity"], pipe_names)
    rows = []
    for pn in pipe_names:
        p = wn_nom.get_link(pn)
        f = flow[pn]
        steady = f.abs().mean() > 0 and abs(f.mean()) / f.abs().mean() > 0.8
        rows.append({"pipe": pn, "start": p.start_node_name, "end": p.end_node_name,
                     "length_m": p.length, "diameter_m": p.diameter, "roughness": p.roughness,
                     "flow_mean_m3s": f.mean(), "abs_flow_mean_m3s": f.abs().mean(),
                     "velocity_mean_ms": vel[pn].mean(), "velocity_min_ms": vel[pn].min(),
                     # +1 if flow goes start->end most of the day, -1 if reversed, 0 if it sloshes
                     "direction": int(np.sign(f.mean())) if steady else 0})
    pipes = pd.DataFrame(rows).set_index("pipe")

    g = hydraulic_graph(wn_nom)
    d = dict(nx.all_pairs_dijkstra_path_length(g, weight="weight"))
    big = 10 * max(max(v.values()) for v in d.values())
    hyd = pd.DataFrame([[d[i].get(j, big) for j in junctions] for i in junctions],
                       index=junctions, columns=junctions, dtype=float)
    coords = {n: wn_nom.get_node(n).coordinates for n in g.nodes}

    return Scenario(name, 0, sample_hour, junctions, None, None, None,
                    age_by_hour, hyd, g, pipes, node_hyd, coords, wn_nom, source_dose, None)
