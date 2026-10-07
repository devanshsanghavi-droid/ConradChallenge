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

import functools
import math
import os
import tempfile
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd
import wntr

from .chemistry import (CHLORAMINE, ER_TRUTH_RANGE_K, THETA_W_TRUTH_RANGE, TOC_REF_MGL, TREF_C, TRUTH_SEED_OFFSET,
                        Chemistry, Warming, apply_water_temperature, arrhenius, remove_epanet_files, set_bulk_kinetics,
                        toc_ratio, use_mg_per_litre)

DAY = 86400
TRUTH_WALL_LAWS = ("theta_w", "arrhenius", "none")
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
    truth_loss_runs: dict | None = None  # the truth rerun with kw = 0, kb = 0 and no decay (truth_loss_split=True)
    truth_initial_share: pd.DataFrame | None = None  # hour x junction, share of the truth's water still the run's
                                                      # starting contents on the last day (truth_age=True)

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
                   truth_age: bool = False, truth_loss_split: bool = False,
                   warming: Warming | None = None, truth_wall_law: str = "theta_w",
                   chloramine_truth=None) -> Scenario:
    """structural_noise: False, True (= "spec") or "persistent"; see apply_structural_noise.
    month_seed: if given, the pipe-level truth (per-pipe wall decay, roughness) comes from `seed` and the
    operating truth (bulk decay, demand, dose) from `month_seed`: the same network in a different month.
    chem: None is the committed truth, draw for draw.  A chemistry.Chemistry gives the truth that chemistry
    (see _chem_truth); it consumes the same rng / rng_m draws in the same order, and every new draw comes from
    its own generator, default_rng(20_000 + seed) for free chlorine or default_rng(40_000 + seed) for chloramine
    (task 12: chloramine.chloramine_truth).
    truth_age: also run EPANET AGE on the truth's perturbed network after all draws (no new draws), and once more
    with every junction's and tank's initial age raised, for truth_initial_share (see initial_water_share).
    truth_loss_split: also rerun the truth's own network with wall decay off, with bulk decay off and with no decay
    (three runs after all draws, no new draws), so age.truth_loss_split can split its chlorine loss between the
    water and the pipe walls.  Neither option changes the chlorine truth.
    warming: a chemistry.Warming (needs chem with a temp_C, the plant temperature): every pipe and tank of the truth
    decays at its own temperature between the plant's and the soil's (see warming_temperatures).  It makes no draw,
    so the truth with warming is paired draw for draw with the truth without it, and with the soil at the plant
    temperature it is that truth bit for bit (a saved check).  Truth only: the model is told the plant temperature.
    truth_wall_law: how the truth's wall rate moves with temperature (needs chem).  'theta_w' (the default, the
    pre-registered task-10 truth) is theta_w^(T - 20); 'arrhenius' is f(T; E_true), the bulk factor (the structure
    the temperature bank M assumes); 'none' is 1 (the structure the bulk-only ablation M1b assumes).  The last two
    exist only for the task-10 mechanism checks run after the stop (seasonal.mechanism_checks); no draw changes.
    chloramine_truth: task 12, with chem=Chemistry(disinfectant='chloramine', kinetics='epa_msx' or 'first'): a
    chloramine.ChloramineTruth (pH range, nitrification stress); None is the default chloramine truth.  The truth is
    total chlorine (mg/L as Cl2); source_dose is its nominal plant dose and kw_m_per_day its wall rate at C = 130 (see
    chloramine.chloramine_truth).  Warming, another wall law and truth_loss_split are refused for chloramine."""
    if warming is not None and (chem is None or chem.temp_C is None):
        raise ValueError("warming needs chem=Chemistry(temp_C=<plant temperature>)")
    if truth_wall_law not in TRUTH_WALL_LAWS:
        raise ValueError(f"truth_wall_law must be one of {TRUTH_WALL_LAWS}")
    if truth_wall_law != "theta_w" and chem is None:
        raise ValueError("truth_wall_law needs chem=Chemistry(...)")
    is_ca = chem is not None and chem.disinfectant == CHLORAMINE
    if chloramine_truth is not None and not is_ca:
        raise ValueError("chloramine_truth needs chem=Chemistry(disinfectant='chloramine', ...)")
    if is_ca and (warming is not None or truth_wall_law != "theta_w" or truth_loss_split):
        raise NotImplementedError("in-network warming, another wall law and the loss split are not built for chloramine")
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
        q, chem_info = _chem_truth(wn, seed, rng, rng_m, chem, source_dose, kb_per_day, kw_m_per_day, warming, name,
                                   truth_wall_law, chloramine_truth)
    junctions = wn.junction_name_list
    truth_by_hour = _last_day(q, junctions).clip(lower=0.0)
    loss_runs = None
    if truth_loss_split:
        scale = 1.0 if chem_info is None else chem_info["quality_scale"]
        loss_runs = _truth_loss_runs(wn, _last_day(q, junctions), scale)
    truth_age_h, truth_share = _truth_age(wn, *age_state) if truth_age else (None, None)

    sc = nominal_scenario(name, sample_hour, source_dose)
    sc.seed, sc.structural = seed, structural
    sc.truth_snapshot, sc.truth_daily_min, sc.truth_by_hour = truth_by_hour.loc[sample_hour], truth_by_hour.min(), truth_by_hour
    sc.truth_age_by_hour_h, sc.chem, sc.truth_loss_runs = truth_age_h, chem_info, loss_runs
    sc.truth_initial_share = truth_share
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
                                   only the direction, wall decay rising with temperature)
    and, for chloramine only (task 12), four more uniforms on [0, 1), drawn after those two so neither moves:
      u_pH, u_cl2n, u_toc, u_alk   mapped by chloramine.truth_chemistry onto the pH, Cl2:N, TOC and alkalinity ranges
                                   (a stress subset maps u_pH onto a lower pH range: the same draw, another range)."""
    rng_seed = TRUTH_SEED_OFFSET[disinfectant] + seed
    rc = np.random.default_rng(rng_seed)
    er_true = float(rc.uniform(*ER_TRUTH_RANGE_K))
    theta_w = float(rc.uniform(*THETA_W_TRUTH_RANGE))
    out = {"rng_seed": rng_seed, "E_true_K": er_true, "theta_w": theta_w}
    if disinfectant == CHLORAMINE:
        out.update({k: float(rc.random()) for k in ("u_pH", "u_cl2n", "u_toc", "u_alk")})
    return out


@functools.lru_cache(maxsize=16)
def _nominal_age_and_flow(name: str, duration_days: int = DURATION_DAYS) -> tuple[dict, dict]:
    """Daily-mean water age (h) at every node (junctions, tanks, reservoirs) and daily-mean flow (m3/s) of every
    pipe on the operator's NOMINAL file, last day of the run.  One EPANET AGE run in a temporary directory, cached
    per process (the file does not change).  Returned as plain dicts; callers must not modify them."""
    wn = load(name, duration_days)
    wn.options.quality.parameter = "AGE"
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as tmp:
        prefix = os.path.join(tmp, "run")
        try:
            res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=prefix)
        finally:
            remove_epanet_files(prefix)
    age = (_last_day(res.node["quality"], wn.node_name_list) / 3600.0).mean()
    flow = _last_day(res.link["flowrate"], wn.pipe_name_list).mean()
    return {k: float(v) for k, v in age.items()}, {k: float(v) for k, v in flow.items()}


def warming_temperatures(name: str, plant_C: float, warming: Warming, wn=None) -> dict:
    """Water temperature in every pipe and tank of the truth under in-network warming (task 10, variant V3):
      pipe p : T_p = T_soil + (T_plant - T_soil) exp(-a_p / tau_pipe), a_p the nominal daily-mean age at the pipe's
               downstream node (by the nominal daily-mean flow; the file's end node when the flow is zero)
      tank t : the same with the tank's own nominal daily-mean age and tau_tank
      mixed  : T_soil + (T_plant - T_soil) x (the |flow|-weighted mean over pipes of exp(-a_p / tau_pipe)), the single
               temperature used for EPANET's global VISCOSITY and DIFFUSIVITY options.
    Written so every temperature is exactly T_plant when T_soil == T_plant.  Ages are the operator's nominal ones (the
    truth's own perturbed network is not used), from _nominal_age_and_flow."""
    age, flow = _nominal_age_and_flow(name)
    wn = wn if wn is not None else load(name)
    pipes, weights, e_pipe = {}, [], []
    for pn, p in wn.pipes():
        down = p.start_node_name if flow.get(pn, 0.0) < 0 else p.end_node_name
        e = math.exp(-max(age.get(down, 0.0), 0.0) / warming.tau_pipe_h)
        pipes[pn] = warming.soil_temp_C + (float(plant_C) - warming.soil_temp_C) * e
        weights.append(abs(flow.get(pn, 0.0)))
        e_pipe.append(e)
    tanks = {tn: warming.temperature(plant_C, age.get(tn, 0.0), warming.tau_tank_h) for tn, _ in wn.tanks()}
    wsum = float(np.sum(weights))
    e_mix = float(np.dot(weights, e_pipe) / wsum) if wsum > 0 else 1.0
    mixed = warming.soil_temp_C + (float(plant_C) - warming.soil_temp_C) * e_mix
    return {"pipes": pipes, "tanks": tanks, "mixed": mixed}


def _chem_truth(wn, seed, rng, rng_m, chem: Chemistry, source_dose, kb_per_day, kw_m_per_day,
                warming: Warming | None = None, name: str | None = None, wall_law: str = "theta_w",
                chloramine_truth=None):
    """The hidden truth under a chemistry condition.  The committed draws are consumed exactly as in the
    default branch of build_scenario (bulk factor, per-pipe wall and roughness, global and per-node demand,
    dose per source); the new ones come from hidden_chem_draws (their own generator).
    At T = 20 C (or None), TOC = TOC_ref (or None) and first order, every factor is exactly 1 and the truth
    equals build_scenario(chem=None) bit for bit (a saved check).  kinetics 'first_si' runs the same first
    order with the corrected unit recipe (chemistry.use_mg_per_litre); 'order2' is not built.  A chloramine condition
    (task 12) goes to chloramine.chloramine_truth, which consumes the same committed draws in the same order.
    kinetics 'clark' (task 11, Clark 1998, J Environ Eng 124(1):16): the bulk and tank reaction is EPANET order 2
    with a limiting potential, dC/dt = -k2 C (C - CL), with CL = d - phi TOC (d the mean of the month's drawn source
    doses, one CL per run, so one plant dose is assumed) and k2 = kb_per_day u f(T; E_true) / (phi TOC_ref): at
    TOC_ref the initial apparent first-order rate k2 (d - CL) equals the committed truth's kb u f.  The decay slows
    as the water ages, the residual responds nonlinearly to the dose, and when CL > 0 the bulk demand runs out.  The
    wall stays first order as committed (it consumes chlorine, not organics: Clark's invariant does not hold at the
    walls, which is acceptable for a truth generator).  Corrected units (use_mg_per_litre); no new draw.  In-network
    warming with Clark kinetics is refused (not built).
    warming (task 10, V3): every pipe gets its own bulk coefficient and wall factor at its own temperature, every
    tank its own bulk coefficient, and the viscosity and diffusivity options are set at the flow-weighted mixed
    temperature (warming_temperatures); no extra draw is made."""
    chem.require_built()
    if chem.disinfectant == CHLORAMINE:      # task 12: its own physics, the same committed draws (chloramine.py)
        from .chloramine import chloramine_truth as _ca_truth
        return _ca_truth(wn, seed, rng, rng_m, chem, source_dose, kw_m_per_day, name, chloramine_truth)
    if chem.kinetics not in ("first", "first_si", "clark"):
        raise NotImplementedError(f"truth kinetics {chem.kinetics!r} is not built")
    if chem.kinetics == "clark" and warming is not None:
        raise NotImplementedError("in-network warming with Clark kinetics is not built")
    hd = hidden_chem_draws(chem.disinfectant, seed)
    er_true, theta_w = hd["E_true_K"], hd["theta_w"]
    T = chem.temp_C
    fb = 1.0 if T is None else arrhenius(T, er_true)
    def wall_factor(t):            # build_scenario's truth_wall_law; 'theta_w' is the pre-registered truth
        if wall_law == "theta_w":
            return theta_w ** (t - TREF_C)
        return arrhenius(t, er_true) if wall_law == "arrhenius" else 1.0

    fw = 1.0 if T is None else wall_factor(T)
    tr = toc_ratio(chem.toc_mgL)
    wt = None
    if warming is not None:
        if T is None or name is None:
            raise ValueError("warming needs the plant temperature (chem.temp_C) and the network name")
        wt = warming_temperatures(name, T, warming, wn)

    u = rng_m.uniform(0.8, 1.2)
    wn.options.reaction.bulk_coeff = -kb_per_day * u * fb * tr / DAY
    wn.options.reaction.wall_coeff = -kw_m_per_day * fw / DAY
    for pn, pipe in wn.pipes():
        fw_p = fw if wt is None else wall_factor(wt["pipes"][pn])
        pipe.wall_coeff = -kw_m_per_day * fw_p * roughness_factor(pipe.roughness, 1.0) * np.exp(rng.normal(0, 0.4)) / DAY
        pipe.roughness = pipe.roughness * np.exp(rng.normal(0, 0.10))
        if wt is not None:
            pipe.bulk_coeff = -kb_per_day * u * arrhenius(wt["pipes"][pn], er_true) * tr / DAY
    if wt is not None:
        for tn, tank in wn.tanks():
            tank.bulk_coeff = -kb_per_day * u * arrhenius(wt["tanks"][tn], er_true) * tr / DAY
    global_mult = rng_m.uniform(0.85, 1.15)
    for _, j in wn.junctions():
        for ts in j.demand_timeseries_list:
            ts.base_value = ts.base_value * global_mult * np.exp(rng_m.normal(0.0, 0.15))
    doses = {res: source_dose * rng_m.uniform(0.9, 1.1) for res in source_nodes(wn)}
    v_ratio, d_ratio = apply_water_temperature(wn, T if wt is None else wt["mixed"])
    if chem.kinetics == "first":
        for res, dose in doses.items():
            wn.add_source(f"src_{res}", res, "CONCEN", dose)
        scale = 1.0
    else:
        scale = use_mg_per_litre(wn, doses)
    clark = None
    if chem.kinetics == "clark":     # task 11: replaces the first-order bulk and tank coefficient set above
        # rounded to the 4 decimals WNTR's .inp writer keeps, so the record is exactly what EPANET runs
        k2 = round(kb_per_day * u * fb / (chem.phi * TOC_REF_MGL), 4)
        cl = round(float(np.mean(list(doses.values()))) - chem.phi * chem.toc_mgL, 4)
        set_bulk_kinetics(wn, 2, -k2, cl)
        clark = {"phi_mg_per_mgC": chem.phi, "k2_L_per_mg_day": k2, "limiting_mgL": cl,
                 "bulk_demand_mgL": chem.phi * chem.toc_mgL}
    q = _run_quality(wn) * scale
    info = {"disinfectant": chem.disinfectant, "species": chem.species, "kinetics": chem.kinetics,
            "temp_C": T, "toc_mgL": chem.toc_mgL, "rng_seed": hd["rng_seed"],
            "E_true_K": er_true, "theta_w": theta_w, "bulk_month_factor": float(u),
            "bulk_temp_factor": fb, "bulk_toc_factor": tr, "wall_temp_factor": fw,
            "viscosity_ratio": v_ratio, "diffusivity_ratio": d_ratio,
            "kb_per_day_effective": kb_per_day * u * fb * tr,
            "source_doses_mgL": {k: float(v) for k, v in doses.items()}, "quality_scale": scale}
    if clark is not None:   # TOC acts through CL, not as a rate factor; the initial apparent rate is k2 phi TOC
        info.update(bulk_toc_factor=None, kb_per_day_effective=clark["k2_L_per_mg_day"] * clark["bulk_demand_mgL"],
                    clark=clark)
    if wall_law != "theta_w":
        info["truth_wall_law"] = wall_law
    if wt is not None:
        tp = np.array(list(wt["pipes"].values()))
        info["warming"] = {"soil_temp_C": warming.soil_temp_C, "tau_pipe_h": warming.tau_pipe_h,
                           "tau_tank_h": warming.tau_tank_h, "mixed_temp_C": wt["mixed"],
                           "pipe_temp_C_min": float(tp.min()), "pipe_temp_C_median": float(np.median(tp)),
                           "pipe_temp_C_max": float(tp.max()), "tank_temp_C": dict(wt["tanks"])}
    return q, info


def _truth_loss_runs(wn, full: pd.DataFrame, scale: float) -> dict:
    """The truth's own network (after all draws) rerun three times, each in a temporary directory: with every wall
    coefficient 0 ('bulk_only'), with the bulk coefficient 0 ('wall_only') and with both 0 ('no_decay', the
    reference that carries each source's own dose).  `full` is the truth run itself (last day, unclipped).  Every
    coefficient is restored afterwards, so a truth_age run that follows sees the network as it was.  Returns
    hour x junction frames in the same units as the truth (multiplied by `scale`)."""
    rx = wn.options.reaction
    bulk0, wall0 = rx.bulk_coeff, rx.wall_coeff
    pipe_wall = {n: p.wall_coeff for n, p in wn.pipes()}
    pipe_bulk = {n: p.bulk_coeff for n, p in wn.pipes()}      # None unless set per pipe (in-network warming)
    tank_bulk = {n: t.bulk_coeff for n, t in wn.tanks()}
    runs = {"full": full}
    try:
        for tag, bulk_on, wall_on in (("bulk_only", True, False), ("wall_only", False, True), ("no_decay", False, False)):
            rx.bulk_coeff = bulk0 if bulk_on else 0.0
            rx.wall_coeff = wall0 if wall_on else 0.0
            for n, p in wn.pipes():
                p.wall_coeff = pipe_wall[n] if wall_on else 0.0
                if pipe_bulk[n] is not None:
                    p.bulk_coeff = pipe_bulk[n] if bulk_on else 0.0
            for n, t in wn.tanks():
                if tank_bulk[n] is not None:
                    t.bulk_coeff = tank_bulk[n] if bulk_on else 0.0
            runs[tag] = _last_day(_run_quality(wn), wn.junction_name_list) * scale
    finally:
        rx.bulk_coeff, rx.wall_coeff = bulk0, wall0
        for n, p in wn.pipes():
            p.wall_coeff = pipe_wall[n]
            p.bulk_coeff = pipe_bulk[n]
        for n, t in wn.tanks():
            t.bulk_coeff = tank_bulk[n]
    return runs


def _truth_age(wn, tolerance, initial_quality) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Water age (hours, hour x junction, last day) on the truth's own perturbed network, and its initial-water
    share (initial_water_share).  The file's tolerance and initial qualities (initial ages) are restored first,
    as the nominal AGE run uses them."""
    wn.options.quality.parameter = "AGE"
    wn.options.quality.tolerance = tolerance
    for n, q0 in initial_quality.items():
        wn.get_node(n).initial_quality = q0
    age = _last_day(_run_quality(wn), wn.junction_name_list) / 3600.0
    return age, initial_water_share(wn, age)


INITIAL_AGE_OFFSET_H = 1000.0      # added to every junction's and tank's initial age by initial_water_share


def initial_water_share(wn, age_h: pd.DataFrame) -> pd.DataFrame:
    """Share of each junction's water on the last day (hour x junction, 0 to 1) that is still the water the
    7-day run started with in its pipes, junctions and tanks, not water that entered from a source during the
    run.  Where it is above zero the run's water age is a lower bound: that water is older than the run.

    Water age mixes linearly, so rerunning AGE with every junction's and tank's initial age raised by
    INITIAL_AGE_OFFSET_H raises each junction's age by exactly that offset times this share (sources keep
    their own age).  wn must be the network, set up for AGE, that gave age_h; its initial ages are restored
    afterwards.  One extra EPANET run, in a temporary directory."""
    nodes = [n for n, _ in wn.junctions()] + [n for n, _ in wn.tanks()]
    q0 = {n: wn.get_node(n).initial_quality for n in nodes}
    try:
        for n in nodes:
            wn.get_node(n).initial_quality = q0[n] + INITIAL_AGE_OFFSET_H * 3600.0
        raised = _last_day(_run_quality(wn), wn.junction_name_list) / 3600.0
    finally:
        for n in nodes:
            wn.get_node(n).initial_quality = q0[n]
    return ((raised - age_h) / INITIAL_AGE_OFFSET_H).clip(0.0, 1.0)


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
