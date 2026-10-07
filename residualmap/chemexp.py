"""
chemexp.py: task 13, the richer-truth audit (iteration 4).  What does today's first-order free-chlorine model (B0, the
committed SimGP24 on the committed 675-member grid) cost against chemistry its structure cannot represent?

The richer truth is Fisher's two-reactant chlorine-organics chemistry (2RA: a fast and a slow pool of organic
reactants; Fisher, Kastl & Sathasivan 2012, doi:10.1016/j.watres.2012.03.017), run in EPANET-MSX on the committed
draws (msx.two_reactant_truth), at fixed conditions (20 C, TOC at its reference).  Its apparent first-order rate falls
as the water ages and the residual responds nonlinearly to the dose, neither of which a first-order grid can do.  The
plan's seasonal variant (2ra_warm) and its cost decomposition were dropped with the temperature models (plan,
addendum 2).  Everything is simulation; no real grab-sample data exists.

    python -m residualmap.chemexp scaling                             # outputs/chem_2ra/scaling.json (truth tuning, s)
    python -m residualmap.experiment Net3 8 --chemistry=first_order   # the paired experiments (experiment.py), fresh
    python -m residualmap.experiment Net3 8 --chemistry=2ra           # seeds 400 to 407 -> outputs/chem_2ra/
    python -m residualmap.experiment Net2 8 --chemistry=2ra --match=96   # Net2 sensitivity: s matched at 96 h
    python -m residualmap.chemexp audit Net3 8 --workers 6            # oracle, age deciles, dose step (Net3), the
                                                                      # opt-in low-kb grid -> outputs/chem/, chem_2ra/
    python -m residualmap.chemexp summarise Net3                      # outputs/chem/summary_audit_Net3.json, figures
    python -m residualmap.chemexp full2r Net3 8 --workers 6           # the gated two-rate model (only if material)
    python -m residualmap.chemexp near Net3 8 --workers 4             # after the review, no bar: junction-days just
                                                                      # above 0.2 mg/L per truth (before summarise)

Parts (the bars are pre-registered in docs/iteration3_journal.md, task 13):
  * materiality: B0's daily-minimum RMSE and recall on the 2RA truth against the committed first-order truth, paired
    seed by seed on the same fresh seeds (the same draws; the random rule's samples are the same junctions, hours and
    reading errors), at n = 8 and 15 under the random and straddle_min rules, from the two results_time CSVs;
  * the oracle: the best single member and dose of a grid against the full daytime truth, scored on the daily minimum
    of every junction (how well the grid's structure can represent each truth, with no discrepancy GP);
  * error by water-age decile: B0 after task 9's 8 demo samples (the app's demo draw), by decile of the nominal
    daily-mean age, against the paired first-order truth and task 9's committed baseline;
  * the dose step (Net3): B0 fitted on 15 demo samples at the 1.2 mg/L dose, its daily-minimum map moved to doses
    0.9 and 1.6 by the exact ln-offset first order implies, scored against each truth rerun at that dose with every
    draw the same (on the first-order truth the offset is exact, so that is the twin);
  * the opt-in low-kb grid (plan addendum 3): B0 with the bulk-rate grid extended below its 0.10 per day floor
    (simgp.GRIDS['full_lowkb'], its own cache name), on every truth, with the committed rules and samples;
  * the gated two-rate model ('full2r'), built only if the cost is material: C = f C(k1) + (1 - f) C(k2) assembled from
    the committed grid with no new run (exact for EPANET's linear first-order operator above the grid's 0.02 mg/L
    floor: the grid stores ln max(C, FLOOR), so the mixture is of floored values, not masked below FLOOR as the plan
    said), plus the committed members.
"""
from __future__ import annotations

import argparse
import contextlib
import json
import math
import os
import time
import warnings

import numpy as np
import pandas as pd

from .experiment import CHEM_FIRST_SEED, NET_CHEM_MATCH_H, NET_TRUTH, chem_outdir, run_scenario_time
from .simgp import (DAY_HOURS, DOSE_GRID, FLOOR, HOURS, KB_GRID, KB_LOWKB, SimGP24, check_grid_order, grid_cache_path,
                    n_hydraulic, simulator_grid_24h)

warnings.filterwarnings("ignore")      # GP optimiser bound warnings, as in experiment.py
THRESHOLD = 0.2
N_SEEDS = 8
SMOKE_SEEDS = (990, 991)               # code checks before the pre-registration; never scored
TRUTHS = {"first_order": None, "2ra": NET_CHEM_MATCH_H, "2ra_m96": 96.0}   # truth -> 2RA match age (None: committed)
NET_TRUTHS = {"Net3": ("first_order", "2ra"), "Net2": ("first_order", "2ra", "2ra_m96")}
PRIMARY = ("first_order", "2ra")       # the materiality pair; 2ra_m96 is the plan's Net2 sensitivity, no bar
DOSE_NETS = ("Net3",)
DOSE_RATIOS = (0.75, 4.0 / 3.0)        # 1.2 mg/L x 0.75 = 0.9 and x 4/3 = 1.6 (the plan's x0.75 and x1.33)
AGE_N, DOSE_N = 8, 15
CELL_RULES, CELL_NS = ("random", "straddle_min"), (8, 15)
MATERIAL_RMSE, MATERIAL_RECALL = 1.10, 0.90
AGE_BAR_MGL = 0.05
COVERAGE_BAND = (0.85, 0.97)
BOOT_N, BOOT_SEED = 2000, 2026
OUT_CHEM = os.path.join("outputs", "chem")
OUT_2RA = os.path.join("outputs", "chem_2ra")
CSV_FLOAT = "%.6g"
SUM_COLS = ("n", "sse", "sae", "se", "in90", "n_low", "tp", "fp")
SINGLE_THREAD_ENV = {k: "1" for k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS",
                                      "MKL_NUM_THREADS")}


# ----------------------------------------------------------------------------- truths
def truth_kwargs(net: str, truth: str) -> dict:
    """build_scenario keywords for a truth on a network: NET_TRUTH's committed settings, plus the 2RA chemistry."""
    if truth not in TRUTHS:
        raise ValueError(f"truth must be one of {tuple(TRUTHS)}")
    kw = dict(NET_TRUTH.get(net, {}))
    if TRUTHS[truth] is not None:
        from .chemistry import Chemistry
        from .msx import TwoReactantTruth
        kw.update(chem=Chemistry(kinetics="2ra"), two_reactant=TwoReactantTruth(match_h=TRUTHS[truth]))
    return kw


def scenario(net: str, seed: int, truth: str, dose_ratio: float = 1.0):
    """The hidden truth (and the operator's side) for one network, seed and truth; dose_ratio scales the nominal plant
    dose (every draw the same, so the per-source dose factors are the same)."""
    from .simulate import build_scenario
    kw = truth_kwargs(net, truth)
    dose = float(kw.pop("source_dose", 1.2)) * float(dose_ratio)
    return build_scenario(net, seed=seed, source_dose=dose, **kw)


def scaling() -> dict:
    """The truth tuning, disclosed: s per network and match age, the reactant loads at u = 1, the batch's apparent
    first-order rate by age at the 1.2 mg/L dose, its local rate in old water, and the batch's ln C shift for the dose
    step's doses against the shift first order implies.  Written to outputs/chem_2ra/scaling.json."""
    from . import msx as M
    from .simulate import build_scenario
    import inspect
    kb_default = inspect.signature(build_scenario).parameters["kb_per_day"].default
    out = {"generated_by": "python -m residualmap.chemexp scaling",
           "about": "Truth tuning for task 13's two-reactant (2RA) truth, disclosed. s is solved by brentq so that the "
                    "2RA batch's apparent first-order rate -ln(C(match)/C0)/match at C0 = 1.2 mg/L and 20 C equals the "
                    "network's committed truth bulk rate (Net3 0.40, Net2 0.10 per day). Batch = plug flow, no wall, "
                    "u = 1 (the truth multiplies both reactant loads by its monthly draw u ~ U(0.8, 1.2)).",
           "constants": {"kfast_L_per_mg_h": M.KF20_L_PER_MG_H, "kslow_L_per_mg_h": M.KS20_L_PER_MG_H,
                         "fast0_greenvale_mgL": M.F0_GREENVALE_MGL, "slow0_greenvale_mgL": M.S0_GREENVALE_MGL,
                         "c0_mgL": M.TWO_RA_C0_MGL, "temp_C": 20.0,
                         "source": "Greenvale water, Fisher, Kastl & Sathasivan 2012 (doi:10.1016/j.watres.2012.03.017), "
                                   "via Walski's Bentley blog (secondary)"},
           "networks": {}}
    ages = (2, 6, 12, 24, 48, 96, 144)
    for net in ("Net3", "Net2"):
        kb = float(NET_TRUTH.get(net, {}).get("kb_per_day", kb_default))
        for truth in ([t for t in NET_TRUTHS[net] if TRUTHS[t] is not None]):
            match = TRUTHS[truth]
            s = M.two_reactant_scale(kb, match)
            f0, s0 = M.F0_GREENVALE_MGL * s, M.S0_GREENVALE_MGL * s
            c, _, _ = M.two_reactant_batch(1.2, f0, s0, [96.0, 120.0, 144.0, 168.0])
            local = (-np.diff(np.log(c)) / 24.0 * 24.0).tolist()
            shift = {}
            for d in (0.9, 1.6):
                a, _, _ = M.two_reactant_batch(d, f0, s0, [6.0, 24.0, 48.0, 96.0])
                b, _, _ = M.two_reactant_batch(1.2, f0, s0, [6.0, 24.0, 48.0, 96.0])
                shift[f"{d:g}"] = {"ages_h": [6, 24, 48, 96], "ln_C_dose_over_C_1.2": np.log(a / b).round(4).tolist(),
                                   "first_order_ln_ratio": round(math.log(d / 1.2), 4)}
            out["networks"].setdefault(net, {})[truth] = {
                "kb_matched_per_day": kb, "match_h": match, "scale_s": s, "fast0_mgL_at_u1": f0, "slow0_mgL_at_u1": s0,
                "apparent_rate_per_day_by_age_h": {str(a): round(M.two_reactant_apparent_rate(s, a), 4) for a in ages},
                "local_rate_per_day_96_to_168h": [round(v, 4) for v in local],
                "dose_shift": shift}
    return out


# ----------------------------------------------------------------------------- scoring helpers
def _sums(t: pd.Series, med: pd.Series, lo: pd.Series | None, hi: pd.Series | None, flag: pd.Series) -> dict:
    """Sums over a set of junctions, so rates pool exactly over seeds (task 12's convention)."""
    e = (med - t).values
    tv = (t < THRESHOLD).values
    fl = flag.values.astype(bool)
    return {"n": int(len(t)), "sse": float((e ** 2).sum()), "sae": float(np.abs(e).sum()), "se": float(e.sum()),
            "in90": float(((t >= lo) & (t <= hi)).sum()) if lo is not None else float("nan"),
            "n_low": int(tv.sum()), "tp": int((tv & fl).sum()), "fp": int((~tv & fl).sum())}


def rates(d: dict) -> dict:
    """RMSE, MAE, bias, coverage90, recall, precision, F1 and false alarms from pooled sums.  A rate with nothing to
    find (no junction below the threshold) or nothing flagged is NaN ('not testable'), never 1.0."""
    n = d["n"]
    rec = d["tp"] / d["n_low"] if d["n_low"] else float("nan")
    flagged = d["tp"] + d["fp"]
    prec = d["tp"] / flagged if flagged else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if (prec == prec and rec == rec and prec + rec > 0) else float("nan")
    return {"rmse": math.sqrt(d["sse"] / n), "mae": d["sae"] / n, "bias": d["se"] / n,
            "coverage90": d["in90"] / n if d["in90"] == d["in90"] else float("nan"), "recall": rec, "precision": prec,
            "f1": f1, "false_alarms": int(d["fp"]), "n_low": int(d["n_low"]), "found": int(d["tp"]), "n": int(n)}


def oracle(sc, params, Z, doses=DOSE_GRID) -> dict:
    """The best single member and dose of a grid against the full daytime truth (every junction, hours 07:00 to
    17:00, ln C with the grid's FLOOR), and its daily-minimum map scored on every junction.  No discrepancy GP: this is
    how well the grid's structure can represent the truth at its best."""
    zt = np.log(np.clip(sc.truth_by_hour.loc[HOURS, sc.junctions].values, FLOOR, None))      # 24 x J
    day = list(DAY_HOURS)
    night = [h for h in HOURS if h not in day]
    offs = np.log(np.asarray(doses, dtype=float))
    best = (None, None, np.inf)
    for d in offs:
        r = np.sqrt(((Z[:, day, :] + d - zt[None, day, :]) ** 2).mean(axis=(1, 2)))
        k = int(np.argmin(r))
        if r[k] < best[2]:
            best = (k, float(d), float(r[k]))
    k, d, rday = best
    e = Z[k] + d - zt
    pred_min = pd.Series(np.exp(Z[k] + d).min(axis=0), index=sc.junctions)
    t = sc.truth_daily_min.loc[sc.junctions]
    s = _sums(t, pred_min, None, None, pred_min < THRESHOLD)
    kb, kw, g, dm, rm = params[k]
    return {**s, "member_kb": kb, "member_kw": kw, "member_gamma": g, "member_demand": dm, "member_rough": rm,
            "member_dose": float(np.exp(d)), "rms_ln_day": rday, "rms_ln_night": float(np.sqrt((e[night] ** 2).mean()))}


# ----------------------------------------------------------------------------- models
def lowkb_model(sc, X, seed, cache_dir, threshold):
    """B0 on the opt-in low-kb grid: everything else (likelihood, dose axis, discrepancy GP, Monte Carlo) is today's."""
    return SimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=threshold, grid="full_lowkb")


def two_rate_members(kb=KB_GRID, fractions=(0.25, 0.5, 0.75)) -> list[tuple]:
    """The decay-rate hypotheses of 'full2r', in grid order: today's single rates (k, None, 1.0), then every pair
    k1 > k2 of today's bulk-rate grid (10 pairs) at each fast-share f, (k1, k2, f)."""
    out = [(float(k), None, 1.0) for k in kb]
    pairs = [(float(k1), float(k2)) for k1 in kb for k2 in kb if k1 > k2]
    for k1, k2 in pairs:
        for f in fractions:
            out.append((k1, k2, float(f)))
    return out


_TWO_RATE_CACHE: dict = {}


def two_rate_grid(params: list[tuple], Z: np.ndarray, key=None) -> tuple[list[tuple], list[tuple], np.ndarray]:
    """'full2r', assembled from the committed 'full' grid with no new EPANET run.  A two-rate member is
    C = f C(k1) + (1 - f) C(k2) for the same (kw, gamma, demand, roughness): the source water is split into two parts
    that decay in the bulk at k1 and k2 and at the walls alike, which is exact for EPANET's linear first-order operator
    (Vieira, Coelho & Loureiro 2004, doi:10.2166/aqua.2004.0036; Al Heboos & Licsko 2017).  The committed grid stores
    ln max(C, FLOOR), so the mixture inherits that floor (it is never below FLOOR).  Members are decay-major and
    hydraulic-minor like every grid here; the first 675 are today's grid (so B0 is nested).  The plan said the members
    are 'masked where C < FLOOR'; here they are not masked: where the fast component sits at FLOOR and the slow one does
    not, the mixture is overstated by at most f x FLOOR (0.005 to 0.015 mg/L), and where both do it is FLOOR, as for a
    single-rate member.  So it is exact only above the floor (a deviation recorded in the journal).
    Returns (members, params,
    Z): members[i] = (k1, k2 or None, f, kw, gamma, demand, rough); params are 5-tuples whose first entry is k for a
    single-rate member and -(1 + index of its (k1, k2, f) hypothesis) for a two-rate one (a code, not a rate)."""
    if key is not None and key in _TWO_RATE_CACHE:
        return _TWO_RATE_CACHE[key]
    P = np.asarray(params, dtype=float)
    kbs = sorted(set(P[:, 0]))
    if kbs != sorted(map(float, KB_GRID)):
        raise ValueError("full2r is assembled from the committed 'full' grid (KB_GRID bulk rates)")
    block = {kb: np.where(P[:, 0] == kb)[0] for kb in kbs}
    nb = len(block[kbs[0]])
    if any(len(v) != nb for v in block.values()):
        raise ValueError("unequal bulk-rate blocks")
    sub = [tuple(p[1:]) for p in P[block[kbs[0]]]]
    for kb in kbs:
        if [tuple(p[1:]) for p in P[block[kb]]] != sub:
            raise ValueError("the bulk-rate blocks do not share one (kw, gamma, demand, rough) order")
    hyps = two_rate_members(kbs)
    members, params2, Zs = [], [], []
    for h_i, (k1, k2, f) in enumerate(hyps):
        if k2 is None:
            Zh = Z[block[k1]]
            code = k1
        else:
            Zh = np.log(f * np.exp(Z[block[k1]].astype(np.float64)) + (1 - f) * np.exp(Z[block[k2]].astype(np.float64)))
            Zh = Zh.astype(np.float32)
            code = -(1.0 + h_i)
        Zs.append(Zh)
        for p in sub:
            members.append((k1, k2, f, *p))
            params2.append((code, *p))
    out = (members, params2, np.concatenate(Zs, axis=0))
    if key is not None:
        _TWO_RATE_CACHE.clear()
        _TWO_RATE_CACHE[key] = out
    return out


class TwoRateSimGP24(SimGP24):
    """The gated 'full2r' model: today's SimGP24 (likelihood, dose axis, discrepancy GP, daily-minimum Monte Carlo,
    uniform prior over members) on two_rate_grid.  Built only if the fixed-conditions cost is material.
    The in-process cache of the assembled grid is keyed on the committed grid's cache file (its path carries the
    network and the plant-dose tag) and that file's modification time, so a scenario at another dose, or a rebuilt
    file, never reuses another grid's mixtures (review fix; every committed full2r run was at the nominal dose)."""

    def __init__(self, sc, X, cache_dir: str = "outputs/cache", **kw):
        kw.pop("grid", None)
        super().__init__(sc, X, grid="full", cache_dir=cache_dir, **kw)
        path = os.path.abspath(grid_cache_path(sc, cache_dir, "full"))
        key = (path, os.path.getmtime(path) if os.path.exists(path) else None, self.Z.shape)
        self.members_, self.params, self.Z = two_rate_grid(self.params, self.Z, key=key)
        check_grid_order(self.params, self.n_hyd)
        self._grp_mean = None


def full2r_model(sc, X, seed, cache_dir, threshold):
    return TwoRateSimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=threshold)


MODEL_FACTORIES = {"lowkb": lowkb_model, "full2r": full2r_model}


# ----------------------------------------------------------------------------- one (network, seed, truth)
def _age_frame(sc, X, seed: int, cache_dir: str, model: str) -> pd.DataFrame:
    """Task 9's protocol: 8 demo samples, the daily-minimum map at the unsampled junctions."""
    from .experiment_chem import demo_samples
    S = demo_samples(sc, seed, AGE_N)
    m = (SimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=THRESHOLD) if model == "B0"
         else lowkb_model(sc, X, seed, cache_dir, THRESHOLD)).fit(S)
    pmin = m.predict_daily_min()
    uns = [j for j in sc.junctions if j not in set(S.junction)]
    return pd.DataFrame({"seed": seed, "junction": uns, "nominal_age_h": sc.age_daily_mean_h.loc[uns].values,
                         "true_age_h": np.nan, "true_min": sc.truth_daily_min.loc[uns].values,
                         "pred_median": pmin.loc[uns, "median"].values, "lo90": pmin.loc[uns, "lo90"].values,
                         "hi90": pmin.loc[uns, "hi90"].values, "flag": (pmin.loc[uns, "p_below"] > 0.5).values,
                         "map_kb": m.map_params_[0]})


def _dose_rows(net: str, seed: int, truth: str, sc, X, cache_dir: str) -> list[dict]:
    """B0 fitted on 15 demo samples at the nominal dose; its daily-minimum map moved to each dose by the exact
    ln-offset (median and band times the ratio, P(daily min x r < 0.2) from the model's own Monte-Carlo draws), scored
    on the unsampled junctions against the truth rerun at that dose (every draw the same)."""
    from .experiment_chem import demo_samples
    S = demo_samples(sc, seed, DOSE_N)
    m = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=THRESHOLD).fit(S)
    pmin = m.predict_daily_min()
    uns = [j for j in sc.junctions if j not in set(S.junction)]
    rows = []
    base_min = sc.truth_daily_min.loc[uns]
    for r in (1.0, *DOSE_RATIOS):
        t = base_min if r == 1.0 else scenario(net, seed, truth, r).truth_daily_min.loc[uns]
        flag = m.p_below_mc(THRESHOLD / r).loc[uns] > 0.5
        s = _sums(t, pmin.loc[uns, "median"] * r, pmin.loc[uns, "lo90"] * r, pmin.loc[uns, "hi90"] * r, flag)
        ok = (t > FLOOR) & (base_min > FLOOR)
        shift = np.log(t[ok] / base_min[ok])
        rows.append({"part": "dose_step", "model": "B0", "dose_ratio": round(r, 6), "dose_mgL": round(1.2 * r, 6), **s,
                     "truth_ln_shift_mean": float(shift.mean()) if len(shift) else float("nan"),
                     "truth_ln_shift_min": float(shift.min()) if len(shift) else float("nan"),
                     "truth_ln_shift_max": float(shift.max()) if len(shift) else float("nan"),
                     "first_order_ln_shift": math.log(r)})
    return rows


def run_task(net: str, seed: int, truth: str, cache_dir: str, parts=("oracle", "age", "dose", "lowkb")) -> dict:
    """Every audit part for one network, seed and truth.  Returns {'audit': rows, 'age': frame, 'lowkb': frame}."""
    t0 = time.time()
    from .features import build_features
    sc = scenario(net, seed, truth)
    X = build_features(sc)
    out = {"audit": [], "age": None, "lowkb": None, "full2r": None}
    base = {"net": net, "truth": truth, "seed": seed}
    if "oracle" in parts:
        for grid in ("full", "full_lowkb"):
            params, Z = simulator_grid_24h(sc, cache_dir, grid)
            out["audit"].append({**base, "part": "oracle", "model": grid, **oracle(sc, params, Z)})
    if "age" in parts:
        out["age"] = pd.concat([_age_frame(sc, X, seed, cache_dir, "B0").assign(truth=truth, model="B0"),
                                _age_frame(sc, X, seed, cache_dir, "lowkb").assign(truth=truth, model="lowkb")],
                               ignore_index=True)
    if "dose" in parts and net in DOSE_NETS and truth in PRIMARY:
        out["audit"] += [{**base, **r} for r in _dose_rows(net, seed, truth, sc, X, cache_dir)]
    for kind in ("lowkb", "full2r"):
        if kind in parts:
            df, _ = run_scenario_time(sc, X, seed=seed, cache_dir=cache_dir, make_model=MODEL_FACTORIES[kind],
                                      baselines=False, extra=True)
            out[kind] = df.assign(model=f"simgp24_{kind}", truth=truth)
    info = {k: v for k, v in (sc.chem or {}).items() if k != "_volatile"}
    out["truth_info"] = {**base, "n_low_daily_min": int((sc.truth_daily_min < THRESHOLD).sum()),
                         "mean_daily_min": float(sc.truth_daily_min.mean()), **{k: info.get(k) for k in
                         ("scale_s", "bulk_month_factor", "fast_mgL", "slow_mgL", "msx_compiler")}}
    out["seconds"] = time.time() - t0
    print(f"{net} {truth} seed {seed}: {', '.join(parts)} in {out['seconds']:.0f} s", flush=True)
    return out


def _task(args):
    return run_task(*args)


def _pool_init(workdir: str) -> None:
    """Worker: one thread per library and a private working directory (the committed truth and nominal_scenario write
    EPANET's temp files to the cwd; MSX changes the cwd for its runs)."""
    warnings.filterwarnings("ignore")
    os.environ.update(SINGLE_THREAD_ENV)
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)


class _single_thread_children:
    def __enter__(self):
        self.saved = {k: os.environ.get(k) for k in SINGLE_THREAD_ENV}
        os.environ.update(SINGLE_THREAD_ENV)

    def __exit__(self, *exc):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


@contextlib.contextmanager
def _scratch_cwd(prefix: str = "rm_2ra_parent_"):
    """Run a block of the parent process in a private temporary working directory and restore the old one: the
    committed truth and nominal_scenario write EPANET's temp.inp, temp.rpt and temp.bin to the cwd (review fix: the
    parent's grid prebuild and nominal_bins had left them in the folder the audit was started from).  Paths passed in
    must be absolute."""
    import tempfile
    old = os.getcwd()
    with tempfile.TemporaryDirectory(prefix=prefix) as tmp:
        os.chdir(tmp)
        try:
            yield tmp
        finally:
            os.chdir(old)


def _pool(tasks: list, workers: int, fn=None) -> list:
    import tempfile
    from concurrent.futures import ProcessPoolExecutor
    fn = fn or _task
    with tempfile.TemporaryDirectory(prefix="rm_2ra_") as tmp:
        if workers <= 1:
            old = os.getcwd()
            try:
                _pool_init(tmp)
                return [fn(t) for t in tasks]
            finally:
                os.chdir(old)
        with _single_thread_children(), ProcessPoolExecutor(max_workers=workers, initializer=_pool_init,
                                                            initargs=(tmp,)) as ex:
            return list(ex.map(fn, tasks))


def check_seeds(seeds, smoke_ok: bool = False) -> None:
    """The audit scores fresh seeds only (400 onwards); the smoke seeds are allowed only into a scratch folder."""
    bad = [s for s in seeds if s < CHEM_FIRST_SEED or (s in SMOKE_SEEDS and not smoke_ok)]
    if bad:
        raise ValueError(f"seeds {bad} are not the audit's fresh seeds")


def _write_csv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False, float_format=CSV_FLOAT)


def nominal_bins(net: str):
    """The nominal daily-mean age of every junction and task 9's lower-bound flag (more than 5% of the run's starting
    water), from the operator's file: the decile bins are task 9's."""
    from .age import hydraulic_age_band
    from .simulate import nominal_scenario
    sc = nominal_scenario(net)
    band = hydraulic_age_band(sc)
    return band.nominal.mean(), band.lower_bound("mean")


def audit(net: str, seeds, workers: int = 4, cache_dir: str = "outputs/cache", out_chem: str = OUT_CHEM,
          out_2ra: str = OUT_2RA, parts=("oracle", "age", "dose", "lowkb"), smoke_ok: bool = False) -> dict:
    """Runs every part for every (seed, truth) of the network in a pool and writes the compact CSVs:
    <out_chem>/audit_<net>.csv (oracle and dose-step sums), <out_chem>/error_by_age_2ra_<net>.csv (deciles, task 9's
    format) and <out_2ra>/results_lowkb_<net>.csv (the opt-in grid under the committed rules)."""
    from .experiment_chem import error_by_age
    check_seeds(seeds, smoke_ok)
    cache_dir = os.path.abspath(cache_dir)
    from .simulate import nominal_scenario
    with _scratch_cwd():                     # build any missing grid once, before the pool (each worker would build it)
        for grid in ("full", "full_lowkb"):
            simulator_grid_24h(nominal_scenario(net), cache_dir, grid)
    truths = PRIMARY if tuple(parts) == ("full2r",) else NET_TRUTHS[net]     # the gate needs the primary pair only
    tasks = [(net, int(s), truth, cache_dir, tuple(parts)) for truth in truths for s in seeds]
    t0 = time.time()
    res = _pool(tasks, workers)
    A = pd.DataFrame([r for x in res for r in x["audit"]])
    if len(A):
        _write_csv(A, os.path.join(out_chem, f"audit_{net}.csv"))
    ages = [x["age"] for x in res if x["age"] is not None]
    if ages:
        J = pd.concat(ages, ignore_index=True)
        with _scratch_cwd():
            nom, lb = nominal_bins(net)
        J = J.assign(truth=J.truth + "|" + J.model)
        _write_csv(error_by_age(J, nom, 10, lb), os.path.join(out_chem, f"error_by_age_2ra_{net}.csv"))
    for kind in ("lowkb", "full2r"):
        fr = [x[kind] for x in res if x[kind] is not None]
        if fr:
            _write_csv(pd.concat(fr, ignore_index=True), os.path.join(out_2ra, f"results_{kind}_{net}.csv"))
    info = pd.DataFrame([x["truth_info"] for x in res])
    secs = {f"{x['truth_info']['truth']}/{x['truth_info']['seed']}": round(x["seconds"], 1) for x in res}
    print(f"{net}: {len(tasks)} tasks in {time.time() - t0:.0f} s on {workers} workers", flush=True)
    return {"truth_info": info, "seconds": secs}


NEAR_BINS = ((0.0, 0.2), (0.2, 0.25), (0.25, 0.3), (0.3, 0.4))     # mg/L, [lo, hi)


def _near_task(args) -> dict:
    """One (network, seed, truth): the truth's daily minimum at every junction counted by band around the threshold."""
    net, seed, truth = args
    sc = scenario(net, seed, truth)
    t = sc.truth_daily_min.loc[sc.junctions].values
    row = {"net": net, "truth": truth, "seed": int(seed), "n_junctions": int(len(t)), "mean_daily_min": float(t.mean())}
    for lo, hi in NEAR_BINS:
        row[f"n_{lo:.2f}_to_{hi:.2f}"] = int(((t >= lo) & (t < hi)).sum())
    return row


def near_threshold(net: str, seeds, workers: int = 4, out_2ra: str = OUT_2RA, smoke_ok: bool = False) -> pd.DataFrame:
    """Added after the review, no bar: how many junction-days each truth puts just above the 0.2 mg/L threshold, where
    a map biased low raises false alarms.  Rebuilds every truth of the network on the audit's seeds (the same
    deterministic truths the scored runs used) and writes <out_2ra>/near_threshold_<net>.csv: per truth and seed, the
    count of junctions whose true daily minimum lies in each band of NEAR_BINS (every junction, not only the unsampled
    ones a scored cell sees)."""
    check_seeds(seeds, smoke_ok)
    tasks = [(net, int(s), truth) for truth in NET_TRUTHS[net] for s in seeds]
    df = pd.DataFrame(_pool(tasks, workers, _near_task))
    _write_csv(df, os.path.join(out_2ra, f"near_threshold_{net}.csv"))
    return df


def _near_section(net: str, root: str = "outputs") -> dict:
    """near_threshold_<net>.csv summed over seeds per truth, and the seeds with more junctions just above the threshold
    (0.20 to 0.25 mg/L) under the two-reactant truth than under the first-order truth.  Reported, no bar."""
    path = os.path.join(root, "chem_2ra", f"near_threshold_{net}.csv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, float_precision="round_trip")
    cols = [c for c in df.columns if c.startswith("n_") and c != "n_junctions"]
    out = {"about": "added after the review, no bar: every junction's true daily minimum (all junctions, 8 seeds) by "
                    "band in mg/L, [lo, hi); a map biased low flags junction-days just above 0.2",
           "file": f"outputs/chem_2ra/near_threshold_{net}.csv", "by_truth": {}}
    for truth, g in df.groupby("truth", sort=False):
        out["by_truth"][truth] = {**{c: int(g[c].sum()) for c in cols}, "junction_days": int(g.n_junctions.sum())}
    fo = df[df.truth == "first_order"].set_index("seed").sort_index()
    for truth in [t for t in df.truth.unique() if t != "first_order"]:
        tr = df[df.truth == truth].set_index("seed").sort_index()
        seeds = sorted(set(fo.index) & set(tr.index))
        c = "n_0.20_to_0.25"
        out["by_truth"][truth]["seeds_more_0.20_to_0.25_than_first_order"] = int((tr.loc[seeds, c] > fo.loc[seeds, c]).sum())
        out["by_truth"][truth]["seeds_fewer_0.20_to_0.25_than_first_order"] = int((tr.loc[seeds, c] < fo.loc[seeds, c]).sum())
    return out


# ----------------------------------------------------------------------------- summary
def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if math.isnan(float(x)) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def results_time_path(net: str, truth: str, root: str = "outputs") -> str:
    """The results_time CSV of one truth's committed-model run (experiment.py --chemistry)."""
    chem = "first_order" if truth == "first_order" else "2ra"
    return os.path.join(chem_outdir(root, chem, TRUTHS[truth] or NET_CHEM_MATCH_H), f"results_time_{net}.csv")


def load_time(net: str, truth: str, root: str = "outputs") -> pd.DataFrame:
    df = pd.read_csv(results_time_path(net, truth, root), float_precision="round_trip")
    return df.assign(truth=truth)


def cell_values(df: pd.DataFrame, model: str, rule: str, n: int) -> pd.DataFrame:
    """Per-seed rows of one model, rule and n, indexed by seed."""
    g = df[(df.model == model) & (df.strategy == rule) & (df.n == n)].set_index("seed").sort_index()
    return g


def pooled_recall(g: pd.DataFrame) -> float:
    """Found over low junction-days, summed over seeds (seeds with no low junction drop out)."""
    tp = g["tp_min"].sum() if "tp_min" in g else (g.recall_min * g.n_true_viol_min).round().sum()
    low = g.n_true_viol_min.sum()
    return float(tp / low) if low else float("nan")


def _boot_idx(n_seeds: int) -> np.ndarray:
    return np.random.default_rng(BOOT_SEED).integers(0, n_seeds, size=(BOOT_N, n_seeds))


def boot_rmse_ratio(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Paired bootstrap over seeds of mean RMSE(a) / mean RMSE(b) (the committed aggregation: the mean of per-seed
    daily-minimum RMSEs)."""
    seeds = sorted(set(a.index) & set(b.index))
    x, y = a.loc[seeds, "rmse_min"].values, b.loc[seeds, "rmse_min"].values
    idx = _boot_idx(len(seeds))
    r = x[idx].mean(1) / y[idx].mean(1)
    lo, hi = np.quantile(r, [0.05, 0.95])
    return {"ratio": float(x.mean() / y.mean()), "lo90": float(lo), "hi90": float(hi), "n_seeds": len(seeds)}


def boot_recall_ratio(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Paired bootstrap over seeds of pooled recall(a) / pooled recall(b); resamples where either is undefined or b's
    is zero are dropped and counted."""
    seeds = sorted(set(a.index) & set(b.index))
    ta, la = a.loc[seeds, "tp_min"].values.astype(float), a.loc[seeds, "n_true_viol_min"].values.astype(float)
    tb, lb = b.loc[seeds, "tp_min"].values.astype(float), b.loc[seeds, "n_true_viol_min"].values.astype(float)
    point = (ta.sum() / la.sum()) / (tb.sum() / lb.sum()) if la.sum() and lb.sum() and tb.sum() else float("nan")
    idx = _boot_idx(len(seeds))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (ta[idx].sum(1) / la[idx].sum(1)) / (tb[idx].sum(1) / lb[idx].sum(1))
    ok = np.isfinite(r)
    lo, hi = np.quantile(r[ok], [0.05, 0.95]) if ok.any() else (float("nan"), float("nan"))
    return {"ratio": float(point), "lo90": float(lo), "hi90": float(hi), "n_seeds": len(seeds),
            "n_resamples_dropped": int((~ok).sum())}


def materiality_cells(net: str, root: str = "outputs", truth: str = "2ra") -> list[dict]:
    """The pre-registered cells: B0 (simgp24) daily-minimum RMSE (mean of per-seed values) and recall (pooled counts)
    on the 2RA truth against the paired first-order truth, at n = 8 and 15, random and straddle_min.  Material when
    the RMSE ratio is above 1.10 or the recall ratio below 0.90 (point values); 'robust' when the paired bootstrap 90%
    interval lies wholly on the worse side of 1.0."""
    fo, tr = load_time(net, "first_order", root), load_time(net, truth, root)
    cells = []
    for rule in CELL_RULES:
        for n in CELL_NS:
            a, b = cell_values(tr, "simgp24", rule, n), cell_values(fo, "simgp24", rule, n)
            br, bc = boot_rmse_ratio(a, b), boot_recall_ratio(a, b)
            for metric, boot, worse in (("rmse", br, br["ratio"] > MATERIAL_RMSE),
                                        ("recall", bc, bc["ratio"] < MATERIAL_RECALL)):
                if metric == "rmse":
                    va, vb = float(a.rmse_min.mean()), float(b.rmse_min.mean())
                    robust = boot["lo90"] > 1.0
                else:
                    va, vb = pooled_recall(a), pooled_recall(b)
                    robust = boot["hi90"] < 1.0
                cells.append({"network": net, "truth": truth, "rule": rule, "n": n, "metric": metric,
                              "first_order": vb, "richer_truth": va, "ratio": boot["ratio"], "lo90": boot["lo90"],
                              "hi90": boot["hi90"], "material": bool(worse),
                              "weight": ("robust" if robust else "within seed noise") if worse else None,
                              "counts": ({"found_richer": int(a.tp_min.sum()), "low_richer": int(a.n_true_viol_min.sum()),
                                          "found_first_order": int(b.tp_min.sum()),
                                          "low_first_order": int(b.n_true_viol_min.sum())} if metric == "recall" else None),
                              **({"n_resamples_dropped": boot["n_resamples_dropped"]} if metric == "recall" else {}),
                              # added after the review, no bar: the same ratio with per-seed RMSEs pooled as a root
                              # mean square instead of the committed mean (the trigger's sensitivity to aggregation)
                              **({"ratio_root_mean_square_no_bar": _rms_ratio(a, b)} if metric == "rmse" else {})})
    return cells


def _rms_ratio(a: pd.DataFrame, b: pd.DataFrame) -> float:
    seeds = sorted(set(a.index) & set(b.index))
    x, y = a.loc[seeds, "rmse_min"].values, b.loc[seeds, "rmse_min"].values
    return float(np.sqrt((x ** 2).mean()) / np.sqrt((y ** 2).mean()))


def false_alarm_cells(net: str, root: str = "outputs", truth: str = "2ra") -> list[dict]:
    """Added after the scored run, no bar: B0's false alarms (flagged junction-days whose true daily minimum is not low)
    summed over seeds on the richer truth against the paired first-order truth, at the materiality cells, with the
    paired bootstrap over seeds of the ratio of the sums, the seeds with more false alarms under each truth, and the
    pooled precision."""
    fo, tr = load_time(net, "first_order", root), load_time(net, truth, root)
    out = []
    for rule in CELL_RULES:
        for n in CELL_NS:
            a, b = cell_values(tr, "simgp24", rule, n), cell_values(fo, "simgp24", rule, n)
            seeds = sorted(set(a.index) & set(b.index))
            fa, fb = a.loc[seeds, "fp_min"].values.astype(float), b.loc[seeds, "fp_min"].values.astype(float)
            idx = _boot_idx(len(seeds))
            with np.errstate(divide="ignore", invalid="ignore"):
                r = fa[idx].sum(1) / fb[idx].sum(1)
            ok = np.isfinite(r)
            lo, hi = np.quantile(r[ok], [0.05, 0.95]) if ok.any() else (float("nan"), float("nan"))
            pa = a.tp_min.sum() / (a.tp_min.sum() + a.fp_min.sum())
            pb = b.tp_min.sum() / (b.tp_min.sum() + b.fp_min.sum())
            out.append({"network": net, "truth": truth, "rule": rule, "n": n, "false_alarms_first_order": int(fb.sum()),
                        "false_alarms_richer": int(fa.sum()), "ratio": float(fa.sum() / fb.sum()) if fb.sum() else float("nan"),
                        "lo90": float(lo), "hi90": float(hi), "seeds_more_under_richer": int((fa > fb).sum()),
                        "seeds_more_under_first_order": int((fa < fb).sum()), "precision_first_order": float(pb),
                        "precision_richer": float(pa)})
    return out


def table(df: pd.DataFrame, models, rules, ns=(3, 8, 15)) -> dict:
    """README columns from results_time rows with the extra columns: per model, rule and n, RMSE and MAE (means of
    per-seed values, as committed), coverage 50/80/90/95 (means), recall, precision and F1 (pooled counts), false
    alarms (sum), low junction-days (sum), and the committed mean-of-seeds recall."""
    out = {}
    for model in models:
        for rule in rules:
            for n in ns:
                g = df[(df.model == model) & (df.strategy == rule) & (df.n == n)]
                if not len(g):
                    continue
                tp, fp, low = int(g.tp_min.sum()), int(g.fp_min.sum()), int(g.n_true_viol_min.sum())
                rec = tp / low if low else float("nan")
                prec = tp / (tp + fp) if tp + fp else float("nan")
                ent = {"rmse": float(g.rmse_min.mean()), "mae": float(g.mae_min.mean()), "bias": float(g.bias_min.mean()),
                       "recall": rec, "precision": prec,
                       "f1": 2 * prec * rec / (prec + rec) if (rec == rec and prec == prec and prec + rec) else float("nan"),
                       "false_alarms": fp, "found": tp, "low": low, "recall_mean_of_seeds": float(g.recall_min.mean()),
                       "n_seeds": int(g.seed.nunique())}
                for q in (50, 80, 90, 95):
                    c = f"coverage{q}_min"
                    if c in g and g[c].notna().any():
                        ent[f"coverage{q}"] = float(g[c].mean())
                if "map_kb" in g:
                    ent["map_kb_mean"] = float(g.map_kb.mean())
                out.setdefault(model, {}).setdefault(rule, {})[str(n)] = ent
    return out


def map_kb_shares(df: pd.DataFrame, model: str, floor: float) -> dict:
    """Share of the fits (every rule and n) whose most probable bulk rate is the grid's floor, and below 0.10."""
    g = df[df.model == model]
    return {"n_fits": int(len(g)), "share_map_kb_at_floor": float((np.isclose(g.map_kb, floor)).mean()),
            "share_map_kb_below_0.10": float((g.map_kb < 0.10 - 1e-12).mean()), "floor_per_day": floor}


def _paired_decile_difference(dec: pd.DataFrame, key_a: str, key_b: str, decile: int = 10) -> dict:
    """Added after the review, no bar: the paired bootstrap over seeds of the oldest decile's pooled bias under truth a
    minus truth b (the same seeds, demo samples and unsampled junctions, so the same n per seed), 90% interval."""
    a = dec[(dec.truth == key_a) & (dec.seed != "all") & (dec.bin == decile)].set_index("seed").sort_index()
    b = dec[(dec.truth == key_b) & (dec.seed != "all") & (dec.bin == decile)].set_index("seed").sort_index()
    seeds = sorted(set(a.index) & set(b.index))
    wa, wb = a.loc[seeds, "n"].values.astype(float), b.loc[seeds, "n"].values.astype(float)
    ba, bb = a.loc[seeds, "bias"].values, b.loc[seeds, "bias"].values
    idx = _boot_idx(len(seeds))
    d = (ba[idx] * wa[idx]).sum(1) / wa[idx].sum(1) - (bb[idx] * wb[idx]).sum(1) / wb[idx].sum(1)
    lo, hi = np.quantile(d, [0.05, 0.95])
    point = float((ba * wa).sum() / wa.sum() - (bb * wb).sum() / wb.sum())
    return {"difference_lo90": float(lo), "difference_hi90": float(hi), "difference_n_seeds": len(seeds),
            "difference_same_n_per_seed": bool(np.array_equal(wa, wb)),
            "difference_inside_bar": bool(abs(point) <= AGE_BAR_MGL)}


def _age_section(net: str, out_chem: str = OUT_CHEM) -> dict:
    dec = pd.read_csv(os.path.join(out_chem, f"error_by_age_2ra_{net}.csv"), dtype={"seed": str})
    allr = dec[dec.seed == "all"]
    out = {"bins": "deciles of the operator's nominal daily-mean age (task 9's bins; 1 = youngest)", "by_truth": {}}
    for key, g in allr.groupby("truth", sort=False):
        g = g.sort_values("bin")
        per_seed = dec[(dec.truth == key) & (dec.seed != "all") & (dec.bin == 10)]
        b = per_seed.bias.values
        w = per_seed.n.values.astype(float)
        idx = _boot_idx(len(b))
        boot = (b[idx] * w[idx]).sum(1) / w[idx].sum(1)
        lo, hi = np.quantile(boot, [0.05, 0.95])
        out["by_truth"][key] = {"bias_by_decile": [round(float(v), 4) for v in g.bias], "rmse_by_decile":
                                [round(float(v), 4) for v in g.rmse], "n_by_decile": [int(v) for v in g.n],
                                "oldest_decile": {"bias": float(g.bias.iloc[-1]), "rmse": float(g.rmse.iloc[-1]),
                                                  "n": int(g.n.iloc[-1]), "bias_lo90": float(lo), "bias_hi90": float(hi),
                                                  "nominal_age_h": [float(g.nominal_age_lo_h.iloc[-1]),
                                                                    float(g.nominal_age_hi_h.iloc[-1])],
                                                  "frac_lower_bound": float(g.frac_lower_bound.iloc[-1])}}
    base = os.path.join(out_chem, f"error_by_age_{net}.csv")
    if os.path.exists(base):
        t9 = pd.read_csv(base, dtype={"seed": str})
        t9 = t9[(t9.truth == "default") & (t9.seed == "all")].sort_values("bin")
        out["task9_committed_baseline"] = {"file": f"outputs/chem/error_by_age_{net}.csv (truth default, seeds 0 to 7)",
                                           "bias_by_decile": [round(float(v), 4) for v in t9.bias],
                                           "oldest_decile_bias": float(t9.bias.iloc[-1])}
    two = out["by_truth"].get("2ra|B0", {}).get("oldest_decile", {})
    fo = out["by_truth"].get("first_order|B0", {}).get("oldest_decile", {})
    if two:
        diff_ci = _paired_decile_difference(dec, "2ra|B0", "first_order|B0") if fo else {}
        out["rule"] = {"bar": f"2RA truth, B0, 8 demo samples, pooled oldest decile: |bias| > {AGE_BAR_MGL} mg/L is recorded "
                              "as the first-order limit and the next step is named, not built",
                       "bias_2ra": two["bias"], "bias_first_order_paired": fo.get("bias"),
                       "difference_2ra_minus_first_order": two["bias"] - fo["bias"] if fo else None,
                       **diff_ci,
                       "first_order_truth_also_past_bar": (bool(abs(fo["bias"]) > AGE_BAR_MGL
                                                                and np.sign(fo["bias"]) == np.sign(two["bias"]))
                                                           if fo else None),
                       "triggered": bool(abs(two["bias"]) > AGE_BAR_MGL),
                       "expectation_negative": bool(two["bias"] < 0),
                       "next_step_named": ("an EPANET order-2 grid with a simulated dose axis (about 5 times the runs)"
                                           if abs(two["bias"]) > AGE_BAR_MGL else None)}
    return out


def _sum_rates(g: pd.DataFrame) -> dict:
    return rates({c: float(g[c].sum()) for c in SUM_COLS})


def _oracle_section(A: pd.DataFrame) -> dict:
    out = {}
    for (truth, grid), g in A[A.part == "oracle"].groupby(["truth", "model"], sort=False):
        out.setdefault(truth, {})[grid] = {**_sum_rates(g), "rms_ln_day_mean": float(g.rms_ln_day.mean()),
                                           "rms_ln_night_mean": float(g.rms_ln_night.mean()),
                                           "member_kb_by_seed": [float(v) for v in g.member_kb],
                                           "member_dose_by_seed": [float(v) for v in g.member_dose]}
    return out


def _dose_section(A: pd.DataFrame) -> dict:
    D = A[A.part == "dose_step"]
    if not len(D):
        return {}
    out = {"protocol": f"B0 fitted on {DOSE_N} demo samples at 1.2 mg/L; its daily-minimum map times the dose ratio "
                       "(the exact ln-offset of first order), P(daily min < 0.2) from its own draws at 0.2 / ratio; "
                       "scored on the unsampled junctions against the truth rerun at that dose, every draw the same",
           "by_truth": {}, "paired_bias_difference": {}}
    for truth, g in D.groupby("truth", sort=False):
        for r, h in g.groupby("dose_ratio"):
            out["by_truth"].setdefault(truth, {})[f"{float(r):.4g}"] = {
                **_sum_rates(h), "dose_mgL": float(h.dose_mgL.iloc[0]),
                "truth_ln_shift_mean": float(h.truth_ln_shift_mean.mean()),
                "truth_ln_shift_range": [float(h.truth_ln_shift_min.min()), float(h.truth_ln_shift_max.max())],
                "first_order_ln_shift": float(h.first_order_ln_shift.iloc[0])}
    if {"first_order", "2ra"} <= set(D.truth):
        for r, h in D.groupby("dose_ratio"):
            a = h[h.truth == "2ra"].set_index("seed").sort_index()
            b = h[h.truth == "first_order"].set_index("seed").sort_index()
            seeds = sorted(set(a.index) & set(b.index))
            da = (a.loc[seeds, "se"] / a.loc[seeds, "n"]).values
            db = (b.loc[seeds, "se"] / b.loc[seeds, "n"]).values
            d = da - db
            idx = _boot_idx(len(seeds))
            lo, hi = np.quantile(d[idx].mean(1), [0.05, 0.95])
            out["paired_bias_difference"][f"{float(r):.4g}"] = {"mean_bias_2ra_minus_first_order_mgL": float(d.mean()),
                                                                "lo90": float(lo), "hi90": float(hi)}
        up = out["paired_bias_difference"].get(f"{DOSE_RATIOS[1]:.4g}", {})
        dn = out["paired_bias_difference"].get(f"{DOSE_RATIOS[0]:.4g}", {})
        out["expectation"] = {"pre_registered": "under 2RA, relative to the first-order twin, the map moved up a dose "
                                                "is conservative (biased low) and moved down a dose optimistic (biased high)",
                              "dose_up_conservative": bool(up.get("mean_bias_2ra_minus_first_order_mgL", 0) < 0),
                              "dose_down_optimistic": bool(dn.get("mean_bias_2ra_minus_first_order_mgL", 0) > 0)}
    return out


def _variant_section(net: str, kind: str, root: str = "outputs") -> dict:
    """An opt-in variant (lowkb, or the gated full2r) against B0 on the same truth: per rule and n (8, 15), RMSE ratio
    and pooled recall ratio with paired bootstrap intervals, coverage, false alarms, and its MAP bulk-rate shares."""
    path = os.path.join(root, "chem_2ra", f"results_{kind}_{net}.csv")
    if not os.path.exists(path):
        return {}
    V = pd.read_csv(path, float_precision="round_trip")
    out = {"by_truth": {}}
    for truth in NET_TRUTHS[net]:
        v = V[V.truth == truth]
        if not len(v):
            continue
        b = load_time(net, truth, root)
        ent = {"cells": {}, "table": table(v, [f"simgp24_{kind}"], ("random", "uncertainty", "straddle", "straddle_min")),
               "table_B0": table(b, ["simgp24"], ("random", "straddle_min"))}
        for rule in CELL_RULES:
            for n in CELL_NS:
                a, c = cell_values(v, f"simgp24_{kind}", rule, n), cell_values(b, "simgp24", rule, n)
                ent["cells"][f"{rule}_{n}"] = {"rmse": boot_rmse_ratio(a, c), "recall": boot_recall_ratio(a, c),
                                               "rmse_variant": float(a.rmse_min.mean()), "rmse_B0": float(c.rmse_min.mean()),
                                               "recall_variant": pooled_recall(a), "recall_B0": pooled_recall(c),
                                               "false_alarms_variant": int(a.fp_min.sum()), "false_alarms_B0": int(c.fp_min.sum())}
        cov = v[(v.strategy == "random") & v.n.isin((3, 8, 15))].groupby("n").coverage90_min.mean()
        covb = b[(b.model == "simgp24") & (b.strategy == "random") & b.n.isin((3, 8, 15))].groupby("n").coverage90_min.mean()
        ent["coverage90_random_mean_over_n"] = float(cov.mean())
        ent["coverage90_random_by_n"] = {str(int(k)): float(x) for k, x in cov.items()}
        ent["coverage90_random_mean_over_n_B0"] = float(covb.mean())
        if kind == "lowkb":
            ent["map_kb_variant"] = map_kb_shares(v, f"simgp24_{kind}", min(KB_LOWKB))
        ent["map_kb_B0"] = map_kb_shares(b, "simgp24", min(KB_GRID))
        out["by_truth"][truth] = ent
    return out


def lowkb_recommendation(net_sections: dict) -> dict:
    """The pre-registered reading of the opt-in low-kb grid (a recommendation for Devansh, never applied here)."""
    r1, r2, r3, notes = True, True, False, []
    for net, sec in net_sections.items():
        for truth, ent in sec.get("by_truth", {}).items():
            for key, c in ent["cells"].items():
                if truth == "first_order":
                    if c["rmse"]["ratio"] > MATERIAL_RMSE or (c["recall"]["ratio"] == c["recall"]["ratio"]
                                                              and c["recall"]["ratio"] < MATERIAL_RECALL):
                        r1 = False
                        notes.append(f"{net} first_order {key}: RMSE ratio {c['rmse']['ratio']:.3f}, recall ratio "
                                     f"{c['recall']['ratio']:.3f}")
                if truth == "2ra" and c["rmse"]["hi90"] < 1.0:
                    r3 = True
            cv = ent["coverage90_random_mean_over_n"]
            if not (COVERAGE_BAND[0] <= cv <= COVERAGE_BAND[1]):
                r2 = False
                notes.append(f"{net} {truth}: random-rule coverage90 mean over n {cv:.3f}")
    ok = r1 and r2 and r3
    return {"R1_no_loss_over_10pct_on_first_order_truths": r1, "R2_coverage90_in_band_every_truth": r2,
            "R3_measurable_rmse_gain_on_a_2ra_cell": r3, "candidate_for_default": ok, "notes": notes,
            "reading": ("a candidate for the default grid; adopting it would move committed numbers, which is Devansh's "
                        "decision" if ok else "tested, not recommended for the default grid"),
            "default_unchanged": True}


def committed_context(net: str, root: str = "outputs") -> dict:
    """The committed first-order outputs on seeds 0 to 7 (outputs/results_time_<net>.csv), B0 at the materiality cells,
    for context only: those seeds are not the audit's (the audit pairs on fresh seeds).  Recall here is pooled from
    recall x low junction-days (the committed file has no counts)."""
    path = os.path.join(root, f"results_time_{net}.csv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, float_precision="round_trip")
    out = {"file": f"outputs/results_time_{net}.csv", "seeds": sorted(int(x) for x in df.seed.unique()), "cells": {}}
    for rule in CELL_RULES:
        for n in CELL_NS:
            g = cell_values(df, "simgp24", rule, n)
            low = g.n_true_viol_min.sum()
            out["cells"][f"{rule}_{n}"] = {"rmse": float(g.rmse_min.mean()),
                                           "recall_pooled": float((g.recall_min * g.n_true_viol_min).round().sum() / low)
                                           if low else float("nan"), "recall_mean_of_seeds": float(g.recall_min.mean())}
    return out


def snapshot_context(net: str, root: str = "outputs") -> dict:
    """The 14:00 snapshot experiment (results_<net>.csv) and the routes (results_routes_<net>.csv) of each truth's run,
    the calibrated-simulator GP at n = 3, 8, 15 and the optimised route at K = 5, 8, 12; reported, no bar."""
    out = {}
    for truth in NET_TRUTHS[net]:
        d = os.path.dirname(results_time_path(net, truth, root))
        snap = pd.read_csv(os.path.join(d, f"results_{net}.csv"))
        rt = pd.read_csv(os.path.join(d, f"results_routes_{net}.csv"))
        ent = {}
        for rule in ("random", "straddle"):
            for n in (3, 8, 15):
                g = snap[(snap.model == "simgp_core") & (snap.strategy == rule) & (snap.n == n)]
                low = g.n_true_viol.sum()
                ent[f"simgp_core_{rule}_{n}"] = {"rmse": float(g.rmse.mean()), "coverage90": float(g.coverage90.mean()),
                                                 "recall_pooled": float((g.recall * g.n_true_viol).round().sum() / low)
                                                 if low else float("nan")}
        for K in (5, 8, 12):
            g = rt[(rt.route == "optimised") & (rt.K == K)]
            ent[f"route_optimised_K{K}"] = {"recall_min": float(g.recall_min.mean()), "rmse_min": float(g.rmse_min.mean())}
        out[truth] = ent
    return out


def gate_section(cells: list[dict], root: str = "outputs") -> dict:
    """The gated two-rate model's adoption bars, per network where a cell is material (pre-registered)."""
    out = {}
    for net in sorted({c["network"] for c in cells if c["material"]}):
        sec = _variant_section(net, "full2r", root)
        if not sec:
            out[net] = {"built": False, "note": "material, full2r not run yet"}
            continue
        two, fo = sec["by_truth"].get("2ra"), sec["by_truth"].get("first_order")
        g1 = []
        for c in [c for c in cells if c["network"] == net and c["material"]]:
            key = f"{c['rule']}_{c['n']}"
            ent = two["cells"][key]
            b0_2ra, b0_fo = c["richer_truth"], c["first_order"]
            if c["metric"] == "rmse":
                target = b0_2ra - 0.5 * (b0_2ra - b0_fo)
                got = ent["rmse_variant"]
                ok = got <= target
            else:
                target = b0_2ra + 0.5 * (b0_fo - b0_2ra)
                got = ent["recall_variant"]
                ok = got >= target
            g1.append({"cell": f"{key}_{c['metric']}", "B0_richer": b0_2ra, "B0_first_order": b0_fo, "target": target,
                       "full2r": got, "recovers_half": bool(ok)})
        g2 = {t: sec["by_truth"][t]["coverage90_random_mean_over_n"] for t in ("2ra", "first_order") if t in sec["by_truth"]}
        g3 = [{"cell": k, "rmse_ratio": v["rmse"]["ratio"], "recall_ratio": v["recall"]["ratio"],
               "ok": bool(v["rmse"]["ratio"] <= MATERIAL_RMSE and not (v["recall"]["ratio"] < MATERIAL_RECALL))}
              for k, v in fo["cells"].items()]
        passed = (all(x["recovers_half"] for x in g1) and all(COVERAGE_BAND[0] <= v <= COVERAGE_BAND[1] for v in g2.values())
                  and all(x["ok"] for x in g3))
        out[net] = {"built": True, "G1_recovers_half_of_each_material_gap": g1, "G2_coverage90_random_mean_over_n": g2,
                    "G2_band": list(COVERAGE_BAND), "G3_first_order_truth_loss_at_most_10pct": g3, "passed": bool(passed),
                    "reading": ("passes its bars: a candidate for the default model, which would move committed numbers "
                                "(Devansh's decision)" if passed else "tried and rejected")}
    return out


def summarise(net: str, root: str = "outputs") -> dict:
    """outputs/chem/summary_audit_<net>.json from the committed CSVs as written (both networks' results_time files feed
    the materiality gate, which spans Net3 and Net2)."""
    out_chem = os.path.join(root, "chem")
    A = pd.read_csv(os.path.join(out_chem, f"audit_{net}.csv"), float_precision="round_trip")
    truths = {t: load_time(net, t, root) for t in NET_TRUTHS[net]}
    cells = []
    for n2 in ("Net3", "Net2"):
        if all(os.path.exists(results_time_path(n2, t, root)) for t in PRIMARY):
            cells += materiality_cells(n2, root)
    material = any(c["material"] for c in cells)
    sens = materiality_cells("Net2", root, "2ra_m96") if (net == "Net2" and os.path.exists(results_time_path("Net2", "2ra_m96", root))) else []
    models = ("simgp24", "simgp_timeblind", "mean_of_samples", "simgp24_g75")
    S = {"generated_by": f"python -m residualmap.chemexp summarise {net}", "network": net,
         "about": "Task 13, the richer-truth audit. Simulation only. Today's model (B0) against Fisher's two-reactant "
                  "chlorine-organics chemistry in EPANET-MSX (2RA) and the committed first-order truth, paired on the same "
                  "fresh seeds and draws. Every number is computed from the CSVs named in 'sources' as written.",
         "sources": {"results_time": {t: results_time_path(net, t, root) for t in NET_TRUTHS[net]},
                     "audit": f"outputs/chem/audit_{net}.csv", "age": f"outputs/chem/error_by_age_2ra_{net}.csv",
                     "lowkb": f"outputs/chem_2ra/results_lowkb_{net}.csv", "scaling": "outputs/chem_2ra/scaling.json"},
         "seeds": sorted(int(s) for s in truths["first_order"].seed.unique()),
         "definitions": {
             "B0": "today's model: SimGP24 on the committed 675-member grid, committed settings (simgp24 in results_time)",
             "rmse, mae, coverage": "means over seeds of per-seed daily-minimum values at unsampled junctions (as committed)",
             "recall, precision, F1, false_alarms": "pooled counts over seeds (found of low junction-days; a seed with no "
                                                   "low junction adds nothing); recall_mean_of_seeds is the committed "
                                                   "convention (1.0 for a seed with nothing to find)",
             "oracle": "best single member and dose against the full daytime truth, daily minimum at every junction",
             "lowkb": "B0 on simgp.GRIDS['full_lowkb'] (bulk rates 0.0125 to 0.70 per day), an opt-in variant"},
         "readme_columns": {t: table(df, models, ("random", "straddle_min")) for t, df in truths.items()},
         "map_kb_B0": {t: map_kb_shares(df, "simgp24", min(KB_GRID)) for t, df in truths.items()},
         "oracle": _oracle_section(A), "age_deciles": _age_section(net, out_chem), "dose_step": _dose_section(A),
         "lowkb": _variant_section(net, "lowkb", root),
         "false_alarms_no_bar": false_alarm_cells(net, root) + (false_alarm_cells(net, root, "2ra_m96") if net == "Net2" else []),
         "near_threshold_no_bar": _near_section(net, root),
         "snapshot_and_routes_no_bar": snapshot_context(net, root),
         "committed_seeds_0_to_7_context": committed_context(net, root)}
    S["acceptance"] = {
        "materiality": {"bar": "material if, at n = 8 or 15 under the random or straddle_min rule, on Net3 or Net2 (24 h "
                               "match), B0's 2RA-truth daily-minimum RMSE is above 1.10 x, or its pooled recall below 0.90 x, "
                               "the paired first-order truth's (point ratios; intervals reported)",
                        "cells": cells, "material": bool(material),
                        "triggered": [f"{c['network']} {c['rule']} n={c['n']} {c['metric']}: {c['ratio']:.3f} "
                                      f"({c['lo90']:.3f} to {c['hi90']:.3f}), {c['weight']}" for c in cells if c["material"]],
                        "reading": ("material: the gated two-rate model is built" if material else
                                    "first order costs at most 10% at fixed conditions")},
        "net2_match96_sensitivity_no_bar": sens,
        "age_decile": S["age_deciles"].get("rule"),
        "dose_step_expectation_no_bar": S["dose_step"].get("expectation"),
        "lowkb_recommendation": lowkb_recommendation({n2: _variant_section(n2, "lowkb", root) for n2 in ("Net3", "Net2")
                                                      if os.path.exists(os.path.join(root, "chem_2ra", f"results_lowkb_{n2}.csv"))}),
        "checks": "outputs/chem/checks_report.json (the task-13 group) and outputs/chem/baseline_reproduction.json"}
    if material:
        S["acceptance"]["full2r"] = gate_section(cells, root)
    return S


def write_summary(net: str, root: str = "outputs") -> dict:
    S = summarise(net, root)
    path = os.path.join(root, "chem", f"summary_audit_{net}.json")
    with open(path, "w") as fh:
        json.dump(_clean(S), fh, indent=1)
        fh.write("\n")
    print(f"wrote {path}: material = {S['acceptance']['materiality']['material']}", flush=True)
    return S


# ----------------------------------------------------------------------------- figures
BIG = {"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 11, "xtick.labelsize": 11,
       "ytick.labelsize": 11}


def plot_cost(net: str, S: dict, root: str = "outputs", out: str | None = None) -> None:
    """Daily-minimum RMSE, recall (pooled), false alarms (summed) and coverage90 against n for B0 on the first-order and
    2RA truths (paired seeds), random and straddle_min, with the opt-in low-kb grid on the 2RA truth (random rule)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = out or os.path.join(root, "chem", f"cost_of_first_order_{net}.png")
    lv = os.path.join(root, "chem_2ra", f"results_lowkb_{net}.csv")
    V = pd.read_csv(lv) if os.path.exists(lv) else None
    style = {"first_order": ("tab:blue", "first-order truth"), "2ra": ("tab:red", "two-reactant truth (2RA)"),
             "2ra_m96": ("tab:orange", "2RA, matched at 96 h")}

    def rec(g):
        return g.apply(lambda x: x.tp_min.sum() / x.n_true_viol_min.sum() if x.n_true_viol_min.sum() else np.nan)

    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 4, figsize=(27, 6.4))
        for truth in NET_TRUTHS[net]:
            df = load_time(net, truth, root)
            c, lab = style[truth]
            for rule, ls in (("random", "--"), ("straddle_min", "-")):
                g = df[(df.model == "simgp24") & (df.strategy == rule)].groupby("n")
                axes[0].plot(g.rmse_min.mean().index, g.rmse_min.mean().values, ls, color=c, marker="o", ms=4,
                             label=f"today's model, {lab}, {rule}")
                r = rec(g)
                axes[1].plot(r.index, r.values, ls, color=c, marker="o", ms=4)
                axes[2].plot(g.fp_min.sum().index, g.fp_min.sum().values, ls, color=c, marker="o", ms=4)
                axes[3].plot(g.coverage90_min.mean().index, g.coverage90_min.mean().values, ls, color=c, marker="o", ms=4)
            if V is not None and truth == "2ra":
                v = V[(V.truth == truth) & (V.strategy == "random")].groupby("n")
                kw = dict(color="tab:purple", marker="s", ms=4)
                axes[0].plot(v.rmse_min.mean().index, v.rmse_min.mean().values, ":", label="low-kb grid (opt-in), 2RA truth, random", **kw)
                r = rec(v)
                axes[1].plot(r.index, r.values, ":", **kw)
                axes[2].plot(v.fp_min.sum().index, v.fp_min.sum().values, ":", **kw)
                axes[3].plot(v.coverage90_min.mean().index, v.coverage90_min.mean().values, ":", **kw)
        axes[0].set(title="Daily-minimum error at unsampled junctions", xlabel="daytime grab samples", ylabel="RMSE (mg/L)")
        axes[1].set(title=f"Recall: low junction-days found\n(daily minimum < {THRESHOLD} mg/L, pooled over seeds)",
                    xlabel="daytime grab samples", ylabel="recall", ylim=(0.8, 1.01))
        axes[2].set(title="False alarms (summed over the 8 seeds)", xlabel="daytime grab samples", ylabel="flagged, not low")
        axes[3].set(title="Truth inside the 90% band (daily minimum)", xlabel="daytime grab samples", ylabel="coverage",
                    ylim=(0.7, 1.01))
        axes[3].axhline(0.9, color="k", lw=0.8, ls=":")
        for ax in axes:
            ax.grid(alpha=0.3)
        axes[0].legend(fontsize=9)
        m = S["acceptance"]["materiality"]
        own = [c for c in m["cells"] if c["network"] == net]
        n_mat = sum(c["material"] for c in own)
        state = (f"{n_mat} of {len(own)} pre-registered cells on {net} more than 10% worse" if n_mat else
                 f"no pre-registered cell on {net} more than 10% worse")
        fig.suptitle(f"{net}, seeds {min(S['seeds'])} to {max(S['seeds'])}: today's model against two-reactant chemistry "
                     f"(simulated, same seeds and draws); {state}", fontsize=15)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_age(net: str, S: dict, root: str = "outputs", out: str | None = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = out or os.path.join(root, "chem", f"error_by_age_2ra_{net}.png")
    a = S["age_deciles"]
    style = {"first_order|B0": ("tab:blue", "-", "B0, first-order truth (same seeds)"),
             "2ra|B0": ("tab:red", "-", "B0, two-reactant truth"),
             "2ra|lowkb": ("tab:purple", ":", "low-kb grid (opt-in), two-reactant truth"),
             "first_order|lowkb": ("tab:cyan", ":", "low-kb grid (opt-in), first-order truth"),
             "2ra_m96|B0": ("tab:orange", "--", "B0, 2RA matched at 96 h")}
    x = np.arange(1, 11)
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(figsize=(11, 6.4))
        for key, ent in a["by_truth"].items():
            if key in style:
                c, ls, lab = style[key]
                ax.plot(x, ent["bias_by_decile"], ls, color=c, marker="o", ms=5, lw=2, label=lab)
        if "task9_committed_baseline" in a:
            ax.plot(x, a["task9_committed_baseline"]["bias_by_decile"], "-", color="gray", lw=1.2, marker=".",
                    label="task 9's committed baseline (first-order truth, seeds 0 to 7)")
        ax.axhspan(-AGE_BAR_MGL, AGE_BAR_MGL, color="green", alpha=0.08)
        ax.axhline(0, color="k", lw=0.8)
        r = a.get("rule", {})
        ax.set(xticks=x, xlabel="nominal water-age decile (1 = youngest)", ylabel="bias of the daily-minimum map (mg/L)",
               title=f"{net}: daily-minimum bias by water age after 8 samples (simulated)\noldest decile under 2RA "
                     f"{r.get('bias_2ra', float('nan')):+.3f} mg/L (first-order truth {r.get('bias_first_order_paired', float('nan')):+.3f}); "
                     f"pre-registered bar {AGE_BAR_MGL} mg/L (shaded)")
        ax.grid(alpha=0.3); ax.legend(fontsize=9)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_dose(net: str, S: dict, root: str = "outputs", out: str | None = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = S["dose_step"]
    if not d:
        return
    out = out or os.path.join(root, "chem", f"dose_step_{net}.png")
    style = {"first_order": ("tab:blue", "first-order truth (the offset is exact)"), "2ra": ("tab:red", "two-reactant truth")}
    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 3, figsize=(20, 5.8))
        for truth, ent in d["by_truth"].items():
            c, lab = style[truth]
            ks = sorted(ent, key=float)
            doses = [ent[k]["dose_mgL"] for k in ks]
            axes[0].plot(doses, [ent[k]["bias"] for k in ks], "-o", color=c, label=lab)
            axes[1].plot(doses, [ent[k]["rmse"] for k in ks], "-o", color=c, label=lab)
            axes[2].plot(doses, [ent[k]["truth_ln_shift_mean"] for k in ks], "-o", color=c, label=f"truth, {lab}")
        axes[2].plot([0.9, 1.2, 1.6], [math.log(0.75), 0.0, math.log(4 / 3)], "k:", label="first order: ln of the dose ratio")
        axes[0].axhline(0, color="k", lw=0.8)
        axes[0].set(xlabel="plant dose (mg/L)", ylabel="bias of the moved map (mg/L)",
                    title="Map fitted at 1.2 mg/L, moved to another dose")
        axes[1].set(xlabel="plant dose (mg/L)", ylabel="RMSE (mg/L)", title="Daily-minimum error at the new dose")
        axes[2].set(xlabel="plant dose (mg/L)", ylabel="mean ln(daily min / daily min at 1.2)",
                    title="How the truth's daily minimum moves with the dose")
        for ax in axes:
            ax.grid(alpha=0.3); ax.legend(fontsize=9)
        pb = d.get("paired_bias_difference", {})
        up, dn = pb.get(f"{DOSE_RATIOS[1]:.4g}", {}), pb.get(f"{DOSE_RATIOS[0]:.4g}", {})
        fig.suptitle(f"{net}: the dose step, 15 samples at 1.2 mg/L (simulated).\nAgainst the first-order twin, under 2RA the map "
                     f"moved to 1.6 mg/L is {abs(up.get('mean_bias_2ra_minus_first_order_mgL', float('nan'))):.3f} mg/L "
                     f"{'lower (conservative)' if up.get('mean_bias_2ra_minus_first_order_mgL', 0) < 0 else 'higher'}, "
                     f"moved to 0.9 mg/L {abs(dn.get('mean_bias_2ra_minus_first_order_mgL', float('nan'))):.3f} mg/L "
                     f"{'higher (optimistic)' if dn.get('mean_bias_2ra_minus_first_order_mgL', 0) > 0 else 'lower'}",
                     fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------- CLI
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("scaling", "audit", "summarise", "replot", "full2r", "near"))
    ap.add_argument("net", nargs="?", default="Net3")
    ap.add_argument("n", nargs="?", type=int, default=N_SEEDS)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=CHEM_FIRST_SEED)
    ap.add_argument("--root", default="outputs", help="outputs root (a scratch folder for a smoke run)")
    ap.add_argument("--parts", default="oracle,age,dose,lowkb")
    a = ap.parse_args(argv)
    if a.command == "scaling":
        S = scaling()
        os.makedirs(os.path.join(a.root, "chem_2ra"), exist_ok=True)
        with open(os.path.join(a.root, "chem_2ra", "scaling.json"), "w") as fh:
            json.dump(_clean(S), fh, indent=1)
            fh.write("\n")
        print(json.dumps(_clean(S)["networks"], indent=1))
        return 0
    smoke = a.root != "outputs"
    if a.seed0 != CHEM_FIRST_SEED and not smoke:
        raise SystemExit("seeds other than the audit's go to a scratch --root")
    seeds = tuple(range(a.seed0, a.seed0 + a.n))
    if a.command in ("audit", "full2r"):
        parts = ("full2r",) if a.command == "full2r" else tuple(p for p in a.parts.split(",") if p)
        res = audit(a.net, seeds, a.workers, cache_dir=os.path.join("outputs", "cache"),
                    out_chem=os.path.join(a.root, "chem"), out_2ra=os.path.join(a.root, "chem_2ra"), parts=parts,
                    smoke_ok=smoke)
        print(res["truth_info"].to_string())
        return 0
    if a.command == "near":
        df = near_threshold(a.net, seeds, a.workers, out_2ra=os.path.join(a.root, "chem_2ra"), smoke_ok=smoke)
        print(df.groupby("truth", sort=False).sum(numeric_only=True).drop(columns="seed").to_string())
        return 0
    if a.command in ("summarise", "replot"):
        S = write_summary(a.net, a.root) if a.command == "summarise" else json.load(
            open(os.path.join(a.root, "chem", f"summary_audit_{a.net}.json")))
        plot_cost(a.net, S, a.root)
        plot_age(a.net, S, a.root)
        plot_dose(a.net, S, a.root)
        return 0
    return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
