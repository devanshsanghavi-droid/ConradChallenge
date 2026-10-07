"""
organics.py: plant TOC as a logged input (iteration 4, journal task 11).

Organic matter in the water consumes chlorine.  Today's model (simgp.SimGP24) has no organics input: their effect is
lumped into the bulk decay rate that is refitted from the last three months of grab samples, so a model fitted in a dry
late summer has no way to know that the first winter storms bring more organics.  Here the plant's monthly TOC (a
plant-log column, toc_mgL) becomes an input that scales the bulk decay rate.

THE MODEL, M_TOC (TocSimGP24, a SeasonalSimGP24 subclass, so SimGP24 and the temperature models are untouched)
  Hypotheses, stacked along the grid's member axis:
    H0     today's 675-member grid, the same every month (TOC has no effect);
    H_TOC  the same grid re-run with every bulk rate times TOC_m / TOC_ref, TOC_ref = 2.0 mg/L, the AWWARF / Powell
           form k = alpha TOC exp(-(E/R)/T) at a fixed temperature (AWWARF 1996, WRF project 815; Powell et al. 2000,
           Water Res 34:117, doi:10.1016/S0043-1354(99)00097-4, via Liu, Reckhow & Li 2014, doi:10.1016/j.watres.
           2014.01.010; refit on another water by Saidan et al. 2017, doi:10.5277/epe170417).  alpha is not taken from
           the literature (it varies about 6 times across waters): it is absorbed into the calibrated kb, read as the
           rate at TOC_ref.  The C0 exponent of the refits is set to 0 on purpose, so the dose stays an exact ln-offset
           (a stated limitation).  The wall rate is not scaled.  Powell's validity range is TOC 1 to 3 mg/L.
  At TOC_ref the H_TOC block IS the committed grid, so with TOC 2.0 every month M_TOC is today's model (a saved check).
  Prior: H0 1/4, H_TOC 3/4 (H0 keeps the 1/4 it had in tasks 10 and 10b); sensitivity priors 1/2 and 0 on H0 are
  reported and never used to choose.  Likelihood, discrepancy GP, daily-minimum Monte Carlo and the 3-month window
  are SeasonalSimGP24's (SimGP24's).  Temperature is not in M_TOC: tasks 10 and 10b were not adopted (journal).

THE TRUTHS (pilot.synthetic_log with a schedule, simulate.build_scenario(chem=...)), all at a fixed 20 C so temperature
cannot confound the TOC test (the truth's E/R and wall-temperature draws are made and multiply by exactly 1):
  O1  matched structure, a debug check: first order, kb_m = kb_net u_m TOC_m / 2 (u_m the committed U(0.8, 1.2) draw);
  O2  Clark (Clark 1998, J Environ Eng 124(1):16), EPANET order 2 with a limiting potential:
      dC/dt = -k2 u_m C (C - CL_m), CL_m = d_m - phi TOC_m, phi = 0.85 mg Cl2 per mg C, k2 = kb_net / (phi TOC_ref)
      (0.235 L/mg/day on Net3, 0.059 on Net2): decay slows as the water ages, the residual responds nonlinearly to the
      dose, and the TOC effect on aged water is steeper than linear;
  O3  low demand (like membrane-filtered water): phi = 0.5, k2 = kb_net / 1.0; CL = +0.2 mg/L at TOC 2.0 and dose
      1.2, so the bulk demand runs out.
  phi, k2 and the TOC schedule are ASSUMPTIONS.  The wall stays first order as committed.  Schedules: S, the assumed
  monthly TOC (chemistry.monthly_toc: 2.5 mg/L January to March, 2.0 April, 1.5 May to October, 3.0 November, the
  first storms, 2.5 December); C, a control at 2.0 every month; D (O2 on Net3, no bar), S with the plant dose stepped
  from 1.2 to 1.5 mg/L in July, the logged dose entering as an exact ln-offset.

THE EXPERIMENT (python -m residualmap.organics <net> <n_seeds>; journal task 11, bars pre-registered there).  12-month
logs as in task 10 (10 route taps and 3 rotating taps a month, hours 07:00 to 17:00, N(0, 0.03) mg/L noise), fresh
seeds (Net3 from 32, Net2 from 16; the seeds of tasks 10 and 10b are refused).  Models: B0 (today's SimGP24, previous 3
months, TOC-blind), M_TOC, the two sensitivity priors (S), M_TOC told no dose (D), persistence and the network mean
(readings, S) and the oracle (the best single member and dose of M_TOC's bank for the target month against the full
daytime truth: the floor for any model built on this first-order grid).  Tests: rolling (fit the previous 3 months,
predict April to December) in four classes by the TOC's course: falling (April to July), steady (August to October),
first_storm (November: TOC doubles from 1.5 to 3.0 after an August-to-October fit) and after_storm (December); and an
extrapolation without refit (fit January to March at 2.5 mg/L, predict July, August and September at 1.5).  Scores as
task 10 (held-out readings at their own junction and hour; the daily-minimum map at every junction not sampled in the
3-month window), plus the map's RMSE and bias by nominal-age tercile.  Outputs: outputs/chem/organics_<net>.csv (the
sums, per truth, schedule, test, seed, class and model, 6 significant digits), organics_windows_<net>.csv (one row per
rolling fit), summary_organics_<net>.json (pooled metrics, paired bootstrap intervals over seeds, the acceptance key),
organics_<net>.png and first_storm_<net>.png.

THE SECONDARY TEST (--seasonal; reported, no bar; addendum 2 of the plan): O2 on task 10's seasonal truth (V1 plant
temperatures, the truth's own bulk E/R, task 10's weak wall), Net3, models B0, task 10b's M2 (temperature, TOC-blind),
M_TOC (TOC, temperature-blind), M2_TOC (M2's 15 bulk and wall pairs with every bulk rate also times TOC_m / TOC_ref;
M2TocSimGP24 over a KeyedWallBank) and the oracle of M2_TOC's bank.  Its 105 new blocks are built in memory, never
cached.  Outputs: outputs/chem/organics_seasonal_<net>.csv and summary_organics_seasonal_<net>.json.

AFTER THE REVIEW (--review; no bar, no scored number changes): outputs/chem/task11_review_<net>.json, a grid-only
control, the out-of-season fit read at the fitting months' TOC, and splits of the committed CSVs (see review()).
Simulation only.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import shutil
import sys
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .chemistry import TOC_REF_MGL, Chemistry, monthly_temperature, monthly_toc
from .seasonal import (BOOT_N, BOOT_SEED, LEVELS, MIN_FREE_GB_RUN, MONTHS, NOISE_SD, PER_MONTH, ROTATING, THRESHOLD,
                       Bank, SeasonalSimGP24, WallBank, WallSeasonalSimGP24, _build_conditions, _clean, _dose,
                       _map_sums, _reading_sums, _round_sig, _tkey, bootstrap_ratio, bootstrap_recall_ratio,
                       oracle_member, wall_bank, wall_pairs)
from .simgp import (DAY_HOURS, KB_GRID, SimGP24, disk_free_gb, disk_preflight, grid_cache_path, grid_edge_mass,
                    simulator_grid_24h)

FIXED_TEMP_C = 20.0
TRUTHS = {"O1": {"kinetics": "first", "phi": None}, "O2": {"kinetics": "clark", "phi": 0.85},
          "O3": {"kinetics": "clark", "phi": 0.5}}
TRUTH_LABELS = {
    "O1": "matched structure: first order, kb_m = kb_net u_m TOC_m / 2 (a debug check)",
    "O2": "Clark: dC/dt = -k2 u_m C (C - CL_m), CL_m = d_m - 0.85 TOC_m, k2 = kb_net / 1.7 (ASSUMPTIONS)",
    "O3": "Clark, low demand: phi = 0.5, k2 = kb_net / 1.0, CL = +0.2 mg/L at TOC 2.0 and dose 1.2 (ASSUMPTIONS)"}
VARIANTS = ("S", "C", "D")
VARIANT_LABELS = {"S": "the assumed TOC schedule (chemistry.monthly_toc)", "C": "control, TOC 2.0 mg/L every month",
                  "D": "S with the plant dose stepped from 1.2 to 1.5 mg/L in July (O2 on Net3 only; no bar)"}
DOSE_STEP = (7, 1.5)                 # month from which the plant dose is 1.5 mg/L in schedule D
CLASSES = {"falling": (4, 5, 6, 7), "steady": (8, 9, 10), "first_storm": (11,), "after_storm": (12,)}
ROLLING_TARGETS = tuple(range(4, 13))
EXTRAPOLATION = {"jan_mar_to_jul_sep": ((1, 2, 3), (7, 8, 9))}
STOP_CLASSES = (("rolling", "falling"), ("rolling", "steady"), ("rolling", "first_storm"), ("rolling", "after_storm"),
                ("jan_mar_to_jul_sep", "all_months"))
TOC_HYPOTHESES = ("TOC",)
TOC_PRIOR = (0.25, 0.75)             # (H0, H_TOC)
SENSITIVITY_PRIORS = {"M_TOC_h0half": (0.5, 0.5), "M_TOC_linear": (0.0, 1.0)}
FIRST_SEED = {"Net3": 32, "Net2": 16}
USED_SEEDS = {"Net3": range(32), "Net2": range(16)}     # scored by tasks 10 (Net3 0-15, Net2 0-7) and 10b (16-31, 8-15)
DOSE_EDGE_BAR, DOSE_EDGE_SHARE = 0.5, 0.90
SUM_COLS = ("n_true_low_all", "r_all_n", "r_all_sse", "r_all_sae", "r_all_low", "r_all_tp", "r_all_fp", "r_all_in50",
            "r_all_in80", "r_all_in90", "r_all_in95", "r_seen_n", "r_seen_in90", "r_new_n", "r_new_in90", "d_n", "d_sse",
            "d_sae", "d_se", "d_in90", "d_low", "d_tp", "d_fp", "n_flag_all",
            "d_n_a1", "d_sse_a1", "d_se_a1", "d_n_a2", "d_sse_a2", "d_se_a2", "d_n_a3", "d_sse_a3", "d_se_a3")
KEY_COLS = ("truth", "variant", "test", "seed", "cls", "model")
CSV_FLOAT = "%.6g"
SECONDARY_MODELS = ("B0", "M2", "M_TOC", "M2_TOC", "oracle")


def toc_class(month: int) -> str:
    for k, ms in CLASSES.items():
        if month in ms:
            return k
    return "winter"


# ----------------------------------------------------------------------------- schedules (ASSUMPTIONS)
def toc_schedule(truth: str, variant: str, months: int = MONTHS, seasonal: bool = False, dose: float = 1.2) -> list:
    """One dict per month for pilot.synthetic_log: the plant temperature (20 C, or task 10's V1 schedule when
    seasonal), the plant TOC (S and D: chemistry.monthly_toc; C: TOC_ref every month), the truth's kinetics and phi,
    and in D the month's plant dose."""
    if truth not in TRUTHS or variant not in VARIANTS:
        raise ValueError(f"truth must be one of {tuple(TRUTHS)} and variant one of {VARIANTS}")
    out = []
    for i in range(months):
        m = i % 12 + 1
        d = {"temp_C": float(monthly_temperature(m)) if seasonal else FIXED_TEMP_C, "soil_temp_C": None,
             "toc_mgL": float(TOC_REF_MGL if variant == "C" else monthly_toc(m))}
        if TRUTHS[truth]["kinetics"] != "first":
            d.update(kinetics=TRUTHS[truth]["kinetics"], phi=TRUTHS[truth]["phi"])
        if variant == "D":
            d["dose_mgL"] = float(DOSE_STEP[1] if m >= DOSE_STEP[0] else dose)
        out.append(d)
    return out


# ----------------------------------------------------------------------------- the TOC bank and model
@dataclass
class TocBank(Bank):
    """Blocks keyed by (plant TOC, hypothesis): blocks[(c, 'TOC')] is the committed grid re-run with every bulk rate
    times c / TOC_ref (chemistry.Chemistry(toc_mgL=c)); at TOC_ref it IS the committed grid.  h0 is hypothesis H0, the
    same at every TOC.  Bank's machinery (stack, blocks_at, merged) works unchanged with TOC in place of temperature."""

    @property
    def hypotheses(self) -> tuple:
        return ("H0",) + tuple(self.ers)

    @property
    def tocs(self) -> list:
        return self.temps


def toc_condition(toc_mgL: float) -> Chemistry:
    """The model condition of a TOC block: TOC scales the bulk rate, nothing else (no temperature)."""
    return Chemistry(toc_mgL=_tkey(toc_mgL))


def toc_bank(sc, tocs, grid: str = "full", cache_dir: str = "outputs/cache", cache: str = "readwrite",
             n_jobs: int | None = None) -> TocBank:
    """M_TOC's bank for the given TOC levels, as seasonal.covariate_bank: TOC_ref is the committed grid; every other
    level is read from its tagged cache file (simgp.grid_cache_path(sc, cache_dir, grid, Chemistry(toc_mgL=c))) or
    built, all in one process-pool call; cache 'readwrite', 'read' (never writes) or 'off' (builds in memory)."""
    if cache not in ("readwrite", "read", "off"):
        raise ValueError("cache must be 'readwrite', 'read' or 'off'")
    params, h0 = simulator_grid_24h(sc, cache_dir, grid)
    blocks, todo = {}, []
    for c in sorted({_tkey(t) for t in tocs}):
        cond = toc_condition(c)
        if cond.is_default:
            blocks[(c, "TOC")] = h0
            continue
        f = grid_cache_path(sc, cache_dir, grid, cond)
        if cache != "off" and os.path.exists(f):
            with open(f, "rb") as fh:
                p, Z = pickle.load(fh)
            if list(p) != list(params) or Z.shape != h0.shape:
                raise ValueError(f"{f} does not hold this network's grid")
            blocks[(c, "TOC")] = Z
        else:
            todo.append((c, cond, f))
    if todo:
        built = _build_conditions(sc, [cond for _, cond, _ in todo], grid, n_jobs, [cache_dir])
        for (c, cond, f), Z in zip(todo, built):
            blocks[(c, "TOC")] = Z
            if cache == "readwrite":
                disk_preflight([cache_dir])
                with open(f, "wb") as fh:
                    pickle.dump((params, Z), fh)
    return TocBank(list(params), h0, blocks, TOC_HYPOTHESES, "toc_linear", grid)


def toc_bank_cache_files(sc, tocs, grid: str = "full", cache_dir: str = "outputs/cache") -> list:
    return [grid_cache_path(sc, cache_dir, grid, toc_condition(c)) for c in sorted({_tkey(t) for t in tocs})
            if not toc_condition(c).is_default]


class TocSimGP24(SeasonalSimGP24):
    """Model M_TOC (see the module docstring).

        model = TocSimGP24(sc, X, bank).fit(samples, target_toc_mgL=3.0)
        model.predict_hours(); model.predict_daily_min()        # the month at TOC 3.0 mg/L
        model.set_target(1.5); model.predict_daily_min()        # the same fit, another month

    samples: junction, hour, y (mg/L), toc_mgL (the plant TOC of the sample's month) and optionally dose_ratio (the
    month's logged plant dose over the model's dose, an exact ln-offset).  prior_mass: one weight per hypothesis
    (H0, H_TOC); default TOC_PRIOR."""

    condition_column = "toc_mgL"

    def __init__(self, sc, X, bank: TocBank, prior_mass=None, **kw):
        if not isinstance(bank, TocBank):
            raise TypeError("TocSimGP24 needs a TocBank")
        if "prior" in kw:
            raise TypeError("M_TOC's prior is TOC_PRIOR or prior_mass, not a named prior")
        super().__init__(sc, X, bank, prior="uniform", **kw)
        p = np.asarray(TOC_PRIOR if prior_mass is None else prior_mass, dtype=float)
        if p.shape != (len(bank.hypotheses),) or (p < 0).any() or abs(float(p.sum()) - 1.0) > 1e-12:
            raise ValueError(f"prior_mass needs {len(bank.hypotheses)} non-negative weights summing to 1")
        self.prior = "toc" if prior_mass is None else "custom"
        self.prior_mass = p
        with np.errstate(divide="ignore"):
            self.log_prior = np.log(p[self.hyp_index] / bank.n_base)

    def fit(self, samples: pd.DataFrame, target_toc_mgL: float, target_dose_ratio: float = 1.0) -> "TocSimGP24":
        return super().fit(samples, target_temp_C=target_toc_mgL, target_dose_ratio=target_dose_ratio)

    def fit_prior(self, target_toc_mgL: float = TOC_REF_MGL, target_dose_ratio: float = 1.0) -> "TocSimGP24":
        out = super().fit_prior(target_toc_mgL, target_dose_ratio)
        self.samples_ = pd.DataFrame(columns=["junction", "hour", "y", "toc_mgL"])
        return out

    def set_target(self, toc_mgL: float, dose_ratio: float = 1.0) -> "TocSimGP24":
        """Predict the month at plant TOC toc_mgL (and the month's dose ratio): the same posterior and GP."""
        return super().set_target(toc_mgL, dose_ratio)

    @property
    def target_toc_mgL(self) -> float:
        return self.target_temp_C

    def kb_ref(self) -> float:
        """Posterior mean bulk rate over the H_TOC members, read at TOC_ref (the calibrated alpha TOC_ref)."""
        return self.rate20(0)

    def p_h0(self) -> float:
        return float(self.hypothesis_posterior()["H0"])


# ----------------------------------------------------------------------------- the secondary: M2 plus TOC
def cond_key(temp_C: float, toc_mgL: float) -> str:
    """The bank key of a (plant temperature, plant TOC) month in the secondary test."""
    return f"T{_tkey(temp_C):g}_TOC{_tkey(toc_mgL):g}"


@dataclass
class KeyedWallBank(WallBank):
    """A WallBank keyed by condition labels (cond_key) instead of temperatures: blocks[(key, (b, w))]."""

    def blocks_at(self, key) -> list:
        k = str(key)
        missing = [p for p in self.ers if (k, p) not in self.blocks]
        if missing:
            raise KeyError(f"the bank has no block at {k} for pairs {missing[:3]}")
        return [self.h0] + [self.blocks[(k, p)] for p in self.ers]

    def stack(self, key) -> np.ndarray:
        return np.concatenate(self.blocks_at(key), axis=0)


class M2TocSimGP24(WallSeasonalSimGP24):
    """M2_TOC, the secondary test's model: task 10b's M2 over a KeyedWallBank whose blocks scale every bulk rate by
    f(T; bulk E/R) TOC / TOC_ref and every wall rate by f(T; wall E/R).  Samples carry cond_key; set_target takes one."""

    condition_column = "cond_key"

    @staticmethod
    def _condition_key(value):
        return str(value)


def m2toc_bank(sc, months_tc, toc_bank_: TocBank, grid: str = "full", cache_dir: str = "outputs/cache",
               n_jobs: int | None = None) -> tuple[KeyedWallBank, WallBank]:
    """(M2_TOC's bank keyed by cond_key, M2's WallBank at the temperatures), for the (temperature, TOC) pairs of the
    months.  At TOC_ref a block is M2's own (read from task 10b's cache); at 20 C a block is the TOC bank's (f = 1, so
    every pair is the same physics); every other block is built in memory and never cached."""
    pairs = wall_pairs()
    temps = sorted({_tkey(t) for t, _ in months_tc})
    m2 = wall_bank(sc, temps, pairs, grid, cache_dir, cache="read", n_jobs=n_jobs)
    blocks, todo = {}, []
    for t, c in sorted({(_tkey(t), _tkey(c)) for t, c in months_tc}):
        key = cond_key(t, c)
        for p in pairs:
            if c == TOC_REF_MGL:
                blocks[(key, p)] = m2.blocks[(t, p)]
            elif t == FIXED_TEMP_C:
                blocks[(key, p)] = toc_bank_.blocks[(c, "TOC")]
            else:
                todo.append(((key, p), Chemistry(temp_C=t, er_K=p[0], wall_er_K=p[1], toc_mgL=c)))
    if todo:
        disk_preflight([cache_dir])
        conds = list({cond: None for _, cond in todo})
        built = dict(zip(conds, _build_conditions(sc, conds, grid, n_jobs, [cache_dir])))
        for k, cond in todo:
            blocks[k] = built[cond]
    return KeyedWallBank(m2.params, m2.h0, blocks, pairs, "wall_er_toc", grid), m2


# ----------------------------------------------------------------------------- one (network, seed, truth, schedule)
def _age_terciles(sc) -> pd.Series:
    """1-based nominal-age tercile of every junction (the operator's daily-mean water age; ties broken by order)."""
    a = sc.age_daily_mean_h.loc[sc.junctions]
    return pd.Series(pd.qcut(a.rank(method="first"), 3, labels=False) + 1, index=a.index)


def _age_sums(median: pd.Series, truth_min: pd.Series, junctions, terc: pd.Series) -> dict:
    out = {}
    e_all = median.loc[junctions].values - truth_min.loc[junctions].values
    out["d_se"] = float(e_all.sum())
    for k in (1, 2, 3):
        js = [j for j in junctions if terc[j] == k]
        e = median.loc[js].values - truth_min.loc[js].values
        out.update({f"d_n_a{k}": int(len(js)), f"d_sse_a{k}": float((e ** 2).sum()), f"d_se_a{k}": float(e.sum())})
    return out


def run_task(net: str, seed: int, truth: str, variant: str, banks: dict, cache_dir: str, months: int = MONTHS,
             targets=ROLLING_TARGETS, models=None, extrapolate: bool = True, seasonal: bool = False,
             want_fig: bool = False) -> dict:
    """Every model and test for one network, seed, truth and schedule.  banks: {'TOC': TocBank} (and, seasonal=True,
    'M2': WallBank, 'M2_TOC': KeyedWallBank).  Returns rows (one per test, month and model), the truth's months, and
    (want_fig) what the first-storm figure needs.  months, targets, models and extrapolate=False run a subset, for the
    saved check that recomputes committed rows; the rows it produces are the same numbers (every fit and draw is
    seeded)."""
    from .experiment import NET_TRUTH
    from .features import build_features
    from .pilot import synthetic_log
    from .simulate import nominal_scenario
    truth_kw = NET_TRUTH.get(net, {})
    dose0 = _dose(net)
    sc = nominal_scenario(net, 14, dose0)
    X = build_features(sc)
    terc = _age_terciles(sc)
    sched = toc_schedule(truth, variant, months, seasonal=seasonal, dose=dose0)
    log, _, plant, truths = synthetic_log(net, months=months, per_month=PER_MONTH, rotating=ROTATING, seed=seed,
                                          schedule=sched, return_truth=True, **truth_kw)
    log = log.assign(mi=[int(m[-2:]) for m in log.month])
    log["dose_ratio"] = log.dose_mgL / dose0 if "dose_mgL" in log else 1.0
    log["cond_key"] = [cond_key(t, c) for t, c in zip(log.temp_C, log.toc_mgL)]
    toc_of = {m: sched[m - 1]["toc_mgL"] for m in range(1, months + 1)}
    T_of = {m: sched[m - 1]["temp_C"] for m in range(1, months + 1)}
    ratio_of = {m: sched[m - 1].get("dose_mgL", dose0) / dose0 for m in range(1, months + 1)}
    key_of = {m: cond_key(T_of[m], toc_of[m]) for m in range(1, months + 1)}
    tmin = {t["month"]: t["truth_daily_min"] for t in truths}
    tbh = {t["month"]: t["truth_by_hour"] for t in truths}
    junctions = list(sc.junctions)
    jidx = {j: i for i, j in enumerate(junctions)}
    if seasonal:
        names = [n for n in ("B0", "M2", "M_TOC", "M2_TOC") if models is None or n in models]
    else:
        names = ["B0", "M_TOC"] + (list(SENSITIVITY_PRIORS) if variant == "S" else []) + (
            ["M_TOC_nodose"] if variant == "D" else [])
        names = [n for n in names if models is None or n in models]
    baselines_on = variant == "S" and not seasonal and (models is None or "persistence" in models)
    oracle_on = models is None or "oracle" in models
    orc_bank, orc_key = (banks["M2_TOC"], key_of) if seasonal else (banks["TOC"], toc_of)
    rows = []

    def win(ms):
        return log[log.mi.isin(ms)]

    def fit(name, ms, m):
        """name's model fitted on months ms, set to predict month m."""
        tr = win(ms)
        kw = dict(seed=seed, cache_dir=cache_dir, threshold=THRESHOLD)
        if name == "B0":
            return SimGP24(sc, X, **kw).fit(tr[["junction", "hour", "y"]])
        if name.startswith("M_TOC"):
            prior = SENSITIVITY_PRIORS.get(name, TOC_PRIOR)
            use_dose = name != "M_TOC_nodose"
            cols = ["junction", "hour", "y", "toc_mgL"] + (["dose_ratio"] if use_dose else [])
            return TocSimGP24(sc, X, banks["TOC"], prior_mass=prior, **kw).fit(
                tr[cols], target_toc_mgL=toc_of[m], target_dose_ratio=ratio_of[m] if use_dose else 1.0)
        if name == "M2":
            return WallSeasonalSimGP24(sc, X, banks["M2"], **kw).fit(tr[["junction", "hour", "y", "temp_C"]],
                                                                     target_temp_C=T_of[m])
        if name == "M2_TOC":
            return M2TocSimGP24(sc, X, banks["M2_TOC"], **kw).fit(tr[["junction", "hour", "y", "cond_key"]],
                                                                  target_temp_C=key_of[m])
        raise ValueError(name)

    def retarget(name, model, m):
        if name.startswith("M_TOC"):
            model.set_target(toc_of[m], ratio_of[m] if name != "M_TOC_nodose" else 1.0)
        elif name == "M2":
            model.set_target(T_of[m])
        elif name == "M2_TOC":
            model.set_target(key_of[m])
        return model

    def base(test_name, m, model_name, window, n_train):
        return {"network": net, "truth": truth, "variant": variant, "test": test_name, "seed": seed, "month": m,
                "cls": toc_class(m) if test_name == "rolling" else "all_months", "model": model_name, "window": window,
                "n_train": n_train, "toc_mgL": toc_of[m], "temp_C": T_of[m], "dose_ratio": ratio_of[m],
                "n_true_low_all": int((tmin[m] < THRESHOLD).sum())}

    def score(test_name, m, name, window, train, pred, model=None, unsampled=None):
        test = log[log.mi == m]
        seen = test.junction.isin(set(train.junction)).values
        z_mu, z_sd, pmin = pred
        ii, hh = np.array([jidx[j] for j in test.junction]), test.hour.values.astype(int)
        r = _reading_sums(z_mu[hh, ii], z_sd[hh, ii], test, seen, THRESHOLD)
        flag = pmin["p_below"] > 0.5
        d = _map_sums(pmin["median"], pmin["lo90"], pmin["hi90"], flag, tmin[m], unsampled, THRESHOLD)
        row = {**base(test_name, m, name, window, len(train)), **r, **d, "n_flag_all": int(flag.sum()),
               **_age_sums(pmin["median"], tmin[m], unsampled, terc)}
        if model is not None:
            P = np.asarray(model.params, float)
            row.update(map_kb=float(model.map_params_[0]), map_kw=float(model.map_params_[1]),
                       map_dose=float(model.map_dose_), post_kb=float(model.w_ @ P[:, 0]),
                       post_kw=float(model.w_ @ P[:, 1]))
            e = grid_edge_mass(model.W_, model.params, model.doses)
            row.update(dose_edge=e["dose"], kb_edge=e["kb"], kw_edge=e["kw"])
            if isinstance(model, SeasonalSimGP24):
                row.update(P_H0=float(model.hypothesis_posterior()["H0"]), map_hypothesis=model.map_hypothesis_)
                if isinstance(model, WallSeasonalSimGP24):
                    row.update(E_wall_post_mean_K=model.wall_er_posterior_mean(), E_post_mean_K=model.e_posterior_mean())
        rows.append(row)

    def score_point(test_name, m, name, window, train, read_pred, map_pred=None, unsampled=None):
        test = log[log.mi == m]
        seen = test.junction.isin(set(train.junction)).values
        r = _reading_sums(np.log(np.clip(read_pred, 1e-9, None)), None, test, seen, THRESHOLD)
        if map_pred is not None:
            flag = map_pred < THRESHOLD
            d = {**_map_sums(map_pred, None, None, flag, tmin[m], unsampled, THRESHOLD),
                 **_age_sums(map_pred, tmin[m], unsampled, terc)}
            nf = int(flag.sum())
        else:
            d, nf = {}, float("nan")
        rows.append({**base(test_name, m, name, window, len(train)), **r, **d, "n_flag_all": nf})

    def baselines(test_name, m, window, train):
        test = log[log.mi == m]
        last = train.sort_values("date").groupby("tap_id").y.last()
        pers = test.tap_id.map(last).fillna(train.y.mean()).values.astype(float)
        score_point(test_name, m, "persistence", window, train, pers)
        score_point(test_name, m, "network_mean", window, train, np.full(len(test), float(train.y.mean())))

    def oracle(test_name, m, window, train, unsampled):
        Z = orc_bank.stack(orc_key[m])
        lr = math.log(ratio_of[m])
        Zs = Z if lr == 0.0 else Z + np.float32(lr)
        k, dmult, rm = oracle_member(Zs, tbh[m], junctions)
        zz = Zs[k].astype(float) + math.log(dmult)
        test = log[log.mi == m]
        ii, hh = np.array([jidx[j] for j in test.junction]), test.hour.values.astype(int)
        dmin = pd.Series(np.exp(zz.min(axis=0)), index=junctions)
        score_point(test_name, m, "oracle", window, train, np.exp(zz[hh, ii]), dmin, unsampled)
        rows[-1].update(oracle_hypothesis=orc_bank.hypotheses[k // orc_bank.n_base], oracle_dose=dmult,
                        oracle_daytime_rmse=rm, map_kb=float(orc_bank.stacked_params[k][0]),
                        map_kw=float(orc_bank.stacked_params[k][1]))

    def predictions(model):
        z_mu, z_sd = model.predict_hours()
        return z_mu, z_sd, model.predict_daily_min()

    def wlabel(ms):
        return f"{min(ms)}-{max(ms)}"

    fig = None
    for m in targets:
        w3 = (m - 3, m - 2, m - 1)
        train = win(w3)
        uns = [j for j in junctions if j not in set(train.junction)]
        fitted = {}
        for name in names:
            mod = fit(name, w3, m)
            pred = predictions(mod)
            score("rolling", m, name, wlabel(w3), train, pred, mod, uns)
            fitted[name] = (mod, pred)
            if want_fig and m == 11 and name in ("B0", "M_TOC"):
                fig = fig or {"truth_min": tmin[11], "toc": {k: toc_of[k] for k in (8, 9, 10, 11)}}
                fig[name] = pred[2]
        if baselines_on:
            baselines("rolling", m, wlabel(w3), train)
        if oracle_on:
            oracle("rolling", m, wlabel(w3), train, uns)
        for test_name, (fit_ms, ex_targets) in (EXTRAPOLATION.items() if extrapolate else ()):
            if w3 != fit_ms:
                continue
            for t in ex_targets:
                for name in names:
                    mod, pred = fitted[name]
                    if name != "B0":
                        pred = predictions(retarget(name, mod, t))
                    score(test_name, t, name, wlabel(w3), train, pred, mod, uns)
                if baselines_on:
                    baselines(test_name, t, wlabel(w3), train)
                if oracle_on:
                    oracle(test_name, t, wlabel(w3), train, uns)
    monthly = []
    for t in truths:
        c, m = t["chem"], t["month"]
        cl = c.get("clark") or {}
        dmean = float(np.mean(list(c["source_doses_mgL"].values())))
        monthly.append({"truth": truth, "variant": variant, "seed": seed, "month": m, "toc_mgL": toc_of[m],
                        "temp_C": T_of[m], "n_true_low_all": int((t["truth_daily_min"] < THRESHOLD).sum()),
                        "n_junctions": len(junctions), "bulk_month_factor": c["bulk_month_factor"],
                        "bulk_temp_factor": c["bulk_temp_factor"], "mean_source_dose_mgL": dmean,
                        "plant_dose_mgL": sched[m - 1].get("dose_mgL", dose0),
                        "k2_L_per_mg_day": cl.get("k2_L_per_mg_day"), "limiting_mgL": cl.get("limiting_mgL"),
                        "kb_per_day_effective": c["kb_per_day_effective"],
                        "mean_daily_min_mgL": float(t["truth_daily_min"].mean()),
                        "mean_daily_min_per_dose": float(t["truth_daily_min"].mean()) / dmean})
    return {"net": net, "seed": seed, "truth": truth, "variant": variant, "rows": rows, "monthly": monthly, "fig": fig}


# ----------------------------------------------------------------------------- worker pool
_W: dict = {}


def _init_worker(net, cache_dir, tocs, workdir, extra_path=None, temps=None):
    warnings.filterwarnings("ignore")
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tb = toc_bank(sc, tocs, cache_dir=cache_dir, cache="read")
    _W.update(net=net, cache_dir=cache_dir, banks={"TOC": tb})
    if extra_path:            # the secondary: M2's cached bank, and M2_TOC's in-memory blocks from the main process
        m2 = wall_bank(sc, temps, cache_dir=cache_dir, cache="read")
        with open(extra_path, "rb") as fh:
            new = pickle.load(fh)
        blocks = {}
        for key, p, src in new["index"]:
            blocks[(key, p)] = (new["arrays"][src[1]] if src[0] == "new" else
                                m2.blocks[src[1]] if src[0] == "m2" else tb.blocks[src[1]])
        _W["banks"].update(M2=m2, M2_TOC=KeyedWallBank(m2.params, m2.h0, blocks, wall_pairs(), "wall_er_toc", m2.grid))


def _task(args):
    net, seed, truth, variant, want_fig, seasonal = args
    t0 = time.time()
    out = run_task(net, seed, truth, variant, _W["banks"], _W["cache_dir"], seasonal=seasonal, want_fig=want_fig)
    out["seconds"] = time.time() - t0
    return out


# ----------------------------------------------------------------------------- compact rows, pooling, bootstrap
def compact_rows(df: pd.DataFrame, models=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(agg, windows): agg sums every count and squared-error column over the target months of each (truth, schedule,
    test, seed, class, model); windows is one row per rolling fit (truth, schedule, seed, month) with the fitted rates,
    M_TOC's P(H0) and its posterior mass at the edges of the dose and kb axes."""
    if models is not None:
        df = df[df.model.isin(models)]
    keys = list(KEY_COLS)
    g = df.groupby(keys, sort=False)
    agg = g[[c for c in SUM_COLS if c in df.columns]].sum(min_count=1)
    agg["n_months"] = g.size()
    agg = agg.reset_index()
    r = df[(df.test == "rolling") & df.model.isin(["B0", "M_TOC", "M2", "M2_TOC"])]
    idx = ["truth", "variant", "seed", "month"]
    win = r[r.model == "B0"][idx + ["cls", "toc_mgL", "temp_C", "dose_ratio"]].set_index(idx)
    want = {"B0": ["map_kb", "post_kb", "dose_edge", "kb_edge"],
            "M_TOC": ["map_kb", "post_kb", "P_H0", "map_hypothesis", "dose_edge", "kb_edge", "kw_edge", "map_dose"],
            "M2": ["post_kb", "P_H0", "dose_edge"], "M2_TOC": ["post_kb", "P_H0", "dose_edge", "E_wall_post_mean_K"]}
    for model, cols in want.items():
        part = r[r.model == model]
        if len(part):
            win = win.join(part.set_index(idx)[cols].add_prefix(f"{model}_"))
    return agg, win.reset_index()


def pooled(df: pd.DataFrame) -> dict:
    """Pooled metrics of a set of aggregated rows (sums over seeds and months)."""
    s = df.sum(numeric_only=True)
    out = {"n_rows": int(len(df)), "n_months": int(s["n_months"])}
    n, low = s["r_all_n"], s["r_all_low"]
    r = {"n": int(n), "rmse": float(np.sqrt(s["r_all_sse"] / n)) if n else float("nan"),
         "mae": float(s["r_all_sae"] / n) if n else float("nan"), "n_low": int(low),
         "recall": float(s["r_all_tp"] / low) if low else float("nan"), "false_alarms": int(s["r_all_fp"])}
    tp, fp = s["r_all_tp"], s["r_all_fp"]
    r["precision"] = float(tp / (tp + fp)) if tp + fp > 0 else float("nan")
    for q in LEVELS:
        col = df[f"r_all_in{q}"]
        r[f"coverage{q}"] = float(col.sum() / n) if n and col.notna().all() else float("nan")
    for tag in ("seen", "new"):
        nn, col = s[f"r_{tag}_n"], df[f"r_{tag}_in90"]
        r[f"n_{tag}"] = int(nn)
        r[f"coverage90_{tag}"] = float(col.sum() / nn) if nn and col.notna().all() else float("nan")
    out["readings"] = r
    if "d_n" in df and df["d_n"].notna().all() and s["d_n"] > 0:
        dn, tp, fp = s["d_n"], s["d_tp"], s["d_fp"]
        dm = {"n": int(dn), "rmse": float(np.sqrt(s["d_sse"] / dn)), "mae": float(s["d_sae"] / dn),
              "bias": float(s["d_se"] / dn),
              "coverage90": float(df.d_in90.sum() / dn) if df.d_in90.notna().all() else float("nan"),
              "n_low": int(s["d_low"]), "recall": float(tp / s["d_low"]) if s["d_low"] else float("nan"),
              "false_alarms": int(fp), "precision": float(tp / (tp + fp)) if tp + fp > 0 else float("nan")}
        dm["by_age_tercile"] = {str(k): {"n": int(s[f"d_n_a{k}"]),
                                         "rmse": float(np.sqrt(s[f"d_sse_a{k}"] / s[f"d_n_a{k}"])) if s[f"d_n_a{k}"] else float("nan"),
                                         "bias": float(s[f"d_se_a{k}"] / s[f"d_n_a{k}"]) if s[f"d_n_a{k}"] else float("nan")}
                                for k in (1, 2, 3)}
        out["daily_min_map"] = dm
    return out


def _cls_rows(df: pd.DataFrame, cls: str) -> pd.DataFrame:
    return df if cls == "all_months" else df[df.cls == cls]


def pooled_tree(agg: pd.DataFrame) -> dict:
    """P[truth][schedule][test][class][model]; class all_months or one of CLASSES (rolling only)."""
    P = {}
    for (truth, v, test), g in agg.groupby(["truth", "variant", "test"], sort=True):
        for cls in ("all_months",) + (tuple(CLASSES) if test == "rolling" else ()):
            gc = _cls_rows(g, cls)
            if not len(gc):
                continue
            node = P.setdefault(truth, {}).setdefault(v, {}).setdefault(test, {}).setdefault(cls, {})
            for model, gm in gc.groupby("model", sort=True):
                node[model] = pooled(gm)
    return P


def _safe_boot(fn, rows, num, den, prefix):
    try:
        return fn(rows, num, den, prefix)
    except (IndexError, ValueError, ZeroDivisionError):
        return {"ratio": float("nan"), "lo90": float("nan"), "hi90": float("nan"), "n_seeds": None}


def bootstrap_tree(agg: pd.DataFrame, num: str = "M_TOC", others=()) -> dict:
    """Paired bootstrap over seeds (2000 resamples, 90% intervals) of pooled RMSE ratios num / B0 (and others / B0)."""
    B = {}
    for (truth, v, test), g in agg.groupby(["truth", "variant", "test"], sort=True):
        for cls in ("all_months",) + (tuple(CLASSES) if test == "rolling" else ()):
            gc = _cls_rows(g, cls)
            node = B.setdefault(truth, {}).setdefault(v, {}).setdefault(test, {}).setdefault(cls, {})
            for name in (num,) + tuple(others):
                if {name, "B0"} <= set(gc.model):
                    tag = "" if name == num else f"_{name}"
                    node[f"readings{tag}"] = _safe_boot(bootstrap_ratio, gc, name, "B0", "r_all")
                    node[f"daily_min_map{tag}"] = _safe_boot(bootstrap_ratio, gc, name, "B0", "d")
    return B


# ----------------------------------------------------------------------------- acceptance (pre-registered)
def _get(P, *keys):
    for k in keys:
        if not isinstance(P, dict) or k not in P:
            return None
        P = P[k]
    return P


def _stop_entry(P_cls: dict, rows: pd.DataFrame, num: str = "M_TOC", den: str = "B0") -> tuple[dict, list, list]:
    """One class of the stop rule: the table entry, its triggers (with counts and bootstrap intervals) and the false
    alarm rises above 10% (rises from zero included)."""
    ent, trig, fa = {}, [], []
    for kind, key, prefix in (("readings", "readings", "r_all"), ("daily_min_map", "daily_min_map", "d")):
        if key not in P_cls.get(num, {}) or key not in P_cls.get(den, {}):
            continue
        m, b = P_cls[num][key], P_cls[den][key]
        rr = m["rmse"] / b["rmse"]
        rc = m["recall"] / b["recall"] if b["recall"] and not math.isnan(b["recall"]) else float("nan")
        far = m["false_alarms"] / b["false_alarms"] if b["false_alarms"] else float("nan")
        s = rows.groupby("model")[[f"{prefix}_tp", f"{prefix}_fp"]].sum()
        # low_flagged: the low readings (or junctions) flagged, the true positives; n_flagged: every one flagged among
        # those scored (true positives plus false alarms).  Before the review both were one key, flagged_<model>,
        # that held the true positives; the numbers are unchanged.
        ent[kind] = {f"rmse_{num}": m["rmse"], "rmse_B0": b["rmse"], "rmse_ratio": rr, f"recall_{num}": m["recall"],
                     "recall_B0": b["recall"], "recall_ratio": rc, "n_low": b["n_low"],
                     f"low_flagged_{num}": int(s.loc[num].iloc[0]), "low_flagged_B0": int(s.loc[den].iloc[0]),
                     f"n_flagged_{num}": int(s.loc[num].sum()), "n_flagged_B0": int(s.loc[den].sum()),
                     f"false_alarms_{num}": m["false_alarms"], "false_alarms_B0": b["false_alarms"],
                     "false_alarm_ratio": far}
        if rr > 1.10:
            bt = bootstrap_ratio(rows, num, den, prefix)
            trig.append({"kind": kind, "metric": "rmse", num: m["rmse"], "B0": b["rmse"], "ratio": rr,
                         "n_scored": int(rows[rows.model == num][f"{prefix}_n"].sum()), "bootstrap": bt,
                         "interval_excludes_1": bool(bt["lo90"] > 1.0), "interval_past_bar": bool(bt["lo90"] > 1.10)})
        if not math.isnan(rc) and rc < 0.90:
            bt = _safe_boot(bootstrap_recall_ratio, rows, num, den, prefix)
            trig.append({"kind": kind, "metric": "recall", num: m["recall"], "B0": b["recall"], "ratio": rc,
                         "n_low": b["n_low"], f"low_flagged_{num}": ent[kind][f"low_flagged_{num}"],
                         "low_flagged_B0": ent[kind]["low_flagged_B0"], f"n_flagged_{num}": ent[kind][f"n_flagged_{num}"],
                         "n_flagged_B0": ent[kind]["n_flagged_B0"], "bootstrap": bt,
                         "interval_excludes_1": bool(bt["hi90"] < 1.0), "interval_past_bar": bool(bt["hi90"] < 0.90)})
        from_zero = not b["false_alarms"] and m["false_alarms"] > 0
        if from_zero or (not math.isnan(far) and far > 1.10):
            fa.append({"kind": kind, num: m["false_alarms"], "B0": b["false_alarms"], "ratio": far,
                       "from_zero": bool(from_zero)})
    return ent, trig, fa


def acceptance(net: str, P: dict, agg: pd.DataFrame, B: dict, win: pd.DataFrame) -> dict:
    """The pre-registered bars of task 11 (journal, task 11, written before the first scored run).  T1 and T2 are
    set for Net3 and reported for every network; T3 (the stop rule) and D apply on every network."""
    A = {"C1": {"bar": "outputs/chem/checks_report.json: every check passes, the task-11 checks included",
                "pass": None, "judged_from": "outputs/chem/checks_report.json after this run"}}
    # T1: matched structure (debug)
    g = _get(P, "O1", "S", "rolling", "all_months")
    t1 = None
    if g:
        m, b = g["M_TOC"]["daily_min_map"], g["B0"]["daily_min_map"]
        bt = _get(B, "O1", "S", "rolling", "all_months", "daily_min_map")
        t1 = {"bar": "truth O1 (first order, kb proportional to TOC), schedule S, rolling test, all target months: "
                     "pooled daily-minimum map RMSE_M_TOC < RMSE_B0; otherwise a bug, debugged before reporting",
              "applies": net == "Net3", "rmse_M_TOC": m["rmse"], "rmse_B0": b["rmse"], "ratio": m["rmse"] / b["rmse"],
              "bootstrap": bt, "pass": bool(m["rmse"] < b["rmse"])}
    A["T1"] = t1
    # T2: the first storm, Clark truth
    g = _get(P, "O2", "S", "rolling", "first_storm")
    t2 = None
    if g:
        m, b = g["M_TOC"]["daily_min_map"], g["B0"]["daily_min_map"]
        parts = {"rmse": m["rmse"] <= 0.85 * b["rmse"],
                 "recall": (m["recall"] >= b["recall"]) if not math.isnan(b["recall"]) else None}
        t2 = {"bar": "truth O2 (Clark), schedule S, first-storm month (November, fitted on August to October): "
                     "daily-minimum map RMSE_M_TOC <= 0.85 RMSE_B0 and recall_M_TOC >= recall_B0; a failure is recorded",
              "applies": net == "Net3", "rmse_M_TOC": m["rmse"], "rmse_B0": b["rmse"], "ratio": m["rmse"] / b["rmse"],
              "recall_M_TOC": m["recall"], "recall_B0": b["recall"], "n_low": b["n_low"],
              "false_alarms_M_TOC": m["false_alarms"], "false_alarms_B0": b["false_alarms"],
              "bootstrap": _get(B, "O2", "S", "rolling", "first_storm", "daily_min_map"), "parts": parts,
              "pass": bool(all(v is not False for v in parts.values()) and parts["rmse"])}
    A["T2"] = t2
    # N: in the control (TOC 2.0 every month) M_TOC is today's model by construction (both blocks are the committed
    # grid); only the daily minimum's Monte-Carlo draws differ (twice as many entries to draw from)
    nest = {}
    for truth in sorted(P):
        g = _get(P, truth, "C", "rolling", "all_months")
        if not g or "M_TOC" not in g:
            continue
        e = {}
        for kind in ("readings", "daily_min_map"):
            m, b = g["M_TOC"][kind], g["B0"][kind]
            e[kind] = {"rmse_ratio": m["rmse"] / b["rmse"], "recall_M_TOC": m["recall"], "recall_B0": b["recall"],
                       "recall_diff": m["recall"] - b["recall"]}
        nest[truth] = e
    A["N"] = {"bar": "control (TOC 2.0 every month), rolling test, every truth: |RMSE_M_TOC / RMSE_B0 - 1| <= 0.02 and "
                     "|recall_M_TOC - recall_B0| <= 0.03, for held-out readings and the daily-minimum map (M_TOC is "
                     "today's model there by construction; a failure is a bug)", "by_truth": nest,
              "pass": bool(nest) and all(abs(x["rmse_ratio"] - 1) <= 0.02 and abs(x["recall_diff"]) <= 0.03
                                         for e in nest.values() for x in e.values())}
    # T3: the stop rule
    trig, fa, table = [], [], {}
    for truth in ("O2", "O3"):
        for v in ("S",):
            for test, cls in STOP_CLASSES:
                Pc = _get(P, truth, v, test, cls)
                if not Pc or "M_TOC" not in Pc:
                    continue
                rows = agg[(agg.truth == truth) & (agg.variant == v) & (agg.test == test)]
                rows = _cls_rows(rows, cls)
                ent, tg, f = _stop_entry(Pc, rows)
                table[f"{truth}_{v}_{test}_{cls}"] = ent
                trig += [{"truth": truth, "variant": v, "test": test, "cls": cls, **x} for x in tg]
                fa += [{"truth": truth, "variant": v, "test": test, "cls": cls, **x} for x in f]
    A["T3"] = {"bar": "stop rule, truths O2 and O3, TOC schedule S, each class (rolling falling, steady, first storm, "
                      "after the storm; January to March predicting July to September): stop if RMSE_M_TOC > 1.10 "
                      "RMSE_B0 or recall_M_TOC < 0.90 recall_B0, for held-out readings or the daily-minimum map. Every "
                      "trigger is reported with its counts and a paired bootstrap 90% interval over seeds; robust = the "
                      "interval lies wholly on M_TOC's worse side of 1.0; both kinds count. False alarms up by more than "
                      "10% are listed as regressions.",
               "by_class": table, "triggered": trig,
               "robust_triggers": [x for x in trig if x["interval_excludes_1"]],
               "false_alarm_rises_over_10pct": fa, "stop": bool(trig)}
    # D: the dose axis
    w = win[win.variant.isin(["S", "C"])]
    de = w["M_TOC_dose_edge"].astype(float)
    share = float((de < DOSE_EDGE_BAR).mean()) if len(de) else float("nan")
    A["D"] = {"bar": f"M_TOC's posterior mass on the dose axis's edges (multipliers 0.90 and 1.10) is below "
                     f"{DOSE_EDGE_BAR} in at least {DOSE_EDGE_SHARE:.0%} of its rolling fits (every truth, schedules S and C)",
              "n_windows": int(len(de)), "share_below": share, "median_edge_mass": float(de.median()),
              "by_truth": {t: float((g_["M_TOC_dose_edge"] < DOSE_EDGE_BAR).mean()) for t, g_ in w.groupby("truth")},
              "pass": bool(share >= DOSE_EDGE_SHARE)}
    # G: gains claimed only from the bootstrap
    gains = {}
    for truth in ("O1", "O2", "O3"):
        for test, cls in (("rolling", "all_months"),) + STOP_CLASSES:
            for kind in ("readings", "daily_min_map"):
                bt = _get(B, truth, "S", test, cls, kind)
                if bt:
                    gains[f"{truth}_S_{test}_{cls}_{kind}"] = {"ratio": bt["ratio"], "lo90": bt["lo90"], "hi90": bt["hi90"],
                                                               "gain_claimed": bool(bt["hi90"] < 1.0)}
    A["G"] = {"bar": "a gain is claimed only where the upper end of the paired bootstrap 90% interval of "
                     "RMSE_M_TOC / RMSE_B0 is below 1.0", "by_case": gains}
    first = (bool(t1 and t1["pass"] and t2 and t2["pass"]) if net == "Net3" else True)
    adopt = bool(first and not A["T3"]["stop"] and A["D"]["pass"] and A["N"]["pass"])
    A["adoption"] = {"rule": "M_TOC becomes an optional logged input (pilot --plant toc_mgL, app, default off) only if "
                             "T1 and T2 pass on Net3, and T3 does not stop and N and D pass on every network",
                     "passes_on_this_network": adopt}
    return A


def dose_step_summary(monthly: pd.DataFrame) -> dict:
    """Schedule D against S for truth O2, months 7 to 12, the same seeds and draws (only the plant dose differs, 1.5
    against 1.2 mg/L): the truth's mean daily minimum per unit of the month's mean source dose.  First-order decay makes
    the two equal; Clark's kinetics do not (the plan's probe: 0.259 at 1.2 and 0.286 at 1.5 on one Net3 seed)."""
    d = monthly[(monthly.variant == "D") & (monthly.truth == "O2") & (monthly.month >= DOSE_STEP[0])]
    s = monthly[(monthly.variant == "S") & (monthly.truth == "O2") & (monthly.month >= DOSE_STEP[0])]
    if not len(d) or not len(s):
        return {}
    j = d.merge(s, on=["seed", "month"], suffixes=("_D", "_S"))
    a, b = float(j.mean_daily_min_per_dose_S.mean()), float(j.mean_daily_min_per_dose_D.mean())
    return {"months": [int(DOSE_STEP[0]), 12], "n_seed_months": int(len(j)),
            "mean_daily_min_per_dose_at_1.2_schedule_S": a, "mean_daily_min_per_dose_at_1.5_schedule_D": b,
            "ratio_1.5_over_1.2": b / a, "first_order_would_give": 1.0}


def _truth_low_by_month(monthly: pd.DataFrame) -> dict:
    out = {}
    for (t, v, m), g in monthly.groupby(["truth", "variant", "month"]):
        out.setdefault(t, {}).setdefault(v, {})[int(m)] = round(float(g.n_true_low_all.mean() / g.n_junctions.iloc[0]), 4)
    return out


def _p_h0_by_class(win: pd.DataFrame) -> dict:
    out = {}
    for (t, v), g in win[win.variant.isin(["S", "C"])].groupby(["truth", "variant"]):
        out.setdefault(t, {})[v] = {c: float(gc.M_TOC_P_H0.mean()) for c, gc in g.groupby("cls")}
    return out


def _settings(net: str, truths, variants) -> dict:
    from .seasonal import _truth_settings
    return {"months": MONTHS, "route_taps": PER_MONTH, "rotating_taps": ROTATING, "noise_sd_mgL": NOISE_SD,
            "sampling_hours": [min(DAY_HOURS), max(DAY_HOURS)], "threshold_mgL": THRESHOLD,
            "fixed_temperature_C": FIXED_TEMP_C, "toc_ref_mgL": TOC_REF_MGL,
            "month_seed": "1000 + 100 seed + month index (pilot.synthetic_log with a schedule)",
            "truths": {t: TRUTH_LABELS[t] for t in truths}, "variants": {v: VARIANT_LABELS[v] for v in variants},
            "toc_schedule_ASSUMPTION": [monthly_toc(m) for m in range(1, 13)], "dose_step": list(DOSE_STEP),
            "M_TOC": "hypotheses H0 (today's grid) and H_TOC (bulk rates times TOC_m / 2.0); prior H0 1/4, H_TOC 3/4",
            "sensitivity_priors": {k: list(v) for k, v in SENSITIVITY_PRIORS.items()},
            "classes": {k: list(v) for k, v in CLASSES.items()}, "rolling_targets": list(ROLLING_TARGETS),
            "extrapolation": {k: [list(a), list(b)] for k, (a, b) in EXTRAPOLATION.items()},
            "stop_classes": [list(x) for x in STOP_CLASSES],
            "map_junctions": "every junction not sampled in the 3-month window, the same set for every model",
            "age_terciles": "by the operator's nominal daily-mean water age over all junctions (1 youngest, 3 oldest)",
            "oracle": "the best single member and dose of M_TOC's bank at the target month's TOC (times the logged dose)",
            "bootstrap": {"resamples": BOOT_N, "rng_seed": BOOT_SEED, "unit": "seed"},
            "csv": {"float_format": CSV_FLOAT, "aggregated_over": "the target months of each class"},
            "truth": _truth_settings(net)}


def summarise(net: str, agg: pd.DataFrame, win: pd.DataFrame, monthly_summary: dict, seeds, settings: dict) -> dict:
    P = pooled_tree(agg)
    B = bootstrap_tree(agg, "M_TOC", ("oracle",) + tuple(SENSITIVITY_PRIORS) + ("M_TOC_nodose",))
    acc = acceptance(net, P, agg, B, win)
    bands = {}
    for t in sorted(P):
        g = _get(P, t, "S", "rolling", "all_months", "M_TOC", "readings")
        if g:
            bands[t] = {"seen_taps": {"coverage90": g["coverage90_seen"], "band": [0.85, 0.95]},
                        "new_taps": {"coverage90": g["coverage90_new"], "band": [0.80, 0.95]}}
            for b in bands[t].values():
                b["in_band"] = bool(b["band"][0] <= b["coverage90"] <= b["band"][1])
    return _round_sig({
        "generated_by": f"python -m residualmap.organics {net} {len(seeds)}", "network": net,
        "seeds": [int(s) for s in seeds],
        "about": "Simulation only (task 11). Fresh seeds. Hidden truths at a fixed 20 C whose bulk decay follows the "
                 "plant TOC: O1 first order proportional to TOC (matched structure), O2 Clark second-order chlorine-"
                 "organics kinetics with a limiting concentration, O3 the same with a low demand that runs out. Today's "
                 "TOC-blind model B0 against M_TOC (bulk rates times TOC / 2.0 as a hypothesis next to today's grid). "
                 "Pooled over seeds and target months; the CSVs hold the sums every number here is computed from. The "
                 "TOC schedule, phi and k2 are ASSUMPTIONS.",
        "settings": settings, "acceptance": acc, "pilot_coverage_bands_M_TOC": bands,
        "p_h0_by_class": _p_h0_by_class(win), "dose_step": monthly_summary.get("dose_step", {}),
        "pooled": P, "bootstrap": B, "truth_share_below_threshold_daily_min_by_month": monthly_summary.get("truth_low", {}),
        "truth_k2_and_limiting": monthly_summary.get("clark", {})})


def _clark_summary(monthly: pd.DataFrame) -> dict:
    out = {}
    m = monthly[monthly.k2_L_per_mg_day.notna()]
    for (t, v), g in m.groupby(["truth", "variant"]):
        out.setdefault(t, {})[v] = {"k2_L_per_mg_day_median": float(g.k2_L_per_mg_day.median()),
                                    "limiting_mgL_by_toc": {f"{c:g}": float(gc.limiting_mgL.mean()) for c, gc in g.groupby("toc_mgL")}}
    return out


def _monthly_summary(monthly: pd.DataFrame) -> dict:
    return {"truth_low": _truth_low_by_month(monthly), "dose_step": dose_step_summary(monthly),
            "clark": _clark_summary(monthly)}


# ----------------------------------------------------------------------------- figures
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_organics(net: str, summ: dict, out: str) -> None:
    """Map RMSE and recall of B0, M_TOC and the oracle in every class, for each truth (schedule S)."""
    plt = _plt()
    from .experiment import BIG
    P = summ["pooled"]
    cls_tests = [("rolling", c) for c in CLASSES] + [("jan_mar_to_jul_sep", "all_months")]
    labels = ["Apr to Jul\n(falling)", "Aug to Oct\n(steady)", "Nov\n(first storm)", "Dec\n(after)", "Jul to Sep\nfrom Jan to Mar"]
    truths = [t for t in ("O1", "O2", "O3") if t in P]
    models = (("B0", "tab:gray", "B0, today (TOC-blind)"), ("M_TOC", "tab:blue", "M_TOC, told the plant TOC"),
              ("oracle", "k", "oracle (best single member)"))
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(2, len(truths), figsize=(8 * len(truths), 11), squeeze=False)
        x = np.arange(len(cls_tests))
        for j, t in enumerate(truths):
            for i, (metric, ylab) in enumerate((("rmse", "daily-minimum map RMSE (mg/L)"), ("recall", "recall of junctions below 0.2"))):
                a = ax[i, j]
                for k, (mname, c, lab) in enumerate(models):
                    vals = [(_get(P, t, "S", test, cls, mname, "daily_min_map", metric) or float("nan")) for test, cls in cls_tests]
                    a.bar(x + (k - 1) * 0.27, [np.nan if v is None else v for v in vals], width=0.27, color=c, label=lab)
                a.set_xticks(x, labels, fontsize=10)
                a.set_ylabel(ylab)
                a.grid(alpha=0.3, axis="y")
                if i == 0:
                    a.set_title(f"truth {t}: {TRUTH_LABELS[t].split(':')[0]}", fontsize=13)
                if metric == "recall":
                    a.set_ylim(0, 1.05)
        ax[0, 0].legend(fontsize=10, loc="upper left")
        a2 = summ["acceptance"]["T2"]
        head = (f"first-storm map RMSE under the Clark truth, M_TOC {a2['rmse_M_TOC']:.3f} vs today's {a2['rmse_B0']:.3f} mg/L "
                f"(ratio {a2['ratio']:.3f})") if a2 else "daily-minimum map by class"
        names = {"falling": "Apr to Jul", "steady": "Aug to Oct", "first_storm": "Nov", "after_storm": "Dec",
                 "all_months": "Jul to Sep from Jan to Mar"}
        trig = summ["acceptance"]["T3"]["triggered"]
        stop = ("stop rule triggered (" + "; ".join(
            f"{x['truth']} {names[x['cls']]} {x['kind'].replace('daily_min_map', 'map')} {x['metric']} {x['ratio']:.3f}"
            + (", robust" if x["interval_excludes_1"] else ", within count noise") for x in trig) + ")") if trig else \
            "stop rule not triggered"
        import textwrap
        fig.suptitle(f"{net}, simulated 12-month logs at 20 C, {len(summ['seeds'])} seeds, assumed TOC schedule: {head}\n"
                     + "\n".join(textwrap.wrap(f"{stop}; the panels show the daily-minimum map only", 150)), fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_first_storm(net: str, fig_in: dict, summ: dict, seed: int, out: str) -> None:
    plt = _plt()
    import wntr
    from .experiment import BIG
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tm = fig_in["truth_min"]
    J = len(sc.junctions)
    size = max(14, int(4000 / J))
    vmax = float(max(tm.max(), fig_in["B0"]["median"].max(), fig_in["M_TOC"]["median"].max()))
    true_low = tm < THRESHOLD
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(1, 3, figsize=(24, 7))
        panels = [("truth", tm, true_low, f"TRUE November daily minimum (simulated, Clark truth O2)\n"
                                          f"{int(true_low.sum())} of {J} junctions below {THRESHOLD} mg/L"),
                  ("B0", fig_in["B0"]["median"], fig_in["B0"]["p_below"] > 0.5, None),
                  ("M_TOC", fig_in["M_TOC"]["median"], fig_in["M_TOC"]["p_below"] > 0.5, None)]
        for a, (name, val, ring, title) in zip(ax, panels):
            if title is None:
                tp, fp = int((ring & true_low).sum()), int((ring & ~true_low).sum())
                lab = "today's model, TOC-blind" if name == "B0" else "told November's TOC is 3.0 mg/L"
                title = f"{name}: {lab}\nfitted on August to October (TOC 1.5); flags {tp} of {int(true_low.sum())}, {fp} false alarms"
            wntr.graphics.plot_network(sc.wn, node_attribute=val.to_dict(), node_size=size, node_cmap="RdYlBu",
                                       node_range=(0, vmax), ax=a, link_width=0.6, add_colorbar=True, title=title)
            xy = np.array([sc.wn.get_node(j).coordinates for j in val.index[ring.values]]) if ring.any() else np.zeros((0, 2))
            if len(xy):
                a.scatter(xy[:, 0], xy[:, 1], s=size * 2.6, facecolors="none", edgecolors="k", linewidths=1.1, zorder=5)
        a2 = summ["acceptance"]["T2"]
        bt = a2.get("bootstrap") or {}
        bar = ("pre-registered bar 0.85 " + ("met" if a2["parts"]["rmse"] else "missed")) if a2["applies"] else \
            "the 0.85 bar is set for Net3"
        fig.suptitle(f"{net}, scenario {seed}: the first storm (plant TOC 1.5 to 3.0 mg/L, an ASSUMED schedule), map from "
                     f"August-to-October grab samples. Rings: below {THRESHOLD} mg/L (truth) or flagged.\nOver all "
                     f"{len(summ['seeds'])} seeds: November map RMSE M_TOC {a2['rmse_M_TOC']:.3f} vs B0 "
                     f"{a2['rmse_B0']:.3f} mg/L (ratio {a2['ratio']:.3f}, 90% interval {bt.get('lo90', float('nan')):.3f} "
                     f"to {bt.get('hi90', float('nan')):.3f}; {bar}), recall {a2['recall_M_TOC']:.2f} vs "
                     f"{a2['recall_B0']:.2f}, false alarms {a2['false_alarms_M_TOC']} vs {a2['false_alarms_B0']}",
                     fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------- the experiment
def _check_seeds(net, seeds):
    used = sorted(set(seeds) & set(USED_SEEDS.get(net, ())))
    if used:
        raise ValueError(f"seeds {used} were scored in tasks 10 or 10b on {net}; task 11 uses fresh seeds only")


def _pool(tasks, net, cache_dir, tocs, workers, extra_path=None, temps=None):
    workdir = tempfile.mkdtemp(prefix="rm_organics_")
    results = []
    try:
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            os.environ[k] = "1"
        n_workers = workers or max(1, min(len(tasks), (os.cpu_count() or 2) - 1))
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker,
                                 initargs=(net, cache_dir, tocs, workdir, extra_path, temps)) as ex:
            for res in ex.map(_task, tasks):
                results.append(res)
                print(f"  {net} {res['truth']} {res['variant']} seed {res['seed']}: {len(res['rows'])} rows in "
                      f"{res['seconds']:.0f} s", flush=True)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    return results


def _frames(results):
    df = pd.DataFrame([r for res in results for r in res["rows"]])
    order = {"rolling": 0, "jan_mar_to_jul_sep": 1}
    df = df.assign(_o=df.test.map(order)).sort_values(["truth", "variant", "_o", "seed", "month", "model"],
                                                     kind="mergesort").drop(columns="_o").reset_index(drop=True)
    monthly = pd.DataFrame([r for res in results for r in res["monthly"]]).sort_values(
        ["truth", "variant", "seed", "month"]).reset_index(drop=True)
    return df, monthly


ALL_TOCS = sorted({monthly_toc(m) for m in range(1, 13)} | {TOC_REF_MGL})


def run(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache", workers: int | None = None,
        truths=tuple(TRUTHS), variants=VARIANTS, tag: str = "") -> dict:
    """Task 11's primary experiment (see the module docstring)."""
    os.makedirs(outdir, exist_ok=True)
    cache_dir = os.path.abspath(cache_dir)
    t0 = time.time()
    _check_seeds(net, seeds)
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free, below {MIN_FREE_GB_RUN} GB")
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tb = time.time()
    n_new = sum(not os.path.exists(f) for f in toc_bank_cache_files(sc, ALL_TOCS, cache_dir=cache_dir))
    toc_bank(sc, ALL_TOCS, cache_dir=cache_dir)
    t_bank = time.time() - tb
    print(f"{net}: TOC bank ready in {t_bank:.0f} s ({len(ALL_TOCS)} TOC levels, {n_new} new cached blocks); "
          f"{disk_free_gb(cache_dir):.2f} GB free", flush=True)
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free after the bank, below {MIN_FREE_GB_RUN} GB")
    tasks = [(net, s, t, v, bool(t == "O2" and v == "S" and s == seeds[0]), False)
             for t in truths for v in variants for s in seeds if v != "D" or (t == "O2" and net == "Net3")]
    results = _pool(tasks, net, cache_dir, ALL_TOCS, workers)
    df, monthly = _frames(results)
    agg, win = compact_rows(df)
    p_agg = os.path.join(outdir, f"organics{tag}_{net}.csv")
    p_win = os.path.join(outdir, f"organics{tag}_windows_{net}.csv")
    agg.to_csv(p_agg, index=False, float_format=CSV_FLOAT)
    win.to_csv(p_win, index=False, float_format=CSV_FLOAT)
    agg = pd.read_csv(p_agg, float_precision="round_trip")          # every number of the summary comes from the files
    win = pd.read_csv(p_win, float_precision="round_trip")
    ms = _monthly_summary(monthly)
    summ = summarise(net, agg, win, ms, seeds, _settings(net, truths, variants))
    summ["monthly_truth_note"] = ("truth_share_below_threshold_daily_min_by_month, dose_step and truth_k2_and_limiting "
                                  "are computed from the truths of the run, not from the CSVs")
    with open(os.path.join(outdir, f"summary_organics{tag}_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    plot_organics(net, summ, os.path.join(outdir, f"organics{tag}_{net}.png"))
    fig_in = next((r["fig"] for r in results if r["fig"] is not None), None)
    if fig_in is not None and "B0" in fig_in and "M_TOC" in fig_in and summ["acceptance"]["T2"]:
        plot_first_storm(net, fig_in, summ, seeds[0], os.path.join(outdir, f"first_storm{tag}_{net}.png"))
    acc = summ["acceptance"]
    print(f"{net}: {time.time() - t0:.0f} s in all (bank {t_bank:.0f} s); T3 stop: {acc['T3']['stop']}", flush=True)
    return summ


# ----------------------------------------------------------------------------- the secondary (seasonal truth)
def secondary_months_tc(months: int = MONTHS) -> list:
    return [(float(monthly_temperature(i % 12 + 1)), float(monthly_toc(i % 12 + 1))) for i in range(months)]


def run_seasonal(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache",
                 workers: int | None = None) -> dict:
    """The secondary test (reported, no bar): O2 on task 10's V1 seasonal truth (weak wall, W1), models B0, M2, M_TOC,
    M2_TOC and the oracle; rolling test and January to March predicting July to September."""
    os.makedirs(outdir, exist_ok=True)
    cache_dir = os.path.abspath(cache_dir)
    t0 = time.time()
    _check_seeds(net, seeds)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tb_ = toc_bank(sc, ALL_TOCS, cache_dir=cache_dir, cache="read")
    mtc = secondary_months_tc()
    temps = sorted({t for t, _ in mtc})
    kb, m2 = m2toc_bank(sc, mtc, tb_, cache_dir=cache_dir)
    t_bank = time.time() - t0
    # the in-memory blocks go to the workers through one temporary file; M2's and the TOC bank's are read there
    arrays, index, seen = [], [], {}
    for (key, p), Z in kb.blocks.items():
        t, c = key[1:].split("_TOC")
        t, c = float(t), float(c)
        if c == TOC_REF_MGL:
            index.append((key, p, ("m2", (_tkey(t), p))))
        elif t == FIXED_TEMP_C:
            index.append((key, p, ("toc", (_tkey(c), "TOC"))))
        else:
            if id(Z) not in seen:
                seen[id(Z)] = len(arrays)
                arrays.append(Z)
            index.append((key, p, ("new", seen[id(Z)])))
    del kb, m2
    path = os.path.join(cache_dir, f"_tmp_organics_m2toc_{net}_{os.getpid()}.pkl")
    try:
        with open(path, "wb") as fh:
            pickle.dump({"index": index, "arrays": arrays}, fh, protocol=pickle.HIGHEST_PROTOCOL)
        n_arr = len(arrays)
        del arrays
        print(f"{net}: secondary banks ready in {t_bank:.0f} s ({n_arr} in-memory blocks, {os.path.getsize(path) / 1e9:.2f} GB "
              f"temporary file); {disk_free_gb(cache_dir):.2f} GB free", flush=True)
        if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
            raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free, below {MIN_FREE_GB_RUN} GB")
        tasks = [(net, s, "O2", "S", False, True) for s in seeds]
        results = _pool(tasks, net, cache_dir, ALL_TOCS, workers, extra_path=path, temps=temps)
    finally:
        if os.path.exists(path):
            os.remove(path)
    df, monthly = _frames(results)
    agg, win = compact_rows(df, models=SECONDARY_MODELS)
    p_agg = os.path.join(outdir, f"organics_seasonal_{net}.csv")
    agg.to_csv(p_agg, index=False, float_format=CSV_FLOAT)
    agg = pd.read_csv(p_agg, float_precision="round_trip")
    P = pooled_tree(agg)
    B = bootstrap_tree(agg, "M2_TOC", ("M2", "M_TOC"))
    summ = _round_sig({
        "generated_by": f"python -m residualmap.organics {net} {len(seeds)} --seasonal", "network": net,
        "seeds": [int(s) for s in seeds],
        "about": "Simulation only (task 11, the secondary test of the plan's addendum 2; reported, no bar). Truth O2 "
                 "(Clark chlorine-organics kinetics) on task 10's seasonal truth: the assumed plant temperature schedule "
                 "(V1), the truth's own bulk Arrhenius E/R, task 10's weak wall response (W1), and the assumed TOC "
                 "schedule. The ladder B0 (today) -> M2 (task 10b's temperature model, TOC-blind; not adopted) -> M2_TOC "
                 "(M2 with every bulk rate also times TOC / 2.0), with M_TOC (TOC only, temperature-blind) beside them.",
        "settings": {**_settings(net, ("O2",), ("S",)), "temperature_schedule_ASSUMPTION": [t for t, _ in mtc],
                     "models": list(SECONDARY_MODELS), "M2_TOC_new_blocks_in_memory": n_arr,
                     "wall_world": "W1, task 10's weak wall theta_w^(T - 20), theta_w ~ U(1.00, 1.07) (an ASSUMPTION)"},
        "pooled": P, "bootstrap": B,
        "truth_share_below_threshold_daily_min_by_month": _truth_low_by_month(monthly)})
    with open(os.path.join(outdir, f"summary_organics_seasonal_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    print(f"{net}: secondary done in {time.time() - t0:.0f} s", flush=True)
    return summ


def resummarise(net: str, outdir: str = "outputs/chem") -> dict:
    """Rewrite summary_organics_<net>.json from the committed CSVs (no rerun); the monthly truth summaries and the
    seeds are read from the summary being replaced, and the figures are not redrawn."""
    path = os.path.join(outdir, f"summary_organics_{net}.json")
    with open(path) as fh:
        old = json.load(fh)
    agg = pd.read_csv(os.path.join(outdir, f"organics_{net}.csv"), float_precision="round_trip")
    win = pd.read_csv(os.path.join(outdir, f"organics_windows_{net}.csv"), float_precision="round_trip")
    ms = {"truth_low": old["truth_share_below_threshold_daily_min_by_month"], "dose_step": old["dose_step"],
          "clark": old["truth_k2_and_limiting"]}
    truths, variants = tuple(old["settings"]["truths"]), tuple(old["settings"]["variants"])
    summ = summarise(net, agg, win, ms, old["seeds"], _settings(net, truths, variants))
    summ["monthly_truth_note"] = old.get("monthly_truth_note")
    with open(path, "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    return summ


def replot_first_storm(net: str, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache") -> dict:
    """Redraw first_storm_<net>.png from the committed summary and the first seed's November, refitted (B0 and M_TOC on
    August to October, truth O2, schedule S; every fit and draw is seeded, so the panels are the scored run's).  Returns
    each model's flags on that seed, for the record."""
    with open(os.path.join(outdir, f"summary_organics_{net}.json")) as fh:
        summ = json.load(fh)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tb = toc_bank(sc, ALL_TOCS, cache_dir=os.path.abspath(cache_dir), cache="read")
    seed = int(summ["seeds"][0])
    workdir = tempfile.mkdtemp(prefix="rm_organics_")
    here = os.getcwd()
    try:
        os.chdir(workdir)
        res = run_task(net, seed, "O2", "S", {"TOC": tb}, os.path.abspath(os.path.join(here, cache_dir)),
                       targets=(11,), models=("B0", "M_TOC"), extrapolate=False, want_fig=True)
    finally:
        os.chdir(here)
        shutil.rmtree(workdir, ignore_errors=True)
    plot_first_storm(net, res["fig"], summ, seed, os.path.join(outdir, f"first_storm_{net}.png"))
    tm = res["fig"]["truth_min"] < THRESHOLD
    return {m: {"tp": int(((res["fig"][m]["p_below"] > 0.5) & tm).sum()),
                "fp": int(((res["fig"][m]["p_below"] > 0.5) & ~tm).sum())} for m in ("B0", "M_TOC")}


# ----------------------------------------------------------------------------- after the review (no bar)
# Supporting numbers for the docs, written after the stop and the reviews; no scored number, interval or bar changes.
# (1) A grid-only control: M_TOC fitted and read with every month labelled TOC 1.5 mg/L, so it carries no TOC
#     information, only H_TOC's 1.5 mg/L block, whose bulk rates (0.075 to 0.525 per day) reach below today's 0.10 floor.
# (2) The January-to-March fit read in July to September at the fitting months' TOC (2.5, told nothing changed) as well
#     as at the logged 1.5, with the signed error of the held-out readings, which the scored outputs do not keep.
# (3) From the committed CSVs only: the map's RMSE split into bias and spread, the stop trigger's seeds, the recall
#     ratios of the out-of-season test under every truth and prior, and where today's fitted bulk rate sits on its grid.
REVIEW_GRID_TOC = 1.5
REVIEW_SUM = ("r_all_n", "r_all_sse", "r_all_se", "r_all_low", "r_all_tp", "r_all_fp", "d_n", "d_sse", "d_se", "d_low",
              "d_tp", "d_fp")
REVIEW_CHECKED = ("r_all_n", "r_all_sse", "r_all_low", "r_all_tp", "r_all_fp", "d_n", "d_sse", "d_se", "d_low", "d_tp",
                  "d_fp")


def review_task(net: str, seed: int, truth: str, banks: dict, cache_dir: str) -> list:
    """Rows (one per test, target month and model) of the review's recomputation for one seed and truth, schedule S.
    B0, M_TOC and M_TOC_linear are fitted exactly as run_task fits them (their sums are compared with the committed CSV);
    M_TOC_grid15 is the grid-only control; the *_at_2.5 rows read the January-to-March fit at TOC 2.5."""
    from .experiment import NET_TRUTH
    from .features import build_features
    from .pilot import synthetic_log
    from .simulate import nominal_scenario
    dose0 = _dose(net)
    sc = nominal_scenario(net, 14, dose0)
    X = build_features(sc)
    terc = _age_terciles(sc)
    sched = toc_schedule(truth, "S", MONTHS, dose=dose0)
    log, _, _, truths = synthetic_log(net, months=MONTHS, per_month=PER_MONTH, rotating=ROTATING, seed=seed,
                                      schedule=sched, return_truth=True, **NET_TRUTH.get(net, {}))
    log = log.assign(mi=[int(m[-2:]) for m in log.month])
    log["dose_ratio"] = log.dose_mgL / dose0 if "dose_mgL" in log else 1.0
    toc_of = {m: sched[m - 1]["toc_mgL"] for m in range(1, MONTHS + 1)}
    ratio_of = {m: sched[m - 1].get("dose_mgL", dose0) / dose0 for m in range(1, MONTHS + 1)}
    tmin = {t["month"]: t["truth_daily_min"] for t in truths}
    junctions = list(sc.junctions)
    jidx = {j: i for i, j in enumerate(junctions)}
    kw = dict(seed=seed, cache_dir=cache_dir, threshold=THRESHOLD)
    cols = ["junction", "hour", "y", "toc_mgL", "dose_ratio"]
    rows = []

    def score(test_name, m, name, train, model, uns):
        test = log[log.mi == m]
        seen = test.junction.isin(set(train.junction)).values
        z_mu, z_sd = model.predict_hours()
        pmin = model.predict_daily_min()
        ii, hh = np.array([jidx[j] for j in test.junction]), test.hour.values.astype(int)
        r = _reading_sums(z_mu[hh, ii], z_sd[hh, ii], test, seen, THRESHOLD)
        d = _map_sums(pmin["median"], pmin["lo90"], pmin["hi90"], pmin["p_below"] > 0.5, tmin[m], uns, THRESHOLD)
        a = _age_sums(pmin["median"], tmin[m], uns, terc)
        rows.append({"truth": truth, "seed": seed, "test": test_name, "month": m,
                     "cls": toc_class(m) if test_name == "rolling" else "all_months", "model": name,
                     "r_all_se": float((np.exp(z_mu[hh, ii]) - test.y.values.astype(float)).sum()),
                     **{k: r[k] for k in ("r_all_n", "r_all_sse", "r_all_low", "r_all_tp", "r_all_fp")},
                     **{k: d[k] for k in ("d_n", "d_sse", "d_low", "d_tp", "d_fp")}, "d_se": a["d_se"]})

    for m in ROLLING_TARGETS:
        w3 = (m - 3, m - 2, m - 1)
        train = log[log.mi.isin(w3)]
        uns = [j for j in junctions if j not in set(train.junction)]
        b0 = SimGP24(sc, X, **kw).fit(train[["junction", "hour", "y"]])
        mt = TocSimGP24(sc, X, banks["TOC"], prior_mass=TOC_PRIOR, **kw).fit(
            train[cols], target_toc_mgL=toc_of[m], target_dose_ratio=ratio_of[m])
        g15 = TocSimGP24(sc, X, banks["TOC"], prior_mass=TOC_PRIOR, **kw).fit(
            train[cols].assign(toc_mgL=REVIEW_GRID_TOC), target_toc_mgL=REVIEW_GRID_TOC, target_dose_ratio=ratio_of[m])
        for name, mod in (("B0", b0), ("M_TOC", mt), ("M_TOC_grid15", g15)):
            score("rolling", m, name, train, mod, uns)
        for test_name, (fit_ms, ex_targets) in EXTRAPOLATION.items():
            if w3 != fit_ms:
                continue
            lin = TocSimGP24(sc, X, banks["TOC"], prior_mass=SENSITIVITY_PRIORS["M_TOC_linear"], **kw).fit(
                train[cols], target_toc_mgL=toc_of[m], target_dose_ratio=ratio_of[m])
            fit_toc = float(np.mean([toc_of[k] for k in fit_ms]))
            for t in ex_targets:
                score(test_name, t, "B0", train, b0, uns)
                for name, mod in (("M_TOC", mt), ("M_TOC_linear", lin)):
                    score(test_name, t, name, train, mod.set_target(toc_of[t], ratio_of[t]), uns)
                    score(test_name, t, f"{name}_at_{fit_toc:g}", train, mod.set_target(fit_toc, ratio_of[t]), uns)
    return rows


def _review_worker(args):
    net, seed, truth = args
    return review_task(net, seed, truth, _W["banks"], _W["cache_dir"])


def _review_pooled(df: pd.DataFrame) -> dict:
    s = df[list(REVIEW_SUM)].sum()
    out = {}
    for p, kind in (("r_all", "readings"), ("d", "daily_min_map")):
        n, low, tp, fp = s[f"{p}_n"], s[f"{p}_low"], s[f"{p}_tp"], s[f"{p}_fp"]
        rmse, bias = math.sqrt(s[f"{p}_sse"] / n), s[f"{p}_se"] / n
        out[kind] = {"n": int(n), "rmse": rmse, "bias": bias, "spread": math.sqrt(max(rmse ** 2 - bias ** 2, 0.0)),
                     "n_low": int(low), "low_flagged": int(tp), "recall": float(tp / low) if low else float("nan"),
                     "false_alarms": int(fp)}
    return out


def _bias_spread(df: pd.DataFrame) -> dict:
    s = df[["d_n", "d_sse", "d_se"]].sum()
    rmse, bias = math.sqrt(s.d_sse / s.d_n), s.d_se / s.d_n
    return {"rmse": rmse, "bias": bias, "spread": math.sqrt(max(rmse ** 2 - bias ** 2, 0.0))}


def _csv_review(net: str, outdir: str) -> dict:
    """The review's numbers that come from the committed CSVs alone."""
    agg = pd.read_csv(os.path.join(outdir, f"organics_{net}.csv"), float_precision="round_trip")
    win = pd.read_csv(os.path.join(outdir, f"organics_windows_{net}.csv"), float_precision="round_trip")
    out = {}
    # the map's RMSE split into bias and spread (spread = sqrt(RMSE^2 - bias^2)), schedules S and C
    bs = {}
    for (truth, v, test), g in agg[agg.variant.isin(["S", "C"])].groupby(["truth", "variant", "test"], sort=True):
        for cls in ("all_months",) + (tuple(CLASSES) if test == "rolling" else ()):
            gc = _cls_rows(g, cls)
            e = {m: _bias_spread(gc[gc.model == m]) for m in ("B0", "M_TOC", "oracle") if (gc.model == m).any()}
            if "B0" in e and "M_TOC" in e:
                e["rmse_ratio_M_TOC_B0"] = e["M_TOC"]["rmse"] / e["B0"]["rmse"]
                e["spread_ratio_M_TOC_B0"] = e["M_TOC"]["spread"] / e["B0"]["spread"]
            bs.setdefault(truth, {}).setdefault(v, {}).setdefault(test, {})[cls] = e
    out["map_bias_and_spread"] = bs
    # the out-of-season reading recall: every truth and prior against B0, and B0's true positives minus each model's
    ex = agg[(agg.variant == "S") & (agg.test == "jan_mar_to_jul_sep")]
    rec = {}
    for truth, g in ex.groupby("truth", sort=True):
        tb = g[g.model == "B0"].set_index("seed").sort_index()
        e = {"n_low": int(tb.r_all_low.sum()), "B0": {"low_flagged": int(tb.r_all_tp.sum()),
                                                      "false_alarms": int(tb.r_all_fp.sum())}}
        for m in ("M_TOC",) + tuple(SENSITIVITY_PRIORS):
            gm = g[g.model == m].set_index("seed").sort_index()
            if not len(gm):
                continue
            d = (tb.r_all_tp - gm.r_all_tp).astype(int)
            e[m] = {"low_flagged": int(gm.r_all_tp.sum()), "false_alarms": int(gm.r_all_fp.sum()),
                    "recall_ratio_bootstrap": _safe_boot(bootstrap_recall_ratio, g, m, "B0", "r_all"),
                    "B0_minus_model_low_flagged_by_seed": {int(k): int(x) for k, x in d.items()},
                    "seeds_where_B0_flags_more": int((d > 0).sum()), "seeds_where_model_flags_more": int((d < 0).sum())}
        rec[truth] = e
    out["jan_mar_to_jul_sep_reading_recall"] = rec
    # where the fitted bulk rate sits: today's MAP kb at the grid's floor, and the kb-axis edge mass
    floor = float(min(KB_GRID))
    fl = {}
    for (truth, v), g in win[win.variant.isin(["S", "C"])].groupby(["truth", "variant"], sort=True):
        for cls, gc in list(g.groupby("cls", sort=True)) + [("all_months", g)]:
            fl.setdefault(truth, {}).setdefault(v, {})[cls] = {
                "n_windows": int(len(gc)), "B0_map_kb_at_floor_share": float((gc.B0_map_kb <= floor + 1e-12).mean()),
                "B0_kb_edge_mass_mean": float(gc.B0_kb_edge.mean()), "M_TOC_kb_edge_mass_mean": float(gc.M_TOC_kb_edge.mean())}
    out["kb_grid_floor"] = {"floor_per_day": floor, "H_TOC_floor_at_1.5_mgL": floor * 1.5 / TOC_REF_MGL,
                            "note": "kb_edge is the posterior mass at the kb axis's two ends (0.40 under a uniform "
                                    "posterior); M_TOC's map_kb is read at TOC_ref, so its effective floor is not shown",
                            "by_truth": fl}
    return out


def review(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache",
           workers: int | None = None) -> dict:
    """Write outputs/chem/task11_review_<net>.json (see the comment above REVIEW_GRID_TOC).  Nothing else is written."""
    cache_dir = os.path.abspath(cache_dir)
    t0 = time.time()
    _check_seeds(net, seeds)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    missing = [f for f in toc_bank_cache_files(sc, ALL_TOCS, cache_dir=cache_dir) if not os.path.exists(f)]
    if missing:
        raise RuntimeError(f"the TOC bank is not cached ({len(missing)} blocks missing); run the experiment first")
    tasks = [(net, s, t) for t in TRUTHS for s in seeds]
    workdir = tempfile.mkdtemp(prefix="rm_organics_")
    results = []
    try:
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            os.environ[k] = "1"
        n_workers = workers or max(1, min(len(tasks), (os.cpu_count() or 2) - 1))
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker,
                                 initargs=(net, cache_dir, ALL_TOCS, workdir)) as ex:
            for rows in ex.map(_review_worker, tasks):
                results += rows
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    df = pd.DataFrame(results)
    # every recomputed B0, M_TOC and M_TOC_linear sum against the committed CSV, at its 6 significant digits
    agg = pd.read_csv(os.path.join(outdir, f"organics_{net}.csv"), float_precision="round_trip")
    keys = ["truth", "test", "seed", "cls", "model"]
    mine = df[df.model.isin(["B0", "M_TOC", "M_TOC_linear"])].groupby(keys)[list(REVIEW_CHECKED)].sum()
    ref = agg[(agg.variant == "S")].set_index(keys)[list(REVIEW_CHECKED)]
    j = mine.join(ref, rsuffix="_csv", how="inner")
    same, worst = 0, 0.0
    for c in REVIEW_CHECKED:
        a = np.array([float(CSV_FLOAT % x) for x in j[c].values])
        b = j[f"{c}_csv"].values.astype(float)
        same += int((a == b).sum())
        worst = max(worst, float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-12))) if len(a) else 0.0)
    cross = {"n_values": int(len(j) * len(REVIEW_CHECKED)), "n_identical_at_6_significant_digits": same,
             "max_relative_difference": worst, "rows_compared": int(len(j)),
             "rows_expected": int(len(mine))}
    # pooled and bootstrapped, per truth and class
    roll, ext = df[df.test == "rolling"], df[df.test != "rolling"]
    ctl = {}
    for truth, g in roll.groupby("truth", sort=True):
        for cls in ("all_months",) + tuple(CLASSES):
            gc = _cls_rows(g, cls)
            sums = gc.groupby(["model", "seed"], as_index=False)[list(REVIEW_SUM)].sum()
            e = {m: _review_pooled(sums[sums.model == m]) for m in ("B0", "M_TOC", "M_TOC_grid15")}
            for m in ("M_TOC", "M_TOC_grid15"):
                e[f"{m}_vs_B0"] = {"readings_rmse_ratio": _safe_boot(bootstrap_ratio, sums, m, "B0", "r_all"),
                                   "map_rmse_ratio": _safe_boot(bootstrap_ratio, sums, m, "B0", "d")}
            ctl.setdefault(truth, {})[cls] = e
    exo = {}
    for truth, g in ext.groupby("truth", sort=True):
        sums = g.groupby(["model", "seed"], as_index=False)[list(REVIEW_SUM)].sum()
        e = {m: _review_pooled(gm) for m, gm in sums.groupby("model", sort=True)}
        for m in sorted(set(sums.model) - {"B0"}):
            e[f"{m}_vs_B0"] = {"readings_rmse_ratio": _safe_boot(bootstrap_ratio, sums, m, "B0", "r_all"),
                               "reading_recall_ratio": _safe_boot(bootstrap_recall_ratio, sums, m, "B0", "r_all"),
                               "map_rmse_ratio": _safe_boot(bootstrap_ratio, sums, m, "B0", "d")}
        exo[truth] = e
    summ = _round_sig({
        "generated_by": f"python -m residualmap.organics {net} {len(seeds)} --review", "network": net,
        "seeds": [int(s) for s in seeds],
        "about": "Simulation only. Task 11's supporting numbers written after the stop rule triggered and after the "
                 "reviews; not a bar, and no scored number, interval or bar changes. Schedule S, the scored seeds and "
                 "logs. grid_only_control: M_TOC_grid15 is M_TOC fitted and read with every month labelled TOC 1.5 "
                 "mg/L, so it carries no TOC information, only H_TOC's 1.5 mg/L block (bulk rates 0.075 to 0.525 per "
                 "day, below today's 0.10 floor). jan_mar_to_jul_sep_read_at_fitting_toc: the January-to-March fit read "
                 "in July to September at the logged TOC 1.5 (the scored M_TOC) and at the fitting months' 2.5 (told "
                 "nothing changed). Readings' bias is the mean signed error of the held-out readings (prediction minus "
                 "reading). The csv_only section is computed from organics_<net>.csv and organics_windows_<net>.csv "
                 "alone. Paired bootstrap over seeds as the summary's (2000 resamples, the same RNG seed).",
        "recomputed_against_committed_csv": cross,
        "grid_only_control": ctl, "jan_mar_to_jul_sep_read_at_fitting_toc": exo, "csv_only": _csv_review(net, outdir)})
    with open(os.path.join(outdir, f"task11_review_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    print(f"{net}: review written in {time.time() - t0:.0f} s; {same} of {cross['n_values']} recomputed values "
          f"identical to the committed CSV", flush=True)
    return summ


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("net")
    ap.add_argument("seeds", type=int, help="number of seeds, from --first-seed")
    ap.add_argument("--first-seed", type=int, default=None, help="default 32 on Net3, 16 on Net2 (fresh seeds)")
    ap.add_argument("--outdir", default="outputs/chem")
    ap.add_argument("--cache", default="outputs/cache")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--truths", default=",".join(TRUTHS))
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--tag", default="", help="suffix for the output names (smoke runs)")
    ap.add_argument("--seasonal", action="store_true", help="the secondary test on the seasonal truth (no bar)")
    ap.add_argument("--resummarise", action="store_true", help="rewrite the summary from the committed CSVs")
    ap.add_argument("--replot", action="store_true", help="redraw organics_<net>.png from the committed summary")
    ap.add_argument("--replot-first-storm", action="store_true",
                    help="redraw first_storm_<net>.png (the first seed's November refitted; the subtitle from the summary)")
    ap.add_argument("--review", action="store_true",
                    help="after the review: write task11_review_<net>.json (grid-only control, out-of-season reading "
                         "at the fitting months' TOC, the CSV-only splits); no scored output changes")
    a = ap.parse_args(argv)
    warnings.filterwarnings("ignore")
    if a.resummarise:
        resummarise(a.net, a.outdir)
        return 0
    if a.replot:
        with open(os.path.join(a.outdir, f"summary_organics_{a.net}.json")) as fh:
            plot_organics(a.net, json.load(fh), os.path.join(a.outdir, f"organics_{a.net}.png"))
        return 0
    if a.replot_first_storm:
        print(f"{a.net}: first seed's flags {replot_first_storm(a.net, a.outdir, a.cache)}")
        return 0
    first = a.first_seed if a.first_seed is not None else FIRST_SEED.get(a.net, 0)
    seeds = tuple(range(first, first + a.seeds))
    if a.review:
        review(a.net, seeds, a.outdir, a.cache, a.workers)
        return 0
    if a.seasonal:
        run_seasonal(a.net, seeds, a.outdir, a.cache, a.workers)
        return 0
    truths = tuple(t for t in a.truths.split(",") if t)
    variants = tuple(v for v in a.variants.split(",") if v)
    s = run(a.net, seeds, a.outdir, a.cache, a.workers, truths, variants, a.tag)
    return 2 if s["acceptance"]["T3"]["stop"] else 0


if __name__ == "__main__":
    sys.exit(main())
