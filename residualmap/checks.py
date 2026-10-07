"""
checks.py: the saved checks for the chemistry work (iteration 4, journal tasks 8 to 14).

    python -m residualmap.checks            # every check, 20 to 40 s on this machine
    python -m residualmap.checks --quick    # skips the fresh grid, the synthetic pilot and the app (a few seconds)

Plain asserts on purpose (pytest is not in .venv).  Exits non-zero if any check fails.  A full run writes
outputs/chem/checks_report.json, which holds only results that do not change from run to run (no dates, timings
or free disk), so rerunning it on an unchanged tree leaves git status unchanged.  Every run (full or --quick)
also writes outputs/chem/checks_run.json, git-ignored, with the date, the timings and the free disk.  Every
EPANET run happens inside a temporary working directory whose outputs/cache is a link to the repo's cache, so
nothing is written to the repo except those two files, and no grid cache is ever written outside outputs/cache.

The baseline-reproduction mode reruns the committed experiments in a scratch directory and compares every file
they write with the committed outputs, file by file and CSV row by CSV row, and runs the app's default demo:

    python -m residualmap.checks --reproduce /tmp/rm_rerun --label post_edit
    python -m residualmap.checks --compare /tmp/rm_rerun --label pre_edit --code "..." --code-root <old tree>

Results go to outputs/chem/baseline_reproduction.json (one entry per label).  It exits non-zero if a rerun
fails or a compared file was not written by the rerun.
"""
from __future__ import annotations

import argparse
import contextlib
import glob
import json
import math
import os
import pickle
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:      # grid workers are spawned processes and the checks run from a temporary cwd
    sys.path.insert(0, REPO)
OUT_DIR = os.path.join(REPO, "outputs", "chem")
CACHE = os.path.join(REPO, "outputs", "cache")
REPORT = os.path.join(OUT_DIR, "checks_report.json")       # committed: run-invariant results only
RUN_LOG = os.path.join(OUT_DIR, "checks_run.json")          # git-ignored: date, timings, free disk
REPRO = os.path.join(OUT_DIR, "baseline_reproduction.json")
PY = sys.executable

# anchors from committed files (CHANGELOG 2026-09-21 task 6 line; outputs/summary_Net3.json; outputs/pilot)
APP_DEMO = {"violations": 41, "junctions": 92, "found": 36, "unsampled_violations": 36, "false_alarms": 1,
            "flagged": 42, "sample_mean_mgL": 0.53}
CHAIN_TOL_MGL = 0.01            # the bar for every analytic check (EPANET's default quality tolerance; the unit
                                # recipe itself runs EPANET at a 1e-5 tolerance)
SI_VS_LEGACY_TOL_MGL = 1e-3
# sha256 of the Net3 'full' grid (params, dtype, shape, values) in outputs/cache/grid24_full_Net3.pkl, built on
# 2026-09-21 16:36 by the committed code.  outputs/cache is git-ignored, so the hash is what is committed: a fresh
# build is compared with it, and with the cached file when that file predates the run.
GRID24_FULL_NET3_SHA256 = "bc849e767dc6e82e41ed3a3bc09e4b673fba8cb8b6698b15eb26efba006c99fb"
CACHE_BEFORE: dict | None = None    # outputs/cache snapshot taken by run_checks before the first check


# ----------------------------------------------------------------------------- tiny harness
CHECKS: list[tuple[str, str, bool, callable]] = []


def check(group: str, quick: bool = True):
    def deco(fn):
        CHECKS.append((fn.__name__, group, quick, fn))
        return fn
    return deco


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return None if math.isnan(x) else x
    if isinstance(x, np.bool_):
        return bool(x)
    return x


@contextlib.contextmanager
def scratch_cwd():
    """A temporary working directory with outputs/cache linked to the repo's cache."""
    old = os.getcwd()
    with tempfile.TemporaryDirectory(prefix="rm_checks_") as tmp:
        os.makedirs(os.path.join(tmp, "outputs"))
        os.symlink(CACHE, os.path.join(tmp, "outputs", "cache"))
        os.chdir(tmp)
        try:
            yield tmp
        finally:
            os.chdir(old)


def _git_status() -> str:
    return subprocess.run(["git", "-C", REPO, "status", "--porcelain", "--untracked-files=all"],
                          capture_output=True, text=True).stdout


def _cache_snapshot() -> dict:
    """(size, mtime_ns) of every file in outputs/cache, which git ignores and so git status cannot watch."""
    out = {}
    for f in sorted(os.listdir(CACHE)) if os.path.isdir(CACHE) else []:
        st = os.stat(os.path.join(CACHE, f))
        out[f] = (st.st_size, st.st_mtime_ns)
    return out


TAGGED_CACHE = re.compile(r"_(free|ca)_.*_[0-9a-f]{10}\.pkl$")   # a chemistry condition's tagged grid name


def _grid_sha256(params, Z) -> str:
    import hashlib
    h = hashlib.sha256()
    h.update(json.dumps([[float(v) for v in p] for p in params]).encode())
    h.update(f"{Z.dtype.str}{Z.shape}".encode())
    h.update(np.ascontiguousarray(Z).tobytes())
    return h.hexdigest()


def _root_scratch_files() -> list[str]:
    return sorted(os.path.basename(f) for f in glob.glob(os.path.join(REPO, "msx*")) + glob.glob(os.path.join(REPO, "en*"))
                  if os.path.isfile(f) and re.fullmatch(r"(msx|en)[A-Za-z0-9]{6}", os.path.basename(f)))


# ----------------------------------------------------------------------------- units and helpers
@check("units")
def arrhenius_values():
    from .chemistry import ER_HYPOTHESES_K, arrhenius
    v10, v25 = arrhenius(10, 8000), arrhenius(25, 8000)
    assert round(v10, 3) == 0.381 and round(v25, 3) == 1.580, (v10, v25)
    assert all(arrhenius(20.0, e) == 1.0 for e in ER_HYPOTHESES_K)
    ratios = {e: arrhenius(20, e) / arrhenius(10, e) for e in ER_HYPOTHESES_K}
    assert [round(r, 2) for r in ratios.values()] == [1.83, 2.62, 4.24], ratios
    return {"f10_E8000": v10, "f25_E8000": v25, "ratio_20_over_10": ratios}


@check("units")
def water_property_ratios():
    from .chemistry import diffusivity_ratio, viscosity_ratio
    got = {"nu10": viscosity_ratio(10), "nu25": viscosity_ratio(25), "D10": diffusivity_ratio(10), "D25": diffusivity_ratio(25)}
    want = {"nu10": 1.299, "nu25": 0.891, "D10": 0.744, "D25": 1.142}
    assert all(round(got[k], 3) == want[k] for k in want), got
    assert viscosity_ratio(20.0) == 1.0 and diffusivity_ratio(20.0) == 1.0
    return got


@check("units")
def seasonal_schedules():
    from .chemistry import monthly_temperature, monthly_toc
    T = [monthly_temperature(m) for m in range(1, 13)]
    assert T == [10.5, 10, 10.5, 12.5, 15, 17.5, 19.5, 20, 19.5, 17.5, 15, 12.5], T
    assert len(set(T)) == 7
    toc = [monthly_toc(m) for m in range(1, 13)]
    assert toc == [2.5, 2.5, 2.5, 2.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 3.0, 2.5], toc
    return {"temp_C": T, "toc_mgL": toc, "label": "ASSUMPTIONS, not measured schedules"}


@check("units")
def chemistry_thresholds_and_species():
    from .chemistry import CHLORAMINE, Chemistry
    free, ca = Chemistry(), Chemistry(disinfectant=CHLORAMINE)
    assert free.threshold == 0.2 and free.species == "free chlorine"
    assert ca.threshold == 0.5 and ca.species == "total chlorine"
    assert "not a California rule" in ca.threshold_note and "detectable" in ca.threshold_note
    assert Chemistry(disinfectant=CHLORAMINE, threshold_mgL=0.3).threshold == 0.3
    for bad in ({"disinfectant": "chlorine dioxide"}, {"kinetics": "order3"}, {"temp_C": 50.0}, {"wall_mode": "x"}):
        try:
            Chemistry(**bad)
        except ValueError:
            continue
        raise AssertionError(f"Chemistry accepted {bad}")
    assert Chemistry().is_default and Chemistry(temp_C=20.0, toc_mgL=2.0).is_default
    assert not Chemistry(temp_C=12.5).is_default and not ca.is_default
    n = Chemistry(temp_C=20.0, toc_mgL=2.0, er_K=12000.0)
    assert n.kb_scale() == 1.0 and n.kw_scale() == 1.0
    return {"free": free.threshold_note, "chloramine": ca.threshold_note}


@check("units")
def cache_tags_and_names():
    from . import simgp
    from .chemistry import CHLORAMINE, Chemistry, cache_tag
    from .simulate import nominal_scenario
    base = Chemistry(temp_C=12.5)
    variants = {"temp": base.with_(temp_C=15.0), "er": base.with_(er_K=5000.0), "toc": base.with_(toc_mgL=3.0),
                "kinetics": base.with_(kinetics="first_si"), "wall": base.with_(wall_mode="mass_transfer_only"),
                "disinfectant": base.with_(disinfectant=CHLORAMINE)}
    t0 = cache_tag(base)
    assert t0 == cache_tag(Chemistry(temp_C=12.5)) and len(t0) == 10
    tags = {k: cache_tag(v) for k, v in variants.items()}
    assert len(set(tags.values()) | {t0}) == len(tags) + 1, tags
    assert cache_tag(base.with_(threshold_mgL=0.4)) == t0          # the threshold does not change a grid
    assert cache_tag(base, "decay") != t0
    saved = simgp.GRIDS["full"]
    try:
        simgp.GRIDS["full"] = (saved[0][:-1],) + saved[1:]
        changed = cache_tag(base)
    finally:
        simgp.GRIDS["full"] = saved
    assert changed != t0, "editing GRIDS must change the tag"
    from . import chemistry
    saved_v = chemistry.CHEM_CACHE_VERSION
    try:
        chemistry.CHEM_CACHE_VERSION = saved_v + 1
        bumped = cache_tag(base)
    finally:
        chemistry.CHEM_CACHE_VERSION = saved_v
    assert bumped != t0, "bumping CHEM_CACHE_VERSION must change the tag"
    # equal conditions share one tag whatever numeric type they were given in
    same = [Chemistry(temp_C=10), Chemistry(temp_C=10.0), Chemistry(temp_C=np.int64(10)), Chemistry(temp_C=np.float32(10.0))]
    assert len({cache_tag(c) for c in same}) == 1 and all(c == same[0] for c in same), [cache_tag(c) for c in same]
    assert cache_tag(Chemistry(temp_C=12.5, toc_mgL=3)) == cache_tag(Chemistry(temp_C=12.5, toc_mgL=3.0))
    assert cache_tag(Chemistry(temp_C=12.5, er_K=8000)) == t0
    sc = nominal_scenario("Net3")
    names = {k: os.path.basename(simgp.grid_cache_path(sc, CACHE, "full", c))
             for k, c in (("none", None), ("neutral", Chemistry(temp_C=20.0, toc_mgL=2.0)), ("T12.5", base))}
    assert names["none"] == names["neutral"] == "grid24_full_Net3.pkl", names
    assert names["T12.5"] == f"grid24_full_Net3_free_T12.5_E8000_{t0}.pkl", names
    return {"tag_T12.5_E8000": t0, "names": names}


@check("units")
def disk_preflight_floor():
    from .simgp import MIN_FREE_GB, disk_preflight
    free = disk_preflight([CACHE, tempfile.gettempdir()])
    try:
        disk_preflight([CACHE], min_gb=1e9)
    except RuntimeError:
        pass
    else:
        raise AssertionError("preflight did not refuse an impossible floor")
    assert free >= MIN_FREE_GB
    return {"floor_gb": MIN_FREE_GB, "_volatile": {"free_gb": free}}


# ----------------------------------------------------------------------------- EPANET against closed forms
def _chain(n_seg=12, length=1000.0, diam=0.3, flow=0.01, days=3, quality_step=60):
    """Reservoir, then n_seg pipes in series, then one demand: plug flow, travel time k L / v to junction k."""
    import wntr
    wn = wntr.network.WaterNetworkModel()
    wn.add_reservoir("R", base_head=100.0)
    prev = "R"
    for k in range(1, n_seg + 1):
        wn.add_junction(f"J{k}", base_demand=(flow if k == n_seg else 0.0), elevation=0.0)
        wn.add_pipe(f"P{k}", prev, f"J{k}", length=length, diameter=diam, roughness=130)
        prev = f"J{k}"
    wn.options.time.duration = int(days * 86400)
    wn.options.time.hydraulic_timestep = 3600
    wn.options.time.report_timestep = 3600
    wn.options.time.quality_timestep = quality_step
    wn.options.quality.parameter = "CHEMICAL"
    wn.options.reaction.wall_coeff = 0.0
    t_days = np.arange(1, n_seg + 1) * length / (flow / (np.pi * diam ** 2 / 4)) / 86400.0
    return wn, t_days


def _run_last(wn, scale):
    from .simulate import _run_quality
    q = _run_quality(wn)
    return q.iloc[-1][[n for n in wn.junction_name_list]].values * scale


CHAIN_CASES = {   # kind: (order, coeff in .inp units, limiting mg/L, analytic k)
    "order1": (1, -0.8, None, 0.8),
    "order2": (2, -0.8 / 1.5, None, 0.8 / 1.5),
    "order0": (0, -0.8 * 1.5, None, 0.8 * 1.5),
    "order1_limit": (1, -0.8, 0.6, 0.8),
    "clark": (2, -0.6, 0.4, 0.6),
    "clark_negative_CL": (2, -0.6, -0.5, 0.6),
}


@check("epanet_analytic")
def chain_kinetics_match_closed_forms():
    """Each EPANET kinetic form on a 12-pipe plug-flow chain, under the corrected unit recipe, against its closed
    form: max |error| at the 12 junctions within 0.01 mg/L."""
    from .chemistry import analytic_concentration, set_bulk_kinetics, use_mg_per_litre
    c0, out = 1.5, {}
    for kind, (order, coeff, cl, k) in CHAIN_CASES.items():
        wn, t = _chain()
        scale = use_mg_per_litre(wn, {"R": c0})
        set_bulk_kinetics(wn, order, coeff, cl)
        sim = _run_last(wn, scale)
        exact = analytic_concentration("order1_limit" if kind == "order1_limit" else
                                       ("clark" if kind.startswith("clark") else kind), c0, k, t, cl or 0.0)
        out[kind] = {"max_abs_err_mgL": float(np.abs(sim - exact).max()), "C_last": float(sim[-1]),
                     "exact_last": float(exact[-1]), "t_last_days": float(t[-1])}
    bad = {k: v for k, v in out.items() if not v["max_abs_err_mgL"] <= CHAIN_TOL_MGL}
    assert not bad, bad
    return out


@check("epanet_analytic")
def chain_legacy_units_break_order2():
    """Negative control: the repo's legacy source (1.5 written as 1500 mg/L) under order 2 is far from the exact
    answer, so the unit recipe is what makes the nonlinear checks pass."""
    from .chemistry import analytic_concentration
    wn, t = _chain()
    wn.add_source("src", "R", "CONCEN", 1.5)
    wn.options.reaction.bulk_order = wn.options.reaction.tank_order = 2
    wn.options.reaction.bulk_coeff = -0.8 / 1.5
    sim = _run_last(wn, 1.0)
    exact = analytic_concentration("order2", 1.5, 0.8 / 1.5, t)
    err = float(abs(sim[-1] - exact[-1]))
    assert err > 0.5, err
    return {"legacy_C_last": float(sim[-1]), "exact_last": float(exact[-1]), "abs_err": err}


@check("epanet_analytic")
def grid_member_leaves_no_files():
    from .simgp import FLOOR, HOURS, _grid_member
    from .simulate import nominal_scenario
    sc = nominal_scenario("Net3")
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "g0")
        z = _grid_member(("Net3", 0.4, 0.7, 1.0, 1.2, 1.0, 1.0, list(sc.junctions), prefix))
        left = os.listdir(tmp)
    assert left == [], left
    assert z.shape == (len(HOURS), len(sc.junctions)) and np.isfinite(z).all() and z.min() >= np.log(FLOOR) - 1e-6
    return {"files_left": left}


# ----------------------------------------------------------------------------- default path anchors
@check("default_path")
def anchor_net3_seed0():
    from .simulate import build_scenario
    summ = json.load(open(os.path.join(REPO, "outputs", "summary_Net3.json")))["scenario0"]
    sc = build_scenario("Net3", 0)
    got = {"below_threshold_at_sampling_hour": int((sc.truth_snapshot < 0.2).sum()),
           "below_threshold_daily_min": int((sc.truth_daily_min < 0.2).sum()),
           "pct_below_by_hour": {str(int(h)): round(float(v) * 100, 1) for h, v in (sc.truth_by_hour < 0.2).mean(axis=1).items()}}
    assert got["below_threshold_at_sampling_hour"] == summ["below_threshold_at_sampling_hour"] == 11, got
    assert got["below_threshold_daily_min"] == summ["below_threshold_daily_min"] == 41, got
    assert got["pct_below_by_hour"] == summ["pct_below_by_hour"]
    return {k: got[k] for k in ("below_threshold_at_sampling_hour", "below_threshold_daily_min")}


def _inp_text(prefix):
    """The .inp WNTR wrote, without its '; Created: <date and time>' header line, which changes every second
    and would make two identical files compare unequal across a second boundary."""
    with open(prefix + ".inp") as fh:
        return "".join(line for line in fh if not line.startswith("; Created:"))


@check("default_path")
def nominal_defaults_write_identical_inp():
    """simulate_nominal_chlorine's new keywords at their defaults (and at 20 C) write the same .inp byte for byte;
    at 10 C they change VISCOSITY and DIFFUSIVITY by the expected ratios and nothing else in [OPTIONS]."""
    from .chemistry import diffusivity_ratio, remove_epanet_files, viscosity_ratio
    from .simulate import simulate_nominal_chlorine
    texts, res = {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        for tag, kw in (("plain", {}), ("explicit", dict(kb_scale=1.0, kw_scale=1.0, temp_C=None)),
                        ("T20", dict(temp_C=20.0)), ("T10", dict(temp_C=10.0))):
            p = os.path.join(tmp, tag)
            res[tag] = simulate_nominal_chlorine("Net3", 0.4, 0.7, 1.0, 1.2, 1.0, 1.0, file_prefix=p, **kw)
            texts[tag] = _inp_text(p)
            remove_epanet_files(p)
    assert texts["plain"] == texts["explicit"] == texts["T20"]
    assert res["plain"].equals(res["explicit"]) and res["plain"].equals(res["T20"])

    def opt(txt, key):
        return float(re.search(rf"^{key}\s+(\S+)", txt, re.M).group(1))
    v, d = opt(texts["T10"], "VISCOSITY"), opt(texts["T10"], "DIFFUSIVITY")
    assert abs(v - viscosity_ratio(10)) < 1e-9 and abs(d - diffusivity_ratio(10)) < 1e-9, (v, d)
    diff = [(a, b) for a, b in zip(texts["plain"].splitlines(), texts["T10"].splitlines()) if a != b]
    assert len(diff) == 2 and all(a.split()[0] in ("VISCOSITY", "DIFFUSIVITY") for a, _ in diff), diff
    shift = float((res["T10"] - res["plain"]).values.mean())
    return {"T10_viscosity": v, "T10_diffusivity": d, "T10_mean_shift_mgL": shift}


@check("default_path")
def ky4_diffusivity_is_multiplied():
    from .chemistry import apply_water_temperature, diffusivity_ratio
    from .simulate import load
    wn = load("ky4")
    d0 = wn.options.quality.diffusivity
    apply_water_temperature(wn, 10.0)
    assert d0 == 0.05 and abs(wn.options.quality.diffusivity - 0.05 * diffusivity_ratio(10.0)) < 1e-15
    return {"file_diffusivity": d0, "at_10C": wn.options.quality.diffusivity}


# ----------------------------------------------------------------------------- chemistry truths
TRUTHS = {"Net3": {}, "Net2": dict(kb_per_day=0.10, kw_m_per_day=0.20)}   # experiment.NET_TRUTH


@check("chemistry_truth")
def neutral_chemistry_is_bit_identical():
    """build_scenario(chem=neutral) equals build_scenario() bit for bit: the chem branch consumes the committed
    draws in the same order and every new factor is exactly 1."""
    from .chemistry import Chemistry
    from .simulate import build_scenario
    out = {}
    cases = [("Net3", 0, {}), ("Net2", 0, {}), ("Net3", 1, dict(structural_noise="persistent")),
             ("Net3", 2, dict(month_seed=123))]
    for net, seed, extra in cases:
        kw = {**TRUTHS[net], **extra}
        ref = build_scenario(net, seed, **kw).truth_by_hour
        for tag, chem in (("T20_TOCref", Chemistry(temp_C=20.0, toc_mgL=2.0)), ("unset", Chemistry())):
            got = build_scenario(net, seed, chem=chem, **kw).truth_by_hour
            same = bool(np.array_equal(ref.values, got.values) and ref.index.equals(got.index) and ref.columns.equals(got.columns))
            out[f"{net}_seed{seed}_{'_'.join(extra) or 'default'}_{tag}"] = same
    assert all(out.values()), out
    return out


@check("chemistry_truth")
def corrected_units_match_legacy_first_order():
    from .chemistry import Chemistry
    from .simulate import build_scenario
    out = {}
    for net in ("Net3", "Net2"):
        a = build_scenario(net, 0, **TRUTHS[net]).truth_by_hour
        b = build_scenario(net, 0, chem=Chemistry(kinetics="first_si"), **TRUTHS[net]).truth_by_hour
        out[net] = float((a - b).abs().values.max())
    assert all(v <= SI_VS_LEGACY_TOL_MGL for v in out.values()), out
    return {"max_abs_diff_mgL": out}


@check("chemistry_truth")
def warm_water_truth_decays_faster_and_draws_are_isolated():
    """Same seed at 10 C and 25 C: every junction's daily-mean chlorine is lower at 25 C; the hidden E/R and
    theta_w are the same at both temperatures; the truth's water age (hydraulics only) is identical at 10 C,
    25 C and on the default path, so temperature consumed no committed draw."""
    from .chemistry import Chemistry
    from .simulate import build_scenario
    cold = build_scenario("Net3", 0, chem=Chemistry(temp_C=10.0), truth_age=True)
    warm = build_scenario("Net3", 0, chem=Chemistry(temp_C=25.0), truth_age=True)
    plain = build_scenario("Net3", 0, truth_age=True)
    ref = build_scenario("Net3", 0)
    assert np.array_equal(plain.truth_by_hour.values, ref.truth_by_hour.values), "truth_age changed the chlorine truth"
    dc, dw = cold.truth_by_hour.mean(), warm.truth_by_hour.mean()
    assert (dw <= dc + 1e-9).all() and float((dc - dw).mean()) > 0.0
    assert np.array_equal(cold.truth_age_by_hour_h.values, warm.truth_age_by_hour_h.values)
    assert np.array_equal(cold.truth_age_by_hour_h.values, plain.truth_age_by_hour_h.values)
    assert cold.chem["E_true_K"] == warm.chem["E_true_K"] and cold.chem["theta_w"] == warm.chem["theta_w"]
    assert cold.chem["rng_seed"] == 20_000
    other = build_scenario("Net3", 1, chem=Chemistry(temp_C=10.0))
    assert other.chem["E_true_K"] != cold.chem["E_true_K"]
    age_corr = float(np.corrcoef(plain.truth_age_by_hour_h.mean(), plain.age_daily_mean_h)[0, 1])
    assert age_corr > 0.8, age_corr
    return {"mean_mgL_10C": float(dc.mean()), "mean_mgL_25C": float(dw.mean()), "E_true_K_seed0": cold.chem["E_true_K"],
            "theta_w_seed0": cold.chem["theta_w"], "truth_vs_nominal_age_corr": age_corr}


@check("chemistry_truth")
def chloramine_refused_until_built():
    """Chloramine has no decay physics before task 12, so a chloramine truth, grid or grid cache lookup raises
    NotImplementedError instead of running free-chlorine physics under a total-chlorine label.  Its hidden draws
    already have their own generator, default_rng(40_000 + seed)."""
    from .chemistry import CHLORAMINE, Chemistry
    from .simgp import grid_cache_path, simulator_grid_24h
    from .simulate import build_scenario, hidden_chem_draws, nominal_scenario
    ca = Chemistry(disinfectant=CHLORAMINE, temp_C=15.0)
    refused = {}
    sc = nominal_scenario("Net3")
    before = os.path.exists(grid_cache_path(sc, CACHE, "decay", ca))
    for name, call in (("truth", lambda: build_scenario("Net3", 3, chem=ca)),
                       ("truth_neutral", lambda: build_scenario("Net3", 3, chem=Chemistry(disinfectant=CHLORAMINE))),
                       ("sim_kwargs", ca.sim_kwargs),
                       ("grid", lambda: simulator_grid_24h(sc, CACHE, "decay", cond=ca))):
        try:
            call()
        except NotImplementedError:
            refused[name] = True
        else:
            refused[name] = False
    assert all(refused.values()), refused
    assert os.path.exists(grid_cache_path(sc, CACHE, "decay", ca)) == before, "a refused grid call wrote a cache file"
    d_ca, d_fc = hidden_chem_draws(CHLORAMINE, 3), hidden_chem_draws("free_chlorine", 3)
    assert d_ca["rng_seed"] == 40_003 and d_fc["rng_seed"] == 20_003 and d_ca["E_true_K"] != d_fc["E_true_K"]
    fc = build_scenario("Net3", 3, chem=Chemistry(temp_C=15.0))
    assert fc.chem["rng_seed"] == 20_003 and fc.chem["E_true_K"] == d_fc["E_true_K"] and fc.chem["species"] == "free chlorine"
    assert ca.species == "total chlorine"
    return {"refused": refused, "chloramine_rng_seed": d_ca["rng_seed"], "free_chlorine_rng_seed": fc.chem["rng_seed"]}


# ----------------------------------------------------------------------------- the model
def _fitted_net3(threshold=0.2, n=8, seed=0):
    from .features import build_features
    from .simgp import DAY_HOURS, SimGP24
    from .simulate import build_scenario
    sc = build_scenario("Net3", seed)
    X = build_features(sc)
    rng = np.random.default_rng(seed)
    js = list(rng.choice(sc.junctions, n, replace=False)); hs = [int(h) for h in rng.choice(DAY_HOURS, n)]
    y = [float(np.clip(sc.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None)) for j, h in zip(js, hs)]
    S = pd.DataFrame({"junction": js, "hour": hs, "y": y})
    return sc, X, S, SimGP24(sc, X, seed=seed, cache_dir=CACHE, threshold=threshold).fit(S)


@check("model")
def threshold_and_monte_carlo_p_below():
    """SimGP24(threshold=0.2) keeps the committed arithmetic; P(daily min < t) is a Monte-Carlo count at any t
    through the model (p_below_mc, or p_below(frame, t, model=m)); the threshold changes only the comparison,
    never the draws; without the model, SimGP24.p_below keeps the committed normal approximation away from the
    frame's own threshold; the frame carries no draws (plain, serialisable attrs); map_params_ stays a 5-tuple."""
    from scipy.stats import norm
    from .simgp import SimGP24
    sc, X, S, m = _fitted_net3()
    pmin = m.predict_daily_min()
    assert np.array_equal(pmin["p_below"].values, (m.zmin_ < np.log(0.2)).mean(axis=0))
    assert SimGP24.p_below(pmin).equals(pmin["p_below"])
    p05 = m.p_below_mc(0.5)
    assert SimGP24.p_below(pmin, 0.5, model=m).equals(p05) and (p05 >= pmin["p_below"]).all()
    normal = SimGP24.p_below(pmin, 0.5)                    # no model: the committed normal approximation
    want = pd.Series(norm.cdf((np.log(0.5) - pmin["z_mu"]) / pmin["z_sd"]), index=pmin.index)
    assert normal.equals(want)
    sub = pmin.loc[sc.junctions[:10]]
    assert SimGP24.p_below(sub, 0.5, model=m).equals(p05.loc[sc.junctions[:10]])
    edited = pmin.copy(); edited["z_mu"] += 5.0            # an edited copy must not read the model's draws
    for bad, exc in ((edited, ValueError), (pmin.reset_index(drop=True), KeyError)):
        try:
            SimGP24.p_below(bad, 0.5, model=m)
        except exc:
            continue
        raise AssertionError(f"p_below accepted a frame its draws did not produce ({exc.__name__} expected)")
    assert SimGP24.p_below(pmin.reset_index(drop=True), 0.5).notna().all()   # no model: works on any frame
    assert set(pmin.attrs) == {"p_below_threshold"} and json.dumps(pmin.attrs) and m.zmin_.flags.writeable
    m5 = SimGP24(sc, X, seed=0, cache_dir=CACHE, threshold=0.5).fit(S)
    p5 = m5.predict_daily_min()
    assert np.array_equal(p5["p_below"].values, p05.values) and SimGP24.p_below(p5, 0.5).equals(p5["p_below"])
    for c in ("median", "lo90", "hi90", "z_mu", "z_sd"):
        assert np.array_equal(p5[c].values, pmin[c].values), c
    assert SimGP24.p_below(p5, 0.2, model=m5).equals(m.p_below_mc(0.2))
    assert SimGP24.p_below(p5, 0.2).equals(pd.Series(norm.cdf((np.log(0.2) - p5["z_mu"]) / p5["z_sd"]), index=p5.index))
    assert len(m.map_params_) == 5
    m.fit(S)
    assert m.zmin_ is None, "a refit must clear the previous daily-minimum draws"
    return {"n_flagged_0.2": int((pmin["p_below"] > 0.5).sum()), "n_flagged_0.5_mc": int((p05 > 0.5).sum()),
            "n_flagged_0.5_normal_approx": int((normal > 0.5).sum()),
            "max_abs_diff_mc_vs_normal_0.5": float((p05 - normal).abs().max())}


@check("model")
def edge_mass_helper():
    from .simgp import DOSE_GRID, GRIDS, grid_edge_mass
    import itertools
    params = list(itertools.product(*GRIDS["full"]))
    W = np.full((len(params), len(DOSE_GRID)), 1.0 / (len(params) * len(DOSE_GRID)))
    e = grid_edge_mass(W, params, DOSE_GRID)
    want = {"kb": 0.4, "kw": 0.4, "gamma": 2 / 3, "demand": 2 / 3, "rough": 2 / 3, "dose": 0.4}
    assert all(abs(e[k] - v) < 1e-12 for k, v in want.items()), e
    W2 = np.zeros_like(W); W2[0, 2] = 1.0              # all mass on kb, kw, gamma, demand, rough minimum, dose 1.0
    e2 = grid_edge_mass(W2, params, DOSE_GRID)
    assert e2["kb_low"] == 1.0 and e2["kb_high"] == 0.0 and e2["dose"] == 0.0
    decay = list(itertools.product(*GRIDS["decay"]))
    e3 = grid_edge_mass(np.full(75, 1 / 75), decay)   # member marginal only
    assert math.isnan(e3["demand"]) and "dose" not in e3
    e4 = grid_edge_mass(np.full((75, 1), 1 / 75), decay, doses=[1.0])   # experiment.py's one-dose decay baseline
    assert math.isnan(e4["dose"]) and abs(e4["kb"] - 0.4) < 1e-12
    for bad in (dict(W=W, params=params), dict(W=W, params=params, doses=DOSE_GRID[:-1]),
                dict(W=np.full((75, 1), 1 / 75), params=decay), dict(W=W[:-1], params=params, doses=DOSE_GRID)):
        try:
            grid_edge_mass(**bad)
        except ValueError:
            continue
        raise AssertionError("grid_edge_mass accepted a posterior that does not match its params or doses")
    return {"uniform": e}


# ----------------------------------------------------------------------------- full-run anchors
@check("anchors", quick=False)
def fresh_grid_equals_cached_grid():
    """A fresh 675-run Net3 grid (built in memory, never pickled) has the sha256 recorded for the 2026-09-21 grid
    (GRID24_FULL_NET3_SHA256), and equals outputs/cache/grid24_full_Net3.pkl, max |diff| 0, when that git-ignored
    file predates this run.  If an earlier check had to build the cached file, only the hash comparison counts
    (the cached file would be the current code compared with itself)."""
    from .simgp import build_grid_24h, disk_preflight
    from .simulate import nominal_scenario
    free = disk_preflight([CACHE, tempfile.gettempdir()])
    sc = nominal_scenario("Net3")
    t0 = time.time()
    p1, Z1 = build_grid_24h(sc, "full", preflight_paths=[CACHE])
    secs = time.time() - t0
    fresh_hash = _grid_sha256(p1, Z1)
    assert fresh_hash == GRID24_FULL_NET3_SHA256, fresh_hash
    f = os.path.join(CACHE, "grid24_full_Net3.pkl")
    before = CACHE_BEFORE.get("grid24_full_Net3.pkl") if CACHE_BEFORE is not None else None
    out = {"sha256_matches_recorded": True, "members": len(p1),
           "cached_file_predates_run": before is not None if CACHE_BEFORE is not None else None}
    if os.path.exists(f):
        with open(f, "rb") as fh:
            p0, Z0 = pickle.load(fh)
        assert p0 == p1 and Z0.shape == Z1.shape and Z0.dtype == Z1.dtype
        d = float(np.abs(Z0 - Z1).max())
        assert d == 0.0 and _grid_sha256(p0, Z0) == GRID24_FULL_NET3_SHA256, d
        out.update(max_abs_diff_vs_cached_file=d, cached_file_mtime=time.strftime(
            "%Y-%m-%d %H:%M", time.localtime(os.stat(f).st_mtime)))
    return {**out, "_volatile": {"seconds": secs, "free_gb_before": free}}


@check("anchors", quick=False)
def tagged_grid_for_a_cold_month():
    """A non-default condition is cached under its own tagged name in outputs/cache, the default cached pickle is not
    touched, and a 10 C grid has at least as much chlorine everywhere as the 20 C one (decay-only grid, 75
    runs).  The tagged pickle is deleted afterwards if this check created it."""
    from .chemistry import Chemistry
    from .simgp import grid_cache_path, simulator_grid_24h
    from .simulate import nominal_scenario
    sc = nominal_scenario("Net3")
    cold = Chemistry(temp_C=10.0)
    f_def, f_cold = grid_cache_path(sc, CACHE, "decay"), grid_cache_path(sc, CACHE, "decay", cold)
    stat0 = os.stat(f_def) if os.path.exists(f_def) else None
    existed = os.path.exists(f_cold)
    try:
        pd_, Zd = simulator_grid_24h(sc, CACHE, "decay")
        pc, Zc = simulator_grid_24h(sc, CACHE, "decay", cond=cold)
        assert os.path.exists(f_cold) and os.path.basename(f_cold) != os.path.basename(f_def)
    finally:
        if not existed and os.path.exists(f_cold):
            os.remove(f_cold)
    if stat0 is not None:
        s1 = os.stat(f_def)
        assert (s1.st_mtime_ns, s1.st_size) == (stat0.st_mtime_ns, stat0.st_size), "the default cached pickle changed"
    assert pd_ == pc
    frac = float((Zc >= Zd - 1e-6).mean())
    assert frac == 1.0 and float((Zc - Zd).mean()) > 0, frac
    return {"tagged_file": os.path.basename(f_cold), "mean_ln_shift": float((Zc - Zd).mean())}


@check("anchors", quick=False)
def synthetic_pilot_net3():
    from .pilot import LEVELS, synthetic_log, validate
    ref = json.load(open(os.path.join(REPO, "outputs", "pilot", "validation_Net3.json")))
    log, inp = synthetic_log("Net3")
    _, summary = validate(inp, log, dose=1.2, threshold=0.2, window_months=3, cache_dir=CACHE)
    got = {}
    for subset, g in summary.groupby("taps"):
        w = g.n_test / g.n_test.sum()
        got[subset] = {"n": int(g.n_test.sum()), "rmse": float((g.rmse * w).sum()),
                       "rmse_persistence": float((g.rmse_persistence * w).sum()),
                       "rmse_network_mean": float((g.rmse_network_mean * w).sum()),
                       **{f"coverage{q}": float((g[f"coverage{q}"] * w).sum()) for q in LEVELS}}
    assert got == ref["mean"], (got, ref["mean"])
    assert round(got["all"]["rmse"], 3) == 0.077 and round(got["all"]["coverage90"], 3) == 0.923
    return {"rmse": got["all"]["rmse"], "coverage90": got["all"]["coverage90"]}


# Self-contained on purpose: --compare runs it against an older copy of the code (PYTHONPATH = that tree), which
# has no residualmap.checks to import.
APP_DEMO_SCRIPT = r'''
import json, re, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
at = AppTest.from_file(sys.argv[1], default_timeout=600).run()
assert not at.exception, [e.value for e in at.exception]
metric = [m for m in at.metric if m.label.startswith("Junctions likely below")][0]
[t for t in at.toggle if t.label.startswith("Show the true daily-minimum map")][0].set_value(True)
at.run()
assert not at.exception, [e.value for e in at.exception]
txt = [m.value for m in at.markdown if "True daily-minimum violations" in m.value][0]
m = re.search(r"\*\*(\d+) of (\d+)\*\* junctions.*found \*\*(\d+) of (\d+)\*\*.*with (\d+) false alarms", txt)
out = dict(zip(("violations", "junctions", "found", "unsampled_violations", "false_alarms"), map(int, m.groups())))
out["flagged"] = int(re.fullmatch(r"(\d+) of (\d+)", metric.value).group(1))
out["flagged_label"] = metric.label
out["sample_mean_mgL"] = float(re.search(r"Mean of your samples: ([0-9.]+) mg/L", txt).group(1))
out["reveal_text"] = txt
print("APP_DEMO_JSON " + json.dumps(out))
'''


def app_demo(code_root: str, cwd: str) -> dict:
    """The app's default demo (Net3, scenario 0, 8 samples) through a headless AppTest, in a subprocess that
    imports the code under `code_root` and runs with `cwd` as its working directory (which needs outputs/cache)."""
    r = subprocess.run([PY, "-c", APP_DEMO_SCRIPT, os.path.join(code_root, "app.py")], cwd=cwd, capture_output=True,
                       text=True, env={**os.environ, "PYTHONPATH": code_root})
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_DEMO_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-2000:]}")
    return json.loads(line[-1][len("APP_DEMO_JSON "):])


@check("anchors", quick=False)
def app_default_demo():
    """Headless AppTest of the app's default demo (Net3, scenario 0, 8 samples), as in the CHANGELOG's task-6
    line: the map flags 42 of 92 junctions at their daily minimum; the truth reveal reads 41 violations, 36 of 36
    unsampled found, 1 false alarm; the mean of the 8 samples is 0.53 mg/L."""
    got = app_demo(REPO, os.getcwd())
    num = {k: got[k] for k in APP_DEMO}
    assert num == APP_DEMO, num
    assert got["flagged_label"] == "Junctions likely below 0.2 mg/L (daily minimum)", got["flagged_label"]
    return num


# ----------------------------------------------------------------------------- runner
def run_checks(quick: bool) -> tuple[dict, dict]:
    """Returns (report, run_log): the report holds only results that do not change from run to run; the run log
    adds the date, the timings and the free disk."""
    global CACHE_BEFORE
    import warnings
    warnings.filterwarnings("ignore")      # GP optimiser bound warnings, as in experiment.py
    status0, CACHE_BEFORE = _git_status(), _cache_snapshot()
    results, timing, t_all = [], [], time.time()
    with scratch_cwd() as tmp:
        for name, group, is_quick, fn in CHECKS:
            if quick and not is_quick:
                results.append({"name": name, "group": group, "pass": None, "skipped": "--quick"})
                continue
            t0 = time.time()
            try:
                detail, ok, err = fn(), True, None
            except Exception as e:  # noqa: BLE001
                detail, ok, err = None, False, f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"
            vol = detail.pop("_volatile", {}) if isinstance(detail, dict) else {}
            results.append({"name": name, "group": group, "pass": ok,
                            "detail": _jsonable(detail), **({"error": err} if err else {})})
            timing.append({"name": name, "seconds": round(time.time() - t0, 2), **_jsonable(vol)})
            print(f"{'PASS' if ok else 'FAIL'}  {group:16s} {name}  ({time.time() - t0:.1f} s)", flush=True)
        left = sorted(os.listdir(tmp))
    status1, root1, cache1 = _git_status(), _root_scratch_files(), _cache_snapshot()
    # the default truth path writes temp.inp/.rpt/.bin to its cwd, as it always has: here that is the scratch cwd
    stray = [f for f in left if f not in ("outputs", "temp.inp", "temp.rpt", "temp.bin")]
    # outputs/cache is git-ignored, so it is watched here: no cached file may change or vanish, and the only new
    # files allowed are default (untagged) grids an earlier check needed and found missing
    changed = sorted(f for f in CACHE_BEFORE if cache1.get(f) != CACHE_BEFORE[f])
    new = sorted(set(cache1) - set(CACHE_BEFORE))
    bad_new = [f for f in new if TAGGED_CACHE.search(f) or not (f.startswith("grid24_") and f.endswith(".pkl"))]
    ok_repo = status0 == status1 and root1 == [] and stray == [] and changed == [] and bad_new == []
    results.append({"name": "repo_untouched", "group": "hygiene", "pass": ok_repo,
                    "detail": {"git_status_unchanged": status0 == status1, "msx_or_en_files_in_root": root1,
                               "unexpected_files_in_scratch_cwd": stray, "cache_files_changed_or_removed": changed,
                               "cache_files_added": new, "cache_files_added_not_allowed": bad_new}})
    print(f"{'PASS' if ok_repo else 'FAIL'}  {'hygiene':16s} repo_untouched", flush=True)
    import wntr
    import sklearn
    import scipy
    counts = {"all_pass": all(r["pass"] is not False for r in results),
              "n_pass": sum(r["pass"] is True for r in results), "n_fail": sum(r["pass"] is False for r in results),
              "n_skipped": sum(r["pass"] is None for r in results)}
    gen = "python -m residualmap.checks" + (" --quick" if quick else "")
    report = {"generated_by": gen,
              "about": "Run-invariant results only; the date, timings and free disk of each run are in "
                       "outputs/chem/checks_run.json (git-ignored).",
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                              "scipy": scipy.__version__, "sklearn": sklearn.__version__, "wntr": wntr.__version__},
              **counts, "checks": results}
    run_log = {"generated_by": gen, "date": time.strftime("%Y-%m-%d %H:%M"), "seconds": round(time.time() - t_all, 1),
               **counts, "timing": timing, "results": results}
    return report, run_log


# ----------------------------------------------------------------------------- baseline reproduction
RERUNS = {   # subdirectory -> commands, run in order with that directory as the cwd
    "net3": [["-m", "residualmap.experiment", "Net3", "8"],
             ["-m", "residualmap.experiment", "Net3", "8", "--structural=persistent"]],
    "net2": [["-m", "residualmap.experiment", "Net2", "8"]],
    "pilot": [["-m", "residualmap.pilot", "--synthetic", "Net3"]],
}
KEY_COLS = ("model", "route", "strategy", "n", "K", "seed", "held_out_month", "taps")


def _prepare_rerun_dir(d: str) -> None:
    os.makedirs(os.path.join(d, "outputs"), exist_ok=True)
    link = os.path.join(d, "outputs", "cache")
    if not os.path.islink(link):
        os.symlink(CACHE, link)


def _seed_structural_summary(d: str) -> None:
    # experiment.plot_stress reads <outdir>/summary_Net3.json BEFORE main() writes it (an ordering bug in the
    # committed code, kept as is); in a fresh directory the structural run would crash, so the committed
    # summary is put where the repo has it, which is what every committed run saw.  The run then overwrites it;
    # reproduce() checks that it did, so the seeded copy is never counted as a reproduced file.
    s = os.path.join(d, "outputs", "structural_persistent", "summary_Net3.json")
    if not os.path.exists(s):
        os.makedirs(os.path.dirname(s), exist_ok=True)
        shutil.copy(os.path.join(REPO, "outputs", "structural_persistent", "summary_Net3.json"), s)


def _output_mtimes(workdir: str) -> dict:
    """mtime_ns of every file the reruns can compare (cache excluded), keyed by '<sub>/<path under outputs>'."""
    out = {}
    for sub in RERUNS:
        root = os.path.join(workdir, sub, "outputs")
        for dirpath, dirnames, fnames in os.walk(root):
            dirnames[:] = [x for x in dirnames if x != "cache"]
            for fn in fnames:
                p = os.path.join(dirpath, fn)
                out[f"{sub}/{os.path.relpath(p, root)}"] = os.stat(p).st_mtime_ns
    return out


def reproduce(workdir: str) -> tuple[dict, list[str]]:
    """Runs RERUNS in workdir.  Returns the runs (cwd, seconds, return code, start and end time) and the files
    under workdir that no run wrote (seeded or left from an earlier run), which compare_outputs then excludes."""
    runs = {}
    env = {**os.environ, "PYTHONPATH": REPO}
    for sub in RERUNS:
        _prepare_rerun_dir(os.path.join(workdir, sub))
    _seed_structural_summary(os.path.join(workdir, "net3"))
    before = _output_mtimes(workdir)
    for sub, cmds in RERUNS.items():
        d = os.path.join(workdir, sub)
        for cmd in cmds:
            t0 = time.time()
            r = subprocess.run([PY, *cmd], cwd=d, env=env, capture_output=True, text=True)
            with open(os.path.join(d, f"log_{'_'.join(c.strip('-') for c in cmd[1:])}.txt"), "w") as fh:
                fh.write(r.stdout + r.stderr)
            runs[" ".join(cmd[1:])] = {"cwd": d, "started": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t0)),
                                       "ended": time.strftime("%Y-%m-%d %H:%M:%S"),
                                       "seconds": round(time.time() - t0, 1), "returncode": r.returncode}
            print(f"{' '.join(cmd[1:])}: {r.returncode} in {time.time() - t0:.0f} s", flush=True)
    after = _output_mtimes(workdir)
    not_written = sorted(k for k, t in after.items() if before.get(k) == t)
    return runs, not_written


def _csv_compare(a: str, b: str) -> dict:
    """Byte identity of the file, then row by row on the parsed values (floats read round-trip exactly, NaN
    equal to NaN): a row drifts if any of its values differs.  Drifting rows are grouped by model / route."""
    same_bytes = open(a, "rb").read() == open(b, "rb").read()
    da, db = pd.read_csv(a, float_precision="round_trip"), pd.read_csv(b, float_precision="round_trip")
    out = {"rows": len(db), "byte_identical": same_bytes}
    if list(da.columns) != list(db.columns) or len(da) != len(db):
        out.update(rows_identical=None, rows_drifting=None,
                   note=f"columns or row count differ ({len(da)} vs {len(db)} rows)")
        return out
    eq = (da == db) | (da.isna() & db.isna())
    bad = ~eq.all(axis=1)
    out.update(rows_identical=int((~bad).sum()), rows_drifting=int(bad.sum()))
    if bad.any():
        group = next((c for c in ("model", "route", "taps") if c in db.columns), None)
        num = [c for c in db.columns if pd.api.types.is_numeric_dtype(db[c]) and c not in KEY_COLS]
        drift = {}
        for g, rows in (db[bad].groupby(group).groups.items() if group else [("all", db.index[bad])]):
            diff = (da.loc[rows, num] - db.loc[rows, num]).abs().max()
            drift[str(g)] = {"rows": len(rows), "max_abs_diff": {c: float(v) for c, v in diff.items() if v > 0}}
        out["drift_by_group"] = drift
    return out


def _json_diff(a, b, path=""):
    if isinstance(b, dict) and isinstance(a, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append({"path": f"{path}/{k}", "missing_in": "rerun" if k not in a else "committed"})
            else:
                out += _json_diff(a[k], b[k], f"{path}/{k}")
        return out
    if isinstance(b, list) and isinstance(a, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += _json_diff(x, y, f"{path}[{i}]")
        return out
    if a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b)):
        return []
    d = {"path": path}
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        d["abs_diff"] = abs(float(a) - float(b))
    return [d]


def compare_outputs(workdir: str, reference: str | None = None, exclude=()) -> dict:
    """Every file the reruns wrote (except the cache) against the same path under the repo's outputs/ (or under
    another rerun's directories when `reference` is given).  Files in `exclude` ('<sub>/<path>' keys from
    reproduce: present in the workdir but not written by its runs) are listed and not compared."""
    files, excluded = {}, []
    for sub in RERUNS:
        root = os.path.join(workdir, sub, "outputs")
        for dirpath, dirnames, fnames in os.walk(root):
            dirnames[:] = [x for x in dirnames if x != "cache"]
            for fn in sorted(fnames):
                p = os.path.join(dirpath, fn)
                rel = os.path.relpath(p, root)
                if f"{sub}/{rel}" in exclude:
                    excluded.append(rel)
                    continue
                ref = os.path.join(reference, sub, "outputs", rel) if reference else os.path.join(REPO, "outputs", rel)
                if not os.path.exists(ref):
                    files[rel] = {"status": "not in reference"}
                    continue
                if fn.endswith(".csv"):
                    files[rel] = _csv_compare(p, ref)
                elif fn.endswith(".json"):
                    same = open(p, "rb").read() == open(ref, "rb").read()
                    files[rel] = {"byte_identical": same}
                    if not same:
                        files[rel]["differences"] = _json_diff(json.load(open(p)), json.load(open(ref)))
                else:
                    files[rel] = {"byte_identical": open(p, "rb").read() == open(ref, "rb").read()}
    drifting = sorted(f"{rel}::{g}" for rel, v in files.items() for g in v.get("drift_by_group", {}))
    not_identical = sorted(rel for rel, v in files.items() if v.get("byte_identical") is False)
    committed = set(subprocess.run(["git", "-C", REPO, "ls-files", "outputs"], capture_output=True, text=True).stdout.split())
    produced = {os.path.join("outputs", r) for r in files}
    return {"files": files, "n_files": len(files), "n_byte_identical": sum(v.get("byte_identical") is True for v in files.values()),
            "not_byte_identical": not_identical, "drifting_rows": drifting, "not_compared_not_written_by_run": sorted(excluded),
            "committed_outputs_not_regenerated": sorted(c for c in committed - produced if not reference)}


def _write_repro(label: str, entry: dict) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    data = json.load(open(REPRO)) if os.path.exists(REPRO) else {}
    data["about"] = ("Reruns of the committed experiments from a scratch working directory, compared file by file "
                     "(and CSV row by row) with the committed outputs/, plus the app's default demo. Written by "
                     "python -m residualmap.checks --compare/--reproduce. The workdir paths are scratch directories of "
                     "the session that made each entry and do not persist; the comparison results are the record.")
    prev = data.setdefault("runs", {}).get(label, {})
    if entry.get("runs") is None and prev.get("workdir") == entry["workdir"] and prev.get("runs"):
        entry["runs"] = prev["runs"]          # a re-comparison of the same rerun keeps its timings
    data["runs"][label] = entry
    pre = data["runs"].get("pre_edit", {}).get("vs_committed", {})
    data["drifting_rows_pre_recorded"] = pre.get("drifting_rows", [])
    data["not_byte_identical_pre_recorded"] = pre.get("not_byte_identical", [])
    with open(REPRO, "w") as fh:
        json.dump(_jsonable(data), fh, indent=1)


def _fmt_ns(t_ns: int) -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t_ns / 1e9))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="skip the fresh grid, the synthetic pilot and the app")
    ap.add_argument("--reproduce", metavar="WORKDIR", help="rerun Net3 8, Net3 8 --structural=persistent, Net2 8 and "
                                                           "the synthetic pilot in WORKDIR, then compare")
    ap.add_argument("--compare", metavar="WORKDIR", help="compare an existing rerun directory (net3/, net2/, pilot/)")
    ap.add_argument("--reference", metavar="WORKDIR", help="also compare against another rerun directory")
    ap.add_argument("--label", default="rerun")
    ap.add_argument("--code", default=None, help="what code the rerun used (default: git HEAD plus local changes)")
    ap.add_argument("--code-root", default=None, help="the tree whose app.py and residualmap/ the app demo runs "
                                                      "(default: this repo)")
    ap.add_argument("--no-app", action="store_true", help="do not run the app's default demo")
    ap.add_argument("--note", default=None, help="free text stored with the entry (how the rerun was made)")
    a = ap.parse_args()
    if a.reproduce or a.compare:
        workdir = os.path.abspath(a.reproduce or a.compare)
        runs, not_written = reproduce(workdir) if a.reproduce else (None, [])
        head = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", REPO, "status", "--porcelain", "residualmap", "app.py"], capture_output=True, text=True).stdout.split("\n")
        written = [t for k, t in _output_mtimes(workdir).items() if k not in not_written]
        entry = {"entry_written": time.strftime("%Y-%m-%d %H:%M"),
                 "outputs_written_between": f"{_fmt_ns(min(written))} and {_fmt_ns(max(written))} (file times)" if written else None,
                 "workdir": workdir,
                 "code": a.code or f"git {head}" + (f" plus local changes to {len([d for d in dirty if d])} files" if any(dirty) else ""),
                 "commands": {sub: [" ".join(c[1:]) for c in cmds] for sub, cmds in RERUNS.items()},
                 "runs": runs, "all_runs_ok": None if runs is None else all(r["returncode"] == 0 for r in runs.values()),
                 "not_written_by_run": not_written if runs is not None else "unknown (--compare of an existing rerun)",
                 **({"note": a.note} if a.note else {}), "vs_committed": compare_outputs(workdir, exclude=set(not_written))}
        if a.reference:
            entry["vs_reference"] = {"reference": os.path.abspath(a.reference),
                                     **compare_outputs(workdir, os.path.abspath(a.reference), exclude=set(not_written))}
        app_ok = True
        if not a.no_app:
            app_dir = os.path.join(workdir, "app")
            _prepare_rerun_dir(app_dir)
            root = os.path.abspath(a.code_root or REPO)
            try:
                got = app_demo(root, app_dir)
                app_ok = {k: got[k] for k in APP_DEMO} == APP_DEMO
                entry["app_default_demo"] = {"code_root": root, "run_at": time.strftime("%Y-%m-%d %H:%M"),
                                             "matches_changelog": app_ok, **got}
            except Exception as e:  # noqa: BLE001
                app_ok = False
                entry["app_default_demo"] = {"code_root": root, "error": f"{type(e).__name__}: {e}"}
        _write_repro(a.label, entry)
        v = entry["vs_committed"]
        print(f"{a.label}: {v['n_byte_identical']} of {v['n_files']} files byte-identical to the committed outputs; "
              f"drifting rows: {v['drifting_rows'] or 'none'}; not identical: {v['not_byte_identical'] or 'none'}; "
              f"not written by the run: {not_written or 'none'}")
        if a.reference:
            r = entry["vs_reference"]
            print(f"vs reference: {r['n_byte_identical']} of {r['n_files']} byte-identical; not identical: {r['not_byte_identical'] or 'none'}")
        if "app_default_demo" in entry:
            print("app default demo:", {k: entry["app_default_demo"].get(k) for k in (*APP_DEMO, "error") if k in entry["app_default_demo"]})
        ok = entry["all_runs_ok"] is not False and not not_written and app_ok
        return 0 if ok else 1
    report, run_log = run_checks(a.quick)
    os.makedirs(OUT_DIR, exist_ok=True)
    if not a.quick:   # the committed report is written by full runs only
        old = open(REPORT, "rb").read() if os.path.exists(REPORT) else None
        new = (json.dumps(_jsonable(report), indent=1) + "\n").encode()
        with open(REPORT, "wb") as fh:
            fh.write(new)
        run_log["report_changed"] = old != new
    with open(RUN_LOG, "w") as fh:
        json.dump(_jsonable(run_log), fh, indent=1)
    dest = RUN_LOG if a.quick else REPORT
    print(f"\n{report['n_pass']} passed, {report['n_fail']} failed, {report['n_skipped']} skipped in {run_log['seconds']} s "
          f"-> {os.path.relpath(dest, REPO)}" + ("" if a.quick else f" ({'changed' if run_log['report_changed'] else 'unchanged'})"))
    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
