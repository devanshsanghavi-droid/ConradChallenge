"""
checks.py: the saved checks for the chemistry work (iteration 4, journal tasks 8 to 14).

    python -m residualmap.checks            # every check, about 400 s on this machine (63 checks, after task 13)
    python -m residualmap.checks --quick    # skips the fresh grid, the synthetic pilot, the app and the slow
                                            # seasonal, organics and chloramine (MSX) checks (about 24 s)

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
def chloramine_built_and_kept_separate():
    """Task 12 replaced task 8's refusal: chloramine has its own physics now, and is kept apart from free chlorine.  A
    chloramine grid condition is first order at 20 C with neutral keywords; a seasonal or TOC-scaled chloramine model
    or truth is refused (NotImplementedError), and so is 'epa_msx' for free chlorine (ValueError); the chloramine grid
    runs only under a chloramine condition and a chloramine condition only on it (ValueError, before any cache is read);
    its cache name carries its own tag.  The truth's hidden draws keep task 8's order: E_true and theta_w first, then the
    four chloramine uniforms, from default_rng(40_000 + seed); free chlorine's draws are unchanged."""
    from .chemistry import CHLORAMINE, Chemistry
    from .simgp import GRIDS, grid_cache_path, simulator_grid_24h
    from .simulate import build_scenario, hidden_chem_draws, nominal_scenario
    ca = Chemistry(disinfectant=CHLORAMINE)
    assert ca.sim_kwargs() == {"kb_scale": 1.0, "kw_scale": 1.0, "temp_C": None}
    assert Chemistry(disinfectant=CHLORAMINE, temp_C=20.0).sim_kwargs() == {"kb_scale": 1.0, "kw_scale": 1.0, "temp_C": 20.0}
    refused = {}
    sc = nominal_scenario("Net3")
    before = sorted(os.listdir(CACHE))
    for name, call, exc in (
            ("seasonal_grid", Chemistry(disinfectant=CHLORAMINE, temp_C=15.0).sim_kwargs, NotImplementedError),
            ("seasonal_truth", lambda: build_scenario("Net3", 3, chem=Chemistry(disinfectant=CHLORAMINE, temp_C=15.0)), NotImplementedError),
            ("toc_input", Chemistry(disinfectant=CHLORAMINE, toc_mgL=2.0).sim_kwargs, NotImplementedError),
            ("clark_chloramine", lambda: Chemistry(disinfectant=CHLORAMINE, kinetics="clark", toc_mgL=2.0, phi=0.85).require_built(), NotImplementedError),
            ("epa_msx_free", Chemistry(kinetics="epa_msx").require_built, ValueError),
            ("epa_msx_grid", Chemistry(disinfectant=CHLORAMINE, kinetics="epa_msx").sim_kwargs, ValueError),
            ("ca_on_free_grid", lambda: simulator_grid_24h(sc, CACHE, "decay", cond=ca), ValueError),
            ("ca_grid_without_ca", lambda: simulator_grid_24h(sc, CACHE, "chloramine"), ValueError),
            ("loss_split_chloramine", lambda: build_scenario("Net3", 3, chem=ca, truth_loss_split=True), NotImplementedError)):
        try:
            call()
        except exc:
            refused[name] = True
        else:
            refused[name] = False
    assert all(refused.values()), refused
    assert sorted(os.listdir(CACHE)) == before, "a refused call wrote a cache file"
    name = os.path.basename(grid_cache_path(nominal_scenario("Net3", 14, 2.0), CACHE, "chloramine", ca))
    assert name.startswith("grid24_chloramine_Net3_d2_ca_") and len(GRIDS["chloramine"][0]) == 8, name
    d_ca, d_fc = hidden_chem_draws(CHLORAMINE, 3), hidden_chem_draws("free_chlorine", 3)
    assert d_ca["rng_seed"] == 40_003 and d_fc["rng_seed"] == 20_003 and d_ca["E_true_K"] != d_fc["E_true_K"]
    rc = np.random.default_rng(40_003)
    e, t = rc.uniform(4660.0, 12104.0), rc.uniform(1.00, 1.07)
    assert (d_ca["E_true_K"], d_ca["theta_w"]) == (float(e), float(t)) and set(d_fc) == {"rng_seed", "E_true_K", "theta_w"}
    assert [d_ca[k] for k in ("u_pH", "u_cl2n", "u_toc", "u_alk")] == [float(rc.random()) for _ in range(4)]
    fc = build_scenario("Net3", 3, chem=Chemistry(temp_C=15.0))
    assert fc.chem["rng_seed"] == 20_003 and fc.chem["E_true_K"] == d_fc["E_true_K"] and fc.chem["species"] == "free chlorine"
    return {"refused": refused, "chloramine_grid_cache_name": name, "chloramine_rng_seed": d_ca["rng_seed"]}


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


# ----------------------------------------------------------------------------- water age and the loss split (task 9)
@check("water_age")
def age_band_members_and_nominal_identity():
    """hydraulic_age_band runs the grid's 9 hydraulic members in the grid's own order; its (1, 1) member is the nominal
    model's age bit for bit (Net3, Net2); min <= median <= max at every junction and hour; no age beyond the run.
    Its initial-water share (the AGE run with the starting ages raised) agrees with an independent measure, one minus
    a no-decay chlorine run with every source at 1 and the starting water at 0, within 0.001 (Net3, Net2)."""
    import itertools
    from .age import hydraulic_age_band, hydraulic_members
    import wntr
    from .simgp import GRIDS, HOURS
    from .simulate import DURATION_DAYS, _last_day, load, nominal_scenario, source_nodes
    want = [tuple(map(float, p[3:])) for p in itertools.product(*GRIDS["full"])][:9]
    assert hydraulic_members() == want, (hydraulic_members(), want)
    out = {}
    for net in ("Net3", "Net2"):
        sc = nominal_scenario(net)
        band = hydraulic_age_band(sc)
        assert band.ages.shape == (9, len(HOURS), len(sc.junctions))
        assert np.array_equal(band.nominal.values, sc.age_by_hour_h.loc[HOURS, sc.junctions].values), net
        lo, med, hi = band.lo().values, band.med().values, band.hi().values
        assert (lo <= med).all() and (med <= hi).all()
        assert band.ages.min() >= 0.0 and band.ages.max() <= DURATION_DAYS * 24 + 24, band.ages.max()
        w = (band.hi() - band.lo()).max()
        sh = band.initial_share
        assert sh.shape == (len(HOURS), len(sc.junctions)) and float(sh.values.min()) >= 0 and float(sh.values.max()) <= 1
        with tempfile.TemporaryDirectory(prefix="rm_check_") as tmp:
            wn = load(net)
            wn.options.quality.parameter = "CHEMICAL"
            wn.options.reaction.bulk_coeff = wn.options.reaction.wall_coeff = 0.0
            for _, pipe in wn.pipes():
                pipe.wall_coeff = 0.0
            for n, node in wn.nodes():
                node.initial_quality = 0.0
            for res in source_nodes(wn):
                wn.add_source(f"src_{res}", res, "CONCEN", 1.0)
            prefix = os.path.join(tmp, "nd")
            c = _last_day(wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=prefix).node["quality"], wn.junction_name_list)
        dev = float((sh - (1.0 - c.loc[HOURS, sc.junctions])).abs().values.max())
        assert dev < 1e-3, (net, dev)
        out[net] = {"members": len(band.members), "max_range_width_h": float(w.max()), "median_range_width_h": float(w.median()),
                    "initial_share_max_abs_diff_vs_no_decay_run": dev,
                    "n_lower_bound_daily_mean": int(band.lower_bound("mean").sum()),
                    "n_lower_bound_daily_max": int(band.lower_bound("max").sum())}
    return out


@check("water_age")
def posterior_age_is_the_weighted_member_age():
    """hydraulic_weights sums the joint (member x dose) posterior over the dose and decay axes; posterior_age is that
    weighted sum of the 9 members' ages; one-hot weights return the member itself with a zero-width band; uniform
    weights return the members' mean; weighted_quantile is the smallest value reaching the cumulative weight."""
    from types import SimpleNamespace
    from .age import hydraulic_age_band, hydraulic_weights, posterior_age, weighted_quantile
    sc, X, S, m = _fitted_net3()
    band = hydraulic_age_band(sc)
    w = hydraulic_weights(m)
    assert w.shape == (9,) and abs(w.sum() - 1.0) < 1e-12
    assert np.allclose(w, m.W_.sum(axis=1).reshape(-1, 9).sum(axis=0), rtol=0, atol=1e-15)
    pa = posterior_age(m, band)
    assert np.allclose(pa.by_hour.values, np.tensordot(w, band.ages, axes=1), rtol=0, atol=1e-12)
    assert (pa.lo90 <= pa.daily_mean + 1e-9).all() and (pa.daily_mean <= pa.hi90 + 1e-9).all()
    for k in (0, 4, 8):
        W = np.zeros_like(m.W_); W[np.arange(k, len(m.params), 9), 2] = 1.0 / (len(m.params) // 9)
        one = posterior_age(SimpleNamespace(params=m.params, W_=W, sc=sc), band)
        assert np.allclose(one.by_hour.values, band.ages[k], rtol=0, atol=1e-12)
        assert np.allclose(one.lo90.values, band.ages[k].mean(axis=0)) and np.allclose(one.hi90.values, band.ages[k].mean(axis=0))
    uni = posterior_age(SimpleNamespace(params=m.params, W_=np.full_like(m.W_, 1.0 / m.W_.size), sc=sc), band)
    assert np.allclose(uni.by_hour.values, band.ages.mean(axis=0), rtol=0, atol=1e-12)
    v = np.array([[3.0], [1.0], [2.0]])
    assert weighted_quantile(v, np.array([0.2, 0.5, 0.3]), 0.05)[0] == 1.0
    assert weighted_quantile(v, np.array([0.2, 0.5, 0.3]), 0.6)[0] == 2.0
    assert weighted_quantile(v, np.array([0.2, 0.5, 0.3]), 0.95)[0] == 3.0
    try:
        hydraulic_weights(SimpleNamespace(params=m.params[:75], W_=m.W_[:75]))
    except ValueError:
        pass
    else:
        raise AssertionError("hydraulic_weights accepted a grid without the full hydraulic block in order")
    return {"posterior_weights": [round(float(x), 4) for x in w]}


@check("water_age")
def loss_split_exact_on_plug_flow_chain():
    """On the 12-pipe plug-flow chain with first-order bulk (0.8 /day) and wall (0.5 m/day) decay, the four-run split
    is exact: L_bulk / travel time = kb, L_wall / travel time is the same at every junction (one pipe size, one flow),
    and (L_bulk + L_wall) / L_tot = 1, each within 0.1%."""
    from .age import split_from_runs
    from .chemistry import DAY, set_bulk_kinetics, use_mg_per_litre
    from .simulate import _run_quality
    kb, kw, runs = 0.8, 0.5, {}
    for tag, b_on, w_on in (("full", 1, 1), ("bulk_only", 1, 0), ("wall_only", 0, 1), ("no_decay", 0, 0)):
        wn, t = _chain()
        scale = use_mg_per_litre(wn, {"R": 1.5})
        set_bulk_kinetics(wn, 1, -kb * b_on)
        for _, p in wn.pipes():
            p.wall_coeff = -kw * w_on / DAY
        runs[tag] = _run_quality(wn).iloc[[-1]][wn.junction_name_list] * scale
    s = split_from_runs(runs)
    kb_hat, kw_eff = s.L_bulk.values / t, s.L_wall.values / t
    assert np.abs(kb_hat / kb - 1).max() < 1e-3, kb_hat
    assert np.abs(kw_eff / kw_eff.mean() - 1).max() < 1e-3, kw_eff
    assert np.abs(s.nonadditivity.values - 1).max() < 1e-3, s.nonadditivity.values
    assert np.allclose(s.wall_share.values, kw_eff.mean() / (kb + kw_eff.mean()), atol=1e-3)
    return {"kb_recovered_per_day": float(kb_hat.mean()), "kw_effective_per_day": float(kw_eff.mean()),
            "wall_share": float(s.wall_share.mean()), "max_abs_nonadditivity_minus_1": float(np.abs(s.nonadditivity - 1).max())}


@check("water_age")
def truth_loss_split_leaves_the_truth_unchanged():
    """build_scenario(truth_loss_split=True) adds three runs after all draws: the chlorine truth and the truth's age are
    bit-identical to the runs without it (every coefficient is restored), on the default and the chemistry path; the
    corrected-unit truth (results x1000) gives the same split as the legacy one; the no-decay run never exceeds the
    highest source dose; the split adds up to within 5% on Net3 seed 0.  The no-decay run falls well below the doses
    where water from the file's initial tank contents (quality 0) is still arriving after 6 days; that is why the
    split's reference is this run and not the source dose (the dose would count that dilution as chlorine loss), and
    the number of such junctions is recorded."""
    from .age import truth_loss_split
    from .chemistry import Chemistry
    from .simulate import build_scenario
    ref = build_scenario("Net3", 0, truth_age=True)
    sc = build_scenario("Net3", 0, truth_age=True, truth_loss_split=True)
    assert np.array_equal(ref.truth_by_hour.values, sc.truth_by_hour.values)
    assert np.array_equal(ref.truth_age_by_hour_h.values, sc.truth_age_by_hour_h.values)
    assert ref.truth_loss_runs is None
    ts = truth_loss_split(sc)
    nd = sc.truth_loss_runs["no_decay"]
    assert 0.0 <= float(nd.values.min()) and float(nd.values.max()) <= 1.2 * 1.1 + 1e-6, (nd.values.min(), nd.values.max())
    n_unflushed = int((nd.mean() < 0.95 * 1.2 * 0.9).sum())
    na = ts.nonadditivity.dropna()
    assert na.between(0.95, 1.05).all(), (na.min(), na.max())
    assert np.allclose(sc.truth_loss_runs["full"].clip(lower=0).values, sc.truth_by_hour.values, rtol=0, atol=0)
    chem = Chemistry(temp_C=12.5)
    c0, c1 = build_scenario("Net3", 0, chem=chem), build_scenario("Net3", 0, chem=chem, truth_loss_split=True)
    assert np.array_equal(c0.truth_by_hour.values, c1.truth_by_hour.values) and c1.chem["quality_scale"] == 1.0
    si = build_scenario("Net3", 0, chem=Chemistry(kinetics="first_si"), truth_loss_split=True)
    assert si.chem["quality_scale"] == 1000.0
    d = float((truth_loss_split(si).wall_share - ts.wall_share).abs().max())
    assert d < 1e-3, d
    return {"wall_share_median": float(ts.wall_share.median()), "nonadditivity_range": [float(na.min()), float(na.max())],
            "corrected_units_max_abs_diff_wall_share": d, "no_decay_min_mgL": float(nd.values.min()),
            "junctions_with_unflushed_initial_water": n_unflushed}


@check("water_age")
def route_reason_default_unchanged_age_opt_in():
    """plan_route's default reasons are the committed text (Net3 scenario 0, K = 8, from the prior, equals
    route_scenario0_K8 in outputs/summary_Net3.json), so the committed experiments reproduce; with_age=True picks the
    same sites, hours and scores, adds the nominal water age at that junction and hour ('at least' where the
    starting water's share is over 5%), and writes no em dash.  At another threshold (0.5) the route's P is the
    model's own P(daily min < 0.5), whether the model was built for 0.5 or for 0.2."""
    from .features import build_features
    from .route import plan_route
    from .simgp import SimGP24
    from .simulate import nominal_scenario
    sc = nominal_scenario("Net3")
    prior = SimGP24(sc, build_features(sc), seed=0, cache_dir=CACHE).fit_prior()
    base = plan_route(prior, 8)
    committed = json.load(open(os.path.join(REPO, "outputs", "summary_Net3.json")))["route_scenario0_K8"]
    assert base.to_dict("records") == committed
    from .age import INITIAL_SHARE_MAX, hydraulic_age_band
    share = hydraulic_age_band(sc).initial_share
    aged = plan_route(prior, 8, with_age=True, age_initial_share=share)
    for c in ("junction", "hour", "p_below", "score"):
        assert aged[c].tolist() == base[c].tolist(), c
    n_at_least = 0
    for j, h, a, r, w in zip(aged.junction, aged.hour, aged.water_age_h, aged.reason, aged.why):
        age = float(sc.age_by_hour_h.loc[h, j])
        word = "at least" if share.loc[h, j] > INITIAL_SHARE_MAX else "about"
        n_at_least += word == "at least"
        assert a == round(age, 1) and f"sample at {h:02d}:00" in r and f"water {word} {age:.0f} h old then" in r, r
        assert r.endswith(w) and "\u2014" not in r and "\u2013" not in r, r   # no em or en dash
    m5 = SimGP24(sc, build_features(sc), seed=0, cache_dir=CACHE, threshold=0.5).fit_prior()
    r5, r5b = plan_route(m5, 8, threshold=0.5), plan_route(prior, 8, threshold=0.5)
    d5 = m5.predict_daily_min()
    assert r5.p_below.tolist() == r5b.p_below.tolist() == [round(float(d5.loc[j, "p_below"]), 2) for j in r5.junction]
    return {"first_reason_with_age": aged.reason.iloc[0], "n_at_least": n_at_least,
            "p_below_at_0.5": r5.p_below.tolist()}


# ----------------------------------------------------------------------------- temperature (task 10)
@check("seasonal")
def seasonal_bank_block_order():
    """The temperature bank keeps every hypothesis block in the grid's own member order.  Full-grid layout (no EPANET
    runs): the stacked parameters pass check_grid_order, a stack shifted by one member fails it, the hypothesis index
    is block by block, and local_hydraulic_var of a stack equals the sum of its blocks' (it is linear in the weights).
    On the 75-member decay grid, built in memory (cache 'off', so no file is added): H0 is the committed grid; at 20 C
    every E/R block IS that grid; a 10 C block equals a direct build_grid_24h of the same condition bit for bit; and at
    10 C every block has at least the 20 C chlorine, more for a larger E/R (smaller rate below 20 C)."""
    import itertools
    from .chemistry import ER_HYPOTHESES_K, Chemistry
    from .seasonal import Bank, covariate_bank
    from .simgp import GRIDS, build_grid_24h, check_grid_order, local_hydraulic_var, n_hydraulic, simulator_grid_24h
    from .simulate import nominal_scenario
    params = list(itertools.product(*GRIDS["full"]))
    n_hyd, nb = n_hydraulic("full"), len(params)
    layout = Bank(params, np.zeros((nb, 1, 1), np.float32), {}, tuple(ER_HYPOTHESES_K), "arrhenius", "full")
    st = layout.stacked_params
    check_grid_order(st, n_hyd)
    assert len(st) == 4 * nb and np.array_equal(layout.hyp_index, np.repeat(np.arange(4), nb))
    try:
        check_grid_order(st[1:] + st[:1], n_hyd)
    except ValueError:
        pass
    else:
        raise AssertionError("check_grid_order accepted a stack shifted by one member")
    rng = np.random.default_rng(0)
    Zs, w = rng.normal(size=(4 * nb, 2, 3)), rng.random(4 * nb)
    w /= w.sum()
    whole = local_hydraulic_var(w, Zs, n_hyd)
    parts = sum(local_hydraulic_var(w[b * nb:(b + 1) * nb], Zs[b * nb:(b + 1) * nb], n_hyd) for b in range(4))
    assert np.allclose(whole, parts, rtol=0, atol=1e-12)
    sc = nominal_scenario("Net3")
    bank = covariate_bank(sc, [10.0, 20.0], grid="decay", cache_dir=CACHE, cache="off")
    p0, Z0 = simulator_grid_24h(sc, CACHE, "decay")
    assert bank.params == list(p0) and np.array_equal(bank.h0, Z0) and bank.temps == [10.0, 20.0]
    assert all(np.array_equal(bank.blocks[(20.0, e)], Z0) for e in bank.ers)
    direct = build_grid_24h(sc, "decay", cond=Chemistry(temp_C=10.0, er_K=8000.0))[1]
    assert np.array_equal(bank.blocks[(10.0, 8000.0)], direct)
    S = bank.stack(10.0)
    assert S.shape == (4 * len(p0),) + Z0.shape[1:]
    assert all(np.array_equal(S[b * len(p0):(b + 1) * len(p0)], blk) for b, blk in enumerate(bank.blocks_at(10.0)))
    z5, z8, z12 = (bank.blocks[(10.0, e)] for e in bank.ers)
    gaps = {"E5000_vs_20C": float((z5 - Z0).min()), "E8000_vs_E5000": float((z8 - z5).min()),
            "E12000_vs_E8000": float((z12 - z8).min())}
    assert all(v >= -1e-6 for v in gaps.values()), gaps
    try:
        bank.stack(15.0)
    except KeyError:
        pass
    else:
        raise AssertionError("a bank served a temperature it does not hold")
    return {"stacked_members_full": len(st), "min_ln_gap_at_10C": gaps,
            "mean_ln_gain_10C_E8000_vs_20C": float((z8 - Z0).mean())}


@check("seasonal")
def seasonal_model_nests_simgp24():
    """At 20 C every hypothesis block is the committed grid, so the temperature-aware model must be today's model:
    on Net3 scenario 0 with 8 samples, its posterior is SimGP24's split evenly over the 4 blocks, its MAP member and
    hourly predictions are SimGP24's (to 1e-10), and kb20 is SimGP24's posterior-mean kb.  The daily minimum agrees to
    Monte-Carlo noise only (the draws are taken over 4 times as many entries).  With identical blocks the samples cannot
    tell the hypotheses apart, so under the 'cejas' and 'h0_0.1' priors the posterior over hypotheses equals the prior.
    The plant dose enters as an exact ln-offset: readings doubled with a dose ratio of 2 give the same posterior and
    predictions shifted by exactly ln 2 (readings above the 0.02 mg/L floor)."""
    from .seasonal import SeasonalSimGP24, covariate_bank, hypothesis_prior
    sc, X, S, m = _fitted_net3()
    bank = covariate_bank(sc, [20.0], cache_dir=CACHE, cache="off")
    S20 = S.assign(temp_C=20.0)
    ms = SeasonalSimGP24(sc, X, bank, seed=0, cache_dir=CACHE).fit(S20, target_temp_C=20.0)
    dW = max(float(np.abs(ms.W_.reshape(4, len(m.params), -1)[b] - m.W_ / 4).max()) for b in range(4))
    assert dW < 1e-15, dW
    assert ms.map_params_ == m.map_params_ and ms.map_dose_ == m.map_dose_
    a, b = m.predict_hours(), ms.predict_hours()
    d_mu, d_sd = float(np.abs(a[0] - b[0]).max()), float(np.abs(a[1] - b[1]).max())
    assert d_mu < 1e-10 and d_sd < 1e-10, (d_mu, d_sd)
    P = np.asarray(m.params, dtype=float)
    assert abs(ms.kb20() - float(m.w_ @ P[:, 0])) < 1e-12
    pa, pb = m.predict_daily_min(), ms.predict_daily_min()
    dp = (pa.p_below - pb.p_below).abs()
    assert float(dp.mean()) < 0.01 and int((pa.p_below > 0.5).sum()) == int((pb.p_below > 0.5).sum())
    for pr in ("cejas", "h0_0.1"):
        mp = SeasonalSimGP24(sc, X, bank, prior=pr, seed=0, cache_dir=CACHE).fit(S20, target_temp_C=20.0)
        assert np.allclose(list(mp.hypothesis_posterior().values()), hypothesis_prior(pr), rtol=0, atol=1e-12), pr
        assert float(np.abs(mp.predict_hours()[0] - b[0]).max()) < 1e-10, pr
    Sa = S20[S20.y >= 0.02]
    m1 = SeasonalSimGP24(sc, X, bank, seed=0, cache_dir=CACHE).fit(Sa, target_temp_C=20.0)
    m2 = SeasonalSimGP24(sc, X, bank, seed=0, cache_dir=CACHE).fit(Sa.assign(y=Sa.y * 2, dose_ratio=2.0),
                                                                    target_temp_C=20.0, target_dose_ratio=2.0)
    h1, h2 = m1.predict_hours(), m2.predict_hours()
    d_dose = float(np.abs(h2[0] - h1[0] - np.log(2.0)).max())
    assert float(np.abs(m2.W_ - m1.W_).max()) < 1e-15 and d_dose < 1e-8 and float(np.abs(h2[1] - h1[1]).max()) < 1e-8, d_dose
    return {"max_abs_diff_W_block_vs_W0_over_4": dW, "max_abs_diff_hourly_z_mu": d_mu, "max_abs_diff_hourly_z_sd": d_sd,
            "daily_min_p_below_mean_abs_diff_mc": float(dp.mean()), "n_flagged": int((pb.p_below > 0.5).sum()),
            "kb20": ms.kb20(), "dose_offset_max_abs_err": d_dose}


@check("seasonal")
def seasonal_truth_warming_and_log():
    """The truth side of task 10.  In-network warming with the soil at the plant temperature is the truth without it,
    bit for bit (Net3 and Net2; so warming makes no draw); with warmer soil every junction has at most the chlorine of
    that truth, pipe temperatures lie between the plant's and the soil's, and the hidden E/R and theta_w are the seed's
    own.  Warming without a plant temperature is refused.  The seasonal synthetic log uses month_seed = 1000 + 100 seed
    + m (distinct for seeds 0 to 15 and months 0 to 11, where the committed 100 + 10 seed + m repeats), logs each month's
    plant temperature with every reading, and returns the plant log."""
    from .chemistry import Chemistry, Warming
    from .pilot import synthetic_log
    from .seasonal import plant_schedule
    from .simulate import build_scenario, hidden_chem_draws
    out = {}
    for net in ("Net3", "Net2"):
        kw = TRUTHS[net]
        a = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=12.5), **kw)
        b = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=12.5), warming=Warming(soil_temp_C=12.5), **kw)
        c = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=12.5), warming=Warming(soil_temp_C=19.0), **kw)
        assert np.array_equal(a.truth_by_hour.values, b.truth_by_hour.values), net
        assert (c.truth_by_hour.values <= a.truth_by_hour.values + 1e-9).all(), net
        w = c.chem["warming"]
        assert 12.5 <= w["pipe_temp_C_min"] <= w["pipe_temp_C_max"] <= 19.0 and 12.5 <= w["mixed_temp_C"] <= 19.0
        hd = hidden_chem_draws("free_chlorine", 0)
        assert c.chem["E_true_K"] == hd["E_true_K"] and c.chem["theta_w"] == hd["theta_w"]
        out[net] = {"mean_shift_mgL_soil_19C": float((c.truth_by_hour - a.truth_by_hour).values.mean()),
                    "mixed_temp_C": w["mixed_temp_C"], "pipe_temp_C_median": w["pipe_temp_C_median"]}
    try:
        build_scenario("Net3", 0, warming=Warming(soil_temp_C=15.0))
    except ValueError:
        pass
    else:
        raise AssertionError("warming without a plant temperature was accepted")
    new = {1000 + 100 * s + m for s in range(16) for m in range(12)}
    old = [100 + 10 * s + m for s in range(16) for m in range(12)]
    assert len(new) == 16 * 12 and len(set(old)) < len(old) and 100 + 10 * 0 + 10 == 100 + 10 * 1 + 0
    sched = plant_schedule("V3")
    log, net, plant, truths = synthetic_log("Net3", months=3, seed=2, schedule=sched, return_truth=True)
    assert [t["month_seed"] for t in truths] == [1200, 1201, 1202]
    assert list(plant.columns) == ["month", "temp_C", "toc_mgL", "dose_mgL"] and plant.temp_C.tolist() == [10.5, 10.0, 10.5]
    assert (log.temp_C == log.month.map(dict(zip(plant.month, plant.temp_C)))).all() and len(log) == 3 * 13
    assert all(t["chem"]["warming"]["soil_temp_C"] == sched[i]["soil_temp_C"] for i, t in enumerate(truths))
    return {**out, "n_old_formula_collisions_16_seeds_12_months": len(old) - len(set(old))}


@check("seasonal", quick=False)
def seasonal_direction():
    """Raising the temperature lowers the predicted chlorine: a SeasonalSimGP24 fitted on Net3 scenario 0's 8 samples
    (logged at 10.5 C) and predicted at 19.5 C instead of 10.5 C has a lower hourly median at every junction and hour
    (equal only where every member sits at the 0.02 mg/L floor), for the bank (M) and for bulk-only Arrhenius (M1b),
    and a lower daily-minimum median at every junction.  The bank blocks are read from outputs/cache when the experiment
    has cached them and otherwise built in memory (never written), so the result is the same either way."""
    from .seasonal import SeasonalSimGP24, covariate_bank
    sc, X, S, _ = _fitted_net3()
    out = {}
    for tag, mode in (("M", "arrhenius"), ("M1b", "mass_transfer_only")):
        bank = covariate_bank(sc, [10.5, 19.5], wall_mode=mode, cache_dir=CACHE, cache="read")
        mod = SeasonalSimGP24(sc, X, bank, seed=0, cache_dir=CACHE).fit(S.assign(temp_C=10.5), target_temp_C=10.5)
        cold_h, cold_d = mod.predict_hours()[0], mod.predict_daily_min()["median"]
        mod.set_target(19.5)
        warm_h, warm_d = mod.predict_hours()[0], mod.predict_daily_min()["median"]
        assert (warm_h <= cold_h + 1e-12).all(), tag
        strictly = (warm_h.mean(axis=0) < cold_h.mean(axis=0))
        assert strictly.all(), (tag, int((~strictly).sum()))
        lower_d = warm_d < cold_d
        assert lower_d.all(), (tag, int((~lower_d).sum()))
        out[tag] = {"mean_ln_drop_hourly": float((cold_h - warm_h).mean()),
                    "median_daily_min_ratio_warm_over_cold": float((warm_d / cold_d).median()),
                    "n_junctions_lower_daily_min": int(lower_d.sum()), "n_junctions": len(lower_d)}
    return out


@check("seasonal", quick=False)
def pilot_plant_log_path():
    """pilot.validate with a plant log: on the example grab log (Net3, January to June) and docs/example_plant_log.csv
    it runs the temperature-aware model and adds the month's temperature, kb20 and the hypothesis posterior to the
    summary; with a plant log at 20 C every month it reproduces the committed temperature-blind path's predictions
    (to 1e-6 mg/L: the two paths sum the same numbers in a different order, and the GP's hyperparameter search turns
    that last-bit difference into about 1e-8 mg/L; the first full run measured 1.3e-8 against a 1e-9 bar set before any
    measurement).  The bank is read from outputs/cache or built in memory, never written."""
    from .pilot import load_log, load_plant, validate
    log = load_log(os.path.join(REPO, "docs", "example_grab_log.csv"), os.path.join(REPO, "docs", "example_tap_map.csv"))
    plant = load_plant(os.path.join(REPO, "docs", "example_plant_log.csv"))
    preds0, summ0 = validate("Net3", log, cache_dir=CACHE)
    p20 = plant.assign(temp_C=20.0)
    preds20, summ20 = validate("Net3", log, cache_dir=CACHE, plant=p20, bank_cache="read")
    d20 = float(np.abs(preds20.pred.values - preds0.pred.values).max())
    assert d20 < 1e-6, d20
    preds, summ = validate("Net3", log, cache_dir=CACHE, plant=plant, bank_cache="read")
    assert {"temp_C", "kb20", "P_H0", "map_hypothesis"} <= set(summ.columns)
    assert summ.groupby("held_out_month").temp_C.first().tolist() == [12.5, 15.0, 17.5]
    a = summ[summ.taps == "all"]
    return {"max_abs_diff_pred_plant_20C_vs_blind_mgL": d20,
            "rmse_by_month_blind": [round(float(x), 6) for x in summ0[summ0.taps == "all"].rmse],
            "rmse_by_month_plant": [round(float(x), 6) for x in a.rmse], "kb20_by_month": [round(float(x), 6) for x in a.kb20]}


SEASON_NETS = ("Net3", "Net2")
SINGLE_THREAD_ENV = {k: "1" for k in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
# numpy here uses Apple's Accelerate, whose thread count only these variables set before numpy loads can fix
# (threadpoolctl does not see it); the experiment's workers run single-threaded, so the recompute runs in a
# subprocess started the same way, or its sums would differ in the last bit.
SEASON_ROWS_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap.seasonal import covariate_bank, run_task
from residualmap.simulate import nominal_scenario
cache = sys.argv[1]
sc = nominal_scenario("Net3")
bank = covariate_bank(sc, [10.0, 10.5, 12.5], cache_dir=cache, cache="read")
res = run_task("Net3", 0, "V1", bank, None, None, cache, months=4, targets=(4,), models={"B0", "M"}, extrapolate=False)
def conv(v):
    return v if isinstance(v, (str, bool)) else float(v)
print("SEASON_ROWS_JSON " + json.dumps([{k: conv(v) for k, v in r.items()} for r in res["rows"]]))
'''


@check("seasonal", quick=False)
def season_outputs_reproduce():
    """Task 10's committed outputs: summary_season_<net>.json for Net3 and Net2 carries the acceptance key A1 to A7 and
    season_<net>.csv exists; recomputing Net3 seed 0, V1, the April rolling fits of B0 and M, in a subprocess with
    single-threaded BLAS as in the experiment's workers, gives the committed CSV rows exactly."""
    out = {}
    for net in SEASON_NETS:
        d = json.load(open(os.path.join(OUT_DIR, f"summary_season_{net}.json")))
        assert all(f"A{i}" in d["acceptance"] for i in range(1, 8)), net
        assert os.path.exists(os.path.join(OUT_DIR, f"season_{net}.csv")), net
        out[net] = {(f"{k}_pass" if "pass" in v else f"{k}_stop"): v.get("pass", v.get("stop")) for k, v in d["acceptance"].items()}
    r = subprocess.run([PY, "-c", SEASON_ROWS_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV})
    line = [x for x in r.stdout.splitlines() if x.startswith("SEASON_ROWS_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"recompute failed (exit {r.returncode}): {r.stderr[-2000:]}")
    rows = json.loads(line[-1][len("SEASON_ROWS_JSON "):])
    committed = pd.read_csv(os.path.join(OUT_DIR, "season_Net3.csv"), float_precision="round_trip")
    diffs, n_compared = {}, 0
    for row in rows:
        want = committed[(committed.variant == "V1") & (committed.test == "rolling") & (committed.seed == 0)
                         & (committed.month == 4) & (committed.model == row["model"])].iloc[0]
        bad = {}
        for k, v in row.items():
            if k not in want.index or isinstance(v, (str, bool)):
                continue
            n_compared += 1
            if not (pd.isna(v) and pd.isna(want[k])) and v != want[k]:
                bad[k] = (v, want[k])
        diffs[row["model"]] = bad
    assert not any(diffs.values()), diffs
    return {"acceptance": out, "net3_seed0_april_rows_identical": sorted(diffs), "n_values_compared": n_compared}


# ----------------------------------------------------------------------------- task 10b: the wall's response, learned
TASK10_TAG_T12_5_E8000 = "72f9c19f47"     # cache_tag(Chemistry(temp_C=12.5)) in the committed task-10 checks report


@check("seasonal_wall")
def wall_condition_canonical_forms():
    """Chemistry's wall E/R (task 10b) leaves every task-10 condition as it was: the committed T12.5/E8000 tag is
    unchanged; a wall E/R equal to the bulk E/R IS task 10's M condition and a wall E/R of 0 IS its M1b condition
    (equal objects, equal tags, equal cache names), so M2's bank reads those cached files instead of rebuilding them;
    any other wall E/R gets its own tag and label, scales the wall by f(T; wall E/R) and the bulk by f(T; E/R); a
    negative wall E/R, or one given with wall_mode 'mass_transfer_only', is refused.  wall_bank_cache_files for the
    full 3 x 5 grid contains every file of task 10's M and M1b banks."""
    from . import simgp
    from .chemistry import Chemistry, arrhenius, cache_tag
    from .seasonal import (WALL_ER_HYPOTHESES_K, all_bank_temps, bank_cache_files, nested_prior, wall_bank_cache_files,
                           wall_pairs, WallBank)
    from .simulate import nominal_scenario
    assert cache_tag(Chemistry(temp_C=12.5)) == TASK10_TAG_T12_5_E8000
    out = {}
    for e in (5000.0, 8000.0, 12000.0):
        m = Chemistry(temp_C=10.0, er_K=e, wall_er_K=e)
        b = Chemistry(temp_C=10.0, er_K=e, wall_er_K=0.0)
        assert m == Chemistry(temp_C=10.0, er_K=e) and cache_tag(m) == cache_tag(Chemistry(temp_C=10.0, er_K=e))
        assert b == Chemistry(temp_C=10.0, er_K=e, wall_mode="mass_transfer_only")
        assert cache_tag(b) == cache_tag(Chemistry(temp_C=10.0, er_K=e, wall_mode="mass_transfer_only"))
    w = Chemistry(temp_C=10.0, er_K=8000.0, wall_er_K=2500.0)
    assert w.kb_scale() == arrhenius(10.0, 8000.0) and w.kw_scale() == arrhenius(10.0, 2500.0)
    assert w.label() == "free_T10_E8000_W2500"
    tags = {cache_tag(Chemistry(temp_C=10.0, er_K=b_, wall_er_K=w_)) for b_, w_ in wall_pairs()}
    assert len(tags) == 15, len(tags)
    for bad in ({"wall_er_K": -1.0}, {"wall_er_K": 2500.0, "wall_mode": "mass_transfer_only"}):
        try:
            Chemistry(temp_C=10.0, **bad)
        except ValueError:
            continue
        raise AssertionError(f"Chemistry accepted {bad}")
    sc = nominal_scenario("Net3")
    temps, _ = all_bank_temps()
    t10 = set(bank_cache_files(sc, temps, cache_dir=CACHE)) | set(bank_cache_files(sc, temps, wall_mode="mass_transfer_only",
                                                                                     cache_dir=CACHE))
    allf = wall_bank_cache_files(sc, temps, cache_dir=CACHE)
    assert t10 <= set(allf) and len(set(allf)) == len(allf) == 6 * 15, (len(t10), len(allf))
    p = nested_prior(WallBank([], np.zeros(1), {}, wall_pairs(), "wall_er", "full"))
    assert p[0] == 0.25 and np.allclose(p[1:], 0.05, rtol=0, atol=1e-15) and abs(p.sum() - 1) < 1e-12
    out.update(tag_T10_E8000_W2500=cache_tag(w), wall_grid_K=list(WALL_ER_HYPOTHESES_K),
               wall_ratio_20C_over_10C=[round(1.0 / arrhenius(10.0, x), 4) for x in WALL_ER_HYPOTHESES_K],
               n_task10_files_reused=len(t10), n_new_files_per_network=len(allf) - len(t10))
    return out


@check("seasonal_wall", quick=False)
def wall_bank_and_model_nest_m_and_m1b():
    """On the 75-member decay grid, built in memory (cache 'off'): the wall bank's (8000, 8000) block at 10.5 C equals a
    direct build of task 10's M condition, its (8000, 0) block a direct build of the M1b condition and its (8000, 2500)
    block a direct build of that condition, bit for bit; at 20 C every pair is the committed grid; at 10.5 C, for a fixed
    bulk E/R, chlorine never falls as the wall E/R rises.  nested('M') and nested('M1b') are task 10's banks over the
    same arrays.  The model: M2 fitted on Net3 scenario 0's 8 samples (logged at 10.5 C, predicting 19.5 C) with all its
    prior mass on H0 and the pairs (b, b) gives M's posterior and hourly predictions; on H0 and the pairs (b, 0), M1b's."""
    from .chemistry import Chemistry
    from .seasonal import SeasonalSimGP24, WallSeasonalSimGP24, wall_bank
    from .simgp import build_grid_24h, simulator_grid_24h
    sc, X, S, _ = _fitted_net3()
    pairs = ((5000.0, 0.0), (5000.0, 5000.0), (8000.0, 0.0), (8000.0, 2500.0), (8000.0, 8000.0), (12000.0, 0.0),
             (12000.0, 12000.0))
    wb = wall_bank(sc, [10.5, 19.5, 20.0], pairs=pairs, grid="decay", cache_dir=CACHE, cache="off")
    p0, Z0 = simulator_grid_24h(sc, CACHE, "decay")
    assert all(np.array_equal(wb.blocks[(20.0, p)], Z0) for p in pairs)
    for p, cond in (((8000.0, 8000.0), Chemistry(temp_C=10.5, er_K=8000.0)),
                    ((8000.0, 0.0), Chemistry(temp_C=10.5, er_K=8000.0, wall_mode="mass_transfer_only")),
                    ((8000.0, 2500.0), Chemistry(temp_C=10.5, er_K=8000.0, wall_er_K=2500.0))):
        assert np.array_equal(wb.blocks[(10.5, p)], build_grid_24h(sc, "decay", cond=cond)[1]), p
    z0, z25, z8 = (wb.blocks[(10.5, (8000.0, w))] for w in (0.0, 2500.0, 8000.0))
    gaps = {"W2500_vs_W0": float((z25 - z0).min()), "W8000_vs_W2500": float((z8 - z25).min())}
    assert all(v >= -1e-6 for v in gaps.values()), gaps
    out = {"min_ln_gap_at_10.5C_bulk8000": gaps}
    H = wb.hypotheses
    s = S.assign(temp_C=10.5)
    for name in ("M", "M1b"):
        nb = wb.nested(name)
        assert nb.ers == (5000.0, 8000.0, 12000.0) and nb.wall_mode == ("arrhenius" if name == "M" else "mass_transfer_only")
        keep = [f"B{b:g}_W{(b if name == 'M' else 0.0):g}" for b in nb.ers]
        assert all(nb.blocks[(T, b)] is wb.blocks[(T, (b, b if name == "M" else 0.0))] for T in nb.temps for b in nb.ers)
        pm = np.array([1.0 if h == "H0" or h in keep else 0.0 for h in H])
        m2 = WallSeasonalSimGP24(sc, X, wb, prior_mass=pm / pm.sum(), seed=0, cache_dir=CACHE).fit(s, target_temp_C=19.5)
        mm = SeasonalSimGP24(sc, X, nb, seed=0, cache_dir=CACHE).fit(s, target_temp_C=19.5)
        idx = [0] + [H.index(k) for k in keep]
        W2 = m2.W_.reshape(len(H), -1, m2.W_.shape[1])[idx].reshape(-1, m2.W_.shape[1])
        dW = float(np.abs(W2 - mm.W_).max())
        a, b = m2.predict_hours(), mm.predict_hours()
        d_mu, d_sd = float(np.abs(a[0] - b[0]).max()), float(np.abs(a[1] - b[1]).max())
        assert dW < 1e-12 and d_mu < 1e-10 and d_sd < 1e-10, (name, dW, d_mu, d_sd)
        assert abs(m2.kb20() - mm.kb20()) < 1e-12
        out[f"M2_as_{name}"] = {"max_abs_diff_W": dW, "max_abs_diff_hourly_z_mu": d_mu, "max_abs_diff_hourly_z_sd": d_sd}
    m2 = WallSeasonalSimGP24(sc, X, wb, seed=0, cache_dir=CACHE).fit(s, target_temp_C=19.5)
    pc = m2.posterior_columns()
    assert abs(pc["P_H0"] + sum(v for k, v in pc.items() if k.startswith("P_W")) - 1) < 1e-12
    assert abs(sum(v for k, v in pc.items() if k.startswith("P_W")) - sum(v for k, v in pc.items() if k.startswith("P_B"))) < 1e-12
    return out


SEASON2_ROWS_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap.seasonal import run_task, wall_bank
from residualmap.simulate import nominal_scenario
cache = sys.argv[1]
sc = nominal_scenario("Net3")
wb = wall_bank(sc, [12.5, 17.5, 19.5, 20.0], cache_dir=cache, cache="read")
res = run_task("Net3", 16, "V1", wb.nested("M"), wb.nested("M1b"), None, cache, targets=(10,),
               models={"B0", "M", "M1b", "M2"}, temp_errors=(), truth_wall_law="arrhenius", bank_m2=wb)
def conv(v):
    return v if isinstance(v, (str, bool)) else float(v)
print("SEASON2_ROWS_JSON " + json.dumps([{k: conv(v) for k, v in r.items()} for r in res["rows"]]))
'''


@check("seasonal_wall", quick=False)
def season2_outputs_reproduce():
    """Task 10b's committed outputs: summary_season2_<net>.json for Net3 and Net2 carries A2 to A7 and R1 for both worlds,
    and its CSVs exist and are under 1 MB; recomputing Net3 seed 16, world W2, V1 (fit July to September, predict
    October and December; B0, M, M1b, M2 and the oracle) in a single-threaded subprocess gives the committed
    July-to-September-to-December rows (one month, so the aggregated row is that month's) and the October fit's row of
    the windows file, to the CSVs' 6 significant digits."""
    from .seasonal import CSV_FLOAT_10B, SUM_COLS_10B
    out = {}
    for net in SEASON_NETS:
        d = json.load(open(os.path.join(OUT_DIR, f"summary_season2_{net}.json")))
        for w in ("W1", "W2"):
            assert all(k in d["acceptance"][w] for k in ("A2", "A3", "A4", "A5", "A6", "A7", "R1")), (net, w)
        sizes = {f: os.path.getsize(os.path.join(OUT_DIR, f)) for f in (f"season2_{net}.csv", f"season2_windows_{net}.csv")}
        assert all(v < 1_000_000 for v in sizes.values()), sizes
        out[net] = {"A7_stop": {w: d["acceptance"][w]["A7"]["stop"] for w in ("W1", "W2")}, "csv_bytes": sizes}
    r = subprocess.run([PY, "-c", SEASON2_ROWS_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV})
    line = [x for x in r.stdout.splitlines() if x.startswith("SEASON2_ROWS_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"recompute failed (exit {r.returncode}): {r.stderr[-2000:]}")
    rows = pd.DataFrame(json.loads(line[-1][len("SEASON2_ROWS_JSON "):]))
    agg = pd.read_csv(os.path.join(OUT_DIR, "season2_Net3.csv"), float_precision="round_trip")
    win = pd.read_csv(os.path.join(OUT_DIR, "season2_windows_Net3.csv"), float_precision="round_trip")

    def same(a, b):
        if pd.isna(a) and pd.isna(b):
            return True
        return float(CSV_FLOAT_10B % float(a)) == float(b)
    bad, n = {}, 0
    for model in ("B0", "M", "M1b", "M2", "oracle"):
        got = rows[(rows.test == "jul_sep_to_dec") & (rows.model == model)]
        want = agg[(agg.world == "W2") & (agg.variant == "V1") & (agg.test == "jul_sep_to_dec") & (agg.seed == 16)
                   & (agg.model == model) & (agg.temp_error == "none")]
        assert len(got) == 1 and len(want) == 1 and int(want.n_months.iloc[0]) == 1, (model, len(got), len(want))
        for c in SUM_COLS_10B:
            n += 1
            if not same(got[c].iloc[0], want[c].iloc[0]):
                bad[f"{model}.{c}"] = (got[c].iloc[0], want[c].iloc[0])
    wrow = win[(win.world == "W2") & (win.variant == "V1") & (win.seed == 16) & (win.month == 10)].iloc[0]
    for model in ("B0", "M", "M1b", "M2"):
        g = rows[(rows.test == "rolling") & (rows.model == model)].iloc[0]
        for c in [c for c in win.columns if c.startswith(f"{model}_")]:
            n += 1
            if not same(g[c[len(model) + 1:]], wrow[c]):
                bad[c] = (g[c[len(model) + 1:]], wrow[c])
    assert not bad, bad
    return {**out, "net3_W2_seed16_V1_values_compared": n, "identical_to_6_significant_digits": not bad}


# ----------------------------------------------------------------------------- task 11: organics (TOC)
@check("organics")
def clark_truth_condition_and_draws():
    """Task 11's truth side.  Chemistry: kinetics 'clark' needs toc_mgL and phi, phi goes with 'clark' only, and an unset
    phi leaves every earlier cache tag as it was (the committed T12.5/E8000 tag).  On Net3 and Net2 (seed 0, month seed
    1000, 20 C, TOC 2.0), a Clark truth consumes the committed draws exactly as the first-order truth does (the same
    bulk month factor and source doses, and the same water age on the truth's network, so no draw moved) and records
    k2 = kb u / (phi TOC_ref) and CL = mean dose - phi TOC, each rounded to the 4 decimals the .inp keeps.  The dose
    response: with every draw the same and the plant dose 1.5 instead of 1.2 mg/L, first-order chlorine scales by
    exactly 1.25 (to 5e-5 mg/L: EPANET's quality tolerance in the legacy units is 1e-5 mg/L; the first run measured
    6.4e-6 mg/L against a 1e-6 bar set before any measurement) and Clark chlorine does not (its mean daily minimum
    moves by a ratio more than 0.01 away from 1.25).  Clark with in-network warming is refused."""
    from .chemistry import TOC_REF_MGL, Chemistry, Warming, cache_tag
    from .simulate import build_scenario
    for bad in ({"kinetics": "clark"}, {"kinetics": "clark", "toc_mgL": 2.0}, {"phi": 0.5},
                {"kinetics": "clark", "toc_mgL": 2.0, "phi": 0.0}):
        try:
            Chemistry(**bad)
        except ValueError:
            continue
        raise AssertionError(f"Chemistry accepted {bad}")
    assert cache_tag(Chemistry(temp_C=12.5)) == TASK10_TAG_T12_5_E8000
    out = {"clark_label": Chemistry(kinetics="clark", toc_mgL=3.0, phi=0.85, temp_C=20.0).label()}
    for net in ("Net3", "Net2"):
        kw = TRUTHS[net]
        kb = kw.get("kb_per_day", 0.40)
        f = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=20.0, toc_mgL=2.0), truth_age=True, **kw)
        c = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=20.0, toc_mgL=2.0, kinetics="clark", phi=0.85),
                           truth_age=True, **kw)
        assert c.chem["bulk_month_factor"] == f.chem["bulk_month_factor"] and c.chem["source_doses_mgL"] == f.chem["source_doses_mgL"]
        assert np.array_equal(c.truth_age_by_hour_h.values, f.truth_age_by_hour_h.values), net
        cl = c.chem["clark"]
        u, d = c.chem["bulk_month_factor"], float(np.mean(list(c.chem["source_doses_mgL"].values())))
        assert cl["k2_L_per_mg_day"] == round(kb * u / (0.85 * TOC_REF_MGL), 4) and cl["limiting_mgL"] == round(d - 0.85 * 2.0, 4)
        ratios = {}
        for kin, phi in (("first", None), ("clark", 0.85)):
            lo = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=20.0, toc_mgL=3.0, kinetics=kin, phi=phi),
                                **kw).truth_by_hour
            hi = build_scenario(net, 0, month_seed=1000, chem=Chemistry(temp_C=20.0, toc_mgL=3.0, kinetics=kin, phi=phi),
                                **{**kw, "source_dose": 1.5}).truth_by_hour
            ratios[kin] = (float((hi - 1.25 * lo).abs().values.max()), float(hi.min().mean() / lo.min().mean()))
        assert ratios["first"][0] < 5e-5, ratios
        assert abs(ratios["clark"][1] - 1.25) > 0.01, ratios
        out[net] = {"k2": cl["k2_L_per_mg_day"], "CL_mgL": cl["limiting_mgL"],
                    "first_order_max_abs_dev_from_x1.25_mgL": ratios["first"][0],
                    "clark_mean_daily_min_ratio_1.5_over_1.2": ratios["clark"][1]}
    try:
        build_scenario("Net3", 0, chem=Chemistry(temp_C=15.0, toc_mgL=2.0, kinetics="clark", phi=0.85),
                       warming=Warming(soil_temp_C=18.0))
    except NotImplementedError:
        pass
    else:
        raise AssertionError("Clark kinetics with in-network warming was accepted")
    return out


@check("organics")
def toc_bank_and_model_nest_simgp24():
    """The TOC bank (decay grid, built in memory): its block at TOC_ref IS the committed grid, its TOC 3.0 block equals a
    direct build of Chemistry(toc_mgL=3.0) bit for bit, and chlorine never rises with TOC (3.0 <= 1.5 everywhere).
    M_TOC (full grid; no TOC block is needed) fitted on Net3 scenario 0's 8 samples logged at TOC 2.0: its posterior is
    SimGP24's split as the prior (H0 1/4, H_TOC 3/4) and stays there, its MAP member and hourly predictions are SimGP24's
    (to 1e-10), under every prior; the logged dose is an exact ln-offset.  Direction (the TOC 3.0 block read from
    outputs/cache when the experiment cached it, else built in memory and never written): told the month's TOC is 3.0
    instead of 2.0, M_TOC predicts less chlorine at every junction and hour."""
    from .organics import TOC_PRIOR, TocBank, TocSimGP24, toc_bank, toc_condition
    from .simgp import build_grid_24h, simulator_grid_24h
    from .simulate import nominal_scenario
    sc0 = nominal_scenario("Net3")
    tb = toc_bank(sc0, [1.5, 2.0, 3.0], grid="decay", cache_dir=CACHE, cache="off")
    p0, Z0 = simulator_grid_24h(sc0, CACHE, "decay")
    assert isinstance(tb, TocBank) and tb.hypotheses == ("H0", "TOC") and tb.tocs == [1.5, 2.0, 3.0]
    assert tb.blocks[(2.0, "TOC")] is tb.h0 and np.array_equal(tb.h0, Z0)
    assert np.array_equal(tb.blocks[(3.0, "TOC")], build_grid_24h(sc0, "decay", cond=toc_condition(3.0))[1])
    gap = float((tb.blocks[(1.5, "TOC")] - tb.blocks[(3.0, "TOC")]).min())
    assert gap >= -1e-6, gap
    sc, X, S, m = _fitted_net3()
    full = toc_bank(sc, [2.0], cache_dir=CACHE, cache="read")
    s2 = S.assign(toc_mgL=2.0)
    out = {"min_ln_gap_toc1.5_minus_toc3.0_decay_grid": gap}
    for prior in (None, (0.5, 0.5), (0.0, 1.0)):
        mt = TocSimGP24(sc, X, full, prior_mass=prior, seed=0, cache_dir=CACHE).fit(s2, target_toc_mgL=2.0)
        p = np.asarray(TOC_PRIOR if prior is None else prior)
        dW = max(float(np.abs(mt.W_.reshape(2, len(m.params), -1)[b] - p[b] * m.W_).max()) for b in range(2))
        a, b = m.predict_hours(), mt.predict_hours()
        d_mu, d_sd = float(np.abs(a[0] - b[0]).max()), float(np.abs(a[1] - b[1]).max())
        assert dW < 1e-15 and d_mu < 1e-10 and d_sd < 1e-10, (prior, dW, d_mu, d_sd)
        assert abs(mt.p_h0() - p[0]) < 1e-12 and mt.map_params_ == m.map_params_
        out[f"prior_{p[0]:g}"] = {"max_abs_diff_W": dW, "max_abs_diff_hourly_z_mu": d_mu, "max_abs_diff_hourly_z_sd": d_sd}
    Sa = s2[s2.y >= 0.02]
    m1 = TocSimGP24(sc, X, full, seed=0, cache_dir=CACHE).fit(Sa, target_toc_mgL=2.0)
    m2 = TocSimGP24(sc, X, full, seed=0, cache_dir=CACHE).fit(Sa.assign(y=Sa.y * 2, dose_ratio=2.0), target_toc_mgL=2.0,
                                                               target_dose_ratio=2.0)
    d_dose = float(np.abs(m2.predict_hours()[0] - m1.predict_hours()[0] - np.log(2.0)).max())
    assert d_dose < 1e-8, d_dose
    tb3 = toc_bank(sc, [2.0, 3.0], cache_dir=CACHE, cache="read")
    mt = TocSimGP24(sc, X, tb3, seed=0, cache_dir=CACHE).fit(s2, target_toc_mgL=2.0)
    h2 = mt.predict_hours()[0]
    h3 = mt.set_target(3.0).predict_hours()[0]
    assert (h3 <= h2 + 1e-12).all() and (h3.mean(axis=0) < h2.mean(axis=0)).all()
    return {**out, "dose_offset_max_abs_err": d_dose, "mean_ln_drop_toc3_vs_toc2": float((h2 - h3).mean())}


@check("organics", quick=False)
def m2toc_bank_nests_m2():
    """The secondary test's bank and model (decay grid, built in memory, never cached): at TOC_ref a (temperature, TOC)
    condition is M2's block (the same array), at 20 C it is the TOC bank's block, and a (10.5 C, TOC 2.5) block equals a
    direct build of Chemistry(temp_C=10.5, er_K=8000, toc_mgL=2.5) bit for bit and has at most the TOC 2.0 block's
    chlorine.  M2_TOC fitted on samples whose months all sit at TOC_ref is M2: the same posterior and hourly predictions."""
    from .chemistry import Chemistry
    from .organics import M2TocSimGP24, cond_key, m2toc_bank, toc_bank
    from .seasonal import WallSeasonalSimGP24
    from .simgp import build_grid_24h
    sc, X, S, _ = _fitted_net3()
    tb = toc_bank(sc, [1.5, 2.0], grid="decay", cache_dir=CACHE, cache="off")
    kb, m2 = m2toc_bank(sc, [(10.5, 2.0), (10.5, 2.5), (20.0, 1.5)], tb, grid="decay", cache_dir=CACHE)
    k20, k25, k_hot = cond_key(10.5, 2.0), cond_key(10.5, 2.5), cond_key(20.0, 1.5)
    assert all(kb.blocks[(k20, p)] is m2.blocks[(10.5, p)] for p in kb.ers)
    assert all(kb.blocks[(k_hot, p)] is tb.blocks[(1.5, "TOC")] for p in kb.ers)
    direct = build_grid_24h(sc, "decay", cond=Chemistry(temp_C=10.5, er_K=8000.0, toc_mgL=2.5))[1]
    assert np.array_equal(kb.blocks[(k25, (8000.0, 8000.0))], direct)
    gap = min(float((kb.blocks[(k20, p)] - kb.blocks[(k25, p)]).min()) for p in kb.ers)
    assert gap >= -1e-6, gap
    s = S.assign(temp_C=10.5, cond_key=k20)
    a = WallSeasonalSimGP24(sc, X, m2, seed=0, cache_dir=CACHE).fit(s[["junction", "hour", "y", "temp_C"]], target_temp_C=10.5)
    b = M2TocSimGP24(sc, X, kb, seed=0, cache_dir=CACHE).fit(s[["junction", "hour", "y", "cond_key"]], target_temp_C=k20)
    dW = float(np.abs(a.W_ - b.W_).max())
    ha, hb = a.predict_hours(), b.predict_hours()
    d_mu, d_sd = float(np.abs(ha[0] - hb[0]).max()), float(np.abs(ha[1] - hb[1]).max())
    assert dW < 1e-15 and d_mu < 1e-10 and d_sd < 1e-10, (dW, d_mu, d_sd)
    return {"min_ln_gap_toc2.0_minus_toc2.5_at_10.5C": gap, "max_abs_diff_W": dW, "max_abs_diff_hourly_z_mu": d_mu}


@check("organics")
def toc_schedules_and_log():
    """The organics schedules and log: S follows chemistry.monthly_toc, C is TOC_ref every month, D steps the plant dose
    from 1.2 to 1.5 mg/L in July, the temperature is 20 C throughout (task 10's V1 schedule for the secondary test);
    O2 and O3 carry kinetics 'clark' and their phi.  A 3-month O2 log under D: every reading carries its month's TOC and
    dose, the plant log logs both, each month's truth is Clark with that TOC and dose, and task 10's schedules still give
    readings with no TOC or dose column."""
    from .chemistry import monthly_temperature, monthly_toc
    from .organics import DOSE_STEP, TRUTHS, toc_schedule
    from .pilot import synthetic_log
    from .seasonal import plant_schedule
    s, c, d = toc_schedule("O1", "S"), toc_schedule("O2", "C"), toc_schedule("O2", "D")
    assert [x["toc_mgL"] for x in s] == [monthly_toc(m) for m in range(1, 13)] and all(x["temp_C"] == 20.0 for x in s)
    assert all(x["toc_mgL"] == 2.0 for x in c) and "kinetics" not in s[0] and c[0]["kinetics"] == "clark"
    assert [x["dose_mgL"] for x in d] == [1.2] * (DOSE_STEP[0] - 1) + [1.5] * (13 - DOSE_STEP[0])
    assert toc_schedule("O3", "S")[0]["phi"] == TRUTHS["O3"]["phi"] == 0.5
    assert [x["temp_C"] for x in toc_schedule("O2", "S", seasonal=True)] == [monthly_temperature(m) for m in range(1, 13)]
    sched = d[5:8]          # June (1.2 mg/L), July and August (1.5 mg/L)
    log, net, plant, truths = synthetic_log("Net3", months=3, seed=2, schedule=sched, return_truth=True)
    assert plant.toc_mgL.tolist() == [1.5, 1.5, 1.5] and plant.dose_mgL.tolist() == [1.2, 1.5, 1.5]
    assert (log.toc_mgL == 1.5).all() and log.groupby("month").dose_mgL.first().tolist() == [1.2, 1.5, 1.5]
    assert all(t["chem"]["kinetics"] == "clark" and t["chem"]["clark"]["phi_mg_per_mgC"] == 0.85 for t in truths)
    doses = [np.mean(list(t["chem"]["source_doses_mgL"].values())) for t in truths]
    assert 1.08 <= doses[0] <= 1.32 and 1.35 <= doses[1] <= 1.65
    log10 = synthetic_log("Net3", months=1, seed=2, schedule=plant_schedule("V1"))[0]
    assert "toc_mgL" not in log10 and "dose_mgL" not in log10
    return {"mean_source_dose_by_month": [round(float(x), 6) for x in doses],
            "CL_by_month": [t["chem"]["clark"]["limiting_mgL"] for t in truths]}


ORGANICS_ROWS_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap.organics import ALL_TOCS, run_task, toc_bank
from residualmap.simulate import nominal_scenario
cache = sys.argv[1]
sc = nominal_scenario("Net3")
tb = toc_bank(sc, ALL_TOCS, cache_dir=cache, cache="read")
res = run_task("Net3", 32, "O2", "S", {"TOC": tb}, cache, targets=(11,), models={"B0", "M_TOC", "oracle"},
               extrapolate=False)
def conv(v):
    return v if isinstance(v, (str, bool)) else float(v)
print("ORGANICS_ROWS_JSON " + json.dumps([{k: conv(v) for k, v in r.items()} for r in res["rows"]]))
'''


@check("organics", quick=False)
def organics_outputs_reproduce():
    """Task 11's committed outputs: summary_organics_<net>.json for Net3 and Net2 carries the acceptance key (T1, T2, N,
    T3, D, G, adoption) and its CSVs exist and are under 1 MB; recomputing Net3 seed 32, truth O2, schedule S, the
    first-storm month (November, fitted on August to October; B0, M_TOC and the oracle) in a single-threaded subprocess
    gives the committed rows (one month, so the aggregated row is that month's) and the November row of the windows file,
    to the CSVs' 6 significant digits."""
    from .organics import CSV_FLOAT, SUM_COLS
    out = {}
    for net in SEASON_NETS:
        d = json.load(open(os.path.join(OUT_DIR, f"summary_organics_{net}.json")))
        assert all(k in d["acceptance"] for k in ("T1", "T2", "N", "T3", "D", "G", "adoption")), net
        sizes = {f: os.path.getsize(os.path.join(OUT_DIR, f)) for f in (f"organics_{net}.csv", f"organics_windows_{net}.csv")}
        assert all(v < 1_000_000 for v in sizes.values()), sizes
        out[net] = {"T3_stop": d["acceptance"]["T3"]["stop"], "csv_bytes": sizes}
    r = subprocess.run([PY, "-c", ORGANICS_ROWS_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV})
    line = [x for x in r.stdout.splitlines() if x.startswith("ORGANICS_ROWS_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"recompute failed (exit {r.returncode}): {r.stderr[-2000:]}")
    rows = pd.DataFrame(json.loads(line[-1][len("ORGANICS_ROWS_JSON "):]))
    agg = pd.read_csv(os.path.join(OUT_DIR, "organics_Net3.csv"), float_precision="round_trip")
    win = pd.read_csv(os.path.join(OUT_DIR, "organics_windows_Net3.csv"), float_precision="round_trip")

    def same(a, b):
        if pd.isna(a) and pd.isna(b):
            return True
        return float(CSV_FLOAT % float(a)) == float(b)
    bad, n = {}, 0
    for model in ("B0", "M_TOC", "oracle"):
        got = rows[rows.model == model]
        want = agg[(agg.truth == "O2") & (agg.variant == "S") & (agg.test == "rolling") & (agg.seed == 32)
                   & (agg.cls == "first_storm") & (agg.model == model)]
        assert len(got) == 1 and len(want) == 1 and int(want.n_months.iloc[0]) == 1, (model, len(got), len(want))
        for c in SUM_COLS:
            n += 1
            if not same(got[c].iloc[0], want[c].iloc[0]):
                bad[f"{model}.{c}"] = (got[c].iloc[0], want[c].iloc[0])
    wrow = win[(win.truth == "O2") & (win.variant == "S") & (win.seed == 32) & (win.month == 11)].iloc[0]
    for model in ("B0", "M_TOC"):
        g = rows[rows.model == model].iloc[0]
        for c in [c for c in win.columns if c.startswith(f"{model}_")]:
            n += 1
            v, w = g[c[len(model) + 1:]], wrow[c]
            if isinstance(v, str) or isinstance(w, str):
                ok = str(v) == str(w)
            else:
                ok = same(v, w)
            if not ok:
                bad[c] = (v, w)
    assert not bad, bad
    return {**out, "net3_O2_S_seed32_november_values_compared": n, "identical_to_6_significant_digits": not bad}


# ----------------------------------------------------------------------------- task 12: the chloramine mode
CA_NETS = ("Net3", "Net2")


def _ca_chain(n_seg=7, diam=0.3, flow=0.01, days=10):
    """Reservoir, then n_seg pipes of one day's travel each, then one demand (plug flow; junction k is k days old).
    Built in code, so its node order differs from EPANET's (msx.run renames the columns; this chain checks that too)."""
    import wntr
    wn = wntr.network.WaterNetworkModel()
    wn.add_reservoir("R", base_head=100.0)
    length = flow / (np.pi * diam ** 2 / 4) * 86400.0
    prev = "R"
    for k in range(1, n_seg + 1):
        wn.add_junction(f"J{k}", base_demand=(flow if k == n_seg else 0.0), elevation=0.0)
        wn.add_pipe(f"P{k}", prev, f"J{k}", length=length, diameter=diam, roughness=130)
        prev = f"J{k}"
    wn.options.time.duration = int(days * 86400)
    wn.options.time.hydraulic_timestep = 3600
    wn.options.time.report_timestep = 3600
    wn.options.time.quality_timestep = 300
    return wn


@check("chloramine")
def batch_port_reproduces_wahman_example():
    """The Python batch port of EPA's model reproduces Wahman 2018's worked example (4 mg/L held 10 days falls to 0.84
    mg/L at pH 7 and 3.2 at pH 9) within 0.15 mg/L, with ASSUMED inputs (Cl2:N 5, 25 C, alkalinity 50, no TOC,
    simultaneous addition); the committed calibration file holds the same numbers.  The prior table's rate falls as the
    pH rises and rises with the Cl2:N ratio."""
    from .chloramine import kb_prior_table, wahman_example
    w = wahman_example()
    assert w["max_abs_diff_mgL"] <= 0.15, w
    cal = json.load(open(os.path.join(REPO, "outputs", "chloramine", "calibration_chloramine.json")))["wahman_example"]
    assert cal["port_mgL"] == w["port_mgL"], (cal["port_mgL"], w["port_mgL"])
    tab = kb_prior_table()
    P = tab.pivot(index="pH", columns="cl2n", values="k_app_per_day")
    assert (np.diff(P.values, axis=0) < 0).all() and (np.diff(P.values, axis=1) > 0).all()
    committed = pd.read_csv(os.path.join(REPO, "outputs", "chloramine", "kb_prior_table.csv"))
    assert np.allclose(committed.k_app_per_day.values, tab.k_app_per_day.values, rtol=1e-5, atol=0)
    return {"port_mgL": w["port_mgL"], "published_mgL": w["published_mgL"], "max_abs_diff_mgL": w["max_abs_diff_mgL"],
            "k_app_range_per_day": [float(P.values.min()), float(P.values.max())]}


@check("chloramine", quick=False)
def msx_port_matches_batch_port():
    """EPA's model in EPANET-MSX (the network truth's reaction model, without the wall) on a plug-flow chain of 1-day
    pipes against the Python batch port at 1, 3 and 7 days, pH 7, 8 and 9 (dose 2.0 mg/L, Cl2:N 4.5, TOC 2,
    alkalinity 100, 20 C): within 0.05 mg/L total chlorine.  No MSX scratch file is left in the working directory."""
    from . import chloramine as C
    from . import msx as M
    out, worst = {}, 0.0
    for ph in (7.0, 8.0, 9.0):
        wn = _ca_chain()
        m = C.build_msx_model(wn, C.rate_constants(ph, 100.0, 20.0))
        M.add_sources(m, wn, {"R": C.source_species(2.0, 4.5, 2.0)})
        res, info = M.run(wn, m)
        last = C.total_chlorine_mgL({sp: res.node[sp] for sp in ("TOTCL", "NH2CL", "NHCL2", "NCL3")}).iloc[-1]
        port = C.batch(ph, 100.0, 20.0, 2.0, 4.5, 2.0, days=8.0, t_eval_days=[1, 3, 7]).total_mgL.values
        msxv = np.array([float(last[f"J{d}"]) for d in (1, 3, 7)])
        d = float(np.abs(msxv - port).max())
        worst = max(worst, d)
        out[str(ph)] = {"msx_mgL": msxv.round(4).tolist(), "port_mgL": port.round(4).tolist(), "max_abs_diff": d,
                        "compiler": info["compiler"], "source_mgL": float(last["R"])}
    assert worst <= 0.05, out
    assert all(abs(v["source_mgL"] - 2.0) < 1e-3 for v in out.values()), "the source column is not the reservoir"
    assert _root_scratch_files() == [] and not [f for f in os.listdir(".") if re.fullmatch(r"(msx|en)[A-Za-z0-9]{6}", f)]
    return {"by_pH": out, "worst_mgL": worst}


_FALLBACK_CHILD = r'''
import json, sys
from residualmap import checks as K, chloramine as C, msx as M
mode = sys.argv[1]
wn = K._ca_chain(n_seg=2, days=3)
m = C.build_msx_model(wn, C.rate_constants(8.0, 100.0, 20.0))
M.add_sources(m, wn, {"R": C.source_species(2.0, 4.5, 2.0)})
out = {"compiler_available": M.compiler_available()}
if mode == "no_gcc":
    try:
        M.run(wn, m, fallback=False)
        out["refused_without_fallback"] = False
    except RuntimeError:
        out["refused_without_fallback"] = True
else:
    res, info = M.run(wn, m)
    last = C.total_chlorine_mgL({sp: res.node[sp] for sp in ("TOTCL", "NH2CL", "NHCL2", "NCL3")}).iloc[-1]
    port = C.batch(8.0, 100.0, 20.0, 2.0, 4.5, 2.0, days=3.0, t_eval_days=[1, 2]).total_mgL.values
    out.update(compiler=info["compiler"], fallback_reason=info["fallback_reason"],
               msx_mgL=[float(last["J1"]), float(last["J2"])], port_mgL=[float(v) for v in port],
               source_mgL=float(last["R"]))
print("RESULT " + json.dumps(out))
'''


@check("chloramine", quick=False)
def msx_falls_back_to_uncompiled():
    """COMPILER GC with no working compiler falls back to COMPILER NONE and gives the same chemistry.  In a subprocess
    whose PATH holds only a `gcc` that fails, MSX's compile fails (MSX error 522), the runner closes the half-open MSX
    project and reruns uncompiled; the result matches the batch port on a 2-pipe plug-flow chain (1 and 2 days old,
    pH 8) within 0.05 mg/L.  With no gcc on the PATH at all, the runner does not try GC, and refuses when the fallback
    is off.  (Found in the task-12 review: the retry used to fail with MSX error 520, a project already open.)"""
    out = {}
    with tempfile.TemporaryDirectory(prefix="rm_nocc_") as tmp:
        fake = os.path.join(tmp, "failing_gcc")
        os.makedirs(fake)
        with open(os.path.join(fake, "gcc"), "w") as fh:
            fh.write("#!/bin/sh\nexit 1\n")
        os.chmod(os.path.join(fake, "gcc"), 0o755)
        for mode, path in (("failing_gcc", fake), ("no_gcc", os.path.join(tmp, "empty"))):
            env = {**os.environ, "PATH": path, "PYTHONPATH": REPO}
            r = subprocess.run([PY, "-c", _FALLBACK_CHILD, mode], cwd=tmp, env=env, capture_output=True, text=True,
                               timeout=900)
            line = [x for x in r.stdout.splitlines() if x.startswith("RESULT ")]
            assert r.returncode == 0 and line, (mode, r.returncode, r.stderr[-2000:])
            out[mode] = json.loads(line[-1][len("RESULT "):])
            out[mode]["loud_message_on_stderr"] = "COMPILER NONE" in r.stderr
    f, n = out["failing_gcc"], out["no_gcc"]
    assert f["compiler_available"] and f["compiler"] == "NONE" and "522" in f["fallback_reason"], f
    diff = float(np.abs(np.array(f["msx_mgL"]) - np.array(f["port_mgL"])).max())
    assert diff <= 0.05 and abs(f["source_mgL"] - 2.0) < 1e-3 and f["loud_message_on_stderr"], f
    assert not n["compiler_available"] and n["refused_without_fallback"], n
    assert _root_scratch_files() == []
    f["fallback_reason"] = f["fallback_reason"].split(" [")[0]          # drop the run-specific pointer text
    return {"failing_gcc": {**f, "max_abs_diff_vs_port_mgL": diff}, "no_gcc": n}


@check("chloramine", quick=False)
def msx_first_order_matches_epanet():
    """One species with EPANET's own first-order bulk (0.40 /day) and mass-transfer-limited wall (0.70 m/day x the
    roughness factor) written as MSX expressions, against EPANET's CHEMICAL quality on Net3 (dose 1.2 mg/L, every
    junction and hour of the last day).  The plan's bar is 0.03 mg/L.  At the repo's 300 s quality step the two engines
    differ by up to about 0.031 mg/L at one junction and hour where a front passes (each engine alone moves by up to
    0.18 mg/L between a 300 s and a 60 s step there), so the bar is checked at a 60 s step in both engines and the 300 s
    difference is recorded."""
    from . import chloramine as C
    from . import msx as M
    from . import simulate as S
    kb, kw, dose = 0.40, 0.70, 1.2
    out, runs = {}, {}
    for step in (300, 60):
        old = S.QUALITY_STEP_S
        S.QUALITY_STEP_S = step
        try:
            ref = S.simulate_nominal_chlorine("Net3", kb, kw, 1.0, dose, file_prefix=os.path.join(os.getcwd(), "fo_ref"))
        finally:
            S.QUALITY_STEP_S = old
            for f in glob.glob("fo_ref.*"):
                os.remove(f)
        wn = S.load("Net3")
        wn.options.time.quality_timestep = step
        m = C.first_order_msx_model(wn, kb, {pn: kw * S.roughness_factor(p.roughness, 1.0) for pn, p in wn.pipes()})
        m.options.timestep = step
        M.add_sources(m, wn, {s: {"CL": dose} for s in S.source_nodes(wn)})
        res, _ = M.run(wn, m)
        q = S._last_day(res.node["CL"], wn.junction_name_list)
        d = (q - ref).abs()
        out[f"step_{step}s"] = {"max_abs_diff_mgL": float(d.values.max()), "mean_abs_diff_mgL": float(d.values.mean()),
                                "p999_abs_diff_mgL": float(np.quantile(d.values, 0.999))}
        runs[step] = (ref, q)
    for i, eng in enumerate(("epanet", "msx")):    # each engine against itself: the quality step's own effect at fronts
        out[f"{eng}_300s_vs_60s"] = {"max_abs_diff_mgL": float((runs[300][i] - runs[60][i]).abs().values.max())}
    assert out["step_60s"]["max_abs_diff_mgL"] <= 0.03, out
    assert out["step_300s"]["mean_abs_diff_mgL"] <= 0.005, out
    return out


@check("chloramine", quick=False)
def chloramine_truth_draws_twin_and_cleanup():
    """The chloramine truth (Net2, seed 900, never scored) consumes the committed draws in their committed order: its
    source doses are the free-chlorine chemistry truth's at the same seed and dose, and the first-order twin has the
    same doses; its chemistry is the seed's own draws on their ranges; the low-pH stress maps the same pH draw onto
    7.0 to 7.5; MSX ran compiled, wn.msx is detached after the run (a later EPANET AGE run of the same network works),
    and no MSX scratch file is left in the working directory."""
    from .chemistry import Chemistry
    from .chloramine import PH_RANGE_STRESS, ChloramineTruth, truth_chemistry, truth_scenario
    from .simulate import build_scenario, hidden_chem_draws
    msx_sc = truth_scenario("Net2", 900, 0.05, "epa_msx")
    twin = truth_scenario("Net2", 900, 0.05, "first")
    free = build_scenario("Net2", 900, source_dose=2.0, kb_per_day=0.10, kw_m_per_day=0.20, chem=Chemistry(temp_C=20.0))
    assert msx_sc.chem["source_doses_mgL"] == twin.chem["source_doses_mgL"] == free.chem["source_doses_mgL"]
    hd = hidden_chem_draws("chloramine", 900)
    want = truth_chemistry(hd)
    assert all(msx_sc.chem[k] == want[k] for k in ("pH", "cl2n", "toc_mgL", "alk_mgL_caco3"))
    lo = truth_chemistry(hd, ChloramineTruth(ph_range=PH_RANGE_STRESS))["pH"]
    assert abs((lo - 7.0) / 0.5 - (want["pH"] - 7.5) / 1.0) < 1e-12 and 7.0 <= lo < 7.5
    assert msx_sc.chem["msx_compiler"] == "GC" and msx_sc.chem["msx_fallback_reason"] is None
    assert msx_sc.wn.msx is None
    dmin, tmin = msx_sc.truth_daily_min, twin.truth_daily_min
    assert (dmin >= 0).all() and float(dmin.max()) <= 2.0 * 1.1 + 1e-6
    assert _root_scratch_files() == [] and not [f for f in os.listdir(".") if re.fullmatch(r"(msx|en)[A-Za-z0-9]{6}", f)]
    return {"pH": msx_sc.chem["pH"], "cl2n": msx_sc.chem["cl2n"], "median_daily_min_msx": float(dmin.median()),
            "median_daily_min_twin": float(tmin.median()), "kb_twin_per_day": twin.chem["kb_twin_per_day"]}


@check("chloramine")
def chloramine_model_prior_and_nesting():
    """Model (c)'s prior: a log prior of zeros gives exactly the uniform posterior (and today's arithmetic when it is
    None); the pH and Cl2:N prior is centred on k_hat and, with no samples, the prior is the posterior.  On Net3's
    chloramine grid (cached) at seed 0's chemistry.  SimGP24's likelihood scale, dose axis and threshold, left unset,
    follow the condition's disinfectant (chloramine's 0.10, DOSES_CA and 0.5 with a chloramine condition)."""
    from .chloramine import ca_condition, k_hat, nominal_ca, prior_log_vector, prior_table_cached
    from .features import build_features
    from .simgp import DOSES_CA, LIK_SD_CA, SimGP24
    sc = nominal_ca("Net3")
    X = build_features(sc)
    def model(lp=None):
        m = SimGP24(sc, X, seed=0, cache_dir=CACHE, grid="chloramine", cond=ca_condition(), lik_sd=LIK_SD_CA,
                    doses=DOSES_CA, threshold=0.5)
        m.log_prior = lp
        return m
    S = pd.DataFrame({"junction": ["15", "123", "247", "60"], "hour": [9, 11, 14, 16], "y": [1.9, 1.1, 0.6, 1.4]})
    m0 = model().fit(S)
    mz = model(np.zeros(len(m0.params))).fit(S)
    assert np.array_equal(m0.W_, mz.W_), float(np.abs(m0.W_ - mz.W_).max())
    lp = prior_log_vector(m0.params, 8.0, 4.5)
    kh = k_hat(8.0, 4.5, prior_table_cached())
    kbs = np.array([p[0] for p in m0.params])
    assert np.isclose(np.exp(np.sum(np.log(kbs) * np.exp(lp)) / np.exp(lp).sum()), kh, rtol=0.6)
    mp = model(lp).fit_prior()
    marg = pd.Series(mp.w_).groupby(kbs).sum()
    assert int(np.argmax(marg.values)) == int(np.argmin(np.abs(np.log(marg.index.values) - np.log(kh))))
    mc = model(lp).fit(S)
    # left unset, the settings follow the condition's disinfectant (review fix): chloramine's with a chloramine
    # condition, the committed free-chlorine ones without one
    from .simgp import DOSE_GRID, LIK_SD
    dflt = SimGP24(sc, X, seed=0, cache_dir=CACHE, grid="chloramine", cond=ca_condition())
    assert (dflt.lik_sd, dflt.doses, dflt.threshold) == (LIK_SD_CA, DOSES_CA, 0.5)
    assert np.array_equal(dflt.fit(S).W_, m0.W_)
    import inspect
    free = {k: v.default for k, v in inspect.signature(SimGP24.__init__).parameters.items()}
    assert all(free[k] is None for k in ("lik_sd", "doses", "threshold")) and (LIK_SD, DOSE_GRID) == (0.35, [0.90, 0.95, 1.00, 1.05, 1.10])
    return {"defaults_follow_disinfectant": True, "k_hat_pH8_cl2n4.5": kh, "prior_mode_kb": float(marg.idxmax()),
            "posterior_geo_kb_uniform": float(np.exp(m0.w_ @ np.log(kbs))),
            "posterior_geo_kb_prior": float(np.exp(mc.w_ @ np.log(kbs)))}


@check("chloramine")
def pilot_log_species_never_mixed():
    """A7: the pilot's grab log takes free_chlorine_mgL or total_chlorine_mgL, never both: a mixed log is refused, and so
    is a log whose species does not match --disinfectant; docs/example_grab_log_total.csv loads as total chlorine and
    the committed free-chlorine example still loads as before.  Since the review: the pilot's default dose follows the
    disinfectant (1.2 free, 2.0 total chlorine), and its pH and Cl2:N prior inputs are refused with free chlorine or
    one without the other."""
    from .pilot import load_log
    taps = os.path.join(REPO, "docs", "example_tap_map.csv")
    free_csv, total_csv = os.path.join(REPO, "docs", "example_grab_log.csv"), os.path.join(REPO, "docs", "example_grab_log_total.csv")
    free = load_log(free_csv, taps)
    total = load_log(total_csv, taps, "chloramine")
    mixed = pd.read_csv(free_csv)
    mixed["total_chlorine_mgL"] = mixed["free_chlorine_mgL"]
    tmp = os.path.join(os.getcwd(), "mixed_log.csv")
    mixed.to_csv(tmp, index=False)
    refused = {}
    try:
        for name, args in (("mixed_free", (tmp, taps, "free_chlorine")), ("mixed_chloramine", (tmp, taps, "chloramine")),
                           ("total_as_free", (total_csv, taps, "free_chlorine")), ("free_as_chloramine", (free_csv, taps, "chloramine"))):
            try:
                load_log(*args)
            except ValueError:
                refused[name] = True
            else:
                refused[name] = False
    finally:
        os.remove(tmp)
    assert all(refused.values()), refused
    assert len(free) == 78 and len(total) == len(free) and list(free.columns) == list(total.columns)
    # review fixes: the pilot's dose defaults follow the disinfectant, and the chloramine prior's pH and Cl2:N go
    # together and only with chloramine (refused before any simulation)
    from .chloramine import CA_DOSE_MGL
    from .pilot import DEFAULT_DOSE_MGL, validate
    assert DEFAULT_DOSE_MGL == {"free_chlorine": 1.2, "chloramine": CA_DOSE_MGL}
    for kw in ({"ph": 8.0, "cl2n": 4.5}, {"ph": 8.0, "disinfectant": "chloramine"}):
        try:
            validate("Net3", free, **kw)
        except ValueError:
            refused[f"validate_{'_'.join(kw)}"] = True
        else:
            refused[f"validate_{'_'.join(kw)}"] = False
    assert all(refused.values()), refused
    return {"refused": refused, "total_log_rows": len(total), "total_log_mean_mgL": float(total.y.mean())}


@check("chloramine")
def not_testable_is_never_one():
    """A6: a rate with no junction below the threshold is null in the summary and printed 'not testable', never 1.0;
    the committed summaries hold no recall of exactly 1.0 for a cell with no true violation."""
    from .chloramine import _fmt, rates
    z = {"n_uns": 10, "sse": 0.1, "sae": 0.5, "sse_ln": 0.2, "se": 0.1, "in50": 5.0, "in80": 8.0, "in90": 9.0,
         "in95": 10.0, "n_true_viol": 0, "tp": 0, "fp": 2, "fn": 0}
    r = rates(z)
    assert math.isnan(r["recall"]) and _fmt(r["recall"]) == "not testable" and r["precision"] == 0.0
    r2 = rates({**z, "fp": 0})
    assert math.isnan(r2["precision"]) and math.isnan(r2["f1"])
    n_null = 0
    for net in CA_NETS:
        d = json.load(open(os.path.join(REPO, "outputs", "chloramine", f"summary_chloramine_{net}.json")))
        for v in d["pooled"].values():
            for m in v.values():
                for rr in m.values():
                    for cell in rr.values():
                        if cell["n_true_viol"] == 0:
                            assert cell["recall"] is None, cell
                            n_null += 1
    return {"recall_no_violation": _fmt(r["recall"]), "summary_cells_not_testable": n_null}


CA_ROWS_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap.chloramine import KW_REF, run_task
res = run_task("Net2", 300, "twin", KW_REF["Net2"], sys.argv[1])
print("CA_ROWS_JSON " + json.dumps([{k: (v if isinstance(v, str) else float(v)) for k, v in r.items()} for r in res["rows"]]))
print("CA_TRUTH_JSON " + json.dumps({k: v for k, v in res["truth"].items() if not isinstance(v, (dict, list))}, default=str))
"""


@check("chloramine", quick=False)
def chloramine_outputs_reproduce():
    """Task 12's committed outputs: each network's summary carries the acceptance key and is what --resummarise computes
    from the committed CSVs; every CSV is under 1 MB; the experimental label follows the bars; and recomputing Net2 seed
    300 under the first-order twin (model (b) and the oracle, both rules, n = 3, 8, 15) in a single-threaded subprocess
    gives the committed rows to the CSV's 6 significant digits."""
    from .chloramine import CSV_FLOAT, OUT_DIR as CA_OUT, _clean, summarise
    out = {}
    for net in CA_NETS:
        df = pd.read_csv(os.path.join(REPO, CA_OUT, f"results_time_{net}.csv"))
        tr = pd.read_csv(os.path.join(REPO, CA_OUT, f"truths_{net}.csv"))
        d = json.load(open(os.path.join(REPO, CA_OUT, f"summary_chloramine_{net}.json")))
        again = json.loads(json.dumps(_clean(summarise(net, df, tr))))
        assert again == d, f"{net}: the committed summary is not what the committed CSVs give"
        A = d["acceptance"]
        bars = [v["pass"] for k, v in A.items() if k[:2] in ("A1", "A2", "A3", "A4")]
        assert A["experimental"] == (not all(bars)), A
        sizes = {f: os.path.getsize(os.path.join(REPO, CA_OUT, f)) for f in (f"results_time_{net}.csv", f"truths_{net}.csv")}
        assert all(v < 1_000_000 for v in sizes.values()), sizes
        out[net] = {"experimental": A["experimental"], "csv_bytes": sizes}
    r = subprocess.run([PY, "-c", CA_ROWS_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV})
    line = [x for x in r.stdout.splitlines() if x.startswith("CA_ROWS_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"recompute failed (exit {r.returncode}): {r.stderr[-2000:]}")
    rows = pd.DataFrame(json.loads(line[-1][len("CA_ROWS_JSON "):]))
    want = pd.read_csv(os.path.join(REPO, CA_OUT, "results_time_Net2.csv"), float_precision="round_trip")
    want = want[(want.variant == "twin") & (want.seed == 300)]
    bad, n = {}, 0
    for _, g in rows.iterrows():
        w = want[(want.model == g.model) & (want.rule == g.rule) & (want.n == g.n)]
        assert len(w) == 1, (g.model, g.rule, g.n)
        for c in want.columns:
            if c in ("net", "variant", "model", "rule") or c not in g:
                continue
            a, b = g[c], w[c].iloc[0]
            n += 1
            if not ((pd.isna(a) and pd.isna(b)) or float(CSV_FLOAT % float(a)) == float(b)):
                bad[f"{g.model}.{g.rule}.{int(g.n)}.{c}"] = (a, b)
    assert len(rows) == len(want) and not bad, (len(rows), len(want), bad)
    return {**out, "net2_seed300_twin_values_compared": n, "identical_to_6_significant_digits": True}


APP_CA_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
at = AppTest.from_file(sys.argv[1], default_timeout=900).run()
assert not at.exception, [e.value for e in at.exception]
at.selectbox[0].set_value([o for o in at.selectbox[0].options if o.startswith("Net2")][0]).run()
[r for r in at.radio if r.label == "Disinfectant"][0].set_value("Chloramine (total chlorine)").run()
assert not at.exception, [e.value for e in at.exception]
metric = [m for m in at.metric if m.label.startswith("Junctions likely below")][0]
res = {"flagged_label": metric.label, "info": [i.value for i in at.info],
       "subheaders": [s.value for s in at.subheader], "captions": [c.value for c in at.caption],
       "number_inputs": [n.label for n in at.number_input]}
[t for t in at.toggle if t.label.startswith("Show the true daily-minimum map")][0].set_value(True)
at.run()
assert not at.exception, [e.value for e in at.exception]
res["reveal"] = [m.value for m in at.markdown if "True daily-minimum violations" in m.value]
[n for n in at.number_input if n.label.startswith("Minimum residual, total chlorine")][0].set_value(0.05)
at.run()
assert not at.exception, [e.value for e in at.exception]
res["reveal_no_violation"] = [m.value for m in at.markdown if "True daily-minimum violations" in m.value]
print("APP_CA_JSON " + json.dumps(res))
"""


@check("chloramine", quick=False)
def app_chloramine_mode():
    """Headless AppTest: on Net2 with the demo on, choosing 'Chloramine (total chlorine)' runs without exceptions; the
    threshold defaults to 0.5 mg/L total chlorine with its label (a common utility operating target, not a California
    rule); the mode's banner says total chlorine and carries the experimental label exactly when the committed summaries
    say so; the nitrification watch is shown with its 'not validated' label; the demo truth's reveal works, and with
    the threshold at 0.05 mg/L (no junction below it) it says recall is not testable instead of printing a rate (A6;
    a review fix).  The demo truth runs in a separate process.  The free chlorine default demo is checked by
    app_default_demo."""
    from .chloramine import THRESHOLD_NOTE
    r = subprocess.run([PY, "-c", APP_CA_SCRIPT, os.path.join(REPO, "app.py")], cwd=os.getcwd(), capture_output=True,
                       text=True, env={**os.environ, "PYTHONPATH": REPO})
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_CA_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-2000:]}")
    res = json.loads(line[-1][len("APP_CA_JSON "):])
    assert res["flagged_label"] == "Junctions likely below 0.5 mg/L (daily minimum)", res["flagged_label"]
    banner = [i for i in res["info"] if i.startswith("Chloramine mode")]
    assert banner and "TOTAL chlorine" in banner[0] and "not a California rule" in banner[0], res["info"]
    exp = any(json.load(open(os.path.join(REPO, "outputs", "chloramine", f"summary_chloramine_{n}.json")))["acceptance"]["experimental"]
              for n in CA_NETS)
    assert ("Experimental" in banner[0]) == exp, (banner[0], exp)
    assert "Nitrification watch (chloramine)" in res["subheaders"]
    assert any("not validated" in c for c in res["captions"])
    assert "Minimum residual, total chlorine (mg/L)" in res["number_inputs"] and "not a California rule" in THRESHOLD_NOTE
    assert res["reveal"], "the demo reveal did not render"
    nv = res["reveal_no_violation"]
    assert nv and "**0 of 35**" in nv[0] and "recall is not testable" in nv[0] and "%" not in nv[0], nv
    return {"flagged_label": res["flagged_label"], "experimental_label_shown": "Experimental" in banner[0],
            "reveal": res["reveal"][0][:160], "reveal_no_violation": nv[0][:160]}


# ----------------------------------------------------------------------------- task 13: the richer-truth audit
TRA_SMOKE_SEED = 990      # chemexp.SMOKE_SEEDS: never scored
TRA_NET2 = dict(kb_per_day=0.10, kw_m_per_day=0.20)     # experiment.NET_TRUTH['Net2']


def _two_ra_chain_run(temp_C: float, er_fast: float, er_slow: float, wall: bool = False):
    """12-pipe plug-flow chain (0.3 m pipes of 1000 m, 0.01 m3/s: 1.96 h per pipe) under the 2RA MSX model at temp_C,
    the rate constants scaled by the Arrhenius factor; returns (MSX last-hour chlorine at J1..J12, ages in hours)."""
    from .chemistry import arrhenius
    from . import msx as M
    wn, t_days = _chain(quality_step=300)
    wn.options.quality.parameter = "NONE"
    kf = M.KF20_L_PER_MG_H * arrhenius(temp_C, er_fast)
    ks = M.KS20_L_PER_MG_H * arrhenius(temp_C, er_slow)
    m = M.two_reactant_model(wn, None, kf, ks)
    f0, s0 = M.F0_GREENVALE_MGL * 0.304, M.S0_GREENVALE_MGL * 0.304
    M.add_sources(m, wn, {"R": {"CL2": 1.2, "FAST": f0, "SLOW": s0}})
    res, info = M.run(wn, m)
    last = res.node["CL2"].iloc[-1]
    return np.array([float(last[f"J{k}"]) for k in range(1, 13)]), t_days * 24.0, (kf, ks, f0, s0), info


@check("two_reactant", quick=False)
def two_reactant_msx_matches_lsoda():
    """The 2RA reactions in EPANET-MSX (msx.two_reactant_model, no wall) on a 12-pipe plug-flow chain against the scipy
    LSODA batch solution (msx.two_reactant_batch) at every junction's travel time, at 10, 20 and 25 C (Greenvale's
    constants times the Arrhenius factor with E/R 6500 K fast and 10000 K slow, inside the plan's ranges; the audit's
    truth runs at 20 C only): within 0.01 mg/L.  The source column reads the 1.2 mg/L dose, compiled reactions ran,
    and no MSX scratch file is left."""
    from . import msx as M
    out, worst = {}, 0.0
    for T in (10.0, 20.0, 25.0):
        sim, ages, (kf, ks, f0, s0), info = _two_ra_chain_run(T, 6500.0, 10000.0)
        ref = M.two_reactant_batch(1.2, f0, s0, ages, kf, ks)[0]
        d = float(np.abs(sim - ref).max())
        worst = max(worst, d)
        out[f"{T:g}C"] = {"max_abs_diff_mgL": d, "C_last_msx": float(sim[-1]), "C_last_lsoda": float(ref[-1]),
                          "age_last_h": float(ages[-1]), "compiler": info["compiler"]}
    assert worst <= CHAIN_TOL_MGL, out
    assert all(v["compiler"] == "GC" for v in out.values()), out
    assert _root_scratch_files() == [] and not [f for f in os.listdir(".") if re.fullmatch(r"(msx|en)[A-Za-z0-9]{6}", f)]
    return {"by_temperature": out, "worst_mgL": worst, "bar_mgL": CHAIN_TOL_MGL}


@check("two_reactant", quick=False)
def first_order_msx_wall_matches_epanet_at_temperatures():
    """The plan's second engine check: one species with EPANET's first-order bulk (0.40 per day) and mass-transfer-
    limited wall (0.70 m/day x the roughness factor) as MSX expressions against EPANET's CHEMICAL quality on Net3, with
    the water's viscosity and diffusivity set for 10 and 25 C in both engines (chemistry.apply_water_temperature; the
    20 C case is task 12's msx_first_order_matches_epanet), every junction and hour of the last day, at a 60 s quality
    step (task 12's finding: at 300 s the two engines differ by about 0.03 mg/L at passing fronts).  Bar 0.03 mg/L.
    Its first run failed at 10 C: MSX computes its Reynolds number with water's 20 C viscosity, whatever the file says.
    msx.add_mass_transfer_terms now divides Re by the file's relative viscosity when that is not 1; at 1 (every
    committed MSX run: Net3, Net2 and ky4) its expressions are as before, checked here.  The uncorrected runner is
    rerun too (its term ReV set back to MSX's own Re), so the gap the fix closes is recorded (before_fix_*), and the
    check asserts that the uncorrected runner fails the bar at 10 C (added after task 13's review)."""
    from . import chloramine as C
    from . import msx as M
    from . import simulate as S
    from .chemistry import apply_water_temperature
    kb, kw, dose, step = 0.40, 0.70, 1.2, 60
    out = {}
    for T in (10.0, 25.0):
        old = S.QUALITY_STEP_S
        S.QUALITY_STEP_S = step
        try:
            ref = S.simulate_nominal_chlorine("Net3", kb, kw, 1.0, dose, file_prefix=os.path.join(os.getcwd(), "fo_t"),
                                              temp_C=T)
        finally:
            S.QUALITY_STEP_S = old
            for f in glob.glob("fo_t.*"):
                os.remove(f)
        diffs = {}
        for fixed in (True, False):
            wn = S.load("Net3")
            wn.options.time.quality_timestep = step
            v, d = apply_water_temperature(wn, T)
            m = C.first_order_msx_model(wn, kb, {pn: kw * S.roughness_factor(p.roughness, 1.0) for pn, p in wn.pipes()})
            m.options.timestep = step
            if not fixed:                  # the runner before the fix: MSX's own Re, whatever the file's viscosity
                m.reaction_system.terms["ReV"].expression = "Re"
            M.add_sources(m, wn, {s: {"CL": dose} for s in S.source_nodes(wn)})
            res, _ = M.run(wn, m)
            q = S._last_day(res.node["CL"], wn.junction_name_list)
            diffs[fixed] = (q - ref).abs()
        out[f"{T:g}C"] = {"viscosity_ratio": v, "diffusivity_ratio": d, "max_abs_diff_mgL": float(diffs[True].values.max()),
                          "mean_abs_diff_mgL": float(diffs[True].values.mean()),
                          "before_fix_max_abs_diff_mgL": float(diffs[False].values.max()),
                          "before_fix_mean_abs_diff_mgL": float(diffs[False].values.mean())}
    assert all(x["max_abs_diff_mgL"] <= 0.03 for x in out.values()), out
    assert out["10C"]["before_fix_max_abs_diff_mgL"] > 0.03, out     # the check catches the fault it found
    wn = S.load("Net3")
    m = C.first_order_msx_model(wn, kb, {pn: kw for pn in wn.pipe_name_list})
    terms = {t: m.reaction_system.terms[t].expression for t in m.term_name_list}
    assert "ReV" not in terms and terms["Yg"] == "D/Len*Re*Sc", terms
    out["relative_viscosity_1_expressions_unchanged"] = True
    return out


@check("two_reactant")
def two_reactant_scale_reproduces():
    """s is what brentq gives (the 2RA batch's apparent first-order rate at the match age equals the committed truth
    rate to 1e-9 per day: Net3 0.40 at 24 h, Net2 0.10 at 24 h and at 96 h), and the committed
    outputs/chem_2ra/scaling.json is exactly what chemexp.scaling() writes."""
    from . import msx as M
    from .chemexp import _clean, scaling
    path = os.path.join(REPO, "outputs", "chem_2ra", "scaling.json")
    committed = json.load(open(path))
    again = json.loads(json.dumps(_clean(scaling())))
    assert again == committed, "outputs/chem_2ra/scaling.json is not what chemexp.scaling() gives"
    out = {}
    for net, ent in committed["networks"].items():
        for truth, e in ent.items():
            s = M.two_reactant_scale(e["kb_matched_per_day"], e["match_h"])
            k = M.two_reactant_apparent_rate(s, e["match_h"])
            assert s == e["scale_s"] and abs(k - e["kb_matched_per_day"]) < 1e-9, (net, truth, s, k)
            out[f"{net}/{truth}"] = {"scale_s": s, "apparent_rate_at_match": k}
    return out


@check("two_reactant", quick=False)
def two_reactant_truth_draws_and_refusals():
    """The 2RA truth (Net2, smoke seed 990, never scored) consumes the committed draws in their committed order: its
    source doses and monthly bulk factor equal the free-chlorine chemistry truth's at the same seed (Chemistry(temp_C=20),
    which consumes them in the default branch's order), and its reactant loads are 1.13 s u and 2.87 s u.  Compiled
    reactions ran, wn.msx is detached, no MSX scratch file is left.  Refused: a 2RA truth at 10 C, with a TOC input,
    as chloramine, with in-network warming or the loss split, as a grid condition, and two_reactant= without kinetics
    '2ra'."""
    from . import msx as M
    from .chemistry import CHLORAMINE, Chemistry, Warming
    from .simulate import build_scenario
    tr = build_scenario("Net2", TRA_SMOKE_SEED, chem=Chemistry(kinetics="2ra"), **TRA_NET2)
    fo = build_scenario("Net2", TRA_SMOKE_SEED, chem=Chemistry(temp_C=20.0), **TRA_NET2)
    assert tr.chem["source_doses_mgL"] == fo.chem["source_doses_mgL"]
    assert tr.chem["bulk_month_factor"] == fo.chem["bulk_month_factor"]
    s = M.two_reactant_scale(0.10)
    u = tr.chem["bulk_month_factor"]
    assert tr.chem["scale_s"] == s and np.isclose(tr.chem["fast_mgL"], 1.13 * s * u, rtol=1e-12)
    assert np.isclose(tr.chem["slow_mgL"], 2.87 * s * u, rtol=1e-12)
    assert tr.chem["msx_compiler"] == "GC" and tr.chem["msx_fallback_reason"] is None and tr.wn.msx is None
    assert (tr.truth_daily_min >= 0).all()
    refused = {}
    cases = {"10C": lambda: Chemistry(kinetics="2ra", temp_C=10.0).require_built(),
             "toc": lambda: Chemistry(kinetics="2ra", toc_mgL=2.5).require_built(),
             "chloramine": lambda: Chemistry(disinfectant=CHLORAMINE, kinetics="2ra").require_built(),
             "grid": lambda: Chemistry(kinetics="2ra").sim_kwargs(),
             "warming": lambda: build_scenario("Net2", TRA_SMOKE_SEED, chem=Chemistry(kinetics="2ra", temp_C=20.0),
                                               warming=Warming(soil_temp_C=18.0), **TRA_NET2),
             "loss_split": lambda: build_scenario("Net2", TRA_SMOKE_SEED, chem=Chemistry(kinetics="2ra"),
                                                  truth_loss_split=True, **TRA_NET2),
             "two_reactant_without_2ra": lambda: build_scenario("Net2", TRA_SMOKE_SEED, two_reactant=M.TwoReactantTruth(),
                                                                **TRA_NET2)}
    for k, fn in cases.items():
        try:
            fn()
            refused[k] = False
        except (ValueError, NotImplementedError):
            refused[k] = True
    assert all(refused.values()), refused
    assert _root_scratch_files() == [] and not [f for f in os.listdir(".") if re.fullmatch(r"(msx|en)[A-Za-z0-9]{6}", f)]
    return {"scale_s": s, "bulk_month_factor": u, "refused": refused,
            "median_daily_min_2ra": float(tr.truth_daily_min.median()), "median_daily_min_first_order": float(fo.truth_daily_min.median())}


_PAIRED_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap import experiment as E
from residualmap.chemexp import scenario
from residualmap.features import build_features
from residualmap.simgp import SimGP24
cache = sys.argv[1]
seen = {}
class Rec(SimGP24):
    def fit(self, S):
        seen.setdefault(KEY, []).append([[str(j), int(h), float(y)] for j, h, y in zip(S.junction, S.hour, S.y)])
        return super().fit(S)
def mk(sc, X, seed, cache_dir, threshold):
    return Rec(sc, X, seed=seed, cache_dir=cache_dir, threshold=threshold)
truths = {}
for truth, bl in (("first_order", True), ("first_order", False), ("2ra", False)):
    sc = scenario("Net2", 990, truth)
    truths[truth] = {f"{h}|{j}": float(sc.truth_by_hour.loc[h, j]) for h in range(24) for j in sc.junctions}
    KEY = f"{truth}|{bl}"
    E.run_scenario_time(sc, build_features(sc), n_max=6, seed=990, cache_dir=cache, make_model=mk, baselines=bl, extra=True)
print("PAIRED_JSON " + json.dumps({"seen": seen, "truths": truths}))
'''


@check("two_reactant", quick=False)
def paired_samples_are_the_same_across_truths():
    """The pairing the materiality test rests on: experiment.run_scenario_time on Net2 smoke seed 990 (n up to 6), once
    on the committed first-order truth and once on the 2RA truth, gives the random rule the same junctions, hours and
    reading errors (reading minus the truth at that junction and hour, where the reading is not clipped at 0.01), and
    the same 3 first samples to every rule; skipping the random-rule baselines (baselines=False, the audit's opt-in
    variants) changes no sample.  Run in a single-threaded subprocess."""
    r = subprocess.run([PY, "-c", _PAIRED_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV}, timeout=900)
    line = [x for x in r.stdout.splitlines() if x.startswith("PAIRED_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"paired-sample run failed (exit {r.returncode}): {r.stderr[-2000:]}")
    d = json.loads(line[-1][len("PAIRED_JSON "):])
    seen, truths = d["seen"], d["truths"]
    fo_b, fo, tr = seen["first_order|True"], seen["first_order|False"], seen["2ra|False"]
    assert fo_b == fo, "baselines=False changed the samples"
    n_fits_rule = len(fo) // 4                                    # four rules, in experiment.TIME_STRATEGIES order
    rnd_fo, rnd_tr = fo[:n_fits_rule], tr[:n_fits_rule]
    errs = []
    for a, b in zip(rnd_fo, rnd_tr):
        assert [x[:2] for x in a] == [x[:2] for x in b], "the random rule's junctions or hours differ between truths"
        for (j, h, ya), (_, _, yb) in zip(a, b):
            if ya > 0.01 and yb > 0.01:
                ea, eb = ya - truths["first_order"][f"{h}|{j}"], yb - truths["2ra"][f"{h}|{j}"]
                errs.append(abs(ea - eb))
    firsts = {tuple(map(tuple, (x[:2] for x in seen[k][i * n_fits_rule][:3]))) for k in ("first_order|False", "2ra|False")
              for i in range(4)}
    assert len(firsts) == 1, "the 3 first samples differ between rules or truths"
    assert errs and max(errs) < 1e-9, max(errs) if errs else None
    return {"random_rule_fits_compared": len(rnd_fo), "reading_errors_compared": len(errs),
            "max_reading_error_difference": max(errs)}


@check("two_reactant", quick=False)
def lowkb_grid_nests_the_committed_grid():
    """The opt-in low-kb grid (simgp.GRIDS['full_lowkb']): its members at today's five bulk rates are the committed
    'full' grid's members bit for bit (Net3 and Net2), in grid order; its cache file carries a tag; the committed grids'
    axes, cache names and the committed condition tag are unchanged; and the default SimGP24 still runs on 'full'."""
    import inspect
    from .chemistry import Chemistry, cache_tag
    from .simgp import GRIDS, KB_GRID, KB_LOWKB, SimGP24, grid_cache_path
    from .simulate import nominal_scenario
    assert GRIDS["full"][0] == [0.10, 0.25, 0.40, 0.55, 0.70] and KB_LOWKB == [0.0125, 0.025, 0.05, 0.075] + KB_GRID
    assert inspect.signature(SimGP24.__init__).parameters["grid"].default == "full"
    assert cache_tag(Chemistry(temp_C=12.5, er_K=8000.0)) == "72f9c19f47"
    out = {}
    for net in ("Net3", "Net2"):
        sc = nominal_scenario(net)
        assert os.path.basename(grid_cache_path(sc)) == f"grid24_full_{net}.pkl"
        lp = grid_cache_path(sc, CACHE, "full_lowkb")
        assert re.search(r"grid24_full_lowkb_%s_[0-9a-f]{10}\.pkl$" % net, lp), lp
        assert os.path.exists(lp), f"{lp} is built by the audit (python -m residualmap.chemexp audit {net} 8)"
        P0, Z0 = pickle.load(open(grid_cache_path(sc, CACHE), "rb"))
        P1, Z1 = pickle.load(open(lp, "rb"))
        keep = [i for i, p in enumerate(P1) if p[0] in KB_GRID]
        assert [tuple(P1[i]) for i in keep] == [tuple(p) for p in P0]
        assert np.array_equal(Z1[keep], Z0), net
        out[net] = {"members": len(P1), "members_at_committed_rates_identical": len(keep), "file": os.path.basename(lp)}
    return out


@check("two_reactant", quick=False)
def two_rate_grid_is_exact_and_nests_today():
    """The gated 'full2r' grid (chemexp.two_rate_grid) from Net3's committed grid: its first 675 members are today's
    grid (B0 nested), it is in grid order (simgp.check_grid_order), and a two-rate member is ln(f C(k1) + (1 - f) C(k2))
    of the committed members.  The exactness it rests on, checked against an engine: an MSX run carrying two chlorine
    parts (dose f x 1.2 at bulk rate k1, (1 - f) x 1.2 at k2, the same mass-transfer-limited wall) on Net3 sums to the
    mixture of two EPANET runs at k1 and k2 within 0.03 mg/L (60 s quality step)."""
    from . import msx as M
    from . import simulate as S
    from .chemexp import two_rate_grid
    from .simgp import check_grid_order, simulator_grid_24h
    from .simulate import nominal_scenario
    sc = nominal_scenario("Net3")
    P, Z = simulator_grid_24h(sc, CACHE, "full")
    members, P2, Z2 = two_rate_grid(P, Z)
    assert len(P2) == 4725 and np.array_equal(Z2[:675], Z) and [tuple(p) for p in P2[:675]] == [tuple(p) for p in P]
    check_grid_order(P2, 9)
    i = members.index((0.55, 0.10, 0.25, 0.60, 1.0, 1.0, 1.0))
    a, b = P.index((0.55, 0.60, 1.0, 1.0, 1.0)), P.index((0.10, 0.60, 1.0, 1.0, 1.0))
    mix = np.log(0.25 * np.exp(Z[a].astype(np.float64)) + 0.75 * np.exp(Z[b].astype(np.float64)))
    assert np.allclose(Z2[i], mix, atol=1e-6)
    # engine check of the linearity
    from wntr.msx import MsxModel
    k1, k2, f, kw, dose, step = 0.55, 0.10, 0.25, 0.60, 1.2, 60
    runs = {}
    for k in (k1, k2):
        old = S.QUALITY_STEP_S
        S.QUALITY_STEP_S = step
        try:
            runs[k] = S.simulate_nominal_chlorine("Net3", k, kw, 1.0, dose, file_prefix=os.path.join(os.getcwd(), "tr_m"))
        finally:
            S.QUALITY_STEP_S = old
            for fn in glob.glob("tr_m.*"):
                os.remove(fn)
    wn = S.load("Net3")
    wn.options.time.quality_timestep = step
    m = MsxModel()
    m.options.rate_units, m.options.area_units, m.options.solver, m.options.timestep = "SEC", "M2", "RK5", step
    for sp in ("CA", "CB"):
        m.add_species(sp, "bulk", units="MG", atol=1e-8, rtol=1e-6)
    m.add_constant("kA", k1 / 86400.0)
    m.add_constant("kB", k2 / 86400.0)
    M.add_mass_transfer_terms(m, wn)
    m.add_parameter("kwp", 0.0)
    M.set_pipe_parameter(m, "kwp", {pn: M.wall_rate_unit(wn, kw * S.roughness_factor(p.roughness, 1.0)) for pn, p in wn.pipes()})
    m.add_reaction("CA", "pipe", "rate", "-kA*CA - (4/D)*kwp*Kf/(kwp+Kf)*CA")
    m.add_reaction("CB", "pipe", "rate", "-kB*CB - (4/D)*kwp*Kf/(kwp+Kf)*CB")
    m.add_reaction("CA", "tank", "rate", "-kA*CA")
    m.add_reaction("CB", "tank", "rate", "-kB*CB")
    M.add_sources(m, wn, {s: {"CA": f * dose, "CB": (1 - f) * dose} for s in S.source_nodes(wn)})
    res, _ = M.run(wn, m)
    tot = S._last_day(res.node["CA"] + res.node["CB"], wn.junction_name_list)
    want = f * runs[k1] + (1 - f) * runs[k2]
    d = (tot - want).abs()
    assert float(d.values.max()) <= 0.03, float(d.values.max())
    return {"members": len(P2), "msx_two_parts_vs_epanet_mixture_max_abs_mgL": float(d.values.max()),
            "mean_abs_mgL": float(d.values.mean())}


@check("two_reactant")
def default_path_untouched_by_task13_hooks():
    """Task 13's hooks keep the committed defaults: run_scenario_time's make_model None, baselines True and extra
    False; _metrics_time with extra False returns exactly the committed keys; experiment.main's chemistry None and
    figures True; the chemistry runs' folders.  (The byte-identical reruns are in baseline_reproduction.json.)"""
    import inspect
    from . import experiment as E
    p = inspect.signature(E.run_scenario_time).parameters
    assert p["make_model"].default is None and p["baselines"].default is True and p["extra"].default is False
    m = inspect.signature(E.main).parameters
    assert m["chemistry"].default is None and m["figures"].default is True and m["outdir"].default == "outputs"
    assert inspect.signature(E._metrics_time).parameters["extra"].default is False
    assert E.chem_outdir("outputs", "2ra") == os.path.join("outputs", "chem_2ra")
    assert E.chem_outdir("outputs", "2ra", 96) == os.path.join("outputs", "chem_2ra", "match96")
    assert E.chem_outdir("outputs", "first_order") == os.path.join("outputs", "chem_2ra", "first_order")
    import types
    sc = types.SimpleNamespace(truth_daily_min=pd.Series([0.1, 0.3, 0.5], index=["a", "b", "c"]))
    pm = pd.DataFrame({"median": [0.15, 0.25, 0.6], "lo90": [0.05, 0.2, 0.4], "hi90": [0.3, 0.4, 0.8],
                       "p_below": [0.9, 0.2, 0.0]}, index=["a", "b", "c"])
    base = E._metrics_time(sc, pm, None, ["a", "b", "c"])
    assert list(base) == ["rmse_min", "precision_min", "recall_min", "f1_min", "n_true_viol_min", "coverage90_min"]
    ext = E._metrics_time(sc, pm, None, ["a", "b", "c"], extra=True)
    assert {k: ext[k] for k in base} == base and ext["tp_min"] == 1 and ext["fp_min"] == 0
    return {"committed_keys": list(base), "extra_keys": [k for k in ext if k not in base]}


_TRA_ROWS_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from residualmap.chemexp import run_task, scenario
from residualmap.experiment import run_scenario_time
from residualmap.features import build_features
cache = sys.argv[1]
sc = scenario("Net2", 400, "2ra")
df, _ = run_scenario_time(sc, build_features(sc), seed=400, cache_dir=cache, extra=True)
res = run_task("Net2", 400, "2ra", cache, ("oracle", "lowkb"))
def rows(d):
    return [{k: (v if isinstance(v, str) else (None if v != v else float(v))) for k, v in r.items()} for r in d]
print("TRA_JSON " + json.dumps({"time": rows(df.to_dict("records")), "audit": rows(res["audit"]),
                                "lowkb": rows(res["lowkb"].to_dict("records"))}))
'''


def _same6(a, b) -> bool:
    if (a is None or (isinstance(a, float) and a != a)) and pd.isna(b):
        return True
    return float("%.6g" % float(a)) == float("%.6g" % float(b))


@check("two_reactant", quick=False)
def two_reactant_outputs_reproduce():
    """Task 13's committed outputs: each network's summary_audit is what chemexp.summarise computes from the committed
    CSVs as written; every task-13 CSV is under 1 MB; and recomputing Net2 seed 400 under the 2RA truth in a
    single-threaded subprocess gives the committed rows of results_time (every model and rule, n = 3 to 15), of the
    oracle and of the opt-in low-kb grid to 6 significant digits."""
    from .chemexp import _clean, summarise
    out = {}
    for net in ("Net3", "Net2"):
        d = json.load(open(os.path.join(REPO, "outputs", "chem", f"summary_audit_{net}.json")))
        old = os.getcwd()
        os.chdir(REPO)                 # summarise reads outputs/ by relative path; it writes nothing
        try:
            again = json.loads(json.dumps(_clean(summarise(net))))
        finally:
            os.chdir(old)
        assert again == d, f"{net}: the committed summary_audit is not what the committed CSVs give"
        out[net] = {"material": d["acceptance"]["materiality"]["material"]}
    sizes = {os.path.relpath(f, REPO): os.path.getsize(f) for f in
             glob.glob(os.path.join(REPO, "outputs", "chem_2ra", "**", "*.csv"), recursive=True)
             + glob.glob(os.path.join(REPO, "outputs", "chem", "*2ra*.csv")) + glob.glob(os.path.join(REPO, "outputs", "chem", "audit_*.csv"))}
    assert sizes and all(v < 1_000_000 for v in sizes.values()), sizes
    r = subprocess.run([PY, "-c", _TRA_ROWS_SCRIPT, CACHE], cwd=os.getcwd(), capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": REPO, **SINGLE_THREAD_ENV}, timeout=1800)
    line = [x for x in r.stdout.splitlines() if x.startswith("TRA_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"recompute failed (exit {r.returncode}): {r.stderr[-2000:]}")
    got = json.loads(line[-1][len("TRA_JSON "):])
    n, bad = 0, {}
    want_t = pd.read_csv(os.path.join(REPO, "outputs", "chem_2ra", "results_time_Net2.csv"), float_precision="round_trip")
    want_t = want_t[want_t.seed == 400]
    assert len(got["time"]) == len(want_t), (len(got["time"]), len(want_t))
    for g, (_, w) in zip(got["time"], want_t.iterrows()):
        assert (g["model"], g["strategy"], g["n"]) == (w.model, w.strategy, w.n)
        for c in want_t.columns:
            if c in ("model", "strategy"):
                continue
            n += 1
            if not _same6(g.get(c), w[c]):
                bad[f"time.{g['model']}.{g['strategy']}.{int(w.n)}.{c}"] = (g.get(c), w[c])
    A = pd.read_csv(os.path.join(REPO, "outputs", "chem", "audit_Net2.csv"), float_precision="round_trip")
    A = A[(A.seed == 400) & (A.truth == "2ra") & (A.part == "oracle")]
    for g in got["audit"]:
        w = A[A.model == g["model"]]
        assert len(w) == 1, g["model"]
        for c in A.columns:
            if c in ("net", "truth", "part", "model") or c not in g:
                continue
            n += 1
            if not _same6(g[c], w[c].iloc[0]):
                bad[f"oracle.{g['model']}.{c}"] = (g[c], w[c].iloc[0])
    L = pd.read_csv(os.path.join(REPO, "outputs", "chem_2ra", "results_lowkb_Net2.csv"), float_precision="round_trip")
    L = L[(L.seed == 400) & (L.truth == "2ra")]
    assert len(got["lowkb"]) == len(L)
    for g, (_, w) in zip(got["lowkb"], L.iterrows()):
        for c in L.columns:
            if c in ("model", "strategy", "truth"):
                continue
            n += 1
            if not _same6(g.get(c), w[c]):
                bad[f"lowkb.{g['strategy']}.{int(w.n)}.{c}"] = (g.get(c), w[c])
    assert not bad, bad
    return {**out, "csv_bytes_max": max(sizes.values()), "n_csv": len(sizes),
            "net2_seed400_2ra_values_compared": n, "identical_to_6_significant_digits": True}


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


AGE_NETS = ("Net3", "Net2", "ky4")


@check("anchors", quick=False)
def water_age_outputs_reproduce():
    """Task 9's committed outputs: water_age_<net>.json for Net3, Net2 and ky4 carry rmse_h, mae_h, band_coverage and
    n_excluded, and call the range a 90% band only if every coverage in the file is at least 0.85; error_by_age_<net>.csv
    exists for Net3 and Net2.  Rerunning Net3 seed 0 under the default truth gives the committed CSV row exactly, and
    that row is the app's default demo: 36 of 36 unsampled violations found, 1 false alarm."""
    from .age import BAND_COVERAGE_BAR, RANGE_LABEL
    from .experiment_chem import run_age_seed
    out = {}
    for net in AGE_NETS:
        d = json.load(open(os.path.join(OUT_DIR, f"water_age_{net}.json")))
        assert all(k in d for k in ("rmse_h", "mae_h", "band_coverage", "n_excluded")), net
        assert all(k in d["by_truth"]["default"] for k in ("converged_only", "lower_bounds")), net
        jf = pd.read_csv(os.path.join(OUT_DIR, f"water_age_junctions_{net}.csv"))
        assert int(jf.lower_bound_daily_max.sum()) == d["nominal_age"]["n_lower_bound_daily_max"], net
        covs = [v[k] for v in d["by_truth"].values() for k in ("band_coverage", "band_coverage_daily_max")]
        assert d["band_label"] == ("90% band" if min(covs) >= BAND_COVERAGE_BAR else RANGE_LABEL), (net, d["band_label"], covs)
        out[net] = {**{k: d[k] for k in ("rmse_h", "band_coverage", "n_excluded", "band_label")},
                    "converged_only_rmse_h": d["by_truth"]["default"]["converged_only"]["rmse_h"]}
    for net in ("Net3", "Net2"):
        assert os.path.exists(os.path.join(OUT_DIR, f"error_by_age_{net}.csv")), net
    committed = pd.read_csv(os.path.join(OUT_DIR, "water_age_Net3.csv"), float_precision="round_trip")
    want = committed[(committed.truth == "default") & (committed.seed == 0)].iloc[0]
    row, _, _ = run_age_seed("Net3", 0, "default", None, CACHE)
    diff = {k: (row[k], want[k]) for k in row if k in want.index and isinstance(row[k], (int, float))
            and not (pd.isna(want[k]) and pd.isna(row[k])) and row[k] != want[k]}
    assert not diff, diff
    assert (row["n_true_viol_min"], row["recall_min"], row["false_alarms_min"]) == (36, 1.0, 1), row
    return {"summaries": out, "net3_seed0_row_identical": True}


APP_AGE_SCRIPT = r'''
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
app = sys.argv[1]
res = {}
for tag in ("Net3", "Net2", "demo_off"):
    at = AppTest.from_file(app, default_timeout=600).run()
    if tag == "Net2":
        at.selectbox[0].set_value([o for o in at.selectbox[0].options if o.startswith("Net2")][0]).run()
    if tag == "demo_off":
        [t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False).run()
    assert not at.exception, (tag, [e.value for e in at.exception])
    subs = [s.value for s in at.subheader]
    caps = [c.value for c in at.caption]
    route_cols = [list(d.value.columns) for d in at.dataframe if "reason" in d.value.columns][0]
    res[tag] = {"water_age_panel": "Water age" in subs, "loss_panel": "Where chlorine is lost" in subs,
                "oldest_water": any("oldest water:" in c for c in caps),
                "tested_coverage_in_caption": any("In simulation on this network it held" in c for c in caps),
                "loss_caption": any(c.startswith("For the single most likely decay rates") for c in caps),
                "lower_bound_note": any("its age is a lower bound" in c for c in caps),
                "needs_sample_note": any("Enter at least one grab sample" in i.value for i in at.info),
                "route_age_column": "water age (h)" in route_cols,
                "pdf_button": len(at.get("download_button")) == 1}
print("APP_AGE_JSON " + json.dumps(res))
'''


@check("anchors", quick=False)
def app_water_age_panels():
    """Headless AppTest of task 9's app changes: on Net3 and Net2 (demo on) and with the demo off (no samples) the app
    runs without exceptions and shows the 'Water age' and 'Where chlorine is lost' panels, the oldest-water line, the
    route table's water-age column and the PDF button; the loss split needs a sample, so with none the panel says so;
    on the example networks the caption quotes the simulated coverage from outputs/chem."""
    r = subprocess.run([PY, "-c", APP_AGE_SCRIPT, os.path.join(REPO, "app.py")], cwd=os.getcwd(), capture_output=True,
                       text=True, env={**os.environ, "PYTHONPATH": REPO})
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_AGE_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-2000:]}")
    res = json.loads(line[-1][len("APP_AGE_JSON "):])
    for tag in ("Net3", "Net2", "demo_off"):
        g = res[tag]
        assert g["water_age_panel"] and g["loss_panel"] and g["oldest_water"] and g["route_age_column"] and g["pdf_button"], (tag, g)
        assert g["lower_bound_note"], (tag, g)
        assert g["tested_coverage_in_caption"], (tag, g)
    assert res["Net3"]["loss_caption"] and res["Net2"]["loss_caption"] and not res["Net3"]["needs_sample_note"]
    assert res["demo_off"]["needs_sample_note"] and not res["demo_off"]["loss_caption"]
    return res


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
