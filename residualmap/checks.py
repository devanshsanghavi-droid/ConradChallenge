"""
checks.py: the saved checks for the chemistry work (iteration 4, journal tasks 8 to 14) and for iteration 5's inputs
and access (task 15 on: groups task15, task15b and task16, files built by residualmap/fixtures.py).

    python -m residualmap.checks            # every check, about 990 s on this machine after task 16 (about 530 s
                                            # after task 14; 69 checks after task 14: 68 in the committed report,
                                            # plus the names check, skipped without the variable; 76 after task 15;
                                            # 102 after task 16 and its review, 101 in the committed report;
                                            # 103 after its final verification, 102 in the report; each run's
                                            # time is in outputs/chem/checks_run.json)
    python -m residualmap.checks --quick    # skips the fresh grid, the synthetic pilot, the app and the slow
                                            # seasonal, organics and chloramine (MSX) checks (about 50 s after task 16)
    RESIDUALMAP_NAMES_CSV=<private csv> python -m residualmap.checks    # also runs the names check (task 14)

Task 14's app check reads the app's PDF with pypdf (in requirements.txt); without pypdf it still requires the PDF
and skips only its PDF-text assertions, which the run log records (pdf_text_checked).  Its names check reads a private
list kept outside the repository (RESIDUALMAP_NAMES_CSV) and is skipped without it; its result goes to the git-ignored
run log only (RUN_LOG_ONLY), so the committed report is the same whether or not the list is on this machine, and a full
run with or without the variable leaves it unchanged.  A names-check failure still fails the run.

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
import types

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
FRESH_GRID_SECONDS: list[float] = []   # the fresh Net3 grid's build time in this run (run_cost_gate compares, task 15b)


# ----------------------------------------------------------------------------- tiny harness
CHECKS: list[tuple[str, str, bool, callable]] = []


class SkipCheck(Exception):
    """Raised by a check that cannot run here (task 14: a private input that is not on this machine); recorded as
    skipped with its reason, like a --quick skip, never as a pass."""


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
    FRESH_GRID_SECONDS.append(secs)
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


# ----------------------------------------------------------------------------- task 14: app, capabilities, docs, numbers
DASH = ("—", "–")       # em dash, en dash
# iteration 5's new files, held to 'no em or en dashes' whole (task 15 on)
IT5_FILES = ("residualmap/importers.py", "residualmap/fixtures.py", ".streamlit/config.toml",
             "residualmap/intake_exp.py",                                        # task 15b
             "pages/1_Get_started.py", "docs/how_to_export_your_model.md",       # task 16
             "docs/templates/grab_log_template.csv", "docs/templates/grab_log_total_template.csv",
             "docs/templates/tap_map_by_id_template.csv", "docs/templates/plant_log_template.csv")
# (em, en) dash counts at commit 1fcf39a (before task 14): task 14 adds none to these files
DASHES_BEFORE_TASK14 = {"docs/pilot_protocol.md": (21, 6), "app.py": (10, 1), "docs/feature_dictionary.md": (3, 1),
                        "README.md": (27, 28)}
ACCOUNTS_ROWS = ["Water age", "Organics and chlorine demand", "Chlorine type", "Temperature"]
NAMES_ENV = "RESIDUALMAP_NAMES_CSV"   # a CSV kept outside the repository with an 'org' column (and optionally
                                      # 'system_name', 'type'): names that must not appear in any tracked text file
# checks recorded in outputs/chem/checks_run.json (git-ignored) only, never in the committed report, with the reason
RUN_LOG_ONLY = {"no_listed_names_in_tracked_files":
                f"reads a private list kept outside the repository ({NAMES_ENV}) and is skipped without it, so its result "
                f"is in outputs/chem/checks_run.json only; a failure still fails the run"}


def _norm_ws(s: str) -> str:
    return " ".join(s.split())


@check("task14")
def capabilities_table_reads_outputs():
    """The app's 'What this model accounts for' table (capabilities.accounts_for): the manager's four points in order,
    every source in capabilities.SOURCES resolves, the numbers it prints are the outputs' values at the printed
    precision, every 'tested' cell says simulation, the adoption words follow the outputs (temperature and organics
    tested and not adopted; chloramine labelled experimental exactly when task 12's summaries say so), and no cell has
    an em or en dash.  The PDF's chemistry lines (capabilities.report_lines) name the disinfectant, the measured
    species, the threshold and its source, say the seasons are absorbed by refitting and 'Validated in simulation
    only', and carry no kb20 or temperature field (no temperature model was adopted)."""
    from .capabilities import COLUMNS, SOURCES, accounts_for, numbers, report_lines, resolve
    for name in SOURCES:
        resolve(name)
    n = numbers()
    t = accounts_for()
    assert list(t.columns) == list(COLUMNS) and list(t["point"]) == ACCOUNTS_ROWS, (list(t.columns), list(t["point"]))
    cells = [c for c in t.values.ravel()]
    assert not any(d in c for c in cells for d in DASH), "dash in the table"
    assert all("simulat" in c for c in t["tested, in simulation only"]), "a tested cell does not say simulation"
    row = {r["point"]: r for _, r in t.iterrows()}
    printed = {"Water age": [f"{n['age_cov_max_Net3']:.3f}", f"{n['age_cov_max_Net2']:.3f}", f"{n['age_cov_max_ky4']:.3f}",
                             f"{n['wall_share_mae_Net3']:.3f}",
                             *(f"{n[f'age_lb_max_{k}']} of {k}'s {n[f'age_n_junctions_{k}']}" for k in ("Net3", "Net2", "ky4"))],
               "Organics and chlorine demand": [f"{n['toc_readings_ratio_Net3']:.3f}", f"{n['toc_first_storm_ratio_Net3']:.3f}",
                                                f"{n['audit_n_within']} of {n['audit_n_cells']}",
                                                f"{abs(n['audit_oldest_bias_Net3']):.3f}"],
               "Chlorine type": [f"{n['ca_found_straddle8_Net3']} of {int(n['ca_low_straddle8_Net3'])}",
                                 f"{n['ca_cov90_random_Net3']:.3f}", f"{n['ca_cov90_random_Net2']:.3f}",
                                 f"fitted to 15 samples, put {n['ca_free_floor_mass_Net3']:.2f} of their posterior",
                                 f"at 8 random samples they raised {int(n['ca_false_alarms_free_Net3'])} false alarms",
                                 f"{n['ca_cov90_straddle8_Net3']:.3f}"],
               "Temperature": [f"{n['temp_m_julsep_rmse_Net3']:.3f}", f"{n['temp_b0_julsep_rmse_Net3']:.3f}",
                               f"{n['temp_m2_ratio_min']:.3f} to {n['temp_m2_ratio_max']:.3f}",
                               f"{n['temp_m2_n_stop']} of the 4"]}
    assert (n["age_n_junctions_Net3"], n["age_n_junctions_Net2"], n["age_n_junctions_ky4"]) == (92, 35, 959)
    for point, strs in printed.items():
        text = " ".join(row[point][c] for c in COLUMNS[1:])
        for s in strs:
            assert s in text, (point, s)
    assert not n["temp_adopted"] and row["Temperature"][COLUMNS[3]].startswith("Tested, not adopted")
    assert not n["toc_adopted"] and row["Organics and chlorine demand"][COLUMNS[3]].startswith("Tested, not adopted")
    assert n["toc_adopted"] == (n["toc_adopted_Net3"] and n["toc_adopted_Net2"]), "TOC adoption follows task 11's own flag"
    # the assumptions behind the temperature and organics numbers are named where the numbers are
    assert "assumed temperature schedule" in row["Temperature"][COLUMNS[2]], "temperature assumptions"
    assert "assumed TOC schedule" in row["Organics and chlorine demand"][COLUMNS[2]], "TOC assumptions"
    assert n["temp_m_stop_Net3"] and n["temp_m_stop_Net2"], "task 10's stop rule is what the table reports"
    assert ("experimental" in row["Chlorine type"][COLUMNS[1]]) == n["ca_experimental"]
    lines = {"free": report_lines(False, 0.2), "chloramine": report_lines(True, 0.5, 8.0, 4.5, 3),
             "custom": report_lines(False, 0.3)}
    for k, ls in lines.items():
        txt = " ".join(ls)
        assert not any(d in txt for d in DASH), k
        assert "Validated in simulation only" in txt and "absorbed into the decay rates" in txt, k
        assert "kb20" not in txt and "kb(T)" not in txt and "temperature this month" not in txt.lower(), k
    assert "Measured species: free chlorine" in lines["free"][0] and "not a California rule" in lines["free"][1]
    assert "Measured species: TOTAL chlorine" in lines["chloramine"][0] and "not a California rule" in lines["chloramine"][1]
    assert "nitrification watch: 3 junctions" in lines["chloramine"][0]
    assert "used only by the nitrification watch" in lines["chloramine"][2] and "No seasonal chloramine" in lines["chloramine"][2]
    assert "tested in simulation and not adopted" in lines["free"][2]
    assert "set by you" in lines["custom"][1] and "detectable residual" in lines["custom"][1]
    return {"rows": list(t["point"]), "numbers": {k: v for k, v in n.items() if not isinstance(v, (list, dict))},
            "pdf_lines_free": lines["free"]}


@check("task14")
def grid_edge_warning_text():
    """capabilities.grid_edge_warning, the app's grid-edge warning, on the real grids with made-up posteriors: silent
    when no edge of the bulk or wall axis holds more than half of the calibration (uniform weights; 0.49 on the edge);
    on the free-chlorine grid with all the weight on the bulk floor it names the rate (0.1 per day) and adds the
    simulations' floor note; on the chloramine grid with the weight on the highest wall rate it names that rate and has
    no floor note; no em or en dash."""
    import itertools
    from .capabilities import EDGE_WARN, grid_edge_warning
    from .simgp import GRIDS, grid_edge_mass
    out = {}
    for grid, axis, side in (("full", 0, "low"), ("chloramine", 1, "high")):
        params = list(itertools.product(*GRIDS[grid]))
        P = np.asarray(params, dtype=float)
        edge_v = P[:, axis].min() if side == "low" else P[:, axis].max()
        uniform = np.full(len(params), 1.0 / len(params))
        assert grid_edge_warning(grid_edge_mass(uniform, params), params, grid == "chloramine") is None, grid
        on = (P[:, axis] == edge_v).astype(float)
        under = EDGE_WARN - 0.01          # just under the trigger on that edge
        near = under * on / on.sum() + (1 - under) * (1 - on) / (1 - on).sum()
        assert abs(grid_edge_mass(near, params)[f"{'kb' if axis == 0 else 'kw'}_{side}"] - under) < 1e-12
        assert grid_edge_warning(grid_edge_mass(near, params), params, grid == "chloramine") is None, grid
        msg = grid_edge_warning(grid_edge_mass(on / on.sum(), params), params, grid == "chloramine")
        assert msg and msg.startswith("Grid edge: 100% of the calibration") and not any(d in msg for d in DASH), msg
        out[grid] = msg
    assert "lowest bulk decay rate the grid allows (0.1 per day)" in out["full"] and "floor of 0.10 per day" in out["full"]
    assert "highest wall decay rate the grid allows (1 m/day)" in out["chloramine"] and "floor" not in out["chloramine"]
    return out


APP_T14_SCRIPT = r'''
import io, json, sys, time, warnings
warnings.filterwarnings("ignore")
import streamlit.runtime.memory_media_file_storage as mm
PDF = {}
_orig = mm.MemoryMediaFileStorage.load_and_get_id
def _rec(self, path_or_data, mimetype, kind, filename=None):
    if mimetype == "application/pdf":
        PDF["last"] = path_or_data
    return _orig(self, path_or_data, mimetype, kind, filename)
mm.MemoryMediaFileStorage.load_and_get_id = _rec
try:      # reads the PDF's text; without pypdf the check still requires a PDF and skips only its text assertions
    from pypdf import PdfReader
except ImportError:
    PdfReader = None
from streamlit.testing.v1 import AppTest
app, inp = sys.argv[1], sys.argv[2]
T = 900

def state(at, t0):
    assert not at.exception, [e.value for e in at.exception]
    acc = [t.value for t in at.table if list(t.value.columns) == ["in the model", "tested, in simulation only",
                                                                  "not done, or tested and not adopted"]]
    metric = [m for m in at.metric if m.label.startswith("Junctions likely below")]
    calib = [m for m in at.metric if m.label.startswith("Calibrated decay")]
    raw = PDF.pop("last", None)
    pdf = PdfReader(io.BytesIO(raw)).pages[0].extract_text() if raw is not None and PdfReader is not None else None
    return {"seconds": round(time.time() - t0, 1),
            "flagged_label": metric[0].label if metric else None, "flagged": metric[0].value if metric else None,
            "calib_value": calib[0].value if calib else None,
            "subheaders": [s.value for s in at.subheader], "warnings": [w.value for w in at.warning],
            "errors": [e.value for e in at.error], "info": [i.value for i in at.info],
            "captions": [c.value for c in at.caption], "accounts_rows": list(acc[0].index) if acc else None,
            "pdf_made": raw is not None, "pdf": pdf,
            "reveal": [m.value for m in at.markdown if "True daily-minimum violations" in m.value]}

def fresh():
    PDF.clear()
    return AppTest.from_file(app, default_timeout=T).run()

def reveal(at):
    [t for t in at.toggle if t.label.startswith("Show the true daily-minimum map")][0].set_value(True)
    PDF.clear()
    return at.run()

def pick(at, label, prefix):
    if label == "selectbox":
        at.selectbox[0].set_value([o for o in at.selectbox[0].options if o.startswith(prefix)][0])
    else:
        [r for r in at.radio if r.label == label][0].set_value(prefix)

res = {}
t0 = time.time(); at = reveal(fresh()); res["default"] = state(at, t0)              # (1) free chlorine default demo
t0 = time.time(); pick(at, "Show", "03:00"); PDF.clear(); at.run(); res["hourly"] = state(at, t0)   # (7) hourly view
t0 = time.time(); at = fresh(); pick(at, "selectbox", "Net2"); PDF.clear(); at.run(); res["net2"] = state(at, t0)  # (4)
t0 = time.time(); at = fresh(); pick(at, "Disinfectant", "Chloramine (total chlorine)"); at.run()   # (2) chloramine, Net3
at = reveal(at); res["chloramine_net3"] = state(at, t0)
t0 = time.time(); at = fresh(); pick(at, "Disinfectant", "Chloramine (total chlorine)"); pick(at, "selectbox", "ky4")
PDF.clear(); at.run(); res["ky4_chloramine"] = state(at, t0)                         # (5) ky4: refused before any run
import residualmap.chloramine as CA                                                 # (8) EPANET-MSX unavailable
def _no_msx(*a, **k):
    raise RuntimeError("EPANET-MSX needs Homebrew's OpenMP runtime, which is missing (simulated by this check)")
_keep, CA.truth_in_subprocess = CA.truth_in_subprocess, _no_msx
try:
    t0 = time.time(); at = fresh(); pick(at, "selectbox", "Net2"); pick(at, "Disinfectant", "Chloramine (total chlorine)")
    at.run(); at = reveal(at); res["msx_fallback"] = state(at, t0)
finally:
    CA.truth_in_subprocess = _keep
t0 = time.time(); at = fresh(); pick(at, "Model", "Upload my .inp"); at.run()      # (6) the upload path, Net3.inp
at.file_uploader[0].set_value(("Net3.inp", open(inp, "rb").read(), "application/octet-stream")); PDF.clear(); at.run()
res["upload"] = state(at, t0)
print("APP_T14_JSON " + json.dumps({"pypdf": PdfReader is not None, "runs": res}))
'''


def _reveal_counts(txt: str) -> dict:
    m = re.search(r"\*\*(\d+) of (\d+)\*\* junctions.*found \*\*(\d+) of (\d+)\*\*.*with (\d+) false alarms", txt)
    out = dict(zip(("violations", "junctions", "found", "unsampled_violations", "false_alarms"), map(int, m.groups())))
    out["sample_mean_mgL"] = float(re.search(r"Mean of your samples: ([0-9.]+) mg/L", txt).group(1))
    return out


@check("task14", quick=False)
def app_task14_runs():
    """Headless AppTest by widget label of task 14's app (the plan's seven runs, with its two temperature-outlook runs,
    dropped with the temperature models, replaced by the ky4 chloramine refusal and the EPANET-MSX fallback):
    (1) the free chlorine default demo (Net3, scenario 0, 8 samples): 42 of 92 flagged; 41 violations, 36 of 36 found,
        1 false alarm, sample mean 0.53 mg/L; the 'What this model accounts for' table with its four rows; the PDF's
        chemistry lines (free chlorine, its threshold note, the seasons absorbed by refitting, 'Validated in
        simulation only') and no kb20 field;
    (2) the chloramine demo on Net3 (EPA's chemistry in EPANET-MSX): the 0.5 mg/L total chlorine label, the TOTAL
        chlorine banner and sample-table caption, a calibrated bulk rate from the chloramine grid, the PDF's chloramine
        lines, no reduced-truth warning;
    (4) Net2, free chlorine; (5) ky4 with the chloramine demo: refused with its message before any simulation;
    (6) the upload path with Net3.inp: a fresh grid built inside the app (finding 1 of task 14; the grid's cache file
        must be created by the run) and the same map as the example network; (7) the hourly view; (8) EPANET-MSX
        unavailable (simulated): the reduced chloramine truth, labelled.  An earlier upload of the same Net3.inp and its
        cached grids (from a crashed run, or from Devansh uploading the file) are moved aside first, so the upload run
        always builds its grid, and put back afterwards; the run's own copy and grid are deleted.  Every run that renders
    the page must produce a PDF; its text is read with pypdf, and without pypdf only the PDF-text assertions are skipped,
    which the run log records (pdf_text_checked), so the committed report does not depend on it."""
    import hashlib
    from .capabilities import report_lines
    from .chemistry import THRESHOLD_NOTE
    from .simgp import GRIDS, KB_GRID
    from .simulate import LIB
    inp = os.path.join(LIB, "Net3.inp")
    up_name = f"{hashlib.sha1(open(inp, 'rb').read()).hexdigest()[:8]}_Net3.inp"    # app.py's name for the upload
    up_dir = os.path.join(REPO, "outputs", "app_uploads")
    had_dir = os.path.isdir(up_dir)
    mine = lambda d: [f for f in (os.listdir(d) if os.path.isdir(d) else []) if f == up_name or f"_{up_name}" in f]
    aside = tempfile.mkdtemp(prefix="rm_t14_aside_")
    moved = [(d, f) for d in (up_dir, CACHE) for f in mine(d)]
    built = []
    try:
        for d, f in moved:
            os.makedirs(os.path.join(aside, os.path.basename(d)), exist_ok=True)
            shutil.move(os.path.join(d, f), os.path.join(aside, os.path.basename(d), f))
        try:
            r = subprocess.run([PY, "-c", APP_T14_SCRIPT, os.path.join(REPO, "app.py"), inp],
                               cwd=os.getcwd(), capture_output=True, text=True, env={**os.environ, "PYTHONPATH": REPO})
        finally:  # the upload test's copy of Net3.inp and its grid (outputs/app_uploads and outputs/cache, git-ignored)
            built = mine(CACHE)
            for d in (up_dir, CACHE):
                for f in mine(d):
                    os.remove(os.path.join(d, f))
            if not had_dir and os.path.isdir(up_dir) and not os.listdir(up_dir):
                os.rmdir(up_dir)
    finally:      # put back what was there before
        for d, f in moved:
            os.makedirs(d, exist_ok=True)
            shutil.move(os.path.join(aside, os.path.basename(d), f), os.path.join(d, f))
        shutil.rmtree(aside, ignore_errors=True)
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_T14_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-3000:]}")
    payload = json.loads(line[-1][len("APP_T14_JSON "):])
    res, pdf_text = payload["runs"], payload["pypdf"]
    if not pdf_text:
        print("      app_task14_runs: pypdf is not installed, so the PDF-text assertions were skipped "
              "(pip install pypdf; requirements.txt lists it)", flush=True)
    d = res["default"]
    got = {**_reveal_counts(d["reveal"][0]), "flagged": int(d["flagged"].split()[0])}
    assert got == APP_DEMO, got
    assert d["flagged_label"] == "Junctions likely below 0.2 mg/L (daily minimum)"
    for k in ("default", "hourly", "net2", "chloramine_net3", "msx_fallback", "upload"):
        assert res[k]["accounts_rows"] == ACCOUNTS_ROWS, (k, res[k]["accounts_rows"])
        assert "What this model accounts for" in res[k]["subheaders"], k
        assert res[k]["pdf_made"], k
        assert res[k]["pdf"] or not pdf_text, k
    if pdf_text:
        pdf = _norm_ws(d["pdf"])
        for ln in report_lines(False, 0.2):
            assert _norm_ws(ln) in pdf, ln
        assert "kb20" not in pdf and THRESHOLD_NOTE["free_chlorine"] in pdf
    assert res["hourly"]["flagged_label"] == "Junctions likely below 0.2 mg/L (03:00)", res["hourly"]["flagged_label"]
    assert res["net2"]["flagged_label"] == "Junctions likely below 0.2 mg/L (daily minimum)"
    c = res["chloramine_net3"]
    assert c["flagged_label"] == "Junctions likely below 0.5 mg/L (daily minimum)", c["flagged_label"]
    assert any(i.startswith("Chloramine mode") and "TOTAL chlorine" in i for i in c["info"]), c["info"]
    assert any("TOTAL chlorine reading" in x for x in c["captions"])
    kb = float(c["calib_value"].split()[0])
    assert kb in GRIDS["chloramine"][0] and kb not in KB_GRID, kb
    if pdf_text:
        cpdf = _norm_ws(c["pdf"])
        assert "Disinfectant: chloramine. Measured species: TOTAL chlorine" in cpdf and THRESHOLD_NOTE["chloramine"] in cpdf
        assert "Validated in simulation only" in cpdf and "kb20" not in cpdf
    assert not any(w.startswith("Reduced demo truth") for w in c["warnings"]), c["warnings"]
    k4 = res["ky4_chloramine"]
    assert any("offered on Net3 and Net2 only" in e for e in k4["errors"]) and k4["accounts_rows"] is None, k4["errors"]
    fb = res["msx_fallback"]
    assert any(w.startswith("Reduced demo truth") and "flatter" in w for w in fb["warnings"]), fb["warnings"]
    assert fb["reveal"] and fb["flagged_label"] == "Junctions likely below 0.5 mg/L (daily minimum)"
    u = res["upload"]
    assert f"grid24_full_{up_name}.pkl" in built, ("the upload run built no grid", built)
    assert u["flagged"] == d["flagged"] and u["calib_value"] == d["calib_value"], (u["flagged"], d["flagged"])
    assert any("has not been tested against a simulated truth for this file" in x for x in u["captions"])
    out = {"default": got, "hourly_flagged": res["hourly"]["flagged"], "net2_flagged": res["net2"]["flagged"],
           "chloramine_net3": {"flagged": c["flagged"], **_reveal_counts(c["reveal"][0]), "calibrated": c["calib_value"]},
           "msx_fallback_net2": {"flagged": fb["flagged"], **_reveal_counts(fb["reveal"][0])},
           "upload_flagged": u["flagged"], "upload_grid_built_in_the_app": True, "ky4_chloramine_refused": True,
           "grid_edge_warnings": {k: sum(w.startswith("Grid edge") for w in v["warnings"]) for k, v in res.items()},
           "_volatile": {"pdf_text_checked": pdf_text, **{f"{k}_seconds": v["seconds"] for k, v in res.items()}}}
    return out


@check("task14")
def numbers_trace_to_outputs():
    """python -m residualmap.numbers: every number in the iteration-4 parts of the README and the journal is found in a
    committed output file or classified (time, date, citation, identifier, run cost, derived count, value written in
    the code, or an allowlisted literature value); none is unexplained; every task-14 headline number in the registry
    is its file's value at its key, printed in its section with the words around it; two wrong entries (another
    network's value, a coarser format) fail; every source of the app's table resolves.  The scan's power (how often a
    one-unit change in the last printed digit would be caught) is recorded, not barred."""
    from . import numbers as N
    res = N.scan()
    N.check_registry(res)
    bad = [m for ok, e, m in res.registry if not ok]
    assert res.sources_ok, res.sources_error
    assert not res.unexplained, res.unexplained[:10]
    assert not bad and len(res.registry) > 0, bad[:10]
    # the registry can fail: Net2's value at the place where the README prints Net3's, and a value one digit off
    keep = N.REGISTRY[:]
    try:
        N.REGISTRY[:] = [("README", N.ONE, "age at {} (Net3)", "age_cov_max_Net2", ".3f"),
                         ("README", N.ONE, "age at {} (Net3)", "age_cov_max_Net3", ".2f")]
        neg = N.Result()
        N.check_registry(neg)
    finally:
        N.REGISTRY[:] = keep
    assert [ok for ok, _, _ in neg.registry] == [False, False], neg.registry
    return {"unexplained": 0, "registry_ok": len(res.registry), "registry_negative_controls_fail": 2,
            "classified": {k: len(v) for k, v in sorted(res.classified.items())},
            "traced": res.traced, "n_tokens": res.n_tokens,
            "power_caught_of_tested": {k: list(v) for k, v in res.power.items()}}


def _it4_texts() -> dict:
    """The iteration-4 texts task 14 holds to 'no em or en dashes': the journal from 'Iteration 4, chemistry' on, the
    README's iteration-4 sections (numbers.sections()), the CHANGELOG's iteration-4 entries and capabilities.py."""
    from . import numbers as N
    out = {f"{s.doc}: {s.title[:50]}": s.text for s in N.sections()}
    ch = open(os.path.join(REPO, "CHANGELOG.md")).read()
    out["CHANGELOG iteration 4"] = ch[:ch.index("## 2026-09-21")]
    out["residualmap/capabilities.py"] = open(os.path.join(REPO, "residualmap", "capabilities.py")).read()
    for f in IT5_FILES:        # iteration 5 (task 15 on): its new modules and files, whole
        out[f] = open(os.path.join(REPO, f), encoding="utf-8").read()
    return out


@check("task14")
def no_dashes_in_iteration4_text():
    """No em or en dash in any iteration-4 text (journal, README sections, CHANGELOG entries, capabilities.py), and task
    14 adds none to the files whose older parts have them (pilot protocol, app.py, feature dictionary, README).  Since
    iteration 5 (task 15) the same holds for the iteration-5 journal sections and CHANGELOG entries (both inside the
    spans above) and for its new files, whole (IT5_FILES)."""
    bad = {k: [t.count(d) for d in DASH] for k, t in _it4_texts().items() if any(d in t for d in DASH)}
    assert not bad, bad
    now = {}
    for f, (em, en) in DASHES_BEFORE_TASK14.items():
        t = open(os.path.join(REPO, f)).read()
        now[f] = (t.count(DASH[0]), t.count(DASH[1]))
        assert now[f][0] <= em and now[f][1] <= en, (f, now[f], (em, en))
    return {"iteration4_texts_checked": len(_it4_texts()), "older_files": now}


@check("task14")
def no_listed_names_in_tracked_files():
    """No tracked text file names an organisation on a private list kept outside the repository (the CSV named by
    RESIDUALMAP_NAMES_CSV, column 'org' and, when present, 'system_name'; rows with a 'type' column other than
    'utility' are left out).  A name matches as a whole sequence of words, case and punctuation ignored, with any
    parenthetical dropped; its leading words before the first generic word ('water', 'district', 'county', ...) are
    matched too when there are at least two of them.  The one allowed match is the public provenance of EPA's Net3
    example in experiment.NET_DESC, which reaches committed figure titles and is left as it is.  Without the variable
    the check is skipped: the list is never committed, and only counts are recorded, never a name."""
    import csv
    path = os.environ.get(NAMES_ENV)
    if not path or not os.path.exists(path):
        raise SkipCheck(f"{NAMES_ENV} not set")
    generic = set("water district county community services service mutual company co inc city town village of the and "
                  "public utility utilities municipal irrigation authority department dept system systems association "
                  "assn corp corporation cwd csd mwd pud wd ccsd ud mwc wsd sd improvement assessment area zone no rural "
                  "cooperative coop board works ca california llc".split())
    words = lambda s: re.findall(r"[a-z0-9]+", s.lower())
    names = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("type") not in (None, "", "utility"):
                continue
            for k in ("org", "system_name"):
                w = words(re.sub(r"\(.*?\)", " ", row.get(k) or ""))
                if w:
                    names.add(tuple(w))
                    stem = []
                    for x in w:
                        if x in generic:
                            break
                        stem.append(x)
                    if len(stem) >= 2:
                        names.add(tuple(stem))
    lengths: dict[str, set] = {}
    for n in names:
        lengths.setdefault(n[0], set()).add(len(n))
    files = subprocess.run(["git", "-C", REPO, "ls-files"], capture_output=True, text=True).stdout.split()
    files += [f for f in subprocess.run(["git", "-C", REPO, "ls-files", "--others", "--exclude-standard"],
                                        capture_output=True, text=True).stdout.split() if f.startswith("residualmap/")]
    allowed, unexpected, n_files = 0, [], 0
    for f in files:
        if f.endswith((".png", ".pkl", ".pdf", ".bin", ".jpg")):
            continue
        try:
            lines = open(os.path.join(REPO, f), encoding="utf-8").read().splitlines()
        except (OSError, UnicodeDecodeError):
            continue
        n_files += 1
        for i, ln in enumerate(lines):
            w = words(ln)
            hit = any(tuple(w[a:a + k]) in names for a in range(len(w)) if w[a] in lengths for k in lengths[w[a]])
            if hit:
                provenance = f == "residualmap/experiment.py" and any("NET_DESC" in x for x in lines[max(0, i - 3):i + 1])
                if provenance:
                    allowed += 1
                else:
                    unexpected.append(f"{f}:{i + 1}")
    assert not unexpected, unexpected
    return {"names_checked": len(names), "text_files_scanned": n_files, "allowed_provenance_lines": allowed,
            "unexpected": 0}


# ----------------------------------------------------------------------------- task 15: read any model file safely
INTAKE_RECORD = os.path.join(REPO, "outputs", "intake", "task15_import.json")


def _intake_record(part: str) -> dict:
    """One part of the committed task-15 record (python -m residualmap.fixtures --record)."""
    with open(INTAKE_RECORD) as fh:
        return json.load(fh)[part]


def _same_as_record(got: dict, part: str) -> dict:
    want = _intake_record(part)
    got = json.loads(json.dumps(got))
    diff = _json_diff(got, want)
    assert not diff, (part, diff[:10])
    return want


def _count_files(root: str) -> int:
    return sum(len(fs) for _, _, fs in os.walk(root))     # does not follow the scratch cwd's cache link


@check("task15")
def sniff_identifies_realistic_files():
    """importers.sniff on the 36 fixtures of fixtures.SNIFF_TABLE (model text in every encoding and wrapping the review
    found, drawings, vendor model files, GIS layers, tables, a PDF, an IFC file and an EPANET results file): every one
    comes back as its expected kind, none 'unknown'; the file name decided something only for the two files whose names
    say they are WaterCAD or WaterGEMS databases, and then the detail says '(from the file name)'; the three .bak files
    are reported as backup copies of what they contain.  Two OLE headers with no directory (TIE_BREAKS) cannot be typed
    by content: the extension breaks the tie, with low confidence, and the result says so.  The review's cases
    (fixtures.REVIEW_SNIFF: a SWMM model and an EPANET-MSX file, which share EPANET headers; CR-only and CR-CR-LF line
    ends; node and link sections after the first 64 KB; zips made on a Mac; a city-size shapefile; malformed and deeply
    nested JSON; a DOCTYPE with an internal subset; from the final verification, singular core headers such as [PIPE]
    and a form feed before a header, both of which WNTR reads, and a small model padded with NUL bytes after [END])
    each come back as their kind, by content, without an
    exception; the inputs that took seconds to hours before the fixes (a run of XML comments or processing
    instructions, a GeoPackage whose gpkg_contents is a slow view) each finish within fixtures.TIME_BOUND_S, and the
    two quote-heavy CSV texts (the 263-byte one from the final verification) within fixtures.CSV_TIME_BOUND_S, a bound
    the old code misses on any machine.  The whole result equals the committed outputs/intake/task15_import.json (part
    'sniff'; timings are recorded as within the bound or not, never as seconds)."""
    from . import fixtures as FX
    from .importers import KINDS
    with tempfile.TemporaryDirectory(prefix="rm_t15_") as d:
        paths = FX.build_all(d, [r[0] for r in FX.SNIFF_TABLE] + [t[0] for t in FX.TIE_BREAKS]
                             + [r[0] for r in FX.REVIEW_SNIFF])
        got = FX.sniff_part(paths)
    wrong = [(r["fixture"], r["kind"], r["by_extension"], r["backup"]) for r in got["rows"] if not r["ok"]]
    assert not wrong, wrong
    assert (got["fixtures"], got["identified"], got["unknown"]) == (36, 36, 0), got
    assert got["decided_with_the_name"] == got["name_use_stated"] == 2, got
    assert got["backups"] == 3 and got["tie_breaks_ok"] == len(FX.TIE_BREAKS) == 2, got
    bad = [(r["fixture"], r["kind"], r["raised"]) for r in got["review_rows"] if not r["ok"]]
    assert not bad and got["review_ok"] == got["review_cases"] == len(FX.REVIEW_SNIFF), bad
    assert got["review_raised"] == 0 and got["review_timed_within_bound"] == got["review_timed"] > 0, got
    assert not any(d in t for v in KINDS.values() for t in v[1:] for d in DASH), "dash in the advice"
    _same_as_record(got, "sniff")
    return {k: got[k] for k in ("fixtures", "identified", "unknown", "decided_with_the_name", "backups", "kinds",
                                "tie_breaks_ok", "review_cases", "review_ok", "review_timed",
                                "review_timed_within_bound")}


@check("task15")
def inp_encodings_normalized():
    """prepare_inp on an .inp with a UTF-8 byte-order mark, one in cp1252 with an accented title and one in UTF-16
    (each of which WNTR cannot read as it is: the control records ENSyntaxError for the first and UnicodeDecodeError for
    the others): each gives a UTF-8 copy, named by its content, that WNTR loads with
    Net3's nodes and links, and lists the decoding as a fix.  Net3.inp as WNTR ships it (CRLF), a CRLF copy and an LF
    copy come back as the same path with the same bytes, and nothing is written.  The review's cases: CR-only and
    CR-CR-LF line ends, a [TITLE] past 64 KB, a large [PATTERNS] block first and a 70 KB comment header also pass
    through untouched, and WNTR reads each as Net3, as do Net3 with singular core headers and Net3 with a form feed
    before each (the final verification), and Net1 padded with NUL bytes after its [END] passes through and reads as
    Net1; a cp1252 file with a byte cp1252 leaves undefined keeps its euro sign and
    curly quotes (only that byte is read as latin-1).  A clean model uploaded as bytes under a 250-character name (the
    final verification: it raised OSError and left its temporary copy) is prepared under a shortened name, with one
    file left and no exception, and a name that fits is kept as it is.  Equals the committed record (part
    'encodings')."""
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t15_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.ENCODING_CASES + FX.PASS_THROUGH_CASES
                             + FX.REVIEW_PASS_THROUGH + FX.REVIEW_ENCODING)
        out = os.path.join(d, "prepared")
        got = FX.encodings_part(paths, out)
        written = sorted(os.listdir(out))
    assert got["normalized"] == got["of"] == 3 and got["passed_through"] == got["pass_through_of"] == 3, got
    assert all(r["wntr_on_the_original"].startswith("raises") for r in got["rows"]), got["rows"]   # the control
    assert got["review_passed_through"] == got["review_pass_through_of"] == len(FX.REVIEW_PASS_THROUGH), got
    gaps = got["review_cp1252_gaps"]
    assert all(g["ok"] and g["euro_sign_and_quotes_kept"] and g["wntr_loads_same_network"] for g in gaps), gaps
    assert written == sorted(r["prepared"] for r in got["rows"] + gaps), written
    long_name, fits = got["review_upload_names"]
    assert long_name["raised"] is None and long_name["ok"] and not long_name["name_kept"], long_name
    assert long_name["files_left"] == [long_name["prepared"]] and long_name["prepared_bytes"] <= 255, long_name
    assert fits["ok"] and fits["name_kept"] and fits["files_left"] == [fits["prepared"]], fits
    _same_as_record(got, "encodings")
    return {"normalized": got["normalized"], "passed_through": got["passed_through"],
            "decoded_as": {r["fixture"]: r["decoded_as"] for r in got["rows"] + gaps},
            "review_passed_through": got["review_passed_through"],
            "long_upload_name": {k: long_name[k] for k in ("name_characters", "prepared_bytes", "files_left")}}


@check("task15")
def inp_duplicate_ids_block():
    """A junction listed twice, a junction that shares a tank's ID, a pipe listed twice and a 32-character ID each block
    in prepare_inp with 'export with IDs, not labels', and WNTR is never called (its model class and .inp reader are
    traps during the four calls).  The control, WNTR on the same files with nothing in front: it keeps 92 of 93
    junction lines and 117 of 118 pipe lines without a word, loads the junction-tank clash and then the simulator
    raises AttributeError, and stops on the long ID with AssertionError.  A 31-character ID runs.  The review's cases
    block the same way, before WNTR: a 31-character ID with two accented letters saved in cp1252 (33 bytes in the
    UTF-8 copy; without the scan WNTR loads it and the simulator fails) and a 33-character pattern ID (without the
    scan WNTR raises AssertionError).  A model blocked after the scan (a byte-order mark, so a copy is written, and a
    pipe to a missing node, so WNTR stops) leaves nothing in the folder, and its message names no folder.  From the
    final verification: a duplicated pipe under a singular [PIPE] header, which WNTR reads as [PIPES] and keeps the
    later line of without a word, blocks the same way; and Net3 with a junction no pipe connects, which WNTR loads and
    EPANET cannot run, is blocked by EPANET's test run of the first hour (given as a path and as bytes, leaving
    nothing either way), naming the lone node.  Equals the committed record (part 'duplicate_ids')."""
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t15_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.DUPLICATE_CASES + ("net3_long_id_31.inp",)
                             + FX.REVIEW_ID_BLOCKS + (FX.REVIEW_WNTR_BLOCK, FX.REVIEW_TEST_RUN_BLOCK))
        out = os.path.join(d, "prepared")
        got = FX.ids_part(paths, out)
        written = sorted(os.listdir(out)) if os.path.isdir(out) else []
    assert got["blocked"] == got["of"] == 4 and got["wntr_calls_while_blocking"] == 0, got
    assert got["id_of_31_characters_runs"] and len(written) == 0, (got, written)    # the 31-character file is clean
    assert got["review_blocked"] == got["review_of"] == len(FX.REVIEW_ID_BLOCKS), got
    assert [r["without_the_scan"] for r in got["review_rows"]] == ["loads, and the simulator raises EpanetException",
                                                                  "raises AssertionError",
                                                                  "loads, and the simulator runs"], got["review_rows"]
    t = got["blocked_by_the_test_run"]
    assert t["without_the_test_run"] == "loads, and the simulator raises EpanetException", t
    for how in ("path", "bytes"):
        assert all(t[how][k] for k in ("blocked", "blocked_by_test_run", "names_the_lone_node", "wntr_read_it",
                                       "nothing_left")), (how, t[how])
    w = got["blocked_by_wntr_leaves_nothing"]
    assert w["blocked_by_wntr"] and w["fixes_it_needed"] and w["prepared"] is None and w["nothing_left"], w
    assert w["no_folder_in_message"], w
    _same_as_record(got, "duplicate_ids")
    return {"blocked": got["blocked"], "wntr_calls_while_blocking": 0,
            "without_the_scan": {c["fixture"]: c["wntr"] for c in got["without_the_scan"]},
            "review_blocked": got["review_blocked"],
            "review_without_the_scan": {r["fixture"]: r["without_the_scan"] for r in got["review_rows"]},
            "blocked_by_wntr_leaves_nothing": w["nothing_left"],
            "blocked_by_the_test_run": {how: t[how]["blocked_by_test_run"] and t[how]["nothing_left"]
                                        for how in ("path", "bytes")}}


@check("task15")
def zip_uploads_are_safe():
    """Zips are read in memory only: a member path that climbs out of the zip (../../evil.inp) or starts at the root
    (/tmp/evil.inp) is refused at the listing, and a model that would unpack to 21 MB at more than 100 to 1 is refused
    when it is about to be read; no member is opened for any of the three (zipfile's open is counted).  A zip with two
    models is blocked ('send one at a time'); the emailed zip's model is read in memory and written once, as a
    prepared .inp; a zipped shapefile whose padded DBF packs far more than 100 to 1 but unpacks to under 10 MB passes;
    a lowered total cap refuses.  The review's cases: a model and a shapefile zipped on a Mac (Finder's __MACOSX/ ._
    copies and a .DS_Store are left out, so the zip holds one model or one layer), a city-size shapefile whose DBF is
    over 10 MB and packs past 100 to 1 (only listed, so it passes) and the model zipped with that layer (only the
    model is read, so it runs).  From the final verification: a zip whose end record puts every member at a negative
    offset lists fine, and reading its model made zipfile raise ValueError (the control), which reached the app as a
    traceback; it is refused now, with a plain message and nothing written.  The file counts of the checks' scratch
    directory and of outputs/ are the same before and after, and only the prepared models are written.  Equals the
    committed record (part 'zips')."""
    from . import fixtures as FX
    n_scratch, n_out = _count_files(os.getcwd()), _count_files(os.path.join(REPO, "outputs"))
    with tempfile.TemporaryDirectory(prefix="rm_t15_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.REFUSED_ZIPS + ("zip_two_models.zip", "emailed_model.zip",
                                                                      "shp_padded_dbf.zip") + FX.REVIEW_ZIPS
                             + FX.REVIEW_BROKEN_ZIPS)
        out = os.path.join(d, "prepared")
        got = FX.zips_part(paths, out)
        written = sorted(os.listdir(out))
    assert (_count_files(os.getcwd()), _count_files(os.path.join(REPO, "outputs"))) == (n_scratch, n_out)
    assert got["refused"] == got["of"] == 3 and got["zip_member_opens_while_refusing"] == 0, got
    assert got["two_models_blocked"] and got["emailed_model"]["ok"], got
    models = [r["prepared"] for r in got["review_rows"] if r["kind"] == "epanet_inp_zip"]
    assert written == sorted([got["emailed_model"]["prepared"]] + models), written
    assert got["padded_dbf"]["kind"] == "shapefile_zip" and got["padded_dbf"]["dbf_ratio"] > 100, got["padded_dbf"]
    assert got["total_cap_refuses_when_lowered"]
    assert got["review_ok"] == got["review_of"] == len(FX.REVIEW_ZIPS), got["review_rows"]
    city = [r for r in got["review_rows"] if r["fixture"] == "city_mains_shp.zip"][0]
    assert city["dbf_mb"] > 10 and city["dbf_ratio"] > 100, city            # what the first rule refused
    assert all(b["kind"] == "epanet_inp_zip" and b["zipfile_alone"] == "raises ValueError" and b["refused"]
               and b["nothing_written"] for b in got["review_broken"]), got["review_broken"]
    _same_as_record(got, "zips")
    return {"refused": got["refused"], "reasons": {r["fixture"]: r["reason"] for r in got["rows"]},
            "padded_dbf_ratio": got["padded_dbf"]["dbf_ratio"], "files_written": len(written),
            "review_ok": got["review_ok"], "city_dbf": {k: city[k] for k in ("records", "dbf_mb", "dbf_ratio")},
            "broken_zips_refused": sum(b["refused"] for b in got["review_broken"])}


APP_T15_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
app, bak, gpkg = sys.argv[1], sys.argv[2], sys.argv[3]

def upload(name, path):
    at = AppTest.from_file(app, default_timeout=900).run()
    [r for r in at.radio if r.label == "Model"][0].set_value("Upload my .inp")
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    up = at.file_uploader[0]
    info = {"allowed_type": up.allowed_type, "label": up.label, "help": up.help}
    up.set_value((name, open(path, "rb").read(), "application/octet-stream"))
    at.run()
    assert not at.exception, (name, [e.value for e in at.exception])
    metric = [m for m in at.metric if m.label.startswith("Junctions likely below")]
    info.update({"flagged": metric[0].value if metric else None, "flagged_label": metric[0].label if metric else None,
                 "expanders": [e.label for e in at.expander],
                 "readiness": [m.value for m in at.markdown if m.value.startswith("**Model readiness")],
                 "warnings": [w.value for w in at.warning], "errors": [e.value for e in at.error],
                 "subheaders": [s.value for s in at.subheader]})
    return info

print("APP_T15_JSON " + json.dumps({"bak": upload("net3_backup.inp.bak", bak), "gpkg": upload("mains.gpkg", gpkg)}))
"""


@check("task15", quick=False)
def uploader_accepts_any_type():
    """Headless AppTest of the app's uploader: it sets no file type (allowed_type is empty; the label 'EPANET .inp file'
    is unchanged, so at.file_uploader[0] still finds it, and its help says any file is identified by its content).
    net3_backup.inp.bak (Net3 under a backup name) runs the Net3 path: a readiness report that identifies an EPANET
    model by its content, says the name marks it as a backup copy, says EPANET ran its first hour (the test run, after
    the final verification) and, since task 15b, gives the physics checks' rows (reaction orders, run cost) where task
    15 listed the physics settings as not checked yet, and the default demo's map, 42 of 92 junctions flagged.  mains.gpkg (a GeoPackage) gets the plain advice that it
    is a GIS layer and what to send instead, and the page stops before any simulation; nothing is written for it
    (measured: no file named after it in outputs/app_uploads or outputs/cache).  The run's copy of the backup and its
    grid (outputs/app_uploads and outputs/cache, git-ignored) are
    deleted afterwards, and an earlier copy is moved aside first and put back.  AppTest's set_value may bypass a
    browser's type filter, so one real-browser upload of a .bak is recorded in the journal (task 15)."""
    import hashlib
    from . import fixtures as FX
    data = FX.FIXTURES["net3_backup.inp.bak"]()
    up_name = f"{hashlib.sha1(data).hexdigest()[:8]}_net3_backup.inp"       # importers.prepare_inp's name for it
    up_dir = os.path.join(REPO, "outputs", "app_uploads")
    had_dir = os.path.isdir(up_dir)
    mine = lambda d: [f for f in (os.listdir(d) if os.path.isdir(d) else []) if f == up_name or f"_{up_name}" in f
                      or "mains" in f]
    aside = tempfile.mkdtemp(prefix="rm_t15_aside_")
    moved = [(d, f) for d in (up_dir, CACHE) for f in mine(d)]
    built = []
    try:
        for d, f in moved:
            os.makedirs(os.path.join(aside, os.path.basename(d)), exist_ok=True)
            shutil.move(os.path.join(d, f), os.path.join(aside, os.path.basename(d), f))
        with tempfile.TemporaryDirectory(prefix="rm_t15_") as fx:
            paths = FX.build_all(fx, ["net3_backup.inp.bak", "mains.gpkg"])
            try:
                r = subprocess.run([PY, "-c", APP_T15_SCRIPT, os.path.join(REPO, "app.py"), paths["net3_backup.inp.bak"],
                                    paths["mains.gpkg"]], cwd=os.getcwd(), capture_output=True, text=True,
                                   env={**os.environ, "PYTHONPATH": REPO})
            finally:
                built = {d: mine(d) for d in (up_dir, CACHE)}
                for d in (up_dir, CACHE):
                    for f in mine(d):
                        os.remove(os.path.join(d, f))
                if not had_dir and os.path.isdir(up_dir) and not os.listdir(up_dir):
                    os.rmdir(up_dir)
    finally:
        for d, f in moved:
            os.makedirs(d, exist_ok=True)
            shutil.move(os.path.join(aside, os.path.basename(d), f), os.path.join(d, f))
        shutil.rmtree(aside, ignore_errors=True)
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_T15_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-3000:]}")
    res = json.loads(line[-1][len("APP_T15_JSON "):])
    b, g = res["bak"], res["gpkg"]
    for k in (b, g):
        assert k["allowed_type"] == [] and k["label"] == "EPANET .inp file", k
        assert "identified by its content" in (k["help"] or ""), k["help"]
    assert b["flagged"] == f"{APP_DEMO['flagged']} of {APP_DEMO['junctions']}", b["flagged"]
    assert "Model readiness report" in b["expanders"] and b["readiness"], b
    assert "the name marks it as a backup copy" in b["readiness"][0] and not b["errors"], b
    assert "| Reaction orders |" in b["readiness"][0] and "| Run cost |" in b["readiness"][0], b["readiness"][0]
    assert "not checked yet" not in b["readiness"][0], b["readiness"][0]          # task 15b checks them
    assert "EPANET runs its first hour" in b["readiness"][0] and "| Test run |" in b["readiness"][0], b["readiness"][0]
    assert built[up_dir] == [up_name] and f"grid24_full_{up_name}.pkl" in built[CACHE], built
    gpkg_written = sorted(f for d in (up_dir, CACHE) for f in built[d] if "mains" in f)        # measured, not assumed
    gpkg_stopped = g["flagged"] is None and not g["expanders"] and "What this model accounts for" not in g["subheaders"]
    assert gpkg_stopped and not gpkg_written, (g, gpkg_written)
    assert any("GIS layer (GeoPackage)" in w and "EPANET export" in w for w in g["warnings"]), g["warnings"]
    return {"allowed_type": b["allowed_type"], "bak_flagged": b["flagged"], "bak_upload_name": "<sha1>[:8]_net3_backup.inp",
            "gpkg_advice": g["warnings"][0][:160], "gpkg_stopped_before_any_run": gpkg_stopped,
            "gpkg_files_written": len(gpkg_written)}


def _pilot_inputs(d: str) -> dict:
    """Two different models both named zone.inp (Net1, and Net1 with one junction's base demand doubled), a four-month
    grab log at five taps and its tap map, all invented, for the pilot's real-data path."""
    from .simulate import LIB
    net1 = open(os.path.join(LIB, "Net1.inp"), "rb").read()
    text = net1.decode("utf-8")
    i = text.index("[JUNCTIONS]")
    j = text.index("\n 22", i)
    line = text[j + 1:text.index("\n", j + 1)]
    tok = line.split()
    other = text[:j + 1] + line.replace(tok[2], str(2 * float(tok[2])), 1) + text[text.index("\n", j + 1):]
    assert other != text
    paths = {}
    for sub, data in (("a", net1), ("b", other.encode("utf-8"))):
        os.makedirs(os.path.join(d, sub))
        paths[sub] = os.path.join(d, sub, "zone.inp")
        with open(paths[sub], "wb") as fh:
            fh.write(data)
    junctions = ["10", "11", "12", "21", "22"]
    rows = ["date,time,tap_id,free_chlorine_mgL"]
    for m in range(1, 5):
        for k, _ in enumerate(junctions):
            rows.append(f"2026-{m:02d}-{3 + 5 * k:02d},{8 + 2 * k:02d}:15,T{k + 1},{0.9 - 0.12 * k + 0.02 * m:.2f}")
    paths["log"] = os.path.join(d, "grab_log.csv")
    with open(paths["log"], "w") as fh:
        fh.write("\n".join(rows) + "\n")
    paths["taps"] = os.path.join(d, "tap_map.csv")
    with open(paths["taps"], "w") as fh:
        fh.write("tap_id,junction_id\n" + "".join(f"T{k + 1},{j}\n" for k, j in enumerate(junctions)))
    return paths


@check("task15", quick=False)
def pilot_inp_is_content_hashed_and_private():
    """The pilot's real-data path, run twice from a fresh working directory with two different models that are both
    named zone.inp: each is copied to the git-ignored outputs/private/inputs/<sha1>[:8]_zone.inp, so the two get two
    grids with different values (before task 15 the second run was served the first one's grid, which is keyed by the
    file name).  A real run's results go to outputs/private/pilot by default, with import_report_zone.json beside them,
    and outputs/pilot (tracked, the synthetic demonstration) is not created; --synthetic keeps outputs/pilot.  The
    repository's .gitignore covers outputs/private/, and git status is unchanged.  After the review: --inp with the name
    of a network bundled with WNTR (Net3) runs by that name as before task 15 (its report, no private copy, so the
    cached grid is used), and a zipped model leaves one private copy, not two.  After the final verification, the
    journal's claim for that run is measured, not assumed: 'pilot --inp Net3' on docs/example_grab_log.csv and
    docs/example_tap_map.csv writes validation_Net3.json, validation_Net3.csv and predictions_Net3.csv byte-identical
    to the committed files in outputs/pilot (whose all-taps mean is rmse 0.077 and coverage90 0.923).  Since task 15b
    the two zone.inp files (Net1, one C value on every pipe) are calibrated on the 'full_fixed_gamma' grid, under its
    own hashed cache name."""
    import pickle as _pickle
    from . import fixtures as FX
    from .pilot import PRIVATE_INPUTS, PRIVATE_OUT, SYNTHETIC_OUT, default_out, prepare_pilot_inp
    status0 = _git_status()
    with tempfile.TemporaryDirectory(prefix="rm_t15_pilot_") as d:
        pin = os.path.join(d, "inputs")
        bundled, brep = prepare_pilot_inp("Net3", inputs_dir=pin)
        assert bundled == "Net3" and brep.ok and brep.same_file and not os.path.exists(pin), (bundled, brep.blocking)
        z = FX.build("emailed_model.zip", d)
        zp, zrep = prepare_pilot_inp(z, inputs_dir=pin)
        assert zrep.ok and sorted(os.listdir(pin)) == [os.path.basename(zp)], os.listdir(pin)
        zipped_copies = len(os.listdir(pin))
    with tempfile.TemporaryDirectory(prefix="rm_t15_pilot_") as d:      # the bundled name, end to end
        os.makedirs(os.path.join(d, "outputs"))
        os.symlink(CACHE, os.path.join(d, "outputs", "cache"))          # Net3's committed grid is read, not rebuilt
        r = subprocess.run([PY, "-m", "residualmap.pilot", "--inp", "Net3",
                            "--samples", os.path.join(REPO, "docs", "example_grab_log.csv"),
                            "--taps", os.path.join(REPO, "docs", "example_tap_map.csv")],
                           cwd=d, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": REPO})
        assert r.returncode == 0, r.stderr[-2000:]
        bundled_same = {}
        for f in ("validation_Net3.json", "validation_Net3.csv", "predictions_Net3.csv"):
            with open(os.path.join(d, PRIVATE_OUT, f), "rb") as a, open(os.path.join(REPO, "outputs", "pilot", f),
                                                                         "rb") as b:
                bundled_same[f] = a.read() == b.read()
    assert all(bundled_same.values()), bundled_same
    assert default_out(True) == SYNTHETIC_OUT == os.path.join("outputs", "pilot")
    assert default_out(False) == PRIVATE_OUT == os.path.join("outputs", "private", "pilot")
    for f in ("outputs/private/pilot/x.csv", "outputs/private/inputs/x.inp"):
        ign = subprocess.run(["git", "-C", REPO, "check-ignore", "-q", f], capture_output=True)
        assert ign.returncode == 0, f"{f} is not git-ignored"
    with tempfile.TemporaryDirectory(prefix="rm_t15_pilot_") as d:
        p = _pilot_inputs(d)
        runs = {}
        for sub in ("a", "b"):
            r = subprocess.run([PY, "-m", "residualmap.pilot", "--inp", p[sub], "--samples", p["log"], "--taps", p["taps"]],
                               cwd=d, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": REPO})
            assert r.returncode == 0, (sub, r.stderr[-2000:])
            runs[sub] = sorted(os.listdir(os.path.join(d, PRIVATE_INPUTS)))
        inputs = runs["b"]
        grids = sorted(f for f in os.listdir(os.path.join(d, "outputs", "cache")) if f.startswith("grid24_full_"))
        assert len(inputs) == 2 and all(f.endswith("_zone.inp") for f in inputs) and len(runs["a"]) == 1, runs
        from .simgp import grid_cache_path      # task 15b: Net1 has one C value, so its grid holds gamma at 0
        want = sorted(os.path.basename(grid_cache_path(types.SimpleNamespace(wn_name=f, source_dose=1.2), "",
                                                       "full_fixed_gamma")) for f in inputs)
        assert grids == want, (grids, want)
        Z = [_pickle.load(open(os.path.join(d, "outputs", "cache", g), "rb"))[1] for g in grids]
        grid_diff = float(np.abs(Z[0] - Z[1]).max())
        assert grid_diff > 0, grid_diff
        out = sorted(os.listdir(os.path.join(d, PRIVATE_OUT)))
        assert {"import_report_zone.json", "predictions_zone.csv", "validation_zone.csv", "validation_zone.json"} <= set(out), out
        assert not os.path.exists(os.path.join(d, "outputs", "pilot")), "a real run wrote to outputs/pilot"
        rep = json.load(open(os.path.join(d, PRIVATE_OUT, "import_report_zone.json")))
        assert rep["ok"] and rep["prepared"] == inputs[[f.startswith(rep["sha1"][:8]) for f in inputs].index(True)], rep
    assert _git_status() == status0, "git status changed"
    return {"private_inputs": len(inputs), "grids": len(grids), "grid_values_differ": grid_diff > 0,
            "default_out_real": PRIVATE_OUT, "default_out_synthetic": SYNTHETIC_OUT, "results": out,
            "bundled_name_runs_by_name": bundled, "zipped_model_private_copies": zipped_copies,
            "bundled_name_results_byte_identical_to_outputs_pilot": sorted(bundled_same)}


STREAMLIT_CFG_SCRIPT = r"""
import json, sys
from streamlit import config
config._main_script_path = sys.argv[1]          # what 'streamlit run <app.py>' sets before it reads the config
config.get_config_options(force_reparse=True)
print("CFG_JSON " + json.dumps({"address": config.get_option("server.address"),
                                "stats": config.get_option("browser.gatherUsageStats")}))
"""


@check("task15")
def streamlit_config_private():
    """.streamlit/config.toml holds server.address = 'localhost' (Streamlit otherwise listens on every network
    interface) and browser.gatherUsageStats = false, and Streamlit's own config loader, set up the way 'streamlit run
    app.py' sets it up and started from another folder, reads both."""
    import tomllib
    path = os.path.join(REPO, ".streamlit", "config.toml")
    with open(path, "rb") as fh:
        cfg = tomllib.load(fh)
    assert cfg.get("server", {}).get("address") == "localhost", cfg
    assert cfg.get("browser", {}).get("gatherUsageStats") is False, cfg
    with tempfile.TemporaryDirectory(prefix="rm_t15_cfg_") as d:
        r = subprocess.run([PY, "-c", STREAMLIT_CFG_SCRIPT, os.path.join(REPO, "app.py")], cwd=d, capture_output=True,
                           text=True, env={**os.environ, "PYTHONPATH": REPO})
    line = [x for x in r.stdout.splitlines() if x.startswith("CFG_JSON ")]
    assert r.returncode == 0 and line, r.stderr[-2000:]
    got = json.loads(line[-1][len("CFG_JSON "):])
    assert got == {"address": "localhost", "stats": False}, got
    return {"file": {"server.address": cfg["server"]["address"], "browser.gatherUsageStats": False},
            "streamlit_reads": got}


# ----------------------------------------------------------------------------- task 15b: an imported model's physics
INTAKE_RECORD_15B = os.path.join(REPO, "outputs", "intake", "task15b_physics.json")


def _record_15b(part: str) -> dict:
    """One part of the committed task-15b record (python -m residualmap.intake_exp record)."""
    with open(INTAKE_RECORD_15B) as fh:
        return json.load(fh)[part]


def _same_as_15b(got: dict, part: str) -> dict:
    want = _record_15b(part)
    diff = _json_diff(json.loads(json.dumps(got)), want)
    assert not diff, (part, diff[:10])
    return want


@contextlib.contextmanager
def _physics_fixtures():
    """The task-15b fixtures built into a temporary folder, and a folder for their prepared copies."""
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t15b_") as d:
        yield FX.build_all(os.path.join(d, "fx"), FX.PHYSICS_FIXTURES), os.path.join(d, "prepared")


@check("task15b", quick=False)
def inp_reaction_orders_fixed():
    """ky10 declares reaction orders of 0 (bulk, wall and tank): the prepared copy sets all three to 1 and lists the fix,
    and its grid member's mean residual differs between bulk rates 0.2 and 0.8 per day, while the original's is the
    same at both (probe d2 saw 0.5748; with order 0 the grid's rates have no effect).  Net6 gets the same fix and the
    size warning (3,323 junctions).  Equals the committed record (outputs/intake/task15b_physics.json, part
    'reaction_orders')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, _):
        got = E.part_reaction_orders(paths)
    k, n6 = got["ky10"], got["Net6"]
    assert set(k["orders_fixed"]) == {"bulk", "wall", "tank"} and k["fixes"], k
    assert k["original_same_at_both"] and abs(k["fixed_difference"]) > 1e-3, k
    assert set(n6["orders_fixed"]) == {"bulk", "wall", "tank"} and n6["large"] and n6["size_warning"], n6
    _same_as_15b(got, "reaction_orders")
    return {"ky10_mean_original": k["mean_original"], "ky10_mean_fixed": k["mean_fixed"],
            "net6_junctions": n6["junctions"], "net6_size_warning": n6["size_warning"]}


@check("task15b")
def inp_bulk_lines_cleared():
    """net3_bulk_lines.inp (BULK lines on 10 pipes, a TANK line on one tank, LIMITING POTENTIAL 0.5): the prepared copy
    clears them and its grid member at kb 0.40, kw 0.70, gamma 1 equals Net3's at every junction-hour (max difference
    0.0); the file as it came differs (the control).  Equals the committed record (part 'bulk_lines')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_bulk_lines(paths, out)
    assert got["ok"] and (got["bulk_lines"], got["tank_lines"], got["limiting_potential_fixed"]) == (10, 1, "0.5"), got
    assert got["max_diff_prepared_vs_net3"] == 0.0 and got["max_diff_original_vs_net3"] > 0, got
    _same_as_15b(got, "bulk_lines")
    return {k: got[k] for k in ("bulk_lines", "tank_lines", "junction_hours", "max_diff_prepared_vs_net3",
                                "max_diff_original_vs_net3")}


@check("task15b")
def inp_headloss_formula():
    """net3_dw.inp (Net3 written in Darcy-Weisbach by the inverse of the stated formula at Net3's own reference
    velocities) gives a prepared Hazen-Williams copy, and the gate's water-age numbers are reported; the pre-registered
    gate (median change at most 1 h, 95th percentile at most 4 h, over the scored day's junction-hours) passes.  The
    114 pipes a positive roughness can match recover their C within 0.5; the 3 C 199 pipes (a tiny roughness in the
    file) are reported, not scored.  net3_cm.inp (Chezy-Manning) is blocked.  The review's case, net3_dw_uniform.inp
    (one Darcy-Weisbach roughness on every pipe): the uniform test reads the file's own roughness, so it gets the
    'full_fixed_gamma' grid and the old-pipe factor is not identifiable, although its converted copy has many C
    values; net3_dw.inp's old-pipe row is tagged derived (its C values come from the conversion).  The second
    review's case, net3_dw_plastic.inp (0.0015 mm and 0.0024 mm on alternating pipes, about 1e-6 m apart): the uniform
    test compares roughness to 6 significant digits, so it finds 2 values and the 'full' grid (rounded to 6 decimals
    of a metre they were one value, and the row named a roughness no pipe has).  Equals the committed record (part
    'headloss')."""
    from . import intake_exp as E
    from .importers import DW_GATE_MEDIAN_H, DW_GATE_P95_H
    with _physics_fixtures() as (paths, out):
        got = E.part_headloss(paths, out)
    dw = got["dw"]
    assert dw["ok"] and dw["headloss_in_copy"] == "H-W", dw
    g = dw["gate"]
    assert g["passes"] and g["age_change_median_h"] <= DW_GATE_MEDIAN_H and g["age_change_p95_h"] <= DW_GATE_P95_H, g
    assert (dw["pipes"], dw["pipes_scored"], len(dw["pipes_reported_not_scored"])) == (117, 114, 3), dw
    assert dw["within_tolerance"] == 114 and dw["max_c_error"] <= 0.5, dw
    assert got["cm"]["blocked"], got["cm"]
    assert dw["old_pipe_row_source"] == "derived", dw            # the converted C values are derived, not the file's
    u = got["dw_uniform"]       # the review's case: one roughness on every pipe, many converted C values
    assert u["ok"] and u["file_roughness_values"] == 1 and u["converted_c_values"] > 1, u
    assert u["uniform_c"] and u["grid"] == "full_fixed_gamma" and u["row"].startswith("not identifiable"), u
    pl = got["dw_plastic"]      # the second review's case: two plastic-pipe roughness values
    assert pl["ok"] and pl["file_roughness_values"] == 2 and not pl["uniform_c"] and pl["grid"] == "full", pl
    assert pl["row"].startswith("learned from"), pl
    _same_as_15b(got, "headloss")
    return {"gate": g, "max_c_error": dw["max_c_error"], "within_0.5": dw["within_tolerance"],
            "unmatched_below_199": dw["max_below_unmatched"], "cm_blocked": True,
            "dw_uniform": {k: u[k] for k in ("converted_c_values", "uniform_c", "grid")},
            "dw_plastic": {k: pl[k] for k in ("file_roughness_values", "uniform_c", "grid")}}


@check("task15b")
def inp_clock_and_weekly():
    """Net2.inp: its 8 h start (START CLOCKTIME) is recorded and its in-file source at the dosed inflow is flagged as
    overridden; net3_clock6.inp gets 6 h; net3_weekly_pattern.inp (a 168-step pattern) gets the weekday warning.  Since
    the review, a start other than midnight also warns that every hour in the app and the pilot is a model hour (Net2
    and net3_clock6 warn, Net3 does not), and the start's HH:MM is rounded to the minute once (5:59:59 is 06:00, not
    05:60).  Since the second review the row and the warning go by that rounded minute too: a start within 30 s of
    midnight (net3_clock_235950.inp, net3_clock_000020.inp, end to end) is '0 h (00:00)' with no model-hours warning,
    and the warning quotes the rounded offset ('Add 6 h' for 5:59:59).  Equals the committed record (part
    'clock')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_clock(paths, out)
    n2 = got["Net2"]
    assert n2["start_clock_h"] == 8.0 and n2["source_flagged"] and n2["sources_overridden"] == ["1"], n2
    assert got["net3_clock6"]["start_clock_h"] == 6.0 and got["net3_clock6"]["ok"], got["net3_clock6"]
    w = got["net3_weekly_pattern"]
    assert w["weekly_patterns"] == ["1"] and w["weekday_warning"] and w["ok"], w
    assert n2["model_hours_warning"] and got["net3_clock6"]["model_hours_warning"], got
    assert got["Net3"]["start_clock_h"] == 0.0 and not got["Net3"]["model_hours_warning"], got["Net3"]
    from .importers import _clock
    edges = {"5:59:59": _clock(5 + 59 / 60 + 59 / 3600), "23:59:50": _clock(23 + 59 / 60 + 50 / 3600),
             "8:00": _clock(8.0), "6:30": _clock(6.5)}           # the review's 05:60 and 23:60
    assert edges == {"5:59:59": "06:00", "23:59:50": "00:00", "8:00": "08:00", "6:30": "06:30"}, edges
    for f, v in got["near_midnight"].items():      # the second review: no 'starts at 00:00, not midnight' warning
        assert v["ok"] and v["row"] == "0 h (00:00)" and not v["model_hours_warning"], (f, v)
    rd = got["rounding"]
    assert [rd[k]["model_hours_warning"] for k in ("23:59:50", "0:00:20", "5:59:59", "6:30")] == [False, False, True, True], rd
    assert (rd["5:59:59"]["row"], rd["5:59:59"]["warning_offset"]) == ("6 h (06:00)", "Add 6 h"), rd["5:59:59"]
    assert (rd["6:30"]["row"], rd["6:30"]["warning_offset"]) == ("6.5 h (06:30)", "Add 6.5 h"), rd["6:30"]
    _same_as_15b(got, "clock")
    return {"net2_start_h": n2["start_clock_h"], "net3_clock6_h": 6.0, "weekly_warning": True,
            "net2_multi_day_patterns": n2["multi_day_patterns"],
            "near_midnight": {k: (v["row"], v["model_hours_warning"]) for k, v in rd.items()}}


@check("task15b", quick=False)
def inp_prv_elevations():
    """ky10.inp (5 PRVs) is flagged for J-16 at elevation 0 and Net6.inp (2 PRVs) for its three; ky10_const_elev.inp
    (every junction at one elevation) is blocked with 'elevations needed at these valves', quotes the committed
    evidence, and lists its reaction fixes as fixes that would be applied, since no copy is kept; Net3.inp (no PRVs, 8
    junctions at 0) is not flagged.  outputs/intake/elevation_prv_ky10.json reproduces
    (ky10 with every junction at the mean elevation: the junctions whose water age moves by more than 2 h at some hour
    of the scored day, and the largest move).  Equals the committed record (part 'prv')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_prv(paths, out)
    assert got["ky10"]["flagged"] and got["ky10"]["zero_elevation_junctions"] == ["J-16"], got["ky10"]
    assert got["Net6"]["flagged"] and len(got["Net6"]["zero_elevation_junctions"]) == 3, got["Net6"]
    c = got["ky10_const_elev"]
    assert c["blocked"] and c["constant_elevation"] and c["quotes_evidence"], c
    assert c["fixes"] and c["fixes_heading"].startswith("**Fixes that would be applied**"), c
    assert not got["Net3"]["flagged"] and len(got["Net3"]["zero_elevation_junctions"]) == 8, got["Net3"]
    _same_as_15b(got, "prv")
    ev = E.elevation_prv()
    with open(E.ELEVATION_PRV) as fh:
        committed = json.load(fh)
    diff = _json_diff(json.loads(json.dumps(ev)), committed)
    assert not diff, diff[:10]
    return {"ky10_zero": got["ky10"]["zero_elevation_junctions"], "net6_zero": len(got["Net6"]["zero_elevation_junctions"]),
            "const_blocked": True, "net3_flagged": False,
            "elevation_evidence": {k: ev[k] for k in ("junctions", "junctions_moved", "max_change_h")}}


@check("task15b", quick=False)
def inp_coordinate_scale():
    """The ratio of each pipe's drawn length (with vertices) to its declared length: Net3, Net2 and Net6 are schematic
    (95th over 5th percentile above 1.5); ky4 and ky10 are to scale, in feet (median 3.281 per metre); net3_no_coords.inp
    has no coordinates.  Equals the committed record (part 'coordinates')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_coordinates(paths, out)
    for n in ("Net3", "Net2", "Net6"):
        assert got[n]["status"] == "schematic", (n, got[n])
    for n in ("ky4", "ky10"):
        assert (got[n]["status"], got[n]["units"]) == ("to scale", "feet"), (n, got[n])
    assert got["net3_no_coords"]["status"] == "no coordinates", got["net3_no_coords"]
    _same_as_15b(got, "coordinates")
    return {n: (v["status"], v.get("units"), v.get("spread")) for n, v in got.items()}


@check("task15b")
def inp_controls_and_structure():
    """net3_no_controls.inp is flagged 'pumps but no controls'; net3_speed_pattern_pumps.inp (no controls, a speed
    pattern on each pump) gets no such warning; net3_zero_demand.inp and net3_two_parts.inp are flagged;
    net3_no_source.inp is blocked; watergems_export_like.inp (connector pipes of 1 ft and 99 in tagged
    AddedForINPConversion) is reported.  The review's cases: net3_rule_premise.inp (no controls, one rule that reads
    pump 10 and sets pipe 20) flags both pumps, because a rule's premise reads a link and only its actions set one;
    Net6's one pipe of the connector size, with no tag, is counted and not reported.  Equals the committed record
    (part 'structure')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_structure(paths, out)
    assert got["net3_no_controls.inp"]["pump_warning"] and got["net3_no_controls.inp"]["ok"], got["net3_no_controls.inp"]
    sp = got["net3_speed_pattern_pumps.inp"]
    assert sp["ok"] and not sp["pump_warning"] and sp["pumps_without_controls"] == [], sp
    assert got["net3_zero_demand.inp"]["zero_total_demand"] and got["net3_two_parts.inp"]["parts"] == [97, 2], got
    assert not got["net3_no_source.inp"]["ok"] and got["net3_no_source.inp"]["blocking"], got["net3_no_source.inp"]
    wg = got["watergems_export_like.inp"]
    assert wg["ok"] and wg["connector_row"] and (wg["connector_pipes"], wg["connector_tags"]) == (2, 2), wg
    rp = got["net3_rule_premise.inp"]      # a rule that reads pump 10 and sets pipe 20: pump 10 is not controlled
    assert rp["ok"] and rp["pump_warning"] and rp["pumps_without_controls"] == ["10", "335"], rp
    n6 = got["Net6"]                       # one pipe of the connector size, no tag: counted, not reported
    assert (n6["connector_pipes"], n6["connector_tags"], n6["connector_row"]) == (1, 0, False), n6
    _same_as_15b(got, "structure")
    return {f: {k: v.get(k) for k in ("ok", "pump_warning", "zero_total_demand", "parts", "connector_pipes")}
            for f, v in got.items()}


@check("task15b", quick=False)
def uniform_c_fixes_gamma():
    """Net1.inp (C 100 on every pipe) is flagged: the old-pipe factor is not identifiable, and the file gets the
    'full_fixed_gamma' grid, whose 225 members each equal the matching gamma-0 member of Net1's 'full' grid bit for bit
    (both built in memory, never cached).  The committed networks keep more than one C value (Net3 5, Net2 2, ky4 2),
    so no default moves.  In chloramine mode (no fixed-gamma chloramine grid exists) Net1 gets the 1440-run chloramine
    grid, and its row says that grid still varies the factor, not that it is held at 0 (the review's case).  The second
    review's case, Net1 written as an export (net1_export_like.inp: a connector of 1 ft and 99 in with C 130 at its
    pump, tagged AddedForINPConversion): the uniform test leaves the recognised connector out, so the file still gets
    'full_fixed_gamma', one C value, and a row that says the connector was left out.  Equals the committed record
    (part 'uniform_c')."""
    from . import intake_exp as E
    got = E.part_uniform_c()
    n1 = got["Net1"]
    assert n1["uniform_c"] and n1["grid"] == "full_fixed_gamma" and n1["runs_per_dose"] == 225, n1
    assert got["fixed_gamma_members"] == got["members_bit_identical_to_full_gamma0"] == got["full_gamma0_members"] == 225
    assert got["gamma_values"] == [0.0] and got["same_order_as_full"], got
    assert all(v > 1 for v in got["committed_networks_c_values"].values()), got["committed_networks_c_values"]
    assert n1["row_says_held_at_0"], n1
    ca = got["Net1_chloramine"]           # no chloramine fixed-gamma grid: the row must not say the factor is held
    assert ca["uniform_c"] and ca["grid"] == "chloramine" and ca["runs_per_dose"] == 1440, ca
    assert not ca["row_says_held_at_0"] and ca["row_says_chloramine_grid_varies_it"], ca
    ex = got["Net1_export_like"]          # the second review: the export's connector is not a second C value
    assert ex["ok"] and ex["connector_row"] and (ex["connector_pipes"], ex["connector_tags"]) == (1, 1), ex
    assert ex["uniform_c"] and ex["c_values"] == 1 and ex["grid"] == "full_fixed_gamma", ex
    assert "connector pipe left out" in ex["row"], ex
    _same_as_15b(got, "uniform_c")
    return {"members": got["fixed_gamma_members"], "bit_identical": got["members_bit_identical_to_full_gamma0"],
            "committed_networks_c_values": got["committed_networks_c_values"],
            "net1_export_like": {k: ex[k] for k in ("c_values", "uniform_c", "grid")}}


@check("task15b", quick=False)
def run_cost_gate():
    """The run cost (task 15b, item 10).  Asserted, because it does not depend on timing: the estimate the product
    reports is the stated formula (one timed 7-day run times the grid's runs over importers.EFFECTIVE_SPEEDUP, from
    the unrounded run time, and the minutes shown are that estimate rounded); Net6 (3,323 junctions) gets the hard
    warning; a 10,000-junction count for chloramine on 11 workers with the free disk stood in at 3 GB is refused (1.38
    GB of grid plus about 1 GB of temporary files would leave under 1.5 GB), and at 6 GB accepted; one run's EPANET
    files on ky4, ky10 and Net6 are as the record says, and TEMP_MB_PER_JUNCTION is ky4's.  Equals the committed record
    (part 'run_cost').
    Not asserted: the pre-registered timing bar (Net3's estimate within a factor of 2 of the fresh grid time).  Its
    scored run missed (2.10, outputs/intake/run_cost_timing.json), and the ratio moves across the bar from run to run
    on an unchanged tree, so as an assertion it made the saved run pass or fail by chance.  The miss stands, recorded
    in the journal; this run's ratio, against the grid fresh_grid_equals_cached_grid built in this run (or, run alone,
    a fresh build here), goes to the run log only, and no passing run counts as a pass of that bar."""
    from . import importers as I
    from . import intake_exp as E
    from .simgp import build_grid_24h
    from .simulate import nominal_scenario
    got = E.part_run_cost()
    assert got["Net6"]["large"] and got["Net6"]["hard_warning"], got["Net6"]
    gate = got["mock"]["gate"]
    assert got["mock"]["workers"] == E.MOCK_WORKERS, got["mock"]
    assert not gate["3"]["accepted"] and gate["3"]["message_says_floor"] and gate["6"]["accepted"], gate
    t = got["temp_files_one_run"]
    assert round(I.TEMP_MB_PER_JUNCTION * t["ky4"]["junctions"], 1) == t["ky4"]["mb"], t
    _same_as_15b(got, "run_cost")
    cost = E.bundled("Net3")[1].run_cost
    estimate = cost["estimate_s"]                                  # the number the product reports, unrounded
    probe = I.run_cost(cost["junctions"], "free_chlorine", seconds_per_run=0.04)     # the formula, on a fixed run time
    assert probe["estimate_s"] == I.estimate_seconds(0.04, 675) == 0.04 * 675 / I.EFFECTIVE_SPEEDUP, probe
    assert abs(estimate * I.EFFECTIVE_SPEEDUP / cost["runs_per_dose"] - cost["seconds_per_run"]) <= 0.0005, cost
    assert cost["minutes"] == round(estimate / 60, 1), cost
    if FRESH_GRID_SECONDS:
        fresh, source = FRESH_GRID_SECONDS[-1], "fresh_grid_equals_cached_grid"
    else:
        t0 = time.time()
        build_grid_24h(nominal_scenario("Net3"), "full")
        fresh, source = time.time() - t0, "a fresh build in this check"
    ratio = fresh / estimate
    return {"estimate_is_the_stated_formula": True, "net6_hard_warning": True, "mock_10000_chloramine": gate,
            "mock_workers": got["mock"]["workers"], "temp_files_one_run": t,
            "timing_bar": "pre-registered, missed in its scored run (outputs/intake/run_cost_timing.json); not "
                          "asserted here, this run's ratio is in the run log only",
            "_volatile": {"net3_seconds_per_run": round(estimate * I.EFFECTIVE_SPEEDUP / cost["runs_per_dose"], 4),
                          "net3_estimate_s": round(estimate, 2), "net3_fresh_grid_s": round(fresh, 2),
                          "fresh_from": source, "ratio": round(ratio, 3),
                          "within_factor_2_this_run": 0.5 <= ratio <= 2.0}}


@check("task15b", quick=False)
def readiness_bundled_reproduces():
    """importers.prepare_inp on the six networks WNTR bundles (Net1, Net2, Net3, Net6, ky4, ky10), rerun, equals the
    committed outputs/intake/readiness_bundled.json (run times and the worker count are left out of it)."""
    from . import intake_exp as E
    got = E.readiness_bundled()
    with open(E.READINESS_BUNDLED) as fh:
        committed = json.load(fh)
    diff = _json_diff(json.loads(json.dumps(got)), committed)
    assert not diff, diff[:10]
    nets = got["networks"]
    return {"networks": sorted(nets), "ok": {n: v["ok"] for n, v in nets.items()},
            "fixes": {n: len(v["fixes"]) for n, v in nets.items()}, "warnings": {n: len(v["warnings"]) for n, v in nets.items()}}


@check("task15b")
def task15_leftovers_fixed():
    """Task 15's four leftovers, each reproduced on the task-15 code first: a SWMM model whose SWMM-only sections come
    after the first 64 KB (whose first 64 KB alone is still typed as EPANET, the control) is typed as SWMM with its
    advice, not blocked with WNTR's syntax error; a zip holding Net3 renamed zone_model.txt is a zip holding a model,
    prepared as Net3; a UTF-8 file with accented IDs and one stray cp1252 byte is read as UTF-8 with that byte as
    cp1252, keeps its IDs and warns of the mixed encoding (read as cp1252 throughout, as task 15 did, WNTR fails: the
    control); a file in cp1252 throughout whose accented capitals before a curly quote are valid UTF-8 pairs by chance
    (net3_cp1252_pair.inp, the review's second control: 4 such pairs, 9 bytes that are not UTF-8) still decodes as
    cp1252 throughout, as task 15 decoded it, with its IDs kept and no mixed-encoding warning; a content-named copy
    that task 15 accepted of a model blocked now (Chezy-Manning) is removed, and another file in the folder is kept.
    The second review's encoding cases (fixtures.ENCODING_CONTROLS), on both sides of the count rule it replaced
    (valid UTF-8 characters against stray bytes, by occurrence): in cp1252 throughout, one ID CITÉ'15 written 3 times
    against 1 stray byte, and five IDs with an en dash (18 against 2), decode as cp1252 with their IDs kept; in UTF-8,
    more stray bytes than accents, one capital-only ID against 6 stray bytes, a pasted block of 10 cp1252 lines, and a
    pasted comment holding a chance pair, decode as UTF-8 with their IDs kept (the pasted 'CITÉ'S' too).  The
    review's fuzz (seed 11: 150 files in cp1252 throughout, 150 in UTF-8 with stray bytes, decode_text only) changes
    no ID.  Equals the committed record (part 'leftovers')."""
    from . import intake_exp as E
    with _physics_fixtures() as (paths, out):
        got = E.part_leftovers(paths, out)
    s = got["swmm_after_64kb"]
    assert s["first_64kb_alone"] == "epanet_inp" and s["kind"] == "swmm_inp" and s["not_a_model"] and s["swmm_advice"], s
    z = got["zip_model_as_txt"]
    assert z["kind"] == "epanet_inp_zip" and z["ok"] and z["same_network_as_net3"], z
    m = got["utf8_stray_cp1252"]
    assert m["ok"] and m["mixed_warning"] and m["accented_ids_kept"] and m["decoded_as_cp1252_throughout"] != "loads", m
    p = got["cp1252_with_utf8_pairs"]      # the review's control: decoded as task 15 decoded it
    assert p["decoded_as"] == "cp1252 (Windows Western European)" and p["ok"] and p["ids_kept"], p
    assert p["same_as_cp1252_throughout"] and not p["mixed_warning"], p
    from .fixtures import ENCODING_CONTROLS
    ec = got["encoding_controls"]
    assert list(ec) == [f for f, _, _ in ENCODING_CONTROLS], list(ec)
    for f, v in ec.items():
        mixed = v["decoded_as"].startswith("utf-8 mixed with")
        assert v["ok"] and v["ids_kept"] and mixed == v["mixed_warning"] == (v["must_decode_as"] == "mixed"), (f, v)
    one = ec["net3_cp1252_one_pair_id.inp"]      # the review's file: the old rule's side it got wrong
    assert one["valid_utf8_chars"] > one["bytes_not_utf8"] and not one["mixed_warning"], one
    few = ec["utf8_one_id_six_stray.inp"]         # and the other side
    assert few["valid_utf8_chars"] < few["bytes_not_utf8"] and few["mixed_warning"], few
    assert ec["utf8_pasted_chance_pair.inp"]["pasted_comment_kept"], ec["utf8_pasted_chance_pair.inp"]
    fz = got["encoding_fuzz"]
    assert fz["ids_changed"] == {"cp1252_throughout": 0, "utf8_with_stray": 0}, fz
    st = got["stale_copy"]
    assert st["blocked"] and st["stale_removed"] and st["other_file_kept"], st
    _same_as_15b(got, "leftovers")
    return {"swmm": s["kind"], "zip_txt": z["kind"], "mixed": m["decoded_as"], "cp1252_with_pairs": p["decoded_as"],
            "stale_removed": True, "encoding_controls": {f: v["decoded_as"] for f, v in ec.items()},
            "encoding_fuzz_ids_changed": fz["ids_changed"]}


@check("task15b")
def numbers_scan_settles():
    """The second review: the numbers scan reads the checks report as a source (every entry but its own), so the run
    that added task 15b's entries committed a scan of the previous run's report, and the next run on the same tree
    rewrote it.  Reproduced on a copy of the report in a temporary folder (the committed one is not touched), with one
    probe entry that holds every number of the README's and journal's scanned sections, standing in for a run's new
    details: the probe changes the scan's result (the review's case), settle_numbers_entry, which the full run applies
    before it writes the report, puts the scan of the new report in the entry, and settling again changes nothing
    (one pass), so one run is enough.  Its first version updated the report's pass counts after the scan, and the
    counts are numbers in the scan's pool too: the full run caught it, when the journal came to quote the new count."""
    import copy
    from . import numbers as N
    with open(REPORT) as fh:
        base = json.load(fh)
    vals = sorted({N.parse_token(m.group(0))[0] for s in N.sections() for m in N.NUM.finditer(s.text)})
    with tempfile.TemporaryDirectory(prefix="rm_settle_") as d:
        path = os.path.join(d, "checks_report.json")
        settled = copy.deepcopy(base)
        settle_numbers_entry(settled, path=path)
        before = copy.deepcopy(next(r for r in settled["checks"] if r["name"] == "numbers_trace_to_outputs"))
        probed = copy.deepcopy(settled)
        probed["checks"].append({"name": "settle_probe", "group": "task15b", "pass": True, "detail": {"values": vals}})
        first = settle_numbers_entry(probed, path=path)
        after = next(r for r in probed["checks"] if r["name"] == "numbers_trace_to_outputs")
        second = settle_numbers_entry(probed, path=path)
    assert before["pass"] and after["pass"], (before.get("error"), after.get("error"))
    assert first["changed"] and before["detail"] != after["detail"], first
    assert second == {"rerun": True, "changed": False, "passes": 1}, second
    moved = sorted(k for k in after["detail"]["traced"] if after["detail"]["traced"][k] != before["detail"]["traced"].get(k))
    return {"probe_changes_the_scan": True, "second_pass_changes_nothing": True, "traced_levels_moved": moved}


@check("task15b")
def route_small_network():
    """plan_route on Net1 (9 junctions) from its fixed-gamma prior: with 8 junctions sampled and 6 sites asked for it
    returns the 1 junction left (before task 15b it stopped with 'attempt to get argmax of an empty sequence', at HEAD
    before task 15 too); with all 9 sampled, an empty route with the usual columns; with none sampled and 3 asked for,
    3 sites.  The committed default route is unchanged (route_reason_default_unchanged_age_opt_in compares Net3's with
    outputs/summary_Net3.json, and the committed reruns compare every route output).  Equals the committed record
    (part 'route')."""
    from . import intake_exp as E
    from .features import build_features
    from .route import plan_route
    from .simgp import SimGP24
    from .simulate import nominal_scenario
    got = E.part_route()
    c = got["cases"]
    assert (c["K6_8_sampled"]["rows"], c["K6_9_sampled"]["rows"], c["K3_none_sampled"]["rows"]) == (1, 0, 3), c
    assert c["K6_9_sampled"]["columns"] == c["K6_8_sampled"]["columns"], c
    sc = nominal_scenario("Net3")
    base = plan_route(SimGP24(sc, build_features(sc), seed=0, cache_dir=CACHE).fit_prior(), 8)
    committed = json.load(open(os.path.join(REPO, "outputs", "summary_Net3.json")))["route_scenario0_K8"]
    assert base.to_dict("records") == committed
    _same_as_15b(got, "route")
    return {"rows": {k: v["rows"] for k, v in c.items()}, "net3_default_route_unchanged": True}


APP_T15B_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
import streamlit.runtime.memory_media_file_storage as mm
PDF = {}
_orig = mm.MemoryMediaFileStorage.load_and_get_id
def _rec(self, path_or_data, mimetype, kind, filename=None):
    if mimetype == "application/pdf":
        PDF["last"] = path_or_data
    return _orig(self, path_or_data, mimetype, kind, filename)
mm.MemoryMediaFileStorage.load_and_get_id = _rec
from streamlit.testing.v1 import AppTest
app, net1, net6 = sys.argv[1], sys.argv[2], sys.argv[3]

def upload(name, path):
    PDF.clear()
    at = AppTest.from_file(app, default_timeout=900).run()
    [r for r in at.radio if r.label == "Model"][0].set_value("Upload my .inp")
    at.run()
    at.file_uploader[0].set_value((name, open(path, "rb").read(), "application/octet-stream"))
    at.run()
    metric = [m for m in at.metric if m.label.startswith("Junctions likely below")]
    calib = [m for m in at.metric if m.label.startswith("Calibrated decay")]
    route = [d.value for d in at.dataframe if "reason" in list(d.value.columns)]
    return {"exception": [e.value for e in at.exception], "flagged": metric[0].value if metric else None,
            "calib_value": calib[0].value if calib else None, "calib_delta": calib[0].delta if calib else None,
            "route_rows": len(route[0]) if route else None, "checkboxes": [c.label for c in at.checkbox],
            "info": [i.value for i in at.info], "errors": [e.value for e in at.error],
            "subheaders": [s.value for s in at.subheader], "pdf_made": "last" in PDF,
            "readiness": [m.value for m in at.markdown if m.value.startswith("**Model readiness")]}

print("APP_T15B_JSON " + json.dumps({"net1": upload("Net1.inp", net1), "net6": upload("Net6.inp", net6)}))
"""


@check("task15b", quick=False)
def app_small_network_upload():
    """Headless AppTest by widget label: Net1.inp uploaded (9 junctions, C 100 on every pipe) runs end to end with the
    default demo (8 samples), no exception: the 'full_fixed_gamma' grid is built inside the app, the calibrated-decay
    metric shows the old-pipe factor as not identifiable, the route has the 1 junction left, and the PDF is made (the
    route planner stopped this page before task 15b).  Net6.inp (3,323 junctions) shows its readiness report with the
    'Build the grid anyway' checkbox, unticked, and the page stops before any grid run (no grid file is written).  The
    runs' copies and grids (outputs/app_uploads and outputs/cache, git-ignored) are deleted afterwards, and earlier
    copies are moved aside first and put back."""
    from .simulate import LIB
    up_dir = os.path.join(REPO, "outputs", "app_uploads")
    had_dir = os.path.isdir(up_dir)
    mine = lambda d: [f for f in (os.listdir(d) if os.path.isdir(d) else []) if "_Net1.inp" in f or "_Net6.inp" in f]
    aside = tempfile.mkdtemp(prefix="rm_t15b_aside_")
    moved = [(d, f) for d in (up_dir, CACHE) for f in mine(d)]
    built = {}
    try:
        for d, f in moved:
            os.makedirs(os.path.join(aside, os.path.basename(d)), exist_ok=True)
            shutil.move(os.path.join(d, f), os.path.join(aside, os.path.basename(d), f))
        try:
            r = subprocess.run([PY, "-c", APP_T15B_SCRIPT, os.path.join(REPO, "app.py"), os.path.join(LIB, "Net1.inp"),
                                os.path.join(LIB, "Net6.inp")], cwd=os.getcwd(), capture_output=True, text=True,
                               env={**os.environ, "PYTHONPATH": REPO})
        finally:
            built = {os.path.basename(d): mine(d) for d in (up_dir, CACHE)}
            for d in (up_dir, CACHE):
                for f in mine(d):
                    os.remove(os.path.join(d, f))
            if not had_dir and os.path.isdir(up_dir) and not os.listdir(up_dir):
                os.rmdir(up_dir)
    finally:
        for d, f in moved:
            os.makedirs(d, exist_ok=True)
            shutil.move(os.path.join(aside, os.path.basename(d), f), os.path.join(d, f))
        shutil.rmtree(aside, ignore_errors=True)
    line = [x for x in r.stdout.splitlines() if x.startswith("APP_T15B_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-3000:]}")
    res = json.loads(line[-1][len("APP_T15B_JSON "):])
    n1, n6 = res["net1"], res["net6"]
    assert not n1["exception"] and not n1["errors"], (n1["exception"], n1["errors"])
    assert n1["flagged"] and n1["flagged"].endswith("of 9"), n1["flagged"]
    assert "not identifiable" in (n1["calib_delta"] or ""), n1["calib_delta"]
    assert n1["route_rows"] == 1 and n1["pdf_made"], n1
    assert "| Old-pipe factor | not identifiable" in n1["readiness"][0], n1["readiness"]
    grids = [f for f in built["cache"] if "_Net1.inp" in f]
    assert len(grids) == 1 and grids[0].startswith("grid24_full_fixed_gamma_"), built
    assert not n6["exception"] and any(c.startswith("Build the grid anyway") for c in n6["checkboxes"]), n6
    assert n6["flagged"] is None and any("tick 'Build the grid anyway'" in i for i in n6["info"]), n6
    assert not [f for f in built["cache"] if "_Net6.inp" in f], built
    assert "Large model:" in n6["readiness"][0], n6["readiness"]
    return {"net1_flagged": n1["flagged"], "net1_calib_delta_says": "old-pipe factor: not identifiable from one C value",
            "net1_route_rows": n1["route_rows"], "net1_grid": "grid24_full_fixed_gamma_<sha1>[:8]_Net1.inp_<tag>.pkl",
            "net6_checkbox": True, "net6_stopped_before_any_run": True}


# ----------------------------------------------------------------------------- task 16: grab logs, clock time, pages
INTAKE_RECORD_16 = os.path.join(REPO, "outputs", "intake", "task16_samples.json")
GET_STARTED_PAGE = os.path.join(REPO, "pages", "1_Get_started.py")


def _record_16(part: str) -> dict:
    """One part of the committed task-16 record (python -m residualmap.intake_exp samples)."""
    with open(INTAKE_RECORD_16) as fh:
        return json.load(fh)[part]


def _same_as_16(got: dict, part: str) -> dict:
    want = _record_16(part)
    diff = _json_diff(json.loads(json.dumps(got)), want)
    assert not diff, (part, diff[:10])
    return want


@contextlib.contextmanager
def _sample_fixtures():
    """The task-16 fixtures built into a temporary folder."""
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t16_") as d:
        yield FX.build_all(os.path.join(d, "fx"), FX.SAMPLE_FIXTURES)


def _app_run(script: str, args: list[str], match, label: str) -> tuple[dict, dict]:
    """Run an AppTest script in a subprocess (PYTHONPATH = the repo, cwd = the checks' scratch folder) and return its
    JSON line (printed after label) and the files it left in outputs/app_uploads and outputs/cache whose names
    match(name) is true for.  Those files (git-ignored) are deleted afterwards; earlier files with such names are
    moved aside first and put back, as the task-15 app checks do."""
    up_dir = os.path.join(REPO, "outputs", "app_uploads")
    had_dir = os.path.isdir(up_dir)
    mine = lambda d: [f for f in (os.listdir(d) if os.path.isdir(d) else []) if match(f)]
    aside = tempfile.mkdtemp(prefix="rm_t16_aside_")
    moved = [(d, f) for d in (up_dir, CACHE) for f in mine(d)]
    built = {}
    try:
        for d, f in moved:
            os.makedirs(os.path.join(aside, os.path.basename(d)), exist_ok=True)
            shutil.move(os.path.join(d, f), os.path.join(aside, os.path.basename(d), f))
        try:
            r = subprocess.run([PY, "-c", script, *args], cwd=os.getcwd(), capture_output=True, text=True,
                               env={**os.environ, "PYTHONPATH": REPO})
        finally:
            built = {os.path.basename(d): mine(d) for d in (up_dir, CACHE)}
            for d in (up_dir, CACHE):
                for f in mine(d):
                    os.remove(os.path.join(d, f))
            if not had_dir and os.path.isdir(up_dir) and not os.listdir(up_dir):
                os.rmdir(up_dir)
    finally:
        for d, f in moved:
            os.makedirs(d, exist_ok=True)
            shutil.move(os.path.join(aside, os.path.basename(d), f), os.path.join(d, f))
        shutil.rmtree(aside, ignore_errors=True)
    line = [x for x in r.stdout.splitlines() if x.startswith(label + " ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-3000:]}")
    return json.loads(line[-1][len(label) + 1:]), built


@check("task16")
def grab_log_column_mapping():
    """lims_one_datetime_us.csv (a sample ID per row, 'Sample Point', one 'Sample Date/Time' column of US dates with
    AM/PM times, 'Total Chlorine (mg/L)', a units column) maps with the guessed columns and an explicit date format
    (MDY, and the same rows with the strptime format), and without a format it is refused with the question, because
    every one of its dates reads two ways; lims_censored.csv ('<0.02', 'ND', '< 0.05', 'BDL', a blank, a decimal comma,
    '0.27 mg/L') gives its 3 readings, lists its 4 censored readings as censored and not fitted and drops the blank
    with its reason; log_semicolon.csv (';' and '0,45') and log_tab.txt give the rows of their table, with the guess
    and with an explicit mapping; ambiguous_dates.csv is refused with 'Which is it?' and read with DMY.  The expected
    rows come from each fixture's own table, never from the reader.  Equals the committed record (part 'mapping')."""
    from . import intake_exp as E
    with _sample_fixtures() as paths:
        got = E.part_mapping(paths)
    lu = got["lims_one_datetime_us"]
    assert lu["rows"] == 40 and lu["rows_as_expected"] and lu["with_strptime_same"], lu
    assert lu["mapping"] == {"tap": "Sample Point", "reading": "Total Chlorine (mg/L)", "date": "Sample Date/Time",
                             "time": None}, lu["mapping"]
    assert lu["without_format_refused"] and lu["without_format_asks"], lu
    lc = got["lims_censored"]
    assert lc["rows"] == 3 and lc["rows_as_expected"] and lc["as_expected"] and lc["censored_reasons_say_so"], lc
    for k in ("log_semicolon", "log_tab"):
        assert got[k]["rows"] == 4 and got[k]["rows_as_expected"], (k, got[k])
    assert got["log_semicolon_explicit_mapping_same"]
    a = got["ambiguous_dates"]
    assert a["refused"] and a["asks"] and a["with_DMY_rows"] == 3, a
    _same_as_16(got, "mapping")
    return {"lims_us_rows": lu["rows"], "lims_us_needs_a_date_format": True, "censored_listed": len(lc["censored"]),
            "censored_raw": lc["censored_raw"], "blank_dropped": len(lc["dropped"]), "semicolon_rows": 4, "tab_rows": 4,
            "ambiguous_refused_with_question": True}


@check("task16")
def grab_log_encodings():
    """grab_log_cp1252.csv (cp1252 throughout, Windows line ends, a tap ID 'Église_3' that the UTF-8 tap map holds)
    and grab_log_bom.csv (a UTF-8 byte-order mark before the first header) give the rows of their tables, the accented
    tap matched, and say how they were decoded.  Equals the committed record (part 'encodings')."""
    from . import intake_exp as E
    with _sample_fixtures() as paths:
        got = E.part_encodings(paths)
    c, b = got["grab_log_cp1252"], got["grab_log_bom"]
    assert c["rows"] == 3 and c["rows_as_expected"] and "cp1252" in c["how"], c
    assert b["rows"] == 2 and b["rows_as_expected"] and "byte-order mark" in b["how"] and b["header_first"] == "date", b
    _same_as_16(got, "encodings")
    return {"cp1252_rows": c["rows"], "bom_rows": b["rows"], "cp1252_how": c["how"], "bom_how": b["how"]}


@check("task16")
def excel_log_equals_csv_log():
    """grab_log_typed.xlsx (typed date and time cells; tap IDs 007 and 012 stored as the numbers 7 and 12, as Excel
    keeps them) gives the same S (junction, hour, reading) and the same rows as its CSV twin, whose IDs are text, with
    a warning that the two taps were matched by number; the CSV gets no warning.  The workbook is written by openpyxl,
    the library that reads it: same-library self-consistency.  Equals the committed record (part 'excel')."""
    from . import intake_exp as E
    with _sample_fixtures() as paths:
        got = E.part_excel(paths)
    assert got["xlsx_rows"] == got["csv_rows"] == 5 and got["same_rows"] and got["same_S"], got
    assert got["rows_as_expected"] and got["leading_zero_warning"] and got["csv_warnings"] == 0, got
    assert got["leading_zero_taps"] == ["007", "012"] and "12 is 012, 7 is 007" in got["warning"], got
    _same_as_16(got, "excel")
    return {"rows": got["xlsx_rows"], "same_S": True, "leading_zero_warning": True, "label": got["label"]}


@check("task16")
def analyzer_thinning():
    """analyzer_15min.csv (two online analyzers, 96 readings a day for 30 days, and three grab samples) keeps one
    reading per site per day, the one at 10:00 (600 minutes), 60 in all, and the three grab samples; every thinned
    reading is in the not-used list with its reason; at --analyzer-hour 8 the kept readings are at 08:00.  Equals the
    committed record (part 'analyzer')."""
    from . import intake_exp as E
    from . import fixtures as FX
    with _sample_fixtures() as paths:
        got = E.part_analyzer(paths)
    n_in = len(FX.ANALYZER_SITES) * FX.ANALYZER_DAYS * 24 * 60 // FX.ANALYZER_STEP_MIN
    assert got["readings_in"] == n_in and got["sites"] == list(FX.ANALYZER_SITES), got
    assert got["kept"] == len(FX.ANALYZER_SITES) * FX.ANALYZER_DAYS and got["days"] == FX.ANALYZER_DAYS, got
    assert got["per_site_per_day_max"] == got["per_site_per_day_min"] == 1 and got["clock_minutes_kept"] == [600], got
    assert got["grab_rows_kept"] == 3 and got["noted"] and got["hour_8"] == [480], got
    assert got["thinned_listed"] == got["readings_in"] - got["kept"], got        # every thinned reading is listed
    _same_as_16(got, "analyzer")
    return {"readings_in": got["readings_in"], "kept": got["kept"], "per_site_per_day": 1, "at": "10:00"}


@check("task16")
def utc_warning():
    """utc_shifted_log.csv (six daytime samples, 10:05 to 15:00, written 8 h late as a UTC export would) warns that the
    times may be UTC; the semicolon log and the committed example log do not.  Equals the committed record (part
    'utc')."""
    from . import intake_exp as E
    with _sample_fixtures() as paths:
        got = E.part_utc(paths)
    assert got["warned"] and "UTC" in got["warning"], got
    assert not got["semicolon_log_warned"] and not got["example_log_warned"], got
    _same_as_16(got, "utc")
    return {"warned": True, "warning_starts": got["warning"][:60]}


@check("task16")
def grab_log_review_cases():
    """The cases the review of task 16 found, each a fixture (part 'review' of the committed record; expected rows
    from each fixture's own table): a title line, a 'Lab:/Exported' block and Excel's 'sep=;' line before the header
    are skipped and said so; a 'Chloride (mg/L)' column before 'Chlorine Residual (mg/L)' is not taken as the reading;
    month-name dates and dotted times are read; after a quoted comment over two lines a row's number is its line;
    '>2.2' and '> 4.0' are kept apart from '<0.02' (kind 'above', never called below detection, in the reader and in
    the pilot's notes); a log in µg/L is divided by 1000, with a warning; 9999 and 450 mg/L are listed, not fitted; a
    long-format LIMS export keeps only its free chlorine rows (its total chlorine rows in chloramine mode) and lists the
    others with their reason; a log with a free and a total column asks which is the reading and reads only the one
    chosen; a junction list is no tap map; an analyzer's typed datetimes from 00:00 in Excel give its CSV twin's
    readings; a workbook's cover sheet is passed over and a stray cell far out is cut; a diurnal grab-sampling study
    keeps every reading; ten analyzers interleaved by time give a short not-used list that still counts every row;
    every reading is kept when asked; and rows before the date window are listed with their reason."""
    from . import intake_exp as E
    from . import fixtures as FX
    with _sample_fixtures() as paths:
        got = E.part_review(paths)
    for k in ("log_title_line", "log_lab_preamble", "log_sep_hint"):
        assert got[k]["rows"] == len(FX.REVIEW_ROWS) and got[k]["rows_as_expected"], (k, got[k])
    assert "skipped" in got["log_title_line"]["how"] and "'sep=;', names the separator" in got["log_sep_hint"]["how"]
    assert got["log_chloride_first"] == {"reading": "Chlorine Residual (mg/L)", "rows_as_expected": True}, got
    assert got["log_month_names"]["rows_as_expected"] and got["log_month_names"]["dropped"] == 0, got["log_month_names"]
    assert got["log_multiline_comment"]["dropped"] == [["T99", 4]], got["log_multiline_comment"]
    o = got["log_overrange"]
    assert o["censored"] == o["as_written"] and o["rows"] == 1, o
    assert any("1 reading below detection" in n for n in o["pilot_notes"]), o["pilot_notes"]
    assert any("2 readings above the measuring range" in n for n in o["pilot_notes"]), o["pilot_notes"]
    assert got["log_ugL"]["rows_as_expected"] and got["log_ugL"]["warned"], got["log_ugL"]
    li = got["log_implausible"]
    assert li["rows_as_expected"] and sorted(t for t, _ in li["dropped"]) == ["T1", "T3"], li
    for dis in ("free_chlorine", "chloramine"):
        lg = got["lims_long"][dis]
        assert lg["rows"] == len(FX.LONG_SAMPLES) and lg["rows_as_expected"] and lg["analyte"] == "Analyte", (dis, lg)
        assert sum(lg["not_used"].values()) == 3 * len(FX.LONG_SAMPLES) and not lg["species_warning"], (dis, lg)
    b = got["both_columns_chosen"]
    assert b["bare_guess_refused_asks"] and b["free_chlorine"]["y"] == [0.45, 0.38], b
    assert b["chloramine"]["y"] == [1.10, 0.95] and b["free_chlorine"]["other_column_noted"], b
    assert got["junction_list"] == {"tap": None, "junction": "Node ID"}, got["junction_list"]
    am = got["analyzer_midnight"]
    assert am["xlsx_readings"] == am["csv_readings"] == am["written"] and am["same_rows"], am
    assert am["no_time_dropped"] == 0, am
    cs = got["log_cover_and_stray"]
    assert cs["rows"] == 2 and len(cs["header"]) == 4 and "'Data'" in cs["how"], cs
    gs = got["grab_study_24h"]
    assert gs["rows"] == gs["written"] and gs["analyzer_taps"] == [] and gs["dropped"] == 0, gs
    ai = got["analyzer_interleaved"]
    assert ai["analyzer_taps"] == len(FX.INTERLEAVED_SITES) and ai["kept"] == len(FX.INTERLEAVED_SITES) * FX.INTERLEAVED_DAYS
    assert ai["thinned_listed"] == ai["readings_in"] - ai["kept"] and ai["lines_capped"] == ai["not_used_lines"], ai
    assert ai["chars"] < 5000, ai
    kr = got["analyzer_keep_every_reading"]
    assert kr["rows"] == len(FX.ANALYZER_SITES) * FX.ANALYZER_DAYS * 96 + 3 and kr["thinned_listed"] == 0, kr
    w = got["window_outside_listed"]
    assert w["n_outside"] == w["listed"] == w["rows"] - w["used"] and "before the window" in w["reason"], w
    _same_as_16(got, "review")
    return {"preamble_logs_read": 3, "overrange_kept_apart": o["censored"], "ugL_converted": True,
            "implausible_listed": len(li["dropped"]), "long_format_rows": got["lims_long"]["free_chlorine"]["rows"],
            "both_columns_chosen_reads_one": True, "grab_study_kept": gs["rows"],
            "interleaved_not_used_chars": ai["chars"], "window_outside_listed": w["listed"]}


@check("task16")
def grab_log_final_review_cases():
    """The cases the final verification of task 16 found, each a fixture (part 'final_review' of the committed
    record; expected rows from each fixture's own table): the species is read inside brackets, so a long-format export
    of 'Chlorine (Free)', 'Chlorine (Total)' and pH reads only the chosen species' rows (its 3 free chlorine rows; 3
    total chlorine rows in chloramine mode), and so do 'Free Cl' and 'Total Cl' with no units column and the LIMS
    codes FCL2 and TCL2 under 'Parameter Code'; two wide columns 'Chlorine (Free) mg/L' and 'Chlorine (Total) mg/L'
    (either order) and 'Chlorine Residual' beside 'Free Chlorine' ask which is the reading and read only the one
    chosen; 'Chloramine (mg/L)' is total chlorine (refused in free mode, read in chloramine mode); analytes that cannot
    pick the rows (codes the reader does not know, two chlorine analytes that name no species, a 'Parameter Code'
    column the guess passes over beside a plain 'Result') are refused with the question and read once the analyte
    value is chosen; µg/L headers written 'ugL', '_ugl', 'Âµg/L', 'ug per L' and 'µg L-1' are divided by 1000, units
    'mg Cl2/L', 'mgCl2/L' and 'mg/litre' are mg/L, and a µg/L log whose header does not say so warns; a comma title
    before a semicolon log and an Excel title row of two cells are skipped and named; an Excel date with no time among
    timed ones gives its CSV twin's rows, listed as 'no time of day'; a time column holding '1899-12-30 09:15:00'
    reads its times; a tap ID with a superscript digit is listed, never a crash; 'Node ID,Elevation' is no tap map;
    a lab qualifier on every reading of ten interleaved analyzers and 200 taps from another system give one not-used
    line each, and the list is capped at DROPPED_MAX_LINES lines; a 9 MB file of 3 million tiny rows is refused, by
    its line count, within 1 s; and the bundled Net2's 08:00 start moves a 09:00 sample to model hour 1 in the pilot
    while the app's bundled Net2 keeps it at 9."""
    from . import intake_exp as E
    from . import fixtures as FX
    from .importers import DROPPED_MAX_LINES
    with _sample_fixtures() as paths:
        got = E.part_final_review(paths)
    for f, res in got["long_species"].items():
        for dis, r in res.items():
            assert r["rows"] == len(FX.LONG_SAMPLES) and r["rows_as_expected"] and not r["species_warning"], (f, dis, r)
            assert r["analyte"] in ("Analyte", "Parameter Code") and r["reading"] == "Result", (f, dis, r)
    for k, r in got["wide_species"].items():
        assert r["bare_guess_refused_asks"] and r["rows_as_expected"] and r["other_column_noted"], (k, r)
    ch = got["chloramine_header"]
    assert ch["free_mode_refused_as_total"] and ch["chloramine_rows_as_expected"] and not ch["species_warning"], ch
    av = got["analyte_value"]
    for k in ("long_unknown_codes", "long_two_unnamed"):
        assert av[k]["refused_asks_value"] and av[k]["rows_as_expected"] and av[k]["warned"], (k, av[k])
    assert av["long_unknown_codes_param"]["refused_names_the_column"] and av["long_unknown_codes_param"]["rows_as_expected"]
    assert av["value_not_in_column_refused"], av
    u = got["units"]
    for h in FX.UGL_HEADERS:
        assert u[h]["y"] == [round(float(v) / 1000, 6) for *_, v in FX.UGL_ROWS] and u[h]["warned"], (h, u[h])
    for x in FX.UNIT_SPELLINGS:
        assert u[x]["rows_as_expected"] and u[x]["dropped"] == 0, (x, u[x])
    assert u["unnamed_ugL_warned"], u
    ti = got["titles"]
    assert ti["comma_title_rows_as_expected"] and "semicolons" in ti["comma_title_how"] and "skipped" in ti["comma_title_how"]
    assert ti["xlsx_title_rows_as_expected"] and "skipped" in ti["xlsx_title_how"], ti
    ed = got["excel_date_only"]
    assert ed["same_rows_as_csv"] and len(ed["xlsx_rows"]) == 2 and ed["xlsx_no_time"] == ed["csv_no_time"] == [["T3", 4]]
    assert got["time_column_epoch_date"] == {"rows_as_expected": True, "dropped": 0}, got["time_column_epoch_date"]
    sp = got["superscript_tap"]
    assert "error" not in sp["log"] and "error" not in sp["tap_map"], sp
    assert sp["log"]["dropped"] == [["1²", "tap 1² is not in the tap map"]] and sp["tap_map"]["rows"] == 2, sp
    assert got["junctions_2col"] == {"tap": None, "junction": "Node ID"}, got["junctions_2col"]
    nl = got["not_used_lists"]
    for k in ("qualifier_interleaved", "unknown_taps_many"):
        assert nl[k]["lines"] == 1 and nl[k]["not_used"] == nl[k]["rows_in"] and nl[k]["chars"] < 300, (k, nl[k])
    assert nl["line_cap"]["lines"] == DROPPED_MAX_LINES and nl["line_cap"]["last"].startswith("and 6 more reasons")
    tr = got["csv_tiny_rows"]
    assert tr["refused"] and "lines, more than" in tr["refused"] and tr["within_1s"], tr
    n2 = got["net2_pilot_hours"]
    assert n2["start_clock_h"] == 8.0 and n2["app_hours"] == n2["clock_hours"], n2
    assert n2["pilot_hours"] == [h - 8 for h in n2["clock_hours"]], n2
    _same_as_16(got, "final_review")
    return {"long_species_rows": {k: v["free_chlorine"]["rows"] for k, v in got["long_species"].items()},
            "wide_both_species_ask": len(got["wide_species"]), "analyte_value_asked_then_read": 3,
            "ugl_headers_divided": len(FX.UGL_HEADERS), "unit_spellings_read": len(FX.UNIT_SPELLINGS),
            "not_used_lines": {k: nl[k]["lines"] for k in ("qualifier_interleaved", "unknown_taps_many")},
            "tiny_rows_refused_by_line_count": True, "net2_pilot_hours": n2["pilot_hours"]}


APP_CLOCK_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
app, net3_log, clock6_inp, clock6_log, taps, net2_log, net2_taps = sys.argv[1:8]
EDITOR = ["junction", "hour", "mg/L"]

def run(at):
    at.run()
    assert not at.exception, [e.value for e in at.exception]
    return at

def put(at, label, name, path):
    [u for u in at.file_uploader if u.label == label][0].set_value((name, open(path, "rb").read(), "text/csv"))

def state(at):
    m = {x.label: x.value for x in at.metric}
    route = [d.value for d in at.dataframe if "reason" in d.value.columns][0]
    every = [d.value for d in at.dataframe if "hour of the daily minimum" in d.value.columns][0]
    ed = [d.value for d in at.dataframe if list(d.value.columns) == EDITOR][0]
    flagged = [k for k in m if k.startswith("Junctions likely below")][0]
    return {"worst": m["Worst hour of the day"], "flagged_label": flagged, "flagged": m[flagged],
            "samples": m["Samples used"], "route_hours": list(route.hour),
            "every_min_hour": dict(zip(every.junction.astype(str), every["hour of the daily minimum"])),
            "editor": [[str(a), int(b), float(c)] for a, b, c in ed.values.tolist()],
            "show_options": [r for r in at.radio if r.label == "Show"][0].options,
            "captions": [c.value for c in at.caption], "info": [i.value for i in at.info],
            "warnings": [w.value for w in at.warning]}

def show(at, opt):
    [r for r in at.radio if r.label == "Show"][0].set_value(opt)
    run(at)
    m = {x.label: x.value for x in at.metric}
    k = [k for k in m if k.startswith("Junctions likely below")][0]
    return {"label": k, "value": m[k]}

res = {}
at = run(AppTest.from_file(app, default_timeout=900))                     # A: Net3, demo off, the log at clock hours
[t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False); run(at)
put(at, "Grab log (CSV or Excel)", "clock_log_net3.csv", net3_log); put(at, "Tap map (CSV)", "task16_tap_map.csv", taps)
run(at); res["net3"] = state(at); res["net3_at_03"] = show(at, "03:00")
at = run(AppTest.from_file(app, default_timeout=900))                     # B: net3_clock6.inp uploaded, demo off
[r for r in at.radio if r.label == "Model"][0].set_value("Upload my .inp"); run(at)
at.file_uploader[0].set_value(("net3_clock6.inp", open(clock6_inp, "rb").read(), "application/octet-stream")); run(at)
res["clock6_demo"] = state(at)                                            # the demo on: its hours by the clock
[t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False); run(at)
res["uploaders"] = [u.label for u in at.file_uploader]
put(at, "Grab log (CSV or Excel)", "clock_log_clock6.csv", clock6_log); put(at, "Tap map (CSV)", "task16_tap_map.csv", taps)
run(at); res["clock6"] = state(at); res["clock6_at_09"] = show(at, "09:00")
res["clock6_readiness"] = [m.value for m in at.markdown if m.value.startswith("**Model readiness")]
at = run(AppTest.from_file(app, default_timeout=900))                     # C: the Net2 example, demo on
at.selectbox[0].set_value([o for o in at.selectbox[0].options if o.startswith("Net2")][0]); run(at)
res["net2"] = state(at)
[t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False); run(at)   # D: and a log at clock hours
put(at, "Grab log (CSV or Excel)", "net2_clock_log.csv", net2_log); put(at, "Tap map (CSV)", "net2_tap_map.csv", net2_taps)
run(at); res["net2_log"] = state(at)
print("APP_CLOCK_JSON " + json.dumps(res))
"""


@check("task16", quick=False)
def clock_conversion():
    """Clock time in and out (critique 10).  The helpers: a 09:00 sample in a file that starts at 06:00 is model hour 3
    (09:59 too: minutes truncated), the daytime window is model hours 1 to 11, shown as 07:00 to 17:00, a start of
    6:30 gives 07:30 to 17:30, and a start at midnight gives simgp.DAY_HOURS and leaves every hour as it is; the
    bundled Net3 and ky4 start at midnight, so nothing of theirs moves, and Net2 (08:00) keeps model hours by lead
    decision 2.  In the app (AppTest): (A) the Net3 example with clock_log_net3.csv, whose times of day equal its
    model hours, and (B) net3_clock6.inp (Net3 starting at 06:00, its grid bit-identical to Net3's) with
    clock_log_clock6.csv, the same readings written 6 h later.  B's sample table holds A's rows exactly (its 09:00
    readings are model hour 3), so the fit is A's: B's worst hour is A's plus 6, mod 24, every junction's hour of the
    daily minimum is A's plus 6, B's map at 09:00 is A's at 03:00, and B's route hours lie within 07:00 to 17:00 by the
    clock; A shows model hours as before.  (C) the Net2 example keeps its model hours, with the caption 'model hour
    (this file starts at 08:00)'.  Since the review of task 16, (B) with the demo on draws its daytime samples at the
    model hours of 07:00 to 17:00 by the clock (1 to 11), not at DAY_HOURS shown as 13:00 to 23:00.  The upload's copy
    and grid are deleted afterwards.  ky4's app run is not repeated here (its free-chlorine grid is not cached); its
    start is midnight, so its hours are unchanged by construction.  Since the final verification of task 16, (D) the
    Net2 example with the demo off and net2_clock_log.csv (09:00 to 15:00) keeps those times as model hours 9 to 15,
    unconverted (lead decision 2), and its caption now says so and that the pilot with --inp Net2 subtracts 8 h; the
    demo-on Net2 page (C) has no such words."""
    from . import intake_exp as E
    from . import fixtures as FX
    from .simulate import LIB
    h = E.part_clock_helpers()
    assert h["09:00_at_6h"] == 3 and h["09:59_at_6h"] == 3 and h["model_3_label_6h"] == "09:00", h
    assert h["daytime_6h"] == list(range(1, 12)) and h["daytime_6h_labels"][0] == "07:00", h
    assert h["daytime_6h_labels"][-1] == "17:00" and h["daytime_6h30_labels"][0] == "07:30", h
    assert h["daytime_0h_is_DAY_HOURS"] and h["identity_at_0h"], h
    assert h["bundled_start_h"] == {"Net3": 0.0, "Net2": 8.0, "ky4": 0.0}, h
    assert h["shows_clock"] == {"0": False, "6": True, "23:59:50": False, "0:00:20": False}, h
    _same_as_16(h, "clock_helpers")
    with _sample_fixtures() as paths:
        res, built = _app_run(APP_CLOCK_SCRIPT, [os.path.join(REPO, "app.py"), paths["clock_log_net3.csv"],
                                                 paths["net3_clock6.inp"], paths["clock_log_clock6.csv"],
                                                 paths["task16_tap_map.csv"], paths["net2_clock_log.csv"],
                                                 paths["net2_tap_map.csv"]],
                              lambda f: "_net3_clock6.inp" in f, "APP_CLOCK_JSON")
    a, b, c = res["net3"], res["clock6"], res["net2"]
    plus6 = lambda lab: f"{(int(lab[:2]) + FX.CLOCK_SHIFT_H) % 24:02d}:00"
    hours = [f"{x:02d}:00" for x in range(24)]
    assert res["uploaders"] == ["EPANET .inp file", "Grab log (CSV or Excel)", "Tap map (CSV)"], res["uploaders"]
    assert b["editor"] == a["editor"] and a["samples"] == b["samples"], (a["editor"][:4], b["editor"][:4])
    first = [r for r in b["editor"] if r[1] == 3]
    assert first and a["show_options"][1:] == hours and b["show_options"][1:] == hours, (first, b["show_options"])
    assert b["worst"] == plus6(a["worst"]), (a["worst"], b["worst"])
    assert all(b["every_min_hour"][j] == plus6(v) for j, v in a["every_min_hour"].items()), "hour of the daily minimum"
    assert res["clock6_at_09"]["value"] == res["net3_at_03"]["value"], (res["clock6_at_09"], res["net3_at_03"])
    assert res["clock6_at_09"]["label"].endswith("(09:00)") and res["net3_at_03"]["label"].endswith("(03:00)")
    assert b["route_hours"] and all(7 <= int(x[:2]) <= 17 for x in b["route_hours"]), b["route_hours"]
    assert all(7 <= int(x[:2]) <= 17 for x in a["route_hours"]), a["route_hours"]
    assert any("Hours are times of day: this file starts at 06:00" in x for x in b["captions"]), b["captions"]
    assert any("MODEL hour" in x and "is hour 3" in x for x in b["captions"]), b["captions"]
    assert not any("times of day" in x for x in a["captions"]), a["captions"]
    assert any("model hour (this file starts at 08:00)" in x for x in c["captions"]), c["captions"]
    assert c["show_options"][1:] == hours and all(7 <= int(x[:2]) <= 17 for x in c["route_hours"]), c["route_hours"]
    assert "Times of day are converted to model hours" in " ".join(b["info"]), b["info"]
    assert any(f"grid24_full_" in f for f in built["cache"]), built
    # the review of task 16: on a converted upload the demo draws its daytime samples at 07:00 to 17:00 by the clock
    demo_h = [r[1] for r in res["clock6_demo"]["editor"]]
    assert demo_h and all(x in h["daytime_6h"] for x in demo_h), demo_h
    # the final verification of task 16: the bundled Net2 with an uploaded log says its times are not converted
    d = res["net2_log"]
    assert [r[1] for r in d["editor"]] == list(range(9, 9 + len(FX.NET2_LOG_JUNCTIONS))), d["editor"]
    assert any("not converted (the pilot with --inp Net2 subtracts 8 h)" in x for x in d["captions"]), d["captions"]
    assert not any("not converted" in x for x in c["captions"]), c["captions"]
    return {"helpers": {k: h[k] for k in ("09:00_at_6h", "daytime_6h_labels", "daytime_0h_is_DAY_HOURS")},
            "net3_worst": a["worst"], "clock6_worst": b["worst"], "net3_at_03": res["net3_at_03"]["value"],
            "clock6_at_09": res["clock6_at_09"]["value"], "clock6_route_hours": b["route_hours"],
            "net3_route_hours": a["route_hours"], "samples_used": int(a["samples"]),
            "junctions_whose_daily_min_hour_moved_by_6": len(a["every_min_hour"]), "net2_caption": True,
            "net2_worst_model_hour": c["worst"],
            "clock6_demo_sample_clock_hours": sorted({(x + FX.CLOCK_SHIFT_H) % 24 for x in demo_h}),
            "net2_upload_hours_unconverted": [r[1] for r in d["editor"]], "net2_upload_caption_says_so": True}


APP_LOG_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
import matplotlib.axes
_scatter = matplotlib.axes.Axes.scatter
MARKS = []
ABOVE = []
def _rec(self, *a, **k):
    if str(k.get("label") or "").startswith("measured below detection"):
        MARKS.append(len(a[0]))
    if str(k.get("label") or "").startswith("measured above the range"):
        ABOVE.append(len(a[0]))
    return _scatter(self, *a, **k)
matplotlib.axes.Axes.scatter = _rec
from streamlit.testing.v1 import AppTest
app, ex_log, ex_taps, fx = sys.argv[1], sys.argv[2], sys.argv[3], json.loads(sys.argv[4])
EDITOR = ["junction", "hour", "mg/L"]

def run(at):
    at.run()
    return at

def put(at, label, name, path):
    [u for u in at.file_uploader if u.label == label][0].set_value((name, open(path, "rb").read(), "application/octet-stream"))

def demo_off(ca=False):
    at = run(AppTest.from_file(app, default_timeout=900))
    [t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False)
    run(at)
    if ca:      # after the demo is off, so the chloramine demo truth (EPANET-MSX) never runs here
        [r for r in at.radio if r.label == "Disinfectant"][0].set_value("Chloramine (total chlorine)")
        run(at)
    return at

def state(at):
    m = {x.label: x.value for x in at.metric}
    ed = [d.value for d in at.dataframe if list(d.value.columns) == EDITOR]
    return {"exception": [e.value for e in at.exception], "samples": m.get("Samples used"),
            "editor": [[str(a), int(b), float(c)] for a, b, c in ed[0].values.tolist()] if ed else None,
            "info": [i.value for i in at.info], "warnings": [w.value for w in at.warning],
            "errors": [e.value for e in at.error], "uploaders": [u.label for u in at.file_uploader],
            "selectboxes": [s.label for s in at.selectbox], "downloads": len(at.get("download_button")),
            "flagged": next((v for k, v in m.items() if k.startswith("Junctions likely below")), None),
            "reading_choice": next((s.value for s in at.selectbox if s.label == "Reading column"), None)}

def upload(log_name, log_path, taps_name, taps_path, ca=False):
    at = demo_off(ca)
    put(at, "Grab log (CSV or Excel)", log_name, log_path)
    put(at, "Tap map (CSV)", taps_name, taps_path)
    return run(at)

res = {}
at = run(AppTest.from_file(app, default_timeout=900))                     # the demo path: no log uploaders
res["demo_on"] = state(at)
at = upload("example_grab_log.csv", ex_log, "example_tap_map.csv", ex_taps)
res["example"] = state(at)
put(at, "Grab log (CSV or Excel)", "multi_month_log.csv", fx["multi_month_log.csv"])      # a changed re-upload
res["reupload_multi_month"] = state(run(at))
res["unknown_tap"] = state(upload("unknown_tap_log.csv", fx["unknown_tap_log.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["mixed"] = state(upload("mixed_species_log.csv", fx["mixed_species_log.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["xlsx"] = state(upload("grab_log_typed.xlsx", fx["grab_log_typed.xlsx"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["xlsx_csv"] = state(upload("grab_log_typed.csv", fx["grab_log_typed.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
MARKS.clear()
res["censored"] = state(upload("lims_censored.csv", fx["lims_censored.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["censored_marks"] = list(MARKS)
res["ambiguous"] = state(upload("ambiguous_dates.csv", fx["ambiguous_dates.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
at = upload("ambiguous_dates.csv", fx["ambiguous_dates.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"])
[s for s in at.selectbox if s.label == "Date format"][0].set_value("day first (DD/MM/YYYY)")
res["ambiguous_dmy"] = state(run(at))
at = demo_off()
put(at, "Grab log (CSV or Excel)", "example_grab_log.csv", ex_log)
res["no_tap_map"] = state(run(at))
MARKS.clear(); ABOVE.clear()                                               # the review of task 16
res["overrange"] = state(upload("log_overrange.csv", fx["log_overrange.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["overrange_marks"] = {"below": list(MARKS), "above": list(ABOVE)}
res["long"] = state(upload("lims_long.csv", fx["lims_long.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
at = upload("mixed_species_log.csv", fx["mixed_species_log.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"])
[s for s in at.selectbox if s.label == "Reading column"][0].set_value("free_chlorine_mgL")
res["mixed_chosen"] = state(run(at))
# the final verification of task 16: the species inside brackets; the question for two species; the analyte value
res["long_paren"] = state(upload("long_paren.csv", fx["long_paren.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"]))
res["wide_paren_ca"] = state(upload("wide_paren.csv", fx["wide_paren.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"],
                                    ca=True))
at = upload("long_unknown_codes.csv", fx["long_unknown_codes.csv"], "task16_tap_map.csv", fx["task16_tap_map.csv"])
res["unknown_codes"] = state(at)
[s for s in at.selectbox if s.label.startswith("Analyte value")][0].set_value("CHLF")
res["unknown_codes_chosen"] = state(run(at))
print("APP_LOG_JSON " + json.dumps(res))
"""


@check("task16", quick=False)
def app_grab_log_upload():
    """Headless AppTest of the grab-log upload, the demo off, the Net3 example (task 16's acceptance: from a grab-log
    file, CSV or Excel, and a tap map by ID to the chlorine map without typing a junction ID):
    docs/example_grab_log.csv with docs/example_tap_map.csv makes 'Samples used' equal load_log's row count inside the
    date window (the last 3 of its 6 months), and says 'using N samples from <date> to <date>' and that within the
    window every sample is one simulated day; a changed re-upload (multi_month_log.csv, 9 months) changes the table
    and uses its last 3 months, and says so; unknown_tap_log.csv lists its unknown tap in a warning; mixed_species_log.csv
    is refused with an error before any map; grab_log_typed.xlsx gives the same table as its CSV twin, with the
    leading-zero warning; lims_censored.csv lists its censored readings with their reason and marks their junctions on
    the map ('measured below detection (not fitted)', counted through matplotlib's scatter); ambiguous_dates.csv
    stops with the question, and runs once 'day first' is chosen; a log with no tap map asks for it.  The demo path
    has no grab-log uploader, and every page keeps one download button.  Since the review of task 16: the 9-month
    log's 24 rows before the window are counted on the page and listed with their reason; log_overrange.csv marks
    its '<0.02' junction as below detection and its two '>' junctions apart, as above the range; lims_long.csv (one
    row per analyte) uses its 3 free chlorine rows, with the analyte column guessed; and mixed_species_log.csv runs
    once its 'Reading column' is chosen.  Since its final verification: long_paren.csv ('Chlorine (Free)',
    'Chlorine (Total)', pH rows) uses its 3 free chlorine rows with no species warning; wide_paren.csv in chloramine
    mode stops with 'Which column is the reading?' and its 'Reading column' starts at '(choose a column)' (it read the
    free column as total chlorine before); and long_unknown_codes.csv stops with 'Which value is the chlorine
    reading?' and uses its 3 CHLF rows once 'Analyte value' is CHLF."""
    from . import intake_exp as E
    from . import fixtures as FX
    want = _record_16("app_inputs")
    with _sample_fixtures() as paths:
        got_inputs = E.part_app_inputs(paths)
        res, _ = _app_run(APP_LOG_SCRIPT, [os.path.join(REPO, "app.py"), os.path.join(REPO, "docs", "example_grab_log.csv"),
                                           os.path.join(REPO, "docs", "example_tap_map.csv"), json.dumps(paths)],
                          lambda f: False, "APP_LOG_JSON")
    assert not _json_diff(json.loads(json.dumps(got_inputs)), want), _json_diff(got_inputs, want)[:5]
    for k, v in res.items():
        if isinstance(v, dict) and "exception" in v:
            assert not v["exception"], (k, v["exception"])
    d0 = res["demo_on"]
    assert not any(u.startswith(("Grab log", "Tap map")) for u in d0["uploaders"]) and d0["downloads"] == 1, d0
    ex = res["example"]
    assert int(ex["samples"]) == want["example_rows_in_window"] == 39, ex["samples"]
    assert any(i.startswith("Grab log: using 39 samples from 2026-04-03 to 2026-06-27") and "one simulated day" in i
               and "the last 3 months of the log" in i for i in ex["info"]), ex["info"]
    assert ex["downloads"] == 1 and not ex["errors"], ex
    mm = res["reupload_multi_month"]
    assert mm["editor"] != ex["editor"] and int(mm["samples"]) == want["multi_month"]["rows_used"], (mm["samples"],)
    assert any("the last 3 months of the log (2026-07, 2026-08, 2026-09" in i for i in mm["info"]), mm["info"]
    u = res["unknown_tap"]
    assert any("tap T99 is not in the tap map: row 3" in w for w in u["warnings"]) and int(u["samples"]) == 3, u
    mx = res["mixed"]
    assert any("both a free chlorine column" in e for e in mx["errors"]) and mx["flagged"] is None, mx
    x, xc = res["xlsx"], res["xlsx_csv"]
    assert x["editor"] == xc["editor"] and int(x["samples"]) == 5, (x["editor"], xc["editor"])
    assert any("leading zeros" in w for w in x["warnings"]) and not any("leading zeros" in w for w in xc["warnings"]), x
    c = res["censored"]
    cens = [w for w in c["warnings"] if w.startswith("Rows of the grab log that are not used")]
    assert cens and cens[0].count("censored") == 4 and "blank reading" in cens[0], cens
    assert int(c["samples"]) == 3 and res["censored_marks"] and set(res["censored_marks"]) == {4}, res["censored_marks"]
    am = res["ambiguous"]
    assert any("Which is it?" in e for e in am["errors"]) and am["flagged"] is None, am["errors"]
    assert int(res["ambiguous_dmy"]["samples"]) == 3 and not res["ambiguous_dmy"]["errors"], res["ambiguous_dmy"]
    assert any("Upload the tap map too" in i for i in res["no_tap_map"]["info"]), res["no_tap_map"]["info"]
    # the review of task 16: rows before the window listed; '>2.2' marked apart from '<0.02'; a long-format export;
    # a log with both species read once its reading column is chosen
    out_n = want["multi_month"]["rows_outside"]
    assert any(f"{out_n} rows dated before 2026-07 are outside the window" in i for i in mm["info"]), mm["info"]
    assert any(w.startswith("Rows of the grab log that are not used") and "dated before the window" in w
               for w in mm["warnings"]), mm["warnings"]
    ov = res["overrange"]
    marks = {k: sorted(set(v)) for k, v in res["overrange_marks"].items()}      # the page's map and the PDF's
    assert int(ov["samples"]) == 1 and marks == {"below": [1], "above": [2]}, res["overrange_marks"]
    lg = res["long"]
    assert int(lg["samples"]) == 3 and not lg["errors"] and "Analyte column (a log with one row per analyte)" in lg["selectboxes"], lg
    assert any("another analyte ('pH')" in w for w in lg["warnings"]), lg["warnings"]
    mc = res["mixed_chosen"]
    assert int(mc["samples"]) == 2 and not mc["errors"] and not mc["exception"], mc
    # the final verification of task 16
    lp = res["long_paren"]
    assert int(lp["samples"]) == len(FX.LONG_SAMPLES) and not lp["errors"], lp
    assert not any("does not say free or total" in w for w in lp["warnings"]), lp["warnings"]
    assert any("total chlorine, not free chlorine ('Chlorine (Total)')" in w for w in lp["warnings"]), lp["warnings"]
    wp = res["wide_paren_ca"]
    assert any("Which column is the reading?" in e for e in wp["errors"]) and wp["flagged"] is None, wp
    assert wp["reading_choice"] == "(choose a column)", wp["reading_choice"]
    uc, ucc = res["unknown_codes"], res["unknown_codes_chosen"]
    assert any("Which value is the chlorine reading?" in e for e in uc["errors"]) and uc["flagged"] is None, uc
    assert any(x.startswith("Analyte value") for x in uc["selectboxes"]), uc["selectboxes"]
    assert int(ucc["samples"]) == len(FX.LONG_SAMPLES) and not ucc["errors"] and not ucc["exception"], ucc
    return {"samples_used_example": int(ex["samples"]), "load_log_rows_in_window": want["example_rows_in_window"],
            "load_log_rows": want["example_rows"], "multi_month_samples_used": int(mm["samples"]),
            "unknown_tap_listed": True, "mixed_species_refused": True, "xlsx_equals_csv": True,
            "censored_listed": 4, "censored_junctions_marked_on_the_map": res["censored_marks"][0],
            "ambiguous_asks_then_runs_with_day_first": True, "download_buttons": ex["downloads"],
            "multi_month_rows_before_the_window_listed": out_n, "overrange_marked_apart": marks,
            "long_format_samples_used": int(lg["samples"]), "both_species_chosen_samples_used": int(mc["samples"]),
            "species_in_brackets_samples_used": int(lp["samples"]), "wide_species_in_brackets_asks": True,
            "analyte_value_asked_then_samples_used": int(ucc["samples"])}


APP_TABLE_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
app, ex_log, ex_taps = sys.argv[1:4]
def table_state(at):
    assert not at.exception, [e.value for e in at.exception]
    every = [d.value for d in at.dataframe if "hour of the daily minimum" in d.value.columns]
    plain = [t.value for t in at.table if "hour of the daily minimum" in t.value.columns]
    m = {x.label: x.value for x in at.metric}
    return {"every_shape": list(every[0].shape) if every else None, "every_columns": list(every[0].columns) if every else None,
            "plain_rows": len(plain[0]) if plain else None, "plain_index": plain[0].index.name if plain else None,
            "plain_top_p": [float(v) for v in plain[0].iloc[:, 3]] if plain else None,
            "downloads": len(at.get("download_button")),
            "flagged": next((v for k, v in m.items() if k.startswith("Junctions likely below")), None),
            "subheaders": [s.value for s in at.subheader], "captions": [c.value for c in at.caption]}
res = {}
at = AppTest.from_file(app, default_timeout=900).run()
res["default"] = table_state(at)
[t for t in at.toggle if t.label.startswith("Demo")][0].set_value(False); at.run()
[u for u in at.file_uploader if u.label == "Grab log (CSV or Excel)"][0].set_value(("example_grab_log.csv", open(ex_log, "rb").read(), "text/csv"))
[u for u in at.file_uploader if u.label == "Tap map (CSV)"][0].set_value(("example_tap_map.csv", open(ex_taps, "rb").read(), "text/csv"))
at.run()
res["demo_off_log"] = table_state(at)
print("APP_TABLE_JSON " + json.dumps(res))
"""


@check("task16", quick=False)
def app_per_junction_table():
    """Headless AppTest (critique 25): with the default demo and with the demo off and the example log uploaded, the
    'Every junction' dataframe has one row per junction (92 on Net3) and the 6 columns junction, median, lo90, hi90
    (mg/L), P(daily min < threshold) and the hour of the daily minimum; the plain table of the junctions most likely
    below the threshold has 20 rows, sorted by that probability; the page keeps exactly one download_button (the
    dataframe's own toolbar offers the CSV); and the default demo still flags 42 of 92.  No claim of WCAG conformance
    is made; the dataframe toolbar's CSV download was checked once in a real browser (journal, task 16)."""
    res, _ = _app_run(APP_TABLE_SCRIPT, [os.path.join(REPO, "app.py"), os.path.join(REPO, "docs", "example_grab_log.csv"),
                                         os.path.join(REPO, "docs", "example_tap_map.csv")], lambda f: False,
                      "APP_TABLE_JSON")
    cols = ["junction", "median (mg/L)", "lo90 (mg/L)", "hi90 (mg/L)", "P(daily min < 0.2)", "hour of the daily minimum"]
    for k in ("default", "demo_off_log"):
        g = res[k]
        assert g["every_shape"] == [APP_DEMO["junctions"], 6] and g["every_columns"] == cols, (k, g)
        assert g["plain_rows"] == 20 and g["plain_index"] == "junction", (k, g)
        assert g["plain_top_p"] == sorted(g["plain_top_p"], reverse=True), (k, g["plain_top_p"])
        assert g["downloads"] == 1 and "Every junction: the daily minimum" in g["subheaders"], (k, g)
        assert any("downloads it as CSV" in c and "screen readers" in c for c in g["captions"]), k
    assert res["default"]["flagged"] == f"{APP_DEMO['flagged']} of {APP_DEMO['junctions']}", res["default"]["flagged"]
    return {"rows": res["default"]["every_shape"][0], "columns": cols, "plain_table_rows": 20, "download_buttons": 1,
            "default_flagged": res["default"]["flagged"]}


PAGE_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
import streamlit.runtime.memory_media_file_storage as mm
FILES = {}
_orig = mm.MemoryMediaFileStorage.load_and_get_id
def _rec(self, path_or_data, mimetype, kind, filename=None):
    FILES[filename] = (mimetype, len(path_or_data) if isinstance(path_or_data, (bytes, bytearray)) else None)
    return _orig(self, path_or_data, mimetype, kind, filename)
mm.MemoryMediaFileStorage.load_and_get_id = _rec
from streamlit.testing.v1 import AppTest
page, fx = sys.argv[1], json.loads(sys.argv[2])
at = AppTest.from_file(page, default_timeout=900).run()
assert not at.exception, [e.value for e in at.exception]
res = {"downloads": [b.label for b in at.get("download_button")], "files": FILES, "choices": {}}
radio = lambda: [r for r in at.radio if r.label == "Which software holds your network?"][0]
for o in radio().options:
    radio().set_value(o)
    at.run()
    assert not at.exception, (o, [e.value for e in at.exception])
    res["choices"][o] = [m.value for m in at.markdown][0]
for name in ("mains.gpkg", "net3_backup.inp.bak", "lims_censored.csv", "task16_tap_map.csv", "ambiguous_dates.csv",
             "junction_list.csv", "log_title_line.csv", "multi_month_log.csv", "junctions_2col.csv", "long_paren.csv"):
    at.file_uploader[0].set_value((name, open(fx[name], "rb").read(), "application/octet-stream"))
    at.run()
    assert not at.exception, (name, [e.value for e in at.exception])
    res[name] = {"markdown": [m.value for m in at.markdown if m.value.startswith(("**This", "**Not sure", "**As",
                                                                                    "**Model readiness"))],
                 "warnings": [w.value for w in at.warning]}
print("PAGE_JSON " + json.dumps(res))
"""


@check("task16", quick=False)
def get_started_page_runs():
    """Headless AppTest of pages/1_Get_started.py (PYTHONPATH set): every software choice of step 1 renders its own
    section of docs/how_to_export_your_model.md (12 choices; the GIS and CAD ones marked 'coming'); step 2 has one
    download button per template (the four CSVs and the Excel grab log, written on the fly by openpyxl); every template
    parses (the grab logs with load_log and the tap map template, the tap map as load_log's tap map, the plant log with
    load_plant, the workbook through the same reader: part 'templates' of the committed record); step 3, 'Check a
    file', gives mains.gpkg its sniff result and advice and net3_backup.inp.bak its sniff result and readiness report,
    reads lims_censored.csv as a grab log (its censored rows listed), task16_tap_map.csv as a tap map, and asks the
    question for ambiguous_dates.csv; since the review of task 16, junction_list.csv (a junction column and no tap
    column, which crashed the page) is said to be a junction list, not a tap map, log_title_line.csv is read past its
    title, and multi_month_log.csv gives the main page's 3-month window; since its final verification,
    junctions_2col.csv ('Node ID,Elevation') is a junction list too, never a tap map whose junction column is
    'Elevation', and long_paren.csv ('Chlorine (Free)', 'Chlorine (Total)', pH rows) gives its 3 free chlorine
    readings only, with no 'does not say free or total' warning.  No grid is built (outputs/cache unchanged),
    and nothing is left in outputs/private/check_a_file or outputs/app_uploads.  The task-19 form parser does not exist
    yet."""
    from . import intake_exp as E
    from . import fixtures as FX
    t = E.part_templates()
    for f, v in t.items():
        assert v.get("rows", v.get("taps", v.get("months", 0))) > 0, (f, v)
    assert t["grab_log_template.csv"]["censored"] == 1 and t["grab_log_template.xlsx"]["same_as_csv"], t
    _same_as_16(t, "templates")
    sections = [ln[3:].strip() for ln in open(os.path.join(REPO, "docs", "how_to_export_your_model.md"),
                                              encoding="utf-8") if ln.startswith("## ")]
    cache0 = _cache_snapshot()
    check_dir = os.path.join(REPO, "outputs", "private", "check_a_file")
    up_dir = os.path.join(REPO, "outputs", "app_uploads")
    ups0 = sorted(os.listdir(up_dir)) if os.path.isdir(up_dir) else []
    with tempfile.TemporaryDirectory(prefix="rm_t16_page_") as d:
        paths = FX.build_all(d, ["mains.gpkg", "net3_backup.inp.bak", "lims_censored.csv", "task16_tap_map.csv",
                                 "ambiguous_dates.csv", "junction_list.csv", "log_title_line.csv",
                                 "multi_month_log.csv", "junctions_2col.csv", "long_paren.csv"])
        r = subprocess.run([PY, "-c", PAGE_SCRIPT, GET_STARTED_PAGE, json.dumps(paths)], cwd=os.getcwd(),
                           capture_output=True, text=True, env={**os.environ, "PYTHONPATH": REPO})
    line = [x for x in r.stdout.splitlines() if x.startswith("PAGE_JSON ")]
    if r.returncode != 0 or not line:
        raise RuntimeError(f"AppTest failed (exit {r.returncode}): {r.stderr[-3000:]}")
    res = json.loads(line[-1][len("PAGE_JSON "):])
    choices = res["choices"]
    assert list(choices) == [s for s in sections if s != "Your grab samples"] and len(choices) == 12, list(choices)
    assert len(set(choices.values())) == len(choices), "two choices show the same steps"
    coming = [c for c in choices if "coming" in c]
    assert len(coming) == 5 and all(any(w in c for w in ("AutoCAD", "MicroStation", "ArcGIS", "QGIS", "Google")) for c in coming)
    assert len(res["downloads"]) == 5, res["downloads"]
    xlsx = [f for f in res["files"] if f and f.endswith(".xlsx")]
    assert xlsx and sorted(f for f in res["files"] if f) == sorted(
        ["grab_log_template.csv", "grab_log_total_template.csv", "tap_map_by_id_template.csv", "plant_log_template.csv",
         "grab_log_template.xlsx"]), res["files"]
    g, b = res["mains.gpkg"], res["net3_backup.inp.bak"]
    assert g["markdown"][0].startswith("**This file is a GIS layer (GeoPackage)") and any("EPANET export" in w for w in g["warnings"]), g
    assert b["markdown"][0].startswith("**This file is a backup copy of an EPANET model") and any(
        m.startswith("**Model readiness: net3_backup.inp.bak") and "| Test run |" in m for m in b["markdown"]), b
    lc = res["lims_censored.csv"]
    assert any(m.startswith("**As a grab log:** 3 readings") for m in lc["markdown"]) and any(
        "censored" in w for w in lc["warnings"]), lc
    assert any(m.startswith("**As a tap map:** 20 taps") for m in res["task16_tap_map.csv"]["markdown"])
    assert any("Which is it?" in w for w in res["ambiguous_dates.csv"]["warnings"])
    # the review of task 16: a junction list crashed the page; a title line before the header; the main page's window
    assert any(m.startswith("**As a table:** it has a junction column ('Node ID')") for m in
               res["junction_list.csv"]["markdown"]), res["junction_list.csv"]
    assert any(m.startswith("**As a grab log:** 4 readings") and "skipped" in m for m in
               res["log_title_line.csv"]["markdown"]), res["log_title_line.csv"]
    mm = _record_16("app_inputs")["multi_month"]
    assert any(f"The main page uses the last 3 months that hold a reading: {mm['rows_used']} readings" in m
               for m in res["multi_month_log.csv"]["markdown"]), res["multi_month_log.csv"]
    # the final verification of task 16
    j2 = res["junctions_2col.csv"]
    assert any(m.startswith("**As a table:** it has a junction column ('Node ID')") for m in j2["markdown"]), j2
    assert not any(m.startswith("**As a tap map") for m in j2["markdown"]), j2
    lp = res["long_paren.csv"]
    assert any(m.startswith(f"**As a grab log:** {len(FX.LONG_SAMPLES)} readings") for m in lp["markdown"]), lp
    assert not any("does not say free or total" in w for w in lp["warnings"]), lp["warnings"]
    cache1 = _cache_snapshot()
    assert cache1 == cache0, sorted(set(cache1) ^ set(cache0))
    assert not os.path.exists(check_dir) or not os.listdir(check_dir), os.listdir(check_dir)
    assert (sorted(os.listdir(up_dir)) if os.path.isdir(up_dir) else []) == ups0
    return {"software_choices": len(choices), "coming": len(coming), "template_downloads": len(res["downloads"]),
            "templates_parse": sorted(t), "check_a_file": ["mains.gpkg: advice", "net3_backup.inp.bak: readiness report",
                                                           "lims_censored.csv: as a grab log", "task16_tap_map.csv: as a tap map",
                                                           "ambiguous_dates.csv: the question",
                                                           "junction_list.csv: a junction list, not a tap map",
                                                           "log_title_line.csv: as a grab log, its title skipped",
                                                           "multi_month_log.csv: the main page's window",
                                                           "junctions_2col.csv: a junction list, not a tap map",
                                                           "long_paren.csv: its free chlorine rows only"],
            "grids_built": 0, "files_left": 0}


APP_ORDER_SCRIPT = r"""
import json, sys, warnings
warnings.filterwarnings("ignore")
from streamlit.testing.v1 import AppTest
import residualmap.chloramine as CA
app, inp = sys.argv[1], sys.argv[2]
def _no_msx(*a, **k):      # the demo truth's EPANET-MSX run is stood in by the reduced truth: no widget depends on it
    raise RuntimeError("EPANET-MSX not run by app_widget_order_unchanged")
CA.truth_in_subprocess = _no_msx
def fresh():
    return AppTest.from_file(app, default_timeout=900).run()
def pick(at, label, value):
    if label == "selectbox":
        at.selectbox[0].set_value([o for o in at.selectbox[0].options if o.startswith(value)][0])
    else:
        [r for r in at.radio if r.label == label][0].set_value(value)
    return at.run()
def toggle(at, prefix, value):
    [t for t in at.toggle if t.label.startswith(prefix)][0].set_value(value)
    return at.run()
def order(at):
    assert not at.exception, [e.value for e in at.exception]
    return {"selectbox0": at.selectbox[0].label if at.selectbox else None,
            "uploaders": [u.label for u in at.file_uploader], "example_shown": any(s.label == "Example" for s in at.selectbox)}
res = {}
at = fresh(); res["t14_default"] = order(at)
res["t14_hourly"] = order(pick(at, "Show", "03:00"))
res["t14_net2"] = order(pick(fresh(), "selectbox", "Net2"))
res["t14_chloramine_net3"] = order(pick(fresh(), "Disinfectant", "Chloramine (total chlorine)"))
at = pick(fresh(), "Disinfectant", "Chloramine (total chlorine)"); res["t14_ky4_chloramine"] = order(pick(at, "selectbox", "ky4"))
at = pick(fresh(), "selectbox", "Net2"); res["t14_msx_fallback"] = order(pick(at, "Disinfectant", "Chloramine (total chlorine)"))
at = pick(fresh(), "Model", "Upload my .inp")
at.file_uploader[0].set_value(("Net3.inp", open(inp, "rb").read(), "application/octet-stream")); at.run()
res["t14_upload"] = order(at)
at = pick(fresh(), "selectbox", "Net2"); res["ca_net2"] = order(pick(at, "Disinfectant", "Chloramine (total chlorine)"))
res["age_net3"] = order(fresh())
res["age_net2"] = order(pick(fresh(), "selectbox", "Net2"))
res["age_demo_off"] = order(toggle(fresh(), "Demo", False))
at = pick(fresh(), "Model", "Upload my .inp")                              # new: upload, demo off
at.file_uploader[0].set_value(("Net3.inp", open(inp, "rb").read(), "application/octet-stream")); at.run()
res["new_upload_demo_off"] = order(toggle(at, "Demo", False))
print("APP_ORDER_JSON " + json.dumps(res))
"""


@check("task16", quick=False)
def app_widget_order_unchanged():
    """The widget order the older app checks rely on (critique 18), replayed by headless AppTest scenario by scenario:
    app_task14_runs (default, hourly view, Net2, chloramine on Net3, ky4 with chloramine, the EPANET-MSX fallback on
    Net2, the Net3.inp upload), app_chloramine_mode (Net2 with chloramine) and app_water_age_panels (Net3, Net2, the
    demo off).  In each, selectbox[0] is 'Example' wherever the example picker is shown, and file_uploader[0] is
    'EPANET .inp file' in every upload scenario; the example network with the demo off has no model uploader, so its
    first uploader is the grab log, by design.  One new scenario (upload, demo off) shows the model, grab-log and
    tap-map uploaders in that order.  The demo truth's EPANET-MSX run is stood in by the reduced truth (raising in its
    place), which changes no widget; the upload's copy and grid are deleted afterwards."""
    import hashlib
    from .simulate import LIB
    inp = os.path.join(LIB, "Net3.inp")
    up_name = f"{hashlib.sha1(open(inp, 'rb').read()).hexdigest()[:8]}_Net3.inp"    # app.py's name for the upload
    res, _ = _app_run(APP_ORDER_SCRIPT, [os.path.join(REPO, "app.py"), inp],
                      lambda f: f == up_name or f"_{up_name}" in f, "APP_ORDER_JSON")
    for k, v in res.items():
        if v["example_shown"]:
            assert v["selectbox0"] == "Example", (k, v)
    for k in ("t14_upload", "new_upload_demo_off"):
        assert res[k]["uploaders"][:1] == ["EPANET .inp file"] and not res[k]["example_shown"], (k, res[k])
    assert res["new_upload_demo_off"]["uploaders"] == ["EPANET .inp file", "Grab log (CSV or Excel)", "Tap map (CSV)"]
    assert res["age_demo_off"]["uploaders"] == ["Grab log (CSV or Excel)", "Tap map (CSV)"], res["age_demo_off"]
    shown = [k for k, v in res.items() if v["example_shown"]]
    assert len(shown) == 10, shown
    return {"scenarios": len(res), "example_picker_first": len(shown), "upload_scenarios_model_uploader_first": 2,
            "new_order": res["new_upload_demo_off"]["uploaders"]}


@check("task16", quick=False)
def pilot_clock_and_log():
    """The pilot's real-data path with task 16's reader, from fresh working directories: (a) lims_one_datetime_us.csv
    on Net3 in chloramine mode with --date-format MDY runs (its model hours are its times of day, Net3 starting at
    midnight) and, without --date-format, stops with the question; (b) clock_log_clock6.csv on net3_clock6.inp (a
    file, so it gets a private copy and a fresh grid) and (c) clock_log_net3.csv, the same readings written 6 h
    earlier, on Net3 by name: the hours of (b)'s predictions are its times of day minus 6 h, and (b) and (c) write
    byte-identical predictions and validation files, because net3_clock6's grid is Net3's bit for bit, and (b) says
    that its outputs' hours are model hours (the review of task 16).
    synthetic_pilot_net3 (unchanged) covers the synthetic path, whose hours are model hours already."""
    from . import fixtures as FX
    out = {}
    with _sample_fixtures() as paths, tempfile.TemporaryDirectory(prefix="rm_t16_pilot_") as d:
        def pilot(sub: str, link: bool, *args) -> subprocess.CompletedProcess:
            wd = os.path.join(d, sub)
            os.makedirs(os.path.join(wd, "outputs"))
            if link:
                os.symlink(CACHE, os.path.join(wd, "outputs", "cache"))      # a committed grid is read, not rebuilt
            return subprocess.run([PY, "-m", "residualmap.pilot", *args], cwd=wd, capture_output=True, text=True,
                                  env={**os.environ, "PYTHONPATH": REPO})
        lims = ["--inp", "Net3", "--samples", paths["lims_one_datetime_us.csv"], "--taps", paths["task16_tap_map.csv"],
                "--disinfectant", "chloramine", "--dose", "2.0"]
        ra = pilot("a", True, *lims, "--date-format", "MDY")
        assert ra.returncode == 0, ra.stderr[-2000:]
        pa = pd.read_csv(os.path.join(d, "a", "outputs", "private", "pilot", "predictions_Net3.csv"), dtype={"tap_id": str})
        want_h = {(f"2026-{mo:02d}-{dd:02d}", t): h for mo, dd, h, mi, t, v in FX.LIMS_US_ROWS}
        assert all(want_h[(r.date, r.tap_id)] == r.hour for r in pa.itertuples()) and len(pa) > 0, pa.head()
        rq = pilot("q", True, *lims)
        assert rq.returncode != 0 and "Which is it?" in (rq.stderr + rq.stdout), (rq.returncode, rq.stderr[-500:])
        rb = pilot("b", False, "--inp", paths["net3_clock6.inp"], "--samples", paths["clock_log_clock6.csv"],
                   "--taps", paths["task16_tap_map.csv"])
        assert rb.returncode == 0, rb.stderr[-2000:]
        assert "every hour in the pilot's outputs (predictions and validation files) is a model hour" in rb.stdout
        rc = pilot("c", True, "--inp", "Net3", "--samples", paths["clock_log_net3.csv"], "--taps",
                   paths["task16_tap_map.csv"])
        assert rc.returncode == 0, rc.stderr[-2000:]
        pb = pd.read_csv(os.path.join(d, "b", "outputs", "private", "pilot", "predictions_net3_clock6.csv"), dtype={"tap_id": str})
        model_h = {(dt_, t): h for dt_, h, mi, t, v in FX.CLOCK_ROWS}
        assert all(model_h[(r.date, r.tap_id)] == r.hour for r in pb.itertuples()) and len(pb) > 0, pb.head()
        same = {}
        for kind in ("predictions", "validation"):
            for ext in ("csv",) if kind == "predictions" else ("csv", "json"):
                fb = os.path.join(d, "b", "outputs", "private", "pilot", f"{kind}_net3_clock6.{ext}")
                fc = os.path.join(d, "c", "outputs", "private", "pilot", f"{kind}_Net3.{ext}")
                same[f"{kind}.{ext}"] = open(fb, "rb").read() == open(fc, "rb").read()
        assert all(same.values()), same
        out = {"lims_us_with_date_format_rows": len(pa), "lims_us_without_date_format_stops_with_the_question": True,
               "clock6_predictions_rows": len(pb), "clock6_hours_are_clock_minus_6": True,
               "clock6_and_net3_byte_identical": sorted(same), "clock6_says_its_outputs_are_model_hours": True}
    return out


# ----------------------------------------------------------------------------- task 16b: tap mapping by coordinates
# Task 16b stopped at its pre-registered bar (ky4, 0.967 of 2,000 taps against 0.98).  Five of the plan's six checks
# are here; pilot_suggested_tap_map is not, because the pilot's --taps by coordinates and --accept-snapped were not
# built.  Each check reruns its part of outputs/intake/tap_snapping.json and compares; tap_coordinates_review_cases
# holds the review's cases.
INTAKE_RECORD_16B = os.path.join(REPO, "outputs", "intake", "tap_snapping.json")


def _record_16b(part: str) -> dict:
    """One part of the committed task-16b record (python -m residualmap.intake_exp snapping)."""
    with open(INTAKE_RECORD_16B) as fh:
        return json.load(fh)[part]


def _same_as_16b(got: dict, part: str) -> dict:
    want = _record_16b(part)
    diff = _json_diff(json.loads(json.dumps(got)), want)
    assert not diff, (part, diff[:10])
    return want


@contextlib.contextmanager
def _tap_fixtures():
    """The task-16b fixtures built into a temporary folder."""
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t16b_") as d:
        yield FX.build_all(os.path.join(d, "fx"), FX.TAP_FIXTURES)


@check("task16b")
def georef_enu_affine():
    """Taps in latitude and longitude through the local east-north plane and the affine fit, against the ellipsoidal
    Lambert conformal conic truth of fixtures.py (which matches Snyder's worked example to 0.1 m), over 20 and 50 km
    zones.  Bar (plan task 16b, met): at most 1 m maximum error with 3 control taps, a fit marked unchecked; with 6
    control taps each one's leave-one-out residual is printed; the direct latitude and longitude affine, kept only in
    the record, is more than 10 times worse; a control tap moved 25 m is flagged, with the largest residual; nearly
    collinear control taps are refused.  The review's cases: a zone across the antimeridian equals the same zone at
    100 W, a narrow valley (8 km by 300 m, 3 m of GPS error on its control taps) is accepted with a warning, control
    taps that fail only the triangle test are refused naming that test alone, and 4,000 control taps are measured in
    under 2 s.  The final verification's cases: the plane's tangent point at the taps' median (the bar passes with it
    and as scored, with it at their mean); one typo tap in a small tap map leaves every test tap within 1 m (with the
    plane at the mean, as first built, it did not); control taps along one avenue for a 4 km square zone flag every
    test tap placed more than 10 m off; control taps at one place or two are refused in those words; the leave-one-out
    residuals in closed form equal one refit per control tap, and 4,000 control taps take under 2 s through
    georeference.  Equals the committed record (part 'georef')."""
    from . import fixtures as FX
    from . import intake_exp as E
    from .importers import control_spread, georeference
    with _tap_fixtures() as paths:
        got = E.part_georef(paths)
    bar = got["bar"]
    assert bar["passed"] and bar["at_most_m"] == 1.0 and bar["result_m"] <= 1.0, bar
    assert bar["result_m_as_scored"] <= 1.0, bar
    assert got["moved_control"]["closed_form_against_a_refit_per_control_tap_max_abs_diff_m"] <= 1e-6, got["moved_control"]
    for km in ("20km", "50km"):
        z = got[km]
        assert not z["3_controls"]["checked"] and z["3_controls"]["loo_m"] == {}, z["3_controls"]
        assert z["6_controls"]["checked"] and len(z["6_controls"]["loo_m"]) == 6, z["6_controls"]
        assert z["3_controls"]["direct_latlon_affine_max_error_m"] > 10 * z["3_controls"]["max_error_m"], z
    mc = got["moved_control"]
    assert mc["moved"] in mc["flagged"] and max(mc["loo_m"], key=mc["loo_m"].get) == mc["moved"], mc
    assert got["collinear"]["refused"], got["collinear"]
    rv = got["review"]
    assert rv["lcc_against_snyder"]["max_abs_diff_m"] <= 0.1, rv["lcc_against_snyder"]
    am = rv["antimeridian"]
    assert am["lon_min"] < -179 and am["lon_max"] > 179, am
    for k in ("3_controls", "6_controls"):
        assert am[k]["max_error_m"] == am[k]["zone_at_100w_max_error_m"] and not am[k]["flagged"], am
        assert am[k]["arithmetic_mean_longitude_max_error_m"] > 1.0, am          # as first built
    nz = rv["narrow_zone"]
    assert nz["accepted"] and nz["relative_tests_alone_refuse"] and nz["warned"] and nz["max_error_m"] <= 5.0, nz
    tr = rv["triangle_only"]
    assert tr["refused"] and tr["names_the_triangle_test_only"], tr
    fin = rv["final"]
    for zone, z in fin["typo_in_a_small_map"].items():
        for label, v in z.items():
            assert v["plane_at_the_median"]["max_error_m"] <= 1.0, (zone, label, v)
            assert v["plane_at_the_median"]["loo_max_m"] < 10.0, (zone, label, v)
    assert max(v["plane_at_the_mean_as_first_built"]["max_error_m"]
               for z in fin["typo_in_a_small_map"].values() for v in z.values()) > 10.0, fin["typo_in_a_small_map"]
    st = fin["strip_controls"]
    assert st["accepted"] and st["narrow_strip"] and st["over_30m"] > 0, st
    assert st["over_10m_not_flagged"] == 0 and st["over_30m_not_flagged"] == 0, st
    for k, words in (("one_place", "at the same place"), ("two_places", "only two places"),
                     ("within_4m", "about the same place")):
        assert fin[k]["refused"] and words in fin[k]["message"] and "inf" not in fin[k]["message"], fin[k]
    t0 = time.time()
    control_spread(np.random.default_rng(1).uniform(-5000, 5000, size=(4000, 2)))
    spread_s = time.time() - t0
    assert spread_s < 2.0, spread_s
    # 4,000 control taps through georeference (positions laid out in degrees: only the time is measured)
    E4 = np.random.default_rng(1).uniform(-5000, 5000, size=(4001, 2))
    lat4 = FX.LCC_CENTRE[0] + E4[:, 1] / 111_320.0
    lon4 = FX.LCC_CENTRE[1] + E4[:, 0] / (111_320.0 * math.cos(math.radians(FX.LCC_CENTRE[0])))
    t0 = time.time()
    _, info4 = georeference(lat4, lon4, list(range(4000)), E4[:4000], 1.0)
    georef_s = time.time() - t0
    assert georef_s < 2.0 and len(info4["loo"]) == 4000, georef_s
    _same_as_16b(got, "georef")
    return {"bar_max_error_3_controls_m": bar["result_m"], "bar_at_most_m": bar["at_most_m"], "bar_passed": True,
            "max_error_m": {km: {k: got[km][k]["max_error_m"] for k in ("3_controls", "6_controls")} for km in ("20km", "50km")},
            "direct_latlon_affine_max_error_m": {km: got[km]["3_controls"]["direct_latlon_affine_max_error_m"]
                                                  for km in ("20km", "50km")},
            "moved_control_flagged": mc["flagged"], "collinear_refused": True,
            "snyder_max_abs_diff_m": rv["lcc_against_snyder"]["max_abs_diff_m"],
            "antimeridian_max_error_m": am["3_controls"]["max_error_m"], "narrow_zone_max_error_m": nz["max_error_m"],
            "triangle_only_refused_naming_that_test": True,
            "bar_result_m_as_scored": bar["result_m_as_scored"],
            "typo_small_map_max_error_m": max(v["plane_at_the_median"]["max_error_m"]
                                              for z in fin["typo_in_a_small_map"].values() for v in z.values()),
            "strip_over_10m_not_flagged": st["over_10m_not_flagged"],
            "_volatile": {"control_spread_4000_s": round(spread_s, 3), "georeference_4000_s": round(georef_s, 3)}}


@check("task16b", quick=False)
def tap_snapping_ky4():
    """The pre-registered scored run of task 16b, rerun: 2,000 taps uniform by length on ky4's pipes under 16 in
    (feet, with vertices), at 5, 10 and 20 m of radial RMS GPS noise, each put on the nearest pipe under 16 in and
    given that pipe's downstream junction by the simulated daytime flow.  Equals the committed record (part 'ky4'),
    and the measurement of its misses made after the stop (part 'ky4_after_the_stop', not a rescoring).
    Not asserted: the bar (at least 0.98 of taps with the right junction at 5 m).  Its scored run missed it, 0.967
    counted strictly and 0.9775 counting a right flagged second choice (open decision 11), so the task stopped; the
    miss stands, recorded in the journal, and no passing run of this check counts as a pass of that bar."""
    from . import intake_exp as E
    got = E.part_snap_ky4()
    after = E.part_ky4_after_the_stop()
    _same_as_16b(got, "ky4")
    _same_as_16b(after, "ky4_after_the_stop")
    bar, five = got["bar"], got["by_noise"]["5"]
    assert (bar["network"], bar["noise_m"], bar["share_at_least"]) == ("ky4", 5.0, 0.98), bar
    assert bar["counting"].startswith("strict") and bar["result"] == five["right_junction"], bar
    assert bar["passed"] == (five["right_junction"] >= 0.98), bar
    assert after["by_noise"]["5"]["wrong_junction"] == round((1 - five["right_junction"]) * got["taps"]), after
    assert all(v["ambiguous_rows_list_two_candidates"] for v in got["by_noise"].values()), got["by_noise"]
    return {"bar": "pre-registered: at least 0.98 of 2,000 taps get the right junction at 5 m; missed in its scored "
                   "run, so task 16b stopped; not asserted here",
            "bar_passed": bar["passed"],
            "right_junction": {k: v["right_junction"] for k, v in got["by_noise"].items()},
            "right_junction_or_flagged_second": {k: v["right_junction_or_flagged_second"] for k, v in got["by_noise"].items()},
            "wrong_junction_unflagged": {k: v["wrong_junction_unflagged"] for k, v in after["by_noise"].items()}}


@check("task16b")
def tap_snapping_urban_grid():
    """The synthetic urban grid (invented, metres, to scale): 2,000 taps on its pipes under 16 in and 100 on its 24 in
    transmission main tagged on_main yes, at 5, 10 and 20 m: no untagged tap is put on the transmission main, the 4 in
    line 20 m beside it included, and every ambiguous row lists two candidates.  The review's variants: with the main
    a metric DN400 (400 mm) or a 16 in main by its internal diameter (388.6 mm), no untagged tap is put on it either
    (with the cut at 16.0 in on the model's diameter, as first built, some were).  The final verification's: a 16 in
    plastic main by its internal diameter under the 15.0 in cut (C905 DR14, 378.8 mm; PE100 SDR11, 327.4 mm) takes
    untagged taps, and every one of them is flagged as on a pipe under the cut.  Equals the committed record (part
    'urban_grid')."""
    from . import intake_exp as E
    with _tap_fixtures() as paths:
        got = E.part_snap_urban(paths)
    for k, v in got["by_noise"].items():
        assert v["untagged_taps_put_on_the_transmission_main"] == 0 and v["four_inch_taps_put_on_the_transmission_main"] == 0, (k, v)
        assert all(v[m]["ambiguous_rows_list_two_candidates"] for m in ("small_mains", "dual_mains", "four_inch_line")), (k, v)
    first_built = under_cut = 0
    for name, z in got["transmission_main_variants"].items():
        for k, v in z["by_noise"].items():
            c15 = v["cut_15in"]
            if z["main_in_model"] >= 15.0:
                assert c15["untagged_taps_put_on_the_transmission_main"] == 0, (name, k, v)
                first_built += v["cut_16in_as_first_built"]["untagged_taps_put_on_the_transmission_main"]
            else:                       # under the cut: the taps put on it are flagged, every one
                n = c15["untagged_taps_put_on_the_transmission_main"]
                assert c15["untagged_taps_on_the_main_with_no_flag_naming_a_main"] == 0, (name, k, v)
                assert c15["untagged_taps_on_the_main_flagged_under_the_cut"] == n, (name, k, v)
                under_cut += n
    assert first_built > 0 and under_cut > 0, got["transmission_main_variants"]
    _same_as_16b(got, "urban_grid")
    return {"right_junction_small_mains": {k: v["small_mains"]["right_junction"] for k, v in got["by_noise"].items()},
            "right_junction_dual_mains": {k: v["dual_mains"]["right_junction"] for k, v in got["by_noise"].items()},
            "untagged_taps_on_the_transmission_main": 0,
            "variants_untagged_taps_on_the_main": {n: {k: (v["cut_15in"]["untagged_taps_put_on_the_transmission_main"],
                                                          v["cut_16in_as_first_built"]["untagged_taps_put_on_the_transmission_main"])
                                                      for k, v in z["by_noise"].items()}
                                                  for n, z in got["transmission_main_variants"].items()}}


@check("task16b", quick=False)
def tap_snapping_direction_invariant():
    """ky4_half_reversed.inp (every second pipe of ky4 drawn the other way, its vertices reversed, so the drawn
    direction and the sign of the simulated flow both flip there) gives the same 2,000 assignments at 5 m as ky4, and
    the flow sign flips on every flipped pipe that carries flow.  Equals the committed record (part 'direction')."""
    from . import intake_exp as E
    with _tap_fixtures() as paths:
        got = E.part_snap_direction(paths)
    assert got["assignments_differ"] == 0 and got["junction_rules_differ"] == 0 and got["same_origin"], got
    assert got["flow_sign_flipped_on"] == got["flipped_with_flow"] > 0, got
    _same_as_16b(got, "direction")
    return {k: got[k] for k in ("taps", "pipes_flipped", "flipped_with_flow", "flow_sign_flipped_on", "assignments_differ")}


@check("task16b")
def tap_snapping_refuses_schematic():
    """Net3 (a schematic drawing) with taps in x and y and in latitude and longitude, and Net3 without coordinates,
    are refused with 'Use junction IDs instead'; no accuracy is computed on Net3.  Equals the committed record (part
    'schematic')."""
    from . import intake_exp as E
    with _tap_fixtures() as paths:
        got = E.part_snap_schematic(paths)
    assert len(got) == 3 and all(v["refused"] and v["says_use_junction_ids"] for v in got.values()), got
    _same_as_16b(got, "schematic")
    return {"refused_with_use_junction_ids": sorted(got)}


@check("task16b")
def tap_coordinates_review_cases():
    """The review of task 16b, each case a fixture (part 'review'), and the coordinates template's suggestion read
    back through load_log as a tap map by ID (part 'template'; the template is a fixture only, not shipped in
    docs/templates).  Which pair of columns places the taps (latitude and longitude left blank, both pairs and no
    control taps: x and y; both pairs with control taps: latitude and longitude, or x and y when chosen, the same
    junctions); '1,197' in x or y not settled by the column is not used, with its reason, and settled by plain values
    of 1,000 or more it is 1197; a latitude with the letter W is not used, with its reason; swapped latitude and
    longitude columns are refused as swapped; a tap listed again elsewhere warns, and with two junctions is refused; a
    position at x = 1e200 is not used and raises no overflow warning; a tap 1 m from the 24 in main is flagged; a
    mapping that names a missing column is refused; a junction column known only through a mapping makes the map 'id',
    and unmapped it is named in a warning; a model drawn in degrees is refused (no traceback), one drawn in miles gets
    1609.344 m per unit and the same junctions as in metres, two nodes drawn at 0,0 leave their 8 pipes out and no
    tap on them, and a model drawn in part warns.  The final verification's cases: '1.197' in a file of decimal
    commas is not used, with its reason, and settled by '1.503,5' it is 1197; longitudes written from 0 to 360 give
    the template's junctions; control taps along one street flag the taps far across it, and the control tap the fit
    cannot be made without is not flagged; one node at -3.4e38 or at nan leaves every tap's suggestion as on the clean
    drawing, with a warning naming it; drawings in degrees with latitude as x, or with one node at 9999,9999, are
    refused.  Equals the committed record."""
    from . import fixtures as FX
    from . import intake_exp as E
    with _tap_fixtures() as paths:
        got = E.part_review_16b(paths)
        tpl = E.part_template(paths)
    c = lambda n: got[n]
    assert c("taps_both_pairs_latlon_blank.csv")["kind"] == "xy" and len(c("taps_both_pairs_latlon_blank.csv")["rows"]) == 3
    assert c("taps_both_pairs_no_controls.csv")["kind"] == "xy"
    bp = c("taps_both_pairs_controls.csv")
    assert bp["guessed"]["kind"] == "latlon" and bp["x_and_y_chosen"]["kind"] == "xy" and bp["same_junctions"], bp
    ta = c("taps_thousands_ambiguous.csv")
    assert [d["tap_id"] for d in ta["not_used"]] == ["T2", "T3"] and all("reads two ways" in d["reason"] for d in ta["not_used"]), ta
    ts = {r["tap_id"]: r["junction_id"] for r in c("taps_thousands_settled.csv")["rows"]}
    want_t2 = {r["tap_id"]: r["junction_id"] for r in c("taps_both_pairs_latlon_blank.csv")["rows"]}["T2"]
    assert ts["T2"] == want_t2, (ts, want_t2)
    ls = c("taps_letters_swapped.csv")
    assert any(d["tap_id"] == "T1" and "has the letter W" in d["reason"] for d in ls["not_used"]), ls["not_used"]
    assert c("taps_latlon_swapped.csv")["refused"] and "look swapped" in c("taps_latlon_swapped.csv")["message"]
    dp = c("taps_duplicate_positions.csv")
    assert any("listed again with another position" in w for w in dp["warnings"]), dp
    assert c("taps_duplicate_junctions.csv")["refused"] and "two different junctions" in c("taps_duplicate_junctions.csv")["message"]
    ab = c("taps_absurd.csv")
    assert ab["runtime_warnings"] == 0 and any(d["tap_id"] == "T1" and "not near the model" in d["reason"]
                                                for d in ab["not_used"]), ab
    bm = {r["tap_id"]: r for r in c("taps_beside_main.csv")["rows"]}
    assert "a main of 16 in or more is 1.0 m away" in bm["T1"]["flags"] and bm["T2"]["pipe"].startswith("TM"), bm
    assert c("taps_mapping_missing.csv")["refused"] and "no column 'node ref'" in c("taps_mapping_missing.csv")["message"]
    jh = c("taps_junction_header_unknown.csv")
    assert jh["kind_guessed"] == "coordinates" and jh["kind_with_mapping"] == "id", jh
    assert jh["read_by_id_with_mapping"] == {"T1": "J0101", "T2": "J0909"}, jh
    assert any("'EPANET_Name' holds junction IDs" in w for w in jh["by_coordinates"]["warnings"]), jh
    m = got["models"]
    deg = m["urban_grid_degrees.inp"]
    assert deg["refused"] and "degrees" in deg["message"] and deg["suggest_tap_map"]["refused"], deg
    mi = m["urban_grid_miles.inp"]
    assert abs(mi["metres_per_unit"] / 1609.344 - 1) < 1e-4 and mi["taps_5m"]["same_junction_as_the_metres_drawing"] == 1.0, mi
    assert abs(mi["metres_per_unit_as_first_built"] / 1609.344 - 1) > 0.01, mi                 # the rounded median
    ph = m["urban_grid_placeholders.inp"]
    assert len(ph["left_out"]) == 8 and ph["taps_5m"]["put_on_a_pipe_left_out"] == 0 and ph["warnings"], ph
    assert any(w.startswith("Only ") for w in m["urban_grid_part_drawn.inp"]["warnings"]), m["urban_grid_part_drawn.inp"]
    assert tpl["checked"] and tpl["accepted_rows_through_load_log"] == 32 and not tpl["not_used"], tpl
    fin = got["final"]
    dt = fin["taps_dot_thousands.csv"]
    assert [r["tap_id"] for r in dt["rows"]] == ["T1"] and [d["tap_id"] for d in dt["not_used"]] == ["T2", "T3"], dt
    assert all("reads two ways ('1." in d["reason"] and "decimal point" in d["reason"] for d in dt["not_used"]), dt
    ds = {r["tap_id"]: r["junction_id"] for r in fin["taps_dot_thousands_settled.csv"]["rows"]}
    clean = {t: j for t, j, _ in fin["models"]["clean drawing"]["rows"]}
    assert ds["T2"] == clean["T2"] and ds["T3"] == clean["T3"], (ds, clean)
    lz = fin["taps_lon_0_to_360.csv"]
    assert lz["same_junctions_as_the_template"] and any("written from 0 to 360" in x for x in lz["notes"]), lz
    sc = {r["tap_id"]: r["flags"] for r in fin["taps_controls_along_one_street.csv"]["rows"]}
    assert all("the fit is uncertain here" in sc[t] for t in ("T2", "T3")) and "uncertain" not in sc["T1"], sc
    assert "cannot be measured" in sc["C4"] and not any("over 10 m" in v for v in sc.values()), sc
    for mname in FX.URBAN_SENTINELS:
        mm = fin["models"][mname]
        assert mm["same_rows_as_the_clean_drawing"] and any("J0505" in w for w in mm["warnings"]), mm
    for mname in ("urban_grid_degrees_lat_as_x.inp", "urban_grid_degrees_placeholder.inp"):
        assert fin["models"][mname]["refused"] and "degrees" in fin["models"][mname]["message"], fin["models"][mname]
    _same_as_16b(got, "review")
    _same_as_16b(tpl, "template")
    return {"cases": len(got) - 2, "models": sorted(m), "template_rows_through_load_log": tpl["accepted_rows_through_load_log"],
            "miles_metres_per_unit": mi["metres_per_unit"], "placeholder_pipes_left_out": len(ph["left_out"]),
            "final_cases": len(fin) - 2, "final_models": sorted(fin["models"])}


# ----------------------------------------------------------------------------- runner
def _counts(rs: list) -> dict:
    return {"all_pass": all(r["pass"] is not False for r in rs),
            "n_pass": sum(r["pass"] is True for r in rs), "n_fail": sum(r["pass"] is False for r in rs),
            "n_skipped": sum(r["pass"] is None for r in rs)}


def settle_numbers_entry(report: dict, run_log: dict | None = None, path: str = REPORT) -> dict:
    """The numbers scan (numbers_trace_to_outputs) reads the checks report as one of its sources, every entry but its
    own, so in a run whose other entries changed (a new check, a new detail) it scanned the previous run's report, and
    the next run on the same tree wrote a different entry (task 15b's second review).  This writes the report to
    path, with its pass counts brought up to date first (they are numbers in the scan's pool too), reruns the scan on
    that file and puts the result in the entry, and repeats while the entry changes: the scan leaves its own entry
    out, so the entry settles on the first pass, unless its pass or fail changed the counts.  Returns what it did,
    for the run log."""
    from . import numbers as N
    entry = next((r for r in report["checks"] if r["name"] == "numbers_trace_to_outputs"), None)
    if entry is None or entry["pass"] is None:
        return {"rerun": False}
    changed, passes = False, 0
    keep = N.PATH_OVERRIDE.get(N.SELF_REPORT)
    N.PATH_OVERRIDE[N.SELF_REPORT] = path
    try:
        for passes in range(1, 4):
            report.update(_counts(report["checks"]))
            with open(path, "wb") as fh:
                fh.write((json.dumps(_jsonable(report), indent=1) + "\n").encode())
            N._VALUES.clear()
            try:
                detail, ok, err = numbers_trace_to_outputs(), True, None
            except Exception as e:  # noqa: BLE001
                detail, ok, err = None, False, f"{type(e).__name__}: {e}"
            detail = _jsonable(detail)
            if detail == entry.get("detail") and ok == entry["pass"] and err == entry.get("error"):
                break
            changed = True
            entry["detail"], entry["pass"] = detail, ok
            entry.pop("error", None)
            if err:
                entry["error"] = err
    finally:
        if keep is None:
            N.PATH_OVERRIDE.pop(N.SELF_REPORT, None)
        else:
            N.PATH_OVERRIDE[N.SELF_REPORT] = keep
        N._VALUES.clear()
    report.update(_counts(report["checks"]))
    if run_log is not None:
        run_log.update(_counts(run_log["results"]))
    return {"rerun": True, "changed": changed, "passes": passes}


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
            except SkipCheck as e:
                results.append({"name": name, "group": group, "pass": None, "skipped": str(e)})
                print(f"SKIP  {group:16s} {name}  ({e})", flush=True)
                continue
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
    counts = _counts
    # a check whose result depends on a private input outside the repository goes to the run log only, so the
    # committed report is the same whether or not that input is on this machine (task 14's review)
    in_report = [r for r in results if r["name"] not in RUN_LOG_ONLY]
    gen = "python -m residualmap.checks" + (" --quick" if quick else "")
    report = {"generated_by": gen,
              "about": "Run-invariant results only; the date, timings and free disk of each run are in "
                       "outputs/chem/checks_run.json (git-ignored).",
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                              "scipy": scipy.__version__, "sklearn": sklearn.__version__, "wntr": wntr.__version__},
              **counts(in_report),
              "run_log_only": {name: why for name, why in RUN_LOG_ONLY.items()},
              "checks": in_report}
    run_log = {"generated_by": gen, "date": time.strftime("%Y-%m-%d %H:%M"), "seconds": round(time.time() - t_all, 1),
               **counts(results), "timing": timing, "results": results}
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
    ap.add_argument("--numbers", action="store_true", help="only run python -m residualmap.numbers: every iteration-4 "
                                                           "number in the README and journal traces to outputs/")
    a = ap.parse_args()
    if a.numbers:
        from .numbers import main as numbers_main
        return numbers_main([])
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
        run_log["numbers_scan_settled"] = settle_numbers_entry(report, run_log)   # one run is enough
        new = (json.dumps(_jsonable(report), indent=1) + "\n").encode()
        with open(REPORT, "wb") as fh:
            fh.write(new)
        run_log["report_changed"] = old != new
    with open(RUN_LOG, "w") as fh:
        json.dump(_jsonable(run_log), fh, indent=1)
    dest = RUN_LOG if a.quick else REPORT
    print(f"\n{report['n_pass']} passed, {report['n_fail']} failed, {report['n_skipped']} skipped in {run_log['seconds']} s "
          f"-> {os.path.relpath(dest, REPO)}" + ("" if a.quick else f" ({'changed' if run_log['report_changed'] else 'unchanged'})"))
    for r in run_log["results"]:
        if r["name"] in RUN_LOG_ONLY:
            state = {True: "passed", False: "FAILED", None: f"skipped ({r.get('skipped')})"}[r["pass"]]
            print(f"{r['name']} (run log only, {os.path.relpath(RUN_LOG, REPO)}): {state}")
    return 0 if report["all_pass"] and run_log["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
