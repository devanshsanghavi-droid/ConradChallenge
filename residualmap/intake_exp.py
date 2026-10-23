"""
intake_exp.py: the iteration-5 intake experiments and their records under outputs/intake/ (task 15b on).

    python -m residualmap.intake_exp              # all three task-15b files below
    python -m residualmap.intake_exp elevation    # only outputs/intake/elevation_prv_ky10.json
    python -m residualmap.intake_exp bundled      # only outputs/intake/readiness_bundled.json
    python -m residualmap.intake_exp record       # only outputs/intake/task15b_physics.json
    python -m residualmap.intake_exp timing       # outputs/intake/run_cost_timing.json (not part of 'all'; add
                                                  # --ky10 for ky10's grid, about 4 minutes)
    python -m residualmap.intake_exp samples      # task 16: outputs/intake/task16_samples.json (not part of 'all')
    python -m residualmap.intake_exp snapping     # task 16b: outputs/intake/tap_snapping.json (not part of 'all';
                                                  # run it from a folder outside the repository, as EPANET writes
                                                  # scratch files into the working directory while it runs)
    python -m residualmap.intake_exp zone         # task 16c: outputs/intake/zone_extract_Net3.json (not part of 'all';
                                                  # run it from outside the repository too: the pilot's validation
                                                  # writes EPANET's temp.inp, .rpt and .bin to the working directory)
    python -m residualmap.intake_exp fidelity [Net3 ky10 ky4]   # task 19: the fidelity ladder, outputs/intake/
                                                  # fidelity_<net>.json (Net3 minutes, ky10 about 17, ky4 about 35;
                                                  # from outside the repository, as EPANET writes scratch files)
    python -m residualmap.intake_exp boundary     # task 19: ky10's boundary valve, picked once (about 2 minutes),
                                                  # outputs/intake/boundary_valve_ky10.json
    python -m residualmap.intake_exp form         # task 19: the operator form's record, outputs/intake/operator_form.json

- elevation_prv_ky10.json: ky10 (five PRVs) with its own junction elevations and with every junction at one elevation
  (the mean of its junctions'), 7-day water-age runs set up as every run is: how many junctions have an hour of the
  scored day whose water age moves by more than 2 h, and the largest move.  The readiness report quotes it when a
  model with pressure-control valves has a junction at elevation 0, or one elevation everywhere.
- readiness_bundled.json: importers.prepare_inp on the six networks WNTR bundles (Net1, Net2, Net3, Net6, ky4, ky10),
  so the journal's numbers on them trace to a file.  The app never prepares the bundled examples; this is a record,
  not a step of any run.  Run times vary, so they stay out (Readiness.bundled_dict) and go to the checks' run log.
- task15b_physics.json: each task-15b check's result on the fixtures (residualmap/fixtures.py) and the bundled
  networks, part by part; the task15b checks in residualmap/checks.py rerun each part and compare.
- run_cost_timing.json: the run cost's time estimate against fresh grid builds, timed, so it varies from run to run
  and no check reruns it.  It holds the pre-registered scored run of run_cost_gate (which missed its bar), the
  review's full run, and the repeats this command times (Net3, and ky10 with --ky10).

- task16_samples.json (task 16): the grab-log reader on the task-16 fixtures (column mapping, censored readings,
  separators and decimal commas, ambiguous dates, encodings, Excel against CSV, analyzer thinning, the UTC warning),
  the clock-conversion helpers, the window the app uses and every template, part by part; the task16 checks rerun each
  part and compare.  The readings are invented, at Net3's junctions.
- tap_snapping.json (task 16b, stopped at its pre-registered bar): tap mapping by coordinates on ky4 (the scored run
  of the bar, which missed: 0.967 of 2,000 taps got the right junction at 5 m against 0.98), the synthetic urban grid
  and its variants, ky4 drawn half the other way, Net3's refusals, the georeference against a Lambert conformal conic
  truth, the review's cases, and a measurement of the ky4 misses made after the stop (not a rescoring); the task16b
  checks rerun each part and compare.

- zone_extract_Net3.json (task 16c): per-inlet residuals on Net3 (two reservoirs): ratios (1, 1) as today's path, a
  fresh grid at them against the cached grid's recorded sha256; the superposition of the two inlets at one member (kb
  0.40, kw 0.70 m/day, the truth's default and off the grid's kw axis; a member on the grid added after the scored run)
  against the pre-registered bars, at the file's quality tolerance and at 1e-9, with the dose axis and the map's
  change for information; leave-3-sites-out validation on a synthetic fixed route (10 taps, no rotating taps), each
  group fit's training size recorded; an online analyzer at a zone inlet read as a boundary condition, with the
  refusals and the review's cases (a declared inlet with no usable reading in the window, roles that mark no inlet);
  the task16c checks rerun each part.

- operator_form.json (task 19, stopped at its pre-registered bar): ky4's clean drawing completed by the operator form
  against the original (part 'complete', the scored run, which missed), the characterization after the stop (part
  'after_the_stop', not a rescoring), ky10's PRVs restored and its boundary valve closed by the form (parts 'prvs' and
  'boundary'), and the refusals (parts 'no_source' and 'bad_rows'); boundary_valve_ky10.json, the one search that
  picked the boundary pipe; fidelity_<net>.json, the fidelity ladder (its protocol is written above its code, before
  the scored runs), Net3 rerun by the saved check, ky4's and ky10's grids checked by sha256.

Every file is built from public models (WNTR's bundled networks, edited by fixtures.py), so the results are converter
self-consistency on public models, never a test on a utility's own files.
"""
from __future__ import annotations

import argparse
import atexit
import contextlib
import functools
import json
import math
import os
import re
import shutil
import sys
import tempfile
import warnings

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "outputs", "intake")
READINESS_BUNDLED = os.path.join(OUT, "readiness_bundled.json")
ELEVATION_PRV = os.path.join(OUT, "elevation_prv_ky10.json")
RECORD_15B = os.path.join(OUT, "task15b_physics.json")
BUNDLED = ("Net1", "Net2", "Net3", "Net6", "ky4", "ky10")
ELEVATION_MOVE_H = 2.0                 # a junction 'moves' when an hour of its scored day changes by more than this
C_RECOVERY_TOL = 0.5                   # inp_headloss_formula: the recovered C of every matched pipe, within this
KB_PAIR = (0.2, 0.8)                   # inp_reaction_orders_fixed: two bulk rates, 1/day (probe d2)
MEMBER = (0.40, 0.70, 1.0)             # inp_bulk_lines_cleared: one grid member, kb 1/day, kw m/day, gamma
MOCK_JUNCTIONS = 10000                 # run_cost_gate: a junction count past every bundled network
MOCK_FREE_GB = (3.0, 6.0)              # ... with the free disk stood in at these values
MOCK_WORKERS = 11                      # ... and this many workers, the plan's feasibility review's machine (12 cores),
                                       # so the record does not depend on the machine that reruns it
TEMP_FILE_NETWORKS = ("ky4", "ky10", "Net6")   # one 7-day run's EPANET files measured (importers.TEMP_MB_PER_JUNCTION)
RUN_COST_TIMING = os.path.join(OUT, "run_cost_timing.json")
TIMING_REPEATS = 5                     # Net3 pairs (one timed run, one fresh grid) per 'timing' run
# The pre-registered scored run of run_cost_gate (task 15b), and the review's full run of the saved checks, as their
# logs give them (kept outside the repository, with the build's other scratch files).  Not rerun: timings vary.
SCORED_RUN = {"what": "run_cost_gate's pre-registered scored run, in the task-15b check run on 18 October 2026, "
                      "before the decision to keep the timing out of the saved check",
              "network": "Net3", "seconds_per_run": 0.032, "runs": 675, "effective_speedup_assumed": 4.0,
              "estimate_s": 5.4, "fresh_grid_s": 11.34, "fresh_from": "fresh_grid_equals_cached_grid, earlier in the "
              "same check run", "ratio": 2.1, "bar": "estimate within a factor of 2 of the fresh grid time "
              "(ratio 0.5 to 2.0)", "passed": False}
REVIEW_FULL_RUN = {"what": "the review's full run of the saved checks on the built tree, 18 October 2026 (its run log)",
                   "network": "Net3", "seconds_per_run": 0.03, "runs": 675, "estimate_s": 5.06,
                   "fresh_grid_s": 9.02, "ratio": 1.783, "passed_the_bar": True,
                   "counts_as": "a characterization after the stop, not a rescoring"}
# Copied in after task 15b's second review, so the journal's numbers about them trace to this file (R3).  Not rerun.
FEASIBILITY_REVIEW = {"what": "the plan's feasibility review, 18 October 2026, on this machine (11 workers): one grid "
                              "of 675 runs and one 7-day run of each network; the effective speedup is 675 times "
                              "the run over the grid time (the private plan, section 'Feasibility review')",
                      "Net3": {"runs": 675, "grid_s": 11.8, "seconds_per_run": 0.041, "effective_speedup": 2.3},
                      "ky10": {"runs": 675, "grid_s": 242.0, "seconds_per_run": 1.77, "effective_speedup": 4.9}}
STANDALONE_BUILDS = {"what": "the implementer's standalone fresh builds of Net3's 675-run grid, 18 October 2026, "
                             "outside any check run; two build logs, kept outside the repository with the build's "
                             "other scratch files",
                     "fresh_grid_s": [10.15, 11.31, 9.35, 9.66, 9.74],
                     "fresh_grid_s_min": 9.35, "fresh_grid_s_max": 11.31, "counts_as": "characterizations, not scored"}


def _quiet():
    warnings.filterwarnings("ignore")


_WORK: list[str] = []


def _work_dir() -> str:
    """One temporary folder per process for the prepared copies of the bundled networks (removed at exit)."""
    if not _WORK:
        d = tempfile.mkdtemp(prefix="rm_t15b_bundled_")
        _WORK.append(d)
        atexit.register(shutil.rmtree, d, True)
    return _WORK[0]


@functools.lru_cache(maxsize=None)
def bundled(name: str):
    """(prepared path, Readiness) of a network WNTR bundles, prepared once per process."""
    _quiet()
    from .importers import prepare_inp
    from .simulate import LIB
    return prepare_inp(os.path.join(LIB, f"{name}.inp"), os.path.join(_work_dir(), "prepared"))


def _prepare(path: str, out: str, name: str | None = None):
    from .importers import prepare_inp
    return prepare_inp(path, out, name=name)


def _mean_member(path: str, kb: float, kw: float = 0.70, gamma: float = 1.0) -> np.ndarray:
    """One grid member (simulate_nominal_chlorine) of a file, hour x junction, run in a temporary folder."""
    from .simulate import simulate_nominal_chlorine
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        return simulate_nominal_chlorine(path, kb, kw, gamma, file_prefix=os.path.join(d, "m")).values


# ----------------------------------------------------------------------------- elevation_prv_ky10.json
def elevation_prv(name: str = "ky10") -> dict:
    """Water age on the scored day with the file's elevations and with every junction at their mean."""
    _quiet()
    import wntr
    from .simulate import _last_day, load
    def age(const: bool):
        wn = load(name)
        wn.options.quality.parameter = "AGE"
        if const:
            e = float(np.mean([j.elevation for _, j in wn.junctions()]))
            for _, j in wn.junctions():
                j.elevation = e
        with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
            r = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "e"))
        return _last_day(r.node["quality"], wn.junction_name_list) / 3600.0, wn
    a0, wn = age(False)
    a1, _ = age(True)
    d = (a1 - a0).abs()
    moved = d.max() > ELEVATION_MOVE_H
    pv = [v for v in wn.valve_name_list if str(wn.get_link(v).valve_type) in ("PRV", "PSV", "PBV")]
    elev_m = float(np.mean([j.elevation for _, j in wn.junctions()]))
    return {"about": "ky10 (public, bundled with WNTR) with its own junction elevations and with every junction at one "
                     "elevation (their mean): water age at every junction and hour of the scored day (the last of a "
                     "7-day run, set up as every run is). Written by python -m residualmap.intake_exp elevation.",
            "network": name, "junctions": len(d.columns), "pressure_valves": len(pv),
            "constant_elevation_m": round(elev_m, 3), "threshold_h": ELEVATION_MOVE_H,
            "junctions_moved": int(moved.sum()), "max_change_h": round(float(d.values.max()), 3),
            "mean_change_h": round(float(d.values.mean()), 4),
            "moved_junctions": sorted(str(j) for j in d.columns[moved.values])}


# ----------------------------------------------------------------------------- readiness_bundled.json
def readiness_bundled() -> dict:
    nets = {}
    for n in BUNDLED:
        path, rep = bundled(n)
        nets[n] = {"same_file": rep.same_file, **rep.bundled_dict()}
    return {"about": "importers.prepare_inp on the six networks WNTR bundles, as a record (the app never prepares the "
                     "bundled examples). Run times and the worker count are left out; they vary by run and machine. "
                     "Written by python -m residualmap.intake_exp bundled.", "networks": nets}


# ----------------------------------------------------------------------------- task15b_physics.json, part by part
def _rows(rep) -> dict:
    return {i: v for i, v, _ in rep.rows}


def part_reaction_orders(paths: dict) -> dict:
    """ky10's orders fixed, and the fixed copy responds to the bulk rate while the original does not; Net6 too, with
    the size warning."""
    from .simulate import LIB
    p10, r10 = bundled("ky10")
    orig = {kb: float(_mean_member(os.path.join(LIB, "ky10.inp"), kb).mean()) for kb in KB_PAIR}
    fixed = {kb: float(_mean_member(p10, kb).mean()) for kb in KB_PAIR}
    p6, r6 = bundled("Net6")
    return {"ky10": {"orders_fixed": r10.flags["reactions"]["orders_fixed"], "fixes": r10.fixes,
                     "mean_original": {str(k): round(v, 4) for k, v in orig.items()},
                     "mean_fixed": {str(k): round(v, 4) for k, v in fixed.items()},
                     "original_same_at_both": abs(orig[KB_PAIR[0]] - orig[KB_PAIR[1]]) < 1e-9,
                     "fixed_difference": round(fixed[KB_PAIR[0]] - fixed[KB_PAIR[1]], 4)},
            "Net6": {"orders_fixed": r6.flags["reactions"]["orders_fixed"], "junctions": r6.run_cost["junctions"],
                     "large": r6.run_cost["large"], "size_warning": any(w.startswith("Large model:") for w in r6.warnings)}}


def part_bulk_lines(paths: dict, out: str) -> dict:
    """net3_bulk_lines.inp: cleared, and its grid member equals Net3's; the original does not (the control)."""
    from .simulate import LIB
    p, rep = _prepare(paths["net3_bulk_lines.inp"], out)
    ref = _mean_member(os.path.join(LIB, "Net3.inp"), *MEMBER)
    got = _mean_member(p, *MEMBER) if p else None
    raw = _mean_member(paths["net3_bulk_lines.inp"], *MEMBER)
    return {"ok": rep.ok, "fixes": rep.fixes, "bulk_lines": rep.flags["reactions"]["bulk_lines"],
            "tank_lines": rep.flags["reactions"]["tank_lines"],
            "limiting_potential_fixed": rep.flags["reactions"]["limiting_potential_fixed"],
            "member": {"kb": MEMBER[0], "kw": MEMBER[1], "gamma": MEMBER[2]},
            "junction_hours": int(ref.size),
            "max_diff_prepared_vs_net3": float(np.abs(got - ref).max()) if got is not None else None,
            "max_diff_original_vs_net3": round(float(np.abs(raw - ref).max()), 4)}


def part_headloss(paths: dict, out: str) -> dict:
    """net3_dw.inp converted, its gate's numbers, and each pipe's recovered C against Net3's; net3_cm.inp blocked."""
    _quiet()
    import wntr
    from . import fixtures as FX
    p, rep = _prepare(paths["net3_dw.inp"], out)
    net3 = FX._wn("Net3")
    _, tiny = FX.net3_dw_roughness()
    got = wntr.network.WaterNetworkModel(p) if p else None
    err, tiny_rows = [], []
    if got is not None:
        for pn in net3.pipe_name_list:
            c0, c1 = float(net3.get_link(pn).roughness), float(got.get_link(pn).roughness)
            if pn in tiny:
                tiny_rows.append({"pipe": pn, "c_net3": c0, "c_recovered": round(c1, 3), "below": round(c0 - c1, 3)})
            else:
                err.append(abs(c1 - c0))
    pc, rc = _prepare(paths["net3_cm.inp"], out)
    pu, ru = _prepare(paths["net3_dw_uniform.inp"], out)          # the review's case: one roughness on every pipe
    pl, rl = _prepare(paths["net3_dw_plastic.inp"], out)          # the second review's: two plastic roughness values
    conv = rep.flags.get("headloss_conversion", {})
    return {"dw": {"ok": rep.ok, "headloss_in_copy": got.options.hydraulic.headloss if got else None,
                   "old_pipe_row_source": next((src for i, _, src in rep.rows if i == "Old-pipe factor"), None),
                   "gate": conv, "pipes": len(net3.pipe_name_list), "pipes_scored": len(err),
                   "max_c_error": round(max(err), 4) if err else None,
                   "within_tolerance": sum(e <= C_RECOVERY_TOL for e in err),
                   "pipes_reported_not_scored": tiny_rows,
                   "max_below_unmatched": round(max(r["below"] for r in tiny_rows), 3) if tiny_rows else None},
            "cm": {"blocked": pc is None and bool(rc.blocking), "message": rc.blocking[:1]},
            "dw_uniform": {"ok": ru.ok, "file_roughness_values": ru.flags.get("c_values"),
                           "converted_c_values": ru.flags.get("converted_c_values"),
                           "uniform_c": ru.flags.get("uniform_c"), "grid": ru.run_cost.get("grid"),
                           "gate": ru.flags.get("headloss_conversion", {}),
                           "row": _rows(ru).get("Old-pipe factor", "")},
            "dw_plastic": {"ok": rl.ok, "roughness_mm": list(FX.DW_PLASTIC_MM),
                           "file_roughness_values": rl.flags.get("c_values"), "uniform_c": rl.flags.get("uniform_c"),
                           "grid": rl.run_cost.get("grid"), "row": _rows(rl).get("Old-pipe factor", "")}}


CLOCK_EDGES = {"23:59:50": 23 + 59 / 60 + 50 / 3600, "0:00:20": 20 / 3600, "5:59:59": 5 + 59 / 60 + 59 / 3600,
               "6:30": 6.5}   # the second review's starts within a minute of midnight, and two controls


def _clock_row(rep) -> str:
    """The start clock row up to its HH:MM: '6.5 h (06:30)'."""
    v = _rows(rep)["Start clock time"]
    return v[:v.index(")") + 1]


def _clock_edge(hours_: float) -> dict:
    """The start clock row and model-hours warning importers._clock_report gives for a start, in hours."""
    from .importers import Readiness, _clock_report
    r = Readiness(name="clock")
    _clock_report(r, hours_)
    w = [x for x in r.warnings if x.startswith("The model starts at")]
    return {"row": _clock_row(r), "model_hours_warning": bool(w), "warning_says": w[0][:w[0].index(",")] if w else None,
            "warning_offset": w[0][w[0].index("Add "):w[0].index(" to a model hour")] if w else None}


def part_clock(paths: dict, out: str) -> dict:
    """Net2's 8 h start and its in-file source; net3_clock6.inp's 6 h; net3_weekly_pattern.inp's weekday warning;
    since the second review, starts within a minute of midnight (net3_clock_235950.inp, net3_clock_000020.inp, and
    the rounding of the row and the warning)."""
    _, r2 = bundled("Net2")
    _, r6 = _prepare(paths["net3_clock6.inp"], out)
    _, rw = _prepare(paths["net3_weekly_pattern.inp"], out)
    hours = lambda r: any(w.startswith("The model starts at") and "model hour" in w for w in r.warnings)
    _, r3 = bundled("Net3")
    near = {}
    for f in ("net3_clock_235950.inp", "net3_clock_000020.inp"):
        _, rn = _prepare(paths[f], out)
        near[f] = {"ok": rn.ok, "start_clock_h": round(rn.flags["start_clock_h"], 6),
                   "row": _clock_row(rn), "model_hours_warning": hours(rn)}
    return {"Net2": {"start_clock_h": r2.flags["start_clock_h"], "sources_overridden": r2.flags["sources_overridden"],
                     "source_flagged": any("own source" in w for w in r2.warnings),
                     "multi_day_patterns": r2.flags["multi_day_patterns"], "model_hours_warning": hours(r2)},
            "net3_clock6": {"start_clock_h": r6.flags["start_clock_h"], "ok": r6.ok, "model_hours_warning": hours(r6)},
            "Net3": {"start_clock_h": r3.flags["start_clock_h"], "model_hours_warning": hours(r3)},
            "net3_weekly_pattern": {"weekly_patterns": rw.flags["weekly_patterns"],
                                    "weekday_warning": any("weekly pattern" in w and "weekday" in w for w in rw.warnings),
                                    "ok": rw.ok},
            "near_midnight": near, "rounding": {k: _clock_edge(v) for k, v in CLOCK_EDGES.items()}}


def part_prv(paths: dict, out: str) -> dict:
    """ky10 and Net6 flagged at their zero-elevation junctions, ky10_const_elev.inp blocked, Net3 not flagged."""
    out_ = {}
    for n in ("ky10", "Net6", "Net3"):
        _, r = bundled(n)
        out_[n] = {"pressure_valves": r.flags["pressure_valves"], "zero_elevation_junctions":
                   r.flags["zero_elevation_junctions"], "flagged": r.flags["zero_elevation_flagged"],
                   "elevations_needed_at": r.flags["elevations_needed_at"], "ok": r.ok}
    p, rc = _prepare(paths["ky10_const_elev.inp"], out)
    out_["ky10_const_elev"] = {"blocked": p is None and any(b.startswith("Elevations needed at these valves")
                                                             for b in rc.blocking),
                               "constant_elevation": rc.flags.get("constant_elevation"),
                               "quotes_evidence": any("elevation_prv_ky10.json" in b for b in rc.blocking),
                               "fixes": len(rc.fixes),
                               "fixes_heading": next((ln for ln in rc.to_markdown().splitlines()
                                                      if ln.startswith("**Fixes")), None)}
    return out_


def part_coordinates(paths: dict, out: str) -> dict:
    res = {}
    for n in ("Net3", "Net2", "Net6", "ky4", "ky10"):
        res[n] = bundled(n)[1].flags["coordinates"]
    res["net3_no_coords"] = _prepare(paths["net3_no_coords.inp"], out)[1].flags["coordinates"]
    return res


def part_structure(paths: dict, out: str) -> dict:
    res = {}
    for f in ("net3_no_controls.inp", "net3_speed_pattern_pumps.inp", "net3_zero_demand.inp", "net3_two_parts.inp",
              "net3_no_source.inp", "watergems_export_like.inp", "net3_rule_premise.inp"):
        p, r = _prepare(paths[f], out)
        fl = r.flags
        res[f] = {"ok": r.ok, "blocking": r.blocking,
                  "pumps_without_controls": fl.get("pumps_without_controls"),
                  "pump_warning": any(w.startswith("Pumps but no controls") for w in r.warnings),
                  "zero_total_demand": fl.get("zero_total_demand"), "parts": fl.get("parts"),
                  "connector_pipes": fl.get("connector_pipes"), "connector_tags": fl.get("connector_tags"),
                  "connector_row": "Export connectors" in _rows(r)}
    _, r6 = bundled("Net6")          # one pipe of the connector size and no tag: counted, not reported
    res["Net6"] = {"connector_pipes": r6.flags["connector_pipes"], "connector_tags": r6.flags["connector_tags"],
                   "connector_row": "Export connectors" in _rows(r6)}
    return res


def part_uniform_c() -> dict:
    """Net1 (C 100 on every pipe) is flagged and gets the fixed-gamma grid, whose 225 members equal the gamma-0
    members of its full grid bit for bit (both grids built in memory, never cached)."""
    _quiet()
    import itertools
    from .simgp import GRIDS, build_grid_24h
    from .simulate import LIB, nominal_scenario
    _, r1 = bundled("Net1")
    sc = nominal_scenario(os.path.join(LIB, "Net1.inp"))
    pf, zf = build_grid_24h(sc, "full")
    pg, zg = build_grid_24h(sc, "full_fixed_gamma")
    with tempfile.TemporaryDirectory(prefix="rm_t15b_") as d:     # chloramine mode: no fixed-gamma grid there
        from . import fixtures as FX
        from .importers import prepare_inp
        _, rc = prepare_inp(os.path.join(LIB, "Net1.inp"), d, "chloramine")
        # the second review's case: Net1 written as an export, with a 1 ft connector of C 130 at its pump
        _, rx = prepare_inp(FX.net1_export_like(), d, name="net1_export_like.inp")
    idx = {tuple(p): k for k, p in enumerate(pf)}
    same = [bool(np.array_equal(zg[k], zf[idx[tuple(p)]])) for k, p in enumerate(pg)]
    return {"Net1": {"uniform_c": r1.flags["uniform_c"], "c_values": r1.flags["c_values"], "grid": r1.run_cost["grid"],
                     "runs_per_dose": r1.run_cost["runs_per_dose"],
                     "row": _rows(r1).get("Old-pipe factor", "")[:60],
                     "row_says_held_at_0": "it is held at 0" in _rows(r1).get("Old-pipe factor", "")},
            "Net1_chloramine": {"uniform_c": rc.flags["uniform_c"], "grid": rc.run_cost["grid"],
                                "runs_per_dose": rc.run_cost["runs_per_dose"],
                                "row_says_held_at_0": "it is held at 0" in _rows(rc).get("Old-pipe factor", ""),
                                "row_says_chloramine_grid_varies_it": "chloramine grid still varies it"
                                                                      in _rows(rc).get("Old-pipe factor", "")},
            "Net1_export_like": {"ok": rx.ok, "connector_pipes": rx.flags["connector_pipes"],
                                 "connector_tags": rx.flags["connector_tags"],
                                 "connector_row": "Export connectors" in _rows(rx), "uniform_c": rx.flags["uniform_c"],
                                 "c_values": rx.flags["c_values"], "grid": rx.run_cost["grid"],
                                 "row": _rows(rx).get("Old-pipe factor", "")[:100]},
            "fixed_gamma_members": len(pg), "gamma_values": sorted({float(p[2]) for p in pg}),
            "members_bit_identical_to_full_gamma0": int(sum(same)),
            "full_gamma0_members": sum(1 for p in pf if p[2] == 0.0),
            "axes": [len(a) for a in GRIDS["full_fixed_gamma"]],
            "same_order_as_full": [tuple(p) for p in pg] == [tuple(p) for p in itertools.product(*GRIDS["full_fixed_gamma"])],
            "committed_networks_c_values": {n: bundled(n)[1].flags["c_values"] for n in ("Net3", "Net2", "ky4")}}


def temp_files_one_run(name: str) -> dict:
    """The EPANET files (.inp, .rpt, .bin) one 7-day grid member of a bundled network writes, in MB (what
    importers.TEMP_MB_PER_JUNCTION scales by the junction count): their size depends on the network only."""
    from .simulate import LIB, load, simulate_nominal_chlorine
    _quiet()
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        simulate_nominal_chlorine(os.path.join(LIB, f"{name}.inp"), 0.40, 0.70, 1.0, file_prefix=os.path.join(d, "m"))
        size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
    return {"junctions": load(name).num_junctions, "mb": round(size / 1e6, 1)}


def part_run_cost() -> dict:
    """Net6's hard warning, the size gate on a 10,000-junction chloramine grid with the free disk stood in at 3 GB
    (refused) and 6 GB (accepted), with MOCK_WORKERS workers, and one run's EPANET files on ky4, ky10 and Net6.  The
    time estimate against a fresh grid is timed, so it is in the checks' run log and run_cost_timing.json only."""
    from . import importers as I
    _, r6 = bundled("Net6")
    cost = I.run_cost(MOCK_JUNCTIONS, "chloramine", workers=MOCK_WORKERS)
    keep = I.disk_free_gb
    gate = {}
    try:
        for gb in MOCK_FREE_GB:
            I.disk_free_gb = lambda path, gb=gb: gb
            ok, msg = I.size_gate(cost)
            gate[f"{gb:g}"] = {"accepted": ok, "message_says_floor": "below the 1.5 GB" in msg}
    finally:
        I.disk_free_gb = keep
    return {"Net6": {"junctions": r6.run_cost["junctions"], "large": r6.run_cost["large"],
                     "hard_warning": any(w.startswith(I.LARGE_WARNING_START) for w in r6.warnings),
                     "grid_file_mb": r6.run_cost["grid_file_mb"]},
            "mock": {"junctions": MOCK_JUNCTIONS, "disinfectant": "chloramine", "runs": cost["runs_per_dose"],
                     "workers": cost["workers"], "grid_file_mb": cost["grid_file_mb"],
                     "temp_files_mb": cost["temp_files_mb"], "gate": gate},
            "kb_per_junction": {"free_chlorine": I.run_cost(1000, "free_chlorine")["grid_file_mb"],
                                "chloramine": I.run_cost(1000, "chloramine")["grid_file_mb"]},
            "temp_files_one_run": {n: temp_files_one_run(n) for n in TEMP_FILE_NETWORKS},
            "temp_mb_per_junction_used": round(I.TEMP_MB_PER_JUNCTION * 1000, 3)}


def _decode_task15(data: bytes) -> str:
    """Task 15's decoding of a file without a byte-order mark that fails strict UTF-8: cp1252 throughout (a baseline
    for the encoding fuzz only; the fuzz files have no byte-order mark and are not UTF-16)."""
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("cp1252", errors="residualmap_cp1252_gaps")


def _decode_first_fix(data: bytes) -> str:
    """The rule task 15b's first review fix used, kept here as a baseline for the encoding fuzz only: UTF-8 with the
    stray bytes as cp1252 when the valid multi-byte characters, counted by occurrence, are at least as many as the
    stray bytes, else cp1252 throughout (replaced after the second review: importers._mixed_utf8)."""
    from .importers import _NON_ASCII, _STRAY, _cp1252_byte
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        pass
    t = data.decode("utf-8", errors="surrogateescape")
    stray = len(_STRAY.findall(t))
    if stray and len(_NON_ASCII.findall(t)) - stray >= max(1, stray):
        return _STRAY.sub(_cp1252_byte, t)
    return data.decode("cp1252", errors="residualmap_cp1252_gaps")


def part_leftovers(paths: dict, out: str) -> dict:
    """Task 15's four leftovers: a SWMM model whose SWMM-only sections come after 64 KB, a zip holding a model renamed
    .txt, a UTF-8 file with one stray cp1252 byte, and a content-named copy an earlier version accepted of a model
    that is blocked now.  Since the second review, the encoding controls (fixtures.ENCODING_CONTROLS) and the
    review's encoding fuzz (fixtures.encoding_fuzz, seed 11: 150 files in cp1252 throughout, 150 in UTF-8 with stray
    bytes; decode_text only)."""
    _quiet()
    import hashlib
    import wntr
    from . import fixtures as FX
    from .importers import HEAD_BYTES, decode_text, prepare_inp, sniff
    from .importers import _sniff_content
    sw = FX.FIXTURES["swmm_after_64kb.inp"]()
    s_sw, p_sw = sniff(sw, "storm.inp"), prepare_inp(sw, out, name="storm.inp")[1]
    head_only = _sniff_content(sw[:HEAD_BYTES], "storm.inp").kind     # the first 64 KB alone: what task 15 saw
    zt = FX.FIXTURES["zip_model_as_txt.zip"]()
    s_zt = sniff(zt, "zone.zip")
    pz, rz = prepare_inp(zt, out, name="zone.zip")
    mixed = FX.FIXTURES["utf8_stray_cp1252.inp"]()
    pm, rm = prepare_inp(mixed, out, name="mixed.inp")
    want = [new for _, new in FX.ACCENTED_RENAMES]
    wn_m = wntr.network.WaterNetworkModel(pm) if pm else None
    with tempfile.TemporaryDirectory(prefix="rm_t15b_") as d:       # the control: cp1252 throughout, as task 15 did
        f = os.path.join(d, "as_cp1252.inp")
        with open(f, "w", encoding="utf-8") as fh:
            fh.write(mixed.decode("cp1252", errors="replace"))
        try:
            wntr.network.WaterNetworkModel(f)
            control = "loads"
        except Exception as e:  # noqa: BLE001
            control = f"raises {type(e).__name__}"
    pair = FX.FIXTURES["net3_cp1252_pair.inp"]()
    pp, rp = prepare_inp(pair, out, name="pair.inp")
    wn_p = wntr.network.WaterNetworkModel(pp) if pp else None
    from .importers import _NON_ASCII, _STRAY
    t_pair = pair.decode("utf-8", errors="surrogateescape")
    n_stray = len(_STRAY.findall(t_pair))
    controls = {}                  # the second review's encoding cases, each on both sides of the old count rule
    for f, how, ids in FX.ENCODING_CONTROLS:
        data = FX.FIXTURES[f]()
        pc_, rc_ = prepare_inp(data, out, name=f)
        wn_c = wntr.network.WaterNetworkModel(pc_) if pc_ else None
        t_c = data.decode("utf-8", errors="surrogateescape")
        controls[f] = {"must_decode_as": how, "decoded_as": _rows(rc_).get("Text encoding"), "ok": rc_.ok,
                       "ids_kept": bool(wn_c) and all(i in wn_c.junction_name_list for i in ids),
                       "mixed_warning": any(w.startswith("The file mixes two text encodings") for w in rc_.warnings),
                       "valid_utf8_chars": len(_NON_ASCII.findall(t_c)) - len(_STRAY.findall(t_c)),
                       "bytes_not_utf8": len(_STRAY.findall(t_c))}
    controls["utf8_pasted_chance_pair.inp"]["pasted_comment_kept"] = (
        FX.PASTED_PAIR in decode_text(FX.FIXTURES["utf8_pasted_chance_pair.inp"]())[0])
    zero = lambda: {"cp1252_throughout": 0, "utf8_with_stray": 0}
    fuzz = {"seed": 11, "files_per_kind": 150, "ids_changed": zero(), "ids_changed_task15_rule": zero(),
            "ids_changed_first_fix_rule": zero()}
    for kind, data, ids in FX.encoding_fuzz(11, 150):
        for key, dec in (("ids_changed", lambda b: decode_text(b)[0]), ("ids_changed_task15_rule", _decode_task15),
                         ("ids_changed_first_fix_rule", _decode_first_fix)):
            text = dec(data)
            fuzz[key][kind] += not all(i in text for i in ids)
    cm = FX.FIXTURES["net3_cm.inp"]()
    with tempfile.TemporaryDirectory(prefix="rm_t15b_") as d:       # a copy task 15 accepted (it only warned on C-M)
        stale = f"{hashlib.sha1(cm).hexdigest()[:8]}_net3_cm.inp"
        for n, b in ((stale, cm), ("other_model.inp", b"unrelated")):
            with open(os.path.join(d, n), "wb") as fh:
                fh.write(b)
        before = sorted(os.listdir(d))
        p_cm, r_cm = prepare_inp(cm, d, name="net3_cm.inp")
        after = sorted(os.listdir(d))
    return {"swmm_after_64kb": {"first_64kb_alone": head_only, "kind": s_sw.kind, "not_a_model": p_sw.not_a_model,
                                "swmm_advice": "SWMM" in p_sw.advice(), "first_swmm_section_at": sw.index(b"[OUTFALLS]")},
            "zip_model_as_txt": {"kind": s_zt.kind, "members": s_zt.members, "ok": rz.ok, "prepared": rz.prepared,
                                 "same_network_as_net3": bool(pz) and FX._same_network(pz, FX._wn("Net3"))},
            "utf8_stray_cp1252": {"decoded_as": _rows(rm).get("Text encoding"), "ok": rm.ok,
                                  "mixed_warning": any(w.startswith("The file mixes two text encodings")
                                                       for w in rm.warnings),
                                  "accented_ids_kept": bool(wn_m) and all(i in wn_m.junction_name_list for i in want),
                                  "decoded_as_cp1252_throughout": control},
            "cp1252_with_utf8_pairs": {"decoded_as": _rows(rp).get("Text encoding"), "ok": rp.ok,
                                       "bytes_not_utf8": n_stray,
                                       "valid_utf8_pairs": len(_NON_ASCII.findall(t_pair)) - n_stray,
                                       "same_as_cp1252_throughout": bool(pp) and open(pp, "rb").read().decode("utf-8")
                                       == pair.decode("cp1252"),
                                       "ids_kept": bool(wn_p) and all(i in wn_p.junction_name_list
                                                                      for _, i in FX.CP1252_PAIR_RENAMES),
                                       "mixed_warning": any(w.startswith("The file mixes two text encodings")
                                                            for w in rp.warnings)},
            "encoding_controls": controls, "encoding_fuzz": fuzz,
            "stale_copy": {"before": before, "blocked": p_cm is None, "after": after,
                           "stale_removed": stale not in after, "other_file_kept": "other_model.inp" in after}}


def part_route() -> dict:
    """plan_route on Net1 (9 junctions) with fewer junctions left than the route asks for, on its fixed-gamma prior
    (a grid built into a temporary cache)."""
    _quiet()
    from .features import build_features
    from .route import plan_route
    from .simgp import SimGP24
    from .simulate import LIB, nominal_scenario
    with tempfile.TemporaryDirectory(prefix="rm_t15b_cache_") as d:
        sc = nominal_scenario(os.path.join(LIB, "Net1.inp"))
        m = SimGP24(sc, build_features(sc), seed=0, cache_dir=d, grid="full_fixed_gamma").fit_prior()
        cases = {}
        for label, K, excl in (("K6_8_sampled", 6, 8), ("K6_9_sampled", 6, 9), ("K3_none_sampled", 3, 0)):
            r = plan_route(m, K, exclude=list(sc.junctions[:excl]), with_age=True)
            cases[label] = {"rows": len(r), "columns": list(r.columns),
                            "junctions": [str(j) for j in r.junction]}
    return {"junctions": len(sc.junctions), "cases": cases}


def record() -> dict:
    _quiet()
    from . import fixtures as FX
    names = [n for n in FX.FIXTURES if n in FX.PHYSICS_FIXTURES]
    with tempfile.TemporaryDirectory(prefix="rm_t15b_rec_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), names)
        out = os.path.join(d, "prepared")
        rec = {"about": "Task 15b's checks of an imported model's physics, on files built by residualmap/fixtures.py "
                        "from WNTR's bundled public networks and on those networks themselves (converter "
                        "self-consistency on public models, not real utility files). Written by python -m "
                        "residualmap.intake_exp record; the task15b checks rerun each part and compare.",
               "reaction_orders": part_reaction_orders(paths), "bulk_lines": part_bulk_lines(paths, out),
               "headloss": part_headloss(paths, out), "clock": part_clock(paths, out), "prv": part_prv(paths, out),
               "coordinates": part_coordinates(paths, out), "structure": part_structure(paths, out),
               "uniform_c": part_uniform_c(), "run_cost": part_run_cost(), "leftovers": part_leftovers(paths, out),
               "route": part_route()}
    return json.loads(json.dumps(rec))


# ----------------------------------------------------------------------------- run_cost_timing.json (timed: not rerun)
def _timing_pair(name: str) -> dict:
    """One timed 7-day grid member (importers._timed_member_run, as prepare_inp times it) and one fresh 'full' grid
    build in memory (simgp.build_grid_24h, as fresh_grid_equals_cached_grid builds it; nothing is cached), with the
    estimate importers.run_cost makes from that run."""
    import time
    from . import importers as I
    from .simgp import build_grid_24h
    from .simulate import LIB, nominal_scenario
    secs = I._timed_member_run(os.path.join(LIB, f"{name}.inp"), 1.0)
    cost = I.run_cost(1, "free_chlorine", seconds_per_run=secs)
    sc = nominal_scenario(name)
    t0 = time.perf_counter()
    build_grid_24h(sc, "full")
    fresh = time.perf_counter() - t0
    runs = cost["runs_per_dose"]
    return {"seconds_per_run": round(secs, 4), "estimate_s": round(cost["estimate_s"], 2), "fresh_grid_s": round(fresh, 2),
            "ratio": round(fresh / cost["estimate_s"], 3), "effective_speedup": round(runs * secs / fresh, 2),
            "within_factor_2": 0.5 <= fresh / cost["estimate_s"] <= 2.0}


def run_cost_timing(repeats: int = TIMING_REPEATS, ky10: bool = False) -> dict:
    """The time estimate against fresh grid builds on this machine, with the scored run and the review's full run."""
    import time
    from . import importers as I
    _quiet()
    rep_ = [_timing_pair("Net3") for _ in range(repeats)]
    ratios = [r["ratio"] for r in rep_]
    out = {"about": "The run cost's time estimate (importers.estimate_seconds: one timed 7-day run times the runs, over "
                    "an assumed effective speedup of 4) against fresh grid builds of 675 runs, in seconds. Timed, so "
                    "it varies from run to run and no check reruns it; written by python -m residualmap.intake_exp "
                    "timing. The pre-registered bar (run_cost_gate, task 15b): Net3's estimate within a factor of 2 "
                    "of the fresh grid time. Its scored run missed; the repeats below are characterizations after "
                    "the stop, not rescorings.",
           "written": time.strftime("%Y-%m-%d %H:%M"), "cpu_count": os.cpu_count(),
           "workers": max(1, (os.cpu_count() or 2) - 1), "effective_speedup_assumed": I.EFFECTIVE_SPEEDUP,
           "scored_run": SCORED_RUN, "review_full_run": REVIEW_FULL_RUN,
           "feasibility_review": FEASIBILITY_REVIEW, "standalone_builds": STANDALONE_BUILDS,
           "net3_repeats": rep_,
           "net3_repeats_summary": {"n": len(rep_), "ratio_min": min(ratios), "ratio_max": max(ratios),
                                    "within_factor_2": sum(r["within_factor_2"] for r in rep_),
                                    "fresh_grid_s_min": min(r["fresh_grid_s"] for r in rep_),
                                    "fresh_grid_s_max": max(r["fresh_grid_s"] for r in rep_),
                                    "effective_speedup_min": min(r["effective_speedup"] for r in rep_),
                                    "effective_speedup_max": max(r["effective_speedup"] for r in rep_)}}
    if ky10:
        out["ky10"] = _timing_pair("ky10")
    return out


def fixed_timing_records() -> dict:
    """The existing run_cost_timing.json with the records copied from logs set to this module's values, in the order
    run_cost_timing writes them; the timed records are kept as they are (nothing is timed)."""
    with open(RUN_COST_TIMING) as fh:
        old = json.load(fh)
    fixed = {"scored_run": SCORED_RUN, "review_full_run": REVIEW_FULL_RUN, "feasibility_review": FEASIBILITY_REVIEW,
             "standalone_builds": STANDALONE_BUILDS}
    out = {k: old[k] for k in ("about", "written", "cpu_count", "workers", "effective_speedup_assumed") if k in old}
    out.update(fixed)
    out.update({k: v for k, v in old.items() if k not in out})
    return out


def _write(path: str, d: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1)
        fh.write("\n")
    print(f"wrote {os.path.relpath(path, REPO)}")


# ----------------------------------------------------------------------------- task 16: grab logs, tap maps, clock time
RECORD_16 = os.path.join(OUT, "task16_samples.json")
TEMPLATES_DIR = os.path.join(REPO, "docs", "templates")
TEMPLATES = {"grab_log_template.csv": ("log", "free_chlorine"), "grab_log_total_template.csv": ("log", "chloramine"),
             "tap_map_by_id_template.csv": ("taps", None), "plant_log_template.csv": ("plant", None)}


def _read_bytes(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def _log(paths: dict, name: str, disinfectant: str = "free_chlorine", taps: str = "task16_tap_map.csv", **kw):
    """importers.read_grab_log on a built fixture, with the task-16 tap map."""
    from .importers import read_grab_log
    return read_grab_log(_read_bytes(paths[name]), name, _read_bytes(paths[taps]), taps, disinfectant, **kw)


def _refusal(paths: dict, name: str, disinfectant: str = "free_chlorine", **kw) -> str | None:
    """The plain message a fixture is refused with, or None when it is read."""
    from .importers import LogRefused
    try:
        _log(paths, name, disinfectant, **kw)
    except LogRefused as e:
        return str(e)
    return None


def _taps() -> dict:
    from . import fixtures as FX
    return dict(FX.T16_TAPS)


def _rows_of(frame) -> list:
    return [[r.junction, int(r.hour), round(float(r.y), 6), r.month, r.date, r.tap_id] for r in frame.itertuples()]


def _expected(rows) -> list:
    """The frame rows a fixture must give, from its own table: (date, hour, minute, tap, reading) in, junction from
    the tap map, model hour = the hour (every model here starts at midnight)."""
    taps = _taps()
    return [[taps[t], int(h), round(float(v), 6), d[:7], d, t] for d, h, _, t, v in rows]


def part_mapping(paths: dict) -> dict:
    """grab_log_column_mapping: the LIMS export with one US datetime column (guessed columns, an explicit date
    format), the censored log, the semicolon and tab logs, and the ambiguous dates refused with the question."""
    from . import fixtures as FX
    out = {}
    g = _log(paths, "lims_one_datetime_us.csv", "chloramine", date_format="MDY")
    want = _expected([(f"2026-{mo:02d}-{d:02d}", h, mi, t, v) for mo, d, h, mi, t, v in FX.LIMS_US_ROWS])
    q = _refusal(paths, "lims_one_datetime_us.csv", "chloramine")
    out["lims_one_datetime_us"] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want,
                                   "mapping": {k: g.mapping[k] for k in ("tap", "reading", "date", "time")},
                                   "species_warning": any("does not say free or total" in w for w in g.warnings),
                                   "without_format_refused": q is not None,
                                   "without_format_asks": bool(q) and "Which is it?" in q,
                                   "with_strptime_same": _rows_of(_log(paths, "lims_one_datetime_us.csv", "chloramine",
                                                                       date_format="%m/%d/%Y %I:%M %p").frame) == want}
    g = _log(paths, "lims_censored.csv")
    used = [r for r in FX.LIMS_CENSORED_ROWS if r[2] in FX.LIMS_CENSORED_EXPECT["used"]]
    want = _expected([(d, int(t[:2]), int(t[3:]), tap, float(v.split()[0].replace(",", "."))) for d, t, tap, v in used])
    out["lims_censored"] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want,
                            "n_censored": len(g.censored), "n_dropped": len(g.dropped),
                            "censored": sorted(g.censored.tap_id), "censored_raw": sorted(g.censored.raw),
                            "censored_reasons_say_so": all(r.startswith("censored") and "not used in the fit" in r
                                                           for r in g.censored.reason),
                            "dropped": [[d["tap_id"], d["reason"]] for d in g.dropped],
                            "as_expected": (sorted(g.censored.tap_id) == sorted(FX.LIMS_CENSORED_EXPECT["censored"])
                                            and [d["tap_id"] for d in g.dropped] == list(FX.LIMS_CENSORED_EXPECT["blank"]))}
    want = _expected([(d, int(t[:2]), int(t[3:]), tap, float(v.replace(",", "."))) for d, t, tap, v in FX.LOG_ROWS])
    for f in ("log_semicolon.csv", "log_tab.txt"):
        g = _log(paths, f)
        out[f.split(".")[0]] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want, "how": g.notes[0]}
    g = _log(paths, "log_semicolon.csv", mapping={"tap": "tap_id", "reading": "free_chlorine_mgL", "date": "date",
                                                   "time": "time"})
    out["log_semicolon_explicit_mapping_same"] = _rows_of(g.frame) == want
    q = _refusal(paths, "ambiguous_dates.csv")
    out["ambiguous_dates"] = {"refused": q is not None, "asks": bool(q) and "Which is it?" in q,
                              "question": q,
                              "with_DMY_rows": len(_log(paths, "ambiguous_dates.csv", date_format="DMY").frame)}
    return out


def part_encodings(paths: dict) -> dict:
    """grab_log_encodings: the cp1252 log (an accented tap ID matched to the UTF-8 tap map) and the log with a UTF-8
    byte-order mark."""
    from . import fixtures as FX
    out = {}
    for f, rows in (("grab_log_cp1252.csv", FX.CP1252_ROWS),
                    ("grab_log_bom.csv", [r for r in FX.CP1252_ROWS if r[2] != "Église_3"])):
        g = _log(paths, f)
        want = _expected([(d, int(t[:2]), int(t[3:]), tap, float(v)) for d, t, tap, v, _ in rows])
        out[f.split(".")[0]] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want, "how": g.notes[0],
                                "header_first": g.mapping["date"]}
    return out


def part_excel(paths: dict) -> dict:
    """excel_log_equals_csv_log: grab_log_typed.xlsx (typed date and time cells, numeric tap IDs that lost their
    leading zeros, written by openpyxl: same-library self-consistency) against its CSV twin."""
    from . import fixtures as FX
    gx, gc = _log(paths, "grab_log_typed.xlsx"), _log(paths, "grab_log_typed.csv")
    want = _expected([(d, h, mi, t, v) for d, h, mi, t, v in FX.TYPED_ROWS])
    lz = [w for w in gx.warnings if "leading zeros" in w]
    return {"xlsx_rows": len(gx.frame), "csv_rows": len(gc.frame), "same_rows": _rows_of(gx.frame) == _rows_of(gc.frame),
            "rows_as_expected": _rows_of(gc.frame) == want, "same_S": gx.frame[["junction", "hour", "y"]].equals(
                gc.frame[["junction", "hour", "y"]]), "leading_zero_warning": bool(lz),
            "leading_zero_taps": sorted({t for _, _, _, t, _ in FX.TYPED_ROWS if t.isdigit() and t.startswith("0")}),
            "warning": lz[0] if lz else None,
            "csv_warnings": len(gc.warnings), "how": gx.notes[0], "label": "same-library self-consistency (openpyxl)"}


def part_analyzer(paths: dict) -> dict:
    """analyzer_thinning: two analyzers at 15-minute steps for 30 days, and three grab samples."""
    from . import fixtures as FX
    from .importers import ANALYZER_HOUR
    g = _log(paths, "analyzer_15min.csv")
    an = g.frame[g.frame.tap_id.isin(FX.ANALYZER_SITES)]
    per = an.groupby(["tap_id", "date"]).size()
    clock = sorted(set(g.clock_min[an.index]))
    return {"readings_in": sum(v["readings"] for v in g.analyzer.values()), "sites": sorted(g.analyzer),
            "kept": len(an), "per_site_per_day_max": int(per.max()), "per_site_per_day_min": int(per.min()),
            "days": int(an.groupby("tap_id").date.nunique().min()), "clock_minutes_kept": clock,
            "at_hour": ANALYZER_HOUR, "grab_rows_kept": int((~g.frame.tap_id.isin(FX.ANALYZER_SITES)).sum()),
            "noted": any("online analyzer" in n for n in g.notes),
            "thinned_listed": sum(d["reason"].startswith("analyzer series at tap") for d in g.dropped),
            "hour_8": sorted(set(_log(paths, "analyzer_15min.csv", analyzer_hour=8.0).clock_min[an.index]))}


def part_utc(paths: dict) -> dict:
    """utc_warning: the log whose daytime samples were written 8 h late warns; the committed example log and the
    semicolon log do not."""
    g = _log(paths, "utc_shifted_log.csv")
    w = [x for x in g.warnings if "UTC" in x]
    plain = _log(paths, "log_semicolon.csv")
    from .importers import read_grab_log
    ex = read_grab_log(_read_bytes(os.path.join(REPO, "docs", "example_grab_log.csv")), "example_grab_log.csv",
                       _read_bytes(os.path.join(REPO, "docs", "example_tap_map.csv")), "example_tap_map.csv")
    return {"warned": bool(w), "warning": w[0] if w else None, "semicolon_log_warned": any("UTC" in x for x in plain.warnings),
            "example_log_warned": any("UTC" in x for x in ex.warnings)}


def part_clock_helpers() -> dict:
    """clock_conversion, its helper half: model_hour, clock_label and daytime_model_hours for a file that starts at
    06:00, at midnight, at 6:30, and within a minute of midnight; the bundled examples' starts."""
    from .importers import clock_label, daytime_model_hours, model_hour, shows_clock
    from .simgp import DAY_HOURS
    bundled = json.load(open(READINESS_BUNDLED))["networks"]
    return {"09:00_at_6h": model_hour(9 * 60, 6.0), "09:59_at_6h": model_hour(9 * 60 + 59, 6.0),
            "03:00_at_6h": model_hour(3 * 60, 6.0), "model_3_label_6h": clock_label(3, 6.0),
            "daytime_6h": daytime_model_hours(6.0), "daytime_6h_labels": [clock_label(h, 6.0) for h in daytime_model_hours(6.0)],
            "daytime_0h_is_DAY_HOURS": daytime_model_hours(0.0) == DAY_HOURS,
            "identity_at_0h": all(model_hour(m, 0.0) == m // 60 for m in range(1440)),
            "daytime_6h30_labels": [clock_label(h, 6.5) for h in daytime_model_hours(6.5)],
            "09:15_at_6h30": model_hour(9 * 60 + 15, 6.5),
            "shows_clock": {k: shows_clock(v) for k, v in (("0", 0.0), ("6", 6.0), ("23:59:50", 23 + 59 / 60 + 50 / 3600),
                                                             ("0:00:20", 20 / 3600))},
            "bundled_start_h": {n: bundled[n]["flags"]["start_clock_h"] for n in ("Net3", "Net2", "ky4")}}


def part_app_inputs(paths: dict) -> dict:
    """What app_grab_log_upload expects, read here without the app: the committed example log's rows in the last three
    months, the unknown tap, the mixed species refused, the 9-month log's window."""
    from .pilot import load_log
    from .importers import read_grab_log
    ex = load_log(os.path.join(REPO, "docs", "example_grab_log.csv"), os.path.join(REPO, "docs", "example_tap_map.csv"),
                  quiet=True)
    months = sorted(set(ex.month))[-3:]
    u = _log(paths, "unknown_tap_log.csv", taps="task16_tap_map.csv")
    mm = _log(paths, "multi_month_log.csv")
    f, _, info = mm.window()
    return {"example_rows": len(ex), "example_rows_in_window": int(ex.month.isin(months).sum()), "window_months": months,
            "unknown_tap_dropped": [[d["tap_id"], d["reason"]] for d in u.dropped],
            "mixed_refused": _refusal(paths, "mixed_species_log.csv") is not None,
            "mixed_refused_chloramine": _refusal(paths, "mixed_species_log.csv", "chloramine") is not None,
            "multi_month": {"months_in_log": info["n_months_in_log"], "months_used": info["months"],
                            "rows_used": info["n_used"], "rows_outside": info["n_outside"]}}


def part_templates() -> dict:
    """Every template parses: the grab logs with load_log (with the tap map template), the tap map as load_log's tap
    map, the plant log with load_plant, and the grab log as a workbook (openpyxl, same-library)."""
    from .importers import read_grab_log, read_tap_map, xlsx_bytes
    from .pilot import load_log, load_plant
    import csv as _csv
    import datetime as _dt
    tp = os.path.join(TEMPLATES_DIR, "tap_map_by_id_template.csv")
    out = {}
    for f, (kind, dis) in TEMPLATES.items():
        p = os.path.join(TEMPLATES_DIR, f)
        if kind == "log":
            g = read_grab_log(_read_bytes(p), f, _read_bytes(tp), "tap_map_by_id_template.csv", dis)
            out[f] = {"parses_with": "load_log", "rows": len(load_log(p, tp, dis, quiet=True)),
                      "censored": len(g.censored), "dropped": len(g.dropped)}
        elif kind == "taps":
            out[f] = {"parses_with": "load_log (its tap map)", "taps": len(read_tap_map(_read_bytes(p), f)[0])}
        else:
            out[f] = {"parses_with": "load_plant", "months": len(load_plant(p))}
    with open(os.path.join(TEMPLATES_DIR, "grab_log_template.csv"), newline="") as fh:
        rows = list(_csv.reader(fh))
    typed = [rows[0]] + [[_dt.date.fromisoformat(r[0]), _dt.time.fromisoformat(r[1]), r[2],
                          float(r[3]) if r[3][0].isdigit() else r[3]] + r[4:] for r in rows[1:]]
    gx = read_grab_log(xlsx_bytes(typed), "grab_log_template.xlsx", _read_bytes(tp), "tap_map_by_id_template.csv")
    gc = read_grab_log(_read_bytes(os.path.join(TEMPLATES_DIR, "grab_log_template.csv")), "grab_log_template.csv",
                       _read_bytes(tp), "tap_map_by_id_template.csv")
    out["grab_log_template.xlsx"] = {"parses_with": "load_log's reader (read_grab_log)", "rows": len(gx.frame),
                                     "censored": len(gx.censored), "same_as_csv": _rows_of(gx.frame) == _rows_of(gc.frame),
                                     "label": "same-library self-consistency (openpyxl)"}
    return out


def part_review(paths: dict) -> dict:
    """The review of task 16, each case it found as a fixture: text logs with a title, a 'Lab:/Exported' block or
    Excel's 'sep=;' line before the header; a chloride column before the chlorine one; month-name dates and dotted
    times; a quoted comment over two lines (row numbers are lines); readings above the measuring range kept apart from
    readings below detection; a log in µg/L; a LIMS missing code and a reading in the wrong unit; a long-format LIMS
    export (one row per analyte); a log with a free and a total column, with the reading column chosen; a junction
    list (no tap column); an analyzer's typed datetimes from 00:00 in Excel against its CSV twin; a workbook with a
    cover sheet and a stray cell far out; a diurnal grab-sampling study (kept whole); ten analyzers interleaved by
    time (the not-used list stays short and still counts every row); and every reading kept when asked.  The
    expected rows come from each fixture's own table."""
    import contextlib as _cl
    import io as _io
    from . import fixtures as FX
    from .importers import ANALYZER_HOUR, dropped_lines, guess_tap_columns, read_table
    from .pilot import load_log
    out = {}
    want = _expected([(d, int(t[:2]), int(t[3:]), tap, float(v)) for d, t, tap, v in FX.REVIEW_ROWS])
    for f in ("log_title_line.csv", "log_lab_preamble.csv", "log_sep_hint.csv"):
        g = _log(paths, f)
        out[f.split(".")[0]] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want, "how": g.notes[0]}
    g = _log(paths, "log_chloride_first.csv")
    out["log_chloride_first"] = {"reading": g.mapping["reading"], "rows_as_expected": _rows_of(g.frame) == want}
    g = _log(paths, "log_month_names.csv")
    out["log_month_names"] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want,
                              "dropped": len(g.dropped)}
    g = _log(paths, "log_multiline_comment.csv")
    out["log_multiline_comment"] = {"rows": len(g.frame), "dropped": [[d["tap_id"], d["row"]] for d in g.dropped]}
    g = _log(paths, "log_overrange.csv")
    buf = _io.StringIO()
    with _cl.redirect_stdout(buf):
        load_log(paths["log_overrange.csv"], paths["task16_tap_map.csv"])
    kinds = {k: sorted(g.censored.tap_id[g.censored.kind == k]) for k in ("below", "above")}
    out["log_overrange"] = {"rows": len(g.frame), "censored": kinds,
                            "as_written": {"below": sorted(t for _, _, t, v in FX.OVERRANGE_ROWS if v.startswith("<")),
                                           "above": sorted(t for _, _, t, v in FX.OVERRANGE_ROWS if v.startswith(">"))},
                            "pilot_notes": [ln for ln in buf.getvalue().splitlines() if ln.startswith("note: ")
                                            and "kept and listed" in ln]}
    g = _log(paths, "log_ugL.csv")
    out["log_ugL"] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == want,
                      "warned": any("µg/L" in w for w in g.warnings)}
    g = _log(paths, "log_implausible.csv")
    keep = [r for r in FX.IMPLAUSIBLE_ROWS if float(r[3]) <= 10]
    out["log_implausible"] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == _expected(
        [(d, int(t[:2]), int(t[3:]), tap, float(v)) for d, t, tap, v in keep]),
        "dropped": [[d["tap_id"], d["reason"]] for d in g.dropped]}
    long_ = {}
    for dis, col in (("free_chlorine", 3), ("chloramine", 4)):
        g = _log(paths, "lims_long.csv", dis)
        exp = _expected([(d, int(t[:2]), int(t[3:]), tap, float(r[col])) for r in FX.LONG_SAMPLES
                         for d, t, tap in [r[:3]]])
        reasons = {}
        for d in g.dropped:
            reasons[d["reason"]] = reasons.get(d["reason"], 0) + 1
        long_[dis] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == exp, "analyte": g.mapping["analyte"],
                      "reading": g.mapping["reading"], "not_used": reasons,
                      "species_warning": any("does not say free or total" in w for w in g.warnings)}
    out["lims_long"] = long_
    both = {"bare_guess_refused_asks": bool(_refusal(paths, "mixed_species_log.csv")) and "Which column is the reading?"
            in _refusal(paths, "mixed_species_log.csv")}
    for dis, col in (("free_chlorine", "free_chlorine_mgL"), ("chloramine", "total_chlorine_mgL")):
        g = _log(paths, "mixed_species_log.csv", dis, mapping={"reading": col})
        both[dis] = {"reading": col, "y": [float(v) for v in g.frame.y],
                     "other_column_noted": any("not read" in n for n in g.notes)}
    out["both_columns_chosen"] = both
    t = read_table(_read_bytes(paths["junction_list.csv"]), "junction_list.csv")
    out["junction_list"] = guess_tap_columns(t.header)
    gx, gc = _log(paths, "analyzer_midnight.xlsx"), _log(paths, "analyzer_midnight.csv")
    out["analyzer_midnight"] = {"xlsx_readings": gx.analyzer["AN-1"]["readings"],
                                "csv_readings": gc.analyzer["AN-1"]["readings"],
                                "written": 96 * FX.ANALYZER_MIDNIGHT_DAYS, "same_rows": _rows_of(gx.frame) == _rows_of(gc.frame),
                                "no_time_dropped": sum(d["reason"] == "no time of day" for d in gx.dropped),
                                "label": "same-library self-consistency (openpyxl)"}
    t = read_table(_read_bytes(paths["log_cover_and_stray.xlsx"]), "log_cover_and_stray.xlsx")
    out["log_cover_and_stray"] = {"rows": len(t.rows), "header": t.header, "how": t.how, "stray_cell": list(FX.STRAY_CELL),
                                  "label": "same-library self-consistency (openpyxl)"}
    g = _log(paths, "grab_study_24h.csv")
    written = len(FX.STUDY_HOURLY_DAYS) * 24 + 12 + 3
    out["grab_study_24h"] = {"rows": len(g.frame), "written": written, "analyzer_taps": sorted(g.analyzer),
                             "dropped": len(g.dropped)}
    g = _log(paths, "analyzer_interleaved.csv")
    lines = dropped_lines(g.dropped)
    n_in = len(FX.INTERLEAVED_SITES) * FX.INTERLEAVED_DAYS * 96
    out["analyzer_interleaved"] = {"readings_in": n_in, "analyzer_taps": len(g.analyzer), "kept": len(g.frame),
                                   "thinned_listed": sum(d["reason"].startswith("analyzer series at tap") for d in g.dropped),
                                   "not_used_lines": len(lines), "lines_capped": sum("more runs of rows" in ln for ln in lines),
                                   "chars": sum(len(ln) for ln in lines)}
    g = _log(paths, "analyzer_15min.csv", thin_analyzers=False)
    out["analyzer_keep_every_reading"] = {"rows": len(g.frame), "thinned_listed": sum(
        d["reason"].startswith("analyzer series at tap") for d in g.dropped), "noted": any("every reading is kept" in n
                                                                                         for n in g.notes)}
    mm = _log(paths, "multi_month_log.csv")
    _, _, info = mm.window()
    out["window_outside_listed"] = {"n_outside": info["n_outside"], "listed": len(info["outside"]),
                                    "reason": info["outside"][0]["reason"] if info["outside"] else None,
                                    "rows": len(mm.frame), "used": info["n_used"]}
    return out


def part_final_review(paths: dict) -> dict:
    """The final verification of task 16, each case it found as a fixture (expected rows from each fixture's own
    table): the species read inside brackets ('Chlorine (Free)', 'Cl2 (total)') in a long-format export and in two
    wide columns; 'Chloramine (mg/L)' as total chlorine; 'Chlorine Residual' beside 'Free Chlorine' asks; 'Free Cl',
    'Total Cl' and LIMS codes (FCL2, TCL2) pick their rows; analytes that do not pick the rows (unknown codes, two
    chlorine analytes with no species, a 'Parameter Code' column the guess passes over) are refused with the
    question, and read once the analyte value is chosen; µg/L headers without a slash or with a mangled µ, and units
    'mg Cl2/L', 'mgCl2/L', 'mg/litre'; a µg/L log whose header does not say so warns; a comma title before a
    semicolon log, and an Excel title row of two cells; an Excel date with no time among timed ones is no 00:00
    sample; a time column holding a full datetime; a superscript digit in a tap ID; a two-column junction export is
    no tap map; a lab qualifier on every reading and a tap map from another system give one not-used line each; a
    small file of 3 million tiny rows is refused before it is parsed; and the bundled Net2's 08:00 start in the
    pilot."""
    import time as _time
    from . import fixtures as FX
    from .importers import (DROPPED_MAX_LINES, LOG_MAX_LINES, LogRefused, dropped_lines, guess_tap_columns, model_hour,
                            read_grab_log, read_table)
    out = {}
    taps = _taps()

    def exp(rows):          # (date, 'HH:MM', tap, reading) -> frame rows
        return _expected([(d, int(t[:2]), int(t[3:5]), tap, float(v)) for d, t, tap, v in rows])

    long_free = [(d, t, tap, fr) for d, t, tap, fr, _ in FX.LONG_SAMPLES]
    long_total = [(d, t, tap, to) for d, t, tap, _, to in FX.LONG_SAMPLES]
    wide_free = list(FX.REVIEW_ROWS)
    wide_total = [(d, t, tap, f"{float(v) + 0.8:.2f}") for d, t, tap, v in FX.REVIEW_ROWS]
    longs = {}
    for f in ("long_paren.csv", "long_short_names.csv", "long_lims_codes.csv"):
        res = {}
        for dis, want in (("free_chlorine", long_free), ("chloramine", long_total)):
            g = _log(paths, f, dis)
            res[dis] = {"rows": len(g.frame), "rows_as_expected": _rows_of(g.frame) == exp(want),
                        "analyte": g.mapping["analyte"], "reading": g.mapping["reading"],
                        "species_warning": any("does not say free or total" in w for w in g.warnings)}
        longs[f.split(".")[0]] = res
    out["long_species"] = longs
    wide = {}
    for f, dis, col, want in (("wide_paren.csv", "free_chlorine", "Chlorine (Free) mg/L", wide_free),
                              ("wide_paren.csv", "chloramine", "Chlorine (Total) mg/L", wide_total),
                              ("wide_paren_total_first.csv", "free_chlorine", "Chlorine (Free) mg/L", wide_free),
                              ("wide_residual_and_free.csv", "free_chlorine", "Free Chlorine", wide_free)):
        q = _refusal(paths, f, dis)
        g = _log(paths, f, dis, mapping={"reading": col})
        wide[f"{f.split('.')[0]}|{dis}"] = {"bare_guess_refused_asks": bool(q) and "Which column is the reading?" in q,
                                            "chosen": col, "rows_as_expected": _rows_of(g.frame) == exp(want),
                                            "other_column_noted": any("not read" in n for n in g.notes)}
    out["wide_species"] = wide
    q = _refusal(paths, "chloramine_header.csv")
    g = _log(paths, "chloramine_header.csv", "chloramine")
    out["chloramine_header"] = {"free_mode_refused_as_total": bool(q) and "is total chlorine" in q,
                                "chloramine_rows_as_expected": _rows_of(g.frame) == exp(wide_total),
                                "species_warning": any("does not say free or total" in w for w in g.warnings)}
    asks = {}
    for f, value in (("long_unknown_codes.csv", "CHLF"), ("long_two_unnamed.csv", "Chlorine Residual")):
        q = _refusal(paths, f)
        g = _log(paths, f, mapping={"analyte_value": value})
        asks[f.split(".")[0]] = {"refused_asks_value": bool(q) and "Which value is the chlorine reading?" in q,
                                 "value": value, "rows_as_expected": _rows_of(g.frame) == exp(long_free),
                                 "warned": any("does not say free or total" in w for w in g.warnings)}
    q = _refusal(paths, "long_unknown_codes_param.csv", mapping={"reading": "Result"})
    g = _log(paths, "long_unknown_codes_param.csv", mapping={"reading": "Result", "analyte": "Parameter Code",
                                                             "analyte_value": "CHLF"})
    asks["long_unknown_codes_param"] = {"refused_names_the_column": bool(q) and "'Parameter Code' looks like" in q,
                                        "value": "CHLF", "rows_as_expected": _rows_of(g.frame) == exp(long_free)}
    q = _refusal(paths, "long_unknown_codes.csv", mapping={"analyte_value": "TCL2X"})
    asks["value_not_in_column_refused"] = bool(q) and "No row's analyte" in q
    out["analyte_value"] = asks
    units = {}
    for k, h in enumerate(FX.UGL_HEADERS):
        g = _log(paths, f"ugl_header_{k}.csv")
        units[h] = {"y": [round(float(v), 6) for v in g.frame.y], "warned": any("µg/L" in w for w in g.warnings)}
    for k, u in enumerate(FX.UNIT_SPELLINGS):
        g = _log(paths, f"units_column_{k}.csv")
        units[u] = {"rows_as_expected": _rows_of(g.frame) == exp(FX.REVIEW_ROWS), "dropped": len(g.dropped)}
    g = _log(paths, "ugl_unnamed.csv")
    units["unnamed_ugL_warned"] = any("a log in µg/L" in w for w in g.warnings)
    out["units"] = units
    g = _log(paths, "log_comma_title_semicolon.csv")
    gx = _log(paths, "log_two_cell_title.xlsx")
    out["titles"] = {"comma_title_rows_as_expected": _rows_of(g.frame) == exp(FX.REVIEW_ROWS), "comma_title_how": g.notes[0],
                     "xlsx_title_rows_as_expected": _rows_of(gx.frame) == exp(FX.REVIEW_ROWS), "xlsx_title_how": gx.notes[0],
                     "label": "the workbook: same-library self-consistency (openpyxl)"}
    gx, gc = _log(paths, "excel_date_only.xlsx"), _log(paths, "excel_date_only.csv")
    out["excel_date_only"] = {"xlsx_rows": _rows_of(gx.frame), "same_rows_as_csv": _rows_of(gx.frame) == _rows_of(gc.frame),
                              "xlsx_no_time": [[d["tap_id"], d["row"]] for d in gx.dropped if d["reason"] == "no time of day"],
                              "csv_no_time": [[d["tap_id"], d["row"]] for d in gc.dropped if d["reason"] == "no time of day"],
                              "label": "same-library self-consistency (openpyxl)"}
    g = _log(paths, "time_column_epoch_date.csv")
    out["time_column_epoch_date"] = {"rows_as_expected": _rows_of(g.frame) == exp(FX.REVIEW_ROWS), "dropped": len(g.dropped)}
    sup = {}
    for name, fn in (("log", lambda: _log(paths, "superscript_tap_log.csv")),
                     ("tap_map", lambda: read_grab_log(_read_bytes(paths["superscript_tap_log.csv"]), "s.csv",
                                                       _read_bytes(paths["superscript_tap_map.csv"]), "t.csv"))):
        try:
            g = fn()
            sup[name] = {"rows": len(g.frame), "dropped": [[d["tap_id"], d["reason"]] for d in g.dropped]}
        except (LogRefused, ValueError) as e:
            sup[name] = {"error": f"{type(e).__name__}: {e}"}
    out["superscript_tap"] = sup
    t = read_table(_read_bytes(paths["junctions_2col.csv"]), "junctions_2col.csv")
    out["junctions_2col"] = guess_tap_columns(t.header)
    lists = {}
    for f in ("qualifier_interleaved.csv", "unknown_taps_many.csv"):
        g = _log(paths, f)
        lines = dropped_lines(g.dropped)
        lists[f.split(".")[0]] = {"rows_in": g.n_rows, "used": len(g.frame), "not_used": len(g.dropped),
                                  "lines": len(lines), "chars": sum(len(x) for x in lines), "first": lines[0][:140]}
    many = [{"row": 2 + i, "reason": f"reason {i}"} for i in range(DROPPED_MAX_LINES + 5)]
    lines = dropped_lines(many)
    lists["line_cap"] = {"reasons": len(many), "lines": len(lines), "last": lines[-1]}
    out["not_used_lists"] = lists
    t0 = _time.perf_counter()
    try:
        read_table(_read_bytes(paths["csv_tiny_rows.csv"]), "csv_tiny_rows.csv")
        tiny = {"refused": None}
    except LogRefused as e:
        tiny = {"refused": str(e)}
    tiny["within_1s"] = _time.perf_counter() - t0 < 1.0
    tiny["rows_written"] = FX.TINY_ROWS
    tiny["line_limit"] = LOG_MAX_LINES
    out["csv_tiny_rows"] = tiny
    start = json.load(open(READINESS_BUNDLED))["networks"]["Net2"]["flags"]["start_clock_h"]
    g = _log(paths, "net2_clock_log.csv", taps="net2_tap_map.csv")
    out["net2_pilot_hours"] = {"start_clock_h": start, "clock_hours": [int(x // 60) for x in g.clock_min],
                               "app_hours": [int(h) for h in g.frame.hour],
                               "pilot_hours": [model_hour(int(m), start) for m in g.clock_min]}
    return out


# the one real-browser run of task 16 (headless Chromium through Playwright, against 'streamlit run app.py' on
# localhost), 19 October 2026: after the review's fixes it ran at 08:41, before a last edit to importers.py (so not
# on the final tree, as this note first said); since the final verification it ran again at 11:48, on the app's
# final code (app.py, the Get started page and every module they import, none edited after it; only this record's
# writer, the checks and the numbers registry were edited afterwards), with the same results; copied here from
# its run log, which is kept outside the repository; the saved checks do not rerun it (it needs a browser and a server)
BROWSER_CHECK = {"what": "the 'Every junction' dataframe's toolbar CSV download in a real browser, with the default demo "
                         "and again with the demo off after the example grab log and tap map were set on the page's own "
                         "file inputs; Chrome's save dialog stood in by a stub that keeps what the page writes",
                 "when": "2026-10-19 11:48, after the final verification of task 16, on the app's final code "
                         "(app.py, the Get started page and every module they import; the run after the review, at "
                         "08:41, was before a last edit to importers.py)",
                 "toolbar_buttons": ["Show/hide columns", "Download as CSV", "Search", "Fullscreen"],
                 "demo_csv": {"rows": 92, "columns": 6}, "upload_csv": {"rows": 92, "columns": 6},
                 "upload_samples_used": 39, "get_started_title": "Bring your data", "get_started_downloads": 5}


# what the final verification of task 16 measured with its own probe scripts, before its fixes ('before', from its
# logs: the reader cases, the sizes and lists, the superscript tap, the Net2 hours, and AppTest of the app and the Get
# started page, on a copy of the tree it checked), and the same scripts run again on the final tree after the last
# code edit ('after', 19 October 2026); copied here from the logs, which are kept outside the repository; the saved
# checks do not rerun them (the AppTest probes and the timed sizes); timings and memory are run costs of this Mac
FINAL_REVIEW_PROBES = {
    "before": {
        "app_example_log": {"samples_used": 39, "flagged": "35 of 92"},
        "app_example_log_as_long_rows_species_in_brackets": {"samples_used": 78, "flagged": "20 of 92"},
        "app_wide_species_in_brackets_chloramine": {"samples_used": 39, "flagged": "87 of 92",
                                                    "reading_column": "Chlorine (Free) mg/L"},
        "reader_long_species_in_brackets_y": [0.45, 1.2, 0.38, 1.05],
        "reader_chloramine_header_free_mode_y": [2.1, 1.9],
        "reader_free_cl_total_cl_ph_no_units_y": [0.45, 1.2, 7.6, 0.38, 1.05, 7.8],
        "csv_tiny_rows": {"1000000": {"peak_mb": 454, "seconds": 2.05},
                          "3000000": {"peak_mb": 1036, "seconds": 7.66, "refused": True},
                          "6000000": {"peak_mb": 1930, "seconds": 443.54, "refused": True}},
        "not_used_lists": {"qualifier_interleaved": {"rows": 57600, "lines": 1451, "chars": 324067},
                           "unknown_taps_interleaved": {"rows": 60000, "lines": 2000, "chars": 476796}},
        "check_a_file_junctions_2col": "As a tap map: 4 taps, tap column 'Node ID', junction column 'Elevation'",
        "superscript_tap": "ValueError: invalid literal for int() with base 10: '1²'",
        "net2_upload": {"app_hours": [9, 10, 11, 12, 13, 14, 15], "pilot_hours": [1, 2, 3, 4, 5, 6, 7],
                        "caption_says_unconverted": False}},
    "after": {
        "app_example_log": {"samples_used": 39, "flagged": "35 of 92"},
        "app_example_log_as_long_rows_species_in_brackets": {"samples_used": 39, "flagged": "35 of 92"},
        "app_wide_species_in_brackets_chloramine": {"refused_asks": True, "reading_column": "(choose a column)"},
        "reader_long_species_in_brackets_y": [0.45, 0.38],
        "reader_chloramine_header_free_mode": "refused: total chlorine in the free-chlorine model",
        "reader_free_cl_total_cl_ph_no_units_y": [0.45, 0.38],
        "csv_tiny_rows": {"1000000": {"peak_mb": 388, "seconds": 1.15},
                          "2000000": {"peak_mb": 247, "seconds": 0.59, "refused": True},
                          "3000000": {"peak_mb": 56, "seconds": 0.01, "refused": True},
                          "6000000": {"peak_mb": 74, "seconds": 0.02, "refused": True},
                          "33000000": {"peak_mb": 236, "seconds": 0.11, "refused": True}},
        "not_used_lists": {"qualifier_interleaved": {"rows": 57600, "lines": 1, "chars": 102},
                           "unknown_taps_interleaved": {"rows": 60000, "lines": 1, "chars": 131}},
        "check_a_file_junctions_2col": "As a table: it has a junction column ('Node ID') but no tap column",
        "superscript_tap": "read; the tap is listed as not in the tap map",
        "net2_upload": {"app_hours": [9, 10, 11, 12, 13, 14, 15], "pilot_hours": [1, 2, 3, 4, 5, 6, 7],
                        "caption_says_unconverted": True}}}


def samples_record() -> dict:
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t16_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.SAMPLE_FIXTURES)
        return {"about": "Task 16's grab-log and clock-time tests on files built by residualmap/fixtures.py (invented "
                         "readings at Net3's junctions; grab_log_typed.xlsx and the Excel template written by openpyxl, "
                         "the library that reads them, so those two are same-library self-consistency). Written by "
                         "python -m residualmap.intake_exp samples; the task16 checks rerun each part and compare. "
                         "The app's own numbers (AppTest) are in outputs/chem/checks_report.json. Part 'review' holds "
                         "the review's cases and part 'final_review' the final verification's; parts "
                         "'final_review_probes' and 'browser_check' are copied from the final verification's probe "
                         "logs and the one real-browser run's log, and the checks do not rerun them.",
                "mapping": part_mapping(paths), "encodings": part_encodings(paths), "excel": part_excel(paths),
                "analyzer": part_analyzer(paths), "utc": part_utc(paths), "clock_helpers": part_clock_helpers(),
                "app_inputs": part_app_inputs(paths), "templates": part_templates(), "review": part_review(paths),
                "final_review": part_final_review(paths), "final_review_probes": FINAL_REVIEW_PROBES,
                "browser_check": BROWSER_CHECK}


# ----------------------------------------------------------------------------- task 16b: taps by coordinates
RECORD_16B = os.path.join(OUT, "tap_snapping.json")
SNAP_TAPS = 2000                       # taps per network, uniform by length on the pipes they may be on
SNAP_NOISE_M = (5.0, 10.0, 20.0)       # radial RMS GPS noise; per axis r / sqrt(2): 3.54, 7.07 and 14.14 m
SNAP_POS_SEED = 0                      # tap positions (numpy default_rng)
SNAP_NOISE_SEEDS = (1, 2, 3)           # the noise at 5, 10 and 20 m
URBAN_MAIN_TAPS, URBAN_MAIN_SEED = 100, 10   # taps on the transmission main, tagged on_main yes (noise seeds 11 to 13)
# Pre-registered before the scored run (plan task 16b; the scratch note prereg.txt of the build holds the same words
# and its time).  The bar is as written (lead decision 7).  Open decision 11, how a right junction shown only as the
# flagged second choice counts, was still open for Devansh: the run counted strictly, with the lenient share reported
# beside it, and both miss 0.98, so the choice between decision 11 and a lower bar is still his.
SNAP_BAR = {"network": "ky4", "noise_m": 5.0, "share_at_least": 0.98,
            "counting": "strict: the suggested junction_id is the right junction (the rule applied to the tap's true "
                        "pipe); a right junction shown only as the flagged second choice does not count (open decision "
                        "11; that share is reported beside it, not scored)",
            "taps": "2,000 uniform by length on ky4's pipes under 16 in that have a junction at one end, "
                    "default_rng(0); noise Gaussian, r / sqrt(2) per axis, default_rng(1), (2), (3) for 5, 10, 20 m",
            "rule": "nearest pipe under 16 in, segments and vertices in metres; its downstream end by the sign of the "
                    "mean simulated flow between 07:00 and 17:00 of the scored day; no flow (every hour at most 1e-7 "
                    "m3/s) takes the end with the older water; a downstream tank or reservoir takes the other end; "
                    "second choice the nearest pipe with another junction, shown within 1.5 times the distance"}
DUAL_PIPES_PREFIX = ("DM", "EM")       # the urban grid's second mains, and the avenue pipes they run beside
FOUR_IN_PREFIX = "H"                   # the 4 in line: H<i>05


def _r(x, k: int = 4):
    return None if x is None else round(float(x), k)


def _seg_pipe(geom):
    n = len(geom.seg_a)
    starts = list(geom.seg_start) + [n]
    return np.repeat(np.arange(len(geom.pipes)), np.diff(starts))


def taps_on_pipes(geom, eligible, n: int, seed: int):
    """n points uniform by length on the pipes where eligible (a mask over geom.pipes) is True: (points in metres
    around geom.origin, each point's pipe index)."""
    sp = _seg_pipe(geom)
    seg_len = np.hypot(*(geom.seg_b - geom.seg_a).T)
    w = np.where(eligible[sp], seg_len, 0.0)
    rng = np.random.default_rng(seed)
    k = rng.choice(len(w), size=n, p=w / w.sum())
    t = rng.uniform(size=n)[:, None]
    return geom.seg_a[k] + t * (geom.seg_b[k] - geom.seg_a[k]), sp[k]


def _coords_of(geom, P, on_main=None, prefix: str = "T"):
    """A TapCoords (x and y in the model's own units) for points P in metres around geom.origin."""
    from .importers import TapCoords
    xy = np.asarray(P) / geom.metres_per_unit + np.array(geom.origin)
    om = [False] * len(xy) if on_main is None else list(on_main)
    return TapCoords(kind="xy", taps=[{"tap_id": f"{prefix}{i + 1:04d}", "row": i + 2, "a": float(x), "b": float(y),
                                       "junction": "", "on_main": bool(m)} for i, ((x, y), m) in enumerate(zip(xy, om))],
                     n_rows=len(xy))


def _noisy(P, r: float, seed: int):
    return P + np.random.default_rng(seed).normal(0.0, r / math.sqrt(2.0), size=np.shape(P))


def _score(sug, geom, rules, true_pipe, mask=None) -> dict:
    """Shares over the taps (mask: a subset): the right junction (strict), the right junction or a right flagged
    second choice, the right pipe; flagged two pipes close, far, flow reversing, no flow; snap distances."""
    rows = sug.rows
    m = np.ones(len(rows), bool) if mask is None else np.asarray(mask, bool)
    tj = np.array([rules[p][0] for p in true_pipe], object)
    gj = np.array([r["junction_id"] for r in rows], object)
    sj = np.array([r["second_choice"] or "" for r in rows], object)
    gp = np.array([r["pipe"] or "" for r in rows], object)
    two = np.array([any(f.startswith("two pipes close") for f in r["flags"]) for r in rows])
    ok = gj == tj
    sel = lambda a: float(np.mean(a[m])) if m.any() else None
    snap = np.array([r["snap_m"] if r["snap_m"] is not None else np.nan for r in rows])
    return {"taps": int(m.sum()), "right_junction": _r(sel(ok)), "right_junction_or_flagged_second": _r(sel(ok | (two & (sj == tj)))),
            "right_pipe": _r(sel(gp == np.array([geom.pipes[p] for p in true_pipe], object))),
            "flagged_two_pipes_close": _r(sel(two)),
            "flagged_far": _r(sel(np.array([any(f.startswith("far:") for f in r["flags"]) for r in rows]))),
            "flagged_flow_reverses": _r(sel(np.array([any("reverses" in f for f in r["flags"]) for r in rows]))),
            "flagged_no_flow": _r(sel(np.array([any(f.startswith("no flow") for f in r["flags"]) for r in rows]))),
            "snap_m_median": _r(np.nanmedian(snap[m]), 2), "snap_m_p95": _r(np.nanpercentile(snap[m], 95), 2),
            "ambiguous_rows_list_two_candidates": bool(all(r["second_choice"] and r["second_m"] is not None
                                                           for r, t, k in zip(rows, two, m) if t and k))}


def _assignments(sug) -> list:
    return [[r["junction_id"], r["pipe"], r["second_choice"] or "", _r(r["snap_m"], 2)] for r in sug.rows]


@functools.lru_cache(maxsize=12)
def _geom_hyd(path: str):
    from .importers import pipe_junctions, tap_geometry, tap_hydraulics
    g = tap_geometry(path)
    h = tap_hydraulics(path)
    return g, h, pipe_junctions(g, h)


def part_snap_ky4() -> dict:
    """tap_snapping_ky4: 2,000 taps on ky4's pipes under 16 in, at 5, 10 and 20 m of radial RMS noise, given in the
    model's x and y (feet); the scored run of SNAP_BAR."""
    from .importers import suggest_from
    g, h, rules = _geom_hyd("ky4")
    eligible = (g.diam_in < 16.0) & np.array([j is not None for j, _, _ in rules])
    P, tp = taps_on_pipes(g, eligible, SNAP_TAPS, SNAP_POS_SEED)
    out = {"network": "ky4", "units": g.scale["units"], "pipes": len(g.pipes), "segments": len(g.seg_a),
           "pipes_16in_or_more": int((g.diam_in >= 16.0).sum()), "taps": SNAP_TAPS,
           "rule_codes": {c: sum(c in cs for _, _, cs in rules) for c in ("reverses", "no_flow", "end_not_junction")},
           "by_noise": {}}
    for r, seed in zip(SNAP_NOISE_M, SNAP_NOISE_SEEDS):
        sug = suggest_from(_coords_of(g, _noisy(P, r, seed)), g, h)
        out["by_noise"][f"{r:g}"] = _score(sug, g, rules, tp)
    got = out["by_noise"]["5"]["right_junction"]
    out["bar"] = dict(SNAP_BAR, result=got, passed=bool(got >= SNAP_BAR["share_at_least"]))
    return out


def part_snap_urban(paths: dict) -> dict:
    """tap_snapping_urban_grid: 2,000 taps on the synthetic urban grid's pipes under 16 in, and 100 taps on its 24 in
    transmission main tagged on_main yes, at 5, 10 and 20 m; the dual mains and the 4 in line beside the transmission
    main scored apart; no untagged tap may be put on the transmission main."""
    from . import fixtures as FX
    from .importers import suggest_from
    g, h, rules = _geom_hyd(paths["synthetic_urban_grid.inp"])
    tm = set(FX._urban_parts()[2]["transmission"])
    has_j = np.array([j is not None for j, _, _ in rules])
    P, tp = taps_on_pipes(g, (g.diam_in < 16.0) & has_j, SNAP_TAPS, SNAP_POS_SEED)
    inside = np.array([p in tm and p != "TM0" for p in g.pipes]) & has_j       # the main inside the grid, TM1 to TM4
    Pm, tpm = taps_on_pipes(g, inside, URBAN_MAIN_TAPS, URBAN_MAIN_SEED)
    names = [g.pipes[p] for p in tp]
    n, nm = len(P), len(Pm)
    small = np.r_[np.ones(n, bool), np.zeros(nm, bool)]
    dual = np.r_[[x.startswith(DUAL_PIPES_PREFIX) or (x.startswith("H") and x.endswith("03")) or x.startswith("V07")
                  for x in names], np.zeros(nm, bool)]
    four = np.r_[[x.startswith(FOUR_IN_PREFIX) and x.endswith("05") for x in names], np.zeros(nm, bool)]
    out = {"network": "synthetic urban grid", "units": g.scale["units"], "pipes": len(g.pipes),
           "segments": len(g.seg_a), "taps": n, "dual_main_taps": int(dual.sum()), "four_inch_taps": int(four.sum()),
           "on_main_taps": nm, "by_noise": {}}
    for k, (r, seed) in enumerate(zip(SNAP_NOISE_M, SNAP_NOISE_SEEDS)):
        Q = np.vstack([_noisy(P, r, seed), _noisy(Pm, r, URBAN_MAIN_SEED + 1 + k)])
        sug = suggest_from(_coords_of(g, Q, [False] * n + [True] * nm), g, h)
        tp_all = np.r_[tp, tpm]
        on_tm = np.array([(row["pipe"] or "") in tm for row in sug.rows])
        right = np.array([row["pipe"] == g.pipes[p] for row, p in zip(sug.rows, tp_all)])
        out["by_noise"][f"{r:g}"] = {
            "small_mains": _score(sug, g, rules, tp_all, small), "dual_mains": _score(sug, g, rules, tp_all, dual),
            "four_inch_line": _score(sug, g, rules, tp_all, four),
            "untagged_taps_put_on_the_transmission_main": int(on_tm[small].sum()),
            "four_inch_taps_put_on_the_transmission_main": int(on_tm[four].sum()),
            "on_main_tagged": {"taps": nm, "put_on_the_transmission_main": _r(on_tm[~small].mean()),
                               "right_pipe": _r(right[~small].mean())}}
    out["transmission_main_variants"] = _urban_tm_variants(paths, P, Pm, tp, tpm, small, four, tm)
    from .importers import MAIN_BAND_RATIO, MAIN_CUT_RATIO, MAIN_MIN_IN
    out["main_cut_in"], out["main_band_from_in"] = _r(MAIN_MIN_IN * MAIN_CUT_RATIO, 3), _r(MAIN_MIN_IN * MAIN_BAND_RATIO, 3)
    return out


def _urban_tm_variants(paths: dict, P, Pm, tp, tpm, small, four, tm) -> dict:
    """The review of task 16b: the same taps on the urban grid with its transmission main a metric DN400 (400 mm,
    15.75 in) and a 16 in main by its internal diameter (388.6 mm, 15.3 in), with today's cut (MAIN_CUT_RATIO of 16 in,
    15.0 in on the model's diameter) and with the cut at 16.0 in as first built: untagged taps put on the
    transmission main, and the 4 in line's taps among them.  The final verification added two 16 in plastic mains by
    internal diameter, under the cut (C905 DR14, 378.8 mm; PE100 SDR11, 327.4 mm): with today's cut, the untagged
    taps put on the main with no flag naming a main (the band below the cut, MAIN_BAND_RATIO, flags them), and those
    whose flag names their true pipe as the nearest smaller one."""
    from . import fixtures as FX
    from .importers import MAIN_CUT_RATIO, MAIN_MIN_IN, suggest_from
    n, nm = len(P), len(Pm)
    out = {}
    for name, mm in FX.URBAN_TM_VARIANTS_MM.items():
        g, h, rules = _geom_hyd(paths[name])
        z = {"main_mm": mm, "main_in_model": _r(g.diam_in[g.pipes.index("TM1")], 3), "by_noise": {}}
        for k, (r, seed) in enumerate(zip(SNAP_NOISE_M, SNAP_NOISE_SEEDS)):
            Q = np.vstack([_noisy(P, r, seed), _noisy(Pm, r, URBAN_MAIN_SEED + 1 + k)])
            c = _coords_of(g, Q, [False] * n + [True] * nm)
            row = {}
            for label, mmin in (("cut_15in", MAIN_MIN_IN), ("cut_16in_as_first_built", MAIN_MIN_IN / MAIN_CUT_RATIO)):
                sug = suggest_from(c, g, h, mmin)
                on_tm = np.array([(x["pipe"] or "") in tm for x in sug.rows])
                row[label] = {"untagged_taps_put_on_the_transmission_main": int(on_tm[small].sum()),
                              "four_inch_taps_put_on_the_transmission_main": int(on_tm[four].sum()),
                              "small_mains_right_junction": _score(sug, g, rules, np.r_[tp, tpm], small)["right_junction"]}
                if label == "cut_15in" and z["main_in_model"] < MAIN_MIN_IN * MAIN_CUT_RATIO:
                    tp_all = np.r_[tp, tpm]
                    on = [i for i in np.flatnonzero(on_tm & small)]
                    named = [re.search(r"the nearest smaller pipe is (\S+) ", " ".join(sug.rows[i]["flags"])) for i in on]
                    row[label].update({
                        "untagged_taps_on_the_main_with_no_flag_naming_a_main":
                            sum(not any("main" in f for f in sug.rows[i]["flags"]) for i in on),
                        "untagged_taps_on_the_main_flagged_under_the_cut":
                            sum(any(f.startswith("a pipe of ") and "under the" in f for f in sug.rows[i]["flags"]) for i in on),
                        "flag_names_the_true_pipe_as_the_nearest_smaller_one":
                            sum(bool(m) and m.group(1) == g.pipes[tp_all[i]] for i, m in zip(on, named))})
            z["by_noise"][f"{r:g}"] = row
        out[name] = z
    return out


def part_snap_direction(paths: dict) -> dict:
    """tap_snapping_direction_invariant: the 2,000 ky4 taps at 5 m on ky4 and on ky4_half_reversed.inp (half the pipes
    drawn the other way, vertices reversed): the same assignments, and the flow sign flipped on every flipped pipe
    that has flow."""
    from . import fixtures as FX
    from .importers import suggest_from
    g1, h1, r1 = _geom_hyd("ky4")
    g2, h2, r2 = _geom_hyd(paths["ky4_half_reversed.inp"])
    eligible = (g1.diam_in < 16.0) & np.array([j is not None for j, _, _ in r1])
    P, _ = taps_on_pipes(g1, eligible, SNAP_TAPS, SNAP_POS_SEED)
    Q = _noisy(P, SNAP_NOISE_M[0], SNAP_NOISE_SEEDS[0])
    a1 = _assignments(suggest_from(_coords_of(g1, Q), g1, h1))
    a2 = _assignments(suggest_from(_coords_of(g2, Q), g2, h2))
    flipped = FX.ky4_flipped_pipes()
    with_flow = [p for p in flipped if max(abs(x) for x in h1["flow"][p]) > 1e-7]
    sign_flipped = sum(all((x > 0) != (y > 0) for x, y in zip(h1["flow"][p], h2["flow"][p]) if abs(x) > 1e-7)
                       for p in with_flow)
    return {"taps": len(a1), "pipes_flipped": len(flipped), "flipped_with_flow": len(with_flow),
            "flow_sign_flipped_on": sign_flipped, "same_origin": list(g1.origin) == list(g2.origin),
            "assignments_differ": sum(x != y for x, y in zip(a1, a2)),
            "junction_rules_differ": sum(x != y for x, y in zip(r1, r2))}


def part_snap_schematic(paths: dict) -> dict:
    """tap_snapping_refuses_schematic: Net3 (a schematic drawing) with taps in x and y and in latitude and longitude,
    and Net3 without coordinates, are refused with 'use junction IDs'.  No accuracy is computed on Net3."""
    from .importers import LogRefused, suggest_tap_map
    out = {}
    for model, f in (("Net3", "net3_taps_xy.csv"), ("Net3", "net3_taps_latlon.csv"),
                     (paths["net3_no_coords.inp"], "net3_taps_xy.csv")):
        key = f"{os.path.basename(model)} with {f}"
        try:
            suggest_tap_map(_read_bytes(paths[f]), f, model)
            out[key] = {"refused": False}
        except LogRefused as e:
            out[key] = {"refused": True, "message": str(e), "says_use_junction_ids": "Use junction IDs" in str(e)}
    return out


def _lcc_table(path: str) -> dict:
    import csv as _csv
    with open(path, newline="") as fh:
        rows = list(_csv.DictReader(fh))
    return {"tap": [r["tap_id"] for r in rows], "role": [r["role"] for r in rows],
            "lat": np.array([float(r["lat"]) for r in rows]), "lon": np.array([float(r["lon"]) for r in rows]),
            "xy": np.array([[float(r["x_m"]), float(r["y_m"])] for r in rows])}


def part_georef(paths: dict) -> dict:
    """georef_enu_affine: the ENU-then-affine georeference against the ellipsoidal Lambert conformal conic truth of
    fixtures.py, over 20 and 50 km zones, with 3 control taps (exact, unchecked) and 6 (checked, leave-one-out
    residuals); the direct affine from latitude and longitude, kept only here, for the error the ENU step removes; a
    control tap moved by 25 m is flagged by its leave-one-out residual; nearly collinear control taps are refused.
    Since the final verification the plane's tangent point is the taps' median position, not their mean: each zone
    also gives the maximum error with the plane at the mean, as first built and as the bar was scored, and the
    moved control tap's leave-one-out residuals (now in closed form) are compared with one refit per control tap."""
    from .importers import GeorefRefused, LOO_FLAG_M, control_spread, enu_from_latlon, fit_affine, apply_affine, georeference
    out = {}
    for km in (20, 50):
        t = _lcc_table(paths[f"lcc_taps_{km}km.csv"])
        ctrl_all = [i for i, r in enumerate(t["role"]) if r == "control"]
        test = np.array([i for i, r in enumerate(t["role"]) if r == "test"])
        z = {"test_taps": len(test)}
        for k in (3, 6):
            ctrl = ctrl_all[:k]
            xy, info = georeference(t["lat"], t["lon"], ctrl, t["xy"][ctrl], 1.0, t["tap"])
            err = np.hypot(*(xy[test] - t["xy"][test]).T)
            xm, _ = georeference(t["lat"], t["lon"], ctrl, t["xy"][ctrl], 1.0, t["tap"], centre="mean")
            M = fit_affine(np.column_stack([t["lat"][ctrl], t["lon"][ctrl]]), t["xy"][ctrl])
            derr = np.hypot(*(apply_affine(M, np.column_stack([t["lat"][test], t["lon"][test]])) - t["xy"][test]).T)
            z[f"{k}_controls"] = {"max_error_m": _r(err.max()), "median_error_m": _r(np.median(err)),
                                  "max_error_m_plane_at_the_mean_as_first_built": _r(np.hypot(*(xm[test] - t["xy"][test]).T).max()),
                                  "checked": info["checked"], "loo_m": {x["tap_id"]: _r(x["residual_m"]) for x in info["loo"]},
                                  "direct_latlon_affine_max_error_m": _r(derr.max()),
                                  "direct_latlon_affine_median_error_m": _r(np.median(derr))}
        out[f"{km}km"] = z
        if km == 20:            # a control tap whose junction is 25 m from where its position was taken
            moved = t["xy"][ctrl_all].copy()
            moved[3] += np.array([25.0, 0.0])
            _, info = georeference(t["lat"], t["lon"], ctrl_all, moved, 1.0, t["tap"])
            refit = _loo_by_refit(t["lat"], t["lon"], ctrl_all, moved)
            out["moved_control"] = {"moved": t["tap"][ctrl_all[3]], "by_m": 25.0, "flagged": info["flagged"],
                                    "loo_m": {x["tap_id"]: _r(x["residual_m"]) for x in info["loo"]},
                                    "flag_above_m": LOO_FLAG_M,
                                    "closed_form_against_a_refit_per_control_tap_max_abs_diff_m":
                                        _r(max(abs(x["residual_m"] - y) for x, y in zip(info["loo"], refit)
                                               if x["residual_m"] is not None), 9)}
    out["bar"] = {"what": "maximum error over the test taps with 3 control taps, both zones", "at_most_m": 1.0,
                  "result_m": max(out["20km"]["3_controls"]["max_error_m"], out["50km"]["3_controls"]["max_error_m"]),
                  "result_m_as_scored": max(out["20km"]["3_controls"]["max_error_m_plane_at_the_mean_as_first_built"],
                                            out["50km"]["3_controls"]["max_error_m_plane_at_the_mean_as_first_built"]),
                  "plane": "tangent at the taps' median position since the final verification (the bar was scored with "
                           "it at their mean: result_m_as_scored)"}
    out["bar"]["passed"] = out["bar"]["result_m"] <= out["bar"]["at_most_m"] and \
        out["bar"]["result_m_as_scored"] <= out["bar"]["at_most_m"]
    c = _lcc_table(paths["collinear_controls.csv"])
    ctrl = [i for i, r in enumerate(c["role"]) if r == "control"]
    enu = enu_from_latlon(c["lat"], c["lon"], float(c["lat"].mean()), float(c["lon"].mean()))
    sp = control_spread(enu[ctrl])
    try:
        georeference(c["lat"], c["lon"], ctrl, c["xy"][ctrl], 1.0, c["tap"])
        out["collinear"] = {"refused": False}
    except GeorefRefused as e:
        out["collinear"] = {"refused": True, "message": str(e), "cond": _r(sp["cond"], 1), "area": _r(sp["area"], 6)}
    out["review"] = _georef_review(paths)
    return out


def _loo_by_refit(lat, lon, controls, control_xy) -> list:
    """Each control tap's leave-one-out residual in metres by one refit of the affine on the other control taps, as
    first built (the check of georeference's closed form), on the same plane."""
    from .importers import apply_affine, enu_from_latlon, fit_affine, plane_centre
    enu = enu_from_latlon(lat, lon, *plane_centre(lat, lon))
    ce, cxy = enu[controls], np.asarray(control_xy, float)
    c0 = cxy.mean(0)
    out = []
    for i in range(len(controls)):
        rest = [j for j in range(len(controls)) if j != i]
        Mi = fit_affine(ce[rest], cxy[rest] - c0)
        out.append(float(np.hypot(*(apply_affine(Mi, ce[[i]])[0] + c0 - cxy[i]))))
    return out


TYPO_TAPS = {"latitude 4.001 for 40.01": (4.001, -100.6), "longitude -10.06 for -100.6": (40.01, -10.06),
             "latitude -40.01 for 40.01": (-40.01, -100.6)}    # the final verification's typo taps


def _georef_final(paths: dict) -> dict:
    """The final verification of task 16b, georeference cases: one typo tap in a small tap map (the control taps, one
    test tap and the typo tap, for each of the 400 test taps) with the plane at the taps' median position and at their
    mean as first built; control taps along one avenue for a 4 km square zone (lcc_strip_controls.csv), the taps
    whose predicted error is above FIT_PREDICTED_FLAG_M flagged; and control taps taken at one place and at two."""
    from .importers import (FIT_GPS_ASSUMED_M, FIT_PREDICTED_FLAG_M, LOO_LEVERAGE_MAX, SAME_PLACE_M, GeorefRefused,
                            georeference)
    from . import fixtures as FX
    out = {"settings": {"same_place_m": SAME_PLACE_M, "loo_leverage_max": LOO_LEVERAGE_MAX,
                        "gps_assumed_m": FIT_GPS_ASSUMED_M, "predicted_flag_above_m": FIT_PREDICTED_FLAG_M},
           "typo_in_a_small_map": {}}
    for km in (20, 50):
        t = _lcc_table(paths[f"lcc_taps_{km}km.csv"])
        ctrl_all = [i for i, r in enumerate(t["role"]) if r == "control"]
        test = [i for i, r in enumerate(t["role"]) if r == "test"]
        for k in (4, 6):
            ctrl = ctrl_all[:k]
            z = {}
            for label, (la, lo) in TYPO_TAPS.items():
                row = {}
                for how in ("median", "mean"):
                    errs, loo = [], []
                    for i in test:
                        lat, lon = np.r_[t["lat"][ctrl + [i]], la], np.r_[t["lon"][ctrl + [i]], lo]
                        xy, info = georeference(lat, lon, list(range(k)), t["xy"][ctrl], 1.0, centre=how)
                        errs.append(float(np.hypot(*(xy[k] - t["xy"][i]))))
                        loo.append(info["loo_max_m"])
                    row[how] = {"max_error_m": _r(max(errs), 3),
                                "loo_max_m": _r(max((x for x in loo if x is not None), default=None), 3)}
                z[label] = {"plane_at_the_median": row["median"], "plane_at_the_mean_as_first_built": row["mean"]}
            out["typo_in_a_small_map"][f"{km}km_{k}_controls"] = z
    st = _lcc_table(paths["lcc_strip_controls.csv"])
    ctrl = [i for i, r in enumerate(st["role"]) if r == "control"]
    test = np.array([i for i, r in enumerate(st["role"]) if r == "test"])
    xy, info = georeference(st["lat"], st["lon"], ctrl, st["xy"][ctrl], 1.0, st["tap"])
    err = np.hypot(*(xy[test] - st["xy"][test]).T)
    flagged = np.array(info["predicted_m"])[test] > FIT_PREDICTED_FLAG_M
    out["strip_controls"] = {"controls": len(ctrl), "control_gps_error_m": FX.STRIP_GPS_M, "zone_m": FX.STRIP_SIDE_M,
                             "accepted": True, "narrow_strip": info["narrow_strip"], "across_m": _r(info["across_m"], 1),
                             "cond": _r(info["cond"], 1), "test_taps": len(test), "max_error_m": _r(err.max(), 2),
                             "over_10m": int((err > 10).sum()), "over_30m": int((err > 30).sum()),
                             "flagged_uncertain": int(flagged.sum()),
                             "over_10m_not_flagged": int(((err > 10) & ~flagged).sum()),
                             "over_30m_not_flagged": int(((err > 30) & ~flagged).sum()),
                             "loo_max_m": _r(info["loo_max_m"], 2), "control_flagged_by_loo": info["flagged"],
                             "flag_above_m": FIT_PREDICTED_FLAG_M}
    x0, y0 = FX.lcc_forward(*FX.LCC_CENTRE)
    for label, pts in (("one_place", [(0.0, 0.0)] * 4), ("two_places", [(0.0, 0.0), (0.0, 0.0), (500.0, 0.0), (500.0, 0.0)]),
                       ("within_4m", [(0.0, 0.0), (2.0, 1.0), (1.0, 3.0), (3.0, 2.0)])):
        ll = [FX.lcc_inverse(x0 + e, y0 + n) for e, n in pts + [(100.0, 100.0)]]
        try:
            georeference(np.array([p[0] for p in ll]), np.array([p[1] for p in ll]), list(range(len(pts))),
                         np.array(pts), 1.0)
            out[label] = {"refused": False}
        except GeorefRefused as e:
            out[label] = {"refused": True, "message": str(e)}
    return out


def _georef_review(paths: dict) -> dict:
    """The review of task 16b, georeference cases: the projection against Snyder's worked example; the 20 km zone
    moved to the antimeridian, its longitudes written from -180 to 180; a narrow valley zone whose
    control taps carry 3 m of GPS error, accepted with a warning (8 km by 300 m: the relative tests alone, as first
    built, refused it); four control taps that fail only the triangle test,
    refused with that test alone named; and the spread measures of 4,000 control taps (default_rng(1), a 10 km
    square), which the pair loop took seconds for."""
    from . import fixtures as FX
    from .importers import GeorefRefused, apply_affine, control_spread, enu_from_latlon, fit_affine, georeference
    out = {}
    sx, sy = FX.lcc_forward(*FX.SNYDER_LCC["point"], FX.SNYDER_LCC["zone"])
    out["lcc_against_snyder"] = {"x_m": _r(sx, 2), "y_m": _r(sy, 2), "snyder_x_m": FX.SNYDER_LCC["xy_m"][0],
                                 "snyder_y_m": FX.SNYDER_LCC["xy_m"][1],
                                 "max_abs_diff_m": _r(max(abs(sx - FX.SNYDER_LCC["xy_m"][0]),
                                                          abs(sy - FX.SNYDER_LCC["xy_m"][1])), 3)}
    a, b = _lcc_table(paths["lcc_taps_20km.csv"]), _lcc_table(paths["lcc_antimeridian.csv"])
    test = np.array([i for i, r in enumerate(a["role"]) if r == "test"])
    am = {"lon_min": _r(b["lon"].min(), 3), "lon_max": _r(b["lon"].max(), 3)}
    for k in (3, 6):
        ctrl = [i for i, r in enumerate(a["role"]) if r == "control"][:k]
        xa, _ = georeference(a["lat"], a["lon"], ctrl, a["xy"][ctrl], 1.0, a["tap"])
        xb, ib = georeference(b["lat"], b["lon"], ctrl, b["xy"][ctrl], 1.0, b["tap"])
        # as first built: the plane's centre at the arithmetic mean of the longitudes as written
        e0 = enu_from_latlon(b["lat"], b["lon"], float(b["lat"].mean()), float(b["lon"].mean()))
        c0 = b["xy"][ctrl].mean(0)
        x0 = apply_affine(fit_affine(e0[ctrl], b["xy"][ctrl] - c0), e0) + c0
        am[f"{k}_controls"] = {"max_error_m": _r(np.hypot(*(xb[test] - b["xy"][test]).T).max()),
                               "zone_at_100w_max_error_m": _r(np.hypot(*(xa[test] - a["xy"][test]).T).max()),
                               "flagged": ib["flagged"],
                               "arithmetic_mean_longitude_max_error_m": _r(np.hypot(*(x0[test] - b["xy"][test]).T).max())}
    out["antimeridian"] = am
    nz = _lcc_table(paths["lcc_narrow_zone.csv"])
    ctrl = [i for i, r in enumerate(nz["role"]) if r == "control"]
    test = np.array([i for i, r in enumerate(nz["role"]) if r == "test"])
    xy, info = georeference(nz["lat"], nz["lon"], ctrl, nz["xy"][ctrl], 1.0, nz["tap"])
    err = np.hypot(*(xy[test] - nz["xy"][test]).T)
    out["narrow_zone"] = {"zone_m": list(FX.NARROW_ZONE_M), "control_gps_error_m": FX.NARROW_CONTROL_NOISE_M,
                          "accepted": True, "relative_tests_alone_refuse": info["narrow_strip"],
                          "warned": info["narrow_strip"], "cond": _r(info["cond"], 1),
                          "area": _r(info["area"], 4), "across_m": _r(info["across_m"], 1),
                          "test_taps": len(test), "max_error_m": _r(err.max(), 2), "median_error_m": _r(np.median(err), 2)}
    tr = _lcc_table(paths["triangle_only_controls.csv"])
    ctrl = [i for i, r in enumerate(tr["role"]) if r == "control"]
    try:
        georeference(tr["lat"], tr["lon"], ctrl, tr["xy"][ctrl], 1.0, tr["tap"])
        out["triangle_only"] = {"refused": False}
    except GeorefRefused as e:
        out["triangle_only"] = {"refused": True, "message": str(e),
                                "names_the_triangle_test_only": "largest triangle" in str(e) and "longer than wide" not in str(e)}
    E = np.random.default_rng(1).uniform(-5000, 5000, size=(4000, 2))
    sp = control_spread(E)
    out["many_controls"] = {"controls": len(E), "cond": _r(sp["cond"], 4), "area": _r(sp["area"], 4),
                            "largest_distance_m": _r(sp["largest_distance_m"], 2), "across_m": _r(sp["across_m"], 2)}
    out["final"] = _georef_final(paths)
    return out


def _case(data: bytes, name: str, g, h, mapping: dict | None = None) -> dict:
    """One tap map of the review read and suggested on the model of g (its flows h): the rows, what is not used, the
    notes and warnings, or the refusal."""
    from .importers import LogRefused, read_tap_coordinates, suggest_from
    try:
        s = suggest_from(read_tap_coordinates(data, name, mapping), g, h)
    except LogRefused as e:
        return {"refused": True, "message": str(e)}
    return {"refused": False, "kind": s.kind, "rows": s.table(), "not_used": s.dropped, "notes": s.notes,
            "warnings": s.warnings}


def part_review_16b(paths: dict) -> dict:
    """The review of task 16b, tap-map and model cases on the synthetic urban grid and its rewrites (fixtures.py,
    REVIEW_16B_FIXTURES): which pair of columns places the taps, thousands commas, hemisphere letters, swapped columns,
    a junction column known only through a mapping, a mapping that names a missing column, duplicate taps, a position
    far beyond the drawing, a tap beside a main, and the models drawn in degrees, in miles, with two nodes at 0,0 and
    drawn in part."""
    import warnings as _w
    from .importers import LogRefused, read_tap_map, suggest_tap_map, tap_geometry, tap_map_kind_of
    g, h, rules = _geom_hyd(paths["synthetic_urban_grid.inp"])
    b = lambda n: _read_bytes(paths[n])
    out = {}
    for n in ("taps_both_pairs_latlon_blank.csv", "taps_both_pairs_no_controls.csv", "taps_thousands_ambiguous.csv",
              "taps_thousands_settled.csv", "taps_decimal_comma.csv", "taps_letters_swapped.csv",
              "taps_latlon_swapped.csv", "taps_duplicate_positions.csv", "taps_duplicate_junctions.csv",
              "taps_beside_main.csv"):
        out[n] = _case(b(n), n, g, h)
    n = "taps_both_pairs_controls.csv"
    out[n] = {"guessed": _case(b(n), n, g, h), "x_and_y_chosen": _case(b(n), n, g, h, {"x": "x", "y": "y"})}
    out[n]["same_junctions"] = ([r["junction_id"] for r in out[n]["guessed"]["rows"]] ==
                                [r["junction_id"] for r in out[n]["x_and_y_chosen"]["rows"]])
    n = "taps_mapping_missing.csv"
    out[n] = _case(b(n), n, g, h, {"junction": "node ref"})
    n = "taps_junction_header_unknown.csv"
    out[n] = {"kind_guessed": tap_map_kind_of(b(n), n), "kind_with_mapping": tap_map_kind_of(b(n), n, {"junction": "EPANET_Name"}),
              "read_by_id_with_mapping": read_tap_map(b(n), n, {"junction": "EPANET_Name"})[0],
              "by_coordinates": _case(b(n), n, g, h)}
    n = "taps_absurd.csv"
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        out[n] = _case(b(n), n, g, h)
    out[n]["runtime_warnings"] = sum(issubclass(w.category, RuntimeWarning) for w in caught)
    models = {}
    for m in ("urban_grid_degrees.inp", "urban_grid_miles.inp", "urban_grid_placeholders.inp", "urban_grid_part_drawn.inp"):
        try:
            gm = tap_geometry(paths[m])
            models[m] = {"refused": False, "units": gm.scale["units"], "metres_per_unit": _r(gm.metres_per_unit, 3),
                         "pipes": len(gm.pipes), "left_out": gm.left_out, "warnings": gm.warnings}
        except LogRefused as e:
            models[m] = {"refused": True, "message": str(e)}
    try:
        suggest_tap_map(b("taps_both_pairs_latlon_blank.csv"), "t.csv", paths["urban_grid_degrees.inp"])
        models["urban_grid_degrees.inp"]["suggest_tap_map"] = {"refused": False}
    except LogRefused as e:
        models["urban_grid_degrees.inp"]["suggest_tap_map"] = {"refused": True, "message": str(e)}
    models["urban_grid_miles.inp"]["metres_per_unit_true"] = 1609.344
    models["urban_grid_miles.inp"]["metres_per_unit_as_first_built"] = _r(        # from the median rounded to 4 places
        1.0 / tap_geometry(paths["urban_grid_miles.inp"]).scale["ratio_median"], 3)
    # the urban part's 2,000 taps at 5 m on the miles drawing and on the drawing with two nodes at 0,0
    from .importers import suggest_from
    eligible = (g.diam_in < 16.0) & np.array([j is not None for j, _, _ in rules])
    P, tp = taps_on_pipes(g, eligible, SNAP_TAPS, SNAP_POS_SEED)
    Q = _noisy(P, SNAP_NOISE_M[0], SNAP_NOISE_SEEDS[0])
    base = suggest_from(_coords_of(g, Q), g, h)
    for m in ("urban_grid_miles.inp", "urban_grid_placeholders.inp"):
        gm, hm, rm = _geom_hyd(paths[m])
        xy = Q + np.array(g.origin)                                   # the true positions in metres, as model x, y
        if m == "urban_grid_miles.inp":
            xy = xy / 1609.344
        P2 = (xy - np.array(gm.origin)) * gm.metres_per_unit
        s = suggest_from(_coords_of(gm, P2), gm, hm)
        same = np.array([x["junction_id"] == y["junction_id"] for x, y in zip(s.rows, base.rows)])
        on_left_out = sum((x["pipe"] or "") in set(gm.left_out) for x in s.rows)
        models[m]["taps_5m"] = {"taps": len(s.rows), "same_junction_as_the_metres_drawing": _r(same.mean()),
                                "put_on_a_pipe_left_out": int(on_left_out)}
    out["models"] = models
    out["final"] = _review_final_16b(paths, g, h)
    return out


FINAL_TAPS = b"tap_id,x,y\nT1,503,403\nT2,1197,1110\nT3,700,1004\nT4,1503,1797\n"   # the sentinel models' taps


def _review_final_16b(paths: dict, g, h) -> dict:
    """The final verification of task 16b, tap-map and model cases (fixtures.FINAL_16B_FIXTURES): a dot thousands
    separator in a file that writes decimal commas, unsettled and settled; longitudes written from 0 to 360 (the
    same junctions as the template); control taps along one street of the urban grid with taps far across it; one
    node at a no-data value (-3.4e38, and nan) through suggest_tap_map, against the clean drawing; and the drawings in
    degrees with latitude as x and with one node at 9999,9999, refused."""
    from . import fixtures as FX
    from .importers import FAR_NODE_SPANS, LogRefused, suggest_tap_map, tap_geometry
    b = lambda n: _read_bytes(paths[n])
    out = {"far_node_spans": FAR_NODE_SPANS}
    for n in ("taps_dot_thousands.csv", "taps_dot_thousands_settled.csv", "taps_controls_along_one_street.csv"):
        out[n] = _case(b(n), n, g, h)
    n = "taps_lon_0_to_360.csv"
    out[n] = _case(b(n), n, g, h)
    tpl = _case(b("tap_map_by_coordinates.csv"), "tap_map_by_coordinates.csv", g, h)
    out[n]["same_junctions_as_the_template"] = ([r["junction_id"] for r in out[n]["rows"]] ==
                                                [r["junction_id"] for r in tpl["rows"]])
    clean = _case(FINAL_TAPS, "t.csv", g, h)
    models = {"clean drawing": {"rows": [[r["tap_id"], r["junction_id"], r["snap_m"]] for r in clean["rows"]]}}
    for m, v in FX.URBAN_SENTINELS.items():
        gm = tap_geometry(paths[m])
        s = suggest_tap_map(FINAL_TAPS, "t.csv", paths[m])
        rows = s.table()
        models[m] = {"node_value": v, "origin": [_r(x, 3) for x in gm.origin], "extent_m": _r(gm.extent_m, 1),
                     "pipes": len(gm.pipes), "warnings": gm.warnings,
                     "rows": [[r["tap_id"], r["junction_id"], r["snap_m"]] for r in rows],
                     "flags": [r["flags"] for r in rows], "not_used": s.dropped,
                     "same_rows_as_the_clean_drawing": [[r["tap_id"], r["junction_id"], r["snap_m"]] for r in rows] ==
                                                       models["clean drawing"]["rows"]}
    for m in ("urban_grid_degrees_lat_as_x.inp", "urban_grid_degrees_placeholder.inp"):
        try:
            tap_geometry(paths[m])
            models[m] = {"refused": False}
        except LogRefused as e:
            models[m] = {"refused": True, "message": str(e)}
    out["models"] = models
    return out


def part_ky4_after_the_stop() -> dict:
    """Measured after the stop, with nothing changed in the scored rule (not a rescoring): where the 66 misses of the
    scored ky4 run at 5 m are, and the share of taps given a wrong junction with no flag at all, at 5, 10 and 20 m,
    with the flags of the code after the review (a main close by is flagged since; ky4's five 16 in pipes)."""
    from .importers import suggest_from
    g, h, rules = _geom_hyd("ky4")
    eligible = (g.diam_in < 16.0) & np.array([j is not None for j, _, _ in rules])
    P, tp = taps_on_pipes(g, eligible, SNAP_TAPS, SNAP_POS_SEED)
    o = np.array(g.origin)
    end_xy = lambda n: (np.array(g.node_xy[n]) - o) * g.metres_per_unit
    d_end = np.array([min(np.hypot(*(P[i] - end_xy(g.start[p]))), np.hypot(*(P[i] - end_xy(g.end[p]))))
                      for i, p in enumerate(tp)])
    tj = [rules[p][0] for p in tp]
    out = {"by_noise": {}}
    for r, seed in zip(SNAP_NOISE_M, SNAP_NOISE_SEEDS):
        sug = suggest_from(_coords_of(g, _noisy(P, r, seed)), g, h)
        miss = np.array([x["junction_id"] != t for x, t in zip(sug.rows, tj)])
        z = {"taps": len(miss), "wrong_junction": int(miss.sum()),
             "wrong_junction_unflagged": int((miss & np.array([not x["flags"] for x in sug.rows])).sum())}
        z["wrong_junction_unflagged_share"] = _r(z["wrong_junction_unflagged"] / len(miss))
        if r == SNAP_NOISE_M[0]:
            wrong_pipe = np.array([x["pipe"] != g.pipes[p] for x, p in zip(sug.rows, tp)])
            two = np.array([any(f.startswith("two pipes close") for f in x["flags"]) for x in sug.rows])
            second_right = np.array([x["second_choice"] == t for x, t in zip(sug.rows, tj)])
            z.update({"misses_on_the_wrong_pipe": int((miss & wrong_pipe).sum()),
                      "misses_within_5m_of_an_end_node_of_the_true_pipe": int((miss & (d_end <= 5.0)).sum()),
                      "taps_within_5m_of_an_end_node": int((d_end <= 5.0).sum()),
                      "misses_flagged_two_pipes_close": int((miss & two).sum()),
                      "misses_with_the_right_junction_as_flagged_second_choice": int((miss & two & second_right).sum())})
        out["by_noise"][f"{r:g}"] = z
    return out


def part_template(paths: dict) -> dict:
    """The coordinates template, a fixture only (fixtures.coordinates_template; not shipped in docs/templates, since
    no app or pilot step reads a tap map by coordinates yet), on the synthetic urban grid: each tap's suggested
    junction, the fit's control taps and residuals, and the suggestion, written as a tap map by ID, read back through
    load_log with urban_coord_log.csv.  The pilot's own step (pilot_suggested_tap_map) is not built."""
    from .importers import suggest_tap_map
    from .pilot import load_log
    data = _read_bytes(paths["tap_map_by_coordinates.csv"])
    sug = suggest_tap_map(data, "tap_map_by_coordinates.csv", paths["synthetic_urban_grid.inp"])
    with tempfile.TemporaryDirectory(prefix="rm_t16b_") as d:
        p = os.path.join(d, "tap_map_suggested.csv")
        with open(p, "wb") as fh:
            fh.write(sug.to_csv())
        log = load_log(paths["urban_coord_log.csv"], p, quiet=True)
    return {"rows": [{k: v for k, v in r.items()} for r in sug.table()],
            "controls": sug.georef.get("controls"), "checked": sug.georef.get("checked"),
            "loo_max_m": _r(sug.georef.get("loo_max_m")), "notes": sug.notes, "warnings": sug.warnings,
            "not_used": sug.dropped, "accepted_rows_through_load_log": len(log),
            "accepted_junctions": sorted(set(log.junction))}


def snapping_record() -> dict:
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t16b_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.TAP_FIXTURES)
        return {"about": "Task 16b's tap mapping by coordinates, stopped at its pre-registered bar, on files built by "
                         "residualmap/fixtures.py: ky4 (a public network, feet, with vertices) and a synthetic urban "
                         "grid (invented, metres) for snapping, an invented Lambert conformal conic zone for the "
                         "georeference. Simulated GPS noise only, never a utility's own taps (R6). Part 'ky4' is the "
                         "pre-registered scored run, which missed its bar (ky4.bar.passed is false); part "
                         "'ky4_after_the_stop' measures its misses afterwards and is not a rescoring; part 'review' and "
                         "the 'review' and 'transmission_main_variants' entries hold the review's cases, and the "
                         "'final' entries (of part 'review', and of part 'georef' in its 'review' entry) the final "
                         "verification's, its plastic mains among the 'transmission_main_variants'. Written by "
                         "python -m residualmap.intake_exp snapping; six task16b checks in residualmap/checks.py (five of the "
                         "plan's six, the pilot's step and its check not being built, and tap_coordinates_review_cases) "
                         "rerun every part and compare.",
                "ky4": part_snap_ky4(), "ky4_after_the_stop": part_ky4_after_the_stop(),
                "urban_grid": part_snap_urban(paths), "direction": part_snap_direction(paths),
                "schematic": part_snap_schematic(paths), "georef": part_georef(paths), "template": part_template(paths),
                "review": part_review_16b(paths)}


# ----------------------------------------------------------------------------- task 16c: zone extracts
RECORD_16C = os.path.join(OUT, "zone_extract_Net3.json")
ZONE_NET = "Net3"                      # two entry points, its reservoirs River and Lake
ZONE_RATIOS = {"River": 1.0, "Lake": 0.5}
SUPERPOSE_BAR = {"default_tolerance": 5e-5, "tolerance_1e-9": 1e-6}   # mg/L, pre-registered (plan task 16c)
TIGHT_TOLERANCE = 1e-9                 # EPANET's quality tolerance for the second bar, in the model's units
DOSE_AXIS_MULT = 1.1                   # the dose axis's largest multiplier (simgp.DOSE_GRID)
ZONE_DOSE = 1.2
HOLDOUT_K = 3                          # holdout_sites_fixed_route: leave 3 sites out
FIXED_ROUTE = {"per_month": 10, "rotating": 0, "seed": 0, "months": 6}   # pilot.synthetic_log, no rotating taps
ZONE_CHANGE_MGL = 0.1                  # the effect part: a junction's daily mean 'moves' by more than this
GRID24_FULL_NET3_SHA256 = "bc849e767dc6e82e41ed3a3bc09e4b673fba8cb8b6698b15eb26efba006c99fb"   # checks.py's record


def _grid_sha256(params, Z) -> str:
    """checks._grid_sha256, the same bytes in the same order (intake_exp never imports the checks)."""
    import hashlib
    h = hashlib.sha256()
    h.update(json.dumps([[float(v) for v in p] for p in params]).encode())
    h.update(f"{Z.dtype.str}{Z.shape}".encode())
    h.update(np.ascontiguousarray(Z).tobytes())
    return h.hexdigest()


ZONE_GRID_MEMBER = (0.40, 0.60, 1.0)  # added after the scored run: MEMBER's kw 0.70 m/day is the truth's default and
                                       # not on the grid's kw axis (0.60 and 1.00 are its neighbours); this one is a
                                       # member of every free-chlorine grid


def _zone_run(ratios, tol=None, dose: float = ZONE_DOSE, keep_inp: bool = False, member=MEMBER):
    """One member of Net3 (kb, kw, gamma; demand and roughness x1) at these inlet ratios and quality tolerance, hour x
    junction (last day), run in a temporary folder; with keep_inp also the .inp WNTR wrote, without its date line."""
    from .simulate import simulate_nominal_chlorine
    with tempfile.TemporaryDirectory(prefix="rm_t16c_") as d:
        p = os.path.join(d, "m")
        c = simulate_nominal_chlorine(ZONE_NET, *member, dose, 1.0, 1.0, file_prefix=p, inlet_ratios=ratios,
                                      quality_tolerance=tol)
        if not keep_inp:
            return c
        with open(p + ".inp") as fh:
            return c, "".join(line for line in fh if not line.startswith("; Created:"))


def part_zone_identity(fresh_grid: bool = True) -> dict:
    """inlet_ratios_identity: ratios (1, 1) are today's path.  They normalize to None, so do two equal residuals; the
    cache name is today's; one member's .inp is byte for byte today's (with ratios given and without); and, with
    fresh_grid, a fresh 675-run grid at ratios (1, 1), built in memory, has the sha256 recorded for the cached grid."""
    from .importers import inlet_ratios_from
    from .simgp import build_grid_24h, grid_cache_path
    from .simulate import nominal_scenario, normalize_inlet_ratios, source_nodes
    sc = nominal_scenario(ZONE_NET)
    entry = source_nodes(sc.wn)
    ones = {e: 1.0 for e in entry}
    c0, t0 = _zone_run(None, keep_inp=True)
    c1, t1 = _zone_run(ones, keep_inp=True)
    out = {"entry_points": list(entry), "ratios": ones, "normalized": normalize_inlet_ratios(ones),
           "equal_residuals_give_ratios": inlet_ratios_from({e: 0.8 for e in entry}, entry)[0],
           "cache_name_with_ratios": os.path.basename(grid_cache_path(sc, "outputs/cache", "full", None, ones)),
           "cache_name_today": os.path.basename(grid_cache_path(sc, "outputs/cache", "full")),
           "member_inp_identical": t0 == t1, "member_values_identical": bool(c0.equals(c1))}
    if fresh_grid:
        p, Z = build_grid_24h(sc, "full", preflight_paths=[os.path.join(REPO, "outputs", "cache")], inlet_ratios=ones)
        out.update(fresh_grid_members=len(p), fresh_grid_sha256=_grid_sha256(p, Z),
                   fresh_grid_sha256_matches_recorded=_grid_sha256(p, Z) == GRID24_FULL_NET3_SHA256)
    return out


def part_zone_superposition() -> dict:
    """inlet_ratios_superpose: Net3 at ratios (1.0, 0.5), one member (MEMBER: kb 0.40, kw 0.70 m/day, the truth's
    default, which is off the grid's kw axis), against C(1, 0) + 0.5 C(0, 1), at the file's quality tolerance and at
    1e-9; the pre-registered bars; the dose axis at fixed ratios (information, no bar); how much the half residual at
    Lake moves the map at this member; and, added after the scored run, the same test at a member on the grid
    (ZONE_GRID_MEMBER)."""
    from .simulate import load
    tol_file = float(load(ZONE_NET).options.quality.tolerance)
    entry = list(ZONE_RATIOS)
    only = {e: {f: (1.0 if f == e else 0.0) for f in entry} for e in entry}
    out = {"network": ZONE_NET, "member": {"kb_per_day": MEMBER[0], "kw_m_per_day": MEMBER[1], "gamma": MEMBER[2],
                                           "demand": 1.0, "rough": 1.0, "dose_mgL": ZONE_DOSE},
           "ratios": dict(ZONE_RATIOS), "file_tolerance": tol_file,
           "file_tolerance_in_grid_units_mgL": tol_file / 1000.0,
           "units_note": "the grid writes each dose to EPANET in the model's units, a thousand times its mg/L, so the "
                         "file's tolerance of 0.01 is 1e-5 mg/L in the grid's units",
           "bar_mgL": dict(SUPERPOSE_BAR), "by_tolerance": {}}
    runs = {}
    for label, tol in (("default_tolerance", None), ("tolerance_1e-9", TIGHT_TOLERANCE)):
        both = _zone_run(ZONE_RATIOS, tol).values
        parts = {e: _zone_run(only[e], tol).values for e in entry}
        summed = sum(ZONE_RATIOS[e] * parts[e] for e in entry)
        err = float(np.abs(both - summed).max())
        dose_up = _zone_run(ZONE_RATIOS, tol, dose=ZONE_DOSE * DOSE_AXIS_MULT).values
        d_err = float(np.abs(dose_up - DOSE_AXIS_MULT * both).max())
        out["by_tolerance"][label] = {"tolerance": tol_file if tol is None else tol, "max_abs_error_mgL": err,
                                      "junction_hours": int(both.size), "max_concentration_mgL": float(both.max()),
                                      "bar_mgL": SUPERPOSE_BAR[label], "passed": err <= SUPERPOSE_BAR[label],
                                      "dose_axis_max_abs_diff_mgL": d_err, "dose_axis_multiplier": DOSE_AXIS_MULT}
        runs[label] = both
    out["passed"] = all(v["passed"] for v in out["by_tolerance"].values())
    # added after the scored run (the bars applied, the scored result above unchanged): the same test at a member that
    # is on the grid's axes
    on = {"member": {"kb_per_day": ZONE_GRID_MEMBER[0], "kw_m_per_day": ZONE_GRID_MEMBER[1],
                     "gamma": ZONE_GRID_MEMBER[2], "demand": 1.0, "rough": 1.0, "dose_mgL": ZONE_DOSE},
          "why": "the scored member's kw (0.70 m/day, the truth's default) is not on the grid's kw axis; this member is "
                 "on every free-chlorine grid; added after the scored run, which stands as scored",
          "by_tolerance": {}}
    for label, tol in (("default_tolerance", None), ("tolerance_1e-9", TIGHT_TOLERANCE)):
        both = _zone_run(ZONE_RATIOS, tol, member=ZONE_GRID_MEMBER).values
        summed = sum(ZONE_RATIOS[e] * _zone_run(only[e], tol, member=ZONE_GRID_MEMBER).values for e in entry)
        err = float(np.abs(both - summed).max())
        on["by_tolerance"][label] = {"max_abs_error_mgL": err, "bar_mgL": SUPERPOSE_BAR[label],
                                     "within_bar": err <= SUPERPOSE_BAR[label]}
    out["on_grid_member"] = on
    today = _zone_run(None).values
    half = runs["default_tolerance"]
    dm = today.mean(axis=0) - half.mean(axis=0)
    out["effect"] = {"what": "this member with Lake's residual half River's, against every inlet at the one dose "
                             "(today's path), daily means at each junction",
                     "junctions": int(today.shape[1]), "moved_more_than_mgL": ZONE_CHANGE_MGL,
                     "junctions_moved": int((np.abs(dm) > ZONE_CHANGE_MGL).sum()),
                     "largest_drop_mgL": float(dm.max()),
                     "junctions_daily_min_below_0.2_today": int((today.min(axis=0) < 0.2).sum()),
                     "junctions_daily_min_below_0.2_half_at_lake": int((half.min(axis=0) < 0.2).sum())}
    return out


def _subset_rows(summary, subset: str) -> list:
    keep = ("held_out_month", "n_train", "n_test", "rmse", "mae", "rmse_network_mean", "coverage50", "coverage80",
            "coverage90", "coverage95", "band90_width_mgL", "n_true_below", "recall_below", "false_alarms")
    s = summary[summary.taps == subset]
    extra = (("sites", "site_groups", "sites_per_group", "n_train_min", "n_train_max") if subset == "held_out_sites"
             else ())
    return [{k: (int(r[k]) if k in ("n_train", "n_test", "n_true_below", "false_alarms") + extra else r[k])
             for k in keep + extra} for _, r in s.iterrows()]


def _subset_mean(summary, subset: str) -> dict:
    from .pilot import LEVELS
    g = summary[summary.taps == subset]
    w = g.n_test / g.n_test.sum()
    viol = g.n_true_below.sum()
    return {"n": int(g.n_test.sum()), "rmse": float((g.rmse * w).sum()),
            "rmse_network_mean": float((g.rmse_network_mean * w).sum()),
            **{f"coverage{q}": float((g[f"coverage{q}"] * w).sum()) for q in LEVELS},
            "n_true_below": int(viol), "found_below": int(round(float((g.recall_below.fillna(0) * g.n_true_below).sum()))),
            "false_alarms": int(g.false_alarms.sum())}


def part_zone_holdout() -> dict:
    """holdout_sites_fixed_route: pilot.validate with holdout_sites=3 on a synthetic Net3 log of a fixed route (10
    taps, no rotating taps, six months; the hold-out is the last three).  Every site of each held-out month is
    predicted once by a fit that never saw it; its table and means, beside the same months' 'all' rows (the fixed
    route's own taps, seen in training: a month-ahead forecast at sampled places)."""
    from .pilot import HOLDOUT_SITES_SEED, synthetic_log, validate
    cache = os.path.join(REPO, "outputs", "cache")
    log, inp = synthetic_log(ZONE_NET, **FIXED_ROUTE)
    preds, summary = validate(inp, log, dose=ZONE_DOSE, threshold=0.2, window_months=3, cache_dir=cache,
                              holdout_sites=HOLDOUT_K)
    ps = preds[preds.validation == "held_out_sites"]
    test = preds[preds.validation == "held_out_month"]
    per = lambda f: f.groupby(["month", "junction"]).size().sort_index()
    groups = ps.groupby(["month", "site_group"]).junction.nunique()
    return {"network": ZONE_NET, "log": dict(FIXED_ROUTE), "taps": sorted(set(log.tap_id)), "samples": int(len(log)),
            "months": sorted(set(log.month)), "held_out_months": sorted(set(summary.held_out_month)),
            "k": HOLDOUT_K, "seed": HOLDOUT_SITES_SEED + FIXED_ROUTE["seed"],
            "every_site_predicted_once": bool(per(ps).equals(per(test))), "held_out_predictions": int(len(ps)),
            "sites_per_group": sorted(int(v) for v in groups.values),
            "by_month": {"held_out_sites": _subset_rows(summary, "held_out_sites"),
                         "all": _subset_rows(summary, "all")},
            "mean": {"held_out_sites": _subset_mean(summary, "held_out_sites"), "all": _subset_mean(summary, "all")},
            "subsets": sorted(set(summary.taps))}


def part_zone_analyzer(paths: dict) -> dict:
    """analyzer_is_boundary_not_data: analyzer_15min.csv read with zone_tap_map.csv, whose role column makes AN-1
    Lake's inlet analyzer.  Its 2,880 readings leave the frame the model is fitted on (S) and set Lake's residual,
    their median, equal to the median of what the fixture wrote; River's residual comes from an inlet residuals file;
    r_Lake = Lake's residual over River's.  The same log with the task-16 tap map (no role column) keeps AN-1 as an
    analyzer, thinned to one reading a day, as before.  And the refusals: an inlet missing, a node that is no entry
    point, a residual below detection, an inlet given twice (the file and an inlet tap)."""
    from . import fixtures as FX
    from .importers import (LogRefused, combine_inlet_residuals, inlet_ratios_from, read_grab_log,
                            read_inlet_residuals)
    rd = lambda n: _read_bytes(paths[n])
    tap, inlet = FX.ZONE_INLET_TAP
    g = read_grab_log(rd("analyzer_15min.csv"), "analyzer_15min.csv", rd("zone_tap_map.csv"), "zone_tap_map.csv")
    g0 = read_grab_log(rd("analyzer_15min.csv"), "analyzer_15min.csv", rd("task16_tap_map.csv"), "task16_tap_map.csv")
    S = g.frame[["junction", "hour", "y"]]
    taps_res = g.inlet_residuals()
    window = g.window()[2]["months"]
    file_res, _ = read_inlet_residuals(rd("zone_inlet_residuals.csv"), "zone_inlet_residuals.csv")
    comb = combine_inlet_residuals(file_res, taps_res)
    ratios, info = inlet_ratios_from(comb, list(ZONE_RATIOS))
    expected = float(np.median(FX.analyzer_values(FX.ANALYZER_SITES.index(tap))))

    def refusal(fn) -> str | None:
        try:
            fn()
        except LogRefused as e:
            return str(e)
        return None
    both = {inlet: 0.5, "River": 1.2}
    review = _zone_review_cases(paths)
    return {"inlet_tap": tap, "inlet": inlet, "inlet_readings": int(len(g.inlets)),
            "inlet_rows_in_S": int((g.frame.tap_id == tap).sum()), "rows_in_S": int(len(S)),
            "S_junctions": sorted(set(S.junction)), "other_analyzer_rows_in_S": int((g.frame.tap_id == "AN-2").sum()),
            "grab_rows_in_S": int((~g.frame.tap_id.isin(FX.ANALYZER_SITES)).sum()),
            "inlet_residual_mgL": taps_res[inlet]["residual_mgL"], "expected_median_mgL": expected,
            "inlet_residual_in_the_window_mgL": g.inlet_residuals(window)[inlet]["residual_mgL"], "window": window,
            "river_from_file_mgL": file_res["River"], "ratios": ratios, "highest_mgL": info["highest_mgL"],
            "r_inlet_expected": round(expected / FX.ZONE_RIVER_MGL, 6),
            "note": next(n for n in g.notes if "zone inlet" in n),
            "without_role_column": {"inlet_tap_rows_in_S": int((g0.frame.tap_id == tap).sum()),
                                    "inlets": g0.inlets, "rows_in_S": int(len(g0.frame))},
            "ugl_file": read_inlet_residuals(rd("zone_inlets_ugl.csv"), "zone_inlets_ugl.csv")[0],
            "refused": {
                "inlet_missing": refusal(lambda: inlet_ratios_from(read_inlet_residuals(rd("zone_inlets_missing.csv"))[0],
                                                                   list(ZONE_RATIOS))),
                "not_an_entry_point": refusal(lambda: inlet_ratios_from(
                    read_inlet_residuals(rd("zone_inlets_unknown.csv"))[0], list(ZONE_RATIOS))),
                "below_detection": refusal(lambda: read_inlet_residuals(rd("zone_inlets_censored.csv"))),
                "given_twice": refusal(lambda: combine_inlet_residuals(both, taps_res))},
            "review": review}


def _zone_review_cases(paths: dict) -> dict:
    """The review of task 16c's cases, read as the app and the pilot read them: a declared inlet tap (INL, Lake's)
    whose 3 readings are all before the grab rows' 3-month window is refused in the app's window and used over the
    pilot's whole log, its readings listed as not used in the window; the same tap read only below detection is
    refused everywhere; role cells that mention an inlet but mark none are read as taps, with a warning; a tap marked
    as an inlet at a junction is refused naming the tap and its role; and three inlet files are refused (a µg/L value
    with a thousands separator, a total chlorine column in free-chlorine mode, a residual under a hundredth of the
    highest), while the µg/L file with decimal commas is read with its two notes."""
    from . import fixtures as FX
    from .importers import (LOG_WINDOW_MONTHS, LogRefused, inlet_ratios_from, inlet_tap_residuals, read_grab_log,
                            read_inlet_residuals)
    from .simulate import load, source_nodes
    rd = lambda n: _read_bytes(paths[n])
    wn = load(ZONE_NET)
    entry, J = source_nodes(wn), wn.junction_name_list

    def refusal(fn) -> str | None:
        try:
            fn()
        except LogRefused as e:
            return str(e)
        return None
    log = lambda name, taps: read_grab_log(rd(name), name, rd(taps), taps).only_junctions(J)
    g = log("zone_inlet_window_log.csv", "zone_window_tap_map.csv")
    months = g.window(LOG_WINDOW_MONTHS)[2]["months"]
    whole = inlet_tap_residuals(g, entry, None)
    gc = log("zone_inlet_censored_log.csv", "zone_window_tap_map.csv")
    gn = log("zone_inlet_window_log.csv", "zone_tap_map_not_inlet.csv")
    gw = log("zone_inlet_window_log.csv", "zone_tap_map_wrong_node.csv")
    return {"outside_window": {"inlet_readings": int(len(g.inlets)), "window": months,
                               "inlet_months": sorted(set(g.inlets.month)),
                               "app_window_refused": refusal(lambda: inlet_tap_residuals(g, entry, months)),
                               "pilot_whole_log_mgL": whole[FX.ZONE_WINDOW_INLET[1]]["residual_mgL"],
                               "expected_median_mgL": float(np.median(FX.ZONE_WINDOW_INLET_MGL)),
                               "listed_not_used_in_window": len(g.inlet_unused(months)),
                               "note": next(n for n in g.notes if "zone inlet" in n)},
            "only_below_detection": {"usable": int(len(gc.inlets)),
                                     "listed_not_used": sum(d["reason"].startswith("zone inlet ") for d in gc.dropped),
                                     "refused": refusal(lambda: inlet_tap_residuals(gc, entry, None)),
                                     "note": next(n for n in gc.notes if "zone inlet" in n)},
            "roles_that_mark_no_inlet": {"roles": list(FX.ZONE_NEAR_MISS_ROLES), "inlet_taps": len(gn.inlet_taps),
                                         "rows_fitted": int(len(gn.frame)), "rows_in_log": gn.n_rows,
                                         "warning": next(w for w in gn.warnings if "mention" in w)},
            "inlet_at_a_junction_refused": refusal(lambda: inlet_tap_residuals(gw, entry, None)),
            "files_refused": {
                "thousands_separator_ugl": refusal(lambda: read_inlet_residuals(rd("zone_inlets_thousands_ugl.csv"))),
                "total_chlorine_column": refusal(lambda: read_inlet_residuals(rd("zone_inlets_total.csv"))),
                "under_a_hundredth": refusal(lambda: inlet_ratios_from(
                    read_inlet_residuals(rd("zone_inlets_tiny.csv"))[0], entry))},
            "ugl_file_notes": read_inlet_residuals(rd("zone_inlets_ugl.csv"), "zone_inlets_ugl.csv")[1]["notes"]}


def zone_record(fresh_grid: bool = True) -> dict:
    from . import fixtures as FX
    _quiet()
    with tempfile.TemporaryDirectory(prefix="rm_t16c_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.ZONE_FIXTURES)
        return {"about": "Task 16c's zone extracts on Net3, a public network with two entry points (its reservoirs "
                         "River and Lake), simulated only (R6): each inlet dosed at its own share of the dose, the "
                         "superposition of the two inlets at one member (kb 0.40 per day, kw 0.70 m/day, the truth's "
                         "default and off the grid's kw axis; a member on the grid added after the scored run) against "
                         "the pre-registered bars (part 'superposition'), ratios (1, 1) as today's path (part "
                         "'identity'), leave-sites-out validation on a synthetic fixed route of 10 taps, with each group "
                         "fit's training size (part 'holdout'), and an online analyzer at a zone inlet read as a "
                         "boundary condition, not data, with the review's cases of declared inlets (part 'analyzer', "
                         "invented readings). Written by python -m residualmap.intake_exp zone; the task16c checks in "
                         "residualmap/checks.py rerun every part and compare.",
                "identity": part_zone_identity(fresh_grid), "superposition": part_zone_superposition(),
                "holdout": part_zone_holdout(), "analyzer": part_zone_analyzer(paths)}


# ----------------------------------------------------------------------------- task 17: gis_roundtrip.json
RECORD_17 = os.path.join(OUT, "gis_roundtrip.json")
GIS_BAR = {"tolerance_ft": 6.0, "pipes_at_least": 1151, "of_pipes": 1156, "components": 1,
           "diameter_share_at_least": 0.99}                     # pre-registered (plan task 17), every scored fixture
GIS_DIAG_FT = (1.0, 3.0, 6.0)          # the scored fixtures' tolerances (the bar at 6 ft; 1 and 3 ft recorded)
DIAMETER_REL_TOL = 0.01                # a rebuilt pipe's diameter matches within 1% (a mm size rounds by under 0.7%)
WNTR_MEMBER = (0.40, 0.70, 1.0)        # gis_wntr_schema_roundtrip's chlorine member: kb 1/day, kw m/day, gamma
NEGATIVE_BAR_MGL = 0.5                 # its negative control (no pattern workaround) must exceed this on Net3


def _gis(data: bytes, name: str, **kw):
    from . import gis_read as G
    return G.build_gis(data, name, G.GisOptions(**kw))


def score_ky4(b, index_of=lambda j: j, man: dict | None = None) -> dict:
    """A rebuilt ky4 mains layer scored by endpoint-pair identity, not counts: each ky4 node stands for the rebuilt
    node most of the line ends drawn at it fell into (so nodes the tolerance joins, the pump node pairs drawn at one
    point and J-702 and J-703, stand for one another), and a ky4 pipe is rebuilt when a rebuilt pipe joins its two
    nodes' rebuilt nodes (each rebuilt pipe used once; one whose two ends are one node is a miss).  Its diameter
    matches within DIAMETER_REL_TOL.  index_of: a record's feature index in the file.  man: another network's
    manifest in ky4_manifest's shape (task 18: ky10's drawing), scored by the same rule."""
    from . import fixtures as FX
    from collections import Counter
    man = man or FX.ky4_manifest()
    net = b.net
    votes: dict = {}
    for j, r in enumerate(man["records"]):
        f = index_of(j)
        for end in (0, 1):
            n = net.end_node.get((f, 0, end))
            if n is not None:
                votes.setdefault(r["ends"][end], Counter())[n] += 1
    rebuilt_of = {u: sorted(c.items(), key=lambda x: (-x[1], x[0]))[0][0] for u, c in votes.items()}
    pool: dict = {}
    for p in net.pipes:
        pool.setdefault(frozenset((p["start"], p["end"])), []).append(b.diam_in[p["id"]])
    matched, d_ok, misses, collapsed = 0, 0, [], []
    for name, (s, e, d) in sorted(man["pipes"].items()):
        a, c = rebuilt_of.get(s), rebuilt_of.get(e)
        if a is None or c is None or a == c:
            misses.append(name)
            collapsed.append(name) if a is not None and a == c else None
            continue
        k = frozenset((a, c))
        if not pool.get(k):
            misses.append(name)
            continue
        ds = pool[k]
        best = min(range(len(ds)), key=lambda i: abs(ds[i] - d))
        got = ds.pop(best)
        matched += 1
        d_ok += abs(got - d) <= DIAMETER_REL_TOL * d
    extra = sum(len(v) for v in pool.values())
    r = net.report
    return {"pipes_rebuilt": matched, "of_pipes": len(man["pipes"]), "missed": misses, "collapsed": collapsed,
            "diameters_matched": d_ok, "diameter_share": d_ok / max(matched, 1), "extra_pipes": extra,
            "components": r["components"], "nodes": r["nodes"], "dead_ends": r["dead_ends"],
            "tees_split": r["tees_split"], "nodes_flagged_wide": r["nodes_flagged_wide"], "nodes_split": r["nodes_split"],
            "dropped_pieces": r["dropped_pieces"], "pipes_built": r["pipes"]}


def _bar_met(s: dict) -> bool:
    return (s["pipes_rebuilt"] >= GIS_BAR["pipes_at_least"] and s["components"] == GIS_BAR["components"]
            and s["diameter_share"] >= GIS_BAR["diameter_share_at_least"])


def part_gis_rebuild(paths: dict) -> dict:
    """gis_utility_style_rebuild: the four ky4 mains fixtures at 1, 3 and 6 ft, the bar at 6 ft."""
    from . import fixtures as FX
    from .topology import FT
    out = {"bar": GIS_BAR, "fixtures": {}}
    for name in FX.GIS_SCORED:
        data = _read_bytes(paths[name])
        idx = FX.ky4_pyshp_index if name == "ky4_mains_pyshp.zip" else (lambda j: j)
        by = {}
        for ft in GIS_DIAG_FT:
            b = _gis(data, name, tolerance_m=ft * FT)
            assert b.net is not None, (name, b.refusal, b.questions)
            by[f"{ft:g}"] = score_ky4(b, idx)
        b = _gis(data, name)
        lay = b.layer
        out["fixtures"][name] = {
            "reader": {"ky4_mains_stdlib_shp.zip": "pyshp, on a file written by the fixtures' own shapefile writer",
                       "ky4_mains_pyshp.zip": "pyshp, on a file pyshp wrote (same-library self-consistency)",
                       "ky4_mains_z.gpkg": "sqlite3 and the WKB parser, on a GeoPackage written by the fixtures' WKB "
                                           "writer",
                       "ky4_mains_lonlat.geojson": "json and the east-north plane, on longitude and latitude from the "
                                                   "fixtures' Lambert conformal conic inverse"}[name],
            "default_tolerance_ft": b.tolerance_m / FT, "lines": lay.n_lines, "crs": {k: lay.crs.get(k) for k in
                                                                                  ("kind", "unit", "how")},
            "dropped": lay.dropped, "by_tolerance_ft": by, "bar_met": _bar_met(by["6"])}
    out["bar_met_on_every_fixture"] = all(v["bar_met"] for v in out["fixtures"].values())
    out["near_node_misses"] = _near_node_misses(out["fixtures"]["ky4_mains_stdlib_shp.zip"]["by_tolerance_ft"]["6"])
    return out


def _near_node_misses(s: dict) -> dict:
    """For each ky4 pipe the 6 ft rebuild missed whose two nodes did not collapse into one: the nearest ky4 node that
    is not on it, and its distance to the pipe in ky4.inp's own geometry, in feet (the review of task 17 asked that
    the journal's 'under a foot' be traced to the record)."""
    from . import fixtures as FX
    pipes, xy, _, _ = FX._ky4_parts()
    out = {}
    for pn in s["missed"]:
        if pn in s["collapsed"]:
            continue
        a, c, pts, _, _ = pipes[pn]
        P = np.asarray(pts, float)
        A, B = P[:-1], P[1:]
        best = (math.inf, None)
        for n, q in xy.items():
            if n in (a, c):
                continue
            q = np.asarray(q, float)
            AB = B - A
            t = np.clip(((q - A) * AB).sum(1) / np.maximum((AB * AB).sum(1), 1e-300), 0.0, 1.0)
            d = float(np.hypot(*(A + t[:, None] * AB - q).T).min())
            if d < best[0]:
                best = (d, n)
        out[pn] = {"node": best[1], "distance_ft": round(best[0], 3)}
    return out


def part_gis_kml(paths: dict) -> dict:
    """gis_kml_handdrawn: ky4 drawn by hand (each line end about 5 m off) at 5, 10, 15 and 20 m; no bar."""
    from . import fixtures as FX
    data = _read_bytes(paths["ky4_googleearth.kmz"])
    out = {"endpoint_error_m": FX.KML_ERROR_M, "default_tolerance_m": None, "by_tolerance_m": {}}
    man = FX.ky4_manifest()
    for tol in FX.KML_TOLERANCES_M:
        b = _gis(data, "ky4_googleearth.kmz", tolerance_m=tol)
        out["default_tolerance_m"] = out["default_tolerance_m"] or _gis(data, "ky4_googleearth.kmz").tolerance_m
        # each placemark is one ky4 pipe, in ky4's order: score with one record per pipe
        names = list(man["pipes"])
        s = _score_pipes(b, names)
        out["by_tolerance_m"][f"{tol:g}"] = s
    out["diagnostic_table_at_default"] = [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}
                                          for row in _gis(data, "ky4_googleearth.kmz").diagnostics]
    return out


def _score_pipes(b, names: list) -> dict:
    """score_ky4 for a layer with one feature per ky4 pipe, in names' order."""
    from . import fixtures as FX
    from collections import Counter
    pipes = FX.ky4_manifest()["pipes"]
    votes: dict = {}
    for f, n in enumerate(names):
        for end, node in ((0, pipes[n][0]), (1, pipes[n][1])):
            r = b.net.end_node.get((f, 0, end))
            if r is not None:
                votes.setdefault(node, Counter())[r] += 1
    of = {u: sorted(c.items(), key=lambda x: (-x[1], x[0]))[0][0] for u, c in votes.items()}
    pool = Counter(frozenset((p["start"], p["end"])) for p in b.net.pipes)
    ok = 0
    for n in names:
        a, c = of.get(pipes[n][0]), of.get(pipes[n][1])
        if a is not None and c is not None and a != c and pool[frozenset((a, c))] > 0:
            pool[frozenset((a, c))] -= 1
            ok += 1
    r = b.net.report
    return {"pipes_rebuilt": ok, "of_pipes": len(names), "components": r["components"], "dead_ends": r["dead_ends"],
            "nodes": r["nodes"], "tees_split": r["tees_split"], "nodes_flagged_wide": r["nodes_flagged_wide"],
            "nodes_split": r["nodes_split"], "dropped_pieces": r["dropped_pieces"]}


def part_gis_quirks(paths: dict) -> dict:
    """gis_spec_quirks: the shapefile and GeoPackage quirks of critique 12 and the CRS and unit quirks of critique 14."""
    from . import fixtures as FX
    out = {}
    for n in ("multipart_two_parts.gpkg", "multipart_two_parts.zip"):
        b = _gis(_read_bytes(paths[n]), n)
        r = b.net.report
        out[n] = {"pipes": r["pipes"], "nodes": r["nodes"], "components": r["components"],
                  "bridge": any(abs(p["xy"][0][0] - p["xy"][-1][0]) < 1e-9 for p in b.net.pipes) and False}
        # a bridge would join the first part's last point (1100, 1000) to the second part's first (1300, 1000)
        ends = {tuple(round(c / 0.3048006096012192 + o, 3) for c, o in zip(p["xy"][k], b.frame["origin"]))
                for p in b.net.pipes for k in (0, -1)}
        out[n]["bridge"] = any({(1100.0, 1000.0), (1300.0, 1000.0)} <= {
            tuple(round(c / b.frame["unit_m"] + o, 3) for c, o in zip(p["xy"][k], b.frame["origin"])) for k in (0, -1)}
            for p in b.net.pipes)
        out[n]["end_points"] = len(ends)
    for n in ("null_shape_mid.zip", "deleted_dbf_record.zip"):
        b = _gis(_read_bytes(paths[n]), n)
        # each rebuilt pipe's diameter against the line it was drawn from (lines in QUIRK_DIAMETERS' order by x)
        got = []
        for p in sorted(b.net.pipes, key=lambda p: min(v[0] for v in p["xy"])):
            got.append(round(b.diam_in[p["id"]], 3))
        out[n] = {"pipes": len(b.net.pipes), "diameters_in_line_order": got,
                  "aligned": got == [float(d) for d in FX.QUIRK_DIAMETERS], "dropped": b.layer.dropped,
                  "notes": b.layer.notes}
    b = _gis(_read_bytes(paths["cpg_utf8.zip"]), "cpg_utf8.zip")
    out["cpg_utf8.zip"] = {"streets": [f.props.get("STREET") for f in b.layer.features],
                           "decoded": [f.props.get("STREET") for f in b.layer.features] == list(FX.CPG_NAMES),
                           "encoding": b.layer.encoding}
    b = _gis(_read_bytes(paths["wkb_z_m_zm.gpkg"]), "wkb_z_m_zm.gpkg")
    out["wkb_z_m_zm.gpkg"] = {"features": len(b.layer.features),
                              "flags": [{"z": f.has_z, "m": f.has_m} for f in b.layer.features],
                              "read_in_2d": [f.parts[0].tolist() for f in b.layer.features] ==
                                            [[list(p) for p in ln] for ln in FX.wkb_z_m_zm_expected()],
                              "pipes": len(b.net.pipes)}
    b = _gis(_read_bytes(paths["prj_stateplane_ftus.zip"]), "prj_stateplane_ftus.zip")
    out["prj_stateplane_ftus.zip"] = {"unit": b.layer.crs["unit"], "unit_m": b.layer.crs["unit_m"],
                                      "how": b.layer.crs["how"]}
    b = _gis(_read_bytes(paths["gpkg_srs_undefined.gpkg"]), "gpkg_srs_undefined.gpkg")
    out["gpkg_srs_undefined.gpkg"] = {"questions": [q["key"] for q in b.questions], "built": b.inp_text is not None,
                                      "question": b.questions[0]["text"] if b.questions else None}
    b2 = _gis(_read_bytes(paths["gpkg_srs_undefined.gpkg"]), "gpkg_srs_undefined.gpkg", unit="ft")
    out["gpkg_srs_undefined.gpkg"]["with_unit_ft"] = {"built": b2.inp_text is not None, "pipes": len(b2.net.pipes)}
    b = _gis(_read_bytes(paths["shape_length_mismatch.zip"]), "shape_length_mismatch.zip")
    out["shape_length_mismatch.zip"] = {**{k: (round(v, 4) if isinstance(v, float) else v)
                                           for k, v in b.length_check.items()}}
    b = _gis(_read_bytes(paths["quantized_4dp.geojson"]), "quantized_4dp.geojson")
    out["quantized_4dp.geojson"] = {"decimals": b.layer.crs.get("quantized_decimals"),
                                    "warned": any("decimals of a degree" in w for w in b.warnings)}
    out.update(_review_quirks(paths))
    return out


def _review_quirks(paths: dict) -> dict:
    """The quirks the review of task 17 found: Web Mercator read as ground metres, Esri JSON's code-only
    spatialReference, GeoJSON positions dropped without a word, a preview KML in a shapefile zip, two shapefiles of one
    name, KML sizes guessed from names and attribute tables, and an empty install-year field."""
    from . import fixtures as FX
    from . import gis_read as G
    from .importers import Readiness, sniff
    out = {}
    for n in ("mercator_3857.geojson", "mercator_esri_prj.zip"):
        data = _read_bytes(paths[n])
        b = _gis(data, n)
        got = sorted(p["length_m"] for p in b.net.pipes)
        want = sorted(FX.MERC_LENGTHS_M)
        naive = sorted(math.hypot(ln[1][0] - ln[0][0], ln[1][1] - ln[0][1]) for ln in FX._mercator_lines())
        out[n] = {"crs_kind": b.layer.crs["kind"], "from_web_mercator": bool(b.layer.crs.get("from_web_mercator")),
                  "lengths_m": [round(x, 4) for x in got], "ground_lengths_m": want,
                  "max_rel_error": max(abs(g - w) / w for g, w in zip(got, want)),
                  "as_mercator_metres_m": [round(x, 2) for x in naive], "notes": b.notes}
    b = _gis(_read_bytes(paths["esri_code_only.json"]), "esri_code_only.json")
    out["esri_code_only.json"] = {"refusal": b.refusal, "epsg": b.layer.crs.get("epsg"), "unit": b.layer.crs.get("unit"),
                                  "pipes": len(b.net.pipes)}
    b = _gis(_read_bytes(paths["geojson_bad_positions.geojson"]), "geojson_bad_positions.geojson")
    out["geojson_bad_positions.geojson"] = {"lines_read": b.layer.n_lines, "features": len(b.layer.features),
                                            "warned": [w for w in b.warnings if "cannot be read" in w]}
    data = _read_bytes(paths["shp_with_preview_kml.zip"])
    b = _gis(data, "shp_with_preview_kml.zip")
    out["shp_with_preview_kml.zip"] = {"sniffed": sniff(data, "shp_with_preview_kml.zip").kind, "layers": b.layers,
                                       "pipes": len(b.net.pipes), "notes": [x for x in b.notes if ".kml" in x]}
    layers, _ = G.read_layers(_read_bytes(paths["two_folders_same_name.zip"]), "two_folders_same_name.zip")
    out["two_folders_same_name.zip"] = {"layers": [[x.name, x.n_lines] for x in layers]}
    for n in ("kml_numbered_paths.kml", "kml_attribute_table.kml"):
        b = _gis(_read_bytes(paths[n]), n)
        out[n] = {"diameter_field": b.fields.get("diameter"), "year_field": b.fields.get("year"),
                  "questions": [q["key"] for q in b.questions],
                  "sizes_in": sorted({round(v, 3) for v in b.diam_in.values()}),
                  "length_check_agrees": b.length_check.get("agrees")}
    b = _gis(_read_bytes(paths["year_all_null.geojson"]), "year_all_null.geojson")
    rows = G.readiness(b, Readiness(name="year_all_null.geojson")).rows
    out["year_all_null.geojson"] = {"year_field": b.fields.get("year"), "years": b.years,
                                    "row": [list(r) for r in rows if r[0] == "Install year"]}
    return out


def _canonical(net) -> set:
    """A network's topology with no ID or coordinate in it: each node named by the smallest (feature, part, end) of
    the line ends it holds, each pipe by its feature, part and piece and its two nodes."""
    label: dict = {}
    for key, n in sorted(net.end_node.items()):
        label.setdefault(n, key)
    # a node no line ends at (the middle of a closed line cut in two, since the review of task 17): the smallest
    # (feature, part, piece, end) of the pieces that end there
    for key, n in sorted(((p["feature"], p["part"], p["piece"], k), p[e]) for p in net.pipes
                         for k, e in ((0, "start"), (1, "end"))):
        label.setdefault(n, (-1,) + key)
    return {(p["feature"], p["part"], p["piece"], tuple(sorted((label.get(p["start"]), label.get(p["end"])))))
            for p in net.pipes}


def part_gis_invariance() -> dict:
    """gis_snapping_invariance: ky4's mains layer (feet) translated by 5e6 ft and rotated 37 degrees, in metres and in
    feet, built at 6 ft through the same frame code as a layer read from a file; and a skeleton written in GPM and feet
    and read back by WNTR."""
    import wntr
    from . import fixtures as FX
    from . import gis_read as G
    from . import topology as T
    recs = FX.ky4_mains_records()

    def layer(transform, unit):
        feats = [G.Feature(index=k, kind="line", parts=[np.array([transform(p) for p in r["pts"]], float)],
                           props={"PIPEDIAM": r["text"]}) for k, r in enumerate(recs)]
        return G.Layer(name="ky4_mains", source="geojson", features=feats, fields=["PIPEDIAM"],
                       crs={"kind": "projected", "unit": unit, "unit_m": T.UNIT_M[unit], "how": "the check"})
    th = math.radians(37.0)
    rot = lambda p: (5e6 + p[0] * math.cos(th) - p[1] * math.sin(th), 5e6 + p[0] * math.sin(th) + p[1] * math.cos(th))
    base = G.build_from_layer(layer(lambda p: p, "ft"))
    moved = G.build_from_layer(layer(rot, "ft"))
    metres = G.build_from_layer(layer(lambda p: (p[0] * T.FT, p[1] * T.FT), "m"))
    out = {"translated_rotated_same_topology": _canonical(base.net) == _canonical(moved.net),
           "pipes": len(base.net.pipes), "metres_same_topology": _canonical(base.net) == _canonical(metres.net)}
    lab_b = {(p["feature"], p["part"], p["piece"]): p for p in base.net.pipes}
    lab_m = {(p["feature"], p["part"], p["piece"]): p for p in metres.net.pipes}
    rel = [abs(lab_b[k]["length_m"] - lab_m[k]["length_m"]) / lab_b[k]["length_m"] for k in lab_b]
    dia = [abs(base.diam_in[lab_b[k]["id"]] - metres.diam_in[lab_m[k]["id"]]) for k in lab_b]
    out["metres_vs_feet_max_rel_length"] = max(rel)
    out["metres_vs_feet_max_diameter_diff_in"] = max(dia)
    with tempfile.TemporaryDirectory(prefix="rm_t17_inp_") as d:
        p = os.path.join(d, "skeleton.inp")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(base.inp_text)
        wn = wntr.network.WaterNetworkModel(p)
    rl = max(abs(wn.get_link(q["id"]).length - q["length_m"]) / q["length_m"] for q in base.net.pipes)
    rd = max(abs(wn.get_link(q["id"]).diameter - base.diam_in[q["id"]] * 0.0254) / (base.diam_in[q["id"]] * 0.0254)
             for q in base.net.pipes)
    out.update(skeleton_units=wn.options.hydraulic.inpfile_units, skeleton_max_rel_length=rl,
               skeleton_max_rel_diameter=rd, skeleton_pipes=len(wn.pipe_name_list),
               skeleton_junctions=len(wn.junction_name_list))
    return out


def part_gis_chain(paths: dict) -> dict:
    """gis_snapping_chain: a 5 ft plus 5 ft chain, a valve cluster and a pump manifold at the 6 ft default."""
    from . import fixtures as FX
    b = _gis(_read_bytes(paths["chain_cluster.geojson"]), "chain_cluster.geojson")
    r = b.net.report
    # which group each node's ends come from
    groups = {}
    feats = {f.index: f.props["PIPEID"].rsplit("_", 1)[0] for f in b.layer.features}
    for (f, _, end), n in b.net.end_node.items():
        if end == 0:
            groups.setdefault(feats[f], set()).add(n)
    tees = {}
    for n in FX.TEE_CASES:           # the review of task 17: tees on a main of dense vertices, an offset cross, a loop
        t = _gis(_read_bytes(paths[n]), n).net.report
        tees[n] = {k: t[k] for k in ("pipes", "nodes", "components", "dead_ends", "tees_split", "loops_cut",
                                     "dropped_pieces", "nodes_flagged_wide")}
    return {"nodes_flagged_wide": r["nodes_flagged_wide"], "flagged_detail": r["flagged_detail"],
            "nodes_split": r["nodes_split"], "split_detail": r["split_detail"],
            "nodes_per_group": {g: len(v) for g, v in sorted(groups.items())},
            "ends_per_group": {g: len(v) for g, v in FX.CHAIN_ENDS.items()},
            "warnings": [w for w in b.warnings if "tolerance" in w], "tees_and_loops": tees}


def part_gis_ids(paths: dict) -> dict:
    """gis_ids_deterministic: the lon/lat ky4 layer shuffled gives the same node and pipe IDs; dirty IDs go through
    the mapping table."""
    a = _gis(_read_bytes(paths["ky4_mains_lonlat.geojson"]), "ky4_mains_lonlat.geojson")
    b = _gis(_read_bytes(paths["ky4_mains_shuffled.geojson"]), "ky4_mains_shuffled.geojson")
    pa = {(p["id"], p["start"], p["end"]) for p in a.net.pipes}
    pb = {(p["id"], p["start"], p["end"]) for p in b.net.pipes}
    d = _gis(_read_bytes(paths["ids_dirty.geojson"]), "ids_dirty.geojson")
    longest = max(len(p["id"].encode("utf-8")) for p in a.net.pipes + d.net.pipes)
    longest_node = max(len(n.encode("utf-8")) for n in list(a.net.nodes) + list(d.net.nodes))
    return {"shuffled_same_node_ids": set(a.net.nodes) == set(b.net.nodes),
            "shuffled_same_pipes": pa == pb, "nodes": len(a.net.nodes), "pipes": len(pa),
            "inp_identical": a.inp_text == b.inp_text, "dirty_mapping": d.id_table,
            "dirty_model_ids": sorted(p["id"] for p in d.net.pipes), "longest_id_bytes": max(longest, longest_node)}


def part_gis_diameters(paths: dict) -> dict:
    """diameter_parsing: every row of diameter_cases.csv, and the layer-level bare-number ranges."""
    import csv as _csv
    from . import fixtures as FX
    from . import topology as T
    rows, wrong = [], []
    with open(paths["diameter_cases.csv"], newline="", encoding="utf-8") as fh:
        for r in _csv.DictReader(fh):
            got = T.parse_diameter(r["text"])
            want_v = float(r["value"]) if r["value"] else None
            ok = (got["value"] == want_v and (got["unit"] or "") == r["unit"] and (got["flag"] or "") == r["flag"]
                  if r["flag"] not in ("empty", "unreadable", "reducer", "two sizes") else
                  got["value"] is None and got["flag"] == r["flag"])
            rows.append({"text": r["text"], "value": got["value"], "unit": got["unit"], "flag": got["flag"], "ok": ok})
            if not ok:
                wrong.append(r["text"])
    bare = [{"values": list(v), "unit": T.bare_unit(v), "ok": T.bare_unit(v) == u} for v, u in FX.BARE_CASES]
    return {"cases": len(rows), "right": len(rows) - len(wrong), "wrong": wrong, "rows": rows, "bare": bare}


# the review of task 17: (fixture, options, read with the feature cap at fixtures.REVIEW_MAX_FEATURES)
REVIEW_REFUSAL_CASES = (("bad_vertices.gpkg", {}, False), ("empty_lines_4326.gpkg", {}, False),
                        ("lonlat_no_prj.zip", {}, False), ("lonlat_no_prj.zip", {"unit": "ft"}, False),
                        ("lonlat_no_prj.zip", {"unit": "deg"}, False), ("lonlat_as_3857.geojson", {}, False),
                        ("stateplane_as_4326.geojson", {}, False), ("stateplane_as_4326.geojson", {"unit": "us_ft"}, False),
                        ("world_mercator.zip", {}, False), ("many_shapefiles.zip", {}, False),
                        ("esri_over_limit.json", {}, True), ("mains_and_many_points.gpkg", {}, True))


@contextlib.contextmanager
def _feature_cap(n: int | None):
    """gis_read.MAX_FEATURES set to n for a while (the over-limit fixtures are small, read with a small cap)."""
    from . import gis_read as G
    old = G.MAX_FEATURES
    G.MAX_FEATURES = n or old
    try:
        yield
    finally:
        G.MAX_FEATURES = old


def _case_key(name: str, opts: dict, capped: bool) -> str:
    from . import fixtures as FX
    return name + "".join(f"|{k}={v}" for k, v in sorted(opts.items())) + (
        f"|cap={FX.REVIEW_MAX_FEATURES}" if capped else "")


def part_gis_refusals(paths: dict) -> dict:
    """gis_readers_refuse_cleanly: each file gets a question or a plain message, and none raises; since the review of
    task 17 also the layers that broke the reader then (REVIEW_REFUSAL_CASES), each read as build_gis reads it and as
    the Get started page does (read_layers, then build_from_layer)."""
    from . import fixtures as FX
    from . import gis_read as G
    out = {}
    for n in ("no_prj.zip", "utility_gdb.zip", "points_only.geojson", "empty_layer.gpkg"):
        try:
            b = _gis(_read_bytes(paths[n]), n)
            out[n] = {"raised": False, "refusal": b.refusal, "questions": [q["key"] for q in b.questions],
                      "built": b.inp_text is not None}
        except Exception as e:  # noqa: BLE001
            out[n] = {"raised": True, "error": f"{type(e).__name__}: {e}"}
    for n, opts, capped in REVIEW_REFUSAL_CASES:
        key = _case_key(n, opts, capped)
        data = _read_bytes(paths[n])
        with _feature_cap(FX.REVIEW_MAX_FEATURES if capped else None):
            try:
                b = _gis(data, n, **opts)
                r = b.net.report if b.net is not None else None
                out[key] = {"raised": False, "refusal": b.refusal, "questions": [q["key"] for q in b.questions],
                            "question_options": [q["options"] for q in b.questions],
                            "question_text": [q["text"] for q in b.questions], "built": b.inp_text is not None,
                            "pipes": r["pipes"] if r else None, "components": r["components"] if r else None,
                            "layers": b.layers, "warnings": b.warnings, "notes": b.notes}
            except Exception as e:  # noqa: BLE001
                out[key] = {"raised": True, "error": f"{type(e).__name__}: {e}"}
            try:                      # the Get started page's own path
                layers, kind = G.read_layers(data, n)
                lay = next((x for x in layers if x.n_lines), layers[0])
                pb = G.build_from_layer(lay, G.GisOptions(layer=lay.name, **opts), n, kind, layers)
                out[key]["page_path"] = {"raised": False, "refusal": pb.refusal, "built": pb.inp_text is not None}
            except G.GisRefused as e:
                out[key]["page_path"] = {"raised": False, "refusal": str(e), "built": False}
            except Exception as e:  # noqa: BLE001
                out[key]["page_path"] = {"raised": True, "error": f"{type(e).__name__}: {e}"}
    return out


# the final verification of task 17: the inlet residuals headers it listed, each read or refused (free chlorine)
FINAL_INLET_HEADERS = ("Cl", "FAC", "Free", "FCl", "Res", "Resid", "mg_l", "mg/L", "ppm", "residual_mgL", "value",
                       "reading", "conc", "concentration", "average", "median", "level", "measured", "dose", "setpoint",
                       "target", "Total", "TRC", "monochloramine", "NH2Cl")
FINAL_SIZE_TEXTS = ("9" * 400, "9" * 400 + " in", '8"x8"x6"', '8"x8"x8"', '6"x6"', "12x8 RED", "6x6x4", "12x12x8 TEE",
                    "6x4x6", "8x8x8", '6" x 6" x 4"')
FINAL_STRICT_TEXTS = ("Line 3 in Oak St", "8 in main", "8 in", "8 in.", '8"', "8 inch main", "200 mm")
FINAL_TEE_SWEEP = {"branches": (2, 3, 4, 5), "spacing_ft": (3.0, 4.0, 5.0, 5.5)}


def _final_case(data: bytes, n: str, opts: dict | None = None, capped: bool = False) -> dict:
    """One fixture read as build_gis reads it (with the hosted demo's unpacked limit) and as the Get started page
    does (read_layers, then build_from_layer): the plain outcome, and whether either path raised."""
    from . import cloud
    from . import fixtures as FX
    from . import gis_read as G
    opts = opts or {}
    lim = cloud.UPLOAD_MB * 1_000_000
    with _feature_cap(FX.REVIEW_MAX_FEATURES if capped else None):
        try:
            b = G.build_gis(data, n, G.GisOptions(**opts), max_unpacked=lim)
            r = b.net.report if b.net is not None else None
            out = {"raised": False, "refusal": b.refusal, "questions": [q["key"] for q in b.questions],
                   "built": b.inp_text is not None, "pipes": r["pipes"] if r else None,
                   "components": r["components"] if r else None, "offer_unit": b.offer_unit, "notes": b.notes,
                   "warnings": b.warnings, "layers": b.layers}
        except Exception as e:  # noqa: BLE001
            return {"raised": True, "error": f"{type(e).__name__}: {e}"}
        try:
            layers, kind = G.read_layers(data, n, lim)
            lay = next((x for x in layers if x.name == opts.get("layer")), None) or next(
                (x for x in layers if x.n_lines), layers[0])        # the page offers only the layers it read
            pb = G.build_from_layer(lay, G.GisOptions(**{**opts, "layer": lay.name}), n, kind, layers)
            out["page_path"] = {"raised": False, "built": pb.inp_text is not None, "refusal": pb.refusal,
                                "layer": lay.name}
        except G.GisRefused as e:
            out["page_path"] = {"raised": False, "built": False, "refusal": str(e)}
        except Exception as e:  # noqa: BLE001
            out["page_path"] = {"raised": True, "error": f"{type(e).__name__}: {e}"}
    return out, b


def _tee_sweep() -> dict:
    """The verifier's sweep: n branches ending 4 ft short of one main, s ft apart along it, fanned out, on a straight
    two-point main and on a main with a vertex every 2 ft, at 6 ft: the connected parts of each."""
    from . import topology as T
    from .fixtures import _branches
    out = {}
    for main_kind, main in (("straight", [(-300.0, 0.0), (300.0, 0.0)]),
                            ("dense", [(float(-300 + 2 * k), 0.0) for k in range(301)])):
        for n in FINAL_TEE_SWEEP["branches"]:
            for sp in FINAL_TEE_SWEEP["spacing_ft"]:
                lines = [main] + _branches(n, sp)
                L = [T.Line(k, 0, np.array(ln, float) * T.FT, np.array(ln, float)) for k, ln in enumerate(lines)]
                r = T.build_network(L, 6 * T.FT).report
                out[f"{main_kind}|n={n}|s={sp:g}"] = r["components"]
    return out


def part_gis_final(paths: dict) -> dict:
    """gis_final_verification: each issue the final verification of task 17 confirmed, as a fixture or a case, read
    as build_gis reads it (with the hosted demo's unpacked limit) and as the page does.  Times are measured by the
    check, not recorded."""
    from . import fixtures as FX
    from . import gis_read as G
    from . import importers as I
    from . import topology as T
    out = {}

    def case(n, opts=None, capped=False, key=None):
        res, b = (lambda r: (r, None) if isinstance(r, dict) else r)(_final_case(_read_bytes(paths[n]), n, opts, capped))
        out[key or _case_key(n, opts or {}, capped)] = res
        return b
    for n in ("kmz_long_whitespace.kmz", "kmz_long_lt.kmz", "kmz_open_rows.kmz", "kmz_doctype_after_blank_lines.kmz",
              "kmz_doctype_after_comment.kmz", "prj_spaces_comma.zip", "prj_spaces_comma_20kb.zip",
              "esri_wkt_spaces_comma.json", "world_mercator_3395.geojson", "esri_wkid_3395.json",
              "kml_line_in_street.kml", "zipped_gpkg.zip", "zipped_geojson.zip", "zipped_two_geojson.zip",
              "geojson_crs_props_string.geojson", "esri_odd_attributes.json", "tee_fan_dense_main.geojson",
              "tee_fan_straight_main.geojson"):
        case(n)
    for n in ("mains_and_service_lines.gpkg", "mains_and_service_lines.zip"):
        case(n, capped=True)
        case(n, {"layer": "mains"}, capped=True)
        case(n, {"layer": "service_lines"}, capped=True)
    for n in ("lonlat_stateplane_prj.zip", "lonlat_as_2263.geojson"):
        case(n)
        case(n, {"unit": "deg"})
    # swapped axes: asked; swapped, the lengths of the same grid written longitude first; kept, the lengths near the pole
    truth = _gis(FX._fc(FX._grid_lines(*FX.FINAL_SWAP_GRID), [{"PIPEDIAM": "8 in"}] * 24,
                        "urn:ogc:def:crs:EPSG::4326"), "truth.geojson")
    want = sorted(p["length_m"] for p in truth.net.pipes)
    for n in ("swapped_latlon_4326.geojson", "swapped_latlon_nocrs.geojson"):
        case(n)
        for sw in (True, False):
            b = case(n, {"swap_axes": sw})
            got = sorted(p["length_m"] for p in b.net.pipes)
            out[_case_key(n, {"swap_axes": sw}, False)].update(
                lengths_m=[round(min(got), 1), round(max(got), 1)], true_lengths_m=[round(min(want), 1),
                                                                                    round(max(want), 1)],
                max_rel_error=max(abs(g - w) / w for g, w in zip(got, want)))
    b = case("loop_joined_through_feed.geojson")
    out["loop_joined_through_feed.geojson"].update({k: b.net.report[k] for k in ("loops_cut", "dropped_pieces")})
    b = case("esri_epoch_dates.json")
    out["esri_epoch_dates.json"]["years"] = b.years
    b = case("gpkg_bad_utf8.gpkg")
    out["gpkg_bad_utf8.gpkg"]["first_size_text"] = str(b.layer.features[0].props.get("PIPEDIAM"))
    b = case("size_overflow.geojson")
    out["size_overflow.geojson"].update(diameters=b.diam_counts, inf_written="inf" in b.inp_text.lower())
    out["kml_line_in_street.kml"]["diameter_field"] = _gis(_read_bytes(paths["kml_line_in_street.kml"]),
                                                           "kml_line_in_street.kml").fields.get("diameter")
    out["wkt_unit_of_the_20kb_run"] = G.wkt_unit(FX._prj_run(FX.FINAL_PRJ_RUN[1]))
    out["size_texts"] = {t[:12] + (f"...({len(t)} characters)" if len(t) > 12 else ""): T.parse_diameter(t)["flag"]
                         or T.parse_diameter(t)["value"] for t in FINAL_SIZE_TEXTS}
    out["strict_in_texts"] = {t: [T.parse_diameter(t, strict_in=True)["unit"], T.parse_diameter(t)["unit"]]
                              for t in FINAL_STRICT_TEXTS}
    out["tee_sweep_components"] = _tee_sweep()
    inlet = {}
    for h in FINAL_INLET_HEADERS:
        try:
            I.read_inlet_residuals(f"inlet,{h}\nRiver,1.20\nLake,0.85\n".encode())
            inlet[h] = "read"
        except I.LogRefused as e:
            inlet[h] = "refused: " + str(e)[:90]
    out["inlet_headers"] = inlet
    return out


def part_gis_valves(paths: dict) -> dict:
    """gis_valve_status: a valve layer's NORMAL_POS gives proposed closed_valve rows, each at its node when one is
    within the tolerance."""
    from . import gis_read as G
    b = _gis(_read_bytes(paths["valves_normal_pos.geojson"]), "valves_normal_pos.geojson")
    v = G.proposed_closed_valves(b.layer, b.net, frame=b.frame)
    rows = lambda v: [{"valve": r["valve"], "status": r["status"], "node": r["node"],
                       "distance_m": None if r["distance_m"] is None else round(r["distance_m"], 4)} for r in v["rows"]]
    # the same mains and valves as two layers of one GeoPackage, the valve layer found beside the line layer
    layers, _ = G.read_layers(_read_bytes(paths["mains_and_valves.gpkg"]), "mains_and_valves.gpkg")
    b2 = _gis(_read_bytes(paths["mains_and_valves.gpkg"]), "mains_and_valves.gpkg")
    vl = G.valve_layer(layers, b2.layer)
    v2 = G.proposed_closed_valves(vl, b2.net, frame=b2.frame)
    return {"status_field": v["status_field"], "id_field": v["id_field"], "proposed": rows(v),
            "open": v["open"], "unread": v["unread"],
            "separate_layer": {"line_layer": b2.layer.name, "valve_layer": vl.name, "proposed": rows(v2)}}


def _wntr_runs(path: str) -> dict:
    """Pressure, flow, tank head and water age (7 days, as every run is set up) and chlorine at WNTR_MEMBER."""
    import wntr
    from .simulate import load, simulate_nominal_chlorine
    with tempfile.TemporaryDirectory(prefix="rm_t17_wntr_") as d:
        wn = load(path)
        wn.options.quality.parameter = "AGE"
        res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "h"))
        c = simulate_nominal_chlorine(path, *WNTR_MEMBER, file_prefix=os.path.join(d, "c"))
    tanks = wn.tank_name_list
    return {"pressure": res.node["pressure"], "flow": res.link["flowrate"],
            "head": res.node["head"][tanks] if tanks else None, "age": res.node["quality"], "chlorine": c}


def _max_diff(a: dict, b: dict) -> dict:
    out = {}
    for k in a:
        if a[k] is None:
            out[k] = None
            continue
        x, y = a[k], b[k][a[k].columns]
        out[k] = float(np.nanmax(np.abs(x.values - y.values))) if x.size else 0.0
    return out


def part_gis_wntr(paths: dict) -> dict:
    """gis_wntr_schema_roundtrip: Net3, Net2 and ky4 in WNTR's GIS schema with WNTR's JSON rebuild with 0.0 difference
    against the bundled file as WNTR reads it and writes it again (the same writer), and the negative control: without
    copying demand_pattern into pattern_name, Net3's chlorine moves by more than NEGATIVE_BAR_MGL."""
    import wntr
    from . import fixtures as FX
    from . import gis_read as G
    from .simulate import LIB
    out = {"member": {"kb_per_day": WNTR_MEMBER[0], "kw_m_per_day": WNTR_MEMBER[1], "gamma": WNTR_MEMBER[2]},
           "networks": {}}
    with tempfile.TemporaryDirectory(prefix="rm_t17_schema_") as d:
        for net in FX.WNTR_SCHEMA_NETS:
            data = _read_bytes(paths[f"wntr_schema_{net}.zip"])
            ref = os.path.join(d, f"{net}_ref.inp")
            wntr.network.write_inpfile(wntr.network.WaterNetworkModel(os.path.join(LIB, f"{net}.inp")), ref)
            got = os.path.join(d, f"{net}_gis.inp")
            with open(got, "w", encoding="utf-8") as fh:
                fh.write(G.wntr_schema_inp(data, fix_patterns=True))
            a, b = _wntr_runs(ref), _wntr_runs(got)
            entry = {"max_abs_diff": _max_diff(a, b)}
            if net == "Net3":
                bad = os.path.join(d, f"{net}_nofix.inp")
                with open(bad, "w", encoding="utf-8") as fh:
                    fh.write(G.wntr_schema_inp(data, fix_patterns=False))
                c = _wntr_runs(bad)
                diff = np.abs(a["chlorine"].values - c["chlorine"][a["chlorine"].columns].values)
                entry["without_pattern_workaround"] = {"chlorine_max_abs_diff_mgL": float(diff.max()),
                                                       "chlorine_mean_abs_diff_mgL": float(diff.mean()),
                                                       "above_bar": float(diff.max()) > NEGATIVE_BAR_MGL}
            out["networks"][net] = entry
    return out


def part_gis_external() -> dict:
    """gis_external_file: the public open-data mains layer was not used (lead decision 8: fixtures are generated
    locally), so the check is skipped and says so."""
    return {"skipped": True, "why": "no external file: the lead's decision for this iteration is that fixtures are "
                                    "generated locally, so no public layer was downloaded"}


def gis_record() -> dict:
    from . import fixtures as FX
    _quiet()
    with tempfile.TemporaryDirectory(prefix="rm_t17_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.GIS_FIXTURES)
        return {"about": "Task 17's GIS importer, converter self-consistency on public models converted to GIS by our "
                         "own writers (residualmap/fixtures.py), never a test on a utility's own layer (R6): ky4's mains "
                         "rebuilt from four utility-style layers against the pre-registered bar (part 'rebuild'), ky4 "
                         "drawn by hand in Google Earth (part 'kml', no bar), the file-format quirks, snapping "
                         "invariance and chains, deterministic IDs, diameter text, clean refusals, a valve layer's "
                         "status, WNTR's own GIS schema with its JSON (part 'wntr'), and the cases of the task's final verification "
                         "(part 'final'). Written by python -m "
                         "residualmap.intake_exp gis; the task17 checks in residualmap/checks.py rerun every part and "
                         "compare.",
                "rebuild": part_gis_rebuild(paths), "kml": part_gis_kml(paths), "quirks": part_gis_quirks(paths),
                "invariance": part_gis_invariance(), "chain": part_gis_chain(paths), "ids": part_gis_ids(paths),
                "diameters": part_gis_diameters(paths), "refusals": part_gis_refusals(paths),
                "valves": part_gis_valves(paths), "wntr": part_gis_wntr(paths), "external": part_gis_external(),
                "final": part_gis_final(paths)}


# ----------------------------------------------------------------------------- task 18: cad_roundtrip.json
RECORD_18 = os.path.join(OUT, "cad_roundtrip.json")
CAD_BAR = {"tolerance_ft": 6.0, "pipes_at_least": 1151, "of_pipes": 1156, "components": 1,
           "diameter_share_at_least": 0.99, "fixture": "ky4_r12_handwritten.dxf"}   # pre-registered (plan task 18)
CAD_SCORED = ("ky4_r12_handwritten.dxf", "ky4_r2018_lwpolyline.dxf")
CAD_READERS = {"ky4_r12_handwritten.dxf": "ezdxf (recover mode), on an R12 drawing written by the fixtures' own "
                                          "writer: cross-implementation",
               "ky4_r2018_lwpolyline.dxf": "ezdxf, on an R2018 drawing ezdxf wrote: same-library self-consistency"}


def _cad(data: bytes, name: str, **kw):
    from . import cad_read as C
    return C.build_cad(data, name, C.CadOptions(**kw))


def _cad_score_keys(s: dict) -> dict:
    return {k: s[k] for k in ("pipes_rebuilt", "of_pipes", "diameters_matched", "diameter_share", "extra_pipes",
                              "components", "nodes", "dead_ends", "tees_split", "nodes_flagged_wide", "nodes_split",
                              "dropped_pieces", "pipes_built", "missed", "collapsed")}


def _labels_elsewhere(b, plan) -> list:
    """The labels of the ky4 drawing that went to a line other than their own: the pipe they label, the pipe they went
    to, their distance from each (ft) and their own line's longest segment (ft)."""
    from . import cad_read as C
    from . import fixtures as FX
    from .topology import FT
    d, recs = b.drawing, FX.ky4_mains_records()
    origin = np.asarray(b.frame["origin"], float)
    f = C.UNIT_M[b.unit]
    lines = [(ln.pts - origin) * f for ln in d.lines]
    texts = [tx for tx in d.texts if C.label_size(tx.text)["value"] is not None]
    got = C.assign_labels(lines, [(tx.box - origin) * f for tx in texts],
                          [b.labels["radius_text_heights"] * tx.height * f for tx in texts])
    owners = [j for j, p in enumerate(plan) if p["label"]]
    out = []
    for g, j, tx in zip(got, owners, texts):
        if g["line"] == j:
            continue
        own = float(C._box_seg_dist(np.repeat(((tx.box - origin) * f)[None], len(lines[j]) - 1, axis=0),
                                    lines[j][:-1], lines[j][1:]).min())
        out.append({"label_of": recs[j]["pipes"][0], "text": tx.text,
                    "went_to": recs[g["line"]]["pipes"][0] if g["line"] is not None else None,
                    "distance_ft": round(g["distance"] / FT, 3) if g["distance"] is not None else None,
                    "own_distance_ft": round(own / FT, 3),
                    "own_longest_segment_ft": round(max(FX._seg_len(recs[j]["pts"])), 1)})
    return out


def part_cad_rebuild(paths: dict) -> dict:
    """dxf_rebuilds_network: ky4's mains drawn as an R12 drawing by the fixtures' own writer (the bar) and as an R2018
    drawing by ezdxf (reported), rebuilt at the 6 ft default and scored with task 17's rule."""
    from . import fixtures as FX
    from .topology import FT
    out = {"bar": CAD_BAR, "fixtures": {}}
    plan = FX.ky4_cad_plan(False)
    out["drawn"] = {"lines": len(plan), "sizes_in_layer_names": sum(1 for p in plan if p["label"] is None),
                    "sizes_as_labels": sum(1 for p in plan if p["label"] is not None),
                    "text_height_ft": FX.CAD_TEXT_H, "label_offsets_heights": list(FX.CAD_OFFSETS_H)}
    for name in CAD_SCORED:
        b = _cad(_read_bytes(paths[name]), name)
        assert b.net is not None and b.inp_text, (name, b.refusal, b.questions)
        s = score_ky4(b)
        d = b.drawing
        out["fixtures"][name] = {
            "reader": CAD_READERS[name], "version": d.counts.get("version_name"), "lines": len(d.lines),
            "texts": len(d.texts), "unit": b.unit, "default_tolerance_ft": b.tolerance_m / FT,
            "score_6ft": _cad_score_keys(s), "size_counts": b.size_counts, "labels": b.labels,
            "c_counts": b.c_counts, "labels_on_another_line": _labels_elsewhere(b, plan if name == CAD_SCORED[0]
                                                                                else FX.ky4_cad_plan(True)),
            "default_sized": [{"pipe": FX.ky4_mains_records()[p["feature"]]["pipes"][0], "piece": p["piece"],
                               "true_in": FX.ky4_mains_records()[p["feature"]]["diameter_in"],
                               "given_in": b.diam_in[p["id"]]}
                              for p in b.net.pipes if b.size_source[p["id"]] == "default"],
            "bar_met": (s["pipes_rebuilt"] >= CAD_BAR["pipes_at_least"]
                                                and s["components"] == CAD_BAR["components"]
                                                and s["diameter_share"] >= CAD_BAR["diameter_share_at_least"])}
    out["bar_met"] = out["fixtures"][CAD_BAR["fixture"]]["bar_met"]
    return out


def part_cad_tolerance(paths: dict) -> dict:
    """dxf_tolerance_table_feet: ky4's and ky10's drawings at 0.5, 1, 3, 6 and 12 ft: pipes rebuilt, connected parts,
    dead ends and real pipes collapsed (both of a pipe's nodes in one rebuilt node), and ky10's node pairs drawn at
    one spot, which no drawing of pipes can separate."""
    from . import cad_read as C
    from . import fixtures as FX
    from .topology import FT
    out = {"feet": list(C.DIAGNOSTIC_FT), "networks": {}}
    for net, name, man in (("ky4", "ky4_r12_handwritten.dxf", None), ("ky10", "ky10_r12.dxf", FX.ky10_manifest())):
        data = _read_bytes(paths[name])
        rows = {}
        for ft in C.DIAGNOSTIC_FT:
            b = _cad(data, name, tolerance_m=ft * FT)
            assert b.net is not None, (name, ft, b.refusal)
            s = score_ky4(b, man=man)
            if ft == 6.0:
                sizes_6 = {"size_counts": b.size_counts, "labels": b.labels, "diameters_matched": s["diameters_matched"],
                           "diameter_share": s["diameter_share"]}
            rows[f"{ft:g}"] = {"pipes_rebuilt": s["pipes_rebuilt"], "of_pipes": s["of_pipes"],
                               "components": s["components"], "dead_ends": s["dead_ends"],
                               "real_pipes_collapsed": len(s["collapsed"]), "collapsed": s["collapsed"],
                               "extra_pipes": s["extra_pipes"], "pipes_built": s["pipes_built"],
                               "dropped_pieces": s["dropped_pieces"], "tees_split": s["tees_split"]}
        out["networks"][net] = {"fixture": name, "by_tolerance_ft": rows, "sizes_at_6ft": sizes_6}
    man = FX.ky10_manifest()
    pairs = man["same_spot_pairs"]
    joined = {}
    data = _read_bytes(paths["ky10_r12.dxf"])
    for ft in C.DIAGNOSTIC_FT:
        b = _cad(data, "ky10_r12.dxf", tolerance_m=ft * FT)
        votes = {}
        for j, r in enumerate(man["records"]):
            for end in (0, 1):
                n = b.net.end_node.get((j, 0, end))
                if n is not None:
                    votes.setdefault(r["ends"][end], {}).setdefault(n, 0)
                    votes[r["ends"][end]][n] += 1
        node = {u: max(c.items(), key=lambda x: (x[1], x[0]))[0] for u, c in votes.items()}
        joined[f"{ft:g}"] = sum(1 for p in pairs if p["nodes"][0] in node and node.get(p["nodes"][0]) ==
                                node.get(p["nodes"][1]))
    out["ky10_same_spot_pairs"] = {"pairs": len(pairs), "kinds": sorted({(p["link"] or "none").split()[0] for p in pairs}),
                                   "list": pairs, "in_one_node_by_tolerance_ft": joined,
                                   "unrecoverable_from_geometry": True}
    return out


def _sizes_of_feature(b, f: int) -> list:
    return [(b.size_source[p["id"]], round(b.diam_in[p["id"]], 3)) for p in b.net.pipes if p["feature"] == f]


def part_cad_labels(paths: dict) -> dict:
    """dxf_labels_and_propagation: labels_offset_cases.dxf, each line's size and where it came from, against what
    the fixture says it should be."""
    from . import fixtures as FX
    b = _cad(_read_bytes(paths["labels_offset_cases.dxf"]), "labels_offset_cases.dxf")
    assert b.net is not None, (b.refusal, b.questions)
    sizes = [v for k, v in b.diam_in.items() if b.size_source[k] != "default"]
    default = max(set(round(v, 2) for v in sizes), key=lambda s: (sum(1 for v in sizes if round(v, 2) == s), -s))
    cases, f = {}, 0
    for c in FX.LABEL_CASES:
        got, right = [], True
        for want_src, want_size in c["want"]:
            g = _sizes_of_feature(b, f)
            f += 1
            want = want_size if want_size is not None else default
            ok = len(g) == 1 and g[0][0] == want_src and abs(g[0][1] - want) <= 0.01 * want
            right &= ok
            got.append({"source": g[0][0] if g else None, "size_in": g[0][1] if g else None,
                        "want_source": want_src, "want_size_in": round(want, 3), "right": ok})
        cases[c["case"]] = {"lines": got, "right": right}
    return {"cases": cases, "right": sum(v["right"] for v in cases.values()), "of": len(cases),
            "size_counts": b.size_counts, "labels": b.labels, "default_size_in": default}


def part_cad_blocks(paths: dict) -> dict:
    """dxf_blocks_and_proxies: xref_block_mains.dxf's mains inside blocks (nested, rotated, scaled, a block array)
    found where AutoCAD draws them, on the layer it shows them on, the external xref named; Civil 3D objects get the
    AECTOACAD message, alone (refused) and beside plain lines (built, warned)."""
    from . import cad_read as C
    from . import fixtures as FX
    d = C.read_drawing(_read_bytes(paths["xref_block_mains.dxf"]), "xref_block_mains.dxf")
    b = C.build_from_drawing(d, C.CadOptions())
    want = FX.xref_block_expected()
    used, found, worst = set(), [], 0.0
    for w in want:
        W = np.asarray(w["pts"], float)
        best = None
        for ln in d.lines:
            if ln.index in used or len(ln.pts) != 2:
                continue
            e = min(float(np.abs(ln.pts - W).max()), float(np.abs(ln.pts[::-1] - W).max()))
            if best is None or e < best[0]:
                best = (e, ln)
        ok = best is not None and best[0] <= 1e-6 and best[1].layer == w["layer"]
        if ok:
            used.add(best[1].index)
            worst = max(worst, best[0])
            sz = _sizes_of_feature(b, best[1].index)
            ok = bool(sz) and all(abs(x[1] - w["size_in"]) <= 1e-9 for x in sz)
        found.append(ok)
    proxy = C.build_cad(_read_bytes(paths["civil3d_proxy.dxf"]), "civil3d_proxy.dxf")
    mixed = C.build_cad(_read_bytes(paths["civil3d_proxy_with_lines.dxf"]), "civil3d_proxy_with_lines.dxf")
    return {"block_lines_expected": len(want), "block_lines_found_in_place": sum(found),
            "worst_position_error_ft": worst, "lines_from_blocks": d.counts.get("lines_from_blocks", 0),
            "inserts_exploded": d.counts.get("inserts_exploded", 0), "xrefs": d.xrefs,
            "xref_warned": any(FX.XREF_PATH in w and "XSITE" in w for w in d.warnings),
            "civil3d_proxy.dxf": {"refused": proxy.refusal is not None, "aectoacad": "AECTOACAD" in (proxy.refusal or ""),
                                  "objects": proxy.drawing.proxies + proxy.drawing.civil3d},
            "civil3d_proxy_with_lines.dxf": {"pipes": len(mixed.net.pipes) if mixed.net else None,
                                             "aectoacad_warned": any("AECTOACAD" in w for w in mixed.warnings)}}


def part_cad_styles(paths: dict) -> dict:
    """dxf_layer_style_listing: abandoned_linetype.dxf lists the dashed abandoned mains as their own layer style, and
    leaving that style out removes those lines."""
    from . import cad_read as C
    d = C.read_drawing(_read_bytes(paths["abandoned_linetype.dxf"]), "abandoned_linetype.dxf")
    listing = {C.style_label(k): {"lines": v["lines"], "texts": v["texts"]} for k, v in d.styles.items()}
    dashed = [k for k in d.styles if k[1].upper() == "DASHED"]
    all_ = C.build_from_drawing(d, C.CadOptions())
    less = C.build_from_drawing(d, C.CadOptions(exclude=dashed))
    return {"styles": listing, "dashed_styles": [C.style_label(k) for k in dashed],
            "pipes_with_all": len(all_.net.pipes), "pipes_without_dashed": len(less.net.pipes),
            "lines_left_out": len(all_.net.pipes) - len(less.net.pipes)}


def part_cad_units(paths: dict) -> dict:
    """dxf_units_and_mtext: a unitless drawing asks its unit (and builds once it is given); metres declared on a
    drawing from an imperial template is flagged; the same with no measurement setting is not (a known limit);
    millimetres declared on a drawing in feet is asked (its median line is implausible); MTEXT with formatting codes
    gives its sizes."""
    from . import fixtures as FX
    out = {}
    for name in ("insunits_zero.dxf", "insunits_wrong.dxf", "insunits_wrong_nomeasure.dxf", "insunits_mm.dxf"):
        b = _cad(_read_bytes(paths[name]), name)
        out[name] = {"questions": [(q["key"], q["blocking"]) for q in b.questions], "built": b.inp_text is not None,
                     "unit": b.unit, "insunits": b.drawing.insunits, "measurement": b.drawing.measurement,
                     "median_ft": round(b.unit_check["median_ft"], 3) if b.unit_check.get("median_ft") else None,
                     "flagged": bool(b.unit_check.get("measurement_disagrees") or b.unit_check.get("implausible")
                                     or any(q["key"] == "unit" for q in b.questions))}
    ft = _cad(_read_bytes(paths["insunits_zero.dxf"]), "insunits_zero.dxf", unit="ft")
    out["insunits_zero.dxf"]["built_in_feet"] = {"pipes": len(ft.net.pipes), "unit": ft.unit}
    m = _cad(_read_bytes(paths["mtext_codes.dxf"]), "mtext_codes.dxf")
    out["mtext_codes.dxf"] = {"texts": [t.text for t in m.drawing.texts],
                              "sizes": [_sizes_of_feature(m, k) for k in range(len(FX.MTEXT_CASES))],
                              "want": [c[4] for c in FX.MTEXT_CASES]}
    return out


def part_cad_bak_binary(paths: dict) -> dict:
    """dxf_bak_and_binary: the ky4 drawing as an ASCII DXF, the same bytes named .bak and the same drawing as binary DXF
    (the fixtures' own binary writer) give the same skeleton, byte for byte."""
    import hashlib
    from .importers import sniff
    out = {}
    for name in ("ky4_r12_handwritten.dxf", "ky4_drawing.dxf.bak", "ky4_binary.dxf"):
        data = _read_bytes(paths[name])
        b = _cad(data, name)
        s = sniff(data, name)
        out[name] = {"sniff": s.kind, "backup": s.backup, "bytes": len(data), "pipes": len(b.net.pipes),
                     "inp_sha1": hashlib.sha1(b.inp_text.encode("utf-8")).hexdigest()}
    out["identical"] = len({v["inp_sha1"] for v in out.values()}) == 1
    return out


def part_cad_guidance(paths: dict) -> dict:
    """dwg_dgn_ifc_get_guidance: a DWG, an AutoCAD backup holding a DWG, DGN V7 and V8 and an IFC file are identified
    by their content and get what to send instead, with no traceback."""
    from . import cad_read as C
    from . import fixtures as FX
    out = {}
    for name in FX.CAD_GUIDANCE:
        d = C.read_drawing(_read_bytes(paths[name]), name)
        r = d.refusal or ""
        out[name] = {"sniff": d.sniff_kind, "first_sentence": r.split(". ")[0] + ".",
                     "says_dxf": "DXF" in r, "oda_link": C.ODA_LINK in r, "rename_to_dwg": "rename it to .dwg" in r}
    return out


def part_cad_hostile(paths: dict) -> dict:
    """cad_reads_in_linear_time: each hostile drawing, read with the local and the hosted demo's limits: what was read
    or why it was refused (the check times each read)."""
    from . import cad_read as C
    from . import cloud
    from . import fixtures as FX
    nb, fb = FX.HOSTILE_NEG_MINSERT
    out = {"sizes": {"mtext_codes_chars": FX.HOSTILE_MTEXT_CHARS, "stacked_mtexts": FX.HOSTILE_STACK_TEXTS[0],
                     "stacked_mtext_chars": FX.HOSTILE_STACK_TEXTS[1], "bomb_depth": FX.HOSTILE_BOMB[0],
                     "bomb_fan": FX.HOSTILE_BOMB[1], "bomb_lines": FX.HOSTILE_BOMB[1] ** (FX.HOSTILE_BOMB[0] - 1),
                     "minsert_copies": FX.HOSTILE_MINSERT[0] * FX.HOSTILE_MINSERT[1],
                     "spline_points": FX.HOSTILE_SPLINE_POINTS, "long_name_chars": FX.HOSTILE_LAYER_CHARS,
                     "upload_limit_bytes": cloud.UPLOAD_MB * 1_000_000,
                     # since the review
                     "tall_label_lines": FX.HOSTILE_TALL[0], "tall_label_points": FX.HOSTILE_TALL[0] * FX.HOSTILE_TALL[1],
                     "tall_labels": FX.HOSTILE_TALL[2], "tall_label_height_ft": FX.HOSTILE_TALL[3],
                     "crowded_labels": FX.HOSTILE_CROWDED[2], "crowded_label_height_ft": FX.HOSTILE_CROWDED[3],
                     "negative_minsert_depth": nb, "negative_minsert_fan": fb,
                     "negative_minsert_lines": fb * (fb + 1) ** (nb - 1),
                     "bulge_vertices": FX.HOSTILE_BULGES[0],
                     "bulge_points": 1 + FX.HOSTILE_BULGES[0] * C.MAX_ARC_SEGMENTS,
                     "cloud_text_limit": cloud.CAD_MAX_TEXTS, "cloud_label_pairs_limit": cloud.CAD_MAX_LABEL_PAIRS,
                     "local_label_pairs_limit": C.MAX_LABEL_PAIRS,
                     # since the final verification
                     "long_handle_hex_digits": FX.HOSTILE_LONG_HANDLE[0], "zero_handle_lines": FX.HOSTILE_LONG_HANDLE[1],
                     "inf_bulge_vertices": FX.HOSTILE_INF_BULGES, "handle_hex_digits_max": 16}}
    for name in FX.HOSTILE_FIXTURES:
        data = _read_bytes(paths[name])
        row = {"bytes": len(data)}
        if name == "hostile_long_handle.dxf":       # the final verification: what renumbering makes of the file
            repaired, info = C._repair_ascii(data)
            row["repaired_bytes"], row["handles_renumbered"] = len(repaired), info["handles_renumbered"]
            row["repaired_over_input"] = round(len(repaired) / len(data), 3)
        for label, lim in (("local", None), ("cloud", cloud.cad_limits())):
            d = C.read_drawing(data, name, lim)
            row[label] = {"lines": len(d.lines), "texts": len(d.texts), "refused": d.refusal is not None,
                          "refusal": (d.refusal or "").split(". ")[0], "warnings": [w.split(". ")[0] for w in d.warnings],
                          "texts_too_long": d.counts.get("texts_too_long", 0)}
            if name == "hostile_inf_bulges.dxf":    # the final verification: left out before its arcs are drawn
                row[label].update(lines_not_finite=d.counts.get("lines_not_finite", 0), points=d.counts.get("points"))
            if name in ("hostile_tall_labels.dxf", "hostile_crowded_labels.dxf"):     # the review: the build's bound
                b = C.build_from_drawing(d, C.CadOptions(to_scale=True), lim)
                row[label]["build"] = {"refused": b.refusal is not None, "refusal": (b.refusal or "").split(": ")[0],
                                       "pipes": len(b.net.pipes) if b.net else None, "size_counts": b.size_counts,
                                       "labels_too_tall": b.labels.get("not_sizes", {}).get("too tall", 0),
                                       "label_pairs_within_reach": b.labels.get("pairs_within_reach")}
        out[name] = row
    return out


def part_cad_review(paths: dict) -> dict:
    """dxf_review_cases: one drawing for each case the review of task 18 found, what the reader and the build now give
    (and, where the case was ezdxf's own, what ezdxf alone gives)."""
    import io
    import logging
    from ezdxf import recover
    from . import cad_read as C
    from . import fixtures as FX
    from .topology import FT
    logging.getLogger("ezdxf").setLevel(logging.ERROR)
    out = {}

    def qs(b):
        return [[q["key"], q["blocking"]] for q in b.questions]

    # inches declared on a drawing made in feet: asked, both medians shown; in feet it is ky4's skeleton; true inches
    data = _read_bytes(paths["ky4_inches_declared.dxf"])
    b, ft, inch = (_cad(data, "ky4_inches_declared.dxf", **kw) for kw in ({}, {"unit": "ft"}, {"unit": "in"}))
    s = score_ky4(ft)
    t = _cad(_read_bytes(paths["inches_true.dxf"]), "inches_true.dxf")
    out["inches"] = {"ky4_inches_declared.dxf": {
        "insunits": b.drawing.insunits, "questions": qs(b), "built": b.inp_text is not None,
        "median_ft_read_as_inches": round(b.unit_check["median_ft"], 1),
        "median_ft_read_as_feet": round(b.unit_check["median_ft_if_feet"], 1),
        "in_feet": {"pipes_rebuilt": s["pipes_rebuilt"], "of_pipes": s["of_pipes"], "components": s["components"],
                    "diameter_share": s["diameter_share"]},
        "in_inches": {"pipes": len(inch.net.pipes), "dropped_pieces": inch.net.report["dropped_pieces"],
                      "median_pipe_ft": round(float(np.median([p["length_m"] for p in inch.net.pipes])) / FT, 1)}},
        "inches_true.dxf": {"questions": qs(t), "unit": t.unit, "median_ft": round(t.unit_check["median_ft"], 1),
                            "median_ft_read_as_feet": round(C._median_ft(t.drawing.lines, FT), 1),
                            "pipes": len(t.net.pipes)}}
    # numbers in layer names that are no sizes
    z = _read_bytes(paths["zone_layers.dxf"])
    runs = {"asked": {}, "not_sizes": {"layer_numbers": "none"}, "inches": {"layer_numbers": "in"},
            "layer_names_off": {"layer_sizes": False}}
    zone = {}
    for k, kw in runs.items():
        bz = _cad(z, "zone_layers.dxf", **kw)
        zone[k] = {"questions": qs(bz), "built": bz.inp_text is not None, "size_counts": bz.size_counts,
                   "sizes_in": sorted({round(v, 3) for v in bz.diam_in.values()}),
                   "labels_win": bz.labels.get("layers_whose_labels_win", []),
                   "warned_labels_win": any("the labels are used on that layer" in w for w in bz.warnings)}
    zone["layer_names"] = {n: {k: C.layer_size(n)[k] for k in ("value", "unit", "cued")}
                           for n in ("WATER-ZONE-2", "WATER-PHASE-2", "SHEET-12", "W-DIST-10", "W_8_PVC", "WATER-8IN",
                                     "WTR-200MM", "MAINS-8-DIA")}
    out["zone_layers.dxf"] = zone
    # bytes invalid in the drawing's encoding (cp1252 in an R2018 file)
    data = _read_bytes(paths["latin1_r2018.dxf"])
    doc, _ = recover.read(io.BytesIO(data))
    d = C.read_drawing(data, "latin1_r2018.dxf")
    names = [k[0] for k in d.styles]
    out["latin1_r2018.dxf"] = {"ezdxf_default_gives_surrogates": any(C._SURROGATES.search(str(x.dxf.layer)) for x in
                                                                     doc.modelspace()),
                               "layer_names": names, "utf8_writable": all(C._SURROGATES.search(n) is None for n in names),
                               "pipes": len(C.build_from_drawing(d, C.CadOptions()).net.pipes)}
    # handles drawn twice, or 0
    for n in ("dup_handles.dxf", "zero_handles.dxf"):
        data = _read_bytes(paths[n])
        try:
            alone = len(recover.read(io.BytesIO(data))[0].modelspace())
        except Exception as e:  # noqa: BLE001  (ezdxf refuses the file)
            alone = type(e).__name__
        bh = _cad(data, n)
        out[n] = {"ezdxf_alone_lines": alone, "lines": len(bh.drawing.lines), "refused": bh.refusal is not None,
                  "handles_renumbered": bh.drawing.counts.get("handles_renumbered", 0),
                  "pipes": len(bh.net.pipes) if bh.net else None,
                  "components": bh.net.report["components"] if bh.net else None,
                  "noted": any("given new ones" in x for x in bh.drawing.notes)}
    # a spline-fit POLYLINE's control frame
    d = C.read_drawing(_read_bytes(paths["spline_fit_polyline.dxf"]), "spline_fit_polyline.dxf")
    ln = d.lines[0]
    curve = np.array(FX.spline_fit_curve(), float)
    with_frame = np.vstack([np.array(FX.SPLINE_FRAME, float), curve])
    length = lambda P: round(float(np.hypot(*(P[1:] - P[:-1]).T).sum()), 1)
    out["spline_fit_polyline.dxf"] = {"points_read": len(ln.pts), "length_read_ft": length(ln.pts),
                                      "curve_length_ft": length(curve), "length_with_frame_ft": length(with_frame)}
    # layers AutoCAD hides
    d = C.read_drawing(_read_bytes(paths["frozen_layers.dxf"]), "frozen_layers.dxf")
    fz = {"hidden_layers": d.hidden_layers, "default_excluded": [C.style_label(k) for k in C.default_excluded(d)],
          "in_frozen_blocks": d.counts.get("in_frozen_blocks", {})}
    for k, kw in (("default", {}), ("frozen_blocks", {"frozen_blocks": True}), ("all_hidden", {"include_hidden": True})):
        bf = C.build_from_drawing(d, C.CadOptions(to_scale=True, **kw))
        fz[k] = {"lines_used": bf.lines_used, "pipes": len(bf.net.pipes), "size_counts": bf.size_counts,
                 "sizes_in": sorted({round(v, 3) for v in bf.diam_in.values()}),
                 "notes": [x.split(";")[0] for x in bf.notes if "frozen" in x]}
    out["frozen_layers.dxf"] = fz
    # files cut short
    for n in [f"ky4_cut_{round(100 * x)}.dxf" for x in FX.CUT_SHORT_SHARES] + ["ky4_no_eof.dxf"]:
        data = _read_bytes(paths[n])
        bc = _cad(data, n)
        out[n] = {"bytes": len(data), "ezdxf_alone_lines": len(recover.read(io.BytesIO(data))[0].modelspace()),
                  "lines": len(bc.drawing.lines), "texts": len(bc.drawing.texts), "pipes": len(bc.net.pipes),
                  "warned_cut_short": any(w.startswith("The file looks cut short") for w in bc.warnings),
                  "refused": bc.refusal is not None}
    # AutoCAD's diameter sign, U+2205
    d = C.read_drawing(_read_bytes(paths["diameter_sign.dxf"]), "diameter_sign.dxf")
    bd = C.build_from_drawing(d, C.CadOptions(to_scale=True))
    out["diameter_sign.dxf"] = {"labels": [[t.text, C.label_size(t.text)["value"], C.label_size(t.text)["marked"]]
                                           for t in d.texts], "size_counts": bd.size_counts,
                                "sizes_in": sorted({round(v, 3) for v in bd.diam_in.values()})}
    # C from a label's material on a line whose layer names its size
    bc = _cad(_read_bytes(paths["c_from_label.dxf"]), "c_from_label.dxf")
    out["c_from_label.dxf"] = {"size_counts": bc.size_counts, "c_counts": bc.c_counts}
    # no size anywhere: asked, built only at the size the operator gives, and said so
    data = _read_bytes(paths["no_sizes.dxf"])
    ask, given = _cad(data, "no_sizes.dxf"), _cad(data, "no_sizes.dxf", default_diameter_in=8.0)
    out["no_sizes.dxf"] = {"questions": qs(ask), "built": ask.inp_text is not None,
                           "given": {"built": given.inp_text is not None, "default_size": given.default_size,
                                     "warned_size_you_gave": any("the size you gave" in w for w in given.warnings)}}
    # a text too long to be read is said
    d = C.read_drawing(_read_bytes(paths["long_label.dxf"]), "long_label.dxf")
    out["long_label.dxf"] = {"texts_too_long": d.counts.get("texts_too_long", 0),
                             "noted": any(x.startswith("1 text longer than") for x in d.notes)}
    # the readiness report's Diameters row: from the file only for layer names; a label is matched by distance
    from .importers import prepare_inp
    rows = {}
    for k, bb in (("layer names only", bc), ("labels only", bd), ("the size you gave", given)):
        with tempfile.TemporaryDirectory(prefix="rm_t18_ready_") as tmp:
            _, rep = prepare_inp(bb.inp_text.encode("utf-8"), tmp, "free_chlorine", name="skeleton.inp")
        row = next(r for r in C.readiness(bb, rep).rows if r[0] == "Diameters")
        rows[k] = {"value": row[1], "source": row[2]}
    out["diameters_row"] = rows
    return out


LAYER_CUE_NAMES = ("W-MAIN-CL-ZONE-2", "WATER-CL-PHASE-3", "C-WATR-PIPE-CL-2", "WATER-ST-12", "W-PE-SHEET-4",
                   "WATER-ZONE-2-PVC", "PVC-ZONE-3", "W-POLY-ZONE-3", "W-AC-DIST-10", "WS-ZONE-2", "W_8_CL",
                   "WATER-ZONE-2", "W_8_PVC", "PVC-8", "W-8PVC", "W-\u22058", "MAINS-8-DIA", "WATER-8IN")
NEAR_LIMIT_TOL_FT = 1.0                       # the near-limit drawing's rows are 2 ft apart: joined within 1 ft


def near_limit_build(path: str) -> dict:
    """The final verification of task 18, in a process of its own (cloud.run_in_folder), so its peak memory is its
    own: label_pairs_near_limit.dxf read, its pairs of a label and a segment counted (a build whose limit is one pair
    refuses with the count), then built with the local limits; the peak resident memory after the read and after the
    build (GB) and the build's seconds, which vary from run to run, beside what the build gave."""
    import resource
    import sys as _sys
    import time as _time
    from . import cad_read as C
    from .topology import FT
    _quiet()
    scale = 1.0 if _sys.platform == "darwin" else 1024.0        # ru_maxrss is in bytes on macOS, KiB on Linux
    peak = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 1e9
    d = C.read_drawing(_read_bytes(path), os.path.basename(path))
    after_read = peak()
    opts = C.CadOptions(to_scale=True, tolerance_m=NEAR_LIMIT_TOL_FT * FT)
    counted = C.build_from_drawing(d, opts, {"max_label_pairs": 1}).labels.get("pairs_within_reach")
    t = _time.perf_counter()
    b = C.build_from_drawing(d, opts)
    secs = _time.perf_counter() - t
    return {"lines": len(d.lines), "texts": len(d.texts), "pairs_within_reach": counted,
            "under_local_limit": counted is not None and counted <= C.MAX_LABEL_PAIRS, "refused": b.refusal is not None,
            "pipes": len(b.net.pipes) if b.net else None, "assigned": b.labels.get("assigned"),
            "size_counts": b.size_counts, "_measured": {"peak_gb_after_read": round(after_read, 3),
                                                        "peak_gb_after_build": round(peak(), 3),
                                                        "build_s": round(secs, 1)}}


def part_cad_final(paths: dict) -> dict:
    """dxf_final_verification: one case for each issue the final verification of task 18 found that the hostile part
    does not hold: which layer names cue their number as a size; a cued number its labels dispute (asked); a byte
    order mark before a drawing whose handles are drawn twice, or that is cut short; a binary DXF whose handles are
    drawn twice (not repaired: a recorded limit); the points a bulge between infinite vertices is charged; and the
    drawing just under the local limit of pairs of a label and a segment, built in a process of its own.  Its key
    '_measured' (peak memory and seconds, which vary from run to run) is not part of the record."""
    import io
    import logging
    import tempfile as _tf
    from ezdxf import recover
    from . import cad_read as C
    from . import cloud
    from . import fixtures as FX
    logging.getLogger("ezdxf").setLevel(logging.ERROR)
    out = {"layer_cues": {n: {k: C.layer_size(n)[k] for k in ("value", "unit", "cued")} for n in LAYER_CUE_NAMES},
           "cue_ambiguous": list(C.CUE_AMBIGUOUS), "dispute_lines": C.DISPUTE_LINES}
    # a number its name cues (W-2-PVC) that the labels on its lines dispute, beside one stray label (W-8-PVC)
    data = _read_bytes(paths["disputed_layer.dxf"])
    disp = {}
    for k, kw in (("asked", {}), ("not_sizes", {"layer_numbers": "none"}), ("inches", {"layer_numbers": "in"})):
        bd = _cad(data, "disputed_layer.dxf", to_scale=True, **kw)
        q = next((x["text"] for x in bd.questions if x["key"] == "layer_numbers"), "")
        disp[k] = {"questions": [[x["key"], x["blocking"]] for x in bd.questions], "built": bd.inp_text is not None,
                   "size_counts": bd.size_counts, "sizes_in": sorted({round(v, 3) for v in bd.diam_in.values()}),
                   "disputed": bd.labels.get("layers_disputed", []),
                   "question_names": [n for n in FX.DISPUTED_LAYERS if n in q],
                   "warned_stray_label": any("with a label of another size near" in w for w in bd.warnings)}
    out["disputed_layer.dxf"] = disp
    # a byte order mark before the first group code
    for n, plain in (("bom_dup_handles.dxf", "dup_handles.dxf"), ("bom_ky4_no_eof.dxf", "ky4_no_eof.dxf")):
        data = _read_bytes(paths[n])
        d = C.read_drawing(data, n)
        try:
            alone = len(recover.read(io.BytesIO(data))[0].modelspace())
        except Exception as e:  # noqa: BLE001  (ezdxf refuses the file)
            alone = type(e).__name__
        out[n] = {"starts_with_bom": data.startswith(C._BOM), "same_as": plain, "lines": len(d.lines),
                  "ezdxf_alone_lines": alone, "handles_renumbered": d.counts.get("handles_renumbered", 0),
                  "cut_short": bool(d.counts.get("cut_short")), "refused": d.refusal is not None,
                  "warned_cut_short": any(w.startswith("The file looks cut short") for w in d.warnings)}
    # a binary DXF whose handles are drawn twice: read as it is, the audit's removals warned (a recorded limit)
    d = C.read_drawing(_read_bytes(paths["binary_dup_handles.dxf"]), "binary_dup_handles.dxf")
    out["binary_dup_handles.dxf"] = {"kind": d.kind, "lines": len(d.lines), "of_lines": 10,
                                     "audit_removed": d.counts.get("audit_removed", 0),
                                     "handles_renumbered": d.counts.get("handles_renumbered", 0),
                                     "warned": any(w.startswith("The DXF repair removed") for w in d.warnings)}
    # a bulge between two vertices that compare equal but whose chord is not 0 (infinite ones): charged what it makes
    inf = float("inf")
    out["bulge_counts"] = [{"vertices": [[str(x), str(y), b] for x, y, b in xyb],
                            "counted": C._bulged_points(xyb, False), "made": len(C._poly_with_bulges(xyb, False))}
                           for xyb in ([(0.0, 0.0, 1.0), (10.0, 0.0, 0.0)], [(inf, inf, 1.0), (inf, inf, 0.0)],
                                       [(inf, 0.0, 1.0), (inf, 0.0, 0.0)], [(5.0, 5.0, 1.0), (5.0, 5.0, 0.0)])]
    n = FX.HOSTILE_INF_BULGES                   # hostile_inf_bulges.dxf's polyline: what flattening it would make
    out["inf_polyline"] = {"vertices": n, "points_its_arcs_make": C._bulged_points([(inf, inf, 1.0)] * n, False)}
    # just under the local limit of pairs: built in a process of its own, its peak memory measured
    from . import intake_exp as this      # so the child unpickles the function under -m too (not as __main__'s)
    with _tf.TemporaryDirectory(prefix="rm_t18_near_") as tmp:
        near = cloud.run_in_folder(tmp, this.near_limit_build, paths["label_pairs_near_limit.dxf"])
    out["label_pairs_near_limit.dxf"] = near
    out["_measured"] = {"label_pairs_near_limit.dxf": near.pop("_measured")}
    return out


RECORD_18_TIMING = os.path.join(OUT, "cad_mtext_timing.json")
MTEXT_TIMING_CHARS = (200_000, 400_000, 800_000)
MTEXT_SHAPES = {"codes_A": "\\A1;", "stack_S": "\\S1^2;", "height_H": "\\H1.5x;", "paragraphs": "\\P",
                "font_f": "\\fArial;"}
LOADER_SHAPES = ("text_spaces", "layer_long", "float_long", "many_lines")


def _timing_dxf(pairs: list) -> bytes:
    head = [(0, "SECTION"), (2, "HEADER"), (9, "$ACADVER"), (1, "AC1015"), (9, "$INSUNITS"), (70, 2), (0, "ENDSEC"),
            (0, "SECTION"), (2, "ENTITIES")]
    return ("\r\n".join(f"{c}\r\n{v}" for c, v in head + pairs + [(0, "ENDSEC"), (0, "EOF")]) + "\r\n").encode("latin-1")


def _timing_mtext(s: str) -> list:
    return [(0, "MTEXT"), (5, "1A"), (100, "AcDbEntity"), (8, "T"), (100, "AcDbMText"), (10, 0.0), (20, 0.0), (40, 2.0),
            (71, 1), (1, s)]


def _loader_case(shape: str, n: int) -> list:
    if shape == "text_spaces":
        return [(0, "TEXT"), (5, "1B"), (100, "AcDbEntity"), (8, "T"), (100, "AcDbText"), (10, 0.0), (20, 0.0),
                (40, 2.0), (1, " " * n + '8"')]
    if shape == "layer_long":
        return [(0, "TEXT"), (5, "1B"), (100, "AcDbEntity"), (8, "L" * n), (100, "AcDbText"), (10, 0.0), (20, 0.0),
                (40, 2.0), (1, "8")]
    if shape == "float_long":
        return [(0, "LINE"), (5, "1D"), (100, "AcDbEntity"), (8, "W"), (100, "AcDbLine"), (10, "1" * n), (20, 0.0),
                (11, 1.0), (21, 0.0)]
    return [p for i in range(n // 60) for p in [(0, "LINE"), (8, "W"), (10, float(i)), (20, 0.0), (11, float(i) + 1),
                                               (21, 0.0)]]


def cad_mtext_timing() -> dict:
    """The probe that found ezdxf's default MTEXT reader quadratic, rerun and kept (the review of task 18): seconds to
    load a one-MTEXT drawing in recover mode and read its text with the default reader (plain_text, its 'fast' mode)
    and with the accurate parser cad_read uses (plain_mtext), for runs of formatting codes of 200,000, 400,000 and
    800,000 characters, and the loader alone on other shapes.  Timed on this Mac, so it varies from run to run and no
    check reruns it; written by python -m residualmap.intake_exp cad_timing."""
    import io
    import logging
    import platform
    import time
    import ezdxf
    from ezdxf import recover
    from ezdxf.tools.text import plain_mtext
    from . import fixtures as FX
    logging.getLogger("ezdxf").setLevel(logging.ERROR)

    def timed(data: bytes, read) -> float:
        t = time.perf_counter()
        doc, _ = recover.read(io.BytesIO(data))
        for e in doc.modelspace():
            if e.dxftype() in ("TEXT", "MTEXT"):
                read(e)
        return round(time.perf_counter() - t, 3)

    def growth(row: dict) -> dict:
        ks = sorted(row, key=int)
        return {f"{b}_over_{a}": round(row[b] / max(row[a], 1e-6), 2) for a, b in zip(ks, ks[1:])}

    out = {"about": " ".join(cad_mtext_timing.__doc__.split("\n\n")[0].split()),
           "written": time.strftime("%Y-%m-%d %H:%M"), "ezdxf": ezdxf.__version__, "python": platform.python_version(),
           "machine": platform.machine(), "chars": list(MTEXT_TIMING_CHARS)}
    fast = {str(n): timed(_timing_dxf(_timing_mtext("\\A1;" * (n // 4))), lambda e: e.plain_text())
            for n in MTEXT_TIMING_CHARS}
    out["default_reader_codes_A_s"] = fast
    out["default_reader_growth"] = growth(fast)
    hi = MTEXT_TIMING_CHARS[-1]
    # an extrapolation, not a run: time grows with the square of the characters, so the hosted demo's largest MTEXT
    # (the hostile fixture's 4,800,000 characters) would take (4,800,000 / 800,000) squared times the 800,000 run
    out["extrapolated_default_reader_s_at_hostile_chars"] = {
        "chars": FX.HOSTILE_MTEXT_CHARS, "seconds": round(fast[str(hi)] * (FX.HOSTILE_MTEXT_CHARS / hi) ** 2),
        "minutes": round(fast[str(hi)] * (FX.HOSTILE_MTEXT_CHARS / hi) ** 2 / 60), "is_extrapolation": True}
    acc = {}
    for shape, code in MTEXT_SHAPES.items():
        acc[shape] = {str(n): timed(_timing_dxf(_timing_mtext(code * (n // len(code)))), lambda e: plain_mtext(e.text))
                      for n in MTEXT_TIMING_CHARS[1:]}
        acc[shape]["growth"] = growth({k: v for k, v in acc[shape].items()})
    out["accurate_parser_s"] = acc
    loader = {}
    for shape in LOADER_SHAPES:
        row = {str(n): timed(_timing_dxf(_loader_case(shape, n)), lambda e: e.plain_text()) for n in MTEXT_TIMING_CHARS[1:]}
        loader[shape] = {**row, "growth": growth(row)}
    out["loader_s"] = loader
    return out


def cad_record() -> dict:
    from . import fixtures as FX
    _quiet()
    with tempfile.TemporaryDirectory(prefix="rm_t18_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.CAD_FIXTURES)
        return {"about": "Task 18's CAD importer, converter self-consistency on public models drawn by our own writers "
                         "(residualmap/fixtures.py), never a test on a real drawing (R6): ky4's mains as an R12 drawing "
                         "written by hand and read by ezdxf against the pre-registered bar, and as an R2018 drawing "
                         "ezdxf wrote (part 'rebuild'), the tolerance table in feet on ky4 and ky10 (part "
                         "'tolerance'), labels and propagation, blocks and Civil 3D objects, layer styles, units and "
                         "MTEXT, the backup and binary copies, DWG, DGN and IFC guidance, hostile drawings (part "
                         "'hostile'), and since the review of task 18 one drawing for each case it found (part "
                         "'review') and since its final verification one for each issue that found (part 'final'). "
                         "Every CAD-built model is topology and sizes only, and needs the operator form. "
                         "Written by python -m residualmap.intake_exp cad; the task18 checks in residualmap/checks.py "
                         "rerun every part and compare.",
                "rebuild": part_cad_rebuild(paths), "tolerance": part_cad_tolerance(paths),
                "labels": part_cad_labels(paths), "blocks": part_cad_blocks(paths), "styles": part_cad_styles(paths),
                "units": part_cad_units(paths), "bak_binary": part_cad_bak_binary(paths),
                "guidance": part_cad_guidance(paths), "hostile": part_cad_hostile(paths),
                "review": part_cad_review(paths), "final": _unmeasured(part_cad_final(paths))}


# ----------------------------------------------------------------------------- task 19: the operator form's record
# operator_form_completes_skeleton's pre-registered bars (plan task 19, written before its scored run): the skeleton
# from ky4_clean_r12.dxf (joined within fixtures.FORM_TOL_FT, 1 ft: the build's own setting, fixed before the scored
# run; the plan named no joining distance, and task 17's record, gis_roundtrip.json's near_node_misses, already listed
# the three ky4 dead ends that 1 ft joins as tees) completed by ky4_form_from_original.csv against the
# original, at every junction of the original and every hour of the scored day (the last of a 7-day run): water age's
# median absolute difference at most 0.5 h and 95th percentile at most 2 h; chlorine (one grid member, kb 0.40 per day,
# kw 0.70 m per day, gamma 1, at ky4's truth dose of 2.0 mg/L) median at most 0.01 mg/L and 95th percentile at most
# 0.03 mg/L; the junctions whose daily minimum is below 0.2 mg/L the same set, Jaccard at least 0.95.  A drawing and the
# form carry no solver or quality options, so 'the original' is ky4 with the skeleton's options (FORM_ALIGN: ky4's own
# file has 100 trials, an accuracy of 0.0001 and a diffusivity of 0.05); the difference those options make on the
# original itself is recorded beside the bars, not barred.
RECORD_19 = os.path.join(OUT, "operator_form.json")
# The scored run itself, as its console output gave it (23 October 2026, run from the build's scratch folder outside
# the repository).  Part 'complete' is that setup rerun after changes that do not touch what it measures, so its
# numbers differ in the fifth significant digit; the verdict, part by part, is the same.  Copied, not rerun.  A first
# attempt stopped before any score: its comparison could not map the three junctions the skeleton had joined to a
# passing main (their nodes sit a fraction of a foot from ky4's), and it then took the nearest node within the joining
# distance, the run scored here.
SCORED_RUN_19 = {"what": "operator_form_completes_skeleton's pre-registered scored run, as its console output gave it",
                 "changes_since": ["the skeleton's zero-demand entries are cleared before the form's demands are added, "
                                   "so each junction carries the form's demand alone (the scored run had a zero entry "
                                   "beside it)",
                                   "only a pump the form adds needs a pump_rule row (ky4's two are both the form's)",
                                   "WNTR's header comments are left out of the written model",
                                   "task 19's review fixed cases ky4's form does not use (part 'review'): the completed "
                                   "model's text is the same byte for byte after the fixes as before them"],
                 "age_h": {"median": 1.471684, "p95": 22.577524, "max": 140.577175, "mean": 4.736041},
                 "chlorine_mgL": {"median": 0.001502, "p95": 0.092337, "max": 1.029479, "mean": 0.019533},
                 "lows": {"original": 897, "completed": 898, "both": 897, "jaccard": 0.998886},
                 "bar_parts_met": {"age_median_h": False, "age_p95_h": False, "chlorine_median_mgL": True,
                                   "chlorine_p95_mgL": False, "lows_jaccard": True}, "bar_met": False}
FORM_BAR = {"age_median_h": 0.5, "age_p95_h": 2.0, "chlorine_median_mgL": 0.01, "chlorine_p95_mgL": 0.03,
            "lows_jaccard_at_least": 0.95}
FORM_MEMBER = (0.40, 0.70, 1.0)        # kb 1/day, kw m/day, gamma
FORM_ALIGN = {"hydraulic": ("trials", "accuracy", "checkfreq", "maxcheck", "damplimit", "unbalanced",
                            "unbalanced_value", "viscosity", "specific_gravity", "emitter_exponent"),
              "quality": ("diffusivity", "tolerance")}
FORM_MOVE_H = 2.0                      # the PRVs-missing record: a junction 'moves' above this, as elevation_prv_ky10


def _form_runs(path: str, dose: float) -> tuple:
    """(water age, chlorine) on the scored day, hour x junction, of one file."""
    from .simulate import simulate_nominal_chlorine
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        cl = simulate_nominal_chlorine(path, *FORM_MEMBER, source_dose=dose, file_prefix=os.path.join(d, "c"))
    return _fid_age(path), cl


def _write_wn(wn, d: str, name: str) -> str:
    import wntr
    p = os.path.join(d, name)
    wntr.network.write_inpfile(wn, p, units="GPM")
    return p


def _diff(a, b) -> dict:
    x = np.abs(np.asarray(a, float) - np.asarray(b, float)).ravel()
    return {"median": float(np.median(x)), "p95": float(np.quantile(x, 0.95)), "max": float(x.max()),
            "mean": float(x.mean())}


def _lows(cl, cols) -> set:
    m = cl[cols].min()
    return set(m.index[m < FID_THRESHOLD])


def _jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if a | b else 1.0


def _form(paths: dict, name: str):
    from .operator_form import read_operator_form
    with open(paths[name], "rb") as fh:
        return read_operator_form(fh.read(), name)


def ky4_clean_skeleton(paths: dict):
    """The CAD skeleton of ky4_clean_r12.dxf, joined within fixtures.FORM_TOL_FT."""
    from . import cad_read as CAD
    from . import fixtures as FX
    from . import topology as T
    with open(paths["ky4_clean_r12.dxf"], "rb") as fh:
        d = CAD.read_drawing(fh.read(), "ky4_clean_r12.dxf")
    return CAD.build_from_drawing(d, CAD.CadOptions(tolerance_m=FX.FORM_TOL_FT * T.FT))


def part_form_complete(paths: dict) -> dict:
    """operator_form_completes_skeleton: the clean ky4 drawing's skeleton and the operator's form against the original
    (with the skeleton's options), on the pre-registered bars."""
    import wntr
    from scipy.spatial import cKDTree
    from .experiment import NET_TRUTH
    from .operator_form import apply_operator_form, skeleton_from_build
    from .simulate import LIB
    _quiet()
    dose = float(NET_TRUTH["ky4"]["source_dose"])
    b = ky4_clean_skeleton(paths)
    res = apply_operator_form(skeleton_from_build(b, "cad"), _form(paths, "ky4_form_from_original.csv"))
    assert res.ok, (res.refused_lines()[:5], res.blocking)
    orig = wntr.network.WaterNetworkModel(os.path.join(LIB, "ky4.inp"))
    comp = res.wn
    aligned = wntr.network.WaterNetworkModel(os.path.join(LIB, "ky4.inp"))
    for sec, keys in FORM_ALIGN.items():
        for k in keys:
            setattr(getattr(aligned.options, sec), k, getattr(getattr(comp.options, sec), k))
    # the original's junctions on the completed model: the nearest node within the joining distance (a node the
    # skeleton joined to a passing line sits where the ends met, a fraction of a foot away), and the two ends of each
    # pump, which share one place, by the pump
    from . import fixtures as FX
    names = comp.node_name_list
    tree = cKDTree(np.array([comp.get_node(n).coordinates for n in names], float))
    node_of = {}
    for j in orig.junction_name_list:
        dist, i = tree.query(orig.get_node(j).coordinates)
        if dist <= FX.FORM_TOL_FT:
            node_of[j] = names[int(i)]
    for pn in orig.pump_name_list:
        p, q = orig.get_link(pn), comp.get_link(pn)
        node_of[p.start_node_name], node_of[p.end_node_name] = q.start_node_name, q.end_node_name
    js = orig.junction_name_list
    assert all(j in node_of for j in js), [j for j in js if j not in node_of][:5]
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        a_age, a_cl = _form_runs(_write_wn(aligned, d, "aligned.inp"), dose)
        o_age, o_cl = _form_runs(os.path.join(LIB, "ky4.inp"), dose)
        c_path = os.path.join(d, "completed.inp")
        with open(c_path, "w", encoding="utf-8") as fh:
            fh.write(res.inp_text)
        c_age, c_cl = _form_runs(c_path, dose)
    cols = [node_of[j] for j in js]
    age, cl = _diff(c_age[cols].values, a_age[js].values), _diff(c_cl[cols].values, a_cl[js].values)
    lows_a = _lows(a_cl, js)
    lows_c = {j for j in js if node_of[j] in _lows(c_cl, cols)}
    jac = _jaccard(lows_a, lows_c)
    bar = {"age_median_h": age["median"] <= FORM_BAR["age_median_h"], "age_p95_h": age["p95"] <= FORM_BAR["age_p95_h"],
           "chlorine_median_mgL": cl["median"] <= FORM_BAR["chlorine_median_mgL"],
           "chlorine_p95_mgL": cl["p95"] <= FORM_BAR["chlorine_p95_mgL"],
           "lows_jaccard": jac >= FORM_BAR["lows_jaccard_at_least"]}
    r = b.net.report
    import hashlib
    return {"completed_sha1": hashlib.sha1(res.inp_text.encode("utf-8")).hexdigest(),
            "skeleton": {"tolerance_ft": FX.FORM_TOL_FT, "pipes": r["pipes"], "nodes": r["nodes"], "components": r["components"],
                         "tees_split": r["tees_split"], "original_pipes": len(orig.pipe_name_list),
                         "original_nodes": len(orig.node_name_list), "c_classes": b.c_counts["classes"]},
            "form": {"rows": res.counts["rows"], **{k: v for k, v in res.counts.items() if k != "rows"}},
            "completed": {"junctions": len(comp.junction_name_list), "tanks": len(comp.tank_name_list),
                          "reservoirs": len(comp.reservoir_name_list), "pumps": len(comp.pump_name_list),
                          "pipes": len(comp.pipe_name_list), "dose_mgL": res.dose_mgL},
            "member": {"kb_per_day": FORM_MEMBER[0], "kw_m_per_day": FORM_MEMBER[1], "gamma": FORM_MEMBER[2],
                       "dose_mgL": dose},
            "age_h": age, "chlorine_mgL": cl,
            "lows": {"original": len(lows_a), "completed": len(lows_c), "both": len(lows_a & lows_c), "jaccard": jac},
            "bar": FORM_BAR, "bar_parts_met": bar, "bar_met": all(bar.values()),
            "original_options_effect": {"age_h": _diff(o_age[js].values, a_age[js].values),
                                        "chlorine_mgL": _diff(o_cl[js].values, a_cl[js].values),
                                        "aligned": {s: list(k) for s, k in FORM_ALIGN.items()},
                                        "ky4_own": {"trials": orig.options.hydraulic.trials,
                                                    "accuracy": orig.options.hydraulic.accuracy,
                                                    "diffusivity": orig.options.quality.diffusivity},
                                        "skeleton_own": {"trials": comp.options.hydraulic.trials,
                                                         "accuracy": comp.options.hydraulic.accuracy,
                                                         "diffusivity": comp.options.quality.diffusivity}}}


AFTER_STOP_TOL_FT = 0.1                # the after-the-stop characterization's joining distance, chosen after the scored
                                       # run, just under the 0.104 ft at which the nearest of the dead ends the 1 ft join
                                       # met passes its main: not a rescoring
PRODUCT_TOL_FT = 6.0                   # the CAD step's default joining distance, for the count of such dead ends there


def _near_touches(wn, within_ft: float) -> list:
    """Junctions of a model that end within within_ft of a pipe they are not on (a dead end drawn touching a passing
    main): what a drawing of ky4 joins as a tee though ky4 does not."""
    out = []
    geo_segs = []
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        pts = [p.start_node.coordinates] + list(p.vertices) + [p.end_node.coordinates]
        for a, b in zip(pts, pts[1:]):
            geo_segs.append((pn, a, b))
    S = np.array([(a[0], a[1], b[0], b[1]) for _, a, b in geo_segs], float)
    owner = [pn for pn, _, _ in geo_segs]
    A, B = S[:, :2], S[:, 2:]
    AB = B - A
    L2 = np.where((AB ** 2).sum(1) > 0, (AB ** 2).sum(1), 1)
    for j in wn.junction_name_list:
        P = np.array(wn.get_node(j).coordinates, float)
        own = set(wn.get_links_for_node(j))
        t = np.clip(((P - A) * AB).sum(1) / L2, 0, 1)
        d = np.hypot(*(A + AB * t[:, None] - P).T)
        for i in np.argsort(d)[:8]:
            if owner[i] not in own and d[i] <= within_ft and 0 < t[i] < 1:
                out.append({"junction": j, "pipe": owner[i], "distance_ft": float(d[i]),
                            "junction_pipes": len(own)})
                break
    return sorted(out, key=lambda x: x["junction"])


def part_form_after_the_stop(paths: dict) -> dict:
    """After operator_form_completes_skeleton missed its bar (a characterization, not a rescoring): ky4's dead ends
    that the clean drawing puts within the 1 ft join of a passing main, which the skeleton joins as tees though ky4
    does not; what the drawing's geometric lengths alone do to the original; and the same form on the same drawing
    joined within AFTER_STOP_TOL_FT, where none of those tees is made, scored the same way."""
    import math as _m
    import wntr
    from . import cad_read as CAD
    from . import fixtures as FX
    from . import topology as T
    from .simulate import LIB
    _quiet()
    orig = wntr.network.WaterNetworkModel(os.path.join(LIB, "ky4.inp"))
    touches = _near_touches(orig, FX.FORM_TOL_FT)
    comp_opts = None
    b = ky4_clean_skeleton(paths)
    from .operator_form import apply_operator_form, skeleton_from_build
    res = apply_operator_form(skeleton_from_build(b, "cad"), _form(paths, "ky4_form_from_original.csv"))
    comp_opts = res.wn.options

    def aligned():
        wn = wntr.network.WaterNetworkModel(os.path.join(LIB, "ky4.inp"))
        for sec, keys in FORM_ALIGN.items():
            for k in keys:
                setattr(getattr(wn.options, sec), k, getattr(getattr(comp_opts, sec), k))
        return wn
    geo = aligned()
    for pn in geo.pipe_name_list:
        p = geo.get_link(pn)
        pts = [p.start_node.coordinates] + list(p.vertices) + [p.end_node.coordinates]
        p.length = sum(_m.hypot(q[0] - o[0], q[1] - o[1]) for o, q in zip(pts, pts[1:])) * T.FT
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        a_age, a_cl = _form_runs(_write_wn(aligned(), d, "aligned.inp"), 2.0)
        g_age, g_cl = _form_runs(_write_wn(geo, d, "geometric.inp"), 2.0)
    js = orig.junction_name_list
    lengths_only = {"age_h": _diff(g_age[js].values, a_age[js].values),
                    "chlorine_mgL": _diff(g_cl[js].values, a_cl[js].values)}
    keep = FX.FORM_TOL_FT
    try:
        FX.FORM_TOL_FT = AFTER_STOP_TOL_FT
        tight = part_form_complete(paths)
    finally:
        FX.FORM_TOL_FT = keep
    return {"counts_as": "a characterization after the stop, not a rescoring: the scored run (part 'complete') "
                         "missed its bar and the miss stands",
            "dead_ends_touching_a_main_within_1_ft": touches,
            "junctions_touching_a_main_within_6_ft": len(_near_touches(orig, PRODUCT_TOL_FT)),
            "geometric_lengths_only": lengths_only,
            "joined_within_0p1_ft": {"tolerance_ft": AFTER_STOP_TOL_FT, "skeleton": tight["skeleton"],
                                     "age_h": tight["age_h"], "chlorine_mgL": tight["chlorine_mgL"],
                                     "lows": tight["lows"],
                                     # the bar's values applied after the stop, at a joining distance chosen after
                                     # the scored run (the review: a key 'bar_met' read as a pass of the bar)
                                     "within_the_bar_values_post_hoc": tight["bar_parts_met"]}}


def _same_model(a_path: str, b_wn) -> dict:
    """Whether a model written by WNTR from b_wn is the file a_path's model writes, and the runs' largest differences."""
    import wntr
    from .operator_form import write_text
    a_wn = wntr.network.WaterNetworkModel(a_path)
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        bp = _write_wn(b_wn, d, "b.inp")
        a_age, a_cl = _form_runs(a_path, 1.2)
        b_age, b_cl = _form_runs(bp, 1.2)
    return {"same_inp_text": write_text(a_wn) == write_text(b_wn),
            "age_max_abs_h": float(np.abs(a_age.values - b_age.values).max()),
            "chlorine_max_abs_mgL": float(np.abs(a_cl.values - b_cl.values).max())}


def _moved(a_path: str, b_path: str) -> dict:
    a_age, a_cl = _form_runs(a_path, 1.2)
    b_age, b_cl = _form_runs(b_path, 1.2)
    d = (a_age - b_age).abs()
    return {"age_h": _diff(a_age.values, b_age.values), "chlorine_mgL": _diff(a_cl.values, b_cl.values),
            "junctions_age_moved": int((d.max() > FORM_MOVE_H).sum()), "move_h": FORM_MOVE_H,
            "junctions": len(d.columns)}


def part_form_prvs(paths: dict) -> dict:
    """operator_form_restores_prvs_ky10: ky10_prv_removed.inp (each PRV an open pipe of 1 ft) and ky10_prv_form.csv
    against the order-fixed original: the same .inp text and 0.0 difference; and what the PRVs' loss moves."""
    from .operator_form import Skeleton, apply_operator_form
    _quiet()
    with open(paths["ky10_prv_removed.inp"], encoding="utf-8") as fh:
        sk = Skeleton(inp_text=fh.read(), name="ky10_prv_removed.inp")
    res = apply_operator_form(sk, _form(paths, "ky10_prv_form.csv"))
    assert res.ok, (res.refused_lines()[:5], res.blocking)
    return {"form_rows": res.counts["rows"], "prvs_restored": res.counts["prv"],
            "restored_vs_original": _same_model(paths["ky10_order_fixed.inp"], res.wn),
            "prvs_missing_vs_original": _moved(paths["ky10_order_fixed.inp"], paths["ky10_prv_removed.inp"])}


def part_form_boundary(paths: dict) -> dict:
    """operator_form_boundary_valve: the order-fixed ky10 (the pinned pipe open, as a skeleton has it) and
    ky10_closed_valve_form.csv (the pipe given by a point on it) against the closed original."""
    from . import fixtures as FX
    from .operator_form import Skeleton, apply_operator_form
    _quiet()
    with open(paths["ky10_order_fixed.inp"], encoding="utf-8") as fh:
        sk = Skeleton(inp_text=fh.read(), name="ky10_boundary_valve_open.inp")
    res = apply_operator_form(sk, _form(paths, "ky10_closed_valve_form.csv"))
    assert res.ok, (res.refused_lines()[:5], res.blocking)
    placed = list(res.placed.values())
    return {"pipe": FX.BOUNDARY_PIPE_KY10, "placed_on": placed[0] if placed else None,
            "closed_vs_closed_original": _same_model(paths["ky10_boundary_valve_closed.inp"], res.wn),
            "open_vs_closed_original": _moved(paths["ky10_order_fixed.inp"], paths["ky10_boundary_valve_closed.inp"])}


def part_form_no_source(paths: dict) -> dict:
    """no_source_never_runs: ky4's clean skeleton stops at the main page's intake (no entry point), and with the form
    less its source row the form stops with NO_SOURCE and writes no model."""
    from .importers import prepare_inp
    from .operator_form import NO_SOURCE, apply_operator_form, skeleton_from_build
    _quiet()
    b = ky4_clean_skeleton(paths)
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        path, rep = prepare_inp(b.inp_text.encode("utf-8"), d, name="ky4_clean_skeleton.inp")
        left = os.listdir(d)
    res = apply_operator_form(skeleton_from_build(b, "cad"), _form(paths, "ky4_form_no_source.csv"))
    return {"skeleton_prepared": path is not None, "skeleton_blocking": rep.blocking, "files_left": len(left),
            "form_ok": res.ok, "form_blocking": res.blocking, "no_source_sentence": NO_SOURCE in res.blocking,
            "model_written": res.inp_text is not None, "form_refused": len(res.refused)}


def part_form_bad_rows(paths: dict) -> dict:
    """form_parser_refuses_bad_rows: every row of operator_form_bad_rows.csv that cannot be used, with its reason, on
    ky4's clean skeleton, and no model."""
    from .operator_form import apply_operator_form, skeleton_from_build
    _quiet()
    b = ky4_clean_skeleton(paths)
    form = _form(paths, "operator_form_bad_rows.csv")
    res = apply_operator_form(skeleton_from_build(b, "cad"), form)
    return {"rows": len(form.rows) + len(form.refused), "refused_when_read": len(form.refused),
            "refused": res.refused, "model_written": res.inp_text is not None, "blocking": res.blocking}


# ----------------------------------------------------------------------------- task 19's review: the form's cases
# Each case the review of task 19 found, on the review's own drawings and layer (fixtures.FORM_REVIEW_*), with what
# the fixed form does now (part 'review' of operator_form.json, operator_form_review_cases).  None of them is in ky4's
# form, so the scored run's model is the same byte for byte (part 'complete', completed_sha1).
REVIEW_SRC = {"kind": "source", "id": "WELL", "location": "0,-800", "head_ft": "200", "dose_mgL": "1.2"}
REVIEW_USE = {"kind": "production", "gallons_per_day": "500000"}
REVIEW_TANK = {"kind": "tank", "id": "TANK", "location": "4000,2000", "base_elev_ft": "120", "diameter_ft": "40",
               "min_level_ft": "5", "max_level_ft": "30", "init_level_ft": "20"}
REVIEW_PUMP = {"kind": "pump", "location": "0,0", "pipe": "0,-400", "power_hp": "30"}
REVIEW_PRV = {"kind": "prv", "pipe": "1500,0", "upstream": "1000,0", "setting_psi": "40", "elevation_ft": "10"}
REVIEW_PATTERN = tuple(round(0.5 + h / 23.0, 4) for h in range(24))   # 24 distinct multipliers, so a shift shows
REVIEW_DEG = (-84.5, 38.05, 288_000.0, 364_000.0)    # a completed model put in degrees: lon0, lat0, ft per degree


def _review_build(paths: dict, name: str):
    """The skeleton of one of the review's drawings, joined within fixtures.FORM_TOL_FT, its layer names' bare numbers
    in inches (W_8_AC: AC cues no size, being also a common CAD abbreviation, so the step asks)."""
    from . import cad_read as CAD
    from . import fixtures as FX
    from . import topology as T
    with open(paths[name], "rb") as fh:
        d = CAD.read_drawing(fh.read(), name)
    b = CAD.build_from_drawing(d, CAD.CadOptions(tolerance_m=FX.FORM_TOL_FT * T.FT, layer_numbers="in"))
    assert b.inp_text, [q["key"] for q in b.questions]
    return b


def _review_form(rows):
    from .operator_form import form_csv, read_operator_form
    return read_operator_form(form_csv(rows), "form.csv")


def _review_run(text: str, hours: int):
    """EPANET hydraulics of a completed model's text for hours, reported hourly, in a temporary folder."""
    import wntr
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        p = os.path.join(d, "m.inp")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        wn = wntr.network.WaterNetworkModel(p)
        wn.options.time.duration = hours * 3600
        wn.options.time.report_timestep = 3600
        wn.options.quality.parameter = "NONE"
        return wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "r"))


def _node_at(wn, xy) -> str:
    return min(wn.node_name_list, key=lambda n: math.dist(wn.get_node(n).coordinates, xy))


def _c_by_material(res, sk) -> dict:
    out = {}
    for pn in sorted(res.wn.pipe_name_list):
        out.setdefault((sk.material or {}).get(pn) or "none", set()).add(float(res.wn.get_link(pn).roughness))
    return {m: sorted(v) for m, v in sorted(out.items())}


def _first_sentence(blocking: list) -> list:
    """Blocking sentences without the measured change a no-rule sentence quotes (read from the ladder's files)."""
    return [b.split(" Measured change")[0] for b in blocking]


def part_form_review(paths: dict) -> dict:
    """Task 19's review (operator_form_review_cases): each case it found, and what the form does now."""
    from collections import Counter
    from . import fixtures as FX
    from . import gis_read as G
    from . import topology as T
    from .operator_form import (GPM_M3S, PATTERN_NAME, FormRefused, Skeleton, _load, apply_operator_form,
                                read_operator_form, skeleton_from_build, write_text)
    _quiet()
    out = {}
    b = _review_build(paths, "form_review_r12.dxf")
    sk = skeleton_from_build(b, "cad")
    base = [REVIEW_SRC, REVIEW_USE]

    def apply(rows, s=sk):
        return apply_operator_form(s, _review_form(rows))

    def gpm(series) -> float:
        return float(series.mean()) / GPM_M3S

    # 1. pipe_c rows by material: DIP and AC share the C table's class, and each keeps its own C
    pc = lambda m, c: {"kind": "pipe_c", "material": m, "hw_c": c}
    rows_c = base + [pc("DIP", "140"), pc("AC", "100"), pc("PVC", "150"), pc("CI", "90")]
    rc = apply(rows_c)
    pvc_plastic = apply(base + [pc("PVC", "150"), pc("plastic", "120")])
    out["pipe_c_by_material"] = {
        "materials": dict(sorted(Counter(v or "none" for v in sk.material.values()).items())),
        "c_by_material": _c_by_material(rc, sk), "pipes": len(rc.wn.pipe_name_list),
        "readiness": next(v for i, v, _ in rc.rows if i == "Hazen-Williams C"),
        "second_row_for_one_material": apply(rows_c + [pc("ductile iron", "130")]).refused,
        "plastic_alone": _c_by_material(apply(base + [pc("plastic", "120")]), sk),
        "pvc_and_plastic": {"c_by_material": _c_by_material(pvc_plastic, sk),
                            "warnings": [w for w in pvc_plastic.warnings if "(pipe_c)" in w]},
        "cast_iron_then_all": _c_by_material(apply(base + [pc("CI", "90"), pc("all", "110")]), sk)}

    # 2. a layer's own C field (GIS, in longitude and latitude): kept by a pipe_c row, and by the age note
    with open(paths["form_review_mains.geojson"], "rb") as fh:
        layers, kind = G.read_layers(fh.read(), "form_review_mains.geojson")
    gb = G.build_from_layer(layers[0], G.GisOptions(tolerance_m=FX.FORM_TOL_FT * T.FT), "form_review_mains.geojson",
                            kind, layers)
    gsk = skeleton_from_build(gb, "gis")
    lon, lat = FX.form_review_lonlat((0, 0))
    gsrc = {"kind": "source", "id": "PLANT", "location": f"{lat:.7f},{lon:.7f}", "head_ft": "150", "dose_mgL": "1.2"}
    guse = {"kind": "production", "gallons_per_day": "300000"}
    g1 = apply([gsrc, guse, pc("cast iron", "100")], gsk)
    g0 = apply([gsrc, guse], gsk)
    g_all = apply([gsrc, guse, pc("all", "110")], gsk)     # the final verification: an 'all' row reaches M1 and M2 too
    out["layer_c_field_kept"] = {
        "c_from_field": sorted(gsk.c_kept), "c_before": {pn: gb.c_of[pn] for pn in sorted(gb.c_of)},
        "c_after_a_cast_iron_row": {pn: g1.wn.get_link(pn).roughness for pn in sorted(g1.wn.pipe_name_list)},
        "c_after_an_all_row": {pn: g_all.wn.get_link(pn).roughness for pn in sorted(g_all.wn.pipe_name_list)},
        "readiness_after_an_all_row": next(v for i, v, _ in g_all.rows if i == "Hazen-Williams C"),
        "readiness": next(v for i, v, _ in g1.rows if i == "Hazen-Williams C"),
        "age_note_without_rows": [w for w in g0.warnings if "Hazen-Williams" in w],
        "longitude_first": apply([dict(gsrc, location=f"{lon:.7f},{lat:.7f}"), guse], gsk).refused}

    # 3. two pumps side by side at one station, the well's main drawn from the station toward the well
    rp = apply(base + [REVIEW_TANK, dict(REVIEW_PUMP, id="PUMP-A"), dict(REVIEW_PUMP, id="PUMP-B"),
                       {"kind": "pump_rule", "pump": "PUMP-A", "on_hours": "always"},
                       {"kind": "pump_rule", "pump": "PUMP-B", "tank": "TANK", "on_below_ft": "10",
                        "off_above_ft": "25"}])
    fl = _review_run(rp.inp_text, 48).link["flowrate"]
    out["pumps_side_by_side"] = {
        "ok": rp.ok, "pumps": {p: [rp.wn.get_link(p).start_node_name, rp.wn.get_link(p).end_node_name]
                               for p in sorted(rp.wn.pump_name_list)},
        "links_at_the_well": sorted(rp.wn.get_links_for_node("WELL")),
        "mean_flow_gpm_48_h": {p: gpm(fl[p]) for p in sorted(rp.wn.pump_name_list)}}

    # 4. a PSV in place of a pipe: its elevation at the node its setting holds, upstream
    rv = apply(base + [dict(REVIEW_TANK, base_elev_ft="100", min_level_ft="2"),
                       {"kind": "psv", "id": "V1", "pipe": "3000,1500", "upstream": "3000,1000", "setting_psi": "40",
                        "elevation_ft": "150"}])
    up, dn = _node_at(rv.wn, (3000, 1000)), _node_at(rv.wn, (3000, 2000))
    run = _review_run(rv.inp_text, 24)
    lvl = run.node["pressure"]["TANK"] / FT_M
    out["psv_elevation_upstream"] = {
        "ok": rv.ok, "valve": [rv.wn.get_link("V1").start_node_name, rv.wn.get_link("V1").end_node_name],
        "upstream_elevation_ft": rv.wn.get_node(up).elevation / FT_M, "downstream_elevation_ft":
            rv.wn.get_node(dn).elevation / FT_M, "psv_mean_flow_gpm_24_h": gpm(run.link["flowrate"]["V1"]),
        "tank_level_ft_first_last": [float(lvl.iloc[0]), float(lvl.iloc[-1])]}

    # 5. a drawing in two parts, one source
    bg = _review_build(paths, "form_review_gap_r12.dxf")
    ru = apply(base, skeleton_from_build(bg, "cad"))
    # the final verification: a tank alone is no source in a part either (a tank in the far part, no source there)
    rt = apply(base + [{"kind": "tank", "id": "T2", "location": "5000,3000", "base_elev_ft": "120", "diameter_ft": "40",
                        "min_level_ft": "5", "max_level_ft": "30", "init_level_ft": "20"}], skeleton_from_build(bg, "cad"))
    out["part_without_source"] = {"components": bg.net.report["components"], "ok": ru.ok, "blocking": ru.blocking,
                                  "model_written": ru.inp_text is not None,
                                  "with_a_tank_there": {"ok": rt.ok, "blocking": _first_sentence(rt.blocking),
                                                        "model_written": rt.inp_text is not None}}

    # 6. an EPANET model whose clock starts at 08:00: a '6-22' pump starts on, and the pattern keeps each clock hour's
    # multiplier
    text8 = re.sub(r"(Start ClockTime\s+)12 am", r"\g<1>8:00 AM", sk.inp_text)
    assert text8 != sk.inp_text
    r8 = apply(base + [REVIEW_TANK, dict(REVIEW_PUMP, id="PUMP-A"),
                       {"kind": "pump_rule", "pump": "PUMP-A", "on_hours": "6-22"},
                       {"kind": "pattern", "multipliers": " ".join(f"{x:g}" for x in REVIEW_PATTERN)}],
               Skeleton(inp_text=text8, name="zone_0800.inp"))
    st = _review_run(r8.inp_text, 30).link["status"]["PUMP-A"]
    out["clock_starts_at_0800"] = {
        "ok": r8.ok, "start_clocktime_h": r8.wn.options.time.start_clocktime / 3600.0,
        "pump_starts": str(r8.wn.get_link("PUMP-A").initial_status),
        "pattern_first_multiplier": r8.wn.get_pattern(PATTERN_NAME).multipliers[0],
        "form_multiplier_at_08": REVIEW_PATTERN[8], "notes": r8.notes,
        "pump_open_by_model_hour": [int(x) for x in st.values]}

    # 7. a form applied to a model the form completed
    rr = apply([{"kind": "closed_valve", "pipe": "1500,0"}], Skeleton(inp_text=rc.inp_text, name="completed.inp"))
    out["applied_again"] = {"ok": rr.ok, "closed": rr.placed, "blocking": rr.blocking}

    # 8. values every earlier check passed that EPANET refuses or that mean nothing
    rule = {"kind": "pump_rule", "pump": "PA", "on_hours": "always"}
    head = dict(REVIEW_PUMP, id="PA", power_hp="")
    values = {"power_hp 0": [dict(REVIEW_PUMP, id="PA", power_hp="0"), rule],
              "three curve points with heads rising": [dict(head, flow_gpm="0 500 1000", head_ft="100 150 200"), rule],
              "a design point at 0 gpm": [dict(head, flow_gpm="0", head_ft="100"), rule],
              "a design point with a negative head": [dict(head, flow_gpm="500", head_ft="-100"), rule],
              "a PRV of 0 in": [dict(REVIEW_PRV, id="V1", diameter_in="0")],
              "a PRV set at -20 psi": [dict(REVIEW_PRV, id="V1", setting_psi="-20")],
              "a tank level below its bottom": [dict(REVIEW_TANK, min_level_ft="-10", init_level_ft="-5")]}
    got = {k: apply(base + v).refused for k, v in values.items()}
    for k, dose in (("a dose of -2 mg/L", "-2"), ("a dose written 1,500", "1,500")):
        got[k] = apply([dict(REVIEW_SRC, dose_mgL=dose), REVIEW_USE]).refused
    out["values_refused"] = got

    # 9. what no row check sees, EPANET's first hour does: two PRVs side by side at one node
    side = {"kind": "prv", "location": "2000,0", "pipe": "1500,0", "setting_psi": "40", "elevation_ft": "10"}
    r9 = apply(base + [dict(side, id="V1"), dict(side, id="V2")])
    out["first_hour_refuses"] = {"ok": r9.ok, "blocking": r9.blocking, "model_written": r9.inp_text is not None}

    # 10. a header with its unit in brackets, and a column the form does not read
    f10 = read_operator_form(b'kind,id,location,head_ft,Dose (mg/L),gallons_per_day,Shift\n'
                             b'source,WELL,"0,-800",200,1.5,,night\nproduction,,,,,500000,\n', "form.csv")
    r10 = apply_operator_form(sk, f10)
    out["headers"] = {"ok": r10.ok, "dose_mgL": r10.dose_mgL, "notes": r10.notes}

    # 11. an unnamed pump's name taken by another pump, and an unnamed pump with no rule
    r11 = apply(base + [REVIEW_TANK, dict(REVIEW_PUMP),
                        {"kind": "pump", "id": "PUMP1", "location": "1000,0", "pipe": "500,0", "power_hp": "30"},
                        {"kind": "pump_rule", "pump": "PUMP1", "tank": "TANK", "on_below_ft": "10",
                         "off_above_ft": "25"}])
    r11b = apply(base + [REVIEW_TANK, dict(REVIEW_PUMP)])
    out["unnamed_pump"] = {"its_name_taken": r11.refused, "no_rule": _first_sentence(r11b.blocking),
                           "no_rule_quotes_ky4_and_ky10": any("Measured change when pump rules are lost" in b_
                                                              and "ky4" in b_ and "ky10" in b_ for b_ in r11b.blocking)}
    # a PRV or PSV with no elevation is refused when the form is read (its setting means nothing without one)
    out["prv_without_elevation"] = _review_form(base + [{k: v for k, v in dict(REVIEW_PRV, id="V1").items()
                                                         if k != "elevation_ft"}]).refused

    # 12. rows that cannot all hold at once, and a well with its pump entered the way the form says
    conflicts = {"two PRVs in place of one pipe": [dict(REVIEW_PRV, id="V1"), dict(REVIEW_PRV, id="V2")],
                 "a PRV and a closed valve on one pipe": [dict(REVIEW_PRV, id="V1"),
                                                          {"kind": "closed_valve", "pipe": "1500,0"}],
                 "a pump splitting the well's node": [{"kind": "pump", "id": "WP", "location": "0,-800",
                                                       "pipe": "0,-400", "power_hp": "40"},
                                                      {"kind": "pump_rule", "pump": "WP", "on_hours": "always"}],
                 "a PRV splitting the tank's node": [REVIEW_TANK, {"kind": "prv", "id": "V1", "location": "4000,2000",
                                                                   "pipe": "3500,2000", "setting_psi": "40",
                                                                   "elevation_ft": "10"}],
                 "a PRV in place of the tank's pipe": [REVIEW_TANK, {"kind": "prv", "id": "V1", "pipe": "3500,2000",
                                                                     "upstream": "3000,2000", "setting_psi": "40",
                                                                     "elevation_ft": "10"}]}
    out["conflicts_refused"] = {k: apply(base + v).refused for k, v in conflicts.items()}
    out["conflicts_refused"]["a pump splitting the model's own reservoir"] = apply(
        [{"kind": "pump", "id": "WP", "location": "WELL", "pipe": "0,-400", "power_hp": "40"},
         {"kind": "pump_rule", "pump": "WP", "on_hours": "always"}],
        Skeleton(inp_text=rc.inp_text, name="completed.inp")).refused
    rw = apply(base + [{"kind": "pump", "id": "WP", "pipe": "0,-400", "upstream": "0,-800", "power_hp": "40"},
                       {"kind": "pump_rule", "pump": "WP", "on_hours": "always"}])
    out["well_and_its_pump"] = {"ok": rw.ok, "pump": [rw.wn.get_link("WP").start_node_name,
                                                      rw.wn.get_link("WP").end_node_name]}

    # 13. points: on a model whose coordinates are degrees, at a tee, beside two pipes
    lon0, lat0, fx, fy = REVIEW_DEG
    wdeg = _load(rc.inp_text)
    for n in wdeg.node_name_list:
        x, y = wdeg.get_node(n).coordinates
        wdeg.get_node(n).coordinates = (lon0 + x / fx, lat0 + y / fy)
    sdeg = Skeleton(inp_text=write_text(wdeg), name="degrees.inp")
    mid = (lon0 + 1500 / fx, lat0)
    cv = lambda pt, s=sk: apply([{"kind": "closed_valve", "pipe": pt}], s) if s is sdeg else \
        apply(base + [{"kind": "closed_valve", "pipe": pt}], s)
    far_deg, on_deg = cv(f"{mid[0] - 1:.7f},{mid[1]:.7f}", sdeg), cv(f"{mid[0]:.7f},{mid[1]:.7f}", sdeg)
    out["points"] = {"degrees_a_degree_off": far_deg.refused, "degrees_on_the_pipe": on_deg.placed,
                     "at_a_tee": cv("1000,0").refused, "beside_two_pipes": cv("999,1").refused,
                     "halfway_along_a_pipe": cv("1500,0").placed}

    # 14. the model's text: an ID given twice (task 15's scan), and a byte-order mark
    lines = sk.inp_text.split("\n")
    i = lines.index("[JUNCTIONS]") + 1
    while lines[i].startswith(";") or not lines[i].strip():
        i += 1
    dup = "\n".join(lines[:i + 1] + [lines[i]] + lines[i + 1:])
    rd = apply(base, Skeleton(inp_text=dup, name="duplicate.inp"))
    rb = apply(base, Skeleton(inp_text="﻿" + sk.inp_text, name="bom.inp"))
    out["model_text"] = {"id_given_twice": {"ok": rd.ok, "blocking": rd.blocking},
                         "byte_order_mark": {"ok": rb.ok, "same_model": rb.inp_text == apply(base, Skeleton(
                             inp_text=sk.inp_text, name="plain.inp")).inp_text}}

    # 15. the form's own words: a PDF given as the form, an elevation row with no elevation
    try:
        read_operator_form(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n1 0 obj\n<<>>\nendobj\n", "form.pdf")
        pdf = None
    except FormRefused as e:
        pdf = str(e)
    out["form_words"] = {"pdf": pdf, "elevation_row": _review_form(base + [{"kind": "elevation",
                                                                             "location": "1000,0"}]).refused}
    return out


def _outcome(r) -> dict:
    """What a form gave, as the final verification's record keeps it."""
    return {"ok": r.ok, "refused": r.refused, "blocking": _first_sentence(r.blocking),
            "model_written": r.inp_text is not None}


FINAL_LONG_HEADER = 6400               # the final verification's form: 'Dose (', this many spaces, 'x' (6,477 bytes)
# The final verifier's grammar fuzzer (seed 11: 3,000 random well-formed forms, on the review's drawing and its gap
# drawing, each completed model run for 24 h in EPANET and its report searched for disconnected nodes), as its console
# gave it before this pass's fixes, after the fix that counts junctions reachable from a source, and after the tank no
# source can fill was refused too.  Copied from the run logs (a scratch script, not part of the repository), not rerun.
FINAL_FUZZ_19 = {"what": "the final verifier's grammar fuzzer, seed 11, 3,000 random well-formed forms, as its console "
                         "gave it; a completed model is counted 'disconnected' when EPANET's 24 h report names a node "
                         "disconnected",
                 "forms": 3000,
                 "before": {"ok": 285, "ok_but_disconnected": 34, "refused": 2457, "blocked": 223, "form_refused": 1},
                 "after_reachability": {"ok": 258, "ok_but_disconnected": 8, "refused": 2521, "blocked": 212,
                                        "form_refused": 1},
                 "after_the_tank_check": {"ok": 255, "ok_but_disconnected": 1, "refused": 2537, "blocked": 206,
                                          "form_refused": 1},
                 "left": "a tank on the only path to a branch, low enough to fill, fills to its maximum and then "
                         "closes, so EPANET cuts the branch off: when a tank fills or empties is not checked"}


def final_long_header_form(n: int = FINAL_LONG_HEADER) -> bytes:
    """The final verification's slow form: one header cell with an opening bracket, n spaces and no closing bracket
    (it took over 100 s to read before the fix)."""
    return ("kind,id,location,head_ft,gallons_per_day,Dose (" + " " * n + "x\n" + 'source,WELL,"0,-800",200,,1.2\n'
            "production,,,,500000,\n").encode()


def part_form_final(paths: dict) -> dict:
    """Task 19's final verification (operator_form_final_verification): each confirmed issue's case, on the review's
    own drawing (fixtures.FORM_REVIEW_*), and what the form does now.  Timings are not kept here (the check times the
    slow form itself)."""
    import io
    import zipfile
    from .operator_form import (FormRefused, Skeleton, _canonical, _Geometry, _load, apply_operator_form,
                                first_hour_error, form_csv, read_model_file, read_operator_form, skeleton_from_build,
                                write_text)
    from . import fixtures as FX
    _quiet()
    out = {}
    b = _review_build(paths, "form_review_r12.dxf")
    sk = skeleton_from_build(b, "cad")
    base = [REVIEW_SRC, REVIEW_USE]

    def apply(rows, s=sk):
        return apply_operator_form(s, _review_form(rows))

    # 1. a header with an opening bracket and a run of spaces (three competing runs of spaces took cubic time)
    f1 = read_operator_form(final_long_header_form(), "form.csv")
    out["bracket_header"] = {
        "form_bytes": len(final_long_header_form()), "rows": len(f1.rows), "refused": f1.refused, "notes": f1.notes,
        "headers_read": {h: _canonical(h) for h in ("Dose (mg/L)", "Head [ft]", "Base elevation (ft)", "Dose(mg/L)",
                                                    "Dose ( mg / L )", "Gallons per day (gpd)", "Setting (psig)",
                                                    "Dose (" + " " * 40 + "x", "Dose ((mg/L))")}}

    # 2. a branch the form itself cuts off (its only pipe, (3000,1000)-(3000,2000), given by its midpoint), and the
    # same rows where they cut nothing off
    pump = {"kind": "pump", "id": "BP", "pipe": "3000,1500", "power_hp": "10"}
    cut = {"a closed_valve row on its only pipe": [{"kind": "closed_valve", "pipe": "3000,1500"}],
           "a check valve pointing out of it": [{"kind": "check_valve", "pipe": "3000,1500", "upstream": "3000,2000"}],
           "a pump whose rule is never": [dict(pump, upstream="3000,1000"),
                                          {"kind": "pump_rule", "pump": "BP", "on_hours": "never"}],
           "a pump pointing out of it": [dict(pump, upstream="3000,2000"),
                                         {"kind": "pump_rule", "pump": "BP", "on_hours": "always"}]}
    out["cut_off_by_the_form"] = {k: _outcome(apply(base + v)) for k, v in cut.items()}
    fed = {"a check valve pointing into it": [{"kind": "check_valve", "pipe": "3000,1500", "upstream": "3000,1000"}],
           "an always-on pump pointing into it": [dict(pump, upstream="3000,1000"),
                                                  {"kind": "pump_rule", "pump": "BP", "on_hours": "always"}],
           "a pump on the clock pointing into it": [dict(pump, upstream="3000,1000"),
                                                    {"kind": "pump_rule", "pump": "BP", "on_hours": "6-22"}],
           "a closed valve on a loop's pipe": [{"kind": "closed_valve", "pipe": "1500,0"}]}
    fed_r = {k: apply(base + v) for k, v in fed.items()}
    out["fed_through_open_or_one_way_links"] = {k: {**_outcome(r), "part_time_warning": [
        w for w in r.warnings if "holds water for those hours" in w]} for k, r in fed_r.items()}

    # 3. clock rules: minutes are refused (WNTR writes 'AT CLOCKTIME 6.5', which its reader refuses); whole hours,
    # '06:00-22:00' among them, run, and at 08:00 and 23:00 starts EPANET's pump status follows each rule hour by hour
    pa = dict(REVIEW_PUMP, id="PA")
    rule = lambda h: {"kind": "pump_rule", "pump": "PA", "on_hours": h}
    out["clock_rules"] = {"refused": {h: apply(base + [REVIEW_TANK, pa, rule(h)]).refused
                                      for h in ("6:30-22:00", "6:15-22", "06:00-22:30")},
                          "runs": {h: apply(base + [REVIEW_TANK, pa, rule(h)]).ok for h in ("06:00-22:00", "6 to 22")}}
    edges = {}
    for start, start_h in (("8:00 AM", 8.0), ("11:00 PM", 23.0)):
        s = Skeleton(inp_text=re.sub(r"(Start ClockTime\s+)12 am", rf"\g<1>{start}", sk.inp_text), name="zone.inp")
        for h in ("0-6", "18-24", "22-6", "6-22"):
            r = apply(base + [REVIEW_TANK, pa, rule(h)], s)
            st = [int(x) for x in _review_run(r.inp_text, 48).link["status"]["PA"].values]
            a, z = (int(x) % 24 for x in h.split("-"))
            want = [int((a <= (start_h + k) % 24 < z) if a < z else ((start_h + k) % 24 >= a or (start_h + k) % 24 < z))
                    for k in range(len(st))]
            edges[f"{start} {h}"] = [k for k, (x, y) in enumerate(zip(st, want)) if x != y]
    out["clock_rules"]["model_hours_where_the_pump_differs_from_its_rule"] = edges

    # 4. messages that named a temporary file: WNTR's reader failing in the first hour and in the last guard
    p0 = sorted(_load(sk.inp_text).pipe_name_list)[0]
    half = sk.inp_text.replace("[CONTROLS]", f"[CONTROLS]\nLINK {p0} OPEN AT CLOCKTIME 6.5", 1)
    assert half != sk.inp_text
    fh_msg = first_hour_error(half)
    guard = apply(base, Skeleton(inp_text=half, name="half.inp"))
    out["messages_name_no_file"] = {
        "first_hour": fh_msg, "last_guard": guard.blocking,
        "any_path": any(x in m for m in [fh_msg or ""] + guard.blocking for x in ("/var/", "/tmp", "rm_form", ".inp'"))}

    # 5. ids that start with '[' (a data line that starts with '[' is a section header to WNTR and EPANET)
    ids = {}
    for label, rows in (("pump [END]", [REVIEW_TANK, dict(REVIEW_PUMP, id="[END]"), {"kind": "pump_rule", "pump": "[END]",
                                                                                     "on_hours": "always"}]),
                        ("pump [P1]", [REVIEW_TANK, dict(REVIEW_PUMP, id="[P1]"), {"kind": "pump_rule", "pump": "[P1]",
                                                                                   "on_hours": "always"}]),
                        ("prv [END]", [dict(REVIEW_PRV, id="[END]")])):
        r = apply(base + rows)
        back = _load(r.inp_text) if r.inp_text else None
        ids[label] = {"ok": r.ok, "named": sorted(set(r.wn.pump_name_list) | set(r.wn.valve_name_list)) if r.wn else None,
                      "links_written_and_read_back": [len(r.wn.link_name_list) if r.wn else None,
                                                      len(back.link_name_list) if back else None]}
    out["ids_starting_with_a_bracket"] = ids

    # 6. a point with an enormous coordinate (the KD-tree found no neighbour, and its index was one past the last node)
    geo = _Geometry(_load(sk.inp_text))
    out["far_points"] = {
        "nearest_node": {f"{y:g}": geo.node((0.0, y))[0] for y in (1e150, 1e160, 1e200)},
        "refused": {loc[:24]: apply(base + [{"kind": "large_user", "location": loc, "gallons_per_day": "100"}]).refused
                    for loc in ("0,1e150", "0,1e160", "2000," + "1" + "0" * 180)},
        "pipe_refused": apply(base + [{"kind": "closed_valve", "pipe": "0,1e160"}]).refused}

    # 7. rows taken with ok that changed nothing or that another row overrode
    done = apply(base)
    completed = Skeleton(inp_text=done.inp_text, name="completed.inp")
    el = lambda loc, ft: {"kind": "elevation", "location": loc, "elevation_ft": ft}
    prv60 = dict(REVIEW_PRV, id="V1", elevation_ft="60")
    out["silent_rows"] = {
        "a pattern row on a model with water use": apply(
            [{"kind": "pattern", "multipliers": " ".join(f"{x:g}" for x in REVIEW_PATTERN)}], completed).refused,
        "two elevation rows for one node": apply(base + [el("1000,0", "120"), el("1000,0", "20")]).refused,
        "an elevation row at the source's node": apply(base + [el("0,-800", "35")]).refused,
        "an elevation row against a PRV's elevation_ft": apply(base + [prv60, el("2000,0", "5")]).refused,
        "an elevation row that agrees with it": _outcome(apply(base + [prv60, el("2000,0", "60")]))["ok"],
        "an elevation row at the PRV's other side": _outcome(apply(base + [prv60, el("1000,0", "5")]))["ok"]}

    # 8. pipe_c rows counted over the pipes still in the model, once a PRV and a pump replace two of them
    r8 = apply(base + [{"kind": "pipe_c", "material": "all", "hw_c": "110"}, dict(REVIEW_PRV, id="V1"),
                       {"kind": "pump", "id": "PA", "pipe": "0,-400", "upstream": "0,-800", "power_hp": "30"},
                       rule("always")])
    out["pipe_c_counted_after_step_4"] = {"ok": r8.ok, "pipes": len(r8.wn.pipe_name_list),
                                          "readiness": next(v for i, v, _ in r8.rows if i == "Hazen-Williams C")}

    # 9. pipe_c rows on a model whose headloss formula is Darcy-Weisbach (its roughness is a height, not a C)
    wd = _load(re.sub(r"(Headloss\s+)H-W", r"\1D-W", sk.inp_text))
    for pn in wd.pipe_name_list:
        wd.get_link(pn).roughness = 0.0015 * FT_M
    sdw = Skeleton(inp_text=write_text(wd), name="dw.inp")
    out["pipe_c_on_darcy_weisbach"] = {"refused": apply(base + [{"kind": "pipe_c", "material": "all", "hw_c": "110"}],
                                                        sdw).refused, "without_the_row_ok": apply(base, sdw).ok}

    # 10. values the row checks had not covered
    tank_rule = lambda on, off: [REVIEW_TANK, pa, {"kind": "pump_rule", "pump": "PA", "tank": "TANK",
                                                   "on_below_ft": on, "off_above_ft": off}]
    out["values_refused"] = {
        "head_ft 490 500": apply([dict(REVIEW_SRC, head_ft="490 500"), REVIEW_USE]).refused,
        "gallons_per_day 1_500": apply([REVIEW_SRC, {"kind": "production", "gallons_per_day": "1_500"}]).refused,
        "on below 2 ft, the tank's minimum 5 ft": apply(base + tank_rule("2", "25")).refused,
        "off above 40 ft, the tank full at 30 ft": apply(base + tank_rule("10", "40")).refused,
        "a tank's base at 646000 ft": apply(base + [dict(REVIEW_TANK, base_elev_ft="646000")]).refused,
        "a source's head at 490500 ft": apply([dict(REVIEW_SRC, head_ft="490500"), REVIEW_USE]).refused}
    out["values_kept"] = {"head_ft 1,200": apply([dict(REVIEW_SRC, head_ft="1,200"), REVIEW_USE]).ok,
                          "a rule inside the tank's range": apply(base + tank_rule("10", "25")).ok}

    # 11. a model in degrees: 50 ft in every direction, on a tangent plane (east and west it had been about 40 ft)
    lon0, lat0, fx, fy = REVIEW_DEG
    wg = _load(sk.inp_text)
    for n in wg.node_name_list:
        x, y = wg.get_node(n).coordinates
        wg.get_node(n).coordinates = (lon0 + x / fx, lat0 + y / fy)
    sdeg = Skeleton(inp_text=write_text(wg), name="degrees.inp")
    src_deg = dict(REVIEW_SRC, location=f"{lon0:.9f},{lat0 - 800 / fy:.9f}")
    at = (lon0 + 2000 / fx, lat0 + 1000 / fy)              # the junction at (2000,1000) ft
    deg = {}
    for label, (dx, dy) in (("45 ft east", (45, 0)), ("45 ft north", (0, 45)), ("55 ft east", (55, 0)),
                            ("55 ft north", (0, 55))):
        r = apply([src_deg, REVIEW_USE, {"kind": "large_user", "location": f"{at[0] + dx / fx:.9f},{at[1] + dy / fy:.9f}",
                                         "gallons_per_day": "1000"}], sdeg)
        deg[label] = {"placed": r.placed.get(4), "refused": [x["reason"] for x in r.refused]}
    out["degrees_both_ways"] = deg

    # 11b. a tank no source can fill (no pump, its water above the source's head), which the grammar fuzzer's last
    # disconnected models were: it drains and never fills, and what only it feeds is then cut off
    hi = dict(REVIEW_TANK, base_elev_ft="300")
    out["tank_no_source_can_fill"] = {"refused": apply(base + [hi]).refused,
                                      "with_a_pump_to_fill_it": _outcome(apply(base + [hi, pa, rule("always")]))["refused"]}
    out["grammar_fuzz"] = FINAL_FUZZ_19

    # 12. other files given as the form, or as the model on the command line
    zb = io.BytesIO()
    with zipfile.ZipFile(zb, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("zone.inp", sk.inp_text)
    words = {}
    for label, (name, data) in (("an EPANET model", ("model.inp", sk.inp_text.encode())),
                                ("a DXF drawing", ("mains.dxf", FX.form_review_r12_dxf())),
                                ("a GeoJSON layer", ("mains.geojson", FX.form_review_mains_geojson()))):
        try:
            read_operator_form(data, name)
            words[label] = None
        except FormRefused as e:
            words[label] = str(e)
    model = {}
    for label, (name, data) in (("a zip holding the model", ("zone.zip", zb.getvalue())),
                                ("a DXF drawing", ("mains.dxf", FX.form_review_r12_dxf())),
                                ("a GeoJSON layer", ("mains.geojson", FX.form_review_mains_geojson())),
                                ("a form", ("form.csv", form_csv(base)))):
        try:
            model[label] = {"read_as_the_model": read_model_file(data, name) == sk.inp_text}
        except FormRefused as e:
            model[label] = {"refused": str(e)}
    out["files_given_as_the_form"] = words
    out["files_given_as_the_model"] = model
    return out


FT_M = 0.3048
# Task 19's review: the public models' own C values and relative diffusivity, which the ladder's material-only C rung
# and every committed ky4 and ky10 chlorine result inherit (part 'public_models').  The wall-rate step is the old-pipe
# factor of C 130 against C 150 at gamma 1, 2^(20/30), applied to every pipe's wall rate with the hydraulics unchanged.
DIFF_KW_STEP = 2.0 ** ((150.0 - 130.0) / 30.0)


def part_form_public_models(paths: dict) -> dict:
    """Net3's, ky4's and ky10's Hazen-Williams C values (how many pipes carry each, and the new-pipe C the ladder's
    c_material rung gives each, intake_exp.material_only_c) and their files' relative diffusivity; and, for ky4 and
    ky10 at their truths (experiment.NET_TRUTH), at the file's diffusivity and at chlorine's own (1.0): how far a wall
    rate DIFF_KW_STEP times higher moves each junction's daily minimum on the scored day, and how many junctions are
    below 0.2 mg/L at the daily minimum."""
    from collections import Counter
    import wntr
    from .experiment import NET_TRUTH
    from .simulate import LIB, simulate_nominal_chlorine
    _quiet()
    out = {"kw_step": DIFF_KW_STEP, "threshold_mgL": FID_THRESHOLD}
    for net, path in (("Net3", os.path.join(LIB, "Net3.inp")), ("ky4", os.path.join(LIB, "ky4.inp")),
                      ("ky10", paths["ky10_order_fixed.inp"])):
        wn = wntr.network.WaterNetworkModel(path)
        cs = Counter(f"{p.roughness:g}" for _, p in wn.pipes())
        row = {"c_counts": dict(sorted(cs.items(), key=lambda kv: -float(kv[0]))),
               "material_only_c": {c: material_only_c(float(c)) for c in sorted(cs, key=float, reverse=True)},
               "file_diffusivity": wn.options.quality.diffusivity, "junctions": len(wn.junction_name_list)}
        if net != "Net3":
            t = {"kb_per_day": 0.40, "kw_m_per_day": 0.70, "source_dose": 1.2, **NET_TRUTH.get(net, {})}
            row["truth"] = t
            with tempfile.TemporaryDirectory(prefix="rm_diff_") as d:
                for label, dif in (("at_file_diffusivity", row["file_diffusivity"]), ("at_chlorine_diffusivity", 1.0)):
                    w = wntr.network.WaterNetworkModel(path)
                    w.options.quality.diffusivity = dif
                    p = _write_wn(w, d, f"{net}_{label}.inp")
                    mins = [simulate_nominal_chlorine(p, t["kb_per_day"], t["kw_m_per_day"], 1.0,
                                                      source_dose=t["source_dose"], kw_scale=k,
                                                      file_prefix=os.path.join(d, f"c{i}")).min()
                            for i, k in enumerate((1.0, DIFF_KW_STEP))]
                    ch = (mins[0] - mins[1]).abs()
                    row[label] = {"diffusivity": dif, "daily_min_change_mgL": {"median": float(ch.median()),
                                                                               "p95": float(ch.quantile(0.95))},
                                  "junctions_low": int((mins[0] < FID_THRESHOLD).sum())}
        out[net] = row
    return out


def form_record() -> dict:
    from . import fixtures as FX
    with tempfile.TemporaryDirectory(prefix="rm_t19_") as d:
        paths = FX.build_all(os.path.join(d, "fx"), FX.FORM_FIXTURES)
        out = {"about": "Task 19's operator form, converter self-consistency on public models (WNTR's bundled ky4 and "
                        "ky10, rewritten by fixtures.py), never a test on a utility's own files: ky4's clean drawing "
                        "completed by a form built from the original, against the pre-registered bars; ky10's PRVs "
                        "restored and its boundary valve closed by the form, against the original; the refusals. "
                        "Written by python -m residualmap.intake_exp form.",
               "complete": part_form_complete(paths), "prvs": part_form_prvs(paths),
               "boundary": part_form_boundary(paths), "no_source": part_form_no_source(paths),
               "bad_rows": part_form_bad_rows(paths), "after_the_stop": part_form_after_the_stop(paths),
               "review": part_form_review(paths), "public_models": part_form_public_models(paths),
               "final": part_form_final(paths), "scored_run": SCORED_RUN_19}
    return _r6(out)


# ----------------------------------------------------------------------------- task 19: the fidelity ladder
# Pre-registered on 23 October 2026, before any scored run (plan task 19 and its lead note on material-only C).
# The ladder measures what each thing a map of the mains lacks costs the product's own model (SimGP24: the 675-run
# grid with its demand and roughness axes and the GP, or its 225-run fixed-gamma slice where every pipe has one C, as
# importers.grid_for gives a file), on Net3, ky4 and ky10 (ky10 through task 15b's reaction-order fix, its prepared
# copy).  Each rung removes one more thing, cumulatively, from the network as the product reads it (rung r0):
#   r1 C 130 on every pipe; r2 plus demand by pipe length (the day's water use kept, one pattern: the one most
#   junctions carry, operator_form.length_shares); r3 plus one constant elevation on every junction (their mean, as
#   task 15b's elevation record); r4 plus every control and rule lost (pumps stay as the file starts them); r5 plus
#   every tank 1.3 times as wide; r6 (ky10 only) plus its 5 PRVs lost (each an open pipe of 1 ft at its diameter).
#   Beside the ladder, c_material (the lead note): from r0, each pipe's C replaced by the new-pipe C of the material its
#   C suggests (145 or more: plastic, topology.C_TABLE's 145; above 130: ductile iron, its 130; 130 or less: old cast
#   iron, its 135), so a skeleton's material-only C, blind to age, is measured on its own.
# The truth is always the full network's (experiment.NET_TRUTH; build_scenario with truth_age), seeds 0 to 7.  Each
# seed's 12 samples are drawn as experiment.run_scenario_time draws its first ones (default_rng(2000 + seed): 12
# junctions without replacement, hours from 07:00 to 17:00, the truth at that hour plus N(0, 0.03) mg/L, clipped at
# 0.01), the same samples on every rung.  Measured on every rung and seed: water-age MAE (h) of the model's nominal
# age against the truth's, over every junction and hour of the scored day (calibration does not change the model's
# age, so it is one number), and its MAE against the full model's own age; then before calibration (the prior,
# fit_prior) and after (the 12 samples), on the junctions not sampled: the daily minimum's MAE (mg/L), low spots
# found (truth's daily minimum below 0.2 mg/L and P(below) above 0.5) and false lows (P(below) above 0.5, truth not
# low).  The summary is the mean over the seeds.  Simulation on public models only.
FIDELITY = {n: os.path.join(OUT, f"fidelity_{n}.json") for n in ("Net3", "ky4", "ky10")}
FID_SEEDS = tuple(range(8))
FID_SAMPLES = 12
FID_NOISE_MGL = 0.03
FID_SAMPLE_HOUR = 14
FID_THRESHOLD = 0.2
FID_SAMPLE_RNG = 2000
FID_TANK_FACTOR = 1.3
FID_PRV_PIPE_FT = 1.0
FID_C_ALL = 130.0
FID_RUNGS = {"Net3": ("r0", "r1", "r2", "r3", "r4", "r5", "c_material"),
             "ky4": ("r0", "r1", "r2", "r3", "r4", "r5", "c_material"),
             "ky10": ("r0", "r1", "r2", "r3", "r4", "r5", "r6", "c_material")}
FID_WHAT = {"r0": "the network as the product reads it (the full model)",
            "r1": "C 130 on every pipe",
            "r2": "plus demand by pipe length (the day's water use kept, one pattern)",
            "r3": "plus one constant elevation on every junction (their mean)",
            "r4": "plus every control and rule lost (pump rules; pumps stay as the file starts them)",
            "r5": "plus every tank 1.3 times as wide",
            "r6": "plus its PRVs lost (each an open pipe of 1 ft)",
            "c_material": "beside the ladder, from the full model: C from the material alone (new-pipe values, blind "
                          "to age)"}


def _fid_c130(wn) -> None:
    for _, p in wn.pipes():
        p.roughness = FID_C_ALL


def _fid_demand(wn) -> None:
    from collections import Counter
    from .operator_form import length_shares
    total, count = 0.0, Counter()
    for _, j in wn.junctions():
        for ts in j.demand_timeseries_list:
            m = ts.pattern.multipliers if ts.pattern is not None else [1.0]
            total += ts.base_value * float(np.mean(m))
            count[ts.pattern_name] += 1
    main = count.most_common(1)[0][0]
    mean_m = float(np.mean(wn.get_pattern(main).multipliers)) if main is not None else 1.0
    shares = length_shares(wn)
    for name, j in wn.junctions():
        j.demand_timeseries_list.clear()
        j.add_demand(total * shares.get(name, 0.0) / mean_m, main)


def _fid_elevation(wn) -> None:
    e = float(np.mean([j.elevation for _, j in wn.junctions()]))
    for _, j in wn.junctions():
        j.elevation = e


def _fid_controls(wn) -> None:
    for c in list(wn.control_name_list):
        wn.remove_control(c)


def _fid_tanks(wn) -> None:
    for _, t in wn.tanks():
        t.diameter = t.diameter * FID_TANK_FACTOR


def _fid_prvs(wn) -> None:
    for v in [v for v in wn.valve_name_list if str(wn.get_link(v).valve_type) == "PRV"]:
        lk = wn.get_link(v)
        s, e, d, ml = lk.start_node_name, lk.end_node_name, lk.diameter, lk.minor_loss
        wn.remove_link(v)
        wn.add_pipe(v, s, e, length=FID_PRV_PIPE_FT * 0.3048, diameter=d, roughness=FID_C_ALL, minor_loss=ml)


def material_only_c(c: float) -> float:
    """The new-pipe C (topology.C_TABLE) of the material a C suggests: plastic for 145 or more, ductile iron above
    130, old cast iron for 130 or less (the lead note's case: an old iron main given cast iron's new-pipe value)."""
    from . import topology as T
    word = "PVC" if c >= 145 else "DIP" if c > 130 else "CI"
    return T.hazen_c(word)[0]


def _fid_material(wn) -> None:
    for _, p in wn.pipes():
        p.roughness = material_only_c(p.roughness)


# The first builds' times, copied from the ladder's run log (23 October 2026, this Mac, 11 workers, other work running at
# times): later runs read every grid from the cache, so their own times are 0.  Not rerun; run costs, not records.
# Rung r0's grids were cached already for Net3 (built 21 September) and ky4 (the committed 'full' grid, about 15
# minutes when it was built); ky10's r0 is its prepared copy's first grid.
FID_FIRST_BUILD = {"what": "seconds per grid in the first build, and for the truths, copied from the ladder's run log",
                   "Net3": {"r0": "cached", "r1": 4.1, "r2": 4.0, "r3": 4.0, "r4": 4.0, "r5": 4.2, "c_material": 8.1,
                            "truths_8_seeds": 1.0},
                   "ky10": {"r0": 218.9, "r1": 67.1, "r2": 80.1, "r3": 75.6, "r4": 73.4, "r5": 75.2, "r6": 73.8,
                            "c_material": 228.3, "truths_8_seeds": 53.9},
                   "ky4": {"r0": "cached", "r1": 270.8, "r2": 316.8, "r3": 237.1, "r4": 193.6, "r5": 170.4,
                           "c_material": 709.4, "truths_8_seeds": 98.1}}
FID_EDITS = {"r1": _fid_c130, "r2": _fid_demand, "r3": _fid_elevation, "r4": _fid_controls, "r5": _fid_tanks,
             "r6": _fid_prvs, "c_material": _fid_material}


def _sha1_text(text: str) -> str:
    import hashlib
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def fidelity_original(net: str, work: str) -> tuple[str, str]:
    """(the name the product reads the full network by, its .inp path): Net3 and ky4 as bundled (so ky4's committed
    grid serves rung r0), ky10 as task 15b's prepared copy (its reaction orders set to 1)."""
    from .simulate import LIB
    if net == "ky10":
        path, rep = _prepare(os.path.join(LIB, "ky10.inp"), os.path.join(work, "prepared"))
        assert path is not None and rep.ok, rep.blocking
        return path, path
    return net, os.path.join(LIB, f"{net}.inp")


def fidelity_rung_files(net: str, work: str) -> dict:
    """Every rung's .inp, written by WNTR into work under a content name (fid_<sha1[:8]>_<net>_<rung>.inp, so a grid
    cache name never serves another file's grid): {rung: {'name': what the product reads, 'path', 'sha1',
    'uniform_c'}}, and the WNTR rewrite of r0 with no edit (rewrite_r0), whose water age is checked against r0's."""
    import wntr
    name0, path0 = fidelity_original(net, work)
    out = {}
    with open(path0, encoding="utf-8", errors="replace") as fh:
        t0 = fh.read()
    wn0 = wntr.network.WaterNetworkModel(path0)
    out["r0"] = {"name": name0, "path": path0, "sha1": _sha1_text(t0), "uniform_c": _uniform_c(wn0)}
    rungs = FID_RUNGS[net]
    for r in [x for x in rungs if x != "r0"] + ["rewrite_r0"]:
        wn = wntr.network.WaterNetworkModel(path0)
        edits = [] if r == "rewrite_r0" else ["c_material"] if r == "c_material" else \
            [x for x in rungs if x not in ("r0", "c_material") and int(x[1:]) <= int(r[1:])]
        for e in edits:
            FID_EDITS[e](wn)
        tmp = os.path.join(work, f"tmp_{net}_{r}.inp")
        wntr.network.write_inpfile(wn, tmp, units="GPM")
        with open(tmp, encoding="utf-8") as fh:
            # WNTR's header comments name the file it read and the time it wrote: left out (operator_form.WNTR_HEADER),
            # so a rung's file, its sha1 and its grid's cache name are the same on every run
            from .operator_form import WNTR_HEADER
            text = "".join(ln for ln in fh if not ln.startswith(WNTR_HEADER))
        h = _sha1_text(text)
        p = os.path.join(work, f"fid_{h[:8]}_{net}_{r}.inp")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.remove(tmp)
        out[r] = {"name": p, "path": p, "sha1": h, "uniform_c": _uniform_c(wn), "edits": edits}
    return out


def _uniform_c(wn) -> bool:
    return len({float(f"{p.roughness:.6g}") for _, p in wn.pipes()}) == 1


def _fid_age(path: str) -> "np.ndarray":
    """Water age (h) of a file on the scored day, hour x junction, run in a temporary folder."""
    import wntr
    from .simulate import _last_day, load
    wn = load(path)
    wn.options.quality.parameter = "AGE"
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        r = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "a"))
    return (_last_day(r.node["quality"], wn.junction_name_list) / 3600.0)


def _fid_samples(truth, seed: int):
    import pandas as pd
    from .simgp import DAY_HOURS
    rng = np.random.default_rng(FID_SAMPLE_RNG + seed)

    def observe(js, hs):
        t = np.array([truth.truth_by_hour.loc[h, j] for j, h in zip(js, hs)])
        return np.clip(t + rng.normal(0, FID_NOISE_MGL, len(js)), 0.01, None)
    js = list(rng.choice(truth.junctions, FID_SAMPLES, replace=False))
    hs = [int(h) for h in rng.choice(DAY_HOURS, FID_SAMPLES)]
    return pd.DataFrame({"junction": js, "hour": hs, "y": list(observe(js, hs))})


def _fid_metrics(pred, truth, mask) -> dict:
    t, m = truth.truth_daily_min.loc[mask], pred.loc[mask]
    tv, pv = t < FID_THRESHOLD, m["p_below"] > 0.5
    return {"chlorine_mae_mgL": float(np.mean(np.abs(m["median"] - t))), "lows_found": int((tv & pv).sum()),
            "false_lows": int((~tv & pv).sum()), "n_true_lows": int(tv.sum()), "flagged": int(pv.sum())}


def _r6(x):
    if isinstance(x, dict):
        return {k: _r6(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r6(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return round(float(x), 6)
    if isinstance(x, np.integer):
        return int(x)
    return x


def fidelity_ladder(net: str, cache_dir: str | None = None, seeds=FID_SEEDS, log=print) -> dict:
    """The ladder on one network (the pre-registered protocol above): every rung's file, grid (its cache file and
    sha256) and scores, per seed and averaged.  cache_dir: where the grids are read or built (default outputs/cache;
    the Net3 rerun of the saved check builds its own in a temporary folder).  Run times go to '_measured' only."""
    import time
    import pandas as pd
    import types
    from .experiment import NET_TRUTH
    from .features import build_features
    from .importers import grid_for
    from .simgp import SimGP24, grid_cache_path, simulator_grid_24h
    from .simulate import build_scenario, nominal_scenario
    _quiet()
    cache_dir = cache_dir or os.path.join(REPO, "outputs", "cache")
    truth_kw = dict(NET_TRUTH.get(net, {}))
    dose = float(truth_kw.get("source_dose", 1.2))
    measured = {"grid_s": {}, "score_s": {}, "truth_s": None,
                "first_build": {"what": FID_FIRST_BUILD["what"], **FID_FIRST_BUILD.get(net, {})}}
    work = tempfile.mkdtemp(prefix=f"rm_fid_{net}_")
    try:
        files = fidelity_rung_files(net, work)
        # the WNTR rewrite of r0 with no edit must not move the water age (every rung is a WNTR rewrite)
        a0, aw = _fid_age(files["r0"]["path"]), _fid_age(files["rewrite_r0"]["path"])
        rewrite = {"max_abs_age_change_h": float((a0 - aw).abs().values.max()),
                   "sha1_r0": files["r0"]["sha1"], "sha1_rewrite": files["rewrite_r0"]["sha1"]}
        t0 = time.time()
        truths = {s: build_scenario(files["r0"]["name"], seed=s, sample_hour=FID_SAMPLE_HOUR, truth_age=True,
                                    **truth_kw) for s in seeds}
        measured["truth_s"] = round(time.time() - t0, 1)
        samples = {s: _fid_samples(truths[s], s) for s in seeds}
        log(f"{net}: truths for {len(seeds)} seeds in {measured['truth_s']} s")
        age_full = None
        rungs = {}
        for r in FID_RUNGS[net]:
            f = files[r]
            sc = nominal_scenario(f["name"], sample_hour=FID_SAMPLE_HOUR, source_dose=dose)
            X = build_features(sc)
            grid = grid_for("free_chlorine", f["uniform_c"])
            gpath = grid_cache_path(types.SimpleNamespace(wn_name=f["name"], source_dose=dose), cache_dir, grid)
            t0 = time.time()
            params, Z = simulator_grid_24h(sc, cache_dir, grid)
            measured["grid_s"][r] = round(time.time() - t0, 1)
            log(f"{net} {r}: grid {grid} ({len(params)} runs) in {measured['grid_s'][r]} s")
            age = sc.age_by_hour_h
            if r == "r0":
                age_full = age
            t0 = time.time()
            per_seed = []
            for s in seeds:
                tr = truths[s]
                assert list(tr.junctions) == list(sc.junctions), (net, r)
                S = samples[s]
                mask = [j for j in sc.junctions if j not in set(S.junction)]
                m0 = SimGP24(sc, X, seed=s, cache_dir=cache_dir, grid=grid, threshold=FID_THRESHOLD,
                             grid_data=(params, Z)).fit_prior()
                b = _fid_metrics(m0.predict_daily_min(), tr, mask)
                m1 = SimGP24(sc, X, seed=s, cache_dir=cache_dir, grid=grid, threshold=FID_THRESHOLD,
                             grid_data=(params, Z)).fit(S)
                a = _fid_metrics(m1.predict_daily_min(), tr, mask)
                per_seed.append({"seed": s, "age_mae_h": float((age - tr.truth_age_by_hour_h[sc.junctions]).abs()
                                                               .values.mean()),
                                 "age_mae_vs_full_model_h": float((age - age_full).abs().values.mean()),
                                 "before": b, "after": a, "map_params": list(m1.map_params_),
                                 "map_dose": m1.map_dose_})
            measured["score_s"][r] = round(time.time() - t0, 1)
            mean = lambda f_: float(np.mean([f_(x) for x in per_seed]))
            summary = {"age_mae_h": mean(lambda x: x["age_mae_h"]),
                       "age_mae_vs_full_model_h": mean(lambda x: x["age_mae_vs_full_model_h"]),
                       "n_true_lows": mean(lambda x: x["after"]["n_true_lows"])}
            for when in ("before", "after"):
                for k in ("chlorine_mae_mgL", "lows_found", "false_lows", "flagged"):
                    summary[f"{when}_{k}"] = mean(lambda x: x[when][k])
            rungs[r] = {"what": FID_WHAT[r], "edits": f.get("edits", []), "file_sha1": f["sha1"],
                        "uniform_c": f["uniform_c"], "grid": grid, "runs": len(params),
                        "grid_file": os.path.basename(gpath), "grid_sha256": _grid_sha256(params, Z),
                        "per_seed": per_seed, "summary": summary}
            log(f"{net} {r}: age MAE {summary['age_mae_h']:.2f} h, after false lows {summary['after_false_lows']:.1f}, "
                f"found {summary['after_lows_found']:.1f} of {summary['n_true_lows']:.1f} ({measured['score_s'][r]} s)")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return _r6({"about": f"Task 19's fidelity ladder on {net}: what each thing a map of the mains lacks costs the "
                         f"product's own model (SimGP24 and its grid), removed one more at each rung, before and after "
                         f"calibration on 12 samples, mean over seeds {seeds[0]} to {seeds[-1]}. Simulation on a public "
                         f"model (WNTR's bundled network), never a test on a utility's own system. Written by python "
                         f"-m residualmap.intake_exp fidelity {net}.",
                "network": net, "truth": {**truth_kw, "seeds": list(seeds)},
                "protocol": {"samples": FID_SAMPLES, "noise_mgL": FID_NOISE_MGL, "sample_hours": "07:00 to 17:00",
                             "sample_rng": f"default_rng({FID_SAMPLE_RNG} + seed)", "threshold_mgL": FID_THRESHOLD,
                             "scored_on": "the junctions not sampled, against the truth's daily minimum",
                             "age_scored_on": "every junction and hour of the scored day",
                             "tank_factor": FID_TANK_FACTOR, "c_everywhere": FID_C_ALL,
                             "prv_pipe_ft": FID_PRV_PIPE_FT},
                "rewrite_check": rewrite, "rungs": rungs, "_measured": measured})


# The boundary valve of operator_form_boundary_valve, picked once (pre-registered 23 October 2026, before the search
# ran): among ky10's closable pipes (simulate.closable_pipes on its prepared copy: closing one keeps the network
# connected), the pipes whose closure keeps every junction fed (no junction whose oldest water on the scored day is
# BOUNDARY_FED_H or more, unless it already was with the pipe open), and of those the one whose closure moves water age
# the most (the largest mean over junctions of the absolute change in the scored day's mean age; ties by pipe ID).  It
# is then pinned by ID in fixtures.py (BOUNDARY_PIPE_KY10), so the check never repeats the search.
BOUNDARY_RECORD = os.path.join(OUT, "boundary_valve_ky10.json")
BOUNDARY_FED_H = 144.0                 # the scored day's start: water at least this old has not been renewed since


def _closed_age(args):
    """One 7-day water-age run with one pipe closed (None: none), in a temporary folder: the scored day's mean and
    maximum age at every junction (module level, for the search's worker processes)."""
    import wntr
    from .simulate import _last_day, load
    path, pipe = args
    wn = load(path)
    wn.options.quality.parameter = "AGE"
    if pipe is not None:
        wn.get_link(pipe).initial_status = wntr.network.LinkStatus.Closed
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d:
        r = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "b"))
    a = _last_day(r.node["quality"], wn.junction_name_list) / 3600.0
    return pipe, a.mean().values.astype(float), a.max().values.astype(float)


def boundary_valve_search(n_jobs: int | None = None) -> dict:
    import multiprocessing
    import time
    from concurrent.futures import ProcessPoolExecutor
    import wntr
    from .simulate import LIB, closable_pipes
    _quiet()
    work = tempfile.mkdtemp(prefix="rm_bv_")
    try:
        path, rep = _prepare(os.path.join(LIB, "ky10.inp"), os.path.join(work, "prepared"))
        cands = closable_pipes(wntr.network.WaterNetworkModel(path))
        t0 = time.time()
        _, mean0, max0 = _closed_age((path, None))
        n_jobs = n_jobs or max(1, (os.cpu_count() or 2) - 1)
        with ProcessPoolExecutor(max_workers=n_jobs, mp_context=multiprocessing.get_context("spawn")) as ex:
            got = list(ex.map(_closed_age, [(path, p) for p in cands], chunksize=4))
        secs = time.time() - t0
    finally:
        shutil.rmtree(work, ignore_errors=True)
    rows = []
    for pipe, mean1, max1 in got:
        starved = int(((max1 >= BOUNDARY_FED_H) & (max0 < BOUNDARY_FED_H)).sum())
        rows.append({"pipe": pipe, "starved_junctions": starved,
                     "mean_abs_change_h": float(np.abs(mean1 - mean0).mean()),
                     "max_abs_change_h": float(np.abs(mean1 - mean0).max())})
    fed = sorted([r for r in rows if r["starved_junctions"] == 0], key=lambda r: (-r["mean_abs_change_h"], r["pipe"]))
    return _r6({"about": "Task 19: the boundary valve of operator_form_boundary_valve, picked once on ky10's prepared "
                         "copy (public, bundled with WNTR) by the pre-registered rule in residualmap/intake_exp.py "
                         "(every closable pipe closed in turn, a 7-day water-age run each; of the closures that keep "
                         "every junction fed, the one that moves the scored day's mean age the most), then pinned in "
                         "fixtures.py. Written by python -m residualmap.intake_exp boundary.",
                "closable_pipes": len(cands), "fed_closures": len(fed), "fed_threshold_h": BOUNDARY_FED_H,
                "picked": fed[0] if fed else None, "next_four": fed[1:5],
                "starving_closures": len(rows) - len(fed), "_measured": {"seconds": round(secs, 1)}})


def fidelity_reproduce_grids(net: str, committed: dict, cache_dir: str | None = None, grids: bool = True) -> dict:
    """For the saved checks (ky4 and ky10, not rerun): every rung's file rebuilt by today's code has the committed
    sha1, and, with grids, every rung's cached grid has the committed sha256 (read from the cache, never built).
    grids=False reads nothing from the cache (fidelity_outputs_reproduce, in the committed report: task 19's final
    verification found it unpickling every git-ignored grid only to drop the sha256, so a grid cut short by an
    interrupted build failed a committed check); a cached grid that cannot be read is 'grid_readable' False."""
    import pickle
    cache_dir = cache_dir or os.path.join(REPO, "outputs", "cache")
    work = tempfile.mkdtemp(prefix=f"rm_fidchk_{net}_")
    try:
        files = fidelity_rung_files(net, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    out = {}
    for r, v in committed["rungs"].items():
        out[r] = {"file_sha1_same": files[r]["sha1"] == v["file_sha1"]}
        if not grids:
            continue
        p = os.path.join(cache_dir, v["grid_file"])
        got, readable = None, True
        if os.path.exists(p):
            try:
                with open(p, "rb") as fh:
                    params, Z = pickle.load(fh)
                got = _grid_sha256(params, Z)
            except Exception:  # noqa: BLE001  (a grid cut short or damaged: reported, never raised)
                readable = False
        out[r].update(grid_cached=os.path.exists(p), grid_readable=readable,
                      grid_sha256_same=got is not None and got == v["grid_sha256"])
    return out


def _unmeasured(part: dict) -> dict:
    """A part without its '_measured' key (peak memory and seconds, which vary from run to run)."""
    return {k: v for k, v in part.items() if k != "_measured"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="?", choices=("all", "elevation", "bundled", "record", "timing", "samples",
                                                "snapping", "zone", "gis", "cad", "cad_timing", "fidelity",
                                                "boundary", "form"),
                    default="all")
    ap.add_argument("nets", nargs="*", help="fidelity: the networks (Net3, ky4, ky10; default all three)")
    ap.add_argument("--ky10", action="store_true", help="timing: also time ky10's grid (about 4 minutes)")
    ap.add_argument("--fixed-only", action="store_true", help="timing: rewrite only the records copied from logs "
                                                              "(scored run, review run, feasibility review, standalone "
                                                              "builds) in the existing file, timing nothing")
    a = ap.parse_args(argv)
    if a.what == "timing" and a.fixed_only:
        _write(RUN_COST_TIMING, fixed_timing_records())
        return 0
    if a.what == "timing":
        _write(RUN_COST_TIMING, run_cost_timing(ky10=a.ky10))
        return 0
    if a.what == "samples":              # task 16 (not part of 'all', which is task 15b's three files)
        _write(RECORD_16, samples_record())
        return 0
    if a.what == "snapping":             # task 16b (not part of 'all')
        _write(RECORD_16B, snapping_record())
        return 0
    if a.what == "zone":                 # task 16c (not part of 'all'); a fresh 675-run Net3 grid, about 15 s
        _write(RECORD_16C, zone_record())
        return 0
    if a.what == "gis":                  # task 17 (not part of 'all'); about a minute
        _write(RECORD_17, gis_record())
        return 0
    if a.what == "cad":                  # task 18 (not part of 'all')
        _write(RECORD_18, cad_record())
        return 0
    if a.what == "cad_timing":           # task 18's review: the MTEXT timing (about two minutes; not part of 'all')
        _write(RECORD_18_TIMING, cad_mtext_timing())
        return 0
    if a.what == "fidelity":             # task 19: the ladder (Net3 minutes; ky10 about 17, ky4 about 35 with grids)
        for net in a.nets or ("Net3", "ky10", "ky4"):
            _write(FIDELITY[net], fidelity_ladder(net))
        return 0
    if a.what == "boundary":             # task 19: the boundary valve's one search (a few minutes)
        _write(BOUNDARY_RECORD, boundary_valve_search())
        return 0
    if a.what == "form":                 # task 19: the operator form's record
        _write(RECORD_19, form_record())
        return 0
    if a.what in ("all", "elevation"):      # first: the readiness text quotes it
        _write(ELEVATION_PRV, elevation_prv())
    if a.what in ("all", "bundled"):
        _write(READINESS_BUNDLED, readiness_bundled())
    if a.what in ("all", "record"):
        _write(RECORD_15B, record())
    return 0


if __name__ == "__main__":
    sys.exit(main())
