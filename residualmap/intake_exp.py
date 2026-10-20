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

Every file is built from public models (WNTR's bundled networks, edited by fixtures.py), so the results are converter
self-consistency on public models, never a test on a utility's own files.
"""
from __future__ import annotations

import argparse
import atexit
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="?", choices=("all", "elevation", "bundled", "record", "timing", "samples",
                                                "snapping"),
                    default="all")
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
    if a.what in ("all", "elevation"):      # first: the readiness text quotes it
        _write(ELEVATION_PRV, elevation_prv())
    if a.what in ("all", "bundled"):
        _write(READINESS_BUNDLED, readiness_bundled())
    if a.what in ("all", "record"):
        _write(RECORD_15B, record())
    return 0


if __name__ == "__main__":
    sys.exit(main())
