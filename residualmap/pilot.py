"""
pilot.py — the real-data path (docs/pilot_protocol.md made executable).

A pilot utility gives us their EPANET .inp and their grab-sample log.  The validation is a hold-out in
time: fit on the earlier months, predict the held-out month's samples at their own junction and hour,
report RMSE, MAE, coverage of the 50/80/90/95% bands, and the violation call, against two baselines an
operator could run by hand (last reading at the same tap; mean of all recent samples).

    python -m residualmap.pilot --inp their_model.inp --samples grab_log.csv --taps tap_map.csv --dose 1.2
    python -m residualmap.pilot --synthetic Net3          # six synthetic months, to see the report format
    python -m residualmap.pilot --inp ... --plant plant_log.csv   # with the plant's monthly water temperature (task 10)
    python -m residualmap.pilot --inp ... --samples total_log.csv --taps ... --disinfectant chloramine --dose 2.0 \
                                --ph 8.0 --cl2n 4.5
                                       # a chloraminated system (task 12): total chlorine, the chloramine grid,
                                       # threshold 0.5 mg/L total chlorine unless --threshold is given; with the
                                       # plant's pH and Cl2:N the bulk rate gets the chloramine mode's prior (model
                                       # (c), the one held to task 12's bars), without them a uniform prior (model (b))

Real data stays private (iteration 5, task 15): --inp goes through residualmap.importers.prepare_inp (identified by
its content, decoded to UTF-8, IDs checked before WNTR parses), and the prepared model is copied to the git-ignored
outputs/private/inputs/<sha1 of the file>[:8]_<stem>.inp, so a revised file sent under the same name gets its own grid.
A real run writes to outputs/private/pilot by default (git-ignored), with import_report_<tag>.json beside the results;
--synthetic still writes to outputs/pilot.  --inp with the name of a network bundled with WNTR (Net3, ky4, ...) runs
that network by its name, as before task 15, with its report and no private copy.

Task 15b: the readiness report checks the model's physics (prepare_inp).  A model of more than 3,000 junctions stops
with the report unless --allow-large is given; a grid build that would leave less than 1.5 GB of disk free is refused
with a plain message; a file with one C value on every pipe is calibrated on the 'full_fixed_gamma' grid (the old-pipe
factor held at 0, since it cannot be learned from one C value).  A model blocked now that an earlier version accepted
leaves no private copy behind.

Task 16: the grab log is read by importers.read_grab_log.  Its columns are guessed from the header (date and time, or
one datetime column; 'Sample Point' or tap_id; 'Free Chlorine (mg/L)' or free_chlorine_mgL), or given with --columns; it
may be in any encoding, separated by semicolons or tabs with decimal commas, or an Excel workbook (.xlsx); a date that
reads two ways stops the run with the question unless --date-format says; '<0.02' and 'ND' are listed and not fitted;
an online analyzer's series keeps one reading a day (--analyzer-hour, default 10:00); and each time of day becomes a
model hour, the time minus the model's START CLOCKTIME, mod 24 h (minutes truncated, as before), so a file that starts
at 06:00 reads a 09:00 sample as model hour 3, and the run prints that its outputs' hours are model hours.  A model
that starts at midnight (the bundled Net3 and ky4, and the committed example logs' runs) gives exactly the rows it gave
before.  --synthetic logs are model hours already.  After the review of task 16: --columns also takes analyte= and
unit= (a long-format export), a log with both species runs once --columns names its reading column,
--keep-every-reading keeps a tap that looks like an analyzer whole, and readings above the measuring range ('>2.2')
are noted apart from those below detection.  After the final verification: --columns also takes analyte_value=, the
one analyte whose rows are read, for a long-format export whose analytes do not name the species in words the reader
knows (a log is refused with that question rather than read whole).

Task 16c, zone extracts: --inlets inlet_residuals.csv (inlet_node_id, residual_mgL) gives each inlet of a zone its
measured residual; each entry point is dosed at r_i x --dose, r_i = residual_i / the largest, so --dose means the
highest residual entering the zone (default: the file's highest).  A tap map's role column can mark a tap as a zone
'inlet' (its junction column naming the inlet's node): the median of its readings over the whole log sets that
inlet's residual, and its rows never enter the calibration.  Every entry point needs a residual, or none.  The grid is
cached under its own name (simgp.grid_cache_path); outputs get '_inlets' in their names.  Not built with --plant.
--holdout-sites K adds leave-sites-out rows ('held_out_sites'): within each held-out month, the model is refitted on
the same earlier months without every sample at K of that month's sites, and those sites are predicted at their own
hours, rotating through all its sites in seeded groups of K; outputs get '_sites<K>' in their names.  Without these
flags the pilot is unchanged.  Its review: an inlet tap with no usable reading, or at a node that is no entry point, is
refused; the inlet file's notes are printed; per-inlet residuals are refused in the chloramine mode (--inlets with
--disinfectant chloramine, and inlet taps); a held-out-sites row's n_train is its smallest group fit's (n_train_min and
n_train_max give the range).

grab_log.csv  : date (YYYY-MM-DD), time (HH:MM), tap_id, free_chlorine_mgL   [optional: method, notes]
                or, for a chloraminated system (--disinfectant chloramine), total_chlorine_mgL instead (never both
                columns: a log that mixes the two species is refused; docs/example_grab_log_total.csv)
tap_map.csv   : tap_id, junction_id                                          [optional: flush_min, notes, role]
inlet_residuals.csv : inlet_node_id, residual_mgL                            [optional: notes]   (task 16c)
plant_log.csv : month (YYYY-MM), temp_C                                      [optional: toc_mgL, dose_mgL]
                With a plant log the model is seasonal.SeasonalSimGP24: every grab sample is explained at its own
                month's water temperature and the held-out month is predicted at its own.  dose_mgL, when given, is
                the month's plant dose; first-order decay makes it an exact offset against --dose.  toc_mgL is carried
                and not used: task 11 tested a TOC-aware model in simulation (residualmap/organics.py) and it stopped
                at its pre-registered stop rule, so it is not wired in (journal, task 11).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import re
import sys
import types

import numpy as np
import pandas as pd
from scipy.stats import norm

from .chemistry import DEFAULT_THRESHOLD_MGL, Chemistry, Warming
from .chloramine import CA_DOSE_MGL, ca_condition
from .contract import CLI as CONTRACT_CLI, GRAB_LOG_MAPPING_KEYS
from .features import build_features
from .importers import (ANALYZER_HOUR, DATE_ORDERS, INLET_APP_ROUTE, INLET_PILOT_ROUTE, INLET_ZONE_NOTE, LARGE_JUNCTIONS,
                        LogRefused,
                        Readiness, _cell_text, _clock, clock_label,
                        combine_inlet_residuals, dropped_lines, inlet_ratios_from, inlet_tap_residuals, prepare_inp,
                        read_grab_log, read_inlet_residuals, read_table, shows_clock, size_gate, stem as inp_stem)
from .simgp import DAY_HOURS, SimGP24, grid_cache_path
from .simulate import build_scenario, load, net_path, nominal_scenario

LEVELS = (50, 80, 90, 95)


SPECIES_COLUMN = {"free_chlorine": "free_chlorine_mgL", "chloramine": "total_chlorine_mgL"}

# iteration 5, task 15: where a real run's inputs and results go (relative to the working directory, like the
# synthetic default; outputs/private/ is git-ignored)
PRIVATE_INPUTS = os.path.join("outputs", "private", "inputs")
PRIVATE_OUT = os.path.join("outputs", "private", "pilot")
SYNTHETIC_OUT = os.path.join("outputs", "pilot")


def default_out(synthetic: bool) -> str:
    """--out when it is not given: outputs/pilot for --synthetic (the committed demonstration), and the git-ignored
    outputs/private/pilot for a utility's real data, so a real run never writes into a tracked folder."""
    return SYNTHETIC_OUT if synthetic else PRIVATE_OUT


def prepare_pilot_inp(inp: str, disinfectant: str = "free_chlorine",
                      inputs_dir: str = PRIVATE_INPUTS) -> tuple[str | None, Readiness]:
    """--inp through importers.prepare_inp, then a private copy named by its content: inputs_dir/<sha1 of the prepared
    file>[:8]_<stem>.inp.  The grid cache is keyed by the model's file name, so two different files sent as zone.inp
    get two grids, and a revised file is never served the old one.  A file that needed a fix (or came in a zip) is
    already written there by prepare_inp, under its content name, and is used as it is (one copy).  The name of a
    network bundled with WNTR (Net3, ky4, ...) that is not a file here runs by that name, as before task 15: it is
    public, so it gets its report but no private copy, and its grid is the cached one; when the report fixes it
    (task 15b: ky10 and Net6 declare reaction orders of 0), the fixed copy in inputs_dir runs instead.  (None, report)
    when the file is not a model or has a blocking problem; nothing is copied then, and a private copy an earlier
    version made of the same file is removed (prepare_inp)."""
    if not os.path.exists(inp) and os.path.exists(net_path(inp)):
        path, rep = prepare_inp(net_path(inp), inputs_dir, disinfectant)
        if path is None:
            return None, rep
        return (inp if rep.same_file else path), rep
    if not os.path.exists(inp):
        raise FileNotFoundError(f"{inp}: no such file, and not the name of a network bundled with WNTR")
    path, rep = prepare_inp(inp, inputs_dir, disinfectant)
    if path is None or not rep.same_file:
        return path, rep
    with open(path, "rb") as fh:
        data = fh.read()
    dest = os.path.join(inputs_dir, f"{hashlib.sha1(data).hexdigest()[:8]}_{inp_stem(inp)}.inp")
    if os.path.abspath(path) != os.path.abspath(dest) and not os.path.exists(dest):
        os.makedirs(inputs_dir, exist_ok=True)
        with open(dest, "wb") as fh:
            fh.write(data)
    rep.prepared = os.path.basename(dest)
    return dest, rep


def read_log(samples: str, taps: str, disinfectant: str = "free_chlorine", date_format: str | None = None,
             start_clock_h: float = 0.0, analyzer_hour: float = ANALYZER_HOUR, columns: dict | None = None,
             thin_analyzers: bool = True):
    """The grab log and tap map read for the model (iteration 5, task 16: importers.read_grab_log): any encoding,
    separator, header and date format the reader recognises, CSV or Excel, censored readings kept and listed, online
    analyzer series thinned, and each time of day converted to a model hour with the file's START CLOCKTIME
    (start_clock_h; 0, a file that starts at midnight, changes no hour).  Returns the GrabLog; raises ValueError
    (importers.LogRefused) when the log cannot be read as it is."""
    with open(samples, "rb") as fh:
        log_data = fh.read()
    with open(taps, "rb") as fh:
        taps_data = fh.read()
    return read_grab_log(log_data, os.path.basename(samples), taps_data, os.path.basename(taps), disinfectant,
                         date_format=date_format, mapping=columns, start_clock_h=start_clock_h,
                         analyzer_hour=analyzer_hour, thin_analyzers=thin_analyzers)


def load_log(samples_csv: str, taps_csv: str, disinfectant: str = "free_chlorine", date_format: str | None = None,
             start_clock_h: float = 0.0, analyzer_hour: float = ANALYZER_HOUR, columns: dict | None = None,
             quiet: bool = False, thin_analyzers: bool = True) -> pd.DataFrame:
    """Join the grab log to the tap map -> junction, hour, y, month, date, tap_id (hour a model hour).  The log's
    reading column names its species: free chlorine (free_chlorine_mgL, 'Free Chlorine (mg/L)') or total chlorine
    (total_chlorine_mgL, a chloraminated system).  A log with both is refused (one likelihood must never mix species),
    and so is a log whose species does not match `disinfectant`.  Since task 16 the columns are guessed from the header
    (or given in columns), the file may be in any encoding, separated by semicolons or tabs with decimal commas, or an
    Excel workbook; a date that reads two ways needs date_format ('MDY', 'DMY', 'YMD' or a strptime format); readings
    below detection ('<0.02', 'ND') are listed and not used; and each time of day becomes a model hour with
    start_clock_h.  The committed example logs give the same rows as before task 16.  Prints what it dropped, censored
    (below detection and above the measuring range, each under its own words), thinned or warns about (quiet=True:
    nothing); thin_analyzers=False keeps every reading of a tap that looks like an online analyzer."""
    if disinfectant not in SPECIES_COLUMN:
        raise ValueError(f"disinfectant must be one of {sorted(SPECIES_COLUMN)}")
    g = read_log(samples_csv, taps_csv, disinfectant, date_format, start_clock_h, analyzer_hour, columns,
                 thin_analyzers)
    if not quiet:
        print_log(g)
    return g.frame


def print_log(g) -> None:
    """What load_log prints about a GrabLog: the rows dropped, censored, thinned or warned about (task 16)."""
    unknown = sorted({d["tap_id"] for d in g.dropped if d["reason"].endswith("is not in the tap map")})
    if unknown:     # at most 20 named (the final verification of task 16: a tap map from another system)
        n = sum(d["reason"].endswith("is not in the tap map") for d in g.dropped)
        print(f"warning: {n} samples at taps with no junction in the tap map are dropped: {unknown[:20]}"
              + (f" and {len(unknown) - 20} more taps" if len(unknown) > 20 else ""))
    for line in dropped_lines([d for d in g.dropped if not d["reason"].endswith("is not in the tap map")]):
        print(f"warning: not used: {line}")
    for kind, words in (("below", "below detection"), ("above", "above the measuring range")):
        c = g.censored[g.censored.kind == kind] if len(g.censored) else g.censored
        if len(c):
            print(f"note: {len(c)} reading{'s' if len(c) != 1 else ''} {words} {'are' if len(c) != 1 else 'is'} "
                  f"kept and listed, not used in the fit: "
                  + "; ".join(f"row {r.row} tap {r.tap_id} {r.raw!r}" for r in c.itertuples()))
    for n_ in g.notes[2:]:
        print(f"note: {n_}")
    for w in g.warnings:
        print(f"warning: {w}")


def synthetic_log(network: str = "Net3", months: int = 6, per_month: int = 10, rotating: int = 3, seed: int = 0,
                  schedule: list[dict] | None = None, return_truth: bool = False, **truth_kw):
    """Six 'months' of daytime grab samples: a fixed route of `per_month` taps plus `rotating` validation
    taps at different junctions each month.  The network's pipe-level truth is fixed; each month re-draws
    the operating truth (bulk decay ±20%, demand ±15% global and ±15% per node, dose ±10%).

    schedule=None is the committed log, draw for draw: month_seed = 100 + 10 seed + m (that formula repeats across
    seeds after ten months: seed 0's month 10 is seed 1's month 0).  A schedule (seasonal.plant_schedule: one dict per
    month, {'temp_C': the plant water temperature, 'soil_temp_C': None or the soil temperature for in-network warming})
    gives the seasonal path (task 10): month_seed = 1000 + 100 seed + m, distinct for any seed and up to 100 months;
    month m's truth is build_scenario(chem=Chemistry(temp_C=...), warming=Warming(soil_temp_C=...) when a soil
    temperature is given); every reading carries its month's logged plant temperature (temp_C); and the plant log
    (month, temp_C, toc_mgL, dose_mgL; TOC not logged here, the dose is the set point) comes back as a third item.
    Task 11: a schedule month may also give 'toc_mgL' (the plant TOC: the truth's Chemistry gets it, every reading of
    the month carries it, and the plant log logs it), 'kinetics' and 'phi' (the truth's kinetics, e.g. 'clark' with its
    demand per mg TOC; default first order) and 'dose_mgL' (the month's plant dose: the truth runs at that dose, and the
    readings and the plant log carry it).  A schedule without these keys gives task 10's log, draw for draw.
    return_truth=True appends a list with each month's truth: month (1-based), month_seed, truth_by_hour,
    truth_daily_min and chem (the truth's chemistry and hidden draws; None on the committed path)."""
    rng = np.random.default_rng(seed)
    if schedule is None:
        sc0 = build_scenario(network, seed=seed, **truth_kw)
        junctions = sc0.junctions
    else:
        if len(schedule) < months:
            raise ValueError(f"the schedule has {len(schedule)} months, {months} asked for")
        junctions = load(network).junction_name_list          # the same list build_scenario returns
    taps = list(rng.choice(junctions, per_month, replace=False))           # the utility's fixed route
    rows, truths, plant = [], [], []
    for m in range(months):
        label = f"{2026 + m // 12}-{m % 12 + 1:02d}"
        if schedule is None:
            ms = 100 + seed * 10 + m
            sc = build_scenario(network, seed=seed, month_seed=ms, **truth_kw)
        else:
            cond = schedule[m]
            ms = 1000 + 100 * seed + m
            warm = Warming(soil_temp_C=cond["soil_temp_C"]) if cond.get("soil_temp_C") is not None else None
            chem = Chemistry(temp_C=cond["temp_C"], toc_mgL=cond.get("toc_mgL"), kinetics=cond.get("kinetics", "first"),
                             phi=cond.get("phi"))
            kw_m = dict(truth_kw, source_dose=float(cond["dose_mgL"])) if cond.get("dose_mgL") is not None else truth_kw
            sc = build_scenario(network, seed=seed, month_seed=ms, chem=chem, warming=warm, **kw_m)
            plant.append({"month": label, "temp_C": float(cond["temp_C"]),
                          "toc_mgL": float(cond["toc_mgL"]) if cond.get("toc_mgL") is not None else float("nan"),
                          "dose_mgL": float(kw_m.get("source_dose", 1.2))})
        others = [j for j in sc.junctions if j not in taps]
        rot = list(rng.choice(others, rotating, replace=False)) if rotating else []
        for k, j in enumerate(taps + rot):
            h = int(rng.choice(DAY_HOURS)); day = int(rng.integers(1, 28))
            y = float(np.clip(sc.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None))
            row = {"junction": j, "hour": h, "y": round(y, 2), "month": label, "date": f"{label}-{day:02d}",
                   "tap_id": f"T{k + 1}" if k < per_month else f"V{m + 1}-{k - per_month + 1}"}
            if schedule is not None:
                row["temp_C"] = float(schedule[m]["temp_C"])
                for key in ("toc_mgL", "dose_mgL"):      # task 11; absent from task 10's schedules
                    if schedule[m].get(key) is not None:
                        row[key] = float(schedule[m][key])
            rows.append(row)
        if return_truth:
            truths.append({"month": m + 1, "month_seed": ms, "truth_by_hour": sc.truth_by_hour,
                           "truth_daily_min": sc.truth_daily_min, "chem": sc.chem})
    out = (pd.DataFrame(rows), network) + ((pd.DataFrame(plant),) if schedule is not None else ())
    return out + ((truths,) if return_truth else ())


PLANT_MONTH = re.compile(r"\d{4}-(?:0[1-9]|1[0-2])")     # a plant log's month: YYYY-MM
PLANT_NUMBERS = ("temp_C", "toc_mgL", "dose_mgL")         # read as numbers (a decimal comma too); the rest as text


def _plant_number(v) -> float | None:
    """A plant log cell as a number: None when blank; ValueError when it is text that is not a number."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return None
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v)
    s = _cell_text(v).strip()
    if not s:
        return None
    if not re.fullmatch(r"-?(?:\d+(?:[.,]\d+)?|[.,]\d+)", s):
        raise ValueError(s)
    return float(s.replace(",", "."))


def _plant_month(v) -> str | None:
    """A plant log's month as YYYY-MM (an Excel date cell gives its year and month), or None."""
    if isinstance(v, (datetime.date, datetime.datetime)):
        return f"{v.year:04d}-{v.month:02d}"
    s = _cell_text(v).strip()
    return s if PLANT_MONTH.fullmatch(s) else None


PLANT_TEMP_C = (0.0, 35.0)      # a plant log's temp_C range: Chemistry's own (chemistry.Chemistry refuses one outside it)


def load_plant(path: str) -> pd.DataFrame:
    """The plant log: month (YYYY-MM), temp_C (required), toc_mgL and dose_mgL (optional), one row per month.  Read
    as importers.read_table reads a grab log (any encoding, separator and decimal comma, or an Excel sheet; the task
    21 review: a plain read_csv stopped the pilot on a semicolon file with a misleading message, and on UTF-16 with a
    traceback).  A text whose first line holds month and temp_C takes that line as its header, and a row that leaves
    out the optional cells at its end is read with them blank, as read_csv read it before (the task 21 final
    verification).  Refused with a ValueError that names the file and the rows: a required column missing, a month
    not written YYYY-MM or listed twice, a temperature blank, not a number or outside PLANT_TEMP_C (Chemistry's range;
    when every one is above it, the message says it looks like Fahrenheit), or an optional value that is not a number
    (a blank one is left empty).  Other columns are kept as text."""
    name = os.path.basename(path)
    need = {"month", "temp_C"}
    with open(path, "rb") as fh:
        data = fh.read()
    try:
        t = read_table(data, name, header_names=tuple(sorted(need)))
    except LogRefused as e:
        raise ValueError(f"the plant log {name} cannot be read: {e}") from None
    if not need <= set(t.header):
        raise ValueError(f"the plant log {name} needs the columns month (YYYY-MM) and temp_C; its first row, read as "
                         f"the header, is {', '.join(t.header[:8])}")
    lo, hi = PLANT_TEMP_C
    cols, bad, temps = {h: [] for h in t.header}, [], []
    for r, row_no in zip(t.rows, t.row_numbers):
        for h, v in zip(t.header, r):
            if h == "month":
                m = _plant_month(v)
                if m is None:
                    bad.append(f"row {row_no}: month {_cell_text(v)!r} is not written YYYY-MM")
                cols[h].append(m)
            elif h in PLANT_NUMBERS:
                try:
                    x = _plant_number(v)
                    if x is None and h == "temp_C":
                        bad.append(f"row {row_no}: temp_C is blank")
                    elif h == "temp_C":
                        temps.append(x)
                        if not lo <= x <= hi:
                            bad.append(f"row {row_no}: temp_C {x:g} is outside {lo:g} to {hi:g} C")
                except ValueError:
                    bad.append(f"row {row_no}: {h} {_cell_text(v)!r} is not a number")
                    x = None
                cols[h].append(np.nan if x is None else x)
            else:
                cols[h].append(_cell_text(v))
    if bad:
        hint = (f" Every temperature is above {hi:g}, so they look like Fahrenheit: give temp_C in degrees Celsius, "
                f"(F - 32) / 1.8." if temps and all(x > hi for x in temps) else "")
        raise ValueError(f"the plant log {name} has {len(bad)} value{'s' if len(bad) != 1 else ''} that cannot be "
                         f"used: {'; '.join(bad[:6])}{'; ...' if len(bad) > 6 else ''}. Each row needs its month "
                         f"(YYYY-MM) and the water temperature in C, {lo:g} to {hi:g} (the range of every rate law "
                         f"used here).{hint}")
    p = pd.DataFrame(cols)
    if p.month.duplicated().any():
        raise ValueError(f"months listed twice in the plant log {name}: {sorted(p.month[p.month.duplicated()])}")
    return p


DEFAULT_DOSE_MGL = {"free_chlorine": 1.2, "chloramine": CA_DOSE_MGL}   # mg/L; typical chloramine doses are 1.5 to 4.0
# the review of task 16c: per-inlet residuals are not built for the chloramine mode until a check covers them
INLETS_NOT_CA = ("per-inlet residuals are not built for the chloramine mode yet: its dose axis (an effective dose of "
                 "0.65 to 1.10 times the dose) is made for the total chlorine leaving a plant, which the plant's fast "
                 "organic demand lowers, and a residual measured where water enters a zone has already passed that "
                 "demand.")
HOLDOUT_SITES_SEED = 50_000      # task 16c: the held-out site groups come from default_rng(HOLDOUT_SITES_SEED + seed)
HELD_OUT_SITES = "held_out_sites"


def validate(inp: str, log: pd.DataFrame, dose: float | None = None, threshold: float | None = None,
             window_months: int = 3, holdout_months: int | None = None, cache_dir: str = "outputs/cache", seed: int = 0,
             plant: pd.DataFrame | None = None, bank_cache: str = "readwrite",
             disinfectant: str = "free_chlorine", ph: float | None = None,
             cl2n: float | None = None, grid: str = "full", inlet_ratios: dict | None = None,
             holdout_sites: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rolling hold-out: for each held-out month, fit on the `window_months` before it and predict its
    samples.  Returns (per-sample predictions, per-month summary).
    plant: the plant log (load_plant).  None is the committed temperature-blind SimGP24, unchanged.  With a plant log,
    seasonal.SeasonalSimGP24 explains each sample at its month's water temperature (and dose, when dose_mgL is given
    for every month of the plant log; otherwise the dose is taken as --dose throughout) and predicts the held-out month
    at its own; the summary gains the month's temperature, kb20 and the posterior over
    the temperature hypotheses.  The bank of grids at the log's temperatures is built once and cached (bank_cache, as
    seasonal.covariate_bank's cache: 'read' never writes a cache file).
    disinfectant 'chloramine' (task 12): the log is total chlorine and the model is the chloramine mode's (the
    chloramine grid, its dose axis and likelihood scale).  With the plant's logged pH and Cl2:N (ph, cl2n) the bulk rate
    gets the chloramine mode's prior, model (c), the one held to task 12's bars; without them the prior is uniform,
    model (b).  A plant log is refused (no seasonal chloramine model is built).  dose and threshold default to the
    disinfectant's (DEFAULT_DOSE_MGL and chemistry.DEFAULT_THRESHOLD_MGL).
    grid: the free-chlorine grid without a plant log, 'full' (the committed path) or 'full_fixed_gamma' for a file
    with one C value on every pipe (task 15b; importers.grid_for).
    inlet_ratios (task 16c): {entry point: r}, a zone whose inlets carry different residuals, each dosed at r times
    dose (importers.inlet_ratios_from; dose is then the highest residual entering the zone); None, or every ratio 1,
    is the committed path.  Refused with a plant log (the seasonal grids are not built at inlet ratios) and, since the
    review of task 16c, with chloramine (INLETS_NOT_CA); checked as given against the model's entry points.
    holdout_sites (task 16c): None is the committed validation, row for row.  k adds leave-sites-out rows ('taps'
    'held_out_sites'): within each held-out month its sites (junctions) are split into groups of k, in an order drawn
    from default_rng(HOLDOUT_SITES_SEED + seed), and each group's samples of that month are predicted, at their own
    junction and hour, by a model refitted on the same earlier months with every sample at those sites removed, so
    every site of the month is predicted once by a fit that never saw it.  This is the map's accuracy at places
    without a tap that a fixed compliance route can give, where the 'new_taps' rows need rotating taps.  Their
    rmse_persistence is NaN (a held-out site has no reading of its own in training) and their map_* columns too (each
    group has its own fit); the predictions frame then gains a 'validation' column and the held-out sites' rows.
    Each group's fit trains on fewer samples than the month's rows (its sites' earlier samples are dropped): their
    n_train is the smallest group fit's, with n_train_min and n_train_max (the review of task 16c; it was the month's
    whole training set).  A group that leaves no training sample is refused, and k at least the month's sites warns."""
    if disinfectant not in DEFAULT_DOSE_MGL:
        raise ValueError(f"disinfectant must be one of {sorted(DEFAULT_DOSE_MGL)}")
    if holdout_sites is not None and (int(holdout_sites) != holdout_sites or holdout_sites < 1):
        raise ValueError("holdout_sites must be a whole number of sites, 1 or more")
    from .simulate import normalize_inlet_ratios
    ratios = normalize_inlet_ratios(inlet_ratios)
    if ratios is not None and plant is not None:
        raise NotImplementedError("per-inlet residuals with a plant log: the seasonal model's grids are not built at "
                                  "inlet ratios")
    if ratios is not None and disinfectant == "chloramine":      # the review of task 16c
        raise NotImplementedError(INLETS_NOT_CA)
    inlet_kw = {} if ratios is None else {"inlet_ratios": ratios}
    if disinfectant == "chloramine" and plant is not None:
        raise NotImplementedError("a plant log (temperature) with chloramine: no seasonal chloramine model is built")
    if (ph is None) != (cl2n is None):
        raise ValueError("give both the plant's pH and its Cl2:N ratio, or neither")
    if ph is not None and disinfectant != "chloramine":
        raise ValueError("pH and Cl2:N set the chloramine mode's prior; they are not used with free chlorine")
    dose = DEFAULT_DOSE_MGL[disinfectant] if dose is None else dose
    threshold = DEFAULT_THRESHOLD_MGL[disinfectant] if threshold is None else threshold
    sc = nominal_scenario(inp, sample_hour=14, source_dose=dose)
    if inlet_ratios:     # the review of task 16c: checked as given, before the all-ones shortcut, as the grid checks them
        from .simulate import check_inlet_ratios, source_nodes
        check_inlet_ratios(inlet_ratios, source_nodes(sc.wn))
    X = build_features(sc)
    unknown = sorted(set(log.junction) - set(sc.junctions))
    if unknown:
        raise ValueError(f"junction IDs not in the .inp: {unknown[:10]}")
    months = sorted(log.month.unique())
    held = months[-holdout_months:] if holdout_months else months[window_months:]
    bank = None
    if plant is not None:
        from .seasonal import SeasonalSimGP24, covariate_bank
        pl = plant.set_index(plant.month.astype(str))
        missing = sorted(set(months) - set(pl.index))
        if missing:
            raise ValueError(f"months in the grab log with no plant-log row: {missing}")
        ratio = (pl.dose_mgL.astype(float) / dose) if "dose_mgL" in pl and pl.dose_mgL.notna().all() else pd.Series(1.0, index=pl.index)
        log = log.assign(temp_C=log.month.map(pl.temp_C.astype(float)), dose_ratio=log.month.map(ratio))
        bank = covariate_bank(sc, sorted(set(log.temp_C)), cache_dir=cache_dir, cache=bank_cache)
    def fit_on(train: pd.DataFrame, m: str, warn: bool = True):
        """The month's model fitted on train (the committed constructors, plus the inlet ratios when given)."""
        if bank is None and disinfectant == "chloramine":
            from .chloramine import ca_condition, prior_log_vector
            from .simgp import DOSES_CA, LIK_SD_CA
            model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid="chloramine", cond=ca_condition(),
                            lik_sd=LIK_SD_CA, doses=DOSES_CA, threshold=threshold, **inlet_kw)
            if ph is not None:
                model.log_prior = prior_log_vector(model.params, float(ph), float(cl2n))
            model.fit(train[["junction", "hour", "y"]])
            if warn and model.map_dose_ in (DOSES_CA[0], DOSES_CA[-1]):
                print(f"warning: held-out month {m}: the fitted effective dose sits at the edge of the chloramine dose "
                      f"axis (x{model.map_dose_:.2f} of {dose:g} mg/L): check --dose, the total chlorine leaving the "
                      f"plant", file=sys.stderr, flush=True)
        elif bank is None and grid == "full":
            model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, **inlet_kw).fit(train[["junction", "hour", "y"]])
        elif bank is None:
            model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid=grid, **inlet_kw).fit(
                train[["junction", "hour", "y"]])
        else:
            model = SeasonalSimGP24(sc, X, bank, seed=seed, cache_dir=cache_dir).fit(
                train[["junction", "hour", "y", "temp_C", "dose_ratio"]], target_temp_C=float(pl.temp_C[m]),
                target_dose_ratio=float(ratio[m]))
        return model

    if holdout_sites:       # task 17's fix of task 16c's leftovers: a group of every site is refused before any fit
        for m in held:
            train_ = log[log.month.isin([t for t in months if t < m][-window_months:])]
            sites_ = sorted(set(map(str, log[log.month == m].junction)))
            if int(holdout_sites) >= len(sites_) and not (~train_.junction.isin(sites_)).any():
                raise ValueError(f"holdout_sites {int(holdout_sites)}: in held-out month {m} a group of sites leaves no "
                                 f"training sample (its {len(sites_)} sites are every site of the earlier months); use "
                                 f"fewer sites per group")
    rng_sites = np.random.default_rng(HOLDOUT_SITES_SEED + seed) if holdout_sites else None
    preds, rows, site_preds = [], [], []
    for m in held:
        train_months = [t for t in months if t < m][-window_months:]
        train, test = log[log.month.isin(train_months)], log[log.month == m]
        model = fit_on(train, m)
        z_mu, z_sd = model.predict_hours()
        jidx = np.array([sc.junctions.index(j) for j in test.junction]); h = test.hour.values
        mu, sd = z_mu[h, jidx], z_sd[h, jidx]
        p = test.assign(pred=np.exp(mu), z_mu=mu, z_sd=sd, p_below=norm.cdf((np.log(threshold) - mu) / sd),
                        train_months=",".join(train_months))
        for q in LEVELS:
            k = norm.ppf(0.5 + q / 200)
            p[f"in{q}"] = (test.y.values >= np.exp(mu - k * sd)) & (test.y.values <= np.exp(mu + k * sd))
        # baselines an operator could run by hand
        last = train.sort_values("date").groupby("tap_id").y.last()
        p["persistence"] = test.tap_id.map(last).fillna(train.y.mean()).values
        p["network_mean"] = train.y.mean()
        # a tap the model has seen in training tests the month-to-month forecast; a tap it has never seen
        # tests the MAP — the claim that matters
        p["seen_tap"] = test.junction.isin(set(train.junction)).values
        preds.append(p)
        for subset, mask in (("all", np.ones(len(p), bool)), ("seen_taps", p.seen_tap.values), ("new_taps", ~p.seen_tap.values)):
            if mask.sum() == 0:
                continue
            t, q_ = test.y.values[mask], p[mask]
            viol = t < threshold; flagged = q_.p_below.values > 0.5
            rows.append({"held_out_month": m, "taps": subset, "trained_on": ",".join(train_months), "n_train": len(train), "n_test": int(mask.sum()),
                         "rmse": float(np.sqrt(np.mean((t - q_.pred.values) ** 2))), "mae": float(np.mean(np.abs(t - q_.pred.values))),
                         "rmse_persistence": float(np.sqrt(np.mean((t - q_.persistence.values) ** 2))),
                         "rmse_network_mean": float(np.sqrt(np.mean((t - q_.network_mean.values) ** 2))),
                         **{f"coverage{q}": float(q_[f"in{q}"].mean()) for q in LEVELS},
                         "band90_width_mgL": float(np.median(np.exp(q_.z_mu + 1.645 * q_.z_sd) - np.exp(q_.z_mu - 1.645 * q_.z_sd))),
                         "n_true_below": int(viol.sum()), "recall_below": float((viol & flagged).sum() / viol.sum()) if viol.sum() else float("nan"),
                         "false_alarms": int((~viol & flagged).sum()),
                         "map_kb": model.map_params_[0], "map_kw": model.map_params_[1], "map_gamma": model.map_params_[2], "map_dose": model.map_dose_})
            if disinfectant == "chloramine":
                from .simgp import DOSES_CA
                rows[-1].update({"model": "(c) pH and Cl2:N prior" if ph is not None else "(b) uniform prior",
                                 "map_dose_at_axis_edge": bool(model.map_dose_ in (DOSES_CA[0], DOSES_CA[-1]))})
            if bank is not None:
                rows[-1].update({"temp_C": float(pl.temp_C[m]), "kb20": model.kb20(), "map_hypothesis": model.map_hypothesis_,
                                 **{f"P_{k}": v for k, v in model.hypothesis_posterior().items()}})
        if holdout_sites:        # task 16c: leave k sites out, rotating through every site of the month
            p["validation"] = "held_out_month"
            sites = sorted(set(map(str, test.junction)))
            order = [sites[i] for i in rng_sites.permutation(len(sites))]
            groups = [order[i:i + int(holdout_sites)] for i in range(0, len(order), int(holdout_sites))]
            # the review of task 16c: each group's fit trains on fewer samples than the month's (those sites' earlier
            # samples are dropped), recorded as n_train_min and n_train_max; a group that leaves no training sample
            # would be a fit on the prior alone, and is refused
            sizes = [int((~train.junction.isin(gr)).sum()) for gr in groups]
            if min(sizes) == 0:
                raise ValueError(f"holdout_sites {int(holdout_sites)}: in held-out month {m} a group of sites leaves no "
                                 f"training sample (its {len(sites)} sites are every site of the earlier months); use "
                                 f"fewer sites per group")
            if int(holdout_sites) >= len(sites):
                print(f"warning: held-out month {m}: --holdout-sites {int(holdout_sites)} is at least its {len(sites)} "
                      f"sites, so they are all held out together, in one group", file=sys.stderr, flush=True)
            parts = []
            for gi, group in enumerate(groups):
                tr_g, te_g = train[~train.junction.isin(group)], test[test.junction.isin(group)]
                mg = fit_on(tr_g, m, warn=False)
                gz_mu, gz_sd = mg.predict_hours()
                gj = np.array([sc.junctions.index(j) for j in te_g.junction]); gh = te_g.hour.values
                mu_g, sd_g = gz_mu[gh, gj], gz_sd[gh, gj]
                pg = te_g.assign(pred=np.exp(mu_g), z_mu=mu_g, z_sd=sd_g,
                                 p_below=norm.cdf((np.log(threshold) - mu_g) / sd_g), train_months=",".join(train_months))
                for q in LEVELS:
                    k_ = norm.ppf(0.5 + q / 200)
                    pg[f"in{q}"] = (te_g.y.values >= np.exp(mu_g - k_ * sd_g)) & (te_g.y.values <= np.exp(mu_g + k_ * sd_g))
                pg["persistence"] = np.nan          # no reading of a held-out site is in its fit
                pg["network_mean"] = tr_g.y.mean() if len(tr_g) else np.nan
                pg["seen_tap"] = False
                pg["validation"], pg["site_group"] = HELD_OUT_SITES, gi
                parts.append(pg)
            ps = pd.concat(parts)
            site_preds.append(ps)
            t = ps.y.values
            viol = t < threshold; flagged = ps.p_below.values > 0.5
            rows.append({"held_out_month": m, "taps": HELD_OUT_SITES, "trained_on": ",".join(train_months),
                         "n_train": min(sizes), "n_test": int(len(ps)),
                         "rmse": float(np.sqrt(np.mean((t - ps.pred.values) ** 2))),
                         "mae": float(np.mean(np.abs(t - ps.pred.values))), "rmse_persistence": float("nan"),
                         "rmse_network_mean": float(np.sqrt(np.mean((t - ps.network_mean.values) ** 2))),
                         **{f"coverage{q}": float(ps[f"in{q}"].mean()) for q in LEVELS},
                         "band90_width_mgL": float(np.median(np.exp(ps.z_mu + 1.645 * ps.z_sd) - np.exp(ps.z_mu - 1.645 * ps.z_sd))),
                         "n_true_below": int(viol.sum()),
                         "recall_below": float((viol & flagged).sum() / viol.sum()) if viol.sum() else float("nan"),
                         "false_alarms": int((~viol & flagged).sum()),
                         "map_kb": float("nan"), "map_kw": float("nan"), "map_gamma": float("nan"),
                         "map_dose": float("nan"), "sites": len(sites), "site_groups": len(groups),
                         "sites_per_group": int(holdout_sites), "n_train_min": min(sizes), "n_train_max": max(sizes)})
    return pd.concat(preds + site_preds, ignore_index=True), pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inp", help="the utility's EPANET model")
    ap.add_argument("--samples", help="grab log CSV: date,time,tap_id,free_chlorine_mgL (or total_chlorine_mgL)")
    ap.add_argument("--taps", help="tap map CSV: tap_id,junction_id")
    # the task 21 final verification: a help that quotes a default or a limit renders it from the code's constant (or
    # takes the contract's words, which are rendered from them), so it cannot go stale
    ap.add_argument("--dose", type=float, default=None,
                    help=f"free chlorine leaving the plant, mg/L (default {DEFAULT_DOSE_MGL['free_chlorine']:.1f}); with "
                         f"chloramine, total chlorine (default {DEFAULT_DOSE_MGL['chloramine']:.1f})")
    ap.add_argument("--threshold", type=float, default=None,
                    help=f"minimum residual, mg/L: default {DEFAULT_THRESHOLD_MGL['free_chlorine']:.1f} free chlorine, or "
                         f"{DEFAULT_THRESHOLD_MGL['chloramine']:.1f} total chlorine with chloramine (a common utility "
                         f"operating target, not a California rule: California requires a detectable residual)")
    ap.add_argument("--disinfectant", choices=sorted(SPECIES_COLUMN), default="free_chlorine",
                    help="free_chlorine (the default) or chloramine (total chlorine; task 12)")
    ap.add_argument("--ph", type=float, default=None, help="chloramine only: the plant water's logged pH (with --cl2n "
                                                           "it sets the bulk-rate prior, model (c))")
    ap.add_argument("--cl2n", type=float, default=None, help="chloramine only: the plant's chlorine to ammonia-N ratio "
                                                             "by mass (with --ph)")
    ap.add_argument("--window", type=int, default=3, help="months of history to fit on")
    ap.add_argument("--synthetic", metavar="NETWORK", help="demonstrate on a synthetic six-month log from this bundled network")
    ap.add_argument("--plant", help="plant log CSV: month,temp_C[,toc_mgL,dose_mgL]; makes the model temperature-aware")
    ap.add_argument("--out", default=None, help="where the results go: default outputs/private/pilot (git-ignored) for "
                                                "real data, outputs/pilot with --synthetic")
    ap.add_argument("--allow-large", action="store_true",
                    help=f"build the grid of a model with more than {LARGE_JUNCTIONS:,} junctions (task 15b: the "
                         f"readiness report gives its run cost)")
    ap.add_argument("--date-format", default=None,
                    help="how the log's dates are written, when they read two ways (task 16): MDY (month first, "
                         "10/03/2026 is October 3), DMY (day first), YMD, or a strptime format such as "
                         "'%%m/%%d/%%Y %%I:%%M %%p'")
    ap.add_argument("--analyzer-hour", type=float, default=ANALYZER_HOUR,
                    help=CONTRACT_CLI["residualmap.pilot"]["--analyzer-hour"])
    ap.add_argument("--keep-every-reading", action="store_true",
                    help="keep every reading of a tap that looks like an online analyzer (a sampling study)")
    ap.add_argument("--columns", default=None,
                    help="the grab log's columns when the guess is wrong (task 16): 'tap=Sample Point;reading=Total "
                         "Chlorine (mg/L);date=Sample Date/Time;time=' (an empty time: the time is in the date column); "
                         "also analyte=<column> and unit=<column> for a log with one row per analyte, and "
                         "analyte_value=<value> for the analyte whose rows are the reading (a LIMS code such as "
                         "CHLF)")
    ap.add_argument("--inlets", default=None,
                    help="task 16c: an inlet residuals CSV (inlet_node_id,residual_mgL; "
                         "docs/templates/inlet_residuals_template.csv) for a zone whose inlets carry different "
                         "residuals: each inlet is dosed in proportion, and --dose means the highest residual entering "
                         "the zone (default: the file's highest). A tap map's role column can mark a tap as a zone "
                         "'inlet' too: the median of its readings sets that inlet's residual, and none is fitted")
    ap.add_argument("--holdout-sites", type=int, default=None, metavar="K",
                    help="task 16c: also leave K sites out: within each held-out month, refit without every sample at "
                         "K of its sites and predict them, rotating through all its sites in seeded groups of K "
                         "(rows 'held_out_sites': the map's accuracy where there is no tap)")
    a = ap.parse_args()
    if a.holdout_sites is not None and a.holdout_sites < 1:
        ap.error("--holdout-sites needs 1 or more sites per group")
    if a.inlets and a.synthetic:
        ap.error("--inlets is for a zone's real data; --synthetic doses every inlet alike")
    if a.inlets and a.disinfectant == "chloramine":      # the review of task 16c
        ap.error(f"--inlets is refused with --disinfectant chloramine, since {INLETS_NOT_CA}")
    dose_given = a.dose is not None      # task 16c: without --dose, inlet residuals set it to their highest
    if a.threshold is None:
        a.threshold = DEFAULT_THRESHOLD_MGL[a.disinfectant]     # task 21: the defaults the contract quotes
    if a.dose is None:
        a.dose = DEFAULT_DOSE_MGL[a.disinfectant]
        if a.disinfectant == "chloramine":
            print(f"--dose not given: taking {a.dose:g} mg/L total chlorine leaving the plant (typical chloramine doses "
                  f"are 1.5 to 4.0); give the plant's own", file=sys.stderr)
    if (a.ph is not None or a.cl2n is not None) and a.disinfectant != "chloramine":
        ap.error("--ph and --cl2n set the chloramine mode's prior; use them with --disinfectant chloramine")
    if (a.ph is None) != (a.cl2n is None):
        ap.error("give both --ph and --cl2n, or neither")
    if a.synthetic and a.disinfectant != "free_chlorine":
        ap.error("--synthetic makes a free-chlorine log; the chloramine mode's simulated test is "
                 "python -m residualmap.experiment <net> --disinfectant=chloramine")
    a.out = a.out or default_out(bool(a.synthetic))
    os.makedirs(a.out, exist_ok=True)
    inlet_ratios = None      # task 16c: per-inlet dose ratios (None: every inlet at the one dose, the committed path)
    if a.synthetic:
        truth_kw = {"Net2": dict(kb_per_day=0.10, kw_m_per_day=0.20), "ky4": dict(source_dose=2.0)}.get(a.synthetic, {})
        log, inp = synthetic_log(a.synthetic, **truth_kw)
        dose = truth_kw.get("source_dose", a.dose)
        log.to_csv(os.path.join(a.out, f"synthetic_log_{a.synthetic}.csv"), index=False)
    else:
        if not (a.inp and a.samples and a.taps):
            ap.error("--inp, --samples and --taps are required (or --synthetic NETWORK)")
        try:
            inp, report = prepare_pilot_inp(a.inp, a.disinfectant)
        except FileNotFoundError as e:
            sys.exit(str(e))
        # task 17's fix of task 16c's leftovers: the readiness warning for several entry points names the app's way to
        # give inlet residuals; the pilot says its own
        report.warnings = [w.replace(INLET_APP_ROUTE, INLET_PILOT_ROUTE) for w in report.warnings]
        with open(os.path.join(a.out, f"import_report_{inp_stem(a.inp)}{'_plant' if a.plant else ''}.json"), "w") as f:
            f.write(report.to_json() + "\n")
        print(report.to_markdown())
        if inp is None:
            sys.exit(f"{a.inp}: {report.summary()} (the report is in {a.out})")
        if report.run_cost.get("large") and not a.allow_large:
            sys.exit(f"{a.inp}: a large model ({report.run_cost['junctions']} junctions): its grid is "
                     f"{report.run_cost['runs_per_dose']} runs, about {report.run_cost.get('minutes', 0):g} min and "
                     f"{report.run_cost['grid_file_mb']:g} MB. Give --allow-large to build it.")
        grid = report.run_cost.get("grid", "full")
        sc_name = types.SimpleNamespace(wn_name=inp, source_dose=float(a.dose))
        if not os.path.exists(grid_cache_path(sc_name, "outputs/cache", grid, ca_condition() if grid == "chloramine"
                                              else None)):
            ok, msg = size_gate(report.run_cost, "outputs/cache")
            if not ok:
                sys.exit(f"{a.inp}: {msg}")
        cols = None
        if a.columns:
            cols = {}
            for part in a.columns.split(";"):
                k, _, v = part.partition("=")
                if k.strip() not in GRAB_LOG_MAPPING_KEYS:      # task 21: the keys come from the contract
                    ap.error(f"--columns: {k.strip()!r} is not one of {', '.join(GRAB_LOG_MAPPING_KEYS)}")
                cols[k.strip()] = v.strip() or None
        start_h = float(report.flags.get("start_clock_h") or 0.0)
        try:            # task 16: each time of day becomes a model hour with the file's START CLOCKTIME
            g_log = read_log(a.samples, a.taps, a.disinfectant, date_format=a.date_format, start_clock_h=start_h,
                             analyzer_hour=a.analyzer_hour, columns=cols, thin_analyzers=not a.keep_every_reading)
        except ValueError as e:
            sys.exit(f"{a.samples}: {e}")
        print_log(g_log)        # what load_log printed
        log = g_log.frame
        # task 16c: per-inlet residuals, from --inlets and from taps the tap map marks as a zone inlet (the median of
        # their readings over the whole log); without either, nothing below runs and the committed path is unchanged.
        # Its review: every tap marked as a zone inlet must give its inlet a usable reading (inlet_tap_residuals), the
        # file is read for the disinfectant's species and its notes are printed, and not in the chloramine mode
        file_res, entry = None, report.flags.get("entry_points") or []
        if a.inlets:
            try:
                with open(a.inlets, "rb") as fh:
                    file_res, f_info = read_inlet_residuals(fh.read(), os.path.basename(a.inlets),
                                                            disinfectant=a.disinfectant)
            except (OSError, LogRefused) as e:
                sys.exit(f"{a.inlets}: {e}")
            for n_ in f_info["notes"]:
                print(f"note: {a.inlets}: {n_}")
        if g_log.inlet_taps and a.disinfectant == "chloramine":
            sys.exit(f"{a.taps}: the tap map marks {', '.join(sorted(g_log.inlet_taps)[:6])} as a zone inlet, and "
                     f"{INLETS_NOT_CA} Clear the role in the tap map, or leave those readings out of the log.")
        try:
            tap_res = inlet_tap_residuals(g_log, entry, None, file_res)
        except LogRefused as e:
            sys.exit(f"{a.taps}: {e}")
        if file_res or tap_res:
            try:
                inlets = combine_inlet_residuals(file_res, tap_res)
                inlet_ratios, inlet_info = inlet_ratios_from(inlets, entry)
            except LogRefused as e:
                sys.exit(f"{a.inlets or a.taps}: {e}")
            print(f"note: {inlet_info['sentence']} From: "
                  + "; ".join(f"{n}, {d['source']}" for n, d in inlets.items()) + ".")
            if tap_res:
                print("note: an inlet tap's residual is the median of its readings over the whole log, the held-out "
                      "months' included: one set of ratios for every held-out month, as a boundary condition known "
                      "when the map is made (the app takes the median over its window instead).")
            top = inlet_info["highest_mgL"]
            if not dose_given:
                a.dose = top
                print(f"note: --dose not given: taking the highest inlet residual, {top:g} mg/L.")
            elif abs(a.dose - top) > 0.1 * top:
                print(f"warning: --dose {a.dose:g} mg/L is the highest residual entering the zone when inlet residuals "
                      f"are given, and the highest given is {top:g} mg/L: check one of them.")
            print(f"note: {INLET_ZONE_NOTE}")
            if a.plant and inlet_ratios is not None:
                sys.exit(f"{a.inlets or a.taps}: per-inlet residuals with --plant are not built (the temperature model's "
                         f"grids are made at one dose for every inlet). Run without one of them.")
            sc_in = types.SimpleNamespace(wn_name=inp, source_dose=float(a.dose))
            if inlet_ratios is not None and not os.path.exists(grid_cache_path(
                    sc_in, "outputs/cache", grid, ca_condition() if grid == "chloramine" else None, inlet_ratios)):
                ok, msg = size_gate(report.run_cost, "outputs/cache")
                if not ok:
                    sys.exit(f"{a.inp}: {msg}")
        if shows_clock(start_h):    # the review of task 16: the outputs' hours are model hours, and the run says so
            print(f"note: {a.inp} starts at {_clock(start_h)} (START CLOCKTIME): the log's times of day were converted "
                  f"to model hours, and every hour in the pilot's outputs (predictions and validation files) is a "
                  f"model hour; model hour 0 is {clock_label(0, start_h)} and model hour 3 is "
                  f"{clock_label(3, start_h)}.")
        dose = a.dose
    try:            # the task 21 review: a plant log that cannot be read is a plain message, not a traceback
        plant = load_plant(a.plant) if a.plant else None
    except (OSError, ValueError) as e:
        sys.exit(f"{a.plant}: {e}")
    grid = "full" if a.synthetic or a.disinfectant == "chloramine" or plant is not None else grid
    try:            # task 17's fix of task 16c's leftovers: a refusal inside validate is a plain message, not a traceback
        preds, summary = validate(inp, log, dose=dose, threshold=a.threshold, window_months=a.window, plant=plant,
                                  disinfectant=a.disinfectant, ph=a.ph, cl2n=a.cl2n, grid=grid,
                                  inlet_ratios=inlet_ratios, holdout_sites=a.holdout_sites)
    except ValueError as e:
        sys.exit(f"{a.samples or a.synthetic}: {e}")
    tag = (a.synthetic or inp_stem(a.inp)) + ("_plant" if plant is not None else "")
    tag += ("_inlets" if inlet_ratios is not None else "") + (f"_sites{a.holdout_sites}" if a.holdout_sites else "")
    preds.to_csv(os.path.join(a.out, f"predictions_{tag}.csv"), index=False)
    summary.to_csv(os.path.join(a.out, f"validation_{tag}.csv"), index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 40)
    print(f"\n== hold-out validation, {tag}: fit on the previous {a.window} months, predict the held-out month's samples at their own junction and hour")
    print(summary[["held_out_month", "taps", "n_train", "n_test", "rmse", "rmse_persistence", "rmse_network_mean", "coverage50", "coverage80", "coverage90", "coverage95",
                   "band90_width_mgL", "n_true_below", "recall_below", "false_alarms", "map_kb", "map_kw"]
                  + (["temp_C", "kb20", "P_H0"] if plant is not None else [])].round(3).to_string(index=False))
    if a.holdout_sites:     # task 16c
        print(f"held_out_sites: within each held-out month, groups of {a.holdout_sites} sites, each predicted by a fit on "
              f"the same months without any of its samples (a held-out site has no reading of its own in that fit, so "
              f"there is no persistence baseline; map_kb and map_kw are per group, not shown)")
    agg = {}
    for subset, g in summary.groupby("taps"):
        w = g.n_test / g.n_test.sum()
        pers = (g.rmse_persistence * w).sum(min_count=1) if subset == HELD_OUT_SITES else (g.rmse_persistence * w).sum()
        agg[subset] = {"n": int(g.n_test.sum()), "rmse": float((g.rmse * w).sum()), "rmse_persistence": float(pers),
                       "rmse_network_mean": float((g.rmse_network_mean * w).sum()), **{f"coverage{q}": float((g[f"coverage{q}"] * w).sum()) for q in LEVELS}}
        print(f"mean over held-out months, {subset}:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in agg[subset].items()})
    with open(os.path.join(a.out, f"validation_{tag}.json"), "w") as f:
        json.dump({"per_month": summary.to_dict("records"), "mean": agg}, f, indent=2, default=float)


if __name__ == "__main__":
    main()
