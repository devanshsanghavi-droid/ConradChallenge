"""
operator_form.py: the operator form, a short table of what a map of the mains cannot give, which turns a skeleton model
into one that runs (iteration 5, task 19).

A skeleton built from a GIS layer (residualmap/gis_read.py, step 4 of the Get started page) or a CAD drawing
(residualmap/cad_read.py, step 5) has pipes, their sizes and where they meet, and nothing else: no source, no tank, no
pump, no valve, no control and no water use.  The form gives those, one row per item, in a plain table an operator can
fill in a spreadsheet (template_csv: the columns of COLUMNS, one row of each kind), and apply_operator_form builds the
full model from the skeleton and the form.  It also fills the same gaps in an EPANET model (a missing PRV, a boundary
valve left open).

    python -m residualmap.cad_read DRAWING.dxf --form FORM.csv --out completed.inp   # a drawing, completed
    python -m residualmap.gis_read LAYER --form FORM.csv --out completed.inp         # a GIS layer, completed
    python -m residualmap.operator_form MODEL.inp FORM.csv [--out completed.inp]     # an EPANET model's gaps filled

The first two build the skeleton and apply the form to the build, so every pipe keeps the material its label, layer
or material field names, and a point is read in the drawing's or the layer's own coordinates.  An .inp file keeps
neither, not even a skeleton .inp written from a map: on one, a location is a node ID or a point in the file's own
coordinates, and a pipe_c row can only give one C to every pipe (the material 'all').  The saved check
operator_form_command_line runs all three.

Not used by the app or the pilot: task 19 stopped at its pre-registered bar (ky4's exact drawing, joined within 1 ft and
completed by this form, missed the water-age bars and the chlorine 95th percentile, because the join connects three
ky4 dead ends that lie a fraction of a foot from a passing main; outputs/intake/operator_form.json), so the Get started
page does not offer the form yet.

The form (one row per item; the 'kind' column says what the row is; a location is a node ID of the skeleton, or the
point 'x,y' in the coordinates of the drawing or layer it came from, or 'latitude,longitude' for a layer in degrees,
such as one traced in Google Earth; a pipe is a pipe ID or a point on it, away from its ends):
- source: location, head_ft (the hydraulic grade a reservoir or clearwell holds: a skeleton cannot run without it),
  and dose_mgL or residual_mgL (the chlorine leaving the plant, or entering the zone at this inlet, 0 to 10 mg/L;
  with several sources their residuals become task 16c's inlet ratios on the main page).
- tank: location, base_elev_ft, diameter_ft, min_level_ft, max_level_ft, init_level_ft.
- pump: where it is, either location (the node at the pump, where a drawing puts both of its ends) with pipe (the pipe
  on its suction side: the node is split in two and the pump joins them), or pipe (a line drawn for the pump) with
  upstream (the node on its suction side: the pump replaces that line); and what it does: power_hp (a constant-power
  pump, as every pump of the public models ky4 and ky10), or flow_gpm with head_ft (one design point, or three curve
  points with flows rising and heads falling, numbers separated by spaces), or nothing (a design point is derived from
  the production, assumed).  Pumps side by side in one station (duty and standby, drawn as one point): one row each,
  with the same location and pipe, and they run in parallel.  A well and its pump drawn at one point: the source row
  at the well's node, and the pump in place of the first pipe out of it (pipe, and upstream the well's node).
- pump_rule: pump (the pump row's id) and either tank, on_below_ft and off_above_ft (the tank's water level, in feet
  above its bottom, both inside the tank's own range: above its minimum level and below its maximum), or on_hours
  ('always', 'never', or whole clock hours such as '6-22', by the model's clock; minutes such as '6:30-22' are
  refused, since the clock rules of the model WNTR writes take whole hours only).
- closed_valve: pipe (a normally closed valve, such as a zone boundary).  check_valve: pipe and upstream (the end flow
  comes from).  prv and psv: placed as a pump is, with setting_psi and elevation_ft (the valve's own elevation, which
  the node where its setting acts takes, below a PRV and above a PSV: a pressure setting means nothing without it)
  and diameter_in (optional; by default the pipe it replaces).  EPANET joins no PRV or PSV straight to a tank or a
  source, so one is refused there.
- elevation: location and elevation_ft (needed at and below PRVs and PSVs; a skeleton puts every junction at 0 ft):
  one row a node, none at a source's or a tank's node, and none against a PRV's or PSV's own elevation_ft there.
- production: gallons_per_day (the average day's water use of the whole system).  large_user: location and
  gallons_per_day (its average day).  The production less the large users is spread over the junctions by pipe
  length (each junction takes half the length of every pipe that ends at it).
- pattern: multipliers, 24 hourly numbers from midnight, separated by spaces; without one, DEFAULT_PATTERN, labelled
  assumed.  On an EPANET model whose clock does not start at midnight they are shifted, so each clock hour keeps its
  multiplier.  The pattern shapes the form's own water use only: a model's own demands keep their own patterns, so a
  pattern row in a form with no production or large_user row is refused.
- pipe_c: material and hw_c, the Hazen-Williams C you use for that material (from your records or a test), for every
  pipe whose label, layer or material field names that material (topology.material_of: DIP and ductile iron are one
  material, AC another, though the C table puts both in one class); 'plastic' gives every plastic pipe no PVC or
  polyethylene row gives, and 'all' every pipe no other row gives.  A C the layer's own C field gives is kept.
  Without pipe_c rows a skeleton's C is the new-pipe value of its material (topology.C_TABLE), which cannot see a
  pipe's age: see AGE_NOTE.  On an EPANET model whose headloss formula is not Hazen-Williams a pipe's roughness is
  no C, so pipe_c rows are refused there.

A number is one number a cell ('490 500' is refused, not read as 490,500); an elevation or a head outside ELEV_RANGE_FT
is refused as a slip.

What it never does: run a model with no water source (no_source refusal), with a part of the network that no
source feeds, or with junctions whose water use no source's water can reach (task 19's final verification: a part
the form itself cuts off, by a closed_valve row, a check valve, PRV, PSV or pump pointing away from it, or a pump
that never runs, had been completed with ok); run pumps that have no rule (a pump the form adds with no pump_rule row
is a blocking problem: the measured change when pump rules are lost is in the readiness text); guess an elevation for
a pressure valve (a prv or psv row without elevation_ft is refused); or write a model EPANET cannot run (the
completed model's first hour is run before it is written).  A row it cannot use (an unknown kind, a missing or
implausible value, a location off the network, a pipe two rows claim, a row another row overrides or that changes
nothing) is refused with its reason, every row's reason is listed, and no model is built until the form is fixed.
Its messages name no temporary file.

The readiness report's 'why we need this' sentences read their numbers from the fidelity ladder's files
(outputs/intake/fidelity_ky4.json and fidelity_ky10.json, residualmap/intake_exp.py) at run time, as
capabilities.SOURCES does, and always quote ky4 and ky10 side by side: they are simulations on public models, not a
test on a utility's own system.

Uses WNTR (its model objects and its .inp writer), numpy, scipy and networkx (which WNTR uses) only.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass, field

import numpy as np

from . import topology as T
from .importers import ASSUMED, DERIVED, FROM_FILE, MISSING, LogRefused, read_table

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KINDS = ("source", "tank", "pump", "pump_rule", "closed_valve", "check_valve", "prv", "psv", "elevation",
         "production", "large_user", "pattern", "pipe_c")
COLUMNS = ("kind", "id", "location", "pipe", "upstream", "head_ft", "dose_mgL", "residual_mgL", "base_elev_ft",
           "diameter_ft", "min_level_ft", "max_level_ft", "init_level_ft", "power_hp", "flow_gpm", "pump", "tank",
           "on_below_ft", "off_above_ft", "on_hours", "setting_psi", "elevation_ft", "diameter_in", "gallons_per_day",
           "multipliers", "material", "hw_c", "note")
ALIASES = {"type": "kind", "item": "kind", "name": "id", "node": "location", "at": "location", "node_id": "location",
           "pipe_id": "pipe", "suction_pipe": "pipe", "from": "upstream", "upstream_node": "upstream",
           "residual_entering_mgl": "residual_mgL", "residual_mgl": "residual_mgL", "dose_mgl": "dose_mgL",
           "gpd": "gallons_per_day", "gallons_per_day_gpd": "gallons_per_day", "c": "hw_c", "c_value": "hw_c",
           "hazen_williams_c": "hw_c", "notes": "note", "comment": "note", "base_elevation_ft": "base_elev_ft",
           "minimum_level_ft": "min_level_ft", "maximum_level_ft": "max_level_ft", "initial_level_ft": "init_level_ft",
           "setting_psig": "setting_psi", "flow_gpm_points": "flow_gpm", "head_ft_points": "head_ft"}
KIND_ALIASES = {"reservoir": "source", "clearwell": "source", "well": "source", "plant": "source",
                "closed valve": "closed_valve", "boundary_valve": "closed_valve", "check valve": "check_valve",
                "cv": "check_valve", "pump rule": "pump_rule", "large user": "large_user", "demand": "large_user",
                "pipe c": "pipe_c", "roughness": "pipe_c"}
SNAP_FT = 50.0                 # a location given as a point must lie within this of a node (or of a pipe, for a pipe)
PIPE_CLEAR_FT = 3.0            # ... and, for a pipe, this much nearer to it than to any other pipe (the review: a
                               # point at a tee is as near to three pipes, and took whichever came first)
FAR_UNITS = 1e9                # a point more than this (feet, or model units) from the model's extent is refused as far
                               # off the network, before any distance is measured (task 19's final verification: a point
                               # at 1e160 reached an internal IndexError, and one at 1e150 quoted a 151-digit distance)
ELEV_RANGE_FT = (-1500.0, 20000.0)   # an elevation or a head outside this is refused as a slip: the lowest ground on
                               # Earth is about -1,400 ft and no town lies above about 17,000 ft (646000 for 646 passed)
HEADER_MAX = 120               # a header longer than this names no column the form reads (it is noted, cut short)
SPLIT_SUFFIX = "_in"           # the suction-side node a split adds: <node>_in
FT = 0.3048
GPM_M3S = 6.30901964e-05       # one US gallon per minute in cubic metres per second (WNTR's own conversion)
GPD_PER_GPM = 1440.0
C_RANGE = (40.0, 200.0)        # a Hazen-Williams C outside this is refused (a unit or a typing slip)
HEADLOSS_WORDS = {"D-W": "Darcy-Weisbach (D-W)", "C-M": "Chezy-Manning (C-M)"}
ROUGHNESS_WORDS = {"D-W": "a roughness height", "C-M": "Manning's n"}
DOSE_RANGE = (0.0, 10.0)       # a dose or residual outside this (mg/L) is refused: a unit slip, or a decimal comma read
                               # as a thousands separator ('1,500' is 1500), as the grab-log reader lists one above 10
DERIVED_HEAD_MARGIN = 1.1      # a derived pump design head: the lift from its source to the highest tank's top, +10%
DERIVED_HEAD_MIN_FT = 10.0
FIRST_HOUR_S = 3600            # the completed model's test run before it is written (hydraulics only)
# The default 24-hour demand pattern, used only when the form gives none (labelled assumed): the diurnal pattern both
# public models ky4 and ky10 carry (the KYPIPE research database of water distribution models, Jolly et al. 2014,
# bundled with WNTR): low at night, a morning rise, an evening peak at 20:00.  No published residential table was
# checked for this build, so it is a public model's own pattern, not a standard.
DEFAULT_PATTERN = (0.33, 0.25, 0.209, 0.209, 0.259, 0.36, 0.529, 0.91, 1.2, 1.299, 1.34, 1.34, 1.32, 1.269, 1.25, 1.25,
                   1.279, 1.37, 1.519, 1.7, 1.75, 1.669, 0.899, 0.479)
DEFAULT_PATTERN_SOURCE = ("the diurnal pattern of the public models ky4 and ky10 (the KYPIPE research database, Jolly "
                          "et al. 2014), not your system's")
PATTERN_NAME = "form_demand"
WNTR_HEADER = ("; Filename:", "; WNTR:", "; Created:")    # the comment lines WNTR's .inp writer starts a file with
NO_SOURCE = ("No water source: the form has no source row, so the model has nowhere for water (and chlorine) to enter "
             "and cannot run. Add a source row: the node at the plant, well or zone inlet, and the hydraulic grade it "
             "holds (head_ft).")
NO_USE = ("No water use: the form gives neither the production (gallons per day) nor any large user, so no water "
          "would move. Add a production row with your average day's production.")
AGE_NOTE = ("{n} of the skeleton take{s} Hazen-Williams C from {their} material alone (new-pipe values, or 130 where the "
            "map names none), so a decades-old cast iron main gets more than a new ductile iron one, and the model's "
            "old-pipe factor, which scales wall decay by C, cannot see the pipe's age. Give your own C for each "
            "material (pipe_c rows) if you have them.")
NO_MATERIALS = ("this model carries no pipe materials (an EPANET file keeps none, not even one written from a drawing "
                "or a GIS layer): use the material 'all' for one C on every pipe, or apply the form to the drawing or "
                "the layer itself, which keeps each pipe's material")


def age_note(n_table: int, n_kept: int = 0) -> str:
    """AGE_NOTE for the pipes whose C is still the table's (none of the form's pipe_c rows gave it, and no C field of
    the layer), with the pipes whose C the layer's own field gave, which keep it."""
    s = AGE_NOTE.format(n=f"{n_table:,} pipe{'s' if n_table != 1 else ''}", s="s" if n_table == 1 else "",
                        their="its" if n_table == 1 else "their")
    if n_kept:
        s += (f" The {n_kept:,} pipe{'s' if n_kept != 1 else ''} whose C the layer's own C field gives "
              f"keep{'s' if n_kept == 1 else ''} it.")
    return s


# The template's rows, one of each kind an operator is most likely to need, with a plain note each (template_csv).  Not
# shipped in docs/templates/: the app does not take the form (task 19 stopped at its bar before the form was wired in),
# so the template is a fixture (fixtures.py) and the command line's example only.
TEMPLATE_ROWS = (
    {"kind": "source", "id": "PLANT", "location": "1000,2000", "head_ft": "490", "dose_mgL": "1.2",
     "note": "where water enters: the node at the plant or well (its ID, or x,y on your map), the hydraulic grade it "
             "holds in feet, and the chlorine leaving it"},
    {"kind": "tank", "id": "TANK-1", "location": "5200,4100", "base_elev_ft": "646", "diameter_ft": "58",
     "min_level_ft": "5", "max_level_ft": "30", "init_level_ft": "20",
     "note": "levels are feet of water above the tank's bottom"},
    {"kind": "pump", "id": "PUMP-1", "location": "1100,2050", "pipe": "1150,2040", "power_hp": "50",
     "note": "a pump drawn as one point: the node there, and a point on the pipe that brings water to it; or flow_gpm "
             "and head_ft for its design point; pumps side by side: one row each, the same node and pipe"},
    {"kind": "pump_rule", "pump": "PUMP-1", "tank": "TANK-1", "on_below_ft": "12", "off_above_ft": "27",
     "note": "the pump starts when the tank falls below on_below_ft and stops above off_above_ft; or on_hours: always, "
             "never, or hours such as 6-22"},
    {"kind": "closed_valve", "pipe": "3000,3500", "note": "a valve kept closed, such as a zone boundary: a point on its "
                                                          "pipe"},
    {"kind": "prv", "location": "4000,1500", "pipe": "3950,1520", "setting_psi": "60", "elevation_ft": "520",
     "note": "a pressure-reducing valve: its node and a point on the pipe on its upstream side, its setting and its "
             "elevation"},
    {"kind": "elevation", "location": "4100,1400", "elevation_ft": "515",
     "note": "a junction's ground elevation, needed at and below PRVs"},
    {"kind": "production", "gallons_per_day": "1500000",
     "note": "your average day's production; what large users do not use is spread over the mains by length"},
    {"kind": "large_user", "location": "2500,2600", "gallons_per_day": "80000",
     "note": "a school, a plant or a hospital, at its node"},
    {"kind": "pattern", "multipliers": " ".join(f"{x:g}" for x in (0.33, 0.25, 0.21, 0.21, 0.26, 0.36, 0.53, 0.91, 1.2,
                                                                  1.3, 1.34, 1.34, 1.32, 1.27, 1.25, 1.25, 1.28, 1.37,
                                                                  1.52, 1.7, 1.75, 1.67, 0.9, 0.48)),
     "note": "optional: 24 hourly multipliers from midnight; without this row a public model's pattern is used, "
             "marked assumed"},
    {"kind": "pipe_c", "material": "PVC", "hw_c": "150",
     "note": "optional: the Hazen-Williams C you use for a material in your map"})


def form_csv(rows) -> bytes:
    """A form as CSV, the columns of COLUMNS."""
    import csv
    import io
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(COLUMNS)
    for r in rows:
        w.writerow([r.get(c, "") for c in COLUMNS])
    return buf.getvalue().encode("utf-8")


def template_csv() -> bytes:
    return form_csv(TEMPLATE_ROWS)


class FormRefused(ValueError):
    """A form that cannot be read as a table at all, in plain words."""


# ----------------------------------------------------------------------------- reading the form
@dataclass
class FormRow:
    line: int                    # its line in the file (or row of the sheet)
    kind: str
    cells: dict                  # canonical column -> text (stripped), only the filled ones


@dataclass
class OperatorForm:
    name: str
    rows: list = field(default_factory=list)          # FormRow
    refused: list = field(default_factory=list)       # {'line', 'kind', 'reason'}
    notes: list = field(default_factory=list)
    header: list = field(default_factory=list)

    def of(self, kind: str) -> list:
        return [r for r in self.rows if r.kind == kind]


def _canonical(h: str) -> str | None:
    """A header as one of COLUMNS: case, spaces and dashes do not matter, and a unit in brackets is the column's unit
    ('Dose (mg/L)' is dose_mgL, 'Head [ft]' head_ft, 'Base elevation (ft)' base_elev_ft; the review: the grab-log
    reader dropped the brackets, this reader dropped the column).  Linear in the header (task 19's final verification:
    the bracket pattern had three runs of spaces competing, so 'Dose (' and 1,600 spaces took 3.6 s and a 6.5 KB
    form over 100 s): a bracket's contents hold no bracket, nothing else competes for a space, and a header longer than
    HEADER_MAX names no column."""
    s = str(h or "").replace("﻿", "").strip()
    if len(s) > HEADER_MAX:
        return None
    s = re.sub(r"[(\[]([^()\[\]]*)[)\]]", lambda m: " " + re.sub(r"[\s/.]+", "", m.group(1)) + " ", s)
    s = re.sub(r"[\s\-/]+", "_", s.strip()).strip("_")
    low = s.lower()
    for c in COLUMNS:
        if low == c.lower():
            return c
    return ALIASES.get(low)


def _short(v, n: int = 40) -> str:
    """A cell or header as a message quotes it: at most n characters, its spaces' runs shown as one."""
    s = re.sub(r"\s+", " ", str(v)).strip()
    return s if len(s) <= n else s[:n] + "..."


def _cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return "" if math.isnan(v) else repr(v)
    return str(v).strip()


def number(text: str) -> float | None:
    """A number written in a cell ('1,200' and '1.5' are read; '1,5' too, as a decimal comma).  One number a cell: a
    space inside it is not read as a thousands separator (task 19's final verification: a head written '490 500' was
    read as 490,500 ft)."""
    s = str(text).strip()
    if "_" in s:                   # Python's float reads '1_500' as 1500: one number, as written, only
        return None
    if re.fullmatch(r"[-+]?\d{1,3}(?:,\d{3})+(?:\.\d+)?", s):
        s = s.replace(",", "")
    elif re.fullmatch(r"[-+]?\d+,\d+", s):
        s = s.replace(",", ".")
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def numbers(text: str) -> list[float] | None:
    """Several numbers in one cell, separated by spaces or semicolons (commas too when no decimal comma is used)."""
    parts = [p for p in re.split(r"[\s;]+" if re.search(r"\d,\d", text or "") and " " in (text or "") else r"[\s;,]+",
                                 str(text or "").strip()) if p]
    vals = [number(p) for p in parts]
    return None if not vals or any(v is None for v in vals) else vals


# the values each kind needs (any of the tuples in a list satisfies it), and the numbers it may carry
NEEDS = {"source": [("location", "head_ft")], "tank": [("location", "base_elev_ft", "diameter_ft", "min_level_ft",
                                                        "max_level_ft", "init_level_ft")],
         "pump": [("location", "pipe"), ("pipe", "upstream")],
         "pump_rule": [("pump", "tank", "on_below_ft", "off_above_ft"), ("pump", "on_hours")],
         "closed_valve": [("pipe",)], "check_valve": [("pipe", "upstream")],
         "prv": [("location", "pipe", "setting_psi", "elevation_ft"), ("pipe", "upstream", "setting_psi", "elevation_ft")],
         "psv": [("location", "pipe", "setting_psi", "elevation_ft"), ("pipe", "upstream", "setting_psi", "elevation_ft")],
         "elevation": [("location", "elevation_ft")], "production": [("gallons_per_day",)],
         "large_user": [("location", "gallons_per_day")], "pattern": [("multipliers",)], "pipe_c": [("material", "hw_c")]}
NUMERIC = ("head_ft", "dose_mgL", "residual_mgL", "base_elev_ft", "diameter_ft", "min_level_ft", "max_level_ft",
           "init_level_ft", "power_hp", "on_below_ft", "off_above_ft", "setting_psi", "elevation_ft", "diameter_in",
           "gallons_per_day", "hw_c")
WORDS = {"head_ft": "head_ft (the hydraulic grade it holds, in feet)", "location": "location (a node ID, or x,y)",
         "pipe": "pipe (a pipe ID, or a point on it as x,y)", "upstream": "upstream (the node on the side water comes from)",
         "setting_psi": "setting_psi", "elevation_ft": "elevation_ft (the valve's elevation)",
         "gallons_per_day": "gallons_per_day", "multipliers": "multipliers (24 hourly numbers)",
         "material": "material", "hw_c": "hw_c", "pump": "pump (the pump row's id)", "tank": "tank",
         "on_below_ft": "on_below_ft", "off_above_ft": "off_above_ft", "on_hours": "on_hours"}
KIND_WORDS = {("elevation", "elevation_ft"): "elevation_ft (the ground elevation there, in feet)"}
# the grab-log reader's refusals say 'log' (importers.read_table); the form's say 'form' (the review)
FORM_WORDS = (("A grab log or tap map is a table", "The operator form is a table"),
              ("a log is read up to", "a form is read up to"),
              ("Send the last few months only.", "An operator form has one row per item."),
              ("at least one reading are needed", "at least one row are needed"))


def _missing_words(kind: str, cells: dict) -> str | None:
    """None when the row has what its kind needs; else what is missing, in plain words."""
    best = None
    for need in NEEDS[kind]:
        miss = [c for c in need if not cells.get(c)]
        if not miss:
            return None
        if best is None or len(miss) < len(best):
            best = miss
    alts = NEEDS[kind]
    say = ", ".join(KIND_WORDS.get((kind, c), WORDS.get(c, c)) for c in best)
    if len(alts) > 1 and kind in ("pump", "prv", "psv") and set(best) & {"location", "pipe", "upstream"}:
        say += (" (a node to split with the pipe on its upstream side, or the pipe it replaces with its upstream "
                "node)")
    elif len(alts) > 1 and kind == "pump_rule":
        say = ("a tank with on_below_ft and off_above_ft, or on_hours ('always', 'never' or hours such as '6-22')"
               + ("" if cells.get("pump") else ", and the pump (the pump row's id)"))
    return f"missing {say}"


def form_advice(s) -> str:
    """What to do when the file given as the form is not a table, in the form's own words (task 19's final
    verification: the sniff's advice is the main page's, for a model upload: an EPANET model given as the form was told
    'It is checked below and then run', a drawing or a GIS layer 'send its EPANET export instead')."""
    if s.category in ("model", "other_model", "results"):
        return ("It is a model, not the form: the form is a table, one row per item (a source, a tank, a pump and so "
                "on), that completes a model. On the command line the model comes first and the form second.")
    if s.category in ("drawing", "gis"):
        return ("It is a map of the mains, not the form: the form is a table, one row per item (a source, a tank, a "
                "pump and so on), that completes the model built from the map (on the command line, cad_read or "
                "gis_read with --form, the map first and the form after --form).")
    return ("Fill in the operator form's template, one row per item (a source, a tank, a pump and so on), in a "
            "spreadsheet and save it as CSV or as an Excel workbook (.xlsx).")


def read_operator_form(data: bytes, name: str = "operator_form.csv", max_unpacked: int | None = None) -> OperatorForm:
    """The form as rows: CSV in any encoding or separator the grab-log reader takes, or an Excel workbook.  Rows of an
    unknown kind, rows missing a value their kind needs and numbers that do not read are refused with their reason
    (OperatorForm.refused); locations are checked against the skeleton later (apply_operator_form).  Raises
    FormRefused for a file that is no table (with form_advice, sniffed before the grab-log reader gives its own), or
    a table with no 'kind' column."""
    from .importers import LOG_MAX_BYTES, sniff
    if len(data) <= LOG_MAX_BYTES:     # a larger file is refused by the reader's own size limit, unsniffed
        s = sniff(data, name)
        if not (s.category == "table" or s.kind in ("text_other", "empty")):     # what read_table refuses as no table
            raise FormRefused(f"{s.sentence()} The operator form is a table: CSV, or an Excel workbook (.xlsx). "
                              f"{form_advice(s)}")
    try:
        t = read_table(data, name, max_unpacked=max_unpacked)
    except LogRefused as e:
        msg = str(e)
        for old, new in FORM_WORDS:
            msg = msg.replace(old, new)
        raise FormRefused(msg) from None
    canon = [_canonical(h) for h in t.header]
    if "kind" not in canon:
        raise FormRefused("The form has no 'kind' column: each row needs one (source, tank, pump, pump_rule, "
                          "closed_valve, check_valve, prv, psv, elevation, production, large_user, pattern or pipe_c). "
                          "Start from the operator form template.")
    form = OperatorForm(name=name, header=list(t.header))
    unknown = [h for h, c in zip(t.header, canon) if c is None and str(h).strip()]
    if unknown:
        form.notes.append(f"Columns not read: {', '.join(_short(h) for h in unknown[:8])}.")
    for line, r in zip(t.row_numbers, t.rows):
        cells = {}
        for c, v in zip(canon, r):
            if c and c != "note" and _cell(v) != "" and c not in cells:
                cells[c] = _cell(v)
        raw_kind = cells.pop("kind", "")
        k = re.sub(r"\s+", " ", raw_kind.strip().lower())
        k = KIND_ALIASES.get(k, k.replace(" ", "_"))
        if not raw_kind:
            if cells:
                form.refused.append({"line": line, "kind": "", "reason": "no kind given"})
            continue
        if k not in KINDS:
            form.refused.append({"line": line, "kind": raw_kind,
                                 "reason": f"unknown kind {raw_kind!r} (the kinds are {', '.join(KINDS)})"})
            continue
        why = _missing_words(k, cells)
        bad = [c for c in NUMERIC if c in cells and number(cells[c]) is None]
        if k in ("pump", "pattern") and cells.get("flow_gpm") and numbers(cells["flow_gpm"]) is None:
            bad.append("flow_gpm")
        if k == "pump" and cells.get("head_ft") and numbers(cells["head_ft"]) is None:
            bad = [b for b in bad if b != "head_ft"] + ["head_ft"]
        elif k == "pump" and "head_ft" in bad and numbers(cells["head_ft"]) is not None:
            bad.remove("head_ft")
        if why:
            form.refused.append({"line": line, "kind": k, "reason": why})
            continue
        if bad:
            form.refused.append({"line": line, "kind": k,
                                 "reason": f"not a number: {', '.join(f'{b} {cells[b]!r}' for b in bad)}"})
            continue
        form.rows.append(FormRow(line=line, kind=k, cells=cells))
    if not form.rows and not form.refused:
        raise FormRefused("The form has no rows: fill in at least a source and the production.")
    return form


# ----------------------------------------------------------------------------- the skeleton
@dataclass
class Skeleton:
    """What apply_operator_form completes: the model's .inp text, and, for a skeleton from a GIS layer or a drawing,
    each pipe's material (topology.material_of's, or None), the pipes whose C the layer's own C field gave (c_kept,
    which pipe_c rows never change), and how a point in the layer's or drawing's own coordinates maps to the
    skeleton's (feet): frame {'kind': 'projected', 'unit_m'} or {'kind': 'enu', 'lat0', 'lon0'} (a layer in degrees).
    For an EPANET model the frame and the materials are None: points are in the model's own coordinates."""
    inp_text: str
    name: str = "skeleton.inp"
    source: str = "inp"                 # 'gis', 'cad' or 'inp'
    material: dict | None = None        # pipe ID -> material or None
    c_kept: frozenset = frozenset()     # pipe IDs whose C came from the layer's own C field
    frame: dict | None = None


def skeleton_from_build(b, source: str) -> Skeleton:
    """A Skeleton from a GisBuild (source 'gis') or a CadBuild ('cad')."""
    fr = dict(getattr(b, "frame", {}) or {})
    frame = None
    if fr.get("kind") == "enu":
        frame = {"kind": "enu", "lat0": float(fr["lat0"]), "lon0": float(fr["lon0"])}
    elif fr.get("kind") == "projected" and fr.get("unit_m"):
        frame = {"kind": "projected", "unit_m": float(fr["unit_m"])}
    return Skeleton(inp_text=b.inp_text, name=getattr(b, "name", "skeleton"), source=source,
                    material=dict(getattr(b, "c_material", {}) or {}),
                    c_kept=frozenset(getattr(b, "c_from_field", ()) or ()), frame=frame)


def read_model_text(data: bytes) -> str:
    """A model file's text, decoded as task 15's intake decodes one (importers.decode_text: a byte-order mark, UTF-16,
    a Windows code page)."""
    from .importers import decode_text
    return decode_text(data)[0]


def _load(text: str):
    import wntr
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        p = os.path.join(d, "m.inp")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        return wntr.network.WaterNetworkModel(p)


def write_text(wn) -> str:
    import wntr
    with tempfile.TemporaryDirectory(prefix="rm_form_") as d:
        p = os.path.join(d, "m.inp")
        wntr.network.write_inpfile(wn, p, units="GPM")
        with open(p, encoding="utf-8") as fh:
            # WNTR's header comments name the file it read (a temporary path) and the time it wrote: left out, so the
            # same model always writes the same text
            return "".join(ln for ln in fh if not ln.startswith(WNTR_HEADER))


def error_words(e: BaseException, n: int = 300) -> str:
    """An error's message as an operator is shown it: the specific error WNTR's reader wraps, when there is one, with
    no temporary file named (task 19's final verification: messages quoted '/var/folders/.../rm_form_*/m.inp') and no
    '%s' left over, at most n characters (cut first, so the patterns below see a short text)."""
    cause = e.__cause__ if isinstance(e.__cause__, Exception) and str(e.__cause__).strip() else None
    s = str(cause if cause is not None else e)[:1000]
    s = re.sub(r"\s*in input file\s*(?:'[^'\n]*'|%s)?", " in the model", s)
    s = re.sub(r"'(?:/|[A-Za-z]:[\\/])[^'\n]*'", "a temporary file", s)
    s = s.replace(" (%s)", "").replace("%s", "").strip()
    return (s or type(e).__name__)[:n]


def first_hour_error(text: str) -> str | None:
    """None when EPANET runs the model's first hour (hydraulics only, in a temporary folder); else a plain sentence
    with EPANET's error (the review: values that pass every row's check gave a model EPANET refuses), naming no
    temporary file (error_words)."""
    import wntr
    with tempfile.TemporaryDirectory(prefix="rm_form_run_") as d:
        p = os.path.join(d, "m.inp")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text)
        try:
            w = wntr.network.WaterNetworkModel(p)
            w.options.time.duration = FIRST_HOUR_S
            w.options.time.report_start = 0
            w.options.quality.parameter = "NONE"
            wntr.sim.EpanetSimulator(w).run_sim(file_prefix=os.path.join(d, "run"))
        except Exception as e:  # noqa: BLE001  (EPANET's own refusal, whatever its kind)
            msg = error_words(e)
            return (f"EPANET cannot run the completed model: its first hour stops with {msg}. Check the form's pumps "
                    f"(their curves and powers), valves (their sizes and settings; EPANET joins no pressure valve "
                    f"straight to a tank, nor two in series or side by side) and tanks (their levels), and that water "
                    f"can reach every junction.")
    return None


def _looks_like_degrees(P) -> bool:
    """An EPANET model's coordinates that read as longitude and latitude of one place (every x within 180, every y
    within 90, all within 2 degrees): the snap distance is then taken in degrees."""
    if not len(P):
        return False
    x, y = P[:, 0], P[:, 1]
    return bool(np.all(np.abs(x) <= 180) and np.all(np.abs(y) <= 90)
                and math.hypot(float(x.max() - x.min()), float(y.max() - y.min())) < 2.0)


def degrees_to_feet(P):
    """For an EPANET model whose coordinates read as degrees (x longitude, y latitude): a function putting such points
    on the tangent plane at the nodes' median place, in feet east and north (importers.enu_from_latlon).  Task 19's
    final verification: the snap had been 50 / 364,000 degrees in every direction, one degree of latitude, so about
    40 ft east and west at latitude 38, and a refusal said 57 ft for a point 45 ft east of a junction."""
    from .importers import enu_from_latlon, plane_centre
    lat0, lon0 = plane_centre(P[:, 1], P[:, 0])

    def to_ft(Q):
        Q = np.asarray(Q, float).reshape(-1, 2)
        return enu_from_latlon(Q[:, 1], Q[:, 0], lat0, lon0) / FT
    return to_ft


class _Geometry:
    """Nearest node and nearest pipes of a point, in the model's coordinates, or, with to_ft (degrees_to_feet), in feet
    on a tangent plane, the point given in the model's own coordinates either way."""

    def __init__(self, wn, to_ft=None):
        from scipy.spatial import cKDTree
        self.to_ft = to_ft
        self.nodes = list(wn.node_name_list)
        P = np.array([wn.get_node(n).coordinates for n in self.nodes], float).reshape(-1, 2)
        self.P = to_ft(P) if to_ft is not None and len(P) else P
        self.tree = cKDTree(self.P) if len(self.P) else None
        segs, owner = [], []
        for pn in wn.pipe_name_list:
            p = wn.get_link(pn)
            pts = [p.start_node.coordinates] + list(p.vertices) + [p.end_node.coordinates]
            for a, b in zip(pts, pts[1:]):
                segs.append((a[0], a[1], b[0], b[1]))
                owner.append(pn)
        S = np.array(segs, float).reshape(-1, 4)
        if to_ft is not None and len(S):
            S = np.hstack([to_ft(S[:, :2]), to_ft(S[:, 2:])])
        self.segs = S
        self.owner = owner
        self.lo, self.hi = (self.P.min(0), self.P.max(0)) if len(self.P) else (np.zeros(2), np.zeros(2))

    def q(self, xy) -> np.ndarray:
        """A point in the units distances are measured in (feet on the tangent plane for a model in degrees)."""
        p = np.asarray(xy, float).reshape(1, 2)
        return (self.to_ft(p) if self.to_ft is not None else p)[0]

    def beyond(self, xy) -> bool:
        """The point is not finite, or more than FAR_UNITS from the box of the model's nodes: far off the network,
        before any distance is measured (task 19's final verification: at 1e160 the KD-tree found no neighbour and
        returned an index one past the last node)."""
        with np.errstate(all="ignore"):
            p = self.q(xy)
            if not np.all(np.isfinite(p)) or self.tree is None:
                return True
            gap = np.maximum(np.maximum(self.lo - p, 0.0), p - self.hi)
            return not float(np.hypot(*gap)) <= FAR_UNITS

    def node(self, xy) -> tuple[str | None, float]:
        if self.tree is None or self.beyond(xy):
            return None, math.inf
        d, i = self.tree.query(self.q(xy))
        if not math.isfinite(float(d)) or int(i) >= len(self.nodes):
            return None, math.inf
        return self.nodes[int(i)], float(d)

    def pipes(self, xy) -> tuple[tuple[str | None, float], tuple[str, float] | None]:
        """(the nearest pipe and its distance, the nearest other pipe and its distance or None)."""
        if not len(self.segs) or self.beyond(xy):
            return (None, math.inf), None
        xy = self.q(xy)
        A, B = self.segs[:, :2], self.segs[:, 2:]
        AB = B - A
        L2 = (AB ** 2).sum(1)
        t = np.clip(((np.asarray(xy) - A) * AB).sum(1) / np.where(L2 > 0, L2, 1), 0, 1)
        d = np.hypot(*(A + AB * t[:, None] - np.asarray(xy)).T)
        order = np.argsort(d, kind="stable")
        first = self.owner[int(order[0])]
        other = next((int(i) for i in order[1:] if self.owner[int(i)] != first), None)
        return (first, float(d[order[0]])), (None if other is None else (self.owner[other], float(d[other])))


def _point(text: str, frame: dict | None) -> tuple[float, float] | None:
    """'x,y' (or 'x y', 'x;y') as a point in the skeleton's coordinates; for a layer in degrees, 'latitude,longitude'
    (as Google Earth copies a place) is put on the skeleton's tangent plane."""
    vals = re.split(r"\s*[,;]\s*|\s+", str(text).strip())
    if len(vals) != 2:
        return None
    a, b = (number(v) for v in vals)
    if a is None or b is None:
        return None
    if frame and frame.get("kind") == "enu":
        if abs(a) <= 90 and abs(b) <= 180:
            from .importers import enu_from_latlon
            e, n = enu_from_latlon(np.array([a]), np.array([b]), frame["lat0"], frame["lon0"])[0]
            return float(e) / FT, float(n) / FT
        return None
    if frame and frame.get("kind") == "projected":
        f = frame["unit_m"] / FT
        a, b = a * f, b * f
    return (a, b) if math.isfinite(a) and math.isfinite(b) else None


# ----------------------------------------------------------------------------- the result
@dataclass
class FormResult:
    ok: bool = False
    inp_text: str | None = None
    wn: object = None
    refused: list = field(default_factory=list)      # {'line', 'kind', 'reason'}: every row that could not be used
    blocking: list = field(default_factory=list)     # plain sentences: what stops the model from running
    rows: list = field(default_factory=list)         # readiness rows: (item, value, source)
    warnings: list = field(default_factory=list)
    notes: list = field(default_factory=list)        # the form reader's notes (columns not read) and the model's
    inlet_residuals: dict | None = None              # {source node: mg/L} with two or more sources (task 16c's inlets)
    dose_mgL: float | None = None                    # one source's dose, or the highest residual given
    placed: dict = field(default_factory=dict)       # form line -> the node or link it was put at (for checks)
    counts: dict = field(default_factory=dict)

    def refused_lines(self) -> list[str]:
        return [f"line {r['line']} ({r['kind'] or 'no kind'}): {r['reason']}" for r in self.refused]


def _safe_id(want: str) -> str:
    """An ID EPANET reads: the form's own when it has no space, ';' or quote, does not start with '[' and fits in 31
    bytes (so '~@Pump-1' stays itself), else topology.clean_id's ('[END]' becomes END).  A data line that starts with
    '[' is a section header to WNTR's reader and EPANET's (task 19's final verification: a pump id '[END]' dropped the
    rest of the written model, and the operator was given an internal error)."""
    s = str(want).strip()
    if s and not re.search(r'[\s;"]', s) and not s.startswith("[") and len(s.encode("utf-8")) <= T.MAX_ID_BYTES:
        return s
    return T.clean_id(s)[0]


def _unique_name(wn, want: str) -> str:
    base = _safe_id(want)
    names = set(wn.node_name_list) | set(wn.link_name_list)
    if base not in names:
        return base
    k = 2
    while f"{base[:27]}_{k}" in names:
        k += 1
    return f"{base[:27]}_{k}"


def _unique_pattern(wn, want: str) -> str:
    """A pattern name the model does not use yet (the review: a form applied to a model the form had completed removed
    the pattern its junctions use)."""
    names = set(wn.pattern_name_list)
    if want not in names:
        return want
    k = 2
    while f"{want}_{k}" in names:
        k += 1
    return f"{want}_{k}"


def _to_si(value: float, param: str) -> float:
    from wntr.epanet.util import FlowUnits, HydParam, to_si
    return float(to_si(FlowUnits.GPM, float(value), getattr(HydParam, param)))


def length_shares(wn) -> dict:
    """Each junction's share of the pipe length: half of every pipe that ends at it, over all such halves (a pipe's
    half at a tank or a reservoir goes to no junction).  Demand by pipe length, in the form and in the fidelity
    ladder's rung 2."""
    half = {j: 0.0 for j in wn.junction_name_list}
    for _, p in wn.pipes():
        for n in (p.start_node_name, p.end_node_name):
            if n in half:
                half[n] += p.length / 2.0
    tot = sum(half.values())
    return {j: v / tot for j, v in half.items()} if tot > 0 else {}


def _convert_node(wn, old: str, make) -> str:
    """Replace junction old by the node make(wn) adds (a tank or reservoir, at old's place): every link moved onto it,
    then old removed.  Returns the new node's name."""
    j = wn.get_node(old)
    new = make(j)
    nn = wn.get_node(new)
    for ln in list(wn.get_links_for_node(old)):
        lk = wn.get_link(ln)
        if lk.start_node_name == old:
            lk.start_node = nn
        if lk.end_node_name == old:
            lk.end_node = nn
    wn.remove_node(old)
    return new


def _split(wn, node: str, pipe: str) -> str:
    """Split node in two at its place: pipe (which ends at node) moves to a new junction <node>_in; returns it."""
    j = wn.get_node(node)
    new = _unique_name(wn, f"{node[:31 - len(SPLIT_SUFFIX)]}{SPLIT_SUFFIX}")
    wn.add_junction(new, base_demand=0.0, elevation=j.elevation, coordinates=tuple(j.coordinates))
    lk = wn.get_link(pipe)
    nn = wn.get_node(new)
    if lk.start_node_name == node:
        lk.start_node = nn
    else:
        lk.end_node = nn
    return new


MINUTES = "minutes"


def _hours(text: str) -> tuple | str | None:
    """'always', 'never', or (on, off) whole clock hours from '6-22', '6 to 22', '06:00-22:00', '22-6' (overnight);
    MINUTES for hours with minutes ('6:30-22'), which the form refuses: WNTR writes a clock rule's time as a decimal
    hour ('AT CLOCKTIME 6.5') that its own reader then refuses (task 19's final verification, where such a form was
    blocked with a message that blamed its pumps, valves and tanks)."""
    s = str(text).strip().lower()
    if s in ("always", "all day", "24/7", "on", "0-24", "0 to 24", "continuous"):
        return "always"
    if s in ("never", "off", "standby", "stand-by", "spare"):
        return "never"
    m = re.fullmatch(r"(\d{1,2})(?::(\d{2}))?\s*(?:-|to)\s*(\d{1,2})(?::(\d{2}))?", s)
    if not m:
        return None
    if int(m.group(2) or 0) >= 60 or int(m.group(4) or 0) >= 60:
        return None
    if int(m.group(2) or 0) or int(m.group(4) or 0):
        return MINUTES
    a = int(m.group(1)) + int(m.group(2) or 0) / 60.0
    b = int(m.group(3)) + int(m.group(4) or 0) / 60.0
    if not (0 <= a <= 24 and 0 <= b <= 24) or a == b:
        return None
    if (a, b) == (0, 24):
        return "always"
    return (a % 24, b % 24)


def _clock(seconds: float) -> str:
    s = int(round(seconds)) % 86400
    return f"{s // 3600:02d}:{s % 3600 // 60:02d}"


def _align_pattern(wn, mult: list, res: FormResult) -> list | None:
    """The form's 24 hourly multipliers (from midnight) as the model's pattern: on a model with no pattern yet (a
    skeleton) the pattern step becomes 1 hour from the run's start; a model's own pattern step and start are kept.
    Either way the multipliers are shifted by the model's START CLOCKTIME less its pattern start (the review: a model
    that starts at 08:00 took midnight's multiplier at 08:00).  None, with a blocking sentence, when the model's
    patterns do not step by the hour."""
    t = wn.options.time
    if not wn.pattern_name_list:
        t.pattern_timestep, t.pattern_start = 3600, 0
    elif int(round(t.pattern_timestep)) != 3600:
        res.blocking.append(f"The model's patterns step every {t.pattern_timestep / 60:g} minutes, and the form's "
                            f"pattern is hourly, so the form's water use was not added: give the water use in the "
                            f"model itself, or leave the model's other patterns out.")
        return None
    off = (float(t.start_clocktime) - float(t.pattern_start)) / 3600.0
    k = int(round(off)) % 24
    if abs(off - round(off)) > 1e-9:
        res.notes.append(f"The model's clock starts at {_clock(t.start_clocktime)}, not on the hour: the form's hourly "
                         f"pattern is aligned to the nearest hour.")
    if k:
        res.notes.append(f"The model's clock starts at {_clock(t.start_clocktime)}: the form's pattern, given from "
                         f"midnight, is shifted so each clock hour keeps its multiplier.")
    return [float(mult[(i + k) % 24]) for i in range(24)]


def _parts(wn) -> list[list[str]]:
    """The model's connected parts (every link joins its two nodes, as importers.prepare_inp counts them), largest
    first."""
    import networkx as nx
    g = nx.Graph()
    g.add_nodes_from(wn.node_name_list)
    g.add_edges_from((lk.start_node_name, lk.end_node_name) for _, lk in wn.links())
    return sorted((sorted(c) for c in nx.connected_components(g)), key=lambda c: (-len(c), c[0]))


def _unfed_parts(wn) -> tuple[int, list[list[str]]]:
    """(the number of connected parts, the parts with a junction and no source: no reservoir and no junction that
    takes water in; a tank alone is no source, as for the whole model)."""
    parts = _parts(wn)
    fed = set(wn.reservoir_name_list) | {j for j, n in wn.junctions()
                                         if any(ts.base_value < 0 for ts in n.demand_timeseries_list)}
    junctions = set(wn.junction_name_list)
    return len(parts), [c for c in parts if not fed & set(c) and junctions & set(c)]


ONE_WAY_VALVES = ("PRV", "PSV", "FCV")       # valves EPANET lets water through in one direction only, start to end


def _controlled_links(wn) -> set:
    """The links some control or rule of the model acts on (one may open a link that starts closed)."""
    out = set()
    for cn in wn.control_name_list:
        for act in wn.get_control(cn).actions():
            name = getattr(act.target()[0], "name", None)
            if name is not None:
                out.add(name)
    return out


def _carries(name: str, lk, controlled: set) -> str | None:
    """How a link can carry water: 'both' ways, 'one' way (start to end: a pump, a check valve, a PRV, PSV or FCV), or
    None: closed with nothing to open it (a closed_valve row, a pump whose rule is 'never', or the model's own link
    closed with no control)."""
    from wntr.network import LinkStatus
    if lk.initial_status == LinkStatus.Closed and name not in controlled:
        return None
    if lk.link_type == "Pump" or (lk.link_type == "Pipe" and lk.check_valve) or (
            lk.link_type == "Valve" and lk.valve_type in ONE_WAY_VALVES):
        return "one"
    return "both"


def _reached(g, src: set) -> set:
    """The nodes of g that water from any of src reaches (src included), in one walk from a node joined to each."""
    import networkx as nx
    top = ("every source",)        # no model node is a tuple
    g.add_edges_from((top, x) for x in src)
    try:
        return nx.descendants(g, top)
    finally:
        g.remove_node(top)


def _cut_off(wn, skip=frozenset()) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    """(the junctions with water use that no source's water can reach, sorted; the links that keep water from them,
    each with 'closed' or 'one way'; the tanks among them, which feed them until they are empty), over the links that
    can carry water (_carries), in the direction each lets it through, and through tanks, which pass water on.  The
    sources are the reservoirs and the junctions that take water in; a tank alone is no source.  Nodes in skip (the parts no link joins to a source, refused on their own) are left
    out.  Task 19's final verification: the parts had been counted over every link, open or not, so a branch the form
    itself cut off (a closed_valve row on its only pipe, a check valve pointing out of it, a pump whose rule is
    'never') was given water use, completed with ok, and EPANET reported its junctions disconnected."""
    import networkx as nx
    controlled = _controlled_links(wn)
    g = nx.DiGraph()
    g.add_nodes_from(wn.node_name_list)
    how = {}
    for name, lk in wn.links():
        a, b = lk.start_node_name, lk.end_node_name
        how[name] = _carries(name, lk, controlled)
        if how[name] is None:
            continue
        g.add_edge(a, b)
        if how[name] == "both":
            g.add_edge(b, a)
    src = set(wn.reservoir_name_list) | {j for j, n in wn.junctions()
                                         if any(ts.base_value < 0 for ts in n.demand_timeseries_list)}
    reached = _reached(g, src)
    cut = sorted(j for j, n in wn.junctions() if j not in reached and j not in skip
                 and any(ts.base_value > 0 for ts in n.demand_timeseries_list))
    if not cut:
        return [], [], []
    # the unreached nodes around the cut-off junctions (joined to them by any link, open or not), and the links between
    # them and the nodes the sources reach: each is closed, or lets water through only out of them
    h = nx.Graph()
    h.add_nodes_from(n for n in wn.node_name_list if n not in reached and n not in skip)
    h.add_edges_from((lk.start_node_name, lk.end_node_name) for _, lk in wn.links()
                     if lk.start_node_name in h and lk.end_node_name in h)
    cut_set = set(cut)
    region = set().union(*(c for c in nx.connected_components(h) if c & cut_set))
    keep = sorted((name, "closed" if how[name] is None else "one way") for name, lk in wn.links()
                  if (lk.start_node_name in region and lk.end_node_name in reached)
                  or (lk.end_node_name in region and lk.start_node_name in reached))
    return cut, keep, sorted(region & set(wn.tank_name_list))


def _cut_off_words(wn, cut: list, keep: list, rows_of: dict, tanks: list = ()) -> str:
    """The blocking sentence for _cut_off's junctions, naming each link that keeps water from them and the form's rows
    that made it so (rows_of: link name -> the FormRows that placed, closed or turned it)."""
    def link(name: str, why: str) -> str:
        lk = wn.get_link(name)
        rows = " and ".join(f"line {r.line}'s {r.kind} row" for r in rows_of.get(name, []))
        kind = {"Pump": "pump", "Valve": str(getattr(lk, "valve_type", "valve"))}.get(lk.link_type, "pipe")
        if why == "closed":
            if lk.link_type == "Pump":
                return f"pump {name}, which never runs ({rows or 'closed in the model, with no control to open it'})"
            return f"{kind} {name}, closed " + (f"by {rows}" if rows else "in the model")
        what = ("which pumps only out of them" if lk.link_type == "Pump" else
                "a check valve that lets water through only out of them" if lk.link_type == "Pipe" else
                "which lets water through only out of them")
        return f"{kind} {name}, {what} ({rows or 'as the model has it'})"
    many = len(cut) != 1
    names = ", ".join(cut[:3]) + (f" and {len(cut) - 3:,} more" if len(cut) > 3 else "")
    ways = "; ".join(link(n, w) for n, w in keep[:4]) + (f"; and {len(keep) - 4:,} more" if len(keep) > 4 else "")
    feed = (f" (only the tank{'s' if len(tanks) != 1 else ''} {', '.join(tanks[:3])} on {'their' if many else 'its'} "
            f"side, which nothing fills: a tank alone is no source)" if tanks else "")
    met = ("Once the tank is empty the water use there is met by nothing" if len(tanks) == 1 else
           "Once the tanks are empty the water use there is met by nothing" if tanks else
           "The water use there would be met by nothing")
    return (f"{len(cut):,} junction{'s' if many else ''} with water use ({names}) can get water from no source{feed}: "
            f"{'they are' if many else 'it is'} joined to the sources only through {ways or 'no link'}. {met} "
            f"(EPANET reports such junctions disconnected), so the model is not written. Open the link or turn it "
            f"round, give the pump a rule that runs it, or add a source row on {'their' if many else 'its'} side.")


def _part_time(wn, skip=frozenset()) -> tuple[list[str], list[str]]:
    """(the junctions with water use that only a link some control or rule turns off at times can feed, with no tank
    on their side to hold water for those hours; those links), each link carrying water as _carries says and the
    controlled ones left out, from the reservoirs, the junctions that take water in and the tanks.  A warning, not a
    refusal (task 19's final verification): such a pump does feed them, while it runs."""
    import networkx as nx
    controlled = _controlled_links(wn)
    g = nx.DiGraph()
    g.add_nodes_from(wn.node_name_list)
    for name, lk in wn.links():
        how = _carries(name, lk, controlled)
        if how is None or name in controlled:
            continue
        g.add_edge(lk.start_node_name, lk.end_node_name)
        if how == "both":
            g.add_edge(lk.end_node_name, lk.start_node_name)
    src = set(wn.reservoir_name_list) | set(wn.tank_name_list) | {
        j for j, n in wn.junctions() if any(ts.base_value < 0 for ts in n.demand_timeseries_list)}
    reached = _reached(g, src)
    alone = sorted(j for j, n in wn.junctions() if j not in reached and j not in skip
                   and any(ts.base_value > 0 for ts in n.demand_timeseries_list))
    feeds = sorted({name for name, lk in wn.links() if name in controlled and _carries(name, lk, controlled)
                    and ((lk.start_node_name in reached) != (lk.end_node_name in reached))}) if alone else []
    return alone, feeds


def use_rows_present(form: OperatorForm) -> bool:
    return any(r.kind in ("production", "large_user") for r in form.rows)


def apply_operator_form(skeleton: Skeleton, form: OperatorForm, snap_ft: float = SNAP_FT) -> FormResult:
    """The full model from a skeleton and the form (FormResult.inp_text, in GPM and feet, written by WNTR, after
    EPANET has run its first hour), or every row it cannot use with its reason and no model (FormResult.refused), or
    what stops the model from running (FormResult.blocking: no source, a part no source feeds, a pump with no rule,
    no water use, a model EPANET refuses).  Never raises for what is in the form or the model: a failure it did not
    foresee comes back as a blocking problem."""
    res = FormResult(refused=list(form.refused), notes=list(form.notes))
    try:
        return _apply(skeleton, form, snap_ft, res)
    except Exception as e:  # noqa: BLE001  (the last guard: a form or model the code did not foresee)
        res.ok, res.inp_text, res.wn = False, None, None
        res.blocking.append(f"The form could not be applied to this model ({type(e).__name__}: {error_words(e, 200)}): "
                            f"a case this tool does not handle yet, in the model or in the form. Check the model and "
                            f"the rows the error names.")
        return res


def _apply(sk: Skeleton, form: OperatorForm, snap_ft: float, res: FormResult) -> FormResult:
    from collections import Counter
    from wntr.network import LinkStatus
    from wntr.network.controls import Control, ControlAction
    from .importers import scan_inp_ids
    # ---- 0. the model's text as task 15's intake reads one: a byte-order mark dropped, and its IDs scanned before
    # WNTR parses it (the review: WNTR keeps the later of two lines with one ID without a word)
    text = sk.inp_text[1:] if sk.inp_text.startswith("﻿") else sk.inp_text
    scan = scan_inp_ids(text)
    if scan.blocking:
        res.blocking.extend(scan.blocking)
        return res
    res.warnings.extend(scan.warnings)
    wn = _load(text)
    P0 = np.array([wn.get_node(n).coordinates for n in wn.node_name_list], float).reshape(-1, 2)
    degrees = sk.source == "inp" and _looks_like_degrees(P0)
    geo = _Geometry(wn, degrees_to_feet(P0) if degrees else None)   # a model in degrees: feet on a tangent plane
    node_names, pipe_names = set(wn.node_name_list), set(wn.pipe_name_list)   # the model as given (step 1 reads it)
    unit = "ft" if sk.source != "inp" or degrees else "model units"
    snap, clear = snap_ft, PIPE_CLEAR_FT
    near_words = f"{snap_ft:g} {unit}" + ("; the model's coordinates read as degrees" if degrees else "")

    def far(d: float, digits: int = 0) -> str:
        return f"{d:,.{digits}f} {unit}"

    def beyond(r: FormRow, col: str, v) -> None:
        refuse(r, f"{col} {_short(v)!r} is far off the network: more than {FAR_UNITS:,.0f} {unit} from every node "
                  f"(check the point's digits)" + swapped(v))

    def refuse(r: FormRow, why: str) -> None:
        res.refused.append({"line": r.line, "kind": r.kind, "reason": why})

    def swapped(v) -> str:
        """For a layer in degrees: a hint when the point, read longitude first, is on the network."""
        parts = re.split(r"\s*[,;]\s*|\s+", str(v).strip())
        if not (sk.frame and sk.frame.get("kind") == "enu") or len(parts) != 2:
            return ""
        xy = _point(f"{parts[1]},{parts[0]}", sk.frame)
        n, d = geo.node(xy) if xy is not None else (None, math.inf)
        return ("; read as longitude,latitude it is on the network: give the latitude first, as Google Earth copies "
                "a place" if n is not None and d <= snap else "")

    def node_of(r: FormRow, col: str) -> str | None:
        v = r.cells.get(col)
        if v in node_names:
            return v
        xy = _point(v, sk.frame)
        if xy is None:
            refuse(r, f"{col} {v!r} is not a node of the model, nor a point x,y"
                      + (" (latitude,longitude for this layer)" if sk.frame and sk.frame.get("kind") == "enu" else ""))
            return None
        if geo.beyond(xy):
            beyond(r, col, v)
            return None
        n, d = geo.node(xy)
        if n is None or d > snap:
            refuse(r, f"{col} {v!r} is off the network: the nearest node is {far(d)} away (more than {near_words})"
                      + swapped(v))
            return None
        return n

    def pipe_of(r: FormRow, col: str = "pipe") -> str | None:
        v = r.cells.get(col)
        if v in pipe_names:
            return v
        xy = _point(v, sk.frame)
        if xy is None:
            refuse(r, f"{col} {v!r} is not a pipe of the model, nor a point x,y on one")
            return None
        if geo.beyond(xy):
            beyond(r, col, v)
            return None
        (p, d), other = geo.pipes(xy)
        if p is None or d > snap:
            refuse(r, f"{col} {v!r} is off the network: the nearest pipe is {far(d)} away (more than {near_words})"
                      + swapped(v))
            return None
        if other is not None and other[1] - d < clear:
            refuse(r, f"{col} {v!r} is about as near to pipe {other[0]} ({far(other[1], 1)} away) as to pipe {p} "
                      f"({far(d, 1)}): give a point along the pipe, away from its ends and other pipes, or its ID")
            return None
        return p

    # ---- 1. resolve every row's places on the model as given (before anything changes it)
    plan = []                      # (row, resolved dict)
    ids = {}                       # form id -> row (pumps and tanks, for the rules)
    used_nodes = {}                # node -> the row that converts it (source, tank)
    c_seen = {}                    # pipe_c material -> its row
    for r in form.rows:
        got = {}
        k = r.kind
        if k in ("source", "tank", "elevation", "large_user"):
            got["node"] = node_of(r, "location")
            if got["node"] is None:
                continue
            if k in ("source", "tank"):
                nd = wn.get_node(got["node"])
                if nd.node_type != "Junction":
                    refuse(r, f"location {r.cells['location']!r} is already a {nd.node_type.lower()} ({got['node']})")
                    continue
                if got["node"] in used_nodes:
                    refuse(r, f"location {r.cells['location']!r} is the node line {used_nodes[got['node']].line} "
                              f"already makes a {used_nodes[got['node']].kind}")
                    continue
                used_nodes[got["node"]] = r
        elif k in ("pump", "prv", "psv"):
            if r.cells.get("location"):
                got["node"] = node_of(r, "location")
                got["pipe"] = pipe_of(r) if got["node"] else None
                if got["node"] is None or got["pipe"] is None:
                    continue
                lk = wn.get_link(got["pipe"])
                if got["node"] not in (lk.start_node_name, lk.end_node_name):
                    refuse(r, f"pipe {r.cells['pipe']!r} ({got['pipe']}) does not end at the node {got['node']}: give "
                              f"a point on a pipe that joins the node, on its upstream side")
                    continue
                got["how"] = "split"
            else:
                got["pipe"] = pipe_of(r)
                got["upstream"] = node_of(r, "upstream") if got["pipe"] else None
                if got["pipe"] is None or got["upstream"] is None:
                    continue
                lk = wn.get_link(got["pipe"])
                if got["upstream"] not in (lk.start_node_name, lk.end_node_name):
                    refuse(r, f"upstream {r.cells['upstream']!r} ({got['upstream']}) is not an end of the pipe "
                              f"{got['pipe']}")
                    continue
                got["how"] = "replace"
            if k == "pump":
                pts = (numbers(r.cells["flow_gpm"]) if r.cells.get("flow_gpm") else None,
                       numbers(r.cells["head_ft"]) if r.cells.get("head_ft") else None)
                if (pts[0] is None) != (pts[1] is None) or (pts[0] and len(pts[0]) != len(pts[1])):
                    refuse(r, "flow_gpm and head_ft must give the same number of points (one design point, or up to "
                              "three curve points)")
                    continue
                if pts[0] and len(pts[0]) not in (1, 3):
                    refuse(r, "give one design point, or three curve points (flow_gpm and head_ft)")
                    continue
                got["curve"] = list(zip(*pts)) if pts[0] else None
        elif k in ("closed_valve", "check_valve"):
            got["pipe"] = pipe_of(r)
            if got["pipe"] is None:
                continue
            if k == "check_valve":
                got["upstream"] = node_of(r, "upstream")
                if got["upstream"] is None:
                    continue
                lk = wn.get_link(got["pipe"])
                if got["upstream"] not in (lk.start_node_name, lk.end_node_name):
                    refuse(r, f"upstream {r.cells['upstream']!r} is not an end of the pipe {got['pipe']}")
                    continue
        elif k == "pattern":
            m = numbers(r.cells["multipliers"])
            if m is None or len(m) != 24 or any(x < 0 for x in m) or sum(m) <= 0:
                refuse(r, "multipliers must be 24 numbers of 0 or more, one for each hour from midnight")
                continue
            got["multipliers"] = m
        elif k == "pipe_c":
            c = number(r.cells["hw_c"])
            hl = str(wn.options.hydraulic.headloss).upper()
            if hl != "H-W":             # task 19's final verification: a C was written as a Darcy-Weisbach roughness
                refuse(r, f"this model's headloss formula is {HEADLOSS_WORDS.get(hl, hl)}, where a pipe's roughness is "
                          f"{ROUGHNESS_WORDS.get(hl, 'another quantity')}, not a Hazen-Williams C: a pipe_c row cannot "
                          f"be used on it (give the roughness in the model itself)")
                continue
            if not C_RANGE[0] <= c <= C_RANGE[1]:
                refuse(r, f"hw_c {c:g} is outside {C_RANGE[0]:g} to {C_RANGE[1]:g}")
                continue
            mat = r.cells["material"].strip()
            if mat.lower() in ("all", "every", "*", "any"):
                got["material"], got["class"], got["rank"] = "all", None, 0
            else:
                _, cls, name = T.material_of(mat)
                if cls is None:
                    refuse(r, f"material {mat!r} is not one the table reads (cast iron, ductile iron, PVC, steel and "
                              f"the others of the template), or 'all'")
                    continue
                if sk.material is None:
                    refuse(r, NO_MATERIALS)
                    continue
                # a class's own word (plastic) gives the pipes of its class no row names more closely
                generic = name == cls and sum(1 for v in T.MATERIAL_CLASS.values() if v == cls) > 1
                got["material"], got["class"], got["rank"] = name, cls, 1 if generic else 2
            if got["material"] in c_seen:
                refuse(r, f"a second pipe_c row for {got['material']} (line {c_seen[got['material']].line} gives one): "
                          f"give one C a material")
                continue
            c_seen[got["material"]] = r
            got["c"] = c
        if k in ("tank", "pump") and r.cells.get("id"):
            if r.cells["id"] in ids:
                refuse(r, f"id {r.cells['id']!r} is used by line {ids[r.cells['id']].line} too")
                continue
            ids[r.cells["id"]] = r
        plan.append((r, got))
    # rules refer to pump and tank rows by id (or to the model's own tanks and pumps)
    for r, got in plan:
        if r.kind != "pump_rule":
            continue
        p = r.cells["pump"]
        if not (p in ids and ids[p].kind == "pump") and p not in wn.pump_name_list:
            refuse(r, f"pump {p!r} is not the id of a pump row (nor a pump of the model)")
            continue
        if r.cells.get("tank"):
            t = r.cells["tank"]
            if not (t in ids and ids[t].kind == "tank") and t not in wn.tank_name_list:
                refuse(r, f"tank {t!r} is not the id of a tank row (nor a tank of the model)")
                continue
        elif _hours(r.cells["on_hours"]) is None:
            refuse(r, f"on_hours {_short(r.cells['on_hours'])!r} is not 'always', 'never' or hours such as '6-22'")
            continue
        elif _hours(r.cells["on_hours"]) == MINUTES:
            refuse(r, f"on_hours {_short(r.cells['on_hours'])!r} gives minutes: the form takes whole clock hours only, "
                      f"such as '6-22' or '06:00-22:00' (round each time to its hour)")
            continue

    def tank_range(t: str) -> tuple[float, float] | None:
        """(min, max) level in feet above its bottom of the tank a rule names: a tank row's, or the model's own."""
        if t in ids and ids[t].kind == "tank":
            lv = (number(ids[t].cells["min_level_ft"]), number(ids[t].cells["max_level_ft"]))
            return lv if None not in lv else None
        if t in wn.tank_name_list:
            tk = wn.get_node(t)
            return tk.min_level / FT, tk.max_level / FT
        return None

    lines_now = {x["line"] for x in res.refused}
    for r, got in plan:     # numbers that must make sense (the review: values EPANET refuses, or no physical value)
        c = {k: number(v) for k, v in r.cells.items() if k in NUMERIC}
        elev = [(k, c[k]) for k in (("head_ft",) if r.kind == "source" else ("base_elev_ft",) if r.kind == "tank" else
                                    ("elevation_ft",) if r.kind in ("prv", "psv", "elevation") else ()) if k in c]
        lv = tank_range(r.cells["tank"]) if r.kind == "pump_rule" and r.cells.get("tank") else None
        if any(not ELEV_RANGE_FT[0] <= v <= ELEV_RANGE_FT[1] for _, v in elev):
            # task 19's final verification: a tank's base at 646000 ft (for 646) was taken
            k, v = next((k, v) for k, v in elev if not ELEV_RANGE_FT[0] <= v <= ELEV_RANGE_FT[1])
            refuse(r, f"{k} {v:,g} ft is outside {ELEV_RANGE_FT[0]:,.0f} to {ELEV_RANGE_FT[1]:,.0f} ft, beyond any "
                      f"ground or water level a system has (a slip in the digits or the unit)")
        elif lv is not None and r.line not in lines_now and c["on_below_ft"] < c["off_above_ft"] and not (
                lv[0] < c["on_below_ft"] and c["off_above_ft"] < lv[1]):
            # task 19's final verification: on below 2 ft on a tank whose minimum is 5 ft never started its pump
            refuse(r, f"on_below_ft {c['on_below_ft']:g} and off_above_ft {c['off_above_ft']:g} must lie inside the "
                      f"tank's range, above its minimum level ({lv[0]:g} ft) and below its maximum ({lv[1]:g} ft): its "
                      f"level never falls below the minimum nor rises above the maximum, so the pump would "
                      f"{'never start' if c['on_below_ft'] <= lv[0] else 'never stop'}")
        elif r.kind == "tank" and min(c["min_level_ft"], c["init_level_ft"], c["max_level_ft"]) < 0:
            refuse(r, "tank levels are feet of water above the tank's bottom, 0 or more")
        elif r.kind == "tank" and not c["min_level_ft"] < c["max_level_ft"]:
            refuse(r, "max_level_ft must be above min_level_ft")
        elif r.kind == "tank" and not (c["min_level_ft"] <= c["init_level_ft"] <= c["max_level_ft"]):
            refuse(r, "init_level_ft must lie between min_level_ft and max_level_ft")
        elif r.kind == "tank" and c["diameter_ft"] <= 0:
            refuse(r, "diameter_ft must be more than 0")
        elif r.kind in ("production", "large_user") and c["gallons_per_day"] < 0:
            refuse(r, "gallons_per_day must be 0 or more")
        elif r.kind == "pump_rule" and r.cells.get("tank") and not c["on_below_ft"] < c["off_above_ft"]:
            refuse(r, "on_below_ft must be below off_above_ft (the pump starts when the tank falls below the first, "
                      "and stops above the second)")
        elif r.kind == "source" and any(x is not None and not DOSE_RANGE[0] <= x <= DOSE_RANGE[1]
                                        for x in (c.get("dose_mgL"), c.get("residual_mgL"))):
            refuse(r, f"dose_mgL and residual_mgL must lie between {DOSE_RANGE[0]:g} and {DOSE_RANGE[1]:g} mg/L (a "
                      f"thousands separator reads '1,500' as 1500: write 1.5)")
        elif r.kind == "pump" and "power_hp" in c and not c["power_hp"] > 0:
            refuse(r, "power_hp must be more than 0")
        elif r.kind == "pump" and got.get("curve") and len(got["curve"]) == 1 and not min(got["curve"][0]) > 0:
            refuse(r, "a design point needs a flow (flow_gpm) and a head (head_ft) above 0")
        elif r.kind == "pump" and got.get("curve") and len(got["curve"]) == 3 and not (
                0 <= got["curve"][0][0] < got["curve"][1][0] < got["curve"][2][0]
                and got["curve"][0][1] > got["curve"][1][1] > got["curve"][2][1] >= 0):
            refuse(r, "three curve points need flows rising from 0 or more and heads falling to 0 or more (EPANET "
                      "reads no other pump curve)")
        elif r.kind in ("prv", "psv") and c["setting_psi"] < 0:
            refuse(r, "setting_psi must be 0 or more")
        elif r.kind in ("prv", "psv") and "diameter_in" in c and not c["diameter_in"] > 0:
            refuse(r, "diameter_in must be more than 0")
    dup = [r for r in form.of("pattern")][1:] + [r for r in form.of("production")][1:]
    for r in dup:
        refuse(r, f"a second {r.kind} row: give one")
    if not use_rows_present(form) and not any(x["kind"] in ("production", "large_user") for x in form.refused):
        for r in form.of("pattern")[:1]:      # task 19's final verification: the row was taken with ok, and changed nothing
            refuse(r, "a pattern row shapes the form's own water use (its production and large_user rows), and the "
                      "form gives none: a model's own water use keeps its own patterns. Leave the row out, or give "
                      "the water use in the form")
    rule_pumps = {}
    for r in form.of("pump_rule"):
        if r.cells["pump"] in rule_pumps:
            refuse(r, f"a second rule for the pump {r.cells['pump']!r} (line {rule_pumps[r.cells['pump']].line} gives "
                      f"one): give one rule a pump")
        rule_pumps.setdefault(r.cells["pump"], r)
    # ---- 1b. rows that cannot all hold at once (the review: each reached an internal error, or a wrong model)
    bad_now = {x["line"] for x in res.refused}
    live = [(r, g) for r, g in plan if r.line not in bad_now]
    src_rows = {g["node"]: r for r, g in live if r.kind == "source"}
    tank_rows = {g["node"]: r for r, g in live if r.kind == "tank"}
    fixed = set(wn.reservoir_name_list) | set(wn.tank_name_list) | set(src_rows) | set(tank_rows)
    # elevation rows that another row overrides or that change nothing (task 19's final verification: each was taken
    # with ok and no word): a second row for one node, a row at a source's or a tank's node, a row against the
    # elevation_ft of a PRV or PSV that puts its own elevation at that node (below a PRV, above a PSV, or a split node)
    valve_at = {}                  # node -> (row, elevation_ft) of the PRV or PSV whose elevation goes there
    for r, g in live:
        if r.kind in ("prv", "psv"):
            if g["how"] == "split":
                at = g["node"]
            else:
                lk = wn.get_link(g["pipe"])
                down = lk.end_node_name if lk.start_node_name == g["upstream"] else lk.start_node_name
                at = down if r.kind == "prv" else g["upstream"]
            valve_at.setdefault(at, (r, number(r.cells["elevation_ft"])))
    elev_row = {}
    for r, g in live:
        if r.kind != "elevation":
            continue
        n, el = g["node"], number(r.cells["elevation_ft"])
        if n in elev_row:
            refuse(r, f"a second elevation row for the node {n} (line {elev_row[n].line} gives "
                      f"{number(elev_row[n].cells['elevation_ft']):g} ft): give one elevation a node")
        elif n in fixed:
            which = (f"the source of line {src_rows[n].line}" if n in src_rows else
                     f"the tank of line {tank_rows[n].line}" if n in tank_rows else
                     f"the model's {wn.get_node(n).node_type.lower()} {n}")
            refuse(r, f"location {_short(r.cells['location'])!r} is {which}, which takes no ground elevation (a "
                      f"source holds its head_ft, a tank stands at its base_elev_ft): leave the row out")
        elif n in valve_at and abs(valve_at[n][1] - el) > 1e-9:
            vr, vel = valve_at[n]
            refuse(r, f"the node {n} takes the elevation_ft of the {vr.kind.upper()} of line {vr.line} ({vel:g} ft), not "
                      f"{el:g} ft: give the node one elevation")
        elev_row.setdefault(n, r)
    # a tank no source can fill: with no pump, and only reservoirs (fixed heads, no pattern) as sources, water rises no
    # higher than the highest source's head, so a tank whose water stands above it drains and never fills again, and
    # what it alone feeds is cut off once it is empty (task 19's final verification: the 8 completed forms its grammar
    # fuzzer still found EPANET reporting disconnected were each such a tank, on the only path from the source)
    heads = [(number(r.cells["head_ft"]), r) for r in src_rows.values()] + [
        (wn.get_node(n).base_head / FT, None) for n in wn.reservoir_name_list]
    pumped = bool(wn.pump_name_list) or any(r.kind == "pump" for r, _ in live) or any(
        ts.base_value < 0 for _, j in wn.junctions() for ts in j.demand_timeseries_list) or any(
        wn.get_node(n).head_pattern_name for n in wn.reservoir_name_list)
    if heads and not pumped:
        top, top_row = max(heads, key=lambda x: x[0])
        for n, r in tank_rows.items():
            floor = number(r.cells["base_elev_ft"]) + number(r.cells["min_level_ft"])
            if floor > top + 1e-6:
                refuse(r, f"the tank's water stands at {floor:g} ft or more (base_elev_ft plus min_level_ft), above the "
                          f"highest head a source holds ({top:g} ft"
                          + (f", line {top_row.line}" if top_row is not None else ", the model's own reservoir")
                          + "), and no pump lifts water to it: it would drain and never fill again. Check both "
                            "elevations (each is a height in feet on one datum), or add the pump that fills it")
    claims = {}                    # pipe -> [(row, 'split' or another kind, node)]
    for r, g in live:
        if r.kind in ("pump", "prv", "psv") and g["how"] == "split":
            node = g["node"]
            if node in src_rows or node in wn.reservoir_name_list:
                which = (f"the source of line {src_rows[node].line}" if node in src_rows else
                         f"the model's reservoir {node}")
                refuse(r, f"location {r.cells['location']!r} is {which}: a source's node cannot be split. For a well "
                          f"and its pump at one point, give the pump in place of the first pipe out of the well (pipe, "
                          f"and upstream the well's node)")
                continue
            if r.kind in ("prv", "psv") and node in fixed:
                refuse(r, f"location {r.cells['location']!r} is a tank or a source: EPANET joins no {r.kind.upper()} "
                          f"straight to one (give it in place of a pipe beside it)")
                continue
            claims.setdefault(g["pipe"], []).append((r, "split", node))
        elif r.kind in ("pump", "prv", "psv"):
            lk = wn.get_link(g["pipe"])
            if r.kind in ("prv", "psv") and {lk.start_node_name, lk.end_node_name} & fixed:
                refuse(r, f"pipe {r.cells['pipe']!r} ({g['pipe']}) ends at a tank or a source: EPANET joins no "
                          f"{r.kind.upper()} straight to one (give a pipe one further along)")
                continue
            claims.setdefault(g["pipe"], []).append((r, r.kind, None))
        elif r.kind in ("closed_valve", "check_valve"):
            claims.setdefault(g["pipe"], []).append((r, r.kind, None))
    for pipe, cl in claims.items():
        r0, role0, node0 = cl[0]
        for r, role, node in cl[1:]:
            if role == role0 == "split" and node == node0:
                continue           # side by side at one node: they share its split (pumps in parallel)
            refuse(r, f"pipe {pipe} is given by line {r0.line} too (a {r0.kind}): one row a pipe, except pumps or "
                      f"valves side by side at one node, which give the same node and pipe")
    pump_rows = [(r, g) for r, g in plan if r.kind == "pump"]
    default_name = {r.line: f"PUMP{i}" for i, (r, g) in enumerate(pump_rows, 1) if not r.cells.get("id")}
    for line, nm in default_name.items():
        if nm in ids and ids[nm].kind == "pump":
            refuse(ids[nm], f"id {nm!r} is the name the form gives the pump of line {line}, which has no id: give "
                            f"that pump an id, or this one another")
    bad_lines = {x["line"] for x in res.refused}
    if res.refused:
        res.refused.sort(key=lambda x: x["line"])
        res.blocking.append(f"{len(res.refused)} row{'s' if len(res.refused) != 1 else ''} of the form cannot be "
                            f"used, so no model was built. Fix {'them' if len(res.refused) != 1 else 'it'} and give the "
                            f"form again.")
        return res
    plan = [(r, g) for r, g in plan if r.line not in bad_lines]
    rows = res.rows
    by = lambda kind: [(r, g) for r, g in plan if r.kind == kind]

    # ---- 2. Hazen-Williams C by material: 'all', then a class's own word, then each material (the most exact wins);
    # a C the layer's own C field gives is kept
    mat_of = sk.material or {}
    given = {}                     # pipe -> (rank, line, C)
    for r, g in by("pipe_c"):
        for pn in pipe_names:
            if pn in sk.c_kept:
                continue
            m = mat_of.get(pn)
            hit = (g["rank"] == 0 or (g["rank"] == 1 and m is not None and T.MATERIAL_CLASS.get(m) == g["class"])
                   or (g["rank"] == 2 and m == g["material"]))
            if hit and (pn not in given or g["rank"] > given[pn][0]):
                given[pn] = (g["rank"], r.line, g["c"])
        res.placed[r.line] = g["material"]
    for pn, (_, _, c) in given.items():
        wn.get_link(pn).roughness = c

    def c_rows() -> None:
        """The readiness row and warnings of step 2, counted over the pipes still in the model once step 4 has put
        pumps and valves in place of some (task 19's final verification: a pipe a PRV or a pump replaced was still
        counted as given the row's C).  Nothing between step 2 and step 4 adds a row or a warning, so their order is
        as before."""
        now = set(wn.pipe_name_list)
        got = {pn: v for pn, v in given.items() if pn in now}
        kept = len(sk.c_kept & now)
        if by("pipe_c"):
            per_line = Counter(line for _, line, _ in got.values())
            per_given = Counter(line for _, line, _ in given.values())
            each = "; ".join(f"{g['material']} {g['c']:g}: {per_line.get(r.line, 0):,}" for r, g in by("pipe_c"))
            rows.append(("Hazen-Williams C", f"{len(got):,} pipe{'s' if len(got) != 1 else ''} given the C of your "
                                             f"form's pipe_c rows ({each})"
                         + (f"; {kept:,} keep the C of the layer's own C field" if kept else ""), FROM_FILE))
            have = Counter(v for v in mat_of.values() if v)
            for r, g in by("pipe_c"):
                if not per_line.get(r.line):
                    listed = ", ".join(f"{m} ({n:,})" for m, n in have.most_common(6)) or "none named"
                    closer = g["rank"] == 1 and any(T.MATERIAL_CLASS.get(m) == g["class"] for m in have)
                    why = ("each pipe it reaches is now a pump or a valve of the form" if per_given.get(r.line) else
                           f"every {g['class']} pipe of the model takes a row that names its material more closely"
                           if closer else f"no pipe of the model is {g['material']}")
                    res.warnings.append(f"Line {r.line} (pipe_c): {why}, so the row gives no pipe its C; the materials "
                                        f"the map names: {listed}.")
        if sk.source in ("gis", "cad"):
            n_table = sum(1 for pn in now if pn not in given and pn not in sk.c_kept)
            if n_table:
                res.warnings.append(age_note(n_table, kept))

    # ---- 3. sources and tanks (junctions become reservoirs and tanks, named by their row's id)
    srcs = []
    rename = {}
    n_kind = {"source": 0, "tank": 0}
    for r, g in by("source") + by("tank"):
        c = r.cells
        n_kind[r.kind] += 1
        want = c.get("id") or f"{r.kind.upper()}{n_kind[r.kind]}"
        if r.kind == "source":
            def make(j, c=c, want=want):
                nm = _unique_name(wn, want)
                wn.add_reservoir(nm, base_head=_to_si(number(c["head_ft"]), "HydraulicHead"),
                                 coordinates=tuple(j.coordinates))
                return nm
        else:
            def make(j, c=c, want=want):
                nm = _unique_name(wn, want)
                wn.add_tank(nm, elevation=_to_si(number(c["base_elev_ft"]), "Elevation"),
                            init_level=_to_si(number(c["init_level_ft"]), "Length"),
                            min_level=_to_si(number(c["min_level_ft"]), "Length"),
                            max_level=_to_si(number(c["max_level_ft"]), "Length"),
                            diameter=_to_si(number(c["diameter_ft"]), "TankDiameter"),
                            coordinates=tuple(j.coordinates))
                return nm
        new = _convert_node(wn, g["node"], make)
        rename[g["node"]] = new
        res.placed[r.line] = new
        if r.kind == "source":
            srcs.append((new, r))
        else:
            ids_tank = c.get("id")
            if ids_tank:
                rename[f"tank:{ids_tank}"] = new
    re_ = lambda n: rename.get(n, n)

    # ---- 4. pumps, PRVs and PSVs: a node split in two (rows side by side at one node share the split), or a pipe
    # replaced
    splits = {}                    # (node, pipe) -> the suction-side node its split made
    halves = {}                    # node -> the suction-side nodes split from it
    form_pumps = {}                # pump row line -> its name in the model
    pump_ids = {}                  # the form's pump id -> its name in the model
    for r, g in by("pump") + by("prv") + by("psv"):
        c = r.cells
        if g["how"] == "split":
            node = re_(g["node"])
            key = (node, g["pipe"])
            if key not in splits:
                lk = wn.get_link(g["pipe"])
                if node not in (lk.start_node_name, lk.end_node_name):     # step 1b keeps this from happening
                    res.blocking.append(f"Line {r.line}: the pipe {g['pipe']} no longer ends at the node {node}, so "
                                        f"the {r.kind} could not be placed; give it the pipe on its own upstream side.")
                    continue
                splits[key] = _split(wn, node, g["pipe"])
                halves.setdefault(node, []).append(splits[key])
            a, b = splits[key], node
            diam = wn.get_link(g["pipe"]).diameter
        else:
            lk = wn.get_link(g["pipe"])
            u = re_(g["upstream"])
            s, e = re_(lk.start_node_name), re_(lk.end_node_name)
            a, b = (s, e) if s == u else (e, s)
            diam, minor = lk.diameter, lk.minor_loss
            wn.remove_link(g["pipe"])
        if r.kind == "pump":
            nm = _unique_name(wn, c.get("id") or default_name[r.line])
            if c.get("power_hp"):
                wn.add_pump(nm, a, b, pump_type="POWER", pump_parameter=_to_si(number(c["power_hp"]), "Power"))
            elif g["curve"]:
                cname = _unique_name(wn, f"{nm}_curve")
                pts = [(_to_si(q, "Flow"), _to_si(h, "HydraulicHead")) for q, h in g["curve"]]
                wn.add_curve(cname, "HEAD", pts)
                wn.add_pump(nm, a, b, pump_type="HEAD", pump_parameter=cname)
            else:
                wn.add_pump(nm, a, b, pump_type="POWER", pump_parameter=1.0)   # set once the production is known
                res.placed.setdefault("derive_pump", []).append(nm)
            form_pumps[r.line] = nm
            if c.get("id"):
                pump_ids[c["id"]] = nm
            res.placed[r.line] = nm
        else:
            nm = _unique_name(wn, c.get("id") or f"{r.kind.upper()}{r.line}")
            if c.get("diameter_in"):
                diam = _to_si(number(c["diameter_in"]), "PipeDiameter")
            wn.add_valve(nm, a, b, diameter=diam, valve_type=r.kind.upper(),
                         minor_loss=0.0 if g["how"] == "split" else minor,
                         initial_setting=_to_si(number(c["setting_psi"]), "Pressure"), initial_status="ACTIVE")
            el = _to_si(number(c["elevation_ft"]), "Elevation")
            # the elevation goes to the node where the setting acts: below a PRV (its downstream node), above a PSV
            # (its upstream node: the review found it on the wrong side); a split node's two halves are at one place,
            # so both take it
            for n in ((a, b) if g["how"] == "split" else (b,) if r.kind == "prv" else (a,)):
                nd = wn.get_node(n)
                if nd.node_type == "Junction":
                    nd.elevation = el
            res.placed[r.line] = nm

    c_rows()

    # ---- 5. closed valves and check valves on pipes
    for r, g in by("closed_valve"):
        wn.get_link(g["pipe"]).initial_status = LinkStatus.Closed
        res.placed[r.line] = g["pipe"]
    for r, g in by("check_valve"):
        lk = wn.get_link(g["pipe"])
        u = re_(g["upstream"])
        if lk.start_node_name != u:
            s, e = lk.start_node, lk.end_node
            lk.start_node, lk.end_node = e, s
            lk.vertices = list(reversed(lk.vertices))
        lk.check_valve = True
        res.placed[r.line] = g["pipe"]

    # ---- 6. elevations (a split node's halves both)
    for r, g in by("elevation"):
        nd = wn.get_node(re_(g["node"]))
        for n in [nd.name] + halves.get(nd.name, []):
            x = wn.get_node(n)
            if x.node_type == "Junction":
                x.elevation = _to_si(number(r.cells["elevation_ft"]), "Elevation")
        res.placed[r.line] = nd.name

    # ---- 7. water use: the pattern, large users at their nodes, the rest by pipe length
    pat_rows = by("pattern")
    mult = list(pat_rows[0][1]["multipliers"]) if pat_rows else list(DEFAULT_PATTERN)
    has_demand = any(ts.base_value != 0 for _, j in wn.junctions() for ts in j.demand_timeseries_list)
    if not has_demand and use_rows_present(form):
        for _, j in wn.junctions():          # a skeleton's zero entries, so each junction keeps the form's alone
            j.demand_timeseries_list.clear()
    use_rows = by("production") + by("large_user")
    mean_m = float(np.mean(mult))
    total_gpd = None
    aligned = _align_pattern(wn, mult, res) if use_rows and not has_demand else None
    if use_rows and has_demand:
        res.blocking.append("The model already has water use at its junctions, so the form's production and large "
                            "users were not added: leave them out of the form, or remove the model's demands.")
    elif use_rows and aligned is not None:
        pname = _unique_pattern(wn, PATTERN_NAME)
        wn.add_pattern(pname, aligned)
        large = 0.0
        for r, g in by("large_user"):
            gpd = number(r.cells["gallons_per_day"])
            node = re_(g["node"])
            nd = wn.get_node(node)
            if nd.node_type != "Junction":
                res.blocking.append(f"Line {r.line}: a large user at {node}, which is a {nd.node_type.lower()}; put it "
                                    f"at a junction.")
                continue
            if gpd > 0:
                nd.add_demand(gpd / GPD_PER_GPM * GPM_M3S / mean_m, pname)
            large += gpd
            res.placed[r.line] = node
        prod = by("production")
        rest = 0.0
        if prod:
            total_gpd = number(prod[0][0].cells["gallons_per_day"])
            rest = total_gpd - large
            if rest < -1e-9 * max(total_gpd, 1.0):
                res.blocking.append(f"The large users use {large:,.0f} gallons a day, more than the production "
                                    f"({total_gpd:,.0f}): check both.")
                rest = 0.0
        else:
            total_gpd = large
        if rest > 0:
            for j, s in length_shares(wn).items():
                if s > 0:
                    wn.get_node(j).add_demand(rest * s / GPD_PER_GPM * GPM_M3S / mean_m, pname)
        n_lu = len(by("large_user"))
        rows.append(("Water use", f"{total_gpd:,.0f} gallons a day: {large:,.0f} at {n_lu:,} large user"
                                  f"{'s' if n_lu != 1 else ''}, {max(rest, 0):,.0f} spread over the junctions by pipe "
                                  f"length", FROM_FILE if not rest else DERIVED))
        rows.append(("Demand pattern", "24 hourly multipliers from your form" if pat_rows else
                     f"24 hourly multipliers, {DEFAULT_PATTERN_SOURCE}", FROM_FILE if pat_rows else ASSUMED))
    elif not use_rows and not has_demand:
        res.blocking.append(NO_USE)

    # ---- 8. pumps with no performance: a design point from the production (assumed)
    for nm in res.placed.get("derive_pump", []):
        q_gpm = (total_gpd or 0.0) / GPD_PER_GPM * max(mult) / mean_m
        heads = [t.elevation + t.max_level for _, t in wn.tanks()]
        suction = [r_.base_head for _, r_ in wn.reservoirs()]
        lift = (max(heads) - max(suction)) / FT if heads and suction else 0.0
        h_ft = max(DERIVED_HEAD_MIN_FT, DERIVED_HEAD_MARGIN * lift)
        p = wn.get_link(nm)
        cname = _unique_name(wn, f"{nm}_curve")
        wn.add_curve(cname, "HEAD", [(_to_si(max(q_gpm, 1.0), "Flow"), _to_si(h_ft, "HydraulicHead"))])
        wn.remove_link(nm)
        wn.add_pump(nm, p.start_node_name, p.end_node_name, pump_type="HEAD", pump_parameter=cname)
        rows.append(("Pump design point", f"{nm}: {max(q_gpm, 1.0):,.0f} gpm at {h_ft:,.0f} ft, derived from the "
                                          f"production's peak hour and the lift to the highest tank (plus 10%)", ASSUMED))
    res.placed.pop("derive_pump", None)

    # ---- 9. pump rules: controls on a tank's level or the clock, and each pump's starting state (at the model's
    # START CLOCKTIME: the review found a '6-22' pump closed at an 08:00 start)
    k_ctl = 0
    ruled = set()
    t0_h = (float(wn.options.time.start_clocktime) / 3600.0) % 24
    for r, g in by("pump_rule"):
        c = r.cells
        pname = pump_ids.get(c["pump"], c["pump"])
        pump = wn.get_link(pname)
        ruled.add(pname)
        if c.get("tank"):
            tname = rename.get(f"tank:{c['tank']}", c["tank"])
            tank = wn.get_node(tname)
            on, off = _to_si(number(c["on_below_ft"]), "HydraulicHead"), _to_si(number(c["off_above_ft"]), "HydraulicHead")
            for oper, thr, st in ((np.less, on, LinkStatus.Open), (np.greater, off, LinkStatus.Closed)):
                k_ctl += 1
                wn.add_control(f"form rule {k_ctl}",
                               Control._conditional_control(tank, "level", oper, thr, ControlAction(pump, "status", st),
                                                            f"form rule {k_ctl}"))
            pump.initial_status = LinkStatus.Open if tank.init_level <= on else LinkStatus.Closed
        else:
            h = _hours(c["on_hours"])
            if h in ("always", "never"):
                pump.initial_status = LinkStatus.Open if h == "always" else LinkStatus.Closed
            else:
                a, b = h
                for t_h, st in ((a, LinkStatus.Open), (b, LinkStatus.Closed)):
                    k_ctl += 1
                    wn.add_control(f"form rule {k_ctl}",
                                   Control._time_control(wn, int(round(t_h * 3600)), "CLOCK_TIME", True,
                                                         ControlAction(pump, "status", st), f"form rule {k_ctl}"))
                on_now = (a <= t0_h < b) if a < b else (t0_h >= a or t0_h < b)
                pump.initial_status = LinkStatus.Open if on_now else LinkStatus.Closed
        res.placed[r.line] = pname
    controlled = set()
    for cn in wn.control_name_list:
        for act in wn.get_control(cn).actions():
            obj = act.target()[0]
            if getattr(obj, "name", None) in wn.pump_name_list:
                controlled.add(obj.name)
    # a pump the form adds needs a rule (a map shows none); a pump the model already had runs as the model says, and
    # the main page's readiness report warns when it has no control.  Every pump the form adds is checked, by its row
    # (the review: an id equal to an unnamed pump's name hid that pump from this check)
    no_rule = sorted(nm for nm in form_pumps.values() if nm not in ruled and nm not in controlled)
    if no_rule:
        unnamed = [nm for line, nm in form_pumps.items() if nm in no_rule and line in default_name]
        res.blocking.append(f"Pump{'s' if len(no_rule) != 1 else ''} with no rule: {', '.join(no_rule[:6])}"
                            f"{f' and {len(no_rule) - 6} more' if len(no_rule) > 6 else ''}. Say when "
                            f"{'each' if len(no_rule) != 1 else 'it'} runs, with a pump_rule row (a tank's levels, its "
                            f"hours, or 'always')"
                            + (", and give the pump row an id for its rule to name" if unnamed else "")
                            + f". {pump_rule_cost()}")

    # ---- 10. a water source for the model and for each of its parts, and the dose or the inlets' residuals
    from .simulate import source_nodes
    inflow = [j for j, n in wn.junctions() if any(ts.base_value < 0 for ts in n.demand_timeseries_list)]
    if not wn.reservoir_name_list and not srcs and not inflow:
        res.blocking.append(NO_SOURCE)
    else:      # (since task 19's final verification on a model fed by inflow junctions alone too)
        n_parts, unfed = _unfed_parts(wn)
        if unfed:
            res.blocking.append(
                f"The model is in {n_parts} parts that no pipe, pump or valve joins (sizes "
                f"{', '.join(str(len(c)) for c in _parts(wn)[:6])}{' and more' if n_parts > 6 else ''}), and "
                f"{len(unfed)} of them {'have' if len(unfed) != 1 else 'has'} no source (for example the part of "
                f"{unfed[0][0]}, {len(unfed[0]):,} node{'s' if len(unfed[0]) != 1 else ''}), so no water reaches "
                f"{'them' if len(unfed) != 1 else 'it'} and the model cannot run. Add a source row in each such part, or "
                f"join it to the rest: the pipe the map lacks, or a larger joining distance where the map's lines "
                f"stop short of each other.")
        # and in the parts a source feeds: junctions with water use that the links the form closes, turns one way or
        # never runs keep every source's water from (task 19's final verification)
        cut, keep, cut_tanks = _cut_off(wn, frozenset().union(*map(frozenset, unfed)) if unfed else frozenset())
        if cut:
            rows_of = {}
            for r, g in plan:
                if r.kind in ("closed_valve", "check_valve", "pump", "prv", "psv", "pump_rule") and r.line in res.placed:
                    rows_of.setdefault(res.placed[r.line], []).append(r)
            res.blocking.append(_cut_off_words(wn, cut, keep, rows_of, cut_tanks))
        alone, feeds = _part_time(wn) if not unfed and not cut else ([], [])
        if alone:                  # fed, but only while a pump (or another controlled link) is on, and no tank there
            many, one = len(alone) != 1, len(feeds) == 1
            names = ", ".join(alone[:3]) + (f" and {len(alone) - 3:,} more" if len(alone) > 3 else "")
            res.warnings.append(
                f"{len(alone):,} junction{'s' if many else ''} with water use ({names}) get{'' if many else 's'} water "
                f"only through {', '.join(feeds[:4])}{' and more' if len(feeds) > 4 else ''}, which "
                f"{'its rule turns' if one else 'their rules turn'} off at times, and no tank on "
                f"{'their' if many else 'its'} side holds water for those hours: while {'it is' if one else 'they are'} "
                f"off, EPANET reports {'them' if many else 'it'} disconnected. Give the pump a rule that keeps it "
                f"running, or add the tank that feeds {'them' if many else 'it'}.")
    vals = {}
    for nm, r in srcs:
        v = r.cells.get("residual_mgL") or r.cells.get("dose_mgL")
        if v:
            vals[nm] = number(v)
    entry = source_nodes(wn) if (wn.reservoir_name_list or srcs) else []
    if vals:
        res.dose_mgL = max(vals.values())
        if len(entry) > 1:
            if set(vals) == set(entry):
                res.inlet_residuals = dict(sorted(vals.items()))
            else:
                res.warnings.append("Not every source has a dose or residual in the form, so every inlet is dosed at "
                                    "one dose, the dose you set on the main page.")
    for nm, r in srcs:
        rows.append(("Source", f"{nm} at the node the form gives, holding {number(r.cells['head_ft']):g} ft"
                               + (f", {vals[nm]:g} mg/L" if nm in vals else ""), FROM_FILE))
    for r, g in by("tank"):
        rows.append(("Tank", f"{res.placed[r.line]}: {number(r.cells['diameter_ft']):g} ft across, levels "
                             f"{number(r.cells['min_level_ft']):g} to {number(r.cells['max_level_ft']):g} ft above its "
                             f"bottom at {number(r.cells['base_elev_ft']):g} ft", FROM_FILE))
    for kind, word in (("pump", "Pump"), ("prv", "Pressure-reducing valve"), ("psv", "Pressure-sustaining valve")):
        for r, g in by(kind):
            if r.line not in res.placed:
                continue
            how = ("at a node split in two, its suction side on the pipe given" if g["how"] == "split" else
                   "in place of the pipe given")
            rows.append((word, f"{res.placed[r.line]}, {how}", FROM_FILE))
    if by("pump_rule"):
        rows.append(("Pump rules", f"{len(by('pump_rule'))} from your form", FROM_FILE))
    if by("closed_valve"):
        rows.append(("Closed valves", f"{len(by('closed_valve'))} pipe{'s' if len(by('closed_valve')) != 1 else ''} "
                                      f"closed", FROM_FILE))
    if by("elevation"):
        rows.append(("Elevations", f"{len(by('elevation'))} from your form; every other junction keeps the model's",
                     FROM_FILE))
    res.counts = {"rows": len(form.rows), "sources": len(srcs), "tanks": len(by("tank")), "pumps": len(by("pump")),
                  "pump_rules": len(by("pump_rule")), "prv": len(by("prv")), "psv": len(by("psv")),
                  "closed_valves": len(by("closed_valve")), "check_valves": len(by("check_valve")),
                  "large_users": len(by("large_user")), "elevations": len(by("elevation")),
                  "pipe_c": len(by("pipe_c")), "controls": len(wn.control_name_list)}
    res.wn = wn
    # ---- 11. a model that cannot run is never written, so it can be neither handed to the main page nor
    # downloaded: what the steps above cannot see, EPANET's first hour does
    if not res.blocking:
        text = write_text(wn)
        err = first_hour_error(text)
        if err:
            res.blocking.append(err)
        else:
            res.inp_text = text
    res.ok = not res.blocking
    return res


# ----------------------------------------------------------------------------- the measured change, from the ladder
FIDELITY_FILES = {n: f"outputs/intake/fidelity_{n}.json" for n in ("Net3", "ky4", "ky10")}
# name -> (file, key path): every number the readiness text quotes, read at run time (as capabilities.SOURCES)
FIDELITY_KEYS = ("age_mae_h", "after_false_lows", "after_lows_found", "n_true_lows", "after_chlorine_mae_mgL")
FIDELITY_FORMAT = {"after_chlorine_mae_mgL": ".3f"}           # every other value is quoted to one decimal
FIDELITY_SOURCES: dict[str, tuple[str, tuple]] = {}
for _n in ("ky4", "ky10"):
    for _r in ("r0", "r1", "r2", "r3", "r4", "r5", "c_material"):
        for _k in FIDELITY_KEYS:
            FIDELITY_SOURCES[f"{_k}_{_r}_{_n}"] = (FIDELITY_FILES[_n], ("rungs", _r, "summary", _k))
for _k in FIDELITY_KEYS:
    FIDELITY_SOURCES[f"{_k}_r6_ky10"] = (FIDELITY_FILES["ky10"], ("rungs", "r6", "summary", _k))


def fidelity_value(name: str, root: str = REPO):
    rel, keys = FIDELITY_SOURCES[name]
    with open(os.path.join(root, rel)) as fh:
        x = json.load(fh)
    for k in keys:
        x = x[k]
    return x


def fidelity_format(name: str) -> str:
    """The format a readiness sentence prints a FIDELITY_SOURCES value in."""
    return next((f for k, f in FIDELITY_FORMAT.items() if name.startswith(k + "_")), ".1f")


def _q(name: str) -> str:
    return format(fidelity_value(name), fidelity_format(name))


def _moved(a: str, b: str, unit: str, where: str) -> str:
    """'rose from a to b h on ky4', 'fell from ...', or 'stayed at a h on ky4', as printed (the review: the sentences
    called a fall a cost)."""
    u = (f" {unit}" if unit else "") + (f" {where}" if where else "")
    if a == b:
        return f"stayed at {a}{u}"
    return f"{'rose' if float(b) > float(a) else 'fell'} from {a} to {b}{u}"


def _change(rung: str, before: str, what: str, chlorine: bool = False, tail: str = "") -> str:
    """'Measured change <what>: in a simulation on two public models (not your system), the model's water-age error
    (or, for C, its daily-minimum chlorine error after 12 samples) rose from a to b on ky4 and ..., and after 12 samples
    its false lows ...'.  None-safe: '' when the files are not there."""
    k = "after_chlorine_mae_mgL" if chlorine else "age_mae_h"
    unit = "mg/L" if chlorine else "h"
    try:
        q = {n: (_q(f"{k}_{before}_{n}"), _q(f"{k}_{rung}_{n}"), _q(f"after_false_lows_{before}_{n}"),
                 _q(f"after_false_lows_{rung}_{n}")) for n in ("ky4", "ky10")}
    except (OSError, KeyError, ValueError, TypeError):
        return ""
    err = ("after 12 samples the model's daily-minimum chlorine error" if chlorine else "the model's water-age error")
    first = " and ".join(_moved(q[n][0], q[n][1], unit, f"on {n}") for n in ("ky4", "ky10"))
    lows = " and ".join(_moved(q[n][2], q[n][3], "", f"on {n}") for n in ("ky4", "ky10"))
    return (f"Measured change {what}: in a simulation on two public models (not your system), {err} {first}, and "
            f"{'its' if chlorine else 'after 12 samples its'} false lows {lows}, on average over the ladder's seeds."
            + (f" {tail}" if tail else ""))


def pump_rule_cost() -> str:
    return _change("r4", "r3", "when pump rules are lost (each pump stays as the file starts it: on ky4 one of its "
                               "two pumps stays off all week, on ky10 its three controlled pumps stay open all week), on top "
                               "of the ladder's steps before it (one C, water use by pipe length, one elevation)")


def demand_cost() -> str:
    return _change("r2", "r1", "when water use is spread by pipe length instead of where it is used, on top of one C "
                               "(130) on every pipe")


def c_cost() -> str:
    return _change("r1", "r0", "when every pipe takes one C (130) instead of its own, the ladder's first step",
                   chlorine=True)


def material_c_cost() -> str:
    return _change("c_material", "r0", "when C comes from the material alone (new-pipe values, blind to age), beside "
                                       "the ladder, from the full model", chlorine=True,
                   tail="These public models carry no pipe ages: on ky10 the step raises its lowest-C pipes to a "
                        "new-pipe value, as a map does for old mains, and ky4's C values are new-pipe ones already.")


def tank_cost() -> str:
    return _change("r5", "r4", "when every tank is 1.3 times too wide, on top of the ladder's steps before it")


def prv_cost() -> str:
    """ky10 only (ky4 has no PRV): its rung 6."""
    try:
        a, b = _q("age_mae_h_r5_ky10"), _q("age_mae_h_r6_ky10")
        f, g = _q("after_false_lows_r5_ky10"), _q("after_false_lows_r6_ky10")
    except (OSError, KeyError, ValueError, TypeError):
        return ""
    return (f"Measured change when PRVs are lost (ky10, a public model with 5 PRVs; ky4 has none), on top of the "
            f"ladder's other steps: in a simulation, not on your system, the model's water-age error "
            f"{_moved(a, b, 'h', '')} and after 12 samples its false lows {_moved(f, g, '', '')}.")


def why_sentences(result: FormResult | None = None, source: str = "cad") -> list[str]:
    """The readiness report's 'why we need this' sentences for a skeleton (or a model the form completed), each with
    the measured change on ky4 and ky10 read from the fidelity files now."""
    out = []
    for fn in (pump_rule_cost, demand_cost, tank_cost, c_cost, material_c_cost, prv_cost):
        s = fn()
        if s:
            out.append(s)
    return out


def completed_readiness(result: FormResult, prepared, built_from: str):
    """The completed model's readiness report: what the drawing or layer gave, what the form gave (each row's source),
    what is assumed, then importers.prepare_inp's rows for the completed model, with the form's blocking problems and
    the measured changes."""
    from .importers import Readiness
    rep = Readiness(name=prepared.name, disinfectant=prepared.disinfectant, sha1=prepared.sha1, sniff=prepared.sniff,
                    not_a_model=prepared.not_a_model, run_cost=prepared.run_cost, prepared=prepared.prepared,
                    same_file=prepared.same_file, flags={**prepared.flags, "operator_form": True})
    rep.add("Built from", built_from, FROM_FILE)
    rep.add("What it is", "a skeleton completed by your operator form: topology and sizes from the map, the rest from "
                          "the form", DERIVED)
    for item, value, src in result.rows:
        rep.add(item, value, src)
    rep.add("Elevations", "the form's where it gives them, else the skeleton's 0 ft: one elevation changes no flow "
                          "unless a pressure valve or a pressure control is there", ASSUMED)
    for item, value, src in prepared.rows:
        if item in ("Text encoding",):
            continue
        if item == "File":
            value = "the model completed from your skeleton and form (EPANET .inp)"
        rep.rows.append((item, value, DERIVED if src == FROM_FILE else src))
    rep.blocking = list(result.blocking) + list(prepared.blocking)
    rep.fixes = list(prepared.fixes)
    rep.warnings = list(result.warnings) + list(prepared.warnings)
    return rep


# ----------------------------------------------------------------------------- command line
def report_and_write(res: FormResult, out: str | None) -> int:
    """Print what the form gave (refused rows, blocking problems, warnings, notes) and write the completed model to
    out: the command line of this module and of cad_read and gis_read (--form)."""
    for ln in res.refused_lines():
        print("refused:", ln)
    for b in res.blocking:
        print("blocking:", b)
    for w in res.warnings:
        print("warning:", w)
    for n in res.notes:
        print("note:", n)
    if res.inp_text and out:
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(res.inp_text)
        print(f"wrote {out}")
    elif res.inp_text:
        print("The model is complete; give --out to write it.")
    return 0 if res.ok else 1


def complete_from_files(sk: Skeleton, form_path: str, out: str | None) -> int:
    """Read the form at form_path, apply it to sk, report and write (a form that is no table is refused in words)."""
    try:
        with open(form_path, "rb") as fh:
            form = read_operator_form(fh.read(), os.path.basename(form_path))
    except FormRefused as e:
        print("refused:", e)
        return 1
    return report_and_write(apply_operator_form(sk, form), out)


def read_model_file(data: bytes, name: str) -> str:
    """The text of the EPANET model the command line is given, identified by its content as task 15's intake
    identifies one (importers.sniff), a model inside a zip read from it in memory; FormRefused in plain words for any
    other file (task 19's final verification: a zipped model, a drawing or a GIS layer had reached WNTR's reader, and
    the operator saw a zip's compressed bytes or a syntax error at line 1), a drawing or a GIS layer pointed to its own
    command line, which takes the form with --form."""
    from .importers import ZipRefused, read_zip_member, sniff
    s = sniff(data, name)
    if s.kind == "epanet_inp_zip":
        if len(s.members) != 1:
            raise FormRefused(f"{s.sentence()} Give one model at a time.")
        try:
            data = read_zip_member(data, s.members[0])
        except ZipRefused as e:
            raise FormRefused(f"The zip is refused: {e}.") from None
        s = sniff(data, s.members[0])
    if s.kind == "epanet_inp":
        return read_model_text(data)
    dxf = s.kind in ("dxf_ascii", "dxf_binary")
    if s.category == "drawing":
        raise FormRefused(f"{s.sentence()} This command completes an EPANET model (.inp). A drawing is completed by its "
                          f"own command, which keeps each pipe's material and the drawing's coordinates: python -m "
                          f"residualmap.cad_read DRAWING.dxf --form FORM.csv --out completed.inp"
                          + ("." if dxf else " (save the drawing as DXF first)."))
    if s.category == "gis":
        raise FormRefused(f"{s.sentence()} This command completes an EPANET model (.inp). A GIS layer is completed by its "
                          f"own command, which keeps each pipe's material and the layer's coordinates: python -m "
                          f"residualmap.gis_read LAYER --form FORM.csv --out completed.inp (a GeoPackage, a zipped "
                          f"shapefile, GeoJSON, KML or KMZ).")
    raise FormRefused(f"{s.sentence()} This command completes an EPANET model (.inp): give the model first and the form "
                      f"second.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="an EPANET .inp file (for a drawing or a GIS layer, use cad_read or gis_read --form)")
    ap.add_argument("form")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    with open(a.model, "rb") as fh:
        data = fh.read()
    try:
        text = read_model_file(data, os.path.basename(a.model))
    except FormRefused as e:
        print("refused:", e)
        return 1
    sk = Skeleton(inp_text=text, name=os.path.basename(a.model))
    return complete_from_files(sk, a.form, a.out)


if __name__ == "__main__":
    sys.exit(main())
