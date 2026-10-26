"""
contract.py: ResidualMap's inputs and outputs as data (iteration 5, task 21).

docs/io_contract.md says the same in plain language; this module is its machine-readable twin.  It holds the columns
each reader takes and the names it accepts for them, the header of every template in docs/templates/, the flags of
every command that takes a utility's files, the label of every input and section on the app's two pages (with the
range and default of each number), the headline numbers and tables the main page shows, the columns of the pilot's
result files, and the limits the readers and the hosted demo apply.

Who reads it:
- the Get started page: its template downloads (TEMPLATE_FILES) and its 'What files it takes' table (FILES_TABLE);
- the main page: the column names that count as the templates' own (GRAB_LOG_TEMPLATE_NAMES);
- the pilot: the keys its --columns flag takes (GRAB_LOG_MAPPING_KEYS);
- the saved check io_contract_matches_code (residualmap/checks.py), which holds this module to the code in both
  directions (every flag, column, template header, label, range, default, limit and output column the code has is
  here, and everything here is in the code) and to docs/io_contract.md: every table of the doc is the one
  doc_tables() renders from this module, under its heading and in its order (any block with a '|' counts as a
  table), the doc's prose quotes a number only in the sentences doc_sentences() renders, every phrase it puts in
  single quotes names a label, section, option or example here, the doc names every name here, and every name it
  puts in backticks is here or in the repository.

    python -m residualmap.contract                  # prints the 'What files it takes' table, checks the templates'
                                                    # headers and says whether the doc's tables are the contract's
    python -m residualmap.contract --write-tables   # rewrites the doc's tables in place from this module

No number in this module is a result: the ranges, defaults and limits are the code's own settings, and the check
compares each with the code.  Since the task 21 final verification that holds for every number the module quotes in
words too: each one is rendered from LIMITS, DEFAULTS, CLI_DEFAULTS or TYPED_SAMPLES_RANGES, whose values the check
compares with the code's constants (and with each add_argument default, read from the source), and the check fails on
a number written into any other string here (labels and column names, which it holds to the code word for word,
aside).  Everything ResidualMap has been tested on is a public model with a simulated truth.
"""
from __future__ import annotations

import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATES_DIR = os.path.join(REPO, "docs", "templates")
DOC = os.path.join(REPO, "docs", "io_contract.md")

# ----------------------------------------------------------------------------- limits (the code's own settings)
LIMITS = {   # name -> (module, constant, value, what it limits); the check compares each value with the constant
    "large model": ("residualmap.importers", "LARGE_JUNCTIONS", 3000,
                    "junctions above which a grid build needs a confirm (the app's checkbox, the pilot's flag)"),
    "date window": ("residualmap.importers", "LOG_WINDOW_MONTHS", 3,
                    "the main page fits the last this many calendar months that hold a reading"),
    "analyzer readings a day": ("residualmap.importers", "ANALYZER_MIN_PER_DAY", 12,
                                "a day of a tap with this many readings or more can be an online analyzer's"),
    "analyzer step": ("residualmap.importers", "ANALYZER_MAX_STEP_MIN", 30,
                      "minutes: such a day whose readings are at most this far apart (the median step) is an "
                      "analyzer day"),
    "analyzer days": ("residualmap.importers", "ANALYZER_MIN_DAYS", 2,
                      "analyzer days that make a tap an online analyzer"),
    "analyzer dense days": ("residualmap.importers", "ANALYZER_DENSE_DAYS", 7,
                            "days with that many readings, at any step, that make a tap an online analyzer too"),
    "analyzer hour": ("residualmap.importers", "ANALYZER_HOUR", 10.0,
                      "an online analyzer keeps the reading nearest this clock hour"),
    "analyzer reach": ("residualmap.importers", "ANALYZER_MAX_GAP_MIN", 60,
                       "minutes: a day keeps its reading nearest the analyzer hour only when it lies within this; a "
                       "day with none keeps none"),
    "highest reading": ("residualmap.importers", "READING_MAX_MGL", 10.0,
                        "a reading above this many mg/L is listed, not used"),
    "hosted upload": ("residualmap.cloud", "UPLOAD_MB", 5, "the hosted demo's largest file, in MB, if its uploads "
                                                           "are turned on"),
    "hosted model": ("residualmap.cloud", "MAX_JUNCTIONS", 300, "the hosted demo's largest model, in junctions"),
}


def _limit(name: str, limits: dict | None = None):
    return (limits or LIMITS)[name][2]


def _hhmm(hours: float) -> str:
    return f"{int(hours):02d}:{int(round((hours - int(hours)) * 60)):02d}"


def analyzer_rule(limits: dict | None = None) -> str:
    """The analyzer rule of importers.read_grab_log in words, from LIMITS (the doc quotes it whole)."""
    v = lambda n: _limit(n, limits)
    return (f"A tap with {v('analyzer readings a day')} or more readings a day at a median step of "
            f"{v('analyzer step')} minutes or less on {v('analyzer days')} or more days, or with "
            f"{v('analyzer readings a day')} or more readings a day on {v('analyzer dense days')} or more days at any "
            f"step, is read as an online analyzer: it keeps one reading a day, the one nearest "
            f"{_hhmm(v('analyzer hour'))} when one lies within {v('analyzer reach')} minutes of it, and a day with "
            f"none keeps no reading and is listed.")


_WORDS = ("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def doc_sentences(limits: dict | None = None) -> tuple[str, ...]:
    """The sentences (or parts of one) of docs/io_contract.md's prose that quote a limit, a default or a setting,
    rendered from LIMITS, DEFAULTS and the other tables here; the doc holds each verbatim, and its prose holds no other
    number (doc_problems), so none of them can go stale.  limits: LIMITS, or another dict of its shape (the check's
    negative controls)."""
    v = lambda n: _limit(n, limits)
    lo, hi = _d("plant log temperatures")
    band = re.search(r"lo(\d+)", EVERY_JUNCTION_COLUMNS[2]).group(1)
    typed = " and ".join([", ".join(f"`{c}` ({w})" for c, w in list(TYPED_SAMPLES_COLUMNS.items())[:-1]),
                          "`{}` ({})".format(*list(TYPED_SAMPLES_COLUMNS.items())[-1])])
    return (f"A model of more than {v('large model'):,} junctions needs a confirm before its grid is built",
            analyzer_rule(limits),
            f"A reading above {v('highest reading'):g} mg/L is listed, not used.",
            f"The main page fits the last {v('date window')} calendar months that hold a reading.",
            f"a {_d('run days')}-day run with hourly hydraulic steps and a {_d('quality step') / 60:g}-minute quality "
            f"step, of which only the last day is scored",
            f"Line ends closer than the joining distance ({_d('GIS joining distance'):g} ft by default, "
            f"{_d('KML joining distance'):g} m for lines drawn by hand in Google Earth) become one node.",
            f"The 'Grab samples' table has {_WORDS[len(TYPED_SAMPLES_COLUMNS)]} columns: {typed}.",
            f"a temperature that is blank, not a number or outside {lo:g} to {hi:g} C (the range of every rate law used "
            f"here; a log whose temperatures are all above it is said to look like Fahrenheit), or a month not written "
            f"YYYY-MM, is refused with the file's name and its rows.",
            f"Its maps are the median in mg/L, the {band}% band's width in mg/L, the probability of being below the "
            f"threshold and the route's sites.")

# ----------------------------------------------------------------------------- defaults (the code's own settings)
# The defaults and settings this module quotes in words: name -> (module, constant, key, value, unit).  The check reads
# the constant (key: an item of a dict, or a function's parameter default; None: the constant itself; module 'app.py':
# a module-level name of the main page, read from its source) and compares it with value, in the code's unit: 'ft'
# means the code holds value feet in metres (topology.FT).  A computed default with no constant behind it is no entry.
DEFAULTS = {
    "free chlorine dose": ("residualmap.pilot", "DEFAULT_DOSE_MGL", "free_chlorine", 1.2, None),
    "chloramine dose": ("residualmap.pilot", "DEFAULT_DOSE_MGL", "chloramine", 2.0, None),
    "free chlorine threshold": ("residualmap.chemistry", "DEFAULT_THRESHOLD_MGL", "free_chlorine", 0.2, None),
    "chloramine threshold": ("residualmap.chemistry", "DEFAULT_THRESHOLD_MGL", "chloramine", 0.5, None),
    "synthetic months": ("residualmap.pilot", "synthetic_log", "months", 6, None),
    "GIS joining distance": ("residualmap.topology", "DEFAULT_TOL_M", "gis", 6, "ft"),
    "KML joining distance": ("residualmap.topology", "DEFAULT_TOL_M", "kml", 10.0, None),
    "chloramine pH": ("app.py", "CA_PH0", None, 8.0, None),
    "chloramine Cl2:N": ("app.py", "CA_CL2N0", None, 4.5, None),
    "dose multipliers": ("residualmap.simgp", "DOSE_GRID", None, (0.90, 0.95, 1.00, 1.05, 1.10), None),
    "chloramine dose multipliers": ("residualmap.simgp", "DOSES_CA", None,
                                    (0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10), None),
    "nitrification watch temperature": ("residualmap.chloramine", "WATCH_MIN_TEMP_C", None, 15.0, None),
    "run days": ("residualmap.simulate", "DURATION_DAYS", None, 7, None),
    "quality step": ("residualmap.simulate", "QUALITY_STEP_S", None, 300, None),
    "plant log temperatures": ("residualmap.pilot", "PLANT_TEMP_C", None, (0.0, 35.0), None),
}
# The add_argument defaults of the command lines (CLI), read from each module's source by the check, both ways: every
# default other than None is here, and nothing else is.
CLI_DEFAULTS = {
    "residualmap.pilot": {"--disinfectant": "free_chlorine", "--window": 3, "--analyzer-hour": 10.0},
    "residualmap.importers": {"--out": "outputs/private/inputs", "--disinfectant": "free_chlorine"},
}


# the main page's 'Grab samples' table: its number columns' ranges (the data editor's column_config, read by the check)
TYPED_SAMPLES_RANGES = {"hour": (0, 23), "mg/L": (0.0, 5.0)}


def _d(name: str):
    return DEFAULTS[name][3]


def _spread(name: str) -> str:
    """A dose axis as the uncertainty it gives: 'by 10%' for multipliers 0.90 to 1.10."""
    lo, hi = min(_d(name)), max(_d(name))
    return f"by {round((hi - 1) * 100):g}%" if abs((1 - lo) - (hi - 1)) < 1e-9 else f"{lo:.2f} to {hi:.2f} times"


# ----------------------------------------------------------------------------- units
UNITS = {
    "chlorine": "mg/L as Cl2: free chlorine, or total chlorine in the chloramine mode (never both in one run)",
    "hour": f"{TYPED_SAMPLES_RANGES['hour'][0]} to {TYPED_SAMPLES_RANGES['hour'][1]}, a model hour counted from the "
            f"model's start time (START CLOCKTIME); shown as a time of day when the model does not start at midnight",
    "water age": "hours",
    "bulk decay": "per day",
    "wall decay": "metres per day",
    "pipe length": "feet or metres on input, as the file says; metres inside the tool",
    "pipe size": "inches or millimetres on input",
    "flow": "the model file's own EPANET units",
}

# ----------------------------------------------------------------------------- the 'What files it takes' table
# The Get started page shows these rows (where files can be given: on the laptop, and in the hosted demo only if its
# uploads are turned on).  Plain words for an operator: no repository paths, no command-line flags, no task numbers.
FILES_TABLE_COLUMNS = ("What", "Formats", "Where it goes", "Needed?")
FILES_TABLE = (
    ("Your network model", "An EPANET model (.inp), as EPANET writes it or as your modelling software exports it; a "
     "zip holding one model is fine", "Main page, '1. Your network': 'Upload my .inp'",
     "Yes, or build a skeleton in step 4 or step 5"),
    ("A GIS layer of your water mains", "GeoPackage (.gpkg), zipped shapefile (.zip), GeoJSON, or a Google Earth KML "
     "or KMZ", "This page, step 4",
     "Only with no hydraulic model: it gives pipes and nodes, and the source, tanks, pumps and demand are added in "
     "EPANET before it runs"),
    ("A CAD drawing of your water mains", "DXF, text or binary (save a DWG as DXF first)", "This page, step 5",
     "Only with no hydraulic model and no GIS layer: topology and sizes only"),
    ("Grab log", "CSV or Excel (.xlsx): the date, the time, the tap ID and one chlorine reading column, free or total "
     "chlorine, never both", "Main page, '2. Your samples', with the demo off",
     "Yes, or type the samples in the main page's table"),
    ("Tap map", "CSV: each tap's ID and the junction ID it is at in the model", "Main page, beside the grab log",
     "Yes, with a grab log"),
    ("Inlet residuals", "CSV: each inlet's node ID in the model and its measured residual in mg/L",
     "Main page, 'My zone's inlets carry different residuals'",
     "Only for a zone whose inlets carry different residuals; free chlorine only"),
    ("Plant log", "CSV: the month and the water temperature (TOC and dose are optional)",
     "The pilot script only; the app does not read it", "Optional"),
    ("Operator form (sources, tanks, pumps, valves, demand)", "A table, one row per item",
     "Not taken by the app yet: it is built and tested, but not offered here",
     "Until it is offered, add those to a skeleton model in EPANET"),
)

# ----------------------------------------------------------------------------- tables a utility gives
# Each reader: the roles it fills from a table's columns, the template's column for each role, and the other names
# it accepts.  'roles' are what the code calls the columns (the keys of its column guess and of a mapping); a reader
# that guesses from the header also takes the common names of lab and field exports (docs/io_contract.md lists some).
GRAB_LOG_ROLES = {          # importers.guess_log_columns and read_grab_log: role -> the template's column (None: none)
    "tap": "tap_id",
    "reading": "free_chlorine_mgL",     # total_chlorine_mgL for a chloraminated system
    "date": "date",
    "time": "time",
    "analyte": None,                    # a long-format export (one row per analyte): its analyte column
    "unit": None,                       # ... and its units column
}
UNIT_DIVISORS = {"µg/L": 1000, "ug/L": 1000, "ppb": 1000}   # a reading in these units is divided by this (the check
#                                                             asks importers.header_unit)
_UNITS_DIVIDED = (", ".join(list(UNIT_DIVISORS)[:-1]) + f" and {list(UNIT_DIVISORS)[-1]} are divided by "
                  f"{max(UNIT_DIVISORS.values())}")
GRAB_LOG_ROLE_MEANING = {   # the doc's roles table: what each role is
    "tap": "the tap or sample point's ID",
    "reading": "the chlorine reading in mg/L; exactly one species, never both",
    "date": "the sample's date (YYYY-MM-DD, or any common form; a date that reads two ways is asked)",
    "time": "the clock time (HH:MM), or none when the time is in the date column",
    "analyte": "a long-format export's analyte column, one row per analyte",
    "unit": f"a long-format export's units column ({_UNITS_DIVIDED})",
}
GRAB_LOG_SPECIES_COLUMNS = {"free_chlorine": "free_chlorine_mgL", "chloramine": "total_chlorine_mgL"}
GRAB_LOG_MAPPING_KEYS = ("tap", "reading", "date", "time", "analyte", "unit", "analyte_value")   # pilot --columns
GRAB_LOG_TEMPLATE_NAMES = frozenset({"date", "time", "tap_id", "free_chlorine_mgL", "total_chlorine_mgL"})
GRAB_LOG_NOT_READ = ("method", "notes")       # template columns kept in the file for the operator, not read
TAP_MAP_ROLES = {"tap": "tap_id", "junction": "junction_id", "role": None}   # role: 'inlet' marks a zone inlet's tap
TAP_MAP_NOT_READ = ("flush_min", "notes")
INLET_ROLES = {"inlet": "inlet_node_id", "residual": "residual_mgL"}          # importers.read_inlet_residuals
INLET_NOT_READ = ("notes",)
PLANT_LOG_REQUIRED = ("month", "temp_C")      # pilot.load_plant refuses a log without them
PLANT_LOG_OPTIONAL = ("toc_mgL", "dose_mgL")  # dose_mgL: the temperature model's monthly dose; toc_mgL: carried only
OPERATOR_FORM_KINDS = ("source", "tank", "pump", "pump_rule", "closed_valve", "check_valve", "prv", "psv",
                       "elevation", "production", "large_user", "pattern", "pipe_c")
OPERATOR_FORM_COLUMNS = ("kind", "id", "location", "pipe", "upstream", "head_ft", "dose_mgL", "residual_mgL",
                         "base_elev_ft", "diameter_ft", "min_level_ft", "max_level_ft", "init_level_ft", "power_hp",
                         "flow_gpm", "pump", "tank", "on_below_ft", "off_above_ft", "on_hours", "setting_psi",
                         "elevation_ft", "diameter_in", "gallons_per_day", "multipliers", "material", "hw_c", "note")
OPERATOR_FORM_ALIASES = {    # other column names the form reader takes, and the column each means
    "type": "kind", "item": "kind", "name": "id", "node": "location", "at": "location", "node_id": "location",
    "pipe_id": "pipe", "suction_pipe": "pipe", "from": "upstream", "upstream_node": "upstream",
    "residual_entering_mgl": "residual_mgL", "residual_mgl": "residual_mgL", "dose_mgl": "dose_mgL",
    "gpd": "gallons_per_day", "gallons_per_day_gpd": "gallons_per_day", "c": "hw_c", "c_value": "hw_c",
    "hazen_williams_c": "hw_c", "notes": "note", "comment": "note",
    "base_elevation_ft": "base_elev_ft", "minimum_level_ft": "min_level_ft", "maximum_level_ft": "max_level_ft",
    "initial_level_ft": "init_level_ft", "setting_psig": "setting_psi", "flow_gpm_points": "flow_gpm",
    "head_ft_points": "head_ft"}
OPERATOR_FORM_KIND_ALIASES = {   # other words the form reader takes in the kind column
    "reservoir": "source", "clearwell": "source", "well": "source", "plant": "source", "closed valve": "closed_valve",
    "boundary_valve": "closed_valve", "check valve": "check_valve", "cv": "check_valve", "pump rule": "pump_rule",
    "large user": "large_user", "demand": "large_user", "pipe c": "pipe_c", "roughness": "pipe_c"}
OPERATOR_FORM_STATUS = "built and tested, not yet in the app (its command lines only)"   # the doc's section 4

# ----------------------------------------------------------------------------- templates (docs/templates/)
TEMPLATES = {     # file -> its header, which io_contract_matches_code compares with the file's first line
    "grab_log_template.csv": ("date", "time", "tap_id", "free_chlorine_mgL", "method", "notes"),
    "grab_log_total_template.csv": ("date", "time", "tap_id", "total_chlorine_mgL", "method", "notes"),
    "tap_map_by_id_template.csv": ("tap_id", "junction_id", "flush_min", "notes"),
    "plant_log_template.csv": ("month", "temp_C", "toc_mgL", "dose_mgL"),
    "inlet_residuals_template.csv": ("inlet_node_id", "residual_mgL", "notes"),
}
TEMPLATE_READER = {"grab_log_template.csv": "grab log", "grab_log_total_template.csv": "grab log",
                   "tap_map_by_id_template.csv": "tap map", "plant_log_template.csv": "plant log",
                   "inlet_residuals_template.csv": "inlet residuals"}
# the Get started page's downloads, in its order (file, button label); the Excel grab log is written on the fly
TEMPLATE_FILES = (("grab_log_template.csv", "Grab log, free chlorine (CSV)"),
                  ("grab_log_total_template.csv", "Grab log, total chlorine for a chloraminated system (CSV)"),
                  ("tap_map_by_id_template.csv", "Tap map by ID (CSV)"),
                  ("plant_log_template.csv", "Plant log, optional, read by the pilot script only (CSV)"))
XLSX_TEMPLATE = ("grab_log_template.xlsx", "Grab log, free chlorine (Excel)")
NOT_ON_THE_PAGE = {"inlet_residuals_template.csv": "named in the pilot's --inlets help; the main page's inlet "
                                                   "uploader says its two columns"}

# ----------------------------------------------------------------------------- command lines
# module -> {flag or positional argument: what it means}.  Only the commands that take a utility's files; the
# experiments, fixtures, checks and numbers scripts are the project's own tools and are documented in the README.
CLI = {
    "residualmap.pilot": {
        "--inp": "the utility's model file, any form the importer takes (a bundled network's name runs it by name)",
        "--samples": "the grab log, CSV or Excel",
        "--taps": "the tap map, CSV",
        "--dose": f"free chlorine leaving the plant, mg/L (default {_d('free chlorine dose'):.1f}); total chlorine "
                  f"with chloramine (default {_d('chloramine dose'):.1f}); with inlet residuals, the highest residual "
                  f"entering the zone",
        "--threshold": f"the minimum residual, mg/L: default {_d('free chlorine threshold'):.1f} free chlorine, "
                       f"{_d('chloramine threshold'):.1f} total chlorine",
        "--disinfectant": f"{CLI_DEFAULTS['residualmap.pilot']['--disinfectant']} (the default) or chloramine",
        "--ph": "chloramine only, with --cl2n: the plant water's logged pH",
        "--cl2n": "chloramine only, with --ph: the chlorine to ammonia-N ratio by mass",
        "--window": f"months of history to fit on (default {CLI_DEFAULTS['residualmap.pilot']['--window']})",
        "--synthetic": f"a bundled network's name: run on a synthetic {_d('synthetic months')}-month log instead of real "
                       f"data",
        "--plant": "the plant log, CSV (free chlorine only): makes the model temperature-aware",
        "--out": "where results go: default the git-ignored private folder for real data, the committed pilot folder "
                 "for --synthetic",
        "--allow-large": f"build the grid of a model with more than {_limit('large model'):,} junctions",
        "--date-format": "how the log's dates are written when they read two ways: MDY, DMY, YMD or a strptime format",
        "--analyzer-hour": f"the clock hour an online analyzer's reading is kept nearest (default "
                           f"{CLI_DEFAULTS['residualmap.pilot']['--analyzer-hour']:g}, for "
                           f"{_hhmm(CLI_DEFAULTS['residualmap.pilot']['--analyzer-hour'])}). {analyzer_rule()}",
        "--keep-every-reading": "keep every reading of a tap that looks like an online analyzer",
        "--columns": "the grab log's columns when the guess is wrong: 'tap=...;reading=...;date=...;time=...', also "
                     "analyte=, unit= and analyte_value=",
        "--inlets": "the inlet residuals CSV, for a zone whose inlets carry different residuals (free chlorine only)",
        "--holdout-sites": "also leave this many sites out of each held-out month's fit, and predict them",
    },
    "residualmap.importers": {
        "files": "one or more files of any kind: each is identified by its content",
        "--prepare": "prepare each EPANET model and print its readiness report",
        "--out": "where a fixed copy goes (default the git-ignored private inputs folder)",
        "--disinfectant": f"{CLI_DEFAULTS['residualmap.importers']['--disinfectant']} (the default) or chloramine, "
                          f"for the readiness report",
        "--json": "print JSON instead of plain text",
    },
    "residualmap.gis_read": {
        "layer": "a GIS layer of the water mains (GeoPackage, zipped shapefile, GeoJSON, KML or KMZ)",
        "--unit": "the coordinates' unit when the layer does not say",
        "--tolerance-ft": f"join line ends closer than this, in feet (default {_d('GIS joining distance'):g} ft; "
                          f"{_d('KML joining distance'):g} m for a KML)",
        "--out": "where the skeleton model (.inp) is written",
        "--form": "an operator form: complete the skeleton with it and write the completed model to --out",
    },
    "residualmap.cad_read": {
        "drawing": "a DXF drawing of the water mains (text or binary, or a .bak holding one)",
        "--unit": "the drawing's unit when it does not say",
        "--tolerance-ft": f"join line ends closer than this, in feet (default {_d('GIS joining distance'):g} ft)",
        "--out": "where the skeleton model (.inp) is written",
        "--form": "an operator form: complete the skeleton with it and write the completed model to --out",
        "--layer-numbers": "numbers in layer names with no unit beside them: none (not sizes), in, mm or m",
        "--default-in": "the size of pipes no layer name or label sizes, in inches",
    },
    "residualmap.operator_form": {
        "model": "an EPANET model whose gaps the form fills",
        "form": "the operator form, CSV or Excel",
        "--out": "where the completed model is written",
    },
}

# ----------------------------------------------------------------------------- the app's inputs and sections
# One entry per place in the code that draws an input: (where, kind, label, (min, max, default) or None, meaning).
# A label with {} is written from values at run time.  A range's item that is a string is computed at run time (the
# string says how); a number is the code's own value.  where: 'sidebar' or 'main' (app.py), 'page' (the Get started
# page).
def _w(where, kind, label, meaning, rng=None):
    return {"where": where, "kind": kind, "label": label, "range": rng, "meaning": meaning}


APP_INPUTS = (
    _w("sidebar", "header", "1. Your network", "the model"),
    _w("sidebar", "radio", "Model", "an example network, your own model file, or a model built on the Get started "
                                    "page (offered while one is kept in this session)"),
    _w("sidebar", "file_uploader", "EPANET .inp file", "your model: any file, identified by its content"),
    _w("sidebar", "selectbox", "Example", "Net3, Net2 or ky4 (Net3 and Net2 only in the hosted demo)"),
    _w("sidebar", "radio", "Disinfectant", "free chlorine, or chloramine (total chlorine); never mixed"),
    _w("sidebar", "number_input", "Highest residual entering the zone (mg/L)",
       f"with inlet residuals: the highest residual entering the zone, treated as uncertain {_spread('dose multipliers')}",
       (0.05, 4.0, "the highest residual given")),
    _w("sidebar", "number_input", "Total chlorine leaving the plant (mg/L as Cl2)",
       f"chloramine: total chlorine at the source (an effective dose of {_spread('chloramine dose multipliers')} "
       f"this)", (0.5, 4.0, 2.0)),
    _w("sidebar", "number_input", "Chlorine dose leaving the plant (mg/L)",
       f"free chlorine at the source, treated as uncertain {_spread('dose multipliers')}", (0.2, 4.0, 1.2)),
    _w("sidebar", "number_input", "Minimum residual, total chlorine (mg/L)",
       "chloramine: the threshold a junction is flagged below (a common utility target, not a California rule)",
       (0.05, 2.0, 0.5)),
    _w("sidebar", "number_input", "Minimum residual (mg/L)", "free chlorine: the threshold a junction is flagged below",
       (0.05, 1.0, 0.2)),
    _w("sidebar", "header", "2. Your samples", "the grab samples"),
    _w("sidebar", "toggle", "Demo: simulate a hidden truth and draw samples from it",
       "on: a simulated network with samples drawn from it; off: your grab log, or samples typed in the table"),
    _w("sidebar", "file_uploader", "Grab log (CSV or Excel)", "your grab-sample log, as your lab or field system "
                                                              "exports it"),
    _w("sidebar", "file_uploader", "Tap map (CSV)", "which junction each tap is at"),
    _w("sidebar", "checkbox", "My zone's inlets carry different residuals", "a zone with several inlets, each with "
                                                                           "its own residual"),
    _w("sidebar", "file_uploader", "Inlet residuals (CSV)", "each inlet's node ID and its measured residual"),
    _w("sidebar", "selectbox", "Tap column", "the grab log's column of tap IDs (guessed; check it)"),
    _w("sidebar", "selectbox", "Reading column", "the grab log's chlorine reading column"),
    _w("sidebar", "selectbox", "Date column", "the grab log's date column"),
    _w("sidebar", "selectbox", "Time column", "the grab log's time column, or the time is in the date column"),
    _w("sidebar", "selectbox", "Analyte column (a log with one row per analyte)", "a long-format export's analyte "
                                                                                 "column"),
    _w("sidebar", "selectbox", "Analyte value (the analyte whose rows are the reading)",
       "which analyte's rows are the chlorine reading"),
    _w("sidebar", "selectbox", "Units column", "a long-format export's units column"),
    _w("sidebar", "selectbox", "Date format", "how the dates are written: guessed, month first, day first, year "
                                              "first, or typed"),
    _w("sidebar", "text_input", "Date format, as strptime writes it (%d day, %b month name, %m month number, %y or "
                                "%Y year, %H:%M time)", "a typed date format"),
    _w("sidebar", "checkbox", "Keep every reading of these taps (a sampling study, not an online analyzer)",
       "keep a tap that looks like an online analyzer whole"),
    _w("sidebar", "number_input", "Online analyzer: keep the reading nearest this hour (0 to 23)",
       f"an online analyzer keeps one reading a day, the one nearest this clock hour when one lies within "
       f"{_limit('analyzer reach')} minutes of it", (0, 23, f"{_limit('analyzer hour'):g}, the reader's hour")),
    _w("sidebar", "slider", "Demo samples (daytime, random junctions)", "the demo's number of samples", (0, 15, 8)),
    _w("sidebar", "number_input", "Demo scenario", "the demo's simulated scenario", (0, 99, 0)),
    _w("sidebar", "number_input", "Plant water pH (logged)", "chloramine: sets the prior centre of the bulk decay rate",
       (6.5, 9.5, f"{_d('chloramine pH'):.1f}, or in the demo its hidden truth's")),
    _w("sidebar", "number_input", "Chlorine to ammonia-N ratio, by mass (logged)", "chloramine: the plant's Cl2:N",
       (3.0, 5.0, f"{_d('chloramine Cl2:N'):.1f}, or in the demo its hidden truth's")),
    _w("sidebar", "number_input", "Water temperature this month (C)",
       f"chloramine: used only by the nitrification watch ({_d('nitrification watch temperature'):g} C or warmer)",
       (0.0, 35.0, 20.0)),
    _w("sidebar", "header", "3. Route", "next month's sampling route"),
    _w("sidebar", "slider", "Sites on next month's route", "how many sites the route picks", (3, 12, 6)),
    _w("main", "page_link", "Bring your data: what to export from your software, and 'Check a file'",
       "to the Get started page, under a file that is not a model"),
    _w("main", "checkbox", "Build the grid anyway: {} runs, about {} min and {} MB (a large model)",
       f"a model of more than {_limit('large model'):,} junctions: confirm before its runs start"),
    _w("main", "radio", "Show", "the map's hour: the daily minimum (the compliance number) or one hour of the day"),
    _w("main", "toggle", "Show the true daily-minimum map and score the flags", "the demo only: reveal the hidden "
                                                                                "truth"),
    _w("main", "download_button", "Download the one-page PDF report", "the PDF report"),
    _w("main", "page_link", "Bring your data: what to export from your software, and templates for the grab log and "
                            "tap map", "to the Get started page (the hosted demo, uploads off)"),
    _w("main", "page_link", "Bring your data: what to export from your software, templates for the grab log and tap "
                            "map, and 'Check a file'", "to the Get started page"),
)
_TR = TYPED_SAMPLES_RANGES
TYPED_SAMPLES_COLUMNS = {"junction": "a junction ID of the model",      # the main page's 'Grab samples' table
                         "hour": f"the model hour, {_TR['hour'][0]:g} to {_TR['hour'][1]:g}",
                         "mg/L": f"the reading, {_TR['mg/L'][0]:g} to {_TR['mg/L'][1]:g} mg/L"}

PAGE_INPUTS = (
    _w("page", "title", "Bring your data", "the Get started page"),
    _w("page", "header", "Step 1. Which software holds your network?", "export steps per software"),
    _w("page", "radio", "Which software holds your network?", "each choice shows its export steps"),
    _w("page", "header", "Step 2. Your grab samples: templates", "the templates"),
    _w("page", "download_button", "Grab log, free chlorine (CSV)", "a template"),
    _w("page", "download_button", "Grab log, total chlorine for a chloraminated system (CSV)", "a template"),
    _w("page", "download_button", "Tap map by ID (CSV)", "a template"),
    _w("page", "download_button", "Plant log, optional, read by the pilot script only (CSV)", "a template"),
    _w("page", "download_button", "Grab log, free chlorine (Excel)", "a template, written on the fly"),
    _w("page", "header", "Step 3. Check a file", "any file: what it is, by its content, and what to send instead"),
    _w("page", "subheader", "What files it takes", "the table of every file ResidualMap takes (above)"),
    _w("page", "radio", "Disinfectant", "the species a grab log or a model is checked for"),
    _w("page", "file_uploader", "Check a file", "any file; nothing is kept"),
    _w("page", "subheader", "What this file is", "the check's answer"),
    _w("page", "header", "Step 4. Build a model from a GIS layer", "a skeleton model from a mains layer"),
    _w("page", "file_uploader", "GIS layer of your water mains", "a GeoPackage, zipped shapefile, GeoJSON, KML or KMZ"),
    _w("page", "selectbox", "Unit of the coordinates", "asked when the layer does not say, or may be wrong"),
    _w("page", "selectbox", "Layer", "which line layer, when the file has several"),
    _w("page", "selectbox", "Order of the coordinates", "asked when latitude may be written first"),
    _w("page", "selectbox", "Pipe size field", "the layer's field with each pipe's size (guessed)"),
    _w("page", "selectbox", "Material field", "the layer's field with each pipe's material (guessed)"),
    _w("page", "selectbox", "Install year field", "the layer's field with each pipe's install year (read, not used)"),
    _w("page", "selectbox", "Declared length field", "the layer's field with each pipe's length (compared with the "
                                                    "drawn length)"),
    _w("page", "selectbox", "Pipe ID field", "the layer's field with each pipe's ID"),
    _w("page", "number_input", "Join line ends closer than (ft)", "the joining distance",
       (0.1, 200.0, f"{_d('GIS joining distance'):g} ft for a GIS layer, {_d('KML joining distance'):g} m for a "
                    f"hand-drawn KML")),
    _w("page", "selectbox", "Unit of the pipe sizes", "asked when sizes are bare numbers"),
    _w("page", "number_input", "Size of every pipe (in)", "asked when no size can be read; empty until given",
       (1.0, 120.0, None)),
    _w("page", "radio", "Is the layer drawn to scale?", "asked for a layer whose unit was chosen"),
    _w("page", "text_area", "Pipes whose true length you know (optional): one per line, 'pipe ID, length in ft'",
       "known lengths, to scale a drawing that is not to scale"),
    _w("page", "download_button", "Download the model (.inp)", "the skeleton model"),
    _w("page", "button", "Open it on the main page", "keeps the skeleton in this session and opens the main page"),
    _w("page", "header", "Step 5. Build a model from a CAD drawing", "a skeleton model from a DXF drawing"),
    _w("page", "file_uploader", "CAD drawing of your water mains (DXF)", "a DXF, text or binary"),
    _w("page", "multiselect", "Layer styles to leave out", "lines and texts that are not mains in service"),
    _w("page", "checkbox", "Read the lines inside blocks inserted on frozen layers", "blocks AutoCAD does not draw"),
    _w("page", "checkbox", "Read pipe sizes from layer names", "a size in a layer's name sizes every line on it"),
    _w("page", "selectbox", "Unit of the drawing", "the drawing's unit, as it says or chosen"),
    _w("page", "number_input", "A label names the nearest line within (text heights)", "a label's reach",
       (0.5, 20.0, 3.0)),
    _w("page", "selectbox", "Unit of the sizes in layer names", "asked when layer names carry bare sizes"),
    _w("page", "selectbox", "Unit of the sizes marked with the diameter symbol", "asked when labels carry bare sizes"),
    _w("page", "selectbox", "Numbers in layer names", "asked: are bare numbers in layer names sizes, and in what unit"),
    _w("page", "number_input", "Join line ends closer than (ft)", "the joining distance", (0.1, 200.0, 6.0)),
    _w("page", "number_input", "Size of every pipe (in)", "asked when no size can be read; empty until given",
       (1.0, 120.0, None)),
    _w("page", "radio", "Is the drawing to scale?", "asked for a drawing whose unit was chosen"),
)

# The options of a choice that this module or the doc names in quotes: option -> (where, kind, label) of the input that
# offers it.  The check reads that input's options from the source (every value the options list can hold).
APP_OPTIONS_NAMED = {"Upload my .inp": ("sidebar", "radio", "Model")}
# Example column headers and readings the doc quotes, each with what the reader makes of it (the check asks the reader)
DOC_EXAMPLES = {"Sample Point": "tap", "Sample Date/Time": "date", "Free Chlorine (mg/L)": "reading",
                "<0.02": "censored", "ND": "censored", "tap=...;reading=...;date=...;time=...": "--columns"}

# ----------------------------------------------------------------------------- what the app shows
APP_TITLE = "ResidualMap"                       # the main page's title (st.title); the check reads every title
APP_SECTIONS = ("Grab samples", "Sample here next", "Every junction: the daily minimum",
                "Nitrification watch (chloramine)", "Water age", "Where chlorine is lost", "When is the network worst?",
                "Demo only: reveal the hidden truth", "What this model accounts for")    # the main page's subheaders
APP_METRICS = ("Samples used", "Junctions likely below {} mg/L ({})", "Worst hour of the day",
               "Calibrated decay (bulk / wall)", "Calibrated decay")
ROUTE_COLUMNS = ("junction", "hour", "P(daily min < {})", "water age (h)", "reason")
EVERY_JUNCTION_COLUMNS = ("junction", "median (mg/L)", "lo90 (mg/L)", "hi90 (mg/L)", "P(daily min < {})",
                          "hour of the daily minimum")
PDF_NAME = "residualmap_{}.pdf"                 # {}: the model file's name, such as Net3.inp

# ----------------------------------------------------------------------------- the pilot's files
PILOT_FILES = ("predictions_{tag}.csv", "validation_{tag}.csv", "validation_{tag}.json", "import_report_{tag}.json",
               "synthetic_log_{network}.csv")
PREDICTIONS_COLUMNS = ("junction", "hour", "y", "month", "date", "tap_id", "pred", "z_mu", "z_sd", "p_below",
                       "train_months", "in50", "in80", "in90", "in95", "persistence", "network_mean", "seen_tap")
VALIDATION_COLUMNS = ("held_out_month", "taps", "trained_on", "n_train", "n_test", "rmse", "mae", "rmse_persistence",
                      "rmse_network_mean", "coverage50", "coverage80", "coverage90", "coverage95", "band90_width_mgL",
                      "n_true_below", "recall_below", "false_alarms", "map_kb", "map_kw", "map_gamma", "map_dose")
SYNTHETIC_LOG_COLUMNS = ("junction", "hour", "y", "month", "date", "tap_id")
COMMITTED_PILOT = {"predictions_Net3.csv": PREDICTIONS_COLUMNS, "validation_Net3.csv": VALIDATION_COLUMNS,
                   "synthetic_log_Net3.csv": SYNTHETIC_LOG_COLUMNS}     # outputs/pilot, the committed synthetic run

# ----------------------------------------------------------------------------- files ResidualMap writes
WRITES = (
    ("outputs/app_uploads/", "the main page's prepared copies of uploaded models (git-ignored)"),
    ("outputs/private/", "real data: the pilot's prepared models and results, the page's checks (git-ignored)"),
    ("outputs/cache/", "the simulator grids, one per model and dose (git-ignored)"),
    ("temp.inp, temp.rpt, temp.bin", "EPANET's scratch files in the working directory, locally (git-ignored)"),
)


def all_labels(entries=APP_INPUTS + PAGE_INPUTS) -> list[str]:
    return [e["label"] for e in entries]


def template_header(name: str) -> str:
    return ",".join(TEMPLATES[name])


def read_template_header(name: str, folder: str = TEMPLATES_DIR) -> str:
    with open(os.path.join(folder, name), encoding="utf-8") as fh:
        return fh.readline().rstrip("\r\n")


# ----------------------------------------------------------------------------- the doc's tables, rendered from here
def _table(head, rows) -> str:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _range_text(e) -> str:
    if e["range"] is None:
        return ""
    lo, hi, d = e["range"]
    d = "empty until given" if d is None else d
    return f"{lo:g} to {hi:g}; default {d if isinstance(d, str) else format(d, 'g')}"


def _labels_table(entries, where) -> str:
    return _table(("label", "kind", "range and default", "what it is"),
                  [(f"'{e['label'].replace('{}', 'N')}'", e["kind"].replace("_", " "), _range_text(e), e["meaning"])
                   for e in entries if e["where"] == where])


def _role_column(role: str) -> str:
    col = GRAB_LOG_ROLES[role]
    if role == "reading":
        return (f"`{GRAB_LOG_SPECIES_COLUMNS['free_chlorine']}`, or `{GRAB_LOG_SPECIES_COLUMNS['chloramine']}` for a "
                f"chloraminated system")
    return f"`{col}`" if col else "(none)"


def doc_tables(limits: dict | None = None) -> list[tuple[str, str]]:
    """Every table of docs/io_contract.md, in its order, as (the '## ' heading it sits under, the table's text),
    rendered from this module.  limits: LIMITS, or another dict of the same shape (the check's negative controls)."""
    lim = limits or LIMITS
    on_page = dict(TEMPLATE_FILES)
    cli = [("8. Command lines", _table(("argument", "what it is"), [(f"`{k}`", v) for k, v in CLI[mod].items()]))
           for mod in CLI]
    return [
        ("Units", _table(("quantity", "unit"), list(UNITS.items()))),
        ("What files it takes", _table(FILES_TABLE_COLUMNS, FILES_TABLE)),
        ("2. Grab samples", _table(("role", "the template's column", "what it is"),
                                   [(f"`{r}`", _role_column(r), GRAB_LOG_ROLE_MEANING[r]) for r in GRAB_LOG_ROLES])),
        ("5. Templates", _table(("template", "header", "read as", "on the Get started page"),
                                [(f"`{f}`", template_header(f), TEMPLATE_READER[f],
                                  f"yes: '{on_page[f]}'" if f in on_page else "no: " + NOT_ON_THE_PAGE.get(f, ""))
                                 for f in TEMPLATES])),
        ("6. The main page's inputs", _labels_table(APP_INPUTS, "sidebar")),
        ("6. The main page's inputs", _labels_table(APP_INPUTS, "main")),
        ("7. The Get started page's inputs", _labels_table(PAGE_INPUTS, "page")),
        *cli,
        ("Outputs", _table(("where", "what"), [(f"`{w}`" if "/" in w else w, why) for w, why in WRITES])),
        ("Limits", _table(("limit", "value", "what it limits"),
                          [(name, f"{v:,g}", why) for name, (_, _, v, why) in lim.items()])),
    ]


HEADING = re.compile(r" {0,3}#{1,6}(?:\s|$)")     # a Markdown heading line (it ends a table, as a blank line does)


def _blocks(text: str) -> list[tuple[int, list[str]]]:
    """The text's blocks: runs of lines that are not blank, as (first line's index, lines); a heading is a block of
    its own."""
    out, cur, first = [], [], 0
    for i, ln in enumerate(text.splitlines() + [""]):
        if HEADING.match(ln):
            if cur:
                out.append((first, cur))
            out.append((i, [ln]))
            cur = []
        elif ln.strip():
            if not cur:
                first = i
            cur.append(ln)
        elif cur:
            out.append((first, cur))
            cur = []
    return out


def tables_in(text: str) -> list[tuple[str | None, str, int]]:
    """The tables of a Markdown text: (the '## ' heading above it, its text, its first line's index).  Since the task
    21 final verification a table is any block (lines between blank lines) with a '|' in any line: GitHub's Markdown
    also renders a row without its end pipes, a row indented up to three spaces, a table with no end pipes at all and
    a line with no pipe right under a table as table rows, and each of those now sits in a block that must be one the
    contract renders, whole, or the doc fails."""
    out, head = [], None
    for first, lines in _blocks(text):
        if lines[0].startswith("## "):
            head = lines[0][3:].strip()
        elif any("|" in ln for ln in lines) and not HEADING.match(lines[0]):
            out.append((head, "\n".join(lines), first))
    return out


QUOTED = re.compile(r"(?<![\w'])'(.+?)'(?![\w])")       # a phrase in single quotes (an apostrophe in a word is none)
PROSE_NUMBER = re.compile(r"(?<![\w.])\d+(?:[.,]\d+)*%?(?![\w])")


def quote_targets() -> tuple[set, list]:
    """What a phrase in single quotes may name: (the exact names, the patterns of names written from values at run
    time, a {} standing for one capital letter such as N): every label, a label without its closing parenthesis, the
    main page's title, sections, metrics and table columns, the options named in APP_OPTIONS_NAMED and the examples
    of DOC_EXAMPLES."""
    names = set()
    for e in APP_INPUTS + PAGE_INPUTS:
        names.add(e["label"])
        names.add(re.sub(r"\s*\([^()]*\)$", "", e["label"]))
    names |= {APP_TITLE, *APP_SECTIONS, *APP_METRICS, *ROUTE_COLUMNS, *EVERY_JUNCTION_COLUMNS, *APP_OPTIONS_NAMED,
              *DOC_EXAMPLES}
    pats = [re.compile("".join("[A-Z]" if p == "{}" else re.escape(p) for p in re.split(r"(\{\})", n)))
            for n in names if "{}" in n]
    return {n for n in names if "{}" not in n}, pats


def quote_problems(text: str, where: str) -> list[str]:
    """The phrases in single quotes in text that name nothing quote_targets() holds (a renamed label, section or
    option leaves its old name behind in prose, or in the 'What files it takes' table the Get started page shows)."""
    names, pats = quote_targets()
    flat = " ".join(re.sub(r"`[^`\n]*`", "``", text).split())
    return [f"{where} quotes {q!r}, which is no label, section, metric, column, named option or example of the "
            f"contract" for q in QUOTED.findall(flat) if q not in names and not any(p.fullmatch(q) for p in pats)]


def step_problems(text: str, where: str) -> list[str]:
    """'step N' in text must be a step of the Get started page (a header 'Step N. ...')."""
    steps = {re.match(r"Step (\d+)\.", e["label"]).group(1) for e in PAGE_INPUTS
             if e["kind"] == "header" and re.match(r"Step \d+\.", e["label"])}
    return [f"{where} names step {n}, and the Get started page has no step {n}"
            for n in re.findall(r"\bstep (\d+)\b", text, re.I) if n not in steps]


def prose_of(text: str) -> str:
    """The doc's prose: its blocks that are no table, without headings, flattened to one line."""
    keep = [ln for first, lines in _blocks(text) if not any("|" in x for x in lines) for ln in lines
            if not ln.startswith("#")]
    return " ".join(" ".join(keep).split())


def doc_problems(text: str, limits: dict | None = None) -> list[str]:
    """What keeps a text of docs/io_contract.md from being the contract's: a table that is not the one doc_tables()
    renders at its place (under its heading, in order; any other table counts), or a limit's sentence
    (doc_sentences()) that is not in it word for word.  Empty when the doc is the contract's."""
    want, got, probs = doc_tables(limits), tables_in(text), []
    if len(got) != len(want):
        probs.append(f"the doc has {len(got)} tables and the contract renders {len(want)}")
    for k, ((w_head, w_text), (g_head, g_text, line)) in enumerate(zip(want, got)):
        if g_head != w_head:
            probs.append(f"table {k + 1} (line {line + 1}) is under '{g_head}', the contract's under '{w_head}'")
        if g_text != w_text:
            gl, wl = g_text.splitlines(), w_text.splitlines()
            j = next((j for j, (a, b) in enumerate(zip(gl, wl)) if a != b), min(len(gl), len(wl)))
            g_ln, w_ln = (gl[j] if j < len(gl) else "(nothing)"), (wl[j] if j < len(wl) else "(nothing)")
            c = max(0, next((c for c, (a, b) in enumerate(zip(g_ln, w_ln)) if a != b), min(len(g_ln), len(w_ln))) - 30)
            probs.append(f"table {k + 1} under '{w_head}' differs at its row {j} (line {line + j + 1}): the doc has "
                         f"{g_ln[c:c + 90]!r}, the contract renders {w_ln[c:c + 90]!r}")
    flat = " ".join(text.split())
    sentences = [" ".join(s.split()) for s in doc_sentences(limits)]
    probs += [f"the doc lacks the sentence {s[:70]!r}" for s in sentences if s not in flat]
    # the task 21 final verification: the prose quotes no number but in a rendered sentence (a second sentence that
    # quotes a limit, or a default written by hand, fails), and every phrase in quotes names something of the contract
    prose = prose_of(text)
    rest = prose
    for s in sentences:
        rest = rest.replace(s, " ")
    rest = re.sub(r"`[^`]*`", " ", rest)
    rest = QUOTED.sub(" ", rest)
    rest = re.sub(r"\bstep \d+\b", " ", rest, flags=re.I)
    probs += [f"the doc's prose quotes the number {m.group()!r} outside the contract's sentences: "
              f"{rest[max(0, m.start() - 50):m.end() + 30]!r}" for m in PROSE_NUMBER.finditer(rest)]
    probs += quote_problems(prose, "the doc's prose") + step_problems(prose, "the doc's prose")
    return probs


def write_tables(text: str) -> str:
    """text with each table replaced, in order, by the one doc_tables() renders; the doc must already have as many
    tables, each under the contract's heading for it (a table added or moved is a change of prose, made by hand)."""
    want, got = doc_tables(), tables_in(text)
    if len(got) != len(want) or any(g[0] != w[0] for g, w in zip(got, want)):
        raise ValueError("the doc's tables are not where the contract puts them: " + "; ".join(
            p for p in doc_problems(text) if "tables" in p or "is under" in p))
    # the task 21 final verification: a table block with a line that is no row of the contract's form (a row without
    # its end pipes, an indented row, prose right under it) is fixed by hand, never written over with the lines next
    # to it
    odd = [first + j + 1 for _, g_text, first in got for j, ln in enumerate(g_text.splitlines())
           if not ln.startswith("|")]
    if odd:
        raise ValueError(f"lines {odd[:6]} of the doc sit in a table and are no row of one; fix them by hand")
    lines = text.splitlines()
    for (_, w_text), (_, g_text, first) in reversed(list(zip(want, got))):
        lines[first:first + len(g_text.splitlines())] = w_text.splitlines()
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="ResidualMap's input and output contract (docs/io_contract.md)")
    ap.add_argument("--write-tables", action="store_true", help="rewrite the doc's tables in place from the contract")
    a = ap.parse_args(argv)
    if a.write_tables:
        with open(DOC, encoding="utf-8") as fh:
            text = fh.read()
        new = write_tables(text)
        with open(DOC, "w", encoding="utf-8") as fh:
            fh.write(new)
        what = "tables rewritten" if new != text else "tables already the contract's"
        print(f"{os.path.relpath(DOC, REPO)}: {what}")
    print(" | ".join(FILES_TABLE_COLUMNS))
    for row in FILES_TABLE:
        print(" | ".join(row))
    bad = 0
    for name in TEMPLATES:
        got = read_template_header(name)
        ok = got == template_header(name)
        bad += not ok
        print(f"{name}: " + ("header as the contract says" if ok else
                              f"header {got!r}, the contract says {template_header(name)!r}"))
    with open(DOC, encoding="utf-8") as fh:
        probs = doc_problems(fh.read())
    bad += bool(probs)
    print(f"{os.path.relpath(DOC, REPO)}: " + ("its tables, numbers and quoted names are the contract's" if not probs else
                                               "; ".join(probs[:6])))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
