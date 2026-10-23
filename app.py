"""
app.py — ResidualMap for an operator.

    streamlit run app.py

Upload the EPANET model you already have, type in the chlorine grab samples you already take, and get:
a chlorine map of every junction with a 90% band, the probability each junction is below the 0.2 mg/L
minimum residual (at any hour, or at its daily minimum), the best next places to sample with a reason
each, and a one-page PDF.  No accounts, no database; everything runs on this laptop.

Demo mode simulates a hidden "true" network from the same file (the experiment's scenario generator)
and draws daytime samples from it, so the app can be shown without real data — and the truth can be
revealed to see how the map did.

Disinfectant (iteration 4, task 12): free chlorine is the default and runs exactly as before.  "Chloramine (total
chlorine)" switches to the chloramine mode: total chlorine readings, the chloramine grid, dose axis, likelihood scale
and prior (residualmap/chloramine.py), a 0.5 mg/L total chlorine default threshold (a common utility operating target,
not a California rule), a nitrification watch (literature thresholds, not validated) and, in demo mode, a hidden truth
from EPA's chloramine chemistry in EPANET-MSX.  Free and total chlorine are never mixed.

Iteration 4, task 14: a 'What this model accounts for' table (water age, organics and chlorine demand, chlorine type,
temperature; every number read from outputs/ through residualmap/capabilities.py), a warning when the calibration sits
on the edge of the decay grid, a labelled fallback for the chloramine demo truth when EPANET-MSX cannot run, and PDF
lines for the disinfectant, the measured species, the threshold and its source, how the seasons are handled, and
'validated in simulation only'.  No temperature or TOC input to the decay model: both were tested in simulation and
not adopted (chloramine mode asks for this month's water temperature for the nitrification watch only).

Iteration 5, task 15: the uploader takes any file and identifies it by its content (residualmap/importers.py).  An
EPANET model goes through prepare_inp: a byte-order mark, UTF-16 or a Windows code page is fixed in a copy, IDs are
checked before WNTR parses, and a 'Model readiness report' says what came from the file, what was derived, what is
assumed and what is missing (task 15 listed the model's physics settings as not checked yet; task 15b checks them,
below).  Any other file gets plain
advice on what to send instead, and the page stops.  A clean UTF-8 .inp is stored exactly as before
(outputs/app_uploads/<sha1>[:8]_<name>), so its grid cache name is unchanged.  .streamlit/config.toml keeps the
server on localhost and usage statistics off, from the next start of the server.

Iteration 5, task 15b: the readiness report checks an uploaded model's physics (reaction orders set to 1, per-pipe and
per-tank bulk lines cleared, Darcy-Weisbach converted to Hazen-Williams behind a water-age gate, Chezy-Manning blocked,
the start clock recorded, pressure-control valves and elevations, coordinates, sources, pumps and controls,
connectivity) and gives its run cost.  A file with one C value (or one Darcy-Weisbach roughness) on every pipe gets
the 'full_fixed_gamma' grid (the old-pipe factor held at 0; in chloramine mode the chloramine grid, which still varies
it) and the old-pipe factor is shown as not identifiable; a model of more than 3,000 junctions
needs the checkbox in its readiness report; a grid that would leave less than 1.5 GB of disk free is refused.  A
network so small that the samples leave fewer junctions than the route asks for gets a shorter route (it stopped the
app before).  The bundled examples never go through prepare_inp and run exactly as before.

Iteration 5, task 16: with the demo off, '2. Your samples' takes the grab log (CSV or Excel) and the tap map as files,
after the model uploader (so file_uploader[0] is still 'EPANET .inp file'), and every new selectbox sits in the sidebar
after the example picker (so selectbox[0] is still 'Example' where it is shown).  importers.read_grab_log reads the
log: its columns guessed from the header and shown for the operator to check, any encoding and separator, decimal
commas, US or European dates (a date that reads two ways stops the page with the question), '<0.02' and 'ND' kept,
listed and marked on the map but not fitted, an online analyzer's 15-minute series thinned to one reading a day, and
free and total chlorine never mixed.  After the review of task 16: a log with a free and a total column asks for its
reading column ('Reading column' starts at '(choose a column)') and reads only that one; a long-format export's
analyte column ('Analyte column') keeps only the chosen species' rows; a units column ('Units column'), µg/L divided
by 1000; a date format can be typed; a tap that looks like an online analyzer can be kept whole (a sampling study);
'>2.2' is marked on the map apart from '<0.02', never as below detection; the rows before the date window are counted
and listed; the expander opens when a guess is not a template name; and on an uploaded file that does not start at
midnight the demo draws its daytime samples at 07:00 to 17:00 by the clock.  The rows of the last three months
prefill the sample table (its key holds a hash of them, so a new upload replaces the table), and every row that is
not used is listed with its reason.  An
uploaded file whose START CLOCKTIME is not midnight has its log's times converted to model hours, and shows the
route's hours, the worst hour, the hour picker and the hour of each junction's daily minimum as times of day, with the
route's window 07:00 to 17:00 by the clock; typed and demo hours are model hours, never converted, and the bundled
examples keep model hours (Net2 with a caption).  A table of every junction's daily minimum (its toolbar downloads the
CSV, so the page keeps one download button) and a plain table of the 20 most likely below the threshold follow the
route, and the page links to the Get started page (pages/1_Get_started.py).  After the final verification of task 16:
the species is read inside brackets too ('Chlorine (Free)'), a second reading column whose header names no species
also asks for the reading column, a long-format export whose analytes do not name the species stops with the
question until 'Analyte value' names the analyte whose rows are read (or '(every row)'), and on the bundled Net2 an
uploaded log's times of day are said to be read as model hours, unconverted.

Iteration 5, task 20 (residualmap/cloud.py): this file is the entrypoint and its own main page.  As the entrypoint it
names the sidebar's pages 'ResidualMap' and 'Get started' (st.navigation) and runs the chosen one; run as the page, it
skips that.  In cloud mode (RESIDUALMAP_MODE=cloud, or the hosted demo's own start, hosted/streamlit_app.py) the Model
radio offers the example networks only (Net3 and Net2), with demo off the samples are typed (no uploader), a privacy
notice opens the page, the chloramine demo goes straight to the reduced truth with its warning, and the last line says
where the numbers are computed.  With RESIDUALMAP_ALLOW_UPLOADS as well, an upload (every file, the grab log and tap
map included, up to cloud.UPLOAD_MB MB, a zip's model too once unpacked, and a model of up to cloud.MAX_JUNCTIONS
junctions, counted in its text before WNTR parses it) is prepared in a per-run temporary folder deleted when the run
ends, and everything made from it (the prepared model, the grid, built in memory, the water age, the demo truth, the
loss split, the grab log) is kept in this browser session only (cloud.session_cached), never in the process-wide
caches, and nothing is written under outputs/.  Since the final verification of task 20 the model's intake runs in a
process of its own with the hosted demo's other limits (links, tanks and reservoirs, lines, time steps, the timed
run) and a time limit (cloud.prepare_upload), every EPANET run on the upload runs in a process whose working
directory is that folder (in_folder), so EPANET's own scratch file goes with it, and an Excel grab log is held to
cloud.WORKBOOK_MB once unpacked (TABLE_LIMITS).  Local mode runs exactly as before.

Iteration 5, task 16c (zone extracts): with the demo off, 'My zone's inlets carry different residuals' (a checkbox
after the tap map uploader, so the uploaders before it are unchanged) offers 'Inlet residuals (CSV)', each inlet's node
ID and measured residual (importers.read_inlet_residuals); a tap map's role column can mark a tap as a zone inlet too,
and the median of its readings over the grab log's window sets that inlet's residual (its readings are never fitted).
Every entry point needs a residual (importers.inlet_ratios_from); each is dosed at r_i times the dose, so the dose
input, written into its place in the sidebar once the inlets are known, is labelled 'Highest residual entering the
zone' and starts at the highest residual given.  The grid at those ratios has its own cache name (and in cloud mode is
kept in the visitor's session, on an example network too); the page and the PDF say which residuals were used and
that water age starts at 0 at every inlet.  Without inlet residuals every path is the one before task 16c.
Its review: a tap the tap map marks as a zone inlet with no usable reading in the window is refused in plain words
(it was left out without a word), its readings outside the window are listed as not used, the inlet file's notes
(decimal commas, µg/L) are shown, the chloramine mode does not offer inlet residuals yet (its dose range is the
plant's), and every stop between the sidebar and the inlets draws the plain dose input first (stop_page).

Iteration 5, task 17: while a model prepared on the Get started page (a skeleton from a GIS layer, or WNTR's own GIS
files) is kept in this browser session, the Model radio offers a third option, 'Model prepared on Get started', and the
page's 'Open it on the main page' picks it; the model goes through the same intake as an upload (its readiness report;
a skeleton stops there, as its report says).  Without one the radio is the one before task 17, so every older check's
widget order holds; with uploads off (the hosted demo) it never appears.  Since task 18 a skeleton built from a CAD
drawing (DXF) on the same page is offered the same way.
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import secrets
import sys
import textwrap
import types
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import wntr
from scipy.stats import norm

from residualmap.age import INITIAL_SHARE_MAX, RANGE_LABEL, hydraulic_age_band, loss_split, oldest_water
from residualmap.capabilities import ACCOUNTS_FOR_NOTE, accounts_for, grid_edge_warning, report_lines
from residualmap.chemistry import FREE_CHLORINE, THRESHOLD_NOTE, cache_tag
from residualmap.features import build_features
from residualmap.importers import (ANALYZER_HOUR, DATE_ORDERS, INLET_ZONE_NOTE, LOG_WINDOW_MONTHS, LogRefused,
                                   clock_label, clock_minute_of, combine_inlet_residuals, daytime_model_hours,
                                   distinct_values, dropped_lines, guess_columns, inlet_ratios_from, inlet_tap_residuals,
                                   mgl_text, prepare_inp, read_grab_log, read_inlet_residuals, read_table, shows_clock,
                                   size_gate)
from residualmap.gis_read import MODEL_RADIO_KEY, PREPARED_KEY, PREPARED_OPEN, PREPARED_OPTION
from residualmap.route import plan_route
from residualmap.simgp import (DAY_HOURS, DOSES_CA, GRIDS, LIK_SD_CA, SimGP24, grid_cache_path,
                               grid_edge_mass, simulator_grid_24h)
from residualmap.simulate import LIB, build_scenario, hidden_chem_draws, nominal_scenario, source_nodes
from residualmap import chloramine as CA
from residualmap import cloud

warnings.filterwarnings("ignore")
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
if not cloud.inside_page():     # task 20: the entrypoint names the pages and runs the chosen one (this file, as the page)
    cloud.run_pages(__file__, os.path.join(APP_ROOT, "pages", "1_Get_started.py"))
    st.stop()
CLOUD = cloud.is_cloud(APP_ROOT)              # task 20: the hosted demo (residualmap/cloud.py); False runs as before
UPLOADS = cloud.uploads_allowed(APP_ROOT)     # always True in local mode
UPLOAD_DIR = os.path.join(APP_ROOT, "outputs", "app_uploads")   # absolute: the app never depends on its working directory
AGE_RESULTS = os.path.join(APP_ROOT, "outputs", "chem", "water_age_{}.json")
CACHE_DIR = os.path.join(APP_ROOT, "outputs", "cache")
CA_RESULTS = os.path.join(APP_ROOT, "outputs", "chloramine", "summary_chloramine_{}.json")
DISINFECTANTS = ["Free chlorine", "Chloramine (total chlorine)"]
PDF_WHY_CHARS = 84     # characters per line of the PDF route table's 'why' column
PDF_CHEM_CHARS = 225   # characters per line of the PDF header's chemistry lines
TEXAS_EXAMPLE = {False: "Texas, for example, requires 0.2 mg/L free chlorine (30 TAC 290.110).",
                 True: "Texas, for example, requires 0.5 mg/L total chlorine in chloraminated systems (30 TAC 290.110)."}
EXAMPLES = {"Net3 (EPANET example, 92 junctions)": "Net3", "Net2 (EPANET example, tank-fed, 35 junctions)": "Net2",
            "ky4 (KYPIPE dataset, 959 junctions; first run takes about 15 min)": "ky4"}
if CLOUD:      # task 20: Net3 and Net2 only; ky4's grid would take hours on the hosted demo's small container
    EXAMPLES = {k: v for k, v in EXAMPLES.items() if v in cloud.CLOUD_EXAMPLES}

st.set_page_config(page_title="ResidualMap", layout="wide")
st.title("ResidualMap")
st.caption("Your EPANET model + the grab samples you already take → a chlorine map of every junction, "
           "with an honest band, the junctions likely below 0.2 mg/L at any hour of the day, and where to sample next.")
if CLOUD:      # task 20: the privacy notice at the top of the page
    st.info(cloud.notice(UPLOADS))
intake_box = st.container()        # the uploaded file's readiness report or advice (iteration 5, task 15)
UPLOAD_HELP = "any file: it is identified by its content"
UPLOAD_MB = cloud.UPLOAD_MB if CLOUD else None     # task 20: the uploaders' own limit in the hosted demo
TABLE_LIMITS = cloud.table_limits() if CLOUD else {}     # task 20: a workbook's size once unpacked, in the hosted demo
# task 16: what to export from each package, templates, 'Check a file'; an absolute path since task 20, so the link
# resolves from the hosted demo's start too (st.navigation registers the page by its path)
GET_STARTED = os.path.join(APP_ROOT, "pages", "1_Get_started.py")
LOG_HELP = ("Your grab-sample log as your lab or field system exports it: CSV (any separator, decimal commas, Windows "
            "encodings) or Excel (.xlsx). The columns are guessed from its header; check them below.")
TAP_HELP = ("Which junction of the model each tap is at: two columns, tap ID and junction ID. The Get started page has "
            "a template.")
CHOOSE = "(choose a column)"
IN_DATE = "(the time is in the date column)"
NO_ANALYTE = "(none: one column per reading)"
NO_UNIT = "(none: the reading column's header says, or mg/L)"
ANALYTE_AUTO = "(the chlorine chosen, found by its name)"      # the final verification of task 16: 'Analyte value'
ANALYTE_EVERY = "(every row: this column does not pick the rows)"
TYPE_FORMAT = "type a format (strptime, such as %d-%b-%y)"
DATE_FORMAT_CHOICES = {"guess (ask when a date reads two ways)": None, **{v: k for k, v in DATE_ORDERS.items()},
                       TYPE_FORMAT: "typed"}
TEMPLATE_NAMES = {"date", "time", "tap_id", "free_chlorine_mgL", "total_chlorine_mgL"}   # docs/templates' headers
KEEP_ALL = "Keep every reading of these taps (a sampling study, not an online analyzer)"
# task 16c: a zone's inlets, each with its own measured residual
INLET_BOX = "My zone's inlets carry different residuals"
INLET_BOX_HELP = ("A pressure zone cut from a larger system has several inlets, each with its own chlorine residual. "
                  "Tick this to give each inlet's residual, so each is dosed in proportion instead of all at one dose.")
INLET_LABEL = "Inlet residuals (CSV)"
INLET_HELP = ("Two columns: each inlet's node ID in your model and its measured residual in mg/L (a typical recent "
              "value). Every inlet needs one. Each inlet is dosed in proportion to its residual, and the dose input then "
              "means the highest residual entering the zone. A tap map can mark a tap as a zone inlet instead, with a "
              "'role' column that says inlet: the median of that tap's readings in the grab log's window sets the "
              "inlet's residual, and its readings are not fitted.")
# task 17: a model prepared on the Get started page that cannot run yet (a skeleton from a GIS layer; since task 18,
# or from a CAD drawing)
PREPARED_NOT_READY = ("This model was built on the Get started page. A model built from a GIS layer or a CAD drawing is a "
                      "skeleton: pipes and nodes only, with no source, tank, pump or demand, so it cannot run until those "
                      "are added: download it on the Get started page and add them in EPANET (a form for them here is "
                      "not built yet).")
# the review of task 16c: inlet residuals are not offered in the chloramine mode until a check covers them
INLET_NOT_CA_WHY = ("Inlet residuals are not available in the chloramine mode yet: its dose range is made for the total "
                    "chlorine leaving a plant, which the plant's fast organic demand lowers, and a residual measured "
                    "where water enters a zone has already passed that demand.")
INLET_NOT_CA = INLET_NOT_CA_WHY + " Every inlet is dosed at the one dose."


@st.cache_data(show_spinner=False, max_entries=16)
def upload_model(data: bytes, name: str, disinfectant: str):
    """importers.prepare_inp on an upload, cached per file, name and disinfectant: (path or None, readiness report).
    A clean UTF-8 .inp is written as outputs/app_uploads/<sha1>[:8]_<name>, as before task 15; a fixed copy as
    <sha1 of the new text>[:8]_<stem>.inp; a file that is not a model, or is blocked, leaves nothing there."""
    return prepare_inp(data, UPLOAD_DIR, disinfectant, name=name)


def session_token() -> str:
    """This browser session's name for its temporary folder (cloud mode, task 20)."""
    return st.session_state.setdefault("_residualmap_scratch", secrets.token_hex(8))


def cloud_upload(data: bytes, name: str, disinfectant: str):
    """upload_model for cloud mode (task 20): prepare_inp in this run's temporary folder (cloud.scratch_dir, deleted
    when the run ends), with the prepared text and the readiness report kept in this browser session only, so a later
    run writes the same file back to the same path and nothing is written under outputs/.  The cloud's limits
    (cloud.intake_limits) refuse a zip's model over the file limit before it is unpacked, and a model over the
    junction limit before WNTR parses it.  Since the final verification of task 20 prepare_inp runs in a process of
    its own, in that folder and with a time limit (cloud.prepare_upload)."""
    folder = cloud.scratch_dir(session_token())

    def prepare_once():
        path, rep = cloud.prepare_upload(data, folder, disinfectant, name, APP_ROOT)
        if path is None:
            return None, rep
        with open(path, "rb") as fh:
            return (os.path.basename(path), fh.read()), rep

    saved, rep = cloud.session_cached(st.session_state, "upload", (hashlib.sha1(data).hexdigest(), name, disinfectant),
                                      prepare_once)
    if saved is None:
        return None, rep
    path = os.path.join(folder, saved[0])
    if not os.path.exists(path):
        with open(path, "wb") as fh:
            fh.write(saved[1])
    return path, rep


def in_folder(fn, *args, **kwargs):
    """fn(*args, **kwargs) for a visitor's upload in cloud mode, in a process of its own whose working directory is
    this run's temporary folder (cloud.run_in_folder), so the scratch file EPANET writes to its working directory is
    deleted with that folder (the final verification of task 20)."""
    return cloud.run_in_folder(cloud.scratch_dir(session_token()), fn, *args, **kwargs)


def run_cached(fn, session: bool, key: tuple, *args, job: tuple | None = None):
    """fn(*args) through its st.cache_data or st.cache_resource decorator, as before task 20 (local mode, and the
    example networks in cloud mode); with session True (a visitor's upload in cloud mode) its undecorated body, kept in
    this browser session only under key (cloud.session_cached).  job: for a function that runs EPANET, the same work
    as an importable function and its arguments, run in this run's folder (in_folder) in its place."""
    if not session:
        return fn(*args)
    if job is not None:
        return cloud.session_cached(st.session_state, fn.__name__, key, in_folder, *job)
    return cloud.session_cached(st.session_state, fn.__name__, key, fn.__wrapped__, *args)


def file_key(*args) -> tuple:
    """A session cache key for arguments that hold a file's bytes: the bytes by their SHA-1."""
    return tuple(hashlib.sha1(a).hexdigest() if isinstance(a, bytes) else a for a in args)


@st.cache_data(show_spinner=False, max_entries=16)
def log_header(data: bytes, name: str, disinfectant: str):
    """The grab log's header and the columns guessed from it (task 16, importers.guess_columns, an analyte column
    found by its values included): (header, guess), or (None, why) when the file cannot be read as a table."""
    try:
        t = read_table(data, name, **TABLE_LIMITS)
    except LogRefused as e:
        return None, str(e)
    except Exception as e:  # noqa: BLE001  (the last guard: a file the reader did not foresee)
        return None, f"The file could not be read ({type(e).__name__}). Save it again as CSV and upload that."
    return t.header, guess_columns(t, disinfectant)


@st.cache_data(show_spinner=False, max_entries=16)
def column_values(data: bytes, name: str, column: str) -> list[str]:
    """The distinct values of one column of the grab log (the analyte column's, for 'Analyte value'), at most 200."""
    try:
        return distinct_values(read_table(data, name, **TABLE_LIMITS), column)
    except Exception:  # noqa: BLE001  (log_header has already said why)
        return []


@st.cache_data(show_spinner=False, max_entries=16)
def grab_log(data: bytes, name: str, taps: bytes, taps_name: str, disinfectant: str, mapping: tuple,
             date_format: str | None, start_clock_h: float, analyzer_hour: float, thin: bool = True):
    """importers.read_grab_log on the uploads, cached per content and choice: (GrabLog, None), or (None, why)."""
    try:
        return read_grab_log(data, name, taps, taps_name, disinfectant, date_format=date_format, mapping=dict(mapping),
                             start_clock_h=start_clock_h, analyzer_hour=analyzer_hour, thin_analyzers=thin,
                             **TABLE_LIMITS), None
    except LogRefused as e:
        return None, str(e)
    except Exception as e:  # noqa: BLE001  (the last guard: a file the reader did not foresee)
        return None, (f"The grab log or the tap map could not be read ({type(e).__name__}). Save them again as CSV and "
                      f"upload those.")


@st.cache_data(show_spinner=False, max_entries=16)
def inlet_file(data: bytes, name: str, disinfectant: str = "free_chlorine"):
    """importers.read_inlet_residuals on the upload (task 16c): (({inlet: mg/L}, notes), None), or (None, why).  Its
    review: read for the disinfectant's species, and how it was read (decimal commas, µg/L) is shown."""
    try:
        res, info = read_inlet_residuals(data, name, disinfectant=disinfectant, **TABLE_LIMITS)
        return (res, info["notes"]), None
    except LogRefused as e:
        return None, str(e)
    except Exception as e:  # noqa: BLE001  (the last guard: a file the reader did not foresee)
        return None, f"The inlet residuals file could not be read ({type(e).__name__}). Save it again as CSV and upload it."


def show_intake(box, rep) -> bool:
    """The upload's advice (not a model), its blocking problems, or its readiness report, in the main area.  For a
    model of more than importers.LARGE_JUNCTIONS junctions (task 15b) the report holds a checkbox: True when the grid
    may be built (always True for a smaller model).  Since task 16 the advice links to the Get started page."""
    if rep.not_a_model:
        box.warning(rep.advice())
        box.page_link(GET_STARTED, label="Bring your data: what to export from your software, and 'Check a file'")
        return False
    if rep.blocking:
        box.error("This model cannot run yet. " + " ".join(rep.blocking))
    large = bool(rep.run_cost.get("large"))
    with box.expander("Model readiness report", expanded=bool(rep.blocking or rep.fixes or rep.warnings)):
        st.markdown(rep.to_markdown())
        if large and not rep.blocking:
            return st.checkbox(f"Build the grid anyway: {rep.run_cost['runs_per_dose']} runs, about "
                               f"{rep.run_cost.get('minutes', 0):g} min and {rep.run_cost['grid_file_mb']:g} MB "
                               f"(a large model)", value=False)
    return True


def old_pipe_text(g: float, uniform_c: bool) -> str:
    """The old-pipe factor for the page and the PDF: its calibrated value, or, for a file with one C value on every
    pipe (task 15b), that it cannot be learned."""
    return "old-pipe factor: not identifiable from one C value" if uniform_c else f"old-pipe factor γ={g:.1f}"


class PreparedModel:
    """A model prepared on the Get started page (task 17: a skeleton from a GIS layer, or WNTR's own GIS files), read
    through the same intake as an upload: name, size and getvalue(), as Streamlit's uploaded file has them."""

    def __init__(self, entry: dict):
        self.name, self._data, self.source = entry["name"], entry["data"], entry.get("from", "the Get started page")
        self.size = len(self._data)

    def getvalue(self) -> bytes:
        return self._data


# ----------------------------------------------------------------------------- inputs
with st.sidebar:
    st.header("1. Your network")
    # task 17: a third option, only while a model prepared on the Get started page is kept in this session; without one
    # the radio is the one before task 17, unkeyed
    prepared = st.session_state.get(PREPARED_KEY) if UPLOADS else None
    model_options = ["Example network"] + (["Upload my .inp"] if UPLOADS else [])
    if prepared:
        model_options.append(PREPARED_OPTION)
        if st.session_state.pop(PREPARED_OPEN, False):
            st.session_state[MODEL_RADIO_KEY] = PREPARED_OPTION
        if st.session_state.get(MODEL_RADIO_KEY) not in model_options:
            st.session_state.pop(MODEL_RADIO_KEY, None)
        src = st.radio("Model", model_options, horizontal=True, key=MODEL_RADIO_KEY)
    else:
        src = st.radio("Model", model_options, horizontal=True)
    up = None
    if src == "Upload my .inp":
        up = st.file_uploader("EPANET .inp file", type=None, help=UPLOAD_HELP, max_upload_size=UPLOAD_MB)
        if up is None:
            st.stop()
    elif src == PREPARED_OPTION:
        up = PreparedModel(prepared)
        st.caption(f"{up.name}, built on the Get started page from {up.source}.")
    else:
        net_path = EXAMPLES[st.selectbox("Example", list(EXAMPLES))]
    ca = st.radio("Disinfectant", DISINFECTANTS, index=0,
                  help="A chloraminated system measures TOTAL chlorine, which decays far more slowly than free chlorine. "
                       "The chloramine mode uses its own decay ranges, dose axis, likelihood scale and threshold; free and "
                       "total chlorine readings are never mixed.") == DISINFECTANTS[1]
    uniform_c = False          # task 15b: one C value on every pipe (an uploaded file only; the examples have several)
    clock_h = 0.0              # task 16: an uploaded file's START CLOCKTIME, in hours (0 for every bundled example)
    if up is not None:
        disinfectant = "chloramine" if ca else "free_chlorine"
        if CLOUD and cloud.too_big(up.size):        # task 20: the hosted demo's limit, before the file is read
            intake_box.error(cloud.too_big(up.size))
            st.stop()
        try:            # prepare_inp turns what is in a file into plain problems; this is the last guard for the rest
            if CLOUD:
                net_path, readiness = cloud_upload(up.getvalue(), up.name, disinfectant)
            else:
                net_path, readiness = upload_model(up.getvalue(), up.name, disinfectant)
                if net_path is not None and not os.path.exists(net_path):     # the stored copy was removed since
                    upload_model.clear()
                    net_path, readiness = upload_model(up.getvalue(), up.name, disinfectant)
        except Exception as e:  # noqa: BLE001
            net_path, readiness = None, None
            intake_box.error(f"This file could not be read ({type(e).__name__}). Send the model's EPANET export "
                             f"(.inp) again, on its own.")
        if CLOUD and readiness is not None and cloud.limit_refusal(readiness):     # before WNTR parses it (task 20)
            intake_box.error(cloud.limit_refusal(readiness))
            st.stop()
        confirmed = show_intake(intake_box, readiness) if readiness is not None else False
        if net_path is None:
            if src == PREPARED_OPTION:         # task 17: a skeleton (GIS layer or drawing) stops here, as its report said
                intake_box.info(PREPARED_NOT_READY)
            st.stop()
        if not confirmed:
            intake_box.info("A large model: tick 'Build the grid anyway' in the readiness report to start its runs.")
            st.stop()
        uniform_c = bool(readiness.flags.get("uniform_c"))
        clock_h = float(readiness.flags.get("start_clock_h") or 0.0)
    # task 16c: the dose input is written into this place once the inlets are known (below the sidebar), so its label
    # can say 'the highest residual entering the zone' when inlet residuals are given; same position, same widget.
    # Its review: every stop before then (stop_page) draws the plain dose input first, as before task 16c
    dose_box = st.container()
    dose_drawn: list = []

    def draw_dose(inlet_info=None) -> float:
        """The sidebar's dose input in its place (dose_box), drawn once per run: the plant's dose, or with inlet
        residuals (task 16c) the highest residual entering the zone, starting at the highest given."""
        if dose_drawn:
            return dose_drawn[0]
        with dose_box:
            if inlet_info is not None:
                top = float(inlet_info["highest_mgL"])
                v = st.number_input("Highest residual entering the zone (mg/L)", 0.05, 4.0,
                                    float(min(max(round(top, 2), 0.05), 4.0)), 0.05,
                                    help="Your inlet residuals dose each inlet in proportion, and this is the highest of "
                                         "them (the highest you gave, unless you change it). The model treats it as "
                                         "uncertain, by ±10%.")
            elif ca:
                v = st.number_input("Total chlorine leaving the plant (mg/L as Cl2)", 0.5, 4.0, 2.0, 0.1,
                                    help="Total chlorine (chloramine) at the source. The model treats it as uncertain and "
                                         "lets the fast organic demand lower it (an effective dose of 0.65 to 1.10 times "
                                         "this).")
            else:
                v = st.number_input("Chlorine dose leaving the plant (mg/L)", 0.2, 4.0, 1.2, 0.1,
                                    help="Free chlorine at the source. The model treats it as ±10% uncertain.")
        dose_drawn.append(v)
        return v

    def stop_page():
        """st.stop, with the sidebar's dose input drawn first (the review of task 16c: a stop between the sidebar and the
        inlets left the page without it)."""
        draw_dose()
        st.stop()
    if ca:
        threshold = st.number_input("Minimum residual, total chlorine (mg/L)", 0.05, 2.0, 0.5, 0.05,
                                    help=f"Default {CA.THRESHOLD_NOTE}. {TEXAS_EXAMPLE[True]} Set your own target here.")
    else:
        threshold = st.number_input("Minimum residual (mg/L)", 0.05, 1.0, 0.2, 0.05,
                                    help=f"Default {THRESHOLD_NOTE[FREE_CHLORINE]}. {TEXAS_EXAMPLE[False]} Set your own "
                                         f"target here.")

    st.header("2. Your samples")
    demo = st.toggle("Demo: simulate a hidden truth and draw samples from it", value=True,
                     help="Off = upload your grab log and tap map, or type your own grab samples in the table."
                     if UPLOADS else cloud.DEMO_HELP_NO_UPLOADS)
    log_up = tap_up = inlet_up = None
    glog, glog_error = None, None
    if not demo and not UPLOADS:      # task 20: the hosted demo, examples only: the samples are typed
        st.caption(cloud.SAMPLES_OFF)
    elif not demo:     # task 16: the grab log and its tap map, after the model uploader (file_uploader[0] is unchanged)
        log_up = st.file_uploader("Grab log (CSV or Excel)", type=None, help=LOG_HELP, max_upload_size=UPLOAD_MB)
        tap_up = st.file_uploader("Tap map (CSV)", type=None, help=TAP_HELP, max_upload_size=UPLOAD_MB)
        # task 16c: a zone's inlets with different residuals; behind a checkbox, so the uploaders above stay the only
        # ones until it is ticked
        if st.checkbox(INLET_BOX, value=False, help=INLET_BOX_HELP):
            if ca:      # the review of task 16c: not built for the chloramine mode (its dose axis is the plant's)
                st.caption(INLET_NOT_CA)
            else:
                inlet_up = st.file_uploader(INLET_LABEL, type=None, help=INLET_HELP, max_upload_size=UPLOAD_MB)
    for f in (log_up, tap_up, inlet_up) if CLOUD else ():   # task 20: the hosted demo's file limit, before it is read
        if f is not None and cloud.too_big(f.size):
            intake_box.error(cloud.too_big(f.size))
            stop_page()
    if log_up is not None:
        log_bytes = log_up.getvalue()
        log_key = hashlib.sha1(log_bytes).hexdigest()[:10]
        species = "chloramine" if ca else "free_chlorine"
        header, guess = run_cached(log_header, CLOUD, file_key(log_bytes, log_up.name, species), log_bytes, log_up.name,
                                   species)
        if header is None:
            glog_error = guess
        else:
            def pick_column(label: str, field: str, extra: list[str]) -> str | None:
                opts = extra + header
                g = guess.get(field)
                if field == "reading" and guess.get("both_species"):
                    g = None        # a free and a total column: the operator chooses (the review of task 16)
                v = st.selectbox(label, opts, index=opts.index(g) if g in opts else 0, key=f"log_{field}_{log_key}")
                return None if v in (CHOOSE, IN_DATE, NO_ANALYTE, NO_UNIT) else v
            # open while a guess is not a template name, or the log needs a choice (both species, an analyte column)
            check_them = bool(guess.get("both_species") or guess.get("analyte")
                              or any(guess.get(k) not in TEMPLATE_NAMES for k in ("tap", "reading", "date", "time")))
            with st.expander("Columns of your grab log", expanded=check_them):
                st.caption("ResidualMap's guesses from the header: check them, and change any that is wrong.")
                mapping = {"tap": pick_column("Tap column", "tap", [CHOOSE]),
                           "reading": pick_column("Reading column", "reading", [CHOOSE]),
                           "date": pick_column("Date column", "date", [CHOOSE]),
                           "time": pick_column("Time column", "time", [IN_DATE]),
                           "analyte": pick_column("Analyte column (a log with one row per analyte)", "analyte",
                                                  [NO_ANALYTE])}
                if mapping["analyte"]:     # the final verification of task 16: which analyte's rows are read
                    v = st.selectbox("Analyte value (the analyte whose rows are the reading)",
                                     [ANALYTE_AUTO, ANALYTE_EVERY] + run_cached(
                                         column_values, CLOUD, file_key(log_bytes, log_up.name, mapping["analyte"]),
                                         log_bytes, log_up.name, mapping["analyte"]),
                                     key=f"log_analyte_value_{log_key}_{mapping['analyte']}")
                    mapping["analyte_value"] = None if v == ANALYTE_AUTO else "*" if v == ANALYTE_EVERY else v
                mapping["unit"] = pick_column("Units column", "unit", [NO_UNIT])
                date_format = DATE_FORMAT_CHOICES[st.selectbox("Date format", list(DATE_FORMAT_CHOICES),
                                                               key=f"log_date_format_{log_key}")]
                if date_format == "typed":
                    typed = st.text_input("Date format, as strptime writes it (%d day, %b month name, %m month "
                                          "number, %y or %Y year, %H:%M time)", value="%d-%b-%y",
                                          key=f"log_date_typed_{log_key}")
                    date_format = typed if "%" in typed else None
                if tap_up is not None:
                    args = (log_bytes, log_up.name, tap_up.getvalue(), tap_up.name, species, tuple(mapping.items()),
                            date_format, clock_h)
                    glog, glog_error = run_cached(grab_log, CLOUD, file_key(*args, ANALYZER_HOUR), *args, ANALYZER_HOUR)
                    if glog is not None and glog.analyzer:
                        if st.checkbox(KEEP_ALL, value=False, key=f"log_keep_all_{log_key}"):
                            glog, glog_error = run_cached(grab_log, CLOUD, file_key(*args, ANALYZER_HOUR, False), *args,
                                                          ANALYZER_HOUR, False)
                        else:
                            a_hour = st.number_input("Online analyzer: keep the reading nearest this hour (0 to 23)", 0,
                                                     23, int(ANALYZER_HOUR), key=f"log_analyzer_{log_key}")
                            if float(a_hour) != ANALYZER_HOUR:
                                glog, glog_error = run_cached(grab_log, CLOUD, file_key(*args, float(a_hour)), *args,
                                                              float(a_hour))
    if demo:
        n_demo = st.slider("Demo samples (daytime, random junctions)", 0, 15, 8)
        demo_seed = st.number_input("Demo scenario", 0, 99, 0)
        demo_kb = 0.40 if os.path.basename(net_path) != "Net2" else 0.10
        demo_kw = 0.70 if os.path.basename(net_path) != "Net2" else 0.20
    if ca:
        # the chloramine mode's prior centre: the plant's logged pH and chlorine-to-ammonia ratio (in demo mode, the
        # hidden truth's own, which an operator would log)
        ph0, r0 = 8.0, 4.5
        if demo:
            ch0 = CA.truth_chemistry(hidden_chem_draws("chloramine", int(demo_seed)))
            ph0, r0 = round(ch0["pH"], 2), round(ch0["cl2n"], 2)
        ph_log = st.number_input("Plant water pH (logged)", 6.5, 9.5, float(ph0), 0.05,
                                 help="Sets the prior centre of the bulk decay rate (EPA's chloramine model at 20 C).")
        cl2n_log = st.number_input("Chlorine to ammonia-N ratio, by mass (logged)", 3.0, 5.0, float(r0), 0.1)
        temp_log = st.number_input("Water temperature this month (C)", 0.0, 35.0, 20.0, 0.5,
                                   help="Used only by the nitrification watch (15 C or warmer).")

    st.header("3. Route")
    K = st.slider("Sites on next month's route", 3, 12, 6)


# ----------------------------------------------------------------------------- heavy lifting, cached
@contextlib.contextmanager
def workers_skip_this_app():
    """Process pools (the grid build, the chloramine demo truth) start their workers with multiprocessing's spawn method,
    which re-runs the main module's file in every worker before its job.  Under Streamlit the main module is this app,
    so each worker reran the whole app in bare mode first, and the grid build's parallel workers collided on EPANET's
    temporary files and crashed (a fresh grid for an uploaded .inp failed this way; found by task 14's AppTest of the
    upload path).  While a pool starts, the main module is an empty stand-in, so workers import only residualmap.  The
    swap is process-wide, like Streamlit's own install of each script run as __main__."""
    main = sys.modules.get("__main__")
    sys.modules["__main__"] = types.ModuleType("__main__")
    try:
        yield
    finally:
        sys.modules["__main__"] = main


@st.cache_resource(show_spinner=False)
def prepare(path: str, dose: float, chloramine: bool = False, chem_tag: str = "", grid: str = "full",
            inlets: tuple = ()):
    """Nominal model, physics features and the simulator grid (cached on disk too): the 675-run free-chlorine grid
    (or, for an uploaded file with one C value on every pipe, its 225-run gamma-0 slice 'full_fixed_gamma', task 15b),
    or the chloramine mode's 1440-run grid under its own cache tag.  The in-memory cache key is (file, dose,
    disinfectant, chemistry tag, grid): chem_tag is the chloramine grid's chemistry.cache_tag ('' for free chlorine),
    so a change to the chloramine grid's definition or physics version is never served a grid from before it.
    inlets (task 16c): the per-inlet dose ratios as sorted (node, ratio) pairs, () for none (today's grid and key)."""
    sc = nominal_scenario(path, sample_hour=14, source_dose=dose)
    X = build_features(sc)
    kw = {"inlet_ratios": dict(inlets)} if inlets else {}
    with workers_skip_this_app():      # n_jobs: None (the default) in local mode; 2 workers in cloud mode (task 20)
        if chloramine:
            simulator_grid_24h(sc, CACHE_DIR, "chloramine", n_jobs=cloud.grid_workers(APP_ROOT), cond=CA.ca_condition(),
                               **kw)
        else:
            simulator_grid_24h(sc, CACHE_DIR, grid, n_jobs=cloud.grid_workers(APP_ROOT), **kw)
    return sc, X


def prepare_in_session(path: str, dose: float, chloramine: bool, grid: str, inlets: tuple = ()):
    """prepare for a visitor's upload in cloud mode (task 20): the same nominal model, features and grid, the grid
    built in memory (build_grid_24h) and never written to outputs/cache.  Returns (sc, X, (params, Z)); the caller
    keeps it in the browser session only.  Since the final verification of task 20 it runs in this run's folder
    (in_folder, cloud.upload_grid), its grid's workers too.  inlets (task 16c): as prepare's."""
    return in_folder(cloud.upload_grid, path, dose, chloramine, grid, cloud.grid_workers(APP_ROOT), dict(inlets) or None)


@st.cache_resource(show_spinner=False)
def water_age(path: str, _sc):
    """Water age of the operator's model at the grid's 9 demand x roughness settings, and the nominal model's share
    of the simulation's starting water (10 EPANET AGE runs, once per network)."""
    return hydraulic_age_band(_sc)


@st.cache_resource(show_spinner=False)
def chlorine_loss(path: str, dose: float, member: tuple, dose_mult: float, _model, chloramine: bool = False,
                  inlets: tuple = ()):
    """Where the calibrated member loses its chlorine: four EPANET runs (as calibrated, wall decay off, bulk decay
    off, no decay), cached per network, dose and member."""
    return loss_split(_model)


def age_test(path: str) -> dict | None:
    """The simulated test of an example network under the default truth (outputs/chem/water_age_<net>.json):
    the age range's label, its coverage of the true daily-mean and daily-max age, and how far the calibrated
    loss split was from the truth's.  An uploaded file has no such test: None."""
    try:
        with open(AGE_RESULTS.format(os.path.basename(path))) as fh:
            d = json.load(fh)
        t = d["by_truth"]["default"]
        return {"label": d["band_label"], "cov": t["band_coverage"], "cov_max": t["band_coverage_daily_max"],
                "split_median_diff": t["loss_split"]["median_abs_diff_mean"],
                "split_median_diff_max": t["loss_split"]["median_abs_diff_max"],
                "split_mae": t["loss_split"]["wall_share_mae"], "n_seeds": t["n_seeds"]}
    except (OSError, KeyError, ValueError, TypeError):
        return None


@st.cache_resource(show_spinner=False)
def demo_truth(path: str, dose: float, seed: int, kb: float, kw: float):
    return build_scenario(path, seed=seed, sample_hour=14, source_dose=dose, kb_per_day=kb, kw_m_per_day=kw)


@st.cache_resource(show_spinner=False)
def demo_truth_ca(path: str, dose: float, seed: int):
    """A hidden chloramine truth: EPA's chloramine chemistry in EPANET-MSX (about 15 s on Net3 or Net2 with compiled
    reactions, minutes without a C compiler), with the wall rate the task-12 experiment set for the example networks
    (0.20 m/day elsewhere).  It runs in a separate process: MSX changes the working directory, which is process-wide,
    and this app serves every browser session from threads of one process.  Returns (truth, None), or (None, the error)
    when EPANET-MSX cannot run here (no OpenMP runtime, a broken MSX library, a failed worker).  The failure is cached
    like a result, so MSX is tried once per network, dose and scenario in this server process: a widget change does
    not start a new worker each time, and the demo's hidden truth cannot switch between the MSX and the reduced truth
    within a session."""
    kw = CA.KW_REF.get(os.path.basename(path), 0.20)
    try:
        with workers_skip_this_app():
            return CA.truth_in_subprocess(path, seed, kw, dose=dose), None
    except Exception as e:  # noqa: BLE001
        return None, f"{type(e).__name__}: {str(e)[:240]}"


@st.cache_resource(show_spinner=False)
def demo_truth_ca_reduced(path: str, dose: float, seed: int):
    """The fallback when EPANET-MSX cannot run here (no OpenMP runtime, a broken library): the same hidden network and
    draws with the reduced truth, EPANET first-order decay at the rate EPA's chemistry gives for this water (task 12's
    first-order twin, chloramine.truth_scenario(kinetics='first')).  It runs in this process and never imports
    EPANET-MSX's library; like demo_truth and prepare, it writes EPANET's temp.inp, .rpt and .bin to the app's working
    directory (no change of working directory).  It shares the model's own kinetics, so the demo's scores flatter the
    model, and the app says so."""
    kw = CA.KW_REF.get(os.path.basename(path), 0.20)
    return CA.truth_scenario(path, seed, kw, kinetics="first", dose=dose)


def chloramine_status() -> dict | None:
    """The chloramine mode's simulated test (outputs/chloramine/summary_chloramine_<net>.json): experimental when the
    mode's model missed a pre-registered bar on either network."""
    out = {}
    for net in ("Net3", "Net2"):
        try:
            with open(CA_RESULTS.format(net)) as fh:
                out[net] = json.load(fh)["acceptance"]
        except (OSError, KeyError, ValueError):
            return None
    return {"experimental": any(a.get("experimental") for a in out.values()), "acceptance": out}


try:
    wn_check = wntr.network.WaterNetworkModel(net_path if os.path.exists(net_path) else os.path.join(LIB, f"{net_path}.inp"))
except Exception as e:  # noqa: BLE001
    st.error(f"Could not read this .inp with EPANET/WNTR: {e}")
    stop_page()
if not source_nodes(wn_check):
    st.error("No reservoir, inflow junction or tank found — the model needs a place where water enters.")
    stop_page()

# task 16c: a zone's inlet residuals (the inlet file, and taps the tap map marks as a zone inlet: the median of their
# readings over the grab log's window), checked against the model's entry points.  None: every inlet at the one dose.
# Its review: every tap the tap map marks as a zone inlet must give its inlet a residual in that window (refused in
# plain words otherwise, where it was left out without a word), its readings outside the window are listed as not
# used, and the inlet file's notes (decimal commas, µg/L) are shown
inlet_ratios, inlet_info, inlet_from, inlet_notes, inlet_unused = None, None, {}, [], []
file_res, tap_res = None, {}
if inlet_up is not None:
    got, inlet_error = run_cached(inlet_file, CLOUD, file_key(inlet_up.getvalue(), inlet_up.name, "free_chlorine"),
                                  inlet_up.getvalue(), inlet_up.name, "free_chlorine")
    if inlet_error:
        st.error(f"The inlet residuals cannot be used as they are. {inlet_error}")
        stop_page()
    file_res, inlet_notes = got
if glog is not None and glog.inlet_taps and tap_up is not None:
    if ca:
        n_ = len(glog.inlet_taps)
        st.error(f"The tap map marks tap{'s' if n_ != 1 else ''} {', '.join(sorted(glog.inlet_taps)[:6])}"
                 f"{f' and {n_ - 6} more' if n_ > 6 else ''} as a zone inlet. {INLET_NOT_CA_WHY} Clear the role in the "
                 f"tap map (the readings are then fitted as a tap's), or leave those readings out of the grab log.")
        stop_page()
    g_in = glog.only_junctions(wn_check.junction_name_list)
    in_months = g_in.window(LOG_WINDOW_MONTHS)[2]["months"] or None       # no grab rows: every month of the log
    try:
        tap_res = inlet_tap_residuals(g_in, source_nodes(wn_check), in_months, file_res)
    except LogRefused as e:
        st.error(f"The inlet residuals cannot be used as they are. {e}")
        stop_page()
    inlet_unused = g_in.inlet_unused(in_months)
if file_res or tap_res:
    try:
        inlet_from = combine_inlet_residuals(file_res, tap_res)
        inlet_ratios, inlet_info = inlet_ratios_from(inlet_from, source_nodes(wn_check))
    except LogRefused as e:
        st.error(f"The inlet residuals cannot be used as they are. {e}")
        stop_page()
inlets_key = tuple(inlet_ratios.items()) if inlet_ratios is not None else ()
dose = draw_dose(inlet_info)      # the sidebar's dose input, in its place before the threshold (task 16c)

if ca and demo and os.path.basename(net_path) == "ky4":
    st.error("The chloramine demo truth (EPA's chloramine chemistry in EPANET-MSX) is offered on Net3 and Net2 only: on "
             "ky4 one simulation takes too long here. Turn the demo off to use your own total chlorine samples.")
    st.stop()
grid = "chloramine" if ca else ("full_fixed_gamma" if uniform_c else "full")
n_runs = int(np.prod([len(a) for a in GRIDS[grid]]))
if up is not None:     # task 15b: refuse a grid build that would leave less than 1.5 GB free (a cached grid needs none)
    cond_ = CA.ca_condition() if ca else None
    if not os.path.exists(grid_cache_path(types.SimpleNamespace(wn_name=net_path, source_dose=float(dose)), CACHE_DIR,
                                          grid, cond_, inlet_ratios)):
        disk_ok, disk_msg = size_gate(readiness.run_cost, CACHE_DIR)
        if not disk_ok:
            st.error(disk_msg)
            st.stop()
# task 20: a visitor's upload, kept in this browser session only; task 16c: so is a grid made from their inlet residuals
SESSION_ONLY = CLOUD and (up is not None or inlet_ratios is not None)
grid_data = None                             # the session's own grid (a cloud upload); None: the cached grid
with st.spinner(f"Running {n_runs} EPANET simulations of your model over the decay and hydraulic-mismatch grid "
                f"(once per network; cached afterwards)…"):
    if SESSION_ONLY:
        sc, X, grid_data = cloud.session_cached(st.session_state, "prepare", (net_path, float(dose), ca, grid, inlets_key),
                                                prepare_in_session, net_path, float(dose), ca, grid, inlets_key)
    elif inlets_key:     # task 16c: the grid at the inlets' ratios, under its own cache name
        sc, X = prepare(net_path, float(dose), ca, cache_tag(CA.ca_condition(), "chloramine") if ca else "", grid,
                        inlets_key)
    else:
        sc, X = (prepare(net_path, float(dose), True, cache_tag(CA.ca_condition(), "chloramine")) if ca else
                 prepare(net_path, float(dose)) if grid == "full" else prepare(net_path, float(dose), False, "", grid))
with st.spinner("Running 10 EPANET water-age simulations of your model (demand and roughness settings; once per network)…"):
    band = run_cached(water_age, SESSION_ONLY, (net_path,), net_path, sc, job=(hydraulic_age_band, sc))
junctions = list(sc.junctions)
# task 16: hours on the page.  An uploaded file whose START CLOCKTIME is not midnight shows times of day (route hours,
# worst hour, hour picker, hour of the daily minimum); every bundled example keeps model hours (Net2's with a caption)
use_clock = shows_clock(clock_h)
NET2_CAPTION = "model hour (this file starts at 08:00)"
is_net2 = up is None and os.path.basename(net_path) == "Net2"


def hlabel(h: int) -> str:
    """The 'HH:00' the page prints for model hour h: a time of day for an uploaded file that does not start at
    midnight (importers.clock_label), the model hour as before otherwise."""
    return clock_label(h, clock_h) if use_clock else f"{h:02d}:00"


def htag(h: int) -> str:
    """The short hour tag next to a sample on the maps ('9h' as before; a time of day for a converted file)."""
    return f"{clock_label(h, clock_h)}" if use_clock else f"{h}h"


# samples table -----------------------------------------------------------------------------
if ca:
    status = chloramine_status()
    exp_note = (" **Experimental:** in simulation the mode missed at least one of its pre-registered bars (README and "
                "journal, task 12)." if status is None or status["experimental"] else "")
    st.info(f"Chloramine mode: every number is TOTAL chlorine (mg/L as Cl2), with the chloramine mode's own decay ranges, "
            f"dose axis and prior (plant pH {ph_log:g}, Cl2:N {cl2n_log:g}); minimum residual {threshold:g} mg/L total "
            f"chlorine ({'a common utility operating target, not a California rule' if abs(threshold - 0.5) < 1e-9 else 'set by you'}; "
            f"California requires a detectable residual). Tested in simulation only.{exp_note}")
if demo:
    if ca and CLOUD:     # task 20: the hosted demo goes straight to the reduced truth, with its warning
        with st.spinner("Simulating a hidden chloraminated network with EPANET (the reduced chloramine truth)…"):
            tr = run_cached(demo_truth_ca_reduced, SESSION_ONLY, (net_path, float(dose), int(demo_seed)), net_path,
                            float(dose), int(demo_seed),
                            job=(cloud.upload_truth_reduced, net_path, float(dose), int(demo_seed)))
        st.warning(cloud.REDUCED_TRUTH)
    elif ca:
        from residualmap.msx import compiler_available
        how_long = ("about 15 s" if compiler_available() else
                    "several minutes: no C compiler was found, so the reactions run uncompiled")
        with st.spinner(f"Simulating a hidden chloraminated network with EPA's chloramine chemistry in EPANET-MSX ({how_long})…"):
            tr, msx_error = demo_truth_ca(net_path, float(dose), int(demo_seed))
        if msx_error:
            with st.spinner("EPANET-MSX could not run here: simulating the reduced chloramine truth with EPANET instead…"):
                tr = demo_truth_ca_reduced(net_path, float(dose), int(demo_seed))
            st.warning(f"Reduced demo truth. EPA's chloramine chemistry could not run here in EPANET-MSX "
                       f"({msx_error}). The hidden truth is its reduced form instead: EPANET "
                       f"first-order decay at the rate EPA's chemistry gives for this water, on the same hidden network "
                       f"(task 12's first-order twin). It has the model's own kinetics, so the demo's scores flatter the "
                       f"model.")
    else:
        tr = run_cached(demo_truth, SESSION_ONLY, (net_path, float(dose), int(demo_seed), demo_kb, demo_kw), net_path,
                        float(dose), int(demo_seed), demo_kb, demo_kw,
                        job=(cloud.upload_truth, net_path, float(dose), int(demo_seed), demo_kb, demo_kw))
    rng = np.random.default_rng(int(demo_seed))
    n_demo = min(int(n_demo), len(junctions))          # task 15b: a network smaller than the slider (Net1 has 9)
    js = list(rng.choice(junctions, int(n_demo), replace=False)) if n_demo else []
    # daytime by the clock: for an uploaded file that does not start at midnight, the model hours of 07:00 to 17:00
    # (the review of task 16); DAY_HOURS for every bundled example, so the default demo draws exactly as before
    hs = [int(h) for h in rng.choice(daytime_model_hours(clock_h) if use_clock else DAY_HOURS, int(n_demo))] if n_demo else []
    ys = [float(np.clip(tr.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None)) for j, h in zip(js, hs)]
    default = pd.DataFrame({"junction": js, "hour": hs, "mg/L": np.round(ys, 2)})
else:
    tr = None
    default = pd.DataFrame({"junction": pd.Series(dtype=str), "hour": pd.Series(dtype=int), "mg/L": pd.Series(dtype=float)})

# task 16: an uploaded grab log prefills the table (rows of the last LOG_WINDOW_MONTHS months, model hours)
log_key_suffix, below_detection, above_range = "", [], []
if log_up is not None:
    if glog_error:
        st.error(f"The grab log cannot be used as it is. {glog_error}")
        st.stop()
    if tap_up is None:
        st.info("Upload the tap map too (tap ID and junction ID, under '2. Your samples'): the grab log's taps are "
                "matched to the model's junctions through it.")
    elif glog is not None:
        gj = glog.only_junctions(junctions)
        rows_w, cens_w, winfo = gj.window(LOG_WINDOW_MONTHS)
        if winfo["n_months_in_log"] > LOG_WINDOW_MONTHS:
            span = (f"the last {LOG_WINDOW_MONTHS} months of the log ({', '.join(winfo['months'])}; the log runs from "
                    f"{winfo['log_first']} to {winfo['log_last']})")
        else:
            span = f"the whole log ({winfo['n_months_in_log']} month{'s' if winfo['n_months_in_log'] != 1 else ''})"
        outside = (f" {winfo['n_outside']} row{'s' if winfo['n_outside'] != 1 else ''} dated before "
                   f"{winfo['months'][0]} {'are' if winfo['n_outside'] != 1 else 'is'} outside the window and not used "
                   f"(listed below)." if winfo["n_outside"] else "")
        if winfo["first"] is None:
            st.warning("The grab log has no reading the model can use: every row is listed below with its reason. "
                       + " ".join(gj.notes))
        else:
            st.info(f"Grab log: using {winfo['n_used']} samples from {winfo['first']} to {winfo['last']}, {span}.{outside} "
                    f"Within the window every sample is treated as one simulated day: a reading is compared with the "
                    f"model's day at its junction and hour, whatever its date. " + " ".join(gj.notes))
        not_used = (gj.dropped + [{"row": int(r.row), "reason": r.reason} for r in cens_w.itertuples()]
                    + winfo["outside"] + inlet_unused)        # task 16c's review: zone inlet readings outside the window
        if not_used:
            st.warning(f"Rows of the grab log that are not used ({len(not_used)}), each with its reason:\n\n"
                       + "\n".join(f"- {ln}" for ln in dropped_lines(sorted(not_used, key=lambda d: d["row"]))))
        for w in gj.warnings:
            st.warning(w)
        if len(cens_w):      # below detection and above the range are marked apart (the review of task 16)
            below_detection = sorted(set(cens_w.junction[cens_w.kind == "below"]))
            above_range = sorted(set(cens_w.junction[cens_w.kind == "above"]))
        default = pd.DataFrame({"junction": rows_w.junction.astype(str), "hour": rows_w.hour.astype(int),
                                "mg/L": rows_w.y.astype(float)})
        log_key_suffix = "_log" + hashlib.sha1(default.to_csv(index=False).encode()).hexdigest()[:10]

if inlet_info is not None:      # task 16c: the inlet residuals in use, where each comes from, and what they do not cover
    st.info(f"{inlet_info['sentence']} From: " + "; ".join(f"{n}, {d['source']}" for n, d in inlet_from.items())
            + f". {INLET_ZONE_NOTE}" + "".join(f" {n_}" for n_ in inlet_notes))

st.subheader("Grab samples")
st.caption("One row per sample: the junction ID from your model, the hour it was taken (0 to 23), the TOTAL chlorine "
           "reading (chloramine)." if ca else
           "One row per sample: the junction ID from your model, the hour it was taken (0–23), the free chlorine reading.")
if use_clock:
    st.caption(f"The hour column is the MODEL hour: this file starts at {clock_label(0, clock_h)} (START CLOCKTIME), so a "
               f"sample taken at {clock_label(3, clock_h)} is hour 3. An uploaded grab log's times are converted for you; "
               f"a typed hour is a model hour. Everywhere else on this page, hours are times of day.")
elif is_net2:
    st.caption(f"Net2: the hour column, and every hour on this page, is a {NET2_CAPTION}."
               + (" Your log's times of day are read as model hours here, not converted (the pilot with --inp Net2 "
                  "subtracts 8 h)." if log_up is not None else ""))
samples = st.data_editor(default, num_rows="dynamic", width='stretch',
                         key=f"samples_{net_path}_{demo}_{demo and (n_demo, demo_seed)}" + ("_ca" if ca else "")
                             + log_key_suffix,
                         column_config={"junction": st.column_config.SelectboxColumn("junction", options=junctions, required=True),
                                        "hour": st.column_config.NumberColumn("hour", min_value=0, max_value=23, step=1),
                                        "mg/L": st.column_config.NumberColumn("mg/L", min_value=0.0, max_value=5.0, step=0.01, format="%.2f")})
samples = samples.dropna(subset=["junction", "hour", "mg/L"])
samples = samples[samples.junction.isin(junctions)]
S = pd.DataFrame({"junction": samples.junction.astype(str), "hour": samples.hour.astype(int), "y": samples["mg/L"].astype(float)})

# model ------------------------------------------------------------------------------------
if ca:
    model = SimGP24(sc, X, seed=0, cache_dir=CACHE_DIR, threshold=float(threshold), grid="chloramine",
                    cond=CA.ca_condition(), lik_sd=LIK_SD_CA, doses=DOSES_CA, grid_data=grid_data,
                    inlet_ratios=inlet_ratios)
    model.log_prior = CA.prior_log_vector(model.params, float(ph_log), float(cl2n_log))
elif grid == "full":
    model = SimGP24(sc, X, seed=0, cache_dir=CACHE_DIR, threshold=float(threshold), grid_data=grid_data,
                    inlet_ratios=inlet_ratios)
else:
    model = SimGP24(sc, X, seed=0, cache_dir=CACHE_DIR, threshold=float(threshold), grid=grid, grid_data=grid_data,
                    inlet_ratios=inlet_ratios)
model.fit(S if len(S) else None)
hourly = model.predict_hours()
pmin = model.predict_daily_min()
if use_clock:     # task 16: the route's daytime window, 07:00 to 17:00, by the clock of this file
    route = plan_route(model, int(K), hours=daytime_model_hours(clock_h), threshold=float(threshold),
                       exclude=list(S.junction), with_age=True, age_initial_share=band.initial_share, hour_label=hlabel)
else:
    route = plan_route(model, int(K), threshold=float(threshold), exclude=list(S.junction), with_age=True,
                       age_initial_share=band.initial_share)

# ----------------------------------------------------------------------------- headline numbers
hour_options = {hlabel(h): h for h in sorted(range(24), key=lambda h: clock_minute_of(h, clock_h) if use_clock else h)}
view = st.radio("Show", ["Daily minimum (the compliance number)"] + list(hour_options), horizontal=True, index=0)
if view.startswith("Daily"):
    med, lo, hi = pmin["median"], pmin["lo90"], pmin["hi90"]
    p_below = pmin["p_below"]          # the model's own draws, at the threshold it was built with
    label = "daily minimum"
else:
    h = hour_options[view]; z_mu, z_sd = hourly[0][h], hourly[1][h]
    med = pd.Series(np.exp(z_mu), index=junctions); lo = pd.Series(np.exp(z_mu - 1.645 * z_sd), index=junctions)
    hi = pd.Series(np.exp(z_mu + 1.645 * z_sd), index=junctions)
    p_below = pd.Series(norm.cdf((np.log(float(threshold)) - z_mu) / z_sd), index=junctions)
    label = view
flag = p_below > 0.5
frac_by_hour = pd.Series([float((hourly[0][h] < np.log(float(threshold))).mean()) for h in range(24)], index=range(24))
worst = int(frac_by_hour.idxmax())

KB_FMT = ".3g" if ca else ".2f"     # the chloramine grid's bulk rates go down to 0.0025 per day
c1, c2, c3, c4 = st.columns(4)
c1.metric("Samples used", len(S))
c2.metric(f"Junctions likely below {threshold:g} mg/L ({label})", f"{int(flag.sum())} of {len(junctions)}")
c3.metric("Worst hour of the day", hlabel(worst), f"{frac_by_hour[worst]:.0%} of junctions below by median")
if model.map_params_ is not None:
    kb, kw, g, dm, rm = model.map_params_
    c4.metric("Calibrated decay (bulk / wall)", f"{kb:{KB_FMT}} /d · {kw:.2f} m/d", f"{old_pipe_text(g, uniform_c)}, demand ×{dm:.2f}, dose ×{model.map_dose_:.2f}")
else:
    c4.metric("Calibrated decay", "no samples yet", "map = your model's physics alone")
if is_net2:
    st.caption(f"Net2: every hour on this page is a {NET2_CAPTION}, and its demand patterns are 55 h long, so an hour "
               f"of its day is an hour of that cycle{'' if CLOUD else ' (README, Net2)'}.")     # task 20: not in cloud
elif use_clock:
    st.caption(f"Hours are times of day: this file starts at {clock_label(0, clock_h)} (START CLOCKTIME), and every model "
               f"hour is shown as its time of day.")


def edge_warning(m) -> str | None:
    """capabilities.grid_edge_warning for the fitted model: a plain warning when more than half of the calibration
    (capabilities.EDGE_WARN) sits on the lowest or highest bulk or wall decay rate of the grid (simgp.grid_edge_mass).
    None with no samples."""
    if m.map_params_ is None:
        return None
    return grid_edge_warning(grid_edge_mass(m.W_, m.params, m.doses), m.params, ca)


edge_msg = edge_warning(model)
if edge_msg:
    st.warning(edge_msg)


# ----------------------------------------------------------------------------- the four panels
def network_panel(ax, series, title, cmap, vrange, marks=None, mark_labels=None, colorbar_label=None, below=None,
                  above=None):
    wntr.graphics.plot_network(sc.wn, node_attribute=series.to_dict(), node_size=max(12, int(4000 / len(junctions))),
                               node_cmap=cmap, node_range=vrange, ax=ax, link_width=0.6, add_colorbar=True, title=title)
    if below:     # task 16: junctions with a reading below detection ('<0.02', 'ND'): shown, not fitted
        ax.scatter([sc.coords[j][0] for j in below], [sc.coords[j][1] for j in below], s=90, marker="x", c="red",
                   linewidths=1.8, zorder=6, label="measured below detection (not fitted)")
    if above:     # a reading above the measuring range ('>2.2'): high chlorine, never called below detection
        ax.scatter([sc.coords[j][0] for j in above], [sc.coords[j][1] for j in above], s=90, marker="+", c="blue",
                   linewidths=1.8, zorder=6, label="measured above the range (not fitted)")
    if below or above:
        ax.legend(loc="lower left", fontsize=8)
    if marks:
        ax.scatter([sc.coords[j][0] for j in marks], [sc.coords[j][1] for j in marks], s=130, facecolors="none",
                   edgecolors="cyan", linewidths=1.8, zorder=5)
        if mark_labels:
            for j, lab in zip(marks, mark_labels):
                ax.annotate(lab, sc.coords[j], fontsize=8, color="darkcyan", xytext=(4, 3), textcoords="offset points")


def route_tag(k: int, h: int) -> str:
    """'1: 09h', a route site's number and hour on the maps (a time of day for a converted file)."""
    return f"{k + 1}: {clock_label(h, clock_h)}" if use_clock else f"{k + 1}: {h:02d}h"


def worst_hour_bars(ax, ylabel: str, legend_label: str):
    """The share of junctions below the threshold at each hour of the day: by model hour as before, or (task 16) at
    each model hour's time of day for a converted file, with the sampling window by the clock."""
    xs = [clock_minute_of(h, clock_h) / 60 if use_clock else h for h in frac_by_hour.index]
    ax.bar(xs, frac_by_hour.values * 100, color=["tab:red" if x >= 18 or x < 6 else "tab:blue" for x in xs])
    ax.axvspan(6.5, 17.5, color="gold", alpha=0.15, label=legend_label)
    ax.set(xlabel=NET2_CAPTION if is_net2 else "hour of day", ylabel=ylabel)


def four_panels(figsize=(24, 6.5)):
    fig, axes = plt.subplots(1, 4, figsize=figsize)
    top = max(1.2, float(dose))
    title0 = (f"Total chlorine map: {label} (mg/L)\nfrom {len(S)} grab samples" if ca else
              f"Chlorine map — {label} (mg/L)\nfrom {len(S)} grab samples")
    network_panel(axes[0], med, title0, "viridis", (0, top),
                  list(S.junction), [htag(h) for h in S.hour], below=below_detection, above=above_range)
    network_panel(axes[1], hi - lo, "Uncertainty — width of the 90% band (mg/L)", "magma", (0, float((hi - lo).quantile(0.98))))
    network_panel(axes[2], p_below, f"P(residual < {threshold:g} mg/L) — {label}\n{int(flag.sum())} junctions likely below", "Reds", (0, 1))
    network_panel(axes[3], p_below, f"Sample here next: {len(route)} sites for next month's route", "Reds", (0, 1),
                  list(route.junction), [route_tag(k, h) for k, h in enumerate(route.hour)])
    fig.tight_layout()
    return fig


st.pyplot(four_panels(), width='stretch')

# ----------------------------------------------------------------------------- next samples with reasons
st.subheader("Sample here next")
st.dataframe(route.assign(hour=[hlabel(h) for h in route.hour])[["junction", "hour", "p_below", "water_age_h", "reason"]]
             .rename(columns={"p_below": f"P(daily min < {threshold:g})", "water_age_h": "water age (h)"}), width='stretch', hide_index=True)

# ----------------------------------------------------------------------------- every junction (task 16)
P_COL = f"P(daily min < {threshold:g})"
every = pd.DataFrame({"junction": junctions, "median (mg/L)": pmin["median"].round(3).values,
                      "lo90 (mg/L)": pmin["lo90"].round(3).values, "hi90 (mg/L)": pmin["hi90"].round(3).values,
                      P_COL: pmin["p_below"].round(2).values,
                      "hour of the daily minimum": [hlabel(int(h)) for h in np.argmin(hourly[0], axis=0)]})
st.subheader("Every junction: the daily minimum")
st.dataframe(every, width='stretch', hide_index=True)
n_top = min(20, len(junctions))
st.caption(f"Each junction's daily minimum (median and 90% band, mg/L), the probability it is below {threshold:g} mg/L "
           f"and the hour of the day its median is lowest. The table's toolbar (top right of the table, with the pointer "
           f"over it) downloads it as CSV. The {n_top} junctions most likely below {threshold:g} mg/L follow as a plain "
           f"table, which screen readers read more easily than the grid above.")
st.table(every.sort_values(P_COL, ascending=False, kind="stable").head(n_top).set_index("junction"))

if ca:
    st.subheader("Nitrification watch (chloramine)")
    watch = CA.watch_flags(model, sc, float(temp_log))
    st.caption(f"Literature thresholds, not validated: a junction is on the watch list when P(daily minimum total chlorine "
               f"< {CA.WATCH_CRIT_MGL:g} mg/L) > {CA.WATCH_P:g} (critical range 0.2 to 0.65 mg/L, Sathasivan, Fisher & Tam "
               f"2008), the water is {CA.WATCH_MIN_TEMP_C:g} C or warmer (this month: {float(temp_log):g} C; nitrification is "
               f"reported from 8 to 26 C, US EPA 2002) and its water age is in the oldest quarter of your model's. It marks "
               f"where to look (nitrite, free ammonia, HPC), not where nitrification is. {int(watch.sum())} of "
               f"{len(junctions)} junctions are on it.")
    if watch.any():
        st.dataframe(pd.DataFrame({"junction": list(watch.index[watch]),
                                   f"P(daily min < {CA.WATCH_CRIT_MGL:g})": model.p_below_mc(CA.WATCH_CRIT_MGL)[watch].round(2).values,
                                   "water age (h, daily mean)": sc.age_by_hour_h.mean()[watch].round(1).values}),
                     width='stretch', hide_index=True)

# ----------------------------------------------------------------------------- water age, and where chlorine is lost
age_max = band.nominal.max()
age_lo, age_hi = band.daily_range("max")
j_old, h_old = oldest_water(band.nominal)
test = age_test(net_path)
range_label = test["label"] if test else RANGE_LABEL
# every age is from a 7-day simulation that starts with the file's water in the pipes and tanks; where more than
# INITIAL_SHARE_MAX of the water at a junction is still that starting water, its age is a lower bound
lb_max, lb_mean = band.lower_bound("max"), band.lower_bound("mean")
share_old = float(band.initial_share[j_old].iloc[int(np.argmax(band.nominal[j_old].values))])
if lb_max[j_old]:
    oldest_line = (f"oldest water: at least {h_old:.0f} h at junction {j_old} (your model; at that hour {share_old:.0%} of "
                   f"the water there is still the 7-day simulation's starting water, so it is older than the simulation can show)")
else:
    oldest_line = (f"oldest water: {h_old:.0f} h at junction {j_old} (your model; {age_lo[j_old]:.0f} to {age_hi[j_old]:.0f} h "
                   f"across demand and roughness errors)")
lower_note = (f"Every age comes from a 7-day simulation that starts with your file's water in the pipes and tanks. Where "
              f"more than {INITIAL_SHARE_MAX:.0%} of a junction's water is still that starting water, its age is a lower bound: "
              f"{int(lb_max.sum())} of {len(junctions)} junctions at their oldest hour, {int(lb_mean.sum())} on the day's average.")
age_col, loss_col = st.columns([2, 1])
with age_col:
    st.subheader("Water age")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    network_panel(axes[0], age_max, "Water age at its oldest hour of the day (h)\nyour EPANET model, 7-day simulation", "YlOrBr",
                  (0, float(age_max.quantile(0.98))))
    network_panel(axes[1], age_hi - age_lo, f"How much that age moves (h): {range_label}\n9 demand and roughness settings",
                  "magma", (0, max(float((age_hi - age_lo).quantile(0.98)), 1.0)))
    fig.tight_layout()
    st.pyplot(fig, width='stretch')
    if test and range_label == RANGE_LABEL:
        tested = (f" In simulation on this network it held the true daily-maximum age at {test['cov_max']:.0%} of junctions "
                  f"(daily mean: {test['cov']:.0%}), so it is shown as a range, not a 90% band: errors in single junctions' "
                  f"demands are not among the 9 settings.")
    elif test:
        tested = f" In simulation on this network it held the true daily-maximum age at {test['cov_max']:.0%} of junctions."
    else:
        tested = " It has not been tested against a simulated truth for this file, so it is shown as a range, not a 90% band."
    st.caption(f"Hours since the water left the source, from your EPANET file: {oldest_line}. {lower_note} Right: how far "
               f"each junction's age moves across the 9 settings the model already weighs (demand x0.85, x1, x1.15; pipe "
               f"roughness x0.9, x1, x1.1).{tested}"
               + (" Water age starts at 0 at every inlet: in a zone cut from a larger system the water entering it is "
                  "already old, so the zone's water is older than this map shows." if inlet_info is not None else ""))
with loss_col:
    st.subheader("Where chlorine is lost")
    if model.map_params_ is None:
        st.info("Enter at least one grab sample: the split uses the decay rates your samples calibrate.")
    else:
        member = (net_path, float(dose), tuple(model.map_params_), float(model.map_dose_))
        split = run_cached(chlorine_loss, SESSION_ONLY, (*member, ca, inlets_key), *member, model, ca, inlets_key,
                           job=(loss_split, model))
        ws = split.wall_share.dropna()
        if ws.empty:
            st.info("Chlorine loss is under 2% at every junction for these decay rates: there is nothing to split.")
        else:
            fig, ax = plt.subplots(figsize=(7, 5.6))
            network_panel(ax, ws, "Share of the chlorine loss at pipe walls\n0 = all in the water, 1 = all at the walls", "RdYlBu_r", (0, 1))
            fig.tight_layout()
            st.pyplot(fig, width='stretch')
            na = split.nonadditivity.dropna()
            n_small = len(junctions) - len(ws)
            accuracy = (f" In simulation on this network ({test['n_seeds']} scenarios, 8 samples each) the calibrated median was "
                        f"off from the truth's by {test['split_median_diff']:.2f} on average and by up to {test['split_median_diff_max']:.2f} "
                        f"in one scenario." if test else " It has not been tested against a simulated truth for this file.")
            gamma_words = ("old-pipe factor not identifiable from one C value" if uniform_c else
                           f"old-pipe factor {model.map_params_[2]:.1f}")
            st.caption(f"For the single most likely decay rates your samples calibrate (bulk {model.map_params_[0]:{KB_FMT}} /day, wall "
                       f"{model.map_params_[1]:.2f} m/day, {gamma_words}): a median {ws.median():.0%} of the "
                       f"chlorine lost on the way to a junction is lost at the pipe walls, the rest in the water (organics and other "
                       f"reactants). It is an estimate: with few samples, bulk and wall decay can trade off against each other."
                       f"{accuracy} Exact along a single path from the source; approximate where flows mix (here the two parts add "
                       f"up to {na.min():.2f} to {na.max():.2f} of the total). Not coloured (under 2% loss): {n_small} "
                       f"junction{'' if n_small == 1 else 's'}. Field studies report wall loss up to 97% of the total in old pipes "
                       f"(Maleki et al. 2023).")

# ----------------------------------------------------------------------------- worst hour + demo truth
left, right = st.columns([1, 1])
with left:
    st.subheader("When is the network worst?")
    fig, ax = plt.subplots(figsize=(8, 3.2))
    worst_hour_bars(ax, f"% junctions below {threshold:g} mg/L (median)", "when you sample"); ax.legend(loc="upper center")
    st.pyplot(fig, width='stretch')
with right:
    if tr is not None:
        st.subheader("Demo only: reveal the hidden truth")
        if st.toggle("Show the true daily-minimum map and score the flags", value=False):
            tv = tr.truth_daily_min < float(threshold)
            uns = [j for j in junctions if j not in set(S.junction)]
            tp = int((tv.loc[uns] & (pmin.loc[uns, "p_below"] > 0.5)).sum()); fn = int((tv.loc[uns] & ~(pmin.loc[uns, "p_below"] > 0.5)).sum())
            fp = int((~tv.loc[uns] & (pmin.loc[uns, "p_below"] > 0.5)).sum())
            found = (f"On the unsampled junctions the map found **{tp} of {tp + fn}** (recall {tp / (tp + fn):.0%}) with {fp} false alarms. "
                     if tp + fn else
                     f"No unsampled junction is truly below the threshold, so recall is not testable here; the map raised {fp} false alarms. ")
            st.write(f"True daily-minimum violations: **{int(tv.sum())} of {len(junctions)}** junctions. {found}"
                     f"Mean of your samples: {S.y.mean():.2f} mg/L — which flags {'everything' if S.y.mean() < threshold else 'nothing'}.")
            fig, ax = plt.subplots(figsize=(8, 5))
            network_panel(ax, tr.truth_daily_min, "TRUE daily minimum (hidden from the model)", "viridis", (0, max(1.2, float(dose))))
            st.pyplot(fig, width='stretch')


# ----------------------------------------------------------------------------- what the model accounts for
st.subheader("What this model accounts for")
try:
    st.table(accounts_for(APP_ROOT).set_index("point"))
    st.caption(cloud.ACCOUNTS_FOR_NOTE if CLOUD else ACCOUNTS_FOR_NOTE)    # task 20: no folder names in the hosted demo
except (OSError, KeyError, IndexError, ValueError, TypeError) as e:     # outputs/ missing or changed on this machine
    st.info(f"The table of what this model accounts for is read from the committed results in outputs/, which could "
            f"not be read here ({type(e).__name__}). The README and docs/iteration3_journal.md carry the same answer.")


# ----------------------------------------------------------------------------- one-page PDF
def pdf_chem_lines() -> list[str]:
    """The chemistry lines of the PDF header (capabilities.report_lines): disinfectant and measured species, the
    threshold and its source, how the seasons are handled, and 'validated in simulation only'."""
    if ca:
        return report_lines(True, float(threshold), float(ph_log), float(cl2n_log), int(watch.sum()))
    return report_lines(False, float(threshold))


def pdf_bytes() -> bytes:
    fig = plt.figure(figsize=(16.5, 11.7))  # A3 landscape-ish, prints fine on A4
    gs = fig.add_gridspec(3, 4, height_ratios=[0.62, 1.5, 1.1])
    ax = fig.add_subplot(gs[0, :]); ax.axis("off")
    name = os.path.basename(net_path)
    pdf_hour_note = (" (a model hour: this file starts at 08:00)" if is_net2 else
                     f" (time of day; the file starts at {clock_label(0, clock_h)})" if use_clock else "")
    ax.text(0, 0.97, f"ResidualMap: monthly chlorine report, {name}", fontsize=18, fontweight="bold", va="top")
    gamma_pdf = "old-pipe factor not identifiable (one C value)" if uniform_c else f"old-pipe factor {model.map_params_[2]:.1f}" if model.map_params_ is not None else ""
    calib = (f"calibrated decay: bulk {model.map_params_[0]:{KB_FMT}} /day, wall {model.map_params_[1]:.2f} m/day, {gamma_pdf}; "
             f"demand ×{model.map_params_[3]:.2f}, roughness ×{model.map_params_[4]:.2f}, dose ×{model.map_dose_:.2f}") if model.map_params_ is not None else "no samples yet: map from the model's physics alone"
    # task 16c: with inlet residuals, the dose is the highest of them, and the residuals and the age caveat are printed
    # one line, so it never reaches the chemistry lines below it (at most 6 inlets named)
    inlet_pdf = "" if inlet_info is None else (
        "\ninlet residuals: " + ", ".join(f"{e} {mgl_text(v)}" for e, v in list(inlet_info["residuals_mgL"].items())[:6])
        + (f" and {len(inlet_info['residuals_mgL']) - 6} more" if len(inlet_info["residuals_mgL"]) > 6 else "")
        + " mg/L, each inlet dosed in proportion; water age starts at 0 at every inlet, so a zone's water is older "
          "than its age map shows")
    dose_words = "highest inlet residual" if inlet_info is not None else "dose"
    ax.text(0, 0.74, f"{len(S)} grab samples{' of TOTAL chlorine (chloramine mode)' if ca else ''} · {dose_words} {dose:g} mg/L · minimum residual {threshold:g} mg/L · "
                     f"{int((pmin['p_below'] > 0.5).sum())} of {len(junctions)} junctions likely below the minimum at their daily minimum · "
                     f"worst hour {hlabel(worst)}{pdf_hour_note}\n{calib}\n{oldest_line}{inlet_pdf}", fontsize=10.5, va="top")
    chem = "\n".join(w for line in pdf_chem_lines() for w in textwrap.wrap(line, PDF_CHEM_CHARS))
    ax.text(0, 0.36, chem, fontsize=9.5, va="top", color="0.15")
    axes = [fig.add_subplot(gs[1, i]) for i in range(4)]
    top = max(1.2, float(dose))
    network_panel(axes[0], med, f"Chlorine — {label} (mg/L)", "viridis", (0, top), list(S.junction), [htag(h) for h in S.hour],
                  below=below_detection, above=above_range)
    network_panel(axes[1], hi - lo, "90% band width (mg/L)", "magma", (0, float((hi - lo).quantile(0.98))))
    network_panel(axes[2], p_below, f"P(residual < {threshold:g} mg/L)", "Reds", (0, 1))
    network_panel(axes[3], p_below, f"Next route: {len(route)} sites", "Reds", (0, 1), list(route.junction), [route_tag(k, h) for k, h in enumerate(route.hour)])
    ax = fig.add_subplot(gs[2, :2]); ax.axis("off")
    ax.set_title("Sample here next", loc="left", fontsize=12, fontweight="bold")
    # the reason without its P(...) part (its own column), wrapped so none of it is cut off
    rows = [[str(j), hlabel(h), f"{p:.2f}", "\n".join(textwrap.wrap(w, PDF_WHY_CHARS))]
            for j, h, p, w in zip(route.junction, route.hour, route.p_below, route.why)]
    if not rows:      # task 15b: every junction is sampled already (a small network)
        rows = [["none", "", "", "every junction of this network is sampled already"]]
    tbl = ax.table(cellText=rows, colLabels=["junction", "hour", f"P(min<{threshold:g})", "why"], loc="upper left", cellLoc="left",
                   colWidths=[0.1, 0.1, 0.12, 0.68])
    n_lines = [1] + [r[3].count("\n") + 1 for r in rows]
    unit = 0.97 / (sum(n_lines) + 0.6 * len(n_lines))       # the table fills the panel's height
    for (r, c), cell in tbl.get_celld().items():
        cell.set_height(unit * (n_lines[r] + 0.6))
    # each text line gets about unit of the table's height: a long route gets a smaller font, never overlapping rows
    line_pt = unit * ax.get_position().height * fig.get_figheight() * 72
    tbl.auto_set_font_size(False); tbl.set_fontsize(min(7.5, 0.95 * line_pt))
    ax = fig.add_subplot(gs[2, 2:])
    worst_hour_bars(ax, f"% junctions below {threshold:g} mg/L", "sampling window")
    ax.set(title="When the network is worst"); ax.legend()
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="pdf"); plt.close(fig)
    return buf.getvalue()


st.download_button("Download the one-page PDF report", data=pdf_bytes(), file_name=f"residualmap_{os.path.basename(net_path)}.pdf", mime="application/pdf")
if CLOUD:      # task 20: the hosted demo says where its numbers are computed
    st.caption(cloud.footer(UPLOADS))
else:
    st.caption("Every number on this page is computed on this laptop from your .inp and your samples. Nothing is uploaded anywhere.")
st.page_link(GET_STARTED, label=cloud.LINK_NO_CHECK if not UPLOADS else
             "Bring your data: what to export from your software, templates for the grab log and tap map, and 'Check a "
             "file'")
