"""
cloud.py: the hosted demo's mode (iteration 5, task 20), so ResidualMap can be tried without installing anything.

    RESIDUALMAP_MODE=cloud streamlit run app.py      # the hosted demo's behaviour, on this laptop
    streamlit run hosted/streamlit_app.py            # the hosted demo's own start (always cloud mode)

Mode.  RESIDUALMAP_MODE = "cloud" (a root-level Streamlit secret, which Streamlit also puts in the environment when
the server starts; a string, as a TOML true is not put there) gives cloud mode, and "local" gives local mode.  Unset,
the mode is cloud when the app runs from Streamlit Community Cloud's checkout under /mount/src/, so the privacy
default fails safe; any other value is cloud mode too.  With neither, it is local mode, exactly as before task 20.
The hosted start (hosted/streamlit_app.py) sets cloud mode itself, whatever the secrets say.

Cloud mode shows example networks only (Net3 and Net2; ky4's grid would take hours on a small shared container), no
file uploader on either page, and a privacy notice at the top of both pages.  Uploads come back only when
RESIDUALMAP_ALLOW_UPLOADS = "true" is set as well (open decision 8: off by default), with limits (UPLOAD_MB for every
uploaded file and for a zip's model once unpacked, a model of at most MAX_JUNCTIONS junctions, counted in its text
before it is parsed: intake_limits, limit_refusal): each upload is processed in memory and in a per-run temporary
folder that is deleted when the run ends (run_pages), its results are kept in this browser session only
(session_cached, in st.session_state, never in the process-wide caches every visitor shares), and nothing is written
under outputs/.  What the example networks compute at a visitor's settings still goes through the shared caches, as
in local mode.

The final verification of task 20 found uploads under both limits that cost far more than their size says, so the
limits also hold a model's links, tanks and reservoirs and lines of model data in its text, its pattern and rule time
steps once WNTR has read it, and its timed 7-day run (MAX_SECONDS_PER_RUN, MAX_RUN_MINUTES for the grid it gives);
a model's intake runs in a process of its own that is stopped after INTAKE_TIMEOUT_S (prepare_upload); and every
EPANET run on an upload runs in a process whose working directory is the run's temporary folder (run_in_folder,
child_main), because EPANET writes a scratch file of its own to the working directory and deletes it only when a
run ends normally.  An Excel grab log is held to WORKBOOK_MB once unpacked (table_limits).

Pages.  run_pages names the sidebar's pages 'ResidualMap' and 'Get started' with st.navigation (without it the main
page is listed as 'app', after app.py).  app.py is the entrypoint and its own main page: as the entrypoint it only
calls run_pages; run as the page, inside_page() is true and it skips that call.

Nothing here imports Streamlit at module level: simulate.py reads the mode, and so do the grid's worker processes.
"""
from __future__ import annotations

import contextlib
import os
import shutil
import tempfile
import threading

MODE_ENV = "RESIDUALMAP_MODE"
UPLOADS_ENV = "RESIDUALMAP_ALLOW_UPLOADS"
CLOUD_CHECKOUT = "/mount/src/"          # where Streamlit Community Cloud checks out the repository it runs
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRUE_WORDS = ("1", "true", "yes", "on")
UPLOAD_MB = 5                           # cloud uploads: the largest file taken, in MB (1e6 bytes)
MAX_JUNCTIONS = 300                     # cloud uploads: the largest model simulated, in junctions
# the final verification of task 20: limits on what a model under those two costs (intake_limits)
MAX_LINKS = 2 * MAX_JUNCTIONS           # cloud uploads: pipes, pumps and valves together
MAX_TANKS_RESERVOIRS = MAX_JUNCTIONS    # cloud uploads: tanks and reservoirs together
MAX_MODEL_LINES = 20_000                # cloud uploads: data lines of the sections EPANET simulates
MIN_PATTERN_STEP_S = 900                # cloud uploads: the shortest pattern time step (EPANET cuts its steps to it)
MIN_RULE_STEP_S = 60                    # cloud uploads with rule-based controls: the shortest rule time step
MAX_SECONDS_PER_RUN = 3.0               # cloud uploads: the longest 7-day run, timed on this server
MAX_RUN_MINUTES = 10                    # cloud uploads: the longest grid build, as the run cost estimates it here
INTAKE_TIMEOUT_S = 60                   # cloud uploads: a model's intake process is stopped after this long
WORKBOOK_MB = 4 * UPLOAD_MB             # cloud uploads: an Excel grab log or tap map once unpacked, in MB
GRID_WORKERS = 2                        # cloud mode: grid worker processes, so a small container keeps its memory
CLOUD_EXAMPLES = ("Net3", "Net2")       # cloud mode: the example networks offered
MAIN_TITLE, GET_STARTED_TITLE = "ResidualMap", "Get started"     # the sidebar's page names
SESSION_KEEP = 3                        # session_cached: results kept per kind in one session


# ----------------------------------------------------------------------------- mode
def mode(app_root: str | None = None, environ=None) -> str:
    """'cloud' or 'local' (module docstring): the variable first, then the checkout path, failing safe to cloud."""
    env = os.environ if environ is None else environ
    v = str(env.get(MODE_ENV, "")).strip().lower()
    if v == "local":
        return "local"
    if v:
        return "cloud"          # "cloud", or a value that is neither: the safe side
    root = os.path.abspath(REPO_ROOT if app_root is None else app_root)
    return "cloud" if (root + os.sep).startswith(CLOUD_CHECKOUT) else "local"


def is_cloud(app_root: str | None = None, environ=None) -> bool:
    return mode(app_root, environ) == "cloud"


def uploads_allowed(app_root: str | None = None, environ=None) -> bool:
    """Local mode always takes uploads; cloud mode only with RESIDUALMAP_ALLOW_UPLOADS set to a true word."""
    env = os.environ if environ is None else environ
    if not is_cloud(app_root, env):
        return True
    return str(env.get(UPLOADS_ENV, "")).strip().lower() in TRUE_WORDS


def grid_workers(app_root: str | None = None, environ=None) -> int | None:
    """Worker processes for a grid build: None (simgp's default, one fewer than the cores) in local mode, so local
    runs are unchanged; GRID_WORKERS in cloud mode, where the core count the container reports can be the host's."""
    return GRID_WORKERS if is_cloud(app_root, environ) else None


# ----------------------------------------------------------------------------- what the pages say
NOTICE_START = ("This is a public demo on Streamlit Community Cloud, run by Snowflake under its privacy policy.")
NOTICE_END = ("Do not use it for a utility's real data unless the utility agrees; the local version keeps everything "
              "on your computer.")


def notice(uploads: bool) -> str:
    """The privacy notice at the top of both pages in cloud mode."""
    if not uploads:
        return f"{NOTICE_START} It shows example networks only. {NOTICE_END}"
    return (f"{NOTICE_START} Uploads are processed in this server's memory and in temporary files deleted when "
            f"processing ends; nothing you upload is saved, logged or kept after your session. Your file does leave "
            f"your computer. It takes files up to {UPLOAD_MB} MB and models of up to {MAX_JUNCTIONS} junctions whose "
            f"grid it can build in about {MAX_RUN_MINUTES} minutes. {NOTICE_END}")


# task 17's review: the page has a fourth step since task 17, for a system with no hydraulic model
GIS_STEP_INTRO = ("A fourth step, for a system with no hydraulic model, builds a skeleton model from a GIS layer of its "
                  "water mains.")
GET_STARTED_LOCAL = ("Three steps: export the network from the software that holds it, put your grab samples in a log "
                     "ResidualMap reads, and check any file before you use it on the main page. " + GIS_STEP_INTRO
                     + " Everything runs on this laptop.")


def get_started_intro(cloud: bool, uploads: bool) -> str:
    """The Get started page's intro: 'Everything runs on this laptop' is true in local mode only."""
    if not cloud:
        return GET_STARTED_LOCAL
    if not uploads:
        return ("Three steps: export the network from the software that holds it, put your grab samples in a log "
                "ResidualMap reads, and check your files with the local version of ResidualMap before you use them. "
                "This public demo shows the steps and the templates; it runs on a server, not on your computer, so "
                "checking a file is turned off here.")
    return ("Three steps: export the network from the software that holds it, put your grab samples in a log "
            "ResidualMap reads, and check any file before you use it on the main page. " + GIS_STEP_INTRO + " This "
            "public demo runs on a server, not on your computer: a file you check here leaves your computer, and "
            "nothing is kept after your session.")


CHECK_OFF = ("Checking a file is turned off in this public demo, so nothing of yours leaves your computer. The local "
             "version of ResidualMap checks any file, a model, a drawing, a GIS layer, a grab log or a tap map, on "
             "your own computer and keeps nothing.")
SAMPLES_OFF = ("In this public demo, type your samples in the 'Grab samples' table. Uploading a grab log and a tap map "
               "needs the local version, which keeps your files on your computer.")
REDUCED_TRUTH = ("Reduced demo truth. This public demo does not run EPA's chloramine chemistry in EPANET-MSX: without a "
                 "C compiler it runs uncompiled and takes minutes per run, and whether this kind of server has a "
                 "compiler and the runtime EPANET-MSX needs is untested. The hidden truth is its reduced form instead: "
                 "EPANET first-order decay at the rate EPA's chemistry gives for this water, on the same hidden "
                 "network. It has the model's own kinetics, so the demo's scores flatter the model.")
LINK_NO_CHECK = "Bring your data: what to export from your software, and templates for the grab log and tap map"
DEMO_HELP_NO_UPLOADS = "Off = type your own grab samples in the table."
# capabilities.ACCOUNTS_FOR_NOTE without the repository's folder names, for the hosted demo's visitors
ACCOUNTS_FOR_NOTE = ("Every test in this table is a simulation, on EPANET's example networks Net3 and Net2 (and, for water "
                     "age, the KYPIPE network ky4), against hidden chemistry with assumed constants and schedules: no "
                     "real grab-sample data has been scored yet. Every result number is read from the project's "
                     "committed results when the app runs (settings and bars are fixed text); the project's journal "
                     "has the full results.")
SHARED_CACHES = ("The samples you type are kept in your session only; results for the example networks at the settings "
                 "you pick are kept on the server, shared by every visitor, until it restarts, and none of them says "
                 "who you are.")


def footer(uploads: bool) -> str:
    """The main page's last line in cloud mode (local mode keeps 'computed on this laptop')."""
    if not uploads:
        return ("Every number on this page is computed on this public demo's server, from an example network and the "
                f"samples in the table. {SHARED_CACHES}")
    return ("Every number on this page is computed on this public demo's server. An uploaded file is processed in "
            "memory and in temporary files deleted when processing ends, and what is made from it is kept in your "
            f"session only. {SHARED_CACHES}")


def too_big(n_bytes: int, unpacked: bool = False) -> str | None:
    """None when a cloud upload (or, unpacked, the model in an uploaded zip) is within UPLOAD_MB, else the plain
    refusal."""
    if n_bytes <= UPLOAD_MB * 1_000_000:
        return None
    what = "The model in this zip is {} MB once unpacked." if unpacked else "This file is {} MB."
    return (f"{what.format(f'{n_bytes / 1e6:.1f}')} The public demo takes files up to {UPLOAD_MB} MB; the local "
            f"version of ResidualMap takes larger files and keeps them on your computer.")


def too_many_junctions(n: int | None) -> str | None:
    """None when a cloud upload's model is within MAX_JUNCTIONS, else the plain refusal (before the grid is built)."""
    if n is None or int(n) <= MAX_JUNCTIONS:
        return None
    return (f"This model has {int(n):,} junctions. The public demo simulates models of up to {MAX_JUNCTIONS} "
            f"junctions; run a larger model, or a zone of up to {MAX_JUNCTIONS} junctions, with the local version, "
            f"which keeps it on your computer.")


def intake_limits(app_root: str | None = None) -> dict:
    """importers.prepare_inp's limits for a cloud upload: a model of at most UPLOAD_MB (a zip's by its declared size,
    before it is unpacked) and MAX_JUNCTIONS junctions (counted in the text before WNTR parses it, so a model over the
    limit costs no load and no EPANET run), and the run cost estimated for the cloud's grid workers.  Since the final
    verification of task 20 also, in the text, at most MAX_LINKS links, MAX_TANKS_RESERVOIRS tanks and reservoirs and
    MAX_MODEL_LINES lines of model data; before EPANET runs it, a pattern time step of at least MIN_PATTERN_STEP_S and,
    with rules, a rule time step of at least MIN_RULE_STEP_S; after its timed run, at most MAX_SECONDS_PER_RUN a run
    and MAX_RUN_MINUTES for its grid."""
    return {"max_junctions": MAX_JUNCTIONS, "max_model_bytes": UPLOAD_MB * 1_000_000,
            "workers": GRID_WORKERS if is_cloud(app_root) else None, "max_links": MAX_LINKS,
            "max_tanks_reservoirs": MAX_TANKS_RESERVOIRS, "max_model_lines": MAX_MODEL_LINES,
            "min_pattern_step_s": MIN_PATTERN_STEP_S, "min_rule_step_s": MIN_RULE_STEP_S,
            "max_seconds_per_run": MAX_SECONDS_PER_RUN, "max_run_minutes": MAX_RUN_MINUTES}


def table_limits() -> dict:
    """importers.read_table's (and read_grab_log's) limits for a cloud upload (the final verification of task 20): an
    Excel workbook of at most WORKBOOK_MB unpacked, and none of its parts packed more than 100 to 1, before openpyxl
    opens it (a 466 KB workbook had unpacked to 190 MB)."""
    return {"max_unpacked": WORKBOOK_MB * 1_000_000}


LOCAL_TAIL = "run it with the local version of ResidualMap, which keeps it on your computer."


def _minutes(seconds: float) -> str:
    """A time step in words: '1 second', '30 seconds', '1 minute', '15 minutes'."""
    s = float(seconds)
    n, unit = (s, "second") if s < 60 else (s / 60, "minute")
    return f"{n:g} {unit}{'' if n == 1 else 's'}"


# what a model over each of intake_limits' limits is told (the refusal names the model's own value)
REFUSALS = {
    "links": lambda v: (f"This model has {v:,} links (pipes, pumps and valves). The public demo simulates models of up "
                        f"to {MAX_LINKS} links;"),
    "tanks_reservoirs": lambda v: (f"This model has {v:,} tanks and reservoirs. The public demo simulates models of up "
                                   f"to {MAX_TANKS_RESERVOIRS} of them;"),
    "model_lines": lambda v: (f"This model has {v:,} lines of model data, its title and drawing aside. The public demo "
                              f"takes models of up to {MAX_MODEL_LINES:,} such lines;"),
    "pattern_step_s": lambda v: (f"This model's demand patterns change every {_minutes(v)}, and EPANET solves the "
                                 f"network at each change. The public demo takes patterns that change every "
                                 f"{_minutes(MIN_PATTERN_STEP_S)} or less often;"),
    "rule_step_s": lambda v: (f"This model checks its rules every {_minutes(v)}. The public demo takes rules checked "
                              f"every {_minutes(MIN_RULE_STEP_S)} or less often;"),
    "seconds_per_run": lambda v: (f"One 7-day simulation of this model took {v:.1f} s on this server. The public demo "
                                  f"takes models that run in {MAX_SECONDS_PER_RUN:g} s or less;"),
    "run_minutes": lambda v: (f"Building this model's grid would take about {v:.0f} minutes on this server. The public "
                              f"demo takes models whose grid builds in {MAX_RUN_MINUTES} minutes or less;"),
    "intake_seconds": lambda v: (f"Checking this model took longer than {v:g} s on this server, so it was stopped and "
                                 f"nothing was kept. The public demo takes models it can check in that time;"),
    "intake_failed": lambda v: ("Checking this model stopped before it finished, and nothing was kept. The public demo "
                                "could not check it;"),
}


def limit_refusal(rep) -> str | None:
    """The plain refusal for a cloud upload's model over a limit (prepare_inp's rep.flags['over_limit'], or, as a
    backstop, its run cost's junctions), else None."""
    over = rep.flags.get("over_limit") or {}
    what = over.get("what")
    if what == "model_bytes":
        return too_big(over["value"], unpacked=bool(over.get("unpacked")))
    if what == "junctions":
        return too_many_junctions(over["value"])
    if what in REFUSALS:
        return f"{REFUSALS[what](over['value'])} {LOCAL_TAIL}"
    return too_many_junctions(rep.run_cost.get("junctions"))


# task 17's review: a GIS layer is read and joined in the Streamlit process, so in the hosted demo a layer over these
# is refused right after it is read, before any line is joined (a 2 MB zipped shapefile of 167,620 lines took 37 s and
# 1.7 GB before it was refused for its junctions)
MAX_GIS_LINES = MAX_LINKS               # cloud uploads: a GIS layer's lines
MAX_GIS_POINTS = 50 * MAX_LINKS         # cloud uploads: a GIS layer's points (its lines' vertices)


def gis_too_big(n_lines: int, n_points: int) -> str | None:
    """None when a GIS layer is within the hosted demo's limits, else the plain refusal."""
    if n_lines <= MAX_GIS_LINES and n_points <= MAX_GIS_POINTS:
        return None
    return (f"This layer has {n_lines:,} lines and {n_points:,} points. The public demo builds layers of up to "
            f"{MAX_GIS_LINES} lines and {MAX_GIS_POINTS:,} points, since it runs models of up to {MAX_JUNCTIONS} "
            f"junctions; build a larger layer with the local version, on your own computer.")



# task 18: a CAD drawing in the hosted demo is read in a process of its own (run_in_folder) with a time limit and
# tighter limits than the local version's, then held to the GIS step's line and point limits, and to a limit on its
# texts, before any line is joined; since the review its build (labels, joining, the tolerance table) runs in a process
# of its own with a time limit too, and measures at most CAD_MAX_LABEL_PAIRS pairs of a label and a line segment
CAD_MAX_ENTITIES = 100_000              # cloud uploads: entities visited in one drawing, exploded blocks included
CAD_MAX_POINTS = 10 * MAX_GIS_POINTS    # cloud uploads: line points read from one drawing
CAD_MAX_INSERT_COPIES = 1_000           # cloud uploads: copies of one block array (MINSERT)
CAD_MAX_TEXTS = 10 * MAX_GIS_LINES      # cloud uploads: texts on the layer styles chosen (ky4's drawing has 557)
CAD_MAX_LABEL_PAIRS = 2_000_000         # cloud uploads: (label, line segment) pairs within reach, counted before any
CAD_READ_TIMEOUT_S = 30                 # cloud uploads: a drawing's read is stopped after this long
CAD_BUILD_TIMEOUT_S = 30                # cloud uploads: and its build


def cad_limits() -> dict:
    """cad_read.read_drawing's and build_from_drawing's limits for a cloud upload."""
    return {"max_entities": CAD_MAX_ENTITIES, "max_points": CAD_MAX_POINTS, "max_insert_copies": CAD_MAX_INSERT_COPIES,
            "max_label_pairs": CAD_MAX_LABEL_PAIRS}


def cad_too_big(n_lines: int, n_points: int, n_texts: int = 0) -> str | None:
    """None when a drawing's lines and texts are within the hosted demo's limits (the GIS step's for lines and
    points, CAD_MAX_TEXTS for texts), else the plain refusal."""
    if n_lines <= MAX_GIS_LINES and n_points <= MAX_GIS_POINTS and n_texts <= CAD_MAX_TEXTS:
        return None
    return (f"This drawing has {n_lines:,} lines, {n_points:,} points and {n_texts:,} texts on the layer styles chosen. "
            f"The public demo builds drawings of up to {MAX_GIS_LINES} lines, {MAX_GIS_POINTS:,} points and "
            f"{CAD_MAX_TEXTS:,} texts, since it runs models of up to {MAX_JUNCTIONS} junctions; leave out the layer "
            f"styles that are not water mains or their labels, or build the drawing with the local version, on your "
            f"own computer.")


def read_cad_upload(data: bytes, folder: str, name: str, timeout: float | None = None):
    """cad_read.read_drawing on a cloud upload, with cad_limits, in a process of its own whose working directory is
    folder (run_in_folder), stopped after CAD_READ_TIMEOUT_S (or timeout).  Returns a cad_read.Drawing; one whose read
    ran past the time limit, or whose process failed, comes back with a refusal that says so."""
    from .cad_read import Drawing, read_drawing
    timeout = CAD_READ_TIMEOUT_S if timeout is None else float(timeout)
    try:
        return run_in_folder(folder, read_drawing, data, name, cad_limits(), timeout=timeout)
    except ChildFailed as e:
        d = Drawing(name=os.path.basename(name))
        d.refusal = ((f"Reading this drawing took longer than {timeout:g} s, so it was stopped; nothing was kept."
                      if e.timed_out else "Reading this drawing stopped before it finished; nothing was kept.")
                     + " Save a drawing of the water mains alone as DXF, or read it with the local version of "
                       "ResidualMap, on your own computer.")
        return d


def build_cad_upload(drawing, opts, folder: str, timeout: float | None = None):
    """cad_read.build_from_drawing on a cloud upload's drawing, with cad_limits, in a process of its own whose working
    directory is folder (run_in_folder), stopped after CAD_BUILD_TIMEOUT_S (or timeout), so no drawing within the
    upload limits holds the server's process (the review measured gigabytes there for 600 lines and thousands of labels).
    Returns a cad_read.CadBuild; one whose build ran past the time limit, or whose process failed, comes back with a
    refusal that says so."""
    from .cad_read import CadBuild, build_from_drawing
    timeout = CAD_BUILD_TIMEOUT_S if timeout is None else float(timeout)
    try:
        return run_in_folder(folder, build_from_drawing, drawing, opts, cad_limits(), timeout=timeout)
    except ChildFailed as e:
        return CadBuild(name=drawing.name, drawing=drawing,
                        refusal=(f"Building this drawing took longer than {timeout:g} s, so it was stopped; nothing was "
                                 f"kept." if e.timed_out else "Building this drawing stopped before it finished; nothing "
                                 "was kept.") + " Leave out the layer styles that are not water mains or their labels, "
                                                "or build it with the local version of ResidualMap, on your own computer.")

# task 17: the GIS sections' last paragraph names the Get started page's step 4, which the hosted demo shows only with
# uploads allowed
GIS_STEP_LOCAL = "Step 4 of the Get started page builds a skeleton model from the layer:"
GIS_STEP_CLOUD = "The local version of ResidualMap builds a skeleton model from the layer, on your own computer:"
# task 18: the CAD sections name step 5 the same way
CAD_STEP_LOCAL = "Step 5 of the Get started page builds a skeleton model from the DXF:"
CAD_STEP_CLOUD = "The local version of ResidualMap builds a skeleton model from the DXF, on your own computer:"
CAD_MSTN_LOCAL = "then use step 5 of the Get started page (previous section)"
CAD_MSTN_CLOUD = "then the local version of ResidualMap reads it (previous section)"
CAD_GUIDE_SWAPS = (("AutoCAD or Civil 3D", CAD_STEP_LOCAL, CAD_STEP_CLOUD, True),
                   ("MicroStation", CAD_MSTN_LOCAL, CAD_MSTN_CLOUD, True))
# The Get started page's guide (docs/how_to_export_your_model.md) in cloud mode: (section, the local words, the cloud
# words, only when uploads are off).  The guide stays as it is for the local version; guide_section swaps the words
# that tell a visitor to upload, or name the repository's folders, and the saved check asserts each one applies.
GUIDE_CLOUD = (
    ("EPANET (.net or .inp)", "1. An `.inp` is ready: upload it as it is.",
     "1. An `.inp` is ready as it is: the local version of ResidualMap reads it.", True),
    ("Not sure", 'Go to step 3 of the Get started page, "Check a file", and give it any file you have: it says what the '
                 'file is, by its content, and what to send instead. If the file is a model, it shows the readiness '
                 'report.',
     "The local version of ResidualMap checks any file you have, on your own computer: it says what the file is, by "
     "its content, and what to send instead, and for a model it shows the readiness report. This public demo does not "
     "take files.", True),
    ("Your grab samples", 'On the main page, turn the demo off and upload two files under "2. Your samples": the grab '
                          'log, with your chlorine readings, and the tap map, which says which junction of the model '
                          'each sampling tap is at.',
     'The local version of ResidualMap takes two files, on its main page under "2. Your samples" with the demo off: the '
     'grab log, with your chlorine readings, and the tap map, which says which junction of the model each sampling tap '
     'is at. In this public demo, type a few samples into the main page\'s "Grab samples" table instead.', True),
    ("Your grab samples", " (in the repository they are in `docs/templates/`)", "", False),
) + tuple((head, GIS_STEP_LOCAL, GIS_STEP_CLOUD, True) for head in ("ArcGIS Pro or ArcGIS Online", "QGIS",
                                                                    "Google Earth")) + CAD_GUIDE_SWAPS


def guide_section(head: str, text: str, cloud: bool, uploads: bool) -> str:
    """A guide section as the Get started page shows it: unchanged in local mode; in cloud mode with GUIDE_CLOUD's
    words."""
    if not cloud:
        return text
    for h, old, new, no_uploads_only in GUIDE_CLOUD:
        if h == head and (uploads is False or not no_uploads_only):
            text = text.replace(old, new)
    return text


# ----------------------------------------------------------------------------- a visitor's EPANET runs, in the folder
class ChildFailed(RuntimeError):
    """run_in_folder's process ended without a result: stopped at its time limit (timed_out) or failed (code)."""

    def __init__(self, message: str, timed_out: bool = False, code: int | None = None, detail: str = ""):
        super().__init__(message)
        self.timed_out, self.code, self.detail = timed_out, code, detail


# the process run_in_folder starts: the repository first on its path (the folder holds a visitor's files, never
# modules), then child_main
CHILD_CODE = "import sys; sys.path[0] = {root!r}; from residualmap.cloud import child_main; child_main(sys.argv[1:])"


def run_in_folder(folder: str, fn, *args, timeout: float | None = None, **kwargs):
    """fn(*args, **kwargs) in a Python process of its own whose working directory is folder (the final verification of
    task 20).  EPANET writes a scratch file of its own (named 'en' and six letters, as large as the run's hydraulics)
    to the working directory, which for the hosted demo's server is the repository's checkout, and deletes it only
    when a run ends normally; with folder as the working directory it goes with the folder, and so do the files of
    the grid's worker processes, which that process starts and which inherit it.  The process's temporary directory
    is folder too (TMPDIR), so the temporary folders it makes are in folder even when it is stopped.  fn must be
    importable (a function of a module); its arguments and its result pass through two pickle files in folder,
    removed afterwards.  timeout: seconds before the process, and every process it started, is stopped (ChildFailed
    with timed_out).  Raises ChildFailed when the process ends without a result, and an error fn raised as ChildFailed
    with its type and message."""
    import pickle
    import signal
    import subprocess
    import sys
    import uuid
    tag = uuid.uuid4().hex[:12]
    job, out = os.path.join(folder, f".rm_job_{tag}.pkl"), os.path.join(folder, f".rm_out_{tag}.pkl")
    try:
        with open(job, "wb") as fh:
            pickle.dump((fn, args, kwargs), fh, protocol=pickle.HIGHEST_PROTOCOL)
        env = {**os.environ, "TMPDIR": folder, "TEMP": folder, "TMP": folder}     # its temporary files there too
        proc = subprocess.Popen([sys.executable, "-c", CHILD_CODE.format(root=REPO_ROOT), job, out], cwd=folder,
                                env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                start_new_session=True)
        try:
            _, err = proc.communicate(timeout=timeout)
        except BaseException as e:      # the time limit, or anything that stops this thread while it waits
            with contextlib.suppress(OSError):
                os.killpg(proc.pid, signal.SIGKILL)       # the process and the workers it started (its own group)
            proc.communicate()
            if isinstance(e, subprocess.TimeoutExpired):
                raise ChildFailed(f"it ran past its time limit of {timeout:g} s", timed_out=True) from None
            raise
        if not os.path.exists(out):
            tail = (err or b"").decode("utf-8", "replace").strip().splitlines()[-1:]
            raise ChildFailed(f"its process ended with code {proc.returncode} and no result", code=proc.returncode,
                              detail=tail[0][:300] if tail else "")
        with open(out, "rb") as fh:
            ok, value = pickle.load(fh)
    finally:
        for f in (job, out, out + ".part"):
            with contextlib.suppress(OSError):
                os.remove(f)
    if not ok:
        raise ChildFailed(f"{value[0]}: {value[1]}")
    return value


def child_main(argv: list[str]) -> None:
    """The process run_in_folder starts: argv is (job file, result file).  The result is written under a temporary
    name and moved into place, so a result file is always whole."""
    import pickle
    job, out = argv[0], argv[1]
    with open(job, "rb") as fh:
        fn, args, kwargs = pickle.load(fh)
    try:
        result = (True, fn(*args, **kwargs))
    except Exception as e:  # noqa: BLE001  (the caller gets the error's type and message)
        result = (False, (type(e).__name__, str(e)[:500]))
    with open(out + ".part", "wb") as fh:
        pickle.dump(result, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(out + ".part", out)


def prepare_upload(data: bytes, folder: str, disinfectant: str, name: str, app_root: str | None = None,
                   timeout: float | None = None, limits: dict | None = None):
    """importers.prepare_inp on a cloud upload, with intake_limits (or limits, for a check), in a process of its own
    whose working directory is folder (run_in_folder), stopped after INTAKE_TIMEOUT_S (or timeout).  Returns what
    prepare_inp returns; a model whose intake ran past the time limit, or whose process failed, comes back blocked,
    with rep.flags['over_limit'] ('intake_seconds' or 'intake_failed') and whatever its process left in folder
    removed."""
    import hashlib
    from .importers import Readiness, prepare_inp
    timeout = INTAKE_TIMEOUT_S if timeout is None else float(timeout)
    lim = {**intake_limits(app_root), **(limits or {})}
    before = set(os.listdir(folder))
    try:
        return run_in_folder(folder, prepare_inp, data, folder, disinfectant, name=name, timeout=timeout, **lim)
    except ChildFailed as e:
        for f in set(os.listdir(folder)) - before:           # the copy it was preparing, EPANET's scratch file
            p = os.path.join(folder, f)
            if os.path.isdir(p) and not os.path.islink(p):
                shutil.rmtree(p, ignore_errors=True)
            else:
                with contextlib.suppress(OSError):
                    os.remove(p)
        rep = Readiness(name=os.path.basename(name), disinfectant=disinfectant, sha1=hashlib.sha1(data).hexdigest())
        if e.timed_out:
            rep.flags["over_limit"] = {"what": "intake_seconds", "value": timeout, "limit": timeout}
            rep.blocking.append(f"Checking this model took longer than {timeout:g} s, so it was stopped; nothing was "
                                f"kept.")
        else:
            rep.flags["over_limit"] = {"what": "intake_failed", "value": e.code, "limit": None}
            rep.blocking.append("Checking this model stopped before it finished; nothing was kept. Send the model's "
                                "EPANET export (.inp) again, on its own.")
        return None, rep


def upload_grid(path: str, dose: float, chloramine: bool, grid: str, workers: int | None,
                inlet_ratios: dict | None = None):
    """app.py's prepare for a visitor's upload in cloud mode, run in the folder (run_in_folder): the nominal model, its
    features and the grid, built in memory (simgp.build_grid_24h) and never cached on disk.  (sc, X, (params, Z)).
    inlet_ratios (task 16c): a visitor's inlet residuals as dose ratios, so their grid stays in their session too (on
    an example network as well); None is the grid before task 16c."""
    from . import chloramine as CA
    from .features import build_features
    from .simgp import build_grid_24h
    from .simulate import nominal_scenario
    sc = nominal_scenario(path, sample_hour=14, source_dose=dose)
    X = build_features(sc)
    kw = {"inlet_ratios": inlet_ratios} if inlet_ratios else {}
    if chloramine:
        data = build_grid_24h(sc, "chloramine", workers, CA.ca_condition(), **kw)
    else:
        data = build_grid_24h(sc, grid, workers, **kw)
    return sc, X, data


def upload_truth(path: str, dose: float, seed: int, kb: float, kw: float):
    """app.py's demo_truth for a visitor's upload in cloud mode, run in the folder (run_in_folder)."""
    from .simulate import build_scenario
    return build_scenario(path, seed=seed, sample_hour=14, source_dose=dose, kb_per_day=kb, kw_m_per_day=kw)


def upload_truth_reduced(path: str, dose: float, seed: int):
    """app.py's demo_truth_ca_reduced for a visitor's upload in cloud mode, run in the folder (run_in_folder)."""
    from . import chloramine as CA
    kw = CA.KW_REF.get(os.path.basename(path), 0.20)
    return CA.truth_scenario(path, seed, kw, kinetics="first", dose=dose)


# ----------------------------------------------------------------------------- pages and the per-run folder
_PAGE = threading.local()       # Streamlit runs each browser session's script in its own thread


def inside_page() -> bool:
    """True while run_pages runs a page in this thread (so app.py, run as its own main page, skips run_pages)."""
    return getattr(_PAGE, "depth", 0) > 0


def run_pages(main_file: str, get_started_file: str) -> None:
    """The entrypoint's whole job: the two pages under their names in the sidebar, and the chosen page run.  The
    per-run temporary folder (scratch_dir) is deleted when the page's run ends, however it ends (st.stop, a rerun,
    an error)."""
    import streamlit as st
    page = st.navigation([st.Page(os.path.abspath(main_file), title=MAIN_TITLE, default=True),
                          st.Page(os.path.abspath(get_started_file), title=GET_STARTED_TITLE)])
    _PAGE.depth = getattr(_PAGE, "depth", 0) + 1
    try:
        page.run()
    finally:
        _PAGE.depth -= 1
        drop_scratch()


def scratch_dir(token: str) -> str:
    """This run's temporary folder for an upload (cloud mode), under the system's temporary directory: the same path
    for every run of one browser session (token), so a model prepared in an earlier run is written back where its
    cached results expect it; created on first use in a run and deleted when the run ends (run_pages)."""
    path = getattr(_PAGE, "scratch", None)
    if path is None:
        path = os.path.join(tempfile.gettempdir(), f"rm_cloud_{token}")
        os.makedirs(path, mode=0o700, exist_ok=True)
        _PAGE.scratch = path
    return path


def drop_scratch() -> None:
    path = getattr(_PAGE, "scratch", None)
    _PAGE.scratch = None
    if path:
        shutil.rmtree(path, ignore_errors=True)


@contextlib.contextmanager
def temporary_folder(prefix: str, local_dir: str | None):
    """A temporary folder that is deleted afterwards: under local_dir in local mode (as before task 20), under the
    system's temporary directory in cloud mode, so nothing is written under outputs/."""
    if local_dir is not None and not is_cloud():
        os.makedirs(local_dir, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory(prefix=prefix, dir=local_dir) as d:
                yield d
        finally:
            with contextlib.suppress(OSError):
                os.rmdir(local_dir)              # only when empty: another check may be running
        return
    with tempfile.TemporaryDirectory(prefix=prefix) as d:
        yield d


def session_cached(state, kind: str, key: tuple, fn, *args, **kwargs):
    """fn(*args, **kwargs), kept in this browser session only (state is st.session_state): the cloud route's cache
    for anything made from a visitor's upload, so it never sits in st.cache_data or st.cache_resource, which every
    visitor of the server shares.  At most SESSION_KEEP results per kind are kept; it all goes with the session."""
    store = state.setdefault("_residualmap_session_results", {})
    k = (kind, key)
    if k in store:
        return store[k]
    value = fn(*args, **kwargs)
    mine = [x for x in store if x[0] == kind]
    for old in mine[:max(0, len(mine) - SESSION_KEEP + 1)]:
        del store[old]
    store[k] = value
    return value
