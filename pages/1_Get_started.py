"""
pages/1_Get_started.py: 'Bring your data', the app's second page (iteration 5, task 16).

    streamlit run app.py          # the page is in the sidebar's page list, under the main page

Step 1, 'Which software holds your network?': each package's export steps, read from docs/how_to_export_your_model.md
(one section per choice, so the guide and the page never disagree).  GIS and CAD are marked 'coming' until their
readers are built (tasks 17 and 18).

Step 2, templates: the grab log (free and total chlorine), the tap map by ID and the plant log, from docs/templates/,
and the grab log as an Excel workbook written on the fly by openpyxl (nothing binary is committed).

Step 3, 'Check a file': any file, identified by its content (importers.sniff).  An EPANET model gets the main page's
readiness report (importers.prepare_inp: WNTR reads it and EPANET runs it, as the main page does before its map; the
map's grid is never built), from a copy in a temporary folder under the git-ignored outputs/private/, deleted at once.
A table is read as a grab log on its own (its columns guessed, its rows counted, the main page's 3-month window
given, every row that would not be used listed with its reason), as a tap map, or as a junction list (a junction
column and no tap column: not a tap map, said so).  Anything else gets its advice.  Nothing is kept.
"""
from __future__ import annotations

import contextlib
import csv
import datetime as dt
import io
import os
import re
import sys
import tempfile

import streamlit as st

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if APP_ROOT not in sys.path:
    sys.path.insert(0, APP_ROOT)

from residualmap.importers import (LOG_WINDOW_MONTHS, LogRefused, dropped_lines, guess_log_columns,  # noqa: E402
                                   guess_tap_columns, have_openpyxl, prepare_inp, read_grab_log, read_table, sniff,
                                   xlsx_bytes)

GUIDE = os.path.join(APP_ROOT, "docs", "how_to_export_your_model.md")
TEMPLATES = os.path.join(APP_ROOT, "docs", "templates")
CHECK_DIR = os.path.join(APP_ROOT, "outputs", "private", "check_a_file")     # git-ignored; emptied after each check
SAMPLES_SECTION = "Your grab samples"
TEMPLATE_FILES = (("grab_log_template.csv", "Grab log, free chlorine (CSV)"),
                  ("grab_log_total_template.csv", "Grab log, total chlorine for a chloraminated system (CSV)"),
                  ("tap_map_by_id_template.csv", "Tap map by ID (CSV)"),
                  ("plant_log_template.csv", "Plant log, optional, read by the pilot script only (CSV)"))
XLSX_TEMPLATE = ("grab_log_template.xlsx", "Grab log, free chlorine (Excel)")
DISINFECTANTS = {"Free chlorine": "free_chlorine", "Chloramine (total chlorine)": "chloramine"}


def guide_sections(path: str = GUIDE) -> dict[str, str]:
    """The guide's '## ' sections: heading -> its text, in the guide's order."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    out = {}
    for part in re.split(r"(?m)^## ", text)[1:]:
        head, _, body = part.partition("\n")
        out[head.strip()] = body.strip()
    return out


def xlsx_template() -> bytes:
    """docs/templates/grab_log_template.csv as a workbook with typed date and time cells (openpyxl)."""
    with open(os.path.join(TEMPLATES, "grab_log_template.csv"), newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    out = [rows[0]]
    for r in rows[1:]:
        d, t, tap, v = r[0], r[1], r[2], r[3]
        try:
            v = float(v)
        except ValueError:
            pass                      # '<0.02' stays text, as a lab writes it
        out.append([dt.date.fromisoformat(d), dt.time.fromisoformat(t), tap, v] + r[4:])
    return xlsx_bytes(out, sheet="grab log")


def check_table(data: bytes, name: str, species: str, s) -> None:
    """'Check a file' for a table: as a tap map when it has a tap and a junction column and no reading column; as a
    junction list (a GIS or EPANET junction export) when it has a junction column but no tap column (the review of
    task 16: this crashed the page); else as a grab log read without a tap map, with the main page's window."""
    try:
        t = read_table(data, name)
    except LogRefused as e:
        st.warning(str(e) + ("" if s.category == "table" else f" {s.advice}"))
        return
    tc = guess_tap_columns(t.header)
    reading = guess_log_columns(t.header, species)["reading"]
    if tc["junction"] and not reading and tc["tap"]:
        n = len({str(r[t.header.index(tc["tap"])]).strip() for r in t.rows})
        st.markdown(f"**As a tap map:** {n} taps, tap column '{tc['tap']}', junction column '{tc['junction']}'. "
                    f"Upload it as the tap map on the main page, with the grab log.")
        return
    if tc["junction"] and not reading:
        st.markdown(f"**As a table:** it has a junction column ('{tc['junction']}') but no tap column and no chlorine "
                    f"reading, so it looks like the model's list of junctions (a GIS or EPANET export), not a tap map "
                    f"or a grab log. A tap map has two columns, the tap ID and the junction ID (see the tap map "
                    f"template in step 2).")
        return
    try:
        g = read_grab_log(data, name, None, disinfectant=species)
    except LogRefused as e:
        st.warning(f"As a grab log: {e}")
        return
    months = sorted(set(g.frame.month) | (set(g.censored.month) if len(g.censored) else set()))
    rows_w, cens_w, w = g.window(LOG_WINDOW_MONTHS)
    use = (f" The main page uses the last {LOG_WINDOW_MONTHS} months that hold a reading: {w['n_used']} readings "
           f"from {w['first']} to {w['last']}." if len(months) > LOG_WINDOW_MONTHS and w["first"] else "")
    st.markdown(f"**As a grab log:** {len(g.frame)} readings can be read, from {g.n_rows} rows"
                + (f", dated {min(g.frame.date)} to {max(g.frame.date)} ({len(months)} month"
                   f"{'s' if len(months) != 1 else ''})" if len(g.frame) else "") + "." + use + " " + " ".join(g.notes))
    not_used = g.dropped + [{"row": int(r.row), "reason": r.reason} for r in g.censored.itertuples()]
    if not_used:
        st.warning("Rows that would not be used, each with its reason:\n\n"
                   + "\n".join(f"- {ln}" for ln in dropped_lines(sorted(not_used, key=lambda x: x["row"]))))
    for x in g.warnings:
        st.warning(x)
    st.caption("Checked without a tap map, so every tap ID stands for itself here; on the main page the tap map "
               "matches each tap to a junction of the model.")


st.set_page_config(page_title="ResidualMap: bring your data", layout="wide")
st.title("Bring your data")
st.caption("Three steps: export the network from the software that holds it, put your grab samples in a log ResidualMap "
           "reads, and check any file before you use it on the main page. Everything runs on this laptop.")

# ----------------------------------------------------------------------------- step 1
st.header("Step 1. Which software holds your network?")
sections = guide_sections()
choices = [k for k in sections if k != SAMPLES_SECTION]
choice = st.radio("Which software holds your network?", choices, index=0)
st.markdown(sections[choice])
st.caption("'Coming soon' means the app already recognises this kind of file and says what to export, but cannot read "
           "it yet. ResidualMap reads EPANET models, and every package above can write one.")

# ----------------------------------------------------------------------------- step 2
st.header("Step 2. Your grab samples: templates")
samples_intro, _, samples_detail = sections[SAMPLES_SECTION].partition("\n\nYour own export")
st.markdown(samples_intro)
if samples_detail:
    with st.expander("What your own export can look like"):
        st.markdown("Your own export" + samples_detail)
for fname, label in TEMPLATE_FILES:
    with open(os.path.join(TEMPLATES, fname), "rb") as fh:
        st.download_button(label, data=fh.read(), file_name=fname, mime="text/csv", key=f"template_{fname}")
if have_openpyxl():
    st.download_button(XLSX_TEMPLATE[1], data=xlsx_template(), file_name=XLSX_TEMPLATE[0],
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key="template_xlsx")
else:
    st.info("The Excel template needs openpyxl, which is not installed here (pip install -r requirements.txt). The CSV "
            "templates above open in Excel too.")

# ----------------------------------------------------------------------------- step 3
st.header("Step 3. Check a file")
st.caption("Give it any file: a model, a drawing, a GIS layer, a grab log, a tap map. It says what the file is, by what "
           "is inside it, and what to send instead. A model gets the same readiness report as on the main page: the "
           "model is opened and run for one hour to make sure it works, but no map is built. Nothing is kept.")
species_label = st.radio("Disinfectant", list(DISINFECTANTS), index=0, horizontal=True, key="check_disinfectant")
up = st.file_uploader("Check a file", type=None, help="any file: it is identified by its content")
if up is not None:
    data, name, species = up.getvalue(), up.name, DISINFECTANTS[species_label]
    s = sniff(data, name)
    st.subheader("What this file is")
    st.markdown(f"**{s.sentence()}**")
    if s.is_model:
        os.makedirs(CHECK_DIR, exist_ok=True)
        try:
            with tempfile.TemporaryDirectory(prefix="check_", dir=CHECK_DIR) as d:
                path, rep = prepare_inp(data, d, species, name=name)
                st.markdown(rep.to_markdown())
        finally:
            with contextlib.suppress(OSError):
                os.rmdir(CHECK_DIR)              # only when empty: another check may be running
        st.caption("This is the readiness report the main page shows when you upload this model there.")
    elif s.category == "table" or s.kind == "text_other":     # text with a title line before its header too
        try:
            check_table(data, name, species, s)
        except Exception as e:  # noqa: BLE001  (the last guard, as on the main page: a file the reader did not foresee)
            st.warning(f"The file could not be checked as a table ({type(e).__name__}). Save it again as CSV.")
    else:
        st.warning(f"ResidualMap reads EPANET models. {s.advice}")
