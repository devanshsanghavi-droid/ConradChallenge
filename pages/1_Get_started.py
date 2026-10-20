"""
pages/1_Get_started.py: 'Bring your data', the app's second page (iteration 5, task 16).

    streamlit run app.py          # the page is in the sidebar's page list, under the main page

Step 1, 'Which software holds your network?': each package's export steps, read from docs/how_to_export_your_model.md
(one section per choice, so the guide and the page never disagree).  CAD is marked 'coming' until its reader is built
(task 18); since task 17 the GIS sections point to step 4.

Step 2, templates: the grab log (free and total chlorine), the tap map by ID and the plant log, from docs/templates/,
and the grab log as an Excel workbook written on the fly by openpyxl (nothing binary is committed).

Step 3, 'Check a file': any file, identified by its content (importers.sniff).  An EPANET model gets the main page's
readiness report (importers.prepare_inp: WNTR reads it and EPANET runs it, as the main page does before its map; the
map's grid is never built), from a copy in a temporary folder under the git-ignored outputs/private/, deleted at once.
A table is read as a grab log on its own (its columns guessed, its rows counted, the main page's 3-month window
given, every row that would not be used listed with its reason), as a tap map, or as a junction list (a junction
column and no tap column: not a tap map, said so).  Anything else gets its advice.  Nothing is kept.

Iteration 5, task 20 (residualmap/cloud.py): in cloud mode a privacy notice opens the page, the intro no longer says
everything runs on this laptop, and 'Check a file' is turned off (no uploader) unless cloud uploads are allowed; then a
file of up to cloud.UPLOAD_MB MB is checked in a temporary folder under the system's temporary directory, not under
outputs/, deleted at once, with the hosted demo's model limits (cloud.intake_limits): since the final verification
of task 20 the model is checked in a process of its own, in that folder and with a time limit (cloud.prepare_upload),
and an Excel table is held to cloud.WORKBOOK_MB once unpacked (TABLE_LIMITS).  The guide's words that tell a
visitor to upload, or name the repository's folders, are swapped for the hosted demo's (cloud.guide_section).  The
sidebar names this page 'Get started' (cloud.run_pages).

Iteration 5, task 17 (residualmap/gis_read.py, residualmap/topology.py): step 4, 'Build a model from a GIS layer', takes
a GeoPackage, a zipped shapefile, GeoJSON or a KML or KMZ (or WNTR's own GIS files with its JSON), asks what the layer
does not say (the coordinates' unit, the unit of bare pipe sizes, whether a layer of a chosen unit is to scale), offers
the guessed fields and the joining distance, and shows the network's counts at half, one and twice that distance, the
warnings, the IDs changed for EPANET, a valve layer's proposed closed valves and the skeleton's readiness report (not
ready: no source, tank, pump, valve, control or demand).  'Download the model' gives the .inp; 'Open it on the main
page' keeps it in this browser session (st.session_state), writes it to outputs/app_uploads/ in local mode (nothing is
written under outputs/ in cloud mode) and opens the main page, whose Model radio then offers 'Model prepared on Get
started'.  Step 4 comes after step 3, so every widget before it keeps its place, and in cloud mode with uploads off the
page stops before both.  Since the review of task 17: nothing in step 4 can raise (gis_read.build_from_layer never
does), the unit choices include degrees (longitude and latitude), a layer is picked by its position in the file (two
may share a name), each layer is built once, or twice when an answer changes it, and in cloud mode a GIS zip is read
within cloud.UPLOAD_MB once unpacked and a layer over cloud.gis_too_big's limits is refused before any line is joined.
Since its final verification: every reader takes time linear in the file (the KML text and coordinate system
definitions were read with patterns that took hours on a KMZ of 5 KB, before any limit of the hosted demo applied);
the coordinates' unit is offered for a layer whose file names a projected system but whose coordinates lie in
longitude and latitude ranges, and after a refusal as crowded, with the file's own unit first (there was no control to
change it); and a layer whose latitude and longitude look swapped asks their order.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import io
import os
import re
import sys

import streamlit as st
from streamlit.errors import StreamlitAPIException

APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if APP_ROOT not in sys.path:
    sys.path.insert(0, APP_ROOT)

from residualmap import cloud  # noqa: E402
from residualmap import gis_read as G  # noqa: E402
from residualmap import topology as TOPO  # noqa: E402
from residualmap.importers import (LOG_WINDOW_MONTHS, LogRefused, dropped_lines, guess_log_columns,  # noqa: E402
                                   guess_tap_columns, have_openpyxl, prepare_inp, prepared_name, read_grab_log,
                                   read_table, sniff, xlsx_bytes)

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
CLOUD = cloud.is_cloud(APP_ROOT)              # task 20: the hosted demo
UPLOADS = cloud.uploads_allowed(APP_ROOT)     # always True in local mode
TABLE_LIMITS = cloud.table_limits() if CLOUD else {}     # task 20: a workbook's size once unpacked, in the hosted demo
MAIN_PAGE = os.path.join(APP_ROOT, "app.py")              # task 17: the hand-off of a model built here
UPLOAD_DIR = os.path.join(APP_ROOT, "outputs", "app_uploads")     # git-ignored; where the main page keeps its models
GIS_DIR = os.path.join(APP_ROOT, "outputs", "private", "gis_check")   # git-ignored; emptied after each readiness check
GIS_HEADER = "Step 4. Build a model from a GIS layer"
GIS_INTRO = ("No hydraulic model? Give the GIS layer of your water mains (a GeoPackage, a zipped shapefile, GeoJSON, or "
             "a Google Earth KML or KMZ). The pipes are joined where their ends meet, their sizes are read, and a "
             "skeleton model is written: pipes and nodes only, with no source, tank, pump, valve, control or demand, "
             "so it is not ready to run until those are added. If a hydraulic model exists, its EPANET export is "
             "always better.")
GIS_HELP = "any file: it is identified by its content"
UNIT_CHOICES = {"(choose the unit)": None, "feet": "ft", "US survey feet": "us_ft", "metres": "m",
                "degrees (longitude and latitude)": "deg"}
GIS_UNPACKED = cloud.UPLOAD_MB * 1_000_000 if CLOUD else None   # the hosted demo: a GIS zip's members once unpacked
DIAMETER_UNITS = {"(choose the unit)": None, "inches": "in", "millimetres": "mm", "metres": "m"}
AXES_CHOICES = {"(choose the order)": None, "latitude first: swap them": True, "longitude first: keep them": False}
NO_FIELD = "(none)"
FIELD_LABELS = (("diameter", "Pipe size field"), ("material", "Material field"), ("year", "Install year field"),
                ("length", "Declared length field"), ("id", "Pipe ID field"))


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
        t = read_table(data, name, **TABLE_LIMITS)
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
        g = read_grab_log(data, name, None, disinfectant=species, **TABLE_LIMITS)
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


def _known_lengths(text: str) -> dict:
    """'pipe ID, length in feet' lines -> {pipe ID: feet}; lines that do not read are left out."""
    out = {}
    for ln in (text or "").splitlines():
        parts = [x.strip() for x in re.split(r"[,;\t]", ln) if x.strip()]
        if len(parts) == 2:
            try:
                out[parts[0]] = float(parts[1].replace(",", "."))
            except ValueError:
                continue
    return out


def gis_readiness(b) -> object:
    """The skeleton's readiness report: importers.prepare_inp on it, in a temporary folder deleted at once (under the
    git-ignored outputs/private/ in local mode; the system's temporary directory and the hosted demo's limits in cloud
    mode), with the GIS layer's rows first (gis_read.readiness)."""
    data, name = b.inp_text.encode("utf-8"), f"{TOPO.clean_id(b.layer.name if b.layer else 'skeleton')[0]}.inp"
    with cloud.temporary_folder("gis_", GIS_DIR) as d:
        _, rep = (cloud.prepare_upload(data, d, "free_chlorine", name, APP_ROOT) if CLOUD else
                  prepare_inp(data, d, "free_chlorine", name=name))
    if CLOUD and cloud.limit_refusal(rep):
        return None, cloud.limit_refusal(rep)
    return (rep if b.wntr_schema else G.readiness(b, rep)), None


def hand_off(b, inp_name: str) -> None:
    """Keep the prepared model for the main page (st.session_state, this browser session only) and, in local mode,
    write it to outputs/app_uploads/ under the main page's own name for it; then open the main page on it."""
    data = b.inp_text.encode("utf-8")
    if not CLOUD:
        os.makedirs(UPLOAD_DIR, exist_ok=True)
        with open(os.path.join(UPLOAD_DIR, prepared_name(hashlib.sha1(data).hexdigest()[:8], inp_name, True)), "wb") as fh:
            fh.write(data)
    st.session_state[G.PREPARED_KEY] = {"name": inp_name, "data": data,
                                        "from": f"the GIS layer {b.layer.name if b.layer else b.name}"
                                                if not b.wntr_schema else f"the WNTR GIS files {b.name}"}
    st.session_state[G.PREPARED_OPEN] = True
    try:
        st.switch_page(MAIN_PAGE)
    except StreamlitAPIException:          # the page run on its own (no page list): say where to go
        st.success("The model is ready on the main page: open 'ResidualMap' in the sidebar and choose "
                   f"'{G.PREPARED_OPTION}'.")


def unit_override(lay) -> str | None:
    """The coordinates' unit for a layer whose file names one that may be wrong (it was refused as crowded, or its
    coordinates lie in longitude and latitude ranges): the file's own unit first, then the others, degrees among them
    (the final verification of task 17)."""
    own = f"as the file says ({G.UNIT_WORDS.get(lay.crs.get('unit'), 'its unit')})"
    choices = {own: None, **{k: v for k, v in UNIT_CHOICES.items() if v is not None}}
    u = st.selectbox("Unit of the coordinates", list(choices), key=f"gis_unit_override_{lay.name}",
                     help=f"The file says {lay.crs.get('how')}, but its coordinates may not be in that unit (all of "
                          f"them in longitude and latitude ranges, or its line ends crowded together): choose another "
                          f"unit to override it.")
    return choices[u]


def gis_step(data: bytes, name: str) -> None:
    """Step 4: the layer read, its questions asked, its skeleton built and its readiness report shown, with the
    download and the hand-off.  Since the review of task 17 nothing here can raise (gis_read.build_from_layer never
    does), the hosted demo reads a zip's members within its upload limit and refuses a layer over its size before the
    lines are joined, and each layer is built once, and again when a question's answer or a unit chosen over the
    file's changes it."""
    unpacked = {"max_unpacked": GIS_UNPACKED} if GIS_UNPACKED else {}
    if data[:2] == b"PK" and G.wntr_schema_members(data, **unpacked) is not None:
        b = G.build_gis(data, name, **unpacked)
    else:
        try:
            layers, kind = G.read_layers(data, name, **unpacked)
        except G.GisRefused as e:
            st.warning(str(e))
            return
        except Exception as e:  # noqa: BLE001  (the last guard: a file the readers did not foresee)
            st.warning(f"The file could not be read as a GIS layer ({type(e).__name__}). Export the layer again.")
            return
        line_layers = [i for i, l in enumerate(layers) if l.n_lines]
        if len(line_layers) > 1:        # by position: two layers may share a name (the review of task 17)
            pick = st.selectbox("Layer", line_layers, index=0, format_func=lambda i: layers[i].name)
            lay = layers[pick]
        else:
            lay = layers[line_layers[0]] if line_layers else layers[0]
        if CLOUD:                        # the hosted demo's limit, before any line is joined
            too_big = cloud.gis_too_big(lay.n_lines, sum(len(p) for f in lay.features for p in f.parts))
            if too_big:
                st.error(too_big)
                return
        opts = G.GisOptions(layer=lay.name)
        unit_shown = False
        if lay.crs.get("axes_suspect"):        # the final verification of task 17: latitude written first?
            ax = st.selectbox("Order of the coordinates", list(AXES_CHOICES), key=f"gis_axes_{lay.name}",
                              help="The layer reads as longitude and latitude only with each coordinate's two "
                                   "numbers swapped: some files write latitude first.")
            opts.swap_axes = AXES_CHOICES[ax]
        if lay.crs.get("kind") == "unknown" or lay.crs.get("unit") not in (*TOPO.UNIT_M, "deg", "other"):
            u = st.selectbox("Unit of the coordinates", list(UNIT_CHOICES),
                             help=f"The layer does not say ({lay.crs.get('how')}). Lines are joined within a distance, "
                                  f"so the unit must be known.")
            opts.unit, unit_shown = UNIT_CHOICES[u], True
        elif lay.crs.get("kind") == "projected" and G.looks_lonlat(lay):
            opts.unit, unit_shown = unit_override(lay), True     # the final verification of task 17
        guess = G.guess_fields(lay)
        tag = hashlib.sha1(data).hexdigest()[:8] + lay.name     # a new file or layer starts from its own guesses
        with st.expander("Fields", expanded=False):
            st.caption("Guessed from the layer's field names: check them.")
            for role, label in FIELD_LABELS:
                options = [NO_FIELD] + list(lay.fields)
                g = guess.get(role)
                pickf = st.selectbox(label, options, index=options.index(g) if g in options else 0,
                                     key=f"gis_{role}_{tag}")
                opts.fields[role] = None if pickf == NO_FIELD else pickf
        kml = lay.source == "kml"
        default_ft = TOPO.DEFAULT_TOL_M["kml" if kml else "gis"] / TOPO.FT
        tol_ft = st.number_input("Join line ends closer than (ft)", 0.1, 200.0, round(default_ft, 1), 0.5,
                                 help="The default is 6 ft for a GIS layer and 10 m (32.8 ft) for lines drawn by hand in "
                                      "Google Earth. The table below shows what half and twice this distance would do.")
        opts.tolerance_m = tol_ft * TOPO.FT
        b = G.build_from_layer(lay, opts, name, kind, layers)
        changed = False
        for q in b.questions:
            if q["key"] == "bare_unit":
                du = st.selectbox("Unit of the pipe sizes", list(DIAMETER_UNITS), help=q["text"])
                if DIAMETER_UNITS[du]:
                    opts.bare_unit, changed = DIAMETER_UNITS[du], True
            elif q["key"] == "default_diameter_in":
                opts.default_diameter_in = st.number_input("Size of every pipe (in)", 1.0, 120.0, 8.0, 1.0,
                                                           help=q["text"])
                changed = True
            elif q["key"] == "to_scale":
                ans = st.radio("Is the layer drawn to scale?", ["not sure", "yes", "no"], horizontal=True,
                               help=q["text"])
                opts.to_scale = {"yes": True, "no": False}.get(ans)
                kl = st.text_area("Pipes whose true length you know (optional): one per line, 'pipe ID, length in ft'",
                                  "", key="gis_known_lengths")
                opts.known_lengths_ft = _known_lengths(kl)
                changed = changed or opts.to_scale is not None or bool(opts.known_lengths_ft)
        if changed:
            b = G.build_from_layer(lay, opts, name, kind, layers)
        if b.offer_unit and not unit_shown:  # refused as crowded: the file's unit may be wrong (the final
            opts.unit = unit_override(lay)    # verification of task 17: no control was offered to change it)
            if opts.unit:
                b = G.build_from_layer(lay, opts, name, kind, layers)
    if b.refusal:
        st.warning(b.refusal)
        return
    for n in b.notes:
        st.caption(n)
    blocking = b.blocking_questions()
    if blocking:
        for q in blocking:
            st.info(q["text"])
        return
    if b.wntr_schema:
        st.markdown(f"**WNTR's own GIS files:** the model is read from them. {' '.join(b.notes[-1:])}")
    else:
        r = b.net.report
        st.markdown(f"**The layer {b.layer.name}:** {b.layer.n_lines:,} lines become {r['pipes']:,} pipes and "
                    f"{r['nodes']:,} nodes, ends joined within {b.tolerance_m / TOPO.FT:.1f} ft "
                    f"({b.tolerance_m:.2f} m): {r['components']} connected part{'s' if r['components'] != 1 else ''}, "
                    f"{r['dead_ends']:,} dead ends, {r['tees_split']:,} cut{'s' if r['tees_split'] != 1 else ''} "
                    f"where another line ends on one.")
        st.dataframe([{"join ends within (ft)": round(x["tolerance_ft"], 1), "pipes": x["pipes"], "nodes": x["nodes"],
                       "connected parts": x["components"], "dead ends": x["dead_ends"],
                       "cuts at a tee": x["tees_split"], "nodes split": x["nodes_split"],
                       "pieces dropped": x["dropped_pieces"]} for x in b.diagnostics], hide_index=True)
        for w in b.warnings:
            st.warning(w)
        if b.id_table:
            with st.expander(f"IDs changed for EPANET ({len(b.id_table)})"):
                st.dataframe([{"in the layer": x["layer_id"], "in the model": x["model_id"], "why": x["why"]}
                              for x in b.id_table], hide_index=True)
        vlay = G.valve_layer(layers, b.layer)
        valves = G.proposed_closed_valves(vlay, b.net, frame=b.frame) if vlay is not None else {"rows": []}
        if valves["rows"]:
            st.markdown(f"**Normally closed valves** (from the field '{valves['status_field']}'"
                        + (f" of the layer {vlay.name}" if vlay is not b.layer else "") + "): proposed for you to "
                        f"confirm. The model does not close them yet.")
            st.dataframe([{"valve": x["valve"], "status": x["status"], "at node": x["node"] or "(no node within the "
                                                                                           "distance)"}
                          for x in valves["rows"]], hide_index=True)
    rep, refusal = gis_readiness(b)
    if refusal:
        st.error(refusal)
        return
    st.markdown(rep.to_markdown())
    inp_name = f"{TOPO.clean_id(os.path.splitext(os.path.basename(name))[0])[0]}_skeleton.inp" if not b.wntr_schema \
        else f"{TOPO.clean_id(os.path.splitext(os.path.basename(name))[0])[0]}.inp"
    st.download_button("Download the model (.inp)", data=b.inp_text.encode("utf-8"), file_name=inp_name,
                       mime="text/plain", key="gis_download")
    if st.button("Open it on the main page", key="gis_hand_off",
                 help="The main page shows its readiness report; a skeleton stops there until its source, tanks, "
                      "pumps and demand are added, in EPANET on the downloaded file (a form for them here is not "
                      "built yet)."):
        hand_off(b, inp_name)


st.set_page_config(page_title="ResidualMap: bring your data", layout="wide")
st.title("Bring your data")
if CLOUD:      # task 20: the privacy notice at the top of the page
    st.info(cloud.notice(UPLOADS))
st.caption(cloud.get_started_intro(CLOUD, UPLOADS))

# ----------------------------------------------------------------------------- step 1
st.header("Step 1. Which software holds your network?")
# task 20: in cloud mode the guide's words that tell a visitor to upload, or name the repository's folders, are swapped
sections = {k: cloud.guide_section(k, v, CLOUD, UPLOADS) for k, v in guide_sections().items()}
choices = [k for k in sections if k != SAMPLES_SECTION]
choice = st.radio("Which software holds your network?", choices, index=0)
st.markdown(sections[choice])
st.caption("'Coming soon' means the app already recognises this kind of file and says what to export, but cannot read "
           "it yet. ResidualMap reads EPANET models, and every modelling package above can write one; "
           + ("step 4 below builds a skeleton model from a GIS layer." if UPLOADS else
              "the local version of ResidualMap builds a skeleton model from a GIS layer."))

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
if not UPLOADS:      # task 20: the hosted demo, examples only
    st.caption(cloud.CHECK_OFF)
    st.stop()
st.caption("Give it any file: a model, a drawing, a GIS layer, a grab log, a tap map. It says what the file is, by what "
           "is inside it, and what to send instead. A model gets the same readiness report as on the main page: the "
           "model is opened and run for one hour to make sure it works, but no map is built. Nothing is kept.")
species_label = st.radio("Disinfectant", list(DISINFECTANTS), index=0, horizontal=True, key="check_disinfectant")
up = st.file_uploader("Check a file", type=None, help="any file: it is identified by its content",
                      max_upload_size=cloud.UPLOAD_MB if CLOUD else None)
if up is not None and CLOUD and cloud.too_big(up.size):      # task 20: the hosted demo's limit
    st.error(cloud.too_big(up.size))
    st.stop()
if up is not None:
    data, name, species = up.getvalue(), up.name, DISINFECTANTS[species_label]
    s = sniff(data, name)
    st.subheader("What this file is")
    st.markdown(f"**{s.sentence()}**")
    if s.is_model:
        # under the git-ignored outputs/private/check_a_file in local mode, as before; in cloud mode under the system's
        # temporary directory, so nothing is written under outputs/ (task 20); deleted at once either way
        # and in cloud mode with the hosted demo's limits (cloud.intake_limits): a zip's model over the file limit is
        # refused before it is unpacked, and a model over the junction limit before WNTR parses it
        with cloud.temporary_folder("check_", CHECK_DIR) as d:
            path, rep = (cloud.prepare_upload(data, d, species, name, APP_ROOT) if CLOUD else
                         prepare_inp(data, d, species, name=name))
            refusal = cloud.limit_refusal(rep) if CLOUD else None
            if refusal:
                st.error(refusal)
            else:
                st.markdown(rep.to_markdown())
        if not refusal:
            st.caption("This is the readiness report the main page shows when you upload this model there.")
    elif s.category == "table" or s.kind == "text_other":     # text with a title line before its header too
        try:
            check_table(data, name, species, s)
        except Exception as e:  # noqa: BLE001  (the last guard, as on the main page: a file the reader did not foresee)
            st.warning(f"The file could not be checked as a table ({type(e).__name__}). Save it again as CSV.")
    else:
        st.warning(f"ResidualMap reads EPANET models. {s.advice}")

# ----------------------------------------------------------------------------- step 4 (task 17)
st.header(GIS_HEADER)
st.caption(GIS_INTRO)
gis_up = st.file_uploader("GIS layer of your water mains", type=None, help=GIS_HELP, key="gis_layer",
                          max_upload_size=cloud.UPLOAD_MB if CLOUD else None)
if gis_up is not None and CLOUD and cloud.too_big(gis_up.size):      # task 20: the hosted demo's limit
    st.error(cloud.too_big(gis_up.size))
    st.stop()
if gis_up is not None:
    gis_step(gis_up.getvalue(), gis_up.name)
