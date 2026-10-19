"""
importers.py: read any model file safely (iteration 5, task 15), and make an imported model mean what the model
assumes (task 15b).

    python -m residualmap.importers FILE [FILE ...]               # what each file is, by its content, and what to do
    python -m residualmap.importers --prepare FILE [--out DIR]    # the model readiness report for an EPANET model

A utility sends what it has: an EPANET .inp saved by a vendor package (a byte-order mark, ';' comment lines first, a
Windows code page, UTF-16), the same file renamed .txt or kept as a .bak, an emailed zip, a CAD drawing, a GIS layer, a
vendor model database or a spreadsheet.  This module identifies a file by its content, never by its name, and turns
an EPANET model into a file WNTR reads exactly as the model means it, with a plain readiness report.

- sniff(data, name) reads the first 64 KB (the whole file for zips, SQLite, OLE and JSON, all in memory; the whole text
  when the first 64 KB looks like the start of an EPANET file but holds none of its node or link sections) and returns
  a Sniff: kind, plain label, detail, confidence, whether the file name decided anything (by_extension), whether the
  name marks a backup copy, and advice in plain words.  The content decides the kind.  The name either adds a note
  (a .wtg.mdb or .wtg.sqlite name says which program uses the database, at high confidence) or breaks a tie the
  content cannot settle (at low confidence), and either way the result says so.  sniff never raises: a file it cannot
  read comes back as 'unknown', with advice.  Every test runs in time linear in the bytes it reads (no backtracking
  pattern over the file, csv.Sniffer on a bounded sample with no quote characters, SQLite queries under a step
  budget).
- decode_text(data) decodes a UTF-8 byte-order mark, then a UTF-16 one, then strict UTF-8, then cp1252 (any of the five
  bytes cp1252 leaves undefined is read as latin-1, the rest as cp1252), and says which.  The BOMs come first: a file
  with a UTF-8 BOM also decodes as plain UTF-8, with a leading U+FEFF that WNTR rejects (ENSyntaxError 201 at line 1).
- scan_inp_ids(text) checks the IDs in the text, before WNTR parses: WNTR keeps the later of two lines with one ID in a
  section without a word, loads a junction that shares a tank's ID and then the simulator fails, and stops on an ID
  of 32 characters with an AssertionError; EPANET itself counts bytes, so an accented ID that fits in a Windows code
  page can be too long once the copy is UTF-8.  Each of these blocks, with 'export with IDs, not labels'.  Pattern and
  curve IDs get the same length rule.
- prepare_inp(src, out_dir, disinfectant) returns (inp_path, Readiness).  A clean UTF-8 file given as a path comes back
  as the same path, untouched (so the app's upload name and grid cache do not change); clean bytes are written as
  <sha1 of the bytes>[:8]_<name>.  A file that needed a fix is written as <sha1 of the new text>[:8]_<stem>.inp in
  out_dir, which callers keep git-ignored (outputs/app_uploads, outputs/private/inputs).  The original is never
  changed.  A file that is not a model, or has a blocking problem, gives (None, Readiness) and leaves nothing in
  out_dir: a copy is written under a temporary name and put in place only after WNTR has read it, an entry point
  was found and EPANET ran the model's first hour (a test run, set up as the app's first run, in a temporary folder;
  WNTR reads a junction no pipe connects, and EPANET then stops).  An upload name too long for a file name once the
  digest is added is shortened (prepared_name).  Nothing in a file makes it raise: any failure after the file is read
  is a blocking problem.
- Task 15b, the model's physics, each a row of the readiness report (fix_reactions_text, dw_to_hw, _physics_report):
  reaction orders other than 1 are set to 1; per-pipe BULK and per-tank TANK lines are cleared and a LIMITING
  POTENTIAL set to 0 (per-pipe WALL lines are ignored, since the grid sets every wall rate); a Darcy-Weisbach model is
  converted to Hazen-Williams (each pipe's C matched to its headloss at its reference velocity) and kept only when the
  water age of the copy stays within DW_GATE_MEDIAN_H and DW_GATE_P95_H of the file's, and Chezy-Manning blocks; the
  start clock time is recorded and a weekly pattern warned about; with pressure-control valves every junction at
  elevation 0 is flagged and one elevation everywhere blocks; coordinates are typed as schematic, to scale in feet or
  metres, or missing; sources at nodes the tool does not dose (boosters), a source the dose overrides, and more than
  one inlet are warned about; pumps with no control, rule or speed pattern, no time patterns, zero demand and
  disconnected parts are warned about, and a WaterGEMS export's connector pipes reported; one C value on every pipe
  means the old-pipe factor cannot be learned (grid_for gives 'full_fixed_gamma'); and the run cost gives the grid's
  runs, its file and its temporary files, and minutes from one timed 7-day run, with a hard warning above
  LARGE_JUNCTIONS and size_gate's refusal when the build would leave under 1.5 GB free.  A model blocked now leaves
  no copy behind, not even a content-named copy an earlier version accepted.
- Task 15b also fixed three of task 15's leftovers in the sniff: a SWMM model whose SWMM-only sections come after the
  first 64 KB (the whole text is searched for them), a zip holding a model under another name (its members' contents
  are typed), and a UTF-8 file with a few stray bytes of a Windows code page (decode_text reads it as UTF-8 with those
  bytes as cp1252, and the report says the file mixes two encodings; since the second review the byte pairs a cp1252
  file makes by chance count for cp1252 and lines are counted, not characters: _utf8_kind, _mixed_utf8).
- Section headers and their names are read the way WNTR's reader reads them: after any whitespace (a form feed
  included), and a singular name such as [PIPE] is the [PIPES] section; the ID scan also stops at [END], as WNTR does.
- Zips are read in memory only, never extracted to disk.  A member with an absolute path or a '..' part is refused,
  and so is a zip whose members add up to more than 200 MB uncompressed.  The member that is read (the model) is
  refused if it is over 10 MB uncompressed and compressed more than 100 to 1; members that are only listed are not
  held to that ratio, because a shapefile's space-padded DBF packs past it at city size.  Files that macOS or Windows
  add to a zip (__MACOSX/, ._ copies, .DS_Store, Thumbs.db, desktop.ini) are left out of the listing.

Nothing here is applied inside simulate.load, so the bundled networks run bit for bit as before.  Uses the standard
library, WNTR and what WNTR itself needs (numpy, networkx).
"""
from __future__ import annotations

import argparse
import codecs
import contextlib
import csv
import hashlib
import io
import json
import math
import os
import re
import sqlite3
import struct
import sys
import tempfile
import time
import warnings
import zipfile
from dataclasses import asdict, dataclass, field

HEAD_BYTES = 64 * 1024                 # how much of a file the text and signature tests read
MAX_ID_BYTES = 31                      # EPANET's longest ID, in bytes of the UTF-8 copy (WNTR counts characters)
ZIP_MAX_TOTAL = 200_000_000           # uncompressed bytes in all members (200 MB)
ZIP_RATIO_MAX = 100                    # compression ratio refused, for the member that is read ...
ZIP_RATIO_MIN_SIZE = 10_000_000        # ... when it is over this many uncompressed bytes (10 MB)
ZIP_MAX_MEMBERS = 10000
JSON_MAX_BYTES = 50 * 1024 * 1024      # JSON files parsed whole up to this size, by their first 64 KB above it
TEXT_SCAN_MAX = 50 * 1024 * 1024       # an EPANET-looking text whose sections start after 64 KB is searched whole up to this
JSON_MAX_DEPTH = 64                    # nesting followed when coordinates are collected from a GeoJSON geometry
CSV_SNIFF_CHARS = 4096                 # csv.Sniffer sees this much, its quote characters taken out (_sniffer_sample)
XML_MAX_SKIP = 1000                    # comments, processing instructions and DOCTYPEs skipped before the root element
SQLITE_MAX_STEPS = 2000                # progress-handler calls (1,000 SQLite steps each) before a query is stopped
MAX_NAME_BYTES = 160                   # a prepared file's name, in UTF-8 bytes (a file name has 255, and the grid
                                       # cache adds up to about 60 for its own words)
TEST_RUN_ROW = "Test run"              # the readiness row of the EPANET run of the first hour
RUN_COST_ROW = "Run cost"              # the readiness row of the grid's size and time (task 15b)
LARGE_WARNING_START = "Large model:"   # the start of the warning for more than LARGE_JUNCTIONS junctions
FLOAT32_BYTES = 4                      # the grid cache stores ln C as float32: runs x 24 hours x junctions
CP1252_UNDEFINED = b"\x81\x8d\x8f\x90\x9d"   # the five bytes cp1252 leaves undefined


def _cp1252_gaps(err):
    """Error handler for cp1252: a byte cp1252 leaves undefined is read as latin-1 (the same code point)."""
    if isinstance(err, UnicodeDecodeError):
        return err.object[err.start:err.end].decode("latin-1"), err.end
    raise err


codecs.register_error("residualmap_cp1252_gaps", _cp1252_gaps)

# provenance of a readiness row
FROM_FILE, DERIVED, ASSUMED, MISSING = "from your file", "derived", "assumed", "missing, please supply"
SOURCES = (FROM_FILE, DERIVED, ASSUMED, MISSING)
DISINFECTANT_WORDS = {"free_chlorine": "free chlorine", "chloramine": "chloramine (total chlorine)"}

INP_SECTIONS = ("TITLE", "JUNCTIONS", "RESERVOIRS", "TANKS", "PIPES", "PUMPS", "VALVES", "TAGS", "DEMANDS", "STATUS",
                "PATTERNS", "CURVES", "CONTROLS", "RULES", "ENERGY", "EMITTERS", "QUALITY", "SOURCES", "REACTIONS",
                "MIXING", "TIMES", "REPORT", "OPTIONS", "COORDINATES", "VERTICES", "LABELS", "BACKDROP", "END",
                "LEAKAGE")
CORE_SECTIONS = {"JUNCTIONS", "PIPES", "RESERVOIRS", "TANKS"}
# section headers that only a SWMM model (sewers, storm drains) or an EPANET-MSX reaction file has; both also use some
# EPANET headers ([JUNCTIONS], [PIPES], [TANKS]), so they are told apart before the EPANET test
# (only names no EPANET file and few other text files use: a generic one like [MAP] or [STORAGE] decides nothing)
SWMM_SECTIONS = ("OUTFALLS", "CONDUITS", "XSECTIONS", "SUBCATCHMENTS", "SUBAREAS", "RAINGAGES", "ORIFICES", "WEIRS",
                 "DIVIDERS", "TRANSECTS", "INFILTRATION", "DWF", "EVAPORATION")
MSX_SECTIONS = ("SPECIES", "COEFFICIENTS", "TERMS", "PARAMETERS", "DIFFUSIVITY")   # [SPECIES], which MSX requires, decides
# a header starts a line, after any whitespace (WNTR strips every kind, a form feed included), once the line ends are
# made '\n' (so a classic Mac CR-only file, or CR-CR-LF line ends from a double conversion, are read like any other);
# it ends at whitespace, a ';' comment or the line end.  The name is matched in Python (_canonical_section), the way
# WNTR's reader does it, so the pattern only finds the bracketed word: linear, with no alternation to backtrack.
_INP_HEADER = re.compile(r"(?m)^[^\S\n]*\[([A-Za-z]+)\](?=[^\S\n]|;|$)")
NODE_SECTIONS = ("JUNCTIONS", "RESERVOIRS", "TANKS")
LINK_SECTIONS = ("PIPES", "PUMPS", "VALVES")
NAMED_SECTIONS = ("PATTERNS", "CURVES")   # their first token is a pattern or curve ID, repeated on continuation lines


def _canonical_section(name: str) -> str:
    """A section name as WNTR's .inp reader reads it: upper case, and a name with a missing or an extra S
    ([JUNCTION], [PIPE], [TANK], [RESERVOIR], [OPTION], [TITLES]) is the EPANET section it stands for.  A SWMM or MSX
    name, or any other, comes back as it is (upper case)."""
    n = name.strip().upper()
    if n in INP_SECTIONS:
        return n
    if n + "S" in INP_SECTIONS:
        return n + "S"
    if n.endswith("S") and n[:-1] in INP_SECTIONS:
        return n[:-1]
    return n


_KNOWN_SECTIONS = frozenset(INP_SECTIONS + SWMM_SECTIONS + MSX_SECTIONS)


def _section_headers(text: str) -> set[str]:
    """The section names (upper case, as WNTR reads them) of every header line in the text, whatever its line ends."""
    flat = text.replace("\r\n", "\n").replace("\r", "\n")
    found = {_canonical_section(m.group(1)) for m in _INP_HEADER.finditer(flat)}
    return found & _KNOWN_SECTIONS

DWG_VERSIONS = {"AC1.2": "R1.2", "AC1.40": "R1.40", "AC1.50": "R2.05", "AC2.10": "R2.10", "AC1001": "R2.22",
                "AC1002": "R2.5", "AC1003": "R2.6", "AC1004": "R9", "AC1006": "R10", "AC1009": "R11 and R12",
                "AC1012": "R13", "AC1014": "R14", "AC1015": "AutoCAD 2000", "AC1018": "AutoCAD 2004",
                "AC1021": "AutoCAD 2007", "AC1024": "AutoCAD 2010", "AC1027": "AutoCAD 2013",
                "AC1032": "AutoCAD 2018"}
OLE_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
SQLITE_MAGIC = b"SQLite format 3\x00"
GPKG_APP_IDS = {0x47504B47: "GPKG", 0x47503130: "GP10", 0x47503131: "GP11"}
BINARY_DXF = b"AutoCAD Binary DXF\r\n\x1a\x00"
EPANET_BIN_MAGIC = 516114521
SHP_FILE_CODE = 9994

_EPANET_EXPORT = ("If a hydraulic model of this system exists, send its EPANET export instead: in WaterCAD or WaterGEMS, "
                  "make the scenario you want active and use File, Export, EPANET, with IDs rather than labels; in "
                  "InfoWater, use its EPANET export; most modelling packages can write an EPANET .inp.")
_SKELETON_LATER = ("Building a skeleton model from a drawing or a GIS layer is planned, not built yet, and even then it "
                   "carries no demands, tanks, pumps or controls.")

# kind -> (category, plain label, advice).  Categories: model, other_model (a model of something else, or part of
# one), results, drawing, bim, vendor_model, gis, table, document, text, archive, refused, unknown.
KINDS: dict[str, tuple[str, str, str]] = {
    "epanet_inp": ("model", "an EPANET model (.inp)",
                   "It is checked below and then run."),
    "epanet_inp_zip": ("model", "a zip holding an EPANET model",
                       "The model is read from the zip in memory (the zip is not unpacked) and saved as a prepared "
                       "copy."),
    "swmm_inp": ("other_model", "a SWMM sewer or stormwater model (.inp)",
                 "SWMM models sewers and storm drains, not a pressurised drinking-water network, although its files "
                 "also end in .inp. Send the water distribution model's EPANET .inp instead."),
    "epanet_msx": ("other_model", "an EPANET-MSX reaction file (.msx)",
                   "It holds water quality reactions, not the pipe network. Send the network's EPANET .inp: "
                   "ResidualMap sets and calibrates the chlorine decay itself."),
    "epanet_bin": ("results", "an EPANET binary results file (.bin)",
                   "It holds the results of a run, not the model. Send the .inp file that produced it."),
    "dxf_ascii": ("drawing", "a CAD drawing (DXF, text)",
                  "A drawing holds lines, not a hydraulic model. " + _EPANET_EXPORT + " " + _SKELETON_LATER),
    "dxf_binary": ("drawing", "a CAD drawing (binary DXF)",
                   "A drawing holds lines, not a hydraulic model. Save it again as a text (ASCII) DXF. " + _EPANET_EXPORT
                   + " " + _SKELETON_LATER),
    "dwg": ("drawing", "an AutoCAD drawing (DWG)",
            "DWG is a closed format that cannot be read here: save it as DXF (AutoCAD: Save As, DXF) or export the water "
            "mains as a GIS layer. " + _EPANET_EXPORT),
    "dgn_v7": ("drawing", "a MicroStation drawing (DGN V7)",
               "Export it as DXF or as a GIS layer. " + _EPANET_EXPORT),
    "dgn_v8": ("drawing", "a MicroStation drawing (DGN V8)",
               "DGN V8 needs MicroStation itself to convert: export it as DXF (File, Export, DXF) or as a GIS layer. In a "
               "WaterGEMS for MicroStation project the drawing is not the model. " + _EPANET_EXPORT),
    "ifc": ("bim", "a BIM model (IFC)",
            "IFC describes buildings and plant equipment, not the hydraulics of a pipe network. " + _EPANET_EXPORT
            + " Otherwise send the water mains as a GIS layer."),
    "pdf": ("document", "a PDF (a map book, drawing set or report)",
            "A PDF map cannot be read as a network. " + _EPANET_EXPORT + " Otherwise send the water mains as a GIS "
            "layer."),
    "access_db": ("vendor_model", "a Microsoft Access (Jet or ACE) database",
                  "Older WaterCAD and WaterGEMS versions keep the model in an Access database (.wtg.mdb), whose layout is "
                  "not published, so it is not read here. In WaterCAD or WaterGEMS use File, Export, EPANET and send "
                  "the .inp."),
    "sqlite_model": ("vendor_model", "a SQLite database that is not a GeoPackage",
                     "If it came from WaterGEMS or WaterCAD (a .wtg.sqlite file) it is the model database, whose layout "
                     "is not published, so it is not read here: in WaterGEMS or WaterCAD use File, Export, EPANET and "
                     "send the .inp."),
    "wtg_settings": ("vendor_model", "a WaterCAD or WaterGEMS project file (.wtg)",
                     "The .wtg file holds project settings; the model is in the .wtg.sqlite database next to it, and "
                     "neither is read here. In WaterCAD or WaterGEMS use File, Export, EPANET and send the .inp."),
    "vendor_model_zip": ("vendor_model", "a zip of WaterCAD or WaterGEMS model files",
                         "These files are not read here. In WaterCAD or WaterGEMS use File, Export, EPANET and send the "
                         ".inp."),
    "geopackage": ("gis", "a GIS layer (GeoPackage)",
                   "A GIS layer holds pipes and their sizes, not a hydraulic model. " + _EPANET_EXPORT + " "
                   + _SKELETON_LATER),
    "geojson": ("gis", "a GIS layer (GeoJSON)",
                "A GIS layer holds pipes and their sizes, not a hydraulic model. " + _EPANET_EXPORT + " "
                + _SKELETON_LATER),
    "esri_json": ("gis", "a GIS layer (Esri JSON)",
                  "A GIS layer holds pipes and their sizes, not a hydraulic model. " + _EPANET_EXPORT + " "
                  + _SKELETON_LATER),
    "gml": ("gis", "a GIS layer (GML)",
            "A GIS layer holds pipes and their sizes, not a hydraulic model. " + _EPANET_EXPORT + " " + _SKELETON_LATER),
    "kml": ("gis", "a Google Earth file (KML)",
            "Lines drawn in Google Earth are not a hydraulic model. " + _EPANET_EXPORT + " " + _SKELETON_LATER),
    "kmz": ("gis", "a Google Earth file (KMZ)",
            "Lines drawn in Google Earth are not a hydraulic model. " + _EPANET_EXPORT + " " + _SKELETON_LATER),
    "shapefile_zip": ("gis", "a zipped shapefile",
                      "A GIS layer holds pipes and their sizes, not a hydraulic model. " + _EPANET_EXPORT + " "
                      + _SKELETON_LATER),
    "shp": ("gis", "a shapefile's geometry file (.shp)",
            "A .shp alone has no attributes: send it zipped with its .shx, .dbf and .prj. A GIS layer is not a hydraulic "
            "model. " + _EPANET_EXPORT),
    "file_gdb_zip": ("gis", "a zipped Esri file geodatabase (.gdb)",
                     "A file geodatabase cannot be read here: export the water mains as a shapefile or GeoPackage. "
                     + _EPANET_EXPORT),
    "wntr_json": ("model", "a WNTR model file (JSON)",
                  "Write it out as an EPANET .inp (wntr.network.write_inpfile) and send that."),
    "json_other": ("text", "a JSON file that is not a GIS layer or a model",
                   "It is not a model, a drawing or a GIS layer that ResidualMap recognises. Tell us which program made "
                   "it."),
    "xml_other": ("text", "an XML file",
                  "It is not a model, a drawing or a GIS layer that ResidualMap recognises. Tell us which program made "
                  "it."),
    "gpx": ("gis", "a GPS track (GPX)",
            "GPS points are useful for sample taps later; they are not a hydraulic model. " + _EPANET_EXPORT),
    "xlsx": ("table", "an Excel workbook (.xlsx)",
             "If it is a grab-sample log, type the samples into the table for now; uploading a sample log is planned. "
             "The pilot path reads grab logs saved as CSV."),
    "xls_legacy": ("table", "a legacy Excel workbook (.xls)",
                   "Save it as .xlsx or CSV. If it is a grab-sample log, the pilot path reads it as CSV."),
    "office_other": ("document", "a Word or PowerPoint file",
                     "It is not a model. " + _EPANET_EXPORT),
    "delimited_text": ("table", "a table in text form (CSV or similar)",
                       "If it is a grab-sample log, type the samples into the table for now; the pilot path reads grab "
                       "logs as CSV. It is not a model."),
    "text_other": ("text", "plain text with no EPANET sections",
                   "It has none of the EPANET section headers ([JUNCTIONS], [PIPES] and so on). If it should be a model, "
                   "open it in EPANET or the modelling package and save it again as .inp."),
    "zip_other": ("archive", "a zip with no model, drawing or GIS layer recognised",
                  "Send the model's EPANET export (.inp) on its own, or tell us what the zip holds."),
    "ole_other": ("unknown", "an old Microsoft Office or MicroStation file (OLE) of unknown type",
                  "Not sure what it is. Tell us which program made it. " + _EPANET_EXPORT),
    "refused_zip": ("refused", "a zip that is refused",
                    "Nothing in it was opened or unpacked. Send the model file on its own."),
    "empty": ("unknown", "an empty file", "The file has no content. Send it again."),
    "unknown": ("unknown", "a file of unknown type",
                "Not sure what this file is: it is not an EPANET model, a drawing, a GIS layer or a table that "
                "ResidualMap recognises. Tell us which program made it. " + _EPANET_EXPORT),
}


class ZipRefused(ValueError):
    """A zip that is not opened: an unsafe member path, or too large when unpacked."""


@dataclass
class Sniff:
    kind: str
    detail: str = ""
    confidence: str = "high"           # high, medium or low; low means 'not sure', and the advice asks the user
    by_extension: bool = False         # the file name decided the kind or added to the detail (the detail says so)
    backup: bool = False               # the name ends in .bak: a backup copy of whatever the content is
    encoding: str | None = None        # text files: how the text was decoded
    members: list[str] = field(default_factory=list)   # zips: the members that decided the kind

    @property
    def category(self) -> str:
        return KINDS[self.kind][0]

    @property
    def label(self) -> str:
        return KINDS[self.kind][1]

    @property
    def advice(self) -> str:
        return KINDS[self.kind][2]

    @property
    def is_model(self) -> bool:
        return self.kind in ("epanet_inp", "epanet_inp_zip")

    def sentence(self) -> str:
        """One plain sentence: what the file is, by its content."""
        what = f"a backup copy of {self.label}" if self.backup else self.label
        s = f"This file is {what}"
        if self.detail:
            s += f": {self.detail}"
        s += "."
        if self.confidence == "low":
            s = "Not sure. " + s
        return s

    def as_dict(self) -> dict:
        d = asdict(self)
        d.update(category=self.category, label=self.label)
        return d


# ----------------------------------------------------------------------------- text
def decode_text(data: bytes) -> tuple[str, str]:
    """(text, how): a UTF-8 byte-order mark (stripped), a UTF-16 one, UTF-16 without one (NUL bytes on every other
    position), strict UTF-8, UTF-8 with a few stray bytes of a Windows code page (task 15b: _mixed_utf8, when the
    file's valid UTF-8 characters are not the byte pairs a cp1252 file makes by chance), then cp1252, where any of the
    five bytes cp1252 leaves undefined is read as latin-1 and every other byte as cp1252 (so a euro sign or a curly
    quote elsewhere in the file still reads right).
    how is 'utf-8' only for plain UTF-8 with no BOM, the one case that needs no rewrite."""
    for bom, codec, how in ((b"\xef\xbb\xbf", "utf-8", "utf-8 with a byte-order mark"),
                            (b"\xff\xfe", "utf-16-le", "utf-16 (little-endian, with a byte-order mark)"),
                            (b"\xfe\xff", "utf-16-be", "utf-16 (big-endian, with a byte-order mark)")):
        if data.startswith(bom):
            try:
                return data[len(bom):].decode(codec), how
            except UnicodeDecodeError:
                return data[len(bom):].decode(codec, errors="replace"), how + ", with undecodable bytes replaced"
    guess = _utf16_without_bom(data[:4096])
    if guess:
        return data.decode(guess, errors="replace"), f"utf-16 ({guess[-2:]}, no byte-order mark)"
    try:
        return data.decode("utf-8"), "utf-8"
    except UnicodeDecodeError:
        pass
    mixed = _mixed_utf8(data)
    if mixed is not None:
        text, n = mixed[:2]
        return text, f"{MIXED_PREFIX} {n} stray byte{'s' if n != 1 else ''} read as cp1252 (Windows Western European)"
    try:
        return data.decode("cp1252"), "cp1252 (Windows Western European)"
    except UnicodeDecodeError:
        n = sum(data.count(bytes([b])) for b in CP1252_UNDEFINED)
        return (data.decode("cp1252", errors="residualmap_cp1252_gaps"),
                f"cp1252 (Windows Western European), with {n} byte{'s' if n != 1 else ''} it leaves undefined read as "
                f"latin-1")


MIXED_PREFIX = "utf-8 mixed with"       # decode_text's word for UTF-8 with a few bytes of a Windows code page
_STRAY = re.compile("[\udc80-\udcff]")   # a byte that is not UTF-8, as surrogateescape decodes it
_NON_ASCII = re.compile("[^\x00-\x7f]")
_LINE_BREAK = re.compile("\r\n|\r|\n")


def _cp1252_byte(m: re.Match) -> str:
    return bytes([ord(m.group(0)) - 0xDC00]).decode("cp1252", errors="residualmap_cp1252_gaps")


def _utf8_kind(t: str, i: int) -> str:
    """How the valid multi-byte UTF-8 character t[i] reads in a file that also holds bytes that are not UTF-8 (task
    15b's second review).  A file in cp1252 throughout holds such characters only by chance, as an accented capital
    (bytes C2 to DF) followed by a cp1252 symbol or punctuation byte (80 to BF), and the character those two bytes make
    in UTF-8 says which it is:
    'strong': a lowercase accented letter (U+00E0 to U+00FF, but not U+00F0), a curly quote, a dash or a bullet
      (U+2010 to U+203A) or the euro sign.  In cp1252 these bytes read as 'Ã' or 'â€' and a symbol ('Ã©', 'â€™'),
      which no text in cp1252 writes, so one of them shows the file holds UTF-8.
    'plausible': another character a model's text uses (U+00A0 to U+017F, the Latin-1 and Latin Extended-A letters,
      the degree sign and the no-break space; the Romanian letters with a comma below; the Greek letters; any other
      character of three or four bytes), or a letter of Cyrillic to NKo (U+0400 to U+07FF) next to another letter of
      those scripts or Greek, as in a word.
    'chance': anything else, the byte pairs a cp1252 file makes by chance: the C1 controls (U+0080 to U+009F, 'Â' and
      a quote or dash), Latin Extended-B, the IPA and the combining marks (U+0180 to U+037F, so 'É', 'È', 'Ç', 'Ì' or
      'Í' and a quote, dash, degree sign or no-break space: 'CITÉ’' is the bytes C9 92, U+0252), and a Cyrillic to NKo
      letter alone between Latin letters ('Ó' or 'Ñ' and a dash)."""
    o = ord(t[i])
    if (0xE0 <= o <= 0xFF and o != 0xF0) or 0x2010 <= o <= 0x203A or o == 0x20AC:
        return "strong"
    if 0xA0 <= o <= 0x17F or 0x218 <= o <= 0x21B or 0x386 <= o <= 0x3CE or o >= 0x800:
        return "plausible"
    if 0x400 <= o <= 0x7FF and any(0 <= k < len(t) and 0x370 <= ord(t[k]) <= 0x7FF for k in (i - 1, i + 1)):
        return "plausible"
    return "chance"


def _mixed_utf8(data: bytes) -> tuple[str, int, list[int], int, int] | None:
    """A UTF-8 file with a few stray bytes of a Windows code page (task 15b; one pasted degree sign or accented letter
    is enough to fail strict UTF-8): (text, the number of bytes read as cp1252, the first three lines that hold them,
    the number of lines that hold them, the number of chance pairs among them), or None to decode the file as cp1252
    throughout, as task 15 did.
    Each valid UTF-8 character is weighed with _utf8_kind.  The file is read as UTF-8 when it holds a 'strong' one, or
    when the lines that hold 'plausible' ones are at least as many as the lines that hold stray bytes or chance pairs.
    Then every stray byte and every chance pair is read as cp1252 (any of the five bytes cp1252 leaves undefined as
    latin-1) and every other UTF-8 character is kept.  Decoding such a file as cp1252 throughout, as task 15 did,
    turned every UTF-8 accent into two characters: an accented ID changed, a letter whose second byte is 0xA0 became a
    no-break space that WNTR splits on, and the model was blocked with WNTR's or EPANET's error rather than a word
    about the encoding.
    Lines, not occurrences, are counted (the second review): an ID is written on its own line, at each pipe end and
    in the coordinates, so counting occurrences let one ID with a chance pair outweigh the stray bytes of a whole file
    ('CITÉ’15' with a degree sign in the title: 3 against 1).  In the two encodings a stray byte is read the same way
    and a chance pair too; only the 'plausible' and 'strong' characters read differently, so the rule weighs only the
    evidence about them.  What it still misreads: a file in cp1252 throughout whose only chance pairs make plausible
    letters (an 'Ã', 'Ä' or 'Å' before a quote or dash: 'KÅ“' reads as 'Kœ') on at least as many lines as its other
    bytes that are not UTF-8; and a UTF-8 file whose only accents are capitals or other 'plausible' characters, with
    stray bytes on more lines than they are, reads as cp1252 throughout."""
    t = data.decode("utf-8", errors="surrogateescape")
    if not _STRAY.search(t):
        return None
    import bisect
    breaks = [m.start() for m in _LINE_BREAK.finditer(t)]
    utf8_lines, other_lines, strong, chance, n_bytes = set(), set(), False, set(), 0
    for m in _NON_ASCII.finditer(t):
        i = m.start()
        line = bisect.bisect_left(breaks, i) + 1
        if "\udc80" <= t[i] <= "\udcff":
            other_lines.add(line)
            n_bytes += 1
            continue
        kind = _utf8_kind(t, i)
        if kind == "chance":
            chance.add(i)
            other_lines.add(line)
            n_bytes += len(t[i].encode("utf-8"))
        else:
            utf8_lines.add(line)
            strong = strong or kind == "strong"
    if not utf8_lines or not (strong or len(utf8_lines) >= len(other_lines)):
        return None

    def cp(m: re.Match) -> str:
        if m.start() in chance:
            return m.group(0).encode("utf-8").decode("cp1252", errors="residualmap_cp1252_gaps")
        return _cp1252_byte(m) if "\udc80" <= m.group(0) <= "\udcff" else m.group(0)
    return _NON_ASCII.sub(cp, t), n_bytes, sorted(other_lines)[:3], len(other_lines), len(chance)


def _utf16_without_bom(head: bytes) -> str | None:
    """'utf-16-le' or 'utf-16-be' when a text file without a byte-order mark has NUL bytes on every other position."""
    h = head[:4096]
    if len(h) < 16 or h.count(0) < len(h) // 4:
        return None
    odd, even = h[1::2].count(0), h[0::2].count(0)
    if odd > 0.9 * (len(h) // 2) and even < 0.05 * len(h):
        return "utf-16-le"
    if even > 0.9 * (len(h) // 2) and odd < 0.05 * len(h):
        return "utf-16-be"
    return None


def _head_text(data: bytes) -> tuple[str, str]:
    """The first HEAD_BYTES as text, with how they were decoded (a multi-byte character cut by the 64 KB limit does not
    turn UTF-8 into cp1252)."""
    head = data[:HEAD_BYTES + 4]
    guess = None if head.startswith((b"\xef\xbb\xbf", b"\xff\xfe", b"\xfe\xff")) else _utf16_without_bom(head)
    if guess:
        return head[:HEAD_BYTES].decode(guess, errors="replace"), f"utf-16 ({guess[-2:]}, no byte-order mark)"
    if len(data) <= HEAD_BYTES:
        return decode_text(data)
    text, how = decode_text(head[:HEAD_BYTES])
    if how.startswith(("cp1252", MIXED_PREFIX)):
        try:
            head[:HEAD_BYTES].decode("utf-8")
        except UnicodeDecodeError as e:
            if e.start >= HEAD_BYTES - 3:           # only the cut at the end fails
                return head[:e.start].decode("utf-8"), "utf-8"
    return text, how


def _is_binary(head: bytes) -> bool:
    h = head[:8192]
    if not h:
        return False
    ctrl = sum(1 for b in h if b < 9 or (13 < b < 32 and b != 27))
    return ctrl > 0.01 * len(h)


# ----------------------------------------------------------------------------- zips, in memory
def _unsafe_member(name: str) -> bool:
    n = name.replace("\\", "/")
    return n.startswith("/") or bool(re.match(r"^[A-Za-z]:", n)) or any(p == ".." for p in n.split("/"))


def zip_junk(name: str) -> bool:
    """A file that macOS or Windows adds to a zip by itself: anything under __MACOSX/ (Finder's Compress), an
    AppleDouble '._' copy, .DS_Store, Thumbs.db or desktop.ini.  It is left out of what a zip is said to hold."""
    parts = [p for p in name.replace("\\", "/").split("/") if p]
    base = parts[-1] if parts else ""
    return "__MACOSX" in parts or base.startswith("._") or base.lower() in (".ds_store", "thumbs.db", "desktop.ini")


def zip_listing(data: bytes, max_total: int = ZIP_MAX_TOTAL) -> list[zipfile.ZipInfo]:
    """The members of a zip held in memory, after the safety rules (no unsafe path, at most ZIP_MAX_MEMBERS members,
    at most max_total bytes unpacked in all); raises ZipRefused.  Nothing is unpacked, and the compression ratio is
    tested only for the member that is read (read_zip_member)."""
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            infos = zf.infolist()
    except Exception as e:  # noqa: BLE001  (BadZipFile, ValueError, struct.error, ...: any of them means not readable)
        raise ZipRefused(f"it is not a readable zip ({type(e).__name__}: {e})") from None
    if len(infos) > ZIP_MAX_MEMBERS:
        raise ZipRefused(f"it has {len(infos)} members, more than {ZIP_MAX_MEMBERS}")
    bad = [i.filename for i in infos if _unsafe_member(i.filename)]
    if bad:
        raise ZipRefused(f"a member's path leaves the zip's own folder ({bad[0]})")
    total = sum(i.file_size for i in infos)
    if total > max_total:
        raise ZipRefused(f"its members add up to {total / 1e6:.0f} MB unpacked, more than {max_total / 1e6:.0f} MB")
    return infos


def read_zip_member(data: bytes, member: str, limit: int = ZIP_MAX_TOTAL, ratio_max: float = ZIP_RATIO_MAX,
                    ratio_min_size: int = ZIP_RATIO_MIN_SIZE) -> bytes:
    """One member of a zip held in memory, read with a cap: a member over ratio_min_size bytes unpacked that is packed
    more than ratio_max to 1 is refused before it is opened, and a member whose real size exceeds its declared size
    (or the cap) is refused, so a zip that lies in its directory cannot fill memory."""
    infos = {i.filename: i for i in zip_listing(data)}
    if member not in infos:
        raise ZipRefused(f"no member {member}")
    i = infos[member]
    if i.file_size > ratio_min_size and i.file_size > ratio_max * max(i.compress_size, 1):
        raise ZipRefused(f"member {member} would unpack to {i.file_size / 1e6:.0f} MB, "
                         f"{i.file_size / max(i.compress_size, 1):.0f} times its packed size (more than "
                         f"{ratio_max:g} to 1)")
    cap = min(i.file_size, limit)
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as zf, zf.open(member) as fh:
            out = fh.read(cap + 1)
    except Exception as e:  # noqa: BLE001  (a corrupt offset raises ValueError 'negative seek value'; any failure refuses)
        raise ZipRefused(f"member {member} cannot be read ({type(e).__name__}: {e})") from None
    if len(out) > cap:
        raise ZipRefused(f"member {member} is larger than its zip directory says")
    return out


def _sniff_zip(data: bytes) -> Sniff:
    try:
        infos = zip_listing(data)
    except ZipRefused as e:
        return Sniff("refused_zip", detail=str(e), confidence="high")
    files = [i.filename for i in infos if not i.is_dir()]
    names = [n for n in files if not zip_junk(n)]
    low = [n.lower() for n in names]
    s = _sniff_zip_names(names, low)
    if s.kind == "zip_other" and names:
        models = _zip_models_by_content(data, names)
        if models:
            s = Sniff("epanet_inp_zip", detail="it holds " + ", ".join(models[:5]) + (" and more" if len(models) > 5
                                                                                     else "")
                      + f", an EPANET model by {'its' if len(models) == 1 else 'their'} content (the name does not "
                        f"end in .inp)", members=models)
    junk = len(files) - len(names)
    if junk:
        s.detail += (f"{'; ' if s.detail else ''}{junk} file{'s' if junk > 1 else ''} that macOS or Windows adds "
                     f"(__MACOSX, ._ copies, .DS_Store, Thumbs.db) left out")
    return s


ZIP_CONTENT_SNIFF_MAX = 50             # members whose content is read when no member's name says what the zip holds


def _zip_models_by_content(data: bytes, names: list[str]) -> list[str]:
    """Members that are EPANET models by their content, read in memory (task 15b: a zip holding a model renamed .txt
    was typed as 'a zip with nothing recognised').  Only when no member's name decided the kind; at most
    ZIP_CONTENT_SNIFF_MAX members, each under the same rules as the member that is read (read_zip_member: the
    ratio rule over 10 MB, the size cap), none over TEXT_SCAN_MAX unpacked, and a member that is itself a zip is not
    opened."""
    sizes = {}
    try:
        sizes = {i.filename: i.file_size for i in zip_listing(data)}
    except ZipRefused:
        return []
    out = []
    for n in names[:ZIP_CONTENT_SNIFF_MAX]:
        if sizes.get(n, TEXT_SCAN_MAX + 1) > TEXT_SCAN_MAX:
            continue
        try:
            b = read_zip_member(data, n)
        except ZipRefused:
            continue
        if b.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
            continue
        try:
            if _sniff_content(b, n).kind == "epanet_inp":
                out.append(n)
        except Exception:  # noqa: BLE001  (a member that cannot be read decides nothing)
            continue
    return out


def _sniff_zip_names(names: list[str], low: list[str]) -> Sniff:
    if not names:
        return Sniff("zip_other", detail="no files", confidence="medium")
    if "[content_types].xml" in low:
        if any(n.startswith("xl/") for n in low):
            return Sniff("xlsx", detail=f"{len(names)} parts, typed by its [Content_Types].xml and xl/ folder")
        return Sniff("office_other", detail="typed by its [Content_Types].xml")
    inp = [n for n, l in zip(names, low) if l.endswith(".inp") or l.endswith(".inp.bak")]
    if inp:
        return Sniff("epanet_inp_zip", detail="it holds " + ", ".join(inp[:5]) + (" and more" if len(inp) > 5 else ""),
                     members=inp)
    kml = [n for n, l in zip(names, low) if l.endswith(".kml")]
    if kml:
        return Sniff("kmz", detail="it holds " + ", ".join(kml[:3]), members=kml)
    gdb = sorted({re.sub(r"(?i)(\.gdb)/.*$", r"\1", n) for n in names if re.search(r"(?i)\.gdb/", n)})
    if gdb and any(l.endswith((".gdbtable", ".gdbtablx")) for l in low):
        return Sniff("file_gdb_zip", detail="geodatabase " + ", ".join(gdb), members=gdb)
    shp = [n for n, l in zip(names, low) if l.endswith(".shp")]
    if shp:
        layers = []
        for n in shp:
            base = n[:-4]
            have = {ext for ext in (".shx", ".dbf", ".prj", ".cpg") if (base + ext).lower() in low}
            miss = [ext for ext in (".shx", ".dbf", ".prj") if ext not in have]
            layers.append(base + (f" (missing {', '.join(miss)})" if miss else ""))
        return Sniff("shapefile_zip", detail=f"{len(shp)} layer{'s' if len(shp) > 1 else ''}: " + "; ".join(layers),
                     members=shp)
    wtg = [n for n, l in zip(names, low) if re.search(r"\.wtg(\.sqlite|\.mdb)?(\.bak)?$", l)]
    if wtg:
        return Sniff("vendor_model_zip", detail="it holds " + ", ".join(wtg[:4]), members=wtg)
    exts = sorted({os.path.splitext(l)[1] or "(no extension)" for l in low})
    return Sniff("zip_other", detail=f"{len(names)} files ({', '.join(exts[:8])})", confidence="medium")


# ----------------------------------------------------------------------------- OLE, SQLite
def ole_names(data: bytes, max_dir_sectors: int = 64) -> list[str]:
    """Names of the storages and streams in an OLE compound file's directory (the first max_dir_sectors directory
    sectors, followed through the first FAT sectors held in data).  [] when the directory is not in data."""
    if len(data) < 512 or not data.startswith(OLE_MAGIC):
        return []
    ssz = 1 << struct.unpack_from("<H", data, 30)[0]
    if ssz not in (512, 4096):
        return []
    first_dir = struct.unpack_from("<I", data, 48)[0]
    fat_secs = [s for s in struct.unpack_from("<109I", data, 76) if s < 0xFFFFFFFA]

    def sector(i):
        off = (i + 1) * ssz
        return data[off:off + ssz] if off + ssz <= len(data) else None

    fat: list[int] = []
    for s in fat_secs:
        b = sector(s)
        if b is None:
            break
        fat.extend(struct.unpack(f"<{ssz // 4}I", b))
    names, seen, s = [], set(), first_dir
    while s < 0xFFFFFFFA and s not in seen and len(seen) < max_dir_sectors:
        seen.add(s)
        b = sector(s)
        if b is None:
            break
        for k in range(0, ssz, 128):
            e = b[k:k + 128]
            nlen, typ = struct.unpack_from("<H", e, 64)[0], e[66]
            if typ in (1, 2, 5) and 2 <= nlen <= 64:
                names.append(e[:nlen - 2].decode("utf-16-le", errors="replace"))
        s = fat[s] if s < len(fat) else 0xFFFFFFFE
    return names


def _sqlite_open(data: bytes) -> sqlite3.Connection:
    """An in-memory copy of an uploaded SQLite database (sqlite3 deserialize; nothing written), with a step budget:
    any query that runs past SQLITE_MAX_STEPS progress calls is stopped (sqlite3.OperationalError 'interrupted'), so
    a crafted file cannot keep a query busy."""
    con = sqlite3.connect(":memory:")
    calls = [0]

    def budget() -> int:
        calls[0] += 1
        return 1 if calls[0] > SQLITE_MAX_STEPS else 0

    try:
        con.deserialize(data)
        con.set_progress_handler(budget, 1000)
    except BaseException:
        con.close()
        raise
    return con


def _sqlite_tables(data: bytes) -> dict[str, tuple[str, str]] | None:
    """{name: (type, sql)} for the tables and views of a SQLite database held in memory, read from sqlite_master only
    (no view is run), or None when it cannot be read."""
    try:
        con = _sqlite_open(data)
    except (sqlite3.Error, AttributeError, OverflowError, ValueError, MemoryError):
        return None
    with contextlib.closing(con):
        try:
            rows = con.execute("select name, type, sql from sqlite_master where type in ('table', 'view') "
                               "limit 10000").fetchall()
        except (sqlite3.Error, OverflowError, ValueError):
            return None
        return {str(n): (str(t), str(q or "")) for n, t, q in rows}


def _gpkg_layers(data: bytes, tables: dict) -> str:
    """'layers: ...' from gpkg_contents, read only when it is a plain table (a view or a virtual table would run code
    from the file), at most 50 rows, under the step budget."""
    typ, sql = tables.get("gpkg_contents", ("", ""))
    if not typ:
        return "layers not listed (no gpkg_contents table)"
    if typ != "table" or sql.lstrip().upper().startswith("CREATE VIRTUAL"):
        what = "a view" if typ == "view" else "a virtual table"
        return f"layers not listed (its gpkg_contents is {what}, not a plain table, so it is not read)"
    try:
        con = _sqlite_open(data)
        with contextlib.closing(con):
            rows = con.execute("select table_name, data_type from gpkg_contents limit 50").fetchall()
    except (sqlite3.Error, AttributeError, OverflowError, ValueError, MemoryError):
        return "layers not readable"
    layers = [f"{t} ({d})" for t, d in rows]
    return f"layers: {', '.join(layers) if layers else 'none listed'}"


# ----------------------------------------------------------------------------- JSON and XML
_GEOJSON_TYPES = {"FeatureCollection", "Feature", "Point", "MultiPoint", "LineString", "MultiLineString", "Polygon",
                  "MultiPolygon", "GeometryCollection"}


def _coords(geom, out: list, cap: int = 2000, max_depth: int = JSON_MAX_DEPTH) -> None:
    """Up to cap (x, y) pairs from a GeoJSON geometry, depth first and left to right, with an explicit stack and a
    nesting cap (so a crafted file cannot exhaust the recursion limit); members of the wrong type are skipped."""
    stack = [(geom, 0)]
    while stack and len(out) < cap:
        g, depth = stack.pop()
        if g is None or depth > max_depth:
            continue
        if isinstance(g, dict):
            subs = g.get("geometries")
            if isinstance(subs, list):
                stack.extend((x, depth + 1) for x in reversed(subs[:cap]))
            if "coordinates" in g:
                stack.append((g["coordinates"], depth + 1))
        elif isinstance(g, (list, tuple)):
            if len(g) >= 2 and all(isinstance(v, (int, float)) for v in g[:2]):
                try:
                    out.append((float(g[0]), float(g[1])))
                except OverflowError:
                    continue
            else:
                stack.extend((x, depth + 1) for x in reversed(g[:cap]))


def _crs_name(d: dict) -> str:
    """The GeoJSON crs member's properties.name, when every part of it has the type the specification gives."""
    crs = d.get("crs")
    props = crs.get("properties") if isinstance(crs, dict) else None
    name = props.get("name") if isinstance(props, dict) else None
    return name if isinstance(name, str) else ""


def _geojson_detail(d: dict) -> str:
    feats = d.get("features") if d.get("type") == "FeatureCollection" else None
    bad_feats = feats is not None and not isinstance(feats, list)
    if bad_feats:
        feats = []
    pts: list = []
    for f in (feats or [d])[:500]:
        _coords(f.get("geometry", f) if isinstance(f, dict) else None, pts)
    crs = _crs_name(d)
    m = re.search(r"EPSG:+(\d+)", crs)
    epsg = m.group(1) if m else ("4326" if "CRS84" in crs else None)
    lonlat = bool(pts) and all(abs(x) <= 180 and abs(y) <= 90 for x, y in pts)
    n = f"{len(feats)} feature{'s' if len(feats) != 1 else ''}, " if feats is not None else ""
    if bad_feats:
        n = "a features member that is not a list, "
    if epsg and epsg != "4326":
        where = f"projected coordinates (EPSG:{epsg}, from its crs member)"
        if lonlat:
            where += ", although the numbers fall in the lon/lat range"
    elif lonlat:
        where = "lon/lat coordinates (WGS84)"
    elif pts:
        where = "coordinates outside the lon/lat range, and no crs member says which system they are in"
    else:
        where = "no coordinates"
    return n + where


def _sniff_json(text: str, complete: bool) -> Sniff | None:
    try:
        d = json.loads(text) if complete else None
    except RecursionError:
        return Sniff("json_other", detail="nested too deeply to read", confidence="low")
    except ValueError:
        d = None
    if d is None:
        if complete:
            return None
        if re.search(r'"type"\s*:\s*"FeatureCollection"', text):
            return Sniff("geojson", detail="too large to read whole here; typed by its first 64 KB", confidence="medium")
        if re.search(r'"geometryType"\s*:\s*"esriGeometry', text):
            return Sniff("esri_json", detail="typed by its first 64 KB", confidence="medium")
        return None
    if isinstance(d, dict):
        keys = set(d)
        if str(d.get("version", "")).lower().startswith("wntr") or {"nodes", "links", "options"} <= keys:
            return Sniff("wntr_json", detail=f"WNTR version {d.get('version', 'not given')}")
        if str(d.get("geometryType", "")).startswith("esriGeometry") or (
                "spatialReference" in keys and "features" in keys):
            sr = d.get("spatialReference")
            if isinstance(sr, dict):
                wkid = sr.get("latestWkid") or sr.get("wkid")
            else:                      # some exports write the code alone
                wkid = sr if isinstance(sr, (int, str)) and not isinstance(sr, bool) else None
            feats = d.get("features")
            n = len(feats) if isinstance(feats, list) else 0
            return Sniff("esri_json", detail=f"{d.get('geometryType', 'features')}, {n} feature{'s' if n != 1 else ''}"
                                             + (f", spatial reference {wkid}" if wkid else ""))
        if d.get("type") in _GEOJSON_TYPES:
            return Sniff("geojson", detail=_geojson_detail(d))
        return Sniff("json_other", detail="keys " + ", ".join(sorted(map(str, keys))[:6]), confidence="medium")
    return Sniff("json_other", detail="a JSON list", confidence="medium")


_XML_NAME = re.compile(r"<([A-Za-z_][\w.\-]*)(?::([A-Za-z_][\w.\-]*))?")
_XML_SPACE = re.compile(r"\s*")


def _xml_root(text: str) -> tuple[str | None, str] | None:
    """(prefix, local name) of an XML text's root element, after any leading processing instructions, comments and
    DOCTYPE (its internal subset in square brackets included), or None.  A linear scan with str.find, at most
    XML_MAX_SKIP of them: a pattern with a repeated group of lazy matches (the first version) backtracks
    exponentially on a run of comments with no element after it."""
    i, n = 0, len(text)
    for _ in range(XML_MAX_SKIP + 1):
        i = _XML_SPACE.match(text, i).end()
        if text.startswith("<?", i):
            j = text.find("?>", i + 2)
            if j < 0:
                return None
            i = j + 2
        elif text.startswith("<!--", i):
            j = text.find("-->", i + 4)
            if j < 0:
                return None
            i = j + 3
        elif text[i:i + 9].upper() == "<!DOCTYPE":
            j, depth = i + 9, 0
            while j < n:
                c = text[j]
                if text.startswith("<!--", j):      # a comment in the internal subset
                    k = text.find("-->", j + 4)
                    if k < 0:
                        return None
                    j = k + 2
                elif c in "\"'":                    # a quoted literal may hold '>' or brackets
                    k = text.find(c, j + 1)
                    if k < 0:
                        return None
                    j = k
                elif c == "[":
                    depth += 1
                elif c == "]":
                    depth -= 1
                elif c == ">" and depth <= 0:
                    break
                j += 1
            if j >= n:
                return None
            i = j + 1
        else:
            m = _XML_NAME.match(text, i)
            if not m:
                return None
            return (m.group(1), m.group(2)) if m.group(2) else (None, m.group(1))
    return None


def _sniff_xml(text: str) -> Sniff | None:
    root = _xml_root(text)
    if root is None:
        return None
    prefix, local = root
    tag = f"{prefix}:{local}" if prefix else local
    low = local.lower()
    if low == "kml":
        return Sniff("kml", detail=f"{len(re.findall(r'<Placemark', text))} placemarks, "
                                   f"{len(re.findall(r'<LineString', text))} lines in the first 64 KB")
    if low == "haestaddatasource":
        v = re.search(r'Version="([^"]+)"', text)
        return Sniff("wtg_settings", detail="root element HaestadDataSource" + (f", version {v.group(1)}" if v else ""))
    if low == "gpx":
        return Sniff("gpx", detail="root element gpx")
    if (prefix or "").lower() in ("gml", "wfs", "ogr") or "www.opengis.net/gml" in text[:4096]:
        return Sniff("gml", detail=f"root element {tag}")
    return Sniff("xml_other", detail=f"root element {tag}", confidence="medium")


# ----------------------------------------------------------------------------- delimited text and DXF
_DELIMS = {",": "commas", ";": "semicolons", "\t": "tabs", "|": "vertical bars"}


def _sniffer_sample(sample: str) -> str:
    """What csv.Sniffer is given: the first CSV_SNIFF_CHARS characters with every quote character taken out.
    Sniffer's quote regexes (a lazy match between two quotes, across lines) are super-linear in the number of quote
    characters: 90 KB of a,"a,"... took 27.7 s, and even inside a 4 KB sample a line of a few hundred quote marks
    took seconds (the final verification), the time growing about as the fourth power of the count.  With no quote
    characters those regexes fail at once and Sniffer falls back to counting each delimiter per line, which is
    linear.  Quoting still matters, but only to the consistency test, which reads the lines with csv.reader and
    decides."""
    return sample[:CSV_SNIFF_CHARS].replace('"', "").replace("'", "")


def _sniff_table(text: str) -> Sniff | None:
    lines = [ln for ln in text.splitlines() if ln.strip()][:25]
    if len(lines) < 2:
        return None
    sample = "\n".join(lines)
    delim = None
    try:            # csv.Sniffer only orders the delimiters tried below; see _sniffer_sample for what it is given
        delim = csv.Sniffer().sniff(_sniffer_sample(sample), delimiters="".join(_DELIMS)).delimiter
    except csv.Error:   # the consistency test below reads all 25 lines with csv.reader's quoting, in linear time
        pass
    best = None
    for d in ([delim] if delim else []) + [x for x in _DELIMS if x != delim]:
        try:
            rows = list(csv.reader(lines, delimiter=d))
        except csv.Error:
            continue
        n0 = len(rows[0])
        if n0 < 2:
            continue
        same = sum(len(r) == n0 for r in rows) / len(rows)
        if same >= 0.8:
            best = (d, rows)
            break
    if best is None:
        return None
    d, rows = best
    cols = [c.strip() for c in rows[0]][:8]
    comma_decimals = d != "," and any(re.fullmatch(r"-?\d+,\d+", c.strip()) for r in rows[1:] for c in r)
    detail = f"separated by {_DELIMS[d]}; first row: {', '.join(cols)}" + ("; decimal commas" if comma_decimals else "")
    return Sniff("delimited_text", detail=detail, confidence="medium")


_DXF_SECTIONS = ("HEADER", "CLASSES", "TABLES", "BLOCKS", "ENTITIES", "OBJECTS", "THUMBNAILIMAGE")


def _is_ascii_dxf(text: str) -> tuple[bool, int]:
    """(is a DXF, number of 999 comment groups skipped): group code and value lines, comments first allowed, the first
    real group a 0 SECTION followed by a 2 with a section name."""
    lines = text.splitlines()[:400]
    i, comments = 0, 0
    while i + 1 < len(lines) and lines[i].strip() == "999":
        i, comments = i + 2, comments + 1
    if i + 3 >= len(lines):
        return False, comments
    ok = (lines[i].strip() == "0" and lines[i + 1].strip().upper() == "SECTION" and lines[i + 2].strip() == "2"
          and lines[i + 3].strip().upper() in _DXF_SECTIONS)
    return ok, comments


# ----------------------------------------------------------------------------- the sniff
def sniff(data: bytes, name: str = "") -> Sniff:
    """What a file is, by its content.  The name adds a note or breaks a tie (and marks a .bak as a backup copy); when
    it was used the result says so (by_extension, and the detail).  Never raises: a file whose reading fails in any
    test comes back as 'unknown', with advice, so the uploader always has something plain to say."""
    try:
        s = _sniff_content(data, name)
    except Exception as e:  # noqa: BLE001  (a crafted or broken file must not reach the user as a traceback)
        s = Sniff("unknown", detail=f"it could not be read ({type(e).__name__})", confidence="low")
    low = os.path.basename(name or "").lower()
    if low.endswith(".bak"):
        s.backup = True
    return s


def _only_comments(text: str) -> bool:
    """Every non-blank line starts with ';' (a vendor's comment header, which can run past 64 KB)."""
    lines = [ln.strip(" \t﻿") for ln in text.replace("\r", "\n").split("\n")]
    lines = [ln for ln in lines if ln]
    return bool(lines) and all(ln.startswith(";") for ln in lines)


def _by_name(s: Sniff, note: str) -> Sniff:
    s.by_extension = True
    s.detail = (s.detail + "; " if s.detail else "") + note + " (from the file name)"
    return s


def _sniff_content(data: bytes, name: str) -> Sniff:
    low = os.path.basename(name or "").lower()
    while low.endswith(".bak"):
        low = low[:-4]
    ext = os.path.splitext(low)[1]
    head = data[:HEAD_BYTES]
    if not data:
        return Sniff("empty")
    # binary signatures, strongest first
    if head.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return _sniff_zip(data)
    if head.startswith(OLE_MAGIC):
        names = ole_names(data)
        if any(n.startswith("Dgn~") for n in names):
            return Sniff("dgn_v8", detail="an OLE compound file with streams named Dgn~")
        if any(n in ("Workbook", "Book") for n in names):
            which = "Workbook" if "Workbook" in names else "Book"
            return Sniff("xls_legacy", detail=f"an OLE compound file with a {which} stream")
        if any(n in ("WordDocument", "PowerPoint Document") for n in names):
            return Sniff("office_other", detail="an OLE compound file with a " +
                         ("WordDocument" if "WordDocument" in names else "PowerPoint Document") + " stream")
        s = Sniff("ole_other", detail="an OLE compound file" + (f" with streams {', '.join(names[:5])}" if names else
                                                                 " whose directory could not be read"),
                  confidence="low")
        if ext == ".xls":
            return _by_name(Sniff("xls_legacy", detail=s.detail, confidence="low"), "the name ends in .xls")
        if ext == ".dgn":
            return _by_name(Sniff("dgn_v8", detail=s.detail, confidence="low"), "the name ends in .dgn")
        return s
    if head.startswith(SQLITE_MAGIC):
        app = struct.unpack_from(">I", head, 68)[0] if len(head) >= 72 else 0
        tables = _sqlite_tables(data)
        if app in GPKG_APP_IDS or (tables and "gpkg_contents" in tables):
            layers = _gpkg_layers(data, tables or {})
            why = (f"application id {GPKG_APP_IDS[app]}" if app in GPKG_APP_IDS else
                   f"a gpkg_contents {(tables or {}).get('gpkg_contents', ('table',))[0]}")
            return Sniff("geopackage", detail=f"{why}; {layers}")
        s = Sniff("sqlite_model", detail=f"{len(tables)} table{'s' if len(tables) != 1 else ''}" if tables is not None
                   else "tables not readable")
        if re.search(r"\.wtg\.sqlite$", low):
            return _by_name(s, "the name ends in .wtg.sqlite, which WaterGEMS and WaterCAD use for the model database")
        return s
    if head[4:19] in (b"Standard Jet DB", b"Standard ACE DB"):
        s = Sniff("access_db", detail=head[4:19].decode())
        if re.search(r"\.wtg\.mdb$", low):
            return _by_name(s, "the name ends in .wtg.mdb, which older WaterCAD and WaterGEMS versions use for the model")
        if low.endswith(".mdb") and "map" in low:
            return _by_name(s, "the name looks like a GIS personal geodatabase")
        return s
    pdf_at = head[:1024].find(b"%PDF-")
    if pdf_at >= 0:
        v = head[pdf_at + 5:pdf_at + 8].decode("latin-1", errors="replace")
        return Sniff("pdf", detail=f"PDF version {v}")
    if re.match(rb"AC(?:1\d{3}|[12]\.\d{1,2})", head[:6]):
        code = head[:6].decode("latin-1").rstrip("\x00")
        ver = DWG_VERSIONS.get(code) or DWG_VERSIONS.get(code[:5]) or "version not listed"
        return Sniff("dwg", detail=f"{ver} format ({code})")
    if head.startswith(BINARY_DXF):
        return Sniff("dxf_binary", detail="it starts with the binary DXF sentinel")
    if len(head) >= 4 and head[0] in (0x08, 0xC8) and head[1:4] == b"\x09\xfe\x02":
        return Sniff("dgn_v7", detail="3D" if head[0] == 0xC8 else "2D")
    if len(head) >= 8 and struct.unpack_from("<i", head, 0)[0] == EPANET_BIN_MAGIC:
        ver, nn, nt, nl = struct.unpack_from("<4i", head, 4) if len(head) >= 20 else (0, 0, 0, 0)
        return Sniff("epanet_bin", detail=f"EPANET version {ver}, {nn} nodes ({nt} tanks and reservoirs), {nl} links")
    if len(head) >= 36 and struct.unpack_from(">i", head, 0)[0] == SHP_FILE_CODE and \
            struct.unpack_from("<i", head, 28)[0] == 1000:
        st = struct.unpack_from("<i", head, 32)[0]
        kind = {1: "points", 3: "lines", 5: "polygons", 8: "multipoints", 11: "points with z", 13: "lines with z",
                15: "polygons with z", 21: "points with m", 23: "lines with m"}.get(st, f"shape type {st}")
        return Sniff("shp", detail=kind)
    # text
    bom = head.startswith((b"\xef\xbb\xbf", b"\xff\xfe", b"\xfe\xff"))
    utf16 = bom and head[:2] in (b"\xff\xfe", b"\xfe\xff")
    guess16 = None if bom else _utf16_without_bom(head)
    # trailing NUL bytes (a copy padded to a block size) do not make a text binary: WNTR stops reading at [END]
    if not (utf16 or guess16) and _is_binary(head.rstrip(b"\x00")):
        return Sniff("unknown", detail="binary content with no signature recognised", confidence="low")
    text, how = _head_text(data)
    found_all = _section_headers(text)
    deciding = CORE_SECTIONS | set(SWMM_SECTIONS) | {"SPECIES"}
    where = "in the first 64 KB"
    if (not found_all & deciding and HEAD_BYTES < len(data) <= TEXT_SCAN_MAX
            and (found_all or _only_comments(text))):
        # it starts like an EPANET file ([TITLE], [PATTERNS] or ';' comment lines) but its node and link sections
        # come later: a long title, a large pattern block first, a long comment header.  The file is in memory.
        full, how_full = decode_text(data)
        found_full = _section_headers(full)
        if found_full & deciding:
            found_all, how = found_full, how_full
            where = "in the whole file (its first 64 KB holds none of [JUNCTIONS], [PIPES], [RESERVOIRS] or [TANKS])"
    other_model = set(SWMM_SECTIONS) | {"SPECIES"}
    if (found_all & CORE_SECTIONS and not found_all & other_model and HEAD_BYTES < len(data) <= TEXT_SCAN_MAX
            and where == "in the first 64 KB"):
        # task 15b: a large SWMM model can hold only shared headers ([JUNCTIONS], [CONDUITS] later) in its first 64 KB,
        # and was typed as EPANET, then blocked with WNTR's syntax error; the whole text is searched for the
        # sections only SWMM or EPANET-MSX use (the file is in memory; the search is linear)
        full, how_full = decode_text(data)
        found_full = _section_headers(full)
        if found_full & other_model:
            found_all, how = found_full, how_full
            where = "in the whole file (after the first 64 KB, which holds only headers EPANET also uses)"
    swmm = sorted(found_all & set(SWMM_SECTIONS))
    if swmm:
        return Sniff("swmm_inp", detail=f"SWMM section{'s' if len(swmm) > 1 else ''} "
                                        f"{', '.join(f'[{x}]' for x in swmm[:4])} found {where}", encoding=how)
    msx = sorted(found_all & set(MSX_SECTIONS))
    if "SPECIES" in found_all and not found_all & {"JUNCTIONS", "RESERVOIRS"}:
        return Sniff("epanet_msx", detail=f"EPANET-MSX section{'s' if len(msx) > 1 else ''} "
                                          f"{', '.join(f'[{x}]' for x in msx[:4])} found {where}, and no "
                                          f"[JUNCTIONS] or [RESERVOIRS]", encoding=how)
    found = found_all & set(INP_SECTIONS)
    if found & CORE_SECTIONS:
        lead = text.lstrip("﻿ \t\r\n")
        extra = "; ';' comment lines come first" if lead.startswith(";") else ""
        n = f"{len(found)} EPANET section{'s' if len(found) > 1 else ''}"
        return Sniff("epanet_inp", detail=f"{n} found {where}, text in {how}{extra}",
                     confidence="high" if len(found) >= 2 else "medium", encoding=how)
    t = text.lstrip("﻿ \t\r\n")
    if t.startswith("ISO-10303-21"):
        m = re.search(r"FILE_SCHEMA\s*\(\s*\(\s*'([^']+)'", text)
        if m and m.group(1).upper().startswith("IFC"):
            return Sniff("ifc", detail=f"schema {m.group(1)}", encoding=how)
        return Sniff("ifc", detail="a STEP file (ISO-10303-21) with no IFC schema in its header", confidence="medium",
                     encoding=how)
    is_dxf, n999 = _is_ascii_dxf(text)
    if is_dxf:
        return Sniff("dxf_ascii", detail="group codes and values" + (f", {n999} comment group{'s' if n999 > 1 else ''} (999) first" if n999 else ""),
                     encoding=how)
    if t.startswith("<"):
        s = _sniff_xml(t)
        if s:
            s.encoding = how
            return s
    if t.startswith(("{", "[")):
        complete = len(data) <= JSON_MAX_BYTES
        s = _sniff_json(decode_text(data)[0].lstrip("﻿ \t\r\n") if complete else t, complete)
        if s:
            s.encoding = how
            return s
    s = _sniff_table(text)
    if s:
        s.encoding = how
        return s
    printable = sum(ch.isprintable() or ch in "\t\r\n" for ch in text[:8192]) / max(1, len(text[:8192]))
    if printable > 0.95:
        hint = " (the name ends in .inp, but the content has no EPANET sections)" if ext == ".inp" else ""
        return Sniff("text_other", detail=f"text in {how}{hint}", confidence="medium", encoding=how)
    return Sniff("unknown", detail="no signature or structure recognised", confidence="low")


# ----------------------------------------------------------------------------- IDs, before WNTR
@dataclass
class IdScan:
    counts: dict = field(default_factory=dict)        # section -> number of data lines
    blocking: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    n_coordinates: int = 0                            # data lines in [COORDINATES]
    coordinate_ids: set = field(default_factory=set)  # the node IDs those lines name


def _where(sections: list[str]) -> str:
    if len(set(sections)) == 1:
        return f"{len(sections)} times in [{sections[0]}]"
    return " and ".join(f"[{s}]" for s in sorted(set(sections)))


def scan_inp_ids(text: str) -> IdScan:
    """IDs in the node sections ([JUNCTIONS], [RESERVOIRS], [TANKS]) and the link sections ([PIPES], [PUMPS], [VALVES]),
    read from the text before WNTR parses it.  Blocking problems, each with 'export with IDs, not labels':
    - a duplicate within either group: WNTR keeps the later of two lines in one section without a word, and loads a
      junction that shares a tank's ID, after which the simulator fails;
    - an ID longer than EPANET allows, here or in [PATTERNS] and [CURVES]: EPANET counts bytes (31), and the prepared
      copy is UTF-8, so an accented ID of 31 characters in a Windows code page can be too long once copied (WNTR,
      which counts characters, stops on 32 with an AssertionError).
    IDs that differ only by upper and lower case are a warning (some EPANET versions treat them as one ID).  Also
    collects the node IDs of the [COORDINATES] lines, for the readiness report.
    The text is read the way WNTR's reader reads it, so no header form WNTR accepts gets past the scan: lines end at
    LF, CR or CR-LF only (str.splitlines would also split at a form feed or a vertical tab, which WNTR reads as
    spaces), a header's name is matched like WNTR's ([PIPE] is [PIPES]: _canonical_section), and [END] ends the
    file, as it does for WNTR."""
    out = IdScan()
    seen = {"node": {}, "link": {}}
    long_ids: list[tuple[str, str]] = []
    section = None
    for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        s = ln.split(";", 1)[0].strip()
        if not s:
            continue
        if s.startswith("["):
            section = _canonical_section(s.split("]", 1)[0].lstrip("["))
            if section == "END":
                break
            continue
        if section is None:
            continue
        tok = s.split()
        if section in NODE_SECTIONS or section in LINK_SECTIONS:
            group = "node" if section in NODE_SECTIONS else "link"
            out.counts[section] = out.counts.get(section, 0) + 1
            seen[group].setdefault(tok[0], []).append(section)
            if len(tok[0].encode("utf-8")) > MAX_ID_BYTES:
                long_ids.append((section, tok[0]))
        elif section in NAMED_SECTIONS:
            if len(tok[0].encode("utf-8")) > MAX_ID_BYTES:
                long_ids.append((section, tok[0]))
        elif section == "COORDINATES":
            out.n_coordinates += 1
            out.coordinate_ids.add(tok[0])
    words = {"node": "nodes (junctions, reservoirs, tanks)", "link": "links (pipes, pumps, valves)"}
    for group in ("node", "link"):
        dup = {i: secs for i, secs in seen[group].items() if len(secs) > 1}
        if dup:
            ex = ", ".join(f"{i} ({_where(secs)})" for i, secs in list(dup.items())[:5])
            harm = []
            if any(len(set(secs)) < len(secs) for secs in dup.values()):
                harm.append("WNTR would keep only the later of two lines in one section, without a word.")
            if any(len(set(secs)) > 1 for secs in dup.values()):
                harm.append("WNTR would load an ID listed in two sections (a junction that shares a tank's ID), and "
                            "the simulation would then fail.")
            out.blocking.append(f"{len(dup)} ID{'s are' if len(dup) > 1 else ' is'} used more than once among the "
                                f"{words[group]}: {ex}. {' '.join(harm)} Export the model with IDs, not labels.")
        fold: dict[str, set] = {}
        for i in seen[group]:
            fold.setdefault(i.upper(), set()).add(i)
        case = [sorted(v) for v in fold.values() if len(v) > 1]
        if case:
            out.warnings.append(f"{len(case)} pair{'s' if len(case) > 1 else ''} of {group} IDs differ only by upper "
                                f"and lower case ({', '.join(' and '.join(c) for c in case[:3])}); some EPANET versions "
                                f"treat them as one ID.")
    if long_ids:
        first = {}
        for sec, i in long_ids:
            first.setdefault(i, sec)
        accented = any(not i.isascii() for i in first)
        size = lambda i: f"{len(i)} characters" if i.isascii() else f"{len(i.encode('utf-8'))} bytes"
        ex = ", ".join(f"{i} ({size(i)}, in [{sec}])" for i, sec in list(first.items())[:3])
        out.blocking.append(f"{len(first)} ID{'s are' if len(first) > 1 else ' is'} longer than EPANET allows "
                            f"({MAX_ID_BYTES} characters"
                            + ("; EPANET counts bytes, and an accented letter takes 2 or more in the UTF-8 copy it reads"
                               if accented else "")
                            + f"): {ex}. Export the model with IDs, not labels.")
    return out


# ----------------------------------------------------------------------------- physics settings (task 15b)
# The decay model assumes first-order reactions, one bulk rate for every pipe and tank, and a Hazen-Williams C on every
# pipe (the old-pipe factor and the roughness axis read it).  A file can say otherwise in ways WNTR reads without a
# word; these settings and checks turn each into a fix, a block or a warning in the readiness report.
G_SI = 9.81                                # m/s2
NU_20C_M2S = 1.1e-5 * 0.3048 ** 2          # EPANET's kinematic viscosity of water at 20 C, 1.1e-5 ft2/s, in m2/s; the file's
                                           # VISCOSITY option is a multiple of it
HW_K, HW_Q_EXP, HW_D_EXP = 10.67, 1.852, 4.8704   # Hazen-Williams in SI, as the plan states it:
                                                  # h/L = 10.67 Q^1.852 / (C^1.852 D^4.8704)
V_REF_FLOOR_MS = 0.3                       # a pipe's reference velocity: the median |V| over the scored day, at least this
DW_GATE_MEDIAN_H, DW_GATE_P95_H = 1.0, 4.0  # the Darcy-Weisbach gate (plan open decision 3, lead decision 3), pre-registered
                                            # before any run: water age may move by at most these over the scored day's
                                            # junction-hours between the file and its converted copy
PRESSURE_VALVES = ("PRV", "PSV", "PBV")    # valves whose setting is a pressure: they need real elevations
SCHEMATIC_SPREAD = 1.5                     # 95th over 5th percentile of drawn over declared pipe length above this: schematic
FEET_PER_METRE = 1 / 0.3048
UNIT_RATIO_TOL = 0.1                       # a median ratio within 10% of 3.281 is feet, of 1.0 is metres
CONNECTOR_LENGTH_M, CONNECTOR_DIAM_M = 0.3048, 99 * 0.0254   # the 1 ft, 99 inch pipes a WaterGEMS export adds
CONNECTOR_TAG = "addedforinpconversion"
CONNECTOR_MIN_PIPES = 3                    # untagged pipes of that size reported only from this many (Net6 has one)
LARGE_JUNCTIONS = 3000                     # above this a grid build needs an explicit confirm (app checkbox, --allow-large)
EFFECTIVE_SPEEDUP = 4.0                    # the plan's assumed effective parallel speedup of a grid build, not the worker
                                           # count: its feasibility review measured 2.3 on Net3 and 4.9 on ky10 (11
                                           # workers).  Task 15b's pre-registered bar on Net3 missed (the build took
                                           # 2.10 times the estimate, against a bar of 2); the measured runs are in
                                           # outputs/intake/run_cost_timing.json.  Not restated: that is a decision.
TEMP_MB_PER_JUNCTION = 9.5 / 959           # one 7-day run's EPANET files (.inp, .rpt, .bin): 9.5 MB on ky4 (959
                                           # junctions); ky4, ky10 and Net6 are measured in outputs/intake/
                                           # task15b_physics.json (part run_cost, temp_files_one_run)
VOLATILE_COST = ("seconds_per_run", "estimate_s", "minutes", "workers", "temp_files_mb")   # vary by run or machine:
                                                                                           # not in the record
REACTION_ORDERS = ("BULK", "WALL", "TANK")
ALL_INP_SECTIONS = frozenset(INP_SECTIONS)


def _lines_keep(text: str) -> list[list[str]]:
    """The text as [line, line end] pairs, split at LF, CR or CR-LF only (as WNTR's reader splits it), so that joining
    them gives the text back exactly."""
    parts = re.split(r"(\r\n|\r|\n)", text)
    return [[parts[i], parts[i + 1] if i + 1 < len(parts) else ""] for i in range(0, len(parts), 2)]


def _data_lines(pairs: list[list[str]]):
    """(index, section, tokens) for every data line of the text (comments taken off), with section names read as
    WNTR reads them; stops at [END]."""
    section = None
    for i, (body, _) in enumerate(pairs):
        s = body.split(";", 1)[0].strip()
        if not s:
            continue
        if s.startswith("["):
            section = _canonical_section(s.split("]", 1)[0].lstrip("["))
            if section == "END":
                return
            continue
        if section is not None:
            yield i, section, s.split()


def _num(tok: str) -> float | None:
    try:
        v = float(tok)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def fix_reactions_text(text: str) -> tuple[str, dict]:
    """The [REACTIONS] settings the decay model cannot use as they are, fixed in the text (the rest of the file is
    left exactly as it is):
    - a reaction order other than 1 (bulk, wall or tank) is set to 1 (ky10 and Net6 declare 0, and then the bulk rate
      has no effect: the grid's members are all the same);
    - per-pipe BULK and per-tank TANK lines are turned into comments: they override the global bulk rate the grid
      sets, so those pipes and tanks would ignore the calibrated rate;
    - a LIMITING POTENTIAL other than 0 is set to 0 (it changes the first-order law).
    Per-pipe WALL lines and a ROUGHNESS CORRELATION are left in the copy: the grid sets every pipe's wall rate itself,
    so they are ignored, and the report says so.  Each changed line keeps a ';ResidualMap:' comment with what the file
    had.  Returns (text, info)."""
    pairs = _lines_keep(text)
    info = {"orders": {}, "orders_fixed": {}, "bulk_lines": [], "tank_lines": [], "wall_lines": 0,
            "limiting_potential": 0.0, "limiting_potential_fixed": None, "roughness_correlation": 0.0}
    for i, sec, tok in list(_data_lines(pairs)):
        if sec != "REACTIONS" or len(tok) < 3:
            continue
        k1, k2 = tok[0].upper(), tok[1].upper()
        if k1 == "ORDER" and k2 in REACTION_ORDERS:
            v = _num(tok[2])
            info["orders"][k2.lower()] = v
            if v != 1.0:
                info["orders_fixed"][k2.lower()] = tok[2]
                pairs[i][0] = f" ORDER {k2} 1 ;ResidualMap: the file had {tok[2]}; the decay model is first order"
        elif k1 in ("BULK", "TANK"):
            info["bulk_lines" if k1 == "BULK" else "tank_lines"].append(tok[1])
            pairs[i][0] = (";" + pairs[i][0].strip() + " ;ResidualMap: cleared, the calibrated bulk rate applies to "
                                                       "every pipe and tank")
        elif k1 == "WALL":
            info["wall_lines"] += 1
        elif k1 == "LIMITING":
            v = _num(tok[2])
            info["limiting_potential"] = v
            if v != 0.0:
                info["limiting_potential_fixed"] = tok[2]
                pairs[i][0] = f" LIMITING POTENTIAL 0 ;ResidualMap: the file had {tok[2]}"
        elif k1 == "ROUGHNESS":
            info["roughness_correlation"] = _num(tok[2])
    return "".join(a + b for a, b in pairs), info


def set_headloss_text(text: str, roughness: dict[str, str]) -> str:
    """The text with [OPTIONS] HEADLOSS set to H-W and each named pipe's roughness field replaced (the sixth field of
    its [PIPES] line); every other byte is left as it is."""
    pairs = _lines_keep(text)
    for i, sec, tok in list(_data_lines(pairs)):
        if sec == "OPTIONS" and tok[0].upper() == "HEADLOSS" and len(tok) >= 2:
            pairs[i][0] = f" Headloss            H-W ;ResidualMap: the file had {tok[1]}"
        elif sec == "PIPES" and len(tok) >= 6 and tok[0] in roughness:
            body, sep, comment = pairs[i][0].partition(";")
            a, b = [m.span() for m in re.finditer(r"\S+", body)][5]
            pairs[i][0] = body[:a] + roughness[tok[0]] + body[b:] + sep + comment
    return "".join(a + b for a, b in pairs)


def hw_equivalent_c(diameter_m: float, roughness_m: float, velocity_ms: float, nu_m2s: float) -> float:
    """The Hazen-Williams C that gives the Darcy-Weisbach headloss of a pipe at velocity V:
    C = [10.67 Q^1.852 (2 g D) / (D^4.8704 f V^2)]^(1/1.852), Q = V pi D^2 / 4, with f from Swamee-Jain,
    f = 0.25 / [log10(e / (3.7 D) + 5.74 / Re^0.9)]^2, Re = V D / nu (what EPANET uses above Re 4000)."""
    D, e, V = diameter_m, roughness_m, velocity_ms
    re_ = V * D / nu_m2s
    f = 0.25 / math.log10(e / (3.7 * D) + 5.74 / re_ ** 0.9) ** 2
    q = V * math.pi * D * D / 4
    return (HW_K * q ** HW_Q_EXP * 2 * G_SI * D / (D ** HW_D_EXP * f * V * V)) ** (1 / HW_Q_EXP)


def _run_7day(path: str, parameter: str = "AGE") -> tuple:
    """A 7-day run of the file, set up as simulate.load sets up every run, in a temporary folder that is removed:
    (water age in hours at every junction, pipe velocity in m/s), each hour of the scored day (the last)."""
    import wntr
    from .simulate import _last_day, load
    wn = load(path)
    wn.options.quality.parameter = parameter
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "run"))
    age = _last_day(res.node["quality"], wn.junction_name_list) / 3600.0
    vel = _last_day(res.link["velocity"], wn.pipe_name_list).abs()
    return age, vel


def dw_to_hw(read_path: str, text: str, work_dir: str) -> tuple[str, dict]:
    """A Darcy-Weisbach model converted to Hazen-Williams, and the gate's numbers.  Each pipe's C is the one that gives
    its Darcy-Weisbach headloss at its reference velocity (the median |V| over the file's own scored day, at least
    V_REF_FLOOR_MS), with nu from the file's VISCOSITY option.  The converted copy is run too, and water age is compared
    at every junction-hour of the scored day: the gate passes when the median change is at most DW_GATE_MEDIAN_H and
    the 95th percentile at most DW_GATE_P95_H.  Exact only at the reference velocity."""
    import numpy as np
    from .simulate import load
    wn = load(read_path)
    age0, vel0 = _run_7day(read_path)
    nu = float(wn.options.hydraulic.viscosity or 1.0) * NU_20C_M2S
    c, v_ref = {}, {}
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        v_ref[pn] = max(float(np.median(vel0[pn].values)), V_REF_FLOOR_MS)
        c[pn] = hw_equivalent_c(float(p.diameter), float(p.roughness), v_ref[pn], nu)
    new_text = set_headloss_text(text, {pn: f"{v:.3f}" for pn, v in c.items()})
    conv = os.path.join(work_dir, "converted.inp")
    with open(conv, "w", encoding="utf-8", newline="") as fh:
        fh.write(new_text)
    age1, _ = _run_7day(conv)
    d = np.abs(age1[age0.columns].values - age0.values).ravel()
    med, p95 = float(np.median(d)), float(np.percentile(d, 95))
    vals = np.array(list(c.values()))
    return new_text, {"pipes": len(c), "c_min": float(vals.min()), "c_max": float(vals.max()),
                      "c": {k: round(v, 3) for k, v in c.items()}, "v_ref": {k: round(v, 4) for k, v in v_ref.items()},
                      "age_change_median_h": med, "age_change_p95_h": p95,
                      "passes": med <= DW_GATE_MEDIAN_H and p95 <= DW_GATE_P95_H, "converted_path": conv}


def control_scan(text: str) -> dict:
    """What the [CONTROLS] and [RULES] act on and read, and the [TAGS] a WaterGEMS export adds, from the text:
    the link IDs a control or rule sets, the node IDs a control or rule reads (any NODE in a simple control; a NODE or
    JUNCTION whose PRESSURE or HEAD a rule's premise reads), and the number of each.  In a rule only the action
    clauses set a link (THEN and ELSE, and the AND lines after them); a link named in a premise (IF, and the AND and OR
    lines after it) is only read, so a pump that a rule reads and never sets is not controlled (found by the review:
    'IF PUMP 10 STATUS IS OPEN THEN PIPE 20 STATUS IS OPEN' left pump 10 closed for 7 days without a warning)."""
    pairs = _lines_keep(text)
    links, nodes, n_controls, n_rules, tags = set(), set(), 0, 0, []
    clause = None                         # in [RULES]: 'premise' or 'action', the clause the current line belongs to
    for _, sec, tok in _data_lines(pairs):
        up = [t.upper() for t in tok]
        if sec == "CONTROLS":
            n_controls += 1
            if up[0] == "LINK" and len(tok) > 1:
                links.add(tok[1])
            nodes.update(tok[k + 1] for k in range(len(up) - 1) if up[k] == "NODE")
        elif sec == "RULES":
            if up[0] == "RULE":
                n_rules += 1
                clause = None
                continue
            if up[0] == "IF":
                clause = "premise"
            elif up[0] in ("THEN", "ELSE"):
                clause = "action"
            elif up[0] == "PRIORITY":
                clause = None
            for k in range(len(up) - 1):
                if clause == "action" and up[k] in ("LINK", "PIPE", "PUMP", "VALVE"):
                    links.add(tok[k + 1])
                elif (clause == "premise" and up[k] in ("NODE", "JUNCTION") and k + 2 < len(up)
                      and up[k + 2] in ("PRESSURE", "HEAD")):
                    nodes.add(tok[k + 1])
        elif sec == "TAGS" and len(tok) >= 3 and any(t.lower() == CONNECTOR_TAG for t in tok[2:]):
            tags.append(tok[1])
    return {"links": links, "nodes": nodes, "n_controls": n_controls, "n_rules": n_rules, "connector_tags": tags}


def coordinate_scale(wn, have_xy: set) -> dict:
    """Drawn pipe length (with vertices) over declared length, for every pipe with both ends drawn: a spread (95th
    percentile over 5th) above SCHEMATIC_SPREAD is a schematic drawing; otherwise the median names the units (3.281
    per metre is feet, 1.0 is metres)."""
    import numpy as np
    if not have_xy:
        return {"status": "no coordinates", "pipes_measured": 0}
    r = []
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        if p.start_node_name not in have_xy or p.end_node_name not in have_xy or not p.length > 0:
            continue
        pts = [p.start_node.coordinates, *p.vertices, p.end_node.coordinates]
        d = sum(math.dist(pts[k], pts[k + 1]) for k in range(len(pts) - 1))
        if d > 0:
            r.append(d / p.length)
    if len(r) < 2:
        return {"status": "too few pipes drawn", "pipes_measured": len(r)}
    p5, p50, p95 = (float(x) for x in np.percentile(r, [5, 50, 95]))
    spread = p95 / p5
    out = {"pipes_measured": len(r), "ratio_p5": round(p5, 4), "ratio_median": round(p50, 4),
           "ratio_p95": round(p95, 4), "spread": round(spread, 2)}
    if spread > SCHEMATIC_SPREAD:
        return {"status": "schematic", "units": "unknown", **out}
    units = ("feet" if abs(p50 / FEET_PER_METRE - 1) <= UNIT_RATIO_TOL else
             "metres" if abs(p50 - 1) <= UNIT_RATIO_TOL else "unknown")
    return {"status": "to scale", "units": units, **out}


def grid_for(disinfectant: str, uniform_c: bool) -> str:
    """The simulator grid a file gets: chloramine's for chloramine; for free chlorine, 'full', or 'full_fixed_gamma'
    (the old-pipe factor held at 0) when every pipe has one C value and the factor cannot be learned."""
    if disinfectant == "chloramine":
        return "chloramine"
    return "full_fixed_gamma" if uniform_c else "full"


def disk_free_gb(path: str) -> float:
    """Free disk in GB (simgp.disk_free_gb); here so a saved check can stand in a value."""
    from .simgp import disk_free_gb as free
    return free(path)


def size_gate(cost: dict, path: str = ".", free_gb: float | None = None) -> tuple[bool, str]:
    """(ok, plain message): a grid build is refused when its file plus the temporary EPANET files of one run per worker
    would leave less than simgp.MIN_FREE_GB free on the volume of `path`."""
    from .simgp import MIN_FREE_GB
    free = disk_free_gb(path) if free_gb is None else float(free_gb)
    need = (cost["grid_file_mb"] + cost["temp_files_mb"]) / 1000
    left = free - need
    msg = (f"The grid needs about {cost['grid_file_mb']:g} MB for its file and up to {cost['temp_files_mb']:g} MB of "
           f"temporary EPANET files while it runs ({cost['workers']} at a time); this disk has {free:.2f} GB free, so "
           f"{left:.2f} GB would be left")
    if left < MIN_FREE_GB:
        return False, (msg + f", below the {MIN_FREE_GB:g} GB this tool keeps free. Free some disk space, or use a "
                             f"smaller zone of the model.")
    return True, msg + "."


# ----------------------------------------------------------------------------- the readiness report
@dataclass
class Readiness:
    name: str
    disinfectant: str = "free_chlorine"
    sha1: str = ""                                   # of the file as received
    sniff: dict = field(default_factory=dict)
    not_a_model: bool = False
    rows: list = field(default_factory=list)         # (item, value, source), source in SOURCES
    blocking: list = field(default_factory=list)
    fixes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    run_cost: dict = field(default_factory=dict)
    prepared: str | None = None                      # the prepared file's name (never its folder)
    same_file: bool = False                          # the prepared file is the file given, untouched
    flags: dict = field(default_factory=dict)        # task 15b: what the physics checks found, by check (bundled_dict)

    def add(self, item: str, value, source: str) -> None:
        assert source in SOURCES, source
        self.rows.append((item, str(value), source))

    @property
    def ok(self) -> bool:
        return not self.not_a_model and not self.blocking

    def advice(self) -> str:
        """For a file that is not a model: what it is and what to send instead, in plain words."""
        s = Sniff(**{k: v for k, v in self.sniff.items() if k in Sniff.__dataclass_fields__}) if self.sniff else None
        if s is None:
            return "This file could not be read."
        return f"{s.sentence()} ResidualMap reads EPANET models. {s.advice}"

    def summary(self) -> str:
        if self.not_a_model:
            return self.advice()
        if self.blocking:
            n = len(self.blocking)
            return f"{n} problem{'s' if n > 1 else ''} stop{'' if n > 1 else 's'} this model from running."
        if any(item == TEST_RUN_ROW for item, _, _ in self.rows):
            bits = ["WNTR reads it and EPANET runs its first hour without an error, and one full 7-day run (timed for "
                    "the run cost); the grid's runs come next. Its physics settings are checked in the rows below "
                    "(reaction orders, bulk lines, headloss formula, start time, valves and elevations, coordinates, "
                    "sources, pumps and controls, connectivity, pipe roughness)."]
        else:                          # prepare_inp(parse=False)
            bits = ["Prepared without a check: WNTR has not read it and EPANET has not run it; only the reaction "
                    "settings in its text were checked."]
        if self.fixes:
            bits.append(f"{len(self.fixes)} fix{'es' if len(self.fixes) > 1 else ''} applied to a copy; your file is "
                        f"not changed.")
        if self.warnings:
            bits.append(f"{len(self.warnings)} warning{'s' if len(self.warnings) > 1 else ''}.")
        return " ".join(bits)

    def to_markdown(self) -> str:
        out = [f"**Model readiness: {self.name}**", "", self.summary(), ""]
        if self.not_a_model:
            return "\n".join(out)
        if self.blocking:
            out += ["**Problems that stop the run:**"] + [f"- {b}" for b in self.blocking] + [""]
        if self.rows:
            out += ["| item | value | source |", "|---|---|---|"]
            out += [f"| {i} | {v.replace('|', '/')} | {s} |" for i, v, s in self.rows] + [""]
        if self.fixes and self.blocking:
            out += ["**Fixes that would be applied** (no copy is kept while a problem stops the run; your file is not "
                    "changed):"] + [f"- {f}" for f in self.fixes] + [""]
        elif self.fixes:
            out += ["**Fixes applied to a copy** (your file is not changed):"] + [f"- {f}" for f in self.fixes] + [""]
        if self.warnings:
            out += ["**Warnings:**"] + [f"- {w}" for w in self.warnings] + [""]
        return "\n".join(out).rstrip() + "\n"

    def as_dict(self) -> dict:
        return {"name": self.name, "disinfectant": self.disinfectant, "sha1": self.sha1, "sniff": self.sniff,
                "not_a_model": self.not_a_model, "ok": self.ok, "summary": self.summary(),
                "rows": [{"item": i, "value": v, "source": s} for i, v, s in self.rows], "blocking": self.blocking,
                "fixes": self.fixes, "warnings": self.warnings, "run_cost": self.run_cost, "prepared": self.prepared,
                "same_file": self.same_file, "flags": self.flags}

    def bundled_dict(self) -> dict:
        """What outputs/intake/readiness_bundled.json keeps of a report: everything that does not vary from run to run
        or machine to machine (the measured run time, the minutes estimated from it and the worker count are left out;
        they go to the checks' run log), and no file name."""
        cost = {k: v for k, v in self.run_cost.items() if k not in VOLATILE_COST}
        flags = {k: v for k, v in self.flags.items() if k != "run_cost"}
        return {"ok": self.ok, "blocking": self.blocking, "fixes": self.fixes, "warnings": self.warnings, "run_cost": cost,
                "flags": flags, "rows": [{"item": i, "value": v, "source": s} for i, v, s in self.rows
                                         if i not in (RUN_COST_ROW, "File", "Text encoding")]}

    def to_json(self) -> str:
        return json.dumps(self.as_dict(), indent=1)


def stem(name: str) -> str:
    """The name without its folder, any .bak and one .inp or .txt (or other) extension, made safe for a file name."""
    b = os.path.basename(name or "model")
    while b.lower().endswith(".bak"):
        b = b[:-4]
    root, ext = os.path.splitext(b)
    b = root if ext else b
    b = re.sub(r"[^A-Za-z0-9._-]+", "_", b).strip("._") or "model"
    return b[:80]


def _sha1(b: bytes) -> str:
    return hashlib.sha1(b).hexdigest()


def run_cost(n_junctions: int, disinfectant: str, uniform_c: bool = False, seconds_per_run: float | None = None,
             workers: int | None = None) -> dict:
    """The simulator grid for one dose (task 15b): its EPANET runs, its cache file (ln C as float32 for every run, hour
    of the last day and junction: 64.8 KB per junction for the 675-run free-chlorine grid, 138.2 KB for chloramine's
    1440 runs), the temporary EPANET files of one run per worker while it builds, and, given one timed 7-day run, the
    minutes it takes (runs times that run, over the measured effective speedup EFFECTIVE_SPEEDUP, not the worker
    count).  'large': more junctions than LARGE_JUNCTIONS, which needs an explicit confirm."""
    from .simgp import GRIDS
    grid = grid_for(disinfectant, uniform_c)
    runs = 1
    for axis in GRIDS[grid]:
        runs *= len(axis)
    mb = runs * 24 * n_junctions * FLOAT32_BYTES / 1e6
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    out = {"grid": grid, "runs_per_dose": runs, "grid_file_mb": round(mb, 1), "junctions": n_junctions,
           "large": n_junctions > LARGE_JUNCTIONS, "workers": workers,
           "temp_files_mb": round(TEMP_MB_PER_JUNCTION * n_junctions * workers, 1)}
    if seconds_per_run is not None:
        est = estimate_seconds(float(seconds_per_run), runs)
        out["seconds_per_run"] = round(float(seconds_per_run), 3)
        out["estimate_s"] = est                   # unrounded: what run_cost_gate scores
        out["minutes"] = round(est / 60, 1)
    return out


def estimate_seconds(seconds_per_run: float, runs: int) -> float:
    """The grid build's time estimate (task 15b, item 10, as the plan states it): one timed 7-day run times the number
    of runs, over the assumed effective parallel speedup EFFECTIVE_SPEEDUP.  On a small model, where starting the
    workers takes most of the time, the build takes longer than this, and on a large one less
    (outputs/intake/run_cost_timing.json: Net3's builds took 1.4 to 2.1 times the estimate, ky10's 0.64)."""
    return seconds_per_run * runs / EFFECTIVE_SPEEDUP


def prepared_name(digest: str, name: str, keep_name: bool) -> str:
    """The prepared copy's file name: <digest>_<name> when keep_name (a clean upload, so the app's upload name and
    grid cache name stay as they were before task 15) and the name ends in .inp and fits in MAX_NAME_BYTES of UTF-8;
    otherwise <digest>_<stem>.inp, whose stem is at most 80 safe characters.  A file name has at most 255 bytes, and the
    grid cache's name adds its own words to this one, so a long upload name is shortened here, never refused."""
    base = os.path.basename(name or "")
    if keep_name and base.lower().endswith(".inp") and len(f"{digest}_{base}".encode("utf-8")) <= MAX_NAME_BYTES:
        return f"{digest}_{base}"
    return f"{digest}_{stem(name)}.inp"


def prepare_inp(src, out_dir: str, disinfectant: str = "free_chlorine", name: str | None = None,
                parse: bool = True) -> tuple[str | None, Readiness]:
    """(inp_path, Readiness) for a model file given as a path or as bytes.

    A clean UTF-8 file given as a path comes back as that path, untouched.  Clean bytes are written to out_dir as
    <sha1 of the bytes>[:8]_<name> (the app's upload name; a name not ending in .inp, or too long, gets <stem>.inp:
    prepared_name).  A file that needed a fix (a byte-order mark, UTF-16, a Windows code page, a model inside a zip)
    is written to out_dir as <sha1 of the new text>[:8]_<stem>.inp; the original is never changed.  A file that is not
    a model, or has a blocking problem, gives (None, Readiness) and leaves nothing in out_dir: a new copy is written
    under a temporary name, read by WNTR from there (after the ID scan, never before), and moved to its name only when
    WNTR read it, the model has an entry point and EPANET ran its first hour; on a blocking problem the temporary file
    is removed.  parse=True does that load and test run, for the report's counts and checks; a failure in either is
    a blocking problem.  With parse=False the copy is moved into place at once, unread.

    Never raises for what is in the file: a failure anywhere after the file is read (a zip that breaks on reading, a
    copy that cannot be saved, an error no test foresaw) comes back as a blocking problem in plain words, with
    nothing left in out_dir.  A wrong disinfectant, or a path that does not exist, still raises."""
    if disinfectant not in DISINFECTANT_WORDS:
        raise ValueError(f"disinfectant must be one of {sorted(DISINFECTANT_WORDS)}")
    if isinstance(src, (bytes, bytearray, memoryview)):
        data, path = bytes(src), None
        name = name or "upload.inp"
    else:
        path = os.fspath(src)
        name = name or os.path.basename(path)
        with open(path, "rb") as fh:
            data = fh.read()
    rep = Readiness(name=os.path.basename(name), disinfectant=disinfectant, sha1=_sha1(data))
    try:
        return _prepare(data, path, name, out_dir, disinfectant, parse, rep)
    except Exception as e:  # noqa: BLE001  (a crafted or broken file must reach the user as a plain problem)
        rep.prepared, rep.same_file = None, False
        rep.blocking.append(f"This file could not be prepared ({type(e).__name__}); nothing was kept. Send the "
                            f"model's EPANET export (.inp) again, on its own.")
        return None, rep


def _stale_names(name: str, contents: list[bytes]) -> set[str]:
    """Every content name a prepared copy of this file can have had, under task 15's scheme and this one: the digest
    of the bytes as received, as decoded, or after the reaction fixes, with the upload name kept or the safe stem.
    The digest makes each name this file's own."""
    out = set()
    for b in contents:
        d = _sha1(b)[:8]
        out.update(prepared_name(d, name, keep) for keep in (True, False))
    return out


def _remove_stale(out_dir: str, names: set[str], keep: str | None) -> None:
    """A blocked model leaves no prepared copy behind (task 15b): a content-named copy that an earlier version (or an
    earlier, less strict check) accepted is removed from out_dir.  The file given (keep) is never removed."""
    for n in names:
        p = os.path.join(out_dir, n)
        if keep is not None and os.path.abspath(p) == os.path.abspath(keep):
            continue
        if os.path.isfile(p) and not os.path.islink(p):
            with contextlib.suppress(OSError):
                os.remove(p)


def _write_temp(out_dir: str, b: bytes) -> str:
    """b in a new temporary file in out_dir (removed again if the write fails)."""
    fd, tmp = tempfile.mkstemp(prefix=".prepare_", suffix=".inp", dir=out_dir)
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(b)
    except BaseException:
        with contextlib.suppress(OSError):
            os.remove(tmp)
        raise
    return tmp


def _decode_fix(rep: Readiness, data: bytes, how: str) -> None:
    """The decoding fix (and, for a file that mixes two encodings, a warning naming the lines)."""
    if how == "utf-8":
        return
    if how.startswith(MIXED_PREFIX):
        _, n, lines, n_lines, pairs = _mixed_utf8(data)
        chance = (f", {2 * pairs if pairs > 1 else 'two'} of them in {pairs if pairs > 1 else 'a'} pair"
                  f"{'s' if pairs > 1 else ''} that {'are' if pairs > 1 else 'is'} valid UTF-8 only by chance (an "
                  f"accented capital before a quote, a dash or a degree sign, as in 'CITÉ’')" if pairs else "")
        rep.fixes.append(f"Read the text as UTF-8, with {n} stray byte{'s' if n != 1 else ''} read as cp1252 (Windows "
                         f"Western European){chance}, and saved it as UTF-8 (WNTR reads UTF-8 only).")
        rep.warnings.append(f"The file mixes two text encodings: it is UTF-8 except for {n} byte{'s' if n != 1 else ''} "
                            f"(on line{'s' if n_lines != 1 else ''} {', '.join(map(str, lines))}"
                            f"{' and more' if n_lines > len(lines) else ''}), which {'were' if n != 1 else 'was'} "
                            f"read as cp1252 (Windows Western European), as a pasted degree sign or accented letter "
                            f"usually is. If a name or an ID there looks wrong, correct it in the modelling package and "
                            f"export again.")
        return
    what = {"utf-8 with a byte-order mark": "Removed the byte-order mark at the start (WNTR stops at line 1 on it)"}
    rep.fixes.append(what.get(how, f"Decoded the text from {how} and saved it as UTF-8 (WNTR reads UTF-8 only)") + ".")


def _ids_list(ids: list[str], k: int = 6) -> str:
    return ", ".join(ids[:k]) + (f" and {len(ids) - k} more" if len(ids) > k else "")


def _reaction_fixes(rep: Readiness, info: dict) -> None:
    """The reaction settings fix_reactions_text changed, as fixes (task 15b, items 1 and 2), kept for the rows."""
    rep.flags["reactions"] = {"orders": info["orders"], "orders_fixed": info["orders_fixed"],
                              "bulk_lines": len(info["bulk_lines"]), "tank_lines": len(info["tank_lines"]),
                              "wall_lines": info["wall_lines"], "limiting_potential": info["limiting_potential"],
                              "limiting_potential_fixed": info["limiting_potential_fixed"],
                              "roughness_correlation": info["roughness_correlation"]}
    of = info["orders_fixed"]
    if of:
        rep.fixes.append("Set the reaction order to 1, as the decay model assumes (first order), where the file had "
                         + " and ".join(f"{k} {v}" for k, v in of.items()) + ". An order of 0 is another law (a fixed "
                         "amount lost per day), so the decay rates would not mean what they say, and the bulk rate can "
                         "have no effect at all.")
    b, t = info["bulk_lines"], info["tank_lines"]
    if b or t:
        parts = ([f"the BULK lines of {len(b)} pipe{'s' if len(b) != 1 else ''} ({_ids_list(b)})"] if b else []) + \
                ([f"the TANK lines of {len(t)} tank{'s' if len(t) != 1 else ''} ({_ids_list(t)})"] if t else [])
        rep.fixes.append("Cleared " + " and ".join(parts) + ": a per-pipe or per-tank bulk rate overrides the one the "
                         "model calibrates, so those pipes and tanks would ignore it. The calibrated bulk rate now "
                         "applies to every pipe and tank.")
    if info["limiting_potential_fixed"] is not None:
        rep.fixes.append(f"Set LIMITING POTENTIAL to 0, where the file had {info['limiting_potential_fixed']}: with a "
                         f"limiting potential the bulk decay is not first order in the chlorine.")


def _reaction_rows(rep: Readiness) -> None:
    r = rep.flags.get("reactions")
    if r is None:
        return
    if r["orders_fixed"]:
        had = ", ".join(f"{k} {v}" for k, v in r["orders_fixed"].items())
        rep.add("Reaction orders", f"the file has {had}; the copy uses 1 (first order) for bulk, wall and tank, as the "
                                   f"decay model assumes", DERIVED)
    else:
        rep.add("Reaction orders", "1 (first order) for bulk, wall and tank, as the decay model assumes"
                + ("" if len(r["orders"]) == 3 else " (EPANET's default where the file sets none)"), FROM_FILE)
    bits = []
    if r["bulk_lines"] or r["tank_lines"]:
        bits.append(f"{r['bulk_lines']} per-pipe BULK and {r['tank_lines']} per-tank TANK lines cleared in the copy, "
                    f"so the calibrated bulk rate applies everywhere")
    if r["limiting_potential_fixed"] is not None:
        bits.append(f"LIMITING POTENTIAL {r['limiting_potential_fixed']} set to 0")
    if r["wall_lines"] or r["roughness_correlation"]:
        bits.append((f"{r['wall_lines']} per-pipe WALL lines" if r["wall_lines"] else "")
                    + (" and " if r["wall_lines"] and r["roughness_correlation"] else "")
                    + (f"a ROUGHNESS CORRELATION of {r['roughness_correlation']:g}" if r["roughness_correlation"] else "")
                    + " ignored: the model sets every pipe's wall rate from its C")
    rep.add("Per-pipe and per-tank rates", "; ".join(bits) if bits else "none: one bulk rate for every pipe and tank, "
                                                                        "as the decay model assumes",
            DERIVED if bits else FROM_FILE)


def _prepare(data: bytes, path: str | None, name: str, out_dir: str, disinfectant: str, parse: bool,
             rep: Readiness) -> tuple[str | None, Readiness]:
    s = sniff(data, name)
    rep.sniff = s.as_dict()
    if s.kind == "epanet_inp_zip":
        if len(s.members) != 1:
            rep.blocking.append(f"The zip holds {len(s.members)} EPANET models ({', '.join(s.members[:5])}): send one "
                                f"model at a time.")
            return None, rep
        try:
            data = read_zip_member(data, s.members[0])
        except ZipRefused as e:
            rep.blocking.append(f"The zip is refused: {e}.")
            return None, rep
        rep.fixes.append(f"Read {s.members[0]} from the zip in memory (the zip was not unpacked); the model is saved as "
                         f"a prepared copy.")
        name, path = s.members[0], None
        inner = sniff(data, name)
        if inner.kind != "epanet_inp":
            rep.not_a_model = True
            rep.sniff = inner.as_dict()
            return None, rep
        s = inner
    elif not s.is_model:
        rep.not_a_model = True
        return None, rep
    rep.add("File", "an EPANET model (.inp), identified by its content"
            + ("; the name marks it as a backup copy" if s.backup else ""), FROM_FILE)
    text, how = decode_text(data)
    rep.add("Text encoding", how, FROM_FILE)
    _decode_fix(rep, data, how)
    ids = scan_inp_ids(text)
    rep.blocking += ids.blocking
    rep.warnings += ids.warnings
    fixed, rinfo = (text, None) if rep.blocking else fix_reactions_text(text)
    stale = _stale_names(name, [data, text.encode("utf-8"), fixed.encode("utf-8")])
    if rep.blocking:
        _remove_stale(out_dir, stale, path)
        return None, rep
    _reaction_fixes(rep, rinfo)
    changed = fixed != text
    # write only when needed, under a temporary name until WNTR has read it, an entry point is found, the physics
    # checks and the test runs passed; a temporary file left by any failure (the save itself included) is removed in
    # the finally, and so is a content-named copy of this file that an earlier version accepted (task 15b)
    tmp, tmp2, accepted, out_path = None, None, False, None
    try:
        if path is not None and how == "utf-8" and not changed:
            read_path = out_path = path
            rep.same_file = True
        else:
            new = data if how == "utf-8" and not changed else fixed.encode("utf-8")
            fname = prepared_name(_sha1(new)[:8], name, keep_name=how == "utf-8" and not rep.fixes)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, fname)
            read_path = tmp = _write_temp(out_dir, new)
        rep.prepared = os.path.basename(out_path)
        replacement = None
        if parse:
            accepted, replacement = _read_and_report(read_path, out_path, rep, ids, disinfectant, fixed)
        else:
            accepted = True
        if accepted and replacement is not None:      # a Darcy-Weisbach model, converted: the copy holds the new text
            new2 = replacement.encode("utf-8")
            out_path = os.path.join(out_dir, prepared_name(_sha1(new2)[:8], name, keep_name=False))
            os.makedirs(out_dir, exist_ok=True)
            tmp2 = _write_temp(out_dir, new2)
            os.replace(tmp2, out_path)
            tmp2 = None
            rep.prepared, rep.same_file = os.path.basename(out_path), False
        elif accepted and tmp is not None:
            os.replace(tmp, out_path)
            tmp = None
    except OSError as e:
        accepted = False
        rep.blocking.append(f"The prepared copy could not be saved ({e.strerror or type(e).__name__}); nothing was "
                            f"kept.")
    finally:
        for t in (tmp, tmp2):
            if t is not None:
                with contextlib.suppress(OSError):
                    os.remove(t)
        if not accepted:
            rep.prepared, rep.same_file = None, False
            _remove_stale(out_dir, stale, path)
    return (out_path if accepted else None), rep


def _test_run(wn) -> str | None:
    """EPANET on the model's first hour, set up the way the app's first run is (simulate.load's time steps, water
    age), in a temporary folder that is removed afterwards.  None when it runs; otherwise the problem in plain words.
    WNTR reading a file does not mean EPANET runs it: a node that no pipe, pump or valve connects reads fine and then
    stops EPANET with its Error 200 (found by the final verification of task 15)."""
    import wntr
    from .simulate import QUALITY_STEP_S
    t = wn.options.time
    t.duration, t.hydraulic_timestep, t.report_timestep, t.quality_timestep = 3600, 3600, 3600, QUALITY_STEP_S
    wn.options.quality.parameter = "AGE"
    try:
        with tempfile.TemporaryDirectory(prefix="rm_test_run_") as d, warnings.catch_warnings():
            warnings.simplefilter("ignore")
            prefix = os.path.join(d, "check")
            try:
                wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=prefix)
            except Exception as e:  # noqa: BLE001
                msg = str(e).replace(prefix, "the test copy").replace(d, "").replace(" %s", "").replace("%s", "")
                linked = set()
                for _, link in wn.links():
                    linked.update((link.start_node_name, link.end_node_name))
                lone = [n for n in wn.node_name_list if n not in linked]
                hint = ""
                if lone:
                    hint = (f" {len(lone)} node{'s' if len(lone) > 1 else ''} in the file "
                            f"{'are' if len(lone) > 1 else 'is'} not connected to any pipe, pump or valve "
                            f"({', '.join(lone[:5])}{' and more' if len(lone) > 5 else ''}), and EPANET stops on such "
                            f"a node.")
                return (f"EPANET could not run the first hour of this model ({type(e).__name__}: {msg[:200]})." + hint
                        + " Open the model in EPANET or the modelling package, fix what it reports, and export it "
                          "again.")
    except OSError as e:          # the temporary folder itself
        return f"The EPANET test run could not start ({e.strerror or type(e).__name__})."
    return None


ELEVATION_EVIDENCE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs", "intake",
                                  "elevation_prv_ky10.json")


def _elevation_evidence() -> str:
    """One sentence from the committed outputs/intake/elevation_prv_ky10.json, or '' when it is not there."""
    try:
        with open(ELEVATION_EVIDENCE) as fh:
            d = json.load(fh)
        return (f" In a test on the public model ky10 (its {d['pressure_valves']} PRVs), setting every junction to one "
                f"elevation moved the water age by more than {d['threshold_h']:g} h at {d['junctions_moved']} of "
                f"{d['junctions']} junctions, by up to {d['max_change_h']:.1f} h (outputs/intake/elevation_prv_ky10.json).")
    except (OSError, KeyError, ValueError, TypeError):
        return ""


def _clock_minute(hours: float) -> int:
    """A time in hours as the minute of the day, rounded once (5:59:59 is minute 360; 23:59:50 and 0:00:20 are 0)."""
    return int(round(hours * 60)) % 1440


def _clock(hours: float) -> str:
    """HH:MM of a time in hours, rounded to the minute once (so 5:59:59 is 06:00, never 05:60)."""
    m = _clock_minute(hours)
    return f"{m // 60:02d}:{m % 60:02d}"


def _clock_report(rep: Readiness, clock_h: float) -> None:
    """The start clock row and, for a start other than midnight, the model-hours warning (task 15b, item 4).  Both
    quote the start rounded to the minute once, and the warning goes by that minute (the second review: a start within
    30 s of midnight is 00:00 and gets no warning, and the offset quoted is the rounded one)."""
    m = _clock_minute(clock_h)
    off = f"{m / 60:g} h"
    rep.add("Start clock time", f"{off} ({_clock(clock_h)}): the model's hour 0 is this time of day; recorded, the file "
                                f"is not changed", FROM_FILE)
    if m:
        rep.warnings.append(f"The model starts at {_clock(clock_h)} (START CLOCKTIME), so its hour 0 is {_clock(clock_h)}, "
                            f"not midnight. Until clock time is converted (planned), every hour in the app and the "
                            f"pilot is a model hour: the map's hours, the worst hour, the route's hours, and the hours "
                            f"of your samples, which are read as model hours too. Add {off} to a model hour for "
                            f"the time of day, and take {off} off a sample's time of day before entering it.")


def _physics_report(rep: Readiness, wn, ctl: dict, have_xy: set, src_nodes: list[str]) -> set[str]:
    """The task-15b checks that read the model WNTR loaded (no EPANET run): headloss formula, start clock and pattern
    lengths, pressure-control valves and elevations, coordinates, sources and entry points, pumps and controls,
    demand, connectivity and export connectors.  Rows, warnings and blocking problems go to rep; what each found goes
    to rep.flags.  Returns the export's connector pipes when the Export connectors row recognises them (the pipes
    tagged AddedForINPConversion, and the pipes of 1 ft and 99 in once the row is shown), which the uniform-C test
    leaves out; else an empty set."""
    import networkx as nx
    from wntr.network import LinkStatus
    f = rep.flags
    _reaction_rows(rep)
    # 3. headloss formula (the conversion itself runs later, after the test run)
    hl = wn.options.hydraulic.headloss
    f["headloss"] = hl
    if hl == "H-W":
        rep.add("Headloss formula", "H-W (Hazen-Williams), as the decay model assumes: the wall decay and the "
                                    "roughness axis read each pipe's roughness as a C", FROM_FILE)
    elif hl == "D-W":
        rep.add("Headloss formula", "D-W (Darcy-Weisbach) in your file; converted to Hazen-Williams in the copy, each "
                                    "pipe's C matched to its headloss at its reference velocity (next rows)", DERIVED)
    else:
        rep.add("Headloss formula", f"{hl} (Chezy-Manning) in your file", FROM_FILE)
        rep.blocking.append("The file uses the Chezy-Manning headloss formula. The wall decay reads each pipe's roughness "
                            "as a Hazen-Williams C, and a Manning n has no such reading, so the decay rates would not "
                            "mean what they say. Export the model with Hazen-Williams, or ask your modeller.")
    # 4. start clock time and pattern lengths
    t = wn.options.time
    clock_h = float(t.start_clocktime or 0) / 3600
    f["start_clock_h"] = clock_h
    _clock_report(rep, clock_h)
    step = float(t.pattern_timestep or 3600)
    weekly, multi = [], []
    for pn in wn.pattern_name_list:
        dur_h = len(wn.get_pattern(pn).multipliers) * step / 3600
        if abs(dur_h - 168) < 1e-9:
            weekly.append(pn)
        elif dur_h > 24 + 1e-9:
            multi.append([pn, dur_h])
    f["weekly_patterns"], f["multi_day_patterns"] = weekly, multi
    if weekly:
        rep.warnings.append(f"{len(weekly)} weekly pattern{'s' if len(weekly) != 1 else ''} ({_ids_list(weekly)}: a "
                            f"week long). Every run is 7 days and only the last day is scored, so the scored day is "
                            f"the pattern's day 7, six days after the weekday it starts on. The file does not record "
                            f"that weekday: tell us which day of the week the pattern starts on, so the scored day can "
                            f"be named.")
    if multi:
        rep.warnings.append(f"{len(multi)} pattern{'s' if len(multi) != 1 else ''} longer than a day and not a week ("
                            + ", ".join(f"{pn}: {d:g} h" for pn, d in multi[:6]) + "): EPANET repeats each one over its "
                            "own length, so the scored day (the last of the 7-day run) is not a repeat of one daily "
                            "cycle.")
    # 5. pressure-control valves and elevations
    junctions = wn.junction_name_list
    kinds: dict[str, list[str]] = {}
    for vn in wn.valve_name_list:
        kinds.setdefault(str(wn.get_link(vn).valve_type), []).append(vn)
    pvalves = [vn for k in PRESSURE_VALVES for vn in kinds.get(k, [])]
    elev = {j: float(wn.get_node(j).elevation) for j in junctions}
    zeros = [j for j in junctions if elev[j] == 0.0]
    constant = len(junctions) > 1 and len(set(elev.values())) == 1
    pnodes = [j for j in junctions if j in ctl["nodes"]]
    need = set(pnodes)
    for vn in pvalves:
        link = wn.get_link(vn)
        need.update((link.start_node_name, link.end_node_name))
    need_j = [j for j in junctions if j in need]
    emitters = [j for j in junctions if (wn.get_node(j).emitter_coefficient or 0) > 0]
    pressure_driven = str(wn.options.hydraulic.demand_model).upper() in ("PDA", "PDD")
    f["pressure_valves"] = {k: len(kinds.get(k, [])) for k in PRESSURE_VALVES if kinds.get(k)}
    f["zero_elevation_junctions"] = zeros
    f["zero_elevation_flagged"] = bool(pvalves and zeros and not constant)
    f["constant_elevation"] = constant
    f["elevations_needed_at"] = len(need_j)
    if pvalves:
        what = ", ".join(f"{len(kinds[k])} {k}" for k in PRESSURE_VALVES if kinds.get(k))
        at_zero = [j for j in need_j if j in set(zeros)]
        rep.add("Elevations", f"{len(pvalves)} pressure-control valve{'s' if len(pvalves) != 1 else ''} ({what}): the "
                              f"elevations change the flows here; real elevations are needed at the valves"
                              + (", downstream of them and in pressure-based controls" if pnodes else
                                 " and downstream of them")
                              + f"; {len(need_j)} junction{'s' if len(need_j) != 1 else ''} "
                              + ("sit on either side of a valve or in such a control" if pnodes else
                                 "sit on either side of a valve")
                              + (f", {len(at_zero)} of them at elevation 0" if at_zero else
                                 ", none of them at elevation 0"), FROM_FILE)
        if constant:
            rep.blocking.append(
                f"Elevations needed at these valves: {_ids_list(pvalves)}. Every junction has the same elevation "
                f"({next(iter(elev.values())):g} m), and a pressure-reducing, -sustaining or -breaker valve works on the "
                f"pressure, the head minus the elevation, so with one made-up elevation the valves set the wrong flows "
                f"and the water age is wrong." + _elevation_evidence() + " Export the model with its real elevations, "
                f"at least at the valves, the junctions downstream of them and the junctions named in pressure-based "
                f"controls.")
        elif zeros:
            rep.warnings.append(
                f"{len(zeros)} junction{'s' if len(zeros) != 1 else ''} at elevation 0 in a model with pressure-control "
                f"valves: {_ids_list(zeros, 10)}. If a 0 means 'not surveyed', the pressure there is wrong, and near a "
                f"valve or a pressure-based control so are the flows." + _elevation_evidence() + " Check these "
                f"elevations in the modelling package.")
    elif pnodes or emitters or pressure_driven:
        users = ([f"{len(pnodes)} junctions read by pressure-based controls"] if pnodes else []) + \
                ([f"{len(emitters)} emitters"] if emitters else []) + (["pressure-driven demands"] if pressure_driven
                                                                       else [])
        rep.add("Elevations", "no pressure-control valves; " + " and ".join(users) + " use the pressure, so the "
                              "elevations there change the flows; read from your file", FROM_FILE)
        if constant:
            rep.warnings.append("Every junction has the same elevation, and the model uses the pressure ("
                                + " and ".join(users) + "): the flows there may be wrong. Export the model with its "
                                "real elevations.")
    else:
        rep.add("Elevations", "no pressure-control valves, pressure-based controls or emitters: in this demand-driven "
                              "run the junction elevations change the pressures, not the flows", FROM_FILE)
    # 6. coordinates
    cs = coordinate_scale(wn, have_xy)
    f["coordinates"] = cs
    later = "tap mapping by coordinates, when it comes (planned), will not be offered for this file; use junction IDs"
    if cs["status"] == "no coordinates":
        rep.add("Coordinate scale", f"no coordinates: the maps will show every node at one point, and {later}", MISSING)
    elif cs["status"] == "too few pipes drawn":
        rep.add("Coordinate scale", f"too few pipes drawn to tell: {later}", DERIVED)
    elif cs["status"] == "schematic":
        rep.add("Coordinate scale", f"schematic: the drawn pipe lengths are not in proportion to the lengths in the "
                                    f"file (95th over 5th percentile of their ratio {cs['spread']:g}, above "
                                    f"{SCHEMATIC_SPREAD:g}), so {later}", DERIVED)
    else:
        unit = {"feet": "in feet", "metres": "in metres"}.get(cs["units"], "in units that are not feet or metres")
        rep.add("Coordinate scale", f"to scale, {unit} (median drawn length {cs['ratio_median']:g} per metre of "
                                    f"declared length; 95th over 5th percentile {cs['spread']:g})", DERIVED)
    # 7. sources in the file and entry points
    dosed = set(src_nodes)
    srcs = [(s.node_name, str(s.source_type)) for _, s in wn.sources()]
    over = [n for n, _ in srcs if n in dosed]
    boost = [n for n, _ in srcs if n not in dosed]
    f["sources_overridden"], f["booster_sources"], f["entry_points"] = over, boost, list(src_nodes)
    if not srcs:
        rep.add("Sources in the file", "none: the dose you set enters at every entry point", FROM_FILE)
    else:
        rep.add("Sources in the file", "; ".join(([f"at {_ids_list(over)}, an entry point: replaced by the dose you "
                                                   f"set"] if over else [])
                                                 + ([f"at {_ids_list(boost)}, not an entry point (a booster): kept in "
                                                     f"the file, but at the tool's dose scale it adds almost nothing"]
                                                    if boost else [])), FROM_FILE)
    if over:
        rep.warnings.append(f"The file's own source at {_ids_list(over)} is overridden: the tool doses every entry point "
                            f"with the dose you set.")
    if boost:
        rep.warnings.append(f"{len(boost)} source{'s' if len(boost) != 1 else ''} at node{'s' if len(boost) != 1 else ''} "
                            f"the tool does not dose ({_ids_list(boost)}): boosters, kept as the file has them. The tool "
                            f"writes its own dose to EPANET in the model's units, a thousand times the file's mg/L, and "
                            f"leaves these as they are, so they stay in the run at about a thousandth of their strength "
                            f"and add almost nothing. How to treat boosters is an open decision; until then the map "
                            f"has almost none of their chlorine.")
    if len(src_nodes) > 1:
        rep.warnings.append(f"{len(src_nodes)} entry points ({_ids_list(list(src_nodes))}): every one is dosed at the "
                            f"same value, the dose you set. If they carry different residuals, the map is off near "
                            f"them; a residual per inlet is planned.")
    # 8. structure: pumps and controls, demand and patterns, connected parts, export connectors
    loose = []
    for pn in wn.pump_name_list:
        p = wn.get_link(pn)
        if getattr(p, "speed_pattern_name", None) or pn in ctl["links"]:
            continue
        loose.append((pn, "closed" if p.initial_status == LinkStatus.Closed else "open"))
    f["pumps_without_controls"] = [pn for pn, _ in loose]
    f["controls"], f["rules"] = ctl["n_controls"], ctl["n_rules"]
    rep.add("Pumps and controls", f"{wn.num_pumps} pump{'s' if wn.num_pumps != 1 else ''}; {ctl['n_controls']} "
                                  f"control{'s' if ctl['n_controls'] != 1 else ''} and {ctl['n_rules']} "
                                  f"rule{'s' if ctl['n_rules'] != 1 else ''}", FROM_FILE)
    if loose:
        opened = [pn for pn, st in loose if st == "open"]
        closed = [pn for pn, st in loose if st == "closed"]
        desc = "; ".join(([f"open all the time: {_ids_list(opened)}"] if opened else [])
                         + ([f"closed all the time: {_ids_list(closed)}"] if closed else []))
        rep.warnings.append(f"Pumps but no controls: {len(loose)} of {wn.num_pumps} pump{'s' if wn.num_pumps != 1 else ''} "
                            f"{'has' if len(loose) == 1 else 'have'} no control, rule or speed pattern ({desc}), so "
                            f"{'it stays' if len(loose) == 1 else 'each stays'} as the file starts it for all 7 days. If "
                            f"the real pumps cycle, the water age and the tanks' chlorine are wrong: export the model "
                            f"with its controls.")
    total = sum(ts.base_value for j in junctions for ts in wn.get_node(j).demand_timeseries_list if ts.base_value > 0)
    f["patterns"], f["zero_total_demand"] = len(wn.pattern_name_list), total == 0
    rep.add("Demand", f"{len(wn.pattern_name_list)} time pattern{'s' if len(wn.pattern_name_list) != 1 else ''}"
                      + ("; no junction takes water" if total == 0 else ""), FROM_FILE)
    if not wn.pattern_name_list:
        rep.warnings.append("No time patterns: every demand is the same all day, so the daily low point of the residual "
                            "is not represented.")
    if total == 0:
        rep.warnings.append("Zero total demand: no junction takes water, so water moves only as the tanks fill and "
                            "drain, and its age is not the network's.")
    g = nx.Graph()
    g.add_nodes_from(wn.node_name_list)
    g.add_edges_from((link.start_node_name, link.end_node_name) for _, link in wn.links())
    parts = sorted((sorted(c) for c in nx.connected_components(g)), key=lambda c: (-len(c), c[0]))
    f["parts"] = [len(c) for c in parts]
    if len(parts) > 1:
        dry = [c for c in parts if not set(c) & (dosed | set(wn.tank_name_list) | set(wn.reservoir_name_list))]
        rep.add("Connected parts", f"{len(parts)} parts that no pipe, pump or valve joins (sizes "
                                   f"{', '.join(str(len(c)) for c in parts[:6])})", DERIVED)
        rep.warnings.append(f"{len(parts)} parts that no pipe, pump or valve joins (sizes "
                            f"{', '.join(str(len(c)) for c in parts[:6])}{' and more' if len(parts) > 6 else ''}): "
                            f"each part needs its own entry point"
                            + (f", and {len(dry)} part{'s have' if len(dry) != 1 else ' has'} none (for example "
                               f"{_ids_list(dry[0], 3)}), so no chlorine reaches it" if dry else "")
                            + ". If the zone was cut from a larger model, check the cut.")
    else:
        rep.add("Connected parts", "one connected network", DERIVED)
    conn = [pn for pn in wn.pipe_name_list
            if abs(wn.get_link(pn).length - CONNECTOR_LENGTH_M) < 0.01 and abs(wn.get_link(pn).diameter
                                                                                 - CONNECTOR_DIAM_M) < 0.01]
    f["connector_pipes"], f["connector_tags"] = len(conn), len(ctl["connector_tags"])
    if ctl["connector_tags"] or len(conn) >= CONNECTOR_MIN_PIPES:   # one pipe of that size alone says nothing (Net6)
        rep.add("Export connectors", f"{len(conn)} pipe{'s' if len(conn) != 1 else ''} of 1 ft and 99 in (the size a "
                                     f"WaterGEMS EPANET export uses for its connector pipes) and "
                                     f"{len(ctl['connector_tags'])} element{'s' if len(ctl['connector_tags']) != 1 else ''} "
                                     f"tagged AddedForINPConversion (the tag it gives them); kept as they are (for "
                                     f"information)", FROM_FILE)
        return set(conn) | (set(ctl["connector_tags"]) & set(wn.pipe_name_list))
    return set()


def _dw_report(rep: Readiness, conv: dict) -> None:
    """Rows, fix or block for a Darcy-Weisbach file and its gate (task 15b, item 3)."""
    from .simulate import roughness_factor
    keep = ("pipes", "c_min", "c_max", "age_change_median_h", "age_change_p95_h", "passes")
    rep.flags["headloss_conversion"] = {k: (round(conv[k], 4) if isinstance(conv[k], float) else conv[k]) for k in keep}
    med, p95 = conv["age_change_median_h"], conv["age_change_p95_h"]
    rep.add("Headloss conversion", f"water age moved by a median of {med:.2f} h and a 95th percentile of {p95:.2f} h "
                                   f"over the scored day's junction-hours between your file and the copy; the copy is "
                                   f"kept up to {DW_GATE_MEDIAN_H:g} h and {DW_GATE_P95_H:g} h", DERIVED)
    if conv["passes"]:
        rep.fixes.append(f"Converted the Darcy-Weisbach headloss formula to Hazen-Williams: each of the {conv['pipes']} "
                         f"pipes gets the C ({conv['c_min']:.0f} to {conv['c_max']:.0f}) that gives its Darcy-Weisbach "
                         f"headloss at its reference velocity (the median over the scored day, at least "
                         f"{V_REF_FLOOR_MS:g} m/s). The wall decay reads each pipe's roughness as a C, and a "
                         f"Darcy-Weisbach roughness (a fraction of a millimetre) read as a C near 0 would multiply every "
                         f"wall rate by about {roughness_factor(0.0):.0f} at the old-pipe factor 1. The conversion is "
                         f"exact at that velocity only, which is why the water age is compared.")
    else:
        rep.blocking.append(f"The Darcy-Weisbach model could not be converted closely enough: water age moved by a "
                            f"median of {med:.2f} h (at most {DW_GATE_MEDIAN_H:g} h is kept) and a 95th percentile of "
                            f"{p95:.2f} h (at most {DW_GATE_P95_H:g} h) between the file and a Hazen-Williams copy. Export "
                            f"the model with Hazen-Williams, or ask your modeller.")


def _timed_member_run(path: str, gamma: float) -> float:
    """Seconds for one 7-day EPANET run of the file exactly as a grid member runs it (simulate_nominal_chlorine at
    kb 0.40, kw 0.70), in a temporary folder that is removed; it raises when EPANET cannot run it."""
    from .simulate import simulate_nominal_chlorine
    with tempfile.TemporaryDirectory(prefix="rm_test_run_") as d, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.perf_counter()
        simulate_nominal_chlorine(path, 0.40, 0.70, gamma, file_prefix=os.path.join(d, "member"))
        return time.perf_counter() - t0


def _read_and_report(read_path: str, out_path: str, rep: Readiness, ids: IdScan, disinfectant: str,
                     text: str) -> tuple[bool, str | None]:
    """WNTR's load of the prepared copy, the report rows that need it, the physics checks (task 15b), EPANET's test run
    of the first hour, the Darcy-Weisbach conversion and its gate, and one timed 7-day run for the run cost.
    (accepted, replacement): accepted is False (with a blocking problem) when WNTR cannot read it, the model has no
    entry point, a physics check blocks or a run fails; replacement is the converted text of a Darcy-Weisbach model,
    which _prepare saves instead of the copy read here."""
    import wntr
    from .simgp import GRIDS
    from .simulate import source_nodes
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wn = wntr.network.WaterNetworkModel(read_path)
    except Exception as e:  # noqa: BLE001
        msg = str(e).replace(read_path, os.path.basename(out_path))     # the file's name, never its folder
        rep.blocking.append(f"EPANET (through WNTR) could not read this model: {type(e).__name__}: {msg[:300]}")
        return False, None
    for item, n in (("Junctions", wn.num_junctions), ("Reservoirs", wn.num_reservoirs), ("Tanks", wn.num_tanks),
                    ("Pipes", wn.num_pipes), ("Pumps", wn.num_pumps), ("Valves", wn.num_valves)):
        rep.add(item, n, FROM_FILE)
    rep.add("Flow units", wn.options.hydraulic.inpfile_units, FROM_FILE)
    n_nodes = wn.num_nodes
    have_xy = ids.coordinate_ids & set(wn.node_name_list)
    if have_xy:
        rep.add("Node coordinates", f"{len(have_xy)} of {n_nodes} nodes", FROM_FILE)
    else:
        rep.add("Node coordinates", "none: the maps need them (EPANET or the modelling package writes them in "
                                    "[COORDINATES])", MISSING)
    src_nodes = source_nodes(wn)
    if not src_nodes:
        rep.blocking.append("No reservoir, inflow junction or tank found: the model needs a place where water enters.")
        return False, None
    how_in = ("reservoirs" if wn.reservoir_name_list else
              "junctions with a negative demand (inflows)" if set(src_nodes) <= set(wn.junction_name_list) else "tanks")
    rep.add("Where water enters", f"{how_in}: {', '.join(src_nodes[:6])}{' and more' if len(src_nodes) > 6 else ''}",
            DERIVED)
    rep.add("Disinfectant", f"{DISINFECTANT_WORDS[disinfectant]} (your setting; every entry point gets the same dose)",
            ASSUMED)
    connectors = _physics_report(rep, wn, control_scan(text), have_xy, src_nodes)
    if rep.blocking:
        return False, None
    headloss = wn.options.hydraulic.headloss
    # the uniform-C test leaves out the export's connector pipes (the second review: a WaterGEMS-style export of a
    # one-C network has a C 130 connector of 1 ft at each pump, which is not a second C value of the network)
    scored = [pn for pn in wn.pipe_name_list if pn not in connectors] or list(wn.pipe_name_list)
    roughness = [float(wn.get_link(pn).roughness) for pn in scored]
    left_out = len(wn.pipe_name_list) - len(scored)
    n_junctions = wn.num_junctions
    problem = _test_run(wn)            # it changes wn's time and quality options: nothing below reads them
    if problem:
        rep.blocking.append(problem)
        return False, None
    replacement, final_path = None, read_path
    with tempfile.TemporaryDirectory(prefix="rm_test_run_") as work:
        if headloss == "D-W":
            try:
                replacement, conv = dw_to_hw(read_path, text, work)
            except Exception as e:  # noqa: BLE001
                rep.blocking.append(f"The Darcy-Weisbach model could not be converted ({type(e).__name__}: "
                                    f"{str(e)[:200]}). Export the model with Hazen-Williams, or ask your modeller.")
                return False, None
            _dw_report(rep, conv)
            if not conv["passes"]:
                return False, None
            final_path = conv["converted_path"]
            keep = set(scored)
            converted = sorted({round(float(v), 6) for k, v in conv["c"].items() if k in keep})
            rep.flags["converted_c_values"] = len(converted)
        # the uniform test reads the file's own roughness, before any conversion: a Darcy-Weisbach file with one
        # roughness on every pipe gets many converted C values, which differ by diameter and velocity, not by age.
        # It compares them to 6 significant digits, not 6 decimals (the second review: a Darcy-Weisbach roughness in
        # metres is about 1e-6 for plastic pipe, so 6 decimals made 0.0015 mm and 0.0024 mm one value)
        own = sorted({float(f"{v:.6g}") for v in roughness})
        uniform = len(own) == 1
        rep.flags["c_values"], rep.flags["uniform_c"] = len(own), uniform
        but = (f" (the export's {left_out} connector pipe{'s' if left_out != 1 else ''} left out)" if left_out else "")
        if uniform:
            what = (f"one C value: every pipe has C {own[0]:g}{but}, so the factor would scale every wall rate alike "
                    f"and cannot be told from the wall rate itself" if headloss != "D-W" else
                    f"one roughness value: every pipe has the Darcy-Weisbach roughness {own[0] * 1000:g} mm in your "
                    f"file{but}, so "
                    f"the C values of the converted copy ({converted[0]:g} to {converted[-1]:g}) differ only by each "
                    f"pipe's diameter and velocity, not by its age, and the factor would read them as age")
            if grid_for(disinfectant, True) == "full_fixed_gamma":
                fixed_runs = math.prod(len(a) for a in GRIDS["full_fixed_gamma"])
                held = f"it is held at 0 (the free-chlorine grid with the factor held at 0 has {fixed_runs} runs)"
            else:
                held = ("the chloramine grid still varies it (a chloramine grid with the factor held at 0 is not built "
                        "yet), so its calibrated value means nothing for this file and the app shows it as not "
                        "identifiable")
            rep.add("Old-pipe factor", f"not identifiable from {what}; {held}", DERIVED)
        elif headloss == "D-W" and converted:
            rep.add("Old-pipe factor", f"learned from the {len(converted)} C values the conversion gives "
                                       f"({converted[0]:g} to {converted[-1]:g}), from your file's {len(own)} "
                                       f"roughness values{but}", DERIVED)
        elif own:
            rep.add("Old-pipe factor", f"learned from your file's {len(own)} C values ({own[0]:g} to "
                                       f"{own[-1]:g}){but}", FROM_FILE)
        try:
            secs = _timed_member_run(final_path, 0.0 if uniform else 1.0)
        except Exception as e:  # noqa: BLE001
            msg = str(e).replace(work, "")
            rep.blocking.append(f"EPANET could not run the 7-day simulation that every grid member is ({type(e).__name__}: "
                                f"{msg[:200]}). Open the model in EPANET or the modelling package, fix what it reports, "
                                f"and export it again.")
            return False, None
    rep.run_cost = run_cost(n_junctions, disinfectant, uniform, secs)
    c = rep.run_cost
    est = c["estimate_s"]
    rep.add(RUN_COST_ROW, f"{c['runs_per_dose']} EPANET runs per dose (the {c['grid']} grid), once per file: about "
                          + (f"{max(1, round(est)):d} s" if est < 60 else f"{c['minutes']:g} min")
                          + f" here, estimated as {c['runs_per_dose']} times one 7-day run ({c['seconds_per_run']:.2f} s) "
                          f"over an assumed parallel speedup of {EFFECTIVE_SPEEDUP:g}. The real speedup grows with the "
                          f"model's size, because starting the workers takes a fixed time: a small model's grid takes "
                          f"longer than this estimate, a large one's less. A grid file of about {c['grid_file_mb']:g} MB "
                          f"and up to "
                          f"{c['temp_files_mb']:g} MB of temporary files while it runs", DERIVED)
    if c["large"]:
        rep.warnings.append(f"{LARGE_WARNING_START} {n_junctions} junctions, more than {LARGE_JUNCTIONS}. The grid is "
                            f"{c['runs_per_dose']} runs and a file of about {c['grid_file_mb']:g} MB (the run cost row "
                            f"has the time), so building it needs an explicit confirm (a checkbox in the app's readiness "
                            f"report, or --allow-large in the pilot). The model is not skeletonized: trimming branches would "
                            f"delete the dead ends, where residuals are lowest.")
    rep.add(TEST_RUN_ROW, "EPANET ran the first hour (hydraulics and water age, set up as the app's first run) "
                          "without an error, and one full 7-day run of the grid's kind; the grid's runs come next",
            DERIVED)
    return True, replacement


# ----------------------------------------------------------------------------- command line
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--prepare", action="store_true", help="prepare each EPANET model and print its readiness report")
    ap.add_argument("--out", default=os.path.join("outputs", "private", "inputs"),
                    help="where a fixed copy goes (default outputs/private/inputs, git-ignored)")
    ap.add_argument("--disinfectant", choices=sorted(DISINFECTANT_WORDS), default="free_chlorine")
    ap.add_argument("--json", action="store_true", help="print JSON instead of plain text")
    a = ap.parse_args(argv)
    bad = 0
    for f in a.files:
        if a.prepare:
            path, rep = prepare_inp(f, a.out, a.disinfectant)
            print(rep.to_json() if a.json else rep.to_markdown())
            if path:
                print(f"prepared: {path}\n")
            bad += path is None
        else:
            with open(f, "rb") as fh:
                s = sniff(fh.read(), os.path.basename(f))
            print(json.dumps({"file": f, **s.as_dict()}) if a.json else f"{os.path.basename(f)}: {s.sentence()} {s.advice}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
