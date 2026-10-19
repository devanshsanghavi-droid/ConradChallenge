"""
importers.py: read any model file safely (iteration 5, task 15).

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
  is a blocking problem.  The report does not check the model's physics yet (headloss formula, reaction orders, bulk
  lines, start time; iteration 5, task 15b): it says so in an 'assumed' row, and warns when the headloss formula is
  not Hazen-Williams or a reaction order is not 1.
- Section headers and their names are read the way WNTR's reader reads them: after any whitespace (a form feed
  included), and a singular name such as [PIPE] is the [PIPES] section; the ID scan also stops at [END], as WNTR does.
- Zips are read in memory only, never extracted to disk.  A member with an absolute path or a '..' part is refused,
  and so is a zip whose members add up to more than 200 MB uncompressed.  The member that is read (the model) is
  refused if it is over 10 MB uncompressed and compressed more than 100 to 1; members that are only listed are not
  held to that ratio, because a shapefile's space-padded DBF packs past it at city size.  Files that macOS or Windows
  add to a zip (__MACOSX/, ._ copies, .DS_Store, Thumbs.db, desktop.ini) are left out of the listing.

Nothing here is applied inside simulate.load, so the bundled networks run bit for bit as before.  Uses the standard
library and WNTR only.
"""
from __future__ import annotations

import argparse
import codecs
import contextlib
import csv
import hashlib
import io
import json
import os
import re
import sqlite3
import struct
import sys
import tempfile
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
    position), strict UTF-8, then cp1252, where any of the five bytes cp1252 leaves undefined is read as latin-1 and
    every other byte as cp1252 (so a euro sign or a curly quote elsewhere in the file still reads right).  how is
    'utf-8' only for plain UTF-8 with no BOM, the one case that needs no rewrite."""
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
    try:
        return data.decode("cp1252"), "cp1252 (Windows Western European)"
    except UnicodeDecodeError:
        n = sum(data.count(bytes([b])) for b in CP1252_UNDEFINED)
        return (data.decode("cp1252", errors="residualmap_cp1252_gaps"),
                f"cp1252 (Windows Western European), with {n} byte{'s' if n != 1 else ''} it leaves undefined read as "
                f"latin-1")


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
    if how.startswith("cp1252"):
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
    junk = len(files) - len(names)
    if junk:
        s.detail += (f"{'; ' if s.detail else ''}{junk} file{'s' if junk > 1 else ''} that macOS or Windows adds "
                     f"(__MACOSX, ._ copies, .DS_Store, Thumbs.db) left out")
    return s


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
            bits = ["WNTR reads it and EPANET runs its first hour without an error; the full runs come next. Its "
                    "physics settings are used as the file has them and are not checked yet (the 'Physics settings' "
                    "row)."]
        else:                          # prepare_inp(parse=False)
            bits = ["Prepared without a check: WNTR has not read it and EPANET has not run it."]
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
        if self.fixes:
            out += ["**Fixes applied to a copy** (your file is not changed):"] + [f"- {f}" for f in self.fixes] + [""]
        if self.warnings:
            out += ["**Warnings:**"] + [f"- {w}" for w in self.warnings] + [""]
        return "\n".join(out).rstrip() + "\n"

    def as_dict(self) -> dict:
        return {"name": self.name, "disinfectant": self.disinfectant, "sha1": self.sha1, "sniff": self.sniff,
                "not_a_model": self.not_a_model, "ok": self.ok, "summary": self.summary(),
                "rows": [{"item": i, "value": v, "source": s} for i, v, s in self.rows], "blocking": self.blocking,
                "fixes": self.fixes, "warnings": self.warnings, "run_cost": self.run_cost, "prepared": self.prepared,
                "same_file": self.same_file}

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


def run_cost(n_junctions: int, disinfectant: str) -> dict:
    """The simulator grid for one dose: its EPANET runs and its cache file (ln C as float32 for every run, hour of the
    last day and junction)."""
    from .simgp import GRIDS
    grid = "chloramine" if disinfectant == "chloramine" else "full"
    runs = 1
    for axis in GRIDS[grid]:
        runs *= len(axis)
    mb = runs * 24 * n_junctions * FLOAT32_BYTES / 1e6
    return {"grid": grid, "runs_per_dose": runs, "grid_file_mb": round(mb, 1)}


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
    if how != "utf-8":
        what = {"utf-8 with a byte-order mark": "Removed the byte-order mark at the start (WNTR stops at line 1 on it)"}
        rep.fixes.append(what.get(how, f"Decoded the text from {how} and saved it as UTF-8 (WNTR reads UTF-8 only)")
                         + ".")
    ids = scan_inp_ids(text)
    rep.blocking += ids.blocking
    rep.warnings += ids.warnings
    if rep.blocking:
        return None, rep
    # write only when needed, under a temporary name until WNTR has read it, an entry point is found and the test
    # run passed; a temporary file left by any failure (the save itself included) is removed in the finally
    tmp, accepted, out_path = None, False, None
    try:
        if path is not None and how == "utf-8":
            out_path, rep.same_file = path, True
        else:
            new = data if how == "utf-8" else text.encode("utf-8")
            fname = prepared_name(_sha1(new)[:8], name, keep_name=how == "utf-8" and not rep.fixes)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, fname)
            if not os.path.exists(out_path):      # content-named: an existing file is this one, already accepted
                fd, tmp = tempfile.mkstemp(prefix=".prepare_", suffix=".inp", dir=out_dir)
                with os.fdopen(fd, "wb") as fh:
                    fh.write(new)
        rep.prepared = os.path.basename(out_path)
        accepted = True if not parse else _read_and_report(tmp or out_path, out_path, rep, ids, disinfectant)
        if accepted and tmp is not None:
            os.replace(tmp, out_path)
            tmp = None
    except OSError as e:
        accepted = False
        rep.blocking.append(f"The prepared copy could not be saved ({e.strerror or type(e).__name__}); nothing was "
                            f"kept.")
    finally:
        if tmp is not None:
            with contextlib.suppress(OSError):
                os.remove(tmp)
        if not accepted:
            rep.prepared, rep.same_file = None, False
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


def _read_and_report(read_path: str, out_path: str, rep: Readiness, ids: IdScan, disinfectant: str) -> bool:
    """WNTR's load of the prepared copy, the report rows that need it and EPANET's test run of the first hour; False
    (with a blocking problem) when WNTR cannot read it, the model has no entry point or the test run fails."""
    import wntr
    from .simulate import source_nodes
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wn = wntr.network.WaterNetworkModel(read_path)
    except Exception as e:  # noqa: BLE001
        msg = str(e).replace(read_path, os.path.basename(out_path))     # the file's name, never its folder
        rep.blocking.append(f"EPANET (through WNTR) could not read this model: {type(e).__name__}: {msg[:300]}")
        return False
    for item, n in (("Junctions", wn.num_junctions), ("Reservoirs", wn.num_reservoirs), ("Tanks", wn.num_tanks),
                    ("Pipes", wn.num_pipes), ("Pumps", wn.num_pumps), ("Valves", wn.num_valves)):
        rep.add(item, n, FROM_FILE)
    rep.add("Flow units", wn.options.hydraulic.inpfile_units, FROM_FILE)
    headloss = wn.options.hydraulic.headloss
    rep.add("Headloss formula", headloss, FROM_FILE)
    n_nodes = wn.num_nodes
    n_xy = len(ids.coordinate_ids & set(wn.node_name_list))
    if n_xy:
        rep.add("Node coordinates", f"{n_xy} of {n_nodes} nodes", FROM_FILE)
    else:
        rep.add("Node coordinates", "none: the maps need them (EPANET or the modelling package writes them in "
                                    "[COORDINATES])", MISSING)
    src_nodes = source_nodes(wn)
    if not src_nodes:
        rep.blocking.append("No reservoir, inflow junction or tank found: the model needs a place where water enters.")
        return False
    how_in = ("reservoirs" if wn.reservoir_name_list else
              "junctions with a negative demand (inflows)" if set(src_nodes) <= set(wn.junction_name_list) else "tanks")
    rep.add("Where water enters", f"{how_in}: {', '.join(src_nodes[:6])}{' and more' if len(src_nodes) > 6 else ''}",
            DERIVED)
    rep.add("Disinfectant", f"{DISINFECTANT_WORDS[disinfectant]} (your setting; every entry point gets the same dose)",
            ASSUMED)
    rep.add("Physics settings", "used as your file has them and not checked yet: the headloss formula, the reaction "
                                "orders, per-pipe and per-tank BULK lines, LIMITING POTENTIAL and START CLOCKTIME. A "
                                "model that does not use Hazen-Williams or first-order reactions runs, but its decay "
                                "rates do not mean what they say; these checks are planned, not built", ASSUMED)
    if headloss != "H-W":
        rep.warnings.append(f"The file uses the {headloss} headloss formula. The wall decay reads each pipe's roughness "
                            f"as a Hazen-Williams C, so on this file the wall rates and the old-pipe factor do not mean "
                            f"what they say. Converting the formula is planned, not built: export the model with "
                            f"Hazen-Williams if the modelling package can.")
    r = wn.options.reaction
    orders = [(w, v) for w, v in (("bulk", r.bulk_order), ("wall", r.wall_order), ("tank", r.tank_order))
              if v is not None and float(v) != 1.0]
    if orders:
        rep.warnings.append("The file sets the " + " and ".join(f"{w} reaction order to {float(v):g}" for w, v in orders)
                            + ". The decay model assumes first order (1) and does not change the file's order yet, so "
                              "its decay rates mean something else on this file.")
    rep.run_cost = run_cost(wn.num_junctions, disinfectant)
    rep.add("Run cost", f"{rep.run_cost['runs_per_dose']} EPANET runs per dose, once per file; a grid file of about "
                        f"{rep.run_cost['grid_file_mb']:g} MB", DERIVED)
    problem = _test_run(wn)            # last: it changes wn's time and quality options
    if problem:
        rep.blocking.append(problem)
        return False
    rep.add(TEST_RUN_ROW, "EPANET ran the first hour (hydraulics and water age, set up as the app's first run) "
                          "without an error; the full runs are longer and come next", DERIVED)
    return True


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
