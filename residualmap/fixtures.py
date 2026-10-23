"""
fixtures.py: the files the iteration-5 importer checks read, built from the bundled public networks (task 15).

    python -m residualmap.fixtures --build DIR     # write every fixture into DIR (for a look, or a browser upload)
    python -m residualmap.fixtures --record        # rerun the task-15 importer tests -> outputs/intake/task15_import.json

Nothing binary is committed: every fixture is built here, from WNTR's bundled networks (Net1, Net3) and from byte
literals for the header-only files (DWG, DGN, Jet, EPANET .bin).  The writers never share code with the readers they
test: the OLE compound file, the zips, the R12 DXF, the shapefile, the GeoPackage geometry, the PDF and the workbook
are written here by hand with the standard library, not by residualmap.importers.  The realistic quirks (a byte-order
mark, ';' comment lines first, a long title, UTF-16, a Windows code page, an emailed zip, a 999 comment before a DXF's
first section, a WNTR JSON, a State Plane GeoJSON with a crs member, Esri JSON, a settings file named .wtg, semicolons
and decimal commas) come from the iteration-5 review's realistic files.  No fixture names a real utility, and the
coordinates are invented.

SNIFF_TABLE lists the 36 fixtures of sniff_identifies_realistic_files with the kind each must come back as, whether
the file name may decide anything, and whether it is a backup copy.  TIE_BREAKS lists two more files whose content
cannot decide (an OLE header with no directory), where the extension must break the tie and the result must say so.
The REVIEW_* tables hold the cases the task-15 review found (a SWMM model and an EPANET-MSX file typed as models,
CR-only line ends and sections after 64 KB missed, zips made on a Mac, a city-size shapefile refused, malformed JSON
that raised, a run of XML comments, a csv.Sniffer sample and a GeoPackage view that took seconds to hours, an accented
ID too long in UTF-8, a long pattern ID, a copy left behind by a model WNTR could not read), and then the cases of
the task's final verification (singular core headers and a form feed before a header, which WNTR reads; a duplicated
pipe under [PIPE]; a junction no pipe connects, which EPANET cannot run; a zip with broken member offsets; a 263-byte
quote run that held csv.Sniffer for seconds; a 250-character upload name; a small model padded with NUL bytes).
They were added after the pre-registered tables above, which are unchanged, and each sits in the check it belongs
to.  The checks in residualmap/checks.py (group task15) rerun each part of record() and compare it with the
committed outputs/intake/task15_import.json.

Task 16 adds the grab-log fixtures (SAMPLE_FIXTURES): a LIMS export with one US datetime column, a log with '<0.02', 'ND'
and blanks, a log whose dates all read two ways, a cp1252 log and one with a byte-order mark, a workbook with typed date
and time cells and numeric tap IDs (written by openpyxl, the library that reads it: same-library self-consistency) with
its CSV twin, two online analyzers' 15-minute series, a log written 8 h late as a UTC export would be, an unknown tap,
both species, nine months, and a log written at Net3's model hours and again 6 h later for net3_clock6.inp, all with
one tap map.  The readings are invented and the junctions are Net3's; residualmap/intake_exp.py writes their record.
The review of task 16 and its final verification add the layouts they found (part 'review' and part 'final_review'
of that record): species in brackets, LIMS codes and analytes that cannot pick the rows, µg/L spellings, titles before
the header, Excel dates with no time, a superscript digit in a tap ID, a two-column junction export, interleaved
lab qualifiers, a tap map from another system, three million tiny rows, and a clock-time log for the bundled Net2.

Task 16c adds the zone-inlet fixtures (ZONE_FIXTURES): the task-16 tap map with a role column that makes one online
analyzer Lake's inlet analyzer, River's residual in an inlet residuals file, both residuals in one (and again in µg/L,
semicolon-separated, under headers to guess), and three files the reader must refuse (an inlet missing, a node that is
no entry point, a residual below detection).  Its review adds a declared inlet tap read outside the grab rows' window
and one read only below detection, role cells that mention an inlet but mark none, a tap marked as an inlet at a
junction, and three more inlet files to refuse (a thousands separator in µg/L, total chlorine, a tiny residual).
"""
from __future__ import annotations

import argparse
import functools
import io
import json
import math
import os
import sqlite3
import struct
import sys
import tempfile
import time
import warnings
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RECORD = os.path.join(REPO, "outputs", "intake", "task15_import.json")
ZIP_TIME = (1980, 1, 1, 0, 0, 0)       # every zip member gets this date, so a fixture's bytes never depend on the clock


def _lib(name: str) -> bytes:
    from .simulate import LIB
    with open(os.path.join(LIB, f"{name}.inp"), "rb") as fh:
        return fh.read()


def _net3() -> bytes:
    return _lib("Net3")


def _net3_text_lf() -> str:
    return _net3().decode("utf-8").replace("\r\n", "\n")


def _wn(name: str = "Net3"):
    import wntr
    from .simulate import LIB
    return wntr.network.WaterNetworkModel(os.path.join(LIB, f"{name}.inp"))


def _pipe_lines(name: str = "Net3", n: int = 12) -> list[tuple[tuple[float, float], tuple[float, float], float]]:
    """The first n pipes of a bundled network as ((x1, y1), (x2, y2), diameter in inches)."""
    wn = _wn(name)
    out = []
    for pn in wn.pipe_name_list[:n]:
        p = wn.get_link(pn)
        a, b = wn.get_node(p.start_node_name).coordinates, wn.get_node(p.end_node_name).coordinates
        out.append(((float(a[0]), float(a[1])), (float(b[0]), float(b[1])), round(p.diameter / 0.0254)))
    return out


def _lonlat(lines, lon0: float = -100.0, lat0: float = 40.0, span: float = 0.04):
    """Invented lon/lat: the drawing scaled into a box of span degrees at (lon0, lat0)."""
    xs = [c for ln in lines for c in (ln[0][0], ln[1][0])]
    ys = [c for ln in lines for c in (ln[0][1], ln[1][1])]
    sx = span / max(1e-9, max(xs) - min(xs))
    sy = span / max(1e-9, max(ys) - min(ys))
    f = lambda p: (round(lon0 + (p[0] - min(xs)) * sx, 6), round(lat0 + (p[1] - min(ys)) * sy, 6))
    return [(f(a), f(b), d) for a, b, d in lines]


# ----------------------------------------------------------------------------- writers (standard library only)
def _zip(members: list[tuple[str, bytes]], method: int = zipfile.ZIP_DEFLATED) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=method) as z:
        for name, data in members:
            zi = zipfile.ZipInfo(name, date_time=ZIP_TIME)
            zi.compress_type = method
            zi.external_attr = 0o644 << 16
            z.writestr(zi, data)
    return buf.getvalue()


def _cfb(entries: list[tuple[str, int, bytes]]) -> bytes:
    """A minimal OLE compound file (version 3, 512-byte sectors): one FAT sector, one directory sector (the root plus
    up to three entries, all children of the root), then each stream's data in whole sectors.  entries: (name, type,
    data), type 1 a storage, 2 a stream."""
    assert len(entries) <= 3
    ssz, ENDCHAIN, FREE, FATSECT, NOSTREAM = 512, 0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFD, 0xFFFFFFFF
    fat = [FATSECT, ENDCHAIN]
    starts, blobs = [], []
    for _, typ, data in entries:
        if typ != 2 or not data:
            starts.append(ENDCHAIN)
            continue
        n = -(-len(data) // ssz)
        first = len(fat)
        fat += [first + k + 1 for k in range(n - 1)] + [ENDCHAIN]
        starts.append(first)
        blobs.append(data + b"\x00" * (n * ssz - len(data)))
    assert len(fat) <= ssz // 4
    fat += [FREE] * (ssz // 4 - len(fat))

    def entry(name, typ, left, right, child, start, size):
        raw = name.encode("utf-16-le") + b"\x00\x00"
        e = raw.ljust(64, b"\x00") + struct.pack("<HBB", len(raw), typ, 1) + struct.pack("<III", left, right, child)
        e += b"\x00" * 16 + b"\x00" * 4 + b"\x00" * 16 + struct.pack("<I", start) + struct.pack("<Q", size)
        assert len(e) == 128
        return e

    dirs = [entry("Root Entry", 5, NOSTREAM, NOSTREAM, 1 if entries else NOSTREAM, ENDCHAIN, 0)]
    for k, ((name, typ, data), start) in enumerate(zip(entries, starts)):
        right = k + 2 if k + 1 < len(entries) else NOSTREAM
        dirs.append(entry(name, typ, NOSTREAM, right, NOSTREAM, start, len(data) if typ == 2 else 0))
    dirs += [b"\x00" * 64 + struct.pack("<HBB", 0, 0, 0) + struct.pack("<III", NOSTREAM, NOSTREAM, NOSTREAM)
             + b"\x00" * 52] * (4 - len(dirs))
    header = (b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" + b"\x00" * 16 + struct.pack("<HHHHH", 0x3E, 3, 0xFFFE, 9, 6)
              + b"\x00" * 6 + struct.pack("<IIIIIIIII", 0, 1, 1, 0, 4096, ENDCHAIN, 0, ENDCHAIN, 0)
              + struct.pack("<109I", 0, *([FREE] * 108)))
    assert len(header) == 512
    return header + struct.pack(f"<{ssz // 4}I", *fat) + b"".join(dirs) + b"".join(blobs)


def _dxf_r12(lines, layer: str = "WATER-MAINS", comments: tuple[str, ...] = ()) -> bytes:
    """An R12 ASCII DXF with one LINE per pipe and its size as a TEXT label, written by hand."""
    g = []
    for c in comments:
        g += ["999", c]
    g += ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1009", "9", "$INSUNITS", "70", "2", "0", "ENDSEC",
          "0", "SECTION", "2", "ENTITIES"]
    for (x1, y1), (x2, y2), d in lines:
        g += ["0", "LINE", "8", layer, "10", f"{x1:.3f}", "20", f"{y1:.3f}", "30", "0.0", "11", f"{x2:.3f}", "21",
              f"{y2:.3f}", "31", "0.0"]
        g += ["0", "TEXT", "8", "WATER-LABELS", "10", f"{(x1 + x2) / 2:.3f}", "20", f"{(y1 + y2) / 2 + 3:.3f}", "30",
              "0.0", "40", "2.5", "1", f'{d}" DIP']
    g += ["0", "ENDSEC", "0", "EOF"]
    return ("\r\n".join(f"{v:>3}" if k % 2 == 0 else v for k, v in enumerate(g)) + "\r\n").encode("ascii")


def _shp_polylines(lines) -> tuple[bytes, bytes]:
    """.shp and .shx for polylines (shape type 3), written by hand from the ESRI shapefile specification."""
    recs = []
    for (x1, y1), (x2, y2), _ in lines:
        box = struct.pack("<4d", min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
        recs.append(struct.pack("<i", 3) + box + struct.pack("<ii", 1, 2) + struct.pack("<i", 0)
                    + struct.pack("<4d", x1, y1, x2, y2))
    xs = [c for ln in lines for c in (ln[0][0], ln[1][0])]
    ys = [c for ln in lines for c in (ln[0][1], ln[1][1])]

    def header(length_words):
        return (struct.pack(">7i", 9994, 0, 0, 0, 0, 0, length_words) + struct.pack("<ii", 1000, 3)
                + struct.pack("<8d", min(xs), min(ys), max(xs), max(ys), 0, 0, 0, 0))

    body, index, offset = [], [], 50
    for i, r in enumerate(recs):
        body.append(struct.pack(">ii", i + 1, len(r) // 2) + r)
        index.append(struct.pack(">ii", offset, len(r) // 2))
        offset += 4 + len(r) // 2
    body, index = b"".join(body), b"".join(index)      # joined once: a city-size layer has 15,000 records
    return header(50 + len(body) // 2) + body, header(50 + len(index) // 2) + index


def _shp_points(points) -> tuple[bytes, bytes]:
    xs, ys = [p[0] for p in points], [p[1] for p in points]

    def header(length_words):
        return (struct.pack(">7i", 9994, 0, 0, 0, 0, 0, length_words) + struct.pack("<ii", 1000, 1)
                + struct.pack("<8d", min(xs), min(ys), max(xs), max(ys), 0, 0, 0, 0))

    body, index, offset = b"", b"", 50
    for i, (x, y) in enumerate(points):
        r = struct.pack("<idd", 1, x, y)
        body += struct.pack(">ii", i + 1, len(r) // 2) + r
        index += struct.pack(">ii", offset, len(r) // 2)
        offset += 4 + len(r) // 2
    return header(50 + len(body) // 2) + body, header(50 + len(index) // 2) + index


def _dbf(fields: list[tuple[str, str, int]], rows: list[list]) -> bytes:
    """A dBASE III table: fields (name, 'C' or 'N', width)."""
    rec_len = 1 + sum(w for _, _, w in fields)
    head = struct.pack("<BBBBIHH20x", 3, 126, 1, 1, len(rows), 32 + 32 * len(fields) + 1, rec_len)
    for name, typ, w in fields:
        head += name.encode("ascii")[:10].ljust(11, b"\x00") + typ.encode() + b"\x00" * 4 + bytes([w, 0]) + b"\x00" * 14
    head += b"\r"
    body = b"".join(b" " + b"".join((str(v).encode("latin-1")[:w].ljust(w) if t == "C" else str(v).encode()[:w].rjust(w))
                                    for v, (_, t, w) in zip(r, fields)) for r in rows)
    return head + body + b"\x1a"


PRJ_FT = ('PROJCS["NAD_1983_StatePlane_Example_FIPS_0000_Feet",GEOGCS["GCS_North_American_1983",DATUM["D_North_'
          'American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",'
          '0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["Standard_Parallel_1",37.0],'
          'PARAMETER["Standard_Parallel_2",39.0],PARAMETER["Central_Meridian",-98.0],PARAMETER["Latitude_Of_Origin",'
          '36.0],UNIT["Foot_US",0.3048006096012192]]')


def _gpkg_line(srs_id: int, pts) -> bytes:
    """A GeoPackage geometry blob: the 'GP' header (version 0, little-endian, no envelope) and a WKB LineString."""
    wkb = struct.pack("<BII", 1, 2, len(pts)) + b"".join(struct.pack("<dd", x, y) for x, y in pts)
    return b"GP" + bytes([0, 0b00000001]) + struct.pack("<i", srs_id) + wkb


def _sqlite(statements: list[tuple[str, tuple]]) -> bytes:
    con = sqlite3.connect(":memory:")
    try:
        for sql, args in statements:
            con.execute(sql, args)
        con.commit()
        return con.serialize()
    finally:
        con.close()


def _pdf() -> bytes:
    objs = [b"<< /Type /Catalog /Pages 2 0 R >>", b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Contents 4 0 R >>",
            b"<< /Length 44 >>\nstream\nBT /F1 12 Tf 72 720 Td (Map book) Tj ET\nendstream"]
    out, offs = b"%PDF-1.4\n", []
    for i, o in enumerate(objs):
        offs.append(len(out))
        out += f"{i + 1} 0 obj\n".encode() + o + b"\nendobj\n"
    x = len(out)
    out += f"xref\n0 {len(objs) + 1}\n0000000000 65535 f \n".encode()
    out += b"".join(f"{o:010d} 00000 n \n".encode() for o in offs)
    out += f"trailer\n<< /Size {len(objs) + 1} /Root 1 0 R >>\nstartxref\n{x}\n%%EOF\n".encode()
    return out


def _xlsx(rows: list[list]) -> bytes:
    """A minimal workbook (one sheet, inline strings and numbers), written by hand from the Office Open XML parts."""
    ct = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Types xmlns="http://schemas.openxmlformats.org/'
          'package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.'
          'relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/'
          'workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
          '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.'
          'spreadsheetml.worksheet+xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.'
            'openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats'
            '.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
    wb = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<workbook xmlns="http://schemas.openxmlformats.org/'
          'spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
          '<sheets><sheet name="grab log" sheetId="1" r:id="rId1"/></sheets></workbook>')
    wbrels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.'
              'openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.'
              'openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
              '</Relationships>')
    cells = []
    for r, row in enumerate(rows, start=1):
        cs = []
        for c, v in enumerate(row):
            ref = f"{chr(65 + c)}{r}"
            if isinstance(v, (int, float)):
                cs.append(f'<c r="{ref}"><v>{v}</v></c>')
            else:
                cs.append(f'<c r="{ref}" t="inlineStr"><is><t>{v}</t></is></c>')
        cells.append(f'<row r="{r}">{"".join(cs)}</row>')
    sheet = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<worksheet xmlns="http://schemas.openxmlformats.'
             f'org/spreadsheetml/2006/main"><sheetData>{"".join(cells)}</sheetData></worksheet>')
    return _zip([("[Content_Types].xml", ct.encode()), ("_rels/.rels", rels.encode()), ("xl/workbook.xml", wb.encode()),
                 ("xl/_rels/workbook.xml.rels", wbrels.encode()), ("xl/worksheets/sheet1.xml", sheet.encode())])


LOG_ROWS = [("2026-03-02", "09:15", "T1", "0,45"), ("2026-03-02", "10:40", "T2", "0,38"),
            ("2026-03-03", "08:55", "T3", "0,52"), ("2026-03-03", "11:20", "T4", "0,21")]


# ----------------------------------------------------------------------------- model text
def net3_as_txt() -> bytes:
    return _net3()


def net3_backup_inp_bak() -> bytes:
    return _net3()


def inp_with_bom() -> bytes:
    return b"\xef\xbb\xbf" + _net3()


def inp_comment_first() -> bytes:
    return b";Exported from a modelling package\r\n;Scenario: Base Average Day\r\n" + _net3()


def inp_long_title() -> bytes:
    n = _net3()
    return b"[TITLE]\r\n" + b"A" * 200 + b"\r\n" + n.split(b"\n", 1)[1]


def inp_utf16() -> bytes:
    return _net3().decode("utf-8").encode("utf-16")      # little-endian, with a byte-order mark


def net3_cp1252_title() -> bytes:
    t = _net3().decode("utf-8")
    return t.replace("[TITLE]\r\n", "[TITLE]\r\nZona de presión norte, cloración 1,2 mg/L, 50 °F\r\n", 1).encode("cp1252")


def emailed_model_zip() -> bytes:
    return _zip([("zone_model/zone7.inp", _net3())])


# ----------------------------------------------------------------------------- drawings
def dxf_with_999_comment() -> bytes:
    return _dxf_r12(_pipe_lines("Net3", 6), comments=("dxfrw 0.6.3",))


def ascii_drawing_dxf_bak() -> bytes:
    return _dxf_r12(_pipe_lines("Net3", 12))


def binary_dxf() -> bytes:
    return b"AutoCAD Binary DXF\r\n\x1a\x00" + b"\x00\x00SECTION\x00\x02\x00HEADER\x00" + b"\x00" * 48


def dwg_ac1027() -> bytes:
    return b"AC1027" + b"\x00" * 5 + b"\x1f\x00\x00\x00\x00\x00" + b"\x00" * 115


def dwg_ac1032() -> bytes:
    return b"AC1032" + b"\x00" * 5 + b"\x1f\x00\x00\x00\x00\x00" + b"\x00" * 115


def drawing_bak() -> bytes:
    return b"AC1024" + b"\x00" * 5 + b"\x1f\x00\x00\x00\x00\x00" + b"\x00" * 115


def dgn_v7() -> bytes:
    return bytes([0x08, 0x09, 0xFE, 0x02]) + b"\x00" * 1532


def ole_dgn_v8() -> bytes:
    return _cfb([("Dgn~H", 1, b""), ("Dgn~S", 2, b"\x00" * 1024)])


# ----------------------------------------------------------------------------- vendor model files
def jet_model_wtg_mdb() -> bytes:
    return (b"\x00\x01\x00\x00Standard Jet DB\x00" + b"\x01\x00\x00\x00").ljust(4096, b"\x00")


def model_wtg_sqlite() -> bytes:
    return _sqlite([("create table ModelElement (Id integer primary key, Label text, Kind integer)", ()),
                    ("create table Scenario (Id integer primary key, Label text)", ()),
                    ("insert into Scenario values (1, 'Base')", ())])


def wtg_settings() -> bytes:
    return (b'<?xml version="1.0" encoding="utf-8"?>\r\n<HaestadDataSource Version="10.4">\r\n'
            b'  <Settings Units="US Customary" />\r\n</HaestadDataSource>\r\n')


# ----------------------------------------------------------------------------- GIS
def mains_gpkg() -> bytes:
    lines = _pipe_lines("Net3", 8)
    st = [("PRAGMA application_id = 1196444487", ()), ("PRAGMA user_version = 10300", ()),
          ("create table gpkg_spatial_ref_sys (srs_name text, srs_id integer primary key, organization text, "
           "organization_coordsys_id integer, definition text, description text)", ()),
          ("insert into gpkg_spatial_ref_sys values ('Undefined cartesian', -1, 'NONE', -1, 'undefined', '')", ()),
          ("create table gpkg_contents (table_name text primary key, data_type text, identifier text, description "
           "text, last_change text, min_x real, min_y real, max_x real, max_y real, srs_id integer)", ()),
          ("insert into gpkg_contents values ('water_mains', 'features', 'water_mains', '', '2026-10-18T00:00:00Z', "
           "null, null, null, null, -1)", ()),
          ("create table gpkg_geometry_columns (table_name text, column_name text, geometry_type_name text, srs_id "
           "integer, z integer, m integer)", ()),
          ("insert into gpkg_geometry_columns values ('water_mains', 'geom', 'LINESTRING', -1, 0, 0)", ()),
          ("create table water_mains (fid integer primary key, geom blob, PIPEDIAM text)", ())]
    for k, (a, b, d) in enumerate(lines):
        st.append(("insert into water_mains values (?, ?, ?)", (k + 1, _gpkg_line(-1, [a, b]), f'{d}"')))
    return _sqlite(st)


def _feature_collection(lines, crs: str | None = None) -> dict:
    fc = {"type": "FeatureCollection", "name": "water_mains"}
    if crs:
        fc["crs"] = {"type": "name", "properties": {"name": crs}}
    fc["features"] = [{"type": "Feature", "properties": {"PIPEID": f"P{k + 1}", "PIPEDIAM": f"{d} in"},
                       "geometry": {"type": "LineString", "coordinates": [list(a), list(b)]}}
                      for k, (a, b, d) in enumerate(lines)]
    return fc


def mains_geojson_as_json() -> bytes:
    return json.dumps(_feature_collection(_lonlat(_pipe_lines("Net3", 12)))).encode()


def geojson_stateplane() -> bytes:
    lines = [((6.0e6 + 50 * a[0], 2.1e6 + 50 * a[1]), (6.0e6 + 50 * b[0], 2.1e6 + 50 * b[1]), d)
             for a, b, d in _pipe_lines("Net3", 12)]
    return json.dumps(_feature_collection(lines, "urn:ogc:def:crs:EPSG::2227")).encode()


def esri_features() -> bytes:
    feats = [{"attributes": {"OBJECTID": k + 1, "DIAM": d},
              "geometry": {"paths": [[[6.0e6 + 50 * a[0], 2.1e6 + 50 * a[1]], [6.0e6 + 50 * b[0], 2.1e6 + 50 * b[1]]]]}}
             for k, (a, b, d) in enumerate(_pipe_lines("Net3", 6))]
    return json.dumps({"displayFieldName": "", "geometryType": "esriGeometryPolyline",
                       "spatialReference": {"wkid": 2227, "latestWkid": 2227}, "fields": [], "features": feats}).encode()


def wntr_sidecar() -> bytes:
    import wntr
    return json.dumps(wntr.network.to_dict(_wn("Net3"))).encode()


def gml_mains() -> bytes:
    members = "".join(
        f'<gml:featureMember><water_mains><geom><gml:LineString srsName="EPSG:2227"><gml:posList>'
        f'{6.0e6 + 50 * a[0]:.1f} {2.1e6 + 50 * a[1]:.1f} {6.0e6 + 50 * b[0]:.1f} {2.1e6 + 50 * b[1]:.1f}'
        f'</gml:posList></gml:LineString></geom><DIAM>{d}</DIAM></water_mains></gml:featureMember>'
        for a, b, d in _pipe_lines("Net3", 6))
    return ('<?xml version="1.0" encoding="utf-8"?>\n<gml:FeatureCollection xmlns:gml="http://www.opengis.net/gml">'
            + members + "</gml:FeatureCollection>\n").encode()


def mains_kml() -> bytes:
    marks = "".join(
        f"<Placemark><name>main {k + 1}</name><description>{d} in</description><LineString><coordinates>"
        f"{a[0]},{a[1]},0 {b[0]},{b[1]},0</coordinates></LineString></Placemark>"
        for k, (a, b, d) in enumerate(_lonlat(_pipe_lines("Net3", 10))))
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
            f"<name>water mains</name>{marks}</Document></kml>\n").encode()


def mains_kmz() -> bytes:
    return _zip([("doc.kml", mains_kml())])


def mains_shp_nested_two_layers_zip() -> bytes:
    lines = _pipe_lines("Net3", 10)
    shp, shx = _shp_polylines(lines)
    dbf = _dbf([("PIPEID", "C", 10), ("PIPEDIAM", "C", 8)], [[f"P{k + 1}", f"{d}"] for k, (_, _, d) in enumerate(lines)])
    pts = [a for a, _, _ in lines[:5]]
    pshp, pshx = _shp_points(pts)
    pdbf = _dbf([("HYDID", "C", 8)], [[f"H{k + 1}"] for k in range(len(pts))])
    root = "export/water/"
    return _zip([(root + "mains.shp", shp), (root + "mains.shx", shx), (root + "mains.dbf", dbf),
                 (root + "mains.prj", PRJ_FT.encode()), (root + "mains.cpg", b"UTF-8"),
                 (root + "hydrants.shp", pshp), (root + "hydrants.shx", pshx), (root + "hydrants.dbf", pdbf),
                 (root + "hydrants.prj", PRJ_FT.encode())])


def utility_gdb_zip() -> bytes:
    root = "utility.gdb/"
    return _zip([(root + "a00000001.gdbtable", b"\x03\x00\x00\x00" + b"\x00" * 36),
                 (root + "a00000001.gdbtablx", b"\x03\x00\x00\x00" + b"\x00" * 12),
                 (root + "a00000009.gdbtable", b"\x03\x00\x00\x00" + b"\x00" * 36), (root + "gdb", b"\x05\x00"),
                 (root + "timestamps", b"\xff" * 16)])


# ----------------------------------------------------------------------------- other
def grab_log_xlsx() -> bytes:
    return _xlsx([["date", "time", "tap_id", "free_chlorine_mgL"]] +
                 [[d, t, k, float(v.replace(",", "."))] for d, t, k, v in LOG_ROWS])


def old_excel_log_xls() -> bytes:
    return _cfb([("Workbook", 2, b"\x09\x08\x10\x00\x00\x06\x05\x00" + b"\x00" * 4088)])


def log_semicolon() -> bytes:
    return ("date;time;tap_id;free_chlorine_mgL\n" + "".join(";".join(r) + "\n" for r in LOG_ROWS)).encode()


def log_tab() -> bytes:
    return ("date\ttime\ttap_id\tfree_chlorine_mgL\n" +
            "".join("\t".join(r[:3] + (r[3].replace(",", "."),)) + "\n" for r in LOG_ROWS)).encode()


def map_book() -> bytes:
    return _pdf()


def plant_ifc() -> bytes:
    return ("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1');\n"
            "FILE_NAME('plant.ifc','2026-10-18T00:00:00',(''),(''),'','','');\nFILE_SCHEMA(('IFC4'));\nENDSEC;\nDATA;\n"
            "#1=IFCPROJECT('0YvctVUKr0kugbFTf53O9L',$,'Treatment plant',$,$,$,$,$,$);\nENDSEC;\n"
            "END-ISO-10303-21;\n").encode()


def results_bin() -> bytes:
    wn = _wn("Net3")
    prolog = struct.pack("<15i", 516114521, 20200, wn.num_nodes, wn.num_tanks + wn.num_reservoirs, wn.num_links,
                         wn.num_pumps, wn.num_valves, 1, 0, 1, 0, 0, 0, 3600, 7 * 86400)
    return prolog + b"\x00" * 1024


# ----------------------------------------------------------------------------- not for the sniff table
def ole_header_only() -> bytes:
    """An OLE header with no directory sector after it (as the review's .xls arrived): its content cannot say what it
    is."""
    return b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" + b"\x00" * 504


def net3_crlf() -> bytes:
    return _net3_text_lf().replace("\n", "\r\n").encode("utf-8")


def net3_lf() -> bytes:
    return _net3_text_lf().encode("utf-8")


def net3_copy() -> bytes:
    return _net3()


def _insert_after_header(text: str, section: str, line: str) -> str:
    """Insert a data line right after a section's header line and its column comment."""
    i = text.index(f"[{section}]\n")
    j = text.index("\n", i) + 1
    if text[j:].startswith(";"):
        j = text.index("\n", j) + 1
    return text[:j] + line + "\n" + text[j:]


def net3_dup_junction() -> bytes:
    return _insert_after_header(_net3_text_lf(), "JUNCTIONS", " 10              \t99          \t0").encode()


def net3_junction_tank_collision() -> bytes:
    return _insert_after_header(_net3_text_lf(), "JUNCTIONS", " 1               \t99          \t0").encode()


def net3_dup_pipe() -> bytes:
    t = _net3_text_lf()
    i = t.index("[PIPES]\n")
    j = t.index("\n", t.index("\n", i) + 1) + 1          # past the header and the column comment
    first = t[j:t.index("\n", j)]
    return (t[:j] + first + "\n" + t[j:]).encode()


def _rename_node(t: str, old: str, new: str) -> str:
    """Rename one node everywhere a node ID can stand (its own line, the link ends, coordinates, demands, quality,
    sources, emitters, mixing, tags, controls and rules), and nowhere else: Net3 also has a pump and option values
    written '10'."""
    first = {"JUNCTIONS", "RESERVOIRS", "TANKS", "COORDINATES", "DEMANDS", "QUALITY", "SOURCES", "EMITTERS", "MIXING"}
    out, sec = [], None
    for ln in t.split("\n"):
        body = ln.split(";", 1)[0]
        tok = body.split()
        if body.strip().startswith("["):
            sec = body.strip().upper()[1:-1]
        elif tok:
            idx = []
            if sec in first:
                idx = [0]
            elif sec in ("PIPES", "PUMPS", "VALVES"):
                idx = [1, 2]
            elif sec == "TAGS" and tok[0].upper() == "NODE":
                idx = [1]
            elif sec in ("CONTROLS", "RULES"):
                idx = [k + 1 for k, w in enumerate(tok[:-1]) if w.upper() in ("NODE", "JUNCTION", "TANK")]
            if any(k < len(tok) and tok[k] == old for k in idx):
                ln = " " + "\t".join(new if (k in idx and w == old) else w for k, w in enumerate(tok))
        out.append(ln)
    return "\n".join(out)


def net3_long_id_32() -> bytes:
    return _rename_node(_net3_text_lf(), "10", "J" * 32).encode()


def net3_long_id_31() -> bytes:
    return _rename_node(_net3_text_lf(), "10", "J" * 31).encode()


def zip_slip() -> bytes:
    return _zip([("../../evil.inp", _net3())])


def zip_absolute() -> bytes:
    return _zip([("/tmp/evil.inp", _net3())])


def zip_ratio_bomb() -> bytes:
    return _zip([("model.inp", b"0" * (20 * 1024 * 1024))])


def zip_two_models() -> bytes:
    return _zip([("a/zone1.inp", _net3()), ("b/zone2.inp", _net3())])


def shp_padded_dbf_zip() -> bytes:
    """A zipped shapefile whose DBF (1,156 records, three blank 254-character text fields) packs far more than 100 to 1
    but unpacks to under 10 MB, so the ratio rule must let it through."""
    lines = _pipe_lines("Net3", 12) * 97
    lines = lines[:1156]
    shp, shx = _shp_polylines(lines)
    dbf = _dbf([("PIPEID", "C", 10), ("NOTES", "C", 254), ("OWNER", "C", 254), ("WORKORDER", "C", 254)],
               [[f"P{k + 1}", "", "", ""] for k in range(len(lines))])
    return _zip([("mains.shp", shp), ("mains.shx", shx), ("mains.dbf", dbf), ("mains.prj", PRJ_FT.encode())])


# ----------------------------------------------------------------------------- the review's cases (task 15)
CITY_RECORDS = 15000                   # a city's mains layer: past the 10 MB line, the DBF still packs past 100 to 1


@functools.lru_cache(maxsize=1)
def _city_mains() -> tuple[tuple[str, bytes], ...]:
    """A city-size mains layer as GIS packages export it: CITY_RECORDS polylines and a DBF with three blank
    254-character text fields."""
    lines = (_pipe_lines("Net3", 12) * (CITY_RECORDS // 12 + 1))[:CITY_RECORDS]
    shp, shx = _shp_polylines(lines)
    dbf = _dbf([("PIPEID", "C", 10), ("NOTES", "C", 254), ("OWNER", "C", 254), ("WORKORDER", "C", 254)],
               [[f"P{k + 1}", "", "", ""] for k in range(len(lines))])
    return (("mains.shp", shp), ("mains.shx", shx), ("mains.dbf", dbf), ("mains.prj", PRJ_FT.encode()))


def city_mains_shp_zip() -> bytes:
    return _zip(list(_city_mains()))


def model_and_gis_zip() -> bytes:
    """The model and its GIS layer sent in one zip: only the model is read."""
    return _zip([("model/zone3.inp", _net3())] + [("gis/" + n, b) for n, b in _city_mains()])


def _appledouble() -> bytes:
    """An AppleDouble header (magic, version 2, the 'Mac OS X' filler, no entries), as in the ._ copies that Finder's
    Compress puts under __MACOSX/ for every file with extended attributes (a downloaded file has a quarantine one)."""
    return struct.pack(">II", 0x00051607, 0x00020000) + b"Mac OS X        " + struct.pack(">H", 0)


def finder_model_zip() -> bytes:
    """A model zipped by Finder's Compress: its folder, the model, and the __MACOSX/ folder with the ._ copy."""
    return _zip([("Zone3/", b""), ("Zone3/zone.inp", _lib("Net1")), ("__MACOSX/", b""), ("__MACOSX/Zone3/", b""),
                 ("__MACOSX/Zone3/._zone.inp", _appledouble())])


def finder_shapefile_zip() -> bytes:
    """A shapefile zipped on a Mac, with __MACOSX/ copies of every part and a .DS_Store in its folder."""
    lines = _pipe_lines("Net3", 8)
    shp, shx = _shp_polylines(lines)
    dbf = _dbf([("PIPEID", "C", 10)], [[f"P{k + 1}"] for k in range(len(lines))])
    parts = [("shp", shp), ("shx", shx), ("dbf", dbf), ("prj", PRJ_FT.encode())]
    return _zip([("mains/", b"")] + [(f"mains/mains.{e}", b) for e, b in parts]
                + [("mains/.DS_Store", b"\x00\x00\x00\x01Bud1" + b"\x00" * 60)]
                + [(f"__MACOSX/mains/._mains.{e}", _appledouble()) for e, _ in parts])


def net3_cr_only() -> bytes:
    return _net3_text_lf().replace("\n", "\r").encode()           # classic Mac line ends


def net3_cr_cr_lf() -> bytes:
    return _net3_text_lf().replace("\n", "\r\r\n").encode()       # CRLF converted to CRLF a second time


def inp_title_over_64kb() -> bytes:
    """A [TITLE] of 1,000 note lines before the first section, so no node or link section is in the first 64 KB."""
    notes = "".join(f"Note {k:04d}: scenario notes kept by the modelling package, one line for each edit.\r\n"
                    for k in range(1000))
    return b"[TITLE]\r\n" + notes.encode() + _net3().split(b"\n", 1)[1]


def inp_patterns_first() -> bytes:
    """A large [PATTERNS] block (1,500 patterns) written before [JUNCTIONS], as some exports order the sections."""
    t = _net3_text_lf()
    i = t.index("[JUNCTIONS]")
    block = ("[PATTERNS]\n;ID              \tMultipliers\n"
             + "".join(f" Q{k:04d}           \t" + "\t".join(["1.0"] * 12) + "\n" for k in range(1500)) + "\n")
    return (t[:i] + block + t[i:]).replace("\n", "\r\n").encode()


def inp_comment_header_70kb() -> bytes:
    head = "".join(f";Exported from a modelling package, scenario notes line {k:04d}, kept as comments\r\n"
                   for k in range(900))
    return head.encode() + _net3()


def swmm_sewer_inp() -> bytes:
    """A small SWMM storm-drain model (invented): it shares [JUNCTIONS] with EPANET and has SWMM's own sections."""
    return ("[TITLE]\r\nStorm drain, invented\r\n\r\n[OPTIONS]\r\nFLOW_UNITS           CFS\r\nINFILTRATION         HORTON\r\n"
            "FLOW_ROUTING         KINWAVE\r\n\r\n[JUNCTIONS]\r\n;;Name  Elevation  MaxDepth  InitDepth  SurDepth  Aponded\r\n"
            "J1  96  4  0  0  0\r\nJ2  90  4  0  0  0\r\n\r\n[OUTFALLS]\r\n;;Name  Elevation  Type  Stage  Gated\r\n"
            "Out1  85  FREE    NO\r\n\r\n[CONDUITS]\r\n;;Name  From  To  Length  Roughness  InOffset  OutOffset\r\n"
            "C1  J1  J2  400  0.01  0  0\r\nC2  J2  Out1  400  0.01  0  0\r\n\r\n[XSECTIONS]\r\n"
            ";;Link  Shape  Geom1  Geom2  Geom3  Geom4  Barrels\r\nC1  CIRCULAR  1  0  0  0  1\r\n"
            "C2  CIRCULAR  1  0  0  0  1\r\n\r\n[COORDINATES]\r\nJ1  0  0\r\nJ2  100  0\r\nOut1  200  0\r\n").encode()


def chlorine_decay_msx() -> bytes:
    """A two-species EPANET-MSX reaction file (bulk and wall chlorine), the kind that sits next to a model's .inp: it
    has [PIPES] and [TANKS] sections too, with RATE lines instead of pipes."""
    return ("[TITLE]\r\nChlorine decay, bulk and wall\r\n\r\n[OPTIONS]\r\n AREA_UNITS FT2\r\n RATE_UNITS DAY\r\n"
            " SOLVER RK5\r\n TIMESTEP 300\r\n\r\n[SPECIES]\r\n BULK CL2 MG 0.01 0.001\r\n WALL CLW MG 0.01 0.001\r\n\r\n"
            "[COEFFICIENTS]\r\n PARAMETER Kb 0.3\r\n PARAMETER Kw 1.0\r\n\r\n[TERMS]\r\n Kf 1.5826e-4 * RE^0.88 / D\r\n\r\n"
            "[PIPES]\r\n RATE CL2 -Kb*CL2 - (4/D)*Kw*Kf/(Kw+Kf)*CL2\r\n RATE CLW 0\r\n\r\n[TANKS]\r\n RATE CL2 -Kb*CL2\r\n"
            "\r\n[QUALITY]\r\n GLOBAL CL2 0.5\r\n\r\n[REPORT]\r\n SPECIES CL2 YES\r\n").encode()


def _small_fc() -> dict:
    return _feature_collection(_lonlat(_pipe_lines("Net3", 4)))


def geojson_crs_string() -> bytes:
    fc = _small_fc()
    fc["crs"] = {"type": "name", "properties": "EPSG:2227"}          # properties should be an object
    return json.dumps(fc).encode()


def geojson_crs_name_number() -> bytes:
    fc = _small_fc()
    fc["crs"] = {"type": "name", "properties": {"name": 2227}}       # name should be a string
    return json.dumps(fc).encode()


def geojson_features_not_list() -> bytes:
    return json.dumps({"type": "FeatureCollection", "features": 3}).encode()


def esri_sr_string() -> bytes:
    return json.dumps({"geometryType": "esriGeometryPolyline", "spatialReference": "102100", "fields": [],
                       "features": []}).encode()


def geojson_deep_geometry() -> bytes:
    return ('{"type": "Feature", "properties": {}, "geometry": {"type": "LineString", "coordinates": '
            + "[" * 3000 + "1, 2" + "]" * 3000 + "}}").encode()


def json_deep() -> bytes:
    return b"[" * 100000 + b"]" * 100000


def kml_doctype_subset() -> bytes:
    """A KML file with a DOCTYPE whose internal subset holds '>' in a quoted value and a comment with a quote mark."""
    return (b'<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE kml [\n  <!ENTITY note "pipe size > 12 in">\n'
            b"  <!-- the team's own entity -->\n]>\n<!-- mains -->\n" + mains_kml().split(b"\n", 1)[1])


XML_COMMENT_RUN = 26                   # the review's case: about 8 s before the linear scan, doubling per comment


def xml_comment_run() -> bytes:
    return b"<!-- x -->" * XML_COMMENT_RUN + b"not an element"


def xml_comment_run_long() -> bytes:
    return b"<!-- x -->" * 6000 + b"not an element"


def xml_pi_run() -> bytes:
    return b'<?xml version="1.0"?>' + b"<?p x?>" * 6000 + b"no root"


def slow_view_gpkg() -> bytes:
    """A 4 KB SQLite file whose gpkg_contents is a view over a recursive count to 1e11 that returns no rows: reading
    the view would keep a query busy for hours."""
    return _sqlite([("create view gpkg_contents as with recursive c(x) as (select 1 union all select x + 1 from c "
                     "where x < 100000000000) select 'mains' as table_name, 'features' as data_type from c where x < 0",
                     ())])


def csv_quotes() -> bytes:
    """Two 30 KB lines of a,"a,"... : csv.Sniffer's quote regexes take seconds on this when given all of it."""
    return (('a,"' * 10000) + "\n").encode() * 2


def net3_cp1252_undefined_byte() -> bytes:
    """cp1252 with a euro sign and curly quotes in the title, and one byte (0x81) that cp1252 leaves undefined."""
    title = "[TITLE]\r\nZona norte, coste 12 € al mes, “revisado”".encode("cp1252") + b" \x81\r\n"
    return title + _net3().split(b"\n", 1)[1]


def net1_accented_id_cp1252() -> bytes:
    """Net1 with junction 22 renamed to a 31-character label with two accented letters, saved in cp1252: 31 bytes
    there, 33 in the UTF-8 copy EPANET reads."""
    return _rename_node(_lib("Net1").decode("utf-8").replace("\r\n", "\n"), "22",
                        "Hydrant_Grünwald_Straße_000022x").encode("cp1252")


def net3_long_pattern_id() -> bytes:
    return _insert_after_header(_net3_text_lf(), "PATTERNS", " " + "P" * 33 + "\t1.0\t1.0\t1.0").encode()


def net3_bom_missing_node() -> bytes:
    """A byte-order mark (so a fixed copy is written) and a pipe to a node the file does not have (so WNTR stops)."""
    return b"\xef\xbb\xbf" + _insert_after_header(_net3_text_lf(), "PIPES",
                                                  " 9999            \t10              \tNOPE            \t100"
                                                  "         \t12          \t100         \t0           \tOpen").encode()


# ----------------------------------------------------------------------------- the final verification's cases (task 15)
CORE_HEADERS = ("JUNCTIONS", "PIPES", "TANKS", "RESERVOIRS")


def net3_singular_headers() -> bytes:
    """Net3 with its four core headers written singular ([JUNCTION], [PIPE], [TANK], [RESERVOIR]): WNTR's reader
    allows a missing or an extra S, so it reads and runs this file, which the sniff had refused as plain text."""
    t = _net3_text_lf()
    for h in CORE_HEADERS:
        t = t.replace(f"[{h}]", f"[{h[:-1]}]")
    return t.replace("\n", "\r\n").encode()


def net3_form_feed_headers() -> bytes:
    """Net3 with a form feed (a page break some exports write) before each core header: WNTR strips it like any
    other whitespace."""
    t = _net3_text_lf()
    for h in CORE_HEADERS:
        t = t.replace(f"[{h}]", f"\x0c[{h}]")
    return t.replace("\n", "\r\n").encode()


def net3_singular_pipe_dup() -> bytes:
    """[PIPE] written singular, and the first pipe listed a second time with its ends swapped: WNTR reads [PIPE] as
    [PIPES] and keeps the later line without a word, so the duplicate-ID rule has to read [PIPE] too."""
    t = _net3_text_lf().replace("[PIPES]", "[PIPE]")
    i = t.index("[PIPE]\n")
    j = t.index("\n", t.index("\n", i) + 1) + 1          # past the header and the column comment
    first = t[j:t.index("\n", j)]
    tok = first.split()
    swapped = " " + "\t".join([tok[0], tok[2], tok[1]] + tok[3:])
    k = t.index("\n", j) + 1
    return (t[:k] + swapped + "\n" + t[k:]).replace("\n", "\r\n").encode()


def net3_isolated_junction() -> bytes:
    """Net3 plus junction 9999, with coordinates and no pipe (a node left alone is common in an exported zone model):
    WNTR reads it, and EPANET then stops with its Error 200, which the app showed as a traceback under a report
    that said the model 'Reads and runs'."""
    t = _insert_after_header(_net3_text_lf(), "JUNCTIONS", " 9999            \t100         \t10          \t"
                                                          "                \t;")
    t = _insert_after_header(t, "COORDINATES", " 9999            \t5.00            \t5.00")
    return t.replace("\n", "\r\n").encode()


ZIP_BAD_OFFSET = 27                    # bytes added to the end record's central-directory offset


def zip_bad_offset() -> bytes:
    """The emailed zip with the central-directory offset in its end record raised by ZIP_BAD_OFFSET bytes: zipfile
    still lists it (it finds the directory from the end), but every member's offset becomes negative, and reading the
    model raised ValueError ('negative seek value') out of prepare_inp, a traceback in the app."""
    z = bytearray(emailed_model_zip())
    end = z.rindex(b"PK\x05\x06")
    struct.pack_into("<I", z, end + 16, struct.unpack_from("<I", z, end + 16)[0] + ZIP_BAD_OFFSET)
    return bytes(z)


CSV_QUOTE_RUN = 250                    # quote marks on the middle line


def csv_quote_run() -> bytes:
    """Three short lines, the middle one CSV_QUOTE_RUN quote marks: 263 bytes on which csv.Sniffer took seconds, all
    in its quote regexes, while holding the interpreter lock (every app session waited)."""
    return (',"x",\n,' + '"' * CSV_QUOTE_RUN + "a\nb,c\n").encode()


def net1_nul_padded() -> bytes:
    """Net1 (6 KB) followed by 200 NUL bytes after its [END], as a copy padded to a block size arrives: the NULs made
    the small text look binary to the sniff, and WNTR, which stops reading at [END], runs it."""
    return _lib("Net1").rstrip() + b"\r\n" + b"\x00" * 200


# ----------------------------------------------------------------------------- task 15b: physics settings
# Each file below is Net3 (or ky10) with one physics setting changed, written by these helpers, which share no code
# with residualmap.importers (its reaction, headloss and structure checks are what they test).
def _section_edit(text: str, section: str, fn) -> str:
    """Apply fn(tokens, line) to every data line of one section of an LF text: it returns the new line, the same
    line, or None to drop it.  Sections are matched by their header as Net3 and ky10 write them."""
    out, sec = [], None
    for ln in text.split("\n"):
        s = ln.split(";", 1)[0].strip()
        if s.startswith("["):
            sec = s.upper()
            out.append(ln)
            continue
        if sec == f"[{section}]" and s:
            new = fn(s.split(), ln)
            if new is not None:
                out.append(new)
            continue
        out.append(ln)
    return "\n".join(out)


def _section_ids(text: str, section: str) -> list[str]:
    ids = []
    _section_edit(text, section, lambda tok, ln: ids.append(tok[0]) or ln)
    return ids


def _crlf(text: str) -> bytes:
    return text.replace("\n", "\r\n").encode("utf-8")


BULK_PIPES = 10                        # pipes given their own BULK line in net3_bulk_lines.inp


def net3_bulk_lines() -> bytes:
    """Net3 with per-pipe BULK lines on its first BULK_PIPES pipes, a per-tank TANK line on tank 1 and LIMITING
    POTENTIAL 0.5 (critique 9a): WNTR reads them, and those pipes and the tank would ignore the bulk rate the grid
    sets, while the limiting potential changes the first-order law."""
    t = _net3_text_lf()
    lines = [f" Bulk            \t{p:<16}\t-0.9" for p in _section_ids(t, "PIPES")[:BULK_PIPES]]
    t = _insert_after_header(t, "REACTIONS", "\n".join(lines + [" Tank            \t1               \t-0.3"]))
    assert " Limiting Potential    \t0.0" in t
    return _crlf(t.replace(" Limiting Potential    \t0.0", " Limiting Potential    \t0.5"))


FX_G, FX_NU_FT2S = 9.81, 1.1e-5         # the fixture's own constants: m/s2, and EPANET's water viscosity at 20 C (ft2/s)
DW_TINY_MFT = 0.001                    # millifeet, for a pipe no positive roughness can match (WNTR refuses 0)


@functools.lru_cache(maxsize=1)
def _net3_reference_velocity() -> dict:
    """Net3's own reference velocity per pipe (m/s): the median |V| over the last day of a 7-day run, at least 0.3."""
    import numpy as np
    import wntr
    wn = _wn("Net3")
    t = wn.options.time
    t.duration, t.hydraulic_timestep, t.report_timestep = 7 * 86400, 3600, 3600
    with tempfile.TemporaryDirectory(prefix="rm_epanet_") as d, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(d, "v")).link["velocity"]
    v = v[(v.index >= 6 * 86400) & (v.index < 7 * 86400)]
    return {pn: max(float(np.median(np.abs(v[pn].values))), 0.3) for pn in wn.pipe_name_list}


@functools.lru_cache(maxsize=1)
def net3_dw_roughness() -> tuple[dict, tuple]:
    """({pipe: Darcy-Weisbach roughness in millifeet}, pipes no positive roughness can match): the inverse of the stated
    formula at Net3's own reference velocities.  The target friction factor is the one whose Darcy-Weisbach headloss
    equals Net3's Hazen-Williams headloss there, f = 10.67 Q^1.852 (2 g D) / (C^1.852 D^4.8704 V^2), and Swamee-Jain
    is inverted for the roughness, e = 3.7 D (10^(-1 / (2 sqrt f)) - 5.74 / Re^0.9).  A pipe smoother than a smooth
    pipe at its velocity (Net3's C 199 pipes) gets DW_TINY_MFT."""
    wn = _wn("Net3")
    v_ref = _net3_reference_velocity()
    nu = FX_NU_FT2S * 0.3048 ** 2 * float(wn.options.hydraulic.viscosity)
    out, tiny = {}, []
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        D, C, V = float(p.diameter), float(p.roughness), v_ref[pn]
        q = V * math.pi * D * D / 4
        f = 10.67 * q ** 1.852 * 2 * FX_G * D / (C ** 1.852 * D ** 4.8704 * V * V)
        e = 3.7 * D * (10 ** (-1 / (2 * math.sqrt(f))) - 5.74 / (V * D / nu) ** 0.9)
        if e <= 0:
            tiny.append(pn)
            out[pn] = DW_TINY_MFT
        else:
            out[pn] = e / 0.3048 * 1000
    return out, tuple(tiny)


def _set_token(tok: list[str], ln: str, k: int, value: str) -> str:
    tok = list(tok)
    tok[k] = value
    return " " + "\t".join(tok) + ("\t;" + ln.split(";", 1)[1] if ";" in ln else "")


def net3_dw() -> bytes:
    """Net3 converted to Darcy-Weisbach by the inverse of the stated formula at Net3's own reference velocities
    (net3_dw_roughness), so a converter that recovers each C at the reference velocity it measures on this file should
    give Net3's C back, except for the pipes no positive roughness matches."""
    rough, _ = net3_dw_roughness()
    t = _net3_text_lf()
    assert " Headloss           \tH-W" in t
    t = t.replace(" Headloss           \tH-W", " Headloss           \tD-W")
    t = _section_edit(t, "PIPES", lambda tok, ln: _set_token(tok, ln, 5, f"{rough[tok[0]]:.6g}"))
    return _crlf(t)


def net3_cm() -> bytes:
    """Net3 with its headloss formula set to Chezy-Manning (its C values left as they are): always blocked."""
    return _crlf(_net3_text_lf().replace(" Headloss           \tH-W", " Headloss           \tC-M"))


def net3_clock6() -> bytes:
    """Net3 starting at 6 am (START CLOCKTIME)."""
    t = _net3_text_lf()
    assert " Start ClockTime    \t12 am" in t
    return _crlf(t.replace(" Start ClockTime    \t12 am", " Start ClockTime    \t6 am"))


WEEK_DAY_FACTORS = (1.0, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8)   # net3_weekly_pattern.inp: a lighter weekend


def net3_weekly_pattern() -> bytes:
    """Net3 with its default demand pattern 1 (24 hourly steps) made weekly: 168 hourly steps, the day's shape times
    WEEK_DAY_FACTORS."""
    t = _net3_text_lf()
    day = []
    _section_edit(t, "PATTERNS", lambda tok, ln: day.extend(tok[1:]) if tok[0] == "1" else None)
    assert len(day) == 24
    week = [f"{float(m) * k:.4g}" for k in WEEK_DAY_FACTORS for m in day]
    t = _section_edit(t, "PATTERNS", lambda tok, ln: None if tok[0] == "1" else ln)
    lines = [" 1               \t" + "\t".join(week[i:i + 6]) for i in range(0, len(week), 6)]
    return _crlf(_insert_after_header(t, "PATTERNS", "\n".join(lines)))


def ky10_const_elev() -> bytes:
    """ky10 with every junction at one elevation, the mean of its junctions' (a skeleton made from a drawing has no
    elevations): ky10 has five PRVs, so this must block."""
    t = _lib("ky10").decode("utf-8").replace("\r\n", "\n")
    elev = []
    _section_edit(t, "JUNCTIONS", lambda tok, ln: elev.append(float(tok[1])) or ln)
    mean = f"{sum(elev) / len(elev):.2f}"
    return _crlf(_section_edit(t, "JUNCTIONS", lambda tok, ln: _set_token(tok, ln, 1, mean)))


def net3_no_coords() -> bytes:
    """Net3 without its [COORDINATES] and [VERTICES] lines."""
    t = _section_edit(_net3_text_lf(), "COORDINATES", lambda tok, ln: None)
    return _crlf(_section_edit(t, "VERTICES", lambda tok, ln: None))


def net3_no_controls() -> bytes:
    """Net3 with its 18 controls removed: pump 335 then runs all the time and pump 10 stays closed."""
    return _crlf(_section_edit(_net3_text_lf(), "CONTROLS", lambda tok, ln: None))


def net3_speed_pattern_pumps() -> bytes:
    """Net3 with its controls removed and each pump given a speed pattern instead (PATTERN in [PUMPS]; pump 10 open
    from 01:00 to 15:00, pump 335 always): a pump with a speed pattern is controlled, and must not be flagged."""
    t = _section_edit(_net3_text_lf(), "CONTROLS", lambda tok, ln: None)
    t = _section_edit(t, "STATUS", lambda tok, ln: None)
    pats = {"10": "SPD10", "335": "SPD335"}
    t = _section_edit(t, "PUMPS", lambda tok, ln: " " + "\t".join(tok + ["PATTERN", pats[tok[0]]]))
    on10 = ["0"] + ["1"] * 14 + ["0"] * 9
    lines = [f" SPD10           \t" + "\t".join(on10[i:i + 6]) for i in range(0, 24, 6)] + [" SPD335          \t1"]
    return _crlf(_insert_after_header(t, "PATTERNS", "\n".join(lines)))


def net3_zero_demand() -> bytes:
    """Net3 with every junction's base demand set to 0."""
    return _crlf(_section_edit(_net3_text_lf(), "JUNCTIONS", lambda tok, ln: _set_token(tok, ln, 2, "0")))


def net3_two_parts() -> bytes:
    """Net3 plus a second part no pipe joins to it: reservoir RX, pipe PX and junction JX (which runs: RX is its
    entry point)."""
    t = _insert_after_header(_net3_text_lf(), "RESERVOIRS", " RX              \t200         \t                \t;")
    t = _insert_after_header(t, "JUNCTIONS", " JX              \t100         \t10          \t                \t;")
    t = _insert_after_header(t, "PIPES", " PX              \tRX              \tJX              \t1000        \t12"
                                         "          \t100         \t0           \tOpen  \t;")
    t = _insert_after_header(t, "COORDINATES", " RX              \t1.00            \t1.00")
    return _crlf(_insert_after_header(t, "COORDINATES", " JX              \t2.00            \t1.00"))


def net3_no_source() -> bytes:
    """Net3 with its reservoirs and tanks written as junctions: no reservoir, no tank and no negative demand, so no
    place where water enters (blocked)."""
    t = _net3_text_lf()
    moved = []
    t = _section_edit(t, "RESERVOIRS", lambda tok, ln: moved.append((tok[0], tok[1])))
    t = _section_edit(t, "TANKS", lambda tok, ln: moved.append((tok[0], tok[1])))
    return _crlf(_insert_after_header(t, "JUNCTIONS", "\n".join(f" {i:<16}\t{e:<12}\t0" for i, e in moved)))


def watergems_export_like() -> bytes:
    """Net3 as a WaterGEMS EPANET export writes a pump: each pump's suction end moved to a new junction joined to the
    old one by a connector pipe of 1 ft and 99 in, tagged AddedForINPConversion in [TAGS]."""
    t = _net3_text_lf()
    pumps = {}
    _section_edit(t, "PUMPS", lambda tok, ln: pumps.__setitem__(tok[0], (tok[1], tok[2])) or ln)
    xy, elev = {}, {}
    _section_edit(t, "COORDINATES", lambda tok, ln: xy.__setitem__(tok[0], (tok[1], tok[2])) or ln)
    _section_edit(t, "JUNCTIONS", lambda tok, ln: elev.__setitem__(tok[0], tok[1]) or ln)
    t = _section_edit(t, "PUMPS", lambda tok, ln: _set_token(tok, ln, 1, f"{tok[0]}_CN"))
    for p, (a, b) in pumps.items():
        t = _insert_after_header(t, "JUNCTIONS", f" {p}_CN{'':<10}\t{elev.get(a, elev[b])}\t0")
        t = _insert_after_header(t, "PIPES", f" {p}_CP{'':<10}\t{a:<16}\t{p}_CN{'':<10}\t1\t99\t130\t0\tOpen\t;")
        t = _insert_after_header(t, "COORDINATES", f" {p}_CN{'':<10}\t{xy[a][0]}\t{xy[a][1]}")
        t = _insert_after_header(t, "TAGS", f" LINK {p}_CP{'':<10}\tAddedForINPConversion")
    return _crlf(t)


# ----------------------------------------------------------------------------- task 15's leftovers, fixed in task 15b
SWMM_JUNCTIONS = 3000                  # junction lines before the first SWMM-only section, past the first 64 KB


def swmm_after_64kb() -> bytes:
    """A SWMM model whose first 64 KB holds only headers EPANET also uses ([TITLE], [OPTIONS], [JUNCTIONS]) and whose
    SWMM-only sections ([OUTFALLS], [CONDUITS], [XSECTIONS]) come after: typed as EPANET by task 15, then blocked
    with WNTR's syntax error instead of the SWMM advice."""
    rows = "".join(f"J{i:<8d} {100 + i % 50:<8.2f} 10 0 0 0\n" for i in range(SWMM_JUNCTIONS))
    return ("[TITLE]\nstorm sewer model\n\n[OPTIONS]\nFLOW_UNITS CFS\nINFILTRATION HORTON\n\n[JUNCTIONS]\n"
            ";Name Elev MaxDepth InitDepth SurDepth Aponded\n" + rows + "\n[OUTFALLS]\nO1 90 FREE NO\n\n[CONDUITS]\n"
            "C1 J0 O1 400 0.01 0 0 0 0\n\n[XSECTIONS]\nC1 CIRCULAR 1 0 0 0 1\n\n[END]\n").encode("ascii")


def zip_model_as_txt() -> bytes:
    """An emailed zip whose one member is Net3 renamed zone_model.txt: typed by task 15 as a zip with nothing
    recognised, because zips were typed by their member names only."""
    return _zip([("zone_model.txt", _net3())])


ACCENTED_RENAMES = (("10", "Válvula_10"), ("15", "Presión_Á"), ("20", "Bóveda_à1"), ("35", "Niño_35"))


def utf8_stray_cp1252() -> bytes:
    """Net3 in UTF-8 with four accented junction IDs (ACCENTED_RENAMES) and one stray cp1252 byte, a degree sign
    pasted into the title: task 15 decoded the whole file as cp1252, so each accent became two characters (the 0xA0
    of 'à' a no-break space WNTR splits on) and WNTR could not read it."""
    t = _net3_text_lf()
    for old, new in ACCENTED_RENAMES:
        t = _rename_node(t, old, new)
    b = _crlf(t.replace("[TITLE]\n", "[TITLE]\nZona norte, 50 °F\n", 1))
    return b.replace("°".encode("utf-8"), b"\xb0")


# ----------------------------------------------------------------------------- task 15b's review cases
DW_UNIFORM_MFT = 0.26 / 0.3048         # 0.26 mm in millifeet: one Darcy-Weisbach roughness on every pipe


def net3_dw_uniform() -> bytes:
    """Net3 in Darcy-Weisbach with one roughness on every pipe (0.26 mm, the plan's example), the Darcy-Weisbach form
    of Net1's one C value (the review's case): the converted copy's C values differ by diameter and velocity only, so
    the old-pipe factor cannot be learned from them and the file must get the fixed-gamma grid."""
    t = _net3_text_lf()
    assert " Headloss           \tH-W" in t
    t = t.replace(" Headloss           \tH-W", " Headloss           \tD-W")
    return _crlf(_section_edit(t, "PIPES", lambda tok, ln: _set_token(tok, ln, 5, f"{DW_UNIFORM_MFT:.6g}")))


def net3_rule_premise() -> bytes:
    """Net3 with its controls removed and one rule that reads pump 10's status and sets pipe 20 (the review's case):
    pump 10 starts closed and nothing opens it, so it must be flagged with pump 335."""
    t = _section_edit(_net3_text_lf(), "CONTROLS", lambda tok, ln: None)
    return _crlf(_insert_after_header(t, "RULES", "RULE 1\nIF PUMP 10 STATUS IS OPEN\nTHEN PIPE 20 STATUS IS OPEN"))


CP1252_PAIR_RENAMES = (("10", "ÉCOLE_10"), ("15", "CITÉ’15"))


def net3_cp1252_pair() -> bytes:
    """Net3 saved in cp1252 throughout (the review's case), with a title and an ID in which an accented capital is
    followed by a curly quote: 'É' then '”' or '’' are the bytes C9 94 and C9 92, valid two-byte UTF-8 by chance, so a
    test that took any valid UTF-8 pair as proof of UTF-8 read the file as UTF-8 with stray bytes and changed the ID.
    It has 4 such pairs and 9 bytes that are not UTF-8: it must decode as cp1252 throughout, with its IDs kept."""
    t = _net3_text_lf().replace("[TITLE]\n", "[TITLE]\nQuartier “CITÉ” près de l’école, 10°C\n", 1)
    for old, new in CP1252_PAIR_RENAMES:
        t = _rename_node(t, old, new)
    return t.replace("\n", "\r\n").encode("cp1252")


# ----------------------------------------------------------------------------- task 15b's second review cases
ONE_PAIR_ID = "CITÉ’15"                                  # the review's ID: É then a curly apostrophe, bytes C9 92
EN_DASH = "\u2013"            # written as an escape: the iteration-5 files hold no dash character (R4)
DASH_IDS = tuple((j, f"CITÉ{EN_DASH}{j}") for j in ("15", "20", "35", "40", "50"))   # É then an en dash, bytes C9 96


def net3_cp1252_one_pair_id() -> bytes:
    """Net3 in cp1252 throughout with one junction renamed CITÉ’15 and a title with one degree sign (the second
    review's file): the chance pair is written wherever the ID is (its junction line, the pipe ends, the coordinates),
    3 times against 1 stray byte, so a rule that counted occurrences read it as UTF-8 and changed the ID."""
    t = _net3_text_lf().replace("[TITLE]\n", "[TITLE]\nZone 10°C\n", 1)
    return _rename_node(t, "15", ONE_PAIR_ID).replace("\n", "\r\n").encode("cp1252")


def net3_cp1252_dash_ids() -> bytes:
    """Net3 in cp1252 throughout with five IDs 'CITÉ', an en dash and the number, and a title 'Réseau 10°C': 18 valid
    UTF-8 pairs by chance against 2 stray bytes (the second review's case A)."""
    t = _net3_text_lf().replace("[TITLE]\n", "[TITLE]\nRéseau 10°C\n", 1)
    for old, new in DASH_IDS:
        t = _rename_node(t, old, new)
    return t.replace("\n", "\r\n").encode("cp1252")


MORE_STRAY_RENAMES = (("15", "Presión_15"), ("20", "Bóveda_20"))
PASTED_CP1252 = "“Zona norte” l’agua 50°F"                 # 4 bytes that are not UTF-8 once saved in cp1252


def utf8_more_stray_than_accents() -> bytes:
    """Net3 in UTF-8 with two accented IDs and a comment pasted in cp1252 with 4 stray bytes, more than the accented
    letters (the second review's case C): it must read as UTF-8 with the 4 bytes as cp1252."""
    t = _net3_text_lf()
    for old, new in MORE_STRAY_RENAMES:
        t = _rename_node(t, old, new)
    b = t.replace("[TITLE]\n", "[TITLE]\n;PASTE\n", 1).replace("\n", "\r\n").encode("utf-8")
    return b.replace(b";PASTE", b";" + PASTED_CP1252.encode("cp1252"))


ONE_CAPITAL_ID = ("35", "CÔTE_35")                         # Ô in UTF-8 is C3 94: 'Ã”' in cp1252, plausible both ways


def utf8_one_id_six_stray() -> bytes:
    """Net3 in UTF-8 with one ID whose only accent is a capital (CÔTE_35) and six stray cp1252 bytes in the title:
    the stray bytes outnumber the ID's occurrences, which read the file as cp1252 throughout and changed the ID
    before the second review (12 of the review's 150 such files did)."""
    t = _rename_node(_net3_text_lf(), *ONE_CAPITAL_ID)
    b = t.replace("[TITLE]\n", "[TITLE]\nZone @@\n", 1).replace("\n", "\r\n").encode("utf-8")
    return b.replace(b"@@", "°’“”°’".encode("cp1252"))


PASTED_BLOCK_RENAMES = (("15", "Bóveda_15"), ("20", "Señal_20"))
PASTED_BLOCK = ("Réseau “Nord”, 10°C", "Château d’eau, 50 °F", "Débit 12 m³/h", f"Presión {EN_DASH} 40 psi", "pH 7,5 ± 0,2",
                "cloración 1,2 mg/L", "Bombeo “Peña”", "Zona de presión norte", "l’eau traitée", "Été 2024")


def utf8_pasted_block() -> bytes:
    """Net3 in UTF-8 with two IDs with a lowercase accent and a block of ten comment lines pasted in cp1252: the stray
    bytes sit on more lines than the IDs, but a lowercase accented letter in UTF-8 ('ó' is C3 B3, 'Ã³' in cp1252) is
    not a pair a cp1252 file makes, so the file reads as UTF-8 with those lines as cp1252."""
    t = _net3_text_lf()
    for old, new in PASTED_BLOCK_RENAMES:
        t = _rename_node(t, old, new)
    b = t.replace("[TITLE]\n", "[TITLE]\n;PASTE\n", 1).replace("\n", "\r\n").encode("utf-8")
    return b.replace(b";PASTE", "\r\n".join(";" + ln for ln in PASTED_BLOCK).encode("cp1252"))


PASTED_PAIR = "Quartier CITÉ’S, 10°C"                     # a chance pair (C9 92) inside a pasted cp1252 comment


def utf8_pasted_chance_pair() -> bytes:
    """Net3 in UTF-8 with the four accented IDs of utf8_stray_cp1252 and a comment pasted in cp1252 that holds a pair
    valid in UTF-8 by chance ('É’'): the file reads as UTF-8, and the pair, like the stray degree sign, is read as
    cp1252, so the comment keeps 'CITÉ’S'."""
    t = _net3_text_lf()
    for old, new in ACCENTED_RENAMES:
        t = _rename_node(t, old, new)
    b = t.replace("[TITLE]\n", "[TITLE]\n;PASTE\n", 1).replace("\n", "\r\n").encode("utf-8")
    return b.replace(b";PASTE", b";" + PASTED_PAIR.encode("cp1252"))


def encoding_fuzz(seed: int = 11, n: int = 150) -> list[tuple[str, bytes, list[str]]]:
    """(kind, file, the IDs it must keep): n files of Net3 in cp1252 throughout, with 1 to 4 junctions renamed to a
    random ID with one accented letter (sometimes followed by a quote, a dash or a degree sign) and a random title,
    then n files in UTF-8 with such IDs and 1 to 6 stray cp1252 bytes in the title.  The second review's generator,
    call for call, so seed 11 gives its 300 files."""
    import random
    letters = "ÀÁÂÃÄÅÇÈÉÊËÌÍÎÏÑÒÓÔÕÖØÙÚÛÜÝàáâãäåçèéêëìíîïñòóôõöøùúûüý"
    punct = "’“”" + EN_DASH + "°"
    base = _net3_text_lf()
    js = ["15", "20", "35", "40", "50", "60", "101", "103", "105", "107"]
    r = random.Random(seed)

    def name(j):
        s = "".join(r.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for _ in range(r.randint(1, 4)))
        pos = r.randrange(len(s) + 1)
        ch = r.choice(letters)
        if r.random() < 0.3:
            ch += r.choice(punct)
        return s[:pos] + ch + s[pos:] + "_" + j
    out = []
    for kind in ("cp1252_throughout", "utf8_with_stray"):
        for _ in range(n):
            t = base
            want = []
            for j in r.sample(js, r.randint(1, 4)):
                nm = name(j)
                t = _rename_node(t, j, nm)
                want.append(nm)
            title = "Zone " + "".join(r.choice(letters + punct + "abc ") for _ in range(r.randint(0, 12)))
            if kind == "cp1252_throughout":
                data = t.replace("[TITLE]\n", "[TITLE]\n" + title + "\n", 1).replace("\n", "\r\n").encode("cp1252")
            else:
                stray = "".join(r.choice("°’“”") for _ in range(r.randint(1, 6)))
                data = (t.replace("[TITLE]\n", "[TITLE]\n@@\n", 1).replace("\n", "\r\n").encode("utf-8")
                        .replace(b"@@", ("Zone " + stray).encode("cp1252")))
            out.append((kind, data, want))
    return out


def _clock_start(clock: str) -> bytes:
    t = _net3_text_lf()
    assert " Start ClockTime    \t12 am" in t
    return _crlf(t.replace(" Start ClockTime    \t12 am", f" Start ClockTime    \t{clock}"))


def net3_clock_235950() -> bytes:
    """Net3 starting at 23:59:50 (START CLOCKTIME): 00:00 to the minute, so no model-hours warning."""
    return _clock_start("23:59:50")


def net3_clock_000020() -> bytes:
    """Net3 starting at 0:00:20: 00:00 to the minute, so no model-hours warning."""
    return _clock_start("0:00:20")


DW_PLASTIC_MM = (0.0015, 0.0024)       # net3_dw_plastic.inp: two Darcy-Weisbach roughness values of plastic pipe


def net3_dw_plastic() -> bytes:
    """Net3 in Darcy-Weisbach with 0.0015 mm on every other pipe and 0.0024 mm on the rest (in millifeet in the
    file): two roughness values about 1e-6 m apart, which a test that rounded metres to 6 decimals read as one."""
    t = _net3_text_lf()
    assert " Headloss           \tH-W" in t
    t = t.replace(" Headloss           \tH-W", " Headloss           \tD-W")
    k = [0]

    def f(tok, ln):
        k[0] += 1
        return _set_token(tok, ln, 5, f"{DW_PLASTIC_MM[k[0] % 2] / 0.3048:.6g}")
    return _crlf(_section_edit(t, "PIPES", f))


def net1_export_like() -> bytes:
    """Net1 (C 100 on every pipe) as watergems_export_like writes an export: pump 9's suction end moved to a new
    junction joined to the old one by a connector of 1 ft and 99 in with C 130, tagged AddedForINPConversion.  The
    connector is not a second C value of the network: the file must still get the fixed-gamma grid."""
    t = _lib("Net1").decode("utf-8").replace("\r\n", "\n")
    pumps, xy, elev = {}, {}, {}
    _section_edit(t, "PUMPS", lambda tok, ln: pumps.__setitem__(tok[0], (tok[1], tok[2])) or ln)
    _section_edit(t, "COORDINATES", lambda tok, ln: xy.__setitem__(tok[0], (tok[1], tok[2])) or ln)
    _section_edit(t, "JUNCTIONS", lambda tok, ln: elev.__setitem__(tok[0], tok[1]) or ln)
    t = _section_edit(t, "PUMPS", lambda tok, ln: _set_token(tok, ln, 1, f"{tok[0]}_CN"))
    if "[TAGS]" not in t:
        t = t.replace("[END]", "[TAGS]\n\n[END]")
    for p, (a, b) in pumps.items():
        t = _insert_after_header(t, "JUNCTIONS", f" {p}_CN\t{elev.get(a, elev.get(b, '0'))}\t0")
        t = _insert_after_header(t, "PIPES", f" {p}_CP\t{a}\t{p}_CN\t1\t99\t130\t0\tOpen\t;")
        t = _insert_after_header(t, "COORDINATES", f" {p}_CN\t{xy[a][0]}\t{xy[a][1]}")
        t = _insert_after_header(t, "TAGS", f" LINK {p}_CP\tAddedForINPConversion")
    return _crlf(t)


# ----------------------------------------------------------------------------- task 16: grab logs and tap maps
# Written here by hand (string formatting and the codecs), never by importers.read_grab_log; grab_log_typed.xlsx is
# written by openpyxl, the library that reads it, so its check is same-library self-consistency.  The readings, taps
# and dates are invented; the junctions are Net3's, as in docs/example_tap_map.csv.
T16_TAPS = (("T1", "225"), ("T2", "239"), ("T3", "193"), ("T4", "173"), ("T5", "127"), ("T6", "35"), ("T7", "15"),
            ("T8", "139"), ("T9", "113"), ("T10", "60"), ("SP-01", "143"), ("SP-02", "261"), ("SP-03", "243"),
            ("SP-04", "61"), ("AN-1", "187"), ("AN-2", "263"), ("007", "219"), ("012", "267"), ("101", "206"),
            ("Église_3", "269"))
# lims_one_datetime_us.csv: (month, day, h, min, tap, total chlorine mg/L); every day is 12 or less and differs from
# its month, so each date reads two ways and the file needs a date format; times are written '9:15 AM'
LIMS_US_ROWS = tuple((mo, d, h, mi, tap, round(1.05 + 0.04 * k + 0.03 * (mo - 6), 2))
                     for mo, days in ((6, (3, 9)), (7, (2, 11)), (8, (4, 10)), (9, (1, 8)), (10, (5, 12)))
                     for d in days
                     for k, (h, mi, tap) in enumerate(((9, 15, "SP-01"), (10, 40, "SP-02"), (13, 5, "SP-03"),
                                                       (14, 30, "SP-04"))))
# lims_censored.csv: (date, time, tap, the reading as written); the expected outcome of each row
LIMS_CENSORED_ROWS = (("2026-03-02", "09:15", "T1", "0.45"), ("2026-03-02", "10:40", "T2", "<0.02"),
                      ("2026-03-03", "08:55", "T3", "ND"), ("2026-03-03", "11:20", "T4", ""),
                      ("2026-03-04", "09:05", "T5", "< 0.05"), ("2026-03-04", "13:30", "T6", "0,31"),
                      ("2026-03-05", "14:10", "T7", "BDL"), ("2026-03-05", "15:45", "T8", "0.27 mg/L"))
LIMS_CENSORED_EXPECT = {"used": ("T1", "T6", "T8"), "censored": ("T2", "T3", "T5", "T7"), "blank": ("T4",)}
AMBIGUOUS_ROWS = (("03/04/2026", "09:15", "T1", "0.45"), ("05/06/2026", "10:40", "T2", "0.38"),
                  ("07/08/2026", "08:55", "T3", "0.52"))
CP1252_ROWS = (("2026-03-02", "09:15", "T1", "0.45", "résidu normal"), ("2026-03-02", "10:40", "Église_3", "0.38",
                                                                            "près de l'école"),
               ("2026-03-03", "08:55", "T3", "0.52", "°C noté"))
# grab_log_typed.xlsx and its CSV twin: (date, h, min, tap as the tap map writes it, reading)
TYPED_ROWS = (("2026-04-06", 9, 10, "007", 0.62), ("2026-04-06", 11, 45, "012", 0.41), ("2026-04-07", 8, 30, "101", 0.55),
              ("2026-04-07", 14, 20, "T2", 0.33), ("2026-04-08", 10, 5, "007", 0.59))
ANALYZER_SITES = ("AN-1", "AN-2")
ANALYZER_DAYS = 30
ANALYZER_STEP_MIN = 15
UTC_SHIFT_H = 8                         # utc_shifted_log.csv: local times 10:00 to 15:00 written 8 h later
UTC_LOCAL = ((10, 5), (11, 20), (12, 40), (13, 15), (14, 30), (15, 0))
MULTI_MONTHS = 9                        # multi_month_log.csv: January to September 2026, 4 readings a month
# the clock conversion logs: (date, model hour, minute, tap, reading); written once with the time of day equal to the
# model hour (for Net3, which starts at midnight) and once 6 h later (for net3_clock6.inp, which starts at 06:00)
# over five months (March to July 2026), so the pilot's hold-out has months to hold out and the app's window takes three
CLOCK_ROWS = tuple((f"2026-{m:02d}-{d + 2 * k:02d}", h, mi, tap, round(v + 0.01 * (m - 3), 2))
                   for m in range(3, 8)
                   for k, (d, h, mi, tap, v) in enumerate(((4, 3, 0, "T1", 0.41), (4, 4, 30, "T2", 0.36),
                                                          (5, 5, 15, "T3", 0.58), (5, 7, 40, "T4", 0.29),
                                                          (6, 9, 5, "T5", 0.47), (6, 10, 20, "T6", 0.52))))
CLOCK_SHIFT_H = 6


def _csv_text(header: str, rows, sep: str = ",") -> str:
    def q(c):
        c = str(c)
        return f'"{c}"' if (sep in c or '"' in c) else c
    return header + "\n" + "".join(sep.join(q(c) for c in r) + "\n" for r in rows)


def task16_tap_map() -> bytes:
    return _csv_text("tap_id,junction_id,notes", [(t, j, "invented") for t, j in T16_TAPS]).encode("utf-8")


def lims_one_datetime_us() -> bytes:
    """A LIMS export: a sample ID per row, 'Sample Point', one 'Sample Date/Time' column with US dates and AM/PM times,
    'Total Chlorine (mg/L)' and a units column, over five months (enough for the pilot's hold-out)."""
    rows = []
    for k, (mo, d, h, mi, tap, v) in enumerate(LIMS_US_ROWS, start=1):
        h12, ap = (h - 12, "PM") if h > 12 else (h, "PM" if h == 12 else "AM")
        rows.append((f"S{k:04d}", tap, f"{mo:02d}/{d:02d}/2026 {h12}:{mi:02d} {ap}", f"{v:.2f}", "mg/L", "lab"))
    return _csv_text("Sample ID,Sample Point,Sample Date/Time,Total Chlorine (mg/L),Units,Analyst", rows).encode()


def lims_censored() -> bytes:
    return _csv_text("Sample Point,Sample Date,Sample Time,Free Chlorine (mg/L),Qualifier",
                     [(t, d, h, v, "") for d, h, t, v in LIMS_CENSORED_ROWS]).encode()


def ambiguous_dates() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL", AMBIGUOUS_ROWS).encode()


def grab_log_cp1252() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL,notes", CP1252_ROWS).replace("\n", "\r\n").encode("cp1252")


def grab_log_bom() -> bytes:
    return b"\xef\xbb\xbf" + _csv_text("date,time,tap_id,free_chlorine_mgL", [r[:4] for r in CP1252_ROWS
                                                                             if r[2] != "Église_3"]).encode()


def grab_log_typed_xlsx() -> bytes:
    """Typed date and time cells, and the numeric tap IDs as numbers (007 is 7), as Excel keeps them: written by
    openpyxl (same-library self-consistency with the reader)."""
    import datetime as _dt
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "samples"
    ws.append(["date", "time", "tap_id", "free_chlorine_mgL"])
    for d, h, mi, tap, v in TYPED_ROWS:
        ws.append([_dt.date.fromisoformat(d), _dt.time(h, mi), int(tap) if tap.isdigit() else tap, v])
        ws.cell(row=ws.max_row, column=1).number_format = "yyyy-mm-dd"
        ws.cell(row=ws.max_row, column=2).number_format = "hh:mm"
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def grab_log_typed_csv() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL",
                     [(d, f"{h:02d}:{mi:02d}", tap, f"{v:.2f}") for d, h, mi, tap, v in TYPED_ROWS]).encode()


def analyzer_15min() -> bytes:
    """Two online analyzers (AN-1, AN-2) at 15-minute steps for 30 days of April 2026, and three grab samples."""
    rows = []
    for day in range(1, ANALYZER_DAYS + 1):
        for k, site in enumerate(ANALYZER_SITES):
            for step in range(0, 24 * 60, ANALYZER_STEP_MIN):
                v = 0.55 + 0.1 * k + 0.08 * math.sin(2 * math.pi * step / 1440) - 0.002 * day
                rows.append((f"2026-04-{day:02d}", f"{step // 60:02d}:{step % 60:02d}", site, f"{v:.3f}"))
    rows += [("2026-04-02", "09:30", "T1", "0.44"), ("2026-04-09", "11:10", "T2", "0.37"),
             ("2026-04-16", "14:25", "T3", "0.51")]
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def utc_shifted_log() -> bytes:
    rows = [(f"2026-03-{2 + i:02d}", f"{(h + UTC_SHIFT_H) % 24:02d}:{mi:02d}", f"T{i + 1}", f"{0.4 + 0.02 * i:.2f}")
            for i, (h, mi) in enumerate(UTC_LOCAL)]
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def unknown_tap_log() -> bytes:
    rows = [("2026-03-02", "09:15", "T1", "0.45"), ("2026-03-02", "10:40", "T99", "0.38"),
            ("2026-03-03", "08:55", "T3", "0.52"), ("2026-03-03", "11:20", "T4", "0.21")]
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def mixed_species_log() -> bytes:
    rows = [("2026-03-02", "09:15", "T1", "0.45", "1.10"), ("2026-03-03", "10:40", "T2", "0.38", "0.95")]
    return _csv_text("date,time,tap_id,free_chlorine_mgL,total_chlorine_mgL", rows).encode()


def multi_month_log() -> bytes:
    rows = []
    for m in range(1, MULTI_MONTHS + 1):
        for k in range(4):
            rows.append((f"2026-{m:02d}-{3 + 6 * k:02d}", f"{8 + 2 * k:02d}:{10 + 5 * k:02d}", f"T{1 + (m + k) % 10}",
                         f"{0.35 + 0.03 * k + 0.01 * m:.2f}"))
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def _clock_log(shift_h: int) -> bytes:
    rows = []
    for d, h, mi, tap, v in CLOCK_ROWS:
        rows.append((d, f"{(h + shift_h) % 24:02d}:{mi:02d}", tap, f"{v:.2f}"))
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def clock_log_net3() -> bytes:
    return _clock_log(0)


def clock_log_clock6() -> bytes:
    return _clock_log(CLOCK_SHIFT_H)


# ----------------------------------------------------------------------------- task 16's review: logs as they arrive
# Each case the review of task 16 found, written by hand (the workbooks by openpyxl: same-library self-consistency).
# (date, time, tap, reading as written) rows; the junctions are Net3's through task16_tap_map.csv.
REVIEW_ROWS = (("2026-03-10", "09:15", "T1", "0.45"), ("2026-03-11", "10:20", "T2", "0.38"),
               ("2026-03-12", "11:05", "T3", "0.52"), ("2026-03-13", "13:40", "T4", "0.29"))
OVERRANGE_ROWS = (("2026-03-10", "09:15", "T1", ">2.2"), ("2026-03-10", "10:20", "T3", "<0.02"),
                  ("2026-03-11", "11:05", "T2", "0.41"), ("2026-03-11", "13:40", "T4", "> 4.0"))
IMPLAUSIBLE_ROWS = (("2026-03-10", "09:15", "T1", "9999"), ("2026-03-10", "10:20", "T2", "0.38"),
                    ("2026-03-11", "11:05", "T3", "450"), ("2026-03-11", "13:40", "T4", "0.29"))
# lims_long.csv: one row per analyte, as a LIMS exports it; each sample (date, time, tap, free, total) gets a free and a
# total chlorine row, a pH row and a temperature row
LONG_SAMPLES = (("2026-03-10", "09:15", "SP-01", "0.45", "1.20"), ("2026-03-11", "10:40", "SP-02", "0.38", "1.05"),
                ("2026-03-12", "13:05", "SP-03", "0.52", "1.31"))
ANALYZER_MIDNIGHT_DAYS = 2              # analyzer_midnight.xlsx: one site, 15-minute data, typed datetimes from 00:00
STUDY_HOURLY_DAYS = ("2026-03-04", "2026-03-05")    # grab_study_24h.csv: T1 hourly on two days, T3 two-hourly on one
INTERLEAVED_SITES = tuple(f"T{i}" for i in range(1, 11))   # analyzer_interleaved.csv: 10 analyzers, by timestamp
INTERLEAVED_DAYS = 2
STRAY_CELL = (20000, 2000)              # log_cover_and_stray.xlsx: one stray cell far down and far right


def _review_log(rows, header: str = "date,time,tap_id,free_chlorine_mgL", sep: str = ",") -> str:
    return _csv_text(header.replace(",", sep), rows, sep)


def junction_list() -> bytes:
    """A GIS or EPANET junction export: a junction column, no tap column, no reading (crashed 'Check a file')."""
    return _csv_text("Node ID,Elevation,Base Demand", [(j, 10 + k, 1.5) for k, (_, j) in enumerate(T16_TAPS[:6])]).encode()


def log_title_line() -> bytes:
    return ("Residuals, March 2026\n" + _review_log(REVIEW_ROWS)).encode()


def log_lab_preamble() -> bytes:
    return ("Lab: Distribution\nExported 2026-04-01 by the LIMS\n" + _review_log(REVIEW_ROWS)).encode()


def log_sep_hint() -> bytes:
    """Excel's 'sep=;' hint line, semicolons and decimal commas, as a European export writes them."""
    rows = [(d, t, tap, v.replace(".", ",")) for d, t, tap, v in REVIEW_ROWS]
    return ("sep=;\n" + _review_log(rows, sep=";")).encode()


def log_chloride_first() -> bytes:
    """A chloride column before the chlorine one, in a LIMS's alphabetical order."""
    rows = [(d, t, tap, f"{35 + k}", v) for k, (d, t, tap, v) in enumerate(REVIEW_ROWS)]
    return _csv_text("date,time,tap_id,Chloride (mg/L),Chlorine Residual (mg/L)", rows).encode()


def log_month_names() -> bytes:
    """Dates with the month's name ('10-Mar-26', '11 March 2026', 'Mar 12, 2026') and times with dots ('09.15')."""
    def mon(d, k):
        y, m, dd = d.split("-")
        name = ("Mar", "March")[k % 2]
        return (f"{int(dd)}-{name[:3]}-{y[2:]}", f"{int(dd)} {name} {y}", f"{name[:3]} {int(dd)}, {y}")[k % 3]
    rows = [(mon(d, k), t.replace(":", "."), tap, v) for k, (d, t, tap, v) in enumerate(REVIEW_ROWS)]
    return _review_log(rows).encode()


def log_multiline_comment() -> bytes:
    """A quoted comment that runs over two lines, then a tap the tap map does not have, on line 4 of the file."""
    return ('date,time,tap_id,free_chlorine_mgL,notes\n2026-03-10,09:15,T1,0.45,"line one\nline two"\n'
            '2026-03-11,10:15,T99,0.55,\n').encode()


def log_overrange() -> bytes:
    return _review_log(OVERRANGE_ROWS).encode()


def log_ugL() -> bytes:
    """The readings in micrograms per litre, as the header says."""
    rows = [(d, t, tap, f"{round(float(v) * 1000):d}") for d, t, tap, v in REVIEW_ROWS]
    return _csv_text("date,time,tap_id,Free Chlorine (ug/L)", rows).encode()


def log_implausible() -> bytes:
    """A LIMS missing code (9999) and a reading in the wrong unit (450) among two real readings."""
    return _review_log(IMPLAUSIBLE_ROWS).encode()


def lims_long() -> bytes:
    rows = []
    for d, t, tap, free, total in LONG_SAMPLES:
        rows += [(tap, d, t, "Free Chlorine", free, "mg/L"), (tap, d, t, "Total Chlorine", total, "mg/L"),
                 (tap, d, t, "pH", "7.6", "SU"), (tap, d, t, "Temperature", "15.2", "C")]
    return _csv_text("Sample Point,Sample Date,Sample Time,Analyte,Result,Units", rows).encode()


def _analyzer_midnight_rows():
    import datetime as _dt
    t0 = _dt.datetime(2026, 3, 1)
    return [(t0 + _dt.timedelta(minutes=15 * k), "AN-1", round(0.5 + 0.001 * (k % 96), 3))
            for k in range(96 * ANALYZER_MIDNIGHT_DAYS)]


def analyzer_midnight_xlsx() -> bytes:
    """One analyzer's 15-minute series with typed datetimes, so each day starts with a cell at exactly 00:00."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["timestamp", "site", "free_chlorine_mgL"])
    for ts, site, v in _analyzer_midnight_rows():
        ws.append([ts, site, v])
        ws.cell(row=ws.max_row, column=1).number_format = "yyyy-mm-dd hh:mm"
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def analyzer_midnight_csv() -> bytes:
    return _csv_text("timestamp,site,free_chlorine_mgL",
                     [(f"{ts:%Y-%m-%d %H:%M}", site, f"{v:.3f}") for ts, site, v in _analyzer_midnight_rows()]).encode()


def grab_study_24h() -> bytes:
    """A diurnal grab-sampling study: T1 every hour on two days, T3 every two hours on one, and three route samples;
    every reading is a grab sample and is kept (it was thinned as an analyzer before the review)."""
    rows = [(d, f"{h:02d}:10", "T1", f"{0.30 + 0.01 * h:.2f}") for d in STUDY_HOURLY_DAYS for h in range(24)]
    rows += [(STUDY_HOURLY_DAYS[0], f"{h:02d}:40", "T3", f"{0.50 - 0.01 * h:.2f}") for h in range(0, 24, 2)]
    rows += [("2026-03-11", "09:30", "T2", "0.44"), ("2026-03-18", "10:45", "T4", "0.31"),
             ("2026-03-25", "14:05", "T5", "0.27")]
    return _review_log(rows).encode()


def analyzer_interleaved() -> bytes:
    """Ten analyzers at 15-minute steps for two days, written row by timestamp as a historian exports them, so each
    site's thinned rows are not consecutive."""
    rows = []
    for day in range(1, INTERLEAVED_DAYS + 1):
        for step in range(0, 24 * 60, 15):
            for k, site in enumerate(INTERLEAVED_SITES):
                rows.append((f"2026-03-{day:02d}", f"{step // 60:02d}:{step % 60:02d}", site, f"{0.4 + 0.01 * k:.2f}"))
    return _review_log(rows).encode()


def log_cover_and_stray_xlsx() -> bytes:
    """A cover sheet first, then the data sheet, which holds one stray cell far down and far right (pandas padded
    such a sheet to rows x columns cells)."""
    import datetime as _dt
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Cover"
    ws.append(["Monthly residuals report"])
    ws.append(["Prepared by the lab"])
    ws2 = wb.create_sheet("Data")
    ws2.append(["date", "time", "tap_id", "free_chlorine_mgL"])
    for d, t, tap, v in REVIEW_ROWS[:2]:
        ws2.append([_dt.date.fromisoformat(d), _dt.time.fromisoformat(t), tap, float(v)])
    ws2.cell(row=STRAY_CELL[0], column=STRAY_CELL[1], value="x")
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()

# ----------------------------------------------------------------------------- task 16's final verification
# Each case the final verification of task 16 found, written by hand (the workbooks by openpyxl: same-library
# self-consistency).  The long-format logs reuse LONG_SAMPLES (free and total chlorine per sample) with the analyte
# names as the verification wrote them; the wide logs reuse REVIEW_ROWS, total chlorine 0.80 mg/L above free.
UGL_HEADERS = ("Free Chlorine ugL", "free_chlorine_ugl", "Free Chlorine (Âµg/L)", "Free Chlorine (ug per L)",
               "Free Chlorine µg L-1")
UGL_ROWS = (("2026-03-10", "09:15", "T1", "450"), ("2026-03-11", "10:20", "T2", "380"),
            ("2026-03-12", "11:05", "T3", "8"), ("2026-03-13", "13:40", "T4", "5"))
UNIT_SPELLINGS = ("mg Cl2/L", "mgCl2/L", "mg/litre")
QUALIFIER_SITES = tuple(f"T{i}" for i in range(1, 11))   # qualifier_interleaved.csv: 10 analyzers, 2 days, '0.452 J'
QUALIFIER_DAYS = 2
UNKNOWN_SITES = 200                                       # unknown_taps_many.csv: 200 taps from another system
TINY_ROWS = 3_000_000                                     # csv_tiny_rows.csv: 'a,b' and 3 million rows '1,'
NET2_LOG_JUNCTIONS = ("3", "7", "11", "15", "19", "23", "27")   # net2_clock_log.csv: 09:00 to 15:00, one a day


def _long_log(names: tuple, header: str = "Analyte", units: bool = True) -> bytes:
    rows = []
    for d, t, tap, free, total in LONG_SAMPLES:
        for nm, v, u in ((names[0], free, "mg/L"), (names[1], total, "mg/L"), ("pH", "7.6", "SU")):
            rows.append((tap, d, t, nm, v) + ((u,) if units else ()))
    return _csv_text(f"Sample Point,Sample Date,Sample Time,{header},Result" + (",Units" if units else ""), rows).encode()


def long_paren() -> bytes:
    """One row per analyte, the species in brackets: 'Chlorine (Free)', 'Chlorine (Total)', 'pH'."""
    return _long_log(("Chlorine (Free)", "Chlorine (Total)"))


def long_short_names() -> bytes:
    """'Free Cl', 'Total Cl' and 'pH', with no units column."""
    return _long_log(("Free Cl", "Total Cl"), units=False)


def long_lims_codes() -> bytes:
    """LIMS codes in a 'Parameter Code' column: FCL2, TCL2, pH."""
    return _long_log(("FCL2", "TCL2"), header="Parameter Code")


def long_unknown_codes() -> bytes:
    """Codes the reader does not know (CHLF, CHLT) under 'Analyte': the analyte value must be chosen."""
    return _long_log(("CHLF", "CHLT"))


def long_unknown_codes_param() -> bytes:
    """The same codes under 'Parameter Code', a header the analyte guess passes over."""
    return _long_log(("CHLF", "CHLT"), header="Parameter Code")


def long_two_unnamed() -> bytes:
    """Two chlorine analytes that name no species: 'Chlorine Residual' and 'Chlorine'."""
    return _long_log(("Chlorine Residual", "Chlorine"))


def _wide(h1: str, h2: str, first_free: bool = True) -> bytes:
    rows = [(d, t, tap, *((v, f"{float(v) + 0.8:.2f}") if first_free else (f"{float(v) + 0.8:.2f}", v)))
            for d, t, tap, v in REVIEW_ROWS]
    return _csv_text(f"date,time,tap_id,{h1},{h2}", rows).encode()


def wide_paren() -> bytes:
    return _wide("Chlorine (Free) mg/L", "Chlorine (Total) mg/L")


def wide_paren_total_first() -> bytes:
    return _wide("Chlorine (Total) mg/L", "Chlorine (Free) mg/L", first_free=False)


def wide_residual_and_free() -> bytes:
    """'Chlorine Residual' (no species) beside 'Free Chlorine' (here the residual column holds the total readings)."""
    return _wide("Free Chlorine", "Chlorine Residual")


def chloramine_header() -> bytes:
    """A 'Chloramine (mg/L)' column, total chlorine as a chloraminated system measures it."""
    return _csv_text("date,time,tap_id,Chloramine (mg/L)", [(d, t, tap, f"{float(v) + 0.8:.2f}")
                                                            for d, t, tap, v in REVIEW_ROWS]).encode()


def ugl_header_log(header: str) -> bytes:
    """UGL_ROWS (µg/L) under one of the UGL_HEADERS spellings, UTF-8."""
    return _csv_text(f"date,time,tap_id,{header}", UGL_ROWS).encode("utf-8")


def units_column_log(unit: str) -> bytes:
    """REVIEW_ROWS with a units column spelled one of the UNIT_SPELLINGS ways."""
    return _csv_text("date,time,tap_id,free_chlorine_mgL,units", [r + (unit,) for r in REVIEW_ROWS]).encode()


def ugl_unnamed() -> bytes:
    """UGL_ROWS under free_chlorine_mgL: a µg/L log whose header does not say so."""
    return _review_log(UGL_ROWS).encode()


def log_comma_title_semicolon() -> bytes:
    """A title with a comma before a semicolon export with decimal commas."""
    rows = [(d, t, tap, v.replace(".", ",")) for d, t, tap, v in REVIEW_ROWS]
    return ("Residuals report, March 2026\n" + _review_log(rows, sep=";")).encode()


def log_two_cell_title_xlsx() -> bytes:
    """A workbook whose first row holds a title and a month, then a blank row, then the header."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["Residuals report", "March 2026"])
    ws.append([])
    ws.append(["date", "time", "tap_id", "free_chlorine_mgL"])
    for d, t, tap, v in REVIEW_ROWS:
        ws.append([d, t, tap, float(v)])
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


DATE_ONLY_ROWS = (((2026, 3, 10, 9, 15), "T1", 0.45), ((2026, 3, 11, 10, 20), "T2", 0.38), ((2026, 3, 12), "T3", 0.21))


def excel_date_only_xlsx() -> bytes:
    """Typed datetimes, one of them a date with no time (a time that was not typed), in one column."""
    import datetime as _dt
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["Sample DateTime", "tap_id", "free_chlorine_mgL"])
    for ts, tap, v in DATE_ONLY_ROWS:
        ws.append([_dt.datetime(*ts), tap, v])
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def excel_date_only_csv() -> bytes:
    return _csv_text("Sample DateTime,tap_id,free_chlorine_mgL",
                     [(f"{ts[0]}-{ts[1]:02d}-{ts[2]:02d}" + (f" {ts[3]:02d}:{ts[4]:02d}" if len(ts) > 3 else ""), tap,
                       f"{v:.2f}") for ts, tap, v in DATE_ONLY_ROWS]).encode()


def time_column_epoch_date() -> bytes:
    """A time column as an Access or SQL time field exports it: '1899-12-30 09:15:00'."""
    return _review_log([(d, f"1899-12-30 {t}:00", tap, v) for d, t, tap, v in REVIEW_ROWS]).encode()


def superscript_tap_log() -> bytes:
    """A tap ID '1²' (a digit and a superscript two: isdigit() is true, int() fails), UTF-8."""
    return _review_log([REVIEW_ROWS[0], ("2026-03-11", "10:20", "1²", "0.38")]).encode("utf-8")


def superscript_tap_map() -> bytes:
    return _csv_text("tap_id,junction_id", [("T1", "225"), ("1²", "239"), ("007", "219")]).encode("utf-8")


def junctions_2col() -> bytes:
    """A two-column junction export: 'Node ID,Elevation'."""
    return _csv_text("Node ID,Elevation", [(j, 10 + k) for k, (_, j) in enumerate(T16_TAPS[:4])]).encode()


def qualifier_interleaved() -> bytes:
    """Ten analyzers at 15-minute steps for two days, by timestamp, every reading with a lab qualifier ('0.452 J')."""
    rows = []
    for day in range(1, QUALIFIER_DAYS + 1):
        for step in range(0, 24 * 60, 15):
            for k, site in enumerate(QUALIFIER_SITES):
                rows.append((f"2026-03-{day:02d}", f"{step // 60:02d}:{step % 60:02d}", site,
                             f"{0.3 + 0.001 * ((step // 15 + 7 * k) % 400):.3f} J"))
    return _review_log(rows).encode()


def unknown_taps_many() -> bytes:
    """200 sampling points, none of them in the task-16 tap map (a tap map from another system)."""
    return _review_log([("2026-03-10", f"10:{s % 60:02d}", f"SITE-{s:05d}", "0.40") for s in range(UNKNOWN_SITES)]).encode()


def csv_tiny_rows() -> bytes:
    """'a,b' and three million rows '1,': 9 MB, far fewer bytes than LOG_MAX_BYTES and far more rows than the cap."""
    return b"a,b\n" + b"1,\n" * TINY_ROWS


def net2_clock_log() -> bytes:
    """Seven samples at 09:00 to 15:00 by the clock, for the bundled Net2 (which starts at 08:00)."""
    return _review_log([(f"2026-03-{10 + k:02d}", f"{9 + k:02d}:00", f"N{k + 1}", f"0.{30 + k}")
                        for k in range(len(NET2_LOG_JUNCTIONS))]).encode()


def net2_tap_map() -> bytes:
    return _csv_text("tap_id,junction_id", [(f"N{k + 1}", j) for k, j in enumerate(NET2_LOG_JUNCTIONS)]).encode()


LONG_UPLOAD_NAME = "Z" * 246 + ".inp"  # a 250-character upload name: legal on APFS, too long once the digest is added


FIXTURES = {
    # the 36 of the sniff table
    "net3_as_txt.txt": net3_as_txt, "net3_backup.inp.bak": net3_backup_inp_bak, "inp_with_bom.inp": inp_with_bom,
    "inp_comment_first.inp": inp_comment_first, "inp_long_title.inp": inp_long_title, "inp_utf16.inp": inp_utf16,
    "net3_cp1252_title.inp": net3_cp1252_title, "emailed_model.zip": emailed_model_zip,
    "dxf_with_999_comment.dxf": dxf_with_999_comment, "ascii_drawing.dxf.bak": ascii_drawing_dxf_bak,
    "binary.dxf": binary_dxf, "dwg_ac1027.dwg": dwg_ac1027, "dwg_ac1032.dwg": dwg_ac1032, "drawing.bak": drawing_bak,
    "dgn_v7.dgn": dgn_v7, "ole_dgn_v8.dgn": ole_dgn_v8,
    "jet_model.wtg.mdb": jet_model_wtg_mdb, "model.wtg.sqlite": model_wtg_sqlite, "wtg_settings.wtg": wtg_settings,
    "mains.gpkg": mains_gpkg, "mains_geojson_as.json": mains_geojson_as_json,
    "geojson_stateplane.geojson": geojson_stateplane, "esri_features.json": esri_features,
    "wntr_sidecar.json": wntr_sidecar, "gml_mains.gml": gml_mains, "mains.kml": mains_kml, "mains.kmz": mains_kmz,
    "mains_shp_nested_two_layers.zip": mains_shp_nested_two_layers_zip, "utility_gdb.zip": utility_gdb_zip,
    "grab_log.xlsx": grab_log_xlsx, "old_excel_log.xls": old_excel_log_xls, "log_semicolon.csv": log_semicolon,
    "log_tab.txt": log_tab, "map_book.pdf": map_book, "plant.ifc": plant_ifc, "results.bin": results_bin,
    # the tie-break cases
    "ole_header_only.xls": ole_header_only, "ole_header_only.dgn": ole_header_only,
    # the other task-15 checks
    "net3_crlf.inp": net3_crlf, "net3_lf.inp": net3_lf, "Net3.inp": net3_copy,
    "net3_dup_junction.inp": net3_dup_junction, "net3_junction_tank_collision.inp": net3_junction_tank_collision,
    "net3_dup_pipe.inp": net3_dup_pipe, "net3_long_id_32.inp": net3_long_id_32,
    "net3_long_id_31.inp": net3_long_id_31,
    "zip_slip.zip": zip_slip, "zip_absolute.zip": zip_absolute, "zip_ratio_bomb.zip": zip_ratio_bomb,
    "zip_two_models.zip": zip_two_models, "shp_padded_dbf.zip": shp_padded_dbf_zip,
    # the review's cases (task 15)
    "swmm_sewer.inp": swmm_sewer_inp, "chlorine_decay.msx": chlorine_decay_msx, "net3_cr_only.inp": net3_cr_only,
    "net3_cr_cr_lf.inp": net3_cr_cr_lf, "inp_title_over_64kb.inp": inp_title_over_64kb,
    "inp_patterns_first.inp": inp_patterns_first, "inp_comment_header_70kb.inp": inp_comment_header_70kb,
    "finder_model.zip": finder_model_zip, "finder_shapefile.zip": finder_shapefile_zip,
    "city_mains_shp.zip": city_mains_shp_zip, "model_and_gis.zip": model_and_gis_zip,
    "geojson_crs_string.geojson": geojson_crs_string, "geojson_crs_name_number.geojson": geojson_crs_name_number,
    "geojson_features_not_list.geojson": geojson_features_not_list, "esri_sr_string.json": esri_sr_string,
    "geojson_deep_geometry.geojson": geojson_deep_geometry, "json_deep.json": json_deep,
    "kml_doctype_subset.kml": kml_doctype_subset, "xml_comment_run.xml": xml_comment_run,
    "xml_comment_run_long.xml": xml_comment_run_long, "xml_pi_run.xml": xml_pi_run, "slow_view.gpkg": slow_view_gpkg,
    "csv_quotes.csv": csv_quotes, "net3_cp1252_undefined_byte.inp": net3_cp1252_undefined_byte,
    "net1_accented_id_cp1252.inp": net1_accented_id_cp1252, "net3_long_pattern_id.inp": net3_long_pattern_id,
    "net3_bom_missing_node.inp": net3_bom_missing_node,
    # the final verification's cases (task 15)
    "net3_singular_headers.inp": net3_singular_headers, "net3_form_feed_headers.inp": net3_form_feed_headers,
    "net3_singular_pipe_dup.inp": net3_singular_pipe_dup, "net3_isolated_junction.inp": net3_isolated_junction,
    "zip_bad_offset.zip": zip_bad_offset, "csv_quote_run.csv": csv_quote_run, "net1_nul_padded.inp": net1_nul_padded,
    # task 15b: physics settings, and task 15's leftovers
    "net3_bulk_lines.inp": net3_bulk_lines, "net3_dw.inp": net3_dw, "net3_cm.inp": net3_cm,
    "net3_clock6.inp": net3_clock6, "net3_weekly_pattern.inp": net3_weekly_pattern,
    "ky10_const_elev.inp": ky10_const_elev, "net3_no_coords.inp": net3_no_coords,
    "net3_no_controls.inp": net3_no_controls, "net3_speed_pattern_pumps.inp": net3_speed_pattern_pumps,
    "net3_zero_demand.inp": net3_zero_demand, "net3_two_parts.inp": net3_two_parts,
    "net3_no_source.inp": net3_no_source, "watergems_export_like.inp": watergems_export_like,
    "swmm_after_64kb.inp": swmm_after_64kb, "zip_model_as_txt.zip": zip_model_as_txt,
    "utf8_stray_cp1252.inp": utf8_stray_cp1252,
    # task 15b's review cases
    "net3_dw_uniform.inp": net3_dw_uniform, "net3_rule_premise.inp": net3_rule_premise,
    "net3_cp1252_pair.inp": net3_cp1252_pair,
    # task 15b's second review cases
    "net3_cp1252_one_pair_id.inp": net3_cp1252_one_pair_id, "net3_cp1252_dash_ids.inp": net3_cp1252_dash_ids,
    "utf8_more_stray_than_accents.inp": utf8_more_stray_than_accents, "utf8_one_id_six_stray.inp": utf8_one_id_six_stray,
    "utf8_pasted_block.inp": utf8_pasted_block, "utf8_pasted_chance_pair.inp": utf8_pasted_chance_pair,
    "net3_clock_235950.inp": net3_clock_235950, "net3_clock_000020.inp": net3_clock_000020,
    "net3_dw_plastic.inp": net3_dw_plastic, "net1_export_like.inp": net1_export_like,
    # task 16: grab logs and tap maps
    "task16_tap_map.csv": task16_tap_map, "lims_one_datetime_us.csv": lims_one_datetime_us,
    "lims_censored.csv": lims_censored, "ambiguous_dates.csv": ambiguous_dates, "grab_log_cp1252.csv": grab_log_cp1252,
    "grab_log_bom.csv": grab_log_bom, "grab_log_typed.xlsx": grab_log_typed_xlsx, "grab_log_typed.csv": grab_log_typed_csv,
    "analyzer_15min.csv": analyzer_15min, "utc_shifted_log.csv": utc_shifted_log, "unknown_tap_log.csv": unknown_tap_log,
    "mixed_species_log.csv": mixed_species_log, "multi_month_log.csv": multi_month_log,
    "clock_log_net3.csv": clock_log_net3, "clock_log_clock6.csv": clock_log_clock6,
    # task 16's review
    "junction_list.csv": junction_list, "log_title_line.csv": log_title_line, "log_lab_preamble.csv": log_lab_preamble,
    "log_sep_hint.csv": log_sep_hint, "log_chloride_first.csv": log_chloride_first,
    "log_month_names.csv": log_month_names, "log_multiline_comment.csv": log_multiline_comment,
    "log_overrange.csv": log_overrange, "log_ugL.csv": log_ugL, "log_implausible.csv": log_implausible,
    "lims_long.csv": lims_long, "analyzer_midnight.xlsx": analyzer_midnight_xlsx,
    "analyzer_midnight.csv": analyzer_midnight_csv, "grab_study_24h.csv": grab_study_24h,
    "analyzer_interleaved.csv": analyzer_interleaved, "log_cover_and_stray.xlsx": log_cover_and_stray_xlsx,
    # task 16's final verification
    "long_paren.csv": long_paren, "long_short_names.csv": long_short_names, "long_lims_codes.csv": long_lims_codes,
    "long_unknown_codes.csv": long_unknown_codes, "long_unknown_codes_param.csv": long_unknown_codes_param,
    "long_two_unnamed.csv": long_two_unnamed, "wide_paren.csv": wide_paren,
    "wide_paren_total_first.csv": wide_paren_total_first, "wide_residual_and_free.csv": wide_residual_and_free,
    "chloramine_header.csv": chloramine_header, "ugl_unnamed.csv": ugl_unnamed,
    "log_comma_title_semicolon.csv": log_comma_title_semicolon, "log_two_cell_title.xlsx": log_two_cell_title_xlsx,
    "excel_date_only.xlsx": excel_date_only_xlsx, "excel_date_only.csv": excel_date_only_csv,
    "time_column_epoch_date.csv": time_column_epoch_date, "superscript_tap_log.csv": superscript_tap_log,
    "superscript_tap_map.csv": superscript_tap_map, "junctions_2col.csv": junctions_2col,
    "qualifier_interleaved.csv": qualifier_interleaved, "unknown_taps_many.csv": unknown_taps_many,
    "csv_tiny_rows.csv": csv_tiny_rows, "net2_clock_log.csv": net2_clock_log, "net2_tap_map.csv": net2_tap_map,
    **{f"ugl_header_{k}.csv": (lambda h=h: ugl_header_log(h)) for k, h in enumerate(UGL_HEADERS)},
    **{f"units_column_{k}.csv": (lambda u=u: units_column_log(u)) for k, u in enumerate(UNIT_SPELLINGS)},
}


# ----------------------------------------------------------------------------- task 16b: taps by coordinates
# The georeference test's truth is an ellipsoidal Lambert conformal conic projection (two standard parallels; Snyder,
# Map Projections: A Working Manual, 1987, equations 15-1 to 15-10), written here and sharing no code with the ENU
# step it tests (importers.enu_from_latlon).  Its zone is invented (standard parallels 38.5 and 41.5, origin 37.0 N,
# central meridian 100.0 W, on the WGS84 ellipsoid), and so is every position below; no fixture is a real place's map.
LCC_ZONE = {"lat1": 38.5, "lat2": 41.5, "lat0": 37.0, "lon0": -100.0}
LCC_CENTRE = (40.0, -100.6)           # every fixture's taps sit around this point, 51 km west of the central meridian
GEOREF_ZONES_KM = (20, 50)            # lcc_taps_20km.csv and lcc_taps_50km.csv: square zones of this side
GEOREF_TEST_TAPS = 400                # test taps per zone, uniform over the square (default_rng(0))
# control taps, as fractions of the zone's half side (east, north): the first three make the 3-tap fit, all six the
# checked fit; the test taps reach the corners, outside the first three's triangle
GEOREF_CONTROLS = ((-0.8, -0.8), (0.8, -0.7), (-0.1, 0.9), (0.85, 0.8), (-0.85, 0.75), (0.05, -0.1))
COLLINEAR_OFFSETS_M = ((-900.0, -2.0), (-300.0, 3.0), (250.0, -1.5), (850.0, 2.5))   # along a street, 2 m off it
_M_PER_DEG = 111_320.0                # only to lay positions out in degrees; the truth is the projection


def _lcc_constants(z: dict = LCC_ZONE) -> tuple:
    a, f = z.get("a", 6378137.0), z.get("f", 1 / 298.257223563)        # WGS84 unless the zone names its ellipsoid
    e = math.sqrt(f * (2 - f))
    m = lambda p: math.cos(p) / math.sqrt(1 - (e * math.sin(p)) ** 2)
    t = lambda p: math.tan(math.pi / 4 - p / 2) / ((1 - e * math.sin(p)) / (1 + e * math.sin(p))) ** (e / 2)
    p1, p2, p0 = (math.radians(z[k]) for k in ("lat1", "lat2", "lat0"))
    n = (math.log(m(p1)) - math.log(m(p2))) / (math.log(t(p1)) - math.log(t(p2)))
    F = m(p1) / (n * t(p1) ** n)
    return a, e, n, F, t, a * F * t(p0) ** n


def lcc_forward(lat: float, lon: float, z: dict = LCC_ZONE) -> tuple[float, float]:
    """(x, y) in metres of a WGS84 latitude and longitude in the invented Lambert conformal conic zone (or in zone z,
    on its own ellipsoid when it names one: SNYDER_LCC's worked example checks the formulas)."""
    a, e, n, F, t, rho0 = _lcc_constants(z)
    rho = a * F * t(math.radians(lat)) ** n
    th = n * math.radians(lon - z["lon0"])
    return rho * math.sin(th), rho0 - rho * math.cos(th)


def lcc_inverse(x: float, y: float, z: dict = LCC_ZONE) -> tuple[float, float]:
    """The latitude and longitude of (x, y) in the zone, by Newton's method on lcc_forward (to 1e-7 m)."""
    lat, lon = LCC_CENTRE
    for _ in range(50):
        fx, fy = lcc_forward(lat, lon, z)
        rx, ry = x - fx, y - fy
        if math.hypot(rx, ry) < 1e-7:
            break
        h = 1e-6
        ax, ay = lcc_forward(lat + h, lon, z)
        bx, by = lcc_forward(lat, lon + h, z)
        j11, j21, j12, j22 = (ax - fx) / h, (ay - fy) / h, (bx - fx) / h, (by - fy) / h
        det = j11 * j22 - j12 * j21
        lat += (j22 * rx - j12 * ry) / det
        lon += (-j21 * rx + j11 * ry) / det
    return lat, lon


def _lcc_taps(side_km: int) -> bytes:
    """tap_id, role (control or test), lat, lon, and the truth x_m, y_m: the projection's x and y minus the centre's.
    Controls at GEOREF_CONTROLS, test taps uniform over the square zone, laid out in degrees around LCC_CENTRE."""
    import numpy as np
    half = side_km * 500.0
    lat_c, lon_c = LCC_CENTRE
    dlat, dlon = half / _M_PER_DEG, half / (_M_PER_DEG * math.cos(math.radians(lat_c)))
    x0, y0 = lcc_forward(lat_c, lon_c)
    rng = np.random.default_rng(0)
    pts = [(f"C{k + 1}", "control", lat_c + fy * dlat, lon_c + fx * dlon) for k, (fx, fy) in enumerate(GEOREF_CONTROLS)]
    u = rng.uniform(-1, 1, size=(GEOREF_TEST_TAPS, 2))
    pts += [(f"T{k + 1:03d}", "test", lat_c + v * dlat, lon_c + w * dlon) for k, (w, v) in enumerate(u)]
    rows = []
    for tap, role, la, lo in pts:
        x, y = lcc_forward(la, lo)
        rows.append((tap, role, f"{la:.9f}", f"{lo:.9f}", f"{x - x0:.4f}", f"{y - y0:.4f}"))
    return _csv_text("tap_id,role,lat,lon,x_m,y_m", rows).encode()


def lcc_taps_20km() -> bytes:
    return _lcc_taps(20)


def lcc_taps_50km() -> bytes:
    return _lcc_taps(50)


def collinear_controls() -> bytes:
    """Four control taps along one street (1.8 km long, at most 3 m off its line), and two taps to place."""
    lat_c, lon_c = LCC_CENTRE
    x0, y0 = lcc_forward(lat_c, lon_c)
    rows = []
    for k, (e, n) in enumerate(COLLINEAR_OFFSETS_M):
        la, lo = lcc_inverse(x0 + e, y0 + n)
        rows.append((f"C{k + 1}", "control", f"{la:.9f}", f"{lo:.9f}", f"{e:.4f}", f"{n:.4f}"))
    for k, (e, n) in enumerate(((0.0, 400.0), (-500.0, -350.0))):
        la, lo = lcc_inverse(x0 + e, y0 + n)
        rows.append((f"T{k + 1}", "test", f"{la:.9f}", f"{lo:.9f}", f"{e:.4f}", f"{n:.4f}"))
    return _csv_text("tap_id,role,lat,lon,x_m,y_m", rows).encode()


# The synthetic urban grid: invented, in metres, to scale (every declared length is its drawn length), 2 km by 2 km.
# Streets every 200 m (junctions J<i><j> at 200 i, 200 j); 8 in mains on every second street, 6 in on the others, and
# the street at y = 1000 m a 4 in line.  Dual mains on two avenues: a second 6 in main 12 m north of y = 600 m
# (junctions D<i>03) and 14 m east of x = 1400 m (junctions E07<j>), tied to the grid at every crossing.  A 24 in
# transmission main 20 m north of the 4 in line, at y = 1020 m, from the reservoir west of the grid (junctions M0 to
# M4 every 500 m), fed into the grid by 12 in pipes at x = 0, 1000 and 2000 m.  One curved 6 in street with vertices.
URBAN_N = 10                         # blocks per side (11 junctions)
URBAN_BLOCK_M = 200.0
URBAN_DUAL = {"D": ("y", 3, 12.0), "E": ("x", 7, 14.0)}   # name: (axis of the avenue, its index, offset in m)
URBAN_FOUR_IN_ROW = 5                # the street at y = 1000 m
URBAN_TM_Y = 1020.0
URBAN_TM_X = (0.0, 500.0, 1000.0, 1500.0, 2000.0)
URBAN_TM_FEEDS = ((0, 0), (2, 5), (4, 10))   # (main junction index, grid column i) at the 4 in row
URBAN_RESERVOIR = ("R1", (-400.0, URBAN_TM_Y), 70.0)
URBAN_CURVE = ("C1", (2, 8), (4, 8), ((480.0, 1650.0), (600.0, 1680.0), (720.0, 1650.0)))
URBAN_PATTERN = (0.6, 0.5, 0.5, 0.5, 0.6, 0.8, 1.1, 1.4, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0, 1.0, 1.1, 1.2, 1.4, 1.5, 1.3,
                 1.1, 0.9, 0.8, 0.7)
IN_MM = 25.4


def _urban_parts() -> tuple[dict, list, dict]:
    """(nodes: id -> (x, y, elevation, demand L/s), pipes: (id, start, end, diameter in, roughness, vertices),
    the transmission main's pipe IDs and the dual mains' pipe IDs)."""
    nodes, pipes = {}, []
    elev = lambda x: round(10.0 + 0.002 * x, 3)
    J = lambda i, j: f"J{i:02d}{j:02d}"
    for i in range(URBAN_N + 1):
        for j in range(URBAN_N + 1):
            x, y = i * URBAN_BLOCK_M, j * URBAN_BLOCK_M
            nodes[J(i, j)] = (x, y, elev(x), 0.4)
    for i in range(URBAN_N + 1):
        for j in range(URBAN_N + 1):
            if i < URBAN_N:
                d = 4 if j == URBAN_FOUR_IN_ROW else (8 if j % 2 == 0 else 6)
                pipes.append((f"H{i:02d}{j:02d}", J(i, j), J(i + 1, j), d, 100.0, ()))
            if j < URBAN_N:
                pipes.append((f"V{i:02d}{j:02d}", J(i, j), J(i, j + 1), 8 if i % 2 == 0 else 6, 100.0, ()))
    dual = []
    for name, (axis, k, off) in URBAN_DUAL.items():
        for m in range(URBAN_N + 1):
            if axis == "y":
                nid, x, y, base = f"{name}{m:02d}{k:02d}", m * URBAN_BLOCK_M, k * URBAN_BLOCK_M + off, J(m, k)
            else:
                nid, x, y, base = f"{name}{k:02d}{m:02d}", k * URBAN_BLOCK_M + off, m * URBAN_BLOCK_M, J(k, m)
            nodes[nid] = (x, y, elev(x), 0.2)
            pipes.append((f"{name}C{m:02d}", base, nid, 6, 100.0, ()))
            if m:
                prev = f"{name}{m - 1:02d}{k:02d}" if axis == "y" else f"{name}{k:02d}{m - 1:02d}"
                pipes.append((f"{name}M{m:02d}", prev, nid, 6, 100.0, ()))
                dual.append(f"{name}M{m:02d}")
    rid, (rx, ry), _ = URBAN_RESERVOIR
    tm = []
    for k, x in enumerate(URBAN_TM_X):
        nodes[f"M{k}"] = (x, URBAN_TM_Y, elev(x), 0.0)
        pipes.append((f"TM{k}", rid if k == 0 else f"M{k - 1}", f"M{k}", 24, 120.0, ()))
        tm.append(f"TM{k}")
    for k, i in URBAN_TM_FEEDS:
        pipes.append((f"F{k}", f"M{k}", J(i, URBAN_FOUR_IN_ROW), 12, 110.0, ()))
    cid, a, b, verts = URBAN_CURVE
    pipes.append((cid, J(*a), J(*b), 6, 100.0, verts))
    return nodes, pipes, {"transmission": tm, "dual": dual}


def _polyline_length(pts) -> float:
    return sum(math.dist(pts[k], pts[k + 1]) for k in range(len(pts) - 1))


def synthetic_urban_grid() -> bytes:
    """The synthetic urban grid as an EPANET .inp (LPS, Hazen-Williams, metres), written by hand here."""
    nodes, pipes, _ = _urban_parts()
    rid, (rx, ry), head = URBAN_RESERVOIR
    xy = {n: (v[0], v[1]) for n, v in nodes.items()}
    xy[rid] = (rx, ry)
    out = ["[TITLE]", "Synthetic urban grid: invented, in metres, to scale (ResidualMap fixture)", "", "[JUNCTIONS]",
           ";ID\tElev\tDemand\tPattern"]
    out += [f" {n}\t{v[2]:.3f}\t{v[3]:.2f}\tP1\t;" for n, v in nodes.items()]
    out += ["", "[RESERVOIRS]", ";ID\tHead\tPattern", f" {rid}\t{head:.1f}\t\t;", "", "[PIPES]",
            ";ID\tNode1\tNode2\tLength\tDiameter\tRoughness\tMinorLoss\tStatus"]
    for pid, a, b, d, c, verts in pipes:
        length = _polyline_length([xy[a], *verts, xy[b]])
        out.append(f" {pid}\t{a}\t{b}\t{length:.3f}\t{d * IN_MM:.1f}\t{c:g}\t0\tOpen\t;")
    out += ["", "[PATTERNS]", ";ID\tMultipliers"]
    for k in range(0, len(URBAN_PATTERN), 6):
        out.append(" P1\t" + "\t".join(f"{v:g}" for v in URBAN_PATTERN[k:k + 6]))
    out += ["", "[TIMES]", " Duration\t24:00", " Hydraulic Timestep\t1:00", " Quality Timestep\t0:05",
            " Pattern Timestep\t1:00", " Report Timestep\t1:00", " Start ClockTime\t12 am", "", "[OPTIONS]",
            " Units\tLPS", " Headloss\tH-W", " Quality\tNone", "", "[COORDINATES]", ";Node\tX-Coord\tY-Coord"]
    out += [f" {n}\t{x:.3f}\t{y:.3f}" for n, (x, y) in xy.items()]
    out += ["", "[VERTICES]", ";Link\tX-Coord\tY-Coord"]
    out += [f" {pid}\t{x:.3f}\t{y:.3f}" for pid, _, _, _, _, verts in pipes for x, y in verts]
    out += ["", "[END]", ""]
    return "\n".join(out).encode()


def ky4_half_reversed() -> bytes:
    """ky4 with every second pipe of [PIPES] drawn the other way: its start and end nodes swapped and its vertices in
    [VERTICES] in reverse order (check-valve pipes, none in ky4, are left as they are), so the drawn direction and
    the sign of the simulated flow both flip on half the pipes.  Nothing else changes."""
    text = _lib("ky4").decode("utf-8")
    lines = text.splitlines(keepends=True)
    sec, k, flipped, out = None, 0, set(), []
    for ln in lines:
        s = ln.strip()
        if s.startswith("["):
            sec = s.upper()
            out.append(ln)
            continue
        if sec == "[PIPES]" and s and not s.startswith(";"):
            tok = ln.split(";")[0].split()
            if k % 2 == 1 and not (len(tok) > 7 and tok[7].upper() == "CV"):
                body, sep, comment = ln.partition(";")
                end = ln[len(ln.rstrip("\r\n")):]
                tok[1], tok[2] = tok[2], tok[1]
                ln = " " + "\t".join(tok) + ("\t;" + comment.rstrip("\r\n") if sep else "") + end
                flipped.add(tok[0])
            k += 1
        out.append(ln)
    # the vertices of each flipped pipe, in reverse order, in the places its lines held
    sec, idx = None, {}
    for i, ln in enumerate(out):
        s = ln.strip()
        if s.startswith("["):
            sec = s.upper()
            continue
        if sec == "[VERTICES]" and s and not s.startswith(";"):
            pid = s.split()[0]
            if pid in flipped:
                idx.setdefault(pid, []).append(i)
    for pid, rows in idx.items():
        vals = [out[i] for i in rows]
        for i, v in zip(rows, reversed(vals)):
            out[i] = v
    return "".join(out).encode("utf-8")


def ky4_flipped_pipes() -> list[str]:
    """The pipes ky4_half_reversed draws the other way: every second pipe of ky4's [PIPES], in order."""
    wn = _wn("ky4")
    return [p for k, p in enumerate(wn.pipe_name_list) if k % 2 == 1
            and not wn.get_link(p).check_valve]


def _net3_some_junctions(k: int = 6) -> list[tuple[str, float, float]]:
    wn = _wn("Net3")
    out = []
    for j in ("10", "15", "20", "35", "40", "50")[:k]:
        x, y = wn.get_node(j).coordinates
        out.append((j, float(x), float(y)))
    return out


def net3_taps_xy() -> bytes:
    """Taps in Net3's own x and y (a schematic drawing): refused, with 'use junction IDs'."""
    return _csv_text("tap_id,x,y", [(f"T{k + 1}", f"{x + 0.5:.2f}", f"{y - 0.5:.2f}")
                                    for k, (_, x, y) in enumerate(_net3_some_junctions())]).encode()


def net3_taps_latlon() -> bytes:
    """Taps in latitude and longitude, four of them control taps at Net3 junctions: refused all the same."""
    lat_c, lon_c = LCC_CENTRE
    rows = []
    for k, (j, x, y) in enumerate(_net3_some_junctions()):
        rows.append((f"T{k + 1}", f"{lat_c + 0.001 * k:.6f}", f"{lon_c + 0.0013 * (k % 3):.6f}", j if k < 4 else ""))
    return _csv_text("tap_id,lat,lon,junction_id", rows).encode()


# The tap map by coordinates template, a fixture only: four control taps at junctions of the synthetic urban grid and
# four taps to place, in latitude and longitude, the grid laid at URBAN_ORIGIN in the invented zone (model x = the
# projection's x minus URBAN_ORIGIN).  It is not shipped in docs/templates: nothing in the app or the pilot reads a
# tap map by coordinates yet (task 16b stopped at its pre-registered bar before those steps were built).
URBAN_ORIGIN = lcc_forward(*LCC_CENTRE)
TEMPLATE_CONTROLS = (("C1", "J0101"), ("C2", "J0901"), ("C3", "J0909"), ("C4", "J0109"))
TEMPLATE_TAPS = (("T1", (503.0, 403.0), ""), ("T2", (1197.0, 1110.0), ""), ("T3", (700.0, 1004.0), ""),
                 ("T4", (1250.0, 1021.0), "yes"))
TEMPLATE_NOTES = {"C1": "control tap: its junction is known and its position was taken at that junction; the "
                        "positions and junction IDs here are from the synthetic test grid: replace them with your own",
                  "C2": "control tap", "C3": "control tap", "C4": "control tap: four or more let the fit be checked",
                  "T1": "a tap to place: ResidualMap suggests its junction for you to check",
                  "T2": "a tap to place", "T3": "a tap to place, beside a large main it is not on",
                  "T4": "on_main yes: this tap is on the large main"}


def _urban_latlon(x: float, y: float) -> tuple[str, str]:
    la, lo = lcc_inverse(URBAN_ORIGIN[0] + x, URBAN_ORIGIN[1] + y)
    return f"{la:.6f}", f"{lo:.6f}"


def coordinates_template() -> bytes:
    nodes, _, _ = _urban_parts()
    rows = []
    for tap, j in TEMPLATE_CONTROLS:
        la, lo = _urban_latlon(nodes[j][0], nodes[j][1])
        rows.append((tap, la, lo, j, "", TEMPLATE_NOTES[tap]))
    for tap, (x, y), main in TEMPLATE_TAPS:
        la, lo = _urban_latlon(x, y)
        rows.append((tap, la, lo, "", main, TEMPLATE_NOTES[tap]))
    return _csv_text("tap_id,lat,lon,junction_id,on_main,notes", rows).encode()


URBAN_LOG_MONTHS = ("2026-01", "2026-02", "2026-03", "2026-04")   # enough for the pilot's hold-out (3 months, then 1)


def urban_coord_log() -> bytes:
    """A grab log at the template's eight taps, one reading a month each for four months (invented readings)."""
    taps = [t for t, _ in TEMPLATE_CONTROLS] + [t for t, _, _ in TEMPLATE_TAPS]
    rows = []
    for m, month in enumerate(URBAN_LOG_MONTHS):
        for k, tap in enumerate(taps):
            rows.append((f"{month}-{5 + k:02d}", f"{9 + k % 6:02d}:{10 + 5 * k:02d}", tap,
                         f"{0.35 + 0.05 * ((k + m) % 9):.2f}"))
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


FIXTURES.update({"synthetic_urban_grid.inp": synthetic_urban_grid, "ky4_half_reversed.inp": ky4_half_reversed,
                 "lcc_taps_20km.csv": lcc_taps_20km, "lcc_taps_50km.csv": lcc_taps_50km,
                 "collinear_controls.csv": collinear_controls, "net3_taps_xy.csv": net3_taps_xy,
                 "net3_taps_latlon.csv": net3_taps_latlon, "tap_map_by_coordinates.csv": coordinates_template,
                 "urban_coord_log.csv": urban_coord_log})
TAP_FIXTURES = ("synthetic_urban_grid.inp", "ky4_half_reversed.inp", "lcc_taps_20km.csv", "lcc_taps_50km.csv",
                "collinear_controls.csv", "net3_taps_xy.csv", "net3_taps_latlon.csv", "tap_map_by_coordinates.csv",
                "urban_coord_log.csv", "net3_no_coords.inp")

# The review of task 16b: its cases as fixtures (the tap maps written as text here, the models rewritten from the
# synthetic urban grid's own text, the zones from the Lambert conformal conic truth), for part 'review' of
# outputs/intake/tap_snapping.json.
# Snyder's worked example of the ellipsoidal Lambert conformal conic (Map Projections: A Working Manual, 1987, p. 296):
# Clarke 1866, standard parallels 33 and 45 N, origin 23 N, 96 W; the point 35 N, 75 W is at these x and y.
SNYDER_LCC = {"zone": {"lat1": 33.0, "lat2": 45.0, "lat0": 23.0, "lon0": -96.0, "a": 6378206.4,
                       "f": 1 - 6356583.8 / 6378206.4},
              "point": (35.0, -75.0), "xy_m": (1894410.9, 1564649.5)}
URBAN_TM_VARIANTS_MM = {"urban_grid_tm_dn400.inp": 400.0, "urban_grid_tm_id388.inp": 388.6,   # the metric 16 in, and a
                        # 16 in main by its internal diameter (15.3 in); the final verification's 16 in plastic mains by
                        # internal diameter, under the 15.0 in cut: AWWA C905 DR14 (14.9 in) and PE100 DN400 SDR11
                        # (12.9 in)
                        "urban_grid_tm_c905_dr14.inp": 378.8, "urban_grid_tm_pe100_sdr11.inp": 327.4}
URBAN_DEG_LAT = 45.0                 # urban_grid_degrees.inp: the grid drawn in degrees around 45 N, 100 W
URBAN_PLACEHOLDERS = ("J0808", "J0604")   # urban_grid_placeholders.inp: drawn at 0,0, as an export does for a node
URBAN_DRAWN_MAX_X = 1200.0           # urban_grid_part_drawn.inp: no position for any node east of this
NARROW_ZONE_M = (8000.0, 300.0)      # lcc_narrow_zone.csv: a valley 8 km long and 300 m wide
NARROW_CONTROL_NOISE_M = 3.0         # radial RMS GPS error on its control taps' positions (default_rng(5))
ANTIMERIDIAN_SHIFT_DEG = 180.0 - LCC_CENTRE[1]   # lcc_antimeridian.csv: the 20 km zone moved to 180 degrees
TRIANGLE_ONLY_OFFSETS_M = ((-1000.0, 0.0), (1000.0, 0.0), (0.0, 60.0), (0.0, -60.0))   # 17 times longer than wide


def _urban_rewrite(text: str, coord=None, tm_mm: float | None = None, drop=None) -> bytes:
    """The urban grid's .inp text with its [COORDINATES] and [VERTICES] rewritten by coord(node or link, x, y) ->
    (x text, y text), the lines for which drop(node, x, y) is true removed, and the transmission main's diameter set
    to tm_mm."""
    out, sec = [], None
    for ln in text.split("\n"):
        s = ln.strip()
        if s.startswith("["):
            sec = s.upper()
        elif s and not s.startswith(";") and sec in ("[COORDINATES]", "[VERTICES]"):
            name, x, y = s.split("\t")
            if sec == "[COORDINATES]" and drop and drop(name, float(x), float(y)):
                continue
            if coord:
                ln = " " + "\t".join((name, *coord(name, float(x), float(y))))
        elif s and not s.startswith(";") and sec == "[PIPES]" and tm_mm is not None and s.split("\t")[0].startswith("TM"):
            tok = s.split("\t")
            tok[4] = f"{tm_mm:.1f}"
            ln = " " + "\t".join(tok)
        out.append(ln)
    return "\n".join(out).encode()


def _urban_text() -> str:
    return synthetic_urban_grid().decode()


def urban_grid_degrees() -> bytes:
    """The urban grid drawn in degrees of longitude and latitude (lengths in the file still in metres), as a GIS
    export that keeps WGS84 writes it."""
    k = 111_320.0 * math.cos(math.radians(URBAN_DEG_LAT))
    return _urban_rewrite(_urban_text(), coord=lambda n, x, y: (f"{-100.0 + x / k:.8f}", f"{URBAN_DEG_LAT + y / 110_574.0:.8f}"))


def urban_grid_miles() -> bytes:
    """The urban grid drawn in miles (lengths in the file still in metres)."""
    return _urban_rewrite(_urban_text(), coord=lambda n, x, y: (f"{x / 1609.344:.7f}", f"{y / 1609.344:.7f}"))


def urban_grid_placeholders() -> bytes:
    """The urban grid with two junctions drawn at 0,0, as an export writes a node it has no position for."""
    return _urban_rewrite(_urban_text(), coord=lambda n, x, y: ("0.000", "0.000") if n in URBAN_PLACEHOLDERS else
                          (f"{x:.3f}", f"{y:.3f}"))


def urban_grid_part_drawn() -> bytes:
    """The urban grid with no position for any node east of URBAN_DRAWN_MAX_X (a model drawn in part)."""
    return _urban_rewrite(_urban_text(), drop=lambda n, x, y: x > URBAN_DRAWN_MAX_X)


def _urban_tm(mm: float):
    return lambda: _urban_rewrite(_urban_text(), tm_mm=mm)


def _urban_ll_rows(points, junction: bool = False) -> list:
    rows = []
    for tap, (x, y), j in points:
        la, lo = _urban_latlon(x, y)
        rows.append((tap, la, lo, f"{x:.1f}", f"{y:.1f}", *((j,) if junction else ())))
    return rows


_REVIEW_TAPS = (("T1", (503.0, 403.0), ""), ("T2", (1197.0, 1110.0), ""), ("T3", (700.0, 1004.0), ""))


def taps_both_pairs_latlon_blank() -> bytes:
    """Latitude and longitude columns left blank, x and y filled (the template with x and y columns added)."""
    return _csv_text("tap_id,lat,lon,x,y", [(t, "", "", f"{x:g}", f"{y:g}") for t, (x, y), _ in _REVIEW_TAPS]).encode()


def taps_both_pairs_no_controls() -> bytes:
    """Both pairs filled on every row and no junction column: x and y place every tap, latitude and longitude would
    need control taps."""
    return _csv_text("tap_id,lat,lon,x,y", _urban_ll_rows(_REVIEW_TAPS)).encode()


def taps_both_pairs_controls() -> bytes:
    """Both pairs filled, with the template's four control taps: latitude and longitude by default, x and y when
    those columns are chosen."""
    nodes, _, _ = _urban_parts()
    pts = [(tap, nodes[j][:2], j) for tap, j in TEMPLATE_CONTROLS] + [(t, p, "") for t, p, _ in _REVIEW_TAPS]
    return _csv_text("tap_id,lat,lon,x,y,junction_id", _urban_ll_rows(pts, junction=True)).encode()


def taps_thousands_ambiguous() -> bytes:
    """x and y with thousands commas ('1,197') and nothing else in the column to say so."""
    return b'tap_id,x,y\nT1,503,403\nT2,"1,197","1,110"\nT3,700,"1,004"\n'


def taps_thousands_settled() -> bytes:
    """The same, with values of 1,000 or more written plainly in both columns, which settle the commas as thousands."""
    return b'tap_id,x,y\nT1,1503.5,1403\nT2,"1,197","1,110"\n'


def taps_decimal_comma() -> bytes:
    """A semicolon file with decimal commas ('503,5'), which settle '1,197' as 1.197."""
    return b"tap_id;x;y\nT1;503,5;403,25\nT2;1,197;1,110\n"


def taps_letters_swapped() -> bytes:
    """The template with T1's latitude written with W and its longitude with N."""
    rows = [r.split(",") for r in coordinates_template().decode().strip().split("\n")]
    for r in rows[1:]:
        if r[0] == "T1":
            r[1], r[2] = f"{r[1]} W", f"{r[2].lstrip('-')} N"
    return ("\n".join(",".join(r) for r in rows) + "\n").encode()


def taps_latlon_swapped() -> bytes:
    """The template with its latitude and longitude values in each other's columns."""
    rows = [r.split(",") for r in coordinates_template().decode().strip().split("\n")]
    for r in rows[1:]:
        r[1], r[2] = r[2], r[1]
    return ("\n".join(",".join(r) for r in rows) + "\n").encode()


def taps_junction_header_unknown() -> bytes:
    """A tap map by ID whose junction column has a header nothing guesses, with x and y too."""
    return b"tap_id,EPANET_Name,x,y\nT1,J0101,203,198\nT2,J0909,1797,1802\n"


def taps_duplicate_positions() -> bytes:
    """Tap T1 listed twice, at two positions 900 m apart."""
    return b"tap_id,x,y\nT1,503,403\nT1,1197,1110\nT2,700,1004\n"


def taps_duplicate_junctions() -> bytes:
    """Tap T1 listed twice, with two different junctions."""
    return b"tap_id,x,y,junction_id\nT1,503,403,J0101\nT1,505,404,J0909\nT2,700,1004,\n"


def taps_absurd() -> bytes:
    """A position far beyond any drawing (x = 1e200) beside a good one."""
    return b"tap_id,x,y\nT1,1e200,5\nT2,503,403\n"


def taps_beside_main() -> bytes:
    """One position 1 m from the 24 in main and 21 m from the 4 in line, untagged and tagged on_main."""
    return b"tap_id,x,y,on_main\nT1,1250,1021,\nT2,1250,1021,yes\n"


def taps_mapping_missing() -> bytes:
    """A junction column 'node_ref', for a mapping that names 'node ref', a column the file does not have."""
    return b"tap_id,x,y,node_ref\nT1,503,403,J0101\nT2,700,1004,\n"


def _lcc_rows(points, controls: int, noise_m: float = 0.0, seed: int = 0) -> bytes:
    """tap_id, role, lat, lon, x_m, y_m for points (east, north) in metres around LCC_CENTRE's projected position,
    the first `controls` of them control taps whose latitude and longitude carry noise_m of radial RMS GPS error."""
    import numpy as np
    x0, y0 = lcc_forward(*LCC_CENTRE)
    rng = np.random.default_rng(seed)
    rows = []
    for k, (e, n) in enumerate(points):
        de, dn = rng.normal(0.0, noise_m / math.sqrt(2.0), 2) if k < controls and noise_m else (0.0, 0.0)
        la, lo = lcc_inverse(x0 + e + de, y0 + n + dn)
        tap = f"C{k + 1}" if k < controls else f"T{k - controls + 1:03d}"
        rows.append((tap, "control" if k < controls else "test", f"{la:.9f}", f"{lo:.9f}", f"{e:.4f}", f"{n:.4f}"))
    return _csv_text("tap_id,role,lat,lon,x_m,y_m", rows).encode()


def lcc_narrow_zone() -> bytes:
    """A narrow valley zone: six control taps at its corners and the middle of its long sides (positions with
    NARROW_CONTROL_NOISE_M of GPS error), 300 test taps inside it (default_rng(4)).  The relative tests alone (more
    than COLLINEAR_MAX_COND times longer than wide, or a triangle under COLLINEAR_MIN_AREA) refuse it."""
    import numpy as np
    L, W = NARROW_ZONE_M
    ctrl = [(-L / 2, -W / 2), (L / 2, -W / 2), (L / 2, W / 2), (-L / 2, W / 2), (0.0, W / 2), (0.0, -W / 2)]
    u = np.random.default_rng(4).uniform(-0.5, 0.5, size=(300, 2)) * np.array([L, W])
    return _lcc_rows(ctrl + [tuple(p) for p in u], len(ctrl), NARROW_CONTROL_NOISE_M, 5)


def lcc_antimeridian() -> bytes:
    """lcc_taps_20km.csv moved east by ANTIMERIDIAN_SHIFT_DEG, so the zone's centre is at 180 degrees, its
    longitudes written from -180 to 180 as a GPS writes them; its x_m and y_m are unchanged (the projection moves
    with its central meridian)."""
    rows = [r.split(",") for r in lcc_taps_20km().decode().strip().split("\n")]
    for r in rows[1:]:
        lo = float(r[3]) + ANTIMERIDIAN_SHIFT_DEG
        r[3] = f"{lo - 360.0 if lo > 180.0 else lo:.9f}"
    return ("\n".join(",".join(r) for r in rows) + "\n").encode()


def triangle_only_controls() -> bytes:
    """Four control taps 14 times longer than wide (under the condition-number limit) whose largest triangle is below
    the area limit, within about 42 m of their line, and two taps to place."""
    return _lcc_rows(list(TRIANGLE_ONLY_OFFSETS_M) + [(300.0, 20.0), (-400.0, -30.0)], len(TRIANGLE_ONLY_OFFSETS_M))


# The final verification of task 16b: its cases as fixtures, for part 'review' of outputs/intake/tap_snapping.json.
URBAN_SENTINEL_NODE = "J0505"        # drawn at a no-data value: the float32 maximum some exports write, or nan
URBAN_SENTINELS = {"urban_grid_sentinel_node.inp": "-3.4028235e38", "urban_grid_nan_node.inp": "nan"}
URBAN_DEG_PLACEHOLDER = ("J0808", "9999")   # urban_grid_degrees_placeholder.inp: one node at 9999,9999
STRIP_CONTROLS_M = ((-2000.0, -55.0), (-1000.0, 55.0), (0.0, -55.0), (1000.0, 55.0), (2000.0, -55.0), (0.0, 55.0))
STRIP_SIDE_M, STRIP_GPS_M = 4000.0, 3.0    # lcc_strip_controls.csv: control taps along one avenue, 3 m of GPS error
                                           # (default_rng(0)), 400 test taps over the 4 km square (default_rng(4))


def _urban_degrees(n, x, y, lat_as_x: bool = False):
    k = 111_320.0 * math.cos(math.radians(URBAN_DEG_LAT))
    lon, lat = f"{-100.0 + x / k:.8f}", f"{URBAN_DEG_LAT + y / 110_574.0:.8f}"
    return (lat, lon) if lat_as_x else (lon, lat)


def urban_grid_degrees_lat_as_x() -> bytes:
    """The urban grid drawn in degrees with latitude as x and longitude as y, the order some exports write."""
    return _urban_rewrite(_urban_text(), coord=lambda n, x, y: _urban_degrees(n, x, y, lat_as_x=True))


def urban_grid_degrees_placeholder() -> bytes:
    """The urban grid drawn in degrees (longitude as x) with one node at the placeholder 9999,9999."""
    j, v = URBAN_DEG_PLACEHOLDER
    return _urban_rewrite(_urban_text(), coord=lambda n, x, y: (v, v) if n == j else _urban_degrees(n, x, y))


def _urban_sentinel(value: str):
    return lambda: _urban_rewrite(_urban_text(), coord=lambda n, x, y: (value, value) if n == URBAN_SENTINEL_NODE else
                                  (f"{x:.3f}", f"{y:.3f}"))


def taps_dot_thousands() -> bytes:
    """A European tap map (semicolons, decimal commas) whose values of 1,000 or more carry a dot thousands separator
    ('1.197' for 1197), with nothing in the columns to settle the dots."""
    return b"tap_id;x;y\nT1;503,5;403,25\nT2;1.197;1.110\nT3;700;1.004\n"


def taps_dot_thousands_settled() -> bytes:
    """The same, with '1.503,5' and '1.403,25' (a dot thousands separator for certain) settling both columns."""
    return b"tap_id;x;y\nT1;1.503,5;1.403,25\nT2;1.197;1.110\nT3;700;1.004\n"


def taps_lon_0_to_360() -> bytes:
    """The coordinates template with 360 added to every longitude (259.4 for 100.6 W), as some receivers write them."""
    rows = [r.split(",") for r in coordinates_template().decode().strip().split("\n")]
    for r in rows[1:]:
        r[2] = f"{float(r[2]) + 360.0:.6f}"
    return ("\n".join(",".join(r) for r in rows) + "\n").encode()


def lcc_strip_controls() -> bytes:
    """Six control taps along one avenue, 4 km long and 55 m to either side of it (STRIP_GPS_M of GPS error on each,
    default_rng(0)), and 400 test taps over the 4 km square around it (default_rng(4)): accepted, as a narrow strip
    more than ACROSS_MIN_M across, and the taps far across the avenue are where the fit is uncertain."""
    import numpy as np
    h = STRIP_SIDE_M / 2
    test = np.random.default_rng(4).uniform(-h, h, size=(400, 2))
    return _lcc_rows(list(STRIP_CONTROLS_M) + [tuple(p) for p in test], len(STRIP_CONTROLS_M), STRIP_GPS_M, 0)


STREET_CONTROLS = (("C1", "J0100"), ("C2", "J0500"), ("C3", "J0900"), ("C4", "J0501"))   # along y = 0, one 200 m off
STREET_TAPS = (("T1", (503.0, 203.0)), ("T2", (1197.0, 1110.0)), ("T3", (700.0, 1797.0)))


def taps_controls_along_one_street() -> bytes:
    """Four control taps on the urban grid, three along its street at y = 0 and one 200 m off it (accepted: not a
    narrow strip), and three taps to place, one beside them and two far across that street, in latitude and
    longitude through the invented zone (no GPS error)."""
    nodes, _, _ = _urban_parts()
    rows = [(tap, *_urban_latlon(*nodes[j][:2]), j) for tap, j in STREET_CONTROLS]
    rows += [(tap, *_urban_latlon(x, y), "") for tap, (x, y) in STREET_TAPS]
    return _csv_text("tap_id,lat,lon,junction_id", rows).encode()


FINAL_16B_FIXTURES = {
    "taps_controls_along_one_street.csv": taps_controls_along_one_street,
    **{n: _urban_sentinel(v) for n, v in URBAN_SENTINELS.items()},
    "urban_grid_degrees_lat_as_x.inp": urban_grid_degrees_lat_as_x,
    "urban_grid_degrees_placeholder.inp": urban_grid_degrees_placeholder,
    "taps_dot_thousands.csv": taps_dot_thousands, "taps_dot_thousands_settled.csv": taps_dot_thousands_settled,
    "taps_lon_0_to_360.csv": taps_lon_0_to_360, "lcc_strip_controls.csv": lcc_strip_controls,
}


REVIEW_16B_FIXTURES = {
    **{n: _urban_tm(mm) for n, mm in URBAN_TM_VARIANTS_MM.items()},
    "urban_grid_degrees.inp": urban_grid_degrees, "urban_grid_miles.inp": urban_grid_miles,
    "urban_grid_placeholders.inp": urban_grid_placeholders, "urban_grid_part_drawn.inp": urban_grid_part_drawn,
    "taps_both_pairs_latlon_blank.csv": taps_both_pairs_latlon_blank,
    "taps_both_pairs_no_controls.csv": taps_both_pairs_no_controls,
    "taps_both_pairs_controls.csv": taps_both_pairs_controls, "taps_thousands_ambiguous.csv": taps_thousands_ambiguous,
    "taps_thousands_settled.csv": taps_thousands_settled, "taps_decimal_comma.csv": taps_decimal_comma,
    "taps_letters_swapped.csv": taps_letters_swapped, "taps_latlon_swapped.csv": taps_latlon_swapped,
    "taps_junction_header_unknown.csv": taps_junction_header_unknown,
    "taps_duplicate_positions.csv": taps_duplicate_positions, "taps_duplicate_junctions.csv": taps_duplicate_junctions,
    "taps_absurd.csv": taps_absurd, "taps_beside_main.csv": taps_beside_main,
    "taps_mapping_missing.csv": taps_mapping_missing, "lcc_narrow_zone.csv": lcc_narrow_zone,
    "lcc_antimeridian.csv": lcc_antimeridian, "triangle_only_controls.csv": triangle_only_controls,
}
REVIEW_16B_FIXTURES.update(FINAL_16B_FIXTURES)
FIXTURES.update(REVIEW_16B_FIXTURES)
TAP_FIXTURES = TAP_FIXTURES + tuple(REVIEW_16B_FIXTURES)


# ----------------------------------------------------------------------------- task 20: the hosted demo's limits
# A square grid model of n x n junctions, invented (schematic coordinates, 100 m apart, Hazen-Williams, LPS), for the
# hosted demo's limit checks (the review of task 20): a model far over the 300-junction limit whose zip is well under
# the 5 MB file limit (n = 180: 32,400 junctions, about 4 MB of text in a zip of about 1 MB), and the same padded with
# comment lines so that its zip unpacks past 5 MB.  Not in FIXTURES, since build_all would then write megabytes for
# every check that builds them all: the task-20 checks call these directly.
def cloud_grid_inp(n: int, pad_bytes: int = 0) -> bytes:
    """An n x n grid model fed from one reservoir; pad_bytes of ';' comment lines after its title."""
    out = ["[TITLE]", f"Synthetic {n} x {n} grid (ResidualMap fixture, invented)"]
    pad = ";" + "padding comment line, ignored by EPANET " * 2 + "\n"
    text_pad = pad * (pad_bytes // len(pad) + 1) if pad_bytes else ""
    out += [text_pad.rstrip("\n")] if text_pad else []
    out += ["", "[JUNCTIONS]"]
    out += [f"J{i}_{j} {100 + (7 * i + 3 * j) % 20} 0.1" for i in range(n) for j in range(n)]
    out += ["", "[RESERVOIRS]", "R1 200", "", "[PIPES]", "P0 R1 J0_0 100 600 130 0 Open"]
    k = 1
    for i in range(n):
        for j in range(n):
            for di, dj in ((0, 1), (1, 0)):
                if i + di < n and j + dj < n:
                    out.append(f"P{k} J{i}_{j} J{i + di}_{j + dj} {50 + (11 * k) % 150} {(150, 200, 300)[k % 3]} "
                               f"{90 + (13 * k) % 50} 0 Open")
                    k += 1
    out += ["", "[COORDINATES]", "R1 -100 -100"]
    out += [f"J{i}_{j} {j * 100} {i * 100}" for i in range(n) for j in range(n)]
    out += ["", "[TIMES]", "Duration 24:00", "Hydraulic Timestep 1:00", "", "[OPTIONS]", "Units LPS",
            "Headloss H-W", "", "[END]", ""]
    return "\n".join(out).encode("ascii")


def cloud_grid_zip(n: int, pad_bytes: int = 0) -> bytes:
    """cloud_grid_inp(n, pad_bytes) as the one member of a zip, as a model is often emailed."""
    return _zip([(f"grid{n}.inp", cloud_grid_inp(n, pad_bytes))])


# The final verification of task 20: models under both of those limits that cost far more than their size says, each
# the 17 x 17 grid (289 junctions) with one thing added, and a small grab-log workbook whose shared strings unpack to
# far more than the file is.  Invented; not in FIXTURES either.
CLOUD_COST_N = 17
CLOUD_PATTERN = "1 " + " ".join(f"{0.8 + 0.4 * ((i * 7) % 10) / 10:.2f}" for i in range(24))


def _cloud_grid_text(n: int = CLOUD_COST_N) -> str:
    return cloud_grid_inp(n).decode("ascii")


def cloud_parallel_pipes_inp(total_bytes: int = 4_700_000, n: int = CLOUD_COST_N) -> bytes:
    """The n x n grid and parallel pipes between neighbouring junctions until the file is about total_bytes (the
    verification's 4.7 MB file had about 125,000 links: 127 s and 3.7 GB to prepare under the old limits)."""
    text = _cloud_grid_text(n)
    jn = [f"J{i}_{j}" for i in range(n) for j in range(n)]
    extra, size, k = [], len(text), 0
    while size < total_bytes:
        a = k % (n * n)
        line = f"Q{k} {jn[a]} {jn[(a + 1) % (n * n)]} 100 150 100 0 Open"
        extra.append(line)
        size += len(line) + 1
        k += 1
    return text.replace("\n\n[COORDINATES]", "\n" + "\n".join(extra) + "\n\n[COORDINATES]", 1).encode("ascii")


def cloud_pattern_step_inp(step: str = "0:01", n: int = CLOUD_COST_N) -> bytes:
    """The n x n grid with a 24-step demand pattern and PATTERN TIMESTEP step (h:mm; the verification's 1-minute step
    took 38 s a 7-day run, and a 10-second step more than 400 s, against 0.7 s at an hour)."""
    return _cloud_grid_text(n).replace("[TIMES]", f"[PATTERNS]\n{CLOUD_PATTERN}\n\n[TIMES]\nPattern Timestep {step}",
                                       1).encode("ascii")


def cloud_rules_inp(n_rules: int = 300, step: str = "0:00:01", n: int = CLOUD_COST_N) -> bytes:
    """The n x n grid with a tank and n_rules rule-based controls checked every step (h:mm:ss; 300 rules every second
    took 34 s a run, every minute 0.6 s)."""
    text = _cloud_grid_text(n).replace("[PIPES]", "[TANKS]\nT0 100 5 0 10 10 0\n\n[PIPES]", 1)
    text = text.replace("\n\n[COORDINATES]", "\nPT0 T0 J0_5 100 150 100 0 Open\n\n[COORDINATES]", 1)
    rules = [f"RULE R{k}\nIF SYSTEM CLOCKTIME >= {k % 12 + 1}:{k % 60:02d} AM\nAND TANK T0 LEVEL ABOVE {k % 7}\n"
             f"THEN PIPE P{k % 500 + 1} STATUS IS OPEN\nELSE PIPE P{k % 500 + 1} STATUS IS CLOSED\n" for k in range(n_rules)]
    text = text.replace("[COORDINATES]", "[RULES]\n" + "\n".join(rules) + "\n\n[COORDINATES]", 1)
    return text.replace("[TIMES]", f"[TIMES]\nRule Timestep {step}", 1).encode("ascii")


def cloud_controls_inp(total_bytes: int = 4_900_000, n: int = CLOUD_COST_N) -> bytes:
    """The n x n grid with timed controls, one a second, until the file is about total_bytes (about 150,000 lines:
    its intake ran past 90 s)."""
    text = _cloud_grid_text(n)
    lines, size, k = [], len(text), 1
    while size < total_bytes:
        line = (f"LINK P{1 + k % 500} {'OPEN' if k % 2 else 'CLOSED'} AT TIME {k // 3600}:{(k // 60) % 60:02d}:"
                f"{k % 60:02d}")
        lines.append(line)
        size += len(line) + 1
        k += 1
    return text.replace("[COORDINATES]", "[CONTROLS]\n" + "\n".join(lines) + "\n\n[COORDINATES]", 1).encode("ascii")


SHARED_STRINGS_TYPE = "application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"


def shared_strings_bomb_xlsx(unpacked_bytes: int = 190_000_000) -> bytes:
    """A 20-row grab log (written by openpyxl) whose workbook also declares a shared-strings part of about
    unpacked_bytes of one repeated string (the verification's: a 466 KB file unpacking to 190 MB, 407 to 1, which
    openpyxl read in 32 s and 1 GB, even in read-only mode)."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["date", "time", "tap_id", "free_chlorine_mgL"])
    for d in range(1, 21):
        ws.append([f"2026-09-{d:02d}", "10:00", f"T{d % 3 + 1}", 0.5])
    buf = io.BytesIO()
    wb.save(buf)
    item = "<si><t>x</t></si>"
    k = unpacked_bytes // len(item)
    strings = (f'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="{k}" '
               f'uniqueCount="{k}">' + item * k + "</sst>").encode("ascii")
    import re
    src = zipfile.ZipFile(io.BytesIO(buf.getvalue()))
    parts = []
    for i in src.infolist():
        data = src.read(i.filename)
        if i.filename == "[Content_Types].xml":
            data = data.replace(b"</Types>", f'<Override PartName="/xl/sharedStrings.xml" '
                                             f'ContentType="{SHARED_STRINGS_TYPE}" /></Types>'.encode("ascii"))
        if i.filename == "docProps/core.xml":        # openpyxl stamps the time of writing: fixed, so the bytes are too
            data = re.sub(rb"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", b"2026-01-01T00:00:00Z", data)
        parts.append((i.filename, data))
    parts.append(("xl/sharedStrings.xml", strings))
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for fname, data in parts:
            zi = zipfile.ZipInfo(fname, date_time=ZIP_TIME)
            zi.external_attr = 0o644 << 16
            z.writestr(zi, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return out.getvalue()


# (fixture, kind it must come back as, the file name may decide something, a backup copy)
SNIFF_TABLE: list[tuple[str, str, bool, bool]] = [
    ("net3_as_txt.txt", "epanet_inp", False, False), ("net3_backup.inp.bak", "epanet_inp", False, True),
    ("inp_with_bom.inp", "epanet_inp", False, False), ("inp_comment_first.inp", "epanet_inp", False, False),
    ("inp_long_title.inp", "epanet_inp", False, False), ("inp_utf16.inp", "epanet_inp", False, False),
    ("net3_cp1252_title.inp", "epanet_inp", False, False), ("emailed_model.zip", "epanet_inp_zip", False, False),
    ("dxf_with_999_comment.dxf", "dxf_ascii", False, False), ("ascii_drawing.dxf.bak", "dxf_ascii", False, True),
    ("binary.dxf", "dxf_binary", False, False), ("dwg_ac1027.dwg", "dwg", False, False),
    ("dwg_ac1032.dwg", "dwg", False, False), ("drawing.bak", "dwg", False, True), ("dgn_v7.dgn", "dgn_v7", False, False),
    ("ole_dgn_v8.dgn", "dgn_v8", False, False),
    ("jet_model.wtg.mdb", "access_db", True, False), ("model.wtg.sqlite", "sqlite_model", True, False),
    ("wtg_settings.wtg", "wtg_settings", False, False),
    ("mains.gpkg", "geopackage", False, False), ("mains_geojson_as.json", "geojson", False, False),
    ("geojson_stateplane.geojson", "geojson", False, False), ("esri_features.json", "esri_json", False, False),
    ("wntr_sidecar.json", "wntr_json", False, False), ("gml_mains.gml", "gml", False, False),
    ("mains.kml", "kml", False, False), ("mains.kmz", "kmz", False, False),
    ("mains_shp_nested_two_layers.zip", "shapefile_zip", False, False), ("utility_gdb.zip", "file_gdb_zip", False, False),
    ("grab_log.xlsx", "xlsx", False, False), ("old_excel_log.xls", "xls_legacy", False, False),
    ("log_semicolon.csv", "delimited_text", False, False), ("log_tab.txt", "delimited_text", False, False),
    ("map_book.pdf", "pdf", False, False), ("plant.ifc", "ifc", False, False), ("results.bin", "epanet_bin", False, False),
]
TIE_BREAKS: list[tuple[str, str]] = [("ole_header_only.xls", "xls_legacy"), ("ole_header_only.dgn", "dgn_v8")]
ENCODING_CASES = ("inp_with_bom.inp", "net3_cp1252_title.inp", "inp_utf16.inp")
PASS_THROUGH_CASES = ("net3_crlf.inp", "net3_lf.inp", "Net3.inp")
DUPLICATE_CASES = ("net3_dup_junction.inp", "net3_junction_tank_collision.inp", "net3_dup_pipe.inp",
                   "net3_long_id_32.inp")
REFUSED_ZIPS = ("zip_slip.zip", "zip_ratio_bomb.zip", "zip_absolute.zip")
# The review's cases, added after the pre-registered tables above (which are unchanged), each in the check it belongs
# to; the final verification's cases follow each table's first rows.  REVIEW_SNIFF: (fixture, kind it must come back
# as or None when only the time bound applies, time bound in seconds or None when untimed); every row must also come
# back without an exception and with by_extension false, and a timed row within its bound.  The bound is a row's own:
# TIME_BOUND_S for the inputs that ran for many seconds to hours before their fix, CSV_TIME_BOUND_S for the two CSV
# texts, whose slow version took 2 to 7 s (so a 2 s bound could miss a return of the old code on a fast machine) and
# whose fixed version takes hundredths of a second.  Only within the bound or not is recorded, never seconds.
TIME_BOUND_S = 2.0
CSV_TIME_BOUND_S = 0.25
REVIEW_SNIFF: list[tuple[str, str | None, float | None]] = [
    ("swmm_sewer.inp", "swmm_inp", None), ("chlorine_decay.msx", "epanet_msx", None),
    ("net3_cr_only.inp", "epanet_inp", None), ("net3_cr_cr_lf.inp", "epanet_inp", None),
    ("inp_title_over_64kb.inp", "epanet_inp", None), ("inp_patterns_first.inp", "epanet_inp", None),
    ("inp_comment_header_70kb.inp", "epanet_inp", None),
    ("finder_model.zip", "epanet_inp_zip", None), ("finder_shapefile.zip", "shapefile_zip", None),
    ("city_mains_shp.zip", "shapefile_zip", None), ("model_and_gis.zip", "epanet_inp_zip", None),
    ("geojson_crs_string.geojson", "geojson", None), ("geojson_crs_name_number.geojson", "geojson", None),
    ("geojson_features_not_list.geojson", "geojson", None), ("esri_sr_string.json", "esri_json", None),
    ("geojson_deep_geometry.geojson", "geojson", None), ("json_deep.json", "json_other", None),
    ("kml_doctype_subset.kml", "kml", None),
    ("xml_comment_run.xml", "text_other", TIME_BOUND_S), ("xml_comment_run_long.xml", "text_other", TIME_BOUND_S),
    ("xml_pi_run.xml", "text_other", TIME_BOUND_S), ("slow_view.gpkg", "geopackage", TIME_BOUND_S),
    ("csv_quotes.csv", None, CSV_TIME_BOUND_S),
    # the final verification: headers WNTR reads (singular names, a form feed first), a 263-byte quote run
    ("net3_singular_headers.inp", "epanet_inp", None), ("net3_form_feed_headers.inp", "epanet_inp", None),
    ("net1_nul_padded.inp", "epanet_inp", None), ("csv_quote_run.csv", None, CSV_TIME_BOUND_S),
]
REVIEW_PASS_THROUGH = ("net3_cr_only.inp", "net3_cr_cr_lf.inp", "inp_title_over_64kb.inp", "inp_patterns_first.inp",
                       "inp_comment_header_70kb.inp",
                       "net3_singular_headers.inp", "net3_form_feed_headers.inp",      # the final verification
                       "net1_nul_padded.inp")
PASS_THROUGH_NETWORK = {"net1_nul_padded.inp": "Net1"}       # the network WNTR must read; Net3 for every other row
REVIEW_ENCODING = ("net3_cp1252_undefined_byte.inp",)
REVIEW_ID_BLOCKS = ("net1_accented_id_cp1252.inp", "net3_long_pattern_id.inp",
                    "net3_singular_pipe_dup.inp")                                        # the final verification
REVIEW_WNTR_BLOCK = "net3_bom_missing_node.inp"
REVIEW_TEST_RUN_BLOCK = "net3_isolated_junction.inp"                                     # the final verification
REVIEW_ZIPS = ("finder_model.zip", "model_and_gis.zip", "finder_shapefile.zip", "city_mains_shp.zip")
REVIEW_BROKEN_ZIPS = ("zip_bad_offset.zip",)                                             # the final verification
# task 15b: the physics fixtures (residualmap/intake_exp.py builds them for outputs/intake/task15b_physics.json), and
# task 15's four leftovers, fixed in task 15b
PHYSICS_FIXTURES = ("net3_bulk_lines.inp", "net3_dw.inp", "net3_cm.inp", "net3_clock6.inp", "net3_weekly_pattern.inp",
                    "ky10_const_elev.inp", "net3_no_coords.inp", "net3_no_controls.inp", "net3_speed_pattern_pumps.inp",
                    "net3_zero_demand.inp", "net3_two_parts.inp", "net3_no_source.inp", "watergems_export_like.inp",
                    "swmm_after_64kb.inp", "zip_model_as_txt.zip", "utf8_stray_cp1252.inp",
                    # task 15b's review cases
                    "net3_dw_uniform.inp", "net3_rule_premise.inp", "net3_cp1252_pair.inp",
                    # task 15b's second review cases
                    "net3_cp1252_one_pair_id.inp", "net3_cp1252_dash_ids.inp", "utf8_more_stray_than_accents.inp",
                    "utf8_one_id_six_stray.inp", "utf8_pasted_block.inp", "utf8_pasted_chance_pair.inp",
                    "net3_clock_235950.inp", "net3_clock_000020.inp", "net3_dw_plastic.inp", "net1_export_like.inp")
# task 16: the grab-log fixtures (residualmap/intake_exp.py builds them for outputs/intake/task16_samples.json), with
# net3_clock6.inp, log_semicolon.csv and log_tab.txt from the tables above
SAMPLE_FIXTURES = ("task16_tap_map.csv", "lims_one_datetime_us.csv", "lims_censored.csv", "ambiguous_dates.csv",
                   "grab_log_cp1252.csv", "grab_log_bom.csv", "grab_log_typed.xlsx", "grab_log_typed.csv",
                   "analyzer_15min.csv", "utc_shifted_log.csv", "unknown_tap_log.csv", "mixed_species_log.csv",
                   "multi_month_log.csv", "clock_log_net3.csv", "clock_log_clock6.csv", "log_semicolon.csv",
                   "log_tab.txt", "net3_clock6.inp", "mains.gpkg", "net3_backup.inp.bak",
                   # task 16's review
                   "junction_list.csv", "log_title_line.csv", "log_lab_preamble.csv", "log_sep_hint.csv",
                   "log_chloride_first.csv", "log_month_names.csv", "log_multiline_comment.csv", "log_overrange.csv",
                   "log_ugL.csv", "log_implausible.csv", "lims_long.csv", "analyzer_midnight.xlsx",
                   "analyzer_midnight.csv", "grab_study_24h.csv", "analyzer_interleaved.csv",
                   "log_cover_and_stray.xlsx",
                   # task 16's final verification
                   "long_paren.csv", "long_short_names.csv", "long_lims_codes.csv", "long_unknown_codes.csv",
                   "long_unknown_codes_param.csv", "long_two_unnamed.csv", "wide_paren.csv",
                   "wide_paren_total_first.csv", "wide_residual_and_free.csv", "chloramine_header.csv",
                   "ugl_unnamed.csv", "log_comma_title_semicolon.csv", "log_two_cell_title.xlsx", "excel_date_only.xlsx",
                   "excel_date_only.csv", "time_column_epoch_date.csv", "superscript_tap_log.csv",
                   "superscript_tap_map.csv", "junctions_2col.csv", "qualifier_interleaved.csv",
                   "unknown_taps_many.csv", "csv_tiny_rows.csv", "net2_clock_log.csv", "net2_tap_map.csv",
                   *(f"ugl_header_{k}.csv" for k in range(len(UGL_HEADERS))),
                   *(f"units_column_{k}.csv" for k in range(len(UNIT_SPELLINGS))))
# the second review's encoding cases: (fixture, how it must decode, the IDs it must keep)
ENCODING_CONTROLS = (("net3_cp1252_one_pair_id.inp", "cp1252", (ONE_PAIR_ID,)),
                     ("net3_cp1252_dash_ids.inp", "cp1252", tuple(n for _, n in DASH_IDS)),
                     ("utf8_more_stray_than_accents.inp", "mixed", tuple(n for _, n in MORE_STRAY_RENAMES)),
                     ("utf8_one_id_six_stray.inp", "mixed", (ONE_CAPITAL_ID[1],)),
                     ("utf8_pasted_block.inp", "mixed", tuple(n for _, n in PASTED_BLOCK_RENAMES)),
                     ("utf8_pasted_chance_pair.inp", "mixed", tuple(n for _, n in ACCENTED_RENAMES)))


# ----------------------------------------------------------------------------- task 16c: zone inlets
# Written here by hand (string formatting), never by importers.read_inlet_residuals or read_grab_log.  Net3's two
# entry points are its reservoirs River and Lake.  zone_tap_map.csv is task16_tap_map.csv with a role column: AN-1, one
# of analyzer_15min.csv's two online analyzers, becomes Lake's inlet analyzer (its junction column names Lake); AN-2
# stays an analyzer in the network.  The residuals are invented.
ZONE_INLET_TAP = ("AN-1", "Lake")      # (tap, the inlet it measures) in zone_tap_map.csv
ZONE_RIVER_MGL = 1.20                  # River's residual in zone_inlet_residuals.csv (Lake's comes from AN-1)
ZONE_RATIOS_FILE = (("River", 1.20), ("Lake", 0.60))   # zone_inlets_net3.csv: Lake at half River's residual


def analyzer_values(site_index: int) -> list[float]:
    """The readings analyzer_15min() writes for one of its sites, as written (3 decimals): the expected inlet median."""
    return [float(f"{0.55 + 0.1 * site_index + 0.08 * math.sin(2 * math.pi * step / 1440) - 0.002 * day:.3f}")
            for day in range(1, ANALYZER_DAYS + 1) for step in range(0, 24 * 60, ANALYZER_STEP_MIN)]


def zone_tap_map() -> bytes:
    rows = [(t, ZONE_INLET_TAP[1], "inlet", "zone inlet analyzer, invented") if t == ZONE_INLET_TAP[0]
            else (t, j, "route" if t.startswith("T") else "", "invented") for t, j in T16_TAPS]
    return _csv_text("tap_id,junction_id,role,notes", rows).encode("utf-8")


def zone_inlet_residuals() -> bytes:
    return _csv_text("inlet_node_id,residual_mgL,notes", [("River", f"{ZONE_RIVER_MGL:.2f}", "invented")]).encode()


def zone_inlets_net3() -> bytes:
    return _csv_text("inlet_node_id,residual_mgL", [(n, f"{v:.2f}") for n, v in ZONE_RATIOS_FILE]).encode()


def zone_inlets_ugl() -> bytes:
    """The same two residuals in µg/L, under headers the reader has to guess ('Inlet', 'Residual (ug/L)'), with a
    semicolon separator and decimal commas."""
    return _csv_text("Inlet;Residual (ug/L)", [(n, f"{v * 1000:.1f}".replace(".", ",")) for n, v in ZONE_RATIOS_FILE],
                     sep=";").encode()


EXAMPLE_INLET_TAP = ("T10", "Lake")    # example_tap_map_inlet.csv: docs/example_tap_map.csv with T10 as Lake's inlet tap


def example_tap_map_inlet() -> bytes:
    """docs/example_tap_map.csv (read with the csv module, not the reader under test) with a role column: T10, a tap
    of the committed example log sampled every month, becomes Lake's inlet tap, so a six-month log has an inlet tap."""
    import csv as _csv
    with open(os.path.join(REPO, "docs", "example_tap_map.csv"), newline="", encoding="utf-8") as fh:
        rows = list(_csv.reader(fh))
    out = [(r[0], EXAMPLE_INLET_TAP[1] if r[0] == EXAMPLE_INLET_TAP[0] else r[1],
            "inlet" if r[0] == EXAMPLE_INLET_TAP[0] else "") for r in rows[1:]]
    return _csv_text("tap_id,junction_id,role", out).encode()


# the review of task 16c: a declared inlet tap whose readings are all before the grab rows' window (three months of
# grab samples at two taps, April to June, and three inlet readings at INL, January to March), the same log with the
# inlet's only readings below detection, role cells that mention an inlet but mark none, a tap marked as an inlet at a
# junction, and inlet files the reader must refuse (a µg/L value written with a thousands separator, a total chlorine
# column in free-chlorine mode, a residual under a hundredth of the highest)
ZONE_WINDOW_TAPS = (("T1", "10"), ("T2", "15"))
ZONE_WINDOW_INLET = ("INL", "Lake")
ZONE_WINDOW_INLET_MGL = (0.90, 0.80, 0.85)      # INL's readings, January to March (median 0.85)


def _zone_window_rows(inlet_values) -> list:
    rows = [(f"2026-{m:02d}-1{k}", "09:00", t, f"{0.50 + 0.05 * k:.2f}")
            for m in (4, 5, 6) for k, (t, _) in enumerate(ZONE_WINDOW_TAPS)]
    rows += [(f"2026-{m:02d}-10", "08:00", ZONE_WINDOW_INLET[0], v) for m, v in zip((1, 2, 3), inlet_values)]
    return sorted(rows)


def zone_inlet_window_log() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL",
                     _zone_window_rows([f"{v:.2f}" for v in ZONE_WINDOW_INLET_MGL])).encode()


def zone_inlet_censored_log() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL", _zone_window_rows(["<0.02", "ND", "<0.02"])).encode()


def zone_window_tap_map() -> bytes:
    return _csv_text("tap_id,junction_id,role", [(t, j, "route") for t, j in ZONE_WINDOW_TAPS]
                     + [(*ZONE_WINDOW_INLET, "inlet")]).encode()


ZONE_NEAR_MISS_ROLES = ("not inlet", "non-inlet", "downstream of inlet")


def zone_tap_map_not_inlet() -> bytes:
    """ZONE_WINDOW_TAPS and INL, with role cells that mention an inlet but mark none (and INL's own junction a
    junction, 20, so a reader that took any cell with the word as an inlet would refuse it as no entry point)."""
    return _csv_text("tap_id,junction_id,role", [("T1", "10", ZONE_NEAR_MISS_ROLES[0]), ("T2", "15", ZONE_NEAR_MISS_ROLES[1]),
                                                 ("INL", "20", ZONE_NEAR_MISS_ROLES[2])]).encode()


def zone_tap_map_wrong_node() -> bytes:
    """INL marked as a zone inlet, but its junction column names junction 20, which is no entry point of Net3."""
    return _csv_text("tap_id,junction_id,role", [(t, j, "route") for t, j in ZONE_WINDOW_TAPS]
                     + [("INL", "20", "Zone inlet")]).encode()


# task 17's fix of task 16c's leftovers (its final verification): a tap marked as an inlet with no node ID; a µg/L grab
# log whose inlet tap's values are 1,000 µg/L or more, written without and with a thousands separator; inlet files
# whose second column is no residual (pH, NH2Cl, a monochloramine residual); an inlet tap read once a day with a date
# and no time of day, beside one row whose date cannot be read, and the same inlet with only unreadable dates; and a
# fixed route of ten taps for the pilot's leave-sites-out refusal
ZONE_UGL_INLET = (1180, 1200, 1190)    # INL's readings in µg/L, April to June (median 1.19 mg/L)
ZONE_DAILY_INLET_MGL = (0.84, 0.86, 0.85, 0.83, 0.87)   # INL's daily values with no time of day (median 0.85)
FIXED_ROUTE_TAPS = 10


def zone_tap_map_inlet_no_node() -> bytes:
    return _csv_text("tap_id,junction_id,role", [(t, j, "route") for t, j in ZONE_WINDOW_TAPS]
                     + [(ZONE_WINDOW_INLET[0], "", "inlet")]).encode()


def _ugl_rows(inlet_cells) -> list:
    rows = [(f"2026-{m:02d}-1{k}", "09:00", t, f"{500 + 50 * k}") for m in (4, 5, 6)
            for k, (t, _) in enumerate(ZONE_WINDOW_TAPS)]
    rows += [(f"2026-{m:02d}-10", "08:00", ZONE_WINDOW_INLET[0], v) for m, v in zip((4, 5, 6), inlet_cells)]
    return sorted(rows)


def zone_inlet_ugl_log() -> bytes:
    return _csv_text("date,time,tap_id,Free Chlorine (ug/L)", _ugl_rows([str(v) for v in ZONE_UGL_INLET])).encode()


def zone_inlet_ugl_thousands_log() -> bytes:
    """The same, with the inlet's values written with a thousands separator ('1,180', as Excel's #,##0 writes them),
    quoted so the separator stays in the cell."""
    rows = _ugl_rows([f"{v:,}" for v in ZONE_UGL_INLET])
    text = "date,time,tap_id,Free Chlorine (ug/L)\n" + "".join(
        ",".join(f'"{c}"' if "," in c else c for c in r) + "\n" for r in rows)
    return text.encode()


def _daily_rows(inlet_dates) -> list:
    rows = [(f"2026-{m:02d}-1{k}", "09:00", t, f"{0.50 + 0.05 * k:.2f}")
            for m in (4, 5, 6) for k, (t, _) in enumerate(ZONE_WINDOW_TAPS)]
    rows += [(d, "", ZONE_WINDOW_INLET[0], f"{v:.2f}") for d, v in zip(inlet_dates, ZONE_DAILY_INLET_MGL)]
    return rows


def zone_inlet_daily_log() -> bytes:
    dates = ["2026-06-01", "2026-06-02", "2026-06-03", "2026-06-04", "2026-06-05"]
    rows = _daily_rows(dates) + [("2026-06-31", "", ZONE_WINDOW_INLET[0], "0.90")]
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def zone_inlet_bad_dates_log() -> bytes:
    return _csv_text("date,time,tap_id,free_chlorine_mgL",
                     _daily_rows(["2026-06-31", "2026-13-02", "June", "2026-02-30", "?"])).encode()


def fixed_route_log() -> bytes:
    """Ten taps sampled once a month for six months (January to June 2026), invented readings: a fixed compliance
    route with no rotating taps, for the pilot's --holdout-sites refusal."""
    rows = [(f"2026-{m:02d}-{5 + k:02d}", "09:30", f"F{k + 1}", f"{0.40 + 0.03 * k + 0.01 * m:.2f}")
            for m in range(1, 7) for k in range(FIXED_ROUTE_TAPS)]
    return _csv_text("date,time,tap_id,free_chlorine_mgL", rows).encode()


def fixed_route_tap_map() -> bytes:
    js = _wn("Net3").junction_name_list[:FIXED_ROUTE_TAPS]
    return _csv_text("tap_id,junction_id", [(f"F{k + 1}", j) for k, j in enumerate(js)]).encode()


LEFTOVER_16C_FIXTURES = {
    "zone_tap_map_inlet_no_node.csv": zone_tap_map_inlet_no_node, "zone_inlet_ugl_log.csv": zone_inlet_ugl_log,
    "zone_inlet_ugl_thousands_log.csv": zone_inlet_ugl_thousands_log,
    "zone_inlets_ph.csv": lambda: _csv_text("inlet,pH", [("River", "7.8"), ("Lake", "7.6")]).encode(),
    "zone_inlets_nh2cl.csv": lambda: _csv_text("inlet_node_id,NH2Cl", [("River", "1.20"), ("Lake", "0.60")]).encode(),
    "zone_inlets_mono_residual.csv": lambda: _csv_text("inlet_node_id,Monochloramine residual",
                                                       [("River", "1.20"), ("Lake", "0.60")]).encode(),
    "zone_inlet_daily_log.csv": zone_inlet_daily_log, "zone_inlet_bad_dates_log.csv": zone_inlet_bad_dates_log,
    "fixed_route_log.csv": fixed_route_log, "fixed_route_tap_map.csv": fixed_route_tap_map}


ZONE_FIXTURES = ("analyzer_15min.csv", "task16_tap_map.csv", "zone_tap_map.csv", "zone_inlet_residuals.csv",
                 "zone_inlets_net3.csv", "zone_inlets_ugl.csv", "zone_inlets_missing.csv", "zone_inlets_unknown.csv",
                 "zone_inlets_censored.csv", "example_tap_map_inlet.csv",
                 "zone_inlet_window_log.csv", "zone_inlet_censored_log.csv", "zone_window_tap_map.csv",
                 "zone_tap_map_not_inlet.csv", "zone_tap_map_wrong_node.csv", "zone_inlets_thousands_ugl.csv",
                 "zone_inlets_total.csv", "zone_inlets_tiny.csv")
FIXTURES.update({"example_tap_map_inlet.csv": example_tap_map_inlet,
    "zone_tap_map.csv": zone_tap_map, "zone_inlet_residuals.csv": zone_inlet_residuals,
    "zone_inlets_net3.csv": zone_inlets_net3, "zone_inlets_ugl.csv": zone_inlets_ugl,
    "zone_inlets_missing.csv": lambda: _csv_text("inlet_node_id,residual_mgL", [("River", "1.20")]).encode(),
    "zone_inlets_unknown.csv": lambda: _csv_text("inlet_node_id,residual_mgL",
                                                 [("River", "1.20"), ("Lake", "0.60"), ("Tank9", "0.50")]).encode(),
    "zone_inlets_censored.csv": lambda: _csv_text("inlet_node_id,residual_mgL",
                                                  [("River", "1.20"), ("Lake", "<0.02")]).encode(),
    "zone_inlet_window_log.csv": zone_inlet_window_log, "zone_inlet_censored_log.csv": zone_inlet_censored_log,
    "zone_window_tap_map.csv": zone_window_tap_map, "zone_tap_map_not_inlet.csv": zone_tap_map_not_inlet,
    "zone_tap_map_wrong_node.csv": zone_tap_map_wrong_node,
    "zone_inlets_thousands_ugl.csv": lambda: _csv_text("inlet_node_id,residual_ugL",
                                                       [("River", "1,200"), ("Lake", "850")]).encode(),
    "zone_inlets_total.csv": lambda: _csv_text("inlet_node_id,total_chlorine_mgL",
                                               [("River", "1.20"), ("Lake", "0.60")]).encode(),
    "zone_inlets_tiny.csv": lambda: _csv_text("inlet_node_id,residual_mgL",
                                              [("River", "1.20"), ("Lake", "0.005")]).encode()})
FIXTURES.update(LEFTOVER_16C_FIXTURES)
ZONE_FIXTURES = ZONE_FIXTURES + tuple(LEFTOVER_16C_FIXTURES)


# ----------------------------------------------------------------------------- task 17: GIS layers
# Every GIS fixture is written here by hand (the shapefile and dBASE writers above and below, a WKB writer, the
# GeoPackage tables, GeoJSON and KML text), except ky4_mains_pyshp.zip, which pyshp writes and pyshp reads (same-library
# self-consistency, labelled so).  The ky4 mains layer is ky4.inp's 1,156 pipes as lines only (no pump, valve, tank or
# reservoir), each pipe's PIPEDIAM in one of three styles, every line end moved by up to +-0.3 ft on each axis, 19 ends
# pulled back 4 ft short of their node, and 43 pairs of pipes of one size drawn as one line through the tee where a
# third pipe ends (GIS_RUNS, GIS_UNDERSHOOTS, GIS_JITTER_FT; seeded, chosen before any scored run).  Longitude and
# latitude come from the ellipsoidal Lambert conformal conic inverse above (task 16b's), with ky4's feet laid out
# around LCC_CENTRE: the reader's east-north plane shares no code with it.  No fixture is a real place's map.
GIS_JITTER_FT = 0.3                    # every line end: uniform in +-0.3 ft on each axis
GIS_UNDERSHOOTS, GIS_UNDERSHOOT_FT = 19, 4.0
GIS_RUNS = 43
GIS_SEED = 17
GIS_DIAMETER_STYLES = ('{d:g}"', "{d:g} in", "{mm} mm")    # record k uses style k % 3
KML_ERROR_M = 5.0                      # ky4_googleearth.kmz: radial RMS error of each line end (per axis r / sqrt(2))
KML_SEED = 18
KML_TOLERANCES_M = (5.0, 10.0, 15.0, 20.0)
PRJ_M = ('PROJCS["NAD_1983_UTM_Example_Zone",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",'
         'SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
         'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["Central_Meridian",-99.0],'
         'PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]')
GPKG_SRS_ID = 100001                   # the GeoPackage fixtures' coordinate system, defined by PRJ_FT's text


@functools.lru_cache(maxsize=1)
def _ky4_parts():
    """ky4's pipes as polylines in its own feet: {pipe: (start node, end node, [(x, y), ...], diameter in, length ft)},
    its node positions and elevations (ft), and each node's pipes."""
    wn = _wn("ky4")
    xy = {n: tuple(map(float, wn.get_node(n).coordinates)) for n in wn.node_name_list}
    elev = {}
    for n in wn.node_name_list:
        nd = wn.get_node(n)
        elev[n] = float(getattr(nd, "elevation", None) or getattr(nd, "base_head", 0.0) or 0.0) / 0.3048
    pipes, inc = {}, {}
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        pts = [xy[p.start_node_name]] + [tuple(map(float, v)) for v in p.vertices] + [xy[p.end_node_name]]
        pipes[pn] = (p.start_node_name, p.end_node_name, pts, round(p.diameter / 0.0254, 3), p.length / 0.3048)
        inc.setdefault(p.start_node_name, []).append(pn)
        inc.setdefault(p.end_node_name, []).append(pn)
    return pipes, xy, elev, inc


def _seg_len(pts) -> list[float]:
    return [math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:])]


@functools.lru_cache(maxsize=1)
def ky4_mains_records() -> tuple:
    """The ky4 mains layer: a tuple of records, each {'pts': [(x, y) ft], 'pipes': (pipe, ...), 'diameter_in',
    'ends': (node at the first point, node at the last), 'undershoot': end index or None}, in ky4's pipe order (a run at
    its first pipe's place).  Seeded (GIS_SEED); the runs, undershoots and jitter are chosen without looking at what
    the importer does with them."""
    import numpy as np
    pipes, xy, _, inc = _ky4_parts()
    rng = np.random.default_rng(GIS_SEED)
    tees = sorted(n for n, ps in inc.items() if len(ps) == 3)
    rng.shuffle(tees)
    used, runs = set(), {}
    for t in tees:
        if sum(v is not None for v in runs.values()) >= GIS_RUNS:
            break
        ps = sorted(inc[t])
        pairs = [(a, b) for i, a in enumerate(ps) for b in ps[i + 1:]
                 if pipes[a][3] == pipes[b][3] and a not in used and b not in used]
        pairs = [(a, b) for a, b in pairs if ({pipes[a][0], pipes[a][1]} - {t}) != ({pipes[b][0], pipes[b][1]} - {t})
                 and len({pipes[a][0], pipes[a][1]}) == 2 and len({pipes[b][0], pipes[b][1]}) == 2]
        if not pairs:
            continue
        a, b = pairs[int(rng.integers(len(pairs)))]
        used |= {a, b}
        runs[a] = (a, b, t)
        runs[b] = None
    recs = []
    for pn, (s, e, pts, d, _) in pipes.items():
        if pn in runs:
            if runs[pn] is None:
                continue
            a, b, t = runs[pn]
            pa = pipes[a][2] if pipes[a][1] == t else pipes[a][2][::-1]
            pb = pipes[b][2] if pipes[b][0] == t else pipes[b][2][::-1]
            ea = pipes[a][0] if pipes[a][1] == t else pipes[a][1]
            eb = pipes[b][1] if pipes[b][0] == t else pipes[b][0]
            recs.append({"pts": list(pa) + list(pb[1:]), "pipes": (a, b), "diameter_in": d, "ends": (ea, eb),
                         "tee": t, "undershoot": None})
        else:
            recs.append({"pts": list(pts), "pipes": (pn,), "diameter_in": d, "ends": (s, e), "tee": None,
                         "undershoot": None})
    # undershoots: single pipes, an end at a node of 2 or more pipes, the end's segment longer than twice the pullback
    cand = []
    for k, r in enumerate(recs):
        if len(r["pipes"]) != 1:
            continue
        for end in (0, 1):
            node = r["ends"][end]
            seg = _seg_len(r["pts"])[0 if end == 0 else -1]
            if len(inc[node]) >= 2 and seg > 2 * GIS_UNDERSHOOT_FT:
                cand.append((k, end))
    pick = rng.choice(len(cand), size=GIS_UNDERSHOOTS, replace=False)
    taken = set()
    for i in sorted(pick):
        k, end = cand[i]
        if k in taken:
            continue
        taken.add(k)
        pts = recs[k]["pts"]
        a, b = (pts[0], pts[1]) if end == 0 else (pts[-1], pts[-2])
        L = math.hypot(b[0] - a[0], b[1] - a[1])
        moved = (a[0] + (b[0] - a[0]) * GIS_UNDERSHOOT_FT / L, a[1] + (b[1] - a[1]) * GIS_UNDERSHOOT_FT / L)
        if end == 0:
            pts[0] = moved
        else:
            pts[-1] = moved
        recs[k]["undershoot"] = end
    # jitter every line end
    for r in recs:
        for end in (0, -1):
            dx, dy = rng.uniform(-GIS_JITTER_FT, GIS_JITTER_FT, 2)
            x, y = r["pts"][end]
            r["pts"][end] = (x + float(dx), y + float(dy))
    for k, r in enumerate(recs):
        d = r["diameter_in"]
        r["text"] = GIS_DIAMETER_STYLES[k % 3].format(d=d, mm=round(d * 25.4))
    return tuple(recs)


def ky4_manifest() -> dict:
    """What the scorer needs: the original pipes (start, end, diameter in), and for every record its pipes and end
    nodes, in the layer's order (records with no ky4 pipe, the pyshp fixture's extras, are not listed)."""
    pipes = _ky4_parts()[0]
    return {"pipes": {p: (v[0], v[1], v[3]) for p, v in pipes.items()},
            "records": [{"pipes": r["pipes"], "ends": r["ends"], "tee": r["tee"]} for r in ky4_mains_records()]}


def _shp_lines(records, shape_type: int = 3) -> tuple[bytes, bytes]:
    """.shp and .shx of polylines with several parts (records of [part, ...], each part [(x, y), ...]; None for a
    null shape), from the ESRI shapefile specification, written by hand."""
    recs = []
    xs, ys = [], []
    for parts in records:
        if parts is None:
            recs.append(struct.pack("<i", 0))
            continue
        pts = [p for part in parts for p in part]
        px, py = [p[0] for p in pts], [p[1] for p in pts]
        xs += px
        ys += py
        starts, k = [], 0
        for part in parts:
            starts.append(k)
            k += len(part)
        recs.append(struct.pack("<i", shape_type) + struct.pack("<4d", min(px), min(py), max(px), max(py))
                    + struct.pack("<ii", len(parts), len(pts)) + struct.pack(f"<{len(parts)}i", *starts)
                    + b"".join(struct.pack("<2d", p[0], p[1]) for p in pts))

    def header(length_words):
        return (struct.pack(">7i", 9994, 0, 0, 0, 0, 0, length_words) + struct.pack("<ii", 1000, shape_type)
                + struct.pack("<8d", min(xs), min(ys), max(xs), max(ys), 0, 0, 0, 0))

    body, index, offset = [], [], 50
    for i, r in enumerate(recs):
        body.append(struct.pack(">ii", i + 1, len(r) // 2) + r)
        index.append(struct.pack(">ii", offset, len(r) // 2))
        offset += 4 + len(r) // 2
    body, index = b"".join(body), b"".join(index)
    return header(50 + len(body) // 2) + body, header(50 + len(index) // 2) + index


def _dbf_text(fields: list[tuple[str, str, int]], rows: list[list], encoding: str = "latin-1",
              deleted: tuple = ()) -> bytes:
    """A dBASE III table whose text is in encoding, rows in `deleted` marked deleted ('*')."""
    rec_len = 1 + sum(w for _, _, w in fields)
    head = struct.pack("<BBBBIHH20x", 3, 126, 1, 1, len(rows), 32 + 32 * len(fields) + 1, rec_len)
    for name, typ, w in fields:
        head += name.encode("ascii")[:10].ljust(11, b"\x00") + typ.encode() + b"\x00" * 4 + bytes([w, 0]) + b"\x00" * 14
    head += b"\r"
    out = []
    for i, r in enumerate(rows):
        cells = b"".join((str(v).encode(encoding)[:w].ljust(w) if t == "C" else str(v).encode()[:w].rjust(w))
                         for v, (_, t, w) in zip(r, fields))
        out.append((b"*" if i in deleted else b" ") + cells)
    return head + b"".join(out) + b"\x1a"


def _shp_zip(stem: str, records, rows, fields, prj: str | None = PRJ_FT, cpg: str | None = "UTF-8",
             encoding: str = "utf-8", deleted: tuple = ()) -> bytes:
    shp, shx = _shp_lines(records)
    members = [(f"{stem}.shp", shp), (f"{stem}.shx", shx), (f"{stem}.dbf", _dbf_text(fields, rows, encoding, deleted))]
    if prj is not None:
        members.append((f"{stem}.prj", prj.encode("ascii")))
    if cpg is not None:
        members.append((f"{stem}.cpg", cpg.encode("ascii")))
    return _zip(members)


def ky4_mains_stdlib_shp_zip() -> bytes:
    recs = ky4_mains_records()
    return _shp_zip("ky4_mains", [[r["pts"]] for r in recs], [[r["text"]] for r in recs], [("PIPEDIAM", "C", 12)])


KY4_PYSHP_NULLS = (100, 500, 900)       # ky4_mains_pyshp.zip: null shapes written before these records
KY4_PYSHP_DELETED_AT = 700              # ... and a deleted record (a stray line between two far ends) before this one


def ky4_pyshp_index(j: int) -> int:
    """The position in ky4_mains_pyshp.zip's .shp of the layer's record j (the extra records come before it)."""
    return j + sum(1 for n in KY4_PYSHP_NULLS if n <= j) + (1 if j >= KY4_PYSHP_DELETED_AT else 0)


def ky4_mains_pyshp_zip() -> bytes:
    """The ky4 mains layer written by pyshp (the library that reads it: same-library self-consistency), with 3 null
    shapes and 1 record whose DBF deletion flag is set after writing, between ky4's records."""
    import shapefile
    recs = ky4_mains_records()
    shp, shx, dbf = io.BytesIO(), io.BytesIO(), io.BytesIO()
    w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf, shapeType=shapefile.POLYLINE)
    w.field("PIPEDIAM", "C", size=12)
    n, deleted_index = 0, None
    far = (recs[0]["pts"][0], recs[-1]["pts"][-1])
    for k, r in enumerate(recs):
        if k in KY4_PYSHP_NULLS:
            w.null()
            w.record('8"')
            n += 1
        if k == KY4_PYSHP_DELETED_AT:
            w.line([[list(far[0]), list(far[1])]])
            w.record('12"')
            deleted_index = n
            n += 1
        w.line([[list(p) for p in r["pts"]]])
        w.record(r["text"])
        n += 1
    w.close()
    d = bytearray(dbf.getvalue())
    head_len, rec_len = struct.unpack("<HH", bytes(d[8:12]))
    d[head_len + deleted_index * rec_len] = ord("*")
    return _zip([("ky4_mains.shp", shp.getvalue()), ("ky4_mains.shx", shx.getvalue()), ("ky4_mains.dbf", bytes(d)),
                 ("ky4_mains.prj", PRJ_FT.encode()), ("ky4_mains.cpg", b"UTF-8")])


def _wkb_line(parts, dims: str = "", multi: bool = False) -> bytes:
    """ISO WKB: a LineString (one part) or a MultiLineString; dims '', 'Z', 'M' or 'ZM' (Z and M filled from the
    third and fourth value of each point, or 0)."""
    off = {"": 0, "Z": 1000, "M": 2000, "ZM": 3000}[dims]
    nd = 2 + len(dims)

    def one(pts):
        body = struct.pack("<BII", 1, 2 + off, len(pts))
        for p in pts:
            v = list(p[:2]) + [float(p[2]) if len(p) > 2 else 0.0, float(p[3]) if len(p) > 3 else 0.0][:nd - 2]
            body += struct.pack(f"<{nd}d", *v)
        return body
    if not multi:
        return one(parts[0])
    return struct.pack("<BII", 1, 5 + off, len(parts)) + b"".join(one(p) for p in parts)


def _gpkg_blob(wkb: bytes, srs_id: int) -> bytes:
    return b"GP" + bytes([0, 0b00000001]) + struct.pack("<i", srs_id) + wkb


def _gpkg_file(layers: list[tuple[str, str, list[tuple[bytes, list]]]], cols: list[tuple[str, str]],
               srs: tuple = (GPKG_SRS_ID, "NONE", GPKG_SRS_ID, PRJ_FT), z: int = 0, m: int = 0) -> bytes:
    """A GeoPackage of feature tables: layers [(table, geometry type, [(blob, [values]), ...])], every table with the
    columns cols [(name, type)], in the coordinate system srs (srs_id, organization, id, definition)."""
    sid, org, oid, defn = srs
    st = [("PRAGMA application_id = 1196444487", ()), ("PRAGMA user_version = 10300", ()),
          ("create table gpkg_spatial_ref_sys (srs_name text, srs_id integer primary key, organization text, "
           "organization_coordsys_id integer, definition text, description text)", ()),
          ("insert into gpkg_spatial_ref_sys values ('Undefined cartesian', -1, 'NONE', -1, 'undefined', '')", ()),
          ("insert into gpkg_spatial_ref_sys values ('Undefined geographic', 0, 'NONE', 0, 'undefined', '')", ())]
    if sid not in (-1, 0):
        st.append(("insert into gpkg_spatial_ref_sys values (?, ?, ?, ?, ?, '')", ("fixture system", sid, org, oid,
                                                                                    defn)))
    st += [("create table gpkg_contents (table_name text primary key, data_type text, identifier text, description "
            "text, last_change text, min_x real, min_y real, max_x real, max_y real, srs_id integer)", ()),
           ("create table gpkg_geometry_columns (table_name text, column_name text, geometry_type_name text, srs_id "
            "integer, z integer, m integer)", ())]
    for layer in layers:
        table, gtype, rows = layer[:3]
        lcols = layer[3] if len(layer) > 3 else cols         # a layer may name its own columns
        extra = "".join(f", {c} {t}" for c, t in lcols)
        st.append(("insert into gpkg_contents values (?, 'features', ?, '', '2026-10-20T00:00:00Z', null, null, null, "
                   "null, ?)", (table, table, sid)))
        st.append(("insert into gpkg_geometry_columns values (?, 'geom', ?, ?, ?, ?)", (table, gtype, sid, z, m)))
        st.append((f"create table {table} (fid integer primary key, geom blob{extra})", ()))
        for k, (blob, vals) in enumerate(rows):
            st.append((f"insert into {table} values (?, ?{', ?' * len(lcols)})", (k + 1, blob, *vals)))
    return _sqlite(st)


def ky4_mains_z_gpkg() -> bytes:
    """The ky4 mains layer as a GeoPackage written by the WKB writer here, with Z (each end at its node's elevation,
    each vertex between them by distance)."""
    _, xy, elev, _ = _ky4_parts()
    rows = []
    for r in ky4_mains_records():
        za, zb = elev[r["ends"][0]], elev[r["ends"][1]]
        L = _seg_len(r["pts"])
        tot = sum(L) or 1.0
        acc, pts = 0.0, []
        for i, p in enumerate(r["pts"]):
            pts.append((p[0], p[1], za + (zb - za) * acc / tot))
            if i < len(L):
                acc += L[i]
        rows.append((_gpkg_blob(_wkb_line([pts], "Z"), GPKG_SRS_ID), [r["text"]]))
    return _gpkg_file([("ky4_mains", "LINESTRING", rows)], [("PIPEDIAM", "text")], z=1)


@functools.lru_cache(maxsize=1)
def _ky4_frame() -> tuple:
    """ky4's feet laid out around LCC_CENTRE: (x0, y0) ky4 feet at the centre, (E0, N0) metres in the LCC zone."""
    _, xy, _, _ = _ky4_parts()
    xs, ys = [p[0] for p in xy.values()], [p[1] for p in xy.values()]
    x0, y0 = round((min(xs) + max(xs)) / 2), round((min(ys) + max(ys)) / 2)
    return (x0, y0), lcc_forward(*LCC_CENTRE)


def ky4_lonlat(p) -> tuple[float, float]:
    """(lon, lat) of a ky4 position in feet: metres in the invented LCC zone around LCC_CENTRE, then lcc_inverse."""
    (x0, y0), (E0, N0) = _ky4_frame()
    lat, lon = lcc_inverse(E0 + (p[0] - x0) * 0.3048, N0 + (p[1] - y0) * 0.3048)
    return lon, lat


def _ky4_lonlat_features() -> list[dict]:
    feats = []
    for r in ky4_mains_records():
        coords = [[round(c, 10) for c in ky4_lonlat(p)] for p in r["pts"]]
        feats.append({"type": "Feature", "properties": {"PIPEDIAM": r["text"]},
                      "geometry": {"type": "LineString", "coordinates": coords}})
    return feats


def ky4_mains_lonlat_geojson() -> bytes:
    return json.dumps({"type": "FeatureCollection", "name": "ky4_mains", "features": _ky4_lonlat_features()}).encode()


KY4_SHUFFLE_SEED = 19


def ky4_mains_shuffled_geojson() -> bytes:
    import numpy as np
    feats = _ky4_lonlat_features()
    order = np.random.default_rng(KY4_SHUFFLE_SEED).permutation(len(feats))
    return json.dumps({"type": "FeatureCollection", "name": "ky4_mains",
                       "features": [feats[i] for i in order]}).encode()


def ky4_googleearth_kmz() -> bytes:
    """ky4's pipes drawn by hand in Google Earth: one placemark per pipe, its size as the name ('8 in'), each line
    end about KML_ERROR_M off (radial RMS), interior vertices as drawn."""
    import numpy as np
    pipes = _ky4_parts()[0]
    rng = np.random.default_rng(KML_SEED)
    marks = []
    for k, (pn, (s, e, pts, d, _)) in enumerate(pipes.items()):
        ll = [list(ky4_lonlat(p)) for p in pts]
        (x0, y0), (E0, N0) = _ky4_frame()
        for end in (0, -1):
            dx, dy = rng.normal(0.0, KML_ERROR_M / math.sqrt(2), 2) / 0.3048
            p = pts[end]
            ll[end] = list(ky4_lonlat((p[0] + dx, p[1] + dy)))
        coords = " ".join(f"{lon:.7f},{lat:.7f},0" for lon, lat in ll)
        marks.append(f"<Placemark><name>{d:g} in</name><LineString><tessellate>1</tessellate><coordinates>{coords}"
                     f"</coordinates></LineString></Placemark>")
    kml = ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
           f"<name>water mains</name><Folder><name>mains</name>{''.join(marks)}</Folder></Document></kml>\n")
    return _zip([("doc.kml", kml.encode("utf-8"))])


def ky4_layer_feet() -> list:
    """The ky4 mains records as plain lines in feet (for the invariance checks, which build a layer in memory)."""
    return [[tuple(p) for p in r["pts"]] for r in ky4_mains_records()]


# the spec quirks (critique 12)
_TWO_PARTS = [[(1000.0, 1000.0), (1100.0, 1000.0)], [(1300.0, 1000.0), (1400.0, 1050.0)]]
QUIRK_DIAMETERS = (4, 6, 8, 10, 12, 16)


def _quirk_lines() -> list:
    return [[(2000.0 + 300 * k, 2000.0), (2000.0 + 300 * k + 200, 2150.0)] for k in range(len(QUIRK_DIAMETERS))]


def multipart_two_parts_gpkg() -> bytes:
    return _gpkg_file([("mains", "MULTILINESTRING", [(_gpkg_blob(_wkb_line(_TWO_PARTS, multi=True), GPKG_SRS_ID),
                                                      ['8"'])])], [("PIPEDIAM", "text")])


def multipart_two_parts_zip() -> bytes:
    return _shp_zip("mains", [_TWO_PARTS], [['8"']], [("PIPEDIAM", "C", 8)])


def null_shape_mid_zip() -> bytes:
    lines = _quirk_lines()
    recs = [[ln] for ln in lines]
    recs.insert(3, None)
    rows = [[f'{d}"'] for d in QUIRK_DIAMETERS]
    rows.insert(3, ['24"'])
    return _shp_zip("mains", recs, rows, [("PIPEDIAM", "C", 8)])


def deleted_dbf_record_zip() -> bytes:
    lines = _quirk_lines()
    recs = [[ln] for ln in lines]
    recs.insert(3, [[(2000.0, 2000.0), (3700.0, 2150.0)]])         # a stray line between two far ends
    rows = [[f'{d}"'] for d in QUIRK_DIAMETERS]
    rows.insert(3, ['24"'])
    return _shp_zip("mains", recs, rows, [("PIPEDIAM", "C", 8)], deleted=(3,))


CPG_NAMES = ("Calle Señora Ávila", "Rue du Château", "Straße Nord")


def cpg_utf8_zip() -> bytes:
    lines = _quirk_lines()[:3]
    return _shp_zip("mains", [[ln] for ln in lines], [[n, '8"'] for n in CPG_NAMES],
                    [("STREET", "C", 30), ("PIPEDIAM", "C", 8)], cpg="UTF-8", encoding="utf-8")


def wkb_z_m_zm_gpkg() -> bytes:
    rows = []
    for k, dims in enumerate(("Z", "M", "ZM")):
        base = [(3000.0 + 300 * k, 3000.0), (3200.0 + 300 * k, 3100.0), (3250.0 + 300 * k, 3300.0)]
        pts = [(x, y, 500.0 + k, 7.5) for x, y in base]
        rows.append((_gpkg_blob(_wkb_line([pts], dims), GPKG_SRS_ID), ['6"']))
    return _gpkg_file([("mains", "LINESTRING", rows)], [("PIPEDIAM", "text")], z=2, m=2)


def wkb_z_m_zm_expected() -> list:
    return [[(3000.0 + 300 * k, 3000.0), (3200.0 + 300 * k, 3100.0), (3250.0 + 300 * k, 3300.0)] for k in range(3)]


def prj_stateplane_ftus_zip() -> bytes:
    return _shp_zip("mains", [[ln] for ln in _quirk_lines()], [[f'{d}"'] for d in QUIRK_DIAMETERS],
                    [("PIPEDIAM", "C", 8)], prj=PRJ_FT)


def gpkg_srs_undefined_gpkg() -> bytes:
    rows = [(_gpkg_blob(_wkb_line([ln]), -1), [f'{d}"']) for ln, d in zip(_quirk_lines(), QUIRK_DIAMETERS)]
    return _gpkg_file([("mains", "LINESTRING", rows)], [("PIPEDIAM", "text")], srs=(-1, "NONE", -1, "undefined"))


def shape_length_mismatch_zip() -> bytes:
    """Coordinates in metres (the .prj says so) and a Shape_Leng field in feet: the declared lengths are 3.28 times
    the geometry's."""
    lines = _quirk_lines()
    rows = [[f'{d}"', f"{math.hypot(ln[1][0] - ln[0][0], ln[1][1] - ln[0][1]) / 0.3048:.3f}"]
            for ln, d in zip(lines, QUIRK_DIAMETERS)]
    return _shp_zip("mains", [[ln] for ln in lines], rows, [("PIPEDIAM", "C", 8), ("Shape_Leng", "N", 16)], prj=PRJ_M)


def quantized_4dp_geojson() -> bytes:
    feats = []
    for r in ky4_mains_records()[:200]:
        coords = [[round(c, 4) for c in ky4_lonlat(p)] for p in r["pts"]]
        feats.append({"type": "Feature", "properties": {"PIPEDIAM": r["text"]},
                      "geometry": {"type": "LineString", "coordinates": coords}})
    return json.dumps({"type": "FeatureCollection", "features": feats}).encode()


# snapping chains (critique 13): ends A, B, C in a row 5 ft apart; a valve cluster of five ends 4 ft apart; a pump
# manifold of four ends 5 ft apart; every end on its own 500 ft line, so no far end meets another (EPSG 2227, US feet)
CHAIN_ENDS = {"chain": [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)],
              "valve_cluster": [(1000.0 + 4 * k, 0.0) for k in range(5)],
              "pump_manifold": [(2000.0, 5.0 * k) for k in range(4)]}


def chain_cluster_geojson() -> bytes:
    """Each end's line leaves its row at right angles, on alternate sides, so no end lies on another line."""
    feats, k = [], 0
    for group, ends in CHAIN_ENDS.items():
        along_x = ends[0][1] == ends[-1][1]
        for i, (x, y) in enumerate(ends):
            side = 1 if i % 2 == 0 else -1
            far = (x, y + side * 500) if along_x else (x + side * 500, y)
            feats.append({"type": "Feature", "properties": {"PIPEID": f"{group}_{i + 1}", "PIPEDIAM": "8 in"},
                          "geometry": {"type": "LineString", "coordinates": [[x, y], list(far)]}})
            k += 1
    return json.dumps({"type": "FeatureCollection", "crs": {"type": "name",
                                                            "properties": {"name": "urn:ogc:def:crs:EPSG::2227"}},
                       "features": feats}).encode()


DIRTY_IDS = ("MAIN 12", "VLV;7", "W-17", "Ünïcode_9", "A" * 40, '"QUOTED"')


def ids_dirty_geojson() -> bytes:
    feats = []
    for k, pid in enumerate(DIRTY_IDS):
        a = (5000.0 + 250 * k, 5000.0)
        b = (5000.0 + 250 * (k + 1), 5000.0)
        feats.append({"type": "Feature", "properties": {"PIPEID": pid, "PIPEDIAM": "6 in"},
                      "geometry": {"type": "LineString", "coordinates": [list(a), list(b)]}})
    return json.dumps({"type": "FeatureCollection", "crs": {"type": "name",
                                                            "properties": {"name": "urn:ogc:def:crs:EPSG::2227"}},
                       "features": feats}).encode()


DIAMETER_CASES = (('8"', "8", "in", ""), ("8 in", "8", "in", ""), ("8in", "8", "in", ""), ("200 mm", "200", "mm", ""),
                  ("200mm PVC", "200", "mm", ""), ("DN150", "150", "mm", ""), ('1965 8" CI', "8", "in", ""),
                  ("L=350' 8\"", "8", "in", ""), ("12x8 RED", "", "", "reducer"), ('24" x 16" reducer', "", "", "reducer"),
                  ("8 & 6", "", "", "two sizes"), ('3/4"', "0.75", "in", ""), ('1-1/2"', "1.5", "in", ""),
                  ("2 1/2 in", "2.5", "in", ""), ("0.2 m", "0.2", "m", ""), ("1,200 mm", "1200", "mm", ""),
                  ("8'' DIP", "8", "in", ""), ("NPS 8", "8", "in", ""), ("6", "6", "", ""), ("6.0", "6", "", ""),
                  ("Ø200", "200", "", ""), ("unknown", "", "", "unreadable"), ("", "", "", "empty"),
                  ('12" DI 1972', "12", "in", ""),
                  # the review of task 17: one size labelled twice, and a coupling of one size
                  ('6" (150mm)', "6", "in", ""), ('150 mm (6")', "6", "in", ""), ('DN 150 (6")', "6", "in", ""),
                  ('6"x6"', "6", "in", ""))
BARE_CASES = (((4, 6, 8, 12), "in"), ((100, 150, 200, 300), "mm"), ((0.1, 0.15, 0.3), "m"), ((1, 8, 100), None),
              ((1.5, 6, 8), None))


def diameter_cases_csv() -> bytes:
    import csv as _csv
    buf = io.StringIO()
    w = _csv.writer(buf, lineterminator="\n")
    w.writerow(["text", "value", "unit", "flag"])
    w.writerows(DIAMETER_CASES)
    return buf.getvalue().encode("utf-8")


def no_prj_zip() -> bytes:
    return _shp_zip("mains", [[ln] for ln in _quirk_lines()], [[f'{d}"'] for d in QUIRK_DIAMETERS],
                    [("PIPEDIAM", "C", 8)], prj=None)


def points_only_geojson() -> bytes:
    feats = [{"type": "Feature", "properties": {"HYDID": f"H{k}"},
              "geometry": {"type": "Point", "coordinates": [-100.6 + 0.001 * k, 40.0]}} for k in range(5)]
    return json.dumps({"type": "FeatureCollection", "features": feats}).encode()


def empty_layer_gpkg() -> bytes:
    return _gpkg_file([("mains", "LINESTRING", [])], [("PIPEDIAM", "text")])


VALVE_STATUS = (("V1", "OPEN"), ("V2", "CLOSED"), ("V3", "Normally Closed"), ("V4", "C"), ("V5", "O"), ("V6", ""),
                ("V7", "closed"))


def valves_normal_pos_geojson() -> bytes:
    """A small mains layer (four lines meeting in a row) and its valves as points in one file: NORMAL_POS gives three
    closed, two open, one blank, and one closed valve away from any node."""
    feats = []
    xs = [6000.0 + 300 * k for k in range(5)]
    for k in range(4):
        feats.append({"type": "Feature", "properties": {"PIPEDIAM": "8 in"},
                      "geometry": {"type": "LineString", "coordinates": [[xs[k], 6000.0], [xs[k + 1], 6000.0]]}})
    spots = [(xs[0], 6000.0), (xs[1] + 1.0, 6000.5), (xs[2], 6000.0), (xs[3], 6000.0), (xs[4], 6000.0),
             (xs[1], 6000.0), (6450.0, 6900.0)]
    for (vid, st), (x, y) in zip(VALVE_STATUS, spots):
        feats.append({"type": "Feature", "properties": {"VALVEID": vid, "NORMAL_POS": st},
                      "geometry": {"type": "Point", "coordinates": [x, y]}})
    return json.dumps({"type": "FeatureCollection", "crs": {"type": "name",
                                                            "properties": {"name": "urn:ogc:def:crs:EPSG::2227"}},
                       "features": feats}).encode()


def mains_and_valves_gpkg() -> bytes:
    """valves_normal_pos.geojson's mains and valves as two layers of one GeoPackage (in its feet), as a utility keeps
    them: a line layer 'mains' and a point layer 'valves' with VALVEID and NORMAL_POS."""
    d = json.loads(valves_normal_pos_geojson())
    lines = [f for f in d["features"] if f["geometry"]["type"] == "LineString"]
    pts = [f for f in d["features"] if f["geometry"]["type"] == "Point"]
    mains = [(_gpkg_blob(_wkb_line([f["geometry"]["coordinates"]]), GPKG_SRS_ID), [f["properties"]["PIPEDIAM"]])
             for f in lines]
    valves = [(_gpkg_blob(struct.pack("<BIdd", 1, 1, *f["geometry"]["coordinates"]), GPKG_SRS_ID),
               [f["properties"]["VALVEID"], f["properties"]["NORMAL_POS"]]) for f in pts]
    return _gpkg_file([("mains", "LINESTRING", mains, [("PIPEDIAM", "text")]),
                       ("valves", "POINT", valves, [("VALVEID", "text"), ("NORMAL_POS", "text")])], [])


# the review of task 17: one fixture for each issue it found, built by the writers above (never by the readers)
def _grid_lines(n: int, step: float, x0: float, y0: float) -> list:
    """An n by n street grid of two-point lines, step apart, from (x0, y0)."""
    out = []
    for i in range(n):
        for j in range(n):
            x, y = x0 + i * step, y0 + j * step
            if i < n - 1:
                out.append([(x, y), (x + step, y)])
            if j < n - 1:
                out.append([(x, y), (x, y + step)])
    return out


def _fc(lines, props, crs: str | None = None) -> bytes:
    d = {"type": "FeatureCollection",
         "features": [{"type": "Feature", "properties": p, "geometry": {"type": "LineString",
                                                                       "coordinates": [list(q) for q in ln]}}
                      for ln, p in zip(lines, props)]}
    if crs:
        d["crs"] = {"type": "name", "properties": {"name": crs}}
    return json.dumps(d).encode()


def _shp_members(stem: str, records, rows, fields, prj: str | None = PRJ_FT, cpg: str | None = "UTF-8") -> list:
    """_shp_zip's members, for a zip that holds more than one layer or other files."""
    shp, shx = _shp_lines(records)
    out = [(f"{stem}.shp", shp), (f"{stem}.shx", shx), (f"{stem}.dbf", _dbf_text(fields, rows, "utf-8"))]
    if prj is not None:
        out.append((f"{stem}.prj", prj.encode("ascii")))
    if cpg is not None:
        out.append((f"{stem}.cpg", cpg.encode("ascii")))
    return out


GEOG_WGS84 = ('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",'
              '0.0],UNIT["Degree",0.0174532925199433]]')
REVIEW_GRID = (4, 400.0, 1.0e6, 2.0e5)            # 24 lines of 400 ft, in PRJ_FT's US feet
REVIEW_BAD_VERTICES = (float("nan"), float("inf"), 1e300)
REVIEW_LONLAT = (30, 0.001, -84.5, 38.0)          # 1,740 lines, blocks of 0.001 degree (about 90 by 110 m)
REVIEW_MAX_FEATURES = 50                          # the feature cap the over-limit fixtures are read with
REVIEW_POINTS = 60                                # ... and their points, over that cap


def bad_vertices_gpkg() -> bytes:
    """REVIEW_GRID's 24 lines and three more, each with one vertex of NaN, infinity or 1e300 (the review: one NaN vertex
    gave the Get started page a traceback)."""
    x0, y0 = REVIEW_GRID[2:]
    lines = _grid_lines(*REVIEW_GRID) + [[(x0, y0), (v, y0 + 400.0)] for v in REVIEW_BAD_VERTICES]
    rows = [(_gpkg_blob(_wkb_line([ln]), GPKG_SRS_ID), ['8"']) for ln in lines]
    return _gpkg_file([("mains", "LINESTRING", rows)], [("PIPEDIAM", "text")])


def empty_lines_4326_gpkg() -> bytes:
    """A GeoPackage in longitude and latitude (EPSG 4326) whose three lines are all LINESTRING EMPTY."""
    rows = [(_gpkg_blob(struct.pack("<BII", 1, 2, 0), 4326), ['8"']) for _ in range(3)]
    return _gpkg_file([("mains", "LINESTRING", rows)], [("PIPEDIAM", "text")], srs=(4326, "EPSG", 4326, GEOG_WGS84))


def lonlat_no_prj_zip() -> bytes:
    """A street grid in longitude and latitude, zipped with no .prj (the review: read as feet, every end fell in one
    node, and joining them took time and memory quadratic in the lines)."""
    lines = _grid_lines(*REVIEW_LONLAT)
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'] for _ in lines], [("PIPEDIAM", "C", 8)], prj=None)


def lonlat_as_3857_geojson() -> bytes:
    """Longitude and latitude in a file whose crs member says EPSG:3857."""
    lines = _grid_lines(8, 0.001, -84.5, 38.0)
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::3857")


def stateplane_as_4326_geojson() -> bytes:
    """State Plane feet in a file whose crs member says EPSG:4326 (a Define Projection mix-up)."""
    lines = _grid_lines(8, 400.0, 1.5e6, 2.0e5)
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::4326")


# Web Mercator (EPSG 3857): an east-west main of MERC_LENGTHS_M[0] and a north-south one of MERC_LENGTHS_M[1] along the
# WGS84 ellipsoid at MERC_LAT, projected by Web Mercator's forward formula (its sphere, the ellipsoid's latitude); in
# Web Mercator metres they measure about 1 / cos(45 degrees) = 1.41 times longer
MERC_LON, MERC_LAT = -93.0, 45.0
MERC_LENGTHS_M = (1000.0, 500.0)
WEB_MERC_PRJ = ('PROJCS["WGS_1984_Web_Mercator_Auxiliary_Sphere",' + GEOG_WGS84 + ',PROJECTION["Mercator_Auxiliary_'
                'Sphere"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",'
                '0.0],PARAMETER["Standard_Parallel_1",0.0],PARAMETER["Auxiliary_Sphere_Type",0.0],UNIT["Meter",1.0]]')
WORLD_MERC_PRJ = ('PROJCS["World_Mercator",' + GEOG_WGS84 + ',PROJECTION["Mercator"],PARAMETER["False_Easting",0.0],'
                  'PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",0.0],PARAMETER["Standard_Parallel_1",'
                  '0.0],UNIT["Meter",1.0]]')


def _mercator_lines() -> list:
    a, e2 = 6378137.0, 6.69437999014e-3
    phi = math.radians(MERC_LAT)
    w = 1.0 - e2 * math.sin(phi) ** 2
    dlon = math.degrees(MERC_LENGTHS_M[0] / (a / math.sqrt(w) * math.cos(phi)))
    dlat = math.degrees(MERC_LENGTHS_M[1] / (a * (1.0 - e2) / w ** 1.5))
    fwd = lambda lon, lat: (6378137.0 * math.radians(lon),
                            6378137.0 * math.log(math.tan(math.pi / 4.0 + math.radians(lat) / 2.0)))
    p0, p1, p2 = (MERC_LON, MERC_LAT), (MERC_LON + dlon, MERC_LAT), (MERC_LON + dlon, MERC_LAT + dlat)
    return [[fwd(*p0), fwd(*p1)], [fwd(*p1), fwd(*p2)]]


def mercator_3857_geojson() -> bytes:
    return _fc(_mercator_lines(), [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "8 in"}], "urn:ogc:def:crs:EPSG::3857")


def mercator_esri_prj_zip() -> bytes:
    lines = _mercator_lines()
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'], ['8"']], [("PIPEDIAM", "C", 8)], prj=WEB_MERC_PRJ)


def world_mercator_zip() -> bytes:
    lines = _mercator_lines()
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'], ['8"']], [("PIPEDIAM", "C", 8)], prj=WORLD_MERC_PRJ)


def mains_and_many_points_gpkg() -> bytes:
    """REVIEW_GRID's mains and a layer of REVIEW_POINTS meters (points), read with the cap at REVIEW_MAX_FEATURES (the
    review: a large meters layer beside the mains refused the whole file)."""
    x0, y0 = REVIEW_GRID[2:]
    mains = [(_gpkg_blob(_wkb_line([ln]), GPKG_SRS_ID), ['8"']) for ln in _grid_lines(*REVIEW_GRID)]
    meters = [(_gpkg_blob(struct.pack("<BIdd", 1, 1, x0 + 7.0 * k, y0 + 13.0), GPKG_SRS_ID), [f"M{k:04d}"])
              for k in range(REVIEW_POINTS)]
    return _gpkg_file([("mains", "LINESTRING", mains, [("PIPEDIAM", "text")]),
                       ("meters", "POINT", meters, [("METER_ID", "text")])], [])


def _esri(lines, sr) -> dict:
    return {"geometryType": "esriGeometryPolyline", "spatialReference": sr,
            "features": [{"attributes": {"DIAMETER": "8"}, "geometry": {"paths": [[list(p) for p in ln]]}}
                         for ln in lines]}


def esri_over_limit_json() -> bytes:
    """Esri JSON of REVIEW_POINTS lines, read with the cap at REVIEW_MAX_FEATURES (the review: cut at the cap without a
    word, where GeoJSON was refused)."""
    lines = _grid_lines(7, 400.0, 1.0e6, 2.0e5)          # 84 lines
    return json.dumps(_esri(lines[:REVIEW_POINTS], {"wkid": 2263})).encode()


def esri_code_only_json() -> bytes:
    """Esri JSON whose spatialReference is the code alone (the review: refused as an AttributeError)."""
    return json.dumps(_esri(_grid_lines(*REVIEW_GRID), 2263)).encode()


def geojson_bad_positions_geojson() -> bytes:
    """REVIEW_GRID's lines, two of them with a position that cannot be read (a null position; 2D and 3D mixed)."""
    lines = [[list(p) for p in ln] for ln in _grid_lines(*REVIEW_GRID)]
    lines[3] = [lines[3][0], [None, None]]
    lines[5] = [lines[5][0] + [10.0], lines[5][1]]
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::2263")


def shp_with_preview_kml_zip() -> bytes:
    """A zipped shapefile that also carries a preview.kml (the review: the zip was read as a KMZ)."""
    lines = _grid_lines(*REVIEW_GRID)
    kml = ('<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><Placemark>'
           '<name>preview</name><LineString><coordinates>-98.0,36.0,0 -97.99,36.0,0</coordinates></LineString>'
           '</Placemark></Document></kml>')
    return _zip(_shp_members("mains", [[ln] for ln in lines], [['8"'] for _ in lines], [("PIPEDIAM", "C", 8)])
                + [("preview.kml", kml.encode())])


def two_folders_same_name_zip() -> bytes:
    """north/mains and south/mains, two shapefiles of one name (the review: both were layer 'mains', and the second
    could not be picked)."""
    north, south = _grid_lines(*REVIEW_GRID), _grid_lines(2, 400.0, 1.2e6, 2.0e5)
    return _zip(_shp_members("north/mains", [[ln] for ln in north], [['8"'] for _ in north], [("PIPEDIAM", "C", 8)])
                + _shp_members("south/mains", [[ln] for ln in south], [['6"'] for _ in south], [("PIPEDIAM", "C", 8)]))


MANY_SHAPEFILES = 101                             # one more than the reader takes from one zip


def many_shapefiles_zip() -> bytes:
    """MANY_SHAPEFILES shapefiles of one line each (the review: every layer's members listed the whole zip again)."""
    members = []
    for j in range(MANY_SHAPEFILES):
        members += _shp_members(f"l{j:04d}", [[[(1.0e6 + 500.0 * j, 2.0e5), (1.0e6 + 500.0 * j + 100.0, 2.0e5)]]],
                                [['8"']], [("PIPEDIAM", "C", 8)])
    return _zip(members)


def _kml_doc(marks: list[str]) -> bytes:
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document>'
            + "".join(marks) + "</Document></kml>\n").encode()


def _kml_path(i: int, name: str, desc: str = "") -> str:
    x0 = -84.5 + 0.001 * i
    return (f"<Placemark><name>{name}</name>" + (f"<description>{desc}</description>" if desc else "")
            + f"<LineString><coordinates>{x0:.6f},38.000000,0 {x0 + 0.001:.6f},38.000000,0</coordinates></LineString>"
              f"</Placemark>")


def kml_numbered_paths_kml() -> bytes:
    """20 paths named by their author 'Pipe 2' to 'Pipe 21' (the review: read as sizes of 2 to 21 in)."""
    return _kml_doc([_kml_path(i, f"Pipe {i + 2}") for i in range(20)])


KML_TABLE_DIAMETER = 8
KML_TABLE_LENGTH_FT = 287.8               # each path's length on the ground, 0.001 degree of longitude at 38 N


def kml_attribute_table_kml() -> bytes:
    """20 placemarks named by facility ID, each with the attribute table ArcGIS's Layer To KML writes in its
    description (the review: the table parsed as two sizes, and the ID became the size)."""
    tbl = lambda i: (f"<![CDATA[<html><body><table><tr><th>FACILITYID</th><td>W-{i}</td></tr><tr><td>DIAMETER</td>"
                     f"<td>{KML_TABLE_DIAMETER}</td></tr><tr><td>INSTALLYEAR</td><td>1965</td></tr><tr><td>"
                     f"Shape_Length</td><td>{KML_TABLE_LENGTH_FT}</td></tr></table></body></html>]]>")
    return _kml_doc([_kml_path(i, f"W-{i}", tbl(i)) for i in range(20)])


def year_all_null_geojson() -> bytes:
    """REVIEW_GRID's lines with an install-date field that is empty everywhere (the review: 'read for 0 pipes (None to
    None)')."""
    lines = _grid_lines(*REVIEW_GRID)
    return _fc(lines, [{"PIPEDIAM": "8 in", "INSTALLDATE": None} for _ in lines], "urn:ogc:def:crs:EPSG::2263")


# tees and loops (the review of task 17), each at the 6 ft default, in US feet (EPSG 2227)
TEE_CASES = ("tee_dense_vertex_main.geojson", "tee_offset_cross.geojson", "loop_main_closed.geojson")


def tee_dense_vertex_main_geojson() -> bytes:
    """A main with a vertex every 10 ft and a branch ending 4 ft short of it, its foot 5 ft from the nearest vertex (the
    review: the main was cut at that vertex, 6.4 ft from the branch end, which then joined nothing)."""
    main = [(float(x), 0.0) for x in range(0, 101, 10)]
    return _fc([main, [(55.0, -204.0), (55.0, -4.0)]], [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "6 in"}],
               "urn:ogc:def:crs:EPSG::2227")


def tee_offset_cross_geojson() -> bytes:
    """A straight main and two branches from either side, ending 4 ft short of it with their feet 5 ft apart (the
    review: one cut, at the second foot, 6.4 ft from the first branch end)."""
    return _fc([[(0.0, 0.0), (100.0, 0.0)], [(50.0, 204.0), (50.0, 4.0)], [(55.0, -204.0), (55.0, -4.0)]],
               [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "6 in"}, {"PIPEDIAM": "6 in"}], "urn:ogc:def:crs:EPSG::2227")


def loop_main_closed_geojson() -> bytes:
    """A feed main and a loop main drawn as one closed line of 24 vertices, starting and ending where the feed ends
    (the review: the loop was dropped, as a pipe from a node to itself)."""
    ring = [(100.0 * math.sin(2 * math.pi * k / 24), 100.0 - 100.0 * math.cos(2 * math.pi * k / 24)) for k in range(24)]
    ring.append(ring[0])
    return _fc([[(0.0, -300.0), (0.0, 0.0)], ring], [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "6 in"}],
               "urn:ogc:def:crs:EPSG::2227")


REVIEW_FIXTURES = ("bad_vertices.gpkg", "empty_lines_4326.gpkg", "lonlat_no_prj.zip", "lonlat_as_3857.geojson",
                   "stateplane_as_4326.geojson", "mercator_3857.geojson", "mercator_esri_prj.zip", "world_mercator.zip",
                   "mains_and_many_points.gpkg", "esri_over_limit.json", "esri_code_only.json",
                   "geojson_bad_positions.geojson", "shp_with_preview_kml.zip", "two_folders_same_name.zip",
                   "many_shapefiles.zip", "kml_numbered_paths.kml", "kml_attribute_table.kml",
                   "year_all_null.geojson") + TEE_CASES


# the final verification of task 17: each issue it confirmed, as a fixture of the check gis_final_verification
FINAL_BIG_CHARS = 4_900_000                # the long run in the KMZ fixtures: within the hosted demo's 5 MB unpacked
FINAL_PRJ_RUN = (3_000, 10_000)            # the runs of spaces around a comma in the .prj fixtures (738 B and 20 KB)
FINAL_EPOCH_YEARS = (1931, 1948, 1965, 1972, 1990, 2004)


def _final_kmz(extra_mark: str) -> bytes:
    """Seven paths of 8 in end to end along 38 N, and one more (extra_mark), as doc.kml in a KMZ."""
    return _zip([("doc.kml", _kml_doc([_kml_path(i, '8"') for i in range(7)] + [extra_mark]))])


def kmz_long_whitespace_kmz() -> bytes:
    """A path whose coordinates hold a run of FINAL_BIG_CHARS spaces between its two points (the final verification:
    a pattern took time quadratic in the run, about 3.5 h for this file in the hosted demo)."""
    x0 = -84.5 + 0.007
    return _final_kmz(f'<Placemark><name>8"</name><LineString><coordinates>{x0:.6f},38.000000,0' + " " * FINAL_BIG_CHARS
                      + f'{x0 + 0.001:.6f},38.000000,0</coordinates></LineString></Placemark>')


def kmz_long_lt_kmz() -> bytes:
    """A path whose description is FINAL_BIG_CHARS / 4 escaped '<' (the final verification: stripping its tags took
    time quadratic in them)."""
    return _final_kmz(_kml_path(7, '8"', "&lt;" * (FINAL_BIG_CHARS // 4)))


def kmz_open_rows_kmz() -> bytes:
    """A path whose description is table rows never closed, about FINAL_BIG_CHARS long (the final verification: the
    attribute-table pattern took time quadratic in them)."""
    return _final_kmz(_kml_path(7, '8"', "<![CDATA[" + "<tr><td>a</td>" * (FINAL_BIG_CHARS // 14) + "]]>"))


def _doctype_kml(lead: str, where: str) -> bytes:
    """A KML whose name uses an internal entity, the declaration after lead (blank lines before the XML declaration)
    or where (a comment after the XML declaration)."""
    body = ('<kml xmlns="http://www.opengis.net/kml/2.2"><Document><Placemark><name>&e; 8"</name><LineString>'
            '<coordinates>-84.5,38.0,0 -84.499,38.0,0</coordinates></LineString></Placemark></Document></kml>')
    return _zip([("doc.kml", (lead + '<?xml version="1.0" encoding="UTF-8"?>\n' + where
                              + '<!DOCTYPE kml [<!ENTITY e "EXPANDED">]>' + body).encode())])


def kmz_doctype_after_blank_lines_kmz() -> bytes:
    return _doctype_kml("\n" * 200_001, "")


def kmz_doctype_after_comment_kmz() -> bytes:
    return _doctype_kml("", "<!--" + "x" * 200_000 + "-->")


def _prj_run(n: int) -> str:
    return 'PROJCS["x",' + " " * n + "," + " " * n


def prj_spaces_comma_zip() -> bytes:
    """REVIEW_GRID's mains with a .prj of 'PROJCS["x",', FINAL_PRJ_RUN[0] spaces, a comma and as many spaces, no
    closing bracket (the final verification: 28 s for this zip of under 1 KB, cubic in the run)."""
    lines = _grid_lines(*REVIEW_GRID)
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'] for _ in lines], [("PIPEDIAM", "C", 8)],
                    prj=_prj_run(FINAL_PRJ_RUN[0]), cpg=None)


def prj_spaces_comma_20kb_zip() -> bytes:
    lines = _grid_lines(*REVIEW_GRID)
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'] for _ in lines], [("PIPEDIAM", "C", 8)],
                    prj=_prj_run(FINAL_PRJ_RUN[1]), cpg=None)


def esri_wkt_spaces_comma_json() -> bytes:
    """Esri JSON whose spatialReference wkt is the same run (wkt_unit reads it as the .prj is read)."""
    return json.dumps(_esri(_grid_lines(*REVIEW_GRID), {"wkt": _prj_run(FINAL_PRJ_RUN[1])})).encode()


def mains_and_service_lines_gpkg() -> bytes:
    """REVIEW_GRID's mains beside a layer of REVIEW_POINTS service lines and one of REVIEW_POINTS hydrants declared
    GEOMETRY, both over the cap of REVIEW_MAX_FEATURES they are read with (the final verification: a large line layer,
    or points declared GEOMETRY, refused the whole GeoPackage, even with the mains layer chosen)."""
    x0, y0 = REVIEW_GRID[2:]
    mains = [(_gpkg_blob(_wkb_line([ln]), GPKG_SRS_ID), ['8"']) for ln in _grid_lines(*REVIEW_GRID)]
    service = [(_gpkg_blob(_wkb_line([[(x0 + 20.0 * k, y0 + 5.0), (x0 + 20.0 * k, y0 + 45.0)]]), GPKG_SRS_ID), ['1"'])
               for k in range(REVIEW_POINTS)]
    hydrants = [(_gpkg_blob(struct.pack("<BIdd", 1, 1, x0 + 9.0 * k, y0 + 3.0), GPKG_SRS_ID), [f"H{k:03d}"])
                for k in range(REVIEW_POINTS)]
    return _gpkg_file([("mains", "LINESTRING", mains, [("PIPEDIAM", "text")]),
                       ("service_lines", "MULTILINESTRING", service, [("PIPEDIAM", "text")]),
                       ("hydrants", "GEOMETRY", hydrants, [("HYD_ID", "text")])], [])


def mains_and_service_lines_zip() -> bytes:
    """The same mains and service lines as two shapefiles in one zip."""
    x0, y0 = REVIEW_GRID[2:]
    mains = _grid_lines(*REVIEW_GRID)
    service = [[(x0 + 20.0 * k, y0 + 5.0), (x0 + 20.0 * k, y0 + 45.0)] for k in range(REVIEW_POINTS)]
    return _zip(_shp_members("mains", [[ln] for ln in mains], [['8"'] for _ in mains], [("PIPEDIAM", "C", 8)])
                + _shp_members("service_lines", [[ln] for ln in service], [['1"'] for _ in service],
                               [("PIPEDIAM", "C", 8)]))


def world_mercator_3395_geojson() -> bytes:
    """Lines in World Mercator, named by its code alone (the final verification: its unit was asked, and metres gave
    lengths 1.41 times too long at 45 N)."""
    return _fc(_mercator_lines(), [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "8 in"}], "urn:ogc:def:crs:EPSG::3395")


def esri_wkid_3395_json() -> bytes:
    return json.dumps(_esri(_mercator_lines(), {"wkid": 3395})).encode()


FINAL_SWAP_GRID = (4, 0.002, -77.03, 38.89)      # a small grid at 38.89 N, 77.03 W, in longitude and latitude


def _swapped_grid() -> list:
    return [[(lat, lon) for lon, lat in ln] for ln in _grid_lines(*FINAL_SWAP_GRID)]


def swapped_latlon_4326_geojson() -> bytes:
    """FINAL_SWAP_GRID written latitude first, with a crs of EPSG::4326 (the final verification: it built near the
    pole, with lengths of 50 to 223 m against 173 to 222 m, and no question)."""
    lines = _swapped_grid()
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::4326")


def swapped_latlon_nocrs_geojson() -> bytes:
    lines = _swapped_grid()
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines])


def lonlat_stateplane_prj_zip() -> bytes:
    """A 6 by 6 grid in longitude and latitude, zipped with a State Plane .prj in US feet (the final verification: it
    was refused as crowded and told to check the unit, with no control to change it)."""
    lines = _grid_lines(6, 0.001, -84.5, 38.0)
    return _shp_zip("mains", [[ln] for ln in lines], [['8"'] for _ in lines], [("PIPEDIAM", "C", 8)], prj=PRJ_FT)


def lonlat_as_2263_geojson() -> bytes:
    lines = _grid_lines(6, 0.001, -84.5, 38.0)
    return _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::2263")


def _branches(n: int, s: float, x0: float = 0.0) -> list:
    """n branches ending 4 ft short of a main along y = 0, s ft apart from x0, fanned out 60 ft apart at 150 ft."""
    return [[(x0 + s * k + 60.0 * (k - n / 2), 150.0), (x0 + s * k, 4.0)] for k in range(n)]


def tee_fan_dense_main_geojson() -> bytes:
    """Four branches 5 ft apart ending 4 ft short of a main drawn with a vertex every 2 ft (the final verification: a
    split of the wide node put two branch ends in a node of their own, apart from the main)."""
    main = [(float(-300 + 2 * k), 0.0) for k in range(301)]
    lines = [main] + _branches(4, 5.0)
    return _fc(lines, [{"PIPEDIAM": "8 in"}] + [{"PIPEDIAM": "6 in"}] * 4, "urn:ogc:def:crs:EPSG::2227")


def tee_fan_straight_main_geojson() -> bytes:
    """Five branches 3 ft apart ending 4 ft short of a straight two-point main (the same, two connected parts)."""
    lines = [[(-300.0, 0.0), (300.0, 0.0)]] + _branches(5, 3.0)
    return _fc(lines, [{"PIPEDIAM": "8 in"}] + [{"PIPEDIAM": "6 in"}] * 5, "urn:ogc:def:crs:EPSG::2227")


def loop_joined_through_feed_geojson() -> bytes:
    """A loop main of 41 vertices whose two ends are 8 ft apart, each 4 ft from the end of a feed main, so the two
    ends join into one node only through the feed's end (the final verification: the loop, 321.8 ft, was dropped)."""
    ring = [(50.0 * math.sin(t), 50.0 - 50.0 * math.cos(t)) for t in (2 * math.pi * k / 40 for k in range(41))]
    ring[0], ring[-1] = (-4.0, 0.0), (4.0, 0.0)
    return _fc([[(0.0, -300.0), (0.0, 0.0)], ring], [{"PIPEDIAM": "8 in"}, {"PIPEDIAM": "6 in"}],
               "urn:ogc:def:crs:EPSG::2227")


def kml_line_in_street_kml() -> bytes:
    """20 paths named 'Line 1 in Oak St' to 'Line 20 in Oak St' (the final verification: 'in' after a number was read
    as inches, 1 to 20 in, with no question)."""
    return _kml_doc([_kml_path(i, f"Line {i + 1} in Oak St") for i in range(20)])


def esri_epoch_dates_json() -> bytes:
    """12 pipes whose INSTALLDATE is epoch milliseconds, as Esri JSON writes a date, for FINAL_EPOCH_YEARS (the final
    verification: digits inside them were read as years, 1984 to 2000)."""
    import datetime as _dt
    lines = _grid_lines(*REVIEW_GRID)[:12]
    ms = [int((_dt.datetime(FINAL_EPOCH_YEARS[k % 6], 6, 1) - _dt.datetime(1970, 1, 1)).total_seconds() * 1000)
          for k in range(12)]
    d = _esri(lines, {"wkid": 2263})
    for f, v in zip(d["features"], ms):
        f["attributes"]["INSTALLDATE"] = v
    return json.dumps(d).encode()


def _grid_gpkg() -> bytes:
    rows = [(_gpkg_blob(_wkb_line([ln]), GPKG_SRS_ID), ['8"']) for ln in _grid_lines(*REVIEW_GRID)]
    return _gpkg_file([("mains", "LINESTRING", rows)], [("PIPEDIAM", "text")])


def zipped_gpkg_zip() -> bytes:
    """A zip holding one GeoPackage of REVIEW_GRID's mains (the final verification: refused with a sentence that named
    GeoPackage among the layers read)."""
    return _zip([("mains.gpkg", _grid_gpkg())])


def zipped_geojson_zip() -> bytes:
    lines = _grid_lines(*REVIEW_GRID)
    return _zip([("mains.geojson", _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::2263"))])


def zipped_two_geojson_zip() -> bytes:
    lines = _grid_lines(*REVIEW_GRID)
    fc = _fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines], "urn:ogc:def:crs:EPSG::2263")
    return _zip([("north.geojson", fc), ("south.geojson", fc)])


def geojson_crs_props_string_geojson() -> bytes:
    """A crs member whose properties are a string (the fuzzer: AttributeError on the page)."""
    lines = _grid_lines(*REVIEW_GRID)
    d = json.loads(_fc(lines, [{"PIPEDIAM": "8 in"} for _ in lines]))
    d["crs"] = {"type": "name", "properties": "EPSG:2263"}
    return json.dumps(d).encode()


def esri_odd_attributes_json() -> bytes:
    """Esri JSON with one feature's attributes a list, one a string, and one feature's paths a number (the fuzzer)."""
    d = _esri(_grid_lines(*REVIEW_GRID), {"wkid": 2263})
    d["features"][0]["attributes"] = ["DIAMETER", "8"]
    d["features"][1]["attributes"] = "8"
    d["features"][2]["geometry"]["paths"] = 7
    return json.dumps(d).encode()


def gpkg_bad_utf8_gpkg() -> bytes:
    """REVIEW_GRID's mains with one size stored as text that is not valid UTF-8 (the fuzzer: UnicodeDecodeError)."""
    con = sqlite3.connect(":memory:")
    try:
        con.deserialize(_grid_gpkg())
        con.execute("update mains set PIPEDIAM = CAST(x'38228a' AS TEXT) where fid = 1")
        con.commit()
        return con.serialize()
    finally:
        con.close()


def size_overflow_geojson() -> bytes:
    """REVIEW_GRID's mains in bare inches, one of them '9' written 400 times (the fuzzer: it parsed as inf, and inf was
    written as a diameter)."""
    lines = _grid_lines(*REVIEW_GRID)
    return _fc(lines, [{"PIPEDIAM": "9" * 400 if k == 0 else "8"} for k in range(len(lines))],
               "urn:ogc:def:crs:EPSG::2263")


FINAL_FIXTURES = ("kmz_long_whitespace.kmz", "kmz_long_lt.kmz", "kmz_open_rows.kmz", "kmz_doctype_after_blank_lines.kmz",
                  "kmz_doctype_after_comment.kmz", "prj_spaces_comma.zip", "prj_spaces_comma_20kb.zip",
                  "esri_wkt_spaces_comma.json", "mains_and_service_lines.gpkg", "mains_and_service_lines.zip",
                  "world_mercator_3395.geojson", "esri_wkid_3395.json", "swapped_latlon_4326.geojson",
                  "swapped_latlon_nocrs.geojson", "lonlat_stateplane_prj.zip", "lonlat_as_2263.geojson",
                  "tee_fan_dense_main.geojson", "tee_fan_straight_main.geojson", "loop_joined_through_feed.geojson",
                  "kml_line_in_street.kml", "esri_epoch_dates.json", "zipped_gpkg.zip", "zipped_geojson.zip",
                  "zipped_two_geojson.zip", "geojson_crs_props_string.geojson", "esri_odd_attributes.json",
                  "gpkg_bad_utf8.gpkg", "size_overflow.geojson")


# WNTR's own GIS schema (critique 24), written here the way WNTR's to_gis lays it out, as far as this writer emulates
# it (each element's base and optional attributes that are plain values, where to_gis keeps every plain column of
# to_dict; coordinates and vertices as the geometry) and WNTR's JSON as the sidecar: never a file WNTR itself wrote
# (to_gis needs geopandas, which is not installed; the review of task 17)
WNTR_SCHEMA_NETS = ("Net3", "Net2", "ky4")


def wntr_schema_zip(net: str) -> bytes:
    import wntr
    from wntr.network.elements import Junction, Pipe, Pump, Reservoir, Tank, Valve
    wn = _wn(net)
    d = wntr.network.to_dict(wn)
    keep = {"Junction": Junction, "Tank": Tank, "Reservoir": Reservoir, "Pipe": Pipe, "Pump": Pump, "Valve": Valve}
    layer_of = {"Junction": "junctions", "Tank": "tanks", "Reservoir": "reservoirs", "Pipe": "pipes", "Pump": "pumps",
                "Valve": "valves"}
    out = {v: [] for v in layer_of.values()}
    xy = {n["name"]: n["coordinates"] for n in d["nodes"]}
    for el in d["nodes"] + d["links"]:
        t = el.get("node_type") or el.get("link_type")
        cls = keep[t]
        names = [a for a in cls._base_attributes + cls._optional_attributes if a not in ("coordinates", "vertices")]
        props = {a: el.get(a) for a in names if isinstance(el.get(a), (str, int, float, bool, type(None)))}
        if "node_type" in el:
            geom = {"type": "Point", "coordinates": list(el["coordinates"])}
        else:
            pts = [xy[el["start_node_name"]]] + list(el.get("vertices") or []) + [xy[el["end_node_name"]]]
            geom = {"type": "LineString", "coordinates": [list(p) for p in pts]}
        out[layer_of[t]].append({"type": "Feature", "properties": props, "geometry": geom})
    members = [(f"{net}_{k}.geojson", json.dumps({"type": "FeatureCollection", "features": v}).encode())
               for k, v in out.items() if v]
    members.append((f"{net}_wn.json", json.dumps(d, default=str).encode()))
    return _zip(members)


GIS_FIXTURES = ("ky4_mains_stdlib_shp.zip", "ky4_mains_pyshp.zip", "ky4_mains_z.gpkg", "ky4_mains_lonlat.geojson",
                "ky4_googleearth.kmz", "ky4_mains_shuffled.geojson", "multipart_two_parts.gpkg",
                "multipart_two_parts.zip", "null_shape_mid.zip", "deleted_dbf_record.zip", "cpg_utf8.zip",
                "wkb_z_m_zm.gpkg", "prj_stateplane_ftus.zip", "gpkg_srs_undefined.gpkg", "shape_length_mismatch.zip",
                "quantized_4dp.geojson", "chain_cluster.geojson", "ids_dirty.geojson", "diameter_cases.csv",
                "no_prj.zip", "utility_gdb.zip", "points_only.geojson", "empty_layer.gpkg",
                "valves_normal_pos.geojson", "mains_and_valves.gpkg") + tuple(f"wntr_schema_{n}.zip"
                                                                             for n in WNTR_SCHEMA_NETS) + REVIEW_FIXTURES \
    + FINAL_FIXTURES
GIS_SCORED = ("ky4_mains_stdlib_shp.zip", "ky4_mains_pyshp.zip", "ky4_mains_z.gpkg", "ky4_mains_lonlat.geojson")
FIXTURES.update({"ky4_mains_stdlib_shp.zip": ky4_mains_stdlib_shp_zip, "ky4_mains_pyshp.zip": ky4_mains_pyshp_zip,
                 "ky4_mains_z.gpkg": ky4_mains_z_gpkg, "ky4_mains_lonlat.geojson": ky4_mains_lonlat_geojson,
                 "ky4_googleearth.kmz": ky4_googleearth_kmz, "ky4_mains_shuffled.geojson": ky4_mains_shuffled_geojson,
                 "multipart_two_parts.gpkg": multipart_two_parts_gpkg, "multipart_two_parts.zip": multipart_two_parts_zip,
                 "null_shape_mid.zip": null_shape_mid_zip, "deleted_dbf_record.zip": deleted_dbf_record_zip,
                 "cpg_utf8.zip": cpg_utf8_zip, "wkb_z_m_zm.gpkg": wkb_z_m_zm_gpkg,
                 "prj_stateplane_ftus.zip": prj_stateplane_ftus_zip, "gpkg_srs_undefined.gpkg": gpkg_srs_undefined_gpkg,
                 "shape_length_mismatch.zip": shape_length_mismatch_zip, "quantized_4dp.geojson": quantized_4dp_geojson,
                 "chain_cluster.geojson": chain_cluster_geojson, "ids_dirty.geojson": ids_dirty_geojson,
                 "diameter_cases.csv": diameter_cases_csv, "no_prj.zip": no_prj_zip,
                 "points_only.geojson": points_only_geojson, "empty_layer.gpkg": empty_layer_gpkg,
                 "valves_normal_pos.geojson": valves_normal_pos_geojson, "mains_and_valves.gpkg": mains_and_valves_gpkg,
                 **{f"wntr_schema_{n}.zip": functools.partial(wntr_schema_zip, n) for n in WNTR_SCHEMA_NETS},
                 **{n: globals()[n.replace(".", "_")] for n in REVIEW_FIXTURES + FINAL_FIXTURES}})


# ----------------------------------------------------------------------------- task 18: CAD drawings
# Every drawing here is written by hand, group code by group code from the DXF reference, as ASCII (R12 and R2000) or
# binary DXF (R12), sharing no code with ezdxf, which reads them; the one exception is ky4_r2018_lwpolyline.dxf, which
# ezdxf writes and ezdxf reads (same-library self-consistency, labelled so).  The ky4 drawing is task 17's mains
# layer (ky4_mains_records: 43 runs drawn through tees, 19 undershoots of 4 ft, every end jittered by up to 0.3 ft),
# drawn as POLYLINEs, half the sizes in layer names and half as labels (pre-registered before any reader code: the
# label 2 ft high, at the middle of its line's longest segment, parallel to it, on a seeded side, the near edge of its
# box 0.5 to 2.5 heights off the line, so the 3-height reach decides none by float error).
CAD_TEXT_H = 2.0                       # label height, ft (the reach, 3 heights, is the plan's original 6 ft)
CAD_OFFSETS_H = (0.5, 2.5)             # the near edge of a label's box off its line, in heights
CAD_SEED = 181                         # which records carry their size in the layer name, label sides and offsets
CAD_KY10_SEED = 182                    # ky10's jitter and label sides and offsets
CAD_JITTER_FT = 0.3
CAD_CHAR_W = 0.6                       # the writer's own width estimate, in heights per character
CAD_LAYER_STYLES = ("WATER-{d:g}IN", "W_{d:g}_{mat}", "WTR-{mm}MM")
CAD_MATERIALS = ("PVC", "DIP", "CI")
CAD_TEXT_STYLES = ('{d:g}" DIP', "{d:g} in", "{mm} mm", "%%c{mm}")
CAD_MTEXT_STYLES = (('\\A1;{{\\fArial|b0;{d:g}" DIP}}', '{d:g}" DIP'), ("{{\\H1.0x;{d:g} in PVC}}", "{d:g} in PVC"))
CAD_MAINS_LAYER, CAD_TEXT_LAYER = "WATER-MAINS", "WATER-TEXT"


def _dxf_v(v) -> str:
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _dxf_ascii(pairs, encoding: str = "cp1252") -> bytes:
    """ASCII DXF: each group code right-aligned in 3 characters on its line, its value on the next."""
    return ("\r\n".join(f"{c:>3}\r\n{_dxf_v(v)}" for c, v in pairs) + "\r\n").encode(encoding)


def _bin_kind(code: int) -> str:
    """A group code's value type in binary DXF (the DXF reference's group code ranges)."""
    if 0 <= code <= 9 or code in (100, 102, 105) or 300 <= code <= 369 or 390 <= code <= 399 or 410 <= code <= 419 \
            or 999 <= code <= 1009:
        return "s"
    if 10 <= code <= 59 or 110 <= code <= 149 or 210 <= code <= 239 or 460 <= code <= 469 or 1010 <= code <= 1059:
        return "d"
    if 60 <= code <= 79 or 170 <= code <= 179 or 270 <= code <= 289 or 370 <= code <= 389 or 400 <= code <= 409 \
            or 1060 <= code <= 1070:
        return "h"
    if 90 <= code <= 99 or 420 <= code <= 429 or 440 <= code <= 449 or code == 1071:
        return "i"
    if 290 <= code <= 299:
        return "b"
    raise ValueError(f"group code {code} not written")


def _dxf_binary_r12(pairs) -> bytes:
    """Binary DXF of R12: the sentinel, then each group code as one byte (255 and two bytes above 254) and its value
    typed by the code (a NUL-ended string, a little-endian double, int16 or int32)."""
    out = [b"AutoCAD Binary DXF\r\n\x1a\x00"]
    for c, v in pairs:
        out.append(bytes([c]) if c < 255 else b"\xff" + struct.pack("<H", c))
        k = _bin_kind(c)
        if k == "s":
            out.append(str(v).encode("cp1252") + b"\x00")
        elif k == "d":
            out.append(struct.pack("<d", float(v)))
        elif k == "h":
            out.append(struct.pack("<h", int(v)))
        elif k == "i":
            out.append(struct.pack("<i", int(v)))
        else:
            out.append(bytes([1 if v else 0]))
    return b"".join(out)


def _r12_header(insunits: int | None = 2, measurement: int | None = None) -> list:
    g = [(0, "SECTION"), (2, "HEADER"), (9, "$ACADVER"), (1, "AC1009")]
    if insunits is not None:
        g += [(9, "$INSUNITS"), (70, insunits)]
    if measurement is not None:
        g += [(9, "$MEASUREMENT"), (70, measurement)]
    return g + [(0, "ENDSEC")]


def _r12_polyline(pts, layer: str, linetype: str | None = None, color: int | None = None) -> list:
    g = [(0, "POLYLINE"), (8, layer)] + ([(6, linetype)] if linetype else []) + ([(62, color)] if color else [])
    g += [(66, 1), (10, 0.0), (20, 0.0), (30, 0.0), (70, 0)]
    for x, y in pts:
        g += [(0, "VERTEX"), (8, layer), (10, float(x)), (20, float(y)), (30, 0.0)]
    return g + [(0, "SEQEND"), (8, layer)]


def _r12_line(a, b, layer: str, linetype: str | None = None, color: int | None = None) -> list:
    return ([(0, "LINE"), (8, layer)] + ([(6, linetype)] if linetype else []) + ([(62, color)] if color else [])
            + [(10, float(a[0])), (20, float(a[1])), (30, 0.0), (11, float(b[0])), (21, float(b[1])), (31, 0.0)])


def _r12_text(ins, h: float, rot: float, text: str, layer: str = CAD_TEXT_LAYER, halign: int = 0, valign: int = 0,
              align=None) -> list:
    g = [(0, "TEXT"), (8, layer), (10, float(ins[0])), (20, float(ins[1])), (30, 0.0), (40, float(h)), (1, text)]
    if rot:
        g.append((50, float(rot)))
    if halign:
        g.append((72, halign))
    if align is not None:
        g += [(11, float(align[0])), (21, float(align[1])), (31, 0.0)]
    if valign:
        g.append((73, valign))
    return g


def _r12_doc(entities: list, insunits: int | None = 2, measurement: int | None = None, tables: list | None = None,
             blocks: list | None = None) -> list:
    g = _r12_header(insunits, measurement)
    if tables:
        g += [(0, "SECTION"), (2, "TABLES")] + tables + [(0, "ENDSEC")]
    if blocks:
        g += [(0, "SECTION"), (2, "BLOCKS")] + blocks + [(0, "ENDSEC")]
    return g + [(0, "SECTION"), (2, "ENTITIES")] + entities + [(0, "ENDSEC"), (0, "EOF")]


def _plain(text: str) -> str:
    """The writer's own plain text of a TEXT value (the %%c code is one character, the diameter symbol)."""
    return text.replace("%%c", "Ø")


def _label_place(pts, n_chars: int, side: int, off_h: float, h: float = CAD_TEXT_H, attach_middle: bool = False):
    """(insertion point, rotation in degrees) of a label at the middle of the line's longest segment, parallel to it,
    reading left to right, on side +1 (left of the reading direction) or -1, the near edge of its box off_h heights
    from the line.  attach_middle: the point is the box's middle left (an MTEXT's attachment point 4), else its lower
    left corner (a TEXT's insertion point)."""
    seg = _seg_len(pts)
    i = int(max(range(len(seg)), key=lambda k: seg[k]))
    (ax, ay), (bx, by) = pts[i], pts[i + 1]
    L = math.hypot(bx - ax, by - ay)
    ux, uy = (bx - ax) / L, (by - ay) / L
    if ux < 0 or (ux == 0 and uy < 0):
        ux, uy = -ux, -uy
    nx, ny = -uy, ux
    w = n_chars * h * CAD_CHAR_W
    sx, sy = (ax + bx) / 2 - ux * w / 2, (ay + by) / 2 - uy * w / 2
    d = off_h * h
    k = (d + h / 2) if attach_middle else (d if side > 0 else d + h)
    k = k if side > 0 else -k
    return (sx + nx * k, sy + ny * k), math.degrees(math.atan2(uy, ux))


@functools.lru_cache(maxsize=2)
def ky4_cad_plan(mtext: bool = False) -> tuple:
    """For each ky4 mains record: {'layer', 'size_in' (the layer's), 'label': None or {'text', 'plain', 'kind',
    'insert', 'rotation', 'side', 'offset_h'}}.  Half the records (floor), chosen by a seeded permutation, carry their
    size in the layer name; the rest are on WATER-MAINS with a label.  mtext: every second label an MTEXT (the R2018
    drawing); the sides and offsets are the same either way."""
    import numpy as np
    recs = ky4_mains_records()
    rng = np.random.default_rng(CAD_SEED)
    perm = rng.permutation(len(recs))
    by_layer = set(int(i) for i in perm[: len(recs) // 2])
    out, kl, kt = [], 0, 0
    for j, r in enumerate(recs):
        d = r["diameter_in"]
        mm = round(d * 25.4)
        if j in by_layer:
            layer = CAD_LAYER_STYLES[kl % 3].format(d=d, mm=mm, mat=CAD_MATERIALS[(kl // 3) % 3])
            kl += 1
            out.append({"layer": layer, "size_in": d, "label": None})
            continue
        side = 1 if int(rng.integers(2)) == 1 else -1
        off = float(rng.uniform(*CAD_OFFSETS_H))
        if mtext and kt % 2 == 1:
            raw, plain = CAD_MTEXT_STYLES[(kt // 2) % 2]
            raw, plain, kind = raw.format(d=d), plain.format(d=d), "MTEXT"
        else:
            raw = CAD_TEXT_STYLES[((kt // 2) if mtext else kt) % 4].format(d=d, mm=mm)
            plain, kind = _plain(raw), "TEXT"
        ins, rot = _label_place(r["pts"], len(plain), side, off, attach_middle=kind == "MTEXT")
        kt += 1
        out.append({"layer": CAD_MAINS_LAYER, "size_in": d,
                    "label": {"text": raw, "plain": plain, "kind": kind, "insert": ins, "rotation": rot, "side": side,
                              "offset_h": off}})
    return tuple(out)


@functools.lru_cache(maxsize=1)
def _ky4_r12_pairs() -> tuple:
    ents = []
    for r, pl in zip(ky4_mains_records(), ky4_cad_plan(False)):
        ents += _r12_polyline(r["pts"], pl["layer"])
        if pl["label"]:
            lb = pl["label"]
            ents += _r12_text(lb["insert"], CAD_TEXT_H, lb["rotation"], lb["text"])
    return tuple(_r12_doc(ents, insunits=2))


def ky4_r12_handwritten_dxf() -> bytes:
    return _dxf_ascii(_ky4_r12_pairs())


def ky4_drawing_dxf_bak() -> bytes:
    """The same bytes under a backup's name."""
    return ky4_r12_handwritten_dxf()


def ky4_binary_dxf() -> bytes:
    """The same drawing as binary DXF, by the fixtures' own binary writer."""
    return _dxf_binary_r12(_ky4_r12_pairs())


def ky4_r2018_lwpolyline_dxf() -> bytes:
    """ky4's mains drawn by ezdxf as R2018 LWPOLYLINEs, every second label an MTEXT (same-library self-consistency)."""
    import ezdxf
    doc = ezdxf.new("R2018")
    doc.header["$INSUNITS"] = 2
    doc.header["$MEASUREMENT"] = 0
    msp = doc.modelspace()
    for r, pl in zip(ky4_mains_records(), ky4_cad_plan(True)):
        msp.add_lwpolyline([tuple(p) for p in r["pts"]], dxfattribs={"layer": pl["layer"]})
        lb = pl["label"]
        if lb and lb["kind"] == "MTEXT":
            msp.add_mtext(lb["text"], dxfattribs={"layer": CAD_TEXT_LAYER, "char_height": CAD_TEXT_H,
                                                  "insert": (*lb["insert"], 0.0), "attachment_point": 4,
                                                  "rotation": lb["rotation"]})
        elif lb:
            msp.add_text(lb["text"], height=CAD_TEXT_H, rotation=lb["rotation"],
                         dxfattribs={"layer": CAD_TEXT_LAYER, "insert": (*lb["insert"], 0.0)})
    buf = io.StringIO()
    doc.write(buf)
    return buf.getvalue().encode("utf-8")


@functools.lru_cache(maxsize=1)
def _ky10_records() -> tuple:
    """ky10's pipes as one line each (its feet), both ends jittered by up to CAD_JITTER_FT, each with a label."""
    import numpy as np
    wn = _wn("ky10")
    rng = np.random.default_rng(CAD_KY10_SEED)
    xy = {n: tuple(map(float, wn.get_node(n).coordinates)) for n in wn.node_name_list}
    recs = []
    for pn in wn.pipe_name_list:
        p = wn.get_link(pn)
        pts = [xy[p.start_node_name]] + [tuple(map(float, v)) for v in p.vertices] + [xy[p.end_node_name]]
        for end in (0, -1):
            dx, dy = rng.uniform(-CAD_JITTER_FT, CAD_JITTER_FT, 2)
            pts[end] = (pts[end][0] + float(dx), pts[end][1] + float(dy))
        d = round(p.diameter / 0.0254, 3)
        side = 1 if int(rng.integers(2)) == 1 else -1
        off = float(rng.uniform(*CAD_OFFSETS_H))
        text = f'{d:g}" DIP'
        ins, rot = _label_place(pts, len(text), side, off)
        recs.append({"pipe": pn, "ends": (p.start_node_name, p.end_node_name), "pts": pts, "diameter_in": d,
                     "label": {"text": text, "insert": ins, "rotation": rot}})
    return tuple(recs)


def ky10_manifest() -> dict:
    """What the scorer needs for ky10's drawing: its pipes (start, end, diameter in) and one record per pipe, and the
    node pairs ky10 draws at one spot (a pump's or a valve's two ends), which no drawing of pipes can separate."""
    from scipy.spatial import cKDTree
    import numpy as np
    wn = _wn("ky10")
    recs = _ky10_records()
    names = wn.node_name_list
    P = np.array([wn.get_node(n).coordinates for n in names], float)
    pairs = sorted(tuple(sorted((names[a], names[b]))) for a, b in cKDTree(P).query_pairs(0.01))
    links = {}
    for ln in wn.pump_name_list + wn.valve_name_list:
        lk = wn.get_link(ln)
        links[tuple(sorted((lk.start_node_name, lk.end_node_name)))] = f"{lk.link_type.lower()} {ln}"
    return {"pipes": {r["pipe"]: (r["ends"][0], r["ends"][1], r["diameter_in"]) for r in recs},
            "records": [{"pipes": (r["pipe"],), "ends": r["ends"], "tee": None} for r in recs],
            "same_spot_pairs": [{"nodes": list(p), "link": links.get(p)} for p in pairs]}


def ky10_r12_dxf() -> bytes:
    ents = []
    for r in _ky10_records():
        ents += _r12_polyline(r["pts"], CAD_MAINS_LAYER)
        lb = r["label"]
        ents += _r12_text(lb["insert"], CAD_TEXT_H, lb["rotation"], lb["text"])
    return _dxf_ascii(_r12_doc(ents, insunits=2))


# labels_offset_cases.dxf: one case a row, 1,000 ft apart; each row's lines, their labels and what each line should
# get ('label', 'propagated' or 'default', and the size).  Offsets are of the near edge of the label's box, in heights.
LABEL_CASES = (
    {"case": "corner_0.5h", "lines": [((0, 0), (300, 0))], "labels": [("corner", 0, 0.5, '4" DIP', 1)],
     "want": [("label", 4.0)]},
    {"case": "corner_2.5h", "lines": [((0, 0), (300, 0))], "labels": [("corner", 0, 2.5, '6" DIP', 1)],
     "want": [("label", 6.0)]},
    {"case": "below_2h", "lines": [((0, 0), (300, 0))], "labels": [("corner", 0, 2.0, "8 in", -1)],
     "want": [("label", 8.0)]},
    {"case": "centred_1h", "lines": [((0, 0), (300, 0))], "labels": [("centre", 0, 1.0, "250 mm", 1)],
     "want": [("label", 250 / 25.4)]},
    {"case": "middle_1.5h", "lines": [((0, 0), (300, 0))], "labels": [("middle", 0, 1.5, '12" PVC', 1)],
     "want": [("label", 12.0)]},
    {"case": "right_2h", "lines": [((0, 0), (300, 0))], "labels": [("right", 0, 2.0, '16"', -1)],
     "want": [("label", 16.0)]},
    {"case": "vertical_1h", "lines": [((150, -150), (150, 150))], "labels": [("corner", 0, 1.0, '3" CI', 1)],
     "want": [("label", 3.0)]},
    {"case": "too_far_3.5h", "lines": [((0, 0), (300, 0))], "labels": [("corner", 0, 3.5, '14" DIP', 1)],
     "want": [("default", None)]},
    {"case": "run_of_three", "lines": [((0, 0), (200, 0)), ((200, 0), (400, 0)), ((400, 0), (600, 0))],
     "labels": [("corner", 0, 1.0, '10" DIP', 1)], "want": [("label", 10.0), ("propagated", 10.0), ("propagated", 10.0)]},
    {"case": "stops_at_a_tee", "lines": [((0, 0), (200, 0)), ((200, 0), (400, 0)), ((200, 0), (200, 200))],
     "labels": [("corner", 0, 1.0, '12" DIP', 1)], "want": [("label", 12.0), ("default", None), ("default", None)]},
    {"case": "nearest_of_two", "lines": [((0, 0), (300, 0)), ((0, -9), (300, -9))],
     "labels": [("corner", 0, 1.0, '6" DIP', 1), ("corner", 1, 1.0, '8" DIP', 1)],
     "want": [("label", 6.0), ("label", 8.0)]},
    {"case": "lone_line", "lines": [((0, 0), (300, 0))], "labels": [], "want": [("default", None)]},
)


def _case_text(kind: str, line, off_h: float, text: str, side: int, h: float = CAD_TEXT_H) -> list:
    """A TEXT label of the line ((x1, y1), (x2, y2)) in one alignment: 'corner' (left, its insertion point the box's
    corner, centred along the line by the writer's estimate), 'centre' (halign 1 at its baseline's middle), 'middle'
    (halign 4 at its box's centre) or 'right' (halign 2 at its baseline's right end)."""
    plain = _plain(text)
    if kind == "corner":
        ins, rot = _label_place(list(line), len(plain), side, off_h, h)
        return _r12_text(ins, h, rot, text)
    (ax, ay), (bx, by) = line
    L = math.hypot(bx - ax, by - ay)
    ux, uy = (bx - ax) / L, (by - ay) / L
    nx, ny = -uy, ux
    mx, my = (ax + bx) / 2, (ay + by) / 2
    d = off_h * h
    w = len(plain) * h * CAD_CHAR_W
    base = (d if side > 0 else d + h) * side          # the baseline's offset from the line
    if kind == "centre":
        p = (mx + nx * base, my + ny * base)
        return _r12_text(p, h, 0.0, text, halign=1, align=p)
    if kind == "middle":
        k = (d + h / 2) * side
        p = (mx + nx * k, my + ny * k)
        return _r12_text(p, h, 0.0, text, halign=4, align=p)
    p = (mx + ux * w / 2 + nx * base, my + uy * w / 2 + ny * base)
    return _r12_text(p, h, 0.0, text, halign=2, align=p)


def labels_offset_cases_dxf() -> bytes:
    ents = []
    for row, c in enumerate(LABEL_CASES):
        y0 = -1000.0 * row
        lines = [((a[0], a[1] + y0), (b[0], b[1] + y0)) for a, b in c["lines"]]
        for a, b in lines:
            ents += _r12_line(a, b, CAD_MAINS_LAYER)
        for kind, li, off, text, side in c["labels"]:
            ents += _case_text(kind, lines[li], off, text, side)
    return _dxf_ascii(_r12_doc(ents, insunits=2))


# xref_block_mains.dxf: mains inside blocks, nested, rotated and scaled, and a block array, and an external xref.
XREF_INSERT = {"at": (5000.0, 3000.0), "rotation": 30.0, "scale": 2.0}
XREF_MINSERT = {"at": (0.0, 0.0), "cols": 2, "rows": 2, "spacing": (500.0, 500.0)}
XREF_MAINS_A = [((0.0, 0.0), (100.0, 0.0)), ((100.0, 0.0), (100.0, 80.0))]
XREF_STUB = ((0.0, 0.0), (0.0, 40.0))
XREF_PATH = "site_mains.dwg"


def xref_block_expected() -> list:
    """The lines a reader must find, in world coordinates (the fixtures' own transform: scale, then rotate, then move),
    each with the layer AutoCAD shows it on and its size from that layer's name."""
    s, a = XREF_INSERT["scale"], math.radians(XREF_INSERT["rotation"])
    X, Y = XREF_INSERT["at"]
    f = lambda p: (X + s * (p[0] * math.cos(a) - p[1] * math.sin(a)), Y + s * (p[0] * math.sin(a) + p[1] * math.cos(a)))
    out = [{"pts": [f(p), f(q)], "layer": "WATER-8IN", "size_in": 8.0} for p, q in XREF_MAINS_A]
    stub0 = f((100.0, 80.0))
    out.append({"pts": [stub0, f((100.0 + XREF_STUB[1][0], 80.0 + XREF_STUB[1][1]))], "layer": "WATER-6IN",
                "size_in": 6.0})
    for r in range(XREF_MINSERT["rows"]):
        for c in range(XREF_MINSERT["cols"]):
            dx, dy = c * XREF_MINSERT["spacing"][0], r * XREF_MINSERT["spacing"][1]
            out.append({"pts": [(dx, dy), (dx + XREF_STUB[1][0], dy + XREF_STUB[1][1])], "layer": "WATER-6IN",
                        "size_in": 6.0})
    out.append({"pts": [(0.0, -100.0), (300.0, -100.0)], "layer": "WATER-8IN", "size_in": 8.0})
    return out


def xref_block_mains_dxf() -> bytes:
    blocks = [(0, "BLOCK"), (8, "0"), (2, "MAINS_A"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, "MAINS_A")]
    for p, q in XREF_MAINS_A:
        blocks += _r12_line(p, q, "0")
    blocks += [(0, "INSERT"), (8, "0"), (2, "STUB"), (10, 100.0), (20, 80.0), (30, 0.0)]
    blocks += [(0, "ENDBLK"), (8, "0")]
    blocks += [(0, "BLOCK"), (8, "0"), (2, "STUB"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, "STUB")]
    blocks += _r12_line(*XREF_STUB, "WATER-6IN") + [(0, "ENDBLK"), (8, "0")]
    blocks += [(0, "BLOCK"), (8, "0"), (2, "XSITE"), (70, 4), (10, 0.0), (20, 0.0), (30, 0.0), (3, "XSITE"),
               (1, XREF_PATH), (0, "ENDBLK"), (8, "0")]
    X, Y = XREF_INSERT["at"]
    ents = [(0, "INSERT"), (8, "WATER-8IN"), (2, "MAINS_A"), (10, X), (20, Y), (30, 0.0),
            (41, XREF_INSERT["scale"]), (42, XREF_INSERT["scale"]), (50, XREF_INSERT["rotation"])]
    mx, my = XREF_MINSERT["at"]
    ents += [(0, "INSERT"), (8, "WATER-4IN"), (2, "STUB"), (10, mx), (20, my), (30, 0.0),
             (70, XREF_MINSERT["cols"]), (71, XREF_MINSERT["rows"]), (44, XREF_MINSERT["spacing"][0]),
             (45, XREF_MINSERT["spacing"][1])]
    ents += [(0, "INSERT"), (8, "WATER-8IN"), (2, "XSITE"), (10, 0.0), (20, 0.0), (30, 0.0)]
    ents += _r12_line((0.0, -100.0), (300.0, -100.0), "WATER-8IN")
    return _dxf_ascii(_r12_doc(ents, insunits=2, blocks=blocks))


def _r2000(entities: list, insunits: int | None = 2, measurement: int | None = None, classes: list | None = None,
           version: str = "AC1015", encoding: str = "cp1252") -> bytes:
    """A minimal R2000 ASCII DXF (header, classes, entities with handles and subclass markers), written by hand; with
    version AC1032 an R2018 one (whose text a reader takes as UTF-8), in the given encoding."""
    g = [(0, "SECTION"), (2, "HEADER"), (9, "$ACADVER"), (1, version), (9, "$HANDSEED"), (5, "FFFFF")]
    if version >= "AC1021":           # R2007 on: UTF-8 by its version, given beside the code page as AutoCAD writes it
        g += [(9, "$DWGCODEPAGE"), (3, "ANSI_1252")]
    if insunits is not None:
        g += [(9, "$INSUNITS"), (70, insunits)]
    if measurement is not None:
        g += [(9, "$MEASUREMENT"), (70, measurement)]
    g += [(0, "ENDSEC")]
    if classes:
        g += [(0, "SECTION"), (2, "CLASSES")] + classes + [(0, "ENDSEC")]
    g += [(0, "SECTION"), (2, "ENTITIES")]
    for k, ent in enumerate(entities):
        g += [ent[0], (5, f"{0x100 + k:X}")] + list(ent[1:])
    return _dxf_ascii(g + [(0, "ENDSEC"), (0, "EOF")], encoding)


def _r2000_line(a, b, layer: str) -> list:
    return [(0, "LINE"), (100, "AcDbEntity"), (8, layer), (100, "AcDbLine"), (10, float(a[0])), (20, float(a[1])),
            (30, 0.0), (11, float(b[0])), (21, float(b[1])), (31, 0.0)]


def _r2000_mtext(ins, h: float, text: str, attach: int, layer: str = CAD_TEXT_LAYER) -> list:
    return [(0, "MTEXT"), (100, "AcDbEntity"), (8, layer), (100, "AcDbMText"), (10, float(ins[0])), (20, float(ins[1])),
            (30, 0.0), (40, float(h)), (71, attach), (1, text)]


CIVIL3D_CLASSES = [(0, "CLASS"), (1, "AECC_PRESSURE_PIPE"), (2, "AeccDbPressurePipe"), (3, "AeccDbPipeNetwork"),
                   (90, 1153), (91, 3), (280, 1), (281, 1)]


def _civil3d_entities() -> list:
    ents = []
    for k in range(3):
        ents.append([(0, "ACAD_PROXY_ENTITY"), (100, "AcDbEntity"), (8, "C-WATR-PIPE"), (100, "AcDbProxyEntity"),
                     (90, 498), (91, 500), (95, 2), (70, 0), (92, 0), (93, 0)])
    for k in range(2):
        ents.append([(0, "AECC_PRESSURE_PIPE"), (100, "AcDbEntity"), (8, "C-WATR-PIPE"), (100, "AeccDbPressurePipe"),
                     (10, 100.0 * k), (20, 0.0)])
    return ents


def civil3d_proxy_dxf() -> bytes:
    """A Civil 3D pressure network saved to DXF as it is: proxy entities and AECC objects only."""
    return _r2000(_civil3d_entities(), classes=CIVIL3D_CLASSES)


def civil3d_proxy_with_lines_dxf() -> bytes:
    """The same with three plain lines beside it (a drawing that mixes them)."""
    lines = [_r2000_line((0.0, -500.0), (300.0, -500.0), "W-8IN"), _r2000_line((300.0, -500.0), (600.0, -500.0), "W-8IN"),
             _r2000_line((300.0, -500.0), (300.0, -800.0), "W-8IN")]
    return _r2000(_civil3d_entities() + lines, classes=CIVIL3D_CLASSES)


# abandoned_linetype.dxf: a small network on WATER-MAINS, with two abandoned mains drawn dashed on the same layer
ABANDONED_ACTIVE = [((0, 0), (400, 0)), ((400, 0), (800, 0)), ((0, 0), (0, 400)), ((0, 400), (400, 400)),
                    ((400, 400), (800, 400)), ((800, 0), (800, 400))]
ABANDONED_DASHED = [((400, 0), (400, 400)), ((0, 400), (0, 800))]


def abandoned_linetype_dxf() -> bytes:
    tables = [(0, "TABLE"), (2, "LTYPE"), (70, 2),
              (0, "LTYPE"), (2, "CONTINUOUS"), (70, 0), (3, "Solid line"), (72, 65), (73, 0), (40, 0.0),
              (0, "LTYPE"), (2, "DASHED"), (70, 0), (3, "Dashed __ __ __"), (72, 65), (73, 2), (40, 0.75), (49, 0.5),
              (49, -0.25), (0, "ENDTAB"),
              (0, "TABLE"), (2, "LAYER"), (70, 2),
              (0, "LAYER"), (2, "0"), (70, 0), (62, 7), (6, "CONTINUOUS"),
              (0, "LAYER"), (2, CAD_MAINS_LAYER), (70, 0), (62, 5), (6, "CONTINUOUS"), (0, "ENDTAB")]
    ents = []
    for a, b in ABANDONED_ACTIVE:
        ents += _r12_line(a, b, CAD_MAINS_LAYER)
    for a, b in ABANDONED_DASHED:
        ents += _r12_line(a, b, CAD_MAINS_LAYER, linetype="DASHED", color=8)
    ents += _case_text("corner", ABANDONED_ACTIVE[0], 1.0, '8" DIP', 1)
    return _dxf_ascii(_r12_doc(ents, insunits=2, tables=tables))


def _units_lines(n: int = 60) -> list:
    """The first n ky4 mains records (feet), each with its size in the layer name (WATER-<d>IN)."""
    ents = []
    for r in ky4_mains_records()[:n]:
        ents += _r12_polyline(r["pts"], CAD_LAYER_STYLES[0].format(d=r["diameter_in"]))
    return ents


def insunits_zero_dxf() -> bytes:
    return _dxf_ascii(_r12_doc(_units_lines(), insunits=0))


def insunits_wrong_dxf() -> bytes:
    """Metres declared, feet drawn: a drawing from an imperial template ($MEASUREMENT 0) whose units say metres."""
    return _dxf_ascii(_r12_doc(_units_lines(), insunits=6, measurement=0))


def insunits_wrong_nomeasure_dxf() -> bytes:
    """The same with no measurement setting: the mistake the reader cannot see (a known limit, recorded)."""
    return _dxf_ascii(_r12_doc(_units_lines(), insunits=6))


def insunits_mm_dxf() -> bytes:
    """Millimetres declared (a metric template's default), feet drawn: the median line is then 0.5 m long."""
    return _dxf_ascii(_r12_doc(_units_lines(), insunits=4, measurement=1))


MTEXT_CASES = ((((0.0, 0.0), (300.0, 0.0)), (100.0, 2.0), 7, '\\A1;{\\fArial|b0;8" DIP}', 8.0),
               (((0.0, -1000.0), (300.0, -1000.0)), (100.0, -1002.0), 1, "{\\H1.0x;6 in PVC}", 6.0),
               (((0.0, -2000.0), (300.0, -2000.0)), (150.0, -1997.0), 5, "\\pxqc;{\\C1;\\L12\\l\" DIP}", 12.0))


def mtext_codes_dxf() -> bytes:
    """Three lines, each labelled by an MTEXT with formatting codes (bottom left, top left and middle centre attachment,
    the near edge of each box 1 height off its line)."""
    ents = []
    for line, ins, attach, text, _ in MTEXT_CASES:
        ents.append(_r2000_line(*line, CAD_MAINS_LAYER))
        ents.append(_r2000_mtext(ins, CAD_TEXT_H, text, attach))
    return _r2000(ents)


# the review of task 18: one drawing for each case it found, written by these writers (none by ezdxf)
def _ky4_entities() -> list:
    """ky4's R12 drawing's entities (its mains and labels), as _ky4_r12_pairs writes them."""
    pairs = list(_ky4_r12_pairs())
    start = next(i for i, p in enumerate(pairs) if p == (2, "ENTITIES")) + 1
    return pairs[start:-2]


def ky4_inches_declared_dxf() -> bytes:
    """ky4's drawing, in feet, whose insertion units say inches (INSUNITS 1, the unit AutoCAD's imperial template
    starts with): read in inches every length is twelve times too short, yet the median line (42 ft) is plausible."""
    return _dxf_ascii(_r12_doc(_ky4_entities(), insunits=1, measurement=0))


REVIEW_INCH_LINES = [((0.0, 0.0), (7200.0, 0.0)), ((7200.0, 0.0), (14400.0, 0.0)), ((7200.0, 0.0), (7200.0, 7200.0)),
                     ((7200.0, 7200.0), (14400.0, 7200.0))]       # 600 ft mains drawn in inches


def inches_true_dxf() -> bytes:
    """Four 600 ft mains drawn in inches (INSUNITS 1): read in feet the median line would be 7,200 ft, implausible, so
    inches are taken without a question."""
    ents = []
    for a, b in REVIEW_INCH_LINES:
        ents += _r12_line(a, b, "WATER-8IN")
    return _dxf_ascii(_r12_doc(ents, insunits=1, measurement=0))


ZONE_LAYERS = ("WATER-ZONE-2", "WATER-ZONE-3", "PVC-ZONE-3", "WATER-ZONE-2-PVC")


def _zone_label(layer: str) -> str:
    return '8" PVC' if "PVC" in layer else '8" DIP'


def zone_layers_dxf() -> bytes:
    """Nineteen mains on four layers whose names carry a zone number (WATER-ZONE-2 and -3, and since the final
    verification PVC-ZONE-3 and WATER-ZONE-2-PVC, whose material word once cued the zone number as a size), in four
    rows joined by three mains up the left side, each main with an 8 in label beside it: the zone number is no size
    (the review's zones.dxf)."""
    ents = []
    for k, lay in enumerate(ZONE_LAYERS):
        y = 600.0 * k
        for i in range(4):
            ents += _r12_line((300.0 * i, y), (300.0 * (i + 1), y), lay)
            ents += _r12_text((300.0 * i + 100, y + 1.0), CAD_TEXT_H, 0.0, _zone_label(lay))
    for k, lay in enumerate(ZONE_LAYERS[:-1]):
        ents += _r12_line((0.0, 600.0 * k), (0.0, 600.0 * (k + 1)), lay)
        ents += _r12_text((1.0, 600.0 * k + 300.0), CAD_TEXT_H, 90.0, _zone_label(lay))
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def latin1_r2018_dxf() -> bytes:
    """An R2018 drawing (UTF-8 by its version) whose writer put cp1252 bytes in a layer name: three mains on
    'TUBERÍA-8IN' as the byte 0xCD, invalid in UTF-8."""
    lay = "TUBERÍA-8IN"
    ents = [_r2000_line(a, b, lay) for a, b in (((0.0, 0.0), (300.0, 0.0)), ((300.0, 0.0), (300.0, 300.0)),
                                                ((300.0, 300.0), (600.0, 300.0)))]
    return _r2000(ents, measurement=0, version="AC1032", encoding="cp1252")


def _handled_line(a, b, handle: str, layer: str = "WATER-8IN") -> list:
    return [(0, "LINE"), (5, handle), (8, layer), (10, float(a[0])), (20, float(a[1])), (30, 0.0), (11, float(b[0])),
            (21, float(b[1])), (31, 0.0)]


def dup_handles_dxf() -> bytes:
    """An R12 chain of ten mains from a writer that reuses handles: each pair of lines shares one (ezdxf's audit
    removed the second of each, half the mains, without a word)."""
    ents = []
    for i in range(10):
        ents += _handled_line((300.0 * i, 0.0), (300.0 * (i + 1), 0.0), f"{0x20 + i // 2:X}")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def zero_handles_dxf() -> bytes:
    """Three mains whose handles are all 0 (refused outright as damaged before the review)."""
    ents = []
    for a, b in (((0.0, 0.0), (300.0, 0.0)), ((300.0, 0.0), (300.0, 300.0)), ((300.0, 300.0), (600.0, 300.0))):
        ents += _handled_line(a, b, "0")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


SPLINE_FRAME = ((0.0, 0.0), (300.0, 300.0), (600.0, 0.0))


def spline_fit_curve() -> list:
    """The fitted vertices of spline_fit_polyline.dxf: a parabola from (0, 0) to (600, 0), 150 ft high."""
    return [(0.0, 0.0)] + [(600.0 * t, 600.0 * t * (1 - t)) for t in (k / 10 for k in range(1, 10))] + [(600.0, 0.0)]


def spline_fit_polyline_dxf() -> bytes:
    """A spline-fit 2D POLYLINE as AutoCAD keeps it (PEDIT Spline): polyline flag 4, its control frame's three
    vertices (flag 16), which AutoCAD does not draw, then the fitted curve's (flag 8), and a plain main after it."""
    g = [(0, "POLYLINE"), (8, "WATER-8IN"), (66, 1), (10, 0.0), (20, 0.0), (30, 0.0), (70, 4), (75, 6)]
    for x, y in SPLINE_FRAME:
        g += [(0, "VERTEX"), (8, "WATER-8IN"), (10, x), (20, y), (30, 0.0), (70, 16)]
    for x, y in spline_fit_curve():
        g += [(0, "VERTEX"), (8, "WATER-8IN"), (10, float(x)), (20, float(y)), (30, 0.0), (70, 8)]
    g += [(0, "SEQEND"), (8, "WATER-8IN")]
    g += _r12_line((600.0, 0.0), (900.0, 0.0), "WATER-8IN")
    return _dxf_ascii(_r12_doc(g, insunits=2, measurement=0))


def _layer_row(name: str, flags: int = 0, color: int = 7) -> list:
    return [(0, "LAYER"), (2, name), (70, flags), (62, color), (6, "CONTINUOUS")]


def frozen_layers_dxf() -> bytes:
    """Layers AutoCAD hides: a main on 'old-mains' (the table names OLD-MAINS, frozen: layer names are
    case-insensitive), an old label '12\" CI' on the frozen OLD-TEXT beside an unsized main, and a block of three 12 in
    mains inserted on the frozen XREF-PROPOSED (a bound proposal, frozen), beside three 8 in mains in service."""
    tables = [(0, "TABLE"), (2, "LAYER"), (70, 5)] + _layer_row("WATER-8IN") + _layer_row("OLD-MAINS", 1) + \
        _layer_row("OLD-TEXT", 1) + _layer_row("XREF-PROPOSED", 1) + _layer_row("PROPOSED$0$WATER-12IN") + \
        [(0, "ENDTAB")]
    blk = [(0, "BLOCK"), (8, "0"), (2, "PROPOSED"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, "PROPOSED")]
    for i in range(3):
        blk += _r12_line((300.0 * i, -500.0), (300.0 * (i + 1), -500.0), "PROPOSED$0$WATER-12IN")
    blk += [(0, "ENDBLK"), (8, "0")]
    ents = []
    for a, b in (((0.0, 0.0), (300.0, 0.0)), ((300.0, 0.0), (600.0, 0.0)), ((0.0, 0.0), (0.0, -500.0))):
        ents += _r12_line(a, b, "WATER-8IN")
    ents += _r12_line((600.0, 0.0), (900.0, 0.0), "old-mains")
    ents += _r12_line((0.0, 0.0), (0.0, 300.0), "WATER")
    ents += _r12_text((2.0, 100.0), CAD_TEXT_H, 90.0, '12" CI', "OLD-TEXT")
    ents += [(0, "INSERT"), (8, "XREF-PROPOSED"), (2, "PROPOSED"), (10, 0.0), (20, 0.0), (30, 0.0)]
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0, tables=tables, blocks=blk))


CUT_SHORT_SHARES = (0.99, 0.9, 0.5)


def ky4_cut_short_dxf(share: float = 0.9) -> bytes:
    """ky4's drawing cut to a share of its bytes, as after an interrupted copy or upload (no '0 EOF')."""
    data = ky4_r12_handwritten_dxf()
    return data[:int(len(data) * share)]


def ky4_no_eof_dxf() -> bytes:
    """ky4's drawing with only its last '0 ENDSEC 0 EOF' missing."""
    data = ky4_r12_handwritten_dxf()
    cut = data.rindex(b"  0\r\nENDSEC")
    return data[:cut]


def diameter_sign_dxf() -> bytes:
    """Labels with AutoCAD's diameter sign U+2205: a TEXT '∅8' and an MTEXT '∅12', each beside its main, in an
    R2018 drawing (UTF-8)."""
    ents = [_r2000_line((0.0, 0.0), (300.0, 0.0), CAD_MAINS_LAYER), _r2000_line((300.0, 0.0), (600.0, 0.0), CAD_MAINS_LAYER),
            [(0, "TEXT"), (100, "AcDbEntity"), (8, CAD_TEXT_LAYER), (100, "AcDbText"), (10, 100.0), (20, 1.0), (30, 0.0),
             (40, CAD_TEXT_H), (1, "\u22058"), (100, "AcDbText")],
            _r2000_mtext((400.0, 3.0), CAD_TEXT_H, "\u220512", 7)]
    return _r2000(ents, measurement=0, version="AC1032", encoding="utf-8")


def c_from_label_dxf() -> bytes:
    """Two mains on WATER-8IN, each with an '8\" CI' label: the layer gives the size, the label the material."""
    ents = []
    for i in range(2):
        ents += _r12_line((300.0 * i, 0.0), (300.0 * (i + 1), 0.0), "WATER-8IN")
        ents += _r12_text((300.0 * i + 100, 1.0), CAD_TEXT_H, 0.0, '8" CI')
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def no_sizes_dxf() -> bytes:
    """Five mains on WATER-MAINS with no size anywhere: the size of every pipe is asked."""
    ents = []
    for a, b in (((0.0, 0.0), (300.0, 0.0)), ((300.0, 0.0), (600.0, 0.0)), ((300.0, 0.0), (300.0, 300.0)),
                 ((300.0, 300.0), (600.0, 300.0)), ((600.0, 0.0), (600.0, 300.0))):
        ents += _r12_line(a, b, CAD_MAINS_LAYER)
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


LONG_LABEL_CHARS = 3_000


def long_label_dxf() -> bytes:
    """A main and a label of 3,000 characters (over cad_read's cap), which is not read: the drawing says so."""
    ents = _r12_line((0.0, 0.0), (300.0, 0.0), "WATER-8IN") + _r12_line((300.0, 0.0), (600.0, 0.0), "WATER-8IN")
    ents += _r12_text((100.0, 1.0), CAD_TEXT_H, 0.0, "x" * (LONG_LABEL_CHARS - 4) + "8 in")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


REVIEW_18_FIXTURES = ("ky4_inches_declared.dxf", "inches_true.dxf", "zone_layers.dxf", "latin1_r2018.dxf",
                      "dup_handles.dxf", "zero_handles.dxf", "spline_fit_polyline.dxf", "frozen_layers.dxf",
                      *(f"ky4_cut_{round(100 * s)}.dxf" for s in CUT_SHORT_SHARES), "ky4_no_eof.dxf",
                      "diameter_sign.dxf", "c_from_label.dxf", "no_sizes.dxf", "long_label.dxf")


# hostile drawings, each within the hosted demo's 5 MB upload limit: each must be read, or refused, in time linear in
# the file (cad_read's text cap, block weights and limits)
HOSTILE_MTEXT_CHARS = 4_800_000        # one MTEXT of formatting codes (ezdxf's fast reader: outputs/intake/cad_mtext_timing.json)
HOSTILE_STACK_TEXTS = (1_000, 1_990)   # MTEXTs of '\\S1^2;' stacks, each just under cad_read's text cap
HOSTILE_BOMB = (9, 10)                 # blocks nested 9 deep, each but the last inserting the next 10 times: 10**8 lines
HOSTILE_MINSERT = (1000, 1000)         # one block array of a million copies
HOSTILE_SPLINE_POINTS = 60_000         # control points of one spline
HOSTILE_LAYER_CHARS = 1_000_000        # a layer name, and a TEXT of as many spaces before its size


def hostile_mtext_codes_dxf() -> bytes:
    ents = [_r2000_line((0.0, 0.0), (300.0, 0.0), CAD_MAINS_LAYER),
            _r2000_mtext((100.0, 2.0), CAD_TEXT_H, "\\A1;" * (HOSTILE_MTEXT_CHARS // 4), 7)]
    return _r2000(ents)


def hostile_stacked_mtext_dxf() -> bytes:
    n, chars = HOSTILE_STACK_TEXTS
    ents = [_r2000_line((0.0, 0.0), (300.0, 0.0), CAD_MAINS_LAYER)]
    ents += [_r2000_mtext((float(k), 2.0), CAD_TEXT_H, "\\S1^2;" * (chars // 6), 7) for k in range(n)]
    return _r2000(ents)


def hostile_block_bomb_dxf() -> bytes:
    depth, fan = HOSTILE_BOMB
    blocks = []
    for k in range(depth):
        blocks += [(0, "BLOCK"), (8, "0"), (2, f"B{k}"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, f"B{k}")]
        if k == depth - 1:
            blocks += _r12_line((0.0, 0.0), (10.0, 0.0), "0")
        else:
            for j in range(fan):
                blocks += [(0, "INSERT"), (8, "0"), (2, f"B{k + 1}"), (10, float(j)), (20, 0.0), (30, 0.0)]
        blocks += [(0, "ENDBLK"), (8, "0")]
    ents = [(0, "INSERT"), (8, CAD_MAINS_LAYER), (2, "B0"), (10, 0.0), (20, 0.0), (30, 0.0)]
    return _dxf_ascii(_r12_doc(ents, insunits=2, blocks=blocks))


def hostile_minsert_dxf() -> bytes:
    cols, rows = HOSTILE_MINSERT
    blocks = [(0, "BLOCK"), (8, "0"), (2, "S"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, "S")]
    blocks += _r12_line((0.0, 0.0), (10.0, 0.0), "0") + [(0, "ENDBLK"), (8, "0")]
    ents = [(0, "INSERT"), (8, CAD_MAINS_LAYER), (2, "S"), (10, 0.0), (20, 0.0), (30, 0.0), (70, cols), (71, rows),
            (44, 20.0), (45, 20.0)]
    ents += _r12_line((0.0, -100.0), (300.0, -100.0), "WATER-8IN") + _r12_line((300.0, -100.0), (600.0, -100.0),
                                                                                  "WATER-8IN")
    return _dxf_ascii(_r12_doc(ents, insunits=2, blocks=blocks))


def hostile_spline_dxf() -> bytes:
    n = HOSTILE_SPLINE_POINTS
    sp = [(0, "SPLINE"), (100, "AcDbEntity"), (8, CAD_MAINS_LAYER), (100, "AcDbSpline"), (70, 8), (71, 3),
          (72, n + 4), (73, n), (74, 0)]
    sp += [(40, float(max(0, min(k - 3, n - 3)))) for k in range(n + 4)]
    for k in range(n):
        sp += [(10, float(k)), (20, float(k % 7)), (30, 0.0)]
    ents = [sp, _r2000_line((0.0, -100.0), (300.0, -100.0), "WATER-8IN"),
            _r2000_line((300.0, -100.0), (600.0, -100.0), "WATER-8IN")]
    return _r2000(ents)


def hostile_long_names_dxf() -> bytes:
    n = HOSTILE_LAYER_CHARS
    ents = _r12_line((0.0, 0.0), (300.0, 0.0), "W" * n) + _r12_line((300.0, 0.0), (600.0, 0.0), "WATER-8IN")
    ents += _r12_text((100.0, 2.0), CAD_TEXT_H, 0.0, " " * n + '8"')
    ents += _r12_text((100.0, 4.0), CAD_TEXT_H, 0.0, "(" * (n // 2) + "8 in")
    return _dxf_ascii(_r12_doc(ents, insunits=2))


# the review of task 18: hostile drawings it found the reader or the build unbounded on
HOSTILE_TALL = (600, 50, 2_000, 1e7)        # 600 lines of 50 points, and 2,000 labels '8 in' of height 10**7 ft
HOSTILE_CROWDED = (600, 50, 2_000, 100.0)   # the same lines 1 ft apart, and 2,000 labels 100 ft high among them
HOSTILE_NEG_MINSERT = (8, 10)               # blocks nested 8 deep, each inserting the next 10 times and once as a
                                            # MINSERT of -10 rows (a negative count once read as negative copies)
HOSTILE_BULGES = (120_000, 100.0)           # one LWPOLYLINE of 120,000 vertices, each segment a near circle (bulge 100)


def _grid_mains(n_lines: int, n_pts: int, dx: float, dy: float, layer: str = "WATER-8IN") -> list:
    ents = []
    for i in range(n_lines):
        ents += _r12_polyline([(float(j * dx), float(i * dy)) for j in range(n_pts)], layer)
    return ents


def hostile_tall_labels_dxf() -> bytes:
    """Within the hosted demo's line and point limits: 600 lines of 50 points, and 2,000 labels of a size whose height
    makes each reach every line (the review measured seconds and gigabytes for a few hundred such labels)."""
    n, p, k, h = HOSTILE_TALL
    ents = _grid_mains(n, p, 10.0, 40.0)
    for _ in range(k):
        ents += _r12_text((0.0, 0.0), h, 0.0, "8 in", "T")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def hostile_crowded_labels_dxf() -> bytes:
    """600 lines of 50 points, 1 ft apart, and 2,000 labels 100 ft high among them, a plausible height: every label
    reaches every segment."""
    n, p, k, h = HOSTILE_CROWDED
    ents = _grid_mains(n, p, 1.0, 1.0)
    for i in range(k):
        ents += _r12_text((float(i % p), float(i % n)), h, 0.0, "8 in", "T")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def hostile_negative_minsert_dxf() -> bytes:
    depth, fan = HOSTILE_NEG_MINSERT
    blocks = []
    for k in range(1, depth + 1):
        blocks += [(0, "BLOCK"), (8, "0"), (2, f"B{k}"), (70, 0), (10, 0.0), (20, 0.0), (30, 0.0), (3, f"B{k}")]
        if k == depth:
            for i in range(fan):
                blocks += _r12_line((0.0, float(i)), (300.0, float(i)), "0")
        else:
            blocks += [(0, "INSERT"), (8, "0"), (2, f"B{k + 1}"), (10, 0.0), (20, 0.0), (30, 0.0), (70, 1),
                       (71, -fan), (44, 1.0), (45, 1.0)]
            for _ in range(fan):
                blocks += [(0, "INSERT"), (8, "0"), (2, f"B{k + 1}"), (10, 0.0), (20, 0.0), (30, 0.0)]
        blocks += [(0, "ENDBLK"), (8, "0")]
    ents = [(0, "INSERT"), (8, "WATER-8IN"), (2, "B1"), (10, 0.0), (20, 0.0), (30, 0.0)]
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0, blocks=blocks))


def hostile_bulges_dxf() -> bytes:
    n, b = HOSTILE_BULGES
    lw = [(0, "LWPOLYLINE"), (100, "AcDbEntity"), (8, "WATER-8IN"), (100, "AcDbPolyline"), (90, n), (70, 0)]
    for k in range(n):
        lw += [(10, float(10 * k)), (20, 0.0), (42, b)]
    return _r2000([lw], measurement=0)


# the final verification of task 18: what its fixes still left open
HOSTILE_LONG_HANDLE = (100_000, 3_000)      # one handle of 100,000 hex digits, and 3,000 lines whose handle is 0
HOSTILE_INF_BULGES = 160_000                # one LWPOLYLINE of 160,000 vertices at infinity, each with bulge 1
NEAR_LIMIT = (190, 10, 10, 4.0, 2.0)        # rows, polylines a row, segments a polyline, segment length, row spacing (ft)
NEAR_LIMIT_LABELS = (1_000, 100.0, 7)       # labels '8 in', their height (ft), the seed of their places
DISPUTED_LAYERS = ("W-2-PVC", "W-8-PVC")    # a cued zone number its labels dispute; a cued size one stray label does not


def hostile_long_handle_dxf() -> bytes:
    """One main whose handle is HOSTILE_LONG_HANDLE[0] hex digits long, one whose handle is the largest a 64-bit handle
    can be, and HOSTILE_LONG_HANDLE[1] whose handle is 0: each renumbered handle was written as long as the long one
    (the final verification: 0.4 MB repaired to 300 MB)."""
    n_hex, n_zero = HOSTILE_LONG_HANDLE
    ents = _handled_line((0.0, 0.0), (300.0, 0.0), "F" * n_hex) + \
        _handled_line((300.0, 0.0), (600.0, 0.0), "FFFFFFFFFFFFFFFF")
    for i in range(n_zero):
        ents += _handled_line((300.0 * (i + 2), 0.0), (300.0 * (i + 3), 0.0), "0")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def hostile_inf_bulges_dxf() -> bytes:
    """An R2000 drawing within the hosted demo's upload limit: one LWPOLYLINE of HOSTILE_INF_BULGES vertices, all at
    (inf, inf) with bulge 1, beside two plain mains (the final verification: its arcs were charged no points, and
    made 5.8 million)."""
    lw = [(0, "LWPOLYLINE"), (100, "AcDbEntity"), (8, "WATER-8IN"), (100, "AcDbPolyline"), (90, HOSTILE_INF_BULGES),
          (70, 0)]
    for _ in range(HOSTILE_INF_BULGES):
        lw += [(10, "inf"), (20, "inf"), (42, 1.0)]
    return _r2000([lw, _r2000_line((0.0, 0.0), (300.0, 0.0), "WATER-8IN"),
                   _r2000_line((300.0, 0.0), (600.0, 0.0), "WATER-8IN")], measurement=0)


HOSTILE_FIXTURES = ("hostile_mtext_codes.dxf", "hostile_stacked_mtext.dxf", "hostile_block_bomb.dxf",
                    "hostile_minsert.dxf", "hostile_spline.dxf", "hostile_long_names.dxf",
                    "hostile_tall_labels.dxf", "hostile_crowded_labels.dxf", "hostile_negative_minsert.dxf",
                    "hostile_bulges.dxf", "hostile_long_handle.dxf", "hostile_inf_bulges.dxf")


def label_pairs_near_limit_dxf() -> bytes:
    """Just under the local limit of pairs of a label and a segment (cad_read.MAX_LABEL_PAIRS): NEAR_LIMIT[0] rows of
    NEAR_LIMIT[1] polylines of NEAR_LIMIT[2] segments each, end to end, and NEAR_LIMIT_LABELS[0] labels '8 in' 100 ft
    high (a height the build takes as a label's) among them, each reaching every segment: 19,000,000 pairs (the final
    verification measured 5.2 GB for as many pairs over separate lines, and 1.3 GB over long polylines)."""
    import numpy as np
    rows, per_row, segs, seg_ft, dy = NEAR_LIMIT
    n, h, seed = NEAR_LIMIT_LABELS
    ents = []
    for j in range(rows):
        for k in range(per_row):
            x0 = k * segs * seg_ft
            ents += _r12_polyline([(x0 + seg_ft * i, dy * j) for i in range(segs + 1)], CAD_MAINS_LAYER)
    rng = np.random.default_rng(seed)
    for x, y in rng.uniform(0.0, 200.0, (n, 2)):
        ents += _r12_text((100.0 + float(x), 90.0 + float(y)), h, 0.0, "8 in", "T")
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def disputed_layer_dxf() -> bytes:
    """Two rows of eight mains joined by one: on W-2-PVC, whose 2 the material beside it cues as a size, three mains
    carry '8\" PVC' labels (a zone or phase number, disputed by its labels: asked); on W-8-PVC one main carries a
    stray '12\" PVC' label (one label is not enough: the layer's 8 in is kept, and warned, as on ky4's drawing)."""
    ents = []
    for k, lay in enumerate(DISPUTED_LAYERS):
        y = 600.0 * k
        for i in range(8):
            ents += _r12_line((300.0 * i, y), (300.0 * (i + 1), y), lay)
        for i in ((0, 3, 6) if k == 0 else (4,)):
            ents += _r12_text((300.0 * i + 100, y + 1.0), CAD_TEXT_H, 0.0, '8" PVC' if k == 0 else '12" PVC')
    ents += _r12_line((0.0, 0.0), (0.0, 600.0), DISPUTED_LAYERS[1])
    return _dxf_ascii(_r12_doc(ents, insunits=2, measurement=0))


def bom_dup_handles_dxf() -> bytes:
    """dup_handles.dxf after a UTF-8 byte order mark (the final verification: the repair skipped such a file)."""
    return b"\xef\xbb\xbf" + dup_handles_dxf()


def bom_ky4_no_eof_dxf() -> bytes:
    """ky4_no_eof.dxf after a UTF-8 byte order mark."""
    return b"\xef\xbb\xbf" + ky4_no_eof_dxf()


def binary_dup_handles_dxf() -> bytes:
    """dup_handles.dxf's ten mains as a binary DXF (the fixtures' own binary writer): the repair is for ASCII DXF only,
    so ezdxf's audit removes half of them, and the reader warns (a recorded limit)."""
    pairs = _r12_header(2, 0) + [(0, "SECTION"), (2, "ENTITIES")]
    for i in range(10):
        pairs += _handled_line((300.0 * i, 0.0), (300.0 * (i + 1), 0.0), f"{0x20 + i // 2:X}")
    return _dxf_binary_r12(pairs + [(0, "ENDSEC"), (0, "EOF")])


FINAL_18_FIXTURES = ("disputed_layer.dxf", "bom_dup_handles.dxf", "bom_ky4_no_eof.dxf", "binary_dup_handles.dxf",
                     "label_pairs_near_limit.dxf")
CAD_FIXTURES = ("ky4_r12_handwritten.dxf", "ky4_r2018_lwpolyline.dxf", "ky4_drawing.dxf.bak", "ky4_binary.dxf",
                "ky10_r12.dxf", "labels_offset_cases.dxf", "xref_block_mains.dxf", "civil3d_proxy.dxf",
                "civil3d_proxy_with_lines.dxf", "abandoned_linetype.dxf", "insunits_zero.dxf", "insunits_wrong.dxf",
                "insunits_wrong_nomeasure.dxf", "insunits_mm.dxf", "mtext_codes.dxf",
                "dwg_ac1027.dwg", "drawing.bak", "dgn_v7.dgn", "ole_dgn_v8.dgn", "plant.ifc") + HOSTILE_FIXTURES \
    + REVIEW_18_FIXTURES + FINAL_18_FIXTURES
CAD_GUIDANCE = ("dwg_ac1027.dwg", "drawing.bak", "dgn_v7.dgn", "ole_dgn_v8.dgn", "plant.ifc")
FIXTURES.update({"ky4_r12_handwritten.dxf": ky4_r12_handwritten_dxf, "ky4_r2018_lwpolyline.dxf": ky4_r2018_lwpolyline_dxf,
                 "ky4_drawing.dxf.bak": ky4_drawing_dxf_bak, "ky4_binary.dxf": ky4_binary_dxf,
                 "ky10_r12.dxf": ky10_r12_dxf, "labels_offset_cases.dxf": labels_offset_cases_dxf,
                 "xref_block_mains.dxf": xref_block_mains_dxf, "civil3d_proxy.dxf": civil3d_proxy_dxf,
                 "civil3d_proxy_with_lines.dxf": civil3d_proxy_with_lines_dxf,
                 "abandoned_linetype.dxf": abandoned_linetype_dxf, "insunits_zero.dxf": insunits_zero_dxf,
                 "insunits_wrong.dxf": insunits_wrong_dxf, "insunits_wrong_nomeasure.dxf": insunits_wrong_nomeasure_dxf,
                 "insunits_mm.dxf": insunits_mm_dxf, "mtext_codes.dxf": mtext_codes_dxf,
                 **{n: globals()[n.replace(".", "_")] for n in HOSTILE_FIXTURES},
                 **{n: globals()[n.replace(".", "_")] for n in REVIEW_18_FIXTURES if not n.startswith("ky4_cut_")},
                 **{n: globals()[n.replace(".", "_")] for n in FINAL_18_FIXTURES},
                 **{f"ky4_cut_{round(100 * s)}.dxf": functools.partial(ky4_cut_short_dxf, s) for s in CUT_SHORT_SHARES}})


def build(name: str, dest: str) -> str:
    path = os.path.join(dest, name)
    with open(path, "wb") as fh:
        fh.write(FIXTURES[name]())
    return path


def build_all(dest: str, names=None) -> dict[str, str]:
    os.makedirs(dest, exist_ok=True)
    return {n: build(n, dest) for n in (names or FIXTURES)}


# ----------------------------------------------------------------------------- the task-15 record, in parts
def _read(path: str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def sniff_part(paths: dict) -> dict:
    from .importers import sniff
    rows = []
    for name, kind, by_ext, backup in SNIFF_TABLE:
        s = sniff(_read(paths[name]), name)
        rows.append({"fixture": name, "expected": kind, "kind": s.kind, "category": s.category,
                     "by_extension": s.by_extension, "expected_by_extension": by_ext, "backup": s.backup,
                     "expected_backup": backup, "confidence": s.confidence, "detail": s.detail,
                     "ok": s.kind == kind and s.by_extension == by_ext and s.backup == backup})
    ties = []
    for name, kind in TIE_BREAKS:
        s = sniff(_read(paths[name]), name)
        ties.append({"fixture": name, "expected": kind, "kind": s.kind, "by_extension": s.by_extension,
                     "confidence": s.confidence, "says_so": "(from the file name)" in s.detail,
                     "ok": s.kind == kind and s.by_extension and s.confidence == "low" and "(from the file name)" in
                     s.detail})
    review = []
    for name, kind, bound in REVIEW_SNIFF:
        data = _read(paths[name])
        t0 = time.perf_counter()
        try:
            s, raised = sniff(data, name), None
        except Exception as e:  # noqa: BLE001  (sniff must not raise; a raise is recorded as a failure)
            s, raised = None, type(e).__name__
        fast = bound is None or time.perf_counter() - t0 < bound
        row = {"fixture": name, "bytes": len(data), "expected": kind, "kind": s.kind if s else None,
               "by_extension": s.by_extension if s else None, "confidence": s.confidence if s else None,
               "detail": s.detail if s else None, "raised": raised}
        if bound is not None:
            row[f"within_{bound:g}_s"] = fast
        row["ok"] = raised is None and (kind is None or s.kind == kind) and not s.by_extension and fast
        review.append(row)
    return {"fixtures": len(rows), "identified": sum(r["ok"] for r in rows),
            "unknown": sum(r["kind"] == "unknown" for r in rows),
            "decided_with_the_name": sum(r["by_extension"] for r in rows),
            "name_use_stated": sum(r["by_extension"] and "(from the file name)" in r["detail"] for r in rows),
            "backups": sum(r["backup"] for r in rows),
            "kinds": len({r["kind"] for r in rows}), "rows": rows,
            "tie_breaks": ties, "tie_breaks_ok": sum(t["ok"] for t in ties),
            "review_cases": len(review), "review_ok": sum(r["ok"] for r in review),
            "review_timed": sum(b is not None for _, _, b in REVIEW_SNIFF),
            "review_timed_within_bound": sum(any(k.startswith("within_") and v for k, v in r.items()) for r in review),
            "review_raised": sum(r["raised"] is not None for r in review),
            "xml_comment_run": XML_COMMENT_RUN, "review_rows": review}


def _same_network(path: str, ref) -> bool:
    import wntr
    wn = wntr.network.WaterNetworkModel(path)
    return (sorted(wn.node_name_list) == sorted(ref.node_name_list)
            and sorted(wn.link_name_list) == sorted(ref.link_name_list))


def encodings_part(paths: dict, out_dir: str) -> dict:
    import hashlib
    from .importers import prepare_inp
    import warnings
    import wntr
    ref = _wn("Net3")
    rows = []
    for name in ENCODING_CASES:
        try:                    # the control: WNTR on the file as it arrived
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wntr.network.WaterNetworkModel(paths[name])
            as_is = "loads"
        except Exception as e:  # noqa: BLE001
            as_is = f"raises {type(e).__name__}"
        p, rep = prepare_inp(paths[name], out_dir)
        rows.append({"fixture": name, "wntr_on_the_original": as_is, "prepared": rep.prepared,
                     "new_file": p is not None and p != paths[name],
                     "decoded_as": dict((i, v) for i, v, _ in rep.rows)["Text encoding"], "fixes": rep.fixes,
                     "wntr_loads_same_network": p is not None and _same_network(p, ref), "ok": rep.ok})
    def through(names):
        out = []
        for name in names:
            before = hashlib.sha1(_read(paths[name])).hexdigest()
            files_before = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
            p, rep = prepare_inp(paths[name], out_dir)
            out.append({"fixture": name, "same_path": p == paths[name], "same_file": rep.same_file,
                        "bytes_unchanged": hashlib.sha1(_read(paths[name])).hexdigest() == before,
                        "nothing_written": (sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []) == files_before,
                        "fixes": rep.fixes, "ok": rep.ok})
        return out

    passes = through(PASS_THROUGH_CASES)
    # the review's cases: line ends WNTR reads (CR only, CR-CR-LF) and files whose sections start after 64 KB pass
    # through untouched, and WNTR reads each as Net3
    review = through(REVIEW_PASS_THROUGH)
    for r in review:
        net = PASS_THROUGH_NETWORK.get(r["fixture"], "Net3")
        r["network"] = net
        r["wntr_loads_same_network"] = r["ok"] and _same_network(paths[r["fixture"]], ref if net == "Net3" else _wn(net))
    gaps = []
    for name in REVIEW_ENCODING:
        p, rep = prepare_inp(paths[name], out_dir)
        text = _read(p).decode("utf-8") if p else ""
        gaps.append({"fixture": name, "prepared": rep.prepared, "decoded_as": dict((i, v) for i, v, _ in rep.rows)
                     ["Text encoding"] if rep.ok else None, "fixes": rep.fixes,
                     "euro_sign_and_quotes_kept": "€" in text and "“revisado”" in text,
                     "wntr_loads_same_network": p is not None and _same_network(p, ref), "ok": rep.ok})
    # the final verification: a clean model uploaded (as bytes, the app's way) under a 250-character name was accepted
    # and then raised OSError at the rename, leaving its temporary copy; the name is shortened now.  A name that fits
    # is kept as it is, so the app's upload and grid cache names do not change.
    from .importers import MAX_NAME_BYTES
    net3 = _net3()
    names = []
    for label, up_name in (("long", LONG_UPLOAD_NAME), ("fits", "Z" * (MAX_NAME_BYTES - 13) + ".inp")):
        with tempfile.TemporaryDirectory(prefix="rm_names_") as nd:
            try:
                p, rep = prepare_inp(net3, nd, name=up_name)
                raised = None
            except Exception as e:  # noqa: BLE001  (recorded as a failure)
                p, rep, raised = None, None, type(e).__name__
            left = sorted(os.listdir(nd))
        names.append({"case": label, "name_characters": len(up_name), "raised": raised, "ok": bool(rep and rep.ok),
                      "prepared": rep.prepared if rep else None,
                      "prepared_bytes": len(rep.prepared.encode("utf-8")) if rep and rep.prepared else None,
                      "name_kept": bool(rep and rep.prepared == f"{hashlib.sha1(net3).hexdigest()[:8]}_{up_name}"),
                      "files_left": left})
    return {"normalized": sum(r["new_file"] and r["wntr_loads_same_network"] and bool(r["fixes"]) for r in rows),
            "of": len(rows), "passed_through": sum(r["same_path"] and r["bytes_unchanged"] and r["nothing_written"]
                                                   for r in passes),
            "pass_through_of": len(passes), "rows": rows, "pass_through": passes,
            "review_passed_through": sum(r["same_path"] and r["bytes_unchanged"] and r["nothing_written"]
                                         and r["wntr_loads_same_network"] for r in review),
            "review_pass_through_of": len(review), "review_pass_through": review,
            "review_cp1252_gaps": gaps, "review_upload_names": names}


def ids_part(paths: dict, out_dir: str) -> dict:
    """Each ID problem blocks before WNTR is called (WNTR's model class and its .inp reader are replaced by a trap for
    the call), and, as a control, what WNTR itself does with each file when nothing checks it first."""
    import warnings
    import wntr
    from . import importers
    calls = []
    keep_wn, keep_read = wntr.network.WaterNetworkModel, wntr.epanet.io.InpFile.read

    def trap(*a, **k):
        calls.append(a[:1])
        raise AssertionError("WNTR was called")

    rows, review = [], []
    try:
        wntr.network.WaterNetworkModel = trap
        wntr.epanet.io.InpFile.read = trap
        for names, into in ((DUPLICATE_CASES, rows), (REVIEW_ID_BLOCKS, review)):
            for name in names:
                p, rep = importers.prepare_inp(paths[name], out_dir)
                into.append({"fixture": name, "blocked": p is None and bool(rep.blocking), "blocking": rep.blocking,
                             "export_advice": all("IDs, not labels" in b for b in rep.blocking)})
    finally:
        wntr.network.WaterNetworkModel, wntr.epanet.io.InpFile.read = keep_wn, keep_read
    # the review's controls: the accented ID loads in WNTR from its UTF-8 copy (31 characters), and EPANET then
    # stops on it (33 bytes); the long pattern ID stops WNTR
    for r in review:
        text, _ = importers.decode_text(_read(paths[r["fixture"]]))
        with tempfile.TemporaryDirectory(prefix="rm_ids_") as tmp:
            utf8 = os.path.join(tmp, "copy.inp")
            with open(utf8, "w", encoding="utf-8") as fh:
                fh.write(text)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    wn = wntr.network.WaterNetworkModel(utf8)
                wn.options.time.duration = 3600
                try:
                    wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(tmp, "run"))
                    said = "loads, and the simulator runs"
                except Exception as e:  # noqa: BLE001
                    said = f"loads, and the simulator raises {type(e).__name__}"
            except Exception as e:  # noqa: BLE001
                said = f"raises {type(e).__name__}"
        r["without_the_scan"] = said
    # a model blocked after the ID scan (WNTR cannot read it) leaves no prepared copy behind
    files_before = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
    pb, repb = importers.prepare_inp(paths[REVIEW_WNTR_BLOCK], out_dir)
    files_after = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
    wntr_block = {"fixture": REVIEW_WNTR_BLOCK, "blocked": pb is None and bool(repb.blocking),
                  "blocked_by_wntr": any(b.startswith("EPANET (through WNTR)") for b in repb.blocking),
                  "fixes_it_needed": repb.fixes, "prepared": repb.prepared, "nothing_left": files_after == files_before,
                  "no_folder_in_message": not any(os.sep in b for b in repb.blocking)}
    # the final verification: a model WNTR reads and EPANET cannot run (a junction no pipe connects) was reported as
    # 'Reads and runs' and then stopped the app; EPANET's test run of the first hour blocks it now, given as a path
    # (nothing is written) and as bytes (the app's way: the temporary copy is removed)
    test_run = {"fixture": REVIEW_TEST_RUN_BLOCK}
    for how, src in (("path", paths[REVIEW_TEST_RUN_BLOCK]), ("bytes", _read(paths[REVIEW_TEST_RUN_BLOCK]))):
        before = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
        pt, rept = importers.prepare_inp(src, out_dir, name=REVIEW_TEST_RUN_BLOCK)
        after = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
        test_run[how] = {"blocked": pt is None and bool(rept.blocking), "blocking": rept.blocking,
                         "blocked_by_test_run": any(b.startswith("EPANET could not run the first hour") for b in
                                                    rept.blocking),
                         "names_the_lone_node": any("(9999)" in b for b in rept.blocking),
                         "wntr_read_it": any(i == "Junctions" for i, _, _ in rept.rows), "nothing_left": after == before}
    with tempfile.TemporaryDirectory(prefix="rm_ids_") as tmp:     # the control: WNTR alone loads it, EPANET stops
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wn = wntr.network.WaterNetworkModel(paths[REVIEW_TEST_RUN_BLOCK])
            wn.options.time.duration = 3600
            try:
                wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(tmp, "run"))
                test_run["without_the_test_run"] = "loads, and the simulator runs"
            except Exception as e:  # noqa: BLE001
                test_run["without_the_test_run"] = f"loads, and the simulator raises {type(e).__name__}"
        except Exception as e:  # noqa: BLE001
            test_run["without_the_test_run"] = f"raises {type(e).__name__}"
    n_lines = {"net3_dup_junction.inp": ("junctions", 93), "net3_junction_tank_collision.inp": ("junctions", 93),
               "net3_dup_pipe.inp": ("pipes", 118), "net3_long_id_32.inp": ("junctions", 92)}
    control = []
    for name in DUPLICATE_CASES:
        what, lines = n_lines[name]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wn = wntr.network.WaterNetworkModel(paths[name])
            got = wn.num_junctions if what == "junctions" else wn.num_pipes
            said = f"loads without a word: {got} {what} from {lines} lines"
            if name == "net3_junction_tank_collision.inp":      # it loads; the simulator is where it breaks
                wn.options.time.duration = 3600
                with tempfile.TemporaryDirectory(prefix="rm_ids_") as tmp:
                    try:
                        wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(tmp, "run"))
                        said += "; the simulator runs"
                    except Exception as e:  # noqa: BLE001
                        said += f"; the simulator raises {type(e).__name__}"
            control.append({"fixture": name, "wntr": said})
        except Exception as e:  # noqa: BLE001
            control.append({"fixture": name, "wntr": f"raises {type(e).__name__}"})
    ok31, rep31 = importers.prepare_inp(paths["net3_long_id_31.inp"], out_dir)
    return {"blocked": sum(r["blocked"] and r["export_advice"] for r in rows), "of": len(rows),
            "wntr_calls_while_blocking": len(calls), "rows": rows, "without_the_scan": control,
            "id_of_31_characters_runs": ok31 is not None and rep31.ok,
            "review_blocked": sum(r["blocked"] and r["export_advice"] for r in review), "review_of": len(review),
            "review_rows": review, "blocked_by_wntr_leaves_nothing": wntr_block,
            "blocked_by_the_test_run": test_run}


def zips_part(paths: dict, out_dir: str) -> dict:
    from .importers import ZipRefused, prepare_inp, sniff, zip_listing
    refused, opens = [], []
    keep_open = zipfile.ZipFile.open

    def counted_open(self, *a, **k):        # a member opened while a refused zip is handled would be counted here
        opens.append(a[:1])
        return keep_open(self, *a, **k)

    try:
        zipfile.ZipFile.open = counted_open
        for name in REFUSED_ZIPS:
            s = sniff(_read(paths[name]), name)
            p, rep = prepare_inp(paths[name], out_dir)
            # refused at the listing (an unsafe path), or when the one member to read is checked (the ratio)
            reason = s.detail if s.kind == "refused_zip" else " ".join(rep.blocking)
            refused.append({"fixture": name, "kind": s.kind, "reason": reason, "prepared": p is not None,
                            "refused_at": "listing" if s.kind == "refused_zip" else "reading the model",
                            "refused": p is None and ("refused" in reason or s.kind == "refused_zip")})
    finally:
        zipfile.ZipFile.open = keep_open
    two_p, two = prepare_inp(paths["zip_two_models.zip"], out_dir)
    em_p, em = prepare_inp(paths["emailed_model.zip"], out_dir)

    def dbf_stats(path):
        with zipfile.ZipFile(path) as z:
            i = [x for x in z.infolist() if x.filename.endswith(".dbf") and "__MACOSX" not in x.filename][0]
        return {"dbf_ratio": round(i.file_size / i.compress_size, 1), "dbf_mb": round(i.file_size / 1e6, 3)}

    pad_s = sniff(_read(paths["shp_padded_dbf.zip"]), "shp_padded_dbf.zip")
    try:
        zip_listing(_read(paths["emailed_model.zip"]), max_total=1000)
        cap_refuses = False
    except ZipRefused:
        cap_refuses = True
    # the review's cases: zips made on a Mac (the __MACOSX/ copies are left out), a city-size shapefile whose DBF is
    # over 10 MB and packs past 100 to 1 (listed, never read, so it passes), and the model sent with that layer
    review = []
    for name in REVIEW_ZIPS:
        s = sniff(_read(paths[name]), name)
        row = {"fixture": name, "kind": s.kind, "detail": s.detail, "members": s.members}
        if s.kind == "epanet_inp_zip":
            p, rep = prepare_inp(paths[name], out_dir)
            row.update(prepared=rep.prepared, fixes=rep.fixes, blocking=rep.blocking, ok=rep.ok)
        else:
            row["ok"] = s.kind == "shapefile_zip" and s.detail.startswith("1 layer:")
        if name in ("city_mains_shp.zip", "model_and_gis.zip"):
            row.update(dbf_stats(paths[name]), records=CITY_RECORDS)
        review.append(row)
    # the final verification: a zip whose member offsets are broken lists fine and then raised ValueError on reading;
    # it is refused now, with nothing written
    broken = []
    for name in REVIEW_BROKEN_ZIPS:
        s = sniff(_read(paths[name]), name)
        before = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
        try:
            p, rep = prepare_inp(paths[name], out_dir)
            raised = None
        except Exception as e:  # noqa: BLE001  (recorded as a failure)
            p, rep, raised = None, None, type(e).__name__
        after = sorted(os.listdir(out_dir)) if os.path.isdir(out_dir) else []
        try:                    # the control: zipfile itself on the model member, in memory as prepare_inp reads it
            with zipfile.ZipFile(io.BytesIO(_read(paths[name]))) as z, z.open(s.members[0]) as fh:
                fh.read()
            zipfile_says = "reads"
        except Exception as e:  # noqa: BLE001
            zipfile_says = f"raises {type(e).__name__}"
        broken.append({"fixture": name, "kind": s.kind, "raised": raised, "blocking": rep.blocking if rep else None,
                       "refused": p is None and raised is None and bool(rep.blocking)
                       and rep.blocking[0].startswith("The zip is refused"),
                       "nothing_written": after == before, "zipfile_alone": zipfile_says})
    return {"refused": sum(r["refused"] for r in refused), "of": len(refused), "rows": refused,
            "zip_member_opens_while_refusing": len(opens),
            "two_models_blocked": two_p is None and bool(two.blocking),
            "emailed_model": {"prepared": em.prepared, "fixes": em.fixes, "ok": em.ok},
            "padded_dbf": {"kind": pad_s.kind, **dbf_stats(paths["shp_padded_dbf.zip"])},
            "total_cap_refuses_when_lowered": cap_refuses,
            "review_ok": sum(r["ok"] for r in review), "review_of": len(review), "review_rows": review,
            "review_broken": broken}


def record() -> dict:
    with tempfile.TemporaryDirectory(prefix="rm_fixtures_") as d:
        paths = build_all(os.path.join(d, "fx"))
        out = {"about": "Task 15's importer tests on files built by residualmap/fixtures.py from WNTR's bundled public "
                        "networks (converter self-consistency on public models, not real utility files). Written by "
                        "python -m residualmap.fixtures --record; the task15 checks rerun each part and compare.",
               "sniff": sniff_part(paths), "encodings": encodings_part(paths, os.path.join(d, "enc")),
               "duplicate_ids": ids_part(paths, os.path.join(d, "ids")),
               "zips": zips_part(paths, os.path.join(d, "zip"))}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", metavar="DIR", help="write every fixture into DIR")
    ap.add_argument("--record", action="store_true", help=f"write {os.path.relpath(RECORD, REPO)}")
    a = ap.parse_args(argv)
    if a.build:
        for n, p in build_all(a.build).items():
            print(p)
    if a.record:
        rec = record()
        os.makedirs(os.path.dirname(RECORD), exist_ok=True)
        with open(RECORD, "w") as fh:
            json.dump(rec, fh, indent=1)
            fh.write("\n")
        s = rec["sniff"]
        print(f"sniff: {s['identified']} of {s['fixtures']} identified, {s['unknown']} unknown; encodings: "
              f"{rec['encodings']['normalized']} of {rec['encodings']['of']} normalized; IDs: "
              f"{rec['duplicate_ids']['blocked']} of {rec['duplicate_ids']['of']} blocked; zips: "
              f"{rec['zips']['refused']} of {rec['zips']['of']} refused -> {os.path.relpath(RECORD, REPO)}")
    if not (a.build or a.record):
        ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
