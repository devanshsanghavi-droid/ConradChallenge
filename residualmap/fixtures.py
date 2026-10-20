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
