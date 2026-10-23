"""
gis_read.py: GIS layers of water mains read into a skeleton model (iteration 5, task 17).

Replies in October 2026 from people who work with small systems said those systems rarely have a hydraulic model, only
a map of their mains, so this is the route for them: a mains layer becomes a skeleton EPANET model (pipes, their sizes
and lengths, nodes where they meet), with a readiness report that says it is not ready until the source, tanks, pumps,
valves, controls and demand are given (in EPANET, on the downloaded file, until the operator form, built in task 19,
is offered on the page; on the command line, --form applies it).
Converter self-consistency on public models converted to GIS by our own writers is all that is tested; no utility's
own layer was read.

    python -m residualmap.gis_read LAYER [--unit ft|us_ft|m|deg] [--tolerance-ft 6] [--out skeleton.inp]
    python -m residualmap.gis_read LAYER --form FORM.csv --out completed.inp   # task 19: completed by the operator
                                    # form (operator_form.py), each pipe's material and the layer's frame kept

Readers, all in memory (a zip is never unpacked to disk; residualmap.importers' zip limits apply):
- GeoPackage: sqlite3 on an in-memory copy (a step budget stops a crafted query), the GeoPackage geometry header and a
  WKB parser that reads Z, M and ZM (kept in 2D) and keeps a MultiLineString's parts apart (the d3 prototype joined
  them into phantom pipes; critique 12).  The unit comes from gpkg_spatial_ref_sys: the definition's last top-level
  UNIT of its PROJCS, or EPSG_UNITS; an undefined system is asked.
- Zipped shapefile: pyshp 3.1.6 on the members' bytes (never a path, so pyshp never opens a URL or a file).  The unit
  comes from the .prj's last top-level UNIT of the PROJCS (a State Plane .prj lists its GEOGCS's degrees first), the
  encoding from the .cpg, else latin-1 with a warning.  Null shapes and deleted DBF records are dropped with their
  attributes, so the rest stay aligned (records(deleted_as_None=True)).  Parts are kept apart.
- GeoJSON (and Esri JSON): a crs member is honoured; otherwise coordinates within longitude and latitude ranges are
  read as WGS84 (RFC 7946), anything else is asked.  Coordinates written with few decimals are detected (4 decimals of
  a degree is about 11 m) and warned.
- KML and KMZ: xml.etree and zipfile, WGS84 only; a document type declaration is refused before parsing (it can hide
  entity expansion), wherever it is.  Sizes come from the name, the description or the ExtendedData.
- Longitude and latitude go through the same WGS84 to east-north conversion as task 16b (importers.enu_from_latlon), on
  a plane at the layer's median position, so pyproj is not needed.

Then residualmap.topology builds the network (tees, snapping, wide nodes, IDs), and fields are mapped from a list of
guesses the operator confirms: the diameter (topology.parse_diameter, bare numbers read by their range), the material
(Hazen-Williams C from topology.C_TABLE, 130 the labelled default), the install year (shown), the declared length
(geometry length in true units is used; the field is a cross-check: a ratio near 3.28 or 0.3048 means wrong units),
the layer's own pipe ID.  A layer with no length field, read in a unit the operator picked, is asked whether it is to
scale, with an optional check against known lengths.  A valve layer's status field (NORMAL_POS and similar) gives
proposed closed-valve rows for the operator form, for confirmation.  Every junction is at one elevation (a skeleton
has no pressure-control valve; task 15b's rule applies once one is added).

A layer in WNTR's own GIS schema, with WNTR's JSON file beside it, is read through wntr.network.from_dict after
demand_pattern is copied into pattern_name (WNTR 1.5.0's from_dict reads pattern_name; without the copy every
junction falls back to the default pattern).  This is for WNTR users only: other tools do not write WNTR's schema.
It is tested on files our own writer lays out as WNTR's to_gis does (fixtures.wntr_schema_zip), never on files WNTR
itself wrote (to_gis needs geopandas, which is not installed).

The review of task 17 (each fix with a fixture in its check):
- Coordinates that are not numbers, or beyond COORD_MAX in any unit, leave their part out with a warning (_finish), and
  build_from_layer never raises: anything unforeseen is a plain refusal (one NaN vertex stopped the page).
- A layer whose unit is not known is asked with 'degrees (longitude and latitude)' among the choices, and told when its
  coordinates all lie in longitude and latitude ranges; a layer that says it is in longitude and latitude but whose
  coordinates mostly lie outside them is asked as mislabelled (a few such parts are left out with a warning).  A layer
  whose ends crowd into a few squares of the tolerance is refused before anything is joined (topology), and one that
  builds no pipe, or loses most of its pieces, is refused, never written as an empty model.
- Web Mercator (EPSG 3857 and its aliases, or a definition that names it) is converted to longitude and latitude by
  the spherical inverse that defines it, since its metres are not ground metres (1.31 at 40 degrees north); any other
  Mercator is refused with what to export instead.
- In a KML, a placemark's name or description is taken as its size only when most of them carry a unit (8 in, 8",
  200 mm, DN200), and an attribute table in the description (ArcGIS's Layer To KML) is read as fields.
- A layer over MAX_FEATURES that holds no lines is skipped with a note, not the whole file refused; a zip's members are
  read with one opening of it, at most MAX_ZIP_LAYERS shapefiles; two shapefiles of one name in different folders
  keep their folder in their layer name; Esri JSON over MAX_FEATURES is refused as GeoJSON is, a spatialReference
  given as its code alone is read, and features whose geometry cannot be read are counted and warned.

The final verification of task 17 (each case in gis_final_verification):
- Every text a reader scans is read in time linear in it: a KML's coordinates (_kml_points, with a budget of
  WKB_MAX_POINTS points a file), its tags (_strip_tags) and its descriptions' attribute tables (description_fields)
  in one pass each, and a coordinate system's definition (_top_level_items); patterns there had taken time
  quadratic or cubic in a run of spaces or of '<', so a KMZ of 5 KB within the hosted demo's 5 MB once unpacked ran
  for hours, holding the interpreter lock, and a .prj zipped in 738 bytes for 28 s.  A KML name or description is
  cut at TEXT_MAX_CHARS and a definition over WKT_MAX_CHARS is not read.
- A document type or entity declaration is refused wherever it is in a KML (an expat pass with only those handlers;
  the guard looked at the first 200,000 characters only).
- A layer over MAX_FEATURES is left unread with a note naming it, whatever its type (a large layer of service lines,
  or points declared GEOMETRY, refused the whole file), and one chosen by name is refused saying to export the mains
  alone.  A zip that holds one GeoPackage, GeoJSON or KML and no shapefile is read as that file.
- A Mercator named by its code alone (OTHER_MERCATOR_EPSG: World Mercator 3395 and others) is refused as one named
  by its definition is, not asked for its unit.  A layer in longitude and latitude whose two numbers look swapped
  (latitude first: its median latitude beyond SWAP_LAT, or most of it out of range, and every coordinate in range
  swapped) asks their order (GisOptions.swap_axes).  A layer refused as crowded, or in longitude and latitude ranges
  while its file names a projected system, sets GisBuild.offer_unit, so the page offers the unit over the file's.
- KML sizes: 'in' after a number is inches only when no word follows it ('Line 3 in Oak St').  Install years
  (year_of): epoch milliseconds are read as dates, and a year in a text only when it stands alone.
- A crs member whose properties are no object, Esri attributes or paths of the wrong type, and text that is not UTF-8
  in a GeoPackage are read in plain words, not raised.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
import os
import re
import sqlite3
import struct
import sys
import zipfile
from dataclasses import dataclass, field

import numpy as np

from . import topology as T
from .importers import (ASSUMED, DERIVED, FROM_FILE, JSON_MAX_BYTES, MISSING, ZIP_MAX_TOTAL, LogRefused, Readiness,
                        ZipRefused, decode_text, enu_from_latlon, plane_centre, read_zip_member, read_zip_members, sniff,
                        zip_junk, zip_listing)

GPKG_MAX_STEPS = 200_000              # progress calls (1,000 SQLite steps each) for one query on a GeoPackage layer
WKB_MAX_DEPTH = 8                     # nesting of geometry collections followed
WKB_MAX_POINTS = 5_000_000            # points read from one file
MAX_FEATURES = 200_000                # features read from one layer
QUANTIZED_DEGREES = 4                 # longitude and latitude written with this many decimals or fewer warn (11 m)
LENGTH_RATIO_TOL = 0.05               # declared over geometry length within 5% of 1, 3.281 or 0.3048
SCALE_TOL = 0.05                      # known lengths: the median ratio within 5% of 1 says the layer is to scale
DEG_TO_M = 111_320.0                  # metres per degree of latitude, for the quantization sentence only
GIS_KINDS = ("geopackage", "shapefile_zip", "geojson", "esri_json", "kml", "kmz")
# a few EPSG codes and their unit (an undefined or unlisted code is asked)
EPSG_UNITS = {4326: "deg", 4269: "deg", 4267: "deg", 4258: "deg", 4283: "deg", 3857: "m", 900913: "m", 102100: "m",
              102113: "m", 3785: "m", 27700: "m", 2193: "m", 3006: "m", 25832: "m", 25833: "m", 2227: "us_ft",
              2229: "us_ft", 2263: "us_ft", 2276: "us_ft", 2278: "us_ft"}
EPSG_UNITS.update({c: "m" for c in list(range(32601, 32661)) + list(range(32701, 32761)) + list(range(26901, 26924))
                   + list(range(26701, 26723))})
UNIT_FACTORS = ((T.FT_US, "us_ft"), (T.FT, "ft"), (1.0, "m"))
COORD_MAX = 1e9                       # a coordinate beyond this, in any unit, is no map position (the review of task 17)
MAX_ZIP_LAYERS = 100                  # shapefiles read from one zip
DROPPED_MAX_SHARE = 0.5               # a layer that would lose more than this share of its pieces is refused
WEB_MERCATOR_EPSG = {3857, 900913, 102100, 102113, 3785}
WEB_MERCATOR_R = 6378137.0            # the sphere of Web Mercator's definition (EPSG 3857), metres
WEB_MERCATOR_MAX = 20037508.342789244 # its x at 180 degrees, metres
UNIT_WORDS = {**T.UNIT_WORDS, "deg": "degrees (longitude and latitude)"}
# the final verification of task 17
WKT_MAX_CHARS = 65_536                # a coordinate system's definition (a real .prj is under 2 KB)
TEXT_MAX_CHARS = 65_536               # a KML name or description, as read (each is cut here before any pattern sees it)
# Mercator projections that are not Web Mercator, by their code alone (a definition that names Mercator is caught by
# its words): World Mercator, PDC Mercator, Mercator 41 (and their deprecated codes), Esri's World_Mercator, Brazil
# and Caspian Sea Mercator, and the old Simple Mercator
OTHER_MERCATOR_EPSG = {3395, 3832, 3994, 54004, 3349, 3752, 5641, 3388, 41001}
SWAP_LAT = 75.0                       # a layer whose median latitude is beyond this, read in longitude and latitude
                                      # that would fit swapped, is asked whether its axes are swapped
_WEB_MERC_WORDS = re.compile(r"mercator[ _]auxiliary[ _]sphere|pseudo[ _-]?mercator|popular[ _]visuali[sz]ation|"
                             r"web[ _-]?mercator|google[ _]maps[ _]global[ _]mercator", re.I)
_OTHER_MERC = re.compile(r'(?:PROJECTION|METHOD)\s*\[\s*"\s*mercator(?:[ _]?\(?(?:1sp|2sp|variant[ _][abc])\)?)?\s*"',
                         re.I)
_MERC_PARAM = re.compile(r'PARAMETER\s*\[\s*"\s*(false[ _]easting|false[ _]northing|central[ _]meridian|'
                         r'longitude[ _]of[ _](?:natural[ _])?origin)\s*"\s*,\s*([-+0-9.eE]+)', re.I)

_DIAMETER_EXACT = ("diameter", "diam", "pipediam", "pipe diam", "pipe diameter", "dia", "size", "pipesize", "pipe size",
                   "nom diam", "nominal diameter", "nomdiam", "nom dia", "diameter in", "diam in", "diameter mm",
                   "diam mm", "dia in", "dia mm", "main size", "mainsize", "pipe dia", "dn")
_MATERIAL_EXACT = ("material", "mat", "pipe mat", "pipemat", "pipe material", "pipematerial", "matl", "pipe matl",
                   "main material", "mainmat")
_YEAR_EXACT = ("installyear", "install year", "year inst", "yearinst", "yr inst", "install date", "installdate",
               "installed", "year installed", "yearinstalled", "inst year", "year", "date installed")
_LENGTH_EXACT = ("shape length", "shape_length", "shape leng", "shape len", "length", "len", "pipe length",
                 "pipelength", "length ft", "length m", "measuredlength", "measured length", "gis length", "leng")
_ID_EXACT = ("facilityid", "facility id", "pipeid", "pipe id", "mainid", "main id", "assetid", "asset id", "gisid",
             "gis id", "valveid", "valve id", "valve no", "hydrantid", "hydrant id", "id", "pipe no", "pipeno", "name")
_C_EXACT = ("c factor", "cfactor", "c_factor", "hw c", "hwc", "hazen c", "roughness", "c value", "cvalue")
_STATUS_EXACT = ("normal pos", "normalpos", "normal position", "norm pos", "normpos", "normally", "status",
                 "valve status", "valvestat", "position", "normal status", "operable position")
_CLOSED_WORDS = {"closed", "close", "c", "nc", "normally closed", "shut", "off", "cl", "clsd"}
_NOT_ID = {"objectid", "fid", "ogc fid", "ogc_fid", "oid", "shape", "geometry", "geom"}


class GisRefused(LogRefused):
    """A GIS file that cannot be read or built into a network, in plain words."""


@dataclass
class Feature:
    index: int
    kind: str                     # 'line', 'point', 'polygon', 'none'
    parts: list                   # each an (n, 2) array of the layer's coordinates (x, y or lon, lat)
    props: dict
    has_z: bool = False
    has_m: bool = False


@dataclass
class Layer:
    name: str
    source: str                   # 'geopackage', 'shapefile', 'geojson', 'esri_json', 'kml'
    features: list
    fields: list
    crs: dict                     # {'kind': 'projected' | 'geographic' | 'unknown', 'unit': 'ft' | 'us_ft' | 'm' |
                                  #  'deg' | None, 'unit_m', 'how', 'epsg'}
    notes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    encoding: str | None = None
    dropped: dict = field(default_factory=dict)   # {'null_shapes', 'deleted_records', ...}
    not_read: list = field(default_factory=list)  # the file's other layers left unread, as over MAX_FEATURES

    @property
    def n_lines(self) -> int:
        return sum(f.kind == "line" for f in self.features)

    @property
    def n_points(self) -> int:
        return sum(f.kind == "point" for f in self.features)


# ----------------------------------------------------------------------------- units
def _unit_of_factor(f: float) -> str | None:
    for v, k in UNIT_FACTORS:
        if abs(f - v) <= 1e-9 * max(1.0, v) or abs(f / v - 1) < 1e-7:
            return k
    return None


_WKT_KEYWORD = re.compile(r"[A-Za-z_]+")


def _top_level_items(body: str) -> list[tuple[str, str]]:
    """(keyword, its bracketed text) for each item directly inside one WKT item's brackets, in one pass (the final
    verification of task 17: a pattern with two optional runs of spaces, tried at every offset, took time cubic in a
    run of spaces and commas, 28 s for a .prj zipped in 738 bytes).  Spaces, commas and a quoted string between items
    are passed over; a word not followed by a bracket is passed over whole."""
    out, i, n = [], 0, len(body)
    while i < n:
        ch = body[i]
        if ch.isspace() or ch == ",":
            i += 1
            continue
        if ch == '"':
            q = body.find('"', i + 1)
            i = n if q < 0 else q + 1
            continue
        m = _WKT_KEYWORD.match(body, i)
        if not m:
            i += 1
            continue
        k = m.end()
        while k < n and body[k].isspace():
            k += 1
        if k >= n or body[k] not in "[(":
            i = k
            continue
        kw, start = m.group(0), k + 1
        depth, j, in_q = 1, start, False
        while j < n and depth:
            ch = body[j]
            if ch == '"':
                in_q = not in_q
            elif not in_q and ch in "[(":
                depth += 1
            elif not in_q and ch in "])":
                depth -= 1
            j += 1
        out.append((kw.upper(), body[start:j - 1]))
        i = j
    return out


def wkt_unit(wkt: str) -> dict:
    """{'kind', 'unit', 'unit_m', 'name'} of a .prj or GeoPackage definition: a PROJCS (or PROJCRS) takes its last
    top-level UNIT (a State Plane .prj names its GEOGCS's degrees first); a GEOGCS alone is degrees.  A definition
    longer than WKT_MAX_CHARS is no coordinate system this tool reads (the final verification of task 17)."""
    s = (wkt or "").strip()
    if len(s) > WKT_MAX_CHARS:
        return {"kind": "unknown", "unit": None, "unit_m": None, "name": None}
    m = re.match(r"([A-Za-z_]+)\s*[\[(]", s)
    if not m:
        return {"kind": "unknown", "unit": None, "unit_m": None, "name": None}
    kw = m.group(1).upper()
    body = s[m.end():s.rfind("]") if "]" in s else len(s)]
    if kw in ("GEOGCS", "GEOGCRS", "GEODCRS", "GEOGRAPHICCRS"):
        return {"kind": "geographic", "unit": "deg", "unit_m": None, "name": "degree"}
    if kw not in ("PROJCS", "PROJCRS", "PROJECTEDCRS"):
        return {"kind": "unknown", "unit": None, "unit_m": None, "name": None}
    units = []
    for k, text in _top_level_items(body):
        if k in ("UNIT", "LENGTHUNIT"):
            parts = [p.strip().strip('"') for p in text.split(",")]
            try:
                units.append((parts[0], float(parts[1])))
            except (IndexError, ValueError):
                continue
        elif k == "CS":
            continue
        elif k in ("COORDINATESYSTEM",):
            continue
    if not units:             # WKT2 puts the unit inside CS axes: take the last LENGTHUNIT anywhere
        for nm, f in re.findall(r'LENGTHUNIT\s*\[\s*"([^"]*)"\s*,\s*([0-9.eE+-]+)', s):
            units.append((nm, float(f)))
    merc = mercator_of_wkt(s)
    if not units:
        return {"kind": "projected", "unit": None, "unit_m": None, "name": None, **merc}
    name, f = units[-1]
    return {"kind": "projected", "unit": _unit_of_factor(f) or "other", "unit_m": f, "name": name, **merc}


def mercator_of_wkt(wkt: str) -> dict:
    """{'mercator': 'web' (Web Mercator, EPSG 3857 by its names), 'other' (any other Mercator, not transverse or
    oblique) or None, 'merc_params': {'fe', 'fn', 'lon0'}} of a projected definition (the review of task 17: Web
    Mercator's metres are not ground metres)."""
    s = wkt or ""
    kind = "web" if _WEB_MERC_WORDS.search(s) else "other" if _OTHER_MERC.search(s) else None
    params = {"fe": 0.0, "fn": 0.0, "lon0": 0.0}
    for nm, v in _MERC_PARAM.findall(s):
        n = nm.lower().replace(" ", "_")
        key = "fe" if n.startswith("false_e") else "fn" if n.startswith("false_n") else "lon0"
        with contextlib.suppress(ValueError):
            params[key] = float(v)
    return {"mercator": kind, "merc_params": params} if kind else {"mercator": None}


def _crs(kind: str, unit: str | None, how: str, unit_m: float | None = None, epsg: int | None = None) -> dict:
    if unit in T.UNIT_M and unit_m is None:
        unit_m = T.UNIT_M[unit]
    return {"kind": kind, "unit": unit, "unit_m": unit_m, "how": how, "epsg": epsg}


def _with_merc(crs: dict, w: dict) -> dict:
    """crs with wkt_unit's Mercator keys (w) copied, when it names a Mercator."""
    if w.get("mercator"):
        crs = {**crs, "mercator": w["mercator"], "merc_params": w.get("merc_params") or {}}
    return crs


def _epsg_crs(code: int, how: str) -> dict:
    if int(code) in OTHER_MERCATOR_EPSG:   # the final verification of task 17: its unit was asked, then read as ground
        return {**_crs("projected", "m", how + f" (EPSG {code}, a Mercator)", epsg=int(code)),   # metres
                "mercator": "other", "merc_params": {}}
    u = EPSG_UNITS.get(int(code))
    if u == "deg":
        return _crs("geographic", "deg", how, epsg=int(code))
    if u:
        crs = _crs("projected", u, how, epsg=int(code))
        if int(code) in WEB_MERCATOR_EPSG:
            crs.update(mercator="web", merc_params={"fe": 0.0, "fn": 0.0, "lon0": 0.0})
        return crs
    return _crs("unknown", None, how + f" (EPSG {code}, a code whose unit is not in this tool's short list)",
                epsg=int(code))


# ----------------------------------------------------------------------------- WKB
class _WkbError(ValueError):
    pass


def parse_wkb(b: bytes, off: int = 0, depth: int = 0, budget: list | None = None):
    """(geometries, offset after) of a WKB geometry (ISO, with Z, M and ZM types, or EWKB flags): each geometry as
    (kind, (n, 2) array), a MultiLineString's parts as separate lines, Z and M read and dropped."""
    budget = budget if budget is not None else [WKB_MAX_POINTS]
    if depth > WKB_MAX_DEPTH:
        raise _WkbError("geometry collections nested too deep")
    if off + 5 > len(b):
        raise _WkbError("truncated geometry")
    order = "<" if b[off] == 1 else ">" if b[off] == 0 else None
    if order is None:
        raise _WkbError("not a WKB geometry")
    (t,) = struct.unpack_from(order + "I", b, off + 1)
    off += 5
    has_z = has_m = False
    if t & 0x20000000:                    # EWKB SRID
        off += 4
    if t & 0x80000000:
        has_z = True
    if t & 0x40000000:
        has_m = True
    t &= 0x0FFFFFFF
    base, dim = t % 1000, t // 1000
    has_z = has_z or dim in (1, 3)
    has_m = has_m or dim in (2, 3)
    nd = 2 + has_z + has_m

    def pts(o):
        if o + 4 > len(b):
            raise _WkbError("truncated geometry")
        (n,) = struct.unpack_from(order + "I", b, o)
        o += 4
        need = n * nd * 8
        if n > budget[0] or o + need > len(b):
            raise _WkbError("a geometry with more points than the file holds")
        budget[0] -= n
        a = np.frombuffer(b, dtype=(order + "f8"), count=n * nd, offset=o).reshape(n, nd)[:, :2].astype(float)
        return a, o + need

    flags = {"z": has_z, "m": has_m}
    if base == 1:
        if off + nd * 8 > len(b):
            raise _WkbError("truncated geometry")
        a = np.frombuffer(b, dtype=(order + "f8"), count=nd, offset=off)[:2].astype(float).reshape(1, 2)
        off += nd * 8
        out = [] if np.isnan(a).any() else [("point", a)]
        return out, off, flags
    if base == 2:
        a, off = pts(off)
        return [("line", a)], off, flags
    if base == 3:
        (nr,) = struct.unpack_from(order + "I", b, off)
        off += 4
        rings = []
        for _ in range(min(nr, 100000)):
            a, off = pts(off)
            rings.append(("polygon", a))
        return rings[:1], off, flags
    if base in (4, 5, 6, 7):
        if off + 4 > len(b):
            raise _WkbError("truncated geometry")
        (n,) = struct.unpack_from(order + "I", b, off)
        off += 4
        if n > budget[0] + 1:
            raise _WkbError("a collection with more parts than the file holds")
        out = []
        for _ in range(n):
            g, off, f2 = parse_wkb(b, off, depth + 1, budget)
            out += g
            flags["z"] |= f2["z"]
            flags["m"] |= f2["m"]
        return out, off, flags
    raise _WkbError(f"geometry type {t} is not read")


def gpkg_geometry(blob: bytes):
    """(geometries, srs_id, flags) of a GeoPackage geometry blob ('GP' header, then WKB), or of plain WKB."""
    if blob is None:
        return [], None, {"z": False, "m": False}
    b = bytes(blob)
    if b[:2] == b"GP":
        if len(b) < 8:
            raise _WkbError("truncated GeoPackage geometry")
        flags = b[3]
        if (flags >> 5) & 1:
            raise _WkbError("an extended GeoPackage geometry type is not read")
        order = "<" if flags & 1 else ">"
        (srs,) = struct.unpack_from(order + "i", b, 4)
        env = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}.get((flags >> 1) & 7)
        if env is None:
            raise _WkbError("a GeoPackage geometry header with an unknown envelope")
        if (flags >> 4) & 1:
            return [], srs, {"z": False, "m": False}
        g, _, f = parse_wkb(b, 8 + env)
        return g, srs, f
    g, _, f = parse_wkb(b, 0)
    return g, None, f


# ----------------------------------------------------------------------------- GeoPackage
def _gpkg_open(data: bytes) -> sqlite3.Connection:
    con = sqlite3.connect(":memory:")
    calls = [0]

    def budget() -> int:
        calls[0] += 1
        return 1 if calls[0] > GPKG_MAX_STEPS else 0

    try:
        con.deserialize(data)
        con.set_progress_handler(budget, 1000)
        # text that is not valid UTF-8 is read with the bad bytes replaced (the final verification of task 17: it
        # raised UnicodeDecodeError, which the page could only name)
        con.text_factory = lambda b: b.decode("utf-8", errors="replace")
    except BaseException:
        con.close()
        raise
    return con


def _q(name: str) -> str:
    return '"' + str(name).replace('"', '""') + '"'


def read_geopackage(data: bytes, name: str = "layer.gpkg") -> list[Layer]:
    try:
        con = _gpkg_open(data)
    except (sqlite3.Error, ValueError, OverflowError, MemoryError) as e:
        raise GisRefused(f"The GeoPackage cannot be opened ({type(e).__name__}). Save it again from your GIS.") from None
    layers, skipped = [], []
    with contextlib.closing(con):
        try:
            master = {n: (t, (q or "")) for n, t, q in con.execute(
                "select name, type, sql from sqlite_master where type in ('table','view') limit 10000")}
            if master.get("gpkg_contents", ("",))[0] != "table" or master.get("gpkg_geometry_columns", ("",))[0] != "table":
                raise GisRefused("The GeoPackage has no gpkg_contents or gpkg_geometry_columns table, so its layers "
                                 "cannot be found. Save the layer again as a GeoPackage from your GIS.")
            contents = con.execute("select table_name, data_type, srs_id from gpkg_contents limit 1000").fetchall()
            geom_cols = {t: (c, s, z, m, str(g or "")) for t, c, g, s, z, m in con.execute(
                "select table_name, column_name, geometry_type_name, srs_id, z, m from gpkg_geometry_columns "
                "limit 1000")}
            srs = {}
            if master.get("gpkg_spatial_ref_sys", ("",))[0] == "table":
                for sid, org, oid, defn in con.execute("select srs_id, organization, organization_coordsys_id, "
                                                       "definition from gpkg_spatial_ref_sys limit 10000"):
                    srs[sid] = (str(org or ""), oid, str(defn or ""))
            for table, dtype, c_srs in contents:
                if str(dtype).lower() != "features" or table not in geom_cols:
                    continue
                typ, sql = master.get(table, ("", ""))
                if typ != "table" or sql.lstrip().upper().startswith("CREATE VIRTUAL"):
                    continue
                gcol, s_id, _, _, gtype = geom_cols[table]
                (n_rows,) = con.execute(f"select count(*) from {_q(table)}").fetchone()
                if n_rows > MAX_FEATURES:       # the review of task 17, and its final verification: a large layer
                    skipped.append((str(table), f"{table} ({n_rows:,} features, {gtype or 'no type'})"))   # of any
                    continue                    # type beside the mains (service lines, meters) refused the file
                cur = con.execute(f"select * from {_q(table)} limit {MAX_FEATURES + 1}")
                cols = [d[0] for d in cur.description]
                rows = cur.fetchall()
                lay = Layer(name=str(table), source="geopackage", features=[], fields=[c for c in cols if c != gcol],
                            crs=_gpkg_crs(srs, s_id if s_id is not None else c_srs))
                if len(rows) > MAX_FEATURES:
                    raise GisRefused(f"The layer {table} has more than {MAX_FEATURES:,} features.")
                gi = cols.index(gcol)
                bad, budget = 0, [WKB_MAX_POINTS]
                for k, r in enumerate(rows):
                    props = {c: v for c, v in zip(cols, r) if c != gcol and not isinstance(v, (bytes, memoryview))}
                    try:
                        geoms, _, fl = gpkg_geometry(r[gi]) if r[gi] is not None else ([], None, {"z": 0, "m": 0})
                    except (_WkbError, struct.error, ValueError):
                        bad += 1
                        geoms, fl = [], {"z": False, "m": False}
                    lay.features.append(_feature(k, geoms, props, fl))
                if bad:
                    lay.warnings.append(f"{bad} feature{'s' if bad != 1 else ''} of {table} could not be read (a broken "
                                        f"or unsupported geometry) and {'are' if bad != 1 else 'is'} left out.")
                layers.append(_finish(lay))
        except sqlite3.OperationalError as e:
            raise GisRefused(f"Reading the GeoPackage was stopped ({e}): a query ran longer than this tool allows. "
                             f"Save the layer again as a GeoPackage from your GIS.") from None
        except sqlite3.Error as e:
            raise GisRefused(f"The GeoPackage cannot be read ({type(e).__name__}: {e}). Save it again from your GIS.") \
                from None
    if not layers and skipped:
        raise GisRefused(_over_cap_refusal("The GeoPackage", skipped))
    if not layers:
        raise GisRefused("The GeoPackage holds no feature layer. Save the water mains layer into it from your GIS.")
    _note_skipped(layers, "GeoPackage's", skipped)
    return layers


def _over_cap_refusal(what: str, skipped: list) -> str:
    """A file whose every layer is over MAX_FEATURES, refused naming them (the final verification of task 17)."""
    return (f"{what}'s layer{'s are' if len(skipped) != 1 else ' is'} each larger than {MAX_FEATURES:,} features, "
            f"the most this tool reads from one layer: {', '.join(t for _, t in skipped[:4])}. Export the water mains "
            f"layer alone (as a GeoPackage or a zipped shapefile), and use that.")


def _note_skipped(layers: list, whose: str, skipped: list) -> None:
    """Every layer read notes the layers left unread as over MAX_FEATURES, and keeps their names (Layer.not_read), so
    a layer chosen by name among them is refused in plain words (the final verification of task 17)."""
    if not skipped:
        return
    for lay in layers:
        lay.notes.append(f"Not read, as larger than {MAX_FEATURES:,} features: the {whose} "
                         f"layer{'s' if len(skipped) != 1 else ''} {', '.join(t for _, t in skipped[:4])}. If the "
                         f"water mains are one of them, export the mains layer alone.")
        lay.not_read = [n for n, _ in skipped]


def _gpkg_crs(srs: dict, sid) -> dict:
    if sid is None or sid not in srs:
        return _crs("unknown", None, "the GeoPackage names no coordinate system for the layer")
    org, oid, defn = srs[sid]
    if sid in (-1, 0) or defn.strip().lower() in ("undefined", ""):
        if sid == 0:
            return _crs("geographic", "deg", "the GeoPackage's undefined geographic system (srs_id 0)")
        return _crs("unknown", None, "the GeoPackage's coordinate system is 'undefined'")
    code = None
    if str(org).upper() in ("EPSG", "ESRI") and oid is not None:
        with contextlib.suppress(TypeError, ValueError):
            code = int(oid)
    if code in OTHER_MERCATOR_EPSG:           # the final verification of task 17
        return _epsg_crs(code, "the GeoPackage's EPSG code")
    w = wkt_unit(defn)
    if w["kind"] != "unknown" and w["unit"] not in (None, "other"):
        code = code if str(org).upper() == "EPSG" else None
        crs = _with_merc(_crs(w["kind"], w["unit"], f"the GeoPackage's coordinate system ({w['name'] or 'degree'})",
                              w["unit_m"], code), w)
        if code in WEB_MERCATOR_EPSG and not crs.get("mercator"):
            crs = _with_merc(crs, {"mercator": "web", "merc_params": {}})
        return crs
    if str(org).upper() == "EPSG" and code is not None:
        return _epsg_crs(code, "the GeoPackage's EPSG code")
    return _crs("unknown", None, "the GeoPackage's coordinate system definition names no unit this tool reads")


def _feature(k: int, geoms: list, props: dict, fl: dict) -> Feature:
    kinds = {g[0] for g in geoms}
    if not geoms:
        kind = "none"
    elif kinds == {"line"}:
        kind = "line"
    elif kinds == {"point"}:
        kind = "point"
    elif "line" in kinds:
        kind = "line"
    else:
        kind = sorted(kinds)[0]
    parts = [g[1] for g in geoms if g[0] == kind]
    return Feature(index=k, kind=kind, parts=parts, props=props, has_z=bool(fl.get("z")), has_m=bool(fl.get("m")))


def _keep_parts(lay: Layer, ok) -> int:
    """Leave out every part of a line or point feature for which ok(part) is false; a feature left with no part has no
    geometry ('none').  Returns how many parts were left out."""
    n = 0
    for f in lay.features:
        if f.kind not in ("line", "point") or not f.parts:
            continue
        keep = [p for p in f.parts if ok(p)]
        if len(keep) != len(f.parts):
            n += len(f.parts) - len(keep)
            f.parts = keep
            if not keep:
                f.kind = "none"
    return n


def _in_lonlat(p: np.ndarray) -> bool:
    return bool(len(p)) and bool(np.all(np.abs(p[:, 0]) <= 180.0)) and bool(np.all(np.abs(p[:, 1]) <= 90.0))


def _line_points(lay: Layer) -> np.ndarray:
    pts = [part for f in lay.features if f.kind == "line" for part in f.parts if len(part)]
    return np.vstack(pts) if pts else np.zeros((0, 2))


def looks_lonlat(lay: Layer) -> bool:
    """Every coordinate of the layer's lines lies within longitude and latitude ranges (the review of task 17: the
    GeoJSON reader's test, applied to a shapefile or GeoPackage whose coordinate system is not known)."""
    a = _line_points(lay)
    return bool(len(a)) and _in_lonlat(a) and bool(np.ptp(a, axis=0).max() > 0)


def _finish(lay: Layer) -> Layer:
    """Every reader's last step (the review of task 17): parts with a coordinate that is not a number, or beyond
    COORD_MAX, or a line part of fewer than two points, are left out with a warning (one NaN vertex stopped cKDTree,
    and the page with it); Web Mercator is converted to longitude and latitude, any other Mercator is refused; a layer
    said to be in longitude and latitude whose coordinates mostly lie outside them is marked mislabelled (its unit is
    asked; a KML is refused), and the few parts outside them of one that mostly lies inside are left out."""
    bad = 0
    for f in lay.features:
        if f.kind not in ("line", "point") or not f.parts:
            continue
        need = 2 if f.kind == "line" else 1
        keep = [p for p in f.parts if np.ndim(p) == 2 and len(p) >= need and bool(np.isfinite(p).all())
                and float(np.abs(p).max()) <= COORD_MAX]
        if len(keep) != len(f.parts):
            bad += len(f.parts) - len(keep)
            f.parts = keep
            if not keep:
                f.kind = "none"
    if bad:
        lay.warnings.append(f"{bad} part{'s' if bad != 1 else ''} of the layer's features "
                            f"{'have' if bad != 1 else 'has'} no usable coordinates (a value that is not a number or "
                            f"beyond any map's range, or a line of fewer than two points) and "
                            f"{'are' if bad != 1 else 'is'} left out.")
    crs = lay.crs or {}
    merc = crs.get("mercator")
    if merc == "other":
        lay.crs = {**crs, "refusal": "The layer is in a Mercator projection, whose metres are not ground metres (they "
                                     "stretch with latitude), so its lengths cannot be read from it. Export it in your "
                                     "State Plane, UTM or other local projected system, or as longitude and latitude "
                                     "(WGS84)."}
    elif merc == "web":
        p = crs.get("merc_params") or {}
        fe, fn, lon0 = float(p.get("fe", 0.0)), float(p.get("fn", 0.0)), float(p.get("lon0", 0.0))
        inside = lambda a: bool(np.all(np.abs(a[:, 0] - fe) <= WEB_MERCATOR_MAX * 1.0001)
                                and np.all(np.abs(a[:, 1] - fn) <= 3.0 * WEB_MERCATOR_MAX))
        parts = [a for f in lay.features if f.kind in ("line", "point") for a in f.parts]
        n_out = sum(not inside(a) for a in parts)
        if parts and n_out > len(parts) / 2:
            lay.crs = _crs("unknown", None, f"{crs.get('how')} says Web Mercator, but its coordinates lie outside it")
            lay.crs["mislabelled"] = True
            return lay
        if n_out:
            _keep_parts(lay, inside)
            lay.warnings.append(f"{n_out} part{'s' if n_out != 1 else ''} of the layer's features lie outside Web "
                                f"Mercator's range and {'are' if n_out != 1 else 'is'} left out.")
        for f in lay.features:
            if f.kind in ("line", "point"):
                f.parts = [np.column_stack([np.degrees((a[:, 0] - fe) / WEB_MERCATOR_R) + lon0,
                                            np.degrees(2.0 * np.arctan(np.exp((a[:, 1] - fn) / WEB_MERCATOR_R))
                                                       - math.pi / 2)]) for a in f.parts]
        lay.crs = _crs("geographic", "deg", f"{crs.get('how')}: Web Mercator, converted to longitude and latitude, "
                                            f"since its metres are not ground metres", epsg=crs.get("epsg"))
        lay.crs["from_web_mercator"] = True
        lay.notes.append("The layer is in Web Mercator (EPSG 3857), whose metres are not ground metres: it was converted "
                         "to longitude and latitude, and its lengths are measured on the ground.")
        return lay
    if crs.get("kind") == "geographic":
        _check_lonlat(lay, crs)
    return lay


def _check_lonlat(lay: Layer, crs: dict, axes: bool = True) -> None:
    """The range test of a layer in longitude and latitude (part of _finish).  First (axes; the final verification of
    task 17), the axis order: a layer whose coordinates would all fit as longitude and latitude with their two numbers
    swapped, and whose median latitude is beyond SWAP_LAT (where almost no water network is) or most of whose parts
    lie outside the ranges, is marked axes_suspect, and its order is asked (latitude written first is a known
    EPSG:4326 mix-up, and a layer built near the pole has wrong lengths and no warning)."""
    parts = [a for f in lay.features if f.kind in ("line", "point") for a in f.parts]
    n_out = sum(not _in_lonlat(a) for a in parts)
    if axes and parts:
        a = np.vstack(parts)
        fits_swapped = bool(np.all(np.abs(a[:, 0]) <= 90.0) and np.all(np.abs(a[:, 1]) <= 180.0))
        med_lat = float(np.median(np.abs(a[:, 1])))
        if fits_swapped and (med_lat > SWAP_LAT or n_out > len(parts) / 2):
            lay.crs = {**crs, "axes_suspect": True, "axes_out_of_range": n_out > len(parts) / 2,
                       "median_lat": float(np.median(a[:, 1])), "median_lon": float(np.median(a[:, 0]))}
            return
    if parts and n_out > len(parts) / 2:
        big = float(np.abs(np.vstack(parts)).max())
        if lay.source == "kml":
            lay.crs = {**crs, "refusal": f"The KML's coordinates run to {big:,.0f}, outside longitude and latitude, "
                                         f"which is all KML holds: it cannot be read. Save it again from Google Earth."}
            return
        lay.crs = _crs("unknown", None, f"{crs.get('how')} says longitude and latitude, but its coordinates run to "
                                        f"{big:,.0f}, outside -180 to 180 and -90 to 90: the coordinate system is "
                                        f"probably mislabelled")
        lay.crs["mislabelled"] = True
    elif n_out:
        _keep_parts(lay, _in_lonlat)
        lay.warnings.append(f"{n_out} part{'s' if n_out != 1 else ''} of the layer's features "
                            f"{'have' if n_out != 1 else 'has'} coordinates outside longitude and latitude (beyond "
                            f"180 or 90 degrees) and {'are' if n_out != 1 else 'is'} left out.")


def swap_axes(lay: Layer) -> Layer:
    """A copy of the layer with each coordinate's two numbers swapped (latitude was written first), range-tested again
    (the final verification of task 17).  The layer itself is not changed."""
    feats = [Feature(index=f.index, kind=f.kind, parts=[np.ascontiguousarray(p[:, ::-1]) for p in f.parts],
                     props=f.props, has_z=f.has_z, has_m=f.has_m) for f in lay.features]
    crs = {k: v for k, v in lay.crs.items() if k not in ("axes_suspect", "axes_out_of_range", "median_lat",
                                                          "median_lon")}
    crs["how"] = f"{crs.get('how')}, with each coordinate's two numbers swapped by you (latitude was written first)"
    crs["axes_swapped"] = True
    out = Layer(name=lay.name, source=lay.source, features=feats, fields=list(lay.fields), crs=crs,
                notes=list(lay.notes), warnings=list(lay.warnings), encoding=lay.encoding, dropped=dict(lay.dropped),
                not_read=list(lay.not_read))
    _check_lonlat(out, crs, axes=False)
    return out


# ----------------------------------------------------------------------------- shapefiles
_CPG = {"utf8": "utf-8", "utf-8": "utf-8", "65001": "utf-8", "1252": "cp1252", "ansi 1252": "cp1252",
        "windows-1252": "cp1252", "88591": "latin-1", "8859-1": "latin-1", "iso-8859-1": "latin-1", "iso88591": "latin-1",
        "latin1": "latin-1", "latin-1": "latin-1", "ascii": "ascii", "us-ascii": "ascii", "936": "gbk", "932": "cp932",
        "sjis": "cp932", "1250": "cp1250", "1251": "cp1251"}


def _encoding_of_cpg(text: str) -> str | None:
    import codecs
    t = text.strip().lower()
    enc = _CPG.get(t)
    if enc is None:
        with contextlib.suppress(LookupError):
            enc = codecs.lookup(t).name
    return enc


def read_shapefile_zip(data: bytes, name: str = "layer.zip", max_unpacked: int = ZIP_MAX_TOTAL) -> list[Layer]:
    import shapefile
    try:
        infos = zip_listing(data, max_unpacked)
    except ZipRefused as e:
        raise GisRefused(f"The zip is refused: {e}. Nothing in it was opened.") from None
    files = {i.filename: i for i in infos if not i.is_dir() and not zip_junk(i.filename)}
    low = {n.lower(): n for n in files}
    if any(n.endswith((".gdbtable", ".gdbtablx")) or "/gdb" in n or n.endswith(".gdb/") for n in low):
        raise GisRefused("This zip holds an Esri file geodatabase (.gdb), which cannot be read here. In ArcGIS Pro, "
                         "export the water mains layer as a shapefile (zip its .shp, .shx, .dbf and .prj together) or "
                         "a GeoPackage, and use that.")
    shps = sorted(n for n in low if n.endswith(".shp"))
    if not shps:
        raise GisRefused("The zip holds no shapefile (.shp). Zip the layer's .shp, .shx, .dbf and .prj together.")
    if len(shps) > MAX_ZIP_LAYERS:      # the review of task 17: every layer's members cost a listing of the whole zip
        raise GisRefused(f"The zip holds {len(shps):,} shapefiles, and at most {MAX_ZIP_LAYERS} are read. Zip the "
                         f"water mains layer alone (its .shp, .shx, .dbf, .prj and .cpg).")
    want = [low[n[:-4] + e] for n in shps for e in (".shp", ".shx", ".dbf", ".prj", ".cpg") if n[:-4] + e in low]
    try:              # one listing and one opening of the zip for every member (read_zip_member's rules for each)
        got = read_zip_members(data, want, limit=max_unpacked, max_total=max_unpacked)
    except ZipRefused as e:
        raise GisRefused(f"The zip is refused: {e}.") from None
    # two shapefiles of one name in different folders keep their folder in their layer name (the review of task 17)
    stems = [os.path.basename(low[n])[:-4] for n in shps]
    names = [os.path.basename(low[n])[:-4] if stems.count(os.path.basename(low[n])[:-4]) == 1
             else low[n][:-4].strip("/") for n in shps]
    kml = sorted(low[n] for n in low if n.endswith(".kml"))
    layers, skipped = [], []
    for shp_low, lname in zip(shps, names):
        stem_low = shp_low[:-4]
        member = lambda ext: got.get(low.get(stem_low + ext))
        shp, shx, dbf, prj, cpg = (member(e) for e in (".shp", ".shx", ".dbf", ".prj", ".cpg"))
        warns, notes = [], []
        if cpg is not None:
            enc = _encoding_of_cpg(decode_text(cpg)[0])
            if enc is None:
                enc = "latin-1"
                warns.append(f"The layer's .cpg names an encoding this tool does not know "
                             f"({decode_text(cpg)[0].strip()[:20]!r}): its text is read as latin-1. Check accented names.")
            else:
                notes.append(f"Text read as {enc}, as the layer's .cpg says.")
        else:
            enc = "latin-1"
            if dbf is not None:
                warns.append("The layer has no .cpg file naming its text encoding: its text is read as latin-1. Check "
                             "accented names, or export the layer with its .cpg.")
        if prj is not None:
            w = wkt_unit(decode_text(prj)[0])
            if w["kind"] == "geographic":
                crs = _crs("geographic", "deg", "the layer's .prj (longitude and latitude)")
            elif w["unit"] in T.UNIT_M:
                crs = _with_merc(_crs("projected", w["unit"], f"the layer's .prj ({w['name']})", w["unit_m"]), w)
            elif w["unit"] == "other":
                crs = _crs("projected", "other", f"the layer's .prj ({w['name']}, {w['unit_m']:g} m)", w["unit_m"])
            else:
                crs = _crs("unknown", None, "the layer's .prj names no unit this tool reads")
        else:
            crs = _crs("unknown", None, "the layer has no .prj, so its coordinates' unit is not known")
        try:
            r = shapefile.Reader(shp=io.BytesIO(shp), shx=io.BytesIO(shx) if shx is not None else None,
                                 dbf=io.BytesIO(dbf) if dbf is not None else None, encoding=enc)
        except Exception as e:  # noqa: BLE001  (pyshp raises its own and struct errors on a broken file)
            raise GisRefused(f"The shapefile {lname} cannot be read ({type(e).__name__}). Export it again.") from None
        try:
            n_records = len(r)
        except Exception:  # noqa: BLE001
            n_records = 0
        if n_records > MAX_FEATURES:     # the review of task 17 and its final verification: a large layer of any
            skipped.append((lname, f"{lname} ({n_records:,} features)"))   # type beside the mains refused the zip
            continue
        try:
            shapes = r.shapes()
            if dbf is not None:
                try:
                    recs = r.records(deleted_as_None=True)
                except UnicodeDecodeError:
                    r = shapefile.Reader(shp=io.BytesIO(shp), shx=io.BytesIO(shx) if shx is not None else None,
                                         dbf=io.BytesIO(dbf), encoding="latin-1")
                    shapes, recs = r.shapes(), r.records(deleted_as_None=True)
                    warns.append(f"The layer's text is not valid {enc}: it is read as latin-1. Check accented names.")
                    enc = "latin-1"
                fnames = [f[0] if isinstance(f, (list, tuple)) else getattr(f, "name", str(f)) for f in r.fields[1:]]
            else:
                recs, fnames = [None] * len(shapes), []
                warns.append("The layer has no .dbf, so it has no attributes (no pipe sizes). Zip its .dbf with it.")
        except Exception as e:  # noqa: BLE001
            raise GisRefused(f"The shapefile {lname} cannot be read ({type(e).__name__}: {str(e)[:80]}). Export it "
                             f"again.") from None
        if len(shapes) > MAX_FEATURES:
            raise GisRefused(f"The layer {lname} has more than {MAX_FEATURES:,} features.")
        lay = Layer(name=lname, source="shapefile", features=[], fields=fnames, crs=crs, notes=notes, warnings=warns,
                    encoding=enc)
        if dbf is not None and len(recs) != len(shapes):
            lay.warnings.append(f"The .shp has {len(shapes)} shapes and the .dbf {len(recs)} records: they are matched "
                                f"by position up to the shorter, and the rest are left out.")
        null = deleted = 0
        for k, (s, rec) in enumerate(zip(shapes, recs if dbf is not None else [None] * len(shapes))):
            if dbf is not None and rec is None:
                deleted += 1
                continue
            props = dict(zip(fnames, list(rec))) if rec is not None else {}
            st = s.shapeType
            if st == 0 or not s.points:
                null += 1
                continue
            P = np.asarray([p[:2] for p in s.points], float)
            starts = list(s.parts) + [len(P)]
            parts = [P[a:b] for a, b in zip(starts, starts[1:]) if b > a]
            kind = {3: "line", 13: "line", 23: "line", 1: "point", 11: "point", 21: "point", 8: "point", 18: "point",
                    28: "point", 5: "polygon", 15: "polygon", 25: "polygon"}.get(st, "other")
            if kind == "point":
                parts = [P[i:i + 1] for i in range(len(P))]
            lay.features.append(Feature(index=k, kind=kind, parts=parts, props=props, has_z=st in (11, 13, 15, 18),
                                        has_m=st in (21, 23, 25, 28)))
        lay.dropped = {"null_shapes": null, "deleted_records": deleted}
        if null:
            lay.notes.append(f"{null} record{'s have' if null != 1 else ' has'} no shape (a null shape) and "
                             f"{'are' if null != 1 else 'is'} left out with {'their' if null != 1 else 'its'} attributes; "
                             f"the others keep theirs.")
        if deleted:
            lay.notes.append(f"{deleted} record{'s are' if deleted != 1 else ' is'} marked deleted in the .dbf and "
                             f"{'are' if deleted != 1 else 'is'} left out with {'their' if deleted != 1 else 'its'} shape.")
        if kml:            # the review of task 17: a preview .kml beside the shapefile made the zip a KMZ
            lay.notes.append(f"The zip also holds {', '.join(kml[:3])}, which was not read: the shapefile was.")
        layers.append(_finish(lay))
    if not layers:
        raise GisRefused(_over_cap_refusal("The zip", skipped))
    _note_skipped(layers, "zip's", skipped)
    return layers


# ----------------------------------------------------------------------------- GeoJSON and Esri JSON
def _coord_array(c, depth: int = 0):
    """An (n, 2) array of a GeoJSON position list, or None when it is not one."""
    try:
        a = np.asarray(c, dtype=float)
    except (TypeError, ValueError):
        return None
    if a.ndim != 2 or a.shape[1] < 2 or not np.isfinite(a[:, :2]).all():
        return None
    return a[:, :2]


def _geojson_geoms(g, depth: int = 0, bad: list | None = None) -> list:
    """The geometries of a GeoJSON geometry, as (kind, (n, 2) array).  bad: a one-item list counting the line and
    point parts whose positions cannot be read (a null position, 2D and 3D positions mixed), which are left out (the
    review of task 17: they were left out without a word)."""
    bad = bad if bad is not None else [0]
    if not isinstance(g, dict) or depth > WKB_MAX_DEPTH:
        return []
    t, c = g.get("type"), g.get("coordinates")

    def lines(cs):
        out = []
        for x in cs:
            a = _coord_array(x)
            if a is None:
                bad[0] += 1
            elif len(a):
                out.append(("line", a))
        return out
    if t == "LineString":
        return lines([c])
    if t == "MultiLineString":
        if not isinstance(c, list):
            bad[0] += 1
            return []
        return lines(c)
    if t == "Point":
        a = _coord_array([c])
        bad[0] += a is None
        return [("point", a)] if a is not None else []
    if t == "MultiPoint" and isinstance(c, list):
        pts = [_coord_array([x]) for x in c]
        bad[0] += sum(a is None for a in pts)
        return [("point", a) for a in pts if a is not None]
    if t in ("Polygon", "MultiPolygon"):
        return [("polygon", np.zeros((0, 2)))]
    if t == "GeometryCollection":
        gs = g.get("geometries")
        gs = gs if isinstance(gs, list) else []
        return [x for gg in gs[:100000] for x in _geojson_geoms(gg, depth + 1, bad)]
    return []


def _warn_unread(lay: Layer, n: int) -> None:
    if n:
        lay.warnings.append(f"{n} line{'s' if n != 1 else ''} or point{'s' if n != 1 else ''} of the layer "
                            f"{'have' if n != 1 else 'has'} positions that cannot be read (a missing coordinate, or "
                            f"2D and 3D positions mixed in one line) and {'are' if n != 1 else 'is'} left out.")


def _crs_member(d: dict) -> dict | None:
    c = d.get("crs")
    if not isinstance(c, dict):
        return None
    props = c.get("properties")
    if not isinstance(props, dict):          # the final verification of task 17: it raised AttributeError
        return _crs("unknown", None, "the file's crs member has no properties this tool reads")
    nm = str(props.get("name") or "")[:200]
    if "CRS84" in nm.upper():
        return _crs("geographic", "deg", "the file's crs member (CRS84, longitude and latitude)")
    m = re.search(r"EPSG:+(?:[0-9.]*:)?(\d+)", nm, re.I)
    if m:
        return _epsg_crs(int(m.group(1)), "the file's crs member")
    return _crs("unknown", None, f"the file's crs member ({nm[:40]!r}) is not one this tool reads")


def _decimals(vals: np.ndarray, most: int = 9) -> int | None:
    """The fewest decimals every value is written with (None: more than most)."""
    v = np.asarray(vals, float)[:20000]
    for d in range(0, most + 1):
        s = v * 10.0 ** d
        if np.all(np.abs(s - np.round(s)) <= 1e-6 * np.maximum(1.0, np.abs(s)) * 1e-3 + 1e-7):
            return d
    return None


def read_geojson(data: bytes, name: str = "layer.geojson") -> list[Layer]:
    if len(data) > JSON_MAX_BYTES:
        raise GisRefused(f"The file is larger than {JSON_MAX_BYTES // (1024 * 1024)} MB. Export the water mains alone, or "
                         f"as a GeoPackage.")
    text, _ = decode_text(data)
    try:
        d = json.loads(text)
    except RecursionError:
        raise GisRefused("The JSON file is nested too deep to read.") from None
    except ValueError as e:
        raise GisRefused(f"The file is not valid JSON ({str(e)[:80]}). Export the layer again.") from None
    if not isinstance(d, dict):
        raise GisRefused("The JSON file is not a GeoJSON layer (no FeatureCollection or Feature).")
    stem = os.path.splitext(os.path.basename(name or "layer"))[0]
    if "geometryType" in d or ("features" in d and isinstance(d.get("features"), list) and d["features"]
                                and isinstance(d["features"][0], dict) and "attributes" in d["features"][0]):
        return [_read_esri(d, stem)]
    if d.get("type") == "FeatureCollection":
        feats = d.get("features")
        if not isinstance(feats, list):
            raise GisRefused("The GeoJSON's features member is not a list.")
    elif d.get("type") == "Feature":
        feats = [d]
    elif d.get("type") in ("LineString", "MultiLineString", "Point", "MultiPoint", "GeometryCollection"):
        feats = [{"type": "Feature", "properties": {}, "geometry": d}]
    else:
        raise GisRefused("The JSON file is not a GeoJSON layer (no FeatureCollection or Feature).")
    if len(feats) > MAX_FEATURES:
        raise GisRefused(f"The layer has more than {MAX_FEATURES:,} features.")
    lay = Layer(name=str(d.get("name") or stem), source="geojson", features=[], fields=[], crs=_crs_member(d) or {})
    fields: dict = {}
    unread = [0]
    for k, f in enumerate(feats):
        if not isinstance(f, dict):
            continue
        props = f.get("properties") if isinstance(f.get("properties"), dict) else {}
        props = {str(a): v for a, v in props.items() if not isinstance(v, (dict, list))}
        if "id" in f and "id" not in props and isinstance(f["id"], (str, int, float)):
            props["id"] = f["id"]
        for a in props:
            fields.setdefault(a, None)
        geoms = _geojson_geoms(f.get("geometry"), 0, unread)
        lay.features.append(_feature(k, geoms, props, {}))
    lay.fields = list(fields)
    _warn_unread(lay, unread[0])
    if not lay.crs:
        xy = [p for ft in lay.features for part in ft.parts for p in part]
        a = np.asarray(xy, float) if xy else np.zeros((0, 2))
        if len(a) and np.all(np.abs(a[:, 0]) <= 180) and np.all(np.abs(a[:, 1]) <= 90):
            lay.crs = _crs("geographic", "deg", "GeoJSON's default (RFC 7946: longitude and latitude, WGS84), since "
                                                 "the file names no other and its coordinates are in range")
        else:
            lay.crs = _crs("unknown", None, "the file names no coordinate system and its coordinates are not "
                                            "longitude and latitude")
    _finish(lay)
    _quantization(lay)
    return [lay]


def _esri_crs(sr) -> dict:
    """The coordinate system of an Esri JSON spatialReference: {'wkid' or 'latestWkid'}, {'wkt'}, or (the review of
    task 17, which the sniffer accepts) the code alone, as a number or a string."""
    if isinstance(sr, bool):
        sr = None
    if isinstance(sr, (int, float)) or (isinstance(sr, str) and sr.strip().isdigit()):
        sr = {"wkid": int(float(sr))}
    sr = sr if isinstance(sr, dict) else {}
    code = sr.get("latestWkid") or sr.get("wkid")
    if isinstance(code, str) and code.strip().isdigit():
        code = int(code)
    if isinstance(code, int) and not isinstance(code, bool):
        return _epsg_crs(code, "the file's spatialReference")
    if isinstance(sr.get("wkt"), str):
        w = wkt_unit(sr["wkt"])
        if w["kind"] == "geographic":
            return _crs("geographic", "deg", "the file's spatialReference (longitude and latitude)")
        if w["unit"] in T.UNIT_M:
            return _with_merc(_crs("projected", w["unit"], "the file's spatialReference", w["unit_m"]), w)
    return _crs("unknown", None, "the file's spatialReference names no unit this tool reads")


def _read_esri(d: dict, stem: str) -> Layer:
    feats = d.get("features") or []
    if not isinstance(feats, list):
        raise GisRefused("The Esri JSON file's features member is not a list.")
    if len(feats) > MAX_FEATURES:        # the review of task 17: the rest were cut off without a word
        raise GisRefused(f"The layer has more than {MAX_FEATURES:,} features.")
    lay = Layer(name=stem, source="esri_json", features=[], fields=[], crs=_esri_crs(d.get("spatialReference")))
    fields: dict = {}
    unread = 0
    for k, f in enumerate(feats):
        if not isinstance(f, dict):
            continue
        attrs = f.get("attributes")         # the final verification of task 17: a list or a string here, or paths
        attrs = attrs if isinstance(attrs, dict) else {}       # that are no list, raised in the reader
        props = {str(a): v for a, v in attrs.items() if not isinstance(v, (dict, list))}
        for a in props:
            fields.setdefault(a, None)
        g = f.get("geometry") or {}
        g = g if isinstance(g, dict) else {}
        geoms = []
        paths = g.get("paths")
        if paths is not None and not isinstance(paths, list):
            unread += 1
            paths = []
        for path in (paths or [])[:100000]:
            a = _coord_array(path)
            if a is None:
                unread += 1
            elif len(a):
                geoms.append(("line", a))
        if "x" in g and "y" in g:
            a = _coord_array([[g["x"], g["y"]]])
            if a is not None:
                geoms.append(("point", a))
            else:
                unread += 1
        lay.features.append(_feature(k, geoms, props, {}))
    lay.fields = list(fields)
    _warn_unread(lay, unread)
    _finish(lay)
    _quantization(lay)
    return lay


def _quantization(lay: Layer) -> None:
    xy = [p for ft in lay.features for part in ft.parts for p in part]
    if not xy:
        return
    a = np.asarray(xy, float)
    d = _decimals(a.ravel())
    if d is None:
        return
    if lay.crs.get("kind") == "geographic" and d <= QUANTIZED_DEGREES:
        step = 10.0 ** -d * DEG_TO_M
        lay.warnings.append(f"Every coordinate is written with {d} decimal{'s' if d != 1 else ''} of a degree or "
                            f"fewer, about {step:.0f} m apart at most: line ends that meet may be rounded apart by that "
                            f"much. Export the layer with more decimals (6 or more), or raise the snapping tolerance and "
                            f"check the result.")
        lay.crs["quantized_decimals"] = d
    elif lay.crs.get("unit_m") and d == 0:
        step = lay.crs["unit_m"]
        if step >= 0.3:
            lay.notes.append(f"Every coordinate is a whole number of {T.UNIT_WORDS.get(lay.crs.get('unit'), 'units')} "
                             f"(rounded to {step:.2f} m).")


# ----------------------------------------------------------------------------- KML and KMZ
def _strip_ns(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _strip_tags(text: str) -> str:
    """text with every tag ('<' up to the next '>', with something between) replaced by a space, in one pass (the
    final verification of task 17: re.sub(r"<[^>]+>") rescanned to the end from every '<' with no '>' after it, 24 s
    for a description of 1.2 MB).  A '<' with no '>' after it leaves the rest as text, as that pattern did."""
    out, i, n = [], 0, len(text)
    while i < n:
        j = text.find("<", i)
        if j < 0:
            out.append(text[i:])
            break
        k = text.find(">", j + 1)
        if k < 0:
            out.append(text[i:])
            break
        if k == j + 1:                       # '<>' is no tag
            out.append(text[i:k + 1])
        else:
            out.append(text[i:j])
            out.append(" ")
        i = k + 1
    return "".join(out)


def _text_of(el) -> str:
    return _strip_tags("".join(el.itertext())[:TEXT_MAX_CHARS]).strip() if el is not None else ""


def _tags(html: str):
    """(tag name in lower case, is an end tag, start, end) of each tag of html, then ('', False, n, n), in one pass;
    the text between two tags is html[end of one:start of the next]."""
    i, n = 0, len(html)
    while i < n:
        j = html.find("<", i)
        if j < 0:
            break
        k = html.find(">", j + 1)
        if k < 0:
            break
        body = html[j + 1:k].lstrip()
        end = body.startswith("/")
        body = body[1:].lstrip() if end else body
        m = 0
        while m < len(body) and body[m].isalnum():
            m += 1
        yield body[:m].lower(), end, j, k + 1
        i = k + 1
    yield "", False, n, n


_FIELD_NAME = re.compile(r"[A-Za-z_][\w .\-]{0,39}")


def description_fields(html: str) -> dict:
    """The 'FIELD value' rows of an attribute table in a placemark's description, as ArcGIS's Layer To KML writes it
    (<tr><td>DIAMETER</td><td>8</td></tr>, the first cell a td or a th), {field: value}; {} when it holds none (the
    review of task 17: such a description parsed as 'two sizes', and the facility ID in the name became the size).
    One pass over the tags (the final verification of task 17: a pattern of nested lazy groups took time quadratic in
    the rows left open, 41 s for a description of 961 KB); a row counts only when it is exactly two cells with nothing
    but spaces between them."""
    import html as _html
    out = {}
    text = (html or "")[:TEXT_MAX_CHARS]
    row = cell = None           # row: [cells, ok]; cell: [tag, text pieces]
    prev_end = 0
    for name, end, a, b in _tags(text):
        between = text[prev_end:a]
        prev_end = b
        if cell is not None:
            cell[1].append(between)
        elif row is not None and between.strip():
            row[1] = False
        if name in ("td", "th"):
            if not end:
                if cell is not None and row is not None:
                    row[1] = False                       # a cell opened inside a cell
                cell = [name, []] if row is not None else None
            elif cell is not None:
                row[0].append((cell[0], "".join(cell[1])))
                cell = None
        elif name == "tr":
            if not end:
                row, cell = [[], True], None
            elif row is not None:
                cells, ok = row
                if ok and cell is None and len(cells) == 2 and cells[1][0] == "td":
                    k = _html.unescape(_strip_tags(cells[0][1])).strip()
                    v = _html.unescape(_strip_tags(cells[1][1])).strip()
                    if _FIELD_NAME.fullmatch(k) and k not in out:
                        out[k] = v
                row = cell = None
        elif cell is not None and name:
            cell[1].append(" ")                          # a tag inside a cell, as a space
    return out


def _kml_points(text: str, budget: list) -> list[tuple[float, float]]:
    """The (longitude, latitude) tuples of a KML coordinates element, in one pass: tuples are separated by spaces,
    and 'lon, lat' written with a space after the comma is one tuple (the review of task 17).  The final verification
    of task 17: rewriting the spaces around each comma with a pattern took time quadratic in a run of spaces, 46 s
    for 300,000 of them.  budget: a one-item list of the points the file may still hold (WKB_MAX_POINTS a file)."""
    pts, cur = [], []

    def flush():
        v = "".join(cur).split(",")
        try:
            pts.append((float(v[0]), float(v[1])))
        except (IndexError, ValueError):
            pass
    for tok in text.split():
        if cur and len(cur) < 8 and (cur[-1].endswith(",") or tok.startswith(",")):
            cur.append(tok)
            continue
        if cur:
            flush()
            if len(pts) > budget[0]:
                break
        cur = [tok]
    if cur:
        flush()
    budget[0] -= len(pts)
    if budget[0] < 0:
        raise GisRefused(f"The KML holds more than {WKB_MAX_POINTS:,} points, the most this tool reads from one file. "
                         f"Export the water mains alone.")
    return pts


class _Doctype(Exception):
    pass


def _refuse_doctype(xml_bytes: bytes) -> None:
    """GisRefused when the document declares a document type or an entity, wherever it is (the final verification of
    task 17: the guard looked at the first 200,000 characters only, so a KMZ whose doc.kml began with 200 KB of blank
    lines had its entity expanded).  An expat pass with only these handlers set stops at the declaration, before any
    entity is expanded."""
    import xml.parsers.expat as X

    def stop(*_a):
        raise _Doctype()
    p = X.ParserCreate()
    p.StartDoctypeDeclHandler = stop
    p.EntityDeclHandler = stop
    try:
        p.Parse(xml_bytes, True)
    except _Doctype:
        raise GisRefused("The KML declares a document type (<!DOCTYPE>), which is not read: it can hide entity "
                         "expansion. Save it again from Google Earth (Save Place As, KMZ).") from None
    except X.ExpatError as e:
        raise GisRefused(f"The KML cannot be read ({str(e)[:80]}). Save it again from Google Earth.") from None


def read_kml(data: bytes, name: str = "layer.kml", max_unpacked: int = ZIP_MAX_TOTAL) -> list[Layer]:
    import xml.etree.ElementTree as ET
    if data[:2] == b"PK":
        try:
            infos = zip_listing(data, max_unpacked)
        except ZipRefused as e:
            raise GisRefused(f"The KMZ is refused: {e}.") from None
        kmls = sorted((i.filename for i in infos if i.filename.lower().endswith(".kml") and not zip_junk(i.filename)),
                      key=lambda n: (n.count("/"), n.lower() != "doc.kml", n))
        if not kmls:
            raise GisRefused("The KMZ holds no .kml document.")
        try:
            data = read_zip_member(data, kmls[0], limit=max_unpacked)
        except ZipRefused as e:
            raise GisRefused(f"The KMZ is refused: {e}.") from None
        name = kmls[0]
    text, _ = decode_text(data)
    text = text.lstrip("\ufeff \t\r\n")        # a KML with a blank line before its declaration (the review of task 17)
    xml_bytes = (text.encode("utf-8") if not text.startswith("<?xml") else
                 re.sub(r"^<\?xml[^>]*\?>", "", text).encode("utf-8"))
    _refuse_doctype(xml_bytes)                 # anywhere in the document (the final verification of task 17)
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as e:
        raise GisRefused(f"The KML cannot be read ({str(e)[:80]}). Save it again from Google Earth.") from None
    lay = Layer(name=os.path.splitext(os.path.basename(name))[0], source="kml", features=[], fields=[],
                crs=_crs("geographic", "deg", "KML is always longitude and latitude (WGS84)"))
    fields: dict = {"name": None, "description": None}
    k = 0
    budget = [WKB_MAX_POINTS]
    for pm in root.iter():
        if _strip_ns(pm.tag) != "Placemark":
            continue
        props = {}
        geoms = []
        for ch in pm.iter():
            t = _strip_ns(ch.tag)
            if t == "name" and "name" not in props:
                props["name"] = (ch.text or "")[:TEXT_MAX_CHARS].strip()
            elif t == "description" and "description" not in props:
                props["description"] = _text_of(ch)
                table = description_fields("".join(ch.itertext())[:TEXT_MAX_CHARS])
                for a, v in table.items():           # an attribute table: its fields, after the KML's own
                    props.setdefault(a, v)
            elif t == "Data" and ch.get("name"):
                v = next((c for c in ch if _strip_ns(c.tag) == "value"), None)
                props[ch.get("name")] = (v.text or "").strip() if v is not None else ""
            elif t == "SimpleData" and ch.get("name"):
                props[ch.get("name")] = (ch.text or "").strip()
            elif t in ("LineString", "Point"):
                c = next((x for x in ch if _strip_ns(x.tag) == "coordinates"), None)
                if c is None or not (c.text or "").strip():
                    continue
                pts = _kml_points(c.text, budget)
                if pts:
                    a = np.asarray(pts, float)
                    if np.isfinite(a).all():
                        geoms.append(("line" if t == "LineString" else "point", a if t == "LineString" else a[:1]))
        for a in props:
            fields.setdefault(a, None)
        lay.features.append(_feature(k, geoms, props, {}))
        k += 1
        if k > MAX_FEATURES:
            raise GisRefused(f"The KML has more than {MAX_FEATURES:,} placemarks.")
    lay.fields = list(fields)
    _finish(lay)
    _quantization(lay)
    return [lay]


# ----------------------------------------------------------------------------- any GIS file
_ZIPPED_LAYER_EXT = (".gpkg", ".geojson", ".json", ".kml")


def _zipped_layer(data: bytes, max_unpacked: int) -> tuple[bytes, str] | None:
    """The one GeoPackage, GeoJSON or KML in a zip that holds no shapefile, as (its bytes, its name), read within the
    zip limits; None when the zip holds none; GisRefused when it holds several (the final verification of task 17: a
    zipped GeoPackage was refused with a sentence that named GeoPackage among the layers read)."""
    try:
        infos = zip_listing(data, max_unpacked)
    except ZipRefused as e:
        raise GisRefused(f"The zip is refused: {e}. Nothing in it was opened.") from None
    names = [i.filename for i in infos if not i.is_dir() and not zip_junk(i.filename)
             and i.filename.lower().endswith(_ZIPPED_LAYER_EXT)]
    if not names:
        return None
    if len(names) > 1:
        raise GisRefused(f"The zip holds {len(names)} GIS files ({', '.join(os.path.basename(n) for n in names[:4])}), "
                         f"and this tool reads one layer: unzip it and give the water mains file itself.")
    try:
        inner = read_zip_member(data, names[0], limit=max_unpacked)
    except ZipRefused as e:
        raise GisRefused(f"The zip is refused: {e}.") from None
    if inner[:2] == b"PK":
        raise GisRefused(f"The zip holds {os.path.basename(names[0])}, itself a zip: unzip it and give the water mains "
                         f"file itself.")
    return inner, os.path.basename(names[0])


def read_layers(data: bytes, name: str = "layer", max_unpacked: int = ZIP_MAX_TOTAL) -> tuple[list[Layer], str]:
    """(layers, kind) of a GIS file identified by its content (importers.sniff); GisRefused for anything else, with
    the reason in plain words.  max_unpacked: a zip's members together, and each, once unpacked (the hosted demo
    passes its upload limit; the review of task 17).  A zip that holds one GeoPackage, GeoJSON or KML and no shapefile
    is read as that file (the final verification of task 17)."""
    s = sniff(data, name)
    kind = s.kind
    if kind == "file_gdb_zip":
        raise GisRefused("This is a zipped Esri file geodatabase (.gdb), which cannot be read here. In ArcGIS Pro, export "
                         "the water mains layer as a shapefile (zip its .shp, .shx, .dbf and .prj together) or a "
                         "GeoPackage, and use that.")
    if kind == "geopackage":
        return read_geopackage(data, name), kind
    if kind == "shapefile_zip":
        return read_shapefile_zip(data, name, max_unpacked), kind
    if kind in ("geojson", "esri_json"):
        return read_geojson(data, name), kind
    if kind in ("kml", "kmz"):
        return read_kml(data, name, max_unpacked), kind
    if kind == "shp":
        raise GisRefused("A .shp alone has no attributes and no coordinate system: zip it with its .shx, .dbf and .prj "
                         "and use the zip.")
    if kind in ("zip_other", "archive") and data[:2] == b"PK":
        inner = _zipped_layer(data, max_unpacked)
        if inner is not None:
            got, k = read_layers(inner[0], inner[1], max_unpacked)
            if k in ("geopackage", "geojson", "esri_json", "kml"):
                for lay in got:
                    lay.notes.append(f"Read from {inner[1]}, the one GIS file in the zip.")
                return got, k
        raise GisRefused("The zip holds no GIS layer this tool reads: a shapefile needs its .shp, .shx and .dbf zipped "
                         "together (with its .prj). A GeoPackage, GeoJSON or KML file can be given as it is, unzipped.")
    if kind == "zip_other" or kind == "archive":
        raise GisRefused("The file is an archive this tool does not open: give the GIS layer itself (a GeoPackage, "
                         "GeoJSON or KML), or a zipped shapefile.")
    raise GisRefused(f"{s.sentence()} It is not a GIS layer this tool reads: use a GeoPackage, a zipped shapefile, "
                     f"GeoJSON or a Google Earth KML or KMZ.")


# ----------------------------------------------------------------------------- fields
def _norm(h: str) -> str:
    return re.sub(r"[\s_\-.()/]+", " ", str(h).strip().lower()).strip()


def _guess(fields: list[str], exact: tuple) -> str | None:
    norm = {f: _norm(f) for f in fields}
    for e in exact:
        for f, n in norm.items():
            if n == e or n.replace(" ", "") == e.replace(" ", ""):
                return f
    return None


def guess_fields(lay: Layer) -> dict:
    """The fields guessed for each role, for the operator to confirm: {'diameter', 'material', 'year', 'length', 'id',
    'c', 'status'}, each a field name or None.  For KML, the name and the description are tried for the size."""
    fields = [f for f in lay.fields if _norm(f) not in _NOT_ID or _norm(f) in ("name",)]
    g = {"diameter": _guess(fields, _DIAMETER_EXACT), "material": _guess(fields, _MATERIAL_EXACT),
         "year": _guess(fields, _YEAR_EXACT), "length": _guess(fields, _LENGTH_EXACT),
         "id": _guess([f for f in fields if _norm(f) not in _NOT_ID], _ID_EXACT), "c": _guess(fields, _C_EXACT),
         "status": _guess(fields, _STATUS_EXACT)}
    if g["diameter"] is None:
        # a KML's name or description is its size only when most of them carry a unit (8 in, 8", 200 mm, DN200): the
        # review of task 17 found paths numbered by their author ('Pipe 2') read as 2 in, and facility IDs ('W-12')
        for f in ("description", "name") if lay.source == "kml" else ():
            vals = [ft.props.get(f) for ft in lay.features if ft.kind == "line"]
            parsed = [T.parse_diameter(v, strict_in=True) for v in vals]   # 'in Oak St' is no unit (the final
                                                                              # verification of task 17)
            if vals and sum(d["value"] is not None and d["unit"] is not None for d in parsed) >= 0.5 * len(vals):
                g["diameter"] = f
                break
    if lay.source == "kml" and g["id"] == "name":
        g["id"] = None                     # a placemark's name is a label, not an ID
    return g


def _layer_lines(lay: Layer, crs: dict) -> tuple[list[T.Line], str, dict]:
    """The layer's lines in metres, in the coordinate system crs settles (the layer's own, or the one the operator
    picked): (lines, id_style, frame) where frame says how metres map back to the coordinates written in the skeleton
    (feet)."""
    lines = []
    if crs.get("kind") == "geographic":
        allp = _line_points(lay)
        if not len(allp):                # the review of task 17: no line point at all raised IndexError here
            return [], "lonlat", {"kind": "enu", "lat0": 0.0, "lon0": 0.0}
        lat0, lon0 = plane_centre(allp[:, 1], allp[:, 0])
        for ft in lay.features:
            if ft.kind != "line":
                continue
            for k, part in enumerate(ft.parts):
                if len(part) < 2:
                    continue
                xy = enu_from_latlon(part[:, 1], part[:, 0], lat0, lon0)
                lines.append(T.Line(ft.index, k, xy, np.array(part, float)))
        return lines, "lonlat", {"kind": "enu", "lat0": lat0, "lon0": lon0}
    f = crs.get("unit_m")
    allp = _line_points(lay)
    origin = np.floor(allp.min(axis=0)) if len(allp) else np.zeros(2)
    for ft in lay.features:
        if ft.kind != "line":
            continue
        for k, part in enumerate(ft.parts):
            if len(part) < 2:
                continue
            raw = np.array(part, float)
            lines.append(T.Line(ft.index, k, (raw - origin) * f, raw))
    return lines, "xy", {"kind": "projected", "origin": origin.tolist(), "unit_m": f}


@dataclass
class GisOptions:
    layer: str | None = None
    unit: str | None = None            # 'ft', 'us_ft' or 'm', for a layer whose unit is not known
    tolerance_m: float | None = None
    fields: dict = field(default_factory=dict)   # role -> field name ('' for none), over the guesses
    bare_unit: str | None = None       # 'in', 'mm' or 'm', for bare diameters the range cannot settle
    default_diameter_in: float | None = None     # for pipes whose size cannot be read (default: the commonest size)
    known_lengths_ft: dict = field(default_factory=dict)   # {pipe ID: its true length in feet}, the scale check
    to_scale: bool | None = None       # the operator's answer, for a layer read in a unit they picked
    swap_axes: bool | None = None      # the operator's answer, for a layer whose axes may be swapped (True: swap them)


@dataclass
class GisBuild:
    name: str
    kind: str
    layers: list
    layer: Layer | None = None
    questions: list = field(default_factory=list)   # {'key', 'text', 'options', 'blocking'}
    refusal: str | None = None
    net: T.Network | None = None
    tolerance_m: float | None = None
    diagnostics: list = field(default_factory=list)
    fields: dict = field(default_factory=dict)
    guess: dict = field(default_factory=dict)
    diam_in: dict = field(default_factory=dict)
    diam_counts: dict = field(default_factory=dict)
    c_of: dict = field(default_factory=dict)
    c_counts: dict = field(default_factory=dict)
    c_material: dict = field(default_factory=dict)  # task 19: pipe ID -> the material its material field names
                                                    # (topology.material_of), or None, as for a pipe whose C came from
                                                    # the layer's own C field (the operator form's pipe_c rows set C by
                                                    # material and never touch a C the layer gives)
    c_from_field: list = field(default_factory=list)   # task 19's review: the pipes whose C the layer's C field gave
    id_table: list = field(default_factory=list)
    length_check: dict = field(default_factory=dict)
    scale_check: dict = field(default_factory=dict)
    years: dict = field(default_factory=dict)
    inp_text: str | None = None
    notes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    frame: dict = field(default_factory=dict)
    wntr_schema: bool = False
    offer_unit: bool = False           # the page shows the coordinates' unit even for a layer that names one: it was
                                       # refused as crowded, or lies in longitude and latitude ranges while its file
                                       # names a projected system (the final verification of task 17)

    @property
    def ok(self) -> bool:
        return self.inp_text is not None

    def blocking_questions(self) -> list:
        return [q for q in self.questions if q.get("blocking")]


def _ask(b: GisBuild, key: str, text: str, options=None, blocking: bool = True) -> None:
    b.questions.append({"key": key, "text": text, "options": list(options or []), "blocking": blocking})


def build_from_layer(lay: Layer, opts: GisOptions | None = None, name: str = "layer", kind: str = "",
                     layers: list | None = None) -> GisBuild:
    """One layer built into a skeleton: its unit settled (or asked), its lines snapped at the tolerance, its fields
    read, its skeleton .inp written (GisBuild.inp_text; None when a blocking question is open or nothing can be
    built).  It never raises (the review of task 17: the Get started page called it bare, and one NaN vertex gave a
    traceback): anything it did not foresee is a plain refusal."""
    try:
        return _build_from_layer(lay, opts, name, kind, layers)
    except Exception as e:  # noqa: BLE001  (the last guard: a layer the builder did not foresee)
        return GisBuild(name=name, kind=kind, layers=[l.name for l in (layers or [lay])], layer=lay,
                        refusal=f"The layer could not be built into pipes ({type(e).__name__}). Export the layer "
                                f"again.")


UNIT_HINT = ("The coordinates' unit or coordinate system is probably wrong (longitude and latitude read as feet or "
             "metres, for example), or the joining distance is far too large for this layer. Check the unit, or "
             "export the layer again with its coordinate system.")


def _build_from_layer(lay: Layer, opts: GisOptions | None, name: str, kind: str, layers: list | None) -> GisBuild:
    opts = opts or GisOptions()
    b = GisBuild(name=name, kind=kind, layers=[l.name for l in (layers or [lay])], layer=lay)
    b.notes += lay.notes
    b.warnings += lay.warnings
    lines_n = lay.n_lines
    if not lay.features:
        b.refusal = f"The layer {lay.name} has no features: it is empty. Export the water mains layer with its lines."
        return b
    if lines_n == 0:
        what = {"point": "points", "polygon": "areas (polygons)", "none": "no geometry"}
        kinds = sorted({what.get(f.kind, f.kind) for f in lay.features})
        b.refusal = (f"The layer {lay.name} holds {', '.join(kinds)}, not lines, so it has no pipes. A water mains layer "
                     f"is made of lines; a layer of valves or hydrants can be given beside it.")
        return b
    crs = dict(lay.crs)
    if crs.get("refusal"):              # a Mercator that is not Web Mercator, or a KML outside longitude and latitude
        b.refusal = crs["refusal"]
        return b
    if crs.get("axes_suspect"):         # the final verification of task 17: latitude written first?
        if opts.swap_axes is None:
            lat, lon = crs.get("median_lat", 0.0), crs.get("median_lon", 0.0)
            where = (f"its median latitude is {abs(lat):.1f} degrees {'north' if lat >= 0 else 'south'}, where almost "
                     f"no water network is" if not crs.get("axes_out_of_range") else
                     "most of its coordinates lie outside longitude and latitude ranges")
            _ask(b, "axes", f"Read as longitude and latitude, {where}; with each coordinate's two numbers swapped, the "
                            f"layer lies at {abs(lon):.2f} degrees {'north' if lon >= 0 else 'south'}, {abs(lat):.2f} "
                            f"degrees {'east' if lat >= 0 else 'west'}. Some files write latitude first (a known "
                            f"EPSG:4326 axis-order mix-up). Swap them?", ["swap", "keep"])
            b.fields, b.guess = {}, guess_fields(lay)
            return b
        if opts.swap_axes:
            lay = swap_axes(lay)
            b.layer, crs = lay, dict(lay.crs)
            b.notes.append("Each coordinate's two numbers were swapped, as you chose: the file writes latitude first.")
        elif crs.get("axes_out_of_range"):    # kept as written: the range test as before the axis question
            keep = Layer(name=lay.name, source=lay.source, features=lay.features, fields=lay.fields,
                         crs={k: v for k, v in crs.items() if k != "axes_suspect"})
            _check_lonlat(keep, keep.crs, axes=False)
            crs = dict(keep.crs)
            if crs.get("refusal"):
                b.refusal = crs["refusal"]
                return b
        else:
            b.warnings.append("The coordinates are read as written, longitude first, as you chose, though the layer's "
                              "median latitude is beyond any but a few water networks: check the map.")
    unit, unit_m = crs.get("unit"), crs.get("unit_m")
    b.offer_unit = crs.get("kind") == "projected" and looks_lonlat(lay)
    if crs.get("kind") != "geographic":
        if opts.unit == "deg":          # the review of task 17: longitude and latitude, chosen
            if not looks_lonlat(lay):
                big = float(np.abs(_line_points(lay)).max()) if len(_line_points(lay)) else 0.0
                b.refusal = (f"The coordinates are not longitude and latitude: they run to {big:,.0f}, outside -180 to "
                             f"180 and -90 to 90. Choose the unit they are in.")
                return b
            crs = _crs("geographic", "deg", "chosen by you (longitude and latitude, WGS84)")
            crs["picked"] = True
            unit, unit_m = "deg", None
        elif opts.unit in T.UNIT_M:
            if unit in T.UNIT_M and opts.unit != unit:
                b.warnings.append(f"The coordinates' unit was set to {T.UNIT_WORDS[opts.unit]}, over "
                                  f"{T.UNIT_WORDS[unit]} from {crs.get('how')}.")
            unit, unit_m = opts.unit, T.UNIT_M[opts.unit]
            crs.update(unit=unit, unit_m=unit_m, how=f"chosen by you ({T.UNIT_WORDS[unit]})", picked=True)
            if looks_lonlat(lay):
                b.warnings.append(f"Every coordinate lies within longitude and latitude ranges, so the layer is "
                                  f"probably in degrees, not {T.UNIT_WORDS[unit]}: if the network below is wrong, "
                                  f"choose degrees.")
        elif unit not in T.UNIT_M and not (unit == "other" and unit_m):
            hint = (" They all lie within longitude and latitude ranges, so they are probably degrees (WGS84)."
                    if looks_lonlat(lay) else "")
            _ask(b, "unit", f"In what unit are the coordinates? ({crs.get('how', 'the layer names none')}.){hint} Lines "
                            f"cannot be joined at a tolerance until the unit is known.", ["ft", "us_ft", "m", "deg"])
            b.fields, b.guess = {}, guess_fields(lay)
            return b
    b.frame = {"crs": crs}
    tol = opts.tolerance_m or (T.DEFAULT_TOL_M["kml"] if lay.source == "kml" else T.DEFAULT_TOL_M["gis"])
    b.tolerance_m = tol
    guess = guess_fields(lay)
    used = dict(guess)
    for k, v in (opts.fields or {}).items():
        if k in used:
            used[k] = v or None
    for k, v in used.items():
        if v is not None and v not in lay.fields:
            b.refusal = f"The field {v!r} chosen for the {k} is not in the layer ({', '.join(lay.fields[:12])})."
            return b
    b.fields, b.guess = used, guess
    lines, id_style, frame = _layer_lines(lay, crs)
    b.frame.update(frame)
    pipe_ids = None
    if used["id"]:
        raw_ids = {f.index: f.props.get(used["id"]) for f in lay.features if f.kind == "line"}
        vals = [v for v in raw_ids.values() if v not in (None, "")]
        if len(vals) == len(raw_ids) and len({str(v) for v in vals}) == len(vals):
            pipe_ids = {}
            for i, v in raw_ids.items():
                c, why = T.clean_id(str(v).strip() if not isinstance(v, float) or not v.is_integer() else str(int(v)))
                pipe_ids[i] = c
                if why:
                    b.id_table.append({"layer_id": str(v), "model_id": c, "why": why})
            if len(set(pipe_ids.values())) != len(pipe_ids):
                b.warnings.append(f"The ID field {used['id']!r} gives two pipes one ID once cleaned for EPANET: pipes "
                                  f"get IDs from their position instead.")
                pipe_ids, b.id_table = None, []
        else:
            b.warnings.append(f"The ID field {used['id']!r} is empty or repeated for some pipes, so pipes get IDs from "
                              f"their position instead (P and 8 characters).")
    degrees_hint = (" Every coordinate lies within longitude and latitude ranges: choose degrees (longitude and "
                    "latitude) as the unit of the coordinates." if crs.get("kind") != "geographic" and looks_lonlat(lay)
                    else "")
    try:
        net = T.build_network(lines, tol, id_style, pipe_ids)
    except T.TopologyRefused as e:
        b.refusal = str(e) + degrees_hint
        b.offer_unit = b.offer_unit or str(e).startswith(T.CROWDED[:40])   # the final verification of task 17
        return b
    r = net.report
    if not net.pipes or r["dropped_pieces"] > DROPPED_MAX_SHARE * r["pieces"]:
        # the review of task 17: a layer in longitude and latitude read as feet put every end in one node and was
        # written as an empty model, offered for download
        b.refusal = ((f"No pipe can be built from the layer's {r['lines']:,} lines" if not net.pipes else
                      f"{r['dropped_pieces']:,} of the layer's {r['pieces']:,} pieces of line would be left out")
                     + f": joining ends within {tol / T.FT:.1f} ft puts both ends of "
                     + ("every one" if not net.pipes else "each") + " in one node. " + UNIT_HINT + degrees_hint)
        b.offer_unit = True
        return b
    b.net = net
    b.warnings += net.flags
    # IDs that changed with a piece suffix are listed too
    if pipe_ids:
        for p in net.pipes:
            base = pipe_ids.get(p["feature"])
            if base is not None and p["id"] != base:
                b.id_table.append({"layer_id": base, "model_id": p["id"],
                                   "why": "a piece of the line, cut where another line ends on it"
                                   if p["pieces"] > 1 else "a part of a line with several parts"})
    feats = {f.index: f for f in lay.features}
    # diameters
    strict = lay.source == "kml" and used["diameter"] in ("name", "description")   # 'Line 3 in Oak St' is no size
    parsed = {p["id"]: T.parse_diameter(feats[p["feature"]].props.get(used["diameter"]), strict_in=strict)
              if used["diameter"] else {"value": None, "unit": None, "flag": "empty", "text": ""} for p in net.pipes}
    bare = [d["value"] for d in parsed.values() if d["value"] is not None and d["unit"] is None]
    bu = opts.bare_unit or (T.bare_unit(bare) if bare else None)
    if bare and bu is None:
        lo, hi = min(bare), max(bare)
        _ask(b, "bare_unit", f"The diameter field {used['diameter']!r} holds numbers with no unit, from {lo:g} to {hi:g}, "
                             f"which fit none of inches (2 to 48), mm (50 to 1,200) or metres (0.05 to 1.2). Which unit "
                             f"are they in?", ["in", "mm", "m"])
    counts = {"read": 0, "bare": len(bare), "bare_unit": bu, "reducer": 0, "two sizes": 0, "unreadable": 0, "empty": 0}
    sizes = []
    for pid, d in parsed.items():
        if d["value"] is not None and (d["unit"] or bu):
            v = T.to_inches(d["value"], d["unit"] or bu)
            b.diam_in[pid] = v
            sizes.append(round(v, 2))
            counts["read"] += 1
        else:
            counts[d["flag"] or "unreadable"] = counts.get(d["flag"] or "unreadable", 0) + 1
    missing = [pid for pid in parsed if pid not in b.diam_in]
    if missing and not (bare and bu is None):
        fill = opts.default_diameter_in
        if fill is None:
            fill = max(set(sizes), key=lambda s: (sizes.count(s), -s)) if sizes else None
        if fill is None:
            _ask(b, "default_diameter_in", "No pipe size can be read from the layer. What size are its pipes (inches)? "
                                           "Every pipe gets that size until you give the real ones.", [])
        else:
            for pid in missing:
                b.diam_in[pid] = float(fill)
            counts["filled_in"] = float(fill)
            ex = [f"{parsed[p]['text']!r}" for p in missing if parsed[p]["text"]][:4]
            b.warnings.append(f"{len(missing)} pipe{'s' if len(missing) != 1 else ''} of {len(parsed)} {'have' if len(missing) != 1 else 'has'} "
                              f"no size that can be read" + (f" ({', '.join(ex)}: " if ex else " (") +
                              f"{counts['reducer']} reducers, {counts['two sizes']} with two sizes, "
                              f"{counts['unreadable']} unreadable, {counts['empty']} empty), and "
                              f"{'are' if len(missing) != 1 else 'is'} given {fill:g} in, "
                              f"{'the commonest size' if opts.default_diameter_in is None else 'the size you gave'}. "
                              f"Give the real sizes in the layer.")
    b.diam_counts = counts
    # Hazen-Williams C
    cc = {"from_field": 0, "from_material": 0, "default": 0}
    classes: dict = {}
    for p in net.pipes:
        pr = feats[p["feature"]].props
        c = None
        if used["c"]:
            try:
                c = float(pr.get(used["c"]))
                c = c if 40 <= c <= 160 else None
            except (TypeError, ValueError):
                c = None
        mat = None
        if c is not None:
            cc["from_field"] += 1
            b.c_from_field.append(p["id"])
        else:
            c, cls, mat = T.material_of(pr.get(used["material"]) if used["material"] else None)
            if cls:
                cc["from_material"] += 1
                classes[cls] = classes.get(cls, 0) + 1
            else:
                cc["default"] += 1
        b.c_of[p["id"]] = c
        b.c_material[p["id"]] = mat
    cc["classes"] = classes
    cc["uniform"] = len(set(b.c_of.values())) <= 1
    b.c_counts = cc
    if used["year"]:
        ys = [feats[p["feature"]].props.get(used["year"]) for p in net.pipes]
        yv = [y for y in (year_of(v) for v in ys) if y is not None]
        b.years = {"read": len(yv), "first": min(yv) if yv else None, "last": max(yv) if yv else None}
    # lengths: geometry, cross-checked against the declared field
    if used["length"]:
        ratios = []
        multipart = {q["feature"] for q in net.pipes if q["part"] > 0}
        geographic = crs.get("kind") == "geographic"
        for p in net.pipes:
            if p["pieces"] > 1 or p["feature"] in multipart:
                continue
            v = feats[p["feature"]].props.get(used["length"])
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            if v > 0 and p["length_m"] > 0:
                geo_units = p["length_m"] / (unit_m if not geographic and unit_m else 1.0)
                ratios.append(v / geo_units)
        if ratios:
            med = float(np.median(ratios))
            near = lambda x: abs(med / x - 1) <= LENGTH_RATIO_TOL
            if geographic and near(1 / T.FT):       # longitude and latitude with a length field in feet: consistent
                med = med * T.FT
            what = ("agree with the geometry" if near(1.0) else
                    "are about 3.28 times the geometry's: the field is in feet and the coordinates in metres, or the "
                    "coordinates' unit is wrong" if near(1 / T.FT) else
                    "are about 0.305 times the geometry's: the field is in metres and the coordinates in feet, or the "
                    "coordinates' unit is wrong" if near(T.FT) else
                    f"are {med:.3g} times the geometry's, which matches no unit: check the field and the coordinates")
            b.length_check = {"field": used["length"], "pipes": len(ratios), "median_ratio": med,
                              "agrees": near(1.0), "flag": None if near(1.0) else what}
            if not near(1.0):
                b.warnings.append(f"The declared lengths in {used['length']!r} {what}. Lengths are taken from the "
                                  f"geometry, in {T.UNIT_WORDS.get(unit, 'the layer unit')}.")
    elif crs.get("picked") and crs.get("kind") != "geographic":     # longitude and latitude are to scale
        if opts.to_scale is None:
            _ask(b, "to_scale", "The layer has no length field and its unit was chosen by you: is it drawn to scale (true "
                                "lengths)? A map traced from a sketch is often not. Give two or three pipes whose true "
                                "length you know to check.", ["yes", "no"], blocking=False)
        elif opts.to_scale is False:
            b.refusal = ("The layer is not drawn to scale, so its pipes' lengths cannot be read from it. Give the real "
                         "lengths in a length field, or trace the mains again over a map that is to scale.")
            return b
    if opts.known_lengths_ft:
        rows = []
        lengths = {p["id"]: p["length_m"] / T.FT for p in net.pipes}
        for pid, true_ft in opts.known_lengths_ft.items():
            if pid in lengths and true_ft and true_ft > 0:
                rows.append({"pipe": pid, "drawn_ft": lengths[pid], "true_ft": float(true_ft),
                             "ratio": lengths[pid] / float(true_ft)})
        med = float(np.median([r["ratio"] for r in rows])) if rows else None
        b.scale_check = {"pipes": rows, "median_ratio": med, "to_scale": med is not None and abs(med - 1) <= SCALE_TOL,
                         "unknown_ids": sorted(set(opts.known_lengths_ft) - set(lengths))}
        if med is not None and abs(med - 1) > SCALE_TOL:
            b.warnings.append(f"The known lengths say the layer is drawn at {med:.3g} times true length: its lengths "
                              f"are not to scale (or the unit is wrong).")
    if b.blocking_questions():
        return b
    # the diagnostic table at half and twice the tolerance (the tolerance's own row is this build's report), only once
    # every blocking question is answered (the review of task 17: the page built each layer up to eight times)
    b.diagnostics = T.diagnostics(lines, tol, id_style, known={1.0: net.report})
    title = f"Skeleton from the GIS layer {lay.name} (ResidualMap): pipes only, not ready to run"
    if frame["kind"] == "enu":
        to_ft = lambda xy: (xy[0] / T.FT, xy[1] / T.FT)
    else:
        o, f = frame["origin"], frame["unit_m"]
        to_ft = lambda xy: ((xy[0] / f + o[0]) * f / T.FT, (xy[1] / f + o[1]) * f / T.FT)
    b.inp_text = T.write_inp(net, b.diam_in, b.c_of, title, 0.0, to_ft)
    return b


_YEAR_TOKEN = re.compile(r"(?<![0-9])(18[5-9][0-9]|19[0-9][0-9]|20[0-2][0-9])(?![0-9])")
_YMD = re.compile(r"\s*(18[5-9][0-9]|19[0-9][0-9]|20[0-2][0-9])(?:0[1-9]|1[0-2])(?:0[1-9]|[12][0-9]|3[01])\s*")
EPOCH_MS_MIN = 1e10                   # a number this large or larger (either sign) in a date field is epoch milliseconds


def year_of(v) -> int | None:
    """The install year of a field value (1850 to 2029), or None (the final verification of task 17: any four digits
    that looked like a year were taken from str(value), so Esri JSON's dates in epoch milliseconds gave wrong years):
    a date's year; a number of EPOCH_MS_MIN or more, either sign, as epoch milliseconds (how Esri JSON and ArcGIS's
    GeoJSON write a date); a whole number from 1850 to 2029 as a year; in a text, a year standing alone (no digit
    either side), the first one; a date written yyyymmdd, as a number or a text, its year."""
    import datetime as _dt
    if v is None or isinstance(v, bool):
        return None
    if isinstance(v, (_dt.date, _dt.datetime)):
        y = v.year
    elif isinstance(v, (int, float, np.integer, np.floating)):
        x = float(v)
        if not math.isfinite(x):
            return None
        if abs(x) >= EPOCH_MS_MIN:
            try:
                y = (_dt.datetime(1970, 1, 1) + _dt.timedelta(milliseconds=x)).year
            except (OverflowError, ValueError):
                return None
        elif x.is_integer():
            y = int(x) // 10000 if 18500101 <= x <= 20291231 else int(x)     # a date written yyyymmdd
        else:
            return None
    else:
        t = str(v)[:200]
        m = _YMD.fullmatch(t) or _YEAR_TOKEN.search(t)
        if not m:
            return None
        y = int(m.group(1))
    return y if 1850 <= y <= 2029 else None


def proposed_closed_valves(lay: Layer, net: T.Network | None = None, status_field: str | None = None,
                           tol_m: float | None = None, frame: dict | None = None) -> dict:
    """A valve layer's normally closed valves as proposed 'closed_valve' rows for the operator form (critique 2), for
    the operator to confirm: {'status_field', 'rows': [{'kind', 'valve', 'status', 'lon_or_x', 'lat_or_y', 'node',
    'distance_m'}], 'open', 'unread'}.  With a network (and its frame), each valve's nearest node within the tolerance
    is named."""
    g = guess_fields(lay)
    sf = status_field or g["status"]
    out = {"status_field": sf, "rows": [], "open": 0, "unread": 0, "id_field": g["id"]}
    if sf is None:
        return out
    pts = [f for f in lay.features if f.kind == "point"]
    node_xy = node_ids = None
    if net is not None and frame is not None and net.nodes:
        node_ids = sorted(net.nodes)
        node_xy = np.array([net.nodes[n]["xy"] for n in node_ids])
    for f in pts:
        v = f.props.get(sf)
        s = _norm(v) if v is not None else ""
        if not s:
            out["unread"] += 1
            continue
        if s not in _CLOSED_WORDS:
            out["open"] += 1
            continue
        x, y = (float(c) for c in f.parts[0][0])
        row = {"kind": "closed_valve", "valve": str(f.props.get(g["id"]) if g["id"] else f.index), "status": str(v),
               "x": x, "y": y, "node": None, "distance_m": None}
        if node_xy is not None:
            if frame.get("kind") == "enu":
                P = enu_from_latlon(np.array([y]), np.array([x]), frame["lat0"], frame["lon0"])[0]
            else:
                P = (np.array([x, y]) - np.array(frame["origin"])) * frame["unit_m"]
            d = np.hypot(*(node_xy - P).T)
            k = int(np.argmin(d))
            if d[k] <= (tol_m or net.tol_m):
                row.update(node=node_ids[k], distance_m=float(d[k]))
        out["rows"].append(row)
    return out


def valve_layer(layers: list, line_layer: Layer) -> Layer | None:
    """The layer whose points give proposed closed valves: the line layer itself when its points have a status field,
    else the first other layer of points (in the same GeoPackage or zip) whose fields include one (a valve layer)."""
    for lay in [line_layer] + [l for l in layers if l is not line_layer]:
        if lay.n_points and guess_fields(lay)["status"]:
            return lay
    return None


# ----------------------------------------------------------------------------- the readiness report
def readiness(b: GisBuild, prepared: Readiness) -> Readiness:
    """The skeleton's readiness report: rows that say what came from the layer, what was derived, assumed and missing,
    before importers.prepare_inp's own rows for the skeleton, which say why it cannot run yet."""
    rep = Readiness(name=prepared.name, disinfectant=prepared.disinfectant, sha1=prepared.sha1, sniff=prepared.sniff,
                    not_a_model=prepared.not_a_model, run_cost=prepared.run_cost, prepared=prepared.prepared,
                    same_file=prepared.same_file, flags={**prepared.flags, "gis_skeleton": True})
    lay, net = b.layer, b.net
    crs = b.frame.get("crs", lay.crs)
    rep.add("Built from", f"the GIS layer {lay.name} ({b.kind or lay.source}), {lay.n_lines:,} lines", FROM_FILE)
    geographic = crs.get("kind") == "geographic"
    rep.add("Coordinates", f"{'longitude and latitude' if geographic else T.UNIT_WORDS.get(crs.get('unit'), 'unknown unit')}: "
                           f"{crs.get('how')}", ASSUMED if crs.get("picked") else FROM_FILE)
    rep.add("Pipes and nodes", f"{len(net.pipes):,} pipes and {len(net.nodes):,} nodes, ends joined within "
                               f"{b.tolerance_m / T.FT:.1f} ft ({b.tolerance_m:.2f} m); {net.report['components']} "
                               f"connected part{'s' if net.report['components'] != 1 else ''}, "
                               f"{net.report['dead_ends']:,} dead ends, {net.report['tees_split']:,} "
                               f"cut{'s' if net.report['tees_split'] != 1 else ''} at a tee",
            DERIVED)
    dc = b.diam_counts
    rep.add("Diameters", f"{dc.get('read', 0):,} of {len(net.pipes):,} read from {b.fields.get('diameter')!r}"
            + (f" (bare numbers read as {dc['bare_unit']})" if dc.get("bare") else "")
            + (f"; {len(net.pipes) - dc.get('read', 0):,} given {dc['filled_in']:g} in" if dc.get("filled_in") else ""),
            FROM_FILE if dc.get("read", 0) == len(net.pipes) else ASSUMED)
    cc = b.c_counts
    if cc.get("from_field") == len(net.pipes):
        rep.add("Hazen-Williams C", f"from the field {b.fields.get('c')!r}", FROM_FILE)
    else:
        rep.add("Hazen-Williams C", f"{cc.get('from_material', 0):,} from the material (new-pipe values of the "
                                    f"{T.C_SOURCE}; {T.C_CLASSES_NOTE}), {cc.get('from_field', 0):,} from a C field, "
                                    f"{cc.get('default', 0):,} at {T.C_DEFAULT:g}, the default", ASSUMED)
    if b.years.get("read"):
        rep.add("Install year", f"read for {b.years['read']:,} pipes ({b.years['first']} to {b.years['last']}), shown only: "
                                f"it does not change C, so no pipe's age reaches the model", FROM_FILE)
    elif b.years:
        rep.add("Install year", f"no year could be read from {b.fields.get('year')!r}", FROM_FILE)
    where = ("longitude and latitude on the WGS84 ellipsoid, so lengths along the ground"
             + (" (converted from Web Mercator, whose metres are not ground metres)" if crs.get("from_web_mercator")
                else "") if geographic else
             f"{T.UNIT_WORDS.get(crs.get('unit'), 'the layer unit')}, " + ("the unit you chose: check it against a few "
                                                                            "known lengths" if crs.get("picked") else
                                                                            f"from {crs.get('how')}"))
    rep.add("Lengths", f"from the geometry, in {where}" + (f"; the declared lengths agree (median ratio "
                                                            f"{b.length_check['median_ratio']:.3f})"
                                                            if b.length_check.get("agrees") else
                                                            f"; {b.length_check['flag']}" if b.length_check else ""),
            DERIVED)
    rep.add("Elevations", "every junction at 0 ft: a skeleton has no pressure-control valve, so one elevation changes "
                          "no flow; give real elevations before adding a PRV", ASSUMED)
    rep.add("Demands", "none: the layer has no customer demand", MISSING)
    rep.add("Source, tanks, pumps, valves, controls", "none: a mains layer carries none of them", MISSING)
    for item, value, source in prepared.rows:          # the skeleton's own checks: derived from the layer, not read
        if item == "Text encoding":
            continue
        if item == "File":
            value = "the skeleton model written from your layer (EPANET .inp)"
        rep.rows.append((item, value, DERIVED if source == FROM_FILE else source))
    rep.blocking = [NOT_READY] + list(prepared.blocking)
    rep.fixes = list(prepared.fixes)
    rep.warnings = list(b.warnings) + list(prepared.warnings)
    return rep


NOT_READY = ("Not ready to run: a mains layer gives pipes, their sizes and where they meet, but no source, tank, pump, "
             "valve, control or demand. To run it, download the model and add those in EPANET (or the modelling "
             "program you use); a form for them here is not built yet.")
# the hand-off from the Get started page to the main page (st.session_state, shared by a session's pages)
PREPARED_KEY = "residualmap_prepared_model"     # {'name', 'data' (the .inp bytes), 'from' (a plain phrase)}
PREPARED_OPEN = "residualmap_prepared_open"     # set by the page's button: the main page picks the prepared model once
PREPARED_OPTION = "Model prepared on Get started"
MODEL_RADIO_KEY = "residualmap_model_source"    # the main page's Model radio, keyed only while a prepared model exists


# ----------------------------------------------------------------------------- WNTR's own GIS schema
WNTR_LAYERS = ("junctions", "tanks", "reservoirs", "pipes", "pumps", "valves")
_WNTR_NODE_TYPE = {"junctions": "Junction", "tanks": "Tank", "reservoirs": "Reservoir"}
_WNTR_LINK_TYPE = {"pipes": "Pipe", "pumps": "Pump", "valves": "Valve"}


def wntr_schema_members(data: bytes, max_unpacked: int = ZIP_MAX_TOTAL) -> dict | None:
    """{'layers': {element: member}, 'sidecar': member} of a zip that holds WNTR's GIS schema (GeoJSON files named
    <prefix>_junctions.geojson and so on, the pipes with start_node_name and end_node_name) and one JSON file with the
    model's options and patterns (WNTR's JSON); None otherwise."""
    try:
        infos = zip_listing(data, max_unpacked)
    except ZipRefused:
        return None
    names = [i.filename for i in infos if not i.is_dir() and not zip_junk(i.filename)]
    layers = {}
    for n in names:
        m = re.search(r"(?:^|[_/-])(junctions|tanks|reservoirs|pipes|pumps|valves)\.(?:geo)?json$", n, re.I)
        if m:
            layers[m.group(1).lower()] = n
    side = [n for n in names if n.lower().endswith(".json") and n not in layers.values()]
    if "pipes" not in layers or "junctions" not in layers or len(side) != 1:
        return None
    try:
        pipes = json.loads(decode_text(read_zip_member(data, layers["pipes"], limit=max_unpacked))[0])
        feats = pipes.get("features") or []
        if not feats or not {"start_node_name", "end_node_name"} <= set(feats[0].get("properties") or {}):
            return None
    except (ValueError, ZipRefused, AttributeError, RecursionError):
        return None
    return {"layers": layers, "sidecar": side[0]}


def wntr_schema_dict(data: bytes, fix_patterns: bool = True, max_unpacked: int = ZIP_MAX_TOTAL) -> dict:
    """The model dict (wntr.network.to_dict's shape) of a WNTR-schema zip: nodes and links from the layers, every other
    table (options, patterns, curves, controls, sources) from the sidecar, whose own nodes and links are not used.
    fix_patterns copies each junction's demand_pattern into pattern_name, which from_dict reads."""
    mem = wntr_schema_members(data, max_unpacked)
    if mem is None:
        raise GisRefused("This zip is not a layer in WNTR's GIS schema with its JSON file.")
    side = json.loads(decode_text(read_zip_member(data, mem["sidecar"], limit=max_unpacked))[0])
    d = {k: v for k, v in side.items() if k not in ("nodes", "links")}
    nodes, links = [], []
    for el, member in mem["layers"].items():
        fc = json.loads(decode_text(read_zip_member(data, member, limit=max_unpacked))[0])
        for f in fc.get("features") or []:
            pr = dict(f.get("properties") or {})
            g = f.get("geometry") or {}
            c = g.get("coordinates")
            if el in _WNTR_NODE_TYPE:
                node = {**pr, "node_type": _WNTR_NODE_TYPE[el], "coordinates": tuple(float(v) for v in c[:2])}
                if fix_patterns and el == "junctions" and "demand_pattern" in node and "pattern_name" not in node:
                    node["pattern_name"] = node["demand_pattern"]
                nodes.append(node)
            else:
                link = {**pr, "link_type": _WNTR_LINK_TYPE[el]}
                if g.get("type") == "LineString" and isinstance(c, list) and len(c) > 2:
                    link["vertices"] = [tuple(float(v) for v in p[:2]) for p in c[1:-1]]
                else:
                    link["vertices"] = []
                links.append(link)
    d["nodes"], d["links"] = nodes, links
    return d


def wntr_schema_inp(data: bytes, fix_patterns: bool = True, max_unpacked: int = ZIP_MAX_TOTAL) -> str:
    """The EPANET .inp text of a WNTR-schema zip, written by WNTR from the model from_dict builds."""
    import tempfile
    import wntr
    d = wntr_schema_dict(data, fix_patterns, max_unpacked)
    wn = wntr.network.from_dict(d)
    with tempfile.TemporaryDirectory(prefix="rm_wntr_schema_") as tmp:
        p = os.path.join(tmp, "model.inp")
        wntr.network.write_inpfile(wn, p)
        with open(p, encoding="utf-8") as fh:
            return "".join(ln for ln in fh if not ln.startswith("; Created:"))


# ----------------------------------------------------------------------------- one call
def build_gis(data: bytes, name: str, opts: GisOptions | None = None, max_unpacked: int = ZIP_MAX_TOTAL) -> GisBuild:
    """Any GIS file to a GisBuild: refused in plain words (GisBuild.refusal), a WNTR-schema zip read as its model, or
    its line layer (opts.layer, else the first layer with lines) built into a skeleton.  max_unpacked: read_layers'."""
    opts = opts or GisOptions()
    if data[:2] == b"PK" and wntr_schema_members(data, max_unpacked) is not None:
        b = GisBuild(name=name, kind="wntr_schema", layers=list(WNTR_LAYERS), wntr_schema=True)
        try:
            b.inp_text = wntr_schema_inp(data, True, max_unpacked)
            b.notes.append("A layer in WNTR's own GIS schema with WNTR's JSON file: the model is read from them "
                           "(nodes and links from the layers, options, patterns, curves and controls from the JSON), "
                           "with each junction's demand pattern carried over. This route is for WNTR users: it is "
                           "tested on files laid out as WNTR's to_gis writes them, by our own writer, not on files "
                           "WNTR itself wrote.")
        except Exception as e:  # noqa: BLE001  (from_dict raises KeyError and others on a broken file)
            b.refusal = f"The WNTR GIS files could not be read as a model ({type(e).__name__}: {str(e)[:80]})."
        return b
    try:
        layers, kind = read_layers(data, name, max_unpacked)
    except GisRefused as e:
        return GisBuild(name=name, kind="", layers=[], refusal=str(e))
    except Exception as e:  # noqa: BLE001  (the last guard: a file the readers did not foresee)
        return GisBuild(name=name, kind="", layers=[], refusal=f"The file could not be read as a GIS layer "
                                                             f"({type(e).__name__}). Export the layer again.")
    pick = next((l for l in layers if l.name == opts.layer), None) if opts.layer else None
    if opts.layer and pick is None and any(opts.layer in l.not_read for l in layers):
        return GisBuild(name=name, kind=kind, layers=[l.name for l in layers],
                        refusal=f"The layer {opts.layer} is larger than {MAX_FEATURES:,} features, the most this tool "
                                f"reads from one layer, so it was not read. Export the water mains layer alone (as a "
                                f"GeoPackage or a zipped shapefile), and use that.")
    pick = pick or next((l for l in layers if l.n_lines), layers[0])
    try:
        b = build_from_layer(pick, opts, name, kind, layers)
    except Exception as e:  # noqa: BLE001
        b = GisBuild(name=name, kind=kind, layers=[l.name for l in layers], layer=pick,
                     refusal=f"The layer could not be built into pipes ({type(e).__name__}).")
    return b


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("layer")
    ap.add_argument("--unit", choices=sorted(T.UNIT_M) + ["deg"])
    ap.add_argument("--tolerance-ft", type=float)
    ap.add_argument("--out")
    ap.add_argument("--form", help="an operator form (CSV or Excel): complete the skeleton with it, each pipe's "
                                   "material and the coordinates kept, and write the completed model to --out "
                                   "(task 19)")
    a = ap.parse_args(argv)
    with open(a.layer, "rb") as fh:
        data = fh.read()
    b = build_gis(data, os.path.basename(a.layer), GisOptions(unit=a.unit, tolerance_m=a.tolerance_ft * T.FT
                                                              if a.tolerance_ft else None))
    if b.refusal:
        print(b.refusal)
        return 1
    for q in b.questions:
        print("Question:", q["text"])
    if b.net is not None:
        print(json.dumps({k: v for k, v in b.net.report.items() if not k.endswith("_detail")}, indent=1))
    for w in b.warnings:
        print("Warning:", w)
    if b.inp_text and a.form:          # task 19's review: the form applied to the build, not to its .inp
        from . import operator_form as OF
        return OF.complete_from_files(OF.skeleton_from_build(b, "gis"), a.form, a.out)
    if b.inp_text and a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(b.inp_text)
        print("Skeleton written to", a.out)
    return 0 if b.inp_text else 1


if __name__ == "__main__":
    sys.exit(main())
