"""
cad_read.py: a CAD drawing of water mains (DXF) read into a skeleton model (iteration 5, task 18).

A town engineer's drawing of the water mains is often the only map a small system has.  This builds a skeleton EPANET
model from it, as residualmap/gis_read.py does from a GIS layer, with the same topology engine (residualmap/topology.py):
pipes, their sizes and lengths, and nodes where they meet.  A drawing gives topology and sizes only: the model needs
the operator form (or EPANET) for its source, tanks, pumps, valves, controls, demand, elevations and roughness, and the
readiness report says so on every CAD-built model.  Tested only as converter self-consistency on public models drawn
by our own writers (residualmap/fixtures.py); no real drawing was read.

    python -m residualmap.cad_read DRAWING.dxf [--unit ft|us_ft|m|in|mm|cm|yd] [--tolerance-ft 6] [--out skeleton.inp]
    python -m residualmap.cad_read DRAWING.dxf --form FORM.csv --out completed.inp   # task 19: completed by the
                                    # operator form (operator_form.py), each pipe's material kept
        [--layer-numbers none|in|mm|m] [--default-in 8]

Reading (ezdxf 1.4.4, imported only here and only when a drawing is read, so the .inp path never loads it):
- The file is identified by its content (importers.sniff): an ASCII DXF is read in ezdxf's recover mode (it repairs
  what it can in a damaged file), a binary DXF through ezdxf's binary tag loader and an audit.  A .bak whose content is
  a DXF is read as that DXF.  Since the review: bytes invalid in the drawing's encoding read as U+FFFD (ezdxf's
  default gave lone surrogates, which the page's table could not write); before ezdxf reads an ASCII DXF, a file cut
  short (no '0 EOF') is closed at its last entity's start and warned, and an entity handle drawn twice, or 0, is
  renumbered (ezdxf's audit removed such entities without a word); anything the audit still removes is warned.
  Since the final verification a handle that is no handle (over 16 hex digits, or not hex) is renumbered too and
  sets no new handle's length, and a file that starts with a byte order mark is repaired as well.  Both repairs are
  for ASCII DXF only: a binary DXF is read as it is, and what its audit removes is warned.  A DWG, a MicroStation
  DGN (V7 or V8, an OLE file), an IFC file or a .bak holding a DWG
  gets task 15's guidance, and the free ODA File Converter is offered as a link only (it is closed source and, per
  ezdxf's documentation, cannot run headless on macOS).  A GIS layer is sent to the GIS step, a model to the main page.
- Model space only.  Lines: LINE, LWPOLYLINE and 2D POLYLINE (bulge arcs flattened), 3D POLYLINE, ARC and SPLINE, the
  arcs at most ARC_STEP_DEG degrees a segment (a chord's length is within 0.04% of its arc's), a spline evaluated at a
  bounded number of points; a spline-fit POLYLINE's control-frame vertices (flag 16), which AutoCAD does not draw, are
  skipped (the review).  Points in an entity's own coordinate system (OCS) are taken to world coordinates, so a
  mirrored arc or polyline lands where AutoCAD draws it.  Texts: TEXT, ATTRIB and MTEXT (MTEXT through ezdxf's
  accurate parser, plain_mtext: '\\A1;{\\fArial|b0;8" DIP}' reads as '8" DIP').
- Blocks: an INSERT is exploded with its transform through ezdxf's virtual_entities, a MINSERT's every copy through
  multi_insert, nested blocks to MAX_BLOCK_DEPTH, so mains inside blocks and bound xrefs are read where AutoCAD draws
  them.  AutoCAD's block rules are applied: an entity on layer 0 inside a block takes the INSERT's layer, and BYBLOCK
  linetype and colour take the INSERT's; everything inside an INSERT on a frozen layer, which AutoCAD does not draw, is
  marked so and left out of the build unless asked for (the review).  A MINSERT's row or column count below 1 counts
  as 1.  An external xref (a block that points to another file) is not resolved, with a warning naming it.
- Civil 3D: ACAD_PROXY_ENTITY and AECC objects (a pressure network saved as it is) cannot be read as lines, and get
  "run Export Civil 3D Drawing to AutoCAD (AECTOACAD) first" (critique 15).
- Linear time and bounded reads (the review of task 17 found patterns that took hours on hostile KML and coordinate
  system text): ezdxf's default MTEXT reader ('fast' mode) is quadratic in a run of formatting codes (its timing is in
  outputs/intake/cad_mtext_timing.json, written by python -m residualmap.intake_exp cad_timing), so MTEXT goes through
  the accurate parser, and no text longer than TEXT_MAX_CHARS is parsed at all (a size label is a few dozen
  characters), which bounds every text's cost, and the drawing's notes say how many were not read; layer names are cut
  at topology.SIZE_TEXT_MAX before the size parser sees them.  A read stops, with a plain refusal, at MAX_ENTITIES
  entities visited (exploded blocks included: a block nested in a block ten times over is a few kilobytes and billions
  of lines) and at MAX_POINTS points (a polyline's counted for the points its bulge arcs become, before they are
  flattened); a MINSERT of more than MAX_INSERT_COPIES copies and a spline of more than MAX_SPLINE_POINTS points are
  left out, each with a warning.  The hosted demo reads and builds with tighter limits (residualmap/cloud.py,
  cad_limits), each in a process of its own with a time limit, and refuses a drawing over its line, point and text
  limits before any line is joined.

Building (build_from_drawing):
- Layer styles: every (layer, linetype, colour) combination is listed with its lines and texts, linetype and colour
  resolved as AutoCAD shows them (BYLAYER from the layer), and any can be left out, so abandoned or proposed mains
  drawn dashed on the mains layer can be excluded, not just layers.  Combinations on layers that are off or frozen in
  the drawing are left out by default (the operator can put them back).
- Units from $INSUNITS (the drawing's insertion units).  A unitless drawing (INSUNITS absent or 0), a unit no water
  network is drawn in, or a median line whose length in the stated unit is outside MEDIAN_LINE_FT is asked, and so,
  since the review, are inches whenever the median line read in feet is plausible too (inches are what AutoCAD's
  imperial template starts with, so a drawing made in feet often says inches, and read so it is twelve times too
  short yet plausible).  A stated unit that disagrees with the drawing's measurement system ($MEASUREMENT, imperial or
  metric) is flagged: a drawing made from an imperial template whose units say metres.  The median length cannot catch
  feet read as metres on a network like ky4 (a 501 ft median read as 501 m is still plausible), so without the
  measurement setting that mistake is not flagged.  Both header variables came after R12 ($MEASUREMENT with R14,
  $INSUNITS with R2000), so a drawing saved as R12 DXF is asked its unit.
- Sizes (critique 14 and 15): from the layer name (WATER-8IN, W_8_PVC, WTR-200MM) through task 17's parser
  (topology.parse_diameter, its separators read as spaces; a lone year passed over).  Since the review a bare number
  in a layer name is a size only beside a size word (DIA), a diameter mark or a material (W_8_PVC), and then read
  together with the others by their range; one with none of these (WATER-ZONE-2, SHEET-12) is asked: not sizes, or
  their unit.  Since the final verification the cue must be right beside the number (a material elsewhere in the
  name cued it: PVC-ZONE-3 gave 3 in), a number right after a numbering word (ZONE, PHASE, SHEET) is never cued
  (WATER-ZONE-2-PVC), CL, ST, PE, AC, WS and POLY, common CAD abbreviations, cue nothing, and a cued number is asked
  too when at least DISPUTE_LINES of its layer's lines carry labels of another size and they outnumber those that
  agree.  The operator can also turn layer-name sizes off.  Or a size comes from the nearest label: a text that
  parses as a size with a unit, or a bare number marked with the diameter symbol (%%c, U+2205 and the other marks,
  read as a bare number of the drawing's labels together); other numbers (addresses, stations, lengths) are no labels,
  and nor is a text taller than LABEL_MAX_HEIGHT_FT or LABEL_HEIGHT_REL times the labels' median height (a title or a
  sheet's label, since the review).  Each label goes to its nearest line, by the distance from the line to the text's
  box (a TEXT's insertion point is a corner of its box, or the point its alignment names; an MTEXT's is its attachment
  point), within LABEL_RADIUS_H times the text's height; the pairs of a label and a segment within reach are counted
  before any is measured, and more than max_label_pairs is a plain refusal; since the final verification they are
  measured in memory that does not grow with them (assign_labels).  Labels are assigned to the drawn lines
  before the lines are split at tees, so a line drawn through a tee gives both pieces its label.  A line whose layer
  names a size takes it (a label of another size on it is counted), unless at least half of that layer's lines carry
  a label and most of those labels give another size: then the labels are used on that layer, with a warning.  Then
  sizes propagate along runs of nodes joining two pipes, each unsized pipe taking the size of the nearest sized pipe
  along its run, and what is left gets the commonest size, or the size asked (nothing is built until it is given).
  The report counts the pipes sized from labels, from layer names, by propagation and by the default.
- Hazen-Williams C from a material in the line's label (on a line sized by its layer, a label of the layer's size,
  since the review), else in the layer name (topology.C_TABLE), 130 the labelled default.
- Topology: task 17's engine at 6 ft by default, with the diagnostic table at 0.5, 1, 3, 6 and 12 ft (critique 3c).
- Scale: a drawing with no geographic location (no GEODATA) is warned that mapping taps by GPS needs control taps;
  whether it is drawn to scale is asked (atlas drawings often are not), with a check against known lengths.
- Output: a skeleton .inp titled 'topology and sizes only, needs the operator form', and its readiness report.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import logging
import math
import os
import re
import sys
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import topology as T
from .importers import ASSUMED, DERIVED, DWG_VERSIONS, FROM_FILE, MISSING, Readiness, sniff

CAD_KINDS = ("dxf_ascii", "dxf_binary")
GUIDANCE_KINDS = ("dwg", "dgn_v7", "dgn_v8", "ifc", "ole_other")
ODA_LINK = "https://www.opendesign.com/guestfiles/oda_file_converter"
# $INSUNITS codes a water network can be drawn in: code -> (unit, metres per unit)
INSUNITS = {1: "in", 2: "ft", 4: "mm", 5: "cm", 6: "m", 10: "yd", 21: "us_ft"}
INSUNITS_OTHER = {0: "unitless", 3: "miles", 7: "kilometres", 8: "microinches", 9: "mils", 11: "angstroms",
                  12: "nanometres", 13: "microns", 14: "decimetres", 15: "decametres", 16: "hectometres",
                  17: "gigametres", 18: "astronomical units", 19: "light years", 20: "parsecs",
                  22: "US survey inches", 23: "US survey yards", 24: "US survey miles"}
UNIT_M = {"in": 0.0254, "ft": T.FT, "us_ft": T.FT_US, "yd": 0.9144, "mm": 0.001, "cm": 0.01, "m": 1.0}
UNIT_WORDS = {"in": "inches", "ft": "feet", "us_ft": "US survey feet", "yd": "yards", "mm": "millimetres",
              "cm": "centimetres", "m": "metres"}
IMPERIAL = ("in", "ft", "us_ft", "yd")
MEDIAN_LINE_FT = (20.0, 5000.0)       # a median drawn line outside this, in feet, is implausible for the stated unit
CHAR_W = 0.6                          # a character's width in text heights (an estimate: it only sets a box's length)
LINE_SPACING = 5.0 / 3.0              # MTEXT line spacing, in text heights (AutoCAD's default)
LABEL_RADIUS_H = 3.0                  # a label reaches lines within this many text heights of its box
DIAGNOSTIC_FT = (0.5, 1.0, 3.0, 6.0, 12.0)    # the tolerance table (critique 3c)
ARC_STEP_DEG = 5.0                    # an arc is drawn as chords of at most this many degrees
MAX_ARC_SEGMENTS = 72
SPLINE_SEGMENTS_PER_POINT = 8         # a spline is evaluated at this many points per control point ...
MAX_SPLINE_SEGMENTS = 512             # ... at most
MAX_SPLINE_POINTS = 2_000             # a spline with more control or fit points than this is not read
MAX_BLOCK_DEPTH = 8                   # blocks nested deeper than this are not exploded
MAX_ENTITIES = 2_000_000              # entities visited in one drawing, exploded blocks included
MAX_POINTS = 5_000_000                # line points read from one drawing (as gis_read.WKB_MAX_POINTS)
MAX_INSERT_COPIES = 10_000            # copies of one MINSERT
TEXT_MAX_CHARS = 2_000                # a TEXT or MTEXT longer than this is no size label and is not parsed
LABEL_MAX_HEIGHT_FT = 200.0           # a text taller than this is no pipe label (a title, a sheet label) and sizes nothing
LABEL_HEIGHT_REL = 10.0               # ... nor one taller than this many times the median height of the sized labels
MAX_LABEL_PAIRS = 20_000_000          # (label, line segment) pairs within reach in one drawing, measured before any is
LABEL_CHUNK_PAIRS = 200_000           # pairs measured at once, and segments in one search tree, which bound the
                                      # measure's memory whatever the pairs (the final verification)
DISPUTE_LINES = 2                     # labelled lines of another size that put a layer name's cued number in question
HIDE_WORDS = ("off", "frozen")
TITLE = "CAD skeleton (ResidualMap): topology and sizes only, needs the operator form"
SIZE_SOURCES = ("label", "layer", "propagated", "default")
_DIAMETER_MARKS = ("Ø", "ø", "⌀", "∅")    # U+2205 is the diameter sign AutoCAD's symbol list gives (\U+2205)
_LAYER_SEPARATORS = re.compile(r"[-_.,;|/\\]+")
_LAYER_SIZE_WORDS = ("dia", "diam", "diameter")   # a bare number beside one of these in a layer name is a size
_SURROGATES = re.compile("[\ud800-\udfff]")
LOCAL_LIMITS = {"max_entities": MAX_ENTITIES, "max_points": MAX_POINTS, "max_insert_copies": MAX_INSERT_COPIES,
                "max_label_pairs": MAX_LABEL_PAIRS}
# what the audit of a loaded drawing removed: an entity it could not keep (a handle another entity has, no owner)
_REMOVED_FIXES = ("REMOVED_INVALID_GRAPHIC_ENTITY", "REMOVED_INVALID_DXF_OBJECT",
                  "REMOVED_ENTITY_WITH_INVALID_OWNER_HANDLE", "REMOVED_STANDALONE_ATTRIB_ENTITY")


class CadRefused(ValueError):
    """A drawing that cannot be read or built, in plain words."""


class _Budget(Exception):
    """A read stopped at one of its limits."""


class _LabelsCrowded(Exception):
    """More (label, line segment) pairs within reach than a build measures."""


def _clean(s: str) -> str:
    """A name or a text with any byte its drawing's encoding could not decode (a lone surrogate) as U+FFFD, so that
    every string the reader keeps can be written as UTF-8 (a page's table, a JSON record)."""
    return _SURROGATES.sub("�", s) if s else s

# ----------------------------------------------------------------------------- what is read
@dataclass
class CadLine:
    """One drawn line, flattened: its points in drawing units (world coordinates, x and y)."""
    index: int
    layer: str
    linetype: str
    color: int
    kind: str
    pts: np.ndarray
    block: str | None = None
    frozen_block: str | None = None    # the frozen layer of the INSERT it came from: AutoCAD does not draw it


@dataclass
class CadText:
    """One text: its plain words, its box (4 corners, drawing units) and its height."""
    text: str
    layer: str
    linetype: str
    color: int
    kind: str
    box: np.ndarray
    height: float
    frozen_block: str | None = None


@dataclass
class Drawing:
    name: str
    kind: str = ""
    version: str = ""
    insunits: int | None = None
    measurement: int | None = None
    lines: list = field(default_factory=list)
    texts: list = field(default_factory=list)
    styles: dict = field(default_factory=dict)     # (layer, linetype, colour) -> {'lines', 'texts', 'hidden'}
    hidden_layers: dict = field(default_factory=dict)   # layer -> 'off' or 'frozen'
    counts: dict = field(default_factory=dict)
    proxies: int = 0                  # ACAD_PROXY_ENTITY
    civil3d: int = 0                  # AECC objects
    xrefs: list = field(default_factory=list)
    geodata: bool = False
    notes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    refusal: str | None = None
    sniff_kind: str = ""
    sentence: str = ""


def style_label(key: tuple) -> str:
    """'WATER-MAINS / DASHED / colour 8': a layer style as the page lists it."""
    return f"{key[0]} / {key[1]} / colour {key[2]}"


def guidance(s, data_kind: str | None = None) -> str:
    """What to send instead of a file that is not a DXF (task 15's advice, and the ODA link for a DWG)."""
    from .importers import KINDS
    kind = data_kind or s.kind
    cat = KINDS[kind][0]
    if cat == "gis":
        return (f"{s.sentence()} A GIS layer is built in the step above, 'Build a model from a GIS layer'; this step "
                f"reads CAD drawings saved as DXF.")
    if s.is_model:
        return (f"{s.sentence()} It is a model already: open it on the main page with 'Upload my .inp'; it needs no "
                f"skeleton.")
    out = f"{s.sentence()} {s.advice}"
    if kind == "dwg":
        if s.backup:
            out += (" A .bak file is AutoCAD's backup copy of a drawing: rename it to .dwg, open it in AutoCAD and "
                    "save it as DXF.")
        out += (f" Without AutoCAD, the free ODA File Converter turns a DWG into a DXF: {ODA_LINK} (a link only: "
                f"ResidualMap does not install or run it).")
    return out


# ----------------------------------------------------------------------------- reading
def _get(e, attr: str, default=None):
    try:
        v = e.dxf.get(attr, default)
    except Exception:  # noqa: BLE001  (an unknown entity has no such attribute)
        return default
    return default if v is None else v


def _type(e) -> str:
    try:
        return str(e.dxftype())
    except Exception:  # noqa: BLE001
        return "?"


def _finite(P: np.ndarray) -> bool:
    return P.ndim == 2 and len(P) > 0 and bool(np.isfinite(P).all())


def _arc_points(cx: float, cy: float, r: float, a0: float, sweep: float) -> np.ndarray:
    """Points of an arc of radius r about (cx, cy) from angle a0 (radians) through sweep (signed)."""
    n = int(min(MAX_ARC_SEGMENTS, max(1, math.ceil(abs(sweep) / math.radians(ARC_STEP_DEG)))))
    a = a0 + sweep * np.arange(n + 1) / n
    return np.column_stack([cx + r * np.cos(a), cy + r * np.sin(a)])


def _bulge_points(p1, p2, b: float) -> np.ndarray:
    """The arc a polyline's bulge b draws from p1 to p2 (b = tan of a quarter of its included angle, positive
    counterclockwise), without p1: its interior points and p2."""
    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    c = math.hypot(x2 - x1, y2 - y1)
    if not b or not math.isfinite(b) or c == 0:
        return np.array([[x2, y2]])
    theta = 4.0 * math.atan(b)
    r = c / (2.0 * abs(math.sin(theta / 2.0)))
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    nx, ny = -(y2 - y1) / c, (x2 - x1) / c             # the chord's left normal
    t = math.tan(theta / 2.0)
    k = (c / 2.0) / t if abs(t) > 1e-300 else 0.0
    cx, cy = mx + nx * k, my + ny * k
    a0 = math.atan2(y1 - cy, x1 - cx)
    pts = _arc_points(cx, cy, r, a0, theta)[1:]
    pts[-1] = (x2, y2)                                  # exactly the vertex
    return pts


def _ocs_to_wcs(e, P2: np.ndarray, elevation: float = 0.0) -> np.ndarray:
    """Points in an entity's own coordinate system (with its elevation) in world x and y."""
    try:
        ocs = e.ocs()
        if not ocs.transform:
            return P2
        from ezdxf.math import Vec3
        out = [ocs.to_wcs(Vec3(float(x), float(y), elevation)) for x, y in P2]
        return np.array([(v.x, v.y) for v in out], float)
    except Exception:  # noqa: BLE001
        return P2


def _poly_with_bulges(xyb: list[tuple[float, float, float]], closed: bool) -> np.ndarray:
    pts = [np.array([[xyb[0][0], xyb[0][1]]])]
    seq = list(xyb) + ([xyb[0]] if closed and len(xyb) > 1 else [])
    for (x1, y1, b), (x2, y2, _) in zip(seq, seq[1:]):
        pts.append(_bulge_points((x1, y1), (x2, y2), b) if b else np.array([[x2, y2]]))
    return np.vstack(pts)


def _bulged_points(xyb: list[tuple[float, float, float]], closed: bool) -> int:
    """How many points _poly_with_bulges gives, counted before it runs (each bulge's arc is up to MAX_ARC_SEGMENTS
    chords), so the point budget is charged for what the polyline becomes, not for its vertices (the review).  A
    segment is charged an arc by _bulge_points' own test (a chord of length 0 is straight): two infinite ends compare
    equal, but their chord is nan, and an arc is drawn there (the final verification)."""
    seq = list(xyb) + ([xyb[0]] if closed and len(xyb) > 1 else [])
    step = math.radians(ARC_STEP_DEG)
    n = 1
    for (x1, y1, b), (x2, y2, _) in zip(seq, seq[1:]):
        if b and math.isfinite(b) and math.hypot(x2 - x1, y2 - y1) != 0:
            n += int(min(MAX_ARC_SEGMENTS, max(1, math.ceil(abs(4.0 * math.atan(b)) / step))))
        else:
            n += 1
    return n


def _vertices(e) -> list:
    """A heavy POLYLINE's vertices AutoCAD draws: a spline-fit polyline (PEDIT Spline) keeps its control frame's
    vertices (flag 16), which AutoCAD does not draw, beside the fitted ones (the review)."""
    return [v for v in e.vertices if not (int(_get(v, "flags", 0) or 0) & 16)]


def _unbulged_if_not_finite(xyb: list) -> np.ndarray | None:
    """A polyline's vertices as they are when one is not finite (the walk counts it and leaves it out), so its bulges
    are never flattened (the final verification: 160,000 vertices at infinity made 5.8 million points), else None."""
    if all(math.isfinite(x) and math.isfinite(y) for x, y, _ in xyb):
        return None
    return np.array([(x, y) for x, y, _ in xyb], float)


def _flatten(e, t: str, ctx) -> np.ndarray | None:
    """An entity's line in world x and y, or None (a kind or a shape that is not read)."""
    if t == "LINE":
        s, f = _get(e, "start"), _get(e, "end")
        if s is None or f is None:
            return None
        return np.array([[s[0], s[1]], [f[0], f[1]]], float)
    if t == "LWPOLYLINE":
        xyb = [(float(x), float(y), float(b)) for x, y, b in e.get_points("xyb")]
        if len(xyb) < 2:
            return None
        bad = _unbulged_if_not_finite(xyb)
        if bad is not None:
            return bad
        ctx.budget_points(_bulged_points(xyb, bool(e.closed)))
        P = _poly_with_bulges(xyb, bool(e.closed))
        return _ocs_to_wcs(e, P, float(_get(e, "elevation", 0.0) or 0.0))
    if t == "POLYLINE":
        if e.is_2d_polyline:
            xyb = []
            for v in _vertices(e):
                loc = v.dxf.location
                xyb.append((float(loc[0]), float(loc[1]), float(_get(v, "bulge", 0.0) or 0.0)))
            if len(xyb) < 2:
                return None
            bad = _unbulged_if_not_finite(xyb)
            if bad is not None:
                return bad
            ctx.budget_points(_bulged_points(xyb, bool(e.is_closed)))
            P = _poly_with_bulges(xyb, bool(e.is_closed))
            elev = _get(e, "elevation")
            return _ocs_to_wcs(e, P, float(elev[2]) if elev is not None and len(elev) > 2 else 0.0)
        if e.is_3d_polyline:
            P = np.array([(float(v.dxf.location[0]), float(v.dxf.location[1])) for v in _vertices(e)], float)
            if len(P) < 2:
                return None
            ctx.budget_points(len(P))
            return np.vstack([P, P[:1]]) if e.is_closed else P
        ctx.not_read["POLYLINE (mesh or polyface)"] += 1
        return None
    if t == "ARC":
        c, r = _get(e, "center"), float(_get(e, "radius", 0.0) or 0.0)
        if c is None or not (r > 0 and math.isfinite(r)):
            return None
        a0, a1 = float(_get(e, "start_angle", 0.0)), float(_get(e, "end_angle", 360.0))
        sweep = (a1 - a0) % 360.0 or 360.0
        P = _arc_points(float(c[0]), float(c[1]), r, math.radians(a0), math.radians(sweep))
        return _ocs_to_wcs(e, P, float(c[2]) if len(c) > 2 else 0.0)
    if t == "SPLINE":
        n_pts = max(len(e.control_points), len(e.fit_points))
        if n_pts > MAX_SPLINE_POINTS:
            ctx.counts["splines_too_long"] += 1
            return None
        ct = e.construction_tool()
        segs = int(min(MAX_SPLINE_SEGMENTS, max(8, SPLINE_SEGMENTS_PER_POINT * n_pts)))
        P = np.array([(v.x, v.y) for v in ct.approximate(segs)], float)
        ctx.budget_points(len(P))
        return P
    return None


def _text_box(x0: float, y0: float, w: float, h: float, rot: float, ref) -> np.ndarray:
    c, s = math.cos(rot), math.sin(rot)
    local = np.array([[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]])
    return np.column_stack([ref[0] + local[:, 0] * c - local[:, 1] * s, ref[1] + local[:, 0] * s + local[:, 1] * c])


def _text_entity(e, t: str) -> tuple[str, np.ndarray, float] | None:
    """(plain text, box in world x and y, height) of a TEXT, ATTRIB or MTEXT, or None.  Raw text over TEXT_MAX_CHARS
    is not parsed."""
    if t == "MTEXT":
        raw = e.text
        if len(raw) > TEXT_MAX_CHARS:
            return ("", None, -1.0)
        from ezdxf.tools.text import plain_mtext
        plain = plain_mtext(raw)
        h = float(_get(e, "char_height", 0.0) or 0.0)
        ins = _get(e, "insert")
        if not (h > 0 and math.isfinite(h)) or ins is None:
            return None
        rows = plain.split("\n") or [""]
        w = max(len(r) for r in rows) * h * CHAR_W
        H = h + (len(rows) - 1) * h * LINE_SPACING * float(_get(e, "line_spacing_factor", 1.0) or 1.0)
        ap = int(_get(e, "attachment_point", 1) or 1)
        col, row = (ap - 1) % 3, (ap - 1) // 3
        x0 = (0.0, -w / 2, -w)[col]
        y0 = (-H, -H / 2, 0.0)[row] if 0 <= row <= 2 else 0.0
        try:
            rot = math.radians(float(e.get_rotation()))
        except Exception:  # noqa: BLE001
            rot = math.radians(float(_get(e, "rotation", 0.0) or 0.0))
        return plain, _text_box(x0, y0, w, h if len(rows) == 1 else H, rot, (float(ins[0]), float(ins[1]))), h
    raw = str(_get(e, "text", "") or "")
    if len(raw) > TEXT_MAX_CHARS:
        return ("", None, -1.0)
    plain = e.plain_text()
    h = float(_get(e, "height", 0.0) or 0.0)
    ins = _get(e, "insert")
    if not (h > 0 and math.isfinite(h)) or ins is None:
        return None
    rot = math.radians(float(_get(e, "rotation", 0.0) or 0.0))
    wf = float(_get(e, "width", 1.0) or 1.0)
    halign, valign = int(_get(e, "halign", 0) or 0), int(_get(e, "valign", 0) or 0)
    al = _get(e, "align_point")
    w = len(plain) * h * CHAR_W * wf
    if halign in (3, 5) and al is not None:          # aligned or fit: between the two points
        dx, dy = float(al[0]) - float(ins[0]), float(al[1]) - float(ins[1])
        w, rot, ref, x0, y0 = math.hypot(dx, dy), math.atan2(dy, dx), (float(ins[0]), float(ins[1])), 0.0, 0.0
    else:
        p = ins if (halign == 0 and valign == 0) or al is None else al
        ref = (float(p[0]), float(p[1]))
        x0 = {1: -w / 2, 2: -w, 4: -w / 2}.get(halign, 0.0)
        y0 = -h / 2 if halign == 4 else {2: -h / 2, 3: -h}.get(valign, 0.0)
    box = _text_box(x0, y0, w, h, rot, ref)
    return plain, _ocs_to_wcs(e, box, float(ins[2]) if len(ins) > 2 else 0.0), h


class _Ctx:
    def __init__(self, doc, limits: dict):
        self.doc = doc
        self.limits = {**LOCAL_LIMITS, **(limits or {})}
        self.visited = 0
        self.points = 0
        self.lines: list[CadLine] = []
        self.texts: list[CadText] = []
        self.not_read: Counter = Counter()
        self.counts: Counter = Counter()
        self.by_type: Counter = Counter()
        self.proxies = 0
        self.civil3d = 0
        self.xrefs: dict[str, str] = {}
        self.layers = {}          # lower-case name -> (name, linetype, colour, 'off', 'frozen' or None): AutoCAD's layer
        for lay in doc.layers:    # names are case-insensitive, so every lookup is by the lower-case name
            try:
                name = _clean(str(lay.dxf.name))
                state = "frozen" if lay.is_frozen() else ("off" if lay.is_off() else None)
                self.layers[name.lower()] = (name, _clean(str(lay.dxf.get("linetype", "Continuous") or "Continuous")),
                                             abs(int(lay.dxf.get("color", 7) or 7)), state)
            except Exception:  # noqa: BLE001
                continue
        self.frozen_inserts: Counter = Counter()      # frozen layer of an INSERT -> lines and texts inside it

    def budget_points(self, n: int) -> None:
        self.points += int(n)
        if self.points > self.limits["max_points"]:
            raise _Budget("points")

    def attrs(self, e, parent) -> tuple[str, str, int]:
        """(layer, linetype, colour) as AutoCAD shows the entity: layer 0 and BYBLOCK inside a block take the INSERT's,
        BYLAYER the layer's."""
        layer = _clean(str(_get(e, "layer", "0") or "0"))
        lt = _clean(str(_get(e, "linetype", "BYLAYER") or "BYLAYER"))
        col = int(_get(e, "color", 256) if _get(e, "color", 256) is not None else 256)
        if parent is not None:
            if layer == "0":
                layer = parent[0]
            if lt.upper() == "BYBLOCK":
                lt = parent[1]
            if col == 0:
                col = parent[2]
        lay = self.layers.get(layer.lower())
        if lt.upper() == "BYLAYER":
            lt = lay[1] if lay else "Continuous"
        if col == 256:
            col = lay[2] if lay else 7
        return layer, lt, int(col)

    def hidden_state(self, layer: str) -> str | None:
        """'off', 'frozen' or None for a layer, by its name in any case."""
        lay = self.layers.get(layer.lower())
        return lay[3] if lay else None

    def walk(self, entities, parent=None, depth: int = 0, block: str | None = None, frozen: str | None = None) -> None:
        for e in entities:
            self.visited += 1
            if self.visited > self.limits["max_entities"]:
                raise _Budget("entities")
            t = _type(e)
            try:
                if t == "INSERT":
                    self.explode(e, parent, depth, frozen)
                elif t in ("LINE", "LWPOLYLINE", "POLYLINE", "ARC", "SPLINE"):
                    P = _flatten(e, t, self)
                    if P is None:
                        continue
                    if not _finite(P):
                        self.counts["lines_not_finite"] += 1
                        continue
                    if t in ("LINE", "ARC"):
                        self.budget_points(len(P))
                    layer, lt, col = self.attrs(e, parent)
                    self.lines.append(CadLine(len(self.lines), layer, lt, col, t, P, block, frozen))
                    self.by_type[t] += 1
                    if block is not None:
                        self.counts["lines_from_blocks"] += 1
                    if frozen is not None:
                        self.frozen_inserts[frozen] += 1
                elif t in ("TEXT", "ATTRIB", "MTEXT"):
                    self.text(e, t, parent, frozen)
                elif t == "ACAD_PROXY_ENTITY":
                    self.proxies += 1
                elif t.upper().startswith("AECC"):
                    self.civil3d += 1
                else:
                    self.not_read[t] += 1
            except _Budget:
                raise
            except Exception:  # noqa: BLE001  (one entity ezdxf could not give: counted, not fatal)
                self.counts["entities_unreadable"] += 1

    def text(self, e, t: str, parent, frozen: str | None = None) -> None:
        got = _text_entity(e, t)
        if got is None:
            return
        plain, box, h = got
        if box is None:
            self.counts["texts_too_long"] += 1
            return
        if not _finite(box):
            return
        layer, lt, col = self.attrs(e, parent)
        self.texts.append(CadText(_clean(plain), layer, lt, col, t, box, h, frozen))
        if frozen is not None:
            self.frozen_inserts[frozen] += 1

    def explode(self, e, parent, depth: int, frozen: str | None = None) -> None:
        raw = str(_get(e, "name", "") or "")
        blk = self.doc.blocks.get(raw) if raw else None
        name = _clean(raw)
        if blk is None:
            self.counts["inserts_of_missing_blocks"] += 1
            return
        if blk.block.is_xref or blk.block.is_xref_overlay:
            self.xrefs[name] = _clean(str(blk.block.dxf.get("xref_path", "") or ""))
            return
        if depth >= MAX_BLOCK_DEPTH:
            self.counts["blocks_too_deep"] += 1
            return
        eff = self.attrs(e, parent)
        if frozen is None and self.hidden_state(eff[0]) == "frozen":
            frozen = eff[0]             # an INSERT on a frozen layer: AutoCAD draws nothing of it (the review)
        copies = _copies(e)
        if copies > self.limits["max_insert_copies"]:
            self.counts["minserts_too_large"] += 1
            return
        self.counts["inserts_exploded"] += 1
        for ins in (e.multi_insert() if copies > 1 else (e,)):
            self.walk(ins.virtual_entities(), eff, depth + 1, name, frozen)
        if parent is None:
            for a in getattr(e, "attribs", []) or []:
                self.visited += 1
                self.text(a, "ATTRIB", None, frozen)


def _copies(e) -> int:
    """The copies a (M)INSERT draws, as ezdxf's mcount counts them, with a row or column count below 1 read as 1 (and
    set so in the entity, so multi_insert agrees): a negative count gave a negative mcount, which the block-bomb count
    took as negative copies while the walk exploded the insert once (the review)."""
    try:
        dxf = e.dxf
        rows, cols = int(dxf.get("row_count", 1) or 1), int(dxf.get("column_count", 1) or 1)
        if rows < 1 or cols < 1:
            dxf.row_count, dxf.column_count = max(1, rows), max(1, cols)
        return max(1, int(e.mcount))
    except Exception:  # noqa: BLE001  (an INSERT ezdxf gives no counts for: one copy)
        return 1


def exploded_size(doc, msp, limits: dict) -> int:
    """The entities a walk of model space visits once its blocks are exploded (each entity once, an INSERT's block
    once per copy, to MAX_BLOCK_DEPTH), counted over the block definitions in one pass each per depth, so a block bomb
    (a few kilobytes that explode to billions of lines) is refused before anything is exploded (_Budget)."""
    cap = limits["max_entities"]
    memo: dict = {}

    def weight(entities, depth: int) -> int:
        total = 0
        for e in entities:
            total += 1
            if _type(e) == "INSERT" and depth < MAX_BLOCK_DEPTH:
                name = str(_get(e, "name", "") or "")
                blk = doc.blocks.get(name) if name else None
                copies = _copies(e)
                if blk is not None and not (blk.block.is_xref or blk.block.is_xref_overlay) and \
                        copies <= limits["max_insert_copies"]:
                    key = (name.lower(), depth + 1)
                    if key not in memo:
                        memo[key] = cap + 1         # a block that inserts itself counts as over the limit
                        memo[key] = weight(blk, depth + 1)
                    total += copies * memo[key]
            if total > cap:
                raise _Budget("entities")
        return total

    return weight(msp, 0)


_END_JUNK = b" \t\r\n\x00\x1a"
_KEEP_HANDLES = (b"HEADER", b"CLASSES", b"TABLES", b"OBJECTS")
_BOM = b"\xef\xbb\xbf"
_HANDLE = re.compile(rb"[0-9A-F]{1,16}")    # a DXF handle is a 64-bit number: 1 to 16 hex digits
HANDLE_MAX = 2 ** 64 - 1


def _new_handles(n: int, top: int, taken: set) -> list[int]:
    """n fresh handle numbers: counting up from above the largest valid handle (top), or, when that would pass a 64-bit
    handle, the smallest numbers no valid handle has (the final verification: each fresh handle is at most 16 hex
    digits, so renumbering cannot make the file much larger, whatever handles it holds)."""
    if top + n <= HANDLE_MAX:
        return [top + 1 + k for k in range(n)]
    out, h = [], 1
    while len(out) < n:
        if h not in taken:
            out.append(h)
        h += 1
    return out


def _repair_ascii(data: bytes) -> tuple[bytes, dict]:
    """Two repairs of an ASCII DXF's group codes before ezdxf reads it (the review), each in one pass over its lines,
    and the bytes unchanged when neither is needed (ASCII DXF only: a binary DXF is read as it is, and whatever its
    audit removes is warned):
    - a file cut short (no '0 EOF' at its end, as after an interrupted copy or upload; ezdxf's recover mode drops an
      unended ENTITIES section whole): it is cut at its last entity's start, which may be incomplete, and the
      sequence, block and section open there are closed, so what came before is read;
    - an entity handle (group code 5) that an earlier object already has, the handle 0, or one that is no handle
      (not 1 to 16 hex digits; the final verification), in the ENTITIES or BLOCKS section (ezdxf's audit removes an
      entity whose handle another has, so a writer that reuses handles loses mains without a word): each is given a
      fresh handle (_new_handles) above every valid handle in the file (handles are optional in R12, and a reference
      to a handle drawn twice is ambiguous anyway).  A handle longer than 16 digits sets no 'top', so one very long
      handle cannot make every fresh handle as long (it made a 0.4 MB file 300 MB).
    A UTF-8 byte order mark before the first group code is kept, and the repairs run on what follows it (the final
    verification: with the mark, neither repair ran).
    Returns (bytes, {'cut_short', 'handles_renumbered', 'entity_left_out'})."""
    info = {"cut_short": False, "handles_renumbered": 0, "entity_left_out": False}
    bom = _BOM if data.startswith(_BOM) else b""
    body = data[len(bom):]
    repaired, info = _repair_body(body, info)
    return (data if repaired is body else bom + repaired), info


def _repair_body(data: bytes, info: dict) -> tuple[bytes, dict]:
    """_repair_ascii on a DXF's bytes from its first group code: (the bytes, the very object when nothing is repaired;
    info)."""
    sep = b"\n" if b"\n" in data else b"\r"
    lines = data.split(sep)
    trailer = None
    if lines and lines[-1].strip(_END_JUNK) == b"":
        trailer = lines.pop()          # what follows the last line end (nothing, or padding)
    if len(lines) < 2:
        return data, info
    codes = [c.strip() for c in lines[0::2]]
    try:
        for c in codes:
            int(c)
    except ValueError:              # not one group code a line: nothing is repaired, ezdxf's recover mode decides
        return data, info
    n_pairs = len(lines) // 2
    ends = False
    section, sec_open, seq_open, blk_open = None, False, False, False
    states = {}                     # pair index of each code 0 -> (section, section open, sequence open, block open)
    seen, taken, renumber, top = set(), set(), [], 0
    expect_name = False
    for i in range(n_pairs):
        c = codes[i]
        if c == b"0":
            v = lines[2 * i + 1].strip(_END_JUNK).upper()
            states[i] = (section, sec_open, seq_open, blk_open)
            if v == b"EOF":
                ends = True
                break
            if v == b"SECTION":
                sec_open, expect_name, section = True, True, None
            elif v == b"ENDSEC":
                sec_open, section, seq_open, blk_open = False, None, False, False
            elif v == b"BLOCK":
                blk_open, seq_open = True, False
            elif v == b"ENDBLK":
                blk_open, seq_open = False, False
            elif v == b"SEQEND":
                seq_open = False
            elif v in (b"VERTEX", b"ATTRIB"):
                seq_open = True
            else:
                seq_open = v == b"POLYLINE"
        elif c == b"2" and expect_name:
            section, expect_name = lines[2 * i + 1].strip(_END_JUNK).upper(), False
        elif c in (b"5", b"105"):
            h = lines[2 * i + 1].strip(_END_JUNK).upper()
            v = int(h, 16) if _HANDLE.fullmatch(h) else None
            if v is not None:
                top = max(top, v)
                taken.add(v)
            if c == b"5" and section not in _KEEP_HANDLES and (h in seen or not v):
                renumber.append(i)
            seen.add(h)
    if not ends and states:
        last = max(states)
        section_at, sec_at, seq_at, blk_at = states[last]
        eol = b"\r\n" if lines[0].endswith(b"\r") else b"\n"
        close = []
        if seq_at:
            close.append(b"SEQEND")
        if blk_at and section_at == b"BLOCKS":
            close.append(b"ENDBLK")
        if sec_at:
            close.append(b"ENDSEC")
        tail = b"".join(b"  0" + eol + w + eol for w in close + [b"EOF"])
        info["cut_short"] = True
        v = lines[2 * last + 1].strip(_END_JUNK).upper()
        info["entity_left_out"] = v not in (b"ENDSEC", b"EOF", b"SECTION", b"ENDBLK", b"SEQEND")
        lines = lines[:2 * last]
        renumber = [i for i in renumber if i < last]
    else:
        tail = None
    for i, new in zip(renumber, _new_handles(len(renumber), top, taken)):
        old = lines[2 * i + 1]
        cr = b"\r" if old.endswith(b"\r") else b""
        lines[2 * i + 1] = f"{new:X}".encode() + cr
    info["handles_renumbered"] = len(renumber)
    if tail is None and not renumber:
        return data, info
    if tail is not None:
        return (sep.join(lines) + sep if lines else b"") + tail, info
    return sep.join(lines + ([trailer] if trailer is not None else [])), info


def _load(data: bytes, kind: str):
    """(ezdxf document, audit errors, audit fixes, what the audit removed, the ASCII repairs) of DXF bytes.  Bytes
    that are not valid in the drawing's encoding are read as U+FFFD (errors='replace'), not as lone surrogates, which
    no page or record can write (the review)."""
    logging.getLogger("ezdxf").setLevel(logging.ERROR)
    repairs = {"cut_short": False, "handles_renumbered": 0, "entity_left_out": False}
    if kind == "dxf_binary":
        from ezdxf.document import Drawing as Doc
        from ezdxf.lldxf.tagger import binary_tags_loader
        doc = Doc.load(binary_tags_loader(data, errors="replace"))
        aud = doc.audit()
    else:
        from ezdxf import recover
        data, repairs = _repair_ascii(data)
        doc, aud = recover.read(io.BytesIO(data), errors="replace")
    removed = sum(1 for f in aud.fixes if getattr(getattr(f, "code", None), "name", "") in _REMOVED_FIXES)
    return doc, len(aud.errors), len(aud.fixes), removed, repairs


def read_drawing(data: bytes, name: str = "drawing.dxf", limits: dict | None = None) -> Drawing:
    """A DXF's lines, texts and layer styles (Drawing), or Drawing.refusal in plain words.  Never raises.  limits:
    max_entities, max_points, max_insert_copies (LOCAL_LIMITS by default; the hosted demo's are tighter)."""
    s = sniff(data, name)
    d = Drawing(name=os.path.basename(name), sniff_kind=s.kind, sentence=s.sentence())
    if s.kind not in CAD_KINDS:
        d.refusal = guidance(s)
        return d
    d.kind = s.kind
    try:
        doc, n_err, n_fix, removed, repairs = _load(data, s.kind)
    except Exception as e:  # noqa: BLE001  (a damaged file ezdxf could not repair)
        d.refusal = (f"The drawing could not be read as a DXF ({type(e).__name__}): it is damaged beyond what the "
                     f"reader repairs. Open it in the CAD program and save it again as DXF.")
        return d
    try:
        return _read_doc(d, doc, n_err, n_fix, limits, removed, repairs)
    except _Budget as b:
        lim = {**LOCAL_LIMITS, **(limits or {})}
        what = {"entities": f"more than {lim['max_entities']:,} entities once its blocks are exploded",
                "points": f"more than {lim['max_points']:,} line points"}[str(b)]
        d.lines, d.texts = [], []
        d.refusal = (f"The drawing holds {what}, more than this reader takes. Save a drawing of the water mains alone "
                     f"(their layers, with blocks exploded) and use that.")
        return d
    except Exception as e:  # noqa: BLE001  (the last guard: a drawing the reader did not foresee)
        d.lines, d.texts = [], []
        d.refusal = f"The drawing could not be read ({type(e).__name__}). Save it again as DXF."
        return d


def _read_doc(d: Drawing, doc, n_err: int, n_fix: int, limits: dict | None, removed: int = 0,
              repairs: dict | None = None) -> Drawing:
    repairs = repairs or {}
    d.version = str(doc.dxfversion)
    d.counts["version_name"] = DWG_VERSIONS.get(d.version, d.version)
    d.counts["audit_errors"], d.counts["audit_fixes"] = n_err, n_fix
    if removed:
        d.counts["audit_removed"] = removed
    for k in ("cut_short", "handles_renumbered"):
        if repairs.get(k):
            d.counts[k] = repairs[k]
    for var, attr in (("$INSUNITS", "insunits"), ("$MEASUREMENT", "measurement")):
        try:
            v = doc.header.get(var)
            setattr(d, attr, int(v) if v is not None else None)
        except (TypeError, ValueError):
            setattr(d, attr, None)
    ctx = _Ctx(doc, limits)
    msp = doc.modelspace()
    ctx.counts["entities_once_exploded"] = exploded_size(doc, msp, ctx.limits)
    ctx.walk(msp)
    d.lines, d.texts = ctx.lines, ctx.texts
    d.proxies, d.civil3d, d.xrefs = ctx.proxies, ctx.civil3d, sorted(ctx.xrefs)
    d.hidden_layers = {v[0]: v[3] for v in ctx.layers.values() if v[3]}
    styles: dict = {}
    for ln in d.lines:
        st = styles.setdefault((ln.layer, ln.linetype, ln.color), {"lines": 0, "texts": 0})
        st["lines"] += 1
        if ln.frozen_block is not None:
            st["in_frozen_blocks"] = st.get("in_frozen_blocks", 0) + 1
    for tx in d.texts:
        st = styles.setdefault((tx.layer, tx.linetype, tx.color), {"lines": 0, "texts": 0})
        st["texts"] += 1
        if tx.frozen_block is not None:
            st["in_frozen_blocks"] = st.get("in_frozen_blocks", 0) + 1
    for (layer, _, _), st in styles.items():
        st["hidden"] = ctx.hidden_state(layer)      # AutoCAD's layer names are case-insensitive (the review)
    d.styles = dict(sorted(styles.items()))
    if ctx.frozen_inserts:
        d.counts["in_frozen_blocks"] = dict(sorted(ctx.frozen_inserts.items()))
    paper = 0
    with contextlib.suppress(Exception):
        for lay in doc.layouts:
            if not lay.is_modelspace:
                paper += len(lay)
    try:
        d.geodata = msp.get_geodata() is not None
    except Exception:  # noqa: BLE001  (an R12 drawing has none)
        d.geodata = False
    d.counts.update({"lines": len(d.lines), "texts": len(d.texts), "by_type": dict(sorted(ctx.by_type.items())),
                     "not_read": dict(sorted(ctx.not_read.items())), "paper_space_entities": paper,
                     "entities_visited": ctx.visited, "points": ctx.points,
                     **{k: int(v) for k, v in sorted(ctx.counts.items())}})
    if d.proxies or d.civil3d:
        d.warnings.append(aectoacad(d.proxies + d.civil3d))
    if d.xrefs:
        d.warnings.append(f"{len(d.xrefs)} external reference{'s' if len(d.xrefs) != 1 else ''} ("
                          + ", ".join(f"{k}" + (f" from {v}" if v else "") for k, v in
                                      sorted(ctx.xrefs.items())[:5]) +
                          f") point{'s' if len(d.xrefs) == 1 else ''} to other drawings, which are not read: bind "
                          f"{'them' if len(d.xrefs) != 1 else 'it'} in AutoCAD (XREF, Bind) and save again, or send "
                          f"{'those drawings' if len(d.xrefs) != 1 else 'that drawing'} too.")
    if ctx.counts.get("blocks_too_deep"):
        d.warnings.append(f"{ctx.counts['blocks_too_deep']} block{'s' if ctx.counts['blocks_too_deep'] != 1 else ''} "
                          f"nested more than {MAX_BLOCK_DEPTH} deep {'were' if ctx.counts['blocks_too_deep'] != 1 else 'was'} "
                          f"not exploded.")
    if ctx.counts.get("minserts_too_large"):
        d.warnings.append(f"{ctx.counts['minserts_too_large']} block array{'s' if ctx.counts['minserts_too_large'] != 1 else ''} "
                          f"(MINSERT) of more than {ctx.limits['max_insert_copies']:,} copies "
                          f"{'were' if ctx.counts['minserts_too_large'] != 1 else 'was'} not exploded.")
    if ctx.counts.get("splines_too_long"):
        d.warnings.append(f"{ctx.counts['splines_too_long']} spline{'s' if ctx.counts['splines_too_long'] != 1 else ''} "
                          f"with more than {MAX_SPLINE_POINTS:,} points {'were' if ctx.counts['splines_too_long'] != 1 else 'was'} "
                          f"not read.")
    if removed:
        d.warnings.append(f"The DXF repair removed {removed:,} entit{'ies' if removed != 1 else 'y'} it could not keep (an "
                          f"entity with no valid owner, or one whose handle another entity has), so lines may be "
                          f"missing: open the drawing in the CAD program, run AUDIT, and save it again as DXF.")
    if repairs.get("handles_renumbered"):
        n = repairs["handles_renumbered"]
        d.notes.append(f"{n:,} entit{'ies' if n != 1 else 'y'} had a handle (the drawing's ID of an entity) that another "
                       f"entity has, the handle 0, or one that is no handle (not 1 to 16 hex digits), and "
                       f"{'were' if n != 1 else 'was'} given new ones, so none is lost; "
                       f"a program that writes DXF this way may make other mistakes.")
    if repairs.get("cut_short"):
        d.warnings.append("The file looks cut short: it has no end-of-file marker (0 EOF), as after an interrupted copy "
                          "or upload. It was read up to where it stops"
                          + (", leaving out the last entity, which may be incomplete" if repairs.get("entity_left_out")
                             else "") + "; send the whole file again if lines are missing.")
    if ctx.counts.get("entities_unreadable"):
        d.warnings.append(f"{ctx.counts['entities_unreadable']} entit{'ies' if ctx.counts['entities_unreadable'] != 1 else 'y'} "
                          f"could not be read and {'were' if ctx.counts['entities_unreadable'] != 1 else 'was'} left out.")
    if ctx.counts.get("texts_too_long"):
        n = ctx.counts["texts_too_long"]
        d.notes.append(f"{n:,} text{'s' if n != 1 else ''} longer than {TEXT_MAX_CHARS:,} characters "
                       f"{'were' if n != 1 else 'was'} not read (a pipe label is a few dozen characters): if one is a "
                       f"label, its pipe gets its size another way.")
    if paper:
        d.notes.append(f"Only model space is read: the {paper:,} entities of the drawing's layouts (paper space) are not.")
    nr = {k: v for k, v in ctx.not_read.items() if k not in ("POINT", "DIMENSION", "HATCH", "VIEWPORT", "SOLID")}
    if nr:
        d.notes.append("Not read as mains: " + ", ".join(f"{v:,} {k}" for k, v in sorted(nr.items())) + ".")
    if not d.lines:
        d.refusal = ((aectoacad(d.proxies + d.civil3d) + " The drawing has no plain line to build pipes from.")
                     if d.proxies or d.civil3d else
                     "The file looks cut short (it has no end-of-file marker, 0 EOF), and what it holds has no line in "
                     "model space: send the whole file again." if repairs.get("cut_short") else
                     "The drawing's model space has no line, polyline, arc or spline, so it has no pipes. Save the "
                     "drawing of the water mains as DXF, with the mains in model space.")
    return d


def aectoacad(n: int) -> str:
    """The plain message for Civil 3D objects (proxies or AECC objects) that cannot be read as lines (critique 15)."""
    return (f"{n:,} Civil 3D object{'s' if n != 1 else ''} (a pipe network saved as it is, most likely) cannot be read "
            f"as lines: in Civil 3D, run Export Civil 3D Drawing to AutoCAD (the AECTOACAD command) first, then save "
            f"that drawing as DXF.")


# ----------------------------------------------------------------------------- sizes
# material words of topology.C_TABLE that are also common CAD abbreviations (CL a centre line, ST a street, PE, AC, WS,
# POLY a polyline): they cue no number in a layer name (the final verification), though they still give C
CUE_AMBIGUOUS = ("cl", "st", "pe", "ac", "ws", "poly")
# words that number things in a layer name: a number right after one is that zone's, phase's or sheet's number,
# whatever follows it (WATER-ZONE-2-PVC), never a size cued by a material beside it (the final verification)
NUMBERING_WORDS = ("zone", "zn", "pz", "phase", "ph", "sheet", "sht", "page", "pg", "district", "dist", "area",
                   "stage", "sector", "sec", "section", "lot", "block", "blk", "rev", "no", "num", "number", "level",
                   "lvl", "grid", "map", "plate", "tile", "unit")
_CUE_WORDS = frozenset(_LAYER_SIZE_WORDS) | frozenset(
    w.split()[0] for _, _, words in T.C_TABLE for w in words if w.split()[0] not in CUE_AMBIGUOUS)
_MARKS = "".join(_DIAMETER_MARKS)
_CUE_TOKEN = re.compile(rf"([{_MARKS}]*)(\d+)([{_MARKS}]*|[a-z]+)", re.I)


def _cue(word: str) -> bool:
    return bool(word) and (word in _CUE_WORDS or all(ch in _MARKS for ch in word))


def _size_cue(s: str, value: float) -> bool:
    """Whether a layer name's bare number (value) is cued as a size: one of its tokens equal to value has a size word
    (DIA), a diameter mark or a material (topology.C_TABLE: PVC, DIP, CI) right beside it, as its own neighbour
    (W_8_PVC, MAINS-8-DIA) or joined to it (8PVC, Ø8), and it does not follow a numbering word (ZONE, PHASE, SHEET).
    Since the final verification a cue must be next to the number, not anywhere in the name, and CUE_AMBIGUOUS
    (CL, ST, PE, AC, WS, POLY) cues nothing."""
    toks = s.split()
    for i, tok in enumerate(toks):
        m = _CUE_TOKEN.fullmatch(tok)
        if not m or float(m.group(2)) != value:
            continue
        before = toks[i - 1].lower() if i else ""
        if before in NUMBERING_WORDS:
            continue                    # ZONE 2, PHASE 3: that zone's or phase's number, whatever follows it
        if m.group(1) or _cue(m.group(3).lower()):
            return True
        after = toks[i + 1].lower() if i + 1 < len(toks) else ""
        if _cue(before) or _cue(after):
            return True
    return False


def layer_size(layer: str) -> dict:
    """A layer name's size (topology.parse_diameter on the name with its separators read as spaces), or a flag; a lone
    year (1850 to 2030) is no size.  'cued' (the review): the size has a unit (WATER-8IN, WTR-200MM), or its bare
    number stands right beside a size word (DIA), a diameter mark or a material (W_8_PVC; _size_cue, since the final
    verification only beside it, never after a numbering word, and never from CL, ST, PE, AC, WS or POLY).  A bare
    number with none of these (WATER-ZONE-2, SHEET-12, W-DIST-10, WATER-ZONE-2-PVC, PVC-ZONE-3) may be a zone,
    phase, sheet or district number, so the build asks."""
    s = _LAYER_SEPARATORS.sub(" ", layer[:T.SIZE_TEXT_MAX])
    d = T.parse_diameter(s)
    if d["value"] is not None and d["unit"] is None and float(d["value"]).is_integer() and 1850 <= d["value"] <= 2030:
        return {"value": None, "unit": None, "flag": "unreadable", "text": layer, "cued": False}
    return {**d, "cued": d["value"] is not None and (d["unit"] is not None or _size_cue(s, float(d["value"])))}


def label_size(text: str) -> dict:
    """A label's size: parse_diameter with a unit, or a bare number marked with the diameter symbol ('marked': True);
    a bare number with no mark is no label ('flag': 'no unit')."""
    d = T.parse_diameter(text)
    marked = any(m in text for m in _DIAMETER_MARKS)
    d["marked"] = marked
    if d["value"] is not None and d["unit"] is None and not marked:
        return {**d, "value": None, "flag": "no unit"}
    return d


def _box_seg_dist(C: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Distance between each box (k, 4, 2) and its segment A-B (k, 2): 0 when they meet."""
    d = np.full(len(A), np.inf)
    for j in range(4):
        d = np.minimum(d, T._seg_dist(C[:, j], A, B)[0])
        E0, E1 = C[:, j], C[:, (j + 1) % 4]
        d = np.minimum(d, T._seg_dist(A, E0, E1)[0])
        d = np.minimum(d, T._seg_dist(B, E0, E1)[0])

    def orient(P, Q, R):
        return (Q[:, 0] - P[:, 0]) * (R[:, 1] - P[:, 1]) - (Q[:, 1] - P[:, 1]) * (R[:, 0] - P[:, 0])
    meet = np.zeros(len(A), bool)
    for j in range(4):
        E0, E1 = C[:, j], C[:, (j + 1) % 4]
        o1, o2, o3, o4 = orient(A, B, E0), orient(A, B, E1), orient(E0, E1, A), orient(E0, E1, B)
        meet |= (o1 * o2 < 0) & (o3 * o4 < 0)
    sides = np.stack([orient(C[:, j], C[:, (j + 1) % 4], A) for j in range(4)], axis=1)
    inside = (sides >= 0).all(1) | (sides <= 0).all(1)
    d[meet | inside] = 0.0
    return d


def assign_labels(lines_m: list[np.ndarray], boxes_m: list[np.ndarray], radii_m: list[float],
                  max_pairs: int | None = None) -> list[dict]:
    """For each label box: {'line' (its nearest line's position in lines_m, or None), 'distance', 'second'} where the
    distance is from the line to the box, within the label's radius; 'second' is the next nearest other line's
    distance (None if none within the radius).  Bounded (the review): the (label, segment) pairs within reach are
    counted first, and more than max_pairs raise _LabelsCrowded before any is measured.  Since the final verification
    memory does not grow with the pairs (it grew by about 265 bytes a pair, 5.2 GB for 19 million, after the review's
    fix too): the segments are searched in trees of at most LABEL_CHUNK_PAIRS each, the labels are queried in
    batches whose pairs within reach add up to at most LABEL_CHUNK_PAIRS (a label alone in a batch reaches at most
    that many in one tree), and after each batch only each label's nearest and next nearest line is kept, which is
    all this returns.  Ties go to the lower line position, as before."""
    from itertools import chain
    from scipy.spatial import cKDTree
    out = [{"line": None, "distance": None, "second": None} for _ in boxes_m]
    if not boxes_m or not lines_m:
        return out
    A = np.concatenate([P[:-1] for P in lines_m])
    B = np.concatenate([P[1:] for P in lines_m])
    owner = np.concatenate([np.full(len(P) - 1, i, dtype=np.int64) for i, P in enumerate(lines_m)])
    mid, half = (A + B) / 2, np.hypot(*(B - A).T) / 2
    boxes = np.array(boxes_m, float)
    centres = boxes.mean(axis=1)
    radii = np.asarray(radii_m, float)
    reach = radii + np.hypot(*(boxes - centres[:, None, :]).transpose(2, 0, 1)).max(axis=1)
    tiny = max(float(np.max(half)) * 1e-9, 1e-12)
    cls = np.ceil(np.log2(np.maximum(half, tiny))).astype(int)
    trees, total = [], 0
    for c in np.unique(cls):            # a class of segment lengths, its segments in trees of LABEL_CHUNK_PAIRS
        idx_c = np.where(cls == c)[0]
        r = reach + 2.0 ** c
        for s0 in range(0, len(idx_c), LABEL_CHUNK_PAIRS):
            idx = idx_c[s0:s0 + LABEL_CHUNK_PAIRS]
            tree = cKDTree(mid[idx])
            n_hits = np.asarray(tree.query_ball_point(centres, r, return_length=True), np.int64)
            total += int(n_hits.sum())
            trees.append((idx, tree, r, n_hits))
        if max_pairs is not None and total > max_pairs:
            raise _LabelsCrowded(total)
    n_lab = len(boxes)
    d1, o1 = np.full(n_lab, np.inf), np.full(n_lab, -1, np.int64)     # each label's nearest line so far
    d2, o2 = np.full(n_lab, np.inf), np.full(n_lab, -1, np.int64)     # ... and its next nearest other line

    def measure(li: np.ndarray, si: np.ndarray) -> None:
        dist = _box_seg_dist(boxes[li], A[si], B[si])
        keep = dist <= radii[li]
        if not keep.any():
            return
        li, o, dist = li[keep], owner[si[keep]], dist[keep]
        ks = np.unique(li)
        for dd, oo in ((d1, o1), (d2, o2)):         # with the two lines each label had so far
            had = ks[oo[ks] >= 0]
            li, o, dist = np.concatenate([li, had]), np.concatenate([o, oo[had]]), np.concatenate([dist, dd[had]])
        order = np.lexsort((dist, o, li))           # each (label, line)'s nearest segment first
        li, o, dist = li[order], o[order], dist[order]
        first = np.ones(len(li), bool)
        first[1:] = (li[1:] != li[:-1]) | (o[1:] != o[:-1])
        li, o, dist = li[first], o[first], dist[first]
        order = np.lexsort((o, dist, li))           # each label's lines by distance, then position
        li, o, dist = li[order], o[order], dist[order]
        pos = np.arange(len(li))
        new = np.ones(len(li), bool)
        new[1:] = li[1:] != li[:-1]
        rank = pos - np.maximum.accumulate(np.where(new, pos, 0))
        o2[ks], d2[ks] = -1, np.inf
        for want, dd, oo in ((0, d1, o1), (1, d2, o2)):
            at = rank == want
            oo[li[at]], dd[li[at]] = o[at], dist[at]

    for idx, tree, r, n_hits in trees:
        hit = np.nonzero(n_hits)[0]
        cum = np.cumsum(n_hits[hit])
        start = 0
        while start < len(hit):
            base = int(cum[start - 1]) if start else 0
            end = max(start + 1, int(np.searchsorted(cum, base + LABEL_CHUNK_PAIRS, side="right")))
            ks = hit[start:end]
            res = tree.query_ball_point(centres[ks], r[ks])
            lens = np.fromiter((len(h) for h in res), np.int64, len(ks))
            measure(np.repeat(ks, lens), idx[np.fromiter(chain.from_iterable(res), np.int64, int(lens.sum()))])
            start = end
    for k in np.nonzero(o1 >= 0)[0].tolist():
        out[k] = {"line": int(o1[k]), "distance": float(d1[k]), "second": float(d2[k]) if o2[k] >= 0 else None}
    return out


def propagate_sizes(net: T.Network, size_of: dict) -> dict:
    """{pipe ID: size} for pipes with no size, from the nearest sized pipe along runs of nodes that join exactly two
    pipes (a breadth-first search from every sized pipe at once, in pipe ID order)."""
    deg = net.degree()
    at: dict = {}
    for p in net.pipes:
        at.setdefault(p["start"], []).append(p["id"])
        at.setdefault(p["end"], []).append(p["id"])
    ends = {p["id"]: (p["start"], p["end"]) for p in net.pipes}
    got: dict = {}
    queue = sorted(pid for pid in size_of if pid in ends)
    head = 0
    known = dict(size_of)
    while head < len(queue):
        pid = queue[head]
        head += 1
        for n in ends[pid]:
            if deg.get(n) != 2:
                continue
            for q in at[n]:
                if q != pid and q not in known:
                    known[q] = known[pid]
                    got[q] = known[pid]
                    queue.append(q)
    return got


# ----------------------------------------------------------------------------- building
@dataclass
class CadOptions:
    unit: str | None = None            # over the drawing's $INSUNITS: a key of UNIT_M
    tolerance_m: float | None = None
    exclude: list = field(default_factory=list)   # layer styles (layer, linetype, colour) left out
    include_hidden: bool = False       # lines and texts on layers that are off or frozen: left out unless True
    frozen_blocks: bool = False        # lines and texts inside an INSERT on a frozen layer: left out unless True
    label_radius_h: float = LABEL_RADIUS_H
    layer_sizes: bool = True           # False: sizes in layer names are not read (labels, runs and the default size)
    layer_numbers: str | None = None   # 'none' (not sizes), or 'in', 'mm' or 'm', for numbers in layer names with no
                                       # unit, size word or material beside them (WATER-ZONE-2), and since the final
                                       # verification for cued ones the labels dispute: asked when None
    layer_bare_unit: str | None = None   # 'in', 'mm' or 'm' for bare sizes in layer names the range cannot settle
    label_bare_unit: str | None = None   # ... and for labels marked with the diameter symbol
    default_diameter_in: float | None = None
    known_lengths_ft: dict = field(default_factory=dict)
    to_scale: bool | None = None


@dataclass
class CadBuild:
    name: str
    drawing: Drawing | None = None
    questions: list = field(default_factory=list)
    refusal: str | None = None
    net: T.Network | None = None
    tolerance_m: float | None = None
    unit: str | None = None
    unit_how: str = ""
    unit_check: dict = field(default_factory=dict)
    excluded: list = field(default_factory=list)
    lines_used: int = 0
    diagnostics: list = field(default_factory=list)
    diam_in: dict = field(default_factory=dict)
    size_source: dict = field(default_factory=dict)
    size_counts: dict = field(default_factory=dict)
    labels: dict = field(default_factory=dict)
    layer_sizes: dict = field(default_factory=dict)
    c_of: dict = field(default_factory=dict)
    c_counts: dict = field(default_factory=dict)
    c_material: dict = field(default_factory=dict)     # task 19: pipe ID -> the material its label or layer names
                                                       # (topology.material_of: 'DIP' ductile iron), or None; the
                                                       # operator form's pipe_c rows set C by it
    default_size: dict = field(default_factory=dict)   # {'size_in', 'from': 'commonest' or 'given'} when one is used
    scale_check: dict = field(default_factory=dict)
    frame: dict = field(default_factory=dict)
    inp_text: str | None = None
    notes: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.inp_text is not None

    def blocking_questions(self) -> list:
        return [q for q in self.questions if q.get("blocking")]


def _ask(b: CadBuild, key: str, text: str, options=None, blocking: bool = True) -> None:
    b.questions.append({"key": key, "text": text, "options": list(options or []), "blocking": blocking})


def default_excluded(d: Drawing) -> list:
    """The layer styles left out unless the operator puts them back: those on layers that are off or frozen, with
    lines or texts (since the review a frozen layer's old labels size nothing, as AutoCAD does not show them)."""
    return [k for k, v in d.styles.items() if v.get("hidden")]


def build_from_drawing(d: Drawing, opts: CadOptions | None = None, limits: dict | None = None) -> CadBuild:
    """A read drawing built into a skeleton (CadBuild.inp_text; None while a blocking question is open or nothing can
    be built).  Never raises: anything it did not foresee is a plain refusal.  limits: max_label_pairs (LOCAL_LIMITS
    by default; the hosted demo's is tighter, cloud.cad_limits)."""
    try:
        return _build(d, opts or CadOptions(), {**LOCAL_LIMITS, **(limits or {})})
    except Exception as e:  # noqa: BLE001  (the last guard)
        return CadBuild(name=d.name, drawing=d, refusal=f"The drawing could not be built into pipes ({type(e).__name__}). "
                                                        f"Save it again as DXF.")


def _layer_numbers_text(plain: dict, disputed: dict) -> str:
    """The layer_numbers question: layer names whose number has no cue beside it (plain: name -> number), and those
    whose cued number the labels on their lines dispute (disputed: name -> (lines with another size, lines labelled,
    number))."""
    parts = []
    if plain:
        listed = ", ".join(f"{n} ({v:g})" for n, v in list(plain.items())[:6]) + (", and more" if len(plain) > 6 else "")
        parts.append(f"{len(plain):,} layer name{'s hold' if len(plain) != 1 else ' holds'} a number with no unit, size "
                     f"word or material beside it: {listed}.")
    if disputed:
        listed = ", ".join(f"{n} ({v:g}: {k:,} of its {m:,} labelled lines give another size)"
                           for n, (k, m, v) in list(sorted(disputed.items()))[:6]) + (
            ", and more" if len(disputed) > 6 else "")
        parts.append(f"{len(disputed):,} layer name{'s hold' if len(disputed) != 1 else ' holds'} a number with a "
                     f"material or size word beside it that the labels on its lines dispute: {listed}.")
    return " ".join(parts) + " Are these numbers pipe sizes? A zone, phase, sheet or district number is not one."


def _drawing_words(text: str) -> str:
    """The topology engine's sentences, which name a layer, said of a drawing."""
    return text.replace("in the layer if", "in the drawing if").replace("The layer's", "The drawing's")


def _median_ft(lines: list[CadLine], unit_m: float) -> float:
    L = [float(np.hypot(*(ln.pts[1:] - ln.pts[:-1]).T).sum()) for ln in lines]
    return float(np.median(L)) * unit_m / T.FT if L else 0.0


def _build(d: Drawing, opts: CadOptions, limits: dict | None = None) -> CadBuild:
    limits = {**LOCAL_LIMITS, **(limits or {})}
    b = CadBuild(name=d.name, drawing=d)
    b.notes += d.notes
    b.warnings += d.warnings
    if d.refusal:
        b.refusal = d.refusal
        return b
    excluded = {tuple(k) for k in opts.exclude}
    hidden = set(default_excluded(d)) if not opts.include_hidden else set()
    out = excluded | hidden
    b.excluded = sorted(out)
    lines = [ln for ln in d.lines if (ln.layer, ln.linetype, ln.color) not in out
             and (ln.frozen_block is None or opts.frozen_blocks)]
    texts = [tx for tx in d.texts if (tx.layer, tx.linetype, tx.color) not in out
             and (tx.frozen_block is None or opts.frozen_blocks)]
    b.lines_used = len(lines)
    if hidden and not opts.include_hidden:
        n, t = sum(d.styles[k]["lines"] for k in hidden), sum(d.styles[k]["texts"] for k in hidden)
        what = " and ".join(x for x in (f"{n:,} line{'s' if n != 1 else ''}" if n else "",
                                        f"{t:,} text{'s' if t != 1 else ''}" if t else "") if x)
        b.notes.append(f"{what} on layers that are off or frozen in the drawing {'are' if n + t != 1 else 'is'} left "
                       f"out; put them back in the layer styles if they are mains or their labels.")
    fz = d.counts.get("in_frozen_blocks") or {}
    if fz and not opts.frozen_blocks:
        n = sum(fz.values())
        b.notes.append(f"{n:,} line{'s' if n != 1 else ''} and text{'s' if n != 1 else ''} inside blocks inserted on "
                       f"frozen layers ({', '.join(list(fz)[:5])}), which AutoCAD does not draw, "
                       f"{'are' if n != 1 else 'is'} left out; read them too if they are mains in service.")
    if not lines:
        b.refusal = "No line is left on the layer styles chosen: choose the styles that hold the water mains."
        return b
    # units
    declared = INSUNITS.get(d.insunits) if d.insunits is not None else None
    if opts.unit in UNIT_M:
        unit, how = opts.unit, f"chosen by you ({UNIT_WORDS[opts.unit]})"
        if declared and declared != opts.unit:
            b.warnings.append(f"The drawing's unit was set to {UNIT_WORDS[opts.unit]}, over {UNIT_WORDS[declared]} from "
                              f"its insertion units (INSUNITS).")
    elif declared is None:
        said = ("does not say its unit (it gives no insertion units: INSUNITS is absent or 0)" if d.insunits in (None, 0) else
                f"says its unit is {INSUNITS_OTHER.get(d.insunits, f'code {d.insunits}')}, which no water network "
                f"is drawn in")
        med = {u: _median_ft(lines, m) for u, m in (("ft", T.FT), ("m", 1.0))}
        _ask(b, "unit", f"The drawing {said}. In what unit is it drawn? (Its median line is {med['ft']:,.0f} ft long "
                        f"in feet, {med['m']:,.0f} ft in metres.) Lines cannot be joined until the unit is known.",
             list(UNIT_M))
        return b
    else:
        unit, how = declared, f"from the drawing's insertion units (INSUNITS {d.insunits})"
        med = _median_ft(lines, UNIT_M[unit])
        lo, hi = MEDIAN_LINE_FT
        if not lo <= med <= hi:
            _ask(b, "unit", f"The drawing says it is drawn in {UNIT_WORDS[unit]} (INSUNITS {d.insunits}), but then its "
                            f"median line is {med:,.1f} ft long, which is implausible for water mains (between "
                            f"{lo:,.0f} and {hi:,.0f} ft is expected). In what unit is it drawn?", list(UNIT_M))
            b.unit_check = {"median_ft": med, "implausible": True}
            return b
        med_ft = _median_ft(lines, T.FT)
        if unit == "in" and lo <= med_ft <= hi:
            # inches are what AutoCAD's imperial template starts with, so a drawing in feet often carries them too, and
            # read in inches such a drawing is twelve times too short yet still plausible (the review): asked
            _ask(b, "unit", f"The drawing's insertion units say inches (INSUNITS 1), the unit AutoCAD's imperial "
                            f"template starts with, so a drawing made in feet often carries them too. Read in inches "
                            f"its median line is {med:,.1f} ft long; read in feet, {med_ft:,.0f} ft. In what unit is "
                            f"it drawn? Lines cannot be joined until the unit is known.", list(UNIT_M))
            b.unit_check = {"median_ft": med, "median_ft_if_feet": med_ft, "inches_or_feet": True}
            return b
    b.unit, b.unit_how = unit, how
    unit_m = UNIT_M[unit]
    med = _median_ft(lines, unit_m)
    b.unit_check = {**b.unit_check, "unit": unit, "median_ft": med}
    if d.measurement in (0, 1) and opts.unit is None:
        imperial = unit in IMPERIAL
        if imperial != (d.measurement == 0):
            system = "imperial" if d.measurement == 0 else "metric"
            text = (f"The drawing's units say {UNIT_WORDS[unit]}, but its measurement setting (MEASUREMENT) is "
                    f"{system}: one of them is wrong (a drawing started from an {system} template whose units were "
                    f"changed, most often). Read in {UNIT_WORDS[unit]} its median line is {med:,.0f} ft long. Check "
                    f"the unit, and choose it if it is not {UNIT_WORDS[unit]}.")
            _ask(b, "unit_check", text, list(UNIT_M), blocking=False)
            b.warnings.append(text)
            b.unit_check["measurement_disagrees"] = True
    origin = np.floor(np.min([ln.pts.min(axis=0) for ln in lines], axis=0))
    b.frame = {"kind": "projected", "origin": origin.tolist(), "unit_m": unit_m, "unit": unit}
    tlines = [T.Line(ln.index, 0, (ln.pts - origin) * unit_m, ln.pts) for ln in lines]
    tol = opts.tolerance_m or T.DEFAULT_TOL_M["gis"]
    b.tolerance_m = tol
    # sizes of the drawn lines, before any line is split at a tee: layer names, then labels
    names = sorted({ln.layer for ln in lines})
    lsz = {n: layer_size(n) for n in names} if opts.layer_sizes else {}
    if not opts.layer_sizes:
        b.notes.append("Sizes in layer names are not read, as you chose: each pipe's size comes from its label, a run "
                       "of pipes or the default size.")
    bare_l = [v["value"] for v in lsz.values() if v["value"] is not None and v["unit"] is None and v["cued"]]
    plain = {n: v["value"] for n, v in lsz.items() if v["value"] is not None and v["unit"] is None and not v["cued"]}
    lbu = opts.layer_bare_unit or (T.bare_unit(bare_l) if bare_l else None)
    if bare_l and lbu is None:
        _ask(b, "layer_bare_unit", f"Layer names hold sizes with no unit, from {min(bare_l):g} to {max(bare_l):g}, which "
                                   f"fit none of inches (2 to 48), mm (50 to 1,200) or metres (0.05 to 1.2). Which unit "
                                   f"are they in?", ["in", "mm", "m"])
    numbers = opts.layer_numbers if opts.layer_numbers in ("none", "in", "mm", "m") else None
    if plain and numbers is None:
        _ask(b, "layer_numbers", _layer_numbers_text(plain, {}), ["none", "in", "mm", "m"])
    layer_in = {}
    for n, v in lsz.items():
        if v["value"] is None:
            continue
        unit_of = v["unit"] or (lbu if v["cued"] else (numbers if numbers != "none" else None))
        if unit_of:
            layer_in[n] = T.to_inches(v["value"], unit_of)
    b.layer_sizes = {n: {"size_in": layer_in.get(n), "text": n, "flag": lsz[n]["flag"] if n in lsz else None}
                     for n in names}
    parsed = [(tx, label_size(tx.text)) for tx in texts]
    # a text far taller than the drawing's labels is a title or a sheet's label, not a pipe's (the review): its reach
    # (3 heights) would take in every line of the drawing
    hts = np.array([tx.height * unit_m / T.FT for tx, p in parsed if p["value"] is not None], float)
    if len(hts):
        cap = min(LABEL_MAX_HEIGHT_FT, LABEL_HEIGHT_REL * float(np.median(hts)))
        tall = [i for i, (tx, p) in enumerate(parsed) if p["value"] is not None and tx.height * unit_m / T.FT > cap]
        for i in tall:
            parsed[i] = (parsed[i][0], {**parsed[i][1], "value": None, "flag": "too tall"})
        if tall:
            b.warnings.append(f"{len(tall):,} text{'s' if len(tall) != 1 else ''} that read as a pipe size "
                              f"{'are' if len(tall) != 1 else 'is'} taller than {cap:,.3g} ft ({LABEL_MAX_HEIGHT_FT:g} "
                              f"ft, or {LABEL_HEIGHT_REL:g} times the labels' median height of "
                              f"{float(np.median(hts)):,.3g} ft) and {'are' if len(tall) != 1 else 'is'} not used as "
                              f"labels: a title or a sheet's label, most often.")
    bare_t = [p["value"] for _, p in parsed if p["value"] is not None and p["unit"] is None]
    tbu = opts.label_bare_unit or (T.bare_unit(bare_t) if bare_t else None)
    if bare_t and tbu is None:
        _ask(b, "label_bare_unit", f"Labels marked with the diameter symbol hold sizes with no unit, from "
                                   f"{min(bare_t):g} to {max(bare_t):g}, which fit none of inches (2 to 48), mm (50 to "
                                   f"1,200) or metres (0.05 to 1.2). Which unit are they in?", ["in", "mm", "m"])
    sized = [(tx, T.to_inches(p["value"], p["unit"] or tbu)) for tx, p in parsed
             if p["value"] is not None and (p["unit"] or tbu)]
    flags = Counter(p["flag"] for _, p in parsed if p["value"] is None and p["flag"] not in (None, "empty"))
    radius_h = opts.label_radius_h if opts.label_radius_h and opts.label_radius_h > 0 else LABEL_RADIUS_H
    try:
        got = assign_labels([tl.xy for tl in tlines], [(tx.box - origin) * unit_m for tx, _ in sized],
                            [radius_h * tx.height * unit_m for tx, _ in sized], limits.get("max_label_pairs"))
    except _LabelsCrowded as e:
        b.labels = {"pairs_within_reach": int(e.args[0])}
        b.refusal = (f"The drawing's labels reach too many lines: {e.args[0]:,} pairs of a label and a line segment "
                     f"within reach, more than the {limits['max_label_pairs']:,} this reader measures. Its texts are "
                     f"probably far taller than its lines are long: leave out the layer styles of texts that are not "
                     f"pipe sizes, or lower the reach of a label.")
        return b
    on_line: dict = {}
    ambiguous = 0
    for (tx, size), g in zip(sized, got):
        if g["line"] is None:
            continue
        on_line.setdefault(lines[g["line"]].index, []).append((g["distance"], size, tx))
        if g["second"] is not None:
            ambiguous += 1
    line_size, line_src, line_text = {}, {}, {}
    conflicts = disagree = 0
    nearest = {ln.index: sorted(on_line.get(ln.index, []), key=lambda x: (x[0], x[1])) for ln in lines}
    # a size read from a bare number its name cues (W_8_PVC, no unit) that the labels on its lines dispute: at least
    # DISPUTE_LINES of them give another size and they outnumber those that agree, so the number may be no size: it
    # is asked as a number with no cue is (the final verification; one stray label, as on ky4's drawing, is not enough)
    disputed: dict = {}
    for n, size in layer_in.items():
        if lsz[n]["unit"] is not None:
            continue
        labelled = [ln for ln in lines if ln.layer == n and nearest[ln.index]]
        other = [ln for ln in labelled if abs(nearest[ln.index][0][1] - size) > 0.01 * size]
        if len(other) >= DISPUTE_LINES and 2 * len(other) > len(labelled):
            disputed[n] = (len(other), len(labelled), lsz[n]["value"])
    if disputed:
        if numbers is None:
            q = next((q for q in b.questions if q["key"] == "layer_numbers"), None)
            text = _layer_numbers_text(plain, disputed)
            if q is None:
                _ask(b, "layer_numbers", text, ["none", "in", "mm", "m"])
            else:
                q["text"] = text
        for n, (_, _, value) in disputed.items():
            layer_in.pop(n, None)
            if numbers not in (None, "none"):
                layer_in[n] = T.to_inches(value, numbers)
    b.layer_sizes = {n: {"size_in": layer_in.get(n), "text": n, "flag": lsz[n]["flag"] if n in lsz else None}
                     for n in names}
    # a layer whose name gives a size while most of its lines carry a label of another size: the name's number is
    # probably not a size (a zone, a phase), so the labels are used on its lines (the review)
    labels_win: dict = {}
    for n, size in layer_in.items():
        on = [ln for ln in lines if ln.layer == n]
        labelled = [ln for ln in on if nearest[ln.index]]
        other = [ln for ln in labelled if abs(nearest[ln.index][0][1] - size) > 0.01 * size]
        if len(labelled) >= 2 and 2 * len(labelled) >= len(on) and 2 * len(other) > len(labelled):
            labels_win[n] = (len(other), len(labelled), size)
    for ln in lines:
        lab = nearest[ln.index]
        if lab and len({round(s, 3) for _, s, _ in lab}) > 1:
            conflicts += 1
        if ln.layer in layer_in and ln.layer not in labels_win:
            size = layer_in[ln.layer]
            line_size[ln.index], line_src[ln.index] = size, "layer"
            if lab and abs(lab[0][1] - size) > 0.01 * size:
                disagree += 1
            elif lab:                      # a label of the layer's size: its material gives C (the review)
                line_text[ln.index] = lab[0][2].text
        elif lab:
            line_size[ln.index], line_src[ln.index], line_text[ln.index] = lab[0][1], "label", lab[0][2].text
    b.labels = {"texts": len(texts), "sizes_read": len(sized), "assigned": sum(len(v) for v in on_line.values()),
                "not_within_reach": sum(1 for g in got if g["line"] is None), "within_reach_of_two_lines": ambiguous,
                "lines_with_two_sizes": conflicts, "layer_and_label_disagree": disagree,
                "bare_unit": tbu, "radius_text_heights": radius_h,
                "not_sizes": dict(sorted(flags.items()))}
    if conflicts:
        b.warnings.append(f"{conflicts} line{'s' if conflicts != 1 else ''} {'have' if conflicts != 1 else 'has'} two "
                          f"labels of different sizes within reach: the nearest is used. Check them on the drawing.")
    if disagree:
        b.warnings.append(f"{disagree} line{'s' if disagree != 1 else ''} {'are' if disagree != 1 else 'is'} on a layer "
                          f"whose name gives a size, with a label of another size near {'them' if disagree != 1 else 'it'}: "
                          f"the layer's size is used.")
    for n, (k, m, size) in sorted(labels_win.items()):
        b.warnings.append(f"On layer {n}, {k:,} of the {m:,} lines with a label carry another size than its name gives "
                          f"({size:g} in): the labels are used on that layer's lines, since the name's number is "
                          f"probably not a size. Choose not to read sizes from layer names if none of them are.")
    if labels_win:
        b.labels["layers_whose_labels_win"] = sorted(labels_win)
    if disputed:
        b.labels["layers_disputed"] = sorted(disputed)
    if any(q["key"] in ("layer_bare_unit", "label_bare_unit", "layer_numbers") for q in b.questions):
        return b
    # the network
    try:
        net = T.build_network(tlines, tol, "xy", None)
    except T.TopologyRefused as e:
        b.refusal = _drawing_words(str(e)).replace("the layer", "the drawing")
        return b
    r = net.report
    if not net.pipes or r["dropped_pieces"] > 0.5 * r["pieces"]:
        b.refusal = ((f"No pipe can be built from the drawing's {r['lines']:,} lines" if not net.pipes else
                      f"{r['dropped_pieces']:,} of the drawing's {r['pieces']:,} pieces of line would be left out")
                     + f": joining ends within {tol / T.FT:.1f} ft puts both ends of "
                     + ("every one" if not net.pipes else "each") + " in one node. The drawing's unit is probably "
                     "wrong, or the joining distance far too large for it.")
        return b
    b.net = net
    b.warnings += [_drawing_words(f) for f in net.flags]
    by_index = {ln.index: ln for ln in lines}
    size_of, src = {}, {}
    for p in net.pipes:
        f = p["feature"]
        if f in line_size:
            size_of[p["id"]], src[p["id"]] = line_size[f], line_src[f]
    prop = propagate_sizes(net, size_of)
    for pid, v in prop.items():
        size_of[pid], src[pid] = v, "propagated"
    missing = [p["id"] for p in net.pipes if p["id"] not in size_of]
    if missing:
        sizes = [round(v, 2) for v in size_of.values()]
        fill = opts.default_diameter_in
        if fill is None and sizes:
            fill = max(set(sizes), key=lambda s: (sizes.count(s), -s))
        if fill is None:
            _ask(b, "default_diameter_in", "No pipe size can be read from the drawing (no layer name or label gives "
                                           "one). What size are its pipes (inches)? Every pipe gets that size until you "
                                           "give the real ones.", [])
            return b
        for pid in missing:
            size_of[pid], src[pid] = float(fill), "default"
        given = opts.default_diameter_in is not None
        b.default_size = {"size_in": float(fill), "from": "given" if given else "commonest"}
        b.warnings.append(f"{len(missing)} pipe{'s' if len(missing) != 1 else ''} of {len(net.pipes)} "
                          f"{'have' if len(missing) != 1 else 'has'} no size from a layer name, a label or a run of "
                          f"pipes, and {'are' if len(missing) != 1 else 'is'} given {fill:g} in, "
                          f"{'the size you gave' if given else 'the commonest size'}.")
    b.diam_in, b.size_source = size_of, src
    cnt = Counter(src.values())
    b.size_counts = {k: int(cnt.get(k, 0)) for k in SIZE_SOURCES}
    b.size_counts["lines_from_labels"] = sum(1 for v in line_src.values() if v == "label")
    b.size_counts["lines_from_layers"] = sum(1 for v in line_src.values() if v == "layer")
    # Hazen-Williams C from a material in the label or the layer name
    cc = Counter()
    classes: Counter = Counter()
    for p in net.pipes:
        f = p["feature"]
        c, cls, mat = T.material_of(line_text.get(f)) if f in line_text else (T.C_DEFAULT, None, None)
        if cls is None:
            c, cls, mat = T.material_of(_LAYER_SEPARATORS.sub(" ", by_index[f].layer[:T.SIZE_TEXT_MAX]))
        if cls:
            cc["from_material"] += 1
            classes[cls] += 1
        else:
            cc["default"] += 1
        b.c_of[p["id"]] = c
        b.c_material[p["id"]] = mat
    b.c_counts = {"from_material": int(cc["from_material"]), "default": int(cc["default"]),
                  "classes": dict(sorted(classes.items())), "uniform": len(set(b.c_of.values())) <= 1}
    # scale and location
    if not d.geodata:
        b.warnings.append("The drawing has no geographic location: its coordinates are its own, so mapping sample taps "
                          "by GPS position needs two or three control taps whose place on the drawing you know.")
    if opts.to_scale is None:
        _ask(b, "to_scale", "Is the drawing to scale (true lengths)? Atlas and schematic drawings often are not. Give "
                            "two or three pipes whose true length you know to check.", ["yes", "no"], blocking=False)
    elif opts.to_scale is False:
        b.refusal = ("The drawing is not to scale, so its pipes' lengths cannot be read from it. Give the real lengths "
                     "in EPANET after the import, or trace the mains again over a map that is to scale.")
        return b
    if opts.known_lengths_ft:
        lengths = {p["id"]: p["length_m"] / T.FT for p in net.pipes}
        rows = [{"pipe": pid, "drawn_ft": lengths[pid], "true_ft": float(v), "ratio": lengths[pid] / float(v)}
                for pid, v in opts.known_lengths_ft.items() if pid in lengths and v and v > 0]
        medr = float(np.median([x["ratio"] for x in rows])) if rows else None
        b.scale_check = {"pipes": rows, "median_ratio": medr,
                         "to_scale": medr is not None and abs(medr - 1) <= 0.05,
                         "unknown_ids": sorted(set(opts.known_lengths_ft) - set(lengths))}
        if medr is not None and abs(medr - 1) > 0.05:
            b.warnings.append(f"The known lengths say the drawing is at {medr:.3g} times true length: it is not to "
                              f"scale (or its unit is wrong).")
    if b.blocking_questions():
        return b
    b.diagnostics = diagnostics(tlines, tol, net.report)
    o, f = b.frame["origin"], unit_m
    to_ft = lambda xy: ((xy[0] / f + o[0]) * f / T.FT, (xy[1] / f + o[1]) * f / T.FT)
    b.inp_text = T.write_inp(net, b.diam_in, b.c_of, TITLE, 0.0, to_ft)
    return b


DIAGNOSTIC_KEYS = ("pipes", "nodes", "components", "dead_ends", "ends_joined", "tees_split", "nodes_split",
                   "dropped_pieces")


def diagnostics(tlines: list, tol_m: float, known: dict | None = None, feet=DIAGNOSTIC_FT) -> list[dict]:
    """The tolerance table: the network's counts at each distance in feet (the build's own report reused at its
    tolerance)."""
    out = []
    for ft in feet:
        tm = ft * T.FT
        try:
            r = known if known is not None and abs(tm - tol_m) < 1e-12 else T.build_network(tlines, tm, "xy").report
        except T.TopologyRefused:
            r = {k: None for k in DIAGNOSTIC_KEYS}
        out.append({"tolerance_ft": ft, **{k: r[k] for k in DIAGNOSTIC_KEYS}})
    return out


# ----------------------------------------------------------------------------- the readiness report
NOT_READY = ("Not ready to run: topology and sizes only, needs the operator form. A drawing gives pipes, their sizes "
             "and where they meet, but no source, tank, pump, valve, control or demand. To run it, download the model "
             "and add those in EPANET (or the modelling program you use); the operator form for them here is not "
             "built yet.")
CANNOT_GIVE = ("its scale (if it is a schematic), demands, elevations, roughness (C values), valve status, PRV settings, "
               "which crossings connect (crossing lines are never joined), pump direction, and anything that changes "
               "with time")


def default_words(b: CadBuild) -> str:
    """'the commonest size', or 'the size you gave' when the operator gave the size of pipes with none."""
    return "the size you gave" if b.default_size.get("from") == "given" else "the commonest size"


def readiness(b: CadBuild, prepared: Readiness) -> Readiness:
    """The skeleton's readiness report: the drawing's rows (what came from it, what was derived, assumed and missing,
    labelled 'topology and sizes only'), then importers.prepare_inp's rows for the skeleton."""
    rep = Readiness(name=prepared.name, disinfectant=prepared.disinfectant, sha1=prepared.sha1, sniff=prepared.sniff,
                    not_a_model=prepared.not_a_model, run_cost=prepared.run_cost, prepared=prepared.prepared,
                    same_file=prepared.same_file, flags={**prepared.flags, "cad_skeleton": True})
    d, net = b.drawing, b.net
    rep.add("Built from", f"the CAD drawing {d.name} (DXF, {d.counts.get('version_name', d.version)}), "
                          f"{b.lines_used:,} of its {len(d.lines):,} lines"
            + (f", {len(b.excluded)} layer style{'s' if len(b.excluded) != 1 else ''} left out" if b.excluded else ""),
            FROM_FILE)
    rep.add("What it is", "topology and sizes only, needs the operator form", DERIVED)
    rep.add("Drawing unit", f"{UNIT_WORDS[b.unit]}: {b.unit_how}", ASSUMED if b.unit_how.startswith("chosen") else FROM_FILE)
    rep.add("Pipes and nodes", f"{len(net.pipes):,} pipes and {len(net.nodes):,} nodes, ends joined within "
                               f"{b.tolerance_m / T.FT:.1f} ft ({b.tolerance_m:.2f} m); {net.report['components']} "
                               f"connected part{'s' if net.report['components'] != 1 else ''}, "
                               f"{net.report['dead_ends']:,} dead ends, {net.report['tees_split']:,} "
                               f"cut{'s' if net.report['tees_split'] != 1 else ''} at a tee", DERIVED)
    sc = b.size_counts
    # a label is matched to its nearest line by distance, so a size from a label is derived, not read (the review)
    rep.add("Diameters", f"{sc.get('label', 0):,} pipes from labels (each label to its nearest line), "
                         f"{sc.get('layer', 0):,} from layer names, {sc.get('propagated', 0):,} along runs of pipes, "
                         f"{sc.get('default', 0):,} at {default_words(b)}",
            ASSUMED if sc.get("default", 0) or sc.get("propagated", 0) else DERIVED if sc.get("label", 0) else FROM_FILE)
    cc = b.c_counts
    rep.add("Hazen-Williams C", f"{cc.get('from_material', 0):,} from a material in a label or layer name (new-pipe "
                                f"values of the {T.C_SOURCE}; {T.C_CLASSES_NOTE}), {cc.get('default', 0):,} at "
                                f"{T.C_DEFAULT:g}, the default", ASSUMED)
    rep.add("Lengths", f"from the drawing, in {UNIT_WORDS[b.unit]}"
            + ("; the known lengths agree" if b.scale_check.get("to_scale") else
               f"; the known lengths say {b.scale_check['median_ratio']:.3g} times true length"
               if b.scale_check.get("median_ratio") else ": check them against a few known lengths if the drawing "
                                                         "may not be to scale"), DERIVED)
    rep.add("Location", "a geographic location is in the drawing" if d.geodata else
            "none: the drawing's own coordinates, so mapping taps by GPS needs control taps", FROM_FILE)
    rep.add("Elevations", "every junction at 0 ft: a skeleton has no pressure-control valve, so one elevation changes "
                          "no flow; give real elevations before adding a PRV", ASSUMED)
    rep.add("What a drawing cannot give", CANNOT_GIVE, MISSING)
    rep.add("Source, tanks, pumps, valves, controls, demand", "none: a drawing of mains carries none of them", MISSING)
    for item, value, source in prepared.rows:
        if item == "Text encoding":
            continue
        if item == "File":
            value = "the skeleton model written from your drawing (EPANET .inp)"
        rep.rows.append((item, value, DERIVED if source == FROM_FILE else source))
    rep.blocking = [NOT_READY] + list(prepared.blocking)
    rep.fixes = list(prepared.fixes)
    rep.warnings = list(b.warnings) + list(prepared.warnings)
    return rep


# ----------------------------------------------------------------------------- one call
def build_cad(data: bytes, name: str, opts: CadOptions | None = None, limits: dict | None = None) -> CadBuild:
    """Any file to a CadBuild: read (read_drawing) and built (build_from_drawing)."""
    d = read_drawing(data, name, limits)
    return build_from_drawing(d, opts)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("drawing")
    ap.add_argument("--unit", choices=sorted(UNIT_M))
    ap.add_argument("--tolerance-ft", type=float)
    ap.add_argument("--out")
    ap.add_argument("--form", help="an operator form (CSV or Excel): complete the skeleton with it, each pipe's "
                                   "material and the coordinates kept, and write the completed model to --out "
                                   "(task 19)")
    ap.add_argument("--layer-numbers", choices=("none", "in", "mm", "m"),
                    help="numbers in layer names with no unit beside them: not sizes (none), or their unit")
    ap.add_argument("--default-in", type=float, help="the size of pipes no layer name or label sizes, in inches")
    a = ap.parse_args(argv)
    with open(a.drawing, "rb") as fh:
        data = fh.read()
    b = build_cad(data, os.path.basename(a.drawing), CadOptions(unit=a.unit, tolerance_m=a.tolerance_ft * T.FT
                                                                if a.tolerance_ft else None, to_scale=True,
                                                                layer_numbers=a.layer_numbers,
                                                                default_diameter_in=a.default_in))
    if b.refusal:
        print(b.refusal)
        return 1
    for q in b.questions:
        print("Question:", q["text"])
    if b.net is not None:
        print(json.dumps({k: v for k, v in b.net.report.items() if not k.endswith("_detail")}, indent=1))
        print("sizes:", b.size_counts)
    for w in b.warnings:
        print("Warning:", w)
    if b.inp_text and a.form:          # task 19's review: the form applied to the build, not to its .inp
        from . import operator_form as OF
        return OF.complete_from_files(OF.skeleton_from_build(b, "cad"), a.form, a.out)
    if b.inp_text and a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write(b.inp_text)
        print("Skeleton written to", a.out)
    return 0 if b.inp_text else 1


if __name__ == "__main__":
    sys.exit(main())
