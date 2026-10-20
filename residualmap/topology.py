"""
topology.py: lines to a pipe network, shared by the GIS importer (iteration 5, task 17) and the CAD importer (task 18).

A GIS mains layer or a drawing gives lines, not a network: each line's two ends must be joined to the ends of the lines
it meets, a line drawn through a tee must be split where the branch ends on it, and every node and pipe needs an ID
EPANET accepts.  This module does that, and only that; residualmap/gis_read.py reads the layers.

- Snapping (critique 13).  Line ends closer than the tolerance are joined: scipy cKDTree.query_pairs over every end,
  then union-find, in float64 and in metres (a layer's own unit converted first; never float32, which loses about half
  a foot at State Plane coordinates of 5e6 ft).  Union-find chains: ends A and C 10 ft apart join through an end B 5 ft
  from each.  So every node whose ends spread wider than the tolerance is flagged, and one wider than SPLIT_SPAN times
  the tolerance (a valve cluster or a pump manifold drawn as short lines) is split again by complete linkage at the
  tolerance, so that every end of each part is within the tolerance of every other, and flagged as split.
- Tees.  A line end that lies within the tolerance of another line's interior (not near either of that line's ends)
  splits that line there: at its nearest vertex when that vertex is within the tolerance of every end that asked for
  the cut, else at the ends' foot on the line (since the review of task 17: a vertex near the foot but farther than
  the tolerance from the branch end cut the main for nothing and left the branch apart).  Ends whose feet lie within
  the tolerance of each other along a line share one cut, at the middle of their feet, only when that point is within
  the tolerance of each of them; otherwise each gets its own cut and the snapping joins the cuts.  A pipe crossing
  another without ending on it is not split (crossing mains are not connected).  Since the final verification of
  task 17, a cut more than the tolerance from both ends of the line it cuts (a true tee, not a junction drawn a little
  off) keeps its two piece ends and the ends that asked for it in one node when a wide node is split (_split_units):
  branch ends 5 ft apart on a main drawn with a vertex every 2 ft were put in a node of their own, apart from it.
- A closed line (a loop main drawn as one line that ends where it starts) longer than twice the tolerance is cut at
  the middle of its length into two pipes through a new node, since EPANET cannot run a pipe from a node to itself
  (the review of task 17: it was dropped), and since its final verification so is any piece longer than that whose
  two ends the snapping puts in one node (a loop whose ends join through a third line's end).  Other pieces whose two
  ends fall into one node (a piece not much longer than the tolerance) are dropped and listed.  The report gives the components, the dead ends, the ends joined, the
  tee cuts, the flagged and split nodes, the loops cut, the dropped pieces and the lines drawn twice.
- Crowding (the review of task 17).  Before anything is joined, the line ends are counted in squares the size of the
  tolerance; a layer with more than MAX_ENDS_PER_SQUARE ends in one square, or whose squares hold far more ends than
  any real network (CROWDING_MAX), is refused: its unit or coordinate system is wrong (longitude and latitude read as
  feet put a whole town inside a foot), and joining it would list every pair of its ends.
- IDs (critique 8).  A node's ID comes from its rounded position (N<x>_<y> in the layer's own unit, or N<lat>_<lon> at
  5 decimals for longitude and latitude), so it does not depend on the order of the features: a re-imported layer whose
  node moved gets a new ID, and a tap map naming the old one is refused, never silently repointed.  A pipe's ID is the
  layer's own ID field when it has one (sanitized), else P plus 8 hex digits of a hash of its two node IDs and its
  length.  Any ID with a space, a ';', a quote or more than 31 bytes goes through the written mapping table.
- Diameters (critique 14).  parse_diameter reads a size in a text: 8", 8 in, 200 mm, DN150, 1965 8" CI (8, not 1965),
  L=350' 8" (8, not 350), and flags 12x8 RED (a reducer) and 8 & 6 (two sizes).  Since the review of task 17 a size
  labelled twice in nominal equivalents (6" (150mm), DN 150 (6")), within DUAL_LABEL_REL of each other, is one size,
  its inch value, and a coupling of one size (6"x6") is that size, not a reducer.  Since its final verification a
  reducing tee of three sizes (8"x8"x6") is a reducer again, a number that overflows a float is no size, a text over
  SIZE_TEXT_MAX characters is unreadable, and for a KML's name or description 'in' is inches only when no word
  follows it.  A bare number has no unit:
  bare_unit reads the layer's bare numbers together by their range (2 to 48 is inches, 50 to 1,200 is mm, 0.05 to 1.2
  is metres); anything else is asked.
- Hazen-Williams C from the material: the new-pipe values of the EPANET 2 Users Manual (Rossman 2000, Table 3.2,
  'Roughness Coefficients for New Pipe'), the middle of each range; 130 is the labelled default for a material it does
  not list.  The table names cast iron, concrete or concrete lined, galvanized iron, plastic, steel and vitrified clay
  only: the words that put a layer's material in one of those classes are this tool's own (ductile iron, asbestos
  cement, PCCP and cement-lined cast iron are read as concrete lined, HDPE, PE and PEX as plastic).  The install year
  is read and shown but does not change C: no aging table was checked for this build, so no pipe's age reaches the
  model (the old-pipe factor scales wall decay by C alone, and new-pipe C ranks cast iron above ductile iron).
- write_inp: the skeleton as an EPANET .inp in GPM and feet (lengths in ft, diameters in inches, 12 significant
  digits), every junction at one elevation and no demand, with no source, tank, pump, valve or control: the readiness
  report says it is not ready until those are given.

Never wntr.gis.connect_lines: it took 17 s on ky4 and was still running after 22 min on Net6 (the plan's probe).
"""
from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass, field

import numpy as np

FT = 0.3048                        # the international foot, in metres
FT_US = 1200.0 / 3937.0            # the US survey foot
UNIT_M = {"ft": FT, "us_ft": FT_US, "m": 1.0}
UNIT_WORDS = {"ft": "feet", "us_ft": "US survey feet", "m": "metres"}
DEFAULT_TOL_M = {"gis": 6 * FT, "kml": 10.0}   # 6 ft for a GIS layer; 10 m for lines drawn by hand in Google Earth
SPAN_FLAG = 1.0                    # a node whose ends spread wider than this many tolerances is flagged
SPLIT_SPAN = 2.0                   # ... and one wider than this many is split again (complete linkage at the tolerance)
BIG_CLUSTER = 2000                 # a node of more ends than this is split on a grid, not by linkage (memory)
MAX_ID_BYTES = 31                  # EPANET's longest ID
SEG_BUCKETS = 32                   # segment length classes for the tee search (powers of 2)
DIAGNOSTIC_FACTORS = (0.5, 1.0, 2.0)   # the diagnostic table: half the tolerance, the tolerance, twice it
C_DEFAULT = 130.0                  # Hazen-Williams C for a material the table does not list (roughness_factor = 1)
DUPLICATE_LENGTH_TOL = 0.01        # two pipes between the same nodes within 1% in length: a line drawn twice?
LOOP_CUT_SPAN = 2.0                # a closed piece longer than this many tolerances is cut at its middle (the review)
MAX_ENDS_PER_SQUARE = 200          # line ends in one tolerance-sized square: more is a unit or CRS mix-up (the review)
CROWDING_MAX = 20.0                # ... and the mean, over the ends, of the ends in each one's square (about 3 in a
                                   # real network: the ends of one junction)
DUAL_LABEL_REL = 0.05              # two sizes within 5% of each other are one size labelled twice (6" and 150 mm)

# EPANET 2 Users Manual (Rossman 2000), Table 3.2, Hazen-Williams C for new pipe: the middle of each range.  The
# words that name each material in a layer's material field (lower case, whole words).
C_TABLE = (
    ("cast iron", 135.0, ("ci", "cip", "cast", "castiron", "cast iron")),
    ("concrete or concrete lined", 130.0, ("di", "dip", "ductile", "ductile iron", "dicl", "cl di", "concrete", "rcp",
                                            "pccp", "cement", "cement lined", "ac", "acp", "asbestos", "asbestos cement",
                                            "transite", "cicl", "cl")),
    ("galvanized iron", 120.0, ("gi", "galv", "galvanized", "galvanised", "gip")),
    ("plastic", 145.0, ("pvc", "hdpe", "pe", "cpvc", "pex", "plastic", "poly", "polyethylene", "mdpe", "pvco")),
    ("steel", 145.0, ("st", "stl", "steel", "ws", "welded steel")),
    ("vitrified clay", 110.0, ("vc", "vcp", "clay", "vitrified clay")),
)
C_SOURCE = "EPANET 2 Users Manual (Rossman 2000), Table 3.2, new pipe, middle of each range"
C_CLASSES_NOTE = ("the table names cast iron, concrete or concrete lined, galvanized iron, plastic, steel and vitrified "
                  "clay only; reading ductile iron, asbestos cement and PCCP as concrete lined (130) and HDPE, PE and "
                  "PEX as plastic (145) is this tool's own choice")


class TopologyRefused(ValueError):
    """Lines that cannot be built into a network (none, or only points), in plain words."""


# ----------------------------------------------------------------------------- diameters
_QUOTES = (("″", '"'), ("“", '"'), ("”", '"'), ("''", '"'), ("’", "'"), ("′", "'"), ("‘", "'"), ("ø", ""), ("Ø", ""),
           ("⌀", ""))
_NUM = r"(\d+(?:\.\d+)?\s*-\s*\d+/\d+|\d+\s+\d+/\d+|\d+/\d+|\d+(?:[.,]\d+)?)"
_TOKEN = re.compile(r"(?<![\w./])(nps\s*|dn\s*|l\s*=\s*|len(?:gth)?\s*=?\s*)?" + _NUM +
                    r"\s*(mm|cm|m(?![a-z])|in(?:ch(?:es)?)?(?![a-z])|\"|'|ft(?![a-z])|feet|foot)?", re.I)
# one run of spaces around an optional quote, so a failed match cannot try every split of a run (the final
# verification of task 17)
_BY = r"\s*(?:\"\s*)?(?:x|×|\*|by)\s*"
_REDUCER = re.compile(r"(?<![\w.])" + _NUM + _BY + _NUM + r"\s*\"?(?![\w.])", re.I)
_REDUCER_NEXT = re.compile(_NUM + _BY + _NUM + r"\s*\"?(?![\w.])", re.I)   # a third size after a coupling's two
SIZE_TEXT_MAX = 500                # a size text longer than this is no size (the final verification of task 17)
# three sizes or more joined by x (a reducing tee, '12x12x8 TEE'): runs of digits, quotes and spaces between the x's,
# so a run cannot be split two ways
_CHAIN = re.compile(r'(?<![\w.])[\d./"\s]+(?:[x×*][\d./"\s]+){2,}', re.I)
_TO_IN = {"in": 1.0, "mm": 1 / 25.4, "cm": 10 / 25.4, "m": 1000 / 25.4}
BARE_RANGES = (("in", 2.0, 48.0), ("mm", 50.0, 1200.0), ("m", 0.05, 1.2))


def _number(s: str) -> float:
    s = s.strip()
    s = s.replace(",", "") if re.fullmatch(r"\d{1,3}(?:,\d{3})+", s) else s.replace(",", ".")   # 1,200 mm; 6,0
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(?:-\s*|\s+)(\d+)/(\d+)", s)
    if m:
        x = float(m.group(1)) + float(m.group(2)) / max(float(m.group(3)), 1e-12)
        if not math.isfinite(x):
            raise ValueError("not a finite number")
        return x
    m = re.fullmatch(r"(\d+)/(\d+)", s)
    if m:
        x = float(m.group(1)) / max(float(m.group(2)), 1e-12)
    else:
        x = float(s)
    if not math.isfinite(x):           # 400 nines overflow to inf (the final verification of task 17)
        raise ValueError("not a finite number")
    return x


def parse_diameter(v, strict_in: bool = False) -> dict:
    """A pipe size as a layer or a drawing writes it: {'value', 'unit' ('in', 'mm', 'cm', 'm', or None for a bare
    number), 'flag' (None, 'empty', 'reducer', 'two sizes', 'unreadable'), 'text'}.  A number with a unit wins over a
    bare one; a number marked in feet (') or after L= is a length and is passed over, and so is a year (1850 to 2030)
    beside another number; 12x8 is a reducer and 8 & 6 two sizes, both flagged with no value.  A coupling of one size
    (6"x6") is that size, and two sizes with units within DUAL_LABEL_REL of each other (6" (150mm)) are one size, its
    inch value (the review of task 17).  Since its final verification: a number that overflows (400 nines) is no
    size; a reducing tee of three sizes (8"x8"x6", 12x12x8 TEE) is a reducer, and a tee of one size (8"x8"x8") that
    size; a text longer than SIZE_TEXT_MAX is unreadable; and strict_in (a KML's name or description)
    takes 'in' as inches only when no word follows it ('Line 3 in Oak St' gives a bare 3)."""
    out = {"value": None, "unit": None, "flag": None, "text": "" if v is None else str(v)}
    if v is None or (isinstance(v, float) and math.isnan(v)):
        out["flag"] = "empty"
        return out
    if isinstance(v, bool):
        out["flag"] = "unreadable"
        return out
    if isinstance(v, (int, float, np.integer, np.floating)):
        x = float(v)
        if x > 0 and math.isfinite(x):
            out["value"] = x
        else:
            out["flag"] = "empty" if x == 0 else "unreadable"
        return out
    s = str(v)
    if len(s) > SIZE_TEXT_MAX:
        out["flag"], out["text"] = "unreadable", s[:SIZE_TEXT_MAX]
        return out
    for a, b in _QUOTES:
        s = s.replace(a, b)
    s = s.strip()
    if not s or s.lower() in ("na", "n/a", "none", "null", "unknown", "unk", "?", "-", "0"):
        out["flag"] = "empty" if s in ("", "0") or s.lower() in ("none", "null") else "unreadable"
        return out
    for run in _CHAIN.finditer(s):        # a tee of three sizes, quoted or not (the final verification of task 17)
        vals = []
        for piece in re.split(r"[x×*]", run.group(0), flags=re.I):
            try:
                vals.append(_number(piece.replace('"', "").strip()))
            except ValueError:
                vals = []
                break
        if len(vals) >= 3 and len(set(vals)) > 1:
            out["flag"] = "reducer"
            return out
    red_word = re.search(r"(?<![a-z])(red|reducer|reducing)(?![a-z])", s, re.I)
    m = _REDUCER.search(s)
    if m and not red_word:
        try:
            same = _number(m.group(1)) == _number(m.group(2))
        except ValueError:
            same = False
        third = _REDUCER_NEXT.match(s, m.start(2))    # a tee of three sizes (8"x8"x6"): the final verification
        if same and third:
            try:
                same = _number(third.group(1)) == _number(third.group(2))
            except ValueError:
                same = False
        if same:                 # a coupling of one size (6"x6"): that size (the review of task 17), not a reducer
            end = third.end() if third else m.end()
            s = s[:m.start()] + m.group(1) + ('"' if '"' in s[m.start():end] else "") + s[end:]
            m = None
    if m or red_word and len(re.findall(r"\d+", s)) >= 2:
        out["flag"] = "reducer"
        return out
    with_unit, bare = [], []
    for m in _TOKEN.finditer(s):
        pre, num, unit = (m.group(1) or "").lower().replace(" ", ""), m.group(2), (m.group(3) or "").lower()
        try:
            x = _number(num)
        except ValueError:
            continue
        if pre.startswith(("l=", "len")) or unit in ("'", "ft", "feet", "foot"):
            continue                                  # a length
        if pre.startswith("dn"):
            with_unit.append((x, "mm"))
        elif strict_in and unit == "in" and s[m.end():m.end() + 50].lstrip()[:1].isalpha():
            bare.append(x)                            # 'in Oak St': a word, not inches (strict_in)
        elif pre.startswith("nps") or unit in ('"',) or unit.startswith("in"):
            with_unit.append((x, "in"))
        elif unit in ("mm", "cm"):
            with_unit.append((x, unit))
        elif unit == "m":
            if x <= 2.0:
                with_unit.append((x, "m"))
        else:
            bare.append(x)
    if with_unit:
        sizes = {round(x * _TO_IN[u], 3) for x, u in with_unit if x > 0}
        if len(sizes) > 1 and max(sizes) > (1 + DUAL_LABEL_REL) * min(sizes):
            out["flag"] = "two sizes"
            return out
        if len(sizes) > 1:       # one size labelled twice, in nominal equivalents (6" (150mm)): its inch value
            with_unit = [w for w in with_unit if w[1] == "in"] or with_unit
        x, u = with_unit[0]
        if x <= 0:
            out["flag"] = "unreadable"
            return out
        out.update(value=x, unit=u)
        return out
    years = [x for x in bare if 1850 <= x <= 2030 and float(x).is_integer()]
    if years and len(bare) > len(years):
        bare = [x for x in bare if x not in years]
    bare = [x for x in bare if x > 0]
    if not bare:
        out["flag"] = "unreadable"
        return out
    if len(set(bare)) > 1:
        out["flag"] = "two sizes"
        return out
    out["value"] = bare[0]
    return out


def bare_unit(values) -> str | None:
    """The unit of a layer's bare diameters, read together by their range: every one in 2 to 48 is inches, 50 to
    1,200 mm, 0.05 to 1.2 metres; None when they fit none of these (asked)."""
    xs = [float(x) for x in values if x is not None and x > 0]
    if not xs:
        return None
    for unit, lo, hi in BARE_RANGES:
        if all(lo <= x <= hi for x in xs):
            return unit
    return None


def to_inches(value: float, unit: str) -> float:
    return float(value) * _TO_IN[unit]


def hazen_c(material) -> tuple[float, str | None]:
    """(C, the table's material class or None) for a material field's value: the new-pipe value of C_TABLE, or
    C_DEFAULT (None) for a value it does not list."""
    if material is None:
        return C_DEFAULT, None
    s = re.sub(r"[^a-z ]+", " ", str(material).lower()).split()
    joined = " ".join(s)
    for cls, c, words in C_TABLE:
        for w in words:
            if (" " in w and w in joined) or w in s:
                return c, cls
    return C_DEFAULT, None


# ----------------------------------------------------------------------------- IDs
_BAD_ID = re.compile(r"[^\w.\-:+/#@]")          # \w takes accented letters, which EPANET reads as UTF-8


def clean_id(raw: str) -> tuple[str, str | None]:
    """(an EPANET-safe ID, why it changed or None): characters other than letters, digits and _ . - : + / # @ become
    '_' (accented letters are kept), and an ID over MAX_ID_BYTES bytes (in UTF-8) is cut and given 8 hex digits of its hash, so two long IDs stay apart."""
    s = str(raw).strip()
    why = []
    t = _BAD_ID.sub("_", s)
    if t != s:
        why.append("a space, ';', a quote or another character EPANET IDs do not take" if s else "empty")
    t = re.sub(r"_+", "_", t).strip("_") or "X"
    if len(t.encode("utf-8")) > MAX_ID_BYTES:
        h = hashlib.sha1(s.encode("utf-8")).hexdigest()[:8]
        head = t
        while len((head + "_" + h).encode("utf-8")) > MAX_ID_BYTES:
            head = head[:-1]
        t = head + "_" + h
        why.append(f"longer than {MAX_ID_BYTES} characters")
    return t, ("; ".join(why) if why else None)


def _coord_text(v: float, decimals: int) -> str:
    t = f"{v:.{decimals}f}"
    return "0" if decimals == 0 and t in ("-0",) else t


def node_id(raw_xy, style: str) -> str:
    """A node's ID from its rounded position: N<x>_<y> in the layer's own unit (style 'xy'), or N<lat>_<lon> at 5
    decimals (about a metre; style 'lonlat', raw_xy given as lon, lat)."""
    x, y = float(raw_xy[0]), float(raw_xy[1])
    if style == "lonlat":
        return f"N{_coord_text(y, 5)}_{_coord_text(x, 5)}"
    return f"N{_coord_text(round(x), 0)}_{_coord_text(round(y), 0)}"


# ----------------------------------------------------------------------------- the network
@dataclass
class Line:
    """One part of one feature: its vertices in metres (a local frame, float64) and in the layer's own coordinates."""
    feature: int
    part: int
    xy: np.ndarray            # (n, 2) metres
    raw: np.ndarray           # (n, 2) the layer's coordinates (x, y or lon, lat)


@dataclass
class Network:
    tol_m: float
    nodes: dict = field(default_factory=dict)        # node ID -> {'xy': (x, y) m, 'raw': (x, y), 'ends': n, 'span_m'}
    pipes: list = field(default_factory=list)        # {'id', 'start', 'end', 'length_m', 'feature', 'part', 'piece',
                                                     #  'pieces', 'xy': [(x, y) m, ...]}
    report: dict = field(default_factory=dict)
    end_node: dict = field(default_factory=dict)     # (feature, part, 0 or 1) -> node ID: each line end's node
    flags: list = field(default_factory=list)        # plain sentences

    def degree(self) -> dict:
        d = {n: 0 for n in self.nodes}
        for p in self.pipes:
            d[p["start"]] += 1
            d[p["end"]] += 1
        return d


class _UF:
    def __init__(self, n: int):
        self.p = np.arange(n)

    def find(self, i: int) -> int:
        p = self.p
        r = i
        while p[r] != r:
            r = p[r]
        while p[i] != r:
            p[i], i = r, p[i]
        return r

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if ra < rb:
                self.p[rb] = ra
            else:
                self.p[ra] = rb


def _clean_line(xy: np.ndarray, raw: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """Consecutive repeated vertices removed: (xy, raw, how many were removed)."""
    if len(xy) < 2:
        return xy, raw, 0
    keep = np.ones(len(xy), bool)
    keep[1:] = np.any(np.diff(xy, axis=0) != 0, axis=1)
    return xy[keep], raw[keep], int((~keep).sum())


def _seg_dist(P, A, B):
    """Distance from each point P (k, 2) to its segment A-B (k, 2), and the position t in [0, 1] of its foot."""
    AB = B - A
    L2 = np.maximum((AB * AB).sum(1), 1e-300)
    t = np.clip(((P - A) * AB).sum(1) / L2, 0.0, 1.0)
    F = A + t[:, None] * AB
    return np.hypot(F[:, 0] - P[:, 0], F[:, 1] - P[:, 1]), t


def _tee_candidates(lines: list[Line], ends: np.ndarray, end_line: np.ndarray, tol: float) -> list[tuple]:
    """(line q, arc position on q, end index, distance) for every line end within tol of another line's interior."""
    from scipy.spatial import cKDTree
    segA, segB, segL, segS, segArc = [], [], [], [], []
    for q, ln in enumerate(lines):
        if len(ln.xy) < 2:
            continue
        a, b = ln.xy[:-1], ln.xy[1:]
        seglen = np.hypot(*(b - a).T)
        arc = np.concatenate([[0.0], np.cumsum(seglen)[:-1]])
        segA.append(a)
        segB.append(b)
        segL.append(np.full(len(a), q))
        segS.append(seglen)
        segArc.append(arc)
    if not segA:
        return []
    A, B = np.concatenate(segA), np.concatenate(segB)
    SL, S, ARC = np.concatenate(segL), np.concatenate(segS), np.concatenate(segArc)
    mid = (A + B) / 2
    half = S / 2
    cls = np.ceil(np.log2(np.maximum(half, tol))).astype(int)
    tree_e = cKDTree(ends)
    found = {}
    for c in np.unique(cls):
        idx = np.where(cls == c)[0]
        r = float(2.0 ** c) + tol
        tree_s = cKDTree(mid[idx])
        hits = tree_s.query_ball_tree(tree_e, r)
        si, ei = [], []
        for k, hs in enumerate(hits):
            if hs:
                si.extend([idx[k]] * len(hs))
                ei.extend(hs)
        if not si:
            continue
        si, ei = np.array(si), np.array(ei)
        own = SL[si] == end_line[ei]
        si, ei = si[~own], ei[~own]
        if not len(si):
            continue
        d, t = _seg_dist(ends[ei], A[si], B[si])
        ok = d <= tol
        for s_, e_, d_, t_ in zip(si[ok], ei[ok], d[ok], t[ok]):
            q = int(SL[s_])
            key = (q, int(e_))
            pos = float(ARC[s_] + t_ * S[s_])
            if key not in found or d_ < found[key][1]:
                found[key] = (pos, float(d_))
    out = []
    for (q, e), (pos, d) in found.items():
        ln = lines[q]
        if (np.hypot(*(ends[e] - ln.xy[0])) <= tol) or (np.hypot(*(ends[e] - ln.xy[-1])) <= tol):
            continue                                   # an end-to-end join, not a tee
        out.append((q, pos, e, d))
    return out


def _arc(xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    seglen = np.hypot(*(xy[1:] - xy[:-1]).T)
    return seglen, np.concatenate([[0.0], np.cumsum(seglen)])


def _at(xy: np.ndarray, seglen: np.ndarray, arc: np.ndarray, pos: float) -> tuple[int, float, np.ndarray]:
    """(segment k, fraction f along it, the point) at arc position pos of a polyline."""
    k = int(np.searchsorted(arc, pos, side="right") - 1)
    k = min(max(k, 0), len(seglen) - 1)
    f = (pos - arc[k]) / max(seglen[k], 1e-300)
    return k, f, xy[k] + f * (xy[k + 1] - xy[k])


def _tee_cuts(ln: Line, cand: list[tuple[float, int]], ends: np.ndarray, tol: float) -> list[tuple[float, list]]:
    """The cuts of one line, as (arc position, the indices of the ends that asked for it): candidates (arc position of
    an end's foot, the end) whose feet lie within tol of each other along the line share one cut at the middle of
    their feet when that point is within tol of each of their ends; otherwise each end gets its own cut (the review of
    task 17: the middle candidate's cut could lie farther than tol from another end, which then joined nothing)."""
    groups: list[list] = []
    for pos, e in sorted(cand):
        if groups and pos - groups[-1][-1][0] <= tol:
            groups[-1].append((pos, e))
        else:
            groups.append([(pos, e)])
    seglen, arc = _arc(ln.xy)
    out = []
    for g in groups:
        if len(g) == 1:
            out.append((g[0][0], [g[0][1]]))
            continue
        E = ends[[e for _, e in g]]
        mid = (g[0][0] + g[-1][0]) / 2
        P = _at(ln.xy, seglen, arc, mid)[2]
        if float(np.hypot(*(E - P).T).max()) <= tol:
            out.append((mid, [e for _, e in g]))
        else:
            out += [(pos, [e]) for pos, e in g]
    return out


def _split_line(ln: Line, cuts_in: list[tuple[float, list]], tol: float, ends: np.ndarray
                ) -> tuple[list[tuple[np.ndarray, np.ndarray]], int, list[list]]:
    """ln cut at the arc positions of cuts_in ((position, the indices of the ends that asked for the cut)), each moved
    to the line's nearest vertex when that vertex is within tol of the foot and of every end that asked for it (else
    the foot is inserted): its pieces as (xy, raw), the number of cuts made, and for each cut, in order along the line
    (the boundary between piece i and piece i + 1), the ends that asked for it (the final verification of task 17:
    the snapping keeps them in the cut's node)."""
    xy, raw = ln.xy, ln.raw
    seglen, arc = _arc(xy)
    total = arc[-1]
    cuts = []
    for pos, E_idx in sorted(cuts_in, key=lambda c: c[0]):
        if pos <= 1e-9 or pos >= total - 1e-9:
            continue
        E = ends[list(E_idx)]
        k, f, foot = _at(xy, seglen, arc, pos)
        vd = np.hypot(*(xy - foot).T)
        vd[0] = vd[-1] = np.inf
        j = int(np.argmin(vd))
        if vd[j] <= tol and float(np.hypot(*(np.asarray(E, float) - xy[j]).T).max()) <= tol:
            P = ("v", j)
        else:
            P = ("p", k, f)
        cuts.append((P, list(E_idx)))
    if not cuts:
        return [(xy, raw)], 0, []
    # build the list of points with the new points inserted, then cut at their indices
    pts, rpts = [xy[0]], [raw[0]]
    inserts: dict = {}
    for c, _ in cuts:
        if c[0] == "p":
            inserts.setdefault(c[1], []).append(c[2])
    at_vertex, at_insert = {0: 0}, {}
    for k in range(len(xy) - 1):
        for f in sorted(set(inserts.get(k, []))):
            pts.append(xy[k] + f * (xy[k + 1] - xy[k]))
            rpts.append(raw[k] + f * (raw[k + 1] - raw[k]))
            at_insert[(k, f)] = len(pts) - 1
        pts.append(xy[k + 1])
        rpts.append(raw[k + 1])
        at_vertex[k + 1] = len(pts) - 1
    pts, rpts = np.array(pts), np.array(rpts)
    asked: dict = {}
    for c, E_idx in cuts:
        i = at_vertex[c[1]] if c[0] == "v" else at_insert[(c[1], c[2])]
        if 0 < i < len(pts) - 1:
            asked.setdefault(i, []).extend(E_idx)
    cut_at = sorted(asked)
    pieces, start = [], 0
    for i in cut_at + [len(pts) - 1]:
        if i > start:
            pieces.append((pts[start:i + 1], rpts[start:i + 1]))
        start = i
    return pieces, len(cut_at), [sorted(set(asked[i])) for i in cut_at]


def _complete_linkage(P: np.ndarray, tol: float) -> np.ndarray:
    """Labels of P's points cut into groups whose every two points are within tol (complete linkage), or on a grid of
    tol cells for a group too large to link."""
    if len(P) > BIG_CLUSTER:
        g = np.floor(P / tol).astype(np.int64)
        _, lab = np.unique(g, axis=0, return_inverse=True)
        return lab.ravel()
    from scipy.cluster.hierarchy import fcluster, linkage
    Z = linkage(P, method="complete")
    return fcluster(Z, t=tol, criterion="distance") - 1


def _split_units(P: np.ndarray, unit: np.ndarray, tol: float) -> np.ndarray:
    """Labels of P's points cut into groups by complete linkage at tol, where the points of one unit (unit[i], the same
    value for points that must stay together) are never parted: the distance of two units is the largest distance
    between their points (the final verification of task 17: complete linkage put two branch ends in a node of their
    own, apart from the cut each had asked for on the main).  With every unit a single point it is _complete_linkage."""
    uniq, inv = np.unique(unit, return_inverse=True)
    inv = inv.ravel()
    if len(uniq) == len(P):
        return _complete_linkage(P, tol)
    if len(uniq) == 1:
        return np.zeros(len(P), dtype=np.int64)
    first = np.full(len(uniq), len(P))
    np.minimum.at(first, inv, np.arange(len(P)))
    if len(P) > BIG_CLUSTER:                 # too large to link: a grid of tol cells, each unit in its first point's
        g = np.floor(P[first] / tol).astype(np.int64)
        _, lab = np.unique(g, axis=0, return_inverse=True)
        return lab.ravel()[inv]
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist, squareform
    order = np.argsort(inv, kind="stable")
    starts = np.searchsorted(inv[order], np.arange(len(uniq)))
    D = squareform(pdist(P[order]))
    DU = np.maximum.reduceat(np.maximum.reduceat(D, starts, axis=0), starts, axis=1)
    np.fill_diagonal(DU, 0.0)
    lab = fcluster(linkage(squareform(DU, checks=False), method="complete"), t=tol, criterion="distance") - 1
    return lab[inv]


def _halves(xy: np.ndarray, raw: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    """A piece cut at the middle of its length into two, through a point at the middle (or the vertex there)."""
    seglen, arc = _arc(xy)
    k, f, P = _at(xy, seglen, arc, arc[-1] / 2)
    if f <= 1e-12:                      # the middle is a vertex (k >= 1, since the length is positive)
        return [(xy[:k + 1], raw[:k + 1]), (xy[k:], raw[k:])]
    R = raw[k] + f * (raw[k + 1] - raw[k])
    return [(np.vstack([xy[:k + 1], P]), np.vstack([raw[:k + 1], R])),
            (np.vstack([P, xy[k + 1:]]), np.vstack([R, raw[k + 1:]]))]


def _span(P: np.ndarray) -> float:
    if len(P) < 2:
        return 0.0
    if len(P) > BIG_CLUSTER:
        return float(np.hypot(*(P.max(0) - P.min(0))))
    from scipy.spatial.distance import pdist
    return float(pdist(P).max())


CROWDED = ("The layer's line ends are crowded together far more than in any water network: {top:,} of its {n:,} line "
           "ends fall in one square {tol_ft:.1f} ft across (the joining distance), and an end's square holds {mean:,.0f} "
           "ends on average. Its coordinates' unit or coordinate system is probably wrong (longitude and latitude "
           "read as feet or metres, for example), so nothing was joined. Check the unit, or the coordinate system the "
           "file names.")


def _refuse_crowded(ends: np.ndarray, tol: float) -> None:
    """TopologyRefused when the line ends are counted in squares of the tolerance and one square holds more than
    MAX_ENDS_PER_SQUARE, or an end's square holds more than CROWDING_MAX on average (the review of task 17: a layer
    in longitude and latitude read as feet put every end within the tolerance of every other, and query_pairs then
    listed every pair, quadratic in time and memory)."""
    if len(ends) < 2:
        return
    g = np.floor((ends - ends.min(axis=0)) / tol).astype(np.int64)
    counts = np.unique(g, axis=0, return_counts=True)[1].astype(float)
    top, mean = int(counts.max()), float((counts ** 2).sum() / len(ends))
    if top > MAX_ENDS_PER_SQUARE or mean > CROWDING_MAX:
        raise TopologyRefused(CROWDED.format(top=top, n=len(ends), tol_ft=tol / FT, mean=mean))


def _cut_loop(xy: np.ndarray, raw: np.ndarray, tol: float) -> list[tuple[np.ndarray, np.ndarray]]:
    """A closed piece (its two ends within tol) longer than LOOP_CUT_SPAN times tol, cut at the middle of its length
    into two (a loop main drawn as one line: EPANET cannot run a pipe from a node to itself); any other piece as it
    is."""
    if len(xy) < 2 or float(np.hypot(*(xy[0] - xy[-1]))) > tol:
        return [(xy, raw)]
    if _arc(xy)[1][-1] <= LOOP_CUT_SPAN * tol:
        return [(xy, raw)]
    return _halves(xy, raw)


def build_network(lines: list[Line], tol_m: float, id_style: str = "xy", pipe_ids: dict | None = None) -> Network:
    """The network of a layer's lines at tolerance tol_m (metres): tees split, ends snapped, wide nodes split, IDs given.
    pipe_ids: {feature index: the layer's own pipe ID (already clean)}; pieces of one feature get _1, _2 ... in order
    along it, and a feature's second part _p2 (none: hashed IDs)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    if tol_m <= 0 or not math.isfinite(tol_m):
        raise TopologyRefused("The snapping tolerance must be a positive distance.")
    removed_dupes = not_finite = 0
    work = []
    for ln in lines:
        xy, raw = np.asarray(ln.xy, float), np.asarray(ln.raw, float)
        if xy.ndim != 2 or len(xy) == 0 or not np.isfinite(xy).all():
            not_finite += 1              # the readers leave these out first; a line of NaN would stop cKDTree
            continue
        xy, raw, k = _clean_line(xy, raw)
        removed_dupes += k
        if len(xy) >= 2:
            work.append(Line(ln.feature, ln.part, xy, raw))
    if not work:
        raise TopologyRefused("The layer has no line with two different points, so there is no pipe to build.")
    ends = np.array([p for ln in work for p in (ln.xy[0], ln.xy[-1])])
    end_line = np.repeat(np.arange(len(work)), 2)
    _refuse_crowded(ends, tol_m)
    # tees: a line end on another line's interior splits that line there
    tees = _tee_candidates(work, ends, end_line, tol_m)
    by_line: dict[int, list] = {}
    for q, pos, e, d in tees:
        by_line.setdefault(q, []).append((pos, e))
    pieces = []          # (line index, piece number, xy, raw)
    n_cuts = loops_cut = 0
    tee_bounds = []      # (piece before the cut, piece after it, the ends that asked for the cut)
    line_pieces = []     # (first piece, last piece) of each line
    for q, ln in enumerate(work):
        cuts_q = _tee_cuts(ln, by_line[q], ends, tol_m) if q in by_line else []
        parts, cuts, asked = _split_line(ln, cuts_q, tol_m, ends) if cuts_q else ([(ln.xy, ln.raw)], 0, [])
        n_cuts += cuts
        k, first_of, last_of = 0, [], []
        for pxy, praw in parts:
            first_of.append(len(pieces))
            for sub in _cut_loop(pxy, praw, tol_m):     # a closed piece becomes two (the review of task 17)
                pieces.append((q, k, *sub))
                k += 1
            last_of.append(len(pieces) - 1)
        loops_cut += k - len(parts)
        line_pieces.append((first_of[0], last_of[-1]))
        tee_bounds += [(last_of[i], first_of[i + 1], e) for i, e in enumerate(asked)]
    # snap every piece end
    P = np.array([p for _, _, pxy, _ in pieces for p in (pxy[0], pxy[-1])])
    R = np.array([p for _, _, _, praw in pieces for p in (praw[0], praw[-1])])
    # a tee cut more than the tolerance from both ends of the line it cuts (a true tee, not a junction drawn a little
    # off) keeps its two piece ends and the ends that asked for it in one node however a wide node is split (the
    # final verification of task 17: complete linkage put branch ends in a node of their own, apart from the main
    # each had asked to cut): each such group is one unit of the split
    must = _UF(len(P))
    for before, after, asked_by in tee_bounds:
        q = pieces[before][0]
        cut = pieces[before][2][-1]
        if min(float(np.hypot(*(cut - work[q].xy[0]))), float(np.hypot(*(cut - work[q].xy[-1])))) <= tol_m:
            continue
        must.union(2 * before + 1, 2 * after)
        for e in asked_by:
            fp, lp = line_pieces[e // 2]
            must.union(2 * before + 1, 2 * fp if e % 2 == 0 else 2 * lp + 1)
    unit = np.array([must.find(i) for i in range(len(P))]) if tee_bounds else np.arange(len(P))
    uf = _UF(len(P))
    for a, b in cKDTree(P).query_pairs(tol_m, output_type="ndarray"):
        uf.union(int(a), int(b))
    roots = np.array([uf.find(i) for i in range(len(P))])
    members: dict[int, list] = {}
    for i, r in enumerate(roots):
        members.setdefault(int(r), []).append(i)
    label = np.empty(len(P), dtype=np.int64)
    clusters = []        # list of index arrays
    flagged, split = [], []
    for r in sorted(members):
        idx = np.array(sorted(members[r]))
        span = _span(P[idx])
        if span > SPLIT_SPAN * tol_m:
            lab = _split_units(P[idx], unit[idx], tol_m)
            subs = [idx[lab == u] for u in np.unique(lab)]
            split.append({"ends": len(idx), "span_m": span, "parts": len(subs)})
            for s in subs:
                label[s] = len(clusters)
                clusters.append(s)
        else:
            if span > SPAN_FLAG * tol_m:
                flagged.append({"ends": len(idx), "span_m": span})
            label[idx] = len(clusters)
            clusters.append(idx)
    # a piece whose two ends fall into one node, longer than LOOP_CUT_SPAN times the tolerance, is cut at its middle
    # through a new node (the final verification of task 17: a loop main whose ends join through a third line's end
    # was dropped; _cut_loop sees only a piece whose own two ends are within the tolerance)
    late_loops = {}
    extra_P, extra_R = [], []
    for k, (q, piece, pxy, praw) in enumerate(pieces):
        if label[2 * k] == label[2 * k + 1] and float(np.hypot(*(pxy[1:] - pxy[:-1]).T).sum()) > LOOP_CUT_SPAN * tol_m:
            halves = _halves(pxy, praw)
            n0 = len(P) + len(extra_P)
            extra_P += [halves[0][0][-1]] * 2
            extra_R += [halves[0][1][-1]] * 2
            clusters.append(np.array([n0, n0 + 1]))
            late_loops[k] = (halves, len(clusters) - 1)
    if extra_P:
        P, R = np.vstack([P, extra_P]), np.vstack([R, extra_R])
        loops_cut += len(late_loops)
    # node positions (exactly rounded sums, so the order of the features cannot change them) and IDs
    pos_m, pos_raw = [], []
    for idx in clusters:
        pos_m.append((math.fsum(P[idx, 0]) / len(idx), math.fsum(P[idx, 1]) / len(idx)))
        pos_raw.append((math.fsum(R[idx, 0]) / len(idx), math.fsum(R[idx, 1]) / len(idx)))
    ids = [node_id(r, id_style) for r in pos_raw]
    order = sorted(range(len(ids)), key=lambda i: (ids[i], pos_raw[i]))
    seen: dict[str, int] = {}
    final = [""] * len(ids)
    for i in order:
        k = seen.get(ids[i], 0)
        seen[ids[i]] = k + 1
        final[i] = ids[i] if k == 0 else f"{ids[i]}_{k + 1}"
    net = Network(tol_m=tol_m)
    span_of = {}
    for ci, idx in enumerate(clusters):
        span_of[ci] = _span(P[idx]) if len(idx) > 1 else 0.0
        net.nodes[final[ci]] = {"xy": pos_m[ci], "raw": pos_raw[ci], "ends": int(len(idx)), "span_m": span_of[ci]}
    dropped, raw_pipes = [], []
    n_pieces = np.bincount([x[0] for x in pieces], minlength=len(work))
    for k, (q, piece, pxy, praw) in enumerate(pieces):
        a, b = label[2 * k], label[2 * k + 1]
        ln = work[q]
        if piece == 0:
            net.end_node[(ln.feature, ln.part, 0)] = final[a]
        if piece == n_pieces[q] - 1:
            net.end_node[(ln.feature, ln.part, 1)] = final[b]
        if a == b and k in late_loops:
            (h1, _), (h2, _) = late_loops[k][0]
            mid = late_loops[k][1]
            segs = [(a, mid, h1, piece), (mid, a, h2, piece + 0.5)]
        elif a == b:
            dropped.append({"feature": ln.feature, "part": ln.part, "piece": piece,
                            "length_m": float(np.hypot(*(pxy[1:] - pxy[:-1]).T).sum()), "node": final[a]})
            continue
        else:
            segs = [(a, b, pxy, piece)]
        for a_, b_, gxy, pc in segs:
            geom = np.vstack([pos_m[a_], gxy[1:-1], pos_m[b_]]) if len(gxy) > 2 else np.array([pos_m[a_], pos_m[b_]])
            length = float(np.hypot(*(geom[1:] - geom[:-1]).T).sum())
            raw_pipes.append({"start": final[a_], "end": final[b_], "length_m": length, "feature": ln.feature,
                              "part": ln.part, "piece": pc, "xy": [tuple(map(float, g)) for g in geom], "line": q})
    pieces_of = {}
    for p in raw_pipes:
        pieces_of.setdefault((p["feature"], p["part"]), []).append(p)
    parts_of = {}
    for (f, part) in pieces_of:
        parts_of.setdefault(f, set()).add(part)
    for p in raw_pipes:
        p["pieces"] = len(pieces_of[(p["feature"], p["part"])])
    # pipe IDs: the layer's own (with piece and part suffixes) or a hash of the two node IDs and the length
    want = []
    for p in raw_pipes:
        if pipe_ids and p["feature"] in pipe_ids:
            base = pipe_ids[p["feature"]]
            parts = sorted(parts_of[p["feature"]])
            if len(parts) > 1 and p["part"] != parts[0]:
                base += f"_p{parts.index(p['part']) + 1}"
            if p["pieces"] > 1:
                rank = sorted(x["piece"] for x in pieces_of[(p["feature"], p["part"])]).index(p["piece"])
                base += f"_{rank + 1}"
        else:
            a, b = sorted((p["start"], p["end"]))
            base = "P" + hashlib.sha1(f"{a}|{b}|{p['length_m']:.1f}".encode()).hexdigest()[:8]
        want.append(clean_id(base)[0])
    order = sorted(range(len(raw_pipes)), key=lambda i: (want[i], sorted((raw_pipes[i]["start"], raw_pipes[i]["end"])),
                                                         raw_pipes[i]["length_m"], raw_pipes[i]["xy"]))
    used: dict[str, int] = {}
    node_names = set(net.nodes)
    for i in order:
        w = want[i]
        k = used.get(w, 0)
        used[w] = k + 1
        pid = w if k == 0 else clean_id(f"{w}_{k + 1}")[0]
        while pid in node_names:
            pid = clean_id(pid + "_P")[0]
        raw_pipes[i]["id"] = pid
    net.pipes = sorted(raw_pipes, key=lambda p: p["id"])
    # the report
    deg = net.degree()
    names = list(net.nodes)
    ix = {n: i for i, n in enumerate(names)}
    if net.pipes:
        r = [ix[p["start"]] for p in net.pipes]
        c = [ix[p["end"]] for p in net.pipes]
        g = coo_matrix((np.ones(len(r)), (r, c)), shape=(len(names), len(names)))
        n_comp, comp = connected_components(g, directed=False)
        used_nodes = {p["start"] for p in net.pipes} | {p["end"] for p in net.pipes}
        isolated = [n for n in names if n not in used_nodes]
        comp_sizes = np.bincount(comp)
        n_comp -= len(isolated)
    else:
        n_comp, comp_sizes, isolated = 0, np.array([]), names
    pairs: dict[tuple, list] = {}
    for p in net.pipes:
        pairs.setdefault(tuple(sorted((p["start"], p["end"]))), []).append((p["length_m"], (p["feature"], p["part"])))
    parallel = sum(1 for v in pairs.values() if len(v) > 1)
    dup = 0
    for v in pairs.values():        # the two halves of one closed line cut in two are not a line drawn twice
        v = sorted(v)
        dup += sum(1 for (a, fa), (b, fb) in zip(v, v[1:]) if abs(a - b) <= DUPLICATE_LENGTH_TOL * max(b, 1e-9)
                   and fa != fb)
    joined = int(sum(len(idx) for idx in clusters if len(idx) > 1))
    net.report = {"tolerance_m": tol_m, "lines": len(work), "line_ends": int(len(ends)), "pieces": len(pieces),
                  "pipes": len(net.pipes), "nodes": len(net.nodes) - len(isolated), "components": int(n_comp),
                  "largest_component_nodes": int(comp_sizes.max()) if len(comp_sizes) else 0,
                  "dead_ends": int(sum(1 for n, d in deg.items() if d == 1)),
                  "ends_joined": joined, "tees_split": int(n_cuts), "tee_ends": len({e for _, _, e, _ in tees}),
                  "nodes_flagged_wide": len(flagged), "nodes_split": len(split),
                  "split_detail": [{**s, "span_m": round(s["span_m"], 3)} for s in split],
                  "flagged_detail": [{**s, "span_m": round(s["span_m"], 3)} for s in flagged],
                  "dropped_pieces": len(dropped), "dropped_detail": dropped[:50],
                  "repeated_vertices_removed": removed_dupes, "parallel_pairs": parallel,
                  "lines_drawn_twice": dup, "loops_cut": int(loops_cut), "lines_not_finite": not_finite}
    for n in isolated:
        del net.nodes[n]
    t_ft = tol_m / FT
    if split:
        net.flags.append(f"{len(split)} place{'s' if len(split) != 1 else ''} where line ends spread wider than "
                         f"{SPLIT_SPAN:g} times the tolerance (a valve cluster or a pump manifold drawn as short lines, "
                         f"most likely) {'were' if len(split) != 1 else 'was'} split into nodes whose ends are within "
                         f"the tolerance of each other, not joined into one node. Check them on a map.")
    if flagged:
        net.flags.append(f"{len(flagged)} node{'s' if len(flagged) != 1 else ''} joined line ends that are farther "
                         f"apart than the tolerance ({t_ft:.1f} ft), through ends in between: check that they really "
                         f"meet.")
    if loops_cut:
        net.flags.append(f"{loops_cut} closed line{'s' if loops_cut != 1 else ''} (a loop main drawn as one line that "
                         f"ends where it starts, or a loop between two ends the joining distance puts in one node) "
                         f"{'were' if loops_cut != 1 else 'was'} cut in two at the middle of its length, through a new "
                         f"node: EPANET cannot run a pipe from a node to itself.")
    if dropped:
        net.flags.append(f"{len(dropped)} piece{'s' if len(dropped) != 1 else ''} of line had both ends at one node (a "
                         f"piece not much longer than the joining distance, whose two ends it joins, most often) and "
                         f"{'were' if len(dropped) != 1 else 'was'} left out: EPANET cannot run a pipe from a node to "
                         f"itself.")
    if not_finite:
        net.flags.append(f"{not_finite} line{'s' if not_finite != 1 else ''} with a coordinate that is not a number "
                         f"{'were' if not_finite != 1 else 'was'} left out.")
    if dup:
        net.flags.append(f"{dup} pipe{'s' if dup != 1 else ''} run between the same two nodes as another pipe with the "
                         f"same length: a line drawn twice? Both are kept; delete the copy in the layer if it is one.")
    return net


DIAGNOSTIC_KEYS = ("pipes", "nodes", "components", "dead_ends", "ends_joined", "tees_split", "nodes_flagged_wide",
                   "nodes_split", "dropped_pieces")


def diagnostics(lines: list[Line], tol_m: float, id_style: str = "xy", factors=DIAGNOSTIC_FACTORS,
                known: dict | None = None) -> list[dict]:
    """The diagnostic table: the network's counts at each tolerance factor times tol_m.  known: {factor: a report
    already built at that factor} (the build's own, so it is not built twice).  A factor whose ends are too crowded
    to join (TopologyRefused) gives a row of None."""
    out = []
    for f in factors:
        try:
            r = (known or {}).get(f) or build_network(lines, tol_m * f, id_style).report
        except TopologyRefused:
            r = {k: None for k in DIAGNOSTIC_KEYS}
        out.append({"tolerance_m": tol_m * f, "tolerance_ft": tol_m * f / FT, **{k: r[k] for k in DIAGNOSTIC_KEYS}})
    return out


# ----------------------------------------------------------------------------- the skeleton .inp
def _g(v: float) -> str:
    return f"{float(v):.12g}"


def write_inp(net: Network, diam_in: dict, c_of: dict, title: str, elevation_ft: float = 0.0,
              coords_ft=None) -> str:
    """The skeleton as an EPANET .inp in GPM and feet: every node a junction at elevation_ft with no demand, every pipe
    with its length (ft), diameter (in) and Hazen-Williams C, its vertices; no source, tank, pump, valve or control.
    diam_in, c_of: {pipe ID: value}.  coords_ft: a function from (x, y) metres to the coordinates written (feet)."""
    to = coords_ft or (lambda xy: (xy[0] / FT, xy[1] / FT))
    out = ["[TITLE]", title.replace("\n", " ")[:79], "", "[JUNCTIONS]", ";ID              \tElev        \tDemand      \tPattern"]
    for n in sorted(net.nodes):
        out.append(f" {n}\t{_g(elevation_ft)}\t0\t;")
    out += ["", "[RESERVOIRS]", ";ID              \tHead        \tPattern", "", "[TANKS]", "", "[PIPES]",
            ";ID              \tNode1           \tNode2           \tLength      \tDiameter    \tRoughness   \tMinorLoss   \tStatus"]
    for p in net.pipes:
        out.append(f" {p['id']}\t{p['start']}\t{p['end']}\t{_g(p['length_m'] / FT)}\t{_g(diam_in[p['id']])}\t"
                   f"{_g(c_of[p['id']])}\t0\tOpen\t;")
    out += ["", "[PUMPS]", "", "[VALVES]", "", "[DEMANDS]", "", "[STATUS]", "", "[PATTERNS]", "", "[CURVES]", "",
            "[CONTROLS]", "", "[RULES]", "", "[ENERGY]", "", "[EMITTERS]", "", "[QUALITY]", "", "[SOURCES]", "",
            "[REACTIONS]", " Order Bulk            \t1", " Order Tank            \t1", " Order Wall            \t1",
            " Global Bulk           \t0", " Global Wall           \t0", "", "[MIXING]", "",
            "[TIMES]", " Duration           \t168:00", " Hydraulic Timestep \t1:00", " Quality Timestep   \t0:05",
            " Pattern Timestep   \t1:00", " Report Timestep    \t1:00", " Start ClockTime    \t12 am", "",
            "[REPORT]", " Status             \tNo", " Summary            \tNo", "",
            "[OPTIONS]", " Units              \tGPM", " Headloss           \tH-W", " Specific Gravity   \t1.0",
            " Viscosity          \t1.0", " Trials             \t40", " Accuracy           \t0.001",
            " Unbalanced         \tContinue 10", " Pattern            \t1", " Demand Multiplier  \t1.0",
            " Emitter Exponent   \t0.5", " Quality            \tNone mg/L", " Diffusivity        \t1.0",
            " Tolerance          \t0.01", "", "[COORDINATES]", ";Node            \tX-Coord         \tY-Coord"]
    for n in sorted(net.nodes):
        x, y = to(net.nodes[n]["xy"])
        out.append(f" {n}\t{_g(x)}\t{_g(y)}")
    out += ["", "[VERTICES]", ";Link            \tX-Coord         \tY-Coord"]
    for p in net.pipes:
        for v in p["xy"][1:-1]:
            x, y = to(v)
            out.append(f" {p['id']}\t{_g(x)}\t{_g(y)}")
    out += ["", "[END]", ""]
    return "\n".join(out)
