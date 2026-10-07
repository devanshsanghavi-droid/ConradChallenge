"""
msx.py: a shared EPANET-MSX runner (iteration 4; task 12's chloramine truths, and task 13's audit), and task 13's
two-reactant chlorine-organics chemistry (2RA: the batch solution, the scale s, the MSX model and the hidden truth on the
committed draws; at the end of this file).

EPANET-MSX (multi-species extension) ships inside WNTR 1.5.  On this Mac it runs, but only with the workarounds the
planning probes found, all of which live here so no caller has to remember them:

  * WNTR's darwin libepanetmsx.dylib links @rpath/libepanet2.dylib and carries no LC_RPATH, and WNTR only extends
    DYLD_FALLBACK_LIBRARY_PATH when that variable is already set.  preload() loads the sibling libepanet2.dylib with
    RTLD_GLOBAL first, so dyld resolves the dependency by its install name.
  * The MSX dylib also links Homebrew's libomp (/opt/homebrew/opt/libomp/lib/libomp.dylib).  A Mac without it fails
    with an unreadable dlopen error, so preload() checks for the file first and raises a clear error.
  * OpenMP threads: one per process (OMP_NUM_THREADS=1, set before the library loads), so a pool of workers does not
    oversubscribe the cores.
  * EPANET-MSX leaves scratch files (msxXXXXXX, enXXXXXX) in the CURRENT directory.  Every run happens inside its own
    temporary working directory, which is deleted afterwards; the process's working directory is restored.  The
    change of directory is process-wide: a module lock lets one MSX run at a time per process, and the app runs its
    demo truth in a separate process (chloramine.truth_in_subprocess), so its own threads never see the change.
  * A CONCEN source at a RESERVOIR silently gives zeros under MSX; add_sources() uses SETPOINT at reservoirs and
    tanks and CONCEN at inflow junctions.
  * WNTR writes per-tank PARAMETER values with the PIPE keyword (MSX error 405), so only per-pipe parameters are set
    (set_pipe_parameter); a tank-wide value must be a CONSTANT.
  * MSX names are case-insensitive and some are reserved (a constant 'kc' collides with the hydraulic variable Kc):
    check_names() refuses the reserved ones before a run.
  * WNTR's MSX binary reader names result columns in wn.node_name_list order, but the file holds them in EPANET's
    index order (junctions, reservoirs, tanks); for a network built in code they differ, so run() renames them.
  * COMPILER GC (the reactions compiled to C by `gcc`; on a Mac, Xcode or its command-line tools) is much faster and
    gave identical results in the probes (about 12 times faster on the saved fallback check's 2-pipe chain here, 1.2 s
    against 14.5 s; the plan's probe measured 15 to 18 s against 328 s for a network run).  run() uses it when compiler_available() finds a compiler and,
    if none is found or the compiled run fails, says so loudly on stderr, records it in the run information and runs
    with COMPILER NONE.  A failed MSXopen leaves MSX's global project open (the retry would fail with MSX error 520),
    so every failed run closes the MSX and EPANET projects first (_close_project).

EPANET's own first-order wall term, k_w k_f / (k_w + k_f) x 4/D with the Sherwood-number mass-transfer coefficient
k_f (EPANET 2.2 manual, Rossman 2000), is available as MSX terms (add_mass_transfer_terms), in the network's own length
units, with the file's relative viscosity and diffusivity.
"""
from __future__ import annotations

import ctypes
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass

import pandas as pd
import wntr

LIBOMP = "/opt/homebrew/opt/libomp/lib/libomp.dylib"
FT = 0.3048
US_FLOW_UNITS = ("CFS", "GPM", "MGD", "IMGD", "AFD")
CHLORINE_DIFFUSIVITY_FT2_S = 1.3e-8    # EPANET's reference molecular diffusivity (chlorine in water, 20 C)
WATER_VISCOSITY_FT2_S = 1.1e-5         # EPANET's reference kinematic viscosity of water (20 C)
RESERVED = {"d", "q", "u", "re", "us", "ff", "av", "kc", "len",     # MSX hydraulic variables
            "abs", "sgn", "sqrt", "log", "exp", "sin", "cos", "tan", "cot", "asin", "acos", "atan", "acot",
            "sinh", "cosh", "tanh", "coth", "log10", "step"}
_PRELOADED = False
_LOCK = threading.Lock()        # os.chdir is process-wide: one MSX run at a time per process


def preload() -> None:
    """Make WNTR's MSX library loadable (idempotent).  Raises RuntimeError when Homebrew's libomp is missing."""
    global _PRELOADED
    if _PRELOADED:
        return
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    if platform.system() == "Darwin":
        if not os.path.exists(LIBOMP):
            raise RuntimeError(f"EPANET-MSX needs Homebrew's OpenMP runtime at {LIBOMP}, which is missing. Install it with "
                               "'brew install libomp' (WNTR's MSX library links that exact path).")
        arch = "darwin-arm" if platform.machine().lower() in ("arm64", "aarch64") else "darwin-x64"
        d = os.path.join(os.path.dirname(wntr.__file__), "epanet", "libepanet", arch)
        ctypes.CDLL(os.path.join(d, "libepanet2.dylib"), mode=ctypes.RTLD_GLOBAL)
    _PRELOADED = True


def msx_library_path() -> str:
    """The MSX shared library WNTR loads (the same file, so ctypes returns the same loaded library)."""
    d = os.path.join(os.path.dirname(wntr.__file__), "epanet", "libepanet")
    if os.name in ("nt", "dos"):
        return os.path.join(d, "windows-x64", "epanetmsx.dll")
    if sys.platform == "darwin":
        return os.path.join(d, "darwin-arm" if "arm" in platform.platform().lower() else "darwin-x64", "libepanetmsx.dylib")
    return os.path.join(d, "linux-x64", "libepanetmsx.so")


def compiler_available() -> bool:
    """True when COMPILER GC can work: MSX compiles with `gcc` from the PATH, and on a Mac /usr/bin/gcc is only a stub
    that asks to install the command-line tools unless Xcode or those tools are installed (xcode-select -p)."""
    if shutil.which("gcc") is None:
        return False
    if platform.system() == "Darwin":
        xs = "/usr/bin/xcode-select"
        if not os.path.exists(xs):
            return False
        try:
            return subprocess.run([xs, "-p"], capture_output=True, timeout=20).returncode == 0
        except (OSError, subprocess.SubprocessError):
            return False
    return True


def _close_project() -> None:
    """Close MSX's global project and its EPANET project after a failed run (both calls are harmless when nothing is
    open).  Without this a failed MSXopen leaves the project open and every later run in the process fails with MSX
    error 520, the compiled-to-uncompiled fallback included."""
    try:
        lib = ctypes.CDLL(msx_library_path())
    except OSError:
        return
    for fn in ("MSXclose", "MSXENclose"):
        try:
            getattr(lib, fn)()
        except Exception:  # noqa: BLE001
            pass


def us_units(wn) -> bool:
    """True when the .inp's flow units are US customary (lengths and diameters in feet inside MSX expressions)."""
    return str(wn.options.hydraulic.inpfile_units).upper() in US_FLOW_UNITS


def metres_per_length_unit(wn) -> float:
    return FT if us_units(wn) else 1.0


def check_names(msx) -> None:
    """Refuse species, constants, parameters and terms whose names MSX reserves (case-insensitive)."""
    names = list(msx.species_name_list) + list(msx.constant_name_list) + list(msx.parameter_name_list) + list(msx.term_name_list)
    bad = sorted(n for n in names if n.lower() in RESERVED)
    if bad:
        raise ValueError(f"MSX names that collide with reserved words (case-insensitive): {bad}")
    low = [n.lower() for n in names]
    dup = sorted({n for n in low if low.count(n) > 1})
    if dup:
        raise ValueError(f"MSX names that differ only in case: {dup}")


def add_mass_transfer_terms(msx, wn) -> None:
    """EPANET 2.2's wall mass-transfer coefficient as MSX terms, Kf in the network's length unit per second
    (rate units SEC): Sh = 2 (Re < 1); 3.65 + 0.0668 y / (1 + 0.04 y^0.667), y = (D / Len) Re Sc (laminar);
    0.0149 Re^0.88 Sc^0.333 (Re >= 2300); Kf = Sh Dm / D.  Dm and the viscosity are EPANET's reference values times the
    file's own relative DIFFUSIVITY and VISCOSITY.  The wall rate of a species C with wall coefficient kwp (same units
    as Kf) is then (4 / D) kwp Kf / (kwp + Kf) C, EPANET's first-order wall reaction.
    MSX computes its hydraulic variable Re with water's 20 C viscosity whatever the file's VISCOSITY says, while EPANET
    uses the file's.  Found by task 13's saved check at 10 C (a 1.30 viscosity ratio moved the wall term enough for a
    0.086 mg/L difference on Net3; the check reruns the uncorrected runner and records it as before_fix_max_abs_diff_mgL
    in outputs/chem/checks_report.json): for a file whose relative VISCOSITY is not 1, Re is divided by it here (a term
    ReV), so both engines see the same Reynolds number.  Every MSX run before task 13 used a file at relative viscosity
    1 (Net3, Net2 and ky4 all are), where the expressions are written exactly as before."""
    unit2 = metres_per_length_unit(wn) ** 2 / FT ** 2          # ft2 -> (length unit)2
    dm = CHLORINE_DIFFUSIVITY_FT2_S * float(wn.options.quality.diffusivity) / unit2
    visc_rel = float(wn.options.hydraulic.viscosity)
    nu = WATER_VISCOSITY_FT2_S * visc_rel / unit2
    msx.add_constant("Dm", dm, note="molecular diffusivity, length unit^2 per s (EPANET reference x file DIFFUSIVITY)")
    msx.add_constant("Sc", nu / dm, note="Schmidt number")
    re = "Re"
    if visc_rel != 1.0:
        msx.add_constant("Vrel", visc_rel, note="the file's relative VISCOSITY (MSX's own Re assumes 1)")
        msx.add_term("ReV", "Re/Vrel")
        re = "ReV"
    msx.add_term("Yg", f"D/Len*{re}*Sc")
    msx.add_term("Sh", f"step(1-{re})*2 + step({re}-1)*step(2300-{re})*(3.65+0.0668*Yg/(1+0.04*Yg^0.667))"
                       f" + step({re}-2300)*0.0149*{re}^0.88*Sc^0.333")
    msx.add_term("Kf", "Sh*Dm/D")


def wall_rate_unit(wn, m_per_day: float) -> float:
    """A wall coefficient in m/day as MSX needs it here: the network's length unit per second."""
    return float(m_per_day) / metres_per_length_unit(wn) / 86400.0


def add_sources(msx, wn, values: dict[str, dict[str, float]]) -> None:
    """values: {node: {species: concentration}}.  SETPOINT at reservoirs and tanks (a CONCEN source at a reservoir
    gives zeros under MSX), CONCEN at junctions (an inflow junction modelled as a negative demand, as in Net2)."""
    fixed = set(wn.reservoir_name_list) | set(wn.tank_name_list)
    for node, sp in values.items():
        kind = "SETPOINT" if node in fixed else "CONCEN"
        for name, v in sp.items():
            msx.add_source(kind, name, node, float(v))


def set_pipe_parameter(msx, name: str, values: dict[str, float]) -> None:
    """Per-pipe values of an MSX PARAMETER (per-tank values are not written correctly by WNTR: use a CONSTANT)."""
    pv = msx.network_data.parameter_values[name].pipe_values
    for pipe, v in values.items():
        pv[pipe] = float(v)


def epanet_order(wn) -> tuple[list[str], list[str]]:
    """Node and link names in EPANET's own index order for the .inp WNTR writes: junctions, then reservoirs, then
    tanks; pipes, then pumps, then valves."""
    return (list(wn.junction_name_list) + list(wn.reservoir_name_list) + list(wn.tank_name_list),
            list(wn.pipe_name_list) + list(wn.pump_name_list) + list(wn.valve_name_list))


def _relabel(res, wn, species) -> None:
    """WNTR 1.5's MSX binary reader names the columns with wn.node_name_list (and link_name_list), but the file holds
    them in EPANET's index order.  The two agree for a network read from an .inp (Net3, Net2), not for one built in
    code (a reservoir added first comes first in node_name_list and last in EPANET's order), so the columns are renamed
    to EPANET's order here."""
    nodes, links = epanet_order(wn)
    for sp in species:
        if sp in res.node and len(res.node[sp].columns) == len(nodes):
            res.node[sp].columns = pd.Index(nodes, name=res.node[sp].columns.name)
        if sp in res.link and len(res.link[sp].columns) == len(links):
            res.link[sp].columns = pd.Index(links, name=res.link[sp].columns.name)


def _run_once(wn, msx, compiler: str):
    msx.options.compiler = compiler
    wn.msx = msx
    with _LOCK:
        old = os.getcwd()
        with tempfile.TemporaryDirectory(prefix="rm_msx_") as tmp:
            os.chdir(tmp)                           # MSX writes msxXXXXXX / enXXXXXX scratch files to the cwd
            try:
                res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=os.path.join(tmp, "run"))
            except BaseException:
                _close_project()                    # a failed run leaves MSX's project open (error 520 next time)
                raise
            finally:
                os.chdir(old)
    _relabel(res, wn, list(msx.species_name_list))
    return res


def run(wn, msx, compiler: str = "GC", fallback: bool = True):
    """One EPANET-MSX run of `wn` with the reaction model `msx`, in its own temporary working directory.
    Returns (results, info): WNTR results (results.node[<species>] is a time x node DataFrame of concentrations in the
    species' units) and {'compiler': the compiler that ran, 'seconds': wall time, 'fallback_reason': None or why GC
    failed}.  With compiler 'GC' and fallback=True, a machine with no compiler (compiler_available) or a failed compiled
    run is reported on stderr and run with 'NONE'; with fallback=False either raises."""
    if compiler not in ("GC", "NONE"):
        raise ValueError("compiler must be 'GC' or 'NONE'")
    preload()
    check_names(msx)
    quality0 = wn.options.quality.parameter
    wn.options.quality.parameter = "NONE"          # MSX carries the chemistry; EPANET's own quality is off
    t0 = time.perf_counter()
    try:
        if compiler == "GC" and not compiler_available():
            reason = "no C compiler: gcc is not on the PATH, or (on a Mac) Xcode's command-line tools are not installed"
            if not fallback:
                raise RuntimeError(f"EPANET-MSX COMPILER GC needs a C compiler: {reason}")
            print(f"\n*** EPANET-MSX: {reason}; running with COMPILER NONE, which is much slower ***\n",
                  file=sys.stderr, flush=True)
            res = _run_once(wn, msx, "NONE")
            return res, {"compiler": "NONE", "seconds": time.perf_counter() - t0, "fallback_reason": reason}
        try:
            res = _run_once(wn, msx, compiler)
            return res, {"compiler": compiler, "seconds": time.perf_counter() - t0, "fallback_reason": None}
        except Exception as e:  # noqa: BLE001
            if compiler == "NONE" or not fallback:
                raise
            reason = f"{type(e).__name__}: {e}"
            print(f"\n*** EPANET-MSX with COMPILER GC failed ({reason[:300]}); rerunning with COMPILER NONE, which is "
                  f"much slower ***\n", file=sys.stderr, flush=True)
            t0 = time.perf_counter()
            res = _run_once(wn, msx, "NONE")
            return res, {"compiler": "NONE", "seconds": time.perf_counter() - t0, "fallback_reason": reason}
    finally:
        wn.msx = None                              # later EPANET runs of this network must not run MSX again
        wn.options.quality.parameter = quality0


# ===================================================================== task 13: two-reactant chlorine-organics (2RA)
# Fisher, Kastl & Sathasivan 2012 (Water Res 46:3293, doi:10.1016/j.watres.2012.03.017) and Fisher et al. 2011
# (doi:10.1016/j.watres.2011.06.032): free chlorine C reacts with a fast (F) and a slow (S) pool of organic reactants,
#     dC/dt = -kF C F - kS C S - wall,   dF/dt = -kF C F,   dS/dt = -kS C S.
# The rate constants and the initial reactant split are the Greenvale water's (Fisher 2012, quoted through Walski's
# Bentley blog: a SECONDARY source).  Every number below is at 20 C; task 13 runs the truth at fixed conditions only
# (the plan's seasonal 2ra_warm variant was dropped with the temperature models, addendum 2).
KF20_L_PER_MG_H = 0.141          # fast reactant, L/mg/h at 20 C
KS20_L_PER_MG_H = 0.00366        # slow reactant, L/mg/h at 20 C
F0_GREENVALE_MGL = 1.13          # the water's initial fast and slow reactant concentrations (mg/L as chlorine demand);
S0_GREENVALE_MGL = 2.87          # the truth uses s x these, with s matched to the network's committed bulk rate
TWO_RA_C0_MGL = 1.2              # the dose at which s is matched (the committed nominal dose)
TWO_RA_MATCH_H = 24.0            # the age at which the apparent first-order rate equals the committed kb (the plan's
                                 # default; a 96 h match is the plan's Net2 sensitivity)
TWO_RA_SPECIES = ("CL2", "FAST", "SLOW")
TWO_RA_SPECIES_TOL = (1e-8, 1e-6)   # MSX absolute (mg/L) and relative tolerances per species
TWO_RA_SOLVER = "RK5"


@dataclass(frozen=True)
class TwoReactantTruth:
    """Settings of a 2RA truth beyond the committed draws (simulate.build_scenario(two_reactant=...)).
    match_h  : the age (hours) at which s matches the network's committed bulk rate (24, the plan's default; 96 is the
               plan's Net2 sensitivity)
    compiler : MSX COMPILER, 'GC' with run()'s loud fallback to 'NONE'"""
    match_h: float = TWO_RA_MATCH_H
    compiler: str = "GC"

    def __post_init__(self):
        object.__setattr__(self, "match_h", float(self.match_h))
        if not (1.0 <= self.match_h <= 168.0):
            raise ValueError(f"match_h {self.match_h} is outside 1 to 168 h")
        if self.compiler not in ("GC", "NONE"):
            raise ValueError("compiler must be 'GC' or 'NONE'")


def two_reactant_batch(c0: float, f0: float, s0: float, t_h, kf_per_h: float = KF20_L_PER_MG_H,
                       ks_per_h: float = KS20_L_PER_MG_H, kbulk_per_h: float = 0.0):
    """The 2RA batch (plug-flow) solution, scipy LSODA at rtol 1e-10: returns (C, F, S) at the ages t_h (hours), each
    an array.  kbulk_per_h adds a first-order chlorine loss (0 for the truth; the saved checks use it)."""
    import numpy as np
    from scipy.integrate import solve_ivp
    t = np.atleast_1d(np.asarray(t_h, dtype=float))

    def rhs(_, y):
        c, f, s = y
        return [-(kf_per_h * c * f + ks_per_h * c * s + kbulk_per_h * c), -kf_per_h * c * f, -ks_per_h * c * s]
    tt = np.unique(np.concatenate([[0.0], t]))
    sol = solve_ivp(rhs, (0.0, float(tt.max())), [float(c0), float(f0), float(s0)], t_eval=tt, rtol=1e-10, atol=1e-13,
                    method="LSODA")
    if not sol.success:
        raise RuntimeError(f"2RA batch solve failed: {sol.message}")
    idx = np.searchsorted(tt, t)
    return sol.y[0][idx], sol.y[1][idx], sol.y[2][idx]


def two_reactant_apparent_rate(s: float, age_h: float, c0: float = TWO_RA_C0_MGL) -> float:
    """-ln(C(age) / C0) / age, per day: the first-order rate that would give the same chlorine at that age."""
    import math
    c = float(two_reactant_batch(c0, F0_GREENVALE_MGL * s, S0_GREENVALE_MGL * s, [age_h])[0][0])
    return -math.log(c / c0) / float(age_h) * 24.0


def two_reactant_scale(kb_per_day: float, match_h: float = TWO_RA_MATCH_H, c0: float = TWO_RA_C0_MGL) -> float:
    """s such that the 2RA water's apparent first-order rate at age match_h, dose c0 and 20 C equals kb_per_day
    (brentq on s in (1e-6, 2): above about 2 the Net3 chlorine is gone before 96 h).  Truth tuning, disclosed: it
    makes the 2RA truth agree with the committed first-order truth's bulk rate at that one age."""
    from scipy.optimize import brentq
    return float(brentq(lambda s: two_reactant_apparent_rate(s, match_h, c0) - float(kb_per_day), 1e-6, 2.0,
                        xtol=1e-12, rtol=1e-12))


def two_reactant_model(wn, kw_pipe_m_day: dict | None, kf_per_h: float = KF20_L_PER_MG_H,
                       ks_per_h: float = KS20_L_PER_MG_H):
    """The 2RA reactions as an MSX model for `wn` (sources are added by the caller with add_sources):
    species CL2, FAST, SLOW in mg/L; rates per second; RK5 at the network's quality step.  kw_pipe_m_day: per-pipe
    first-order chlorine wall coefficient (m/day) through EPANET's own mass-transfer-limited form
    (add_mass_transfer_terms), as in the committed truth; None: no wall (the saved chain check).  Tanks react in the
    bulk only, as in EPANET.  The constants are named kFast and kSlow: 'kF' would collide with the term Kf (MSX names
    are case-insensitive)."""
    from wntr.msx import MsxModel
    m = MsxModel()
    m.options.rate_units, m.options.area_units, m.options.solver = "SEC", "M2", TWO_RA_SOLVER
    m.options.coupling = "NONE"
    m.options.timestep = int(wn.options.time.quality_timestep)
    for sp in TWO_RA_SPECIES:
        m.add_species(sp, "bulk", units="MG", atol=TWO_RA_SPECIES_TOL[0], rtol=TWO_RA_SPECIES_TOL[1])
    m.add_constant("kFast", float(kf_per_h) / 3600.0, note="fast reactant, L/mg/s")
    m.add_constant("kSlow", float(ks_per_h) / 3600.0, note="slow reactant, L/mg/s")
    bulk = "-kFast*CL2*FAST - kSlow*CL2*SLOW"
    pipe = bulk
    if kw_pipe_m_day is not None:
        add_mass_transfer_terms(m, wn)
        m.add_parameter("kwp", 0.0)
        set_pipe_parameter(m, "kwp", {p: wall_rate_unit(wn, v) for p, v in kw_pipe_m_day.items()})
        pipe = bulk + " - (4/D)*kwp*Kf/(kwp+Kf)*CL2"
    for where, expr in (("pipe", pipe), ("tank", bulk)):
        m.add_reaction("CL2", where, "rate", expr)
        m.add_reaction("FAST", where, "rate", "-kFast*CL2*FAST")
        m.add_reaction("SLOW", where, "rate", "-kSlow*CL2*SLOW")
    return m


def two_reactant_truth(wn, seed: int, rng, rng_m, source_dose: float, kb_per_day: float, kw_m_per_day: float,
                       match_h: float = TWO_RA_MATCH_H, compiler: str = "GC"):
    """The hidden 2RA truth (task 13; called by simulate._chem_truth for kinetics '2ra').  The committed draws are
    consumed exactly as in build_scenario's default branch: the monthly bulk factor u ~ U(0.8, 1.2), per pipe the wall
    factor exp(N(0, 0.4)) (with the file's own roughness) and the roughness error exp(N(0, 0.10)), the global and
    per-node demand, and one dose U(0.9, 1.1) per source; no new draw is made.
      * reactant loads at every source: FAST = 1.13 s u, SLOW = 2.87 s u mg/L, with s = two_reactant_scale(kb_per_day,
        match_h) (the plan's F0 = 1.13 s (TOC/2) u at TOC_ref; proportional to TOC and to u is an ASSUMPTION);
      * chlorine at every source: its drawn dose; zero initial quality everywhere (the committed truth's file values
        are 1/1200 of a dose or zero);
      * wall: the committed per-pipe coefficient kw 2^(-(C - 130)/30) exp(N(0, 0.4)), EPANET's mass-transfer-limited
        first-order form, the file's viscosity and diffusivity (20 C);
      * the network's quality step (300 s), RK5, compiled reactions (COMPILER GC) with msx.run's loud fallback.
    Returns (chlorine, time x node in mg/L; info)."""
    import numpy as np
    from .simulate import roughness_factor, source_nodes
    u = float(rng_m.uniform(0.8, 1.2))
    kw_pipe = {}
    for pn, pipe in wn.pipes():
        kw_pipe[pn] = kw_m_per_day * roughness_factor(pipe.roughness, 1.0) * np.exp(rng.normal(0, 0.4))
        pipe.roughness = pipe.roughness * np.exp(rng.normal(0, 0.10))
    global_mult = rng_m.uniform(0.85, 1.15)
    for _, j in wn.junctions():
        for ts in j.demand_timeseries_list:
            ts.base_value = ts.base_value * global_mult * np.exp(rng_m.normal(0.0, 0.15))
    doses = {s: source_dose * rng_m.uniform(0.9, 1.1) for s in source_nodes(wn)}
    s = two_reactant_scale(kb_per_day, match_h)
    f0, s0 = F0_GREENVALE_MGL * s * u, S0_GREENVALE_MGL * s * u
    m = two_reactant_model(wn, kw_pipe)
    add_sources(m, wn, {node: {"CL2": d, "FAST": f0, "SLOW": s0} for node, d in doses.items()})
    res, run_info = run(wn, m, compiler=compiler)
    q = res.node["CL2"]
    info = {"disinfectant": "free_chlorine", "species": "free chlorine", "kinetics": "2ra", "temp_C": 20.0,
            "toc_mgL": None, "kb_per_day_matched": float(kb_per_day), "match_h": float(match_h), "scale_s": s,
            "bulk_month_factor": u, "fast_mgL": f0, "slow_mgL": s0,
            "kfast_L_per_mg_h": KF20_L_PER_MG_H, "kslow_L_per_mg_h": KS20_L_PER_MG_H,
            "source_doses_mgL": {k: float(v) for k, v in doses.items()}, "quality_scale": 1.0,
            "msx_compiler": run_info["compiler"], "msx_fallback_reason": run_info["fallback_reason"],
            "_volatile": {"msx_seconds": run_info["seconds"]}}
    return q, info
