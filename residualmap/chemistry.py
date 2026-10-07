"""
chemistry.py: the water-chemistry layer added in iteration 4 (journal tasks 8 to 14).

Nothing in this module is used by the default path (free chlorine, first-order decay, the file's own water
properties).  It holds:

  * temperature physics: the Arrhenius multiplier relative to 20 C and the relative viscosity and molecular
    diffusivity of water (EPANET's VISCOSITY and DIFFUSIVITY options, which move only the wall mass-transfer
    coefficient on a Hazen-Williams network);
  * the seasonal schedules used by the synthetic experiments.  They are ASSUMPTIONS, not measurements: no
    coastal California plant series of water temperature or TOC was available when they were written;
  * `Chemistry`, the description of one chemistry condition (disinfectant, kinetics, temperature, TOC,
    threshold, E/R, wall mode), and `cache_tag`, a short hash that keeps grids built under different
    conditions in different cache files;
  * a truth-only helper for EPANET's nonlinear kinetics (order 0 and 2, limiting concentration).  WNTR takes
    concentrations in kg/m3 and does not convert bulk coefficients of order other than 1, so the repo's
    habit of passing the dose as 1.2 (written to the .inp as 1200 mg/L) is harmless for first order but breaks
    every nonlinear option.  `use_mg_per_litre` applies the verified recipe;
  * the closed-form solutions the saved checks compare EPANET against.

Sources for the numbers are given where they are used.  Free chlorine and chloramine are never mixed: a
`Chemistry` carries one disinfectant, its measured species and its own threshold.
"""
from __future__ import annotations

import glob
import hashlib
import json
import math
import os
from dataclasses import asdict, dataclass, replace

DAY = 86400.0
KELVIN = 273.15
TREF_C = 20.0                      # every calibrated rate in the repo is read as a 20 C value
TOC_REF_MGL = 2.0                  # TOC at which the bulk rate is the calibrated one (task 11)
ER_DEFAULT_K = 8000.0              # prior centre of Cejas, Diaz & Gonzalez 2026 (8000 +/- 2000 K)
ER_HYPOTHESES_K = (5000.0, 8000.0, 12000.0)   # low end of Powell et al. 2000; Cejas 2026; about Blokker et al. 2014
ER_TRUTH_RANGE_K = (4660.0, 12104.0)          # compiled bulk span quoted by Cejas 2026, used to draw hidden truths
THETA_W_TRUTH_RANGE = (1.00, 1.07)            # ASSUMPTION: direction only (Lee et al. 2014); no published wall theta

FREE_CHLORINE, CHLORAMINE = "free_chlorine", "chloramine"
DISINFECTANTS = (FREE_CHLORINE, CHLORAMINE)
SPECIES = {FREE_CHLORINE: "free chlorine", CHLORAMINE: "total chlorine"}
DEFAULT_THRESHOLD_MGL = {FREE_CHLORINE: 0.2, CHLORAMINE: 0.5}
THRESHOLD_NOTE = {
    FREE_CHLORINE: "0.2 mg/L free chlorine: a common operating minimum, not a California rule "
                   "(California requires a detectable residual)",
    CHLORAMINE: "0.5 mg/L total chlorine: a common utility operating target, not a California rule "
                "(California requires a detectable residual)",
}
KINETICS = ("first", "first_si", "order2", "clark")   # first = the repo's legacy-unit first order
WALL_MODES = ("arrhenius", "mass_transfer_only")      # mass_transfer_only = ablation M1b (wall chemistry not scaled)
TRUTH_SEED_OFFSET = {FREE_CHLORINE: 20_000, CHLORAMINE: 40_000}   # same pattern as default_rng(10_000 + seed)
BUILT_DISINFECTANTS = (FREE_CHLORINE,)   # chloramine decay physics arrive with task 12: until then a chloramine
                                         # truth or grid is refused rather than run with free-chlorine physics
CHEM_CACHE_VERSION = 1   # part of every cache_tag: bump it whenever the physics a condition runs changes
                         # (simulate_nominal_chlorine's condition keywords, Chemistry.sim_kwargs), so no grid
                         # pickled under the old physics is ever served for the new one


# ----------------------------------------------------------------------------- temperature physics
def arrhenius(temp_C: float, er_K: float, tref_C: float = TREF_C) -> float:
    """Rate multiplier relative to tref: f(T; E/R) = exp[(E/R)(T - Tref) / ((Tref + 273.15)(T + 273.15))].
    Fisher, Kastl & Sathasivan 2012 (Water Res 46:3293); Cejas, Diaz & Gonzalez 2026 (Water 18:1390).
    f(10; 8000) = 0.381, f(25; 8000) = 1.580, f(20; any) = 1 exactly."""
    return math.exp(er_K * (temp_C - tref_C) / ((tref_C + KELVIN) * (temp_C + KELVIN)))


def kinematic_viscosity(temp_C: float) -> float:
    """nu(T) = 497e-6 / (T + 42.5)^1.5 m2/s (Blokker, Vreeburg & Speight 2014, Procedia Eng 70:172, Table 1)."""
    return 497e-6 / (temp_C + 42.5) ** 1.5


def viscosity_ratio(temp_C: float, tref_C: float = TREF_C) -> float:
    """nu(T) / nu(Tref): 1.299 at 10 C, 0.891 at 25 C, 1 exactly at 20 C."""
    return kinematic_viscosity(temp_C) / kinematic_viscosity(tref_C)


def diffusivity_ratio(temp_C: float, tref_C: float = TREF_C) -> float:
    """D(T) / D(Tref) with D(T) = 1.21e-9 ((T + 273.15)/293.15) nu(20)/nu(T) (Blokker 2014, Table 1):
    0.744 at 10 C, 1.142 at 25 C, 1 exactly at 20 C."""
    return ((temp_C + KELVIN) / (tref_C + KELVIN)) * (kinematic_viscosity(tref_C) / kinematic_viscosity(temp_C))


def apply_water_temperature(wn, temp_C: float | None) -> tuple[float, float]:
    """Set EPANET's relative VISCOSITY and DIFFUSIVITY for water at temp_C by MULTIPLYING the file's own values
    (ky4.inp sets DIFFUSIVITY 0.05).  None leaves the file untouched.  With Hazen-Williams headloss (Net3, Net2,
    ky4) the hydraulics do not change; only kf = Sh D / d inside EPANET's wall term moves.  Returns the ratios."""
    if temp_C is None:
        return 1.0, 1.0
    v, d = viscosity_ratio(temp_C), diffusivity_ratio(temp_C)
    wn.options.hydraulic.viscosity = wn.options.hydraulic.viscosity * v
    wn.options.quality.diffusivity = wn.options.quality.diffusivity * d
    return v, d


# ----------------------------------------------------------------------------- seasonal schedules (ASSUMPTIONS)
def monthly_temperature(month: int) -> float:
    """ASSUMPTION, coastal plant water: T_m = 15 - 5 cos(2 pi (m - 2)/12), rounded to 0.5 C.  Jan to Dec:
    10.5, 10, 10.5, 12.5, 15, 17.5, 19.5, 20, 19.5, 17.5, 15, 12.5 (7 distinct values).  Inside the 3.5 to
    28 C range of Fisher 2012 and below the 30 C limit of Cejas 2026."""
    return round(2.0 * (15.0 - 5.0 * math.cos(2.0 * math.pi * (month - 2) / 12.0))) / 2.0


def soil_temperature(month: int) -> float:
    """ASSUMPTION: T_soil,m = 16 + 6 sin(2 pi (m - 5)/12) C, the temperature water warms toward in the pipes."""
    return 16.0 + 6.0 * math.sin(2.0 * math.pi * (month - 5) / 12.0)


WARMING_TAU_PIPE_H = 12.0          # ASSUMPTION: e-folding time of a pipe's water toward soil temperature (Blokker &
                                   # Pieterse-Quirijns 2013, JAWWA 105(1):E19, say only that heating is faster than
                                   # residence)
WARMING_TAU_TANK_H = 72.0          # ASSUMPTION: the same for a tank (much less wall contact per volume)


@dataclass(frozen=True)
class Warming:
    """In-network warming, for hidden truths only (task 10, variant V3): water leaves the plant at the plant temperature
    T_m and moves toward the soil temperature T_soil as it ages.  Per pipe, T_p = T_soil + (T_m - T_soil)
    exp(-a_p / tau_pipe), with a_p the nominal daily-mean water age at the pipe's downstream node; per tank the same
    with the tank's own age and tau_tank.  The model never sees this: it is told only the plant temperature.  Context:
    Blokker et al. 2014 found that water entering at 10 C and warming toward 25 C soil moved the share of customers
    below 0.2 mg/L from 0.4% to 33%."""
    soil_temp_C: float
    tau_pipe_h: float = WARMING_TAU_PIPE_H
    tau_tank_h: float = WARMING_TAU_TANK_H

    def __post_init__(self):
        object.__setattr__(self, "soil_temp_C", float(self.soil_temp_C))
        if not (0.0 <= self.soil_temp_C <= 35.0):
            raise ValueError(f"soil_temp_C {self.soil_temp_C} is outside 0 to 35 C")
        if not (self.tau_pipe_h > 0 and self.tau_tank_h > 0):
            raise ValueError("warming time constants must be positive")

    def temperature(self, plant_C: float, age_h: float, tau_h: float) -> float:
        """T_soil + (T_plant - T_soil) exp(-age / tau): exactly T_plant when the soil is at the plant temperature."""
        return self.soil_temp_C + (float(plant_C) - self.soil_temp_C) * math.exp(-max(float(age_h), 0.0) / tau_h)


TOC_SCHEDULE_MGL = {1: 2.5, 2: 2.5, 3: 2.5, 4: 2.0, 5: 1.5, 6: 1.5, 7: 1.5, 8: 1.5, 9: 1.5, 10: 1.5,
                    11: 3.0, 12: 2.5}


def monthly_toc(month: int) -> float:
    """ASSUMPTION, creek-fed surface water with wet-season organic matter: 1.5 mg/L May to Oct, 3.0 in Nov
    (first storms), 2.5 Dec to Mar, 2.0 in Apr; inside Powell et al. 2000's 1 to 3 mg/L."""
    return TOC_SCHEDULE_MGL[int(month)]


def toc_ratio(toc_mgL: float | None) -> float:
    """TOC / TOC_ref (AWWARF / Powell form k proportional to TOC); 1 when TOC is not logged."""
    return 1.0 if toc_mgL is None else float(toc_mgL) / TOC_REF_MGL


# ----------------------------------------------------------------------------- one chemistry condition
@dataclass(frozen=True)
class Chemistry:
    """One chemistry condition, for a hidden truth (simulate.build_scenario(chem=...)) or for a grid
    (simgp.simulator_grid_24h(cond=...)).

    disinfectant  : 'free_chlorine' (measured as free chlorine) or 'chloramine' (measured as total chlorine);
                    a chloramine truth or grid is refused (NotImplementedError) until task 12 builds its physics
    kinetics      : 'first' (the repo's path), 'first_si' (first order with the corrected unit recipe),
                    'order2' and 'clark' (EPANET order 2, truth only; built in later tasks)
    temp_C        : water temperature; None = the file's own properties and the calibrated 20 C rates
    toc_mgL       : plant TOC; None = not logged (bulk rate unscaled)
    threshold_mgL : compliance threshold; None = 0.2 free chlorine / 0.5 total chlorine (see THRESHOLD_NOTE)
    er_K          : E/R used by a MODEL condition (a grid); hidden truths draw their own E/R per seed
    wall_mode     : 'arrhenius' scales the wall rate by the same factor as the bulk rate (an ASSUMPTION);
                    'mass_transfer_only' leaves wall chemistry at 20 C (ablation M1b)
    wall_er_K     : task 10b, a MODEL condition's own wall E/R under 'arrhenius': the wall rate is scaled by
                    f(T; wall_er_K) while the bulk rate keeps f(T; er_K).  None (the default) means the wall follows
                    er_K, as before.  It is stored in canonical form, so a condition that is physically one of the
                    task-10 conditions IS that condition, with its cache tag and its cached grid: wall_er_K equal to
                    er_K is stored as None, and wall_er_K = 0 (f = 1 at every T) as wall_mode 'mass_transfer_only'.
    """
    disinfectant: str = FREE_CHLORINE
    kinetics: str = "first"
    temp_C: float | None = None
    toc_mgL: float | None = None
    threshold_mgL: float | None = None
    er_K: float = ER_DEFAULT_K
    wall_mode: str = "arrhenius"
    wall_er_K: float | None = None

    def __post_init__(self):
        # numbers are stored as Python floats, so equal conditions (10, 10.0, numpy 10) share one cache_tag
        for f in ("temp_C", "toc_mgL", "threshold_mgL", "er_K"):
            v = getattr(self, f)
            if v is not None:
                object.__setattr__(self, f, float(v))
        if not self.er_K > 0.0:
            raise ValueError(f"er_K must be positive, got {self.er_K}")
        if self.disinfectant not in DISINFECTANTS:
            raise ValueError(f"disinfectant must be one of {DISINFECTANTS}, got {self.disinfectant!r}")
        if self.kinetics not in KINETICS:
            raise ValueError(f"kinetics must be one of {KINETICS}, got {self.kinetics!r}")
        if self.wall_mode not in WALL_MODES:
            raise ValueError(f"wall_mode must be one of {WALL_MODES}, got {self.wall_mode!r}")
        if self.wall_er_K is not None:
            w = float(self.wall_er_K)
            if not (w >= 0.0 and math.isfinite(w)):
                raise ValueError(f"wall_er_K must be zero or positive, got {self.wall_er_K}")
            if self.wall_mode != "arrhenius":
                raise ValueError("wall_er_K is a wall E/R under wall_mode 'arrhenius'; 'mass_transfer_only' already "
                                 "fixes the wall chemistry at 20 C")
            if w == self.er_K:          # canonical forms (task 10b): the task-10 conditions keep their cache tags
                w = None
            elif w == 0.0:
                object.__setattr__(self, "wall_mode", "mass_transfer_only")
                w = None
            object.__setattr__(self, "wall_er_K", w)
        if self.temp_C is not None and not (0.0 <= self.temp_C <= 35.0):
            raise ValueError(f"temp_C {self.temp_C} is outside 0 to 35 C, beyond every rate law used here")
        if self.toc_mgL is not None and not (0.0 < self.toc_mgL <= 20.0):
            raise ValueError(f"toc_mgL {self.toc_mgL} is outside (0, 20] mg/L")
        if self.threshold_mgL is not None and not (0.0 < self.threshold_mgL < 5.0):
            raise ValueError(f"threshold_mgL {self.threshold_mgL} is outside (0, 5) mg/L")

    @property
    def species(self) -> str:
        return SPECIES[self.disinfectant]

    @property
    def threshold(self) -> float:
        return float(self.threshold_mgL) if self.threshold_mgL is not None else DEFAULT_THRESHOLD_MGL[self.disinfectant]

    @property
    def threshold_note(self) -> str:
        if self.threshold_mgL is None:
            return THRESHOLD_NOTE[self.disinfectant]
        return f"{self.threshold:g} mg/L {self.species}, set by the user (California requires a detectable residual)"

    @property
    def is_default(self) -> bool:
        """True when this condition reproduces today's free-chlorine first-order path exactly: the grid it
        describes is the committed one and needs no cache tag."""
        return (self.disinfectant == FREE_CHLORINE and self.kinetics == "first"
                and (self.temp_C is None or self.temp_C == TREF_C)
                and (self.toc_mgL is None or self.toc_mgL == TOC_REF_MGL))

    def kb_scale(self) -> float:
        """Bulk multiplier for a MODEL condition: f(T; E) (TOC / TOC_ref).  Exactly 1 at 20 C and TOC_ref."""
        f = 1.0 if self.temp_C is None else arrhenius(self.temp_C, self.er_K)
        return f * toc_ratio(self.toc_mgL)

    def kw_scale(self) -> float:
        """Wall multiplier for a MODEL condition: f(T; E) under 'arrhenius' (f(T; wall_er_K) when a wall E/R of its
        own is set, task 10b), 1 under 'mass_transfer_only'."""
        if self.temp_C is None or self.wall_mode == "mass_transfer_only":
            return 1.0
        return arrhenius(self.temp_C, self.er_K if self.wall_er_K is None else self.wall_er_K)

    def require_built(self) -> None:
        """Refuse a disinfectant whose decay physics is not built yet (chloramine until task 12)."""
        if self.disinfectant not in BUILT_DISINFECTANTS:
            raise NotImplementedError(f"{self.disinfectant} decay physics is not built yet (task 12); running it "
                                      "now would give free-chlorine physics under a chloramine label")

    def sim_kwargs(self) -> dict:
        """Keyword arguments for simulate.simulate_nominal_chlorine under this condition."""
        self.require_built()
        if self.kinetics != "first":
            raise ValueError("the calibrated grid is first order in the repo's units; "
                             f"kinetics {self.kinetics!r} is for hidden truths only")
        return {"kb_scale": self.kb_scale(), "kw_scale": self.kw_scale(), "temp_C": self.temp_C}

    def label(self) -> str:
        """Short human-readable part of a cache file name, e.g. 'free_T12.5_E8000'."""
        parts = ["free" if self.disinfectant == FREE_CHLORINE else "ca"]
        if self.kinetics != "first":
            parts.append(self.kinetics)
        if self.temp_C is not None:
            parts.append(f"T{self.temp_C:g}")
            parts.append(f"E{self.er_K:g}")
        if self.toc_mgL is not None:
            parts.append(f"TOC{self.toc_mgL:g}")
        if self.wall_mode != "arrhenius":
            parts.append("M1b")
        if self.wall_er_K is not None:
            parts.append(f"W{self.wall_er_K:g}")
        return "_".join(parts)

    def with_(self, **kw) -> "Chemistry":
        return replace(self, **kw)


def cache_tag(cond: Chemistry | None, grid: str = "full", n_chars: int = 10) -> str:
    """Short sha1 over everything that changes a cached grid: CHEM_CACHE_VERSION, the GRIDS contents of `grid`,
    FLOOR, HOURS, the quality step and run length, and the condition's kinetics, temperature, TOC, E/R, wall mode,
    wall E/R and disinfectant.  The threshold is left out on purpose: it does not change a simulated grid.  An unset
    wall E/R (None, the wall following the bulk E/R) is left out too, so every task-10 tag is unchanged."""
    from . import simgp, simulate      # lazy: simgp imports this module
    payload = {"version": CHEM_CACHE_VERSION, "grid": grid,
               "grid_axes": [list(map(float, ax)) for ax in simgp.GRIDS[grid]],
               "floor": simgp.FLOOR, "hours": list(simgp.HOURS),
               "quality_step_s": simulate.QUALITY_STEP_S, "duration_days": simulate.DURATION_DAYS,
               "cond": None if cond is None else {k: v for k, v in asdict(cond).items()
                                                  if k != "threshold_mgL" and not (k == "wall_er_K" and v is None)}}
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:n_chars]


# ----------------------------------------------------------------------------- EPANET files and units
EPANET_SUFFIXES = (".inp", ".rpt", ".bin", ".hyd", ".msx", ".msx-rpt", ".msx-bin", ".check.msx")


def remove_epanet_files(prefix: str) -> int:
    """Delete the files one EPANET run left behind for `prefix` (WNTR writes prefix.inp, .rpt, .bin).
    Returns how many were removed.  Missing files are fine."""
    n = 0
    for f in glob.glob(glob.escape(prefix) + ".*"):
        if f.endswith(EPANET_SUFFIXES):
            try:
                os.remove(f)
                n += 1
            except FileNotFoundError:
                pass
    return n


MGL_PER_KGM3 = 1000.0
NONLINEAR_TOLERANCE = 1e-5     # WNTR writes TOLERANCE unconverted (mg/L in the .inp); 1e-5 at a 1.2 mg/L scale is
                               # the same relative tolerance as the default 0.01 at the legacy 1200 mg/L scale


def use_mg_per_litre(wn, doses: dict[str, float], tolerance: float = NONLINEAR_TOLERANCE) -> float:
    """The verified unit recipe for EPANET kinetics that are not linear in concentration (truth only):
      * sources at dose/1000 kg/m3 (WNTR's SI unit), so the .inp says CONCEN <dose> mg/L, not <dose x 1000>;
      * every node.initial_quality set to 0 (Net2.inp starts every node at 1.0 mg/L, which becomes material
        once the dose is no longer x1000);
      * quality tolerance 1e-5.
    Bulk coefficients of order 0 or 2 and the limiting potential must then be given RAW in .inp units
    (mg/L/day, L/mg/day, mg/L): see set_bulk_kinetics.  Returns the factor that turns EPANET's quality
    results back into mg/L (1000)."""
    wn.options.quality.tolerance = tolerance
    for _, node in wn.nodes():
        node.initial_quality = 0.0
    for node, dose in doses.items():
        wn.add_source(f"src_{node}", node, "CONCEN", float(dose) / MGL_PER_KGM3)
    return MGL_PER_KGM3


def set_bulk_kinetics(wn, order: int, coeff: float, limiting_mgL: float | None = None) -> None:
    """Global bulk (and tank) reaction of the given order.  coeff is NEGATIVE for decay and in .inp units:
    1/day for order 1 (WNTR converts it, so it is divided by DAY here), mg/L/day for order 0 and L/mg/day for
    order 2 (WNTR writes these unconverted).  limiting_mgL is EPANET's LIMITING POTENTIAL in mg/L: with
    order 2 and coeff = -k2 it gives Clark's single-reactant form dC/dt = -k2 C (C - CL).  Only valid after
    use_mg_per_litre.  The .inp writer keeps 4 decimals, so |coeff| below 5e-5 would be written as 0."""
    if order not in (0, 1, 2):
        raise ValueError("order must be 0, 1 or 2")
    if order != 1 and 0 < abs(coeff) < 5e-5:
        raise ValueError(f"bulk coefficient {coeff} would be written as 0.0000 in the .inp")
    wn.options.reaction.bulk_order = order
    wn.options.reaction.tank_order = order
    wn.options.reaction.bulk_coeff = coeff / DAY if order == 1 else coeff
    if limiting_mgL is not None:
        wn.options.reaction.limiting_potential = float(limiting_mgL)


# ----------------------------------------------------------------------------- closed forms for the checks
def analytic_concentration(kind: str, c0: float, k: float, t_days, cl: float = 0.0):
    """Closed-form batch (plug-flow) solutions, t in days:
      order1        C0 exp(-k t)                         k in 1/day
      order2        C0 / (1 + k C0 t)                    k in L/mg/day
      order0        max(C0 - k t, 0)                     k in mg/L/day
      order1_limit  CL + (C0 - CL) exp(-k t)             first order toward a limiting concentration CL
      clark         CL / (1 + (CL/C0 - 1) exp(-k CL t))  dC/dt = -k C (C - CL), Clark 1998 (J Environ Eng
                                                          124(1):16); tends to order2 as CL -> 0"""
    import numpy as np
    t = np.asarray(t_days, dtype=float)
    if kind == "order1":
        return c0 * np.exp(-k * t)
    if kind == "order2":
        return c0 / (1.0 + k * c0 * t)
    if kind == "order0":
        return np.maximum(c0 - k * t, 0.0)
    if kind == "order1_limit":
        return cl + (c0 - cl) * np.exp(-k * t)
    if kind == "clark":
        if cl == 0.0:
            return c0 / (1.0 + k * c0 * t)
        return cl / (1.0 + (cl / c0 - 1.0) * np.exp(-k * cl * t))
    raise ValueError(kind)
