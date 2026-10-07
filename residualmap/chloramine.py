"""
chloramine.py: the chloramine mode (iteration 4, journal task 12), kept strictly separate from free chlorine.

A chloraminated system is a different measurement and a different chemistry: grab samples read TOTAL chlorine (mostly
monochloramine), which decays one to two orders of magnitude more slowly than free chlorine, through an inorganic
autodecomposition that depends on pH and the chlorine-to-ammonia ratio, plus reactions with organic matter and the pipe
walls, and, where the residual runs low in warm old water, nitrification.  This module holds:

  * EPA's unified chloramine model, the constants copied from the EPA web application's app.R (Wahman 2018, J AWWA
    110(11):E43, doi:10.1002/awwa.1146; inorganic chemistry of Jafvert & Valentine 1992, doi:10.1021/es00027a022, and
    Vikesland et al. 2001, doi:10.1016/S0043-1354(00)00406-1; two-site organic matter of Duirk et al. 2005, Water Res
    39:3418, doi:10.1016/j.watres.2005.06.003): 8 species (TOTNH, TOTCL, NH2CL, NHCL2, NCL3, I, DOC1, DOC2), 14
    inorganic and 2 organic reactions, constant pH.  Two implementations of the same equations: a Python batch port
    (scipy solve_ivp, the reference) and the EPANET-MSX model that runs on a network (build_msx_model); a saved check
    holds them within 0.05 mg/L of each other.  app.R converts with 71000 mg/mol Cl2 and 14000 mg/mol N; this module
    uses 70906 and 14007 (0.13% and 0.05% apart), as the plan says.
  * the kb prior table: the apparent first-order rate of total chlorine from the batch port at 20 C over pH 7 to 9 and
    Cl2:N 3 to 5 (mass), and k_hat(pH, Cl2:N) interpolated from it (model (c)'s prior centre);
  * the hidden chloramine truth on a network (truth_quality): the committed per-pipe, demand and dose draws of
    simulate.build_scenario, the EPA model in MSX with a first-order, mass-transfer-limited wall term on NH2CL, and an
    optional nitrification stress; and its first-order twin (EPANET, kb = the batch port's apparent rate at the seed's
    chemistry), which isolates what first-order kinetics cost;
  * the nitrification watch: a rule, not a fit, with literature thresholds that are not validated here.

What is assumed and what is from the literature is said where it is used.  No chloramine number is ever put on the same
axis, table or figure as a free-chlorine number: a free-chlorine MODEL may be scored against a chloramine TRUTH (the
conflation test, model (a)), but every number is total chlorine.

THE EXPERIMENT (journal task 12; bars pre-registered there):
    python -m residualmap.experiment Net3 --disinfectant=chloramine --workers 6     # also Net2; about 2 min and 30 s
    python -m residualmap.experiment --disinfectant=chloramine --calibrate           # the kw_ref and LIK_SD rules
    python -m residualmap.experiment Net3 --disinfectant=chloramine --resummarise    # summary from the committed CSVs
    python -m residualmap.experiment Net3 --disinfectant=chloramine --replot         # figures (recomputes seed 300)
Models (a) to (g) (MODELS) on fresh seeds (Net3 300 to 307 plus the low-pH stress seeds 308 to 311; Net2 300 to 307),
daytime samples n = 3, 8, 15 under a random and a straddle-on-daily-minimum rule, scored on the daily-minimum total
chlorine at unsampled junctions; the nitrification sweep (Fm 3, 10, 30) on Net3; the first-order twin on both.
Outputs in outputs/chloramine/: results_time_<net>.csv (sums per variant, seed, model, rule and n, 6 significant
digits), truths_<net>.csv, summary_chloramine_<net>.json (pooled rates, the acceptance key, the sweep), the conflation
figure chloramine_mode_vs_free_grid_<net>.png, nitrification_stress_Net3.png, kb_prior_table.csv and
calibration_chloramine.json.  Simulation only.
"""
from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd

MW_CL2 = 70906.0            # mg per mol Cl2 (app.R: 71000)
MW_N = 14007.0              # mg per mol N (app.R: 14000)
MW_C = 12000.0              # mg per mol C (app.R)
S1_FAST = 0.020             # Duirk 2005 / EPA app default: fraction of TOC in the site that reacts with NH2Cl (DOC1)
S2_SLOW = 0.65              # the fraction that reacts with HOCl (DOC2)
D_FAST_PER_TOC = MW_CL2 / MW_C * S1_FAST   # 0.118 mg Cl2 per mg C: the fast demand DOC1 exerts on NH2Cl
SPECIES = ("TOTNH", "TOTCL", "NH2CL", "NHCL2", "NCL3", "I", "DOC1", "DOC2")
TREF_C = 20.0
DAY = 86400.0

# ------------------------------------------------------------------- the truth's draws and settings (task 12 plan)
CA_DOSE_MGL = 2.0                     # nominal plant dose, mg/L as Cl2 (typical chloramine doses are 1.5 to 4.0)
PH_RANGE = (7.5, 8.5)                 # validity range of Roy, Sathasivan & Kastl 2020 (doi:10.1016/j.scitotenv.2020.140410)
PH_RANGE_STRESS = (7.0, 7.5)          # the low-pH stress subset (batch first-order fits 17 to 22% low at 7 days)
CL2N_RANGE = (4.0, 5.0)               # mass ratio; monochloramine dominates below 5.07:1 (Wahman 2018)
TOC_RANGE = (1.0, 3.0)                # mg C/L, ASSUMPTION inside the planning probe's 0 to 6
ALK_RANGE = (50.0, 150.0)             # mg/L as CaCO3, inside the probe's 50 to 200
TRUTH_TEMP_C = 20.0                   # a seasonal chloramine test is not run
WATCH_CRIT_MGL = 0.4                  # nitrification critical residual (range 0.2 to 0.65; Sathasivan, Fisher & Tam 2008)
WATCH_P = 0.5
WATCH_MIN_TEMP_C = 15.0               # nitrification observed from 8 to 26 C, optimum 25 to 30 C (US EPA 2002)
WATCH_AGE_QUANTILE = 0.75
FM_SWEEP = (0.0, 3.0, 10.0, 30.0)     # nitrification stress multiplier: a sweep only, no typical value is published
                                      # (Sathasivan, Fisher & Kastl 2005, doi:10.1021/es048300u)
FIT_WINDOW_D = (0.25, 3.0)            # apparent first-order fits of the batch port: from 6 h (after the fast organic
                                      # demand) to 3 days
PRIOR_TABLE_PH = tuple(np.round(np.arange(7.0, 9.0001, 0.25), 2))
PRIOR_TABLE_CL2N = (3.0, 3.5, 4.0, 4.5, 5.0)
PRIOR_TABLE_SETTINGS = {"temp_C": 20.0, "dose_mgL": CA_DOSE_MGL, "alk_mgL_caco3": 100.0, "toc_mgL": 0.0}
PRIOR_SD_LN = 1.0                     # ASSUMPTION: SD of ln kb20 around k_hat in model (c)


# ===================================================================== EPA unified model: constants and batch port
def equilibria(pH: float, alk_mgL_caco3: float, temp_C: float) -> dict:
    """Acid-base equilibria of app.R (temperature-adjusted pKa values), species fractions and carbonate species (mol/L)."""
    T = temp_C + 273.15
    KHOCl = 10 ** (-(1.18e-4 * T ** 2 - 7.86e-2 * T + 20.5))
    KNH4 = 10 ** (-(1.03e-4 * T ** 2 - 9.21e-2 * T + 27.6))
    KH2CO3 = 10 ** (-(1.48e-4 * T ** 2 - 9.39e-2 * T + 21.2))
    KHCO3 = 10 ** (-(1.19e-4 * T ** 2 - 7.99e-2 * T + 23.6))
    KW = 10 ** (-(1.5e-4 * T ** 2 - 1.23e-1 * T + 37.3))
    H = 10.0 ** -pH
    OH = KW / H
    a0CO = 1 / (1 + KH2CO3 / H + KH2CO3 * KHCO3 / H ** 2)
    a1CO = 1 / (1 + H / KH2CO3 + KHCO3 / H)
    a2CO = 1 / (1 + H / KHCO3 + H ** 2 / (KH2CO3 * KHCO3))
    TOTCO = (alk_mgL_caco3 / 50000 + H - OH) / (a1CO + 2 * a2CO)
    return {"H": H, "OH": OH, "a0Cl": 1 / (1 + KHOCl / H), "a1Cl": 1 / (1 + H / KHOCl), "a1NH": 1 / (1 + H / KNH4),
            "H2CO3": a0CO * TOTCO, "HCO3": a1CO * TOTCO, "CO3": a2CO * TOTCO}


def rate_constants(pH: float, alk_mgL_caco3: float, temp_C: float) -> dict:
    """app.R's rate constants (mol/L and seconds), temperature-adjusted where app.R adjusts them, plus the equilibria."""
    e = equilibria(pH, alk_mgL_caco3, temp_C)
    T = temp_C + 273.15
    k5 = (1.05e7 * math.exp(-2169 / T)) * e["H"] + (4.2e31 * math.exp(-22144 / T)) * e["HCO3"] \
        + (8.19e6 * math.exp(-4026 / T)) * e["H2CO3"]
    return {**e, "k1": 6.6e8 * math.exp(-1510 / T), "k2": 1.38e8 * math.exp(-8800 / T), "k3": 3.0e5 * math.exp(-2010 / T),
            "k4": 6.5e-7, "k5": k5, "k6": 6.0e4, "k7": 1.1e2, "k8": 2.8e4, "k9": 8.3e3, "k10": 1.5e-2,
            "k11p": 3.28e9 * e["OH"] + 6.0e6 * e["CO3"], "k11OCl": 9e4, "k12": 5.56e10, "k13": 1.39e9, "k14": 2.31e2,
            "kDOC1": 5.4, "kDOC2": 180.0}


def initial_state(dose_mgL: float, cl2n: float, toc_mgL: float, mode: str = "preformed") -> np.ndarray:
    """Initial concentrations (mol/L) in SPECIES order.  'preformed': chloramine formed at the plant, all of the dose as
    NH2Cl and the rest of the ammonia free (TOTNH = N added - N in NH2Cl), no free chlorine: what enters a network.
    'simadd': app.R's simultaneous addition of free chlorine and ammonia at the Cl2:N mass ratio (Wahman's examples)."""
    cl = float(dose_mgL) / MW_CL2
    n_total = cl * (MW_CL2 / (MW_N * float(cl2n)))       # mol N per L at a Cl2:N mass ratio
    doc1, doc2 = toc_mgL * S1_FAST / MW_C, toc_mgL * S2_SLOW / MW_C
    if mode == "preformed":
        if n_total < cl:
            raise ValueError(f"Cl2:N {cl2n} is past the monochloramine stoichiometry ({MW_CL2 / MW_N:.2f}:1): no free ammonia")
        return np.array([n_total - cl, 0.0, cl, 0.0, 0.0, 0.0, doc1, doc2])
    if mode == "simadd":
        return np.array([n_total, cl, 0.0, 0.0, 0.0, 0.0, doc1, doc2])
    raise ValueError(mode)


def _rhs(k: dict):
    H, OH, a0Cl, a1Cl, a1NH = k["H"], k["OH"], k["a0Cl"], k["a1Cl"], k["a1NH"]

    def f(t, y):
        TOTNH, TOTCl, NH2Cl, NHCl2, NCl3, I, DOC1, DOC2 = y
        HOCl, OCl, NH3 = a0Cl * TOTCl, a1Cl * TOTCl, a1NH * TOTNH
        r1 = k["k1"] * HOCl * NH3
        r2 = k["k2"] * NH2Cl
        r3 = k["k3"] * HOCl * NH2Cl
        r4 = k["k4"] * NHCl2
        r5 = k["k5"] * NH2Cl ** 2
        r6 = k["k6"] * NHCl2 * NH3 * H
        r7 = k["k7"] * NHCl2 * OH
        r8 = k["k8"] * I * NHCl2
        r9 = k["k9"] * I * NH2Cl
        r10 = k["k10"] * NH2Cl * NHCl2
        r11 = (k["k11p"] + k["k11OCl"] * OCl) * HOCl * NHCl2
        r12 = k["k12"] * NHCl2 * NCl3 * OH
        r13 = k["k13"] * NH2Cl * NCl3 * OH
        r14 = k["k14"] * NHCl2 * OCl
        rD1 = k["kDOC1"] * NH2Cl * DOC1
        rD2 = k["kDOC2"] * HOCl * DOC2
        return [-r1 + r2 + r5 - r6 + rD1,
                -r1 + r2 - r3 + r4 + r8 - r11 + 2 * r12 + r13 - 2 * r14 - rD2,
                r1 - r2 - r3 + r4 - 2 * r5 + 2 * r6 - r9 - r10 - r13 - rD1,
                r3 - r4 + r5 - r6 - r7 - r8 - r10 - r11 - r12 - r14,
                r11 - r12 - r13,
                r7 - r8 - r9,
                -rD1, -rD2]
    return f


def total_chlorine_mgL(Y) -> np.ndarray:
    """Total chlorine as Cl2, mg/L, from species rows (SPECIES order) or a dict/DataFrame of species."""
    g = (lambda s: Y[s]) if not isinstance(Y, np.ndarray) else (lambda s: Y[SPECIES.index(s)])
    return MW_CL2 * (g("TOTCL") + g("NH2CL") + 2 * g("NHCL2") + 3 * g("NCL3"))


def batch(pH: float, alk_mgL_caco3: float, temp_C: float, dose_mgL: float, cl2n: float, toc_mgL: float,
          days: float = 10.0, t_eval_days=None, mode: str = "preformed") -> pd.DataFrame:
    """The Python batch port (ideal plug flow, no wall, no nitrification), as the EPA application solves it.  Returns
    t_days, total_mgL (total chlorine as Cl2), mono_mgL (NH2Cl as Cl2) and free_nh3_mgN.  LSODA, rtol 1e-8, atol 1e-14."""
    from scipy.integrate import solve_ivp
    k = rate_constants(pH, alk_mgL_caco3, temp_C)
    y0 = initial_state(dose_mgL, cl2n, toc_mgL, mode)
    if t_eval_days is None:
        t_eval_days = np.concatenate([np.linspace(0, 1 / 24, 61), np.linspace(2 / 24, days, 480)])
    t = np.asarray(t_eval_days, dtype=float) * DAY
    sol = solve_ivp(_rhs(k), (0.0, float(t.max())), y0, method="LSODA", t_eval=t, rtol=1e-8, atol=1e-14)
    if not sol.success:
        raise RuntimeError(f"batch port did not integrate: {sol.message}")
    Y = sol.y
    return pd.DataFrame({"t_days": sol.t / DAY, "total_mgL": total_chlorine_mgL(Y), "mono_mgL": Y[2] * MW_CL2,
                         "free_nh3_mgN": Y[0] * MW_N})


def apparent_rate(df: pd.DataFrame, window_days=FIT_WINDOW_D) -> tuple[float, float]:
    """Least-squares fit of ln total chlorine = a - k t over the window: (k per day, C(0) = exp(a) in mg/L)."""
    m = (df.t_days >= window_days[0]) & (df.t_days <= window_days[1]) & (df.total_mgL > 1e-3)
    A = np.column_stack([np.ones(int(m.sum())), -df.t_days[m].values])
    a, k = np.linalg.lstsq(A, np.log(df.total_mgL[m].values), rcond=None)[0]
    return float(k), float(math.exp(a))


def kb_prior_table(pHs=PRIOR_TABLE_PH, ratios=PRIOR_TABLE_CL2N, settings=PRIOR_TABLE_SETTINGS) -> pd.DataFrame:
    """Apparent first-order rate of total chlorine (per day, 20 C) from the batch port over pH x Cl2:N, at the stated
    dose, alkalinity and TOC 0 (inorganic decay only: the fast organic demand is the dose axis's job), fitted over
    FIT_WINDOW_D.  Written to outputs/chloramine/kb_prior_table.csv by the experiment."""
    rows = []
    for ph in pHs:
        for r in ratios:
            df = batch(ph, settings["alk_mgL_caco3"], settings["temp_C"], settings["dose_mgL"], r, settings["toc_mgL"], days=4.0)
            k, c0 = apparent_rate(df)
            rows.append({"pH": float(ph), "cl2n": float(r), "k_app_per_day": k, "c0_fit_mgL": c0,
                         "total_1d_mgL": float(np.interp(1.0, df.t_days, df.total_mgL)),
                         "total_3d_mgL": float(np.interp(3.0, df.t_days, df.total_mgL))})
    return pd.DataFrame(rows)


def k_hat(pH: float, cl2n: float, table: pd.DataFrame) -> float:
    """Bilinear interpolation of ln k_app over the prior table's (pH, Cl2:N) grid, clamped to its edges."""
    P = table.pivot(index="pH", columns="cl2n", values="k_app_per_day").sort_index().sort_index(axis=1)
    ph_v, r_v = P.index.values.astype(float), P.columns.values.astype(float)
    L = np.log(P.values)
    x = float(np.clip(pH, ph_v[0], ph_v[-1])); y = float(np.clip(cl2n, r_v[0], r_v[-1]))
    i = int(np.clip(np.searchsorted(ph_v, x) - 1, 0, len(ph_v) - 2)); j = int(np.clip(np.searchsorted(r_v, y) - 1, 0, len(r_v) - 2))
    tx = (x - ph_v[i]) / (ph_v[i + 1] - ph_v[i]); ty = (y - r_v[j]) / (r_v[j + 1] - r_v[j])
    v = (1 - tx) * (1 - ty) * L[i, j] + tx * (1 - ty) * L[i + 1, j] + (1 - tx) * ty * L[i, j + 1] + tx * ty * L[i + 1, j + 1]
    return float(math.exp(v))


# ===================================================================== the same model in EPANET-MSX
TERMS = {
    "r1": "k1*a0Cl*TOTCL*a1NH*TOTNH", "r2": "k2*NH2CL", "r3": "k3*a0Cl*TOTCL*NH2CL", "r4": "k4*NHCL2",
    "r5": "k5*NH2CL*NH2CL", "r6": "k6*NHCL2*a1NH*TOTNH*Hc", "r7": "k7*NHCL2*OHc", "r8": "k8*I*NHCL2",
    "r9": "k9*I*NH2CL", "r10": "k10*NH2CL*NHCL2", "r11": "(k11p + k11OCl*a1Cl*TOTCL)*a0Cl*TOTCL*NHCL2",
    "r12": "k12*NHCL2*NCL3*OHc", "r13": "k13*NH2CL*NCL3*OHc", "r14": "k14*NHCL2*a1Cl*TOTCL",
    "rD1": "kD1*NH2CL*DOC1", "rD2": "kD2*a0Cl*TOTCL*DOC2",
}
RATES = {
    "TOTNH": "-r1 + r2 + r5 - r6 + rD1",
    "TOTCL": "-r1 + r2 - r3 + r4 + r8 - r11 + 2*r12 + r13 - 2*r14 - rD2",
    "NH2CL": "r1 - r2 - r3 + r4 - 2*r5 + 2*r6 - r9 - r10 - r13 - rD1",
    "NHCL2": "r3 - r4 + r5 - r6 - r7 - r8 - r10 - r11 - r12 - r14",
    "NCL3": "r11 - r12 - r13", "I": "r7 - r8 - r9", "DOC1": "-rD1", "DOC2": "-rD2",
}
MSX_CONSTANTS = ("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11p", "k11OCl", "k12", "k13", "k14")
MSX_OPTIONS = {"rate_units": "SEC", "area_units": "M2", "solver": "ROS2", "coupling": "NONE", "timestep": 300,
               "atol": 1e-12, "rtol": 1e-4}
SPECIES_TOL = (1e-14, 1e-4)


def build_msx_model(wn, k: dict, kw_pipe_m_day: dict | None = None, fm: float = 0.0, fm_pipes=(), kn_per_day: float = 0.0):
    """The EPA model as an MSX reaction model for `wn` (sources are added by the caller).
    kw_pipe_m_day: per-pipe first-order NH2Cl wall coefficient (m/day); the wall rate is EPANET's mass-transfer-limited
    (4/D) kw Kf / (kw + Kf) NH2Cl, with chlorine's molecular diffusivity (ASSUMED for monochloramine).  None: no wall.
    fm > 0 adds the nitrification stress: an extra NH2Cl loss fm k_n NH2Cl while total chlorine is below 0.4 mg/L, on
    the pipes in fm_pipes (a per-pipe PARAMETER) and in every tank (a CONSTANT; WNTR cannot write per-tank values);
    k_n = kn_per_day (the batch apparent rate at the water's pH and Cl2:N)."""
    from wntr.msx import MsxModel
    from . import msx as M
    m = MsxModel()
    for key, v in MSX_OPTIONS.items():
        setattr(m.options, key, v)
    for sp in SPECIES:
        m.add_species(sp, "bulk", units="MOL", atol=SPECIES_TOL[0], rtol=SPECIES_TOL[1])
    for name in MSX_CONSTANTS:
        m.add_constant(name, float(k[name]))
    for name, key in (("kD1", "kDOC1"), ("kD2", "kDOC2"), ("Hc", "H"), ("OHc", "OH"), ("a0Cl", "a0Cl"),
                      ("a1Cl", "a1Cl"), ("a1NH", "a1NH")):
        m.add_constant(name, float(k[key]))
    for name, expr in TERMS.items():
        m.add_term(name, expr)
    pipe_nh2cl, tank_nh2cl = RATES["NH2CL"], RATES["NH2CL"]
    if kw_pipe_m_day is not None:
        M.add_mass_transfer_terms(m, wn)
        m.add_parameter("kwp", 0.0)
        M.set_pipe_parameter(m, "kwp", {p: M.wall_rate_unit(wn, v) for p, v in kw_pipe_m_day.items()})
        pipe_nh2cl += " - (4/D)*kwp*Kf/(kwp+Kf)*NH2CL"
    if fm > 0:
        m.add_constant("MWc", MW_CL2)
        m.add_constant("TCcrit", WATCH_CRIT_MGL)
        m.add_constant("kn", float(kn_per_day) / DAY)
        m.add_constant("FmT", float(fm))
        m.add_term("TCmg", "MWc*(TOTCL + NH2CL + 2*NHCL2 + 3*NCL3)")
        m.add_parameter("fmp", 0.0)
        M.set_pipe_parameter(m, "fmp", {p: float(fm) for p in fm_pipes})
        pipe_nh2cl += " - fmp*kn*NH2CL*step(TCcrit - TCmg)"
        tank_nh2cl += " - FmT*kn*NH2CL*step(TCcrit - TCmg)"
    for sp, expr in RATES.items():
        m.add_reaction(sp, "pipe", "rate", pipe_nh2cl if sp == "NH2CL" else expr)
        m.add_reaction(sp, "tank", "rate", tank_nh2cl if sp == "NH2CL" else expr)
    return m


def first_order_msx_model(wn, kb_per_day: float, kw_pipe_m_day: dict):
    """One species, EPANET's own first-order bulk and mass-transfer-limited wall reactions written as MSX expressions
    (a saved check runs it against EPANET's CHEMICAL quality on Net3: it validates the wall term and the runner)."""
    from wntr.msx import MsxModel
    from . import msx as M
    m = MsxModel()
    for key, v in MSX_OPTIONS.items():
        setattr(m.options, key, v)
    m.add_species("CL", "bulk", units="MG", atol=1e-8, rtol=1e-6)
    m.add_constant("kbulk", float(kb_per_day) / DAY)
    M.add_mass_transfer_terms(m, wn)
    m.add_parameter("kwp", 0.0)
    M.set_pipe_parameter(m, "kwp", {p: M.wall_rate_unit(wn, v) for p, v in kw_pipe_m_day.items()})
    m.add_reaction("CL", "pipe", "rate", "-kbulk*CL - (4/D)*kwp*Kf/(kwp+Kf)*CL")
    m.add_reaction("CL", "tank", "rate", "-kbulk*CL")
    return m


def source_species(dose_mgL: float, cl2n: float, toc_mgL: float) -> dict:
    """Preformed chloramine at a source: {species: mol/L} (NH2Cl, free ammonia, the two organic sites)."""
    y = initial_state(dose_mgL, cl2n, toc_mgL, "preformed")
    return {sp: float(v) for sp, v in zip(SPECIES, y) if sp in ("NH2CL", "TOTNH", "DOC1", "DOC2") and v > 0}


# ===================================================================== the hidden chloramine truth on a network
@dataclass(frozen=True)
class ChloramineTruth:
    """Settings of a chloramine truth beyond the committed draws (simulate.build_scenario(chloramine_truth=...)).
    ph_range : the range the seed's u_pH draw is mapped onto (PH_RANGE; PH_RANGE_STRESS for the low-pH subset)
    fm       : nitrification stress multiplier (0 = none; FM_SWEEP)
    compiler : MSX COMPILER, 'GC' with a loud fallback to 'NONE' (msx.run)"""
    ph_range: tuple = PH_RANGE
    fm: float = 0.0
    compiler: str = "GC"

    def __post_init__(self):
        lo, hi = map(float, self.ph_range)
        if not (6.5 <= lo < hi <= 9.5):
            raise ValueError(f"ph_range {self.ph_range} is outside 6.5 to 9.5")
        object.__setattr__(self, "ph_range", (lo, hi))
        object.__setattr__(self, "fm", float(self.fm))
        if self.fm < 0:
            raise ValueError("fm must be zero or positive")
        if self.compiler not in ("GC", "NONE"):
            raise ValueError("compiler must be 'GC' or 'NONE'")


def truth_chemistry(draws: dict, truth: ChloramineTruth | None = None) -> dict:
    """The seed's water chemistry from its raw uniforms (simulate.hidden_chem_draws): pH, Cl2:N (mass), TOC (mg C/L)
    and alkalinity (mg/L as CaCO3), each lo + (hi - lo) u on its range; temperature 20 C."""
    t = truth or ChloramineTruth()
    def m(u, r):
        return float(r[0] + (r[1] - r[0]) * u)
    return {"pH": m(draws["u_pH"], t.ph_range), "cl2n": m(draws["u_cl2n"], CL2N_RANGE),
            "toc_mgL": m(draws["u_toc"], TOC_RANGE), "alk_mgL_caco3": m(draws["u_alk"], ALK_RANGE),
            "temp_C": TRUTH_TEMP_C}


def aged_pipes(name: str, wn, quantile: float = WATCH_AGE_QUANTILE) -> tuple[list[str], float]:
    """Pipes whose downstream node's nominal daily-mean water age is at least the given quantile of the junctions'
    nominal daily-mean ages (the file's end node when the mean flow is zero), and that age threshold in hours."""
    from .simulate import _nominal_age_and_flow
    age, flow = _nominal_age_and_flow(name)
    thr = float(np.quantile([age[j] for j in wn.junction_name_list], quantile))
    out = []
    for pn, p in wn.pipes():
        down = p.start_node_name if flow.get(pn, 0.0) < 0 else p.end_node_name
        if age.get(down, 0.0) >= thr:
            out.append(pn)
    return out, thr


def twin_rate(chem: dict, dose_mgL: float) -> tuple[float, float]:
    """The first-order twin's bulk rate (per day) and effective source factor: the batch port at the seed's chemistry
    and the nominal dose, fitted over FIT_WINDOW_D; the factor is C(0) of the fit over the dose (the fast organic demand
    removed)."""
    df = batch(chem["pH"], chem["alk_mgL_caco3"], chem["temp_C"], dose_mgL, chem["cl2n"], chem["toc_mgL"], days=4.0)
    k, c0 = apparent_rate(df)
    return k, c0 / float(dose_mgL)


def chloramine_truth(wn, seed: int, rng, rng_m, chem, source_dose: float, kw_ref_m_day: float, name: str | None,
                     truth: ChloramineTruth | None = None):
    """The hidden chloramine truth (called by simulate._chem_truth; build_scenario(chem=Chemistry(disinfectant=
    'chloramine', kinetics=...), source_dose=2.0, kw_m_per_day=<kw_ref>)).  The committed draws are consumed exactly
    as in build_scenario's default branch: the monthly bulk factor u (consumed so every later draw stays aligned; EPA's
    mechanistic bulk chemistry has no such factor, so it is recorded and not used), per pipe the wall factor
    exp(N(0, 0.4)) and the roughness error exp(N(0, 0.10)), the global and per-node demand, and one dose U(0.9, 1.1)
    per source.  The new draws are simulate.hidden_chem_draws('chloramine', seed) (default_rng(40_000 + seed)).
      kinetics 'epa_msx': EPA's unified model in EPANET-MSX (build_msx_model), preformed chloramine at every source
        (dose source_dose x U, the seed's Cl2:N, TOC split into the two organic sites), zero initial quality, constant
        pH, 20 C; wall: NH2Cl first order, kw_p = kw_ref 2^(-(C_p - 130)/30) exp(N(0, 0.4)) (the committed form and
        draw), mass-transfer limited with chlorine's diffusivity (ASSUMED); optional nitrification stress (fm > 0) on
        the pipes downstream of the oldest quarter of the network's water and in every tank.
      kinetics 'first': the first-order twin: EPANET first order with kb = the batch port's apparent rate at the seed's
        chemistry (twin_rate), the same per-pipe wall coefficients as EPANET wall coefficients, and every source at its
        dose times the batch fit's C(0) / dose.  It isolates what first-order kinetics cost.
    Returns (node quality, total chlorine in mg/L as Cl2, time x node; info)."""
    from . import msx as M
    from .simulate import _run_quality, roughness_factor, source_nodes
    t = truth or ChloramineTruth()
    from .simulate import hidden_chem_draws
    hd = hidden_chem_draws("chloramine", seed)
    ch = truth_chemistry(hd, t)
    u = rng_m.uniform(0.8, 1.2)                     # committed monthly bulk factor: consumed, not used
    kw_pipe = {}
    for pn, pipe in wn.pipes():
        kw_pipe[pn] = kw_ref_m_day * roughness_factor(pipe.roughness, 1.0) * np.exp(rng.normal(0, 0.4))
        pipe.roughness = pipe.roughness * np.exp(rng.normal(0, 0.10))
    global_mult = rng_m.uniform(0.85, 1.15)
    for _, j in wn.junctions():
        for ts in j.demand_timeseries_list:
            ts.base_value = ts.base_value * global_mult * np.exp(rng_m.normal(0.0, 0.15))
    doses = {s: source_dose * rng_m.uniform(0.9, 1.1) for s in source_nodes(wn)}
    info = {"disinfectant": "chloramine", "species": "total chlorine", "kinetics": chem.kinetics, "rng_seed": hd["rng_seed"],
            **{k: v for k, v in ch.items()}, "ph_range": list(t.ph_range), "kw_ref_m_per_day": float(kw_ref_m_day),
            "bulk_month_factor_unused": float(u), "source_doses_mgL": {k: float(v) for k, v in doses.items()},
            "k_hat_per_day": None, "fm": t.fm, "quality_scale": 1.0}
    tab = prior_table_cached()
    info["k_hat_per_day"] = k_hat(ch["pH"], ch["cl2n"], tab)
    if chem.kinetics == "first":
        kb, factor = twin_rate(ch, source_dose)
        wn.options.quality.parameter = "CHEMICAL"
        wn.options.reaction.bulk_coeff = -kb / DAY
        wn.options.reaction.wall_coeff = -kw_ref_m_day / DAY
        for pn, pipe in wn.pipes():
            pipe.wall_coeff = -kw_pipe[pn] / DAY
        for s, d in doses.items():
            wn.add_source(f"src_{s}", s, "CONCEN", d * factor)
        info.update(kb_twin_per_day=kb, twin_source_factor=factor)
        return _run_quality(wn), info
    k = rate_constants(ch["pH"], ch["alk_mgL_caco3"], ch["temp_C"])
    fm_pipes, age_thr = ([], None)
    if t.fm > 0:
        if name is None:
            raise ValueError("the nitrification stress needs the network name (nominal water age)")
        fm_pipes, age_thr = aged_pipes(name, wn)
    m = build_msx_model(wn, k, kw_pipe, fm=t.fm, fm_pipes=fm_pipes, kn_per_day=info["k_hat_per_day"])
    M.add_sources(m, wn, {s: source_species(d, ch["cl2n"], ch["toc_mgL"]) for s, d in doses.items()})
    res, run_info = M.run(wn, m, compiler=t.compiler)
    q = total_chlorine_mgL({sp: res.node[sp] for sp in ("TOTCL", "NH2CL", "NHCL2", "NCL3")})
    last = q.loc[q.index >= 6 * DAY]
    info["tank_last_day_mean_mgL"] = {tn: float(last[tn].mean()) for tn in wn.tank_name_list}
    info.update(msx_compiler=run_info["compiler"], msx_fallback_reason=run_info["fallback_reason"],
                n_fm_pipes=len(fm_pipes), fm_age_threshold_h=age_thr,
                _volatile={"msx_seconds": run_info["seconds"]})
    return q, info


_PRIOR_TABLE = None


def prior_table_cached() -> pd.DataFrame:
    """kb_prior_table() computed once per process (0.2 s)."""
    global _PRIOR_TABLE
    if _PRIOR_TABLE is None:
        _PRIOR_TABLE = kb_prior_table()
    return _PRIOR_TABLE


# ===================================================================== truths in a pool, and the two settings rules
KW_LADDER = (0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.90)
KW_SHARE_BAND, KW_SHARE_TARGET = (0.10, 0.30), 0.20
TUNE_SEED = 0                       # the kw_ref rule's seed (never scored)
CALIB_SEEDS = (100, 101, 102, 103)  # LIK_SD_CA's seeds (never scored)
SINGLE_THREAD_ENV = {k: "1" for k in ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}


def _pool_init(workdir: str) -> None:
    """Worker: one thread per library, a private working directory (nominal_scenario writes EPANET files to the cwd)."""
    import warnings
    warnings.filterwarnings("ignore")
    os.environ.update(SINGLE_THREAD_ENV)
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)


def truth_scenario(net: str, seed: int, kw_ref: float, kinetics: str = "epa_msx", truth: ChloramineTruth | None = None,
                   dose: float = CA_DOSE_MGL, sample_hour: int = 14, month_seed: int | None = None):
    """build_scenario for a chloramine truth (the scenario's nominal side is the operator's file at the chloramine
    dose)."""
    from .chemistry import CHLORAMINE, Chemistry
    from .simulate import build_scenario
    return build_scenario(net, seed, sample_hour=sample_hour, source_dose=dose, kw_m_per_day=kw_ref, month_seed=month_seed,
                          chem=Chemistry(disinfectant=CHLORAMINE, kinetics=kinetics), chloramine_truth=truth)


def _truth_job(args):
    net, seed, kw_ref, kinetics, truth = args
    sc = truth_scenario(net, seed, kw_ref, kinetics, truth)
    info = dict(sc.chem); vol = info.pop("_volatile", {})
    return {"net": net, "seed": seed, "kw_ref": kw_ref, "kinetics": kinetics, "fm": truth.fm if truth else 0.0,
            "truth_by_hour": sc.truth_by_hour, "chem": info, "seconds": vol.get("msx_seconds")}


def run_truths(jobs: list, workers: int = 4, workdir: str | None = None) -> list:
    """Chloramine truths in a process pool (jobs: (net, seed, kw_ref, kinetics, ChloramineTruth or None)); returns the
    _truth_job dicts in job order.  Each worker gets a private cwd under a temporary directory, deleted afterwards."""
    import tempfile
    from concurrent.futures import ProcessPoolExecutor
    with tempfile.TemporaryDirectory(prefix="rm_ca_pool_", dir=workdir) as tmp:
        if workers <= 1:
            old = os.getcwd()
            try:
                _pool_init(tmp)
                return [_truth_job(j) for j in jobs]
            finally:
                os.chdir(old)
        with _single_thread_children(), ProcessPoolExecutor(max_workers=workers, initializer=_pool_init,
                                                            initargs=(tmp,)) as ex:
            return list(ex.map(_truth_job, jobs))


def _demo_truth_job(args):
    net, seed, kw_ref, dose = args
    sc = truth_scenario(net, seed, kw_ref, "epa_msx", None, dose=dose)
    info = dict(sc.chem)
    info.pop("_volatile", None)
    return {"truth_by_hour": sc.truth_by_hour, "truth_daily_min": sc.truth_daily_min, "chem": info}


def truth_in_subprocess(net: str, seed: int, kw_ref: float, dose: float = CA_DOSE_MGL):
    """One chloramine truth (EPA's chemistry in MSX) computed in a separate process with its own working directory,
    for the app's demo.  MSX needs a change of working directory, which is process-wide; in a separate process the
    app's own threads (one per browser session) never see it.  `net` is a bundled name or a file path.  Returns a
    namespace with truth_by_hour, truth_daily_min and chem, the same values truth_scenario gives."""
    import tempfile
    import types
    from concurrent.futures import ProcessPoolExecutor
    if os.path.exists(net):
        net = os.path.abspath(net)              # the worker runs in its own directory
    with tempfile.TemporaryDirectory(prefix="rm_ca_demo_") as tmp:
        with ProcessPoolExecutor(max_workers=1, initializer=_pool_init, initargs=(tmp,)) as ex:
            out = ex.submit(_demo_truth_job, (net, int(seed), float(kw_ref), float(dose))).result()
    return types.SimpleNamespace(**out)


def choose_kw_ref(shares: dict) -> float | None:
    """The plan's rule made exact (written down before the tuning run): among ladder values whose seed-0 share of
    junctions with a daily minimum below 0.5 mg/L is in [0.10, 0.30], the one closest to 0.20 (ties: the smaller)."""
    ok = [(abs(s - KW_SHARE_TARGET), kw) for kw, s in shares.items() if KW_SHARE_BAND[0] <= s <= KW_SHARE_BAND[1]]
    return min(ok)[1] if ok else None


def tune_kw_ref(nets=("Net3", "Net2"), workers: int = 4) -> dict:
    """The kw_ref rule on TUNE_SEED for each network.  Returns {net: {'shares': {kw: share}, 'kw_ref': chosen}}."""
    jobs = [(net, TUNE_SEED, kw, "epa_msx", None) for net in nets for kw in KW_LADDER]
    res = run_truths(jobs, workers)
    out = {}
    for r in res:
        dm = r["truth_by_hour"].min()
        out.setdefault(r["net"], {"shares": {}, "n_junctions": int(len(dm))})["shares"][r["kw_ref"]] = float((dm < 0.5).mean())
    for net in out:
        out[net]["kw_ref"] = choose_kw_ref(out[net]["shares"])
    return out


def ca_condition():
    """The chloramine grid's chemistry condition (first order, 20 C): its own cache tag, its own grid."""
    from .chemistry import CHLORAMINE, Chemistry
    return Chemistry(disinfectant=CHLORAMINE)


def nominal_ca(net: str, dose: float = CA_DOSE_MGL, sample_hour: int = 14):
    from .simulate import nominal_scenario
    return nominal_scenario(net, sample_hour, dose)


def ca_grid(net: str, cache_dir: str = "outputs/cache", dose: float = CA_DOSE_MGL):
    """The chloramine grid of a network at the chloramine dose: (params, Z [members x 24 x J]), cached in cache_dir."""
    from .simgp import simulator_grid_24h
    return simulator_grid_24h(nominal_ca(net, dose), cache_dir, "chloramine", cond=ca_condition())


def best_member(Z: np.ndarray, truth_by_hour: pd.DataFrame, doses, hours=None) -> dict:
    """The single (member, dose) of a grid with the smallest RMS ln mismatch against the truth over `hours` (the
    daytime window by default) and every junction: its index, dose, RMS ln mismatch, and its ln C for all 24 hours."""
    from .simgp import DAY_HOURS, FLOOR
    hours = list(DAY_HOURS if hours is None else hours)
    zt = np.log(np.clip(truth_by_hour.loc[hours].values, FLOOR, None))
    E = Z[:, hours, :].astype(float) - zt[None]
    m1, m2 = E.mean(axis=(1, 2)), (E ** 2).mean(axis=(1, 2))
    offs = np.log(np.asarray(doses, dtype=float))
    ms = m2[:, None] + 2 * offs[None, :] * m1[:, None] + offs[None, :] ** 2          # members x doses
    k, d = np.unravel_index(int(np.argmin(ms)), ms.shape)
    return {"member": int(k), "dose": float(doses[d]), "rms_ln": float(math.sqrt(ms[k, d])), "z": Z[k] + offs[d]}


def calibrate_lik_sd(kw_ref: dict, nets=("Net3", "Net2"), seeds=CALIB_SEEDS, cache_dir: str = "outputs/cache",
                     workers: int = 4) -> dict:
    """The LIK_SD rule for chloramine: on the calibration seeds (never scored), the daytime RMS ln mismatch of the best
    single (member, dose) of the chloramine grid against the main chloramine truth; the value is set from Net3 (as
    0.35 was for free chlorine): the mean over the seeds, rounded up to a multiple of 0.05.  Net2 is reported."""
    from .simgp import DOSES_CA, GRIDS
    jobs = [(net, s, kw_ref[net], "epa_msx", None) for net in nets for s in seeds]
    res = run_truths(jobs, workers)
    out = {}
    for net in nets:
        params, Z = ca_grid(net, cache_dir)
        rows = []
        for r in (x for x in res if x["net"] == net):
            b = best_member(Z, r["truth_by_hour"], DOSES_CA)
            p = params[b["member"]]
            rows.append({"seed": r["seed"], "rms_ln": b["rms_ln"], "dose": b["dose"],
                         **dict(zip(("kb", "kw", "gamma", "demand", "rough"), map(float, p))),
                         "pH": r["chem"]["pH"], "cl2n": r["chem"]["cl2n"], "toc_mgL": r["chem"]["toc_mgL"],
                         "k_hat_per_day": r["chem"]["k_hat_per_day"]})
        mean = float(np.mean([x["rms_ln"] for x in rows]))
        out[net] = {"seeds": rows, "mean_rms_ln": mean, "rounded_up_0.05": math.ceil(mean / 0.05 - 1e-9) * 0.05,
                    "grid_members": len(params), "grid_axes": [list(map(float, a)) for a in GRIDS["chloramine"]]}
    out["LIK_SD_CA"] = out[nets[0]]["rounded_up_0.05"]
    return out


# ===================================================================== the experiment (journal task 12)
# python -m residualmap.experiment <net> --disinfectant=chloramine [--workers 6]; bars pre-registered in the journal
OUT_DIR = os.path.join("outputs", "chloramine")
THRESHOLD = 0.5                       # mg/L total chlorine: a common utility operating target, not a California rule
THRESHOLD_NOTE = ("0.5 mg/L total chlorine: a common utility operating target, not a California rule (California "
                  "requires a detectable residual)")
FIRST_SEED = 300                      # fresh: no chloramine truth was scored before; 0 is the kw_ref rule's seed,
N_MAIN = 8                            # 100 to 103 the LIK_SD rule's, 900 the smoke run's
STRESS_SEEDS = {"Net3": tuple(range(308, 312))}
FM_NETS = ("Net3",)
SMOKE_SEED = 900
UNSCORED = frozenset({TUNE_SEED, *CALIB_SEEDS, SMOKE_SEED})
KW_REF = {"Net3": 0.20, "Net2": 0.05}  # by the kw_ref rule on seed 0 (outputs/chloramine/calibration_chloramine.json)
NS = (3, 8, 15)
N0, NMAX = 3, 15
RULES = ("random", "straddle_min")
NOISE_SD = 0.03
AFFECTED_DROP_MGL = 0.05              # a junction is nitrification-affected when its true daily minimum is at least this
                                      # much below its daily minimum without the stress (same seed and draws)
MODELS = {"a": "free-chlorine grid and settings (the conflation test)", "b": "chloramine grid, uniform prior",
          "c": "chloramine grid, pH and Cl2:N prior (the chloramine mode's model)", "d": "mean of samples",
          "e": "nearest sampled junction", "f": "best single grid member and dose (oracle)",
          "g": "model (b) on the first-order twin truth"}
GRID_MODELS = ("a", "b", "c", "g")
SUM_COLS = ("n_uns", "sse", "sae", "sse_ln", "se", "in50", "in80", "in90", "in95", "n_true_viol", "tp", "fp", "fn")
CSV_FLOAT = "%.6g"


def check_seeds(seeds) -> None:
    bad = sorted(set(seeds) & UNSCORED)
    if bad:
        raise ValueError(f"seeds {bad} were used to set kw_ref, LIK_SD_CA or the smoke run; they are never scored")


def prior_log_vector(params, pH: float, cl2n: float, table: pd.DataFrame | None = None) -> np.ndarray:
    """Model (c)'s log prior over grid members: ln kb20 ~ N(ln k_hat(pH, Cl2:N), PRIOR_SD_LN^2), uniform elsewhere."""
    kh = k_hat(pH, cl2n, prior_table_cached() if table is None else table)
    lkb = np.log(np.asarray([p[0] for p in params], dtype=float))
    return -0.5 * ((lkb - math.log(kh)) / PRIOR_SD_LN) ** 2


def make_model(kind: str, sc, X, seed: int, cache_dir: str, chem: dict | None = None):
    from .simgp import DOSE_GRID, DOSES_CA, LIK_SD, LIK_SD_CA, SimGP24
    if kind == "a":
        return SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid="full", lik_sd=LIK_SD, doses=DOSE_GRID,
                       threshold=THRESHOLD)
    m = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid="chloramine", cond=ca_condition(), lik_sd=LIK_SD_CA,
                doses=DOSES_CA, threshold=THRESHOLD)
    if kind == "c":
        m.log_prior = prior_log_vector(m.params, chem["pH"], chem["cl2n"])
    return m


def _score(truth_min: pd.Series, median: pd.Series, uns: list, flags: pd.Series, bands: pd.DataFrame | None,
           deciles: pd.Series) -> dict:
    """Sums and rates on the unsampled junctions (the daily minimum).  Recall is NaN (printed 'not testable') when no
    unsampled junction is below the threshold, precision NaN when nothing is flagged; never 1.0 by default."""
    from .simgp import FLOOR
    t, m, f = truth_min.loc[uns], median.loc[uns], flags.loc[uns].astype(bool)
    e = m - t
    eln = np.log(np.clip(m, FLOOR, None)) - np.log(np.clip(t, FLOOR, None))
    tv = t < THRESHOLD
    tp, fp, fn = int((tv & f).sum()), int((~tv & f).sum()), int((tv & ~f).sum())
    out = {"n_uns": len(uns), "sse": float((e ** 2).sum()), "sae": float(e.abs().sum()), "sse_ln": float((eln ** 2).sum()),
           "se": float(e.sum()), "n_true_viol": int(tv.sum()), "tp": tp, "fp": fp, "fn": fn}
    for q in (50, 80, 90, 95):
        out[f"in{q}"] = (float(((t >= bands.loc[uns, f"lo{q}"]) & (t <= bands.loc[uns, f"hi{q}"])).sum())
                         if bands is not None else np.nan)
    dec = deciles.loc[uns]
    for k in range(10):
        mk = (dec == k).values
        out[f"age{k + 1}_n"] = int(mk.sum())
        out[f"age{k + 1}_se"] = float(e.values[mk].sum())
    return out


def rates(d: dict) -> dict:
    """Rates from sums (one row or pooled): rmse, mae, rms_ln, bias, coverage, recall, precision, f1, false alarms."""
    n = d["n_uns"]
    rec = d["tp"] / (d["tp"] + d["fn"]) if d["tp"] + d["fn"] else float("nan")
    prec = d["tp"] / (d["tp"] + d["fp"]) if d["tp"] + d["fp"] else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if (prec == prec and rec == rec and prec + rec > 0) else float("nan")
    out = {"rmse": math.sqrt(d["sse"] / n), "mae": d["sae"] / n, "rms_ln": math.sqrt(d["sse_ln"] / n), "bias": d["se"] / n,
           "recall": rec, "precision": prec, "f1": f1, "false_alarms": d["fp"], "n_true_viol": d["n_true_viol"]}
    for q in (50, 80, 90, 95):
        out[f"coverage{q}"] = d[f"in{q}"] / n if d[f"in{q}"] == d[f"in{q}"] else float("nan")
    return out


def _edge(model) -> dict:
    from .simgp import grid_edge_mass
    e = grid_edge_mass(model.W_, model.params, model.doses)
    lkb = np.log(np.asarray([p[0] for p in model.params], dtype=float))
    return {"kb_edge": e["kb"], "kb_low": e["kb_low"], "kb_high": e["kb_high"], "kw_edge": e["kw"], "kw_low": e["kw_low"],
            "kw_high": e["kw_high"], "dose_edge": e["dose"], "dose_low": e["dose_low"], "dose_high": e["dose_high"],
            "post_kb_geo": float(np.exp(model.w_ @ lkb)), "map_kb": float(model.map_params_[0]),
            "map_kw": float(model.map_params_[1]), "map_gamma": float(model.map_params_[2]), "map_dose": float(model.map_dose_)}


def watch_flags(model, sc, temp_C: float = TRUTH_TEMP_C) -> pd.Series:
    """The nitrification watch (a rule, not a fit; literature thresholds, not validated): a junction is flagged when
    P(daily-minimum total chlorine < 0.4 mg/L) > 0.5 (critical range 0.2 to 0.65, Sathasivan, Fisher & Tam 2008), the
    month's water temperature is at least 15 C (nitrification observed from 8 to 26 C, optimum 25 to 30 C, US EPA 2002)
    and the junction's nominal daily-mean water age is at least the network's 75th percentile.  model: a SimGP24 after
    predict_daily_min."""
    age = sc.age_by_hour_h.mean()
    old = age >= age.quantile(WATCH_AGE_QUANTILE)
    p = model.p_below_mc(WATCH_CRIT_MGL)
    return (p > WATCH_P) & old & (temp_C >= WATCH_MIN_TEMP_C)


def _watch_cols(flags: pd.Series, affected: pd.Series | None) -> dict:
    n_flag = int(flags.sum())
    out = {"watch_n_flagged": n_flag}
    if affected is not None:
        n_aff, hit = int(affected.sum()), int((flags & affected).sum())
        out.update(watch_n_affected=n_aff, watch_affected_on_list=hit, watch_false_flags=n_flag - hit)
    return out


def _variant_truth(variant: str):
    if variant == "stress":
        return "epa_msx", ChloramineTruth(ph_range=PH_RANGE_STRESS)
    if variant.startswith("fm"):
        return "epa_msx", ChloramineTruth(fm=float(variant[2:]))
    if variant == "twin":
        return "first", None
    return "epa_msx", None


def run_task(net: str, seed: int, variant: str, kw_ref: float, cache_dir: str, main_min: pd.Series | None = None,
             keep_fig: bool = False) -> dict:
    """Every model and rule for one network, seed and truth variant ('main', 'stress' (low pH), 'fm3'/'fm10'/'fm30'
    (nitrification stress), 'twin' (first-order twin)).  Returns rows (sums per model, rule and n), the truth's
    summary and, with keep_fig, the first seed's figure data."""
    from .features import build_features
    from .simgp import DAY_HOURS, DOSES_CA
    from .surrogate import acquire_time, baseline_mean, baseline_nearest
    kin, truth = _variant_truth(variant)
    sc = truth_scenario(net, seed, kw_ref, kin, truth)
    X = build_features(sc)
    J = list(sc.junctions); jidx = {j: i for i, j in enumerate(J)}
    tbh, tmin = sc.truth_by_hour, sc.truth_daily_min
    chem = {k: v for k, v in sc.chem.items() if k != "_volatile"}
    age_mean = sc.age_by_hour_h.mean()
    deciles = pd.Series(pd.qcut(age_mean.rank(method="first"), 10, labels=False), index=age_mean.index)
    rng = np.random.default_rng(50_000 + seed)
    noise = rng.normal(0.0, NOISE_SD, size=(24, len(J)))
    def read(j, h):
        return float(np.clip(tbh.loc[h, j] + noise[h, jidx[j]], 0.01, None))
    init_j = [str(x) for x in rng.choice(J, N0, replace=False)]
    init_h = [int(h) for h in rng.choice(DAY_HOURS, N0)]
    rest = [j for j in J if j not in init_j]
    rand_j = init_j + [str(x) for x in rng.permutation(rest)[:NMAX - N0]]
    rand_h = init_h + [int(h) for h in rng.choice(DAY_HOURS, NMAX - N0)]
    def frame(js, hs):
        return pd.DataFrame({"junction": js, "hour": hs, "y": [read(j, h) for j, h in zip(js, hs)]})
    kinds = {"main": ("a", "b", "c"), "stress": ("a", "b", "c"), "twin": ("b",)}.get(variant, ("c",))
    models = {k: make_model(k, sc, X, seed, cache_dir, chem) for k in kinds}
    affected = None
    if main_min is not None:
        affected = (main_min.loc[J] - tmin.loc[J]) >= AFFECTED_DROP_MGL
    rows, fig = [], {}
    base = {"net": net, "variant": variant, "seed": seed}

    def grid_row(kind, rule, n, model, pmin, S):
        uns = [j for j in J if j not in set(S.junction)]
        r = {**base, "model": "g" if variant == "twin" else kind, "rule": rule, "n": n,
             **_score(tmin, pmin["median"], uns, pmin["p_below"] > 0.5, pmin, deciles), **_edge(model)}
        if kind == "c":
            r.update(_watch_cols(watch_flags(model, sc), affected))
        return r

    def baseline_rows(rule, n, S):
        uns = [j for j in J if j not in set(S.junction)]
        out = []
        for kind, med in (("d", baseline_mean(sc, list(S.junction), S.y.values)),
                          ("e", baseline_nearest(sc, list(S.junction), S.y.values))):
            out.append({**base, "model": kind, "rule": rule, "n": n,
                        **_score(tmin, med, uns, med < THRESHOLD, None, deciles)})
        return out

    orc = best_member(models[kinds[-1]].Z, tbh, DOSES_CA) if variant in ("main", "stress", "twin") else None
    if orc is not None:
        o_med = pd.Series(np.exp(orc["z"].min(axis=0)), index=J)
    # random rule: one sample set for every model
    for n in NS:
        S = frame(rand_j[:n], rand_h[:n])
        for kind, m in models.items():
            m.fit(S); pmin = m.predict_daily_min()
            rows.append(grid_row(kind, "random", n, m, pmin, S))
            if keep_fig and n == 8 and kind in ("a", "c"):
                fig[f"random_{kind}"] = {"median": pmin["median"].copy(), "p_below": pmin["p_below"].copy(), "S": S.copy()}
                if kind == "c":
                    fig["random_c"]["watch"] = watch_flags(m, sc)
        if variant in ("main", "stress"):
            rows += baseline_rows("random", n, S)
        if orc is not None:
            uns = [j for j in J if j not in set(S.junction)]
            rows.append({**base, "model": "f", "rule": "random", "n": n,
                         **_score(tmin, o_med, uns, o_med < THRESHOLD, None, deciles)})
    # straddle on the daily minimum: every grid model chooses its own samples from the same three first ones
    for kind, m in models.items():
        S = frame(init_j, init_h)
        for n in range(N0, NMAX + 1):
            m.fit(S); hourly = m.predict_hours(); pmin = m.predict_daily_min()
            if n in NS:
                rows.append(grid_row(kind, "straddle_min", n, m, pmin, S))
                if kind == kinds[-1] and variant in ("main", "stress"):
                    rows += baseline_rows("straddle_min", n, S)      # the baselines read the main model's samples
                if kind == kinds[-1] and orc is not None:
                    uns = [j for j in J if j not in set(S.junction)]
                    rows.append({**base, "model": "f", "rule": "straddle_min", "n": n,
                                 **_score(tmin, o_med, uns, o_med < THRESHOLD, None, deciles)})
            if n == NMAX:
                break
            uns = [j for j in J if j not in set(S.junction)]
            hd = (pd.DataFrame(hourly[0], columns=J), pd.DataFrame(hourly[1], columns=J), pd.DataFrame(m.z_sd_acq_, columns=J))
            j, h = acquire_time("straddle_min", hd, pmin, uns, DAY_HOURS, None, THRESHOLD)
            S = pd.concat([S, frame([j], [h])], ignore_index=True)
    tsum = {**base, **{k: chem.get(k) for k in ("pH", "cl2n", "toc_mgL", "alk_mgL_caco3", "k_hat_per_day",
                                                "kb_twin_per_day", "twin_source_factor", "fm", "n_fm_pipes",
                                                "fm_age_threshold_h", "msx_compiler")},
            "kw_ref_m_per_day": kw_ref, "dose_mean_mgL": float(np.mean(list(chem["source_doses_mgL"].values()))),
            "n_junctions": len(J), "n_below_0.5_min": int((tmin < THRESHOLD).sum()),
            "n_below_0.4_min": int((tmin < WATCH_CRIT_MGL).sum()), "n_below_0.2_min": int((tmin < 0.2).sum()),
            "median_min_mgL": float(tmin.median()), "min_min_mgL": float(tmin.min()),
            "n_affected": int(affected.sum()) if affected is not None else None,
            "tank_min_last_day_mean_mgL": (min(chem["tank_last_day_mean_mgL"].values())
                                           if chem.get("tank_last_day_mean_mgL") else None),
            "n_tanks_below_0.4": (int(sum(v < WATCH_CRIT_MGL for v in chem["tank_last_day_mean_mgL"].values()))
                                  if chem.get("tank_last_day_mean_mgL") else None),
            "oracle_rms_ln_day": orc["rms_ln"] if orc is not None else None}
    if keep_fig:
        fig.update(truth_min=tmin.copy(), main_min=main_min.copy() if main_min is not None else None,
                   affected=affected.copy() if affected is not None else None, age_mean=age_mean.copy())
    return {"rows": rows, "truth": tsum, "fig": fig, "truth_min": tmin.copy(),
            "seconds": sc.chem.get("_volatile", {}).get("msx_seconds")}


def _task(args):
    return run_task(*args)


class _single_thread_children:
    """One BLAS and OpenMP thread in every worker: the variables must be in the environment the workers start with
    (Apple's Accelerate reads them when numpy loads, before a pool initializer runs).  Restored afterwards."""
    def __enter__(self):
        self.saved = {k: os.environ.get(k) for k in SINGLE_THREAD_ENV}
        os.environ.update(SINGLE_THREAD_ENV)

    def __exit__(self, *exc):
        for k, v in self.saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _pool(tasks, workers: int):
    import tempfile
    from concurrent.futures import ProcessPoolExecutor
    with tempfile.TemporaryDirectory(prefix="rm_ca_exp_") as tmp:
        if workers <= 1:
            old = os.getcwd()
            try:
                _pool_init(tmp)
                return [_task(t) for t in tasks]
            finally:
                os.chdir(old)
        with _single_thread_children(), ProcessPoolExecutor(max_workers=workers, initializer=_pool_init,
                                                            initargs=(tmp,)) as ex:
            return list(ex.map(_task, tasks))


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if math.isnan(float(x)) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def _fmt(v, nd=3):
    """A rate for the docs: 'not testable' when it is undefined (no junction below the threshold, or nothing flagged)."""
    return "not testable" if v is None or (isinstance(v, float) and math.isnan(v)) else round(float(v), nd)


def pooled(df: pd.DataFrame) -> dict:
    """Rates from sums pooled over the rows' seeds (and age-decile biases)."""
    s = {c: float(df[c].sum(min_count=1)) if df[c].notna().any() else float("nan") for c in SUM_COLS}
    out = rates(s)
    out.update(n_seeds=int(df.seed.nunique()), n_junction_days=int(s["n_uns"]))
    for c in ("kb_edge", "kb_low", "kb_high", "kw_edge", "kw_low", "kw_high", "dose_edge", "post_kb_geo"):
        if c in df and df[c].notna().any():
            out[f"mean_{c}"] = float(df[c].mean())
    for c in ("watch_n_flagged", "watch_n_affected", "watch_affected_on_list", "watch_false_flags"):
        if c in df and df[c].notna().any():
            out[c] = int(df[c].sum())
    if "watch_n_affected" in out:
        out["watch_share_affected_on_list"] = (out["watch_affected_on_list"] / out["watch_n_affected"]
                                               if out["watch_n_affected"] else float("nan"))
        out["watch_share_false_flags"] = (out["watch_false_flags"] / out["watch_n_flagged"]
                                          if out["watch_n_flagged"] else float("nan"))
    return out


def age_bias(df: pd.DataFrame) -> list:
    """Mean signed error (mg/L) of the daily-minimum prediction by nominal water-age decile (1 = youngest)."""
    out = []
    for k in range(1, 11):
        n = float(df[f"age{k}_n"].sum())
        out.append(float(df[f"age{k}_se"].sum()) / n if n else float("nan"))
    return out


def summary_tree(df: pd.DataFrame) -> dict:
    tree = {}
    for (v, m, r, n), g in df.groupby(["variant", "model", "rule", "n"]):
        tree.setdefault(v, {}).setdefault(m, {}).setdefault(r, {})[str(int(n))] = pooled(g)
    return tree


def acceptance(net: str, df: pd.DataFrame, P: dict) -> dict:
    """The pre-registered bars (journal, task 12), from the results CSV as written."""
    main = df[df.variant == "main"]
    A = {}
    if net == "Net3":
        rec = P["main"]["c"]["straddle_min"]["8"]["recall"]
        A["A1_recall_c_straddle_n8"] = {"value": _fmt(rec), "bar": ">= 0.80",
                                        "pass": bool(rec == rec and rec >= 0.80),
                                        "counts": f"{int(P['main']['c']['straddle_min']['8']['n_true_viol'])} low junction-days"}
    cov = [P["main"]["c"]["random"][str(n)]["coverage90"] for n in NS]
    cv = float(np.mean(cov))
    A["A2_coverage90_c_random_mean_over_n"] = {"value": round(cv, 3), "by_n": [round(c, 3) for c in cov],
                                               "bar": "in [0.85, 0.97]", "pass": bool(0.85 <= cv <= 0.97)}
    e15 = main[(main.model == "c") & (main.n == 15)]
    kb, kw = float(e15.kb_edge.mean()), float(e15.kw_edge.mean())
    A["A3_edge_mass_c_n15"] = {"kb": round(kb, 3), "kb_low": round(float(e15.kb_low.mean()), 3),
                               "kb_high": round(float(e15.kb_high.mean()), 3), "kw": round(kw, 3),
                               "kw_low": round(float(e15.kw_low.mean()), 3), "kw_high": round(float(e15.kw_high.mean()), 3),
                               "bar": "each <= 0.5", "pass": bool(kb <= 0.5 and kw <= 0.5)}
    if net == "Net3":
        per = {}
        for s, g in main.groupby("seed"):
            rc = float(np.mean([rates(r)["rmse"] for _, r in g[g.model == "c"].iterrows()]))
            ra = float(np.mean([rates(r)["rmse"] for _, r in g[g.model == "a"].iterrows()]))
            per[str(int(s))] = {"rmse_c": round(rc, 4), "rmse_a": round(ra, 4), "c_lower": bool(rc < ra)}
        k = sum(v["c_lower"] for v in per.values())
        A["A4_c_beats_free_grid_rmse_seeds"] = {"value": f"{k} of {len(per)}", "bar": ">= 6 of 8", "pass": bool(k >= 6),
                                                "per_seed": per}
    a15 = main[(main.model == "a") & (main.n == 15)]
    cov_a = float(np.mean([P["main"]["a"]["random"][str(n)]["coverage90"] for n in NS]))
    A["A5_free_grid_reported"] = {"kb_floor_mass_n15": round(float(a15.kb_low.mean()), 3),
                                  "kb_edge_n15": round(float(a15.kb_edge.mean()), 3), "coverage90_random": round(cov_a, 3),
                                  "expected_failure_seen": bool(float(a15.kb_low.mean()) > 0.5 or cov_a < 0.7),
                                  "bar": "none (reported either way)"}
    model_bars = [v["pass"] for k, v in A.items() if k[:2] in ("A1", "A2", "A3", "A4")]
    A["experimental"] = not all(model_bars)
    return A


def acceptance_model_b(net: str, df: pd.DataFrame) -> dict:
    """Reported, no bar (added after the review): the bars, where they apply, applied to model (b), the chloramine grid
    with a uniform prior, which the pilot runs when the plant's pH and Cl2:N are not given.  From the same CSV."""
    d = pd.concat([df[df.model != "c"], df[df.model == "b"].assign(model="c")], ignore_index=True)
    A = acceptance(net, d, summary_tree(d))
    out = {"about": "model (b), uniform prior: the pre-registered bars applied for information only; they were set "
                    "for model (c)"}
    for k, v in A.items():
        if k[:2] in ("A1", "A2", "A3", "A4"):
            out[k.replace("_c_", "_b_")] = {kk: vv for kk, vv in v.items() if kk != "per_seed"}
    out["all_applicable_met"] = all(v["pass"] for k, v in out.items() if k[:2] in ("A1", "A2", "A3", "A4"))
    return out


def _settings(net: str) -> dict:
    from .simgp import DOSES_CA, GRIDS, LIK_SD, LIK_SD_CA
    return {"threshold_mgL": THRESHOLD, "threshold_note": THRESHOLD_NOTE, "species": "total chlorine (mg/L as Cl2)",
            "main_seeds": list(range(FIRST_SEED, FIRST_SEED + N_MAIN)), "stress_seeds": list(STRESS_SEEDS.get(net, ())),
            "fm_sweep": list(FM_SWEEP) if net in FM_NETS else [0.0], "kw_ref_m_per_day": KW_REF[net],
            "nominal_dose_mgL": CA_DOSE_MGL, "truth_ranges": {"pH": PH_RANGE, "pH_stress": PH_RANGE_STRESS,
                                                              "cl2n": CL2N_RANGE, "toc_mgL": TOC_RANGE,
                                                              "alk_mgL_caco3": ALK_RANGE, "temp_C": TRUTH_TEMP_C},
            "grid_chloramine": [list(map(float, a)) for a in GRIDS["chloramine"]], "doses_chloramine": DOSES_CA,
            "lik_sd": {"a": LIK_SD, "b_c_g": LIK_SD_CA}, "prior_sd_ln": PRIOR_SD_LN, "ns": list(NS), "rules": list(RULES),
            "noise_sd_mgL": NOISE_SD, "models": MODELS, "affected_drop_mgL": AFFECTED_DROP_MGL,
            "watch": {"p_below_0.4": WATCH_P, "min_temp_C": WATCH_MIN_TEMP_C, "age_quantile": WATCH_AGE_QUANTILE,
                      "label": "literature thresholds, not validated"}}


def summarise(net: str, df: pd.DataFrame, truths: pd.DataFrame) -> dict:
    P = summary_tree(df)
    out = {"network": net, "about": "Task 12, chloramine mode, simulation only. Every number is total chlorine "
           "(mg/L as Cl2); model (a) is the free-chlorine grid and settings run on the chloraminated network (the "
           "conflation test). Rates are pooled over seeds from the sums in results_time_<net>.csv; a rate with no "
           "junction below the threshold (or nothing flagged) is null and printed 'not testable'.",
           "settings": _settings(net), "acceptance": acceptance(net, df, P),
           "acceptance_model_b_no_bar": acceptance_model_b(net, df), "pooled": P}
    out["age_decile_bias_mgL"] = {f"{v}/{m}": [round(x, 4) if x == x else None for x in age_bias(g)]
                                  for (v, m), g in df[df.n == 8].groupby(["variant", "model"])}
    tr = truths.copy()
    out["truths"] = {v: {"n": int(len(g)), "pH": [round(float(g.pH.min()), 2), round(float(g.pH.max()), 2)],
                         "share_below_0.5_min": round(float(g["n_below_0.5_min"].sum() / g.n_junctions.sum()), 3),
                         "share_below_0.2_min": round(float(g["n_below_0.2_min"].sum() / g.n_junctions.sum()), 3),
                         "median_daily_min_mgL": round(float(g.median_min_mgL.median()), 3)}
                     for v, g in tr.groupby("variant")}
    if (df.variant == "twin").any():
        out["first_order_cost"] = {r: {str(n): {"b_on_msx": {k: _fmt(P["main"]["b"][r][str(n)][k]) for k in ("rmse", "rms_ln", "recall", "coverage90")},
                                                "b_on_twin": {k: _fmt(P["twin"]["g"][r][str(n)][k]) for k in ("rmse", "rms_ln", "recall", "coverage90")}}
                                       for n in NS} for r in RULES}
    if any(str(v).startswith("fm") for v in df.variant.unique()):
        nit = {}
        for v in ["main"] + [f"fm{int(f)}" for f in FM_SWEEP if f > 0]:
            if v not in P:
                continue
            nit[v] = {r: {str(n): {k: _fmt(P[v]["c"][r][str(n)].get(k), 3) for k in
                                   ("recall", "coverage90", "rmse", "n_true_viol", "watch_n_flagged", "watch_n_affected",
                                    "watch_affected_on_list", "watch_share_affected_on_list", "watch_false_flags",
                                    "watch_share_false_flags")} for n in NS} for r in RULES}
        out["nitrification_sweep"] = {"label": "literature thresholds, not validated; no bar", "by_fm": nit}
    return out


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


BIG = {"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 11, "xtick.labelsize": 11,
       "ytick.labelsize": 11}


def plot_conflation(net: str, summ: dict, fig_data: dict, wn, out: str, seed: int | None = None) -> None:
    """Free-chlorine settings on a chloraminated network against the chloramine mode: error, band and recall against
    the number of daytime samples, and one scenario's map.  Every number is total chlorine."""
    import wntr
    plt = _plt()
    P = summ["pooled"]["main"]
    style = {"a": ("tab:orange", "free-chlorine grid and settings (wrong mode)"),
             "b": ("tab:blue", "chloramine grid, uniform prior"),
             "c": ("tab:purple", "chloramine mode (pH and Cl2:N prior)"),
             "d": ("gray", "mean of samples")}
    with plt.rc_context(BIG):
        fig = plt.figure(figsize=(31, 6.8))
        gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1.1, 1.1, 1.1])
        axes = [fig.add_subplot(gs[i]) for i in range(3)]
        for m, (col, lab) in style.items():
            for rule, ls in (("random", "--"), ("straddle_min", "-")):
                if m == "d" and rule == "straddle_min":
                    continue
                g = P[m][rule]
                y0 = [g[str(n)]["rmse"] for n in NS]
                axes[0].plot(NS, y0, ls, color=col, marker="o", label=f"{lab}, {'random' if rule == 'random' else 'straddle'}")
                if m != "d":
                    axes[1].plot(NS, [g[str(n)]["coverage90"] for n in NS], ls, color=col, marker="o")
                axes[2].plot(NS, [g[str(n)]["recall"] if g[str(n)]["recall"] is not None else np.nan for n in NS], ls,
                             color=col, marker="o")
        axes[1].axhspan(0.85, 0.97, color="green", alpha=0.08); axes[1].axhline(0.9, color="k", lw=0.8, ls=":")
        axes[2].axhline(0.8, color="k", lw=0.8, ls=":")
        axes[0].set(title="Daily-minimum error, unsampled junctions", xlabel="daytime samples", ylabel="RMSE (mg/L total chlorine)")
        axes[1].set(title="Truth inside the 90% band", xlabel="daytime samples", ylabel="coverage", ylim=(0, 1.02))
        axes[2].set(title=f"Recall of junctions below {THRESHOLD} mg/L", xlabel="daytime samples", ylabel="recall", ylim=(0, 1.02))
        for ax in axes:
            ax.grid(alpha=0.3); ax.set_xticks(NS)
        axes[0].legend(fontsize=9, loc="upper right")
        ra, rc = fig_data.get("random_a"), fig_data.get("random_c")
        tmin = fig_data["truth_min"]
        uns = [j for j in tmin.index if j not in set(ra["S"].junction)]
        tv = tmin.loc[uns] < THRESHOLD
        tag = f"scenario {seed}, " if seed is not None else ""
        for k, (key, title) in enumerate((("truth", "TRUE daily minimum (hidden)"), ("a", "FREE-chlorine settings"),
                                          ("c", "CHLORAMINE mode"))):
            ax = fig.add_subplot(gs[3 + k])
            if key == "truth":
                ser = tmin
                t = f"{title}\n{tag}{int((tmin < THRESHOLD).sum())} junctions below {THRESHOLD} mg/L total chlorine"
            else:
                d = ra if key == "a" else rc
                ser = d["median"]
                e = float(np.sqrt(np.mean((d["median"].loc[uns] - tmin.loc[uns]) ** 2)))
                fl = d["p_below"].loc[uns] > 0.5
                t = (f"{title}\n{tag}8 random daytime samples:\nRMSE {e:.2f} mg/L, found {int((fl & tv).sum())} of "
                     f"{int(tv.sum())} low, {int((fl & ~tv).sum())} false alarms")
            wntr.graphics.plot_network(wn, node_attribute=ser.to_dict(), node_size=45, node_cmap="viridis",
                                       node_range=(0, 2.2), ax=ax, link_width=0.6, add_colorbar=True, title=t)
            ax.set_title(t, fontsize=12)
            if key != "truth":
                S = ra["S"]
                ax.scatter([wn.get_node(j).coordinates[0] for j in S.junction], [wn.get_node(j).coordinates[1] for j in S.junction],
                           s=130, facecolors="none", edgecolors="cyan", linewidths=1.8, zorder=5)
        A = summ["acceptance"]
        fig.suptitle(f"{net}, chloraminated (EPA chloramine chemistry in EPANET-MSX, simulated): free-chlorine settings vs the "
                     f"chloramine mode; all values total chlorine, threshold {THRESHOLD} mg/L (a common utility operating "
                     f"target, not a California rule)" + ("" if not A.get("experimental") else "; chloramine mode: experimental")
                     + ("\nleft: pooled over the scenarios; right: the first scenario only (chosen before the run), which can "
                        "differ from the pooled result" if seed is not None else ""), fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_nitrification(summ: dict, fig_data: dict, wn, out: str, seed: int | None = None) -> None:
    """The nitrification sweep (no bar): recall and coverage of the chloramine mode at 8 samples against the stress
    multiplier, the watch list against the affected junctions, and one scenario's map at the strongest stress."""
    import wntr
    plt = _plt()
    nit = summ["nitrification_sweep"]["by_fm"]
    keys = [k for k in ["main", "fm3", "fm10", "fm30"] if k in nit]
    fms = [0.0 if k == "main" else float(k[2:]) for k in keys]
    def get(k, r, f):
        v = nit[k][r]["8"][f]
        return np.nan if isinstance(v, str) else v
    with plt.rc_context(BIG):
        fig = plt.figure(figsize=(24, 6.6))
        gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.2, 1.2])
        ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        xs = np.arange(len(keys))
        for r, ls in (("random", "--"), ("straddle_min", "-")):
            ax0.plot(xs, [get(k, r, "recall") for k in keys], ls, color="tab:purple", marker="o", label=f"recall, {r}")
            ax0.plot(xs, [get(k, r, "coverage90") for k in keys], ls, color="tab:green", marker="s", label=f"90% coverage, {r}")
        ax0.set(xticks=xs, xticklabels=[f"Fm {f:g}" for f in fms], ylim=(0, 1.02), ylabel="rate",
                title="Chloramine mode at 8 samples\nas the nitrification stress grows")
        ax0.grid(alpha=0.3); ax0.legend(fontsize=9, loc="lower left")
        w = 0.38
        aff = [get(k, "straddle_min", "watch_n_affected") for k in keys]
        hit = [get(k, "straddle_min", "watch_affected_on_list") for k in keys]
        fl = [get(k, "straddle_min", "watch_n_flagged") for k in keys]
        ax1.bar(xs - w / 2, aff, w, color="tab:red", label="affected junctions (truth)")
        ax1.bar(xs + w / 2, fl, w, color="tab:gray", label="on the watch list")
        ax1.bar(xs + w / 2, hit, w, color="tab:olive", label="on the list and affected")
        ax1.set(xticks=xs, xticklabels=[f"Fm {f:g}" for f in fms], ylabel="junction-days, 8 seeds",
                title="Nitrification watch (literature thresholds,\nnot validated), straddle rule, 8 samples")
        ax1.grid(alpha=0.3, axis="y"); ax1.legend(fontsize=9)
        f = fig_data
        if f.get("main_min") is not None and f.get("random_c") is not None:
            drop = (f["main_min"] - f["truth_min"]).clip(lower=0)
            ax2 = fig.add_subplot(gs[2])
            wntr.graphics.plot_network(wn, node_attribute=drop.to_dict(), node_size=45, node_cmap="Reds",
                                       node_range=(0, max(0.1, float(drop.max()))), ax=ax2, link_width=0.6, add_colorbar=True,
                                       title=(f"scenario {seed}, " if seed is not None else "") +
                                             f"Fm {fms[-1]:g}: fall of the true daily minimum (mg/L)\n"
                                             f"{int(f['affected'].sum())} junctions affected (fall >= {AFFECTED_DROP_MGL} mg/L)")
            ax3 = fig.add_subplot(gs[3])
            wl = f["random_c"]["watch"]
            cat = pd.Series(0.15, index=wl.index)
            cat[wl] = 0.85
            wntr.graphics.plot_network(wn, node_attribute=cat.to_dict(), node_size=45, node_cmap="Purples", node_range=(0, 1),
                                       ax=ax3, link_width=0.6, add_colorbar=False,
                                       title=(f"scenario {seed}, " if seed is not None else "") +
                                             f"Fm {fms[-1]:g}: watch list (dark),\n8 random samples: "
                                             f"{int(wl.sum())} flagged, {int((wl & f['affected']).sum())} affected\n"
                                             f"red rings: all {int(f['affected'].sum())} affected junctions")
            aff = [j for j in wl.index if f["affected"].get(j, False)]
            ax3.scatter([wn.get_node(j).coordinates[0] for j in aff], [wn.get_node(j).coordinates[1] for j in aff], s=110,
                        facecolors="none", edgecolors="tab:red", linewidths=1.4, zorder=5)
        fig.suptitle("Nitrification stress sweep (a sweep, not a fit: no typical Fm is published); total chlorine; "
                     "simulated", fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def write_outputs(net: str, df: pd.DataFrame, truths: pd.DataFrame, outdir: str, fig_main: dict | None = None,
                  fig_fm: dict | None = None, seed: int | None = None) -> dict:
    from .simulate import load
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"results_time_{net}.csv"), index=False, float_format=CSV_FLOAT)
    truths.to_csv(os.path.join(outdir, f"truths_{net}.csv"), index=False, float_format=CSV_FLOAT)
    # the summary is computed from the files as written
    df = pd.read_csv(os.path.join(outdir, f"results_time_{net}.csv"))
    truths = pd.read_csv(os.path.join(outdir, f"truths_{net}.csv"))
    summ = summarise(net, df, truths)
    with open(os.path.join(outdir, f"summary_chloramine_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    wn = load(net)
    if fig_main:
        plot_conflation(net, summ, fig_main, wn, os.path.join(outdir, f"chloramine_mode_vs_free_grid_{net}.png"), seed)
    if fig_fm:
        plot_nitrification(summ, fig_fm, wn, os.path.join(outdir, f"nitrification_stress_{net}.png"), seed)
    return summ


def replot(net: str, outdir: str = OUT_DIR, cache_dir: str = "outputs/cache") -> dict:
    """Redraw the figures from the committed summary: the first scored seed's main task (and, on Net3, its Fm 30 task)
    is recomputed for the map panels, and its rows must equal the committed CSV's to its 6 significant digits."""
    cache_dir = os.path.abspath(cache_dir)
    from .simulate import load
    s0 = FIRST_SEED
    with open(os.path.join(outdir, f"summary_chloramine_{net}.json")) as fh:
        summ = json.load(fh)
    want = pd.read_csv(os.path.join(outdir, f"results_time_{net}.csv"), float_precision="round_trip")
    m = run_task(net, s0, "main", KW_REF[net], cache_dir, None, True)
    tasks = [("main", m)]
    f = None
    if net in FM_NETS:
        f = run_task(net, s0, f"fm{int(max(FM_SWEEP))}", KW_REF[net], cache_dir, m["truth_min"], True)
        tasks.append((f"fm{int(max(FM_SWEEP))}", f))
    n_cmp = 0
    for variant, res in tasks:
        got = pd.DataFrame(res["rows"])
        w = want[(want.variant == variant) & (want.seed == s0)]
        if len(got) != len(w):
            raise RuntimeError(f"{variant}: {len(got)} rows recomputed, {len(w)} committed")
        for (_, a), (_, b) in zip(got.iterrows(), w.iterrows()):
            for c in w.columns:
                va, vb = a.get(c), b[c]
                if isinstance(vb, str) or isinstance(va, str):
                    ok = str(va) == str(vb)
                else:
                    ok = (pd.isna(va) and pd.isna(vb)) or float(CSV_FLOAT % float(va)) == float(vb)
                if not ok:
                    raise RuntimeError(f"{variant} row {a.model}/{a.rule}/{a.n} column {c}: {va} recomputed, {vb} committed")
                n_cmp += 1
    wn = load(net)
    plot_conflation(net, summ, m["fig"], wn, os.path.join(outdir, f"chloramine_mode_vs_free_grid_{net}.png"), s0)
    if f is not None:
        plot_nitrification(summ, f["fig"], wn, os.path.join(outdir, f"nitrification_stress_{net}.png"), s0)
    return {"values_compared": n_cmp}


def run_experiment(net: str, outdir: str = OUT_DIR, cache_dir: str = "outputs/cache", workers: int = 4,
                   seeds=None, smoke: bool = False) -> dict:
    """The task-12 experiment for one network.  smoke=True runs SMOKE_SEED only (every variant) and is never scored."""
    from .simgp import disk_preflight
    cache_dir = os.path.abspath(cache_dir)        # the workers run in private working directories
    disk_preflight([cache_dir, outdir])
    seeds = [SMOKE_SEED] if smoke else list(seeds or range(FIRST_SEED, FIRST_SEED + N_MAIN))
    if not smoke:
        check_seeds(seeds)
    stress = [] if smoke else list(STRESS_SEEDS.get(net, ()))
    if smoke and net in STRESS_SEEDS:
        stress = [SMOKE_SEED + 1]
    kw = KW_REF[net]
    ca_grid(net, cache_dir)                       # build or read both grids once, before the workers start
    from .simgp import simulator_grid_24h
    simulator_grid_24h(nominal_ca(net), cache_dir, "full")
    tasks = [(net, s, "main", kw, cache_dir, None, s == seeds[0]) for s in seeds]
    tasks += [(net, s, "stress", kw, cache_dir, None, False) for s in stress]
    tasks += [(net, s, "twin", kw, cache_dir, None, False) for s in seeds]
    t0 = time.time()
    res1 = _pool(tasks, workers)
    print(f"{net}: {len(tasks)} main, stress and twin tasks in {time.time() - t0:.0f} s", flush=True)
    main_min = {r["truth"]["seed"]: r["truth_min"] for r in res1 if r["truth"]["variant"] == "main"}
    res2 = []
    if net in FM_NETS:
        t1 = time.time()
        tasks2 = [(net, s, f"fm{int(f)}", kw, cache_dir, main_min[s], s == seeds[0] and f == max(FM_SWEEP))
                  for f in FM_SWEEP if f > 0 for s in seeds]
        res2 = _pool(tasks2, workers)
        print(f"{net}: {len(tasks2)} nitrification tasks in {time.time() - t1:.0f} s", flush=True)
    res = res1 + res2
    df = pd.DataFrame([r for x in res for r in x["rows"]])
    truths = pd.DataFrame([x["truth"] for x in res])
    fig_main = next((x["fig"] for x in res1 if x["fig"] and x["truth"]["variant"] == "main"), None)
    fig_fm = next((x["fig"] for x in res2 if x["fig"]), None)
    summ = write_outputs(net, df, truths, outdir, fig_main, fig_fm, seeds[0])
    secs = [x["seconds"] for x in res if x["seconds"]]
    print(f"{net}: done in {time.time() - t0:.0f} s; MSX runs {len(secs)}, {np.mean(secs) if secs else 0:.1f} s each on "
          f"average", flush=True)
    return summ


def write_calibration(outdir: str = OUT_DIR, cache_dir: str = "outputs/cache", workers: int = 4) -> dict:
    """Regenerates outputs/chloramine/kb_prior_table.csv and calibration_chloramine.json: the kw_ref rule on seed 0 and
    the LIK_SD rule on seeds 100 to 103 (the settings the experiment uses, KW_REF and simgp.LIK_SD_CA)."""
    os.makedirs(outdir, exist_ok=True)
    kb_prior_table().to_csv(os.path.join(outdir, "kb_prior_table.csv"), index=False, float_format=CSV_FLOAT)
    kw = tune_kw_ref(workers=workers)
    lk = calibrate_lik_sd({n: kw[n]["kw_ref"] for n in kw}, cache_dir=cache_dir, workers=workers)
    out = {"about": "Task 12 settings rules, run on seeds never scored: kw_ref on seed 0 (share of junctions with a daily "
                    "minimum below 0.5 mg/L total chlorine in [0.10, 0.30], closest to 0.20, ties to the smaller); "
                    "LIK_SD_CA on Net3 seeds 100 to 103 (best member's daytime RMS ln mismatch, mean, rounded up to a "
                    "multiple of 0.05). Simulation only.",
           "kw_ladder_m_per_day": list(KW_LADDER), "kw_ref": kw, "lik_sd": lk,
           "wahman_example": wahman_example()}
    with open(os.path.join(outdir, "calibration_chloramine.json"), "w") as fh:
        json.dump(_clean(out), fh, indent=1)
    return out


WAHMAN_INPUTS = {"dose_mgL": 4.0, "cl2n": 5.0, "temp_C": 25.0, "alk_mgL_caco3": 50.0, "toc_mgL": 0.0, "mode": "simadd",
                 "days": 10.0}
WAHMAN_PUBLISHED = {7.0: 0.84, 9.0: 3.2}


def wahman_example() -> dict:
    """Wahman 2018's worked example (4 mg/L held for 10 days falls to 0.84 mg/L at pH 7 and 3.2 at pH 9), reproduced
    by the batch port with ASSUMED inputs (Cl2:N 5, 25 C, alkalinity 50, no TOC, simultaneous addition)."""
    w = WAHMAN_INPUTS
    out = {"inputs_assumed": w, "published_mgL": {str(k): v for k, v in WAHMAN_PUBLISHED.items()}, "port_mgL": {}}
    for ph in WAHMAN_PUBLISHED:
        df = batch(ph, w["alk_mgL_caco3"], w["temp_C"], w["dose_mgL"], w["cl2n"], w["toc_mgL"], days=w["days"],
                   t_eval_days=[0.0, w["days"]], mode=w["mode"])
        out["port_mgL"][str(ph)] = float(df.total_mgL.iloc[-1])
    out["max_abs_diff_mgL"] = max(abs(out["port_mgL"][str(k)] - v) for k, v in WAHMAN_PUBLISHED.items())
    return out


def cli(argv) -> int:
    """python -m residualmap.experiment <net> --disinfectant=chloramine [--workers N] [--smoke DIR] [--resummarise]
    python -m residualmap.experiment --disinfectant=chloramine --calibrate [--workers N]"""
    import argparse
    ap = argparse.ArgumentParser(prog="python -m residualmap.experiment", description=cli.__doc__)
    ap.add_argument("net", nargs="?", default="Net3")
    ap.add_argument("--disinfectant", required=True, choices=["chloramine", "free_chlorine"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--outdir", default=OUT_DIR)
    ap.add_argument("--cache", default=os.path.join("outputs", "cache"))
    ap.add_argument("--smoke", metavar="DIR", help="run the unscored smoke seed only, outputs to DIR")
    ap.add_argument("--calibrate", action="store_true", help="rerun the kw_ref and LIK_SD rules and the prior table")
    ap.add_argument("--resummarise", action="store_true", help="rewrite the summary from the committed CSVs (no figures)")
    ap.add_argument("--replot", action="store_true", help="redraw the figures (recomputes the first seed's map panels)")
    a = ap.parse_args(argv)
    if a.disinfectant != "chloramine":
        ap.error("without --disinfectant=chloramine run python -m residualmap.experiment <net> <seeds> as before")
    if a.calibrate:
        out = write_calibration(a.outdir, a.cache, a.workers)
        print(json.dumps(_clean({"kw_ref": {k: v["kw_ref"] for k, v in out["kw_ref"].items()},
                                 "LIK_SD_CA": out["lik_sd"]["LIK_SD_CA"]}), indent=1))
        return 0
    if a.replot:
        print(json.dumps(replot(a.net, a.outdir, a.cache)))
        return 0
    if a.resummarise:
        df = pd.read_csv(os.path.join(a.outdir, f"results_time_{a.net}.csv"))
        tr = pd.read_csv(os.path.join(a.outdir, f"truths_{a.net}.csv"))
        summ = summarise(a.net, df, tr)
        with open(os.path.join(a.outdir, f"summary_chloramine_{a.net}.json"), "w") as fh:
            json.dump(_clean(summ), fh, indent=1)
        print(json.dumps(_clean(summ["acceptance"]), indent=1))
        return 0
    summ = run_experiment(a.net, a.smoke or a.outdir, a.cache, a.workers, smoke=bool(a.smoke))
    print(json.dumps(_clean(summ["acceptance"]), indent=1))
    return 0


def example_total_log(grab_csv: str = "docs/example_grab_log.csv", taps_csv: str = "docs/example_tap_map.csv",
                      out_csv: str = "docs/example_grab_log_total.csv", net: str = "Net3", seed: int = SMOKE_SEED) -> pd.DataFrame:
    """docs/example_grab_log_total.csv, a FORMAT example for a chloraminated system: the free-chlorine example's taps,
    dates, times and notes, with total_chlorine_mgL readings from a simulated chloramine truth on Net3 (seed 900, never
    scored; one truth per month, month seed 1000 + month index) plus an N(0, 0.03) mg/L reading error.  Not data."""
    log = pd.read_csv(grab_csv, dtype={"tap_id": str})
    taps = pd.read_csv(taps_csv, dtype={"tap_id": str, "junction_id": str}).set_index("tap_id").junction_id
    month = log.date.str[:7]
    rng = np.random.default_rng(60_000 + seed)
    vals = np.full(len(log), np.nan)
    for k, m in enumerate(sorted(month.unique())):
        sc = truth_scenario(net, seed, KW_REF[net], "epa_msx", month_seed=1000 + k)
        for i in np.where(month.values == m)[0]:
            h = int(str(log.time.iloc[i])[:2])
            vals[i] = sc.truth_by_hour.loc[h, taps[log.tap_id.iloc[i]]]
    out = log.drop(columns=["free_chlorine_mgL"])
    out.insert(3, "total_chlorine_mgL", np.round(np.clip(vals + rng.normal(0, NOISE_SD, len(vals)), 0.01, None), 2))
    out.to_csv(out_csv, index=False)
    return out
