"""
seasonal.py: water temperature as a monthly input (iteration 4, journal task 10).

Today's model (simgp.SimGP24) calibrates one set of decay rates from the last three months of grab samples and uses
it for the month ahead whatever the water temperature, so its rates ride up and down with the seasons and a model
fitted in winter has no way to know that summer water decays faster.  Here the plant's monthly water temperature
(a plant log: month, temp_C) becomes an input: each grab sample is explained at its own month's temperature, and
any month, including a coming summer, is predicted at that month's temperature.

THE MODEL (SeasonalSimGP24, a subclass, so SimGP24 itself is untouched)
  Hypotheses, stacked along the grid's member axis next to today's grid:
    H0   today's 675-member grid, the same in every month (temperature has no effect)
    H_E  for E/R in {5000, 8000, 12000} K: member k = (kb_k, kw_k, gamma_k, dm_k, rm_k) is read as 20 C values and
         simulated at the month's temperature T with bulk rate kb_k f(T; E), per-pipe wall rate
         kw_k roughness_factor(C_p, gamma_k) f(T; E), and EPANET's VISCOSITY and DIFFUSIVITY at T
         (chemistry.Chemistry(temp_C=T, er_K=E).sim_kwargs()).
    f(T; E) = exp[E (T - 20) / (293.15 (T + 273.15))] (Fisher, Kastl & Sathasivan 2012, Water Res 46:3293; Cejas,
    Diaz & Gonzalez 2026, Water 18:1390).  E/R values: 5000 K, the low end of the compiled bulk span (Powell et al.
    2000, 4660 to 9560 K, Water Res 34:117; AWWARF 6050 K via Liu, Reckhow & Li 2014); 8000 K, the prior centre of
    Cejas 2026 (8000 +/- 2000 K); 12000 K, about Blokker et al. 2014's 12,104 K.  Their 10-to-20 C rate ratios are
    1.83, 2.62 and 4.24 (1.8 to 4.3 times per 10 C).
    Applying f to the wall rate is an ASSUMPTION (Cejas 2026: no wall model includes temperature; Lee et al. 2014
    give only the direction).  The ablation M1b (Chemistry wall_mode 'mass_transfer_only') leaves the wall chemistry
    at 20 C, so the wall moves only through the mass-transfer coefficient.
  Prior: 1/4 on each hypothesis, then uniform over the 675 members and 5 doses ('uniform').  Sensitivity priors,
  reported and never used to choose: 'cejas' (H0 1/4, the H_E in proportion to N(E; 8000, 2000)) and 'h0_0.1'
  (H0 0.1, each H_E 0.3).
  Posterior: W[k, d] proportional to prior(k) prod_s t3((ln y_s - Z_{T_s}[k, h_s, j_s] - ln d) / 0.35), SimGP24's
  Student-t (nu 3) and scale.  The dose stays an exact ln-offset because the model is still first order; a plant
  log's monthly dose enters the same way (dose_ratio).
  Prediction for a month at temperature T*: moments from the stacked blocks at T*, the local hydraulic variance as
  committed, the discrepancy GP fitted on the residuals against each sample's own month mean (SimGP24's inputs and
  kernel), the daily minimum by SimGP24's Monte Carlo.  Every block keeps the grid's (decay-major, hydraulic-minor)
  order, so grid_dose_weights, posterior_moments, local_hydraulic_var, check_grid_order and the dose offset work
  unchanged.  Reported: kb20, the posterior mean of kb over the H_E members (a 20 C rate).
  At 20 C every H_E block is the committed grid itself, and the model reproduces SimGP24 (a saved check).
  Limits: E/R and kb20 are confounded inside one window (Cejas 2026 needs at least two temperatures 3 to 4 C apart;
  January to March spans 0.5 C), so the posterior over E stays near its prior; it is reported against the truth's
  E and never as identified.  The posterior mass on H0 is not a seasonality test.

THE BANK (covariate_bank): for each distinct temperature, the committed grid re-run under every H_E, all conditions
in one process-pool call, cached per condition in outputs/cache under chemistry.cache_tag names (simgp's tagged
grid names), or built in memory only (cache='off').

THE EXPERIMENT (python -m residualmap.seasonal <net> <seeds>; journal task 10).  Simulation only.  12-month synthetic
logs (pilot.synthetic_log with a schedule): a fixed route of 10 taps plus 3 rotating taps a month, random daytime
hours 07:00 to 17:00, N(0, 0.03) mg/L noise.  Hidden truths (simulate.build_scenario(chem=..., warming=...)):
  V1 seasonal: the plant temperature schedule (chemistry.monthly_temperature, an ASSUMPTION: 10 to 20 C), the
     truth's own bulk E/R ~ U(4660, 12104) K and wall factor theta_w^(T - 20), theta_w ~ U(1.00, 1.07) per seed
     (simulate.hidden_chem_draws), on top of the committed monthly bulk draw U(0.8, 1.2);
  V2 control: 15 C every month;
  V3 in-network warming: V1, plus every pipe and tank decaying at its own temperature between the plant's and the
     soil's (chemistry.Warming, soil curve and time constants ASSUMPTIONS); the model is told the plant temperature.
Models: B0 (today's SimGP24, previous 3 months, temperature-blind), B1 (the same, previous month only), M (the bank,
3 months), M6 (the bank, up to 6 months), M1b (bulk-only Arrhenius), M under the two sensitivity priors (V1 only),
persistence and the network mean of the window (readings only), and the oracle (the best single member and dose of
the target month's bank against the full daytime truth: the floor for any first-order model).
Tests: (1) rolling, fit the previous months, predict April to December; (2) extrapolation without refit, fit January
to March and predict July, August and September, and fit July to September and predict December; (3) in test (2) on
V1, the logged temperature off by +2 C and -2 C, in the forecast month only and in every logged month.
Scores: held-out readings at their own junction and hour (RMSE, MAE, recall below the threshold and false alarms for
all, seen and new taps; coverage 50/80/90/95 for all readings and 90 for seen and new taps) and the daily-minimum map at every junction not sampled in the 3-month
window (the same set for every model): RMSE, MAE, coverage90, recall (P > 0.5), false alarms, precision.  Pooled over
seeds and months, split into warming (April to August) and cooling (September to December); paired bootstrap over
seeds (2000 resamples) of RMSE(M) / RMSE(B0).
Outputs: outputs/chem/season_<net>.csv (one row per variant, test, temperature error, seed, month and model, with
sums), summary_season_<net>.json (pooled metrics, bootstrap intervals, the pre-registered acceptance key), and the
figures season_<net>.png, kb_by_month_<net>.png and summer_forecast_<net>.png.

TASK 10b (python -m residualmap.seasonal <net> <n> --wall; the section near the end of this file).  M stopped at its
stop rule, mostly because it assumes the wall responds to temperature as steeply as the bulk.  M2
(WallSeasonalSimGP24 over a WallBank) crosses the bulk E/R {5000, 8000, 12000} K with a wall E/R of its own
{0, 2500, 5000, 8000, 12000} K, so the samples weigh the wall's response; M (wall = bulk) and M1b (wall 0) are nested in
it and their cached blocks are reused.  Tested on fresh seeds in two truth worlds (W1: task 10's weak wall; W2: the
wall follows the truth's own bulk factor) against task 10's bars plus R1 (M2 not more than 10% worse than the better of
M and M1b); outputs/chem/season2_<net>.csv, season2_windows_<net>.csv, summary_season2_<net>.json, season2_<net>.png and
wall_posterior_<net>.png.  Tested, not adopted (journal, task 10b).  After the review, with no rerun: --wall
--resummarise rewrites the summary from the committed CSVs, and --wall --grid-edges writes
outputs/chem/task10b_grid_edges_<net>.json (M2's posterior mass at the edges of the kb and kw grids).
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import pickle
import shutil
import sys
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from .chemistry import (ER_DEFAULT_K, ER_HYPOTHESES_K, Chemistry, arrhenius, monthly_temperature, soil_temperature)
from .simgp import (DAY_HOURS, DOSE_GRID, FLOOR, GRIDS, HOURS, LIK_SD, SimGP24, _grid_member, check_grid_order,
                    disk_free_gb, disk_preflight, grid_cache_path, grid_dose_weights, grid_loglik, local_hydraulic_var,
                    posterior_moments, simulator_grid_24h)

CONTROL_TEMP_C = 15.0
VARIANTS = ("V1", "V2", "V3")
VARIANT_LABELS = {"V1": "seasonal", "V2": "control, 15 C every month", "V3": "seasonal with in-network warming"}
PRIORS = ("uniform", "cejas", "h0_0.1")
MONTHS = 12
PER_MONTH, ROTATING, NOISE_SD = 10, 3, 0.03
THRESHOLD = 0.2
LEVELS = (50, 80, 90, 95)
WARMING_MONTHS = (4, 5, 6, 7, 8)        # season classes of the target month (1-based)
COOLING_MONTHS = (9, 10, 11, 12)
ROLLING_TARGETS = tuple(range(4, 13))
EXTRAPOLATION = {"jan_mar_to_jul_sep": ((1, 2, 3), (7, 8, 9)), "jul_sep_to_dec": ((7, 8, 9), (12,))}
TEMP_ERRORS_C = (2.0, -2.0)
BOOT_N, BOOT_SEED = 2000, 2026
MIN_FREE_GB_RUN = 1.5


def season_of(month: int) -> str:
    return "warming" if month in WARMING_MONTHS else ("cooling" if month in COOLING_MONTHS else "winter")


def _tkey(temp_C) -> float:
    """Bank key of a temperature (the schedules use 0.5 C steps; rounding guards float noise)."""
    return round(float(temp_C), 3)


# ----------------------------------------------------------------------------- schedules (ASSUMPTIONS)
def plant_schedule(variant: str, months: int = MONTHS) -> list[dict]:
    """Per month (list index m - 1): {'temp_C': the plant water temperature the operator logs, 'soil_temp_C': the soil
    temperature the water warms toward in the truth, or None}.  V1 seasonal: chemistry.monthly_temperature (10 to 20 C);
    V2 control: 15 C every month; V3: V1's plant temperature plus chemistry.soil_temperature for in-network warming.
    All three are ASSUMPTIONS (no coastal California plant series was found), labelled as such wherever used."""
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {VARIANTS}")
    out = []
    for i in range(months):
        m = i % 12 + 1
        T = CONTROL_TEMP_C if variant == "V2" else monthly_temperature(m)
        out.append({"temp_C": float(T), "soil_temp_C": float(soil_temperature(m)) if variant == "V3" else None})
    return out


# ----------------------------------------------------------------------------- the bank
@dataclass
class Bank:
    """Hypothesis blocks of the grid at each temperature.  blocks[(T, E)] is the grid (members, 24, J) re-run at T under
    E/R = E; h0 is the committed grid (hypothesis H0, the same at every T).  Every block is in the grid's own member
    order (params), so stack(T) = [H0, E_1, E_2, ...] keeps the (decay-major, hydraulic-minor) order block by block."""
    params: list
    h0: np.ndarray
    blocks: dict
    ers: tuple
    wall_mode: str
    grid: str

    @property
    def hypotheses(self) -> tuple:
        return ("H0",) + tuple(f"E{e:g}" for e in self.ers)

    @property
    def n_base(self) -> int:
        return len(self.params)

    @property
    def stacked_params(self) -> list:
        return list(self.params) * len(self.hypotheses)

    @property
    def hyp_index(self) -> np.ndarray:
        return np.repeat(np.arange(len(self.hypotheses)), self.n_base)

    @property
    def temps(self) -> list:
        return sorted({t for t, _ in self.blocks})

    def blocks_at(self, temp_C) -> list:
        T = _tkey(temp_C)
        missing = [e for e in self.ers if (T, e) not in self.blocks]
        if missing:
            raise KeyError(f"the bank has no block at {T:g} C for E/R {missing}; it holds {self.temps}")
        return [self.h0] + [self.blocks[(T, e)] for e in self.ers]

    def stack(self, temp_C) -> np.ndarray:
        return np.concatenate(self.blocks_at(temp_C), axis=0)

    def merged(self, other: "Bank") -> "Bank":
        if (other.params != self.params or other.ers != self.ers or other.wall_mode != self.wall_mode
                or other.grid != self.grid or not np.array_equal(other.h0, self.h0)):
            raise ValueError("banks built on different grids or hypotheses cannot be merged")
        return type(self)(self.params, self.h0, {**self.blocks, **other.blocks}, self.ers, self.wall_mode, self.grid)


def _conditions(temps, ers, wall_mode):
    return [(_tkey(T), float(e), Chemistry(temp_C=_tkey(T), er_K=float(e), wall_mode=wall_mode))
            for T in sorted({_tkey(t) for t in temps}) for e in ers]


def _build_conditions(sc, conds, grid: str, n_jobs: int | None, preflight_paths) -> list[np.ndarray]:
    """The grid under several chemistry conditions, every run of every condition in ONE process-pool call.  Each run is
    exactly a build_grid_24h run (same arguments to simgp._grid_member, which deletes its EPANET files as soon as it
    is read).  Refused below simgp.MIN_FREE_GB free disk."""
    disk_preflight(list(preflight_paths) + [tempfile.gettempdir()])
    params = list(itertools.product(*GRIDS[grid]))
    n_jobs = n_jobs or max(1, (os.cpu_count() or 2) - 1)
    with tempfile.TemporaryDirectory(prefix="rm_bank_") as tmp:
        jobs = [(sc.wn_name, kb, kw, g, sc.source_dose, dm, rm, list(sc.junctions), os.path.join(tmp, f"c{c}_g{i}"),
                 cond.sim_kwargs())
                for c, cond in enumerate(conds) for i, (kb, kw, g, dm, rm) in enumerate(params)]
        if n_jobs > 1 and len(jobs) > 8:
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                rows = list(ex.map(_grid_member, jobs, chunksize=8))
        else:
            rows = [_grid_member(j) for j in jobs]
    Z = np.stack(rows).reshape(len(conds), len(params), len(HOURS), len(sc.junctions))
    return [Z[c] for c in range(len(conds))]


def covariate_bank(sc, temps, ers=ER_HYPOTHESES_K, wall_mode: str = "arrhenius", grid: str = "full",
                   cache_dir: str = "outputs/cache", cache: str = "readwrite", n_jobs: int | None = None) -> Bank:
    """The hypothesis bank for the given temperatures.  H0 is the committed grid (simgp.simulator_grid_24h).  At 20 C
    every factor is exactly 1, so every H_E block there IS the committed grid.  Every other (T, E) block is read from
    its tagged cache file in cache_dir (simgp.grid_cache_path(sc, cache_dir, grid, Chemistry(temp_C=T, er_K=E,
    wall_mode=wall_mode))) or built; all the blocks to build are built in one process-pool call.
    cache: 'readwrite' reads and writes the tagged files; 'read' reads them but builds missing blocks in memory only;
    'off' builds every non-20 C block in memory (the plus or minus 2 C banks of the experiment, and the saved checks,
    which may not add tagged files to the cache)."""
    if cache not in ("readwrite", "read", "off"):
        raise ValueError("cache must be 'readwrite', 'read' or 'off'")
    params, h0 = simulator_grid_24h(sc, cache_dir, grid)
    blocks, todo = {}, []
    for T, e, cond in _conditions(temps, ers, wall_mode):
        if cond.is_default:
            blocks[(T, e)] = h0
            continue
        f = grid_cache_path(sc, cache_dir, grid, cond)
        if cache != "off" and os.path.exists(f):
            with open(f, "rb") as fh:
                p, Z = pickle.load(fh)
            if list(p) != list(params) or Z.shape != h0.shape:
                raise ValueError(f"{f} does not hold this network's grid")
            blocks[(T, e)] = Z
        else:
            todo.append((T, e, cond, f))
    if todo:
        built = _build_conditions(sc, [c for _, _, c, _ in todo], grid, n_jobs, [cache_dir])
        for (T, e, cond, f), Z in zip(todo, built):
            blocks[(T, e)] = Z
            if cache == "readwrite":
                disk_preflight([cache_dir])
                with open(f, "wb") as fh:
                    pickle.dump((params, Z), fh)
    return Bank(list(params), h0, blocks, tuple(float(e) for e in ers), wall_mode, grid)


def bank_cache_files(sc, temps, ers=ER_HYPOTHESES_K, wall_mode="arrhenius", grid="full", cache_dir="outputs/cache"):
    """The tagged cache files a covariate_bank call would read or write (20 C blocks excluded: they are the committed
    grid)."""
    return [grid_cache_path(sc, cache_dir, grid, c) for _, _, c in _conditions(temps, ers, wall_mode) if not c.is_default]


# ----------------------------------------------------------------------------- the model
def hypothesis_prior(prior: str, ers=ER_HYPOTHESES_K) -> np.ndarray:
    """Prior mass on (H0, H_E...) for the named prior."""
    n = len(ers)
    if prior == "uniform":
        return np.full(n + 1, 1.0 / (n + 1))
    if prior == "h0_0.1":
        return np.array([0.1] + [0.9 / n] * n)
    if prior == "cejas":     # H_E in proportion to N(E; 8000, 2000) (Cejas, Diaz & Gonzalez 2026), H0 1/4
        g = np.exp(-0.5 * ((np.asarray(ers, float) - ER_DEFAULT_K) / 2000.0) ** 2)
        return np.concatenate([[0.25], 0.75 * g / g.sum()])
    raise ValueError(f"prior must be one of {PRIORS}")


def bank_dose_weights(z_sim, z_obs, lik_sd, lik, nu, doses, log_prior) -> tuple[np.ndarray, np.ndarray]:
    """grid_dose_weights with a log prior over members (uniform over doses)."""
    offs = np.log(np.asarray(doses, dtype=float))
    ll = np.stack([grid_loglik(z_sim, z_obs - d, lik_sd, lik, nu) for d in offs], axis=1) + log_prior[:, None]
    return np.exp(ll - logsumexp(ll)), offs


def discrepancy_gp(Xs: np.ndarray, r: np.ndarray, seed: int, smooth_hours: bool = True) -> GaussianProcessRegressor:
    """The discrepancy GP of SimGP24.fit, kernel and bounds unchanged (the nesting check compares the two models)."""
    lo = np.full(Xs.shape[1], 0.1)
    if smooth_hours:
        lo[[0, -2, -1]] = 1.0
    k = (ConstantKernel(0.1, (1e-3, 5.0))
         * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=list(zip(lo, np.full(Xs.shape[1], 20.0))), nu=1.5)
         + WhiteKernel(0.01, (1e-4, 0.5)))
    return GaussianProcessRegressor(kernel=k, n_restarts_optimizer=4, random_state=seed).fit(Xs, r)


class SeasonalSimGP24(SimGP24):
    """SimGP24 with the plant's monthly water temperature as an input (see the module docstring).

        model = SeasonalSimGP24(sc, X, bank).fit(samples, target_temp_C=19.5)
        model.predict_hours(); model.predict_daily_min()        # for the month at 19.5 C
        model.set_target(12.5); model.predict_daily_min()       # the same fit, another month

    samples: junction, hour, y (mg/L), temp_C (the plant temperature of the sample's month) and optionally
    dose_ratio (the month's logged plant dose over the model's dose; first order makes it an exact ln-offset).
    Subclasses may read another column (condition_column) with another key (_condition_key): task 11's TOC model
    (organics.TocSimGP24) reads toc_mgL, and its seasonal ablation a condition label.  For this class both are
    unchanged from task 10 (temp_C, rounded by _tkey)."""

    condition_column = "temp_C"      # the samples' column naming each sample's bank condition

    @staticmethod
    def _condition_key(value):
        return _tkey(value)

    def __init__(self, sc, X, bank: Bank, prior: str = "uniform", seed: int = 0, lik_sd: float = LIK_SD,
                 cache_dir: str = "outputs/cache", n_draws: int = 1024, lik: str = "t", nu: float = 3.0,
                 doses=DOSE_GRID, smooth_hours: bool = True, threshold: float = 0.2):
        super().__init__(sc, X, seed=seed, lik_sd=lik_sd, cache_dir=cache_dir, n_draws=n_draws, lik=lik, nu=nu,
                         grid=bank.grid, local_hydraulic=True, doses=doses, smooth_hours=smooth_hours, threshold=threshold)
        if list(self.params) != list(bank.params) or self.Z.shape != bank.h0.shape:
            raise ValueError("the bank was not built on this network's grid")
        self.bank, self.prior = bank, prior
        self.base_params = list(self.params)
        self.params = bank.stacked_params
        if self.n_hyd > 1:
            check_grid_order(self.params, self.n_hyd)
        self.hyp_index = bank.hyp_index
        p = hypothesis_prior(prior, bank.ers)
        self.prior_mass = p
        self.log_prior = None if prior == "uniform" else np.log(p[self.hyp_index] / bank.n_base)
        self.Z = None                    # the stacked blocks of the target month, set by set_target
        self.target_temp_C = None

    # ---- fitting
    def _z_at_samples(self, temps: np.ndarray, h: np.ndarray, idx: np.ndarray) -> np.ndarray:
        n0, out = self.bank.n_base, np.empty((len(self.params), len(h)))
        for T in np.unique(temps):
            mk = temps == T
            for b, Zb in enumerate(self.bank.blocks_at(T)):
                out[b * n0:(b + 1) * n0, mk] = Zb[:, h[mk], idx[mk]]
        return out

    def fit(self, samples: pd.DataFrame, target_temp_C: float, target_dose_ratio: float = 1.0) -> "SeasonalSimGP24":
        self.zmin_ = None
        if samples is None or len(samples) == 0:
            return self.fit_prior(target_temp_C, target_dose_ratio)
        col = self.condition_column
        if col not in samples:
            raise ValueError(f"samples need a {col} column: " + ("the plant water temperature of each sample's month"
                                                                  if col == "temp_C" else "each sample's bank condition"))
        idx = np.array([self.jidx[j] for j in samples.junction])
        h = samples.hour.astype(int).values
        temps = np.array([self._condition_key(t) for t in samples[col].values])
        ratio = samples["dose_ratio"].astype(float).values if "dose_ratio" in samples else np.ones(len(samples))
        if not (np.isfinite(ratio).all() and (ratio > 0).all()):
            raise ValueError("dose_ratio must be positive")
        z_obs = np.log(np.clip(samples.y.values, FLOOR, None)) - np.log(ratio)
        z_sim = self._z_at_samples(temps, h, idx)
        if self.log_prior is None:
            W, offs = grid_dose_weights(z_sim, z_obs, self.lik_sd, self.lik, self.nu, self.doses)
        else:
            W, offs = bank_dose_weights(z_sim, z_obs, self.lik_sd, self.lik, self.nu, self.doses, self.log_prior)
        w = W.sum(axis=1)
        self.W_, self.w_, self.offs_base_ = W, w, offs
        k_, d_ = np.unravel_index(int(np.argmax(W)), W.shape)
        self.map_params_, self.map_dose_ = self.params[k_], float(np.exp(offs[d_]))
        self.map_hypothesis_ = self.bank.hypotheses[self.hyp_index[k_]]
        m_s = w @ z_sim + W.sum(axis=0) @ offs          # posterior mean of ln C at each sample, at its month's temperature
        self.gp = discrepancy_gp(self._design(idx, h), z_obs - m_s, self.seed, self.smooth_hours)
        self.samples_ = samples.copy()
        return self.set_target(target_temp_C, target_dose_ratio)

    def fit_prior(self, target_temp_C: float = 20.0, target_dose_ratio: float = 1.0) -> "SeasonalSimGP24":
        """No samples: the prior over hypotheses, members and doses, no discrepancy GP."""
        pm = self.prior_mass[self.hyp_index] / self.bank.n_base
        W = np.repeat(pm[:, None] / len(self.doses), len(self.doses), axis=1)
        self.W_, self.w_ = W, W.sum(axis=1)
        self.offs_base_ = np.log(np.asarray(self.doses, dtype=float))
        self.map_params_ = self.map_dose_ = self.map_hypothesis_ = None
        self.gp = None
        self.samples_ = pd.DataFrame(columns=["junction", "hour", "y", "temp_C"])
        self.zmin_ = None
        return self.set_target(target_temp_C, target_dose_ratio)

    def set_target(self, temp_C: float, dose_ratio: float = 1.0) -> "SeasonalSimGP24":
        """Predict the month at temp_C (with the month's plant dose over the model's dose): the same posterior and GP,
        the bank's blocks at that temperature."""
        self.target_temp_C, self.target_dose_ratio = self._condition_key(temp_C), float(dose_ratio)
        self.Z = self.bank.stack(temp_C)
        self._grp_mean, self.zmin_ = None, None
        self.offs_ = self.offs_base_ + math.log(dose_ratio)
        self.m_, self.v_ = posterior_moments(self.W_, self.offs_, self.Z)
        self.hv_ = local_hydraulic_var(self.w_, self.Z, self.n_hyd) if self.n_hyd > 1 else np.zeros_like(self.v_)
        return self

    # ---- what the posterior says
    def hypothesis_posterior(self) -> dict:
        return dict(zip(self.bank.hypotheses, self.w_.reshape(len(self.bank.hypotheses), -1).sum(axis=1).tolist()))

    def posterior_columns(self) -> dict:
        """The posterior columns of an experiment row: P_<hypothesis> for every hypothesis."""
        return {f"P_{k}": v for k, v in self.hypothesis_posterior().items()}

    def rate20(self, which: int = 0) -> float:
        """Posterior mean of a 20 C rate (0 = kb, 1 = kw) over the H_E members, normalised within them."""
        wE = self.w_.reshape(len(self.bank.hypotheses), -1)[1:].sum(axis=0)
        P = np.asarray(self.base_params, dtype=float)
        return float(wE @ P[:, which] / wE.sum()) if wE.sum() > 0 else float("nan")

    def kb20(self) -> float:
        return self.rate20(0)

    def e_posterior_mean(self) -> float:
        hp = self.hypothesis_posterior()
        pe = np.array([hp[h] for h in self.bank.hypotheses[1:]])
        return float(pe @ np.asarray(self.bank.ers) / pe.sum()) if pe.sum() > 0 else float("nan")


# ----------------------------------------------------------------------------- scoring
def _reading_sums(z_mu: np.ndarray, z_sd: np.ndarray | None, test: pd.DataFrame, seen: np.ndarray, thr: float) -> dict:
    """Held-out readings at their own junction and hour.  z_sd None: a point prediction (baselines, oracle), with
    no coverage.  A reading is flagged when P(reading < thr) > 0.5 (the median below thr)."""
    y = test.y.values.astype(float)
    pred = np.exp(z_mu)
    flag = (norm.cdf((np.log(thr) - z_mu) / z_sd) > 0.5) if z_sd is not None else (pred < thr)
    out = {}
    for tag, mask in (("all", np.ones(len(y), bool)), ("seen", seen), ("new", ~seen)):
        yy, pp, ff = y[mask], pred[mask], flag[mask]
        low = yy < thr
        out.update({f"r_{tag}_n": int(mask.sum()), f"r_{tag}_sse": float(((yy - pp) ** 2).sum()),
                    f"r_{tag}_sae": float(np.abs(yy - pp).sum()), f"r_{tag}_low": int(low.sum()),
                    f"r_{tag}_tp": int((low & ff).sum()), f"r_{tag}_fp": int((~low & ff).sum())})
        for q in (LEVELS if tag == "all" else (90,)):        # every level for all readings, 90% for seen and new taps
            if z_sd is None:
                out[f"r_{tag}_in{q}"] = float("nan")
            else:
                k = norm.ppf(0.5 + q / 200)
                mu, sd = z_mu[mask], z_sd[mask]
                out[f"r_{tag}_in{q}"] = int(((yy >= np.exp(mu - k * sd)) & (yy <= np.exp(mu + k * sd))).sum())
    return out


def _map_sums(median: pd.Series, lo90, hi90, flag: pd.Series, truth_min: pd.Series, junctions, thr: float) -> dict:
    """Daily minimum at the given junctions (those not sampled in the window).  lo90/hi90 None: a point prediction."""
    t, p, f = truth_min.loc[junctions].values, median.loc[junctions].values, flag.loc[junctions].values.astype(bool)
    low = t < thr
    cov = (int(((t >= lo90.loc[junctions].values) & (t <= hi90.loc[junctions].values)).sum())
           if lo90 is not None else float("nan"))
    return {"d_n": int(len(t)), "d_sse": float(((p - t) ** 2).sum()), "d_sae": float(np.abs(p - t).sum()),
            "d_in90": cov, "d_low": int(low.sum()), "d_tp": int((low & f).sum()), "d_fp": int((~low & f).sum())}


def _nan_map() -> dict:
    return {k: float("nan") for k in ("d_n", "d_sse", "d_sae", "d_in90", "d_low", "d_tp", "d_fp")}


def oracle_member(Z: np.ndarray, truth_by_hour: pd.DataFrame, junctions, doses=DOSE_GRID, hours=DAY_HOURS):
    """The best single member and dose of a block stack against the full daytime truth (every junction, hours 07:00 to
    17:00), by RMSE in mg/L.  Returns (member, dose multiplier, rmse)."""
    C = np.exp(Z[:, hours, :].astype(float))                                  # members x hours x J
    T = truth_by_hour.loc[hours, list(junctions)].values.astype(float)
    A = (C ** 2).mean(axis=(1, 2))
    B = (C * T[None]).mean(axis=(1, 2))
    tt = float((T ** 2).mean())
    d = np.asarray(doses, float)
    mse = d[None, :] ** 2 * A[:, None] - 2 * d[None, :] * B[:, None] + tt
    k, j = np.unravel_index(int(np.argmin(mse)), mse.shape)
    return int(k), float(d[j]), float(np.sqrt(max(mse[k, j], 0.0)))


# ----------------------------------------------------------------------------- one (network, seed, variant)
_W: dict = {}      # per worker process: the banks and settings (set by _init_worker)


def _init_worker(net, cache_dir, temps_m, temps_pm2_path, workdir):
    """Pool initializer: a private working directory (nominal_scenario writes EPANET files to the cwd) and the banks read
    once from outputs/cache.  One BLAS thread per worker comes from the environment variables run() sets before the
    pool starts (this numpy uses Apple's Accelerate, which threadpoolctl cannot limit), so the numbers do not depend on
    threading; threadpool_limits is kept for builds that use OpenBLAS."""
    warnings.filterwarnings("ignore")
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    _W.update(net=net, cache_dir=cache_dir,
              M=covariate_bank(sc, temps_m, cache_dir=cache_dir, cache="read"),
              M1b=covariate_bank(sc, temps_m, wall_mode="mass_transfer_only", cache_dir=cache_dir, cache="read"))
    with open(temps_pm2_path, "rb") as fh:
        _W["pm2"] = pickle.load(fh)
    _W["Mx"] = _W["M"].merged(_W["pm2"])


def _dose(net: str) -> float:
    from .experiment import NET_TRUTH
    return float(NET_TRUTH.get(net, {}).get("source_dose", 1.2))


def _task(args):
    net, seed, variant = args
    t0 = time.time()
    out = run_task(net, seed, variant, _W["M"], _W["M1b"], _W["Mx"], _W["cache_dir"])
    out["seconds"] = time.time() - t0
    return out


def run_task(net: str, seed: int, variant: str, bank_m: Bank, bank_m1b: Bank | None, bank_x: Bank | None, cache_dir: str,
             months: int = MONTHS, targets=ROLLING_TARGETS, models=None, extrapolate: bool = True,
             temp_errors=None, truth_wall_law: str = "theta_w", bank_m2: "WallBank | None" = None) -> dict:
    """Every model and test for one network, seed and variant.  Returns the rows (one per test, temperature error,
    month and model), the truth's monthly low counts, and (seed 0, V1) what the summer-forecast figure needs.
    months, targets, models (a set of model names) and extrapolate=False run a subset, for the saved check that
    recomputes committed rows; the rows it does produce are the same numbers (every fit and draw is seeded).
    temp_errors (None = TEMP_ERRORS_C; () skips the plus or minus 2 C test, which needs bank_x) and truth_wall_law
    (simulate.build_scenario's; 'theta_w' is the pre-registered truth) are for mechanism_checks and task 10b.
    bank_m2 (task 10b): the bulk E/R x wall E/R bank (WallBank, merged with its plus or minus 2 C blocks when that test
    runs).  When it is given, model M2 (WallSeasonalSimGP24) is fitted beside M and M1b in every test, the plus or minus
    2 C test runs on M2 instead of M, and the oracle is the best member of M2's bank (which holds M's and M1b's).
    Without it every row is the task-10 row."""
    from .experiment import NET_TRUTH
    from .features import build_features
    from .pilot import synthetic_log
    from .simulate import nominal_scenario
    truth_kw = NET_TRUTH.get(net, {})
    sc = nominal_scenario(net, 14, _dose(net))
    X = build_features(sc)
    sched = plant_schedule(variant)
    errs = TEMP_ERRORS_C if temp_errors is None else tuple(temp_errors)
    wall_kw = {} if truth_wall_law == "theta_w" else {"truth_wall_law": truth_wall_law}
    log, _, plant, truths = synthetic_log(net, months=months, per_month=PER_MONTH, rotating=ROTATING, seed=seed,
                                          schedule=sched, return_truth=True, **truth_kw, **wall_kw)
    want = (lambda name: True) if models is None else (lambda name: name in models)
    log = log.assign(mi=[int(m[-2:]) for m in log.month])
    T_of = {m: sched[m - 1]["temp_C"] for m in range(1, months + 1)}
    tmin = {t["month"]: t["truth_daily_min"] for t in truths}
    tbh = {t["month"]: t["truth_by_hour"] for t in truths}
    chem = {t["month"]: t["chem"] for t in truths}
    junctions = list(sc.junctions)
    jidx = {j: i for i, j in enumerate(junctions)}
    kb_net = float(truth_kw.get("kb_per_day", 0.40))
    # M's bank: the M blocks plus the plus or minus 2 C blocks when given (identical blocks at every M temperature, so
    # M's numbers do not depend on it; the forecast-month error test needs the extra temperatures)
    bank_m = bank_x if bank_x is not None else bank_m
    # task 10b: with M2's bank, M2 carries the plus or minus 2 C test and the oracle is taken over its (larger) bank
    err_bank, err_name = (bank_m2, "M2") if bank_m2 is not None else (bank_x, "M")
    orc_bank = bank_m2 if bank_m2 is not None else bank_m
    rows = []

    def win(ms):
        return log[log.mi.isin(ms)]

    def blind(ms):
        return SimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=THRESHOLD).fit(win(ms)[["junction", "hour", "y"]])

    def banked(ms, bank, target_T, prior="uniform", shift=0.0):
        s = win(ms)[["junction", "hour", "y", "temp_C"]].copy()
        s["temp_C"] = s.temp_C + shift
        if isinstance(bank, WallBank):      # task 10b: M2, with its own nested prior
            return WallSeasonalSimGP24(sc, X, bank, seed=seed, cache_dir=cache_dir,
                                       threshold=THRESHOLD).fit(s, target_temp_C=target_T)
        return SeasonalSimGP24(sc, X, bank, prior=prior, seed=seed, cache_dir=cache_dir,
                               threshold=THRESHOLD).fit(s, target_temp_C=target_T)

    def predictions(model):
        z_mu, z_sd = model.predict_hours()
        return z_mu, z_sd, model.predict_daily_min()

    def base(test_name, err, m, model_name, window, n_train):
        c = chem[m]
        return {"network": net, "variant": variant, "test": test_name, "temp_error": err, "seed": seed, "month": m,
                "season": season_of(m), "plant_temp_C": T_of[m], "model": model_name, "window": window,
                "n_train": n_train, "E_true_K": c["E_true_K"], "theta_w": c["theta_w"],
                "truth_bulk_month_factor": c["bulk_month_factor"],
                "truth_kb_eff_per_day": kb_net * c["bulk_month_factor"] * c["bulk_temp_factor"],
                "truth_mixed_temp_C": c.get("warming", {}).get("mixed_temp_C", T_of[m]),
                "n_true_low_all": int((tmin[m] < THRESHOLD).sum())}

    def score_model(test_name, err, m, name, window, train, pred, model=None, unsampled=None):
        test = log[log.mi == m]
        seen = test.junction.isin(set(train.junction)).values
        z_mu, z_sd, pmin = pred
        ii, hh = np.array([jidx[j] for j in test.junction]), test.hour.values.astype(int)
        r = _reading_sums(z_mu[hh, ii], z_sd[hh, ii], test, seen, THRESHOLD)
        flag = pmin["p_below"] > 0.5
        d = _map_sums(pmin["median"], pmin["lo90"], pmin["hi90"], flag, tmin[m], unsampled, THRESHOLD)
        row = {**base(test_name, err, m, name, window, len(train)), **r, **d, "n_flag_all": int(flag.sum())}
        if model is not None:
            P = np.asarray(model.base_params if isinstance(model, SeasonalSimGP24) else model.params, float)
            row.update(map_kb=float(model.map_params_[0]), map_kw=float(model.map_params_[1]),
                       map_gamma=float(model.map_params_[2]), map_dose=float(model.map_dose_))
            if isinstance(model, SeasonalSimGP24):
                row.update(post_kb=model.kb20(), post_kw=model.rate20(1), map_hypothesis=model.map_hypothesis_,
                           E_post_mean_K=model.e_posterior_mean(),
                           **model.posterior_columns(), target_temp_C=model.target_temp_C)
            else:
                row.update(post_kb=float(model.w_ @ P[:, 0]), post_kw=float(model.w_ @ P[:, 1]))
        rows.append(row)

    def score_point(test_name, err, m, name, window, train, read_pred, map_pred=None, unsampled=None):
        """Point predictions: read_pred mg/L per test reading (persistence, network mean), map_pred mg/L per junction."""
        test = log[log.mi == m]
        seen = test.junction.isin(set(train.junction)).values
        r = _reading_sums(np.log(np.clip(read_pred, 1e-9, None)), None, test, seen, THRESHOLD)
        if map_pred is not None:
            flag = map_pred < THRESHOLD
            d = _map_sums(map_pred, None, None, flag, tmin[m], unsampled, THRESHOLD)
            nf = int(flag.sum())
        else:
            d, nf = _nan_map(), float("nan")
        rows.append({**base(test_name, err, m, name, window, len(train)), **r, **d, "n_flag_all": nf})

    def baselines(test_name, err, m, window, train, unsampled):
        test = log[log.mi == m]
        last = train.sort_values("date").groupby("tap_id").y.last()
        pers = test.tap_id.map(last).fillna(train.y.mean()).values.astype(float)
        score_point(test_name, err, m, "persistence", window, train, pers)
        score_point(test_name, err, m, "network_mean", window, train, np.full(len(test), float(train.y.mean())))

    def oracle(test_name, err, m, window, train, unsampled):
        Z = orc_bank.stack(T_of[m])
        k, dmult, rm = oracle_member(Z, tbh[m], junctions)
        zz = Z[k].astype(float) + math.log(dmult)                       # 24 x J
        test = log[log.mi == m]
        ii, hh = np.array([jidx[j] for j in test.junction]), test.hour.values.astype(int)
        read = np.exp(zz[hh, ii])
        dmin = pd.Series(np.exp(zz.min(axis=0)), index=junctions)
        score_point(test_name, err, m, "oracle", window, train, read, dmin, unsampled)
        rows[-1].update(oracle_member=k, oracle_hypothesis=orc_bank.hypotheses[k // orc_bank.n_base],
                        oracle_dose=dmult, oracle_daytime_rmse=rm, map_kb=float(orc_bank.stacked_params[k][0]),
                        map_kw=float(orc_bank.stacked_params[k][1]))

    def wlabel(ms):
        return f"{min(ms)}-{max(ms)}"

    fig = None
    for m in targets:
        w3, w1, w6 = (m - 3, m - 2, m - 1), (m - 1,), tuple(range(max(1, m - 6), m))
        train3 = win(w3)
        uns = [j for j in junctions if j not in set(train3.junction)]
        b0 = blind(w3); p_b0 = predictions(b0)
        score_model("rolling", "none", m, "B0", wlabel(w3), train3, p_b0, b0, uns)
        if want("B1"):
            b1 = blind(w1)
            score_model("rolling", "none", m, "B1", wlabel(w1), win(w1), predictions(b1), b1, uns)
        mm = banked(w3, bank_m, T_of[m])
        score_model("rolling", "none", m, "M", wlabel(w3), train3, predictions(mm), mm, uns)
        if want("M6"):
            m6 = banked(w6, bank_m, T_of[m])
            score_model("rolling", "none", m, "M6", wlabel(w6), win(w6), predictions(m6), m6, uns)
        if want("M1b"):
            m1b = banked(w3, bank_m1b, T_of[m])
            score_model("rolling", "none", m, "M1b", wlabel(w3), train3, predictions(m1b), m1b, uns)
        m2 = None
        if bank_m2 is not None and want("M2"):
            m2 = banked(w3, bank_m2, T_of[m])
            score_model("rolling", "none", m, "M2", wlabel(w3), train3, predictions(m2), m2, uns)
        if variant == "V1":
            for pr in PRIORS[1:]:
                if want(f"M_prior_{pr}"):
                    mp = banked(w3, bank_m, T_of[m], prior=pr)
                    score_model("rolling", "none", m, f"M_prior_{pr}", wlabel(w3), train3, predictions(mp), mp, uns)
        if want("persistence"):
            baselines("rolling", "none", m, wlabel(w3), train3, uns)
        if want("oracle"):
            oracle("rolling", "none", m, wlabel(w3), train3, uns)
        # extrapolation without refit: the rolling fits on January to March and July to September predict later months
        for test_name, (fit_ms, ex_targets) in (EXTRAPOLATION.items() if extrapolate else ()):
            if w3 != fit_ms:
                continue
            # the plus or minus 2 C error in every logged month needs its own fit (once per error; set_target per month)
            shifted = ({dT: banked(w3, err_bank, T_of[ex_targets[0]] + dT, shift=dT) for dT in errs}
                       if variant == "V1" else {})
            err_model = m2 if bank_m2 is not None else mm
            for t in ex_targets:
                score_model(test_name, "none", t, "B0", wlabel(w3), train3, p_b0, b0, uns)
                for name, mod in (("M", mm), ("M1b", m1b)) + ((("M2", m2),) if m2 is not None else ()):
                    mod.set_target(T_of[t])
                    pred = predictions(mod)
                    score_model(test_name, "none", t, name, wlabel(w3), train3, pred, mod, uns)
                    if name == "M" and variant == "V1" and seed == 0 and t == 7:
                        fig = {"truth_min": tmin[7], "B0": p_b0[2], "M": pred[2], "plant_temp_C": T_of[7],
                               "fit_temps_C": [T_of[x] for x in w3], "E_true_K": chem[7]["E_true_K"]}
                baselines(test_name, "none", t, wlabel(w3), train3, uns)
                oracle(test_name, "none", t, wlabel(w3), train3, uns)
                if variant == "V1":
                    for dT in errs:
                        err_model.set_target(T_of[t] + dT)
                        score_model(test_name, f"{dT:+g}C_forecast_month", t, err_name, wlabel(w3), train3,
                                    predictions(err_model), err_model, uns)
                        ma = shifted[dT].set_target(T_of[t] + dT)
                        score_model(test_name, f"{dT:+g}C_every_month", t, err_name, wlabel(w3), train3, predictions(ma), ma, uns)
    monthly = [{"variant": variant, "seed": seed, "month": t["month"], "plant_temp_C": T_of[t["month"]],
                "n_true_low_all": int((t["truth_daily_min"] < THRESHOLD).sum()), "n_junctions": len(junctions),
                "E_true_K": t["chem"]["E_true_K"], "theta_w": t["chem"]["theta_w"],
                "bulk_month_factor": t["chem"]["bulk_month_factor"], "bulk_temp_factor": t["chem"]["bulk_temp_factor"],
                "mixed_temp_C": t["chem"].get("warming", {}).get("mixed_temp_C", T_of[t["month"]])} for t in truths]
    return {"net": net, "seed": seed, "variant": variant, "rows": rows, "monthly": monthly, "fig": fig,
            "plant": plant.to_dict("records")}


# ----------------------------------------------------------------------------- pooling and acceptance
def pooled(df: pd.DataFrame) -> dict:
    """Pooled metrics over a set of rows (sums over seeds and months)."""
    s = df.sum(numeric_only=True)
    out = {"n_rows": int(len(df))}

    def ratio(a, b):
        return float(s[a] / s[b]) if s[b] > 0 else float("nan")
    for tag in ("all", "seen", "new"):
        n = s[f"r_{tag}_n"]
        d = {"n": int(n), "rmse": float(np.sqrt(s[f"r_{tag}_sse"] / n)) if n else float("nan"),
             "mae": float(s[f"r_{tag}_sae"] / n) if n else float("nan"),
             "n_low": int(s[f"r_{tag}_low"]), "recall": ratio(f"r_{tag}_tp", f"r_{tag}_low"),
             "false_alarms": int(s[f"r_{tag}_fp"])}
        for q in (LEVELS if tag == "all" else (90,)):
            col = df[f"r_{tag}_in{q}"]
            d[f"coverage{q}"] = float(col.sum() / n) if n and col.notna().all() else float("nan")
        out[f"readings_{tag}"] = d
    if df["d_n"].notna().all() and s["d_n"] > 0:
        n = s["d_n"]
        tp, fp = s["d_tp"], s["d_fp"]
        out["daily_min_map"] = {"n": int(n), "rmse": float(np.sqrt(s["d_sse"] / n)), "mae": float(s["d_sae"] / n),
                                "coverage90": float(df.d_in90.sum() / n) if df.d_in90.notna().all() else float("nan"),
                                "n_low": int(s["d_low"]), "recall": ratio("d_tp", "d_low"), "false_alarms": int(fp),
                                "precision": float(tp / (tp + fp)) if tp + fp > 0 else float("nan")}
    return out


def _class_rows(df, cls):
    return df if cls == "all_months" else df[df.season == cls]


def bootstrap_ratio(df: pd.DataFrame, num: str, den: str, prefix: str, rng_seed: int = BOOT_SEED, n: int = BOOT_N) -> dict:
    """Paired bootstrap over seeds of pooled RMSE(num) / RMSE(den); prefix 'r_all' (readings) or 'd' (daily-min map)."""
    a = df[df.model == num].groupby("seed")[[f"{prefix}_sse", f"{prefix}_n"]].sum()
    b = df[df.model == den].groupby("seed")[[f"{prefix}_sse", f"{prefix}_n"]].sum()
    seeds = sorted(set(a.index) & set(b.index))
    a, b = a.loc[seeds].values, b.loc[seeds].values
    point = math.sqrt(a[:, 0].sum() / a[:, 1].sum()) / math.sqrt(b[:, 0].sum() / b[:, 1].sum())
    idx = np.random.default_rng(rng_seed).integers(0, len(seeds), size=(n, len(seeds)))
    r = np.sqrt(a[idx, 0].sum(1) / a[idx, 1].sum(1)) / np.sqrt(b[idx, 0].sum(1) / b[idx, 1].sum(1))
    lo, hi = np.quantile(r, [0.05, 0.95])
    return {"ratio": point, "lo90": float(lo), "hi90": float(hi), "n_seeds": len(seeds), "n_resamples": n}


def _get(P, *keys):
    for k in keys:
        P = P[k]
    return P


def acceptance(net: str, P: dict, boot: dict, kbsd: dict, checks: dict | None) -> dict:
    """The pre-registered bars (journal, Task 10, written before the first scored run)."""
    A = {}
    A["A1"] = {"bar": "outputs/chem/checks_report.json: every check passes, the task-10 checks included",
               "pass": None, "judged_from": "outputs/chem/checks_report.json after this run (its season_outputs_reproduce "
                                            "check reads this file, so the verdict cannot be written here)"}
    # A2: V2 control, rolling
    c = {}
    for kind, key in (("readings", "readings_all"), ("daily_min_map", "daily_min_map")):
        m, b = _get(P, "V2", "rolling", "all_months", "M", key), _get(P, "V2", "rolling", "all_months", "B0", key)
        c[kind] = {"rmse_ratio_M_B0": m["rmse"] / b["rmse"], "recall_M": m["recall"], "recall_B0": b["recall"],
                   "recall_diff": m["recall"] - b["recall"]}
    ok2 = all(abs(v["rmse_ratio_M_B0"] - 1) <= 0.05 and abs(v["recall_diff"]) <= 0.03 for v in c.values())
    A["A2"] = {"bar": "V2 control, rolling: |RMSE_M/RMSE_B0 - 1| <= 0.05 and |recall_M - recall_B0| <= 0.03, "
                      "for held-out readings and for the daily-minimum map", **c, "pass": bool(ok2)}
    # A3: V1 rolling, readings
    m, b = _get(P, "V1", "rolling", "all_months", "M", "readings_all"), _get(P, "V1", "rolling", "all_months", "B0", "readings_all")
    bt = boot["V1"]["rolling"]["all_months"]["readings"]
    A["A3"] = {"bar": "V1 rolling: pooled held-out-reading RMSE_M <= RMSE_B0; a gain is claimed only if the bootstrap 90% "
                      "upper bound of RMSE_M/RMSE_B0 is below 1.0",
               "rmse_M": m["rmse"], "rmse_B0": b["rmse"], "bootstrap": bt, "pass": bool(m["rmse"] <= b["rmse"]),
               "gain_claimed": bool(bt["hi90"] < 1.0),
               "readme_wording": "measurable gain" if bt["hi90"] < 1.0 else f"no measurable gain on {net}"}
    # A4: V1 extrapolation, Net3 bar (computed for every network, the bar applies to Net3)
    e = _get(P, "V1", "jan_mar_to_jul_sep", "none", "all_months")
    mM, mB = e["M"]["daily_min_map"], e["B0"]["daily_min_map"]
    r = _get(P, "V1", "jul_sep_to_dec", "none", "all_months")
    fa_M, fa_B = r["M"]["daily_min_map"]["false_alarms"], r["B0"]["daily_min_map"]["false_alarms"]
    parts = {"rmse": mM["rmse"] <= 0.85 * mB["rmse"], "recall": mM["recall"] >= mB["recall"] + 0.10,
             "coverage90": mM["coverage90"] >= 0.80, "reverse_false_alarms": fa_M <= fa_B}
    A["A4"] = {"bar": "V1, fit January to March, predict July to September, daily-minimum map: RMSE_M <= 0.85 RMSE_B0, "
                      "recall_M >= recall_B0 + 0.10, coverage90_M >= 0.80; fit July to September, predict December: "
                      "false alarms_M <= false alarms_B0 (the bar is set for Net3)",
               "applies": net == "Net3", "rmse_M": mM["rmse"], "rmse_B0": mB["rmse"], "rmse_ratio": mM["rmse"] / mB["rmse"],
               "recall_M": mM["recall"], "recall_B0": mB["recall"], "coverage90_M": mM["coverage90"],
               "coverage90_B0": mB["coverage90"], "false_alarms_M_jul_sep": mM["false_alarms"],
               "false_alarms_B0_jul_sep": mB["false_alarms"], "reverse_false_alarms_M": fa_M,
               "reverse_false_alarms_B0": fa_B, "parts": parts, "pass": bool(all(parts.values()))}
    # A5
    A["A5"] = {"bar": "V1 rolling: median over seeds of the SD of M's kb20 across the 9 windows < the same for B0's MAP kb",
               **kbsd, "pass": bool(kbsd["median_sd_M_kb20"] < kbsd["median_sd_B0_map_kb"])}
    # A6
    base_rec = e["M"]["daily_min_map"]["recall"]
    errs = {}
    for dT in TEMP_ERRORS_C:
        for mode in ("forecast_month", "every_month"):
            k = f"{dT:+g}C_{mode}"
            rec = _get(P, "V1", "jan_mar_to_jul_sep", k, "all_months", "M", "daily_min_map")["recall"]
            errs[k] = {"recall": rec, "cost": base_rec - rec}
    A["A6"] = {"bar": "in the A4 test, a logged temperature off by +2 C or -2 C (forecast month only, or every logged month) "
                      "costs at most 0.10 of M's daily-minimum recall", "recall_no_error": base_rec, "errors": errs,
               "worst_cost": max(v["cost"] for v in errs.values()),
               "pass": bool(all(v["cost"] <= 0.10 for v in errs.values()))}
    # A7: stop rule, rolling, every variant and season class
    trig, fa = [], []
    table = {}
    for v in [x for x in VARIANTS if x in P]:
        for cls in ("warming", "cooling"):
            g = _get(P, v, "rolling", cls)
            ent = {}
            for kind, key in (("readings", "readings_all"), ("daily_min_map", "daily_min_map")):
                m, b = g["M"][key], g["B0"][key]
                rr = m["rmse"] / b["rmse"]
                rc = m["recall"] / b["recall"] if b["recall"] and not math.isnan(b["recall"]) else float("nan")
                far = m["false_alarms"] / b["false_alarms"] if b["false_alarms"] else float("nan")
                ent[kind] = {"rmse_M": m["rmse"], "rmse_B0": b["rmse"], "rmse_ratio": rr, "recall_M": m["recall"],
                             "recall_B0": b["recall"], "recall_ratio": rc, "false_alarms_M": m["false_alarms"],
                             "false_alarms_B0": b["false_alarms"], "false_alarm_ratio": far}
                if rr > 1.10:
                    trig.append(f"{v} {cls} {kind}: RMSE ratio {rr:.3f}")
                if not math.isnan(rc) and rc < 0.90:
                    trig.append(f"{v} {cls} {kind}: recall ratio {rc:.3f}")
                if not math.isnan(far) and far > 1.10:
                    fa.append(f"{v} {cls} {kind}: false alarms {m['false_alarms']} vs {b['false_alarms']} (x{far:.2f})")
            table[f"{v}_{cls}"] = ent
    A["A7"] = {"bar": "stop rule, rolling test, each variant and season class: stop if RMSE_M > 1.10 RMSE_B0 or "
                      "recall_M < 0.90 recall_B0, for held-out readings or the daily-minimum map; false alarms up by more "
                      "than 10% are listed as regressions", "by_variant_and_season": table, "triggered": trig,
               "false_alarm_rises_over_10pct": fa, "stop": bool(trig)}
    return A


def kb_sd_summary(df: pd.DataFrame) -> dict:
    """SD across the 9 rolling windows, per seed, of M's kb20 and B0's MAP kb (A5), plus like-for-like extras."""
    r = df[(df.variant == "V1") & (df.test == "rolling")]
    per = {}
    for name, model, col in (("M_kb20", "M", "post_kb"), ("B0_map_kb", "B0", "map_kb"), ("B0_post_kb", "B0", "post_kb"),
                             ("M_map_kb", "M", "map_kb"), ("M_kw20", "M", "post_kw"), ("B0_post_kw", "B0", "post_kw")):
        per[name] = r[r.model == model].groupby("seed")[col].std(ddof=1)
    return {"median_sd_M_kb20": float(per["M_kb20"].median()), "median_sd_B0_map_kb": float(per["B0_map_kb"].median()),
            "like_for_like_median_sd_B0_posterior_mean_kb": float(per["B0_post_kb"].median()),
            "median_sd_M_map_kb": float(per["M_map_kb"].median()),
            "median_sd_M_kw20": float(per["M_kw20"].median()), "median_sd_B0_posterior_mean_kw": float(per["B0_post_kw"].median()),
            "sd_by_seed": {k: [round(float(x), 4) for x in v.values] for k, v in per.items()}}


def summarise(net: str, df: pd.DataFrame, monthly: pd.DataFrame, seeds, checks: dict | None, settings: dict) -> dict:
    P, boot = {}, {}
    for v in sorted(df.variant.unique()):
        dv = df[df.variant == v]
        for test in sorted(dv.test.unique()):
            dt = dv[dv.test == test]
            for err in sorted(dt.temp_error.unique()):
                de = dt[dt.temp_error == err]
                key = test if err == "none" and test == "rolling" else None
                for cls in ("all_months", "warming", "cooling"):
                    dc = _class_rows(de, cls)
                    if not len(dc):
                        continue
                    node = P.setdefault(v, {}).setdefault(test, {})
                    node = node if test == "rolling" else node.setdefault(err, {})
                    node = node.setdefault(cls, {})
                    for model in sorted(dc.model.unique()):
                        node[model] = pooled(dc[dc.model == model])
                    if err == "none" and {"M", "B0"} <= set(dc.model):
                        b = boot.setdefault(v, {}).setdefault(test, {}).setdefault(cls, {})
                        b["readings"] = bootstrap_ratio(dc, "M", "B0", "r_all")
                        b["daily_min_map"] = bootstrap_ratio(dc, "M", "B0", "d")
                        if "M1b" in set(dc.model):
                            b["readings_M1b"] = bootstrap_ratio(dc, "M1b", "B0", "r_all")
                            b["daily_min_map_M1b"] = bootstrap_ratio(dc, "M1b", "B0", "d")
    kbsd = kb_sd_summary(df)
    acc = acceptance(net, P, boot, kbsd, checks)
    # pilot coverage bands for M's readings (V1 rolling)
    g = P["V1"]["rolling"]["all_months"]["M"]
    bands = {"seen_taps": {"coverage90": g["readings_seen"]["coverage90"], "band": [0.85, 0.95]},
             "new_taps": {"coverage90": g["readings_new"]["coverage90"], "band": [0.80, 0.95]}}
    for v in bands.values():
        v["in_band"] = bool(v["band"][0] <= v["coverage90"] <= v["band"][1])
    # hypothesis posterior against the truth's E (never used as evidence)
    r = df[(df.test == "rolling") & (df.model == "M")]
    hyp = {}
    for v in sorted(r.variant.unique()):
        rv = r[r.variant == v]
        terc = pd.qcut(rv.E_true_K.rank(method="first"), 3, labels=["low", "middle", "high"])
        hyp[v] = {"mean_posterior": {h: float(rv[f"P_{h}"].mean()) for h in ("H0", "E5000", "E8000", "E12000")},
                  "by_E_true_tercile": {str(t): {"E_true_K_range": [float(g_.E_true_K.min()), float(g_.E_true_K.max())],
                                                 **{h: float(g_[f"P_{h}"].mean()) for h in ("H0", "E5000", "E8000", "E12000")}}
                                        for t, g_ in rv.groupby(terc, observed=True)},
                  "corr_E_post_mean_vs_E_true": float(np.corrcoef(rv.E_post_mean_K, rv.E_true_K)[0, 1])}
    truth_low = (monthly.groupby(["variant", "month"]).n_true_low_all.mean() / monthly.n_junctions.iloc[0]).round(4)
    return {"generated_by": f"python -m residualmap.seasonal {net} {len(seeds)}", "network": net,
            "seeds": [int(s) for s in seeds],
            "about": "Simulation only. 12-month synthetic grab logs against hidden truths with their own temperature "
                     "sensitivity (V1), a constant-temperature control (V2) and in-network warming the model cannot see "
                     "(V3); today's temperature-blind model (B0) against the temperature-aware bank (M) and the other "
                     "models. Pooled over seeds and target months. The temperature schedules are ASSUMPTIONS.",
            "settings": settings, "acceptance": acc, "pilot_coverage_bands": bands, "pooled": P, "bootstrap": boot,
            "hypothesis_posterior_vs_truth": hyp,
            "truth_share_below_threshold_daily_min_by_month": {v: {int(m): float(x) for (vv, m), x in truth_low.items() if vv == v}
                                                              for v in sorted(monthly.variant.unique())}}


# ----------------------------------------------------------------------------- figures
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_season(net, df, monthly, summ, out):
    plt = _plt()
    from .experiment import BIG
    r = df[(df.variant == "V1") & (df.test == "rolling")]
    mo = monthly[monthly.variant == "V1"]
    J = int(mo.n_junctions.iloc[0])
    months = list(range(1, 13))
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(1, 3, figsize=(24, 6.4))
        a = ax[0]
        a.plot(months, [plant_schedule("V1")[m - 1]["temp_C"] for m in months], "o-", color="tab:red", label="plant water (logged)")
        a.plot(months, [plant_schedule("V3")[m - 1]["soil_temp_C"] for m in months], "s--", color="tab:brown",
               label="soil (V3 truth only)")
        a.set(xticks=months, xlabel="month", ylabel="temperature (C)", title="Assumed coastal schedule (ASSUMPTION)")
        a.grid(alpha=0.3); a.legend()
        a = ax[1]
        tl = mo.groupby("month").n_true_low_all.mean() / J
        a.plot(tl.index, tl.values, "k-o", lw=2.5, label="truth (mean over seeds)")
        for model, c, lab in (("B0", "tab:gray", "B0, temperature-blind (today)"), ("M", "tab:blue", "M, temperature bank")):
            f = r[r.model == model].groupby("month").n_flag_all.mean() / J
            a.plot(f.index, f.values, "o--", color=c, label=f"{lab}: flagged")
        a.set(xticks=months, xlabel="month", ylabel="share of junctions", ylim=(0, 1),
              title=f"Junctions below {THRESHOLD} mg/L at the daily minimum\n(rolling: fit the previous 3 months)")
        a.grid(alpha=0.3); a.legend(fontsize=10)
        a = ax[2]
        for model, c in (("B0", "tab:gray"), ("B1", "tab:orange"), ("M", "tab:blue"), ("M1b", "tab:green"), ("oracle", "k")):
            g = r[r.model == model].groupby("month")[["r_all_sse", "r_all_n"]].sum()
            a.plot(g.index, np.sqrt(g.r_all_sse / g.r_all_n), "o-", color=c, label=model)
        a.set(xticks=months, xlabel="held-out month", ylabel="RMSE of held-out readings (mg/L)",
              title="Held-out readings, pooled over seeds")
        a.grid(alpha=0.3); a.legend(fontsize=10)
        bt = summ["bootstrap"]["V1"]["rolling"]["all_months"]["readings"]
        fig.suptitle(f"{net}, simulated 12-month logs, {len(summ['seeds'])} seeds: held-out reading RMSE ratio M / B0 "
                     f"{bt['ratio']:.3f} (90% interval {bt['lo90']:.3f} to {bt['hi90']:.3f})", fontsize=15)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_kb_by_month(net, df, monthly, out):
    plt = _plt()
    from .experiment import BIG
    r = df[(df.variant == "V1") & (df.test == "rolling")]
    mo = monthly[monthly.variant == "V1"]
    kbn = _kb_net(net)
    xs = list(ROLLING_TARGETS)
    with plt.rc_context(BIG):
        fig, a = plt.subplots(figsize=(13, 6.6))
        for model, col, c, lab in (("B0", "map_kb", "tab:gray", "B0 (today): most likely kb, refitted each window"),
                                   ("B0", "post_kb", "tab:orange", "B0: posterior-mean kb"),
                                   ("M", "post_kb", "tab:blue", "M: kb20, temperature-normalised (posterior mean)")):
            g = r[r.model == model].groupby("month")[col]
            med, lo, hi = g.median(), g.quantile(0.25), g.quantile(0.75)
            a.plot(med.index, med.values, "o-", color=c, lw=2.2, label=lab)
            a.fill_between(med.index, lo.values, hi.values, color=c, alpha=0.15)
        # the truth's bulk rate over the window's 3 months: at the water's temperature, and normalised to 20 C
        eff, t20 = [], []
        for m in xs:
            w = mo[mo.month.isin([m - 3, m - 2, m - 1])]
            eff.append(kbn * float((w.bulk_month_factor * w.bulk_temp_factor).groupby(w.seed).mean().median()))
            t20.append(kbn * float(w.bulk_month_factor.groupby(w.seed).mean().median()))
        a.plot(xs, eff, "k--", lw=1.6, label="truth: bulk kb at the water's temperature (window mean, median over seeds)")
        a.plot(xs, t20, "k:", lw=1.6, label="truth: bulk kb at 20 C (window mean, median over seeds)")
        a.set(xticks=xs, xlabel="held-out month (window = the 3 months before it)", ylabel="bulk decay kb (1/day)",
              title=f"{net}: the blind refit's kb and the temperature-normalised kb20\n(median and quartiles over seeds; simulated)")
        a2 = a.twinx()
        a2.plot(xs, [np.mean([plant_schedule("V1")[x - 1]["temp_C"] for x in (m - 3, m - 2, m - 1)]) for m in xs],
                "-", color="tab:red", alpha=0.35, lw=6)
        a2.set_ylabel("window mean plant temperature (C, red band)")
        a.grid(alpha=0.3); a.legend(fontsize=10, loc="upper left")
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def _kb_net(net):
    from .experiment import NET_TRUTH
    return float(NET_TRUTH.get(net, {}).get("kb_per_day", 0.40))


def plot_summer_forecast(net, fig_in, out, summ):
    plt = _plt()
    import wntr
    from .experiment import BIG
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    tm = fig_in["truth_min"]
    J = len(sc.junctions)
    size = max(14, int(4000 / J))
    vmax = float(max(tm.max(), fig_in["B0"]["median"].max(), fig_in["M"]["median"].max()))
    true_low = tm < THRESHOLD
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(1, 3, figsize=(24, 7))
        panels = [("truth", tm, true_low, f"TRUE July daily minimum (simulated)\n{int(true_low.sum())} of {J} junctions below {THRESHOLD} mg/L"),
                  ("B0", fig_in["B0"]["median"], fig_in["B0"]["p_below"] > 0.5, None),
                  ("M", fig_in["M"]["median"], fig_in["M"]["p_below"] > 0.5, None)]
        for a, (name, val, ring, title) in zip(ax, panels):
            if title is None:
                f = ring
                tp, fp = int((f & true_low).sum()), int((f & ~true_low).sum())
                lab = ("today's model, temperature-blind" if name == "B0" else
                       f"temperature bank, told July is {fig_in['plant_temp_C']:g} C")
                title = (f"{name}: {lab}\nfitted on January to March only; flags {tp} of {int(true_low.sum())}, "
                         f"{fp} false alarms")
            wntr.graphics.plot_network(sc.wn, node_attribute=val.to_dict(), node_size=size, node_cmap="RdYlBu",
                                       node_range=(0, vmax), ax=a, link_width=0.6, add_colorbar=True, title=title)
            xy = np.array([sc.wn.get_node(j).coordinates for j in val.index[ring.values]]) if ring.any() else np.zeros((0, 2))
            if len(xy):
                a.scatter(xy[:, 0], xy[:, 1], s=size * 2.6, facecolors="none", edgecolors="k", linewidths=1.1, zorder=5)
        A4 = summ["acceptance"]["A4"]
        fig.suptitle(f"{net}, scenario 0: a July map predicted from January-to-March grab samples (logged plant water "
                     f"{', '.join(f'{t:g}' for t in fig_in['fit_temps_C'])} C). Rings: below {THRESHOLD} mg/L (truth) or flagged.\n"
                     f"Over all seeds, pooled over July, August and September: map RMSE M {A4['rmse_M']:.3f} vs B0 "
                     f"{A4['rmse_B0']:.3f} mg/L, recall {A4['recall_M']:.2f} vs {A4['recall_B0']:.2f}", fontsize=14)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------- the experiment
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


def all_bank_temps(variants=VARIANTS) -> tuple[list, list]:
    """(temperatures of the M and M1b banks, temperatures needed only by the plus or minus 2 C test)."""
    temps = sorted({_tkey(c["temp_C"]) for v in variants for c in plant_schedule(v)})
    need = set()
    if "V1" in variants:
        s = plant_schedule("V1")
        for fit_ms, targets in EXTRAPOLATION.values():
            for dT in TEMP_ERRORS_C:
                need |= {_tkey(s[m - 1]["temp_C"] + dT) for m in fit_ms + targets}
    return temps, sorted(need - set(temps))


def run(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache", workers: int | None = None,
        variants=VARIANTS) -> dict:
    os.makedirs(outdir, exist_ok=True)
    cache_dir = os.path.abspath(cache_dir)
    t0 = time.time()
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free, below {MIN_FREE_GB_RUN} GB")
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    temps_m, temps_pm2 = all_bank_temps(variants)
    tb = time.time()
    covariate_bank(sc, temps_m, cache_dir=cache_dir)
    covariate_bank(sc, temps_m, wall_mode="mass_transfer_only", cache_dir=cache_dir)
    pm2 = covariate_bank(sc, temps_pm2, cache_dir=cache_dir, cache="off")
    t_bank = time.time() - tb
    print(f"{net}: banks ready in {t_bank:.0f} s ({len(temps_m)} temperatures x {len(ER_HYPOTHESES_K)} E/R, twice; "
          f"{len(temps_pm2)} more for the +/-2 C test, not cached); {disk_free_gb(cache_dir):.2f} GB free", flush=True)
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free after the banks, below {MIN_FREE_GB_RUN} GB")
    workdir = tempfile.mkdtemp(prefix="rm_season_")
    pm2_path = os.path.join(cache_dir, f"_tmp_season_pm2_{net}_{os.getpid()}.pkl")
    results = []
    try:
        with open(pm2_path, "wb") as fh:
            pickle.dump(pm2, fh)
        del pm2
        tasks = [(net, s, v) for s in seeds for v in variants]
        n_workers = workers or max(1, min(len(tasks), (os.cpu_count() or 2) - 1))
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            os.environ[k] = "1"
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_worker,
                                 initargs=(net, cache_dir, temps_m, pm2_path, workdir)) as ex:
            for res in ex.map(_task, tasks):
                results.append(res)
                print(f"  {net} seed {res['seed']} {res['variant']}: {len(res['rows'])} rows in {res['seconds']:.0f} s", flush=True)
    finally:
        if os.path.exists(pm2_path):
            os.remove(pm2_path)
        shutil.rmtree(workdir, ignore_errors=True)
    df = pd.DataFrame([r for res in results for r in res["rows"]])
    order = {"rolling": 0, "jan_mar_to_jul_sep": 1, "jul_sep_to_dec": 2}
    df = df.assign(_o=df.test.map(order)).sort_values(["variant", "_o", "temp_error", "seed", "month", "model"],
                                                     kind="mergesort").drop(columns="_o").reset_index(drop=True)
    monthly = pd.DataFrame([r for res in results for r in res["monthly"]]).sort_values(["variant", "seed", "month"]).reset_index(drop=True)
    checks = None
    settings = {"months": MONTHS, "route_taps": PER_MONTH, "rotating_taps": ROTATING, "noise_sd_mgL": NOISE_SD,
                "sampling_hours": [min(DAY_HOURS), max(DAY_HOURS)], "threshold_mgL": THRESHOLD,
                "month_seed": "1000 + 100 seed + month index (pilot.synthetic_log with a schedule)",
                "schedules_ASSUMPTION": {v: plant_schedule(v) for v in variants}, "variants": VARIANT_LABELS,
                "hypotheses": ["H0"] + [f"E{e:g}" for e in ER_HYPOTHESES_K],
                "priors": {p: [float(x) for x in hypothesis_prior(p)] for p in PRIORS},
                "bank_temperatures_C": temps_m, "pm2_temperatures_C_not_cached": temps_pm2,
                "season_classes": {"warming": list(WARMING_MONTHS), "cooling": list(COOLING_MONTHS)},
                "rolling_targets": list(ROLLING_TARGETS), "extrapolation": {k: [list(a), list(b)] for k, (a, b) in EXTRAPOLATION.items()},
                "map_junctions": "every junction not sampled in the 3-month window, the same set for every model",
                "bootstrap": {"resamples": BOOT_N, "rng_seed": BOOT_SEED, "unit": "seed"},
                "truth": _truth_settings(net)}
    summ = summarise(net, df, monthly, seeds, checks, settings)
    df.to_csv(os.path.join(outdir, f"season_{net}.csv"), index=False)
    with open(os.path.join(outdir, f"summary_season_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    plot_season(net, df, monthly, summ, os.path.join(outdir, f"season_{net}.png"))
    plot_kb_by_month(net, df, monthly, os.path.join(outdir, f"kb_by_month_{net}.png"))
    fig_in = next((r["fig"] for r in results if r["fig"] is not None), None)
    if fig_in is not None:
        plot_summer_forecast(net, fig_in, os.path.join(outdir, f"summer_forecast_{net}.png"), summ)
    print(f"{net}: {time.time() - t0:.0f} s in all (banks {t_bank:.0f} s); A7 stop: {summ['acceptance']['A7']['stop']} "
          f"{summ['acceptance']['A7']['triggered']}", flush=True)
    return summ


def _truth_settings(net):
    from .experiment_chem import truth_settings
    from .chemistry import ER_TRUTH_RANGE_K, THETA_W_TRUTH_RANGE, WARMING_TAU_PIPE_H, WARMING_TAU_TANK_H
    return {**truth_settings(net), "E_true_K_range": list(ER_TRUTH_RANGE_K), "theta_w_range_ASSUMPTION": list(THETA_W_TRUTH_RANGE),
            "warming_tau_h_ASSUMPTION": {"pipe": WARMING_TAU_PIPE_H, "tank": WARMING_TAU_TANK_H},
            "chemistry_draws": "default_rng(20000 + seed): E_true, theta_w (simulate.hidden_chem_draws)"}


# ----------------------------------------------------------------------------- after the stop: what caused it
# Not pre-registered and never used to choose a model.  Run after the A7 stop, at the reviewers' request, to test the
# explanation the journal gives (task 10, finding 2) directly and to put counts and an interval on every A7 trigger.
# python -m residualmap.seasonal Net3 8 --mechanism -> outputs/chem/task10_mechanism_checks.json (no timings in it).
WALL_WORLDS = {"theta_w": "the pre-registered V1 truth: wall factor theta_w^(T - 20), theta_w ~ U(1.00, 1.07) per seed",
               "arrhenius": "the truth's wall factor is f(T; E_true), its own bulk factor: the structure M assumes",
               "none": "the truth's wall has no temperature factor (1 at every T): the structure M1b assumes"}
MECH_MODELS = ("B0", "M", "M1b")
MECH_COMPARE_COLS = ("r_all_sse", "r_all_tp", "r_all_fp", "r_all_in90", "d_sse", "d_tp", "d_fp", "d_in90", "post_kb",
                     "P_H0", "P_E12000")


@dataclass(frozen=True)
class OracleWallCondition:
    """A bank condition for the oracle-wall check: bulk rate times f(T; E) as in M, wall rate times the truth's own
    theta_w^(T - 20) instead of f(T; E), mass transfer at T (EPANET viscosity and diffusivity), as in M."""
    temp_C: float
    er_K: float
    theta_w: float

    def sim_kwargs(self) -> dict:
        return {"kb_scale": arrhenius(self.temp_C, self.er_K), "kw_scale": self.theta_w ** (self.temp_C - 20.0),
                "temp_C": self.temp_C}


def bootstrap_recall_ratio(df: pd.DataFrame, num: str, den: str, prefix: str, rng_seed: int = BOOT_SEED,
                           n: int = BOOT_N) -> dict:
    """Paired bootstrap over seeds of pooled recall(num) / recall(den); prefix 'r_all' (readings) or 'd' (daily-min
    map).  The same resampled seeds as bootstrap_ratio.  Resamples where a recall is undefined or den's is zero are
    dropped and counted."""
    cols = [f"{prefix}_tp", f"{prefix}_low"]
    a = df[df.model == num].groupby("seed")[cols].sum()
    b = df[df.model == den].groupby("seed")[cols].sum()
    seeds = sorted(set(a.index) & set(b.index))
    a, b = a.loc[seeds].values.astype(float), b.loc[seeds].values.astype(float)
    point = (a[:, 0].sum() / a[:, 1].sum()) / (b[:, 0].sum() / b[:, 1].sum())
    idx = np.random.default_rng(rng_seed).integers(0, len(seeds), size=(n, len(seeds)))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (a[idx, 0].sum(1) / a[idx, 1].sum(1)) / (b[idx, 0].sum(1) / b[idx, 1].sum(1))
    ok = np.isfinite(r)
    lo, hi = np.quantile(r[ok], [0.05, 0.95])
    return {"ratio": float(point), "lo90": float(lo), "hi90": float(hi), "n_seeds": len(seeds), "n_resamples": n,
            "n_resamples_dropped": int((~ok).sum())}


def a7_trigger_intervals(net: str, outdir: str = "outputs/chem") -> list[dict]:
    """Every A7 trigger of summary_season_<net>.json with its counts and a paired bootstrap 90% interval over seeds,
    from season_<net>.csv (the RMSE intervals equal the summary's own).  interval_excludes_1: the interval lies wholly
    on M's worse side of 1.0; interval_past_bar: wholly past the stop rule's 1.10 (RMSE) or 0.90 (recall)."""
    with open(os.path.join(outdir, f"summary_season_{net}.json")) as fh:
        summ = json.load(fh)
    df = pd.read_csv(os.path.join(outdir, f"season_{net}.csv"), float_precision="round_trip")
    roll = df[(df.test == "rolling") & (df.temp_error == "none")]
    out = []
    for key, ent in summ["acceptance"]["A7"]["by_variant_and_season"].items():
        v, cls = key.split("_")
        rows = roll[(roll.variant == v) & (roll.season == cls)]
        for kind, prefix in (("readings", "r_all"), ("daily_min_map", "d")):
            e = ent[kind]
            for metric in ("rmse", "recall"):
                r = e[f"{metric}_ratio"]
                if r is None or not ((metric == "rmse" and r > 1.10) or (metric == "recall" and r < 0.90)):
                    continue
                if metric == "rmse":
                    bt = bootstrap_ratio(rows, "M", "B0", prefix)
                    counts = {"n_scored": int(rows[rows.model == "M"][f"{prefix}_n"].sum())}
                    excl, past = bt["lo90"] > 1.0, bt["lo90"] > 1.10
                else:
                    bt = bootstrap_recall_ratio(rows, "M", "B0", prefix)
                    s = rows.groupby("model")[[f"{prefix}_tp", f"{prefix}_low"]].sum()
                    counts = {"n_low": int(s.loc["B0", f"{prefix}_low"]), "flagged_M": int(s.loc["M", f"{prefix}_tp"]),
                              "flagged_B0": int(s.loc["B0", f"{prefix}_tp"])}
                    excl, past = bt["hi90"] < 1.0, bt["hi90"] < 0.90
                out.append({"variant": v, "season": cls, "kind": kind, "metric": metric, "M": e[f"{metric}_M"],
                            "B0": e[f"{metric}_B0"], **counts, "bootstrap": bt, "interval_excludes_1": bool(excl),
                            "interval_past_bar": bool(past)})
    return out


def _brief(df: pd.DataFrame) -> dict:
    """Pooled readings and daily-minimum map numbers of one model's rows, with the counts."""
    s = df.sum(numeric_only=True)
    rn, dn = s["r_all_n"], s["d_n"]
    return {"readings": {"n": int(rn), "rmse": float(np.sqrt(s["r_all_sse"] / rn)), "n_low": int(s["r_all_low"]),
                         "flagged_low": int(s["r_all_tp"]), "false_alarms": int(s["r_all_fp"])},
            "daily_min_map": {"n": int(dn), "rmse": float(np.sqrt(s["d_sse"] / dn)), "coverage90": float(s["d_in90"] / dn),
                              "n_low": int(s["d_low"]), "flagged_low": int(s["d_tp"]), "false_alarms": int(s["d_fp"])}}


def _world_summary(rows: pd.DataFrame) -> dict:
    sel = {"jan_mar_to_jul_sep": rows.test == "jan_mar_to_jul_sep", "jul_sep_to_dec": rows.test == "jul_sep_to_dec",
           "rolling_all_months": rows.test == "rolling",
           "rolling_warming": (rows.test == "rolling") & (rows.season == "warming"),
           "rolling_cooling": (rows.test == "rolling") & (rows.season == "cooling")}
    return {k: {m: _brief(rows[mask & (rows.model == m)]) for m in MECH_MODELS} for k, mask in sel.items()}


def _init_mech_worker(net, cache_dir, temps_m, workdir):
    """Pool initializer of mechanism_checks: as _init_worker, without the plus or minus 2 C blocks."""
    warnings.filterwarnings("ignore")
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    _W.update(net=net, cache_dir=cache_dir,
              M=covariate_bank(sc, temps_m, cache_dir=cache_dir, cache="read"),
              M1b=covariate_bank(sc, temps_m, wall_mode="mass_transfer_only", cache_dir=cache_dir, cache="read"))


def _mech_task(args):
    kind, net, seed, law, extra = args
    if kind == "world":
        out = run_task(net, seed, "V1", _W["M"], _W["M1b"], None, _W["cache_dir"], models=set(MECH_MODELS),
                       temp_errors=(), truth_wall_law=law)
        return {"kind": "world", "law": law, "seed": seed, "rows": [r for r in out["rows"] if r["model"] in MECH_MODELS]}
    return {"kind": "oracle", "seed": seed, **_oracle_wall_fit(net, seed, *extra)}


def _oracle_wall_fit(net: str, seed: int, blocks: dict, theta_w: float, fit_ms=(1, 2, 3), target: int = 7) -> dict:
    """V1, one seed: B0, M, M1b and M with the oracle-wall bank, fitted on fit_ms and predicting the target month's
    daily-minimum map at the junctions not sampled in fit_ms.  Per H block of the banked models: posterior mass, the
    bias of exp(the block's posterior-mean ln daily minimum) against the truth (no GP, no Monte Carlo), and the share
    of the block's mass at the top of the kw grid."""
    from .experiment import NET_TRUTH
    from .features import build_features
    from .pilot import synthetic_log
    from .simulate import nominal_scenario
    cd = _W["cache_dir"]
    sc = nominal_scenario(net, 14, _dose(net))
    X = build_features(sc)
    sched = plant_schedule("V1")
    log, _, _, truths = synthetic_log(net, months=target, per_month=PER_MONTH, rotating=ROTATING, seed=seed,
                                      schedule=sched, return_truth=True, **NET_TRUTH.get(net, {}))
    log = log.assign(mi=[int(m[-2:]) for m in log.month])
    train = log[log.mi.isin(fit_ms)]
    uns = [j for j in sc.junctions if j not in set(train.junction)]
    jj = np.array([list(sc.junctions).index(j) for j in uns])
    truth = truths[target - 1]["truth_daily_min"].loc[uns].values.astype(float)
    low = truth < THRESHOLD
    T_t = sched[target - 1]["temp_C"]
    bm = _W["M"]
    bank_o = Bank(list(bm.params), bm.h0, blocks, bm.ers, "oracle_wall", bm.grid)
    models = {"B0": SimGP24(sc, X, seed=seed, cache_dir=cd, threshold=THRESHOLD).fit(train[["junction", "hour", "y"]])}
    s = train[["junction", "hour", "y", "temp_C"]]
    for name, bank in (("M", bm), ("M1b", _W["M1b"]), ("M_oracle_wall", bank_o)):
        models[name] = SeasonalSimGP24(sc, X, bank, seed=seed, cache_dir=cd, threshold=THRESHOLD).fit(s, target_temp_C=T_t)
    res = {}
    for name, mod in models.items():
        pm = mod.predict_daily_min()
        p = pm["median"].loc[uns].values.astype(float)
        f = (pm["p_below"] > 0.5).loc[uns].values.astype(bool)
        ent = {"n": int(len(truth)), "sse": float(((p - truth) ** 2).sum()), "rmse": float(np.sqrt(((p - truth) ** 2).mean())),
               "bias_mean_pred_minus_truth": float((p - truth).mean()), "n_low": int(low.sum()),
               "flagged_low": int((low & f).sum()), "false_alarms": int((~low & f).sum())}
        if isinstance(mod, SeasonalSimGP24):
            W, n0 = mod.W_, mod.bank.n_base
            P = np.asarray(mod.base_params, dtype=float)
            top = P[:, 1] == P[:, 1].max()
            zmin = mod.Z.min(axis=1)
            blk = {}
            for b, h in enumerate(mod.bank.hypotheses):
                sl = slice(b * n0, (b + 1) * n0)
                wb = W[sl].sum(axis=1)
                mass = float(wb.sum())
                mz = (wb @ zmin[sl] + W[sl].sum(axis=0) @ mod.offs_) / mass
                blk[h] = {"posterior_mass": mass, "bias_block_mean": float((np.exp(mz[jj]) - truth).mean()),
                          "mass_share_at_top_kw": float(wb[top].sum() / mass), "kw_top_m_per_day": float(P[top, 1][0])}
            ent.update(kb20=mod.kb20(), kw20=mod.rate20(1), blocks=blk)
        res[name] = ent
    return {"theta_w_true": theta_w, "fit_months": list(fit_ms), "target_month": target, "target_temp_C": T_t,
            "fit_temps_C": [sched[m - 1]["temp_C"] for m in fit_ms], "models": res}


def mechanism_checks(net: str = "Net3", n_seeds: int = 8, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache",
                     workers: int | None = None) -> dict:
    """After the stop (not pre-registered): (1) the counts and bootstrap interval of every A7 trigger, both networks;
    (2) the same experiment code (run_task, V1, B0, M and M1b, no plus or minus 2 C test) on seeds 0 to n_seeds - 1 in
    three worlds that differ only in the truth's wall temperature law (WALL_WORLDS); the pre-registered world must
    give the committed rows; (3) on seed 0, M with the truth's own wall law in its bank (OracleWallCondition), in
    memory only; (4) from the committed rows, theta_w against M's out-of-season gap and H0's posterior mass after a
    January-to-March fit.  Writes outdir/task10_mechanism_checks.json."""
    from .chemistry import FREE_CHLORINE
    from .simulate import hidden_chem_draws, nominal_scenario
    t0 = time.time()
    cache_dir = os.path.abspath(cache_dir)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[k] = "1"
    sc = nominal_scenario(net, 14, _dose(net))
    temps_m, _ = all_bank_temps(VARIANTS)
    seeds = list(range(n_seeds))
    sched = plant_schedule("V1")
    theta0 = hidden_chem_draws(FREE_CHLORINE, 0)["theta_w"]
    o_temps = sorted({_tkey(sched[m - 1]["temp_C"]) for m in (1, 2, 3, 7)})
    conds = [OracleWallCondition(T, float(e), theta0) for T in o_temps for e in ER_HYPOTHESES_K]
    blocks = {(c.temp_C, c.er_K): Z for c, Z in zip(conds, _build_conditions(sc, conds, "full", None, [cache_dir]))}
    print(f"{net}: oracle-wall blocks built in {time.time() - t0:.0f} s", flush=True)
    tasks = ([("oracle", net, 0, "theta_w", (blocks, theta0))]
             + [("world", net, s, law, None) for law in WALL_WORLDS for s in seeds])
    workdir = tempfile.mkdtemp(prefix="rm_mech_")
    try:
        with ProcessPoolExecutor(max_workers=workers or max(1, (os.cpu_count() or 2) - 1), initializer=_init_mech_worker,
                                 initargs=(net, cache_dir, temps_m, workdir)) as ex:
            results = list(ex.map(_mech_task, tasks))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    rows = pd.DataFrame([{**r, "wall_world": res["law"]} for res in results if res["kind"] == "world" for r in res["rows"]])
    oracle = next(res for res in results if res["kind"] == "oracle")
    committed = pd.read_csv(os.path.join(outdir, f"season_{net}.csv"), float_precision="round_trip")
    committed = committed[(committed.variant == "V1") & (committed.temp_error == "none") & committed.seed.isin(seeds)
                          & committed.model.isin(MECH_MODELS)]
    keys = ["test", "seed", "month", "model"]
    mine = rows[rows.wall_world == "theta_w"]
    both = mine.merge(committed, on=keys, suffixes=("", "_c"), how="outer", indicator=True)
    diff = 0.0                       # NaN on both sides counts as equal, on one side as a mismatch (inf)
    for c in MECH_COMPARE_COLS:
        a, b = both[c].astype(float).values, both[f"{c}_c"].astype(float).values
        d = np.where(np.isnan(a) & np.isnan(b), 0.0, np.abs(a - b))
        d = np.where(np.isnan(d), np.inf, d)
        diff = max(diff, float(d.max()) if len(d) else 0.0)
    worlds = {law: {"truth_wall": desc, **_world_summary(rows[rows.wall_world == law])} for law, desc in WALL_WORLDS.items()}
    # from the committed rows: per seed, theta_w against M's out-of-season map RMSE gap, and H0's mass after the fit
    full = pd.read_csv(os.path.join(outdir, f"season_{net}.csv"), float_precision="round_trip")
    e = full[(full.variant == "V1") & (full.test == "jan_mar_to_jul_sep") & (full.temp_error == "none")]
    g = e.groupby(["seed", "model"])[["d_sse", "d_n"]].sum()
    rm = np.sqrt(g.d_sse / g.d_n).unstack()
    th = e.groupby("seed").theta_w.first()
    gap = rm["M"] - rm["B0"]
    ph0 = e[e.model == "M"].groupby("seed").P_H0.first()
    out = {"generated_by": f"python -m residualmap.seasonal {net} {n_seeds} --mechanism",
           "about": "Simulation only. Checks run after the task-10 A7 stop, not pre-registered and never used to choose "
                    "a model: the counts and paired bootstrap interval of every A7 trigger, and direct tests of the "
                    "explanation for the stop (the truth's pipe wall responds to temperature more weakly than M assumes).",
           "a7_triggers": {n: a7_trigger_intervals(n, outdir) for n in ("Net3", "Net2")
                           if os.path.exists(os.path.join(outdir, f"season_{n}.csv"))},
           "wall_worlds": {"network": net, "variant": "V1", "seeds": seeds, "models": list(MECH_MODELS),
                           "same_code": "seasonal.run_task with temp_errors=() and simulate.build_scenario's truth_wall_law; "
                                        "every draw, log and model is the same in the three worlds, only the truth's wall "
                                        "temperature law differs",
                           "pre_registered_world_equals_committed_rows": {
                               "rows_matched": int((both._merge == "both").sum()),
                               "rows_unmatched": int((both._merge != "both").sum()),
                               "max_abs_diff": diff, "columns": list(MECH_COMPARE_COLS)},
                           "worlds": worlds},
           "oracle_wall_seed0": {"network": net, "variant": "V1",
                                 "bank": "M's bank with every H_E block's wall rate times the truth's own theta_w^(T - 20) "
                                         "instead of f(T; E); bulk f(T; E) and mass transfer at T as in M; built in memory",
                                 "committed_M_sse_this_seed_and_month": float(full[(full.variant == "V1") & (full.test == "jan_mar_to_jul_sep")
                                                                                   & (full.temp_error == "none") & (full.seed == 0)
                                                                                   & (full.month == 7) & (full.model == "M")].d_sse.iloc[0]),
                                 **{k: v for k, v in oracle.items() if k not in ("kind", "seed")}},
           "from_committed_rows": {
               "theta_w_vs_jan_mar_to_jul_sep_map_rmse_gap": {
                   "n_seeds": int(len(th)), "corr_theta_w_vs_rmse_M_minus_rmse_B0": float(np.corrcoef(th, gap.loc[th.index])[0, 1]),
                   "seeds_where_M_beats_B0": [int(x) for x in gap[gap < 0].index],
                   "per_seed": {int(s): {"theta_w": float(th[s]), "rmse_M": float(rm.loc[s, "M"]), "rmse_B0": float(rm.loc[s, "B0"])}
                                for s in th.index}},
               "P_H0_after_jan_mar_fit": {"mean": float(ph0.mean()), "min": float(ph0.min()), "max": float(ph0.max()),
                                          "prior": 0.25, "n_seeds": int(len(ph0))}}}
    path = os.path.join(outdir, "task10_mechanism_checks.json")
    with open(path, "w") as fh:
        json.dump(_clean(out), fh, indent=1)
    print(f"{net}: mechanism checks written to {path} in {time.time() - t0:.0f} s; pre-registered world vs committed "
          f"rows: {out['wall_worlds']['pre_registered_world_equals_committed_rows']}", flush=True)
    return out


# ============================================================================= task 10b: the wall's response, learned
# Task 10's M scaled each member's wall rate by the same Arrhenius factor as its bulk rate and stopped at its stop rule,
# mostly because the test's truth has a weaker wall response (journal, task 10, finding 2).  M2 assumes neither M's
# wall law nor M1b's: its hypotheses cross the bulk E/R with a wall E/R of its own, so the grab samples weigh how
# strongly the wall responds.  M (wall E/R = bulk E/R) and M1b (wall E/R 0) are nested in it, block for block.
# python -m residualmap.seasonal <net> <n> --wall: journal task 10b, outputs/chem/summary_season2_<net>.json.
WALL_ER_HYPOTHESES_K = (0.0, 2500.0, 5000.0, 8000.0, 12000.0)
# 0 = no wall chemistry response (M1b's structure); 5000, 8000 and 12000 K = the bulk hypotheses (M's structure on the
# diagonal); 2500 K fills the gap.  Their 10-to-20 C wall ratios are 1.00, 1.35, 1.83, 2.62 and 4.24, about even steps
# on a log scale (0.30 to 0.48 in ln).  No published wall temperature response exists (Cejas, Diaz & Gonzalez 2026;
# Lee et al. 2014 give only the direction), so the grid spans no response to the steepest bulk value.
WALL_WORLDS_10B = {"W1": "theta_w", "W2": "arrhenius"}     # simulate.build_scenario's truth_wall_law
WALL_WORLD_LABELS = {
    "W1": "task 10's truth: wall factor theta_w^(T - 20), theta_w ~ U(1.00, 1.07) per seed (an ASSUMPTION, a weak "
          "wall response; about 0 to 5600 K as a wall E/R)",
    "W2": "the truth's wall follows its own bulk Arrhenius factor f(T; E_true), E_true ~ U(4660, 12104) K"}
MODELS_10B = ("B0", "M", "M1b", "M2", "oracle")
FIRST_SEED_10B = {"Net3": 16, "Net2": 8}           # fresh seeds; task 10 scored Net3 0 to 15 and Net2 0 to 7
TASK10_SEEDS = {"Net3": range(16), "Net2": range(8)}
A7_CLASSES_10B = (("rolling", "warming"), ("rolling", "cooling"), ("jan_mar_to_jul_sep", "all_months"),
                  ("jul_sep_to_dec", "all_months"))
SUM_COLS_10B = ("n_true_low_all", "r_all_n", "r_all_sse", "r_all_sae", "r_all_low", "r_all_tp", "r_all_fp",
                "r_all_in50", "r_all_in80", "r_all_in90", "r_all_in95", "r_seen_n", "r_seen_in90", "r_new_n",
                "r_new_in90", "d_n", "d_sse", "d_sae", "d_in90", "d_low", "d_tp", "d_fp", "n_flag_all")
KEY_COLS_10B = ("world", "variant", "test", "temp_error", "seed", "season", "model")
CSV_FLOAT_10B = "%.6g"           # the compact CSVs keep 6 significant digits; the summary is computed from them


@dataclass
class WallBank(Bank):
    """A Bank whose hypotheses are (bulk E/R, wall E/R) pairs: ers holds the pairs, blocks[(T, (b, w))] is the grid
    re-run at T with bulk rates times f(T; b), wall rates times f(T; w) and EPANET's viscosity and diffusivity at T
    (chemistry.Chemistry(temp_C=T, er_K=b, wall_er_K=w)).  A pair (b, b) is task 10's M block and (b, 0) its M1b
    block: the same condition, the same cache file, the same array."""

    @property
    def hypotheses(self) -> tuple:
        return ("H0",) + tuple(f"B{b:g}_W{w:g}" for b, w in self.ers)

    @property
    def bulk_ers(self) -> tuple:
        return tuple(sorted({b for b, _ in self.ers}))

    @property
    def wall_ers(self) -> tuple:
        return tuple(sorted({w for _, w in self.ers}))

    def nested(self, which: str) -> Bank:
        """Task 10's M ('M': wall E/R equal to the bulk E/R) or M1b ('M1b': wall E/R 0) as a Bank over the same
        arrays, so a worker holds them once."""
        if which not in ("M", "M1b"):
            raise ValueError("which must be 'M' or 'M1b'")
        pick = {b: ((b, b) if which == "M" else (b, 0.0)) for b in self.bulk_ers}
        missing = [p for p in pick.values() if p not in self.ers]
        if missing:
            raise KeyError(f"the wall bank has no pairs {missing}")
        blocks = {(T, b): self.blocks[(T, pick[b])] for T in self.temps for b in self.bulk_ers}
        return Bank(self.params, self.h0, blocks, self.bulk_ers, "arrhenius" if which == "M" else "mass_transfer_only",
                    self.grid)


def wall_pairs(bulk_ers=ER_HYPOTHESES_K, wall_ers=WALL_ER_HYPOTHESES_K) -> tuple:
    """Every (bulk E/R, wall E/R) pair, bulk-major."""
    return tuple((float(b), float(w)) for b in bulk_ers for w in wall_ers)


def _wall_conditions(temps, pairs):
    return [(_tkey(T), p, Chemistry(temp_C=_tkey(T), er_K=p[0], wall_er_K=p[1]))
            for T in sorted({_tkey(t) for t in temps}) for p in pairs]


def wall_bank(sc, temps, pairs=None, grid: str = "full", cache_dir: str = "outputs/cache", cache: str = "readwrite",
              n_jobs: int | None = None) -> WallBank:
    """The bulk E/R x wall E/R bank for the given temperatures (pairs: default wall_pairs(), 3 x 5 = 15).  As
    covariate_bank: H0 is the committed grid, every block at 20 C is the committed grid, every other block is read from
    its tagged cache file or built (all in one process-pool call), and cache is 'readwrite', 'read' or 'off'.  The
    pairs (b, b) and (b, 0) are task 10's M and M1b conditions (chemistry.Chemistry stores them in that canonical
    form), so their cached task-10 files are read, never rebuilt."""
    if cache not in ("readwrite", "read", "off"):
        raise ValueError("cache must be 'readwrite', 'read' or 'off'")
    pairs = wall_pairs() if pairs is None else tuple((float(b), float(w)) for b, w in pairs)
    if len(set(pairs)) != len(pairs):
        raise ValueError("a (bulk, wall) pair is listed twice")
    params, h0 = simulator_grid_24h(sc, cache_dir, grid)
    blocks, todo = {}, {}
    for T, p, cond in _wall_conditions(temps, pairs):
        if cond.is_default:
            blocks[(T, p)] = h0
            continue
        f = grid_cache_path(sc, cache_dir, grid, cond)
        if cache != "off" and os.path.exists(f):
            with open(f, "rb") as fh:
                pp, Z = pickle.load(fh)
            if list(pp) != list(params) or Z.shape != h0.shape:
                raise ValueError(f"{f} does not hold this network's grid")
            blocks[(T, p)] = Z
        else:
            todo.setdefault(f, (cond, []))[1].append((T, p))
    if todo:
        built = _build_conditions(sc, [c for c, _ in todo.values()], grid, n_jobs, [cache_dir])
        for (f, (cond, keys)), Z in zip(todo.items(), built):
            for k in keys:
                blocks[k] = Z
            if cache == "readwrite":
                disk_preflight([cache_dir])
                with open(f, "wb") as fh:
                    pickle.dump((params, Z), fh)
    return WallBank(list(params), h0, blocks, pairs, "wall_er", grid)


def wall_bank_cache_files(sc, temps, pairs=None, grid="full", cache_dir="outputs/cache") -> list:
    """The tagged cache files a wall_bank call reads or writes (20 C excluded), task-10 files included."""
    pairs = wall_pairs() if pairs is None else pairs
    return [grid_cache_path(sc, cache_dir, grid, c) for _, _, c in _wall_conditions(temps, pairs) if not c.is_default]


def nested_prior(bank: WallBank) -> np.ndarray:
    """M2's prior over (H0, pairs...): H0 1/4, as in M; the other 3/4 split evenly over the pairs, so with the full
    3 x 5 grid each bulk E/R keeps M's 1/4 and every wall E/R gets the same share of it."""
    n = len(bank.ers)
    return np.concatenate([[0.25], np.full(n, 0.75 / n)])


class WallSeasonalSimGP24(SeasonalSimGP24):
    """Model M2 (task 10b): SeasonalSimGP24 over a WallBank, the wall's temperature response a hypothesis like the
    bulk's.  Prior: nested_prior (H0 1/4, the pairs evenly), or prior_mass, one weight per hypothesis (zeros allowed:
    all the mass on the pairs (b, b) and H0 is M; on (b, 0) and H0 it is M1b; a saved check).  Everything else
    (likelihood, discrepancy GP, set_target, daily minimum) is SeasonalSimGP24's."""

    def __init__(self, sc, X, bank: WallBank, prior_mass=None, **kw):
        if not isinstance(bank, WallBank):
            raise TypeError("WallSeasonalSimGP24 needs a WallBank")
        if "prior" in kw:
            raise TypeError("M2's prior is nested_prior or prior_mass, not a named prior")
        super().__init__(sc, X, bank, prior="uniform", **kw)
        p = nested_prior(bank) if prior_mass is None else np.asarray(prior_mass, dtype=float)
        if p.shape != (len(bank.hypotheses),) or (p < 0).any() or abs(float(p.sum()) - 1.0) > 1e-12:
            raise ValueError(f"prior_mass needs {len(bank.hypotheses)} non-negative weights summing to 1")
        self.prior = "nested" if prior_mass is None else "custom"
        self.prior_mass = p
        with np.errstate(divide="ignore"):
            self.log_prior = np.log(p[self.hyp_index] / bank.n_base)

    def _block_mass(self) -> np.ndarray:
        return self.w_.reshape(len(self.bank.hypotheses), -1).sum(axis=1)

    def _by(self, axis: int) -> dict:
        m = self._block_mass()[1:]
        vals = self.bank.bulk_ers if axis == 0 else self.bank.wall_ers
        return {v: float(sum(mi for mi, p in zip(m, self.bank.ers) if p[axis] == v)) for v in vals}

    def wall_posterior(self) -> dict:
        """Posterior mass on each wall E/R, summed over the bulk E/R (with P(H0) it sums to 1)."""
        return self._by(1)

    def bulk_posterior(self) -> dict:
        return self._by(0)

    def _mean_er(self, axis: int) -> float:
        m = self._block_mass()[1:]
        v = np.array([p[axis] for p in self.bank.ers])
        return float(m @ v / m.sum()) if m.sum() > 0 else float("nan")

    def e_posterior_mean(self) -> float:
        """Posterior mean bulk E/R over the pairs."""
        return self._mean_er(0)

    def wall_er_posterior_mean(self) -> float:
        """Posterior mean wall E/R over the pairs (the prior's is the mean of the wall grid, 5500 K)."""
        return self._mean_er(1)

    def kw_top_share(self) -> dict:
        """Within each wall E/R, the share of its posterior mass on members at the top of the kw grid (kw20 = 1.5 m/day).
        The grid is read as 20 C rates, so at 10 C a wall E/R of 12000 K reaches at most 1.5 x 0.236 = 0.35 m/day; a
        share near 1 says that hypothesis is held at the grid's edge (task 10, finding 11)."""
        n0 = self.bank.n_base
        P = np.asarray(self.base_params, dtype=float)
        top = P[:, 1] == P[:, 1].max()
        Wb = self.w_.reshape(len(self.bank.hypotheses), n0)[1:]
        out = {}
        for w in self.bank.wall_ers:
            g = Wb[[i for i, p in enumerate(self.bank.ers) if p[1] == w]]
            tot = float(g.sum())
            out[w] = float(g[:, top].sum() / tot) if tot > 0 else float("nan")
        return out

    def edge_shares(self) -> dict:
        """Added after the review of task 10b, not in the experiment's rows: within each wall E/R, the share of its
        posterior mass at the top and at the bottom of the kw grid and at the bottom of the kb grid (kb20 = 0.10), and
        the share of all the pairs' mass at the kb floor.  kw_top equals kw_top_share()."""
        P = np.asarray(self.base_params, dtype=float)
        at = {"kw_top": P[:, 1] == P[:, 1].max(), "kw_bottom": P[:, 1] == P[:, 1].min(),
              "kb_bottom": P[:, 0] == P[:, 0].min()}
        Wb = self.w_.reshape(len(self.bank.hypotheses), self.bank.n_base)[1:]
        out = {k: {} for k in at}
        for w in self.bank.wall_ers:
            g = Wb[[i for i, p in enumerate(self.bank.ers) if p[1] == w]]
            tot = float(g.sum())
            for k, sel in at.items():
                out[k][w] = float(g[:, sel].sum() / tot) if tot > 0 else float("nan")
        tot = float(Wb.sum())
        out["kb_bottom_all_pairs"] = float(Wb[:, at["kb_bottom"]].sum() / tot) if tot > 0 else float("nan")
        return out

    def posterior_columns(self) -> dict:
        out = {"P_H0": float(self._block_mass()[0])}
        out.update({f"P_W{w:g}": v for w, v in self.wall_posterior().items()})
        out.update({f"P_B{b:g}": v for b, v in self.bulk_posterior().items()})
        out["E_wall_post_mean_K"] = self.wall_er_posterior_mean()
        out.update({f"kwtop_W{w:g}": v for w, v in self.kw_top_share().items()})
        return out


TREF_K = 293.15


def wall_effective_er(world: str, theta_w: float, E_true: float) -> float:
    """The truth's wall response as a wall E/R: in W2 its own E_true; in W1 the E/R with the same 10-to-20 C ratio as
    theta_w^10, ln(theta_w) x 293.15 x 283.15 K (1.07 gives 5616 K)."""
    if world == "W2":
        return float(E_true)
    return float(math.log(theta_w) * (TREF_K * (TREF_K - 10.0)))


TRUTH_BULK_MONTH_FACTOR = (0.8, 1.2)    # simulate._chem_truth: the truth's kb is kb_net x U(0.8, 1.2) x f(T; E_true)


def kw_cap_table(net: str) -> dict:
    """The reviewer's question (task 10, finding 11): the grid's kw axis is read as 20 C rates, so at a cold month a
    steep wall E/R cannot reach a fast wall.  Per wall E/R and plant temperature: the largest and smallest nominal wall
    rate the model can represent (max and min of KW_GRID times f(T; wall E/R)) against the truth's nominal wall rate
    range at T in each world (the network's kw times the wall factor over the truth's draw range; each pipe also has its
    own roughness factor and a lognormal draw, as the members have their roughness factor).  The cap binds where the
    top is below the truth's range (some or all of it), the floor where the bottom is above it.  The same for the bulk
    rate (kb_floor; added after the review): min(KB_GRID) times f(T; bulk E/R) against the truth's kb, the network's kb
    times U(0.8, 1.2) times f(T; E_true) over E_true's range, the same in both worlds.  Nominal rates only: the
    posterior's mass at the grid's edges is measured by grid_edge_checks."""
    from .chemistry import ER_TRUTH_RANGE_K, THETA_W_TRUTH_RANGE
    from .experiment import NET_TRUTH
    from .simgp import KB_GRID, KW_GRID
    kw_net = float(NET_TRUTH.get(net, {}).get("kw_m_per_day", 0.70))
    kb_net = float(NET_TRUTH.get(net, {}).get("kb_per_day", 0.40))
    temps = sorted({_tkey(c["temp_C"]) for c in plant_schedule("V1")})
    rows, kb_rows = [], []
    for T in temps:
        w1 = sorted(kw_net * th ** (T - 20.0) for th in THETA_W_TRUTH_RANGE)
        w2 = sorted(kw_net * arrhenius(T, e) for e in ER_TRUTH_RANGE_K)
        for w in WALL_ER_HYPOTHESES_K:
            f = arrhenius(T, w)
            top, bot = max(KW_GRID) * f, min(KW_GRID) * f
            rows.append({"temp_C": T, "wall_E_K": w, "model_kw_top": top, "model_kw_bottom": bot,
                         "W1_truth_kw_range": w1, "W2_truth_kw_range": w2,
                         "W1_cap_binds_some": bool(top < w1[1]), "W1_cap_binds_all": bool(top < w1[0]),
                         "W2_cap_binds_some": bool(top < w2[1]), "W2_cap_binds_all": bool(top < w2[0]),
                         "W1_floor_binds_some": bool(bot > w1[0]), "W1_floor_binds_all": bool(bot > w1[1]),
                         "W2_floor_binds_some": bool(bot > w2[0]), "W2_floor_binds_all": bool(bot > w2[1])})
        kb = sorted(kb_net * u * arrhenius(T, e) for u in TRUTH_BULK_MONTH_FACTOR for e in ER_TRUTH_RANGE_K)
        kb = [kb[0], kb[-1]]
        for b in ER_HYPOTHESES_K:
            bot = min(KB_GRID) * arrhenius(T, b)
            kb_rows.append({"temp_C": T, "bulk_E_K": b, "model_kb_bottom": bot, "truth_kb_range": kb,
                            "floor_binds_some": bool(bot > kb[0]), "floor_binds_all": bool(bot > kb[1])})
    binds = [r for r in rows if r["W1_cap_binds_some"] or r["W2_cap_binds_some"]]
    floors = [r for r in rows if r["W1_floor_binds_some"] or r["W2_floor_binds_some"]]
    return {"kw_grid_m_per_day_at_20C": list(KW_GRID), "truth_kw_nominal_m_per_day": kw_net, "rows": rows,
            "where_the_cap_binds": [{k: r[k] for k in ("temp_C", "wall_E_K", "model_kw_top", "W1_cap_binds_some",
                                                         "W1_cap_binds_all", "W2_cap_binds_some", "W2_cap_binds_all")}
                                    for r in binds],
            "where_the_floor_binds": [{k: r[k] for k in ("temp_C", "wall_E_K", "model_kw_bottom", "W1_floor_binds_some",
                                                           "W1_floor_binds_all", "W2_floor_binds_some",
                                                           "W2_floor_binds_all")} for r in floors],
            "kb_floor": {"kb_grid_per_day_at_20C": list(KB_GRID), "truth_kb_nominal_per_day": kb_net,
                         "truth_bulk_month_factor": list(TRUTH_BULK_MONTH_FACTOR), "rows": kb_rows,
                         "where_the_floor_binds": [r for r in kb_rows if r["floor_binds_some"]]},
            "note": "The cap binds only for wall E/R steeper than the truth's own response: a truth whose wall factor "
                    "matches a hypothesis needs kw20 = the network's kw, inside the grid. The kw floor binds only for "
                    "wall E/R flatter than the truth's, and only where the network's kw is within a factor of about 2 "
                    "of the grid's 0.10 (Net2, 0.20 m/day). The kb floor binds where the truth's bulk rate, read at "
                    "20 C, is below 0.10 per day (Net2, 0.10 per day nominal). Nominal rates only; the posterior's "
                    "mass at the edges is in outputs/chem/task10b_grid_edges_<net>.json (added after the review). "
                    "Extending the axes would change every block (a new grid) and is not done; the task-10 blocks are "
                    "reused."}


# ---- the experiment
def _init_wall_worker(net, cache_dir, temps_m, pm2_path, workdir):
    """Pool initializer of run_wall: as _init_worker, with M2's bank (its plus or minus 2 C blocks merged in when
    pm2_path is given) and M and M1b taken from it (the same arrays)."""
    warnings.filterwarnings("ignore")
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass
    d = os.path.join(workdir, f"w{os.getpid()}")
    os.makedirs(d, exist_ok=True)
    os.chdir(d)
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    wb = wall_bank(sc, temps_m, cache_dir=cache_dir, cache="read")
    M, M1b = wb.nested("M"), wb.nested("M1b")
    if pm2_path:
        with open(pm2_path, "rb") as fh:
            wb = wb.merged(pickle.load(fh))
    _W.update(net=net, cache_dir=cache_dir, M=M, M1b=M1b, M2=wb)


def _wall_task(args):
    net, seed, variant, world = args
    t0 = time.time()
    out = run_task(net, seed, variant, _W["M"], _W["M1b"], None, _W["cache_dir"], models=set(MODELS_10B),
                   truth_wall_law=WALL_WORLDS_10B[world], bank_m2=_W["M2"])
    return {"net": net, "seed": seed, "variant": variant, "world": world, "seconds": time.time() - t0,
            "rows": [{"world": world, **r} for r in out["rows"]],
            "monthly": [{"world": world, **r} for r in out["monthly"]]}


def compact_rows(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(agg, windows): agg sums every count and squared-error column over the target months of each (world, variant,
    test, temperature error, seed, season class, model), which is all that pooling and the bootstrap over seeds need;
    windows is one row per rolling fit (world, variant, seed, target month) with the fitted rates and M2's posterior,
    for A5 and R2."""
    df = df[df.model.isin(MODELS_10B)]      # run_task also scores persistence and the network mean out of season
    keys = list(KEY_COLS_10B)
    g = df.groupby(keys, sort=False)
    agg = g[list(SUM_COLS_10B)].sum(min_count=1)
    agg["n_months"] = g.size()
    agg[["E_true_K", "theta_w"]] = g[["E_true_K", "theta_w"]].first()
    agg = agg.reset_index()
    r = df[(df.test == "rolling") & df.model.isin(["B0", "M", "M1b", "M2"])]
    want = {"B0": ["map_kb", "post_kb", "post_kw"], "M": ["post_kb", "post_kw", "P_H0", "E_post_mean_K"],
            "M1b": ["post_kb", "post_kw", "P_H0"],
            "M2": ["post_kb", "post_kw", "P_H0", "E_post_mean_K", "E_wall_post_mean_K"]
                  + [c for c in df.columns if c.startswith(("P_W", "P_B", "kwtop_W"))]}
    idx = ["world", "variant", "seed", "month"]
    win = r[r.model == "B0"][idx + ["plant_temp_C", "E_true_K", "theta_w"]].set_index(idx)
    for model, cols in want.items():
        part = r[r.model == model].set_index(idx)[cols]
        win = win.join(part.add_prefix(f"{model}_"))
    win = win.reset_index()
    win["wall_E_eff_K"] = [wall_effective_er(w, th, e) for w, th, e in zip(win.world, win.theta_w, win.E_true_K)]
    return agg, win


def pooled_10b(df: pd.DataFrame) -> dict:
    """Pooled metrics of a set of aggregated rows (sums over seeds and months)."""
    s = df.sum(numeric_only=True)
    out = {"n_rows": int(len(df)), "n_months": int(s["n_months"])}
    n, low = s["r_all_n"], s["r_all_low"]
    r = {"n": int(n), "rmse": float(np.sqrt(s["r_all_sse"] / n)) if n else float("nan"),
         "mae": float(s["r_all_sae"] / n) if n else float("nan"), "n_low": int(low),
         "recall": float(s["r_all_tp"] / low) if low else float("nan"), "false_alarms": int(s["r_all_fp"])}
    for q in LEVELS:
        col = df[f"r_all_in{q}"]
        r[f"coverage{q}"] = float(col.sum() / n) if n and col.notna().all() else float("nan")
    for tag in ("seen", "new"):
        nn, col = s[f"r_{tag}_n"], df[f"r_{tag}_in90"]
        r[f"n_{tag}"] = int(nn)
        r[f"coverage90_{tag}"] = float(col.sum() / nn) if nn and col.notna().all() else float("nan")
    out["readings"] = r
    if df["d_n"].notna().all() and s["d_n"] > 0:
        dn, tp, fp = s["d_n"], s["d_tp"], s["d_fp"]
        out["daily_min_map"] = {"n": int(dn), "rmse": float(np.sqrt(s["d_sse"] / dn)), "mae": float(s["d_sae"] / dn),
                                "coverage90": float(df.d_in90.sum() / dn) if df.d_in90.notna().all() else float("nan"),
                                "n_low": int(s["d_low"]),
                                "recall": float(tp / s["d_low"]) if s["d_low"] else float("nan"),
                                "false_alarms": int(fp), "precision": float(tp / (tp + fp)) if tp + fp > 0 else float("nan")}
    return out


def _cls_rows(df: pd.DataFrame, cls: str) -> pd.DataFrame:
    return df if cls == "all_months" else df[df.season == cls]


def pooled_tree_10b(agg: pd.DataFrame) -> dict:
    """P[world][variant][test][temperature error][class][model]; class all_months, warming or cooling."""
    P = {}
    for (world, v, test, err), g in agg.groupby(["world", "variant", "test", "temp_error"], sort=True):
        for cls in ("all_months", "warming", "cooling"):
            gc = _cls_rows(g, cls)
            if not len(gc):
                continue
            node = P.setdefault(world, {}).setdefault(v, {}).setdefault(test, {}).setdefault(err, {}).setdefault(cls, {})
            for model, gm in gc.groupby("model", sort=True):
                node[model] = pooled_10b(gm)
    return P


BOOT_PAIRS_10B = (("M2", "B0", "r_all", "readings_M2_B0"), ("M2", "B0", "d", "map_M2_B0"),
                  ("M", "B0", "r_all", "readings_M_B0"), ("M", "B0", "d", "map_M_B0"),
                  ("M1b", "B0", "r_all", "readings_M1b_B0"), ("M1b", "B0", "d", "map_M1b_B0"),
                  ("M2", "M", "d", "map_M2_M"), ("M2", "M1b", "d", "map_M2_M1b"))


def _safe_recall_boot(rows, num, den, prefix):
    try:
        return bootstrap_recall_ratio(rows, num, den, prefix)
    except (IndexError, ValueError, ZeroDivisionError):
        return {"ratio": float("nan"), "lo90": float("nan"), "hi90": float("nan"), "n_seeds": None}


def bootstrap_tree_10b(agg: pd.DataFrame) -> dict:
    """Paired bootstrap over seeds (2000 resamples, 90% intervals) of pooled RMSE ratios, temperature error 'none'."""
    B = {}
    a = agg[agg.temp_error == "none"]
    for (world, v, test), g in a.groupby(["world", "variant", "test"], sort=True):
        for cls in ("all_months", "warming", "cooling"):
            gc = _cls_rows(g, cls)
            if not len(gc):
                continue
            node = B.setdefault(world, {}).setdefault(v, {}).setdefault(test, {}).setdefault(cls, {})
            for num, den, prefix, key in BOOT_PAIRS_10B:
                if {num, den} <= set(gc.model):
                    node[key] = bootstrap_ratio(gc, num, den, prefix)
    return B


def _a7_entry(P_cls: dict, rows: pd.DataFrame, num: str = "M2", den: str = "B0") -> tuple[dict, list, list]:
    """One class of the stop rule: the table entry, its triggers (with counts and bootstrap intervals) and the false
    alarm rises above 10%, rises from zero included (ratio null, from_zero true; added after the review, the scored
    run's summary left them out)."""
    ent, trig, fa = {}, [], []
    for kind, key, prefix in (("readings", "readings", "r_all"), ("daily_min_map", "daily_min_map", "d")):
        if key not in P_cls[num] or key not in P_cls[den]:
            continue
        m, b = P_cls[num][key], P_cls[den][key]
        rr = m["rmse"] / b["rmse"]
        rc = m["recall"] / b["recall"] if b["recall"] and not math.isnan(b["recall"]) else float("nan")
        far = m["false_alarms"] / b["false_alarms"] if b["false_alarms"] else float("nan")
        ent[kind] = {"rmse_M2": m["rmse"], "rmse_B0": b["rmse"], "rmse_ratio": rr, "recall_M2": m["recall"],
                     "recall_B0": b["recall"], "recall_ratio": rc, "n_low": b["n_low"], "flagged_M2": None,
                     "flagged_B0": None, "false_alarms_M2": m["false_alarms"], "false_alarms_B0": b["false_alarms"],
                     "false_alarm_ratio": far}
        s = rows.groupby("model")[[f"{prefix}_tp"]].sum()
        ent[kind]["flagged_M2"], ent[kind]["flagged_B0"] = int(s.loc[num].iloc[0]), int(s.loc[den].iloc[0])
        if rr > 1.10:
            bt = bootstrap_ratio(rows, num, den, prefix)
            trig.append({"kind": kind, "metric": "rmse", "M2": m["rmse"], "B0": b["rmse"], "ratio": rr,
                         "n_scored": int(rows[rows.model == num][f"{prefix}_n"].sum()), "bootstrap": bt,
                         "interval_excludes_1": bool(bt["lo90"] > 1.0), "interval_past_bar": bool(bt["lo90"] > 1.10)})
        if not math.isnan(rc) and rc < 0.90:
            bt = _safe_recall_boot(rows, num, den, prefix)
            trig.append({"kind": kind, "metric": "recall", "M2": m["recall"], "B0": b["recall"], "ratio": rc,
                         "n_low": b["n_low"], "flagged_M2": ent[kind]["flagged_M2"], "flagged_B0": ent[kind]["flagged_B0"],
                         "bootstrap": bt, "interval_excludes_1": bool(bt["hi90"] < 1.0),
                         "interval_past_bar": bool(bt["hi90"] < 0.90)})
        from_zero = not b["false_alarms"] and m["false_alarms"] > 0      # a rise from none has no ratio; still a rise
        if from_zero or (not math.isnan(far) and far > 1.10):
            fa.append({"kind": kind, "M2": m["false_alarms"], "B0": b["false_alarms"], "ratio": far,
                       "from_zero": bool(from_zero)})
    return ent, trig, fa


def acceptance_10b(net: str, world: str, Pw: dict, aggw: pd.DataFrame, Bw: dict, kbsd: dict) -> dict:
    """The pre-registered bars of task 10b for one truth world (journal, task 10b, written before the first scored
    run): task 10's A2 to A6 with M2 in place of M, A7 on M2 against B0 over the rolling and out-of-season classes, R1."""
    A = {"A1": {"bar": "outputs/chem/checks_report.json: every check passes, the task-10b checks included",
                "pass": None, "judged_from": "outputs/chem/checks_report.json after this run"}}
    c = {}
    for kind, key in (("readings", "readings"), ("daily_min_map", "daily_min_map")):
        m, b = Pw["V2"]["rolling"]["none"]["all_months"]["M2"][key], Pw["V2"]["rolling"]["none"]["all_months"]["B0"][key]
        c[kind] = {"rmse_ratio_M2_B0": m["rmse"] / b["rmse"], "recall_M2": m["recall"], "recall_B0": b["recall"],
                   "recall_diff": m["recall"] - b["recall"]}
    ok2 = all(abs(v["rmse_ratio_M2_B0"] - 1) <= 0.05 and abs(v["recall_diff"]) <= 0.03 for v in c.values())
    A["A2"] = {"bar": "V2 control, rolling: |RMSE_M2/RMSE_B0 - 1| <= 0.05 and |recall_M2 - recall_B0| <= 0.03, for "
                      "held-out readings and for the daily-minimum map", **c, "pass": bool(ok2)}
    m, b = Pw["V1"]["rolling"]["none"]["all_months"]["M2"]["readings"], Pw["V1"]["rolling"]["none"]["all_months"]["B0"]["readings"]
    bt = Bw["V1"]["rolling"]["all_months"]["readings_M2_B0"]
    A["A3"] = {"bar": "V1 rolling: pooled held-out-reading RMSE_M2 <= RMSE_B0; a gain is claimed only if the bootstrap "
                      "90% upper bound of RMSE_M2/RMSE_B0 is below 1.0",
               "rmse_M2": m["rmse"], "rmse_B0": b["rmse"], "bootstrap": bt, "pass": bool(m["rmse"] <= b["rmse"]),
               "gain_claimed": bool(bt["hi90"] < 1.0),
               "readme_wording": "measurable gain" if bt["hi90"] < 1.0 else f"no measurable gain on {net}"}
    e = Pw["V1"]["jan_mar_to_jul_sep"]["none"]["all_months"]
    mM, mB = e["M2"]["daily_min_map"], e["B0"]["daily_min_map"]
    r = Pw["V1"]["jul_sep_to_dec"]["none"]["all_months"]
    fa_M, fa_B = r["M2"]["daily_min_map"]["false_alarms"], r["B0"]["daily_min_map"]["false_alarms"]
    parts = {"rmse": mM["rmse"] <= 0.85 * mB["rmse"], "recall": mM["recall"] >= mB["recall"] + 0.10,
             "coverage90": mM["coverage90"] >= 0.80, "reverse_false_alarms": fa_M <= fa_B}
    A["A4"] = {"bar": "V1, fit January to March, predict July to September, daily-minimum map: RMSE_M2 <= 0.85 RMSE_B0, "
                      "recall_M2 >= recall_B0 + 0.10, coverage90_M2 >= 0.80; fit July to September, predict December: "
                      "false alarms_M2 <= false alarms_B0 (the bar is set for Net3)",
               "applies": net == "Net3", "rmse_M2": mM["rmse"], "rmse_B0": mB["rmse"], "rmse_ratio": mM["rmse"] / mB["rmse"],
               "recall_M2": mM["recall"], "recall_B0": mB["recall"], "coverage90_M2": mM["coverage90"],
               "coverage90_B0": mB["coverage90"], "false_alarms_M2_jul_sep": mM["false_alarms"],
               "false_alarms_B0_jul_sep": mB["false_alarms"], "reverse_false_alarms_M2": fa_M,
               "reverse_false_alarms_B0": fa_B, "parts": parts, "pass": bool(all(parts.values()))}
    A["A5"] = {"bar": "V1 rolling: median over seeds of the SD of M2's kb20 across the 9 windows < the same for B0's MAP kb",
               **kbsd, "pass": bool(kbsd["median_sd_M2_kb20"] < kbsd["median_sd_B0_map_kb"])}
    base_rec = mM["recall"]
    errs = {}
    for dT in TEMP_ERRORS_C:
        for mode in ("forecast_month", "every_month"):
            k = f"{dT:+g}C_{mode}"
            rec = Pw["V1"]["jan_mar_to_jul_sep"][k]["all_months"]["M2"]["daily_min_map"]["recall"]
            errs[k] = {"recall": rec, "cost": base_rec - rec,
                       "false_alarms": Pw["V1"]["jan_mar_to_jul_sep"][k]["all_months"]["M2"]["daily_min_map"]["false_alarms"]}
    A["A6"] = {"bar": "in the A4 test, a logged temperature off by +2 C or -2 C (forecast month only, or every logged month) "
                      "costs at most 0.10 of M2's daily-minimum recall", "recall_no_error": base_rec,
               "false_alarms_no_error": mM["false_alarms"], "errors": errs,
               "worst_cost": max(v["cost"] for v in errs.values()),
               "pass": bool(all(v["cost"] <= 0.10 for v in errs.values()))}
    trig, fa, table = [], [], {}
    for v in [x for x in VARIANTS if x in Pw]:
        for test, cls in A7_CLASSES_10B:
            P_cls = Pw[v][test]["none"][cls]
            rows = _cls_rows(aggw[(aggw.variant == v) & (aggw.test == test) & (aggw.temp_error == "none")], cls)
            ent, t, f = _a7_entry(P_cls, rows)
            key = f"{v}|{test}|{cls}"
            table[key] = ent
            trig += [{"variant": v, "test": test, "season": cls, **x} for x in t]
            fa += [{"variant": v, "test": test, "season": cls, **x} for x in f]
    roll = [x for x in trig if x["test"] == "rolling"]
    A["A7"] = {"bar": "stop rule, M2 against B0, each variant (V1, V2, V3) and class (rolling warming, rolling cooling, "
                      "January-to-March fit predicting July to September, July-to-September fit predicting December): "
                      "stop if RMSE_M2 > 1.10 RMSE_B0 or recall_M2 < 0.90 recall_B0, for held-out readings or the "
                      "daily-minimum map; false alarms up by more than 10% are listed as regressions",
               "by_class": table, "triggered": trig, "stop": bool(trig),
               "rolling_only_as_task10": {"triggered": roll, "stop": bool(roll)},
               "robust_triggers": [x for x in trig if x["interval_excludes_1"]],
               "triggers_within_noise": [x for x in trig if not x["interval_excludes_1"]],
               "false_alarm_rises_over_10pct": fa}
    r1, fails = {}, []
    for v in [x for x in VARIANTS if x in Pw]:
        for test, cls in A7_CLASSES_10B:
            g = Pw[v][test]["none"][cls]
            d = {k: g[k]["daily_min_map"] for k in ("M2", "M", "M1b")}
            best_rmse = min(d["M"]["rmse"], d["M1b"]["rmse"])
            rec = [x for x in (d["M"]["recall"], d["M1b"]["recall"]) if not math.isnan(x)]
            best_rec = max(rec) if rec else float("nan")
            ok_r = d["M2"]["rmse"] <= 1.10 * best_rmse
            ok_c = math.isnan(best_rec) or math.isnan(d["M2"]["recall"]) or d["M2"]["recall"] >= 0.90 * best_rec
            key = f"{v}|{test}|{cls}"
            r1[key] = {"rmse_M2": d["M2"]["rmse"], "rmse_M": d["M"]["rmse"], "rmse_M1b": d["M1b"]["rmse"],
                       "rmse_ratio_to_best": d["M2"]["rmse"] / best_rmse, "recall_M2": d["M2"]["recall"],
                       "recall_M": d["M"]["recall"], "recall_M1b": d["M1b"]["recall"],
                       "recall_ratio_to_best": d["M2"]["recall"] / best_rec if best_rec else float("nan"),
                       "n_low": d["M2"]["n_low"], "rmse_ok": bool(ok_r), "recall_ok": bool(ok_c)}
            if not ok_r:
                fails.append(f"{key}: map RMSE {d['M2']['rmse']:.4f} vs best of M and M1b {best_rmse:.4f}")
            if not ok_c:
                fails.append(f"{key}: map recall {d['M2']['recall']:.4f} vs best of M and M1b {best_rec:.4f}")
    A["R1"] = {"bar": "M2 is never more than 10% worse than the better of M and M1b on daily-minimum map RMSE or recall, in "
                      "any variant and class of A7: RMSE_M2 <= 1.10 min(RMSE_M, RMSE_M1b), recall_M2 >= 0.90 "
                      "max(recall_M, recall_M1b)", "by_class": r1, "failures": fails, "pass": not fails}
    return A


def kb_sd_summary_10b(win: pd.DataFrame) -> dict:
    """A5: SD across the 9 rolling windows, per seed, of M2's kb20 and B0's MAP kb (V1), with like-for-like extras."""
    r = win[win.variant == "V1"]
    per = {k: r.groupby("seed")[c].std(ddof=1) for k, c in (("M2_kb20", "M2_post_kb"), ("B0_map_kb", "B0_map_kb"),
                                                             ("B0_post_kb", "B0_post_kb"), ("M_kb20", "M_post_kb"),
                                                             ("M2_kw20", "M2_post_kw"), ("B0_post_kw", "B0_post_kw"))}
    return {"median_sd_M2_kb20": float(per["M2_kb20"].median()), "median_sd_B0_map_kb": float(per["B0_map_kb"].median()),
            "like_for_like_median_sd_B0_posterior_mean_kb": float(per["B0_post_kb"].median()),
            "median_sd_M_kb20": float(per["M_kb20"].median()),
            "median_sd_M2_kw20": float(per["M2_kw20"].median()),
            "median_sd_B0_posterior_mean_kw": float(per["B0_post_kw"].median())}


def _window_span(m: int, variant: str = "V1") -> float:
    s = plant_schedule(variant)
    t = [s[x - 1]["temp_C"] for x in (m - 3, m - 2, m - 1)]
    return float(max(t) - min(t))


def r2_wall_posterior(win: pd.DataFrame, wall_ers=WALL_ER_HYPOTHESES_K) -> dict:
    """R2 (reported, no bar): does M2's posterior over the wall E/R move toward the truth's own wall response
    (wall_effective_er)?  Per world and variant, for all rolling windows, the January-to-March and July-to-September
    windows (0.5 C of contrast) and the windows spanning at least 4 C: the mean posterior mass on each wall E/R within
    the pairs (H0 excluded), the posterior mean wall E/R against the truth's, the share of fits whose posterior mean is
    closer to the truth's than the prior mean (5500 K), and their correlation.  Also, per variant, the posterior mean
    in W2 minus in W1 on the same seeds and windows: positive when the data move it the right way."""
    prior_mean = float(np.mean(wall_ers))
    sets = {"all_rolling": lambda mm: np.ones(len(mm), bool), "jan_mar_window": lambda mm: mm == 4,
            "jul_sep_window": lambda mm: mm == 10,
            "windows_spanning_4C_or_more": lambda mm: np.array([_window_span(int(x)) >= 4.0 for x in mm])}
    out = {"prior_mean_wall_E_K": prior_mean, "by_world": {}, "W2_minus_W1": {}}
    for world in sorted(win.world.unique()):
        for v in sorted(win.variant.unique()):
            w = win[(win.world == world) & (win.variant == v)]
            node = out["by_world"].setdefault(world, {}).setdefault(v, {})
            for name, sel in sets.items():
                g = w[sel(w.month.values)]
                if not len(g):
                    continue
                pe = 1.0 - g["M2_P_H0"]
                mass = {f"W{x:g}": float((g[f"M2_P_W{x:g}"] / pe).mean()) for x in wall_ers}
                post, eff = g["M2_E_wall_post_mean_K"].values, g["wall_E_eff_K"].values
                closer = np.abs(post - eff) < np.abs(prior_mean - eff)
                node[name] = {"n_fits": int(len(g)), "mean_mass_within_pairs": mass,
                              "mean_posterior_wall_E_K": float(post.mean()), "mean_truth_wall_E_eff_K": float(eff.mean()),
                              "share_closer_than_prior": float(closer.mean()),
                              "corr_posterior_vs_truth": float(np.corrcoef(post, eff)[0, 1]) if np.std(eff) > 0 and np.std(post) > 0 else float("nan"),
                              "mean_P_H0": float(g["M2_P_H0"].mean())}
    if {"W1", "W2"} <= set(win.world):
        for v in sorted(win.variant.unique()):
            a = win[(win.world == "W1") & (win.variant == v)].set_index(["seed", "month"])["M2_E_wall_post_mean_K"]
            b = win[(win.world == "W2") & (win.variant == v)].set_index(["seed", "month"])["M2_E_wall_post_mean_K"]
            d = (b - a).dropna()
            node = {"all_rolling": {"mean": float(d.mean()), "share_positive": float((d > 0).mean()), "n": int(len(d))}}
            for name in ("jan_mar_window", "windows_spanning_4C_or_more"):
                mm = d.index.get_level_values("month").values
                dd = d[sets[name](mm)]
                node[name] = {"mean": float(dd.mean()), "share_positive": float((dd > 0).mean()), "n": int(len(dd))}
            out["W2_minus_W1"][v] = node
    return out


def p_h0_summary(win: pd.DataFrame) -> dict:
    """P(H0) of M2 after the January-to-March and July-to-September fits (V1) next to the 15 C control (V2, every
    window), per world (task 10's open item: H0 kept 0.104 after a winter fit against the control's 0.136)."""
    out = {}
    for world in sorted(win.world.unique()):
        w = win[win.world == world]
        out[world] = {"V1_jan_mar_fit": float(w[(w.variant == "V1") & (w.month == 4)].M2_P_H0.mean()),
                      "V1_jul_sep_fit": float(w[(w.variant == "V1") & (w.month == 10)].M2_P_H0.mean()),
                      "V2_control_all_windows": float(w[w.variant == "V2"].M2_P_H0.mean()), "prior": 0.25}
    return out


def kw_top_observed(win: pd.DataFrame, wall_ers=WALL_ER_HYPOTHESES_K) -> dict:
    """Mean share of each wall E/R's posterior mass at the top of the kw grid, after the January-to-March fit (V1) and
    over all V1 rolling fits, per world."""
    out = {}
    for world in sorted(win.world.unique()):
        w = win[(win.world == world) & (win.variant == "V1")]
        out[world] = {name: {f"W{x:g}": float(g[f"M2_kwtop_W{x:g}"].mean()) for x in wall_ers}
                      for name, g in (("jan_mar_fit", w[w.month == 4]), ("all_rolling", w))}
    return out


def _round_sig(x, sig: int = 6):
    if isinstance(x, dict):
        return {k: _round_sig(v, sig) for k, v in x.items()}
    if isinstance(x, list):
        return [_round_sig(v, sig) for v in x]
    if isinstance(x, float) and math.isfinite(x) and x != 0.0:
        return float(f"{x:.{sig}g}")
    return x


def _truth_low_by_month(monthly: pd.DataFrame) -> dict:
    """Share of junctions whose true daily minimum is below the threshold, mean over seeds, per world, variant, month."""
    truth_low = {}
    for (w, v, m), g in monthly.groupby(["world", "variant", "month"]):
        truth_low.setdefault(w, {}).setdefault(v, {})[int(m)] = round(float(g.n_true_low_all.mean() / g.n_junctions.iloc[0]), 4)
    return truth_low


def summarise_wall(net: str, agg: pd.DataFrame, win: pd.DataFrame, truth_low: dict, seeds, settings: dict) -> dict:
    settings = dict(settings)
    kw_cap = settings.pop("kw_cap")
    P = pooled_tree_10b(agg)
    B = bootstrap_tree_10b(agg)
    kbsd = {w: kb_sd_summary_10b(win[win.world == w]) for w in sorted(win.world.unique())}
    acc = {w: acceptance_10b(net, w, P[w], agg[agg.world == w], B[w], kbsd[w]) for w in sorted(P)}
    a7 = {w: {"stop": acc[w]["A7"]["stop"], "rolling_only_stop": acc[w]["A7"]["rolling_only_as_task10"]["stop"],
              "n_triggers": len(acc[w]["A7"]["triggered"]), "n_robust": len(acc[w]["A7"]["robust_triggers"])} for w in acc}
    bands = {}
    for w in P:
        g = P[w]["V1"]["rolling"]["none"]["all_months"]["M2"]["readings"]
        bands[w] = {"seen_taps": {"coverage90": g["coverage90_seen"], "band": [0.85, 0.95]},
                    "new_taps": {"coverage90": g["coverage90_new"], "band": [0.80, 0.95]}}
        for v in bands[w].values():
            v["in_band"] = bool(v["band"][0] <= v["coverage90"] <= v["band"][1])
    return _round_sig({
        "generated_by": f"python -m residualmap.seasonal {net} {len(seeds)} --wall", "network": net,
        "seeds": [int(s) for s in seeds],
        "about": "Simulation only (task 10b). Fresh seeds. Two truth worlds that differ only in how the truth's pipe wall "
                 "decay responds to temperature (W1: task 10's weak theta_w law, an assumption; W2: the truth's own bulk "
                 "Arrhenius factor). Models: today's temperature-blind model B0, task 10's bank M (wall E/R = bulk E/R), "
                 "its bulk-only ablation M1b (wall E/R 0), M2 (bulk E/R x wall E/R, the wall's response learned) and the "
                 "oracle (best single member of M2's bank). Pooled over seeds and target months; the CSVs hold the sums "
                 "every number here is computed from. The temperature schedules are ASSUMPTIONS.",
        "settings": settings, "acceptance": acc, "a7_summary": a7,
        "adoption": {"rule": "M2 becomes the temperature model only if it passes A7 in both worlds (on every network)",
                     "passes_A7_both_worlds_on_this_network": bool(all(not x["stop"] for x in a7.values()) and len(a7) == 2)},
        "pilot_coverage_bands_M2": bands, "r2_wall_posterior": r2_wall_posterior(win), "p_h0": p_h0_summary(win),
        "kw_cap": {**kw_cap, "observed_share_at_kw_top": kw_top_observed(win)},
        "pooled": P, "bootstrap": B, "truth_share_below_threshold_daily_min_by_month": truth_low})


def run_wall(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache", workers: int | None = None,
             variants=VARIANTS, worlds=tuple(WALL_WORLDS_10B)) -> dict:
    """Task 10b's experiment (see the section comment).  One task per world, seed and variant; outputs
    outdir/season2_<net>.csv (aggregated sums), season2_windows_<net>.csv (one row per rolling fit),
    summary_season2_<net>.json and the figures season2_<net>.png and wall_posterior_<net>.png."""
    os.makedirs(outdir, exist_ok=True)
    cache_dir = os.path.abspath(cache_dir)
    t0 = time.time()
    used = sorted(set(seeds) & set(TASK10_SEEDS.get(net, ())))
    if used:
        raise ValueError(f"seeds {used} were scored in task 10 on {net}; task 10b uses fresh seeds only")
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free, below {MIN_FREE_GB_RUN} GB")
    from .simulate import nominal_scenario
    sc = nominal_scenario(net, 14, _dose(net))
    temps_m, temps_pm2 = all_bank_temps(variants)
    tb = time.time()
    n_new = sum(not os.path.exists(f) for f in wall_bank_cache_files(sc, temps_m, cache_dir=cache_dir))
    wall_bank(sc, temps_m, cache_dir=cache_dir)
    pm2 = wall_bank(sc, temps_pm2, cache_dir=cache_dir, cache="off") if temps_pm2 else None
    t_bank = time.time() - tb
    print(f"{net}: wall banks ready in {t_bank:.0f} s ({len(temps_m)} temperatures x {len(wall_pairs())} pairs, "
          f"{n_new} new cached blocks; {len(temps_pm2)} more temperatures for the +/-2 C test, not cached); "
          f"{disk_free_gb(cache_dir):.2f} GB free", flush=True)
    if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
        raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free after the banks, below {MIN_FREE_GB_RUN} GB")
    workdir = tempfile.mkdtemp(prefix="rm_season2_")
    pm2_path = os.path.join(cache_dir, f"_tmp_season2_pm2_{net}_{os.getpid()}.pkl") if pm2 is not None else None
    results = []
    try:
        if pm2 is not None:
            with open(pm2_path, "wb") as fh:
                pickle.dump(pm2, fh)
            del pm2
            if disk_free_gb(cache_dir) < MIN_FREE_GB_RUN:
                raise RuntimeError(f"stopping: {disk_free_gb(cache_dir):.2f} GB free with the +/-2 C blocks on disk")
        tasks = [(net, s, v, w) for w in worlds for s in seeds for v in variants]
        n_workers = workers or max(1, min(len(tasks), (os.cpu_count() or 2) - 1))
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            os.environ[k] = "1"
        with ProcessPoolExecutor(max_workers=n_workers, initializer=_init_wall_worker,
                                 initargs=(net, cache_dir, temps_m, pm2_path, workdir)) as ex:
            for res in ex.map(_wall_task, tasks):
                results.append(res)
                print(f"  {net} {res['world']} seed {res['seed']} {res['variant']}: {len(res['rows'])} rows in "
                      f"{res['seconds']:.0f} s", flush=True)
    finally:
        if pm2_path and os.path.exists(pm2_path):
            os.remove(pm2_path)
        shutil.rmtree(workdir, ignore_errors=True)
    df = pd.DataFrame([r for res in results for r in res["rows"]])
    order = {"rolling": 0, "jan_mar_to_jul_sep": 1, "jul_sep_to_dec": 2}
    df = df.assign(_o=df.test.map(order)).sort_values(["world", "variant", "_o", "temp_error", "seed", "month", "model"],
                                                     kind="mergesort").drop(columns="_o").reset_index(drop=True)
    monthly = pd.DataFrame([r for res in results for r in res["monthly"]]).sort_values(
        ["world", "variant", "seed", "month"]).reset_index(drop=True)
    agg, win = compact_rows(df)
    p_agg, p_win = os.path.join(outdir, f"season2_{net}.csv"), os.path.join(outdir, f"season2_windows_{net}.csv")
    agg.to_csv(p_agg, index=False, float_format=CSV_FLOAT_10B)
    win.to_csv(p_win, index=False, float_format=CSV_FLOAT_10B)
    agg = pd.read_csv(p_agg, float_precision="round_trip")       # every number of the summary comes from the files
    win = pd.read_csv(p_win, float_precision="round_trip")
    summ = summarise_wall(net, agg, win, _truth_low_by_month(monthly), seeds, _wall_settings(net, variants, worlds))
    with open(os.path.join(outdir, f"summary_season2_{net}.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    plot_season2(net, summ, os.path.join(outdir, f"season2_{net}.png"))
    plot_wall_posterior(net, win, summ, os.path.join(outdir, f"wall_posterior_{net}.png"))
    stops = {w: summ["acceptance"][w]["A7"]["stop"] for w in summ["acceptance"]}
    print(f"{net}: {time.time() - t0:.0f} s in all (banks {t_bank:.0f} s); A7 stop by world: {stops}", flush=True)
    for w in summ["acceptance"]:
        for x in summ["acceptance"][w]["A7"]["triggered"]:
            print(f"  {w} {x['variant']} {x['test']} {x['season']} {x['kind']} {x['metric']}: ratio {x['ratio']:.3f} "
                  f"(90% {x['bootstrap']['lo90']:.3f} to {x['bootstrap']['hi90']:.3f})", flush=True)
    return summ


def _wall_settings(net: str, variants, worlds) -> dict:
    """The settings block of summary_season2_<net>.json (run_wall and resummarise_wall)."""
    temps_m, temps_pm2 = all_bank_temps(variants)
    bank_hyp = ("H0",) + tuple(f"B{b:g}_W{w:g}" for b, w in wall_pairs())
    return {"months": MONTHS, "route_taps": PER_MONTH, "rotating_taps": ROTATING, "noise_sd_mgL": NOISE_SD,
            "sampling_hours": [min(DAY_HOURS), max(DAY_HOURS)], "threshold_mgL": THRESHOLD,
            "month_seed": "1000 + 100 seed + month index (pilot.synthetic_log with a schedule)",
            "schedules_ASSUMPTION": {v: plant_schedule(v) for v in variants}, "variants": VARIANT_LABELS,
            "worlds": {w: WALL_WORLD_LABELS[w] for w in worlds}, "models": list(MODELS_10B),
            "M2_hypotheses": list(bank_hyp), "wall_E_grid_K": list(WALL_ER_HYPOTHESES_K),
            "bulk_E_grid_K": list(ER_HYPOTHESES_K),
            "wall_ratio_20C_over_10C": {f"{w:g}": 1.0 / arrhenius(10.0, w) for w in WALL_ER_HYPOTHESES_K},
            "M2_prior": "H0 1/4 (as M); 3/4 split evenly over the 15 (bulk, wall) pairs (0.05 each)",
            "bank_temperatures_C": temps_m, "pm2_temperatures_C_not_cached": temps_pm2,
            "season_classes": {"warming": list(WARMING_MONTHS), "cooling": list(COOLING_MONTHS)},
            "rolling_targets": list(ROLLING_TARGETS),
            "extrapolation": {k: [list(a), list(b)] for k, (a, b) in EXTRAPOLATION.items()},
            "a7_classes": [list(x) for x in A7_CLASSES_10B],
            "map_junctions": "every junction not sampled in the 3-month window, the same set for every model",
            "oracle": "the best single member and dose of M2's bank (which holds M's and M1b's) for the target month",
            "bootstrap": {"resamples": BOOT_N, "rng_seed": BOOT_SEED, "unit": "seed"},
            "csv": {"float_format": CSV_FLOAT_10B, "aggregated_over": "the target months of each season class"},
            "truth": _truth_settings(net), "kw_cap": kw_cap_table(net)}


def resummarise_wall(net: str, outdir: str = "outputs/chem") -> dict:
    """Rewrite summary_season2_<net>.json from the committed season2_<net>.csv and season2_windows_<net>.csv, with no
    rerun (after the review of task 10b: rises of false alarms from zero listed, the kw and kb floors in kw_cap).
    Every pooled number, interval and bar comes from the two CSVs as in run_wall; the truth's monthly low shares, which
    are not in the CSVs, and the seeds are read from the summary being replaced.  The figures are not redrawn."""
    path = os.path.join(outdir, f"summary_season2_{net}.json")
    with open(path) as fh:
        old = json.load(fh)
    agg = pd.read_csv(os.path.join(outdir, f"season2_{net}.csv"), float_precision="round_trip")
    win = pd.read_csv(os.path.join(outdir, f"season2_windows_{net}.csv"), float_precision="round_trip")
    variants = tuple(old["settings"]["schedules_ASSUMPTION"])
    worlds = tuple(old["settings"]["worlds"])
    summ = summarise_wall(net, agg, win, old["truth_share_below_threshold_daily_min_by_month"], old["seeds"],
                          _wall_settings(net, variants, worlds))
    with open(path, "w") as fh:
        json.dump(_clean(summ), fh, indent=1)
    return summ


def grid_edge_checks(net: str, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache") -> dict:
    """After the review of task 10b (not pre-registered, never used to choose a model): how much of M2's posterior sits
    at the edges of the kb and kw grids.  Refits M2 on every V1 rolling window of the committed seeds in both worlds as
    run_task does (the same log, window and target temperature; the bank without its plus or minus 2 C blocks, which no
    V1 rolling fit reads), checks each refit against the committed row of season2_windows_<net>.csv (every M2 column,
    to the CSV's 6 significant digits), and writes per world the mean share of each wall E/R's mass at the top and the
    bottom of the kw grid and at the bottom of the kb grid, after the January-to-March fit and over all rolling fits,
    with the share of windows where B0's most likely kb is the grid's lowest (from the windows file).
    -> outdir/task10b_grid_edges_<net>.json (no timings in it)."""
    from .experiment import NET_TRUTH
    from .features import build_features
    from .pilot import synthetic_log
    from .simgp import KB_GRID, KW_GRID
    from .simulate import nominal_scenario
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass
    outdir, cache_dir = os.path.abspath(outdir), os.path.abspath(cache_dir)
    with open(os.path.join(outdir, f"summary_season2_{net}.json")) as fh:
        seeds = [int(s) for s in json.load(fh)["seeds"]]
    win = pd.read_csv(os.path.join(outdir, f"season2_windows_{net}.csv"), float_precision="round_trip")
    sched = plant_schedule("V1")
    truth_kw = NET_TRUTH.get(net, {})
    prior_mean = float(np.mean(WALL_ER_HYPOTHESES_K))

    def same6(a, b):
        if pd.isna(a) and pd.isna(b):
            return True
        return float(CSV_FLOAT_10B % float(a)) == float(b)

    rows, n_cmp, bad = [], 0, {}
    workdir = tempfile.mkdtemp(prefix="rm_edges_")
    here = os.getcwd()
    try:
        os.chdir(workdir)                     # EPANET's scratch files stay out of the repo
        sc = nominal_scenario(net, 14, _dose(net))
        X = build_features(sc)
        temps_m, _ = all_bank_temps()
        wb = wall_bank(sc, temps_m, cache_dir=cache_dir, cache="read")
        for world in [w for w in WALL_WORLDS_10B if w in set(win.world)]:
            wall_kw = {} if WALL_WORLDS_10B[world] == "theta_w" else {"truth_wall_law": WALL_WORLDS_10B[world]}
            for seed in seeds:
                log, _, _, _ = synthetic_log(net, months=MONTHS, per_month=PER_MONTH, rotating=ROTATING, seed=seed,
                                             schedule=sched, return_truth=True, **truth_kw, **wall_kw)
                log = log.assign(mi=[int(m[-2:]) for m in log.month])
                for m in ROLLING_TARGETS:
                    s = log[log.mi.isin((m - 3, m - 2, m - 1))][["junction", "hour", "y", "temp_C"]].copy()
                    fit = WallSeasonalSimGP24(sc, X, wb, seed=seed, cache_dir=cache_dir,
                                              threshold=THRESHOLD).fit(s, target_temp_C=sched[m - 1]["temp_C"])
                    vals = {"post_kb": fit.kb20(), "post_kw": fit.rate20(1), "E_post_mean_K": fit.e_posterior_mean(),
                            **fit.posterior_columns()}
                    ref = win[(win.world == world) & (win.variant == "V1") & (win.seed == seed) & (win.month == m)]
                    if len(ref) != 1:
                        raise ValueError(f"no single committed window row for {world} seed {seed} month {m}")
                    for c in [c for c in win.columns if c.startswith("M2_")]:
                        n_cmp += 1
                        if not same6(vals[c[3:]], ref[c].iloc[0]):
                            bad[f"{world}|{seed}|{m}|{c}"] = (vals[c[3:]], float(ref[c].iloc[0]))
                    e = fit.edge_shares()
                    rows.append({"world": world, "seed": seed, "month": m,
                                 "E_wall_post_mean_K": vals["E_wall_post_mean_K"],
                                 "wall_E_eff_K": float(ref.wall_E_eff_K.iloc[0]), "P_H0": vals["P_H0"],
                                 **{f"{k}_W{w:g}": v for k in ("kw_top", "kw_bottom", "kb_bottom")
                                    for w, v in e[k].items()},
                                 "kb_bottom_all_pairs": e["kb_bottom_all_pairs"]})
    finally:
        os.chdir(here)
        shutil.rmtree(workdir, ignore_errors=True)
    df = pd.DataFrame(rows)
    share_cols = [c for c in df.columns if c.startswith(("kw_top_", "kw_bottom_", "kb_bottom_"))]
    out = {"generated_by": f"python -m residualmap.seasonal {net} 0 --wall --grid-edges", "network": net,
           "about": "After the review of task 10b; not pre-registered and never used to choose a model. M2 refitted on "
                    "every V1 rolling window of the scored seeds in both worlds; shares are of each wall E/R's posterior "
                    "mass (H0 excluded) on members at the grid's edge, read as 20 C rates: kw_top kw20 = "
                    f"{max(KW_GRID):g} m/day, kw_bottom kw20 = {min(KW_GRID):g} m/day, kb_bottom kb20 = "
                    f"{min(KB_GRID):g} per day; kb_bottom_all_pairs over all 15 pairs. Simulated.",
           "seeds": seeds,
           "reproduces_committed_windows_rows": {"values_compared": n_cmp, "identical_to_6_significant_digits":
                                                 n_cmp - len(bad), "mismatches": dict(list(bad.items())[:20])},
           "by_world": {}}
    for world, g in df.groupby("world", sort=True):
        node = {}
        for name, gg in (("jan_mar_fit", g[g.month == 4]), ("all_rolling", g)):
            post, eff = gg.E_wall_post_mean_K.values, gg.wall_E_eff_K.values
            node[name] = {"n_fits": int(len(gg)), "mean_shares": {c: float(gg[c].mean()) for c in share_cols},
                          "mean_posterior_wall_E_K": float(post.mean()), "mean_truth_wall_E_eff_K": float(eff.mean()),
                          "share_closer_than_prior": float((np.abs(post - eff) < np.abs(prior_mean - eff)).mean())}
        node["jan_mar_fit_per_seed"] = g[g.month == 4].drop(columns=["world", "month"]).to_dict("records")
        w = win[win.world == world]
        node["B0_map_kb_at_floor_share_of_windows"] = {v: float((gv.B0_map_kb == min(KB_GRID)).mean())
                                                       for v, gv in w.groupby("variant", sort=True)}
        node["M2_post_kb20_mean_V1"] = float(w[w.variant == "V1"].M2_post_kb.mean())
        out["by_world"][world] = node
    out = _clean(_round_sig(out))
    with open(os.path.join(outdir, f"task10b_grid_edges_{net}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ---- task 10b figures (light surface; categorical slots validated with the dataviz validator, adjacent pairs)
FIG_COLORS_10B = {"B0": "#8a8984", "M": "#eb6834", "M1b": "#1baf7a", "M2": "#2a78d6", "oracle": "#fcfcfb"}
WORLD_COLORS_10B = {"W1": "#eb6834", "W2": "#2a78d6"}
CLASS_TITLES_10B = {("rolling", "warming"): "rolling, warming months (Apr to Aug)",
                    ("rolling", "cooling"): "rolling, cooling months (Sep to Dec)",
                    ("jan_mar_to_jul_sep", "all_months"): "fit Jan to Mar, predict Jul to Sep",
                    ("jul_sep_to_dec", "all_months"): "fit Jul to Sep, predict Dec"}


def plot_season2(net: str, summ: dict, out: str) -> None:
    """Daily-minimum map RMSE (top) and recall (bottom) of every model, V1, in the four A7 classes, both worlds."""
    plt = _plt()
    from .experiment import BIG
    P = summ["pooled"]
    worlds = [w for w in ("W1", "W2") if w in P]
    models = list(MODELS_10B)
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(2, 4, figsize=(25, 10.5))
        for j, (test, cls) in enumerate(A7_CLASSES_10B):
            for i, (metric, lab) in enumerate((("rmse", "map RMSE (mg/L)"), ("recall", "map recall"))):
                a = ax[i, j]
                width = 0.16
                for k, model in enumerate(models):
                    xs, ys = [], []
                    for wi, w in enumerate(worlds):
                        d = P[w]["V1"][test]["none"][cls].get(model, {}).get("daily_min_map")
                        if d is None or d.get(metric) is None:
                            continue
                        xs.append(wi + (k - 2) * width)
                        ys.append(d[metric])
                    edge = "#0b0b0b" if model == "oracle" else FIG_COLORS_10B[model]
                    a.bar(xs, ys, width * 0.88, color=FIG_COLORS_10B[model], edgecolor=edge, linewidth=1.0,
                          hatch="//" if model == "oracle" else None, label=model if (i, j) == (0, 0) else None)
                    if model in ("B0", "M2"):
                        for x, y in zip(xs, ys):
                            a.text(x, y, f"{y:.3f}" if metric == "rmse" else f"{y:.2f}", ha="center", va="bottom",
                                   fontsize=10, color="#0b0b0b")
                if not any(P[w]["V1"][test]["none"][cls].get(mm_, {}).get("daily_min_map", {}).get(metric) is not None
                           for w in worlds for mm_ in models):
                    a.text(0.5, 0.5, "no low junctions in this class", transform=a.transAxes, ha="center", color="#52514e")
                a.set_xticks(range(len(worlds)))
                a.set_xticklabels([f"{w}: {'weak wall (task 10)' if w == 'W1' else 'wall follows bulk'}" for w in worlds])
                a.set_ylabel(lab)
                if metric == "recall":
                    a.set_ylim(0, 1.08)
                a.grid(axis="y", alpha=0.25)
                a.spines[["top", "right"]].set_visible(False)
                if i == 0:
                    a.set_title(CLASS_TITLES_10B[(test, cls)])
        fig.legend(*ax[0, 0].get_legend_handles_labels(), loc="upper center", ncol=5, frameon=False,
                   bbox_to_anchor=(0.5, 0.935))
        stops = ", ".join(f"{w}: {'triggered' if summ['acceptance'][w]['A7']['stop'] else 'not triggered'}" for w in worlds)
        ns = len(summ["seeds"])
        fig.suptitle(f"{net}, simulated 12-month logs, {ns} fresh seed{'s' if ns != 1 else ''}, seasonal truth (V1): "
                     f"the daily-minimum map by model (stop rule A7 for M2, every variant and class: {stops})\n"
                     f"B0 today's model; M wall E/R = bulk E/R (task 10); M1b no wall response; M2 the wall's response "
                     f"learned; oracle the best single member of M2's bank (hatched)", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.91))
        fig.savefig(out, dpi=100)
        plt.close(fig)


def plot_wall_posterior(net: str, win: pd.DataFrame, summ: dict, out: str) -> None:
    """R2: M2's posterior over the wall E/R in the two worlds (V1), against the prior and the truth's own response."""
    plt = _plt()
    from .experiment import BIG
    r2 = summ["r2_wall_posterior"]["by_world"]
    worlds = [w for w in ("W1", "W2") if w in r2]
    labels = [f"{w:g}" for w in WALL_ER_HYPOTHESES_K]
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(1, 2, figsize=(20, 7))
        a = ax[0]
        width = 0.2
        sets = (("jan_mar_window", "Jan-to-Mar window (0.5 C of contrast)", 0.45), ("windows_spanning_4C_or_more",
                                                                                     "windows spanning 4 C or more", 1.0))
        k = 0
        for w in worlds:
            for name, lab, alpha in sets:
                m = r2[w]["V1"][name]["mean_mass_within_pairs"]
                xs = np.arange(len(labels)) + (k - 1.5) * width
                a.bar(xs, [m[f"W{x}"] for x in labels], width * 0.88, color=WORLD_COLORS_10B[w], alpha=alpha,
                      label=f"{w}, {lab}")
                k += 1
        a.axhline(1.0 / len(labels), color="#52514e", ls="--", lw=1.5, label="prior (each wall E/R)")
        a.set_ylim(0, max(1.0 / len(labels), max(b.get_height() for b in a.patches)) * 1.5)
        a.set_xticks(range(len(labels)))
        a.set_xticklabels([f"{x} K\n(x{1.0 / arrhenius(10.0, float(x)):.2f} per 10 C)" for x in labels])
        a.set(ylabel="posterior mass within the pairs (H0 excluded)", xlabel="wall E/R hypothesis",
              title="M2's posterior over the wall's temperature response (V1, mean over fits)")
        a.grid(axis="y", alpha=0.25)
        a.spines[["top", "right"]].set_visible(False)
        a.legend(fontsize=10, frameon=False)
        a = ax[1]
        r = win[(win.variant == "V1") & np.array([_window_span(int(m)) >= 4.0 for m in win.month])]
        g = r.groupby(["world", "seed"])[["M2_E_wall_post_mean_K", "wall_E_eff_K"]].mean().reset_index()
        for w in worlds:
            gw = g[g.world == w]
            a.scatter(gw.wall_E_eff_K, gw.M2_E_wall_post_mean_K, s=70, color=WORLD_COLORS_10B[w], edgecolor="#fcfcfb",
                      linewidth=1.5, label=f"{w}: one seed (mean over the windows spanning 4 C or more)", zorder=3)
        lim = (-300, 12600)
        a.plot(lim, lim, color="#0b0b0b", lw=1.2, label="posterior = truth")
        a.axhline(float(np.mean(WALL_ER_HYPOTHESES_K)), color="#52514e", ls="--", lw=1.5, label="prior mean, 5500 K")
        a.set(xlim=lim, ylim=lim, xlabel="truth's own wall response as a wall E/R (K)",
              ylabel="M2's posterior mean wall E/R (K)", title="Does the posterior move toward the truth's wall? (V1)")
        a.grid(alpha=0.25)
        a.spines[["top", "right"]].set_visible(False)
        a.legend(fontsize=10, frameon=False, loc="upper left")
        ns = len(summ["seeds"])
        fig.suptitle(f"{net}, simulated, {ns} fresh seed{'s' if ns != 1 else ''}: the wall's temperature response, learned "
                     f"(report only, no bar). W1: weak wall (task 10's truth, an assumption); W2: wall follows the bulk",
                     fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(out, dpi=100)
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("net")
    ap.add_argument("seeds", type=int)
    ap.add_argument("--outdir", default="outputs/chem")
    ap.add_argument("--cache", default="outputs/cache")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--mechanism", action="store_true",
                    help="after the stop: write outdir/task10_mechanism_checks.json (mechanism_checks) instead of the run")
    ap.add_argument("--wall", action="store_true",
                    help="task 10b: M2 (bulk E/R x wall E/R) in worlds W1 and W2 on fresh seeds (run_wall); seeds = "
                         "first-seed .. first-seed + seeds - 1")
    ap.add_argument("--first-seed", type=int, default=None,
                    help="task 10b: first seed (default 16 on Net3, 8 on Net2; task 10's seeds are refused)")
    ap.add_argument("--worlds", default=",".join(WALL_WORLDS_10B), help="task 10b: truth worlds, W1 and/or W2")
    ap.add_argument("--resummarise", action="store_true",
                    help="with --wall, after the review: rewrite outdir/summary_season2_<net>.json from the committed "
                         "CSVs (resummarise_wall; no rerun; seeds is ignored)")
    ap.add_argument("--grid-edges", action="store_true",
                    help="with --wall, after the review: write outdir/task10b_grid_edges_<net>.json (grid_edge_checks; "
                         "refits M2 on the committed seeds' V1 rolling windows; seeds is ignored)")
    a = ap.parse_args(argv)
    if a.resummarise or a.grid_edges:
        if not a.wall:
            ap.error("--resummarise and --grid-edges go with --wall")
        warnings.filterwarnings("ignore")
        if a.resummarise:
            resummarise_wall(a.net, a.outdir)
        if a.grid_edges:
            grid_edge_checks(a.net, a.outdir, a.cache)
        return 0
    if a.wall:
        warnings.filterwarnings("ignore")
        variants = tuple(v for v in a.variants.split(",") if v)
        if "V1" not in variants or "V2" not in variants:
            ap.error("V1 and V2 are needed for the acceptance key")
        first = a.first_seed if a.first_seed is not None else FIRST_SEED_10B.get(a.net, 0)
        worlds = tuple(w for w in a.worlds.split(",") if w)
        s = run_wall(a.net, tuple(range(first, first + a.seeds)), a.outdir, a.cache, a.workers, variants, worlds)
        return 2 if any(x["A7"]["stop"] for x in s["acceptance"].values()) else 0
    if a.mechanism:
        warnings.filterwarnings("ignore")
        mechanism_checks(a.net, a.seeds, a.outdir, a.cache, a.workers)
        return 0
    variants = tuple(v for v in a.variants.split(",") if v)
    if "V1" not in variants or "V2" not in variants:
        ap.error("V1 and V2 are needed for the acceptance key")
    warnings.filterwarnings("ignore")
    s = run(a.net, tuple(range(a.seeds)), a.outdir, a.cache, a.workers, variants)
    return 2 if s["acceptance"]["A7"]["stop"] else 0


if __name__ == "__main__":
    sys.exit(main())
