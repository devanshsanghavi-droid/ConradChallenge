"""
simgp.py — calibrated-simulator GP ("the EPANET model is the prior").

The strongest physics we have is the operator's own EPANET model.  What it lacks is the decay
chemistry: bulk decay kb, wall decay kw, and how much old (rough) pipe accelerates wall decay
(gamma).  Those three numbers are calibrated from grab samples:

    1. simulate chlorine on the NOMINAL model over a grid of (kb, kw, gamma)   [cached]
    2. weight each grid member by how well it explains the samples (Bayesian model averaging)
    3. mean function  m(node) = posterior-weighted mean of ln C_sim(node)
       prior variance v(node) = posterior-weighted variance  (parameter uncertainty -> map uncertainty)
    4. a GP on the residual ln y - m learns what the simulator still gets wrong
       (demand mismatch, pipe-level heterogeneity, dose drift) from the same features as PhysicsGP
    5. predictive:  z ~ N(m + r_mu, v + r_sd^2)

This is a Kennedy-O'Hagan calibration-plus-discrepancy model with a grid posterior instead of MCMC.
It uses every physical fact in the .inp — mixing at junctions, tank turnover, pipe lengths,
diameters, roughness, demand patterns — because EPANET does.

Why this instead of a PINN: with 3-15 samples a neural network has nothing to learn from, and the
transport PDE is already solved exactly by EPANET.  A PINN would re-learn a physics we can compute.
See pinn.py for the graph-PINN baseline that tests this claim rather than asserting it.
"""
from __future__ import annotations

import itertools
import os
import pickle
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.special import logsumexp
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from .chemistry import DEFAULT_THRESHOLD_MGL, Chemistry, cache_tag, remove_epanet_files
from .features import CORE
from .simulate import simulate_nominal_chlorine

FLOOR = 0.02
KB_GRID = [0.10, 0.25, 0.40, 0.55, 0.70]        # 1/day
KW_GRID = [0.10, 0.30, 0.60, 1.00, 1.50]        # m/day
GAMMA_GRID = [0.0, 0.5, 1.0]                    # old-pipe sensitivity
DEMAND_GRID = [0.85, 1.0, 1.15]                 # global demand multiplier (hydraulic mismatch)
ROUGH_GRID = [0.9, 1.0, 1.1]                    # global Hazen-Williams C multiplier (hydraulic mismatch)
# chloramine (task 12): total chlorine, first order on its own ranges.  KB_CA spans the batch port's apparent rates
# (0.012 to 0.16 per day over pH 7.5 to 8.5 and Cl2:N 4 to 5; up to 0.22 at the low-pH stress pH 7.0 to 7.5) with
# margin at both ends (0.32 was added to the plan's 0.0025 to 0.16, per the plan's addendum 3); the compiled literature
# range is 0.004 to 0.16 per day (Odimayomi et al. 2026, secondary).  KW_CA is an ASSUMPTION, a log span from plastic to
# old iron.  GAMMA_CA is widened upward (wall rates at least 4x bulk in iron and cement pipe).
KB_CA = [0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32]   # 1/day at 20 C
KW_CA = [0.01, 0.03, 0.10, 0.30, 1.0]                          # m/day
GAMMA_CA = [0.0, 0.5, 1.0, 2.0]
# task 13, an OPT-IN variant (plan addendum 3): today's grid with four bulk rates below its 0.10 per day floor, which
# binds in most fits under task 11's second-order truths and on Net2.  0.0125 lies inside, not below, the two-reactant
# truth's local decay rate in Net2's oldest water (0.0138, 0.0126 and 0.0115 per day over 96 to 120, 120 to 144 and
# 144 to 168 h at the 24 h match, outputs/chem_2ra/scaling.json), so the new floor does not bracket it.  Its own
# cache name (grid_cache_path adds a tag); the default 'full' grid, its cache file and every committed number are
# untouched.  Adopting it into the default path would move committed numbers: that is Devansh's decision, not the
# task's.
KB_LOWKB = [0.0125, 0.025, 0.05, 0.075] + KB_GRID   # 1/day
GRIDS = {"decay": (KB_GRID, KW_GRID, GAMMA_GRID, [1.0], [1.0]),                    # iteration 2: 75 runs
         "full": (KB_GRID, KW_GRID, GAMMA_GRID, DEMAND_GRID, ROUGH_GRID),         # iteration 3: 675 runs
         "chloramine": (KB_CA, KW_CA, GAMMA_CA, DEMAND_GRID, ROUGH_GRID),         # task 12: 1440 runs, chloramine only
         "full_lowkb": (KB_LOWKB, KW_GRID, GAMMA_GRID, DEMAND_GRID, ROUGH_GRID)}   # task 13 opt-in: 1215 runs
COMMITTED_GRID_NAMES = ("decay", "full", "chloramine")   # cache names without a grid tag (kept as committed)
DOSE_GRID = [0.90, 0.95, 1.00, 1.05, 1.10]     # source-dose multiplier: first-order decay is linear in
                                                # concentration, so this axis is an exact ln-offset, no runs
# chloramine's effective-dose axis also absorbs the fast organic demand (5.92 S1 TOC mg/L, Duirk et al. 2005; about
# 0.24 mg/L at TOC 2, d about 0.88 at a 2.0 mg/L dose), so it mixes dose error with that demand.  0.65 and 0.70 were
# added below the plan's 0.75 for margin (TOC 3 and a dose draw of 0.9 give about 0.74).
DOSES_CA = [0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10]
DOSES = {"free_chlorine": DOSE_GRID, "chloramine": DOSES_CA}


HOURS = list(range(24))
DAY_HOURS = list(range(7, 18))                  # 07:00-17:00: when an operator can take a grab sample
STATIC = ["emb0", "emb1", "emb2", "path_wall_index", "dist_src_km"]   # CORE minus the hour-dependent age


MIN_FREE_GB = 1.5   # a grid build is refused below this much free disk (an EPANET run that fills its temp dir
                    # fails with error 308; a fresh ky4 grid needs about 6 GB of temp files without cleanup)


def disk_free_gb(path: str) -> float:
    """Free space, in GB (1e9 bytes), on the volume holding `path` (or its nearest existing parent)."""
    p = os.path.abspath(path)
    while not os.path.exists(p):
        p = os.path.dirname(p)
    return shutil.disk_usage(p).free / 1e9


def disk_preflight(paths=None, min_gb: float = MIN_FREE_GB) -> float:
    """Refuse to start a grid build with less than `min_gb` GB free on the cache or temp volume.  Returns the
    smallest free space found (GB)."""
    paths = list(paths) if paths else [tempfile.gettempdir()]
    free = min(disk_free_gb(p) for p in paths)
    if free < min_gb:
        raise RuntimeError(f"refusing to build an EPANET grid: {free:.2f} GB free on the cache or temp volume, "
                           f"below the {min_gb:g} GB floor.  Free some disk space first.")
    return free


def _grid_member(args):
    """One EPANET run of the grid (module-level so it can run in a worker process).  Its .inp, .rpt and .bin
    are deleted as soon as the result is read (about 0.9 MB per Net3 run and 9 MB per ky4 run used to stay on
    disk until the pool closed).  An optional 10th element holds extra simulate_nominal_chlorine keywords
    (a chemistry condition's kb_scale, kw_scale, temp_C)."""
    name, kb, kw, g, dose, dm, rm, junctions, prefix, *extra = args
    try:
        c = simulate_nominal_chlorine(name, kb, kw, g, dose, dm, rm, file_prefix=prefix, **(extra[0] if extra else {}))
    finally:
        remove_epanet_files(prefix)
    return np.log(np.clip(c.loc[HOURS, junctions].values, FLOOR, None)).astype(np.float32)


def _non_default(cond: Chemistry | None) -> bool:
    return cond is not None and not cond.is_default


def grid_cache_path(sc, cache_dir: str = "outputs/cache", grid: str = "full", cond: Chemistry | None = None) -> str:
    """Cache file of a grid.  Default calls (cond None, or a condition that reproduces today's path) keep the
    committed names, so existing pickles stay valid; any other condition appends its label and
    chemistry.cache_tag (a hash of the grid definition and the condition).  A grid that is not one of
    COMMITTED_GRID_NAMES (task 13's opt-in 'full_lowkb') always carries chemistry.cache_tag, which hashes its axes,
    so a change to its axes can never be served an old file."""
    dose_tag = "" if abs(sc.source_dose - 1.2) < 1e-9 else f"_d{sc.source_dose:g}"
    cond_tag = f"_{cond.label()}_{cache_tag(cond, grid)}" if _non_default(cond) else ""
    if not cond_tag and grid not in COMMITTED_GRID_NAMES:
        cond_tag = f"_{cache_tag(None, grid)}"
    return os.path.join(cache_dir, f"grid24_{grid}_{os.path.basename(sc.wn_name)}{dose_tag}{cond_tag}.pkl")


def check_grid_disinfectant(grid: str, cond: Chemistry | None) -> None:
    """Free chlorine and chloramine never share a grid: GRIDS['chloramine'] runs only under a chloramine condition,
    and a chloramine condition only on GRIDS['chloramine'].  (A free-chlorine grid may still be scored against a
    chloramine truth: that is task 12's conflation test, run with cond=None.)"""
    ca = cond is not None and cond.disinfectant == "chloramine"
    if (grid == "chloramine") != ca:
        raise ValueError(f"grid {grid!r} with disinfectant {cond.disinfectant if cond is not None else 'free_chlorine'!r}: "
                         "the chloramine grid needs Chemistry(disinfectant='chloramine') and that condition needs the "
                         "chloramine grid")


def simulator_grid_24h(sc, cache_dir: str = "outputs/cache", grid: str = "full",
                       n_jobs: int | None = None, cond: Chemistry | None = None) -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C for every hour of the last day: (params, array [members, 24, junctions]).
    params are (kb, kw, gamma, demand_mult, rough_mult).

    Each EPANET run already produces the whole day; storing all of it is what lets the time-aware
    model (SimGP24) calibrate on daytime samples and predict the night.  grid="full" adds the two
    hydraulic-mismatch axes (675 runs); grid="decay" is the iteration-2 grid (75 runs).
    cond: a chemistry.Chemistry; params stay the 20 C values, the runs use cond's rate multipliers and water
    properties.  Cached under grid_cache_path; a build (not a cache hit) is refused below MIN_FREE_GB free disk.
    A condition the grid cannot run (a seasonal or TOC-scaled chloramine grid, kinetics other than 'first', a
    disinfectant on the other disinfectant's grid) is refused before the cache is read, so a file under its name is
    never served."""
    check_grid_disinfectant(grid, cond)
    if _non_default(cond):
        cond.sim_kwargs()      # raises for a condition that is not built
    os.makedirs(cache_dir, exist_ok=True)
    f = grid_cache_path(sc, cache_dir, grid, cond)
    if os.path.exists(f):
        with open(f, "rb") as fh:
            return pickle.load(fh)
    out = build_grid_24h(sc, grid, n_jobs, cond, preflight_paths=[cache_dir])
    with open(f, "wb") as fh:
        pickle.dump(out, fh)
    return out


def build_grid_24h(sc, grid: str = "full", n_jobs: int | None = None, cond: Chemistry | None = None,
                   preflight_paths=None) -> tuple[list[tuple], np.ndarray]:
    """The grid itself, without the cache (simulator_grid_24h caches it).  Every run's EPANET files are deleted
    as soon as it is read, and the whole build is refused below MIN_FREE_GB free disk."""
    check_grid_disinfectant(grid, cond)
    extra = (cond.sim_kwargs(),) if _non_default(cond) else ()
    disk_preflight(list(preflight_paths or []) + [tempfile.gettempdir()])
    params = list(itertools.product(*GRIDS[grid]))
    n_jobs = n_jobs or max(1, (os.cpu_count() or 2) - 1)
    with tempfile.TemporaryDirectory() as tmp:
        # EPANET writes temp.inp/.rpt/.bin per run: give every run its own prefix so runs can go in parallel
        jobs = [(sc.wn_name, kb, kw, g, sc.source_dose, dm, rm, list(sc.junctions), os.path.join(tmp, f"g{i}")) + extra
                for i, (kb, kw, g, dm, rm) in enumerate(params)]
        if n_jobs > 1 and len(jobs) > 8:
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                rows = list(ex.map(_grid_member, jobs, chunksize=4))
        else:
            rows = [_grid_member(j) for j in jobs]
    return params, np.stack(rows)


def simulator_grid(sc, cache_dir: str = "outputs/cache", grid: str = "full") -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C at the sampling hour: (params list, array [n_members, n_junctions])."""
    params, Z = simulator_grid_24h(sc, cache_dir, grid)
    return params, Z[:, sc.sample_hour, :]


def grid_loglik(z_sim: np.ndarray, z_obs: np.ndarray, lik_sd: float, lik: str = "gauss",
                nu: float = 3.0) -> np.ndarray:
    """Unnormalised log-likelihood of each grid member given log-samples.  z_sim: members x samples.

    lik="gauss" is iteration 2.  lik="t" (Student-t, nu degrees of freedom) is robust: a junction where
    the operator's model is structurally wrong (a front that sits elsewhere, a tank zone that turns over
    at a different hour) can be 3x off for EVERY grid member; under a Gaussian one such reading drags
    the whole calibration, under a Student-t it is discounted."""
    scales = np.atleast_1d(np.asarray(lik_sd, dtype=float))
    lls = []
    for sd in scales:
        e = (z_sim - z_obs[None, :]) / sd
        if lik == "gauss":
            ll = -0.5 * (e ** 2).sum(axis=1) - z_sim.shape[1] * np.log(sd)
        elif lik == "t":
            ll = (-(nu + 1) / 2 * np.log1p(e ** 2 / nu)).sum(axis=1) - z_sim.shape[1] * np.log(sd)
        else:
            raise ValueError(lik)
        lls.append(ll)
    # several scales = the model-error scale is unknown too: marginalise it (uniform prior over the list)
    return logsumexp(np.vstack(lls), axis=0)


def grid_weights(z_sim, z_obs, lik_sd, lik="gauss", nu=3.0) -> np.ndarray:
    ll = grid_loglik(z_sim, z_obs, lik_sd, lik, nu)
    return np.exp(ll - logsumexp(ll))


def grid_dose_weights(z_sim: np.ndarray, z_obs: np.ndarray, lik_sd: float, lik: str, nu: float,
                      doses=DOSE_GRID, log_prior: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Joint posterior over (grid member, dose multiplier): W [members x doses], and the ln-offsets.
    A dose multiplier m shifts every simulated ln C by ln m, so member k at dose d predicts Z_k + offs_d,
    which is scored against z_obs as Z_k against (z_obs - offs_d).
    log_prior: an optional log prior over the members (task 12's model (c)); None is the committed uniform prior and
    the committed arithmetic."""
    offs = np.log(np.asarray(doses, dtype=float))
    ll = np.stack([grid_loglik(z_sim, z_obs - d, lik_sd, lik, nu) for d in offs], axis=1)
    if log_prior is not None:
        ll = ll + np.asarray(log_prior, dtype=float)[:, None]
    return np.exp(ll - logsumexp(ll)), offs


def posterior_moments(W: np.ndarray, offs: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean and variance of ln C under the joint (member, dose) posterior W, for Z of shape
    (members, ...) — per junction for SimGP, per (hour, junction) for SimGP24.
        E[Z + d]       = sum_k w_k Z_k + sum_d wd_d offs_d
        E[(Z + d)^2]   = sum_k w_k Z_k^2 + 2 sum_k (sum_d W_kd offs_d) Z_k + sum_d wd_d offs_d^2
    with w, wd the marginals.  Variance floored at 1e-6."""
    w, wd = W.sum(axis=1), W.sum(axis=0)
    wo = W @ offs
    m = np.tensordot(w, Z, axes=1) + wd @ offs
    v = np.tensordot(w, Z ** 2, axes=1) + 2 * np.tensordot(wo, Z, axes=1) + wd @ offs ** 2 - m ** 2
    return m, np.clip(v, 1e-6, None)


LIK_SD = 0.35   # log-space scale of the calibration likelihood = the day-time RMS mismatch between the truth
                # and the best grid member on Net3 (0.33), i.e. model error, not grab-sample noise (0.03 mg/L)
LIK_SD_CA = 0.10   # task 12, the same rule on the chloramine truth: the daytime RMS ln mismatch of the best member of
                   # GRIDS['chloramine'] on Net3 calibration seeds 100 to 103 (never scored) is 0.0955 on average, rounded
                   # up to a multiple of 0.05 (outputs/chloramine/calibration_chloramine.json; Net2's would be 0.20)
LIK_SD_BY = {"free_chlorine": LIK_SD, "chloramine": LIK_SD_CA}


def n_hydraulic(grid: str) -> int:
    return len(GRIDS[grid][3]) * len(GRIDS[grid][4])


def check_grid_order(params: list[tuple], n_hyd: int) -> None:
    """The reshapes in local_hydraulic_var and SimGP24.predict_daily_min assume itertools.product order:
    consecutive blocks of n_hyd members share one decay triple and run through the hydraulic pairs."""
    P = np.asarray(params, dtype=float).reshape(-1, n_hyd, 5)
    if not ((P[:, :, :3] == P[:, :1, :3]).all() and (P[:, :, 3:] == P[:1, :, 3:]).all()):
        raise ValueError("grid members are not in (decay-major, hydraulic-minor) order; GRIDS changed?")


PARAM_NAMES = ("kb", "kw", "gamma", "demand", "rough")


def grid_edge_mass(W: np.ndarray, params: list[tuple], doses=None) -> dict:
    """Posterior mass on the edges of each grid axis: for every parameter, the weight of members at its
    smallest value (`<name>_low`), at its largest (`<name>_high`) and their sum (`<name>`).  A large edge mass
    says the data want a value outside the grid.  W is the joint (member x dose) posterior or the member
    marginal; with the joint, `doses` must be the model's own dose multipliers (model.doses, in W's column
    order) and the dose axis is reported too (its smallest and largest multiplier).  Axes with a single value
    give NaN."""
    W = np.asarray(W, dtype=float)
    if W.ndim not in (1, 2) or W.shape[0] != len(params):
        raise ValueError(f"W has shape {W.shape}; expected ({len(params)},) or ({len(params)}, n_doses)")
    if W.ndim == 2:
        if doses is None:
            raise ValueError("a joint (member x dose) posterior needs doses=model.doses")
        if len(doses) != W.shape[1]:
            raise ValueError(f"{len(doses)} doses for a posterior with {W.shape[1]} dose columns")
    w = W.sum(axis=1) if W.ndim == 2 else W
    P = np.asarray(params, dtype=float)
    out = {}
    for i, name in enumerate(PARAM_NAMES):
        vals = np.unique(P[:, i])
        if len(vals) < 2:
            out[f"{name}_low"] = out[f"{name}_high"] = out[name] = float("nan")
            continue
        lo, hi = float(w[P[:, i] == vals[0]].sum()), float(w[P[:, i] == vals[-1]].sum())
        out[f"{name}_low"], out[f"{name}_high"], out[name] = lo, hi, lo + hi
    if W.ndim == 2:
        d = np.asarray(doses, dtype=float)
        wd = W.sum(axis=0)
        if len(d) < 2:
            out["dose_low"] = out["dose_high"] = out["dose"] = float("nan")
        else:
            lo, hi = float(wd[d == d.min()].sum()), float(wd[d == d.max()].sum())
            out["dose_low"], out["dose_high"], out["dose"] = lo, hi, lo + hi
    return out


def local_hydraulic_var(w: np.ndarray, Z: np.ndarray, n_hyd: int) -> np.ndarray:
    """Variance of the simulator across the hydraulic axes, averaged over the decay-parameter posterior.

    Calibrating the GLOBAL demand and roughness multipliers does not remove the operator's LOCAL errors
    (per-node demand, per-pipe roughness), which are of the same size as the axes.  The grid's own
    spread along those axes at each (hour, junction) is used as the variance of that local error."""
    nd = len(w) // n_hyd
    Zr = Z.reshape(nd, n_hyd, *Z.shape[1:])
    wd = w.reshape(nd, n_hyd).sum(axis=1)
    return np.tensordot(wd, Zr.var(axis=1), axes=1)


class SimGP:
    def __init__(self, sc, seed: int = 0, columns: list[str] | None = None,
                 lik_sd: float = LIK_SD, cache_dir: str = "outputs/cache", lik: str = "t", nu: float = 3.0,
                 grid: str = "full", local_hydraulic: bool = True, doses=DOSE_GRID):
        self.sc = sc
        self.seed = seed
        self.columns = columns or CORE
        self.lik_sd, self.lik, self.nu = lik_sd, lik, nu   # log-space tolerance / family when scoring grid members
        self.doses = doses
        self.params, self.Z = simulator_grid(sc, cache_dir, grid)   # Z: members x junctions
        self.n_hyd = n_hydraulic(grid) if local_hydraulic else 0
        if self.n_hyd > 1:
            check_grid_order(self.params, self.n_hyd)
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}

    def _calibrate(self, nodes: list[str], y: np.ndarray) -> None:
        idx = [self.jidx[n] for n in nodes]
        z_obs = np.log(np.clip(y, FLOOR, None))
        W, offs = grid_dose_weights(self.Z[:, idx], z_obs, self.lik_sd, self.lik, self.nu, self.doses)
        w = W.sum(axis=1)                                 # marginal over members
        self.w_, self.W_, self.offs_ = w, W, offs
        self.m_, self.v_ = posterior_moments(W, offs, self.Z)          # per junction
        self.hv_ = local_hydraulic_var(w, self.Z, self.n_hyd) if self.n_hyd > 1 else np.zeros_like(self.v_)
        k, d = np.unravel_index(int(np.argmax(W)), W.shape)
        self.map_params_ = self.params[k]
        self.map_dose_ = float(np.exp(offs[d]))

    def fit(self, X: pd.DataFrame, y: np.ndarray) -> "SimGP":
        nodes = list(X.index)
        self._calibrate(nodes, y)
        idx = [self.jidx[n] for n in nodes]
        r = np.log(np.clip(y, FLOOR, None)) - self.m_[idx]
        Xc = X[self.columns]
        self.mu_, self.sd_ = Xc.mean(), Xc.std().replace(0, 1.0)
        Xs = ((Xc - self.mu_) / self.sd_).values
        k = (ConstantKernel(0.1, (1e-3, 5.0))
             * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=(0.1, 20.0), nu=1.5)
             + WhiteKernel(0.01, (1e-4, 0.5)))
        self.gp = GaussianProcessRegressor(kernel=k, n_restarts_optimizer=4, random_state=self.seed)
        self.gp.fit(Xs, r)
        return self

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        idx = [self.jidx[n] for n in X.index]
        Xc = X[self.columns]
        Xs = ((Xc - self.mu_) / self.sd_).values
        r_mu, r_sd = self.gp.predict(Xs, return_std=True)
        noise = self.gp.kernel_.k2.noise_level
        r_var = np.clip(r_sd ** 2 - noise, 1e-4, None)
        z_mu = self.m_[idx] + r_mu
        z_sd = np.sqrt(self.v_[idx] + self.hv_[idx] + r_var)
        out = pd.DataFrame(index=X.index)
        out["median"] = np.exp(z_mu)
        out["lo90"] = np.exp(z_mu - 1.645 * z_sd)
        out["hi90"] = np.exp(z_mu + 1.645 * z_sd)
        out["z_mu"], out["z_sd"] = z_mu, z_sd
        # what a grab sample can still reduce: parameter spread + discrepancy GP, not the local hydraulic
        # error (which stays whatever we sample) — acquisition rules use this one
        out["z_sd_acq"] = np.sqrt(self.v_[idx] + r_var)
        return out

    def predict_prior(self, X: pd.DataFrame) -> pd.DataFrame:
        """The calibrated simulator alone (no discrepancy GP): what the grid posterior says by itself."""
        idx = [self.jidx[n] for n in X.index]
        z_mu, z_sd = self.m_[idx], np.sqrt(self.v_[idx] + self.hv_[idx])
        out = pd.DataFrame(index=X.index)
        out["median"] = np.exp(z_mu)
        out["lo90"] = np.exp(z_mu - 1.645 * z_sd)
        out["hi90"] = np.exp(z_mu + 1.645 * z_sd)
        out["z_mu"], out["z_sd"] = z_mu, z_sd
        return out

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)


# ----------------------------------------------------------------------------- time-aware model
class SimGP24:
    """Time-aware calibrated-simulator GP.  Samples are (junction, hour, mg/L); predictions are the full
    24-h profile at every junction and, from it, the daily minimum — the compliance number.

    Same three steps as SimGP, in time:
      1. grid members are scored on the simulated value at each sample's OWN hour;
      2. the discrepancy GP sees (age at that hour, hydraulic embedding, wall index, distance to source,
         sin/cos of hour) so it can be told a 09:00 and a 16:00 reading at the same junction differ;
      3. the daily minimum is a Monte-Carlo draw: pick a grid member by its weight (whole-day profile),
         add a joint GP draw over the 24 hours of that junction, take the minimum.
    Away from the sampled hours the GP reverts to the simulator, with a wider band — so night predictions
    rest on the calibrated physics, and the band says so.
    """

    def __init__(self, sc, X: pd.DataFrame, seed: int = 0, lik_sd: float | None = None,
                 cache_dir: str = "outputs/cache", n_draws: int = 1024, lik: str = "t", nu: float = 3.0,
                 grid: str = "full", local_hydraulic: bool = True, doses=None, smooth_hours: bool = True,
                 threshold: float | None = None, cond: Chemistry | None = None, log_prior: np.ndarray | None = None):
        # lik_sd, doses and threshold left as None follow the condition's disinfectant (task 12 review): free chlorine
        # (cond None) gets the committed LIK_SD, DOSE_GRID and 0.2, exactly as before; a chloramine condition gets
        # LIK_SD_CA, DOSES_CA and 0.5 (chemistry.DEFAULT_THRESHOLD_MGL), so the chloramine grid never runs on free
        # chlorine's settings by omission
        dis = cond.disinfectant if cond is not None else "free_chlorine"
        lik_sd = LIK_SD_BY[dis] if lik_sd is None else lik_sd
        doses = DOSES[dis] if doses is None else doses
        threshold = DEFAULT_THRESHOLD_MGL[dis] if threshold is None else threshold
        self.sc, self.seed, self.lik_sd, self.n_draws = sc, seed, lik_sd, n_draws
        self.lik, self.nu, self.doses, self.smooth_hours = lik, nu, doses, smooth_hours
        self.threshold = float(threshold)   # compliance threshold of predict_daily_min's p_below column (mg/L)
        # task 12: the grid's chemistry condition (None: the committed free-chlorine grids; a chloramine condition with
        # grid='chloramine') and an optional log prior over the grid members (None: uniform, the committed arithmetic)
        self.cond, self.log_prior = cond, log_prior
        self.params, self.Z = simulator_grid_24h(sc, cache_dir, grid, cond=cond)  # members x 24 x J
        self.n_hyd = n_hydraulic(grid) if local_hydraulic else 0
        if self.n_hyd > 1:
            check_grid_order(self.params, self.n_hyd)
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}
        self.J = len(sc.junctions)
        self._grp_mean = None                                            # lazily: mean over hydraulic siblings
        self.zmin_ = None             # S x J draws of ln(daily min) from the last predict_daily_min; a refit clears them
        self.age = sc.age_by_hour_h.loc[HOURS, sc.junctions].values     # 24 x J, nominal model
        self.static = X.loc[sc.junctions, STATIC].values                # J x |STATIC|
        # design rows for every (junction, hour), junction-major; scaling is fixed by this full grid so
        # the hour features are not re-scaled by whichever hours happened to be sampled
        jj = np.repeat(np.arange(self.J), 24); hh = np.tile(np.arange(24), self.J)
        self.Xall_ = self._raw(jj, hh)
        self.mu_, self.sd_ = self.Xall_.mean(0), self.Xall_.std(0)
        self.sd_[self.sd_ == 0] = 1.0

    def _raw(self, jidx: np.ndarray, hours: np.ndarray) -> np.ndarray:
        ang = 2 * np.pi * np.asarray(hours) / 24.0
        return np.column_stack([self.age[hours, jidx], self.static[jidx], np.sin(ang), np.cos(ang)])

    def _design(self, jidx, hours) -> np.ndarray:
        return (self._raw(np.asarray(jidx), np.asarray(hours)) - self.mu_) / self.sd_

    def fit_prior(self) -> "SimGP24":
        """No samples yet: uniform weights over grid members and doses, no discrepancy GP.  This is what
        the operator gets from the .inp alone, and what the first route is planned from."""
        W = np.full((len(self.params), len(self.doses)), 1.0 / (len(self.params) * len(self.doses)))
        if self.log_prior is not None:
            p = np.exp(np.asarray(self.log_prior, dtype=float) - logsumexp(self.log_prior))
            W = np.repeat(p[:, None] / len(self.doses), len(self.doses), axis=1)
        offs = np.log(np.asarray(self.doses, dtype=float))
        w = W.sum(axis=1)
        self.w_, self.W_, self.offs_ = w, W, offs
        self.m_, self.v_ = posterior_moments(W, offs, self.Z)
        self.hv_ = local_hydraulic_var(w, self.Z, self.n_hyd) if self.n_hyd > 1 else np.zeros_like(self.v_)
        self.map_params_, self.map_dose_ = None, None
        self.gp = None
        self.samples_ = pd.DataFrame(columns=["junction", "hour", "y"])
        self.zmin_ = None
        return self

    def fit(self, samples: pd.DataFrame) -> "SimGP24":
        """samples: columns junction, hour (int 0-23), y (mg/L)."""
        self.zmin_ = None
        if samples is None or len(samples) == 0:
            return self.fit_prior()
        idx = np.array([self.jidx[j] for j in samples.junction])
        h = samples.hour.astype(int).values
        z_obs = np.log(np.clip(samples.y.values, FLOOR, None))
        W, offs = grid_dose_weights(self.Z[:, h, idx], z_obs, self.lik_sd, self.lik, self.nu, self.doses, self.log_prior)
        w = W.sum(axis=1)
        self.w_, self.W_, self.offs_ = w, W, offs
        self.m_, self.v_ = posterior_moments(W, offs, self.Z)          # 24 x J
        self.hv_ = local_hydraulic_var(w, self.Z, self.n_hyd) if self.n_hyd > 1 else np.zeros_like(self.v_)
        k_, d_ = np.unravel_index(int(np.argmax(W)), W.shape)
        self.map_params_, self.map_dose_ = self.params[k_], float(np.exp(offs[d_]))
        r = z_obs - self.m_[h, idx]
        Xs = self._design(idx, h)
        # the discrepancy is driven by the diurnal demand pattern: it cannot change within an hour.  Inputs
        # that vary with the hour (age at that hour, sin, cos) get a length-scale floor of ~3 h, so the
        # joint 24-h draws behind the daily minimum are coherent rather than white noise (whose minimum
        # over 24 values is biased low by ~2 sd)
        lo = np.full(Xs.shape[1], 0.1)
        if self.smooth_hours:
            lo[[0, -2, -1]] = 1.0
        k = (ConstantKernel(0.1, (1e-3, 5.0))
             * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=list(zip(lo, np.full(Xs.shape[1], 20.0))), nu=1.5)
             + WhiteKernel(0.01, (1e-4, 0.5)))
        self.gp = GaussianProcessRegressor(kernel=k, n_restarts_optimizer=4, random_state=self.seed)
        self.gp.fit(Xs, r)
        self.samples_ = samples.copy()
        return self

    def predict_hours(self) -> tuple[np.ndarray, np.ndarray]:
        """Predictive ln C for every (hour, junction): (z_mu, z_sd), each 24 x J."""
        if self.gp is None:                                             # prior: simulator only
            self.z_sd_acq_ = np.sqrt(self.v_)
            return self.m_.copy(), np.sqrt(self.v_ + self.hv_)
        Xs = (self.Xall_ - self.mu_) / self.sd_
        r_mu, r_sd = self.gp.predict(Xs, return_std=True)
        noise = self.gp.kernel_.k2.noise_level
        r_var = np.clip(r_sd ** 2 - noise, 1e-4, None)
        z_mu = self.m_ + r_mu.reshape(self.J, 24).T
        z_sd = np.sqrt(self.v_ + self.hv_ + r_var.reshape(self.J, 24).T)
        self.z_sd_acq_ = np.sqrt(self.v_ + r_var.reshape(self.J, 24).T)   # reducible part, for acquisition
        return z_mu, z_sd

    def predict_hour(self, hour: int, hourly: tuple[np.ndarray, np.ndarray] | None = None) -> pd.DataFrame:
        """Same columns as SimGP.predict, for one hour of the day."""
        z_mu, z_sd = hourly if hourly is not None else self.predict_hours()
        out = pd.DataFrame(index=self.sc.junctions)
        out["median"] = np.exp(z_mu[hour]); out["lo90"] = np.exp(z_mu[hour] - 1.645 * z_sd[hour])
        out["hi90"] = np.exp(z_mu[hour] + 1.645 * z_sd[hour]); out["z_mu"], out["z_sd"] = z_mu[hour], z_sd[hour]
        return out

    def _gp_blocks(self):
        """Posterior mean and 24x24 covariance blocks of the discrepancy GP for every junction, without
        forming the full (24J)^2 matrix:  C_j = k(X_j, X_j) - V_j^T V_j,  V = L^{-1} k(X_train, X_*).
        Returns (r_mu [J, 24], C [J, 24, 24]) — the covariance of the LATENT field (noise removed)."""
        gp = self.gp
        Xs = (self.Xall_ - self.mu_) / self.sd_
        K_trans = gp.kernel_(Xs, gp.X_train_)                                # (24J) x n
        r_mu = (K_trans @ gp.alpha_).reshape(self.J, 24)
        V = solve_triangular(gp.L_, K_trans.T, lower=True).reshape(-1, self.J, 24)   # n x J x 24
        VtV = np.einsum("nja,njb->jab", V, V)
        k1, k2 = gp.kernel_.k1, gp.kernel_.k2                                # Constant * Matern, White
        Xb = Xs.reshape(self.J, 24, -1) / k1.k2.length_scale
        d = np.sqrt(np.maximum(((Xb[:, :, None, :] - Xb[:, None, :, :]) ** 2).sum(-1), 0.0))
        if k1.k2.nu == 1.5:
            Kb = k1.k1.constant_value * (1.0 + np.sqrt(3.0) * d) * np.exp(-np.sqrt(3.0) * d)
        else:                                                                # any other kernel: per-block fallback
            Kb = np.stack([k1(Xs[24 * j:24 * j + 24]) for j in range(self.J)])
        C = Kb - VtV
        return r_mu, 0.5 * (C + C.transpose(0, 2, 1))

    def predict_daily_min(self) -> pd.DataFrame:
        """Per junction: median / 90% band of the daily minimum, P(daily min < self.threshold), and the mean/sd
        of ln(daily min) for acquisition.  Monte Carlo over grid members and joint 24-h GP draws.  The draws stay
        on the model (self.zmin_), so P(daily min < t) at any other t is a Monte-Carlo count through
        self.p_below_mc(t) or SimGP24.p_below(frame, t, model=self).  The frame's attrs hold only the threshold
        its p_below column was computed at (a float)."""
        rng = np.random.default_rng(self.seed)
        S = self.n_draws
        flat = rng.choice(self.W_.size, size=S, p=self.W_.ravel())
        members, doses = np.unravel_index(flat, self.W_.shape)
        Zk = self.Z[members] + self.offs_[doses][:, None, None]
        if self.n_hyd > 1:
            # local hydraulic error: add the whole-day deviation of a random hydraulic sibling of the
            # drawn member from its decay group's mean (same variance as local_hydraulic_var, but as a
            # correlated 24-h profile so the minimum is taken over a coherent day)
            nd = len(self.w_) // self.n_hyd
            if self._grp_mean is None:
                self._grp_mean = self.Z.reshape(nd, self.n_hyd, 24, self.J).mean(axis=1)
            grp = members // self.n_hyd
            sib = grp * self.n_hyd + rng.integers(self.n_hyd, size=S)
            Zk = Zk + self.Z[sib] - self._grp_mean[grp]
        R = np.zeros((S, 24, self.J))
        if self.gp is not None:
            r_mu, C = self._gp_blocks()
            # about 5 of the 24 eigenvalues per block sit below the 1e-6 floor; their eigenvectors are
            # arbitrary, so a draw built from them changes with any 1e-16 change upstream.  Rebuild the
            # floored matrix (basis-independent) and take its Cholesky factor (unique) instead.
            w_, Vv = np.linalg.eigh(C)                                       # batched over junctions
            C_psd = (Vv * np.clip(w_, 1e-6, None)[:, None, :]) @ Vv.transpose(0, 2, 1)
            L = np.linalg.cholesky(0.5 * (C_psd + C_psd.transpose(0, 2, 1)))
            eps = rng.standard_normal((self.J, 24, S))
            R = (r_mu[:, :, None] + L @ eps).transpose(2, 1, 0)              # S x 24 x J
        zmin = (Zk + R).min(axis=1)                                     # S x J  ln(daily min)
        out = pd.DataFrame(index=self.sc.junctions)
        out["median"] = np.exp(np.median(zmin, axis=0))
        for q in (50, 80, 90, 95):
            a = (1 - q / 100) / 2
            out[f"lo{q}"] = np.exp(np.quantile(zmin, a, axis=0))
            out[f"hi{q}"] = np.exp(np.quantile(zmin, 1 - a, axis=0))
        out["z_mu"], out["z_sd"] = zmin.mean(axis=0), zmin.std(axis=0)
        out["p_below"] = (zmin < np.log(self.threshold)).mean(axis=0)
        out["argmin_hour"] = np.argmin(self.m_, axis=0)
        self.zmin_ = zmin
        out.attrs["p_below_threshold"] = self.threshold
        return out

    def p_below_mc(self, threshold: float) -> pd.Series:
        """Monte-Carlo P(daily min < threshold) per junction, from the draws of the last predict_daily_min."""
        if getattr(self, "zmin_", None) is None:
            raise RuntimeError("call predict_daily_min first: there are no daily-minimum draws yet")
        return pd.Series((self.zmin_ < np.log(threshold)).mean(axis=0), index=self.sc.junctions)

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2, model: "SimGP24 | None" = None) -> pd.Series:
        """P(value < threshold).  A frame from predict_daily_min returns its own p_below column at the threshold
        it was computed for (attrs['p_below_threshold'], 0.2 when absent).  At any other threshold:
          * with model=<the SimGP24 that produced the frame>: a Monte-Carlo count over that model's draws of the
            daily minimum, for the frame's junctions (the frame's z_mu must still be the one those draws give,
            so an edited or foreign frame raises instead of reading stale draws);
          * without a model (and for every other frame, e.g. an hourly prediction): the normal approximation on
            ln C, as in the committed code."""
        if "p_below" in pred and threshold == pred.attrs.get("p_below_threshold", 0.2):
            return pred["p_below"]
        if model is not None:
            zmin = getattr(model, "zmin_", None)
            if zmin is None:
                raise RuntimeError("model has no daily-minimum draws: call model.predict_daily_min() first")
            cols = pd.Index(model.sc.junctions).get_indexer(pred.index)
            if (cols < 0).any():
                raise KeyError("the frame's index has junctions the model's draws do not cover")
            if "z_mu" not in pred or not np.array_equal(pred["z_mu"].to_numpy(), zmin.mean(axis=0)[cols]):
                raise ValueError("the frame's z_mu is not the one the model's last draws give: it was edited, or it "
                                 "comes from another model or an earlier predict_daily_min call")
            return pd.Series((zmin[:, cols] < np.log(threshold)).mean(axis=0), index=pred.index)
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)
