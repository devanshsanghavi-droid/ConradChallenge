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

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from .features import CORE
from .simulate import simulate_nominal_chlorine

FLOOR = 0.02
KB_GRID = [0.10, 0.25, 0.40, 0.55, 0.70]        # 1/day
KW_GRID = [0.10, 0.30, 0.60, 1.00, 1.50]        # m/day
GAMMA_GRID = [0.0, 0.5, 1.0]                    # old-pipe sensitivity
DEMAND_GRID = [0.85, 1.0, 1.15]                 # global demand multiplier (hydraulic mismatch)
ROUGH_GRID = [0.9, 1.0, 1.1]                    # global Hazen-Williams C multiplier (hydraulic mismatch)
GRIDS = {"decay": (KB_GRID, KW_GRID, GAMMA_GRID, [1.0], [1.0]),                    # iteration 2: 75 runs
         "full": (KB_GRID, KW_GRID, GAMMA_GRID, DEMAND_GRID, ROUGH_GRID)}         # iteration 3: 675 runs
DOSE_GRID = [0.90, 0.95, 1.00, 1.05, 1.10]     # source-dose multiplier: first-order decay is linear in
                                                # concentration, so this axis is an exact ln-offset, no runs


HOURS = list(range(24))
DAY_HOURS = list(range(7, 18))                  # 07:00-17:00: when an operator can take a grab sample
STATIC = ["emb0", "emb1", "emb2", "path_wall_index", "dist_src_km"]   # CORE minus the hour-dependent age


def simulator_grid_24h(sc, cache_dir: str = "outputs/cache", grid: str = "full") -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C for every hour of the last day: (params, array [members, 24, junctions]).
    params are (kb, kw, gamma, demand_mult, rough_mult).

    Each EPANET run already produces the whole day; storing all of it is what lets the time-aware
    model (SimGP24) calibrate on daytime samples and predict the night.  grid="full" adds the two
    hydraulic-mismatch axes (675 runs); grid="decay" is the iteration-2 grid (75 runs)."""
    os.makedirs(cache_dir, exist_ok=True)
    f = os.path.join(cache_dir, f"grid24_{grid}_{os.path.basename(sc.wn_name)}.pkl")
    if os.path.exists(f):
        with open(f, "rb") as fh:
            return pickle.load(fh)
    params, rows = [], []
    for kb, kw, g, dm, rm in itertools.product(*GRIDS[grid]):
        c = simulate_nominal_chlorine(sc.wn_name, kb, kw, g, sc.source_dose, dm, rm)
        params.append((kb, kw, g, dm, rm))
        rows.append(np.log(np.clip(c.loc[HOURS, sc.junctions].values, FLOOR, None)).astype(np.float32))
    out = (params, np.stack(rows))
    with open(f, "wb") as fh:
        pickle.dump(out, fh)
    return out


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
                      doses=DOSE_GRID) -> tuple[np.ndarray, np.ndarray]:
    """Joint posterior over (grid member, dose multiplier): W [members x doses], and the ln-offsets."""
    offs = np.log(np.asarray(doses, dtype=float))
    ll = np.stack([grid_loglik(z_sim, z_obs - d, lik_sd, lik, nu) for d in offs], axis=1)
    return np.exp(ll - logsumexp(ll)), offs


LIK_SD = 0.35   # log-space scale of the calibration likelihood = the day-time RMS mismatch between the truth
                # and the best grid member on Net3 (0.33), i.e. model error, not grab-sample noise (0.03 mg/L)


def n_hydraulic(grid: str) -> int:
    return len(GRIDS[grid][3]) * len(GRIDS[grid][4])


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
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}

    def _calibrate(self, nodes: list[str], y: np.ndarray) -> None:
        idx = [self.jidx[n] for n in nodes]
        z_obs = np.log(np.clip(y, FLOOR, None))
        W, offs = grid_dose_weights(self.Z[:, idx], z_obs, self.lik_sd, self.lik, self.nu, self.doses)
        w, wd = W.sum(axis=1), W.sum(axis=0)              # marginals over members / over dose
        self.w_, self.W_, self.offs_ = w, W, offs
        self.m_ = w @ self.Z + wd @ offs                  # per junction
        self.v_ = w @ self.Z ** 2 + 2 * (w * (W @ offs) / np.maximum(w, 1e-300)) @ self.Z + wd @ offs ** 2 - self.m_ ** 2
        self.v_ = np.clip(self.v_, 1e-6, None)
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

    def __init__(self, sc, X: pd.DataFrame, seed: int = 0, lik_sd: float = LIK_SD,
                 cache_dir: str = "outputs/cache", n_draws: int = 256, lik: str = "t", nu: float = 3.0,
                 grid: str = "full", local_hydraulic: bool = True, doses=DOSE_GRID, smooth_hours: bool = True):
        self.sc, self.seed, self.lik_sd, self.n_draws = sc, seed, lik_sd, n_draws
        self.lik, self.nu, self.doses, self.smooth_hours = lik, nu, doses, smooth_hours
        self.params, self.Z = simulator_grid_24h(sc, cache_dir, grid)  # members x 24 x J
        self.n_hyd = n_hydraulic(grid) if local_hydraulic else 0
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}
        self.J = len(sc.junctions)
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

    def fit(self, samples: pd.DataFrame) -> "SimGP24":
        """samples: columns junction, hour (int 0-23), y (mg/L)."""
        idx = np.array([self.jidx[j] for j in samples.junction])
        h = samples.hour.astype(int).values
        z_obs = np.log(np.clip(samples.y.values, FLOOR, None))
        W, offs = grid_dose_weights(self.Z[:, h, idx], z_obs, self.lik_sd, self.lik, self.nu, self.doses)
        w, wd = W.sum(axis=1), W.sum(axis=0)
        self.w_, self.W_, self.offs_ = w, W, offs
        self.m_ = np.tensordot(w, self.Z, axes=1) + wd @ offs           # 24 x J
        # E[(Z+d)^2] - m^2 with (member, dose) jointly distributed
        wo = W @ offs                                                   # per member: sum_d W[k,d] offs[d]
        self.v_ = (np.tensordot(w, self.Z ** 2, axes=1) + 2 * np.tensordot(wo, self.Z, axes=1)
                   + wd @ offs ** 2 - self.m_ ** 2)
        self.v_ = np.clip(self.v_, 1e-6, None)
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

    def predict_daily_min(self) -> pd.DataFrame:
        """Per junction: median / 90% band of the daily minimum, P(daily min < 0.2), and the mean/sd of
        ln(daily min) for acquisition.  Monte Carlo over grid members and joint 24-h GP draws."""
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
            Zr = self.Z.reshape(nd, self.n_hyd, 24, self.J)
            grp = members // self.n_hyd
            sib = grp * self.n_hyd + rng.integers(self.n_hyd, size=S)
            Zk = Zk + self.Z[sib] - Zr.mean(axis=1)[grp]
        Xs = (self.Xall_ - self.mu_) / self.sd_
        r_mu, r_cov = self.gp.predict(Xs, return_cov=True)
        noise = self.gp.kernel_.k2.noise_level
        R = np.empty((S, 24, self.J))
        for j in range(self.J):
            sl = slice(24 * j, 24 * j + 24)
            C = r_cov[sl, sl] - noise * np.eye(24)
            C = 0.5 * (C + C.T)
            w_, V = np.linalg.eigh(C)
            L = V * np.sqrt(np.clip(w_, 1e-6, None))
            R[:, :, j] = r_mu[sl] + (L @ rng.standard_normal((24, S))).T
        zmin = (Zk + R).min(axis=1)                                     # S x J  ln(daily min)
        out = pd.DataFrame(index=self.sc.junctions)
        out["median"] = np.exp(np.median(zmin, axis=0))
        for q in (50, 80, 90, 95):
            a = (1 - q / 100) / 2
            out[f"lo{q}"] = np.exp(np.quantile(zmin, a, axis=0))
            out[f"hi{q}"] = np.exp(np.quantile(zmin, 1 - a, axis=0))
        out["z_mu"], out["z_sd"] = zmin.mean(axis=0), zmin.std(axis=0)
        out["p_below"] = (zmin < np.log(0.2)).mean(axis=0)
        out["argmin_hour"] = np.argmin(self.m_, axis=0)
        return out

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
        if "p_below" in pred and threshold == 0.2:
            return pred["p_below"]
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)
