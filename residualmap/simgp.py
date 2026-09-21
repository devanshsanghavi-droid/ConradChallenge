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


HOURS = list(range(24))
DAY_HOURS = list(range(7, 18))                  # 07:00-17:00: when an operator can take a grab sample
STATIC = ["emb0", "emb1", "emb2", "path_wall_index", "dist_src_km"]   # CORE minus the hour-dependent age


def simulator_grid_24h(sc, cache_dir: str = "outputs/cache") -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C for every hour of the last day: (params, array [members, 24, junctions]).

    Each EPANET run already produces the whole day; storing all of it is what lets the time-aware
    model (SimGP24) calibrate on daytime samples and predict the night."""
    os.makedirs(cache_dir, exist_ok=True)
    f = os.path.join(cache_dir, f"grid24_{os.path.basename(sc.wn_name)}.pkl")
    if os.path.exists(f):
        with open(f, "rb") as fh:
            return pickle.load(fh)
    params, rows = [], []
    for kb, kw, g in itertools.product(KB_GRID, KW_GRID, GAMMA_GRID):
        c = simulate_nominal_chlorine(sc.wn_name, kb, kw, g, sc.source_dose)
        params.append((kb, kw, g))
        rows.append(np.log(np.clip(c.loc[HOURS, sc.junctions].values, FLOOR, None)))
    out = (params, np.stack(rows))
    with open(f, "wb") as fh:
        pickle.dump(out, fh)
    return out


def simulator_grid(sc, cache_dir: str = "outputs/cache") -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C at the sampling hour: (params list, array [n_members, n_junctions])."""
    params, Z = simulator_grid_24h(sc, cache_dir)
    return params, Z[:, sc.sample_hour, :]


def grid_weights(z_sim: np.ndarray, z_obs: np.ndarray, lik_sd: float, lik: str = "gauss",
                 nu: float = 3.0) -> np.ndarray:
    """Posterior weight of each grid member given log-samples.  z_sim: members x samples.

    lik="gauss" is iteration 2.  lik="t" (Student-t, nu degrees of freedom) is robust: a junction where
    the operator's model is structurally wrong (a front that sits elsewhere, a tank zone that turns over
    at a different hour) can be 3x off for EVERY grid member; under a Gaussian one such reading drags
    the whole calibration, under a Student-t it is discounted."""
    e = (z_sim - z_obs[None, :]) / lik_sd
    if lik == "gauss":
        ll = -0.5 * (e ** 2).sum(axis=1)
    elif lik == "t":
        ll = (-(nu + 1) / 2 * np.log1p(e ** 2 / nu)).sum(axis=1)
    else:
        raise ValueError(lik)
    return np.exp(ll - logsumexp(ll))


class SimGP:
    def __init__(self, sc, seed: int = 0, columns: list[str] | None = None,
                 lik_sd: float = 0.25, cache_dir: str = "outputs/cache", lik: str = "gauss", nu: float = 3.0):
        self.sc = sc
        self.seed = seed
        self.columns = columns or CORE
        self.lik_sd, self.lik, self.nu = lik_sd, lik, nu   # log-space tolerance / family when scoring grid members
        self.params, self.Z = simulator_grid(sc, cache_dir)   # Z: members x junctions
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}

    def _calibrate(self, nodes: list[str], y: np.ndarray) -> None:
        idx = [self.jidx[n] for n in nodes]
        z_obs = np.log(np.clip(y, FLOOR, None))
        w = grid_weights(self.Z[:, idx], z_obs, self.lik_sd, self.lik, self.nu)
        self.w_ = w
        self.m_ = w @ self.Z                              # per junction
        self.v_ = w @ (self.Z - self.m_[None, :]) ** 2    # per junction
        self.map_params_ = self.params[int(np.argmax(w))]

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
        z_sd = np.sqrt(self.v_[idx] + r_var)
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

    def __init__(self, sc, X: pd.DataFrame, seed: int = 0, lik_sd: float = 0.25,
                 cache_dir: str = "outputs/cache", n_draws: int = 256, lik: str = "t", nu: float = 3.0):
        self.sc, self.seed, self.lik_sd, self.n_draws = sc, seed, lik_sd, n_draws
        self.lik, self.nu = lik, nu
        self.params, self.Z = simulator_grid_24h(sc, cache_dir)        # members x 24 x J
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
        w = grid_weights(self.Z[:, h, idx], z_obs, self.lik_sd, self.lik, self.nu)
        self.w_ = w
        self.m_ = np.tensordot(w, self.Z, axes=1)                       # 24 x J
        self.v_ = np.tensordot(w, (self.Z - self.m_[None]) ** 2, axes=1)
        self.map_params_ = self.params[int(np.argmax(w))]
        r = z_obs - self.m_[h, idx]
        Xs = self._design(idx, h)
        k = (ConstantKernel(0.1, (1e-3, 5.0))
             * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=(0.1, 20.0), nu=1.5)
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
        z_sd = np.sqrt(self.v_ + r_var.reshape(self.J, 24).T)
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
        members = rng.choice(len(self.w_), size=S, p=self.w_)
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
        zmin = (self.Z[members] + R).min(axis=1)                        # S x J  ln(daily min)
        out = pd.DataFrame(index=self.sc.junctions)
        out["median"] = np.exp(np.median(zmin, axis=0))
        out["lo90"] = np.exp(np.quantile(zmin, 0.05, axis=0))
        out["hi90"] = np.exp(np.quantile(zmin, 0.95, axis=0))
        out["z_mu"], out["z_sd"] = zmin.mean(axis=0), zmin.std(axis=0)
        out["p_below"] = (zmin < np.log(0.2)).mean(axis=0)
        out["argmin_hour"] = np.argmin(self.m_, axis=0)
        return out

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
        if "p_below" in pred and threshold == 0.2:
            return pred["p_below"]
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)
