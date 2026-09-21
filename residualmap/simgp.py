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


def simulator_grid(sc, cache_dir: str = "outputs/cache") -> tuple[list[tuple], np.ndarray]:
    """All grid members' ln C at the sampling hour: (params list, array [n_members, n_junctions])."""
    os.makedirs(cache_dir, exist_ok=True)
    f = os.path.join(cache_dir, f"grid_{os.path.basename(sc.wn_name)}_h{sc.sample_hour}.pkl")
    if os.path.exists(f):
        with open(f, "rb") as fh:
            return pickle.load(fh)
    params, rows = [], []
    for kb, kw, g in itertools.product(KB_GRID, KW_GRID, GAMMA_GRID):
        c = simulate_nominal_chlorine(sc.wn_name, kb, kw, g, sc.source_dose)
        params.append((kb, kw, g))
        rows.append(np.log(np.clip(c.loc[sc.sample_hour, sc.junctions].values, FLOOR, None)))
    out = (params, np.vstack(rows))
    with open(f, "wb") as fh:
        pickle.dump(out, fh)
    return out


class SimGP:
    def __init__(self, sc, seed: int = 0, columns: list[str] | None = None,
                 lik_sd: float = 0.25, cache_dir: str = "outputs/cache"):
        self.sc = sc
        self.seed = seed
        self.columns = columns or CORE
        self.lik_sd = lik_sd             # log-space tolerance when scoring grid members
        self.params, self.Z = simulator_grid(sc, cache_dir)   # Z: members x junctions
        self.jidx = {j: i for i, j in enumerate(sc.junctions)}

    def _calibrate(self, nodes: list[str], y: np.ndarray) -> None:
        idx = [self.jidx[n] for n in nodes]
        z_obs = np.log(np.clip(y, FLOOR, None))
        ll = -0.5 * (((self.Z[:, idx] - z_obs[None, :]) / self.lik_sd) ** 2).sum(axis=1)
        w = np.exp(ll - logsumexp(ll))
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
