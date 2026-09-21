"""
surrogate.py — physics-informed GP for chlorine residual, baselines, and acquisition rules.

Model (in log space so predictions stay positive and errors are relative):

    z = ln C
    z = a + b * age + f(x) + eps,     f ~ GP(0, Matern_ARD),  x = [age, MDS(hydraulic distance)]

* a, b : ridge-regularised log-linear decay fit (first-order decay says ln C is linear in age)
* f    : learns what the decay law misses — pipe-wall heterogeneity, mixing zones, tank effects
* eps  : grab-sample noise (learned WhiteKernel)

Outputs per junction: median C, 90% band, and P(C < threshold) — the compliance question.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

FLOOR = 0.02  # mg/L, below this a chlorine reading is effectively "none"


from .features import CORE, build_features  # noqa: F401  (re-exported for callers)


# ----------------------------------------------------------------------------- model
class PhysicsGP:
    """Log-linear decay mean function + GP residual."""

    def __init__(self, prior_b: float = -0.02, ridge: float = 0.5, seed: int = 0,
                 columns: list[str] | None = None):
        self.prior_b = prior_b     # prior slope: ln C loses 0.02 per hour of age (~38 %/day)
        self.ridge = ridge         # shrinkage of the fitted slope toward prior_b
        self.seed = seed
        self.columns = columns or CORE
        self.gp = None

    # mean function -----------------------------------------------------------
    def _fit_mean(self, age: np.ndarray, z: np.ndarray) -> None:
        # ridge regression of (z - prior_b*age) on age, shrinking slope toward prior
        A = np.column_stack([np.ones_like(age), age])
        target = z - self.prior_b * age
        reg = np.diag([0.0, self.ridge])
        coef = np.linalg.solve(A.T @ A + reg, A.T @ target)
        self.a_, self.b_ = coef[0], coef[1] + self.prior_b

    def _mean(self, age: np.ndarray) -> np.ndarray:
        return self.a_ + self.b_ * age

    # fit / predict -----------------------------------------------------------
    def fit(self, X: pd.DataFrame, y: np.ndarray) -> "PhysicsGP":
        X = X[self.columns]
        z = np.log(np.clip(y, FLOOR, None))
        age = X["age_h"].values
        self._fit_mean(age, z)
        r = z - self._mean(age)
        self.mu_, self.sd_ = X.mean(), X.std().replace(0, 1.0)
        Xs = ((X - self.mu_) / self.sd_).values
        k = (ConstantKernel(0.3, (1e-3, 10.0))
             * Matern(length_scale=np.ones(Xs.shape[1]), length_scale_bounds=(0.1, 20.0), nu=1.5)
             + WhiteKernel(0.01, (1e-4, 0.5)))
        self.gp = GaussianProcessRegressor(kernel=k, normalize_y=False, n_restarts_optimizer=4,
                                           random_state=self.seed)
        self.gp.fit(Xs, r)
        return self

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        X = X[self.columns]
        Xs = ((X - self.mu_) / self.sd_).values
        r_mu, r_sd = self.gp.predict(Xs, return_std=True)
        z_mu = self._mean(X["age_h"].values) + r_mu
        # sklearn's predictive variance includes the WhiteKernel noise; remove it so the
        # band describes the latent residual field, not a future noisy grab sample
        noise = self.gp.kernel_.k2.noise_level
        z_sd = np.sqrt(np.clip(r_sd**2 - noise, 1e-4, None))
        out = pd.DataFrame(index=X.index)
        out["median"] = np.exp(z_mu)
        out["lo90"] = np.exp(z_mu - 1.645 * z_sd)
        out["hi90"] = np.exp(z_mu + 1.645 * z_sd)
        out["z_mu"], out["z_sd"] = z_mu, z_sd
        return out

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)


# ----------------------------------------------------------------------------- baselines
def baseline_mean(sc, sampled: list[str], y: np.ndarray) -> pd.Series:
    """What an operator implicitly assumes: the whole system looks like my samples."""
    return pd.Series(float(np.mean(y)), index=sc.junctions)


def baseline_nearest(sc, sampled: list[str], y: np.ndarray) -> pd.Series:
    """Copy the reading from the hydraulically nearest sampled node."""
    d = sc.hyd_dist.loc[:, sampled].values
    return pd.Series(np.asarray(y)[d.argmin(axis=1)], index=sc.junctions)


def baseline_decay_only(sc, X: pd.DataFrame, sampled: list[str], y: np.ndarray) -> pd.Series:
    """First-order decay law fitted to samples, no GP (the mean function alone)."""
    m = PhysicsGP()
    m._fit_mean(X.loc[sampled, "age_h"].values, np.log(np.clip(y, FLOOR, None)))
    return pd.Series(np.exp(m._mean(X["age_h"].values)), index=sc.junctions)


# ----------------------------------------------------------------------------- acquisition
def acquire(strategy: str, pred: pd.DataFrame, candidates: list[str], rng,
            threshold: float = 0.2) -> str:
    """Pick the next junction to sample."""
    p = pred.loc[candidates]
    if strategy == "random":
        return str(rng.choice(candidates))
    if strategy == "uncertainty":
        return str(p["z_sd"].idxmax())
    if strategy == "straddle":
        # level-set estimation (Bryan et al. 2005): most ambiguous w.r.t. the threshold
        score = 1.96 * p["z_sd"] - (p["z_mu"] - np.log(threshold)).abs()
        return str(score.idxmax())
    raise ValueError(strategy)
