"""
pinn.py — graph physics-informed neural network (PINN) baseline.

A PINN adds a physics residual to the data loss.  For a pipe network the natural physics is not a
continuous PDE (EPANET already discretises it) but a per-pipe balance: for a pipe that flows one way
all day from node i to node j, first-order decay says

        ln C_j - ln C_i  =  -k * (age_j - age_i)     (k unknown, learned)

so the physics loss is the mean squared violation of that relation over all steadily-directed pipes,
evaluated at ALL nodes (labels not needed) — that is what makes it "physics-informed".

Architecture: tiny MLP (features -> 16 -> 8 -> 1) in plain numpy, trained by L-BFGS.  A deep
ensemble (5 restarts) supplies a rough uncertainty.  No torch dependency so it runs anywhere.

Purpose: answer "would a PINN help?" with a number instead of a guess.  Expectation going in:
with 3-15 samples the network is data-starved and the GP with a simulator prior wins.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

from .features import rich_columns

FLOOR = 0.02


def _unpack(theta, sizes):
    ws, i = [], 0
    for a, b in zip(sizes[:-1], sizes[1:]):
        W = theta[i:i + a * b].reshape(a, b); i += a * b
        bvec = theta[i:i + b]; i += b
        ws.append((W, bvec))
    k = theta[i]
    return ws, k


def _forward(ws, X):
    h = X
    for W, b in ws[:-1]:
        h = np.tanh(h @ W + b)
    W, b = ws[-1]
    return (h @ W + b).ravel()


class GraphPINN:
    def __init__(self, sc, columns: list[str] | None = None, hidden=(16, 8), lam_phys: float = 1.0,
                 lam_l2: float = 1e-3, n_ensemble: int = 5, seed: int = 0):
        self.sc, self.hidden, self.lam_phys, self.lam_l2 = sc, hidden, lam_phys, lam_l2
        self.n_ensemble, self.seed = n_ensemble, seed
        self.columns = columns
        p = sc.pipes[sc.pipes.direction != 0]
        # unidirectional pipes as (upstream, downstream) junction pairs
        pairs = [(r.start, r.end) if r.direction > 0 else (r.end, r.start) for r in p.itertuples()]
        self.pairs = [(a, b) for a, b in pairs if a in sc.junctions and b in sc.junctions]

    def fit(self, X_all: pd.DataFrame, sampled: list[str], y: np.ndarray) -> "GraphPINN":
        cols = self.columns or rich_columns(X_all)
        Xc = X_all[cols]
        self.mu_, self.sd_ = Xc.mean(), Xc.std().replace(0, 1.0)
        Xs = ((Xc - self.mu_) / self.sd_).values
        jidx = {j: i for i, j in enumerate(X_all.index)}
        s_idx = np.array([jidx[j] for j in sampled])
        z_obs = np.log(np.clip(y, FLOOR, None))
        age = X_all["age_h"].values
        up = np.array([jidx[a] for a, _ in self.pairs]); dn = np.array([jidx[b] for _, b in self.pairs])
        dage = np.clip(age[dn] - age[up], 0.0, None)  # age can only grow downstream (else 0 -> no constraint)
        sizes = [Xs.shape[1], *self.hidden, 1]
        n_par = sum(a * b + b for a, b in zip(sizes[:-1], sizes[1:])) + 1
        z_mean = z_obs.mean()

        def loss(theta):
            ws, k = _unpack(theta, sizes)
            z = _forward(ws, Xs) + z_mean
            data = np.mean((z[s_idx] - z_obs) ** 2)
            phys = np.mean((z[dn] - z[up] + np.exp(k) * dage) ** 2) if len(up) else 0.0
            l2 = self.lam_l2 * np.sum(theta[:-1] ** 2)
            return data + self.lam_phys * phys + l2

        rng = np.random.default_rng(self.seed)
        self.members_ = []
        for _ in range(self.n_ensemble):
            th0 = np.concatenate([rng.normal(0, 0.3, n_par - 1), [np.log(0.02)]])
            res = minimize(loss, th0, method="L-BFGS-B", options={"maxiter": 400})
            self.members_.append(res.x)
        self.sizes_, self.z_mean_, self.Xs_ = sizes, z_mean, Xs
        self.k_ = float(np.mean([np.exp(_unpack(m, sizes)[1]) for m in self.members_]))
        return self

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        cols = list(self.mu_.index)
        Xs = ((X[cols] - self.mu_) / self.sd_).values
        preds = np.vstack([_forward(_unpack(m, self.sizes_)[0], Xs) + self.z_mean_ for m in self.members_])
        z_mu, z_sd = preds.mean(axis=0), np.clip(preds.std(axis=0), 0.05, None)
        out = pd.DataFrame(index=X.index)
        out["median"] = np.exp(z_mu)
        out["lo90"] = np.exp(z_mu - 1.645 * z_sd)
        out["hi90"] = np.exp(z_mu + 1.645 * z_sd)
        out["z_mu"], out["z_sd"] = z_mu, z_sd
        return out

    @staticmethod
    def p_below(pred: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
        return pd.Series(norm.cdf((np.log(threshold) - pred["z_mu"]) / pred["z_sd"]), index=pred.index)
