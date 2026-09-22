"""
route.py — a monthly sampling route, chosen all at once.

An operator does not pick one grab sample at a time; they drive a route.  plan_route picks K
(junction, hour) pairs from the current model — before any sample of the month has been taken, that
model is the operator's EPANET file alone (SimGP24.fit_prior) — greedily on the daily-minimum straddle
score (level-set estimation: where the model is least sure whether the daily minimum is above or below
the limit), with a repulsion in hydraulic-distance space so the K sites spread over the network, and at
most ceil(K / len(hours)) sites per hour so one person can actually drive it.

Baselines for the comparison: K random junctions, and the K highest-demand junctions (a plausible
operator heuristic: sample where the customers are).  Both on an evenly spread daytime schedule.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .simgp import DAY_HOURS


def spread_hours(K: int, hours: list[int] = DAY_HOURS) -> list[int]:
    """Route order 08:00 -> 16:00: K visits spread evenly over the working day."""
    lo, hi = hours[1] if len(hours) > 2 else hours[0], hours[-2] if len(hours) > 2 else hours[-1]
    return [int(round(h)) for h in np.linspace(lo, hi, K)]


def plan_route(model, K: int, hours: list[int] = DAY_HOURS, threshold: float = 0.2,
               repulsion: float = 1.0, length_scale_m: float | None = None,
               exclude: list[str] | None = None) -> pd.DataFrame:
    """K (junction, hour) pairs from a fitted (or prior) SimGP24.  Returns junction, hour, p_below,
    score and a one-sentence reason per site."""
    sc = model.sc
    z_mu, z_sd = model.predict_hours()
    sd_acq = pd.DataFrame(model.z_sd_acq_, columns=sc.junctions).loc[hours]
    dmin = model.predict_daily_min()
    D = sc.hyd_dist
    if length_scale_m is None:
        finite = D.values[np.isfinite(D.values) & (D.values > 0)]
        # 5% of the network's hydraulic diameter; with repulsion = 1.0 the route's mean nearest-neighbour
        # spacing matches a random route's on Net3 (2.2 / 1.8 / 1.3 km at K = 5 / 8 / 12)
        length_scale_m = 0.05 * float(np.max(finite))
    base = 1.96 * dmin["z_sd"] - (dmin["z_mu"] - np.log(threshold)).abs()
    candidates = [j for j in sc.junctions if not exclude or j not in set(exclude)]
    cap = math.ceil(K / len(hours))
    chosen, used = [], {h: 0 for h in hours}
    rows = []
    for _ in range(K):
        score = base.loc[candidates].copy()
        if chosen:
            dnear = D.loc[candidates, chosen].min(axis=1)
            score = score - repulsion * np.exp(-dnear / length_scale_m)
        j = str(score.idxmax())
        free = [h for h in hours if used[h] < cap]
        h = int(sd_acq.loc[free, j].idxmax())                       # the daytime hour that constrains it most
        used[h] += 1
        p = float(dmin.loc[j, "p_below"])
        near = f"{D.loc[j, chosen].min() / 1000:.1f} km from the nearest other route site" if chosen else "first site"
        why = ("most uncertain violation call" if abs(p - 0.5) < 0.25 else
               "likely violation, confirm it" if p >= 0.75 else "likely fine, widest band")
        rows.append({"junction": j, "hour": h, "p_below": round(p, 2), "score": round(float(base[j]), 3),
                     "reason": f"P(daily min < {threshold}) = {p:.2f} — {why}; sample at {h:02d}:00 where the band is widest by day; {near}."})
        chosen.append(j); candidates.remove(j)
    return pd.DataFrame(rows)


def random_route(sc, K: int, rng, hours: list[int] = DAY_HOURS) -> pd.DataFrame:
    js = [str(j) for j in rng.choice(sc.junctions, K, replace=False)]
    return pd.DataFrame({"junction": js, "hour": spread_hours(K, hours)})


def demand_route(sc, K: int, hours: list[int] = DAY_HOURS) -> pd.DataFrame:
    """The K highest base-demand junctions: sample where the customers are."""
    dem = pd.Series({j: sum(ts.base_value for ts in sc.wn.get_node(j).demand_timeseries_list) for j in sc.junctions})
    js = list(dem.sort_values(ascending=False).index[:K])
    return pd.DataFrame({"junction": js, "hour": spread_hours(K, hours)})
