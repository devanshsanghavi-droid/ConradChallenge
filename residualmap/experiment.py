"""
experiment.py — three questions, answered on junctions that were NOT sampled:

  Q1  Does a physics-informed surrogate beat what an operator does today?
  Q2  Which physics helps at 3-15 samples: the decay law, the full simulator, more features, a PINN?
  Q3  Does smart sampling (uncertainty / straddle) beat random sampling?

Plus the finding that drives iteration 3: the network is worst at hours nobody samples.
"""
from __future__ import annotations

import json
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wntr

from .features import CORE, build_features, rich_columns
from .pinn import GraphPINN
from .route import demand_route, plan_route, random_route
from .simgp import DAY_HOURS, SimGP, SimGP24
from .simulate import build_scenario
from .surrogate import (PhysicsGP, acquire, acquire_time, baseline_decay_only, baseline_mean,
                        baseline_nearest)

warnings.filterwarnings("ignore")  # GP optimizer bound warnings are expected with 3-15 points
STRATEGIES = ["random", "uncertainty", "straddle"]
THRESHOLD = 0.2
PINN_AT = (3, 5, 8, 10, 12, 15)
TIME_STRATEGIES = ["random", "uncertainty", "straddle", "straddle_min"]
ROUTE_KS = (5, 8, 12)
TIME_MAIN = "straddle_min"      # rule whose maps are drawn
NIGHT_HOUR = 22
BIG = {"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "legend.fontsize": 11,
       "xtick.labelsize": 11, "ytick.labelsize": 11}   # figures must read in a judged video


def _metrics(truth, pred_median, p_viol, lo, hi, mask) -> dict:
    t, m = truth.loc[mask], pred_median.loc[mask]
    out = {"rmse": float(np.sqrt(np.mean((t - m) ** 2))), "mae": float(np.mean(np.abs(t - m)))}
    tv = t < THRESHOLD
    pv = (p_viol.loc[mask] > 0.5) if p_viol is not None else (m < THRESHOLD)
    tp = int((tv & pv).sum()); fp = int((~tv & pv).sum()); fn = int((tv & ~pv).sum())
    prec = tp / (tp + fp) if tp + fp else 1.0
    rec = tp / (tp + fn) if tp + fn else 1.0
    out.update(precision=prec, recall=rec, f1=(2 * prec * rec / (prec + rec)) if prec + rec else 0.0,
               n_true_viol=int(tv.sum()))
    if lo is not None:
        out["coverage90"] = float(((t >= lo.loc[mask]) & (t <= hi.loc[mask])).mean())
    return out


LEVELS = (50, 80, 90, 95)   # nominal band levels for the reliability diagram


def _coverage_levels(truth, pred, mask, suffix="") -> dict:
    """Empirical coverage of the central 50/80/95% bands (90 is already in the metrics)."""
    from scipy.stats import norm
    t = truth.loc[mask]; z_mu, z_sd = pred.loc[mask, "z_mu"], pred.loc[mask, "z_sd"]
    out = {}
    for q in LEVELS:
        if q == 90:
            continue
        if f"lo{q}" in pred:
            lo, hi = pred.loc[mask, f"lo{q}"], pred.loc[mask, f"hi{q}"]
        else:
            k = norm.ppf(0.5 + q / 200); lo, hi = np.exp(z_mu - k * z_sd), np.exp(z_mu + k * z_sd)
        out[f"coverage{q}{suffix}"] = float(((t >= lo) & (t <= hi)).mean())
    return out


def _row(model_name, strategy, n, seed, truth, pred, unsampled, model_cls):
    pv = model_cls.p_below(pred, THRESHOLD)
    return {"model": model_name, "strategy": strategy, "n": n, "seed": seed,
            **_metrics(truth, pred["median"], pv, pred["lo90"], pred["hi90"], unsampled),
            **_coverage_levels(truth, pred, unsampled)}, pv


def run_scenario(sc, n_seed=3, n_max=15, noise_sd=0.03, seed=0, cache_dir="outputs/cache"):
    rng = np.random.default_rng(1000 + seed)
    X = build_features(sc)
    rich = rich_columns(X)
    truth = sc.truth_snapshot

    def observe(nodes):
        return np.clip(truth.loc[nodes].values + rng.normal(0, noise_sd, len(nodes)), 0.01, None)

    seed_nodes = list(rng.choice(sc.junctions, n_seed, replace=False))
    rows, snapshots = [], {}

    # ---- main model (SimGP) under each acquisition rule; other models ride along on 'random'
    for strat in STRATEGIES:
        sampled, y = list(seed_nodes), list(observe(seed_nodes))
        for n in range(n_seed, n_max + 1):
            unsampled = [j for j in sc.junctions if j not in sampled]
            Xs, ys = X.loc[sampled], np.array(y)
            main = SimGP(sc, seed=seed, cache_dir=cache_dir).fit(Xs, ys)
            pred = main.predict(X)
            r, pv = _row("simgp_core", strat, n, seed, truth, pred, unsampled, SimGP)
            r["map_kb"], r["map_kw"], r["map_gamma"], r["map_demand"], r["map_rough"] = main.map_params_
            r["map_dose"] = main.map_dose_
            # operational recall: a violation is found if its grab sample read below the limit OR the
            # model flags it — the unsampled-only recall above penalises a rule for sampling the violators
            found = int(((truth.loc[sampled] < THRESHOLD) & (np.array(y) < THRESHOLD)).sum()) + \
                    int(((truth.loc[unsampled] < THRESHOLD) & (pv.loc[unsampled] > 0.5)).sum())
            r["recall_all"] = found / max(int((truth < THRESHOLD).sum()), 1)
            rows.append(r)
            if strat == "random":
                # calibrated simulator alone — does the discrepancy GP absorb what the grid cannot?
                rows.append(_row("simgp_core_nogp", strat, n, seed, truth, main.predict_prior(X), unsampled, SimGP)[0])
                # iteration-2 grid (decay parameters only) for the before/after calibration comparison
                p0 = SimGP(sc, seed=seed, cache_dir=cache_dir, grid="decay", lik_sd=0.25, doses=[1.0],
                           lik="gauss").fit(Xs, ys).predict(X)   # exactly the iteration-2 model
                rows.append(_row("simgp_core_g75", strat, n, seed, truth, p0, unsampled, SimGP)[0])
                p1 = PhysicsGP(seed=seed, columns=CORE).fit(Xs, ys).predict(X)
                rows.append(_row("physgp_core", strat, n, seed, truth, p1, unsampled, PhysicsGP)[0])
                p2 = PhysicsGP(seed=seed, columns=rich).fit(Xs, ys).predict(X)
                rows.append(_row("physgp_rich", strat, n, seed, truth, p2, unsampled, PhysicsGP)[0])
                p3 = SimGP(sc, seed=seed, columns=rich, cache_dir=cache_dir).fit(Xs, ys).predict(X)
                rows.append(_row("simgp_rich", strat, n, seed, truth, p3, unsampled, SimGP)[0])
                if n in PINN_AT:
                    p4 = GraphPINN(sc, n_ensemble=3, seed=seed).fit(X, sampled, ys).predict(X)
                    rows.append(_row("pinn_rich", strat, n, seed, truth, p4, unsampled, GraphPINN)[0])
                for name, bp in [("mean_of_samples", baseline_mean(sc, sampled, y)),
                                 ("nearest_sample", baseline_nearest(sc, sampled, y)),
                                 ("decay_only", baseline_decay_only(sc, X, sampled, y))]:
                    rows.append({"model": name, "strategy": strat, "n": n, "seed": seed,
                                 **_metrics(truth, bp, None, None, None, unsampled)})
            if strat == "straddle" and n in (n_seed, 8, n_max):
                snapshots[n] = {"sampled": list(sampled), "pred": pred.copy(), "p_viol": pv.copy()}
            if n == n_max:
                break
            nxt = acquire(strat, pred, unsampled, rng, THRESHOLD)
            sampled.append(nxt); y.append(float(observe([nxt])[0]))
    return pd.DataFrame(rows), snapshots, X


# ----------------------------------------------------------------------------- time-aware (task 1)
def _metrics_time(sc, pmin: pd.DataFrame, p_night: pd.DataFrame | None, mask) -> dict:
    """Scores against the DAILY MINIMUM (the compliance number) and the night snapshot, on unsampled
    junctions.  pmin needs columns median / lo90 / hi90 / p_below."""
    t, m = sc.truth_daily_min.loc[mask], pmin.loc[mask]
    tv, pv = t < THRESHOLD, m["p_below"] > 0.5
    tp = int((tv & pv).sum()); fp = int((~tv & pv).sum()); fn = int((tv & ~pv).sum())
    prec = tp / (tp + fp) if tp + fp else 1.0
    rec = tp / (tp + fn) if tp + fn else 1.0
    out = {"rmse_min": float(np.sqrt(np.mean((t - m["median"]) ** 2))), "precision_min": prec, "recall_min": rec,
           "f1_min": (2 * prec * rec / (prec + rec)) if prec + rec else 0.0, "n_true_viol_min": int(tv.sum()),
           "coverage90_min": float(((t >= m["lo90"]) & (t <= m["hi90"])).mean())}
    if "lo50" in pmin:
        out.update(_coverage_levels(sc.truth_daily_min, pmin, mask, "_min"))
    if p_night is not None:
        tn, mn = sc.truth_by_hour.loc[NIGHT_HOUR, mask], p_night.loc[mask]
        tvn, pvn = tn < THRESHOLD, mn["p_below"] > 0.5
        out["rmse_night"] = float(np.sqrt(np.mean((tn - mn["median"]) ** 2)))
        out["recall_night"] = float((tvn & pvn).sum() / tvn.sum()) if tvn.sum() else 1.0
        out["coverage90_night"] = float(((tn >= mn["lo90"]) & (tn <= mn["hi90"])).mean())
    return out


def run_scenario_time(sc, X, n_seed=3, n_max=15, noise_sd=0.03, seed=0, cache_dir="outputs/cache"):
    """Task 1: samples are (junction, hour) with hour in the operator's 07:00-17:00 window; the model
    predicts all 24 h and the daily minimum.  Baselines: the time-blind iteration-2 model (every sample
    treated as a 14:00 sample) and the mean of samples."""
    rng = np.random.default_rng(2000 + seed)

    def observe(js, hs):
        t = np.array([sc.truth_by_hour.loc[h, j] for j, h in zip(js, hs)])
        return np.clip(t + rng.normal(0, noise_sd, len(js)), 0.01, None)

    seed_js = list(rng.choice(sc.junctions, n_seed, replace=False))
    seed_hs = [int(h) for h in rng.choice(DAY_HOURS, n_seed)]
    seed_y = list(observe(seed_js, seed_hs))
    rows, snapshots = [], {}
    for strat in TIME_STRATEGIES:
        S = pd.DataFrame({"junction": seed_js, "hour": seed_hs, "y": seed_y})
        for n in range(n_seed, n_max + 1):
            unsampled = [j for j in sc.junctions if j not in set(S.junction)]
            m = SimGP24(sc, X, seed=seed, cache_dir=cache_dir).fit(S)
            hourly = m.predict_hours()
            pmin = m.predict_daily_min()
            pnight = m.predict_hour(NIGHT_HOUR, hourly); pnight["p_below"] = SimGP24.p_below(pnight)
            r = {"model": "simgp24", "strategy": strat, "n": n, "seed": seed, **_metrics_time(sc, pmin, pnight, unsampled)}
            r["map_kb"], r["map_kw"], r["map_gamma"], r["map_demand"], r["map_rough"] = m.map_params_
            r["map_dose"] = m.map_dose_
            rows.append(r)
            if strat == "random":
                m0 = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid="decay", lik_sd=0.25, doses=[1.0],
                             smooth_hours=False).fit(S)   # the task-1 model, for the before/after comparison
                rows.append({"model": "simgp24_g75", "strategy": strat, "n": n, "seed": seed,
                             **_metrics_time(sc, m0.predict_daily_min(), None, unsampled)})
                tb = SimGP(sc, seed=seed, cache_dir=cache_dir, lik="t").fit(X.loc[S.junction], S.y.values).predict(X)
                tb["p_below"] = SimGP.p_below(tb)          # its 14:00 flags, scored against the daily minimum
                rows.append({"model": "simgp_timeblind", "strategy": strat, "n": n, "seed": seed,
                             **_metrics_time(sc, tb, tb, unsampled)})
                mu = float(np.mean(S.y))
                bm = pd.DataFrame({"median": mu, "lo90": mu, "hi90": mu, "p_below": float(mu < THRESHOLD)}, index=sc.junctions)
                rows.append({"model": "mean_of_samples", "strategy": strat, "n": n, "seed": seed,
                             **_metrics_time(sc, bm, bm, unsampled)})
            if strat == TIME_MAIN and n in (n_seed, 8, n_max):
                snapshots[n] = {"samples": S.copy(), "pmin": pmin.copy(), "pnight": pnight.copy(),
                                "hourly": (hourly[0].copy(), hourly[1].copy())}
            if n == n_max:
                break
            cols = sc.junctions
            hourly_df = (pd.DataFrame(hourly[0], columns=cols), pd.DataFrame(hourly[1], columns=cols),
                         pd.DataFrame(m.z_sd_acq_, columns=cols))
            j, h = acquire_time(strat, hourly_df, pmin, unsampled, DAY_HOURS, rng, THRESHOLD)
            S = pd.concat([S, pd.DataFrame({"junction": [j], "hour": [h], "y": observe([j], [h])})], ignore_index=True)
    return pd.DataFrame(rows), snapshots


# ----------------------------------------------------------------------------- routes (task 4)
def run_routes(sc, X, seed=0, noise_sd=0.03, Ks=ROUTE_KS, cache_dir="outputs/cache"):
    """Task 4: a fixed monthly route of K sites chosen at once from the .inp alone (prior model), against
    K random sites and the K highest-demand sites.  Scored on the daily minimum at unsampled junctions."""
    rng = np.random.default_rng(3000 + seed)
    prior = SimGP24(sc, X, seed=seed, cache_dir=cache_dir).fit_prior()
    rows, routes = [], {}
    for K in Ks:
        for name, route in (("optimised", plan_route(prior, K)), ("highest_demand", demand_route(sc, K)),
                            ("random", random_route(sc, K, rng))):
            t = np.array([sc.truth_by_hour.loc[h, j] for j, h in zip(route.junction, route.hour)])
            S = pd.DataFrame({"junction": list(route.junction), "hour": [int(h) for h in route.hour],
                              "y": np.clip(t + rng.normal(0, noise_sd, K), 0.01, None)})
            m = SimGP24(sc, X, seed=seed, cache_dir=cache_dir).fit(S)
            uns = [j for j in sc.junctions if j not in set(S.junction)]
            hourly = m.predict_hours(); pnight = m.predict_hour(NIGHT_HOUR, hourly); pnight["p_below"] = SimGP24.p_below(pnight)
            pmin = m.predict_daily_min()
            # a daytime sample does not reveal a junction's night minimum, so the daily-minimum flags are
            # also scored over ALL junctions (the operator's view); the unsampled-only score stays primary
            allj = _metrics_time(sc, pmin, None, sc.junctions)
            rows.append({"route": name, "K": K, "seed": seed, **_metrics_time(sc, pmin, pnight, uns),
                         **{k + "_all": allj[k] for k in ("recall_min", "precision_min", "f1_min")}})
            routes[(name, K)] = {"route": route, "pmin": pmin}
    return pd.DataFrame(rows), routes


def plot_routes(df_r, df_t, sc, routes, out, K_map=8):
    df_r = df_r.copy(); df_r["missed"] = ((1 - df_r.recall_min) * df_r.n_true_viol_min).round(0)
    agg = df_r.groupby(["route", "K"])[["recall_min", "precision_min", "f1_min", "rmse_min", "coverage90_min", "recall_min_all", "f1_min_all"]].mean().reset_index()
    missed = df_r.groupby(["route", "K"]).missed.sum()
    n_months = df_r.seed.nunique()
    seq = df_t[(df_t.model == "simgp24") & (df_t.strategy == "straddle_min") & df_t.n.isin(ROUTE_KS)].groupby("n")[["recall_min", "f1_min"]].mean()
    style = {"optimised": ("tab:purple", "optimised route (from the .inp alone)"),
             "highest_demand": ("tab:orange", "K highest-demand sites (operator heuristic)"),
             "random": ("gray", "K random sites")}
    with plt.rc_context(BIG):
        fig = plt.figure(figsize=(22, 6.4)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25])
        ax0, ax1, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
        w = 0.26
        for i, (name, (c, lab)) in enumerate(style.items()):
            g = agg[agg.route == name].set_index("K")
            xs = np.arange(len(ROUTE_KS)) + (i - 1) * w
            ms = [missed[(name, k)] for k in ROUTE_KS]
            ax0.bar(xs, ms, w, color=c, label=lab)
            for x, v in zip(xs, ms):
                ax0.text(x, v + 0.15, f"{int(v)}", ha="center", fontsize=12, fontweight="bold")
            ax1.bar(xs, g.loc[list(ROUTE_KS), "f1_min"], w, color=c, label=lab)
        ax1.plot(np.arange(len(ROUTE_KS)), seq.loc[list(ROUTE_KS), "f1_min"], "k_", ms=22, mew=2, label="adaptive, one sample at a time (reference)")
        ax0.set(xticks=range(len(ROUTE_KS)), xticklabels=[f"K = {k} sites" for k in ROUTE_KS], ylabel="missed night violations",
                title=f"Night violations MISSED, summed over {n_months} scenario-months (lower is better)"); ax0.grid(alpha=0.3, axis="y")
        ax1.set(xticks=range(len(ROUTE_KS)), xticklabels=[f"K = {k} sites" for k in ROUTE_KS], ylim=(0, 1.05), ylabel="F1", title="F1 of the daily-minimum flags"); ax1.grid(alpha=0.3, axis="y")
        ax1.legend(loc="lower right", fontsize=10)
        r = routes[("optimised", K_map)]; route, pmin = r["route"], r["pmin"]
        wntr.graphics.plot_network(sc.wn, node_attribute=pmin["p_below"].to_dict(), node_size=55, node_cmap="Reds", node_range=(0, 1),
                                   ax=ax2, link_width=0.7, add_colorbar=True, title=f"The {K_map}-site route on scenario {sc.seed}: colour = P(daily min < {THRESHOLD}) before any sample")
        ax2.scatter([sc.coords[j][0] for j in route.junction], [sc.coords[j][1] for j in route.junction], s=170, facecolors="none", edgecolors="tab:purple", linewidths=2.2, zorder=5)
        for k, (j, h) in enumerate(zip(route.junction, route.hour)):
            ax2.annotate(f"{k + 1}: {h:02d}h", sc.coords[j], fontsize=10, color="tab:purple", xytext=(5, 4), textcoords="offset points", fontweight="bold")
        fig.suptitle(f"A 5-site route planned from the EPANET file alone misses {int(missed[('optimised', 5)])} night violation(s) in {n_months} scenario-months; "
                     f"the 5 highest-demand sites miss {int(missed[('highest_demand', 5)])}, 5 random sites {int(missed[('random', 5)])}", fontsize=15)
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)
    return agg


# ----------------------------------------------------------------------------- figures
def plot_maps(sc, snap, out, n):
    truth, pred, pv, sampled = sc.truth_snapshot, snap["pred"], snap["p_viol"], snap["sampled"]
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.4))
    panels = [("True chlorine (hidden from model)", truth, "viridis", (0, 1.2)),
              (f"Predicted median from {n} grab samples", pred["median"], "viridis", (0, 1.2)),
              ("Uncertainty: 90% band width (mg/L)", pred["hi90"] - pred["lo90"], "magma", None),
              (f"P(residual < {THRESHOLD} mg/L)", pv, "Reds", (0, 1))]
    for ax, (title, series, cmap, rng_) in zip(axes, panels):
        kw = dict(node_attribute=series.to_dict(), node_size=48, node_cmap=cmap, ax=ax, title=title,
                  link_width=0.6, add_colorbar=True)
        if rng_:
            kw["node_range"] = rng_
        wntr.graphics.plot_network(sc.wn, **kw)
        ax.scatter([sc.coords[j][0] for j in sampled], [sc.coords[j][1] for j in sampled], s=120,
                   facecolors="none", edgecolors="cyan", linewidths=1.8, zorder=5, label="grab samples")
    axes[1].legend(loc="lower left", fontsize=9)
    fig.suptitle(f"{sc.wn_name} (EPANET example based on North Marin Water District, CA) — "
                 f"{sc.sample_hour}:00 snapshot, scenario {sc.seed}, calibrated-simulator GP", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def plot_day_night(sc, out):
    frac = (sc.truth_by_hour < THRESHOLD).mean(axis=1)
    fig = plt.figure(figsize=(18, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1])
    ax0 = fig.add_subplot(gs[0])
    ax0.bar(frac.index, frac.values * 100, color=["tab:red" if h >= 18 or h < 5 else "tab:blue" for h in frac.index])
    ax0.axvspan(6.5, 16.5, color="gold", alpha=0.15, label="typical sampling window")
    ax0.set(xlabel="hour of day", ylabel=f"% of junctions below {THRESHOLD} mg/L",
            title="The network is worst when nobody samples")
    ax0.legend(loc="upper center")
    for ax, h in [(fig.add_subplot(gs[1]), 14), (fig.add_subplot(gs[2]), 22)]:
        wntr.graphics.plot_network(sc.wn, node_attribute=sc.truth_by_hour.loc[h].to_dict(), node_size=48,
                                   node_cmap="viridis", node_range=(0, 1.2), ax=ax, link_width=0.6,
                                   add_colorbar=True, title=f"True chlorine at {h}:00 "
                                   f"({int((sc.truth_by_hour.loc[h] < THRESHOLD).sum())} junctions below {THRESHOLD})")
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def plot_day_night_predicted(sc, snap, out, n):
    """The pitch in one figure: the true night map next to what n DAYTIME samples predict for it,
    and the daily-minimum violations next to the model's P(daily min < 0.2)."""
    S, pmin, pnight = snap["samples"], snap["pmin"], snap["pnight"]
    uns = [j for j in sc.junctions if j not in set(S.junction)]
    met = _metrics_time(sc, pmin, pnight, uns)
    tv = sc.truth_daily_min < THRESHOLD
    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 4, figsize=(24, 6.2))
        panels = [(f"TRUE chlorine at {NIGHT_HOUR}:00 (hidden)\n{int((sc.truth_by_hour.loc[NIGHT_HOUR] < THRESHOLD).sum())} junctions below {THRESHOLD} mg/L",
                   sc.truth_by_hour.loc[NIGHT_HOUR], "viridis", (0, 1.2)),
                  (f"PREDICTED {NIGHT_HOUR}:00 from {n} daytime samples\nrecall {met['recall_night']:.2f}, RMSE {met['rmse_night']:.2f} mg/L",
                   pnight["median"], "viridis", (0, 1.2)),
                  (f"TRUE daily minimum (hidden)\n{int(tv.sum())} junctions below {THRESHOLD} mg/L at some hour",
                   sc.truth_daily_min, "viridis", (0, 1.2)),
                  (f"P(daily minimum < {THRESHOLD} mg/L) from daytime samples\nrecall {met['recall_min']:.2f} on unsampled junctions",
                   pmin["p_below"], "Reds", (0, 1))]
        for ax, (title, series, cmap, rng_) in zip(axes, panels):
            wntr.graphics.plot_network(sc.wn, node_attribute=series.to_dict(), node_size=60, node_cmap=cmap,
                                       node_range=rng_, ax=ax, title=title, link_width=0.7, add_colorbar=True)
        ax = axes[1]
        ax.scatter([sc.coords[j][0] for j in S.junction], [sc.coords[j][1] for j in S.junction], s=150,
                   facecolors="none", edgecolors="cyan", linewidths=2.0, zorder=5, label="grab samples (hour)")
        for j, h in zip(S.junction, S.hour):
            ax.annotate(f"{h}h", sc.coords[j], fontsize=9, color="cyan", xytext=(4, 4), textcoords="offset points")
        ax.legend(loc="lower left")
        fig.suptitle(f"{sc.wn_name}: {n} grab samples taken between 07:00 and 17:00 predict the night — "
                     f"scenario {sc.seed}, time-aware calibrated-simulator GP ({TIME_MAIN} rule)", fontsize=16)
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def plot_curves_time(df, out):
    """Daily-minimum error, recall and coverage vs number of daytime samples."""
    agg = df.groupby(["model", "strategy", "n"]).agg(rmse=("rmse_min", "mean"), recall=("recall_min", "mean"),
                                                     cov90=("coverage90_min", "mean"), recall_night=("recall_night", "mean")).reset_index()
    style = {("simgp24", "random"): ("tab:red", "--", "time-aware GP, random daytime samples"),
             ("simgp24", "uncertainty"): ("tab:red", ":", "time-aware GP, max-uncertainty rule"),
             ("simgp24", "straddle"): ("tab:red", "-", "time-aware GP, straddle rule"),
             ("simgp24", "straddle_min"): ("tab:purple", "-", "time-aware GP, straddle on daily min"),
             ("simgp24_g75", "random"): ("tab:pink", "--", "time-aware GP, decay-only grid, random"),
             ("simgp_timeblind", "random"): ("tab:orange", "-", "time-blind GP (iteration 2, every sample = 14:00)"),
             ("mean_of_samples", "random"): ("gray", "--", "mean of samples (today's practice)")}
    n15 = agg[(agg.n == agg.n.max()) & (agg.model == "simgp24") & (agg.strategy == TIME_MAIN)].iloc[0]
    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 3, figsize=(21, 5.6))
        for (mname, strat), g in agg.groupby(["model", "strategy"]):
            if (mname, strat) not in style:
                continue
            c, ls, lab = style[(mname, strat)]
            axes[0].plot(g.n, g.rmse, ls, color=c, marker="o", ms=4, label=lab)
            axes[1].plot(g.n, g.recall, ls, color=c, marker="o", ms=4, label=lab)
            if mname.startswith("simgp24"):
                axes[2].plot(g.n, g.cov90, ls, color=c, marker="o", ms=4, label=lab)
        axes[1].axhline(0.8, color="k", lw=0.8, ls=":"); axes[1].text(3.1, 0.81, "target 0.80", fontsize=10)
        axes[2].axhline(0.9, color="k", lw=0.8, ls=":"); axes[2].text(3.1, 0.905, "target 0.90", fontsize=10)
        axes[0].set(title="Daily-minimum error at unsampled junctions", xlabel="daytime grab samples", ylabel="RMSE of daily minimum (mg/L)")
        axes[1].set(title=f"Recall of junctions whose daily minimum < {THRESHOLD} mg/L", xlabel="daytime grab samples", ylabel="recall", ylim=(0, 1.02))
        axes[2].set(title="Calibration: truth inside the 90% band (daily minimum)", xlabel="daytime grab samples", ylabel="coverage", ylim=(0, 1.02))
        for ax in axes:
            ax.grid(alpha=0.3)
        axes[1].legend(loc="lower right", fontsize=9)
        fig.suptitle(f"Daytime samples predict the night: recall {n15.recall:.2f} of daily-minimum violations at "
                     f"{int(n15.n)} daytime samples, {TIME_MAIN} rule (time-blind model: "
                     f"{agg[(agg.model == 'simgp_timeblind') & (agg.n == n15.n)].recall.iloc[0]:.2f})", fontsize=15)
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)
    return agg


def plot_reliability(df, df_t, out):
    """Nominal vs empirical coverage at 50/80/90/95%, averaged over seeds and n (random samples)."""
    def curve(d, suffix=""):
        return [d[f"coverage{q}{suffix}"].mean() for q in LEVELS]
    snap = df[df.strategy == "random"]; tm = df_t[df_t.strategy == "random"]
    series = [("14:00 snapshot — iteration 2 (decay grid, 75 runs)", curve(snap[snap.model == "simgp_core_g75"]), "tab:orange", "--"),
              ("14:00 snapshot — iteration 3 (+ demand, roughness, dose axes; local hydraulic error)", curve(snap[snap.model == "simgp_core"]), "tab:red", "-"),
              ("daily minimum — task 1 (decay grid)", curve(tm[tm.model == "simgp24_g75"], "_min"), "tab:blue", "--"),
              ("daily minimum — iteration 3", curve(tm[tm.model == "simgp24"], "_min"), "tab:purple", "-")]
    before, after = series[0][1][2], series[1][1][2]
    with plt.rc_context(BIG):
        fig, ax = plt.subplots(figsize=(8.5, 7.5))
        ax.plot([0.4, 1.0], [0.4, 1.0], "k:", lw=1.2, label="perfect calibration")
        ax.axhspan(0.88, 0.95, xmin=0.78, xmax=0.86, color="green", alpha=0.12)
        for lab, ys, c, ls in series:
            ax.plot([q / 100 for q in LEVELS], ys, ls, color=c, marker="o", ms=7, lw=2, label=lab)
        ax.set(xlabel="nominal band level", ylabel="empirical coverage (unsampled junctions)", xlim=(0.4, 1.0), ylim=(0.2, 1.0),
               title=f"90% band now holds the truth {after:.2f} of the time (was {before:.2f}); target 0.88-0.95")
        ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=10)
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def plot_stress(outdir, structural_dir, out, network="Net3"):
    """Task 3: the same model scored against a truth with structural errors the grid cannot represent."""
    a, b = pd.read_csv(os.path.join(outdir, f"results_{network}.csv")), pd.read_csv(os.path.join(structural_dir, f"results_{network}.csv"))
    at, bt = pd.read_csv(os.path.join(outdir, f"results_time_{network}.csv")), pd.read_csv(os.path.join(structural_dir, f"results_time_{network}.csv"))
    def g(d, model, strat, cols):
        return d[(d.model == model) & (d.strategy == strat)].groupby("n")[cols].mean()
    series = [("14:00 map, GP, straddle rule", g(a, "simgp_core", "straddle", ["rmse", "recall", "coverage90"]), g(b, "simgp_core", "straddle", ["rmse", "recall", "coverage90"]), "tab:red"),
              ("14:00 map, GP, random", g(a, "simgp_core", "random", ["rmse", "recall", "coverage90"]), g(b, "simgp_core", "random", ["rmse", "recall", "coverage90"]), "tab:orange"),
              ("14:00 map, simulator alone (no GP), random", g(a, "simgp_core_nogp", "random", ["rmse", "recall", "coverage90"]), g(b, "simgp_core_nogp", "random", ["rmse", "recall", "coverage90"]), "tab:gray")]
    tm = [("daily minimum, straddle on daily min", g(at, "simgp24", "straddle_min", ["rmse_min", "recall_min", "coverage90_min"]), g(bt, "simgp24", "straddle_min", ["rmse_min", "recall_min", "coverage90_min"]), "tab:purple")]
    r15 = (g(a, "simgp_core", "straddle", ["recall"]).loc[15, "recall"], g(b, "simgp_core", "straddle", ["recall"]).loc[15, "recall"])
    c = (a[(a.model == "simgp_core") & (a.strategy == "random")].coverage90.mean(), b[(b.model == "simgp_core") & (b.strategy == "random")].coverage90.mean())
    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 3, figsize=(21, 5.8))
        for lab, x, y, col in series + tm:
            for ax, k in zip(axes, range(3)):
                ax.plot(x.index, x.iloc[:, k], "-", color=col, marker="o", ms=4, label=f"{lab} — file correct")
                ax.plot(y.index, y.iloc[:, k], "--", color=col, marker="s", ms=4, label=f"{lab} — file wrong")
        axes[0].set(title="Error at unsampled junctions", xlabel="grab samples", ylabel="RMSE (mg/L)")
        axes[1].set(title=f"Recall of junctions below {THRESHOLD} mg/L", xlabel="grab samples", ylabel="recall", ylim=(0.4, 1.02))
        axes[1].axhline(0.7, color="k", lw=0.8, ls=":"); axes[1].text(3.1, 0.71, "stop line 0.70", fontsize=10)
        axes[2].set(title="Truth inside the 90% band", xlabel="grab samples", ylabel="coverage", ylim=(0.4, 1.02))
        axes[2].axhline(0.9, color="k", lw=0.8, ls=":")
        for ax in axes:
            ax.grid(alpha=0.3)
        axes[2].legend(fontsize=8, loc="lower left")
        info = json.load(open(os.path.join(structural_dir, f"summary_{network}.json"))).get("structural_noise", {})
        n_closed, n_seeds = sum(1 for v in info.values() if v.get("closed_pipe")), max(len(info), 1)
        fig.suptitle(f"Stress test: a tank with half the volume the operator's file says ({n_seeds} of {n_seeds} scenarios) and a closed pipe "
                     f"({n_closed} of {n_seeds}) — straddle-rule recall at 15 samples {r15[0]:.2f} → {r15[1]:.2f}, 90% coverage {c[0]:.2f} → {c[1]:.2f}", fontsize=13)
        fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def plot_curves(df, out):
    d = df[df.strategy == "random"]
    agg = d.groupby(["model", "n"]).agg(rmse=("rmse", "mean"), f1=("f1", "mean"), recall=("recall", "mean"),
                                        cov90=("coverage90", "mean")).reset_index()
    style = {"simgp_core": ("tab:red", "-", "calibrated-simulator GP (core)"),
             "simgp_core_nogp": ("tab:red", "--", "calibrated simulator alone (no discrepancy GP)"),
             "simgp_core_g75": ("tab:pink", "-", "calibrated-simulator GP (core, decay-only grid)"),
             "simgp_rich": ("tab:red", ":", "calibrated-simulator GP (rich)"),
             "physgp_core": ("tab:green", "-", "decay-law GP (core)"),
             "physgp_rich": ("tab:green", ":", "decay-law GP (rich)"),
             "pinn_rich": ("tab:brown", "-.", "graph-PINN (rich)"),
             "mean_of_samples": ("gray", "--", "mean of samples (today)"),
             "nearest_sample": ("tab:orange", "--", "nearest sampled node"),
             "decay_only": ("tab:purple", "--", "decay law only")}
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.8))
    for m, g in agg.groupby("model"):
        c, ls, lab = style[m]
        axes[0].plot(g.n, g.rmse, ls, color=c, marker="o", ms=3, label=lab)
        axes[1].plot(g.n, g.recall, ls, color=c, marker="o", ms=3, label=lab)
        if g.cov90.notna().any():
            axes[2].plot(g.n, g.cov90, ls, color=c, marker="o", ms=3, label=lab)
    axes[2].axhline(0.9, color="k", lw=0.8, ls=":"); axes[2].text(3.1, 0.905, "target 0.90", fontsize=8)
    axes[0].set(title="Q1/Q2  error at UNSAMPLED junctions (random samples)", xlabel="grab samples", ylabel="RMSE (mg/L)")
    axes[1].set(title=f"recall of junctions below {THRESHOLD} mg/L", xlabel="grab samples", ylabel="recall")
    axes[2].set(title="calibration: fraction of truth inside the 90% band", xlabel="grab samples", ylabel="coverage", ylim=(0, 1))
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)

    s = df[df.model == "simgp_core"]
    agg2 = s.groupby(["strategy", "n"]).agg(rmse=("rmse", "mean"), recall=("recall", "mean"), f1=("f1", "mean"),
                                             recall_all=("recall_all", "mean"), cov90=("coverage90", "mean")).reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.6))
    for st, g in agg2.groupby("strategy"):
        for ax, col in zip(axes, ["rmse", "recall", "f1"]):
            ax.plot(g.n, g[col], marker="o", ms=3, label=st)
    axes[0].set(title="Q3  sampling rule vs error", xlabel="grab samples", ylabel="RMSE (mg/L)")
    axes[1].set(title=f"Q3  sampling rule vs recall of nodes < {THRESHOLD}", xlabel="grab samples", ylabel="recall")
    axes[2].set(title="Q3  sampling rule vs F1", xlabel="grab samples", ylabel="F1")
    for ax in axes:
        ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(out.replace("curves_", "strategies_"), dpi=130); plt.close(fig)
    return agg, agg2


def main(network="Net3", seeds=tuple(range(8)), n_seed=3, n_max=15, outdir="outputs", sample_hour=14,
         structural_noise=False, cache=None):
    """structural_noise=True is the task-3 stress test: the truth gets a closed pipe / low tank the
    operator's file does not have; outputs go to <outdir>/structural, the grid cache is shared."""
    cache = cache or os.path.join(outdir, "cache")
    if structural_noise:
        outdir = os.path.join(outdir, "structural" if structural_noise is True else f"structural_{structural_noise}")
    os.makedirs(outdir, exist_ok=True)
    frames, frames_t, frames_r, summary = [], [], [], {}
    for s in seeds:
        sc = build_scenario(network, seed=s, sample_hour=sample_hour, structural_noise=structural_noise)
        if sc.structural:
            summary.setdefault("structural_noise", {})[str(s)] = sc.structural
            print(f"seed {s}: structural noise -> {sc.structural}", flush=True)
        df, snaps, X = run_scenario(sc, n_seed=n_seed, n_max=n_max, seed=s, cache_dir=cache)
        frames.append(df)
        df_t, snaps_t = run_scenario_time(sc, X, n_seed=n_seed, n_max=n_max, seed=s, cache_dir=cache)
        frames_t.append(df_t)
        df_r, routes = run_routes(sc, X, seed=s, cache_dir=cache)
        frames_r.append(df_r)
        if s == seeds[0]:
            sc0, routes0 = sc, routes
            for n, snap in snaps.items():
                plot_maps(sc, snap, os.path.join(outdir, f"map_{network}_n{n}.png"), n)
            plot_day_night(sc, os.path.join(outdir, f"day_vs_night_{network}.png"))
            plot_day_night_predicted(sc, snaps_t[n_max], os.path.join(outdir, f"day_vs_night_predicted_{network}.png"), n_max)
            X.to_csv(os.path.join(outdir, f"features_{network}_seed{s}.csv"))
            summary["scenario0"] = {
                "junctions": len(sc.junctions),
                "below_threshold_at_sampling_hour": int((sc.truth_snapshot < THRESHOLD).sum()),
                "below_threshold_daily_min": int((sc.truth_daily_min < THRESHOLD).sum()),
                "pct_below_by_hour": {int(h): round(float(v) * 100, 1) for h, v in (sc.truth_by_hour < THRESHOLD).mean(axis=1).items()},
            }
        print(f"seed {s}: {int((sc.truth_snapshot < THRESHOLD).sum())}/{len(sc.junctions)} junctions below "
              f"{THRESHOLD} mg/L at {sample_hour}:00", flush=True)
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(os.path.join(outdir, f"results_{network}.csv"), index=False)
    agg, agg2 = plot_curves(df, os.path.join(outdir, f"curves_{network}.png"))
    summary["models_random_sampling"] = {
        str(n): agg[agg.n == n].set_index("model")[["rmse", "recall", "f1", "cov90"]].round(3).to_dict("index")
        for n in (n_seed, 8, n_max)}
    summary["strategies_simgp"] = {
        str(n): agg2[agg2.n == n].set_index("strategy")[["rmse", "recall", "f1", "recall_all", "cov90"]].round(3).to_dict("index")
        for n in (n_seed, 8, n_max)}
    mp = df[(df.model == "simgp_core") & (df.strategy == "random")].groupby("n")[["map_kb", "map_kw", "map_gamma"]].mean()
    summary["map_params_mean_by_n"] = {str(n): mp.loc[n].round(2).to_dict() for n in (n_seed, 8, n_max)}
    df_t = pd.concat(frames_t, ignore_index=True)
    df_t.to_csv(os.path.join(outdir, f"results_time_{network}.csv"), index=False)
    agg_t = plot_curves_time(df_t, os.path.join(outdir, f"curves_time_{network}.png"))
    plot_reliability(df, df_t, os.path.join(outdir, f"reliability_{network}.png"))
    if structural_noise:
        base = os.path.dirname(outdir)
        if os.path.exists(os.path.join(base, f"results_{network}.csv")):
            plot_stress(base, outdir, os.path.join(base, f"stress_test_{network}.png"), network)
    agg_t["key"] = agg_t.model + "/" + agg_t.strategy
    summary["time_aware_daily_min"] = {
        str(n): agg_t[agg_t.n == n].set_index("key")[["rmse", "recall", "cov90", "recall_night"]].round(3).to_dict("index")
        for n in (n_seed, 8, n_max)}
    df_r = pd.concat(frames_r, ignore_index=True)
    df_r.to_csv(os.path.join(outdir, f"results_routes_{network}.csv"), index=False)
    agg_r = plot_routes(df_r, df_t, sc0, routes0, os.path.join(outdir, f"route_comparison_{network}.png"))
    summary["routes_daily_min"] = {
        str(K): agg_r[agg_r.K == K].set_index("route")[["recall_min", "precision_min", "f1_min", "rmse_min", "coverage90_min", "recall_min_all", "f1_min_all"]].round(3).to_dict("index")
        for K in ROUTE_KS}
    summary["route_scenario0_K8"] = routes0[("optimised", 8)]["route"].to_dict("records")
    with open(os.path.join(outdir, f"summary_{network}.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return df, agg, agg2, summary, agg_t


if __name__ == "__main__":
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    net = args[0] if args else "Net3"
    seeds = tuple(range(int(args[1]))) if len(args) > 1 else tuple(range(8))
    structural = "persistent" if "--structural=persistent" in sys.argv else ("--structural" in sys.argv)
    _, agg, agg2, _, agg_t = main(net, seeds=seeds, structural_noise=structural)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30)
    for col in ("rmse", "recall", "cov90"):
        print(f"\n== models under random sampling: {col}")
        print(agg[agg.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="model", values=col).round(3))
    print("\n== sampling rules with calibrated-simulator GP: recall on unsampled junctions")
    print(agg2[agg2.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="strategy", values="recall").round(3))
    print("\n== sampling rules: recall counting violations found by the sample itself")
    print(agg2[agg2.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="strategy", values="recall_all").round(3))
    for col in ("rmse", "recall", "cov90", "recall_night"):
        print(f"\n== time-aware: DAILY MINIMUM from daytime samples, {col}")
        print(agg_t[agg_t.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="key", values=col).round(3))
    print("\n== routes: K sites chosen at once, daily-minimum recall / precision / F1")
    print(pd.read_csv(os.path.join("outputs" if not structural else f"outputs/structural{'_' + structural if isinstance(structural, str) else ''}", f"results_routes_{net}.csv"))
          .groupby(["K", "route"])[["recall_min", "precision_min", "f1_min", "rmse_min", "coverage90_min", "recall_min_all", "f1_min_all"]].mean().round(3))
