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
from .simgp import SimGP
from .simulate import build_scenario
from .surrogate import PhysicsGP, acquire, baseline_decay_only, baseline_mean, baseline_nearest

warnings.filterwarnings("ignore")  # GP optimizer bound warnings are expected with 3-15 points
STRATEGIES = ["random", "uncertainty", "straddle"]
THRESHOLD = 0.2
PINN_AT = (3, 5, 8, 10, 12, 15)


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


def _row(model_name, strategy, n, seed, truth, pred, unsampled, model_cls):
    pv = model_cls.p_below(pred, THRESHOLD)
    return {"model": model_name, "strategy": strategy, "n": n, "seed": seed,
            **_metrics(truth, pred["median"], pv, pred["lo90"], pred["hi90"], unsampled)}, pv


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
            r["map_kb"], r["map_kw"], r["map_gamma"] = main.map_params_
            rows.append(r)
            if strat == "random":
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


def plot_curves(df, out):
    d = df[df.strategy == "random"]
    agg = d.groupby(["model", "n"]).agg(rmse=("rmse", "mean"), f1=("f1", "mean"), recall=("recall", "mean"),
                                        cov90=("coverage90", "mean")).reset_index()
    style = {"simgp_core": ("tab:red", "-", "calibrated-simulator GP (core)"),
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
    agg2 = s.groupby(["strategy", "n"]).agg(rmse=("rmse", "mean"), recall=("recall", "mean"), f1=("f1", "mean")).reset_index()
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


def main(network="Net3", seeds=tuple(range(8)), n_seed=3, n_max=15, outdir="outputs", sample_hour=14):
    os.makedirs(outdir, exist_ok=True)
    cache = os.path.join(outdir, "cache")
    frames, summary = [], {}
    for s in seeds:
        sc = build_scenario(network, seed=s, sample_hour=sample_hour)
        df, snaps, X = run_scenario(sc, n_seed=n_seed, n_max=n_max, seed=s, cache_dir=cache)
        frames.append(df)
        if s == seeds[0]:
            for n, snap in snaps.items():
                plot_maps(sc, snap, os.path.join(outdir, f"map_{network}_n{n}.png"), n)
            plot_day_night(sc, os.path.join(outdir, f"day_vs_night_{network}.png"))
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
        str(n): agg2[agg2.n == n].set_index("strategy")[["rmse", "recall", "f1"]].round(3).to_dict("index")
        for n in (n_seed, 8, n_max)}
    with open(os.path.join(outdir, f"summary_{network}.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return df, agg, agg2, summary


if __name__ == "__main__":
    import sys
    net = sys.argv[1] if len(sys.argv) > 1 else "Net3"
    seeds = tuple(range(int(sys.argv[2]))) if len(sys.argv) > 2 else tuple(range(8))
    _, agg, agg2, _ = main(net, seeds=seeds)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30)
    for col in ("rmse", "recall", "cov90"):
        print(f"\n== models under random sampling: {col}")
        print(agg[agg.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="model", values=col).round(3))
    print("\n== sampling rules with calibrated-simulator GP: recall")
    print(agg2[agg2.n.isin([3, 5, 8, 10, 12, 15])].pivot(index="n", columns="strategy", values="recall").round(3))
