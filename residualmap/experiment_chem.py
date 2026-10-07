"""
experiment_chem.py: the iteration-4 chemistry experiments (journal tasks 9 to 13), one subcommand per task.

    python -m residualmap.experiment_chem age Net3 8     # task 9: water age and where chlorine is lost
    python -m residualmap.experiment_chem age Net2 8
    python -m residualmap.experiment_chem age ky4 2

Everything is simulation: the "truth" is a hidden perturbed copy of the same EPANET file (simulate.build_scenario).

age (task 9).  For each truth (the default one, and the persistent structural stress test: a tank with 0.7 x the
file's diameter and, in about half the scenarios, a closed pipe) and each seed:
  * water age: the operator's nominal daily-mean age against the truth's own (build_scenario(truth_age=True)),
    RMSE / MAE / bias in hours; coverage of the 9-member range of hydraulic_age_band (does the truth's daily-mean age
    lie between the smallest and largest member's?); the same for the posterior-weighted age and its central 90%
    after 8 random daytime samples.  Junctions whose TRUE daily-max age exceeds age.AGE_UNCONVERGED_H (160 h) are
    excluded from the age metrics and counted (the plan's rule).  That rule misses most unconverged ages: every run
    starts with the file's water in the pipes and tanks, and where some of it still arrives on the scored day the
    age is a lower bound.  So the same metrics are also given on the junctions where both the operator's model and
    the truth hold at most age.INITIAL_SHARE_MAX of that starting water (simulate.initial_water_share), as
    'converged_only' in the JSON and conv_* in the CSV, and outputs/chem/water_age_junctions_<net>.csv gives the
    share per junction of the operator's file.
  * chlorine by age: today's SimGP24 at the same 8 samples, scored on the unsampled junctions' daily minimum
    (RMSE, MAE, bias, 90% coverage, recall below the threshold with P > 0.5, false alarms), by nominal-age tercile in
    the JSON and by decile in outputs/chem/error_by_age_<net>.csv.  This is the baseline task 13 compares against.
  * where chlorine is lost: age.loss_split for the model's MAP member against age.truth_loss_split on the truth.

The 8 samples are drawn exactly as the app's demo draws them (default_rng(seed): 8 distinct junctions, 8 hours in
07:00-17:00, N(0, 0.03) mg/L noise, clipped at 0.01 and rounded to 0.01 mg/L like a field reading), so Net3 seed 0
under the default truth is the app's default demo.
"""
from __future__ import annotations

import argparse
import inspect
import json
import os
import time
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wntr

from .age import (AGE_UNCONVERGED_H, BAND_COVERAGE_BAR, INITIAL_SHARE_MAX, MIN_LOSS, RANGE_LABEL, band_label,
                  hydraulic_age_band, loss_split, oldest_water, posterior_age, truth_loss_split)
from .experiment import BIG, NET_TRUTH
from .features import build_features
from .simgp import DAY_HOURS, SimGP24
from .simulate import build_scenario

warnings.filterwarnings("ignore")   # GP optimiser bound warnings, as in experiment.py
THRESHOLD = 0.2
N_SAMPLES = 8
NOISE_SD = 0.03
AGE_TRUTHS = {"default": False, "persistent": "persistent"}
TRUTH_DRAWS = ("per scenario: bulk rate x U(0.8, 1.2); per pipe: wall rate x roughness factor (gamma = 1) x lognormal(0, 0.4) "
               "and Hazen-Williams C x lognormal(0, 0.10); demand x U(0.85, 1.15) globally and x lognormal(0, 0.15) per "
               "junction; dose x U(0.9, 1.1) per source (simulate.build_scenario)")


def truth_settings(net: str) -> dict:
    """The truth's nominal rates and dose for this network: build_scenario's defaults with experiment.NET_TRUTH's
    overrides, so every number the docs quote about the truth is in the output file."""
    sig = inspect.signature(build_scenario).parameters
    out = {k: sig[k].default for k in ("source_dose", "kb_per_day", "kw_m_per_day")}
    out.update(NET_TRUTH.get(net, {}))
    return {**out, "hidden_draws": TRUTH_DRAWS}


def demo_samples(sc, seed: int, n: int = N_SAMPLES, noise_sd: float = NOISE_SD) -> pd.DataFrame:
    """n random daytime samples, drawn exactly as app.py's demo draws them (so Net3 seed 0 is the app's default demo):
    default_rng(seed) picks n distinct junctions, then n hours from DAY_HOURS, then one N(0, noise_sd) error per
    sample; readings are clipped at 0.01 mg/L and rounded to 0.01 mg/L."""
    rng = np.random.default_rng(seed)
    js = list(rng.choice(sc.junctions, n, replace=False)) if n else []
    hs = [int(h) for h in rng.choice(DAY_HOURS, n)] if n else []
    ys = [float(np.clip(sc.truth_by_hour.loc[h, j] + rng.normal(0, noise_sd), 0.01, None)) for j, h in zip(js, hs)]
    return pd.DataFrame({"junction": [str(j) for j in js], "hour": hs, "y": np.round(ys, 2)})


def _chlorine_metrics(t: pd.Series, med: pd.Series, lo: pd.Series, hi: pd.Series, flag: pd.Series) -> dict:
    """Daily-minimum scores on a set of junctions.  Recall is NaN when there is no violation to find (not 1.0)."""
    tv = t < THRESHOLD
    tp, fp = int((tv & flag).sum()), int((~tv & flag).sum())
    e = med - t
    return {"n": int(len(t)), "rmse": float(np.sqrt((e ** 2).mean())), "mae": float(e.abs().mean()),
            "bias": float(e.mean()), "coverage90": float(((t >= lo) & (t <= hi)).mean()),
            "n_true_viol": int(tv.sum()), "n_found": tp, "recall": float(tp / tv.sum()) if tv.sum() else float("nan"),
            "false_alarms": fp, "n_flagged": int(flag.sum())}


def _age_metrics(true_mean: pd.Series, est: pd.Series, lo: pd.Series, hi: pd.Series, inc: pd.Series) -> dict:
    e = (est - true_mean)[inc]
    return {"rmse_h": float(np.sqrt((e ** 2).mean())), "mae_h": float(e.abs().mean()), "bias_h": float(e.mean()),
            "coverage": float(((true_mean >= lo) & (true_mean <= hi))[inc].mean()),
            "width_mean_h": float((hi - lo)[inc].mean())}


def run_age_seed(net: str, seed: int, truth: str, band, cache_dir: str) -> tuple[dict, pd.DataFrame, dict]:
    """One (network, truth, seed): the per-seed row, the per-junction frame (unsampled junctions, for the age bins)
    and what the figure needs, including the age band.  band=None builds it from this scenario's nominal model: the
    band belongs to the operator's file, so one per network serves every seed and truth."""
    sc = build_scenario(net, seed=seed, structural_noise=AGE_TRUTHS[truth], truth_age=True, truth_loss_split=True,
                        **NET_TRUTH.get(net, {}))
    band = band if band is not None else hydraulic_age_band(sc)
    X = build_features(sc)
    # ---- water age: nominal and posterior against the truth's own
    ta = sc.truth_age_by_hour_h
    t_mean, t_max = ta.mean(), ta.max()
    inc = t_max <= AGE_UNCONVERGED_H
    nom = sc.age_daily_mean_h.loc[sc.junctions]
    lo, hi = band.daily_range("mean")
    lo_x, hi_x = band.daily_range("max")
    # where the run's starting water is still arriving, both ages are lower bounds: the converged subsets
    tsh_mean, tsh_max = sc.truth_initial_share.mean(), sc.truth_initial_share.max()
    nsh_mean, nsh_max = band.share("mean"), band.share("max")
    conv = inc & (nsh_mean <= INITIAL_SHARE_MAX) & (tsh_mean <= INITIAL_SHARE_MAX)
    conv_x = inc & (nsh_max <= INITIAL_SHARE_MAX) & (tsh_max <= INITIAL_SHARE_MAX)
    S = demo_samples(sc, seed)
    model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, threshold=THRESHOLD).fit(S)
    post = posterior_age(model, band)
    a_nom = _age_metrics(t_mean, nom, lo, hi, inc)
    a_post = _age_metrics(t_mean, post.daily_mean, post.lo90, post.hi90, inc)
    cov_max = float(((t_max >= lo_x) & (t_max <= hi_x))[inc].mean())
    c_nom = _age_metrics(t_mean, nom, lo, hi, conv)
    c_post = _age_metrics(t_mean, post.daily_mean, post.lo90, post.hi90, conv)
    c_cov_max = float(((t_max >= lo_x) & (t_max <= hi_x))[conv_x].mean()) if conv_x.any() else float("nan")
    # ---- chlorine at the daily minimum, unsampled junctions
    pmin = model.predict_daily_min()
    uns = [j for j in sc.junctions if j not in set(S.junction)]
    flag = pmin["p_below"] > 0.5
    chl = _chlorine_metrics(sc.truth_daily_min.loc[uns], pmin.loc[uns, "median"], pmin.loc[uns, "lo90"],
                            pmin.loc[uns, "hi90"], flag.loc[uns])
    # ---- where chlorine is lost: the MAP member against the truth
    ls, ts = loss_split(model), truth_loss_split(sc)
    both = ls.wall_share.notna() & ts.wall_share.notna()
    d = (ls.wall_share - ts.wall_share)[both]
    kb, kw, g, dm, rm = model.map_params_
    row = {"truth": truth, "seed": seed, "n_junctions": len(sc.junctions), "n_excluded": int((~inc).sum()),
           "rmse_h": a_nom["rmse_h"], "mae_h": a_nom["mae_h"], "bias_h": a_nom["bias_h"],
           "band_coverage": a_nom["coverage"], "band_coverage_daily_max": cov_max, "band_width_mean_h": a_nom["width_mean_h"],
           "post_rmse_h": a_post["rmse_h"], "post_mae_h": a_post["mae_h"], "post_bias_h": a_post["bias_h"],
           "post_band_coverage": a_post["coverage"], "post_band_width_mean_h": a_post["width_mean_h"],
           "true_age_mean_h": float(t_mean[inc].mean()), "nominal_age_mean_h": float(nom[inc].mean()),
           # lower bounds: junctions holding more than INITIAL_SHARE_MAX of the run's starting water (daily mean, or
           # at any hour for the daily max), and the age metrics on the junctions where neither side does
           "n_lower_bound_nominal": int((nsh_mean > INITIAL_SHARE_MAX).sum()),
           "n_lower_bound_truth": int((tsh_mean > INITIAL_SHARE_MAX).sum()),
           "n_lower_bound_nominal_daily_max": int((nsh_max > INITIAL_SHARE_MAX).sum()),
           "n_lower_bound_truth_daily_max": int((tsh_max > INITIAL_SHARE_MAX).sum()),
           "initial_share_nominal_median": float(nsh_mean.median()), "initial_share_truth_median": float(tsh_mean.median()),
           "n_converged": int(conv.sum()), "n_converged_daily_max": int(conv_x.sum()),
           "conv_rmse_h": c_nom["rmse_h"], "conv_mae_h": c_nom["mae_h"], "conv_bias_h": c_nom["bias_h"],
           "conv_band_coverage": c_nom["coverage"], "conv_band_coverage_daily_max": c_cov_max,
           "conv_post_rmse_h": c_post["rmse_h"], "conv_post_band_coverage": c_post["coverage"],
           **{f"{k}_min": v for k, v in chl.items()},
           "map_kb": kb, "map_kw": kw, "map_gamma": g, "map_demand": dm, "map_rough": rm, "map_dose": model.map_dose_,
           "post_weight_on_map_hydraulics": float(post.weights[band.members.index((float(dm), float(rm)))]),
           "wall_share_model_median": float(ls.wall_share.median()), "wall_share_truth_median": float(ts.wall_share.median()),
           "wall_share_mae": float(d.abs().mean()), "wall_share_bias": float(d.mean()),
           "wall_share_corr": float(np.corrcoef(ls.wall_share[both], ts.wall_share[both])[0, 1]) if both.sum() > 2 else float("nan"),
           "n_share_defined": int(both.sum()),
           "nonadditivity_model_min": float(ls.nonadditivity.min()), "nonadditivity_model_max": float(ls.nonadditivity.max()),
           "nonadditivity_truth_min": float(ts.nonadditivity.min()), "nonadditivity_truth_max": float(ts.nonadditivity.max()),
           "frac_wall_dominated_model": float((ls.wall_share[both] > 0.5).mean()),
           "frac_wall_dominated_truth": float((ts.wall_share[both] > 0.5).mean()),
           # start-up transient: the share of each junction's water on the scored day that is still the file's initial
           # contents (quality 0), not water that entered from a source during the 7-day run.  Exact for the model's
           # calibrated member (one dose at every source); for the truth, a count of junctions whose no-decay chlorine
           # is under 95% of the lowest possible source dose (so at least 5% initial water)
           "initial_water_share_median": float((1 - ls.source_water_share).median()),
           "initial_water_share_max": float((1 - ls.source_water_share).max()),
           "n_initial_water_over_5pct": int(((1 - ls.source_water_share) > 0.05).sum()),
           "n_unflushed_truth": int((sc.truth_loss_runs["no_decay"].mean() < 0.95 * 0.9 * sc.source_dose).sum())}
    if sc.structural:
        row["closed_pipe"] = sc.structural.get("closed_pipe") or ""
        row["tank"] = sc.structural.get("tank") or ""
    J = pd.DataFrame({"truth": truth, "seed": seed, "junction": uns, "nominal_age_h": nom.loc[uns].values,
                      "true_age_h": t_mean.loc[uns].values, "nominal_initial_share": nsh_mean.loc[uns].values,
                      "true_initial_share": tsh_mean.loc[uns].values, "true_min": sc.truth_daily_min.loc[uns].values,
                      "pred_median": pmin.loc[uns, "median"].values, "lo90": pmin.loc[uns, "lo90"].values,
                      "hi90": pmin.loc[uns, "hi90"].values, "flag": flag.loc[uns].values})
    fig = {"sc": sc, "X": X, "split": ls, "truth_split": ts, "row": row, "band": band}
    return row, J, fig


def _bins(nominal_mean_age: pd.Series, q: int) -> pd.Series:
    """1-based nominal-age quantile bin of every junction of the network (q = 3 terciles, 10 deciles)."""
    return pd.Series(pd.qcut(nominal_mean_age.rank(method="first"), q, labels=False) + 1, index=nominal_mean_age.index)


def error_by_age(J: pd.DataFrame, nominal_mean_age: pd.Series, q: int, lower_bound: pd.Series | None = None) -> pd.DataFrame:
    """Daily-minimum chlorine error of the unsampled junction-seeds by nominal-age bin, per seed and pooled over
    seeds ('all').  Bins are fixed by the network's nominal daily-mean age (every junction), ties broken by order.
    lower_bound (per junction, True where the nominal daily-mean age is a lower bound) gives each bin's
    frac_lower_bound: in the oldest bins the order of ages, not their size, is what the bins rest on."""
    b = _bins(nominal_mean_age, q)
    J = J.assign(bin=b.loc[J.junction].values)
    rows = []
    for truth in pd.unique(J.truth):                     # rows in run order: each seed, then all seeds pooled
        Jt = J[J.truth == truth]
        for seed_key, g0 in [(s, Jt[Jt.seed == s]) for s in pd.unique(Jt.seed)] + [("all", Jt)]:
            for k, g in g0.groupby("bin"):
                m = _chlorine_metrics(g.true_min, g.pred_median, g.lo90, g.hi90, g.flag.astype(bool))
                ages = nominal_mean_age[b == k]
                rows.append({"truth": truth, "seed": str(seed_key), "bin": int(k), "n_bins": q,
                             "nominal_age_lo_h": float(ages.min()), "nominal_age_hi_h": float(ages.max()),
                             "frac_lower_bound": float(lower_bound[b == k].mean()) if lower_bound is not None else float("nan"),
                             "nominal_age_mean_h": float(g.nominal_age_h.mean()), "true_age_mean_h": float(g.true_age_h.mean()), **m})
    return pd.DataFrame(rows)


def _pooled(rows: pd.DataFrame, J: pd.DataFrame) -> dict:
    """Pooled over seeds: age metrics weighted by the junctions each seed scored, chlorine metrics on every
    unsampled junction-seed, loss-split agreement as the mean of the per-seed numbers."""
    n = rows.n_junctions - rows.n_excluded
    w = n / n.sum()

    def rms(col):
        return float(np.sqrt((rows[col] ** 2 * w).sum()))

    def avg(col):
        return float((rows[col] * w).sum())

    def conv_avg(col, ncol="n_converged", square=False):
        n = rows[ncol]
        m = n > 0
        if not m.any():
            return float("nan")
        v = rows.loc[m, col] ** 2 if square else rows.loc[m, col]
        out = float((v * n[m]).sum() / n[m].sum())
        return float(np.sqrt(out)) if square else out
    chl = _chlorine_metrics(J.true_min, J.pred_median, J.lo90, J.hi90, J.flag.astype(bool))
    med_diff = (rows.wall_share_model_median - rows.wall_share_truth_median).abs()
    return {"n_seeds": int(len(rows)), "n_junction_seeds": int(rows.n_junctions.sum()), "n_excluded": int(rows.n_excluded.sum()),
            "rmse_h": rms("rmse_h"), "mae_h": avg("mae_h"), "bias_h": avg("bias_h"),
            "band_coverage": avg("band_coverage"), "band_coverage_daily_max": avg("band_coverage_daily_max"),
            "band_width_mean_h": avg("band_width_mean_h"),
            "band_coverage_by_seed": [round(float(v), 3) for v in rows.band_coverage],
            "posterior": {"rmse_h": rms("post_rmse_h"), "mae_h": avg("post_mae_h"), "bias_h": avg("post_bias_h"),
                          "band_coverage": avg("post_band_coverage"), "band_width_mean_h": avg("post_band_width_mean_h")},
            "rmse_h_by_seed": [round(float(v), 3) for v in rows.rmse_h],
            "post_rmse_h_by_seed": [round(float(v), 3) for v in rows.post_rmse_h],
            "n_seeds_posterior_worse": int((rows.post_rmse_h > rows.rmse_h).sum()),
            "lower_bounds": {"n_lower_bound_nominal": int(rows.n_lower_bound_nominal.iloc[0]),
                             "n_lower_bound_nominal_daily_max": int(rows.n_lower_bound_nominal_daily_max.iloc[0]),
                             "n_lower_bound_truth_mean": float(rows.n_lower_bound_truth.mean()),
                             "n_lower_bound_truth_daily_max_mean": float(rows.n_lower_bound_truth_daily_max.mean()),
                             "initial_share_nominal_median": float(rows.initial_share_nominal_median.iloc[0]),
                             "initial_share_truth_median_mean": float(rows.initial_share_truth_median.mean())},
            "converged_only": {"n_junction_seeds": int(rows.n_converged.sum()),
                               "n_junction_seeds_daily_max": int(rows.n_converged_daily_max.sum()),
                               "rmse_h": conv_avg("conv_rmse_h", square=True), "mae_h": conv_avg("conv_mae_h"),
                               "bias_h": conv_avg("conv_bias_h"), "band_coverage": conv_avg("conv_band_coverage"),
                               "band_coverage_daily_max": conv_avg("conv_band_coverage_daily_max", "n_converged_daily_max"),
                               "posterior": {"rmse_h": conv_avg("conv_post_rmse_h", square=True),
                                             "band_coverage": conv_avg("conv_post_band_coverage")}},
            "chlorine_daily_min_n8": chl,
            "loss_split": {"wall_share_model_median": float(rows.wall_share_model_median.mean()),
                           "wall_share_truth_median": float(rows.wall_share_truth_median.mean()),
                           "wall_share_mae": float(rows.wall_share_mae.mean()), "wall_share_bias": float(rows.wall_share_bias.mean()),
                           "median_abs_diff_mean": float(med_diff.mean()), "median_abs_diff_max": float(med_diff.max()),
                           "wall_share_corr_mean": float(rows.wall_share_corr.mean()),
                           "frac_wall_dominated_model": float(rows.frac_wall_dominated_model.mean()),
                           "frac_wall_dominated_truth": float(rows.frac_wall_dominated_truth.mean()),
                           "nonadditivity_model_range": [float(rows.nonadditivity_model_min.min()), float(rows.nonadditivity_model_max.max())],
                           "nonadditivity_truth_range": [float(rows.nonadditivity_truth_min.min()), float(rows.nonadditivity_truth_max.max())]},
            "start_up_transient": {"initial_water_share_median": float(rows.initial_water_share_median.mean()),
                                   "initial_water_share_max": float(rows.initial_water_share_max.max()),
                                   "n_initial_water_over_5pct_mean": float(rows.n_initial_water_over_5pct.mean()),
                                   "n_unflushed_truth_mean": float(rows.n_unflushed_truth.mean())}}


def plot_water_age(net: str, band, fig_in: dict, out: str, label: str, coverage_max: float) -> None:
    """Age map, the width of the 9-member range of the daily-max age (titled with that range's coverage of the true
    daily-max age), where chlorine is lost, and true daily-minimum chlorine against age."""
    sc, X, ls, row = fig_in["sc"], fig_in["X"], fig_in["split"], fig_in["row"]
    J = len(sc.junctions)
    size = max(12, int(4000 / J))
    amax = band.nominal.max()
    lo_x, hi_x = band.daily_range("max")
    j_old, h_old = oldest_water(band.nominal)
    n_lb = int(band.lower_bound("max").sum())
    oldest = (f"oldest: at least {h_old:.0f} h, junction {j_old}\n{n_lb} of {J} are lower bounds "
              f"(over {INITIAL_SHARE_MAX:.0%} starting water)" if band.lower_bound("max")[j_old]
              else f"oldest: {h_old:.0f} h, junction {j_old}")
    with plt.rc_context(BIG):
        fig, axes = plt.subplots(1, 4, figsize=(26, 6.6))
        wntr.graphics.plot_network(sc.wn, node_attribute=amax.to_dict(), node_size=size, node_cmap="YlOrBr",
                                   node_range=(0, float(amax.quantile(0.98))), ax=axes[0], link_width=0.6, add_colorbar=True,
                                   title=f"Water age, daily max, operator's model (h)\n7-day run; {oldest}")
        wntr.graphics.plot_network(sc.wn, node_attribute=(hi_x - lo_x).to_dict(), node_size=size, node_cmap="magma",
                                   node_range=(0, float((hi_x - lo_x).quantile(0.98))), ax=axes[1], link_width=0.6, add_colorbar=True,
                                   title=f"Width of the {label}, daily max (h)\n9 demand x roughness settings; holds the true\n"
                                         f"daily-max age at {coverage_max:.0%} of junctions (simulated)")
        ws = ls.wall_share.dropna()
        wntr.graphics.plot_network(sc.wn, node_attribute=ws.to_dict(), node_size=size, node_cmap="RdYlBu_r", node_range=(0, 1),
                                   ax=axes[2], link_width=0.6, add_colorbar=True,
                                   title=f"Where chlorine is lost: share at the pipe walls\ncalibrated member, {N_SAMPLES} samples, "
                                         f"scenario {sc.seed} (0 = all in the water)")
        ax = axes[3]
        sca = ax.scatter(sc.age_daily_mean_h.loc[sc.junctions], sc.truth_daily_min.loc[sc.junctions], c=X.loc[sc.junctions, "path_wall_index"],
                         cmap="viridis", s=28 if J < 200 else 8)
        plt.colorbar(sca, ax=ax, label="path_wall_index (old-pipe contact)")
        ax.axhline(THRESHOLD, color="tab:red", ls=":", lw=1.2)
        ax.text(ax.get_xlim()[0], THRESHOLD, f" {THRESHOLD} mg/L", color="tab:red", va="bottom", fontsize=10)
        ax.set(xlabel="nominal daily-mean water age (h)", ylabel="TRUE daily-minimum chlorine (mg/L)",
               title=f"Older water has less chlorine (scenario {sc.seed}, simulated truth)")
        ax.grid(alpha=0.3)
        fig.suptitle(f"{net} ({J} junctions): the operator's 7-day-run water age is off by {row['rmse_h']:.1f} h RMSE against the simulated "
                     f"truth (scenario {sc.seed}); walls take a median {ls.wall_share.median():.0%} of the chlorine loss "
                     f"(truth {fig_in['truth_split'].wall_share.median():.0%})", fontsize=15)
        fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def age_main(net: str, seeds, outdir: str = "outputs/chem", cache_dir: str = "outputs/cache", truths=tuple(AGE_TRUTHS)) -> dict:
    if not len(seeds):
        raise ValueError("age needs at least one seed")
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    rows, frames, fig_in, band = [], [], None, None
    for truth in truths:
        for s in seeds:
            ts = time.time()
            row, J, f = run_age_seed(net, s, truth, band, cache_dir)
            band = f["band"]
            rows.append(row); frames.append(J)
            if fig_in is None or (truth == "default" and s == seeds[0] and fig_in["row"]["truth"] != "default"):
                fig_in = f     # the figure shows the default truth's first seed (or the first run without it)
            print(f"{net} {truth} seed {s}: age RMSE {row['rmse_h']:.2f} h, range coverage {row['band_coverage']:.2f}, "
                  f"posterior {row['post_rmse_h']:.2f} h / {row['post_band_coverage']:.2f}; chlorine RMSE {row['rmse_min']:.3f}, "
                  f"recall {row['recall_min']:.2f}, {row['false_alarms_min']} false alarms; wall share {row['wall_share_model_median']:.2f} "
                  f"(truth {row['wall_share_truth_median']:.2f}) ({time.time() - ts:.1f} s)", flush=True)
    R = pd.DataFrame(rows)
    J = pd.concat(frames, ignore_index=True)
    nom = band.nominal.mean()
    R.to_csv(os.path.join(outdir, f"water_age_{net}.csv"), index=False)
    lb_mean, lb_max = band.lower_bound("mean"), band.lower_bound("max")
    dec = error_by_age(J, nom, 10, lb_mean)
    dec.to_csv(os.path.join(outdir, f"error_by_age_{net}.csv"), index=False)
    ter = error_by_age(J, nom, 3, lb_mean)
    lo_m, hi_m = band.daily_range("mean")
    lo_x0, hi_x0 = band.daily_range("max")
    pd.DataFrame({"junction": band.junctions, "nominal_age_daily_mean_h": nom.values,
                  "nominal_age_daily_max_h": band.nominal.max().values,
                  "range_daily_mean_lo_h": lo_m.values, "range_daily_mean_hi_h": hi_m.values,
                  "range_daily_max_lo_h": lo_x0.values, "range_daily_max_hi_h": hi_x0.values,
                  "initial_share_daily_mean": band.share("mean").values, "initial_share_daily_max": band.share("max").values,
                  "lower_bound_daily_mean": lb_mean.values, "lower_bound_daily_max": lb_max.values}
                 ).to_csv(os.path.join(outdir, f"water_age_junctions_{net}.csv"), index=False)
    by_truth = {t: _pooled(R[R.truth == t], J[J.truth == t]) for t in truths}
    covs = [by_truth[t][k] for t in truths for k in ("band_coverage", "band_coverage_daily_max")]
    label = band_label(min(covs))
    j_old, h_old = oldest_water(band.nominal)
    lo_x, hi_x = band.daily_range("max")
    d = by_truth[truths[0]]
    summary = {
        "generated_by": f"python -m residualmap.experiment_chem age {net} {len(seeds)}",
        "network": net, "seeds": list(map(int, seeds)), "truths": list(truths),
        "about": "Simulation only. Water age of the operator's nominal model and of 9 demand x roughness variants of it "
                 "(EPANET AGE), scored against the hidden truth's own water age; chlorine error by nominal-age group for "
                 "SimGP24 at 8 random daytime samples; the chlorine loss split between the water and the pipe walls for "
                 "the calibrated member and for the truth. Top-level numbers are for the first truth listed, pooled over "
                 "seeds; by_truth has both.",
        "definitions": {
            "age": "daily-mean water age (h) over the last day of the 7-day run, per junction",
            "rmse_h, mae_h, bias_h": "nominal (operator's file) minus true daily-mean age, over junctions not excluded",
            "band_coverage": "share of junctions whose true daily-mean age lies within [min, max] of the 9 members' daily-mean ages",
            "band_coverage_daily_max": "the same for the daily-maximum age (what the app's map shows)",
            "posterior": "ages weighted by the SimGP24 posterior over the 9 hydraulic members after 8 samples; band = its central 90%",
            "n_excluded": f"junction-seeds whose TRUE daily-max age exceeds {AGE_UNCONVERGED_H:g} h (the plan's rule; it "
                          f"misses most junctions whose age is a lower bound, see lower_bounds)",
            "lower bound": f"every age is a 7-day-run age; where more than {INITIAL_SHARE_MAX:g} of a junction's water is "
                           f"still the run's starting water (simulate.initial_water_share: daily mean, or any hour for the "
                           f"daily max) the age is a lower bound. lower_bounds counts them in the operator's model and "
                           f"the truth; converged_only repeats the age metrics on the junctions where neither is one",
            "chlorine_daily_min_n8": "SimGP24 at 8 random daytime samples, unsampled junctions, daily minimum; recall NaN when nothing to find",
            "wall_share": f"1 - L_bulk / L_tot per junction (mean over hours of ln C_ref / C); undefined below a loss of {MIN_LOSS}",
            "nonadditivity": "(L_bulk + L_wall) / L_tot; 1 on a single plug-flow path",
            "start_up_transient": "share of a junction's water on the scored day that is still the file's initial contents "
                                  "(1 - no-decay chlorine / dose) for the calibrated member: median over junctions averaged "
                                  "over seeds, and the largest; n_unflushed_truth counts truth junctions whose no-decay "
                                  "chlorine is under 95% of the lowest possible source dose"},
        "settings": {"n_samples": N_SAMPLES, "noise_sd_mgL": NOISE_SD, "threshold_mgL": THRESHOLD,
                     "sampling": "as the app's demo: default_rng(seed), 8 distinct junctions, hours in 07:00-17:00, readings rounded to 0.01 mg/L",
                     "truth_settings": truth_settings(net), "band_members": [list(m) for m in band.members],
                     "unconverged_age_h": AGE_UNCONVERGED_H, "initial_share_max": INITIAL_SHARE_MAX,
                     "band_coverage_bar": BAND_COVERAGE_BAR},
        "rmse_h": d["rmse_h"], "mae_h": d["mae_h"], "bias_h": d["bias_h"], "band_coverage": d["band_coverage"],
        "band_coverage_daily_max": d["band_coverage_daily_max"], "n_excluded": d["n_excluded"],
        "n_junction_seeds": d["n_junction_seeds"],
        "band_label": label,
        "band_label_rule": f"'90% band' only if every band_coverage and band_coverage_daily_max here is at least {BAND_COVERAGE_BAR}; "
                           f"otherwise '{RANGE_LABEL}'. Stricter than the plan, which asks only band_coverage (daily mean, "
                           f"default truth) >= {BAND_COVERAGE_BAR}: the app's map shows the daily-max range",
        "nominal_age": {"oldest_junction": j_old, "oldest_daily_max_h": h_old,
                        "oldest_range_h": [float(lo_x[j_old]), float(hi_x[j_old])],
                        "oldest_initial_share_max": float(band.share("max")[j_old]),
                        "oldest_initial_share_at_oldest_hour": float(band.initial_share[j_old].iloc[int(np.argmax(band.nominal[j_old].values))]),
                        "oldest_is_lower_bound": bool(lb_max[j_old]),
                        "n_lower_bound_daily_mean": int(lb_mean.sum()), "n_lower_bound_daily_max": int(lb_max.sum()),
                        "initial_share_daily_mean_median": float(band.share("mean").median()),
                        "daily_mean_median_h": float(nom.median()), "daily_mean_max_h": float(nom.max()),
                        "n_nominal_daily_max_above_unconverged": int((band.nominal.max() > AGE_UNCONVERGED_H).sum())},
        "by_truth": by_truth,
        "chlorine_by_age_tercile": {t: ter[(ter.truth == t) & (ter.seed == "all")].drop(columns=["truth", "seed"]).set_index("bin")
                                    .to_dict("index") for t in truths},
    }   # no timings in the file, so a rerun on an unchanged tree leaves it unchanged
    with open(os.path.join(outdir, f"water_age_{net}.json"), "w") as fh:
        json.dump(_clean(summary), fh, indent=1)
    plot_water_age(net, band, fig_in, os.path.join(outdir, f"water_age_{net}.png"), label, d["band_coverage_daily_max"])
    print(f"{net}: {time.time() - t0:.0f} s; band label: {label}", flush=True)
    return summary


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return None if np.isnan(x) else float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("age", help="task 9: water age, its range over the file's hydraulic errors, chlorine error by "
                                   "age, and where chlorine is lost")
    a.add_argument("net")
    a.add_argument("seeds", type=int, nargs="?", default=8)
    a.add_argument("--outdir", default="outputs/chem")
    a.add_argument("--cache", default="outputs/cache")
    a.add_argument("--truths", default=",".join(AGE_TRUTHS), help="comma list of: " + ", ".join(AGE_TRUTHS))
    args = ap.parse_args(argv)
    if args.cmd == "age":
        if args.seeds < 1:
            ap.error("seeds must be at least 1")
        truths = tuple(t for t in args.truths.split(",") if t)
        bad = [t for t in truths if t not in AGE_TRUTHS]
        if bad:
            ap.error(f"unknown truth(s) {bad}")
        age_main(args.net, tuple(range(args.seeds)), args.outdir, args.cache, truths)


if __name__ == "__main__":
    main()
