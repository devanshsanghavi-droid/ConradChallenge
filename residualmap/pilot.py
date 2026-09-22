"""
pilot.py — the real-data path (docs/pilot_protocol.md made executable).

A pilot utility gives us their EPANET .inp and their grab-sample log.  The validation is a hold-out in
time: fit on the earlier months, predict the held-out month's samples at their own junction and hour,
report RMSE, MAE, coverage of the 50/80/90/95% bands, and the violation call, against two baselines an
operator could run by hand (last reading at the same tap; mean of all recent samples).

    python -m residualmap.pilot --inp their_model.inp --samples grab_log.csv --taps tap_map.csv --dose 1.2
    python -m residualmap.pilot --synthetic Net3          # six synthetic months, to see the report format

grab_log.csv : date (YYYY-MM-DD), time (HH:MM), tap_id, free_chlorine_mgL   [optional: method, notes]
tap_map.csv  : tap_id, junction_id                                          [optional: flush_min, notes]
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

from .features import build_features
from .simgp import DAY_HOURS, SimGP24
from .simulate import build_scenario, nominal_scenario

LEVELS = (50, 80, 90, 95)


def load_log(samples_csv: str, taps_csv: str) -> pd.DataFrame:
    """Join the grab log to the tap map -> junction, hour, y, month, date."""
    log = pd.read_csv(samples_csv, dtype={"tap_id": str})
    taps = pd.read_csv(taps_csv, dtype={"tap_id": str, "junction_id": str})
    df = log.merge(taps[["tap_id", "junction_id"]], on="tap_id", how="left")
    missing = df.junction_id.isna()
    if missing.any():
        print(f"warning: {int(missing.sum())} samples at taps with no junction in the tap map are dropped: "
              f"{sorted(df.loc[missing, 'tap_id'].unique())}")
        df = df[~missing]
    ts = pd.to_datetime(df.date.astype(str) + " " + df.time.astype(str))
    return pd.DataFrame({"junction": df.junction_id.astype(str), "hour": ts.dt.hour.astype(int),
                         "y": df.free_chlorine_mgL.astype(float), "month": ts.dt.to_period("M").astype(str),
                         "date": ts.dt.date.astype(str), "tap_id": df.tap_id.astype(str)}).reset_index(drop=True)


def synthetic_log(network: str = "Net3", months: int = 6, per_month: int = 10, rotating: int = 3, seed: int = 0,
                  **truth_kw) -> tuple[pd.DataFrame, str]:
    """Six 'months' of daytime grab samples: a fixed route of `per_month` taps plus `rotating` validation
    taps at different junctions each month.  The network's pipe-level truth is fixed; each month re-draws
    the operating truth (bulk decay ±20%, demand ±15% global and ±15% per node, dose ±10%)."""
    rng = np.random.default_rng(seed)
    sc0 = build_scenario(network, seed=seed, **truth_kw)
    taps = list(rng.choice(sc0.junctions, per_month, replace=False))           # the utility's fixed route
    rows = []
    for m in range(months):
        sc = build_scenario(network, seed=seed, month_seed=100 + seed * 10 + m, **truth_kw)
        others = [j for j in sc.junctions if j not in taps]
        rot = list(rng.choice(others, rotating, replace=False)) if rotating else []
        for k, j in enumerate(taps + rot):
            h = int(rng.choice(DAY_HOURS)); day = int(rng.integers(1, 28))
            y = float(np.clip(sc.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None))
            rows.append({"junction": j, "hour": h, "y": round(y, 2), "month": f"2026-{m + 1:02d}", "date": f"2026-{m + 1:02d}-{day:02d}",
                         "tap_id": f"T{k + 1}" if k < per_month else f"V{m + 1}-{k - per_month + 1}"})
    return pd.DataFrame(rows), network


def validate(inp: str, log: pd.DataFrame, dose: float = 1.2, threshold: float = 0.2, window_months: int = 3,
             holdout_months: int | None = None, cache_dir: str = "outputs/cache", seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rolling hold-out: for each held-out month, fit on the `window_months` before it and predict its
    samples.  Returns (per-sample predictions, per-month summary)."""
    sc = nominal_scenario(inp, sample_hour=14, source_dose=dose)
    X = build_features(sc)
    unknown = sorted(set(log.junction) - set(sc.junctions))
    if unknown:
        raise ValueError(f"junction IDs not in the .inp: {unknown[:10]}")
    months = sorted(log.month.unique())
    held = months[-holdout_months:] if holdout_months else months[window_months:]
    preds, rows = [], []
    for m in held:
        train_months = [t for t in months if t < m][-window_months:]
        train, test = log[log.month.isin(train_months)], log[log.month == m]
        model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir).fit(train[["junction", "hour", "y"]])
        z_mu, z_sd = model.predict_hours()
        jidx = np.array([sc.junctions.index(j) for j in test.junction]); h = test.hour.values
        mu, sd = z_mu[h, jidx], z_sd[h, jidx]
        p = test.assign(pred=np.exp(mu), z_mu=mu, z_sd=sd, p_below=norm.cdf((np.log(threshold) - mu) / sd),
                        train_months=",".join(train_months))
        for q in LEVELS:
            k = norm.ppf(0.5 + q / 200)
            p[f"in{q}"] = (test.y.values >= np.exp(mu - k * sd)) & (test.y.values <= np.exp(mu + k * sd))
        # baselines an operator could run by hand
        last = train.sort_values("date").groupby("tap_id").y.last()
        p["persistence"] = test.tap_id.map(last).fillna(train.y.mean()).values
        p["network_mean"] = train.y.mean()
        # a tap the model has seen in training tests the month-to-month forecast; a tap it has never seen
        # tests the MAP — the claim that matters
        p["seen_tap"] = test.junction.isin(set(train.junction)).values
        preds.append(p)
        for subset, mask in (("all", np.ones(len(p), bool)), ("seen_taps", p.seen_tap.values), ("new_taps", ~p.seen_tap.values)):
            if mask.sum() == 0:
                continue
            t, q_ = test.y.values[mask], p[mask]
            viol = t < threshold; flagged = q_.p_below.values > 0.5
            rows.append({"held_out_month": m, "taps": subset, "trained_on": ",".join(train_months), "n_train": len(train), "n_test": int(mask.sum()),
                         "rmse": float(np.sqrt(np.mean((t - q_.pred.values) ** 2))), "mae": float(np.mean(np.abs(t - q_.pred.values))),
                         "rmse_persistence": float(np.sqrt(np.mean((t - q_.persistence.values) ** 2))),
                         "rmse_network_mean": float(np.sqrt(np.mean((t - q_.network_mean.values) ** 2))),
                         **{f"coverage{q}": float(q_[f"in{q}"].mean()) for q in LEVELS},
                         "band90_width_mgL": float(np.median(np.exp(q_.z_mu + 1.645 * q_.z_sd) - np.exp(q_.z_mu - 1.645 * q_.z_sd))),
                         "n_true_below": int(viol.sum()), "recall_below": float((viol & flagged).sum() / viol.sum()) if viol.sum() else float("nan"),
                         "false_alarms": int((~viol & flagged).sum()),
                         "map_kb": model.map_params_[0], "map_kw": model.map_params_[1], "map_gamma": model.map_params_[2], "map_dose": model.map_dose_})
    return pd.concat(preds, ignore_index=True), pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inp", help="the utility's EPANET model")
    ap.add_argument("--samples", help="grab log CSV: date,time,tap_id,free_chlorine_mgL")
    ap.add_argument("--taps", help="tap map CSV: tap_id,junction_id")
    ap.add_argument("--dose", type=float, default=1.2, help="free chlorine leaving the plant, mg/L")
    ap.add_argument("--threshold", type=float, default=0.2)
    ap.add_argument("--window", type=int, default=3, help="months of history to fit on")
    ap.add_argument("--synthetic", metavar="NETWORK", help="demonstrate on a synthetic six-month log from this bundled network")
    ap.add_argument("--out", default="outputs/pilot")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.synthetic:
        truth_kw = {"Net2": dict(kb_per_day=0.10, kw_m_per_day=0.20), "ky4": dict(source_dose=2.0)}.get(a.synthetic, {})
        log, inp = synthetic_log(a.synthetic, **truth_kw)
        dose = truth_kw.get("source_dose", a.dose)
        log.to_csv(os.path.join(a.out, f"synthetic_log_{a.synthetic}.csv"), index=False)
    else:
        if not (a.inp and a.samples and a.taps):
            ap.error("--inp, --samples and --taps are required (or --synthetic NETWORK)")
        log, inp, dose = load_log(a.samples, a.taps), a.inp, a.dose
    preds, summary = validate(inp, log, dose=dose, threshold=a.threshold, window_months=a.window)
    tag = a.synthetic or os.path.splitext(os.path.basename(inp))[0]
    preds.to_csv(os.path.join(a.out, f"predictions_{tag}.csv"), index=False)
    summary.to_csv(os.path.join(a.out, f"validation_{tag}.csv"), index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 40)
    print(f"\n== hold-out validation, {tag}: fit on the previous {a.window} months, predict the held-out month's samples at their own junction and hour")
    print(summary[["held_out_month", "taps", "n_train", "n_test", "rmse", "rmse_persistence", "rmse_network_mean", "coverage50", "coverage80", "coverage90", "coverage95",
                   "band90_width_mgL", "n_true_below", "recall_below", "false_alarms", "map_kb", "map_kw"]].round(3).to_string(index=False))
    agg = {}
    for subset, g in summary.groupby("taps"):
        w = g.n_test / g.n_test.sum()
        agg[subset] = {"n": int(g.n_test.sum()), "rmse": float((g.rmse * w).sum()), "rmse_persistence": float((g.rmse_persistence * w).sum()),
                       "rmse_network_mean": float((g.rmse_network_mean * w).sum()), **{f"coverage{q}": float((g[f"coverage{q}"] * w).sum()) for q in LEVELS}}
        print(f"mean over held-out months, {subset}:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in agg[subset].items()})
    with open(os.path.join(a.out, f"validation_{tag}.json"), "w") as f:
        json.dump({"per_month": summary.to_dict("records"), "mean": agg}, f, indent=2, default=float)


if __name__ == "__main__":
    main()
