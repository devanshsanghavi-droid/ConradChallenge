"""
pilot.py — the real-data path (docs/pilot_protocol.md made executable).

A pilot utility gives us their EPANET .inp and their grab-sample log.  The validation is a hold-out in
time: fit on the earlier months, predict the held-out month's samples at their own junction and hour,
report RMSE, MAE, coverage of the 50/80/90/95% bands, and the violation call, against two baselines an
operator could run by hand (last reading at the same tap; mean of all recent samples).

    python -m residualmap.pilot --inp their_model.inp --samples grab_log.csv --taps tap_map.csv --dose 1.2
    python -m residualmap.pilot --synthetic Net3          # six synthetic months, to see the report format
    python -m residualmap.pilot --inp ... --plant plant_log.csv   # with the plant's monthly water temperature (task 10)
    python -m residualmap.pilot --inp ... --samples total_log.csv --taps ... --disinfectant chloramine --dose 2.0 \
                                --ph 8.0 --cl2n 4.5
                                       # a chloraminated system (task 12): total chlorine, the chloramine grid,
                                       # threshold 0.5 mg/L total chlorine unless --threshold is given; with the
                                       # plant's pH and Cl2:N the bulk rate gets the chloramine mode's prior (model
                                       # (c), the one held to task 12's bars), without them a uniform prior (model (b))

grab_log.csv  : date (YYYY-MM-DD), time (HH:MM), tap_id, free_chlorine_mgL   [optional: method, notes]
                or, for a chloraminated system (--disinfectant chloramine), total_chlorine_mgL instead (never both
                columns: a log that mixes the two species is refused; docs/example_grab_log_total.csv)
tap_map.csv   : tap_id, junction_id                                          [optional: flush_min, notes]
plant_log.csv : month (YYYY-MM), temp_C                                      [optional: toc_mgL, dose_mgL]
                With a plant log the model is seasonal.SeasonalSimGP24: every grab sample is explained at its own
                month's water temperature and the held-out month is predicted at its own.  dose_mgL, when given, is
                the month's plant dose; first-order decay makes it an exact offset against --dose.  toc_mgL is carried
                and not used: task 11 tested a TOC-aware model in simulation (residualmap/organics.py) and it stopped
                at its pre-registered stop rule, so it is not wired in (journal, task 11).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

from .chemistry import Chemistry, Warming
from .chloramine import CA_DOSE_MGL
from .features import build_features
from .simgp import DAY_HOURS, SimGP24
from .simulate import build_scenario, load, nominal_scenario

LEVELS = (50, 80, 90, 95)


SPECIES_COLUMN = {"free_chlorine": "free_chlorine_mgL", "chloramine": "total_chlorine_mgL"}


def load_log(samples_csv: str, taps_csv: str, disinfectant: str = "free_chlorine") -> pd.DataFrame:
    """Join the grab log to the tap map -> junction, hour, y, month, date.  The log's reading column names its species:
    free_chlorine_mgL (free chlorine) or total_chlorine_mgL (total chlorine, a chloraminated system).  A log with both
    columns is refused (one likelihood must never mix species), and so is a log whose species does not match
    `disinfectant`."""
    if disinfectant not in SPECIES_COLUMN:
        raise ValueError(f"disinfectant must be one of {sorted(SPECIES_COLUMN)}")
    log = pd.read_csv(samples_csv, dtype={"tap_id": str})
    have = [c for c in SPECIES_COLUMN.values() if c in log.columns]
    if len(have) > 1:
        what = ("the free-chlorine model reads free chlorine only: drop or rename the total_chlorine_mgL column"
                if disinfectant == "free_chlorine" else
                "the chloramine mode reads total chlorine only: drop or rename the free_chlorine_mgL column")
        what += " (or, if the log mixes a free-chlorine and a chloraminated system, split it by disinfectant)"
        raise ValueError(f"the grab log has both {have[0]} and {have[1]}: free and total chlorine are different "
                         f"measurements and one model never mixes them; {what}")
    want = SPECIES_COLUMN[disinfectant]
    if have != [want]:
        found = have[0] if have else "neither column"
        raise ValueError(f"--disinfectant {disinfectant} needs a {want} column; the grab log has {found}")
    taps = pd.read_csv(taps_csv, dtype={"tap_id": str, "junction_id": str})
    df = log.merge(taps[["tap_id", "junction_id"]], on="tap_id", how="left")
    missing = df.junction_id.isna()
    if missing.any():
        print(f"warning: {int(missing.sum())} samples at taps with no junction in the tap map are dropped: "
              f"{sorted(df.loc[missing, 'tap_id'].unique())}")
        df = df[~missing]
    ts = pd.to_datetime(df.date.astype(str) + " " + df.time.astype(str))
    return pd.DataFrame({"junction": df.junction_id.astype(str), "hour": ts.dt.hour.astype(int),
                         "y": df[want].astype(float), "month": ts.dt.to_period("M").astype(str),
                         "date": ts.dt.date.astype(str), "tap_id": df.tap_id.astype(str)}).reset_index(drop=True)


def synthetic_log(network: str = "Net3", months: int = 6, per_month: int = 10, rotating: int = 3, seed: int = 0,
                  schedule: list[dict] | None = None, return_truth: bool = False, **truth_kw):
    """Six 'months' of daytime grab samples: a fixed route of `per_month` taps plus `rotating` validation
    taps at different junctions each month.  The network's pipe-level truth is fixed; each month re-draws
    the operating truth (bulk decay ±20%, demand ±15% global and ±15% per node, dose ±10%).

    schedule=None is the committed log, draw for draw: month_seed = 100 + 10 seed + m (that formula repeats across
    seeds after ten months: seed 0's month 10 is seed 1's month 0).  A schedule (seasonal.plant_schedule: one dict per
    month, {'temp_C': the plant water temperature, 'soil_temp_C': None or the soil temperature for in-network warming})
    gives the seasonal path (task 10): month_seed = 1000 + 100 seed + m, distinct for any seed and up to 100 months;
    month m's truth is build_scenario(chem=Chemistry(temp_C=...), warming=Warming(soil_temp_C=...) when a soil
    temperature is given); every reading carries its month's logged plant temperature (temp_C); and the plant log
    (month, temp_C, toc_mgL, dose_mgL; TOC not logged here, the dose is the set point) comes back as a third item.
    Task 11: a schedule month may also give 'toc_mgL' (the plant TOC: the truth's Chemistry gets it, every reading of
    the month carries it, and the plant log logs it), 'kinetics' and 'phi' (the truth's kinetics, e.g. 'clark' with its
    demand per mg TOC; default first order) and 'dose_mgL' (the month's plant dose: the truth runs at that dose, and the
    readings and the plant log carry it).  A schedule without these keys gives task 10's log, draw for draw.
    return_truth=True appends a list with each month's truth: month (1-based), month_seed, truth_by_hour,
    truth_daily_min and chem (the truth's chemistry and hidden draws; None on the committed path)."""
    rng = np.random.default_rng(seed)
    if schedule is None:
        sc0 = build_scenario(network, seed=seed, **truth_kw)
        junctions = sc0.junctions
    else:
        if len(schedule) < months:
            raise ValueError(f"the schedule has {len(schedule)} months, {months} asked for")
        junctions = load(network).junction_name_list          # the same list build_scenario returns
    taps = list(rng.choice(junctions, per_month, replace=False))           # the utility's fixed route
    rows, truths, plant = [], [], []
    for m in range(months):
        label = f"{2026 + m // 12}-{m % 12 + 1:02d}"
        if schedule is None:
            ms = 100 + seed * 10 + m
            sc = build_scenario(network, seed=seed, month_seed=ms, **truth_kw)
        else:
            cond = schedule[m]
            ms = 1000 + 100 * seed + m
            warm = Warming(soil_temp_C=cond["soil_temp_C"]) if cond.get("soil_temp_C") is not None else None
            chem = Chemistry(temp_C=cond["temp_C"], toc_mgL=cond.get("toc_mgL"), kinetics=cond.get("kinetics", "first"),
                             phi=cond.get("phi"))
            kw_m = dict(truth_kw, source_dose=float(cond["dose_mgL"])) if cond.get("dose_mgL") is not None else truth_kw
            sc = build_scenario(network, seed=seed, month_seed=ms, chem=chem, warming=warm, **kw_m)
            plant.append({"month": label, "temp_C": float(cond["temp_C"]),
                          "toc_mgL": float(cond["toc_mgL"]) if cond.get("toc_mgL") is not None else float("nan"),
                          "dose_mgL": float(kw_m.get("source_dose", 1.2))})
        others = [j for j in sc.junctions if j not in taps]
        rot = list(rng.choice(others, rotating, replace=False)) if rotating else []
        for k, j in enumerate(taps + rot):
            h = int(rng.choice(DAY_HOURS)); day = int(rng.integers(1, 28))
            y = float(np.clip(sc.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None))
            row = {"junction": j, "hour": h, "y": round(y, 2), "month": label, "date": f"{label}-{day:02d}",
                   "tap_id": f"T{k + 1}" if k < per_month else f"V{m + 1}-{k - per_month + 1}"}
            if schedule is not None:
                row["temp_C"] = float(schedule[m]["temp_C"])
                for key in ("toc_mgL", "dose_mgL"):      # task 11; absent from task 10's schedules
                    if schedule[m].get(key) is not None:
                        row[key] = float(schedule[m][key])
            rows.append(row)
        if return_truth:
            truths.append({"month": m + 1, "month_seed": ms, "truth_by_hour": sc.truth_by_hour,
                           "truth_daily_min": sc.truth_daily_min, "chem": sc.chem})
    out = (pd.DataFrame(rows), network) + ((pd.DataFrame(plant),) if schedule is not None else ())
    return out + ((truths,) if return_truth else ())


def load_plant(path: str) -> pd.DataFrame:
    """The plant log: month (YYYY-MM), temp_C (required), toc_mgL and dose_mgL (optional), one row per month."""
    p = pd.read_csv(path, dtype={"month": str})
    if not {"month", "temp_C"} <= set(p.columns):
        raise ValueError("the plant log needs the columns month (YYYY-MM) and temp_C")
    if p.month.duplicated().any():
        raise ValueError(f"months listed twice in the plant log: {sorted(p.month[p.month.duplicated()])}")
    return p


DEFAULT_DOSE_MGL = {"free_chlorine": 1.2, "chloramine": CA_DOSE_MGL}   # mg/L; typical chloramine doses are 1.5 to 4.0


def validate(inp: str, log: pd.DataFrame, dose: float | None = None, threshold: float | None = None,
             window_months: int = 3, holdout_months: int | None = None, cache_dir: str = "outputs/cache", seed: int = 0,
             plant: pd.DataFrame | None = None, bank_cache: str = "readwrite",
             disinfectant: str = "free_chlorine", ph: float | None = None,
             cl2n: float | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rolling hold-out: for each held-out month, fit on the `window_months` before it and predict its
    samples.  Returns (per-sample predictions, per-month summary).
    plant: the plant log (load_plant).  None is the committed temperature-blind SimGP24, unchanged.  With a plant log,
    seasonal.SeasonalSimGP24 explains each sample at its month's water temperature (and dose, when dose_mgL is given
    for every month of the plant log; otherwise the dose is taken as --dose throughout) and predicts the held-out month
    at its own; the summary gains the month's temperature, kb20 and the posterior over
    the temperature hypotheses.  The bank of grids at the log's temperatures is built once and cached (bank_cache, as
    seasonal.covariate_bank's cache: 'read' never writes a cache file).
    disinfectant 'chloramine' (task 12): the log is total chlorine and the model is the chloramine mode's (the
    chloramine grid, its dose axis and likelihood scale).  With the plant's logged pH and Cl2:N (ph, cl2n) the bulk rate
    gets the chloramine mode's prior, model (c), the one held to task 12's bars; without them the prior is uniform,
    model (b).  A plant log is refused (no seasonal chloramine model is built).  dose and threshold default to the
    disinfectant's (1.2 mg/L free chlorine and 0.2; 2.0 mg/L total chlorine and 0.5)."""
    if disinfectant not in DEFAULT_DOSE_MGL:
        raise ValueError(f"disinfectant must be one of {sorted(DEFAULT_DOSE_MGL)}")
    if disinfectant == "chloramine" and plant is not None:
        raise NotImplementedError("a plant log (temperature) with chloramine: no seasonal chloramine model is built")
    if (ph is None) != (cl2n is None):
        raise ValueError("give both the plant's pH and its Cl2:N ratio, or neither")
    if ph is not None and disinfectant != "chloramine":
        raise ValueError("pH and Cl2:N set the chloramine mode's prior; they are not used with free chlorine")
    dose = DEFAULT_DOSE_MGL[disinfectant] if dose is None else dose
    threshold = (0.5 if disinfectant == "chloramine" else 0.2) if threshold is None else threshold
    sc = nominal_scenario(inp, sample_hour=14, source_dose=dose)
    X = build_features(sc)
    unknown = sorted(set(log.junction) - set(sc.junctions))
    if unknown:
        raise ValueError(f"junction IDs not in the .inp: {unknown[:10]}")
    months = sorted(log.month.unique())
    held = months[-holdout_months:] if holdout_months else months[window_months:]
    bank = None
    if plant is not None:
        from .seasonal import SeasonalSimGP24, covariate_bank
        pl = plant.set_index(plant.month.astype(str))
        missing = sorted(set(months) - set(pl.index))
        if missing:
            raise ValueError(f"months in the grab log with no plant-log row: {missing}")
        ratio = (pl.dose_mgL.astype(float) / dose) if "dose_mgL" in pl and pl.dose_mgL.notna().all() else pd.Series(1.0, index=pl.index)
        log = log.assign(temp_C=log.month.map(pl.temp_C.astype(float)), dose_ratio=log.month.map(ratio))
        bank = covariate_bank(sc, sorted(set(log.temp_C)), cache_dir=cache_dir, cache=bank_cache)
    preds, rows = [], []
    for m in held:
        train_months = [t for t in months if t < m][-window_months:]
        train, test = log[log.month.isin(train_months)], log[log.month == m]
        if bank is None and disinfectant == "chloramine":
            from .chloramine import ca_condition, prior_log_vector
            from .simgp import DOSES_CA, LIK_SD_CA
            model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir, grid="chloramine", cond=ca_condition(),
                            lik_sd=LIK_SD_CA, doses=DOSES_CA, threshold=threshold)
            if ph is not None:
                model.log_prior = prior_log_vector(model.params, float(ph), float(cl2n))
            model.fit(train[["junction", "hour", "y"]])
            if model.map_dose_ in (DOSES_CA[0], DOSES_CA[-1]):
                print(f"warning: held-out month {m}: the fitted effective dose sits at the edge of the chloramine dose "
                      f"axis (x{model.map_dose_:.2f} of {dose:g} mg/L): check --dose, the total chlorine leaving the "
                      f"plant", file=sys.stderr, flush=True)
        elif bank is None:
            model = SimGP24(sc, X, seed=seed, cache_dir=cache_dir).fit(train[["junction", "hour", "y"]])
        else:
            model = SeasonalSimGP24(sc, X, bank, seed=seed, cache_dir=cache_dir).fit(
                train[["junction", "hour", "y", "temp_C", "dose_ratio"]], target_temp_C=float(pl.temp_C[m]),
                target_dose_ratio=float(ratio[m]))
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
            if disinfectant == "chloramine":
                from .simgp import DOSES_CA
                rows[-1].update({"model": "(c) pH and Cl2:N prior" if ph is not None else "(b) uniform prior",
                                 "map_dose_at_axis_edge": bool(model.map_dose_ in (DOSES_CA[0], DOSES_CA[-1]))})
            if bank is not None:
                rows[-1].update({"temp_C": float(pl.temp_C[m]), "kb20": model.kb20(), "map_hypothesis": model.map_hypothesis_,
                                 **{f"P_{k}": v for k, v in model.hypothesis_posterior().items()}})
    return pd.concat(preds, ignore_index=True), pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inp", help="the utility's EPANET model")
    ap.add_argument("--samples", help="grab log CSV: date,time,tap_id,free_chlorine_mgL (or total_chlorine_mgL)")
    ap.add_argument("--taps", help="tap map CSV: tap_id,junction_id")
    ap.add_argument("--dose", type=float, default=None,
                    help="free chlorine leaving the plant, mg/L (default 1.2); with chloramine, total chlorine (default 2.0)")
    ap.add_argument("--threshold", type=float, default=None,
                    help="minimum residual, mg/L: default 0.2 free chlorine, or 0.5 total chlorine with chloramine (a "
                         "common utility operating target, not a California rule: California requires a detectable residual)")
    ap.add_argument("--disinfectant", choices=sorted(SPECIES_COLUMN), default="free_chlorine",
                    help="free_chlorine (the default) or chloramine (total chlorine; task 12)")
    ap.add_argument("--ph", type=float, default=None, help="chloramine only: the plant water's logged pH (with --cl2n "
                                                           "it sets the bulk-rate prior, model (c))")
    ap.add_argument("--cl2n", type=float, default=None, help="chloramine only: the plant's chlorine to ammonia-N ratio "
                                                             "by mass (with --ph)")
    ap.add_argument("--window", type=int, default=3, help="months of history to fit on")
    ap.add_argument("--synthetic", metavar="NETWORK", help="demonstrate on a synthetic six-month log from this bundled network")
    ap.add_argument("--plant", help="plant log CSV: month,temp_C[,toc_mgL,dose_mgL]; makes the model temperature-aware")
    ap.add_argument("--out", default="outputs/pilot")
    a = ap.parse_args()
    if a.threshold is None:
        a.threshold = 0.5 if a.disinfectant == "chloramine" else 0.2
    if a.dose is None:
        a.dose = DEFAULT_DOSE_MGL[a.disinfectant]
        if a.disinfectant == "chloramine":
            print(f"--dose not given: taking {a.dose:g} mg/L total chlorine leaving the plant (typical chloramine doses "
                  f"are 1.5 to 4.0); give the plant's own", file=sys.stderr)
    if (a.ph is not None or a.cl2n is not None) and a.disinfectant != "chloramine":
        ap.error("--ph and --cl2n set the chloramine mode's prior; use them with --disinfectant chloramine")
    if (a.ph is None) != (a.cl2n is None):
        ap.error("give both --ph and --cl2n, or neither")
    if a.synthetic and a.disinfectant != "free_chlorine":
        ap.error("--synthetic makes a free-chlorine log; the chloramine mode's simulated test is "
                 "python -m residualmap.experiment <net> --disinfectant=chloramine")
    os.makedirs(a.out, exist_ok=True)
    if a.synthetic:
        truth_kw = {"Net2": dict(kb_per_day=0.10, kw_m_per_day=0.20), "ky4": dict(source_dose=2.0)}.get(a.synthetic, {})
        log, inp = synthetic_log(a.synthetic, **truth_kw)
        dose = truth_kw.get("source_dose", a.dose)
        log.to_csv(os.path.join(a.out, f"synthetic_log_{a.synthetic}.csv"), index=False)
    else:
        if not (a.inp and a.samples and a.taps):
            ap.error("--inp, --samples and --taps are required (or --synthetic NETWORK)")
        log, inp, dose = load_log(a.samples, a.taps, a.disinfectant), a.inp, a.dose
    plant = load_plant(a.plant) if a.plant else None
    preds, summary = validate(inp, log, dose=dose, threshold=a.threshold, window_months=a.window, plant=plant,
                              disinfectant=a.disinfectant, ph=a.ph, cl2n=a.cl2n)
    tag = (a.synthetic or os.path.splitext(os.path.basename(inp))[0]) + ("_plant" if plant is not None else "")
    preds.to_csv(os.path.join(a.out, f"predictions_{tag}.csv"), index=False)
    summary.to_csv(os.path.join(a.out, f"validation_{tag}.csv"), index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_columns", 40)
    print(f"\n== hold-out validation, {tag}: fit on the previous {a.window} months, predict the held-out month's samples at their own junction and hour")
    print(summary[["held_out_month", "taps", "n_train", "n_test", "rmse", "rmse_persistence", "rmse_network_mean", "coverage50", "coverage80", "coverage90", "coverage95",
                   "band90_width_mgL", "n_true_below", "recall_below", "false_alarms", "map_kb", "map_kw"]
                  + (["temp_C", "kb20", "P_H0"] if plant is not None else [])].round(3).to_string(index=False))
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
