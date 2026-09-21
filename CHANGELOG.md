# Changelog

## 2026-09-21 — iteration 3, task 1: time-aware model

- `simgp.simulator_grid_24h`: the simulator grid now keeps all 24 hours of every EPANET run (cache `outputs/cache/grid24_<net>.pkl`; `simulator_grid` slices its hour from it, so iteration-2 numbers are unchanged).
- `simgp.SimGP24`: samples are (junction, hour, mg/L). Grid members are scored on the simulated value at each sample's own hour; the discrepancy GP takes age-at-that-hour, the hydraulic embedding, wall index, distance to source and sin/cos(hour). The daily minimum and P(daily min < 0.2) are Monte-Carlo: draw a grid member by weight (whole-day profile) + a joint 24-h GP draw per junction, take the minimum (256 draws).
- `simgp.grid_weights`: calibration likelihood is now Student-t (ν=3, scale 0.25) for the time-aware model. Diagnosis on Net3 seed 0: the straddle rule samples junctions near 0.2 mg/L, which are the fronts where the operator's model is structurally wrong (junctions 20, 131, 151: truth 3× every grid member); under a Gaussian four such readings pulled the calibration to kb=0.1 and daily-min recall fell from 0.88 (n=4) to 0.62 (n=15). Under Student-t, on the straddle loop over 4 seeds at n=15: RMSE 0.11 vs 0.19, precision 0.85 vs 0.64, coverage 0.75 vs 0.29. `SimGP` keeps `lik="gauss"` by default so iteration-2 results stand; the option is shared.
- `surrogate.acquire_time`: (junction, hour) acquisition restricted to 07:00–17:00: random, max-uncertainty, hourly straddle, and `straddle_min` (straddle on the daily minimum picks the junction; the daytime hour with the largest predictive sd picks the hour).
- `experiment.run_scenario_time` + figures `outputs/curves_time_Net3.png`, `outputs/day_vs_night_predicted_Net3.png`; numbers in `outputs/results_time_Net3.csv` and `summary_Net3.json["time_aware_daily_min"]`. Baselines: mean of samples; the iteration-2 model treating every sample as a 14:00 sample ("time-blind"), scored against the daily minimum.
- Net3, 8 seeds, daytime samples only, scored on unsampled junctions, target = daily minimum: `straddle_min` recall 0.86 (n=3) → 0.94 (n=8) → 0.93 (n=15), precision 0.94 at n=15, F1 0.93; RMSE of the daily minimum 0.08 → 0.12 → 0.13 mg/L; recall at the 22:00 snapshot 0.78 → 0.95 → 0.85. Random daytime samples: recall 0.86 → 0.80 → 0.89. Time-blind model: 0.37 → 0.27 → 0.38. Mean of samples: 0. Per-seed recall at n=15 (`straddle_min`): 1.0, 0.75, 0.88, 0.82, 1.0, 1.0, 1.0, 1.0. **Acceptance (recall ≥ 0.80 at 15 daytime samples): met.**
- Calibration recovers the physics from daytime samples: posterior mode at n=15 averages kb 0.31 /day, kw 0.75 m/day, γ 0.94 (truth 0.40, 0.70, 1.0).
- Iteration-2 snapshot numbers re-verified in this environment: identical to 4 decimals (Python 3.13, WNTR 1.5.0, scikit-learn 1.9.1). PINN and rich-feature GP rows drift (L-BFGS / 23-dim kernel optimiser), nothing else.
- Not met yet, carried to task 2: 90% band on the daily minimum covers only 0.81 (n=3) → 0.65 (n=15), and the daily-minimum RMSE grows with n (0.08 → 0.13). With more daytime samples the posterior concentrates on one grid member (effective members ≈ 2–6); the best member's night profile is biased (truth − sim ≈ −0.4 in ln C at 00:00–03:00 for the true-parameter member: hydraulic mismatch and per-pipe decay noise are not on the grid) and the discrepancy GP, trained by day, reverts to zero at night. The demand/roughness axes of task 2 are the intended fix.
- Fixed: `experiment.__main__` printed `cov` instead of `cov90` and crashed after writing all outputs.

## 2026-09-21 — iteration 2

- Added `simgp.py`: calibrated-simulator GP. Operator's EPANET model run over a 5×5×3 grid of (bulk decay, wall decay, old-pipe sensitivity); Bayesian weights from grab samples; discrepancy GP on residuals; parameter spread feeds uncertainty.
- Added `features.py`: 24 physics features from the `.inp` + one hydraulic run (age, path materials, diameters, wall-exposure integral, tank distance, velocities, sloshing fraction, pressure, elevation, demand, degree). Data dictionary in `docs/feature_dictionary.md`.
- Added `pinn.py`: graph physics-informed NN baseline to test whether a PINN helps at 3–15 samples.
- Truth now includes ±10% pipe-roughness mismatch between the operator's model and reality.
- Full last-day profile retained; day-vs-night figure added.
- Net3, 8 seeds, random samples, unsampled junctions: calibrated-simulator GP RMSE 0.14 (n=3) → 0.09 (n=15) mg/L, recall 0.80 → 0.87, coverage 0.75–0.82. Straddle rule: recall 0.96, F1 0.91 at n=15. Iteration-1 GP for reference: RMSE 0.27 → 0.17, recall 0.66 → 0.59. PINN: RMSE 0.24 → 0.14, recall 0.05 → 0.45, coverage 0.38.
- Verdict: physics helps through the simulator, not the feature list; PINN not justified at this sample count.

## 2026-09-21 — iteration 1

- WNTR/EPANET truth on Net3 with hidden per-pipe wall decay; decay-law GP in log space on water age + hydraulic-distance embedding; random / uncertainty / straddle acquisition; baselines.
- Net3, 8 seeds: GP+uncertainty RMSE 0.136 mg/L at n=15 vs mean-of-samples 0.270; GP+straddle recall 0.89 on nodes < 0.2 mg/L vs 0.00 for mean-of-samples. 90% coverage only 0.32–0.66.
- Finding: 10% of junctions below 0.2 mg/L at 08:00, 37% at 22:00.
