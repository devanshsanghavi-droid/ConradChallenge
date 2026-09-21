# Changelog

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
