# ResidualMap

**One sentence:** a small water utility takes a handful of chlorine grab samples a month; ResidualMap turns them into a chlorine map of the whole pipe network with calibrated uncertainty, flags every junction likely to be below the 0.2 mg/L minimum, and tells the operator where to sample next.

Conrad Challenge 2026–27, Water Challenge. Software only. Working name — rename freely.

## The problem in one paragraph

Most of the ~50,000 US community water systems are small and run by one or two operators with no sensors and no SCADA. They keep the required minimum chlorine residual by grabbing a few samples and assuming the rest of the network looks the same. Every experiment in this repo shows that assumption flags **zero** low-chlorine nodes — because the average of a few samples sits above 0.2 mg/L while 10–25% of the network is under it at midday, and up to 40% at night.

## What the model does (iteration 2)

The operator's own EPANET model is the prior. EPANET already solves the transport physics exactly — mixing at junctions, tank turnover, travel time through every pipe with its real length, diameter and roughness. What it lacks is the decay chemistry. ResidualMap calibrates that from grab samples and models what's left:

1. **Simulator grid.** Chlorine is simulated on the operator's model over a grid of the three unknowns — bulk decay `kb`, wall decay `kw`, and `gamma`, how much old (rough) pipe accelerates wall decay. 75 EPANET runs, ~13 s, cached per network.
2. **Bayesian calibration.** Each grid member is weighted by how well it explains the samples. The weighted mean of log-chlorine is the mean function; the weighted spread is the parameter uncertainty that feeds the map's error bars.
3. **Discrepancy GP.** A Matérn-ARD Gaussian process on water age, a hydraulic-distance embedding, the path wall-exposure index and distance to source learns what the simulator still gets wrong (demand mismatch, pipe-level heterogeneity, dose drift).
4. **Outputs per junction:** median, 90% band, `P(residual < 0.2 mg/L)`.
5. **Next sample:** three rules — random, max-uncertainty, and *straddle* (level-set estimation: sample where the model is least sure whether a node is above or below 0.2).

This is a Kennedy–O'Hagan calibration-plus-discrepancy model with a grid posterior instead of MCMC.

Three other surrogates are kept for comparison: the iteration-1 decay-law GP (first-order decay on water age as the mean), the same GP on the full 24-feature "deep dive" set (`docs/feature_dictionary.md`), and a graph physics-informed neural network (`residualmap/pinn.py`).

## Results — Net3 (EPANET example based on North Marin Water District, CA; 92 junctions), 14:00 snapshot, 8 scenario weeks

Hidden truth per scenario: per-pipe wall decay tied to roughness plus lognormal noise, ±15% global demand, per-node demand noise, ±10% source dose, ±10% pipe roughness (hydraulic model mismatch). Grab samples carry ±0.03 mg/L noise. **All scores are on junctions that were never sampled.**

**Q1/Q2 — which physics helps at 3–15 random samples?**

| model | n=3 RMSE / recall | n=8 RMSE / recall | n=15 RMSE / recall | 90% coverage (n=8) |
|---|---|---|---|---|
| mean of samples (today's practice) | 0.31 / 0.00 | 0.28 / 0.00 | 0.28 / 0.00 | – |
| nearest sampled node | 0.27 / 0.06 | 0.23 / 0.35 | 0.18 / 0.27 | – |
| decay law only (no GP) | 0.27 / 0.66 | 0.25 / 0.41 | 0.25 / 0.29 | – |
| decay-law GP, core features (iter 1) | 0.27 / 0.66 | 0.22 / 0.67 | 0.17 / 0.59 | 0.43 |
| decay-law GP, rich features | 0.27 / 0.66 | 0.23 / 0.59 | 0.20 / 0.60 | 0.51 |
| graph-PINN, rich features | 0.24 / 0.05 | 0.19 / 0.35 | 0.14 / 0.45 | 0.38 |
| **calibrated-simulator GP, core** | **0.14 / 0.80** | **0.12 / 0.82** | **0.09 / 0.87** | 0.75 |
| calibrated-simulator GP, rich | 0.13 / 0.82 | 0.11 / 0.81 | 0.10 / 0.79 | 0.82 |

RMSE in mg/L; recall = share of junctions truly below 0.2 mg/L that the model flags (P > 0.5).

**Q3 — does smart sampling beat random? (calibrated-simulator GP)**

| sampling rule | n=8 RMSE / recall / F1 | n=12 RMSE / recall / F1 | n=15 RMSE / recall / F1 |
|---|---|---|---|
| random | 0.12 / 0.82 / 0.73 | 0.10 / 0.86 / 0.82 | 0.09 / 0.87 / 0.84 |
| uncertainty | 0.13 / 0.82 / 0.78 | 0.12 / 0.89 / 0.81 | 0.11 / 0.90 / 0.80 |
| **straddle** | 0.11 / 0.88 / 0.83 | 0.11 / **0.96** / 0.90 | 0.10 / **0.96** / **0.91** |

**What the numbers say**

- **Three grab samples plus the calibrated simulator beat every other model at fifteen.** RMSE 0.14 mg/L at n=3 vs 0.17 for the iteration-1 GP at n=15.
- **Eighty percent of the low-chlorine junctions are found from three samples;** 96% with twelve samples chosen by the straddle rule.
- **The calibration recovers the physics.** Averaged over seeds, the posterior mode lands at `kw` ≈ 0.7 m/day and `gamma` ≈ 0.9 — the truth was 0.7 and 1.0. From eight grab samples the model works out that old pipes are eating the chlorine.
- **More features did not help; more physics did.** The 24-feature "deep dive" set slightly hurts the decay-law GP at n ≥ 12 (over-fitting) and is a wash for the simulator GP. Physics pays off when it enters through the simulator, not the feature list. The features remain useful for explaining *why* a node is low.
- **PINN verdict.** The graph-PINN (first-order decay enforced along every steadily-directed pipe, trained on all nodes' features with the sample labels) beats the decay-law GP on RMSE at n ≥ 10 but is far worse at finding violations (recall 0.05–0.45), is badly over-confident (coverage 0.38), and is 10× slower. At 3–15 samples a neural network has nothing to learn from that EPANET does not already compute. A PINN would earn its place in iteration 3+ when months of samples exist and the target is the pipe-level decay field.
- **Uncertainty is now usable but still under-covers:** 0.75–0.82 for a 90% band (was 0.43). Remaining gap is hydraulic model mismatch not represented in the grid — see tasks.

**The finding that drives the pitch — `outputs/day_vs_night_Net3.png`:** in scenario 0, 11 of 92 junctions are below 0.2 mg/L at 14:00; 31 are at 22:00, 40% at midnight. Operators sample in the day. The network is worst at night.

Figures: `outputs/map_Net3_n{3,8,15}.png`, `outputs/curves_Net3.png`, `outputs/strategies_Net3.png`, `outputs/day_vs_night_Net3.png`. Numbers: `outputs/results_Net3.csv`, `outputs/summary_Net3.json`, `outputs/features_Net3_seed0.csv`.

## Honest limitations

- **Truth shares the nominal model's topology.** Real EPANET files have closed valves that are open, missing pipes, wrong tank levels. Roughness mismatch is modelled (±10%); structural mismatch is not yet.
- **Snapshot, not daily minimum.** The model predicts residual at the sampling hour. The compliance number is the daily minimum, which happens at night.
- **Synthetic truth.** No real grab-sample data yet. `docs/pilot_protocol.md` (to be written) is the path to it.
- **Grid posterior is coarse** (5×5×3). Fine for three parameters; MCMC or an emulator if more are added.
- **Net3 is mid-size.** ky4/ky10 (≈950 junctions, bundled with WNTR) run through the same code but the simulator grid takes minutes per network.

## Run

```bash
pip install -r requirements.txt
python -m residualmap.experiment Net3 8      # network, number of scenario seeds; ~5 min; writes outputs/
python -m residualmap.experiment ky4 2       # bigger real network, fewer seeds
```

## Layout

```
residualmap/simulate.py    truth scenario (hidden decay, demand, dose, roughness noise) + nominal model hydraulics, age, pipe table, nominal-chlorine simulator
residualmap/features.py    24 physics features from the .inp and one hydraulic run; CORE subset; FEATURE_DOCS
residualmap/simgp.py       calibrated-simulator GP (main model): grid posterior over kb/kw/gamma + discrepancy GP
residualmap/surrogate.py   decay-law GP (iteration 1), baselines, acquisition rules
residualmap/pinn.py        graph-PINN baseline (numpy, L-BFGS, deep ensemble)
residualmap/experiment.py  sequential-sampling loop, metrics, all figures
docs/feature_dictionary.md what every feature means physically
CHANGELOG.md               dated results per iteration
CLAUDE_CODE_PROMPT.md      the prompt to start iteration 3 in Claude Code
```
