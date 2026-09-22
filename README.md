# ResidualMap

**One sentence:** a small water utility takes a handful of chlorine grab samples a month; ResidualMap turns them into a chlorine map of the whole pipe network with calibrated uncertainty, flags every junction likely to be below the 0.2 mg/L minimum, and tells the operator where to sample next.

Conrad Challenge 2026–27, Water Challenge. Software only. Working name — rename freely.

## The problem in one paragraph

Most of the ~50,000 US community water systems are small and run by one or two operators with no sensors and no SCADA. They keep the required minimum chlorine residual by grabbing a few samples and assuming the rest of the network looks the same. Every experiment in this repo shows that assumption flags **zero** low-chlorine nodes — because the average of a few samples sits above 0.2 mg/L while 10–25% of the network is under it at midday, and up to 40% at night.

## What the model does (iteration 3)

The operator's own EPANET model is the prior. EPANET already solves the transport physics exactly — mixing at junctions, tank turnover, travel time through every pipe with its real length, diameter and roughness. What it lacks is the decay chemistry and the honest size of its own hydraulic errors. ResidualMap calibrates the first from grab samples and models the second:

1. **Simulator grid.** Chlorine is simulated on the operator's model over a grid of the three decay unknowns — bulk decay `kb`, wall decay `kw`, and `gamma`, how much old (rough) pipe accelerates wall decay — times two hydraulic-mismatch axes, demand ×{0.85, 1, 1.15} and pipe roughness ×{0.9, 1, 1.1}. 675 EPANET runs, ~20 s on Net3, cached per network; all 24 hours of every run are kept. A source-dose axis ×{0.9 … 1.1} costs nothing: first-order decay is linear in concentration, so it is an exact offset in log space.
2. **Bayesian calibration.** Each grid member is weighted by how well it explains the samples (Student-t likelihood, scale 0.35 in ln C — the model error, not the grab-sample noise — so one junction where the operator's model is structurally wrong cannot hijack the fit). The weighted mean of log-chlorine is the mean function; the weighted spread is the parameter uncertainty that feeds the map's error bars.
3. **Local hydraulic error.** Calibrating the global demand and roughness multipliers does not remove the operator's per-node and per-pipe errors of the same size; the grid's own spread along the hydraulic axes at each (hour, junction) is carried into the band as their variance.
4. **Discrepancy GP.** A Matérn-ARD Gaussian process on water age, a hydraulic-distance embedding, the path wall-exposure index, distance to source and — in the time-aware model — the hour of day learns what the simulator still gets wrong.
5. **Outputs per junction:** median, 90% band, `P(residual < 0.2 mg/L)` at any hour; the daily minimum and `P(daily minimum < 0.2)` by Monte Carlo over grid members, local hydraulic deviations and joint 24-h GP draws.
6. **Next sample:** random, max-uncertainty, *straddle* (level-set estimation: sample where the model is least sure whether a node is above or below 0.2) and *straddle on the daily minimum* (pick the junction whose daily-minimum verdict is least sure, sample it at the daytime hour that constrains it most). Rules use the reducible part of the uncertainty (parameters + GP), not the local hydraulic term that no sample can shrink.

This is a Kennedy–O'Hagan calibration-plus-discrepancy model with a grid posterior instead of MCMC.

Three other surrogates are kept for comparison: the iteration-1 decay-law GP (first-order decay on water age as the mean), the same GP on the full 24-feature "deep dive" set (`docs/feature_dictionary.md`), and a graph physics-informed neural network (`residualmap/pinn.py`).

## Results — Net3 (EPANET example based on North Marin Water District, CA; 92 junctions), 14:00 snapshot, 8 scenario weeks

Hidden truth per scenario: per-pipe wall decay tied to roughness plus lognormal noise, ±15% global demand, per-node demand noise, ±10% source dose, ±10% pipe roughness (hydraulic model mismatch). Grab samples carry ±0.03 mg/L noise. **All scores are on junctions that were never sampled.**

**Q1/Q2 — which physics helps at 3–15 random samples?**

| model | n=3 RMSE / recall | n=8 RMSE / recall | n=15 RMSE / recall | 90% coverage (n=8) |
|---|---|---|---|---|
| mean of samples (today's practice) | 0.31 / 0.00 | 0.28 / 0.00 | 0.28 / 0.00 | – |
| nearest sampled node | 0.27 / 0.07 | 0.23 / 0.35 | 0.18 / 0.27 | – |
| decay law only (no GP) | 0.27 / 0.66 | 0.25 / 0.41 | 0.25 / 0.29 | – |
| decay-law GP, core features (iter 1) | 0.27 / 0.66 | 0.22 / 0.67 | 0.17 / 0.59 | 0.43 |
| decay-law GP, rich features | 0.27 / 0.66 | 0.23 / 0.59 | 0.20 / 0.60 | 0.51 |
| graph-PINN, rich features | 0.24 / 0.05 | 0.19 / 0.33 | 0.14 / 0.38 | 0.37 |
| calibrated-simulator GP, core, iteration-2 grid (75 runs) | 0.14 / 0.80 | 0.12 / 0.82 | 0.09 / 0.87 | 0.75 |
| **calibrated-simulator GP, core, iteration 3** | **0.12 / 0.74** | **0.11 / 0.87** | **0.09 / 0.84** | **0.94** |
| calibrated-simulator GP, rich, iteration 3 | 0.12 / 0.75 | 0.11 / 0.76 | 0.09 / 0.73 | 0.92 |

RMSE in mg/L; recall = share of junctions truly below 0.2 mg/L that the model flags (P > 0.5). Precision of the iteration-3 flags: 0.84 / 0.82 / 0.87 (iteration-2 grid: 0.68 / 0.72 / 0.83).

**Q3 — does smart sampling beat random? (calibrated-simulator GP, iteration 3)**

| sampling rule | n=8 RMSE / recall / F1 | n=12 RMSE / recall / F1 | n=15 RMSE / recall / F1 | found by sample or flag, n=8 / 12 / 15 |
|---|---|---|---|---|
| random | 0.11 / 0.87 / 0.81 | 0.09 / 0.83 / 0.81 | 0.08 / 0.84 / 0.83 | 0.88 / 0.83 / 0.83 |
| uncertainty | 0.10 / 0.74 / 0.73 | 0.10 / 0.72 / 0.64 | 0.10 / 0.70 / 0.65 | 0.77 / 0.81 / 0.82 |
| **straddle** | 0.09 / **0.93** / **0.94** | 0.10 / **0.95** / **0.93** | 0.09 / **0.92** / **0.93** | 0.91 / 0.88 / 0.84 |

The last column counts a violation as found when its own grab sample read below 0.2 mg/L *or* the model flags it — the operational number; the unsampled-only recall penalises a rule for sampling the violators.

**Calibration — `outputs/reliability_Net3.png`.** Nominal vs empirical coverage, mean over seeds and n = 3…15, random samples: iteration-2 grid 0.45 / 0.67 / 0.76 / 0.81 at the 50 / 80 / 90 / 95% levels; iteration 3 **0.68 / 0.90 / 0.95 / 0.97**. The 90% band now holds the truth 0.95 of the time (target 0.88–0.95).

**What the numbers say**

- **Three grab samples plus the calibrated simulator beat every other model at fifteen.** RMSE 0.12 mg/L at n=3 vs 0.17 for the iteration-1 GP at n=15.
- **Nine in ten of the low-chlorine junctions are found with eight samples chosen by the straddle rule** (recall 0.93, F1 0.94); random samples find 0.74 at three and 0.84 at fifteen.
- **The band is honest now.** Iteration 2 under-covered (0.75 for a 90% band); with demand, roughness and dose on the grid and the local hydraulic error carried into the band, coverage is on the diagonal at 80–95% and slightly wide at 50%.
- **The calibration recovers the physics.** Averaged over seeds, the posterior mode lands at `kw` ≈ 0.75 m/day and `gamma` ≈ 0.94 — the truth was 0.7 and 1.0. From eight grab samples the model works out that old pipes are eating the chlorine.
- **More features did not help; more physics did.** The 24-feature "deep dive" set slightly hurts the decay-law GP at n ≥ 12 (over-fitting) and is a wash or worse for the simulator GP. Physics pays off when it enters through the simulator, not the feature list. The features remain useful for explaining *why* a node is low.
- **PINN verdict.** The graph-PINN (first-order decay enforced along every steadily-directed pipe, trained on all nodes' features with the sample labels) beats the decay-law GP on RMSE at n ≥ 10 but is far worse at finding violations (recall 0.05–0.38), is badly over-confident (coverage 0.37), and is 10× slower. At 3–15 samples a neural network has nothing to learn from that EPANET does not already compute. A PINN would earn its place when months of samples exist and the target is the pipe-level decay field.
- **The max-uncertainty rule is no longer useful** (recall 0.70 at n=15, was 0.90): it spends its samples on the violators themselves and then flags only 70% of the harder ones left; use the straddle rule.

**The finding that drives the pitch — `outputs/day_vs_night_Net3.png`:** in scenario 0, 11 of 92 junctions are below 0.2 mg/L at 14:00; 31 are at 22:00, 40% at midnight. Operators sample in the day. The network is worst at night.

Figures: `outputs/map_Net3_n{3,8,15}.png`, `outputs/curves_Net3.png`, `outputs/strategies_Net3.png`, `outputs/day_vs_night_Net3.png`. Numbers: `outputs/results_Net3.csv`, `outputs/summary_Net3.json`, `outputs/features_Net3_seed0.csv`.

## Iteration 3 — daytime samples predict the night

The compliance number is the **daily minimum**, which happens at night; operators sample by day. `residualmap/simgp.py::SimGP24` takes samples as (junction, hour, mg/L), scores every grid member on the simulated value at each sample's own hour, gives the discrepancy GP the hour (sin/cos) and the water age at that hour — with a length-scale floor of about three hours, because a discrepancy driven by the diurnal demand pattern cannot change in fifteen minutes — and predicts the whole 24-h profile at every junction. The daily minimum and `P(daily min < 0.2)` come from Monte-Carlo draws (grid member and dose by weight, a random local hydraulic deviation, a joint 24-h GP draw per junction; take the minimum). Sampling hours are restricted to 07:00–17:00.

**Net3, 8 seeds, samples only between 07:00 and 17:00, target = daily minimum, scored on unsampled junctions**

| model / sampling rule | n=3 RMSE / recall | n=8 RMSE / recall | n=15 RMSE / recall | 90% coverage of daily min (n=8) |
|---|---|---|---|---|
| mean of samples (today's practice) | 0.32 / 0.00 | 0.32 / 0.00 | 0.34 / 0.00 | – |
| time-blind calibrated-simulator GP (every sample treated as 14:00) | 0.22 / 0.27 | 0.22 / 0.30 | 0.21 / 0.35 | 0.43 |
| time-aware GP, task-1 state (decay grid only), random | 0.08 / 0.86 | 0.13 / 0.80 | 0.12 / 0.89 | 0.65 |
| time-aware GP, random daytime samples | 0.07 / 0.94 | 0.07 / 0.95 | 0.08 / 0.98 | 0.94 |
| time-aware GP, hourly straddle | 0.07 / 0.94 | 0.08 / 0.90 | 0.15 / 0.97 | 0.87 |
| time-aware GP, max-uncertainty | 0.07 / 0.94 | 0.10 / 0.94 | 0.14 / 0.98 | 0.84 |
| **time-aware GP, straddle on daily minimum** | **0.07 / 0.94** | **0.07 / 0.99** | **0.12 / 0.99** | **0.92** |

RMSE in mg/L on the daily minimum; recall = share of junctions whose true daily minimum is below 0.2 mg/L that the model flags (`P(daily min < 0.2) > 0.5`). Precision of the flags at n=15: 0.88 straddle-on-daily-minimum (F1 0.91), 0.93 random (F1 0.95). Recall at the 22:00 snapshot from the same daytime samples, straddle-on-daily-minimum: 0.74 → 0.91 → 0.95. Coverage of the daily minimum, mean over n at the 50 / 80 / 90 / 95% levels: 0.56 / 0.85 / 0.93 / 0.96 random, 0.53 / 0.80 / 0.88 / 0.93 straddle-on-daily-minimum (task-1 state: 0.34 / 0.55 / 0.65 / 0.72).

- **Eight daytime grab samples find 99% of the junctions that go below 0.2 mg/L at night** (per seed at n=15: 0.97, 1.0, 1.0, 0.96, 1.0, 1.0, 1.0, 1.0). The same samples fed to the time-blind model find 30%; the mean of samples finds none.
- **The calibration recovers the physics from daytime data:** posterior mode at n=8 averages kb 0.31 /day, kw 0.75 m/day, γ 0.94 against a truth of 0.40 / 0.70 / 1.0.
- **The band on the daily minimum is honest on average** (0.93 for a 90% band with random samples) but still narrows with n: 0.96 at n=3, 0.87 at n=15 (0.79 with the straddle-on-daily-minimum rule, which deliberately samples the borderline junctions).

Figures: `outputs/day_vs_night_predicted_Net3.png` (true 22:00 map next to the prediction from 15 daytime samples; true daily-minimum violations next to `P(daily min < 0.2)`), `outputs/curves_time_Net3.png`, `outputs/reliability_Net3.png`. Numbers: `outputs/results_time_Net3.csv`, `outputs/summary_Net3.json["time_aware_daily_min"]`.

## Honest limitations

- **Truth shares the nominal model's topology.** Real EPANET files have closed valves that are open, missing pipes, wrong tank levels. Global demand and roughness mismatch are on the grid and local (per-node, per-pipe) mismatch is in the band; structural mismatch is not yet.
- **The band is honest on average, not at every n.** Snapshot: 0.95 for a 90% band over n = 3–15, but the 50% band over-covers (0.68) — slightly too wide in the middle. Daily minimum: 0.93 on average, narrowing from 0.96 at n=3 to 0.87 at n=15 (0.79 with the straddle-on-daily-minimum rule). Treat `P(daily min < 0.2)` from a straddle-chosen set of 15 samples as slightly over-confident.
- **Recall traded for precision on the 14:00 snapshot.** The better-calibrated posterior shrinks the low junctions' medians toward the grid centre: recall under random samples is about 5% lower than with the iteration-2 grid (precision 10–16 points higher, F1 higher). The max-uncertainty rule lost 20% recall and should not be used.
- **Synthetic truth.** No real grab-sample data yet. `docs/pilot_protocol.md` (to be written) is the path to it.
- **Grid posterior is coarse** (5×5×3 decay × 3×3 hydraulic × 5 dose). Fine for six parameters with a 20-second cache; MCMC or an emulator if more are added.
- **Net3 is mid-size.** ky4/ky10 (≈950 junctions, bundled with WNTR) run through the same code but the 675-run simulator grid takes minutes per network and ~190 MB in memory at float32.

## Run

```bash
pip install -r requirements.txt
python -m residualmap.experiment Net3 8      # network, number of scenario seeds; ~3 min; writes outputs/
python -m residualmap.experiment ky4 2       # bigger real network, fewer seeds
```

## Layout

```
residualmap/simulate.py    truth scenario (hidden decay, demand, dose, roughness noise) + nominal model hydraulics, age, pipe table, nominal-chlorine simulator
residualmap/features.py    24 physics features from the .inp and one hydraulic run; CORE subset; FEATURE_DOCS
residualmap/simgp.py       calibrated-simulator GP (main model): grid posterior over kb/kw/gamma + discrepancy GP;
                           SimGP24 = the time-aware version (samples at any hour, 24-h profile, daily minimum)
residualmap/surrogate.py   decay-law GP (iteration 1), baselines, acquisition rules (per junction, and per junction-hour)
residualmap/pinn.py        graph-PINN baseline (numpy, L-BFGS, deep ensemble)
residualmap/experiment.py  sequential-sampling loops (snapshot and time-aware), metrics, all figures
docs/feature_dictionary.md what every feature means physically
CHANGELOG.md               dated results per iteration
CLAUDE_CODE_PROMPT.md      the prompt to start iteration 3 in Claude Code
```
