# Claude Code — introductory prompt for `ConradChallenge` (ResidualMap)

Paste everything below the line as your first message to Claude Code, from the repo root.

---

You are the coding agent on **ResidualMap**, a Conrad Challenge 2026–27 entry (Water Challenge). I plan and review; you write the code. Before writing anything, read in this order: `README.md`, `CHANGELOG.md`, `docs/feature_dictionary.md`, then every file in `residualmap/`, then `outputs/summary_Net3.json`. Then run `python -m residualmap.experiment Net3 2` once to confirm the environment works (about 90 seconds) and report the numbers you got against the ones in `README.md`.

## What the product is

An operator of a small water utility uploads the EPANET model they already have and types in the chlorine grab samples they already take. They get back: a chlorine map of every junction with a 90% band, the probability each junction is below the 0.2 mg/L minimum residual, the best next place to sample, and a one-page monthly report. Software only. The only hardware ever contemplated is a phone and, optionally, a Raspberry Pi reading a commodity chlorine probe as an extra data source. Never propose purpose-built hardware.

## Why it matters — the evidence base

- Most US community water systems are small; their operators do everything: O&M, management, sampling, customer service (AWWA small-systems page).
- Many small utilities lack the sensors, SCADA or digital meters an AI tool would need (Environmental Finance Center Network, Jan 2026). The product therefore runs on grab samples plus an EPANET file and nothing else.
- Failing to maintain a 0.2 mg/L chlorine residual is a standard small-system enforcement problem (EPA "special concerns for small water utilities").
- Compliance-paperwork tools for small systems exist already (1water.ai builds the Consumer Confidence Report from lab PDFs). ResidualMap must stay a physics-and-inference product, not a paperwork product — originality is 30% of the Conrad score.
- Judging: innovation 30%, storytelling 20%, practicality/proof-of-concept 20%, marketing 20%, finances 10%. Lean Canvas due Oct 29 2026, full submission Jan 7 2027. Every task should end in something a judge can see in a 2–3 minute video.

## What exists and what it scores (iteration 2, verified)

- `simulate.py` — hidden truth on WNTR-bundled networks: per-pipe wall decay tied to Hazen-Williams roughness plus lognormal noise, ±15% global demand, per-node demand noise, ±10% source dose, ±10% roughness mismatch. Nominal model gives age, hydraulics, pipe table, and a nominal-chlorine simulator.
- `features.py` — 24 physics features; `CORE` = age, 3-D hydraulic-distance embedding, path wall-exposure index, distance to source.
- `simgp.py` — the main model. Grid of 75 EPANET runs over (kb, kw, gamma), Bayesian weights from samples, discrepancy GP, parameter spread in the error bars.
- `surrogate.py` — iteration-1 decay-law GP, baselines, acquisition rules (random / uncertainty / straddle).
- `pinn.py` — graph-PINN baseline. Tested, not adopted: at 3–15 samples it loses to the calibrated simulator on every metric that matters.
- Net3, 8 seeds, scored on unsampled junctions: calibrated-simulator GP RMSE 0.14 mg/L at 3 samples, 0.09 at 15; recall of junctions below 0.2 mg/L 0.80 → 0.87 (0.96 with the straddle rule); 90% coverage 0.75–0.82. "Mean of samples" — today's practice — flags zero low nodes in every scenario.
- Finding for the pitch: 11 of 92 junctions below 0.2 mg/L at 14:00, 31 at 22:00, 40% at midnight. Operators sample by day; the network is worst at night.

## Iteration-3 tasks, in order

Do them in order. After each: run the full Net3 experiment, append a dated entry to `CHANGELOG.md` with the new numbers, update the limitations section of `README.md`. Stop and report if a task makes RMSE or recall worse than the numbers above by more than 10%.

### Task 1 — Time-aware model (the pitch)
Predict the full 24-h profile and the daily minimum, from daytime samples only.
- Extend `simulator_grid` to keep all 24 hours (the sims already produce them; only the sampling hour is stored now).
- Samples become (junction, hour, mg/L). Calibration weights use the simulated value at each sample's own hour. The discrepancy GP gets `hour` as sin/cos inputs plus age at that hour.
- New evaluation targets: **daily minimum** per junction and **P(daily min < 0.2)**.
- Acquisition chooses (junction, hour) pairs restricted to 07:00–17:00.
- Acceptance: with 15 daytime samples, recall on junctions whose daily minimum is below 0.2 mg/L ≥ 0.80; figure `outputs/day_vs_night_predicted_Net3.png` showing the true 22:00 map next to the prediction from daytime samples.

### Task 2 — Close the calibration gap
Coverage is 0.75–0.82 for a 90% band; target 0.88–0.95.
- Add hydraulic mismatch to the grid: a demand multiplier axis (0.85, 1.0, 1.15) and a roughness multiplier axis (0.9, 1.0, 1.1). Cache grows to 675 runs (~2 min on Net3); keep the cache.
- Report a reliability diagram `outputs/reliability_Net3.png` (nominal vs empirical coverage at 50/80/90/95%).
- Acceptance: mean coverage across seeds and n in [0.88, 0.95] with RMSE within 10% of the current numbers.

### Task 3 — Structural mismatch stress test
Real EPANET files are wrong in ways the grid cannot represent.
- Add a `structural_noise` option to `build_scenario`: with probability 0.5 close one random non-bridge pipe (check connectivity with networkx first), and multiply one random tank's initial level by 0.7.
- Run the full comparison with it on. Report how much the calibrated-simulator GP degrades and whether the discrepancy GP absorbs it. If recall falls below 0.7 at n=15, propose a fix in `CHANGELOG.md`; do not hide it.

### Task 4 — Route optimiser
Operators pick a fixed monthly route; give them a better one.
- Batch acquisition: choose K (junction, hour) pairs at once — greedy on straddle score with a diversity penalty in hydraulic-distance space.
- Compare against random K sites and against the K highest-demand sites (a plausible operator heuristic) for K = 5, 8, 12, over 8 seeds.
- Acceptance: optimised route wins on recall of daily-minimum violations at every K; figure `outputs/route_comparison_Net3.png`.

### Task 5 — Runs on small and large real networks
Run everything on `Net2` (35 junctions, tank-fed — set a tank source if there is no reservoir) and `ky4` (959 junctions, bundled with WNTR). Produce the same figures. If the grid is too slow on ky4, sub-sample the grid or fit a quick emulator; say what you did.

### Task 6 — Operator-facing demo
A single Streamlit app (`app.py`): upload `.inp`, enter samples as rows, see the four panels (map, uncertainty, P(below 0.2), "sample here next" with a one-sentence reason each), download a one-page PDF. No accounts, no database. Runs on a laptop.

### Task 7 — Real-data path
Write `docs/pilot_protocol.md`: what a pilot utility gives us (their `.inp`, six months of grab logs with timestamps and sample-tap locations), how sample taps map to junctions, what they get back monthly, and the exact validation — hold out the latest month, predict it, report RMSE and coverage. Target pilot: Purissima Hills Water District (Los Altos Hills, CA) or any district that will share an EPANET file.

## Constraints

- Keep every existing entry point runnable; extend, don't rewrite.
- Ask before adding any dependency other than `gpytorch`/`botorch`/`streamlit`. If you move the GP to GPyTorch/BoTorch, match the author's existing stack (BoTorch `ModelListGP`, `qNEHVI`).
- No claim in a figure title, docstring or README that is not backed by a number in `outputs/`.
- Every figure readable in a judged video: large fonts, one message per panel, the finding in the title.
- Do not build a PINN unless a task above asks for it. The verdict in `README.md` stands until real multi-month data exists.

## Report-back format

After each task: task number, files changed, the metric table in the same columns as `README.md`, coverage, and one sentence on what a judge can now see that they could not before.
