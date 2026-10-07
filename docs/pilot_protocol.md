# Pilot protocol — the real-data path

Everything in this repository is scored against a simulated truth. This document is how that changes: what a pilot utility gives us, how their sample taps become model junctions, what they get back every month, and the exact validation we will run and report — before we claim anything to anyone.

**Target pilot:** a local district that will share an EPANET file. The ask is small: one file they already have, one spreadsheet they already keep, and twenty minutes on the phone.

The validation described here runs as `python -m residualmap.pilot --inp <their .inp> --samples <grab log> --taps <tap map>`; `python -m residualmap.pilot --synthetic Net3` shows the report format on a synthetic six-month log (`docs/example_grab_log.csv`, `docs/example_tap_map.csv` are in that format). With the plant's monthly water temperature, add `--plant <plant log>` (`docs/example_plant_log.csv` is in that format; the temperatures in it are an assumed coastal schedule, not measurements).

## 1. What the utility gives us

| item | what exactly | why we need it |
|---|---|---|
| **EPANET model** | The `.inp` they run (or their consultant runs), with demand patterns, tank geometry and levels, pump curves and the current control setpoints. Any version from the last few years is fine — the stress test shows the model absorbs a closed pipe or a tank with half the real volume. | It is the prior: every pipe's length, diameter and roughness, every tank's turnover, every junction's demand pattern. We never build a model; we calibrate theirs. |
| **Grab-sample log, ≥ 6 months** | Date, time of day, sample location, free chlorine (mg/L), or total chlorine (mg/L) for a chloraminated system, never both in one log, method (DPD colorimeter / test strip / lab). A spreadsheet export or photos of the log sheet both work. Twelve months is better (seasons). | The only measurement the method uses. Time of day matters: the model scores each reading against the simulated value at that hour. |
| **Sample-tap list** | For each location that appears in the log: what it is (hydrant, blow-off, dedicated sampling station, customer tap), its address or coordinates, and how long they flush before sampling. | To map each tap to a junction (section 2) and to know which readings carry premise-plumbing bias. |
| **Plant dose** | Free chlorine leaving the plant or entering the system: the daily log if they have it, the monthly average if not. Also the disinfectant (free chlorine / chloramine) and any booster stations. | The dose is a calibrated axis (±10% by default); if it is logged we fix it and the band narrows. Chloramine decays differently and is read as total chlorine: a chloraminated system runs the separate chloramine mode (`python -m residualmap.pilot ... --disinfectant chloramine --dose <total chlorine> --ph <pH> --cl2n <ratio>`, journal task 12). With the plant's pH and chlorine-to-ammonia ratio it runs the model held to task 12's bars (model (c)); without them it runs the same grid with a uniform prior (model (b)), reported in task 12 with no bar. |
| **Plant water temperature, monthly** | The water temperature leaving the plant or at the entry point, as a monthly mean, in a plant log with the columns `month,temp_C,toc_mgL,dose_mgL` (`docs/example_plant_log.csv`; TOC and dose may be left empty). | Chlorine decays faster in warm water: 1.8 to 4.3 times faster per 10 C in the literature (section 5). With the log, each grab sample is explained at its own month's temperature and the next month is predicted at its own (`--plant plant_log.csv`; tested in simulation and not adopted, see section 5); without it the model refits one temperature-blind rate every three months, which is still the default. |
| **Events** | Main breaks, flushing programmes, dose changes, tank cleaning, source changes, with dates. | Those days are excluded from calibration and from the score. |
| **Twenty minutes** | Where the water comes from, what changed in the network since the model was built, which taps they trust, when in the day they sample. | Most structural errors in a `.inp` are known to the operator and unknown to the file. |

What we do not need: SCADA, online analysers, smart meters, customer data, or anyone's time beyond the sampling they already do. Data stays on our laptop; nothing is uploaded to a service. The utility is named in nothing without written agreement.

## 2. How sample taps become junctions

Every reading has to be attached to a junction ID in the `.inp`. The tap map (`tap_id, junction_id, flush_min, notes`) is built once, together with the operator, and versioned with the data.

1. **Hydrants, blow-offs, sampling stations.** Locate on the map, take the nearest junction on the main the tap is connected to (not the nearest junction by straight line — a tap on a 4-inch line is not the 12-inch transmission main 20 m away). Hydrant leads are short; no correction.
2. **Customer taps.** The reading is after a service line and premise plumbing that may hold water for hours: it reads low relative to the main. Flushing until the temperature steadies (typically 2–5 min) brings it close to the main. Taps with no or unknown flushing are kept in the log but flagged (`notes`) and excluded from the score in the first pilot; they are shown to the model as readings with the extra noise they carry.
3. **Taps that fall between junctions.** Attach to the downstream junction of the pipe (water arrives from upstream; the downstream node is closer in age). If the pipe is long (> 300 m) and the model is fine-grained enough to add a node, add a junction to the `.inp` at the tap's chainage instead — a two-line edit.
4. **Taps with no junction.** Some taps sit on pipes the model does not have (a new subdivision, a private line). These readings are logged, dropped from calibration with a warning (`load_log` prints them), and become the first thing to fix in the model.
5. **Hour.** The log's time of day, rounded to the hour. If only the date is known, the operator's usual sampling hour is used and the reading is flagged.

The map is checked once by running the model: a tap mapped to a junction whose modelled water age is very different from its neighbours' usually means the wrong pipe was chosen.

## 3. What the utility gets back every month

One page (the app's PDF, `app.py`), the first working day of the month, from last month's readings and the `.inp`:

- **The chlorine map at the daily minimum** — every junction, median and 90% band, with last month's samples marked.
- **Junctions likely below the minimum residual** — `P(daily minimum < 0.2 mg/L) > 0.5`, listed by junction with the hour the minimum is expected. This is the list that matters for compliance and for where to look for a problem.
- **When the network is worst** — the hour-of-day bar chart. On Net3 the worst hour is midnight and nobody samples at midnight; on Net2 it is 04:00–08:00 and 18:00.
- **Next month's route** — K sites (their usual number) with an hour and a one-sentence reason each, chosen so the readings tell the model the most about the night; plus 2–3 **rotating validation sites** off the usual route (section 4).
- **What the model learned** — the calibrated decay rates and the demand and dose multipliers, in plain words ("your old pipes are eating chlorine about 40% faster than the new ones"; "actual demand runs about 10% above the model").
- **The score so far** — the validation table of section 4, updated every month. The utility sees the same numbers we do.

The operator changes nothing about how they sample, except that the route is suggested rather than habitual and includes the two or three validation sites.

## 4. The validation — exactly

The claim to test is not "the map looks plausible"; it is that the model predicts readings it has not seen, at places it has not seen, with a band that is right about how often it is wrong.

**Hold-out in time.** For each month *m* from the fourth month on: fit on the three months before *m* only, predict every reading of month *m* at its own junction and hour, then compare. Nothing from month *m* or later touches the fit. Reported per month and pooled.

**Two kinds of held-out readings.**
- *Seen taps* — the fixed route. These test the month-to-month forecast: can the model, given last quarter's readings, say what this tap will read this month?
- *New taps* — the 2–3 rotating validation sites, junctions the model has never had a reading from. These test the map: the claim that we can say what the chlorine is where nobody sampled. This is the number that matters, and it is why the protocol adds rotating sites: a fixed route alone can never validate the map.

**Metrics** (all in `outputs/pilot/validation_<name>.csv`):
- **RMSE and MAE** of the predicted median against the reading, mg/L.
- **Coverage** of the 50 / 80 / 90 / 95% bands: the share of readings inside each band. A calibrated model gives ≈ 0.50 / 0.80 / 0.90 / 0.95. Under-coverage means the band is too narrow (the model is over-confident); over-coverage means it is too wide to be useful.
- **Band width**: median width of the 90% band, mg/L — coverage is cheap if the band is a mile wide.
- **The violation call**: of the readings below the minimum residual, the share the model had flagged (`P(below) > 0.5`) — recall — and the number of false alarms.
- **Baselines an operator could run by hand**: *persistence* (the last reading at the same tap) and the *network mean* (mean of the last three months' readings). For a new tap, persistence is undefined and falls back to the network mean — which is exactly today's practice.

**Success criteria for the pilot**, agreed in advance and reported whether or not they are met:

| | seen taps | new taps |
|---|---|---|
| RMSE | ≤ 0.10 mg/L and better than persistence | ≤ 0.12 mg/L and better than the network mean |
| 90% band coverage | 0.85–0.95 | 0.80–0.95 |
| violation recall | ≥ 0.8 of readings below 0.2 mg/L flagged | ≥ 0.8 |

**What we exclude and say so**: readings on event days (section 1), readings from unflushed customer taps, readings with no time of day. The count of exclusions is in the report.

**What the synthetic version gives**, so the format is fixed before real data arrives — `python -m residualmap.pilot --synthetic Net3`, six months of ten route samples plus three rotating sites, the same network with a different operating truth every month (bulk decay ±20%, demand ±15% globally and per node, dose ±10%), fit on the previous three months: pooled over the three held-out months, 39 readings — RMSE **0.077** mg/L (persistence 0.149, network mean 0.288); coverage 0.54 / 0.87 / 0.92 / 1.00 at the 50 / 80 / 90 / 95% levels. On the seven readings at taps the model had never seen: RMSE **0.088** (persistence and network mean 0.202), 90% coverage 0.86. Of the seven held-out readings below 0.2 mg/L, six were flagged, with two false alarms. These are synthetic numbers; the point of the pilot is to replace them.

## 5. What could go wrong, and what we do about it

- **Chloraminated system.** Different chemistry and a different measurement (total chlorine). Since task 12 the chloramine mode has its own decay grid, dose axis and likelihood scale, a 0.5 mg/L total chlorine default threshold (a common utility operating target, not a California rule; California requires a detectable residual) and, in the app, a nitrification watch labelled as literature thresholds (`python -m residualmap.pilot ... --disinfectant chloramine`, `docs/example_grab_log_total.csv` for the log format). In simulation against EPA's chloramine chemistry its model with the pH and Cl2:N prior (model (c), `--ph` and `--cl2n`) met its pre-registered bars where they apply on Net3 and Net2 (journal, task 12; the uniform-prior model (b) the pilot runs without them met the same bars when they were applied to it for information), but it has never seen real data, its band is too narrow under the straddle rule, and two things it cannot represent should be asked about: where chloraminated water blends with chlorinated wells, the mix can cross breakpoint (Cl2:N 7.6 to 9.5:1), and an active nitrification episode is flagged by the watch, not modelled. A free-chlorine pilot is still the simpler first test.
- **Seasonal decay.** Bulk chlorine decay rises about 1.8 to 4.3 times per 10 C in the literature: an Arrhenius E/R of about 4660 to 12,100 K, from the compiled range of Powell et al. 2000 (Water Res 34:117, 4660 to 9560 K), the prior centre of Cejas, Diaz & Gonzalez 2026 (Water 18:1390, 8000 K) and Blokker et al. 2014 (about 12,100 K). The model tested in task 10 uses E/R 5000, 8000 and 12000 K, 1.8 to 4.2 times per 10 C. How pipe-wall decay responds to temperature is not established (Cejas 2026: no wall model includes temperature; Lee et al. 2014 give only the direction). A three-month window that ignores temperature follows the seasons with a lag, and the calibrated `kb` by month is itself a result worth showing. With the plant's monthly water temperature (`--plant`), the model can explain each reading at its own month's temperature, but in simulation (journal, task 10) that model was not better than the temperature-blind one on held-out readings. Fitted in winter and asked for summer it was worse; fitted in summer and asked for December it had the better map RMSE but lower recall. The main cause is that it assumed wall decay rises with temperature as steeply as bulk decay. It stopped at its pre-registered stop rule and is not the default. Ask for the plant temperature log anyway: it is what any later temperature model will need, and a pilot's own seasons are the real test. For a ky4-sized network, bin temperatures to 1 C and use E/R 8000 K only, or the 75-run decay grid, to keep the bank affordable.
- **The model is badly out of date.** The stress test says the method survives a closed pipe and a wrong tank; it does not survive a `.inp` with the wrong source or a missing pressure zone. The twenty-minute call catches that. If the file is beyond repair, the pilot becomes "help the district get a working model", which is also worth doing.
- **Readings cluster at one or two hours.** Fine for the forecast; weak for the night claim. The route's hour suggestions spread them.
- **Too few readings below 0.2 mg/L to score recall.** Good news for the district; we then report the calibration numbers and the daily-minimum map as the deliverable, and the violation recall stays "not testable this quarter".
- **Anything real that contradicts a synthetic result here** is reported as such. The README's limitations section is where it goes.

## 6. First contact — the checklist

1. Do you have an EPANET model of the system, and could you share the `.inp`? Who built it and when?
2. How many chlorine grab samples a month, at which locations, at what time of day? Can we have six to twelve months of the log?
3. What is the dose leaving the plant, and is it logged?
4. Free chlorine or chloramine? If chloramine: is the log free or total chlorine (total is what the chloramine mode reads), what are the plant's pH and chlorine-to-ammonia ratio, and does chloraminated water blend with chlorinated water anywhere? Any boosters?
5. What has changed in the system since the model was built?
6. Would you be willing to add two or three sampling sites a month that we suggest, for three months?
7. Who should receive the one-page report, and does it help to have it before the monthly operations report?
