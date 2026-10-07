# Feature dictionary

Every input the surrogate may use. None of them requires a chlorine measurement; all come from the operator's EPANET `.inp` plus one hydraulic/age run of it.

**CORE set (used by the GP at 3–15 samples):** `age_h`, `emb0`, `emb1`, `emb2`, `path_wall_index`, `dist_src_km`

| feature | physical meaning |
|---|---|
| `age_h` | Water age at the sampling hour (nominal model). First-order decay makes ln C ~ linear in age. |
| `age_mean_h` | Mean water age over the last simulated day. |
| `age_max_h` | Max water age over the last day — stagnation flag for dead ends and tank-fed zones. |
| `emb0/emb1/emb2` | 3-D classical-MDS embedding of pipe-length-weighted hydraulic distance; puts hydraulically close nodes close. |
| `dist_src_km` | Hydraulic distance to nearest reservoir/source. |
| `dist_tank_km` | Hydraulic distance to nearest tank (tanks release old, low-chlorine water). |
| `path_len_km` | Length of the shortest pipe path from the nearest source. |
| `path_mean_rough` | Length-weighted mean Hazen-Williams C along that path (low C = rough = old = more wall decay). |
| `path_min_rough` | Roughest pipe on the path. |
| `path_mean_diam_m` | Length-weighted mean diameter on the path (wall decay ~ 4*kw/D: small pipes lose chlorine faster). |
| `path_min_diam_m` | Narrowest pipe on the path. |
| `path_sum_L_over_D` | Sum of length/diameter along the path — the wall-contact exposure integral for first-order wall decay. |
| `path_wall_index` | Sum of (L/D) * roughness_factor(C) — same integral, weighted by the old-pipe hypothesis. |
| `elevation_m` | Junction elevation. |
| `demand_lps` | Base demand (people served proxy; high-demand nodes pull fresh water). |
| `degree` | Number of connected links (1 = dead end). |
| `pressure_mean_m` | Mean nominal pressure over the day. |
| `pressure_min_m` | Min nominal pressure over the day. |
| `inc_vel_mean_ms` | Mean velocity of incident pipes (low = stagnant). |
| `inc_vel_min_ms` | Min velocity of incident pipes. |
| `inc_abs_flow_m3s` | Mean |flow| of incident pipes. |
| `inc_frac_sloshing` | Fraction of incident pipes whose flow reverses during the day (mixing zones). |

## What the calibrated simulator uses on top of this
Everything in the `.inp` that EPANET itself uses: pipe lengths, diameters, roughness, junction demands and their daily patterns, tank geometry and levels, pump curves and controls, valve settings, source locations. The three unknowns it calibrates from grab samples are bulk decay `kb` (1/day), wall decay `kw` (m/day) and `gamma`, the strength of the 'rougher pipe = faster wall decay' hypothesis.

Roughness → wall-decay hypothesis: factor = 2^(-gamma·(C−130)/30). With gamma=1: C=110 → 1.6×, C=130 → 1×, C=199 → 0.2×.

## Water age and where chlorine is lost (iteration 4, task 9; `residualmap/age.py`)
Outputs shown to the operator, not new GP inputs: the GP's inputs are unchanged. All come from EPANET runs of the operator's own `.inp`, plus the grab samples for the last three. Tested against a simulated truth only (`outputs/chem/water_age_<net>.json`).

| quantity | physical meaning |
|---|---|
| `age band` (`hydraulic_age_band`) | Water age (hours since the water left a source) from EPANET AGE runs at the grid's 9 hydraulic settings: global demand x0.85, x1, x1.15 times Hazen-Williams C x0.9, x1, x1.1. Min, median and max per junction and hour. The (x1, x1) member is the nominal age, `age_h` above. Like every age here, it is the age on the last day of a 7-day run that starts with the file's water at age 0. |
| `initial_share` (`simulate.initial_water_share`) | The share of a junction's water on the last day (per hour) that is still the water the 7-day run started with in its pipes, junctions and tanks, not water that entered from a source during the run. From a second AGE run with every junction's and tank's starting age raised by 1000 h: ages mix linearly, so each junction's age rises by 1000 h times this share. Above 5% the junction's age is a lower bound (`AgeBand.lower_bound`: on the day's average for the daily-mean age, at any hour for the daily-max age). Per junction of each example file in `outputs/chem/water_age_junctions_<net>.csv`. |
| range width | Max minus min of the 9 members' daily-maximum age at a junction (the app's second water-age map). Called a "range across demand and roughness errors", not a 90% band, unless a simulated test shows at least 0.85 coverage; errors in single junctions' demands are not among the 9 settings. |
| oldest water | The junction with the largest daily-maximum nominal age. Where more than 5% of its water at that hour is still the run's starting water, the age is a lower bound and is reported as "at least X h" (on every example network). The plan's 160 h cut (an age past 160 h on the last day of the run) catches only the extreme cases of this. |
| `posterior age` (`posterior_age`) | The 9 members' ages weighted by the calibrated model's posterior over the demand and roughness settings (its weights summed over the decay and dose axes): what the grab samples say about the hydraulics, in hours. Its band is the central 90% of that 9-point posterior. |
| `wall_share`, `bulk_share` (`loss_split`) | Of each junction's chlorine loss, the share lost at the pipe walls and in the water, for the calibrated member: L = mean over the day of ln(C_ref / C), with C_ref the same EPANET run with no decay; bulk share = L(wall decay off) / L(as calibrated), wall share = 1 minus that. Exact along a single plug-flow path (both decay terms are first order), approximate where flows mix. Undefined below a 2% loss (junctions next to a source). |
| `nonadditivity` | (L(wall decay off) + L(bulk decay off)) / L(as calibrated): 1 on a single path; its distance from 1 says how approximate the split is where flows mix. |
