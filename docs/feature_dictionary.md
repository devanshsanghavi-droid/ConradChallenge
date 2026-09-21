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