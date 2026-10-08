"""
capabilities.py: what the model accounts for, with every result number read from the committed outputs (iteration 4, task 14).

A water-district general manager asked whether the model accounts for water age, organics that impart chlorine
demand, the type of disinfectant and water temperature.  This module holds the answer in one place:

  * accounts_for()  the app's 'What this model accounts for' table: one row per point, saying what the model does,
                    how it was tested (in simulation only) and what it does not do or what was tested and not adopted;
  * report_lines()  the PDF report's chemistry lines: disinfectant, measured species, threshold and its source,
                    how the seasons are handled, and 'validated in simulation only';
  * grid_edge_warning()  the app's warning when the calibration sits on the edge of the decay grid;
  * SOURCES         every result number the table quotes, as (file under the repo, key path into its JSON).  The
                    table is filled from these at run time, and `python -m residualmap.numbers` checks that each one
                    still resolves, so the app never quotes a result that is not in outputs/.  Settings and bars (the
                    0.85 bar, the 5% starting-water rule, the 0.10 per day floor) are fixed text.

Wording rules (plan, task 14): 'accounts for' only where there is a term in the model and a test; 'absorbed into the
calibrated decay rates' where there is not; 'not modelled' for nitrification as a process, blending and pH.  Nothing
here claims real-data validation: every test is a simulation.
"""
from __future__ import annotations

import json
import os

import pandas as pd

from .chemistry import CHLORAMINE, DEFAULT_THRESHOLD_MGL, FREE_CHLORINE, SPECIES, THRESHOLD_NOTE

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# name -> (file under the repo, key path).  An int in a key path indexes a list; ('count_material', ...) and the
# other derived values are computed in numbers() from the keys listed here.
SOURCES: dict[str, tuple[str, tuple]] = {
    # water age (task 9): default truth, pooled over seeds
    "age_cov_max_Net3": ("outputs/chem/water_age_Net3.json", ("by_truth", "default", "band_coverage_daily_max")),
    "age_cov_max_Net2": ("outputs/chem/water_age_Net2.json", ("by_truth", "default", "band_coverage_daily_max")),
    "age_cov_max_ky4": ("outputs/chem/water_age_ky4.json", ("by_truth", "default", "band_coverage_daily_max")),
    "age_band_label_Net3": ("outputs/chem/water_age_Net3.json", ("band_label",)),
    "wall_share_mae_Net3": ("outputs/chem/water_age_Net3.json", ("by_truth", "default", "loss_split", "wall_share_mae")),
    # the file's daily-maximum ages that are lower bounds (more than 5% starting water at their oldest hour), and the
    # junction count (junction-seeds over seeds)
    **{f"age_lb_max_{n}": (f"outputs/chem/water_age_{n}.json", ("nominal_age", "n_lower_bound_daily_max"))
       for n in ("Net3", "Net2", "ky4")},
    **{f"age_n_js_{n}": (f"outputs/chem/water_age_{n}.json", ("by_truth", "default", "n_junction_seeds"))
       for n in ("Net3", "Net2", "ky4")},
    **{f"age_n_seeds_{n}": (f"outputs/chem/water_age_{n}.json", ("by_truth", "default", "n_seeds"))
       for n in ("Net3", "Net2", "ky4")},
    # organics (task 11): second-order (Clark) truth O2, TOC schedule S, Net3
    "toc_readings_ratio_Net3": ("outputs/chem/summary_organics_Net3.json",
                                ("acceptance", "G", "by_case", "O2_S_rolling_all_months_readings", "ratio")),
    "toc_first_storm_ratio_Net3": ("outputs/chem/summary_organics_Net3.json", ("acceptance", "T2", "ratio")),
    "toc_first_storm_pass_Net3": ("outputs/chem/summary_organics_Net3.json", ("acceptance", "T2", "pass")),
    "toc_adopted_Net3": ("outputs/chem/summary_organics_Net3.json", ("acceptance", "adoption", "passes_on_this_network")),
    "toc_adopted_Net2": ("outputs/chem/summary_organics_Net2.json", ("acceptance", "adoption", "passes_on_this_network")),
    # two-reactant audit (task 13)
    "audit_cells": ("outputs/chem/summary_audit_Net3.json", ("acceptance", "materiality", "cells")),
    "audit_oldest_bias_Net3": ("outputs/chem/summary_audit_Net3.json", ("acceptance", "age_decile", "bias_2ra")),
    "lowkb_candidate": ("outputs/chem/summary_audit_Net3.json", ("acceptance", "lowkb_recommendation", "candidate_for_default")),
    # chloramine (task 12): model (c), the mode's model
    "ca_recall_straddle8_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                                 ("pooled", "main", "c", "straddle_min", "8", "recall")),
    "ca_low_straddle8_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                              ("pooled", "main", "c", "straddle_min", "8", "n_true_viol")),
    "ca_cov90_straddle8_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                                ("pooled", "main", "c", "straddle_min", "8", "coverage90")),
    "ca_cov90_random_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                             ("acceptance", "A2_coverage90_c_random_mean_over_n", "value")),
    "ca_cov90_random_Net2": ("outputs/chloramine/summary_chloramine_Net2.json",
                             ("acceptance", "A2_coverage90_c_random_mean_over_n", "value")),
    # free-chlorine settings on the chloraminated network: the posterior on their lowest bulk rate after 15 samples
    # (pooled over both sampling rules), and their false alarms at 8 random samples against the mode's
    "ca_free_floor_mass_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                                ("acceptance", "A5_free_grid_reported", "kb_floor_mass_n15")),
    "ca_false_alarms_mode_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                                  ("pooled", "main", "c", "random", "8", "false_alarms")),
    "ca_false_alarms_free_Net3": ("outputs/chloramine/summary_chloramine_Net3.json",
                                  ("pooled", "main", "a", "random", "8", "false_alarms")),
    "ca_experimental_Net3": ("outputs/chloramine/summary_chloramine_Net3.json", ("acceptance", "experimental")),
    "ca_experimental_Net2": ("outputs/chloramine/summary_chloramine_Net2.json", ("acceptance", "experimental")),
    # temperature (task 10, model M; task 10b, model M2)
    "temp_m_julsep_rmse_Net3": ("outputs/chem/summary_season_Net3.json", ("acceptance", "A4", "rmse_M")),
    "temp_b0_julsep_rmse_Net3": ("outputs/chem/summary_season_Net3.json", ("acceptance", "A4", "rmse_B0")),
    "temp_m_stop_Net3": ("outputs/chem/summary_season_Net3.json", ("acceptance", "A7", "stop")),
    "temp_m_stop_Net2": ("outputs/chem/summary_season_Net2.json", ("acceptance", "A7", "stop")),
    "temp_m2_ratio_Net3_W1": ("outputs/chem/summary_season2_Net3.json", ("acceptance", "W1", "A3", "bootstrap", "ratio")),
    "temp_m2_ratio_Net3_W2": ("outputs/chem/summary_season2_Net3.json", ("acceptance", "W2", "A3", "bootstrap", "ratio")),
    "temp_m2_ratio_Net2_W1": ("outputs/chem/summary_season2_Net2.json", ("acceptance", "W1", "A3", "bootstrap", "ratio")),
    "temp_m2_ratio_Net2_W2": ("outputs/chem/summary_season2_Net2.json", ("acceptance", "W2", "A3", "bootstrap", "ratio")),
    "temp_m2_stop_Net3_W1": ("outputs/chem/summary_season2_Net3.json", ("a7_summary", "W1", "stop")),
    "temp_m2_stop_Net3_W2": ("outputs/chem/summary_season2_Net3.json", ("a7_summary", "W2", "stop")),
    "temp_m2_stop_Net2_W1": ("outputs/chem/summary_season2_Net2.json", ("a7_summary", "W1", "stop")),
    "temp_m2_stop_Net2_W2": ("outputs/chem/summary_season2_Net2.json", ("a7_summary", "W2", "stop")),
    "temp_m2_adopted_Net3": ("outputs/chem/summary_season2_Net3.json", ("adoption", "passes_A7_both_worlds_on_this_network")),
    "temp_m2_adopted_Net2": ("outputs/chem/summary_season2_Net2.json", ("adoption", "passes_A7_both_worlds_on_this_network")),
}


def resolve(name: str, root: str = REPO):
    """The value SOURCES[name] points at (KeyError, IndexError or OSError when it no longer resolves)."""
    rel, keys = SOURCES[name]
    with open(os.path.join(root, rel)) as fh:
        x = json.load(fh)
    for k in keys:
        x = x[k]
    return x


def numbers(root: str = REPO) -> dict:
    """Every SOURCES value plus the derived ones the table quotes.  Raises if a source no longer resolves."""
    v = {name: resolve(name, root) for name in SOURCES}
    cells = v.pop("audit_cells")
    v["audit_n_cells"] = len(cells)
    v["audit_n_within"] = sum(1 for c in cells if not c["material"])
    v["ca_found_straddle8_Net3"] = int(round(v["ca_recall_straddle8_Net3"] * v["ca_low_straddle8_Net3"]))
    m2 = [v[f"temp_m2_ratio_{n}_{w}"] for n in ("Net3", "Net2") for w in ("W1", "W2")]
    v["temp_m2_ratio_min"], v["temp_m2_ratio_max"] = min(m2), max(m2)
    v["temp_m2_ratio_min_world"] = min((v[f"temp_m2_ratio_{n}_{w}"], w) for n in ("Net3", "Net2") for w in ("W1", "W2"))[1]
    v["temp_m2_n_stop"] = sum(bool(v[f"temp_m2_stop_{n}_{w}"]) for n in ("Net3", "Net2") for w in ("W1", "W2"))
    v["temp_adopted"] = bool(v["temp_m2_adopted_Net3"] and v["temp_m2_adopted_Net2"])
    v["toc_adopted"] = bool(v["toc_adopted_Net3"] and v["toc_adopted_Net2"])     # task 11's own adoption flag
    for n in ("Net3", "Net2", "ky4"):
        v[f"age_n_junctions_{n}"] = int(round(v[f"age_n_js_{n}"] / v[f"age_n_seeds_{n}"]))
    v["ca_experimental"] = bool(v["ca_experimental_Net3"] or v["ca_experimental_Net2"])
    return v


COLUMNS = ("point", "in the model", "tested, in simulation only", "not done, or tested and not adopted")


def accounts_for(root: str = REPO) -> pd.DataFrame:
    """The app's 'What this model accounts for' table: water age, organics and chlorine demand, chlorine type and
    temperature, in that order (the manager's four points), every number from outputs/ via SOURCES."""
    n = numbers(root)
    ca_label = "experimental, " if n["ca_experimental"] else ""
    rows = [
        ("Water age",
         "Yes, through EPANET: every simulation moves the water through each pipe and tank of your file hour by hour, "
         "so water age is inside every prediction. Shown since iteration 4: a water-age map with its range across "
         "demand and roughness errors, and where each junction's chlorine is lost (in the water or at the pipe walls).",
         f"Against a simulated truth's own water age, the range held the true daily-maximum age at "
         f"{n['age_cov_max_Net3']:.3f} of Net3's junctions, {n['age_cov_max_Net2']:.3f} of Net2's and "
         f"{n['age_cov_max_ky4']:.3f} of ky4's, so it is called a {n['age_band_label_Net3']}, not a 90% band. The "
         f"calibrated wall share of the loss was off from the truth's by {n['wall_share_mae_Net3']:.3f} per junction "
         f"on average on Net3.",
         f"Ages come from a 7-day run: where more than 5% of a junction's water is still the run's starting water, the "
         f"age is a lower bound and is shown as 'at least'. At their oldest hour that holds at "
         f"{n['age_lb_max_Net3']} of Net3's {n['age_n_junctions_Net3']} junctions, {n['age_lb_max_Net2']} of Net2's "
         f"{n['age_n_junctions_Net2']} and {n['age_lb_max_ky4']} of ky4's {n['age_n_junctions_ky4']}, so the coverage "
         f"beside this is measured mostly on lower bounds. A longer warm-up is proposed, not built."),
        ("Organics and chlorine demand",
         "Absorbed into the calibrated decay rates: decay in the water (where organics use up chlorine) and at the pipe "
         "walls is refitted from your last three months of samples. Plant TOC is not an input.",
         f"A model given the plant's monthly TOC was tested on 12 simulated months against chlorine-organics chemistry "
         f"first order cannot represent, with an assumed TOC schedule (the first storm is an assumed November doubling) "
         f"and assumed constants in the truth: held-out readings were slightly better on Net3 (RMSE "
         f"{n['toc_readings_ratio_Net3']:.3f} of today's), but at the first storm its map RMSE was "
         f"{n['toc_first_storm_ratio_Net3']:.3f} of today's against a bar of 0.85, and it stopped at its "
         f"pre-registered stop rule. Against two-reactant chlorine-organics chemistry (constants from one water, via a "
         f"secondary source), today's model stayed within 10% in {n['audit_n_within']} of {n['audit_n_cells']} test "
         f"cells.",
         f"{'Tested, not adopted' if not n['toc_adopted'] else 'Adopted'}: "
         f"{'TOC is not an input' if not n['toc_adopted'] else 'TOC is an optional logged input'}. Part of the TOC "
         f"model's gain came from decay rates below today's grid floor (0.10 per day); a grid with lower rates "
         f"{'is a candidate for the default, not adopted yet' if n['lowkb_candidate'] else 'was tested'}. Under "
         f"two-reactant chemistry today's map reads {abs(n['audit_oldest_bias_Net3']):.3f} mg/L low in Net3's oldest "
         f"tenth of junctions."),
        ("Chlorine type",
         f"Free chlorine (the default) or chloramine ({ca_label}a separate mode; the two are never mixed). The "
         f"chloramine mode reads total chlorine and has its own decay ranges, dose axis, likelihood scale and a 0.5 mg/L "
         f"default threshold, plus a nitrification watch (literature thresholds, not validated). Gas chlorine and "
         f"hypochlorite give the same free chlorine; they differ mainly through pH, which is not modelled and is "
         f"absorbed into the calibrated decay rates.",
         f"In simulation, against EPA's chloramine chemistry in EPANET-MSX with an assumed wall term: with 8 "
         f"straddle-chosen samples on Net3 the mode found {n['ca_found_straddle8_Net3']} of "
         f"{int(n['ca_low_straddle8_Net3'])} low junction-days, and its random-sample 90% band held the truth "
         f"{n['ca_cov90_random_Net3']:.3f} (Net3) and {n['ca_cov90_random_Net2']:.3f} (Net2) of the time. Free-chlorine "
         f"settings on the same network, fitted to 15 samples, put {n['ca_free_floor_mass_Net3']:.2f} of their "
         f"posterior on their lowest decay rate; at 8 random samples they raised "
         f"{int(n['ca_false_alarms_free_Net3'])} false alarms against the mode's {int(n['ca_false_alarms_mode_Net3'])}.",
         f"Under the straddle rule the chloramine band is too narrow ({n['ca_cov90_straddle8_Net3']:.3f} at 8 samples "
         f"on Net3). Not modelled: blending with chlorinated water (breakpoint), nitrification itself, and a seasonal "
         f"chloramine test."),
        ("Temperature",
         "Absorbed into the calibrated decay rates, refitted from your last three months of samples, so the model "
         "follows the seasons with a lag. Water temperature is not an input of the decay model (in chloramine mode the "
         "nitrification watch uses this month's temperature only to mark warm water).",
         f"Two temperature-aware models were fitted in one season and asked for another over 12 simulated months, with "
         f"an assumed temperature schedule (10 to 20 C) and an assumed pipe-wall temperature response in the truth. The "
         f"first assumed pipe-wall decay speeds up with temperature as steeply as decay in the water: its July-to-"
         f"September map from January-to-March samples had RMSE {n['temp_m_julsep_rmse_Net3']:.3f} mg/L against today's "
         f"{n['temp_b0_julsep_rmse_Net3']:.3f} on Net3. The second learned the wall's response and predicted held-out "
         f"readings better ({n['temp_m2_ratio_min']:.3f} to {n['temp_m2_ratio_max']:.3f} of today's RMSE; the low end "
         f"in {_WORLD[n['temp_m2_ratio_min_world']]}), but its stop rule triggered in {n['temp_m2_n_stop']} of the 4 "
         f"simulated network and wall-response cases. Neither was adopted.",
         f"{'Tested, not adopted' if not n['temp_adopted'] else 'Adopted'}. The deciding unknown is how real pipe walls "
         f"respond to temperature; a year of real grab samples with plant water temperatures would settle it."),
    ]
    return pd.DataFrame(rows, columns=list(COLUMNS))


# task 10b's two simulated truth worlds, as the table names them
_WORLD = {"W1": "the world with a weak wall response, an assumption",
          "W2": "a constructed world where the walls respond to temperature as the water does"}

ACCOUNTS_FOR_NOTE = ("Every test in this table is a simulation, on EPANET's example networks Net3 and Net2 (and, for water "
                     "age, the KYPIPE network ky4), against hidden chemistry with assumed constants and schedules: no real "
                     "grab-sample data has been scored yet. Every result number is read from outputs/chem and "
                     "outputs/chloramine when the app runs (settings and bars are fixed text); the journal "
                     "(docs/iteration3_journal.md) has the full results.")

SEASON_LINE = ("Season: water temperature and organics (TOC) are not inputs; their effect is absorbed into the decay rates, "
               "refitted from the last three months of samples (temperature and TOC models were tested in simulation and "
               "not adopted).")
SEASON_LINE_CA = ("Season: water temperature and organics (TOC) do not enter the decay model (this month's water temperature "
                  "is used only by the nitrification watch); their effect is absorbed into the decay rates, refitted from "
                  "the last three months of samples. No seasonal chloramine test has been run.")
SIMULATION_ONLY_LINE = "Validated in simulation only: no real grab-sample data has been scored yet."
EDGE_WARN = 0.5        # posterior mass on one edge of the bulk or wall decay axis that triggers the grid-edge warning


def threshold_source(chloramine: bool, threshold: float) -> str:
    """The threshold and where it comes from: the default's note (a common operating minimum or target, not a
    California rule), or 'set by you' for any other value."""
    dis = CHLORAMINE if chloramine else FREE_CHLORINE
    if abs(float(threshold) - DEFAULT_THRESHOLD_MGL[dis]) < 1e-9:
        return THRESHOLD_NOTE[dis]
    return f"{float(threshold):g} mg/L {SPECIES[dis]}: set by you (California requires a detectable residual)"


def report_lines(chloramine: bool, threshold: float, ph: float | None = None, cl2n: float | None = None,
                 n_watch: int | None = None) -> list[str]:
    """The PDF report's chemistry lines, one string per line, no em or en dashes.  No temperature or kb20 field: no
    temperature model was adopted (tasks 10 and 10b), so the report says how the seasons are handled instead."""
    if chloramine:
        prior = f", prior from plant pH {ph:g} and Cl2:N {cl2n:g}" if ph is not None and cl2n is not None else ""
        watch = (f"; nitrification watch: {n_watch} junction{'' if n_watch == 1 else 's'} (literature thresholds, not "
                 f"validated)" if n_watch is not None else "")
        dis = (f"Disinfectant: chloramine. Measured species: TOTAL chlorine (mg/L as Cl2); every reading and map value is "
               f"total chlorine, in the chloramine mode (its own decay ranges{prior}){watch}.")
    else:
        dis = ("Disinfectant: free chlorine. Measured species: free chlorine (mg/L as Cl2); every reading and map value "
               "is free chlorine.")
    return [dis, f"Minimum residual: {threshold_source(chloramine, threshold)}.",
            SEASON_LINE_CA if chloramine else SEASON_LINE, SIMULATION_ONLY_LINE]


def grid_edge_warning(edge: dict, params, chloramine: bool, share: float = EDGE_WARN) -> str | None:
    """The app's grid-edge warning: a plain sentence when more than `share` of the calibration (simgp.grid_edge_mass of
    the fitted model, keys kb_low, kb_high, kw_low, kw_high) sits on the lowest or highest bulk or wall decay rate of
    the grid in `params` (the model's members, (kb, kw, gamma, demand, roughness) tuples).  The water may need a rate
    the grid does not have.  For free chlorine at the bulk floor it adds what the simulations found there (tasks 10b,
    11 and 13: the fit often sat there and today's map was biased low, with the cause not settled).  None when no edge
    holds more than `share`."""
    parts = []
    for i, key, name, unit in ((0, "kb", "bulk decay rate", "per day"), (1, "kw", "wall decay rate", "m/day")):
        vals = sorted({float(p[i]) for p in params})
        for side, v in (("low", vals[0]), ("high", vals[-1])):
            mass = float(edge.get(f"{key}_{side}", 0.0))
            if mass > share:
                parts.append(f"{mass:.0%} of the calibration sits on the {'lowest' if side == 'low' else 'highest'} "
                             f"{name} the grid allows ({v:g} {unit})")
    if not parts:
        return None
    floor_note = (" In simulation the free-chlorine fit often sat at this floor of 0.10 per day (journal, tasks 10b, 11 "
                  "and 13), and today's map was biased low in those simulations; the cause is not settled. A grid with "
                  "lower rates was tested and is a candidate, not yet the default."
                  if not chloramine and float(edge.get("kb_low", 0.0)) > share else "")
    return (f"Grid edge: {'; '.join(parts)}. The water may need a rate outside the grid, so treat the map with more "
            f"caution there.{floor_note}")
