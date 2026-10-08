"""
numbers.py: every iteration-4 number quoted in the README and the journal traces to a committed file under outputs/
(iteration 4, task 14).

    python -m residualmap.numbers            # exits 0 when every number is traced or classified, 1 otherwise
    python -m residualmap.numbers --list     # also prints every number that is not traced to a file, with its class
    python -m residualmap.checks --numbers   # the same, from the checks script

Two layers.

1. The registry (REGISTRY).  Each headline number of task 14's own text (the journal's 'Capabilities after iteration
   4' and 'Iteration 4 summary' sections and the README's 'Iteration 4 in one page', their lists of regressions
   included) is tied to one file and key under outputs/, or to one rule over a file (DERIVED, for a range or a share
   the text quotes), and one format, with a few words of the text around it.  The check: the value, formatted and put
   into those words, is in that section.  The app's 'What this model accounts for' table needs no registry: it reads its numbers
   from outputs/ at run time (residualmap/capabilities.py, SOURCES), and every source must resolve.

2. The scan.  Every number in the iteration-4 parts of the README and the journal (the journal from 'Iteration 4,
   chemistry' to its end; the README's iteration-4 sections, its app section and the limitations bullets that name an
   iteration-4 task) is looked up, at its printed precision, among the values in the committed outputs files that its
   section cites by name; if it is not there, among every committed file in outputs/chem, outputs/chem_2ra and
   outputs/chloramine; then among all committed files in outputs/.  A number found in none of them is classified:
   a time of day, a date or year, a citation (volume:page, doi), an identifier (task, addendum, finding, section,
   seed), a run cost (seconds, minutes, MB, GB, cores: these come from run logs and df readings, and the text says so),
   a count derived from a stored rate ('97 of 109' where 97/109 is a stored recall), a setting or literature value
   written in the code (residualmap/*.py, app.py), or an entry of the explicit allowlist ALLOW with its reason.
   Anything left is 'unexplained' and fails the run.

What the scan can and cannot show.  It catches a number that appears in no output at all (a typo, a scratch number, a
number from an earlier run).  It cannot prove that a number came from the right key: the outputs hold over a million
values, so a three-digit decimal often matches some value by chance, and a cited file's pool is smaller but still
large.  An integer is matched at plus or minus 0.5, and an 'X of Y' count is accepted when X/Y is a stored rate, so
its power on counts is low too; the script prints both powers.  The spans it blanks (dates, identifiers, times,
citations, seed ranges) are guarded on both sides so they never take part of a decimal, and a blanked span that
touches a digit is reported as unexplained.  The registry is the exact check, for the numbers task 14 puts in front
of a reader.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOURNAL = os.path.join(REPO, "docs", "iteration3_journal.md")
README = os.path.join(REPO, "README.md")
JOURNAL_START = "## Iteration 4, chemistry (tasks 8 to 14)"
IT4_DIRS = ("outputs/chem/", "outputs/chem_2ra/", "outputs/chloramine/")
SELF_REPORT = "outputs/chem/checks_report.json"      # the saved checks' report: its numbers-scan entry is left out
IT4_TASK = re.compile(r"\b(?:tasks?|Task) (?:8|9|10b?|11|12|13|14)\b|\biteration 4\b|\bIteration 4\b")

MONTHS = ("January|February|March|April|May|June|July|August|September|October|November|December|"
          "Jan|Feb|Mar|Apr|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec")
NUM = re.compile(r"(?<![\w.])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:e[-+]?\d+)?%?(?![\w])")
# Every pattern below that takes digits is guarded on both sides (_L, _R), so it never takes part of a decimal: 'July
# 0.949' leaves 0.949 to be read, '2.5 December' leaves 2.5, 'model 12.5' leaves 12.5.  scan() also reports any blanked
# span that still touches a digit as unexplained.
_L, _R = r"(?<![\d.])", r"(?!\.?\d)"
BLANKS = [   # (class, pattern): spans removed before numbers are read; each match is counted under its class
    ("code or path", re.compile(r"`[^`\n]*`")),
    ("link target", re.compile(r"\]\([^)]*\)")),
    ("time of day", re.compile(rf"{_L}\b\d{{1,2}}:\d{{2}}(?::\d{{2}})?{_R}")),
    ("citation", re.compile(rf"doi:\S+|arXiv:\S+|arxiv:\S+|{_L}\b\d+(?:\(\d+\))?:\d+(?:-\d+)?{_R}|"
                            rf"{_L}\b\d+ (?:CCR|TAC|CFR) \d+(?:\.\d+)*{_R}")),
    ("date", re.compile(rf"{_L}\b\d{{1,2}} (?:{MONTHS})(?: \d{{4}}{_R})?\b|\b(?:{MONTHS}) \d{{1,2}}{_R}(?:, \d{{4}}{_R})?|"
                        rf"{_L}\b\d{{4}}-\d{{2}}(?:-\d{{2}})?{_R}|{_L}\b(?:19|20)\d{{2}}{_R}")),
    ("identifier", re.compile(r"\b(?:[Tt]asks?|[Aa]ddend(?:um|a)|[Ff]indings?|[Ss]ections?|[Ss]teps?|[Ii]terations?|"
                              r"[Tt]able|[Ii]tems?|[Qq]uestions?|[Dd]ay|[Cc]hapters?|[Rr]ows?|[Cc]olumns?|"
                              r"[Ss]eeds?|[Ss]cenarios?|[Mm]odels?|[Cc]ells?|label|labels)"
                              rf"(?: (?:\d+[a-z]?{_R}|\(\w\))(?:(?:, | and | to | or |-)(?:\d+[a-z]?{_R}))*)+")),
    ("seed range", re.compile(rf"\b(?:Net3|Net2|ky4) \d+ to \d+{_R}|{_L}\b\d+ to \d+ \((?:Net3|Net2)\)")),
]
COST = re.compile(r"^\s*(?:s|seconds?|min|minutes?|ms|MB|GB|GiB|KB|kB|cores?|workers?|CPU|EPANET runs)\b")
OF = re.compile(r"(?<![\w.])(\d+) of (\d+)(?![\w.])")


# ----------------------------------------------------------------------------- the registry (task 14's own text)
# (doc, section heading prefix, text, source, format).  text: the words around the number as printed, with {} where
# the number goes.  source: a capabilities.numbers() name, a DERIVED name (a rule over a file, for a range or a share
# the text quotes), or (file, key path); an int in a key path indexes a list.  format: a format spec applied to the
# value ('.3f', 'd', '.0%', 'abs.3f' for |value|, 's' for a string).
CAP, SUM, ONE = "Capabilities after iteration 4", "Iteration 4 summary", "Iteration 4 in one page"
_WA = "outputs/chem/water_age_{}.json"
_S3, _S2 = "outputs/chem/summary_season_Net3.json", "outputs/chem/summary_season_Net2.json"
_O3 = "outputs/chem/summary_organics_Net3.json"
_A3, _A2 = "outputs/chem/summary_audit_Net3.json", "outputs/chem/summary_audit_Net2.json"
_C3, _C2 = "outputs/chloramine/summary_chloramine_Net3.json", "outputs/chloramine/summary_chloramine_Net2.json"
_BR = "outputs/chem/baseline_reproduction.json"
_T14 = ("runs", "post_task14")
_T2 = ("acceptance", "T2")
_MAT14 = ("acceptance", "materiality", "cells", 14)       # Net2, 2ra, straddle_min, 15 samples, RMSE: the material cell
_FA = ("false_alarms_no_bar", 0)                          # Net3, 2ra, random, 8 samples


def _age(net: str, *keys) -> tuple:
    return (_WA.format(net), keys)


def _json(rel: str):
    with open(os.path.join(REPO, rel)) as fh:
        return json.load(fh)


def _pct_range(ratios) -> str:
    """'10 to 14%': the smallest and largest excess of a set of ratios over 1, in whole percent, as the text prints it."""
    lo, hi = (round(100 * (r - 1)) for r in (min(ratios), max(ratios)))
    return f"{lo} to {hi}%"


# Rules over a file, for numbers the text quotes as a range or a share: name -> (files it reads, function).
DERIVED: dict[str, tuple[tuple[str, ...], callable]] = {
    # task 12, low-pH stress, random rule: model (c)'s RMSE over the uniform prior (b)'s at every n
    "ca_lowph_prior_cost": ((_C3,), lambda: _pct_range(
        [_json(_C3)["pooled"]["stress"]["c"]["random"][n]["rmse"] / _json(_C3)["pooled"]["stress"]["b"]["random"][n]["rmse"]
         for n in ("3", "8", "15")])),
    # task 12, Net2: model (b) on the MSX truth against (b) on its first-order twin, RMSE, every random n
    "ca_first_order_cost_Net2": ((_C2,), lambda: _pct_range(
        [_json(_C2)["first_order_cost"]["random"][n]["b_on_msx"]["rmse"] /
         _json(_C2)["first_order_cost"]["random"][n]["b_on_twin"]["rmse"] for n in ("3", "8", "15")])),
    # task 13, Net2: the largest share of fits whose most probable bulk rate sits at the low-rate grid's own floor
    "lowkb_floor_share_max_Net2": ((_A2,), lambda: max(v["map_kb_variant"]["share_map_kb_at_floor"]
                                                       for v in _json(_A2)["lowkb"]["by_truth"].values())),
    # task 13, Net3: today's 90% band in the oldest tenth of junctions (bin 10 of 10), two-reactant truth, all seeds
    "oldest_tenth_cov90_2ra_Net3": (("outputs/chem/error_by_age_2ra_Net3.csv",), lambda: float(
        pd.read_csv(os.path.join(REPO, "outputs/chem/error_by_age_2ra_Net3.csv"), dtype={"seed": str})
        .query("truth == '2ra|B0' and seed == 'all' and bin == 10 and n_bins == 10")["coverage90"].item())),
}


REGISTRY: list[tuple] = [
    # journal, 'Capabilities after iteration 4'
    ("journal", CAP, "daily-maximum age at {}, ", "age_cov_max_Net3", ".3f"),
    ("journal", CAP, "0.791, {} and 0.811 of junctions", "age_cov_max_Net2", ".3f"),
    ("journal", CAP, "0.779 and {} of junctions (Net3, Net2, ky4)", "age_cov_max_ky4", ".3f"),
    ("journal", CAP, "off by {}, 5.62 and 8.11 h RMSE", _age("Net3", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", CAP, "5.82, {} and 8.11 h RMSE", _age("Net2", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", CAP, "5.62 and {} h RMSE", _age("ky4", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", CAP, "wall share is off by {} per junction on Net3", "wall_share_mae_Net3", ".3f"),
    ("journal", CAP, "at their oldest hour {} of 92 junctions on Net3", "age_lb_max_Net3", "d"),
    ("journal", CAP, "87 of {} junctions on Net3", "age_n_junctions_Net3", "d"),
    ("journal", CAP, "{} of 35 on Net2", "age_lb_max_Net2", "d"),
    ("journal", CAP, "29 of {} on Net2", "age_n_junctions_Net2", "d"),
    ("journal", CAP, "{} of 959 on ky4", "age_lb_max_ky4", "d"),
    ("journal", CAP, "799 of {} on ky4", "age_n_junctions_ky4", "d"),
    ("journal", CAP, "(on the day's average {}, 28 and 745)", _age("Net3", "nominal_age", "n_lower_bound_daily_mean"), "d"),
    ("journal", CAP, "average 19, {} and 745)", _age("Net2", "nominal_age", "n_lower_bound_daily_mean"), "d"),
    ("journal", CAP, "19, 28 and {})", _age("ky4", "nominal_age", "n_lower_bound_daily_mean"), "d"),
    ("journal", CAP, "holds the truth at {} on Net3", _age("Net3", "by_truth", "default", "converged_only",
                                                           "band_coverage_daily_max"), ".3f"),
    ("journal", CAP, "on Net3 ({} junction-seeds)", _age("Net3", "by_truth", "default", "converged_only",
                                                         "n_junction_seeds_daily_max"), "d"),
    ("journal", CAP, "{} on Net2 (48)", _age("Net2", "by_truth", "default", "converged_only", "band_coverage_daily_max"),
     ".3f"),
    ("journal", CAP, "0.500 on Net2 ({})", _age("Net2", "by_truth", "default", "converged_only",
                                                "n_junction_seeds_daily_max"), "d"),
    ("journal", CAP, "{} on ky4 (262)", _age("ky4", "by_truth", "default", "converged_only", "band_coverage_daily_max"),
     ".3f"),
    ("journal", CAP, "0.870 on ky4 ({})", _age("ky4", "by_truth", "default", "converged_only",
                                               "n_junction_seeds_daily_max"), "d"),
    ("journal", CAP, "January-to-March samples {} against", "temp_m_julsep_rmse_Net3", ".3f"),
    ("journal", CAP, "against today's {} mg/L on Net3", "temp_b0_julsep_rmse_Net3", ".3f"),
    ("journal", CAP, "held-out readings {} to", "temp_m2_ratio_min", ".3f"),
    ("journal", CAP, "to {} of today's RMSE", "temp_m2_ratio_max", ".3f"),
    ("journal", CAP, "triggered in {} of 4 network", "temp_m2_n_stop", "d"),
    ("journal", CAP, "held-out readings {} of today's RMSE", "toc_readings_ratio_Net3", ".3f"),
    ("journal", CAP, "first-storm map {} against a 0.85 bar", "toc_first_storm_ratio_Net3", ".3f"),
    ("journal", CAP, "in {} of", "audit_n_within", "d"),
    ("journal", CAP, "of {} cells", "audit_n_cells", "d"),
    ("journal", CAP, "oldest tenth of junctions {} mg/L too low", "audit_oldest_bias_Net3", "abs.3f"),
    ("journal", CAP, "recall {} (", "ca_recall_straddle8_Net3", ".3f"),
    ("journal", CAP, "({} of 109", "ca_found_straddle8_Net3", "d"),
    ("journal", CAP, "of {} low junction-days", "ca_low_straddle8_Net3", ".0f"),
    ("journal", CAP, "90% coverage {} (Net3)", "ca_cov90_random_Net3", ".3f"),
    ("journal", CAP, "and {} (Net2)", "ca_cov90_random_Net2", ".3f"),
    ("journal", CAP, "fitted to 15 samples: {} of the posterior", "ca_free_floor_mass_Net3", ".2f"),
    ("journal", CAP, "{} false alarms against", "ca_false_alarms_free_Net3", ".0f"),
    ("journal", CAP, "against the mode's {}", "ca_false_alarms_mode_Net3", ".0f"),
    ("journal", CAP, "too narrow ({} at 8 samples", "ca_cov90_straddle8_Net3", ".3f"),
    # journal, 'Iteration 4 summary' (its key-numbers table and its lists)
    ("journal", SUM, "| {} h (Net3), 5.62", _age("Net3", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", SUM, "(Net3), {} h (Net2)", _age("Net2", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", SUM, "(Net2), {} h (ky4)", _age("ky4", "by_truth", "default", "rmse_h"), ".2f"),
    ("journal", SUM, "| {}, 0.779, 0.811 of junctions", "age_cov_max_Net3", ".3f"),
    ("journal", SUM, "0.791, {}, 0.811 of junctions", "age_cov_max_Net2", ".3f"),
    ("journal", SUM, "0.779, {} of junctions", "age_cov_max_ky4", ".3f"),
    ("journal", SUM, "| {} (Net3, interval", (_S3, ("acceptance", "A3", "bootstrap", "ratio")), ".3f"),
    ("journal", SUM, "interval {} to 1.013", (_S3, ("acceptance", "A3", "bootstrap", "lo90")), ".3f"),
    ("journal", SUM, "0.937 to {}), 1.007", (_S3, ("acceptance", "A3", "bootstrap", "hi90")), ".3f"),
    ("journal", SUM, "), {} (Net2)", (_S2, ("acceptance", "A3", "bootstrap", "ratio")), ".3f"),
    ("journal", SUM, "| {} against 0.101 mg/L", (_S3, ("acceptance", "A4", "rmse_M")), ".3f"),
    ("journal", SUM, "against {} mg/L; 554", (_S3, ("acceptance", "A4", "rmse_B0")), ".3f"),
    ("journal", SUM, "mg/L; {} against 18", (_S3, ("acceptance", "A4", "false_alarms_M_jul_sep")), "d"),
    ("journal", SUM, "554 against {} |", (_S3, ("acceptance", "A4", "false_alarms_B0_jul_sep")), "d"),
    ("journal", SUM, "| {} to 0.916; triggered", "temp_m2_ratio_min", ".3f"),
    ("journal", SUM, "0.774 to {}; triggered", "temp_m2_ratio_max", ".3f"),
    ("journal", SUM, "triggered in {} of 4 cases", "temp_m2_n_stop", "d"),
    ("journal", SUM, "| {} (interval 0.874", (_O3, _T2 + ("ratio",)), ".3f"),
    ("journal", SUM, "(interval {} to 1.016)", (_O3, _T2 + ("bootstrap", "lo90")), ".3f"),
    ("journal", SUM, "0.874 to {}); 0.917", (_O3, _T2 + ("bootstrap", "hi90")), ".3f"),
    ("journal", SUM, "); {} against 0.861", (_O3, _T2 + ("recall_M_TOC",)), ".3f"),
    ("journal", SUM, "against {}; 38", (_O3, _T2 + ("recall_B0",)), ".3f"),
    ("journal", SUM, "; {} against 27", (_O3, _T2 + ("false_alarms_M_TOC",)), "d"),
    ("journal", SUM, "38 against {} |", (_O3, _T2 + ("false_alarms_B0",)), "d"),
    ("journal", SUM, "| {} against 49 of 61", (_O3, ("acceptance", "T3", "robust_triggers", 0, "low_flagged_M_TOC")), "d"),
    ("journal", SUM, "44 against {} of 61", (_O3, ("acceptance", "T3", "robust_triggers", 0, "low_flagged_B0")), "d"),
    ("journal", SUM, "49 of {} |", (_O3, ("acceptance", "T3", "robust_triggers", 0, "n_low")), "d"),
    ("journal", SUM, "| {} (97 of 109)", "ca_recall_straddle8_Net3", ".3f"),
    ("journal", SUM, "0.890 ({} of 109)", "ca_found_straddle8_Net3", "d"),
    ("journal", SUM, "| {} (Net3), 0.887 (Net2)", "ca_cov90_random_Net3", ".3f"),
    ("journal", SUM, "0.896 (Net3), {} (Net2)", "ca_cov90_random_Net2", ".3f"),
    ("journal", SUM, "| {} (Net3), 0.96 (Net2)", "ca_free_floor_mass_Net3", ".2f"),
    ("journal", SUM, "0.60 (Net3), {} (Net2)", (_C2, ("acceptance", "A5_free_grid_reported", "kb_floor_mass_n15")), ".2f"),
    ("journal", SUM, "| {} (interval 0.838", (_A2, _MAT14 + ("ratio",)), ".3f"),
    ("journal", SUM, "(interval {} to 1.463)", (_A2, _MAT14 + ("lo90",)), ".3f"),
    ("journal", SUM, "0.838 to {})", (_A2, _MAT14 + ("hi90",)), ".3f"),
    ("journal", SUM, "| {} / -0.031 mg/L", (_A3, ("acceptance", "age_decile", "bias_2ra")), ".3f"),
    ("journal", SUM, "-0.076 / {} mg/L", (_A3, ("acceptance", "age_decile", "bias_first_order_paired")), ".3f"),
    ("journal", SUM, "| {} / 32 |", (_A3, _FA + ("false_alarms_richer",)), "d"),
    ("journal", SUM, "65 / {} |", (_A3, _FA + ("false_alarms_first_order",)), "d"),
    ("journal", SUM, "| {} of 50 byte-identical", (_BR, _T14 + ("vs_committed", "n_byte_identical")), "d"),
    ("journal", SUM, "of {} byte-identical", (_BR, _T14 + ("vs_committed", "n_files")), "d"),
    ("journal", SUM, "; {} of 92 flagged", (_BR, _T14 + ("app_default_demo", "flagged")), "d"),
    ("journal", SUM, "{} of 36 found", (_BR, _T14 + ("app_default_demo", "found")), "d"),
    ("journal", SUM, "36 of {} found", (_BR, _T14 + ("app_default_demo", "unsampled_violations")), "d"),
    ("journal", SUM, "found, {} false alarm", (_BR, _T14 + ("app_default_demo", "false_alarms")), "d"),
    ("journal", SUM, "off by up to {} in one Net2 scenario", (_WA.format("Net2"), ("by_truth", "default", "loss_split",
                                                                                 "median_abs_diff_max")), ".3f"),
    ("journal", SUM, "a measurable {} times", (_A2, ("acceptance", "materiality", "cells", 11, "ratio")), ".3f"),
    ("journal", SUM, "lower bounds ({} of 92,", "age_lb_max_Net3", "d"),
    ("journal", SUM, "lower bounds (87 of {},", "age_n_junctions_Net3", "d"),
    ("journal", SUM, "{} of 35 and", "age_lb_max_Net2", "d"),
    ("journal", SUM, "29 of {} and", "age_n_junctions_Net2", "d"),
    ("journal", SUM, "and {} of 959 junctions", "age_lb_max_ky4", "d"),
    ("journal", SUM, "799 of {} junctions", "age_n_junctions_ky4", "d"),
    ("journal", SUM, "file's own age on {} of 8 seeds on Net3", _age("Net3", "by_truth", "default", "n_seeds_posterior_worse"),
     "d"),
    ("journal", SUM, "own age on 2 of {} seeds on Net3", _age("Net3", "by_truth", "default", "n_seeds"), "d"),
    ("journal", SUM, "{} of 8 seeds on Net3 and on Net2", _age("Net2", "by_truth", "default", "n_seeds_posterior_worse"), "d"),
    ("journal", SUM, "2 of {} seeds on Net3 and on Net2", _age("Net2", "by_truth", "default", "n_seeds"), "d"),
    ("journal", SUM, "made the mode {} worse in RMSE", "ca_lowph_prior_cost", "s"),
    ("journal", SUM, "kinetics cost {} of RMSE on Net2", "ca_first_order_cost_Net2", "s"),
    ("journal", SUM, "new floor in up to {} of Net2's fits", "lowkb_floor_share_max_Net2", ".0%"),
    ("journal", SUM, "oldest tenth's truth at {} on Net3", "oldest_tenth_cov90_2ra_Net3", ".3f"),
    ("journal", SUM, "biased low by about {} mg/L", (_A3, ("dose_step", "by_truth", "first_order", "1", "bias")), "abs.2f"),
    ("journal", SUM, "biased low by about {} mg/L", (_A3, ("dose_step", "by_truth", "2ra", "1", "bias")), "abs.2f"),
    # README, 'Iteration 4 in one page'
    ("README", ONE, "age at {} (Net3)", "age_cov_max_Net3", ".3f"),
    ("README", ONE, "{} (Net2)", "age_cov_max_Net2", ".3f"),
    ("README", ONE, "{} (ky4) of junctions", "age_cov_max_ky4", ".3f"),
    ("README", ONE, "off by {} per junction on Net3", "wall_share_mae_Net3", ".3f"),
    ("README", ONE, "of the water at {} of 92,", "age_lb_max_Net3", "d"),
    ("README", ONE, "of the water at 87 of {},", "age_n_junctions_Net3", "d"),
    ("README", ONE, "{} of 35 and", "age_lb_max_Net2", "d"),
    ("README", ONE, "29 of {} and", "age_n_junctions_Net2", "d"),
    ("README", ONE, "and {} of 959 junctions", "age_lb_max_ky4", "d"),
    ("README", ONE, "799 of {} junctions", "age_n_junctions_ky4", "d"),
    ("README", ONE, "file's own on {} of 8 seeds on Net3", _age("Net3", "by_truth", "default", "n_seeds_posterior_worse"),
     "d"),
    ("README", ONE, "{} of 8 seeds on Net3 and on Net2", _age("Net2", "by_truth", "default", "n_seeds_posterior_worse"), "d"),
    ("README", ONE, "off by up to {} in one Net2 scenario", (_WA.format("Net2"), ("by_truth", "default", "loss_split",
                                                                                "median_abs_diff_max")), ".3f"),
    ("README", ONE, "made the mode {} worse in RMSE", "ca_lowph_prior_cost", "s"),
    ("README", ONE, "kinetics cost {} of RMSE on Net2", "ca_first_order_cost_Net2", "s"),
    ("README", ONE, "truth at {} there", "oldest_tenth_cov90_2ra_Net3", ".3f"),
    ("README", ONE, "a measurable {} times", (_A2, ("acceptance", "materiality", "cells", 11, "ratio")), ".3f"),
    ("README", ONE, "new floor in up to {} of Net2's fits", "lowkb_floor_share_max_Net2", ".0%"),
    ("README", ONE, "biased low by about {} mg/L", (_A3, ("dose_step", "by_truth", "first_order", "1", "bias")), "abs.2f"),
    ("README", ONE, "biased low by about {} mg/L", (_A3, ("dose_step", "by_truth", "2ra", "1", "bias")), "abs.2f"),
    ("README", ONE, "January-to-March samples {} against", "temp_m_julsep_rmse_Net3", ".3f"),
    ("README", ONE, "against today's {} mg/L on Net3", "temp_b0_julsep_rmse_Net3", ".3f"),
    ("README", ONE, "held-out readings {} to", "temp_m2_ratio_min", ".3f"),
    ("README", ONE, "to {} of today's RMSE", "temp_m2_ratio_max", ".3f"),
    ("README", ONE, "triggered in {} of 4 cases", "temp_m2_n_stop", "d"),
    ("README", ONE, "map RMSE {} of today's against a 0.85 bar", "toc_first_storm_ratio_Net3", ".3f"),
    ("README", ONE, "within 10% in {} of", "audit_n_within", "d"),
    ("README", ONE, "of {} cells", "audit_n_cells", "d"),
    ("README", ONE, "oldest tenth of junctions {} mg/L too low", "audit_oldest_bias_Net3", "abs.3f"),
    ("README", ONE, "recall {} (", "ca_recall_straddle8_Net3", ".3f"),
    ("README", ONE, "({} of 109", "ca_found_straddle8_Net3", "d"),
    ("README", ONE, "of {} low junction-days", "ca_low_straddle8_Net3", ".0f"),
    ("README", ONE, "90% coverage {} (Net3)", "ca_cov90_random_Net3", ".3f"),
    ("README", ONE, "and {} (Net2)", "ca_cov90_random_Net2", ".3f"),
    ("README", ONE, "put {} of their posterior", "ca_free_floor_mass_Net3", ".2f"),
    ("README", ONE, "too narrow under the straddle rule ({} at 8 samples", "ca_cov90_straddle8_Net3", ".3f"),
    ("README", ONE, "write {} of 50 files", (_BR, _T14 + ("vs_committed", "n_byte_identical")), "d"),
    ("README", ONE, "of {} files byte-identical", (_BR, _T14 + ("vs_committed", "n_files")), "d"),
    ("README", ONE, "flags {} of 92", (_BR, _T14 + ("app_default_demo", "flagged")), "d"),
    ("README", ONE, "of {} junctions", (_BR, _T14 + ("app_default_demo", "junctions")), "d"),
    ("README", ONE, "finds {} of 36", (_BR, _T14 + ("app_default_demo", "found")), "d"),
    ("README", ONE, "with {} false alarm", (_BR, _T14 + ("app_default_demo", "false_alarms")), "d"),
]


@dataclass
class Section:
    doc: str
    title: str
    text: str
    cited: list[str] = field(default_factory=list)


def tracked_outputs() -> list[str]:
    out = subprocess.run(["git", "-C", REPO, "ls-files", "outputs"], capture_output=True, text=True).stdout.split()
    return [f for f in out if f.endswith((".json", ".csv"))]


_VALUES: dict[str, np.ndarray] = {}


def file_values(rel: str) -> np.ndarray:
    """Sorted unique absolute values of every number in a committed JSON or CSV file (numbers inside strings too)."""
    if rel in _VALUES:
        return _VALUES[rel]
    vals: list[float] = []

    def from_str(s: str):
        for m in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", s):
            try:
                vals.append(float(m))
            except ValueError:
                pass

    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                from_str(str(k))
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, bool) or x is None:
            return
        elif isinstance(x, (int, float)):
            if not (isinstance(x, float) and math.isnan(x)):
                vals.append(float(x))
        elif isinstance(x, str):
            from_str(x)

    path = os.path.join(REPO, rel)
    if rel.endswith(".json"):
        with open(path) as fh:
            x = json.load(fh)
        if rel == SELF_REPORT and isinstance(x, dict):
            # leave out this scan's own entry in the checks report, which would otherwise feed its own counts back
            # into the pool and change the committed report on every run
            x = {**x, "checks": [c for c in x.get("checks", []) if c.get("name") != "numbers_trace_to_outputs"]}
        walk(x)
    else:
        df = pd.read_csv(path, float_precision="round_trip")
        for c in df.columns:
            from_str(str(c))
            num = pd.to_numeric(df[c], errors="coerce")
            vals.extend(num.dropna().astype(float).tolist())
            if df[c].dtype == object:
                for s in df[c].dropna().astype(str).unique():
                    from_str(s)
    a = np.unique(np.abs(np.asarray(vals, dtype=float)))
    _VALUES[rel] = a[np.isfinite(a)]
    return _VALUES[rel]


def pool(files) -> np.ndarray:
    arrs = [file_values(f) for f in files]
    return np.unique(np.concatenate(arrs)) if arrs else np.zeros(0)


def code_values() -> np.ndarray:
    """Every numeric literal written in residualmap/*.py and app.py (code, comments and docstrings): the settings, the
    plan's bars and the literature values the code cites."""
    vals = []
    files = [os.path.join(REPO, "app.py")] + [os.path.join(REPO, "residualmap", f)
                                              for f in sorted(os.listdir(os.path.join(REPO, "residualmap")))
                                              if f.endswith(".py") and f != "numbers.py"]
    for p in files:
        with open(p) as fh:
            for m in re.findall(r"(?<![\w.])\d+(?:_\d{3})*(?:\.\d+)?(?:e-?\d+)?", fh.read()):
                try:
                    vals.append(float(m.replace("_", "")))
                except ValueError:
                    pass
    return np.unique(np.abs(np.asarray(vals)))


def parse_token(tok: str) -> tuple[float, float, bool]:
    """(value, half-width of its printed precision, is a percent)."""
    pct = tok.endswith("%")
    s = tok.rstrip("%").replace(",", "").lstrip("+")
    if "e" in s:
        mant, exp = s.split("e")
        d = len(mant.split(".")[1]) if "." in mant else 0
        half = 0.5 * 10 ** (int(exp) - d)
    else:
        d = len(s.split(".")[1]) if "." in s else 0
        half = 0.5 * 10 ** (-d)
    return abs(float(s)), half, pct


def found(values: np.ndarray, v: float, half: float, pct: bool) -> bool:
    if len(values) == 0:
        return False
    tol = 1e-9 * max(1.0, v)
    cands = [(v - half - tol, v + half + tol)]
    if pct:
        cands.append(((v - half) / 100 - tol, (v + half) / 100 + tol))
    for lo, hi in cands:
        i = np.searchsorted(values, lo)
        if i < len(values) and values[i] <= hi:
            return True
    return False


# ----------------------------------------------------------------------------- the sections scanned
def _cited(text: str, tracked: list[str]) -> list[str]:
    names = set(re.findall(r"([\w\-<>*]+\.(?:json|csv))", text))
    out = set()
    for n in names:
        rx = re.compile("^" + re.escape(n).replace(re.escape("<net>"), r"[A-Za-z0-9]+")
                        .replace(re.escape("<n>"), r"\d+").replace(re.escape("*"), r".*") + "$")
        out.update(f for f in tracked if rx.match(os.path.basename(f)))
    return sorted(out)


def sections() -> list[Section]:
    tracked = tracked_outputs()
    out = []
    with open(JOURNAL) as fh:
        j = fh.read()
    part = j[j.index(JOURNAL_START):]
    for s in re.split(r"(?m)^(?=##+ )", part):
        if s.strip():
            out.append(Section("journal", s.splitlines()[0].lstrip("# ").strip(), s, _cited(s, tracked)))
    with open(README) as fh:
        r = fh.read()
    for s in re.split(r"(?m)^(?=## )", r):
        head = s.splitlines()[0] if s.strip() else ""
        title = head.lstrip("# ").strip()
        if "iteration 4" in head.lower() or title.startswith("The app"):
            out.append(Section("README", title, s, _cited(s, tracked)))
        elif title == "Honest limitations":
            bullets = [b for b in re.split(r"(?m)^(?=- )", s) if b.startswith("- ") and IT4_TASK.search(b)]
            text = "\n".join(bullets)
            out.append(Section("README", "Honest limitations (iteration-4 bullets)", text, _cited(text, tracked)))
    return out


# ----------------------------------------------------------------------------- explicit allowlist
# (section title prefix or '' for any, the printed token, the reason).  Every entry is a number that is not in any
# output file and is not a time, date, citation, identifier, run cost or a value written in the code.
ALLOW: list[tuple[str, str, str]] = [
    ("Task 10, temperature", "12,100", "literature value: Blokker et al. 2014's E/R, about 12,100 K (cited with its "
                                       "source in the pilot protocol, section 5)"),
]


@dataclass
class Result:
    traced: dict = field(default_factory=dict)        # level -> count
    classified: dict = field(default_factory=dict)    # class -> list of (section, token, context)
    unexplained: list = field(default_factory=list)   # (section, token, context)
    registry: list = field(default_factory=list)      # (ok, entry, message)
    sources_ok: bool = True
    sources_error: str | None = None
    n_tokens: int = 0
    # [caught, tested] per kind of traced number: would a one-unit change in its last printed digit be caught by the
    # pool that traced it?  Integers (counts, 'X of Y') are matched at plus or minus 0.5, so the scan has almost no power
    # on them; the registry is the exact check for the numbers it covers.
    power: dict = field(default_factory=lambda: {"decimal": [0, 0], "integer": [0, 0]})


def scan(verbose: bool = False) -> Result:
    res = Result()
    tracked = tracked_outputs()
    it4 = [f for f in tracked if f.startswith(IT4_DIRS)]
    p_it4, p_all, p_code = pool(it4), pool(tracked), code_values()
    for sec in sections():
        p_sec = pool(sec.cited)
        # paragraph-level pools: the files a paragraph cites (a table also takes its caption paragraph's), tried first
        blocks, pos = [], 0
        for b in re.split(r"(\n\s*\n)", sec.text):
            if b.strip():
                blocks.append([pos, pos + len(b), _cited(b, tracked), b.lstrip().startswith("|")])
            pos += len(b)
        for i, b in enumerate(blocks):
            if b[3] and i > 0:
                b[2] = sorted(set(b[2]) | set(blocks[i - 1][2]))
        p_par = {i: pool(b[2]) for i, b in enumerate(blocks)}
        text = sec.text
        for cls, rx in BLANKS:
            for m in rx.finditer(text):
                res.classified.setdefault(cls, []).append((sec.title, m.group(0)[:40], ""))
                a, b = m.start(), m.end()
                if cls not in ("code or path", "link target") and (
                        (a > 0 and (text[a - 1].isdigit() or (text[a - 1] == "." and a > 1 and text[a - 2].isdigit())))
                        or re.match(r"\.?\d", text[b:b + 2])):
                    # a blanked span that takes part of a number: that number would never be read
                    res.unexplained.append((f"{sec.doc}: {sec.title[:60]}", m.group(0)[:40],
                                            f"blanked as {cls} but touches a digit: " + text[max(0, a - 30):b + 30]))
            text = rx.sub(lambda m: " " * len(m.group(0)), text)
        text = re.sub(r"(?<![\w])[x×](?=\d)", " ", text)          # x0.85, ×1.05: the multiplier's number
        ofs = {m.start(1): (int(m.group(1)), int(m.group(2))) for m in OF.finditer(text)}
        for m in NUM.finditer(text):
            tok = m.group(0)
            res.n_tokens += 1
            v, half, pct = parse_token(tok)
            bi = next((i for i, b in enumerate(blocks) if b[0] <= m.start() < b[1]), None)
            hit = None
            for level, p in (("file cited in its paragraph", p_par.get(bi, np.zeros(0))),
                             ("file cited in its section", p_sec), ("iteration-4 outputs", p_it4),
                             ("outputs", p_all)):
                if found(p, v, half, pct):
                    hit = (level, p)
                    break
            if hit:
                res.traced[hit[0]] = res.traced.get(hit[0], 0) + 1
                # power: would a one-unit error in the last printed digit be caught by the pool that traced it?
                if "e" not in tok and not pct:
                    kind = "decimal" if "." in tok else "integer"
                    for v2 in (v + 2 * half, v - 2 * half):
                        if v2 > 0:
                            res.power[kind][1] += 1
                            res.power[kind][0] += not found(hit[1], v2, half, False)
                continue
            ctx = text[max(0, m.start() - 50):m.end() + 30].replace("\n", " ")
            cls = None
            if COST.match(text[m.end():m.end() + 12]):
                cls = "run cost (logs, df)"
            elif m.start() in ofs:
                x, y = ofs[m.start()]
                if y and any(found(p, x / y, 0.0005, False) for p in (p_sec, p_it4)):
                    cls = "count derived from a stored rate"
            if cls is None and found(p_code, v, half, False):
                cls = "setting or literature value written in the code"
            if cls is None:
                for title, t, why in ALLOW:
                    if t == tok and sec.title.startswith(title):
                        cls = f"allowlist: {why}"
                        break
            if cls is None:
                res.unexplained.append((f"{sec.doc}: {sec.title[:60]}", tok, ctx))
            else:
                res.classified.setdefault(cls, []).append((sec.title[:60], tok, ctx))
    return res


def check_registry(res: Result) -> None:
    from .capabilities import numbers as cap_numbers
    try:
        cap = cap_numbers()
    except Exception as e:  # noqa: BLE001
        res.sources_ok, res.sources_error = False, f"{type(e).__name__}: {e}"
        cap = {}
    secs = {(s.doc, s.title): s.text for s in sections()}
    for entry in _registry():
        doc, head, words, source, fmt = entry
        if "{}" not in words:
            words = "{}"
        text = next((t for (d, ti), t in secs.items() if d == doc and ti.startswith(head)), None)
        try:
            if isinstance(source, str) and source in DERIVED:
                val = DERIVED[source][1]()
            elif isinstance(source, str):
                val = cap[source]
            else:
                rel, keys = source
                with open(os.path.join(REPO, rel)) as fh:
                    val = json.load(fh)
                for k in keys:
                    val = val[k]
            if fmt.startswith("abs"):
                val, fmt = abs(val), fmt[3:]
            got = format(val, fmt) if fmt != "s" else str(val)
        except Exception as e:  # noqa: BLE001
            res.registry.append((False, entry, f"source does not resolve: {type(e).__name__}: {e}"))
            continue
        want = words.format(got)
        if text is None:
            res.registry.append((False, entry, f"section {head!r} not found in {doc}"))
        elif not _contains(_flat(text), _flat(want)):
            res.registry.append((False, entry, f"outputs give {got}, and {want!r} is not in section {head!r}"))
        else:
            res.registry.append((True, entry, "ok"))


def _flat(t: str) -> str:
    return " ".join(t.split())


def _contains(text: str, want: str) -> bool:
    """want in text, not as part of a longer number (so '0.79' does not match inside '0.791')."""
    for m in re.finditer(re.escape(want), text):
        a, b = m.start(), m.end()
        if (a > 0 and want[0].isdigit() and (text[a - 1].isdigit() or text[a - 1] == ".")) or \
           (b < len(text) and want[-1].isdigit() and (text[b].isdigit() or (text[b] == "." and b + 1 < len(text)
                                                                            and text[b + 1].isdigit()))):
            continue
        return True
    return False


def _registry() -> list[tuple]:
    return REGISTRY


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="print every number not traced to a file, with its class")
    a = ap.parse_args(argv)
    res = scan()
    check_registry(res)
    n_traced = sum(res.traced.values())
    n_class = sum(len(v) for k, v in res.classified.items() if k not in ("code or path", "link target"))
    print(f"iteration-4 numbers read: {res.n_tokens}; traced to a committed output file: {n_traced} "
          f"({', '.join(f'{k} {v}' for k, v in res.traced.items())})")
    print(f"classified, not results ({n_class} spans and numbers):")
    for k, v in sorted(res.classified.items()):
        if k in ("code or path", "link target"):
            continue
        print(f"  {k}: {len(v)}")
        if a.list:
            for sec, tok, ctx in v:
                print(f"      [{sec[:40]}] {tok!r}  ...{ctx.strip()[:110]}...")
    for kind, (caught, tested) in res.power.items():
        if tested:
            print(f"power of the scan on {kind}s: a one-unit change in the last printed digit of a traced {kind} would "
                  f"be caught {caught} of {tested} times ({caught / tested:.0%})")
    print("  (an integer is matched at plus or minus 0.5, and an 'X of Y' count found in no file is accepted when X/Y "
          "is a stored rate; the registry below is the exact check for the numbers it covers)")
    print(f"registry: {sum(ok for ok, _, _ in res.registry)} of {len(res.registry)} task-14 headline numbers match "
          f"their file and key")
    for ok, entry, msg in res.registry:
        if not ok:
            print(f"  MISMATCH {entry[:3]}: {msg}")
    print(f"the app's table: every source in capabilities.SOURCES resolves: {res.sources_ok}"
          + (f" ({res.sources_error})" if res.sources_error else ""))
    print(f"unexplained: {len(res.unexplained)}")
    for sec, tok, ctx in res.unexplained:
        print(f"  [{sec}] {tok!r}  ...{ctx.strip()[:120]}...")
    ok = not res.unexplained and res.sources_ok and all(ok for ok, _, _ in res.registry)
    print("OK: every iteration-4 number traces to outputs/ or is classified" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
