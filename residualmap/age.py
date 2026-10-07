"""
age.py: water age made visible, and where chlorine is lost (iteration 4, journal task 9).

Water age has always been inside the model: every EPANET run carries water through every pipe and tank hour by
hour, and the nominal age at each hour is the first input of the discrepancy GP.  What the operator never saw is
the age itself, how far the file's own demand and roughness errors could move it, and how much of each junction's
chlorine loss happens in the water and how much at the pipe walls.  This module computes those three things.

  hydraulic_age_band(sc)  EPANET AGE runs of the operator's nominal model at the 9 hydraulic members of the grid
                          (simgp.DEMAND_GRID x simgp.ROUGH_GRID, the committed axes, in the grid's own order): the
                          min, median and max age per junction and hour.  The (1, 1) member is the nominal age
                          itself, bit for bit (a saved check).
  posterior_age(model)    the 9 ages weighted by a fitted SimGP24's posterior over those members (its weights
                          summed over the decay and dose axes): what the grab samples say about the hydraulics,
                          turned into hours.  Posterior age = sum_k w_k a_jh(dm_k, rm_k).
  loss_split(model)       EPANET runs of the MAP member with wall decay off, with bulk decay off and with both off,
                          giving each junction's chlorine loss in the water and at the pipe walls.
  truth_loss_split(sc)    the same split on the hidden truth (build_scenario(truth_loss_split=True)), so the
                          model's split can be scored.

The loss split (per junction, each term the mean over the 24 hours of the last day):
    L_tot  = ln(C_ref / C)              C     : the member as calibrated
    L_bulk = ln(C_ref / C at kw = 0)    C_ref : the same run with no decay at all; it equals the source dose
    L_wall = ln(C_ref / C at kb = 0)            wherever source water has replaced the file's initial water,
    wall share = 1 - L_bulk / L_tot             and it also handles sources at different doses (the truth)
    bulk share = L_bulk / L_tot
On a single plug-flow path the split is exact: EPANET's first-order bulk term and its mass-transfer-limited
first-order wall term are both linear in C (Rossman, Clark & Grayman 1994, J Environ Eng 120(4):803), so ln C falls
linearly with age at rate kb + kw_eff and L_bulk + L_wall = L_tot.  Where flows of different ages mix it is
approximate, and (L_bulk + L_wall) / L_tot is reported as a non-additivity diagnostic.  For context, wall loss is
reported at up to 97% of the total in old pipes, and bulk loss at up to 35% in newer pipes (Maleki et al. 2023,
Water Supply 23(2):657).

Every age here is a 7-day-run age.  The run starts with the file's water in every pipe, junction and tank, at
the file's initial age (0 in the example files), so wherever some of that water is still arriving on the last
day, the age EPANET reports is lower than the real age of the water: a lower bound.  simulate.initial_water_share
measures that share per junction and hour (one more AGE run with the starting ages raised); a junction holding
more than INITIAL_SHARE_MAX of it is reported as a lower bound, in the outputs and in the app.  The plan's
160 h cut (AGE_UNCONVERGED_H, the true daily-max age) only catches the extreme cases: ages on the last day cannot
exceed 168 h plus the file's initial age.
"""
from __future__ import annotations

import itertools
import os
import tempfile
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .chemistry import remove_epanet_files
from .simgp import DEMAND_GRID, HOURS, ROUGH_GRID, check_grid_order, n_hydraulic
from .simulate import _last_day, _run_quality, initial_water_share, load, simulate_nominal_chlorine

AGE_UNCONVERGED_H = 160.0      # the plan's exclusion: true daily-max age above this on the last day of the 7-day run
INITIAL_SHARE_MAX = 0.05       # above this share of the run's starting water, a junction's age is a lower bound
BAND_COVERAGE_BAR = 0.85       # the 9-member range is called a 90% band only if it holds the truth this often
RANGE_LABEL = "range across demand and roughness errors"
MIN_LOSS = 0.02                # ln units: below a 2% loss the shares are undefined (junctions next to a source)
C_FLOOR = 1e-6                 # concentrations are floored here only to keep the logarithms finite


def band_label(coverage: float | None) -> str:
    """'90% band' only when the measured coverage of the 9-member range is at least BAND_COVERAGE_BAR; otherwise the
    honest name of what it is: the range across the file's demand and roughness errors."""
    return "90% band" if coverage is not None and coverage >= BAND_COVERAGE_BAR else RANGE_LABEL


# ----------------------------------------------------------------------------- the 9 hydraulic members
def hydraulic_members() -> list[tuple[float, float]]:
    """(demand multiplier, roughness multiplier) in the order of the grid's hydraulic block (itertools.product)."""
    return [(float(dm), float(rm)) for dm, rm in itertools.product(DEMAND_GRID, ROUGH_GRID)]


def simulate_nominal_age(name: str, demand_mult: float = 1.0, rough_mult: float = 1.0, with_initial_share: bool = False):
    """Water age (hours, hour x junction, last day) on the operator's NOMINAL model with the global demand and
    Hazen-Williams C scaled as in simulate_nominal_chlorine.  The run happens in a temporary directory.  At (1, 1) it
    writes the same .inp as nominal_scenario's AGE run, so the result is the nominal age bit for bit.
    with_initial_share=True returns (age, initial_water_share), the second from one more run."""
    wn = load(name)
    wn.options.quality.parameter = "AGE"
    for _, pipe in wn.pipes():
        pipe.roughness = pipe.roughness * rough_mult
    if demand_mult != 1.0:
        for _, j in wn.junctions():
            for ts in j.demand_timeseries_list:
                ts.base_value = ts.base_value * demand_mult
    age = _last_day(_run_quality(wn), wn.junction_name_list) / 3600.0
    return (age, initial_water_share(wn, age)) if with_initial_share else age


@dataclass
class AgeBand:
    """Water age of the 9 hydraulic members: ages[k, hour, junction] in hours, members[k] = (dm, rm), and the
    nominal member's initial-water share (hour x junction): where it is above INITIAL_SHARE_MAX the age is a
    lower bound."""
    members: list[tuple[float, float]]
    ages: np.ndarray
    junctions: list[str]
    initial_share: pd.DataFrame | None = None

    def _frame(self, a: np.ndarray) -> pd.DataFrame:
        return pd.DataFrame(a, index=HOURS, columns=self.junctions)

    @property
    def nominal_index(self) -> int:
        return self.members.index((1.0, 1.0))

    @property
    def nominal(self) -> pd.DataFrame:
        return self._frame(self.ages[self.nominal_index])

    def lo(self) -> pd.DataFrame:
        return self._frame(self.ages.min(axis=0))

    def med(self) -> pd.DataFrame:
        return self._frame(np.median(self.ages, axis=0))

    def hi(self) -> pd.DataFrame:
        return self._frame(self.ages.max(axis=0))

    def daily(self, stat: str = "mean") -> pd.DataFrame:
        """Each member's daily mean (or max) age per junction: members x junction."""
        a = self.ages.mean(axis=1) if stat == "mean" else self.ages.max(axis=1)
        return pd.DataFrame(a, index=pd.MultiIndex.from_tuples(self.members, names=["demand", "rough"]),
                            columns=self.junctions)

    def daily_range(self, stat: str = "mean") -> tuple[pd.Series, pd.Series]:
        """[min, max] over the 9 members of the daily mean (or max) age, per junction."""
        d = self.daily(stat)
        return d.min(), d.max()

    def share(self, stat: str = "mean") -> pd.Series:
        """The nominal model's initial-water share per junction: its daily mean, or its largest hour ('max')."""
        return self.initial_share.mean() if stat == "mean" else self.initial_share.max()

    def lower_bound(self, stat: str = "mean") -> pd.Series:
        """True where the nominal daily-mean (or daily-max) age is a lower bound: more than INITIAL_SHARE_MAX of the
        junction's water (on the day's average, or at any hour for 'max') is still the run's starting water."""
        return self.share(stat) > INITIAL_SHARE_MAX


def hydraulic_age_band(sc) -> AgeBand:
    """EPANET AGE for the 9 hydraulic members on the nominal model of scenario sc (0.03 s a run on Net3 and Net2,
    a few seconds on ky4), plus one run for the nominal member's initial-water share."""
    members = hydraulic_members()
    ages, share = [], None
    for dm, rm in members:
        if (dm, rm) == (1.0, 1.0):
            a, sh = simulate_nominal_age(sc.wn_name, dm, rm, with_initial_share=True)
            share = sh.loc[HOURS, sc.junctions]
        else:
            a = simulate_nominal_age(sc.wn_name, dm, rm)
        ages.append(a.loc[HOURS, sc.junctions].values)
    return AgeBand(members, np.stack(ages).astype(float), list(sc.junctions), share)


# ----------------------------------------------------------------------------- posterior age
def hydraulic_weights(model) -> np.ndarray:
    """A fitted SimGP24's posterior over the 9 hydraulic members: its joint (member x dose) weights summed over the
    dose axis and over the decay triples.  Raises unless the model's grid has the full hydraulic block."""
    n_hyd = n_hydraulic("full")
    if len(model.params) % n_hyd:
        raise ValueError("posterior_age needs a model on the 'full' grid (demand and roughness axes)")
    check_grid_order(model.params, n_hyd)
    if [tuple(map(float, p[3:])) for p in model.params[:n_hyd]] != hydraulic_members():
        raise ValueError("the model's hydraulic members are not DEMAND_GRID x ROUGH_GRID in grid order")
    w = np.asarray(model.W_, dtype=float)
    w = w.sum(axis=1) if w.ndim == 2 else w
    return w.reshape(-1, n_hyd).sum(axis=0)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> np.ndarray:
    """Quantile q of the discrete distribution putting weights[k] on values[k, j], per column j: the smallest value
    whose cumulative weight reaches q."""
    order = np.argsort(values, axis=0, kind="stable")
    v = np.take_along_axis(values, order, axis=0)
    cw = np.cumsum(weights[order], axis=0) / weights.sum()
    idx = np.argmax(cw >= q - 1e-12, axis=0)
    return v[idx, np.arange(values.shape[1])]


@dataclass
class PosteriorAge:
    weights: np.ndarray            # (9,) posterior over the hydraulic members, in AgeBand.members order
    by_hour: pd.DataFrame          # hour x junction, posterior-weighted age (h)
    daily_mean: pd.Series          # per junction, posterior-weighted daily-mean age (h)
    lo90: pd.Series                # central 90% of the posterior over the 9 members' daily-mean ages
    hi90: pd.Series


def posterior_age(model, band: AgeBand | None = None) -> PosteriorAge:
    """Posterior age = sum_k w_k a_jh(dm_k, rm_k), with w the model's weights summed over the decay and dose axes.
    The band is the central 90% of that discrete posterior over the 9 members' daily-mean ages (5% and 95% weighted
    quantiles), so it narrows as the samples pin the hydraulics down."""
    band = band if band is not None else hydraulic_age_band(model.sc)
    w = hydraulic_weights(model)
    by_hour = pd.DataFrame(np.tensordot(w, band.ages, axes=1), index=HOURS, columns=band.junctions)
    dm = band.ages.mean(axis=1)                                   # 9 x J daily-mean ages
    lo = pd.Series(weighted_quantile(dm, w, 0.05), index=band.junctions)
    hi = pd.Series(weighted_quantile(dm, w, 0.95), index=band.junctions)
    return PosteriorAge(w, by_hour, pd.Series(w @ dm, index=band.junctions), lo, hi)


# ----------------------------------------------------------------------------- where chlorine is lost
SPLIT_RUNS = ("full", "bulk_only", "wall_only", "no_decay")


def split_from_runs(runs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per junction: L_tot, L_bulk, L_wall (mean over the hours of ln(C_ref / C)), wall_share, bulk_share and
    nonadditivity = (L_bulk + L_wall) / L_tot.  runs maps SPLIT_RUNS to hour x junction chlorine of one network
    under the same hydraulics: as calibrated, with kw = 0, with kb = 0, and with no decay (C_ref).  The shares are
    NaN where L_tot < MIN_LOSS (no material loss to split)."""
    ref = np.clip(runs["no_decay"].values, C_FLOOR, None)
    L = {k: np.log(ref / np.clip(runs[k].values, C_FLOOR, None)).mean(axis=0) for k in ("full", "bulk_only", "wall_only")}
    out = pd.DataFrame({"L_tot": L["full"], "L_bulk": L["bulk_only"], "L_wall": L["wall_only"]},
                       index=runs["full"].columns)
    ok = out.L_tot >= MIN_LOSS
    out["bulk_share"] = (out.L_bulk / out.L_tot).where(ok)
    out["wall_share"] = 1.0 - out["bulk_share"]
    out["nonadditivity"] = ((out.L_bulk + out.L_wall) / out.L_tot).where(ok)
    return out


def loss_split(model) -> pd.DataFrame:
    """The split for a fitted SimGP24's MAP member (kb, kw, gamma, demand, roughness, dose): four EPANET runs of the
    nominal model, as calibrated, with kw = 0, with kb = 0 and with no decay, each in a temporary directory.  Also
    returns source_water_share: the no-decay run over the dose, i.e. the share of each junction's water (mean over the
    last day) that entered from a source during the 7-day run.  Raises when the model has no samples (fit_prior has
    no MAP member)."""
    if getattr(model, "map_params_", None) is None:
        raise ValueError("the loss split needs a calibrated member: fit the model on at least one sample")
    cond = getattr(model, "cond", None)
    if cond is not None and not cond.is_default and cond.sim_kwargs() != {"kb_scale": 1.0, "kw_scale": 1.0, "temp_C": None}:
        # these runs use the member's rates as they are: right for today's grids and for task 12's chloramine grid
        # (first order at 20 C, the neutral keywords), wrong for a condition that rescales them (a temperature or TOC)
        raise NotImplementedError("loss_split runs the 20 C member only; a model fitted under a chemistry "
                                  "condition needs that condition's rates here")
    sc = model.sc
    kb, kw, g, dm, rm = model.map_params_
    dose = float(sc.source_dose) * float(model.map_dose_)
    rates = {"full": (kb, kw), "bulk_only": (kb, 0.0), "wall_only": (0.0, kw), "no_decay": (0.0, 0.0)}
    runs = {}
    with tempfile.TemporaryDirectory(prefix="rm_split_") as tmp:
        for k, (b, w) in rates.items():
            prefix = os.path.join(tmp, k)
            try:
                runs[k] = simulate_nominal_chlorine(sc.wn_name, b, w, g, dose, dm, rm, file_prefix=prefix).loc[HOURS, sc.junctions]
            finally:
                remove_epanet_files(prefix)
    out = split_from_runs(runs)
    # every source carries the same dose here, so the no-decay run over the dose is the share of the water at the
    # junction that entered from a source during the run; the rest is still the file's initial contents (quality 0
    # in the example files), a start-up transient the truth shares
    out["source_water_share"] = (runs["no_decay"].mean() / dose).values
    out.attrs["member"] = {"kb": kb, "kw": kw, "gamma": g, "demand": dm, "rough": rm, "dose_mult": model.map_dose_}
    return out


def truth_loss_split(sc) -> pd.DataFrame:
    """The same split on the hidden truth: needs build_scenario(..., truth_loss_split=True)."""
    if getattr(sc, "truth_loss_runs", None) is None:
        raise ValueError("build the scenario with truth_loss_split=True")
    return split_from_runs({k: sc.truth_loss_runs[k].loc[HOURS, sc.junctions] for k in SPLIT_RUNS})


def oldest_water(age_by_hour: pd.DataFrame) -> tuple[str, float]:
    """(junction, hours) of the largest daily-maximum age."""
    m = age_by_hour.max()
    j = str(m.idxmax())
    return j, float(m[j])
