"""
app.py — ResidualMap for an operator.

    streamlit run app.py

Upload the EPANET model you already have, type in the chlorine grab samples you already take, and get:
a chlorine map of every junction with a 90% band, the probability each junction is below the 0.2 mg/L
minimum residual (at any hour, or at its daily minimum), the best next places to sample with a reason
each, and a one-page PDF.  No accounts, no database; everything runs on this laptop.

Demo mode simulates a hidden "true" network from the same file (the experiment's scenario generator)
and draws daytime samples from it, so the app can be shown without real data — and the truth can be
revealed to see how the map did.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import textwrap
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import wntr
from scipy.stats import norm

from residualmap.age import INITIAL_SHARE_MAX, RANGE_LABEL, hydraulic_age_band, loss_split, oldest_water
from residualmap.features import build_features
from residualmap.route import plan_route
from residualmap.simgp import DAY_HOURS, GRIDS, SimGP24, simulator_grid_24h
from residualmap.simulate import LIB, build_scenario, nominal_scenario, source_nodes

warnings.filterwarnings("ignore")
UPLOAD_DIR = "outputs/app_uploads"
AGE_RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "chem", "water_age_{}.json")
CACHE_DIR = "outputs/cache"
PDF_WHY_CHARS = 84     # characters per line of the PDF route table's 'why' column
EXAMPLES = {"Net3 (EPANET example, 92 junctions)": "Net3", "Net2 (EPANET example, tank-fed, 35 junctions)": "Net2",
            "ky4 (KYPIPE dataset, 959 junctions; first run takes about 15 min)": "ky4"}

st.set_page_config(page_title="ResidualMap", layout="wide")
st.title("ResidualMap")
st.caption("Your EPANET model + the grab samples you already take → a chlorine map of every junction, "
           "with an honest band, the junctions likely below 0.2 mg/L at any hour of the day, and where to sample next.")


# ----------------------------------------------------------------------------- inputs
with st.sidebar:
    st.header("1. Your network")
    src = st.radio("Model", ["Example network", "Upload my .inp"], horizontal=True)
    if src == "Upload my .inp":
        up = st.file_uploader("EPANET .inp file", type=["inp"])
        if up is None:
            st.stop()
        os.makedirs(UPLOAD_DIR, exist_ok=True)
        digest = hashlib.sha1(up.getvalue()).hexdigest()[:8]
        net_path = os.path.join(UPLOAD_DIR, f"{digest}_{os.path.basename(up.name)}")
        if not os.path.exists(net_path):
            with open(net_path, "wb") as fh:
                fh.write(up.getvalue())
    else:
        net_path = EXAMPLES[st.selectbox("Example", list(EXAMPLES))]
    dose = st.number_input("Chlorine dose leaving the plant (mg/L)", 0.2, 4.0, 1.2, 0.1,
                           help="Free chlorine at the source. The model treats it as ±10% uncertain.")
    threshold = st.number_input("Minimum residual (mg/L)", 0.05, 1.0, 0.2, 0.05)

    st.header("2. Your samples")
    demo = st.toggle("Demo: simulate a hidden truth and draw samples from it", value=True,
                     help="Off = type your own grab samples in the table.")
    if demo:
        n_demo = st.slider("Demo samples (daytime, random junctions)", 0, 15, 8)
        demo_seed = st.number_input("Demo scenario", 0, 99, 0)
        demo_kb = 0.40 if os.path.basename(net_path) != "Net2" else 0.10
        demo_kw = 0.70 if os.path.basename(net_path) != "Net2" else 0.20

    st.header("3. Route")
    K = st.slider("Sites on next month's route", 3, 12, 6)


# ----------------------------------------------------------------------------- heavy lifting, cached
@st.cache_resource(show_spinner=False)
def prepare(path: str, dose: float):
    """Nominal model, physics features and the 675-run simulator grid (cached on disk too)."""
    sc = nominal_scenario(path, sample_hour=14, source_dose=dose)
    X = build_features(sc)
    simulator_grid_24h(sc, CACHE_DIR, "full")
    return sc, X


@st.cache_resource(show_spinner=False)
def water_age(path: str, _sc):
    """Water age of the operator's model at the grid's 9 demand x roughness settings, and the nominal model's share
    of the simulation's starting water (10 EPANET AGE runs, once per network)."""
    return hydraulic_age_band(_sc)


@st.cache_resource(show_spinner=False)
def chlorine_loss(path: str, dose: float, member: tuple, dose_mult: float, _model):
    """Where the calibrated member loses its chlorine: four EPANET runs (as calibrated, wall decay off, bulk decay
    off, no decay), cached per network, dose and member."""
    return loss_split(_model)


def age_test(path: str) -> dict | None:
    """The simulated test of an example network under the default truth (outputs/chem/water_age_<net>.json):
    the age range's label, its coverage of the true daily-mean and daily-max age, and how far the calibrated
    loss split was from the truth's.  An uploaded file has no such test: None."""
    try:
        with open(AGE_RESULTS.format(os.path.basename(path))) as fh:
            d = json.load(fh)
        t = d["by_truth"]["default"]
        return {"label": d["band_label"], "cov": t["band_coverage"], "cov_max": t["band_coverage_daily_max"],
                "split_median_diff": t["loss_split"]["median_abs_diff_mean"],
                "split_median_diff_max": t["loss_split"]["median_abs_diff_max"],
                "split_mae": t["loss_split"]["wall_share_mae"], "n_seeds": t["n_seeds"]}
    except (OSError, KeyError, ValueError, TypeError):
        return None


@st.cache_resource(show_spinner=False)
def demo_truth(path: str, dose: float, seed: int, kb: float, kw: float):
    return build_scenario(path, seed=seed, sample_hour=14, source_dose=dose, kb_per_day=kb, kw_m_per_day=kw)


try:
    wn_check = wntr.network.WaterNetworkModel(net_path if os.path.exists(net_path) else os.path.join(LIB, f"{net_path}.inp"))
except Exception as e:  # noqa: BLE001
    st.error(f"Could not read this .inp with EPANET/WNTR: {e}")
    st.stop()
if not source_nodes(wn_check):
    st.error("No reservoir, inflow junction or tank found — the model needs a place where water enters.")
    st.stop()

n_runs = int(np.prod([len(a) for a in GRIDS["full"]]))
with st.spinner(f"Running {n_runs} EPANET simulations of your model over the decay and hydraulic-mismatch grid "
                f"(once per network; cached afterwards)…"):
    sc, X = prepare(net_path, float(dose))
with st.spinner("Running 10 EPANET water-age simulations of your model (demand and roughness settings; once per network)…"):
    band = water_age(net_path, sc)
junctions = list(sc.junctions)

# samples table -----------------------------------------------------------------------------
if demo:
    tr = demo_truth(net_path, float(dose), int(demo_seed), demo_kb, demo_kw)
    rng = np.random.default_rng(int(demo_seed))
    js = list(rng.choice(junctions, int(n_demo), replace=False)) if n_demo else []
    hs = [int(h) for h in rng.choice(DAY_HOURS, int(n_demo))] if n_demo else []
    ys = [float(np.clip(tr.truth_by_hour.loc[h, j] + rng.normal(0, 0.03), 0.01, None)) for j, h in zip(js, hs)]
    default = pd.DataFrame({"junction": js, "hour": hs, "mg/L": np.round(ys, 2)})
else:
    tr = None
    default = pd.DataFrame({"junction": pd.Series(dtype=str), "hour": pd.Series(dtype=int), "mg/L": pd.Series(dtype=float)})

st.subheader("Grab samples")
st.caption("One row per sample: the junction ID from your model, the hour it was taken (0–23), the free chlorine reading.")
samples = st.data_editor(default, num_rows="dynamic", width='stretch', key=f"samples_{net_path}_{demo}_{demo and (n_demo, demo_seed)}",
                         column_config={"junction": st.column_config.SelectboxColumn("junction", options=junctions, required=True),
                                        "hour": st.column_config.NumberColumn("hour", min_value=0, max_value=23, step=1),
                                        "mg/L": st.column_config.NumberColumn("mg/L", min_value=0.0, max_value=5.0, step=0.01, format="%.2f")})
samples = samples.dropna()
samples = samples[samples.junction.isin(junctions)]
S = pd.DataFrame({"junction": samples.junction.astype(str), "hour": samples.hour.astype(int), "y": samples["mg/L"].astype(float)})

# model ------------------------------------------------------------------------------------
model = SimGP24(sc, X, seed=0, cache_dir=CACHE_DIR, threshold=float(threshold))
model.fit(S if len(S) else None)
hourly = model.predict_hours()
pmin = model.predict_daily_min()
route = plan_route(model, int(K), threshold=float(threshold), exclude=list(S.junction), with_age=True,
                   age_initial_share=band.initial_share)

# ----------------------------------------------------------------------------- headline numbers
view = st.radio("Show", ["Daily minimum (the compliance number)"] + [f"{h:02d}:00" for h in range(24)], horizontal=True, index=0)
if view.startswith("Daily"):
    med, lo, hi = pmin["median"], pmin["lo90"], pmin["hi90"]
    p_below = pmin["p_below"]          # the model's own draws, at the threshold it was built with
    label = "daily minimum"
else:
    h = int(view[:2]); z_mu, z_sd = hourly[0][h], hourly[1][h]
    med = pd.Series(np.exp(z_mu), index=junctions); lo = pd.Series(np.exp(z_mu - 1.645 * z_sd), index=junctions)
    hi = pd.Series(np.exp(z_mu + 1.645 * z_sd), index=junctions)
    p_below = pd.Series(norm.cdf((np.log(float(threshold)) - z_mu) / z_sd), index=junctions)
    label = view
flag = p_below > 0.5
frac_by_hour = pd.Series([float((hourly[0][h] < np.log(float(threshold))).mean()) for h in range(24)], index=range(24))
worst = int(frac_by_hour.idxmax())

c1, c2, c3, c4 = st.columns(4)
c1.metric("Samples used", len(S))
c2.metric(f"Junctions likely below {threshold:g} mg/L ({label})", f"{int(flag.sum())} of {len(junctions)}")
c3.metric("Worst hour of the day", f"{worst:02d}:00", f"{frac_by_hour[worst]:.0%} of junctions below by median")
if model.map_params_ is not None:
    kb, kw, g, dm, rm = model.map_params_
    c4.metric("Calibrated decay (bulk / wall)", f"{kb:.2f} /d · {kw:.2f} m/d", f"old-pipe factor γ={g:.1f}, demand ×{dm:.2f}, dose ×{model.map_dose_:.2f}")
else:
    c4.metric("Calibrated decay", "no samples yet", "map = your model's physics alone")


# ----------------------------------------------------------------------------- the four panels
def network_panel(ax, series, title, cmap, vrange, marks=None, mark_labels=None, colorbar_label=None):
    wntr.graphics.plot_network(sc.wn, node_attribute=series.to_dict(), node_size=max(12, int(4000 / len(junctions))),
                               node_cmap=cmap, node_range=vrange, ax=ax, link_width=0.6, add_colorbar=True, title=title)
    if marks:
        ax.scatter([sc.coords[j][0] for j in marks], [sc.coords[j][1] for j in marks], s=130, facecolors="none",
                   edgecolors="cyan", linewidths=1.8, zorder=5)
        if mark_labels:
            for j, lab in zip(marks, mark_labels):
                ax.annotate(lab, sc.coords[j], fontsize=8, color="darkcyan", xytext=(4, 3), textcoords="offset points")


def four_panels(figsize=(24, 6.5)):
    fig, axes = plt.subplots(1, 4, figsize=figsize)
    top = max(1.2, float(dose))
    network_panel(axes[0], med, f"Chlorine map — {label} (mg/L)\nfrom {len(S)} grab samples", "viridis", (0, top),
                  list(S.junction), [f"{h}h" for h in S.hour])
    network_panel(axes[1], hi - lo, "Uncertainty — width of the 90% band (mg/L)", "magma", (0, float((hi - lo).quantile(0.98))))
    network_panel(axes[2], p_below, f"P(residual < {threshold:g} mg/L) — {label}\n{int(flag.sum())} junctions likely below", "Reds", (0, 1))
    network_panel(axes[3], p_below, f"Sample here next: {len(route)} sites for next month's route", "Reds", (0, 1),
                  list(route.junction), [f"{k + 1}: {h:02d}h" for k, h in enumerate(route.hour)])
    fig.tight_layout()
    return fig


st.pyplot(four_panels(), width='stretch')

# ----------------------------------------------------------------------------- next samples with reasons
st.subheader("Sample here next")
st.dataframe(route.assign(hour=[f"{h:02d}:00" for h in route.hour])[["junction", "hour", "p_below", "water_age_h", "reason"]]
             .rename(columns={"p_below": f"P(daily min < {threshold:g})", "water_age_h": "water age (h)"}), width='stretch', hide_index=True)

# ----------------------------------------------------------------------------- water age, and where chlorine is lost
age_max = band.nominal.max()
age_lo, age_hi = band.daily_range("max")
j_old, h_old = oldest_water(band.nominal)
test = age_test(net_path)
range_label = test["label"] if test else RANGE_LABEL
# every age is from a 7-day simulation that starts with the file's water in the pipes and tanks; where more than
# INITIAL_SHARE_MAX of the water at a junction is still that starting water, its age is a lower bound
lb_max, lb_mean = band.lower_bound("max"), band.lower_bound("mean")
share_old = float(band.initial_share[j_old].iloc[int(np.argmax(band.nominal[j_old].values))])
if lb_max[j_old]:
    oldest_line = (f"oldest water: at least {h_old:.0f} h at junction {j_old} (your model; at that hour {share_old:.0%} of "
                   f"the water there is still the 7-day simulation's starting water, so it is older than the simulation can show)")
else:
    oldest_line = (f"oldest water: {h_old:.0f} h at junction {j_old} (your model; {age_lo[j_old]:.0f} to {age_hi[j_old]:.0f} h "
                   f"across demand and roughness errors)")
lower_note = (f"Every age comes from a 7-day simulation that starts with your file's water in the pipes and tanks. Where "
              f"more than {INITIAL_SHARE_MAX:.0%} of a junction's water is still that starting water, its age is a lower bound: "
              f"{int(lb_max.sum())} of {len(junctions)} junctions at their oldest hour, {int(lb_mean.sum())} on the day's average.")
age_col, loss_col = st.columns([2, 1])
with age_col:
    st.subheader("Water age")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    network_panel(axes[0], age_max, "Water age at its oldest hour of the day (h)\nyour EPANET model, 7-day simulation", "YlOrBr",
                  (0, float(age_max.quantile(0.98))))
    network_panel(axes[1], age_hi - age_lo, f"How much that age moves (h): {range_label}\n9 demand and roughness settings",
                  "magma", (0, max(float((age_hi - age_lo).quantile(0.98)), 1.0)))
    fig.tight_layout()
    st.pyplot(fig, width='stretch')
    if test and range_label == RANGE_LABEL:
        tested = (f" In simulation on this network it held the true daily-maximum age at {test['cov_max']:.0%} of junctions "
                  f"(daily mean: {test['cov']:.0%}), so it is shown as a range, not a 90% band: errors in single junctions' "
                  f"demands are not among the 9 settings.")
    elif test:
        tested = f" In simulation on this network it held the true daily-maximum age at {test['cov_max']:.0%} of junctions."
    else:
        tested = " It has not been tested against a simulated truth for this file, so it is shown as a range, not a 90% band."
    st.caption(f"Hours since the water left the source, from your EPANET file: {oldest_line}. {lower_note} Right: how far "
               f"each junction's age moves across the 9 settings the model already weighs (demand x0.85, x1, x1.15; pipe "
               f"roughness x0.9, x1, x1.1).{tested}")
with loss_col:
    st.subheader("Where chlorine is lost")
    if model.map_params_ is None:
        st.info("Enter at least one grab sample: the split uses the decay rates your samples calibrate.")
    else:
        split = chlorine_loss(net_path, float(dose), tuple(model.map_params_), float(model.map_dose_), model)
        ws = split.wall_share.dropna()
        if ws.empty:
            st.info("Chlorine loss is under 2% at every junction for these decay rates: there is nothing to split.")
        else:
            fig, ax = plt.subplots(figsize=(7, 5.6))
            network_panel(ax, ws, "Share of the chlorine loss at pipe walls\n0 = all in the water, 1 = all at the walls", "RdYlBu_r", (0, 1))
            fig.tight_layout()
            st.pyplot(fig, width='stretch')
            na = split.nonadditivity.dropna()
            n_small = len(junctions) - len(ws)
            accuracy = (f" In simulation on this network ({test['n_seeds']} scenarios, 8 samples each) the calibrated median was "
                        f"off from the truth's by {test['split_median_diff']:.2f} on average and by up to {test['split_median_diff_max']:.2f} "
                        f"in one scenario." if test else " It has not been tested against a simulated truth for this file.")
            st.caption(f"For the single most likely decay rates your samples calibrate (bulk {model.map_params_[0]:.2f} /day, wall "
                       f"{model.map_params_[1]:.2f} m/day, old-pipe factor {model.map_params_[2]:.1f}): a median {ws.median():.0%} of the "
                       f"chlorine lost on the way to a junction is lost at the pipe walls, the rest in the water (organics and other "
                       f"reactants). It is an estimate: with few samples, bulk and wall decay can trade off against each other."
                       f"{accuracy} Exact along a single path from the source; approximate where flows mix (here the two parts add "
                       f"up to {na.min():.2f} to {na.max():.2f} of the total). Not coloured (under 2% loss): {n_small} "
                       f"junction{'' if n_small == 1 else 's'}. Field studies report wall loss up to 97% of the total in old pipes "
                       f"(Maleki et al. 2023).")

# ----------------------------------------------------------------------------- worst hour + demo truth
left, right = st.columns([1, 1])
with left:
    st.subheader("When is the network worst?")
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.bar(frac_by_hour.index, frac_by_hour.values * 100, color=["tab:red" if h >= 18 or h < 6 else "tab:blue" for h in frac_by_hour.index])
    ax.axvspan(6.5, 17.5, color="gold", alpha=0.15, label="when you sample")
    ax.set(xlabel="hour of day", ylabel=f"% junctions below {threshold:g} mg/L (median)"); ax.legend(loc="upper center")
    st.pyplot(fig, width='stretch')
with right:
    if tr is not None:
        st.subheader("Demo only: reveal the hidden truth")
        if st.toggle("Show the true daily-minimum map and score the flags", value=False):
            tv = tr.truth_daily_min < float(threshold)
            uns = [j for j in junctions if j not in set(S.junction)]
            tp = int((tv.loc[uns] & (pmin.loc[uns, "p_below"] > 0.5)).sum()); fn = int((tv.loc[uns] & ~(pmin.loc[uns, "p_below"] > 0.5)).sum())
            fp = int((~tv.loc[uns] & (pmin.loc[uns, "p_below"] > 0.5)).sum())
            st.write(f"True daily-minimum violations: **{int(tv.sum())} of {len(junctions)}** junctions. On the unsampled junctions the map "
                     f"found **{tp} of {tp + fn}** (recall {tp / max(tp + fn, 1):.0%}) with {fp} false alarms. "
                     f"Mean of your samples: {S.y.mean():.2f} mg/L — which flags {'everything' if S.y.mean() < threshold else 'nothing'}.")
            fig, ax = plt.subplots(figsize=(8, 5))
            network_panel(ax, tr.truth_daily_min, "TRUE daily minimum (hidden from the model)", "viridis", (0, max(1.2, float(dose))))
            st.pyplot(fig, width='stretch')


# ----------------------------------------------------------------------------- one-page PDF
def pdf_bytes() -> bytes:
    fig = plt.figure(figsize=(16.5, 11.7))  # A3 landscape-ish, prints fine on A4
    gs = fig.add_gridspec(3, 4, height_ratios=[0.35, 1.5, 1.1])
    ax = fig.add_subplot(gs[0, :]); ax.axis("off")
    name = os.path.basename(net_path)
    ax.text(0, 0.9, f"ResidualMap — monthly chlorine report — {name}", fontsize=18, fontweight="bold", va="top")
    calib = (f"calibrated decay: bulk {model.map_params_[0]:.2f} /day, wall {model.map_params_[1]:.2f} m/day, old-pipe factor {model.map_params_[2]:.1f}; "
             f"demand ×{model.map_params_[3]:.2f}, roughness ×{model.map_params_[4]:.2f}, dose ×{model.map_dose_:.2f}") if model.map_params_ is not None else "no samples yet: map from the model's physics alone"
    ax.text(0, 0.45, f"{len(S)} grab samples · dose {dose:g} mg/L · minimum residual {threshold:g} mg/L · "
                     f"{int((pmin['p_below'] > 0.5).sum())} of {len(junctions)} junctions likely below the minimum at their daily minimum · "
                     f"worst hour {worst:02d}:00\n{calib}\n{oldest_line}", fontsize=10.5, va="top")
    axes = [fig.add_subplot(gs[1, i]) for i in range(4)]
    top = max(1.2, float(dose))
    network_panel(axes[0], med, f"Chlorine — {label} (mg/L)", "viridis", (0, top), list(S.junction), [f"{h}h" for h in S.hour])
    network_panel(axes[1], hi - lo, "90% band width (mg/L)", "magma", (0, float((hi - lo).quantile(0.98))))
    network_panel(axes[2], p_below, f"P(residual < {threshold:g} mg/L)", "Reds", (0, 1))
    network_panel(axes[3], p_below, f"Next route: {len(route)} sites", "Reds", (0, 1), list(route.junction), [f"{k + 1}: {h:02d}h" for k, h in enumerate(route.hour)])
    ax = fig.add_subplot(gs[2, :2]); ax.axis("off")
    ax.set_title("Sample here next", loc="left", fontsize=12, fontweight="bold")
    # the reason without its P(...) part (its own column), wrapped so none of it is cut off
    rows = [[str(j), f"{h:02d}:00", f"{p:.2f}", "\n".join(textwrap.wrap(w, PDF_WHY_CHARS))]
            for j, h, p, w in zip(route.junction, route.hour, route.p_below, route.why)]
    tbl = ax.table(cellText=rows, colLabels=["junction", "hour", f"P(min<{threshold:g})", "why"], loc="upper left", cellLoc="left",
                   colWidths=[0.1, 0.1, 0.12, 0.68])
    n_lines = [1] + [r[3].count("\n") + 1 for r in rows]
    unit = 0.97 / (sum(n_lines) + 0.6 * len(n_lines))       # the table fills the panel's height
    for (r, c), cell in tbl.get_celld().items():
        cell.set_height(unit * (n_lines[r] + 0.6))
    # each text line gets about unit of the table's height: a long route gets a smaller font, never overlapping rows
    line_pt = unit * ax.get_position().height * fig.get_figheight() * 72
    tbl.auto_set_font_size(False); tbl.set_fontsize(min(7.5, 0.95 * line_pt))
    ax = fig.add_subplot(gs[2, 2:])
    ax.bar(frac_by_hour.index, frac_by_hour.values * 100, color=["tab:red" if h >= 18 or h < 6 else "tab:blue" for h in frac_by_hour.index])
    ax.axvspan(6.5, 17.5, color="gold", alpha=0.15, label="sampling window")
    ax.set(xlabel="hour of day", ylabel=f"% junctions below {threshold:g} mg/L", title="When the network is worst"); ax.legend()
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="pdf"); plt.close(fig)
    return buf.getvalue()


st.download_button("Download the one-page PDF report", data=pdf_bytes(), file_name=f"residualmap_{os.path.basename(net_path)}.pdf", mime="application/pdf")
st.caption("Every number on this page is computed on this laptop from your .inp and your samples. Nothing is uploaded anywhere.")
