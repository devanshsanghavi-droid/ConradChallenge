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
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import wntr
from scipy.stats import norm

from residualmap.features import build_features
from residualmap.route import plan_route
from residualmap.simgp import DAY_HOURS, GRIDS, SimGP24, simulator_grid_24h
from residualmap.simulate import LIB, build_scenario, nominal_scenario, source_nodes

warnings.filterwarnings("ignore")
UPLOAD_DIR = "outputs/app_uploads"
CACHE_DIR = "outputs/cache"
EXAMPLES = {"Net3 (North Marin, CA — 92 junctions)": "Net3", "Net2 (tank-fed — 35 junctions)": "Net2",
            "ky4 (Kentucky — 959 junctions; first run takes ~15 min)": "ky4"}

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
model = SimGP24(sc, X, seed=0, cache_dir=CACHE_DIR)
model.fit(S if len(S) else None)
hourly = model.predict_hours()
pmin = model.predict_daily_min()
route = plan_route(model, int(K), threshold=float(threshold), exclude=list(S.junction))

# ----------------------------------------------------------------------------- headline numbers
view = st.radio("Show", ["Daily minimum (the compliance number)"] + [f"{h:02d}:00" for h in range(24)], horizontal=True, index=0)
if view.startswith("Daily"):
    med, lo, hi = pmin["median"], pmin["lo90"], pmin["hi90"]
    p_below = pmin["p_below"] if abs(float(threshold) - 0.2) < 1e-9 else SimGP24.p_below(pmin, float(threshold))
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
st.dataframe(route.assign(hour=[f"{h:02d}:00" for h in route.hour])[["junction", "hour", "p_below", "reason"]]
             .rename(columns={"p_below": f"P(daily min < {threshold:g})"}), width='stretch', hide_index=True)

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
                     f"worst hour {worst:02d}:00\n{calib}", fontsize=10.5, va="top")
    axes = [fig.add_subplot(gs[1, i]) for i in range(4)]
    top = max(1.2, float(dose))
    network_panel(axes[0], med, f"Chlorine — {label} (mg/L)", "viridis", (0, top), list(S.junction), [f"{h}h" for h in S.hour])
    network_panel(axes[1], hi - lo, "90% band width (mg/L)", "magma", (0, float((hi - lo).quantile(0.98))))
    network_panel(axes[2], p_below, f"P(residual < {threshold:g} mg/L)", "Reds", (0, 1))
    network_panel(axes[3], p_below, f"Next route: {len(route)} sites", "Reds", (0, 1), list(route.junction), [f"{k + 1}: {h:02d}h" for k, h in enumerate(route.hour)])
    ax = fig.add_subplot(gs[2, :2]); ax.axis("off")
    ax.set_title("Sample here next", loc="left", fontsize=12, fontweight="bold")
    rows = [[str(j), f"{h:02d}:00", f"{p:.2f}", r[:110]] for j, h, p, r in zip(route.junction, route.hour, route.p_below, route.reason)]
    tbl = ax.table(cellText=rows, colLabels=["junction", "hour", f"P(min<{threshold:g})", "why"], loc="upper left", cellLoc="left",
                   colWidths=[0.1, 0.1, 0.12, 0.68])
    tbl.auto_set_font_size(False); tbl.set_fontsize(7.5); tbl.scale(1, 1.25)
    ax = fig.add_subplot(gs[2, 2:])
    ax.bar(frac_by_hour.index, frac_by_hour.values * 100, color=["tab:red" if h >= 18 or h < 6 else "tab:blue" for h in frac_by_hour.index])
    ax.axvspan(6.5, 17.5, color="gold", alpha=0.15, label="sampling window")
    ax.set(xlabel="hour of day", ylabel=f"% junctions below {threshold:g} mg/L", title="When the network is worst"); ax.legend()
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="pdf"); plt.close(fig)
    return buf.getvalue()


st.download_button("Download the one-page PDF report", data=pdf_bytes(), file_name=f"residualmap_{os.path.basename(net_path)}.pdf", mime="application/pdf")
st.caption("Every number on this page is computed on this laptop from your .inp and your samples. Nothing is uploaded anywhere.")
