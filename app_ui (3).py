"""
Streamlit interface for the Spud–SRI leg penetration calculator.

This app implements the V50 logic for leg penetration analysis using the
`lpa_v50` computation engine.  It provides a simple web interface
that allows engineers to describe a spudcan and a soil profile
directly in the browser (without uploading Excel files) and then
visualise the resulting penetration curves.  Capacities are
calculated for undrained (clay) and drained (sand) cases where
relevant soil data exist; the controlling capacity is taken as the
minimum of the available modes at each penetration depth.  A single
penetration prediction is then derived from this controlling curve.

Key features:
  • Input spudcan geometry (diameter, area, tip offset and preload).
  • Define an arbitrary number of soil layers with depths, soil
    classification, and constant values for submerged unit weight
    (γ′), friction angle (φ′) and undrained shear strength (Su).  If
    either φ′ or Su is omitted for a layer, the corresponding sand or
    clay bearing capacity will not be computed in that layer.
  • Choose conservative options: use the minimum of point and B/2‑average
    Su, apply a 5° reduction to φ′, multiply capacities by 0.8 to
    account for a windward leg, and enforce a geometric trigger for
    squeezing failures.  These flags mirror the options available in
    the original VBA macro.
  • Compute and display the penetration curve, highlighting the
    controlling capacity and the governing mode (clay or sand) at
    each depth.  A Plotly graph plots the controlling curve along
    with a vertical line representing the preload.  The plot size has
    been reduced for better presentation in a browser.
  • Report the predicted penetration depth at which the controlling
    capacity equals the preload and the corresponding tip penetration
    after adding the tip offset.  Unlike earlier versions, only a
    single penetration value is reported (rather than separate clay
    and sand bounds) because the controlling capacity already
    incorporates whichever mode governs at each depth.

Run this script with `streamlit run app_ui.py` to launch the app.
"""
# app_ui.py
# Streamlit UI for spud-SRI / Leg Penetration (SNAME-style), no Plotly required.

import io
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from lpa_v50 (3) import (
    Spudcan, SoilPoint, SoilLayer,
    compute_envelopes, penetration_results,
    USE_MIN_CU_POINT_AVG_DEFAULT, APPLY_PHI_REDUCTION_DEFAULT,
    APPLY_WINDWARD_FACTOR_DEFAULT, APPLY_SQUEEZE_TRIGGER_DEFAULT,
)

st.set_page_config(
    page_title="spud-SRI / Leg Penetration (SNAME)",
    page_icon="🦀",
    layout="centered",
)

st.title("spud-SRI · Leg Penetration (SNAME-style)")

with st.sidebar:
    st.subheader("Analysis switches")
    use_min_cu = st.checkbox("Use min(Su point, Su B/2 avg)", value=USE_MIN_CU_POINT_AVG_DEFAULT,
                             help="Conservative default for undrained strength.")
    phi_reduce = st.checkbox("Apply 5° reduction to ϕ′", value=APPLY_PHI_REDUCTION_DEFAULT,
                             help="Optional conservatism for sands.")
    windward80 = st.checkbox("Windward factor 0.8 on REAL", value=APPLY_WINDWARD_FACTOR_DEFAULT,
                             help="Applies 0.8 to the governing REAL capacity only.")
    squeeze_trig = st.checkbox("Enforce squeezing geometric trigger", value=APPLY_SQUEEZE_TRIGGER_DEFAULT)
    dz = st.number_input("Depth step Δz (m)", value=0.25, min_value=0.05, max_value=2.0, step=0.05)
    dmax = st.number_input("Max analysis depth (m)", value=30.0, min_value=5.0, max_value=200.0, step=1.0)

st.markdown("#### Spudcan inputs")
cols = st.columns(5)
rig = cols[0].text_input("Rig name", "Rig-1")
B   = cols[1].number_input("Diameter B (m)", value=8.0, min_value=0.1, step=0.1)
A   = cols[2].number_input("Area A (m²)", value=float(np.pi*(B**2)/4.0), min_value=0.01, step=0.1,
                           help="Projected area of widest section. Defaults to πB²/4.")
tip = cols[3].number_input("Tip elevation (m)", value=1.5, min_value=0.0, step=0.1,
                           help="Distance from tip to widest section; added to analysis depth for tip penetration.")
Pmn = cols[4].number_input("Preload per leg (MN)", value=80.0, min_value=1.0, step=1.0)

spud = Spudcan(rig_name=rig, B=B, A=A, tip_elev=tip, preload_MN=Pmn)

st.markdown("#### Soil profile (mimics Excel sheets)")
st.caption("Add layers from seabed downward. Provide data only where it exists; the app will not fabricate a second envelope.")

# Layer builder
if "layers" not in st.session_state:
    st.session_state.layers = []

def _add_layer():
    st.session_state.layers.append({
        "name": f"Layer {len(st.session_state.layers)+1}",
        "z_top": 0.0 if not st.session_state.layers else st.session_state.layers[-1]["z_bot"],
        "z_bot": (0.0 if not st.session_state.layers else st.session_state.layers[-1]["z_bot"]) + 2.0,
        "type": "clay",
        "gamma_pairs": "0,10.0; 2,10.0",
        "su_pairs":    "0,30;   2,35",
        "phi_pairs":   "",
    })

def _parse_pairs(s: str):
    s = (s or "").strip()
    if not s:
        return []
    pts = []
    for tok in s.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        d, v = tok.split(",")
        pts.append(SoilPoint(float(d.strip()), float(v.strip())))
    return pts

st.button("➕ Add layer", on_click=_add_layer)
for i, L in enumerate(st.session_state.layers):
    with st.expander(f"Layer {i+1}", expanded=True):
        c1, c2, c3, c4 = st.columns([1.2,1.2,1.2,1.2])
        L["name"]  = c1.text_input("Name", L["name"], key=f"name{i}")
        L["z_top"] = c2.number_input("z_top (m)", value=float(L["z_top"]), key=f"z1{i}", step=0.1)
        L["z_bot"] = c3.number_input("z_bot (m)", value=float(L["z_bot"]), key=f"z2{i}", step=0.1)
        L["type"]  = c4.selectbox("Type", ["clay","sand","silt","unknown"], index=["clay","sand","silt","unknown"].index(L["type"]), key=f"type{i}")
        L["gamma_pairs"] = st.text_input("γ′ pairs 'z,val; z,val; ...'  (kN/m³)", L["gamma_pairs"], key=f"g{i}")
        L["su_pairs"]    = st.text_input("Su pairs 'z,val; z,val; ...'  (kPa)", L["su_pairs"],    key=f"su{i}")
        L["phi_pairs"]   = st.text_input("ϕ′ pairs 'z,val; z,val; ...' (deg)", L["phi_pairs"],   key=f"phi{i}")

st.divider()

do_run = st.button("Run analysis", type="primary")
if do_run:
    # Build SoilLayer list
    layers = []
    for L in st.session_state.layers:
        layers.append(
            SoilLayer(
                name=L["name"],
                z_top=float(L["z_top"]),
                z_bot=float(L["z_bot"]),
                soil_type=L["type"],
                gamma=_parse_pairs(L["gamma_pairs"]),
                su=_parse_pairs(L["su_pairs"]),
                phi=_parse_pairs(L["phi_pairs"]),
            )
        )

    if not layers:
        st.error("Please add at least one soil layer.")
        st.stop()

    # Compute
    df = compute_envelopes(
        spud=spud,
        layers=layers,
        max_depth=dmax,
        dz=dz,
        use_min_cu=use_min_cu,
        phi_reduction=phi_reduce,
        windward_factor=windward80,
        squeeze_trigger=squeeze_trig,
        meyerhof_table=None,  # could be made user-editable later
    )

    # Penetration results (analysis depth and tip depth)
    pen = penetration_results(spud, df)

    # --- Results summary ---
    with st.container():
        st.subheader("Results")
        cols = st.columns(3)
        cols[0].metric("Preload per leg", f"{spud.preload_MN:.2f} MN")
        if pen["tip_range_min"] is not None and pen["tip_range_max"] is not None and pen["tip_range_min"] != pen["tip_range_max"]:
            cols[1].metric("Tip penetration range", f"{pen['tip_range_min']:.2f} – {pen['tip_range_max']:.2f} m")
        elif pen["tip_range_min"] is not None:
            cols[1].metric("Tip penetration", f"{pen['tip_range_min']:.2f} m")
        else:
            cols[1].metric("Tip penetration", "—")
        gov_at_P = ""
        # try to infer which envelope governs near the interpolated depths
        cols[2].metric("Note", "REAL is the min of valid envelopes at each depth.")

    # --- Compact, high-DPI chart (Matplotlib) ---
    fig, ax = plt.subplots(figsize=(4.2, 6.2), dpi=200)   # small footprint, crisp
    ax.plot(df["idle_clay_MN"], df["depth"], lw=0.8, ls="-", color="0.6", label="Idle Clay")
    ax.plot(df["idle_sand_MN"], df["depth"], lw=0.8, ls="--", color="0.6", label="Idle Sand")
    ax.plot(df["real_MN"], df["depth"], lw=1.8, color="#0033cc", label="REAL (governing)")

    # preload line
    ax.axvline(spud.preload_MN, ymin=0, ymax=1, color="red", lw=1.2, ls="--", label="Preload")

    ax.set_xlabel("Leg load (MN)")
    ax.set_ylabel("Penetration of widest section (m)")
    ax.invert_yaxis()
    ax.grid(True, ls=":", lw=0.6, color="0.7")
    ax.set_xlim(left=0)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    st.pyplot(fig, clear_figure=True)

    # --- Table & downloads ---
    st.subheader("Detailed table")
    st.dataframe(df, use_container_width=True, height=320)

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("Download CSV", data=csv, file_name=f"{spud.rig_name}_results.csv", mime="text/csv")

    # small textual summary for the report
    st.caption(
        "Notes: Idle lines are computed only where the corresponding soil property exists. "
        "REAL equals the minimum of valid candidates (clay/sand), after applying special-mode reductions "
        "(squeezing or punch-through) when their geometric/material triggers are satisfied. "
        "Backflow sets p′=0 in clay where the Meyerhof stability condition is exceeded."
    )
