
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from lpa_v50 import (
    Spudcan, SoilLayer, SoilDataPoint,
    compute_envelopes, penetration_at_load_MN
)

st.set_page_config(page_title="Spud-SRI | Leg Penetration (SNAME)", layout="wide")
st.title("Spud-SRI • Leg Penetration based on SNAME")

st.caption("Manual input UI mimicking the Excel sheets. Define layers, enter global property profiles for γ′, φ′, and Su, and compute SNAME-compliant penetration envelopes.")

with st.sidebar:
    use_min_cu = st.checkbox("Use min(cu_point, cu_avg) for clay", True)
    phi_reduction = st.checkbox("Apply 5° reduction to φ′", False)
    windward = st.checkbox("Apply windward 0.8 factor", False)
    squeeze_trigger_info = st.caption("Geometric squeezing trigger is enabled in engine.")
    dz = st.number_input("Depth step Δz (m)", value=0.25, min_value=0.05, step=0.05)
    dmax = st.number_input("Max analysis depth (m)", value=30.0, min_value=5.0, step=1.0)

col1, col2 = st.columns([1.1, 1.0])

with col1:
    st.subheader("Spud data")
    rig = st.text_input("Rig Name", "Rig")
    B = st.number_input("Diameter B (m)", min_value=0.0, value=8.0, step=0.1)
    A_default = 0.25*math.pi*B*B if B>0 else 0.0
    A = st.number_input("Area A (m²)", min_value=0.0, value=A_default, step=0.1)
    tip = st.number_input("Tip Elevation from widest section (m)", value=0.0, step=0.1)
    preloadMN = st.number_input("Max Preload per leg (MN)", value=80.0, min_value=0.0, step=1.0)

    st.subheader("Layers (A,B,C,E)")
    layers_df = st.data_editor(pd.DataFrame({
        "Layer": ["Layer 1", "Layer 2"],
        "Top": [0.0, 6.0],
        "Bottom": [6.0, 20.0],
        "Description": ["stiff clay", "medium dense sand"]
    }), use_container_width=True, num_rows="dynamic")

with col2:
    st.subheader("Global profiles")
    gamma_df = st.data_editor(pd.DataFrame({
        "z_gamma": [0.0, 3.0, 6.0, 10.0, 16.0],
        "gamma_kN_m3": [8.5, 9.0, 9.5, 10.0, 10.5]
    }), use_container_width=True, num_rows="dynamic")
    phi_df = st.data_editor(pd.DataFrame({
        "z_phi": [0.0, 6.0, 6.1, 12.0, 18.0],
        "phi_deg": [0.0, 0.0, 32.0, 33.0, 34.0]
    }), use_container_width=True, num_rows="dynamic")
    su_df = st.data_editor(pd.DataFrame({
        "z_su": [0.0, 3.0, 5.5],
        "su_kPa": [25.0, 35.0, 45.0]
    }), use_container_width=True, num_rows="dynamic")
    st.subheader("Meyerhof N vs D/B (optional)")
    mey_df = st.data_editor(pd.DataFrame({
        "D_over_B": [0.2, 1.0],
        "N": [1.5, 6.0]
    }), use_container_width=True, num_rows="dynamic")

if st.button("Compute"):
    if len(layers_df) == 0:
        st.error("Please define at least one layer.")
        st.stop()

    soil_layers = []
    for _, r in layers_df.iterrows():
        name = str(r.get("Layer","Layer"))
        top = float(r.get("Top", 0.0))
        bot = float(r.get("Bottom", max(top+0.5, top)))
        desc = str(r.get("Description","")).lower()
        if "clay" in desc:
            soil_type = "clay"
        elif "silt" in desc:
            soil_type = "silt"
        elif "sand" in desc:
            soil_type = "sand"
        else:
            soil_type = "unknown"
        gamma_prof, phi_prof, su_prof = [], [], []
        for _, g in gamma_df.iterrows():
            try:
                d = float(g["z_gamma"]); v = float(g["gamma_kN_m3"])
                if d >= top and d < bot:
                    gamma_prof.append(SoilDataPoint(d, v))
            except: pass
        for _, p in phi_df.iterrows():
            try:
                d = float(p["z_phi"]); v = float(p["phi_deg"])
                if d >= top and d < bot:
                    phi_prof.append(SoilDataPoint(d, v))
            except: pass
        for _, s in su_df.iterrows():
            try:
                d = float(s["z_su"]); v = float(s["su_kPa"])
                if d >= top and d < bot and v>0:
                    su_prof.append(SoilDataPoint(d, v))
            except: pass
        soil_layers.append(SoilLayer(name=name, top=top, bottom=bot, soil_type=soil_type,
                                     gamma_profile=sorted(gamma_prof, key=lambda x:x.depth),
                                     su_profile=sorted(su_prof, key=lambda x:x.depth),
                                     phi_profile=sorted(phi_prof, key=lambda x:x.depth)))
    spud = Spudcan(rig, B, A, tip, preloadMN*1000.0)

    mey = None
    try:
        mm = mey_df.dropna()
        if {"D_over_B","N"}.issubset(set(mm.columns)):
            mey = pd.DataFrame({"D_over_B": pd.to_numeric(mm["D_over_B"], errors="coerce"),
                                "N": pd.to_numeric(mm["N"], errors="coerce")}).dropna()
    except Exception:
        mey = None

    df = compute_envelopes(spud, soil_layers, mey, dz=dz, dmax=dmax,
                           use_min_cu_point_avg=use_min_cu,
                           apply_phi_reduction=phi_reduction,
                           windward_factor=windward)

    pMN = spud.max_preload_kN/1000.0
    pen_clay = penetration_at_load_MN(df, pMN, "real_clay_MN")
    pen_sand = penetration_at_load_MN(df, pMN, "real_sand_MN")
    tip_clay = spud.tip_elev + pen_clay if not math.isnan(pen_clay) else float('nan')
    tip_sand = spud.tip_elev + pen_sand if not math.isnan(pen_sand) else float('nan')

    st.success(f"Tip penetration range at {pMN:.2f} MN preload: {tip_clay:.2f} m – {tip_sand:.2f} m")

    import matplotlib
    matplotlib.use("Agg")
    fig, ax = plt.subplots(figsize=(7,9), dpi=120)
    ax.plot(df["real_clay_MN"], df["depth"], label="REAL - Clay")
    ax.plot(df["real_sand_MN"], df["depth"], label="REAL - Sand", linestyle="--")
    ax.plot([pMN, pMN], [0, df["depth"].max()], label="Max Preload", color="red", linestyle=":")
    ax.set_xlabel("Leg Load (MN)")
    ax.set_ylabel("Penetration (m)")
    ax.invert_yaxis()
    ax.grid(True, linestyle="--", linewidth=0.6)
    ax.legend()
    st.pyplot(fig, clear_figure=True)

    st.subheader("Results table")
    st.dataframe(df, use_container_width=True)
    st.download_button("Download CSV", df.to_csv(index=False).encode("utf-8"),
                       file_name=f"{spud.rig_name}_lpa_results.csv", mime="text/csv")
