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

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

from typing import List

from lpa_v50 import (
    SoilLayer,
    SoilDataPoint,
    Spudcan,
    calculate_penetration_curve,
    interpolate_penetration,
)


def build_manual_layers(num_layers: int) -> List[SoilLayer]:
    """Collect soil layer definitions from user inputs.

    Each layer is assumed to have constant soil properties across its
    thickness.  Users provide top and bottom depths, choose a soil
    type (clay, sand, silt or unknown) and optionally specify
    numerical values for submerged unit weight (gamma′), friction
    angle (φ′) and undrained shear strength (Su).  If a value is
    omitted (left blank), the corresponding property is not used in
    the calculation for that layer.  Internally, each property is
    stored as a profile with two points (start and end) to simplify
    interpolation.
    """
    layers: List[SoilLayer] = []
    for idx in range(num_layers):
        st.subheader(f"Layer {idx + 1}")
        name = st.text_input(f"Name for layer {idx + 1}", value=f"Layer{idx + 1}", key=f"name_{idx}")
        col1, col2, col3 = st.columns(3)
        with col1:
            top_depth = st.number_input(
                f"Top depth (m) – Layer {idx + 1}", min_value=0.0, step=0.1, format="%.2f", key=f"top_{idx}"
            )
        with col2:
            bottom_depth = st.number_input(
                f"Bottom depth (m) – Layer {idx + 1}", min_value=0.0, step=0.1, format="%.2f", key=f"bot_{idx}"
            )
        with col3:
            soil_type = st.selectbox(
                f"Soil type – Layer {idx + 1}",
                options=["clay", "sand", "silt", "unknown"],
                index=0 if idx == 0 else 1,
                key=f"soil_{idx}",
            )
        # Ensure bottom depth is greater than top depth
        if bottom_depth <= top_depth:
            st.error(f"For layer {idx + 1}, bottom depth must exceed top depth.")
        col4, col5, col6 = st.columns(3)
        with col4:
            gamma = st.number_input(
                f"Submerged unit weight γ′ (kN/m³) – Layer {idx + 1}",
                min_value=0.0,
                step=0.1,
                format="%.2f",
                key=f"gamma_{idx}",
            )
        with col5:
            phi = st.text_input(
                f"Friction angle φ′ (°) – Layer {idx + 1} (leave blank if none)",
                value="",
                key=f"phi_{idx}",
            )
        with col6:
            su = st.text_input(
                f"Undrained shear strength Su (kPa) – Layer {idx + 1} (leave blank if none)",
                value="",
                key=f"su_{idx}",
            )
        # Build property profiles
        unit_weight_profile = []
        phi_profile = []
        su_profile = []
        # Use two points at top and bottom so interpolation returns a constant value
        # across the layer.  If the user leaves phi or su blank, the profile
        # remains empty and that property is ignored for this layer.
        unit_weight_profile = [
            SoilDataPoint(depth=top_depth, value=gamma),
            SoilDataPoint(depth=bottom_depth, value=gamma),
        ]
        try:
            phi_val = float(phi) if phi.strip() != "" else None
        except Exception:
            phi_val = None
        if phi_val is not None and phi_val > 0:
            phi_profile = [
                SoilDataPoint(depth=top_depth, value=phi_val),
                SoilDataPoint(depth=bottom_depth, value=phi_val),
            ]
        try:
            su_val = float(su) if su.strip() != "" else None
        except Exception:
            su_val = None
        if su_val is not None and su_val > 0:
            su_profile = [
                SoilDataPoint(depth=top_depth, value=su_val),
                SoilDataPoint(depth=bottom_depth, value=su_val),
            ]
        layers.append(
            SoilLayer(
                name=name,
                top_depth=top_depth,
                bottom_depth=bottom_depth,
                soil_type=soil_type,
                unit_weight_profile=unit_weight_profile,
                su_profile=su_profile,
                phi_profile=phi_profile,
            )
        )
    return layers


def main() -> None:
    st.set_page_config(page_title="Spud-SRI | Leg Penetration (SNAME)", layout="centered")
    st.title("Spud-SRI | Leg Penetration (SNAME)")

    st.markdown(
        "This tool computes the leg penetration behaviour of a spudcan "
        "based on the SNAME guidelines.  Enter the geometry of your spudcan "
        "and define the soil profile layer by layer.  The app will calculate "
        "bearing capacities for clay and sand where appropriate, apply "
        "squeezing and punch-through criteria, and return a single controlling "
        "penetration estimate."
    )

    # Spudcan input section
    st.header("Spudcan parameters")
    rig_name = st.text_input("Rig name", value="Rig1")
    colA, colB, colC, colD = st.columns(4)
    with colA:
        diameter = st.number_input(
            "Diameter B (m)", min_value=0.0, step=0.1, format="%.2f", value=10.0
        )
    with colB:
        area = st.number_input(
            "Area A (m²)", min_value=0.0, step=1.0, format="%.2f", value=78.5
        )
    with colC:
        tip_elevation = st.number_input(
            "Tip offset from widest section (m)",
            help="Distance between the widest section and the tip. "
            "Positive values indicate that the tip lies below the widest section.",
            format="%.2f",
            value=0.0,
        )
    with colD:
        preload_MN = st.number_input(
            "Maximum preload per leg (MN)", min_value=0.0, step=1.0, format="%.2f", value=50.0
        )
    # Number of soil layers
    st.header("Soil profile")
    num_layers = st.number_input(
        "Number of soil layers", min_value=1, max_value=20, step=1, value=3, format="%d"
    )
    layers = build_manual_layers(int(num_layers))

    # Calculation options
    st.header("Calculation options")
    use_min_cu = st.checkbox(
        "Use minimum of point Su and B/2‑average Su (conservative)",
        value=True,
    )
    apply_phi_red = st.checkbox(
        "Apply 5° reduction to φ′ for sand layers", value=False
    )
    windward_factor = st.checkbox(
        "Apply 0.8 windward capacity factor", value=False
    )
    apply_trigger = st.checkbox(
        "Apply squeezing geometric trigger", value=True
    )
    depth_step = st.number_input(
        "Depth increment for analysis (m)", min_value=0.1, max_value=5.0,
        step=0.1, format="%.2f", value=0.25
    )
    # Compute maximum depth as the deepest bottom among layers
    max_depth = max([layer.bottom_depth for layer in layers]) if layers else 0.0

    if st.button("Compute penetration"):
        # Build spudcan object
        spud = Spudcan(
            rig_name=rig_name,
            diameter=diameter,
            area=area,
            tip_elevation=tip_elevation,
            max_preload=preload_MN * 1000.0,  # convert MN to kN
        )
        # Perform calculation (empty Meyerhof profile for manual input)
        meyerhof_profile: List[SoilDataPoint] = []
        df_results = calculate_penetration_curve(
            spud=spud,
            layers=layers,
            depth_step=depth_step,
            max_depth=max_depth,
            meyerhof=meyerhof_profile,
            use_min_cu=use_min_cu,
            apply_phi_red=apply_phi_red,
            windward_factor=windward_factor,
            apply_trigger=apply_trigger,
        )
        # Interpolate penetration at preload using controlling capacity
        ctrl_pen_depth = interpolate_penetration(df_results, preload_MN, col="REAL Capacity (MN)")
        tip_pen_depth = ctrl_pen_depth + tip_elevation
        st.subheader("Results")
        st.write(f"Controlling penetration depth: **{ctrl_pen_depth:.2f} m**")
        st.write(f"Tip penetration depth: **{tip_pen_depth:.2f} m**")
        # Plot controlling curve
        fig = go.Figure()
        # Controlling capacity curve
        fig.add_trace(
            go.Scatter(
                x=df_results["REAL Capacity (MN)"],
                y=df_results["Depth (m)"],
                mode="lines",
                name="Controlling Capacity",
                line=dict(color="blue", width=3),
            )
        )
        # Vertical preload line
        fig.add_trace(
            go.Scatter(
                x=[preload_MN, preload_MN],
                y=[0, max_depth],
                mode="lines",
                name="Max Preload",
                line=dict(color="red", width=2, dash="dash"),
            )
        )
        fig.update_layout(
            xaxis_title="Leg Load (MN)",
            yaxis_title="Penetration depth (m)",
            yaxis_autorange="reversed",
            width=700,
            height=500,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1.0,
            ),
            margin=dict(l=50, r=50, t=50, b=50),
            template="plotly_white",
        )
        st.plotly_chart(fig, use_container_width=False)
        # Show table of results
        with st.expander("Show detailed results table"):
            st.dataframe(df_results)
            # Provide CSV download
            csv_data = df_results.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="Download results as CSV",
                data=csv_data,
                file_name=f"{rig_name}_penetration_results.csv",
                mime="text/csv",
            )


if __name__ == "__main__":
    main()