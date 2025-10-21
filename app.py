"""
Streamlit application for Leg Penetration Analysis (LPA) based on the
V50 VBA macro provided by the user.  This app allows users to upload
an Excel workbook containing two sheets—`Details‑SPC` and
`Soil Information`—or to manually input data directly via the web UI.

The calculations performed mirror the logic of the VBA macro:
  * A per‑rig penetration curve is computed for a range of depths
    using undrained (clay) and drained (sand) bearing capacity
    formulas from the SNAME guidelines.
  * Squeezing and punch‑through modes are assessed.  An optional
    geometric trigger can be enforced for squeezing.
  * Backflow stability is checked using a Meyerhof N chart if
    available.
  * The engineer can choose conservative or less conservative
    settings: using the minimum of point and B/2‑average shear
    strengths (Su) for clay, applying a 5° reduction to φ for sand,
    multiplying the resulting envelopes by 0.8 for windward legs, and
    enabling/disabling the squeezing trigger.

Usage:
  Run this script with `streamlit run app.py`.  The UI presents
  options to upload a spreadsheet or to enter data manually.  After
  selecting settings and providing input, press the "Compute" button
  to view penetration curves, capacity tables, and summary results.
"""

import io
import math
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go


@dataclass
class SoilDataPoint:
    depth: float
    value: float


@dataclass
class SoilLayer:
    name: str
    top_depth: float
    bottom_depth: float
    soil_type: str  # "clay", "sand", "silt", or "unknown"
    unit_weight_profile: List[SoilDataPoint]
    su_profile: List[SoilDataPoint]
    phi_profile: List[SoilDataPoint]


@dataclass
class Spudcan:
    rig_name: str
    diameter: float
    area: float
    tip_elevation: float
    max_preload: float  # kN


# Default calculation parameters analogous to the VBA flags.
DEFAULT_USE_MIN_CU_POINT_AVG = True
DEFAULT_APPLY_PHI_REDUCTION = False
DEFAULT_APPLY_WINDWARD_CAPACITY_FACTOR = False
DEFAULT_APPLY_SQUEEZING_TRIGGER = True


# Helper functions for interpolation and averaging
def interpolate_profile(target_depth: float, profile: List[SoilDataPoint]) -> float:
    """Interpolate a value from a depth‑value profile using linear interpolation.

    If the depth is shallower than the first point, the first value
    is returned.  If deeper than the last point, the last value is
    returned.
    """
    if not profile:
        return 0.0
    # If below first point
    if target_depth <= profile[0].depth:
        return profile[0].value
    # Find interval for interpolation
    for i in range(1, len(profile)):
        if target_depth <= profile[i].depth:
            d1, v1 = profile[i - 1].depth, profile[i - 1].value
            d2, v2 = profile[i].depth, profile[i].value
            if d2 == d1:
                return v1
            return v1 + (target_depth - d1) * (v2 - v1) / (d2 - d1)
    # Beyond last point
    return profile[-1].value


def avg_over_window_b2(start_depth: float, B: float, layers: List[SoilLayer], prop: str) -> float:
    """Compute the average property value from start_depth to start_depth + B/2.

    This function integrates (via a simple stepwise sum) the soil
    property (Su or gamma) across the relevant layers.  It mirrors
    the VBA `AvgOverWindow_B2` approach, using an incremental step of
    0.2 m to approximate the integral.
    """
    step = 0.2
    end_depth = start_depth + B / 2.0
    if end_depth <= start_depth:
        return 0.0
    # Sum the property across layers by stepping through depth
    total = 0.0
    count = 0
    d = start_depth
    while d <= end_depth:
        idx = get_layer_index_at_depth(d, layers)
        layer = layers[idx]
        if prop == "su":
            val = interpolate_profile(d, layer.su_profile)
        else:
            # gamma'
            val = interpolate_profile(d, layer.unit_weight_profile)
        total += val
        count += 1
        d += step
    if count == 0:
        return 0.0
    return total / count


def get_overburden_pressure(depth: float, layers: List[SoilLayer]) -> float:
    """Compute effective overburden pressure (sum of gamma' × thickness) down to a depth.

    This is the submerged unit weight times thickness of each layer up
    to the specified depth.  It mirrors the VBA GetOverburdenPressure.
    """
    total = 0.0
    for layer in layers:
        if depth > layer.top_depth:
            if depth < layer.bottom_depth:
                thick = depth - layer.top_depth
            else:
                thick = layer.bottom_depth - layer.top_depth
            if thick > 0:
                # Average gamma across this thickness at midpoints
                midpoint = layer.top_depth + 0.5 * thick
                gamma_prime = interpolate_profile(midpoint, layer.unit_weight_profile)
                total += gamma_prime * thick
    return total


def get_layer_index_at_depth(depth: float, layers: List[SoilLayer]) -> int:
    """Return the index of the soil layer at a given depth.

    If the depth is deeper than all layer bottoms, the last index is
    returned.
    """
    for i, layer in enumerate(layers):
        if depth >= layer.top_depth and depth < layer.bottom_depth:
            return i
    return len(layers) - 1


def backflow_check(depth: float, B: float, layers: List[SoilLayer], meyerhof: List[SoilDataPoint]) -> bool:
    """Check for backflow using Meyerhof stability N chart.

    Backflow is deemed to occur if d > (N × cu_avg) / gamma_avg,
    where N is interpolated from the Meyerhof chart (function of
    d/B), cu_avg is the B/2‑average undrained shear strength above the
    layer, and gamma_avg is the B/2‑average submerged unit weight.
    If no Meyerhof data is provided or if the current layer is not
    clay/silt, False is returned.
    """
    if depth <= 0 or B <= 0 or not meyerhof:
        return False
    idx = get_layer_index_at_depth(depth, layers)
    soil_type = layers[idx].soil_type
    # Only clay or silt triggers backflow check
    if soil_type not in ("clay", "silt"):
        return False
    # Interpolate N from Meyerhof data using dimensionless ratio depth/B
    n_ratio = depth / B
    N_val = interpolate_profile(n_ratio, meyerhof)
    cu_avg = avg_over_window_b2(depth, B, layers, "su")
    gamma_avg = avg_over_window_b2(depth, B, layers, "gamma")
    if gamma_avg <= 0:
        return False
    return depth > (N_val * cu_avg) / gamma_avg


def calculate_clay_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool,
                            use_min_cu: bool) -> float:
    """Calculate idle (general shear) bearing capacity in clay.

    Implements the clay capacity calculation with options for
    conservative Su selection.  The result (kN) is returned.
    """
    if spud.diameter <= 0:
        return 0.0
    idx = get_layer_index_at_depth(depth, layers)
    # Point Su and average Su over B/2
    cu_point = interpolate_profile(depth, layers[idx].su_profile)
    cu_avg = avg_over_window_b2(depth, spud.diameter, layers, "su")
    cu_eff: float
    if use_min_cu:
        if cu_point > 0 and cu_avg > 0:
            cu_eff = min(cu_point, cu_avg)
        else:
            cu_eff = max(cu_point, cu_avg)
    else:
        cu_eff = cu_avg
    if cu_eff <= 0:
        return 0.0
    # Shape factors for circular footing
    Nc = 5.14
    sc = 1.2
    # Depth factor dc per SNAME
    d_over_B = depth / spud.diameter if spud.diameter > 0 else 0
    if d_over_B <= 1:
        dc = 1.0 + 0.4 * d_over_B
    else:
        dc = 1.0 + 0.4 * math.atan(d_over_B)
    # Effective overburden
    po_prime = 0.0 if has_backflow else get_overburden_pressure(depth, layers)
    Fv = (cu_eff * Nc * sc * dc + po_prime) * spud.area
    return Fv


def calculate_sand_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer],
                            apply_phi_red: bool) -> float:
    """Calculate idle (general shear) bearing capacity in sand.

    Returns total vertical force in kN.  Applies a 5° reduction to φ if
    requested.
    """
    if spud.diameter <= 0:
        return 0.0
    idx = get_layer_index_at_depth(depth, layers)
    gamma_prime = get_average_property_value(depth, layers, prop="gamma")
    phi = interpolate_profile(depth, layers[idx].phi_profile)
    if apply_phi_red:
        phi = phi - 5.0
    phi = max(phi, 0.0)
    po_prime = get_overburden_pressure(depth, layers)
    if phi <= 0.0:
        return po_prime * spud.area
    phi_rad = math.radians(phi)
    Nq = math.exp(math.pi * math.tan(phi_rad)) * (math.tan(math.pi / 4.0 + phi_rad / 2.0) ** 2)
    Ngamma = 2.0 * (Nq + 1.0) * math.tan(phi_rad)
    d_over_B = depth / spud.diameter if spud.diameter > 0 else 0
    dq = 1.0 + 2.0 * math.tan(phi_rad) * (1.0 - math.sin(phi_rad)) ** 2 * d_over_B
    sq = 1.0 + math.tan(phi_rad)
    sgamma = 0.6
    # Depth factor for Ngamma usually equals 1 for B < 15 m (simplified)
    Fv = (0.5 * gamma_prime * spud.diameter * Ngamma * sgamma + po_prime * Nq * sq * dq) * spud.area
    return Fv


def calculate_squeezing_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool,
                                 apply_trigger: bool) -> Optional[float]:
    """Calculate squeezing capacity for clay‑over‑clay sequences.

    Returns the vertical capacity (kN) or None if squeezing does not
    apply.  The geometric trigger B >= 3.45 * T * (1 + 1.1 * D/B) can
    be enforced by the apply_trigger flag.
    """
    idx = get_layer_index_at_depth(depth, layers)
    if idx + 1 >= len(layers):
        return None
    topL = layers[idx]
    botL = layers[idx + 1]
    # Only consider if top layer is clay
    if topL.soil_type != "clay":
        return None
    B = spud.diameter
    T = topL.bottom_depth - depth
    if T <= 0 or B <= 0:
        return None
    cu_t = avg_over_window_b2(depth, B, layers, "su")
    cu_b = avg_over_window_b2(topL.bottom_depth, B, layers, "su")
    if cu_t <= 0 or cu_b <= 0:
        return None
    # Lower layer must be at least 1.5 times stronger
    if cu_b <= 1.5 * cu_t:
        return None
    # Geometric trigger
    if apply_trigger:
        if B < 3.45 * T * (1.0 + 1.1 * (depth / B)):
            return None
    # Overburden pressure at current depth
    po_prime = 0.0 if has_backflow else get_overburden_pressure(depth, layers)
    # SNAME recommended constants a = 5.0, b = 0.33, c = 1.2
    a = 5.0
    b = 0.33
    c = 1.2
    Fv = spud.area * ((a + b * (B / T) + c * (depth / B)) * cu_t + po_prime)
    return Fv


def calculate_punch_through_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool) -> Optional[float]:
    """Calculate punch‑through capacity for clay/clay or sand/clay layering.

    Returns the vertical capacity (kN) or None if punch‑through does not
    apply.
    """
    idx = get_layer_index_at_depth(depth, layers)
    if idx + 1 >= len(layers):
        return None
    topL = layers[idx]
    botL = layers[idx + 1]
    B = spud.diameter
    H = topL.bottom_depth - depth
    if H <= 0 or B <= 0:
        return None
    # Clay over clay
    if topL.soil_type == "clay" and botL.soil_type == "clay":
        cu_t = avg_over_window_b2(depth, B, layers, "su")
        cu_b = avg_over_window_b2(topL.bottom_depth, B, layers, "su")
        if cu_t <= 0 or cu_b <= 0 or cu_t <= cu_b:
            return None
        po_p = 0.0 if has_backflow else get_overburden_pressure(depth + H, layers)
        d_over_B = (depth + H) / B
        dc = 1.0 + 0.2 * d_over_B
        Fv = spud.area * (3.0 * (H / B) * cu_t + 5.14 * 1.2 * dc * cu_b + po_p)
        return Fv
    # Sand over clay
    if topL.soil_type == "sand" and botL.soil_type == "clay":
        Fv_b = calculate_clay_capacity(spud, depth + H, layers, False, True)
        # gamma at midpoint of sand layer
        gamma_s = get_average_property_value(topL.bottom_depth, layers, "gamma", depth)
        po_s = get_overburden_pressure(depth, layers)
        cu_c = avg_over_window_b2(topL.bottom_depth, B, layers, "su")
        # Equivalent Ks * tan(phi) factor: 3*cu_c/(B*gamma_s)
        KsTanPhi = 0.0
        if B > 0 and gamma_s > 0:
            KsTanPhi = (3.0 * cu_c) / (B * gamma_s)
        Fv = Fv_b - spud.area * H * gamma_s + 2.0 * (H / B) * (H * gamma_s + 2.0 * po_s) * KsTanPhi * spud.area
        return Fv
    return None


def get_average_property_value(end_depth: float, layers: List[SoilLayer], prop: str,
                               start_depth: float = 0.0) -> float:
    """Average a soil property (su or gamma) between start_depth and end_depth.

    If end_depth <= start_depth, returns the value at end_depth.
    Mirrors the VBA GetAveragePropertyValue with a 0.2 m step.
    """
    if end_depth <= start_depth:
        idx = get_layer_index_at_depth(end_depth, layers)
        layer = layers[idx]
        if prop == "su":
            return interpolate_profile(end_depth, layer.su_profile)
        else:
            return interpolate_profile(end_depth, layer.unit_weight_profile)
    step = 0.2
    total = 0.0
    count = 0
    d = start_depth
    while d <= end_depth:
        idx = get_layer_index_at_depth(d, layers)
        layer = layers[idx]
        if prop == "su":
            val = interpolate_profile(d, layer.su_profile)
        else:
            val = interpolate_profile(d, layer.unit_weight_profile)
        total += val
        count += 1
        d += step
    if count == 0:
        return 0.0
    return total / count


def calculate_penetration_curve(spud: Spudcan, layers: List[SoilLayer], depth_step: float,
                                max_depth: float, meyerhof: List[SoilDataPoint],
                                use_min_cu: bool, apply_phi_red: bool,
                                windward_factor: bool, apply_trigger: bool) -> pd.DataFrame:
    """Compute the penetration curve for a given spudcan and soil profile.

    Returns a DataFrame with columns: Depth, Idle_Clay (MN), Idle_Sand (MN),
    Real_Clay (MN), Real_Sand (MN), Backflow (bool), Clay_Mode, Sand_Mode,
    Squeeze_Cap (MN), PT_Cap (MN).  All capacities are returned in MN.
    """
    depths = np.arange(0, max_depth + depth_step, depth_step)
    # Prepare arrays for results
    results = {
        "Depth (m)": [],
        "Idle Clay (MN)": [],
        "Idle Sand (MN)": [],
        "REAL Clay (MN)": [],
        "REAL Sand (MN)": [],
        "Backflow": [],
        "Clay Mode": [],
        "Sand Mode": [],
        "Squeeze Cap (MN)": [],
        "PT Cap (MN)": []
    }
    for d in depths:
        idx = get_layer_index_at_depth(d, layers)
        has_backflow = backflow_check(d, spud.diameter, layers, meyerhof)
        idle_clay_kN = calculate_clay_capacity(spud, d, layers, has_backflow, use_min_cu)
        idle_sand_kN = calculate_sand_capacity(spud, d, layers, apply_phi_red)
        squeeze_cap_kN = calculate_squeezing_capacity(spud, d, layers, has_backflow, apply_trigger)
        pt_cap_kN = calculate_punch_through_capacity(spud, d, layers, has_backflow)
        # Determine real clay and real sand capacities
        # Start with idle; then apply special capacities if smaller
        real_clay_kN = idle_clay_kN
        clay_mode = "General"
        real_sand_kN = idle_sand_kN
        sand_mode = "General"
        # Squeeze or PT only apply if result is not None
        if squeeze_cap_kN is not None:
            if squeeze_cap_kN < real_clay_kN or real_clay_kN == 0:
                real_clay_kN = squeeze_cap_kN
                clay_mode = "Squeezing"
        if pt_cap_kN is not None:
            # If top layer is sand, PT affects sand; else clay
            if layers[idx].soil_type == "sand":
                if pt_cap_kN < real_sand_kN or real_sand_kN == 0:
                    real_sand_kN = pt_cap_kN
                    sand_mode = "Punch‑Through"
            else:
                if pt_cap_kN < real_clay_kN or real_clay_kN == 0:
                    real_clay_kN = pt_cap_kN
                    clay_mode = "Punch‑Through"
        # Apply windward factor if requested
        if windward_factor:
            real_clay_kN *= 0.8
            real_sand_kN *= 0.8
        # Convert to MN
        idle_clay_MN = idle_clay_kN / 1000.0 if idle_clay_kN > 0 else None
        idle_sand_MN = idle_sand_kN / 1000.0 if idle_sand_kN > 0 else None
        real_clay_MN = real_clay_kN / 1000.0 if real_clay_kN > 0 else None
        real_sand_MN = real_sand_kN / 1000.0 if real_sand_kN > 0 else None
        results["Depth (m)"].append(d)
        results["Idle Clay (MN)"].append(idle_clay_MN)
        results["Idle Sand (MN)"].append(idle_sand_MN)
        results["REAL Clay (MN)"].append(real_clay_MN)
        results["REAL Sand (MN)"].append(real_sand_MN)
        results["Backflow"].append(has_backflow)
        results["Clay Mode"].append(clay_mode)
        results["Sand Mode"].append(sand_mode)
        results["Squeeze Cap (MN)"].append(squeeze_cap_kN / 1000.0 if squeeze_cap_kN is not None else None)
        results["PT Cap (MN)"].append(pt_cap_kN / 1000.0 if pt_cap_kN is not None else None)
    df = pd.DataFrame(results)
    return df


def interpolate_penetration(df: pd.DataFrame, preload_MN: float, col: str) -> float:
    """Find penetration (depth) at which the real clay/sand capacity equals the preload.

    Uses linear interpolation between data points in the DataFrame.  If
    preload exceeds the maximum capacity, the maximum depth is returned.
    """
    depths = df["Depth (m)"].values
    capacities = df[col].values
    # If all capacities are None or NaN, return deepest depth
    if np.all(pd.isnull(capacities)):
        return depths[-1]
    # Convert None to NaN for interpolation
    capacities = np.array([c if c is not None else np.nan for c in capacities], dtype=float)
    # Loop through segments
    for i in range(1, len(capacities)):
        x1, x2 = capacities[i - 1], capacities[i]
        if np.isnan(x1) or np.isnan(x2):
            continue
        if (x1 <= preload_MN <= x2) or (x2 <= preload_MN <= x1):
            d1, d2 = depths[i - 1], depths[i]
            if x2 == x1:
                return d2
            # Linear interpolation on capacity vs depth
            return d1 + (preload_MN - x1) * (d2 - d1) / (x2 - x1)
    # If preload is outside range, return deepest depth
    return depths[-1]


# Data parsing from Excel
def parse_excel_details_spc(df: pd.DataFrame) -> List[Spudcan]:
    """Parse spudcan details from the Details‑SPC sheet.

    The DataFrame `df` is assumed to have at least five columns:
    Column 0: RigName, 1: Diameter (m), 2: Area (m2), 3: Tip Elevation (m), 4: Max Preload (MN).
    Columns 11 and 12 (L and M) may contain Meyerhof d/B and N values.
    Returns a list of Spudcan objects and the Meyerhof profile as a list
    of SoilDataPoint (depth ratio vs N).
    """
    spuds = []
    meyerhof_profile: List[SoilDataPoint] = []
    # Determine number of rows; assume header row first or no header; ignore empty lines
    for _, row in df.iterrows():
        rig = str(row.get(0, "")).strip()
        if rig == "" or pd.isna(rig):
            continue
        try:
            diameter = float(row.get(1, 0))
        except Exception:
            diameter = 0.0
        try:
            area = float(row.get(2, 0))
        except Exception:
            area = 0.0
        try:
            tip_elev = float(row.get(3, 0))
        except Exception:
            tip_elev = 0.0
        try:
            preload_MN = float(row.get(4, 0))
        except Exception:
            preload_MN = 0.0
        spuds.append(Spudcan(rig_name=rig, diameter=diameter, area=area,
                              tip_elevation=tip_elev, max_preload=preload_MN * 1000.0))
    # Parse Meyerhof N chart if provided (columns 11 and 12, 0‑indexed 10 and 11)
    try:
        depth_ratio_col = df.columns[10]
        n_col = df.columns[11]
        for _, row in df.iterrows():
            depth_ratio = row.get(depth_ratio_col, None)
            n_value = row.get(n_col, None)
            if depth_ratio is not None and n_value is not None and pd.notna(depth_ratio) and pd.notna(n_value):
                try:
                    dr = float(depth_ratio)
                    nv = float(n_value)
                    meyerhof_profile.append(SoilDataPoint(depth=dr, value=nv))
                except Exception:
                    pass
    except Exception:
        # If columns not present, ignore
        pass
    # Sort Meyerhof profile by depth ratio
    meyerhof_profile.sort(key=lambda p: p.depth)
    return spuds, meyerhof_profile


def parse_soil_information(df: pd.DataFrame) -> Tuple[List[SoilLayer], float]:
    """Parse the Soil Information sheet into a list of SoilLayer objects.

    The expected format mirrors the user’s Excel template:
      * Each soil layer has a row with: A: name, B: top depth, C: bottom depth,
        E: description (text containing clay/silt/sand), F: depth column for gamma,
        G: gamma values, H: depth for phi, I: phi values, J: depth for Su, K: Su values,
        L: total depth (present in cell L3).

    This parser reads each layer header and collects the depth‑value
    profiles under columns F–G (gamma), H–I (phi), and J–K (Su)
    between the layer’s top and bottom depths.  Zero or missing Su
    values are filtered out.
    Returns the list of SoilLayer objects and the maximum analysis
    depth (total depth) if provided.
    """
    layers: List[SoilLayer] = []
    # Total depth may be stored in L3 (0‑based row index 2, column 11)
    max_depth = 0.0
    try:
        max_depth = float(df.iloc[2, 11])
    except Exception:
        max_depth = 0.0
    # Identify layer header rows: those with non‑empty value in column 0 (A)
    layer_indices = []
    for i in range(len(df)):
        val = df.iloc[i, 0]
        if not pd.isna(val) and str(val).strip() != "":
            layer_indices.append(i)
    # Build each layer from header row to next header or end
    for idx, row_index in enumerate(layer_indices):
        name = str(df.iloc[row_index, 0])
        try:
            top_depth = float(df.iloc[row_index, 1])
        except Exception:
            top_depth = 0.0
        try:
            bottom_depth = float(df.iloc[row_index, 2])
        except Exception:
            bottom_depth = top_depth  # default if missing
        soil_desc = str(df.iloc[row_index, 4]).lower() if not pd.isna(df.iloc[row_index, 4]) else ""
        if "clay" in soil_desc:
            soil_type = "clay"
        elif "silt" in soil_desc:
            soil_type = "silt"
        elif "sand" in soil_desc:
            soil_type = "sand"
        else:
            soil_type = "unknown"
        # Collect profile rows within this layer (from top to bottom)
        # Depth columns: F (index 5) with gamma in G (index 6), H (7) phi in I (8), J (9) Su in K (10)
        unit_profile: List[SoilDataPoint] = []
        phi_profile: List[SoilDataPoint] = []
        su_profile: List[SoilDataPoint] = []
        # Determine end row for this layer (start of next header or end of sheet)
        end_row = layer_indices[idx + 1] if idx + 1 < len(layer_indices) else len(df)
        for j in range(row_index + 1, end_row):
            try:
                d_gamma = df.iloc[j, 5]
                v_gamma = df.iloc[j, 6]
                if pd.notna(d_gamma) and pd.notna(v_gamma):
                    d_gamma = float(d_gamma)
                    v_gamma = float(v_gamma)
                    if top_depth <= d_gamma < bottom_depth:
                        unit_profile.append(SoilDataPoint(depth=d_gamma, value=v_gamma))
            except Exception:
                pass
            try:
                d_phi = df.iloc[j, 7]
                v_phi = df.iloc[j, 8]
                if pd.notna(d_phi) and pd.notna(v_phi):
                    d_phi = float(d_phi)
                    v_phi = float(v_phi)
                    if top_depth <= d_phi < bottom_depth:
                        phi_profile.append(SoilDataPoint(depth=d_phi, value=v_phi))
            except Exception:
                pass
            try:
                d_su = df.iloc[j, 9]
                v_su = df.iloc[j, 10]
                if pd.notna(d_su) and pd.notna(v_su):
                    d_su = float(d_su)
                    v_su = float(v_su)
                    if top_depth <= d_su < bottom_depth and v_su > 0:
                        su_profile.append(SoilDataPoint(depth=d_su, value=v_su))
            except Exception:
                pass
        # Sort profiles by depth
        unit_profile.sort(key=lambda p: p.depth)
        phi_profile.sort(key=lambda p: p.depth)
        su_profile.sort(key=lambda p: p.depth)
        layers.append(SoilLayer(name=name, top_depth=top_depth, bottom_depth=bottom_depth,
                                soil_type=soil_type, unit_weight_profile=unit_profile,
                                su_profile=su_profile, phi_profile=phi_profile))
    return layers, max_depth


def display_penetration_curve(df: pd.DataFrame, spud: Spudcan):
    """Display the penetration curve using Plotly for interactive features."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["Idle Clay (MN)"], y=df["Depth (m)"], mode="lines", name="Idle Clay",
                             line=dict(color="black", width=1.2)))
    fig.add_trace(go.Scatter(x=df["Idle Sand (MN)"], y=df["Depth (m)"], mode="lines", name="Idle Sand",
                             line=dict(color="black", width=1.2, dash="dash")))
    fig.add_trace(go.Scatter(x=df["REAL Clay (MN)"], y=df["Depth (m)"], mode="lines", name="REAL Clay",
                             line=dict(color="blue", width=2.4)))
    fig.add_trace(go.Scatter(x=df["REAL Sand (MN)"], y=df["Depth (m)"], mode="lines", name="REAL Sand",
                             line=dict(color="blue", width=2.0, dash="dash")))
    # Preload line
    preload_x = [spud.max_preload / 1000.0, spud.max_preload / 1000.0]
    preload_y = [0, df["Depth (m)"].max()]
    fig.add_trace(go.Scatter(x=preload_x, y=preload_y, mode="lines", name="Max Preload",
                             line=dict(color="red", dash="dash")))
    fig.update_yaxes(autorange="reversed", title="Penetration (m)")
    fig.update_xaxes(title="Leg Load (MN)")
    fig.update_layout(title=f"Leg Penetration Analysis: {spud.rig_name}", legend=dict(orientation="h"))
    st.plotly_chart(fig, use_container_width=True)


def main():
    st.set_page_config(page_title="Leg Penetration Analysis (V50)", layout="wide")
    st.title("Leg Penetration Analysis (V50)")
    st.write("This application implements the V50 leg penetration analysis described in the provided VBA macro.")
    # Sidebar settings
    st.sidebar.header("Settings")
    use_min_cu = st.sidebar.checkbox("Use min(point, B/2-average) for Su (conservative)",
                                     value=DEFAULT_USE_MIN_CU_POINT_AVG)
    apply_phi_red = st.sidebar.checkbox("Apply 5° reduction to φ (sand)",
                                        value=DEFAULT_APPLY_PHI_REDUCTION)
    windward_factor = st.sidebar.checkbox("Multiply REAL capacities by 0.8 (windward leg)",
                                         value=DEFAULT_APPLY_WINDWARD_CAPACITY_FACTOR)
    apply_trigger = st.sidebar.checkbox("Apply geometric trigger for squeezing", value=DEFAULT_APPLY_SQUEEZING_TRIGGER)

    st.sidebar.header("Input Method")
    input_method = st.sidebar.radio("Select data input method", ["Upload Excel", "Manual Input"], index=0)

    spud_list: List[Spudcan] = []
    soil_layers: List[SoilLayer] = []
    meyerhof_data: List[SoilDataPoint] = []
    max_depth = 0.0

    if input_method == "Upload Excel":
        uploaded_file = st.file_uploader("Upload Excel file (.xlsx) with Details-SPC and Soil Information sheets", type=["xlsx"])
        if uploaded_file is not None:
            try:
                xls = pd.ExcelFile(uploaded_file)
                # Read sheets without headers (raw cells)
                df_details = pd.read_excel(xls, sheet_name="Details-SPC", header=None)
                df_soil = pd.read_excel(xls, sheet_name="Soil Information", header=None)
                spud_list, meyerhof_data = parse_excel_details_spc(df_details)
                soil_layers, max_depth = parse_soil_information(df_soil)
                st.success(f"Loaded {len(spud_list)} rig(s) and {len(soil_layers)} soil layer(s) from file.")
            except Exception as e:
                st.error(f"Error reading Excel file: {e}")
    else:
        st.subheader("Manual Input")
        # Manual rig input
        rig_name = st.text_input("Rig Name", value="Rig1")
        diameter = st.number_input("Spudcan Diameter (m)", min_value=0.0, value=20.0, step=0.1)
        area = st.number_input("Spudcan Area (m²)", min_value=0.0, value=314.0, step=0.1)
        tip_elevation = st.number_input("Tip Elevation above base (m)", value=0.0, step=0.1)
        max_preload = st.number_input("Max Preload (MN)", min_value=0.0, value=80.0, step=1.0)
        spud_list = [Spudcan(rig_name=rig_name, diameter=diameter, area=area,
                             tip_elevation=tip_elevation, max_preload=max_preload * 1000.0)]
        # Manual soil input (simplified): number of layers
        st.markdown("### Soil Layers (manual input)")
        num_layers = st.number_input("Number of layers", min_value=1, max_value=10, value=2, step=1)
        soil_layers = []
        for i in range(int(num_layers)):
            st.markdown(f"**Layer {i + 1}**")
            layer_name = st.text_input(f"Layer {i + 1} name", value=f"Layer {i + 1}")
            top_depth = st.number_input(f"Layer {i + 1} top depth (m)", value=float(i * 5), step=0.1)
            bottom_depth = st.number_input(f"Layer {i + 1} bottom depth (m)", value=float((i + 1) * 5), step=0.1)
            soil_type = st.selectbox(f"Layer {i + 1} soil type", ["clay", "sand", "silt", "unknown"], index=0)
            # Uniform properties across layer for manual input
            gamma_prime = st.number_input(f"Layer {i + 1} γ' (kN/m³)", value=10.0, step=0.1)
            phi_val = st.number_input(f"Layer {i + 1} φ (degrees)", value=0.0, step=0.1)
            su_val = st.number_input(f"Layer {i + 1} Su (kPa)", value=0.0, step=0.1)
            # Build profiles with a single data point at layer mid depth
            mid_depth = (top_depth + bottom_depth) / 2.0
            unit_profile = [SoilDataPoint(depth=mid_depth, value=gamma_prime)]
            phi_profile = [SoilDataPoint(depth=mid_depth, value=phi_val)]
            su_profile = []
            if su_val > 0:
                su_profile = [SoilDataPoint(depth=mid_depth, value=su_val)]
            soil_layers.append(SoilLayer(name=layer_name, top_depth=top_depth,
                                         bottom_depth=bottom_depth, soil_type=soil_type,
                                         unit_weight_profile=unit_profile,
                                         su_profile=su_profile, phi_profile=phi_profile))
        # Prompt for maximum depth
        max_depth = st.number_input("Maximum analysis depth (m)", min_value=0.0, value=float(num_layers * 5), step=0.5)
        meyerhof_data = []  # No backflow chart in manual mode

    if spud_list and soil_layers:
        # Depth resolution and compute button
        depth_step = st.sidebar.number_input("Depth increment (m)", min_value=0.1, max_value=5.0, value=0.25, step=0.05)
        # If max_depth not given from file, use bottom of last layer
        if max_depth <= 0:
            max_depth = soil_layers[-1].bottom_depth
        if st.sidebar.button("Compute"):
            for spud in spud_list:
                st.header(f"Results for {spud.rig_name}")
                df_results = calculate_penetration_curve(spud, soil_layers, depth_step,
                                                         max_depth, meyerhof_data, use_min_cu,
                                                         apply_phi_red, windward_factor, apply_trigger)
                # Calculate penetration depths at preload (analysis depth)
                preload_MN = spud.max_preload / 1000.0
                pen_clay = interpolate_penetration(df_results, preload_MN, "REAL Clay (MN)")
                pen_sand = interpolate_penetration(df_results, preload_MN, "REAL Sand (MN)")
                tip_pen_clay = pen_clay + spud.tip_elevation
                tip_pen_sand = pen_sand + spud.tip_elevation
                st.write(f"**Predicted tip penetration (Clay, REAL):** {tip_pen_clay:.2f} m")
                st.write(f"**Predicted tip penetration (Sand, REAL):** {tip_pen_sand:.2f} m")
                range_low = min(tip_pen_clay, tip_pen_sand)
                range_high = max(tip_pen_clay, tip_pen_sand)
                st.write(f"**Penetration range (tip):** {range_low:.2f} – {range_high:.2f} m")
                # Display table
                st.dataframe(df_results)
                # Plot chart
                display_penetration_curve(df_results, spud)
    else:
        st.info("Please provide spud and soil information.")


if __name__ == "__main__":
    main()