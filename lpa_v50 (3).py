"""
lpa_v50.py
===============

This module contains the core computation logic for the Leg Penetration
Analysis (LPA) as requested by the user.  It is a pure‑Python port of
the VBA macro logic provided previously, with additional refinements to
handle cases where either undrained shear strength (Su) or drained
friction angle (φ) is absent.  The functions herein do not depend on
Streamlit or any UI, making them reusable for both web and desktop
applications.

Key Features
------------

* **Data classes** for Spudcan, SoilLayer and SoilDataPoint to hold
  input data cleanly.
* **Interpolation and averaging** routines that mirror the VBA
  behaviour (linear interpolation between depth‑value pairs and
  average over a B/2 window for Su and γ').
* **Bearing capacity calculations** for general shear in clay and sand
  following SNAME §6.2.2–6.2.5, including options to apply
  conservative Su selection (min of point and B/2 average), reduce φ
  by 5°, and apply a windward factor of 0.8.
* **Special failure modes** (squeezing and punch‑through) for
  layered soils in accordance with SNAME §6.2.6, with the
  recommended geometric trigger for squeezing and simplified
  punch‑through formulas.
* **Backflow check** using a Meyerhof N chart (optional).
* **Penetration curve generation** that returns a pandas DataFrame
  containing per‑depth idle capacities, real capacities (with special
  modes applied), and a new **controlling capacity** column.  The
  controlling capacity is the minimum of the available real clay and
  real sand capacities at that depth.  If either Su or φ is not
  provided at a depth, only the corresponding capacity is computed,
  avoiding unrealistic upper/lower bound comparisons.

Usage
-----

Import this module and call `calculate_penetration_curve` with a
`Spudcan` object, a list of `SoilLayer` objects, calculation
parameters and (optionally) a Meyerhof chart.  The returned DataFrame
includes the depths and capacities in mega‑newtons (MN).  See
`app_ui.py` for an example of how to integrate this module into a
Streamlit interface.

"""

from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

import math
import numpy as np
import pandas as pd


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


# Default calculation flags – these mirror the VBA constants but
# remain configurable in the UI layer.
DEFAULT_USE_MIN_CU_POINT_AVG = True
DEFAULT_APPLY_PHI_REDUCTION = False
DEFAULT_APPLY_WINDWARD_CAPACITY_FACTOR = False
DEFAULT_APPLY_SQUEEZING_TRIGGER = True


def interpolate_profile(target_depth: float, profile: List[SoilDataPoint]) -> float:
    """Interpolate a value from a depth‑value profile using linear interpolation.

    If the profile is empty, returns 0.  If the depth is shallower
    than the first point, the first value is returned.  If deeper
    than the last point, the last value is returned.  Otherwise
    performs linear interpolation between the two bounding points.
    """
    if not profile:
        return 0.0
    # If shallower than first point
    if target_depth <= profile[0].depth:
        return profile[0].value
    # Walk through to find correct interval
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

    This integrates the soil property (Su or gamma') across the
    relevant layers by sampling at 0.2 m increments.  Mirrors the VBA
    `AvgOverWindow_B2` logic.
    """
    if B <= 0:
        return 0.0
    step = 0.2
    end_depth = start_depth + B / 2.0
    if end_depth <= start_depth:
        return 0.0
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
    return total / count if count > 0 else 0.0


def get_overburden_pressure(depth: float, layers: List[SoilLayer]) -> float:
    """Compute effective overburden pressure (sum of γ' × thickness) down to a depth.

    Integrates the submerged unit weight across all layers down to
    `depth`.  Mirrors the VBA `GetOverburdenPressure`.
    """
    total = 0.0
    for layer in layers:
        if depth > layer.top_depth:
            if depth < layer.bottom_depth:
                thick = depth - layer.top_depth
            else:
                thick = layer.bottom_depth - layer.top_depth
            if thick > 0:
                midpoint = layer.top_depth + 0.5 * thick
                gamma_prime = interpolate_profile(midpoint, layer.unit_weight_profile)
                total += gamma_prime * thick
    return total


def get_layer_index_at_depth(depth: float, layers: List[SoilLayer]) -> int:
    """Return the index of the soil layer at a given depth.

    If depth is deeper than the bottoms of all layers, returns the
    last index.
    """
    for i, layer in enumerate(layers):
        if depth >= layer.top_depth and depth < layer.bottom_depth:
            return i
    return len(layers) - 1


def backflow_check(depth: float, B: float, layers: List[SoilLayer], meyerhof: List[SoilDataPoint]) -> bool:
    """Check for backflow using a Meyerhof stability N chart.

    Backflow is deemed to occur if d > (N × cu_avg) / γ_avg, where N
    is interpolated from the Meyerhof chart as a function of d/B.
    Returns False if there is no Meyerhof data or if the current
    layer is not clay or silt.  Mirrors the VBA `backflow_check`.
    """
    if depth <= 0 or B <= 0 or not meyerhof:
        return False
    idx = get_layer_index_at_depth(depth, layers)
    soil_type = layers[idx].soil_type
    if soil_type not in ("clay", "silt"):
        return False
    ratio = depth / B
    N_val = interpolate_profile(ratio, meyerhof)
    cu_avg = avg_over_window_b2(depth, B, layers, "su")
    gamma_avg = avg_over_window_b2(depth, B, layers, "gamma")
    if gamma_avg <= 0:
        return False
    return depth > (N_val * cu_avg) / gamma_avg


def calculate_clay_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool,
                            use_min_cu: bool) -> float:
    """Calculate idle (general shear) bearing capacity in clay.

    Returns the total vertical force in kN.  If either the spud
    diameter is non‑positive or the effective undrained shear strength
    is zero, returns 0.  Mirrors the VBA `CalculateClayCapacity`.
    """
    if spud.diameter <= 0:
        return 0.0
    idx = get_layer_index_at_depth(depth, layers)
    cu_point = interpolate_profile(depth, layers[idx].su_profile)
    cu_avg = avg_over_window_b2(depth, spud.diameter, layers, "su")
    # Select effective Su
    if use_min_cu:
        if cu_point > 0 and cu_avg > 0:
            cu_eff = min(cu_point, cu_avg)
        else:
            cu_eff = max(cu_point, cu_avg)
    else:
        cu_eff = cu_avg
    if cu_eff <= 0:
        return 0.0
    # Shape and depth factors for circular footing
    Nc = 5.14
    sc = 1.2
    d_over_B = depth / spud.diameter if spud.diameter > 0 else 0.0
    if d_over_B <= 1.0:
        dc = 1.0 + 0.4 * d_over_B
    else:
        dc = 1.0 + 0.4 * math.atan(d_over_B)
    po_prime = 0.0 if has_backflow else get_overburden_pressure(depth, layers)
    Fv = (cu_eff * Nc * sc * dc + po_prime) * spud.area
    return Fv


def calculate_sand_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer],
                            apply_phi_red: bool) -> float:
    """Calculate idle (general shear) bearing capacity in sand.

    Returns total vertical force in kN.  If φ is zero or undefined,
    returns the overburden pressure contribution only (po' × area).
    Mirrors the VBA `CalculateSandCapacity`.
    """
    if spud.diameter <= 0:
        return 0.0
    idx = get_layer_index_at_depth(depth, layers)
    gamma_prime = get_average_property_value(depth, layers, prop="gamma")
    phi = interpolate_profile(depth, layers[idx].phi_profile)
    if apply_phi_red:
        phi -= 5.0
    phi = max(phi, 0.0)
    po_prime = get_overburden_pressure(depth, layers)
    if phi <= 0.0:
        return po_prime * spud.area
    phi_rad = math.radians(phi)
    Nq = math.exp(math.pi * math.tan(phi_rad)) * (math.tan(math.pi / 4.0 + phi_rad / 2.0) ** 2)
    Ngamma = 2.0 * (Nq + 1.0) * math.tan(phi_rad)
    d_over_B = depth / spud.diameter if spud.diameter > 0 else 0.0
    dq = 1.0 + 2.0 * math.tan(phi_rad) * ((1.0 - math.sin(phi_rad)) ** 2) * d_over_B
    sq = 1.0 + math.tan(phi_rad)
    sgamma = 0.6
    Fv = (0.5 * gamma_prime * spud.diameter * Ngamma * sgamma + po_prime * Nq * sq * dq) * spud.area
    return Fv


def calculate_squeezing_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool,
                                 apply_trigger: bool) -> Optional[float]:
    """Calculate squeezing capacity for clay‑over‑clay layering.

    Returns vertical capacity in kN or None if squeezing does not
    apply.  Implements SNAME §6.2.6.1 and respects the geometric
    trigger (if apply_trigger is True).
    """
    idx = get_layer_index_at_depth(depth, layers)
    if idx + 1 >= len(layers):
        return None
    topL = layers[idx]
    botL = layers[idx + 1]
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
    if cu_b <= 1.5 * cu_t:
        return None
    if apply_trigger:
        if B < 3.45 * T * (1.0 + 1.1 * (depth / B)):
            return None
    po_prime = 0.0 if has_backflow else get_overburden_pressure(depth, layers)
    # SNAME recommended constants a = 5.0, b = 0.33, c = 1.2
    a = 5.0
    b = 0.33
    c = 1.2
    Fv = spud.area * ((a + b * (B / T) + c * (depth / B)) * cu_t + po_prime)
    return Fv


def calculate_punch_through_capacity(spud: Spudcan, depth: float, layers: List[SoilLayer], has_backflow: bool) -> Optional[float]:
    """Calculate punch‑through capacity for layered soils.

    Supports clay over clay and sand over clay sequences.  Returns
    vertical capacity in kN or None if punch‑through does not apply.
    Implements the simplified forms of SNAME §6.2.6.2–6.2.6.3.  If
    Su or φ data are missing, punch‑through is not evaluated.
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
        # Underlying clay capacity at depth d+H (no backflow for PT)
        Fv_b = calculate_clay_capacity(spud, depth + H, layers, False, True)
        gamma_s = get_average_property_value(topL.bottom_depth, layers, prop="gamma", start_depth=depth)
        po_s = get_overburden_pressure(depth, layers)
        cu_c = avg_over_window_b2(topL.bottom_depth, B, layers, "su")
        KsTanPhi = 0.0
        if B > 0 and gamma_s > 0:
            KsTanPhi = (3.0 * cu_c) / (B * gamma_s)
        Fv = Fv_b - spud.area * H * gamma_s + 2.0 * (H / B) * (H * gamma_s + 2.0 * po_s) * KsTanPhi * spud.area
        return Fv
    return None


def get_average_property_value(end_depth: float, layers: List[SoilLayer], prop: str,
                               start_depth: float = 0.0) -> float:
    """Average a soil property (Su or γ') between start_depth and end_depth.

    Mirrors the VBA `GetAveragePropertyValue`.  If end_depth <=
    start_depth, returns the property at end_depth.
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
    return total / count if count > 0 else 0.0


def calculate_penetration_curve(spud: Spudcan, layers: List[SoilLayer], depth_step: float,
                                max_depth: float, meyerhof: List[SoilDataPoint],
                                use_min_cu: bool, apply_phi_red: bool,
                                windward_factor: bool, apply_trigger: bool) -> pd.DataFrame:
    """Compute the penetration curve for a given spudcan and soil profile.

    Returns a DataFrame with the following columns (all capacities in MN):

      * Depth (m)
      * Idle Clay (MN) – general shear capacity in clay (or None if not applicable)
      * Idle Sand (MN) – general shear capacity in sand (or None if not applicable)
      * REAL Clay (MN) – governing clay capacity including squeezing/punch‑through (or None)
      * REAL Sand (MN) – governing sand capacity including punch‑through (or None)
      * REAL Capacity (MN) – minimum of REAL Clay and REAL Sand where available
      * Governing Mode – "clay" or "sand" depending on which capacity controls
      * Backflow – True/False flag
      * Squeeze Cap (MN) – squeezing capacity if computed
      * PT Cap (MN) – punch‑through capacity if computed

    The controlling capacity column facilitates a single penetration
    prediction instead of comparing clay and sand bounds when only one
    soil parameter is available.
    """
    depths = np.arange(0.0, max_depth + depth_step, depth_step)
    results: Dict[str, List] = {
        "Depth (m)": [],
        "Idle Clay (MN)": [],
        "Idle Sand (MN)": [],
        "REAL Clay (MN)": [],
        "REAL Sand (MN)": [],
        "REAL Capacity (MN)": [],
        "Governing Mode": [],
        "Backflow": [],
        "Squeeze Cap (MN)": [],
        "PT Cap (MN)": [],
    }
    for d in depths:
        idx = get_layer_index_at_depth(d, layers)
        has_backflow = backflow_check(d, spud.diameter, layers, meyerhof)
        # Determine if Su or φ data exist at this depth
        layer = layers[idx]
        su_val = interpolate_profile(d, layer.su_profile)
        phi_val = interpolate_profile(d, layer.phi_profile)
        su_valid = bool(layer.su_profile) and su_val > 0.0
        phi_valid = bool(layer.phi_profile) and phi_val > 0.0
        # Compute idle capacities only if data valid
        idle_clay_kN = calculate_clay_capacity(spud, d, layers, has_backflow, use_min_cu) if su_valid else 0.0
        idle_sand_kN = calculate_sand_capacity(spud, d, layers, apply_phi_red) if phi_valid else 0.0
        # Special modes apply only if Su data exists (they rely on clay strength)
        squeeze_cap_kN = calculate_squeezing_capacity(spud, d, layers, has_backflow, apply_trigger) if su_valid else None
        pt_cap_kN = calculate_punch_through_capacity(spud, d, layers, has_backflow) if su_valid else None
        # Determine real clay and sand capacities
        real_clay_kN: Optional[float] = None
        real_sand_kN: Optional[float] = None
        clay_mode: Optional[str] = None
        sand_mode: Optional[str] = None
        if su_valid:
            real_clay_kN = idle_clay_kN
            clay_mode = "General"
            if squeeze_cap_kN is not None:
                if real_clay_kN <= 0 or squeeze_cap_kN < real_clay_kN:
                    real_clay_kN = squeeze_cap_kN
                    clay_mode = "Squeezing"
            if pt_cap_kN is not None and (layer.soil_type != "sand"):
                if real_clay_kN <= 0 or pt_cap_kN < real_clay_kN:
                    real_clay_kN = pt_cap_kN
                    clay_mode = "Punch‑Through"
        # Sand capacities
        if phi_valid:
            real_sand_kN = idle_sand_kN
            sand_mode = "General"
            if pt_cap_kN is not None and (layer.soil_type == "sand"):
                if real_sand_kN <= 0 or pt_cap_kN < real_sand_kN:
                    real_sand_kN = pt_cap_kN
                    sand_mode = "Punch‑Through"
        # Apply windward factor to real capacities
        if windward_factor:
            if real_clay_kN is not None:
                real_clay_kN *= 0.8
            if real_sand_kN is not None:
                real_sand_kN *= 0.8
        # Convert to MN (set None if capacity <= 0)
        idle_clay_MN = idle_clay_kN / 1000.0 if idle_clay_kN > 0 else None
        idle_sand_MN = idle_sand_kN / 1000.0 if idle_sand_kN > 0 else None
        real_clay_MN = real_clay_kN / 1000.0 if real_clay_kN and real_clay_kN > 0 else None
        real_sand_MN = real_sand_kN / 1000.0 if real_sand_kN and real_sand_kN > 0 else None
        # Determine controlling capacity and governing mode
        ctrl_MN = None
        ctrl_mode = None
        # Consider only positive capacities for control
        candidates: List[Tuple[float, str]] = []
        if real_clay_MN is not None:
            candidates.append((real_clay_MN, "clay"))
        if real_sand_MN is not None:
            candidates.append((real_sand_MN, "sand"))
        if candidates:
            ctrl_MN, ctrl_mode = min(candidates, key=lambda x: x[0])
        results["Depth (m)"].append(d)
        results["Idle Clay (MN)"].append(idle_clay_MN)
        results["Idle Sand (MN)"].append(idle_sand_MN)
        results["REAL Clay (MN)"].append(real_clay_MN)
        results["REAL Sand (MN)"].append(real_sand_MN)
        results["REAL Capacity (MN)"].append(ctrl_MN)
        results["Governing Mode"].append(ctrl_mode)
        results["Backflow"].append(has_backflow)
        results["Squeeze Cap (MN)"].append(squeeze_cap_kN / 1000.0 if squeeze_cap_kN else None)
        results["PT Cap (MN)"].append(pt_cap_kN / 1000.0 if pt_cap_kN else None)
    return pd.DataFrame(results)


def interpolate_penetration(df: pd.DataFrame, preload_MN: float, col: str) -> float:
    """Find the penetration (depth) at which a capacity column equals the preload.

    Performs linear interpolation on the specified DataFrame column
    (which may contain None/NaN values).  If all values are NaN,
    returns the deepest depth.  If preload exceeds the maximum
    capacity in the column, returns the deepest depth.  Mirrors the
    VBA `InterpolatePenetration`.
    """
    depths = df["Depth (m)"].values
    capacities = df[col].values
    # Convert None to NaN for interpolation
    capacities = np.array([c if c is not None else np.nan for c in capacities], dtype=float)
    if np.all(np.isnan(capacities)):
        return depths[-1]
    for i in range(1, len(capacities)):
        x1, x2 = capacities[i - 1], capacities[i]
        if np.isnan(x1) or np.isnan(x2):
            continue
        if (x1 <= preload_MN <= x2) or (x2 <= preload_MN <= x1):
            d1, d2 = depths[i - 1], depths[i]
            if x2 == x1:
                return d2
            return d1 + (preload_MN - x1) * (d2 - d1) / (x2 - x1)
    return depths[-1]


def parse_excel_details_spc(df: pd.DataFrame) -> Tuple[List[Spudcan], List[SoilDataPoint]]:
    """Parse spudcan details and Meyerhof N chart from a raw DataFrame.

    The DataFrame `df` corresponds to the `Details‑SPC` sheet with
    no header.  Columns 0–4 contain rig name, diameter (m), area
    (m²), tip elevation (m) and max preload (MN).  Columns 10 and 11
    may contain depth ratio (d/B) and N values for the Meyerhof
    backflow chart.  Returns a list of `Spudcan` objects (max
    preload converted to kN) and a sorted Meyerhof profile.
    """
    spuds: List[Spudcan] = []
    meyerhof_profile: List[SoilDataPoint] = []
    for _, row in df.iterrows():
        rig = str(row.get(0, "")).strip()
        if not rig:
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
    # Parse Meyerhof chart if present
    try:
        dr_col = df.columns[10]
        n_col = df.columns[11]
        for _, row in df.iterrows():
            depth_ratio = row.get(dr_col, None)
            n_value = row.get(n_col, None)
            if depth_ratio is not None and n_value is not None and pd.notna(depth_ratio) and pd.notna(n_value):
                try:
                    dr = float(depth_ratio)
                    nv = float(n_value)
                    meyerhof_profile.append(SoilDataPoint(depth=dr, value=nv))
                except Exception:
                    pass
    except Exception:
        pass
    meyerhof_profile.sort(key=lambda p: p.depth)
    return spuds, meyerhof_profile


def parse_soil_information(df: pd.DataFrame) -> Tuple[List[SoilLayer], float]:
    """Parse the `Soil Information` sheet into layers.

    Expects a raw DataFrame with no header.  Each soil layer begins
    with a row where column 0 (A) is non‑empty.  Columns 1 and 2
    contain top and bottom depths, column 4 (E) contains a text
    description used to classify the soil type.  Depth/value pairs for
    gamma' are in columns F (5) and G (6), φ in H (7) and I (8), and
    Su in J (9) and K (10).  The total analysis depth may appear
    in L3 (0‑based row index 2, column 11).  Returns the list of
    layers and the max depth.
    """
    layers: List[SoilLayer] = []
    max_depth = 0.0
    try:
        max_depth = float(df.iloc[2, 11])
    except Exception:
        max_depth = 0.0
    # Identify header rows
    header_indices = []
    for i in range(len(df)):
        val = df.iloc[i, 0]
        if pd.notna(val) and str(val).strip() != "":
            header_indices.append(i)
    for idx, row_idx in enumerate(header_indices):
        name = str(df.iloc[row_idx, 0])
        try:
            top_depth = float(df.iloc[row_idx, 1])
        except Exception:
            top_depth = 0.0
        try:
            bottom_depth = float(df.iloc[row_idx, 2])
        except Exception:
            bottom_depth = top_depth
        soil_desc = str(df.iloc[row_idx, 4]).lower() if pd.notna(df.iloc[row_idx, 4]) else ""
        if "clay" in soil_desc:
            soil_type = "clay"
        elif "silt" in soil_desc:
            soil_type = "silt"
        elif "sand" in soil_desc:
            soil_type = "sand"
        else:
            soil_type = "unknown"
        # End row for this layer
        end_row = header_indices[idx + 1] if idx + 1 < len(header_indices) else len(df)
        unit_profile: List[SoilDataPoint] = []
        phi_profile: List[SoilDataPoint] = []
        su_profile: List[SoilDataPoint] = []
        for j in range(row_idx + 1, end_row):
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


def SDPArrayHasData(arr: List[SoilDataPoint]) -> bool:
    """Return True if the SoilDataPoint list has data."""
    return bool(arr)


def SDPArrayHasData_2D(arr: List[List[SoilDataPoint]]) -> bool:
    return any(arr)