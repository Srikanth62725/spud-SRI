
from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np
import pandas as pd

USE_MIN_CU_POINT_AVG = True
APPLY_PHI_REDUCTION = False
APPLY_WINDWARD_CAPACITY_FACTOR = False
APPLY_SQUEEZING_TRIGGER = True

@dataclass
class SoilDataPoint:
    depth: float
    value: float

@dataclass
class SoilLayer:
    name: str
    top: float
    bottom: float
    soil_type: str
    gamma_profile: List[SoilDataPoint]
    su_profile: List[SoilDataPoint]
    phi_profile: List[SoilDataPoint]

@dataclass
class Spudcan:
    rig_name: str
    diameter: float
    area: float
    tip_elev: float
    max_preload_kN: float

def _interp_profile(z: float, prof: List[SoilDataPoint]) -> float:
    if not prof:
        return 0.0
    prof = sorted(prof, key=lambda p: p.depth)
    if z <= prof[0].depth:
        return prof[0].value
    for i in range(1, len(prof)):
        if z <= prof[i].depth:
            z1, v1 = prof[i-1].depth, prof[i-1].value
            z2, v2 = prof[i].depth, prof[i].value
            if z2 == z1:
                return v1
            return v1 + (z - z1) * (v2 - v1) / (z2 - z1)
    return prof[-1].value

def _get_layer_index(z: float, layers: List[SoilLayer]) -> int:
    for i, L in enumerate(layers):
        if z >= L.top and z < L.bottom:
            return i
    return len(layers) - 1

def _get_prop_at(z: float, layers: List[SoilLayer], prop: str) -> float:
    if not layers:
        return 0.0
    i = _get_layer_index(z, layers)
    L = layers[i]
    if prop == "su":
        return max(0.0, _interp_profile(z, L.su_profile))
    if prop == "gamma":
        return _interp_profile(z, L.gamma_profile)
    if prop == "phi":
        return _interp_profile(z, L.phi_profile)
    return 0.0

def _avg_between(z1: float, z2: float, layers: List[SoilLayer], prop: str) -> float:
    if z2 < z1:
        z1, z2 = z2, z1
    if abs(z2 - z1) < 1e-9:
        return _get_prop_at(z2, layers, prop)
    n = max(5, int((z2 - z1) / 0.05))
    zz = np.linspace(z1, z2, n)
    vals = [_get_prop_at(z, layers, prop) for z in zz]
    return float(np.mean(vals)) if len(vals) else 0.0

def _avg_over_B2(start_z: float, B: float, layers: List[SoilLayer], prop: str) -> float:
    return _avg_between(start_z, start_z + 0.5*B, layers, prop)

def _overburden(z: float, layers: List[SoilLayer]) -> float:
    if z <= 0:
        return 0.0
    dz = 0.05
    s = 0.0
    x = 0.0
    while x < z:
        g = _get_prop_at(x + 0.5*dz, layers, "gamma")
        s += g * dz
        x += dz
    return s

def _layer_type(z: float, layers: List[SoilLayer]) -> str:
    if not layers:
        return "unknown"
    return layers[_get_layer_index(z, layers)].soil_type

def _interp_meyerhof(db_over_B: float, meyerhof_df: Optional[pd.DataFrame]) -> float:
    if meyerhof_df is None or meyerhof_df.empty:
        return 0.0
    df = meyerhof_df.sort_values(by="D_over_B")
    x = df["D_over_B"].to_numpy()
    y = df["N"].to_numpy()
    if db_over_B <= x[0]:
        return float(y[0])
    if db_over_B >= x[-1]:
        return float(y[-1])
    i = np.searchsorted(x, db_over_B)
    x1,x2 = x[i-1], x[i]
    y1,y2 = y[i-1], y[i]
    return float(y1 + (db_over_B - x1)*(y2 - y1)/(x2 - x1))

def backflow_flag(d: float, B: float, layers: List[SoilLayer], meyerhof_df: Optional[pd.DataFrame]) -> bool:
    if d <= 0 or B <= 0:
        return False
    if _layer_type(d, layers) not in ("clay","silt"):
        return False
    N = _interp_meyerhof(d/B, meyerhof_df)
    cu_avg = _avg_over_B2(d, B, layers, "su")
    gamma_avg = _avg_over_B2(d, B, layers, "gamma")
    if gamma_avg <= 0:
        return False
    return d > (N * cu_avg) / gamma_avg

def clay_capacity_kN(B: float, A: float, d: float, layers: List[SoilLayer],
                     has_backflow: bool, use_min_cu_point_avg: bool = USE_MIN_CU_POINT_AVG) -> float:
    cu_point = _get_prop_at(d, layers, "su")
    cu_avg = _avg_over_B2(d, B, layers, "su")
    if use_min_cu_point_avg and cu_point > 0 and cu_avg > 0:
        cu_eff = min(cu_point, cu_avg)
    else:
        cu_eff = cu_avg if cu_avg > 0 else cu_point
    if cu_eff <= 0:
        return 0.0
    Nc = 5.14
    sc = 1.2
    dc = 1.0 + 0.4 * (d / B) if B > 0 else 1.0
    p0 = 0.0 if has_backflow else _overburden(d, layers)
    Fv = (cu_eff * Nc * sc * dc + p0) * A
    return float(max(0.0, Fv))

def sand_capacity_kN(B: float, A: float, d: float, layers: List[SoilLayer],
                     apply_phi_reduction: bool = APPLY_PHI_REDUCTION) -> float:
    gamma_prime = _avg_over_B2(d, B, layers, "gamma")
    phi = _get_prop_at(d, layers, "phi")
    if apply_phi_reduction:
        phi = max(0.0, phi - 5.0)
    if phi <= 0.0:
        return _overburden(d, layers) * A
    phi_rad = np.deg2rad(phi)
    Nq = np.exp(np.pi*np.tan(phi_rad)) * (np.tan(np.pi/4 + phi_rad/2)**2)
    Ngamma = 2.0 * (Nq + 1.0) * np.tan(phi_rad)
    sq = 1.0 + np.tan(phi_rad)
    sgamma = 0.6
    dq = 1.0 + 2*np.tan(phi_rad) * (1 - np.sin(phi_rad))**2 * (d/B) if B>0 else 1.0
    p0 = _overburden(d, layers)
    Fv = (0.5 * gamma_prime * B * Ngamma * sgamma + p0 * Nq * sq * dq) * A
    return float(max(0.0, Fv))

def squeezing_capacity_kN(B: float, A: float, d: float, layers: List[SoilLayer], has_backflow: bool) -> float | None:
    i = _get_layer_index(d, layers)
    if i+1 >= len(layers):
        return None
    topL = layers[i]
    if topL.soil_type != "clay":
        return None
    cu_t = _avg_over_B2(d, B, layers, "su")
    cu_b = _avg_over_B2(topL.bottom, B, layers, "su")
    if cu_t <= 0 or cu_b <= 0:
        return None
    if not (cu_b > 1.5 * cu_t):
        return None
    T = topL.bottom - d
    if T <= 0:
        return None
    if APPLY_SQUEEZING_TRIGGER:
        if B <= 0:
            return None
        if B < 3.45 * T * (1.0 + 1.1 * (d / B)):
            return None
    p0 = 0.0 if has_backflow else _overburden(d, layers)
    Fv = A * ((5.0 + 0.33 * (B / T) + 1.2 * (d / B)) * cu_t + p0)
    return float(max(0.0, Fv))

def punch_through_capacity_kN(B: float, A: float, d: float, layers: List[SoilLayer], has_backflow: bool) -> float | None:
    i = _get_layer_index(d, layers)
    if i+1 >= len(layers):
        return None
    topL = layers[i]
    botL = layers[i+1]
    H = topL.bottom - d
    if H <= 0:
        return None
    if topL.soil_type == "clay" and botL.soil_type == "clay":
        cu_t = _avg_over_B2(d, B, layers, "su")
        cu_b = _avg_over_B2(topL.bottom, B, layers, "su")
        if not (cu_t > cu_b and B > 0):
            return None
        Nc = 5.14
        sc = 1.2
        p0 = 0.0 if has_backflow else _overburden(d + H, layers)
        Fv = (3.0 * (H / B) * cu_t + Nc * sc * (1 + 0.2 * (d + H) / B) * cu_b + p0) * A
        return float(max(0.0, Fv))
    if topL.soil_type == "sand" and botL.soil_type == "clay":
        Fv_b = clay_capacity_kN(B, A, d + H, layers, has_backflow=False)
        gamma_sand = _avg_between(d, topL.bottom, layers, "gamma")
        p0_sand = _overburden(d, layers)
        cu_clay = _avg_over_B2(topL.bottom, B, layers, "su")
        KsTanPhi = 0.0 if (B<=0 or gamma_sand<=0) else (3.0 * cu_clay) / (B * gamma_sand)
        Fv = Fv_b - A * H * gamma_sand + 2.0 * (H / B) * (H * gamma_sand + 2.0 * p0_sand) * KsTanPhi * A
        return float(max(0.0, Fv))
    return None

def compute_envelopes(spud: Spudcan, layers: List[SoilLayer], meyerhof_df: Optional[pd.DataFrame],
                      dz: float, dmax: float, use_min_cu_point_avg: bool,
                      apply_phi_reduction: bool, windward_factor: bool) -> pd.DataFrame:
    depths = np.arange(0.0, dmax + 1e-9, dz)
    rows = []
    for d in depths:
        has_bk = backflow_flag(d, spud.diameter, layers, meyerhof_df)
        idle_clay = clay_capacity_kN(spud.diameter, spud.area, d, layers, has_bk, use_min_cu_point_avg)
        idle_sand = sand_capacity_kN(spud.diameter, spud.area, d, layers, apply_phi_reduction)
        sq = squeezing_capacity_kN(spud.diameter, spud.area, d, layers, has_bk)
        pt = punch_through_capacity_kN(spud.diameter, spud.area, d, layers, has_bk)

        real_clay = idle_clay
        clay_mode = "General"
        if sq is not None and sq < real_clay:
            real_clay = float(sq); clay_mode = "Squeezing"
        if pt is not None and (_layer_type(d, layers) != "sand") and pt < real_clay:
            real_clay = float(pt); clay_mode = "Punch-Through"

        if np.isnan(idle_sand):
            real_sand = np.nan; sand_mode = "NA"
        else:
            real_sand = idle_sand; sand_mode = "General"
            if pt is not None and _layer_type(d, layers) == "sand" and pt < real_sand:
                real_sand = float(pt); sand_mode = "Punch-Through"

        if windward_factor:
            if not np.isnan(real_clay):
                real_clay *= 0.8
            if not np.isnan(real_sand):
                real_sand *= 0.8

        rows.append(dict(depth=d,
                         idle_clay_MN=idle_clay/1000.0 if idle_clay>0 else np.nan,
                         idle_sand_MN=idle_sand/1000.0 if idle_sand>0 else np.nan,
                         real_clay_MN=real_clay/1000.0 if real_clay>0 else np.nan,
                         real_sand_MN=real_sand/1000.0 if real_sand>0 else np.nan,
                         backflow="Yes" if has_bk else "No",
                         clay_mode=clay_mode, sand_mode=sand_mode,
                         squeeze_MN=np.nan if sq is None else sq/1000.0,
                         pt_MN=np.nan if pt is None else pt/1000.0))
    return pd.DataFrame(rows)

def penetration_at_load_MN(df: pd.DataFrame, target_MN: float, col: str) -> float:
    sub = df[['depth', col]].dropna()
    if sub.empty:
        return float('nan')
    x = sub[col].to_numpy()
    z = sub['depth'].to_numpy()
    idx = np.where(x >= target_MN)[0]
    if len(idx) == 0:
        return float(z[-1])
    i = idx[0]
    if i == 0:
        return float(z[0])
    x1, x2 = x[i-1], x[i]
    z1, z2 = z[i-1], z[i]
    if x2 == x1:
        return float(z1)
    frac = (target_MN - x1) / (x2 - x1)
    return float(z1 + frac * (z2 - z1))
