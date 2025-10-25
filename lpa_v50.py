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
# lpa_v50.py
# Leg Penetration Analysis (V50 streamlit build)
# Author: Srikanth & ChatGPT
# Units: depth in m, Su in kPa, gamma' in kN/m3, forces in kN; plotting in MN

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
import numpy as np
import pandas as pd

# ---------------- Version 50 switches (user-togglable from UI) ----------------
USE_MIN_CU_POINT_AVG_DEFAULT = True      # use min(point, B/2-average)
APPLY_PHI_REDUCTION_DEFAULT   = False    # subtract 5 degrees from phi' if True
APPLY_WINDWARD_FACTOR_DEFAULT = False    # multiply REAL by 0.8 if True
APPLY_SQUEEZE_TRIGGER_DEFAULT = True     # enforce B >= 3.45 T (1+1.1 D/B)

# ---------------- Data models ----------------
@dataclass
class Spudcan:
    rig_name: str
    B: float                 # diameter (m)
    A: float                 # area (m2)
    tip_elev: float          # distance from tip to widest section (m)
    preload_MN: float        # per leg, MN

@dataclass
class SoilPoint:
    z: float
    v: float

@dataclass
class SoilLayer:
    name: str
    z_top: float
    z_bot: float
    soil_type: str           # "clay", "sand", "silt", "unknown"
    gamma: List[SoilPoint] = field(default_factory=list)  # kN/m3 vs depth
    su:    List[SoilPoint] = field(default_factory=list)  # kPa vs depth
    phi:   List[SoilPoint] = field(default_factory=list)  # deg vs depth

# ---------------- Helpers: profiles ----------------
def _interp(depth: float, prof: List[SoilPoint]) -> float:
    """Linear interpolation with edge hold; NaN for empty."""
    if not prof:
        return np.nan
    prof = sorted(prof, key=lambda p: p.z)
    if depth <= prof[0].z:
        return prof[0].v
    for i in range(1, len(prof)):
        if depth <= prof[i].z:
            z1, v1 = prof[i-1].z, prof[i-1].v
            z2, v2 = prof[i].z, prof[i].v
            if z2 == z1:
                return v1
            return v1 + (depth - z1) * (v2 - v1) / (z2 - z1)
    return prof[-1].v

def _avg_over(z1: float, z2: float, prof: List[SoilPoint], dz: float = 0.05) -> float:
    if not prof or z2 <= z1:
        return np.nan
    zs = np.arange(z1, z2 + 1e-9, dz)
    vals = np.array([_interp(z, prof) for z in zs], dtype=float)
    vals = vals[~np.isnan(vals)]
    return float(vals.mean()) if vals.size else np.nan

def _layer_index(z: float, layers: List[SoilLayer]) -> int:
    for i, L in enumerate(layers):
        if L.z_top <= z < L.z_bot:
            return i
    return len(layers) - 1

def _has_su_here(z: float, layers: List[SoilLayer]) -> bool:
    i = _layer_index(z, layers)
    val = _interp(z, layers[i].su)
    return np.isfinite(val) and val > 0

def _has_phi_here(z: float, layers: List[SoilLayer]) -> bool:
    i = _layer_index(z, layers)
    val = _interp(z, layers[i].phi)
    return np.isfinite(val) and val > 0

def _gamma_prime(z: float, layers: List[SoilLayer]) -> float:
    i = _layer_index(z, layers)
    return _interp(z, layers[i].gamma)

def _overburden(z: float, layers: List[SoilLayer], dz: float = 0.1) -> float:
    """Vertical effective stress at depth z from seabed; integrates gamma'."""
    if z <= 0:
        return 0.0
    zs = np.arange(0.0, z + 1e-9, dz)
    gammas = np.array([_gamma_prime(zi, layers) for zi in zs], dtype=float)
    gammas[np.isnan(gammas)] = 0.0
    # trapezoid
    return float(np.trapz(gammas, zs))

# ---------------- Meyerhof stability N(d/B) ----------------
# Default curve (smooth, conservative); users can override via UI
_DEFAULT_MEYERHOF = pd.DataFrame({
    "D_over_B": [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0],
    "N":        [0.0,  2.0,  3.0,  3.6,  4.0,  4.7,  5.1],
})

def _meyerhof_N(d_over_B: float, table: Optional[pd.DataFrame]) -> float:
    df = table if table is not None else _DEFAULT_MEYERHOF
    x = np.asarray(df["D_over_B"], dtype=float)
    y = np.asarray(df["N"], dtype=float)
    if d_over_B <= x[0]:
        return float(y[0])
    if d_over_B >= x[-1]:
        return float(y[-1])
    j = np.searchsorted(x, d_over_B)
    x1, x2 = x[j-1], x[j]
    y1, y2 = y[j-1], y[j]
    return float(y1 + (d_over_B - x1) * (y2 - y1) / (x2 - x1))

# ---------------- Sand factors ----------------
def _Nq(phi_rad: float) -> float:
    return np.exp(np.pi * np.tan(phi_rad)) * np.tan(np.pi/4 + phi_rad/2)**2

def _Ngamma(phi_rad: float) -> float:
    return 2.0 * (_Nq(phi_rad) + 1.0) * np.tan(phi_rad)

# ---------------- Capacities (kN) ----------------
def clay_capacity(spud: Spudcan, z: float, layers: List[SoilLayer],
                  use_min_cu: bool, backflow_zero: bool) -> Optional[float]:
    B, A = spud.B, spud.A
    if B <= 0 or A <= 0:
        return None
    i = _layer_index(z, layers)
    if not _has_su_here(z, layers):
        return None
    cu_point = _interp(z, layers[i].su)
    cu_avg   = _avg_over(z, z + B/2.0, layers[i].su)
    cu_eff   = np.nanmin([cu_point, cu_avg]) if use_min_cu else cu_avg
    if not np.isfinite(cu_eff) or cu_eff <= 0:
        return None
    Nc, sc = 5.14, 1.2
    if B > 0:
        d_over_B = z / B
        dc = 1.0 + 0.2 * d_over_B if d_over_B <= 1.0 else 1.0 + 0.2 * np.arctan(d_over_B)
    else:
        dc = 1.0
    p0 = 0.0 if backflow_zero else _overburden(z, layers)
    Fv = (cu_eff * Nc * sc * dc + p0) * A
    return float(Fv)

def sand_capacity(spud: Spudcan, z: float, layers: List[SoilLayer],
                  apply_phi_reduction: bool) -> Optional[float]:
    B, A = spud.B, spud.A
    if B <= 0 or A <= 0:
        return None
    i = _layer_index(z, layers)
    if not _has_phi_here(z, layers):
        return None
    phi = _interp(z, layers[i].phi)
    if not np.isfinite(phi) or phi <= 0:
        return None
    if apply_phi_reduction:
        phi = max(0.0, phi - 5.0)
    phi_rad = np.deg2rad(phi)
    Nq = _Nq(phi_rad)
    Ng = _Ngamma(phi_rad)
    gamma_p = _gamma_prime(z, layers)
    p0 = _overburden(z, layers)
    # shape/depth factors (simple)
    sq, sg, dq, dg = 1.0 + np.tan(phi_rad), 0.6, 1.0 + 2.0*np.tan(phi_rad)*(1-np.sin(phi_rad))**2*(z/max(B,1e-6)), 1.0
    Fv = (0.5 * gamma_p * B * Ng * sg * dg + p0 * Nq * sq * dq) * A
    return float(max(Fv, 0.0))

# ---------------- Special modes (kN) ----------------
def squeeze_capacity(spud: Spudcan, z: float, layers: List[SoilLayer],
                     enforce_trigger: bool, backflow_zero: bool) -> Optional[float]:
    """Strong over weak clay; requires next layer exists and is stronger."""
    idx = _layer_index(z, layers)
    if idx + 1 >= len(layers):
        return None
    top, bot = layers[idx], layers[idx+1]
    if top.soil_type not in ("clay", "silt"):
        return None
    B, A = spud.B, spud.A
    cu_t = _avg_over(z, z + B/2.0, top.su)
    cu_b = _avg_over(top.z_bot, top.z_bot + B/2.0, bot.su)
    if not np.isfinite(cu_t) or not np.isfinite(cu_b):
        return None
    if cu_b <= 1.5 * cu_t:   # needs pronounced contrast
        return None
    T = top.z_bot - z
    if T <= 0:
        return None
    if enforce_trigger:
        if B < 3.45 * T * (1.0 + 1.1 * (z / max(B,1e-6))):
            return None
    p0 = 0.0 if backflow_zero else _overburden(z, layers)
    # lower-bound expression consistent with SNAME commentary intent (simple conservative form)
    # Fv = A * [(5 + 0.33*(B/T) + 1.2*(z/B)) * cu_t + p0]
    Fv = A * ((5.0 + 0.33*(B/T) + 1.2*(z/max(B,1e-6))) * cu_t + p0)
    return float(Fv)

def punchthrough_capacity(spud: Spudcan, z: float, layers: List[SoilLayer],
                          backflow_zero: bool) -> Optional[float]:
    """Two common cases: CLAY/CLAY with stronger below; SAND over CLAY."""
    idx = _layer_index(z, layers)
    if idx + 1 >= len(layers):
        return None
    top, bot = layers[idx], layers[idx+1]
    H = top.z_bot - z
    if H <= 0:
        return None
    B, A = spud.B, spud.A

    # clay over clay: top weaker than bottom? (punch risk when top is stronger and fails)
    if top.soil_type in ("clay","silt") and bot.soil_type in ("clay","silt"):
        cu_t = _avg_over(z, z + B/2.0, top.su)
        cu_b = _avg_over(top.z_bot, top.z_bot + B/2.0, bot.su)
        if not np.isfinite(cu_t) or not np.isfinite(cu_b):
            return None
        if cu_t <= cu_b:
            return None
        Nc, sc = 5.14, 1.2
        p0b = 0.0 if backflow_zero else _overburden(z + H, layers)
        Fv = A * (3.0*(H/max(B,1e-6))*cu_t + Nc*sc*(1.0 + 0.2*((z+H)/max(B,1e-6))) * cu_b + p0b)
        return float(Fv)

    # sand over clay
    if top.soil_type == "sand" and bot.soil_type in ("clay","silt"):
        # effective capacity at lower interface, minus sand self-weight effect + shear transfer
        Fv_b = clay_capacity(spud, z + H, layers, use_min_cu=True, backflow_zero=False)
        if Fv_b is None:
            return None
        gamma_s = _gamma_prime(top.z_bot, layers)
        p0_s   = _overburden(z, layers)
        cu_c   = _avg_over(top.z_bot, top.z_bot + B/2.0, bot.su)
        if not np.isfinite(gamma_s) or not np.isfinite(cu_c):
            return None
        KsTanPhi = (3.0 * cu_c) / (max(B,1e-6) * max(gamma_s,1e-6))
        Fv = Fv_b - A * H * gamma_s + 2.0 * (H/max(B,1e-6)) * (H*gamma_s + 2.0*p0_s) * KsTanPhi * A
        return float(Fv)

    return None

# ---------------- Master sweep ----------------
def compute_envelopes(
    spud: Spudcan,
    layers: List[SoilLayer],
    max_depth: float,
    dz: float = 0.25,
    use_min_cu: bool = USE_MIN_CU_POINT_AVG_DEFAULT,
    phi_reduction: bool = APPLY_PHI_REDUCTION_DEFAULT,
    windward_factor: bool = APPLY_WINDWARD_FACTOR_DEFAULT,
    squeeze_trigger: bool = APPLY_SQUEEZE_TRIGGER_DEFAULT,
    meyerhof_table: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Returns a dataframe with per-depth results. Capacities stored in MN."""
    n = int(np.floor(max_depth / dz)) + 1
    depths = np.round(np.linspace(0.0, n*dz, n+1)[:n], 6)  # avoid fp drift
    out: Dict[str, list] = {
        "depth": [],
        "idle_clay_MN": [],
        "idle_sand_MN": [],
        "real_MN": [],
        "gov": [],
        "backflow": [],
        "squeeze_MN": [],
        "punch_MN": [],
        "real_clay_only_MN": [],
        "real_sand_only_MN": [],
    }

    for z in depths:
        # backflow check using Meyerhof
        N = _meyerhof_N(z / max(spud.B,1e-6), meyerhof_table)
        cu_avg = np.nan
        gamma_avg = np.nan
        if _has_su_here(z, layers):
            i = _layer_index(z, layers)
            cu_avg = _avg_over(z, z + spud.B/2.0, layers[i].su)
            gamma_avg = _avg_over(z, z + spud.B/2.0, layers[i].gamma)
        backflow = False
        if np.isfinite(cu_avg) and np.isfinite(gamma_avg) and gamma_avg > 0:
            backflow = z > (N * cu_avg) / gamma_avg

        # idle envelopes (compute only where property exists)
        Fc = clay_capacity(spud, z, layers, use_min_cu=use_min_cu, backflow_zero=backflow)
        Fs = sand_capacity(spud, z, layers, apply_phi_reduction=phi_reduction)

        # specials (only meaningful where they return a finite number)
        Fsq = squeeze_capacity(spud, z, layers, enforce_trigger=squeeze_trigger, backflow_zero=backflow)
        Fpt = punchthrough_capacity(spud, z, layers, backflow_zero=backflow)

        # Compose candidates for REAL at this depth:
        # 1) If clay is viable, it may be reduced by the least applicable special.
        # 2) If sand is viable, it may be reduced by punch-through when applicable.
        real_clay = None
        if Fc is not None:
            real_clay = Fc
            if Fsq is not None:
                real_clay = min(real_clay, Fsq)
            if Fpt is not None and layers[_layer_index(z, layers)].soil_type != "sand":
                real_clay = min(real_clay, Fpt)

        real_sand = None
        if Fs is not None:
            real_sand = Fs
            if Fpt is not None and layers[_layer_index(z, layers)].soil_type == "sand":
                real_sand = min(real_sand, Fpt)

        # Apply windward factor to REAL candidates only (if enabled)
        if windward_factor:
            if real_clay is not None:
                real_clay *= 0.8
            if real_sand is not None:
                real_sand *= 0.8

        # Select governing REAL
        gov = ""
        real = None
        if (real_clay is not None) and (real_sand is not None):
            if real_clay <= real_sand:
                real = real_clay; gov = "Clay-governed"
            else:
                real = real_sand; gov = "Sand-governed"
        elif real_clay is not None:
            real = real_clay; gov = "Clay-only"
        elif real_sand is not None:
            real = real_sand; gov = "Sand-only"

        # write MN
        out["depth"].append(z)
        out["idle_clay_MN"].append(np.nan if Fc is None else Fc/1000.0)
        out["idle_sand_MN"].append(np.nan if Fs is None else Fs/1000.0)
        out["real_MN"].append(np.nan if real is None else real/1000.0)
        out["gov"].append(gov if gov else "NA")
        out["backflow"].append("Yes" if backflow else "No")
        out["squeeze_MN"].append(np.nan if Fsq is None else Fsq/1000.0)
        out["punch_MN"].append(np.nan if Fpt is None else Fpt/1000.0)
        out["real_clay_only_MN"].append(np.nan if real_clay is None else real_clay/1000.0)
        out["real_sand_only_MN"].append(np.nan if real_sand is None else real_sand/1000.0)

    return pd.DataFrame(out)

# ---------------- Penetration utility ----------------
def _penetration_for_load_MN(df: pd.DataFrame, col: str, load_MN: float) -> Optional[float]:
    """Linear interpolate depth at which capacity == load; returns None if no crossing."""
    x = df[col].to_numpy(dtype=float)
    z = df["depth"].to_numpy(dtype=float)
    mask = np.isfinite(x)
    x = x[mask]; z = z[mask]
    if x.size < 2:
        return None
    # need first z where x>=load; assume monotonic nondecreasing is typical
    if load_MN > x.max():
        return None
    j = np.where(x >= load_MN)[0]
    if j.size == 0:
        return None
    j = j[0]
    if j == 0:
        return float(z[0])
    x1, x2 = x[j-1], x[j]
    z1, z2 = z[j-1], z[j]
    if x2 == x1:
        return float(z2)
    return float(z1 + (load_MN - x1) * (z2 - z1) / (x2 - x1))

def penetration_results(spud: Spudcan, df: pd.DataFrame) -> Dict[str, Optional[float]]:
    """Returns analysis-depth and tip-penetration results (m). Range is reported when both envelopes exist."""
    P = spud.preload_MN
    # independently compute penetrations from clay-only and sand-only REALs
    z_clay = _penetration_for_load_MN(df, "real_clay_only_MN", P)
    z_sand = _penetration_for_load_MN(df, "real_sand_only_MN", P)

    res: Dict[str, Optional[float]] = {
        "z_clay": None if z_clay is None else float(z_clay),
        "z_sand": None if z_sand is None else float(z_sand),
        "tip_clay": None if z_clay is None else float(z_clay + spud.tip_elev),
        "tip_sand": None if z_sand is None else float(z_sand + spud.tip_elev),
        "z_range_min": None,
        "z_range_max": None,
        "tip_range_min": None,
        "tip_range_max": None,
    }
    candidates = [v for v in [z_clay, z_sand] if v is not None]
    if candidates:
        res["z_range_min"] = min(candidates)
        res["z_range_max"] = max(candidates)
        res["tip_range_min"] = res["z_range_min"] + spud.tip_elev
        res["tip_range_max"] = res["z_range_max"] + spud.tip_elev
    return res
