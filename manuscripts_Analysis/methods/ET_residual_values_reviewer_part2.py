#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q3_ET_partitioning_adjustment_diagnostic.py

Quantifies how often the ET partitioning adjustment (when modeled evaporation
exceeds observed ET) occurred in the final analytical dataset, and assesses
the monthly impact of the exceedance.

This script replicates the exact production PM evaporation calculation
on the full site time series (preserving wet‑mask history and site‑wide
radiation‑coverage decisions), then restricts to the final site‑year‑month
domain. It captures evaporation immediately before the final clip(upper=ET_pos).

Outputs summary statistics to the console.
"""

import os
import math
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS (adjust if necessary)
# ============================================================================
MONTHLY_DATA = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
SITE_INFO_CSV = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
AMERI_INPUT_DIR = r"M:\Research\WUE_CUE\ameri_data\ameri_fill__lai_precip"

# ============================================================================
# COAST CLASSIFIER (copy from Study_area_august.py)
# ============================================================================
EARTH_RADIUS_KM = 6371.0088
COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9),
        (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2),
        (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0),
        (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7),
        (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1),
        (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6),
        (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3),
        (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9),
        (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)
    ],
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
    lat0 = math.radians(lat)
    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y])
    p = np.array([0.0, 0.0])
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)
    ab = b - a
    denom = np.dot(ab, ab)
    if denom == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, np.dot(p - a, ab) / denom))
    return float(np.linalg.norm(p - (a + t * ab)))

def _distance_to_polyline(lat, lon, polyline):
    return min(
        _seg_dist(lat, lon, *polyline[i], *polyline[i+1])
        for i in range(len(polyline)-1)
    )

def assign_coast_region(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return np.nan
    lat, lon = float(lat), float(lon)
    if lat > 50:
        return "AK Coast"
    dists = {
        coast: _distance_to_polyline(lat, lon, polyline)
        for coast, polyline in COAST_POLYLINES.items()
    }
    return min(dists, key=dists.get)

# ============================================================================
# PRODUCTION PARTITIONING FUNCTIONS (copied from 16-ET_partioning_may_2026.py)
# ============================================================================

PRINT = "[PM_Partitioning]"

def _to_num(x):
    return pd.to_numeric(x, errors="coerce")

def _median_dt_seconds(index, fallback=1800.0):
    try:
        di = pd.Index(pd.to_datetime(index)).diff().total_seconds()
        mode = pd.Series(di).dropna().mode()
        return float(mode.iloc[0]) if len(mode) else float(fallback)
    except Exception:
        return float(fallback)

def _solar_declination(doy):
    return 0.409 * np.sin(2.0 * np.pi * (doy - 81.0) / 365.0)

def _inv_rel_dist(doy):
    return 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)

def _sunset_hour_angle(lat_rad, delta):
    cosw = -np.tan(lat_rad) * np.tan(delta)
    return np.arccos(np.clip(cosw, -1.0, 1.0))

def _Ra_MJ_m2_day(lat_deg, doy):
    lat = np.radians(lat_deg)
    delta = _solar_declination(doy)
    dr = _inv_rel_dist(doy)
    ws = _sunset_hour_angle(lat, delta)
    Gsc = 0.0820
    return (24.0 * 60.0 / np.pi) * Gsc * dr * (ws*np.sin(lat)*np.sin(delta) + np.cos(lat)*np.cos(delta)*np.sin(ws))

def _doy_from_index(idx):
    try:
        return pd.DatetimeIndex(idx).dayofyear
    except Exception:
        return pd.Series(182, index=idx)

def rn_available_fao56(df, lat_deg, elev_m, albedo=0.23):
    """FAO-56 net radiation (W m-2) from Rg and Tair."""
    idx = df.index
    Rg = _to_num(df.get("Rg_f"))
    Tair = _to_num(df.get("Tair_f"))
    if Rg is None or Tair is None or Rg.isna().all() or Tair.isna().all():
        return pd.Series(np.nan, index=idx, name="Rn")

    Tk = Tair + 273.15
    doy = _doy_from_index(idx)
    Ra = _Ra_MJ_m2_day(lat_deg, doy)
    Ra_Wm2 = Ra * 1e6 / (24.0*3600.0)

    Rso = (0.75 + 2e-5*float(elev_m)) * Ra_Wm2
    Rso = pd.Series(Rso, index=idx).clip(lower=1.0)

    Rns = (1.0 - albedo) * Rg

    es_kPa = 0.6108 * np.exp(17.27*Tair / (Tair + 237.3))
    VPD_hPa = _to_num(df.get("VPD_f"))
    VPD_kPa = (VPD_hPa/10.0) if VPD_hPa is not None else pd.Series(np.nan, index=idx)
    ea_kPa = (es_kPa - VPD_kPa).clip(lower=0.05)

    sigma = 5.67e-8
    Rs_Rso = (Rg / Rso).clip(0.3, 1.0)
    cloud = (1.35*Rs_Rso - 0.35).clip(0.05, 1.0)
    emiss = (0.34 - 0.14*np.sqrt(ea_kPa)).clip(0.10, 0.34)
    Rnl = sigma * Tk**4 * emiss * cloud

    Rn = (Rns - Rnl).clip(lower=0.0)
    return Rn.rename("Rn")

def rn_available(df, site, meta):
    """
    Exact production rn_available().
    Returns net radiation (W m-2) using measured NETRAD if coverage >= 50%,
    otherwise FAO-56, with 0.45*Rg fallback.
    """
    idx = df.index
    total_rows = len(idx)
    lat = meta.get(site, {}).get("lat", None)
    elev = meta.get(site, {}).get("elevation_m", 0.0)

    # 1. Check measured Rn columns
    measured_col = None
    measured_Rn = None
    measured_coverage = 0.0
    for c in ("NETRAD_f", "NETRAD", "Rn", "Rn_f"):
        if c in df.columns:
            s = _to_num(df[c])
            if s.notna().any():
                measured_col = c
                measured_Rn = s.clip(lower=0.0)
                measured_coverage = measured_Rn.notna().sum() / total_rows
                break

    # 2. FAO-56 estimate
    fao_Rn = None
    if lat is not None and "Rg_f" in df.columns:
        try:
            fao_Rn = rn_available_fao56(df, float(lat), float(elev))
        except Exception:
            fao_Rn = None

    # 3. Decision – EXACT production logic (no filling of measured)
    COVERAGE_THRESHOLD = 0.50

    if measured_Rn is not None and measured_coverage >= COVERAGE_THRESHOLD:
        # Use measured Rn as‑is; do NOT fill missing values
        Rn = measured_Rn
    elif measured_Rn is not None and measured_coverage < COVERAGE_THRESHOLD:
        # Ignore measured, use FAO or fallback
        if fao_Rn is not None:
            Rn = fao_Rn
        else:
            Rg = _to_num(df.get("Rg_f"))
            if Rg is not None and Rg.notna().any():
                Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            else:
                Rn = pd.Series(np.nan, index=idx, name="Rn")
    else:
        # No measured column
        if fao_Rn is not None:
            Rn = fao_Rn
        else:
            Rg = _to_num(df.get("Rg_f"))
            if Rg is not None and Rg.notna().any():
                Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            else:
                Rn = pd.Series(np.nan, index=idx, name="Rn")

    return Rn.rename("Rn")

def guess_ra(idx, h_m, WS, ra_floor=15.0, kB1=2.3):
    try:
        h = float(h_m)
        if not np.isfinite(h) or h <= 0:
            h = 0.3
    except Exception:
        h = 0.3
    z = max(3.0*h, 2.0)
    d = 0.67*h
    z0m = 0.1*h
    z0h = 0.1*h * np.exp(-float(kB1))

    WS = _to_num(WS)
    if WS is None:
        return pd.Series(float(ra_floor), index=idx, name="ra_s_m")
    WS = pd.Series(WS, index=idx).clip(lower=0.5)

    num_m = (z - d) / max(z0m, 1e-6)
    num_h = (z - d) / max(z0h, 1e-6)
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = (np.log(num_m) * np.log(num_h)) / (0.4**2 * WS)
    return pd.Series(ra, index=idx, name="ra_s_m").replace([np.inf, -np.inf], np.nan).fillna(float(ra_floor)).clip(lower=float(ra_floor))

def beer_partition(LAI, idx, k_beer=0.5):
    L = pd.Series(0.0, index=idx) if LAI is None else _to_num(LAI).reindex(idx)
    L = L.fillna(0.0).clip(lower=0.0)
    f_c = (1.0 - np.exp(-float(k_beer) * L)).clip(0.0, 1.0)
    return f_c.rename("f_canopy"), (1.0 - f_c).rename("f_soil")

def soil_resistance(LAI, idx, rsoil_min=150.0, rsoil_max=600.0, k_lai=0.5):
    L = pd.Series(0.0, index=idx) if LAI is None else _to_num(LAI).reindex(idx)
    L = L.fillna(0.0).clip(lower=0.0)
    rs = float(rsoil_min) * np.exp(float(k_lai) * L)
    return pd.Series(rs, index=idx, name="rs_soil").clip(rsoil_min, rsoil_max)

def std_pressure_from_elev_kpa(elev_m):
    z = float(elev_m) if pd.notna(elev_m) else 0.0
    return 101.325 * (1.0 - 2.25577e-5 * z) ** 5.25588

def pressure_for_equation(site, df, meta):
    idx = df.index
    P = _to_num(df.get("PA_f"))
    if P is not None and (P > 0).any():
        elev = meta.get(site, {}).get("elevation_m", 0.0)
        Pfallback = std_pressure_from_elev_kpa(elev)
        return pd.Series(P, index=idx).where(P > 0, Pfallback).rename("PA_kPa")
    P = _to_num(df.get("PA"))
    if P is not None and (P > 0).any():
        elev = meta.get(site, {}).get("elevation_m", 0.0)
        Pfallback = std_pressure_from_elev_kpa(elev)
        return pd.Series(P, index=idx).where(P > 0, Pfallback).rename("PA_kPa")
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    return pd.Series(std_pressure_from_elev_kpa(elev), index=idx, name="PA_kPa")

def wet_mask_from(df, site, params, meta):
    idx = df.index
    precip = _to_num(df.get("precip_mm_per_30min"))
    if precip is None:
        precip = _to_num(df.get("precip_mm"))
    precip = pd.Series(0.0, index=idx) if precip is None else precip.reindex(idx).fillna(0.0)

    thr = float(params.get("precip_threshold_mm", 0.6))
    ext_h = float(params.get("wet_extension_hours", 2.0))

    wet = (precip >= thr)

    dt_s = _median_dt_seconds(idx)
    spph = max(1, int(round(3600.0 / dt_s)))
    win = max(1, int(round(ext_h * spph)))
    wet = wet.rolling(window=win, min_periods=1).max().astype(bool)

    Rn_min = params.get("require_Rn_min", None)
    if Rn_min is not None:
        Rn = rn_available(df, site, meta)
        low = (Rn < float(Rn_min)) & (Rn > 0)
        wet = wet & (~low)

    return pd.Series(wet, index=idx, name="wet_mask")

def parse_canopy_height(height_str):
    if pd.isna(height_str):
        return 0.3
    if isinstance(height_str, (int, float)):
        return float(height_str)
    height_str = str(height_str)
    if '–' in height_str or '-' in height_str:
        height_str = height_str.replace('–', '-')
        parts = height_str.split('-')
        try:
            low = float(parts[0].strip())
            high = float(parts[1].strip())
            return (low + high) / 2.0
        except:
            return 0.3
    try:
        return float(height_str)
    except:
        return 0.3

# Parameter sets (copied from production)
PARAMETER_SETS = {
    'WET': {
        'wet_extension_hours': 2.0,
        'precip_threshold_mm': 0.8,
        'k_beer': 0.48,
        'rsoil_min': 250.0,
        'rsoil_max': 600.0,
        'G_frac': 0.12,
        'ra_floor': 18.0,
        'kB1': 2.6,
        'rs_wet_canopy': 70.0,
        'require_Rn_min': 50.0,
        'wet_E_cap_frac': 0.60,
        'night_T_frac': 0.05,
    },
    'SHRUB_GRASS_CROP': {
        'wet_extension_hours': 1.0,
        'precip_threshold_mm': 0.6,
        'k_beer': 0.58,
        'rsoil_min': 250.0,
        'rsoil_max': 500.0,
        'G_frac': 0.13,
        'ra_floor': 15.0,
        'kB1': 2.3,
        'rs_wet_canopy': 60.0,
        'require_Rn_min': 30.0,
        'wet_E_cap_frac': 0.50,
        'night_T_frac': 0.05,
    },
    'BSV': {
        'wet_extension_hours': 0.5,
        'precip_threshold_mm': 0.8,
        'k_beer': 0.35,
        'rsoil_min': 350.0,
        'rsoil_max': 800.0,
        'G_frac': 0.15,
        'ra_floor': 25.0,
        'kB1': 3.0,
        'rs_wet_canopy': 100.0,
        'require_Rn_min': 80.0,
        'wet_E_cap_frac': 0.45,
        'night_T_frac': 0.03,
    },
    'FOREST': {
        'wet_extension_hours': 3.0,
        'precip_threshold_mm': 0.2,
        'k_beer': 0.55,
        'rsoil_min': 140.0,
        'rsoil_max': 400.0,
        'G_frac': 0.06,
        'ra_floor': 12.0,
        'kB1': 2.0,
        'rs_wet_canopy': 50.0,
        'require_Rn_min': None,
        'wet_E_cap_frac': 0.62,
        'night_T_frac': 0.08,
    },
}

def biome_defaults(biome):
    b = str(biome).upper()
    if b == "WET":
        return PARAMETER_SETS['WET'].copy()
    elif b in ["GRA", "CSH", "OSH", "CRO"]:
        return PARAMETER_SETS['SHRUB_GRASS_CROP'].copy()
    elif b == "BSV":
        return PARAMETER_SETS['BSV'].copy()
    elif b in ["ENF", "DBF", "MF"]:
        return PARAMETER_SETS['FOREST'].copy()
    else:
        return PARAMETER_SETS['SHRUB_GRASS_CROP'].copy()

def load_site_meta(site_info_csv):
    encodings = ['utf-8', 'latin1', 'cp1252', 'iso-8859-1', 'windows-1252']
    info = None
    for encoding in encodings:
        try:
            info = pd.read_csv(site_info_csv, encoding=encoding)
            break
        except:
            continue
    if info is None:
        raise ValueError(f"Could not read {site_info_csv}")
    meta = {}
    for _, r in info.iterrows():
        site = str(r["site_name"]).strip()
        h_str = r.get("canopy_height_m", 0.3)
        h_canopy = parse_canopy_height(h_str)
        elev = r.get("elevation_m", 0.0)
        if pd.isna(elev):
            elev = 0.0
        else:
            elev = float(elev)
        lat = r.get("lat", None)
        if lat is not None and not pd.isna(lat):
            lat = float(lat)
        else:
            lat = 40.0
        meta[site] = dict(
            lat=lat,
            long=float(r.get("long", 0.0)),
            elevation_m=elev,
            biome=str(r.get("IGBP", "WET")),
            canopy_height_m=h_canopy
        )
    return meta

def _params_for_site(site, meta):
    site_info = meta.get(site, {})
    biome = site_info.get("biome", "WET")
    h_canopy = float(site_info.get("canopy_height_m", 0.3) or 0.3)
    params = biome_defaults(biome)
    return params, h_canopy, biome

# ============================================================================
# DIAGNOSTIC PARTITIONING FUNCTION (returns pre-clip evaporation)
# ============================================================================
def pm_partition_diagnostic(df, site, h_canopy_m, params, meta):
    """
    Exact production pm_partition(), but returns evaporation immediately
    after the wet-period cap and before the final clip(upper=ET_pos).
    Also returns the final clipped evaporation for comparison (optional).
    """
    idx = df.index
    ET = _to_num(df.get("ET"))
    if ET is None or ET.isna().all():
        return pd.Series(np.nan, index=idx, name="Evap_raw"), pd.Series(np.nan, index=idx, name="Evap_final")

    Tair = _to_num(df.get("Tair_f"))
    VPD_h = _to_num(df.get("VPD_f"))
    WS = _to_num(df.get("WS"))
    LAI = _to_num(df.get("lai"))
    if LAI is None or LAI.isna().all():
        LAI = _to_num(df.get("LAI_daily"))

    dt_s = _median_dt_seconds(idx, 1800.0)
    conv = dt_s / 2.45e6

    PkPa = pressure_for_equation(site, df, meta)
    gamma = 0.000665 * PkPa
    VPD_kP = (VPD_h / 10.0) if VPD_h is not None else pd.Series(np.nan, index=idx)
    VPD_kP = VPD_kP.clip(lower=0.0)

    if Tair is not None and Tair.notna().any():
        es = 0.6108 * np.exp(17.27 * Tair / (Tair + 237.3))
        delta = 4098.0 * es / ((Tair + 237.3) ** 2)
    else:
        delta = pd.Series(0.15, index=idx)

    R = 8.314462618
    M_air = 0.029
    Cp = 1010.0
    rho = ((PkPa * 1000.0) / (R * (Tair + 273.15))) * M_air

    Rn = rn_available(df, site, meta)
    f_c, f_s = beer_partition(LAI, idx, params["k_beer"])
    Rnc = f_c * Rn
    Rns = f_s * Rn
    G = (float(params["G_frac"]) * Rns).clip(upper=0.30 * Rns)
    Rns_eff = (Rns - G).clip(lower=0.0)

    ra = guess_ra(idx, h_canopy_m, WS, params["ra_floor"], params["kB1"])
    rs_soil = soil_resistance(LAI, idx, params["rsoil_min"], params["rsoil_max"], 0.5)
    rs_wet = pd.Series(float(params["rs_wet_canopy"]), index=idx)
    wet = wet_mask_from(df, site, params, meta)

    Evap_raw = pd.Series(np.nan, index=idx, name="Evap_raw")
    Evap_final = pd.Series(np.nan, index=idx, name="Evap_final")

    # Wet canopy evaporation
    m_wet = wet & Rnc.notna() & (Rnc > 0)
    if m_wet.any():
        num_energy = delta[m_wet] * Rnc[m_wet] * conv
        num_aero = (rho[m_wet] * Cp * VPD_kP[m_wet] / ra[m_wet]) * conv
        denom = (delta[m_wet] + gamma[m_wet] * (1.0 + rs_wet[m_wet] / ra[m_wet])).replace(0, np.nan)
        Evap_raw.loc[m_wet] = ((num_energy + num_aero) / denom).clip(lower=0.0)
        # Apply wet-period cap (this is the value we want to capture)
        cap = float(params.get("wet_E_cap_frac", 0.65))
        for i in idx[m_wet]:
            if pd.notna(ET.loc[i]) and ET.loc[i] > 0:
                Evap_raw.loc[i] = min(Evap_raw.loc[i], cap * ET.loc[i])

    # Dry soil evaporation
    m_dry = (~wet) & Rns_eff.notna() & (Rns_eff > 0)
    if m_dry.any():
        num_energy = delta[m_dry] * Rns_eff[m_dry] * conv
        num_aero = (rho[m_dry] * Cp * VPD_kP[m_dry] / ra[m_dry]) * conv
        denom = (delta[m_dry] + gamma[m_dry] * (1.0 + rs_soil[m_dry] / ra[m_dry])).replace(0, np.nan)
        Evap_raw.loc[m_dry] = ((num_energy + num_aero) / denom).clip(lower=0.0)

    # Final clipped version (same as production)
    ET_pos = ET.clip(lower=0)
    pos = ET_pos.notna() & (ET_pos > 0)
    Evap_final = Evap_raw.copy()
    Evap_final.loc[~pos] = np.nan
    Evap_final = Evap_final.clip(lower=0.0)
    Evap_final = Evap_final.where(ET_pos.notna(), np.nan).clip(upper=ET_pos)

    # Nighttime fallback (same as production)
    night_T_frac = float(params.get("night_T_frac", 0.05))
    missing_evap = (ET_pos > 0) & Evap_final.isna()
    if missing_evap.any():
        Trans_fallback = night_T_frac * ET_pos
        Evap_final.loc[missing_evap] = (ET_pos.loc[missing_evap] - Trans_fallback.loc[missing_evap]).clip(lower=0.0)

    return Evap_raw.rename("Evap_raw"), Evap_final.rename("Evap_final")

# ============================================================================
# STEP 1: RECONSTRUCT FINAL ANALYTICAL DOMAIN
# ============================================================================
print("Loading monthly data and filtering to final analytical sample...")
monthly = pd.read_csv(MONTHLY_DATA)
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

required = ['site_name', 'Year', 'month', 'water_class', 'lat', 'long',
            'Trans_ratio', 'WUE_tra']
for c in required:
    if c not in monthly.columns:
        raise ValueError(f"Column '{c}' not found in monthly data")

df = monthly.dropna(subset=required).copy()
df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])].copy()
df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) &
        np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) &
        (df['Trans_ratio'] <= 1) & np.isfinite(df['WUE_tra'])].copy()

df['coast_region'] = df.apply(lambda r: assign_coast_region(r['lat'], r['long']), axis=1)
df = df.dropna(subset=['coast_region']).copy()

final_domain = df[['site_name', 'Year', 'month']].drop_duplicates().copy()
final_domain['coast_region'] = final_domain['site_name'].map(
    df.groupby('site_name')['coast_region'].first()
)

site_coast = final_domain[['site_name', 'coast_region']].drop_duplicates().set_index('site_name')['coast_region'].to_dict()

print(f"Final analytical domain: {len(final_domain)} site-months, {final_domain['site_name'].nunique()} sites")

# ============================================================================
# STEP 2: RUN DIAGNOSTIC PARTITIONING ON FULL SITE, THEN RESTRICT
# ============================================================================
print("\nLoading site metadata...")
meta = load_site_meta(SITE_INFO_CSV)

# IMPORTANT: reset all containers to avoid cross-run duplication
all_events = []
monthly_impact = []
site_summaries = []   # <-- ADDED: reset from previous Spyder run

for site in final_domain['site_name'].unique():
    input_file = os.path.join(AMERI_INPUT_DIR, f"{site}.csv")
    if not os.path.exists(input_file):
        print(f"  WARNING: Input file for {site} not found at {input_file}. Skipping.")
        continue

    # ---- Robust CSV reading ----
    try:
        df_site = pd.read_csv(input_file, engine='python', on_bad_lines='skip')
    except Exception as e:
        print(f"  ERROR reading {site}: {e}")
        continue

    # Normalise datetime column
    dt_candidates = [col for col in df_site.columns if 'datetime' in col.lower()]
    if dt_candidates:
        dt_col = dt_candidates[0]
        df_site['DateTime'] = pd.to_datetime(df_site[dt_col], errors='coerce')
    else:
        if all(c in df_site.columns for c in ['Year', 'month', 'day', 'Hour']):
            df_site['DateTime'] = pd.to_datetime(
                df_site['Year'].astype(str) + '-' +
                df_site['month'].astype(str) + '-' +
                df_site['day'].astype(str) + ' ' +
                df_site['Hour'].astype(str) + ':00',
                errors='coerce'
            )
        else:
            print(f"  ERROR: {site} missing datetime column and cannot construct from components.")
            continue

    df_site = df_site.dropna(subset=['DateTime'])
    if len(df_site) == 0:
        print(f"  {site}: No valid datetime records.")
        continue

    df_site = df_site.set_index('DateTime')

    # 1. Run partitioning on the FULL site
    params, h_canopy, biome = _params_for_site(site, meta)
    Evap_raw_full, Evap_final_full = pm_partition_diagnostic(df_site, site, h_canopy, params, meta)

    # 2. Build diagnostic dataframe on full site
    diag = pd.DataFrame({
        'ET': _to_num(df_site.get('ET')),
        'Evap_raw': Evap_raw_full,
        'Evap_final': Evap_final_full
    }, index=df_site.index)
    diag['Year'] = diag.index.year
    diag['month'] = diag.index.month

    # 3. Restrict to final analytical domain
    domain_site = (
        final_domain[final_domain['site_name'] == site][['Year', 'month']]
        .drop_duplicates()
        .copy()
    )
    domain_site['Year'] = domain_site['Year'].astype(int)
    domain_site['month'] = domain_site['month'].astype(int)

    diag = (
        diag.reset_index()
        .merge(domain_site, on=['Year', 'month'], how='inner')
        .set_index('DateTime')
    )

    if len(diag) == 0:
        print(f"  {site}: No data in final domain.")
        continue

    ET = diag['ET']
    Evap_raw = diag['Evap_raw']

    # 4. Exceedance detection (reviewer-relevant denominator)
    valid = ET.notna() & (ET > 0) & Evap_raw.notna()
    mask_exceed = valid & (Evap_raw > ET)

    n_exceed = int(mask_exceed.sum())
    total_valid = int(valid.sum())

    if n_exceed > 0:
        excess_vals = Evap_raw[mask_exceed] - ET[mask_exceed]
        mean_excess = excess_vals.mean()
        max_excess = excess_vals.max()
        sum_excess = excess_vals.sum()
        for val in excess_vals:
            all_events.append({
                'site': site,
                'coast': site_coast.get(site, 'Unknown'),
                'excess_mm': val
            })
    else:
        mean_excess = np.nan
        max_excess = np.nan
        sum_excess = 0.0

    site_summaries.append({
        'site': site,
        'coast': site_coast.get(site, 'Unknown'),
        'total_valid': total_valid,
        'n_exceed': n_exceed,
        'frac_exceed': n_exceed / total_valid if total_valid > 0 else 0,
        'mean_excess': mean_excess,
        'max_excess': max_excess,
        'sum_excess': sum_excess,
    })

    # ---- NEW: Aggregate impact by site-month ----
    monthly_group = diag.groupby(['Year', 'month'])

    for (yr, mo), group in monthly_group:

        # Same valid definition used in the half-hourly diagnostic
        valid_month = (
            group['ET'].notna()
            & (group['ET'] > 0)
            & group['Evap_raw'].notna()
        )

        exceed_month = (
            valid_month
            & (group['Evap_raw'] > group['ET'])
        )

        # Sum ET only over valid positive-ET observations
        monthly_ET = group.loc[valid_month, 'ET'].sum()

        # Sum only genuine pre-clip exceedances
        monthly_excess = (
            group.loc[exceed_month, 'Evap_raw']
            - group.loc[exceed_month, 'ET']
        ).sum()

        n_valid_month = int(valid_month.sum())
        n_exceed_month = int(exceed_month.sum())

        if monthly_ET > 0:
            pct_exceed = 100.0 * monthly_excess / monthly_ET
        else:
            pct_exceed = np.nan

        monthly_impact.append({
            'site': site,
            'coast': site_coast.get(site, 'Unknown'),
            'Year': int(yr),
            'month': int(mo),
            'n_valid': n_valid_month,
            'n_exceed': n_exceed_month,
            'frac_exceed': (
                n_exceed_month / n_valid_month
                if n_valid_month > 0
                else np.nan
            ),
            'monthly_ET': monthly_ET,
            'monthly_excess': monthly_excess,
            'pct_exceed': pct_exceed
        })

# ============================================================================
# STEP 3: AGGREGATE AND PRINT RESULTS
# ============================================================================
print("\n" + "="*80)
print("ET PARTITIONING ADJUSTMENT DIAGNOSTIC (BEFORE FINAL CLIP)")
print("="*80)

if not site_summaries:
    print("No data processed.")
    exit()

df_site_sum = pd.DataFrame(site_summaries)

# Overall statistics (correctly calculated from all events)
total_valid_all = df_site_sum['total_valid'].sum()
total_exceed_all = df_site_sum['n_exceed'].sum()
frac_exceed_all = total_exceed_all / total_valid_all if total_valid_all > 0 else 0

if len(all_events) > 0:
    df_events = pd.DataFrame(all_events)
    overall_mean_excess = df_events['excess_mm'].mean()
    overall_max_excess = df_events['excess_mm'].max()
    overall_sum_excess = df_events['excess_mm'].sum()
    overall_median_excess = df_events['excess_mm'].median()
    overall_q25_excess = df_events['excess_mm'].quantile(0.25)
    overall_q75_excess = df_events['excess_mm'].quantile(0.75)
    overall_p95_excess = df_events['excess_mm'].quantile(0.95)
else:
    overall_mean_excess = np.nan
    overall_max_excess = np.nan
    overall_sum_excess = 0.0
    overall_median_excess = np.nan
    overall_q25_excess = np.nan
    overall_q75_excess = np.nan
    overall_p95_excess = np.nan

print("\nOVERALL (all sites, final domain):")
print(f"  Total valid half-hours (ET > 0 and Evap_raw available): {total_valid_all:,.0f}")
print(f"  Number of exceedances (Evap_raw > ET): {total_exceed_all:,.0f}")
print(f"  Percentage of valid records: {frac_exceed_all*100:.2f}%")
print(f"  Sum of excess (Evap_raw - ET) over all exceedances: {overall_sum_excess:.2f} mm")
print(f"  Mean excess per exceedance: {overall_mean_excess:.4f} mm")
print(f"  Median excess per exceedance: {overall_median_excess:.4f} mm")
print(f"  IQR: {overall_q25_excess:.4f} to {overall_q75_excess:.4f} mm")
print(f"  95th percentile: {overall_p95_excess:.4f} mm")
print(f"  Maximum excess: {overall_max_excess:.4f} mm")

# By coast
coast_summary = df_site_sum.groupby('coast').agg({
    'total_valid': 'sum',
    'n_exceed': 'sum',
    'sum_excess': 'sum',
}).reset_index()
coast_summary['frac_exceed'] = coast_summary['n_exceed'] / coast_summary['total_valid']
coast_summary['mean_excess'] = coast_summary['sum_excess'] / coast_summary['n_exceed']

print("\nBY COAST:")
print(coast_summary.to_string(index=False, float_format='%.3f'))

# Sites with highest fraction
print("\nTOP 5 SITES WITH HIGHEST EXCEEDANCE FRACTION:")
top5 = df_site_sum.nlargest(5, 'frac_exceed')[['site', 'coast', 'frac_exceed', 'n_exceed', 'total_valid']]
print(top5.to_string(index=False, float_format='%.3f'))

# Sites with zero exceedances
zero_sites = df_site_sum[df_site_sum['n_exceed'] == 0]['site'].tolist()
print(f"\nSites with zero exceedances: {len(zero_sites)} sites")
if len(zero_sites) > 0:
    print(f"  {', '.join(sorted(zero_sites))}")

# ============================================================================
# NEW: MONTHLY IMPACT SUMMARY (using the corrected aggregation)
# ============================================================================
print("\n" + "="*80)
print("MONTHLY IMPACT OF EXCEEDANCE (summed excess as % of monthly ET)")
print("="*80)

if monthly_impact:
    df_monthly = pd.DataFrame(monthly_impact)
    # Remove NaN percentages (months with zero ET)
    df_monthly_valid = df_monthly.dropna(subset=['pct_exceed'])

    print(f"Number of site-months with positive ET: {len(df_monthly_valid)}")
    if len(df_monthly_valid) > 0:
        med = df_monthly_valid['pct_exceed'].median()
        q25 = df_monthly_valid['pct_exceed'].quantile(0.25)
        q75 = df_monthly_valid['pct_exceed'].quantile(0.75)
        p95 = df_monthly_valid['pct_exceed'].quantile(0.95)
        max_pct = df_monthly_valid['pct_exceed'].max()
        # Find the site-month with the maximum percentage
        max_row = df_monthly_valid.loc[df_monthly_valid['pct_exceed'].idxmax()]
        print(f"  Median percentage of monthly ET exceeded: {med:.3f}%")
        print(f"  IQR: {q25:.3f}% to {q75:.3f}%")
        print(f"  95th percentile: {p95:.3f}%")
        print(f"  Maximum percentage: {max_pct:.3f}%")
        print(f"  Site-month with max: {max_row['site']} {int(max_row['Year'])}-{int(max_row['month']):02d} (coast: {max_row['coast']})")
    else:
        print("  No valid monthly percentages.")
else:
    print("  No monthly impact data collected.")

print("\nDiagnostic complete.")