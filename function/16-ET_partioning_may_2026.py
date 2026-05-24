# -*- coding: utf-8 -*-
"""
PM-based evaporation plus residual transpiration for ET partitioning.

CORE EQUATIONS:
    Evaporation (E_est) computed using Penman-Monteith formulation:
        E = [Δ(Rn - G) + ρ·Cp·(VPD/ra)] / [Δ + γ·(1 + rs/ra)]
    
    Transpiration is NOT modeled directly:
        T = ET_obs - E_est
    
    This is a residual partitioning approach.

Input files require: ET, Tair_f, VPD_f, precip_mm_per_30min (or precip_mm), lai

Output: Evap_pen and Trans_pen columns added to each site file.
"""

import os
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

PRINT = "[PM_Partitioning]"

# ============================================================================
# TUNABLE PARAMETERS - ADJUST THESE VALUES FOR YOUR STUDY SYSTEM
# ============================================================================

PARAMETER_SETS = {
    # Wetland group (WET) - coastal marshes, mangroves
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
        # Empirical constraint: during wet periods, evaporation is limited to a fraction of ET.
        # This acts as a regularization term in the residual partitioning framework,
        # preventing the model from attributing all ET to wet-surface evaporation,
        # which would suppress transpiration estimates due to the residual formulation.
        # This constraint is particularly important in coastal wetlands where
        # shallow water and wet surfaces can otherwise dominate the PM estimate.
        'wet_E_cap_frac': 0.60,
        # Nocturnal transpiration fraction (conservative lower bound from literature)
        # Wu et al. 2023: 5.5-24% in mangroves; using 5% as conservative fallback
        'night_T_frac': 0.05,
    },
    
    # Grass / Shrub / Crop group (GRA, CSH, OSH, CRO)
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
    
    # Barren/Sparse (BSV)
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
        'night_T_frac': 0.03,  # Lower for barren/sparse
    },
    
    # Forest group (ENF, DBF, MF)
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
        'night_T_frac': 0.08,  # Forests may have higher nocturnal T
    },
}


def biome_defaults(biome):
    """Get parameters based on biome/IGBP class."""
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
        print(PRINT, f"  Warning: Unknown biome '{b}', using shrub/grass defaults")
        return PARAMETER_SETS['SHRUB_GRASS_CROP'].copy()


def export_parameters_to_csv(output_dir):
    """Save parameter sets to CSV for reproducibility."""
    import pandas as pd
    rows = []
    for biome, params in PARAMETER_SETS.items():
        row = {'biome_group': biome}
        row.update(params)
        rows.append(row)
    df = pd.DataFrame(rows)
    out_path = Path(output_dir) / 'parameter_sets_used.csv'
    df.to_csv(out_path, index=False)
    print(PRINT, f"Parameter table saved to {out_path}")
    return out_path


# ============================================================================
# DIAGNOSTIC FUNCTION
# ============================================================================

# Set to True to save diagnostic variables (wet_mask, ra, rs_soil, f_canopy)
SAVE_DIAGNOSTICS = False


def create_diagnostic_report(annual_df, output_dir):
    """Create diagnostic plots and summary."""
    try:
        diag_dir = Path(output_dir) / 'diagnostics'
        diag_dir.mkdir(parents=True, exist_ok=True)
        
        # Plot 1: T/ET Distribution by Biome
        fig, ax = plt.subplots(figsize=(12, 6))
        biomes = annual_df.groupby('IGBP')['Trans_Fraction'].agg(['mean', 'std', 'count']).sort_values('mean')
        biomes = biomes[biomes['count'] >= 5]
        
        ax.barh(range(len(biomes)), biomes['mean'], xerr=biomes['std'], 
                capsize=3, alpha=0.7, color='steelblue')
        ax.set_yticks(range(len(biomes)))
        ax.set_yticklabels(biomes.index)
        ax.set_xlabel('Transpiration Fraction (T/ET)')
        ax.set_title('T/ET by Biome (with standard deviation)')
        ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(diag_dir / 't_et_by_biome.png', dpi=150)
        plt.close()
        
        # Plot 2: Distribution of T/ET values
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(annual_df['Trans_Fraction'], bins=30, alpha=0.7, color='steelblue', edgecolor='black')
        ax.axvline(x=annual_df['Trans_Fraction'].mean(), color='red', 
                  linestyle='--', label=f"Mean: {annual_df['Trans_Fraction'].mean():.3f}")
        ax.axvline(x=annual_df['Trans_Fraction'].median(), color='green', 
                  linestyle='--', label=f"Median: {annual_df['Trans_Fraction'].median():.3f}")
        ax.set_xlabel('Transpiration Fraction (T/ET)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Annual T/ET Across All Sites')
        ax.legend()
        plt.tight_layout()
        plt.savefig(diag_dir / 't_et_distribution.png', dpi=150)
        plt.close()
        
        print(PRINT, f"Diagnostic plots saved to {diag_dir}")
        return True
    except Exception as e:
        print(PRINT, f"Warning: Could not create diagnostic plots: {e}")
        return False


# ============================================================================
# CORE PHYSICAL FUNCTIONS
# ============================================================================

def _to_num(x):
    """Convert to numeric, coercing errors to NaN."""
    return pd.to_numeric(x, errors="coerce")


def _median_dt_seconds(index, fallback=1800.0):
    """Median timestep in seconds; robust to gaps."""
    try:
        di = pd.Index(pd.to_datetime(index)).diff().total_seconds()
        mode = pd.Series(di).dropna().mode()
        return float(mode.iloc[0]) if len(mode) else float(fallback)
    except Exception:
        return float(fallback)


def _solar_declination(doy):
    """Solar declination angle (radians)."""
    return 0.409 * np.sin(2.0 * np.pi * (doy - 81.0) / 365.0)


def _inv_rel_dist(doy):
    """Inverse relative Earth-Sun distance."""
    return 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)


def _sunset_hour_angle(lat_rad, delta):
    """Sunset hour angle (radians)."""
    cosw = -np.tan(lat_rad) * np.tan(delta)
    return np.arccos(np.clip(cosw, -1.0, 1.0))


def _Ra_MJ_m2_day(lat_deg, doy):
    """Extraterrestrial radiation (MJ m-2 day-1)."""
    lat = np.radians(lat_deg)
    delta = _solar_declination(doy)
    dr = _inv_rel_dist(doy)
    ws = _sunset_hour_angle(lat, delta)
    Gsc = 0.0820
    return (24.0 * 60.0 / np.pi) * Gsc * dr * (ws*np.sin(lat)*np.sin(delta) + np.cos(lat)*np.cos(delta)*np.sin(ws))


def _doy_from_index(idx):
    """Extract day of year from datetime index."""
    try:
        return pd.DatetimeIndex(idx).dayofyear
    except Exception:
        return pd.Series(182, index=idx)


def rn_available_fao56(df, lat_deg, elev_m, albedo=0.23):
    """Net radiation Rn (W m-2) via FAO-56."""
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
    Get net radiation with coverage-based decision. Does NOT fill/modify measured Rn.
    
    Strategy:
        1. Check measured Rn columns in order: NETRAD_f, NETRAD, Rn, Rn_f
        2. Calculate measured Rn coverage = valid measured values / total rows
        3. If measured Rn coverage >= 0.50:
           - use measured Rn as-is (accept missing values, do NOT fill)
        4. If measured Rn coverage < 0.50:
           - ignore measured Rn for this site (too sparse)
           - use FAO-56 Rn from Rg_f for the whole site
        5. If no measured Rn exists:
           - use FAO-56 Rn from Rg_f
        6. If FAO-56 fails but Rg_f exists:
           - use Rn = 0.45 * Rg_f as final fallback
        7. Return Series named "Rn" (never modifies original data)
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input data frame with meteorological columns
    site : str
        Site name (for metadata lookup)
    meta : dict
        Site metadata containing latitude and elevation
    
    Returns
    -------
    Rn : pandas.Series
        Net radiation (W m-2)
    """
    idx = df.index
    total_rows = len(idx)
    lat = meta.get(site, {}).get("lat", None)
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    
    # ================================================================
    # Step 1: Check for measured Rn columns
    # ================================================================
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
                print(PRINT, f"    Rn: found measured column '{c}' for {site} (coverage = {measured_coverage*100:.1f}%)")
                break
    
    # ================================================================
    # Step 2: Get FAO-56 Rn from Rg_f if possible
    # ================================================================
    fao_Rn = None
    if lat is not None and "Rg_f" in df.columns:
        try:
            fao_Rn = rn_available_fao56(df, float(lat), float(elev))
            fao_coverage = fao_Rn.notna().sum() / total_rows
            print(PRINT, f"    Rn: FAO-56 estimate available for {site} (coverage = {fao_coverage*100:.1f}%)")
        except Exception as e:
            print(PRINT, f"    Rn: FAO-56 failed for {site}: {e}")
            fao_Rn = None
    
    # ================================================================
    # Step 3: Decision based on measured coverage (NO FILLING)
    # ================================================================
    COVERAGE_THRESHOLD = 0.50
    
    # Case 1: measured column exists with good coverage (>= 50%)
    if measured_Rn is not None and measured_coverage >= COVERAGE_THRESHOLD:
        print(PRINT, f"    Rn: using measured '{measured_col}' (coverage {measured_coverage*100:.1f}% >= {COVERAGE_THRESHOLD*100}%)")
        print(PRINT, f"    Rn: NOT filling missing values (accepting {measured_coverage*100:.1f}% coverage as-is)")
        Rn = measured_Rn
    
    # Case 2: measured column exists but coverage too low (< 50%)
    elif measured_Rn is not None and measured_coverage < COVERAGE_THRESHOLD:
        print(PRINT, f"    Rn: measured '{measured_col}' has LOW coverage ({measured_coverage*100:.1f}% < {COVERAGE_THRESHOLD*100}%)")
        print(PRINT, f"    Rn: IGNORING measured Rn for this site")
        
        if fao_Rn is not None:
            print(PRINT, f"    Rn: using FAO-56 instead")
            Rn = fao_Rn
        else:
            print(PRINT, f"    Rn: no FAO-56 available, falling back to 0.45*Rg_f")
            Rg = _to_num(df.get("Rg_f"))
            if Rg is not None and Rg.notna().any():
                Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            else:
                Rn = pd.Series(np.nan, index=idx, name="Rn")
    
    # Case 3: no measured column exists
    elif fao_Rn is not None:
        print(PRINT, f"    Rn: no measured Rn, using FAO-56")
        Rn = fao_Rn
    
    # Case 4: nothing worked, final fallback
    else:
        Rg = _to_num(df.get("Rg_f"))
        if Rg is not None and Rg.notna().any():
            Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            print(PRINT, f"    Rn: using fallback 0.45*Rg_f")
        else:
            Rn = pd.Series(np.nan, index=idx, name="Rn")
            print(PRINT, f"    Rn: WARNING - no radiation data available for {site}")
    
    # ================================================================
    # Step 4: Final diagnostic
    # ================================================================
    final_count = Rn.notna().sum()
    final_coverage = 100 * final_count / total_rows
    print(PRINT, f"    Rn: FINAL valid = {final_count}/{total_rows} ({final_coverage:.1f}%)")
    
    return Rn.rename("Rn")
    """
    Get net radiation with coverage-based intelligent fallback.
    
    Strategy:
        1. Check measured Rn columns in order: NETRAD_f, NETRAD, Rn, Rn_f
        2. Calculate measured Rn coverage = valid measured values / total rows
        3. If measured Rn coverage >= 0.50:
           - use measured Rn where available
           - fill missing measured Rn values with FAO-56 Rn from Rg_f
        4. If measured Rn coverage < 0.50:
           - ignore measured Rn for this site (too sparse)
           - use FAO-56 Rn from Rg_f for the whole site
        5. If no measured Rn exists:
           - use FAO-56 Rn from Rg_f
        6. If FAO-56 fails but Rg_f exists:
           - use Rn = 0.45 * Rg_f as final fallback
        7. Return complete Series named "Rn" whenever possible
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input data frame with meteorological columns
    site : str
        Site name (for metadata lookup)
    meta : dict
        Site metadata containing latitude and elevation
    
    Returns
    -------
    Rn : pandas.Series
        Net radiation (W m-2)
    """
    idx = df.index
    total_rows = len(idx)
    lat = meta.get(site, {}).get("lat", None)
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    
    # ================================================================
    # Step 1: Check for measured Rn columns
    # ================================================================
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
                print(PRINT, f"    Rn: found measured column '{c}' for {site} (coverage = {measured_coverage*100:.1f}%)")
                break
    
    # ================================================================
    # Step 2: Get FAO-56 Rn from Rg_f if possible
    # ================================================================
    fao_Rn = None
    if lat is not None and "Rg_f" in df.columns:
        try:
            fao_Rn = rn_available_fao56(df, float(lat), float(elev))
            fao_coverage = fao_Rn.notna().sum() / total_rows
            print(PRINT, f"    Rn: FAO-56 estimate available for {site} (coverage = {fao_coverage*100:.1f}%)")
        except Exception as e:
            print(PRINT, f"    Rn: FAO-56 failed for {site}: {e}")
            fao_Rn = None
    
    # ================================================================
    # Step 3: Decision based on measured coverage
    # ================================================================
    COVERAGE_THRESHOLD = 0.50
    
    # Case 1: measured column exists with good coverage (>= 50%)
    if measured_Rn is not None and measured_coverage >= COVERAGE_THRESHOLD:
        measured_count = measured_Rn.notna().sum()
        print(PRINT, f"    Rn: using measured '{measured_col}' (coverage {measured_coverage*100:.1f}% >= {COVERAGE_THRESHOLD*100}%)")
        
        if fao_Rn is not None:
            # Fill missing measured values with FAO-56
            missing_count = (~measured_Rn.notna() & fao_Rn.notna()).sum()
            Rn = measured_Rn.combine_first(fao_Rn)
            print(PRINT, f"    Rn: filled {missing_count} missing values with FAO-56")
            print(PRINT, f"    Rn: final coverage = {Rn.notna().sum()/total_rows*100:.1f}%")
        else:
            Rn = measured_Rn
            print(PRINT, f"    Rn: no FAO-56 available, using measured only")
    
    # Case 2: measured column exists but coverage too low (< 50%)
    elif measured_Rn is not None and measured_coverage < COVERAGE_THRESHOLD:
        print(PRINT, f"    Rn: measured '{measured_col}' has LOW coverage ({measured_coverage*100:.1f}% < {COVERAGE_THRESHOLD*100}%)")
        print(PRINT, f"    Rn: IGNORING measured Rn for this site")
        
        if fao_Rn is not None:
            print(PRINT, f"    Rn: using FAO-56 instead (coverage = {fao_Rn.notna().sum()/total_rows*100:.1f}%)")
            Rn = fao_Rn
        else:
            print(PRINT, f"    Rn: no FAO-56 available, falling back to 0.45*Rg_f")
            Rg = _to_num(df.get("Rg_f"))
            if Rg is not None and Rg.notna().any():
                Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            else:
                Rn = pd.Series(np.nan, index=idx, name="Rn")
    
    # Case 3: no measured column exists
    elif fao_Rn is not None:
        print(PRINT, f"    Rn: no measured Rn, using FAO-56 (coverage = {fao_Rn.notna().sum()/total_rows*100:.1f}%)")
        Rn = fao_Rn
    
    # Case 4: nothing worked, final fallback
    else:
        Rg = _to_num(df.get("Rg_f"))
        if Rg is not None and Rg.notna().any():
            Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            print(PRINT, f"    Rn: using fallback 0.45*Rg_f (coverage = {Rn.notna().sum()/total_rows*100:.1f}%)")
        else:
            Rn = pd.Series(np.nan, index=idx, name="Rn")
            print(PRINT, f"    Rn: WARNING - no radiation data available for {site}")
    
    # ================================================================
    # Step 4: Final diagnostic and return
    # ================================================================
    final_count = Rn.notna().sum()
    final_coverage = 100 * final_count / total_rows
    print(PRINT, f"    Rn: FINAL valid = {final_count}/{total_rows} ({final_coverage:.1f}%)")
    
    return Rn.rename("Rn")
    """
    Get net radiation with intelligent fallback.
    
    Strategy:
        1. Look for measured Rn columns in order: NETRAD_f, NETRAD, Rn, Rn_f
        2. If measured Rn exists, keep its valid values
        3. Fill missing measured Rn values with FAO-56 estimate from Rg_f
        4. If no measured Rn exists, use FAO-56 Rn from Rg_f
        5. Final fallback: Rn = 0.45 * Rg_f if Rg_f exists
        6. Return complete Series named "Rn" whenever possible
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input data frame with meteorological columns
    site : str
        Site name (for metadata lookup)
    meta : dict
        Site metadata containing latitude and elevation
    
    Returns
    -------
    Rn : pandas.Series
        Net radiation (W m-2)
    """
    idx = df.index
    lat = meta.get(site, {}).get("lat", None)
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    
    # ================================================================
    # Step 1: Check for measured Rn columns
    # ================================================================
    measured_col = None
    measured_Rn = None
    for c in ("NETRAD_f", "NETRAD", "Rn", "Rn_f"):
        if c in df.columns:
            s = _to_num(df[c])
            if s.notna().any():
                measured_col = c
                measured_Rn = s.clip(lower=0.0)
                print(PRINT, f"    Rn: found measured column '{c}' for {site}")
                break
    
    # ================================================================
    # Step 2: Get FAO-56 Rn from Rg_f if possible
    # ================================================================
    fao_Rn = None
    if lat is not None and "Rg_f" in df.columns:
        try:
            fao_Rn = rn_available_fao56(df, float(lat), float(elev))
            print(PRINT, f"    Rn: FAO-56 estimate available for {site} ({fao_Rn.notna().sum()} valid values)")
        except Exception as e:
            print(PRINT, f"    Rn: FAO-56 failed for {site}: {e}")
            fao_Rn = None
    
    # ================================================================
    # Step 3: Special case - measured column exists with values
    # ================================================================
    if measured_Rn is not None:
        measured_count = measured_Rn.notna().sum()
        print(PRINT, f"    Rn: using measured '{measured_col}' ({measured_count} valid values)")
        
        if fao_Rn is not None:
            # Fill missing measured values with FAO-56
            filled_count = (~measured_Rn.notna() & fao_Rn.notna()).sum()
            Rn = measured_Rn.combine_first(fao_Rn)
            print(PRINT, f"    Rn: filled {filled_count} missing values with FAO-56")
        else:
            Rn = measured_Rn
            print(PRINT, f"    Rn: no FAO-56 available, using measured only")
    
    # ================================================================
    # Step 4: No measured Rn, but FAO-56 available
    # ================================================================
    elif fao_Rn is not None:
        print(PRINT, f"    Rn: no measured Rn, using FAO-56 ({fao_Rn.notna().sum()} valid values)")
        Rn = fao_Rn
    
    # ================================================================
    # Step 5: Final fallback - simple Rg-based estimate
    # ================================================================
    else:
        Rg = _to_num(df.get("Rg_f"))
        if Rg is not None and Rg.notna().any():
            Rn = (0.45 * Rg).clip(lower=0.0).rename("Rn")
            print(PRINT, f"    Rn: using fallback 0.45*Rg_f ({Rn.notna().sum()} valid values)")
        else:
            Rn = pd.Series(np.nan, index=idx, name="Rn")
            print(PRINT, f"    Rn: WARNING - no radiation data available for {site}")
    
    # ================================================================
    # Step 6: Final diagnostic and return
    # ================================================================
    final_count = Rn.notna().sum()
    total_count = len(Rn)
    print(PRINT, f"    Rn: final valid = {final_count}/{total_count} ({100*final_count/total_count:.1f}%)")
    
    return Rn.rename("Rn")


def guess_ra(idx, h_m, WS, ra_floor=15.0, kB1=2.3):
    """Aerodynamic resistance ra [s m-1] (log-law with kB^-1 for heat)."""
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
    """
    Beer's law partitioning of radiation between canopy and soil.
    
    PHYSICAL COMPONENT: Radiation partitioning.
    """
    L = pd.Series(0.0, index=idx) if LAI is None else _to_num(LAI).reindex(idx)
    L = L.fillna(0.0).clip(lower=0.0)
    f_c = (1.0 - np.exp(-float(k_beer) * L)).clip(0.0, 1.0)
    return f_c.rename("f_canopy"), (1.0 - f_c).rename("f_soil")


def soil_resistance(LAI, idx, rsoil_min=150.0, rsoil_max=600.0, k_lai=0.5):
    """
    Soil surface resistance as function of LAI.
    
    EMPIRICAL / PARAMETERIZED COMPONENT: Soil resistance parameterization.
    """
    L = pd.Series(0.0, index=idx) if LAI is None else _to_num(LAI).reindex(idx)
    L = L.fillna(0.0).clip(lower=0.0)
    rs = float(rsoil_min) * np.exp(float(k_lai) * L)
    return pd.Series(rs, index=idx, name="rs_soil").clip(rsoil_min, rsoil_max)


def std_pressure_from_elev_kpa(elev_m):
    """Standard atmospheric pressure from elevation (kPa)."""
    z = float(elev_m) if pd.notna(elev_m) else 0.0
    return 101.325 * (1.0 - 2.25577e-5 * z) ** 5.25588


def pressure_for_equation(site, df, meta):
    """
    Get air pressure for calculations.
    
    Priority: PA_f column > PA column > elevation-based estimate.
    """
    idx = df.index
    
    # Check for PA_f first (standard AmeriFlux naming)
    P = _to_num(df.get("PA_f"))
    if P is not None and (P > 0).any():
        elev = meta.get(site, {}).get("elevation_m", 0.0)
        Pfallback = std_pressure_from_elev_kpa(elev)
        return pd.Series(P, index=idx).where(P > 0, Pfallback).rename("PA_kPa")
    
    # Check for PA (alternative naming)
    P = _to_num(df.get("PA"))
    if P is not None and (P > 0).any():
        elev = meta.get(site, {}).get("elevation_m", 0.0)
        Pfallback = std_pressure_from_elev_kpa(elev)
        return pd.Series(P, index=idx).where(P > 0, Pfallback).rename("PA_kPa")
    
    # Fallback to elevation-based pressure
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    return pd.Series(std_pressure_from_elev_kpa(elev), index=idx, name="PA_kPa")


def wet_mask_from(df, site, params, meta):
    """
    Determine wet periods based on precipitation.
    
    EMPIRICAL / PARAMETERIZED COMPONENT: Precipitation-based wet mask.
    """
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
    """Parse canopy height from range string like '0.10–0.50' or return float."""
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


# ============================================================================
# MAIN PARTITIONING FUNCTION
# ============================================================================

def pm_partition(df, site, h_canopy_m, params, meta):
    """
    Penman-Monteith based partitioning: E_est followed by T = ET_obs - E_est.
    
    CORE EQUATIONS:
        Evaporation (E_est) computed using Penman-Monteith formulation:
            E = [Δ(Rn - G) + ρ·Cp·(VPD/ra)] / [Δ + γ·(1 + rs/ra)]
        
        Transpiration is NOT modeled directly:
            T = ET_obs - E_est
        
        This is a residual partitioning approach.
    
    PHYSICAL COMPONENTS:
        - Penman-Monteith evaporation equation
        - Aerodynamic resistance (ra) via log-law with kB1 for heat
        - Radiation partitioning via Beer-Lambert law (f_c, f_s)
        - Net radiation (Rn) from measured or FAO-56
    
    EMPIRICAL / PARAMETERIZED COMPONENTS:
        - Soil resistance (rs_soil) as exponential function of LAI
        - Wet canopy resistance (rs_wet_canopy) fixed per biome
        - Wet-period evaporation cap (wet_E_cap_frac) - regularization constraint
        - Wet mask (precipitation-based proxy, not soil moisture)
        - Nighttime transpiration fallback (night_T_frac) - conservative lower bound
          from literature (Wu et al. 2023: 5.5-24% in mangroves)
    
    UNCERTAINTY NOTE:
        Radiation partitioning depends on LAI. Uncertainty in LAI propagates
        into both evaporation and transpiration estimates, as it controls the
        distribution of energy between canopy and soil surfaces.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input data with datetime index and required variables
    site : str
        Site name
    h_canopy_m : float
        Canopy height in meters
    params : dict
        Parameter set for this biome
    meta : dict
        Site metadata
    
    Returns
    -------
    Trans_pen, Evap_pen : pandas.Series
        Transpiration and evaporation fluxes
    """
    idx = df.index
    
    # Get observed ET
    ET = _to_num(df.get("ET"))
    if ET is None or ET.isna().all():
        print(PRINT, f"  Warning: No ET data for {site}, skipping partitioning")
        return pd.Series(np.nan, index=idx, name="Trans_pen"), pd.Series(np.nan, index=idx, name="Evap_pen")
    
    # Meteorological drivers
    Tair = _to_num(df.get("Tair_f"))
    VPD_h = _to_num(df.get("VPD_f"))
    WS = _to_num(df.get("WS"))
    
    # Leaf area index
    LAI = _to_num(df.get("lai"))
    if LAI is None or LAI.isna().all():
        LAI = _to_num(df.get("LAI_daily"))
    
    # Time step conversion (mm per half-hour from W m-2)
    dt_s = _median_dt_seconds(idx, 1800.0)
    conv = dt_s / 2.45e6  # converts W m-2 to mm per timestep
    
    # Psychrometric constants
    PkPa = pressure_for_equation(site, df, meta)
    gamma = 0.000665 * PkPa  # kPa K-1
    
    # Vapor pressure deficit
    VPD_kP = (VPD_h / 10.0) if VPD_h is not None else pd.Series(np.nan, index=idx)
    VPD_kP = VPD_kP.clip(lower=0.0)
    
    # Slope of saturation vapor pressure curve
    if Tair is not None and Tair.notna().any():
        es = 0.6108 * np.exp(17.27 * Tair / (Tair + 237.3))
        delta = 4098.0 * es / ((Tair + 237.3) ** 2)
    else:
        delta = pd.Series(0.15, index=idx)
    
    # Air density
    R = 8.314462618
    M_air = 0.029
    Cp = 1010.0
    rho = ((PkPa * 1000.0) / (R * (Tair + 273.15))) * M_air
    
    # Net radiation
    Rn = rn_available(df, site, meta)
    
    # PHYSICAL COMPONENT: Radiation partitioning via Beer-Lambert
    f_c, f_s = beer_partition(LAI, idx, params["k_beer"])
    Rnc = f_c * Rn    # canopy net radiation
    Rns = f_s * Rn    # soil net radiation
    
    # Soil heat flux
    G = (float(params["G_frac"]) * Rns).clip(upper=0.30 * Rns)
    Rns_eff = (Rns - G).clip(lower=0.0)
    
    # PHYSICAL COMPONENT: Aerodynamic resistance
    ra = guess_ra(idx, h_canopy_m, WS, params["ra_floor"], params["kB1"])
    
    # EMPIRICAL COMPONENT: Soil and wet canopy resistances
    rs_soil = soil_resistance(LAI, idx, params["rsoil_min"], params["rsoil_max"], 0.5)
    rs_wet = pd.Series(float(params["rs_wet_canopy"]), index=idx)
    
    # EMPIRICAL COMPONENT: Wet mask (precipitation-based)
    wet = wet_mask_from(df, site, params, meta)
    
    # Optional diagnostic storage
    if SAVE_DIAGNOSTICS:
        diag_dir = Path(output_dir) / 'diagnostics' / site
        diag_dir.mkdir(parents=True, exist_ok=True)
        # Would save wet_mask, ra, rs_soil, f_c here if implemented
    
    # Initialize evaporation array
    Evap = pd.Series(np.nan, index=idx, name="Evap_pen")
    
    # ================================================================
    # PART 1: Wet canopy evaporation (when canopy is wet)
    # ================================================================
    m_wet = wet & Rnc.notna() & (Rnc > 0)
    if m_wet.any():
        # Penman-Monteith equation for wet canopy
        num_energy = delta[m_wet] * Rnc[m_wet] * conv
        num_aero = (rho[m_wet] * Cp * VPD_kP[m_wet] / ra[m_wet]) * conv
        denom = (delta[m_wet] + gamma[m_wet] * (1.0 + rs_wet[m_wet] / ra[m_wet])).replace(0, np.nan)
        Evap.loc[m_wet] = ((num_energy + num_aero) / denom).clip(lower=0.0)
        
        # EMPIRICAL REGULARIZATION CONSTRAINT:
        # During wet periods, evaporation is limited to a fraction of ET.
        # This acts as a regularization term in the residual partitioning framework,
        # preventing the model from attributing all ET to wet-surface evaporation,
        # which would suppress transpiration estimates due to the residual formulation.
        # This constraint is particularly important in coastal wetlands where
        # shallow water and wet surfaces can otherwise dominate the PM estimate.
        cap = float(params.get("wet_E_cap_frac", 0.65))
        for i in idx[m_wet]:
            if pd.notna(ET.loc[i]) and ET.loc[i] > 0:
                Evap.loc[i] = min(Evap.loc[i], cap * ET.loc[i])
    
    # ================================================================
    # PART 2: Dry soil evaporation (when canopy is dry)
    # ================================================================
    m_dry = (~wet) & Rns_eff.notna() & (Rns_eff > 0)
    if m_dry.any():
        # Penman-Monteith equation for dry soil
        num_energy = delta[m_dry] * Rns_eff[m_dry] * conv
        num_aero = (rho[m_dry] * Cp * VPD_kP[m_dry] / ra[m_dry]) * conv
        denom = (delta[m_dry] + gamma[m_dry] * (1.0 + rs_soil[m_dry] / ra[m_dry])).replace(0, np.nan)
        Evap.loc[m_dry] = ((num_energy + num_aero) / denom).clip(lower=0.0)
    
    # ================================================================
    # POST-PROCESSING: Clipping and residual transpiration
    # ================================================================
    ET_pos = ET.clip(lower=0)
    
    # Only calculate where ET > 0
    pos = ET_pos.notna() & (ET_pos > 0)
    Evap.loc[~pos] = np.nan
    
    # Physical constraint: Evaporation cannot exceed total ET
    Evap = Evap.clip(lower=0.0)
    Evap = Evap.where(ET_pos.notna(), np.nan).clip(upper=ET_pos)
    
    # ================================================================
    # FALLBACK FOR MISSING EVAPORATION ESTIMATES (e.g., nighttime / low radiation)
    # ================================================================
    # When Evap is NaN but ET > 0, the residual Trans becomes NaN.
    # This occurs primarily at night or under low-radiation conditions
    # where the PM equations may not produce valid E estimates.
    # Based on literature (Wu et al. 2023, Frontiers in Plant Science),
    # nocturnal transpiration in coastal wetlands ranges from 5.5% to 24% of daily total.
    # We apply a conservative lower bound (5% for wetlands) as a fallback.
    
    night_T_frac = float(params.get("night_T_frac", 0.05))
    missing_evap = (ET_pos > 0) & Evap.isna()
    
    if missing_evap.any():
        Trans_fallback = night_T_frac * ET_pos
        Evap.loc[missing_evap] = (
            ET_pos.loc[missing_evap] - Trans_fallback.loc[missing_evap]
        ).clip(lower=0.0)
    
    # RESIDUAL TRANSPIRATION: T = ET_obs - E_est
    Trans = (ET_pos - Evap).where(ET_pos.notna(), np.nan).clip(lower=0.0)
    
    return Trans.rename("Trans_pen"), Evap.rename("Evap_pen")


# ============================================================================
# SITE PROCESSING FUNCTIONS
# ============================================================================

def load_site_meta(site_info_csv):
    """Load site metadata from CSV file with proper encoding handling."""
    encodings = ['utf-8', 'latin1', 'cp1252', 'iso-8859-1', 'windows-1252']
    
    info = None
    for encoding in encodings:
        try:
            info = pd.read_csv(site_info_csv, encoding=encoding)
            print(PRINT, f"Successfully read CSV with encoding: {encoding}")
            break
        except UnicodeDecodeError:
            continue
        except Exception as e:
            print(PRINT, f"Error with encoding {encoding}: {e}")
            continue
    
    if info is None:
        raise ValueError(f"Could not read {site_info_csv} with any common encoding")
    
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
    print(PRINT, f"Loaded site metadata: {len(info)} sites")
    return meta


def _params_for_site(site: str, meta: dict):
    """Get parameters, canopy height, and biome for a site."""
    site_info = meta.get(site, {})
    biome = site_info.get("biome", "WET")
    h_canopy = float(site_info.get("canopy_height_m", 0.3) or 0.3)
    
    params = biome_defaults(biome)
    
    required = ("k_beer", "G_frac", "ra_floor", "kB1", "rsoil_min", "rsoil_max", "rs_wet_canopy")
    for k in required:
        if k not in params:
            raise KeyError(f"params missing '{k}' for site {site}")
    
    return params, h_canopy, biome


# Desired output column order - UNCHANGED
DESIRED_ORDER = [
    "DateTime", "Year", "month", "day", "DOY", "Hour",
    "NEE", "GPP", "Reco", "NEP", "ET",
    "NEE_f", "LE_f", "H_f", "Tair_f", "PA", "RH", "WS", "WD", "VPD_f", "Rg_f",
    "PAR_f", "NETRAD_f", "lambda", "precip_mm_per_30min", "lai",
    "avg_sos", "avg_eos",
    "Evap_pen", "Trans_pen"
]


def normalize_datetime_column(df):
    """Ensure DateTime is a proper column with datetime64 type."""
    df_clean = df.copy()
    df_clean = df_clean.reset_index(drop=True)
    
    bad = [c for c in df_clean.columns if str(c).startswith("Unnamed:")]
    if bad:
        df_clean = df_clean.drop(columns=bad)
    
    if "DateTime" not in df_clean.columns:
        datetime_cols = [col for col in df_clean.columns if 'date' in col.lower() or 'time' in col.lower()]
        if datetime_cols:
            df_clean = df_clean.rename(columns={datetime_cols[0]: "DateTime"})
        else:
            df_clean.insert(0, "DateTime", pd.NaT)
    
    df_clean["DateTime"] = pd.to_datetime(df_clean["DateTime"], errors="coerce")
    df_clean = df_clean.dropna(subset=["DateTime"])
    df_clean.index.name = None
    
    return df_clean


def run_site_partitioning(df, site, meta):
    """Run PM partitioning for one site."""
    params, h_canopy, biome = _params_for_site(site, meta)
    print(PRINT, f"  Site: {site} | Biome: {biome} | Canopy height: {h_canopy:.2f}m | Elevation: {meta[site]['elevation_m']:.1f}m")
    
    df_temp = df.copy()
    df_temp = df_temp.reset_index(drop=True)
    
    if "DateTime" not in df_temp.columns:
        raise ValueError(f"DateTime column missing for {site}")
    
    df_temp["DateTime"] = pd.to_datetime(df_temp["DateTime"], errors="coerce")
    df_temp = df_temp.dropna(subset=["DateTime"])
    
    df_temp = df_temp.set_index("DateTime")
    
    Trans_pen, Evap_pen = pm_partition(df_temp, site, h_canopy, params, meta)
    
    result = pd.DataFrame({
        "DateTime": df_temp.index,
        "Trans_pen": Trans_pen,
        "Evap_pen": Evap_pen
    }).reset_index(drop=True)
    
    return result


def process_and_save_all(input_dir, output_dir, site_info_csv, run_diagnostics=False):
    """Main function to process all site files."""
    meta = load_site_meta(site_info_csv)
    
    # Create output directory FIRST
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # THEN export parameter table (after directory exists)
    export_parameters_to_csv(output_dir)
    
    files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    print(PRINT, f"Found {len(files)} files | output → {output_dir}")
    
    success_count = 0
    all_annual_data = []
    
    for f in files:
        site = Path(f).stem
        try:
            print(f"\n{PRINT} Processing {site}...")
            raw = pd.read_csv(f, engine="python")
            print(f"  Original rows: {len(raw)}")
            
            raw = normalize_datetime_column(raw)
            print(f"  After normalization: {len(raw)} rows")
            
            required_cols = ["ET", "Tair_f", "VPD_f"]
            missing = [c for c in required_cols if c not in raw.columns]
            
            if "lai" not in raw.columns and "LAI_daily" not in raw.columns:
                missing.append("lai or LAI_daily")
                print(f"  WARNING: No LAI column found (looking for 'lai' or 'LAI_daily')")
            else:
                if "LAI_daily" in raw.columns and "lai" not in raw.columns:
                    raw = raw.rename(columns={"LAI_daily": "lai"})
                    print(f"  Renamed 'LAI_daily' to 'lai'")
            
            if missing:
                print(f"  SKIPPING {site} - missing required data: {missing}")
                continue
            
            part = run_site_partitioning(raw.copy(), site, meta)
            
            for old in ["LE_tr", "Trans_dir", "Evap_pen", "Trans_pen"]:
                if old in raw.columns:
                    raw = raw.drop(columns=[old])
            
            merged = pd.merge(
                raw,
                part[["DateTime", "Evap_pen", "Trans_pen"]],
                on="DateTime",
                how="left",
                validate="m:1"
            )
            
            for col in DESIRED_ORDER:
                if col not in merged.columns:
                    merged[col] = pd.NA
            
            merged = merged[DESIRED_ORDER]
            
            # OUTPUT PATH AND FILENAME UNCHANGED
            out_path = os.path.join(output_dir, f"{site}.csv")
            merged.to_csv(out_path, index=False)
            print(f"  ✓ Saved → {os.path.basename(out_path)}")
            print(f"    Evap_pen valid: {merged['Evap_pen'].notna().sum()}/{len(merged)}")
            print(f"    Trans_pen valid: {merged['Trans_pen'].notna().sum()}/{len(merged)}")
            if 'avg_sos' in merged.columns:
                print(f"    avg_sos present: {merged['avg_sos'].notna().sum()} non-null values")
            if 'avg_eos' in merged.columns:
                print(f"    avg_eos present: {merged['avg_eos'].notna().sum()} non-null values")
            
            if run_diagnostics:
                if 'Year' in merged.columns and 'Evap_pen' in merged.columns and 'Trans_pen' in merged.columns:
                    annual = merged.groupby('Year').agg({
                        'Evap_pen': 'sum',
                        'Trans_pen': 'sum'
                    }).reset_index()
                    annual['site'] = site
                    annual['IGBP'] = meta.get(site, {}).get('biome', 'unknown')
                    annual['Trans_Fraction'] = annual['Trans_pen'] / (annual['Trans_pen'] + annual['Evap_pen'])
                    all_annual_data.append(annual)
            
            success_count += 1
            
        except Exception as e:
            print(f"  ✗ ERROR processing {site}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{PRINT} {'='*50}")
    print(f"Completed! Successfully processed {success_count}/{len(files)} files.")
    print(f"{'='*50}")
    
    if run_diagnostics and all_annual_data:
        annual_df = pd.concat(all_annual_data, ignore_index=True)
        create_diagnostic_report(annual_df, output_dir)
        print(PRINT, f"Diagnostic report saved to {output_dir}/diagnostics/")


if __name__ == "__main__":
    # INPUT AND OUTPUT PATHS - UNCHANGED
    site_info_csv = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
    input_dir = r"M:\Research\WUE_CUE\ameri_data\ameri_fill__lai_precip"
    output_dir = r"M:\Research\WUE_CUE\ameri_data\ET_partitioning"
    
    print(PRINT, "Starting ET partitioning: PM-based evaporation + residual transpiration")
    print(PRINT, "Method: E_est = PM(wet canopy or dry soil); T = ET_obs - E_est")
    print(PRINT, f"Input: {input_dir}")
    print(PRINT, f"Output: {output_dir}")
    print(PRINT, "\nParameter sets applied:")
    print(PRINT, f"  WET: wet_cap=0.60, night_T_frac=0.05, rsoil_min=250, k_beer=0.48")
    print(PRINT, f"  Shrub/Grass: wet_cap=0.50, night_T_frac=0.05, rsoil_min=250, k_beer=0.58")
    print(PRINT, f"  BSV: wet_cap=0.45, night_T_frac=0.03, rsoil_min=350")
    print(PRINT, f"  Forest: wet_cap=0.62, night_T_frac=0.08, rsoil_min=140")
    print(PRINT, "\nConstraints applied:")
    print(PRINT, "  - Wet-period evaporation cap (empirical regularization for coastal wetlands)")
    print(PRINT, "  - Evaporation <= ET (physical constraint)")
    print(PRINT, "  - Nighttime fallback (night_T_frac) for missing Evap estimates")
    print(PRINT, "    Literature basis: Wu et al. 2023 (5.5-24% in mangroves)")
    print(PRINT, "\nUncertainty note: LAI uncertainty propagates into both E and T")
    
    process_and_save_all(input_dir, output_dir, site_info_csv, run_diagnostics=True)
    print(PRINT, "DONE.")