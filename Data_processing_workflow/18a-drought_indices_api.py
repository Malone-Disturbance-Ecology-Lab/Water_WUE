# -*- coding: utf-8 -*-
"""
Daily PET via Penman–Monteith (site-aware) - CRITICALLY REVIEWED VERSION
- Uses IGBP classification + canopy height for aerodynamic and canopy resistance
- FIXED: Canopy resistance scaling (rc decreases with higher LAI)
- FIXED: En dash parsing for canopy height ranges
- FIXED: IGBP grouping with proper fallbacks
- FIXED: Encoding issues with metadata CSV
- FIXED: Duplicate DateTime handling in merge
- Handles missing LAI_daily by using reference LAI values
"""

import os, glob
from pathlib import Path
import numpy as np
import pandas as pd
import warnings

PRINT = "[PET]"

# ---------------------------
# Small utilities
# ---------------------------
def _to_num(x):
    return pd.to_numeric(x, errors="coerce") if x is not None else None

def _median_dt_seconds(index, fallback=86400.0):
    """Robust modal timestep in seconds; fallback to 86400 (daily) if unclear."""
    try:
        di = pd.Index(pd.to_datetime(index)).diff().total_seconds()
        mode = pd.Series(di).dropna().mode()
        return float(mode.iloc[0]) if len(mode) else float(fallback)
    except Exception:
        return float(fallback)

def lambda_Jkg(Tair_C):
    """Latent heat of vaporization (J kg-1) from °C."""
    return (2.501 - 0.002361 * pd.to_numeric(Tair_C, errors="coerce")) * 1e6

# ---------------------------
# Solar/FAO-56 helpers (net radiation fallback)
# ---------------------------
def _solar_declination(doy):
    return 0.409 * np.sin(2.0 * np.pi * (doy - 81.0) / 365.0)

def _inv_rel_dist(doy):
    return 1.0 + 0.033 * np.cos(2.0 * np.pi * doy / 365.0)

def _sunset_hour_angle(lat_rad, delta):
    cosw = -np.tan(lat_rad) * np.tan(delta)
    return np.arccos(np.clip(cosw, -1.0, 1.0))

def _Ra_MJ_m2_day(lat_deg, doy):
    """Daily extraterrestrial radiation (MJ m-2 d-1)."""
    lat = np.radians(lat_deg)
    delta = _solar_declination(doy)
    dr = _inv_rel_dist(doy)
    ws = _sunset_hour_angle(lat, delta)
    Gsc = 0.0820
    return (24.0 * 60.0 / np.pi) * Gsc * dr * (
        ws*np.sin(lat)*np.sin(delta) + np.cos(lat)*np.cos(delta)*np.sin(ws)
    )

def _doy_from_index(idx):
    try:
        return pd.DatetimeIndex(idx).dayofyear
    except Exception:
        return pd.Series(182, index=idx)

def rn_available_fao56(df, lat_deg, elev_m, albedo=0.23):
    """Net radiation via FAO-56 (Allen et al., 1998)."""
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
    """Return net radiation (W m-2); prefer measured NETRAD_f else FAO-56."""
    for c in ("NETRAD_f","NETRAD","Rn","Rn_f"):
        if c in df.columns:
            s = _to_num(df[c])
            if s is not None and s.notna().any():
                return s.clip(lower=0.0).rename("Rn")
    
    lat = meta.get(site, {}).get("lat", None)
    elev = meta.get(site, {}).get("elevation_m", None)
    if lat is not None and elev is not None and "Rg_f" in df.columns:
        return rn_available_fao56(df, float(lat), float(elev))
    
    Rg = _to_num(df.get("Rg_f"))
    return (0.45*Rg).rename("Rn") if Rg is not None else pd.Series(np.nan, index=df.index, name="Rn")

# ---------------------------
# Aerodynamics / pressure
# ---------------------------
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
    elev = meta.get(site, {}).get("elevation_m", 0.0)
    return pd.Series(std_pressure_from_elev_kpa(elev), index=idx, name="PA_kPa")

def guess_ra(idx, h_m, WS, ra_floor=15.0, kB1=2.3):
    """Aerodynamic resistance ra [s m-1] using log-law with kB^-1 for heat."""
    try:
        h = float(h_m)
        if not np.isfinite(h) or h <= 0:
            h = 0.3
    except Exception:
        h = 0.3
    
    z = max(3.0*h, 2.0)
    d = 0.67*h
    z0m = 0.1*h
    z0h = 0.1*h*np.exp(-float(kB1))

    WS = _to_num(WS)
    if WS is None:
        return pd.Series(float(ra_floor), index=idx, name="ra_s_m")
    WS = pd.Series(WS, index=idx).clip(lower=0.5)

    num_m = (z - d) / max(z0m, 1e-6)
    num_h = (z - d) / max(z0h, 1e-6)
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = (np.log(num_m)*np.log(num_h)) / (0.4**2 * WS)
    return pd.Series(ra, index=idx, name="ra_s_m") \
             .replace([np.inf, -np.inf], np.nan) \
             .fillna(float(ra_floor)).clip(lower=float(ra_floor))

# ---------------------------
# IGBP-based parameter mapping
# ---------------------------
def get_params_from_igbp(igbp_code):
    """PET parameters based on IGBP classification."""
    if igbp_code is None:
        code = "UNKNOWN"
    else:
        code = str(igbp_code).strip().upper()
    
    # Forest group (ENF, DBF, MF)
    if code in ["ENF", "DBF", "MF"]:
        return {
            'G_frac': 0.06, 'ra_floor': 15.0, 'kB1': 2.3,
            'rc_ref': 90.0, 'LAI_ref': 3.0, 'rc_min': 50.0, 'rc_max': 200.0
        }
    # Grass/shrub/crop group (GRA, CSH, OSH, CRO)
    elif code in ["GRA", "CSH", "OSH", "CRO"]:
        return {
            'G_frac': 0.10, 'ra_floor': 15.0, 'kB1': 2.3,
            'rc_ref': 70.0, 'LAI_ref': 2.0, 'rc_min': 40.0, 'rc_max': 160.0
        }
    # Barren/sparse group (BSV)
    elif code == "BSV":
        return {
            'G_frac': 0.15, 'ra_floor': 20.0, 'kB1': 2.6,
            'rc_ref': 120.0, 'LAI_ref': 0.5, 'rc_min': 80.0, 'rc_max': 250.0
        }
    # Wetlands (WET)
    elif code == "WET":
        return {
            'G_frac': 0.12, 'ra_floor': 18.0, 'kB1': 2.6,
            'rc_ref': 85.0, 'LAI_ref': 2.5, 'rc_min': 45.0, 'rc_max': 180.0
        }
    # Unknown - use grass/shrub defaults with warning
    else:
        warnings.warn(f"Unknown IGBP code '{code}' - using grass/shrub defaults")
        return {
            'G_frac': 0.10, 'ra_floor': 15.0, 'kB1': 2.3,
            'rc_ref': 70.0, 'LAI_ref': 2.0, 'rc_min': 40.0, 'rc_max': 160.0
        }

def parse_canopy_height(height_value):
    """Parse canopy height from various formats."""
    if pd.isna(height_value):
        return 0.3
    
    h_str = str(height_value).strip()
    
    # Try direct numeric conversion first
    try:
        return float(h_str)
    except ValueError:
        pass
    
    # Check for range separators
    for separator in ['–', '-', '—']:
        if separator in h_str:
            parts = h_str.split(separator)
            if len(parts) >= 2:
                try:
                    low = float(parts[0].strip())
                    high = float(parts[1].strip())
                    return (low + high) / 2.0
                except ValueError:
                    continue
    
    # Extract any number
    import re
    numbers = re.findall(r'[\d.]+', h_str)
    if numbers:
        try:
            return float(numbers[0])
        except ValueError:
            pass
    
    warnings.warn(f"Could not parse canopy height '{height_value}', using 0.3m")
    return 0.3

def load_site_meta(site_info_csv):
    """Load site metadata with IGBP classification."""
    encodings = ['utf-8', 'latin1', 'cp1252', 'ISO-8859-1', 'utf-16']
    
    info = None
    for encoding in encodings:
        try:
            info = pd.read_csv(site_info_csv, encoding=encoding)
            print(PRINT, f"Successfully read metadata with {encoding} encoding")
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
        except Exception as e:
            continue
    
    if info is None:
        try:
            info = pd.read_csv(site_info_csv, encoding_errors='ignore')
        except Exception as e:
            raise Exception(f"Could not read CSV file with any encoding. Last error: {e}")
    
    print(PRINT, f"Metadata columns found: {list(info.columns)}")
    
    meta = {}
    
    for _, r in info.iterrows():
        site = str(r["site_name"]).strip()
        
        lat = float(r["lat"]) if "lat" in info.columns else np.nan
        lon = float(r["long"]) if "long" in info.columns else np.nan
        elev = float(r["elevation_m"]) if "elevation_m" in info.columns else 0.0
        igbp = str(r["IGBP"]).strip().upper() if "IGBP" in info.columns else "UNKNOWN"
        h = parse_canopy_height(r["canopy_height_m"]) if "canopy_height_m" in info.columns else 0.3
        
        meta[site] = {
            'lat': lat,
            'long': lon,
            'elevation_m': elev,
            'igbp': igbp,
            'canopy_height_m': h
        }
        
        if len(meta) <= 10:
            print(PRINT, f"Loaded site: {site}, IGBP: {igbp}, height: {h:.2f}m")
    
    print(PRINT, f"Loaded site metadata: {len(meta)} sites")
    return meta

def _params_for_site(site: str, meta: dict):
    """Get PET parameters for a site based on IGBP classification."""
    site_info = meta.get(site, {})
    
    if not site_info:
        warnings.warn(f"Site {site} not found in metadata, using defaults")
        igbp = "UNKNOWN"
        h_canopy = 0.3
    else:
        igbp = site_info.get('igbp', 'UNKNOWN')
        h_canopy = float(site_info.get('canopy_height_m', 0.3))
    
    params = get_params_from_igbp(igbp)
    
    required_keys = ['G_frac', 'ra_floor', 'kB1', 'rc_ref', 'LAI_ref', 'rc_min', 'rc_max']
    for key in required_keys:
        if key not in params:
            raise KeyError(f"Parameters missing '{key}' for site {site} (IGBP: {igbp})")
    
    return params, h_canopy, igbp

# ---------------------------
# Core PET (Penman–Monteith)
# ---------------------------
def pet_pm(df, site, h_canopy_m, params, meta):
    """Penman–Monteith PET (mm per day)."""
    idx = df.index
    
    Tair = _to_num(df.get("Tair_f"))
    VPD_h = _to_num(df.get("VPD_f"))
    WS = _to_num(df.get("WS"))
    LAI = _to_num(df.get("LAI_daily"))
    
    if Tair is None or VPD_h is None:
        print(PRINT, f"Warning: Missing Tair or VPD for {site}")
        return pd.Series(np.nan, index=idx, name="PET")
    
    dt_s = _median_dt_seconds(idx, 86400.0)
    lam = lambda_Jkg(Tair)
    conv = dt_s / lam
    
    PkPa = pressure_for_equation(site, df, meta)
    gamma = 0.000665 * PkPa
    VPD_kPa = VPD_h / 10.0
    
    es = 0.6108 * np.exp(17.27 * Tair / (Tair + 237.3))
    delta = 4098.0 * es / ((Tair + 237.3) ** 2)
    
    R_specific = 287.04
    Cp = 1004.0
    rho = (PkPa * 1000.0) / (R_specific * (Tair + 273.15))
    
    Rn = rn_available(df, site, meta)
    G = params['G_frac'] * Rn
    G = G.clip(upper=0.30 * Rn)
    Rn_eff = (Rn - G).clip(lower=0.0)
    
    ra = guess_ra(idx, h_canopy_m, WS, params['ra_floor'], params['kB1'])
    
    if LAI is None or LAI.isna().all():
        rc = pd.Series(params['rc_ref'], index=idx)
    else:
        LAI_eff = pd.Series(LAI, index=idx).fillna(params['LAI_ref']).clip(lower=0.5)
        rc = params['rc_ref'] * (params['LAI_ref'] / LAI_eff)
        rc = rc.clip(params['rc_min'], params['rc_max'])
    
    num_energy = delta * Rn_eff
    num_aero = rho * Cp * VPD_kPa / ra
    denom = delta + gamma * (1.0 + rc / ra)
    denom = denom.replace(0, np.nan)
    
    PET = ((num_energy + num_aero) / denom) * conv
    return pd.Series(PET, index=idx, name="PET").clip(lower=0.0)

# ---------------------------
# I/O helpers and driver
# ---------------------------
def _normalize_datetime_column(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure DateTime is a clean column and remove duplicates."""
    df_clean = df.copy()
    df_clean = df_clean.reset_index(drop=True)
    
    # Remove unnamed columns
    bad = [c for c in df_clean.columns if str(c).startswith("Unnamed:")]
    if bad:
        df_clean = df_clean.drop(columns=bad)
    
    # Handle duplicate columns
    if df_clean.columns.duplicated().any():
        df_clean = df_clean.loc[:, ~df_clean.columns.duplicated()]
    
    # Find or create DateTime column
    if "DateTime" not in df_clean.columns:
        candidates = [c for c in df_clean.columns if "date" in c.lower() or "time" in c.lower()]
        if candidates:
            df_clean = df_clean.rename(columns={candidates[0]: "DateTime"})
        else:
            df_clean.insert(0, "DateTime", pd.NaT)
    
    df_clean["DateTime"] = pd.to_datetime(df_clean["DateTime"], errors="coerce")
    df_clean = df_clean.dropna(subset=["DateTime"])
    
    # CRITICAL FIX: Remove duplicate DateTime rows
    # Keep first occurrence of each timestamp
    df_clean = df_clean.drop_duplicates(subset=["DateTime"], keep='first')
    
    # Sort by DateTime
    df_clean = df_clean.sort_values("DateTime")
    
    df_clean.index.name = None
    return df_clean

def run_site_pet(df, site, meta):
    """Compute PET for one site."""
    params, h_canopy, igbp = _params_for_site(site, meta)
    
    if "DateTime" not in df.columns:
        raise ValueError("DateTime column missing")
    
    print(PRINT, f"Site {site}: IGBP={igbp}, height={h_canopy:.2f}m, rc_ref={params['rc_ref']}")
    
    df_i = df.set_index("DateTime")
    df_i.index = pd.to_datetime(df_i.index)
    
    # Check for duplicate indices and remove them
    if df_i.index.duplicated().any():
        print(PRINT, f"  Warning: Found {df_i.index.duplicated().sum()} duplicate timestamps, removing duplicates")
        df_i = df_i[~df_i.index.duplicated(keep='first')]
    
    PET = pet_pm(df_i, site, h_canopy, params, meta)
    
    out = pd.DataFrame({"DateTime": df_i.index, "PET": PET}).reset_index(drop=True)
    return out

def process_and_save_all(input_dir, output_dir, site_info_csv):
    """Process all site files and save with PET column."""
    meta = load_site_meta(site_info_csv)
    
    # Look for CSV files (both with and without US- prefix)
    files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    if not files:
        print(PRINT, f"Warning: No CSV files found in {input_dir}")
        return
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(PRINT, f"Found {len(files)} files | output → {output_dir}")
    
    for f in files:
        site = Path(f).stem
        try:
            print(PRINT, f"Reading {f}...")
            raw = pd.read_csv(f, engine='python')
            raw = _normalize_datetime_column(raw)
            
            if len(raw) == 0:
                print(PRINT, f"  Warning: No valid data after cleaning for {site}")
                continue
            
            print(PRINT, f"Processing {site} | rows={len(raw)}")
            pet_df = run_site_pet(raw.copy(), site, meta)
            
            # Check for duplicates in pet_df before merging
            if pet_df['DateTime'].duplicated().any():
                print(PRINT, f"  Warning: PET DataFrame has duplicates, removing them")
                pet_df = pet_df.drop_duplicates(subset=["DateTime"], keep='first')
            
            # Merge PET back - using left join without validation to avoid errors
            # Then check for issues manually
            merged = pd.merge(raw, pet_df, on="DateTime", how="left")
            
            # Check if merge created duplicate rows
            if len(merged) != len(raw):
                print(PRINT, f"  Warning: Merge changed row count from {len(raw)} to {len(merged)}")
                # Try alternative merge strategy
                merged = raw.copy()
                pet_dict = dict(zip(pet_df['DateTime'], pet_df['PET']))
                merged['PET'] = merged['DateTime'].map(pet_dict)
            
            out_path = os.path.join(output_dir, f"{site}.csv")
            merged.to_csv(out_path, index=False)
            
            if 'PET' in merged.columns and not merged['PET'].isna().all():
                print(PRINT, f"  Saved → {os.path.basename(out_path)} (PET range: {merged['PET'].min():.2f}-{merged['PET'].max():.2f} mm/day)")
            else:
                print(PRINT, f"  Saved → {os.path.basename(out_path)} (PET all NaN - check input data)")
            
        except Exception as e:
            print(PRINT, f"[ERROR] {site}: {e}")
            import traceback
            traceback.print_exc()

# ---------------------------
# Main
# ---------------------------
if __name__ == "__main__":
    site_info_csv = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_lat_long.csv"
    input_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\input_data"
    output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\PET_drought"
    
    print(PRINT, "Starting PET writeout (daily) - CRITICALLY REVIEWED VERSION")
    print(PRINT, "Key fixes: IGBP classification, en dash parsing, corrected rc scaling")
    print(PRINT, f"Metadata file: {site_info_csv}")
    print(PRINT, f"Input directory: {input_dir}")
    print(PRINT, f"Output directory: {output_dir}")
    
    process_and_save_all(input_dir, output_dir, site_info_csv)
    print(PRINT, "DONE.")
######################################################################################################################
# PET plots

  ## PET graphs - UPDATED for IGBP classification (FIXED)
# -*- coding: utf-8 -*-
"""
Make exactly three plots from PET outputs in one folder:
  1) PET_mean_annual_by_site.png        (all sites in one plot)
  2) PET_mean_annual_by_IGBP.png        (changed from biome to IGBP)
  3) PET_mean_annual_by_climate.png

Reads all US-*.csv in PET_drought (expects DateTime, PET).
Joins IGBP+climate from site_lat_long.csv (site_name matches CSV stem).
"""

from pathlib import Path
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import textwrap

PRINT = "[PET-PLOTS]"

PET_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\PET_drought"
META_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info"

def _read_pet_long(pet_dir):
    """Read all site CSVs and return long df: site, DateTime, PET."""
    pet_dir = Path(pet_dir)
    files = sorted(glob.glob(str(pet_dir / "US-*.csv")))
    if not files:
        print(PRINT, f"Warning: No US-*.csv files found in {pet_dir}")
        # Try without US- prefix
        files = sorted(glob.glob(str(pet_dir / "*.csv")))
    
    rows = []
    for f in files:
        site = Path(f).stem
        try:
            df = pd.read_csv(f)
            if "DateTime" not in df.columns or "PET" not in df.columns:
                print(PRINT, f"[WARN] Skipping {site}: missing DateTime or PET")
                continue
            dt = pd.to_datetime(df["DateTime"], errors="coerce")
            pet = pd.to_numeric(df["PET"], errors="coerce")
            m = dt.notna() & pet.notna()
            if not m.any():
                print(PRINT, f"[WARN] Skipping {site}: no valid PET timepoints")
                continue
            rows.append(pd.DataFrame({"site": site, "DateTime": dt[m], "PET": pet[m]}))
        except Exception as e:
            print(PRINT, f"[ERROR] {site}: {e}")
    
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=["site","DateTime","PET"])

def _annual_totals(df_long):
    """Annual total PET per site (mm/year)."""
    d = df_long.copy()
    d["year"] = pd.to_datetime(d["DateTime"]).dt.year
    return (d.groupby(["site","year"], as_index=False)["PET"]
              .sum()
              .rename(columns={"PET":"PET_total_mm"}))

def _mean_annual_by_site(annual_df):
    """Mean of annual totals per site across available years."""
    return (annual_df.groupby("site", as_index=False)
            .agg(mean_annual_PET_mm=("PET_total_mm","mean"),
                 years_n=("year","nunique")))

def _load_meta(meta_dir):
    """
    Load site metadata from site_lat_long.csv.
    UPDATED: Uses IGBP instead of biome.
    Expected columns: site_name, IGBP, climate (optional)
    """
    csv_path = Path(meta_dir) / "site_lat_long.csv"
    
    # Try multiple encodings
    encodings = ['utf-8', 'latin1', 'cp1252', 'ISO-8859-1']
    meta = None
    for encoding in encodings:
        try:
            meta = pd.read_csv(csv_path, encoding=encoding)
            print(PRINT, f"Successfully read metadata with {encoding} encoding")
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
        except Exception as e:
            print(PRINT, f"Error with {encoding}: {e}")
            continue
    
    if meta is None:
        raise ValueError(f"Could not read {csv_path} with any encoding")
    
    # Print available columns for debugging
    print(PRINT, f"Metadata columns found: {list(meta.columns)}")
    
    # Check for required columns (flexible naming)
    # Look for site name column
    site_col = None
    for col in ['site_name', 'site', 'Site', 'SITE', 'name']:
        if col in meta.columns:
            site_col = col
            break
    
    if site_col is None:
        raise ValueError(f"No site name column found. Available: {list(meta.columns)}")
    
    # Look for IGBP column (instead of biome)
    igbp_col = None
    for col in ['IGBP', 'igbp', 'IGBP_class', 'Vegetation', 'vegetation', 'biome']:
        if col in meta.columns:
            igbp_col = col
            break
    
    if igbp_col is None:
        print(PRINT, "[WARN] No IGBP column found. Available columns:", list(meta.columns))
        # Fall back to creating a default column
        meta['IGBP'] = 'UNKNOWN'
        igbp_col = 'IGBP'
    
    # Look for climate column (optional)
    climate_col = None
    for col in ['climate', 'Climate', 'CLIMATE', 'water_class']:
        if col in meta.columns:
            climate_col = col
            break
    
    # Prepare metadata
    cols_to_keep = [site_col, igbp_col]
    if climate_col:
        cols_to_keep.append(climate_col)
    
    meta = meta[cols_to_keep].copy()
    meta = meta.rename(columns={site_col: 'site_name', igbp_col: 'IGBP'})
    if climate_col:
        meta = meta.rename(columns={climate_col: 'climate'})
    else:
        meta['climate'] = 'unknown'
    
    # Clean strings
    meta["site_name"] = meta["site_name"].astype(str).str.strip()
    meta["IGBP"] = meta["IGBP"].astype(str).str.strip().str.upper()
    meta["climate"] = meta["climate"].astype(str).str.strip()
    
    # Remove duplicates
    meta = meta.drop_duplicates(subset=['site_name'])
    
    print(PRINT, f"Loaded {len(meta)} sites with IGBP classes: {sorted(meta['IGBP'].unique())}")
    return meta

def _wrap_labels_auto(labels, n_items, min_width=8, max_width=20, max_lines=3):
    """
    Choose a wrap width from number of categories, then wrap each label.
    Heuristic: more items -> narrower width.
    """
    est_width = max(min_width, min(max_width, int(70 / max(1, n_items))))
    wrapped = []
    max_obs = 1
    for lab in map(str, labels):
        # Handle potential NaN
        if pd.isna(lab) or lab == 'nan':
            lab = 'Unknown'
        parts = textwrap.wrap(lab, width=est_width, break_long_words=False, break_on_hyphens=True)
        if not parts:
            parts = [lab]
        parts = parts[:max_lines]
        max_obs = max(max_obs, len(parts))
        wrapped.append("\n".join(parts))
    return wrapped, max_obs, est_width

def _plot_bar(x_labels, values, title, ylabel, out_path, rotate=0, wrap=True, wrap_lines=3, sort_values=True):
    """
    Bar plot with multi-line x tick labels that don't overlap neighboring ticks.
    - Auto wrap width based on number of categories.
    - Horizontal, centered, multi-line labels.
    - Extra bottom margin and x tick padding.
    """
    # Sort by values if requested
    if sort_values:
        sorted_idx = np.argsort(values)[::-1]  # descending
        x_labels = [x_labels[i] for i in sorted_idx]
        values = [values[i] for i in sorted_idx]
    
    n = len(x_labels)
    if wrap:
        x_labels_display, max_lines, used_width = _wrap_labels_auto(
            x_labels, n_items=n, min_width=10, max_width=18, max_lines=wrap_lines
        )
    else:
        x_labels_display, max_lines = list(map(str, x_labels)), 1

    # Wider figure for many bars; extra height for multi-line labels
    fig_w = max(10, min(20, 0.42 * n + 4))  # Cap at 20 inches wide
    fig_h = 6 + 0.5 * (max_lines - 1)
    plt.figure(figsize=(fig_w, fig_h))

    pos = np.arange(n)
    plt.bar(pos, values, color='steelblue', edgecolor='black', alpha=0.7)

    ax = plt.gca()
    ax.set_xticks(pos)
    if rotate and rotate != 0:
        ax.set_xticklabels(x_labels_display, rotation=rotate, ha="right", fontsize=9)
    else:
        ax.set_xticklabels(x_labels_display, rotation=0, ha="center", fontsize=9)

    # Padding so wrapped lines don't touch the axis
    ax.tick_params(axis='x', pad=8)
    # Small side margins so edge labels don't get clipped
    ax.margins(x=0.02)

    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.grid(axis='y', alpha=0.3, linestyle='--')

    # Extra bottom margin for multi-line labels
    bottom_margin = min(0.55, 0.22 + 0.10 * (max_lines - 1))
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom=bottom_margin)

    plt.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close()
    print(PRINT, f"Saved: {out_path}")

def make_three_plots(pet_dir=PET_DIR, meta_dir=META_DIR):
    """Main function to generate three PET plots."""
    pet_dir = Path(pet_dir)
    meta_dir = Path(meta_dir)

    print(PRINT, "Reading PET files…")
    df_long = _read_pet_long(pet_dir)
    if df_long.empty:
        print(PRINT, "No PET data found. Exiting.")
        return

    print(PRINT, f"Read {df_long['site'].nunique()} sites with {len(df_long)} observations")
    
    print(PRINT, "Computing annual totals and site means…")
    annual = _annual_totals(df_long)
    by_site = _mean_annual_by_site(annual)
    by_site = by_site.sort_values("mean_annual_PET_mm", ascending=False)

    # Merge metadata for IGBP/climate aggregation
    try:
        meta = _load_meta(meta_dir)
        by_site_meta = by_site.merge(meta, left_on="site", right_on="site_name", how="left")
        print(PRINT, f"Merged {len(by_site_meta)} sites with metadata")
    except Exception as e:
        print(PRINT, f"[WARN] Could not attach meta: {e}")
        by_site_meta = by_site.copy()
        by_site_meta["IGBP"] = "UNKNOWN"
        by_site_meta["climate"] = "unknown"

    # 1) All sites in one plot
    print(PRINT, "Plot 1: PET_mean_annual_by_site.png")
    _plot_bar(
        x_labels=by_site_meta["site"].tolist(),
        values=by_site_meta["mean_annual_PET_mm"].values,
        title="Mean Annual Total PET by Site",
        ylabel="PET (mm year⁻¹)",
        out_path=pet_dir / "PET_mean_annual_by_site.png",
        rotate=45 if len(by_site_meta) > 20 else 0,  # Rotate if many sites
        wrap=False,
        sort_values=True
    )

    # 2) By IGBP (was biome) - FIXED: Properly handle aggregation
    print(PRINT, "Plot 2: PET_mean_annual_by_IGBP.png")
    
    # Filter out rows with missing IGBP
    igbp_data = by_site_meta.dropna(subset=["IGBP"]).copy()
    
    # Calculate statistics by IGBP class
    igbp_group = igbp_data.groupby("IGBP")["mean_annual_PET_mm"].agg(['mean', 'std', 'count']).reset_index()
    igbp_group.columns = ['IGBP', 'mean_annual_PET_mm', 'std_PET', 'n_sites']
    igbp_group = igbp_group.sort_values("mean_annual_PET_mm", ascending=False)
    
    # Add count to labels for context
    labels = [f"{row['IGBP']}\n(n={int(row['n_sites'])})" 
              for _, row in igbp_group.iterrows()]
    
    _plot_bar(
        x_labels=labels,
        values=igbp_group["mean_annual_PET_mm"].values,
        title="Mean Annual Total PET by IGBP Class",
        ylabel="PET (mm year⁻¹)",
        out_path=pet_dir / "PET_mean_annual_by_IGBP.png",
        rotate=0,
        wrap=True,
        wrap_lines=3,
        sort_values=False  # Already sorted
    )
    
    # Print summary statistics by IGBP
    print(PRINT, "\nPET by IGBP class:")
    for _, row in igbp_group.iterrows():
        print(PRINT, f"  {row['IGBP']:4s}: {row['mean_annual_PET_mm']:.0f} ± {row['std_PET']:.0f} mm/yr (n={int(row['n_sites'])})")

    # 3) By climate (if available) - FIXED: Properly handle aggregation
    print(PRINT, "Plot 3: PET_mean_annual_by_climate.png")
    
    # Check if we have climate data
    if 'climate' in by_site_meta.columns and by_site_meta['climate'].nunique() > 1:
        climate_data = by_site_meta.dropna(subset=["climate"]).copy()
        
        if len(climate_data) > 0 and climate_data['climate'].nunique() > 1:
            # Calculate statistics by climate class
            climate_group = climate_data.groupby("climate")["mean_annual_PET_mm"].agg(['mean', 'std', 'count']).reset_index()
            climate_group.columns = ['climate', 'mean_annual_PET_mm', 'std_PET', 'n_sites']
            climate_group = climate_group.sort_values("mean_annual_PET_mm", ascending=False)
            
            labels = [f"{row['climate']}\n(n={int(row['n_sites'])})" 
                      for _, row in climate_group.iterrows()]
            
            _plot_bar(
                x_labels=labels,
                values=climate_group["mean_annual_PET_mm"].values,
                title="Mean Annual Total PET by Climate Class",
                ylabel="PET (mm year⁻¹)",
                out_path=pet_dir / "PET_mean_annual_by_climate.png",
                rotate=0,
                wrap=True,
                wrap_lines=3,
                sort_values=False
            )
            
            print(PRINT, "\nPET by climate class:")
            for _, row in climate_group.iterrows():
                print(PRINT, f"  {row['climate']}: {row['mean_annual_PET_mm']:.0f} ± {row['std_PET']:.0f} mm/yr (n={int(row['n_sites'])})")
        else:
            print(PRINT, "[WARN] Not enough climate data for plot 3")
            _create_placeholder_plot(pet_dir / "PET_mean_annual_by_climate.png", "Climate")
    else:
        print(PRINT, "[WARN] No climate column found or insufficient data for plot 3")
        _create_placeholder_plot(pet_dir / "PET_mean_annual_by_climate.png", "Climate")

    print(PRINT, f"\nAll plots saved in: {pet_dir}")

def _create_placeholder_plot(out_path, category):
    """Create a placeholder plot when data is missing."""
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, f'Insufficient {category} data\nfor PET by {category} plot', 
            ha='center', va='center', fontsize=14, transform=plt.gca().transAxes)
    plt.title(f"PET Mean Annual by {category} (Data Unavailable)")
    plt.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close()
    print(PRINT, f"Created placeholder plot: {out_path}")

# Run directly
if __name__ == "__main__":
    make_three_plots()      

######################################################################################################
####  drought indices

# SPEI batch processor with robust data handling and multi-scale support
# UPDATED: CORRECT log-logistic distribution with MONTHLY-SPECIFIC fitting
# Fits separate distribution for each calendar month (Jan-Dec) to account for seasonality
# Following the standard SPEI package methodology
####################################################################################################

import os
import glob
import warnings
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, List
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

# ============================ CONFIGURATION ============================
IN_DIR  = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\PET_drought"
OUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\drought"

DATE_COL   = "DateTime"
PRECIP_COL = "precip_mm"
PET_COL    = "PET"

# SPEI timescales (months) - adjust based on your data length requirements
SCALES = (1, 3, 6, 12, 24, 36, 48)

# Standard reference period following WMO recommendations
CALIB_YEARS = (1991, 2020)

# Data quality thresholds
MIN_MONTHLY_COVERAGE   = 0.70   # Keep months with ≥70% coverage
LOW_COVERAGE_THRESHOLD = 0.50   # Scale months with 50-70% coverage

# ============================ CORE FUNCTIONS ============================

def minimum_required_months(scale: int, min_total: int = 60) -> int:
    """
    Calculate minimum months required for reliable SPEI calculation.
    
    Parameters:
    -----------
    scale : int
        SPEI timescale in months
    min_total : int
        Absolute minimum months required regardless of scale
        
    Returns:
    --------
    int : Minimum months required
    """
    return max(min_total, 5 * scale)


class DataQualityChecker:
    """Helper class for data quality assessment and reporting."""
    
    @staticmethod
    def assess_data_sufficiency(valid_months: int, scales: Tuple[int, ...]) -> Dict[str, Any]:
        """Assess which scales can be reliably computed."""
        feasible_scales = []
        marginal_scales = []
        infeasible_scales = []
        
        for scale in scales:
            min_required = minimum_required_months(scale)
            if valid_months >= min_required:
                feasible_scales.append(scale)
            elif valid_months >= scale * 12:  # Absolute minimum
                marginal_scales.append(scale)
            else:
                infeasible_scales.append(scale)
                
        return {
            'feasible': feasible_scales,
            'marginal': marginal_scales,
            'infeasible': infeasible_scales,
            'total_months': valid_months
        }
    
    @staticmethod
    def print_quality_report(assessment: Dict[str, Any]):
        """Print formatted data quality assessment."""
        print("📊 DATA QUALITY ASSESSMENT:")
        print(f"   Total valid months: {assessment['total_months']}")
        
        if assessment['feasible']:
            scales_str = ", ".join(f"SPEI-{s}" for s in assessment['feasible'])
            print(f"   ✅ Reliable scales: {scales_str}")
        
        if assessment['marginal']:
            scales_str = ", ".join(f"SPEI-{s}" for s in assessment['marginal'])
            print(f"   ⚠️  Marginal scales: {scales_str} (interpret with caution)")
        
        if assessment['infeasible']:
            scales_str = ", ".join(f"SPEI-{s}" for s in assessment['infeasible'])
            print(f"   ❌ Infeasible scales: {scales_str}")


def determine_data_frequency(df: pd.DataFrame, date_col: str) -> str:
    """
    Determine data frequency automatically.
    
    Returns:
    --------
    str : 'half_hourly', 'daily', or 'unknown'
    """
    if len(df) < 2:
        return 'unknown'
    
    df_temp = df.copy()
    df_temp[date_col] = pd.to_datetime(df_temp[date_col], errors='coerce')
    df_temp = df_temp.dropna(subset=[date_col]).sort_values(date_col)
    
    time_diffs = df_temp[date_col].diff().dt.total_seconds().dropna()
    
    if len(time_diffs) == 0:
        return 'unknown'
    
    mode_diff = time_diffs.mode().iloc[0] if not time_diffs.mode().empty else time_diffs.iloc[0]
    
    if 1500 <= mode_diff <= 2100:  # ~30-35 minutes
        return 'half_hourly'
    elif 82800 <= mode_diff <= 90000:  # ~23-25 hours
        return 'daily'
    else:
        return f'unknown_{mode_diff:.0f}s'


def aggregate_to_monthly(df: pd.DataFrame, date_col: str, precip_col: str, 
                        pet_col: str, data_frequency: str) -> pd.DataFrame:
    """
    Aggregate sub-monthly data to monthly sums with coverage tracking.
    """
    df_clean = df.copy()
    df_clean[date_col] = pd.to_datetime(df_clean[date_col], errors='coerce')
    df_clean = df_clean.dropna(subset=[date_col]).sort_values(date_col)
    
    # Add counter for coverage calculation
    df_clean["__ones__"] = 1.0
    
    # Resample to monthly
    monthly = (
        df_clean.set_index(date_col)
        .resample("MS")
        .agg({precip_col: "sum", pet_col: "sum", "__ones__": "sum"})
        .rename(columns={precip_col: "P", pet_col: "PET", "__ones__": "n_present"})
    )
    
    # Calculate expected records based on data frequency
    days_in_month = pd.Index(monthly.index).days_in_month
    
    if data_frequency == 'half_hourly':
        monthly["n_expected"] = 48 * days_in_month
    elif data_frequency == 'daily':
        monthly["n_expected"] = days_in_month
    else:
        # Conservative default
        monthly["n_expected"] = days_in_month
        print(f"   ⚠️  Unknown frequency '{data_frequency}', using daily assumption")
    
    monthly["coverage"] = monthly["n_present"] / monthly["n_expected"]
    
    return monthly


def apply_coverage_correction(monthly: pd.DataFrame, low_threshold: float, 
                            high_threshold: float) -> pd.DataFrame:
    """
    Apply coverage-based corrections to monthly data.
    
    Note: Precipitation is clipped at 0 (cannot be negative).
    PET is left unclipped as negative values can occur in cold months
    and are physically meaningful.
    """
    monthly_corrected = monthly.copy()
    coverage = monthly["coverage"].to_numpy()
    
    P_raw = monthly["P"].to_numpy(dtype=float)
    PET_raw = monthly["PET"].to_numpy(dtype=float)
    
    P_adj = np.full_like(P_raw, np.nan)
    PET_adj = np.full_like(P_raw, np.nan)
    
    # Case 1: High coverage - keep raw values
    mask_keep = coverage >= high_threshold
    P_adj[mask_keep] = P_raw[mask_keep]
    PET_adj[mask_keep] = PET_raw[mask_keep]
    
    # Case 2: Medium coverage - scale values
    mask_scale = (coverage >= low_threshold) & (coverage < high_threshold)
    scale_factors = np.clip(1.0 / coverage[mask_scale], 1.0, 1.5)
    P_adj[mask_scale] = P_raw[mask_scale] * scale_factors
    PET_adj[mask_scale] = PET_raw[mask_scale] * scale_factors
    
    # Case 3: Low coverage - already NaN
    
    monthly_corrected["P"] = P_adj
    monthly_corrected["PET"] = PET_adj
    
    # Remove negative precipitation (physically impossible)
    monthly_corrected["P"] = monthly_corrected["P"].clip(lower=0)
    
    # Note: PET is NOT clipped - negative values can occur in cold months
    # and represent conditions where potential evapotranspiration is near zero
    
    return monthly_corrected


def get_calibration_period(monthly_index: pd.DatetimeIndex, 
                          requested_calib: Optional[Tuple[int, int]] = None) -> Tuple[int, int]:
    """
    Determine effective calibration period considering data availability.
    """
    data_start = int(monthly_index.min().year)
    data_end = int(monthly_index.max().year)
    
    if requested_calib is None:
        return data_start, data_end
    
    req_start, req_end = requested_calib
    eff_start = max(req_start, data_start)
    eff_end = min(req_end, data_end)
    
    if eff_start > eff_end:
        print(f"   ⚠️  Calibration {req_start}-{req_end} outside data span, using {data_start}-{data_end}")
        return data_start, data_end
    elif (eff_start, eff_end) != (req_start, req_end):
        print(f"   ⚠️  Calibration clipped to {eff_start}-{eff_end}")
    
    return eff_start, eff_end


# ============================ MONTHLY-SPECIFIC LOG-LOGISTIC SPEI ============================

def loglogistic_cdf(x: float, alpha: float, beta: float, gamma: float) -> float:
    """
    Log-logistic cumulative distribution function.
    
    Following Vicente-Serrano et al. (2010), Eq. 2:
    F(x) = 1 / (1 + (α/(x-γ))^β)
    
    Parameters:
    -----------
    x : float
        Value at which to evaluate CDF
    alpha : float
        Scale parameter (>0)
    beta : float
        Shape parameter (>0)
    gamma : float
        Location parameter (can be negative)
        
    Returns:
    --------
    float : Cumulative probability (0 to 1)
    """
    if x <= gamma:
        return 0.0
    
    try:
        term = (alpha / (x - gamma)) ** beta
        cdf = 1.0 / (1.0 + term)
        return np.clip(cdf, 0.0, 1.0)
    except (OverflowError, FloatingPointError):
        if x > gamma:
            return 1.0
        else:
            return 0.0


def fit_loglogistic_3param_monthly(data_by_month: Dict[int, np.ndarray]) -> Dict[int, Tuple[float, float, float]]:
    """
    Fit 3-parameter log-logistic distribution SEPARATELY for each calendar month.
    
    This follows the standard SPEI methodology where distribution parameters
    are estimated for each month (Jan-Dec) to account for seasonality.
    
    Parameters:
    -----------
    data_by_month : Dict[int, np.ndarray]
        Dictionary mapping month (1-12) to array of accumulated values
        for that month during the calibration period
        
    Returns:
    --------
    Dict[int, Tuple[float, float, float]] : Month -> (alpha, beta, gamma) parameters
    """
    params_by_month = {}
    
    for month in range(1, 13):
        data = data_by_month.get(month, np.array([]))
        data_clean = data[np.isfinite(data)]
        
        if len(data_clean) < 10:  # Need at least 10 years of data for reliable fitting
            print(f"      ⚠️  Month {month}: insufficient data ({len(data_clean)} values), using global fit")
            params_by_month[month] = (np.nan, np.nan, np.nan)
            continue
        
        # Negative log-likelihood function for 3-parameter log-logistic
        def neg_log_likelihood(params):
            alpha, beta, gamma = params
            
            # Parameter constraints
            if alpha <= 0 or beta <= 0:
                return 1e10
            
            # Log-likelihood
            ll = 0.0
            for x in data_clean:
                if x <= gamma:
                    return 1e10  # Invalid if x <= gamma (probability zero)
                
                try:
                    term = ((x - gamma) / alpha) ** beta
                    # Log of PDF: ln(β/α) + (β-1)*ln((x-γ)/α) - 2*ln(1 + ((x-γ)/α)^β)
                    log_pdf = np.log(beta) - np.log(alpha) + (beta - 1) * np.log((x - gamma) / alpha) - 2 * np.log(1 + term)
                    ll += log_pdf
                except (OverflowError, FloatingPointError, ValueError):
                    return 1e10
            
            return -ll
        
        # Initial parameter guesses using method of moments
        data_sorted = np.sort(data_clean)
        gamma_init = np.min(data_sorted) - 0.01
        
        # Estimate alpha and beta from percentiles
        p25 = np.percentile(data_sorted, 25)
        p75 = np.percentile(data_sorted, 75)
        
        if p75 > p25:
            beta_init = np.log(9) / (np.log(p75 - gamma_init) - np.log(p25 - gamma_init))
            alpha_init = (p75 - gamma_init) / (3 ** (1/beta_init))
        else:
            beta_init = 2.0
            alpha_init = np.std(data_clean)
        
        beta_init = max(beta_init, 0.5)
        alpha_init = max(alpha_init, 0.1)
        
        # Optimize
        try:
            result = minimize(
                neg_log_likelihood,
                [alpha_init, beta_init, gamma_init],
                method='L-BFGS-B',
                bounds=[(1e-6, None), (1e-6, None), (None, np.min(data_clean) - 1e-6)]
            )
            
            if result.success:
                alpha, beta, gamma = result.x
                params_by_month[month] = (alpha, beta, gamma)
            else:
                print(f"      ⚠️  Month {month}: optimization failed, using global fit")
                params_by_month[month] = (np.nan, np.nan, np.nan)
                
        except Exception as e:
            print(f"      ⚠️  Month {month}: fitting error ({e}), using global fit")
            params_by_month[month] = (np.nan, np.nan, np.nan)
    
    return params_by_month


def fit_loglogistic_3param_global(data: np.ndarray) -> Tuple[float, float, float]:
    """
    Fit 3-parameter log-logistic distribution to all data pooled (fallback method).
    
    Used when monthly-specific fitting has insufficient data.
    
    Parameters:
    -----------
    data : np.ndarray
        All accumulated values during calibration period
        
    Returns:
    --------
    Tuple[float, float, float] : (alpha, beta, gamma) parameters
    """
    data_clean = data[np.isfinite(data)]
    
    if len(data_clean) < 24:
        return (np.nan, np.nan, np.nan)
    
    def neg_log_likelihood(params):
        alpha, beta, gamma = params
        
        if alpha <= 0 or beta <= 0:
            return 1e10
        
        ll = 0.0
        for x in data_clean:
            if x <= gamma:
                return 1e10
            
            try:
                term = ((x - gamma) / alpha) ** beta
                log_pdf = np.log(beta) - np.log(alpha) + (beta - 1) * np.log((x - gamma) / alpha) - 2 * np.log(1 + term)
                ll += log_pdf
            except (OverflowError, FloatingPointError, ValueError):
                return 1e10
        
        return -ll
    
    data_sorted = np.sort(data_clean)
    gamma_init = np.min(data_sorted) - 0.01
    
    p25 = np.percentile(data_sorted, 25)
    p75 = np.percentile(data_sorted, 75)
    
    if p75 > p25:
        beta_init = np.log(9) / (np.log(p75 - gamma_init) - np.log(p25 - gamma_init))
        alpha_init = (p75 - gamma_init) / (3 ** (1/beta_init))
    else:
        beta_init = 2.0
        alpha_init = np.std(data_clean)
    
    beta_init = max(beta_init, 0.5)
    alpha_init = max(alpha_init, 0.1)
    
    try:
        result = minimize(
            neg_log_likelihood,
            [alpha_init, beta_init, gamma_init],
            method='L-BFGS-B',
            bounds=[(1e-6, None), (1e-6, None), (None, np.min(data_clean) - 1e-6)]
        )
        
        if result.success:
            return result.x
        else:
            return (np.nan, np.nan, np.nan)
            
    except Exception:
        return (np.nan, np.nan, np.nan)


def cdf_to_spei(cdf: float) -> float:
    """
    Convert cumulative probability to standardized SPEI value.
    
    Uses inverse normal (Gaussian) transformation:
    SPEI = Φ^(-1)(F(x))
    where Φ is the standard normal CDF.
    
    Parameters:
    -----------
    cdf : float
        Cumulative probability from log-logistic distribution
        
    Returns:
    --------
    float : Standardized SPEI value (mean ≈ 0, std ≈ 1)
    """
    cdf_clipped = np.clip(cdf, 1e-6, 1 - 1e-6)
    
    try:
        spei = stats.norm.ppf(cdf_clipped)
        return spei
    except Exception:
        return np.nan


def compute_spei_loglogistic_correct(
    water_balance: np.ndarray,
    scale: int,
    calib_start: int,
    calib_end: int,
    dates: pd.DatetimeIndex
) -> np.ndarray:
    """
    Compute SPEI using MONTHLY-SPECIFIC 3-parameter log-logistic distribution.
    
    Following the standard SPEI methodology:
    1. Calculate water balance D = P - PET
    2. Accumulate D over specified timescale (FULL WINDOW ONLY)
    3. Fit separate 3-parameter log-logistic for each calendar month
    4. Convert probabilities to standardized normal values using month-specific parameters
    
    Parameters:
    -----------
    water_balance : np.ndarray
        Monthly water balance (P - PET) values
    scale : int
        SPEI timescale in months
    calib_start : int
        Start year for calibration period
    calib_end : int
        End year for calibration period
    dates : pd.DatetimeIndex
        Monthly dates corresponding to water_balance
        
    Returns:
    --------
    np.ndarray : SPEI values for each month (first scale-1 months are NaN)
    """
    n_months = len(water_balance)
    spei_values = np.full(n_months, np.nan)
    
    # Calculate accumulated water balance over FULL WINDOW only
    accumulated = np.full(n_months, np.nan)
    
    for i in range(scale - 1, n_months):
        start_idx = i - scale + 1
        acc_values = water_balance[start_idx:i+1]
        
        # Only compute if we have exactly 'scale' months of valid data
        if len(acc_values) == scale and np.all(np.isfinite(acc_values)):
            accumulated[i] = np.nansum(acc_values)
    
    # Prepare data for monthly-specific fitting
    # Group accumulated values by calendar month for calibration period
    calib_mask = (dates.year >= calib_start) & (dates.year <= calib_end)
    
    data_by_month = {month: [] for month in range(1, 13)}
    for i in range(n_months):
        if calib_mask[i] and np.isfinite(accumulated[i]):
            month = dates[i].month
            data_by_month[month].append(accumulated[i])
    
    # Convert to numpy arrays
    for month in range(1, 13):
        data_by_month[month] = np.array(data_by_month[month])
    
    # Fit monthly-specific distributions
    print(f"   Fitting monthly distributions...", end=" ")
    params_by_month = fit_loglogistic_3param_monthly(data_by_month)
    
    # Check if any monthly fits succeeded
    months_with_params = sum(1 for m in range(1, 13) if not np.any(np.isnan(params_by_month[m])))
    
    if months_with_params < 6:
        print(f"⚠️  Only {months_with_params}/12 months have valid fits, using global fit")
        # Fall back to global fit
        all_calib_data = np.concatenate([data_by_month[m] for m in range(1, 13)])
        alpha_g, beta_g, gamma_g = fit_loglogistic_3param_global(all_calib_data)
        
        if np.any(np.isnan([alpha_g, beta_g, gamma_g])):
            print(f"      ❌ Global fit also failed")
            return spei_values
        
        # Use global parameters for all months
        for month in range(1, 13):
            params_by_month[month] = (alpha_g, beta_g, gamma_g)
        print(f"   Using global parameters")
    else:
        print(f"✓ {months_with_params}/12 months fitted")
        
        # For months with failed fits, use the average of successful months
        successful_params = [params_by_month[m] for m in range(1, 13) 
                            if not np.any(np.isnan(params_by_month[m]))]
        if successful_params:
            avg_alpha = np.mean([p[0] for p in successful_params])
            avg_beta = np.mean([p[1] for p in successful_params])
            avg_gamma = np.mean([p[2] for p in successful_params])
            
            for month in range(1, 13):
                if np.any(np.isnan(params_by_month[month])):
                    params_by_month[month] = (avg_alpha, avg_beta, avg_gamma)
                    print(f"      Month {month}: using average parameters")
    
    # Calculate SPEI using month-specific parameters
    valid_count = 0
    for i in range(n_months):
        if np.isfinite(accumulated[i]):
            month = dates[i].month
            alpha, beta, gamma = params_by_month[month]
            
            if not np.any(np.isnan([alpha, beta, gamma])):
                prob = loglogistic_cdf(accumulated[i], alpha, beta, gamma)
                spei_values[i] = cdf_to_spei(prob)
                valid_count += 1
    
    # Clip extreme values
    spei_values = np.clip(spei_values, -5.0, 5.0)
    
    print(f"({valid_count} valid months with full {scale}-month windows)", end=" ")
    
    return spei_values


def compute_spei_safely(
    P: np.ndarray, 
    PET: np.ndarray, 
    scale: int,
    calib_start: int, 
    calib_end: int, 
    dates: pd.DatetimeIndex
) -> np.ndarray:
    """
    Compute SPEI with monthly-specific log-logistic implementation.
    
    This implements the standard SPEI method (Vicente-Serrano et al., 2010)
    with proper 3-parameter log-logistic distribution, full window accumulation,
    AND monthly-specific distribution fitting to account for seasonality.
    """
    try:
        # Calculate water balance (P - PET)
        water_balance = P - PET
        
        # Compute SPEI using correct log-logistic distribution
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            spei_values = compute_spei_loglogistic_correct(
                water_balance, scale, calib_start, calib_end, dates
            )
        
        spei_array = np.asarray(spei_values, dtype=float)
        
        # Validate results
        finite_mask = np.isfinite(spei_array)
        if finite_mask.any():
            finite_values = spei_array[finite_mask]
            
            # Check for extreme values
            extreme_count = np.sum((finite_values < -5) | (finite_values > 5))
            if extreme_count > len(finite_values) * 0.3:
                print(f"⚠️  Excessive extremes ({extreme_count}/{len(finite_values)}), likely fitting issue")
                spei_array[:] = np.nan
            elif extreme_count > 0:
                print(f"ℹ️  {extreme_count} extreme values (beyond ±5)")
        
        return spei_array
        
    except Exception as e:
        print(f"❌ Computation failed: {e}")
        import traceback
        traceback.print_exc()
        return np.full_like(P, np.nan, dtype=float)


# ============================ MAIN SPEI PIPELINE ============================

def spei_from_dataframe(
    df: pd.DataFrame,
    date_col: str = "DateTime",
    precip_col: str = "precip_mm",
    pet_col: str = "PET",
    scales: Tuple[int, ...] = (1, 3, 6, 12, 24, 36, 48),
    calib_years: Optional[Tuple[int, int]] = None,
    min_monthly_coverage: float = 0.70,
    low_coverage_threshold: float = 0.50,
) -> pd.DataFrame:
    """
    Calculate SPEI from sub-monthly precipitation and PET data.
    Uses MONTHLY-SPECIFIC 3-parameter LOG-LOGISTIC distribution (standard SPEI method).
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with datetime, precipitation, and PET columns
    date_col : str
        Name of datetime column
    precip_col : str
        Name of precipitation column (mm)
    pet_col : str
        Name of PET column (mm)
    scales : tuple
        SPEI timescales to compute (months)
    calib_years : tuple
        Calibration period (start_year, end_year)
    min_monthly_coverage : float
        Minimum coverage for keeping data (0-1)
    low_coverage_threshold : float
        Threshold for scaling vs dropping data (0-1)
        
    Returns:
    --------
    pd.DataFrame with monthly SPEI values
    """
    # ==================== INPUT VALIDATION ====================
    if date_col not in df.columns:
        raise ValueError(f"Missing datetime column '{date_col}'")
    if precip_col not in df.columns:
        raise ValueError(f"Missing precipitation column '{precip_col}'")
    if pet_col not in df.columns:
        raise ValueError(f"Missing PET column '{pet_col}'")
    
    print(f"📁 Processing dataset: {len(df):,} records")
    print(f"📅 Date range: {df[date_col].min()} to {df[date_col].max()}")
    
    # ==================== DATA FREQUENCY DETECTION ====================
    data_freq = determine_data_frequency(df, date_col)
    print(f"📊 Data frequency: {data_freq}")
    
    # ==================== MONTHLY AGGREGATION ====================
    monthly = aggregate_to_monthly(df, date_col, precip_col, pet_col, data_freq)
    print(f"📈 Aggregated to {len(monthly)} monthly periods")
    
    # ==================== EDGE MONTH HANDLING ====================
    if len(monthly) >= 2:
        edge_indices = [monthly.index[0], monthly.index[-1]]
        monthly.loc[edge_indices, ["P", "PET", "n_present", "coverage"]] = np.nan
        print("   ⏹️  Masked edge months (often incomplete)")
    
    # ==================== COVERAGE CORRECTION ====================
    monthly = apply_coverage_correction(monthly, low_coverage_threshold, min_monthly_coverage)
    
    # Report coverage statistics
    valid_coverage = monthly["coverage"].dropna()
    if not valid_coverage.empty:
        cov_stats = (valid_coverage.min(), valid_coverage.max(), 
                    (monthly["coverage"] >= min_monthly_coverage).sum())
        print(f"📋 Coverage: {cov_stats[0]:.3f}-{cov_stats[1]:.3f}, {cov_stats[2]} months ≥{min_monthly_coverage}")
    
    # ==================== DATA QUALITY CHECK ====================
    valid_mask = monthly[["P", "PET"]].notna().all(axis=1)
    valid_count = valid_mask.sum()
    print(f"✅ Valid months: {valid_count}/{len(monthly)}")
    
    # Check for negative PET values (normal in cold months)
    pet_neg_count = (monthly["PET"] < 0).sum()
    if pet_neg_count > 0:
        print(f"   ℹ️  {pet_neg_count} months with negative PET (normal in cold conditions)")
    
    if valid_count == 0:
        print("❌ No valid data after quality control")
        for scale in scales:
            monthly[f"SPEI_{scale}"] = np.nan
        return monthly.drop(columns=["n_present", "n_expected", "coverage"])
    
    # Extract final arrays
    P = monthly["P"].to_numpy()
    PET = monthly["PET"].to_numpy()
    dates = monthly.index
    
    # ==================== CALIBRATION PERIOD ====================
    if calib_years is None:
        calib_years = CALIB_YEARS
    
    calib_start, calib_end = get_calibration_period(monthly.index, calib_years)
    calib_mask = ((dates.year >= calib_start) & (dates.year <= calib_end))
    calib_valid = (valid_mask & calib_mask).sum()
    
    print(f"🎯 Reference period: {calib_start}-{calib_end} (following WMO recommendations)")
    print(f"   Valid months in reference period: {calib_valid}")
    print(f"📊 Distribution: MONTHLY-SPECIFIC 3-PARAMETER LOG-LOGISTIC")
    print(f"   Following Vicente-Serrano et al. (2010) and standard SPEI package")
    print(f"   ✓ Separate distribution for each calendar month (Jan-Dec)")
    print(f"   ✓ Location parameter NOT fixed to zero (handles negative values)")
    print(f"   ✓ Full window accumulation only (no partial months)")
    
    if calib_valid < 48:
        print(f"   ⚠️  Short reference period - long scales may be unstable")
    
    # ==================== DATA ASSESSMENT ====================
    assessment = DataQualityChecker.assess_data_sufficiency(valid_count, scales)
    DataQualityChecker.print_quality_report(assessment)
    
    # ==================== SPEI COMPUTATION ====================
    print("🔄 Computing SPEI scales (monthly-specific log-logistic)...")
    
    for scale in scales:
        min_required = minimum_required_months(scale)
        
        if valid_count < min_required:
            status = "marginal" if valid_count >= scale * 12 else "infeasible"
            print(f"   ⚠️  SPEI-{scale:2d}: {valid_count:3d} months ({status})")
            monthly[f"SPEI_{scale}"] = np.nan
            continue
        
        print(f"   🔄 SPEI-{scale:2d}: computing (monthly-specific)...", end=" ")
        spei_array = compute_spei_safely(P, PET, scale, calib_start, calib_end, dates)
        valid_spei = np.isfinite(spei_array).sum()
        monthly[f"SPEI_{scale}"] = spei_array
        print(f"✓ {valid_spei:3d} valid")
        
        # Report first valid date for this scale
        first_valid_idx = np.where(np.isfinite(spei_array))[0]
        if len(first_valid_idx) > 0:
            first_date = dates[first_valid_idx[0]]
            print(f"      (First valid: {first_date.strftime('%Y-%m')})")
    
    return monthly.drop(columns=["n_present", "n_expected", "coverage"])


# ============================ BATCH PROCESSING ============================

def read_csv_robust(path: str) -> pd.DataFrame:
    """Robust CSV reading with fallback options."""
    try:
        return pd.read_csv(path, low_memory=False)
    except Exception as e:
        print(f"   ⚠️  Standard reader failed: {e}, trying Python engine")
        return pd.read_csv(path, engine="python", on_bad_lines="skip")


def process_one_file(in_csv: str, out_dir: str):
    """Process a single CSV file through the SPEI pipeline."""
    name = Path(in_csv).name
    
    try:
        df = read_csv_robust(in_csv)
        print(f"📖 {name}: {len(df):,} rows")
    except Exception as e:
        print(f"❌ READ-ERROR {name}: {e}")
        return
    
    try:
        monthly = spei_from_dataframe(
            df,
            date_col=DATE_COL,
            precip_col=PRECIP_COL,
            pet_col=PET_COL,
            scales=SCALES,
            calib_years=CALIB_YEARS,
            min_monthly_coverage=MIN_MONTHLY_COVERAGE,
            low_coverage_threshold=LOW_COVERAGE_THRESHOLD,
        )
    except Exception as e:
        print(f"❌ SPEI-ERROR {name}: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Merge back to original temporal resolution
    try:
        out = df.copy()
        out[DATE_COL] = pd.to_datetime(out[DATE_COL], errors="coerce")
        out = out.dropna(subset=[DATE_COL]).sort_values(DATE_COL)
        
        # Create month period for merging
        out["__month__"] = out[DATE_COL].dt.to_period("M")
        monthly["__month__"] = monthly.index.to_period("M")
        
        # Merge SPEI columns
        spei_cols = [f"SPEI_{s}" for s in SCALES]
        out = out.merge(monthly[["__month__"] + spei_cols], on="__month__", how="left")
        out = out.drop(columns="__month__")
        
        print(f"🔗 Merged {len(spei_cols)} SPEI columns (monthly-specific log-logistic)")
        
    except Exception as e:
        print(f"❌ MERGE-ERROR {name}: {e}")
        return
    
    # Write output
    try:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_path = Path(out_dir) / name
        out.to_csv(out_path, index=False)
        print(f"💾 Saved: {out_path}")
    except Exception as e:
        print(f"❌ WRITE-ERROR {name}: {e}")


def main():
    """Main batch processing function."""
    files = sorted(glob.glob(os.path.join(IN_DIR, "*.csv")))
    
    if not files:
        print(f"❌ No CSV files found in {IN_DIR}")
        return
    
    print(f"🚀 Starting SPEI batch processing")
    print(f"📂 Input: {IN_DIR}")
    print(f"💾 Output: {OUT_DIR}")
    print(f"📊 Timescales: {SCALES} months")
    print(f"🎯 Reference period: {CALIB_YEARS[0]}-{CALIB_YEARS[1]} (WMO recommendation)")
    print(f"📊 Distribution: MONTHLY-SPECIFIC 3-PARAMETER LOG-LOGISTIC")
    print(f"   Following Vicente-Serrano et al. (2010) and standard SPEI package")
    print(f"   ✓ Separate distribution for each calendar month")
    print(f"   ✓ Location parameter NOT fixed to zero")
    print(f"   ✓ Handles negative accumulated water balance")
    print(f"   ✓ Full window accumulation only")
    print(f"📋 Coverage: keep ≥{MIN_MONTHLY_COVERAGE}, scale {LOW_COVERAGE_THRESHOLD}-{MIN_MONTHLY_COVERAGE}")
    print(f"💧 PET handling: Negative values allowed (normal in cold months)")
    print(f"📈 Files to process: {len(files)}\n")
    
    successful = 0
    for i, file_path in enumerate(files, 1):
        print(f"\n{'='*60}")
        print(f"📦 [{i}/{len(files)}] {Path(file_path).name}")
        print(f"{'='*60}")
        
        try:
            process_one_file(file_path, OUT_DIR)
            successful += 1
        except Exception as e:
            print(f"💥 Unexpected error: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'='*60}")
    print(f"🎉 BATCH COMPLETE: {successful}/{len(files)} files processed successfully")
    print(f"🎉 Method: MONTHLY-SPECIFIC 3-parameter log-logistic SPEI")
    print(f"🎉 Reference period: {CALIB_YEARS[0]}-{CALIB_YEARS[1]}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
#####################################################################################################