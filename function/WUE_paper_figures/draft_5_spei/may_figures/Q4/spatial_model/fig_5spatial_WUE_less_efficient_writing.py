# -*- coding: utf-8 -*-
"""
DIAGNOSTIC REGIONAL SUMMARY FOR FIG. 1 (WUE_T only)
Reads agg_FULLPERIOD_WUE_tra.nc and summarizes "frac_less" by coastal regions.
Uses the SAME coastal regions as Figure 2 for consistency.

Interpretation:
  frac_less = fraction of monthly time steps (Mar 2000–Sep 2025; 182 months)
              where the pixel is classified as "Less efficient" (Class=1).

Higher frac_less  -> more frequent less-efficient response
Lower frac_less   -> more frequent non-decrease response (Increase or NoChange)

Outputs:
  1) Printed summary table
  2) CSV saved next to the NetCDF
  3) Detailed breakdown by region
"""

import os
import numpy as np
import pandas as pd
import xarray as xr

# =============================================================================
# CONFIG
# =============================================================================
# Using WUE_T input file
INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_time_aggregated\agg_FULLPERIOD_WUE_tra.nc"

# If your ds lon is 0–360, convert to -180–180
SHIFT_LON_0_360 = True

# Descriptive thresholds for reporting (no ecological classification)
# These are used only to report percentages; no categorical labels are assigned.
REPORT_THRESH_HIGH = 0.50    # "frac_less ≥ 0.50"
REPORT_THRESH_VHIGH = 0.66   # "frac_less ≥ 0.66"
REPORT_THRESH_LOW = 0.20     # "frac_less ≤ 0.20"

# Map extent (must match your Figure 1)
CONUS_EXTENT = [-126, -66, 24, 50]
ALASKA_EXTENT = [-165, -130, 52, 72]

# =============================================================================
# COASTAL REGION DEFINITIONS (Matching Figure 2 - 5 regions)
# =============================================================================
COASTAL_REGIONS = {
    "Pacific": {
        "lon_bounds": [-130, -116],
        "lat_bounds": [32, 49],
        "states": ["California", "Oregon", "Washington"]
    },
    "Gulf": {
        "lon_bounds": [-98, -80],
        "lat_bounds": [24, 31],
        "states": ["Texas", "Louisiana", "Mississippi", "Alabama", "Florida"]
    },
    "Southeast Atlantic": {
        "lon_bounds": [-82, -75],
        "lat_bounds": [30, 37],
        "states": ["Georgia", "South Carolina", "North Carolina"]
    },
    "Atlantic North": {
        "lon_bounds": [-78, -66],
        "lat_bounds": [37, 46],
        "states": ["Virginia", "Maryland", "Delaware", "New Jersey", 
                   "New York", "Connecticut", "Rhode Island", 
                   "Massachusetts", "New Hampshire", "Maine"]
    },
    "Alaska": {
        "lon_bounds": ALASKA_EXTENT[:2],
        "lat_bounds": ALASKA_EXTENT[2:],
        "states": ["Alaska"]
    },
    "CONUS Total": {
        "lon_bounds": CONUS_EXTENT[:2],
        "lat_bounds": CONUS_EXTENT[2:],
        "states": ["All CONUS states"]
    }
}

# =============================================================================
# HELPERS
# =============================================================================
def shift_lon_0_360_to_minus180_180(ds, lon_name="lon"):
    lon = ds[lon_name].values
    if np.nanmax(lon) > 180:
        new_lon = ((lon + 180) % 360) - 180
        ds = ds.assign_coords({lon_name: new_lon}).sortby(lon_name)
    return ds

def get_lon_lat_grids(ds):
    """Return lon2d, lat2d arrays regardless of 1D/2D coords."""
    lon = ds["lon"].values
    lat = ds["lat"].values
    if lon.ndim == 1 and lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lon2d, lat2d = lon, lat
    return lon2d, lat2d

def area_weights_from_lat(lat2d):
    """Approx area weights proportional to cos(lat)."""
    w = np.cos(np.deg2rad(lat2d))
    w = np.where(np.isfinite(w), w, np.nan)
    return w

def create_region_mask(lon2d, lat2d, region_def):
    """Create mask for a region based on longitude/latitude bounds."""
    mask = (
        (lon2d >= region_def["lon_bounds"][0]) & 
        (lon2d <= region_def["lon_bounds"][1]) &
        (lat2d >= region_def["lat_bounds"][0]) & 
        (lat2d <= region_def["lat_bounds"][1])
    )
    return mask

def summarize_region(name, region_def, frac_less, count_valid, w_area):
    """Compute unweighted + area-weighted summaries for a region."""
    lon2d, lat2d = region_def["_lon2d"], region_def["_lat2d"]
    mask = create_region_mask(lon2d, lat2d, region_def)
    
    # Valid pixels: must have finite frac_less and positive count_valid
    valid = mask & np.isfinite(frac_less) & np.isfinite(count_valid) & (count_valid > 0)

    n = int(np.sum(valid))
    if n == 0:
        return {
            "Region": name,
            "States": ", ".join(region_def["states"]),
            "n_pixels": 0,
            "mean_frac_less": np.nan,
            "median_frac_less": np.nan,
            "iqr_frac_less": np.nan,
            f"frac_less_ge_{REPORT_THRESH_HIGH:.2f}": np.nan,
            f"frac_less_ge_{REPORT_THRESH_VHIGH:.2f}": np.nan,
            f"frac_less_le_{REPORT_THRESH_LOW:.2f}": np.nan,
            "area_wt_mean": np.nan,
        }

    vals = frac_less[valid]
    mean_ = float(np.nanmean(vals))
    med_ = float(np.nanmedian(vals))
    q25, q75 = np.nanpercentile(vals, [25, 75])
    iqr_ = float(q75 - q25)

    pct_high = float(100.0 * np.mean(vals >= REPORT_THRESH_HIGH))
    pct_vhigh = float(100.0 * np.mean(vals >= REPORT_THRESH_VHIGH))
    pct_low = float(100.0 * np.mean(vals <= REPORT_THRESH_LOW))

    # area-weighted mean
    ww = w_area[valid]
    if np.all(~np.isfinite(ww)) or np.nansum(ww) == 0:
        aw_mean = np.nan
    else:
        aw_mean = float(np.nansum(vals * ww) / np.nansum(ww))

    return {
        "Region": name,
        "States": ", ".join(region_def["states"]),
        "n_pixels": n,
        "mean_frac_less": mean_,
        "median_frac_less": med_,
        "iqr_frac_less": iqr_,
        f"frac_less_ge_{REPORT_THRESH_HIGH:.2f}": pct_high,
        f"frac_less_ge_{REPORT_THRESH_VHIGH:.2f}": pct_vhigh,
        f"frac_less_le_{REPORT_THRESH_LOW:.2f}": pct_low,
        "area_wt_mean": aw_mean,
    }

def print_detailed_summary(df):
    """Print detailed interpretation of results using neutral language."""
    print("\n" + "="*90)
    print("DETAILED REGIONAL SUMMARY (WUE$_T$)")
    print("="*90)
    
    # Find regions with highest and lowest mean frac_less (excluding CONUS Total)
    df_regions = df[~df['Region'].isin(['CONUS Total', 'Alaska'])].copy()
    highest_region = df_regions.loc[df_regions['mean_frac_less'].idxmax()]
    lowest_region = df_regions.loc[df_regions['mean_frac_less'].idxmin()]
    
    print(f"\n1. REGION WITH HIGHEST mean frac_less: {highest_region['Region']}")
    print(f"   • Mean frac_less: {highest_region['mean_frac_less']:.3f}")
    print(f"   • Median: {highest_region['median_frac_less']:.3f}")
    print(f"   • {highest_region[f'frac_less_ge_{REPORT_THRESH_HIGH:.2f}']:.1f}% of pixels have frac_less ≥ {REPORT_THRESH_HIGH}")
    print(f"   • States: {highest_region['States']}")
    
    print(f"\n2. REGION WITH LOWEST mean frac_less: {lowest_region['Region']}")
    print(f"   • Mean frac_less: {lowest_region['mean_frac_less']:.3f}")
    print(f"   • Median: {lowest_region['median_frac_less']:.3f}")
    print(f"   • {lowest_region[f'frac_less_le_{REPORT_THRESH_LOW:.2f}']:.1f}% of pixels have frac_less ≤ {REPORT_THRESH_LOW}")
    print(f"   • States: {lowest_region['States']}")
    
    # Regional comparisons
    print("\n3. REGIONAL COMPARISONS:")
    
    pacific_value = df[df['Region'] == 'Pacific']['mean_frac_less'].values[0]
    atlantic_north_value = df[df['Region'] == 'Atlantic North']['mean_frac_less'].values[0]
    
    print(f"   • Pacific Coast mean frac_less: {pacific_value:.3f}")
    print(f"   • Atlantic North Coast mean frac_less: {atlantic_north_value:.3f}")
    if pacific_value > atlantic_north_value:
        print(f"   → Pacific Coast shows {((pacific_value/atlantic_north_value)-1)*100:.1f}% higher mean frac_less")
    else:
        print(f"   → Atlantic North Coast shows {((atlantic_north_value/pacific_value)-1)*100:.1f}% higher mean frac_less")
    
    gulf_value = df[df['Region'] == 'Gulf']['mean_frac_less'].values[0]
    southeast_value = df[df['Region'] == 'Southeast Atlantic']['mean_frac_less'].values[0]
    print(f"   • Gulf of America mean frac_less: {gulf_value:.3f}")
    print(f"   • Southeast Atlantic mean frac_less: {southeast_value:.3f}")
    
    alaska_value = df[df['Region'] == 'Alaska']['mean_frac_less'].values[0]
    conus_value = df[df['Region'] == 'CONUS Total']['mean_frac_less'].values[0]
    print(f"   • Alaska mean frac_less: {alaska_value:.3f}")
    print(f"   • CONUS average mean frac_less: {conus_value:.3f}")
    
    print("\n4. INTERPRETATION NOTES:")
    print("   • frac_less = proportion of monthly time steps with 'Less efficient' response (Class=1)")
    print("   • Higher values = more frequent less-efficient response")
    print("   • Lower values = more frequent non-decrease response (Increase or NoChange)")
    print("   • Values are based on modeled classification frequency from point-based logistic regression")
    print("   • IQR shows within-region variability of frac_less")

# =============================================================================
# MAIN
# =============================================================================
def main():
    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(f"Input file not found:\n{INPUT_FILE}")

    ds = xr.open_dataset(INPUT_FILE)

    if SHIFT_LON_0_360:
        ds = shift_lon_0_360_to_minus180_180(ds, lon_name="lon")

    # Required variable for Fig 1
    if "frac_less" not in ds:
        raise ValueError("Expected variable 'frac_less' not found in dataset.")

    frac_less = ds["frac_less"].values
    count_valid = ds["count_valid"].values if "count_valid" in ds else np.where(np.isfinite(frac_less), 1, np.nan)

    lon2d, lat2d = get_lon_lat_grids(ds)
    w_area = area_weights_from_lat(lat2d)
    
    # Add grid data to region definitions
    for region_name, region_def in COASTAL_REGIONS.items():
        region_def["_lon2d"] = lon2d
        region_def["_lat2d"] = lat2d

    # -------------------------------------------------------------------------
    # SUMMARIZE EACH REGION
    # -------------------------------------------------------------------------
    rows = []
    for name, region_def in COASTAL_REGIONS.items():
        rows.append(summarize_region(name, region_def, frac_less, count_valid, w_area))

    df = pd.DataFrame(rows)

    # Sort by mean frac_less descending
    df = df.sort_values("mean_frac_less", ascending=False)

    # -------------------------------------------------------------------------
    # SAVE + PRINT
    # -------------------------------------------------------------------------
    out_csv = os.path.splitext(INPUT_FILE)[0] + "_COASTAL_REGIONS_DIAGNOSTICS.csv"
    df.to_csv(out_csv, index=False)

    pd.set_option("display.max_columns", 200)
    pd.set_option("display.width", 200)

    print("\n" + "="*90)
    print("COASTAL REGION DIAGNOSTICS (WUE$_T$ - Matching Figure 2)")
    print("="*90)
    print(f"Input: {INPUT_FILE}")
    print(f"Saved: {out_csv}")
    print("\nREGIONS ANALYZED (5 coastal regions from Figure 2):")
    for region_name, region_def in COASTAL_REGIONS.items():
        if region_name not in ["CONUS Total"]:
            print(f"  • {region_name}: {', '.join(region_def['states'])}")
    
    print("\nVARIABLE DEFINITION:")
    print("  frac_less = proportion of monthly time steps classified as 'Less efficient' (Class=1)")
    print("  Higher values = more frequent less-efficient response")
    print("  Lower values = more frequent non-decrease response (Increase or NoChange)")
    
    print(f"\nREPORTING THRESHOLDS (descriptive only, no ecological classification):")
    print(f"  • Pixels with frac_less ≥ {REPORT_THRESH_HIGH:.2f}")
    print(f"  • Pixels with frac_less ≥ {REPORT_THRESH_VHIGH:.2f}")
    print(f"  • Pixels with frac_less ≤ {REPORT_THRESH_LOW:.2f}")
    
    print("\n" + "="*90)
    print("SUMMARY TABLE (sorted by mean frac_less, highest first)")
    print("="*90 + "\n")
    print(df.to_string(index=False, float_format=lambda x: f"{x:0.3f}" if isinstance(x, float) else str(x)))
    
    # Print detailed summary
    print_detailed_summary(df)
    
    # Additional statistics
    print("\n" + "="*90)
    print("OVERALL STATISTICS (WUE$_T$)")
    print("="*90)
    
    conus_data = df[df['Region'] == 'CONUS Total'].iloc[0]
    print(f"CONUS Total (all regions combined):")
    print(f"  • Mean frac_less: {conus_data['mean_frac_less']:.3f}")
    print(f"  • Pixels with frac_less ≥ {REPORT_THRESH_HIGH:.2f}: {conus_data[f'frac_less_ge_{REPORT_THRESH_HIGH:.2f}']:.1f}%")
    print(f"  • Pixels with frac_less ≤ {REPORT_THRESH_LOW:.2f}: {conus_data[f'frac_less_le_{REPORT_THRESH_LOW:.2f}']:.1f}%")
    
    # Ranking by mean frac_less
    print(f"\nRANKING BY MEAN FRAC_LESS (1 = highest mean):")
    df_rank = df[~df['Region'].isin(['CONUS Total', 'Alaska'])].copy()
    df_rank['Rank'] = df_rank['mean_frac_less'].rank(ascending=False).astype(int)
    for _, row in df_rank.iterrows():
        print(f"  {row['Rank']}. {row['Region']}: {row['mean_frac_less']:.3f}")

    ds.close()
    return df, out_csv

if __name__ == "__main__":
    df, out_csv = main()