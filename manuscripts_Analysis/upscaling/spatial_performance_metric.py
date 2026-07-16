#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diagnose_pixel_level_spatial_resistance_analogue.py

Console-only diagnostic for pixel-level spatial drought-resistance analogue.
Uses the exact proximity-based coast classifier from the Q4 workflow.
Reads dry_mean_pct_change raster, computes resistance = 1 + pct/100,
assigns coasts, and prints detailed per-coast summaries including
block bootstrap CIs.

ADDED: Saves pixel-level CSV, coast summary CSV, and boxplot figures.
"""

import sys
import math
import warnings
from pathlib import Path
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
INPUT_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
OUTPUT_DIR = INPUT_DIR / "spatial_performance_metric"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DRY_MEAN_FILE = INPUT_DIR / "EDI_response_v2_dry_mean_pct_change_upland.nc"
# Optional files
OPTIONAL_FILES = {
    "dry_freq": INPUT_DIR / "EDI_response_v2_dry_month_frequency_upland.nc",
    "inc_freq": INPUT_DIR / "EDI_response_v2_dry_freq_strong_increase_upland.nc",
    "dec_freq": INPUT_DIR / "EDI_response_v2_dry_freq_strong_decrease_upland.nc",
    "net_adj": INPUT_DIR / "EDI_response_v2_drought_net_adjustment_minus_impairment_upland.nc",
}

BLOCK_DEGREES = 3.0
N_BOOTSTRAP = 10000
RANDOM_SEED = 42
PRACTICAL_THRESHOLD = 0.005
EARTH_RADIUS_KM = 6371.0088

# ----------------------------------------------------------------------
# Exact proximity-based coast classifier (from 24-EDI_spatial.py)
# ----------------------------------------------------------------------
REGION_NAMES = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]

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

def _segment_distance_grid(LAT, LON, a_lat, a_lon, b_lat, b_lon):
    """
    Vectorized minimum distance from every grid cell to coastline segment AB.
    Uses local equirectangular projection.
    """
    lat0 = np.radians(LAT)

    ax = EARTH_RADIUS_KM * np.radians(a_lon - LON) * np.cos(lat0)
    ay = EARTH_RADIUS_KM * np.radians(a_lat - LAT)

    bx = EARTH_RADIUS_KM * np.radians(b_lon - LON) * np.cos(lat0)
    by = EARTH_RADIUS_KM * np.radians(b_lat - LAT)

    abx = bx - ax
    aby = by - ay
    denom = abx * abx + aby * aby

    with np.errstate(invalid="ignore", divide="ignore"):
        t = -(ax * abx + ay * aby) / denom

    t = np.where(denom == 0, 0.0, t)
    t = np.clip(t, 0.0, 1.0)

    closest_x = ax + t * abx
    closest_y = ay + t * aby

    return np.sqrt(closest_x * closest_x + closest_y * closest_y)

def _distance_to_polyline_grid(LAT, LON, polyline):
    """Minimum distance from every grid cell to one coastline polyline."""
    min_dist = np.full(LAT.shape, np.inf, dtype=np.float64)
    for i in range(len(polyline) - 1):
        d = _segment_distance_grid(LAT, LON, *polyline[i], *polyline[i+1])
        min_dist = np.minimum(min_dist, d)
    return min_dist

def assign_coast_regions_to_grid(LAT, LON):
    """
    Fast vectorized coast assignment using the exact proximity-based rule:
      lat >= 50 -> AK Coast
      otherwise -> nearest Pacific/Gulf/Atlantic polyline
    """
    coast_id = np.full(LAT.shape, -1, dtype=np.int8)
    region_to_id = {name: i for i, name in enumerate(REGION_NAMES)}

    valid = np.isfinite(LAT) & np.isfinite(LON)

    # Alaska first
    ak = valid & (LAT >= 50)
    coast_id[ak] = region_to_id["AK Coast"]

    # Remaining pixels
    candidate = valid & (~ak)

    best_dist = np.full(LAT.shape, np.inf, dtype=np.float64)

    for name, polyline in COAST_POLYLINES.items():
        d = _distance_to_polyline_grid(LAT, LON, polyline)
        update = candidate & (d < best_dist)
        coast_id[update] = region_to_id[name]
        best_dist[update] = d[update]

    return coast_id

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def read_netcdf_var(file_path, var_name=None):
    """
    Read a NetCDF file and return data, lat, lon, and the variable name used.
    Handles masked arrays and _FillValue safely.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    ds = nc.Dataset(file_path, 'r')

    # Identify variable: if var_name given, use it; else take first non-coordinate var
    if var_name is None:
        all_vars = list(ds.variables.keys())
        coord_vars = set(ds.dimensions.keys())
        data_vars = [v for v in all_vars if v not in coord_vars and ds.variables[v].ndim >= 2]
        if not data_vars:
            raise ValueError(f"No 2D+ variable found in {file_path}")
        var_name = data_vars[0]

    # Read data, fill value, lat, lon before closing
    data = ds.variables[var_name][:]
    fill_val = getattr(ds.variables[var_name], '_FillValue', None)

    # Get lat/lon
    lat = ds.variables.get('lat', ds.variables.get('latitude', None))
    lon = ds.variables.get('lon', ds.variables.get('longitude', None))
    if lat is None or lon is None:
        raise ValueError("Could not find lat/lon variables.")
    lat = lat[:]
    lon = lon[:]

    ds.close()

    # Convert masked array to ndarray with NaN for masked values
    if isinstance(data, np.ma.MaskedArray):
        data = np.ma.filled(data, np.nan).astype(float)
    else:
        data = data.astype(float)

    # Apply fill value
    if fill_val is not None:
        data = np.where(data == fill_val, np.nan, data)
    # Treat huge values as missing
    data = np.where(np.abs(data) > 1e30, np.nan, data)

    return data, lat, lon, var_name

# ----------------------------------------------------------------------
# Main diagnostic
# ----------------------------------------------------------------------
def main():
    print("\n" + "=" * 80)
    print("PIXEL-LEVEL SPATIAL DROUGHT-RESISTANCE ANALOGUE DIAGNOSTIC")
    print("=" * 80)

    # 1. Read dry mean percent change
    print("\nReading dry_mean_pct_change raster...")
    if not DRY_MEAN_FILE.exists():
        raise FileNotFoundError(f"Required input file missing: {DRY_MEAN_FILE}")

    pct_data, lat1d, lon1d, var_name = read_netcdf_var(DRY_MEAN_FILE)
    print(f"  Variable used: {var_name}")
    print(f"  Shape: {pct_data.shape} (lat x lon)")

    # Ensure lat/lon are 1D
    if lat1d.ndim == 2:
        lat1d = lat1d[:, 0]
    if lon1d.ndim == 2:
        lon1d = lon1d[0, :]
    nlat, nlon = pct_data.shape

    # Build 2D lat/lon grids
    lat_grid, lon_grid = np.meshgrid(lat1d, lon1d, indexing='ij')
    # Area weight: cos(lat in radians)
    lat_rad = np.deg2rad(lat_grid)
    area_weight = np.cos(lat_rad)

    # 2. Compute resistance
    print("Computing resistance = 1 + pct_change/100 ...")
    resistance = 1.0 + pct_data / 100.0

    # Read optional rasters (we'll store them in dict for later use)
    opt_data = {}
    domain_mask = np.isfinite(pct_data)  # fallback

    # Dry frequency
    dry_freq_data = None
    if OPTIONAL_FILES["dry_freq"].exists():
        print("\nFound dry_month_frequency raster; using it to define the upland domain.")
        dry_freq_data, _, _, _ = read_netcdf_var(OPTIONAL_FILES["dry_freq"])
        if dry_freq_data.shape != resistance.shape:
            print("  WARNING: dry_freq raster has different shape; cannot use for domain mask.")
            dry_freq_data = None
        else:
            # domain_mask: any pixel that is finite in dry_freq (i.e., part of the upland domain)
            domain_mask = np.isfinite(dry_freq_data)
            opt_data["dry_freq"] = dry_freq_data
    else:
        print("\nWARNING: dry_month_frequency raster not found. Cannot reliably count total domain pixels.")
        print("         'Pixels excluded (no dry data)' will be based on resistance validity only.")
        domain_mask = np.isfinite(resistance)

    # Read other optional rasters if available and matching shape
    for key, fpath in OPTIONAL_FILES.items():
        if key == "dry_freq":
            continue  # already handled
        if fpath.exists():
            data, _, _, _ = read_netcdf_var(fpath)
            if data.shape == resistance.shape:
                opt_data[key] = data
                print(f"  Loaded {key}")
            else:
                print(f"  WARNING: {key} raster shape mismatch; skipping.")

    # Prepare valid mask: finite resistance and (if dry_freq present) dry_freq > 0
    valid_mask = np.isfinite(resistance)
    if dry_freq_data is not None:
        dry_freq_clean = np.where(np.isfinite(dry_freq_data), dry_freq_data, 0)
        dry_present = dry_freq_clean > 0
        valid_mask = valid_mask & dry_present
    else:
        dry_present = None

    # 3. Assign coast regions using exact proximity-based classifier
    print("\nAssigning coast regions using exact proximity-based classifier...")
    coast_ids = assign_coast_regions_to_grid(lat_grid, lon_grid)
    unassigned = np.sum((coast_ids == -1) & domain_mask)
    if unassigned > 0:
        print(f"  Warning: {unassigned} domain pixels could not be assigned a coast; will be excluded.")

    # ------------------------------------------------------------
    # COLLECT PIXEL-LEVEL DATA FOR CSV AND PLOTTING
    # ------------------------------------------------------------
    # We'll iterate over valid pixels to build a list of dicts
    pixel_data = []
    coast_summary = []  # will be filled later

    # We'll compute block IDs for each valid pixel
    lat_block_flat = np.floor(lat_grid / BLOCK_DEGREES).astype(int)
    lon_block_flat = np.floor(lon_grid / BLOCK_DEGREES).astype(int)

    # We'll use masks to get indices
    valid_indices = np.where(valid_mask & (coast_ids >= 0))
    # Flatten arrays for faster iteration? We'll loop over valid indices.
    for i, j in zip(valid_indices[0], valid_indices[1]):
        coast_id = coast_ids[i, j]
        coast_name = REGION_NAMES[coast_id]
        lat_val = lat_grid[i, j]
        lon_val = lon_grid[i, j]
        res_val = resistance[i, j]
        pct_val = pct_data[i, j]
        aw_val = area_weight[i, j]
        l_block = lat_block_flat[i, j]
        lon_block = lon_block_flat[i, j]

        row = {
            "coast_region": coast_name,
            "lat": lat_val,
            "lon": lon_val,
            "dry_mean_WUE_T_pct_change": pct_val,
            "resistance_analogue": res_val,
            "area_weight": aw_val,
            "lat_block": l_block,
            "lon_block": lon_block,
        }
        # Add optional fields if available
        if "dry_freq" in opt_data:
            row["dry_month_frequency"] = opt_data["dry_freq"][i, j]
        if "inc_freq" in opt_data:
            row["strong_increase_frequency"] = opt_data["inc_freq"][i, j]
        if "dec_freq" in opt_data:
            row["strong_decrease_frequency"] = opt_data["dec_freq"][i, j]
        if "net_adj" in opt_data:
            row["net_adjustment_minus_impairment"] = opt_data["net_adj"][i, j]

        pixel_data.append(row)

    # Convert to DataFrame
    df_pixels = pd.DataFrame(pixel_data)
    # Ensure coast order
    df_pixels["coast_region"] = pd.Categorical(
        df_pixels["coast_region"], categories=REGION_NAMES, ordered=True
    )
    df_pixels = df_pixels.sort_values(["coast_region", "lat", "lon"])

    # Save pixel-level CSV
    pixel_csv_path = OUTPUT_DIR / "pixel_level_spatial_resistance_analogue_by_coast.csv"
    df_pixels.to_csv(pixel_csv_path, index=False)
    print(f"\nPixel-level CSV saved: {pixel_csv_path}")

    # ------------------------------------------------------------
    # PER-COAST SUMMARIES (also used for console output)
    # ------------------------------------------------------------
    # We'll compute per-coast summaries using the same logic as before,
    # but now we can also use the df_pixels grouped by coast.

    # We'll re-use the existing loop structure but compute summary stats
    # and fill the coast_summary list.
    # We'll also use the df_pixels for each coast to compute stats.

    print("\n" + "-" * 80)
    print("PER-COAST SUMMARIES")
    print("-" * 80)

    # We'll store summary rows in a list of dicts
    summary_rows = []

    for coast_id, coast_name in enumerate(REGION_NAMES):
        mask_coast = (coast_ids == coast_id)
        n_total = np.sum(mask_coast & domain_mask)
        # Valid pixels: finite resistance, dry_present if applicable
        mask_valid = valid_mask & mask_coast
        n_valid = np.sum(mask_valid)

        if n_valid == 0:
            print(f"\n{coast_name}: No valid pixels.")
            if dry_freq_data is not None:
                print(f"  Total domain pixels: {n_total}")
            # Still add a summary row with NaNs for completeness
            summary_rows.append({
                "coast_region": coast_name,
                "n_total_domain_pixels": n_total,
                "n_valid_pixels": 0,
                "n_excluded_no_dry_data": n_total,
                "pct_pixels_represented": 0.0,
                # all other fields NaN
            })
            continue

        # Use df_pixels for this coast
        df_coast = df_pixels[df_pixels["coast_region"] == coast_name]
        res_coast = df_coast["resistance_analogue"].values
        weights_coast = df_coast["area_weight"].values
        pct_coast = df_coast["dry_mean_WUE_T_pct_change"].values

        # Basic stats
        mean_unweighted = np.mean(res_coast)
        median = np.median(res_coast)
        std = np.std(res_coast, ddof=1)
        min_val = np.min(res_coast)
        max_val = np.max(res_coast)
        p01 = np.percentile(res_coast, 1)
        p05 = np.percentile(res_coast, 5)
        p25 = np.percentile(res_coast, 25)
        p75 = np.percentile(res_coast, 75)
        p95 = np.percentile(res_coast, 95)
        p99 = np.percentile(res_coast, 99)

        # Area-weighted mean
        mean_weighted = np.average(res_coast, weights=weights_coast)

        # Percent categories
        pct_greater = np.mean(res_coast > 1) * 100
        pct_less = np.mean(res_coast < 1) * 100
        pct_near = np.mean((res_coast >= 0.995) & (res_coast <= 1.005)) * 100
        pct_greater_meaningful = np.mean(res_coast > 1.005) * 100
        pct_less_meaningful = np.mean(res_coast < 0.995) * 100

        aw_pct_greater = np.average((res_coast > 1).astype(float), weights=weights_coast) * 100
        aw_pct_less = np.average((res_coast < 1).astype(float), weights=weights_coast) * 100
        aw_pct_near = np.average(((res_coast >= 0.995) & (res_coast <= 1.005)).astype(float), weights=weights_coast) * 100
        aw_pct_greater_meaningful = np.average((res_coast > 1.005).astype(float), weights=weights_coast) * 100
        aw_pct_less_meaningful = np.average((res_coast < 0.995).astype(float), weights=weights_coast) * 100

        mean_pct = np.mean(pct_coast)
        median_pct = np.median(pct_coast)

        # Block bootstrap for area-weighted mean
        # Use block IDs from df_pixels
        lat_block = df_coast["lat_block"].values
        lon_block = df_coast["lon_block"].values
        block_ids = lat_block * 10000 + lon_block
        unique_blocks, block_inverse = np.unique(block_ids, return_inverse=True)
        n_blocks = len(unique_blocks)

        block_means = np.zeros(n_blocks)
        block_weights = np.zeros(n_blocks)
        for b in range(n_blocks):
            idx = (block_inverse == b)
            res_b = res_coast[idx]
            w_b = weights_coast[idx]
            block_weights[b] = np.sum(w_b)
            if block_weights[b] > 0:
                block_means[b] = np.average(res_b, weights=w_b)
            else:
                block_means[b] = np.nan
        valid_blocks = ~np.isnan(block_means) & (block_weights > 0)
        block_means = block_means[valid_blocks]
        block_weights = block_weights[valid_blocks]
        n_blocks_eff = len(block_means)

        if n_blocks_eff > 1:
            rng = np.random.default_rng(seed=RANDOM_SEED)
            boot_means = np.zeros(N_BOOTSTRAP)
            for i in range(N_BOOTSTRAP):
                sampled_indices = rng.choice(n_blocks_eff, size=n_blocks_eff, replace=True)
                sampled_means = block_means[sampled_indices]
                sampled_weights = block_weights[sampled_indices]
                boot_means[i] = np.average(sampled_means, weights=sampled_weights)
            ci_low = np.percentile(boot_means, 2.5)
            ci_high = np.percentile(boot_means, 97.5)
        else:
            ci_low = ci_high = np.nan
            print(f"  Warning: only {n_blocks_eff} block(s); bootstrap CI not reliable.")

        # Classification
        if np.isnan(ci_low) or np.isnan(ci_high):
            classification = "cannot determine"
        else:
            effect = mean_weighted - 1.0
            if ci_low > 1.0 and effect >= PRACTICAL_THRESHOLD:
                classification = "meaningfully greater than 1"
            elif ci_high < 1.0 and -effect >= PRACTICAL_THRESHOLD:
                classification = "meaningfully less than 1"
            else:
                if (ci_low > 1.0 or ci_high < 1.0) and abs(effect) < PRACTICAL_THRESHOLD:
                    classification = "statistically different but ecologically negligible"
                else:
                    classification = "near-neutral / ecologically negligible"

        # Store summary row
        summary_row = {
            "coast_region": coast_name,
            "n_total_domain_pixels": n_total,
            "n_valid_pixels": n_valid,
            "n_excluded_no_dry_data": n_total - n_valid,
            "pct_pixels_represented": n_valid / n_total * 100 if n_total > 0 else 0.0,
            "mean_dry_mean_WUE_T_pct_change": mean_pct,
            "median_dry_mean_WUE_T_pct_change": median_pct,
            "mean_resistance_unweighted": mean_unweighted,
            "mean_resistance_area_weighted": mean_weighted,
            "median_resistance": median,
            "std_resistance": std,
            "min_resistance": min_val,
            "max_resistance": max_val,
            "p01_resistance": p01,
            "p05_resistance": p05,
            "p25_resistance": p25,
            "p75_resistance": p75,
            "p95_resistance": p95,
            "p99_resistance": p99,
            "pct_pixels_gt_1": pct_greater,
            "pct_pixels_lt_1": pct_less,
            "pct_pixels_near_neutral": pct_near,
            "pct_pixels_meaningfully_gt_1": pct_greater_meaningful,
            "pct_pixels_meaningfully_lt_1": pct_less_meaningful,
            "pct_pixels_gt_1_area_weighted": aw_pct_greater,
            "pct_pixels_lt_1_area_weighted": aw_pct_less,
            "pct_pixels_near_neutral_area_weighted": aw_pct_near,
            "pct_pixels_meaningfully_gt_1_area_weighted": aw_pct_greater_meaningful,
            "pct_pixels_meaningfully_lt_1_area_weighted": aw_pct_less_meaningful,
            "n_spatial_blocks": n_blocks_eff,
            "bootstrap_ci_lower": ci_low,
            "bootstrap_ci_upper": ci_high,
            "classification": classification,
        }
        summary_rows.append(summary_row)

        # Console output (as before)
        print(f"\n{coast_name}:")
        print(f"  Total domain pixels assigned: {n_total}")
        print(f"  Pixels with finite resistance and dry data: {n_valid}")
        print(f"  Pixels excluded (no dry data): {n_total - n_valid}")
        if n_total > 0:
            print(f"  Represented: {n_valid/n_total*100:.2f}% of coast")
        else:
            print("  Represented: N/A (no domain pixels)")
        print(f"  Dry-month WUE_T % change: mean = {mean_pct:.6f}, median = {median_pct:.6f}")
        print(f"  Resistance (unweighted): mean = {mean_unweighted:.6f}, median = {median:.6f}, std = {std:.6f}")
        print(f"  Resistance (area-weighted): mean = {mean_weighted:.6f}")
        print(f"  Min = {min_val:.6f}, Max = {max_val:.6f}")
        print(f"  Percentiles: p01={p01:.6f}, p05={p05:.6f}, p25={p25:.6f}, "
              f"p75={p75:.6f}, p95={p95:.6f}, p99={p99:.6f}")
        print(f"  % pixels > 1: {pct_greater:.2f}%  (area-w: {aw_pct_greater:.2f}%)")
        print(f"  % pixels < 1: {pct_less:.2f}%  (area-w: {aw_pct_less:.2f}%)")
        print(f"  % pixels near-neutral [0.995,1.005]: {pct_near:.2f}%  (area-w: {aw_pct_near:.2f}%)")
        print(f"  % meaningfully >1 (res>1.005): {pct_greater_meaningful:.2f}%  (area-w: {aw_pct_greater_meaningful:.2f}%)")
        print(f"  % meaningfully <1 (res<0.995): {pct_less_meaningful:.2f}%  (area-w: {aw_pct_less_meaningful:.2f}%)")
        print(f"  Spatial blocks (3°): {n_blocks_eff}")
        if not np.isnan(ci_low):
            print(f"  95% block bootstrap CI for area-weighted mean: [{ci_low:.6f}, {ci_high:.6f}]")
        else:
            print("  95% block bootstrap CI: not available (too few blocks)")
        print(f"  Classification: {classification}")

    # Save summary CSV
    df_summary = pd.DataFrame(summary_rows)
    # Ensure coast order
    df_summary["coast_region"] = pd.Categorical(
        df_summary["coast_region"], categories=REGION_NAMES, ordered=True
    )
    df_summary = df_summary.sort_values("coast_region")
    summary_csv_path = OUTPUT_DIR / "spatial_resistance_analogue_summary_by_coast.csv"
    df_summary.to_csv(summary_csv_path, index=False)
    print(f"\nCoast summary CSV saved: {summary_csv_path}")

    # ------------------------------------------------------------
    # PLOTTING: Boxplot
    # ------------------------------------------------------------
    print("\nGenerating boxplot...")
    # Use df_pixels for plotting
    # Prepare coast order
    coast_order = REGION_NAMES
    # Colors
    coast_colors = {
        "AK Coast": "brown",
        "Pacific Coast": "crimson",
        "Gulf Coast": "teal",
        "Atlantic Coast": "green"
    }

    fig, ax = plt.subplots(figsize=(9, 6))
    # Boxplot: we can create a list of data for each coast
    box_data = []
    coast_labels = []
    for coast in coast_order:
        vals = df_pixels[df_pixels["coast_region"] == coast]["resistance_analogue"].dropna().values
        box_data.append(vals)
        coast_labels.append(coast)

    # Create boxplot with no fliers
    bp = ax.boxplot(box_data, labels=coast_labels, showfliers=False,
                    patch_artist=True, boxprops=dict(linewidth=1.5),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    medianprops=dict(linewidth=2, color='black'),
                    widths=0.6)

    # Color boxes
    for patch, coast in zip(bp['boxes'], coast_order):
        patch.set_facecolor(coast_colors.get(coast, 'lightgray'))
        patch.set_alpha(0.7)

    # Add horizontal lines
    ax.axhline(y=1.0, color='black', linestyle='-', linewidth=1.5, alpha=0.8)
    ax.axhline(y=0.995, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(y=1.005, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

    # Add area-weighted mean as black points
    for i, coast in enumerate(coast_order):
        # Get area-weighted mean from summary
        row = df_summary[df_summary["coast_region"] == coast]
        if not row.empty:
            aw_mean = row["mean_resistance_area_weighted"].values[0]
            if not np.isnan(aw_mean):
                ax.scatter(i+1, aw_mean, color='black', s=50, zorder=5, label='Area-weighted mean' if i==0 else "")

    # Add unweighted mean as small colored points
    for i, coast in enumerate(coast_order):
        row = df_summary[df_summary["coast_region"] == coast]
        if not row.empty:
            uw_mean = row["mean_resistance_unweighted"].values[0]
            if not np.isnan(uw_mean):
                ax.scatter(i+1, uw_mean, color='dodgerblue', s=30, zorder=5, label='Unweighted mean' if i==0 else "")

    # Axis labels and title
    ax.set_ylabel("Pixel-level spatial drought-resistance analogue", fontsize=12)
    ax.set_xlabel("Coast region", fontsize=12)
    ax.set_title("Pixel-level spatial drought-resistance analogue by coast", fontsize=14, fontweight='bold')
    ax.tick_params(axis='both', labelsize=11)
    ax.grid(axis='y', linestyle='--', alpha=0.3)

    # Move legend outside if needed
    # We'll just show a simple legend if we have labels
    # Since we have only two types of points, we can add a legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc='upper right', fontsize=10, frameon=True)

    fig.tight_layout()

    # Save PNG and PDF
    png_path = OUTPUT_DIR / "spatial_resistance_analogue_boxplot_by_coast.png"
    pdf_path = OUTPUT_DIR / "spatial_resistance_analogue_boxplot_by_coast.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    print(f"Boxplot PNG saved: {png_path}")
    print(f"Boxplot PDF saved: {pdf_path}")

    # ------------------------------------------------------------
    # OPTIONAL: Boxplot with jittered points (sample)
    # ------------------------------------------------------------
    print("\nGenerating boxplot with jittered points (sample)...")
    # Sample up to 5000 points per coast to avoid clutter
    df_sample = df_pixels.groupby('coast_region', group_keys=False).apply(
        lambda x: x.sample(min(len(x), 5000), random_state=42)
    )
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    bp2 = ax2.boxplot(box_data, labels=coast_labels, showfliers=False,
                      patch_artist=True, boxprops=dict(linewidth=1.5),
                      whiskerprops=dict(linewidth=1.5),
                      capprops=dict(linewidth=1.5),
                      medianprops=dict(linewidth=2, color='black'),
                      widths=0.6)
    for patch, coast in zip(bp2['boxes'], coast_order):
        patch.set_facecolor(coast_colors.get(coast, 'lightgray'))
        patch.set_alpha(0.5)

    # Add jittered points
    for i, coast in enumerate(coast_order):
        vals = df_sample[df_sample["coast_region"] == coast]["resistance_analogue"].values
        if len(vals) > 0:
            # Add jitter
            x_jitter = np.random.normal(i+1, 0.04, size=len(vals))
            ax2.scatter(x_jitter, vals, s=5, alpha=0.2, color=coast_colors.get(coast, 'gray'), rasterized=True)

    # Horizontal lines
    ax2.axhline(y=1.0, color='black', linestyle='-', linewidth=1.5, alpha=0.8)
    ax2.axhline(y=0.995, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.axhline(y=1.005, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

    ax2.set_ylabel("Pixel-level spatial drought-resistance analogue", fontsize=12)
    ax2.set_xlabel("Coast region", fontsize=12)
    ax2.set_title("Pixel-level spatial drought-resistance analogue by coast\n(with jittered points)", fontsize=14, fontweight='bold')
    ax2.tick_params(axis='both', labelsize=11)
    ax2.grid(axis='y', linestyle='--', alpha=0.3)

    fig2.tight_layout()
    png2_path = OUTPUT_DIR / "spatial_resistance_analogue_boxplot_with_points_by_coast.png"
    fig2.savefig(png2_path, dpi=300, bbox_inches='tight')
    print(f"Boxplot with points PNG saved: {png2_path}")

    # ------------------------------------------------------------
    # Additional console messages
    # ------------------------------------------------------------
    print("\n" + "-" * 80)
    print("OUTPUT FILES SUMMARY")
    print("-" * 80)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Pixel-level CSV: {pixel_csv_path}")
    print(f"Coast summary CSV: {summary_csv_path}")
    print(f"Boxplot PNG: {png_path}")
    print(f"Boxplot PDF: {pdf_path}")
    print(f"Boxplot with points PNG: {png2_path}")

    # Continue with rest of console output (ranking, optional rasters, interpretation)
    # ------------------------------------------------------------
    # Ranking by area-weighted mean (from summary)
    # ------------------------------------------------------------
    if not df_summary.empty:
        print("\n" + "-" * 80)
        print("RANKING OF COASTS BY AREA-WEIGHTED MEAN RESISTANCE")
        sorted_df = df_summary.sort_values("mean_resistance_area_weighted", ascending=False)
        for rank, (idx, row) in enumerate(sorted_df.iterrows(), 1):
            print(f"  {rank}. {row['coast_region']}: {row['mean_resistance_area_weighted']:.6f}")

    # ------------------------------------------------------------
    # Supplementary optional rasters summaries (as before)
    # ------------------------------------------------------------
    print("\n" + "-" * 80)
    print("SUPPLEMENTARY SUMMARIES FROM OPTIONAL RASTERS")
    print("-" * 80)

    if "dry_freq" in opt_data:
        print("\nDry-month frequency (per coast, only pixels with dry data):")
        for coast in REGION_NAMES:
            df_coast = df_pixels[df_pixels["coast_region"] == coast]
            if df_coast.empty:
                continue
            freq_vals = df_coast["dry_month_frequency"].values
            weights = df_coast["area_weight"].values
            mean_freq = np.mean(freq_vals)
            median_freq = np.median(freq_vals)
            p05_freq = np.percentile(freq_vals, 5)
            p95_freq = np.percentile(freq_vals, 95)
            aw_mean_freq = np.average(freq_vals, weights=weights)
            print(f"  {coast}: mean={mean_freq:.4f}, median={median_freq:.4f}, "
                  f"p05={p05_freq:.4f}, p95={p95_freq:.4f}, area-w mean={aw_mean_freq:.4f}")

    if "inc_freq" in opt_data:
        print("\nStrong increase frequency (% dry months with pct_change >= +5%):")
        for coast in REGION_NAMES:
            df_coast = df_pixels[df_pixels["coast_region"] == coast]
            if df_coast.empty:
                continue
            vals = df_coast["strong_increase_frequency"].values
            weights = df_coast["area_weight"].values
            mean_inc = np.mean(vals)
            aw_mean_inc = np.average(vals, weights=weights)
            print(f"  {coast}: mean={mean_inc:.4f}, area-w mean={aw_mean_inc:.4f}")

    if "dec_freq" in opt_data:
        print("\nStrong decrease frequency (% dry months with pct_change <= -5%):")
        for coast in REGION_NAMES:
            df_coast = df_pixels[df_pixels["coast_region"] == coast]
            if df_coast.empty:
                continue
            vals = df_coast["strong_decrease_frequency"].values
            weights = df_coast["area_weight"].values
            mean_dec = np.mean(vals)
            aw_mean_dec = np.average(vals, weights=weights)
            print(f"  {coast}: mean={mean_dec:.4f}, area-w mean={aw_mean_dec:.4f}")

    if "net_adj" in opt_data:
        print("\nNet adjustment minus impairment (adj_freq - impair_freq):")
        for coast in REGION_NAMES:
            df_coast = df_pixels[df_pixels["coast_region"] == coast]
            if df_coast.empty:
                continue
            vals = df_coast["net_adjustment_minus_impairment"].values
            weights = df_coast["area_weight"].values
            mean_net = np.mean(vals)
            aw_mean_net = np.average(vals, weights=weights)
            print(f"  {coast}: mean={mean_net:.4f}, area-w mean={aw_mean_net:.4f}")

    # ------------------------------------------------------------
    # Interpretation paragraph (same as before)
    # ------------------------------------------------------------
    print("\n" + "=" * 80)
    print("INTERPRETATION NOTE")
    print("=" * 80)
    print("This pixel-level analysis is the closest spatial analogue to Q2 resistance because")
    print("Q2 calculates one resistance value per site after pooling drought and near-normal months,")
    print("whereas this analysis calculates one model-derived resistance analogue per pixel after")
    print("pooling dry-month exposure across 2000–2025. However, it is not true observed Q2")
    print("resistance because spatial upscaling does not contain observed absolute")
    print("WUE_T_drought / WUE_T_near-normal values. It should be interpreted as a model-derived")
    print("spatial drought-resistance analogue.")
    print("\nHierarchy of resistance analogues:")
    print("  Main spatial resistance analogue:")
    print("    pixel-level resistance analogue from dry_mean_pct_change raster.")
    print("  Supporting diagnostic:")
    print("    annual coast-level resistance analogue from monthly coast summaries.")
    print("\nIMPORTANT: The bootstrap CIs are descriptive uncertainty intervals only—they are")
    print("not formal significance tests and do not account for spatial autocorrelation.")
    print("Pixels are not independent; use these values for exploratory spatial pattern")
    print("identification, not for hypothesis testing.")
    if dry_freq_data is None:
        print("\nNOTE: dry_month_frequency raster was not found; total pixel counts are based on")
        print("      pixels with finite resistance only, which may underestimate total domain.")
    print("=" * 80)
    print("\nAll requested output files have been saved.")

# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()