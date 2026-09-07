#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Q3_ERA5_reviewer_chunk1_processing.py

Processes the merged and blended data, computes all metrics, slopes, fill fractions,
and saves intermediate results to a cache folder for fast figure regeneration.

Run this once before running Chunk 2.
"""

import os
import math
import numpy as np
import pandas as pd
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS AND OUTPUT FOLDER
# ============================================================================
MONTHLY_DATA = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
MERGED_DIR = r"M:\Research\WUE_CUE\ameri_data\reddy_gaps\merged_ameri_era5"
BLENDED_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps\blended_gaps"
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods\reviewer_ERA5_correction"
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# ============================================================================
# COAST CLASSIFIER (EXACT COPY FROM Study_area_august.py)
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
# HELPER: compute metrics
# ============================================================================
def compute_metrics(obs, pred):
    """Return R², bias (pred - obs), RMSE, n."""
    mask = np.isfinite(obs) & np.isfinite(pred)
    if mask.sum() < 2:
        return np.nan, np.nan, np.nan, 0
    obs = obs[mask]; pred = pred[mask]
    slope, intercept, r_value, _, _ = linregress(obs, pred)
    r2 = r_value**2
    bias = np.mean(pred - obs)
    rmse = np.sqrt(np.mean((pred - obs)**2))
    return r2, bias, rmse, len(obs)

# ============================================================================
# STEP 1: RECONSTRUCT FINAL ANALYTICAL SAMPLE
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

n_sites = final_domain['site_name'].nunique()
n_site_years = final_domain.groupby(['site_name', 'Year']).ngroups
n_site_months = len(final_domain)

print("\n=== FINAL ANALYTICAL SAMPLE ===")
print(f"Number of sites: {n_sites}")
print(f"Number of site-years: {n_site_years}")
print(f"Number of site-months: {n_site_months}")
print("\nSite counts by coast (unique sites):")
coast_counts = (
    final_domain[['site_name','coast_region']]
    .drop_duplicates()['coast_region']
    .value_counts()
)
for coast in ['AK Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
    print(f"  {coast}: {coast_counts.get(coast, 0)}")
print("\nYear range and site-years by coast:")
for coast in ['AK Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
    sub = final_domain[final_domain['coast_region'] == coast]
    if not sub.empty:
        ymin = sub['Year'].min()
        ymax = sub['Year'].max()
        nyears = sub.groupby('site_name')['Year'].nunique().sum()
        print(f"  {coast}: {ymin}–{ymax} (site-years: {nyears})")

final_sites = list(final_domain['site_name'].unique())
print(f"\nFinal sites ({len(final_sites)}): {sorted(final_sites)}")

# ============================================================================
# STEP 2: LOAD EXISTING MERGED AND BLENDED DATA, RESTRICT TO FINAL DOMAIN
# ============================================================================
print("\nLoading merged and blended files for each site...")

# Containers for results
results_tair = []
results_rg = []
results_vpd = []
# We'll also collect pooled scatter arrays (will be saved as .npy)
pooled = {
    'tower_tair': [], 'raw_tair': [], 'corr_tair': [],
    'tower_rg': [], 'raw_rg': [], 'corr_rg': [],
    'tower_vpd': [], 'raw_vpd': [], 'corr_vpd': []
}

vpd_correction_sites = ['US-Atq','US-Brw','US-EvM','US-Snd','US-S02','US-S03',
                        'US-S04','US-StS','US-MRf']

for site in final_sites:
    merged_file = os.path.join(MERGED_DIR, f"merged_{site}.csv")
    blended_file = os.path.join(BLENDED_DIR, f"gaps_blend_{site}.csv")
    if not os.path.exists(merged_file) or not os.path.exists(blended_file):
        print(f"  WARNING: Missing merged or blended file for {site}. Skipping.")
        continue

    # Read data
    merged = pd.read_csv(merged_file)
    blended = pd.read_csv(blended_file)
    merged['DateTime'] = pd.to_datetime(merged['DateTime'])
    blended['DateTime'] = pd.to_datetime(blended['DateTime'])

    # Preserve original tower observations
    for col in ['Tair', 'Rg', 'VPD']:
        if col not in merged.columns:
            merged[col] = np.nan

    if site in vpd_correction_sites:
        merged['VPD'] = merged['VPD'] * 10

    # Ensure raw ERA5 columns exist
    if 'Tair_era' not in merged.columns:
        if 't2m' in merged.columns:
            merged['Tair_era'] = merged['t2m'] - 273.15
        else:
            merged['Tair_era'] = np.nan
    if 'Rg_era' not in merged.columns:
        if 'ssrd' in merged.columns:
            merged['Rg_era'] = merged['ssrd'] / 3600.0
        else:
            merged['Rg_era'] = np.nan
    if 'VPD_era' not in merged.columns:
        if 'd2m' in merged.columns and 't2m' in merged.columns:
            es = 0.6108 * np.exp((17.27 * (merged['t2m']-273.15)) / ((merged['t2m']-273.15) + 237.3))
            ea = 0.6108 * np.exp((17.27 * (merged['d2m']-273.15)) / ((merged['d2m']-273.15) + 237.3))
            merged['VPD_era'] = (es - ea) * 10.0
        else:
            merged['VPD_era'] = np.nan

    for var in ['Tair', 'Rg', 'VPD']:
        col_c = f'{var}_era_c'
        if col_c not in blended.columns:
            blended[col_c] = merged[f'{var}_era'].copy()

    # Restrict to final domain
    merged['year_month'] = merged['DateTime'].dt.to_period('M')
    blended['year_month'] = blended['DateTime'].dt.to_period('M')
    domain_years_months = final_domain[final_domain['site_name'] == site][['Year','month']]
    domain_years_months['year_month'] = pd.to_datetime(
        domain_years_months['Year'].astype(str) + '-' + domain_years_months['month'].astype(str) + '-01'
    ).dt.to_period('M')
    domain_set = set(domain_years_months['year_month'])

    merged_dom = merged[merged['year_month'].isin(domain_set)].copy()
    blended_dom = blended[blended['year_month'].isin(domain_set)].copy()
    merged_dom.drop('year_month', axis=1, inplace=True)
    blended_dom.drop('year_month', axis=1, inplace=True)

    if merged_dom.empty or blended_dom.empty:
        print(f"  {site}: No data in final domain.")
        continue

    merged_dom = merged_dom.set_index('DateTime').sort_index()
    blended_dom = blended_dom.set_index('DateTime').sort_index()
    common_idx = merged_dom.index.intersection(blended_dom.index)
    merged_dom = merged_dom.loc[common_idx]
    blended_dom = blended_dom.loc[common_idx]

    # Combine
    combined = pd.DataFrame(index=common_idx)
    combined['Tair_tower'] = merged_dom['Tair']
    combined['Rg_tower'] = merged_dom['Rg']
    combined['VPD_tower'] = merged_dom['VPD']
    combined['Tair_era'] = merged_dom['Tair_era']
    combined['Rg_era'] = merged_dom['Rg_era']
    combined['VPD_era'] = merged_dom['VPD_era']
    combined['Tair_era_c'] = blended_dom['Tair_era_c']
    combined['Rg_era_c'] = blended_dom['Rg_era_c']
    combined['VPD_era_c'] = blended_dom['VPD_era_c']
    combined['Tair_final'] = blended_dom['Tair']
    if 'RH' in merged_dom.columns:
        combined['RH'] = merged_dom['RH']
    else:
        combined['RH'] = np.nan
    combined['Year'] = combined.index.year
    combined['month'] = combined.index.month

    # Collect pooled scatter data
    mask_tair = combined['Tair_tower'].notna() & combined['Tair_era'].notna() & combined['Tair_era_c'].notna()
    if mask_tair.sum() > 0:
        pooled['tower_tair'].extend(combined.loc[mask_tair, 'Tair_tower'].values)
        pooled['raw_tair'].extend(combined.loc[mask_tair, 'Tair_era'].values)
        pooled['corr_tair'].extend(combined.loc[mask_tair, 'Tair_era_c'].values)

    mask_rg = combined['Rg_tower'].notna() & combined['Rg_era'].notna() & combined['Rg_era_c'].notna()
    if mask_rg.sum() > 0:
        pooled['tower_rg'].extend(combined.loc[mask_rg, 'Rg_tower'].values)
        pooled['raw_rg'].extend(combined.loc[mask_rg, 'Rg_era'].values)
        pooled['corr_rg'].extend(combined.loc[mask_rg, 'Rg_era_c'].values)

    mask_vpd = combined['VPD_tower'].notna() & combined['VPD_era'].notna() & combined['VPD_era_c'].notna()
    if mask_vpd.sum() > 0:
        pooled['tower_vpd'].extend(combined.loc[mask_vpd, 'VPD_tower'].values)
        pooled['raw_vpd'].extend(combined.loc[mask_vpd, 'VPD_era'].values)
        pooled['corr_vpd'].extend(combined.loc[mask_vpd, 'VPD_era_c'].values)

    # Compute slopes from raw vs corrected
    # Tair
    raw = combined['Tair_era'].dropna()
    corr = combined['Tair_era_c'].dropna()
    if len(raw) >= 2 and len(corr) >= 2:
        idx = raw.index.intersection(corr.index)
        if len(idx) >= 2:
            slope, intercept, r_value, _, _ = linregress(raw.loc[idx], corr.loc[idx])
            tair_slope, tair_intercept = slope, intercept
        else:
            tair_slope, tair_intercept = np.nan, np.nan
    else:
        tair_slope, tair_intercept = np.nan, np.nan

    # Rg
    raw = combined['Rg_era'].dropna()
    corr = combined['Rg_era_c'].dropna()
    if len(raw) >= 2 and len(corr) >= 2:
        idx = raw.index.intersection(corr.index)
        if len(idx) >= 2:
            slope, intercept, r_value, _, _ = linregress(raw.loc[idx], corr.loc[idx])
            rg_slope, rg_intercept = slope, intercept
        else:
            rg_slope, rg_intercept = np.nan, np.nan
    else:
        rg_slope, rg_intercept = np.nan, np.nan

    # VPD slopes per site-year (store in dict for now, but will be added to results_vpd later)
    vpd_slopes = {}
    vpd_intercepts = {}
    for yr in combined['Year'].unique():
        sub = combined[combined['Year'] == yr]
        raw = sub['VPD_era'].dropna()
        corr = sub['VPD_era_c'].dropna()
        if len(raw) >= 2 and len(corr) >= 2:
            idx = raw.index.intersection(corr.index)
            if len(idx) >= 2:
                slope, intercept, r_value, _, _ = linregress(raw.loc[idx], corr.loc[idx])
                vpd_slopes[(site, yr)] = slope
                vpd_intercepts[(site, yr)] = intercept
            else:
                vpd_slopes[(site, yr)] = np.nan
                vpd_intercepts[(site, yr)] = np.nan
        else:
            vpd_slopes[(site, yr)] = np.nan
            vpd_intercepts[(site, yr)] = np.nan

    coast = final_domain[final_domain['site_name'] == site]['coast_region'].iloc[0]

    # ---- Tair metrics ----
    obs = combined['Tair_tower']
    raw_era = combined['Tair_era']
    corr_era = combined['Tair_era_c']
    r2_before, bias_before, rmse_before, n_before = compute_metrics(obs, raw_era)
    r2_after, bias_after, rmse_after, n_after = compute_metrics(obs, corr_era)
    results_tair.append({
        'site': site, 'coast': coast,
        'slope': tair_slope, 'intercept': tair_intercept,
        'r2_before': r2_before, 'bias_before': bias_before, 'rmse_before': rmse_before,
        'r2_after': r2_after, 'bias_after': bias_after, 'rmse_after': rmse_after,
        'n': n_before
    })

    # ---- Rg metrics ----
    obs = combined['Rg_tower']
    raw_era = combined['Rg_era']
    corr_era = combined['Rg_era_c']
    r2_before, bias_before, rmse_before, n_before = compute_metrics(obs, raw_era)
    r2_after, bias_after, rmse_after, n_after = compute_metrics(obs, corr_era)
    results_rg.append({
        'site': site, 'coast': coast,
        'slope': rg_slope, 'intercept': rg_intercept,
        'r2_before': r2_before, 'bias_before': bias_before, 'rmse_before': rmse_before,
        'r2_after': r2_after, 'bias_after': bias_after, 'rmse_after': rmse_after,
        'n': n_before
    })

    # ---- VPD metrics (site-year) ----
    for yr in combined['Year'].unique():
        sub = combined[combined['Year'] == yr]
        obs = sub['VPD_tower']
        raw_era = sub['VPD_era']
        corr_era = sub['VPD_era_c']
        r2_before, bias_before, rmse_before, n_before = compute_metrics(obs, raw_era)
        r2_after, bias_after, rmse_after, n_after = compute_metrics(obs, corr_era)
        results_vpd.append({
            'site': site, 'year': yr, 'coast': coast,
            'slope': vpd_slopes.get((site, yr), np.nan),
            'intercept': vpd_intercepts.get((site, yr), np.nan),
            'r2_before': r2_before, 'bias_before': bias_before, 'rmse_before': rmse_before,
            'r2_after': r2_after, 'bias_after': bias_after, 'rmse_after': rmse_after,
            'n': n_before
        })

# Convert to DataFrames
df_tair = pd.DataFrame(results_tair)
df_rg = pd.DataFrame(results_rg)
df_vpd = pd.DataFrame(results_vpd)
df_tair['variable'] = 'Tair'
df_rg['variable'] = 'Rg'
df_vpd['variable'] = 'VPD'

# ============================================================================
# STEP 3: COMPUTE ERA5 FILL FRACTIONS
# ============================================================================
print("\nComputing ERA5 fill fractions...")

fill_data = []
# We need to re‑read the per‑site combined data? We didn't store it. We need to recompute fill fractions from the original merged/blended? 
# Actually we can compute fill fractions from the per‑site combined data we already built. But we didn't store it; we only stored metrics.
# However, we can recompute quickly by re‑reading and processing just the fill counts, without redoing all metrics.
# To save time, we can compute fill fractions in the same loop above. We should have accumulated fill counts per site.
# Let's restructure: In the loop above, we can also accumulate fill counts.
# Since we already finished the loop, we need to recompute fill fractions by re‑reading the data.
# But this adds extra time. Since this is a one‑time processing, it's acceptable.
# We'll re‑read the merged and blended files again, but only for the final domain, to count fills.

print("Re‑reading data for fill fractions...")
fill_data = []
for site in final_sites:
    merged_file = os.path.join(MERGED_DIR, f"merged_{site}.csv")
    blended_file = os.path.join(BLENDED_DIR, f"gaps_blend_{site}.csv")
    if not os.path.exists(merged_file) or not os.path.exists(blended_file):
        continue
    merged = pd.read_csv(merged_file)
    blended = pd.read_csv(blended_file)
    merged['DateTime'] = pd.to_datetime(merged['DateTime'])
    blended['DateTime'] = pd.to_datetime(blended['DateTime'])
    # Ensure columns
    for col in ['Tair', 'Rg', 'VPD']:
        if col not in merged.columns:
            merged[col] = np.nan
    if site in vpd_correction_sites:
        merged['VPD'] = merged['VPD'] * 10
    # Ensure ERA5 columns in merged
    if 'Tair_era' not in merged.columns:
        if 't2m' in merged.columns:
            merged['Tair_era'] = merged['t2m'] - 273.15
        else:
            merged['Tair_era'] = np.nan
    if 'Rg_era' not in merged.columns:
        if 'ssrd' in merged.columns:
            merged['Rg_era'] = merged['ssrd'] / 3600.0
        else:
            merged['Rg_era'] = np.nan
    if 'VPD_era' not in merged.columns:
        if 'd2m' in merged.columns and 't2m' in merged.columns:
            es = 0.6108 * np.exp((17.27 * (merged['t2m']-273.15)) / ((merged['t2m']-273.15) + 237.3))
            ea = 0.6108 * np.exp((17.27 * (merged['d2m']-273.15)) / ((merged['d2m']-273.15) + 237.3))
            merged['VPD_era'] = (es - ea) * 10.0
        else:
            merged['VPD_era'] = np.nan
    # Ensure corrected in blended
    for var in ['Tair', 'Rg', 'VPD']:
        col_c = f'{var}_era_c'
        if col_c not in blended.columns:
            blended[col_c] = merged[f'{var}_era'].copy()
    # Restrict to final domain
    merged['year_month'] = merged['DateTime'].dt.to_period('M')
    blended['year_month'] = blended['DateTime'].dt.to_period('M')
    domain_years_months = final_domain[final_domain['site_name'] == site][['Year','month']]
    domain_years_months['year_month'] = pd.to_datetime(
        domain_years_months['Year'].astype(str) + '-' + domain_years_months['month'].astype(str) + '-01'
    ).dt.to_period('M')
    domain_set = set(domain_years_months['year_month'])
    merged_dom = merged[merged['year_month'].isin(domain_set)].copy()
    blended_dom = blended[blended['year_month'].isin(domain_set)].copy()
    if merged_dom.empty or blended_dom.empty:
        continue
    merged_dom = merged_dom.set_index('DateTime').sort_index()
    blended_dom = blended_dom.set_index('DateTime').sort_index()
    common_idx = merged_dom.index.intersection(blended_dom.index)
    merged_dom = merged_dom.loc[common_idx]
    blended_dom = blended_dom.loc[common_idx]
    combined = pd.DataFrame(index=common_idx)
    combined['Tair_tower'] = merged_dom['Tair']
    combined['Rg_tower'] = merged_dom['Rg']
    combined['VPD_tower'] = merged_dom['VPD']
    combined['Tair_era'] = merged_dom['Tair_era']
    combined['Rg_era'] = merged_dom['Rg_era']
    combined['VPD_era'] = merged_dom['VPD_era']
    combined['Tair_era_c'] = blended_dom['Tair_era_c']
    combined['Rg_era_c'] = blended_dom['Rg_era_c']
    combined['VPD_era_c'] = blended_dom['VPD_era_c']
    combined['Tair_final'] = blended_dom['Tair']
    if 'RH' in merged_dom.columns:
        combined['RH'] = merged_dom['RH']
    else:
        combined['RH'] = np.nan

    total = len(combined)
    coast = final_domain[final_domain['site_name'] == site]['coast_region'].iloc[0]

    # Tair fill
    tair_missing = combined['Tair_tower'].isna()
    tair_era_avail = combined['Tair_era'].notna() | combined['Tair_era_c'].notna()
    tair_fill = tair_missing & tair_era_avail
    tair_fill_count = tair_fill.sum()

    # Rg fill
    rg_missing = combined['Rg_tower'].isna()
    rg_era_avail = combined['Rg_era'].notna() | combined['Rg_era_c'].notna()
    rg_fill = rg_missing & rg_era_avail
    rg_fill_count = rg_fill.sum()

    # VPD fill
    vpd_missing = combined['VPD_tower'].isna()
    tair_avail = combined['Tair_final'].notna()
    rh_avail = combined['RH'].notna()
    tair_rh_derived = vpd_missing & tair_avail & rh_avail
    vpd_era_avail = combined['VPD_era'].notna() | combined['VPD_era_c'].notna()
    vpd_era_fill = vpd_missing & ~tair_rh_derived & vpd_era_avail
    vpd_fill_count = vpd_era_fill.sum()

    fill_data.append({
        'site': site,
        'coast': coast,
        'total_records': total,
        'Tair_fill_count': tair_fill_count,
        'Rg_fill_count': rg_fill_count,
        'VPD_fill_count': vpd_fill_count,
        'Tair_fill_frac': tair_fill_count / total if total > 0 else np.nan,
        'Rg_fill_frac': rg_fill_count / total if total > 0 else np.nan,
        'VPD_fill_frac': vpd_fill_count / total if total > 0 else np.nan,
    })

df_fill = pd.DataFrame(fill_data)

# Merge fill fractions into metrics dataframes
df_tair = df_tair.merge(df_fill[['site', 'Tair_fill_frac', 'Tair_fill_count', 'total_records']],
                        left_on='site', right_on='site', how='left')
df_tair.rename(columns={'Tair_fill_frac': 'fill_frac', 'Tair_fill_count': 'era_fill_count'}, inplace=True)

df_rg = df_rg.merge(df_fill[['site', 'Rg_fill_frac', 'Rg_fill_count', 'total_records']],
                    left_on='site', right_on='site', how='left')
df_rg.rename(columns={'Rg_fill_frac': 'fill_frac', 'Rg_fill_count': 'era_fill_count'}, inplace=True)

# ============================================================================
# SAVE INTERMEDIATE DATA TO CACHE
# ============================================================================
print("\nSaving intermediate data to cache...")
# Save DataFrames
df_tair.to_csv(os.path.join(CACHE_DIR, 'df_tair.csv'), index=False)
df_rg.to_csv(os.path.join(CACHE_DIR, 'df_rg.csv'), index=False)
df_vpd.to_csv(os.path.join(CACHE_DIR, 'df_vpd.csv'), index=False)
df_fill.to_csv(os.path.join(CACHE_DIR, 'df_fill.csv'), index=False)
final_domain.to_csv(os.path.join(CACHE_DIR, 'final_domain.csv'), index=False)

# Save pooled scatter arrays as numpy .npy
for key, arr in pooled.items():
    np.save(os.path.join(CACHE_DIR, f'{key}.npy'), np.array(arr))

# Also save coast counts and other QC info for the figure? Not needed.

print("Chunk 1 processing complete. You can now run Chunk 2 to generate the figure.")