# -*- coding: utf-8 -*-
"""
Q3 WORKFLOW: WUE_T vs SPEI-6 Timescale Analysis
UPDATED: Uses ALL SPEI-6 values (not just NN) - only requires WUE_T and SPEI_6 columns
UPDATED: Strict triple intersection (57 sites)
FIXED: Minimum 5 observations per site (any SPEI value, not just NN)
FIXED: No filtering by SPEI category - uses full SPEI gradient
FIXED: Drop duplicates from metadata before merge
FIXED: Uses len(site_data) not drop_duplicates() for observation counting
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
import os

# =============================================================================
# CONFIGURATION
# =============================================================================

# File paths
BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
MONTHLY_FILE = os.path.join(BASE_DIR, "monthly_data_after_outlier_removal.csv")
NN_MEDIANS_FILE = os.path.join(BASE_DIR, "wue_site_level_NN_medians_SPEI_1.csv")
METADATA_FILE = os.path.join(BASE_DIR, "site_metadata_with_salinity_SPEIinfo.csv")

# Output directory
OUTPUT_DIR = os.path.join(BASE_DIR, "Q3_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output files
OUTPUT_SLOPE_CSV = os.path.join(OUTPUT_DIR, "Q3_WUET_SPEI6_slopes.csv")
OUTPUT_FIGURE_PNG = os.path.join(OUTPUT_DIR, "Q3_WUET_SPEI6_scatter.png")
OUTPUT_FIGURE_PDF = os.path.join(OUTPUT_DIR, "Q3_WUET_SPEI6_scatter.pdf")

# SPEI timescale for Q3
SPEI_TIMESCALE = 'SPEI_6'

# Minimum observations required per site for slope calculation
MIN_OBSERVATIONS = 5

# =============================================================================
# STEP 1: GET STRICT TRIPLE INTERSECTION SITES (57 sites)
# =============================================================================

print("="*80)
print("Q3 WORKFLOW: WUE_T vs SPEI-6 Timescale Analysis")
print("Using ALL SPEI-6 values (full gradient, not just NN)")
print(f"Minimum observations per site: {MIN_OBSERVATIONS}")
print("="*80)

# Load NN medians file
nn_medians = pd.read_csv(NN_MEDIANS_FILE)
print(f"\nLoaded NN medians file: {len(nn_medians)} rows")

# Get sites with ALL THREE metrics (strict triple intersection)
wue_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())

shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
print(f"Strict triple intersection sites: {len(shared_sites)}")

# =============================================================================
# STEP 2: LOAD MONTHLY DATA
# =============================================================================

print("\n" + "="*60)
print("LOADING MONTHLY DATA")
print("="*60)

df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"Loaded monthly data: {len(df_monthly):,} rows")

# IMPORTANT: Do NOT filter by SPEI category - use ALL SPEI values
# Only filter to strict triple intersection sites
df_filtered = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"After strict triple intersection filter: {len(df_filtered):,} rows")
print(f"Unique sites: {df_filtered['site_name'].nunique()} (expected: {len(shared_sites)})")

# =============================================================================
# STEP 3: APPLY STRICT MONTH FILTER (ALL THREE METRICS MUST HAVE DATA)
# =============================================================================

print("\n" + "="*60)
print("APPLYING STRICT MONTH FILTER")
print("="*60)

strict_mask = (
    df_filtered['WUE'].notna() &
    df_filtered['WUE_eva'].notna() &
    df_filtered['WUE_tra'].notna()
)

df_strict = df_filtered[strict_mask].copy()
print(f"After strict month filter: {len(df_strict):,} rows")
print(f"Rows removed: {len(df_filtered) - len(df_strict)} ({((len(df_filtered) - len(df_strict))/len(df_filtered)*100):.1f}%)")

# =============================================================================
# STEP 4: PREPARE DATA FOR Q3 ANALYSIS (WUE_T vs SPEI_6)
# =============================================================================

print("\n" + "="*60)
print("PREPARING Q3 DATA")
print("="*60)

# Rename columns for clarity
df_q3 = df_strict.rename(columns={
    'WUE_tra': 'WUE_T',
    'SPEI_6': 'SPEI_6'
})

# CRITICAL: Only require WUE_T and SPEI_6 for slope calculation
# Do NOT drop rows just because Trans_ratio, lat, or long are missing
df_q3_clean = df_q3.dropna(subset=['WUE_T', 'SPEI_6']).copy()
print(f"Rows after dropping missing WUE_T or SPEI_6: {len(df_q3_clean):,}")

# Convert numeric columns (coerce errors to NaN, but don't drop yet)
for col in ['WUE_T', 'SPEI_6', 'Trans_ratio', 'lat', 'long']:
    if col in df_q3_clean.columns:
        df_q3_clean[col] = pd.to_numeric(df_q3_clean[col], errors='coerce')

# =============================================================================
# STEP 5: LOAD METADATA AND REMOVE DUPLICATES
# =============================================================================

print("\n" + "="*60)
print("LOADING METADATA")
print("="*60)

metadata = pd.read_csv(METADATA_FILE)

# CRITICAL FIX: Remove duplicate site_name entries before merging
if metadata['site_name'].duplicated().any():
    print(f"⚠️ Found {metadata['site_name'].duplicated().sum()} duplicate site_name entries")
    metadata = metadata.drop_duplicates(subset=['site_name'], keep='first')
    print(f"   Removed duplicates, now {len(metadata)} unique sites")

# Keep only relevant columns
metadata_cols = ['site_name', 'Salinity_Category', 'water_class', 'IGBP', 'climate', 'lat', 'long']
available_cols = [col for col in metadata_cols if col in metadata.columns]
metadata = metadata[available_cols].copy()
print(f"Metadata loaded: {len(metadata)} unique sites")

# =============================================================================
# STEP 6: CALCULATE SITE-LEVEL SLOPES
# =============================================================================

print("\n" + "="*60)
print(f"CALCULATING SITE-LEVEL SLOPES (WUE_T vs SPEI_6)")
print(f"Minimum observations required: {MIN_OBSERVATIONS}")
print("="*60)

site_results = []
sites_excluded = []

for site in shared_sites:
    # Get data for this site
    site_data = df_q3_clean[df_q3_clean['site_name'] == site].copy()
    
    if len(site_data) == 0:
        sites_excluded.append({'site': site, 'reason': 'No data available', 'n_obs': 0})
        continue
    
    # Count observations (FIXED: use len, not drop_duplicates)
    n_obs = len(site_data)
    
    # Calculate median T:ET ratio (use dropna, but don't filter entire site if missing)
    median_TET = site_data['Trans_ratio'].dropna().median()
    
    # Extract WUE_T and SPEI_6 for slope calculation
    wue_t_values = site_data['WUE_T'].values
    spei_values = site_data['SPEI_6'].values
    
    # Calculate slope using linear regression (requires minimum observations)
    if len(wue_t_values) >= MIN_OBSERVATIONS:
        slope, intercept, r_value, p_value, std_err = linregress(spei_values, wue_t_values)
        r_squared = r_value**2
        
        # Calculate Spearman correlation
        spearman_rho, spearman_p = spearmanr(spei_values, wue_t_values)
        
        site_results.append({
            'site_name': site,
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'pearson_p': p_value,
            'spearman_rho': spearman_rho,
            'spearman_p': spearman_p,
            'n_observations': n_obs,
            'median_TET': median_TET
        })
        
        print(f"  {site}: n={n_obs}, slope={slope:.4f}, r²={r_squared:.3f}, ρ={spearman_rho:.3f}")
    else:
        sites_excluded.append({'site': site, 'reason': f'Insufficient data (<{MIN_OBSERVATIONS} obs)', 'n_obs': n_obs})
        print(f"  {site}: EXCLUDED - n={n_obs} (<{MIN_OBSERVATIONS})")

# =============================================================================
# STEP 7: CREATE RESULTS DATAFRAME AND MERGE METADATA
# =============================================================================

results_df = pd.DataFrame(site_results)
print(f"\nTotal sites with valid slopes: {len(results_df)}")
print(f"Total sites excluded: {len(sites_excluded)}")

# Merge metadata (safe merge since duplicates removed)
if len(results_df) > 0:
    results_df = results_df.merge(metadata, on='site_name', how='left')
    
    # Reorder columns for better readability
    col_order = ['site_name', 'Salinity_Category', 'slope', 'r_squared', 'spearman_rho', 
                 'pearson_p', 'spearman_p', 'n_observations', 'median_TET', 'intercept',
                 'water_class', 'IGBP', 'climate', 'lat', 'long']
    available_cols = [col for col in col_order if col in results_df.columns]
    results_df = results_df[available_cols]

# Save results
results_df.to_csv(OUTPUT_SLOPE_CSV, index=False)
print(f"\n✅ Saved slopes to: {OUTPUT_SLOPE_CSV}")

# Save excluded sites for reference
if sites_excluded:
    excluded_df = pd.DataFrame(sites_excluded)
    excluded_path = os.path.join(OUTPUT_DIR, "Q3_excluded_sites.csv")
    excluded_df.to_csv(excluded_path, index=False)
    print(f"✅ Saved excluded sites to: {excluded_path}")

# =============================================================================
# STEP 8: PRINT SUMMARY STATISTICS
# =============================================================================

print("\n" + "="*80)
print("Q3 SUMMARY STATISTICS")
print("="*80)

if len(results_df) > 0:
    print(f"\nSites analyzed: {len(results_df)} / {len(shared_sites)} ({len(results_df)/len(shared_sites)*100:.1f}%)")
    print(f"  Sites with positive slope: {(results_df['slope'] > 0).sum()}")
    print(f"  Sites with negative slope: {(results_df['slope'] < 0).sum()}")
    print(f"  Sites with zero slope: {(results_df['slope'] == 0).sum()}")
    
    print(f"\nSlope statistics:")
    print(f"  Mean slope: {results_df['slope'].mean():.4f}")
    print(f"  Median slope: {results_df['slope'].median():.4f}")
    print(f"  Std slope: {results_df['slope'].std():.4f}")
    print(f"  Min slope: {results_df['slope'].min():.4f}")
    print(f"  Max slope: {results_df['slope'].max():.4f}")
    
    print(f"\nR² statistics:")
    print(f"  Mean R²: {results_df['r_squared'].mean():.3f}")
    print(f"  Median R²: {results_df['r_squared'].median():.3f}")
    
    print(f"\nSpearman correlation statistics:")
    print(f"  Mean ρ: {results_df['spearman_rho'].mean():.3f}")
    print(f"  Median ρ: {results_df['spearman_rho'].median():.3f}")
    
    print(f"\nObservation counts:")
    print(f"  Total observations across all sites: {results_df['n_observations'].sum()}")
    print(f"  Mean observations per site: {results_df['n_observations'].mean():.1f}")
    print(f"  Median observations per site: {results_df['n_observations'].median():.0f}")
    print(f"  Min observations per site: {results_df['n_observations'].min()}")
    print(f"  Max observations per site: {results_df['n_observations'].max()}")
    
    # Breakdown by salinity category
    if 'Salinity_Category' in results_df.columns:
        print(f"\nBreakdown by salinity category:")
        for cat in ['Freshwater', 'Saline', 'Upland', 'Brackish']:
            subset = results_df[results_df['Salinity_Category'] == cat]
            if len(subset) > 0:
                print(f"  {cat}: {len(subset)} sites, median slope={subset['slope'].median():.4f}")
else:
    print("\n⚠️ No sites met the minimum observation requirement!")

# Print excluded sites summary
if sites_excluded:
    print(f"\nExcluded sites summary:")
    excluded_by_reason = {}
    for ex in sites_excluded:
        reason = ex['reason']
        excluded_by_reason[reason] = excluded_by_reason.get(reason, 0) + 1
    for reason, count in excluded_by_reason.items():
        print(f"  {reason}: {count} sites")

# =============================================================================
# STEP 9: CREATE SCATTER PLOT (OPTIONAL)
# =============================================================================

print("\n" + "="*80)
print("CREATING SCATTER PLOT")
print("="*80)

fig, ax = plt.subplots(figsize=(10, 8))

if len(results_df) > 0:
    # Color by slope magnitude
    slopes = results_df['slope'].values
    norm = plt.Normalize(vmin=min(slopes), vmax=max(slopes))
    cmap = plt.cm.RdYlBu_r
    
    scatter = ax.scatter(results_df['spearman_rho'], results_df['slope'],
                        c=slopes, cmap=cmap, norm=norm,
                        s=100, alpha=0.7, edgecolors='black', linewidth=1)
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Slope (WUE_T / SPEI-6)', fontsize=14)
    
    # Add reference lines
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Labels
    ax.set_xlabel('Spearman Correlation (ρ)', fontsize=16)
    ax.set_ylabel('Linear Slope (WUE_T vs SPEI-6)', fontsize=16)
    ax.set_title('Q3: WUE_T Sensitivity to SPEI-6 (Full Gradient)', fontsize=18, fontweight='bold')
    
    # Add text annotation for summary
    ax.text(0.05, 0.95, f'n = {len(results_df)} sites\n(min {MIN_OBSERVATIONS} obs/site)',
            transform=ax.transAxes, fontsize=14, va='top')
    
    # Save figure
    plt.tight_layout()
    plt.savefig(OUTPUT_FIGURE_PNG, dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_FIGURE_PDF, bbox_inches='tight')
    print(f"\n✅ Saved scatter plot to: {OUTPUT_FIGURE_PNG}")
    print(f"✅ Saved scatter plot to: {OUTPUT_FIGURE_PDF}")
    
    plt.show()
else:
    print("\n⚠️ No data to plot")
    plt.close()

# =============================================================================
# STEP 10: FINAL SUMMARY
# =============================================================================

print("\n" + "="*80)
print("Q3 WORKFLOW COMPLETE")
print("="*80)
print("\n✅ KEY UPDATES:")
print("  1. Uses ALL SPEI-6 values (full gradient, not filtered by NN)")
print("  2. Strict triple intersection sites: 57 sites")
print(f"  3. Minimum {MIN_OBSERVATIONS} observations per site (any SPEI value)")
print("  4. Only drops rows missing WUE_T or SPEI_6 (slope calculation)")
print("  5. Trans_ratio missing values handled gracefully (dropped only for median)")
print("  6. Latitude/longitude not required - sites kept even if missing")
print("  7. Fixed metadata duplicate issue - removed duplicates before merge")
print("  8. Month counting uses len(site_data) (not drop_duplicates)")
print(f"\n📊 Results saved to: {OUTPUT_DIR}")
print(f"   - Slope data: {os.path.basename(OUTPUT_SLOPE_CSV)}")
print(f"   - Excluded sites: Q3_excluded_sites.csv")
print(f"   - Scatter plot: {os.path.basename(OUTPUT_FIGURE_PNG)}")