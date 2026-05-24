# -*- coding: utf-8 -*-
"""
DIAGNOSTIC ANALYSIS FOR FIGURE 3: PRE vs POST vs Δ (POST-PRE)
Analyzes frac_less values before and after breakpoint.
Produces statistical summaries with appropriate tests for spatial data.

Interpretation:
  PRE (≤ Breakpoint) : Pre-breakpoint period
  POST (> Breakpoint): Post-breakpoint period  
  Δ = POST - PRE    : Change in frequency of WUE_T decline

Positive Δ  → Increased frequency of WUE_T decline
Negative Δ  → Decreased frequency of WUE_T decline

IMPORTANT STATISTICAL NOTE:
- Pixel-level tests treat each spatial pixel as independent
- Spatial autocorrelation exists (p-values may be optimistic)
- Results should be interpreted as descriptive patterns
- Wilcoxon signed-rank test is more robust than t-test for these data
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION - UPDATED FOR WUE_T
# =============================================================================
PRE_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_pre_post_breakpoint\agg_PRE_WUE_tra.nc"
POST_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_pre_post_breakpoint\agg_POST_WUE_tra.nc"

# Statistical significance threshold
ALPHA = 0.05

# Minimum pixels for statistical tests
MIN_PIXELS_FOR_TESTS = 20

# Change magnitude thresholds
THRESH_LARGE_POSITIVE = 0.15   # Substantial increase
THRESH_LARGE_NEGATIVE = -0.15  # Substantial decrease
THRESH_SMALL_CHANGE = 0.05     # Minor/noise-level change

# Map extent (matches Figure 3 - UPDATED)
CONUS_EXTENT = [-126, -66, 24, 50]
ALASKA_EXTENT = [-165, -130, 52, 72]

# =============================================================================
# COASTAL REGION DEFINITIONS - UPDATED TO 5 REGIONS (Matching Figure 2)
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
    "CONUS Coastal Domain": {
        "lon_bounds": CONUS_EXTENT[:2],
        "lat_bounds": CONUS_EXTENT[2:],
        "states": ["All coastal pixels"]
    }
}

# =============================================================================
# DATA LOADING AND PROCESSING
# =============================================================================
def load_and_align_data(pre_file, post_file):
    """Load PRE and POST data and ensure spatial alignment."""
    print("Loading and aligning data...")
    
    # Load PRE data
    ds_pre = xr.open_dataset(pre_file)
    
    # Validate required variable
    if "frac_less" not in ds_pre.data_vars:
        raise ValueError("'frac_less' not found in PRE dataset")
    
    if np.nanmax(ds_pre.lon.values) > 180:
        ds_pre = ds_pre.assign_coords(lon=((ds_pre.lon + 180) % 360) - 180).sortby('lon')
    
    # Load POST data
    ds_post = xr.open_dataset(post_file)
    
    # Validate required variable
    if "frac_less" not in ds_post.data_vars:
        raise ValueError("'frac_less' not found in POST dataset")
    
    if np.nanmax(ds_post.lon.values) > 180:
        ds_post = ds_post.assign_coords(lon=((ds_post.lon + 180) % 360) - 180).sortby('lon')
    
    # Ensure consistent dimensions
    if ds_pre.lat.shape != ds_post.lat.shape or ds_pre.lon.shape != ds_post.lon.shape:
        print("  Interpolating POST to PRE grid...")
        ds_post = ds_post.interp(lat=ds_pre.lat, lon=ds_pre.lon)
    
    # Extract data
    pre_data = ds_pre["frac_less"].values
    post_data = ds_post["frac_less"].values
    lats = ds_pre.lat.values
    lons = ds_pre.lon.values
    
    # Create 2D grids if needed
    if lats.ndim == 1 and lons.ndim == 1:
        lon_grid, lat_grid = np.meshgrid(lons, lats)
    else:
        lon_grid, lat_grid = lons, lats
    
    ds_pre.close()
    ds_post.close()
    
    print(f"  Data shape: {pre_data.shape}")
    print(f"  Valid PRE pixels: {np.sum(np.isfinite(pre_data)):,}")
    print(f"  Valid POST pixels: {np.sum(np.isfinite(post_data)):,}")
    
    return pre_data, post_data, lat_grid, lon_grid

def calculate_delta(pre_data, post_data):
    """Calculate change (Δ = POST - PRE)."""
    mask = np.isfinite(pre_data) & np.isfinite(post_data)
    delta = np.full_like(pre_data, np.nan)
    delta[mask] = post_data[mask] - pre_data[mask]
    
    valid_count = np.sum(mask)
    print(f"  Valid Δ pixels (both periods): {valid_count:,}")
    print(f"    ({100*valid_count/pre_data.size:.1f}% of total grid)")
    
    return delta, mask

def create_region_mask(lon_grid, lat_grid, region_def):
    """Create mask for a region."""
    mask = (
        (lon_grid >= region_def["lon_bounds"][0]) & 
        (lon_grid <= region_def["lon_bounds"][1]) &
        (lat_grid >= region_def["lat_bounds"][0]) & 
        (lat_grid <= region_def["lat_bounds"][1])
    )
    return mask

# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================
def analyze_region(region_name, region_def, pre_data, post_data, delta_data, mask_both):
    """Comprehensive analysis for a single region."""
    lon_grid = region_def["_lon_grid"]
    lat_grid = region_def["_lat_grid"]
    
    # Create region mask
    region_mask = create_region_mask(lon_grid, lat_grid, region_def)
    
    # Combine with data validity mask
    valid_mask = region_mask & mask_both
    
    n_pixels = int(np.sum(valid_mask))
    
    if n_pixels == 0:
        return None
    
    # Extract valid data
    pre_vals = pre_data[valid_mask]
    post_vals = post_data[valid_mask]
    delta_vals = delta_data[valid_mask]
    
    # Basic statistics
    stats_dict = {
        "Region": region_name,
        "States": ", ".join(region_def["states"]),
        "n_pixels": n_pixels,
        
        # PRE statistics
        "PRE_mean": float(np.nanmean(pre_vals)),
        "PRE_median": float(np.nanmedian(pre_vals)),
        "PRE_std": float(np.nanstd(pre_vals)),
        "PRE_q25": float(np.nanpercentile(pre_vals, 25)),
        "PRE_q75": float(np.nanpercentile(pre_vals, 75)),
        "PRE_IQR": float(np.nanpercentile(pre_vals, 75) - np.nanpercentile(pre_vals, 25)),
        
        # POST statistics
        "POST_mean": float(np.nanmean(post_vals)),
        "POST_median": float(np.nanmedian(post_vals)),
        "POST_std": float(np.nanstd(post_vals)),
        "POST_q25": float(np.nanpercentile(post_vals, 25)),
        "POST_q75": float(np.nanpercentile(post_vals, 75)),
        "POST_IQR": float(np.nanpercentile(post_vals, 75) - np.nanpercentile(post_vals, 25)),
        
        # Δ statistics
        "Δ_mean": float(np.nanmean(delta_vals)),
        "Δ_median": float(np.nanmedian(delta_vals)),
        "Δ_std": float(np.nanstd(delta_vals)),
        "Δ_q25": float(np.nanpercentile(delta_vals, 25)),
        "Δ_q75": float(np.nanpercentile(delta_vals, 75)),
        "Δ_IQR": float(np.nanpercentile(delta_vals, 75) - np.nanpercentile(delta_vals, 25)),
        "Δ_min": float(np.nanmin(delta_vals)),
        "Δ_max": float(np.nanmax(delta_vals)),
    }
    
    # Statistical tests
    if n_pixels >= MIN_PIXELS_FOR_TESTS:
        # Paired t-test
        try:
            t_stat, p_value_t = stats.ttest_rel(post_vals, pre_vals, nan_policy='omit')
            stats_dict["t_statistic"] = float(t_stat)
            stats_dict["p_value_t"] = float(p_value_t)
            stats_dict["significant_t"] = p_value_t < ALPHA
        except:
            stats_dict["t_statistic"] = np.nan
            stats_dict["p_value_t"] = np.nan
            stats_dict["significant_t"] = False
        
        # Wilcoxon signed-rank test
        finite_mask = np.isfinite(pre_vals) & np.isfinite(post_vals)
        x = pre_vals[finite_mask]
        y = post_vals[finite_mask]
        
        if len(x) >= MIN_PIXELS_FOR_TESTS:
            try:
                w_stat, p_value_w = stats.wilcoxon(y, x, zero_method='wilcox')
                stats_dict["wilcoxon_stat"] = float(w_stat)
                stats_dict["p_value_wilcoxon"] = float(p_value_w)
                stats_dict["significant_wilcoxon"] = p_value_w < ALPHA
            except:
                stats_dict["wilcoxon_stat"] = np.nan
                stats_dict["p_value_wilcoxon"] = np.nan
                stats_dict["significant_wilcoxon"] = False
        else:
            stats_dict["wilcoxon_stat"] = np.nan
            stats_dict["p_value_wilcoxon"] = np.nan
            stats_dict["significant_wilcoxon"] = False
    else:
        stats_dict["t_statistic"] = np.nan
        stats_dict["p_value_t"] = np.nan
        stats_dict["wilcoxon_stat"] = np.nan
        stats_dict["p_value_wilcoxon"] = np.nan
        stats_dict["significant_t"] = False
        stats_dict["significant_wilcoxon"] = False
        stats_dict["statistical_note"] = f"Insufficient pixels (n={n_pixels})"
    
    # Direction and magnitude of change
    stats_dict["Δ_direction"] = "increase" if stats_dict["Δ_mean"] > 0 else "decrease"
    
    # Proportion of pixels with significant changes
    large_pos = np.sum(delta_vals >= THRESH_LARGE_POSITIVE)
    large_neg = np.sum(delta_vals <= THRESH_LARGE_NEGATIVE)
    small_change = np.sum(np.abs(delta_vals) < THRESH_SMALL_CHANGE)
    
    stats_dict["pct_large_pos"] = float(100 * large_pos / n_pixels)
    stats_dict["pct_large_neg"] = float(100 * large_neg / n_pixels)
    stats_dict["pct_small_change"] = float(100 * small_change / n_pixels)
    
    # Calculate percentage change
    if stats_dict["PRE_mean"] != 0 and np.isfinite(stats_dict["PRE_mean"]):
        pct_change = 100 * stats_dict["Δ_mean"] / stats_dict["PRE_mean"]
        stats_dict["pct_change"] = float(pct_change)
    else:
        stats_dict["pct_change"] = np.nan
    
    return stats_dict

# =============================================================================
# REPORT GENERATION
# =============================================================================
def print_statistical_caveats():
    """Print important caveats about statistical tests."""
    print("\n" + "="*100)
    print("STATISTICAL METHODS AND LIMITATIONS")
    print("="*100)
    print("\nMETHODS:")
    print("1. Paired t-test: Assumes normal distribution of differences")
    print("2. Wilcoxon signed-rank test: Non-parametric, more robust to non-normality")
    print("3. Tests performed only for regions with ≥20 valid pixels")
    print("4. Percentage changes calculated relative to PRE period means")
    
    print("\nIMPORTANT LIMITATIONS:")
    print("1. Spatial autocorrelation: Pixels are not fully independent,")
    print("   p-values may be optimistic (Type I error risk)")
    print("2. Multiple testing: No correction for multiple comparisons across regions")
    print("3. Results are descriptive patterns, not causal evidence")
    print("4. Statistical significance (p < 0.05) ≠ ecological significance")
    
    print("\nRECOMMENDED INTERPRETATION:")
    print("• Focus on effect sizes (Δ means, percentage changes)")
    print("• Consider Wilcoxon results more robust than t-test")
    print("• Emphasize patterns over p-values in ecological interpretation")
    print()

def print_overall_summary(df):
    """Print overall summary of results."""
    print("\n" + "="*100)
    print("OVERALL SUMMARY: PRE vs POST BREAKPOINT")
    print("="*100)
    
    conus_data = df[df['Region'] == 'CONUS Coastal Domain'].iloc[0]
    
    print(f"\nCONUS COASTAL DOMAIN CHANGES:")
    print(f"  • PRE period (≤Breakpoint):  {conus_data['PRE_mean']:.3f} (IQR: {conus_data['PRE_IQR']:.3f})")
    print(f"  • POST period (>Breakpoint): {conus_data['POST_mean']:.3f} (IQR: {conus_data['POST_IQR']:.3f})")
    print(f"  • Absolute change (Δ): {conus_data['Δ_mean']:+.3f} (IQR: {conus_data['Δ_IQR']:.3f})")
    print(f"  • Relative change: {conus_data['pct_change']:+.1f}%")
    
    if pd.notnull(conus_data['p_value_wilcoxon']):
        print(f"  • Wilcoxon test: p = {conus_data['p_value_wilcoxon']:.4f} " +
              f"{'(significant)' if conus_data['significant_wilcoxon'] else '(not significant)'}")
    
    print(f"\nSPATIAL PATTERNS OF CHANGE:")
    print(f"  • Pixels with large increase (Δ ≥ {THRESH_LARGE_POSITIVE}): {conus_data['pct_large_pos']:.1f}%")
    print(f"  • Pixels with large decrease (Δ ≤ {THRESH_LARGE_NEGATIVE}): {conus_data['pct_large_neg']:.1f}%")
    print(f"  • Pixels with minimal change (|Δ| < {THRESH_SMALL_CHANGE}): {conus_data['pct_small_change']:.1f}%")

def print_regional_comparison(df):
    """Print detailed regional comparison."""
    print("\n" + "="*100)
    print("REGIONAL COMPARISONS")
    print("="*100)
    
    # Filter out CONUS Total and Alaska for ranking
    df_regions = df[~df['Region'].isin(['CONUS Coastal Domain', 'Alaska'])].copy()
    df_regions = df_regions.sort_values('Δ_mean', ascending=False)
    df_regions['Rank'] = range(1, len(df_regions) + 1)
    
    # Most and least changed regions
    most_increased = df_regions.iloc[0]
    most_decreased = df_regions.iloc[-1]
    
    print(f"\n1. GREATEST INCREASE IN WUE_T DECLINE FREQUENCY:")
    print(f"   Region: {most_increased['Region']}")
    print(f"   • Δ = {most_increased['Δ_mean']:+.3f} ({most_increased['pct_change']:+.1f}% change)")
    print(f"   • PRE: {most_increased['PRE_mean']:.3f} → POST: {most_increased['POST_mean']:.3f}")
    print(f"   • Pixels analyzed: {most_increased['n_pixels']:,}")
    print(f"   • {most_increased['pct_large_pos']:.1f}% of pixels show large increases")
    
    if pd.notnull(most_increased['p_value_wilcoxon']):
        sig_status = "significant" if most_increased['significant_wilcoxon'] else "not significant"
        print(f"   • Wilcoxon test: p = {most_increased['p_value_wilcoxon']:.4f} ({sig_status})")
    
    print(f"\n2. GREATEST DECREASE IN WUE_T DECLINE FREQUENCY:")
    print(f"   Region: {most_decreased['Region']}")
    print(f"   • Δ = {most_decreased['Δ_mean']:+.3f} ({most_decreased['pct_change']:+.1f}% change)")
    print(f"   • PRE: {most_decreased['PRE_mean']:.3f} → POST: {most_decreased['POST_mean']:.3f}")
    print(f"   • Pixels analyzed: {most_decreased['n_pixels']:,}")
    print(f"   • {most_decreased['pct_large_neg']:.1f}% of pixels show large decreases")
    
    if pd.notnull(most_decreased['p_value_wilcoxon']):
        sig_status = "significant" if most_decreased['significant_wilcoxon'] else "not significant"
        print(f"   • Wilcoxon test: p = {most_decreased['p_value_wilcoxon']:.4f} ({sig_status})")
    
    # Coastal comparisons
    print(f"\n3. COASTAL COMPARISONS:")
    
    # Pacific vs Atlantic North
    pacific_data = df_regions[df_regions['Region'] == 'Pacific']
    atlantic_north_data = df_regions[df_regions['Region'] == 'Atlantic North']
    
    if len(pacific_data) > 0:
        pacific_value = pacific_data['Δ_mean'].iloc[0]
        print(f"   • Pacific Coast Δ: {pacific_value:+.3f}")
    
    if len(atlantic_north_data) > 0:
        atlantic_north_value = atlantic_north_data['Δ_mean'].iloc[0]
        print(f"   • Atlantic North Coast Δ: {atlantic_north_value:+.3f}")
    
    # Gulf vs Southeast Atlantic
    gulf_data = df_regions[df_regions['Region'] == 'Gulf']
    southeast_data = df_regions[df_regions['Region'] == 'Southeast Atlantic']
    
    if len(gulf_data) > 0:
        gulf_value = gulf_data['Δ_mean'].iloc[0]
        print(f"\n   • Gulf of America Δ: {gulf_value:+.3f}")
    
    if len(southeast_data) > 0:
        southeast_value = southeast_data['Δ_mean'].iloc[0]
        print(f"   • Southeast Atlantic Δ: {southeast_value:+.3f}")

def print_nature_style_narrative(df):
    """Generate Nature-style narrative summary of results."""
    print("\n" + "="*100)
    print("NARRATIVE SUMMARY FOR RESEARCH PAPER")
    print("="*100)
    
    df_regions = df[~df['Region'].isin(['CONUS Coastal Domain', 'Alaska'])].copy()
    df_regions = df_regions.sort_values('Δ_mean', ascending=False)
    
    conus_data = df[df['Region'] == 'CONUS Coastal Domain'].iloc[0]
    most_increased = df_regions.iloc[0]
    most_decreased = df_regions.iloc[-1]
    
    # Pacific vs Atlantic North
    pacific_data = df_regions[df_regions['Region'] == 'Pacific']
    atlantic_north_data = df_regions[df_regions['Region'] == 'Atlantic North']
    
    pacific_value = pacific_data['Δ_mean'].iloc[0] if len(pacific_data) > 0 else None
    atlantic_north_value = atlantic_north_data['Δ_mean'].iloc[0] if len(atlantic_north_data) > 0 else None
    
    gulf_data = df_regions[df_regions['Region'] == 'Gulf']
    southeast_data = df_regions[df_regions['Region'] == 'Southeast Atlantic']
    
    gulf_value = gulf_data['Δ_mean'].iloc[0] if len(gulf_data) > 0 else None
    southeast_value = southeast_data['Δ_mean'].iloc[0] if len(southeast_data) > 0 else None
    
    print(f"\nThe frequency of WUE_T decline changed heterogeneously")
    print(f"across U.S. coastal ecosystems following the breakpoint. Mean frequencies")
    print(f"changed by Δ = {conus_data['Δ_mean']:.3f} ({conus_data['pct_change']:+.1f}%)")
    print(f"across the coastal domain, with {conus_data['pct_large_pos']:.1f}% of pixels")
    print(f"showing substantial increases (Δ ≥ {THRESH_LARGE_POSITIVE}) and")
    print(f"{conus_data['pct_large_neg']:.1f}% showing substantial decreases.")
    
    print(f"\nRegional changes varied markedly. The {most_increased['Region']} exhibited")
    print(f"the greatest increase (Δ = {most_increased['Δ_mean']:+.3f}, " +
          f"{most_increased['pct_change']:+.1f}%), rising from")
    print(f"{most_increased['PRE_mean']:.2f} to {most_increased['POST_mean']:.2f}. In contrast,")
    print(f"the {most_decreased['Region']} showed a decrease of Δ = {most_decreased['Δ_mean']:+.3f}")
    print(f"({most_decreased['pct_change']:+.1f}%), from {most_decreased['PRE_mean']:.2f} to")
    print(f"{most_decreased['POST_mean']:.2f}.")
    
    if pacific_value is not None and atlantic_north_value is not None:
        print(f"\nPacific Coast ecosystems experienced Δ = {pacific_value:+.3f}")
        print(f"while Atlantic North Coast regions showed Δ = {atlantic_north_value:+.3f}.")
    
    if gulf_value is not None:
        print(f"The Gulf of America showed Δ = {gulf_value:+.3f}.")
    
    if southeast_value is not None:
        print(f"The Southeast Atlantic showed Δ = {southeast_value:+.3f}.")
    
    # List regions with statistically significant changes
    sig_regions = df_regions[pd.notnull(df_regions['p_value_wilcoxon']) & 
                            df_regions['significant_wilcoxon']]['Region'].tolist()
    if sig_regions:
        print(f"\nNon-parametric Wilcoxon tests indicated statistically significant changes")
        print(f"(p < 0.05) in: {', '.join(sig_regions)}.")
    else:
        print(f"\nNo regions showed statistically significant changes at p < 0.05.")
    
    print(f"\nThese patterns demonstrate spatially heterogeneous post-breakpoint shifts in")
    print(f"WUE_T decline frequency across U.S. coastal ecosystems.")

def print_detailed_statistics(df):
    """Print detailed statistical tables."""
    print("\n" + "="*100)
    print("DETAILED STATISTICAL TABLES")
    print("="*100)
    
    print("\nTABLE 1: Regional Summary Statistics (sorted by Δ)")
    print("-" * 100)
    
    display_cols = ['Region', 'n_pixels', 'PRE_mean', 'POST_mean', 'Δ_mean', 
                   'Δ_IQR', 'pct_change', 'p_value_wilcoxon', 'significant_wilcoxon']
    
    display_df = df[display_cols].copy()
    display_df = display_df.sort_values('Δ_mean', ascending=False)
    
    def format_float(x):
        if isinstance(x, float):
            return f"{x:.3f}"
        return str(x)
    
    print(display_df.to_string(index=False, formatters={
        'n_pixels': lambda x: f"{x:,}",
        'PRE_mean': format_float,
        'POST_mean': format_float,
        'Δ_mean': lambda x: f"{x:+.3f}",
        'Δ_IQR': format_float,
        'pct_change': lambda x: f"{x:+.1f}%" if pd.notnull(x) else "N/A",
        'p_value_wilcoxon': lambda x: f"{x:.4f}" if pd.notnull(x) else "N/A",
        'significant_wilcoxon': lambda x: "Yes" if x else "No"
    }))
    
    print("\n\nTABLE 2: Change Magnitude Distribution")
    print("-" * 100)
    
    change_cols = ['Region', 'pct_large_pos', 'pct_large_neg', 'pct_small_change']
    change_df = df[change_cols].copy()
    change_df = change_df.sort_values('pct_large_pos', ascending=False)
    
    print(change_df.to_string(index=False, formatters={
        'pct_large_pos': lambda x: f"{x:.1f}%",
        'pct_large_neg': lambda x: f"{x:.1f}%",
        'pct_small_change': lambda x: f"{x:.1f}%"
    }))

# =============================================================================
# MAIN ANALYSIS
# =============================================================================
def main():
    print("="*100)
    print("DIAGNOSTIC ANALYSIS FOR FIGURE 3: PRE vs POST BREAKPOINT (WUE$_T$)")
    print("="*100)
    print("\nAnalyzing changes in WUE_T decline frequency")
    print("before and after breakpoint...\n")
    
    # Print caveats first
    print_statistical_caveats()
    
    # Load data
    pre_data, post_data, lat_grid, lon_grid = load_and_align_data(PRE_FILE, POST_FILE)
    
    # Calculate delta
    delta_data, mask_both = calculate_delta(pre_data, post_data)
    
    # Add grids to region definitions
    for region_def in COASTAL_REGIONS.values():
        region_def["_lon_grid"] = lon_grid
        region_def["_lat_grid"] = lat_grid
    
    # Analyze each region
    print("\nAnalyzing coastal regions...")
    results = []
    for region_name, region_def in COASTAL_REGIONS.items():
        region_stats = analyze_region(region_name, region_def, pre_data, post_data, delta_data, mask_both)
        if region_stats is not None:
            results.append(region_stats)
            print(f"  ✓ {region_name}: {results[-1]['n_pixels']:,} pixels, Δ = {results[-1]['Δ_mean']:+.3f}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Generate reports
    print_overall_summary(df)
    print_regional_comparison(df)
    print_detailed_statistics(df)
    print_nature_style_narrative(df)
    
    # Save results
    output_dir = os.path.dirname(PRE_FILE)
    output_csv = os.path.join(output_dir, "FIGURE3_PRE_POST_DELTA_DIAGNOSTICS_WUE_T.csv")
    df.to_csv(output_csv, index=False)
    print(f"\n\nResults saved to: {output_csv}")
    
    # Additional summary statistics
    print("\n" + "="*100)
    print("ADDITIONAL DIAGNOSTIC STATISTICS")
    print("="*100)
    
    # Overall delta statistics
    n_valid = int(np.sum(mask_both))
    valid_delta = delta_data[mask_both]
    print(f"\nOverall Δ Statistics (all valid pixels):")
    print(f"  • Pixels analyzed: {n_valid:,}")
    print(f"  • Mean Δ: {np.nanmean(valid_delta):+.4f}")
    print(f"  • Median Δ: {np.nanmedian(valid_delta):+.4f}")
    print(f"  • Std Δ: {np.nanstd(valid_delta):.4f}")
    print(f"  • IQR: {np.nanpercentile(valid_delta, 75) - np.nanpercentile(valid_delta, 25):.4f}")
    print(f"  • Range: [{np.nanmin(valid_delta):+.4f}, {np.nanmax(valid_delta):+.4f}]")
    
    # Direction summary
    pos_changes = np.sum(valid_delta > 0)
    neg_changes = np.sum(valid_delta < 0)
    no_changes = np.sum(valid_delta == 0)
    
    print(f"\nChange Direction Distribution:")
    print(f"  • Increased (Δ > 0): {pos_changes:,} pixels ({100*pos_changes/len(valid_delta):.1f}%)")
    print(f"  • Decreased (Δ < 0): {neg_changes:,} pixels ({100*neg_changes/len(valid_delta):.1f}%)")
    print(f"  • No change (Δ = 0): {no_changes:,} pixels ({100*no_changes/len(valid_delta):.1f}%)")
    
    # Extreme changes
    extreme_pos = np.sum(valid_delta >= 0.2)
    extreme_neg = np.sum(valid_delta <= -0.2)
    
    print(f"\nExtreme Changes:")
    print(f"  • Extreme increases (Δ ≥ 0.20): {extreme_pos:,} pixels ({100*extreme_pos/len(valid_delta):.1f}%)")
    print(f"  • Extreme decreases (Δ ≤ -0.20): {extreme_neg:,} pixels ({100*extreme_neg/len(valid_delta):.1f}%)")
    
    # Note about spatial autocorrelation
    print(f"\nNOTE: Spatial autocorrelation exists among these {len(valid_delta):,} pixels.")
    print("      Statistical tests should be interpreted with caution.")
    print(f"      Wilcoxon tests performed only for regions with ≥{MIN_PIXELS_FOR_TESTS} pixels.")
    
    return df

if __name__ == "__main__":
    df = main()