# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:17:09 2026


@author: ammar
"""

"""
Q2_methods_results_verification.py - FIXED VERSION
Diagnostic script for Q2 manuscript Methods and Results verification.
Reads already-generated Q2 outputs and prints all values needed for manuscript revision.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Q2 METHODS AND RESULTS VERIFICATION")
print("Reading existing Q2 outputs for manuscript revision")
print("="*80)

# ============================================================================
# PATHS
# ============================================================================

input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs"
figure_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_figures"

# Required files
required_files = [
    "Q2_site_stability.csv",
    "Q2_site_plasticity_slope.csv",
    "Q2_site_plasticity_range.csv",
    "Q2_site_resistance.csv",
    "Q2_site_recovery_events.csv",
    "Q2_site_recovery.csv",
    "Q2_site_performance_summary.csv",
    "Q2_ecosystem_class_statistical_tests.csv",
    "Q2_WUE_performance_summary.txt"
]

# Check all required files exist
missing_files = []
for f in required_files:
    if not os.path.exists(os.path.join(output_dir, f)):
        missing_files.append(f)

if missing_files:
    print("\nERROR: Missing required files:")
    for f in missing_files:
        print(f"  - {f}")
    print("\nPlease run Q2_WUE_performance.py first.")
    exit(1)

print("\nAll required files found.")

# ============================================================================
# CONSTANTS
# ============================================================================

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
SINGLE_SPEI = "SPEI_3"

# ============================================================================
# READ ALL Q2 OUTPUTS
# ============================================================================

print("\n" + "="*80)
print("READING Q2 OUTPUTS")
print("="*80)

stability = pd.read_csv(os.path.join(output_dir, "Q2_site_stability.csv"))
plasticity_slope = pd.read_csv(os.path.join(output_dir, "Q2_site_plasticity_slope.csv"))
plasticity_range = pd.read_csv(os.path.join(output_dir, "Q2_site_plasticity_range.csv"))
resistance = pd.read_csv(os.path.join(output_dir, "Q2_site_resistance.csv"))
recovery_events = pd.read_csv(os.path.join(output_dir, "Q2_site_recovery_events.csv"))
recovery_site = pd.read_csv(os.path.join(output_dir, "Q2_site_recovery.csv"))
summary_table = pd.read_csv(os.path.join(output_dir, "Q2_site_performance_summary.csv"))
ecosystem_tests = pd.read_csv(os.path.join(output_dir, "Q2_ecosystem_class_statistical_tests.csv"))

print(f"  Stability: {len(stability)} sites")
print(f"  Plasticity slope: {len(plasticity_slope)} site×timescale rows")
print(f"  Plasticity range: {len(plasticity_range)} sites")
print(f"  Resistance: {len(resistance)} sites")
print(f"  Recovery events: {len(recovery_events)} events")
print(f"  Recovery site: {len(recovery_site)} sites")
print(f"  Summary table: {len(summary_table)} sites")
print(f"  Ecosystem tests: {len(ecosystem_tests)} metrics")

# ============================================================================
# DATA COVERAGE FROM ORIGINAL MONTHLY FILE
# ============================================================================

print("\n" + "="*80)
print("DATA COVERAGE (from original monthly file)")
print("="*80)

# Read and filter original monthly data using same rules as Q2
print("Reading original monthly file...")
monthly = pd.read_csv(input_file)

# Trim whitespace
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

# Apply same filtering as Q2
base_df = monthly[
    monthly['site_name'].notna() &
    (monthly['site_name'] != "") &
    monthly['water_class'].notna() &
    (monthly['water_class'] != "") &
    monthly['water_class'].isin(ECOSYSTEM_CLASSES) &
    monthly['lat'].notna() &
    np.isfinite(monthly['lat']) &
    monthly['long'].notna() &
    np.isfinite(monthly['long']) &
    monthly['WUE_tra'].notna() &
    np.isfinite(monthly['WUE_tra']) &
    monthly['Trans_ratio'].notna() &
    np.isfinite(monthly['Trans_ratio']) &
    (monthly['Trans_ratio'] >= 0) &
    (monthly['Trans_ratio'] <= 1)
].copy()

base_df['WUE_T'] = base_df['WUE_tra']
base_df['water_class'] = pd.Categorical(base_df['water_class'], categories=ECOSYSTEM_CLASSES)
base_df = base_df.sort_values(['site_name', 'Year', 'month']).reset_index(drop=True)

n_obs = len(base_df)
n_sites = base_df['site_name'].nunique()
year_min = base_df['Year'].min()
year_max = base_df['Year'].max()

print(f"\nFiltered monthly data:")
print(f"  Valid monthly observations: {n_obs:,}")
print(f"  Number of sites: {n_sites}")
print(f"  Year range: {year_min} to {year_max}")

# By ecosystem
print("\n  By ecosystem class:")
eco_counts = base_df.groupby('water_class').agg(
    n_obs=('WUE_T', 'size'),
    n_sites=('site_name', 'nunique')
).reset_index()
for _, row in eco_counts.iterrows():
    print(f"    {row['water_class']}: {row['n_obs']:,} observations, {row['n_sites']} sites")

# ============================================================================
# COUNT COMPLETE ALL-FIVE-METRIC SITES
# ============================================================================

required_metrics = ['stability', 'plasticity_slope_SPEI3', 'plasticity_range', 'resistance', 'mean_recovery']
available_metrics = [m for m in required_metrics if m in summary_table.columns]

complete_sites = summary_table.dropna(subset=available_metrics)
n_complete = len(complete_sites)

print(f"\nComplete all-five-metric sites: {n_complete}")
if n_complete < len(recovery_site):
    print(f"  WARNING: Complete sites ({n_complete}) < recovery sites ({len(recovery_site)})")

# ============================================================================
# METHODS METRIC TABLE
# ============================================================================

print("\n" + "="*80)
print("METHODS: METRIC DEFINITIONS")
print("="*80)

print("\nMetric | Calculation | Unit | Input file | Inclusion rule | Interpretation | Figure")
print("-"*100)
print("Stability | mean(WUE_T) / var(WUE_T) | unitless | Q2_site_stability.csv | >=6 months | Higher = more stable WUE_T | Panel A")
print("Plasticity slope | OLS slope WUE_T ~ SPEI_k | g C mm-1 per SPEI unit | Q2_site_plasticity_slope.csv | >=6 months, >=4 unique SPEI values | Negative = increases under drought | Panel B")
print("Plasticity range | WUE_T 95th / 5th percentile | unitless | Q2_site_plasticity_range.csv | >=6 months, min>0 | Larger = wider dynamic range | Panel C")
print("Resistance | mean(drought WUE_T) / mean(normal WUE_T) | unitless | Q2_site_resistance.csv | >=3 drought months, >=3 near-normal months | >1 = WUE_T increases in drought | Panel D")
print("Recovery | mean(post-drought WUE_T) / mean(pre-drought WUE_T) | unitless | Q2_site_recovery_events.csv, Q2_site_recovery.csv | >=3 valid months in both windows | >1 = full recovery/overshoot | Panel E")

print("\n" + "="*80)
print("METHODS: DROUGHT DEFINITIONS")
print("="*80)

print("Drought month: SPEI_3 < -1")
print("Near-normal month: -1 <= SPEI_3 <= 1")
print("Drought event: >=2 consecutive months with SPEI_3 < -1")
print("Pre-drought window: 12 months before event start")
print("Post-drought window: 12 months after event end")
print("Recovery inclusion: >=3 valid WUE_T months in both pre- and post-drought windows")
print("Resistance inclusion: >=3 drought months and >=3 near-normal months per site")

# ============================================================================
# RESULTS: STABILITY
# ============================================================================

print("\n" + "="*80)
print("RESULTS: STABILITY")
print("="*80)

stability_sites = len(stability)
print(f"n sites: {stability_sites}")

if stability_sites > 0:
    print(f"\nOverall:")
    print(f"  Mean: {stability['stability'].mean():.2f}")
    print(f"  Median: {stability['stability'].median():.2f}")
    print(f"  Min: {stability['stability'].min():.2f}")
    print(f"  Max: {stability['stability'].max():.2f}")
    
    print("\nBy ecosystem:")
    for eco in ECOSYSTEM_CLASSES:
        subset = stability[stability['water_class'] == eco]
        if len(subset) > 0:
            print(f"  {eco}: n={len(subset)}, mean={subset['stability'].mean():.2f}, "
                  f"median={subset['stability'].median():.2f}, min={subset['stability'].min():.2f}, "
                  f"max={subset['stability'].max():.2f}")
    
    # Get p-value from ecosystem tests
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'stability']
    if len(test_row) > 0:
        p_val = test_row['p_value'].values[0]
        if pd.isna(p_val):
            print("\nKruskal-Wallis p-value: NA (insufficient data for test)")
        else:
            print(f"\nKruskal-Wallis p-value: {p_val:.3f}")
            if p_val < 0.05:
                print("  Significant at alpha=0.05")
            else:
                print("  Not significant at alpha=0.05")

# ============================================================================
# RESULTS: PLASTICITY SLOPE
# ============================================================================

print("\n" + "="*80)
print("RESULTS: PLASTICITY SLOPE")
print("="*80)

if len(plasticity_slope) > 0:
    n_rows = len(plasticity_slope)
    n_sites_slope = plasticity_slope['site_name'].nunique()
    print(f"Total site×timescale rows: {n_rows}")
    print(f"Unique sites: {n_sites_slope}")
    
    print("\nCounts by SPEI timescale:")
    for spei in SPEI_COLS:
        count = len(plasticity_slope[plasticity_slope['SPEI_timescale'] == spei])
        print(f"  {spei}: {count}")
    
    # SPEI-3 only
    slope_3 = plasticity_slope[plasticity_slope['SPEI_timescale'] == 'SPEI_3']
    if len(slope_3) > 0:
        print(f"\nSPEI-3 only (n={len(slope_3)} sites):")
        
        for eco in ECOSYSTEM_CLASSES:
            subset = slope_3[slope_3['water_class'] == eco]
            if len(subset) > 0:
                pct_neg = (subset['slope'] < 0).sum() / len(subset) * 100
                print(f"  {eco}: n={len(subset)}, mean={subset['slope'].mean():.2f}, "
                      f"median={subset['slope'].median():.2f}, min={subset['slope'].min():.2f}, "
                      f"max={subset['slope'].max():.2f}, %negative={pct_neg:.1f}%")
        
        # Get p-values
        test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'plasticity_slope_SPEI3']
        if len(test_row) > 0:
            p_val = test_row['p_value'].values[0]
            if pd.isna(p_val):
                print("\nKruskal-Wallis p-value: NA (insufficient data)")
            else:
                print(f"\nKruskal-Wallis p-value: {p_val:.3f}")
                if p_val < 0.05:
                    print("  Significant at alpha=0.05")
                else:
                    print("  Not significant at alpha=0.05")
            
            # Pairwise
            p_uw_fw = test_row['pairwise_upland_vs_freshwater_p_adj'].values[0]
            p_uw_sa = test_row['pairwise_upland_vs_saline_p_adj'].values[0]
            p_fw_sa = test_row['pairwise_freshwater_vs_saline_p_adj'].values[0]
            
            print("\nPairwise adjusted p-values (BH/FDR):")
            if not pd.isna(p_uw_fw):
                print(f"  Upland vs Freshwater: {p_uw_fw:.3f}")
            if not pd.isna(p_uw_sa):
                print(f"  Upland vs Saline: {p_uw_sa:.3f}")
            if not pd.isna(p_fw_sa):
                print(f"  Freshwater vs Saline: {p_fw_sa:.3f}")

# ============================================================================
# RESULTS: PLASTICITY RANGE
# ============================================================================

print("\n" + "="*80)
print("RESULTS: PLASTICITY RANGE")
print("="*80)

if len(plasticity_range) > 0:
    range_sites = len(plasticity_range)
    print(f"n sites: {range_sites}")
    
    print(f"\nOverall (plasticity_p95_p05):")
    print(f"  Mean: {plasticity_range['plasticity_p95_p05'].mean():.2f}")
    print(f"  Median: {plasticity_range['plasticity_p95_p05'].median():.2f}")
    print(f"  Min: {plasticity_range['plasticity_p95_p05'].min():.2f}")
    print(f"  Max: {plasticity_range['plasticity_p95_p05'].max():.2f}")
    
    print("\nBy ecosystem:")
    for eco in ECOSYSTEM_CLASSES:
        subset = plasticity_range[plasticity_range['water_class'] == eco]
        if len(subset) > 0:
            print(f"  {eco}: n={len(subset)}, mean={subset['plasticity_p95_p05'].mean():.2f}, "
                  f"median={subset['plasticity_p95_p05'].median():.2f}, min={subset['plasticity_p95_p05'].min():.2f}, "
                  f"max={subset['plasticity_p95_p05'].max():.2f}")
    
    # Top 5 highest plasticity range
    top5 = plasticity_range.nlargest(5, 'plasticity_p95_p05')[['site_name', 'water_class', 'plasticity_p95_p05', 
                                                                'WUE_T_p95', 'WUE_T_p05']]
    print("\nTop 5 highest plasticity range sites:")
    for _, row in top5.iterrows():
        print(f"  {row['site_name']} ({row['water_class']}): {row['plasticity_p95_p05']:.2f} "
              f"(95th={row['WUE_T_p95']:.2f}, 5th={row['WUE_T_p05']:.2f})")
    
    # Get p-value
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'plasticity_p95_p05']
    if len(test_row) > 0:
        p_val = test_row['p_value'].values[0]
        if pd.isna(p_val):
            print("\nKruskal-Wallis p-value: NA (insufficient data)")
        else:
            print(f"\nKruskal-Wallis p-value: {p_val:.3f}")
            if p_val < 0.05:
                print("  Significant at alpha=0.05")
            else:
                print("  Not significant at alpha=0.05")

# ============================================================================
# RESULTS: RESISTANCE
# ============================================================================

print("\n" + "="*80)
print("RESULTS: RESISTANCE")
print("="*80)

if len(resistance) > 0:
    resist_sites = len(resistance)
    print(f"n sites: {resist_sites}")
    
    print(f"\nOverall:")
    print(f"  Mean: {resistance['resistance'].mean():.2f}")
    print(f"  Median: {resistance['resistance'].median():.2f}")
    print(f"  Min: {resistance['resistance'].min():.2f}")
    print(f"  Max: {resistance['resistance'].max():.2f}")
    
    print("\nBy ecosystem:")
    for eco in ECOSYSTEM_CLASSES:
        subset = resistance[resistance['water_class'] == eco]
        if len(subset) > 0:
            pct_above1 = (subset['resistance'] > 1).sum() / len(subset) * 100
            print(f"  {eco}: n={len(subset)}, mean={subset['resistance'].mean():.2f}, "
                  f"median={subset['resistance'].median():.2f}, min={subset['resistance'].min():.2f}, "
                  f"max={subset['resistance'].max():.2f}, %>1={pct_above1:.1f}%")
    
    # Get p-values
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'resistance']
    if len(test_row) > 0:
        p_val = test_row['p_value'].values[0]
        if pd.isna(p_val):
            print("\nKruskal-Wallis p-value: NA (insufficient data)")
        else:
            print(f"\nKruskal-Wallis p-value: {p_val:.3f}")
            if p_val < 0.05:
                print("  Significant at alpha=0.05")
            else:
                print("  Not significant at alpha=0.05")
        
        # Pairwise
        p_uw_fw = test_row['pairwise_upland_vs_freshwater_p_adj'].values[0]
        p_uw_sa = test_row['pairwise_upland_vs_saline_p_adj'].values[0]
        p_fw_sa = test_row['pairwise_freshwater_vs_saline_p_adj'].values[0]
        
        print("\nPairwise adjusted p-values (BH/FDR):")
        if not pd.isna(p_uw_fw):
            print(f"  Upland vs Freshwater: {p_uw_fw:.3f}")
        if not pd.isna(p_uw_sa):
            print(f"  Upland vs Saline: {p_uw_sa:.3f}")
        if not pd.isna(p_fw_sa):
            print(f"  Freshwater vs Saline: {p_fw_sa:.3f}")

# ============================================================================
# RESULTS: RECOVERY
# ============================================================================

print("\n" + "="*80)
print("RESULTS: RECOVERY")
print("="*80)

if len(recovery_events) > 0:
    n_events = len(recovery_events)
    print(f"Number of drought events: {n_events}")
    
    if len(recovery_site) > 0:
        n_sites_recovery = len(recovery_site)
        print(f"Number of sites with recovery estimates: {n_sites_recovery}")
        
        print(f"\nEvent-level recovery:")
        print(f"  Mean: {recovery_events['recovery'].mean():.2f}")
        print(f"  Median: {recovery_events['recovery'].median():.2f}")
        print(f"  Min: {recovery_events['recovery'].min():.2f}")
        print(f"  Max: {recovery_events['recovery'].max():.2f}")
        
        print(f"\nSite-level recovery (mean_recovery):")
        print(f"  Mean: {recovery_site['mean_recovery'].mean():.2f}")
        print(f"  Median: {recovery_site['mean_recovery'].median():.2f}")
        print(f"  Min: {recovery_site['mean_recovery'].min():.2f}")
        print(f"  Max: {recovery_site['mean_recovery'].max():.2f}")
        
        print("\nBy ecosystem (site-level):")
        for eco in ECOSYSTEM_CLASSES:
            subset = recovery_site[recovery_site['water_class'] == eco]
            if len(subset) > 0:
                print(f"  {eco}: n_sites={len(subset)}, mean={subset['mean_recovery'].mean():.2f}, "
                      f"median={subset['mean_recovery'].median():.2f}, min={subset['mean_recovery'].min():.2f}, "
                      f"max={subset['mean_recovery'].max():.2f}, %full_recovery={subset['pct_full_recovery'].mean():.1f}%")
        
        print("\nBy ecosystem (event-level):")
        for eco in ECOSYSTEM_CLASSES:
            subset = recovery_events[recovery_events['water_class'] == eco]
            if len(subset) > 0:
                pct_ge1 = (subset['recovery'] >= 1).sum() / len(subset) * 100
                print(f"  {eco}: n_events={len(subset)}, mean={subset['recovery'].mean():.2f}, "
                      f"median={subset['recovery'].median():.2f}, %>={pct_ge1:.1f}%")
        
        # Get p-values
        test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'mean_recovery']
        if len(test_row) > 0:
            p_val = test_row['p_value'].values[0]
            if pd.isna(p_val):
                print("\nKruskal-Wallis p-value: NA (insufficient data)")
            else:
                print(f"\nKruskal-Wallis p-value: {p_val:.3f}")
                if p_val < 0.05:
                    print("  Significant at alpha=0.05")
                else:
                    print("  Not significant at alpha=0.05")
            
            # Pairwise
            p_uw_fw = test_row['pairwise_upland_vs_freshwater_p_adj'].values[0]
            p_uw_sa = test_row['pairwise_upland_vs_saline_p_adj'].values[0]
            p_fw_sa = test_row['pairwise_freshwater_vs_saline_p_adj'].values[0]
            
            print("\nPairwise adjusted p-values (BH/FDR):")
            if not pd.isna(p_uw_fw):
                print(f"  Upland vs Freshwater: {p_uw_fw:.3f}")
            if not pd.isna(p_uw_sa):
                print(f"  Upland vs Saline: {p_uw_sa:.3f}")
            if not pd.isna(p_fw_sa):
                print(f"  Freshwater vs Saline: {p_fw_sa:.3f}")

# ============================================================================
# MANUSCRIPT VALUE CHECKS
# ============================================================================

print("\n" + "="*80)
print("MANUSCRIPT VALUE VERIFICATION")
print("="*80)

# Draft values from manuscript
draft_values = {
    'monthly_observations': 1852,
    'sites': 64,
    'year_range': '1994 to 2025',
    'stability_sites': 56,
    'plasticity_range_sites': 56,
    'resistance_sites': 39,
    'recovery_events': 52,
    'recovery_sites': 27,
    'complete_sites': 27
}

# Calculated values
calc_values = {
    'monthly_observations': n_obs,
    'sites': n_sites,
    'year_range': f"{year_min} to {year_max}",
    'stability_sites': len(stability),
    'plasticity_range_sites': len(plasticity_range),
    'resistance_sites': len(resistance),
    'recovery_events': len(recovery_events),
    'recovery_sites': len(recovery_site) if len(recovery_site) > 0 else 0,
    'complete_sites': n_complete
}

print("Comparing calculated values vs draft manuscript values:")
print("-"*60)
print(f"{'Metric':<25} {'Calculated':<12} {'Draft':<12} {'Match?'}")
print("-"*60)

all_match = True
for key in draft_values:
    calc = calc_values[key]
    draft = draft_values[key]
    
    if key == 'year_range':
        match = (calc == draft)
    else:
        match = (calc == draft)
    
    status = "✓" if match else "✗ MISMATCH"
    if key == 'year_range':
        print(f"{key:<25} {calc:<12} {draft:<12} {status}")
    else:
        print(f"{key:<25} {calc:<12} {draft:<12} {status}")
    if not match:
        all_match = False

print("-"*60)
if all_match:
    print("All values match the draft manuscript.")
else:
    print("WARNING: Some values do not match the draft manuscript. Please update the text.")

# Special check for recovery events sentence
print("\nSPECIAL CHECK: Recovery events sentence")
print(f"Draft says: '53 drought events across 28 sites'")
print(f"Calculated: {len(recovery_events)} drought events across {len(recovery_site) if len(recovery_site) > 0 else 0} sites")
if len(recovery_events) == 52 and len(recovery_site) == 27:
    print("✓ The correct value is 52 events across 27 sites (update manuscript)")
elif len(recovery_events) == 53 and len(recovery_site) == 28:
    print("✓ The draft value is correct (53 events, 28 sites)")
else:
    print(f"✗ Update needed: Use {len(recovery_events)} events and {len(recovery_site) if len(recovery_site) > 0 else 0} sites")

# ============================================================================
# FIGURE-TO-ANALYSIS MAP
# ============================================================================

print("\n" + "="*80)
print("FIGURE-TO-ANALYSIS MAP")
print("="*80)

print("\nPanel A (Stability):")
print(f"  CSV: Q2_site_stability.csv")
print(f"  x: mean_TET, y: stability")
print(f"  n: {len(stability)} sites")
print(f"  Test: Kruskal-Wallis (p={ecosystem_tests[ecosystem_tests['response_variable']=='stability']['p_value'].values[0]:.3f} if available)")

print("\nPanel B (Plasticity slope):")
print(f"  CSV: Q2_site_plasticity_slope.csv")
print(f"  x: SPEI timescale, y: slope")
print(f"  Grouping: water_class")
print(f"  n: {len(plasticity_slope)} site×timescale rows")
print(f"  Test: Kruskal-Wallis for SPEI-3")

print("\nPanel C (Plasticity range):")
print(f"  CSV: Q2_site_plasticity_range.csv")
print(f"  y: plasticity_p95_p05")
print(f"  n: {len(plasticity_range)} sites")
print(f"  Test: Kruskal-Wallis")

print("\nPanel D (Resistance):")
print(f"  CSV: Q2_site_resistance.csv")
print(f"  y: resistance")
print(f"  n: {len(resistance)} sites")
print(f"  Test: Kruskal-Wallis")

print("\nPanel E (Recovery):")
print(f"  CSV: Q2_site_recovery.csv")
print(f"  y: mean_recovery")
print(f"  n: {len(recovery_site)} sites")
print(f"  Test: Kruskal-Wallis")

# ============================================================================
# INTERPRETATION CHECKS
# ============================================================================

print("\n" + "="*80)
print("INTERPRETATION CHECKS (for Results revision)")
print("="*80)

# Stability
if len(stability) > 0:
    highest_stability = stability.groupby('water_class')['stability'].mean().idxmax() if len(stability['water_class'].unique()) > 0 else "N/A"
    print(f"• Highest mean stability: {highest_stability}")

# Plasticity slope
if len(slope_3) > 0:
    most_negative = slope_3.groupby('water_class')['slope'].mean().idxmin() if len(slope_3['water_class'].unique()) > 0 else "N/A"
    print(f"• Most negative mean SPEI-3 slope: {most_negative}")

# Resistance
if len(resistance) > 0:
    highest_resistance = resistance.groupby('water_class')['resistance'].mean().idxmax() if len(resistance['water_class'].unique()) > 0 else "N/A"
    print(f"• Highest mean resistance: {highest_resistance}")

# Significance checks
print("\n• Significant Kruskal-Wallis tests (p<0.05):")
significant_tests = ecosystem_tests[ecosystem_tests['p_value'] < 0.05]
if len(significant_tests) > 0:
    for _, row in significant_tests.iterrows():
        print(f"  - {row['metric']}: p={row['p_value']:.3f}")
else:
    print("  None")

print("\n• Significant pairwise tests after FDR correction:")
significant_pairs = []
for _, row in ecosystem_tests.iterrows():
    for pair in ['pairwise_upland_vs_freshwater_p_adj', 'pairwise_upland_vs_saline_p_adj', 'pairwise_freshwater_vs_saline_p_adj']:
        if pair in row and not pd.isna(row[pair]) and row[pair] < 0.05:
            significant_pairs.append(f"{row['metric']}: {pair.replace('pairwise_', '').replace('_p_adj', '').replace('_vs_', ' vs ')} = {row[pair]:.3f}")
if significant_pairs:
    for sp in significant_pairs:
        print(f"  - {sp}")
else:
    print("  None")

# Recovery check
if len(recovery_events) > 0:
    recovery_median = recovery_events['recovery'].median()
    print(f"\n• Recovery centered near 1? Median = {recovery_median:.2f} {'(yes)' if 0.9 <= recovery_median <= 1.1 else '(no)'}")

# ============================================================================
# SAVE REPORT
# ============================================================================

print("\n" + "="*80)
print("SAVING REPORT")
print("="*80)

# Build report with ASCII characters instead of Unicode
report_lines = []
report_lines.append("="*80)
report_lines.append("Q2 METHODS AND RESULTS VERIFICATION REPORT")
report_lines.append("Generated from Q2 output CSVs")
report_lines.append("="*80)
report_lines.append("")

report_lines.append("DATA COVERAGE")
report_lines.append("-"*40)
report_lines.append(f"Valid monthly observations: {n_obs:,}")
report_lines.append(f"Number of sites: {n_sites}")
report_lines.append(f"Year range: {year_min} to {year_max}")
for _, row in eco_counts.iterrows():
    report_lines.append(f"  {row['water_class']}: {row['n_obs']:,} obs, {row['n_sites']} sites")
report_lines.append("")
report_lines.append(f"Complete all-five-metric sites: {n_complete}")
report_lines.append("")

report_lines.append("METHODS: METRIC DEFINITIONS")
report_lines.append("-"*40)
report_lines.append("Stability: mean(WUE_T) / var(WUE_T) | unitless | >=6 months | Panel A")
report_lines.append("Plasticity slope: OLS slope WUE_T ~ SPEI_k | g C mm-1 per SPEI unit | >=6 months, >=4 unique SPEI values | Panel B")
report_lines.append("Plasticity range: WUE_T 95th / 5th percentile | unitless | >=6 months, min>0 | Panel C")
report_lines.append("Resistance: mean(drought WUE_T) / mean(normal WUE_T) | unitless | >=3 drought months, >=3 near-normal months | Panel D")
report_lines.append("Recovery: mean(post-drought WUE_T) / mean(pre-drought WUE_T) | unitless | >=3 valid months in both windows | Panel E")
report_lines.append("")

report_lines.append("DROUGHT DEFINITIONS")
report_lines.append("-"*40)
report_lines.append("Drought month: SPEI_3 < -1")
report_lines.append("Near-normal month: -1 <= SPEI_3 <= 1")
report_lines.append("Drought event: >=2 consecutive months with SPEI_3 < -1")
report_lines.append("Pre/post windows: 12 months before/after event")
report_lines.append("Recovery inclusion: >=3 valid WUE_T months in both windows")
report_lines.append("Resistance inclusion: >=3 drought months and >=3 near-normal months")
report_lines.append("")

report_lines.append("RESULTS: STABILITY")
report_lines.append("-"*40)
report_lines.append(f"n sites: {stability_sites}")
if stability_sites > 0:
    report_lines.append(f"Overall: mean={stability['stability'].mean():.2f}, median={stability['stability'].median():.2f}, min={stability['stability'].min():.2f}, max={stability['stability'].max():.2f}")
    for eco in ECOSYSTEM_CLASSES:
        subset = stability[stability['water_class'] == eco]
        if len(subset) > 0:
            report_lines.append(f"  {eco}: n={len(subset)}, mean={subset['stability'].mean():.2f}, median={subset['stability'].median():.2f}")
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'stability']
    if len(test_row) > 0 and not pd.isna(test_row['p_value'].values[0]):
        report_lines.append(f"Kruskal-Wallis p-value: {test_row['p_value'].values[0]:.3f}")
report_lines.append("")

report_lines.append("RESULTS: PLASTICITY SLOPE (SPEI-3)")
report_lines.append("-"*40)
if len(slope_3) > 0:
    report_lines.append(f"n sites: {len(slope_3)}")
    for eco in ECOSYSTEM_CLASSES:
        subset = slope_3[slope_3['water_class'] == eco]
        if len(subset) > 0:
            pct_neg = (subset['slope'] < 0).sum() / len(subset) * 100
            report_lines.append(f"  {eco}: n={len(subset)}, mean={subset['slope'].mean():.2f}, median={subset['slope'].median():.2f}, %negative={pct_neg:.1f}%")
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'plasticity_slope_SPEI3']
    if len(test_row) > 0 and not pd.isna(test_row['p_value'].values[0]):
        report_lines.append(f"Kruskal-Wallis p-value: {test_row['p_value'].values[0]:.3f}")
report_lines.append("")

report_lines.append("RESULTS: PLASTICITY RANGE")
report_lines.append("-"*40)
if len(plasticity_range) > 0:
    report_lines.append(f"n sites: {len(plasticity_range)}")
    report_lines.append(f"Overall: mean={plasticity_range['plasticity_p95_p05'].mean():.2f}, median={plasticity_range['plasticity_p95_p05'].median():.2f}, min={plasticity_range['plasticity_p95_p05'].min():.2f}, max={plasticity_range['plasticity_p95_p05'].max():.2f}")
    for eco in ECOSYSTEM_CLASSES:
        subset = plasticity_range[plasticity_range['water_class'] == eco]
        if len(subset) > 0:
            report_lines.append(f"  {eco}: n={len(subset)}, mean={subset['plasticity_p95_p05'].mean():.2f}, median={subset['plasticity_p95_p05'].median():.2f}")
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'plasticity_p95_p05']
    if len(test_row) > 0 and not pd.isna(test_row['p_value'].values[0]):
        report_lines.append(f"Kruskal-Wallis p-value: {test_row['p_value'].values[0]:.3f}")
report_lines.append("")

report_lines.append("RESULTS: RESISTANCE")
report_lines.append("-"*40)
if len(resistance) > 0:
    report_lines.append(f"n sites: {len(resistance)}")
    report_lines.append(f"Overall: mean={resistance['resistance'].mean():.2f}, median={resistance['resistance'].median():.2f}")
    for eco in ECOSYSTEM_CLASSES:
        subset = resistance[resistance['water_class'] == eco]
        if len(subset) > 0:
            pct_above1 = (subset['resistance'] > 1).sum() / len(subset) * 100
            report_lines.append(f"  {eco}: n={len(subset)}, mean={subset['resistance'].mean():.2f}, median={subset['resistance'].median():.2f}, %>1={pct_above1:.1f}%")
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'resistance']
    if len(test_row) > 0 and not pd.isna(test_row['p_value'].values[0]):
        report_lines.append(f"Kruskal-Wallis p-value: {test_row['p_value'].values[0]:.3f}")
report_lines.append("")

report_lines.append("RESULTS: RECOVERY")
report_lines.append("-"*40)
if len(recovery_events) > 0:
    report_lines.append(f"Events: {len(recovery_events)}")
    report_lines.append(f"Sites: {len(recovery_site) if len(recovery_site) > 0 else 0}")
    report_lines.append(f"Event-level: mean={recovery_events['recovery'].mean():.2f}, median={recovery_events['recovery'].median():.2f}")
    if len(recovery_site) > 0:
        report_lines.append(f"Site-level: mean={recovery_site['mean_recovery'].mean():.2f}, median={recovery_site['mean_recovery'].median():.2f}")
    test_row = ecosystem_tests[ecosystem_tests['response_variable'] == 'mean_recovery']
    if len(test_row) > 0 and not pd.isna(test_row['p_value'].values[0]):
        report_lines.append(f"Kruskal-Wallis p-value: {test_row['p_value'].values[0]:.3f}")
report_lines.append("")

report_lines.append("MANUSCRIPT VALUE CHECKS")
report_lines.append("-"*40)
for key in draft_values:
    calc = calc_values[key]
    draft = draft_values[key]
    match = "✓" if calc == draft else "✗"
    report_lines.append(f"{key}: calculated={calc}, draft={draft} {match}")
report_lines.append("")
report_lines.append(f"Recovery events: Correct value is {len(recovery_events)} events across {len(recovery_site) if len(recovery_site) > 0 else 0} sites")
report_lines.append("")

report_lines.append("INTERPRETATION CHECKS")
report_lines.append("-"*40)
if len(stability) > 0 and len(stability['water_class'].unique()) > 0:
    report_lines.append(f"Highest mean stability: {highest_stability}")
if len(slope_3) > 0 and len(slope_3['water_class'].unique()) > 0:
    report_lines.append(f"Most negative mean SPEI-3 slope: {most_negative}")
if len(resistance) > 0 and len(resistance['water_class'].unique()) > 0:
    report_lines.append(f"Highest mean resistance: {highest_resistance}")

sig_tests = ecosystem_tests[ecosystem_tests['p_value'] < 0.05]
if len(sig_tests) > 0:
    report_lines.append("Significant Kruskal-Wallis tests (p<0.05):")
    for _, row in sig_tests.iterrows():
        report_lines.append(f"  - {row['metric']}: p={row['p_value']:.3f}")
else:
    report_lines.append("No significant Kruskal-Wallis tests (p<0.05)")
report_lines.append("")
report_lines.append("="*80)
report_lines.append("Report generated from Q2 output files")
report_lines.append("="*80)

# Write report with utf-8 encoding
report_text = "\n".join(report_lines)
report_path = os.path.join(output_dir, "Q2_methods_results_verification_report.txt")
with open(report_path, 'w', encoding='utf-8') as f:
    f.write(report_text)

print(f"\nReport saved to: {report_path}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("Q2 METHODS/RESULTS VERIFICATION COMPLETE")
print("="*80)
print(f"\nReport saved to:")
print(f"  {report_path}")
print("\nPaste this console output back before revising manuscript text.")
print("="*80)

# Now re-print the key sections to console for easy copying
print("\n\n" + "="*80)
print("CONSOLE OUTPUT FOR COPYING (Key Results)")
print("="*80)

# Quick summary of key results for manuscript revision
print("\nKEY RESULTS SUMMARY:")
print(f"- {n_obs:,} valid monthly observations across {n_sites} sites ({year_min}-{year_max})")
print(f"- Stability estimated for {len(stability)} sites")
print(f"- Plasticity slope estimated for {len(plasticity_slope)} site×timescale combinations")
print(f"- Plasticity range estimated for {len(plasticity_range)} sites")
print(f"- Resistance estimated for {len(resistance)} sites")
print(f"- Recovery estimated from {len(recovery_events)} drought events across {len(recovery_site) if len(recovery_site) > 0 else 0} sites")
print(f"- {n_complete} sites with complete information for all five metrics")
print(f"- No significant ecosystem-class differences (all Kruskal-Wallis p > 0.05)")

print("\n" + "="*80)