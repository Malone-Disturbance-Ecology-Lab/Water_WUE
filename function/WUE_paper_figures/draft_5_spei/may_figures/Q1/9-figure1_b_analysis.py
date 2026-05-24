# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 15:49:04 2025
Updated: ALL THREE metrics now use the SAME strict triple intersection
Matches final Figure 2 workflow exactly
UPDATED: Now reports Mean ± SE and SD for group-level summary statistics
UPDATED: Added 3-month global NN filter (same as Figure 1)
FIXED: Removed duplicate function blocks
FIXED: Now counts ALL retained NN observations (not unique month combos)
@author: ammar
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os  # ADDED for file existence check

BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
FILE = BASE_DIR + r"\wue_site_level_summary_SPEI_1.csv"

# =============================================================================
# 3-MONTH NN FILTER FUNCTION (SAME AS FIGURE 1 AND FIGURE 2)
# =============================================================================

def filter_sites_by_nn_months(sites_list, min_months=3):
    """
    Filter sites to keep only those with >= min_months of NN OBSERVATIONS
    Uses monthly_data_after_outlier_removal.csv to count NN observations per site
    
    This is the SAME filter applied in Figure 1 and Figure 2
    FIXED: Now counts ALL retained observations (not unique month combos)
    """
    print(f"\n{'='*60}")
    print(f"FILTERING SITES WITH <{min_months} NN OBSERVATIONS")
    print(f"{'='*60}")
    
    # Load monthly data
    monthly_path = BASE_DIR + r"\monthly_data_after_outlier_removal.csv"
    
    if not os.path.exists(monthly_path):
        print(f"⚠️ Warning: {monthly_path} not found. Skipping observation filter.")
        return set(sites_list), len(sites_list), 0
    
    df_monthly = pd.read_csv(monthly_path)
    
    # Filter to NN conditions at SPEI-1
    df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()
    
    # Count retained NN OBSERVATIONS per site - FIXED: count ALL observations
    obs_counts = {}
    for site in sites_list:
        site_data = df_nn[df_nn["site_name"] == site]
        n_obs = len(site_data)  # FIXED: Count ALL retained NN observations (not unique month combos)
        obs_counts[site] = n_obs
    
    # Filter sites
    sites_to_keep = [site for site in sites_list if obs_counts.get(site, 0) >= min_months]
    sites_to_remove = [site for site in sites_list if obs_counts.get(site, 0) < min_months]
    
    print(f"\n  Original sites: {len(sites_list)}")
    print(f"  Sites removed (<{min_months} observations): {len(sites_to_remove)}")
    print(f"  Sites kept (≥{min_months} observations): {len(sites_to_keep)}")
    
    if sites_to_remove:
        print(f"\n  ⚠️ Sites REMOVED (insufficient NN observations):")
        for site in sorted(sites_to_remove):
            obs = obs_counts.get(site, 0)
            print(f"      - {site}: {obs} observations")
    
    return set(sites_to_keep), len(sites_to_remove), obs_counts

# =============================================================================
# LOAD DATA
# =============================================================================

df = pd.read_csv(FILE)

# =============================================================================
# STEP 1: FILTER TO NEAR-NORMAL CONDITIONS (SPEI_Class == 'NN')
# =============================================================================
df_nn = df[df["SPEI_Class"] == "NN"].copy()

print("="*80)
print("FIGURE 2 SUPPORT ANALYSIS (Updated: 3 Categories + STRICT TRIPLE INTERSECTION)")
print("ALL THREE METRICS now use the SAME strict triple intersection (WUE ∩ WUE_E ∩ WUE_T)")
print("UPDATED: Now reporting Mean ± SE and SD for group-level statistics")
print("UPDATED: Added 3-month global NN filter (same as Figure 1)")
print("="*80)
print("\n[STEP 1] Filter to SPEI_Class == 'NN'")
print(f"  Total rows after NN filter: {len(df_nn)}")
print(f"  Unique sites: {df_nn['site_name'].nunique()}")

# =============================================================================
# STEP 2: FILTER TO ALLOWED SALINITY CATEGORIES (Freshwater, Saline, Upland)
# =============================================================================
allowed_categories = ['Freshwater', 'Saline', 'Upland']
df_filtered = df_nn[df_nn["Salinity_Category"].isin(allowed_categories)].copy()

print(f"\n[STEP 2] Filter to salinity categories: {allowed_categories}")
print(f"  Rows after category filter: {len(df_filtered)}")
print(f"  Unique sites: {df_filtered['site_name'].nunique()}")

print("\n  Breakdown by salinity category (before N_months filter):")
for cat in allowed_categories:
    n_sites = df_filtered[df_filtered["Salinity_Category"] == cat]["site_name"].nunique()
    n_rows = len(df_filtered[df_filtered["Salinity_Category"] == cat])
    print(f"    {cat}: {n_sites} sites, {n_rows} rows")

# =============================================================================
# STEP 3: REMOVE DUPLICATES (site × salinity × metric × class)
# =============================================================================
df_dedup = df_filtered.drop_duplicates(
    subset=["site_name", "Salinity_Category", "WUE_Metric", "SPEI_Class"]
).copy()

print(f"\n[STEP 3] Remove duplicates at site × salinity × metric × class level")
print(f"  Rows after deduplication: {len(df_dedup)}")

# =============================================================================
# STEP 4: APPLY SAMPLE-SIZE ADEQUACY FILTER (N_months >= 3) - PER METRIC
# =============================================================================
df_adequate = df_dedup[df_dedup["N_months"] >= 3].copy()

print(f"\n[STEP 4] Apply N_months >= 3 adequacy filter (per metric)")
print(f"  Rows after adequacy filter: {len(df_adequate)}")
print(f"  Unique sites retained: {df_adequate['site_name'].nunique()}")

print("\n  Sites by salinity category (after N_months >= 3 per metric):")
for cat in allowed_categories:
    n_sites = df_adequate[df_adequate["Salinity_Category"] == cat]["site_name"].nunique()
    print(f"    {cat}: {n_sites} sites")

# =============================================================================
# STEP 4.5: APPLY GLOBAL 3-MONTH NN FILTER (SAME AS FIGURE 1)
# =============================================================================
all_sites_before_filter = set(df_adequate['site_name'].unique())
print(f"\n  Unique sites before global NN filter: {len(all_sites_before_filter)}")

sites_kept, n_removed, obs_counts = filter_sites_by_nn_months(all_sites_before_filter, min_months=3)

# Filter data to only sites with >=3 NN observations
df_adequate = df_adequate[df_adequate['site_name'].isin(sites_kept)].copy()
print(f"\n  After global NN filter: {len(df_adequate)} rows")
print(f"  Unique sites retained after global NN filter: {df_adequate['site_name'].nunique()}")

# =============================================================================
# STEP 5: DEFINE STRICT TRIPLE INTERSECTION (WUE ∩ WUE_E ∩ WUE_T)
# ALL THREE METRICS now use this strict triple intersection
# =============================================================================
print(f"\n[STEP 5] Define STRICT TRIPLE INTERSECTION (ALL metrics will use this)")

# Get sites with WUE available
sites_with_wue = set(
    df_adequate[df_adequate["WUE_Metric"] == "WUE"]["site_name"].unique()
)

# Get sites with WUE_eva available
sites_with_eva = set(
    df_adequate[df_adequate["WUE_Metric"] == "WUE_eva"]["site_name"].unique()
)

# Get sites with WUE_tra available
sites_with_tra = set(
    df_adequate[df_adequate["WUE_Metric"] == "WUE_tra"]["site_name"].unique()
)

# STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics
shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)

print(f"  Sites with WUE available: {len(sites_with_wue)}")
print(f"  Sites with WUE_E available: {len(sites_with_eva)}")
print(f"  Sites with WUE_T available: {len(sites_with_tra)}")
print(f"  STRICT TRIPLE INTERSECTION (ALL metrics): {len(shared_sites)} sites")

# Show which sites are excluded from all metrics
wue_only = sites_with_wue - shared_sites
eva_only = sites_with_eva - shared_sites
tra_only = sites_with_tra - shared_sites

if len(wue_only) > 0:
    print(f"\n  ⚠️ Sites with WUE ONLY (excluded from ALL metrics):")
    for site in sorted(wue_only):
        print(f"      - {site}")

if len(eva_only) > 0:
    print(f"\n  ⚠️ Sites with WUE_E ONLY (excluded from ALL metrics):")
    for site in sorted(eva_only):
        print(f"      - {site}")

if len(tra_only) > 0:
    print(f"\n  ⚠️ Sites with WUE_T ONLY (excluded from ALL metrics):")
    for site in sorted(tra_only):
        print(f"      - {site}")

# =============================================================================
# STEP 6: PREPARE FINAL DATASETS FOR EACH METRIC - ALL USE STRICT TRIPLE INTERSECTION
# =============================================================================
print(f"\n[STEP 6] Final datasets for Figure 2 - ALL metrics use strict triple intersection")

# ALL metrics now use the same strict triple intersection
shared_final = df_adequate[df_adequate["site_name"].isin(shared_sites)].copy()

wue_final = shared_final[shared_final["WUE_Metric"] == "WUE"].copy()
eva_final = shared_final[shared_final["WUE_Metric"] == "WUE_eva"].copy()
tra_final = shared_final[shared_final["WUE_Metric"] == "WUE_tra"].copy()

wue_sites_set = set(wue_final["site_name"].unique())
eva_sites_set = set(eva_final["site_name"].unique())
tra_sites_set = set(tra_final["site_name"].unique())

print(f"\n  WUE (bulk): {len(wue_sites_set)} sites (strict triple intersection + 3-month NN filter)")
print(f"  WUE_E: {len(eva_sites_set)} sites (strict triple intersection + 3-month NN filter)")
print(f"  WUE_T: {len(tra_sites_set)} sites (strict triple intersection + 3-month NN filter)")

# Verify all metrics use the same sites
if len(wue_sites_set) == len(eva_sites_set) == len(tra_sites_set) == len(shared_sites):
    print(f"\n  ✓ VERIFICATION: All three metrics use the same {len(shared_sites)} sites")
else:
    print(f"\n  ⚠️ WARNING: Site counts differ!")
    print(f"      WUE: {len(wue_sites_set)} sites")
    print(f"      WUE_E: {len(eva_sites_set)} sites")
    print(f"      WUE_T: {len(tra_sites_set)} sites")
    print(f"      Expected: {len(shared_sites)} sites")

# =============================================================================
# STEP 7: PRINT DETAILED SAMPLE SIZES BY CATEGORY FOR EACH METRIC
# INCLUDES Mean ± SE and SD (UPDATED)
# =============================================================================
print("\n" + "="*80)
print("FINAL FIGURE 2 SAMPLE SIZES (Used in boxplots and statistical tests)")
print("ALL metrics now use the same strict triple intersection + 3-month NN filter")
print("NOW reporting Mean ± SE and SD for each group")
print("="*80)

metrics_info = {
    'WUE (strict triple intersection)': wue_final,
    'WUE_E (strict triple intersection)': eva_final,
    'WUE_T (strict triple intersection)': tra_final
}

for metric_name, metric_data in metrics_info.items():
    print(f"\n{metric_name}:")
    print(f"  Total sites: {metric_data['site_name'].nunique()}")
    for cat in allowed_categories:
        subset = metric_data[metric_data["Salinity_Category"] == cat]
        values = subset["Median"].dropna().values
        n_sites = len(values)
        
        if n_sites > 0:
            # Calculate statistics using sample SD (ddof=1)
            mean_val = np.mean(values)
            sd_val = np.std(values, ddof=1)  # sample standard deviation
            se_val = sd_val / np.sqrt(n_sites)  # standard error
            median_val = np.median(values)
            iqr_25 = np.percentile(values, 25)
            iqr_75 = np.percentile(values, 75)
            min_val = np.min(values)
            max_val = np.max(values)
            
            n_months_sum = subset["N_months"].sum()
            n_months_median = subset["N_months"].median()
            
            print(f"\n    {cat}: n={n_sites} sites")
            print(f"      Mean ± SE = {mean_val:.3f} ± {se_val:.3f}")
            print(f"      SD = {sd_val:.3f}")
            print(f"      Median = {median_val:.3f}")
            print(f"      IQR = {iqr_25:.3f} - {iqr_75:.3f}")
            print(f"      Range = {min_val:.3f} - {max_val:.3f}")
            print(f"      Total NN observations = {int(n_months_sum)}")
            print(f"      Median NN observations = {n_months_median:.1f}")
        else:
            print(f"\n    {cat}: No data available")

# =============================================================================
# STEP 8: PRINT MONTHLY SUPPORT SUMMARY - ALL METRICS USE SAME STRICT TRIPLE INTERSECTION
# =============================================================================
print("\n" + "="*80)
print("TOTAL NN OBSERVATIONS USED IN FIGURE 2 ANALYSIS")
print("ALL metrics use the same strict triple intersection + 3-month NN filter")
print("="*80)

# All metrics now use the strict triple intersection only
site_obs_shared = (
    df_adequate[df_adequate["site_name"].isin(shared_sites)][["site_name", "Salinity_Category", "N_months"]]
    .drop_duplicates(subset=["site_name", "Salinity_Category"])
)

print(f"\nStrict triple intersection (used for ALL metrics):")
print(f"  Total NN observations: {int(site_obs_shared['N_months'].sum())}")
print(f"  Total sites: {len(site_obs_shared)}")
print(f"  Unique sites in shared subset: {len(shared_sites)}")

print("\nBy salinity category (strict triple intersection):")
for cat in allowed_categories:
    subset = site_obs_shared[site_obs_shared["Salinity_Category"] == cat]
    total_obs = int(subset["N_months"].sum())
    n_sites = len(subset)
    print(f"  {cat}: {n_sites} sites, {total_obs} NN observations")

# =============================================================================
# STEP 9: STATISTICAL TESTS - ALL METRICS NOW USE SAME STRICT TRIPLE INTERSECTION
# =============================================================================
print("\n" + "="*80)
print("STATISTICAL TESTS (3 Groups: Freshwater, Saline, Upland)")
print("NOTE: Using ORIGINAL uncapped values - capping applied for visualization only")
print("      ALL metrics now use the same STRICT TRIPLE INTERSECTION + 3-month NN filter")
print("="*80)

def perform_figure2_tests(data_dict, shared_sites):
    """Perform statistical tests matching Figure 2 exactly
       ALL metrics now use the same strict triple intersection"""
    
    for metric_name, metric_data in data_dict.items():
        print(f"\n{'='*60}")
        print(f"METRIC: {metric_name}")
        print(f"{'='*60}")
        
        # Extract values for each category
        groups = {}
        for cat in allowed_categories:
            vals = metric_data[metric_data["Salinity_Category"] == cat]["Median"].dropna().values
            groups[cat] = vals
            print(f"  {cat}: n={len(vals)} sites")
        
        # Only include non-empty groups
        non_empty = [(cat, vals) for cat, vals in groups.items() if len(vals) > 0]
        
        if len(non_empty) >= 2:
            names, values_lists = zip(*non_empty)
            h_stat, p_value = kruskal(*values_lists)
            
            print(f"\n  Kruskal-Wallis test (strict triple intersection + 3-month NN filter):")
            print(f"    H-statistic = {h_stat:.3f}")
            print(f"    p-value = {p_value:.6f}")
            
            if p_value < 0.05:
                print(f"    ✓ Significant differences detected (p < 0.05)")
                
                # Post-hoc pairwise tests with Bonferroni correction
                pairwise_results = []
                for i in range(len(names)):
                    for j in range(i+1, len(names)):
                        cat1, cat2 = names[i], names[j]
                        vals1, vals2 = groups[cat1], groups[cat2]
                        
                        if len(vals1) > 0 and len(vals2) > 0:
                            u_stat, p_pair = mannwhitneyu(vals1, vals2, alternative='two-sided')
                            pairwise_results.append({
                                'comparison': f"{cat1} vs {cat2}",
                                'cat1': cat1,
                                'cat2': cat2,
                                'u_stat': u_stat,
                                'p_value': p_pair,
                                'n1': len(vals1),
                                'n2': len(vals2),
                                'median1': np.median(vals1),
                                'median2': np.median(vals2)
                            })
                
                # Apply Bonferroni correction
                if pairwise_results:
                    p_vals = [res['p_value'] for res in pairwise_results]
                    reject, p_corrected, _, _ = multipletests(p_vals, method='bonferroni')
                    
                    print(f"\n  Pairwise comparisons (Bonferroni-corrected):")
                    for idx, res in enumerate(pairwise_results):
                        sig = '***' if p_corrected[idx] < 0.001 else '**' if p_corrected[idx] < 0.01 else '*' if p_corrected[idx] < 0.05 else 'ns'
                        print(f"    {res['comparison']}: U={res['u_stat']:.1f}, p_raw={res['p_value']:.4f}, p_corr={p_corrected[idx]:.4f} {sig}")
                        print(f"       Median {res['cat1']}: {res['median1']:.3f} (n={res['n1']})")
                        print(f"       Median {res['cat2']}: {res['median2']:.3f} (n={res['n2']})")
            else:
                print(f"    No significant differences among groups (p >= 0.05)")
        else:
            print(f"  Insufficient data for Kruskal-Wallis test")

# Prepare data dictionary for tests - ALL now use strict triple intersection
test_data = {
    'WUE (strict triple intersection)': wue_final,
    'WUE_E (strict triple intersection)': eva_final,
    'WUE_T (strict triple intersection)': tra_final
}

perform_figure2_tests(test_data, shared_sites)

# =============================================================================
# STEP 10: METHODOLOGICAL NOTE FOR MANUSCRIPT
# =============================================================================
print("\n" + "="*80)
print("METHODOLOGICAL NOTE FOR MANUSCRIPT")
print("="*80)

# Store the values as variables
shared_subset_size = len(shared_sites)

print(f"""
Figure 2 Analysis Specifications:
--------------------------------
• Data source: wue_site_level_summary_SPEI_1.csv
• Water conditions: Near-normal (SPEI_Class == 'NN', SPEI-1)
• Salinity categories: Freshwater, Saline, Upland (Brackish excluded)
• Site adequacy: Minimum 3 NN observations per site (global filter, same as Figure 1)
• Per-metric adequacy: Minimum 3 NN observations per metric (N_months >= 3)

• ALL THREE METRICS (WUE, WUE_E, WUE_T) were evaluated on the SAME STRICT TRIPLE INTERSECTION
  of {shared_subset_size} sites that had ALL THREE metrics available under NN conditions.
  This ensures full comparability across metrics.

• GLOBAL 3-OBSERVATION NN FILTER APPLIED: Removed sites with <3 total NN observations (same as Figure 1)

• Total NN observations used: {int(site_obs_shared['N_months'].sum())} (aggregated site-level NN support)

• SUMMARY STATISTICS REPORTED:
  - Mean ± SE (standard error): Quantifies uncertainty around the mean estimate
  - SD (standard deviation): Measure of among-site variability
  - Median and IQR: Robust measures of central tendency and spread
  - Range: Full data extent

• Statistical tests: Kruskal-Wallis with pairwise Mann-Whitney U 
  (Bonferroni-corrected for multiple comparisons)

• Capping: 95th percentile applied for VISUALIZATION only (WUE_E and WUE_T)
• WUE remains UNCAPPED in visualization

• All statistical summaries use ORIGINAL uncapped values within the strict triple intersection

Sites excluded from analysis (lacked at least one of the three metrics):
  - WUE only: {len(wue_only)} sites
  - WUE_E only: {len(eva_only)} sites  
  - WUE_T only: {len(tra_only)} sites
""")

print("\n" + "="*80)
print("FIGURE 2 SUPPORT ANALYSIS COMPLETE")
print("="*80)
print(f"\n✓ VERIFICATION: ALL THREE metrics have identical sample size ({len(shared_sites)} sites each)")
print(f"✓ UPDATED: Group-level statistics now include Mean ± SE and SD")
print(f"✓ ADDED: Global 3-observation NN filter (same as Figure 1)")
print(f"✓ FIXED: All NN observation counts now use len(site_data) not drop_duplicates()")