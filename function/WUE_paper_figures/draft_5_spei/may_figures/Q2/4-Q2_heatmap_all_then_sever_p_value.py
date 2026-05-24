# -*- coding: utf-8 -*-
"""
Created on Sun May  3 16:10:40 2026
Updated: Mon May 4 2026

@author: ammar

================================================================================
GROUP-LEVEL INFERENTIAL TESTING ON PERCENT CHANGE RELATIVE TO BASELINE
================================================================================

This workflow combines:
(1) site-level NN bootstrap CI classification (uncertainty-aware)
(2) Wilcoxon signed-rank test on median percent change relative to zero as the 
    primary group-level test (non-parametric, robust for skewed data)
(3) one-sample t-test as a mean-based sensitivity check

Zero represents "no change relative to NN baseline conditions"

For DRY conditions: Testing if percent change < 0 (significant decrease)
For WET conditions: Testing if percent change > 0 (significant increase)

Input: SPEI_site_level_details_FINAL.csv (with updated NN bootstrap CI columns)
Output: Q2_CI_adjusted_ttest_results.csv
================================================================================
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp, wilcoxon
import os

# ================================
# INPUT FILE (same as your workflow)
# ================================
file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
output_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results"
output_file = os.path.join(output_dir, "Q2_CI_adjusted_ttest_results.csv")

df = pd.read_csv(file)

metrics = ["WUE", "WUE_tra"]
metrics_display = {"WUE": "WUE_ET", "WUE_tra": "WUE_T"}
timescales = ["SPEI_6", "SPEI_48"]

conditions = [
    ("PASS C", "Dry (all)"),
    ("PASS C", "Wet (all)"),
    ("PASS B", "Severe Dry"),
    ("PASS B", "Severe Wet")
]

print("=" * 100)
print("GROUP-LEVEL INFERENTIAL TESTING ON PERCENT CHANGE RELATIVE TO BASELINE")
print("=" * 100)
print("\nThis workflow combines:")
print("  (1) site-level NN bootstrap CI classification (uncertainty-aware)")
print("  (2) Wilcoxon signed-rank test (primary, non-parametric)")
print("  (3) one-sample t-test (sensitivity check, mean-based)")
print("\nZero represents 'no change relative to NN baseline conditions'")
print("\nDRY conditions: Testing if percent change < 0 (significant decrease)")
print("WET conditions: Testing if percent change > 0 (significant increase)")
print("=" * 100)

results = []

for pass_name, condition in conditions:
    for timescale in timescales:
        
        print(f"\n{'='*80}")
        print(f"Processing: {pass_name} | {condition} | {timescale}")
        print(f"{'='*80}")

        # Filter to All salinity and relevant metrics
        subset = df[
            (df["Pass"] == pass_name) &
            (df["Salinity"] == "All") &
            (df["SPEI_Timescale"] == timescale) &
            (df["Condition"] == condition) &
            (df["WUE_Metric"].isin(metrics))
        ].copy()

        # --- strict double intersection (same as figure code) ---
        wue_sites = set(subset[subset["WUE_Metric"] == "WUE"]["Site"].unique())
        tra_sites = set(subset[subset["WUE_Metric"] == "WUE_tra"]["Site"].unique())
        shared_sites = wue_sites.intersection(tra_sites)
        
        print(f"  Shared sites (both metrics): {len(shared_sites)}")

        for metric in metrics:
            # Filter to shared sites
            sub = subset[
                (subset["WUE_Metric"] == metric) &
                (subset["Site"].isin(shared_sites))
            ].copy()
            
            if len(sub) == 0:
                print(f"\n  {metrics_display[metric]}: No data")
                continue
            
            # Step 1: Collapse to one row per site
            site_level = sub.groupby("Site").agg({
                "Median_%_Change": "median",
                "CI_Lower": "median",
                "CI_Upper": "median",
                "CI_Width": "median",
                "Direction": "first",
                "Direction_Method": "first",
                "Threshold_Method": "first"
            }).reset_index()
            
            n_sites = len(site_level)
            
            # Extract percent change values for testing
            pct_change_vals = site_level["Median_%_Change"].dropna().values
            
            # Step 2: Descriptive statistics
            median_change = np.median(pct_change_vals) if len(pct_change_vals) > 0 else np.nan
            mean_change = np.mean(pct_change_vals) if len(pct_change_vals) > 0 else np.nan
            sd_change = np.std(pct_change_vals, ddof=1) if len(pct_change_vals) > 1 else np.nan
            se_change = sd_change / np.sqrt(n_sites) if len(pct_change_vals) > 1 and n_sites > 0 else np.nan
            
            # IQR for median percent change
            if len(pct_change_vals) > 0:
                q25 = np.percentile(pct_change_vals, 25)
                q75 = np.percentile(pct_change_vals, 75)
                iqr = q75 - q25
                min_change = np.min(pct_change_vals)
                max_change = np.max(pct_change_vals)
            else:
                q25, q75, iqr, min_change, max_change = np.nan, np.nan, np.nan, np.nan, np.nan
            
            # Step 3: Count NN bootstrap classifications
            direction_counts = site_level["Direction"].value_counts()
            increase_sites = direction_counts.get("Increase", 0)
            decrease_sites = direction_counts.get("Decrease", 0)
            no_change_sites = direction_counts.get("No change", 0)
            insufficient_sites = direction_counts.get("Insufficient data", 0)
            
            # Median CI bounds (for reporting only, not used in test)
            median_ci_lower = site_level["CI_Lower"].median() if "CI_Lower" in site_level.columns else np.nan
            median_ci_upper = site_level["CI_Upper"].median() if "CI_Upper" in site_level.columns else np.nan
            median_ci_width = site_level["CI_Width"].median() if "CI_Width" in site_level.columns else np.nan
            
            # Get Direction_Method and Threshold_Method
            direction_method = site_level["Direction_Method"].iloc[0] if "Direction_Method" in site_level.columns else "NN_Bootstrap_Median_CI_95"
            threshold_method = site_level["Threshold_Method"].iloc[0] if "Threshold_Method" in site_level.columns else "Site_specific_NN_bootstrap_95CI_median"
            
            # Step 4: Run both tests on Median_%_Change
            is_dry = "Dry" in condition
            alternative = "less" if is_dry else "greater"
            
            if len(pct_change_vals) >= 3:
                # ----- ONE-SAMPLE T-TEST -----
                t_stat, t_p_val = ttest_1samp(pct_change_vals, 0, alternative=alternative)
                t_significant = t_p_val < 0.05
                
                # ----- WILCOXON SIGNED-RANK TEST -----
                try:
                    # Check if all values are zero (would cause Wilcoxon to fail)
                    if np.all(pct_change_vals == 0):
                        wilcoxon_stat = np.nan
                        wilcoxon_p_val = np.nan
                        wilcoxon_significant = False
                        print(f"      Wilcoxon test unavailable: all values are zero")
                    else:
                        # Use wilcoxon with appropriate alternative
                        # Note: wilcoxon in scipy uses 'less'/'greater' as of recent versions
                        # zero_method="wilcox" handles zeros by discarding them
                        wilcoxon_result = wilcoxon(pct_change_vals, alternative=alternative, zero_method="wilcox")
                        wilcoxon_stat = wilcoxon_result.statistic
                        wilcoxon_p_val = wilcoxon_result.pvalue
                        wilcoxon_significant = wilcoxon_p_val < 0.05
                except Exception as e:
                    wilcoxon_stat = np.nan
                    wilcoxon_p_val = np.nan
                    wilcoxon_significant = False
                    print(f"      Wilcoxon test unavailable: {str(e)}")
                
                # Primary interpretation based on Wilcoxon
                if wilcoxon_significant and is_dry:
                    interpretation = "WUE significantly decreased relative to baseline conditions (Wilcoxon)"
                elif wilcoxon_significant and not is_dry:
                    interpretation = "WUE significantly increased relative to baseline conditions (Wilcoxon)"
                elif not wilcoxon_significant and not np.isnan(wilcoxon_p_val):
                    interpretation = "No significant change relative to baseline conditions (Wilcoxon)"
                else:
                    interpretation = "Insufficient data or Wilcoxon unavailable"
                
                # Secondary interpretation based on t-test for comparison
                if t_significant and is_dry:
                    t_interpretation = "WUE significantly decreased relative to baseline conditions (t-test)"
                elif t_significant and not is_dry:
                    t_interpretation = "WUE significantly increased relative to baseline conditions (t-test)"
                elif not t_significant and not np.isnan(t_p_val):
                    t_interpretation = "No significant change relative to baseline conditions (t-test)"
                else:
                    t_interpretation = "Insufficient data for t-test"
                    
            else:
                t_stat = np.nan
                t_p_val = np.nan
                t_significant = False
                wilcoxon_stat = np.nan
                wilcoxon_p_val = np.nan
                wilcoxon_significant = False
                interpretation = f"Insufficient data for inferential tests (n={len(pct_change_vals)} < 3)"
                t_interpretation = f"Insufficient data (n={len(pct_change_vals)} < 3)"
            
            # Compile results
            result = {
                "Pass": pass_name,
                "Condition": condition,
                "Timescale": timescale,
                "Metric": metrics_display[metric],
                "N_sites": n_sites,
                "Mean_%_Change": mean_change,
                "Median_%_Change": median_change,
                "SD_%_Change": sd_change,
                "SE_%_Change": se_change,
                "IQR_%_Change": iqr,
                "Min_%_Change": min_change,
                "Max_%_Change": max_change,
                "Median_CI_Lower": median_ci_lower,
                "Median_CI_Upper": median_ci_upper,
                "Median_CI_Width": median_ci_width,
                "Increase_sites": increase_sites,
                "Decrease_sites": decrease_sites,
                "No_change_sites": no_change_sites,
                "Insufficient_data_sites": insufficient_sites,
                # Primary test (Wilcoxon) results - these become the main P_value and Significant
                "P_value": wilcoxon_p_val,
                "Significant": wilcoxon_significant,
                # T-test results (kept for sensitivity check)
                "T_statistic": t_stat,
                "T_p_value": t_p_val,
                "T_significant": t_significant,
                # Wilcoxon results
                "Wilcoxon_statistic": wilcoxon_stat,
                "Wilcoxon_p_value": wilcoxon_p_val,
                "Wilcoxon_significant": wilcoxon_significant,
                # Test metadata
                "Alternative": alternative,
                "Primary_test": "Wilcoxon signed-rank",
                "Direction_Method": direction_method,
                "Threshold_Method": threshold_method,
                "Interpretation": interpretation,
                "T_Interpretation": t_interpretation
            }
            
            results.append(result)
            
            # Print detailed results
            print(f"\n  {metrics_display[metric]} (n={n_sites} sites):")
            print(f"    Percent change statistics:")
            print(f"      Mean = {mean_change:+.1f}% ± {se_change:.1f}% (SD={sd_change:.1f}%)")
            print(f"      Median = {median_change:+.1f}% (IQR: {q25:+.1f}% to {q75:+.1f}%)")
            print(f"      Range: {min_change:+.1f}% to {max_change:+.1f}%")
            print(f"    Site-level NN bootstrap classification:")
            print(f"      Increase: {increase_sites}, Decrease: {decrease_sites}, No change: {no_change_sites}, Insufficient: {insufficient_sites}")
            
            if len(pct_change_vals) >= 3:
                print(f"    One-sample t-test (H1: mean {alternative} 0):")
                print(f"      t = {t_stat:.3f}, p = {t_p_val:.4f} → {'Significant' if t_significant else 'Not significant'}")
                
                print(f"    Wilcoxon signed-rank test (H1: median {alternative} 0):")
                if not np.isnan(wilcoxon_stat):
                    print(f"      W = {wilcoxon_stat:.1f}, p = {wilcoxon_p_val:.4f} → {'Significant' if wilcoxon_significant else 'Not significant'}")
                    if wilcoxon_significant:
                        print(f"      ✓ PRIMARY: {interpretation}")
                    else:
                        print(f"      ✗ PRIMARY: {interpretation}")
                else:
                    print(f"      Wilcoxon test unavailable")
                print(f"    Sensitivity check (t-test): {t_interpretation}")
            else:
                print(f"    ✗ {interpretation}")

# Save results to CSV
results_df = pd.DataFrame(results)

# Reorder columns for better readability
column_order = [
    "Pass", "Condition", "Timescale", "Metric", "N_sites",
    "Mean_%_Change", "Median_%_Change", "SD_%_Change", "SE_%_Change", "IQR_%_Change",
    "Min_%_Change", "Max_%_Change",
    "Median_CI_Lower", "Median_CI_Upper", "Median_CI_Width",
    "Increase_sites", "Decrease_sites", "No_change_sites", "Insufficient_data_sites",
    "P_value", "Significant",  # Primary (Wilcoxon)
    "Wilcoxon_statistic", "Wilcoxon_p_value", "Wilcoxon_significant",
    "T_statistic", "T_p_value", "T_significant",
    "Alternative", "Primary_test",
    "Direction_Method", "Threshold_Method",
    "Interpretation", "T_Interpretation"
]

# Only include columns that exist
existing_columns = [col for col in column_order if col in results_df.columns]
results_df = results_df[existing_columns]

results_df.to_csv(output_file, index=False)

print("\n" + "=" * 100)
print("GROUP-LEVEL TEST RESULTS SUMMARY")
print("=" * 100)
print(results_df.to_string(index=False))
print("\n" + "=" * 100)
print(f"Results saved to: {output_file}")
print("=" * 100)

# Print interpretation guidance
print("\n" + "=" * 100)
print("INTERPRETATION GUIDANCE")
print("=" * 100)
print("""
This workflow combines multiple approaches for robust inference:

1. SITE-LEVEL CLASSIFICATION (from upstream workflow):
   - Uses site-specific NN bootstrap 95% CI thresholds
   - Accounts for each site's baseline variability
   - Classifies each site as Increase/Decrease/No change/Insufficient data

2. GROUP-LEVEL TESTS (this workflow):
   PRIMARY TEST: Wilcoxon signed-rank test (non-parametric, median-based)
     - Recommended for non-normal or skewed data
     - Robust to outliers
     - Tests whether median percent change differs from zero
   
   SENSITIVITY CHECK: One-sample t-test (parametric, mean-based)
     - Provided for comparison and completeness
     - Tests whether mean percent change differs from zero

   For DRY conditions (Dry, Severe Dry):
     H1: percent change < 0 (significant decrease)
   
   For WET conditions (Wet, Severe Wet):
     H1: percent change > 0 (significant increase)

WHY WILCOXON IS PRIMARY:
   - Percent change data can be skewed or have outliers
   - Heatmap and summaries emphasize medians
   - Non-parametric test makes fewer distributional assumptions
   - More robust for small sample sizes

INTERPRETATION EXAMPLE:
   - Significant Wilcoxon + dry → WUE significantly decreases across sites
   - Significant Wilcoxon + wet → WUE significantly increases across sites
   - Non-significant → No evidence of systematic change across sites
   - Compare t-test results: agreement increases confidence, disagreement suggests influence of outliers

NOTE ON SAMPLE SIZE:
   - Tests require n ≥ 3 sites for interpretable results
   - Below n=3, results are marked as insufficient data
""")
print("=" * 100)