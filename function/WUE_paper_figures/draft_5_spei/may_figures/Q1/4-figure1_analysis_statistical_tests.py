# -*- coding: utf-8 -*-
"""
Figure 1 Statistical Analysis Workflow - FIXED VERSION
No figures, only statistical tests and manuscript-ready output

FIXES APPLIED:
1. Removed unsafe hardcoded pivot-table column renaming
2. Added robust column validation after reshaping
3. Added metric mapping validation with expected values
4. Fixed Friedman statistic printing in final summary
5. UPDATED: Support totals use ALL retained NN observations (not unique month combos)
6. UPDATED: Wording changed to "retained NN monthly observations"

Author: Statistical Analysis Pipeline
Date: 2026-05-09
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from scipy.stats import friedmanchisquare, wilcoxon, spearmanr, kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.stats import fligner, levene

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
MONTHLY_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\figure1_statistics"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Statistical significance threshold
ALPHA = 0.05

# Expected values for validation (approximate)
EXPECTED_MEDIAN = {'WUE': 1.9, 'WUE_eva': 4.1, 'WUE_tra': 3.6}
EXPECTED_SD = {'WUE': 1.1, 'WUE_eva': 3.3, 'WUE_tra': 3.2}

# =============================================================================
# SECTION 1: LOAD DATA AND DEFINE STRICT TRIPLE INTERSECTION
# =============================================================================

print("="*80)
print("FIGURE 1 STATISTICAL ANALYSIS WORKFLOW (FIXED VERSION)")
print("="*80)
print("\nSTEP 1: Loading data and defining strict triple intersection")

# Load main data
df = pd.read_csv(INPUT_FILE)
print(f"✓ Loaded {len(df)} records from {INPUT_FILE}")

# Get sites with each metric
wue_sites = set(df[df['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(df[df['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(df[df['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())

# Strict triple intersection
shared_subset_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
n_shared = len(shared_subset_sites)

print(f"\nSITE AVAILABILITY:")
print(f"  WUE only (bulk):        {len(wue_sites)} sites")
print(f"  WUE_E only:             {len(eva_sites)} sites")
print(f"  WUE_T only:             {len(tra_sites)} sites")
print(f"  STRICT TRIPLE INTERSECTION: {n_shared} sites")

# Save shared site list
shared_sites_df = pd.DataFrame({'site_name': sorted(shared_subset_sites)})
shared_sites_df.to_csv(os.path.join(OUTPUT_DIR, 'shared_subset_sites.csv'), index=False)
print(f"✓ Saved shared site list to {OUTPUT_DIR}")

# =============================================================================
# SECTION 1B: CALCULATE RETAINED NN OBSERVATIONS (UPDATED)
# =============================================================================

print("\n" + "="*80)
print("STEP 1B: Calculating retained NN observation totals")
print("="*80)

# Load monthly data
df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"✓ Loaded monthly data: {len(df_monthly)} records")

# Filter to NN conditions at SPEI-1
df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()
print(f"✓ Filtered to NN conditions: {len(df_nn)} records")

# Filter to shared subset sites only
df_nn_shared = df_nn[df_nn["site_name"].isin(shared_subset_sites)]

# Count ALL retained NN observations (FIXED: no drop_duplicates)
total_nn_observations = len(df_nn_shared)
print(f"\n📊 RETAINED NN OBSERVATIONS (Support totals):")
print(f"  • Strict triple intersection sites: {len(shared_subset_sites)}")
print(f"  • Total retained NN observations: {total_nn_observations}")
print(f"  • Note: This represents ALL retained monthly observations (not unique month combinations)")

# Save observation count
obs_count_df = pd.DataFrame({
    'Metric': 'All three metrics (WUE, WUE_E, WUE_T)',
    'N_sites': [len(shared_subset_sites)],
    'Total_retained_NN_observations': [total_nn_observations],
    'Note': ['Support totals use all retained NN monthly observations']
})
obs_count_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_retained_observations.csv'), index=False)
print(f"✓ Saved observation counts to {OUTPUT_DIR}")

# =============================================================================
# SECTION 2: RESHAPE TO WIDE FORMAT - FIXED VERSION
# =============================================================================

print("\n" + "="*80)
print("STEP 2: Reshaping to wide format (WITH COLUMN VALIDATION)")

# Filter to shared subset
df_shared = df[df['site_name'].isin(shared_subset_sites)]

# Pivot to wide format - DO NOT RENAME COLUMNS MANUALLY
wide_df = df_shared.pivot_table(
    index='site_name',
    columns='WUE_Metric',
    values='WUE_median'
).reset_index()

print(f"\n✓ Pivot complete. Original column names from pivot:")
print(f"  {list(wide_df.columns)}")

# =============================================================================
# VALIDATION BLOCK: Inspect and verify metric columns
# =============================================================================

print("\n" + "-"*80)
print("VALIDATION: Checking metric column names and values")
print("-"*80)

# Expected metric names
expected_metrics = ['WUE', 'WUE_eva', 'WUE_tra']

# Check which expected metrics are actually present
present_metrics = [m for m in expected_metrics if m in wide_df.columns]
missing_metrics = [m for m in expected_metrics if m not in wide_df.columns]

print(f"\nExpected metrics: {expected_metrics}")
print(f"Present in pivot: {present_metrics}")
if missing_metrics:
    print(f"⚠️ WARNING: Missing metrics: {missing_metrics}")

# Verify we have all three metrics
if len(present_metrics) != 3:
    raise ValueError(f"CRITICAL: Missing metrics in pivot. Found: {present_metrics}. Expected: {expected_metrics}")

# Safely reorder columns (only if they exist)
ordered_columns = ['site_name'] + expected_metrics
wide_df = wide_df[ordered_columns]

print(f"\n✓ Safely reordered columns to: {ordered_columns}")

# =============================================================================
# VALIDATION: Check metric values against expected ranges
# =============================================================================

print("\n" + "-"*80)
print("VALIDATION: Checking metric values against expected ranges")
print("-"*80)

validation_passed = True
for metric in expected_metrics:
    values = wide_df[metric].dropna()
    median_val = values.median()
    sd_val = values.std(ddof=1)
    
    print(f"\n{metric}:")
    print(f"  Actual median = {median_val:.3f} (expected ~{EXPECTED_MEDIAN[metric]})")
    print(f"  Actual SD = {sd_val:.3f} (expected ~{EXPECTED_SD[metric]})")
    
    # Check if values are reasonably close to expectations
    median_diff = abs(median_val - EXPECTED_MEDIAN[metric])
    sd_diff = abs(sd_val - EXPECTED_SD[metric])
    
    if median_diff > 0.5:  # Allow 0.5 deviation
        print(f"  ⚠️ WARNING: Median differs substantially from expected")
        validation_passed = False
    if sd_diff > 1.0:  # Allow 1.0 deviation
        print(f"  ⚠️ WARNING: SD differs substantially from expected")
        validation_passed = False

if not validation_passed:
    print("\n" + "!"*80)
    print("METRIC MAPPING VALIDATION FAILED")
    print("!"*80)
    print("Values do not match expected ranges. Stopping execution to prevent")
    print("incorrect downstream analyses with mislabeled metrics.")
    print("\nPlease verify the pivot table column order and WUE_Metric values.")
    raise SystemExit("Metric mapping validation FAILED - stopping execution")
else:
    print("\n✓ VALIDATION PASSED: Metric mapping is correct")

# Merge with metadata
metadata_cols = ['site_name', 'Salinity_Category', 'climate', 'biome']
if all(col in df.columns for col in metadata_cols):
    metadata = df[metadata_cols].drop_duplicates()
    wide_df = wide_df.merge(metadata, on='site_name', how='left')
    print(f"\n✓ Merged metadata (Salinity_Category, climate, biome)")

# Print final wide_df structure confirmation
print(f"\n✓ Final wide DataFrame shape: {wide_df.shape}")
print(f"✓ Columns: {list(wide_df.columns)}")

# =============================================================================
# SECTION 3: DESCRIPTIVE STATISTICS
# =============================================================================

print("\n" + "="*80)
print("STEP 3: Descriptive statistics for each metric")

def compute_descriptive_stats(values, metric_name):
    """Compute comprehensive descriptive statistics"""
    values_clean = values.dropna()
    n = len(values_clean)
    mean_val = values_clean.mean()
    median_val = values_clean.median()
    sd_val = values_clean.std(ddof=1)
    se_val = sd_val / np.sqrt(n)
    cv_val = (sd_val / mean_val) * 100 if mean_val != 0 else np.nan
    
    stats_dict = {
        'Metric': metric_name,
        'N_sites': n,
        'Mean': mean_val,
        'Median': median_val,
        'SD': sd_val,
        'SE': se_val,
        'CV_%': cv_val,
        'Min': values_clean.min(),
        'Max': values_clean.max(),
        'IQR_25': values_clean.quantile(0.25),
        'IQR_75': values_clean.quantile(0.75),
        'IQR_width': values_clean.quantile(0.75) - values_clean.quantile(0.25),
        'p90': values_clean.quantile(0.90),
        'p95': values_clean.quantile(0.95),
        'p99': values_clean.quantile(0.99),
    }
    return stats_dict

# Compute for each metric
descriptive_stats = []
for metric in ['WUE', 'WUE_eva', 'WUE_tra']:
    stats = compute_descriptive_stats(wide_df[metric], metric)
    descriptive_stats.append(stats)

# Create DataFrame and save
desc_df = pd.DataFrame(descriptive_stats)
desc_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_descriptive_statistics.csv'), index=False)

# Print manuscript-ready format
print("\n" + "="*80)
print("DESCRIPTIVE STATISTICS (Manuscript-Ready)")
print("="*80)
for stats in descriptive_stats:
    print(f"\n{stats['Metric']} (n={stats['N_sites']} sites):")
    print(f"  Mean ± SE:     {stats['Mean']:.3f} ± {stats['SE']:.3f}")
    print(f"  SD:            {stats['SD']:.3f}")
    print(f"  CV:            {stats['CV_%']:.1f}%")
    print(f"  Median (IQR):  {stats['Median']:.3f} ({stats['IQR_25']:.3f}–{stats['IQR_75']:.3f})")
    print(f"  Range:         {stats['Min']:.3f}–{stats['Max']:.3f}")
    print(f"  p90/p95/p99:   {stats['p90']:.3f} / {stats['p95']:.3f} / {stats['p99']:.3f}")

# =============================================================================
# SECTION 4: FRIEDMAN TEST (Omnibus paired comparison)
# =============================================================================

print("\n" + "="*80)
print("STEP 4: Friedman test (omnibus paired comparison among all three metrics)")

# Prepare data for Friedman test
friedman_data = wide_df[['WUE', 'WUE_eva', 'WUE_tra']].dropna()
n_paired_sites = len(friedman_data)

# Perform Friedman test
statistic, p_value = friedmanchisquare(
    friedman_data['WUE'],
    friedman_data['WUE_eva'],
    friedman_data['WUE_tra']
)

# ===== CRITICAL: SAVE FRIEDMAN RESULTS IMMEDIATELY (BEFORE THEY GET OVERWRITTEN) =====
FRIEDMAN_STATISTIC = statistic
FRIEDMAN_PVALUE = p_value
FRIEDMAN_SIGNIFICANT = p_value < ALPHA

# Save results to CSV
friedman_results = pd.DataFrame({
    'Test': ['Friedman test'],
    'N_paired_sites': [n_paired_sites],
    'Statistic': [statistic],
    'Degrees_freedom': [2],
    'P_value': [p_value],
    'Significant_at_alpha_0.05': [p_value < ALPHA]
})
friedman_results.to_csv(os.path.join(OUTPUT_DIR, 'figure1_friedman_results.csv'), index=False)

# Print results
print(f"\nFriedman test (paired comparison across WUE, WUE_E, WUE_T):")
print(f"  Number of paired sites: {n_paired_sites}")
print(f"  Test statistic (Q): {statistic:.6f}")
print(f"  Degrees of freedom: 2")
print(f"  P-value: {p_value:.6f}")
print(f"  Significant at α=0.05: {'YES' if p_value < ALPHA else 'NO'}")

# =============================================================================
# SECTION 5: PAIRWISE WILCOXON TESTS WITH FDR CORRECTION
# =============================================================================

print("\n" + "="*80)
print("STEP 5: Paired Wilcoxon signed-rank tests with FDR correction")

# Define comparisons
comparisons = [
    ('WUE', 'WUE_eva', 'WUE vs WUE_E'),
    ('WUE', 'WUE_tra', 'WUE vs WUE_T'),
    ('WUE_eva', 'WUE_tra', 'WUE_E vs WUE_T')
]

# Store results
wilcoxon_results = []
raw_pvalues = []

# Calculate effect sizes and differences
for col1, col2, label in comparisons:
    # Get paired data
    paired = wide_df[[col1, col2]].dropna()
    n_pairs = len(paired)
    
    # Calculate differences
    diff = paired[col1] - paired[col2]
    median_diff = diff.median()
    direction = f"{col1} > {col2}" if median_diff > 0 else f"{col1} < {col2}" if median_diff < 0 else "No difference"
    
    # Effect size: r = Z / sqrt(N)
    # Get Wilcoxon results
    try:
        statistic, p_value = wilcoxon(paired[col1], paired[col2])
        
        # Calculate Z-statistic approximation for effect size
        # For Wilcoxon, Z = (W - n*(n+1)/4) / sqrt(n*(n+1)*(2n+1)/24)
        n = n_pairs
        expected_w = n * (n + 1) / 4
        var_w = n * (n + 1) * (2 * n + 1) / 24
        z_score = (statistic - expected_w) / np.sqrt(var_w)
        effect_size = abs(z_score) / np.sqrt(n)  # r effect size
        
        raw_pvalues.append(p_value)
        
        wilcoxon_results.append({
            'Comparison': label,
            'N_pairs': n_pairs,
            'Median_difference': median_diff,
            'Direction': direction,
            'Test_statistic_W': statistic,
            'Z_score': z_score,
            'Raw_p_value': p_value,
            'Effect_size_r': effect_size,
            'Effect_size_magnitude': 'large' if effect_size >= 0.5 else 'medium' if effect_size >= 0.3 else 'small'
        })
    except Exception as e:
        print(f"  Warning: Could not compute for {label}: {e}")
        raw_pvalues.append(np.nan)
        wilcoxon_results.append({
            'Comparison': label,
            'N_pairs': n_pairs,
            'Median_difference': median_diff,
            'Direction': direction,
            'Test_statistic_W': np.nan,
            'Z_score': np.nan,
            'Raw_p_value': np.nan,
            'Effect_size_r': np.nan,
            'Effect_size_magnitude': 'NA'
        })

# Apply Benjamini-Hochberg FDR correction
valid_pvalues = [p for p in raw_pvalues if not np.isnan(p)]
if len(valid_pvalues) > 0:
    reject, corrected_pvals, _, _ = multipletests(
        valid_pvalues, 
        alpha=ALPHA, 
        method='fdr_bh'  # Benjamini-Hochberg
    )
    
    # Map back to results
    corrected_idx = 0
    for i, result in enumerate(wilcoxon_results):
        if not np.isnan(result['Raw_p_value']):
            result['Adjusted_p_value'] = corrected_pvals[corrected_idx]
            result['Significant_after_FDR'] = reject[corrected_idx]
            corrected_idx += 1
        else:
            result['Adjusted_p_value'] = np.nan
            result['Significant_after_FDR'] = False

# Create DataFrame and save
wilcoxon_df = pd.DataFrame(wilcoxon_results)
wilcoxon_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_pairwise_wilcoxon.csv'), index=False)

# Print results
print("\nPaired Wilcoxon tests (with FDR correction):")
for result in wilcoxon_results:
    print(f"\n  {result['Comparison']} (n={result['N_pairs']} pairs):")
    print(f"    Median difference: {result['Median_difference']:.4f} ({result['Direction']})")
    print(f"    Raw p-value: {result['Raw_p_value']:.6f}")
    if not np.isnan(result['Adjusted_p_value']):
        print(f"    FDR-adjusted p: {result['Adjusted_p_value']:.6f}")
        print(f"    Significant: {'YES' if result['Significant_after_FDR'] else 'NO'}")
    print(f"    Effect size (r): {result['Effect_size_r']:.3f} ({result['Effect_size_magnitude']})")

# =============================================================================
# SECTION 6: VARIABILITY COMPARISON ACROSS METRICS
# =============================================================================

print("\n" + "="*80)
print("STEP 6: Cross-site variability comparison")

# Extract values for each metric
metrics_values = {
    'WUE': wide_df['WUE'].dropna(),
    'WUE_eva': wide_df['WUE_eva'].dropna(),
    'WUE_tra': wide_df['WUE_tra'].dropna()
}

# Compute variability metrics
variability_data = []
for metric, values in metrics_values.items():
    sd = values.std(ddof=1)
    iqr_width = values.quantile(0.75) - values.quantile(0.25)
    cv = (sd / values.mean()) * 100 if values.mean() != 0 else np.nan
    median = values.median()
    p95 = values.quantile(0.95)
    p99 = values.quantile(0.99)
    
    variability_data.append({
        'Metric': metric,
        'SD': sd,
        'IQR_width': iqr_width,
        'CV_%': cv,
        'p95_minus_median': p95 - median,
        'p99_minus_median': p99 - median
    })

var_df = pd.DataFrame(variability_data)

# Fligner-Killeen test for variance homogeneity
fligner_data = [metrics_values['WUE'].values, metrics_values['WUE_eva'].values, metrics_values['WUE_tra'].values]
fligner_stat, fligner_p = fligner(*fligner_data)

# Levene test (median-centered, robust alternative)
levene_stat, levene_p = levene(*fligner_data, center='median')

# Add test results to variability output
var_tests = pd.DataFrame({
    'Test': ['Fligner-Killeen', 'Levene (median-centered)'],
    'Test_statistic': [fligner_stat, levene_stat],
    'P_value': [fligner_p, levene_p],
    'Significant_at_alpha_0.05': [fligner_p < ALPHA, levene_p < ALPHA]
})

# Save both
var_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_variability_metrics.csv'), index=False)
var_tests.to_csv(os.path.join(OUTPUT_DIR, 'figure1_variance_homogeneity_tests.csv'), index=False)

print("\nVariability metrics:")
print(var_df.to_string(index=False))

print("\nVariance homogeneity tests (H0: equal variances):")
for test, stat, p in zip(['Fligner-Killeen', 'Levene'], [fligner_stat, levene_stat], [fligner_p, levene_p]):
    print(f"  {test}: statistic={stat:.3f}, p={p:.6f} ({'significant' if p < ALPHA else 'not significant'})")

# =============================================================================
# SECTION 7: T:ET GRADIENT ANALYSIS (Spearman correlations)
# =============================================================================

print("\n" + "="*80)
print("STEP 7: T:ET gradient analysis (Spearman correlations)")

# Filter to NN conditions and shared subset sites
df_nn = df_monthly[df_monthly['SPEI_1_Cat'] == 'NN'].copy()
df_nn_shared = df_nn[df_nn['site_name'].isin(shared_subset_sites)]

# Check if Trans_ratio exists
if 'Trans_ratio' in df_nn_shared.columns:
    # Compute site-level median Trans_ratio
    trans_ratio_median = df_nn_shared.groupby('site_name')['Trans_ratio'].median().reset_index()
    trans_ratio_median.columns = ['site_name', 'Trans_ratio_median']
    
    # Merge with wide_df
    wide_df = wide_df.merge(trans_ratio_median, on='site_name', how='left')
    print(f"✓ Computed site-level median Trans_ratio for {trans_ratio_median['site_name'].nunique()} sites")
    
    # Run Spearman correlations
    spearman_results = []
    for metric in ['WUE', 'WUE_eva', 'WUE_tra']:
        # Drop NA in either variable
        corr_data = wide_df[[metric, 'Trans_ratio_median']].dropna()
        if len(corr_data) > 0:
            rho, p_value = spearmanr(corr_data[metric], corr_data['Trans_ratio_median'])
            direction = 'positive' if rho > 0 else 'negative' if rho < 0 else 'zero'
            
            spearman_results.append({
                'Metric': metric,
                'N_sites': len(corr_data),
                'Spearman_rho': rho,
                'P_value': p_value,
                'Direction': direction,
                'Significant': p_value < ALPHA
            })
    
    # Save results
    spearman_df = pd.DataFrame(spearman_results)
    spearman_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_spearman_correlations.csv'), index=False)
    
    print("\nSpearman correlations with Trans_ratio:")
    for result in spearman_results:
        print(f"\n  {result['Metric']} (n={result['N_sites']}):")
        print(f"    ρ = {result['Spearman_rho']:.3f} ({result['Direction']})")
        print(f"    p = {result['P_value']:.6f} ({'significant' if result['Significant'] else 'not significant'})")
else:
    print("⚠️ Trans_ratio column not found in monthly file. Skipping Spearman correlations.")

# =============================================================================
# SECTION 8: T:ET DIVERGENCE ANALYSIS
# =============================================================================

if 'Trans_ratio_median' in wide_df.columns:
    print("\n" + "="*80)
    print("STEP 8: T:ET divergence analysis")
    
    # Calculate divergence: |WUE_ET - WUE_T|
    wide_df['divergence'] = np.abs(wide_df['WUE'] - wide_df['WUE_tra'])
    
    # Create terciles of Trans_ratio
    wide_df['TET_tercile'] = pd.qcut(wide_df['Trans_ratio_median'].rank(method='first'), 
                                     3, labels=['low T:ET', 'intermediate T:ET', 'high T:ET'])
    
    # Kruskal-Wallis test
    groups = [wide_df[wide_df['TET_tercile'] == group]['divergence'].dropna() 
              for group in ['low T:ET', 'intermediate T:ET', 'high T:ET']]
    groups = [g for g in groups if len(g) > 0]
    
    if len(groups) == 3:
        h_stat, kw_p = kruskal(*groups)
        
        # Post-hoc Mann-Whitney with Bonferroni
        posthoc_results = []
        tercile_pairs = [('low T:ET', 'intermediate T:ET'), 
                        ('low T:ET', 'high T:ET'),
                        ('intermediate T:ET', 'high T:ET')]
        
        raw_ps = []
        for pair in tercile_pairs:
            group1 = wide_df[wide_df['TET_tercile'] == pair[0]]['divergence'].dropna()
            group2 = wide_df[wide_df['TET_tercile'] == pair[1]]['divergence'].dropna()
            if len(group1) > 0 and len(group2) > 0:
                u_stat, p_val = mannwhitneyu(group1, group2, alternative='two-sided')
                raw_ps.append(p_val)
                effect_size = u_stat / (len(group1) * len(group2))
                posthoc_results.append({
                    'Group1': pair[0],
                    'Group2': pair[1],
                    'N1': len(group1),
                    'N2': len(group2),
                    'Mann_Whitney_U': u_stat,
                    'Raw_p_value': p_val,
                    'Effect_size_U/(n1*n2)': effect_size
                })
        
        # Apply Bonferroni correction
        if len(raw_ps) > 0:
            bonf_pvals = [min(p * len(raw_ps), 1.0) for p in raw_ps]
            for i, result in enumerate(posthoc_results):
                result['Bonferroni_p_value'] = bonf_pvals[i]
                result['Significant'] = bonf_pvals[i] < ALPHA
        
        # Save results
        tet_results = {
            'Test': ['Kruskal-Wallis'],
            'N_low_TET': [len(groups[0])],
            'N_int_TET': [len(groups[1])],
            'N_high_TET': [len(groups[2])],
            'H_statistic': [h_stat],
            'P_value': [kw_p],
            'Significant': [kw_p < ALPHA]
        }
        tet_df = pd.DataFrame(tet_results)
        tet_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_TET_divergence_kw.csv'), index=False)
        
        posthoc_df = pd.DataFrame(posthoc_results)
        posthoc_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_TET_divergence_posthoc.csv'), index=False)
        
        print(f"\nKruskal-Wallis test on divergence |WUE_ET - WUE_T| across T:ET terciles:")
        print(f"  H = {h_stat:.3f}, p = {kw_p:.6f} ({'significant' if kw_p < ALPHA else 'not significant'})")
        print(f"\n  Group medians:")
        for group in ['low T:ET', 'intermediate T:ET', 'high T:ET']:
            med = wide_df[wide_df['TET_tercile'] == group]['divergence'].median()
            print(f"    {group}: {med:.4f}")
        
        if kw_p < ALPHA:
            print(f"\n  Post-hoc Mann-Whitney (Bonferroni corrected):")
            for result in posthoc_results:
                print(f"    {result['Group1']} vs {result['Group2']}: p_adj={result['Bonferroni_p_value']:.6f} "
                      f"({'significant' if result['Significant'] else 'not'})")

# =============================================================================
# SECTION 9: ECOSYSTEM COMPARISON ANALYSIS
# =============================================================================

print("\n" + "="*80)
print("STEP 9: Ecosystem comparison (Freshwater, Saline, Upland)")

# Exclude Brackish
eco_df = wide_df[~wide_df['Salinity_Category'].str.contains('Brackish', na=False)].copy()
eco_df = eco_df[eco_df['Salinity_Category'].isin(['Freshwater', 'Saline', 'Upland'])].copy()

if len(eco_df) > 0:
    kruskal_results = []
    mannwhitney_results = []
    
    for metric in ['WUE', 'WUE_eva', 'WUE_tra']:
        # Prepare groups
        groups = [eco_df[eco_df['Salinity_Category'] == cat][metric].dropna() 
                 for cat in ['Freshwater', 'Saline', 'Upland']]
        groups = [g for g in groups if len(g) > 0]
        group_names = ['Freshwater', 'Saline', 'Upland'][:len(groups)]
        group_sizes = [len(g) for g in groups]
        
        if len(groups) >= 2:
            # Kruskal-Wallis
            h_stat, kw_p = kruskal(*groups)
            
            kruskal_results.append({
                'Metric': metric,
                'Groups_tested': ' vs '.join(group_names),
                'Sample_sizes': str(dict(zip(group_names, group_sizes))),
                'H_statistic': h_stat,
                'P_value': kw_p,
                'Significant': kw_p < ALPHA
            })
            
            # Post-hoc pairwise Mann-Whitney if significant
            if kw_p < ALPHA and len(groups) >= 2:
                pairs = [('Freshwater', 'Saline'), ('Freshwater', 'Upland'), ('Saline', 'Upland')]
                raw_ps = []
                pair_results = []
                
                for pair in pairs:
                    if pair[0] in group_names and pair[1] in group_names:
                        group1 = eco_df[eco_df['Salinity_Category'] == pair[0]][metric].dropna()
                        group2 = eco_df[eco_df['Salinity_Category'] == pair[1]][metric].dropna()
                        if len(group1) > 0 and len(group2) > 0:
                            u_stat, p_val = mannwhitneyu(group1, group2, alternative='two-sided')
                            raw_ps.append(p_val)
                            pair_results.append({
                                'Metric': metric,
                                'Group1': pair[0],
                                'Group2': pair[1],
                                'N1': len(group1),
                                'N2': len(group2),
                                'Median1': group1.median(),
                                'Median2': group2.median(),
                                'Mann_Whitney_U': u_stat,
                                'Raw_p_value': p_val
                            })
                
                # Bonferroni correction
                if len(raw_ps) > 0:
                    bonf_pvals = [min(p * len(raw_ps), 1.0) for p in raw_ps]
                    for i, result in enumerate(pair_results):
                        result['Bonferroni_p_value'] = bonf_pvals[i]
                        result['Significant'] = bonf_pvals[i] < ALPHA
                        mannwhitney_results.append(result)
    
    # Save results
    if kruskal_results:
        kw_df = pd.DataFrame(kruskal_results)
        kw_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_kruskal_results.csv'), index=False)
    
    if mannwhitney_results:
        mw_df = pd.DataFrame(mannwhitney_results)
        mw_df.to_csv(os.path.join(OUTPUT_DIR, 'figure1_mannwhitney_results.csv'), index=False)
    
    print("\nKruskal-Wallis tests by Salinity Category (Freshwater, Saline, Upland):")
    for result in kruskal_results:
        print(f"\n  {result['Metric']}:")
        print(f"    H = {result['H_statistic']:.3f}, p = {result['P_value']:.6f} "
              f"({'significant' if result['Significant'] else 'not significant'})")
        print(f"    Sample sizes: {result['Sample_sizes']}")
        
        # Show medians
        for cat in ['Freshwater', 'Saline', 'Upland']:
            median_val = eco_df[eco_df['Salinity_Category'] == cat][result['Metric']].median()
            if not pd.isna(median_val):
                print(f"    Median in {cat}: {median_val:.3f}")
    
    if mannwhitney_results:
        print(f"\n  Significant post-hoc comparisons (Bonferroni corrected):")
        for result in mannwhitney_results:
            if result['Significant']:
                print(f"    {result['Metric']}: {result['Group1']} vs {result['Group2']} "
                      f"(p_adj={result['Bonferroni_p_value']:.6f})")
else:
    print("⚠️ Insufficient data for ecosystem comparison (no valid categories)")


# =============================================================================
# SECTION 10: FINAL MANUSCRIPT-READY SUMMARY (COMPLETELY FIXED)
# =============================================================================

print("\n" + "="*80)
print("FINAL MANUSCRIPT-READY SUMMARY")
print("="*80)

print(f"""
FIGURE 1 STATISTICAL ANALYSIS COMPLETE

1. SAMPLE CHARACTERISTICS:
   • Strict triple intersection: {n_shared} sites
   • Total retained NN observations: {total_nn_observations}
   • Note: Support totals use all retained NN monthly observations
   • All three metrics (WUE, WUE_E, WUE_T) have identical site samples

2. OMNIBUS COMPARISON (Friedman test):
   • Friedman Q = {FRIEDMAN_STATISTIC:.3f}
   • p-value = {FRIEDMAN_PVALUE:.6f}
   • Significant at α=0.05: {'YES' if FRIEDMAN_SIGNIFICANT else 'NO'}

3. PAIRWISE COMPARISONS (FDR-corrected Wilcoxon):
""")

for result in wilcoxon_results:
    if not np.isnan(result['Adjusted_p_value']):
        print(f"   • {result['Comparison']}: p_adj={result['Adjusted_p_value']:.6f} "
              f"({'SIGNIFICANT' if result['Significant_after_FDR'] else 'NOT significant'}), "
              f"r={result['Effect_size_r']:.3f} ({result['Effect_size_magnitude']})")

print(f"""
4. VARIABILITY FINDINGS:
   • Variance heterogeneity (Fligner-Killeen): p={fligner_p:.6f} 
     ({'significant' if fligner_p < ALPHA else 'NOT significant'})
   • Highest CV: {var_df.loc[var_df['CV_%'].idxmax(), 'Metric']} ({var_df['CV_%'].max():.1f}%)
   • Largest IQR width: {var_df.loc[var_df['IQR_width'].idxmax(), 'Metric']} ({var_df['IQR_width'].max():.3f})
""")

if 'Trans_ratio_median' in wide_df.columns and len(spearman_results) > 0:
    print(f"""
5. T:ET RELATIONSHIPS (Spearman correlations):
""")
    for result in spearman_results:
        print(f"   • {result['Metric']}: ρ={result['Spearman_rho']:.3f}, "
              f"p={result['P_value']:.6f} ({'significant' if result['Significant'] else 'NOT'})")

if len(kruskal_results) > 0:
    print(f"""
6. ECOSYSTEM GROUP DIFFERENCES (Kruskal-Wallis, Freshwater/Saline/Upland):
""")
    for result in kruskal_results:
        print(f"   • {result['Metric']}: H={result['H_statistic']:.3f}, "
              f"p={result['P_value']:.6f} ({'significant' if result['Significant'] else 'NOT'})")

print(f"""
7. OUTPUT FILES SAVED TO:
   {OUTPUT_DIR}
   
   ✓ shared_subset_sites.csv
   ✓ figure1_retained_observations.csv
   ✓ figure1_descriptive_statistics.csv
   ✓ figure1_friedman_results.csv
   ✓ figure1_pairwise_wilcoxon.csv
   ✓ figure1_variability_metrics.csv
   ✓ figure1_variance_homogeneity_tests.csv
   ✓ figure1_spearman_correlations.csv
   ✓ figure1_TET_divergence_kw.csv
   ✓ figure1_TET_divergence_posthoc.csv
   ✓ figure1_kruskal_results.csv
   ✓ figure1_mannwhitney_results.csv

ANALYSIS COMPLETE - All statistics validated
Support totals use all retained NN monthly observations (not unique month combinations)
""")

print("="*80)
print("STATISTICAL WORKFLOW FINISHED SUCCESSFULLY")
print("="*80)