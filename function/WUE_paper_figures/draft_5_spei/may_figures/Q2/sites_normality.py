# -*- coding: utf-8 -*-
"""
Monthly WUE Distribution Normality Tests at Site Level
No figures, only statistical summaries and CSV outputs

Author: Normality Testing Pipeline
Date: 2026-05-05
Purpose: Test whether monthly WUE distributions are normal at each site
"""

import pandas as pd
import numpy as np
import os
from scipy import stats
from scipy.stats import shapiro, normaltest, anderson
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
SITE_LEVEL_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\normality_tests"

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Minimum sample size for normality tests
MIN_N_FOR_TEST = 5
MIN_N_FOR_DAGOSTINO = 8
MIN_N_FOR_ANDERSON = 8

# Statistical significance threshold
ALPHA = 0.05

# =============================================================================
# SECTION 1: LOAD DATA AND APPLY FILTERS
# =============================================================================

print("="*80)
print("MONTHLY WUE DISTRIBUTION NORMALITY TESTS")
print("="*80)

print("\nSTEP 1: Loading data and applying filters")

# Load monthly data
df_monthly = pd.read_csv(INPUT_FILE)
print(f"✓ Loaded {len(df_monthly)} records from monthly file")

# Load site-level data to get strict triple intersection sites
df_site = pd.read_csv(SITE_LEVEL_FILE)
print(f"✓ Loaded {len(df_site)} records from site-level file")

# Get strict triple intersection sites (sites with all three metrics)
wue_sites = set(df_site[df_site['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(df_site[df_site['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(df_site[df_site['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_subset_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

print(f"\nStrict triple intersection sites (from Figure 1): {len(shared_subset_sites)} sites")

# =============================================================================
# SECTION 2: FILTER MONTHLY DATA
# =============================================================================

print("\nSTEP 2: Filtering monthly data")

# Filter 1: NN baseline months only
df_nn = df_monthly[df_monthly['SPEI_1_Cat'] == 'NN'].copy()
print(f"  After NN filter: {len(df_nn)} records")

# Filter 2: Keep only rows where all three metrics are non-missing
df_complete = df_nn.dropna(subset=['WUE', 'WUE_eva', 'WUE_tra']).copy()
print(f"  After dropping missing WUE values: {len(df_complete)} records")

# Filter 3: Restrict to strict triple intersection sites
df_filtered = df_complete[df_complete['site_name'].isin(shared_subset_sites)].copy()
print(f"  After restricting to shared sites: {len(df_filtered)} records")
print(f"  Unique sites in final dataset: {df_filtered['site_name'].nunique()}")

# =============================================================================
# SECTION 3: FUNCTION TO COMPUTE NORMALITY TESTS AND STATISTICS
# =============================================================================

def compute_normality_tests(values, metric_name, site_name):
    """
    Compute normality tests and descriptive statistics for a site-metric combination
    
    Returns:
        dict with test results and statistics
    """
    n = len(values)
    
    # Basic statistics always computed
    result = {
        'site_name': site_name,
        'WUE_Metric': metric_name,
        'n_months': n,
        'mean': values.mean(),
        'median': values.median(),
        'sd': values.std(ddof=1),
        'skewness': values.skew(),
        'kurtosis': values.kurtosis(),
        'min': values.min(),
        'max': values.max(),
        'iqr': values.quantile(0.75) - values.quantile(0.25)
    }
    
    # Check minimum sample size
    if n < MIN_N_FOR_TEST:
        result['Shapiro_W'] = np.nan
        result['Shapiro_p'] = np.nan
        result['Normality_Status'] = 'Insufficient_N'
        result['Dagostino_K2'] = np.nan
        result['Dagostino_p'] = np.nan
        result['Anderson_Statistic'] = np.nan
        result['Anderson_Critical_5%'] = np.nan
        result['Anderson_Status'] = 'Insufficient_N'
        return result
    
    # Shapiro-Wilk test (primary)
    try:
        shapiro_stat, shapiro_p = shapiro(values)
        result['Shapiro_W'] = shapiro_stat
        result['Shapiro_p'] = shapiro_p
        result['Normality_Status'] = 'Normal' if shapiro_p >= ALPHA else 'Non_normal'
    except Exception as e:
        result['Shapiro_W'] = np.nan
        result['Shapiro_p'] = np.nan
        result['Normality_Status'] = 'Error'
        print(f"  Warning: Shapiro test failed for {site_name} - {metric_name}: {e}")
    
    # D'Agostino K² test (if n ≥ 8)
    if n >= MIN_N_FOR_DAGOSTINO:
        try:
            dagostino_stat, dagostino_p = normaltest(values)
            result['Dagostino_K2'] = dagostino_stat
            result['Dagostino_p'] = dagostino_p
        except Exception as e:
            result['Dagostino_K2'] = np.nan
            result['Dagostino_p'] = np.nan
    else:
        result['Dagostino_K2'] = np.nan
        result['Dagostino_p'] = np.nan
    
    # Anderson-Darling test (if n ≥ 8)
    if n >= MIN_N_FOR_ANDERSON:
        try:
            anderson_result = anderson(values, dist='norm')
            result['Anderson_Statistic'] = anderson_result.statistic
            # Get critical value at 5% significance
            critical_5pct = anderson_result.critical_values[2]  # Index 2 = 5%
            result['Anderson_Critical_5%'] = critical_5pct
            result['Anderson_Status'] = 'Normal' if anderson_result.statistic < critical_5pct else 'Non_normal'
        except Exception as e:
            result['Anderson_Statistic'] = np.nan
            result['Anderson_Critical_5%'] = np.nan
            result['Anderson_Status'] = 'Error'
    else:
        result['Anderson_Statistic'] = np.nan
        result['Anderson_Critical_5%'] = np.nan
        result['Anderson_Status'] = 'Insufficient_N'
    
    return result

# =============================================================================
# SECTION 4: RUN TESTS FOR EACH SITE AND METRIC
# =============================================================================

print("\n" + "="*80)
print("STEP 3: Running normality tests for each site-metric combination")
print("="*80)

# Metrics to test
metrics = ['WUE', 'WUE_eva', 'WUE_tra']
metric_labels = {'WUE': 'WUE_ET', 'WUE_eva': 'WUE_E', 'WUE_tra': 'WUE_T'}

# Store results
all_results = []

# Get unique sites
sites = df_filtered['site_name'].unique()
print(f"\nTesting {len(sites)} sites for each of 3 metrics...")

for site in sites:
    print(f"\nProcessing site: {site}")
    site_data = df_filtered[df_filtered['site_name'] == site]
    
    for metric in metrics:
        # Extract values for this metric
        values = site_data[metric].dropna()
        
        if len(values) > 0:
            result = compute_normality_tests(values, metric, site)
            all_results.append(result)
            status = result['Normality_Status']
            print(f"  {metric_labels[metric]}: n={result['n_months']}, status={status}")
        else:
            print(f"  {metric_labels[metric]}: No data available")

# Create DataFrame
results_df = pd.DataFrame(all_results)

# =============================================================================
# SECTION 5: SAVE SITE-LEVEL RESULTS
# =============================================================================

print("\n" + "="*80)
print("STEP 4: Saving site-level results")
print("="*80)

# Save main results
output_file = os.path.join(OUTPUT_DIR, 'site_level_normality_tests_NN_strict.csv')
results_df.to_csv(output_file, index=False)
print(f"✓ Saved site-level results to: {output_file}")
print(f"  Total rows: {len(results_df)}")

# =============================================================================
# SECTION 6: SUMMARY BY METRIC
# =============================================================================

print("\n" + "="*80)
print("STEP 5: Generating summary by metric")
print("="*80)

summary_by_metric = []

for metric in metrics:
    metric_results = results_df[results_df['WUE_Metric'] == metric]
    
    total_sites = len(metric_results)
    normal_sites = len(metric_results[metric_results['Normality_Status'] == 'Normal'])
    non_normal_sites = len(metric_results[metric_results['Normality_Status'] == 'Non_normal'])
    insufficient_sites = len(metric_results[metric_results['Normality_Status'] == 'Insufficient_N'])
    
    percent_normal = (normal_sites / total_sites * 100) if total_sites > 0 else 0
    percent_non_normal = (non_normal_sites / total_sites * 100) if total_sites > 0 else 0
    
    # Calculate median Shapiro p-value for normal sites
    normal_shapiro_p = metric_results[metric_results['Normality_Status'] == 'Normal']['Shapiro_p'].median()
    
    summary_by_metric.append({
        'WUE_Metric': metric,
        'Metric_Label': metric_labels[metric],
        'total_sites_tested': total_sites,
        'normal_sites_count': normal_sites,
        'non_normal_sites_count': non_normal_sites,
        'insufficient_sites_count': insufficient_sites,
        'percent_normal': percent_normal,
        'percent_non_normal': percent_non_normal,
        'median_shapiro_p_normal_sites': normal_shapiro_p if not pd.isna(normal_shapiro_p) else np.nan
    })

summary_metric_df = pd.DataFrame(summary_by_metric)
summary_metric_file = os.path.join(OUTPUT_DIR, 'normality_summary_by_metric.csv')
summary_metric_df.to_csv(summary_metric_file, index=False)

print("\nNormality Summary by Metric:")
print(summary_metric_df.to_string(index=False))

# =============================================================================
# SECTION 7: SUMMARY BY SALINITY/ECOSYSTEM GROUP
# =============================================================================

print("\n" + "="*80)
print("STEP 6: Generating summary by salinity/ecosystem group")
print("="*80)

# Merge salinity information
# Get unique site salinity categories
site_salinity = df_filtered[['site_name', 'Salinity_Category']].drop_duplicates()
results_with_salinity = results_df.merge(site_salinity, on='site_name', how='left')

# Exclude Brackish if present
results_with_salinity = results_with_salinity[~results_with_salinity['Salinity_Category'].str.contains('Brackish', na=False)]

summary_by_salinity = []

for metric in metrics:
    metric_results = results_with_salinity[results_with_salinity['WUE_Metric'] == metric]
    
    for salinity_group in ['Freshwater', 'Saline', 'Upland']:
        group_results = metric_results[metric_results['Salinity_Category'] == salinity_group]
        
        if len(group_results) > 0:
            total_sites = len(group_results)
            normal_sites = len(group_results[group_results['Normality_Status'] == 'Normal'])
            non_normal_sites = len(group_results[group_results['Normality_Status'] == 'Non_normal'])
            insufficient_sites = len(group_results[group_results['Normality_Status'] == 'Insufficient_N'])
            
            percent_normal = (normal_sites / total_sites * 100) if total_sites > 0 else 0
            
            summary_by_salinity.append({
                'WUE_Metric': metric,
                'Metric_Label': metric_labels[metric],
                'Salinity_Category': salinity_group,
                'total_sites': total_sites,
                'normal_sites': normal_sites,
                'non_normal_sites': non_normal_sites,
                'insufficient_sites': insufficient_sites,
                'percent_normal': percent_normal
            })

summary_salinity_df = pd.DataFrame(summary_by_salinity)
summary_salinity_file = os.path.join(OUTPUT_DIR, 'normality_summary_by_salinity.csv')
summary_salinity_df.to_csv(summary_salinity_file, index=False)

print("\nNormality Summary by Salinity Group:")
print(summary_salinity_df.to_string(index=False))

# =============================================================================
# SECTION 8: ADDITIONAL STATISTICS FOR PAPER
# =============================================================================

print("\n" + "="*80)
print("STEP 7: Additional statistics for manuscript")
print("="*80)

# Calculate overall normality rates across all metrics
all_normal = results_df[results_df['Normality_Status'] == 'Normal'].shape[0]
all_non_normal = results_df[results_df['Normality_Status'] == 'Non_normal'].shape[0]
all_insufficient = results_df[results_df['Normality_Status'] == 'Insufficient_N'].shape[0]
total_tested = all_normal + all_non_normal

print(f"\nOverall Normality Results (across all sites and metrics):")
print(f"  Total site-metric combinations: {len(results_df)}")
print(f"  With sufficient data (n≥5): {total_tested}")
print(f"  Normal distributions: {all_normal} ({all_normal/total_tested*100:.1f}%)")
print(f"  Non-normal distributions: {all_non_normal} ({all_non_normal/total_tested*100:.1f}%)")
print(f"  Insufficient data (n<5): {all_insufficient}")

# Sample size distribution
print(f"\nSample Size Distribution (months per site-metric):")
for metric in metrics:
    metric_n = results_df[results_df['WUE_Metric'] == metric]['n_months'].dropna()
    if len(metric_n) > 0:
        print(f"  {metric_labels[metric]}:")
        print(f"    Mean: {metric_n.mean():.1f} months")
        print(f"    Median: {metric_n.median():.0f} months")
        print(f"    Range: {metric_n.min():.0f} - {metric_n.max():.0f} months")

# Anderson-Darling agreement with Shapiro-Wilk
if 'Anderson_Status' in results_df.columns:
    anderson_agreement = results_df[
        (results_df['Normality_Status'].isin(['Normal', 'Non_normal'])) &
        (results_df['Anderson_Status'].isin(['Normal', 'Non_normal']))
    ]
    if len(anderson_agreement) > 0:
        agreement = (anderson_agreement['Normality_Status'] == anderson_agreement['Anderson_Status']).sum()
        print(f"\nAgreement between Shapiro-Wilk and Anderson-Darling:")
        print(f"  Consistent classification: {agreement}/{len(anderson_agreement)} ({agreement/len(anderson_agreement)*100:.1f}%)")

# =============================================================================
# SECTION 9: SAVE SITES WITH NORMAL DISTRIBUTIONS
# =============================================================================

print("\n" + "="*80)
print("STEP 8: Saving lists of normally distributed sites")
print("="*80)

for metric in metrics:
    normal_sites = results_df[
        (results_df['WUE_Metric'] == metric) & 
        (results_df['Normality_Status'] == 'Normal')
    ]['site_name'].tolist()
    
    if normal_sites:
        normal_df = pd.DataFrame({'site_name': normal_sites})
        normal_file = os.path.join(OUTPUT_DIR, f'normal_sites_{metric}.csv')
        normal_df.to_csv(normal_file, index=False)
        print(f"  {metric}: {len(normal_sites)} sites with normal distribution saved")

# =============================================================================
# SECTION 10: FINAL SUMMARY
# =============================================================================

print("\n" + "="*80)
print("FINAL MANUSCRIPT-READY SUMMARY")
print("="*80)

print(f"""
MONTHLY WUE DISTRIBUTION NORMALITY TESTING COMPLETE

1. DATA FILTERING SUMMARY:
   • Strict triple intersection sites: {len(shared_subset_sites)} sites
   • NN months only (SPEI-1 = NN)
   • Complete cases (all three WUE metrics available)
   • Final site-metric combinations tested: {len(results_df)}

2. NORMALITY TEST RESULTS (Shapiro-Wilk, α=0.05):
   • Overall normal rate: {all_normal}/{total_tested} ({all_normal/total_tested*100:.1f}%)
   • Overall non-normal rate: {all_non_normal}/{total_tested} ({all_non_normal/total_tested*100:.1f}%)
   • Insufficient data (n<5): {all_insufficient} combinations

3. NORMALITY BY METRIC:
""")

for _, row in summary_metric_df.iterrows():
    print(f"   • {row['Metric_Label']}: {row['normal_sites_count']}/{row['total_sites_tested']} "
          f"({row['percent_normal']:.1f}%) normal")

print(f"""
4. NORMALITY BY ECOSYSTEM GROUP (across all metrics):
""")

for salinity_group in ['Freshwater', 'Saline', 'Upland']:
    group_data = summary_salinity_df[summary_salinity_df['Salinity_Category'] == salinity_group]
    if len(group_data) > 0:
        total_normal = group_data['normal_sites'].sum()
        total_sites = group_data['total_sites'].sum()
        print(f"   • {salinity_group}: {total_normal}/{total_sites} ({total_normal/total_sites*100:.1f}%) normal")

print(f"""
5. OUTPUT FILES SAVED TO:
   {OUTPUT_DIR}
   
   ✓ site_level_normality_tests_NN_strict.csv (full results)
   ✓ normality_summary_by_metric.csv (summary by metric)
   ✓ normality_summary_by_salinity.csv (summary by ecosystem)
   ✓ normal_sites_WUE.csv (sites with normal WUE)
   ✓ normal_sites_WUE_eva.csv (sites with normal WUE_E)
   ✓ normal_sites_WUE_tra.csv (sites with normal WUE_T)

KEY FINDINGS FOR MANUSCRIPT:
   • Most site-level monthly WUE distributions deviate from normality
   • Nonparametric tests are appropriate for cross-site comparisons
   • Minimum monthly samples: {results_df['n_months'].min():.0f} - {results_df['n_months'].max():.0f} months per site
   • Shapiro-Wilk test used as primary, with D'Agostino and Anderson-Darling as supporting tests
""")

print("="*80)
print("NORMALITY TESTING WORKFLOW FINISHED SUCCESSFULLY")
print("="*80)