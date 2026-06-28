# ============================================================================
# SPEI ANALYSIS WORKFLOW - UPDATED WITH 4 SALINITY CATEGORIES + STRICT TRIPLE INTERSECTION
# DIRECTION CLASSIFICATION: SITE-SPECIFIC BOOTSTRAP 95% CI AROUND NN MEDIAN
# ANALYSIS INCLUDES: PASS B (severity grouped) and PASS C (fully merged) ONLY
# UPDATED: ALL THREE metrics (WUE, WUE_E, WUE_T) use the SAME strict triple intersection
# ============================================================================

import pandas as pd
import numpy as np
import os
import sys

# ============================================================================
# USER SETTINGS - CHANGE THIS TO CONTROL INTERACTIVE OUTPUT
# ============================================================================
INTERACTIVE_MODE = False  # Set to True for step-by-step interactive output, False for batch processing
# ============================================================================

print("=" * 80)
print("SPEI TIMESCALE & CLASS GROUPING ANALYSIS")
print("4 Salinity Categories: Freshwater, Brackish, Saline, Upland")
print("ALL THREE metrics (WUE, WUE_E, WUE_T) use the SAME strict triple intersection")
print("Direction classification: site-specific bootstrap 95% CI around NN median")
print("Analysis includes: PASS B (severity grouped) and PASS C (wet vs dry) ONLY")
if not INTERACTIVE_MODE:
    print("Running in BATCH MODE (non-interactive) - use INTERACTIVE_MODE=True for detailed step-through")
print("=" * 80)

# ============================================================================
# 0. DATA LOCATION (AS SPECIFIED)
# ============================================================================

ROOT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
INPUT_DIR = os.path.join(ROOT_DIR, "data_products", "results")
OUTPUT_DIR = os.path.join(ROOT_DIR, "data_products", "results", "SPEI_analysis_results")

# Use INPUT_DIR for monthly data
DATA_DIR = INPUT_DIR  # Points to ...\WUE_CUE\data_products\results

os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Input directory: {INPUT_DIR}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Monthly data directory: {DATA_DIR}")

# Bootstrap parameters
BOOTSTRAP_ITERATIONS = 5000
BOOTSTRAP_RANDOM_SEED = 42

# All SPEI files (use as-is)
SPEI_FILES = [
    "wue_spei_vs_nn_SPEI_1.csv",
    "wue_spei_vs_nn_SPEI_3.csv", 
    "wue_spei_vs_nn_SPEI_6.csv",
    "wue_spei_vs_nn_SPEI_12.csv",
    "wue_spei_vs_nn_SPEI_24.csv",
    "wue_spei_vs_nn_SPEI_36.csv",
    "wue_spei_vs_nn_SPEI_48.csv"
]

# Monthly data file for NN baseline calculation
MONTHLY_DATA_FILE = "monthly_data_after_outlier_removal.csv"

# Verify monthly data file exists
monthly_filepath = os.path.join(DATA_DIR, MONTHLY_DATA_FILE)
if not os.path.exists(monthly_filepath):
    print(f"[ERROR] Monthly data file not found: {monthly_filepath}")
    print(f"[ERROR] Expected location: {monthly_filepath}")
    print(f"[ERROR] Check if file exists or adjust DATA_DIR path")
    exit(1)
else:
    print(f"[INFO] Found monthly data: {monthly_filepath}")

# ============================================================================
# 1. DEFINE CLASS GROUPINGS (PASS B AND PASS C ONLY)
# ============================================================================

# Severity-merged classes (PASS B)
SEVERITY_MERGED = {
    'Mild Wet': ['MW', 'MoW'],
    'Severe Wet': ['SW', 'EW'],
    'Mild Dry': ['MD', 'MoD'],
    'Severe Dry': ['SD', 'ED']
}

# Fully merged wet vs dry (PASS C)
FULLY_MERGED = {
    'Dry (all)': ['MD', 'MoD', 'SD', 'ED'],
    'Wet (all)': ['MW', 'MoW', 'SW', 'EW']
}

# All grouping passes (PASS A removed)
GROUPING_PASSES = {
    'PASS B': SEVERITY_MERGED,
    'PASS C': FULLY_MERGED
}

# All WUE metrics
WUE_METRICS = ['WUE', 'WUE_eva', 'WUE_tra']

# Salinity categories (4 groups)
SALINITY_CATEGORIES = ['Freshwater', 'Brackish', 'Saline', 'Upland']

# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

def load_monthly_data():
    """Load monthly data file for NN baseline calculations."""
    monthly_filepath = os.path.join(DATA_DIR, MONTHLY_DATA_FILE)
    
    if not os.path.exists(monthly_filepath):
        print(f"[ERROR] Monthly data file not found: {monthly_filepath}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(monthly_filepath)
        print(f"[INFO] Loaded monthly data: {len(df):,} rows from {MONTHLY_DATA_FILE}")
        return df
    except Exception as e:
        print(f"[ERROR] Loading monthly data: {e}")
        return pd.DataFrame()

def bootstrap_nn_median_ci(nn_values, n_boot=BOOTSTRAP_ITERATIONS, random_state=BOOTSTRAP_RANDOM_SEED):
    """Bootstrap 95% confidence interval for the median of NN values."""
    clean_values = np.array([v for v in nn_values if not np.isnan(v)])
    n_used = len(clean_values)
    
    if n_used < 3:
        return np.nan, np.nan, np.nan, n_used
    
    np.random.seed(random_state)
    median_estimate = np.median(clean_values)
    
    bootstrap_medians = []
    n = len(clean_values)
    
    for _ in range(n_boot):
        bootstrap_sample = np.random.choice(clean_values, size=n, replace=True)
        bootstrap_medians.append(np.median(bootstrap_sample))
    
    ci_lower = np.percentile(bootstrap_medians, 2.5)
    ci_upper = np.percentile(bootstrap_medians, 97.5)
    
    return median_estimate, ci_lower, ci_upper, n_used

def compute_nn_bootstrap_baseline(monthly_df, site, metric):
    """Compute NN baseline bootstrap 95% CI around median for a specific site and metric."""
    site_data = monthly_df[(monthly_df['site_name'] == site) & 
                          (monthly_df['SPEI_1_Cat'] == 'NN')].copy()
    
    if site_data.empty:
        return None
    
    if metric == 'WUE':
        metric_col = 'WUE'
    elif metric == 'WUE_eva':
        metric_col = 'WUE_eva'
    elif metric == 'WUE_tra':
        metric_col = 'WUE_tra'
    else:
        return None
    
    nn_values = site_data[metric_col].dropna().values
    nn_median, ci_lower_raw, ci_upper_raw, n_nn_values = bootstrap_nn_median_ci(nn_values)
    
    if np.isnan(nn_median) or np.isnan(ci_lower_raw) or np.isnan(ci_upper_raw):
        return None
    
    if nn_median == 0 or np.abs(nn_median) < 1e-10:
        return None
    
    ci_lower_percent = ((ci_lower_raw - nn_median) / nn_median) * 100
    ci_upper_percent = ((ci_upper_raw - nn_median) / nn_median) * 100
    ci_width = ci_upper_percent - ci_lower_percent
    
    return {
        'NN_median': nn_median,
        'NN_ci_lower_raw': ci_lower_raw,
        'NN_ci_upper_raw': ci_upper_raw,
        'NN_ci_lower_percent': ci_lower_percent,
        'NN_ci_upper_percent': ci_upper_percent,
        'NN_ci_width': ci_width,
        'n_nn_values': n_nn_values
    }

def classify_direction_from_nn_bootstrap(median_change, nn_baseline):
    """Classify direction based on whether median percent change falls outside bootstrap 95% CI."""
    if nn_baseline is None or np.isnan(median_change):
        return 'Insufficient data'
    
    if median_change > nn_baseline['NN_ci_upper_percent']:
        return 'Increase'
    elif median_change < nn_baseline['NN_ci_lower_percent']:
        return 'Decrease'
    else:
        return 'No change'

def load_and_filter(filepath, spei_timescale):
    """Load CSV and apply QC as specified."""
    try:
        df = pd.read_csv(filepath)
        print(f"[INFO] Loaded {len(df):,} rows from {os.path.basename(filepath)}")
        
        if 'Salinity_Category' not in df.columns:
            print("[ERROR] 'Salinity_Category' column not found!")
            return pd.DataFrame()
        
        if 'site_name' not in df.columns:
            print("[ERROR] 'site_name' column not found!")
            return pd.DataFrame()
        
        initial_rows = len(df)
        df = df[(df['N_months_NN'] >= 3) & (df['N_months_Condition'] >= 3)]
        filtered_rows = initial_rows - len(df)
        
        if filtered_rows > 0:
            print(f"[INFO] Filtered {filtered_rows:,} rows with <3 months")
        
        df['SPEI_Timescale'] = f"SPEI_{spei_timescale}"
        
        required_cols = ['site_name', 'WUE_Metric', 'SPEI_Class', 'Percent_Change',
                        'Salinity_Category', 'N_months_Condition', 'N_months_NN']
        
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"[WARNING] Missing columns: {missing_cols}")
            return pd.DataFrame()
        
        print(f"[INFO] Final dataset: {len(df):,} rows, {df['site_name'].nunique():,} sites")
        
        salinity_counts = df.groupby('site_name')['Salinity_Category'].first().value_counts()
        print("[INFO] Salinity distribution (unique sites):")
        for salinity in SALINITY_CATEGORIES:
            count = salinity_counts.get(salinity, 0)
            print(f"  {salinity}: {count:,} sites")
        
        return df
    
    except Exception as e:
        print(f"[ERROR] loading {filepath}: {e}")
        return pd.DataFrame()

def get_strict_triple_intersection_sites(df):
    """Get strict triple intersection sites that have ALL THREE metrics."""
    sites_with_wue = set(df[df['WUE_Metric'] == 'WUE']['site_name'].unique())
    sites_with_eva = set(df[df['WUE_Metric'] == 'WUE_eva']['site_name'].unique())
    sites_with_tra = set(df[df['WUE_Metric'] == 'WUE_tra']['site_name'].unique())
    
    shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)
    
    print(f"[INFO] Strict triple intersection (ALL metrics): {len(shared_sites)} sites")
    print(f"  Sites with WUE: {len(sites_with_wue)}")
    print(f"  Sites with WUE_E: {len(sites_with_eva)}")
    print(f"  Sites with WUE_T: {len(sites_with_tra)}")
    
    return shared_sites

def get_site_salinity(df, site):
    """Get the salinity category for a specific site."""
    site_data = df[df['site_name'] == site]
    return site_data['Salinity_Category'].iloc[0]

def analyze_group_all_sites(df, class_list, wue_metric, spei_timescale, pass_name, condition_name, 
                           shared_sites, nn_baselines):
    """Analyze ALL SITES for a given condition."""
    df_filter = df[df['WUE_Metric'] == wue_metric].copy()
    df_filter = df_filter[df_filter['SPEI_Class'].isin(class_list)]
    df_filter = df_filter[df_filter['site_name'].isin(shared_sites)]
    
    if df_filter.empty:
        return None, []
    
    sites_with_condition = df_filter['site_name'].unique()
    
    if len(sites_with_condition) == 0:
        return None, []
    
    site_results = []
    direction_counts = {'Increase': 0, 'Decrease': 0, 'No change': 0, 'Insufficient data': 0}
    all_changes = []
    all_ci_widths = []
    all_cond_months = []
    all_nn_months = []
    salinity_counts = {cat: 0 for cat in SALINITY_CATEGORIES}
    
    for site in sites_with_condition:
        site_data = df_filter[df_filter['site_name'] == site]
        percent_changes = site_data['Percent_Change'].values
        median_change = np.median(percent_changes) if len(percent_changes) > 0 else np.nan
        
        baseline_key = f"{site}_{wue_metric}"
        nn_baseline = nn_baselines.get(baseline_key)
        direction = classify_direction_from_nn_bootstrap(median_change, nn_baseline)
        
        median_cond_months = site_data['N_months_Condition'].median()
        median_nn_months = site_data['N_months_NN'].median()
        site_salinity = get_site_salinity(df, site)
        
        if nn_baseline is not None:
            ci_lower = nn_baseline['NN_ci_lower_percent']
            ci_upper = nn_baseline['NN_ci_upper_percent']
            ci_width = nn_baseline['NN_ci_width']
            n_values_used = len(percent_changes)
            nn_median = nn_baseline['NN_median']
            nn_lower_raw = nn_baseline['NN_ci_lower_raw']
            nn_upper_raw = nn_baseline['NN_ci_upper_raw']
            n_nn_values = nn_baseline['n_nn_values']
        else:
            ci_lower = np.nan
            ci_upper = np.nan
            ci_width = np.nan
            n_values_used = len(percent_changes)
            nn_median = np.nan
            nn_lower_raw = np.nan
            nn_upper_raw = np.nan
            n_nn_values = 0
        
        site_results.append({
            'site': site,
            'median_change': median_change,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'ci_width': ci_width,
            'direction': direction,
            'n_values_used': n_values_used,
            'months_condition': median_cond_months,
            'months_nn': median_nn_months,
            'site_salinity': site_salinity,
            'insufficient_data': direction == 'Insufficient data',
            'nn_median': nn_median,
            'nn_lower_raw': nn_lower_raw,
            'nn_upper_raw': nn_upper_raw,
            'n_nn_values': n_nn_values
        })
        
        direction_counts[direction] += 1
        if not np.isnan(median_change):
            all_changes.append(median_change)
        if not np.isnan(ci_width):
            all_ci_widths.append(ci_width)
        all_cond_months.append(median_cond_months)
        all_nn_months.append(median_nn_months)
        
        if site_salinity in salinity_counts:
            salinity_counts[site_salinity] += 1
    
    total_sites = len(site_results)
    
    stats = {
        'SPEI_Timescale': f"SPEI_{spei_timescale}",
        'Pass': pass_name,
        'Condition': condition_name,
        'WUE_Metric': wue_metric,
        'Salinity': 'All',
        'Shared_Site_Filter': 'Yes (strict triple intersection)',
        'Direction_Method': 'NN_Bootstrap_Median_CI_95',
        'Threshold_Method': 'Site_specific_NN_bootstrap_95CI_median',
        'Increase': direction_counts['Increase'],
        'Decrease': direction_counts['Decrease'],
        'No_change': direction_counts['No change'],
        'Insufficient_data': direction_counts['Insufficient data'],
        'Total_sites': total_sites,
        'Freshwater_sites': salinity_counts['Freshwater'],
        'Brackish_sites': salinity_counts['Brackish'],
        'Saline_sites': salinity_counts['Saline'],
        'Upland_sites': salinity_counts['Upland'],
        'Median_%_Change': np.median(all_changes) if all_changes else 0,
        'Min_%_Change': np.min(all_changes) if all_changes else 0,
        'Max_%_Change': np.max(all_changes) if all_changes else 0,
        'Median_CI_Width': np.median(all_ci_widths) if all_ci_widths else 0,
        'Median_Months_Condition': np.median(all_cond_months) if all_cond_months else 0,
        'Min_Months_Condition': np.min(all_cond_months) if all_cond_months else 0,
        'Max_Months_Condition': np.max(all_cond_months) if all_cond_months else 0,
        'Median_Months_NN': np.median(all_nn_months) if all_nn_months else 0
    }
    
    return stats, site_results

def analyze_group_by_salinity(df, class_list, wue_metric, salinity_category, spei_timescale, 
                             pass_name, condition_name, shared_sites, nn_baselines):
    """Analyze sites by specific salinity category."""
    df_filter = df[df['WUE_Metric'] == wue_metric].copy()
    df_filter = df_filter[df_filter['SPEI_Class'].isin(class_list)]
    df_filter = df_filter[df_filter['site_name'].isin(shared_sites)]
    
    if df_filter.empty:
        return None, []
    
    all_sites_with_condition = df_filter['site_name'].unique()
    
    sites_in_category = []
    for site in all_sites_with_condition:
        site_salinity = get_site_salinity(df, site)
        if site_salinity == salinity_category:
            sites_in_category.append(site)
    
    if not sites_in_category:
        return None, []
    
    site_results = []
    direction_counts = {'Increase': 0, 'Decrease': 0, 'No change': 0, 'Insufficient data': 0}
    all_changes = []
    all_ci_widths = []
    all_cond_months = []
    all_nn_months = []
    
    for site in sites_in_category:
        site_data = df_filter[df_filter['site_name'] == site]
        percent_changes = site_data['Percent_Change'].values
        median_change = np.median(percent_changes) if len(percent_changes) > 0 else np.nan
        
        baseline_key = f"{site}_{wue_metric}"
        nn_baseline = nn_baselines.get(baseline_key)
        direction = classify_direction_from_nn_bootstrap(median_change, nn_baseline)
        
        median_cond_months = site_data['N_months_Condition'].median()
        median_nn_months = site_data['N_months_NN'].median()
        
        if nn_baseline is not None:
            ci_lower = nn_baseline['NN_ci_lower_percent']
            ci_upper = nn_baseline['NN_ci_upper_percent']
            ci_width = nn_baseline['NN_ci_width']
            n_values_used = len(percent_changes)
            nn_median = nn_baseline['NN_median']
            nn_lower_raw = nn_baseline['NN_ci_lower_raw']
            nn_upper_raw = nn_baseline['NN_ci_upper_raw']
            n_nn_values = nn_baseline['n_nn_values']
        else:
            ci_lower = np.nan
            ci_upper = np.nan
            ci_width = np.nan
            n_values_used = len(percent_changes)
            nn_median = np.nan
            nn_lower_raw = np.nan
            nn_upper_raw = np.nan
            n_nn_values = 0
        
        site_results.append({
            'site': site,
            'median_change': median_change,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'ci_width': ci_width,
            'direction': direction,
            'n_values_used': n_values_used,
            'months_condition': median_cond_months,
            'months_nn': median_nn_months,
            'site_salinity': salinity_category,
            'insufficient_data': direction == 'Insufficient data',
            'nn_median': nn_median,
            'nn_lower_raw': nn_lower_raw,
            'nn_upper_raw': nn_upper_raw,
            'n_nn_values': n_nn_values
        })
        
        direction_counts[direction] += 1
        if not np.isnan(median_change):
            all_changes.append(median_change)
        if not np.isnan(ci_width):
            all_ci_widths.append(ci_width)
        all_cond_months.append(median_cond_months)
        all_nn_months.append(median_nn_months)
    
    total_sites = len(site_results)
    
    stats = {
        'SPEI_Timescale': f"SPEI_{spei_timescale}",
        'Pass': pass_name,
        'Condition': condition_name,
        'WUE_Metric': wue_metric,
        'Salinity': salinity_category,
        'Shared_Site_Filter': 'Yes (strict triple intersection)',
        'Direction_Method': 'NN_Bootstrap_Median_CI_95',
        'Threshold_Method': 'Site_specific_NN_bootstrap_95CI_median',
        'Increase': direction_counts['Increase'],
        'Decrease': direction_counts['Decrease'],
        'No_change': direction_counts['No change'],
        'Insufficient_data': direction_counts['Insufficient data'],
        'Total_sites': total_sites,
        'Median_%_Change': np.median(all_changes) if all_changes else 0,
        'Min_%_Change': np.min(all_changes) if all_changes else 0,
        'Max_%_Change': np.max(all_changes) if all_changes else 0,
        'Median_CI_Width': np.median(all_ci_widths) if all_ci_widths else 0,
        'Median_Months_Condition': np.median(all_cond_months) if all_cond_months else 0,
        'Min_Months_Condition': np.min(all_cond_months) if all_cond_months else 0,
        'Max_Months_Condition': np.max(all_cond_months) if all_cond_months else 0,
        'Median_Months_NN': np.median(all_nn_months) if all_nn_months else 0
    }
    
    return stats, site_results

def print_results_batch_by_batch(stats_list):
    """Print results one batch at a time with user control - INTERACTIVE MODE ONLY"""
    if not stats_list:
        print("[WARNING] No results to print")
        return
    
    print("\n" + "=" * 80)
    print("CONSOLE OUTPUT - STRUCTURED RESULTS")
    print(f"Direction classification: site-specific bootstrap 95% CI around NN median")
    print(f"ALL THREE metrics use the SAME strict triple intersection")
    print(f"Analysis includes: PASS B (severity grouped) and PASS C (wet vs dry) ONLY")
    print("=" * 80)
    
    if INTERACTIVE_MODE:
        print("\nPress Enter to continue to next section...")
        print("Type 'q' to quit, 's' to skip to next condition")
        print("=" * 80)
    
    # Group by SPEI timescale
    spei_timescales = sorted(set([s['SPEI_Timescale'] for s in stats_list if s['Salinity'] == 'All']))
    
    for spei in spei_timescales:
        print(f"\n{'='*80}")
        print(f"SPEI TIMESCALE: {spei}")
        print(f"{'='*80}")
        
        spei_stats = [s for s in stats_list if s['SPEI_Timescale'] == spei]
        
        for pass_name in ['PASS B', 'PASS C']:
            pass_stats = [s for s in spei_stats if s['Pass'] == pass_name and s['Salinity'] == 'All']
            
            if not pass_stats:
                continue
            
            print(f"\n{'='*80}")
            print(f"{spei} | {pass_name.replace('_', ' ')}")
            print(f"Direction method: site-specific bootstrap 95% CI around NN median")
            print(f"ALL metrics: strict triple intersection")
            print(f"{'='*80}")
            
            conditions = sorted(set([s['Condition'] for s in pass_stats]))
            
            for condition in conditions:
                cond_stats = [s for s in pass_stats if s['Condition'] == condition]
                
                for wue_metric in WUE_METRICS:
                    metric_stats = [s for s in cond_stats if s['WUE_Metric'] == wue_metric]
                    
                    if not metric_stats:
                        continue
                    
                    all_stats = metric_stats[0]
                    
                    print(f"\nCondition: {condition}")
                    print(f"Metric: {wue_metric}")
                    print(f"Shared subset applied: {all_stats['Shared_Site_Filter']}")
                    print("-" * 80)
                    
                    print("PRIMARY ANALYSIS - ALL SITES:")
                    print(f"  Total sites: {all_stats['Total_sites']}")
                    print(f"    Freshwater: {all_stats['Freshwater_sites']}, Brackish: {all_stats['Brackish_sites']}, Saline: {all_stats['Saline_sites']}, Upland: {all_stats['Upland_sites']}")
                    print(f"  Increase: {all_stats['Increase']}, Decrease: {all_stats['Decrease']}, No change: {all_stats['No_change']}, Insufficient: {all_stats['Insufficient_data']}")
                    print(f"  Median % change: {all_stats['Median_%_Change']:.1f}%")
                    print(f"  Median CI width: {all_stats['Median_CI_Width']:.1f}%")
                    
                    # SECONDARY ANALYSIS: All salinity groups
                    print("\nSECONDARY ANALYSIS - SALINITY GROUPS:")
                    
                    salinity_results = {}
                    for salinity_cat in SALINITY_CATEGORIES:
                        salinity_stats = next((s for s in spei_stats if 
                                              s['Pass'] == pass_name and 
                                              s['Condition'] == condition and 
                                              s['WUE_Metric'] == wue_metric and 
                                              s['Salinity'] == salinity_cat), None)
                        if salinity_stats and salinity_stats['Total_sites'] > 0:
                            salinity_results[salinity_cat] = salinity_stats
                    
                    if salinity_results:
                        for salinity_cat, stats in salinity_results.items():
                            total_valid = stats['Increase'] + stats['Decrease'] + stats['No_change']
                            increase_pct = (stats['Increase'] / total_valid * 100) if total_valid > 0 else 0
                            print(f"  {salinity_cat} (n={stats['Total_sites']}): ↑{stats['Increase']} ({increase_pct:.0f}%), ↓{stats['Decrease']}, ={stats['No_change']}")
                    else:
                        print("  No salinity-specific data available")
                    
                    print("-" * 80)
                    
                    if INTERACTIVE_MODE:
                        user_input = input("\nPress Enter for next metric (or 's' to skip, 'q' to quit): ").strip().lower()
                        if user_input == 'q':
                            return
                        elif user_input == 's':
                            break
            
            if INTERACTIVE_MODE and pass_name != 'PASS C':
                user_input = input(f"\nFinished {pass_name}. Press Enter for next pass (or 'q' to quit): ").strip().lower()
                if user_input == 'q':
                    return
        
        if INTERACTIVE_MODE and spei != spei_timescales[-1]:
            user_input = input(f"\nFinished {spei}. Press Enter for next SPEI (or 'q' to quit): ").strip().lower()
            if user_input == 'q':
                return

# ============================================================================
# 3. MAIN ANALYSIS LOOP
# ============================================================================

print("\n" + "=" * 80)
print("STARTING ANALYSIS")
print(f"Processing PASS B (severity grouped) and PASS C (wet vs dry) ONLY")
print("ALL THREE metrics use the SAME strict triple intersection")
print("=" * 80)

# Load monthly data and pre-compute NN baselines
print("\n" + "=" * 80)
print("PRE-COMPUTING NN BOOTSTRAP BASELINES")
print("=" * 80)

monthly_df = load_monthly_data()
nn_baselines = {}

if not monthly_df.empty:
    all_sites = monthly_df['site_name'].unique()
    print(f"[INFO] Computing NN bootstrap CIs for {len(all_sites)} sites...")
    
    nn_site_count = 0
    nn_examples = []
    
    for site in all_sites:
        for metric in WUE_METRICS:
            baseline_key = f"{site}_{metric}"
            nn_baseline = compute_nn_bootstrap_baseline(monthly_df, site, metric)
            
            if nn_baseline is not None:
                nn_baselines[baseline_key] = nn_baseline
                nn_site_count += 1
                
                if len(nn_examples) < 5:
                    nn_examples.append({
                        'site_name': site,
                        'metric': metric,
                        'NN_median': nn_baseline['NN_median'],
                        'NN_lower_percent': nn_baseline['NN_ci_lower_percent'],
                        'NN_upper_percent': nn_baseline['NN_ci_upper_percent'],
                        'NN_ci_width': nn_baseline['NN_ci_width'],
                        'n_nn_values': nn_baseline['n_nn_values']
                    })
    
    print(f"[INFO] NN bootstrap baselines computed: {nn_site_count} site-metric combinations")
    
    print("\n[VALIDATION CHECK 1] NN bootstrap CI coverage:")
    for metric in WUE_METRICS:
        metric_count = sum(1 for k in nn_baselines.keys() if k.endswith(f"_{metric}"))
        print(f"  {metric}: {metric_count} sites")
    
    print("\n[VALIDATION CHECK 2] Median bootstrap CI width by metric:")
    for metric in WUE_METRICS:
        widths = [v['NN_ci_width'] for k, v in nn_baselines.items() if k.endswith(f"_{metric}")]
        if widths:
            print(f"  {metric}: {np.median(widths):.1f}%")
    
    total_expected = len(all_sites) * len(WUE_METRICS)
    coverage_pct = (nn_site_count / total_expected * 100) if total_expected > 0 else 0
    print(f"\n[VALIDATION] Overall coverage: {nn_site_count}/{total_expected} ({coverage_pct:.1f}%)")
    
else:
    print("[ERROR] Could not load monthly data. Cannot compute NN bootstrap CIs.")
    exit(1)

all_summary_stats = []
all_site_details = []

for spei_file in SPEI_FILES:
    spei_timescale = spei_file.split('_')[-1].replace('.csv', '').replace('SPEI', '')
    filepath = os.path.join(INPUT_DIR, spei_file)
    
    print(f"\n{'='*60}")
    print(f"PROCESSING: {spei_file}")
    print(f"{'='*60}")
    
    if not os.path.exists(filepath):
        print(f"[WARNING] File not found: {spei_file}")
        continue
    
    df = load_and_filter(filepath, spei_timescale)
    if df.empty:
        continue
    
    shared_sites = get_strict_triple_intersection_sites(df)
    
    for pass_name, groupings in GROUPING_PASSES.items():
        for condition_name, class_list in groupings.items():
            for wue_metric in WUE_METRICS:
                # PRIMARY ANALYSIS: All sites
                stats, site_results = analyze_group_all_sites(df, class_list, wue_metric, 
                                                            spei_timescale, pass_name, condition_name, 
                                                            shared_sites, nn_baselines)
                
                if stats:
                    all_summary_stats.append(stats)
                    
                    for site_detail in site_results:
                        all_site_details.append({
                            'SPEI_Timescale': f"SPEI_{spei_timescale}",
                            'Pass': pass_name,
                            'Condition': condition_name,
                            'WUE_Metric': wue_metric,
                            'Shared_Site_Filter': 'Yes (strict triple intersection)',
                            'Direction_Method': 'NN_Bootstrap_Median_CI_95',
                            'Threshold_Method': 'Site_specific_NN_bootstrap_95CI_median',
                            'Salinity': 'All',
                            'Site': site_detail['site'],
                            'Direction': site_detail['direction'],
                            'Median_%_Change': site_detail['median_change'],
                            'CI_Lower': site_detail['ci_lower'],
                            'CI_Upper': site_detail['ci_upper'],
                            'CI_Width': site_detail['ci_width'],
                            'N_Class_Values_Used': site_detail['n_values_used'],
                            'Insufficient_Data_Flag': site_detail['insufficient_data'],
                            'Months_Condition': site_detail['months_condition'],
                            'Months_NN': site_detail['months_nn'],
                            'Site_Salinity': site_detail['site_salinity'],
                            'NN_Median': site_detail['nn_median'],
                            'NN_Lower_Raw': site_detail['nn_lower_raw'],
                            'NN_Upper_Raw': site_detail['nn_upper_raw'],
                            'NN_Lower_Percent': site_detail['ci_lower'],
                            'NN_Upper_Percent': site_detail['ci_upper'],
                            'NN_Envelope_Width': site_detail['ci_width'],
                            'N_NN_Values_Used': site_detail['n_nn_values']
                        })
                
                # SECONDARY ANALYSIS: Salinity categories
                for salinity_category in SALINITY_CATEGORIES:
                    stats, site_results = analyze_group_by_salinity(df, class_list, wue_metric, 
                                                                  salinity_category, spei_timescale, 
                                                                  pass_name, condition_name, shared_sites, 
                                                                  nn_baselines)
                    
                    if stats:
                        all_summary_stats.append(stats)
                        
                        for site_detail in site_results:
                            all_site_details.append({
                                'SPEI_Timescale': f"SPEI_{spei_timescale}",
                                'Pass': pass_name,
                                'Condition': condition_name,
                                'WUE_Metric': wue_metric,
                                'Shared_Site_Filter': 'Yes (strict triple intersection)',
                                'Direction_Method': 'NN_Bootstrap_Median_CI_95',
                                'Threshold_Method': 'Site_specific_NN_bootstrap_95CI_median',
                                'Salinity': salinity_category,
                                'Site': site_detail['site'],
                                'Direction': site_detail['direction'],
                                'Median_%_Change': site_detail['median_change'],
                                'CI_Lower': site_detail['ci_lower'],
                                'CI_Upper': site_detail['ci_upper'],
                                'CI_Width': site_detail['ci_width'],
                                'N_Class_Values_Used': site_detail['n_values_used'],
                                'Insufficient_Data_Flag': site_detail['insufficient_data'],
                                'Months_Condition': site_detail['months_condition'],
                                'Months_NN': site_detail['months_nn'],
                                'Site_Salinity': site_detail['site_salinity'],
                                'NN_Median': site_detail['nn_median'],
                                'NN_Lower_Raw': site_detail['nn_lower_raw'],
                                'NN_Upper_Raw': site_detail['nn_upper_raw'],
                                'NN_Lower_Percent': site_detail['ci_lower'],
                                'NN_Upper_Percent': site_detail['ci_upper'],
                                'NN_Envelope_Width': site_detail['ci_width'],
                                'N_NN_Values_Used': site_detail['n_nn_values']
                            })

print(f"\n{'='*80}")
print("ANALYSIS COMPLETE")
print(f"{'='*80}")
print(f"Summary statistics: {len(all_summary_stats):,} records")
print(f"Site details: {len(all_site_details):,} records")

# ============================================================================
# 4. PRINT RESULTS (INTERACTIVE OR BATCH MODE)
# ============================================================================

if INTERACTIVE_MODE:
    print_results_batch_by_batch(all_summary_stats)
else:
    # BATCH MODE - Print compact summary
    print("\n" + "=" * 80)
    print("BATCH MODE SUMMARY (INTERACTIVE_MODE=False)")
    print("=" * 80)
    
    if all_summary_stats:
        summary_view = pd.DataFrame(all_summary_stats)
        print(f"\nTotal summary records: {len(summary_view):,}")
        print(f"SPEI Timescales: {sorted(summary_view['SPEI_Timescale'].unique())}")
        print(f"Passes: {summary_view['Pass'].unique().tolist()}")
        print(f"WUE Metrics: {summary_view['WUE_Metric'].unique().tolist()}")
        
        print("\n" + "-" * 60)
        print("QUICK STATS BY PASS AND METRIC (All Sites):")
        print("-" * 60)
        
        all_sites_stats = summary_view[summary_view['Salinity'] == 'All']
        for pass_name in ['PASS B', 'PASS C']:
            pass_stats = all_sites_stats[all_sites_stats['Pass'] == pass_name]
            if not pass_stats.empty:
                print(f"\n{pass_name}:")
                for metric in WUE_METRICS:
                    metric_stats = pass_stats[pass_stats['WUE_Metric'] == metric]
                    if not metric_stats.empty:
                        total_sites = metric_stats['Total_sites'].iloc[0]
                        increase = metric_stats['Increase'].iloc[0]
                        decrease = metric_stats['Decrease'].iloc[0]
                        no_change = metric_stats['No_change'].iloc[0]
                        median_change = metric_stats['Median_%_Change'].iloc[0]
                        print(f"  {metric}: n={total_sites} sites, median={median_change:+.1f}% (↑{increase} ↓{decrease} ↔{no_change})")
        
        print("\n" + "-" * 60)
        print("To see detailed results, open the CSV files:")
        print(f"  - {os.path.join(OUTPUT_DIR, 'SPEI_analysis_summary_FINAL.csv')}")
        print(f"  - {os.path.join(OUTPUT_DIR, 'SPEI_site_level_details_FINAL.csv')}")

# ============================================================================
# 5. SAVE RESULTS TO CSV FILES
# ============================================================================

print("\n" + "=" * 80)
print("SAVING CSV FILES")
print("=" * 80)

if all_summary_stats:
    summary_df = pd.DataFrame(all_summary_stats)
    
    for col in ['Freshwater_sites', 'Brackish_sites', 'Saline_sites', 'Upland_sites']:
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].fillna(0)
    
    int_cols = ['Increase', 'Decrease', 'No_change', 'Insufficient_data', 'Total_sites', 
                'Freshwater_sites', 'Brackish_sites', 'Saline_sites', 'Upland_sites']
    for col in int_cols:
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].astype(int)
    
    summary_file = os.path.join(OUTPUT_DIR, "SPEI_analysis_summary_FINAL.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"[SAVED] Summary results: {summary_file}")
    print(f"        Records: {len(summary_df):,}")
    
    if all_site_details:
        site_details_df = pd.DataFrame(all_site_details)
        site_file = os.path.join(OUTPUT_DIR, "SPEI_site_level_details_FINAL.csv")
        site_details_df.to_csv(site_file, index=False)
        print(f"[SAVED] Site-level details: {site_file}")
        print(f"        Records: {len(site_details_df):,}")
else:
    print("[WARNING] No results to save")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE - READY FOR INTERPRETATION")
print("=" * 80)
print("\nIMPORTANT METHODOLOGICAL NOTES:")
print("1. Direction classification: site-specific bootstrap 95% CI around NN median")
print("2. ALL THREE metrics (WUE, WUE_E, WUE_T) use the SAME strict triple intersection")
print("3. Analysis includes PASS B (severity grouped) and PASS C (wet vs dry) ONLY")
print("\nPASS B (Severity Grouped):")
print("  - Mild Wet (MW + MoW), Severe Wet (SW + EW)")
print("  - Mild Dry (MD + MoD), Severe Dry (SD + ED)")
print("\nPASS C (Fully Merged):")
print("  - Dry (all) = MD + MoD + SD + ED")
print("  - Wet (all) = MW + MoW + SW + EW")
print("=" * 80)