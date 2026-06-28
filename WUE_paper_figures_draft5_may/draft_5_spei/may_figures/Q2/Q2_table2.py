# -*- coding: utf-8 -*-
"""
TABLE S4: Site-level median WUE responses under hydroclimatic anomalies
Supports Figure 4 (magnitude/direction of median responses)

Structure:
- Timescale | Hydroclimatic condition | Metric | N | Median % change | IQR (%) | Wilcoxon p

Formatting:
- Median % change: 2 significant digits (e.g., -4.1, +3.1)
- IQR: 2 significant digits with dash (e.g., -10.3-5.3)
- N: whole number
- Wilcoxon p: 2 significant digits (e.g., 0.18, 0.020)
- Severe/Extreme conditions: labeled as "Severe/Extreme Dry" and "Severe/Extreme Wet"

@author: WUE Analysis Pipeline
Date: 2026-05-18
"""

import pandas as pd
import numpy as np
import os
from scipy.stats import wilcoxon

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results\SPEI_site_level_details_FINAL.csv"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\SPEI_analysis_results"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "Table_S4_WUE_responses_under_hydroclimatic_anomalies.csv")

# Metrics
METRICS = ['WUE', 'WUE_tra']
METRIC_DISPLAY = {
    'WUE': 'WUE_ET',
    'WUE_tra': 'WUE_T'
}

# Timescales
TIMESCALES = ['SPEI_6', 'SPEI_48']

# Conditions for each Pass
PASS_C_CONDITIONS = ['Dry (all)', 'Wet (all)']
PASS_B_CONDITIONS = ['Severe Dry', 'Severe Wet']


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def format_percent_2sig(value):
    """Format percent change with 2 significant digits, include + sign for positive"""
    if pd.isna(value):
        return '--'
    # Round to 2 significant digits
    if value == 0:
        return '0.0'
    elif abs(value) >= 10:
        return f"{value:+.0f}"
    elif abs(value) >= 1:
        return f"{value:+.1f}"
    else:
        return f"{value:+.2f}".rstrip('0').rstrip('.')


def format_iqr_2sig_dash(q1, q3):
    """Format IQR range with 2 significant digits using dash (e.g., -10.3-5.3)"""
    if pd.isna(q1) or pd.isna(q3):
        return '--'
    
    def sig2(x):
        if abs(x) >= 10:
            return f"{x:.0f}"
        elif abs(x) >= 1:
            return f"{x:.1f}"
        else:
            return f"{x:.2f}".rstrip('0').rstrip('.')
    
    return f"{sig2(q1)}-{sig2(q3)}"


def format_pvalue_2sig(p_val):
    """Format p-value with 2 significant digits"""
    if pd.isna(p_val):
        return '--'
    if p_val < 0.001:
        return '<0.001'
    elif p_val < 0.01:
        return f"{p_val:.3f}"
    else:
        return f"{p_val:.2f}"


def compute_median_iqr(values):
    """Compute median and IQR (Q1, Q3)"""
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    median = np.median(values)
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    return median, q1, q3


def compute_wilcoxon_p(values, alternative):
    """
    Compute Wilcoxon signed-rank test p-value
    alternative: 'less' for dry (testing decrease), 'greater' for wet (testing increase)
    """
    if len(values) < 3:
        return np.nan
    
    # Check if all values are zero
    if np.all(values == 0):
        return np.nan
    
    try:
        result = wilcoxon(values, alternative=alternative, zero_method='wilcox')
        return result.pvalue
    except Exception:
        return np.nan


def format_condition_display(condition, pass_name):
    """Format condition for display with Severe/Extreme label"""
    if pass_name == 'PASS B':
        if condition == 'Severe Dry':
            return 'Severe/Extreme Dry'
        elif condition == 'Severe Wet':
            return 'Severe/Extreme Wet'
    else:
        return condition.replace(' (all)', '')
    return condition


# =============================================================================
# MAIN TABLE GENERATION
# =============================================================================

def generate_table_s4():
    """Generate Table S4 from SPEI_site_level_details_FINAL.csv"""
    
    print("=" * 100)
    print("TABLE S4: Site-level median WUE responses under hydroclimatic anomalies")
    print("=" * 100)
    
    # Load data
    print(f"\nLoading data from: {INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df):,} rows")
    print(f"Columns available: {df.columns.tolist()}")
    
    # Verify required columns exist
    required_cols = ['Pass', 'Salinity', 'SPEI_Timescale', 'Condition', 'WUE_Metric', 
                     'Site', 'Median_%_Change']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"⚠️ WARNING: Missing columns: {missing_cols}")
        print("Available columns:", df.columns.tolist())
        return None
    
    results = []
    
    print("\n" + "-" * 80)
    print("PROCESSING ALL CONDITIONS (PASS C)")
    print("-" * 80)
    
    # Process PASS C (All conditions)
    for timescale in TIMESCALES:
        for condition in PASS_C_CONDITIONS:
            for metric in METRICS:
                process_combination(df, results, 'PASS C', timescale, condition, metric)
    
    print("\n" + "-" * 80)
    print("PROCESSING SEVERE/EXTREME CONDITIONS (PASS B)")
    print("-" * 80)
    
    # Process PASS B (Severe/Extreme conditions)
    for timescale in TIMESCALES:
        for condition in PASS_B_CONDITIONS:
            for metric in METRICS:
                process_combination(df, results, 'PASS B', timescale, condition, metric)
    
    # Create DataFrame
    table_df = pd.DataFrame(results)
    
    # Reorder columns
    column_order = [
        'Timescale', 'Hydroclimatic condition', 'Metric', 'N sites', 
        'Median % change', 'IQR (%)', 'Wilcoxon p'
    ]
    table_df = table_df[column_order]
    
    # Save to CSV
    table_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✅ Table S4 saved to: {OUTPUT_FILE}")
    
    # Print table
    print_table_s4(table_df)
    
    return table_df


def process_combination(df, results, pass_name, timescale, condition, metric):
    """
    Process one combination of pass, timescale, condition, metric
    Uses strict double intersection (sites with both WUE_ET and WUE_T)
    """
    
    # Filter data
    subset = df[
        (df['Pass'] == pass_name) &
        (df['Salinity'] == 'All') &
        (df['SPEI_Timescale'] == timescale) &
        (df['Condition'] == condition) &
        (df['WUE_Metric'].isin(METRICS))
    ].copy()
    
    if len(subset) == 0:
        return
    
    # Get sites for each metric (for double intersection)
    wue_sites = set(subset[subset['WUE_Metric'] == 'WUE']['Site'].unique())
    tra_sites = set(subset[subset['WUE_Metric'] == 'WUE_tra']['Site'].unique())
    
    # Strict double intersection
    shared_sites = wue_sites.intersection(tra_sites)
    
    # Filter to current metric and shared sites
    metric_data = subset[
        (subset['WUE_Metric'] == metric) &
        (subset['Site'].isin(shared_sites))
    ].copy()
    
    if len(metric_data) == 0:
        return
    
    # Aggregate to site level (median per site)
    site_level = metric_data.groupby('Site')['Median_%_Change'].median().dropna()
    values = site_level.values
    n_sites = len(values)
    
    if n_sites == 0:
        return
    
    # Compute statistics
    median_val, q1, q3 = compute_median_iqr(values)
    
    # Determine alternative for Wilcoxon
    is_dry = 'Dry' in condition
    alternative = 'less' if is_dry else 'greater'
    
    # Compute Wilcoxon p-value
    wilcoxon_p = compute_wilcoxon_p(values, alternative)
    
    # Format values
    median_formatted = format_percent_2sig(median_val)
    iqr_formatted = format_iqr_2sig_dash(q1, q3)
    p_formatted = format_pvalue_2sig(wilcoxon_p)
    
    # Condition display name (with Severe/Extreme for PASS B)
    cond_display = format_condition_display(condition, pass_name)
    
    # Print progress
    metric_display = METRIC_DISPLAY[metric]
    print(f"  {timescale} | {cond_display} | {metric_display}: n={n_sites}, median={median_formatted}%, IQR={iqr_formatted}, p={p_formatted}")
    
    # Store result
    results.append({
        'Timescale': timescale.replace('_', '-'),
        'Hydroclimatic condition': cond_display,
        'Metric': metric_display,
        'N sites': n_sites,
        'Median % change': median_formatted,
        'IQR (%)': iqr_formatted,
        'Wilcoxon p': p_formatted
    })


def print_table_s4(df):
    """Print table in readable format"""
    
    print("\n" + "=" * 120)
    print("TABLE S4: Site-level median WUE responses under hydroclimatic anomalies")
    print("=" * 120)
    print("\nValues are site-level medians of percent change relative to NN baseline conditions.")
    print("IQR = interquartile range (Q1-Q3).")
    print("Wilcoxon signed-rank test: H1: median < 0 for dry conditions, H1: median > 0 for wet conditions.")
    print("Bold p-values indicate statistical significance at α = 0.05.\n")
    
    # Print header
    header = f"{'Timescale':<14} {'Condition':<22} {'Metric':<10} {'N':<6} {'Median % change':<18} {'IQR (%)':<20} {'Wilcoxon p':<12}"
    print(header)
    print("-" * 110)
    
    for _, row in df.iterrows():
        timescale = row['Timescale']
        condition = row['Hydroclimatic condition']
        metric = row['Metric']
        n = row['N sites']
        median = row['Median % change']
        iqr = row['IQR (%)']
        p = row['Wilcoxon p']
        
        # Add asterisk for significant p-values
        p_display = p
        if p not in ['--', '<0.001']:
            try:
                p_val = float(p)
                if p_val < 0.05:
                    p_display = f"{p}*"
            except:
                pass
        elif p == '<0.001':
            p_display = '<0.001*'
        
        print(f"{timescale:<14} {condition:<22} {metric:<10} {n:<6} {median:<18} {iqr:<20} {p_display:<12}")
    
    print("-" * 110)
    print("\n* p < 0.05 (Wilcoxon signed-rank test)")
    print("  Dry conditions: testing decrease (median < 0)")
    print("  Wet conditions: testing increase (median > 0)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    
    print("\n" + "=" * 100)
    print("GENERATING TABLE S4")
    print("=" * 100)
    
    table = generate_table_s4()
    
    if table is not None:
        print(f"\n✅ Table S4 complete!")
        print(f"   Output file: {OUTPUT_FILE}")
        print(f"   Rows: {len(table)}")
    else:
        print("\n❌ Failed to generate table. Please check file path and columns.")