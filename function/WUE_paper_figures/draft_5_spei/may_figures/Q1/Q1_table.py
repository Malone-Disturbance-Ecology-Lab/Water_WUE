# -*- coding: utf-8 -*-
"""
COMPLETE TABLE GENERATOR - CORRECTED VERSION
Uses TWO separate data sources matching Figure 1 and Figure 2 exactly

OVERALL section (57 sites): wue_site_level_NN_medians_SPEI_1.csv (already NN, no salinity filter)
Ecosystem section (52 sites): wue_site_level_summary_SPEI_1.csv (NN + Freshwater/Saline/Upland only)

@author: WUE Analysis Pipeline
Date: 2026-05-18
"""

import pandas as pd
import numpy as np
import os
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
OUTPUT_DIR = BASE_DIR

# Data sources (CORRECTED)
FIGURE1_FILE = os.path.join(BASE_DIR, "wue_site_level_NN_medians_SPEI_1.csv")
FIGURE2_FILE = os.path.join(BASE_DIR, "wue_site_level_summary_SPEI_1.csv")

ECOSYSTEM_TYPES = ['Freshwater', 'Saline', 'Upland']
METRICS = ['WUE', 'WUE_eva', 'WUE_tra']
METRIC_DISPLAY = {
    'WUE': 'WUE_ET',
    'WUE_eva': 'WUE_E',
    'WUE_tra': 'WUE_T'
}

# =============================================================================
# DATA LOADING FUNCTIONS (CORRECTED)
# =============================================================================

def load_overall_data():
    """
    Load Figure 1 dataset (57 sites, includes Brackish, already NN-filtered)
    File: wue_site_level_NN_medians_SPEI_1.csv
    - No SPEI_Class filter (already NN)
    - No salinity filter (keep all sites including Brackish)
    - Apply strict triple intersection
    """
    print("=" * 80)
    print("LOADING OVERALL DATA (57 sites - includes Brackish)")
    print("=" * 80)
    
    df = pd.read_csv(FIGURE1_FILE)
    print(f"Loaded {len(df)} rows from {FIGURE1_FILE}")
    print(f"Columns: {df.columns.tolist()}")
    
    # NO SPEI_Class filter - file is already NN-only
    # NO salinity filter - keep all sites
    
    # Remove duplicates if any
    df_dedup = df.drop_duplicates(
        subset=["site_name", "WUE_Metric"]
    ).copy()
    print(f"After deduplication: {len(df_dedup)} rows")
    
    # Strict triple intersection (sites with all three metrics)
    sites_with_wue = set(df_dedup[df_dedup["WUE_Metric"] == "WUE"]["site_name"].unique())
    sites_with_eva = set(df_dedup[df_dedup["WUE_Metric"] == "WUE_eva"]["site_name"].unique())
    sites_with_tra = set(df_dedup[df_dedup["WUE_Metric"] == "WUE_tra"]["site_name"].unique())
    
    shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)
    print(f"\nStrict triple intersection: {len(shared_sites)} sites")
    
    final_data = df_dedup[df_dedup["site_name"].isin(shared_sites)].copy()
    
    for metric in METRICS:
        n_sites = final_data[final_data["WUE_Metric"] == metric]["site_name"].nunique()
        print(f"  {METRIC_DISPLAY[metric]}: {n_sites} sites")
    
    return final_data, shared_sites


def load_ecosystem_data():
    """
    Load Figure 2 dataset (52 sites, Freshwater/Saline/Upland only)
    File: wue_site_level_summary_SPEI_1.csv
    - Filter SPEI_Class == 'NN'
    - Filter Salinity_Category in ['Freshwater', 'Saline', 'Upland']
    - Apply strict triple intersection
    """
    print("\n" + "=" * 80)
    print("LOADING ECOSYSTEM DATA (52 sites - Freshwater/Saline/Upland only)")
    print("=" * 80)
    
    df = pd.read_csv(FIGURE2_FILE)
    print(f"Loaded {len(df)} rows from {FIGURE2_FILE}")
    
    # Filter to near-normal conditions
    df_nn = df[df["SPEI_Class"] == "NN"].copy()
    print(f"After SPEI_Class == 'NN': {len(df_nn)} rows")
    
    # Filter to allowed ecosystem types
    df_filtered = df_nn[df_nn["Salinity_Category"].isin(ECOSYSTEM_TYPES)].copy()
    print(f"After filtering to {ECOSYSTEM_TYPES}: {len(df_filtered)} rows")
    
    # Remove duplicates
    df_dedup = df_filtered.drop_duplicates(
        subset=["site_name", "Salinity_Category", "WUE_Metric"]
    ).copy()
    print(f"After deduplication: {len(df_dedup)} rows")
    
    # Strict triple intersection
    sites_with_wue = set(df_dedup[df_dedup["WUE_Metric"] == "WUE"]["site_name"].unique())
    sites_with_eva = set(df_dedup[df_dedup["WUE_Metric"] == "WUE_eva"]["site_name"].unique())
    sites_with_tra = set(df_dedup[df_dedup["WUE_Metric"] == "WUE_tra"]["site_name"].unique())
    
    shared_sites = sites_with_wue.intersection(sites_with_eva).intersection(sites_with_tra)
    print(f"\nStrict triple intersection (Freshwater/Saline/Upland only): {len(shared_sites)} sites")
    
    final_data = df_dedup[df_dedup["site_name"].isin(shared_sites)].copy()
    
    for metric in METRICS:
        n_sites = final_data[final_data["WUE_Metric"] == metric]["site_name"].nunique()
        print(f"  {METRIC_DISPLAY[metric]}: {n_sites} sites")
    
    return final_data, shared_sites


# =============================================================================
# STATISTICAL FUNCTIONS
# =============================================================================

def compute_median_iqr(values):
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    median = np.median(values)
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    return median, q1, q3


def compute_range_percentile(values):
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    val_min = np.min(values)
    val_max = np.max(values)
    p90 = np.percentile(values, 90)
    return val_min, val_max, p90


def compare_wue_metrics_globally(data):
    """
    Compare WUE_ET, WUE_E, WUE_T across all sites (57 sites)
    Using Kruskal-Wallis test
    """
    values = []
    for metric in METRICS:
        metric_data = data[data["WUE_Metric"] == metric]
        vals = metric_data["WUE_median"].dropna().values  # Note: column name is WUE_median in Figure 1 file
        values.append(vals)
    
    if len(values) >= 2 and all(len(v) > 0 for v in values):
        h_stat, p_value = kruskal(*values)
        return f"WUE metrics differ (p = {p_value:.3f})" if p_value < 0.05 else f"WUE metrics do not differ (p = {p_value:.3f})"
    return "Insufficient data"


def ecosystem_comparison_for_metric(data, metric):
    """
    Compare Freshwater, Saline, Upland for a specific metric (52 sites)
    """
    metric_data = data[data["WUE_Metric"] == metric]
    
    groups = {}
    for eco in ECOSYSTEM_TYPES:
        vals = metric_data[metric_data["Salinity_Category"] == eco]["Median"].dropna().values
        if len(vals) > 0:
            groups[eco] = vals
    
    if len(groups) < 2:
        return "Insufficient data"
    
    # Kruskal-Wallis
    h_stat, p_kw = kruskal(*groups.values())
    
    # Pairwise Mann-Whitney with Bonferroni
    eco_names = list(groups.keys())
    pairs = []
    p_vals = []
    
    for i in range(len(eco_names)):
        for j in range(i+1, len(eco_names)):
            eco1, eco2 = eco_names[i], eco_names[j]
            u_stat, p_pair = mannwhitneyu(groups[eco1], groups[eco2], alternative='two-sided')
            pairs.append(f"{eco1} vs {eco2}")
            p_vals.append(p_pair)
    
    if p_vals:
        reject, p_corrected, _, _ = multipletests(p_vals, method='bonferroni')
        sig_pairs = [pairs[idx] for idx, p_val in enumerate(p_corrected) if p_val < 0.05]
        
        if sig_pairs:
            pairs_text = "; ".join([f"{pair} (p={p_corrected[idx]:.3f})" 
                                    for idx, pair in enumerate(pairs) if p_corrected[idx] < 0.05])
            return f"Kruskal-Wallis p = {p_kw:.3f}\nSignificant: {pairs_text}"
        else:
            return f"Kruskal-Wallis p = {p_kw:.3f}\nNo significant differences"
    else:
        return f"Kruskal-Wallis p = {p_kw:.3f}"


# =============================================================================
# TABLE GENERATION
# =============================================================================

def generate_table(overall_data, ecosystem_data):
    """
    Generate final table using TWO separate data sources:
    - overall_data: 57 sites for OVERALL section (uses 'WUE_median' column)
    - ecosystem_data: 52 sites for Ecosystem section (uses 'Median' column)
    """
    
    rows = []
    
    # ========================================================================
    # SECTION 1: OVERALL (All 57 sites from Figure 1)
    # ========================================================================
    overall_comparison = compare_wue_metrics_globally(overall_data)
    total_sites_57 = overall_data["site_name"].nunique()
    
    for idx, metric in enumerate(METRICS):
        metric_display = METRIC_DISPLAY[metric]
        # Use 'WUE_median' column name (Figure 1 file)
        metric_data = overall_data[overall_data["WUE_Metric"] == metric]
        values = metric_data["WUE_median"].dropna().values
        
        median, q1, q3 = compute_median_iqr(values)
        val_min, val_max, p90 = compute_range_percentile(values)
        
        median_iqr = f"{median:.2f}\n{q1:.2f}-{q3:.2f}" if not np.isnan(median) else "--"
        range_text = f"{val_min:.2f}-{val_max:.2f}" if not np.isnan(val_min) else "--"
        p90_text = f"{p90:.2f}" if not np.isnan(p90) else "--"
        
        # Only show comparison text in first metric row
        group_comp = overall_comparison if idx == 0 else ""
        
        rows.append({
            "Section": f"OVERALL (All {total_sites_57} sites)",
            "Metric": metric_display,
            "Median / IQR": median_iqr,
            "Range (min-max)": range_text,
            "90th percentile": p90_text,
            "Group Comparison": group_comp
        })
    
    # Add spacer row
    rows.append({
        "Section": "", "Metric": "", "Median / IQR": "", 
        "Range (min-max)": "", "90th percentile": "", "Group Comparison": ""
    })
    
    # ========================================================================
    # SECTION 2: ECOSYSTEM TYPES (52 sites from Figure 2)
    # ========================================================================
    # Pre-compute ecosystem comparisons for each metric
    ecosystem_comparisons = {}
    for metric in METRICS:
        ecosystem_comparisons[metric] = ecosystem_comparison_for_metric(ecosystem_data, metric)
    
    for ecosystem in ECOSYSTEM_TYPES:
        # Get n (number of sites) for this ecosystem
        eco_sites = ecosystem_data[ecosystem_data["Salinity_Category"] == ecosystem]["site_name"].unique()
        n_sites = len(eco_sites)
        
        for idx, metric in enumerate(METRICS):
            metric_display = METRIC_DISPLAY[metric]
            metric_data = ecosystem_data[(ecosystem_data["WUE_Metric"] == metric) & 
                                          (ecosystem_data["Salinity_Category"] == ecosystem)]
            values = metric_data["Median"].dropna().values  # Note: 'Median' column in Figure 2 file
            
            median, q1, q3 = compute_median_iqr(values)
            val_min, val_max, p90 = compute_range_percentile(values)
            
            median_iqr = f"{median:.2f}\n{q1:.2f}-{q3:.2f}" if not np.isnan(median) else "--"
            range_text = f"{val_min:.2f}-{val_max:.2f}" if not np.isnan(val_min) else "--"
            p90_text = f"{p90:.2f}" if not np.isnan(p90) else "--"
            
            # Show ecosystem name only in first row of each ecosystem
            section_label = f"{ecosystem} (n={n_sites})" if idx == 0 else ""
            
            # Show comparison only in first ecosystem row (Freshwater)
            if ecosystem == "Freshwater":
                group_comp = ecosystem_comparisons[metric]
            else:
                group_comp = ""
            
            rows.append({
                "Section": section_label,
                "Metric": metric_display,
                "Median / IQR": median_iqr,
                "Range (min-max)": range_text,
                "90th percentile": p90_text,
                "Group Comparison": group_comp
            })
    
    return pd.DataFrame(rows)


# =============================================================================
# OUTPUT FUNCTIONS
# =============================================================================

def print_table(df):
    """Print table in readable format"""
    print("\n" + "=" * 140)
    print("FINAL TABLE: Overall WUE (57 sites) + Ecosystem Comparison (52 sites)")
    print("=" * 140)
    print("\nNOTE: Overall statistics from Figure 1 dataset (includes Brackish, 57 sites)")
    print("      Ecosystem statistics from Figure 2 dataset (Freshwater/Saline/Upland only, 52 sites)")
    print("\n" + "-" * 140)
    
    # Print header
    header = f"{'Section':<35} {'Metric':<12} {'Median / IQR':<20} {'Range (min-max)':<18} {'90th percentile':<16} {'Group Comparison':<45}"
    print(header)
    print("-" * 145)
    
    for _, row in df.iterrows():
        section = str(row["Section"]) if pd.notna(row["Section"]) else ""
        metric = str(row["Metric"]) if pd.notna(row["Metric"]) else ""
        median_iqr = str(row["Median / IQR"]) if pd.notna(row["Median / IQR"]) else ""
        range_text = str(row["Range (min-max)"]) if pd.notna(row["Range (min-max)"]) else ""
        p90 = str(row["90th percentile"]) if pd.notna(row["90th percentile"]) else ""
        group_comp = str(row["Group Comparison"]) if pd.notna(row["Group Comparison"]) else ""
        
        # Replace newlines for cleaner printing
        median_iqr_display = median_iqr.replace('\n', ' / ')
        group_comp_display = group_comp.replace('\n', ' | ')
        
        # Truncate long text
        if len(group_comp_display) > 45:
            group_comp_display = group_comp_display[:42] + "..."
        
        print(f"{section:<35} {metric:<12} {median_iqr_display:<20} {range_text:<18} {p90:<16} {group_comp_display:<45}")
    
    print("-" * 145)
    
    # Add explanatory note
    print("\n" + "=" * 140)
    print("TABLE NOTES:")
    print("  - Overall statistics: 57-site strict triple-intersection dataset (includes Brackish and other coastal systems)")
    print("  - Ecosystem statistics: 52-site subset containing only Freshwater, Saline, and Upland ecosystems")
    print("  - Kruskal-Wallis and pairwise comparisons use the 52-site ecosystem dataset")
    print("  - Range = min-max of site-level medians")
    print("  - 90th percentile = 90th percentile of site-level medians")
    print("=" * 140)


def save_table(df, output_path):
    """Save table to CSV"""
    df.to_csv(output_path, index=False, encoding='utf-8')
    print(f"\n✅ Table saved to: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    print("\n" + "=" * 80)
    print("GENERATING TABLE: Overall (57 sites) + Ecosystem Comparison (52 sites)")
    print("=" * 80)
    
    # Load both data sources
    overall_data, shared_sites_57 = load_overall_data()
    ecosystem_data, shared_sites_52 = load_ecosystem_data()
    
    # Verify site counts
    print(f"\n✓ Figure 1 (OVERALL): {len(shared_sites_57)} sites")
    print(f"✓ Figure 2 (Ecosystem): {len(shared_sites_52)} sites")
    
    # Generate table
    table_df = generate_table(overall_data, ecosystem_data)
    
    # Print and save
    print_table(table_df)
    
    output_path = os.path.join(OUTPUT_DIR, "Table_Overall_57sites_Ecosystem_52sites_CORRECTED.csv")
    save_table(table_df, output_path)
    
    print("\n" + "=" * 80)
    print("TABLE GENERATION COMPLETE")
    print("=" * 80)