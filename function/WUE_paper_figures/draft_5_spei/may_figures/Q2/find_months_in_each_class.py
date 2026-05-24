# -*- coding: utf-8 -*-
"""
STANDALONE SCRIPT: Count observations (monthly records) per class
Run this AFTER the comprehensive SPEI workflow

This script reads:
1. The original monthly data file (has all monthly observations per site per SPEI class)
2. The site-level details CSV (to know which sites pass strict triple intersection)

Outputs: Number of monthly observations per condition class for each panel
"""

import pandas as pd
import numpy as np
import os

# ============================================================================
# FILE PATHS (same as your main workflow)
# ============================================================================
ROOT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
DATA_DIR = os.path.join(ROOT_DIR, "data_products", "results")
OUTPUT_DIR = os.path.join(ROOT_DIR, "data_products", "results", "SPEI_analysis_results")

# Input files
MONTHLY_DATA_FILE = os.path.join(DATA_DIR, "monthly_data_after_outlier_removal.csv")
SITE_DETAILS_FILE = os.path.join(OUTPUT_DIR, "SPEI_site_level_details_FINAL.csv")

# ============================================================================
# DEFINE CLASS GROUPINGS (matching your workflow)
# ============================================================================

# PASS B: Severity grouped
SEVERITY_MERGED = {
    'Mild Wet': ['MW', 'MoW'],
    'Severe Wet': ['SW', 'EW'],
    'Mild Dry': ['MD', 'MoD'],
    'Severe Dry': ['SD', 'ED']
}

# PASS C: Fully merged wet vs dry
FULLY_MERGED = {
    'Dry (all)': ['MD', 'MoD', 'SD', 'ED'],
    'Wet (all)': ['MW', 'MoW', 'SW', 'EW']
}

GROUPING_PASSES = {
    'PASS B': SEVERITY_MERGED,
    'PASS C': FULLY_MERGED
}

# SPEI timescales (from your monthly data - these are columns)
# Note: Your monthly data has columns like 'SPEI_1_Cat', 'SPEI_3_Cat', etc.
SPEI_TIMESCALES = [1, 3, 6, 12, 24, 36, 48]

# WUE metrics
WUE_METRICS = ['WUE', 'WUE_eva', 'WUE_tra']

# ============================================================================
# FUNCTION: Get strict triple intersection sites from site-level details
# ============================================================================

def get_triple_intersection_sites(site_details_df, spei_timescale, pass_name, condition, metric):
    """
    Get sites that passed strict triple intersection for a specific combination
    Uses the site-level details CSV which already has the strict intersection applied
    """
    filtered = site_details_df[
        (site_details_df['SPEI_Timescale'] == f"SPEI_{spei_timescale}") &
        (site_details_df['Pass'] == pass_name) &
        (site_details_df['Condition'] == condition) &
        (site_details_df['WUE_Metric'] == metric) &
        (site_details_df['Salinity'] == 'All')  # All sites combined
    ]
    
    if filtered.empty:
        return set()
    
    # Get sites that had valid data (not insufficient)
    valid_sites = filtered[~filtered['Insufficient_Data_Flag']]['Site'].unique()
    return set(valid_sites)

# ============================================================================
# FUNCTION: Count monthly observations per class
# ============================================================================

def count_observations_per_class(monthly_df, site_list, spei_timescale, class_list):
    """
    Count number of monthly observations for given sites and SPEI classes
    """
    # Get the appropriate SPEI category column
    spei_col = f'SPEI_{spei_timescale}_Cat'
    
    if spei_col not in monthly_df.columns:
        print(f"  Warning: {spei_col} not found in monthly data")
        return 0
    
    # Filter: sites in our list AND SPEI class in class_list
    filtered = monthly_df[
        (monthly_df['site_name'].isin(site_list)) &
        (monthly_df[spei_col].isin(class_list))
    ]
    
    return len(filtered)

# ============================================================================
# FUNCTION: Count observations for all conditions in a panel
# ============================================================================

def count_panel_observations(monthly_df, site_details_df, pass_name, groupings, spei_timescale):
    """
    Count monthly observations for all conditions in a panel
    """
    results = []
    
    for condition, class_list in groupings.items():
        for metric in WUE_METRICS:
            # Get sites that passed strict triple intersection for this combo
            sites = get_triple_intersection_sites(
                site_details_df, spei_timescale, pass_name, condition, metric
            )
            
            if len(sites) == 0:
                continue
            
            # Count observations
            n_obs = count_observations_per_class(
                monthly_df, sites, spei_timescale, class_list
            )
            
            # Also get site count from site_details
            site_count = len(sites)
            
            results.append({
                'SPEI_Timescale': f"SPEI_{spei_timescale}",
                'Pass': pass_name,
                'Condition': condition,
                'WUE_Metric': metric,
                'N_Sites': site_count,
                'N_Monthly_Observations': n_obs,
                'Mean_Obs_Per_Site': round(n_obs / site_count, 1) if site_count > 0 else 0
            })
    
    return results

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("=" * 80)
    print("STANDALONE SCRIPT: Monthly Observation Counts Per Class")
    print("=" * 80)
    
    # Check if files exist
    if not os.path.exists(MONTHLY_DATA_FILE):
        print(f"\n❌ ERROR: Monthly data file not found:")
        print(f"   {MONTHLY_DATA_FILE}")
        print(f"\n   Please run the main SPEI workflow first!")
        return
    
    if not os.path.exists(SITE_DETAILS_FILE):
        print(f"\n❌ ERROR: Site details file not found:")
        print(f"   {SITE_DETAILS_FILE}")
        print(f"\n   Please run the main SPEI workflow first!")
        return
    
    # Load data
    print("\n📂 Loading files...")
    monthly_df = pd.read_csv(MONTHLY_DATA_FILE)
    site_details_df = pd.read_csv(SITE_DETAILS_FILE)
    
    print(f"   Monthly data: {len(monthly_df):,} rows")
    print(f"   Site details: {len(site_details_df):,} rows")
    
    # Get unique sites in monthly data
    all_monthly_sites = set(monthly_df['site_name'].unique())
    print(f"   Unique sites in monthly data: {len(all_monthly_sites)}")
    
    # Collect results for all combinations
    all_results = []
    
    print("\n" + "=" * 80)
    print("COUNTING OBSERVATIONS PER CLASS")
    print("=" * 80)
    
    for spei_timescale in SPEI_TIMESCALES:
        print(f"\n{'─'*60}")
        print(f"SPEI-{spei_timescale}")
        print(f"{'─'*60}")
        
        for pass_name, groupings in GROUPING_PASSES.items():
            print(f"\n  {pass_name}:")
            
            results = count_panel_observations(
                monthly_df, site_details_df, pass_name, groupings, spei_timescale
            )
            
            if results:
                for r in results:
                    print(f"    {r['Condition']} | {r['WUE_Metric']}: "
                          f"{r['N_Sites']} sites, "
                          f"{r['N_Monthly_Observations']:,} obs "
                          f"({r['Mean_Obs_Per_Site']:.1f} avg/site)")
                    all_results.append(r)
            else:
                print(f"    No data found")
    
    # Create DataFrame and save
    if all_results:
        results_df = pd.DataFrame(all_results)
        
        # Save to CSV
        output_file = os.path.join(OUTPUT_DIR, "Monthly_Observations_Per_Class.csv")
        results_df.to_csv(output_file, index=False)
        print(f"\n{'='*80}")
        print(f"✅ Results saved to: {output_file}")
        
        # Print summary table
        print("\n" + "=" * 80)
        print("SUMMARY TABLE - Panel A (All conditions) vs Panel B (Severe/Extreme)")
        print("=" * 80)
        
        # Panel A: PASS C conditions (Dry all, Wet all)
        print("\n📊 PANEL A - All conditions (Dry/Wet):")
        panel_a = results_df[results_df['Pass'] == 'PASS C']
        for metric in WUE_METRICS:
            metric_data = panel_a[panel_a['WUE_Metric'] == metric]
            if not metric_data.empty:
                for _, row in metric_data.iterrows():
                    print(f"  {row['SPEI_Timescale']} | {row['Condition']} | {metric}: "
                          f"{row['N_Sites']} sites, {row['N_Monthly_Observations']:,} obs")
        
        # Panel B: PASS B severe conditions
        print("\n📊 PANEL B - Severe/Extreme conditions:")
        panel_b = results_df[
            (results_df['Pass'] == 'PASS B') & 
            (results_df['Condition'].isin(['Severe Dry', 'Severe Wet']))
        ]
        for metric in WUE_METRICS:
            metric_data = panel_b[panel_b['WUE_Metric'] == metric]
            if not metric_data.empty:
                for _, row in metric_data.iterrows():
                    print(f"  {row['SPEI_Timescale']} | {row['Condition']} | {metric}: "
                          f"{row['N_Sites']} sites, {row['N_Monthly_Observations']:,} obs")
        
        # Also show mild conditions for reference
        print("\n📊 (Reference) Mild conditions:")
        mild = results_df[
            (results_df['Pass'] == 'PASS B') & 
            (results_df['Condition'].isin(['Mild Dry', 'Mild Wet']))
        ]
        for metric in WUE_METRICS:
            metric_data = mild[mild['WUE_Metric'] == metric]
            if not metric_data.empty:
                for _, row in metric_data.iterrows():
                    print(f"  {row['SPEI_Timescale']} | {row['Condition']} | {metric}: "
                          f"{row['N_Sites']} sites, {row['N_Monthly_Observations']:,} obs")
    
    else:
        print("\n❌ No results generated")
    
    print("\n" + "=" * 80)
    print("SCRIPT COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
    
#########################################################################################################

# -*- coding: utf-8 -*-
"""
PRECISE OBSERVATION COUNTS FOR TWO-PANEL FIGURE
Matches EXACTLY the categories used in Combined_Figure_TwoPanel.png

Panel A (All conditions): SPEI-6/48 × Dry(all)/Wet(all)
Panel B (Severe/Extreme): SPEI-6/48 × Severe Dry/Severe Wet
Metrics: WUE_ET (WUE) and WUE_T (WUE_tra) only
Uses strict double intersection (sites with BOTH metrics)
"""

import pandas as pd
import numpy as np
import os

# ============================================================================
# FILE PATHS
# ============================================================================
ROOT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
DATA_DIR = os.path.join(ROOT_DIR, "data_products", "results")
OUTPUT_DIR = os.path.join(ROOT_DIR, "data_products", "results", "SPEI_analysis_results")

# Input files (from your workflow)
MONTHLY_DATA_FILE = os.path.join(DATA_DIR, "monthly_data_after_outlier_removal.csv")
SITE_DETAILS_FILE = os.path.join(OUTPUT_DIR, "SPEI_site_level_details_FINAL.csv")

# ============================================================================
# EXACT FIGURE PARAMETERS (matching your figure workflow)
# ============================================================================

# Panel A conditions (All conditions)
PANEL_A_CONDITIONS = {
    'Dry (all)': ['MD', 'MoD', 'SD', 'ED'],      # All dry classes
    'Wet (all)': ['MW', 'MoW', 'SW', 'EW']       # All wet classes
}

# Panel B conditions (Severe/Extreme)
PANEL_B_CONDITIONS = {
    'Severe Dry': ['SD', 'ED'],                   # Severe + Extreme dry
    'Severe Wet': ['SW', 'EW']                    # Severe + Extreme wet
}

# SPEI timescales used in figure
TIMESCALES = [6, 48]  # SPEI-6 and SPEI-48 only

# Metrics used in figure (double intersection)
FIGURE_METRICS = ['WUE', 'WUE_tra']  # WUE_ET and WUE_T

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_double_intersection_sites(site_details_df, spei_timescale, condition, pass_name):
    """
    Get sites that passed strict DOUBLE intersection (both WUE and WUE_tra)
    Matches exactly the logic in load_and_recompute_all_stats() from figure workflow
    """
    spei_col = f"SPEI_{spei_timescale}"
    
    # Filter to get sites for this specific combination
    filtered = site_details_df[
        (site_details_df['SPEI_Timescale'] == spei_col) &
        (site_details_df['Pass'] == pass_name) &
        (site_details_df['Condition'] == condition) &
        (site_details_df['Salinity'] == 'All') &
        (site_details_df['WUE_Metric'].isin(FIGURE_METRICS))
    ]
    
    if filtered.empty:
        return set()
    
    # Get sites for each metric
    wue_sites = set(filtered[filtered['WUE_Metric'] == 'WUE']['Site'].unique())
    tra_sites = set(filtered[filtered['WUE_Metric'] == 'WUE_tra']['Site'].unique())
    
    # STRICT DOUBLE INTERSECTION (both metrics must exist)
    shared_sites = wue_sites.intersection(tra_sites)
    
    return shared_sites


def count_observations_for_condition(monthly_df, sites, spei_timescale, class_list):
    """
    Count monthly observations for given sites and SPEI classes
    Each row in monthly data = one monthly observation
    """
    spei_col = f'SPEI_{spei_timescale}_Cat'
    
    if spei_col not in monthly_df.columns:
        print(f"  WARNING: {spei_col} not found in monthly data")
        return 0
    
    # Filter: sites in our list AND SPEI class in class_list
    filtered = monthly_df[
        (monthly_df['site_name'].isin(sites)) &
        (monthly_df[spei_col].isin(class_list))
    ]
    
    return len(filtered)


def analyze_panel(monthly_df, site_details_df, conditions_dict, pass_name, panel_name):
    """
    Analyze a panel and return formatted results
    """
    print(f"\n{'='*70}")
    print(f"{panel_name}")
    print(f"Pass: {pass_name}")
    print(f"{'='*70}")
    
    results = []
    
    for spei in TIMESCALES:
        print(f"\n  SPEI-{spei}:")
        
        for condition, class_list in conditions_dict.items():
            # Get sites with double intersection
            sites = get_double_intersection_sites(
                site_details_df, spei, condition, pass_name
            )
            
            if len(sites) == 0:
                print(f"    {condition}: No sites found")
                continue
            
            # Count observations
            n_obs = count_observations_for_condition(
                monthly_df, sites, spei, class_list
            )
            
            n_sites = len(sites)
            obs_per_site = n_obs / n_sites if n_sites > 0 else 0
            
            # Store result
            results.append({
                'Panel': panel_name,
                'SPEI': f'SPEI-{spei}',
                'Condition': condition,
                'N_Sites': n_sites,
                'N_Observations': n_obs,
                'Mean_Obs_Per_Site': round(obs_per_site, 1)
            })
            
            # Print formatted output
            print(f"    {condition}: N_sites={n_sites}, n_obs={n_obs:,} ({obs_per_site:.1f} obs/site)")
    
    return results


def print_compact_summary(all_results):
    """
    Print a compact, publication-ready summary table
    """
    print("\n" + "="*90)
    print("PUBLICATION-READY SUMMARY: Observation Counts per Figure Panel")
    print("="*90)
    
    # Panel A summary
    print("\n📊 PANEL A: All conditions (Dry/Wet)")
    print("-" * 60)
    panel_a = [r for r in all_results if r['Panel'] == 'PANEL A']
    
    print(f"{'SPEI':<8} {'Condition':<15} {'Metric':<12} {'Sites (N)':<12} {'Observations (n)':<20} {'Mean/site':<12}")
    print("-" * 80)
    
    for r in panel_a:
        # For each condition, same sites for both metrics (double intersection)
        print(f"{r['SPEI']:<8} {r['Condition']:<15} {'WUE_ET & WUE_T':<12} {r['N_Sites']:<12} {r['N_Observations']:<20,} {r['Mean_Obs_Per_Site']:<12}")
    
    # Panel B summary
    print("\n📊 PANEL B: Severe/Extreme conditions")
    print("-" * 60)
    panel_b = [r for r in all_results if r['Panel'] == 'PANEL B']
    
    print(f"{'SPEI':<8} {'Condition':<15} {'Metric':<12} {'Sites (N)':<12} {'Observations (n)':<20} {'Mean/site':<12}")
    print("-" * 80)
    
    for r in panel_b:
        print(f"{r['SPEI']:<8} {r['Condition']:<15} {'WUE_ET & WUE_T':<12} {r['N_Sites']:<12} {r['N_Observations']:<20,} {r['Mean_Obs_Per_Site']:<12}")


def print_vertical_table(all_results):
    """
    Print vertical format matching your figure's N values
    """
    print("\n" + "="*70)
    print("VERTICAL FORMAT (matches figure y-axis labels)")
    print("="*70)
    
    print("\nPANEL A (All conditions):")
    print("─" * 50)
    for r in [x for x in all_results if x['Panel'] == 'PANEL A']:
        cond_short = r['Condition'].replace(' (all)', '')
        print(f"  {r['SPEI']} | {cond_short:<10} | N_sites={r['N_Sites']}, n_obs={r['N_Observations']:,}")
    
    print("\nPANEL B (Severe/Extreme):")
    print("─" * 50)
    for r in [x for x in all_results if x['Panel'] == 'PANEL B']:
        cond_short = r['Condition'].replace('Severe ', 'Ext. ')
        print(f"  {r['SPEI']} | {cond_short:<10} | N_sites={r['N_Sites']}, n_obs={r['N_Observations']:,}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("="*80)
    print("PRECISE OBSERVATION COUNTS FOR TWO-PANEL FIGURE")
    print("Matching exactly: Combined_Figure_TwoPanel.png")
    print("="*80)
    
    # Check files
    if not os.path.exists(MONTHLY_DATA_FILE):
        print(f"\n❌ ERROR: Monthly data not found: {MONTHLY_DATA_FILE}")
        return
    
    if not os.path.exists(SITE_DETAILS_FILE):
        print(f"\n❌ ERROR: Site details not found: {SITE_DETAILS_FILE}")
        print("   Please run the main SPEI workflow first!")
        return
    
    # Load data
    print("\n📂 Loading data...")
    monthly_df = pd.read_csv(MONTHLY_DATA_FILE)
    site_details_df = pd.read_csv(SITE_DETAILS_FILE)
    
    print(f"   Monthly data: {len(monthly_df):,} rows")
    print(f"   Site details: {len(site_details_df):,} rows")
    print(f"   SPEI timescales in monthly data: {[c for c in monthly_df.columns if 'SPEI_' in c and '_Cat' in c]}")
    
    # Verify monthly data has required SPEI columns
    for spei in TIMESCALES:
        col = f'SPEI_{spei}_Cat'
        if col not in monthly_df.columns:
            print(f"\n⚠️ WARNING: {col} not found in monthly data!")
            print(f"   Available SPEI columns: {[c for c in monthly_df.columns if 'SPEI_' in c]}")
    
    # Analyze Panel A (PASS C - fully merged)
    print("\n" + "="*70)
    print("ANALYZING PANEL A: All conditions (Dry/Wet)")
    print("Using PASS C (fully merged wet vs dry)")
    print(f"Strict double intersection: sites with BOTH {FIGURE_METRICS}")
    print("="*70)
    
    results_a = analyze_panel(
        monthly_df, site_details_df, 
        PANEL_A_CONDITIONS, 'PASS C', 'PANEL A'
    )
    
    # Analyze Panel B (PASS B - severe/extreme only)
    print("\n" + "="*70)
    print("ANALYZING PANEL B: Severe/Extreme conditions")
    print("Using PASS B (severity grouped)")
    print(f"Strict double intersection: sites with BOTH {FIGURE_METRICS}")
    print("="*70)
    
    results_b = analyze_panel(
        monthly_df, site_details_df,
        PANEL_B_CONDITIONS, 'PASS B', 'PANEL B'
    )
    
    # Combine results
    all_results = results_a + results_b
    
    # Print formatted outputs
    print_compact_summary(all_results)
    print_vertical_table(all_results)
    
    # Save to CSV
    if all_results:
        results_df = pd.DataFrame(all_results)
        output_file = os.path.join(OUTPUT_DIR, "Figure_Observation_Counts.csv")
        results_df.to_csv(output_file, index=False)
        print(f"\n✅ Saved to: {output_file}")
    
    # Also print a quick reference table matching figure rows exactly
    print("\n" + "="*70)
    print("FIGURE ROW REFERENCE (exactly as in Combined_Figure_TwoPanel.png)")
    print("="*70)
    
    print("\nPANEL A (Left panel) - All conditions:")
    print("┌─────────────┬──────────────────┬────────────┬──────────────┐")
    print("│ Row         │ Condition        │ N_sites    │ n_obs        │")
    print("├─────────────┼──────────────────┼────────────┼──────────────┤")
    
    row_order_a = [
        ('SPEI-6', 'Dry (all)'),
        ('SPEI-48', 'Dry (all)'),
        ('SPEI-6', 'Wet (all)'),
        ('SPEI-48', 'Wet (all)')
    ]
    
    for spei, cond in row_order_a:
        match = next((r for r in results_a if r['SPEI'] == spei and r['Condition'] == cond), None)
        if match:
            print(f"│ {spei:<10} │ {cond:<16} │ {match['N_Sites']:<10} │ {match['N_Observations']:<12,} │")
    
    print("└─────────────┴──────────────────┴────────────┴──────────────┘")
    
    print("\nPANEL B (Right panel) - Severe/Extreme:")
    print("┌─────────────┬──────────────────┬────────────┬──────────────┐")
    print("│ Row         │ Condition        │ N_sites    │ n_obs        │")
    print("├─────────────┼──────────────────┼────────────┼──────────────┤")
    
    row_order_b = [
        ('SPEI-6', 'Severe Dry'),
        ('SPEI-48', 'Severe Dry'),
        ('SPEI-6', 'Severe Wet'),
        ('SPEI-48', 'Severe Wet')
    ]
    
    for spei, cond in row_order_b:
        match = next((r for r in results_b if r['SPEI'] == spei and r['Condition'] == cond), None)
        if match:
            cond_short = cond.replace('Severe ', 'Ext. ')
            print(f"│ {spei:<10} │ {cond_short:<16} │ {match['N_Sites']:<10} │ {match['N_Observations']:<12,} │")
    
    print("└─────────────┴──────────────────┴────────────┴──────────────┘")
    
    print("\n" + "="*80)
    print("✅ SCRIPT COMPLETE")
    print("="*80)
    print("\nNOTES:")
    print(f"  • N_sites = Number of sites with BOTH WUE_ET and WUE_T (strict double intersection)")
    print(f"  • n_obs = Total monthly observations across all sites for the specified SPEI classes")
    print(f"  • Each observation = one month of data at a site under the given SPEI condition")
    print("="*80)


if __name__ == "__main__":
    main()    
    
    
    
    
