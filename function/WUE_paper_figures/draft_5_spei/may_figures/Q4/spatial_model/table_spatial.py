# -*- coding: utf-8 -*-
"""
COMPLETE TABLE GENERATION FOR Q4: Spatial Frequency, PRE/POST Breakpoint, and SPEI Shifts

Generates three tables dynamically from existing CSV outputs:
1. Table X: Regional spatial WUE_T decline frequency (supports Figure 11)
2. Table Y: Regional PRE vs POST breakpoint changes (supports Figure 13)
3. Table Sx: Breakpoint timing and SPEI shifts (supplementary)

FORMATTING RULES:
- Mean frequency, Median, IQR: 2 significant digits
- PRE mean, POST mean, Δ (POST-PRE): 2 significant digits
- Relative change (%): 1 decimal place
- Pixels >=0.50 (%): 1 decimal place (using >=, not unicode)
- Pixels with Δ <= -0.15 (%): 1 decimal place
- PRE SPEI, POST SPEI, ΔSPEI: 2 significant digits
- Pixels (n): whole number with commas

@author: WUE Analysis Pipeline
Date: 2026-05-19
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input file paths
SPATIAL_DIAGNOSTIC_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_time_aggregated\agg_FULLPERIOD_WUE_tra_COASTAL_REGIONS_DIAGNOSTICS.csv"

PREPOST_DIAGNOSTIC_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_pre_post_breakpoint\FIGURE3_PRE_POST_DELTA_DIAGNOSTICS_WUE_T.csv"

BREAKPOINT_SPEI_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\regional_SPEI48_pre_post_summary.csv"

# Output directory
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\tables_Q4"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output file paths
TABLE_X_OUTPUT = os.path.join(OUTPUT_DIR, "Table_X_Regional_Spatial_WUET_Decline_Frequency.csv")
TABLE_Y_OUTPUT = os.path.join(OUTPUT_DIR, "Table_Y_Regional_PRE_POST_Breakpoint_Changes.csv")
TABLE_SX_OUTPUT = os.path.join(OUTPUT_DIR, "Table_Sx_Breakpoint_Timing_SPEI_Shifts.csv")

# Region order
REGION_ORDER = ["Pacific", "Gulf", "Southeast Atlantic", "Atlantic North", "Alaska"]

# =============================================================================
# FORMATTING FUNCTIONS
# =============================================================================

def format_2sig(value):
    """Format with 2 significant digits."""
    if pd.isna(value):
        return '--'
    if value == 0:
        return '0.0'
    
    # For very small values (p-values), keep 3 decimals
    if 0 < abs(value) < 0.01:
        if abs(value) < 0.001:
            return "<0.001"
        else:
            return f"{value:.3f}"
    
    from math import floor, log10
    try:
        magnitude = floor(log10(abs(value)))
        if magnitude >= 0:
            if value >= 100:
                return f"{int(round(value))}"
            elif value >= 10:
                return f"{round(value, 1):.1f}"
            else:
                return f"{round(value, 2):.2f}"
        else:
            rounded = round(value, -magnitude + 1)
            formatted = f"{rounded:.{max(2, -magnitude + 2)}f}"
            if '.' in formatted:
                formatted = formatted.rstrip('0').rstrip('.')
            return formatted
    except:
        return f"{value:.2f}"


def format_percent_1dec(value):
    """Format percentage with 1 decimal place."""
    if pd.isna(value):
        return '--'
    return f"{value:.1f}"


def format_whole_number(value):
    """Format whole number with commas."""
    if pd.isna(value):
        return '--'
    return f"{int(value):,}"


def format_delta_with_sign(value):
    """Format delta with + sign (2 significant digits)."""
    if pd.isna(value):
        return '--'
    if value >= 0:
        return f"+{format_2sig(value)}"
    else:
        return f"{format_2sig(value)}"


def format_p_value(value):
    """Format p-value."""
    if pd.isna(value):
        return '--'
    if value < 0.001:
        return "<0.001"
    elif value < 0.01:
        return f"{value:.3f}"
    else:
        return f"{value:.2f}"


# =============================================================================
# TABLE X: Spatial Frequency (from FULLPERIOD diagnostic)
# =============================================================================

def generate_table_x():
    """Generate Table X from spatial diagnostic CSV."""
    print("\n" + "=" * 80)
    print("TABLE X: Regional Spatial WUE_T Decline Frequency")
    print("=" * 80)
    
    if not os.path.exists(SPATIAL_DIAGNOSTIC_FILE):
        print(f"❌ File not found: {SPATIAL_DIAGNOSTIC_FILE}")
        return None
    
    df = pd.read_csv(SPATIAL_DIAGNOSTIC_FILE)
    print(f"✅ Loaded from: {SPATIAL_DIAGNOSTIC_FILE}")
    
    results = []
    
    # Regional rows
    for region in REGION_ORDER:
        row = df[df['Region'] == region]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        
        results.append({
            'Region': region,
            'Pixels (n)': format_whole_number(row['n_pixels']),
            'Mean frequency': format_2sig(row['mean_frac_less']),
            'Median': format_2sig(row['median_frac_less']),
            'IQR': format_2sig(row['iqr_frac_less']),
            'Pixels >=0.50 (%)': format_percent_1dec(row['frac_less_ge_0.50'])
        })
    
    # Overall row: "All Valid Pixels"
    all_row = df[df['Region'] == 'All Valid Pixels']
    if len(all_row) > 0:
        all_data = all_row.iloc[0]
        results.append({
            'Region': 'All valid coastal pixels',
            'Pixels (n)': format_whole_number(all_data['n_pixels']),
            'Mean frequency': format_2sig(all_data['mean_frac_less']),
            'Median': format_2sig(all_data['median_frac_less']),
            'IQR': format_2sig(all_data['iqr_frac_less']),
            'Pixels >=0.50 (%)': format_percent_1dec(all_data['frac_less_ge_0.50'])
        })
    
    df_table = pd.DataFrame(results)
    df_table.to_csv(TABLE_X_OUTPUT, index=False, encoding='utf-8-sig')
    
    # Print to console
    print("\n" + "-" * 110)
    print("TABLE X OUTPUT:")
    print("-" * 110)
    print(df_table.to_string(index=False))
    print("-" * 110)
    print(f"✅ Saved: {TABLE_X_OUTPUT}")
    
    return df_table


# =============================================================================
# TABLE Y: PRE vs POST Changes (from PRE/POST diagnostic)
# =============================================================================

def generate_table_y():
    """Generate Table Y from PRE/POST diagnostic CSV."""
    print("\n" + "=" * 80)
    print("TABLE Y: Regional PRE vs POST Breakpoint Changes")
    print("=" * 80)
    
    if not os.path.exists(PREPOST_DIAGNOSTIC_FILE):
        print(f"❌ File not found: {PREPOST_DIAGNOSTIC_FILE}")
        return None
    
    df = pd.read_csv(PREPOST_DIAGNOSTIC_FILE)
    print(f"✅ Loaded from: {PREPOST_DIAGNOSTIC_FILE}")
    
    results = []
    
    # Regional rows
    for region in REGION_ORDER:
        row = df[df['Region'] == region]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        
        results.append({
            'Region': region,
            'PRE mean': format_2sig(row['PRE_mean']),
            'POST mean': format_2sig(row['POST_mean']),
            'Δ (POST-PRE)': format_delta_with_sign(row['Δ_mean']),
            'Relative change (%)': f"{row['pct_change']:+.1f}" if not pd.isna(row['pct_change']) else '--',
            'Pixels with Δ <= -0.15 (%)': format_percent_1dec(row['pct_large_neg']) if not pd.isna(row['pct_large_neg']) else '--'
        })
    
    # Overall row: "All Valid Pixels"
    all_row = df[df['Region'] == 'All Valid Pixels']
    if len(all_row) > 0:
        all_data = all_row.iloc[0]
        results.append({
            'Region': 'All valid coastal pixels',
            'PRE mean': format_2sig(all_data['PRE_mean']),
            'POST mean': format_2sig(all_data['POST_mean']),
            'Δ (POST-PRE)': format_delta_with_sign(all_data['Δ_mean']),
            'Relative change (%)': f"{all_data['pct_change']:+.1f}" if not pd.isna(all_data['pct_change']) else '--',
            'Pixels with Δ <= -0.15 (%)': format_percent_1dec(all_data['pct_large_neg']) if not pd.isna(all_data['pct_large_neg']) else '--'
        })
    
    df_table = pd.DataFrame(results)
    df_table.to_csv(TABLE_Y_OUTPUT, index=False, encoding='utf-8-sig')
    
    # Print to console
    print("\n" + "-" * 110)
    print("TABLE Y OUTPUT:")
    print("-" * 110)
    print(df_table.to_string(index=False))
    print("-" * 110)
    print(f"✅ Saved: {TABLE_Y_OUTPUT}")
    
    return df_table


# =============================================================================
# TABLE Sx: Breakpoint Timing and SPEI Shifts
# =============================================================================

def generate_table_sx():
    """Generate Table Sx from breakpoint SPEI summary CSV."""
    print("\n" + "=" * 80)
    print("TABLE Sx: Breakpoint Timing and SPEI Shifts")
    print("=" * 80)
    
    if not os.path.exists(BREAKPOINT_SPEI_FILE):
        print(f"❌ File not found: {BREAKPOINT_SPEI_FILE}")
        print("\nPlease add CSV saving to your merged breakpoint workflow.")
        return None
    
    df = pd.read_csv(BREAKPOINT_SPEI_FILE)
    print(f"✅ Loaded from: {BREAKPOINT_SPEI_FILE}")
    
    # Warn if Pacific breakpoint is still old
    pacific_row = df[df['Region'] == 'Pacific']
    if len(pacific_row) > 0:
        bp = str(pacific_row.iloc[0]['Breakpoint'])
        if '2021' in bp:
            print(f"   ⚠️ WARNING: Pacific breakpoint is '{bp}' (should be Apr 2014)")
            print("   Please update your merged breakpoint workflow to save correct CSV.")
    
    results = []
    
    for region in REGION_ORDER:
        row = df[df['Region'] == region]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        
        # Format breakpoint date nicely
        breakpoint = row['Breakpoint']
        if isinstance(breakpoint, str):
            parts = breakpoint.split()
            if len(parts) >= 2:
                if '-' in parts[0]:
                    from datetime import datetime
                    try:
                        dt = datetime.strptime(parts[0], '%Y-%m-%d')
                        breakpoint = dt.strftime('%b %Y')
                    except:
                        pass
                else:
                    breakpoint = f"{parts[0]} {parts[1]}"
        
        results.append({
            'Region': region,
            'Breakpoint': breakpoint,
            'PRE SPEI': format_2sig(row['PRE_mean_SPEI48']),
            'POST SPEI': format_2sig(row['POST_mean_SPEI48']),
            'ΔSPEI': format_delta_with_sign(row['Delta_SPEI48']),
            'Pettitt p-value': format_p_value(row['MWU_pvalue'])
        })
    
    df_table = pd.DataFrame(results)
    df_table.to_csv(TABLE_SX_OUTPUT, index=False, encoding='utf-8-sig')
    
    # Print to console
    print("\n" + "-" * 110)
    print("TABLE Sx OUTPUT:")
    print("-" * 110)
    print(df_table.to_string(index=False))
    print("-" * 110)
    print(f"✅ Saved: {TABLE_SX_OUTPUT}")
    
    return df_table


# =============================================================================
# PRINT ALL TABLES SUMMARY
# =============================================================================

def print_all_tables_summary(table_x, table_y, table_sx):
    """Print a compact summary of all three tables for quick cross-checking."""
    
    print("\n" + "=" * 110)
    print("QUICK CROSS-CHECK SUMMARY (All Tables)")
    print("=" * 110)
    
    if table_x is not None:
        print("\n📊 TABLE X - Regional Spatial WUE_T Decline Frequency:")
        print("   (Long-term persistence, NO breakpoint)")
        print("-" * 90)
        for _, row in table_x.iterrows():
            region = row['Region']
            mean_freq = row['Mean frequency']
            pct_ge_50 = row['Pixels >=0.50 (%)']
            print(f"   {region:<30} Mean = {mean_freq:<8}  >=0.50% = {pct_ge_50}")
    
    if table_y is not None:
        print("\n📊 TABLE Y - Regional PRE vs POST Breakpoint Changes:")
        print("   (Domain-wide Sept 2014 breakpoint)")
        print("-" * 90)
        for _, row in table_y.iterrows():
            region = row['Region']
            pre = row['PRE mean']
            post = row['POST mean']
            delta = row['Δ (POST-PRE)']
            print(f"   {region:<30} PRE = {pre:<8} POST = {post:<8} Δ = {delta}")
    
    if table_sx is not None:
        print("\n📊 TABLE Sx - Breakpoint Timing and SPEI Shifts:")
        print("   (Region-specific SPEI-derived breakpoints)")
        print("-" * 90)
        for _, row in table_sx.iterrows():
            region = row['Region']
            bp = row['Breakpoint']
            delta_spei = row['ΔSPEI']
            print(f"   {region:<30} Breakpoint = {bp:<14} ΔSPEI = {delta_spei}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    print("=" * 80)
    print("Q4 TABLE GENERATION: Spatial Frequency, PRE/POST, SPEI Shifts")
    print("=" * 80)
    print("\n📁 Output directory:", OUTPUT_DIR)
    
    # Check input files
    print("\n" + "-" * 40)
    print("INPUT FILE STATUS")
    print("-" * 40)
    
    files = {
        "Spatial diagnostic (Table X)": SPATIAL_DIAGNOSTIC_FILE,
        "PRE/POST diagnostic (Table Y)": PREPOST_DIAGNOSTIC_FILE,
        "Breakpoint SPEI (Table Sx)": BREAKPOINT_SPEI_FILE
    }
    
    for name, path in files.items():
        if os.path.exists(path):
            print(f"  ✅ {name}: Found")
        else:
            print(f"  ❌ {name}: NOT FOUND")
    
    # Generate all three tables
    table_x = generate_table_x()
    table_y = generate_table_y()
    table_sx = generate_table_sx()
    
    # Print quick summary
    print_all_tables_summary(table_x, table_y, table_sx)
    
    print("\n" + "=" * 80)
    print("TABLE GENERATION COMPLETE")
    print("=" * 80)
    print(f"\n📁 Output files saved to: {OUTPUT_DIR}")
    print("   • Table_X_Regional_Spatial_WUET_Decline_Frequency.csv")
    print("   • Table_Y_Regional_PRE_POST_Breakpoint_Changes.csv")
    print("   • Table_Sx_Breakpoint_Timing_SPEI_Shifts.csv")
    print("\n💡 Column headers use '>=' and '<=' instead of unicode symbols")
    print("   to avoid encoding issues in CSV files.")


if __name__ == "__main__":
    main()