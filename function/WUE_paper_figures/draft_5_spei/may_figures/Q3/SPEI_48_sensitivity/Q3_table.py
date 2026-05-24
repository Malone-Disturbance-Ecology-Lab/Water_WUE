# -*- coding: utf-8 -*-
"""
Created on Mon May 18 23:12:03 2026

@author: ammar
"""
# -*- coding: utf-8 -*-
"""
TABLE Q3: Coastal Region Sensitivity to SPEI-48 and Hydroclimatic Exposure
Supplementary Table for Manuscript

Columns:
- Coastal region
- N sites
- Monthly observations (N)
- Median slope (2 decimal places)
- Slope IQR (Q1–Q3 with en dash, 2 decimal places)
- Significant slopes (count/total (percentage))
- Median SPEI-48 (2 decimal places)
- Dry months (%) - SPEI < -0.5
- Neutral months (%) - -0.5 ≤ SPEI ≤ 0.5
- Wet months (%) - SPEI > 0.5

@author: WUE Analysis Pipeline
Date: 2026-05-18
"""

import pandas as pd
import numpy as np
import pickle
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
pickle_file = output_dir / 'processed_data_SPEI48_filtered.pkl'
monthly_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"

# Output files
csv_output = output_dir / 'Table_Q3_CoastalSensitivity_Hydroclimate.csv'
excel_output = output_dir / 'Table_Q3_CoastalSensitivity_Hydroclimate.xlsx'


# =============================================================================
# HELPER FUNCTIONS FOR 2 DECIMAL PLACES (CONSISTENT)
# =============================================================================

def format_2decimals(value):
    """
    Format a number with EXACTLY 2 decimal places.
    Examples:
        0.020 -> 0.02
        0.045 -> 0.05
        0.066 -> 0.07
        0.152 -> 0.15
        0.629 -> 0.63
        -0.618 -> -0.62
        0.210 -> 0.21
        1.029 -> 1.03
    """
    if pd.isna(value):
        return '--'
    if value == 0:
        return '0.00'
    
    # Round to 2 decimal places
    rounded = round(value, 2)
    
    # Format with exactly 2 decimal places
    if rounded < 0:
        return f"{rounded:.2f}"
    else:
        return f"{rounded:.2f}"


def format_slope_2dec(value):
    """Format slope with exactly 2 decimal places"""
    if pd.isna(value):
        return '--'
    if value == 0:
        return '0.00'
    
    # Round to 2 decimal places
    rounded = round(value, 2)
    
    # Format with exactly 2 decimal places
    if rounded < 0:
        return f"{rounded:.2f}"
    else:
        return f"{rounded:.2f}"


def format_iqr_2dec_en_dash(q1, q3):
    """
    Format IQR range with exactly 2 decimal places using EN DASH.
    Example: -0.62 – 0.21
    """
    if pd.isna(q1) or pd.isna(q3):
        return '--'
    
    q1_rounded = round(q1, 2)
    q3_rounded = round(q3, 2)
    
    q1_str = f"{q1_rounded:.2f}" if q1_rounded < 0 else f"{q1_rounded:.2f}"
    q3_str = f"{q3_rounded:.2f}" if q3_rounded < 0 else f"{q3_rounded:.2f}"
    
    # Use EN DASH (U+2013) with spaces
    return f"{q1_str} – {q3_str}"


def format_significant_slopes(n_sig, n_sites):
    """Format as 'count/total (percentage%)' e.g., '7/12 (58%)'"""
    if n_sites == 0:
        return '--'
    pct = (n_sig / n_sites) * 100
    pct_rounded = int(round(pct))
    return f"{n_sig}/{n_sites} ({pct_rounded}%)"


def format_percent_1dec(value):
    """Format percentage with 1 decimal place"""
    if pd.isna(value):
        return '--'
    return f"{value:.1f}"


def format_median_spei_2dec(value):
    """Format median SPEI-48 with exactly 2 decimal places"""
    if pd.isna(value):
        return '--'
    if value == 0:
        return '0.00'
    
    # Round to 2 decimal places
    rounded = round(value, 2)
    
    # Format with exactly 2 decimal places
    if rounded < 0:
        return f"{rounded:.2f}"
    else:
        return f"{rounded:.2f}"


# =============================================================================
# LOAD DATA
# =============================================================================

print("=" * 80)
print("TABLE Q3: Coastal Region Sensitivity to SPEI-48")
print("=" * 80)

# Load processed sensitivity data
print(f"\n1. Loading sensitivity data from: {pickle_file}")
with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df = data['df']
coast_order_ascending = data['coast_order_ascending']
coast_full_names = data['coast_full_names']
coast_short_labels = data['coast_short_labels']

print(f"   Loaded {len(df)} site records")
print(f"   Coastal regions: {coast_order_ascending}")

# Load monthly data
print(f"\n2. Loading monthly data from: {monthly_file}")
df_monthly = pd.read_csv(monthly_file)
print(f"   Loaded {len(df_monthly):,} rows")

# =============================================================================
# APPLY FILTERS (SAME AS PANEL B AND C)
# =============================================================================

print("\n3. Applying filters...")

# Step 1: Restrict to valid sites only
valid_sites = df['site_name'].unique()
df_monthly = df_monthly[df_monthly['site_name'].isin(valid_sites)]
print(f"   After valid sites filter: {len(df_monthly):,} rows")

# Step 2: Apply Trans_ratio > 0 filter
df_monthly = df_monthly[df_monthly['Trans_ratio'] > 0]
print(f"   After Trans_ratio > 0: {len(df_monthly):,} rows")

# Step 3: Map coastal regions
site_to_coast = dict(zip(df['site_name'], df['coast_region_analysis']))
df_monthly['coast_region_analysis'] = df_monthly['site_name'].map(site_to_coast)

# Step 4: Drop rows without coast assignment
df_monthly = df_monthly.dropna(subset=['coast_region_analysis'])
print(f"   After coast mapping: {len(df_monthly):,} rows")

# =============================================================================
# COMPUTE STATISTICS PER COASTAL REGION
# =============================================================================

print("\n4. Computing statistics per coastal region...")

results = []

for coast in coast_order_ascending:
    coast_name = coast_full_names.get(coast, coast)
    coast_short = coast_short_labels.get(coast, coast)
    
    # ===== SENSITIVITY DATA (from df) =====
    coast_df = df[df['coast_region_analysis'] == coast]
    n_sites = len(coast_df)
    
    if n_sites == 0:
        print(f"   WARNING: No sites for {coast_name}")
        continue
    
    # Median slope (2 decimal places)
    median_slope = coast_df['slope_theilsen'].median()
    median_slope_formatted = format_slope_2dec(median_slope)
    
    # IQR with EN DASH (2 decimal places)
    q25 = coast_df['slope_theilsen'].quantile(0.25)
    q75 = coast_df['slope_theilsen'].quantile(0.75)
    iqr_formatted = format_iqr_2dec_en_dash(q25, q75)
    
    # Significant slopes: count/total (percentage)
    n_sig = coast_df['is_significant'].sum()
    sig_formatted = format_significant_slopes(n_sig, n_sites)
    
    # Median SPEI-48 (2 decimal places) - using median, not mean
    median_spei = coast_df['spei48_median'].median()
    median_spei_formatted = format_median_spei_2dec(median_spei)
    
    # ===== MONTHLY HYDROCLIMATE DATA (from df_monthly) =====
    coast_monthly = df_monthly[df_monthly['coast_region_analysis'] == coast]
    n_months = len(coast_monthly)
    
    if n_months > 0:
        spei_values = coast_monthly['SPEI_48'].dropna().values
        
        # Classify SPEI values
        dry_mask = spei_values < -0.5
        wet_mask = spei_values > 0.5
        neutral_mask = (spei_values >= -0.5) & (spei_values <= 0.5)
        
        pct_dry = (np.sum(dry_mask) / n_months) * 100
        pct_neutral = (np.sum(neutral_mask) / n_months) * 100
        pct_wet = (np.sum(wet_mask) / n_months) * 100
        
        dry_formatted = format_percent_1dec(pct_dry)
        neutral_formatted = format_percent_1dec(pct_neutral)
        wet_formatted = format_percent_1dec(pct_wet)
    else:
        n_months = 0
        dry_formatted = '--'
        neutral_formatted = '--'
        wet_formatted = '--'
    
    # Format monthly observations with comma
    n_months_formatted = f"{n_months:,}" if n_months > 0 else '0'
    
    results.append({
        'Coastal region': coast_name,
        'N sites': n_sites,
        'Monthly observations (N)': n_months_formatted,
        'Median slope': median_slope_formatted,
        'Slope IQR (Q1–Q3)': iqr_formatted,
        'Significant slopes': sig_formatted,
        'Median SPEI-48': median_spei_formatted,
        'Dry months (%)': dry_formatted,
        'Neutral months (%)': neutral_formatted,
        'Wet months (%)': wet_formatted
    })
    
    print(f"\n   {coast_name} ({coast_short}):")
    print(f"     N sites = {n_sites}")
    print(f"     Monthly observations = {n_months:,}")
    print(f"     Median slope = {median_slope_formatted} (IQR: {iqr_formatted})")
    print(f"     Significant slopes = {sig_formatted}")
    print(f"     Median SPEI-48 = {median_spei_formatted}")
    print(f"     Hydroclimate: Dry={dry_formatted}%, Neutral={neutral_formatted}%, Wet={wet_formatted}%")

# =============================================================================
# CREATE DATAFRAME AND SAVE
# =============================================================================

table_df = pd.DataFrame(results)

# Reorder columns as requested
column_order = [
    'Coastal region',
    'N sites',
    'Monthly observations (N)',
    'Median slope',
    'Slope IQR (Q1–Q3)',
    'Significant slopes',
    'Median SPEI-48',
    'Dry months (%)',
    'Neutral months (%)',
    'Wet months (%)'
]
table_df = table_df[column_order]

print("\n" + "=" * 80)
print("FINAL TABLE Q3")
print("=" * 80)
print("\n" + table_df.to_string(index=False))

# Save to CSV
table_df.to_csv(csv_output, index=False)
print(f"\n✅ Saved CSV: {csv_output}")

# Save to Excel with formatting
try:
    with pd.ExcelWriter(excel_output, engine='openpyxl') as writer:
        table_df.to_excel(writer, sheet_name='Table_Q3', index=False)
        
        # Auto-adjust column widths
        worksheet = writer.sheets['Table_Q3']
        for column in worksheet.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = min(max_length + 2, 35)
            worksheet.column_dimensions[column_letter].width = adjusted_width
    
    print(f"✅ Saved Excel: {excel_output}")
except Exception as e:
    print(f"⚠️ Could not save Excel file: {e}")

# =============================================================================
# PRINT MANUSCRIPT-READY TABLE
# =============================================================================

print("\n" + "=" * 80)
print("MANUSCRIPT-READY TABLE")
print("=" * 80)
print("\nSupplementary Table SX. Regional differences in WUE$_T$ sensitivity to SPEI-48 and hydroclimatic exposure")
print("\nValues: Median slope and IQR shown with 2 decimal places. IQR uses en dash (–).\n")

# Print formatted table
header = f"{'Coastal region':<20} {'N sites':<10} {'Monthly obs':<15} {'Median slope':<15} {'Slope IQR (Q1–Q3)':<22} {'Significant slopes':<20} {'Median SPEI-48':<15} {'Dry (%)':<10} {'Neutral (%)':<12} {'Wet (%)':<10}"
print(header)
print("-" * 150)

for _, row in table_df.iterrows():
    print(f"{row['Coastal region']:<20} {row['N sites']:<10} {row['Monthly observations (N)']:<15} {row['Median slope']:<15} {row['Slope IQR (Q1–Q3)']:<22} {row['Significant slopes']:<20} {row['Median SPEI-48']:<15} {row['Dry months (%)']:<10} {row['Neutral months (%)']:<12} {row['Wet months (%)']:<10}")

print("-" * 150)

print("\nNote: Dry = SPEI-48 < -0.5, Neutral = -0.5 to 0.5, Wet = SPEI-48 > 0.5")
print("      Sensitivity slopes represent ΔWUE$_T$ / ΔSPEI-48 from Theil-Sen regression.")
print("      Significant slopes: p < 0.05 from linear regression between WUE$_T$ and SPEI-48.")
print("      Median SPEI-48: site-level median of SPEI-48 values across the study period.")
print("      IQR shown as Q1–Q3 with en dash.")

print("\n" + "=" * 80)
print("TABLE Q3 COMPLETE")
print("=" * 80)
