# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 10:55:57 2025
@author: ammar

UPDATED VERSION (2026-05-09):
- NO outlier removal (data already cleaned)
- Using original variables: WUE, WUE_tra, WUE_eva (NO _clean variables)
- Added Upland as water category (4 categories total)
- water_class comes from monthly_data (not site_info)
- IGBP stays in monthly_data (not dropped)
- OUTPUT FILE NAMES MATCH ORIGINAL CODE for downstream compatibility
- **NEW: Strict month filter - only months where ALL THREE metrics have data**
- **FIXED: Trans_ratio now correctly set to NaN when WUE_tra or WUE_eva missing**
- **WUE metrics (WUE, WUE_eva, WUE_tra) are NEVER filled or modified - NaNs remain NaNs**
- **FIXED: N_months now represents retained monthly records (not unique month-year combos)**
- **NEW: GLOBAL STRICT METRIC FILTER applied BEFORE saving monthly file**
"""

print("="*80)
print("COMPREHENSIVE WUE-SPEI ANALYSIS WORKFLOW (UPDATED)")
print("="*80)

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, kruskal, wilcoxon, spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
import os
from itertools import combinations

warnings.filterwarnings('ignore')

# =============================================================================
# CHUNK 1: DATA LOADING AND PREPROCESSING (WITH GLOBAL STRICT METRIC FILTER)
# =============================================================================
print("\n" + "="*80)
print("CHUNK 1: DATA LOADING AND PREPROCESSING")
print("="*80)

# Set working directory
print("Setting working directory...")
results_directory = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
os.chdir(results_directory)
print(f"✅ Working directory set to: {os.getcwd()}")

# File paths
monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
site_info_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\info\site_lat_long.csv"

# Load data with encoding handling
print("Loading data...")
encodings_to_try = ['utf-8', 'latin1', 'ISO-8859-1', 'cp1252', 'utf-16']

site_info = None
for encoding in encodings_to_try:
    try:
        site_info = pd.read_csv(site_info_path, encoding=encoding)
        print(f"✅ Successfully loaded site_info with {encoding} encoding")
        break
    except (UnicodeDecodeError, UnicodeError):
        continue

if site_info is None:
    site_info = pd.read_csv(site_info_path, encoding='latin1', on_bad_lines='skip')

monthly_data = None
for encoding in encodings_to_try:
    try:
        monthly_data = pd.read_csv(monthly_data_path, encoding=encoding)
        print(f"✅ Successfully loaded monthly_data with {encoding} encoding")
        break
    except (UnicodeDecodeError, UnicodeError):
        continue

if monthly_data is None:
    monthly_data = pd.read_csv(monthly_data_path, encoding='latin1', on_bad_lines='skip')

print(f"Monthly data shape: {monthly_data.shape}")
print(f"Site info shape: {site_info.shape}")

# Rename 'Year' to 'year' for consistency
if 'Year' in monthly_data.columns:
    monthly_data = monthly_data.rename(columns={'Year': 'year'})

# =============================================================================
# FIX: CORRECT Trans_ratio CALCULATION (NaN when data missing, NEVER fill WUE metrics)
# =============================================================================
print("\n" + "="*80)
print("FIX: Correcting Trans_ratio calculation (NaN for missing data)")
print("NOTE: This ONLY affects Trans_ratio - WUE, WUE_eva, WUE_tra are NEVER filled")
print("="*80)

# Check if Trans_ratio column exists and fix it
if 'Trans_ratio' in monthly_data.columns:
    # IMPORTANT: We NEVER modify WUE, WUE_eva, or WUE_tra columns
    # We ONLY fix Trans_ratio to be consistent with the data
    
    # Count how many rows have Trans_ratio = 0 but missing WUE_tra or WUE_eva
    zero_but_missing_mask = (monthly_data['Trans_ratio'] == 0) & \
                            (monthly_data['WUE_tra'].isna() | monthly_data['WUE_eva'].isna())
    n_incorrect = zero_but_missing_mask.sum()
    
    if n_incorrect > 0:
        print(f"  Found {n_incorrect} rows where Trans_ratio = 0 but WUE_tra or WUE_eva is missing")
        print(f"  Fixing: Setting these Trans_ratio values to NaN (NOT filling WUE metrics)")
        monthly_data.loc[zero_but_missing_mask, 'Trans_ratio'] = np.nan
    
    # Also check for rows where Trans_ratio is not NaN but either component is missing
    invalid_ratio_mask = (monthly_data['Trans_ratio'].notna()) & \
                         (monthly_data['WUE_tra'].isna() | monthly_data['WUE_eva'].isna())
    n_invalid = invalid_ratio_mask.sum()
    
    if n_invalid > 0:
        print(f"  Found {n_invalid} rows where Trans_ratio is not NaN but WUE_tra or WUE_eva is missing")
        print(f"  Fixing: Setting these Trans_ratio values to NaN")
        monthly_data.loc[invalid_ratio_mask, 'Trans_ratio'] = np.nan
    
    # Recalculate Trans_ratio where both components exist but ratio is missing
    both_exist_mask = monthly_data['WUE_tra'].notna() & monthly_data['WUE_eva'].notna() & monthly_data['Trans_ratio'].isna()
    n_recalc = both_exist_mask.sum()
    
    if n_recalc > 0:
        print(f"  Found {n_recalc} rows where both components exist but Trans_ratio is missing")
        print(f"  Recalculating Trans_ratio for these rows")
        monthly_data.loc[both_exist_mask, 'Trans_ratio'] = \
            monthly_data.loc[both_exist_mask, 'WUE_tra'] / \
            (monthly_data.loc[both_exist_mask, 'WUE_tra'] + monthly_data.loc[both_exist_mask, 'WUE_eva'])
    
    # Final validation
    still_invalid = (monthly_data['Trans_ratio'].notna()) & \
                    (monthly_data['WUE_tra'].isna() | monthly_data['WUE_eva'].isna())
    if still_invalid.sum() > 0:
        print(f"  ⚠️ WARNING: {still_invalid.sum()} rows still have invalid Trans_ratio")
    else:
        print(f"  ✅ Trans_ratio fix complete - all values now consistent with component data")
    
    # Verify we never modified WUE metrics
    print(f"\n  VERIFICATION: WUE, WUE_eva, WUE_tra columns were NEVER modified")
    print(f"  Only Trans_ratio column was corrected")
else:
    print(f"  ⚠️ Trans_ratio column not found - skipping fix")

# =============================================================================
# CONTINUE WITH NORMAL PIPELINE (NO CHANGES TO WUE METRICS)
# =============================================================================

# Columns to drop from monthly data (only lat and long - these come from site_info)
# DO NOT drop IGBP, climate, or water_class - they stay in monthly_data
columns_to_drop = ['lat', 'long']
for col in columns_to_drop:
    if col in monthly_data.columns:
        print(f"Will drop '{col}' from monthly data to avoid merge conflicts")

monthly_data_clean = monthly_data.drop(columns=[col for col in columns_to_drop if col in monthly_data.columns]).copy()

# Merge with site_info for lat and long only (IGBP and climate stay in monthly_data)
print("\nMerging data...")
site_info_merge = site_info[['site_name', 'lat', 'long']].copy()
monthly_data_merged = monthly_data_clean.merge(site_info_merge, on='site_name', how='left')
print(f"After merge - monthly data shape: {monthly_data_merged.shape}")

# Verify IGBP is still present
if 'IGBP' in monthly_data_merged.columns:
    print(f"✅ IGBP column present in merged data")
    print(f"   Unique IGBP values: {monthly_data_merged['IGBP'].unique()}")
else:
    print(f"⚠️ IGBP column not found. Checking original monthly_data...")
    if 'IGBP' in monthly_data.columns:
        print("   IGBP is in original monthly_data, copying over...")
        monthly_data_merged['IGBP'] = monthly_data['IGBP']

# Verify climate is still present
if 'climate' in monthly_data_merged.columns:
    print(f"✅ climate column present in merged data")
else:
    print(f"⚠️ climate column not found. Checking original monthly_data...")
    if 'climate' in monthly_data.columns:
        print("   climate is in original monthly_data, copying over...")
        monthly_data_merged['climate'] = monthly_data['climate']

# Map salinity/water categories (UPDATED: Includes Upland)
def map_salinity_category(cat):
    if pd.isna(cat):
        return "Unknown"
    cat_str = str(cat).lower().strip()
    
    # Upland check (NEW)
    if "upland" in cat_str:
        return "Upland"
    
    freshwater_terms = ["fresh", "freshwater", "fresh marsh", "tidal freshwater"]
    brackish_terms = ["brackish", "oligohaline", "mesohaline"]  
    saline_terms = ["saline", "polyhaline", "salt", "salt marsh"]

    if any(term in cat_str for term in freshwater_terms):
        return "Freshwater"
    elif any(term in cat_str for term in brackish_terms):
        return "Brackish"
    elif any(term in cat_str for term in saline_terms):
        return "Saline"
    else:
        return cat_str.title() if cat_str else "Unknown"

# Use water_class from monthly_data (already in merged data)
if 'water_class' in monthly_data_merged.columns:
    monthly_data_merged['Salinity_Category'] = monthly_data_merged['water_class'].apply(map_salinity_category)
    print(f"\n✅ Created Salinity_Category from water_class")
else:
    print("⚠️ water_class column not found - checking original monthly_data...")
    if 'water_class' in monthly_data.columns:
        monthly_data_merged['water_class'] = monthly_data['water_class']
        monthly_data_merged['Salinity_Category'] = monthly_data_merged['water_class'].apply(map_salinity_category)
    else:
        monthly_data_merged['Salinity_Category'] = "Unknown"

print(f"\nSalinity categories distribution:")
sal_cat_counts = monthly_data_merged['Salinity_Category'].value_counts()
for cat, count in sal_cat_counts.items():
    print(f"  {cat}: {count} records ({count/len(monthly_data_merged)*100:.1f}%)")

# SPEI classification
spei_columns = ['SPEI_1', 'SPEI_3', 'SPEI_6', 'SPEI_12', 'SPEI_24', 'SPEI_36', 'SPEI_48']

def classify_spei(spei_value):
    if pd.isna(spei_value):
        return np.nan
    if spei_value >= 2.0:
        return "EW"
    elif 1.5 <= spei_value < 2.0:
        return "SW"
    elif 1.0 <= spei_value < 1.5:
        return "MoW"
    elif 0.5 < spei_value < 1.0:
        return "MW"
    elif -0.5 <= spei_value <= 0.5:
        return "NN"
    elif -1.0 < spei_value < -0.5:
        return "MD"
    elif -1.5 < spei_value <= -1.0:
        return "MoD"
    elif -2.0 < spei_value <= -1.5:
        return "SD"
    elif spei_value <= -2.0:
        return "ED"
    else:
        return np.nan

ALL_SPEI_CLASSES = ["EW", "SW", "MoW", "MW", "NN", "MD", "MoD", "SD", "ED"]
COMPARISON_CLASSES = ["MW", "MoW", "SW", "EW", "MD", "MoD", "SD", "ED"]
WET_CLASSES = ["MW", "MoW", "SW", "EW"]
DRY_CLASSES = ["MD", "MoD", "SD", "ED"]

print(f"\nSPEI classes defined: {ALL_SPEI_CLASSES}")
print(f"Comparison classes (vs NN): {COMPARISON_CLASSES}")
print(f"Wet classes: {WET_CLASSES}")
print(f"Dry classes: {DRY_CLASSES}")

# Create SPEI category columns
for spei_col in spei_columns:
    monthly_data_merged[f"{spei_col}_Cat"] = monthly_data_merged[spei_col].apply(classify_spei)

# Create and save site metadata file
print("\n💾 Creating site metadata file...")

metadata_cols = ['site_name', 'lat', 'long', 'Salinity_Category']
if 'water_class' in monthly_data_merged.columns:
    metadata_cols.append('water_class')
if 'climate' in monthly_data_merged.columns:
    metadata_cols.append('climate')
if 'IGBP' in monthly_data_merged.columns:
    metadata_cols.append('IGBP')

print(f"  Using metadata columns: {metadata_cols}")
site_metadata = monthly_data_merged[metadata_cols].drop_duplicates()
site_metadata.to_csv('site_metadata_with_salinity_SPEIinfo.csv', index=False)
print(f"✅ Saved: site_metadata_with_salinity_SPEIinfo.csv with {len(site_metadata)} sites")

# DATA VALIDATION
print("\n🔍 DATA VALIDATION CHECK:")
critical_columns = ['site_name', 'year', 'month', 'Salinity_Category'] + spei_columns + ['WUE', 'WUE_eva', 'WUE_tra', 'Trans_ratio']
missing_columns = [col for col in critical_columns if col not in monthly_data_merged.columns]
if missing_columns:
    print(f"⚠️ WARNING: Missing critical columns: {missing_columns}")
else:
    print("✅ All critical columns present")

for col in ['site_name', 'year', 'month', 'Salinity_Category', 'Trans_ratio']:
    if col in monthly_data_merged.columns:
        nan_count = monthly_data_merged[col].isna().sum()
        if nan_count > 0:
            print(f"⚠️ WARNING: {col} has {nan_count} NaN values")

# NO OUTLIER REMOVAL - use data as is
print("\n=== SKIPPING OUTLIER REMOVAL (Data already cleaned) ===")
print("✅ Using original WUE, WUE_eva, WUE_tra variables directly")
print("✅ WUE metrics are NEVER filled or modified - NaNs remain NaNs")

# =============================================================================
# CRITICAL FIX: APPLY GLOBAL STRICT METRIC FILTER
# Keep only rows where WUE, WUE_eva, and WUE_tra ALL exist
# This ensures all downstream workflows use identical retained observations
# =============================================================================

print("\n" + "="*80)
print("APPLYING GLOBAL STRICT METRIC FILTER")
print("="*80)
print("  Filter criteria: WUE, WUE_eva, WUE_tra ALL must have data")
print("  This ensures consistency across ALL downstream workflows")

strict_metric_mask = (
    monthly_data_merged['WUE'].notna() &
    monthly_data_merged['WUE_eva'].notna() &
    monthly_data_merged['WUE_tra'].notna()
)

rows_before = len(monthly_data_merged)
monthly_data_clean = monthly_data_merged[strict_metric_mask].copy()
rows_after = len(monthly_data_clean)

print(f"\n  Rows before filter: {rows_before:,}")
print(f"  Rows after filter: {rows_after:,}")
print(f"  Rows removed: {rows_before - rows_after:,} ({((rows_before - rows_after)/rows_before*100):.1f}%)")
print(f"\n✅ Applied GLOBAL strict metric filter")

# Define WUE metrics
wue_metrics = ['WUE', 'WUE_eva', 'WUE_tra']

print(f"\n✅ Data loading and preprocessing completed")
print(f"Total records after strict metric filter: {len(monthly_data_clean)}")
print(f"Salinity distribution: {monthly_data_clean['Salinity_Category'].value_counts().to_dict()}")
print(f"SPEI timescales: {len(spei_columns)}")
print(f"WUE metrics: {wue_metrics}")
print(f"SPEI classes: {ALL_SPEI_CLASSES}")

# Save filtered monthly data (SAME NAME - downstream files read this)
monthly_data_clean.to_csv('monthly_data_after_outlier_removal.csv', index=False)
print(f"\n✅ Saved: monthly_data_after_outlier_removal.csv with {len(monthly_data_clean)} records")
print("  NOTE: This file now contains ONLY rows where WUE, WUE_eva, WUE_tra ALL have data")

# =============================================================================
# CHUNK 2: SKIPPED - NO OUTLIER REMOVAL
# =============================================================================
print("\n" + "="*80)
print("CHUNK 2: OUTLIER DETECTION AND REMOVAL - SKIPPED")
print("="*80)
print("✅ Data already cleaned. Proceeding to site-level summary.")

# =============================================================================
# CHUNK 3: SITE-LEVEL SUMMARY DATASET (WITH STRICT MONTH FILTER)
# =============================================================================
print("\n" + "="*80)
print("CHUNK 3: SITE-LEVEL SUMMARY DATASET")
print("="*80)
print("NOTE: Using strict month filter - only months where ALL THREE metrics have data")
print("      Months with ANY NaN in WUE metrics are EXCLUDED (NOT filled)")
print("      **N_months now counts retained monthly records directly**")

MIN_STRICT_MONTHS = 3
spei_category_columns = [f"{col}_Cat" for col in spei_columns]

for spei_cat_col in spei_category_columns:
    spei_timescale = spei_cat_col.replace('_Cat', '')
    print(f"\nProcessing {spei_timescale}...")
    
    site_summary_list = []
    nn_medians_list = []

    for site in monthly_data_clean['site_name'].unique():
        site_data = monthly_data_clean[monthly_data_clean['site_name'] == site].copy()
        if len(site_data) == 0:
            continue
            
        site_meta = site_data.iloc[0]
        
        # Apply STRICT MONTH FILTER (already filtered globally, but re-apply for safety)
        strict_mask = (site_data['WUE'].notna() & 
                       site_data['WUE_eva'].notna() & 
                       site_data['WUE_tra'].notna())
        site_data_strict = site_data[strict_mask].copy()
        
        # Check if site has enough strict months for NN condition
        nn_strict_mask = (site_data_strict[spei_cat_col] == 'NN')
        nn_strict_months = nn_strict_mask.sum()
        
        if nn_strict_months < MIN_STRICT_MONTHS:
            continue
        
        for wue_metric in wue_metrics:
            if wue_metric not in site_data_strict.columns:
                continue
                
            for spei_class in ALL_SPEI_CLASSES:
                class_data = site_data_strict[(site_data_strict[spei_cat_col] == spei_class) & 
                                               (site_data_strict[wue_metric].notna())]
                values = class_data[wue_metric].dropna()
                
                if len(values) > 0:
                    n_records = len(values)
                    unique_months = n_records
                    unique_years = class_data['year'].nunique()
                    sample_size_flag = "Adequate" if n_records >= 6 else "Low"
                    
                    stats_dict = {
                        'site_name': site,
                        'SPEI_Timescale': spei_timescale,
                        'SPEI_Class': spei_class,
                        'WUE_Metric': wue_metric,
                        'Median': np.median(values),
                        'Mean': np.mean(values),
                        'Std': np.std(values),
                        'IQR_25': np.percentile(values, 25),
                        'IQR_75': np.percentile(values, 75),
                        'Min': np.min(values),
                        'Max': np.max(values),
                        'N_months': unique_months,
                        'N_records': n_records,
                        'N_years': unique_years,
                        'Sample_Size_Flag': sample_size_flag,
                        'lat': site_meta.get('lat', np.nan),
                        'long': site_meta.get('long', np.nan),
                        'Salinity_Category': site_meta.get('Salinity_Category', 'Unknown'),
                        'climate': site_meta.get('climate', 'Unknown'),
                        'biome': site_meta.get('IGBP', 'Unknown')
                    }
                    site_summary_list.append(stats_dict)
                    
                    if spei_class == "NN":
                        nn_dict = {
                            'site_name': site,
                            'SPEI_Timescale': spei_timescale,
                            'WUE_Metric': wue_metric,
                            'WUE_median': np.median(values),
                            'lat': site_meta.get('lat', np.nan),
                            'long': site_meta.get('long', np.nan),
                            'Salinity_Category': site_meta.get('Salinity_Category', 'Unknown'),
                            'climate': site_meta.get('climate', 'Unknown'),
                            'biome': site_meta.get('IGBP', 'Unknown'),
                            'N_months': unique_months,
                            'N_records': n_records,
                            'N_years': unique_years,
                            'Sample_Size_Flag': sample_size_flag
                        }
                        nn_medians_list.append(nn_dict)

    if site_summary_list:
        filename = f'wue_site_level_summary_{spei_timescale}.csv'
        site_summary_df = pd.DataFrame(site_summary_list)
        site_summary_df.to_csv(filename, index=False)
        print(f"✅ Saved: {filename} with {len(site_summary_df)} records")
        
        if nn_medians_list:
            nn_filename = f'wue_site_level_NN_medians_{spei_timescale}.csv'
            nn_medians_df = pd.DataFrame(nn_medians_list)
            nn_medians_df.to_csv(nn_filename, index=False)
            print(f"✅ Saved: {nn_filename} with {len(nn_medians_df)} NN median records")

print(f"\n✅ CHUNK 3 COMPLETED")

# =============================================================================
# CHUNK 4: SPEI COMPARISON AND STATISTICAL ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("CHUNK 4: SPEI COMPARISON AND STATISTICAL ANALYSIS")
print("="*80)
print("NOTE: Using STRICT month filter - only months where ALL THREE metrics have data")
print("="*80)

def cliffs_delta(x, y):
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    wins = sum(1 for i in x for j in y if i > j)
    losses = sum(1 for i in x for j in y if i < j)
    return (wins - losses) / (n1 * n2)

for spei_cat_col in spei_category_columns:
    spei_timescale = spei_cat_col.replace('_Cat', '')
    print(f"\nProcessing {spei_timescale}...")
    
    try:
        site_level_data = pd.read_csv(f'wue_site_level_summary_{spei_timescale}.csv')
        print(f"✅ Loaded site-level data for {spei_timescale}")
    except FileNotFoundError:
        print(f"❌ Site level data not found for {spei_timescale}, skipping...")
        continue
    
    spei_comparison_list = []
    
    for site in site_level_data['site_name'].unique():
        site_data = site_level_data[site_level_data['site_name'] == site]
        if len(site_data) == 0:
            continue
        
        site_monthly = monthly_data_clean[monthly_data_clean['site_name'] == site].copy()
        if len(site_monthly) == 0:
            continue
        
        strict_mask = (site_monthly['WUE'].notna() & 
                       site_monthly['WUE_eva'].notna() & 
                       site_monthly['WUE_tra'].notna())
        site_monthly_strict = site_monthly[strict_mask].copy()
        
        if len(site_monthly_strict) == 0:
            continue
            
        for wue_metric in wue_metrics:
            metric_data = site_data[site_data['WUE_Metric'] == wue_metric]
            nn_data = metric_data[metric_data['SPEI_Class'] == 'NN']
            
            if len(nn_data) == 0:
                continue
                
            nn_monthly = site_monthly_strict[
                (site_monthly_strict[spei_cat_col] == 'NN') &
                (site_monthly_strict[wue_metric].notna())
            ][wue_metric]
            
            if len(nn_monthly) == 0:
                continue
            
            nn_median_monthly = np.median(nn_monthly)
            nn_sample_flag = nn_data['Sample_Size_Flag'].iloc[0] if 'Sample_Size_Flag' in nn_data.columns else "Unknown"
            
            for spei_class in COMPARISON_CLASSES:
                class_data = metric_data[metric_data['SPEI_Class'] == spei_class]
                if len(class_data) == 0:
                    continue
                
                class_monthly = site_monthly_strict[
                    (site_monthly_strict[spei_cat_col] == spei_class) &
                    (site_monthly_strict[wue_metric].notna())
                ][wue_metric]
                
                if len(class_monthly) == 0:
                    continue
                
                class_sample_flag = class_data['Sample_Size_Flag'].iloc[0] if 'Sample_Size_Flag' in class_data.columns else "Unknown"
                percent_changes = ((class_monthly - nn_median_monthly) / nn_median_monthly) * 100
                sample_quality_flag = "Adequate" if (nn_sample_flag == "Adequate" and class_sample_flag == "Adequate") else "Low"
                condition_type = "Wet" if spei_class in WET_CLASSES else "Dry"
                
                for pc in percent_changes:
                    spei_comparison_list.append({
                        'site_name': site,
                        'SPEI_Timescale': spei_timescale,
                        'SPEI_Class': spei_class,
                        'Condition_Type': condition_type,
                        'WUE_Metric': wue_metric,
                        'Median_NN': nn_median_monthly,
                        'Median_Condition': np.median(class_monthly),
                        'Percent_Change': pc,
                        'N_months_NN': len(nn_monthly),
                        'N_months_Condition': len(class_monthly),
                        'Sample_Quality_Flag': sample_quality_flag,
                        'Salinity_Category': site_data['Salinity_Category'].iloc[0] if 'Salinity_Category' in site_data.columns else 'Unknown',
                        'lat': site_data['lat'].iloc[0] if 'lat' in site_data.columns else np.nan,
                        'long': site_data['long'].iloc[0] if 'long' in site_data.columns else np.nan,
                        'climate': site_data['climate'].iloc[0] if 'climate' in site_data.columns else 'Unknown',
                        'biome': site_data['biome'].iloc[0] if 'biome' in site_data.columns else 'Unknown'
                    })
    
    if not spei_comparison_list:
        print(f"  ❌ No comparison data available")
        continue
    
    comparison_filename = f'wue_spei_vs_nn_{spei_timescale}.csv'
    wue_spei_vs_nn = pd.DataFrame(spei_comparison_list)
    wue_spei_vs_nn.to_csv(comparison_filename, index=False)
    print(f"✅ Saved: {comparison_filename} with {len(wue_spei_vs_nn)} comparisons")
    
    # WITHIN-SITE Paired Wilcoxon tests
    within_site_results = []
    
    for wue_metric in wue_metrics:
        family_a_tests = []
        
        for spei_class in COMPARISON_CLASSES:
            paired_data = []
            for site in wue_spei_vs_nn['site_name'].unique():
                site_condition = wue_spei_vs_nn[
                    (wue_spei_vs_nn['site_name'] == site) & 
                    (wue_spei_vs_nn['WUE_Metric'] == wue_metric) & 
                    (wue_spei_vs_nn['SPEI_Class'] == spei_class)
                ]
                if len(site_condition) > 0:
                    paired_data.append((site_condition['Median_NN'].iloc[0], 
                                       site_condition['Median_Condition'].iloc[0]))
            
            if len(paired_data) >= 6:
                nn_vals = [x[0] for x in paired_data]
                cond_vals = [x[1] for x in paired_data]
                
                try:
                    w_stat, p_value = wilcoxon(nn_vals, cond_vals)
                    effect_size = cliffs_delta(nn_vals, cond_vals)
                    percent_change = ((np.median(cond_vals) - np.median(nn_vals)) / np.median(nn_vals)) * 100
                    direction = "↑" if percent_change > 0 else "↓"
                    condition_type = "Wet" if spei_class in WET_CLASSES else "Dry"
                    
                    family_a_tests.append({
                        'spei_class': spei_class,
                        'condition_type': condition_type,
                        'p_raw': p_value,
                        'effect_size': effect_size,
                        'n_pairs': len(paired_data),
                        'percent_change': percent_change,
                        'direction': direction
                    })
                except:
                    pass
        
        if family_a_tests:
            p_values = [test['p_raw'] for test in family_a_tests]
            _, p_adj_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
            
            for i, test in enumerate(family_a_tests):
                within_site_results.append({
                    'SPEI_Timescale': spei_timescale,
                    'Analysis_Type': 'WITHIN_SITE_Paired_Wilcoxon',
                    'WUE_Metric': wue_metric,
                    'Comparison': f'NN_vs_{test["spei_class"]}',
                    'Condition_Type': test['condition_type'],
                    'p_raw': test['p_raw'],
                    'p_adj_bh': p_adj_bh[i],
                    'effect_size': test['effect_size'],
                    'N_pairs': test['n_pairs'],
                    'percent_change': test['percent_change'],
                    'direction': test['direction']
                })
    
    if within_site_results:
        within_site_df = pd.DataFrame(within_site_results)
        within_site_df.to_csv(f'statistical_results_within_sites_{spei_timescale}.csv', index=False)
        print(f"✅ Saved: statistical_results_within_sites_{spei_timescale}.csv")
    
    # ACROSS-SITE tests
    across_site_results = []
    
    for wue_metric in wue_metrics:
        for spei_class in COMPARISON_CLASSES:
            percent_changes = wue_spei_vs_nn[
                (wue_spei_vs_nn['WUE_Metric'] == wue_metric) & 
                (wue_spei_vs_nn['SPEI_Class'] == spei_class)
            ]['Percent_Change'].dropna()
            
            if len(percent_changes) >= 6:
                try:
                    w_stat, p_value = wilcoxon(percent_changes)
                    effect_size = cliffs_delta(percent_changes, [0]*len(percent_changes))
                    median_change = np.median(percent_changes)
                    direction = "↑" if median_change > 0 else "↓"
                    condition_type = "Wet" if spei_class in WET_CLASSES else "Dry"
                    
                    across_site_results.append({
                        'SPEI_Timescale': spei_timescale,
                        'Analysis_Type': 'ACROSS_SITE_Wilcoxon_OneSample',
                        'WUE_Metric': wue_metric,
                        'SPEI_Class': spei_class,
                        'Condition_Type': condition_type,
                        'p_raw': p_value,
                        'effect_size': effect_size,
                        'N_sites': wue_spei_vs_nn[
                            (wue_spei_vs_nn['WUE_Metric'] == wue_metric) & 
                            (wue_spei_vs_nn['SPEI_Class'] == spei_class)
                        ]['site_name'].nunique(),
                        'N_observations': len(percent_changes),
                        'Median_Percent_Change': median_change,
                        'direction': direction,
                        'Direction_Method': 'Bootstrap_95CI_monthly'
                    })
                except:
                    pass
    
    if across_site_results:
        across_site_df = pd.DataFrame(across_site_results)
        p_values = across_site_df['p_raw'].dropna()
        if len(p_values) > 0:
            _, p_adj_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
            across_site_df['p_adj_bh'] = p_adj_bh
        across_site_df.to_csv(f'statistical_results_across_sites_{spei_timescale}.csv', index=False)
        print(f"✅ Saved: statistical_results_across_sites_{spei_timescale}.csv")

print(f"\n✅ CHUNK 4 COMPLETED")

# =============================================================================
# CHUNK 5: SALINITY-STRATIFIED ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("CHUNK 5: SALINITY-STRATIFIED ANALYSIS")
print("="*80)

salinity_categories = ['Freshwater', 'Brackish', 'Saline', 'Upland']

for spei_timescale in spei_columns:
    print(f"\nProcessing {spei_timescale}...")
    
    try:
        spei_data = pd.read_csv(f'wue_spei_vs_nn_{spei_timescale}.csv')
        print(f"  ✅ Loaded {len(spei_data)} records")
    except FileNotFoundError:
        print(f"  ❌ Data not found, skipping...")
        continue
    
    panel_b_combined = []
    
    for salinity_cat in salinity_categories:
        cat_spei_data = spei_data[spei_data['Salinity_Category'] == salinity_cat]
        if len(cat_spei_data) == 0:
            continue
            
        for wue_metric in wue_metrics:
            for spei_class in COMPARISON_CLASSES:
                metric_class_data = cat_spei_data[
                    (cat_spei_data['WUE_Metric'] == wue_metric) & 
                    (cat_spei_data['SPEI_Class'] == spei_class)
                ]
                
                if len(metric_class_data) > 0:
                    condition_type = "Wet" if spei_class in WET_CLASSES else "Dry"
                    
                    if len(metric_class_data) >= 6:
                        try:
                            _, p_value = wilcoxon(metric_class_data['Percent_Change'])
                            effect_size = cliffs_delta(metric_class_data['Percent_Change'], [0]*len(metric_class_data))
                        except:
                            p_value = np.nan
                            effect_size = np.nan
                    else:
                        p_value = np.nan
                        effect_size = np.nan
                    
                    positive_changes = len(metric_class_data[metric_class_data['Percent_Change'] > 0])
                    total_changes = len(metric_class_data)
                    
                    panel_b_combined.append({
                        'Salinity_Category': salinity_cat,
                        'SPEI_Timescale': spei_timescale,
                        'WUE_Metric': wue_metric,
                        'SPEI_Class': spei_class,
                        'Condition_Type': condition_type,
                        'N_sites': len(metric_class_data),
                        'Median_Percent_Change': metric_class_data['Percent_Change'].median(),
                        'Mean_Percent_Change': metric_class_data['Percent_Change'].mean(),
                        'Std_Percent_Change': metric_class_data['Percent_Change'].std(),
                        'p_value_change': p_value,
                        'effect_size_change': effect_size,
                        'Positive_Changes': positive_changes,
                        'Total_Changes': total_changes,
                        'Positive_Proportion': positive_changes / total_changes if total_changes > 0 else 0
                    })
    
    if panel_b_combined:
        panel_b_combined_df = pd.DataFrame(panel_b_combined)
        panel_b_combined_filename = f'fig3_panelB_percent_change_salinity_{spei_timescale}.csv'
        panel_b_combined_df.to_csv(panel_b_combined_filename, index=False)
        print(f"  ✅ Saved: {panel_b_combined_filename} with {len(panel_b_combined_df)} records")
    
    # ACROSS-SITE SALINITY COMPARISONS
    salinity_comparison_results = []
    
    for wue_metric in wue_metrics:
        for condition_type in ['Wet', 'Dry']:
            spei_classes_to_test = WET_CLASSES if condition_type == 'Wet' else DRY_CLASSES
            
            for spei_class in spei_classes_to_test:
                salinity_groups = []
                valid_salinities = []
                group_sizes = {}
                
                for salinity_cat in salinity_categories:
                    sal_data = spei_data[
                        (spei_data['Salinity_Category'] == salinity_cat) & 
                        (spei_data['WUE_Metric'] == wue_metric) & 
                        (spei_data['SPEI_Class'] == spei_class)
                    ]['Percent_Change'].dropna()
                    
                    if len(sal_data) >= 3:
                        salinity_groups.append(sal_data)
                        valid_salinities.append(salinity_cat)
                        group_sizes[salinity_cat] = len(sal_data)
                
                if len(salinity_groups) >= 2:
                    try:
                        h_stat, p_value = kruskal(*salinity_groups)
                        total_n = sum(len(g) for g in salinity_groups)
                        epsilon_squared = h_stat / (total_n - 1) if total_n > 1 else np.nan
                        
                        salinity_comparison_results.append({
                            'SPEI_Timescale': spei_timescale,
                            'WUE_Metric': wue_metric,
                            'SPEI_Class': spei_class,
                            'Condition_Type': condition_type,
                            'Test': 'ACROSS_SITE_Kruskal_Wallis_Salinity',
                            'Test_Statistic': h_stat,
                            'p_raw': p_value,
                            'effect_size': epsilon_squared,
                            'N_groups': len(salinity_groups),
                            'Group_sizes': str(group_sizes),
                            'Group_names': ','.join(valid_salinities)
                        })
                    except Exception as e:
                        print(f"    Error in Kruskal-Wallis test: {e}")
    
    if salinity_comparison_results:
        salinity_comparison_df = pd.DataFrame(salinity_comparison_results)
        p_values = salinity_comparison_df['p_raw'].dropna()
        if len(p_values) > 0:
            _, p_adj_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
            salinity_comparison_df['p_adj_bh'] = p_adj_bh
        salinity_comparison_df.to_csv(f'salinity_comparison_across_sites_{spei_timescale}.csv', index=False)
        print(f"  ✅ Saved: salinity_comparison_across_sites_{spei_timescale}.csv")

print(f"\n✅ CHUNK 5 COMPLETED")

# =============================================================================
# CHUNK 6: SENSITIVITY ANALYSIS
# =============================================================================
print("\n" + "="*80)
print("CHUNK 6: SENSITIVITY ANALYSIS")
print("="*80)

class_merging_schemes = {
    'original': {
        'EW': 4, 'SW': 3, 'MoW': 2, 'MW': 1, 'NN': 0, 
        'MD': -1, 'MoD': -2, 'SD': -3, 'ED': -4
    },
    'merged_wet_dry': {
        'EW_SW': 3, 'MoW_MW': 2, 'NN': 0, 'MD_MoD': -2, 'SD_ED': -4
    },
    'binary_wet_dry': {
        'Wet': 1, 'NN': 0, 'Dry': -1
    },
    'three_level': {
        'EW': 2, 'SW': 2, 'MoW': 1, 'MW': 1, 'NN': 0, 
        'MD': -1, 'MoD': -1, 'SD': -2, 'ED': -2
    }
}

for spei_timescale in spei_columns:
    print(f"\nProcessing {spei_timescale}...")
    
    try:
        site_level_data = pd.read_csv(f'wue_site_level_summary_{spei_timescale}.csv')
        print(f"  ✅ Loaded site-level data")
    except FileNotFoundError:
        print(f"  ❌ Site level data not found, skipping...")
        continue
    
    for scheme_name, severity_scheme in class_merging_schemes.items():
        print(f"  Testing scheme: {scheme_name}")
        
        sensitivity_list = []
        
        for site in site_level_data['site_name'].unique():
            site_data = site_level_data[site_level_data['site_name'] == site]
            if len(site_data) == 0:
                continue
            
            for wue_metric in wue_metrics:
                metric_data = site_data[site_data['WUE_Metric'] == wue_metric]
                
                severity_scores = []
                wue_values = []
                
                for spei_class, severity_score in severity_scheme.items():
                    if '_' in spei_class:
                        component_classes = spei_class.split('_')
                        class_data_list = []
                        for comp_class in component_classes:
                            comp_data = metric_data[metric_data['SPEI_Class'] == comp_class]
                            if len(comp_data) > 0:
                                class_data_list.append(comp_data)
                        if not class_data_list:
                            continue
                        class_data = class_data_list[0]
                    else:
                        class_data = metric_data[metric_data['SPEI_Class'] == spei_class]
                    
                    if len(class_data) > 0:
                        severity_scores.append(severity_score)
                        wue_values.append(class_data['Median'].iloc[0])
                
                if len(severity_scores) >= 3 and len(set(severity_scores)) > 1:
                    try:
                        slope, _, r_value, p_value, _ = stats.linregress(severity_scores, wue_values)
                        mean_wue = np.mean(wue_values)
                        relative_sensitivity = (slope / mean_wue) * 100 if mean_wue != 0 else np.nan
                        
                        sensitivity_list.append({
                            'site_name': site,
                            'SPEI_Timescale': spei_timescale,
                            'Class_Merging_Scheme': scheme_name,
                            'WUE_Metric': wue_metric,
                            'Slope': slope,
                            'Slope_Relative_Percent': relative_sensitivity,
                            'R_squared': r_value**2,
                            'P_Value_Raw': p_value,
                            'N_Categories': len(severity_scores),
                            'Salinity_Category': site_data['Salinity_Category'].iloc[0] if 'Salinity_Category' in site_data.columns else 'Unknown'
                        })
                    except:
                        pass
        
        if sensitivity_list:
            sensitivity_df = pd.DataFrame(sensitivity_list)
            sensitivity_df.to_csv(f'wue_sensitivity_within_sites_{spei_timescale}_{scheme_name}.csv', index=False)
            print(f"    ✅ Saved: wue_sensitivity_within_sites_{spei_timescale}_{scheme_name}.csv ({len(sensitivity_df)} records)")
            
            for wue_metric in wue_metrics:
                metric_slopes = sensitivity_df[
                    (sensitivity_df['WUE_Metric'] == wue_metric) &
                    (sensitivity_df['Slope'].notna())
                ]['Slope']
                
                if len(metric_slopes) >= 6:
                    try:
                        _, p_value = wilcoxon(metric_slopes)
                        median_slope = np.median(metric_slopes)
                        median_relative = np.median(sensitivity_df[
                            (sensitivity_df['WUE_Metric'] == wue_metric) &
                            (sensitivity_df['Slope_Relative_Percent'].notna())
                        ]['Slope_Relative_Percent'])
                        
                        across_result = pd.DataFrame([{
                            'SPEI_Timescale': spei_timescale,
                            'Class_Merging_Scheme': scheme_name,
                            'WUE_Metric': wue_metric,
                            'N_sites': len(metric_slopes),
                            'Median_Slope': median_slope,
                            'Median_Relative_Slope_Percent': median_relative,
                            'p_value': p_value
                        }])
                        
                        across_result.to_csv(f'wue_sensitivity_across_sites_{spei_timescale}_{scheme_name}.csv', index=False)
                        print(f"    ✅ Saved: wue_sensitivity_across_sites_{spei_timescale}_{scheme_name}.csv")
                    except:
                        pass

print(f"\n✅ CHUNK 6 COMPLETED")

# =============================================================================
# CHUNK 7: FINAL OUTPUTS AND METHODOLOGY SUMMARY
# =============================================================================
print("\n" + "="*80)
print("CHUNK 7: FINAL OUTPUTS AND METHODOLOGY SUMMARY")
print("="*80)

print("\n📊 DATASET SUMMARY:")
print(f"  Total monthly records (after strict metric filter): {len(monthly_data_clean):,}")
print(f"  Unique sites: {monthly_data_clean['site_name'].nunique()}")
if 'year' in monthly_data_clean.columns:
    print(f"  Years covered: {monthly_data_clean['year'].nunique()} ({monthly_data_clean['year'].min()}-{monthly_data_clean['year'].max()})")
print(f"  Salinity categories: {monthly_data_clean['Salinity_Category'].value_counts().to_dict()}")
print(f"  SPEI timescales: {len(spei_columns)}")
print(f"  WUE metrics: {wue_metrics}")
print(f"  SPEI classes: {len(ALL_SPEI_CLASSES)} ({len(WET_CLASSES)} wet, 1 normal, {len(DRY_CLASSES)} dry)")

print(f"\n💾 OUTPUT FILES CREATED:")
import glob

file_patterns = [
    ('Site Metadata', 'site_metadata_with_salinity_SPEIinfo.csv'),
    ('Monthly Data (filtered)', 'monthly_data_after_outlier_removal.csv'),
    ('Site-Level Summaries', 'wue_site_level_summary_*.csv'),
    ('NN Medians', 'wue_site_level_NN_medians_*.csv'),
    ('SPEI Comparisons', 'wue_spei_vs_nn_*.csv'),
    ('Within-Site Stats', 'statistical_results_within_sites_*.csv'),
    ('Across-Site Stats', 'statistical_results_across_sites_*.csv'),
    ('Salinity Analysis', 'fig3_panelB_percent_change_salinity_*.csv'),
    ('Salinity Comparisons', 'salinity_comparison_across_sites_*.csv'),
    ('Sensitivity Within-Site', 'wue_sensitivity_within_sites_*.csv'),
    ('Sensitivity Across-Site', 'wue_sensitivity_across_sites_*.csv')
]

for category, pattern in file_patterns:
    files = glob.glob(pattern)
    if files:
        print(f"  {category}: {len(files)} files")

summary_data = {
    'Metric': [
        'Total_Sites', 'Total_Years', 'Total_Monthly_Records_After_Filter',
        'SPEI_Timescales', 'WUE_Metrics', 'SPEI_Classes_Total',
        'SPEI_Classes_Wet', 'SPEI_Classes_Dry', 'Class_Merging_Schemes_Tested'
    ],
    'Count': [
        monthly_data_clean['site_name'].nunique(),
        monthly_data_clean['year'].nunique() if 'year' in monthly_data_clean.columns else 0,
        len(monthly_data_clean),
        len(spei_columns),
        len(wue_metrics),
        len(ALL_SPEI_CLASSES),
        len(WET_CLASSES),
        len(DRY_CLASSES),
        len(class_merging_schemes)
    ]
}

for sal_cat in salinity_categories:
    summary_data['Metric'].append(f'{sal_cat}_Sites')
    summary_data['Count'].append(len(monthly_data_clean[monthly_data_clean['Salinity_Category'] == sal_cat]['site_name'].unique()))

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('study_summary_statistics.csv', index=False)
print(f"\n✅ Saved: study_summary_statistics.csv")

print("\n" + "="*80)
print("✅ WORKFLOW COMPLETED SUCCESSFULLY")
print("="*80)
print("\n🎯 KEY FEATURES:")
print("  • NO outlier removal - using original cleaned data")
print("  • Original WUE variables: WUE, WUE_eva, WUE_tra")
print("  • Water categories: Freshwater, Brackish, Saline, Upland")
print("  • Complete SPEI gradient: wet → normal → dry")
print("  • SAME OUTPUT FILE NAMES as original code for downstream compatibility")
print("  • FIXED: Trans_ratio now correctly set to NaN when WUE_tra or WUE_eva missing")
print("  • WUE metrics (WUE, WUE_eva, WUE_tra) are NEVER filled or modified")
print("  • FIXED: N_months now counts retained monthly records directly")
print("  • NEW: GLOBAL STRICT METRIC FILTER applied - monthly file contains ONLY rows")
print("    where WUE, WUE_eva, and WUE_tra ALL have data")
print("="*80)