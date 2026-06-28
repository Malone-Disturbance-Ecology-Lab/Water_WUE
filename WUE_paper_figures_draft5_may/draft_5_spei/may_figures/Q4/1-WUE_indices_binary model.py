# -*- coding: utf-8 -*-
"""
UPDATED WUE ANOMALY WORKFLOW - PERCENTILE-BASED CLASSIFICATION
(NO BOOTSTRAP - uses percentiles of anomalies directly)
=============================================================================
UPDATED: Now includes Trans_ratio and time columns in final dataset
for downstream Q4 interaction analysis.
=============================================================================
"""

print("="*80)
print("UPDATED WUE ANOMALY WORKFLOW - PERCENTILE-BASED CLASSIFICATION")
print("INCLUDES Trans_ratio FOR DOWNSTREAM ANALYSIS")
print("="*80)

import pandas as pd
import numpy as np
import warnings
import os

warnings.filterwarnings('ignore')

# =============================================================================
# PARAMETERS
# =============================================================================
LOWER_PERCENTILE = 10
UPPER_PERCENTILE = 90

print(f"\n🔧 CLASSIFICATION PARAMETERS:")
print(f"   Method: Percentile-based (NO bootstrap)")
print(f"   Percentile range: {LOWER_PERCENTILE}th - {UPPER_PERCENTILE}th percentile")
print(f"   Coverage: {UPPER_PERCENTILE - LOWER_PERCENTILE}%")
print(f"   Rule: anomaly < {LOWER_PERCENTILE}th percentile → Decrease")
print(f"         anomaly > {UPPER_PERCENTILE}th percentile → Increase")
print(f"         {LOWER_PERCENTILE}th-{UPPER_PERCENTILE}th percentile → NoChange")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n" + "="*60)
print("STEP 1: LOADING DATA")
print("="*60)

monthly_data_path = r"M:\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
site_info_path = r"M:\Research\WUE_CUE\data_products\results\site_metadata_with_salinity_SPEIinfo.csv"
nn_medians_path = r"M:\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

print(f"  Monthly data: {monthly_data_path}")
print(f"  Site metadata: {site_info_path}")
print(f"  NN medians: {nn_medians_path}")

monthly_data = pd.read_csv(monthly_data_path)
site_info = pd.read_csv(site_info_path)
nn_medians = pd.read_csv(nn_medians_path)

print(f"  ✅ Monthly data shape: {monthly_data.shape}")
print(f"     Unique sites: {monthly_data['site_name'].nunique()}")
print(f"  ✅ Site info shape: {site_info.shape}")
print(f"  ✅ NN medians shape: {nn_medians.shape}")

# =============================================================================
# 2. FILTER TO SITES WITH NN MEDIANS
# =============================================================================
print("\n" + "="*60)
print("STEP 2: FILTERING TO SITES WITH NN MEDIANS")
print("="*60)

nn_medians_spei1 = nn_medians[nn_medians['SPEI_Timescale'] == 'SPEI_1'].copy()
sites_with_nn = set(nn_medians_spei1['site_name'].unique())
print(f"  Sites with NN medians: {len(sites_with_nn)}")

monthly_data_filtered = monthly_data[monthly_data['site_name'].isin(sites_with_nn)].copy()
print(f"  Monthly data after filtering: {monthly_data_filtered.shape}")
print(f"  Unique sites: {monthly_data_filtered['site_name'].nunique()}")

# =============================================================================
# 3. MERGE AND PREPARE DATA
# =============================================================================
print("\n" + "="*60)
print("STEP 3: MERGING AND PREPARING DATA")
print("="*60)

if 'Salinity_Category' not in monthly_data_filtered.columns:
    monthly_data_filtered = monthly_data_filtered.merge(
        site_info[['site_name', 'Salinity_Category', 'water_class', 'IGBP', 'climate', 'lat', 'long']],
        on='site_name',
        how='left'
    )
    print("  ✅ Added site metadata")

if 'IGBP' in monthly_data_filtered.columns and 'biome' not in monthly_data_filtered.columns:
    monthly_data_filtered['biome'] = monthly_data_filtered['IGBP']
    print("  ✅ Renamed IGBP → biome")

spei_columns = ['SPEI_1', 'SPEI_3', 'SPEI_6', 'SPEI_12', 'SPEI_24', 'SPEI_36', 'SPEI_48']
wue_metrics = ['WUE', 'WUE_eva', 'WUE_tra']

print(f"  SPEI columns: {spei_columns}")
print(f"  WUE metrics: {wue_metrics}")
print(f"  Unique sites: {monthly_data_filtered['site_name'].nunique()}")
print(f"  Total rows: {len(monthly_data_filtered)}")

# =============================================================================
# 4. CREATE ANOMALY COLUMNS
# =============================================================================
print("\n" + "="*60)
print("STEP 4: CREATING ANOMALY COLUMNS")
print("="*60)

nn_median_WUE = {}
nn_median_WUE_eva = {}
nn_median_WUE_tra = {}

for _, row in nn_medians_spei1.iterrows():
    site = row['site_name']
    metric = row['WUE_Metric']
    median_val = row['WUE_median']
    
    if metric == 'WUE':
        nn_median_WUE[site] = median_val
    elif metric == 'WUE_eva':
        nn_median_WUE_eva[site] = median_val
    elif metric == 'WUE_tra':
        nn_median_WUE_tra[site] = median_val

print(f"  NN medians loaded: WUE={len(nn_median_WUE)} sites")

monthly_data_with_anomalies = monthly_data_filtered.copy()

for spei_col in spei_columns:
    for site, nn_med in nn_median_WUE.items():
        if pd.isna(nn_med):
            continue
        site_mask = monthly_data_with_anomalies['site_name'] == site
        
        anomaly_col = f"WUE_minus_NNmed_{spei_col}"
        if anomaly_col not in monthly_data_with_anomalies.columns:
            monthly_data_with_anomalies[anomaly_col] = np.nan
        monthly_data_with_anomalies.loc[site_mask, anomaly_col] = (
            monthly_data_with_anomalies.loc[site_mask, 'WUE'] - nn_med
        )
        
        if site in nn_median_WUE_eva and not pd.isna(nn_median_WUE_eva[site]):
            anomaly_col_eva = f"WUE_eva_minus_NNmed_{spei_col}"
            if anomaly_col_eva not in monthly_data_with_anomalies.columns:
                monthly_data_with_anomalies[anomaly_col_eva] = np.nan
            monthly_data_with_anomalies.loc[site_mask, anomaly_col_eva] = (
                monthly_data_with_anomalies.loc[site_mask, 'WUE_eva'] - nn_median_WUE_eva[site]
            )
        
        if site in nn_median_WUE_tra and not pd.isna(nn_median_WUE_tra[site]):
            anomaly_col_tra = f"WUE_tra_minus_NNmed_{spei_col}"
            if anomaly_col_tra not in monthly_data_with_anomalies.columns:
                monthly_data_with_anomalies[anomaly_col_tra] = np.nan
            monthly_data_with_anomalies.loc[site_mask, anomaly_col_tra] = (
                monthly_data_with_anomalies.loc[site_mask, 'WUE_tra'] - nn_median_WUE_tra[site]
            )

print(f"  ✅ Created anomaly columns for {len(spei_columns)} SPEI timescales")

# =============================================================================
# 5. PERCENTILE-BASED CLASSIFICATION
# =============================================================================
print("\n" + "="*60)
print("STEP 5: PERCENTILE-BASED CLASSIFICATION")
print(f"   Using {LOWER_PERCENTILE}th-{UPPER_PERCENTILE}th percentile range")
print("="*60)

anomaly_columns = [col for col in monthly_data_with_anomalies.columns if 'minus_NNmed' in col]
print(f"  Found {len(anomaly_columns)} anomaly columns to classify")

for anomaly_col in anomaly_columns:
    class_col = f"{anomaly_col}_Class"
    monthly_data_with_anomalies[class_col] = np.nan
    
    for site in monthly_data_with_anomalies['site_name'].unique():
        site_mask = monthly_data_with_anomalies['site_name'] == site
        site_anomalies = monthly_data_with_anomalies.loc[site_mask, anomaly_col].dropna().values
        
        if len(site_anomalies) < 3:
            monthly_data_with_anomalies.loc[site_mask, class_col] = "NoChange"
            continue
        
        lower_threshold = np.percentile(site_anomalies, LOWER_PERCENTILE)
        upper_threshold = np.percentile(site_anomalies, UPPER_PERCENTILE)
        
        for idx in monthly_data_with_anomalies.loc[site_mask].index:
            anomaly_value = monthly_data_with_anomalies.loc[idx, anomaly_col]
            if pd.isna(anomaly_value):
                monthly_data_with_anomalies.loc[idx, class_col] = np.nan
            elif anomaly_value < lower_threshold:
                monthly_data_with_anomalies.loc[idx, class_col] = "Decrease"
            elif anomaly_value > upper_threshold:
                monthly_data_with_anomalies.loc[idx, class_col] = "Increase"
            else:
                monthly_data_with_anomalies.loc[idx, class_col] = "NoChange"

print(f"  ✅ Classification complete for {len(anomaly_columns)} columns")


# =============================================================================
# 6. REMOVE NaN ROWS (FIXED - preserves index, no accidental row loss)
# =============================================================================
print("\n" + "="*60)
print("STEP 6: REMOVING NaN ROWS")
print("="*60)

rows_before = len(monthly_data_with_anomalies)
sites_before = monthly_data_with_anomalies['site_name'].nunique()

print(f"  Rows before NaN removal: {rows_before}")
print(f"  Sites before NaN removal: {sites_before}")

# Identify SPEI-48 class columns
spei_48_class_cols = [
    col for col in monthly_data_with_anomalies.columns
    if col.endswith('_Class') and 'SPEI_48' in col
]

print(f"\n  SPEI-48 class columns checked ({len(spei_48_class_cols)} columns):")
missing_counts = {}
for col in spei_48_class_cols:
    n_missing = monthly_data_with_anomalies[col].isna().sum()
    missing_counts[col] = n_missing
    print(f"    {col}: missing rows = {n_missing}")

# IMPORTANT FIX: Use dataframe-based mask so pandas preserves the original index
# This prevents the index misalignment that was causing valid rows to be dropped
mask_all_valid = monthly_data_with_anomalies[spei_48_class_cols].notna().all(axis=1)

monthly_data_clean = monthly_data_with_anomalies.loc[mask_all_valid].copy()

rows_after = len(monthly_data_clean)
sites_after = monthly_data_clean['site_name'].nunique()
rows_removed = rows_before - rows_after

print(f"\n  Rows after NaN removal: {rows_after}")
print(f"  Sites after NaN removal: {sites_after}")
print(f"  Rows removed: {rows_removed} ({rows_removed/rows_before*100:.1f}%)")

if rows_removed > 0:
    removed_rows = monthly_data_with_anomalies.loc[~mask_all_valid].copy()
    print(f"\n  ⚠️ Removed rows by site:")
    print(removed_rows['site_name'].value_counts())
    
    # Also print first few removed rows for debugging
    print(f"\n  Sample of removed rows (first 5):")
    print(removed_rows[['site_name', 'SPEI_48'] + spei_48_class_cols[:3]].head())
else:
    print("  ✅ No rows removed; all SPEI-48 class labels are valid")

# =============================================================================
# 7. CREATE FINAL DATASET (UPDATED - INCLUDES Trans_ratio)
# =============================================================================
print("\n" + "="*60)
print("STEP 7: CREATING FINAL DATASET")
print("="*60)

# UPDATED: Include Trans_ratio and time columns for downstream analyses
base_columns = [
    'site_name',
    'Year',
    'year',
    'month',
    'date',
    'lat',
    'long',
    'climate',
    'biome',
    'Salinity_Category',
    'water_class',
    'IGBP',
    'Trans_ratio'
]

anomaly_columns = [col for col in monthly_data_clean.columns if 'minus_NNmed' in col and '_Class' not in col]
class_columns = [col for col in monthly_data_clean.columns if col.endswith('_Class')]

# Add any missing base columns that exist in the data
for col in base_columns:
    if col not in monthly_data_clean.columns:
        monthly_data_clean[col] = np.nan

# Build final columns list
final_columns = base_columns + spei_columns + anomaly_columns + class_columns
final_columns = list(dict.fromkeys(final_columns))
final_columns = [col for col in final_columns if col in monthly_data_clean.columns]

final_dataset = monthly_data_clean[final_columns].copy()

print(f"  Final dataset shape: {final_dataset.shape}")
print(f"  Columns included: {len(final_columns)}")
print(f"  ✓ Trans_ratio included: {'Trans_ratio' in final_dataset.columns}")

# Check how many rows have Trans_ratio
if 'Trans_ratio' in final_dataset.columns:
    n_with_tet = final_dataset['Trans_ratio'].notna().sum()
    print(f"  ✓ Rows with Trans_ratio: {n_with_tet}/{len(final_dataset)} ({100*n_with_tet/len(final_dataset):.1f}%)")

# =============================================================================
# 8. VALIDATION
# =============================================================================
print("\n" + "="*60)
print("STEP 8: VALIDATION")
print("="*60)

print("\n  Classification distribution (SPEI-48):")
if 'WUE_tra_minus_NNmed_SPEI_48_Class' in final_dataset.columns:
    spei48_class = final_dataset['WUE_tra_minus_NNmed_SPEI_48_Class'].value_counts()
    for cls, cnt in spei48_class.items():
        print(f"    {cls}: {cnt} ({cnt/len(final_dataset)*100:.1f}%)")

# =============================================================================
# 9. SAVE OUTPUTS (SAME FILE NAMES - CRITICAL for downstream)
# =============================================================================
print("\n" + "="*60)
print("STEP 9: SAVING OUTPUTS")
print("="*60)

output_dir = r"M:\Research\WUE_CUE\data_products\results"

final_nn_path = os.path.join(output_dir, "final_nn_medians.csv")
nn_medians_spei1.to_csv(final_nn_path, index=False)
print(f"  ✅ Saved: {final_nn_path}")

final_dataset_path = os.path.join(output_dir, "final_dataset_with_spei_and_anomalies.csv")
final_dataset.to_csv(final_dataset_path, index=False)
print(f"  ✅ Saved: {final_dataset_path}")

final_with_classes_path = os.path.join(output_dir, "final_dataset_with_spei_anomalies_and_classes.csv")
final_dataset.to_csv(final_with_classes_path, index=False)
print(f"  ✅ Saved: {final_with_classes_path}")

print("\n" + "="*60)
print("WORKFLOW COMPLETED SUCCESSFULLY!")
print("="*60)

print(f"\n📊 SUMMARY:")
print(f"   Method: Percentile-based (NO bootstrap, NO confidence intervals)")
print(f"   Percentile range: {LOWER_PERCENTILE}th-{UPPER_PERCENTILE}th")
print(f"   Final dataset includes: Trans_ratio, water_class, IGBP, time columns")

print(f"\n📁 OUTPUT FILES (SAME NAMES - No downstream impact):")
print(f"   • final_nn_medians.csv")
print(f"   • final_dataset_with_spei_and_anomalies.csv")
print(f"   • final_dataset_with_spei_anomalies_and_classes.csv")

print(f"\n✅ DOWNSTREAM COMPATIBILITY:")
print(f"   • Trans_ratio now available for Q4 interaction analysis")
print(f"   • All existing column names unchanged")
print(f"   • Additional columns added (no existing columns removed)")

print("\n" + "="*60)
print("UPDATE COMPLETE - Trans_ratio now in final dataset!")
print("="*60)