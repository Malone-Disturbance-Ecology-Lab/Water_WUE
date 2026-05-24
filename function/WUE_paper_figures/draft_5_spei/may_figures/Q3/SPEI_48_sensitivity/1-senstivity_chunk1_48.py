# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:01:47 2026

@author: ammar
"""

"""
CHUNK 1: Data Loading and Analysis for SPEI-48 (with ratio filtering)
- Loads monthly data
- Filters for Trans_ratio > 0 (no NaN, no zero)
- Computes Theil-Sen slopes for SPEI-48 (NEW - not from old Q3 file)
- Computes site-level SPEI-48 statistics (min, max, median)
- Computes p-values via linear regression
- Applies global capping for visualization
- Prints detailed data loss tracking
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kruskal
from sklearn.linear_model import TheilSenRegressor
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# SETUP OUTPUT DIRECTORY
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CHUNK 1: Data Loading and Analysis (SPEI-48 with Ratio Filtering)")
print("="*70)

# ============================================
# STEP 1: LOAD DATA
# ============================================

print("\n📂 Loading data...")

# Load monthly data (NOT using Q3 file for slopes anymore)
monthly_filtered_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
nn_medians_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

df_monthly = pd.read_csv(monthly_filtered_path)
nn_data = pd.read_csv(nn_medians_path)

print(f"  Initial monthly data: {len(df_monthly):,} rows")

# ============================================
# STEP 2: GET STRICT TRIPLE INTERSECTION SITES
# ============================================

print("\n📊 Getting strict triple intersection sites...")

wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"  After triple intersection: {len(df_monthly):,} rows")
print(f"  Unique sites: {len(shared_sites)}")

# ============================================
# STEP 3: APPLY RATIO FILTERING (Trans_ratio > 0)
# ============================================

print("\n📊 Applying ratio filtering (Trans_ratio > 0, no NaN)...")

initial_rows = len(df_monthly)
initial_sites = df_monthly['site_name'].nunique()

# Filter out NaN or <= 0 Trans_ratio
df_monthly = df_monthly[df_monthly['Trans_ratio'] > 0].copy()

rows_after_ratio = len(df_monthly)
sites_after_ratio = df_monthly['site_name'].nunique()

print(f"  Rows removed: {initial_rows - rows_after_ratio} ({(1 - rows_after_ratio/initial_rows)*100:.1f}%)")
print(f"  Sites remaining: {sites_after_ratio} (from {initial_sites})")

# ============================================
# STEP 4: COMPUTE SITE-LEVEL SPEI-48 STATISTICS AND SLOPES
# ============================================

print("\n📊 Computing site-level SPEI-48 statistics and Theil-Sen slopes...")

site_stats = []
site_slopes = {}
site_pvalues = {}

for site in df_monthly['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site].dropna(subset=['SPEI_48', 'WUE_tra'])
    
    if len(site_data) >= 5:
        # Get SPEI-48 stats
        spei48_min = site_data['SPEI_48'].min()
        spei48_max = site_data['SPEI_48'].max()
        spei48_median = site_data['SPEI_48'].median()
        spei48_range = spei48_max - spei48_min
        n_months = len(site_data)
        
        # === THEIL-SEN SLOPE for SPEI-48 ===
        X = site_data['SPEI_48'].values.reshape(-1, 1)
        y = site_data['WUE_tra'].values
        
        try:
            ts = TheilSenRegressor(random_state=42).fit(X, y)
            slope_theilsen = ts.coef_[0]
        except:
            slope_theilsen = np.nan
        
        # === Linear regression p-value ===
        try:
            _, _, _, p_value, _ = stats.linregress(site_data['SPEI_48'].values, y)
        except:
            p_value = np.nan
        
        # Store results
        site_stats.append({
            'site_name': site,
            'spei48_min': spei48_min,
            'spei48_max': spei48_max,
            'spei48_median': spei48_median,
            'spei48_range': spei48_range,
            'n_months': n_months,
            'slope_theilsen': slope_theilsen,
            'p_value': p_value
        })
        site_slopes[site] = slope_theilsen
        site_pvalues[site] = p_value
    else:
        print(f"  Site {site} has only {len(site_data)} months (<5) - excluding")

# Create summary dataframe
spei_summary_df = pd.DataFrame(site_stats)

print(f"  Computed SPEI-48 statistics for {len(spei_summary_df)} sites")
print(f"    SPEI-48 min: {spei_summary_df['spei48_min'].min():.2f} to {spei_summary_df['spei48_min'].max():.2f}")
print(f"    SPEI-48 max: {spei_summary_df['spei48_max'].min():.2f} to {spei_summary_df['spei48_max'].max():.2f}")
print(f"    SPEI-48 median: {spei_summary_df['spei48_median'].min():.2f} to {spei_summary_df['spei48_median'].max():.2f}")
print(f"    Months per site: min={spei_summary_df['n_months'].min()}, max={spei_summary_df['n_months'].max()}, mean={spei_summary_df['n_months'].mean():.1f}")

# ============================================
# STEP 5: LOAD COASTAL REGION DATA (from Q3 file only for site metadata)
# ============================================

print("\n📊 Loading coastal region data for sites...")

# Load Q3 file ONLY for site metadata (coast_region_analysis)
q3_file = output_dir / 'Q3_final_dataset.csv'
if q3_file.exists():
    df_q3 = pd.read_csv(q3_file)
    # Keep only site metadata columns
    df_metadata = df_q3[['site_name', 'coast_region_analysis', 'lat', 'long', 'ecosystem_type']].drop_duplicates(subset=['site_name'])
    print(f"  Loaded metadata for {len(df_metadata)} sites")
else:
    print("  ⚠️ Q3_final_dataset.csv not found - creating empty metadata")
    df_metadata = pd.DataFrame(columns=['site_name', 'coast_region_analysis', 'lat', 'long', 'ecosystem_type'])

# Merge SPEI results with metadata
df = spei_summary_df.merge(df_metadata, on='site_name', how='left')

# ============================================
# STEP 6: ADD SIGNIFICANCE FLAGS
# ============================================

print("\n📊 Adding significance flags...")

df['is_significant'] = (df['p_value'] < 0.05) & (df['p_value'].notna())

# Drop sites with missing slopes
df = df.dropna(subset=['slope_theilsen'])
print(f"  Final sites after all filtering: {len(df)}")
print(f"  Significant sites (p < 0.05): {df['is_significant'].sum()} / {len(df)} ({df['is_significant'].sum()/len(df)*100:.1f}%)")

# ============================================
# STEP 7: PRINT DETAILED DATA LOSS SUMMARY
# ============================================

print("\n" + "="*70)
print("DATA LOSS SUMMARY")
print("="*70)
print(f"\n  Step 1 - Initial shared sites: {len(shared_sites)}")
print(f"  Step 2 - After Trans_ratio > 0 filter: {sites_after_ratio} sites retained")
print(f"  Step 3 - After min 5 obs & slope computation: {len(spei_summary_df)} sites")
print(f"  Step 4 - Final sites in analysis: {len(df)}")
print(f"\n  Total months in filtered dataset: {len(df_monthly)}")
print(f"  Date range in data: {df_monthly['year'].min()}-{df_monthly['month'].min()} to {df_monthly['year'].max()}-{df_monthly['month'].max()}")

# ============================================
# STEP 8: APPLY GLOBAL CAPPING FOR VISUALIZATION
# ============================================

print("\n📊 Applying global capping for visualization...")

global_pos_95 = df['slope_theilsen'].quantile(0.95)
global_neg_95 = df['slope_theilsen'].quantile(0.05)

print(f"  Lower bound (5th percentile): {global_neg_95:.4f}")
print(f"  Upper bound (95th percentile): {global_pos_95:.4f}")

df['slope_capped'] = df['slope_theilsen'].clip(lower=global_neg_95, upper=global_pos_95)
df['is_extreme_low'] = df['slope_theilsen'] < global_neg_95
df['is_extreme_high'] = df['slope_theilsen'] > global_pos_95

print(f"  Below lower bound: {df['is_extreme_low'].sum()} sites")
print(f"  Above upper bound: {df['is_extreme_high'].sum()} sites")

# Sort by capped slope for Panel A
df_sorted = df.sort_values('slope_capped').reset_index(drop=True)

# ============================================
# STEP 9: COMPUTE COAST ORDER FOR PANEL B
# ============================================

# Only compute coast order if coast_region_analysis exists and has values
if 'coast_region_analysis' in df.columns and df['coast_region_analysis'].notna().any():
    coast_medians = df.groupby('coast_region_analysis')['slope_theilsen'].median().sort_values()
    coast_order_ascending = coast_medians.index.tolist()
    
    coast_full_names = {
        'AK_coast': 'Alaska',
        'West_Coast': 'Pacific',
        'Gulf_of_America': 'Gulf',
        'Southeast_Atlantic': 'Southeast Atlantic',
        'Atlantic_Coast_North': 'Atlantic North'
    }
    
    coast_short_labels = {
        'AK_coast': 'AK',
        'West_Coast': 'Pacific',
        'Gulf_of_America': 'Gulf',
        'Southeast_Atlantic': 'SE Atl',
        'Atlantic_Coast_North': 'Atl N'
    }
    
    print(f"\n  Coast order based on median slope (ascending):")
    for coast in coast_order_ascending:
        median_val = coast_medians[coast]
        n_sites = len(df[df['coast_region_analysis'] == coast])
        print(f"    {coast_full_names.get(coast, coast)}: n={n_sites}, median={median_val:.4f}")
else:
    coast_order_ascending = []
    coast_short_labels = {}
    coast_full_names = {}
    print(f"\n  ⚠️ No coast region data available")

# ============================================
# STEP 10: SETUP COLOR SCALE (DIVERGING, CENTERED AT 0)
# ============================================

print("\n🎨 Setting up diverging color scale (centered at 0)...")

from matplotlib.colors import TwoSlopeNorm

# Use SPEI-48 median as the color variable
color_var = "spei48_median"

# Diverging normalization centered at 0
norm_spei = TwoSlopeNorm(
    vmin=df["spei48_min"].min(),
    vcenter=0,
    vmax=df["spei48_max"].max()
)

# Use RdBu (Red for negative/dry, Blue for positive/wet)
custom_cmap = plt.cm.RdBu

print(f"  Color variable: {color_var}")
print(f"  SPEI-48 min: {df['spei48_min'].min():.2f}")
print(f"  SPEI-48 max: {df['spei48_max'].max():.2f}")
print(f"  Center at: 0")
print(f"  Colormap: RdBu (RED for negative/dry SPEI, BLUE for positive/wet SPEI)")

# Calculate y-axis range for star placement
y_range = global_pos_95 - global_neg_95
y_margin = y_range * 0.05

# ============================================
# SAVE PROCESSED DATA
# ============================================

processed_data = {
    'df': df,
    'df_sorted': df_sorted,
    'coast_order_ascending': coast_order_ascending,
    'coast_short_labels': coast_short_labels,
    'coast_full_names': coast_full_names,
    'global_neg_95': global_neg_95,
    'global_pos_95': global_pos_95,
    'y_range': y_range,
    'y_margin': y_margin,
    'norm_spei': norm_spei,
    'custom_cmap': custom_cmap,
    'color_var': color_var,
    'spei_version': 'SPEI-48',
    'filtering_criteria': 'Trans_ratio > 0',
    'n_sites_initial': len(shared_sites),
    'n_sites_final': len(df),
    'n_months_total': len(df_monthly)
}

print("\n" + "="*70)
print("CHUNK 1 COMPLETE - Data ready for panels (SPEI-48 with ratio filtering)")
print("="*70)
print(f"\n  ✅ Theil-Sen slopes computed from SPEI-48 (not from old Q3 file)")
print(f"  ✅ Filtering criteria: Trans_ratio > 0 (no NaN, no zero)")
print(f"  ✅ Sites: {len(df)} retained from {len(shared_sites)} initial")
print(f"  ✅ Months: {len(df_monthly)} rows after filtering")
print(f"  ✅ Method: Theil-Sen slopes + linear regression p-values")
print(f"  ✅ SPEI metric: SPEI-48")

# Save processed data to pickle
import pickle
pickle_file = output_dir / 'processed_data_SPEI48_filtered.pkl'
with open(pickle_file, 'wb') as f:
    pickle.dump(processed_data, f)
print(f"\n✅ Saved processed data to: {pickle_file}")

# Print final summary for manuscript
print("\n" + "="*70)
print("MANUSCRIPT SUMMARY (SPEI-48 Analysis)")
print("="*70)
print(f"\n  Sample size: {len(df)} sites")
print(f"  Total months analyzed: {len(df_monthly)}")
print(f"  Temporal extent: {df_monthly['year'].min()}-{df_monthly['month'].min()} to {df_monthly['year'].max()}-{df_monthly['month'].max()}")
print(f"  Filter criteria: Trans_ratio > 0 (excluded invalid transpiration estimates)")
print(f"  Significance: {df['is_significant'].sum()} sites (p < 0.05, linear regression)")
print(f"  Slope method: Theil-Sen regression")

print("\n" + "="*70)
print("READY FOR CHUNK 2 (Panel A with SPEI-48)")
print("="*70)


"""
CREATE SPEI-48 SLOPES FILE
Extract slopes from processed SPEI-48 data for diagnostic analysis
"""

import pandas as pd
from pathlib import Path

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
pickle_file = output_dir / 'processed_data_SPEI48_filtered.pkl'

# Load processed data
import pickle
with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df = data['df']

# Create slopes dataframe
df_slopes = df[['site_name', 'slope_theilsen', 'p_value', 'is_significant',
                'spei48_min', 'spei48_max', 'spei48_median', 'spei48_range',
                'coast_region_analysis', 'n_months']].copy()

# Save to Q3_analysis folder
q3_analysis_dir = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\Q3_analysis")
q3_analysis_dir.mkdir(parents=True, exist_ok=True)

output_file = q3_analysis_dir / 'Q3_WUET_SPEI48_slopes.csv'
df_slopes.to_csv(output_file, index=False)

print(f"✅ Saved SPEI-48 slopes to: {output_file}")
print(f"   Total sites: {len(df_slopes)}")
print(f"   Columns: {df_slopes.columns.tolist()}")


import pandas as pd
from pathlib import Path

slopes_file = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\Q3_analysis\Q3_WUET_SPEI48_slopes.csv")

df = pd.read_csv(slopes_file)

print("Columns:", df.columns.tolist())

n_col = "n_months"

q25 = df[n_col].quantile(0.25)
median = df[n_col].median()
q75 = df[n_col].quantile(0.75)

print(f"N sites = {len(df)}")
print(f"Total months = {df[n_col].sum()}")
print(f"Min = {df[n_col].min()}")
print(f"Max = {df[n_col].max()}")
print(f"Mean = {df[n_col].mean():.1f}")
print(f"Median = {median:.0f}")
print(f"IQR = {q25:.0f}–{q75:.0f}")







