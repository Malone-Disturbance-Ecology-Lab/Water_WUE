"""
CHUNK 1: Data Loading and Analysis (Shared between Panel A and B)
- Loads monthly data and Q3 final dataset
- Computes site-level SPEI statistics (min, max, median)
- Computes p-values and significance flags
- Applies global capping for visualization
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kruskal
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# SETUP OUTPUT DIRECTORY
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CHUNK 1: Data Loading and Analysis")
print("="*70)

# ============================================
# STEP 1: LOAD DATA
# ============================================

print("\n📂 Loading data...")

# Load the Q3 final dataset (has slopes and coastal regions)
input_file = output_dir / 'Q3_final_dataset.csv'
df = pd.read_csv(input_file)

# Load monthly data to compute SPEI statistics and p-values
monthly_filtered_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
nn_medians_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"

df_monthly = pd.read_csv(monthly_filtered_path)
nn_data = pd.read_csv(nn_medians_path)

print(f"  Loaded Q3 data: {len(df)} sites")
print(f"  Loaded monthly data: {len(df_monthly):,} rows")

# ============================================
# STEP 2: GET STRICT TRIPLE INTERSECTION SITES
# ============================================

print("\n📊 Getting strict triple intersection sites...")

wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"  Strict triple intersection sites: {len(shared_sites)}")

# ============================================
# STEP 3: COMPUTE SITE-LEVEL SPEI STATISTICS
# ============================================

print("\n📊 Computing site-level SPEI statistics...")

spei_summary_df = (
    df_monthly
    .dropna(subset=["SPEI_6", "WUE_tra"])
    .groupby("site_name")
    .agg(
        spei6_min=("SPEI_6", "min"),
        spei6_max=("SPEI_6", "max"),
        spei6_median=("SPEI_6", "median"),
        spei6_range=("SPEI_6", lambda x: x.max() - x.min())
    )
    .reset_index()
)

print(f"  Computed SPEI statistics for {len(spei_summary_df)} sites")
print(f"    SPEI min: {spei_summary_df['spei6_min'].min():.2f} to {spei_summary_df['spei6_min'].max():.2f}")
print(f"    SPEI max: {spei_summary_df['spei6_max'].min():.2f} to {spei_summary_df['spei6_max'].max():.2f}")
print(f"    SPEI median: {spei_summary_df['spei6_median'].min():.2f} to {spei_summary_df['spei6_median'].max():.2f}")

# Merge SPEI statistics into main dataframe
df = df.merge(spei_summary_df, on='site_name', how='left')
print(f"  Merged SPEI data: {len(df)} sites")

# ============================================
# STEP 4: COMPUTE P-VALUES AND SIGNIFICANCE
# ============================================

print("\n📊 Computing p-values and significance flags...")

site_pvalues = {}
for site in df['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site].dropna(subset=['SPEI_6', 'WUE_tra'])
    if len(site_data) >= 5:
        X = site_data['SPEI_6'].values
        y = site_data['WUE_tra'].values
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)
            site_pvalues[site] = p_value
        except:
            site_pvalues[site] = np.nan
    else:
        site_pvalues[site] = np.nan

df['p_value'] = df['site_name'].map(site_pvalues)
df['is_significant'] = (df['p_value'] < 0.05) & (df['p_value'].notna())

# Drop missing slopes
df = df.dropna(subset=['slope_theilsen'])
print(f"  Total sites: {len(df)}")
print(f"  Significant sites (p < 0.05): {df['is_significant'].sum()} / {len(df)} ({df['is_significant'].sum()/len(df)*100:.1f}%)")

# ============================================
# STEP 5: APPLY GLOBAL CAPPING
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
# STEP 6: COMPUTE COAST ORDER FOR PANEL B
# ============================================

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
    print(f"    {coast_full_names.get(coast, coast)}: {median_val:.4f}")

# ============================================
# STEP 7: SETUP COLOR SCALE (DIVERGING, CENTERED AT 0)
# ============================================

print("\n🎨 Setting up diverging color scale (centered at 0)...")

from matplotlib.colors import TwoSlopeNorm

# Use SPEI median as the color variable (hydroclimatic tendency)
color_var = "spei6_median"

# Diverging normalization centered at 0
norm_spei = TwoSlopeNorm(
    vmin=df["spei6_min"].min(),
    vcenter=0,
    vmax=df["spei6_max"].max()
)

# Use RdBu (Red for negative/dry, Blue for positive/wet)
# RdBu: Red = negative values, Blue = positive values
custom_cmap = plt.cm.RdBu

print(f"  Color variable: {color_var}")
print(f"  SPEI min: {df['spei6_min'].min():.2f}")
print(f"  SPEI max: {df['spei6_max'].max():.2f}")
print(f"  Center at: 0")
print(f"  Colormap: RdBu (RED for negative/dry SPEI, BLUE for positive/wet SPEI)")

# Calculate y-axis range for star placement
y_range = global_pos_95 - global_neg_95
y_margin = y_range * 0.05

# Save processed data for panels
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
    'color_var': color_var
}

print("\n" + "="*70)
print("CHUNK 1 COMPLETE - Data ready for panels")
print("="*70)

# Save processed data to pickle for panels to use
import pickle
pickle_file = output_dir / 'processed_data.pkl'
with open(pickle_file, 'wb') as f:
    pickle.dump(processed_data, f)
print(f"\n✅ Saved processed data to: {pickle_file}")