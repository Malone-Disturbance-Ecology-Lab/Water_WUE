# -*- coding: utf-8 -*-
"""
REVIEWER RESPONSE: Diagnostic Analysis - SPEI-48 Range vs WUE_T Sensitivity Magnitude
===============================================================================
Test whether site-level WUE_T sensitivity magnitude is related to the 
hydroclimatic variability (SPEI-48 range) experienced by each site.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# File paths
BASE_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results")
MONTHLY_FILE = BASE_DIR / "monthly_data_after_outlier_removal.csv"
NN_MEDIANS_FILE = BASE_DIR / "wue_site_level_NN_medians_SPEI_1.csv"

# Slopes file location
SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI48_slopes.csv"

# Output directory
OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "diagnostic_SPEI48_range"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Output files
OUTPUT_CSV = OUTPUT_DIR / "site_level_SPEI48_range_analysis.csv"
OUTPUT_PNG = OUTPUT_DIR / "Figure_R1_SPEI48_range_vs_sensitivity.png"
OUTPUT_PDF = OUTPUT_DIR / "Figure_R1_SPEI48_range_vs_sensitivity.pdf"

# Minimum observations per site
MIN_OBSERVATIONS = 5

print("="*80)
print("REVIEWER RESPONSE: SPEI-48 Range vs WUE_T Sensitivity Analysis")
print("="*80)
print("\nHypothesis: Does wider SPEI-48 range explain stronger sensitivity slopes?")
print(f"Minimum observations per site: {MIN_OBSERVATIONS}")
print("="*80)

# ============================================================================
# STEP 1: LOAD Q3 SLOPES DATA (SPEI-48)
# ============================================================================

print("\n[STEP 1] Loading Q3 slopes data (SPEI-48)...")

if not SLOPES_FILE.exists():
    print(f"  ERROR: Slopes file not found: {SLOPES_FILE}")
    exit(1)

df_slopes = pd.read_csv(SLOPES_FILE)
print(f"  Loaded {len(df_slopes)} sites from slopes file")
print(f"  Columns: {df_slopes.columns.tolist()}")

# Get slope column
slope_col = 'slope_theilsen'
print(f"  Using column '{slope_col}' for slope values")

# ============================================================================
# STEP 2: LOAD MONTHLY DATA AND COMPUTE SPEI-48 RANGE
# ============================================================================

print("\n[STEP 2] Loading monthly data and computing SPEI-48 range per site...")

if not MONTHLY_FILE.exists():
    print(f"  ERROR: Monthly file not found: {MONTHLY_FILE}")
    exit(1)

df_monthly = pd.read_csv(MONTHLY_FILE)
print(f"  Loaded monthly data: {len(df_monthly):,} rows")

# Get strict triple intersection sites
print("\n[STEP 3] Getting strict triple intersection sites...")

nn_medians = pd.read_csv(NN_MEDIANS_FILE)
wue_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())
eva_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())
tra_sites = set(nn_medians[nn_medians['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
print(f"  Strict triple intersection sites: {len(shared_sites)}")

# Filter monthly data
df_monthly = df_monthly[df_monthly['site_name'].isin(shared_sites)].copy()
print(f"  Monthly data after site filter: {len(df_monthly):,} rows")

# Apply strict month filter
strict_mask = (
    df_monthly['WUE'].notna() &
    df_monthly['WUE_eva'].notna() &
    df_monthly['WUE_tra'].notna()
)
df_monthly = df_monthly[strict_mask].copy()
print(f"  After strict month filter: {len(df_monthly):,} rows")

# Apply Trans_ratio > 0 filter
if 'Trans_ratio' in df_monthly.columns:
    initial_rows = len(df_monthly)
    df_monthly = df_monthly[df_monthly['Trans_ratio'] > 0].copy()
    print(f"  After Trans_ratio > 0 filter: {len(df_monthly):,} rows (removed {initial_rows - len(df_monthly)})")

# Keep only rows with valid SPEI_48 and WUE_tra
df_monthly = df_monthly.dropna(subset=['SPEI_48', 'WUE_tra'])
print(f"  After dropping missing SPEI_48/WUE_tra: {len(df_monthly):,} rows")

# ============================================================================
# STEP 4: COMPUTE SPEI-48 RANGE PER SITE
# ============================================================================

print("\n[STEP 4] Computing SPEI-48 range per site...")

spei_range_list = []
for site in df_slopes['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site]
    
    if len(site_data) >= MIN_OBSERVATIONS:
        spei_min = site_data['SPEI_48'].min()
        spei_max = site_data['SPEI_48'].max()
        spei_range = spei_max - spei_min
        n_obs = len(site_data)
        
        spei_range_list.append({
            'site_name': site,
            'spei48_min': spei_min,
            'spei48_max': spei_max,
            'spei48_range': spei_range,  # This is the column name
            'n_months_used': n_obs
        })
    else:
        print(f"  WARNING: {site} has insufficient data ({len(site_data)} obs)")

df_spei_range = pd.DataFrame(spei_range_list)
print(f"  Computed SPEI-48 range for {len(df_spei_range)} sites")
print(f"  SPEI range columns: {df_spei_range.columns.tolist()}")

# ============================================================================
# STEP 5: MERGE WITH SLOPES DATA
# ============================================================================

print("\n[STEP 5] Merging slope and SPEI-48 range data...")

# Merge on site_name
df_merged = df_slopes.merge(df_spei_range, on='site_name', how='inner')
print(f"  Merged data: {len(df_merged)} sites")
print(f"  Merged columns: {df_merged.columns.tolist()}")

# Calculate absolute slope
df_merged['abs_slope'] = df_merged[slope_col].abs()

# Identify correct SPEI range column (might be from slopes file or merged)
if 'spei48_range' in df_merged.columns:
    spei_range_col = 'spei48_range'
elif 'spei48_range_x' in df_merged.columns:
    spei_range_col = 'spei48_range_x'
elif 'spei48_range_y' in df_merged.columns:
    spei_range_col = 'spei48_range_y'
else:
    # Check if there's any column with 'range' in the name
    range_cols = [col for col in df_merged.columns if 'range' in col.lower()]
    if range_cols:
        spei_range_col = range_cols[0]
        print(f"  Found range column: {spei_range_col}")
    else:
        print(f"  ERROR: No SPEI range column found. Available columns: {df_merged.columns.tolist()}")
        exit(1)

print(f"\n  Using SPEI range column: {spei_range_col}")

print(f"\n  Summary statistics:")
print(f"    SPEI-48 range: min={df_merged[spei_range_col].min():.2f}, "
      f"max={df_merged[spei_range_col].max():.2f}, "
      f"mean={df_merged[spei_range_col].mean():.2f}")
print(f"    Absolute slope: min={df_merged['abs_slope'].min():.3f}, "
      f"max={df_merged['abs_slope'].max():.3f}, "
      f"mean={df_merged['abs_slope'].mean():.3f}")

# ============================================================================
# STEP 6: STATISTICAL TESTS
# ============================================================================

print("\n" + "="*80)
print("STATISTICAL RESULTS (SPEI-48)")
print("="*80)

# Spearman rank correlation
spearman_rho, spearman_p = spearmanr(df_merged[spei_range_col], df_merged['abs_slope'])

print(f"\n📊 SPEARMAN RANK CORRELATION (Primary Test):")
print(f"   rho = {spearman_rho:.4f}")
print(f"   p-value = {spearman_p:.4f}")

if spearman_p < 0.001:
    p_text = "p < 0.001"
elif spearman_p < 0.01:
    p_text = f"p = {spearman_p:.3f}"
elif spearman_p < 0.05:
    p_text = f"p = {spearman_p:.3f}"
else:
    p_text = f"p = {spearman_p:.3f} (not significant)"

print(f"\n📈 INTERPRETATION:")

if spearman_rho > 0 and spearman_p < 0.05:
    print(f"   → POSITIVE correlation: Sites with wider SPEI-48 range tend to have")
    print(f"     stronger (larger absolute) WUE_T sensitivity slopes.")
elif spearman_rho < 0 and spearman_p < 0.05:
    print(f"   → NEGATIVE correlation: Sites with wider SPEI-48 range tend to have")
    print(f"     weaker (smaller absolute) WUE_T sensitivity slopes.")
else:
    print(f"   → NO significant correlation (p > 0.05):")
    print(f"     SPEI-48 range does NOT explain the magnitude of WUE_T sensitivity.")

# Linear regression
slope_reg, intercept, r_value, p_reg, std_err = linregress(
    df_merged[spei_range_col], df_merged['abs_slope']
)
r_squared = r_value**2

print(f"\n📊 LINEAR REGRESSION (Visualization Only):")
print(f"   slope = {slope_reg:.4f}")
print(f"   R² = {r_squared:.4f}")
print(f"   p-value = {p_reg:.4f}")

# ============================================================================
# STEP 7: CREATE SCATTER PLOT
# ============================================================================

print("\n[STEP 7] Creating scatter plot...")

# Coastal region colors (consistent with Panel A and B)
coast_colors = {
    'AK_coast': '#9E0142',
    'West_Coast': '#F46D43',
    'Gulf_of_America': '#74ADD1',
    'Atlantic_Coast': '#2B83BA',
    'Atlantic_Coast_North': '#2B83BA',
    'Southeast_Atlantic': '#ABD9E9'
}

# Check if coast_region_analysis column exists
if 'coast_region_analysis' in df_merged.columns:
    use_coast_colors = True
    print("  Coloring points by coastal region")
else:
    use_coast_colors = False
    print("  Warning: coast_region_analysis not found, using single color")

# Create figure
fig, ax = plt.subplots(1, 1, figsize=(12, 10), dpi=300)

# Plot points
if use_coast_colors:
    for coast in df_merged['coast_region_analysis'].unique():
        coast_data = df_merged[df_merged['coast_region_analysis'] == coast]
        color = coast_colors.get(coast, '#9E9E9E')
        ax.scatter(coast_data[spei_range_col], coast_data['abs_slope'],
                  color=color, s=150, alpha=0.8, edgecolors='black', linewidth=1.5,
                  label=coast)
else:
    ax.scatter(df_merged[spei_range_col], df_merged['abs_slope'],
              color='#2C7FB8', s=150, alpha=0.7, edgecolors='black', linewidth=1.5)

# Add regression line
x_range = np.linspace(df_merged[spei_range_col].min(), df_merged[spei_range_col].max(), 100)
y_pred = slope_reg * x_range + intercept
ax.plot(x_range, y_pred, color='red', linewidth=2.5, linestyle='--', 
        label=f'Linear fit (R²={r_squared:.2f})')

# Labels
ax.set_xlabel('SPEI-48 Range (max - min)', fontsize=28, fontweight='bold')
ax.set_ylabel('|WUE$_T$ Sensitivity Slope|', fontsize=28, fontweight='bold')
ax.set_title('Diagnostic: SPEI-48 Range vs WUE$_T$ Sensitivity Magnitude', 
             fontsize=24, fontweight='bold', pad=15)

# Annotation box
if spearman_p < 0.001:
    p_annotation = "p < 0.001"
else:
    p_annotation = f"p = {spearman_p:.3f}"

annotation_text = (f"Spearman ρ = {spearman_rho:.3f}\n"
                  f"{p_annotation}\n"
                  f"N = {len(df_merged)} sites")

ax.text(0.95, 0.95, annotation_text, transform=ax.transAxes,
        fontsize=22, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# Styling
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_axisbelow(True)
ax.tick_params(axis='both', labelsize=24, width=2, length=8)

for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')

# Legend
if use_coast_colors:
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(handles, labels, loc='lower right', fontsize=16, 
                 framealpha=0.9, edgecolor='gray', title='Coastal Region',
                 title_fontsize=18)

# ============================================================================
# SAVE FIGURE
# ============================================================================

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_PDF, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n  ✓ Saved figure: {OUTPUT_PNG}")
print(f"  ✓ Saved figure: {OUTPUT_PDF}")

plt.show()

# ============================================================================
# STEP 8: SAVE CSV
# ============================================================================

print("\n[STEP 8] Saving site-level data...")

output_cols = ['site_name', slope_col, 'abs_slope', 
               spei_range_col, 'n_months_used']

if 'coast_region_analysis' in df_merged.columns:
    output_cols.append('coast_region_analysis')

df_output = df_merged[output_cols].copy()
df_output = df_output.rename(columns={slope_col: 'slope_theilsen', spei_range_col: 'spei48_range'})

df_output.to_csv(OUTPUT_CSV, index=False)
print(f"  ✓ Saved CSV: {OUTPUT_CSV}")
print(f"    Records: {len(df_output)}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("MANUSCRIPT-READY SUMMARY (SPEI-48)")
print("="*80)

print(f"""
Diagnostic Analysis: SPEI-48 Range vs |WUE_T Sensitivity|

DATA SUMMARY:
  • Sites analyzed: {len(df_merged)}
  • SPEI-48 range: {df_merged[spei_range_col].min():.2f} – {df_merged[spei_range_col].max():.2f} (mean={df_merged[spei_range_col].mean():.2f})
  • |Slope| range: {df_merged['abs_slope'].min():.3f} – {df_merged['abs_slope'].max():.3f} (mean={df_merged['abs_slope'].mean():.3f})

SPEARMAN RANK CORRELATION:
  • ρ = {spearman_rho:.4f}
  • p = {spearman_p:.4f}
  • Interpretation: {'SIGNIFICANT' if spearman_p < 0.05 else 'NOT SIGNIFICANT'}

CONCLUSION:
  {('A significant positive correlation exists between SPEI-48 range and |WUE_T sensitivity|.')
    if spearman_rho > 0 and spearman_p < 0.05 else
   ('A significant negative correlation exists between SPEI-48 range and |WUE_T sensitivity|.')
    if spearman_rho < 0 and spearman_p < 0.05 else
   ('No significant correlation was detected between SPEI-48 range and |WUE_T sensitivity|.')}
""")

print("="*80)
print("DIAGNOSTIC ANALYSIS COMPLETE (SPEI-48)")
print("="*80)