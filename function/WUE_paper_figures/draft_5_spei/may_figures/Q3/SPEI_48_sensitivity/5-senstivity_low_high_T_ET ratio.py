# -*- coding: utf-8 -*-
"""
DIAGNOSTIC ANALYSIS: T:ET Ratio vs WUE_T Sensitivity (SPEI-48)
====================================================
Test whether sites with higher T:ET ratios (more transpiration-dominated)
show stronger or weaker WUE_T sensitivity to SPEI-48.

Research Questions:
1. Does median T:ET ratio correlate with WUE_T sensitivity slope?
2. Does median T:ET ratio correlate with absolute WUE_T sensitivity?

Data Sources:
- Q3_WUET_SPEI48_slopes.csv (has SPEI-48 slopes)
- TET_gradient_site_NN_medians.csv (has median_Trans_ratio for T:ET)
====================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, linregress
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results")

# Input files
SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI48_slopes.csv"
TET_FILE = BASE_DIR / "T_ET_ratio" / "TET_gradient_site_NN_medians.csv"

# Output directory
OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "diagnostic_TET_analysis_SPEI48"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("DIAGNOSTIC ANALYSIS: T:ET Ratio vs WUE_T Sensitivity (SPEI-48)")
print("="*80)

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[STEP 1] Loading data...")

# Load SPEI-48 slopes file
df_slopes = pd.read_csv(SLOPES_FILE)
print(f"  Loaded SPEI-48 slopes file: {len(df_slopes)} sites")
print(f"  Columns: {df_slopes.columns.tolist()}")

# Load TET file (has median_Trans_ratio)
df_tet = pd.read_csv(TET_FILE)
print(f"\n  Loaded TET file: {len(df_tet)} sites")
print(f"  Columns: {df_tet.columns.tolist()}")

# Check for T:ET ratio column
if 'median_Trans_ratio' in df_tet.columns:
    tet_col = 'median_Trans_ratio'
    print(f"  Using T:ET column: {tet_col}")
else:
    print(f"  ERROR: No T:ET ratio column found. Available: {df_tet.columns.tolist()}")
    exit(1)

# Merge on site_name
df_merged = df_slopes.merge(df_tet[['site_name', tet_col, 'Salinity_Category', 'IGBP', 'climate']], 
                            on='site_name', how='inner')
print(f"\n  Merged data: {len(df_merged)} sites")

# Get slope column
slope_col = 'slope_theilsen'
print(f"  Using slope column: {slope_col}")

# ============================================================================
# STEP 2: CHECK DATA QUALITY
# ============================================================================

print("\n[STEP 2] Data quality check...")

# Rename T:ET column for clarity
df_merged = df_merged.rename(columns={tet_col: 'TET_ratio'})

# Check for missing values
print(f"\n  Missing values:")
print(f"    {slope_col}: {df_merged[slope_col].isna().sum()}")
print(f"    TET_ratio: {df_merged['TET_ratio'].isna().sum()}")

# Drop any rows with missing values
df_clean = df_merged.dropna(subset=[slope_col, 'TET_ratio']).copy()
print(f"\n  After dropping missing: {len(df_clean)} sites")

# Summary statistics
print(f"\n  Summary statistics (SPEI-48):")
print(f"    T:ET ratio (median_Trans_ratio):")
print(f"      Min: {df_clean['TET_ratio'].min():.3f}")
print(f"      Max: {df_clean['TET_ratio'].max():.3f}")
print(f"      Mean: {df_clean['TET_ratio'].mean():.3f}")
print(f"      Median: {df_clean['TET_ratio'].median():.3f}")

print(f"\n    Slope (WUE_T sensitivity to SPEI-48):")
print(f"      Min: {df_clean[slope_col].min():.3f}")
print(f"      Max: {df_clean[slope_col].max():.3f}")
print(f"      Mean: {df_clean[slope_col].mean():.3f}")
print(f"      Median: {df_clean[slope_col].median():.3f}")

# Calculate absolute slope
df_clean['abs_slope'] = df_clean[slope_col].abs()

print(f"\n    Absolute slope:")
print(f"      Min: {df_clean['abs_slope'].min():.3f}")
print(f"      Max: {df_clean['abs_slope'].max():.3f}")
print(f"      Mean: {df_clean['abs_slope'].mean():.3f}")
print(f"      Median: {df_clean['abs_slope'].median():.3f}")

# ============================================================================
# STEP 3: STATISTICAL TESTS - T:ET vs SLOPE (SPEI-48)
# ============================================================================

print("\n" + "="*80)
print("STATISTICAL RESULTS (SPEI-48)")
print("="*80)

# Test 1: T:ET vs Slope (direction of sensitivity)
rho1, p1 = spearmanr(df_clean['TET_ratio'], df_clean[slope_col])

print(f"\n📊 TEST 1: T:ET Ratio vs WUE_T Sensitivity Slope (SPEI-48)")
print(f"   Research question: Do transpiration-dominated sites show more")
print(f"   positive or negative sensitivity to wetter conditions?")
print(f"   ===========================================================")
print(f"   Spearman ρ = {rho1:.4f}")
print(f"   p-value = {p1:.4f}")

if p1 < 0.001:
    p_text1 = "p < 0.001"
elif p1 < 0.01:
    p_text1 = f"p = {p1:.3f}"
elif p1 < 0.05:
    p_text1 = f"p = {p1:.3f}"
else:
    p_text1 = f"p = {p1:.3f} (not significant)"

if rho1 > 0 and p1 < 0.05:
    interpretation1 = ("POSITIVE correlation: Sites with higher T:ET ratios "
                      "(more transpiration-dominated) tend to have MORE POSITIVE "
                      "sensitivity slopes, meaning they increase WUE_T more "
                      "strongly under wetter conditions (SPEI-48).")
elif rho1 < 0 and p1 < 0.05:
    interpretation1 = ("NEGATIVE correlation: Sites with higher T:ET ratios "
                      "(more transpiration-dominated) tend to have MORE NEGATIVE "
                      "sensitivity slopes, meaning they decrease WUE_T under "
                      "wetter conditions (or increase under drier conditions).")
else:
    interpretation1 = ("NO significant correlation: T:ET ratio does not predict "
                      "the direction of WUE_T sensitivity to SPEI-48.")

print(f"\n   Interpretation: {interpretation1}")

# ============================================================================
# TEST 2: T:ET vs Absolute Slope (magnitude of sensitivity)
# ============================================================================

rho2, p2 = spearmanr(df_clean['TET_ratio'], df_clean['abs_slope'])

print(f"\n📊 TEST 2: T:ET Ratio vs |WUE_T Sensitivity| (Magnitude) - SPEI-48")
print(f"   Research question: Are transpiration-dominated sites more or less")
print(f"   sensitive overall (regardless of direction) to SPEI-48?")
print(f"   ===========================================================")
print(f"   Spearman ρ = {rho2:.4f}")
print(f"   p-value = {p2:.4f}")

if p2 < 0.001:
    p_text2 = "p < 0.001"
elif p2 < 0.01:
    p_text2 = f"p = {p2:.3f}"
elif p2 < 0.05:
    p_text2 = f"p = {p2:.3f}"
else:
    p_text2 = f"p = {p2:.3f} (not significant)"

if rho2 > 0 and p2 < 0.05:
    interpretation2 = ("POSITIVE correlation: Transpiration-dominated sites "
                      "(higher T:ET) show STRONGER overall sensitivity "
                      "(larger absolute slopes) to SPEI-48.")
elif rho2 < 0 and p2 < 0.05:
    interpretation2 = ("NEGATIVE correlation: Transpiration-dominated sites "
                      "(higher T:ET) show WEAKER overall sensitivity "
                      "(smaller absolute slopes) to SPEI-48.")
else:
    interpretation2 = ("NO significant correlation: T:ET ratio does not predict "
                      "the magnitude of WUE_T sensitivity to SPEI-48.")

print(f"\n   Interpretation: {interpretation2}")

# ============================================================================
# STEP 4: CREATE SCATTER PLOTS
# ============================================================================

print("\n[STEP 4] Creating scatter plots (SPEI-48)...")

# Colors for salinity categories
salinity_colors = {
    'Freshwater': '#0000FF',
    'Saline': '#FFA500',
    'Upland': '#800080',
    'Brackish': '#008080'
}

# FIGURE 1: T:ET vs Slope (SPEI-48)
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 8), dpi=300)

# Color by salinity if available
if 'Salinity_Category' in df_clean.columns:
    for salinity in df_clean['Salinity_Category'].unique():
        subset = df_clean[df_clean['Salinity_Category'] == salinity]
        color = salinity_colors.get(salinity, '#9E9E9E')
        ax1.scatter(subset['TET_ratio'], subset[slope_col],
                   color=color, s=120, alpha=0.7, edgecolors='black', linewidth=1,
                   label=salinity)
else:
    ax1.scatter(df_clean['TET_ratio'], df_clean[slope_col],
               color='#2C7FB8', s=120, alpha=0.7, edgecolors='black', linewidth=1)

# Add regression line
x_range = np.linspace(df_clean['TET_ratio'].min(), df_clean['TET_ratio'].max(), 100)
slope_reg1, intercept1, r_val1, p_reg1, _ = linregress(df_clean['TET_ratio'], df_clean[slope_col])
y_pred1 = slope_reg1 * x_range + intercept1
ax1.plot(x_range, y_pred1, color='red', linewidth=2, linestyle='--',
        label=f'Linear fit (R²={r_val1**2:.2f})')

# Labels and formatting
ax1.set_xlabel('Median T:ET Ratio (Transpiration/Evapotranspiration)', fontsize=14, fontweight='bold')
ax1.set_ylabel('WUE$_T$ Sensitivity to SPEI-48 (ΔWUE$_T$/ΔSPEI-48)', fontsize=14, fontweight='bold')
ax1.set_title('T:ET Ratio vs WUE$_T$ Sensitivity to SPEI-48', fontsize=16, fontweight='bold')

# Add annotation
annotation1 = (f"Spearman ρ = {rho1:.3f}\n"
              f"{p_text1}\n"
              f"N = {len(df_clean)} sites")
ax1.text(0.95, 0.05, annotation1, transform=ax1.transAxes,
        fontsize=12, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_axisbelow(True)

if 'Salinity_Category' in df_clean.columns:
    ax1.legend(loc='upper left', fontsize=10, framealpha=0.9)

# FIGURE 2: T:ET vs Absolute Slope (SPEI-48)
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 8), dpi=300)

if 'Salinity_Category' in df_clean.columns:
    for salinity in df_clean['Salinity_Category'].unique():
        subset = df_clean[df_clean['Salinity_Category'] == salinity]
        color = salinity_colors.get(salinity, '#9E9E9E')
        ax2.scatter(subset['TET_ratio'], subset['abs_slope'],
                   color=color, s=120, alpha=0.7, edgecolors='black', linewidth=1,
                   label=salinity)
else:
    ax2.scatter(df_clean['TET_ratio'], df_clean['abs_slope'],
               color='#2C7FB8', s=120, alpha=0.7, edgecolors='black', linewidth=1)

# Add regression line
slope_reg2, intercept2, r_val2, p_reg2, _ = linregress(df_clean['TET_ratio'], df_clean['abs_slope'])
y_pred2 = slope_reg2 * x_range + intercept2
ax2.plot(x_range, y_pred2, color='red', linewidth=2, linestyle='--',
        label=f'Linear fit (R²={r_val2**2:.2f})')

ax2.set_xlabel('Median T:ET Ratio (Transpiration/Evapotranspiration)', fontsize=14, fontweight='bold')
ax2.set_ylabel('|WUE$_T$ Sensitivity to SPEI-48|', fontsize=14, fontweight='bold')
ax2.set_title('T:ET Ratio vs |WUE$_T$ Sensitivity| to SPEI-48', fontsize=16, fontweight='bold')

# Add annotation
annotation2 = (f"Spearman ρ = {rho2:.3f}\n"
              f"{p_text2}\n"
              f"N = {len(df_clean)} sites")
ax2.text(0.95, 0.95, annotation2, transform=ax2.transAxes,
        fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_axisbelow(True)

if 'Salinity_Category' in df_clean.columns:
    ax2.legend(loc='upper left', fontsize=10, framealpha=0.9)

# ============================================================================
# SAVE FIGURES AND DATA
# ============================================================================

print("\n[STEP 5] Saving outputs...")

# Save figures
fig1.tight_layout()
fig1.savefig(OUTPUT_DIR / "Figure_TET_vs_slope_SPEI48.png", dpi=300, bbox_inches='tight')
fig1.savefig(OUTPUT_DIR / "Figure_TET_vs_slope_SPEI48.pdf", bbox_inches='tight')

fig2.tight_layout()
fig2.savefig(OUTPUT_DIR / "Figure_TET_vs_abs_slope_SPEI48.png", dpi=300, bbox_inches='tight')
fig2.savefig(OUTPUT_DIR / "Figure_TET_vs_abs_slope_SPEI48.pdf", bbox_inches='tight')

print(f"  ✓ Saved: Figure_TET_vs_slope_SPEI48.png/pdf")
print(f"  ✓ Saved: Figure_TET_vs_abs_slope_SPEI48.png/pdf")

# Save CSV with all data
output_cols = ['site_name', 'TET_ratio', slope_col, 'abs_slope']
if 'Salinity_Category' in df_clean.columns:
    output_cols.append('Salinity_Category')
if 'p_value' in df_clean.columns:
    output_cols.append('p_value')
if 'is_significant' in df_clean.columns:
    output_cols.append('is_significant')
if 'coast_region_analysis' in df_clean.columns:
    output_cols.append('coast_region_analysis')

df_output = df_clean[output_cols].copy()
df_output = df_output.rename(columns={slope_col: 'slope_theilsen_SPEI48'})
df_output.to_csv(OUTPUT_DIR / "TET_sensitivity_analysis_data_SPEI48.csv", index=False)
print(f"  ✓ Saved CSV: TET_sensitivity_analysis_data_SPEI48.csv")

plt.show()

# ============================================================================
# FINAL SUMMARY FOR MANUSCRIPT
# ============================================================================

print("\n" + "="*80)
print("MANUSCRIPT-READY SUMMARY (SPEI-48)")
print("="*80)

print(f"""
T:ET Ratio vs WUE_T Sensitivity to SPEI-48 Analysis

DATA SUMMARY:
  • Sites analyzed: {len(df_clean)}
  • T:ET ratio range: {df_clean['TET_ratio'].min():.3f} – {df_clean['TET_ratio'].max():.3f}
  • SPEI-48 Slope range: {df_clean[slope_col].min():.3f} – {df_clean[slope_col].max():.3f}
  • |SPEI-48 Slope| range: {df_clean['abs_slope'].min():.3f} – {df_clean['abs_slope'].max():.3f}

TEST 1: T:ET vs Sensitivity Direction (Slope to SPEI-48)
  • Spearman ρ = {rho1:.4f}
  • p = {p1:.4f}
  • {interpretation1}

TEST 2: T:ET vs Sensitivity Magnitude (|Slope| to SPEI-48)
  • Spearman ρ = {rho2:.4f}
  • p = {p2:.4f}
  • {interpretation2}

CONCLUSION:
  {'The T:ET ratio shows a significant relationship with WUE_T sensitivity to SPEI-48, '
   'suggesting that evaporation vs transpiration dominance influences '
   'how coastal wetlands respond to multi-year hydroclimatic variability.' if (p1 < 0.05 or p2 < 0.05) else
   'No significant relationship was detected between T:ET ratio and '
   'WUE_T sensitivity to SPEI-48, suggesting that evaporation vs transpiration '
   'dominance does not strongly influence sensitivity patterns at 48-month timescales.'}
""")

print("\n" + "="*80)
print("DIAGNOSTIC ANALYSIS COMPLETE (SPEI-48)")
print("="*80)
print(f"\n📁 Output files saved to: {OUTPUT_DIR}")
print("="*80)