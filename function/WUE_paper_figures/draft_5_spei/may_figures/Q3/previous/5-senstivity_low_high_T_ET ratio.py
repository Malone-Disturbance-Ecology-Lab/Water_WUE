# -*- coding: utf-8 -*-
"""
Created on Sun May 10 11:17:03 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
DIAGNOSTIC ANALYSIS: T:ET Ratio vs WUE_T Sensitivity
====================================================
Test whether sites with higher T:ET ratios (more transpiration-dominated)
show stronger or weaker WUE_T sensitivity.

Research Questions:
1. Does median T:ET ratio correlate with WUE_T sensitivity slope?
   (positive slope means wetter conditions increase WUE_T)

2. Does median T:ET ratio correlate with absolute WUE_T sensitivity?
   (transpiration-dominated sites more or less sensitive overall)

Data Sources:
- Q3_WUET_SPEI6_slopes.csv (has slope and median_TET)
- TET_gradient_site_NN_medians.csv (has median_Trans_ratio as alternative)
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
SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI6_slopes.csv"
TET_FILE = BASE_DIR / "T_ET_ratio" / "TET_gradient_site_NN_medians.csv"
ALTERNATIVE_FILE = BASE_DIR / "senstivity_april" / "Q3_final_dataset.csv"

# Output directory
OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "diagnostic_TET_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("DIAGNOSTIC ANALYSIS: T:ET Ratio vs WUE_T Sensitivity")
print("="*80)

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[STEP 1] Loading data...")

# Load slopes file (has median_TET already)
df_slopes = pd.read_csv(SLOPES_FILE)
print(f"  Loaded slopes file: {len(df_slopes)} sites")
print(f"  Columns available: {df_slopes.columns.tolist()}")

# Load TET gradient file for alternative T:ET measure
df_tet = pd.read_csv(TET_FILE)
print(f"\n  Loaded TET file: {len(df_tet)} sites")
print(f"  Columns: {df_tet.columns.tolist()}")

# Merge to ensure consistency (should already have 57 sites)
df_merged = df_slopes.merge(df_tet[['site_name', 'median_Trans_ratio']], on='site_name', how='inner')
print(f"\n  Merged data: {len(df_merged)} sites")

# ============================================================================
# STEP 2: CHECK DATA QUALITY
# ============================================================================

print("\n[STEP 2] Data quality check...")

# Check for missing values
print(f"\n  Missing values:")
print(f"    slope: {df_merged['slope'].isna().sum()}")
print(f"    median_TET (from slopes): {df_merged['median_TET'].isna().sum()}")
print(f"    median_Trans_ratio (from TET file): {df_merged['median_Trans_ratio'].isna().sum()}")

# Summary statistics
print(f"\n  Summary statistics:")
print(f"    T:ET ratio (median_TET):")
print(f"      Min: {df_merged['median_TET'].min():.3f}")
print(f"      Max: {df_merged['median_TET'].max():.3f}")
print(f"      Mean: {df_merged['median_TET'].mean():.3f}")
print(f"      Median: {df_merged['median_TET'].median():.3f}")

print(f"\n    Slope (WUE_T sensitivity):")
print(f"      Min: {df_merged['slope'].min():.3f}")
print(f"      Max: {df_merged['slope'].max():.3f}")
print(f"      Mean: {df_merged['slope'].mean():.3f}")
print(f"      Median: {df_merged['slope'].median():.3f}")

# Calculate absolute slope
df_merged['abs_slope'] = df_merged['slope'].abs()

print(f"\n    Absolute slope:")
print(f"      Min: {df_merged['abs_slope'].min():.3f}")
print(f"      Max: {df_merged['abs_slope'].max():.3f}")
print(f"      Mean: {df_merged['abs_slope'].mean():.3f}")
print(f"      Median: {df_merged['abs_slope'].median():.3f}")

# ============================================================================
# STEP 3: STATISTICAL TESTS - T:ET vs SLOPE
# ============================================================================

print("\n" + "="*80)
print("STATISTICAL RESULTS")
print("="*80)

# Test 1: T:ET vs Slope (direction of sensitivity)
rho1, p1 = spearmanr(df_merged['median_TET'], df_merged['slope'])

print(f"\n📊 TEST 1: T:ET Ratio vs WUE_T Sensitivity Slope")
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
                      "strongly under wetter conditions.")
elif rho1 < 0 and p1 < 0.05:
    interpretation1 = ("NEGATIVE correlation: Sites with higher T:ET ratios "
                      "(more transpiration-dominated) tend to have MORE NEGATIVE "
                      "sensitivity slopes, meaning they decrease WUE_T under "
                      "wetter conditions (or increase under drier conditions).")
else:
    interpretation1 = ("NO significant correlation: T:ET ratio does not predict "
                      "the direction of WUE_T sensitivity to wetter conditions.")

print(f"\n   Interpretation: {interpretation1}")

# ============================================================================
# TEST 2: T:ET vs Absolute Slope (magnitude of sensitivity)
# ============================================================================

rho2, p2 = spearmanr(df_merged['median_TET'], df_merged['abs_slope'])

print(f"\n📊 TEST 2: T:ET Ratio vs |WUE_T Sensitivity| (Magnitude)")
print(f"   Research question: Are transpiration-dominated sites more or less")
print(f"   sensitive overall (regardless of direction)?")
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
                      "(larger absolute slopes).")
elif rho2 < 0 and p2 < 0.05:
    interpretation2 = ("NEGATIVE correlation: Transpiration-dominated sites "
                      "(higher T:ET) show WEAKER overall sensitivity "
                      "(smaller absolute slopes).")
else:
    interpretation2 = ("NO significant correlation: T:ET ratio does not predict "
                      "the magnitude of WUE_T sensitivity.")

print(f"\n   Interpretation: {interpretation2}")

# ============================================================================
# STEP 4: CREATE SCATTER PLOTS
# ============================================================================

print("\n[STEP 3] Creating scatter plots...")

# Define colors for salinity categories (if available)
salinity_colors = {
    'Freshwater': '#0000FF',
    'Saline': '#FFA500',
    'Upland': '#800080',
    'Brackish': '#008080'
}

# FIGURE 1: T:ET vs Slope
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 8), dpi=300)

# Color by salinity if available
if 'Salinity_Category' in df_merged.columns:
    for salinity in df_merged['Salinity_Category'].unique():
        subset = df_merged[df_merged['Salinity_Category'] == salinity]
        color = salinity_colors.get(salinity, '#9E9E9E')
        ax1.scatter(subset['median_TET'], subset['slope'],
                   color=color, s=120, alpha=0.7, edgecolors='black', linewidth=1,
                   label=salinity)
else:
    ax1.scatter(df_merged['median_TET'], df_merged['slope'],
               color='#2C7FB8', s=120, alpha=0.7, edgecolors='black', linewidth=1)

# Add regression line
x_range = np.linspace(df_merged['median_TET'].min(), df_merged['median_TET'].max(), 100)
slope_reg1, intercept1, r_val1, p_reg1, _ = linregress(df_merged['median_TET'], df_merged['slope'])
y_pred1 = slope_reg1 * x_range + intercept1
ax1.plot(x_range, y_pred1, color='red', linewidth=2, linestyle='--',
        label=f'Linear fit (R²={r_val1**2:.2f})')

# Labels and formatting
ax1.set_xlabel('Median T:ET Ratio', fontsize=14, fontweight='bold')
ax1.set_ylabel('WUE$_T$ Sensitivity Slope', fontsize=14, fontweight='bold')
ax1.set_title('T:ET Ratio vs WUE$_T$ Sensitivity', fontsize=16, fontweight='bold')

# Add annotation
annotation1 = (f"Spearman ρ = {rho1:.3f}\n"
              f"{p_text1}\n"
              f"N = {len(df_merged)} sites")
ax1.text(0.95, 0.05, annotation1, transform=ax1.transAxes,
        fontsize=12, verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax1.grid(True, alpha=0.3, linestyle='--')
ax1.set_axisbelow(True)

if 'Salinity_Category' in df_merged.columns:
    ax1.legend(loc='upper left', fontsize=10, framealpha=0.9)

# FIGURE 2: T:ET vs Absolute Slope
fig2, ax2 = plt.subplots(1, 1, figsize=(10, 8), dpi=300)

if 'Salinity_Category' in df_merged.columns:
    for salinity in df_merged['Salinity_Category'].unique():
        subset = df_merged[df_merged['Salinity_Category'] == salinity]
        color = salinity_colors.get(salinity, '#9E9E9E')
        ax2.scatter(subset['median_TET'], subset['abs_slope'],
                   color=color, s=120, alpha=0.7, edgecolors='black', linewidth=1,
                   label=salinity)
else:
    ax2.scatter(df_merged['median_TET'], df_merged['abs_slope'],
               color='#2C7FB8', s=120, alpha=0.7, edgecolors='black', linewidth=1)

# Add regression line
slope_reg2, intercept2, r_val2, p_reg2, _ = linregress(df_merged['median_TET'], df_merged['abs_slope'])
y_pred2 = slope_reg2 * x_range + intercept2
ax2.plot(x_range, y_pred2, color='red', linewidth=2, linestyle='--',
        label=f'Linear fit (R²={r_val2**2:.2f})')

ax2.set_xlabel('Median T:ET Ratio', fontsize=14, fontweight='bold')
ax2.set_ylabel('|WUE$_T$ Sensitivity Slope|', fontsize=14, fontweight='bold')
ax2.set_title('T:ET Ratio vs |WUE$_T$ Sensitivity|', fontsize=16, fontweight='bold')

# Add annotation
annotation2 = (f"Spearman ρ = {rho2:.3f}\n"
              f"{p_text2}\n"
              f"N = {len(df_merged)} sites")
ax2.text(0.95, 0.95, annotation2, transform=ax2.transAxes,
        fontsize=12, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.grid(True, alpha=0.3, linestyle='--')
ax2.set_axisbelow(True)

if 'Salinity_Category' in df_merged.columns:
    ax2.legend(loc='upper left', fontsize=10, framealpha=0.9)

# ============================================================================
# SAVE FIGURES AND DATA
# ============================================================================

print("\n[STEP 4] Saving outputs...")

# Save figures
fig1.tight_layout()
fig1.savefig(OUTPUT_DIR / "Figure_TET_vs_slope.png", dpi=300, bbox_inches='tight')
fig1.savefig(OUTPUT_DIR / "Figure_TET_vs_slope.pdf", bbox_inches='tight')

fig2.tight_layout()
fig2.savefig(OUTPUT_DIR / "Figure_TET_vs_abs_slope.png", dpi=300, bbox_inches='tight')
fig2.savefig(OUTPUT_DIR / "Figure_TET_vs_abs_slope.pdf", bbox_inches='tight')

print(f"  ✓ Saved: Figure_TET_vs_slope.png/pdf")
print(f"  ✓ Saved: Figure_TET_vs_abs_slope.png/pdf")

# Save CSV with all data
output_cols = ['site_name', 'median_TET', 'slope', 'abs_slope']
if 'Salinity_Category' in df_merged.columns:
    output_cols.append('Salinity_Category')
if 'median_Trans_ratio' in df_merged.columns:
    output_cols.append('median_Trans_ratio')

df_output = df_merged[output_cols].copy()
df_output.to_csv(OUTPUT_DIR / "TET_sensitivity_analysis_data.csv", index=False)
print(f"  ✓ Saved CSV: TET_sensitivity_analysis_data.csv")

plt.show()

# ============================================================================
# FINAL SUMMARY FOR MANUSCRIPT
# ============================================================================

print("\n" + "="*80)
print("MANUSCRIPT-READY SUMMARY")
print("="*80)

print(f"""
T:ET Ratio vs WUE_T Sensitivity Analysis

DATA SUMMARY:
  • Sites analyzed: {len(df_merged)}
  • T:ET ratio range: {df_merged['median_TET'].min():.3f} – {df_merged['median_TET'].max():.3f}
  • Slope range: {df_merged['slope'].min():.3f} – {df_merged['slope'].max():.3f}
  • |Slope| range: {df_merged['abs_slope'].min():.3f} – {df_merged['abs_slope'].max():.3f}

TEST 1: T:ET vs Sensitivity Direction (Slope)
  • Spearman ρ = {rho1:.4f}
  • p = {p1:.4f}
  • {interpretation1}

TEST 2: T:ET vs Sensitivity Magnitude (|Slope|)
  • Spearman ρ = {rho2:.4f}
  • p = {p2:.4f}
  • {interpretation2}

CONCLUSION:
  {'The T:ET ratio shows a significant relationship with WUE_T sensitivity, '
   'suggesting that evaporation vs transpiration dominance influences '
   'how wetlands respond to hydroclimatic variability.' if (p1 < 0.05 or p2 < 0.05) else
   'No significant relationship was detected between T:ET ratio and '
   'WUE_T sensitivity, suggesting that evaporation vs transpiration '
   'dominance does not strongly influence sensitivity patterns.'}
""")

print("\n" + "="*80)
print("DIAGNOSTIC ANALYSIS COMPLETE")
print("="*80)
print(f"\n📁 Output files saved to: {OUTPUT_DIR}")
print("="*80)