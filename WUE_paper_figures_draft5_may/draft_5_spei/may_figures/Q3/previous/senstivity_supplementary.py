# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 12:08:11 2026

@author: ammar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import kruskal, wilcoxon
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# SETUP OUTPUT DIRECTORY
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("SUPPLEMENTARY FIGURE: Dry-Wet & Ecosystem Sensitivity")
print("="*70)

# ============================================
# COLOR MAPPING
# ============================================

# Ecosystem colors (matches Q1)
ecosystem_colors = {
    'freshwater': '#0000FF',   # Blue
    'saline': '#FFA500',       # Orange
    'upland': '#800080'        # Purple
}

ecosystem_full_names = {
    'freshwater': 'Freshwater',
    'saline': 'Saline',
    'upland': 'Upland'
}

# Coast colors (same as Figure 2)
coast_colors = {
    'AK_coast': '#1E88E5',      # Blue
    'West_Coast': '#F57C00',    # Orange
    'Gulf_of_America': '#7B1FA2', # Purple
    'Atlantic_Coast': '#E91E63'   # Deep pink
}

coast_full_names = {
    'AK_coast': 'AK Coast',
    'West_Coast': 'West Coast',
    'Gulf_of_America': 'Gulf of America',
    'Atlantic_Coast': 'Atlantic Coast'
}

# ============================================
# LOAD DATA
# ============================================

print("\n📂 Loading data...")

# Load dry-wet data
dry_wet_file = output_dir / 'Q3_site_sensitivity_dry_wet.csv'
if dry_wet_file.exists():
    df_dw = pd.read_csv(dry_wet_file)
    # Calculate diff column
    df_dw['diff'] = df_dw['slope_wet_theilsen'] - df_dw['slope_dry_theilsen']
    print(f"  Loaded {len(df_dw)} dry-wet sites")
    print(f"  Calculated diff = wet - dry")
    
    # Calculate dry-wet statistics BEFORE any manipulation
    n_wet_dom = (df_dw['diff'] > 0).sum()
    n_dry_dom = (df_dw['diff'] < 0).sum()
    median_diff = df_dw['diff'].median()
    
    print(f"  Wet-dominated (>0): {n_wet_dom} sites")
    print(f"  Dry-dominated (<0): {n_dry_dom} sites")
    print(f"  Median difference: {median_diff:.4f}")
else:
    df_dw = pd.DataFrame()
    print(f"  WARNING: Dry-wet file not found")
    n_wet_dom = 0
    n_dry_dom = 0
    median_diff = 0

# Load ecosystem data
input_file = output_dir / 'Q3_final_dataset.csv'
df_eco = pd.read_csv(input_file)
df_eco = df_eco.dropna(subset=['slope_theilsen'])
df_eco = df_eco[df_eco['ecosystem_type'].isin(['freshwater', 'saline', 'upland'])]
print(f"  Loaded {len(df_eco)} ecosystem sites")

# Order ecosystems by median (ascending)
eco_medians = df_eco.groupby('ecosystem_type')['slope_theilsen'].median().sort_values()
ecosystem_order = eco_medians.index.tolist()
print(f"  Ecosystem order by median: {ecosystem_order}")

# ============================================
# CAPPING FOR VISUALIZATION ONLY (95%)
# ============================================

# For Panel A (Dry-Wet) - cap diff values at 95%
if len(df_dw) > 0:
    diff_lower = df_dw['diff'].quantile(0.05)
    diff_upper = df_dw['diff'].quantile(0.95)
    df_dw['diff_capped'] = df_dw['diff'].clip(lower=diff_lower, upper=diff_upper)
    df_dw['is_extreme_low'] = df_dw['diff'] < diff_lower
    df_dw['is_extreme_high'] = df_dw['diff'] > diff_upper
    print(f"\n  Panel A capping (95%):")
    print(f"    Lower: {diff_lower:.4f}, Upper: {diff_upper:.4f}")
    print(f"    Extreme low: {df_dw['is_extreme_low'].sum()}, Extreme high: {df_dw['is_extreme_high'].sum()}")

# For Panel B (Ecosystem) - cap slope values at 95%
slope_lower = df_eco['slope_theilsen'].quantile(0.05)
slope_upper = df_eco['slope_theilsen'].quantile(0.95)
df_eco['slope_capped'] = df_eco['slope_theilsen'].clip(lower=slope_lower, upper=slope_upper)
df_eco['is_extreme_low'] = df_eco['slope_theilsen'] < slope_lower
df_eco['is_extreme_high'] = df_eco['slope_theilsen'] > slope_upper
print(f"\n  Panel B capping (95%):")
print(f"    Lower: {slope_lower:.4f}, Upper: {slope_upper:.4f}")
print(f"    Extreme low: {df_eco['is_extreme_low'].sum()}, Extreme high: {df_eco['is_extreme_high'].sum()}")

# ============================================
# CREATE FIGURE WITH TWO PANELS
# ============================================

print("\n🎨 Creating figure...")

fig, axes = plt.subplots(1, 2, figsize=(28, 14))

# ============================================
# PANEL A: DRY VS WET SENSITIVITY
# ============================================

ax1 = axes[0]

if len(df_dw) > 0:
    # Sort by diff
    df_dw_sorted = df_dw.sort_values('diff').reset_index(drop=True)
    diff_values_capped = df_dw_sorted['diff_capped'].values
    diff_values_true = df_dw_sorted['diff'].values
    is_extreme_low = df_dw_sorted['is_extreme_low'].values
    is_extreme_high = df_dw_sorted['is_extreme_high'].values
    
    # Get coast colors for each site
    site_colors = [coast_colors.get(coast, '#9E9E9E') for coast in df_dw_sorted['coast_region_analysis']]
    
    n_negative = (diff_values_true < 0).sum()
    n_positive = (diff_values_true > 0).sum()
    
    # Plot vertical lines (BLACK only - no green/gray)
    for i, diff_capped in enumerate(diff_values_capped):
        ax1.plot([i, i], [0, diff_capped], color='black', linewidth=2, alpha=0.5)
    
    # Plot points
    for i, (diff_capped, color, is_low, is_high) in enumerate(zip(diff_values_capped, site_colors, is_extreme_low, is_extreme_high)):
        if is_low:
            ax1.scatter(i, diff_lower, color=color, s=200, marker='v', zorder=5, edgecolors='none', alpha=0.85)
        elif is_high:
            ax1.scatter(i, diff_upper, color=color, s=200, marker='^', zorder=5, edgecolors='none', alpha=0.85)
        else:
            ax1.scatter(i, diff_capped, color=color, s=200, zorder=5, edgecolors='none', alpha=0.85)
    
    # Horizontal zero line
    ax1.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.7)
    
    # Vertical separator
    if n_negative > 0 and n_positive > 0:
        separator_pos = n_negative - 0.5
        ax1.axvline(x=separator_pos, color='black', linestyle='--', linewidth=2, alpha=0.6)
        
        # Add labels
        ax1.text(n_negative/2, ax1.get_ylim()[1] * 0.95, 'Dry-dominated\n(Wet < Dry)', 
                ha='center', va='top', fontsize=24, fontweight='bold', color='#666666')
        ax1.text(n_negative + n_positive/2.5, ax1.get_ylim()[1] * 0.99, 'Wet-dominated\n(Wet > Dry)', 
                ha='center', va='top', fontsize=24, fontweight='bold', color='#666666')
    
    # Axes formatting - INCREASED FONT SIZES
    ax1.set_xlabel('Sites ordered by dry–wet sensitivity difference', fontsize=32, fontweight='bold')
    ax1.set_ylabel('Δ Sensitivity (Wet − Dry)', fontsize=32, fontweight='bold')
    ax1.set_xticks([])
    ax1.tick_params(axis='y', labelsize=28, width=2, length=8)
    ax1.grid(True, alpha=0.3, axis='y', linestyle='--')
    ax1.set_axisbelow(True)
    
    # Make axis spines darker
    ax1.spines['bottom'].set_color('black')
    ax1.spines['bottom'].set_linewidth(3)
    ax1.spines['left'].set_color('black')
    ax1.spines['left'].set_linewidth(3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Panel label
    ax1.text(-0.08, 1.02, '(a)', transform=ax1.transAxes, fontsize=38, fontweight='bold', va='bottom')
    
    # Statistics - 2 significant digits, bottom right
    n_sites = len(df_dw)
    
    # Wilcoxon paired test using TRUE uncapped values
    wilcoxon_stat, wilcoxon_p = wilcoxon(df_dw['slope_wet_theilsen'], df_dw['slope_dry_theilsen'])
    
    # Format p-value with 2 significant digits
    if wilcoxon_p < 0.001:
        p_text = '<0.001'
    else:
        p_text = f'{wilcoxon_p:.2g}'
    
    stats_text = f"n = {n_sites}\n"
    stats_text += f"Wet > Dry: {n_wet_dom} ({n_wet_dom/n_sites*100:.0f}%)\n"
    stats_text += f"Dry > Wet: {n_dry_dom} ({n_dry_dom/n_sites*100:.0f}%)\n"
    stats_text += f"Median Δ = {median_diff:.2f}\n"
    stats_text += f"p = {p_text}"
    
    ax1.text(0.98, 0.05, stats_text, transform=ax1.transAxes, fontsize=24,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# PANEL B: ECOSYSTEM-LEVEL SENSITIVITY
# ============================================

ax2 = axes[1]

ecosystem_positions = range(len(ecosystem_order))

# Plot jittered points and median bars
for i, ecosystem in enumerate(ecosystem_order):
    eco_df = df_eco[df_eco['ecosystem_type'] == ecosystem]
    eco_slopes_true = eco_df['slope_theilsen'].values
    eco_slopes_capped = eco_df['slope_capped'].values
    eco_extreme_low = eco_df['is_extreme_low'].values
    eco_extreme_high = eco_df['is_extreme_high'].values
    
    if len(eco_slopes_true) > 0:
        # Calculate median from TRUE values (uncapped for stats)
        median_val = np.median(eco_slopes_true)
        
        # Add jitter
        jitter = np.random.normal(0, 0.12, size=len(eco_slopes_true))
        x_positions = np.ones(len(eco_slopes_true)) * i + jitter
        
        # Plot points with capping for visualization only
        for j, (x_pos, slope_capped, is_low, is_high) in enumerate(zip(
                x_positions, eco_slopes_capped, eco_extreme_low, eco_extreme_high)):
            
            if is_low:
                ax2.scatter(x_pos, slope_lower, color=ecosystem_colors[ecosystem], s=200, 
                           marker='v', alpha=0.7, edgecolors='none', zorder=3)
            elif is_high:
                ax2.scatter(x_pos, slope_upper, color=ecosystem_colors[ecosystem], s=200, 
                           marker='^', alpha=0.7, edgecolors='none', zorder=3)
            else:
                ax2.scatter(x_pos, slope_capped, color=ecosystem_colors[ecosystem], s=200, 
                           alpha=0.7, edgecolors='none', zorder=3)
        
        # Median bar (black)
        ax2.hlines(y=median_val, xmin=i - 0.35, xmax=i + 0.35, 
                  color='black', linewidth=7, zorder=4, alpha=0.9)

# Zero line
ax2.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.7)

# Axes formatting - INCREASED FONT SIZES
ax2.set_xlabel('Ecosystem type (ordered by median sensitivity)', fontsize=32, fontweight='bold')
ax2.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-6)', fontsize=32, fontweight='bold')
ax2.set_xticks(ecosystem_positions)
ax2.set_xticklabels([ecosystem_full_names[eco] for eco in ecosystem_order], fontsize=28)
ax2.set_ylim(slope_lower - 0.3, slope_upper + 0.3)
ax2.tick_params(axis='both', labelsize=28, width=2, length=8)
ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
ax2.set_axisbelow(True)

# Make axis spines darker
ax2.spines['bottom'].set_color('black')
ax2.spines['bottom'].set_linewidth(3)
ax2.spines['left'].set_color('black')
ax2.spines['left'].set_linewidth(3)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Panel label
ax2.text(-0.08, 1.02, '(b)', transform=ax2.transAxes, fontsize=38, fontweight='bold', va='bottom')

# Statistics - Kruskal-Wallis using TRUE uncapped values, 2 significant digits
ecosystem_groups = [df_eco[df_eco['ecosystem_type'] == eco]['slope_theilsen'].values 
                    for eco in ecosystem_order if len(df_eco[df_eco['ecosystem_type'] == eco]) > 0]
kw_stat, kw_p = kruskal(*ecosystem_groups)

# Format p-value with 2 significant digits
if kw_p < 0.001:
    p_text_b = '<0.001'
else:
    p_text_b = f'{kw_p:.2g}'

stats_text_b = f"Kruskal-Wallis p = {p_text_b}"
ax2.text(0.98, 0.05, stats_text_b, transform=ax2.transAxes, fontsize=26,
         verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# LEGEND FOR PANEL B (Ecosystem) - moved down
# ============================================

legend_elements_b = []
for ecosystem in ecosystem_order:
    count = len(df_eco[df_eco['ecosystem_type'] == ecosystem])
    legend_elements_b.append(
        Line2D([0], [0], marker='o', color='w', 
               label=f"{ecosystem_full_names[ecosystem]} (n={count})",
               markerfacecolor=ecosystem_colors[ecosystem], markersize=24,
               markeredgecolor='none')
    )

legend_b = ax2.legend(handles=legend_elements_b, loc='lower center', bbox_to_anchor=(0.5, -0.33),
                      ncol=3, fontsize=24, title='Ecosystem Types', title_fontsize=28, 
                      framealpha=0.95, edgecolor='gray', handlelength=2, handleheight=3, 
                      handletextpad=0.5)

# Color legend text
for i, ecosystem in enumerate(ecosystem_order):
    if i < len(legend_b.get_texts()):
        legend_b.get_texts()[i].set_color(ecosystem_colors[ecosystem])
        legend_b.get_texts()[i].set_fontweight('bold')

# ============================================
# LEGEND FOR PANEL A (Coast) - 2 rows, moved further down
# ============================================

legend_elements_a = []
for coast in coast_colors.keys():
    count = len(df_dw[df_dw['coast_region_analysis'] == coast]) if len(df_dw) > 0 else 0
    if count > 0:
        legend_elements_a.append(
            Line2D([0], [0], marker='o', color='w', 
                   label=f"{coast_full_names[coast]} (n={count})",
                   markerfacecolor=coast_colors[coast], markersize=24,
                   markeredgecolor='none')
        )

if len(legend_elements_a) > 0:
    legend_a = ax1.legend(handles=legend_elements_a, loc='lower center', bbox_to_anchor=(0.5, -0.35),
                          ncol=2, fontsize=22, title='Coastal Regions', title_fontsize=26, 
                          framealpha=0.95, edgecolor='gray', handlelength=2, handleheight=3,
                          handletextpad=0.5)
    
    # Color legend text
    for i, coast in enumerate(coast_colors.keys()):
        if i < len(legend_a.get_texts()):
            legend_a.get_texts()[i].set_color(coast_colors[coast])

# ============================================
# ADJUST LAYOUT
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.28, wspace=0.35)

# ============================================
# SAVE OUTPUTS
# ============================================

print("\n💾 Saving outputs...")

for ext in ['png', 'pdf']:
    output_file = output_dir / f'FigureS_Supplementary.{ext}'
    plt.savefig(output_file, dpi=300 if ext == 'png' else 150, bbox_inches='tight', facecolor='white')
    print(f"  ✓ Saved: {output_file}")

plt.show()

# ============================================
# FINAL SUMMARY
# ============================================

print("\n" + "="*70)
print("SUPPLEMENTARY FIGURE COMPLETE")
print("="*70)

print("\n📊 PANEL A (Dry-Wet Sensitivity):")
if len(df_dw) > 0:
    print(f"  • Total sites: {len(df_dw)}")
    print(f"  • Wet-dominated: {n_wet_dom} ({n_wet_dom/len(df_dw)*100:.1f}%)")
    print(f"  • Dry-dominated: {n_dry_dom} ({n_dry_dom/len(df_dw)*100:.1f}%)")
    print(f"  • Median difference: {median_diff:.4f}")
    print(f"  • Wilcoxon p-value: {p_text} (using uncapped values)")

print("\n📊 PANEL B (Ecosystem Sensitivity):")
for ecosystem in ecosystem_order:
    count = len(df_eco[df_eco['ecosystem_type'] == ecosystem])
    median = df_eco[df_eco['ecosystem_type'] == ecosystem]['slope_theilsen'].median()
    print(f"  • {ecosystem_full_names[ecosystem]}: n={count}, median={median:.4f}")
print(f"  • Kruskal-Wallis p-value: {p_text_b} (using uncapped values)")

# Add ecosystem note based on p-value
if kw_p > 0.05:
    print("\n📝 NOTE: Differences among ecosystem types were not statistically significant (p > 0.05).")
else:
    print("\n📝 NOTE: Differences among ecosystem types were statistically significant (p < 0.05).")

print("\n✅ Key features:")
print("  • Panel A: BLACK vertical lines, dots colored by coast region")
print("  • Panel B: Ecosystems ordered by median sensitivity")
print("  • Both panels: Capped at 95% for visualization (triangles for extremes)")
print("  • Statistical tests use UNCAPPED original values")
print("  • p-values shown with 2 significant digits")
print("  • Legends moved down to avoid x-axis overlap")
print("  • Increased font sizes for all labels and ticks")

print("\n" + "="*70)
print("✅ SUPPLEMENTARY FIGURE COMPLETE")
print(f"📁 Files saved to: {output_dir}")
print("="*70)