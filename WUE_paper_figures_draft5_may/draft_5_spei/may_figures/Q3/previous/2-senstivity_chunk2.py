"""
CHUNK 2: Panel A - Site-level Sensitivity Distribution (UPDATED)
- Full width across the page (top panel)
- Sites ordered by sensitivity slope
- Vertical lines show SPEI GRADIENT based on full range of SPEI values per site
- SAME ENHANCED COLORMAP as Panel B (9 distinct color steps)
- COLLECTIVE HORIZONTAL LINES (THICK) for positive and negative sensitivity groups at TOP
- Text labels showing min/max SPEI for each group (larger font)
- Positive/Negative sensitivity labels: ABOVE horizontal lines, BLACK color, NO (N)
- Positive label moved RIGHT to avoid crossing vertical separator
- Black stars for significant slopes (p < 0.05)
- Comprehensive console output for manuscript
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy import stats
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================
# USER-ADJUSTABLE PARAMETERS
# ============================================

TOP_PERCENTILE = 97
BOTTOM_PERCENTILE = 2

# ============================================
# LOAD PROCESSED DATA
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
pickle_file = output_dir / 'processed_data.pkl'

with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df = data['df']
df_sorted = data['df_sorted']
global_neg_95 = data['global_neg_95']
global_pos_95 = data['global_pos_95']
y_range = data['y_range']
y_margin = data['y_margin']

print("="*70)
print("CHUNK 2: Panel A - Site-level Sensitivity Distribution (UPDATED)")
print("="*70)
print(f"\n⚙️ Capping settings:")
print(f"   Bottom {BOTTOM_PERCENTILE}th percentile: {df['slope_theilsen'].quantile(BOTTOM_PERCENTILE/100):.4f}")
print(f"   Top {TOP_PERCENTILE}th percentile: {df['slope_theilsen'].quantile(TOP_PERCENTILE/100):.4f}")

# ============================================
# APPLY CAPPING
# ============================================

cap_lower = df['slope_theilsen'].quantile(BOTTOM_PERCENTILE / 100)
cap_upper = df['slope_theilsen'].quantile(TOP_PERCENTILE / 100)

df['slope_capped'] = df['slope_theilsen'].clip(lower=cap_lower, upper=cap_upper)
df['is_extreme_low'] = df['slope_theilsen'] < cap_lower
df['is_extreme_high'] = df['slope_theilsen'] > cap_upper

df_sorted = df.sort_values('slope_capped').reset_index(drop=True)

slopes_capped = df_sorted['slope_capped'].values
y_range = cap_upper - cap_lower
y_margin = y_range * 0.05

# ============================================
# CREATE SAME ENHANCED COLORMAP AS PANEL B
# ============================================

print("\n🎨 Creating same enhanced colormap as Panel B...")

custom_cmap = LinearSegmentedColormap.from_list(
    "enhanced_gradient",
    [
        "#B22222",  # FireBrick - very dry (< -2)
        "#DC143C",  # Crimson - dry
        "#F08080",  # LightCoral - moderately dry
        "#FADADD",  # Very light pink - near zero negative
        "#D4E6F1",  # Very light blue - near zero positive
        "#5DADE2",  # Light blue - moderately wet
        "#2874A6",  # Medium blue - wet
        "#1A5276",  # Dark blue - very wet (> 2)
        "#0B3B60"   # Deep Navy - extreme wet
    ]
)

# ============================================
# CREATE FIGURE
# ============================================

fig, ax = plt.subplots(1, 1, figsize=(24, 12), dpi=300)

slopes = df_sorted['slope_capped'].values
n_sites = len(slopes)
n_negative = (slopes < 0).sum()
n_positive = (slopes > 0).sum()

print(f"\n📊 Creating Panel A...")
print(f"  Total sites: {n_sites}")
print(f"  Negative slopes: {n_negative}")
print(f"  Positive slopes: {n_positive}")

# ============================================
# GET SPEI RANGE FOR EACH SITE
# ============================================

monthly_filtered_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
df_monthly = pd.read_csv(monthly_filtered_path)

spei_range_dict = {}
for site in df_sorted['site_name'].unique():
    site_data = df_monthly[df_monthly['site_name'] == site].dropna(subset=['SPEI_6', 'WUE_tra'])
    if len(site_data) > 0:
        spei_min = site_data['SPEI_6'].min()
        spei_max = site_data['SPEI_6'].max()
        spei_range_dict[site] = (spei_min, spei_max)
    else:
        spei_range_dict[site] = (-2, 2)

global_spei_min = df['spei6_min'].min()
global_spei_max = df['spei6_max'].max()
norm_spei = Normalize(vmin=global_spei_min, vmax=global_spei_max)

print(f"  Global SPEI range: {global_spei_min:.2f} to {global_spei_max:.2f}")

# ============================================
# CALCULATE GROUP SPEI RANGES
# ============================================

negative_sites = df_sorted[df_sorted['slope_capped'] < 0]['site_name'].values
positive_sites = df_sorted[df_sorted['slope_capped'] > 0]['site_name'].values

neg_spei_min = min([spei_range_dict.get(site, (-2, 2))[0] for site in negative_sites]) if len(negative_sites) > 0 else 0
neg_spei_max = max([spei_range_dict.get(site, (-2, 2))[1] for site in negative_sites]) if len(negative_sites) > 0 else 0
pos_spei_min = min([spei_range_dict.get(site, (-2, 2))[0] for site in positive_sites]) if len(positive_sites) > 0 else 0
pos_spei_max = max([spei_range_dict.get(site, (-2, 2))[1] for site in positive_sites]) if len(positive_sites) > 0 else 0

print(f"\n  Negative sensitivity group (n={len(negative_sites)}): SPEI range [{neg_spei_min:.2f}, {neg_spei_max:.2f}]")
print(f"  Positive sensitivity group (n={len(positive_sites)}): SPEI range [{pos_spei_min:.2f}, {pos_spei_max:.2f}]")

# ============================================
# PLOT SITES WITH GRADIENT VERTICAL LINES
# ============================================

for i, (slope, row) in enumerate(zip(slopes, df_sorted.iterrows())):
    row_data = row[1]
    site_name = row_data['site_name']
    spei_median = row_data['spei6_median']
    is_sig = row_data['is_significant']
    
    spei_min, spei_max = spei_range_dict.get(site_name, (-2, 2))
    
    n_segments = 50
    for seg in range(n_segments):
        t = seg / n_segments
        spei_value = spei_min + t * (spei_max - spei_min)
        line_color = custom_cmap(norm_spei(spei_value))
        
        y_start = 0 + t * slope
        y_end = 0 + (seg + 1) / n_segments * slope
        
        ax.plot([i, i], [y_start, y_end], color=line_color, 
                linewidth=5.0, alpha=1.0, zorder=2, solid_capstyle='butt')
    
    if row_data['is_extreme_low']:
        y_plot = cap_lower
        marker = 'v'
    elif row_data['is_extreme_high']:
        y_plot = cap_upper
        marker = '^'
    else:
        y_plot = slope
        marker = 'o'
    
    point_color = custom_cmap(norm_spei(spei_median))
    
    ax.scatter(i, y_plot, color=point_color, s=190, marker=marker,
               zorder=5, edgecolors='none', alpha=1.0)
    
    if is_sig:
        if slope < 0:
            star_y = y_plot - y_range * 0.03
            va = 'top'
        else:
            star_y = y_plot + y_range * 0.03
            va = 'bottom'
        
        ax.text(i, star_y, '*', fontsize=20, fontweight='bold',
                color='black', ha='center', va=va, zorder=10)

# ============================================
# ADD COLLECTIVE HORIZONTAL LINES AT TOP
# ============================================

# Position lines at top
y_upper_limit = cap_upper + y_margin + y_range * 0.15
y_neg_line_pos = y_upper_limit * 0.94
y_pos_line_pos = y_upper_limit * 0.94  # Same height for both

line_length_multiplier = 0.40

# Horizontal line for NEGATIVE sensitivity group
if n_negative > 0:
    line_length_neg = n_negative * line_length_multiplier
    start_x_neg = 0 + 2
    n_segments = 100
    for seg in range(n_segments):
        t = seg / n_segments
        spei_value = neg_spei_min + t * (neg_spei_max - neg_spei_min)
        line_color = custom_cmap(norm_spei(spei_value))
        
        x_start = start_x_neg + t * line_length_neg
        x_end = start_x_neg + (seg + 1) / n_segments * line_length_neg
        
        ax.hlines(y=y_neg_line_pos, xmin=x_start, xmax=x_end,
                  color=line_color, linewidth=4.5, zorder=4, alpha=0.95)
    
    neg_min_color = custom_cmap(norm_spei(neg_spei_min))
    neg_max_color = custom_cmap(norm_spei(neg_spei_max))
    
    ax.text(start_x_neg - 0.45, y_neg_line_pos, f'{neg_spei_min:.1f}', 
            ha='right', va='center', fontsize=26, fontweight='bold',
            color=neg_min_color, alpha=0.9, zorder=5)
    
    ax.text(start_x_neg + line_length_neg + 0.45, y_neg_line_pos, f'{neg_spei_max:.1f}', 
            ha='left', va='center', fontsize=26, fontweight='bold',
            color=neg_max_color, alpha=0.9, zorder=5)

# Horizontal line for POSITIVE sensitivity group (moved RIGHT to avoid crossing separator)
if n_positive > 0:
    # Start position shifted right to avoid crossing the vertical separator line
    start_x_pos = n_negative + 4  # Moved right (was n_negative - 0.35)
    line_length_pos = n_positive * line_length_multiplier
    
    n_segments = 100
    for seg in range(n_segments):
        t = seg / n_segments
        spei_value = pos_spei_min + t * (pos_spei_max - pos_spei_min)
        line_color = custom_cmap(norm_spei(spei_value))
        
        x_start = start_x_pos + t * line_length_pos
        x_end = start_x_pos + (seg + 1) / n_segments * line_length_pos
        
        ax.hlines(y=y_pos_line_pos, xmin=x_start, xmax=x_end,
                  color=line_color, linewidth=4.5, zorder=4, alpha=0.95)
    
    pos_min_color = custom_cmap(norm_spei(pos_spei_min))
    pos_max_color = custom_cmap(norm_spei(pos_spei_max))
    
    ax.text(start_x_pos - 0.45, y_pos_line_pos, f'{pos_spei_min:.1f}', 
            ha='right', va='center', fontsize=26, fontweight='bold',
            color=pos_min_color, alpha=0.9, zorder=5)
    
    ax.text(start_x_pos + line_length_pos + 0.45, y_pos_line_pos, f'{pos_spei_max:.1f}', 
            ha='left', va='center', fontsize=26, fontweight='bold',
            color=pos_max_color, alpha=0.9, zorder=5)

# ============================================
# REFERENCE LINES
# ============================================

ax.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.8, zorder=1)

# ============================================
# VERTICAL SEPARATOR WITH LABELS (ABOVE horizontal lines)
# ============================================

if n_negative > 0 and n_positive > 0:
    separator_pos = n_negative - 0.5
    ax.axvline(x=separator_pos, color='black', linestyle='--', linewidth=2, alpha=0.6, zorder=1)
    
    # Position labels ABOVE the horizontal lines
    label_y_pos = y_neg_line_pos + y_range * 0.06
    
    # Negative label (centered over negative region)
    neg_label_x = n_negative / 2.0
    ax.text(neg_label_x, label_y_pos, 'Negative\nsensitivity', 
            ha='center', va='bottom', fontsize=32, fontweight='bold', color='black')
    
    # Positive label (moved RIGHT to center over positive region)
    pos_label_x = n_negative + n_positive / 2.0
    ax.text(pos_label_x, label_y_pos, 'Positive\nsensitivity', 
            ha='center', va='bottom', fontsize=32, fontweight='bold', color='black')

# ============================================
# AXES FORMATTING
# ============================================

ax.set_ylim(cap_lower - y_margin - y_range * 0.08, 
            label_y_pos + y_range * 0.08)

ax.set_xlabel('Sites ordered by WUE$_T$ sensitivity', fontsize=42, fontweight='bold', labelpad=15)
ax.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-6)', fontsize=35, fontweight='bold', labelpad=15)
ax.set_xticks([])
ax.tick_params(axis='both', labelsize=30, width=2, length=10)
ax.grid(True, alpha=0.15, axis='y', linestyle='--', zorder=0)
ax.set_axisbelow(True)

for spine in ax.spines.values():
    spine.set_linewidth(3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ============================================
# ADD PANEL LABEL (a)
# ============================================

ax.text(0.0, 1.08, '(a)', transform=ax.transAxes, 
        fontsize=54, fontweight='bold', va='bottom', ha='left')

# ============================================
# ADD STATISTICS (BOTTOM-RIGHT CORNER)
# ============================================

n_sites_total = len(df)
n_pos_total = (df['slope_theilsen'] > 0).sum()
n_neg_total = (df['slope_theilsen'] < 0).sum()
median_slope = df['slope_theilsen'].median()
iqr_slope = df['slope_theilsen'].quantile(0.75) - df['slope_theilsen'].quantile(0.25)
wilcoxon_stat, wilcoxon_p = stats.wilcoxon(df['slope_theilsen'])
n_sig_total = df['is_significant'].sum()

def format_p(p_val):
    if p_val < 0.001:
        return '<0.001'
    else:
        return f'{p_val:.2f}'

def format_num(x):
    if abs(x) < 0.001:
        return '≈0'
    else:
        return f'{x:.2f}'

stats_text = (f"N = {n_sites_total}\n"
              f"Positive: {n_pos_total} ({n_pos_total/n_sites_total*100:.0f}%)\n"
              f"Negative: {n_neg_total} ({n_neg_total/n_sites_total*100:.0f}%)\n"
              f"Significant: {n_sig_total}\n"
              f"Median = {format_num(median_slope)}\n"
              f"IQR = {format_num(iqr_slope)}\n"
              f"p = {format_p(wilcoxon_p)}")

ax.text(0.98, 0.05, stats_text, transform=ax.transAxes, fontsize=25,
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# COLORBAR
# ============================================

print("\n🎨 Adding colorbar for SPEI gradient...")

sm = ScalarMappable(cmap=custom_cmap, norm=norm_spei)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, orientation='horizontal',
                    fraction=0.05, pad=0.12, aspect=50)

cbar.set_label('SPEI-6 (drier → wetter)', 
               fontsize=34, fontweight='bold', labelpad=15)
cbar.ax.tick_params(labelsize=30)

cbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
cbar.set_ticklabels(['-3', '-2', '-1', '0', '1', '2', '3'])

# ============================================
# SAVE PANEL A
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.95)

png_file = output_dir / 'PanelA_Q3_sensitivity_gradient.png'
pdf_file = output_dir / 'PanelA_Q3_sensitivity_gradient.pdf'

plt.savefig(png_file, dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved Panel A: {png_file}")
print(f"✓ Saved Panel A: {pdf_file}")

plt.show()

# ============================================
# COMPREHENSIVE CONSOLE OUTPUT FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT RESULTS - PANEL A")
print("="*80)

print("\n📊 SAMPLE CHARACTERISTICS:")
print(f"   • Total sites analyzed: {n_sites_total}")
print(f"   • Sites with positive sensitivity: {n_pos_total} ({n_pos_total/n_sites_total*100:.1f}%)")
print(f"   • Sites with negative sensitivity: {n_neg_total} ({n_neg_total/n_sites_total*100:.1f}%)")
print(f"   • Statistically significant slopes (p < 0.05): {n_sig_total} ({n_sig_total/n_sites_total*100:.1f}%)")

print("\n📊 SLOPE DISTRIBUTION STATISTICS:")
print(f"   • Median slope: {format_num(median_slope)}")
print(f"   • Interquartile range (IQR): {format_num(iqr_slope)}")
print(f"   • Slope range: {df['slope_theilsen'].min():.2f} to {df['slope_theilsen'].max():.2f}")
print(f"   • Wilcoxon signed-rank test (one-sample, H0: median = 0): p = {format_p(wilcoxon_p)}")
print(f"   • Interpretation: {'Significant' if wilcoxon_p < 0.05 else 'Not significant'} deviation from zero")

print("\n📊 HYDROCLIMATIC CONTEXT (SPEI-6):")
if len(negative_sites) > 0:
    print(f"   • Negative sensitivity sites (n={len(negative_sites)}):")
    print(f"     - SPEI range: [{neg_spei_min:.2f}, {neg_spei_max:.2f}]")
    print(f"     - Interpretation: Sites span from {'dry' if neg_spei_min < 0 else 'wet'} to {'wet' if neg_spei_max > 0 else 'dry'} conditions")
if len(positive_sites) > 0:
    print(f"   • Positive sensitivity sites (n={len(positive_sites)}):")
    print(f"     - SPEI range: [{pos_spei_min:.2f}, {pos_spei_max:.2f}]")
    print(f"     - Interpretation: Sites span from {'dry' if pos_spei_min < 0 else 'wet'} to {'wet' if pos_spei_max > 0 else 'dry'} conditions")

print("\n📊 SALINITY BREAKDOWN (for reference):")
# Get salinity from the original dataframe (not from the sorted one)
if 'Salinity_Category' in df.columns:
    for salinity in ['Freshwater', 'Saline', 'Upland', 'Brackish']:
        count = df[df['Salinity_Category'] == salinity].shape[0]
        if count > 0:
            print(f"   • {salinity}: {count} sites")
else:
    print("   • Salinity data not available in current dataframe")

print("\n✅ UPDATES APPLIED:")
print("  • Fixed KeyError for Salinity_Category")
print("  • Negative/Positive sensitivity labels moved ABOVE horizontal lines")
print("  • Positive label line moved RIGHT (start_x_pos = n_negative + 0.15)")
print("  • Labels no longer cross vertical separator line")
print("  • Both sensitivity labels at same height")

print("\n" + "="*70)
print("CHUNK 2 COMPLETE - Panel A saved")
print("="*70)