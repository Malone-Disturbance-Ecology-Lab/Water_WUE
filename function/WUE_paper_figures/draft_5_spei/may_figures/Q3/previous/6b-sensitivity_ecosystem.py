"""
CHUNK 2: WUE_T Sensitivity by Ecosystem/Salinity Class - Figure Generation
Purpose: Create single publication-quality figure with slope values only
- Sites ordered by median slope within ecosystem
- Legend at bottom right (single legend with stats)
- Capped at 2nd and 97th percentiles for visualization
- Stars placed with padding to avoid overlap
- No panel label (a)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pickle
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# CONFIGURATION
# ============================================

OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Figure dimensions - WIDER to prevent site name overlap
FIGURE_WIDTH = 20  # Increased from 14
FIGURE_HEIGHT = 10  # Slightly taller
DPI = 300

# Font sizes
AXIS_LABEL_FONT = 18
AXIS_TICK_FONT = 10
ECOSYSTEM_LABEL_FONT = 16
LEGEND_FONT = 11

# Boxplot styling
BOX_WIDTH = 0.6
BOX_ALPHA = 0.7
BOX_LINEWIDTH = 1.5
MEDIAN_LINEWIDTH = 2.5
WHISKER_WIDTH = 1.5
POINT_SIZE = 45
POINT_ALPHA = 0.6
JITTER_WIDTH = 0.15

# Capping percentiles for visualization
CAP_LOWER_PERCENTILE = 2
CAP_UPPER_PERCENTILE = 97

# ============================================
# LOAD PREPARED DATA
# ============================================

print("="*80)
print("CHUNK 2: Figure Generation (Single Panel)")
print("="*80)

pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data.pkl'
with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df_clean = data['df_clean']
ECOSYSTEM_COLORS = data['ECOSYSTEM_COLORS']
ECOSYSTEM_ORDER = data['ECOSYSTEM_ORDER']
groups_present = data['groups_present']
stats = data['statistics']

print(f"\n📊 Loaded data for {len(df_clean)} sites")
print(f"   Groups present: {groups_present}")

# ============================================
# APPLY CAPPING FOR VISUALIZATION ONLY
# ============================================

print("\n" + "="*60)
print("APPLYING CAPPING FOR VISUALIZATION")
print("="*60)

cap_lower = df_clean['slope'].quantile(CAP_LOWER_PERCENTILE / 100)
cap_upper = df_clean['slope'].quantile(CAP_UPPER_PERCENTILE / 100)

print(f"   Lower cap ({CAP_LOWER_PERCENTILE}th percentile): {cap_lower:.3f}")
print(f"   Upper cap ({CAP_UPPER_PERCENTILE}th percentile): {cap_upper:.3f}")
print(f"   Sites below cap: {(df_clean['slope'] < cap_lower).sum()}")
print(f"   Sites above cap: {(df_clean['slope'] > cap_upper).sum()}")

# Create capped version for visualization only
df_clean['slope_capped'] = df_clean['slope'].clip(lower=cap_lower, upper=cap_upper)

# ============================================
# ORDER SITES BY MEDIAN SLOPE WITHIN ECOSYSTEM
# ============================================

print("\n" + "="*60)
print("ORDERING SITES FOR VISUALIZATION")
print("="*60)

# Calculate median slope per site (should already have one per site)
site_slopes = df_clean.groupby(['site_name', 'Salinity_Category'])['slope_capped'].first().reset_index()

# Order within each ecosystem
ordered_sites = []
ordered_ecosystems = []

for eco in ECOSYSTEM_ORDER:
    eco_sites = site_slopes[site_slopes['Salinity_Category'] == eco].copy()
    # Sort by capped slope value
    eco_sites_sorted = eco_sites.sort_values('slope_capped')
    ordered_sites.extend(eco_sites_sorted['site_name'].tolist())
    ordered_ecosystems.extend([eco] * len(eco_sites_sorted))
    print(f"   {eco}: {len(eco_sites_sorted)} sites")

print(f"\nTotal ordered sites: {len(ordered_sites)}")

# Create positions with increased spacing
positions = np.arange(len(ordered_sites)) * 1.2  # 20% more space between sites

# ============================================
# CREATE FIGURE
# ============================================

print("\n" + "="*60)
print("CREATING FIGURE")
print("="*60)

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

# Colors for each site
box_colors = [ECOSYSTEM_COLORS[eco] for eco in ordered_ecosystems]

# Create data lists for boxplot (each site = one value since site-level data)
site_slope_values = []
for site in ordered_sites:
    slope_val = df_clean[df_clean['site_name'] == site]['slope_capped'].iloc[0]
    site_slope_values.append([slope_val])

# Create boxplot
bp = ax.boxplot(site_slope_values, positions=positions, widths=BOX_WIDTH,
                patch_artist=True, showfliers=False, vert=True)

# Color boxes
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(box_colors[i])
    patch.set_alpha(BOX_ALPHA)
    patch.set_edgecolor('black')
    patch.set_linewidth(BOX_LINEWIDTH)

# Style medians
for median in bp['medians']:
    median.set_color('red')
    median.set_linewidth(MEDIAN_LINEWIDTH)

# Style whiskers and caps
for whisker in bp['whiskers']:
    whisker.set_color('black')
    whisker.set_linewidth(WHISKER_WIDTH)
for cap in bp['caps']:
    cap.set_color('black')
    cap.set_linewidth(WHISKER_WIDTH)

# ============================================
# ADD JITTERED POINTS WITH CAPPED VALUES
# ============================================

# Calculate y-range for star placement
y_min, y_max = cap_lower, cap_upper
y_range_plot = y_max - y_min
star_offset = y_range_plot * 0.03  # 3% offset for stars

for i, (site, eco) in enumerate(zip(ordered_sites, ordered_ecosystems)):
    slope_val = df_clean[df_clean['site_name'] == site]['slope_capped'].iloc[0]
    original_slope = df_clean[df_clean['site_name'] == site]['slope'].iloc[0]
    
    # Determine if point is extreme (capped)
    is_capped_low = original_slope < cap_lower
    is_capped_high = original_slope > cap_upper
    
    jitter = np.random.normal(0, JITTER_WIDTH)
    x_pos = positions[i] + jitter
    
    # Different markers for capped vs normal
    if is_capped_low:
        marker = 'v'  # Downward triangle for low outliers
        size = POINT_SIZE
    elif is_capped_high:
        marker = '^'  # Upward triangle for high outliers
        size = POINT_SIZE
    else:
        marker = 'o'
        size = POINT_SIZE
    
    point_color = ECOSYSTEM_COLORS[eco]
    ax.scatter(x_pos, slope_val, color=point_color, s=size, 
              marker=marker, alpha=POINT_ALPHA, edgecolors='black', 
              linewidth=0.8, zorder=3)
    
    # Add asterisk for significant slopes with padding to avoid overlap
    is_sig = False
    if 'spearman_p' in df_clean.columns:
        p_val = df_clean[df_clean['site_name'] == site]['spearman_p'].iloc[0]
        is_sig = p_val < 0.05
    elif 'pearson_p' in df_clean.columns:
        p_val = df_clean[df_clean['site_name'] == site]['pearson_p'].iloc[0]
        is_sig = p_val < 0.05
    
    if is_sig:
        # Position star above or below point with padding
        if slope_val >= 0:
            star_y = slope_val + star_offset
            va = 'bottom'
        else:
            star_y = slope_val - star_offset
            va = 'top'
        
        ax.text(x_pos, star_y, '*', fontsize=16, fontweight='bold',
               ha='center', va=va, color='black', zorder=4)

# Add horizontal line at zero
ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.7, zorder=1)

# ============================================
# ECOSYSTEM LABELS AND SEPARATORS
# ============================================

# Calculate ecosystem boundaries (accounting for increased spacing)
ecosystem_boundaries = []
current_eco = None
start_idx = 0

# Convert to indices based on ordered lists
for i, eco in enumerate(ordered_ecosystems):
    if current_eco is None:
        current_eco = eco
        start_idx = i
    elif eco != current_eco:
        ecosystem_boundaries.append((current_eco, start_idx, i - 1))
        current_eco = eco
        start_idx = i
if current_eco:
    ecosystem_boundaries.append((current_eco, start_idx, len(ordered_sites) - 1))

# Add ecosystem labels at top
y_top = ax.get_ylim()[1]
y_text_pos = y_top + (y_range_plot * 0.08)  # 8% above top

for eco, start_idx, end_idx in ecosystem_boundaries:
    center_x = (positions[start_idx] + positions[end_idx]) / 2
    ax.text(center_x, y_text_pos, eco, ha='center', va='bottom',
           fontsize=ECOSYSTEM_LABEL_FONT, fontweight='bold', 
           color=ECOSYSTEM_COLORS[eco])

# Add vertical separator lines
for eco, start_idx, end_idx in ecosystem_boundaries:
    if start_idx > 0:
        separator_x = (positions[start_idx] + positions[start_idx - 1]) / 2
        ax.axvline(x=separator_x, color='black', linestyle='--', 
                  alpha=0.5, linewidth=1.5, zorder=2)

# ============================================
# FORMATTING
# ============================================

# X-axis: site names (rotated for readability with increased spacing)
ax.set_xticks(positions)
ax.set_xticklabels(ordered_sites, rotation=45, ha='right', 
                  fontsize=AXIS_TICK_FONT)

# Y-axis
ax.set_ylabel('WUE$_T$ Sensitivity Slope\n(ΔWUE$_T$ / ΔSPEI-6)', 
             fontsize=AXIS_LABEL_FONT, fontweight='bold', labelpad=15)
ax.tick_params(axis='y', labelsize=AXIS_TICK_FONT + 2)

# Adjust y-limits to accommodate stars and labels
y_min_plot = cap_lower - (y_range_plot * 0.1)
y_max_plot = y_text_pos + (y_range_plot * 0.05)
ax.set_ylim(y_min_plot, y_max_plot)

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Grid
ax.grid(True, linestyle=':', alpha=0.3, axis='y', zorder=0)
ax.set_axisbelow(True)

# ============================================
# SINGLE LEGEND AT BOTTOM RIGHT
# ============================================

# Create legend elements
legend_elements = [
    # Ecosystem colors (simple rectangles without text labels since labels are on plot)
    Rectangle((0, 0), 1, 1, facecolor=ECOSYSTEM_COLORS['Upland'], 
             edgecolor='black', alpha=BOX_ALPHA, label='_nolegend_'),
    Rectangle((0, 0), 1, 1, facecolor=ECOSYSTEM_COLORS['Freshwater'], 
             edgecolor='black', alpha=BOX_ALPHA, label='_nolegend_'),
    Rectangle((0, 0), 1, 1, facecolor=ECOSYSTEM_COLORS['Brackish'], 
             edgecolor='black', alpha=BOX_ALPHA, label='_nolegend_'),
    Rectangle((0, 0), 1, 1, facecolor=ECOSYSTEM_COLORS['Saline'], 
             edgecolor='black', alpha=BOX_ALPHA, label='_nolegend_'),
    
    # Point markers
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
              markeredgecolor='black', markersize=8, label='Normal value'),
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='gray',
              markeredgecolor='black', markersize=8, label=f'< {CAP_LOWER_PERCENTILE}th pctl'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='gray',
              markeredgecolor='black', markersize=8, label=f'> {CAP_UPPER_PERCENTILE}th pctl'),
    
    # Other elements
    plt.Line2D([0], [0], color='red', linestyle='--', linewidth=2, 
              label='Zero line'),
    plt.Line2D([0], [0], color='black', marker='*', markersize=12, 
              linestyle='None', label='Significant (p < 0.05)'),
]

# Add statistics to legend
if stats['kw_p'] is not None:
    if stats['kw_p'] < 0.001:
        kw_text = f"Kruskal-Wallis p < 0.001"
    elif stats['kw_p'] < 0.01:
        kw_text = f"Kruskal-Wallis p = {stats['kw_p']:.3f}"
    else:
        kw_text = f"Kruskal-Wallis p = {stats['kw_p']:.3f}"
    
    # Add as text element in legend
    legend_elements.append(plt.Line2D([0], [0], color='none', label=kw_text))

# Add sample size info
total_n = len(df_clean)
legend_elements.append(plt.Line2D([0], [0], color='none', label=f'N = {total_n} sites'))

# Place legend at bottom right
legend = ax.legend(handles=legend_elements, loc='lower right', 
                  fontsize=LEGEND_FONT, frameon=True, fancybox=True, 
                  shadow=True, ncol=2, 
                  bbox_to_anchor=(1.0, -0.15))  # Position below x-axis

# ============================================
# ADD SAMPLE SIZES BELOW ECOSYSTEM GROUPS
# ============================================

# Position sample sizes below x-axis
for eco, start_idx, end_idx in ecosystem_boundaries:
    n_sites = end_idx - start_idx + 1
    center_x = (positions[start_idx] + positions[end_idx]) / 2
    y_pos = ax.get_ylim()[0] - (y_range_plot * 0.05)
    ax.text(center_x, y_pos, f'n = {n_sites}', ha='center', va='top', 
           fontsize=AXIS_TICK_FONT + 1, fontweight='bold', color='black')

# ============================================
# ADJUST LAYOUT AND SAVE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.18, top=0.92, left=0.08, right=0.95)

# Save figure
png_file = OUTPUT_DIR / 'Fig_Diagnostic_Sensitivity_by_Ecosystem.png'
pdf_file = OUTPUT_DIR / 'Fig_Diagnostic_Sensitivity_by_Ecosystem.pdf'

plt.savefig(png_file, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {png_file}")
print(f"✅ Saved: {pdf_file}")

plt.show()

# ============================================
# PRINT CAPTION FOR FIGURE
# ============================================

print("\n" + "="*80)
print("FIGURE CAPTION")
print("="*80)

print("\nFigure X: Site-level WUE_T sensitivity slopes to SPEI-6 grouped by ecosystem/salinity class.")
print(f"Sites (n={len(ordered_sites)}) are ordered by increasing sensitivity within each group. ")
print(f"Boxes show the interquartile range with median (red line), whiskers extend to the full range ")
print(f"of non-capped values. Points are jittered individual site values, with triangles indicating ")
print(f"values below the {CAP_LOWER_PERCENTILE}th percentile (▼) or above the {CAP_UPPER_PERCENTILE}th percentile (▲). ")
print(f"Black asterisks denote statistically significant slopes (p < 0.05). ")
print(f"Vertical dashed lines separate ecosystem groups. ")

if stats['kw_p'] is not None:
    if stats['kw_p'] < 0.001:
        print(f"Kruskal-Wallis test: H = {stats['kw_stat']:.2f}, p < 0.001, ε² = {stats['epsilon_sq']:.3f}")
    else:
        print(f"Kruskal-Wallis test: H = {stats['kw_stat']:.2f}, p = {stats['kw_p']:.4f}, ε² = {stats['epsilon_sq']:.3f}")

print("\n" + "="*80)
print("CHUNK 2 COMPLETE")
print("="*80)