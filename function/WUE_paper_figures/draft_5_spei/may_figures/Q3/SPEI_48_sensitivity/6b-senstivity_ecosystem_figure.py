# -*- coding: utf-8 -*-
"""
CHUNK 2: WUE_T Sensitivity by Ecosystem/Salinity Class - Figure Generation (SPEI-48)
Purpose: Create publication-quality figure with points only (no boxplots since one value per site)
- Sites ordered by median slope within ecosystem
- Points colored by ecosystem (fill and edge color match ecosystem)
- Black stars for significant slopes
- Horizontal gradient lines showing SPEI-48 range (min to max) for each ecosystem
- Small vertical colorbar showing SPEI gradient (dry → wet)
- Sample size shown next to ecosystem label
- Only p-value and N in legend box
- Vertical dashed lines between ecosystem groups
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import pickle
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================
# CONFIGURATION
# ============================================

OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april\diagnostics_SPEI48")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Figure dimensions
FIGURE_WIDTH = 30
FIGURE_HEIGHT = 15
DPI = 600

# Font sizes
AXIS_LABEL_FONT = 40
AXIS_TICK_FONT = 20
ECOSYSTEM_LABEL_FONT = 40
LEGEND_FONT = 30
GRADIENT_LABEL_FONT = 24
COLORBAR_FONT = 32

# Point styling
POINT_SIZE = 300
POINT_ALPHA = 0.9
JITTER_WIDTH = 0.2

# Star offset
STAR_OFFSET_FACTOR = 0.015

# Capping percentiles for visualization (internal only)
CAP_LOWER_PERCENTILE = 2
CAP_UPPER_PERCENTILE = 97

# Number of significant digits
SIG_DIGITS = 2

# Line length multipliers (shorter for Brackish to prevent overlap)
LINE_LENGTH_MULTIPLIER = {
    'Upland': 0.7,
    'Freshwater': 0.6,
    'Brackish': 0.4,
    'Saline': 0.6
}

# Horizontal offset for line centering (negative = left, positive = right)
LINE_CENTER_OFFSET = {
    'Upland': 0,
    'Freshwater': -0.8,
    'Brackish': 0,
    'Saline': 0.8
}

# ============================================
# LOAD PREPARED DATA (SPEI-48)
# ============================================

print("="*80)
print("CHUNK 2: Figure Generation (SPEI-48) - Points Only with SPEI Gradient")
print("="*80)

pickle_file = OUTPUT_DIR / 'ecosystem_diagnostic_data_SPEI48.pkl'
with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df_clean = data['df_clean']
ECOSYSTEM_COLORS = data['ECOSYSTEM_COLORS']
ECOSYSTEM_ORDER = data['ECOSYSTEM_ORDER']
groups_present = data['groups_present']
stats = data['statistics']
spei_version = data.get('spei_version', 'SPEI-48')

print(f"\n📊 Loaded data for {len(df_clean)} sites ({spei_version})")
print(f"   Groups present: {groups_present}")

# ============================================
# CREATE SPEI GRADIENT COLORMAP (dry to wet)
# ============================================

print("\n🎨 Creating SPEI gradient colormap (dry → wet)...")

spei_gradient_cmap = LinearSegmentedColormap.from_list(
    "spei_gradient",
    [
        "#B22222",  # FireBrick - dry (< -2)
        "#DC143C",  # Crimson - dry
        "#F08080",  # LightCoral - moderately dry
        "#FADADD",  # Very light pink - near zero negative
        "#D4E6F1",  # Very light blue - near zero positive
        "#5DADE2",  # Light blue - moderately wet
        "#2874A6",  # Medium blue - wet
        "#1A5276",  # Dark blue - very wet (> 2)
    ]
)

# ============================================
# APPLY CAPPING FOR VISUALIZATION (internal only)
# ============================================

cap_lower = df_clean['slope'].quantile(CAP_LOWER_PERCENTILE / 100)
cap_upper = df_clean['slope'].quantile(CAP_UPPER_PERCENTILE / 100)

# Create capped version for visualization only
df_clean['slope_capped'] = df_clean['slope'].clip(lower=cap_lower, upper=cap_upper)

print(f"\n   Capping applied ({CAP_LOWER_PERCENTILE}-{CAP_UPPER_PERCENTILE}th percentile)")

# ============================================
# CALCULATE SPEI-48 RANGE FOR EACH ECOSYSTEM
# ============================================

print("\n📊 Calculating SPEI-48 range per ecosystem...")

ecosystem_spei_ranges = {}
for eco in ECOSYSTEM_ORDER:
    eco_df = df_clean[df_clean['Salinity_Category'] == eco]
    if len(eco_df) > 0:
        spei_min = eco_df['spei48_min'].min()
        spei_max = eco_df['spei48_max'].max()
        ecosystem_spei_ranges[eco] = (spei_min, spei_max)
        print(f"   {eco}: SPEI-48 [{spei_min:.2f}, {spei_max:.2f}]")

# Global SPEI range for colorbar
global_spei_min = df_clean['spei48_min'].min()
global_spei_max = df_clean['spei48_max'].max()
norm_spei = Normalize(vmin=global_spei_min, vmax=global_spei_max)

print(f"   Global SPEI-48 range: [{global_spei_min:.2f}, {global_spei_max:.2f}]")

# ============================================
# ORDER SITES BY MEDIAN SLOPE WITHIN ECOSYSTEM
# ============================================

print("\n" + "="*60)
print("ORDERING SITES FOR VISUALIZATION")
print("="*60)

# Order within each ecosystem
ordered_sites = []
ordered_ecosystems = []
site_slope_dict = {}
site_sig_dict = {}

for eco in ECOSYSTEM_ORDER:
    eco_sites = df_clean[df_clean['Salinity_Category'] == eco].copy()
    eco_sites_sorted = eco_sites.sort_values('slope_capped')
    
    for _, row in eco_sites_sorted.iterrows():
        ordered_sites.append(row['site_name'])
        ordered_ecosystems.append(eco)
        site_slope_dict[row['site_name']] = row['slope_capped']
        
        is_sig = False
        if 'is_significant' in df_clean.columns:
            is_sig = bool(row['is_significant'])
        site_sig_dict[row['site_name']] = is_sig
    
    print(f"   {eco}: {len(eco_sites_sorted)} sites")

print(f"\nTotal ordered sites: {len(ordered_sites)}")

# Create positions with increased spacing
positions = np.arange(len(ordered_sites)) * 1.3

# ============================================
# CREATE FIGURE
# ============================================

print("\n" + "="*60)
print("CREATING FIGURE (SPEI-48 - Points with Gradient Lines)")
print("="*60)

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))

# ============================================
# ADD JITTERED POINTS
# ============================================

y_range_all = []

for i, (site, eco) in enumerate(zip(ordered_sites, ordered_ecosystems)):
    slope_val = site_slope_dict[site]
    y_range_all.append(slope_val)
    
    jitter = np.random.normal(0, JITTER_WIDTH)
    x_pos = positions[i] + jitter
    
    point_color = ECOSYSTEM_COLORS[eco]
    
    ax.scatter(x_pos, slope_val, color=point_color, s=POINT_SIZE, 
              marker='o', alpha=POINT_ALPHA, edgecolors=point_color, 
              linewidth=0, zorder=3)

# ============================================
# ADD BLACK STARS FOR SIGNIFICANT SLOPES
# ============================================

y_min_vals = min(y_range_all) if y_range_all else cap_lower
y_max_vals = max(y_range_all) if y_range_all else cap_upper
y_range_plot_star = y_max_vals - y_min_vals
star_offset = y_range_plot_star * STAR_OFFSET_FACTOR

for i, (site, eco) in enumerate(zip(ordered_sites, ordered_ecosystems)):
    slope_val = site_slope_dict[site]
    is_sig = site_sig_dict[site]
    
    if is_sig:
        jitter = np.random.normal(0, JITTER_WIDTH)
        x_pos = positions[i] + jitter
        
        if slope_val >= 0:
            star_y = slope_val + star_offset
            va = 'bottom'
        else:
            star_y = slope_val - star_offset
            va = 'top'
        
        ax.text(x_pos, star_y, '*', fontsize=24, fontweight='bold',
               ha='center', va=va, color='black', zorder=4)

# ============================================
# ADD HORIZONTAL ZERO LINE
# ============================================

ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.8, zorder=1)

# ============================================
# ECOSYSTEM LABELS AND HORIZONTAL GRADIENT LINES
# ============================================

# Calculate y-position for labels (top of figure)
y_max_vals_adj = max(y_range_all) if y_range_all else cap_upper
y_range_plot = y_max_vals_adj - cap_lower
y_text_pos = y_max_vals_adj + (y_range_plot * 0.12)

# Calculate ecosystem boundaries
ecosystem_boundaries = []
current_eco = None
start_idx = 0

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

# Define y-position for gradient lines (slightly below ecosystem label)
y_line_pos = y_text_pos - (y_range_plot * 0.04)

# Add ecosystem labels and horizontal gradient lines with adjusted lengths/positions
for eco, start_idx, end_idx in ecosystem_boundaries:
    n_sites = end_idx - start_idx + 1
    ecosystem_center_x = (positions[start_idx] + positions[end_idx]) / 2
    
    # Get SPEI range for this ecosystem
    spei_min, spei_max = ecosystem_spei_ranges.get(eco, (-2, 2))
    
    # Get line length multiplier and center offset for this ecosystem
    line_mult = LINE_LENGTH_MULTIPLIER.get(eco, 0.7)
    center_offset = LINE_CENTER_OFFSET.get(eco, 0)
    
    # Calculate line length based on ecosystem width and multiplier
    ecosystem_width = positions[end_idx] - positions[start_idx]
    line_length = ecosystem_width * line_mult
    
    # Calculate starting position with offset
    start_x = ecosystem_center_x + center_offset - (line_length / 2)
    
    # Draw horizontal gradient line
    n_segments = 100
    
    for seg in range(n_segments):
        t = seg / n_segments
        spei_value = spei_min + t * (spei_max - spei_min)
        line_color = spei_gradient_cmap(norm_spei(spei_value))
        
        x_start_seg = start_x + t * line_length
        x_end_seg = start_x + (seg + 1) / n_segments * line_length
        
        ax.hlines(y=y_line_pos, xmin=x_start_seg, xmax=x_end_seg,
                  color=line_color, linewidth=4.5, zorder=4, alpha=0.95)
    
    # Add min/max SPEI labels on the line ends
    min_color = spei_gradient_cmap(norm_spei(spei_min))
    max_color = spei_gradient_cmap(norm_spei(spei_max))
    
    label_offset = 0.25 if line_mult > 0.5 else 0.15
    
    ax.text(start_x - label_offset, y_line_pos, f'{spei_min:.1f}', 
            ha='right', va='center', fontsize=GRADIENT_LABEL_FONT, 
            fontweight='bold', color=min_color, alpha=0.9, zorder=5)
    
    ax.text(start_x + line_length + label_offset, y_line_pos, f'{spei_max:.1f}', 
            ha='left', va='center', fontsize=GRADIENT_LABEL_FONT, 
            fontweight='bold', color=max_color, alpha=0.9, zorder=5)
    
    # Add ecosystem label (above the gradient line)
    label_text = f"{eco}\n(N={n_sites})"
    ax.text(ecosystem_center_x, y_text_pos, label_text, ha='center', va='bottom',
           fontsize=ECOSYSTEM_LABEL_FONT, fontweight='bold', 
           color=ECOSYSTEM_COLORS[eco])

# ============================================
# ADD VERTICAL DASHED LINES BETWEEN ECOSYSTEM GROUPS
# ============================================

for eco, start_idx, end_idx in ecosystem_boundaries:
    if start_idx > 0:
        separator_x = (positions[start_idx] + positions[start_idx - 1]) / 2
        ax.axvline(x=separator_x, color='black', linestyle='--', 
                  alpha=0.6, linewidth=1.5, zorder=2)

# ============================================
# ADD SMALL VERTICAL COLORBAR (SPEI LEGEND)
# ============================================

print("\n🎨 Adding small vertical colorbar for SPEI legend...")

# Create a small inset axes for the colorbar
cbar_ax = fig.add_axes([0.92, 0.25, 0.025, 0.45])  # [left, bottom, width, height]

# Create colorbar
sm = ScalarMappable(cmap=spei_gradient_cmap, norm=norm_spei)
sm.set_array([])

cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
cbar.set_label('SPEI-48\n(drier → wetter)', 
               fontsize=COLORBAR_FONT, fontweight='bold', labelpad=10)
cbar.ax.tick_params(labelsize=COLORBAR_FONT - 1)

# Set reasonable ticks
tick_values = np.linspace(global_spei_min, global_spei_max, 5)
cbar.set_ticks(tick_values)
cbar.set_ticklabels([f'{t:.1f}' for t in tick_values])

print(f"   Colorbar range: {global_spei_min:.1f} to {global_spei_max:.1f}")

# ============================================
# FORMATTING
# ============================================

# X-axis: site names
ax.set_xticks(positions)
ax.set_xticklabels(ordered_sites, rotation=45, ha='right', 
                  fontsize=AXIS_TICK_FONT)

# Y-axis
ax.set_ylabel(f'WUE$_T$ Sensitivity to {spei_version}\n(ΔWUE$_T$ / Δ{spei_version})', 
             fontsize=AXIS_LABEL_FONT, fontweight='bold', labelpad=15)

# Format y-axis ticks
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.{SIG_DIGITS}f}'))
ax.tick_params(axis='y', labelsize=AXIS_TICK_FONT + 10, width=1.5, length=8)

# Adjust y-limits
y_min_plot = cap_lower - (y_range_plot * 0.1)
y_max_plot = y_text_pos + (y_range_plot * 0.08)
ax.set_ylim(y_min_plot, y_max_plot)

# X-axis ticks formatting
ax.tick_params(axis='x', labelsize=AXIS_TICK_FONT, rotation=45)

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Grid
ax.grid(True, linestyle=':', alpha=0.3, axis='y', zorder=0)
ax.set_axisbelow(True)

# ============================================
# SIMPLIFIED LEGEND (p-value and N only)
# ============================================

legend_elements = []

if stats['kw_p'] is not None:
    if stats['kw_p'] < 0.001:
        kw_text = f"Kruskal-Wallis p < 0.001"
    elif stats['kw_p'] < 0.01:
        kw_text = f"p = {stats['kw_p']:.2f}"
    else:
        kw_text = f"p = {stats['kw_p']:.2f}"
    
    legend_elements.append(plt.Line2D([0], [0], color='none', label=kw_text))

total_n = len(df_clean)
legend_elements.append(plt.Line2D([0], [0], color='none', label=f'N = {total_n} sites'))

legend = ax.legend(handles=legend_elements, loc='lower right', 
                  fontsize=LEGEND_FONT, frameon=True, fancybox=True, 
                  shadow=True, ncol=1)

# ============================================
# ADJUST LAYOUT AND SAVE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.90, left=0.08, right=0.91)  # Adjusted right margin for colorbar

# Save figure
png_file = OUTPUT_DIR / 'Fig_Diagnostic_Sensitivity_by_Ecosystem_SPEI48_with_gradient.png'
pdf_file = OUTPUT_DIR / 'Fig_Diagnostic_Sensitivity_by_Ecosystem_SPEI48_with_gradient.pdf'

plt.savefig(png_file, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=DPI, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {png_file}")
print(f"✅ Saved: {pdf_file}")

plt.show()

# ============================================
# PRINT CAPTION FOR FIGURE
# ============================================

print("\n" + "="*80)
print("FIGURE CAPTION (SPEI-48)")
print("="*80)

print(f"\nFigure X: Site-level WUE_T sensitivity slopes to {spei_version} grouped by ecosystem/salinity class.")
print(f"Sites (n={len(ordered_sites)}) are ordered by increasing sensitivity within each group. ")
print(f"Points represent individual site values, colored by ecosystem class. ")
print(f"Black asterisks (*) indicate statistically significant slopes (p < 0.05). ")
print(f"Horizontal gradient lines below each ecosystem label show the range of SPEI-48 values")
print(f"observed at sites within that ecosystem, with red indicating drier conditions and")
print(f"blue indicating wetter conditions. The small vertical colorbar on the right shows the")
print(f"SPEI-48 color mapping. The dashed red line indicates zero sensitivity. ")

if stats['kw_p'] is not None:
    if stats['kw_p'] < 0.001:
        print(f"Kruskal-Wallis test: H = {stats['kw_stat']:.2f}, p < 0.001, ε² = {stats['epsilon_sq']:.2f}")
    else:
        print(f"Kruskal-Wallis test: H = {stats['kw_stat']:.2f}, p = {stats['kw_p']:.4f}, ε² = {stats['epsilon_sq']:.2f}")

print("\n" + "="*80)
print("CHUNK 2 COMPLETE (SPEI-48 - Points with Gradient Lines + Colorbar)")
print("="*80)
print(f"\n📁 Output saved to: {OUTPUT_DIR}")
print(f"   - PNG: {png_file.name}")
print(f"   - PDF: {pdf_file.name}")
print("="*80)