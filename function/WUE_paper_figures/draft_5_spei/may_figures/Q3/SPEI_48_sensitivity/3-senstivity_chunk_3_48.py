# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:25:03 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 12 10:47:53 2026

@author: ammar
"""

"""
CHUNK 3: Panel B - Coastal Region Sensitivity Distribution (SPEI-48)
- Full width across the page (bottom panel)
- SAME figure dimensions as Panel A (24, 12)
- X-axis grouped by coastal region (ordered by median sensitivity)
- POINTS for each site colored by SPEI-48 median gradient
- GRADIENT MEDIAN LINES: transition from region's min SPEI-48 to max SPEI-48
- TEXT LABELS on median line ends showing min and max SPEI-48 values
- BALANCED, ZERO-CENTERED COLORMAP (consistent with Panel A)
- SOFTER COLORS - less dominant blue on wet side
- REDUCED LINE WIDTH (thinner median lines)
- INCREASED DISTANCE between coastal regions (prevents overlap)
- Console table output for manuscript results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import kruskal
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================
# LOAD PROCESSED DATA (SPEI-48 FILTERED)
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
pickle_file = output_dir / 'processed_data_SPEI48_filtered.pkl'

with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df = data['df']
coast_order_ascending = data['coast_order_ascending']
coast_short_labels = data['coast_short_labels']
coast_full_names = data['coast_full_names']

# Get capping values (2nd and 97th percentiles - consistent with Panel A)
TOP_PERCENTILE = 97
BOTTOM_PERCENTILE = 2

cap_lower = df['slope_theilsen'].quantile(BOTTOM_PERCENTILE / 100)
cap_upper = df['slope_theilsen'].quantile(TOP_PERCENTILE / 100)

print("="*70)
print("CHUNK 3: Panel B - Coastal Region Sensitivity Distribution (SPEI-48)")
print("="*70)
print(f"\n⚙️ Filtering criteria: {data.get('filtering_criteria', 'Trans_ratio > 0')}")
print(f"⚙️ Capping settings (2nd/97th percentile - same as Panel A):")
print(f"   Bottom {BOTTOM_PERCENTILE}th percentile: {cap_lower:.4f}")
print(f"   Top {TOP_PERCENTILE}th percentile: {cap_upper:.4f}")

# Apply capping
df['slope_capped'] = df['slope_theilsen'].clip(lower=cap_lower, upper=cap_upper)
df['is_extreme_low'] = df['slope_theilsen'] < cap_lower
df['is_extreme_high'] = df['slope_theilsen'] > cap_upper

y_range = cap_upper - cap_lower
y_margin = y_range * 0.05

# ============================================
# CREATE SOFT, BALANCED COLORMAP (SAME AS PANEL A)
# ============================================

print("\n🎨 Creating balanced, zero-centered colormap for SPEI-48 (same as Panel A)...")

# Same colormap as Panel A for consistency
custom_cmap = LinearSegmentedColormap.from_list(
    "balanced_spei_gradient",
    [
        "#9E0142",  # stronger dry red
        "#D53E4F",
        "#F46D43",
        "#F7F7F7",  # near normal (white/gray)
        "#ABD9E9",
        "#74ADD1",
        "#2B83BA"
    ]
)

print("  Colormap: Balanced red-white-blue (centered at zero)")

# ============================================
# SETUP SPEI-48 NORMALIZATION (ZERO-CENTERED - SAME AS PANEL A)
# ============================================

# Get global SPEI-48 range from processed data
global_spei_min = df['spei48_min'].min()
global_spei_max = df['spei48_max'].max()

# Zero-centered normalization (scientifically justified for SPEI)
spei_absmax = max(abs(global_spei_min), abs(global_spei_max))

norm_spei = TwoSlopeNorm(
    vmin=-spei_absmax,
    vcenter=0,
    vmax=spei_absmax
)

print(f"  Global SPEI-48 range: {global_spei_min:.2f} to {global_spei_max:.2f}")
print(f"  Symmetric range for colorbar: -{spei_absmax:.2f} to +{spei_absmax:.2f}")
print(f"  Center at 0 (ecologically meaningful - drier/wetter than normal)")

# ============================================
# COMPUTE COASTAL REGION SPEI-48 RANGES
# ============================================

print("\n📊 Computing coastal region SPEI-48 ranges...")

coast_spei_ranges = {}
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    if len(coast_df) > 0:
        spei_min = coast_df['spei48_min'].min()
        spei_max = coast_df['spei48_max'].max()
        coast_spei_ranges[coast] = (spei_min, spei_max)
        print(f"  {coast_full_names.get(coast, coast)}: SPEI-48 [{spei_min:.2f}, {spei_max:.2f}]")

# ============================================
# CREATE FIGURE WITH INCREASED SPACING
# ============================================

fig, ax = plt.subplots(1, 1, figsize=(24, 12), dpi=300)

# Spacing between coastal regions (consistent with original)
coast_positions = [i * 2.5 for i in range(len(coast_order_ascending))]

print(f"\n📊 Creating Panel B (SPEI-48)...")
print(f"  Coastal regions: {len(coast_order_ascending)}")
print(f"  Region positions: {coast_positions}")
print(f"  Spacing increased to prevent overlap")

# ============================================
# PLOT POINTS AND GRADIENT MEDIAN LINES
# ============================================

for idx, coast in enumerate(coast_order_ascending):
    i = coast_positions[idx]
    coast_df = df[df['coast_region_analysis'] == coast]
    coast_slopes_true = coast_df['slope_theilsen'].values
    coast_is_sig = coast_df['is_significant'].values
    coast_extreme_low = coast_df['is_extreme_low'].values
    coast_extreme_high = coast_df['is_extreme_high'].values
    coast_spei_median = coast_df['spei48_median'].values  # Changed to spei48_median
    
    if len(coast_slopes_true) > 0:
        median_val = np.median(coast_slopes_true)
        jitter = np.random.normal(0, 0.12, size=len(coast_slopes_true))
        x_positions = np.ones(len(coast_slopes_true)) * i + jitter
        
        for j, (x_pos, slope_true, is_sig, is_low, is_high, spei_median) in enumerate(zip(
                x_positions, coast_slopes_true, coast_is_sig, 
                coast_extreme_low, coast_extreme_high, coast_spei_median)):
            
            point_color = custom_cmap(norm_spei(spei_median))
            
            if is_low:
                y_plot = cap_lower
                marker = 'v'
                size = 190
            elif is_high:
                y_plot = cap_upper
                marker = '^'
                size = 190
            else:
                y_plot = slope_true
                marker = 'o'
                size = 190
            
            ax.scatter(x_pos, y_plot, color=point_color, s=size, marker=marker,
                       alpha=0.9, edgecolors='none', zorder=3)
            
            if is_sig:
                if slope_true < 0:
                    star_y = y_plot - y_range * 0.05
                    va = 'top'
                else:
                    star_y = y_plot + y_range * 0.05
                    va = 'bottom'
                
                ax.text(x_pos, star_y, '*', fontsize=20, fontweight='bold',
                        color='black', ha='center', va=va, zorder=10)
        
        # ============================================
        # GRADIENT MEDIAN LINE (THINNER)
        # ============================================
        spei_min, spei_max = coast_spei_ranges.get(coast, (-2, 2))
        
        n_segments = 100
        line_length = 0.9
        start_x = i - 0.45
        
        for seg in range(n_segments):
            t = seg / n_segments
            spei_value = spei_min + t * (spei_max - spei_min)
            line_color = custom_cmap(norm_spei(spei_value))
            
            x_start = start_x + t * line_length
            x_end = start_x + (seg + 1) / n_segments * line_length
            
            ax.hlines(y=median_val, xmin=x_start, xmax=x_end,
                      color=line_color, linewidth=2.5, zorder=4, alpha=0.95)
        
        # ============================================
        # TEXT LABELS (MIN/MAX SPEI-48)
        # ============================================
        min_color = custom_cmap(norm_spei(spei_min))
        ax.text(i - 0.48, median_val, f'{spei_min:.1f}', 
                ha='right', va='center', fontsize=30, fontweight='bold',
                color=min_color, alpha=0.9, zorder=5)
        
        max_color = custom_cmap(norm_spei(spei_max))
        ax.text(i + 0.48, median_val, f'{spei_max:.1f}', 
                ha='left', va='center', fontsize=30, fontweight='bold',
                color=max_color, alpha=0.9, zorder=5)
        
        # Sample size annotation
        n_sites_region = len(coast_df)
        y_text_pos = cap_upper + y_range * 0.22
        ax.text(i, y_text_pos, f'N = {n_sites_region}',
               ha='center', va='bottom', fontsize=34, fontweight='bold',
               color='black', alpha=0.8)

# ============================================
# REFERENCE LINES
# ============================================

ax.axhline(y=0, color='#666666', linestyle=':', linewidth=2.5, alpha=0.8, zorder=1)

# ============================================
# AXES FORMATTING
# ============================================

ax.set_ylim(cap_lower - y_margin - y_range * 0.08, 
            cap_upper + y_margin + y_range * 0.12)

ax.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-48)', 
              fontsize=30, fontweight='bold', labelpad=20)

ax.set_xlabel('Coastal region (ordered by median sensitivity)', 
              fontsize=42, fontweight='bold', labelpad=15)

ax.set_xticks(coast_positions)
ax.set_xticklabels([coast_short_labels[coast] for coast in coast_order_ascending], 
                   fontsize=40, fontweight='bold')

for tick in ax.get_xticklabels():
    tick.set_color('black')

ax.tick_params(axis='both', labelsize=34, width=2, length=10)
ax.grid(True, alpha=0.15, axis='y', linestyle='--', zorder=0)
ax.set_axisbelow(True)

for spine in ax.spines.values():
    spine.set_linewidth(3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ============================================
# ADD PANEL LABEL (b)
# ============================================

ax.text(-0.02, 1.08, '(b)', transform=ax.transAxes, 
        fontsize=54, fontweight='bold', va='bottom', ha='left')

# ============================================
# ADD STATISTICS (Kruskal-Wallis)
# ============================================

coast_groups = [df[df['coast_region_analysis'] == coast]['slope_theilsen'].values 
                for coast in coast_order_ascending if len(df[df['coast_region_analysis'] == coast]) > 0]
kw_stat, kw_p = kruskal(*coast_groups)

def format_p(p_val):
    if p_val < 0.001:
        return '<0.001'
    else:
        return f'{p_val:.2f}'

stats_text = f"Kruskal-Wallis p = {format_p(kw_p)}"

ax.text(0.98, 0.05, stats_text, transform=ax.transAxes, fontsize=32,
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=2))

# ============================================
# COLORBAR (SYMMETRIC, ZERO-CENTERED - SAME AS PANEL A)
# ============================================

print("\n🎨 Adding symmetric, zero-centered colorbar for SPEI-48 (same as Panel A)...")

sm = ScalarMappable(cmap=custom_cmap, norm=norm_spei)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, orientation='horizontal',
                    fraction=0.05, pad=0.20, aspect=50)

cbar.set_label('SPEI-48 (drier ← normal → wetter)', 
               fontsize=34, fontweight='bold', labelpad=15)
cbar.ax.tick_params(labelsize=30)

# Use standard ticks (consistent with Panel A)
max_tick = min(3, np.ceil(spei_absmax))
ticks = np.arange(-max_tick, max_tick + 0.1, 1)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{t:.0f}' for t in ticks])

print(f"  Colorbar symmetric range: -{spei_absmax:.2f} to +{spei_absmax:.2f}")
print(f"  Colorbar ticks: {ticks}")

# ============================================
# SAVE FIGURE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.95, left=0.08, right=0.95)

png_file = output_dir / 'PanelB_SPEI48_sensitivity_points.png'
pdf_file = output_dir / 'PanelB_SPEI48_sensitivity_points.pdf'

plt.savefig(png_file, dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved Panel B: {png_file}")
print(f"✓ Saved Panel B: {pdf_file}")

plt.show()

# ============================================
# CONSOLE TABLE OUTPUT FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT TABLE: Coastal Region Sensitivity Results (SPEI-48)")
print("="*80)

results_table = []
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    n_total = len(coast_df)
    n_sig = coast_df['is_significant'].sum()
    median_val = coast_df['slope_theilsen'].median()
    q25_val = coast_df['slope_theilsen'].quantile(0.25)
    q75_val = coast_df['slope_theilsen'].quantile(0.75)
    mean_spei = coast_df['spei48_median'].mean()  # Changed to spei48_median
    spei_min, spei_max = coast_spei_ranges.get(coast, (0, 0))
    coast_name = coast_full_names.get(coast, coast)
    
    results_table.append({
        'Coast': coast_name,
        'N': n_total,
        'Significant (p<0.05)': n_sig,
        'Median Slope': f"{median_val:.3f}",
        'IQR': f"{q25_val:.3f}–{q75_val:.3f}",
        'Mean SPEI-48': f"{mean_spei:.2f}",  # Changed label
        'SPEI-48 Range': f"[{spei_min:.2f}, {spei_max:.2f}]"  # Changed label
    })

results_df = pd.DataFrame(results_table)
print("\n" + results_df.to_string(index=False))

table_csv = output_dir / 'PanelB_SPEI48_manuscript_results_table.csv'
results_df.to_csv(table_csv, index=False)
print(f"\n✅ Saved manuscript table to: {table_csv}")

# ============================================
# PRINT SUMMARY STATISTICS
# ============================================

print("\n" + "="*70)
print("PANEL B SUMMARY STATISTICS (SPEI-48)")
print("="*70)

print("\n  Coast regions (ordered by median sensitivity):")
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    n_total = len(coast_df)
    n_sig = coast_df['is_significant'].sum()
    median_val = coast_df['slope_theilsen'].median()
    mean_spei = coast_df['spei48_median'].mean()  # Changed to spei48_median
    spei_min, spei_max = coast_spei_ranges.get(coast, (0, 0))
    coast_name = coast_full_names.get(coast, coast)
    pct_sig = (n_sig / n_total * 100) if n_total > 0 else 0
    print(f"    {coast_name}: N={n_total}, sig={n_sig} ({pct_sig:.0f}%), median={median_val:.3f}, mean SPEI-48={mean_spei:.2f}, SPEI-48 range=[{spei_min:.2f}, {spei_max:.2f}]")

# ============================================
# CAPTION NOTE (printed to console, NOT on figure)
# ============================================

print("\n" + "="*80)
print("CAPTION NOTE FOR MANUSCRIPT (do not put on figure)")
print("="*80)
print("\nLine colors show the observed SPEI-48 range for each coastal region,")
print("with red indicating drier periods, blue indicating wetter periods,")
print("and near-white indicating near-normal conditions; colors provide")
print("hydroclimatic context and do not define the direction of WUE_T sensitivity.")

print("\n✅ UPDATES FOR SPEI-48:")
print("  • Converted from SPEI-6 to SPEI-48")
print("  • Using filtered data (Trans_ratio > 0)")
print("  • Same zero-centered colormap as Panel A (balanced red-white-blue)")
print("  • Same TwoSlopeNorm normalization as Panel A")
print("  • Symmetric colorbar centered at 0")
print("  • All SPEI ranges now from SPEI-48 data")
print("  • Updated labels to SPEI-48")
print(f"  • Filtering criteria: {data.get('filtering_criteria', 'Trans_ratio > 0')}")

print("\n" + "="*70)
print("CHUNK 3 COMPLETE - Panel B saved (SPEI-48 with consistent colors)")
print("="*70)