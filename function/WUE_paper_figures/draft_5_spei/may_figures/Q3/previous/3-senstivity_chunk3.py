"""
CHUNK 3: Panel B - Coastal Region Sensitivity Distribution (UPDATED)
- Full width across the page (bottom panel)
- SAME figure dimensions as Panel A (24, 12)
- X-axis grouped by coastal region (ordered by median sensitivity)
- POINTS for each site colored by SPEI median gradient
- GRADIENT MEDIAN LINES: transition from region's min SPEI to max SPEI
- TEXT LABELS on median line ends showing min and max SPEI values (closer to lines)
- REDUCED LINE WIDTH (thinner median lines)
- INCREASED DISTANCE between coastal regions (prevents overlap)
- Console table output for manuscript results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import kruskal
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================
# LOAD PROCESSED DATA
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")
pickle_file = output_dir / 'processed_data.pkl'

with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

df = data['df']
coast_order_ascending = data['coast_order_ascending']
coast_short_labels = data['coast_short_labels']
coast_full_names = data['coast_full_names']

# Get capping values
TOP_PERCENTILE = 97
BOTTOM_PERCENTILE = 2

cap_lower = df['slope_theilsen'].quantile(BOTTOM_PERCENTILE / 100)
cap_upper = df['slope_theilsen'].quantile(TOP_PERCENTILE / 100)

print("="*70)
print("CHUNK 3: Panel B - Coastal Region Sensitivity Distribution (UPDATED)")
print("="*70)
print(f"\n⚙️ Capping settings:")
print(f"   Bottom {BOTTOM_PERCENTILE}th percentile: {cap_lower:.4f}")
print(f"   Top {TOP_PERCENTILE}th percentile: {cap_upper:.4f}")

# Apply capping
df['slope_capped'] = df['slope_theilsen'].clip(lower=cap_lower, upper=cap_upper)
df['is_extreme_low'] = df['slope_theilsen'] < cap_lower
df['is_extreme_high'] = df['slope_theilsen'] > cap_upper

y_range = cap_upper - cap_lower
y_margin = y_range * 0.05

# ============================================
# CREATE ENHANCED COLORMAP (HIGHER CONTRAST)
# ============================================

print("\n🎨 Creating enhanced colormap with higher contrast...")

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

print("  Colormap: FireBrick → Crimson → Light Pink → Light Blue → Deep Navy")

# ============================================
# SETUP SPEI NORMALIZATION
# ============================================

global_spei_min = df['spei6_min'].min()
global_spei_max = df['spei6_max'].max()
norm_spei = Normalize(vmin=global_spei_min, vmax=global_spei_max)

print(f"  Global SPEI range: {global_spei_min:.2f} to {global_spei_max:.2f}")

# ============================================
# COMPUTE COASTAL REGION SPEI RANGES
# ============================================

print("\n📊 Computing coastal region SPEI ranges...")

coast_spei_ranges = {}
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    if len(coast_df) > 0:
        spei_min = coast_df['spei6_min'].min()
        spei_max = coast_df['spei6_max'].max()
        coast_spei_ranges[coast] = (spei_min, spei_max)
        print(f"  {coast_full_names.get(coast, coast)}: SPEI [{spei_min:.2f}, {spei_max:.2f}]")

# ============================================
# CREATE FIGURE WITH INCREASED SPACING
# ============================================

fig, ax = plt.subplots(1, 1, figsize=(24, 12), dpi=300)

# INCREASED spacing between coastal regions
# Original: x_positions = i (0,1,2,3,4) with default spacing
# Now: Space them further apart by using i * 1.5
coast_positions = [i * 2.5 for i in range(len(coast_order_ascending))]

print(f"\n📊 Creating Panel B...")
print(f"  Coastal regions: {len(coast_order_ascending)}")
print(f"  Region positions: {coast_positions}")
print(f"  Spacing increased by 50% to prevent overlap")

# ============================================
# PLOT POINTS AND GRADIENT MEDIAN LINES
# ============================================

for idx, coast in enumerate(coast_order_ascending):
    i = coast_positions[idx]  # Use spaced position
    coast_df = df[df['coast_region_analysis'] == coast]
    coast_slopes_true = coast_df['slope_theilsen'].values
    coast_is_sig = coast_df['is_significant'].values
    coast_extreme_low = coast_df['is_extreme_low'].values
    coast_extreme_high = coast_df['is_extreme_high'].values
    coast_spei_median = coast_df['spei6_median'].values
    
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
        # GRADIENT MEDIAN LINE (REDUCED LINE WIDTH)
        # ============================================
        spei_min, spei_max = coast_spei_ranges.get(coast, (-2, 2))
        
        # REDUCED line width from 4.5 to 2.5
        n_segments = 100
        line_length = 0.9  # Slightly longer line
        start_x = i - 0.45
        
        for seg in range(n_segments):
            t = seg / n_segments
            spei_value = spei_min + t * (spei_max - spei_min)
            line_color = custom_cmap(norm_spei(spei_value))
            
            x_start = start_x + t * line_length
            x_end = start_x + (seg + 1) / n_segments * line_length
            
            ax.hlines(y=median_val, xmin=x_start, xmax=x_end,
                      color=line_color, linewidth=2.5, zorder=4, alpha=0.95)  # Reduced from 4.5 to 2.5
        
        # ============================================
        # TEXT LABELS (MOVED CLOSER TO LINE ENDS)
        # ============================================
        # Min SPEI label - moved closer (offset from 0.52 to 0.48)
        min_color = custom_cmap(norm_spei(spei_min))
        ax.text(i - 0.48, median_val, f'{spei_min:.1f}', 
                ha='right', va='center', fontsize=30, fontweight='bold',
                color=min_color, alpha=0.9, zorder=5)
        
        # Max SPEI label - moved closer (offset from 0.52 to 0.48)
        max_color = custom_cmap(norm_spei(spei_max))
        ax.text(i + 0.48, median_val, f'{spei_max:.1f}', 
                ha='left', va='center', fontsize=30, fontweight='bold',
                color=max_color, alpha=0.9, zorder=5)
        
        # Sample size annotation (slightly adjusted position)
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

ax.set_ylabel('Sensitivity slope (ΔWUE$_T$ / ΔSPEI-6)', 
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
# ADD STATISTICS
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
# COLORBAR
# ============================================

print("\n🎨 Adding colorbar for SPEI gradient...")

sm = ScalarMappable(cmap=custom_cmap, norm=norm_spei)
sm.set_array([])

cbar = fig.colorbar(sm, ax=ax, orientation='horizontal',
                    fraction=0.05, pad=0.20, aspect=50)

cbar.set_label('SPEI-6 (drier → wetter)', 
               fontsize=34, fontweight='bold', labelpad=15)
cbar.ax.tick_params(labelsize=30)

cbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
cbar.set_ticklabels(['-3', '-2', '-1', '0', '1', '2', '3'])

# ============================================
# SAVE FIGURE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.95, left=0.08, right=0.95)  # Adjusted right margin

png_file = output_dir / 'PanelB_Q3_sensitivity_points.png'
pdf_file = output_dir / 'PanelB_Q3_sensitivity_points.pdf'

plt.savefig(png_file, dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved Panel B: {png_file}")
print(f"✓ Saved Panel B: {pdf_file}")

plt.show()

# ============================================
# CONSOLE TABLE OUTPUT FOR MANUSCRIPT
# ============================================

print("\n" + "="*80)
print("MANUSCRIPT TABLE: Coastal Region Sensitivity Results")
print("="*80)

results_table = []
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    n_total = len(coast_df)
    n_sig = coast_df['is_significant'].sum()
    median_val = coast_df['slope_theilsen'].median()
    q25_val = coast_df['slope_theilsen'].quantile(0.25)
    q75_val = coast_df['slope_theilsen'].quantile(0.75)
    mean_spei = coast_df['spei6_median'].mean()
    spei_min, spei_max = coast_spei_ranges.get(coast, (0, 0))
    coast_name = coast_full_names.get(coast, coast)
    
    results_table.append({
        'Coast': coast_name,
        'N': n_total,
        'Significant (p<0.05)': n_sig,
        'Median Slope': f"{median_val:.3f}",
        'IQR': f"{q25_val:.3f}–{q75_val:.3f}",
        'Mean SPEI': f"{mean_spei:.2f}",
        'SPEI Range': f"[{spei_min:.2f}, {spei_max:.2f}]"
    })

results_df = pd.DataFrame(results_table)
print("\n" + results_df.to_string(index=False))

table_csv = output_dir / 'PanelB_manuscript_results_table.csv'
results_df.to_csv(table_csv, index=False)
print(f"\n✅ Saved manuscript table to: {table_csv}")

# ============================================
# PRINT SUMMARY STATISTICS
# ============================================

print("\n" + "="*70)
print("PANEL B SUMMARY STATISTICS")
print("="*70)

print("\n  Coast regions (ordered by median sensitivity):")
for coast in coast_order_ascending:
    coast_df = df[df['coast_region_analysis'] == coast]
    n_total = len(coast_df)
    n_sig = coast_df['is_significant'].sum()
    median_val = coast_df['slope_theilsen'].median()
    mean_spei = coast_df['spei6_median'].mean()
    spei_min, spei_max = coast_spei_ranges.get(coast, (0, 0))
    coast_name = coast_full_names.get(coast, coast)
    pct_sig = (n_sig / n_total * 100) if n_total > 0 else 0
    print(f"    {coast_name}: N={n_total}, sig={n_sig} ({pct_sig:.0f}%), median={median_val:.3f}, mean SPEI={mean_spei:.2f}, SPEI range=[{spei_min:.2f}, {spei_max:.2f}]")

print("\n✅ UPDATES APPLIED:")
print("  • INCREASED COASTAL REGION SPACING: 50% wider (positions: 0, 1.5, 3.0, 4.5, 6.0)")
print("  • REDUCED MEDIAN LINE WIDTH: from 4.5 to 2.5 (thinner lines)")
print("  • TEXT LABELS MOVED CLOSER: offset from 0.52 to 0.48")
print("  • PREVENTED OVERLAP: numbers now closer to respective line ends")
print("  • ENHANCED COLORMAP: 9 distinct color steps")
print("  • Console table output saved as CSV")

print("\n" + "="*70)
print("CHUNK 3 COMPLETE - Panel B saved")
print("="*70)