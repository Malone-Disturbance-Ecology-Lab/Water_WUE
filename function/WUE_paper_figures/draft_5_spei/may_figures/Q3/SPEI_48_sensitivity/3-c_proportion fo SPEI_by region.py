# -*- coding: utf-8 -*-
"""
PANEL C: SPEI-48 Distribution (Violin + Density Bar Plot) - CORRECTED
Uses ACTUAL MONTHLY SPEI-48 observations from monthly_data_after_outlier_removal.csv
Dry = SPEI < -0.5, Wet = SPEI > 0.5, Neutral = -0.5 to 0.5
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Polygon
from scipy.stats import gaussian_kde
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

# ============================================
# LOAD MONTHLY SPEI-48 DATA (NOT site summaries!)
# ============================================

output_dir = Path(r"M:\Research\WUE_CUE\data_products\results\senstivity_april")

# Load the MONTHLY data file (has actual time series)
monthly_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
df_monthly = pd.read_csv(monthly_file)

# Load site metadata (coastal region assignments, etc.)
pickle_file = output_dir / 'processed_data_SPEI48_filtered.pkl'
with open(pickle_file, 'rb') as f:
    data = pickle.load(f)

# Get the list of sites that passed filtering
valid_sites = data['df']['site_name'].unique()
coast_order_ascending = data['coast_order_ascending']
coast_short_labels = data['coast_short_labels']
coast_full_names = data['coast_full_names']

# Filter monthly data to only valid sites
df_monthly = df_monthly[df_monthly['site_name'].isin(valid_sites)].copy()

# Apply Trans_ratio > 0 filter (same as upstream)
df_monthly = df_monthly[df_monthly['Trans_ratio'] > 0].copy()

# Get coastal region mapping from processed data
site_to_coast = dict(zip(data['df']['site_name'], data['df']['coast_region_analysis']))
df_monthly['coast_region_analysis'] = df_monthly['site_name'].map(site_to_coast)

# Drop rows without coast assignment
df_monthly = df_monthly.dropna(subset=['coast_region_analysis'])

print("="*70)
print("PANEL C: SPEI-48 Distribution (MONTHLY Observations)")
print("="*70)
print(f"\n📊 Monthly data shape: {df_monthly.shape}")
print(f"📊 Unique sites: {df_monthly['site_name'].nunique()}")
print(f"📊 Total monthly observations: {len(df_monthly):,}")
print(f"📊 SPEI-48 column: SPEI_48")
print(f"📊 Date range: {df_monthly['year'].min()}-{df_monthly['month'].min()} to {df_monthly['year'].max()}-{df_monthly['month'].max()}")

# ============================================
# BALANCED, ZERO-CENTERED COLORMAP
# ============================================

custom_cmap = LinearSegmentedColormap.from_list(
    "balanced_spei_gradient",
    ["#9E0142", "#D53E4F", "#F46D43", "#F7F7F7", "#ABD9E9", "#74ADD1", "#2B83BA"]
)

# ============================================
# SPEI-48 NORMALIZATION (ZERO-CENTERED)
# ============================================

global_spei_min = df_monthly['SPEI_48'].min()
global_spei_max = df_monthly['SPEI_48'].max()
spei_absmax = max(abs(global_spei_min), abs(global_spei_max))

norm_spei = TwoSlopeNorm(vmin=-spei_absmax, vcenter=0, vmax=spei_absmax)

print(f"\n  Global SPEI-48 range: {global_spei_min:.2f} to {global_spei_max:.2f}")
print(f"  Symmetric range: -{spei_absmax:.2f} to +{spei_absmax:.2f}")

# ============================================
# PREPARE DATA WITH MONTHLY SPEI VALUES
# ============================================

violin_data = []
region_stats = []

for coast in coast_order_ascending:
    # Get monthly SPEI values for this coast
    coast_df = df_monthly[df_monthly['coast_region_analysis'] == coast].copy()
    spei_values = coast_df['SPEI_48'].dropna().values
    
    if len(spei_values) > 0:
        violin_data.append(spei_values)
        
        # Count sites and months
        n_sites = coast_df['site_name'].nunique()
        n_months = len(spei_values)
        
        # Classify SPEI values
        dry_mask = spei_values < -0.5
        wet_mask = spei_values > 0.5
        neutral_mask = (spei_values >= -0.5) & (spei_values <= 0.5)
        
        pct_dry = (np.sum(dry_mask) / n_months) * 100
        pct_wet = (np.sum(wet_mask) / n_months) * 100
        pct_neutral = (np.sum(neutral_mask) / n_months) * 100
        
        region_stats.append({
            'short': coast_short_labels.get(coast, coast),
            'name': coast_full_names.get(coast, coast),
            'n_sites': n_sites,
            'n_months': n_months,
            'pct_dry': pct_dry,
            'pct_wet': pct_wet,
            'pct_neutral': pct_neutral,
            'values': spei_values
        })
        
        print(f"\n  {coast_full_names.get(coast, coast)} ({coast_short_labels.get(coast, coast)}):")
        print(f"    Sites (n) = {n_sites}")
        print(f"    Monthly observations (N) = {n_months:,}")
        print(f"    Dry (SPEI < -0.5): {pct_dry:.1f}%")
        print(f"    Neutral (-0.5 to 0.5): {pct_neutral:.1f}%")
        print(f"    Wet (SPEI > 0.5): {pct_wet:.1f}%")
        
    else:
        violin_data.append(np.array([0]))
        region_stats.append({
            'short': coast_short_labels.get(coast, coast),
            'n_sites': 0,
            'n_months': 0,
            'pct_dry': 0,
            'pct_wet': 0,
            'pct_neutral': 0,
            'values': np.array([0])
        })

n_regions = len(violin_data)
positions = np.arange(n_regions)

print("\n📊 Right panel: Density bar plot shows SPEI bins (0.5 unit intervals)")
print("   Bar color: Red (dry) → White (neutral) → Blue (wet)")
print("   Bar height: Proportion of MONTHLY observations in each SPEI bin")

# ============================================
# CREATE FIGURE WITH TWO SUBPLOTS
# ============================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 12), dpi=300)
fig.subplots_adjust(wspace=1)

# ============================================
# SUBPLOT 1: VIOLIN PLOT
# ============================================

print("\n📊 Creating violin plot (monthly observations)...")

def gradient_filled_violin(ax, data, positions, cmap, norm, width=0.8):
    for i, (pos, values) in enumerate(zip(positions, data)):
        if len(values) < 2:
            continue
            
        try:
            # Remove any NaN or inf values
            values = values[np.isfinite(values)]
            if len(values) < 2:
                continue
                
            kde = gaussian_kde(values)
            y_vals = np.linspace(values.min(), values.max(), 300)
            density = kde(y_vals)
            
            max_density = density.max()
            if max_density > 0:
                x_left = pos - (density / max_density) * (width / 2)
                x_right = pos + (density / max_density) * (width / 2)
            
            n_segments = 80
            y_min = values.min()
            y_max = values.max()
            
            for seg in range(n_segments):
                y_start = y_min + (seg / n_segments) * (y_max - y_min)
                y_end = y_min + ((seg + 1) / n_segments) * (y_max - y_min)
                y_mid = (y_start + y_end) / 2
                
                color = cmap(norm(y_mid))
                
                mask = (y_vals >= y_start) & (y_vals <= y_end)
                if not mask.any():
                    continue
                
                x_left_at_y = x_left[mask]
                x_right_at_y = x_right[mask]
                y_at_seg = y_vals[mask]
                
                if len(x_left_at_y) > 0 and len(x_right_at_y) > 0:
                    verts = []
                    for xl, y in zip(x_left_at_y, y_at_seg):
                        verts.append([xl, y])
                    for xr, y in zip(reversed(x_right_at_y), reversed(y_at_seg)):
                        verts.append([xr, y])
                    
                    if len(verts) >= 3:
                        polygon = Polygon(verts, closed=True, facecolor=color, 
                                         edgecolor='none', alpha=0.85)
                        ax.add_patch(polygon)
            
            median_val = np.median(values)
            ax.plot([pos - width/3, pos + width/3], [median_val, median_val], 
                   'k-', linewidth=2.5, alpha=0.9, zorder=10)
            
            ax.plot(x_left, y_vals, 'k-', linewidth=0.8, alpha=0.3, zorder=5)
            ax.plot(x_right, y_vals, 'k-', linewidth=0.8, alpha=0.3, zorder=5)
            
        except Exception as e:
            continue

gradient_filled_violin(ax1, violin_data, positions, custom_cmap, norm_spei, width=0.8)

# Customize violin plot
ax1.set_xticks(positions)
ax1.set_xticklabels([region_stats[i]['short'] for i in range(n_regions)], 
                    fontsize=32, fontweight='bold')
ax1.set_ylabel('SPEI-48', fontsize=34, fontweight='bold', labelpad=20)
ax1.set_ylim(-spei_absmax - 0.5, spei_absmax + 0.5)
ax1.axhline(y=0, color='#666666', linestyle='--', linewidth=2.5, alpha=0.7)
ax1.axhline(y=-0.5, color='#9E0142', linestyle=':', linewidth=1.5, alpha=0.5)
ax1.axhline(y=0.5, color='#2B83BA', linestyle=':', linewidth=1.5, alpha=0.5)
ax1.tick_params(axis='both', labelsize=32, width=3, length=12)
ax1.grid(True, alpha=0.15, axis='y', linestyle='--')

# ============================================
# ADD TWO N LABELS: small n (sites) and big N (monthly observations)
# ============================================

for idx, stats in enumerate(region_stats):
    if stats['n_months'] > 0:
        # Big N (monthly observations) - placed higher (closer to top of violin)
        y_pos_bigN = spei_absmax + 0.05
        ax1.text(idx, y_pos_bigN, f'n = {stats["n_months"]:,}', 
                ha='center', va='bottom', fontsize=24, fontweight='bold', 
                color='black', alpha=0.8, zorder=20)
        
        # Small n (sites) - placed below big N
        y_pos_smalln = spei_absmax + 0.05 - (spei_absmax * 0.08)
        ax1.text(idx, y_pos_smalln, f'N = {stats["n_sites"]}', 
                ha='center', va='bottom', fontsize=20, fontweight='normal',
                color='black', alpha=0.7, zorder=20)

# Dry/wet percentage labels
for idx, stats in enumerate(region_stats):
    if stats['n_months'] > 0:
        y_pos_dry = -spei_absmax - 1.3
        y_pos_wet = -spei_absmax - 1.7
        
        ax1.text(idx, y_pos_dry, f'{stats["pct_dry"]:.0f}% dry', 
                ha='center', va='top', fontsize=28, color='#9E0142', fontweight='bold')
        ax1.text(idx, y_pos_wet, f'{stats["pct_wet"]:.0f}% wet', 
                ha='center', va='top', fontsize=28, color='#2B83BA', fontweight='bold')

# Add panel label (a)
ax1.text(-0.12, 1.02, '(a)', transform=ax1.transAxes, 
         fontsize=48, fontweight='bold', va='bottom', ha='left')

# ============================================
# SUBPLOT 2: DENSITY BAR PLOT
# ============================================

print("📊 Creating density bar plot (monthly observations)...")

bins = np.arange(-3.5, 3.6, 0.5)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_colors = [custom_cmap(norm_spei(center)) for center in bin_centers]

# Find global max density
max_density_all = 0
for stats in region_stats:
    if stats['n_months'] > 0:
        hist_counts, _ = np.histogram(stats['values'], bins=bins)
        density = hist_counts / stats['n_months']
        max_density_all = max(max_density_all, density.max())

if max_density_all == 0:
    max_density_all = 0.5

# Create density bar plot
for idx, stats in enumerate(region_stats):
    if stats['n_months'] == 0:
        continue
    
    spei_values = stats['values']
    x_pos = idx + 0.5
    
    hist_counts, _ = np.histogram(spei_values, bins=bins)
    density = hist_counts / stats['n_months']
    
    for bin_idx, dens in enumerate(density):
        if dens > 0:
            bar_x = x_pos + (bin_centers[bin_idx] * 0.06)
            ax2.bar(bar_x, dens, width=0.06, color=bin_colors[bin_idx], 
                   edgecolor='white', linewidth=0.5, alpha=0.9)
    
    # Add density curve
    try:
        kde = gaussian_kde(spei_values)
        y_vals = np.linspace(-spei_absmax, spei_absmax, 200)
        density_curve = kde(y_vals)
        density_curve_scaled = density_curve / density_curve.max() * max_density_all * 0.8
        ax2.plot(x_pos + (y_vals * 0.06), density_curve_scaled, 
                'k-', linewidth=1.5, alpha=0.5, zorder=10)
    except:
        pass

# Customize density bar plot
ax2.set_xticks(positions + 0.5)
ax2.set_xticklabels([region_stats[i]['short'] for i in range(n_regions)], 
                    fontsize=32, fontweight='bold')
ax2.set_ylabel('Density (proportion)', fontsize=34, fontweight='bold', labelpad=20)
ax2.set_ylim(0, max_density_all * 1.15)
ax2.tick_params(axis='both', labelsize=32, width=3, length=12)
ax2.grid(True, alpha=0.15, axis='y', linestyle='--')
ax2.set_axisbelow(True)

# Add panel label (b)
ax2.text(-0.12, 1.02, '(b)', transform=ax2.transAxes, 
         fontsize=48, fontweight='bold', va='bottom', ha='left')

# ============================================
# ADD COLORBAR
# ============================================

sm = ScalarMappable(cmap=custom_cmap, norm=norm_spei)
sm.set_array([])

cbar_ax = fig.add_axes([0.3, 0.02, 0.4, 0.02])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cbar.set_label('SPEI-48 (drier ← normal → wetter)', 
               fontsize=28, fontweight='bold', labelpad=12)
cbar.ax.tick_params(labelsize=24)

max_tick = min(3, np.ceil(spei_absmax))
ticks = np.arange(-max_tick, max_tick + 0.1, 1)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f'{t:.0f}' for t in ticks])

# ============================================
# SAVE FIGURE
# ============================================

plt.tight_layout()
plt.subplots_adjust(bottom=0.22, top=0.95, left=0.08, right=0.95)

png_file = output_dir / 'PanelC_SPEI48_violin_density_bars.png'
pdf_file = output_dir / 'PanelC_SPEI48_violin_density_bars.pdf'

plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')

print(f"\n✓ Saved Panel C: {png_file}")
print(f"✓ Saved Panel C: {pdf_file}")

# ============================================
# FINAL SUMMARY
# ============================================
print("\n" + "="*70)
print("FINAL RESULTS (Monthly SPEI-48 Observations)")
print("Threshold: SPEI < -0.5 = dry, -0.5 to 0.5 = neutral, > 0.5 = wet")
print("="*70)
for stats in region_stats:
    if stats['n_months'] > 0:
        print(f"\n  {stats['name']} ({stats['short']}):")
        print(f"    Sites (n) = {stats['n_sites']}")
        print(f"    Monthly observations (N) = {stats['n_months']:,}")
        print(f"    Dry (SPEI < -0.5): {stats['pct_dry']:.1f}%")
        print(f"    Neutral (-0.5 to 0.5): {stats['pct_neutral']:.1f}%")
        print(f"    Wet (SPEI > 0.5): {stats['pct_wet']:.1f}%")

plt.show()

print("\n" + "="*70)
print("PANEL C COMPLETE - Using MONTHLY SPEI-48 observations")
print("="*70)