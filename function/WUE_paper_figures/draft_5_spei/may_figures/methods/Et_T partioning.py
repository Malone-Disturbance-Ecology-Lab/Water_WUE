# -*- coding: utf-8 -*-
"""
THREE-PANEL T/ET PLOT using SAME sites as Figure 1 (strict triple intersection)
(a) T/ET distributions by Vegetation (IGBP) - boxplot (green shades by LAI)
(b) T/ET vs LAI - scatter plot (VERY DIFFERENT colors per IGBP)
(c) T/ET distributions by Climate - boxplot (bluish colors)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# ============================================
# PATHS
# ============================================
RESULTS_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
NN_MEDIANS_FILE = Path(RESULTS_DIR) / "wue_site_level_NN_medians_SPEI_1.csv"
MONTHLY_FILE = Path(RESULTS_DIR) / "monthly_data_after_outlier_removal.csv"
OUTPUT_DIR = Path(RESULTS_DIR) / "figures"
OUTPUT_FILE = OUTPUT_DIR / "Figure_T_ET_three_panel.png"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*60)
print("LOADING DATA...")
print("="*60)

# ============================================
# STEP 1: Get shared subset sites
# ============================================
nn_data = pd.read_csv(NN_MEDIANS_FILE)
wue_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE']['site_name'].dropna())
eva_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_eva']['site_name'].dropna())
tra_sites = set(nn_data[nn_data['WUE_Metric'] == 'WUE_tra']['site_name'].dropna())
shared_subset_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)
print(f"Shared subset sites: {len(shared_subset_sites)}")

# ============================================
# STEP 2: Load monthly data and STANDARDIZE climate names
# ============================================
df = pd.read_csv(MONTHLY_FILE)
df = df[df['site_name'].isin(shared_subset_sites)].copy()
df['T_ET'] = df['Trans_ratio']
df = df[(df['T_ET'] > 0) & (df['T_ET'].notna())].copy()

# ============================================
# STEP 3: STANDARDIZE climate names (remove 'with', use commas, KEEP FULL NAMES)
# ============================================
print("\n" + "="*60)
print("STANDARDIZING CLIMATE NAMES (remove 'with', use commas)")
print("="*60)

def standardize_climate(name):
    """Standardize climate names - replace 'with' with commas, keep full name"""
    if pd.isna(name):
        return "Unknown"
    
    name = str(name)
    
    # Replace ' with ' with comma and space
    name = name.replace(' with ', ', ')
    name = name.replace(' with', ',')
    
    # Remove trailing parentheses
    name = name.replace(')', '').strip()
    
    # Fix any double commas
    name = name.replace(', ,', ',')
    
    return name

# Apply standardization
df['climate_std'] = df['climate'].apply(standardize_climate)

# ============================================
# STEP 4: Aggregate to site-level with standardized climate
# ============================================
site_medians = df.groupby(['site_name', 'IGBP', 'climate_std']).agg({
    'T_ET': 'median',
    'lai': 'median'
}).reset_index()
site_medians.rename(columns={'climate_std': 'climate'}, inplace=True)

site_medians_lai = site_medians.dropna(subset=['lai']).copy()
print(f"\nSites: {len(site_medians)} | With LAI: {len(site_medians_lai)}")

# Calculate median LAI per vegetation
veg_lai = site_medians.groupby('IGBP')['lai'].median().sort_values()

# Print final unique climates
print(f"\nFinal unique climate classes: {len(site_medians['climate'].unique())}")

# ============================================
# STEP 5: Function to wrap labels (4+ rows for long names)
# ============================================
def wrap_label(text, max_chars=18):
    """Wrap text to multiple lines - can be 4 or more rows for long names"""
    words = text.split()
    lines = []
    current_line = []
    current_length = 0
    
    for word in words:
        if current_length + len(word) + 1 <= max_chars:
            current_line.append(word)
            current_length += len(word) + 1
        else:
            if current_line:
                lines.append(' '.join(current_line))
            current_line = [word]
            current_length = len(word) + 1
    
    if current_line:
        lines.append(' '.join(current_line))
    
    return '\n'.join(lines)

# ============================================
# STEP 6: VERY DIFFERENT COLORS for panel b (no blue shades)
# ============================================
very_different_colors = {
    'ENF': '#E41A1C',  # Bright Red
    'WET': '#377EB8',  # Blue
    'GRA': '#4DAF4A',  # Green
    'CSH': '#984EA3',  # Purple
    'BSV': '#FF7F00',  # Orange
    'DBF': '#FFFF33',  # Yellow
    'CRO': '#A65628',  # Brown
    'MF': '#F781BF',   # Pink
    'OSH': '#999999'   # Gray
}

# Different markers for each IGBP
markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'H', 'X']

# ============================================
# STEP 7: Create three-panel plot (INCREASED HEIGHT)
# ============================================
print("\nCreating three-panel plot...")

# INCREASED figure height from 20 to 28
fig = plt.figure(figsize=(24, 28))
fig.subplots_adjust(hspace=0.35, left=0.08, right=0.95)

# --------------------------------------------
# PANEL (a): T/ET by Vegetation
# --------------------------------------------
ax1 = fig.add_subplot(3, 1, 1)

veg_order = veg_lai.sort_values(ascending=True).index.tolist()
box_data_veg = [site_medians[site_medians['IGBP'] == v]['T_ET'].values for v in veg_order]

n_veg = len(veg_order)
green_colors = plt.cm.Greens(np.linspace(0.4, 0.9, n_veg))

bp1 = ax1.boxplot(box_data_veg, labels=veg_order, patch_artist=True, 
                  widths=0.5, showfliers=False)

for patch, color in zip(bp1['boxes'], green_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.5)

means1 = [site_medians[site_medians['IGBP'] == v]['T_ET'].mean() for v in veg_order]
ax1.scatter(range(1, len(veg_order)+1), means1, color='red', marker='D', s=80, zorder=3)
ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)

ax1.set_ylabel('T/ET', fontsize=32, fontweight='bold')
ax1.set_title('(a) T/ET Distributions by Vegetation (IGBP)', fontsize=40, fontweight='bold', pad=30)

ax1.grid(True, alpha=0.3, axis='y')
ax1.set_ylim(-0.05, 1.05)

# INCREASE y-axis ticks to 5-6
ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax1.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')

plt.setp(ax1.get_xticklabels(), rotation=45, ha='right', fontsize=32, fontweight='bold')

for i, veg in enumerate(veg_order, 1):
    n = len(site_medians[site_medians['IGBP'] == veg])
    ax1.text(i, 0.96, f'N={n}', ha='center', va='top', fontsize=28,
             fontweight='bold', transform=ax1.get_xaxis_transform())

# --------------------------------------------
# PANEL (b): T/ET vs LAI (VERY DIFFERENT COLORS)
# --------------------------------------------
ax2 = fig.add_subplot(3, 1, 2)

if len(site_medians_lai) >= 5:
    x = site_medians_lai['lai'].values
    y = site_medians_lai['T_ET'].values
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r2 = r_value**2
    
    # Format p-value for display (p < 0.0001)
    if p_value < 0.0001:
        p_text = 'p < 0.001'
    else:
        p_text = f'p = {p_value:.3f}'
    
    biomes = site_medians_lai['IGBP'].unique()
    marker_map = {biome: markers[i % len(markers)] for i, biome in enumerate(biomes)}
    
    for biome in biomes:
        subset = site_medians_lai[site_medians_lai['IGBP'] == biome]
        color = very_different_colors.get(biome, '#333333')
        marker = marker_map.get(biome, 'o')
        
        ax2.scatter(subset['lai'], subset['T_ET'], 
                   label=biome, alpha=0.8, s=120,
                   facecolor=color, edgecolor=color, linewidth=0.5,
                   marker=marker)
    
    x_line = np.array([x.min(), x.max()])
    y_line = slope * x_line + intercept
    ax2.plot(x_line, y_line, 'k-', linewidth=2.5)
    
    ax2.set_xlabel('Leaf Area Index (LAI)', fontsize=32, fontweight='bold')
    ax2.set_ylabel('T/ET', fontsize=32, fontweight='bold')
    ax2.set_title('(b) T/ET vs LAI Relationship', fontsize=40, fontweight='bold', pad=30)
    
    ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
    ax2.grid(True, alpha=0.3)
    
    # TWO-COLUMN LEGEND (ncol=2)
    legend = ax2.legend(loc='lower right', fontsize=26, ncol=2, frameon=True, fancybox=True,
                        edgecolor='black')
    
    # Add R² and p-value as text above legend or in title
    ax2.text(0.75, 0.05, f'R² = {r2:.2f}, {p_text}', 
             transform=ax2.transAxes, fontsize=26, fontweight='bold',
             ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    
    ax2.set_ylim(-0.05, 1.05)
    
    # INCREASE y-axis ticks to 5-6
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')
    ax2.tick_params(axis='x', labelsize=32)
    
    # N label
    ax2.text(0.02, 0.97, f'N = {len(site_medians_lai)} sites', transform=ax2.transAxes,
            fontsize=28, fontweight='bold', verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))

# --------------------------------------------
# PANEL (c): T/ET by Climate (FULL NAMES, 4+ rows)
# --------------------------------------------
ax3 = fig.add_subplot(3, 1, 3)

# Get unique climates (standardized, full names)
climate_order = site_medians.groupby('climate')['T_ET'].median().sort_values(ascending=False).index.tolist()

box_data_climate = [site_medians[site_medians['climate'] == c]['T_ET'].values for c in climate_order]

# Wrap labels (can be 4 or more rows for long names)
wrapped_labels = [wrap_label(c, max_chars=18) for c in climate_order]

bp3 = ax3.boxplot(box_data_climate, labels=wrapped_labels, patch_artist=True, 
                  widths=0.5, showfliers=False)

for patch in bp3['boxes']:
    patch.set_facecolor('#ADD8E6')
    patch.set_alpha(0.8)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.5)

means3 = [site_medians[site_medians['climate'] == c]['T_ET'].mean() for c in climate_order]
ax3.scatter(range(1, len(climate_order)+1), means3, color='red', marker='D', s=80, zorder=3)
ax3.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)

ax3.set_xlabel('Climate Zone', fontsize=32, fontweight='bold')
ax3.set_ylabel('T/ET', fontsize=32, fontweight='bold')
ax3.set_title('(c) T/ET Distributions by Climate', fontsize=40, fontweight='bold', pad=30)

ax3.grid(True, alpha=0.3, axis='y')
ax3.set_ylim(-0.05, 1.05)

# INCREASE y-axis ticks to 5-6 for panel c
ax3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax3.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')

# SLIGHTLY REDUCED x-tick font for panel c (from 24 to 22 to fit better)
ax3.tick_params(axis='x', which='major', labelsize=22)
plt.setp(ax3.get_xticklabels(), ha='center', va='top', fontweight='bold')

# N labels
for i, climate in enumerate(climate_order, 1):
    n = len(site_medians[site_medians['climate'] == climate])
    ax3.text(i, 0.96, f'N={n}', ha='center', va='top', fontsize=28,
             fontweight='bold', transform=ax3.get_xaxis_transform())

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=600, bbox_inches='tight')
print(f"\n✅ Plot saved: {OUTPUT_FILE} (600 DPI)")
plt.show()

# ============================================
# PRINT STATISTICS
# ============================================
print("\n" + "="*60)
print("STATISTICS")
print("="*60)

print("\n--- Vegetation (IGBP) ---")
for veg in veg_order:
    data = site_medians[site_medians['IGBP'] == veg]['T_ET']
    lai_val = site_medians[site_medians['IGBP'] == veg]['lai'].median()
    q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
    print(f"{veg} (LAI={lai_val:.2f}): N={len(data)}, Mean={data.mean():.3f}±{data.std():.3f}, Median={data.median():.3f}, IQR={q1:.3f}-{q3:.3f}")

print("\n--- Climate (FULL NAMES, commas instead of 'with') ---")
for climate in climate_order:
    data = site_medians[site_medians['climate'] == climate]['T_ET']
    q1, q3 = np.percentile(data, 25), np.percentile(data, 75)
    print(f"{climate}: N={len(data)}, Mean={data.mean():.3f}±{data.std():.3f}, Median={data.median():.3f}, IQR={q1:.3f}-{q3:.3f}")

if len(site_medians_lai) >= 5:
    print(f"\n--- LAI ---")
    print(f"N={len(site_medians_lai)}, Slope={slope:.4f}, R²={r2:.3f}, p={p_value:.2e}")

print("\n" + "="*60)
print("DONE!")
print("="*60)