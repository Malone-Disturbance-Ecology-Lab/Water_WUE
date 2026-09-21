# -*- coding: utf-8 -*-
"""
THREE-PANEL T/ET PLOT (ENHANCED VISIBILITY)
- Borders darker and thicker
- Ticks longer and thicker
- Boxplot edges dark
- Uses 64 sites from monthly data (no SPEI filter)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

# ============================================
# PATHS
# ============================================
MONTHLY_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FILE = OUTPUT_DIR / "Figure_T_ET_three_panel_64sites.png"

print("="*60)
print("LOADING DATA (64 sites from monthly file)")
print("="*60)

# ============================================
# STEP 1: Load and filter monthly data (exact Study_area filter)
# ============================================
df = pd.read_csv(MONTHLY_FILE)
print(f"Loaded {len(df)} rows")

for col in df.select_dtypes(include="object").columns:
    df[col] = df[col].where(df[col].isna(), df[col].astype(str).str.strip())
    df[col] = df[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

required = ['site_name', 'Year', 'month', 'water_class', 'lat', 'long', 'Trans_ratio', 'WUE_tra']
for c in required + ['IGBP', 'climate', 'lai']:
    if c not in df.columns:
        print(f"Warning: Column '{c}' not found; will handle missing.")

df = df.dropna(subset=required).copy()
df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])].copy()
df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) &
        np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) & (df['Trans_ratio'] <= 1) &
        np.isfinite(df['WUE_tra'])].copy()

print(f"After filtering: {len(df)} rows, {df['site_name'].nunique()} unique sites")

site_medians = df.groupby(['site_name', 'IGBP', 'climate']).agg({
    'Trans_ratio': 'median',
    'lai': 'median'
}).reset_index()
site_medians.rename(columns={'Trans_ratio': 'T_ET'}, inplace=True)

print(f"Sites with T/ET: {len(site_medians)}")
site_medians_lai = site_medians.dropna(subset=['lai']).copy()
print(f"Sites with LAI: {len(site_medians_lai)}")

# ============================================
# STEP 2: Standardize climate names (remove 'with', use commas)
# ============================================
def standardize_climate(name):
    if pd.isna(name):
        return "Unknown"
    name = str(name)
    name = name.replace(' with ', ', ')
    name = name.replace(' with', ',')
    name = name.replace(')', '').strip()
    name = name.replace(', ,', ',')
    return name

site_medians['climate'] = site_medians['climate'].apply(standardize_climate)

veg_lai = site_medians.groupby('IGBP')['lai'].median().sort_values()
veg_order = veg_lai.sort_values(ascending=True).index.tolist()
climate_order = site_medians.groupby('climate')['T_ET'].median().sort_values(ascending=False).index.tolist()

print(f"Vegetation types: {len(veg_order)}")
print(f"Climate types: {len(climate_order)}")

# ============================================
# STEP 3: Wrap labels for climate (max_chars=16)
# ============================================
def wrap_label(text, max_chars=16):
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
# STEP 4: Colors and markers
# ============================================
very_different_colors = {
    'ENF': '#E41A1C',  'WET': '#377EB8',  'GRA': '#4DAF4A',
    'CSH': '#984EA3',  'BSV': '#FF7F00',  'DBF': '#FFFF33',
    'CRO': '#A65628',  'MF': '#F781BF',   'OSH': '#999999'
}
markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'H', 'X']

# ============================================
# STEP 5: Create three-panel plot with enhanced visibility
# ============================================
print("\nCreating three-panel plot (enhanced borders and ticks)...")
fig = plt.figure(figsize=(26, 28))
fig.subplots_adjust(hspace=0.35, left=0.08, right=0.95)

# --- Helper function to style each subplot ---
def style_axes(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
        spine.set_color('black')
    ax.tick_params(axis='both', which='major', length=10, width=2.5, colors='black')
    ax.tick_params(axis='both', which='minor', length=6, width=1.5, colors='black')
    ax.grid(True, alpha=0.3, axis='y', linewidth=1.2)

# --- Panel (a) ---
ax1 = fig.add_subplot(3, 1, 1)
box_data_veg = [site_medians[site_medians['IGBP'] == v]['T_ET'].values for v in veg_order]
n_veg = len(veg_order)
green_colors = plt.cm.Greens(np.linspace(0.4, 0.9, n_veg))
bp1 = ax1.boxplot(box_data_veg, labels=veg_order, patch_artist=True, 
                  widths=0.5, showfliers=False)
for patch, color in zip(bp1['boxes'], green_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.8)
means1 = [site_medians[site_medians['IGBP'] == v]['T_ET'].mean() for v in veg_order]
ax1.scatter(range(1, len(veg_order)+1), means1, color='red', marker='D', s=100, zorder=3, edgecolor='black')
ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=2)
ax1.set_ylabel('T:ET', fontsize=32, fontweight='bold')
ax1.set_title('(a) T:ET Distributions by Vegetation (IGBP)', fontsize=40, fontweight='bold', pad=30)
ax1.set_ylim(-0.05, 1.05)
ax1.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax1.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')
plt.setp(ax1.get_xticklabels(), rotation=45, ha='right', fontsize=32, fontweight='bold')
for i, veg in enumerate(veg_order, 1):
    n = len(site_medians[site_medians['IGBP'] == veg])
    ax1.text(i, 0.96, f'N={n}', ha='center', va='top', fontsize=28,
             fontweight='bold', transform=ax1.get_xaxis_transform())
style_axes(ax1)

# --- Panel (b) ---
ax2 = fig.add_subplot(3, 1, 2)
if len(site_medians_lai) >= 5:
    x = site_medians_lai['lai'].values
    y = site_medians_lai['T_ET'].values
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r2 = r_value**2
    p_text = 'p < 0.001' if p_value < 0.0001 else f'p = {p_value:.3f}'
    biomes = site_medians_lai['IGBP'].unique()
    marker_map = {biome: markers[i % len(markers)] for i, biome in enumerate(biomes)}
    for biome in biomes:
        subset = site_medians_lai[site_medians_lai['IGBP'] == biome]
        color = very_different_colors.get(biome, '#333333')
        marker = marker_map.get(biome, 'o')
        ax2.scatter(subset['lai'], subset['T_ET'], label=biome, alpha=0.8, s=140,
                    facecolor=color, edgecolor='black', linewidth=1.5, marker=marker)
    x_line = np.array([x.min(), x.max()])
    y_line = slope * x_line + intercept
    ax2.plot(x_line, y_line, 'k-', linewidth=3)
    ax2.set_xlabel('Leaf Area Index (LAI)', fontsize=32, fontweight='bold')
    ax2.set_ylabel('T:ET', fontsize=32, fontweight='bold')
    ax2.set_title('(b) T:ET vs LAI Relationship', fontsize=40, fontweight='bold', pad=30)
    ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=2)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='lower right', fontsize=26, ncol=2, frameon=True, fancybox=True, edgecolor='black')
    ax2.text(0.75, 0.05, f'R² = {r2:.2f}, {p_text}', 
             transform=ax2.transAxes, fontsize=26, fontweight='bold',
             ha='right', va='bottom',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')
    ax2.tick_params(axis='x', labelsize=32)
    ax2.text(0.02, 0.97, f'N = {len(site_medians_lai)} sites', transform=ax2.transAxes,
            fontsize=28, fontweight='bold', verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
style_axes(ax2)

# --- Panel (c) ---
ax3 = fig.add_subplot(3, 1, 3)
box_data_climate = [site_medians[site_medians['climate'] == c]['T_ET'].values for c in climate_order]
wrapped_labels = [wrap_label(c, max_chars=16) for c in climate_order]
bp3 = ax3.boxplot(box_data_climate, labels=wrapped_labels, patch_artist=True, 
                  widths=0.5, showfliers=False)
for patch in bp3['boxes']:
    patch.set_facecolor('#ADD8E6')
    patch.set_alpha(0.8)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.8)
means3 = [site_medians[site_medians['climate'] == c]['T_ET'].mean() for c in climate_order]
ax3.scatter(range(1, len(climate_order)+1), means3, color='red', marker='D', s=100, zorder=3, edgecolor='black')
ax3.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7, linewidth=2)
ax3.set_xlabel('Climate Zone', fontsize=32, fontweight='bold')
ax3.set_ylabel('T:ET', fontsize=32, fontweight='bold')
ax3.set_title('(c) T:ET Distributions by Climate', fontsize=40, fontweight='bold', pad=30)
ax3.set_ylim(-0.05, 1.05)
ax3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax3.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=32, fontweight='bold')
ax3.tick_params(axis='x', which='major', labelsize=20)
plt.setp(ax3.get_xticklabels(), rotation=30, ha='right', fontweight='bold')
for i, climate in enumerate(climate_order, 1):
    n = len(site_medians[site_medians['climate'] == climate])
    ax3.text(i, 0.96, f'N={n}', ha='center', va='top', fontsize=28,
             fontweight='bold', transform=ax3.get_xaxis_transform())
style_axes(ax3)

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

print("\n--- Climate ---")
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