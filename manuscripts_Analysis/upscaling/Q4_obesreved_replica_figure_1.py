#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 2, Panel B – Observed SPEI‑3 distribution by coastline.
Uses all available months (data are already growing‑season filtered per site).
"""

import os
import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import PercentFormatter

# ============================================================================
# 1. COASTLINE CLASSIFICATION FUNCTIONS (copied verbatim)
# ============================================================================
EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9),
        (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2),
        (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0),
        (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7),
        (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1),
        (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6),
        (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3),
        (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9),
        (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)
    ],
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
    lat0 = math.radians(lat)
    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y])
    p = np.array([0.0, 0.0])
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)
    ab = b - a
    denom = np.dot(ab, ab)
    if denom == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, np.dot(p - a, ab) / denom))
    return float(np.linalg.norm(p - (a + t * ab)))

def _distance_to_polyline(lat, lon, polyline):
    return min(
        _seg_dist(lat, lon, *polyline[i], *polyline[i + 1])
        for i in range(len(polyline) - 1)
    )

def assign_coast_region(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return np.nan
    lat = float(lat); lon = float(lon)
    if lat > 50:
        return "AK Coast"
    dists = {
        coast: _distance_to_polyline(lat, lon, polyline)
        for coast, polyline in COAST_POLYLINES.items()
    }
    return min(dists, key=dists.get)

# ============================================================================
# 2. DATA LOADING & FILTERING (no month restriction)
# ============================================================================
DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

print("\n" + "="*70)
print("STEP 1: Loading and filtering observational data")
print("="*70)

try:
    df_raw = pd.read_csv(DATA_PATH)
except FileNotFoundError:
    print(f"ERROR: File not found at:\n{DATA_PATH}")
    sys.exit(1)

# Clean strings
for col in df_raw.select_dtypes(include='object').columns:
    df_raw[col] = df_raw[col].map(lambda x: x.strip() if isinstance(x, str) else x)

# Required columns
req_cols = ['site_name', 'Year', 'month', 'water_class',
            'lat', 'long', 'Trans_ratio', 'WUE_tra']
df = df_raw.dropna(subset=req_cols)

# Retain only Upland, Freshwater, Saline
df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])]

# Finite lat/long, Trans_ratio in [0,1], finite WUE_tra
df = df[np.isfinite(df['lat']) & np.isfinite(df['long'])]
df = df[np.isfinite(df['Trans_ratio'])]
df = df[(df['Trans_ratio'] >= 0) & (df['Trans_ratio'] <= 1)]
df = df[np.isfinite(df['WUE_tra'])]

print(f"Retained rows: {len(df)}")
print(f"Unique sites: {df['site_name'].nunique()}")

if len(df) != 1852 or df['site_name'].nunique() != 64:
    print("ERROR: Expected 1,852 rows and 64 sites. Stopping.")
    sys.exit(1)
else:
    print("✓ Checkpoint passed: 1,852 rows, 64 sites.\n")

# ============================================================================
# 3. ASSIGN COASTLINES
# ============================================================================
print("="*70)
print("STEP 2: Assigning coastlines to sites")
print("="*70)

sites = df[['site_name', 'lat', 'long']].drop_duplicates()
sites['coast_raw'] = sites.apply(
    lambda row: assign_coast_region(row['lat'], row['long']), axis=1
)
sites['coastline'] = sites['coast_raw'].replace({'AK Coast': 'Alaska Coast'})

df = df.merge(sites[['site_name', 'coastline']], on='site_name', how='left')
if df['coastline'].isnull().any():
    print("ERROR: Some sites have no coastline.")
    sys.exit(1)

coast_counts = df.groupby('coastline')['site_name'].nunique().sort_index()
expected = {'Alaska Coast': 12, 'Pacific Coast': 23,
            'Gulf Coast': 7, 'Atlantic Coast': 22}

print("Unique sites by coastline:")
for c in ['Alaska Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
    print(f"  {c}: {coast_counts.get(c, 0)}")

if coast_counts.to_dict() != expected:
    print("ERROR: Coastline counts do not match expected.")
    print(f"Expected: {expected}")
    print(f"Observed: {coast_counts.to_dict()}")
    sys.exit(1)
else:
    print("✓ Coastline assignment verified.\n")

# ============================================================================
# 4. SELECT DATA FOR PANEL B – NO MONTH FILTER, NO YEAR FILTER
# ============================================================================
print("="*70)
print("STEP 3: Selecting data with valid SPEI‑3 (all months, all years)")
print("="*70)

# Use all available months – data are already growing‑season specific per site
df_b = df[np.isfinite(df['SPEI_3'])].copy()

# Check for exact duplicates (site-year-month)
dups = df_b.duplicated(subset=['site_name', 'Year', 'month'], keep=False)
if dups.any():
    print("WARNING: Exact repeated site-year-month records found:")
    print(df_b[dups][['site_name', 'Year', 'month']].drop_duplicates())
else:
    print("No exact repeated site-year-month records.")

total_obs = len(df_b)
total_sites = df_b['site_name'].nunique()
print(f"Total observed site‑months: {total_obs}")
print(f"Unique sites represented: {total_sites}")

obs_by_coast = df_b.groupby('coastline').size()
sites_by_coast = df_b.groupby('coastline')['site_name'].nunique()
print("\nObserved site‑months by coastline:")
print(obs_by_coast)
print("\nUnique sites by coastline:")
print(sites_by_coast)

if total_sites != 64:
    missing = set(df['site_name'].unique()) - set(df_b['site_name'].unique())
    print(f"WARNING: Missing sites: {missing}")
else:
    print("✓ All 64 sites retained.\n")

# ============================================================================
# 5. CLASSIFY SPEI‑3 INTO FIVE CLASSES
# ============================================================================
print("="*70)
print("STEP 4: Classifying SPEI‑3 values")
print("="*70)

def classify_spei(x):
    if x <= -1.5:
        return '≤ -1.5'
    elif x <= -1.0:
        return '-1.5 to -1.0'
    elif x < 1.0:
        return '-1.0 to 1.0'
    elif x < 1.5:
        return '1.0 to 1.5'
    else:
        return '≥ 1.5'

df_b['spei_class'] = df_b['SPEI_3'].apply(classify_spei)

class_order = ['≤ -1.5', '-1.5 to -1.0', '-1.0 to 1.0', '1.0 to 1.5', '≥ 1.5']
coastline_order = ['Alaska Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']

percentages = {}
print("Percentages per coastline (sum to 100%):")
for coast in coastline_order:
    subset = df_b[df_b['coastline'] == coast]
    total = len(subset)
    if total == 0:
        pct = {c: 0 for c in class_order}
    else:
        counts = subset['spei_class'].value_counts()
        pct = {c: (counts.get(c, 0) / total) * 100 for c in class_order}
    percentages[coast] = pct
    total_pct = sum(pct.values())
    print(f"  {coast}: {total_pct:.1f}% (OK)" if abs(total_pct-100)<0.01 else f"  {coast}: {total_pct:.1f}% (WARNING)")

print()

# ============================================================================
# 6. CREATE PANEL B – WITH CORRECTED LABEL COLOUR LOGIC
# ============================================================================
print("="*70)
print("STEP 5: Creating Panel B figure")
print("="*70)

# Coast colours
coast_colors = {
    'Alaska Coast':  '#009E73',
    'Pacific Coast': '#0072B2',
    'Gulf Coast':    '#B3B300',
    'Atlantic Coast':'#CC79A7'
}

# SPEI‑3 colours
spei_colors = {
    '≤ -1.5':      '#b2182b',   # dark red
    '-1.5 to -1.0': '#e08214',  # orange
    '-1.0 to 1.0':  '#d9d9d9',  # light gray
    '1.0 to 1.5':   '#67a9cf',  # light blue
    '≥ 1.5':        '#2166ac'   # dark blue
}

# Prepare data
plot_data = pd.DataFrame(percentages).T[class_order].loc[coastline_order]
n_sites = df_b.groupby('coastline')['site_name'].nunique().reindex(coastline_order)

# Wider figure
fig, ax = plt.subplots(figsize=(14, 5.5))
plt.subplots_adjust(left=0.22, bottom=0.25)
ax.set_xlim(0, 105)

bottom = np.zeros(len(plot_data))
stagger_count = np.zeros(len(plot_data), dtype=int)

for i, cls in enumerate(class_order):
    vals = plot_data[cls].values
    ax.barh(plot_data.index, vals, left=bottom,
            color=spei_colors[cls], edgecolor='white', linewidth=0.8)

    for j, val in enumerate(vals):
        if val <= 0:
            continue

        center_x = bottom[j] + val / 2
        label = f'{val:.1f}%'

        # White for dark red, orange, and dark blue
        # Black for gray and light blue
        if cls in ['≤ -1.5', '-1.5 to -1.0', '≥ 1.5']:
            text_color = 'white'
        else:
            text_color = 'black'

        if val >= 5:
            ax.text(
                center_x, j, label,
                ha='center', va='center',
                fontsize=11,
                fontweight='bold',
                color=text_color,
                clip_on=False,
                zorder=5
            )
        else:
            # Narrow segment – stagger vertically
            offset_y = 10 if stagger_count[j] % 2 == 0 else -10
            vertical_align = 'bottom' if stagger_count[j] % 2 == 0 else 'top'
            stagger_count[j] += 1

            ax.annotate(
                label,
                xy=(center_x, j),
                xytext=(0, offset_y),
                textcoords='offset points',
                ha='center',
                va=vertical_align,
                fontsize=9,
                fontweight='bold',
                color=text_color,          # use text_color (not fixed black)
                clip_on=False,
                zorder=6,
                arrowprops=dict(
                    arrowstyle='-',
                    linewidth=0.4,
                    color='gray',
                    shrinkA=2,
                    shrinkB=2
                )
            )

    bottom += vals

# Coast labels
for j, coast in enumerate(plot_data.index):
    label = f"{coast} (N = {int(n_sites[coast])})"
    ax.text(-2.5, j, label, ha='right', va='center',
            fontsize=12, fontweight='bold', color=coast_colors[coast])

# Clean spines
for spine in ['top', 'right', 'left']:
    ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_visible(True)

ax.set_xlabel('Observed growing‑season site‑months (%)', fontsize=13, fontweight='bold')
ax.xaxis.set_major_formatter(PercentFormatter())
ax.set_ylabel('')
ax.tick_params(axis='y', labelsize=0)

# Panel label
ax.text(-0.04, 1.02, 'b)', transform=ax.transAxes,
        fontsize=16, fontweight='bold', va='bottom', ha='left')

# Legend
legend_patches = [mpatches.Patch(color=spei_colors[cls], label=cls) for cls in class_order]
legend = ax.legend(handles=legend_patches, loc='lower center',
                   bbox_to_anchor=(0.5, -0.45), ncol=5,
                   fontsize=11, title='SPEI-3')
legend.get_title().set_fontsize(12)
legend.get_title().set_fontweight('bold')

plt.subplots_adjust(bottom=0.28)

# ============================================================================
# SAVE PNG ONLY TO THE SPECIFIED DIRECTORY
# ============================================================================
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling\aug_figure"
os.makedirs(output_dir, exist_ok=True)

fig.savefig(
    os.path.join(output_dir, 'Figure2_PanelB.png'),
    dpi=600,
    bbox_inches='tight',
    facecolor='white',
    transparent=False
)

print(f"Panel B saved to {os.path.join(output_dir, 'Figure2_PanelB.png')}\n")
plt.show()

# ============================================================================
# 7. PRINT SUMMARY STATISTICS
# ============================================================================
print("="*70)
print("STEP 6: Summary statistics for manuscript")
print("="*70)

summary_list = []
for coast in coastline_order:
    subset = df_b[df_b['coastline'] == coast]
    n_sites = subset['site_name'].nunique()
    n_obs = len(subset)
    if n_obs == 0:
        continue
    mean_sp = subset['SPEI_3'].mean()
    median_sp = subset['SPEI_3'].median()
    min_sp = subset['SPEI_3'].min()
    max_sp = subset['SPEI_3'].max()
    pct_le_minus1 = (subset['SPEI_3'] <= -1.0).mean() * 100
    pct_le_minus1_5 = (subset['SPEI_3'] <= -1.5).mean() * 100
    pct_le_minus2 = (subset['SPEI_3'] <= -2.0).mean() * 100
    pct_ge_1 = (subset['SPEI_3'] >= 1.0).mean() * 100
    pct_ge_1_5 = (subset['SPEI_3'] >= 1.5).mean() * 100
    sites_dry1 = subset[subset['SPEI_3'] <= -1.0]['site_name'].nunique()
    pct_sites_dry1 = (sites_dry1 / n_sites) * 100 if n_sites else 0
    sites_dry1_5 = subset[subset['SPEI_3'] <= -1.5]['site_name'].nunique()
    pct_sites_dry1_5 = (sites_dry1_5 / n_sites) * 100 if n_sites else 0
    pct_classes = percentages[coast]

    summary_dict = {
        'coastline': coast,
        'n_sites': n_sites,
        'n_obs': n_obs,
        'mean_SPEI3': mean_sp,
        'median_SPEI3': median_sp,
        'min_SPEI3': min_sp,
        'max_SPEI3': max_sp,
        'pct_<=-1.0': pct_le_minus1,
        'pct_<=-1.5': pct_le_minus1_5,
        'pct_<=-2.0': pct_le_minus2,
        'pct_>=1.0': pct_ge_1,
        'pct_>=1.5': pct_ge_1_5,
        'sites_dry1': sites_dry1,
        'pct_sites_dry1': pct_sites_dry1,
        'sites_dry1.5': sites_dry1_5,
        'pct_sites_dry1.5': pct_sites_dry1_5,
        **{f'pct_{cls}': pct_classes[cls] for cls in class_order}
    }
    summary_list.append(summary_dict)

summary_df = pd.DataFrame(summary_list)
print(summary_df.to_string(index=False))

summary_df.to_csv(os.path.join(output_dir, 'PanelB_summary_statistics.csv'), index=False)
print(f"\nSummary saved to {os.path.join(output_dir, 'PanelB_summary_statistics.csv')}")

print("\n" + "="*70)
print("SCRIPT FINISHED SUCCESSFULLY")
print("="*70)