#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 2 – Combined panels a & b.

Panel a: Distribution of observed SPEI‑3 classes by coastline (stacked bar).
Panel b: Site‑level drought‑event exposure (lollipop chart).

Both panels use the same observational dataset:
- 1,852 initial rows, 64 sites
- All years, all months (growing‑season already filtered per site)
- Finite SPEI‑3 values
"""

import os
import sys
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter

# ============================================================================
# 1. COASTLINE CLASSIFICATION FUNCTIONS (copied from Panel B)
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
# 2. CONSTANTS
# ============================================================================
COAST_COLORS = {
    "Alaska Coast":   "#009E73",
    "Pacific Coast":  "#0072B2",
    "Gulf Coast":     "#B3B300",
    "Atlantic Coast": "#CC79A7"
}
COAST_ORDER = ["Alaska Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]

# SPEI‑3 class colours
SPEI_COLORS = {
    '≤ -1.5':      '#b2182b',
    '-1.5 to -1.0': '#e08214',
    '-1.0 to 1.0':  '#d9d9d9',
    '1.0 to 1.5':   '#67a9cf',
    '≥ 1.5':        '#2166ac'
}
CLASS_ORDER = ['≤ -1.5', '-1.5 to -1.0', '-1.0 to 1.0', '1.0 to 1.5', '≥ 1.5']

# ============================================================================
# 3. DATA LOADING & FILTERING (identical for both panels)
# ============================================================================
DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

print("\n" + "="*70)
print("Figure 2 – Combined panels: loading and filtering data")
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

print(f"Retained rows after initial filters: {len(df)}")
print(f"Unique sites: {df['site_name'].nunique()}")

if len(df) != 1852 or df['site_name'].nunique() != 64:
    print("ERROR: Expected 1,852 rows and 64 sites. Stopping.")
    sys.exit(1)
else:
    print("✓ Checkpoint passed: 1,852 rows, 64 sites.\n")

# ============================================================================
# 4. ASSIGN COASTLINES
# ============================================================================
print("Assigning coastlines...")
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
for c in COAST_ORDER:
    print(f"  {c}: {coast_counts.get(c, 0)}")
if coast_counts.to_dict() != expected:
    print("ERROR: Coastline counts do not match expected.")
    sys.exit(1)
else:
    print("✓ Coastline assignment verified.\n")

# ============================================================================
# 5. COMMON ANALYSIS DATASET (no month/year filter)
# ============================================================================
print("Creating analysis dataset for both panels...")
df_common = df[np.isfinite(df['SPEI_3'])].copy()

# Check duplicates
dups = df_common.duplicated(subset=['site_name', 'Year', 'month'], keep=False)
if dups.any():
    print("ERROR: Exact repeated site-year-month records found.")
    sys.exit(1)
else:
    print("No exact repeated site-year-month records.")

total_obs = len(df_common)
total_sites = df_common['site_name'].nunique()
print(f"Total observed site‑months: {total_obs}")
print(f"Unique sites: {total_sites}")

if total_sites != 64:
    print(f"ERROR: Only {total_sites} sites represented.")
    sys.exit(1)
else:
    print("✓ All 64 sites retained.\n")

# ============================================================================
# 6. PANEL a DATA: SPEI‑3 CLASSES
# ============================================================================
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

df_common['spei_class'] = df_common['SPEI_3'].apply(classify_spei)

percentages = {}
for coast in COAST_ORDER:
    subset = df_common[df_common['coastline'] == coast]
    total = len(subset)
    if total == 0:
        pct = {c: 0 for c in CLASS_ORDER}
    else:
        counts = subset['spei_class'].value_counts()
        pct = {c: (counts.get(c, 0) / total) * 100 for c in CLASS_ORDER}
    percentages[coast] = pct

# ============================================================================
# 7. PANEL b DATA: DROUGHT EVENT EXPOSURE
# ============================================================================
df_common['date'] = pd.to_datetime(df_common[['Year', 'month']].assign(day=1))

site_results = []
for site_name, group in df_common.groupby('site_name'):
    group = group.sort_values('date').copy()
    n_total = len(group)
    spei = group['SPEI_3'].to_numpy()
    periods = group['date'].dt.to_period('M').to_numpy()
    dry = spei < -1.0

    runs = []
    i = 0
    while i < len(dry):
        if not dry[i]:
            i += 1
            continue
        start = i
        j = i + 1
        while j < len(dry):
            if dry[j] and (periods[j] == periods[j-1] + 1):
                j += 1
            else:
                break
        if j - start >= 2:
            runs.append((start, j))
        i = j

    event_months = np.zeros(len(dry), dtype=bool)
    for start, end in runs:
        event_months[start:end] = True
    n_event_months = int(event_months.sum())
    n_events = len(runs)
    event_durations = [end - start for start, end in runs]
    mean_duration = np.mean(event_durations) if event_durations else 0
    max_duration = max(event_durations) if event_durations else 0
    exposure_pct = (n_event_months / n_total) * 100 if n_total > 0 else 0

    site_results.append({
        'site_name': site_name,
        'coastline': group['coastline'].iloc[0],
        'n_valid_SPEI3_months': n_total,
        'n_drought_events': n_events,
        'n_drought_event_months': n_event_months,
        'mean_event_duration_months': mean_duration,
        'maximum_event_duration_months': max_duration,
        'drought_event_exposure_percent': exposure_pct
    })

site_df = pd.DataFrame(site_results)

# ============================================================================
# 8. CREATE COMBINED FIGURE (panels a & b) – FIXED OVERLAP
# ============================================================================
print("Creating combined Figure 2...")

fig = plt.figure(figsize=(12, 12))   # increased height for more space

# ----- Panel a (stacked bar) -----
ax_a = plt.subplot(2, 1, 1)
plot_data_a = pd.DataFrame(percentages).T[CLASS_ORDER].loc[COAST_ORDER]
n_sites_a = df_common.groupby('coastline')['site_name'].nunique().reindex(COAST_ORDER)

bottom = np.zeros(len(COAST_ORDER))
stagger_count = np.zeros(len(COAST_ORDER), dtype=int)

for i, cls in enumerate(CLASS_ORDER):
    vals = plot_data_a[cls].values
    ax_a.barh(COAST_ORDER, vals, left=bottom,
              color=SPEI_COLORS[cls], edgecolor='white', linewidth=0.8)
    for j, val in enumerate(vals):
        if val <= 0:
            continue
        center_x = bottom[j] + val / 2
        label = f'{val:.1f}%'
        if cls in ['≤ -1.5', '-1.5 to -1.0', '≥ 1.5']:
            txt_col = 'white'
        else:
            txt_col = 'black'
        if val >= 5:
            ax_a.text(center_x, j, label, ha='center', va='center',
                      fontsize=10, fontweight='bold', color=txt_col, clip_on=False)
        else:
            off_y = 10 if stagger_count[j] % 2 == 0 else -10
            va = 'bottom' if stagger_count[j] % 2 == 0 else 'top'
            stagger_count[j] += 1
            ax_a.annotate(label, xy=(center_x, j), xytext=(0, off_y),
                          textcoords='offset points', ha='center', va=va,
                          fontsize=8, fontweight='bold', color=txt_col,
                          arrowprops=dict(arrowstyle='-', linewidth=0.4, color='gray'),
                          clip_on=False)
    bottom += vals

# Formatting Panel a
ax_a.set_xlim(0, 105)
ax_a.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
ax_a.set_xlabel('Observed growing‑season site‑months (%)', fontsize=14, fontweight='bold')
ax_a.set_ylabel('')
ax_a.tick_params(axis='y', labelsize=14)
ax_a.tick_params(axis='x', labelsize=14)
ax_a.set_yticks(range(len(COAST_ORDER)))
ax_a.set_yticklabels(
    [f"{coast} (N = {int(n_sites_a[coast])})" for coast in COAST_ORDER],
    fontsize=14, fontweight='bold'
)
for tick, coast in zip(ax_a.get_yticklabels(), COAST_ORDER):
    tick.set_color(COAST_COLORS[coast])
ax_a.invert_yaxis()

ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.spines['left'].set_visible(False)
ax_a.spines['bottom'].set_visible(True)

# Legend – moved further down to avoid x‑axis overlap
legend_patches = [mpatches.Patch(color=SPEI_COLORS[cls], label=cls) for cls in CLASS_ORDER]
leg_a = ax_a.legend(handles=legend_patches, loc='lower center',
                    bbox_to_anchor=(0.5, -0.45), ncol=5,
                    fontsize=12, title='SPEI-3')
leg_a.get_title().set_fontsize(13)
leg_a.get_title().set_fontweight('bold')

# Panel label
ax_a.text(-0.02, 1.02, 'a)', transform=ax_a.transAxes,
          fontsize=18, fontweight='bold', va='bottom', ha='left')

# ----- Panel b (lollipop) -----
ax_b = plt.subplot(2, 1, 2)
coast_y = {coast: i for i, coast in enumerate(COAST_ORDER)}

for coast in COAST_ORDER:
    sub = site_df[site_df['coastline'] == coast].sort_values('drought_event_exposure_percent')
    base_y = coast_y[coast]
    x_vals = sub['drought_event_exposure_percent'].values
    unique_vals, counts = np.unique(x_vals, return_counts=True)
    y_offsets = np.zeros(len(sub))
    for uv, cnt in zip(unique_vals, counts):
        if cnt == 1:
            offsets = [0]
        else:
            offsets = np.linspace(-0.15, 0.15, cnt)
        idx = np.where(x_vals == uv)[0]
        y_offsets[idx] = offsets

    for x, yoff in zip(x_vals, y_offsets):
        y_pos = base_y + yoff
        ax_b.plot([x, x], [base_y, y_pos], color='gray', linewidth=0.8, alpha=0.6, zorder=1)
        ax_b.scatter(x, y_pos, s=80, color=COAST_COLORS[coast],
                     edgecolor='white', linewidth=0.5, zorder=3)

    ax_b.axhline(y=base_y, color='black', linewidth=0.5, alpha=0.3, zorder=0)
    mean_val = sub['drought_event_exposure_percent'].mean()
    ax_b.scatter(mean_val, base_y, s=120, marker='D', color='black', zorder=5)

# Formatting Panel b
ax_b.set_yticks(range(len(COAST_ORDER)))
ax_b.set_yticklabels(
    [f"{coast} (N = {len(site_df[site_df['coastline']==coast])})" for coast in COAST_ORDER],
    fontsize=14, fontweight='bold'
)
for tick, coast in zip(ax_b.get_yticklabels(), COAST_ORDER):
    tick.set_color(COAST_COLORS[coast])
ax_b.set_ylabel('')
ax_b.set_ylim(-0.45, len(COAST_ORDER) - 1 + 0.45)
ax_b.invert_yaxis()

ax_b.set_xlim(0, 100)
ax_b.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
ax_b.set_xlabel('Observed months in drought events (%)', fontsize=14, fontweight='bold')
ax_b.tick_params(axis='x', labelsize=14)
ax_b.tick_params(axis='y', labelsize=14)
ax_b.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.25)
ax_b.set_axisbelow(True)

ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_visible(False)
ax_b.spines['bottom'].set_visible(True)
ax_b.tick_params(axis='y', length=0, pad=8)
ax_b.tick_params(axis='x', length=5, width=1)

# Mean labels
for coast in COAST_ORDER:
    sub = site_df[site_df['coastline'] == coast]
    mean_val = sub['drought_event_exposure_percent'].mean()
    y_pos = coast_y[coast]
    ax_b.text(1.015, y_pos, f'Mean = {mean_val:.1f}%',
              transform=ax_b.get_yaxis_transform(),
              va='center', ha='left', fontsize=13, fontweight='bold',
              color=COAST_COLORS[coast], clip_on=False)

ax_b.text(-0.02, 1.02, 'b)', transform=ax_b.transAxes,
          fontsize=18, fontweight='bold', va='bottom', ha='left')

# Adjust layout – increased bottom margin and hspace to separate panels
plt.subplots_adjust(
    left=0.25,
    right=0.84,
    bottom=0.15,       # more bottom space
    top=0.92,
    hspace=0.6         # increased spacing between panels
)

# ============================================================================
# 9. SAVE COMBINED FIGURE
# ============================================================================
output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling\aug_figure"
os.makedirs(output_dir, exist_ok=True)

output_png = os.path.join(output_dir, 'Figure2_combined.png')
fig.savefig(output_png, dpi=600, bbox_inches='tight', facecolor='white')

print(f"\nCombined figure saved to: {output_png}")
plt.show()

print("\n" + "="*70)
print("SCRIPT FINISHED SUCCESSFULLY")
print("="*70)