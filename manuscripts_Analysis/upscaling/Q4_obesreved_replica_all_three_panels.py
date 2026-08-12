#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 2 – Three panels:

a) Time‑series of cumulative upscaled SPEI‑3 by coastline (4 subplots).
b) Distribution of observed SPEI‑3 classes by coastline (stacked bar).
c) Site‑level drought‑event exposure (lollipop chart).

All panels use the same coastline order and colours.
"""

import os
import sys
import math
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter, PercentFormatter, MaxNLocator
import matplotlib.gridspec as gridspec

# ============================================================================
# 1. COASTLINE CLASSIFICATION FUNCTIONS
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
# Coast colours – identical to study‑area figure
COAST_COLORS = {
    "AK Coast":       "#009E73",   # bluish green
    "Pacific Coast":  "#0072B2",   # blue
    "Gulf Coast":     "#4D4D4D",   # dark charcoal
    "Atlantic Coast": "#A50F15"    # dark red
}

COAST_ORDER = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]
COAST_LABELS = {
    "AK Coast": "Alaska Coast",
    "Pacific Coast": "Pacific Coast",
    "Gulf Coast": "Gulf Coast",
    "Atlantic Coast": "Atlantic Coast"
}

# SPEI‑class colours – revised for distinct blues and a pinkish orange
SPEI_COLORS = {
    '≤ -1.5':      '#D7191C',    # bright red (different from Atlantic)
    '-1.5 to -1.0': '#D35E7D',   # pinkish‑coral (colour‑blind friendly)
    '-1.0 to 1.0':  '#d9d9d9',   # neutral grey
    '1.0 to 1.5':   '#6A8DC4',   # greyish‑blue (distinct from Pacific)
    '≥ 1.5':        '#3A4E7A'    # dark navy‑purple (clearly different from Pacific)
}
CLASS_ORDER = ['≤ -1.5', '-1.5 to -1.0', '-1.0 to 1.0', '1.0 to 1.5', '≥ 1.5']

# ============================================================================
# 3. PATHS
# ============================================================================
BASE_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
MONTHLY_CUM_FILE = BASE_DIR / "Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv"
ANNUAL_CUM_FILE  = BASE_DIR / "Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv"
OBS_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = BASE_DIR / "aug_figure"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 4. LOAD UPSCALED CUMULATIVE DATA
# ============================================================================
print("Loading upscaled cumulative SPEI data...")
try:
    monthly_cum = pd.read_csv(MONTHLY_CUM_FILE)
    annual_cum = pd.read_csv(ANNUAL_CUM_FILE)
    monthly_cum["date"] = pd.to_datetime(
        monthly_cum["year"].astype(str) + "-" +
        monthly_cum["month"].astype(str).str.zfill(2) + "-15"
    )
    annual_cum["date"] = pd.to_datetime(
        annual_cum["year"].astype(str) + "-07-01"
    )
    monthly_cum = monthly_cum[monthly_cum["coast_region"].isin(COAST_ORDER)].copy()
    annual_cum = annual_cum[annual_cum["coast_region"].isin(COAST_ORDER)].copy()
    print("Upscaled cumulative data loaded successfully.")
except FileNotFoundError as e:
    print(f"WARNING: Could not load upscaled data: {e}")
    monthly_cum = None
    annual_cum = None

# ============================================================================
# 5. LOAD OBSERVATIONAL DATA
# ============================================================================
print("Loading observational data...")
try:
    df_raw = pd.read_csv(OBS_DATA_PATH)
except FileNotFoundError:
    print(f"ERROR: File not found at:\n{OBS_DATA_PATH}")
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

# Assign coastlines
sites = df[['site_name', 'lat', 'long']].drop_duplicates()
sites['coast_raw'] = sites.apply(
    lambda row: assign_coast_region(row['lat'], row['long']), axis=1
)
sites['coastline'] = sites['coast_raw']
df = df.merge(sites[['site_name', 'coastline']], on='site_name', how='left')
if df['coastline'].isnull().any():
    print("ERROR: Some sites have no coastline.")
    sys.exit(1)

df_common = df[np.isfinite(df['SPEI_3'])].copy()
print(f"Observational dataset: {len(df_common)} rows, {df_common['site_name'].nunique()} sites")

# ============================================================================
# 6. PREPARE DATA FOR PANEL b (SPEI classes)
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
    sub = df_common[df_common['coastline'] == coast]
    total = len(sub)
    if total == 0:
        pct = {c: 0 for c in CLASS_ORDER}
    else:
        counts = sub['spei_class'].value_counts()
        pct = {c: (counts.get(c, 0) / total) * 100 for c in CLASS_ORDER}
    percentages[coast] = pct

# ============================================================================
# 7. PREPARE DATA FOR PANEL c (drought exposure)
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
    exposure_pct = (n_event_months / n_total) * 100 if n_total > 0 else 0

    site_results.append({
        'site_name': site_name,
        'coastline': group['coastline'].iloc[0],
        'drought_event_exposure_percent': exposure_pct
    })

site_df = pd.DataFrame(site_results)

# ============================================================================
# 8. CREATE THREE‑PANEL FIGURE WITH SPACER ROWS AND FIXED PERCENTAGE LABELS
# ============================================================================
print("Creating three‑panel figure with spacer rows and hidden‑label fixes...")

# ----- Font-size constants -----
TICK_FS = 26
AXIS_LABEL_FS = 30
SUBPLOT_TITLE_FS = 27
PANEL_LABEL_FS = 33
SHARED_YLABEL_FS = 30
PERCENT_FS = 20
SMALL_PERCENT_FS = 15
LEGEND_FS = 24
LEGEND_TITLE_FS = 26
COAST_LABEL_FS = 25
MEAN_LABEL_FS = 24

# ----- Figure setup with blank spacer rows -----
fig = plt.figure(figsize=(14, 22))  # increased height from 20 to 22

# Eight rows:
# rows 0–3 = Panel A (now with increased height ratios)
# row 4    = blank spacer
# row 5    = Panel B
# row 6    = blank spacer
# row 7    = Panel C
gs = gridspec.GridSpec(
    8, 1,
    height_ratios=[
        1.25, 1.25, 1.25, 1.25,   # Panel A – increased from 1.0 to 1.25 each
        0.50,                     # spacer
        2.15,                     # Panel B
        0.55,                     # spacer
        2.65                      # Panel C
    ],
    hspace=0.30,                  # increased from 0.20 to 0.30
    left=0.18,
    right=0.84,
    bottom=0.06,
    top=0.97
)

# ---- Panel A: time‑series subplots (rows 0-3) ----
ax_time = []
for i, coast in enumerate(COAST_ORDER):
    ax = fig.add_subplot(gs[i])
    ax_time.append(ax)

    if monthly_cum is not None and annual_cum is not None:
        sub_m = monthly_cum[monthly_cum["coast_region"] == coast]
        sub_a = annual_cum[annual_cum["coast_region"] == coast]
        if not sub_m.empty:
            ax.plot(sub_m["date"], sub_m["cumulative_SPEI3"],
                    color=COAST_COLORS[coast], linewidth=1.2, alpha=0.88)
        if not sub_a.empty:
            ax.plot(sub_a["date"], sub_a["cumulative_annual_SPEI3"],
                    color=COAST_COLORS[coast], linewidth=2.8)
        ax.axhline(0, color="0.35", linewidth=0.7)
        if i == 0:
            date_min = monthly_cum["date"].min()
            date_max = monthly_cum["date"].max()
        ax.set_xlim(date_min, date_max)
    else:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)

    ax.set_title(COAST_LABELS[coast], fontweight='bold', fontsize=SUBPLOT_TITLE_FS)
    ax.tick_params(axis='both', labelsize=TICK_FS)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.grid(axis="y", color="0.90", linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis='x', length=8, width=2)
    ax.tick_params(axis='y', length=6, width=1.5)

    if i < 3:
        ax.tick_params(axis='x', labelbottom=False)
    else:
        ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.set_xlabel('Year', fontsize=AXIS_LABEL_FS, fontweight='bold', labelpad=8)

# ---- Draw canvas so subplot positions are final ----
fig.canvas.draw()

# ---- Shared y‑axis label for Panel A (moved further left) ----
panel_a_positions = [ax.get_position() for ax in ax_time]
panel_a_top = max(pos.y1 for pos in panel_a_positions)
panel_a_bottom = min(pos.y0 for pos in panel_a_positions)
panel_a_left = min(pos.x0 for pos in panel_a_positions)
panel_a_center_y = (panel_a_top + panel_a_bottom) / 2

fig.text(
    panel_a_left - 0.090,          # was -0.060, now further left
    panel_a_center_y,
    'Cumulative mean SPEI-3',
    rotation='vertical',
    va='center',
    ha='center',
    fontsize=SHARED_YLABEL_FS,
    fontweight='bold'
)

# ---- Panel label 'a)' ----
ax_time[0].text(-0.10, 1.02, 'a)', transform=ax_time[0].transAxes,
                fontsize=PANEL_LABEL_FS, fontweight='bold', va='bottom', ha='left')

# ---- Panel B (stacked bar) – now on row 5 (index 5) ----
ax_b = fig.add_subplot(gs[5])
plot_data_b = pd.DataFrame(percentages).T[CLASS_ORDER].loc[COAST_ORDER]
n_sites_b = df_common.groupby('coastline')['site_name'].nunique().reindex(COAST_ORDER)

bottom = np.zeros(len(COAST_ORDER))
stagger_count = np.zeros(len(COAST_ORDER), dtype=int)

for i, cls in enumerate(CLASS_ORDER):
    vals = plot_data_b[cls].values
    ax_b.barh(COAST_ORDER, vals, left=bottom,
              color=SPEI_COLORS[cls], edgecolor='white', linewidth=0.8)
    for j, val in enumerate(vals):
        if val <= 0:
            continue

        center_x = bottom[j] + val / 2
        # ************** MODIFICATION: round to whole number **************
        label = f'{val:.0f}%'
        # ***************************************************************

        coast_for_bar = COAST_ORDER[j]

        # ----- Which segments must be placed inside the bar (even if narrow) -----
        force_inside_large = (
            (coast_for_bar == "Atlantic Coast" and cls in ['-1.5 to -1.0', '1.0 to 1.5']) or
            (coast_for_bar == "Gulf Coast" and cls == '1.0 to 1.5')
        )

        force_inside_small = (
            (coast_for_bar == "Atlantic Coast" and cls in ['≤ -1.5', '≥ 1.5']) or
            (coast_for_bar == "Gulf Coast" and cls == '≥ 1.5') or
            (coast_for_bar == "AK Coast" and cls == '-1.5 to -1.0')
        )

        # ----- Determine text colour and font size for inside placement -----
        if force_inside_small:
            use_fontsize = SMALL_PERCENT_FS
            inside_color = 'white'
        elif force_inside_large:
            use_fontsize = PERCENT_FS
            inside_color = 'white'
        else:
            use_fontsize = PERCENT_FS
            if cls in ['≤ -1.5', '-1.5 to -1.0', '1.0 to 1.5', '≥ 1.5']:
                inside_color = 'white'
            else:
                inside_color = 'black'

        # ----- Place label inside if it's forced or bar width is ≥5% -----
        if (val >= 5) or force_inside_large or force_inside_small:
            ax_b.text(
                center_x, j, label,
                ha='center', va='center',
                fontsize=use_fontsize,
                fontweight='bold',
                color=inside_color,
                clip_on=False,
                zorder=5
            )
        else:
            # ----- Narrow segment: place outside with annotation -----
            off_y = 18 if stagger_count[j] % 2 == 0 else -18
            vertical_alignment = 'bottom' if stagger_count[j] % 2 == 0 else 'top'
            stagger_count[j] += 1

            ax_b.annotate(
                label,
                xy=(center_x, j),
                xytext=(0, off_y),
                textcoords='offset points',
                ha='center',
                va=vertical_alignment,
                fontsize=SMALL_PERCENT_FS,
                fontweight='bold',
                color='black',
                annotation_clip=False,
                clip_on=False,
                zorder=10,
                bbox=dict(
                    facecolor='white',
                    edgecolor='none',
                    alpha=0.92,
                    pad=0.5
                ),
                arrowprops=dict(
                    arrowstyle='-',
                    linewidth=0.7,
                    color='0.35',
                    shrinkA=1,
                    shrinkB=1
                )
            )
    bottom += vals

ax_b.set_xlim(0, 100)
ax_b.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
ax_b.set_xlabel('Observed site-months (%)', fontsize=AXIS_LABEL_FS, fontweight='bold')
ax_b.set_ylabel('')
ax_b.tick_params(axis='y', labelsize=TICK_FS)
ax_b.tick_params(axis='x', labelsize=TICK_FS)
ax_b.set_yticks(range(len(COAST_ORDER)))
ax_b.set_yticklabels(
    [f"{COAST_LABELS[coast]} (N = {int(n_sites_b[coast])})" for coast in COAST_ORDER],
    fontsize=COAST_LABEL_FS, fontweight='bold'
)
for tick, coast in zip(ax_b.get_yticklabels(), COAST_ORDER):
    tick.set_color(COAST_COLORS[coast])
ax_b.invert_yaxis()

# Extra vertical room for outside labels on Alaska/Atlantic bars
ax_b.set_ylim(len(COAST_ORDER) - 0.35, -0.70)
ax_b.tick_params(axis='y', pad=12)
ax_b.xaxis.labelpad = 12

ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)
ax_b.spines['left'].set_visible(False)
ax_b.spines['bottom'].set_visible(True)

# ---- Legend on the right side ----
legend_patches = [mpatches.Patch(color=SPEI_COLORS[cls], label=cls) for cls in CLASS_ORDER]
leg_b = ax_b.legend(handles=legend_patches,
                    loc='center left',
                    bbox_to_anchor=(1.02, 0.5),
                    ncol=1,
                    fontsize=LEGEND_FS,
                    title='SPEI-3')
leg_b.get_title().set_fontsize(LEGEND_TITLE_FS)
leg_b.get_title().set_fontweight('bold')

plt.subplots_adjust(right=0.84)

ax_b.text(-0.02, 1.08, 'b)', transform=ax_b.transAxes,
          fontsize=PANEL_LABEL_FS, fontweight='bold', va='bottom', ha='left', clip_on=False)

# ---- Panel C (lollipop) – now on row 7 (index 7) ----
ax_c = fig.add_subplot(gs[7])
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
        ax_c.plot([x, x], [base_y, y_pos], color='gray', linewidth=0.8, alpha=0.6, zorder=1)
        ax_c.scatter(x, y_pos, s=120, color=COAST_COLORS[coast],
                     edgecolor='white', linewidth=0.5, zorder=3)

    ax_c.axhline(y=base_y, color='black', linewidth=0.5, alpha=0.3, zorder=0)
    mean_val = sub['drought_event_exposure_percent'].mean()
    ax_c.scatter(mean_val, base_y, s=180, marker='D', color='black', zorder=5)

# ---- Formatting Panel C (no N= in labels) ----
ax_c.set_yticks(range(len(COAST_ORDER)))
ax_c.set_yticklabels(
    [COAST_LABELS[coast] for coast in COAST_ORDER],
    fontsize=COAST_LABEL_FS, fontweight='bold'
)
for tick, coast in zip(ax_c.get_yticklabels(), COAST_ORDER):
    tick.set_color(COAST_COLORS[coast])
ax_c.set_ylabel('')
ax_c.set_ylim(-0.45, len(COAST_ORDER) - 1 + 0.45)
ax_c.invert_yaxis()
ax_c.set_xlim(0, 100)
ax_c.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
ax_c.set_xlabel('Observed months in drought events (%)', fontsize=AXIS_LABEL_FS, fontweight='bold')
ax_c.tick_params(axis='x', labelsize=TICK_FS)
ax_c.tick_params(axis='y', labelsize=TICK_FS)
ax_c.grid(axis='x', linestyle='--', linewidth=0.6, alpha=0.25)
ax_c.set_axisbelow(True)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.spines['left'].set_visible(False)
ax_c.spines['bottom'].set_visible(True)
ax_c.tick_params(axis='y', length=0, pad=8)
ax_c.tick_params(axis='x', length=5, width=1)

# ---- Mean labels in Panel C (moved further right) ----
for coast in COAST_ORDER:
    sub = site_df[site_df['coastline'] == coast]
    mean_val = sub['drought_event_exposure_percent'].mean()
    y_pos = coast_y[coast]
    ax_c.text(1.08, y_pos, f'Mean = {mean_val:.1f}%',   # was 1.015, now 1.08
              transform=ax_c.get_yaxis_transform(),
              va='center', ha='left', fontsize=MEAN_LABEL_FS, fontweight='bold',
              color=COAST_COLORS[coast], clip_on=False)

ax_c.text(-0.02, 1.08, 'c)', transform=ax_c.transAxes,
          fontsize=PANEL_LABEL_FS, fontweight='bold', va='bottom', ha='left', clip_on=False)

# ============================================================================
# FINAL FONT-SIZE ENFORCEMENT
# ============================================================================

# Panel A tick labels
for ax in ax_time:
    ax.tick_params(
        axis='y',
        which='major',
        labelsize=TICK_FS
    )

# Only the final Panel A subplot displays x-axis labels
ax_time[-1].tick_params(
    axis='x',
    which='major',
    labelsize=TICK_FS
)

# Panel B x-axis ticks
ax_b.tick_params(
    axis='x',
    which='major',
    labelsize=TICK_FS
)

# Panel B coastline labels
for label in ax_b.get_yticklabels():
    label.set_fontsize(COAST_LABEL_FS)
    label.set_fontweight('bold')

# Panel B legend
for legend_text in leg_b.get_texts():
    legend_text.set_fontsize(LEGEND_FS)

leg_b.get_title().set_fontsize(LEGEND_TITLE_FS)
leg_b.get_title().set_fontweight('bold')

# Panel C x-axis ticks
ax_c.tick_params(
    axis='x',
    which='major',
    labelsize=TICK_FS
)

# Panel C coastline labels
for label in ax_c.get_yticklabels():
    label.set_fontsize(COAST_LABEL_FS)
    label.set_fontweight('bold')

# Ensure the final text changes are rendered before saving
fig.canvas.draw()

# ============================================================================
# 9. SAVE PNG, PDF, and SVG
# ============================================================================
output_png = OUTPUT_DIR / "Figure2_three_panels.png"
output_pdf = OUTPUT_DIR / "Figure2_three_panels.pdf"
output_svg = OUTPUT_DIR / "Figure2_three_panels.svg"

fig.savefig(
    output_png,
    dpi=600,
    bbox_inches='tight',
    pad_inches=0.15,
    facecolor='white',
    transparent=False
)

fig.savefig(
    output_pdf,
    bbox_inches='tight',
    pad_inches=0.15,
    facecolor='white'
)

fig.savefig(
    output_svg,
    bbox_inches='tight',
    pad_inches=0.15,
    facecolor='white'
)

print(f"PNG saved to: {output_png.resolve()}")
print(f"PDF saved to: {output_pdf.resolve()}")
print(f"SVG saved to: {output_svg.resolve()}")
plt.show()

print("\n" + "="*70)
print("SCRIPT FINISHED SUCCESSFULLY")
print("="*70)