# -*- coding: utf-8 -*-
"""
STUDY AREA MAP – FINAL WITH PAIRED T:ET & SPEI ROWS
- T:ET and SPEI-3 are computed from the same site-month rows
- Alaska time‑series moved higher to avoid overlap
- All fonts, markers, and labels enlarged
- UPDATED: thicker lines, PDF export, markers added
"""

import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS – UPDATED OUTPUT DIRECTORY
# ============================================================================
MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_PNG = os.path.join(OUTPUT_DIR, "study_area_figure.png")
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "study_area_figure.pdf")   # NEW: PDF export

# ============================================================================
# COAST CLASSIFIER (unchanged)
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
# COAST COLOURS & SHORT NAMES (unchanged)
# ============================================================================
COAST_COLORS = {
    'Atlantic Coast': '#2E8B57',
    'Gulf Coast':     '#00CED1',
    'Pacific Coast':  '#DC143C',
    'AK Coast':       '#8B4513'
}
COAST_SHORT = {
    'Atlantic Coast': 'Atlantic',
    'Pacific Coast':  'Pacific',
    'Gulf Coast':     'Gulf',
    'AK Coast':       'Alaska'
}
COAST_ORDER = ['Atlantic Coast', 'Pacific Coast', 'Gulf Coast', 'AK Coast']

# ============================================================================
# ECOSYSTEM MARKERS (unchanged)
# ============================================================================
ECOSYSTEM_MARKERS = {
    'Upland': 'o',
    'Freshwater': 's',
    'Saline': '^'
}
LEGEND_ORDER = ['Upland', 'Freshwater', 'Saline']
MARKER_EDGE_WIDTH = 0.5

# ============================================================================
# MAP EXTENTS (unchanged)
# ============================================================================
CONUS_EXTENT = [-130, -62, 22, 52]
ALASKA_EXTENT = [-170, -128, 49, 74]
ALASKA_INSET_POSITION = [0.11, 0.61, 0.25, 0.28]

# ============================================================================
# DATA LOADING & FILTERING (unchanged)
# ============================================================================
def load_and_filter_monthly_data(filepath):
    print(f"📂 Loading monthly data from {filepath} ...")
    df = pd.read_csv(filepath)
    
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].where(df[col].isna(), df[col].astype(str).str.strip())
        df[col] = df[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})
    
    required = ['site_name', 'Year', 'month', 'water_class', 'lat', 'long',
                'Trans_ratio', 'WUE_tra']
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found")
    
    spei_candidates = ['SPEI_3', 'spei3', 'SPEI3']
    spei_col = None
    for cand in spei_candidates:
        if cand in df.columns:
            spei_col = cand
            break
    if spei_col is not None:
        df.rename(columns={spei_col: 'spei3'}, inplace=True)
        print(f"✅ Found SPEI column: {spei_col} → renamed to 'spei3'")
    else:
        df['spei3'] = np.nan
        print("⚠️  No 3‑month SPEI column found")
    
    df = df.dropna(subset=required).copy()
    df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])].copy()
    df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) &
            np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) &
            (df['Trans_ratio'] <= 1) & np.isfinite(df['WUE_tra'])].copy()
    
    df['spei3'] = pd.to_numeric(df['spei3'], errors='coerce')
    
    print(f"  Filtered to {len(df)} monthly rows, {df['site_name'].nunique()} unique sites")
    return df

def get_site_summary(df_monthly):
    summary = df_monthly.groupby('site_name').agg(
        lat=('lat', 'first'),
        long=('long', 'first'),
        water_class=('water_class', 'first'),
        Trans_ratio_median=('Trans_ratio', 'median')
    ).reset_index()
    return summary

# ============================================================================
# PLOTTING FUNCTIONS (map – unchanged)
# ============================================================================
def plot_sites_by_coast(ax, df_sites):
    for coast in COAST_ORDER:
        sub = df_sites[df_sites['coast_region'] == coast]
        if sub.empty:
            continue
        color = COAST_COLORS[coast]
        for eco, marker in ECOSYSTEM_MARKERS.items():
            eco_sub = sub[sub['water_class'] == eco]
            if eco_sub.empty:
                continue
            ax.scatter(eco_sub['long'], eco_sub['lat'],
                       c=[color], marker=marker, s=350,
                       edgecolor=color, linewidth=MARKER_EDGE_WIDTH,
                       zorder=6,
                       transform=ccrs.PlateCarree())

def add_coast_label(ax, coast, lon, lat, n_sites):
    label_text = f'{COAST_SHORT[coast]}\nN = {n_sites}'
    ax.text(lon, lat, label_text,
            fontsize=32, fontweight='bold', ha='center', va='center',
            color='black',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      alpha=0.9, edgecolor='black', linewidth=1.2),
            transform=ccrs.PlateCarree(), zorder=5, linespacing=1.5)

def add_alaska_inset(ax_main, df_alaska, n_sites):
    ax_inset = ax_main.inset_axes(ALASKA_INSET_POSITION,
                                  projection=ccrs.PlateCarree())
    ax_inset.set_extent(ALASKA_EXTENT, crs=ccrs.PlateCarree())
    ax_inset.coastlines(resolution='50m', linewidth=1.5,
                        color='darkgray', zorder=2)
    ax_inset.add_feature(cfeature.OCEAN, facecolor='#B8D4E8',
                         alpha=0.3, zorder=0)
    ax_inset.add_feature(cfeature.LAND, facecolor='#f0f0f0',
                         alpha=0.3, zorder=0)
    
    for coast in COAST_ORDER:
        sub = df_alaska[df_alaska['coast_region'] == coast]
        if sub.empty:
            continue
        color = COAST_COLORS[coast]
        for eco, marker in ECOSYSTEM_MARKERS.items():
            eco_sub = sub[sub['water_class'] == eco]
            if eco_sub.empty:
                continue
            ax_inset.scatter(eco_sub['long'], eco_sub['lat'],
                             c=[color], marker=marker, s=220,
                             edgecolor=color, linewidth=MARKER_EDGE_WIDTH,
                             zorder=6, transform=ccrs.PlateCarree())
    
    ax_inset.text(0.99, 0.95, f'Alaska\nN = {n_sites}', transform=ax_inset.transAxes,
                  fontsize=32, fontweight='bold', ha='right', va='center',
                  color='black',
                  bbox=dict(boxstyle='round,pad=0.22', facecolor='white',
                            alpha=0.9, edgecolor='black', linewidth=1.2),
                  zorder=4)
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)
    return ax_inset

def add_ecosystem_legend(ax, ecosystem_counts):
    leg_elements = []
    for eco in LEGEND_ORDER:
        if eco in ECOSYSTEM_MARKERS:
            count = ecosystem_counts.get(eco, 0)
            leg_elements.append(
                Line2D([0], [0], marker=ECOSYSTEM_MARKERS[eco],
                       color='none', markerfacecolor='black',
                       markeredgecolor='black', markersize=26,
                       linestyle='none',
                       label=f'{eco} (N={count})')
            )
    legend = ax.legend(handles=leg_elements, loc='upper right',
                       bbox_to_anchor=(0.90, 0.995),
                       fontsize=30,
                       framealpha=0.95, title='Ecosystem Type',
                       title_fontsize=34)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.2)
    return legend

# ============================================================================
# STACKED TIME‑SERIES PANELS – UPDATED WITH THICKER LINES, MARKERS
# ============================================================================
def add_time_series_panels(fig, df_monthly, coast, rect):
    left, bottom, width, height = rect
    half_height = height / 2.0
    gap = 0.03 * height
    
    ax_top = fig.add_axes([left, bottom + half_height + gap/2, width, half_height - gap/2])
    ax_bottom = fig.add_axes([left, bottom, width, half_height - gap/2])
    
    sub = df_monthly[df_monthly['coast_region'] == coast].copy()
    if sub.empty:
        for ax in [ax_top, ax_bottom]:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=16)
            ax.set_facecolor('white')
            for spine in ax.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1.0)
        return
    
    sub['date'] = pd.to_datetime(sub['Year'].astype(str) + '-' +
                                 sub['month'].astype(str) + '-01')
    
    # Retain only rows with both T:ET and SPEI-3
    paired_sub = sub.dropna(subset=['Trans_ratio', 'spei3']).copy()
    
    if paired_sub.empty:
        for ax in [ax_top, ax_bottom]:
            ax.text(0.5, 0.5, 'No paired data', ha='center', va='center',
                    transform=ax.transAxes, fontsize=16, color='gray')
            ax.set_facecolor('white')
            for spine in ax.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1.0)
        return
    
    # Aggregate by month – both means from the same rows
    monthly_paired = (
        paired_sub.groupby('date')
        .agg(
            Trans_ratio=('Trans_ratio', 'mean'),
            spei3=('spei3', 'mean')
        )
        .sort_index()
    )
    
    # Create full monthly date range from first to last available month
    date_range = pd.date_range(
        start=monthly_paired.index.min(),
        end=monthly_paired.index.max(),
        freq='MS'
    )
    monthly_paired = monthly_paired.reindex(date_range)
    
    # ---- DIAGNOSTIC CHECK (uncomment to verify) ----
    # gap_mismatch = (
    #     monthly_paired['Trans_ratio'].isna()
    #     != monthly_paired['spei3'].isna()
    # )
    # print(coast, "months with mismatched gaps:", gap_mismatch.sum())
    # assert not gap_mismatch.any(), f"{coast}: T:ET and SPEI have different missing months"
    
    monthly_te = monthly_paired['Trans_ratio']
    monthly_spei = monthly_paired['spei3']
    
    # ---- X limits ----
    x_start = pd.Timestamp(date_range.min().year, 1, 1)
    x_end = pd.Timestamp('2025-12-31')
    
    # ---- Dynamic y‑limits ----
    te_valid = monthly_te.dropna()
    if len(te_valid) > 0:
        te_min, te_max = te_valid.min(), te_valid.max()
        te_pad = (te_max - te_min) * 0.1 if te_max > te_min else 0.05
        te_ylim = (max(0, te_min - te_pad), min(1, te_max + te_pad))
    else:
        te_ylim = (0.10, 0.82)
    
    spei_valid = monthly_spei.dropna()
    if len(spei_valid) > 0:
        spei_min, spei_max = spei_valid.min(), spei_valid.max()
        spei_pad = (spei_max - spei_min) * 0.1 if spei_max > spei_min else 0.2
        spei_ylim = (spei_min - spei_pad, spei_max + spei_pad)
    else:
        spei_ylim = (-2.5, 2.5)
    
    # ---- Plot T:ET – linewidth 3.0, with markers ----
    ax_top.plot(monthly_te.index, monthly_te.values,
                color=COAST_COLORS[coast], linestyle='-', linewidth=3.0,
                marker='o', markersize=4, markeredgewidth=0, markerfacecolor=COAST_COLORS[coast])
    ax_top.set_ylim(te_ylim)
    ax_top.yaxis.set_major_locator(MaxNLocator(5, prune='both'))
    ax_top.set_ylabel('T:ET', fontsize=26, color=COAST_COLORS[coast], weight='bold')
    ax_top.tick_params(axis='y', labelsize=22, colors=COAST_COLORS[coast])
    ax_top.set_xlim(x_start, x_end)
    ax_top.grid(True, alpha=0.25, linewidth=0.5)
    
    # ---- Plot SPEI – linewidth 3.0, with markers ----
    if not monthly_spei.isna().all():
        ax_bottom.plot(monthly_spei.index, monthly_spei.values,
                       color='black', linestyle='-', linewidth=3.0,
                       marker='o', markersize=4, markeredgewidth=0, markerfacecolor='black')
    else:
        ax_bottom.text(0.5, 0.5, 'SPEI N/A', transform=ax_bottom.transAxes,
                       ha='center', va='center', fontsize=16, color='gray')
    ax_bottom.set_ylim(spei_ylim)
    ax_bottom.yaxis.set_major_locator(MaxNLocator(5, prune='both'))
    ax_bottom.set_ylabel('SPEI', fontsize=26, color='black', weight='bold')
    ax_bottom.tick_params(axis='y', labelsize=22, colors='black')
    ax_bottom.set_xlim(x_start, x_end)
    ax_bottom.grid(True, alpha=0.25, linewidth=0.5)
    
    # ---- Shared x‑axis formatting ----
    for ax in [ax_top, ax_bottom]:
        ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.YearLocator(base=1))
        ax.tick_params(axis='x', which='major', labelsize=22, length=6, width=1.2, pad=8)
        ax.tick_params(axis='x', which='minor', length=4, width=0.8)
        plt.setp(ax.get_xticklabels(), rotation=25, ha='right', rotation_mode='anchor', weight='bold')
        ax.set_facecolor('white')
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(1.2)
    
    plt.setp(ax_top.get_yticklabels(), weight='bold')
    plt.setp(ax_bottom.get_yticklabels(), weight='bold')
    
    ax_top.set_xlabel('')
    ax_top.tick_params(axis='x', labelbottom=False)
    ax_bottom.set_xlabel('')
    
    return ax_top, ax_bottom

# ============================================================================
# MAIN FIGURE (unchanged except for output format)
# ============================================================================
def create_figure(df_monthly, df_sites, coast_counts):
    fig = plt.figure(figsize=(30, 20))
    
    # ---- Map axes ----
    map_ax = fig.add_axes([0.12, 0.06, 0.72, 0.92], projection=ccrs.PlateCarree())
    map_ax.set_extent(CONUS_EXTENT, crs=ccrs.PlateCarree())
    map_ax.coastlines(resolution='50m', linewidth=1.2, color='#666666', zorder=2)
    map_ax.add_feature(cfeature.OCEAN, facecolor='#E0F0FA', alpha=0.5, zorder=0)
    map_ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', alpha=0.05, zorder=0)
    map_ax.add_feature(cfeature.STATES, linewidth=0.8, edgecolor='#999999', alpha=0.4, zorder=1)
    
    # ---- Coast labels ----
    label_positions = {
        'Pacific Coast': (-127.0, 40.5),
        'Gulf Coast': (-88, 26.5),
        'Atlantic Coast': (-75, 30.5)
    }
    for coast in ['Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
        if coast not in coast_counts:
            continue
        lon, lat = label_positions[coast]
        add_coast_label(map_ax, coast, lon, lat, coast_counts[coast])
    
    plot_sites_by_coast(map_ax, df_sites)
    
    # ---- Alaska inset ----
    alaska_sites = df_sites[df_sites['lat'] >= 50]
    alaska_ax = None
    if not alaska_sites.empty:
        alaska_ax = add_alaska_inset(map_ax, alaska_sites, coast_counts.get('AK Coast', 0))
    
    # ---- Ecosystem legend ----
    eco_counts = df_sites['water_class'].value_counts().to_dict()
    add_ecosystem_legend(map_ax, eco_counts)
    
    # ---- Gridlines ----
    gl = map_ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.3, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 26}
    gl.ylabel_style = {'size': 26}
    
    fig.canvas.draw()
    map_box = map_ax.get_position()
    alaska_box = alaska_ax.get_position() if alaska_ax else None
    
    # ---- Panel sizes ----
    ts_width = 0.18
    ts_height = 0.22
    margin = 0.015
    
    # ---- Pacific – fixed left position ----
    pacific_x = 0.005
    pacific_y = map_box.y0 + 0.01 * map_box.height
    pacific_rect = [pacific_x, pacific_y, ts_width, ts_height]
    
    # ---- Alaska – higher ----
    if alaska_box:
        alaska_x = alaska_box.x0 + 0.12 * alaska_box.width
        alaska_y = min(1.01 - ts_height, alaska_box.y1 + 0.25)
        alaska_rect = [alaska_x, alaska_y, ts_width, ts_height]
    else:
        alaska_rect = [0.08, 0.75, ts_width, ts_height]
    
    # ---- Gulf – below the map ----
    gulf_x = map_box.x0 + 0.45 * map_box.width
    gulf_y = map_box.y0 - ts_height - 0.025
    gulf_rect = [gulf_x, gulf_y, ts_width, ts_height]
    
    # ---- Atlantic – right of the map ----
    atlantic_x = min(0.90, map_box.x1 + margin + 0.5)
    atlantic_y = map_box.y0 + 0.15 * map_box.height
    atlantic_rect = [atlantic_x, atlantic_y, ts_width, ts_height]
    
    # ---- Add panels ----
    add_time_series_panels(fig, df_monthly, 'Pacific Coast', pacific_rect)
    add_time_series_panels(fig, df_monthly, 'AK Coast', alaska_rect)
    add_time_series_panels(fig, df_monthly, 'Gulf Coast', gulf_rect)
    add_time_series_panels(fig, df_monthly, 'Atlantic Coast', atlantic_rect)
    
    return fig

# ============================================================================
# MAIN – SAVES BOTH PNG AND PDF, DISPLAYS FIGURE
# ============================================================================
def main():
    print("="*80)
    print("STUDY AREA MAP – PAIRED T:ET & SPEI ROWS")
    print("Saving to:", OUTPUT_DIR)
    print("="*80)
    
    df_monthly = load_and_filter_monthly_data(MONTHLY_DATA_PATH)
    
    print("Assigning coast regions ...")
    df_monthly['coast_region'] = df_monthly.apply(
        lambda row: assign_coast_region(row['lat'], row['long']), axis=1
    )
    df_monthly = df_monthly.dropna(subset=['coast_region']).copy()
    print(f"  {df_monthly['site_name'].nunique()} sites assigned to coasts")
    
    df_sites = get_site_summary(df_monthly)
    site_coast = df_monthly.groupby('site_name')['coast_region'].first().reset_index()
    df_sites = df_sites.merge(site_coast, on='site_name', how='left')
    
    coast_counts = df_sites['coast_region'].value_counts().to_dict()
    print("\nSite counts per coast:")
    for coast, count in coast_counts.items():
        print(f"  {coast}: {count}")
    
    print("\nCreating figure ...")
    fig = create_figure(df_monthly, df_sites, coast_counts)
    
    # ---- Save both PNG and PDF ----
    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, bbox_inches='tight', facecolor='white')
    print(f"\n✅ Saved PNG:  {OUTPUT_PNG}")
    print(f"✅ Saved PDF:  {OUTPUT_PDF}")
    
    # ---- Display in console (Spyder / interactive) ----
    plt.show()

if __name__ == "__main__":
    main()