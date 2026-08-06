# -*- coding: utf-8 -*-
"""
STUDY AREA FIGURE – ULTRA ENLARGED (TIME‑SERIES ORDER: ALASKA, PACIFIC, GULF, ATLANTIC)
- Figure: 56 x 26
- Marker size: 800 / 400
- Tick labels: 40 (time‑series), 36 (map)
- T:ET label: 34
- Pie diameter: 4.0
- Alaska pie at lon=-169.8 (leftmost inside inset)
- Time‑series order: Alaska, Pacific, Gulf, Atlantic (top to bottom)
- PNG and PDF export
"""

import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS
# ============================================================================
MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_PNG = os.path.join(OUTPUT_DIR, "study_area_figure.png")
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "study_area_figure.pdf")

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
# COLOUR PALETTES & SHORT NAMES
# ============================================================================
# ============================================================================
# COLOUR PALETTES & SHORT NAMES
# ============================================================================
COAST_COLORS = {
    'Pacific Coast':  '#0072B2',  # blue
    'Gulf Coast':     '#4D4D4D',  # dark charcoal
    'Atlantic Coast': '#A50F15',  # dark red
    'AK Coast':       '#009E73'   # bluish green
}

ECOSYSTEM_COLORS = {
    'Upland':      '#800080',
    'Freshwater':  '#0000FF',
    'Saline':      '#FFA500'
}
# Map order (site markers, coast labels, pies)
COAST_ORDER = ['Pacific Coast', 'Gulf Coast', 'Atlantic Coast', 'AK Coast']
# Time‑series order (new: Alaska, Pacific, Gulf, Atlantic)
TS_ORDER = ['AK Coast', 'Pacific Coast', 'Gulf Coast', 'Atlantic Coast']
ECOSYSTEM_ORDER = ['Upland', 'Freshwater', 'Saline']
COAST_SHORT = {
    'Pacific Coast': 'Pacific',
    'Gulf Coast': 'Gulf',
    'Atlantic Coast': 'Atlantic',
    'AK Coast': 'Alaska'
}

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
    
    df = df.dropna(subset=required).copy()
    df = df[df['water_class'].isin(['Upland', 'Freshwater', 'Saline'])].copy()
    df = df[np.isfinite(df['lat']) & np.isfinite(df['long']) &
            np.isfinite(df['Trans_ratio']) & (df['Trans_ratio'] >= 0) &
            (df['Trans_ratio'] <= 1) & np.isfinite(df['WUE_tra'])].copy()
    
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
# PLOTTING FUNCTIONS (with enlarged sizes)
# ============================================================================

def add_site_markers(ax, df_sites, inset=False):
    size = 800 if not inset else 400
    for coast in COAST_ORDER:
        sub = df_sites[df_sites['coast_region'] == coast]
        if sub.empty:
            continue
        color = COAST_COLORS[coast]
        ax.scatter(
            sub['long'], sub['lat'],
            c=[color], marker='o', s=size,
            edgecolor='white', linewidth=1.0,
            zorder=6, transform=ccrs.PlateCarree()
        )

def add_coast_label(ax, coast, lon, lat, n_sites):
    short = COAST_SHORT[coast]
    ax.text(
        lon, lat,
        f'{short}\nN = {n_sites}',
        fontsize=44, fontweight='bold',
        ha='center', va='center', color='black',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.92,
                  edgecolor='black', linewidth=1.2),
        transform=ccrs.PlateCarree(), zorder=15, linespacing=1.2
    )

def add_alaska_inset(ax_main, df_alaska, n_sites):
    alaska_inset_pos = [0.22, 0.70, 0.20, 0.25]
    ax_inset = ax_main.inset_axes(alaska_inset_pos, projection=ccrs.PlateCarree())
    ax_inset.set_extent([-170, -128, 49, 74], crs=ccrs.PlateCarree())
    ax_inset.coastlines(resolution='50m', linewidth=2.0, color='#666666', zorder=2)
    ax_inset.add_feature(cfeature.OCEAN, facecolor='#E0F0FA', zorder=0)
    ax_inset.add_feature(cfeature.LAND, facecolor='#F5F5F5', zorder=0)
    add_site_markers(ax_inset, df_alaska, inset=True)
    ax_inset.text(
        0.97, 0.99,
        f'Alaska\nN = {n_sites}',
        transform=ax_inset.transAxes,
        fontsize=44, fontweight='bold',
        ha='right', va='top', color='black',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.92,
                  edgecolor='black', linewidth=1.2),
        zorder=15
    )
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    return ax_inset

def add_geo_pie(fig, geo_ax, lon, lat, counts, diameter_in=4.0):
    values = np.array([counts.get(eco, 0) for eco in ECOSYSTEM_ORDER], dtype=float)
    if values.sum() <= 0:
        return None
    colors = [ECOSYSTEM_COLORS[eco] for eco in ECOSYSTEM_ORDER]

    x_projected, y_projected = geo_ax.projection.transform_point(lon, lat, ccrs.PlateCarree())
    x_display, y_display = geo_ax.transData.transform((x_projected, y_projected))
    x_figure, y_figure = fig.transFigure.inverted().transform((x_display, y_display))

    pie_width = diameter_in / fig.get_figwidth()
    pie_height = diameter_in / fig.get_figheight()

    pie_ax = fig.add_axes(
        [x_figure - pie_width/2, y_figure - pie_height/2, pie_width, pie_height],
        zorder=30
    )

    def show_percentage(pct):
        return f'{pct:.0f}%' if pct >= 7 else ''

    wedges, texts, autotexts = pie_ax.pie(
        values,
        colors=colors,
        startangle=90,
        counterclock=False,
        autopct=show_percentage,
        pctdistance=0.55,
        textprops={'fontsize': 38, 'fontweight': 'bold', 'color': 'white'},
        wedgeprops={'edgecolor': 'none', 'linewidth': 0}
    )

    pie_ax.set_aspect('equal')
    pie_ax.set_axis_off()
    return pie_ax

def add_ecosystem_legend(ax):
    handles = [
        Patch(facecolor=ECOSYSTEM_COLORS[eco], edgecolor='none', linewidth=0)
        for eco in ECOSYSTEM_ORDER
    ]
    legend = ax.legend(
        handles=handles,
        labels=ECOSYSTEM_ORDER,
        loc='upper center',
        bbox_to_anchor=(0.50, -0.075),
        ncol=3,
        fontsize=44,
        title='Ecosystem Type',
        title_fontsize=48,
        framealpha=0.95,
        edgecolor='black',
        frameon=True,
        columnspacing=3.0,
        handletextpad=1.2,
        handlelength=2.5
    )
    legend.get_frame().set_linewidth(1.5)
    return legend

def add_time_series_panel(ax, df_monthly, coast, common_xlim, common_ylim=(0,1)):
    sub = df_monthly[df_monthly['coast_region'] == coast].copy()
    if sub.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes, fontsize=44, color='gray')
        return ax

    sub['date'] = pd.to_datetime(dict(year=sub['Year'].astype(int), month=sub['month'].astype(int), day=1))
    monthly = sub.groupby('date')['Trans_ratio'].mean().sort_index()
    full_range = pd.date_range(start=common_xlim[0], end=common_xlim[1], freq='MS')
    monthly = monthly.reindex(full_range)

    ax.plot(monthly.index, monthly.values, color=COAST_COLORS[coast], linewidth=4.0)
    ax.set_xlim(common_xlim)
    ax.set_ylim(common_ylim)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.grid(True, axis='y', alpha=0.25, linewidth=0.8)
    ax.tick_params(axis='both', labelsize=50)
    ax.set_facecolor('white')
    ax.margins(x=0)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)

    ax.text(0.018, 0.93, COAST_SHORT[coast], transform=ax.transAxes, fontsize=44, fontweight='bold', va='top', ha='left')

    ax.xaxis.set_major_locator(mdates.YearLocator(base=5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.YearLocator(base=1))
    ax.tick_params(axis='x', which='major', length=12, width=2.0, pad=12)
    ax.tick_params(axis='x', which='minor', length=8, width=1.2)
    return ax

# ============================================================================
# MAIN FIGURE CREATION
# ============================================================================
def create_figure(df_monthly, df_sites, coast_counts, eco_counts_per_coast):
    fig = plt.figure(figsize=(56, 26), facecolor='white')

    outer_gs = fig.add_gridspec(
        nrows=1, ncols=2,
        width_ratios=[2.0, 1.0],
        left=0.04, right=0.98, bottom=0.11, top=0.94,
        wspace=0.12
    )

    # ---- Map ----
    map_ax = fig.add_subplot(outer_gs[0, 0], projection=ccrs.PlateCarree())
    map_ax.set_extent([-130, -62, 22, 52], crs=ccrs.PlateCarree())
    map_ax.set_anchor('C')
    map_ax.coastlines(resolution='50m', linewidth=2.5, color='#666666', zorder=2)
    map_ax.add_feature(cfeature.OCEAN, facecolor='#E0F0FA', zorder=0)
    map_ax.add_feature(cfeature.LAND, facecolor='#F5F5F5', zorder=0)
    map_ax.add_feature(cfeature.STATES, linewidth=1.5, edgecolor='#999999', alpha=0.55, zorder=1)

    gl = map_ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.30, linestyle='-')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 46}
    gl.ylabel_style = {'size': 46}

    label_positions = {
        'Pacific Coast': (-124.0, 41.0),
        'Gulf Coast': (-87.0, 28.0),
        'Atlantic Coast': (-68.0, 39.0)
    }
    for coast in ['Pacific Coast', 'Gulf Coast', 'Atlantic Coast']:
        lon, lat = label_positions[coast]
        add_coast_label(map_ax, coast, lon, lat, coast_counts.get(coast, 0))

    add_site_markers(map_ax, df_sites, inset=False)

    alaska_ax = None
    alaska_sites = df_sites[df_sites['coast_region'] == 'AK Coast']
    if not alaska_sites.empty:
        alaska_ax = add_alaska_inset(map_ax, alaska_sites, coast_counts.get('AK Coast', 0))

    # ---- Time series (order: Alaska, Pacific, Gulf, Atlantic) ----
    right_gs = outer_gs[0, 1].subgridspec(nrows=4, ncols=1, hspace=0.28)
    ts_axes = [fig.add_subplot(right_gs[i, 0]) for i in range(4)]

    all_dates = pd.to_datetime(dict(year=df_monthly['Year'].astype(int), month=df_monthly['month'].astype(int), day=1))
    common_xlim = (pd.Timestamp(all_dates.min().year, 1, 1), pd.Timestamp(all_dates.max().year, 12, 31))
    common_ylim = (0, 1)

    for i, coast in enumerate(TS_ORDER):
        ax = ts_axes[i]
        add_time_series_panel(ax, df_monthly, coast, common_xlim, common_ylim)
        ax.set_ylabel('T:ET', fontsize=44, fontweight='bold')
        if i == len(TS_ORDER) - 1:
            ax.set_xlabel('Year', fontsize=44, fontweight='bold')
        else:
            ax.set_xlabel('')
            ax.tick_params(axis='x', labelbottom=False)

    fig.canvas.draw()

    # ---- Pie charts ----
    PIE_DIAMETER = 4.0
    add_geo_pie(fig, map_ax, lon=-126.5, lat=33.0, counts=eco_counts_per_coast['Pacific Coast'], diameter_in=PIE_DIAMETER)
    add_geo_pie(fig, map_ax, lon=-91.5, lat=24.5, counts=eco_counts_per_coast['Gulf Coast'], diameter_in=PIE_DIAMETER)
    add_geo_pie(fig, map_ax, lon=-70.0, lat=33.5, counts=eco_counts_per_coast['Atlantic Coast'], diameter_in=PIE_DIAMETER)
    if alaska_ax is not None:
        # Alaska pie remains at far left inside inset (lon=-169.8)
        add_geo_pie(fig, alaska_ax, lon=-176, lat=69.5, counts=eco_counts_per_coast['AK Coast'], diameter_in=PIE_DIAMETER)

    add_ecosystem_legend(map_ax)

    map_ax.text(-0.03, 1.03, '(a)', transform=map_ax.transAxes, fontsize=48, fontweight='bold', va='bottom', ha='left')
    ts_axes[0].text(-0.08, 1.07, '(b)', transform=ts_axes[0].transAxes, fontsize=48, fontweight='bold', va='bottom', ha='left')

    return fig

# ============================================================================
# ECOSYSTEM COUNTS
# ============================================================================
def get_ecosystem_counts(df_sites):
    eco_counts = {}
    for coast in COAST_ORDER:
        sub = df_sites[df_sites['coast_region'] == coast]
        counts = sub['water_class'].value_counts().to_dict()
        for eco in ECOSYSTEM_ORDER:
            counts[eco] = counts.get(eco, 0)
        eco_counts[coast] = counts
    return eco_counts

# ============================================================================
# MAIN
# ============================================================================
def main():
    print("="*80)
    print("STUDY AREA FIGURE – ULTRA ENLARGED (TIME‑SERIES ORDER: ALASKA, PACIFIC, GULF, ATLANTIC)")
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
    for coast in COAST_ORDER:
        print(f"  {coast}: {coast_counts.get(coast, 0)}")
    
    eco_counts_per_coast = get_ecosystem_counts(df_sites)
    
    print("\nCreating figure ...")
    fig = create_figure(df_monthly, df_sites, coast_counts, eco_counts_per_coast)
    
    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, bbox_inches='tight', facecolor='white')
    print(f"\n✅ Saved PNG:  {OUTPUT_PNG}")
    print(f"✅ Saved PDF:  {OUTPUT_PDF}")
    
    plt.show()

if __name__ == "__main__":
    main()