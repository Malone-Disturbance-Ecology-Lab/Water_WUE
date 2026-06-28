# -*- coding: utf-8 -*-
"""
STUDY AREA MAP - 57 Coastal Sites by Ecosystem Type
UPDATED: Horizontal bars and global colorbars use MONTHLY min/max (not site medians)
- Dots: Site-level MEDIAN (Trans_ratio_median)
- Region bars: MONTHLY min/max across all sites in region
- Global colorbars: MONTHLY min/max across ALL sites
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# FILE PATHS
# ============================================================================

MONTHLY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
SUMMARY_DATA_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\site_summary_with_SPEI48.csv"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)


INPUT_FILE = SUMMARY_DATA_PATH  # This defines INPUT_FILE for the load_site_data function


OUTPUT_PNG = os.path.join(OUTPUT_DIR, "Figure_StudyArea_57sites_TET_gradient_with_SPEI48.png")
OUTPUT_PDF = os.path.join(OUTPUT_DIR, "Figure_StudyArea_57sites_TET_gradient_with_SPEI48.pdf")

# ============================================================================
# ⚙️ EASY-TO-EDIT PARAMETERS
# ============================================================================

TET_BAR_ABOVE_LABEL = 1.8
SPEI_BAR_BELOW_LABEL = 2
FONT_SIZE_BAR_VALUES = 20
FONT_SIZE_CBAR_TICKS = 20
FONT_SIZE_CBAR_LABEL = 20
PACIFIC_BAR_LENGTH = 5.0

FIGURE_WIDTH = 22
FIGURE_HEIGHT = 14
CONUS_EXTENT = [-130, -62, 22, 52]
ALASKA_EXTENT = [-170, -128, 49, 74]

MARKER_SIZE = 150
BAR_HEIGHT = 8

# ============================================================================
# COLOR GRADIENTS
# ============================================================================
TET_CMAP = LinearSegmentedColormap.from_list(
    "TET_gradient",
    [
        (0.0, "#D7A900"),
        (0.33, "#9B8A71"),
        (0.66, "#4169E1"),
        (1.0, "#9B00B4")
    ]
)

SPEI_CMAP = LinearSegmentedColormap.from_list(
    "balanced_spei_gradient",
    [
        "#9E0142",
        "#D53E4F",
        "#F46D43",
        "#F7F7F7",
        "#ABD9E9",
        "#74ADD1",
        "#2B83BA"
    ]
)

# ============================================================================
# ECOSYSTEM MARKERS
# ============================================================================
ECOSYSTEM_MARKERS = {
    'Upland': 'o',
    'Freshwater': 's',
    'Brackish': 'D',
    'Saline': '^'
}

LEGEND_ORDER = ['Upland', 'Freshwater', 'Brackish', 'Saline']
MARKER_EDGE_WIDTH = 0.5

# ============================================================================
# COASTAL REGIONS
# ============================================================================
COASTAL_REGIONS = {
    'Alaska': {
        'label_pos': [-148, 60],
        'tet_bar_offset': TET_BAR_ABOVE_LABEL,
        'spei_bar_offset': -SPEI_BAR_BELOW_LABEL,
        'bar_length': 6.0,
        'is_inset': True
    },
    'Pacific': {
        'label_pos': [-124, 42],
        'tet_bar_offset': TET_BAR_ABOVE_LABEL,
        'spei_bar_offset': -SPEI_BAR_BELOW_LABEL,
        'bar_length': PACIFIC_BAR_LENGTH,
        'is_inset': False
    },
    'Gulf Coast': {
        'label_pos': [-88, 27],
        'tet_bar_offset': TET_BAR_ABOVE_LABEL-0.1,
        'spei_bar_offset': -SPEI_BAR_BELOW_LABEL-0.3,
        'bar_length': 6.5,
        'is_inset': False
    },
    'Southeast Atlantic': {
        'label_pos': [-75, 30.5],
    'bar_lon_offset': 0.8,
        'tet_bar_offset': TET_BAR_ABOVE_LABEL+0.1,
        'spei_bar_offset': -SPEI_BAR_BELOW_LABEL-0.3,
        'bar_length': 6.0,
        'is_inset': False
    },
    'Atlantic North': {
        'label_pos': [-67, 40],
        'bar_lon_offset': -1.50,
        'tet_bar_offset': TET_BAR_ABOVE_LABEL+ 0.2,
        'spei_bar_offset': -SPEI_BAR_BELOW_LABEL- 0.3,
        'bar_length': 6.0,
        'is_inset': False
    }
}

ALASKA_INSET_POSITION = [0.13, 0.56, 0.36, 0.36]
ALASKA_INSET_LABEL_POS = [0.80, 0.75]
LEGEND_BBOX_X = 0.88

# ============================================================================
# DATA LOADING WITH MONTHLY MIN/MAX
# ============================================================================

def load_monthly_data():
    """Load monthly data for min/max calculations"""
    print(f"📂 Loading monthly data...")
    df_monthly = pd.read_csv(MONTHLY_DATA_PATH)
    print(f"  ✓ Loaded {len(df_monthly):,} rows")
    return df_monthly

def load_site_data(filepath, df_monthly):
    """Load site summary data and merge with monthly min/max for regions"""
    print(f"📂 Loading site data...")
    df = pd.read_csv(filepath)
    print(f"  ✓ Loaded {len(df)} sites")
    
    if 'Salinity_Category' in df.columns:
        df['Ecosystem_Type'] = df['Salinity_Category'].astype(str).str.strip()
    
    if 'coast_region' in df.columns:
        df['coast_region'] = df['coast_region'].astype(str).str.strip()
        df['coast_region'] = df['coast_region'].replace('Gulf of Mexico', 'Gulf Coast')
    
    # Get monthly min/max for T:ET across ALL sites (for global colorbar)
    if 'Trans_ratio' in df_monthly.columns:
        global_tet_min = df_monthly['Trans_ratio'].min()
        global_tet_max = df_monthly['Trans_ratio'].max()
        print(f"\n  🌍 Global T:ET range from monthly data: {global_tet_min:.3f} - {global_tet_max:.3f}")
    else:
        global_tet_min = df['Trans_ratio_median'].min()
        global_tet_max = df['Trans_ratio_median'].max()
        print(f"  ⚠️ Using site medians for global T:ET range")
    
    # Get monthly min/max for SPEI across ALL sites (for global colorbar)
    if 'SPEI_6' in df_monthly.columns:
        global_spei_min = df_monthly['SPEI_6'].min()
        global_spei_max = df_monthly['SPEI_6'].max()
        print(f"  🌍 Global SPEI range from monthly data: {global_spei_min:.2f} - {global_spei_max:.2f}")
    else:
        global_spei_min = df['spei_min'].min()
        global_spei_max = df['spei_max'].max()
        print(f"  ⚠️ Using site min/max for global SPEI range")
    
    return df, global_tet_min, global_tet_max, global_spei_min, global_spei_max

def get_coast_ranges_with_monthly_minmax(df, df_monthly):
    """
    Calculate regional T:ET and SPEI ranges using MONTHLY min/max
    across all sites in each region (not site medians)
    """
    coast_ranges = {}
    
    # Map coast region names to match monthly data
    region_mapping = {
        'Alaska': 'AK_coast',
        'Pacific': 'West_Coast',
        'Gulf Coast': 'Gulf_of_America',
        'Southeast Atlantic': 'Southeast_Atlantic',
        'Atlantic North': 'Atlantic_Coast_North'
    }
    
    for region_name, region_info in COASTAL_REGIONS.items():
        # Get sites in this region from summary data
        region_sites = df[df['coast_region'] == region_name]['site_name'].tolist()
        
        if len(region_sites) == 0:
            continue
        
        # Get monthly data for these sites
        region_monthly = df_monthly[df_monthly['site_name'].isin(region_sites)]
        
        # Calculate T:ET range from MONTHLY data (not site medians)
        if 'Trans_ratio' in region_monthly.columns and len(region_monthly) > 0:
            tet_min = region_monthly['Trans_ratio'].min()
            tet_max = region_monthly['Trans_ratio'].max()
        else:
            # Fallback to site medians
            tet_min = df[df['coast_region'] == region_name]['Trans_ratio_median'].min()
            tet_max = df[df['coast_region'] == region_name]['Trans_ratio_median'].max()
        
        # Calculate SPEI range from MONTHLY data
        if 'SPEI_6' in region_monthly.columns and len(region_monthly) > 0:
            spei_min = region_monthly['SPEI_6'].min()
            spei_max = region_monthly['SPEI_6'].max()
        else:
            # Fallback to site min/max
            spei_min = df[df['coast_region'] == region_name]['spei_min'].min()
            spei_max = df[df['coast_region'] == region_name]['spei_max'].max()
        
        coast_ranges[region_name] = {
            'n_sites': len(region_sites),
            'tet_min': tet_min,
            'tet_max': tet_max,
            'spei_min': spei_min,
            'spei_max': spei_max,
        }
        print(f"  {region_name}: T:ET [{tet_min:.3f}, {tet_max:.3f}]")
        print(f"            SPEI [{spei_min:.2f}, {spei_max:.2f}]")
    
    return coast_ranges

# ============================================================================
# COLOR BAR FUNCTIONS
# ============================================================================

def add_tet_horizontal_bar(ax, lon, lat, value_min, value_max, norm, bar_length):
    """Add T:ET bar ABOVE label box using monthly min/max"""
    start_lon = lon - bar_length / 2
    end_lon = lon + bar_length / 2
    
    n_segments = 100
    for seg in range(n_segments):
        t = seg / n_segments
        value = value_min + t * (value_max - value_min)
        bar_color = TET_CMAP(norm(value))
        
        x_start = start_lon + t * bar_length
        x_end = start_lon + (seg + 1) / n_segments * bar_length
        
        ax.hlines(y=lat, xmin=x_start, xmax=x_end,
                  color=bar_color, linewidth=BAR_HEIGHT, zorder=6, alpha=0.95,
                  transform=ccrs.PlateCarree())
    
    min_color = TET_CMAP(norm(value_min))
    max_color = TET_CMAP(norm(value_max))
    
    offset = 0.35 if lon < -120 else 0.25
    
    ax.text(start_lon - offset, lat, f'{value_min:.2f}',
            ha='right', va='center', fontsize=FONT_SIZE_BAR_VALUES,
            fontweight='bold', color=min_color, transform=ccrs.PlateCarree(), zorder=7)
    
    ax.text(end_lon + offset, lat, f'{value_max:.2f}',
            ha='left', va='center', fontsize=FONT_SIZE_BAR_VALUES,
            fontweight='bold', color=max_color, transform=ccrs.PlateCarree(), zorder=7)

def add_spei_horizontal_bar(ax, lon, lat, value_min, value_max, norm, bar_length):
    """Add SPEI bar BELOW label box using monthly min/max"""
    start_lon = lon - bar_length / 2
    end_lon = lon + bar_length / 2
    
    n_segments = 100
    for seg in range(n_segments):
        t = seg / n_segments
        value = value_min + t * (value_max - value_min)
        bar_color = SPEI_CMAP(norm(value))
        
        x_start = start_lon + t * bar_length
        x_end = start_lon + (seg + 1) / n_segments * bar_length
        
        ax.hlines(y=lat, xmin=x_start, xmax=x_end,
                  color=bar_color, linewidth=BAR_HEIGHT, zorder=6, alpha=0.95,
                  transform=ccrs.PlateCarree())
    
    min_color = SPEI_CMAP(norm(value_min))
    max_color = SPEI_CMAP(norm(value_max))
    
    offset = 0.35 if lon < -120 else 0.25
    
    ax.text(start_lon - offset, lat, f'{value_min:.1f}',
            ha='right', va='center', fontsize=FONT_SIZE_BAR_VALUES,
            fontweight='bold', color=min_color, transform=ccrs.PlateCarree(), zorder=7)
    
    ax.text(end_lon + offset, lat, f'{value_max:.1f}',
            ha='left', va='center', fontsize=FONT_SIZE_BAR_VALUES,
            fontweight='bold', color=max_color, transform=ccrs.PlateCarree(), zorder=7)

def add_two_row_coastal_label(ax, region_name, label_pos, n_sites):
    """Add coastal label box"""
    lon, lat = label_pos
    label_text = f'{region_name}\nN = {n_sites}'
    
    ax.text(lon, lat, label_text,
            fontsize=17, fontweight='bold', ha='center', va='center', color='black',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.95,
                     edgecolor='black', linewidth=1.5),
            transform=ccrs.PlateCarree(), zorder=5, linespacing=1.5)

# ============================================================================
# MAPPING FUNCTIONS
# ============================================================================

def add_coastal_labels_with_bars(ax, extent, coast_ranges, tet_norm, spei_norm):
    """Add labels with T:ET bar ABOVE and SPEI bar BELOW using monthly ranges"""
    for region_name, region_info in COASTAL_REGIONS.items():
        if region_name == 'Alaska':
            continue
        
        if region_name not in coast_ranges:
            continue
        
        label_pos = region_info['label_pos']
        lon, lat = label_pos
        
        if extent[0] <= lon <= extent[1] and extent[2] <= lat <= extent[3]:
            ranges = coast_ranges[region_name]
            
            add_two_row_coastal_label(ax, region_name, label_pos, ranges['n_sites'])
            
            # T:ET bar ABOVE (using monthly min/max for region)
            tet_bar_lat = lat + region_info.get('tet_bar_offset', 1.4)
            bar_length = region_info.get('bar_length', 6.0)
            bar_lon = lon + region_info.get('bar_lon_offset', 0)
            add_tet_horizontal_bar(ax, bar_lon, tet_bar_lat, ranges['tet_min'], 
                                   ranges['tet_max'], tet_norm, bar_length)
            
            # SPEI bar BELOW (using monthly min/max for region)
            spei_bar_lat = lat + region_info.get('spei_bar_offset', -1.6)
            add_spei_horizontal_bar(ax, bar_lon, spei_bar_lat, ranges['spei_min'],
                                     ranges['spei_max'], spei_norm, bar_length)
            
            print(f"  {region_name:20} | Label:{lat:.1f} | T:ET bar:{tet_bar_lat:.1f} | SPEI bar:{spei_bar_lat:.1f}")

def plot_sites(ax, df, tet_cmap, tet_norm):
    """Plot sites with color based on site MEDIAN (not monthly)"""
    ecosystem_counts = df['Ecosystem_Type'].value_counts()
    for ecosystem_type, marker in ECOSYSTEM_MARKERS.items():
        ecosystem_sites = df[df['Ecosystem_Type'] == ecosystem_type]
        for _, site in ecosystem_sites.iterrows():
            tet_value = site['Trans_ratio_median']  # Using MEDIAN for dots
            color = tet_cmap(tet_norm(tet_value))
            ax.scatter(site['long'], site['lat'], c=[color], marker=marker,
                      s=MARKER_SIZE, edgecolor=color, linewidth=MARKER_EDGE_WIDTH,
                      zorder=5, transform=ccrs.PlateCarree())
    return ecosystem_counts

def add_alaska_inset(ax_main, df_alaska, tet_cmap, tet_norm, spei_norm, coast_ranges):
    """Alaska inset with bars using monthly min/max"""
    ax_inset = ax_main.inset_axes(ALASKA_INSET_POSITION, projection=ccrs.PlateCarree())
    ax_inset.set_extent(ALASKA_EXTENT, crs=ccrs.PlateCarree())
    ax_inset.coastlines(resolution='50m', linewidth=1.5, color='darkgray', zorder=2)
    ax_inset.add_feature(cfeature.OCEAN, facecolor='#B8D4E8', alpha=0.3, zorder=0)
    ax_inset.add_feature(cfeature.LAND, facecolor='#f0f0f0', alpha=0.3, zorder=0)
    
    # Plot sites using MEDIAN
    for _, site in df_alaska.iterrows():
        ecosystem = site['Ecosystem_Type']
        tet_value = site['Trans_ratio_median']
        marker = ECOSYSTEM_MARKERS.get(ecosystem, 'o')
        color = tet_cmap(tet_norm(tet_value))
        ax_inset.scatter(site['long'], site['lat'], c=[color], marker=marker,
                        s=MARKER_SIZE * 0.7, edgecolor=color, linewidth=MARKER_EDGE_WIDTH,
                        zorder=5, transform=ccrs.PlateCarree())
    
    if 'Alaska' in coast_ranges:
        ranges = coast_ranges['Alaska']
        
        label_x = ALASKA_INSET_LABEL_POS[0]
        label_y = ALASKA_INSET_LABEL_POS[1]
        label_text = f'Alaska\nN = {ranges["n_sites"]}'
        
        ax_inset.text(label_x, label_y, label_text, transform=ax_inset.transAxes,
                     fontsize=17, fontweight='bold', ha='center', va='center', color='black',
                     bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.97,
                              edgecolor='black', linewidth=1.5), zorder=5, linespacing=1.5)
        
        # T:ET bar ABOVE (using monthly min/max)
        tet_bar_y = label_y + 0.18
        bar_length_axes = 0.65
        bar_x_center = label_x - 0.30
        
        for seg in range(100):
            t = seg / 100
            tet_value = ranges['tet_min'] + t * (ranges['tet_max'] - ranges['tet_min'])
            bar_color = TET_CMAP(tet_norm(tet_value))
            x_start = bar_x_center - bar_length_axes/2 + t * bar_length_axes
            x_end = bar_x_center - bar_length_axes/2 + (seg + 1) / 100 * bar_length_axes
            ax_inset.hlines(y=tet_bar_y, xmin=x_start, xmax=x_end, color=bar_color,
                          linewidth=BAR_HEIGHT * 0.8, zorder=6, alpha=0.95,
                          transform=ax_inset.transAxes)
        
        min_color = TET_CMAP(tet_norm(ranges['tet_min']))
        max_color = TET_CMAP(tet_norm(ranges['tet_max']))
        ax_inset.text(bar_x_center - bar_length_axes/2 - 0.001, tet_bar_y, f'{ranges["tet_min"]:.2f}',
                     ha='right', va='center', fontsize=FONT_SIZE_BAR_VALUES,
                     fontweight='bold', color=min_color, transform=ax_inset.transAxes, zorder=7)
        ax_inset.text(bar_x_center + bar_length_axes/2 + 0.001, tet_bar_y, f'{ranges["tet_max"]:.2f}',
                     ha='left', va='center', fontsize=FONT_SIZE_BAR_VALUES,
                     fontweight='bold', color=max_color, transform=ax_inset.transAxes, zorder=7)
        
        # SPEI bar BELOW (using monthly min/max)
        spei_bar_y = label_y - 0.22
        
        for seg in range(100):
            t = seg / 100
            spei_value = ranges['spei_min'] + t * (ranges['spei_max'] - ranges['spei_min'])
            bar_color = SPEI_CMAP(spei_norm(spei_value))
            x_start = bar_x_center - bar_length_axes/2 + t * bar_length_axes
            x_end = bar_x_center - bar_length_axes/2 + (seg + 1) / 100 * bar_length_axes
            ax_inset.hlines(y=spei_bar_y, xmin=x_start, xmax=x_end, color=bar_color,
                          linewidth=BAR_HEIGHT * 0.8, zorder=6, alpha=0.95,
                          transform=ax_inset.transAxes)
        
        min_color_spei = SPEI_CMAP(spei_norm(ranges['spei_min']))
        max_color_spei = SPEI_CMAP(spei_norm(ranges['spei_max']))
        ax_inset.text(bar_x_center - bar_length_axes/2 - 0.001, spei_bar_y, f'{ranges["spei_min"]:.1f}',
                     ha='right', va='center', fontsize=FONT_SIZE_BAR_VALUES,
                     fontweight='bold', color=min_color_spei, transform=ax_inset.transAxes, zorder=7)
        ax_inset.text(bar_x_center + bar_length_axes/2 + 0.001, spei_bar_y, f'{ranges["spei_max"]:.1f}',
                     ha='left', va='center', fontsize=FONT_SIZE_BAR_VALUES,
                     fontweight='bold', color=max_color_spei, transform=ax_inset.transAxes, zorder=7)
    
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)
    
    return ax_inset

def add_ecosystem_legend(ax, ecosystem_counts):
    legend_elements = []
    for ecosystem_type in LEGEND_ORDER:
        if ecosystem_type in ECOSYSTEM_MARKERS:
            marker = ECOSYSTEM_MARKERS[ecosystem_type]
            count = ecosystem_counts.get(ecosystem_type, 0)
            legend_elements.append(
                plt.Line2D([0], [0], marker=marker, color='none',
                          markerfacecolor='black', markeredgecolor='black',
                          markersize=12, linestyle='none',
                          label=f'{ecosystem_type} (N={count})')
            )
    legend = ax.legend(handles=legend_elements, loc='upper right',
                      bbox_to_anchor=(LEGEND_BBOX_X, 0.98), fontsize=16, framealpha=0.95,
                      title='Ecosystem Type', title_fontsize=18)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.5)
    return legend

# ============================================================================
# MAIN FIGURE
# ============================================================================

def create_study_area_map(df, df_monthly, coast_ranges, global_tet_min, global_tet_max, 
                          global_spei_min, global_spei_max):
    
    conus_sites = df[(df['long'] >= CONUS_EXTENT[0]) & (df['long'] <= CONUS_EXTENT[1]) & 
                     (df['lat'] >= CONUS_EXTENT[2]) & (df['lat'] <= CONUS_EXTENT[3])].copy()
    alaska_sites = df[(df['long'] >= ALASKA_EXTENT[0]) & (df['long'] <= ALASKA_EXTENT[1]) & 
                      (df['lat'] >= ALASKA_EXTENT[2]) & (df['lat'] <= ALASKA_EXTENT[3])].copy()
    
    # Use GLOBAL monthly min/max for colorbar normalization
    tet_norm = Normalize(vmin=global_tet_min, vmax=global_tet_max)
    spei_norm = Normalize(vmin=global_spei_min, vmax=global_spei_max)
    
    print(f"\n  Global T:ET range (monthly): {global_tet_min:.3f} - {global_tet_max:.3f}")
    print(f"  Global SPEI range (monthly): {global_spei_min:.2f} - {global_spei_max:.2f}")
    
    fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT),
                          subplot_kw={'projection': ccrs.PlateCarree()})
    
    ax.set_extent(CONUS_EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m', linewidth=1.2, color='#666666', zorder=2)
    ax.add_feature(cfeature.OCEAN, facecolor='#E0F0FA', alpha=0.5, zorder=0)
    ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', alpha=0.05, zorder=0)
    ax.add_feature(cfeature.STATES, linewidth=0.8, edgecolor='#999999', alpha=0.4, zorder=1)
    
    print("\n📍 Adding coastal labels with bars (using monthly min/max):")
    add_coastal_labels_with_bars(ax, CONUS_EXTENT, coast_ranges, tet_norm, spei_norm)
    
    ecosystem_counts = plot_sites(ax, conus_sites, TET_CMAP, tet_norm)
    
    if len(alaska_sites) > 0:
        add_alaska_inset(ax, alaska_sites, TET_CMAP, tet_norm, spei_norm, coast_ranges)
        for _, site in alaska_sites.iterrows():
            eco = site['Ecosystem_Type']
            ecosystem_counts[eco] = ecosystem_counts.get(eco, 0) + 1
    
    add_ecosystem_legend(ax, ecosystem_counts)
    
    # GLOBAL COLORBARS (using monthly min/max)
    sm_tet = ScalarMappable(norm=tet_norm, cmap=TET_CMAP)
    sm_tet.set_array(np.linspace(global_tet_min, global_tet_max, 100))
    cbar_tet = plt.colorbar(sm_tet, ax=ax, orientation='vertical',
                           fraction=0.025, pad=0.06, shrink=0.55)
    cbar_tet.set_label('T:ET Ratio', fontsize=FONT_SIZE_CBAR_LABEL, fontweight='bold', labelpad=10)
    cbar_tet.ax.tick_params(labelsize=FONT_SIZE_CBAR_TICKS)
    
    sm_spei = ScalarMappable(norm=spei_norm, cmap=SPEI_CMAP)
    sm_spei.set_array(np.linspace(global_spei_min, global_spei_max, 100))
    cbar_spei = plt.colorbar(sm_spei, ax=ax, orientation='vertical',
                            fraction=0.025, pad=0.01, shrink=0.55)
    cbar_spei.set_label('SPEI-48', fontsize=FONT_SIZE_CBAR_LABEL, fontweight='bold', labelpad=10)
    cbar_spei.ax.tick_params(labelsize=FONT_SIZE_CBAR_TICKS)
    pos_tet = cbar_tet.ax.get_position()
    cbar_spei.ax.set_position([pos_tet.x0 + 0.06, pos_tet.y0, pos_tet.width, pos_tet.height])
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.3, linestyle='--', zorder=1)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 16}
    gl.ylabel_style = {'size': 16}
    
    plt.subplots_adjust(right=0.92, left=0.05, top=0.95, bottom=0.05)
    
    return fig, ax

def save_figures(fig, png_path, pdf_path):
    fig.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"✅ Saved PNG: {png_path}")
    fig.savefig(pdf_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"✅ Saved PDF: {pdf_path}")

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("="*80)
    print("STUDY AREA MAP - USING MONTHLY MIN/MAX FOR BARS")
    print("  - Dots: Site MEDIAN (Trans_ratio_median)")
    print("  - Region bars: MONTHLY min/max across sites in region")
    print("  - Global colorbars: MONTHLY min/max across ALL sites")
    print("="*80)
    
    try:
        # Load monthly data first (for min/max calculations)
        df_monthly = load_monthly_data()
        
        # Load site data and get global monthly ranges
        df, global_tet_min, global_tet_max, global_spei_min, global_spei_max = load_site_data(INPUT_FILE, df_monthly)
        
        # Filter to valid sites
        df = df.dropna(subset=['Trans_ratio_median', 'spei_min', 'spei_max'])
        print(f"  Filtered to {len(df)} sites")
        
        # Get coastal region ranges using MONTHLY min/max
        print(f"\n📊 Coastal region ranges (from monthly data):")
        coast_ranges = get_coast_ranges_with_monthly_minmax(df, df_monthly)
        
        print("\n🎨 Creating map...")
        fig, ax = create_study_area_map(df, df_monthly, coast_ranges, 
                                        global_tet_min, global_tet_max,
                                        global_spei_min, global_spei_max)
        
        print("\n💾 Saving...")
        save_figures(fig, OUTPUT_PNG, OUTPUT_PDF)
        
        print("\n✅ COMPLETED!")
        plt.show()
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()