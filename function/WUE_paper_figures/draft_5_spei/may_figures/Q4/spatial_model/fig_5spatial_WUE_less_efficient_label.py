"""
FIGURE 1 — SPATIAL FREQUENCY OF WUE_T CLASSIFICATION (2000-2025)
Single panel with smooth green-to-grey gradient for Less Efficient Response.
Using pcolormesh for pixel-perfect display of gridded NetCDF data.
UPDATED: WUE_T only, coastal regions from Figure 2
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION - WUE_T ONLY
# ============================================================================

# INPUT: WUE_T aggregate file
INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_time_aggregated\agg_FULLPERIOD_WUE_tra.nc"

# Output file paths
OUTPUT_PNG = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure1_WUE_FullPeriod_SpatialFrequency.png"
OUTPUT_PDF = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure1_WUE_FullPeriod_SpatialFrequency.pdf"

# Variables to plot - ONLY ONE PANEL
VARS_TO_PLOT = ['frac_less']

# ============================================================================
# CUSTOM COLOR PALETTE - SMOOTH GREEN TO GREY GRADIENT
# ============================================================================

green_to_grey_colors = np.array([
    [34/255, 139/255, 34/255],    # 0. Forest Green
    [0.20, 0.60, 0.20],           # 1. Dark Green
    [0.35, 0.68, 0.35],           # 2. Medium Green
    [0.55, 0.75, 0.55],           # 3. Light Green
    [0.70, 0.78, 0.70],           # 4. Very Light Green
    [0.75, 0.75, 0.75],           # 5. Light Grey
    [0.55, 0.55, 0.55],           # 6. Medium Grey
    [0.30, 0.30, 0.30]            # 7. Dark Grey
])

COLORBAR_LABELS = {
    'frac_less': 'Frequency of WUE$_T$ Decline'
}

# Map extents
CONUS_EXTENT = [-126, -66, 24, 50]
ALASKA_EXTENT = [-165, -130, 52, 72]

# ============================================================================
# COASTAL REGIONS - MOVED POSITIONS, DISTINCT COLORS, INCREASED FONT SIZE
# ============================================================================

# Colors for each region (distinct, NOT green/grey)
# Using vibrant colors that stand out from the map
COASTAL_REGIONS = {
    'Alaska': {
        'label_pos': [-155, 60],           # Position in Alaska inset
        'fontsize': 18,                    # Increased font size
        'color': '#1ABC9C',                # Deep Pink/Maroon (kept from Figure 2)
        'is_inset': True
    },
    'Pacific': {
        'label_pos': [-122, 44],            # Moved slightly east (was -124)
        'fontsize': 18,                    # Increased font size
        'color': '#E67E22'                 # Bright Orange (distinct, not green/grey)
    },
    'Gulf': {
        'label_pos': [-94, 26.5],           # Same position
        'fontsize': 18,                    # Increased font size
        'color': '#F1C40F'                   # Purple (distinct)
    },
    'Southeast Atlantic': {
        'label_pos': [-76, 30.5],           # Moved east (was -78)
        'fontsize': 18,                    # Increased font size
        'color':  '#9B59B6'                 # Bright Red (distinct)
    },
    'Atlantic North': {
        'label_pos': [-69, 40],             # Moved east (was -72.5)
        'fontsize': 18,                    # Increased font size
        'color': '#3498DB'                 # Bright Blue (distinct)
    }
}

# Alaska inset position - MOVED TO THE RIGHT
# [left, bottom, width, height] in figure coordinates
# Original was [0.01, 0.70, 0.25, 0.25]
# Moved to right: [0.14, 0.70, 0.25, 0.25] (further right to uncover Pacific label)
ALASKA_INSET_POSITION = [0.14, 0.70, 0.25, 0.25]  # Moved right by 0.13

# Alaska label position in inset
ALASKA_INSET_LABEL_POS = [0.5, 0.99]  # Keep same

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def shift_lon_0_360_to_minus180_180(ds, lon_name="lon"):
    lon = ds[lon_name].values
    if np.nanmax(lon) > 180:
        new_lon = ((lon + 180) % 360) - 180
        ds = ds.assign_coords({lon_name: new_lon}).sortby(lon_name)
    return ds

def calculate_quantile_bins(data_vars):
    print("\n📊 Calculating quantile-based bins...")

    quantiles = [0.00, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 1.00]

    bin_edges = {}
    for var in VARS_TO_PLOT:
        valid = data_vars[var][np.isfinite(data_vars[var])]
        edges = np.quantile(valid, quantiles)
        
        for i in range(1, len(edges)):
            if edges[i] <= edges[i-1]:
                edges[i] = edges[i-1] + 1e-10
        
        bin_edges[var] = edges
        edges_print = np.round(edges, 3)
        print(f"  {var} edges: {[f'{e:.3f}' for e in edges_print]}")

    colormaps = {
        'frac_less': ListedColormap(green_to_grey_colors, name='GreenToGrey')
    }

    norms = {
        var: BoundaryNorm(bin_edges[var], colormaps[var].N, extend='neither')
        for var in VARS_TO_PLOT
    }

    return bin_edges, colormaps, norms

def add_coastal_region_labels(ax, extent):
    """Add coastal region labels with distinct colors."""
    for region_name, region_info in COASTAL_REGIONS.items():
        if region_name == 'Alaska':
            continue
            
        label_pos = region_info['label_pos']
        lon, lat = label_pos
        
        if extent[0] <= lon <= extent[1] and extent[2] <= lat <= extent[3]:
            ax.text(lon, lat, region_name,
                   fontsize=region_info['fontsize'],
                   fontweight='bold',
                   ha='center',
                   va='center',
                   color=region_info['color'],
                   bbox=dict(boxstyle='round,pad=0.4',
                           facecolor='white',
                           alpha=0.9,
                           edgecolor=region_info['color'],
                           linewidth=2),
                   zorder=4,
                   transform=ccrs.PlateCarree())
            
            print(f"  Label: {region_name:20} Position: [{lon:6.1f}, {lat:5.1f}] Color: {region_info['color']}")

def add_alaska_inset(ax_main, alaska_data, lons_alaska, lats_alaska, 
                     colormap, norm):
    """Add Alaska inset - MOVED TO RIGHT to uncover Pacific label."""
    # Use the new position for Alaska inset
    ax_inset = ax_main.inset_axes(ALASKA_INSET_POSITION,
                                   projection=ccrs.PlateCarree())
    
    if len(alaska_data.shape) == 1:
        unique_lons = np.unique(lons_alaska)
        unique_lats = np.unique(lats_alaska)
        lon_grid, lat_grid = np.meshgrid(unique_lons, unique_lats)
        
        data_2d = np.empty((len(unique_lats), len(unique_lons)))
        data_2d[:] = np.nan
        
        for i, lat in enumerate(unique_lats):
            for j, lon in enumerate(unique_lons):
                idx = np.where((lons_alaska == lon) & (lats_alaska == lat))[0]
                if len(idx) > 0:
                    data_2d[i, j] = alaska_data[idx[0]]
        
        sc = ax_inset.pcolormesh(
            lon_grid, lat_grid, data_2d,
            cmap=colormap,
            norm=norm,
            transform=ccrs.PlateCarree(),
            zorder=3
        )
    else:
        sc = ax_inset.pcolormesh(
            lons_alaska, lats_alaska, alaska_data,
            cmap=colormap,
            norm=norm,
            transform=ccrs.PlateCarree(),
            zorder=3
        )
    
    ax_inset.set_extent(ALASKA_EXTENT, crs=ccrs.PlateCarree())
    ax_inset.coastlines(resolution='50m', linewidth=1.8, color='darkgray', zorder=2)
    
    # Alaska label with matching color
    ax_inset.text(ALASKA_INSET_LABEL_POS[0], ALASKA_INSET_LABEL_POS[1], 'Alaska',
                 transform=ax_inset.transAxes,
                 fontsize=18,
                 fontweight='bold',
                 ha='center',
                 va='top',
                 color=COASTAL_REGIONS['Alaska']['color'],
                 bbox=dict(boxstyle='round,pad=0.4',
                         facecolor='white',
                         alpha=0.97,
                         edgecolor=COASTAL_REGIONS['Alaska']['color'],
                         linewidth=2),
                 zorder=5)
    
    for spine in ax_inset.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2.5)
    
    return sc

def add_light_ocean(ax):
    ax.add_feature(cfeature.OCEAN, facecolor='#1E88E5', alpha=0.25, zorder=0)
    return ax

def load_and_validate_data(filepath):
    print(f"📂 Loading data from: {filepath}")
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"❌ File not found: {filepath}")
    
    ds = xr.open_dataset(filepath)
    ds = shift_lon_0_360_to_minus180_180(ds, lon_name="lon")
    
    data_vars = {}
    for var_name in VARS_TO_PLOT:
        if var_name in ds.data_vars:
            data_vars[var_name] = ds[var_name].values
            print(f"  ✓ Loaded {var_name}")
        else:
            raise ValueError(f"Variable '{var_name}' not found")
    
    lats = ds['lat'].values if 'lat' in ds.coords else None
    lons = ds['lon'].values if 'lon' in ds.coords else None
    
    if lats is None or lons is None:
        for coord in ds.coords:
            if 'lat' in coord.lower():
                lats = ds[coord].values
            if 'lon' in coord.lower():
                lons = ds[coord].values
    
    if lats.ndim == 1 and lons.ndim == 1:
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        print(f"  Created 2D grid: {lon_grid.shape}")
    else:
        lon_grid, lat_grid = lons, lats
        print(f"  Using existing 2D grid: {lon_grid.shape}")
    
    ds.close()
    
    return data_vars, lat_grid, lon_grid

# ============================================================================
# FIGURE CREATION
# ============================================================================

def create_single_panel_figure(data_vars, lats_2d, lons_2d, bin_edges, colormaps, norms):
    
    FIG_SIZE = (22, 13)
    
    fig, ax = plt.subplots(1, 1, figsize=FIG_SIZE,
                          subplot_kw={'projection': ccrs.PlateCarree()})
    
    fig.suptitle('Spatial Frequency of Declining WUE$_T$ (2000-2025)',
                fontsize=22,
                fontweight='bold',
                y=0.85)
    
    var_name = VARS_TO_PLOT[0]
    data_2d = data_vars[var_name]
    
    sc_main = ax.pcolormesh(
        lons_2d, lats_2d, data_2d,
        cmap=colormaps[var_name],
        norm=norms[var_name],
        transform=ccrs.PlateCarree(),
        zorder=3
    )
    
    ax.set_extent(CONUS_EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m', linewidth=1.5, color='#2f4f4f', zorder=2)
    add_light_ocean(ax)
    
    print("\n📍 Coastal region labels:")
    add_coastal_region_labels(ax, CONUS_EXTENT)
    
    ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', alpha=0.05, zorder=0)
    
    # Add Alaska inset
    alaska_mask = (lons_2d >= ALASKA_EXTENT[0]) & (lons_2d <= ALASKA_EXTENT[1]) & \
                  (lats_2d >= ALASKA_EXTENT[2]) & (lats_2d <= ALASKA_EXTENT[3])
    
    if np.any(alaska_mask):
        alaska_lons = lons_2d[alaska_mask]
        alaska_lats = lats_2d[alaska_mask]
        alaska_data = data_2d[alaska_mask]
        
        sc_inset = add_alaska_inset(
            ax, alaska_data, alaska_lons, alaska_lats,
            colormaps[var_name], norms[var_name]
        )
    
    cbar = plt.colorbar(sc_main, ax=ax, orientation='vertical',
                       fraction=0.04, pad=0.02, shrink=0.8)
    
    cbar.set_label(COLORBAR_LABELS[var_name], 
                  fontsize=20,
                  fontweight='bold',
                  labelpad=15)
    
    edges = bin_edges[var_name]
    tick_positions = edges[:-1] + np.diff(edges)/2
    tick_labels = [f'{pos:.2f}' for pos in tick_positions]
    
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(labelsize=20)
    
    plt.subplots_adjust(left=0.05, right=0.90, top=0.92, bottom=0.05)
    
    return fig, ax

def save_figures(fig, png_path, pdf_path):
    png_dir = os.path.dirname(png_path)
    pdf_dir = os.path.dirname(pdf_path)
    
    for dir_path in [png_dir, pdf_dir]:
        if dir_path and not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print(f"📁 Created directory: {dir_path}")
    
    fig.savefig(png_path, dpi=400, bbox_inches='tight', facecolor='white')
    print(f"✅ Saved PNG (400 DPI): {png_path}")
    
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"✅ Saved PDF: {pdf_path}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("=" * 80)
    print("FIGURE 1: SINGLE PANEL - Frequency of WUE$_T$ decline")
    print("=" * 80)
    print("🎨 COLOR SCHEME: Smooth Forest Green to Dark Grey Gradient")
    print("📐 DISPLAY: Using pcolormesh() for pixel-perfect gridded data")
    print("📍 COASTAL REGIONS (distinct colors, increased font size):")
    print("  • Alaska: #F06292 (Deep Pink) - moved inset to reveal Pacific label")
    print("  • Pacific: #FF6B35 (Bright Orange) - moved east")
    print("  • Gulf: #9B59B6 (Purple)")
    print("  • Southeast Atlantic: #E74C3C (Bright Red) - moved east")
    print("  • Atlantic North: #3498DB (Bright Blue) - moved east")
    print("📝 Note: Labels have colored borders, increased font size (18pt)")
    print("🌊 Ocean colored, lakes remain transparent")
    print()
    
    try:
        data_vars, lats_2d, lons_2d = load_and_validate_data(INPUT_FILE)
        bin_edges, colormaps, norms = calculate_quantile_bins(data_vars)
        
        print("\n🎨 Creating figure...")
        fig, ax = create_single_panel_figure(data_vars, lats_2d, lons_2d, 
                                            bin_edges, colormaps, norms)
        
        print("\n💾 Saving outputs...")
        save_figures(fig, OUTPUT_PNG, OUTPUT_PDF)
        
        var_name = VARS_TO_PLOT[0]
        data_2d = data_vars[var_name]
        valid_data = data_2d[np.isfinite(data_2d)]
        
        print("\n📊 DATA DISTRIBUTION (WUE$_T$):")
        print(f"  Grid shape: {data_2d.shape}")
        print(f"  Valid pixels: {len(valid_data):,}")
        print(f"  Data range: [{valid_data.min():.3f}, {valid_data.max():.3f}]")
        print(f"  Mean: {valid_data.mean():.3f}, Std: {valid_data.std():.3f}")
        
        print("\n✅ Figure creation completed!")
        print(f"\n📋 Files saved:")
        print(f"  1. {OUTPUT_PNG}")
        print(f"  2. {OUTPUT_PDF}")
        print(f"\n💡 Key features:")
        print(f"   • WUE$_T$ (transpiration) only")
        print(f"   • Alaska inset moved to right (position: {ALASKA_INSET_POSITION})")
        print(f"   • Pacific label now visible")
        print(f"   • All coastal labels have distinct colors (no green/grey)")
        print(f"   • Atlantic labels moved further east")
        print(f"   • Font size increased to 18pt")
        
        plt.show()
        
    except FileNotFoundError as e:
        print(f"❌ FILE NOT FOUND ERROR: {e}")
        print(f"\n💡 The file 'agg_FULLPERIOD_WUE_tra.nc' may not exist yet.")
        print(f"   Please run the Concise Spatial Analyzer first to create it.")
        
    except Exception as e:
        print(f"❌ ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()