# -*- coding: utf-8 -*-
"""
FIGURE 3 — PRE vs POST vs Δ (POST−PRE) FOR frac_less
Efficient workflow matching Figure 1 style exactly.
UPDATED: WUE_T only, 5 coastal regions from Figure 2, distinct colors
Coastal labels ONLY in Panel (a) - not repeated in (b) and (c)
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
# CONFIGURATION - UPDATED FOR WUE_T
# ============================================================================

# File paths - UPDATED to WUE_T files
PRE_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_pre_post_breakpoint\agg_PRE_WUE_tra.nc"
POST_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\rasters_pre_post_breakpoint\agg_POST_WUE_tra.nc"

# Output paths
OUTPUT_PNG = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure3_WUE_T_PRE_POST_DELTA_frac_less.png"
OUTPUT_PDF = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure3_WUE_T_PRE_POST_DELTA_frac_less.pdf"

VAR_NAME = "frac_less"

# ============================================================================
# STYLING (IDENTICAL TO UPDATED FIGURE 1)
# ============================================================================

# Color palette from Figure 1 (green to grey)
GREEN_TO_GREY = np.array([
    [34/255, 139/255, 34/255],    # 0. Forest Green (lowest values)
    [0.20, 0.60, 0.20],           # 1. Dark Green
    [0.35, 0.68, 0.35],           # 2. Medium Green
    [0.55, 0.75, 0.55],           # 3. Light Green
    [0.70, 0.78, 0.70],           # 4. Very Light Green
    [0.75, 0.75, 0.75],           # 5. Light Grey (transition)
    [0.55, 0.55, 0.55],           # 6. Medium Grey
    [0.30, 0.30, 0.30]            # 7. Dark Grey (highest values)
])

# Map extents
CONUS_EXTENT = [-126, -66, 24, 50]
ALASKA_EXTENT = [-165, -130, 52, 72]

# Alaska inset position - MOVED TO THE RIGHT (same as Figure 1)
ALASKA_INSET_POSITION = [0.14, 0.70, 0.25, 0.25]

# ============================================================================
# COASTAL REGIONS - DISTINCT COLORS (UPDATED WITH MOVED LABELS)
# ============================================================================

COASTAL_REGIONS = {
    'Alaska': {
        'label_pos': [-155, 60],
        'fontsize': 16,
        'color': '#1ABC9C',                # Turquoise/Teal
        'is_inset': True
    },
    'Pacific': {
        'label_pos': [-126, 46],            # MOVED: was [-122, 44] - now further west/north
        'fontsize': 16,
        'color': '#E67E22',                 # Orange
        'is_inset': False
    },
    'Gulf': {
        'label_pos': [-92, 28],             # MOVED: was [-94, 26.5] - now further north/east
        'fontsize': 16,
        'color': '#F1C40F',                 # Yellow/Gold
        'is_inset': False
    },
    'Southeast Atlantic': {
        'label_pos': [-77, 31],             # MOVED: was [-76, 30.5] - now further west/north
        'fontsize': 16,
        'color': '#9B59B6',                 # Purple
        'is_inset': False
    },
    'Atlantic North': {
        'label_pos': [-67, 43],             # MOVED: was [-69, 40] - now further east/north
        'fontsize': 16,
        'color': '#3498DB',                 # Bright Blue
        'is_inset': False
    }
}

# Quantiles for PRE/POST
QUANTILES = [0.00, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 1.00]

# ============================================================================
# CORE DATA FUNCTIONS
# ============================================================================

def load_data_and_grid(filepath, var_name):
    """Load NetCDF data and create 2D coordinate grids."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    ds = xr.open_dataset(filepath)
    
    if np.nanmax(ds.lon.values) > 180:
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180).sortby('lon')
    
    data = ds[var_name].values
    lats = ds.lat.values
    lons = ds.lon.values
    
    if lats.ndim == 1 and lons.ndim == 1:
        lon_grid, lat_grid = np.meshgrid(lons, lats)
    else:
        lon_grid, lat_grid = lons, lats
    
    ds.close()
    return data, lat_grid, lon_grid

def validate_data_masks(pre_data, post_data):
    """Validate that PRE and POST have similar data masks."""
    print("\n🔍 VALIDATING DATA MASKS:")
    pre_valid = np.sum(np.isfinite(pre_data))
    post_valid = np.sum(np.isfinite(post_data))
    both_valid = np.sum(np.isfinite(pre_data) & np.isfinite(post_data))
    
    print(f"  PRE valid pixels:  {pre_valid:,}")
    print(f"  POST valid pixels: {post_valid:,}")
    print(f"  BOTH valid pixels: {both_valid:,}")
    print(f"  Mask agreement: {100*both_valid/min(pre_valid, post_valid):.1f}%")
    
    return both_valid

def calculate_shared_bins(pre_data, post_data):
    """Calculate shared quantile bins from combined PRE+POST data."""
    pre_valid = pre_data[np.isfinite(pre_data)]
    post_valid = post_data[np.isfinite(post_data)]
    combined = np.concatenate([pre_valid, post_valid])
    
    edges = np.quantile(combined, QUANTILES)
    
    for i in range(1, len(edges)):
        if edges[i] <= edges[i-1]:
            edges[i] = edges[i-1] + 1e-10
    
    return edges

def calculate_delta_bins(delta_data, mask):
    """Calculate symmetric bins for delta map centered at 0."""
    delta_valid = delta_data[mask]
    max_abs = np.nanpercentile(np.abs(delta_valid), 99)
    edges = np.linspace(-max_abs, max_abs, 9)
    edges[4] = 0.0
    return edges

def create_delta_colormap():
    """Create a discrete diverging colormap for delta."""
    base_cmap = plt.cm.RdBu_r
    colors = [base_cmap(i / 9.0) for i in range(10)]
    return ListedColormap(colors, name='DiscreteRdBu_r')

# ============================================================================
# LABELING FUNCTIONS (UPDATED FOR TWO-ROW LABELS)
# ============================================================================

def add_coastal_region_labels(ax, extent):
    """Add coastal region labels with distinct colors and two-row names."""
    for region_name, region_info in COASTAL_REGIONS.items():
        if region_info.get('is_inset', False):
            continue
            
        label_pos = region_info['label_pos']
        lon, lat = label_pos
        
        if extent[0] <= lon <= extent[1] and extent[2] <= lat <= extent[3]:
            # Create two-row labels for multi-word names
            if region_name == 'Southeast Atlantic':
                display_name = 'Southeast\nAtlantic'  # Two rows
            elif region_name == 'Atlantic North':
                display_name = 'Atlantic\nNorth'      # Two rows
            else:
                display_name = region_name            # One row (Alaska, Pacific, Gulf)
            
            ax.text(lon, lat, display_name,
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

def add_alaska_inset(ax, data, lons, lats, cmap, norm):
    """Add Alaska inset with Alaska label."""
    ax_inset = ax.inset_axes(ALASKA_INSET_POSITION,
                            projection=ccrs.PlateCarree())
    
    alaska_mask = ((lons >= ALASKA_EXTENT[0]) & (lons <= ALASKA_EXTENT[1]) &
                   (lats >= ALASKA_EXTENT[2]) & (lats <= ALASKA_EXTENT[3]))
    alaska_data = np.where(alaska_mask, data, np.nan)
    
    ax_inset.pcolormesh(lons, lats, alaska_data,
                       cmap=cmap, norm=norm,
                       transform=ccrs.PlateCarree(), zorder=3)
    
    ax_inset.set_extent(ALASKA_EXTENT, crs=ccrs.PlateCarree())
    ax_inset.coastlines(resolution='50m', linewidth=1.8, color='darkgray', zorder=2)
    
    # Alaska label
    ax_inset.text(0.5, 0.99, 'Alaska',
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

def add_light_ocean(ax):
    """Add light ocean color."""
    ax.add_feature(cfeature.OCEAN, facecolor='#1E88E5', alpha=0.25, zorder=0)
    return ax

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def create_panel(ax, data, lons, lats, cmap, norm, panel_label_with_text, 
                 cbar_label, is_delta=False, show_coastal_labels=False):
    """Create a single panel with styling."""
    
    mesh = ax.pcolormesh(lons, lats, data,
                         cmap=cmap, norm=norm,
                         transform=ccrs.PlateCarree(),
                         zorder=3)
    
    ax.set_extent(CONUS_EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='50m', linewidth=1.5, color='#2f4f4f', zorder=2)
    
    # Add coastal region labels ONLY for Panel (a)
    if show_coastal_labels:
        add_coastal_region_labels(ax, CONUS_EXTENT)
    
    add_light_ocean(ax)
    ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', alpha=0.05, zorder=0)
    
    # Add Alaska inset to ALL panels
    add_alaska_inset(ax, data, lons, lats, cmap, norm)
    
    # Panel label with descriptive text (like original)
    ax.text(0.5, 0.98, panel_label_with_text,
           transform=ax.transAxes,
           fontsize=16, fontweight='bold',
           ha='center', va='top',
           bbox=dict(boxstyle='round,pad=0.3',
                   facecolor='white', alpha=0.9,
                   edgecolor='black', linewidth=1.0),
           zorder=10)
    
    # Colorbar
    cbar = plt.colorbar(mesh, ax=ax, orientation='vertical',
                       fraction=0.04, pad=0.03, shrink=0.8)
    cbar.set_label(cbar_label, fontsize=18, fontweight='bold', labelpad=15)
    
    edges = norm.boundaries
    tick_positions = edges[:-1] + np.diff(edges)/2
    
    if is_delta:
        tick_labels = [f"{t:+.2f}" for t in tick_positions]
    else:
        tick_labels = [f"{t:.2f}" for t in tick_positions]
    
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(labelsize=20)
    
    return mesh

# ============================================================================
# MAIN WORKFLOW FUNCTION
# ============================================================================

def create_figure3():
    """Main function to create Figure 3."""
    print("=" * 80)
    print("FIGURE 3: PRE vs POST vs Δ FOR frac_less (WUE$_T$)")
    print("=" * 80)
    print("📊 Loading data...")
    
    # Load PRE and POST data
    pre_data, lats, lons = load_data_and_grid(PRE_FILE, VAR_NAME)
    post_data, _, _ = load_data_and_grid(POST_FILE, VAR_NAME)
    
    # Validate data masks
    both_valid = validate_data_masks(pre_data, post_data)
    
    # Calculate delta
    mask_both = np.isfinite(pre_data) & np.isfinite(post_data)
    delta_data = np.full_like(pre_data, np.nan)
    delta_data[mask_both] = post_data[mask_both] - pre_data[mask_both]
    
    # Calculate bins
    print("\n📈 Calculating shared bins for PRE/POST...")
    shared_edges = calculate_shared_bins(pre_data, post_data)
    print(f"  Shared edges: {[f'{e:.3f}' for e in shared_edges]}")
    
    print("\n📈 Calculating symmetric bins for Δ...")
    delta_edges = calculate_delta_bins(delta_data, mask_both)
    print(f"  Δ edges: {[f'{e:.3f}' for e in delta_edges]}")
    
    # Create colormaps and norms
    main_cmap = ListedColormap(GREEN_TO_GREY, name='GreenToGrey')
    main_norm = BoundaryNorm(shared_edges, main_cmap.N, extend='neither')
    
    delta_cmap = create_delta_colormap()
    delta_norm = BoundaryNorm(delta_edges, delta_cmap.N, extend='both')
    
    # Create figure
    print("\n🎨 Creating figure...")
    fig, axes = plt.subplots(3, 1, figsize=(22, 18),
                            subplot_kw={'projection': ccrs.PlateCarree()},
                            gridspec_kw={'hspace': 0.02})
    
    fig.suptitle('Spatial Change in Declining WUE$_T$ (PRE vs POST Breakpoint)',
                fontsize=22, fontweight='bold', y=0.98)
    
    # Panel (a): PRE - WITH coastal labels, WITH descriptive text
    print("  Creating panel (a): PRE...")
    create_panel(axes[0], pre_data, lons, lats, 
                main_cmap, main_norm,
                "(a) PRE (≤ Breakpoint)", "Frequency of WUE$_T$ Decline", 
                is_delta=False, show_coastal_labels=True)
    
    # Panel (b): POST - NO coastal labels, WITH descriptive text
    print("  Creating panel (b): POST...")
    create_panel(axes[1], post_data, lons, lats,
                main_cmap, main_norm,
                "(b) POST (> Breakpoint)", "Frequency of WUE$_T$ Decline",
                is_delta=False, show_coastal_labels=False)
    
    # Panel (c): Δ - NO coastal labels, WITH descriptive text
    print("  Creating panel (c): Δ...")
    create_panel(axes[2], delta_data, lons, lats,
                delta_cmap, delta_norm,
                "(c) Δ (POST − PRE)", "Change in Frequency of WUE$_T$ Decline",
                is_delta=True, show_coastal_labels=False)
    
    # Adjust layout
    plt.subplots_adjust(left=0.05, right=0.92, top=0.94, bottom=0.03)
    
    # Save figures
    print("\n💾 Saving outputs...")
    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)
    
    fig.savefig(OUTPUT_PNG, dpi=400, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, bbox_inches='tight', facecolor='white')
    
    print(f"✅ PNG saved: {OUTPUT_PNG}")
    print(f"✅ PDF saved: {OUTPUT_PDF}")
    
    # Print statistics
    print("\n📊 COMPREHENSIVE STATISTICS (WUE$_T$):")
    
    pre_valid_data = pre_data[np.isfinite(pre_data)]
    post_valid_data = post_data[np.isfinite(post_data)]
    
    print(f"\n  PRE statistics:")
    print(f"    Mean: {np.mean(pre_valid_data):.3f}")
    print(f"    Std: {np.std(pre_valid_data):.3f}")
    
    print(f"\n  POST statistics:")
    print(f"    Mean: {np.mean(post_valid_data):.3f}")
    print(f"    Std: {np.std(post_valid_data):.3f}")
    
    delta_valid = delta_data[mask_both]
    print(f"\n  Δ STATISTICS (POST − PRE):")
    print(f"    Mean Δ: {np.mean(delta_valid):.3f}")
    print(f"    Median Δ: {np.median(delta_valid):.3f}")
    
    print(f"\n  Δ CHANGE DIRECTION:")
    print(f"    Δ > 0 (increased less-efficient): {np.sum(delta_valid > 0):,} ({100*np.sum(delta_valid > 0)/len(delta_valid):.1f}%)")
    print(f"    Δ < 0 (decreased less-efficient): {np.sum(delta_valid < 0):,} ({100*np.sum(delta_valid < 0)/len(delta_valid):.1f}%)")
    
    plt.show()
    
    print("\n✅ Figure 3 creation complete!")
    print("📍 Coastal labels ONLY in Panel (a)")
    print("📝 Two-row labels for: Southeast Atlantic → 'Southeast\\nAtlantic', Atlantic North → 'Atlantic\\nNorth'")
    print("📍 Labels moved away from coast: Pacific [-126,46], Gulf [-92,28], Southeast Atlantic [-80,33], Atlantic North [-67,43]")
    print("📝 Panel labels with descriptive text restored: (a) PRE (≤ Breakpoint), (b) POST (> Breakpoint), (c) Δ (POST − PRE)")
    print("🎨 Colors: Alaska: #1ABC9C (Teal), Pacific: #E67E22 (Orange), Gulf: #F1C40F (Yellow), Southeast Atlantic: #9B59B6 (Purple), Atlantic North: #3498DB (Blue)")
    
    return fig, axes

# ============================================================================
# EXECUTION
# ============================================================================

if __name__ == "__main__":
    create_figure3()