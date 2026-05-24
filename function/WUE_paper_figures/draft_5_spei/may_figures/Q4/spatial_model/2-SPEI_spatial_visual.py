"""
SPLIT SPEI48 COASTAL DATA & CONVERT TO POLYGON SHAPEFILES
WITH SEPARATE FIGURES FOR EACH REGION
UPDATED FOR 80km BUFFER - OCEAN COAST ONLY (NO GREAT LAKES)
FAST VERSION - USING VECTORIZED OPERATIONS
NOW READS DIRECTLY FROM PREVIOUS STEP'S NETCDF OUTPUTS
KEEPS ALASKA AND CONUS SEPARATE
"""

import xarray as xr
import geopandas as gpd
import numpy as np
import pandas as pd
import glob
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("SPLIT SPEI48 & CONVERT TO POLYGON SHAPEFILES")
print("80km BUFFER - OCEAN COAST ONLY (NO GREAT LAKES)")
print("FAST VERSION - VECTORIZED OPERATIONS")
print("KEEPING ALASKA AND CONUS SEPARATE")
print("="*80)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ====================
# 1. PATHS
# ====================
input_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\spei_subsetting"
output_shape_dir = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\shape"

# Create output directory
output_path = Path(output_shape_dir)
output_path.mkdir(parents=True, exist_ok=True)
print(f"Input directory: {input_dir}")
print(f"Output shapefile directory: {output_path}")

# ====================
# 2. FIND PROCESSED SPEI48 FILES FROM PREVIOUS STEP
# ====================
print("\n" + "="*80)
print("FINDING PROCESSED SPEI48 FILES (80km OCEAN COAST)")
print("="*80)

# Updated pattern to match the output from previous step (FINAL FIX version)
# Your previous code saved files as: SPEI48_genlogistic_global_era5_moda_ref1991to2020_202503_US_OCEAN_COAST_80km.nc
processed_files = glob.glob(os.path.join(input_dir, "*_US_OCEAN_COAST_80km.nc"))
print(f"Found {len(processed_files)} files with pattern '*_US_OCEAN_COAST_80km.nc'")

if not processed_files:
    # Try alternative pattern from earlier version
    processed_files = glob.glob(os.path.join(input_dir, "*_US_OCEAN_COAST_*.nc"))
    print(f"Found {len(processed_files)} files with pattern '*_US_OCEAN_COAST_*.nc'")
    
    if not processed_files:
        # Try any file with ocean coast in name
        processed_files = glob.glob(os.path.join(input_dir, "*OCEAN_COAST*.nc"))
        print(f"Found {len(processed_files)} files with 'OCEAN_COAST' in name")

if not processed_files:
    print("\nERROR: No processed SPEI48 ocean coast files found!")
    print(f"Please check directory: {input_dir}")
    print("\nFiles in directory:")
    all_files = glob.glob(os.path.join(input_dir, "*.nc"))
    for f in all_files[:10]:
        print(f"  {os.path.basename(f)}")
    exit()

# Sort files by date (extract date from filename)
def extract_date_from_filename(filename):
    basename = os.path.basename(filename)
    # Look for YYYYMM pattern
    match = re.search(r'_(\d{4})(\d{2})_', basename)
    if match:
        year = int(match.group(1))
        month = int(match.group(2))
        return year * 100 + month  # Returns YYYYMM as integer for sorting
    return 0

processed_files.sort(key=extract_date_from_filename)

# Show available files
print("\nAvailable processed files:")
for i, filepath in enumerate(processed_files, 1):
    basename = os.path.basename(filepath)
    date_val = extract_date_from_filename(basename)
    if date_val > 0:
        year = date_val // 100
        month = date_val % 100
        print(f"  {i:2d}. {basename} ({year}-{month:02d})")
    else:
        print(f"  {i:2d}. {basename}")

# Select the most recent file for processing (or let user choose)
selected_file = processed_files[-1]
basename = os.path.basename(selected_file)

# Extract date info
date_val = extract_date_from_filename(basename)
if date_val > 0:
    year = date_val // 100
    month = date_val % 100
    date_str = f"{year}-{month:02d}"
else:
    year = "unknown"
    month = "unknown"
    date_str = "unknown"

print(f"\n{'='*80}")
print(f"SELECTED FILE FOR PROCESSING:")
print(f"  File: {basename}")
print(f"  Date: {date_str}")
print(f"{'='*80}")

# ====================
# 3. LOAD THE NETCDF FILE (ALREADY HAS ALASKA AND CONUS SEPARATED?)
# ====================
print("\n" + "="*80)
print("LOADING NETCDF FILE")
print("="*80)

try:
    ds = xr.open_dataset(selected_file)
    print(f"✓ Successfully opened: {basename}")
    
    # Print dataset info
    print(f"\nDataset dimensions: {dict(ds.dims)}")
    print(f"Data variables: {list(ds.data_vars)}")
    print(f"Coordinates: {list(ds.coords)}")
    
    # Find the SPEI variable
    spei_var = None
    for var_name in ds.data_vars:
        if 'SPEI' in var_name.upper() or 'spei' in var_name.lower():
            spei_var = var_name
            break
    
    if not spei_var:
        spei_var = list(ds.data_vars)[0]
        print(f"  Warning: Using variable '{spei_var}' as SPEI")
    else:
        print(f"  Found SPEI variable: {spei_var}")
    
    # Get coordinates
    lat_var = None
    lon_var = None
    
    for var in ['lat', 'latitude', 'Lat', 'Latitude']:
        if var in ds.coords:
            lat_var = var
            break
    
    for var in ['lon', 'longitude', 'Lon', 'Longitude']:
        if var in ds.coords:
            lon_var = var
            break
    
    if not lat_var or not lon_var:
        print(f"ERROR: Could not find lat/lon coordinates")
        print(f"Available coordinates: {list(ds.coords)}")
        exit()
    
    print(f"  Latitude variable: {lat_var}")
    print(f"  Longitude variable: {lon_var}")
    
    # Get data (assuming no time dimension or take first time step)
    if 'time' in ds.dims:
        data = ds[spei_var].isel(time=0).values
        print(f"  Using first time step from {len(ds.time)} available")
    else:
        data = ds[spei_var].values
    
    lats = ds[lat_var].values
    lons = ds[lon_var].values
    
    print(f"\nData shape: {data.shape}")
    print(f"Grid size: {len(lats)} lats × {len(lons)} lons = {len(lats)*len(lons):,} cells")
    print(f"Latitude range: {lats.min():.2f}° to {lats.max():.2f}°")
    print(f"Longitude range: {lons.min():.2f}° to {lons.max():.2f}°")
    
    # Count valid (non-NaN) cells
    valid_cells = np.sum(~np.isnan(data))
    print(f"Valid cells (non-NaN): {valid_cells:,} ({valid_cells/len(lats)/len(lons)*100:.1f}%)")
    
except Exception as e:
    print(f"ERROR loading NetCDF file: {str(e)}")
    import traceback
    traceback.print_exc()
    exit()

# ====================
# 4. SEPARATE ALASKA AND CONUS USING COORDINATES (NOT STATE BOUNDARIES)
# ====================
print("\n" + "="*80)
print("SEPARATING ALASKA AND CONUS USING COORDINATES")
print("="*80)

# Since the NetCDF already has the coastal band, we can separate based on longitude
# Alaska is typically west of -130° longitude
# CONUS is east of -130° longitude (with some overlap in Pacific NW)

# Define separation longitude
separator_lon = -130  # Approximate boundary between Alaska and CONUS

# Create masks based on longitude
lon_grid, lat_grid = np.meshgrid(lons, lats)
alaska_mask = (lon_grid < separator_lon) & ~np.isnan(data)
conus_mask = (lon_grid >= separator_lon) & ~np.isnan(data)

# Also filter by latitude to remove any Arctic ocean points that might be included
alaska_mask = alaska_mask & (lat_grid >= 51)  # Alaska south of ~72°N, north of 51°N
conus_mask = conus_mask & (lat_grid >= 24) & (lat_grid <= 50)  # CONUS lat range

# Apply masks to data
alaska_data = np.where(alaska_mask, data, np.nan)
conus_data = np.where(conus_mask, data, np.nan)

# Count cells
alaska_valid = np.sum(~np.isnan(alaska_data))
conus_valid = np.sum(~np.isnan(conus_data))

print(f"Alaska cells: {alaska_valid:,}")
print(f"CONUS cells: {conus_valid:,}")
print(f"Total: {alaska_valid + conus_valid:,}")

# Check if separation worked
if alaska_valid == 0:
    print("WARNING: No Alaska cells found! Trying different longitude threshold...")
    # Try a different threshold for Alaska
    alaska_mask = (lon_grid < -120) & (lat_grid >= 51) & ~np.isnan(data)
    alaska_data = np.where(alaska_mask, data, np.nan)
    alaska_valid = np.sum(~np.isnan(alaska_data))
    print(f"  Alaska cells with new threshold: {alaska_valid:,}")

if conus_valid == 0:
    print("WARNING: No CONUS cells found!")

# ====================
# 5. CREATE METADATA
# ====================
metadata = {
    'original_file': basename,
    'year': year,
    'month': month,
    'buffer_km': '80',
    'coast_type': 'Ocean coast only (no Great Lakes)',
    'source': 'SPEI48 coastal subsetting workflow'
}

# ====================
# 6. CONVERT TO POLYGON SHAPEFILES
# ====================
print("\n" + "="*80)
print("CONVERTING TO POLYGON SHAPEFILES")
print("="*80)

def array_to_polygon_shapefile(data_array, lats, lons, region_name, metadata, output_path):
    """
    Convert a 2D array with NaN values to a polygon shapefile
    Each non-NaN cell becomes a polygon
    """
    print(f"\nConverting {region_name} to polygon shapefile...")
    
    try:
        # Calculate grid resolution
        if len(lats) > 1:
            lat_res = abs(lats[1] - lats[0])
        else:
            lat_res = 0.25
        
        if len(lons) > 1:
            lon_res = abs(lons[1] - lons[0])
        else:
            lon_res = 0.25
        
        print(f"  Grid resolution: {lat_res:.4f}° lat, {lon_res:.4f}° lon")
        
        # Find valid cells
        valid_indices = np.where(~np.isnan(data_array))
        total_valid = len(valid_indices[0])
        print(f"  Valid cells to convert: {total_valid:,}")
        
        if total_valid == 0:
            print(f"  WARNING: No valid data points found for {region_name}")
            return None, None
        
        from shapely.geometry import Polygon
        
        polygons = []
        spei_values = []
        cell_ids = []
        lat_centers = []
        lon_centers = []
        
        for idx in range(total_valid):
            i = valid_indices[0][idx]  # latitude index
            j = valid_indices[1][idx]  # longitude index
            value = data_array[i, j]
            
            # Calculate polygon bounds
            lat_center = lats[i]
            lon_center = lons[j]
            
            lat_min = lat_center - lat_res/2
            lat_max = lat_center + lat_res/2
            lon_min = lon_center - lon_res/2
            lon_max = lon_center + lon_res/2
            
            # Create polygon
            polygon = Polygon([
                (lon_min, lat_min),
                (lon_min, lat_max),
                (lon_max, lat_max),
                (lon_max, lat_min),
                (lon_min, lat_min)
            ])
            
            polygons.append(polygon)
            spei_values.append(float(value))
            cell_ids.append(f"{region_name[:3]}_{idx:08d}")
            lat_centers.append(lat_center)
            lon_centers.append(lon_center)
            
            # Progress indicator
            if (idx + 1) % 10000 == 0:
                print(f"    Processed {idx+1:,}/{total_valid:,} cells")
        
        print(f"  Created {len(polygons):,} polygons")
        
        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame({
            'cell_id': cell_ids,
            'spei48': spei_values,
            'lat_center': lat_centers,
            'lon_center': lon_centers,
            'region': region_name,
            'year': str(metadata.get('year', 'unknown')),
            'month': str(metadata.get('month', 'unknown')),
            'buffer_km': metadata.get('buffer_km', '80'),
            'coast_type': metadata.get('coast_type', 'ocean_only'),
            'geometry': polygons
        }, crs="EPSG:4326")
        
        # Create output filename
        output_filename = f"SPEI48_{region_name}_{metadata.get('year', 'unknown')}{metadata.get('month', 'unknown'):02d}_80km_ocean_polygons.shp"
        output_file = output_path / output_filename
        
        # Save shapefile
        gdf.to_file(output_file, driver='ESRI Shapefile', encoding='utf-8')
        
        print(f"  ✓ Saved shapefile to: {output_file}")
        
        return gdf, output_file
        
    except Exception as e:
        print(f"  ERROR converting {region_name}: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, None

# Convert Alaska data
alaska_gdf, alaska_shapefile = array_to_polygon_shapefile(
    alaska_data, lats, lons, 'Alaska', metadata, output_path
)

# Convert CONUS data
conus_gdf, conus_shapefile = array_to_polygon_shapefile(
    conus_data, lats, lons, 'Contiguous_US', metadata, output_path
)

if alaska_gdf is None and conus_gdf is None:
    print("ERROR: Failed to create any shapefiles!")
    exit()

# ====================
# 7. CREATE FIGURES FOR EACH REGION
# ====================
print("\n" + "="*80)
print("CREATING FIGURES FOR EACH REGION")
print("="*80)

def create_region_figure(gdf, region_name, metadata, output_path):
    """
    Create a figure for a region's shapefile
    """
    print(f"Creating {region_name} figure...")
    
    if gdf is None or len(gdf) == 0:
        print(f"  No data for {region_name}, skipping figure")
        return None
    
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot 1: Colored by SPEI48
        ax1 = axes[0]
        gdf.plot(column='spei48', ax=ax1, cmap='RdBu', 
                vmin=-3, vmax=3, legend=True,
                legend_kwds={'label': "SPEI48 Value", 'shrink': 0.6},
                edgecolor='black', linewidth=0.1, alpha=0.9)
        ax1.set_title(f"{region_name} - SPEI48 Values\n{len(gdf):,} polygons", 
                     fontsize=12, fontweight='bold')
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Histogram of SPEI48 values
        ax2 = axes[1]
        spei_values = gdf['spei48'].values
        ax2.hist(spei_values, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
        ax2.axvline(x=np.mean(spei_values), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean: {np.mean(spei_values):.3f}')
        ax2.axvline(x=np.median(spei_values), color='green', linestyle='--',
                   linewidth=2, label=f'Median: {np.median(spei_values):.3f}')
        ax2.set_title(f"{region_name} - SPEI48 Distribution\n{metadata.get('year', 'unknown')}-{metadata.get('month', 'unknown'):02d}", 
                     fontsize=12, fontweight='bold')
        ax2.set_xlabel("SPEI48 Value")
        ax2.set_ylabel("Frequency")
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Add statistics text box
        stats_text = (
            f"STATISTICS:\n"
            f"Polygons: {len(gdf):,}\n"
            f"Range: {gdf['spei48'].min():.3f} to {gdf['spei48'].max():.3f}\n"
            f"Mean: {gdf['spei48'].mean():.3f}\n"
            f"Std Dev: {gdf['spei48'].std():.3f}\n"
            f"Buffer: {metadata.get('buffer_km', '80')}km\n"
            f"Coast: Ocean only"
        )
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
                fontsize=8, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.suptitle(f"{region_name.upper()} SPEI48 POLYGON SHAPEFILE - 80km OCEAN COAST", 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        fig_filename = f"{region_name}_80km_ocean_shapefile_{metadata.get('year', 'unknown')}{metadata.get('month', 'unknown'):02d}.png"
        fig_path = output_path / fig_filename
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        print(f"  ✓ Figure saved: {fig_path}")
        
        plt.show()
        return fig_path
        
    except Exception as e:
        print(f"  ERROR creating figure for {region_name}: {str(e)}")
        return None

# Create separate figures
alaska_fig = create_region_figure(alaska_gdf, 'Alaska', metadata, output_path)
conus_fig = create_region_figure(conus_gdf, 'Contiguous_US', metadata, output_path)

# ====================
# 8. CREATE COMBINED VERIFICATION FIGURE
# ====================
print("\n" + "="*80)
print("CREATING COMBINED VERIFICATION FIGURE")
print("="*80)

if alaska_gdf is not None and conus_gdf is not None:
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        # Alaska
        if len(alaska_gdf) > 0:
            alaska_gdf.plot(ax=ax1, column='spei48', cmap='RdBu',
                          vmin=-3, vmax=3, legend=False,
                          edgecolor='black', linewidth=0.1, alpha=0.9)
            ax1.set_title(f"Alaska - {len(alaska_gdf):,} polygons", fontweight='bold')
        else:
            ax1.text(0.5, 0.5, "No Alaska data", ha='center', va='center')
            ax1.set_title("Alaska - No Data", fontweight='bold')
        ax1.set_xlabel("Longitude")
        ax1.set_ylabel("Latitude")
        ax1.set_xlim([-180, -130])
        ax1.set_ylim([51, 72])
        ax1.grid(True, alpha=0.3)
        
        # CONUS
        if len(conus_gdf) > 0:
            conus_gdf.plot(ax=ax2, column='spei48', cmap='RdBu',
                          vmin=-3, vmax=3, legend=True,
                          legend_kwds={'label': "SPEI48 Value", 'shrink': 0.6},
                          edgecolor='black', linewidth=0.1, alpha=0.9)
            ax2.set_title(f"Contiguous US - {len(conus_gdf):,} polygons", fontweight='bold')
        else:
            ax2.text(0.5, 0.5, "No CONUS data", ha='center', va='center')
            ax2.set_title("CONUS - No Data", fontweight='bold')
        ax2.set_xlabel("Longitude")
        ax2.set_ylabel("Latitude")
        ax2.set_xlim([-130, -65])
        ax2.set_ylim([24, 50])
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle(f"SPEI48 SHAPEFILES - 80km OCEAN COAST (KEPT SEPARATE)\n{metadata.get('year', 'unknown')}-{metadata.get('month', 'unknown'):02d}", 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        combined_fig_path = output_path / f"Combined_Verification_{metadata.get('year', 'unknown')}{metadata.get('month', 'unknown'):02d}.png"
        plt.savefig(combined_fig_path, dpi=150, bbox_inches='tight')
        print(f"✓ Combined verification figure saved: {combined_fig_path}")
        
        plt.show()
        
    except Exception as e:
        print(f"ERROR creating combined figure: {str(e)}")

# ====================
# 9. CREATE SUMMARY REPORT
# ====================
print("\n" + "="*80)
print("CREATING SUMMARY REPORT")
print("="*80)

summary_path = output_path / f"Shapefile_Conversion_Summary_{metadata.get('year', 'unknown')}{metadata.get('month', 'unknown'):02d}.txt"

with open(summary_path, 'w', encoding='utf-8') as f:
    f.write("="*70 + "\n")
    f.write("SPEI48 NETCDF TO POLYGON SHAPEFILE CONVERSION\n")
    f.write("80km BUFFER - OCEAN COAST ONLY\n")
    f.write("ALASKA AND CONUS KEPT SEPARATE\n")
    f.write("="*70 + "\n\n")
    
    f.write(f"Processing date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input file: {basename}\n")
    f.write(f"Date: {date_str}\n")
    f.write(f"Buffer: {metadata.get('buffer_km', '80')}km\n")
    f.write(f"Coast type: Ocean only (no Great Lakes)\n\n")
    
    f.write("-"*70 + "\n")
    f.write("CREATED SHAPEFILES:\n")
    f.write("-"*70 + "\n")
    
    if alaska_gdf is not None:
        f.write(f"\n1. ALASKA SHAPEFILE\n")
        f.write(f"   File: {alaska_shapefile.name if alaska_shapefile else 'N/A'}\n")
        f.write(f"   Polygons: {len(alaska_gdf):,}\n")
        if len(alaska_gdf) > 0:
            f.write(f"   SPEI48 range: {alaska_gdf['spei48'].min():.3f} to {alaska_gdf['spei48'].max():.3f}\n")
            f.write(f"   SPEI48 mean: {alaska_gdf['spei48'].mean():.3f}\n")
            f.write(f"   SPEI48 std dev: {alaska_gdf['spei48'].std():.3f}\n")
            f.write(f"   Spatial extent:\n")
            f.write(f"     Longitude: {alaska_gdf['lon_center'].min():.2f}° to {alaska_gdf['lon_center'].max():.2f}°\n")
            f.write(f"     Latitude: {alaska_gdf['lat_center'].min():.2f}° to {alaska_gdf['lat_center'].max():.2f}°\n")
    
    if conus_gdf is not None:
        f.write(f"\n2. CONTIGUOUS US SHAPEFILE\n")
        f.write(f"   File: {conus_shapefile.name if conus_shapefile else 'N/A'}\n")
        f.write(f"   Polygons: {len(conus_gdf):,}\n")
        if len(conus_gdf) > 0:
            f.write(f"   SPEI48 range: {conus_gdf['spei48'].min():.3f} to {conus_gdf['spei48'].max():.3f}\n")
            f.write(f"   SPEI48 mean: {conus_gdf['spei48'].mean():.3f}\n")
            f.write(f"   SPEI48 std dev: {conus_gdf['spei48'].std():.3f}\n")
            f.write(f"   Spatial extent:\n")
            f.write(f"     Longitude: {conus_gdf['lon_center'].min():.2f}° to {conus_gdf['lon_center'].max():.2f}°\n")
            f.write(f"     Latitude: {conus_gdf['lat_center'].min():.2f}° to {conus_gdf['lat_center'].max():.2f}°\n")
    
    total_polygons = (len(alaska_gdf) if alaska_gdf is not None else 0) + \
                     (len(conus_gdf) if conus_gdf is not None else 0)
    f.write(f"\n{'='*70}\n")
    f.write(f"TOTAL POLYGONS CREATED: {total_polygons:,}\n")
    f.write(f"{'='*70}\n")

print(f"✓ Summary saved to: {summary_path}")

# ====================
# 10. FINAL VERIFICATION
# ====================
print("\n" + "="*80)
print("FINAL VERIFICATION")
print("="*80)

print(f"\n✓ SHAPEFILES CREATED SUCCESSFULLY:")
if alaska_gdf is not None:
    print(f"  • Alaska: {len(alaska_gdf):,} polygons (80km ocean coast)")
if conus_gdf is not None:
    print(f"  • Contiguous US: {len(conus_gdf):,} polygons (80km ocean coast, no Great Lakes)")

print(f"\nOutput directory: {output_path}")
print(f"Shapefiles and figures saved to: {output_shape_dir}")

print("\n" + "="*80)
print("WORKFLOW COMPLETE")
print("="*80)
print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")