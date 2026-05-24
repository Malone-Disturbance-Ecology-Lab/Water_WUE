# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 23:05:35 2025

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
CORRECTED WORKFLOW: SPEI Extraction with MODIS Forest/Wetland Mask - FINAL VERSION
File: spei_extraction_final_fixed.py
"""
import xarray as xr
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin, rowcol
from pathlib import Path
import warnings
import time
import os
import re
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directories
SPEI_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\spei_subsetting")
CONUS_MASK_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\processed\conus_masks_500m_FIXED")
ALASKA_MASK_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\processed\alaska_masks_500m_FIXED")

# Output
OUTPUT_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\extraction_final_fixed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Parameters
FW_THRESHOLD = 0.10  # 10% forest/wetland fraction threshold # this threshold will work on other classes too not just forest
MODIS_NODATA = -9999.0

# ============================================================================
# 1. MODIS YEAR RESOLVER
# ============================================================================

def get_available_modis_years():
    """
    Get list of available MODIS years from both CONUS and Alaska directories.
    """
    print("📅 Scanning available MODIS years...")
    
    conus_years = set()
    alaska_years = set()
    
    # Get CONUS years
    for file in CONUS_MASK_DIR.glob("*.tif"):
        match = re.search(r'conus_mask_(\d{4})_500m\.tif', file.name)
        if match:
            conus_years.add(int(match.group(1)))
    
    # Get Alaska years
    for file in ALASKA_MASK_DIR.glob("*.tif"):
        match = re.search(r'alaska_mask_(\d{4})_500m_360\.tif', file.name)
        if match:
            alaska_years.add(int(match.group(1)))
    
    all_years = sorted(conus_years.union(alaska_years))
    
    print(f"  CONUS years: {sorted(conus_years)}")
    print(f"  Alaska years: {sorted(alaska_years)}")
    print(f"  All available MODIS years: {all_years}")
    
    return all_years

def resolve_modis_year(spei_year, available_years):
    """
    Resolve which MODIS year to use for a given SPEI year.
    """
    if spei_year in available_years:
        return spei_year
    
    # Special edge cases
    if spei_year == 2000:
        return 2001  # Explicit rule
    
    # Find nearest year
    closest = min(available_years, key=lambda y: abs(y - spei_year))
    
    # For future years, use most recent
    if spei_year > max(available_years):
        return max(available_years)
    
    return closest

# ============================================================================
# 2. SPEI FILE PARSING
# ============================================================================

def parse_spei_filename(filename):
    """
    Extract year and month from SPEI filename.
    Expected pattern: *_YYYYMM_*.nc
    """
    # Try pattern: _YYYYMM_
    match = re.search(r'_(\d{6})_', filename.name)
    if match:
        yyyymm = match.group(1)
        year = int(yyyymm[:4])
        month = int(yyyymm[4:6])
        return year, month
    
    # Try pattern: YYYYMM
    match = re.search(r'(\d{6})', filename.name)
    if match:
        yyyymm = match.group(1)
        if len(yyyymm) == 6:
            year = int(yyyymm[:4])
            month = int(yyyymm[4:6])
            return year, month
    
    # Try pattern: _YYYY_MM_ or _YYYY-MM_
    match = re.search(r'_(\d{4})[_-](\d{2})_', filename.name)
    if match:
        year = int(match.group(1))
        month = int(match.group(2))
        return year, month
    
    # Last resort: find any 4-digit year
    match = re.search(r'(\d{4})', filename.name)
    if match:
        year = int(match.group(1))
        month = 1  # Default if not found
        return year, month
    
    raise ValueError(f"Cannot parse year from filename: {filename.name}")

# ============================================================================
# 3. AGGREGATION FUNCTIONS
# ============================================================================

def aggregate_conus_to_spei(conus_path, spei_info):
    """
    Aggregate CONUS MODIS mask to SPEI grid using fraction approach.
    """
    print("  Processing CONUS...")
    
    try:
        with rasterio.open(conus_path) as src:
            # Read MODIS data
            modis_data = src.read(1)
            modis_transform = src.transform
            
            # Get bounds
            left, bottom, right, top = src.bounds
            print(f"    CONUS bounds: [{left:.2f}, {bottom:.2f}, {right:.2f}, {top:.2f}]")
            
            # SPEI grid info
            spei_lons = spei_info['lons']
            spei_lats = spei_info['lats']
            lon_res = spei_info['lon_res']
            lat_res = spei_info['lat_res']
            
            # Create target grid transform for SPEI grid
            west = spei_lons[0] - lon_res/2
            east = spei_lons[-1] + lon_res/2
            south = min(spei_lats) - lat_res/2
            north = max(spei_lats) + lat_res/2
            
            width = len(spei_lons)
            height = len(spei_lats)
            
            # Calculate pixel sizes
            x_res = (east - west) / width
            y_res = (south - north) / height
            
            target_transform = from_origin(west, north, x_res, abs(y_res))
            
            # Create masks for aggregation
            forest_mask = (modis_data == 1.0).astype(np.float32)
            valid_mask = (modis_data != MODIS_NODATA).astype(np.float32)
            
            # Initialize output arrays
            forest_sum = np.zeros((height, width), dtype=np.float32)
            valid_sum = np.zeros((height, width), dtype=np.float32)
            
            # Use Proj4 string
            safe_crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
            
            print(f"    Target grid: {height} rows × {width} columns")
            
            # Aggregate
            reproject(
                source=forest_mask,
                destination=forest_sum,
                src_transform=modis_transform,
                src_crs=safe_crs,
                dst_transform=target_transform,
                dst_crs=safe_crs,
                resampling=Resampling.sum
            )
            
            reproject(
                source=valid_mask,
                destination=valid_sum,
                src_transform=modis_transform,
                src_crs=safe_crs,
                dst_transform=target_transform,
                dst_crs=safe_crs,
                resampling=Resampling.sum
            )
            
            # Calculate fraction
            with np.errstate(divide='ignore', invalid='ignore'):
                fraction_grid = np.zeros_like(forest_sum)
                mask = valid_sum > 0
                fraction_grid[mask] = forest_sum[mask] / valid_sum[mask]
                fraction_grid[~mask] = 0.0
            
            # Statistics
            cells_with_data = (valid_sum > 0).sum()
            total_cells = width * height
            print(f"    Cells with MODIS data: {cells_with_data:,}/{total_cells:,}")
            print(f"    Max fraction: {fraction_grid.max():.3f}")
            print(f"    Mean fraction (where >0): {fraction_grid[fraction_grid>0].mean():.3f}")
            
            return {
                'success': True,
                'fraction_grid': fraction_grid,
                'forest_sum': forest_sum,
                'valid_sum': valid_sum
            }
            
    except Exception as e:
        print(f"    ❌ Error: {e}")
        return {'success': False, 'error': str(e)}

def aggregate_alaska_to_spei(alaska_path, spei_info):
    """
    Process Alaska by manual grid intersection.
    """
    print("  Processing Alaska (manual method)...")
    
    try:
        with rasterio.open(alaska_path) as src:
            # Alaska data in 0-360° space
            alaska_data = src.read(1)
            alaska_transform = src.transform
            alaska_height, alaska_width = alaska_data.shape
            
            # Alaska bounds in 0-360°
            left_360, bottom, right_360, top = src.bounds
            print(f"    Alaska bounds (0-360°): [{left_360:.2f}, {bottom:.2f}, {right_360:.2f}, {top:.2f}]")
            
            # SPEI grid info
            spei_lons = spei_info['lons']
            spei_lats = spei_info['lats']
            lon_res = spei_info['lon_res']
            lat_res = spei_info['lat_res']
            
            # Initialize fraction grid
            height = len(spei_lats)
            width = len(spei_lons)
            fraction_grid = np.zeros((height, width), dtype=np.float32)
            
            # Process each SPEI cell that overlaps Alaska
            cells_processed = 0
            alaska_cells_found = 0
            alaska_above_threshold = 0
            
            for lat_idx in range(height):
                lat = spei_lats[lat_idx]
                
                # SPEI cell latitude bounds
                lat_north = lat + lat_res/2
                lat_south = lat - lat_res/2
                
                # Skip if outside Alaska latitude range
                if lat_north < bottom or lat_south > top:
                    continue
                
                for lon_idx in range(width):
                    lon = spei_lons[lon_idx]
                    
                    # Convert SPEI lon to 0-360° for comparison
                    lon_360 = lon if lon >= 0 else lon + 360
                    
                    # SPEI cell longitude bounds in 0-360°
                    lon_west_360 = lon_360 - lon_res/2
                    lon_east_360 = lon_360 + lon_res/2
                    
                    # Handle wrap-around
                    if lon_east_360 > 360:
                        lon_east_360 -= 360
                    
                    # Check if cell overlaps Alaska
                    overlaps = False
                    if (lon_west_360 >= left_360 and lon_west_360 <= right_360) or \
                       (lon_east_360 >= left_360 and lon_east_360 <= right_360):
                        overlaps = True
                    
                    if lon_west_360 > lon_east_360:  # Cell wraps around 360°
                        if (0 >= left_360 and 0 <= right_360) or \
                           (360 >= left_360 and 360 <= right_360):
                            overlaps = True
                    
                    if not overlaps:
                        continue
                    
                    cells_processed += 1
                    
                    try:
                        # Convert geographic bounds to MODIS pixel coordinates
                        row_north, col_west = rowcol(alaska_transform, [lon_west_360], [lat_north])
                        row_south, col_east = rowcol(alaska_transform, [lon_east_360], [lat_south])
                        
                        # Ensure row_north < row_south
                        row_start = min(row_north[0], row_south[0])
                        row_end = max(row_north[0], row_south[0]) + 1
                        col_start = min(col_west[0], col_east[0])
                        col_end = max(col_west[0], col_east[0]) + 1
                        
                        # Clamp to MODIS grid bounds
                        row_start = max(0, row_start)
                        row_end = min(alaska_height, row_end)
                        col_start = max(0, col_start)
                        col_end = min(alaska_width, col_end)
                        
                        if row_end > row_start and col_end > col_start:
                            # Extract MODIS data for this cell
                            cell_data = alaska_data[row_start:row_end, col_start:col_end]
                            
                            # Calculate statistics
                            valid_pixels = cell_data[cell_data != MODIS_NODATA]
                            
                            if len(valid_pixels) > 0:
                                forest_pixels = (valid_pixels == 1.0).sum()
                                fraction = forest_pixels / len(valid_pixels)
                                fraction_grid[lat_idx, lon_idx] = fraction
                                alaska_cells_found += 1
                                
                                # Count cells above threshold
                                if fraction >= FW_THRESHOLD:
                                    alaska_above_threshold += 1
                                
                    except Exception as cell_error:
                        continue
            
            print(f"    SPEI cells checked: {cells_processed:,}")
            print(f"    Alaska cells found (fraction > 0): {alaska_cells_found:,}")
            print(f"    Alaska cells above threshold: {alaska_above_threshold:,}")
            
            if alaska_cells_found > 0:
                print(f"    Max fraction: {fraction_grid.max():.3f}")
                print(f"    Mean fraction: {fraction_grid[fraction_grid>0].mean():.3f}")
            else:
                print(f"    ⚠️ No Alaska cells found in SPEI grid")
            
            return {
                'success': True,
                'fraction_grid': fraction_grid,
                'alaska_cells_found': alaska_cells_found,
                'alaska_above_threshold': alaska_above_threshold
            }
            
    except Exception as e:
        print(f"    ❌ Error: {e}")
        return {'success': False, 'error': str(e)}

# ============================================================================
# 4. COMBINE REGIONS AND CREATE MASK
# ============================================================================

def combine_regions_and_create_mask(conus_result, alaska_result, spei_info):
    """
    Combine CONUS and Alaska fraction grids, create binary mask at threshold.
    """
    print("  Combining regions...")
    
    # Start with CONUS
    combined_fraction = conus_result['fraction_grid'].copy()
    
    # Add Alaska if available and successful
    alaska_added = 0
    alaska_above_threshold = 0
    if alaska_result and alaska_result['success']:
        alaska_fraction = alaska_result['fraction_grid']
        
        # Find where Alaska has data
        alaska_mask = alaska_fraction > 0
        
        # Add Alaska data (shouldn't overlap with CONUS)
        if alaska_mask.any():
            combined_fraction[alaska_mask] = alaska_fraction[alaska_mask]
            alaska_added = alaska_mask.sum()
            alaska_above_threshold = alaska_result.get('alaska_above_threshold', 0)
            print(f"    Alaska cells added: {alaska_added:,}")
            print(f"    Alaska cells above threshold: {alaska_above_threshold:,}")
        else:
            print(f"    ⚠️ Alaska processed but no cells found in overlap")
    else:
        print(f"    No Alaska data to combine")
    
    # Create binary mask at threshold
    binary_mask = (combined_fraction >= FW_THRESHOLD).astype(np.int8)
    
    # Statistics
    total_cells = combined_fraction.size
    cells_with_forest = (combined_fraction > 0).sum()
    cells_above_threshold = binary_mask.sum()
    
    print(f"    Total cells with forest/wetland: {cells_with_forest:,}")
    print(f"    Cells above {FW_THRESHOLD*100:.0f}% threshold: {cells_above_threshold:,}")
    print(f"    Percentage of SPEI grid: {cells_above_threshold/total_cells*100:.2f}%")
    
    # Calculate CONUS contribution
    conus_above_threshold = cells_above_threshold - alaska_above_threshold
    print(f"    CONUS cells above threshold: {conus_above_threshold:,}")
    
    return {
        'combined_fraction': combined_fraction,
        'binary_mask': binary_mask,
        'stats': {
            'total_cells': total_cells,
            'cells_with_forest': cells_with_forest,
            'cells_above_threshold': cells_above_threshold,
            'conus_above_threshold': conus_above_threshold,
            'alaska_above_threshold': alaska_above_threshold,
            'alaska_added': alaska_added
        }
    }

# ============================================================================
# 5. APPLY MASK TO SPEI DATA WITH SANITY CHECKS
# ============================================================================

def apply_mask_to_spei(spei_data, binary_mask, fill_value):
    """
    Apply binary mask to SPEI data with proper NaN handling and sanity checks.
    """
    print("  Applying mask to SPEI...")
    
    # Convert to float for NaN handling
    data_float = spei_data.astype(np.float32)
    
    # Identify valid SPEI pixels (finite numbers AND not fill value)
    valid_mask = np.isfinite(data_float)  # This excludes NaN and inf
    if fill_value is not None:
        valid_mask = valid_mask & (data_float != fill_value)
    
    # Count original valid pixels
    original_valid = valid_mask.sum()
    
    # Apply forest/wetland mask
    masked_data = np.full_like(data_float, fill_value, dtype=np.float32)
    extraction_mask = (binary_mask == 1) & valid_mask
    
    # Copy valid SPEI values where mask is 1
    masked_data[extraction_mask] = data_float[extraction_mask]
    
    # Count masked valid pixels
    masked_valid = extraction_mask.sum()
    
    # SANITY CHECK 1: Candidate pixels (fw_mask=1 AND SPEI valid)
    candidate_mask = (binary_mask == 1) & valid_mask
    candidate_count = candidate_mask.sum()
    
    # SANITY CHECK 2: Masked_valid should equal candidate_count
    sanity_check1 = masked_valid == candidate_count
    sanity_check2 = masked_valid <= original_valid
    
    if not sanity_check1:
        print(f"    ⚠️  SANITY CHECK FAILED: masked_valid ({masked_valid}) != candidate_count ({candidate_count})")
    if not sanity_check2:
        print(f"    ⚠️  SANITY CHECK FAILED: masked_valid ({masked_valid}) > original_valid ({original_valid})")
    
    # Calculate range of extracted values
    if masked_valid > 0:
        extracted_values = data_float[extraction_mask]
        masked_range = [float(np.min(extracted_values)), float(np.max(extracted_values))]
        masked_mean = float(np.mean(extracted_values))
        masked_std = float(np.std(extracted_values))
        
        # Also calculate range of original valid values for comparison
        original_values = data_float[valid_mask]
        original_range = [float(np.min(original_values)), float(np.max(original_values))]
        original_mean = float(np.mean(original_values))
        original_std = float(np.std(original_values))
    else:
        masked_range = [np.nan, np.nan]
        masked_mean = np.nan
        masked_std = np.nan
        original_range = [np.nan, np.nan]
        original_mean = np.nan
        original_std = np.nan
    
    # Calculate percentages
    if original_valid > 0:
        retention = masked_valid / original_valid * 100
        reduction = (original_valid - masked_valid) / original_valid * 100
    else:
        retention = 0.0
        reduction = 100.0
    
    print(f"    Original valid SPEI pixels: {original_valid:,}")
    print(f"    Masked valid SPEI pixels: {masked_valid:,}")
    print(f"    Candidate pixels (fw_mask=1 & SPEI valid): {candidate_count:,}")
    print(f"    Retention: {retention:.1f}% (Extracted/Original)")
    print(f"    Reduction: {reduction:.1f}% (Removed/Original)")
    print(f"    Sanity checks: {'PASS' if sanity_check1 and sanity_check2 else 'FAIL'}")
    
    if masked_valid > 0:
        print(f"    Original SPEI range: [{original_range[0]:.2f}, {original_range[1]:.2f}]")
        print(f"    Original SPEI mean: {original_mean:.2f}")
        print(f"    Masked SPEI range: [{masked_range[0]:.2f}, {masked_range[1]:.2f}]")
        print(f"    Masked SPEI mean: {masked_mean:.2f}")
    
    return {
        'masked_data': masked_data,
        'original_valid': int(original_valid),
        'masked_valid': int(masked_valid),
        'candidate_count': int(candidate_count),
        'sanity_check1': bool(sanity_check1),
        'sanity_check2': bool(sanity_check2),
        'retention_percent': float(retention),
        'reduction_percent': float(reduction),
        'original_range': original_range,
        'original_mean': float(original_mean),
        'original_std': float(original_std),
        'masked_range': masked_range,
        'masked_mean': float(masked_mean),
        'masked_std': float(masked_std)
    }

# ============================================================================
# 6. SAVE RESULTS WITH PROPER FILENAME (FIXED)
# ============================================================================

def save_results(spei_data, masked_data, combined_results, extraction_stats, 
                 spei_info, spei_filename, spei_year, spei_month, modis_year, timestamp):
    """
    Save all results to NetCDF file with proper naming.
    """
    print("  Saving results to NetCDF...")
    
    # Create dataset
    ds_out = xr.Dataset(
        {
            'SPEI48_original': (['lat', 'lon'], spei_data),
            'SPEI48_masked': (['lat', 'lon'], masked_data),
            'fw_fraction': (['lat', 'lon'], combined_results['combined_fraction']),
            'fw_mask': (['lat', 'lon'], combined_results['binary_mask']),
        },
        coords={
            'lon': spei_info['lons'],
            'lat': spei_info['lats'],
        }
    )
    
    # Add comprehensive attributes
    # Convert boolean to int for NetCDF compatibility
    sanity_passed = extraction_stats['sanity_check1'] and extraction_stats['sanity_check2']
    
    ds_out.attrs.update({
        'title': 'SPEI Extracted for Forest/Wetland Areas',
        'description': 'SPEI48 values extracted only for pixels with ≥10% forest/wetland cover',
        'SPEI_file': spei_filename,
        'SPEI_year': str(spei_year),
        'SPEI_month': str(spei_month),
        'MODIS_year': str(modis_year),
        'threshold_used': float(FW_THRESHOLD),
        'threshold_description': f'Cells with ≥{FW_THRESHOLD*100:.0f}% forest/wetland cover',
        'processing_date': timestamp,
        'grid_resolution': f"{spei_info['lon_res']:.3f}° × {spei_info['lat_res']:.3f}°",
        'extraction_method': 'Fraction-based aggregation from 500m MODIS to 0.25° SPEI grid',
        'valid_spei_pixels': int(extraction_stats['original_valid']),
        'extracted_spei_pixels': int(extraction_stats['masked_valid']),
        'retention_percentage': float(extraction_stats['retention_percent']),
        'original_spei_mean': float(extraction_stats['original_mean']),
        'extracted_spei_mean': float(extraction_stats['masked_mean']),
        'sanity_check_passed': int(sanity_passed)  # Convert bool to int
    })
    
    # Variable-specific attributes
    ds_out['SPEI48_original'].attrs.update({
        'long_name': 'Standardized Precipitation Evapotranspiration Index (48-month)',
        'units': 'standardized units',
        '_FillValue': float(spei_info['fill_value']),
        'description': 'Original SPEI48 data before masking',
        'valid_pixel_count': int(extraction_stats['original_valid'])
    })
    
    ds_out['SPEI48_masked'].attrs.update({
        'long_name': 'SPEI48 for forest/wetland areas only',
        'units': 'standardized units',
        '_FillValue': float(spei_info['fill_value']),
        'description': 'SPEI48 values extracted only where forest/wetland fraction ≥ 10%',
        'extracted_pixel_count': int(extraction_stats['masked_valid'])
    })
    
    ds_out['fw_fraction'].attrs.update({
        'long_name': 'Forest/Wetland fraction per SPEI cell',
        'units': 'fraction (0-1)',
        '_FillValue': 0.0,
        'description': 'Fraction of MODIS 500m pixels within each 0.25° SPEI cell that are forest/wetland'
    })
    
    ds_out['fw_mask'].attrs.update({
        'long_name': 'Forest/Wetland binary mask',
        'units': 'binary (0/1)',
        '_FillValue': 0,
        'description': f'1 = cell has ≥{FW_THRESHOLD*100:.0f}% forest/wetland cover, 0 = otherwise'
    })
    
    # Save to file with proper naming
    base_name = Path(spei_filename).stem
    nc_filename = f"{base_name}_fw_processed.nc"
    nc_path = OUTPUT_DIR / nc_filename
    
    ds_out.to_netcdf(nc_path)
    
    file_size_mb = nc_path.stat().st_size / (1024*1024)
    print(f"    ✅ NetCDF saved: {nc_filename} ({file_size_mb:.1f} MB)")
    
    return nc_path

# ============================================================================
# 7. MAIN PROCESSING PIPELINE
# ============================================================================

def process_single_spei_file(spei_path, modis_year):
    """
    Process a single SPEI file with all fixes and sanity checks.
    """
    print("="*70)
    print(f"PROCESSING: {spei_path.name}")
    print(f"Using MODIS year: {modis_year}")
    print("="*70)
    
    start_time = time.time()
    
    try:
        # Parse SPEI year and month from filename
        spei_year, spei_month = parse_spei_filename(spei_path)
        
        # Step 1: Load SPEI data
        print("\n1. Loading SPEI data...")
        with xr.open_dataset(spei_path) as ds:
            spei_info = {
                'lons': ds.lon.values,
                'lats': ds.lat.values,
                'lon_res': abs(ds.lon.values[1] - ds.lon.values[0]),
                'lat_res': abs(ds.lat.values[1] - ds.lat.values[0]),
                'fill_value': float(ds['SPEI48'].attrs.get('_FillValue', -9999)),
                'spei_data': ds['SPEI48'].isel(time=0).values
            }
        
        print(f"   Grid: {len(spei_info['lats'])}×{len(spei_info['lons'])}")
        print(f"   Resolution: {spei_info['lon_res']:.3f}° × {spei_info['lat_res']:.3f}°")
        print(f"   Fill value: {spei_info['fill_value']}")
        
        # Step 2: Load MODIS masks
        print("\n2. Loading MODIS masks...")
        
        # CONUS mask
        conus_file = CONUS_MASK_DIR / f"conus_mask_{modis_year}_500m.tif"
        if not conus_file.exists():
            # Try alternative naming
            alt_files = list(CONUS_MASK_DIR.glob(f"*{modis_year}*.tif"))
            if alt_files:
                conus_file = alt_files[0]
            else:
                raise FileNotFoundError(f"CONUS mask not found for year {modis_year}")
        
        # Alaska mask
        alaska_file = ALASKA_MASK_DIR / f"alaska_mask_{modis_year}_500m_360.tif"
        if not alaska_file.exists():
            # Try alternative naming
            alt_files = list(ALASKA_MASK_DIR.glob(f"*{modis_year}*.tif"))
            if alt_files:
                alaska_file = alt_files[0]
            else:
                print(f"   ⚠️  Alaska mask not found for year {modis_year}")
                alaska_file = None
        
        # Step 3: Aggregate MODIS to SPEI grid
        print("\n3. Aggregating MODIS to SPEI grid...")
        conus_result = aggregate_conus_to_spei(conus_file, spei_info)
        
        if not conus_result['success']:
            raise RuntimeError(f"CONUS aggregation failed: {conus_result.get('error', 'Unknown error')}")
        
        alaska_result = None
        if alaska_file:
            alaska_result = aggregate_alaska_to_spei(alaska_file, spei_info)
            if not alaska_result['success']:
                print(f"   ⚠️  Alaska aggregation failed: {alaska_result.get('error', 'Unknown error')}")
                alaska_result = None
        
        # Step 4: Combine regions and create mask
        print("\n4. Creating forest/wetland mask...")
        combined_results = combine_regions_and_create_mask(
            conus_result, alaska_result, spei_info
        )
        
        # Step 5: Apply mask to SPEI with sanity checks
        print("\n5. Extracting SPEI data...")
        extraction_stats = apply_mask_to_spei(
            spei_info['spei_data'],
            combined_results['binary_mask'],
            spei_info['fill_value']
        )
        
        # Check sanity checks
        if not extraction_stats['sanity_check1'] or not extraction_stats['sanity_check2']:
            print(f"   ⚠️  WARNING: Sanity checks failed! Investigate before proceeding.")
        
        # Step 6: Save results
        print("\n6. Saving results...")
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        
        # Save NetCDF with proper naming
        nc_path = save_results(
            spei_info['spei_data'],
            extraction_stats['masked_data'],
            combined_results,
            extraction_stats,
            spei_info,
            spei_path.name,
            spei_year,
            spei_month,
            modis_year,
            timestamp
        )
        
        # Step 7: Final summary
        elapsed_time = time.time() - start_time
        
        print("\n" + "="*70)
        print("✅ PROCESSING COMPLETE")
        print("="*70)
        print(f"Processing time: {elapsed_time:.1f} seconds")
        
        # Calculate percentages correctly
        total_grid_cells = spei_info['spei_data'].size
        
        print(f"\n📊 COMPREHENSIVE STATISTICS:")
        print(f"  Grid size: {total_grid_cells:,} cells")
        print(f"  Forest/Wetland cells (≥{FW_THRESHOLD*100:.0f}%): {combined_results['stats']['cells_above_threshold']:,}")
        print(f"    - CONUS: {combined_results['stats']['conus_above_threshold']:,}")
        print(f"    - Alaska: {combined_results['stats']['alaska_above_threshold']:,}")
        print(f"  Forest % of grid: {combined_results['stats']['cells_above_threshold']/total_grid_cells*100:.2f}%")
        
        print(f"\n🌡️ SPEI STATISTICS:")
        print(f"  Valid SPEI pixels: {extraction_stats['original_valid']:,}")
        print(f"    - % of grid: {extraction_stats['original_valid']/total_grid_cells*100:.2f}%")
        print(f"  Extracted SPEI pixels: {extraction_stats['masked_valid']:,}")
        print(f"    - % of valid SPEI: {extraction_stats['retention_percent']:.1f}%")
        print(f"    - % of grid: {extraction_stats['masked_valid']/total_grid_cells*100:.2f}%")
        
        print(f"\n📈 VALUE RANGES:")
        print(f"  Original SPEI: [{extraction_stats['original_range'][0]:.2f}, {extraction_stats['original_range'][1]:.2f}]")
        print(f"  Extracted SPEI: [{extraction_stats['masked_range'][0]:.2f}, {extraction_stats['masked_range'][1]:.2f}]")
        print(f"  Original mean: {extraction_stats['original_mean']:.2f}")
        print(f"  Extracted mean: {extraction_stats['masked_mean']:.2f}")
        
        print(f"\n🔍 SANITY CHECKS:")
        print(f"  Check 1 (masked_valid == candidate): {'PASS' if extraction_stats['sanity_check1'] else 'FAIL'}")
        print(f"  Check 2 (masked_valid ≤ original_valid): {'PASS' if extraction_stats['sanity_check2'] else 'FAIL'}")
        
        # Interpretation
        retention = extraction_stats['retention_percent']
        if retention < 20:
            interpretation = "Very selective extraction"
        elif retention < 40:
            interpretation = "Moderate extraction"
        elif retention < 60:
            interpretation = "Substantial extraction"
        else:
            interpretation = "High extraction"
        
        print(f"\n🔍 INTERPRETATION: {interpretation} ({retention:.1f}% of valid SPEI retained)")
        
        print(f"\n💾 OUTPUT FILE:")
        print(f"  {nc_path.name}")
        
        return {
            'success': True,
            'nc_file': nc_path,
            'spei_year': spei_year,
            'spei_month': spei_month,
            'modis_year': modis_year,
            'stats': {
                'forest_cells': combined_results['stats']['cells_above_threshold'],
                'conus_forest': combined_results['stats']['conus_above_threshold'],
                'alaska_forest': combined_results['stats']['alaska_above_threshold'],
                'original_valid': extraction_stats['original_valid'],
                'masked_valid': extraction_stats['masked_valid'],
                'retention_percent': extraction_stats['retention_percent'],
                'original_mean': extraction_stats['original_mean'],
                'masked_mean': extraction_stats['masked_mean'],
                'sanity_check1': extraction_stats['sanity_check1'],
                'sanity_check2': extraction_stats['sanity_check2'],
                'processing_time': elapsed_time
            }
        }
        
    except Exception as e:
        print(f"\n❌ PROCESSING FAILED: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}

# ============================================================================
# 8. BATCH PROCESSING (ALL FILES)
# ============================================================================

# ============================================================================
# 8. BATCH PROCESSING (ALL FILES) - COMPLETED
# ============================================================================

def process_all_files():
    """
    Process ALL SPEI files with proper MODIS year resolution.
    """
    print("="*70)
    print("BATCH PROCESSING: ALL SPEI FILES")
    print("="*70)
    
    # Get available MODIS years
    available_years = get_available_modis_years()
    if not available_years:
        print("❌ No MODIS years found!")
        return []
    
    # Get all SPEI files
    spei_files = sorted(SPEI_DIR.glob("*.nc"))
    if not spei_files:
        print("❌ No SPEI files found!")
        return []
    
    print(f"\nFound {len(spei_files)} SPEI files")
    print(f"Available MODIS years: {available_years}")
    
    all_results = []
    failed_files = []
    
    # Process each SPEI file
    for i, spei_file in enumerate(spei_files, 1):
        print(f"\n{'='*60}")
        print(f"FILE {i}/{len(spei_files)}: {spei_file.name}")
        print(f"{'='*60}")
        
        try:
            # Parse year from filename
            spei_year, spei_month = parse_spei_filename(spei_file)
            
            # Resolve MODIS year
            modis_year = resolve_modis_year(spei_year, available_years)
            
            print(f"  SPEI date: {spei_year}-{spei_month:02d}")
            print(f"  Using MODIS year: {modis_year}")
            
            # Process the file
            result = process_single_spei_file(spei_file, modis_year)
            
            if result['success']:
                all_results.append(result)
                print(f"  ✅ Success")
            else:
                failed_files.append((spei_file.name, result.get('error', 'Unknown error')))
                print(f"  ❌ Failed")
                
        except Exception as e:
            failed_files.append((spei_file.name, str(e)))
            print(f"  ❌ Failed with error: {e}")
    
    # Save batch summary
    save_batch_summary(all_results, failed_files)
    
    # Final summary
    print("\n" + "="*70)
    print("BATCH PROCESSING COMPLETE")
    print("="*70)
    
    successful = len(all_results)
    total = len(spei_files)
    
    print(f"Successfully processed: {successful}/{total} files ({successful/total*100:.1f}%)")
    
    if successful > 0:
        print(f"\n📊 OVERALL STATISTICS:")
        
        # Calculate average statistics
        retention_percents = [r['stats']['retention_percent'] for r in all_results]
        forest_counts = [r['stats']['forest_cells'] for r in all_results]
        
        print(f"  Average retention: {np.mean(retention_percents):.1f}%")
        print(f"  Min retention: {np.min(retention_percents):.1f}%")
        print(f"  Max retention: {np.max(retention_percents):.1f}%")
        print(f"  Average forest cells: {np.mean(forest_counts):.0f}")
        
        # Check sanity checks
        sanity_check1_passed = sum(1 for r in all_results if r['stats']['sanity_check1'])
        sanity_check2_passed = sum(1 for r in all_results if r['stats']['sanity_check2'])
        
        print(f"\n🔍 SANITY CHECK RESULTS:")
        print(f"  Check 1 passed: {sanity_check1_passed}/{successful}")
        print(f"  Check 2 passed: {sanity_check2_passed}/{successful}")
        
        # Check for red flags
        low_retention = sum(1 for p in retention_percents if p < 10)
        high_retention = sum(1 for p in retention_percents if p > 90)
        
        if low_retention > 0:
            print(f"  ⚠️  {low_retention} files have very low retention (<10%)")
        if high_retention > 0:
            print(f"  ⚠️  {high_retention} files have very high retention (>90%)")
    
    if failed_files:
        print(f"\n⚠️ FAILED FILES ({len(failed_files)}):")
        for filename, error in failed_files[:10]:  # Show first 10
            print(f"  {filename}: {error[:80]}...")
        if len(failed_files) > 10:
            print(f"  ... and {len(failed_files)-10} more")
    
    print(f"\n💾 OUTPUT DIRECTORY: {OUTPUT_DIR}")
    print(f"   Contains {successful} NetCDF files")
    
    return all_results

# ============================================================================
# 9. TEST MULTIPLE FILES FUNCTION
# ============================================================================

def test_multiple_files(num_files=3):
    """
    Test the pipeline on multiple files.
    """
    print("="*70)
    print(f"TESTING {num_files} FILES")
    print("="*70)
    
    # Get available MODIS years
    available_years = get_available_modis_years()
    
    # Get SPEI files
    spei_files = sorted(SPEI_DIR.glob("*.nc"))
    if not spei_files:
        print("❌ No SPEI files found!")
        return []
    
    # Limit to requested number
    test_files = spei_files[:min(num_files, len(spei_files))]
    
    print(f"Testing {len(test_files)} files:")
    for f in test_files:
        print(f"  - {f.name}")
    
    all_results = []
    failed_files = []
    
    # Process each test file
    for i, spei_file in enumerate(test_files, 1):
        print(f"\n{'='*60}")
        print(f"TEST FILE {i}/{len(test_files)}: {spei_file.name}")
        print(f"{'='*60}")
        
        try:
            # Parse year from filename
            spei_year, spei_month = parse_spei_filename(spei_file)
            
            # Resolve MODIS year
            modis_year = resolve_modis_year(spei_year, available_years)
            
            print(f"  SPEI date: {spei_year}-{spei_month:02d}")
            print(f"  Using MODIS year: {modis_year}")
            
            # Process the file
            result = process_single_spei_file(spei_file, modis_year)
            
            if result['success']:
                all_results.append(result)
                print(f"  ✅ Success")
            else:
                failed_files.append((spei_file.name, result.get('error', 'Unknown error')))
                print(f"  ❌ Failed")
                
        except Exception as e:
            failed_files.append((spei_file.name, str(e)))
            print(f"  ❌ Failed with error: {e}")
    
    # Test summary
    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70)
    
    successful = len(all_results)
    total = len(test_files)
    
    print(f"Successfully processed: {successful}/{total} test files ({successful/total*100:.1f}%)")
    
    if successful > 0:
        print(f"\n📊 TEST STATISTICS:")
        
        # Show results for each file
        for result in all_results:
            stats = result['stats']
            print(f"\n  {result['nc_file'].name}:")
            print(f"    Forest cells: {stats['forest_cells']:,} (CONUS: {stats['conus_forest']:,}, Alaska: {stats['alaska_forest']:,})")
            print(f"    Retention: {stats['retention_percent']:.1f}%")
            print(f"    Sanity checks: {'PASS' if stats['sanity_check1'] and stats['sanity_check2'] else 'FAIL'}")
    
    if failed_files:
        print(f"\n⚠️ FAILED TEST FILES:")
        for filename, error in failed_files:
            print(f"  {filename}: {error}")
    
    print(f"\n💾 Output saved to: {OUTPUT_DIR}")
    
    return all_results

# ============================================================================
# 10. BATCH SUMMARY
# ============================================================================

def save_batch_summary(all_results, failed_files):
    """
    Save comprehensive batch processing summary.
    """
    summary_path = OUTPUT_DIR / "batch_processing_summary.txt"
    
    with open(summary_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SPEI EXTRACTION BATCH PROCESSING SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Total files processed: {len(all_results)}\n")
        f.write(f"Failed files: {len(failed_files)}\n\n")
        
        # Summary statistics
        if all_results:
            f.write("OVERALL STATISTICS:\n")
            f.write("-" * 80 + "\n")
            
            retention_percents = [r['stats']['retention_percent'] for r in all_results]
            forest_counts = [r['stats']['forest_cells'] for r in all_results]
            conus_counts = [r['stats']['conus_forest'] for r in all_results]
            alaska_counts = [r['stats']['alaska_forest'] for r in all_results]
            sanity_check1_passed = sum(1 for r in all_results if r['stats']['sanity_check1'])
            sanity_check2_passed = sum(1 for r in all_results if r['stats']['sanity_check2'])
            
            f.write(f"Retention (% of valid SPEI extracted):\n")
            f.write(f"  Mean: {np.mean(retention_percents):.1f}%\n")
            f.write(f"  Min: {np.min(retention_percents):.1f}%\n")
            f.write(f"  Max: {np.max(retention_percents):.1f}%\n")
            f.write(f"  Std: {np.std(retention_percents):.1f}%\n\n")
            
            f.write(f"Forest/Wetland cells (≥{FW_THRESHOLD*100:.0f}%):\n")
            f.write(f"  Mean total: {np.mean(forest_counts):.0f}\n")
            f.write(f"  Mean CONUS: {np.mean(conus_counts):.0f}\n")
            f.write(f"  Mean Alaska: {np.mean(alaska_counts):.0f}\n\n")
            
            f.write(f"Sanity Checks:\n")
            f.write(f"  Check 1 (masked_valid == candidate): {sanity_check1_passed}/{len(all_results)} passed\n")
            f.write(f"  Check 2 (masked_valid ≤ original_valid): {sanity_check2_passed}/{len(all_results)} passed\n\n")
        
        # File-by-file details
        f.write("FILE-BY-FILE DETAILS:\n")
        f.write("="*80 + "\n")
        
        for result in all_results:
            stats = result['stats']
            f.write(f"\n{result['nc_file'].name}:\n")
            f.write(f"  SPEI date: {result['spei_year']}-{result['spei_month']:02d}\n")
            f.write(f"  MODIS year: {result['modis_year']}\n")
            f.write(f"  Forest cells: {stats['forest_cells']:,} (CONUS: {stats['conus_forest']:,}, Alaska: {stats['alaska_forest']:,})\n")
            f.write(f"  Valid SPEI: {stats['original_valid']:,}\n")
            f.write(f"  Extracted SPEI: {stats['masked_valid']:,}\n")
            f.write(f"  Retention: {stats['retention_percent']:.1f}%\n")
            f.write(f"  Sanity checks: {'PASS' if stats['sanity_check1'] and stats['sanity_check2'] else 'FAIL'}\n")
        
        # Failed files
        if failed_files:
            f.write("\n\nFAILED FILES:\n")
            f.write("="*80 + "\n")
            for filename, error in failed_files:
                f.write(f"\n{filename}:\n  {error}\n")
    
    print(f"  ✅ Batch summary saved: {summary_path.name}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """
    Main execution function.
    """
    print("SPEI EXTRACTION WORKFLOW - FINAL FIXED VERSION")
    print("="*70)
    print("Extract SPEI drought data for forest/wetland areas only")
    print(f"Threshold: ≥{FW_THRESHOLD*100:.0f}% forest/wetland cover")
    print(f"Output directory: {OUTPUT_DIR}")
    print("-" * 70)
    
    print("\nOptions:")
    print("1. Process all SPEI files (full batch)")
    print("2. Test single file")
    print("3. Test 3 files")
    print("4. List available MODIS years")
    print("5. Exit")
    
    choice = input("\nEnter choice (1-5): ").strip()
    
    if choice == "1":
        print("\nStarting full batch processing...")
        results = process_all_files()
        
    elif choice == "2":
        print("\nTesting single file...")
        # Get first file
        spei_files = sorted(SPEI_DIR.glob("*.nc"))
        if spei_files:
            available_years = get_available_modis_years()
            test_file = spei_files[0]
            spei_year, spei_month = parse_spei_filename(test_file)
            modis_year = resolve_modis_year(spei_year, available_years)
            result = process_single_spei_file(test_file, modis_year)
        else:
            print("❌ No SPEI files found!")
            
    elif choice == "3":
        print("\nTesting 3 files...")
        results = test_multiple_files(num_files=3)
        
    elif choice == "4":
        available_years = get_available_modis_years()
        print(f"\nAvailable MODIS years: {available_years}")
        
    else:
        print("\nExiting.")

# ============================================================================

if __name__ == "__main__":
    main()