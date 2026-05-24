# -*- coding: utf-8 -*-
"""
Created on Wed May  6 21:43:18 2026

@author: ammar
"""

"""
STEP 1B: Process Alaska MODIS CSV to 500m Binary Masks - FIXED
File: step2_alaska_FIXED.py
UPDATED: TARGET_CLASSES based on 57 sites in final figure (strict triple intersection)
"""
import pandas as pd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from pathlib import Path
import os
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION - ALASKA
# ============================================================================

# MODIS Land Cover Classes for FINAL FIGURE SITES (57 sites in strict triple intersection)
# Based on actual IGBP biomes from the 57 sites:
# ENF, DBF, MF, CSH, OSH, GRA, WET, CRO, BSV
TARGET_CLASSES = [1, 4, 5, 6, 7, 10, 11, 12, 16]
# Class details for reference:
#   1  = Evergreen Needleleaf Forests (ENF)
#   4  = Deciduous Broadleaf Forests (DBF)
#   5  = Mixed Forests (MF)
#   6  = Closed Shrublands (CSH)
#   7  = Open Shrublands (OSH)
#   10 = Grasslands (GRA)
#   11 = Permanent Wetlands (WET)
#   12 = Croplands (CRO)
#   16 = Barren (BSV)

# Alaska paths (SAME - DO NOT CHANGE)
ALASKA_INPUT_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\MODIS_landcover\alaska"
ALASKA_OUTPUT_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\processed\alaska_masks_500m_FIXED")

MODIS_RES = 0.0045

# CRS for Alaska in 0-360° space
CRS_360_PROJ4 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=360 +no_defs"

ALASKA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# FUNCTION: Create Alaska Mask - FIXED
# ============================================================================

def create_alaska_mask(csv_path, output_tif_path):
    """
    Create Alaska binary mask in 0-360° coordinate space.
    """
    filename = Path(csv_path).name
    print(f"\n📁 Processing Alaska: {filename}")
    
    try:
        # 1. Read CSV
        print("  📖 Reading CSV...")
        df = pd.read_csv(csv_path)
        
        # Clean columns
        for col in ['.geo', 'system:index']:
            if col in df.columns:
                df = df.drop(columns=[col])
        
        print(f"    Points: {len(df):,}")
        
        # 2. **CONVERT TO 0-360° for Alaska**
        print("  🔄 Converting to 0-360° coordinate space...")
        df['lon_360'] = df['longitude'].where(df['longitude'] >= 0, df['longitude'] + 360)
        
        # 3. Get bounds in 0-360°
        lat_min, lat_max = df['latitude'].min(), df['latitude'].max()
        lon_min_360, lon_max_360 = df['lon_360'].min(), df['lon_360'].max()
        
        print(f"    Alaska bounds (0-360°):")
        print(f"      Lat: [{lat_min:.3f}, {lat_max:.3f}]")
        print(f"      Lon: [{lon_min_360:.3f}, {lon_max_360:.3f}]")
        
        # 4. Create grid in 0-360°
        lat_min_grid = np.floor(lat_min / MODIS_RES) * MODIS_RES
        lat_max_grid = np.ceil(lat_max / MODIS_RES) * MODIS_RES
        lon_min_grid_360 = np.floor(lon_min_360 / MODIS_RES) * MODIS_RES
        lon_max_grid_360 = np.ceil(lon_max_360 / MODIS_RES) * MODIS_RES
        
        height = int(np.round((lat_max_grid - lat_min_grid) / MODIS_RES))
        width = int(np.round((lon_max_grid_360 - lon_min_grid_360) / MODIS_RES))
        
        print(f"  📐 Alaska grid: {height} rows × {width} columns")
        
        # 5. **CRITICAL FIX: Transform with POSITIVE y-resolution**
        transform = from_origin(lon_min_grid_360, lat_max_grid, MODIS_RES, MODIS_RES)
        
        # 6. Create binary mask
        mask = np.full((height, width), -9999, dtype=np.float32)
        
        print("  🗺️  Creating mask in 0-360° space...")
        
        # Calculate array indices
        lat_idx = ((lat_max_grid - df['latitude'].values) / MODIS_RES).astype(int)
        lon_idx = ((df['lon_360'].values - lon_min_grid_360) / MODIS_RES).astype(int)
        
        # Filter valid indices
        valid_mask = (lat_idx >= 0) & (lat_idx < height) & (lon_idx >= 0) & (lon_idx < width)
        
        if valid_mask.any():
            lc_values = df['LC_Type1'].values[valid_mask]
            is_target = np.isin(lc_values, TARGET_CLASSES)
            mask_values = np.where(is_target, 1.0, 0.0)
            mask[lat_idx[valid_mask], lon_idx[valid_mask]] = mask_values
        
        # 7. **VERIFICATION**
        print("  🔍 Verifying Alaska transform...")
        top_left_lon, top_left_lat = transform * (0, 0)
        bottom_right_lon, bottom_right_lat = transform * (width, height)
        
        print(f"    Top-left: ({top_left_lon:.3f}°, {top_left_lat:.3f}°)")
        print(f"    Bottom-right: ({bottom_right_lon:.3f}°, {bottom_right_lat:.3f}°)")
        print(f"    Transform.e coefficient: {transform.e:.6f}")
        
        # 8. Write GeoTIFF with 0-360° CRS
        print("  💾 Writing Alaska GeoTIFF...")
        with rasterio.open(
            output_tif_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=1,
            dtype='float32',
            crs=CRS_360_PROJ4,
            transform=transform,
            nodata=-9999.0,
            compress='lzw'
        ) as dst:
            dst.write(mask, 1)
            
            dst.update_tags(
                Source_File=filename,
                Region='Alaska',
                Coordinate_System='Longitude_0_360',
                Processing='Step1_Alaska_FIXED'
            )
        
        # 9. **FINAL VERIFICATION**
        print("  ✅ Final verification...")
        with rasterio.open(output_tif_path) as verify:
            bounds = verify.bounds
            print(f"    Bounds (0-360°):")
            print(f"      West:  {bounds.left:.3f}°")
            print(f"      East:  {bounds.right:.3f}°")
            print(f"      South: {bounds.bottom:.3f}°")
            print(f"      North: {bounds.top:.3f}°")
            
            # Critical checks
            if bounds.bottom < bounds.top:
                print(f"    ✅ Correct: South < North")
            else:
                print(f"    ❌ ERROR: South > North")
            
            # Alaska should be in 0-360°
            if bounds.right > 180:
                print(f"    ✅ 0-360° coordinate system (correct for Alaska)")
            else:
                print(f"    ⚠️  -180+180° system")
        
        # Statistics
        valid_pixels = (mask != -9999).sum()
        target_pixels = (mask == 1.0).sum()
        
        print(f"  📊 Statistics:")
        print(f"    Valid pixels: {valid_pixels:,}")
        print(f"    Target land cover types: {target_pixels:,}")
        print(f"    Target MODIS classes: {TARGET_CLASSES}")
        
        file_size_mb = output_tif_path.stat().st_size / (1024 * 1024)
        print(f"  ✅ Saved: {output_tif_path.name} ({file_size_mb:.1f} MB)")
        
        return {
            'success': True,
            'year': Path(csv_path).name.split('_')[2],
            'bounds_360': bounds
        }
        
    except Exception as e:
        print(f"  ❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return {'success': False}

# ============================================================================
# BATCH PROCESS ALASKA
# ============================================================================

def process_alaska():
    """
    Process all Alaska files.
    """
    print("="*70)
    print("STEP 1B: PROCESS ALASKA MODIS FILES")
    print("="*70)
    print("FIX 1: transform uses MODIS_RES (positive) in from_origin()")
    print("FIX 2: Processed in 0-360° space to avoid dateline wrap")
    print(f"TARGET_CLASSES (based on 57 final figure sites): {TARGET_CLASSES}")
    print("-" * 70)
    
    # Get Alaska CSV files
    csv_files = [f for f in os.listdir(ALASKA_INPUT_PATH) 
                 if f.endswith('.csv')]
    
    if not csv_files:
        print("❌ No Alaska CSV files found!")
        return False
    
    print(f"Found {len(csv_files)} Alaska files")
    print("Sample files:", csv_files[:3])
    
    results = []
    
    # Process all files
    for csv_file in tqdm(sorted(csv_files), desc="Processing Alaska"):
        # Extract year
        try:
            year = csv_file.split('_')[2]  # Adjust based on your filename
        except:
            year = "unknown"
        
        # Create output path
        input_path = os.path.join(ALASKA_INPUT_PATH, csv_file)
        output_file = ALASKA_OUTPUT_DIR / f"alaska_mask_{year}_500m_360.tif"
        
        # Process file
        result = create_alaska_mask(input_path, output_file)
        results.append(result)
    
    # Summary
    print("\n" + "="*70)
    print("STEP 1B COMPLETE - ALASKA SUMMARY")
    print("="*70)
    
    successful = [r for r in results if r.get('success', False)]
    failed = [r for r in results if not r.get('success', False)]
    
    print(f"✅ Successful: {len(successful)} files")
    print(f"❌ Failed: {len(failed)} files")
    
    if successful:
        print(f"\n📁 Output location: {ALASKA_OUTPUT_DIR}")
        print(f"📊 MODIS classes retained: {TARGET_CLASSES}")
        print(f"   (Filtered to match 57 sites in final figure)")
        
        # Show first file bounds for verification
        first_file = ALASKA_OUTPUT_DIR / f"alaska_mask_{successful[0]['year']}_500m_360.tif"
        if first_file.exists():
            with rasterio.open(first_file) as src:
                b = src.bounds
                print(f"\n📋 First file verification:")
                print(f"  West (0-360°):  {b.left:.3f}°")
                print(f"  East (0-360°):  {b.right:.3f}°")
                print(f"  South: {b.bottom:.3f}°")
                print(f"  North: {b.top:.3f}°")
                print(f"  South < North: {b.bottom < b.top}")
                print(f"  Right > 180° (0-360°): {b.right > 180}")
    
    return len(successful) > 0

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("PROCESSING ALASKA MODIS FILES")
    print("-" * 70)
    print("This script creates 500m binary masks for Alaska only")
    print("Output: alaska_mask_YYYY_500m_360.tif files (0-360° coordinates)")
    print("=" * 70)
    
    success = process_alaska()
    
    if success:
        print("\n" + "="*70)
        print("✅ ALASKA PROCESSING COMPLETE")
        print("="*70)
        print("Next: Verify outputs, then proceed to drought extraction")
    else:
        print("\n❌ Alaska processing failed.")