"""
STEP 1A: Process CONUS (non-Alaska) MODIS CSV to 500m Binary Masks - FIXED
File: step1_conus_FIXED.py
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
# CONFIGURATION - CONUS (NON-ALASKA)
# ============================================================================

# MODIS Land Cover Classes for FINAL FIGURE SITES (57 sites in strict triple intersection)
# Based on actual IGBP biomes from the 57 sites: 
# ENF, DBF, MF, CSH, OSH, GRA, WET, CRO, BSV
TARGET_CLASSES = [1, 4, 5, 6, 7, 10, 11, 12, 16]
# Class details for reference:
#   1  = Evergreen Needleleaf Forests (ENF) - e.g., US-HB2, US-HB3, US-MRf, US-ONA, US-xSB
#   4  = Deciduous Broadleaf Forests (DBF) - e.g., US-Slt, US-xSE
#   5  = Mixed Forests (MF) - e.g., US-Dix
#   6  = Closed Shrublands (CSH) - e.g., US-Ced, US-KS2, US-SO2, US-SO3, US-SO4
#   7  = Open Shrublands (OSH) - e.g., US-SCs
#   10 = Grasslands (GRA) - e.g., US-BWb, US-CGG, US-NGC, US-PAS, US-SCg, US-Snd
#   11 = Permanent Wetlands (WET) - e.g., Most of the 57 sites (35+ sites)
#   12 = Croplands (CRO) - e.g., US-Tw3, US-Twt
#   16 = Barren (BSV) - e.g., US-A03, US-A10

# Paths (SAME - DO NOT CHANGE)
CONUS_INPUT_PATH = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\MODIS_landcover\non_alaska"
CONUS_OUTPUT_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\modis_landcover_shape\processed\conus_masks_500m_FIXED")

# Resolution
MODIS_RES = 0.0045  # 500m at mid-latitudes

# CRS for CONUS (-180 to +180°)
CRS_PROJ4 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

CONUS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# FUNCTION: Create CONUS Mask - FIXED
# ============================================================================

def create_conus_mask(csv_path, output_tif_path):
    """
    Create CONUS binary mask with CORRECT transform.
    """
    filename = Path(csv_path).name
    print(f"\n📁 Processing CONUS: {filename}")
    
    try:
        # 1. Read CSV
        print("  📖 Reading CSV...")
        df = pd.read_csv(csv_path)
        
        # Clean columns
        for col in ['.geo', 'system:index']:
            if col in df.columns:
                df = df.drop(columns=[col])
        
        print(f"    Points: {len(df):,}")
        
        # 2. Get bounds
        lat_min, lat_max = df['latitude'].min(), df['latitude'].max()
        lon_min, lon_max = df['longitude'].min(), df['longitude'].max()
        
        print(f"    Raw bounds: Lat [{lat_min:.3f}, {lat_max:.3f}]")
        print(f"                Lon [{lon_min:.3f}, {lon_max:.3f}]")
        
        # 3. Create grid
        lat_min_grid = np.floor(lat_min / MODIS_RES) * MODIS_RES
        lat_max_grid = np.ceil(lat_max / MODIS_RES) * MODIS_RES
        lon_min_grid = np.floor(lon_min / MODIS_RES) * MODIS_RES
        lon_max_grid = np.ceil(lon_max / MODIS_RES) * MODIS_RES
        
        height = int(np.round((lat_max_grid - lat_min_grid) / MODIS_RES))
        width = int(np.round((lon_max_grid - lon_min_grid) / MODIS_RES))
        
        print(f"  📐 Grid: {height} rows × {width} columns")
        
        # 4. **CRITICAL FIX: Transform with POSITIVE y-resolution**
        transform = from_origin(lon_min_grid, lat_max_grid, MODIS_RES, MODIS_RES)
        
        # 5. Create binary mask
        mask = np.full((height, width), -9999, dtype=np.float32)
        
        print("  🗺️  Creating mask...")
        
        # Calculate array indices
        lat_idx = ((lat_max_grid - df['latitude'].values) / MODIS_RES).astype(int)
        lon_idx = ((df['longitude'].values - lon_min_grid) / MODIS_RES).astype(int)
        
        # Filter valid indices
        valid_mask = (lat_idx >= 0) & (lat_idx < height) & (lon_idx >= 0) & (lon_idx < width)
        
        if valid_mask.any():
            lc_values = df['LC_Type1'].values[valid_mask]
            is_target = np.isin(lc_values, TARGET_CLASSES)
            mask_values = np.where(is_target, 1.0, 0.0)
            mask[lat_idx[valid_mask], lon_idx[valid_mask]] = mask_values
        
        # 6. **VERIFICATION BEFORE WRITING**
        print("  🔍 Verifying transform...")
        top_left_lon, top_left_lat = transform * (0, 0)
        bottom_right_lon, bottom_right_lat = transform * (width, height)
        
        print(f"    Top-left: ({top_left_lon:.3f}°, {top_left_lat:.3f}°)")
        print(f"    Bottom-right: ({bottom_right_lon:.3f}°, {bottom_right_lat:.3f}°)")
        print(f"    Transform.e coefficient: {transform.e:.6f}")
        
        # 7. Write GeoTIFF
        print("  💾 Writing GeoTIFF...")
        with rasterio.open(
            output_tif_path,
            'w',
            driver='GTiff',
            height=height,
            width=width,
            count=1,
            dtype='float32',
            crs=CRS_PROJ4,
            transform=transform,
            nodata=-9999.0,
            compress='lzw'
        ) as dst:
            dst.write(mask, 1)
            
            dst.update_tags(
                Source_File=filename,
                Region='CONUS',
                Processing='Step1_CONUS_FIXED'
            )
        
        # 8. **FINAL VERIFICATION**
        print("  ✅ Final verification...")
        with rasterio.open(output_tif_path) as verify:
            bounds = verify.bounds
            print(f"    Bounds:")
            print(f"      West:  {bounds.left:.3f}°")
            print(f"      East:  {bounds.right:.3f}°")
            print(f"      South: {bounds.bottom:.3f}°")
            print(f"      North: {bounds.top:.3f}°")
            
            # Critical checks
            if bounds.bottom < bounds.top:
                print(f"    ✅ Correct: South < North")
            else:
                print(f"    ❌ ERROR: South > North")
            
            if bounds.right > 180:
                print(f"    ⚠️  0-360° detected (unexpected for CONUS)")
            else:
                print(f"    ✅ -180+180° coordinate system")
        
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
            'bounds': bounds
        }
        
    except Exception as e:
        print(f"  ❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return {'success': False}

# ============================================================================
# BATCH PROCESS CONUS
# ============================================================================

def process_conus():
    """
    Process all CONUS files.
    """
    print("="*70)
    print("STEP 1A: PROCESS CONUS (NON-ALASKA) MODIS FILES")
    print("="*70)
    print("FIX: transform uses MODIS_RES (positive) in from_origin()")
    print(f"TARGET_CLASSES (based on 57 final figure sites): {TARGET_CLASSES}")
    print("-" * 70)
    
    # Get CONUS CSV files
    csv_files = [f for f in os.listdir(CONUS_INPUT_PATH) 
                 if f.endswith('.csv')]
    
    if not csv_files:
        print("❌ No CONUS CSV files found!")
        return False
    
    print(f"Found {len(csv_files)} CONUS files")
    print("Sample files:", csv_files[:3])
    
    results = []
    
    # Process all files
    for csv_file in tqdm(sorted(csv_files), desc="Processing CONUS"):
        # Extract year
        try:
            year = csv_file.split('_')[2]  # Adjust based on your filename
        except:
            year = "unknown"
        
        # Create output path
        input_path = os.path.join(CONUS_INPUT_PATH, csv_file)
        output_file = CONUS_OUTPUT_DIR / f"conus_mask_{year}_500m.tif"
        
        # Process file
        result = create_conus_mask(input_path, output_file)
        results.append(result)
    
    # Summary
    print("\n" + "="*70)
    print("STEP 1A COMPLETE - CONUS SUMMARY")
    print("="*70)
    
    successful = [r for r in results if r.get('success', False)]
    failed = [r for r in results if not r.get('success', False)]
    
    print(f"✅ Successful: {len(successful)} files")
    print(f"❌ Failed: {len(failed)} files")
    
    if successful:
        print(f"\n📁 Output location: {CONUS_OUTPUT_DIR}")
        print(f"📊 MODIS classes retained: {TARGET_CLASSES}")
        print(f"   (Filtered to match 57 sites in final figure)")
        
        # Show first file bounds for verification
        first_file = CONUS_OUTPUT_DIR / f"conus_mask_{successful[0]['year']}_500m.tif"
        if first_file.exists():
            with rasterio.open(first_file) as src:
                b = src.bounds
                print(f"\n📋 First file verification:")
                print(f"  West:  {b.left:.3f}°")
                print(f"  East:  {b.right:.3f}°")
                print(f"  South: {b.bottom:.3f}°")
                print(f"  North: {b.top:.3f}°")
                print(f"  South < North: {b.bottom < b.top}")
    
    return len(successful) > 0

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("PROCESSING CONUS (NON-ALASKA) MODIS FILES")
    print("-" * 70)
    print("This script creates 500m binary masks for CONUS only")
    print("Output: conus_mask_YYYY_500m.tif files")
    print("=" * 70)
    
    success = process_conus()
    
    if success:
        print("\n" + "="*70)
        print("✅ CONUS PROCESSING COMPLETE")
        print("="*70)
        print("Next: Run Alaska script (step2_alaska_FIXED.py)")
    else:
        print("\n❌ CONUS processing failed.")