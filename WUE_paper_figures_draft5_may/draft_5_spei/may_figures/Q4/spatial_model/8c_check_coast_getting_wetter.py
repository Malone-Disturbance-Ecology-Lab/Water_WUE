"""
EXTRACT SPEI-48 REGIONAL TIME SERIES (IDENTICAL STRUCTURE TO WUE_T)
Creates all_regions_monthly_spei48.csv matching the format expected by 
the SPEI diagnostic analysis.

Uses the SAME regional bounds as your existing WUE_T extraction.
"""

import xarray as xr
import numpy as np
import pandas as pd
import os
import glob
import re
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION - SAME AS YOUR EXISTING EXTRACTION
# ============================================================================

# Input: SPEI NetCDF files (same source used in spatial emulation)
# These should be the original SPEI-48 files, NOT the prediction outputs
SPEI_INPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\extraction_final_fixed"

# Output directory (same as where all_regions_monthly_frac_less.csv lives)
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\regional_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# SPEI variable name in your NetCDF files
SPEI_VAR_NAME = "SPEI48"  # or "spei48" or "SPEI_48" - check your files

# Region bounds - IDENTICAL to your WUE_T extraction
COASTAL_REGIONS = {
    "Pacific": {"lon_bounds": [-130, -116], "lat_bounds": [32, 49]},
    "Gulf": {"lon_bounds": [-98, -80], "lat_bounds": [24, 31]},
    "Southeast Atlantic": {"lon_bounds": [-82, -75], "lat_bounds": [30, 37]},
    "Atlantic North": {"lon_bounds": [-78, -66], "lat_bounds": [37, 46]},
    "Alaska": {"lon_bounds": [-165, -130], "lat_bounds": [52, 72]}
}

# ============================================================================
# FUNCTION TO CHECK WHAT VARIABLES ARE AVAILABLE
# ============================================================================

def inspect_spei_file():
    """Check SPEI file structure to find correct variable name."""
    files = sorted(glob.glob(os.path.join(SPEI_INPUT_DIR, "*.nc")))
    if not files:
        print(f"❌ No NetCDF files found in {SPEI_INPUT_DIR}")
        return None
    
    test_file = files[0]
    print(f"\n🔍 Inspecting: {os.path.basename(test_file)}")
    
    ds = xr.open_dataset(test_file)
    print(f"\nAvailable variables: {list(ds.data_vars.keys())}")
    print(f"Coordinates: {list(ds.coords.keys())}")
    
    # Check for SPEI variable
    spei_candidates = [v for v in ds.data_vars if 'spei' in v.lower()]
    print(f"\nSPEI candidates: {spei_candidates}")
    
    ds.close()
    return spei_candidates[0] if spei_candidates else None

# ============================================================================
# EXTRACT SPEI-48 FOR ONE REGION
# ============================================================================

def extract_spei_region(region_name, lon_bounds, lat_bounds, spei_var):
    """Extract monthly regional mean SPEI-48 values."""
    
    # Get all NetCDF files (sorted by date)
    nc_files = sorted(glob.glob(os.path.join(SPEI_INPUT_DIR, "*.nc")))
    
    if not nc_files:
        raise FileNotFoundError(f"No NetCDF files found in {SPEI_INPUT_DIR}")
    
    print(f"\n📂 Extracting SPEI-48 for {region_name} from {len(nc_files)} files...")
    
    dates = []
    spei_values = []
    
    for nc_file in nc_files:
        # Extract date from filename
        basename = os.path.basename(nc_file)
        # Look for pattern like _202003_ or 202003
        match = re.search(r'_(\d{6})_', basename)
        if not match:
            match = re.search(r'(\d{6})', basename)
        
        if match:
            yyyymm = match.group(1)
            year = int(yyyymm[:4])
            month = int(yyyymm[4:6])
            date = datetime(year, month, 1)
        else:
            print(f"  ⚠️ Could not parse date from: {basename}")
            continue
        
        try:
            ds = xr.open_dataset(nc_file)
            
            # Check if SPEI variable exists
            if spei_var not in ds.data_vars:
                print(f"  ⚠️ {spei_var} not in {basename}")
                ds.close()
                continue
            
            # Get SPEI data
            spei_data = ds[spei_var].values
            
            # Get coordinates - handle both 1D and 2D grids
            if 'lon' in ds.coords and 'lat' in ds.coords:
                lons = ds.lon.values
                lats = ds.lat.values
                
                # Convert longitude if needed (0-360 to -180-180)
                if np.nanmax(lons) > 180:
                    lons = ((lons + 180) % 360) - 180
                    # Regrid if necessary
                    if spei_data.ndim == 2:
                        # Need to reorder if lon was changed
                        sort_idx = np.argsort(lons)
                        lons = lons[sort_idx]
                        spei_data = spei_data[sort_idx, :] if spei_data.shape[0] == len(sort_idx) else spei_data[:, sort_idx]
                
                # Create 2D grid for masking
                if spei_data.ndim == 2:
                    lon_grid, lat_grid = np.meshgrid(lons, lats)
                    spei_grid = spei_data
                else:
                    # Already 2D
                    lon_grid, lat_grid = lons, lats
                    spei_grid = spei_data
                    
            elif 'x' in ds.coords and 'y' in ds.coords:
                # Alternative coordinate names
                lons = ds.x.values
                lats = ds.y.values
                lon_grid, lat_grid = np.meshgrid(lons, lats)
                spei_grid = spei_data
            else:
                print(f"  ⚠️ Unknown coordinate system in {basename}")
                ds.close()
                continue
            
            # Create region mask
            mask = ((lon_grid >= lon_bounds[0]) & (lon_grid <= lon_bounds[1]) &
                    (lat_grid >= lat_bounds[0]) & (lat_grid <= lat_bounds[1]))
            
            # Extract SPEI values for this region
            region_spei = spei_grid[mask]
            
            # Remove NaN values
            valid_spei = region_spei[np.isfinite(region_spei)]
            
            if len(valid_spei) > 0:
                mean_spei = np.mean(valid_spei)
                spei_values.append(mean_spei)
                dates.append(date)
                
                # Debug first few files
                if len(dates) <= 3:
                    print(f"    {date.strftime('%Y-%m')}: {len(valid_spei)} pixels, SPEI={mean_spei:.3f}")
            else:
                if len(dates) <= 3:
                    print(f"    {date.strftime('%Y-%m')}: No valid SPEI pixels")
            
            ds.close()
            
        except Exception as e:
            print(f"  Error processing {basename}: {e}")
            continue
    
    # Create DataFrame (matches structure of all_regions_monthly_frac_less.csv)
    df = pd.DataFrame({
        'date': dates, 
        'mean_SPEI48': spei_values, 
        'region': region_name
    })
    df = df.sort_values('date').reset_index(drop=True)
    
    print(f"\n  ✅ Extracted {len(df)} months for {region_name}")
    if len(df) > 0:
        print(f"     Date range: {df['date'].min().strftime('%Y-%m')} to {df['date'].max().strftime('%Y-%m')}")
        print(f"     Mean SPEI-48: {df['mean_SPEI48'].mean():.3f}")
    
    return df

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def extract_all_spei_regions():
    """Extract SPEI-48 time series for all coastal regions."""
    
    print("="*80)
    print("EXTRACT SPEI-48 REGIONAL TIME SERIES")
    print("="*80)
    print(f"Input directory: {SPEI_INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("-"*80)
    
    # First, check what variable name to use
    spei_var = inspect_spei_file()
    if spei_var is None:
        print("\n❌ Could not identify SPEI variable in files.")
        print("Please check the SPEI NetCDF files and update SPEI_VAR_NAME")
        return None
    
    print(f"\n✅ Using SPEI variable: {spei_var}")
    
    # Extract for each region
    all_results = []
    
    for region_name, bounds in COASTAL_REGIONS.items():
        df = extract_spei_region(
            region_name,
            bounds['lon_bounds'],
            bounds['lat_bounds'],
            spei_var
        )
        
        if len(df) > 0:
            # Save individual region file (optional)
            individual_file = os.path.join(OUTPUT_DIR, f"{region_name}_monthly_spei48.csv")
            df.to_csv(individual_file, index=False)
            print(f"  💾 Saved: {individual_file}")
            all_results.append(df)
        else:
            print(f"  ⚠️ No data for {region_name}")
    
    # Create combined file (THIS IS WHAT YOU NEED)
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        output_file = os.path.join(OUTPUT_DIR, "all_regions_monthly_spei48.csv")
        combined_df.to_csv(output_file, index=False)
        
        print("\n" + "="*80)
        print("✅ SUCCESS!")
        print("="*80)
        print(f"Combined file saved: {output_file}")
        print(f"\nFile structure:")
        print(f"  - Columns: date, mean_SPEI48, region")
        print(f"  - Total records: {len(combined_df)}")
        print(f"  - Regions: {combined_df['region'].unique().tolist()}")
        print(f"  - Date range: {combined_df['date'].min()} to {combined_df['date'].max()}")
        
        # Quick summary
        print("\n📊 REGIONAL SPEI-48 SUMMARY:")
        for region in combined_df['region'].unique():
            region_df = combined_df[combined_df['region'] == region]
            print(f"  {region}: {len(region_df)} months, mean SPEI = {region_df['mean_SPEI48'].mean():.3f}")
        
        return combined_df
    else:
        print("\n❌ No data extracted for any region!")
        return None

# ============================================================================
# VERIFICATION FUNCTION
# ============================================================================

def verify_output():
    """Verify that the output file matches expected format."""
    output_file = os.path.join(OUTPUT_DIR, "all_regions_monthly_spei48.csv")
    
    if not os.path.exists(output_file):
        print(f"❌ File not found: {output_file}")
        return False
    
    df = pd.read_csv(output_file)
    df['date'] = pd.to_datetime(df['date'])
    
    print("\n" + "="*80)
    print("VERIFICATION")
    print("="*80)
    print(f"✅ File exists: {output_file}")
    print(f"   Shape: {df.shape}")
    print(f"   Columns: {df.columns.tolist()}")
    print(f"   Date range: {df['date'].min()} to {df['date'].max()}")
    print(f"   Unique regions: {df['region'].unique().tolist()}")
    
    # Check for missing values
    missing_spei = df['mean_SPEI48'].isna().sum()
    if missing_spei > 0:
        print(f"   ⚠️ Warning: {missing_spei} missing SPEI values")
    else:
        print(f"   ✅ No missing SPEI values")
    
    return True

if __name__ == "__main__":
    # Extract SPEI-48 data
    df = extract_all_spei_regions()
    
    # Verify output
    if df is not None:
        verify_output()
        
        print("\n" + "="*80)
        print("NEXT STEP")
        print("="*80)
        print("Now run the SPEI-48 diagnostic analysis script.")
        print("It will automatically find: all_regions_monthly_spei48.csv")