"""
STEP 1: Extract Regional Time Series from Monthly NetCDF Files (CORRECTED)
Uses Class_WUE_tra from netcdf_outputs folder.
Includes proper debugging to verify region overlap.
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
# CONFIGURATION
# ============================================================================

NETCDF_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs"
OUTPUT_DIR = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\regional_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

VAR_NAME = "Class_WUE_tra"  # 1 = Less efficient (Decrease), 0 = Increase/NoChange
FILL_VALUE = -9999

# Region bounds (must match the grid of your NetCDF files)
COASTAL_REGIONS = {
    "Pacific": {"lon_bounds": [-130, -116], "lat_bounds": [32, 49], "color": "#E67E22"},
    "Gulf": {"lon_bounds": [-98, -80], "lat_bounds": [24, 31], "color": "#F1C40F"},
    "Southeast Atlantic": {"lon_bounds": [-82, -75], "lat_bounds": [30, 37], "color": "#9B59B6"},
    "Atlantic North": {"lon_bounds": [-78, -66], "lat_bounds": [37, 46], "color": "#3498DB"},
    "Alaska": {"lon_bounds": [-165, -130], "lat_bounds": [52, 72], "color": "#1ABC9C"}
}

# ============================================================================
# FUNCTION TO CHECK GRID OVERLAP
# ============================================================================

def check_region_overlap(region_name, lon_bounds, lat_bounds, test_file):
    """Debug function to check if region bounds overlap with grid."""
    ds = xr.open_dataset(test_file)
    
    # Convert longitude if needed
    if np.nanmax(ds.lon.values) > 180:
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180).sortby('lon')
    
    lons = ds.lon.values
    lats = ds.lat.values
    
    print(f"\n🔍 Checking {region_name} bounds against grid:")
    print(f"  Region lon: [{lon_bounds[0]}, {lon_bounds[1]}]")
    print(f"  Grid lon range: [{lons.min():.2f}, {lons.max():.2f}]")
    print(f"  Region lat: [{lat_bounds[0]}, {lat_bounds[1]}]")
    print(f"  Grid lat range: [{lats.min():.2f}, {lats.max():.2f}]")
    
    # Check overlap
    lon_overlap = (lon_bounds[0] <= lons.max()) and (lon_bounds[1] >= lons.min())
    lat_overlap = (lat_bounds[0] <= lats.max()) and (lat_bounds[1] >= lats.min())
    
    print(f"  Lon overlap: {lon_overlap}")
    print(f"  Lat overlap: {lat_overlap}")
    
    ds.close()
    return lon_overlap and lat_overlap

# ============================================================================
# MAIN EXTRACTION FUNCTION
# ============================================================================

def extract_region_time_series(region_name, lon_bounds, lat_bounds):
    """Extract monthly regional mean frac_less from Class_WUE_tra."""
    
    # Get all NetCDF files
    nc_files = sorted(glob.glob(os.path.join(NETCDF_DIR, "WUE_predictions_linearlogistic_*.nc")))
    
    if not nc_files:
        raise FileNotFoundError(f"No NetCDF files found in {NETCDF_DIR}")
    
    # Check overlap using first file
    check_region_overlap(region_name, lon_bounds, lat_bounds, nc_files[0])
    
    print(f"\n📂 Extracting {region_name} from {len(nc_files)} files...")
    
    dates = []
    region_frac_less = []
    
    for nc_file in nc_files:
        # Extract date
        basename = os.path.basename(nc_file)
        match = re.search(r'_(\d{6})\.nc', basename)
        if match:
            yyyymm = match.group(1)
            year = int(yyyymm[:4])
            month = int(yyyymm[4:6])
            date = datetime(year, month, 1)
        else:
            continue
        
        try:
            ds = xr.open_dataset(nc_file)
            
            # Convert longitude from 0-360 to -180-180 if needed
            if np.nanmax(ds.lon.values) > 180:
                ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180).sortby('lon')
            
            # Check if variable exists
            if VAR_NAME not in ds.data_vars:
                print(f"  Warning: {VAR_NAME} not in {basename}")
                ds.close()
                continue
            
            # Get class data
            class_data = ds[VAR_NAME].values
            
            # Create coordinate grids
            lons = ds.lon.values
            lats = ds.lat.values
            if lons.ndim == 1 and lats.ndim == 1:
                lon_grid, lat_grid = np.meshgrid(lons, lats)
            else:
                lon_grid, lat_grid = lons, lats
            
            # Create region mask
            mask = ((lon_grid >= lon_bounds[0]) & (lon_grid <= lon_bounds[1]) &
                    (lat_grid >= lat_bounds[0]) & (lat_grid <= lat_bounds[1]))
            
            # Extract region data with proper valid mask
            region_data = class_data[mask]
            # Remove fill values and NaN
            valid_data = region_data[(region_data != FILL_VALUE) & np.isfinite(region_data)]
            
            if len(valid_data) > 0:
                # frac_less = proportion of pixels with Class=1 (Less efficient)
                frac_less = np.mean(valid_data)
                region_frac_less.append(frac_less)
                dates.append(date)
                # Debug first few files
                if len(dates) <= 3:
                    print(f"    {date.strftime('%Y-%m')}: {len(valid_data)} valid pixels, frac_less={frac_less:.3f}")
            else:
                # No valid data in this region for this month
                if len(dates) <= 3:
                    print(f"    {date.strftime('%Y-%m')}: No valid pixels (mask sum={mask.sum()})")
            
            ds.close()
            
        except Exception as e:
            print(f"  Error processing {basename}: {e}")
            continue
    
    # Create DataFrame
    df = pd.DataFrame({'date': dates, 'frac_less': region_frac_less, 'region': region_name})
    df = df.sort_values('date').reset_index(drop=True)
    
    print(f"\n  Extracted {len(df)} months of data for {region_name}")
    if len(df) > 0:
        print(f"  Date range: {df['date'].min().strftime('%Y-%m')} to {df['date'].max().strftime('%Y-%m')}")
        print(f"  Mean frac_less: {df['frac_less'].mean():.3f}")
    
    return df

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def extract_all_regions():
    """Extract time series for all coastal regions."""
    all_results = []
    
    for region_name, info in COASTAL_REGIONS.items():
        df = extract_region_time_series(
            region_name, 
            info['lon_bounds'], 
            info['lat_bounds']
        )
        
        if len(df) > 0:
            output_file = os.path.join(OUTPUT_DIR, f"{region_name}_monthly_frac_less.csv")
            df.to_csv(output_file, index=False)
            print(f"  Saved: {output_file}")
            all_results.append(df)
        else:
            print(f"  ⚠️ No data for {region_name}")
    
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_file = os.path.join(OUTPUT_DIR, "all_regions_monthly_frac_less.csv")
        combined_df.to_csv(combined_file, index=False)
        print(f"\n✅ Combined data saved: {combined_file}")
        return combined_df
    else:
        print("\n❌ No data extracted for any region!")
        return None

if __name__ == "__main__":
    print("="*80)
    print("EXTRACT REGIONAL TIME SERIES FROM MONTHLY NETCDF FILES")
    print("="*80)
    print(f"Source: {NETCDF_DIR}")
    print(f"Variable: {VAR_NAME}")
    print(f"  Class=1 = Less efficient (Decrease in WUE_T)")
    print("-"*80)
    
    df = extract_all_regions()
    
    if df is not None:
        print(f"\n📊 SUMMARY:")
        print(f"  Total records: {len(df)}")
        print(f"  Regions: {df['region'].unique().tolist()}")
        print(f"  Date range: {df['date'].min()} to {df['date'].max()}")
        
        print("\n📈 REGIONAL STATISTICS:")
        for region in df['region'].unique():
            region_df = df[df['region'] == region]
            print(f"  {region}: {len(region_df)} months, mean frac_less = {region_df['frac_less'].mean():.3f}")