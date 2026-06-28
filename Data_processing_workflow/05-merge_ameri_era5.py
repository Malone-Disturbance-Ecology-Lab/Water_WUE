# -*- coding: utf-8 -*-
"""
Created for merging AmeriFlux data with ERA5 point data
Includes timezone conversion from UTC to local time with proper error handling
"""

import os
import pandas as pd
import zipfile
import numpy as np
from datetime import datetime, timedelta
import re
import pytz
from timezonefinder import TimezoneFinder

def get_site_timezone(latitude, longitude):
    """
    Get timezone for a site based on coordinates
    Returns timezone string and UTC offset in hours
    
    Parameters:
        latitude (float): Site latitude
        longitude (float): Site longitude
    
    Returns:
        tuple: (timezone_str, offset_hours) or (None, None) if failed
    """
    try:
        tf = TimezoneFinder()
        timezone_str = tf.timezone_at(lng=longitude, lat=latitude)
        
        if timezone_str is None:
            print(f"    Warning: Could not find timezone for coordinates ({latitude}, {longitude})")
            return None, None
        
        # Get timezone object
        local_zone = pytz.timezone(timezone_str)
        
        # Get UTC offset for standard time (no DST)
        # Use a date in January to avoid DST
        jan_1 = datetime(2020, 1, 1)
        utc_time = pytz.UTC.localize(jan_1)
        local_time = utc_time.astimezone(local_zone)
        
        # Calculate offset in hours
        offset_hours = local_time.utcoffset().total_seconds() / 3600
        
        print(f"    Timezone detected: {timezone_str}, UTC offset: {offset_hours} hours")
        return timezone_str, offset_hours
        
    except Exception as e:
        print(f"    Error getting timezone: {e}")
        return None, None

def convert_utc_to_local_stable(df, latitude_col='latitude', longitude_col='longitude', time_col='valid_time'):
    """
    Convert UTC timestamps to local standard time based on site coordinates
    Uses a single timezone lookup per site for efficiency
    
    Parameters:
        df (pd.DataFrame): DataFrame with UTC timestamps and coordinates
        latitude_col (str): Name of latitude column
        longitude_col (str): Name of longitude column
        time_col (str): Name of timestamp column
    
    Returns:
        pd.Series: Series with local timestamps, or None if conversion fails
    """
    # Get unique coordinates (should be same for entire site)
    unique_lats = df[latitude_col].unique()
    unique_lons = df[longitude_col].unique()
    
    if len(unique_lats) == 0 or len(unique_lons) == 0:
        print(f"    Error: No latitude/longitude data found")
        return None
    
    # Use the first valid coordinate pair
    lat = unique_lats[0]
    lon = unique_lons[0]
    
    print(f"    Site coordinates: lat={lat}, lon={lon}")
    
    # Get timezone and offset for this site
    timezone_str, offset_hours = get_site_timezone(lat, lon)
    
    if timezone_str is None or offset_hours is None:
        print(f"    Error: Could not determine timezone for site")
        return None
    
    print(f"    Converting UTC to local time (offset: {offset_hours} hours)")
    
    # Convert all timestamps using the fixed offset
    def convert_single_time(utc_time):
        try:
            if not isinstance(utc_time, pd.Timestamp):
                utc_time = pd.to_datetime(utc_time)
            
            if utc_time.tzinfo is None:
                utc_time = utc_time.tz_localize('UTC')
            
            # Apply fixed offset to get local standard time
            local_time = utc_time + timedelta(hours=offset_hours)
            return local_time.replace(tzinfo=None)
            
        except Exception as e:
            print(f"    Error converting timestamp {utc_time}: {e}")
            return pd.NaT
    
    # Apply conversion to all timestamps
    local_times = df[time_col].apply(convert_single_time)
    
    # Check if conversion was successful
    if local_times.isna().all():
        print(f"    Error: All timestamp conversions failed")
        return None
    
    return local_times

def read_era5_from_zip(zip_path):
    """
    Read ERA5 data from zip file containing CSV with proper error handling
    
    Parameters:
        zip_path (str): Path to the zip file
    
    Returns:
        pd.DataFrame: ERA5 data with datetime index, or None if failed
    """
    try:
        # Open zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Find the CSV file inside the zip
            csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv')]
            if not csv_files:
                print(f"    Error: No CSV file found in {zip_path}")
                return None
            
            # Read the first CSV file
            with zip_ref.open(csv_files[0]) as csv_file:
                era5_df = pd.read_csv(csv_file)
        
        print(f"    ERA5 columns: {list(era5_df.columns)}")
        
        # Check required columns
        required_cols = ['valid_time', 'latitude', 'longitude']
        missing_cols = [col for col in required_cols if col not in era5_df.columns]
        if missing_cols:
            print(f"    Error: Missing required columns: {missing_cols}")
            return None
        
        # Parse datetime
        era5_df['TIMESTAMP_UTC'] = pd.to_datetime(era5_df['valid_time'])
        
        # Check if we have valid timestamps
        if era5_df['TIMESTAMP_UTC'].isna().all():
            print(f"    Error: No valid timestamps found in ERA5 data")
            return None
        
        # Convert UTC to local time
        print(f"    Converting ERA5 timestamps from UTC to local time...")
        local_times = convert_utc_to_local_stable(era5_df, 'latitude', 'longitude', 'TIMESTAMP_UTC')
        
        if local_times is None:
            print(f"    Error: Failed to convert timestamps to local time")
            return None
        
        era5_df['TIMESTAMP'] = local_times
        
        # Remove any rows with invalid timestamps
        initial_rows = len(era5_df)
        era5_df = era5_df.dropna(subset=['TIMESTAMP'])
        if len(era5_df) < initial_rows:
            print(f"    Removed {initial_rows - len(era5_df)} rows with invalid timestamps")
        
        if era5_df.empty:
            print(f"    Error: No valid timestamps after conversion")
            return None
        
        # Drop the UTC timestamp column
        era5_df = era5_df.drop(columns=['TIMESTAMP_UTC'])
        
        # Set timestamp as index
        era5_df.set_index('TIMESTAMP', inplace=True)
        
        # Resample to 30-minute frequency and interpolate - better approach
        # Resample the entire dataframe, not just numeric columns
        if len(era5_df) < 2:
            print(f"    Warning: Too few data points for resampling ({len(era5_df)} rows)")
        else:
            # Resample to 30-minute frequency
            era5_df = era5_df.resample('30T').interpolate(method='linear')
            print(f"    Resampled to 30-minute frequency")
        
        # Reset index to have TIMESTAMP as column again
        era5_df.reset_index(inplace=True)
        
        # Print time range for verification
        time_min = era5_df['TIMESTAMP'].min()
        time_max = era5_df['TIMESTAMP'].max()
        print(f"    ERA5 data loaded: {len(era5_df)} rows")
        print(f"    Time range (local): {time_min} to {time_max}")
        
        return era5_df
        
    except Exception as e:
        print(f"    Error reading ERA5 zip file: {e}")
        import traceback
        traceback.print_exc()
        return None

def merge_ameri_era5(ameri_path, era5_path, output_path):
    """
    Merges all AmeriFlux data with corresponding ERA5 data
    
    Parameters:
        ameri_path (str): Path to the AmeriFlux data directory (with gaps_*.csv files)
        era5_path (str): Path to the ERA5 point data directory (with zip files)
        output_path (str): Path to save merged data
    
    Returns:
        dict: Dictionary with site names as keys and merged DataFrames as values
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)
    
    # Get all AmeriFlux files
    ameri_files = [f for f in os.listdir(ameri_path) 
                   if f.startswith('gaps_') and f.endswith('.csv')]
    
    if not ameri_files:
        print(f"No AmeriFlux files found in {ameri_path}")
        return {}
    
    print(f"Found {len(ameri_files)} AmeriFlux files to process")
    print(f"{'='*60}")
    
    merged_results = {}
    successful = 0
    failed = 0
    
    for ameri_file in ameri_files:
        # Extract site name from filename (gaps_US-A03.csv -> US-A03)
        site_name = ameri_file.replace('gaps_', '').replace('.csv', '')
        print(f"\n{'='*50}")
        print(f"Processing site: {site_name}")
        print(f"{'='*50}")
        
        # Find corresponding ERA5 zip file
        era5_zip = None
        for zip_file in os.listdir(era5_path):
            if zip_file.endswith('.zip') and site_name in zip_file:
                era5_zip = os.path.join(era5_path, zip_file)
                break
        
        if not era5_zip:
            print(f"  ✗ Warning: No ERA5 zip file found for site {site_name}")
            failed += 1
            continue
        
        print(f"  ✓ Found ERA5 file: {os.path.basename(era5_zip)}")
        
        try:
            # Read AmeriFlux data
            ameri_file_path = os.path.join(ameri_path, ameri_file)
            ameri_df = pd.read_csv(ameri_file_path)
            
            # Parse datetime for AmeriFlux data
            if 'DateTime' in ameri_df.columns:
                ameri_df['TIMESTAMP'] = pd.to_datetime(ameri_df['DateTime'])
            elif 'datetime' in ameri_df.columns:
                ameri_df['TIMESTAMP'] = pd.to_datetime(ameri_df['datetime'])
            else:
                print(f"  ✗ Error: No datetime column found in AmeriFlux file")
                failed += 1
                continue
            
            # Check AmeriFlux time range
            ameri_min = ameri_df['TIMESTAMP'].min()
            ameri_max = ameri_df['TIMESTAMP'].max()
            print(f"  AmeriFlux time range: {ameri_min} to {ameri_max}")
            print(f"  AmeriFlux rows: {len(ameri_df)}")
            
            # Read ERA5 data (with timezone conversion)
            era5_df = read_era5_from_zip(era5_zip)
            if era5_df is None:
                print(f"  ✗ Error: Could not read ERA5 data")
                failed += 1
                continue
            
            # Check ERA5 time range
            era5_min = era5_df['TIMESTAMP'].min()
            era5_max = era5_df['TIMESTAMP'].max()
            print(f"  ERA5 time range (local): {era5_min} to {era5_max}")
            print(f"  ERA5 rows: {len(era5_df)}")
            
            # Check time overlap
            overlap_start = max(ameri_min, era5_min)
            overlap_end = min(ameri_max, era5_max)
            if overlap_start < overlap_end:
                overlap_days = (overlap_end - overlap_start).days
                print(f"  Time overlap: {overlap_days} days ({overlap_start} to {overlap_end})")
            else:
                print(f"  ⚠ Warning: No time overlap between AmeriFlux and ERA5 data!")
            
            # Ensure TIMESTAMP columns are datetime
            ameri_df['TIMESTAMP'] = pd.to_datetime(ameri_df['TIMESTAMP'])
            era5_df['TIMESTAMP'] = pd.to_datetime(era5_df['TIMESTAMP'])
            
            # Merge dataframes on TIMESTAMP
            merged_df = pd.merge(ameri_df, era5_df, on='TIMESTAMP', how='left')
            
            # Check merge quality
            era5_cols = [col for col in era5_df.columns if col != 'TIMESTAMP']
            if merged_df[era5_cols].notna().any().any():
                merge_success_rate = merged_df[era5_cols[0]].notna().mean() * 100 if era5_cols else 0
                print(f"  ✓ Successfully merged with ERA5 data")
                print(f"  Merge success rate for {era5_cols[0]}: {merge_success_rate:.1f}%")
                print(f"  Merged data shape: {merged_df.shape}")
            else:
                print(f"  ⚠ Warning: No ERA5 data was merged (all NaN)")
            
            # Save merged data
            output_file = f"merged_{site_name}.csv"
            output_file_path = os.path.join(output_path, output_file)
            merged_df.to_csv(output_file_path, index=False)
            print(f"  ✓ Saved merged data to: {output_file_path}")
            
            # Store in results dictionary
            merged_results[site_name] = merged_df
            successful += 1
            
            # Print data completeness summary for key ERA5 variables
            print(f"  Data completeness (% non-missing):")
            key_vars = ['t2m', 'd2m', 'sp', 'ssrd']
            for var in key_vars:
                if var in merged_df.columns:
                    pct = merged_df[var].count() / len(merged_df) * 100
                    print(f"    {var}: {pct:.1f}%")
            
        except Exception as e:
            print(f"  ✗ Error processing site {site_name}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
            continue
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Processing Complete!")
    print(f"Successfully processed: {successful} sites")
    print(f"Failed: {failed} sites")
    print(f"Output saved to: {output_path}")
    print(f"{'='*60}")
    
    return merged_results

def merge_ameri_era5_test(ameri_path, era5_path, output_path, test_sites):
    """
    Test merge for specific sites before processing all
    
    Parameters:
        ameri_path (str): Path to the AmeriFlux data directory
        era5_path (str): Path to the ERA5 point data directory
        output_path (str): Path to save merged data
        test_sites (list): List of site names to test
    
    Returns:
        dict: Dictionary with test site names as keys and merged DataFrames as values
    """
    print("="*60)
    print("TEST MODE: Processing only specified sites")
    print(f"Test sites: {test_sites}")
    print("="*60)
    
    # Create output directory for test results
    test_output_path = os.path.join(output_path, 'test_results')
    os.makedirs(test_output_path, exist_ok=True)
    
    results = {}
    
    for site_name in test_sites:
        # Find AmeriFlux file
        ameri_file = f"gaps_{site_name}.csv"
        ameri_file_path = os.path.join(ameri_path, ameri_file)
        
        if not os.path.exists(ameri_file_path):
            print(f"\n✗ Test site {site_name}: AmeriFlux file not found")
            continue
        
        # Find ERA5 zip file
        era5_zip = None
        for zip_file in os.listdir(era5_path):
            if zip_file.endswith('.zip') and site_name in zip_file:
                era5_zip = os.path.join(era5_path, zip_file)
                break
        
        if not era5_zip:
            print(f"\n✗ Test site {site_name}: ERA5 zip file not found")
            continue
        
        print(f"\n{'='*50}")
        print(f"Testing site: {site_name}")
        print(f"{'='*50}")
        
        try:
            # Read AmeriFlux data
            ameri_df = pd.read_csv(ameri_file_path)
            
            if 'DateTime' in ameri_df.columns:
                ameri_df['TIMESTAMP'] = pd.to_datetime(ameri_df['DateTime'])
            else:
                print(f"  ✗ No DateTime column found")
                continue
            
            # Read ERA5 data
            era5_df = read_era5_from_zip(era5_zip)
            if era5_df is None:
                continue
            
            # Merge
            merged_df = pd.merge(ameri_df, era5_df, on='TIMESTAMP', how='left')
            
            # Save test result
            output_file = f"test_merged_{site_name}.csv"
            output_file_path = os.path.join(test_output_path, output_file)
            merged_df.to_csv(output_file_path, index=False)
            
            results[site_name] = merged_df
            print(f"  ✓ Test passed! Saved to: {output_file_path}")
            
        except Exception as e:
            print(f"  ✗ Test failed: {e}")
            continue
    
    print(f"\n{'='*60}")
    print(f"Test Complete! Results saved to: {test_output_path}")
    print(f"{'='*60}")
    
    return results

# Main execution
if __name__ == "__main__":
    # Define paths
    ameri_path = r'\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\ameri_data\reddy_gaps'
    era5_path = r'M:\Research\WUE_CUE\ameri_data\era5_point_data'
    output_path = r'M:\Research\WUE_CUE\ameri_data\reddy_gaps\merged_ameri_era5'
    
    # Optional: Run test mode first (comment out if you want to skip testing)
    print("="*60)
    print("OPTIONAL TEST MODE")
    print("="*60)
    print("Do you want to run a test on 3 representative sites first?")
    print("Test sites: US-HB2 (Eastern), US-A03 (Alaska), US-Tw1 (Pacific)")
    response = input("Run test mode? (yes/no): ").lower()
    
    if response in ['yes', 'y']:
        test_sites = ['US-HB2', 'US-A03', 'US-Tw1']
        test_results = merge_ameri_era5_test(ameri_path, era5_path, output_path, test_sites)
        
        print("\n" + "="*60)
        print("TEST RESULTS SUMMARY")
        print("="*60)
        print(f"Successfully tested: {len(test_results)}/{len(test_sites)} sites")
        
        if len(test_results) < len(test_sites):
            print("⚠ Some tests failed. Check the output above for details.")
            response = input("\nContinue with processing all sites anyway? (yes/no): ").lower()
            if response not in ['yes', 'y']:
                print("Exiting. Fix issues before processing all sites.")
                exit()
        else:
            print("✓ All tests passed!")
            print("\nProceeding to process all sites...")
    
    # Process all sites
    print("\n" + "="*60)
    print("PROCESSING ALL SITES")
    print("="*60)
    
    merged_data_all = merge_ameri_era5(ameri_path, era5_path, output_path)
    
    print("\nAll processing completed!")
    
    
    
    