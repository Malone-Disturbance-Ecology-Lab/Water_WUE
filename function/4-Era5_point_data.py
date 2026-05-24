# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:23:47 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
ERA5 point downloader for a single site
"""

import cdsapi
import os


def round_to_quarter(value):
    return round(value * 4) / 4


def fetch_site_csv(site_name, lat, lon, start_date, end_date, variables):
    lat_grid = round_to_quarter(lat)
    lon_grid = round_to_quarter(lon)

    client = cdsapi.Client()
    dataset = "reanalysis-era5-single-levels-timeseries"

    output_file = f"{site_name}_{start_date}_{end_date}.csv".replace(":", "-")

    request = {
        "variable": variables,
        "location": {
            "longitude": lon_grid,
            "latitude": lat_grid
        },
        "date": [f"{start_date}/{end_date}"],
        "data_format": "csv"
    }

    client.retrieve(dataset, request).download(target=output_file)

    print(f"{site_name} done")
    print(f"  original: lat={lat}, lon={lon}")
    print(f"  snapped : lat={lat_grid}, lon={lon_grid}")
    print(f"  saved   : {output_file}")


# Set credentials once
os.environ['CDSAPI_URL'] = 'https://cds.climate.copernicus.eu/api'
os.environ['CDSAPI_KEY'] = '0a2b0eec-ed62-4a42-8fb8-83ca726ef382'

# ============================================
# EDIT THESE VALUES FOR YOUR SITE
# ============================================
site_name = "US-xSE"      # Change this to your site name
latitude = 38.8901 # Change this to your latitude
longitude = -76.56  # Change this to your longitude
start_date = "2017-01-01" # Change this to your start date
end_date = "2025-12-31"   # Change this to your end date
# ===========================================



variables = [
    "2m_dewpoint_temperature",
    "surface_pressure",
    "surface_solar_radiation_downwards",
    "2m_temperature", "surface_net_solar_radiation",
        "surface_net_thermal_radiation"
]

# Download the data
fetch_site_csv(
    site_name=site_name,
    lat=latitude,
    lon=longitude,
    start_date=start_date,
    end_date=end_date,
    variables=variables
)

##############################################################################################
## check if ERA5 point data doenloaded correctly 

# -*- coding: utf-8 -*-
"""
Automated ERA5 Data Validator - Automatically detects and reports problem sites
"""

import pandas as pd
import zipfile
import os
from pathlib import Path
from datetime import datetime
import numpy as np


def round_to_quarter(value):
    """Round coordinates to nearest 0.25 degrees (ERA5 grid resolution)"""
    if pd.isna(value):
        return None
    return round(value * 4) / 4


def check_era5_zip_file(zip_path, site_name, expected_lat, expected_lon, expected_start_year, expected_end_year):
    """
    Check if the ERA5 zip file contains correct data for the site
    Returns: (is_valid, issues_list, actual_coords, actual_dates)
    """
    issues = []
    actual_lat = None
    actual_lon = None
    actual_start = None
    actual_end = None
    
    # Check if zip file exists
    if not os.path.exists(zip_path):
        issues.append("MISSING")
        return False, issues, (actual_lat, actual_lon), (actual_start, actual_end)
    
    try:
        # Open zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_file:
            # Find CSV file in zip
            csv_files = [f for f in zip_file.namelist() if f.endswith('.csv')]
            
            if not csv_files:
                issues.append("NO_CSV")
                return False, issues, (actual_lat, actual_lon), (actual_start, actual_end)
            
            # Read CSV data
            with zip_file.open(csv_files[0]) as f:
                df = pd.read_csv(f)
            
            # Check if required columns exist
            required_cols = ['latitude', 'longitude', 'valid_time']
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                issues.append(f"MISSING_COLUMNS: {missing_cols}")
                return False, issues, (actual_lat, actual_lon), (actual_start, actual_end)
            
            # Get actual coordinates
            actual_lat = df['latitude'].iloc[0]
            actual_lon = df['longitude'].iloc[0]
            
            # Check coordinates
            expected_lat_snapped = round_to_quarter(expected_lat)
            expected_lon_snapped = round_to_quarter(expected_lon)
            
            lat_match = abs(actual_lat - expected_lat_snapped) < 0.01
            lon_match = abs(actual_lon - expected_lon_snapped) < 0.01
            
            if not lat_match or not lon_match:
                issues.append(f"COORD_MISMATCH")
            
            # Parse date range
            df['valid_time'] = pd.to_datetime(df['valid_time'])
            actual_start = df['valid_time'].min()
            actual_end = df['valid_time'].max()
            
            # Check date range if specified
            if not pd.isna(expected_start_year) and not pd.isna(expected_end_year):
                if actual_start.year != int(expected_start_year) or actual_end.year != int(expected_end_year):
                    issues.append(f"DATE_MISMATCH")
            
            # Check data completeness
            if not pd.isna(expected_start_year) and not pd.isna(expected_end_year):
                expected_start = datetime(int(expected_start_year), 1, 1)
                expected_end = datetime(int(expected_end_year), 12, 31)
                expected_hours = (expected_end - expected_start).days * 24 + 24
                actual_hours = len(df)
                
                if actual_hours < expected_hours * 0.95:
                    issues.append(f"INCOMPLETE_DATA")
            
            # Check for missing variables
            required_vars = ['d2m', 't2m', 'sp', 'ssrd']
            missing_vars = [var for var in required_vars if var not in df.columns]
            if missing_vars:
                issues.append(f"MISSING_VARS: {missing_vars}")
    
    except Exception as e:
        issues.append(f"ERROR: {str(e)}")
    
    is_valid = len(issues) == 0
    return is_valid, issues, (actual_lat, actual_lon), (actual_start, actual_end)


def validate_all_sites():
    """
    Validate all sites and automatically identify problem sites
    """
    # Paths
    csv_path = r"M:\Research\WUE_CUE\data_products\info\site_lat_long.csv"
    era5_base_path = r"M:\Research\WUE_CUE\ameri_data\era5_data"
    
    # Read CSV
    print("="*80)
    print("AUTOMATED ERA5 DATA VALIDATION")
    print("="*80)
    
    df = pd.read_csv(csv_path)
    
    # Remove rows with NaN in required columns
    df_clean = df.dropna(subset=['site_name', 'lat', 'long'])
    
    print(f"\n📊 Total sites in CSV: {len(df)}")
    print(f"📊 Valid sites to check: {len(df_clean)}")
    print(f"📊 Skipped invalid rows: {len(df) - len(df_clean)}")
    
    # Store results
    valid_sites = []
    problem_sites = []
    
    print("\n" + "-"*80)
    print("CHECKING EACH SITE...")
    print("-"*80)
    
    for idx, row in df_clean.iterrows():
        site_name = row['site_name']
        expected_lat = row['lat']
        expected_lon = row['long']
        
        # Handle date columns
        expected_start = int(row['start_date']) if not pd.isna(row['start_date']) else np.nan
        expected_end = int(row['end_date']) if not pd.isna(row['end_date']) else np.nan
        
        # Look for zip file
        zip_path = os.path.join(era5_base_path, f"{site_name}.zip")
        
        # Validate
        is_valid, issues, actual_coords, actual_dates = check_era5_zip_file(
            zip_path, site_name, expected_lat, expected_lon, expected_start, expected_end
        )
        
        if is_valid:
            valid_sites.append(site_name)
            print(f"   ✓ {site_name}")
        else:
            actual_lat, actual_lon = actual_coords
            actual_start, actual_end = actual_dates
            
            problem_sites.append({
                'site': site_name,
                'issues': issues,
                'expected_lat': expected_lat,
                'expected_lon': expected_lon,
                'expected_snapped_lat': round_to_quarter(expected_lat),
                'expected_snapped_lon': round_to_quarter(expected_lon),
                'actual_lat': actual_lat,
                'actual_lon': actual_lon,
                'expected_start': expected_start,
                'expected_end': expected_end,
                'actual_start': actual_start.year if actual_start else None,
                'actual_end': actual_end.year if actual_end else None
            })
            print(f"   ✗ {site_name} - PROBLEM DETECTED")
    
    # Print detailed problem report
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    print(f"\n✅ VALID SITES: {len(valid_sites)}")
    print(f"❌ PROBLEM SITES: {len(problem_sites)}")
    
    if problem_sites:
        print("\n" + "-"*80)
        print("DETAILED PROBLEM REPORT")
        print("-"*80)
        
        for i, site_info in enumerate(problem_sites, 1):
            print(f"\n{i}. {site_info['site']}")
            print(f"   Issues detected: {', '.join(site_info['issues'])}")
            print(f"   Expected location (snapped): ({site_info['expected_snapped_lat']}, {site_info['expected_snapped_lon']})")
            print(f"   Actual location in ERA5: ({site_info['actual_lat']}, {site_info['actual_lon']})")
            
            # Calculate offset if coordinates are wrong
            if site_info['actual_lat'] and site_info['expected_snapped_lat']:
                lat_offset = abs(site_info['actual_lat'] - site_info['expected_snapped_lat'])
                lon_offset = abs(site_info['actual_lon'] - site_info['expected_snapped_lon'])
                if lat_offset > 0 or lon_offset > 0:
                    print(f"   Offset: {lat_offset:.2f}° lat, {lon_offset:.2f}° lon")
            
            # Show date range issues
            if site_info['expected_start'] and not pd.isna(site_info['expected_start']):
                print(f"   Expected date range: {site_info['expected_start']} to {site_info['expected_end']}")
                print(f"   Actual date range: {site_info['actual_start']} to {site_info['actual_end']}")
            
            print(f"   Action needed: Re-download with correct coordinates")
    
    # Create summary file
    report_path = os.path.join(era5_base_path, f"validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
    with open(report_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("ERA5 DATA VALIDATION REPORT (AUTOMATED)\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Total sites processed: {len(df_clean)}\n")
        f.write(f"Valid sites: {len(valid_sites)}\n")
        f.write(f"Problem sites: {len(problem_sites)}\n\n")
        
        if valid_sites:
            f.write("VALID SITES:\n")
            f.write("-"*40 + "\n")
            for site in valid_sites:
                f.write(f"  {site}\n")
        
        if problem_sites:
            f.write("\nPROBLEM SITES:\n")
            f.write("-"*40 + "\n")
            for site_info in problem_sites:
                f.write(f"\n{site_info['site']}:\n")
                f.write(f"  Issues: {', '.join(site_info['issues'])}\n")
                f.write(f"  Expected snapped: ({site_info['expected_snapped_lat']}, {site_info['expected_snapped_lon']})\n")
                f.write(f"  Actual coordinates: ({site_info['actual_lat']}, {site_info['actual_lon']})\n")
                if site_info['expected_start'] and not pd.isna(site_info['expected_start']):
                    f.write(f"  Expected dates: {site_info['expected_start']}-{site_info['expected_end']}\n")
                    f.write(f"  Actual dates: {site_info['actual_start']}-{site_info['actual_end']}\n")
    
    print(f"\n💾 Full report saved to: {report_path}")
    
    # Return lists for further use if needed
    return valid_sites, problem_sites


if __name__ == "__main__":
    valid, problems = validate_all_sites()
    
    # Final summary
    print("\n" + "="*80)
    print("AUTOMATED VALIDATION COMPLETE")
    print("="*80)
    
    if problems:
        print(f"\n⚠️ {len(problems)} SITE(S) NEED ATTENTION:")
        for site_info in problems:
            print(f"   • {site_info['site']}")
        
        print("\n📋 RECOMMENDATION:")
        print("   Delete the problematic zip files and re-download using the single-site code")
        print("   with the correct coordinates from your CSV file")
    else:
        print("\n✅ ALL SITES HAVE CORRECT DATA!")
    
    print(f"\n📊 Statistics:")
    print(f"   Valid: {len(valid)} sites")
    print(f"   Problematic: {len(problems)} sites")
    print(f"   Success rate: {len(valid)/(len(valid)+len(problems))*100:.1f}%")

####################################################################################################


































