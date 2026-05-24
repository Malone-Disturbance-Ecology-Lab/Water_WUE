# -*- coding: utf-8 -*-
"""
Focus only on identifying sites with invalid LAI data
No file saving, no plots, just reporting
"""

import pandas as pd
import numpy as np
import os
import glob

def check_invalid_lai_sites():
    """
    Only identifies and prints sites with invalid LAI data
    Does NOT save any files or create plots
    """
    # Define paths
    lai_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai"
    phenology_folder = r"M:\Research\WUE_CUE\ameri_data\ameri_phenofit"
    
    print("="*70)
    print("CHECKING FOR SITES WITH INVALID LAI DATA")
    print("="*70)
    
    # Read all phenology files
    phenology_files = glob.glob(os.path.join(phenology_folder, "*_phenofit*.csv"))
    print(f"\nFound {len(phenology_files)} phenology files")
    
    # Create a dictionary to store phenology data for each site
    phenology_dict = {}
    for pheno_file in phenology_files:
        filename = os.path.basename(pheno_file)
        site_name = filename.split('_')[0]
        
        try:
            pheno_df = pd.read_csv(pheno_file)
            if 'avg_sos' in pheno_df.columns and 'avg_eos' in pheno_df.columns:
                avg_sos_site = round(pheno_df['avg_sos'].mean())
                avg_eos_site = round(pheno_df['avg_eos'].mean())
                phenology_dict[site_name] = {
                    'avg_sos': avg_sos_site,
                    'avg_eos': avg_eos_site
                }
        except Exception:
            continue
    
    print(f"Loaded phenology data for {len(phenology_dict)} sites")
    
    # Get all CSV files in LAI folder
    csv_files = glob.glob(os.path.join(lai_folder, "*.csv"))
    print(f"Found {len(csv_files)} LAI files to process")
    
    # Track site status
    sites_with_valid_data = set()
    sites_appearing_in_lai = set()
    missing_phenology_sites = set()
    
    # First, identify all sites that appear in LAI data
    for file_path in csv_files:
        try:
            lai_df = pd.read_csv(file_path)
        except Exception:
            continue
        
        # Check for required columns
        if not all(col in lai_df.columns for col in ['Latitude', 'Longitude', 'Date', 'ID']):
            continue
        
        # Check for LAI column
        lai_column = None
        for possible_col in ['MOD15A2H_061_Lai_500m', 'MCD15A2H_061_Lai_500m']:
            if possible_col in lai_df.columns:
                lai_column = possible_col
                break
        
        if lai_column is None:
            continue
        
        # Add all sites from this file
        unique_sites = lai_df['ID'].unique()
        for site_name in unique_sites:
            sites_appearing_in_lai.add(site_name)
    
    # Now check for valid data in sites that appear in LAI
    for file_path in csv_files:
        try:
            lai_df = pd.read_csv(file_path)
        except Exception:
            continue
        
        # Check for required columns
        if not all(col in lai_df.columns for col in ['Latitude', 'Longitude', 'Date', 'ID']):
            continue
        
        # Check for LAI column
        lai_column = None
        for possible_col in ['MOD15A2H_061_Lai_500m', 'MCD15A2H_061_Lai_500m']:
            if possible_col in lai_df.columns:
                lai_column = possible_col
                break
        
        if lai_column is None:
            continue
        
        # Get unique sites in this LAI file
        unique_sites = lai_df['ID'].unique()
        
        for site_name in unique_sites:
            # Skip if already found valid data
            if site_name in sites_with_valid_data:
                continue
            
            # Check if site has phenology data
            if site_name not in phenology_dict:
                missing_phenology_sites.add(site_name)
                continue
            
            # Filter data for this site
            site_data = lai_df[lai_df['ID'] == site_name].copy()
            
            if len(site_data) == 0:
                continue
            
            # Keep only required columns
            site_data = site_data[['Date', lai_column]].copy()
            site_data.rename(columns={lai_column: 'lai'}, inplace=True)
            
            # Remove extreme outliers (outside 0-10)
            site_data['lai'] = site_data['lai'].apply(lambda x: np.nan if (x < 0 or x > 10) else x)
            
            # Parse dates
            date_parsed = False
            for date_format in ['%m/%d/%Y', '%Y-%m-%d', '%m/%d/%y', '%Y/%m/%d']:
                try:
                    site_data['Date'] = pd.to_datetime(site_data['Date'], format=date_format)
                    date_parsed = True
                    break
                except:
                    continue
            
            if not date_parsed:
                try:
                    site_data['Date'] = pd.to_datetime(site_data['Date'])
                    date_parsed = True
                except Exception:
                    continue
            
            # Add DOY
            site_data['DOY'] = site_data['Date'].dt.dayofyear
            
            # Get phenology data
            avg_sos_site = phenology_dict[site_name]['avg_sos']
            avg_eos_site = phenology_dict[site_name]['avg_eos']
            
            # Filter to growing season with 2-day buffer
            sos_buffer = 2
            eos_buffer = 2
            
            growing_season_mask = ((site_data['DOY'] >= (avg_sos_site - sos_buffer)) & 
                                  (site_data['DOY'] <= (avg_eos_site + eos_buffer)))
            filtered_df = site_data[growing_season_mask].copy()
            
            if filtered_df.empty:
                continue
            
            # Check if we have any valid LAI values
            valid_lai_count = filtered_df['lai'].count()
            
            if valid_lai_count > 0:
                sites_with_valid_data.add(site_name)
    
    # Determine problematic sites
    sites_with_phenology = set(phenology_dict.keys())
    
    # Sites that have phenology but never appear in any LAI file
    sites_missing_from_lai = sites_with_phenology - sites_appearing_in_lai
    
    # Sites that appear in LAI but have no valid data after filtering
    sites_with_no_valid_lai = (sites_with_phenology & sites_appearing_in_lai) - sites_with_valid_data
    
    # Sites that appear in LAI but have no phenology
    missing_phenology_sites = sites_appearing_in_lai - sites_with_phenology
    
    # Print report
    print("\n" + "="*70)
    print("INVALID LAI SITES REPORT")
    print("="*70)
    
    if missing_phenology_sites:
        print(f"\n❌ SITES WITH NO PHENOLOGY DATA: {len(missing_phenology_sites)}")
        print("-"*70)
        for site in sorted(missing_phenology_sites):
            print(f"   • {site}")
    
    if sites_missing_from_lai:
        print(f"\n🔍 SITES WITH PHENOLOGY BUT NEVER APPEAR IN LAI DATA: {len(sites_missing_from_lai)}")
        print("-"*70)
        for site in sorted(sites_missing_from_lai):
            print(f"   • {site}")
    
    if sites_with_no_valid_lai:
        print(f"\n⚠️  SITES IN LAI DATA BUT NO VALID LAI VALUES IN GROWING SEASON: {len(sites_with_no_valid_lai)}")
        print("-"*70)
        for site in sorted(sites_with_no_valid_lai):
            print(f"   • {site}")
    
    if not missing_phenology_sites and not sites_missing_from_lai and not sites_with_no_valid_lai:
        print("\n✅ All sites have valid data!")
    
    print("\n" + "="*70)
    print("SUMMARY STATISTICS:")
    print("="*70)
    print(f"Total sites with phenology data: {len(sites_with_phenology)}")
    print(f"Total sites appearing in LAI data: {len(sites_appearing_in_lai)}")
    print(f"Total sites with valid LAI data: {len(sites_with_valid_data)}")
    print(f"Total invalid sites: {len(missing_phenology_sites) + len(sites_missing_from_lai) + len(sites_with_no_valid_lai)}")
    print("="*70 + "\n")

# Run only the check function
if __name__ == "__main__":
    check_invalid_lai_sites()