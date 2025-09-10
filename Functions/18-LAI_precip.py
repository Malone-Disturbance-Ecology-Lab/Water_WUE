# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 16:08:20 2025

@author: ammar
"""

import pandas as pd
import numpy as np
import os
import glob
from datetime import datetime
import warnings
import matplotlib.pyplot as plt
import seaborn as sns

def lai_processing():
    # Define paths
    lai_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai"
    phenofit_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\phenofit_growing_season.csv"
    output_folder = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Read phenofit growing season data
    try:
        phenofit_df = pd.read_csv(phenofit_file)
    except FileNotFoundError:
        print(f"Error: Phenofit file not found at {phenofit_file}")
        return
    
    # Get all CSV files in LAI folder
    csv_files = glob.glob(os.path.join(lai_folder, "*.csv"))
    
    processed_sites = 0
    all_sites_data = []  # Store data for plotting
    
    for file_path in csv_files:
        # Extract site name from filename (first two parts before second hyphen)
        filename = os.path.basename(file_path)
        # Split by hyphen and take first two parts for site code (e.g., "US-A03")
        parts = filename.split('-')
        if len(parts) >= 2:
            site_name = f"{parts[0]}-{parts[1]}"
        else:
            print(f"Skipping {filename}: Cannot extract site name")
            continue
        
        print(f"Processing site: {site_name}")
        
        # Read LAI data
        try:
            lai_df = pd.read_csv(file_path)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue
        
        # Check for both possible LAI column names
        lai_column = None
        for possible_col in ['MOD15A2H_061_Lai_500m', 'MCD15A2H_061_Lai_500m']:
            if possible_col in lai_df.columns:
                lai_column = possible_col
                break
        
        if lai_column is None:
            print(f"Skipping {filename}: No LAI column found. Available columns: {list(lai_df.columns)}")
            continue
        
        # Keep only required columns and rename
        required_cols = ['Latitude', 'Longitude', 'Date', lai_column]
        if not all(col in lai_df.columns for col in required_cols):
            print(f"Skipping {filename}: Missing required columns. Available: {list(lai_df.columns)}")
            continue
        
        lai_df = lai_df[required_cols].copy()
        lai_df.rename(columns={lai_column: 'lai'}, inplace=True)
        
        # Handle outliers - range: 0-10
        lai_df['lai'] = lai_df['lai'].apply(lambda x: np.nan if (x < 0 or x > 10) else x)
        
        # Parse Date to datetime - handle multiple date formats
        try:
            # First try parsing with the format you specified (M/D/YYYY)
            lai_df['Date'] = pd.to_datetime(lai_df['Date'], format='%m/%d/%Y')
        except ValueError:
            try:
                # If that fails, try parsing with ISO format (YYYY-MM-DD)
                lai_df['Date'] = pd.to_datetime(lai_df['Date'], format='%Y-%m-%d')
            except Exception as e:
                print(f"Skipping {filename}: Error parsing dates - {e}")
                continue
        
        # Add DOY column
        lai_df['DOY'] = lai_df['Date'].dt.dayofyear
        
        # Interpolate LAI over time - set Date as index first
        lai_df = lai_df.sort_values('Date')
        lai_df_indexed = lai_df.set_index('Date')
        lai_df_indexed['lai'] = lai_df_indexed['lai'].interpolate(method='time')
        
        # Reset index to get Date back as a column
        lai_df = lai_df_indexed.reset_index()
        
        # Check if site exists in phenofit data
        site_phenofit = phenofit_df[phenofit_df['site_name'] == site_name]
        if site_phenofit.empty:
            warnings.warn(f"No phenofit data found for site {site_name}. Skipping.")
            continue
        
        # Calculate site-wide averages
        avg_sos_site = round(site_phenofit['avg_sos'].mean())
        avg_eos_site = round(site_phenofit['avg_eos'].mean())
        
        # Filter to growing season
        # Find closest DOY at/after SOS and at/before EOS
        all_doy = sorted(lai_df['DOY'].unique())
        
        # Find closest DOY at/after SOS
        sos_doy = min(all_doy, key=lambda x: (x < avg_sos_site, abs(x - avg_sos_site)))
        
        # Find closest DOY at/before EOS
        eos_doy = min(all_doy, key=lambda x: (x > avg_eos_site, abs(x - avg_eos_site)))
        
        # Filter rows within growing season
        growing_season_mask = (lai_df['DOY'] >= sos_doy) & (lai_df['DOY'] <= eos_doy)
        filtered_df = lai_df[growing_season_mask].copy()
        
        if filtered_df.empty:
            warnings.warn(f"No data within growing season for site {site_name}. SOS: {sos_doy}, EOS: {eos_doy}")
            continue
        
        # Add constant columns
        filtered_df['site'] = site_name
        filtered_df['avg_sos'] = avg_sos_site
        filtered_df['avg_eos'] = avg_eos_site
        
        # Reorder columns
        final_cols = ['site', 'Latitude', 'Longitude', 'Date', 'DOY', 'lai', 'avg_sos', 'avg_eos']
        filtered_df = filtered_df[final_cols]
        
        # Save to CSV
        output_file = os.path.join(output_folder, f"{site_name}.csv")
        filtered_df.to_csv(output_file, index=False)
        
        # Store data for plotting - calculate statistics
        site_mean_lai = filtered_df['lai'].mean()
        site_min_lai = filtered_df['lai'].min()
        site_max_lai = filtered_df['lai'].max()
        site_std_lai = filtered_df['lai'].std()
        
        all_sites_data.append({
            'site': site_name,
            'mean_lai': site_mean_lai,
            'min_lai': site_min_lai,
            'max_lai': site_max_lai,
            'std_lai': site_std_lai,
            'latitude': filtered_df['Latitude'].iloc[0],
            'longitude': filtered_df['Longitude'].iloc[0],
            'num_observations': len(filtered_df),
            'sos': avg_sos_site,
            'eos': avg_eos_site
        })
        
        processed_sites += 1
        print(f"Processed {site_name}: {len(filtered_df)} rows, LAI range: {site_min_lai:.2f}-{site_max_lai:.2f}, mean: {site_mean_lai:.3f}")
    
    # Create summary DataFrame for plotting
    if all_sites_data:
        summary_df = pd.DataFrame(all_sites_data)
        
        # Plot 1: Range of mean LAI by site with error bars
        plt.figure(figsize=(16, 10))
        
        # Sort by mean LAI for better visualization
        summary_df = summary_df.sort_values('mean_lai', ascending=False)
        
        # Create the plot with error bars showing range
        plt.errorbar(summary_df['mean_lai'], range(len(summary_df)), 
                    xerr=[summary_df['mean_lai'] - summary_df['min_lai'], 
                         summary_df['max_lai'] - summary_df['mean_lai']],
                    fmt='o', color='darkgreen', alpha=0.7, 
                    capsize=5, capthick=2, elinewidth=2,
                    label='LAI Range (min-max)')
        
        # Add mean points
        plt.scatter(summary_df['mean_lai'], range(len(summary_df)), 
                   s=100, color='lightgreen', edgecolor='darkgreen', 
                   linewidth=2, zorder=5, label='Mean LAI')
        
        # Customize the plot
        plt.yticks(range(len(summary_df)), summary_df['site'])
        plt.xlabel('Leaf Area Index (LAI)', fontsize=12)
        plt.ylabel('Site', fontsize=12)
        plt.title('Range of Mean LAI by Site (Growing Season)', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend()
        
        # Add value annotations
        for i, (mean_val, min_val, max_val) in enumerate(zip(summary_df['mean_lai'], 
                                                           summary_df['min_lai'], 
                                                           summary_df['max_lai'])):
            plt.text(mean_val + 0.1, i, f'{mean_val:.2f}', va='center', fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        plot1_path = os.path.join(output_folder, "lai_range_by_site_errorbars.png")
        plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Plot 2: Box plot showing distribution of LAI ranges
        plt.figure(figsize=(14, 8))
        
        # Prepare data for box plot
        plot_data = []
        site_labels = []
        for _, row in summary_df.iterrows():
            plot_data.append([row['min_lai'], row['mean_lai'], row['max_lai']])
            site_labels.append(row['site'])
        
        # Create box plot
        box = plt.boxplot(plot_data, vert=False, labels=site_labels, patch_artist=True)
        
        # Customize box plot colors
        for patch in box['boxes']:
            patch.set_facecolor('lightgreen')
            patch.set_alpha(0.7)
        
        for whisker in box['whiskers']:
            whisker.set(color='darkgreen', linewidth=2)
        
        for cap in box['caps']:
            cap.set(color='darkgreen', linewidth=2)
        
        for median in box['medians']:
            median.set(color='red', linewidth=2)
        
        plt.xlabel('Leaf Area Index (LAI)', fontsize=12)
        plt.ylabel('Site', fontsize=12)
        plt.title('Distribution of LAI Ranges by Site', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save the plot
        plot2_path = os.path.join(output_folder, "lai_distribution_by_site_boxplot.png")
        plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Plot 3: Bar plot of mean LAI values
        plt.figure(figsize=(14, 8))
        
        # Create bar plot
        bars = plt.bar(range(len(summary_df)), summary_df['mean_lai'], 
                      color='lightgreen', edgecolor='darkgreen', alpha=0.7)
        
        # Add error bars showing standard deviation
        plt.errorbar(range(len(summary_df)), summary_df['mean_lai'], 
                    yerr=summary_df['std_lai'], fmt='none', 
                    color='darkgreen', capsize=5, capthick=2, elinewidth=2)
        
        # Customize the plot
        plt.xticks(range(len(summary_df)), summary_df['site'], rotation=45, ha='right')
        plt.xlabel('Site', fontsize=12)
        plt.ylabel('Mean Leaf Area Index (LAI)', fontsize=12)
        plt.title('Mean LAI by Site with Standard Deviation', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for i, (mean_val, std_val) in enumerate(zip(summary_df['mean_lai'], summary_df['std_lai'])):
            plt.text(i, mean_val + 0.05, f'{mean_val:.2f} ± {std_val:.2f}', 
                    ha='center', va='bottom', fontsize=8, rotation=45)
        
        plt.tight_layout()
        
        # Save the plot
        plot3_path = os.path.join(output_folder, "mean_lai_by_site_bars.png")
        plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        # Print summary statistics
        print(f"\n{'='*60}")
        print("SUMMARY STATISTICS:")
        print(f"{'='*60}")
        print(f"Total sites processed: {len(summary_df)}")
        print(f"Overall mean LAI: {summary_df['mean_lai'].mean():.3f}")
        print(f"Overall LAI range: {summary_df['min_lai'].min():.2f} - {summary_df['max_lai'].max():.2f}")
        print(f"Site with highest mean LAI: {summary_df.loc[summary_df['mean_lai'].idxmax(), 'site']} ({summary_df['mean_lai'].max():.3f})")
        print(f"Site with lowest mean LAI: {summary_df.loc[summary_df['mean_lai'].idxmin(), 'site']} ({summary_df['mean_lai'].min():.3f})")
        print(f"\nPlots saved to:")
        print(f"1. {plot1_path}")
        print(f"2. {plot2_path}")
        print(f"3. {plot3_path}")
        print(f"{'='*60}")
    
    # Print success message
    if processed_sites == 27:
        print("Success: 27 site CSVs were saved!")
    else:
        print(f"Processed {processed_sites} sites (expected 27)")

# Run the function
if __name__ == "__main__":
    lai_processing()


#################################################################################################################






























write a python function

\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai

this folder has mutiple csv files. for example name of one file US-A03-MOD15A2H-061-results
here name name of site is this US-A03

for each site 

I am only interested in these column

Latitude	Longitude	Date	MOD15A2H_061_Lai_500m

remove outliers, any value greater than 1 or less than 0 in column MOD15A2H_061_Lai_500m is outliers
replace outlier with nans

fill nans by interpolating

calculate day of year (DOY) using Date column
this is how date look like 

7/4/2021
7/12/2021
7/20/2021
7/28/2021

Then use DOY to calulate growing season
to caluclaute growing season, you need data from another folder.
here is the folder

\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products
here is name of csv file in this folder that has info for growing season 
phenofit_growing_season


site_name	year	avg_length	avg_sos	avg_eos

match site name for lai data with growing seaosn data using site_name column in
phenofit_growing_season

for example name of lai csv is US-A03-MOD15A2H-061-results and extract site name from name of file:
    
for example US-A03 is site name and match it with data in site_name column:
    in phenofit_growing_season. that column has mutiple site name that match with lai
    mutiple csv files
    
    when desired data is extracted use these columns avg_sos	avg_eos
    to calculate growign season These two columns has Doy info. and growing seaosn
    start with avg_sos and end with avg_eos. if DoY does not exactly matches th sos and eos then use closest Doy
    
    if new DOY column in lai data does not match with growing seaosn data then use the closest DOY.:
        this is because lai data is weekly. now growing seaon should be same for all years for a particular
        site there for take average of growing season for all years and use same growing season for all years
        for a particualr site:
            Howveer growing seaosn should be different across sites
                       
            
 I am intesrted in only one columsn US-A03-MOD15A2H-061-results
 that should be renamed beforeto lai 
 
 save scv files with same  ane as site name such as US-A03   and add date columns lai and avg_sos avg_eos columns

 refrence data in this column is this folder is half hourly so convert lai to half hourly, using sami lai for a given 
                        
            
            
merge folder
\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\canopy_conduc\gs-calculated\2-evaporation_transpiration_split\3-growing_season_transp_Evap

read csv files in this folder and merge columns lai column . I am intesrted in only one columsn US-A03-MOD15A2H-061-results
that should be renamed beforeto lai before merging

refrence data in this column is this folder is half hourly so convert lai to half hourly, using sami lai for a given 
save folder for csv file

 \\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\drivers\ameri_drivers\lai\growing_season_LAI           
            
 print message if all 27 files are saved
           