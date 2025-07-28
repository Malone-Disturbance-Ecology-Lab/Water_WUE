# Water_WUE

 The purpose of my project is to understand, how much WUE, CUE changes within ecosystems. The specific objective of my research includes, 1) whether differences in WUE, CUE across ecosystems resemble ecosystem variability and 2) determine if WUE, CUE in coastal ecosystems denotes salinity stress. As part of this project, I downloaded data from Ameri flux website.


 There are three steps to execute the analysis for this project:

 1) Preprocess Ameriflux raw data from coastal sites in the United States and prepare a data frame to be used for the reddyproc package for further processing
 2) Gap-fill pre proceeded ameriflux dataset using reddyproc package
 3) Calculate GPP and Reco using the light response curve method
 4) Calculate carbon use efficiency and water use efficiency 


    

**1-ameri_api**  This script uses the amerifluxr R package to programmatically download BASE-BADM metadata files for selected AmeriFlux sites. The amf_download_base() function retrieves data for sites related to water use efficiency (WUE) and salinity impact studies. Metadata is saved locally to a shared network directory. All downloads comply with CCBY4.0 licensing via user agreement.
Function Used: amf_download_base()
Output: BADM .csv files saved to ameri_data directory' 

**2-ameri_un_zip.py** This script defines the function unzip_ameriflux_data(zip_file_path, extracted_folder) which extracts high-frequency (HH or BASE_HH) data files from a downloaded AmeriFlux ZIP archive.
The function:
Automatically creates the target output folder if it doesn't exist.
Scans the ZIP archive and extracts only files containing 'HH' or 'BASE_HH' in their names.
Supports extraction to network drives and prints progress for debugging.
Use case: Prepares high-frequency flux data for analysis by extracting only relevant files from bulk AmeriFlux downloads.. 

**3-ameri_preprocess.py** Purpose of this script is
To process half-hourly AmeriFlux .csv files by:
Handling variable headers and delimiters
Parsing timestamps into usable datetime indices
Reindexing to a continuous 30-minute timestep
Filling gaps with NaN and replacing quality control flags (-9999, -6999)
Creating derived meteorological and flux variables (e.g., NEE, VPD)
Selecting and saving a standardized set of columns to a cleaned CSV
Optionally visualizing time series of key variables

**4-era5_api.py**
This script defines the function fetch_cds_data(area, year_range, month_range, day_range, time_range, output_file) to programmatically download hourly ERA5 single-level reanalysis data (e.g., surface pressure, solar radiation) from the Copernicus Climate Data Store (CDS) using the cdsapi Python client.

**5-merge_ameri_era5.py**
Purpose:
To merge AmeriFlux half-hourly CSV data with ERA5 NetCDF data by:
Reading AmeriFlux data and converting timestamps
Loading and combining ERA5 hourly data from accumulated and instantaneous streams
Converting ERA5 timestamps from UTC to the site's local time zone
Interpolating ERA5 to match the AmeriFlux 30-minute resolution
Merging both datasets on their aligned timestamps

**6-blending_ameri_era**
Key Steps in blending_ameri_era() and blended_save()
Convert ERA5 variables to standard units (e.g., temperature, radiation, VPD)
Fill missing observed Tair and Rg values using ERA5
Apply linear regression to generate corrected ERA5 estimates
Fill missing data based on seasonal completeness thresholds
Recalculate VPD from available Tair and RH if needed
Fill missing VPD using corrected or raw ERA5 values
Remove unrealistic values for VPD and Rg
Generate regression plots for original vs. ERA5 variables
Save the blended output using a standardized filename format
