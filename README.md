# Water_WUE

 The purpose of my project is to understand, how much WUE, CUE changes within ecosystems. The specific objective of my research includes, 1) whether differences in WUE, CUE across ecosystems resemble ecosystem variability and 2) determine if WUE, CUE in coastal ecosystems denotes salinity stress. As part of this project, I downloaded data from Ameri flux website.


 There are three steps to execute the analysis for this project:

 1) Preprocess Ameriflux raw data from coastal sites in the United States and prepare a data frame to be used for the reddyproc package for further processing
 2) Gap-fill pre proceeded ameriflux dataset using reddyproc package
 3) Calculate GPP and Reco using the light response curve method
 4) Calculate carbon use efficiency and water use efficiency 

Please see details 
    

**1-ameri_api**  

This script uses the `amerifluxr` R package to programmatically download BASE-BADM metadata files for selected AmeriFlux sites.

- Uses the `amf_download_base()` function to retrieve metadata.
- Focuses on sites related to water use efficiency (WUE) and salinity impact studies.
- Saves downloaded metadata locally to a shared network directory.
- Ensures all downloads comply with CCBY4.0 licensing via user agreement.

Function Used: amf_download_base()
Output: BADM .csv files saved to ameri_data directory' 

**2-ameri_un_zip.py** This script defines the function unzip_ameriflux_data(zip_file_path, extracted_folder) which extracts high-frequency (HH or BASE_HH) data files from a downloaded AmeriFlux ZIP archive.
The function:

- Automatically creates the target output folder if it doesn't exist.
- Scans the ZIP archive and extracts only files containing `'HH'` or `'BASE_HH'` in their names.
- Supports extraction to network drives.
- Prints progress messages for debugging.
- Prepares high-frequency flux data for analysis by extracting only relevant files from bulk AmeriFlux downloads.

**3-ameri_preprocess.py** Purpose of this script is
To process half-hourly AmeriFlux .csv files by:

- Processes half-hourly AmeriFlux `.csv` files.
- Handles variable headers and delimiters.
- Parses timestamps into usable datetime indices.
- Reindexes to a continuous 30-minute timestep.
- Fills gaps with NaN and replaces quality control flags (-9999, -6999).
- Creates derived meteorological and flux variables (e.g., NEE, VPD).
- Selects and saves a standardized set of columns to a cleaned CSV.
- Optionally visualizes time series of key variables.


**4-era5_api.py**
This script defines the function fetch_cds_data(area, year_range, month_range, day_range, time_range, output_file) to programmatically download hourly ERA5 single-level reanalysis data (e.g., surface pressure, solar radiation) from the Copernicus Climate Data Store (CDS) using the cdsapi Python client.

**5-merge_ameri_era5.py**
Purpose:

- Merges AmeriFlux half-hourly CSV data with ERA5 NetCDF data.
- Reads AmeriFlux data and converts timestamps.
- Loads and combines ERA5 hourly data from accumulated and instantaneous streams.
- Converts ERA5 timestamps from UTC to the site's local time zone.
- Interpolates ERA5 data to match the AmeriFlux 30-minute resolution.
- Merges both datasets on their aligned timestamps.

**6-blending_ameri_era.py**
Key Steps in blending_ameri_era() and blended_save()

- Converts ERA5 variables to standard units (e.g., temperature, radiation, VPD).
- Fills missing observed Tair and Rg values using ERA5.
- Applies linear regression to generate corrected ERA5 estimates.
- Fills missing data based on seasonal completeness thresholds.
- Recalculates VPD from available Tair and RH if needed.
- Fills missing VPD using corrected or raw ERA5 values.
- Removes unrealistic values for VPD and Rg.
- Generates regression plots for original vs. ERA5 variables.
- Saves the blended output using a standardized filename format.



**7-long_gaps.py**
- Creates the output directory if it doesn't exist.
- Loops through all `.csv` files in the input directory.
- Reads each CSV and ensures the `Month` column is present.
- Computes daily-hourly mean values of LE and NEE across the dataset.
- For each year:
  - Checks if the May–August growing season has ≥50% valid data.
  - If so, identifies April and September gaps (≥7 days long).
  - Fills those long gaps using mean values by DoY and Hour.
- Rounds all filled LE and NEE values to 3 decimal places.
- Saves the processed DataFrame to the specified output folder.
- Prints success or error messages for each file.


**8-Loop_gap_fill_gpp.R**

- Loads AmeriFlux-ERA5 blended CSV files.
- Converts date and time columns to POSIX format.
- Initializes REddyProc with available site variables.
- Estimates uStar threshold:
  - Uses default method.
  - If it fails, applies custom control parameters.
  - If still NA, assigns fallback uStar = 0.1.
- Identifies and removes problematic years with invalid uStar.
- Reinitializes REddyProc after year removal.
- Gap-fills:
  - NEE and LE using uStar filtering.
  - Rg, Tair, and VPD using MDS (without uStar).
- Sets site latitude, longitude, and timezone using metadata.
- Performs flux partitioning:
  - Nighttime: Reichstein method.
  - Daytime: Lasslop method.
- Exports filled variables: NEE_f, LE_f, Tair_f, VPD_f, Rg_f, GPP_DT, Reco_DT, GPP_nt, Reco_nt.
- Saves output to `*_fill.csv` per site.



**9-data_merging.py**

Step-by-Step Workflow
- Converts units and computes ET, GPP, Reco, and NPP.  
  - Fluxes are converted from μmol/m²/s to gC/m² and mm H₂O.  
  - Net primary productivity (NPP) is calculated.

- Merges site metadata.  
  - Site-level attributes such as salinity, location, and biome are joined from a reference metadata file.

- Integrates growing season phenology.  
  - Start (sos) and end (eos) of growing seasons per site-year are merged using PhenoFit-derived values.  
  - Falls back on nearby or average values if missing.

- Generates growing-season monthly summaries.  
  - For each site-year-month within the growing season, computes monthly totals for GPP, ET, NEE, Reco, and derived metrics:  
    - WUE = GPP / ET  
    - CUE = NEP / GPP

- Generates growing-season yearly summaries.  
  - Performs yearly aggregation of carbon and water fluxes.  
  - Filters for quality and stores with metadata and climate averages.

- Visualizes results.  
  - Final plots of WUE, GPP, and ET by site are generated using boxplots to support comparison and interpretation.


**10-growing_season_phenofit.R**

- Converts units and computes ET, GPP, Reco, and NPP.
  - Fluxes are converted from μmol/m²/s to gC/m² and mm H₂O.
  - Net primary productivity (NPP) is calculated.

- Merges site metadata.
  - Site-level attributes such as salinity, location, and biome are joined from a reference metadata file.

- Integrates growing season phenology.
  - Start (sos) and end (eos) of growing seasons per site-year are merged using PhenoFit-derived values.
  - Falls back on nearby or average values if missing.

- Generates growing-season monthly summaries.
  - For each site-year-month within the growing season, computes monthly totals for GPP, ET, NEE, Reco, and derived metrics:
    - WUE = GPP / ET
    - CUE = NEP / GPP

- Generates growing-season yearly summaries.
  - Performs yearly aggregation of carbon and water fluxes.
  - Filters for quality and stores with metadata and climate averages.

- Visualizes results.
  - Final plots of WUE, GPP, and ET by site are generated using boxplots to support comparison and interpretation.

**11-continous_salinity_analysis.py**

This workflow merges continuous salinity records with AmeriFlux carbon and water flux data, focusing on growing-season periods. It performs statistical analysis and generates publication-ready visualizations to explore salinity impacts on ecosystem function.

Main Steps
Merge and preprocess data:

Combine continuous salinity data from multiple sites.

Join site-level salinity with growing season dates (SOS/EOS) from PhenoFit.

Merge growing-season monthly salinity with AmeriFlux flux data (ET, GPP, NEE, Reco, etc.).

Data normalization and cleanup:

Normalize ET and GPP within each site to a [-1, 1] scale to account for magnitude differences.

Handle missing or extreme values (e.g., CUE clipping, negative salinity filtering).

Statistical analysis and plotting:

Perform site-level and combined linear regressions of:

ET vs. salinity

GPP vs. salinity

WUE vs. salinity

CUE vs. salinity

Calculate and display slope, $R^2$, and p-values for each relationship.

Export regression plots for individual sites and combined datasets.

GPP time series visualization:

Plot monthly GPP trends by site, identifying gaps in temporal coverage.





