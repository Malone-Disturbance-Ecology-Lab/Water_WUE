## 🌿 Water_WUE: Coastal Ecosystem Water Use Efficiency and Drought Analysis

This project investigates spatial and temporal variability in ecosystem Water Use Efficiency (WUE), transpiration-based WUE (WUE_T), and evapotranspiration partitioning across U.S. coastal ecosystems using AmeriFlux observations, remote sensing products, and hydroclimatic drought indices. The workflow integrates flux tower, ERA5, MODIS, phenology, salinity, and SPEI datasets to evaluate how prolonged drought, precipitation variability, and salinity gradients influence ecosystem carbon–water coupling.

### 🔍 Project Objectives

- **Q1:** Quantify differences in WUE_ET, WUEₜ, and evapotranspiration partitioning across U.S. coastal ecosystems under near-normal hydroclimatic conditions  

- **Q2:** Evaluate how short- and long-term hydroclimatic anomalies alter WUE_ET and WUEₜ responses across coastal ecosystems using multi-timescale SPEI  

- **Q3:** Assess site-level WUEₜ sensitivity to persistent multi-year moisture anomalies and determine whether sensitivity patterns differ among ecosystem types, climate–biome groups, salinity gradients, and coastal regions  

- **Q4:** Determine how persistent multi-year drought influences the probability, spatial persistence, and regional patterns of WUEₜ decline across U.S. coastal regions  

- Integrate AmeriFlux observations, ET partitioning, ERA5 climate reanalysis, MODIS vegetation products, and SPEI datasets within a unified coastal carbon–water analysis framework  

### ⚙️ Code Modules Overview

- `1-ameri_api.R`  
  Downloads AmeriFlux BASE and BADM metadata

- `2-ameri_un_zip.py`  
  Extracts and organizes raw AmeriFlux half-hourly datasets

- `3-ameri_preprocess.py`  
  Cleans and standardizes AmeriFlux flux and meteorological observations

- `4-Era5_point_data.py`  
  Downloads ERA5 meteorological reanalysis data for flux tower locations

- `5-merge_ameri_era5.py`  
  Merges AmeriFlux and ERA5 datasets

- `6-blending_ameri_era.py`  
  Performs meteorological gap-filling and ERA5 blending corrections

- `7-long_gaps.py`  
  Fills extended gaps in flux observations using climatological approaches

- `8a-Loop_gap_fill_gpp.R`  
  Performs REddyProc-based gap-filling and flux partitioning

- `8b-Loop_gap_fill_gpp_unique_sites.R`  
  Applies gap-filling workflow for unique-site processing

- `9-gap_filling_check.py`  
  Evaluates and visualizes gap-filling performance

- `10-data_merging.py`  
  Aggregates processed datasets and computes WUE metrics

- `11-growing_season_phenofit.R`  
  Extracts phenology-based growing season metrics

- `11b-growing_season_phenofit_unique_sites.R`  
  Growing season analysis for unique-site datasets

- `12-LAI_growing_season.py`  
  Processes LAI dynamics during growing seasons

- `13-ameri_flux_lai_precip_growingS.py`  
  Integrates flux, LAI, and precipitation datasets

- `14-Find_elevation.py`  
  Extracts elevation information for site-level analyses

- `15-elevation_driver_data.R`  
  Processes elevation-based environmental drivers

- `16-ET_partioning_may_2026.py`  
  Performs evapotranspiration partitioning analyses

- `16-precip_driver.R`  
  Generates precipitation-based hydroclimatic drivers

- `17-add_biome_plot_partitioning.py`  
  Adds biome classifications and visualization outputs

- `18-drought_indices_api.py`  
  Downloads and processes drought indices including SPEI

- `18-drought_indices_input_data.R`  
  Prepares drought datasets for analysis workflows

- `19-merge_ameri_WUE_indices.py`  
  Merges WUE metrics with hydroclimatic drought indices

- `20-WUE_cleaning_plot.py`  
  Cleans datasets and generates WUE visualization outputs

- `21-coastal buffer for modis.py`  
  Creates MODIS coastal buffer grids around flux tower sites

- `22-SPEI_spatial_visual.py`  
  Generates spatial drought and SPEI visualization products

### 📁 Outputs

- Gap-filled and partitioned flux datasets (`*_fill.csv`)  
- Seasonal and annual WUE/WUEₜ summaries  
- ET partitioning outputs and drought sensitivity analyses  
- SPEI and hydroclimatic driver datasets  
- Growing season phenology summaries  
- Spatial drought visualizations and coastal MODIS products  
- Publication-ready figures and diagnostic plots  
