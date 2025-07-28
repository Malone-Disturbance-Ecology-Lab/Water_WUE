## 🌿 Water_WUE: Ecosystem Efficiency and Salinity Stress Analysis

This project investigates spatial and temporal variations in Water Use Efficiency (WUE) and Carbon Use Efficiency (CUE) across U.S. AmeriFlux sites, with a focus on how salinity stress affects coastal ecosystems. The workflow integrates AmeriFlux and MODIS data, applies gap-filling and partitioning methods, and performs trend and salinity-based statistical analyses.

### 🔍 Project Objectives

- Evaluate how WUE and CUE vary across ecosystems and over time  
- Assess whether coastal ecosystem WUE and CUE indicate salinity stress  
- Integrate remote sensing and flux tower data to analyze ecosystem function  

### ⚙️ Code Modules Overview

- `1-ameri_api.R`  
  - Downloads BASE and BADM metadata from AmeriFlux using `amerifluxr`

- `2-ameri_un_zip.py`  
  - Unzips and filters high-frequency AmeriFlux (HH) data from raw archives

- `3-ameri_preprocess.py`  
  - Processes half-hourly AmeriFlux data: cleans, standardizes, and prepares for gap-filling

- `4-era5_api.py`  
  - Downloads ERA5 reanalysis data (e.g., radiation, temperature) using the Copernicus CDS API

- `5-merge_ameri_era5.py`  
  - Merges AmeriFlux and ERA5 datasets, aligning them to 30-min resolution

- `6-blending_ameri_era.py`  
  - Fills missing meteorological data using ERA5 corrections and linear regression

- `7-long_gaps.py`  
  - Fills long seasonal gaps in LE and NEE using climatological hourly means

- `8-Loop_gap_fill_gpp.R`  
  - Gap-fills flux data using REddyProc and performs GPP/Reco partitioning

- `9-data_merging.py`  
  - Aggregates AmeriFlux data by season/year, computes WUE, CUE, and merges with site metadata

- `10-growing_season_phenofit.R`  
  - Integrates PhenoFit-based growing season metrics and summarizes seasonal fluxes

- `11-continous_salinity_analysis.py`  
  - Combines continuous salinity records with flux data; analyzes and visualizes salinity–WUE/GPP/CUE relationships

- `12-canopy_conduc_gs.R`  
  - Estimates canopy conductance during growing season using FG and iPM methods

- `13-modis_grid.R`  
  - Generates 3×3 MODIS sinusoidal grids (500 m resolution) around flux tower locations

- `14-modis_analysis.py`  
  - Processes MODIS ET and GPP data, filters by phenological growing season, and compares MODIS vs. AmeriFlux trends

### 📁 Outputs

- Cleaned and gap-filled flux data (`*_fill.csv`)  
- Growing season summaries and trends (WUE, CUE, ET, GPP)  
- Long-term trend plots and salinity regressions (`.png`)  
- MODIS grid shapefiles and comparative MODIS–AmeriFlux scatterplots
