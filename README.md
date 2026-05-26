# 🌿 Water_WUE: Coastal Ecosystem Water Use Efficiency and Drought Analysis

This project investigates spatial and temporal variability in ecosystem Water Use Efficiency (WUE), transpiration-based WUE (WUEₜ), and evapotranspiration partitioning across U.S. coastal ecosystems using AmeriFlux observations, remote sensing products, and hydroclimatic drought indices. The workflow integrates flux tower observations, ERA5 climate reanalysis, MODIS products, ET partitioning, phenology, salinity, and SPEI drought datasets to evaluate how prolonged hydroclimatic stress alters ecosystem carbon–water coupling across coastal regions.

---

# 🔍 Project Objectives

- **Q1:** Quantify differences in WUE_ET, WUEₜ, and evapotranspiration partitioning across U.S. coastal ecosystems under near-normal hydroclimatic conditions  

- **Q2:** Evaluate how short- and long-term hydroclimatic anomalies alter WUE_ET and WUEₜ responses across coastal ecosystems using multi-timescale SPEI  

- **Q3:** Assess site-level WUEₜ sensitivity to persistent multi-year moisture anomalies and determine whether sensitivity patterns differ among ecosystem types, climate–biome groups, salinity gradients, and coastal regions  

- **Q4:** Determine how persistent multi-year drought influences the probability, spatial persistence, and regional patterns of WUEₜ decline across U.S. coastal regions  

- Integrate AmeriFlux observations, ET partitioning, ERA5 climate reanalysis, MODIS vegetation products, and SPEI datasets within a unified coastal carbon–water analysis framework  

---

# 📂 Data Directory Structure

Raw, intermediate, and processed datasets are organized outside the repository on the Malone Lab server:

```text
M:\Research\WUE_CUE\
```

---

## Raw AmeriFlux input data

```text
M:\Research\WUE_CUE\ameri_data\
```

Contains raw AmeriFlux `.zip` archives downloaded from AmeriFlux.

Processing begins with:

```text
2-ameri_un_zip.py
```

This script extracts `HH` / `BASE_HH` half-hourly AmeriFlux files from raw archives.

Output directory:

```text
M:\Research\WUE_CUE\ameri_data\ameri_gaps\
```

Extracted half-hourly AmeriFlux files used for preprocessing, gap-filling, and ET partitioning workflows.

---

## ERA5 meteorological forcing data

Generated using:

```text
4-Era5_point_data.py
```

Output directory:

```text
M:\Research\WUE_CUE\era5_point_data\
```

Contains site-level ERA5 meteorological forcing datasets for each AmeriFlux site.

---

## Final merged WUE datasets

Final integration and metric calculations are performed using:

```text
19-merge_ameri_WUE_indices.py
```

Inputs:
- ET partitioning outputs
- AmeriFlux processed datasets
- ERA5 meteorological drivers
- site metadata
- SPEI drought indices

Main input directories:

```text
M:\Research\WUE_CUE\ameri_data\ET_partitioning\
```

```text
M:\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\drought\
```

```text
M:\Research\WUE_CUE\data_products\info\
```

Final merged outputs:

```text
M:\Research\WUE_CUE\data_products\
```

Key output files:

```text
WUE_CUE_monthly_merged_indices.csv
WUE_CUE_yearly_merged_indices.csv
```

---

## Final quality-control and outlier filtering

Final cleaning and outlier filtering are performed using:

```text
20-WUE_cleaning_plot.py
```

Input files:

```text
WUE_CUE_monthly_merged_indices.csv
WUE_CUE_yearly_merged_indices.csv
```

Final cleaned outputs:

```text
WUE_CUE_monthly_merged_indices_clean.csv
WUE_CUE_yearly_merged_indices_clean.csv
```

Additional outputs include:
- outlier audit summaries
- IQR filtering summaries
- diagnostic plots
- quality-control logs

These cleaned datasets are the primary inputs used for downstream manuscript analyses and figure-generation workflows.

---

# 📊 Manuscript Analysis and Figure Workflows

The scripts above generate the cleaned WUE, WUEₜ, ET partitioning, drought, and hydroclimatic datasets used for downstream manuscript analyses.

Reviewers interested only in the final analyses for manuscript figures and statistical workflows do **not** need to rerun the complete preprocessing pipeline.

The primary cleaned datasets used for Q1–Q4 analyses are generated using:

```text
19-merge_ameri_WUE_indices.py
```

and final quality-control filtering is performed using:

```text
20-WUE_cleaning_plot.py
```

Primary cleaned analysis directory:

```text
M:\Research\WUE_CUE\data_products\
```

Key cleaned analysis files:

```text
WUE_CUE_monthly_merged_indices_clean.csv
WUE_CUE_yearly_merged_indices_clean.csv
```

These files contain:
- WUE_ET
- WUEₜ
- ET partitioning metrics
- SPEI drought indices
- climate and biome classifications
- site metadata
- hydroclimatic drivers
- growing season metrics

and serve as the primary inputs for Q1–Q4 analyses.

---

# 📁 Q1–Q4 Manuscript Analysis Directories

Downstream statistical analyses, figure-generation workflows, and manuscript-specific modeling scripts are organized by research question within the repository:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/
```

Subdirectories include:

```text
Q1/
```

Near-normal hydroclimatic analyses, ET partitioning comparisons, and baseline WUE/WUEₜ assessments.

```text
Q2/
```

Short- and long-term drought response analyses using multi-timescale SPEI.

```text
Q3/
```

Site-level WUEₜ sensitivity analyses across ecosystem types, climate–biome groups, salinity gradients, and coastal regions.

```text
Q4/
```

Persistent multi-year drought analyses, spatial persistence workflows, breakpoint analyses, and regional WUEₜ decline probability assessments.

Additional workflow directories:

```text
climate_biome/
```

Climate-zone and biome-level comparison analyses.

```text
methods/
```

Supporting methodological workflows, diagnostics, validation scripts, and supplementary analyses.

---

# 🌎 Q4 Spatial Modeling and Regional Drought Analysis

Spatial drought and regional WUEₜ probability analyses for Q4 are organized within:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/Q4/spatial_model/
```

Key spatial analysis scripts include:

```text
1-coastal buffer for modis.py
```

Generates coastal MODIS spatial buffers and regional grids.

```text
2-SPEI_spatial_visual.py
```

Creates SPEI spatial visualization products and regional drought layers.

```text
5-spatial_prob_model.py
```

Builds spatial probability models for WUEₜ drought-response classification.

```text
6-spatial_prob_model_output_analysis.py
```

Processes modeled spatial probability outputs, breakpoint analyses, and regional aggregation products used in Q4 figures and interpretation. :contentReference[oaicite:0]{index=0}

---

## Spatial analysis data directories

Primary spatial drought datasets:

```text
M:\Research\WUE_CUE\spatial_SPEI\
```

Important spatial model inputs:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs\
```

Contains monthly NetCDF spatial probability outputs generated from the WUEₜ logistic drought-response models.

Key spatial outputs generated for Q4 analyses:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\
```

This directory contains:
- breakpoint analysis outputs
- regional aggregation tables
- monthly probability summaries
- time-aggregated raster products
- pre/post-breakpoint spatial comparisons
- NetCDF spatial probability layers used for manuscript figures

Important generated subdirectories include:

```text
diagnostics/
tables_breakpoint/
rasters_time_aggregated/
rasters_pre_post_breakpoint/
```

These outputs are used directly for:
- Q4 spatial persistence analyses
- regional WUEₜ decline probability mapping
- breakpoint detection
- temporal aggregation analyses
- manuscript figure generation
- supplementary spatial diagnostics

---

# ⚙️ Core Processing Workflow

| Step | Script | Main Purpose |
|---|---|---|
| 1 | `1-ameri_api.R` | Download AmeriFlux metadata and site information |
| 2 | `2-ameri_un_zip.py` | Extract AmeriFlux half-hourly datasets |
| 3 | `3-ameri_preprocess.py` | Preprocess and standardize flux observations |
| 4 | `4-Era5_point_data.py` | Download ERA5 climate reanalysis data |
| 5 | `5-merge_ameri_era5.py` | Merge AmeriFlux and ERA5 datasets |
| 6 | `6-blending_ameri_era.py` | Meteorological gap-filling and blending |
| 7 | `7-long_gaps.py` | Correct long observational gaps |
| 8 | `8a-Loop_gap_fill_gpp.R` | REddyProc gap-filling and partitioning |
| 9 | `10-data_merging.py` | Aggregate datasets and compute WUE metrics |
| 10 | `16-ET_partioning_may_2026.py` | ET partitioning workflows |
| 11 | `18-drought_indices_api.py` | Download and process SPEI drought indices |
| 12 | `19-merge_ameri_WUE_indices.py` | Merge WUE metrics with drought indices |
| 13 | `20-WUE_cleaning_plot.py` | Final quality control and outlier filtering |

---

# 📁 Outputs

Primary outputs include:
- gap-filled AmeriFlux datasets
- ET partitioning products
- WUE and WUEₜ datasets
- SPEI drought datasets
- spatial drought analyses
- publication-ready figures
- diagnostic and validation outputs
- climate–biome comparison analyses
- regional drought sensitivity assessments

> Large AmeriFlux and ERA5 raw datasets are stored on the Malone Lab server and are not fully tracked in the GitHub repository because of file size limitations.
