# 🌿 Changes in Water-Use Efficiency Signal Critical Declines Across U.S. Coastal Zones

This project investigates effect of hydroclimatic extremes on ecosystem Water Use Efficiency (WUE), transpiration-based WUE<sub>T</sub>, and evapotranspiration partitioning across U.S. coastal ecosystems using AmeriFlux observations, remote sensing products, and hydroclimatic indices.

The workflow integrates flux tower observations, ERA5 climate reanalysis, MODIS products, ET partitioning, phenology, and SPEI datasets to evaluate how hydroclimatic extremes alters ecosystem carbon–water coupling across coastal regions.

---

# 🔍 Project Objectives

- **Q1:** Quantify differences in WUE<sub>ET</sub>, WUE<sub>T</sub>, and evapotranspiration partitioning across U.S. coastal ecosystems under near-normal hydroclimatic conditions  

- **Q2:** Evaluate how short- and long-term hydroclimatic anomalies alter WUE<sub>ET</sub> and WUE<sub>T</sub> responses across coastal ecosystems using multi-timescale SPEI gradients  

- **Q3:** Assess site-level WUE<sub>T</sub> sensitivity to persistent multi-year moisture anomalies and determine whether sensitivity patterns differ among ecosystem types, climate–biome groups, salinity gradients, and coastal regions  

- **Q4:** Determine how persistent multi-year hydroclimatic anomalies influence the probability, spatial persistence, and regional patterns of WUE<sub>T</sub> decline across U.S. coastal regions  

- Integrate AmeriFlux observations, ET partitioning, ERA5 climate reanalysis, MODIS vegetation products, and SPEI datasets within a unified coastal carbon–water analysis framework  

---

# 📂 Data Directory Structure

All raw, intermediate, and processed datasets are organized on the Malone Lab server:

```text
M:\Research\WUE_CUE\
```

---

## Raw AmeriFlux ZIP archives

```text
M:\Research\WUE_CUE\ameri_data\
```

Processed using:

```text
2-ameri_un_zip.py
```

Extracted half-hourly AmeriFlux files are generated in:

```text
M:\Research\WUE_CUE\ameri_data\ameri_gaps\
```

---

## ERA5 meteorological forcing data

ERA5 processed datasets are located in:

```text
M:\Research\WUE_CUE\era5_point_data\
```

Processed using:

```text
4-Era5_point_data.py
```

---

## Final merged and cleaned datasets

Merged datasets are generated using:

```text
19-merge_ameri_WUE_indices.py
```

Main input directories:

```text
M:\Research\WUE_CUE\ameri_data\ET_partitioning\
M:\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\drought\
M:\Research\WUE_CUE\data_products\info\
```

Output directory:

```text
M:\Research\WUE_CUE\data_products\
```

Merged datasets:

```text
WUE_CUE_monthly_merged_indices.csv
WUE_CUE_yearly_merged_indices.csv
```

Final cleaned datasets are generated using:

```text
20-WUE_cleaning_plot.py
```

Cleaned datasets:

```text
WUE_CUE_monthly_merged_indices_clean.csv
WUE_CUE_yearly_merged_indices_clean.csv
```

GitHub copies of merged and cleaned datasets are available in:

```text
data/
```

Repository:

:contentReference[oaicite:1]{index=1}

---

# 📊 Q1–Q4 Manuscript Analysis Workflows

Reviewers interested only in manuscript analyses and figure-generation workflows do not need to rerun the full preprocessing pipeline.

Primary cleaned datasets used throughout the manuscript are located in:

```text
M:\Research\WUE_CUE\data_products\
```

Main downstream analysis scripts in the GitHub repository are organized in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/
```

Subdirectories:

```text
Q1/
Q2/
Q3/
Q4/
```

These folders contain manuscript-specific statistical analyses, sensitivity analyses, regional analyses, and figure-generation workflows.

---

# 🌎 Q4 Spatial Modeling and Regional Analysis

Q4 spatial-modeling scripts in the GitHub repository are located in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/Q4/spatial_model/
```

Key scripts:

```text
5-spatial_prob_model.py
```

Builds spatial probability models for WUE<sub>T</sub> response classification across coastal regions.

```text
6-spatial_prob_model_output_analysis.py
```

Processes NetCDF probability outputs, breakpoint analyses, and regional spatial products used in Q4 manuscript figures.

Primary spatial-analysis datasets are located in:

```text
M:\Research\WUE_CUE\spatial_SPEI\
```

Important NetCDF model outputs:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs\
```

Spatial outputs and breakpoint-analysis products used for Q4 figures are generated in:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\
```

GitHub scripts for spatial analyses are available in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/Q4/spatial_model/
```

---

# ⚙️ Processing Workflow

| Step | Script | Main Purpose |
|---|---|---|
| 1 | `1-ameri_api.R` | Download AmeriFlux metadata and site information |
| 2 | `2-ameri_un_zip.py` | Extract AmeriFlux half-hourly ZIP archives |
| 3 | `3-ameri_preprocess.py` | Preprocess and standardize flux observations |
| 4 | `4-Era5_point_data.py` | ERA5 data processing |
| 5 | `5-merge_ameri_era5.py` | Merge AmeriFlux and ERA5 datasets |
| 6 | `6-blending_ameri_era.py` | Meteorological gap-filling and blending |
| 7 | `7-long_gaps.py` | Correct long observational gaps |
| 8 | `8a-Loop_gap_fill_gpp.R` | REddyProc gap-filling and flux partitioning |
| 9 | `10-data_merging.py` | Aggregate datasets and compute WUE metrics |
| 10 | `16-ET_partioning_may_2026.py` | ET partitioning workflows |
| 11 | `18-drought_indices_api.py` | Download and process SPEI indices |
| 12 | `19-merge_ameri_WUE_indices.py` | Generate merged datasets before outlier removal |
| 13 | `20-WUE_cleaning_plot.py` | Generate final cleaned datasets for analyses |

---

# 📁 Outputs

Primary outputs include:
- gap-filled AmeriFlux datasets
- ET partitioning products
- WUE and WUE<sub>T</sub> datasets
- SPEI hydroclimatic datasets
- NetCDF spatial probability products
- regional spatial analyses
- publication-ready figures

> Large AmeriFlux and ERA5 datasets are stored on the Malone Lab server and are not fully tracked in the GitHub repository because of file size limitations.
