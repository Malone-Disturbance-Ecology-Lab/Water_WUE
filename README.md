# 🌿 Water_WUE: Coastal Ecosystem Water Use Efficiency and Drought Analysis

This project investigates spatial and temporal variability in ecosystem Water Use Efficiency (WUE), transpiration-based WUE_T, and evapotranspiration partitioning across U.S. coastal ecosystems using AmeriFlux observations, remote sensing products, and hydroclimatic drought indices. The workflow integrates flux tower observations, ERA5 climate reanalysis, MODIS products, ET partitioning, phenology, salinity, and SPEI drought datasets to evaluate how prolonged hydroclimatic stress alters ecosystem carbon–water coupling across coastal regions.

---

# 🔍 Project Objectives

- **Q1:** Quantify differences in WUE_ET, WUE_T, and evapotranspiration partitioning across U.S. coastal ecosystems under near-normal hydroclimatic conditions  

- **Q2:** Evaluate how short- and long-term hydroclimatic anomalies alter WUE_ET and WUE_T responses across coastal ecosystems using multi-timescale SPEI  

- **Q3:** Assess site-level WUE_T sensitivity to persistent multi-year moisture anomalies and determine whether sensitivity patterns differ among ecosystem types, climate–biome groups, salinity gradients, and coastal regions  

- **Q4:** Determine how persistent multi-year drought influences the probability, spatial persistence, and regional patterns of WUE_T decline across U.S. coastal regions  

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

## Final merged datasets before outlier removal

Generated using:

```text
19-merge_ameri_WUE_indices.py
```

This script merges:
- ET partitioning outputs
- AmeriFlux processed datasets
- ERA5 meteorological drivers
- SPEI drought indices
- climate and biome metadata
- growing season datasets

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

Output directory:

```text
M:\Research\WUE_CUE\data_products\
```

Key output files:

```text
WUE_CUE_monthly_merged_indices.csv
WUE_CUE_yearly_merged_indices.csv
```

---

## Final cleaned datasets after outlier removal

Generated using:

```text
20-WUE_cleaning_plot.py
```

Final cleaned outputs:

```text
WUE_CUE_monthly_merged_indices_clean.csv
WUE_CUE_yearly_merged_indices_clean.csv
```

These files are the primary datasets used for Q1–Q4 analyses and manuscript figure generation.

---

# 📊 Q1–Q4 Manuscript Analysis Workflows

Reviewers interested only in manuscript analyses and figure-generation workflows do not need to rerun the full preprocessing pipeline.

The primary cleaned datasets used throughout the manuscript are located in:

```text
M:\Research\WUE_CUE\data_products\
```

Key cleaned files:

```text
WUE_CUE_monthly_merged_indices_clean.csv
WUE_CUE_yearly_merged_indices_clean.csv
```

Downstream analysis scripts are organized by research question within the repository:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/
```

Subdirectories include:

```text
Q1/
Q2/
Q3/
Q4/
```

Additional directories:

```text
climate_biome/
methods/
```

These folders contain manuscript-specific:
- statistical analyses
- sensitivity analyses
- climate–biome comparisons
- regional analyses
- figure-generation workflows
- supplementary diagnostics

---

# 🌎 Q4 Spatial Modeling and Regional Drought Analysis

Q4 spatial modeling workflows are located in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/Q4/spatial_model/
```

Key scripts include:

```text
5-spatial_prob_model.py
```

Builds spatial probability models for WUE_T drought-response classification.

```text
6-spatial_prob_model_output_analysis.py
```

Processes NetCDF spatial probability outputs, breakpoint analyses, and regional aggregation products used for Q4 manuscript figures and spatial persistence analyses. :contentReference[oaicite:0]{index=0}

Primary spatial analysis data directory:

```text
M:\Research\WUE_CUE\spatial_SPEI\
```

Important NetCDF model outputs used in Q4 analyses:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs\
```

Spatial outputs and breakpoint-analysis products used for manuscript figures are generated in:

```text
M:\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\
```

These directories contain:
- NetCDF spatial probability outputs
- breakpoint-analysis summaries
- regional aggregation tables
- raster products
- spatial persistence outputs
- time-aggregated regional analyses

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
| 13 | `20-WUE_cleaning_plot.py` | Final outlier removal and cleaned dataset generation |

---

# 📁 Outputs

Primary outputs include:
- gap-filled AmeriFlux datasets
- ET partitioning products
- WUE and WUE_T datasets
- SPEI drought datasets
- spatial drought analyses
- NetCDF probability products
- publication-ready figures
- climate–biome comparison analyses
- regional drought sensitivity assessments

> Large AmeriFlux and ERA5 raw datasets are stored on the Malone Lab server and are not fully tracked in the GitHub repository because of file size limitations.
