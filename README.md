# 🌿 Changes in Water-Use Efficiency Signal Critical Declines Across U.S. Coastal Zones

## Project Overview

This repository contains workflows, processed datasets, statistical analyses, and spatial modeling products used to investigate long-term changes in ecosystem water-use efficiency across U.S. coastal ecosystems using AmeriFlux observations, ERA5 climate reanalysis, ET partitioning, SPEI hydroclimatic indices, and spatial probability modeling.

---

# 🔍 Research Objectives

| Question | Objective |
|---|---|
| **Q1** | Quantify differences in WUE<sub>ET</sub>, WUE<sub>T</sub>, and evapotranspiration partitioning across U.S. coastal ecosystems under near-normal hydroclimatic conditions. |
| **Q2** | Evaluate how short- and long-term hydroclimatic anomalies alter WUE<sub>ET</sub>, WUE<sub>T</sub>, and ET partitioning behavior using multi-timescale SPEI gradients. |
| **Q3** | Assess site-level WUE<sub>T</sub> sensitivity to persistent multi-year moisture anomalies across ecosystem types, climate–biome groups, salinity gradients, and coastal regions. |
| **Q4** | Determine how persistent multi-year hydroclimatic anomalies influence probability, spatial persistence, and regional patterns of WUE<sub>T</sub> decline across U.S. coastal regions. |

---

---

# 📂 Final Datasets (Quick Access)

| Dataset | Description | Server Path | GitHub Path |
|---|---|---|---|
| `WUE_CUE_monthly_merged_indices_clean.csv` | Final cleaned monthly ecosystem dataset | `\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\` | `data_products/` |
| `WUE_CUE_yearly_merged_indices_clean.csv` | Final cleaned yearly ecosystem dataset | `\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\` | `data_products/` |
| `WUE_CUE_monthly_merged_indices.csv` | Monthly merged dataset before filtering | `\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\` | `data_products/` |
| `WUE_CUE_yearly_merged_indices.csv` | Yearly merged dataset before filtering | `\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\` | `data_products/` |

---

# 📂 Data Directory Structure

All raw, intermediate, and processed datasets are organized on the Malone Lab server:

```text
M:\Research\WUE_CUE\
```

## AmeriFlux Data

| Item | Path / Script |
|---|---|
| Raw AmeriFlux ZIP archives | `M:\Research\WUE_CUE\ameri_data\` |
| ZIP extraction workflow | `2-ameri_un_zip.py` |
| Extracted half-hourly files | `M:\Research\WUE_CUE\ameri_data\ameri_gaps\` |

---

## ERA5 Climate Data

| Item | Path / Script |
|---|---|
| ERA5 processed datasets | `M:\Research\WUE_CUE\era5_point_data\` |
| ERA5 processing workflow | `4-Era5_point_data.py` |

---

## Final Merged and Cleaned Datasets

| Item | Path / Script |
|---|---|
| Merge workflow | `19-merge_ameri_WUE_indices.py` |
| Main ET partitioning input | `M:\Research\WUE_CUE\ameri_data\ET_partitioning\` |
| SPEI drought input | `M:\Research\WUE_CUE\drivers\ameri_drivers\PET_drought\drought\` |
| Site metadata input | `M:\Research\WUE_CUE\data_products\info\` |
| Output directory | `M:\Research\WUE_CUE\data_products\` |
| Monthly merged dataset | `WUE_CUE_monthly_merged_indices.csv` |
| Yearly merged dataset | `WUE_CUE_yearly_merged_indices.csv` |
| Final cleaning workflow | `20-WUE_cleaning_plot.py` |
| Monthly cleaned dataset | `WUE_CUE_monthly_merged_indices_clean.csv` |
| Yearly cleaned dataset | `WUE_CUE_yearly_merged_indices_clean.csv` |
| GitHub dataset directory | `data_products/` |

---

# ⚙️ Workflow Overview

```text
AmeriFlux + ERA5 Integration
            ↓
ET Partitioning + WUE Metrics
            ↓
Quality Control + Outlier Filtering
            ↓
SPEI Sensitivity Analysis
            ↓
Spatial Probability Modeling
            ↓
Breakpoint + Regional Aggregation Analysis
```

---

# 🌎 Q1–Q4 Manuscript Analysis Workflows

Main manuscript-analysis scripts in the GitHub repository are organized in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/
```

| Directory | Purpose |
|---|---|
| `Q1/` | Baseline ecosystem variability analyses |
| `Q2/` | Multi-timescale SPEI response analyses |
| `Q3/` | Long-term SPEI-48 sensitivity analyses |
| `Q4/` | Spatial probability and regional analyses |

Primary manuscript datasets are located in:

```text
M:\Research\WUE_CUE\data_products\
```

---

# 🌍 Spatial Modeling Workflow (Q4)

Spatial-modeling scripts in the GitHub repository are located in:

```text
function/WUE_paper_figures/draft_5_spei/may_figures/Q4/spatial_model/
```

## Spatial Probability Projection

| Item | Description |
|---|---|
| `5-spatial_prob_model.py` | Builds spatial probability models for WUE<sub>T</sub> response classification using gridded SPEI datasets |
| Main outputs | `WUE_predictions_linearlogistic_YYYYMM.nc` |
| Additional outputs | Monthly and yearly summary tables |

---

## Temporal Aggregation and Breakpoint Analysis

| Item | Description |
|---|---|
| `6-spatial_prob_model_output_analysis.py` | Processes NetCDF outputs into breakpoint analyses, regional summaries, and manuscript figure products |
| Spatial analysis directory | `M:\Research\WUE_CUE\spatial_SPEI\` |
| NetCDF outputs | `M:\Research\WUE_CUE\spatial_SPEI\logistic_model\netcdf_outputs\` |
| Breakpoint outputs | `M:\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\` |

---

# ⚙️ Core Processing Workflow

| Step | Script | Main Purpose |
|---|---|---|
| 1 | `1-ameri_api.R` | Download AmeriFlux metadata |
| 2 | `2-ameri_un_zip.py` | Extract AmeriFlux ZIP archives |
| 3 | `3-ameri_preprocess.py` | Preprocess flux observations |
| 4 | `4-Era5_point_data.py` | ERA5 data processing |
| 5 | `5-merge_ameri_era5.py` | Merge AmeriFlux and ERA5 datasets |
| 6 | `6-blending_ameri_era.py` | Meteorological gap-filling |
| 7 | `7-long_gaps.py` | Correct long observational gaps |
| 8 | `8a-Loop_gap_fill_gpp.R` | REddyProc gap-filling and partitioning |
| 9 | `10-data_merging.py` | Compute WUE metrics |
| 10 | `16-ET_partioning_may_2026.py` | ET partitioning workflows |
| 11 | `18-drought_indices_api.py` | Process SPEI indices |
| 12 | `19-merge_ameri_WUE_indices.py` | Generate merged datasets |
| 13 | `20-WUE_cleaning_plot.py` | Generate final cleaned datasets |

---

# 📁 Repository Structure

```text
Water_WUE/
│
├── data_products/
├── function/
│   ├── WUE_paper_figures/
│   ├── climate_biome/
│   ├── methods/
│   └── spatial_model/
│
├── spatial_SPEI/
├── figures/
├── diagnostics/
└── results/
```

---

# 📦 Primary Outputs

- ET partitioning products
- merged WUE datasets
- cleaned ecosystem datasets
- SPEI sensitivity analyses
- NetCDF spatial probability products
- breakpoint analyses
- regional aggregation products
- publication-ready figures

---

# 👤 Contact

Ammara Talib  
Yale University  
Malone Lab  
