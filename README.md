# Coastal Ecosystems Show Regionally Distinct Resilience to Hydroclimatic Extremes Across U.S. Coastlines

## Project Overview

This repository contains the data-processing and manuscript-analysis workflows used to investigate ecosystem water-use efficiency and hydroclimatic resilience across U.S. coastal ecosystems.

The study uses eddy-covariance observations from 64 AmeriFlux sites representing upland, freshwater, and saline ecosystems. AmeriFlux observations are combined with ERA5 meteorological data, evapotranspiration partitioning, multi-timescale Standardized Precipitation Evapotranspiration Index (SPEI), generalized additive models, and spatial SPEI products.

The analysis compares alternative water-use efficiency metrics, quantifies site-level stability and drought resilience, evaluates nonlinear WUE<sub>T</sub> responses across SPEI timescales and coastlines, and spatially evaluates coastal WUE<sub>T</sub> responses from 2000–2025.

---

## Workflow Overview

```text
AmeriFlux Observations + ERA5 Meteorology
                    ↓
Meteorological Gap Filling and GPP Processing
                    ↓
ET Partitioning and WUE Metric Calculation
                    ↓
SPEI Calculation Across Multiple Timescales
                    ↓
Dataset Merging, Quality Control, and Filtering
                    ↓
Comparison of WUE Metrics Along the T:ET Gradient
                    ↓
Site-Level Stability, Plasticity, Resistance, and Recovery
                    ↓
Ecosystem and Coast-Specific GAM Analysis
                    ↓
Ecological Threshold Identification
                    ↓
Spatial Upscaling Across U.S. Coastlines, 2000–2025
```

---

## Research Questions

| Question | Objective                                                                                                                                                                                        |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Q1**   | Determine how WUE<sub>ET</sub>, WUE<sub>E</sub>, and WUE<sub>T</sub> differ along the transpiration fraction gradient under near-normal hydroclimatic conditions.                                |
| **Q2**   | Evaluate how WUE<sub>T</sub> performance varies across coastal ecosystems using site-level stability, plasticity, drought resistance, and post-drought recovery metrics.                         |
| **Q3**   | Identify when meteorological extremes produce ecologically meaningful WUE<sub>T</sub> responses across SPEI timescales and U.S. coastlines and spatially upscale these responses from 2000–2025. |

---

## Main Workflow Directories

### Dataset Creation

Scripts used to download, process, merge, gap-fill, partition, and clean the AmeriFlux, ERA5, ET, WUE, and SPEI datasets are located in:

```text
M:\Research\WUE_CUE\Water_WUE\Data_processing_workflow\
```

### Manuscript Analyses

Scripts used for the Q1–Q3 manuscript analyses, statistical models, spatial analyses, and manuscript figures are located in:

```text
M:\Research\WUE_CUE\Water_WUE\manuscripts_Analysis\
```

---

## Final Datasets

| Dataset                                    | Purpose                                                    |
| ------------------------------------------ | ---------------------------------------------------------- |
| `WUE_CUE_monthly_merged_indices.csv`       | Monthly merged dataset before final filtering              |
| `WUE_CUE_yearly_merged_indices.csv`        | Yearly merged dataset before final filtering               |
| `WUE_CUE_monthly_merged_indices_clean.csv` | Final cleaned monthly dataset used for manuscript analyses |
| `WUE_CUE_yearly_merged_indices_clean.csv`  | Final cleaned yearly dataset                               |

Final dataset directory:

```text
M:\Research\WUE_CUE\data_products\
```

GitHub dataset directory:

```text
data_products/
```

Repository:

```text
https://github.com/Malone-Disturbance-Ecology-Lab/Water_WUE
```

---

## Core Data-Processing Scripts

| Step | Script                          | Purpose                                                           |
| ---- | ------------------------------- | ----------------------------------------------------------------- |
| 1    | `1-ameri_api.R`                 | Download AmeriFlux metadata                                       |
| 2    | `2-ameri_un_zip.py`             | Extract AmeriFlux ZIP archives                                    |
| 3    | `3-ameri_preprocess.py`         | Preprocess flux and meteorological observations                   |
| 4    | `4-Era5_point_data.py`          | Process ERA5 meteorological data                                  |
| 5    | `5-merge_ameri_era5.py`         | Merge AmeriFlux and ERA5 datasets                                 |
| 6    | `6-blending_ameri_era.py`       | Fill missing meteorological observations                          |
| 7    | `7-long_gaps.py`                | Identify and manage long observational gaps                       |
| 8    | `8a-Loop_gap_fill_gpp.R`        | Perform REddyProc gap filling and flux partitioning               |
| 9    | `10-data_merging.py`            | Merge processed site-level variables                              |
| 10   | `16-ET_partioning_may_2026.py`  | Partition ET and calculate WUE metrics                            |
| 11   | `18-drought_indices_api.py`     | Calculate and process SPEI indices                                |
| 12   | `19-merge_ameri_WUE_indices.py` | Generate merged monthly and yearly datasets                       |
| 13   | `20-WUE_cleaning_plot.py`       | Apply final quality control and filtering                         |
| 14   | `24-EDI_spatial.py`             | Spatially upscale coast-specific SPEI-3 WUE<sub>T</sub> responses |

---

## Repository Structure

```text
Water_WUE/
│
├── Data_processing_workflow/
├── manuscripts_Analysis/
├── data_products/
├── spatial_SPEI/
├── figures/
├── diagnostics/
└── results/
```

---

## Primary Outputs

* Processed AmeriFlux and ERA5 time series
* ET partitioning products
* WUE<sub>ET</sub>, WUE<sub>E</sub>, and WUE<sub>T</sub> datasets
* Monthly and yearly cleaned datasets
* Site-level WUE<sub>T</sub> performance metrics
* Multi-timescale SPEI response analyses
* Ecosystem- and coast-specific GAM analyses
* Spatial WUE<sub>T</sub> response products
* Publication-ready manuscript figures

---

## Study Citation

Manuscript in preparation:

**Talib, A., Wang, J., Zhang, B., and Malone, S.**
*Coastal Ecosystems Show Regionally Distinct Resilience to Hydroclimatic Extremes Across U.S. Coastlines.*

Yale School of the Environment, Yale University.

---

## Contact

**Ammara Talib**
Yale School of the Environment
Yale University
Malone Lab
