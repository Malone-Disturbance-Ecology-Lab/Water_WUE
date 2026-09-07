# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 12:52:26 2026

@author: ammar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Q3_Delta_exclusion_sensitivity_diagnostic.py

Standalone diagnostic script that replicates the exact Q3 coast-threshold GAM
analysis (from Q3.WUE_sensitivity_SPEI.py) for two scenarios:
    FULL      – complete dataset
    DELTA_EXCLUDED – removes 10 Sacramento–San Joaquin Delta sites

All results are printed to the console. No permanent files are saved.
Temporary files (R scripts, CSVs) are automatically deleted.

The GAM is fitted in R using mgcv, exactly as in the production workflow.
"""

import os
import sys
import subprocess
import tempfile
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Q3 DELTA EXCLUSION SENSITIVITY DIAGNOSTIC")
print("="*80)

# ============================================================================
# CONSTANTS – EXACTLY AS IN Q3.WUE_sensitivity_SPEI.py
# ============================================================================

SPEI_COLS = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
SELECTED_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]

# ---------------------------------------------------------------------------
# Coast classifier (exact copy from production)
# ---------------------------------------------------------------------------
EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9),
        (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2),
        (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0),
        (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7),
        (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1),
        (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6),
        (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3),
        (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9),
        (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)
    ],
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
    import math
    lat0 = math.radians(lat)
    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y])
    p = np.array([0.0, 0.0])
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)
    ab = b - a
    denom = np.dot(ab, ab)
    if denom == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, np.dot(p - a, ab) / denom))
    return float(np.linalg.norm(p - (a + t * ab)))

def _distance_to_polyline(lat, lon, polyline):
    distances = []
    for i in range(len(polyline) - 1):
        d = _seg_dist(lat, lon, *polyline[i], *polyline[i+1])
        distances.append(d)
    return min(distances) if distances else float('inf')

def classify_coast_region(lat, long):
    if pd.isna(lat) or pd.isna(long):
        return "Other/Check"
    if lat > 50:
        return "AK Coast"
    lat, lon = float(lat), float(long)
    dists = {
        name: _distance_to_polyline(lat, lon, poly)
        for name, poly in COAST_POLYLINES.items()
    }
    return min(dists, key=dists.get)

# ---------------------------------------------------------------------------
# Delta sites – exactly 10 documented sites
# ---------------------------------------------------------------------------
DELTA_SITES = {
    "US-Dmg", "US-Myb", "US-Snd", "US-Sne", "US-Tw1",
    "US-Tw2", "US-Tw3", "US-Tw4", "US-Tw5", "US-Twt"
}

# ============================================================================
# INPUT DATA PATH – same as production
# ============================================================================
input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# ============================================================================
# STEP 1 — REPRODUCE ORIGINAL Q3 MODEL-BASE DATA
# ============================================================================
print("\nSTEP 1: Loading and filtering data (production copy)")

monthly = pd.read_csv(input_file)

# Trim whitespace from character columns
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

# Required columns (same as production)
required_cols = ["site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra"] + SPEI_COLS

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Filtering – identical to production
model_base = monthly.dropna(subset=required_cols).copy()
model_base = model_base[
    model_base['site_name'].notna() &
    (model_base['site_name'] != "") &
    model_base['water_class'].notna() &
    (model_base['water_class'] != "") &
    model_base['water_class'].isin(ECOSYSTEM_CLASSES) &
    np.isfinite(model_base['lat']) &
    np.isfinite(model_base['long']) &
    np.isfinite(model_base['WUE_tra'])
].copy()

# Add coast region (same function)
model_base['coast_region'] = model_base.apply(
    lambda row: classify_coast_region(row['lat'], row['long']), axis=1
)

# Counts for FULL
n_full_rows = len(model_base)
n_full_sites = model_base['site_name'].nunique()
n_full_pacific = len(model_base[model_base['coast_region'] == "Pacific Coast"])
n_full_pacific_sites = model_base[model_base['coast_region'] == "Pacific Coast"]['site_name'].nunique()
n_delta_in_pacific = len(model_base[model_base['site_name'].isin(DELTA_SITES)])
n_delta_sites = len(set(model_base['site_name']) & DELTA_SITES)

print(f"  FULL model-base rows    : {n_full_rows} (expected 1,872)")
print(f"  FULL sites              : {n_full_sites} (expected 64)")
print(f"  Pacific Coast rows      : {n_full_pacific} (expected 670)")
print(f"  Pacific Coast sites     : {n_full_pacific_sites} (expected 23)")
print(f"  Delta sites in Pacific  : {n_delta_in_pacific} rows, {n_delta_sites} sites")

# Validate against expected
if n_full_rows != 1872:
    raise ValueError(f"FULL rows mismatch: {n_full_rows} vs 1872")
if n_full_sites != 64:
    raise ValueError(f"FULL sites mismatch: {n_full_sites} vs 64")
if n_full_pacific != 670:
    raise ValueError(f"Pacific rows mismatch: {n_full_pacific} vs 670")
if n_full_pacific_sites != 23:
    raise ValueError(f"Pacific sites mismatch: {n_full_pacific_sites} vs 23")
if n_delta_in_pacific != 400:
    raise ValueError(f"Delta rows in Pacific mismatch: {n_delta_in_pacific} vs 400")
if n_delta_sites != 10:
    raise ValueError(f"Delta sites count mismatch: {n_delta_sites} vs 10")

print("  Counts match production expectations.")

# ============================================================================
# STEP 2 — CREATE TWO SCENARIOS
# ============================================================================
print("\nSTEP 2: Creating FULL and DELTA_EXCLUDED scenarios")

# Scenario 1: FULL (already model_base)
full_data = model_base.copy()

# Scenario 2: DELTA_EXCLUDED (remove rows with site_name in DELTA_SITES)
delta_excluded_data = model_base[~model_base['site_name'].isin(DELTA_SITES)].copy()

# Validate counts for DELTA_EXCLUDED
n_delta_excl_rows = len(delta_excluded_data)
n_delta_excl_sites = delta_excluded_data['site_name'].nunique()
n_delta_excl_pacific = len(delta_excluded_data[delta_excluded_data['coast_region'] == "Pacific Coast"])
n_delta_excl_pacific_sites = delta_excluded_data[delta_excluded_data['coast_region'] == "Pacific Coast"]['site_name'].nunique()

print(f"  DELTA_EXCLUDED rows    : {n_delta_excl_rows} (expected 1,472)")
print(f"  DELTA_EXCLUDED sites   : {n_delta_excl_sites} (expected 54)")
print(f"  Pacific rows           : {n_delta_excl_pacific} (expected 270)")
print(f"  Pacific sites          : {n_delta_excl_pacific_sites} (expected 13)")

if n_delta_excl_rows != 1472:
    raise ValueError(f"DELTA_EXCLUDED rows mismatch: {n_delta_excl_rows} vs 1472")
if n_delta_excl_sites != 54:
    raise ValueError(f"DELTA_EXCLUDED sites mismatch: {n_delta_excl_sites} vs 54")
if n_delta_excl_pacific != 270:
    raise ValueError(f"DELTA_EXCLUDED Pacific rows mismatch: {n_delta_excl_pacific} vs 270")
if n_delta_excl_pacific_sites != 13:
    raise ValueError(f"DELTA_EXCLUDED Pacific sites mismatch: {n_delta_excl_pacific_sites} vs 13")

print("  Counts match expected after exclusion.")

# ============================================================================
# STEP 3 — REBUILD LONG-FORMAT DATA FOR EACH SCENARIO
# ============================================================================
print("\nSTEP 3: Reshaping to long format for both scenarios")

def build_long(df, scenario_name):
    """Identical to production melt + factor setup"""
    # Keep all columns except SPEI columns
    id_vars = [c for c in df.columns if c not in SPEI_COLS]
    long = pd.melt(
        df,
        id_vars=id_vars,
        value_vars=SPEI_COLS,
        var_name='SPEI_timescale',
        value_name='SPEI_value'
    )
    # Filter finite
    long = long[np.isfinite(long['SPEI_value']) & np.isfinite(long['WUE_tra'])]

    # Add factors (as in production)
    long['SPEI_timescale'] = pd.Categorical(long['SPEI_timescale'], categories=SPEI_COLS)
    long['site_name'] = long['site_name'].astype('category')
    long['water_class'] = pd.Categorical(long['water_class'], categories=ECOSYSTEM_CLASSES)
    long['coast_region'] = pd.Categorical(long['coast_region'], categories=COAST_REGION_LEVELS)
    long['month_f'] = pd.Categorical(long['month'], categories=range(1, 13))
    long['SPEI_z'] = (long['SPEI_value'] - long['SPEI_value'].mean()) / long['SPEI_value'].std()
    long['spei_ecosystem'] = long['SPEI_timescale'].astype(str) + "__" + long['water_class'].astype(str)
    long['spei_coast'] = long['SPEI_timescale'].astype(str) + "__" + long['coast_region'].astype(str)
    long['WUE_T'] = long['WUE_tra']   # alias

    # Additional counts for Pacific (these are long-format rows, not original observations)
    # We do not need these for the final table, but we print them for diagnostics.
    pacific_long = long[long['coast_region'] == "Pacific Coast"]
    pacific_sites = pacific_long['site_name'].nunique()
    pacific_rows = len(pacific_long)
    # By ecosystem
    eco_counts = pacific_long.groupby('water_class').size().to_dict()

    return long, pacific_rows, pacific_sites, eco_counts

full_long, full_pac_rows, full_pac_sites, full_eco = build_long(full_data, "FULL")
delta_long, delta_pac_rows, delta_pac_sites, delta_eco = build_long(delta_excluded_data, "DELTA_EXCLUDED")

print(f"\nFULL long format:")
print(f"  Total rows     : {len(full_long)}")
print(f"  Pacific rows   : {full_pac_rows} (long-format; original Pacific observations = {n_full_pacific})")
print(f"  Pacific sites  : {full_pac_sites} (unique sites, matches {n_full_pacific_sites})")
for eco in ECOSYSTEM_CLASSES:
    print(f"    {eco}: {full_eco.get(eco, 0)} rows")

print(f"\nDELTA_EXCLUDED long format:")
print(f"  Total rows     : {len(delta_long)}")
print(f"  Pacific rows   : {delta_pac_rows} (long-format; original Pacific observations = {n_delta_excl_pacific})")
print(f"  Pacific sites  : {delta_pac_sites} (unique sites, matches {n_delta_excl_pacific_sites})")
for eco in ECOSYSTEM_CLASSES:
    print(f"    {eco}: {delta_eco.get(eco, 0)} rows")

# ----------------------------------------------------------------------------
# Observed Pacific SPEI-3 support (for later reference)
# ----------------------------------------------------------------------------
def get_pacific_spei3_support(long_df):
    sub = long_df[(long_df['coast_region'] == "Pacific Coast") & (long_df['SPEI_timescale'] == "SPEI_3")]
    vals = sub['SPEI_value'].dropna()
    if len(vals) == 0:
        return None
    return vals

full_spei3_vals = get_pacific_spei3_support(full_long)
delta_spei3_vals = get_pacific_spei3_support(delta_long)

# ============================================================================
# STEP 4 — PREPARE TEMPORARY DIRECTORY AND R SCRIPT
# ============================================================================
print("\nSTEP 4: Preparing R script for GAM fitting and predictions")

# We'll use a single temporary directory
temp_dir = tempfile.TemporaryDirectory(prefix="Q3_delta_diag_")
temp_path = temp_dir.name
print(f"  Temporary directory: {temp_path}")

# Write the long data to CSV files for R
full_csv = os.path.join(temp_path, "full_long.csv")
delta_csv = os.path.join(temp_path, "delta_long.csv")
full_long.to_csv(full_csv, index=False)
delta_long.to_csv(delta_csv, index=False)

# Precompute R-safe paths for compatibility
full_csv_r = full_csv.replace("\\", "/")
delta_csv_r = delta_csv.replace("\\", "/")
temp_path_r = temp_path.replace("\\", "/")

# Determine common observed support for fixed SPEI values (Pacific, SPEI_3)
common_min = max(full_spei3_vals.min(), delta_spei3_vals.min())
common_max = min(full_spei3_vals.max(), delta_spei3_vals.max())
requested_spei = [-1.0, -1.5, -2.0, -2.5, -3.0]
valid_requested = [x for x in requested_spei if common_min <= x <= common_max]
print(f"  Common observed SPEI-3 range: [{common_min:.2f}, {common_max:.2f}]")
print(f"  Valid fixed SPEI values for comparison: {valid_requested}")

# Prepare R script content with fixes
r_script_content = f'''
# Temporary R script for Delta exclusion sensitivity diagnostic
# This replicates the coast-threshold GAM from Q3.WUE_sensitivity_SPEI.py

library(mgcv)
library(dplyr)

# ---- helper to write results ----
write_table <- function(x, file) {{
    write.csv(x, file, row.names = FALSE)
}}

# ---- load data ----
full_data <- read.csv("{full_csv_r}")
delta_data <- read.csv("{delta_csv_r}")

# ---- function to prepare factors (same as production) ----
prepare_data <- function(data) {{
    data$site_name <- as.factor(data$site_name)
    data$water_class <- factor(data$water_class, levels = c("Upland", "Freshwater", "Saline"))
    data$coast_region <- factor(data$coast_region, levels = c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"))
    data$month_f <- as.factor(data$month_f)
    data$SPEI_timescale <- factor(data$SPEI_timescale,
                                  levels = c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"))
    data$spei_coast <- interaction(data$SPEI_timescale, data$coast_region, sep = "__", drop = TRUE)
    # Ensure spei_coast factor levels match
    data$spei_coast <- factor(data$spei_coast, levels = levels(data$spei_coast))
    return(data)
}}

full_data <- prepare_data(full_data)
delta_data <- prepare_data(delta_data)

# ---- function to fit model and extract results ----
run_analysis <- function(data, scenario_name, out_prefix) {{
    cat("\\nFitting model for", scenario_name, "\\n")
    
    # Fit the coast-threshold GAM
    model <- gam(
        WUE_T ~
            SPEI_timescale * coast_region +
            water_class +
            month_f +
            s(SPEI_value, by = spei_coast, k = 6) +
            s(site_name, bs = "re"),
        data = data,
        method = "ML",
        select = TRUE
    )
    
    # ---- Model summary ----
    smry <- summary(model)
    model_info <- data.frame(
        scenario = scenario_name,
        AIC = AIC(model),
        R_sq = smry$r.sq,
        dev_expl = smry$dev.expl,
        n = nrow(data)
    )
    write_table(model_info, paste0(out_prefix, "_model_summary.csv"))
    
    # ---- Smooth terms table ----
    s_tab <- as.data.frame(smry$s.table)
    # Robust p-value column detection
    p_col <- grep("^p", names(s_tab), value = TRUE)[1]
    if (is.na(p_col) || length(p_col) == 0) {{
        stop("Could not identify p-value column in GAM smooth table.")
    }}
    s_tab$p_value <- s_tab[[p_col]]
    s_tab$smooth_term <- rownames(s_tab)
    rownames(s_tab) <- NULL
    s_tab <- s_tab[, c("smooth_term", "edf", "Ref.df", "F", "p_value")]
    s_tab$scenario <- scenario_name
    write_table(s_tab, paste0(out_prefix, "_smooth_terms.csv"))
    
    # ---- Predictions for Pacific Coast, SPEI_3, Upland, July ----
    # Use same SPEI sequence as production: quantile(0.02, 0.98, length=120)
    spei_seq <- seq(quantile(data$SPEI_value, 0.02, na.rm=TRUE),
                    quantile(data$SPEI_value, 0.98, na.rm=TRUE),
                    length.out = 120)
    
    pred_grid <- expand.grid(
        SPEI_timescale = factor("SPEI_3", levels = levels(data$SPEI_timescale)),
        coast_region = factor("Pacific Coast", levels = levels(data$coast_region)),
        water_class = factor("Upland", levels = levels(data$water_class)),
        month_f = factor("7", levels = levels(data$month_f)),
        SPEI_value = spei_seq,
        site_name = levels(data$site_name)[1]  # dummy site
    )
    pred_grid$spei_coast <- interaction(pred_grid$SPEI_timescale, pred_grid$coast_region, sep = "__", drop = TRUE)
    pred_grid$spei_coast <- factor(pred_grid$spei_coast, levels = levels(data$spei_coast))
    
    pred <- predict(model, newdata = pred_grid, type = "link", se.fit = TRUE, exclude = "s(site_name)")
    pred_grid$predicted_WUE_T <- as.numeric(pred$fit)
    pred_grid$predicted_se <- as.numeric(pred$se.fit)
    pred_grid$predicted_lower <- pred_grid$predicted_WUE_T - 1.96 * pred_grid$predicted_se
    pred_grid$predicted_upper <- pred_grid$predicted_WUE_T + 1.96 * pred_grid$predicted_se
    
    # Near-normal baseline (mean of predictions for SPEI in [-1,1])
    baseline_data <- pred_grid[pred_grid$SPEI_value >= -1 & pred_grid$SPEI_value <= 1, ]
    baseline <- mean(baseline_data$predicted_WUE_T, na.rm = TRUE)
    
    pred_grid$pct_change <- 100 * (pred_grid$predicted_WUE_T - baseline) / baseline
    pred_grid$pct_lower <- 100 * (pred_grid$predicted_lower - baseline) / baseline
    pred_grid$pct_upper <- 100 * (pred_grid$predicted_upper - baseline) / baseline
    
    # Save predictions
    pred_out <- pred_grid[, c("SPEI_value", "predicted_WUE_T", "predicted_se",
                              "predicted_lower", "predicted_upper", "pct_change",
                              "pct_lower", "pct_upper")]
    write_table(pred_out, paste0(out_prefix, "_predictions.csv"))
    
    # ---- Dry-side thresholds (5%, 10%, 20%) ----
    thresholds <- c(5, 10, 20)
    threshold_results <- data.frame()
    for (thr in thresholds) {{
        # Decrease (pct <= -thr) in dry side (SPEI < -1)
        dec <- pred_grid[pred_grid$SPEI_value < -1 & pred_grid$pct_change <= -thr, ]
        thr_dec <- if (nrow(dec) > 0) max(dec$SPEI_value, na.rm = TRUE) else NA
        # Increase (pct >= thr) in dry side
        inc <- pred_grid[pred_grid$SPEI_value < -1 & pred_grid$pct_change >= thr, ]
        thr_inc <- if (nrow(inc) > 0) max(inc$SPEI_value, na.rm = TRUE) else NA
        
        temp <- data.frame(
            scenario = scenario_name,
            threshold_pct = thr,
            direction = c("decrease", "increase"),
            SPEI_threshold = c(thr_dec, thr_inc)
        )
        threshold_results <- rbind(threshold_results, temp)
    }}
    write_table(threshold_results, paste0(out_prefix, "_thresholds.csv"))
    
    # ---- Fixed SPEI predictions for requested values ----
    fixed_spei <- c({",".join(map(str, valid_requested))})
    if (length(fixed_spei) > 0) {{
        fixed_grid <- expand.grid(
            SPEI_timescale = factor("SPEI_3", levels = levels(data$SPEI_timescale)),
            coast_region = factor("Pacific Coast", levels = levels(data$coast_region)),
            water_class = factor("Upland", levels = levels(data$water_class)),
            month_f = factor("7", levels = levels(data$month_f)),
            SPEI_value = fixed_spei,
            site_name = levels(data$site_name)[1]
        )
        fixed_grid$spei_coast <- interaction(fixed_grid$SPEI_timescale, fixed_grid$coast_region, sep = "__", drop = TRUE)
        fixed_grid$spei_coast <- factor(fixed_grid$spei_coast, levels = levels(data$spei_coast))
        
        pred_fixed <- predict(model, newdata = fixed_grid, type = "link", se.fit = TRUE, exclude = "s(site_name)")
        fixed_grid$predicted_WUE_T <- as.numeric(pred_fixed$fit)
        fixed_grid$predicted_se <- as.numeric(pred_fixed$se.fit)
        fixed_grid$predicted_lower <- fixed_grid$predicted_WUE_T - 1.96 * fixed_grid$predicted_se
        fixed_grid$predicted_upper <- fixed_grid$predicted_WUE_T + 1.96 * fixed_grid$predicted_se
        
        # Baseline is the same near-normal mean computed above
        fixed_grid$baseline_WUE_T <- baseline
        fixed_grid$pct_change <- 100 * (fixed_grid$predicted_WUE_T - baseline) / baseline
        fixed_grid$pct_lower <- 100 * (fixed_grid$predicted_lower - baseline) / baseline
        fixed_grid$pct_upper <- 100 * (fixed_grid$predicted_upper - baseline) / baseline
        
        fixed_out <- fixed_grid[, c("SPEI_value", "predicted_WUE_T", "predicted_lower",
                                    "predicted_upper", "baseline_WUE_T", "pct_change")]
        write_table(fixed_out, paste0(out_prefix, "_fixed_predictions.csv"))
    }}
    
    # ---- Return baseline for later use ----
    return(list(baseline = baseline, model = model))
}}

# ---- Run for FULL ----
full_res <- run_analysis(full_data, "FULL", file.path("{temp_path_r}", "full"))
# ---- Run for DELTA_EXCLUDED ----
delta_res <- run_analysis(delta_data, "DELTA_EXCLUDED", file.path("{temp_path_r}", "delta"))

cat("\\nR analysis complete.\\n")
'''

# Write R script
r_script_path = os.path.join(temp_path, "run_diagnostic.R")
with open(r_script_path, 'w') as f:
    f.write(r_script_content)

# ============================================================================
# STEP 5 — RUN R SCRIPT
# ============================================================================
print("\nSTEP 5: Running R script (mgcv GAM fitting)")

# Locate R executable
r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if r_path not in os.environ['PATH']:
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + r_path

result = subprocess.run(
    ["Rscript", r_script_path],
    capture_output=True,
    text=True
)

if result.returncode != 0:
    print("R script failed:")
    print(result.stderr)
    sys.exit(1)
else:
    print("R script completed successfully.")
    if result.stdout:
        print(result.stdout)

# ============================================================================
# STEP 6 — READ AND PRINT RESULTS
# ============================================================================
print("\n" + "="*80)
print("DIAGNOSTIC RESULTS")
print("="*80)

# Helper to read CSV if exists
def read_csv_safe(path):
    if os.path.exists(path):
        return pd.read_csv(path)
    else:
        return None

# ---- 6a. Model summaries ----
full_summary = read_csv_safe(os.path.join(temp_path, "full_model_summary.csv"))
delta_summary = read_csv_safe(os.path.join(temp_path, "delta_model_summary.csv"))

print("\nOVERALL MODEL INFORMATION")
print("-------------------------")
if full_summary is not None:
    print(f"FULL:            n={full_summary['n'].iloc[0]}, R²={full_summary['R_sq'].iloc[0]:.4f}, "
          f"Dev.expl={full_summary['dev_expl'].iloc[0]:.3f}, AIC={full_summary['AIC'].iloc[0]:.1f}")
if delta_summary is not None:
    print(f"DELTA_EXCLUDED:  n={delta_summary['n'].iloc[0]}, R²={delta_summary['R_sq'].iloc[0]:.4f}, "
          f"Dev.expl={delta_summary['dev_expl'].iloc[0]:.3f}, AIC={delta_summary['AIC'].iloc[0]:.1f}")

# ---- 6b. Pacific Coast smooth terms (all SPEI timescales) ----
full_smooth = read_csv_safe(os.path.join(temp_path, "full_smooth_terms.csv"))
delta_smooth = read_csv_safe(os.path.join(temp_path, "delta_smooth_terms.csv"))

def extract_pacific_smooth(df, scenario_name):
    if df is None:
        return None
    # Identify rows with "SPEI_" and "Pacific Coast" in smooth_term
    pac = df[df['smooth_term'].str.contains("SPEI_") & df['smooth_term'].str.contains("Pacific Coast")].copy()
    if len(pac) == 0:
        return None
    # Extract SPEI timescale
    pac['SPEI_timescale'] = pac['smooth_term'].str.extract(r'(SPEI_\d+)')[0]
    pac['scenario'] = scenario_name
    return pac[['scenario', 'SPEI_timescale', 'edf', 'Ref.df', 'F', 'p_value']]

full_pac_smooth = extract_pacific_smooth(full_smooth, "FULL")
delta_pac_smooth = extract_pacific_smooth(delta_smooth, "DELTA_EXCLUDED")

print("\nPACIFIC COAST SPEI SMOOTH STATISTICS")
print("-------------------------------------")
if full_pac_smooth is not None:
    print("FULL:")
    print(full_pac_smooth.to_string(index=False))
if delta_pac_smooth is not None:
    print("\nDELTA_EXCLUDED:")
    print(delta_pac_smooth.to_string(index=False))

# ---- 6c. Pacific SPEI-3 comparison (compact table) ----
def get_spei3_row(df, scenario_name, pac_sites, pac_rows):
    if df is None:
        return None
    row = df[df['SPEI_timescale'] == "SPEI_3"]
    if len(row) == 0:
        return None
    row = row.iloc[0]
    return {
        'scenario': scenario_name,
        'Pacific_sites': pac_sites,
        'Pacific_rows': pac_rows,   # This is original site-month observations, not long-format
        'edf': row['edf'],
        'Ref.df': row['Ref.df'],
        'F': row['F'],
        'p_value': row['p_value'],
        'significant': row['p_value'] < 0.05
    }

# Use the original (non-long) site and observation counts for the final table
full_spei3 = get_spei3_row(full_pac_smooth, "FULL", n_full_pacific_sites, n_full_pacific)
delta_spei3 = get_spei3_row(delta_pac_smooth, "DELTA_EXCLUDED", n_delta_excl_pacific_sites, n_delta_excl_pacific)

print("\nPACIFIC SPEI-3 COMPARISON")
print("--------------------------")
if full_spei3:
    print(f"FULL:            sites={full_spei3['Pacific_sites']}, rows={full_spei3['Pacific_rows']}, "
          f"edf={full_spei3['edf']:.2f}, F={full_spei3['F']:.2f}, p={full_spei3['p_value']:.4f}, "
          f"significant={full_spei3['significant']}")
if delta_spei3:
    print(f"DELTA_EXCLUDED:  sites={delta_spei3['Pacific_sites']}, rows={delta_spei3['Pacific_rows']}, "
          f"edf={delta_spei3['edf']:.2f}, F={delta_spei3['F']:.2f}, p={delta_spei3['p_value']:.4f}, "
          f"significant={delta_spei3['significant']}")

# ---- 6d. Observed Pacific SPEI-3 support ----
print("\nOBSERVED PACIFIC SPEI-3 SUPPORT")
print("-------------------------------")
def print_support(vals, label, n_sites):
    if vals is not None and len(vals) > 0:
        q = np.percentile(vals, [0, 2, 5, 25, 50, 75, 95, 98, 100])
        print(f"{label}:")
        print(f"  sites={n_sites}, observations={len(vals)}")
        print(
            f"  min={q[0]:.2f}, 2nd={q[1]:.2f}, 5th={q[2]:.2f}, "
            f"25th={q[3]:.2f}, med={q[4]:.2f}, 75th={q[5]:.2f}, "
            f"95th={q[6]:.2f}, 98th={q[7]:.2f}, max={q[8]:.2f}"
        )
    else:
        print(f"{label}: no data")

print_support(full_spei3_vals, "FULL", n_full_pacific_sites)
print_support(delta_spei3_vals, "DELTA_EXCLUDED", n_delta_excl_pacific_sites)

# ---- 6e. Thresholds ----
full_thresh = read_csv_safe(os.path.join(temp_path, "full_thresholds.csv"))
delta_thresh = read_csv_safe(os.path.join(temp_path, "delta_thresholds.csv"))

def fmt_threshold(x):
    return "NA" if pd.isna(x) else f"{x:.2f}"

def print_thresholds(df, label):
    if df is not None:
        print(f"{label}:")
        for thr in [5, 10, 20]:
            dec = df[(df['threshold_pct'] == thr) & (df['direction'] == 'decrease')]['SPEI_threshold'].values
            inc = df[(df['threshold_pct'] == thr) & (df['direction'] == 'increase')]['SPEI_threshold'].values
            dec_val = dec[0] if len(dec) > 0 else np.nan
            inc_val = inc[0] if len(inc) > 0 else np.nan
            print(f"  {thr}% decrease: {fmt_threshold(dec_val)}, "
                  f"{thr}% increase: {fmt_threshold(inc_val)}")
    else:
        print(f"{label}: no data")

print("\nDRY-SIDE THRESHOLDS (Pacific Coast, SPEI-3)")
print("-------------------------------------------")
print_thresholds(full_thresh, "FULL")
print_thresholds(delta_thresh, "DELTA_EXCLUDED")

# ---- 6f. Fixed SPEI comparison ----
full_fixed = read_csv_safe(os.path.join(temp_path, "full_fixed_predictions.csv"))
delta_fixed = read_csv_safe(os.path.join(temp_path, "delta_fixed_predictions.csv"))

print("\nFIXED SPEI COMPARISON (Pacific Coast, SPEI-3, Upland, July)")
print("-------------------------------------------------------------")
if full_fixed is not None:
    print("FULL:")
    print(full_fixed[['SPEI_value', 'predicted_WUE_T', 'predicted_lower', 'predicted_upper',
                      'baseline_WUE_T', 'pct_change']].to_string(index=False))
if delta_fixed is not None:
    print("\nDELTA_EXCLUDED:")
    print(delta_fixed[['SPEI_value', 'predicted_WUE_T', 'predicted_lower', 'predicted_upper',
                       'baseline_WUE_T', 'pct_change']].to_string(index=False))

# ---- 6g. Response direction ----
print("\nRESPONSE DIRECTION (Pacific dry-side, SPEI-3)")
print("---------------------------------------------")
# Determine from predictions: use the sign of the fitted curve at dry end (e.g., at SPEI = -2 or similar)
def get_direction(pred_df):
    if pred_df is None:
        return "NA"
    # Use the mean of pct_change for SPEI < -1 (if any)
    dry = pred_df[pred_df['SPEI_value'] < -1]
    if len(dry) == 0:
        return "essentially flat"
    mean_pct = dry['pct_change'].mean()
    if mean_pct > 1:
        return "WUE_T increase"
    elif mean_pct < -1:
        return "WUE_T decrease"
    else:
        return "essentially flat"

full_pred = read_csv_safe(os.path.join(temp_path, "full_predictions.csv"))
delta_pred = read_csv_safe(os.path.join(temp_path, "delta_predictions.csv"))

full_dir = get_direction(full_pred)
delta_dir = get_direction(delta_pred)
print(f"FULL:            {full_dir}")
print(f"DELTA_EXCLUDED:  {delta_dir}")
print(f"Direction unchanged after Delta exclusion: {'YES' if full_dir == delta_dir else 'NO'}")

# ---- 6h. Final compact table ----
print("\nFINAL COMPACT TABLE")
print("--------------------")
# Build a row for each scenario
def build_row(scenario, spei3_info, thresh_df, fixed_df):
    row = {
        'scenario': scenario,
        'Pacific_sites': spei3_info['Pacific_sites'] if spei3_info else None,
        'Pacific_rows': spei3_info['Pacific_rows'] if spei3_info else None,
        'Pacific_SPEI3_edf': spei3_info['edf'] if spei3_info else None,
        'Pacific_SPEI3_F': spei3_info['F'] if spei3_info else None,
        'Pacific_SPEI3_p': spei3_info['p_value'] if spei3_info else None,
        'Pacific_SPEI3_significant': spei3_info['significant'] if spei3_info else None,
    }
    # Add thresholds
    for thr in [5,10,20]:
        for dir in ['decrease', 'increase']:
            val = None
            if thresh_df is not None:
                sub = thresh_df[(thresh_df['threshold_pct']==thr) & (thresh_df['direction']==dir)]
                if len(sub) > 0:
                    val = sub['SPEI_threshold'].iloc[0]
            row[f'dry_{thr}pct_{dir}_threshold'] = val
    # Add pct_change at specific SPEI values
    if fixed_df is not None:
        for spei_val in requested_spei:
            sub = fixed_df[fixed_df['SPEI_value'] == spei_val]
            if len(sub) > 0:
                row[f'pct_change_at_SPEI_minus_{abs(spei_val):.1f}'] = sub['pct_change'].iloc[0]
            else:
                row[f'pct_change_at_SPEI_minus_{abs(spei_val):.1f}'] = np.nan
    else:
        for spei_val in requested_spei:
            row[f'pct_change_at_SPEI_minus_{abs(spei_val):.1f}'] = np.nan
    return row

full_row = build_row("FULL", full_spei3, full_thresh, full_fixed)
delta_row = build_row("DELTA_EXCLUDED", delta_spei3, delta_thresh, delta_fixed)

# Print as DataFrame
compact_df = pd.DataFrame([full_row, delta_row])
# Reorder columns for clarity
base_cols = ['scenario', 'Pacific_sites', 'Pacific_rows', 'Pacific_SPEI3_edf',
             'Pacific_SPEI3_F', 'Pacific_SPEI3_p', 'Pacific_SPEI3_significant']
threshold_cols = [f'dry_{t}pct_{d}_threshold' for t in [5,10,20] for d in ['decrease','increase']]
pct_cols = [f'pct_change_at_SPEI_minus_{abs(x):.1f}' for x in requested_spei]
order = base_cols + threshold_cols + pct_cols
compact_df = compact_df[order]
print(compact_df.to_string(index=False))

# ---- 6i. Final statements ----
print("\nFINAL NUMERICAL CHECKS")
print("-----------------------")
def get_significant(spei3_info):
    return "YES" if (spei3_info and spei3_info['significant']) else "NO"

print(f"Pacific SPEI-3 significant in FULL model:             {get_significant(full_spei3)}")
print(f"Pacific SPEI-3 significant after Delta exclusion:    {get_significant(delta_spei3)}")
print(f"Direction of Pacific dry-side response unchanged after Delta exclusion: {'YES' if full_dir == delta_dir else 'NO'}")

# Detectability of thresholds
def threshold_detectable(thresh_df, thr, direction):
    if thresh_df is None:
        return "NO"
    sub = thresh_df[(thresh_df['threshold_pct']==thr) & (thresh_df['direction']==direction)]
    if len(sub) > 0 and not np.isnan(sub['SPEI_threshold'].iloc[0]):
        return "YES"
    else:
        return "NO"

for thr in [5,10,20]:
    full_dec = threshold_detectable(full_thresh, thr, 'decrease')
    full_inc = threshold_detectable(full_thresh, thr, 'increase')
    delta_dec = threshold_detectable(delta_thresh, thr, 'decrease')
    delta_inc = threshold_detectable(delta_thresh, thr, 'increase')
    print(f"{thr}% dry-side threshold remains detectable after Delta exclusion: "
          f"decrease: {delta_dec} (FULL: {full_dec}), increase: {delta_inc} (FULL: {full_inc})")

# ============================================================================
# CLEANUP (automatically via tempdir)
# ============================================================================
print("\nDiagnostic complete. Temporary files will be deleted.")
temp_dir.cleanup()