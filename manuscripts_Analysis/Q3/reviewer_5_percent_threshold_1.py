# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 16:44:30 2026

@author: ammar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5pct_threshold_justification.py

Standalone diagnostic to evaluate and justify the 5% ecological-response threshold
within the graded WUE_T response framework (5–75%).

Uses the saved coast-threshold GAM and proper lpmatrix contrasts to calculate
uncertainty of the difference between predicted WUE_T at a given SPEI and the
near-normal baseline. Excludes the site random effect exactly as in the main Q3 workflow.

All R code is contained in a temporary R script, executed via subprocess.
Outputs are then read back into Python for figure creation.
"""

import os
import sys
import subprocess
import tempfile
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================
BASE_OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\5_percent_threshold_justification"
os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

Q3_OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"
Q3_AUGUST_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_august_update"

MODEL_RDS = os.path.join(Q3_OUTPUT_DIR, "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds")
MODEL_DATA = os.path.join(Q3_OUTPUT_DIR, "Q3_WUE_T_SPEI_model_data_long.csv")
EXISTING_PRED = os.path.join(Q3_AUGUST_DIR, "Q3_reviewer_coast_threshold_predictions_all_ecosystems_all_months.csv")
EXISTING_CURVES_AVG = os.path.join(Q3_AUGUST_DIR, "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv")

# ============================================================================
# CONSTANTS
# ============================================================================
ECOSYSTEMS = ["Upland", "Freshwater", "Saline"]
COASTS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]
MONTHS = list(range(1, 13))
THRESHOLDS = [5, 10, 15, 20, 25, 30, 35, 40, 50, 75]

# R path
R_PATHS = [
    r"C:\Program Files\R\R-4.4.2\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.4.0\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.3\bin\x64\Rscript.exe",
]
RSCRIPT_EXE = None
for p in R_PATHS:
    if os.path.exists(p):
        RSCRIPT_EXE = p
        break
if RSCRIPT_EXE is None:
    RSCRIPT_EXE = shutil.which("Rscript")
if RSCRIPT_EXE is None or not os.path.exists(RSCRIPT_EXE):
    raise RuntimeError("Rscript not found. Please install R and ensure Rscript.exe is in PATH.")

def r_path(path):
    return os.path.abspath(path).replace("\\", "/")

# ============================================================================
# CREATE R SCRIPT
# ============================================================================
print("Creating R script for GAM predictions and lpmatrix contrasts...")

r_script_content = f'''
# 5pct_threshold_justification_diagnostic.R
# Runs in R via subprocess from Python

library(mgcv)
library(dplyr)
library(tidyr)

# ---- Constants (defined inside R) ----
threshold_levels <- c(5, 10, 15, 20, 25, 30, 35, 40, 50, 75)
coast_levels <- c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast")
ecosystem_levels <- c("Upland", "Freshwater", "Saline")
timescale_levels <- c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48")

# ---- Helper functions ----
write_table <- function(df, path) {{
    write.csv(df, path, row.names = FALSE)
}}

# ---- Paths ----
model_rds <- "{r_path(MODEL_RDS)}"
model_data <- "{r_path(MODEL_DATA)}"
existing_pred <- "{r_path(EXISTING_PRED)}"
existing_curves_avg <- "{r_path(EXISTING_CURVES_AVG)}"
output_dir <- "{r_path(BASE_OUTPUT_DIR)}"

# ---- Load data and model ----
cat("Loading model and data...\\n")
data <- read.csv(model_data, stringsAsFactors = FALSE)
model <- readRDS(model_rds)

# ---- Set factors exactly as in original model ----
data$site_name <- factor(data$site_name)
data$water_class <- factor(data$water_class, levels = ecosystem_levels)
data$coast_region <- factor(data$coast_region, levels = coast_levels)
data$month_f <- as.factor(data$month_f)
data$SPEI_timescale <- factor(data$SPEI_timescale, levels = timescale_levels)
data$spei_coast <- interaction(data$SPEI_timescale, data$coast_region, sep = "__", drop = TRUE)
data$spei_coast <- factor(data$spei_coast, levels = levels(data$spei_coast))

# ---- Observational coverage: count unique site-year-month records ----
coverage <- data |>
    distinct(site_name, Year, month, coast_region, water_class) |>
    group_by(coast_region, water_class) |>
    summarise(
        n_sites = n_distinct(site_name),
        n_site_months = n(),
        .groups = "drop"
    )
write_table(coverage, file.path(output_dir, "observational_coverage_by_coast_ecosystem.csv"))

# ---- Verify that we can reproduce existing predictions ----
cat("Verifying predictions against existing Q3 output...\\n")
existing <- read.csv(existing_pred, stringsAsFactors = FALSE)

# Use a subset for verification: Pacific Coast, Upland, SPEI_3, July
verify_subset <- existing[
    existing$SPEI_timescale == "SPEI_3" & 
    existing$coast_region == "Pacific Coast" &
    existing$water_class == "Upland" &
    existing$month_f == "7", 
    c("SPEI_value", "predicted_WUE_T")
]
verify_subset <- verify_subset[order(verify_subset$SPEI_value), ]

# Build prediction grid for the same points
verify_grid <- data.frame(
    SPEI_timescale = factor("SPEI_3", levels = timescale_levels),
    coast_region = factor("Pacific Coast", levels = coast_levels),
    water_class = factor("Upland", levels = ecosystem_levels),
    month_f = factor("7", levels = levels(data$month_f)),
    SPEI_value = verify_subset$SPEI_value,
    site_name = factor(levels(data$site_name)[1], levels = levels(data$site_name))
)
verify_grid$spei_coast <- interaction(verify_grid$SPEI_timescale, verify_grid$coast_region, sep = "__", drop = TRUE)
verify_grid$spei_coast <- factor(verify_grid$spei_coast, levels = levels(data$spei_coast))

# Predict
pred_verify <- predict(model, newdata = verify_grid, type = "link", se.fit = TRUE, exclude = "s(site_name)")
verify_grid$predicted <- as.numeric(pred_verify$fit)

# Compare
max_diff <- max(abs(verify_grid$predicted - verify_subset$predicted_WUE_T), na.rm = TRUE)
if (max_diff > 1e-6) {{
    stop(paste("Predictions do not match existing output. Max diff =", max_diff))
}} else {{
    cat("Verification passed. Max diff =", max_diff, "\\n")
}}

# ---- Define baseline: coast × ecosystem × month × timescale ----
cat("Computing near-normal baselines...\\n")

# Use full observed SPEI range (min to max)
spei_range <- range(data$SPEI_value, na.rm = TRUE)
spei_seq <- seq(spei_range[1], spei_range[2], length.out = 200)

# Create grid for baseline (all combinations)
nn_grid <- expand.grid(
    SPEI_timescale = timescale_levels,
    coast_region = coast_levels,
    water_class = ecosystem_levels,
    month_f = levels(data$month_f),
    SPEI_value = spei_seq,
    site_name = levels(data$site_name)[1],
    stringsAsFactors = FALSE
)
nn_grid$SPEI_timescale <- factor(nn_grid$SPEI_timescale, levels = timescale_levels)
nn_grid$coast_region <- factor(nn_grid$coast_region, levels = coast_levels)
nn_grid$water_class <- factor(nn_grid$water_class, levels = ecosystem_levels)
nn_grid$month_f <- factor(nn_grid$month_f, levels = levels(data$month_f))
nn_grid$site_name <- factor(nn_grid$site_name, levels = levels(data$site_name))
nn_grid$spei_coast <- interaction(nn_grid$SPEI_timescale, nn_grid$coast_region, sep = "__", drop = TRUE)
nn_grid$spei_coast <- factor(nn_grid$spei_coast, levels = levels(data$spei_coast))

# Get lpmatrix for baseline grid
Xp_nn <- predict(model, newdata = nn_grid, type = "lpmatrix", exclude = "s(site_name)")
# Compute predictions
coef_vec <- coef(model)
nn_grid$predicted <- as.numeric(Xp_nn %*% coef_vec)

# Baseline: mean predicted over SPEI in [-1,1] for each group
nn_grid$group_id <- interaction(nn_grid[, c("SPEI_timescale", "coast_region", "water_class", "month_f")], drop = TRUE)
# Compute mean Xp per group for later contrast use
Xp_nn_grouped <- split(as.data.frame(Xp_nn), nn_grid$group_id)

mean_Xp_list <- list()
for (gid in names(Xp_nn_grouped)) {{
    mat <- as.matrix(Xp_nn_grouped[[gid]])
    nn_grid_sub <- nn_grid[nn_grid$group_id == gid, ]
    nn_idx <- nn_grid_sub$SPEI_value >= -1 & nn_grid_sub$SPEI_value <= 1
    if (sum(nn_idx) > 0) {{
        mean_Xp_list[[gid]] <- colMeans(mat[nn_idx, , drop = FALSE], na.rm = TRUE)
    }}
}}

# Baseline table (for later merge)
baseline_table <- nn_grid[
    nn_grid$SPEI_value >= -1 & nn_grid$SPEI_value <= 1, 
] |>
    group_by(SPEI_timescale, coast_region, water_class, month_f) |>
    summarise(baseline_WUE_T = mean(predicted, na.rm = TRUE), .groups = "drop")

# ---- Build full prediction grid for dry-side analysis ----
cat("Building full prediction grid for dry-side analysis...\\n")
full_grid <- expand.grid(
    SPEI_timescale = timescale_levels,
    coast_region = coast_levels,
    water_class = ecosystem_levels,
    month_f = levels(data$month_f),
    SPEI_value = spei_seq,
    site_name = levels(data$site_name)[1],
    stringsAsFactors = FALSE
)
full_grid$SPEI_timescale <- factor(full_grid$SPEI_timescale, levels = timescale_levels)
full_grid$coast_region <- factor(full_grid$coast_region, levels = coast_levels)
full_grid$water_class <- factor(full_grid$water_class, levels = ecosystem_levels)
full_grid$month_f <- factor(full_grid$month_f, levels = levels(data$month_f))
full_grid$site_name <- factor(full_grid$site_name, levels = levels(data$site_name))
full_grid$spei_coast <- interaction(full_grid$SPEI_timescale, full_grid$coast_region, sep = "__", drop = TRUE)
full_grid$spei_coast <- factor(full_grid$spei_coast, levels = levels(data$spei_coast))

# Get lpmatrix for full grid
Xp_full <- predict(model, newdata = full_grid, type = "lpmatrix", exclude = "s(site_name)")

# ---- Compute contrasts using lpmatrix ----
cat("Computing contrasts from near-normal baseline...\\n")

full_grid$group_id <- interaction(full_grid[, c("SPEI_timescale", "coast_region", "water_class", "month_f")], drop = TRUE)
unique_groups <- unique(full_grid$group_id)

# Pre-allocate results
contrast <- rep(NA_real_, nrow(full_grid))
contrast_se <- rep(NA_real_, nrow(full_grid))
contrast_lower <- rep(NA_real_, nrow(full_grid))
contrast_upper <- rep(NA_real_, nrow(full_grid))

Vp <- model$Vp

for (gid in unique_groups) {{
    idx <- which(full_grid$group_id == gid)
    if (length(idx) == 0) next
    if (!(gid %in% names(mean_Xp_list))) {{
        next
    }}
    mean_Xp <- mean_Xp_list[[gid]]
    Xp_rows <- Xp_full[idx, , drop = FALSE]
    C <- Xp_rows - matrix(mean_Xp, nrow = length(idx), ncol = length(mean_Xp), byrow = TRUE)
    contrast[idx] <- C %*% coef_vec
    for (i in 1:length(idx)) {{
        ci <- C[i, , drop = FALSE]  # 1 x p matrix
        # Correct variance: ci (1xp) %*% Vp (pxp) %*% t(ci) (px1) -> 1x1
        var_i <- ci %*% Vp %*% t(ci)
        se_i <- sqrt(max(var_i[1, 1], 0))
        contrast_se[idx[i]] <- se_i
        contrast_lower[idx[i]] <- contrast[idx[i]] - 1.96 * se_i
        contrast_upper[idx[i]] <- contrast[idx[i]] + 1.96 * se_i
    }}
}}

full_grid$contrast <- contrast
full_grid$contrast_se <- contrast_se
full_grid$contrast_lower <- contrast_lower
full_grid$contrast_upper <- contrast_upper

# Merge baseline table
full_grid <- full_grid |>
    left_join(baseline_table, by = c("SPEI_timescale", "coast_region", "water_class", "month_f"))

# Calculate month-specific percent changes (manuscript definition)
full_grid$pct_change <- 100 * full_grid$contrast / full_grid$baseline_WUE_T
full_grid$pct_change_lower <- 100 * full_grid$contrast_lower / full_grid$baseline_WUE_T
full_grid$pct_change_upper <- 100 * full_grid$contrast_upper / full_grid$baseline_WUE_T

# ---- Month-averaged effect-size curve: mean of month-specific pct_change ----
cat("Computing month-averaged effect-size curve (mean of month-specific percent changes)...\\n")
avg_pct <- full_grid |>
    group_by(SPEI_timescale, coast_region, water_class, SPEI_value) |>
    summarise(
        mean_pct = mean(pct_change, na.rm = TRUE),
        .groups = "drop"
    )

# ---- Verification against existing month-averaged curves ----
cat("Verifying month-averaged percent-change curve against existing Q3 file...\\n")
existing_avg <- read.csv(existing_curves_avg, stringsAsFactors = FALSE)

# Compare a subset: Pacific Coast, Upland, SPEI_3
verify_avg <- existing_avg[
    existing_avg$SPEI_timescale == "SPEI_3" &
    existing_avg$coast_region == "Pacific Coast" &
    existing_avg$water_class == "Upland",
    c("SPEI_value", "mean_pct")
]
verify_avg <- verify_avg[order(verify_avg$SPEI_value), ]

diag_avg <- avg_pct[
    avg_pct$SPEI_timescale == "SPEI_3" &
    avg_pct$coast_region == "Pacific Coast" &
    avg_pct$water_class == "Upland",
    c("SPEI_value", "mean_pct")
]
diag_avg <- diag_avg[order(diag_avg$SPEI_value), ]

# Interpolate to common SPEI grid for comparison
common_spei <- intersect(verify_avg$SPEI_value, diag_avg$SPEI_value)
if (length(common_spei) == 0) {{
    stop("Could not directly verify the month-averaged curve against the existing Q3 output because the SPEI grids do not match.")
}}
verify_sub <- verify_avg[verify_avg$SPEI_value %in% common_spei, ]
diag_sub <- diag_avg[diag_avg$SPEI_value %in% common_spei, ]
max_avg_diff <- max(abs(verify_sub$mean_pct - diag_sub$mean_pct), na.rm = TRUE)
if (max_avg_diff > 1e-6) {{
    stop(paste("Month-averaged curves do not match. Max diff =", max_avg_diff))
}} else {{
    cat("Month-averaged curve verification passed. Max diff =", max_avg_diff, "\\n")
}}

# ---- Functions to find threshold crossings ----
find_crossing <- function(group_data, direction, threshold_pct) {{
    sub <- group_data[group_data$SPEI < -1, ]
    if (nrow(sub) == 0) return(NA_real_)
    if (direction == "increase") {{
        crossed <- sub[sub$mean_pct >= threshold_pct, ]
    }} else {{
        crossed <- sub[sub$mean_pct <= -threshold_pct, ]
    }}
    if (nrow(crossed) == 0) return(NA_real_)
    return(max(crossed$SPEI, na.rm = TRUE))
}}

# ---- Month-averaged threshold crossings ----
cat("Computing month-averaged threshold crossings...\\n")

# Compute month-averaged contrast and SE using lpmatrix (proper averaging of contrast vectors)
cat("Computing month-averaged contrast and SE using lpmatrix...\\n")

avg_contrast_df <- data.frame()
for (coast in coast_levels) {{
    for (eco in ecosystem_levels) {{
        for (ts in timescale_levels) {{
            rows <- full_grid[full_grid$coast_region == coast & full_grid$water_class == eco & full_grid$SPEI_timescale == ts, ]
            if (nrow(rows) == 0) next
            for (sp in spei_seq) {{
                sub_rows <- rows[rows$SPEI_value == sp, ]
                if (nrow(sub_rows) == 0) next
                gids <- sub_rows$group_id
                Xp_avg <- rep(0, ncol(Xp_full))
                mean_Xp_avg <- rep(0, ncol(Xp_full))
                valid_count <- 0
                for (gid in gids) {{
                    idx <- which(full_grid$group_id == gid & full_grid$SPEI_value == sp)
                    if (length(idx) == 0) next
                    if (gid %in% names(mean_Xp_list)) {{
                        Xp_row <- Xp_full[idx, , drop = FALSE]
                        Xp_avg <- Xp_avg + Xp_row[1, ]
                        mean_Xp_avg <- mean_Xp_avg + mean_Xp_list[[gid]]
                        valid_count <- valid_count + 1
                    }}
                }}
                if (valid_count == 0) next
                Xp_avg <- Xp_avg / valid_count
                mean_Xp_avg <- mean_Xp_avg / valid_count
                C_avg <- as.numeric(Xp_avg - mean_Xp_avg)  # now a vector of length p
                contrast_avg <- sum(C_avg * coef_vec)
                # Variance: t(C_avg) %*% Vp %*% C_avg
                var_avg <- t(C_avg) %*% Vp %*% C_avg
                se_avg <- sqrt(max(var_avg[1, 1], 0))
                lower_avg <- contrast_avg - 1.96 * se_avg
                upper_avg <- contrast_avg + 1.96 * se_avg

                # Baseline for this group (average across months)
                baseline_avg_val <- baseline_table |>
                    filter(SPEI_timescale == ts, coast_region == coast, water_class == eco) |>
                    summarise(baseline_WUE_T = mean(baseline_WUE_T, na.rm = TRUE)) |>
                    pull(baseline_WUE_T)

                avg_contrast_df <- rbind(avg_contrast_df, data.frame(
                    coast = coast,
                    ecosystem = eco,
                    timescale = ts,
                    SPEI = sp,
                    contrast = contrast_avg,
                    se = se_avg,
                    lower = lower_avg,
                    upper = upper_avg,
                    baseline = baseline_avg_val,
                    pct_change = 100 * contrast_avg / baseline_avg_val,
                    pct_lower = 100 * lower_avg / baseline_avg_val,
                    pct_upper = 100 * upper_avg / baseline_avg_val
                ))
            }}
        }}
    }}
}}

# Now avg_pct contains the month-averaged effect-size curve (mean of month-specific pct changes)
# and avg_contrast_df contains the proper averaged contrast and its SE.

avg_5pct_results <- list()
avg_threshold_results <- list()

for (coast in coast_levels) {{
    for (eco in ecosystem_levels) {{
        for (ts in timescale_levels) {{
            # Effect-size curve (avg_pct)
            curve_sub <- avg_pct[
                avg_pct$coast_region == coast &
                avg_pct$water_class == eco &
                avg_pct$SPEI_timescale == ts,
                c("SPEI_value", "mean_pct")
            ]
            colnames(curve_sub) <- c("SPEI", "mean_pct")

            # Statistical support curve (avg_contrast_df)
            stat_sub <- avg_contrast_df[
                avg_contrast_df$coast == coast &
                avg_contrast_df$ecosystem == eco &
                avg_contrast_df$timescale == ts,
                c("SPEI", "contrast", "lower", "upper", "pct_change", "pct_lower", "pct_upper")
            ]

            # For each direction
            for (dir in c("increase", "decrease")) {{
                # Effect-size crossings for all thresholds (5-75%)
                for (thr in threshold_levels) {{
                    spei_cross <- find_crossing(curve_sub, dir, thr)
                    if (!is.na(spei_cross)) {{
                        row <- curve_sub[which.min(abs(curve_sub$SPEI - spei_cross)), ]
                        # Check statistical support at this crossing from stat_sub
                        stat_row <- stat_sub[which.min(abs(stat_sub$SPEI - spei_cross)), ]
                        if (nrow(stat_row) == 0) {{
                            stat_support <- NA
                        }} else {{
                            if (dir == "increase") {{
                                stat_support <- stat_row$pct_lower > 0
                            }} else {{
                                stat_support <- stat_row$pct_upper < 0
                            }}
                        }}
                        avg_threshold_results <- c(avg_threshold_results, list(list(
                            coast = coast,
                            ecosystem = eco,
                            timescale = ts,
                            direction = dir,
                            threshold = thr,
                            spei_crossing = spei_cross,
                            pct_change_at_crossing = row$mean_pct,
                            stat_support_reached = stat_support
                        )))
                    }}
                }}

                # 5% specific with statistical support threshold
                spei_5 <- find_crossing(curve_sub, dir, 5)
                if (!is.na(spei_5)) {{
                    # Get CI at the 5% crossing directly from stat_sub
                    stat_row_5 <- stat_sub[which.min(abs(stat_sub$SPEI - spei_5)), ]
                    if (nrow(stat_row_5) == 0) {{
                        ci_supports <- FALSE
                    }} else {{
                        if (dir == "increase") {{
                            ci_supports <- stat_row_5$pct_lower[1] > 0
                        }} else {{
                            ci_supports <- stat_row_5$pct_upper[1] < 0
                        }}
                    }}

                    # Statistical support threshold (first SPEI where CI excludes zero)
                    sub_stat <- stat_sub[stat_sub$SPEI < -1, ]
                    if (dir == "increase") {{
                        sig <- sub_stat[sub_stat$contrast > 0 & sub_stat$lower > 0, ]
                    }} else {{
                        sig <- sub_stat[sub_stat$contrast < 0 & sub_stat$upper < 0, ]
                    }}
                    spei_stat <- if (nrow(sig) > 0) max(sig$SPEI) else NA_real_

                    diff_sp <- spei_5 - spei_stat
                    if (is.na(spei_stat)) {{
                        status <- "no_stat_support"
                    }} else if (ci_supports) {{
                        status <- "supported_at_5pct"
                    }} else if (spei_stat > spei_5) {{
                        status <- "support_before_5pct"
                    }} else {{
                        status <- "support_after_5pct"
                    }}
                    avg_5pct_results <- c(avg_5pct_results, list(list(
                        coast = coast,
                        ecosystem = eco,
                        timescale = ts,
                        direction = dir,
                        spei_5 = spei_5,
                        spei_stat = spei_stat,
                        diff = diff_sp,
                        ci_supports = ci_supports,
                        status = status,
                        abs_diff = abs(diff_sp),
                        agreement_0.25 = abs(diff_sp) <= 0.25
                    )))
                }}
            }}
        }}
    }}
}}

avg_5pct_df <- do.call(rbind, lapply(avg_5pct_results, as.data.frame))
avg_threshold_df <- do.call(rbind, lapply(avg_threshold_results, as.data.frame))

# ---- Save outputs ----
cat("Saving outputs...\\n")

# Month-averaged 5% summary
avg_5pct_out <- avg_5pct_df[, c("coast", "ecosystem", "timescale", "direction", "spei_5", "spei_stat", "diff", "ci_supports", "status", "agreement_0.25")]
write_table(avg_5pct_out, file.path(output_dir, "5pct_threshold_month_averaged_summary.csv"))

# Effect-size sensitivity (5-75%)
higher_summary <- avg_threshold_df |>
    group_by(coast, ecosystem, timescale, direction, threshold) |>
    summarise(
        n_crossings = n(),
        mean_spei = mean(spei_crossing, na.rm = TRUE),
        median_spei = median(spei_crossing, na.rm = TRUE),
        pct_stat_supported = mean(stat_support_reached, na.rm = TRUE) * 100,
        .groups = "drop"
    )
write_table(higher_summary, file.path(output_dir, "effect_size_threshold_sensitivity_5_to_75pct.csv"))

# Supplementary table (concise 5% only)
supp_table <- avg_5pct_out[
    !is.na(avg_5pct_out$spei_5), 
    c("coast", "ecosystem", "timescale", "direction", "spei_5", "spei_stat", "diff", "ci_supports")
]
colnames(supp_table) <- c("Coast", "Ecosystem", "SPEI timescale", "Response direction",
                          "SPEI at 5% response", "SPEI at statistical-support threshold",
                          "Difference in SPEI", "95% CI supports departure from near-normal at 5% crossing")
write_table(supp_table, file.path(output_dir, "Table_Sx_5pct_threshold_justification.csv"))

# ---- Console output ----
# Use paste(rep("=", 70), collapse = "") for separator lines
sep_line <- paste(rep("=", 70), collapse = "")
cat("\n", sep_line, "\n", sep="")
cat("DIAGNOSTIC SUMMARY\n")
cat(sep_line, "\n", sep="")

n_5pct_avg <- nrow(avg_5pct_out[!is.na(avg_5pct_out$spei_5), ])
n_supp_5pct <- sum(avg_5pct_out$ci_supports, na.rm = TRUE)
n_support_before <- sum(avg_5pct_out$status == "support_before_5pct", na.rm = TRUE)
n_support_after <- sum(avg_5pct_out$status == "support_after_5pct", na.rm = TRUE)
n_no_stat <- sum(avg_5pct_out$status == "no_stat_support", na.rm = TRUE)
median_diff <- median(abs(avg_5pct_out$diff), na.rm = TRUE)

cat("\n5% CRITERION (month-averaged):\n")
cat("  Total dry-side 5% crossings detected:", n_5pct_avg, "\n")
cat("  Already statistically supported at 5% crossing:", n_supp_5pct, 
    sprintf("(%.1f%%)", 100 * n_supp_5pct / max(1, n_5pct_avg)), "\n")
cat("  Statistical support occurs before 5% crossing:", n_support_before, "\n")
cat("  Statistical support occurs after 5% crossing:", n_support_after, "\n")
cat("  No statistical support threshold found:", n_no_stat, "\n")
cat("  Median absolute SPEI difference between thresholds:", 
    ifelse(is.finite(median_diff), sprintf("%.3f", median_diff), "NA"), "\n")

cat("\nHIGHER THRESHOLDS (5-75%, month-averaged):\n")
for (thr in threshold_levels) {{
    sub <- higher_summary[higher_summary$threshold == thr, ]
    n_cross <- nrow(sub)
    n_supp <- sum(sub$pct_stat_supported > 0, na.rm = TRUE)
    cat(sprintf("  %d%%: %d crossings, %.1f%% groups with any support\n", thr, n_cross, 
                100 * n_supp / max(1, n_cross)))
}}

cat("\nBY COAST (month-averaged 5%):\n")
for (cst in coast_levels) {{
    sub <- avg_5pct_out[avg_5pct_out$coast == cst & !is.na(avg_5pct_out$spei_5), ]
    n_cross <- nrow(sub)
    n_supp <- sum(sub$ci_supports, na.rm = TRUE)
    cat(sprintf("  %s: %d crossings, %.1f%% supported\n", cst, n_cross,
                100 * n_supp / max(1, n_cross)))
}}

cat("\nBY ECOSYSTEM (month-averaged 5%):\n")
for (eco in ecosystem_levels) {{
    sub <- avg_5pct_out[avg_5pct_out$ecosystem == eco & !is.na(avg_5pct_out$spei_5), ]
    n_cross <- nrow(sub)
    n_supp <- sum(sub$ci_supports, na.rm = TRUE)
    cat(sprintf("  %s: %d crossings, %.1f%% supported\n", eco, n_cross,
                100 * n_supp / max(1, n_cross)))
}}

cat("\nKEY MANUSCRIPT COMBINATIONS (month-averaged):\n")
key_combos <- list(
    c("Pacific Coast", "Upland", "SPEI_3"),
    c("Atlantic Coast", "Upland", "SPEI_3"),
    c("AK Coast", "Upland", "SPEI_48"),
    c("Gulf Coast", "Upland", "SPEI_3")
)
for (combo in key_combos) {{
    cst <- combo[1]; eco <- combo[2]; ts <- combo[3]
    sub <- avg_5pct_out[avg_5pct_out$coast == cst & 
                        avg_5pct_out$ecosystem == eco & 
                        avg_5pct_out$timescale == ts, ]
    if (nrow(sub) > 0) {{
        cat(sprintf("\n  %s | %s | %s:\n", cst, eco, ts))
        for (dir in c("increase", "decrease")) {{
            sub2 <- sub[sub$direction == dir, ]
            if (nrow(sub2) > 0) {{
                sp5 <- sub2$spei_5[1]
                sp_stat <- sub2$spei_stat[1]
                supp <- sub2$ci_supports[1]
                status <- sub2$status[1]
                cat(sprintf("    %s: 5%% at SPEI=%.3f, stat support at SPEI=%.3f, supported=%s, status=%s\n",
                           dir, sp5, sp_stat, supp, status))
            }}
        }}
    }}
}}

cat("\nDiagnostic complete. Outputs saved to:", output_dir, "\n")
'''

# ============================================================================
# WRITE R SCRIPT TO TEMP AND RUN
# ============================================================================
temp_r_file = os.path.join(tempfile.gettempdir(), "5pct_threshold_diagnostic.R")
with open(temp_r_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print("Running R script...")
result = subprocess.run([RSCRIPT_EXE, temp_r_file], capture_output=True, text=True)

if result.returncode != 0:
    print("R script failed with error:")
    print(result.stderr)
    sys.exit(1)
else:
    print(result.stdout)
    if result.stderr:
        print("R warnings:")
        print(result.stderr)

# Clean up temp file
try:
    os.remove(temp_r_file)
except:
    pass

# ============================================================================
# CREATE FIGURE IN PYTHON (with fixed indexing)
# ============================================================================
print("\nCreating figure...")

# Read the necessary CSVs
avg_5pct = pd.read_csv(os.path.join(BASE_OUTPUT_DIR, "5pct_threshold_month_averaged_summary.csv"))
higher_summary = pd.read_csv(os.path.join(BASE_OUTPUT_DIR, "effect_size_threshold_sensitivity_5_to_75pct.csv"))
coverage = pd.read_csv(os.path.join(BASE_OUTPUT_DIR, "observational_coverage_by_coast_ecosystem.csv"))

# Figure: Panel a - 5% SPEI vs statistical support SPEI
# Panel b - threshold sensitivity (SPEI crossing vs response magnitude) with direction preserved

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel a
ax1 = axes[0]
# Proper pandas indexing: no trailing comma
mask = (~avg_5pct['spei_5'].isna()) & (~avg_5pct['spei_stat'].isna())
plot_df = avg_5pct[mask].copy()

if len(plot_df) > 0:
    colors = {'Upland': '#800080', 'Freshwater': '#0000FF', 'Saline': '#FFA500'}
    for eco, color in colors.items():
        sub = plot_df[plot_df['ecosystem'] == eco]
        if len(sub) > 0:
            ax1.scatter(sub['spei_5'], sub['spei_stat'], 
                       color=color, alpha=0.6, s=40, label=eco)
    # 1:1 line
    min_val = min(plot_df['spei_5'].min(), plot_df['spei_stat'].min())
    max_val = max(plot_df['spei_5'].max(), plot_df['spei_stat'].max())
    lims = [min_val - 0.5, max_val + 0.5]
    ax1.plot([lims[0], lims[1]], [lims[0], lims[1]], 'k--', alpha=0.5, linewidth=1)
    ax1.set_xlabel('SPEI at 5% response threshold', fontsize=12, fontweight='bold')
    ax1.set_ylabel('SPEI at statistical-support threshold', fontsize=12, fontweight='bold')
    ax1.set_title('(a) 5% Threshold vs Statistical Support', fontsize=13, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
else:
    ax1.text(0.5, 0.5, 'No valid data for panel (a)', ha='center', va='center', transform=ax1.transAxes)

# Panel b: threshold sensitivity, preserving direction
ax2 = axes[1]
for dir in ['increase', 'decrease']:
    # No trailing comma
    sub = higher_summary[higher_summary['direction'] == dir]
    if len(sub) == 0: continue
    # Aggregate across ecosystems for each coast and threshold
    for cst in COASTS:
        # No trailing comma
        sub2 = sub[sub['coast'] == cst]
        if len(sub2) == 0: continue
        med_by_thresh = sub2.groupby('threshold')['median_spei'].mean()
        linestyle = '-' if dir == 'increase' else '--'
        color = {'Pacific Coast':'#0072B2', 'Gulf Coast':'#B3B300', 'Atlantic Coast':'#CC79A7', 'AK Coast':'#009E73'}.get(cst, 'gray')
        ax2.plot(med_by_thresh.index, med_by_thresh.values, 
                marker='o', linewidth=1.5, linestyle=linestyle, color=color, 
                label=f"{cst} {dir}")

ax2.axhline(y=-1, color='gray', linestyle='--', alpha=0.7, linewidth=1)
ax2.set_xlabel('Response magnitude threshold (%)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Median SPEI crossing', fontsize=12, fontweight='bold')
ax2.set_title('(b) Threshold sensitivity (5–75%)', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=9, ncol=2)
ax2.grid(alpha=0.3)
ax2.invert_xaxis()

plt.tight_layout()
fig_path = os.path.join(BASE_OUTPUT_DIR, "Figure_Sx_5pct_threshold_justification.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Figure saved to: {fig_path}")
print("\nAll outputs saved to:", BASE_OUTPUT_DIR)
print("="*70)