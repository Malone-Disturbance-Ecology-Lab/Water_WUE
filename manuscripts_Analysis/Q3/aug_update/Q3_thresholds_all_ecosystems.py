"""
Q3_WUE_T_SPEI_reviewer_all_ecosystem_thresholds.py
==================================================
Lean reviewer-output script.

- Loads the saved coast-threshold GAM model (RDS) and long data (CSV).
- Generates predictions for all coast_region × water_class × month_f × SPEI_timescale × SPEI_value.
- Writes four output files:
  1. Full prediction file (month‑specific)
  2. Month‑averaged prediction curves
  3. Month‑specific threshold markers (5,10,15,20,25,30,35,40,50,75%)
  4. Month‑averaged threshold markers (from averaged curve)

All outputs go to:
    M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_august_update

No old Upland/July files are created.
"""

import os
import subprocess
import warnings
warnings.filterwarnings('ignore')
print("="*60)
print("Q3: REVIEWER – All Ecosystems, All Months (Lean)")
print("="*60)

# ============================================================================
# PATHS
# ============================================================================
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_dir, "Q3_august_update")
temp_dir = os.path.join(output_dir, "temp")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)

# Input files (from original outputs)
data_long_path = os.path.join(base_dir, "Q3_WUE_T_SPEI_sensitivity_outputs",
                              "Q3_WUE_T_SPEI_model_data_long.csv")
model_rds_path = os.path.join(base_dir, "Q3_WUE_T_SPEI_sensitivity_outputs",
                              "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds")

print(f"\nData long path:   {data_long_path}")
print(f"Model RDS path:   {model_rds_path}")
print(f"Output directory: {output_dir}")

# ---- R‑safe path variables ----
data_long_path_r = data_long_path.replace("\\", "/")
model_rds_path_r = model_rds_path.replace("\\", "/")
output_dir_r = output_dir.replace("\\", "/")

# ============================================================================
# CREATE R SCRIPT
# ============================================================================
r_script_content = f'''
# Reviewer script – all ecosystems, all months
# Loads saved model and data, generates predictions & thresholds.

library(mgcv)
library(dplyr)
library(tidyr)

# ---- paths ----
data_long_path <- "{data_long_path_r}"
model_rds_path <- "{model_rds_path_r}"
output_dir <- "{output_dir_r}"

write_table <- function(x, filename) {{
    write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}}

# ---- load data and model ----
cat("Loading data...\\n")
data <- read.csv(data_long_path)
cat("Loading model...\\n")
smooth_coast_model <- readRDS(model_rds_path)

# ---- ensure factors ----
data$site_name <- factor(data$site_name)
data$water_class <- factor(data$water_class, levels = c("Upland", "Freshwater", "Saline"))
data$coast_region <- factor(data$coast_region, levels = c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"))
data$month_f <- as.factor(data$month_f)
data$SPEI_timescale <- factor(data$SPEI_timescale,
                              levels = c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"))
data$spei_coast <- interaction(data$SPEI_timescale, data$coast_region, sep = "__", drop = TRUE)

# ---- prediction grid: all combinations ----
cat("Creating prediction grid...\\n")
spei_lower <- min(data$SPEI_value, na.rm = TRUE)
spei_upper <- max(data$SPEI_value, na.rm = TRUE)

cat("Prediction SPEI range:", round(spei_lower, 3), "to", round(spei_upper, 3), "\n")

spei_sequence <- seq(
    spei_lower,
    spei_upper,
    length.out = 200
)

# Get all levels
water_classes <- levels(data$water_class)
coast_regions <- levels(data$coast_region)
spei_timescales <- levels(data$SPEI_timescale)
months <- levels(data$month_f)   # all months present (already growing‑season filtered)

prediction_grid <- expand.grid(
    SPEI_timescale = spei_timescales,
    coast_region = coast_regions,
    water_class = water_classes,
    month_f = months,
    SPEI_value = spei_sequence,
    site_name = levels(data$site_name)[1]   # dummy site for random effect
)

# ---- explicitly set factor levels ----
prediction_grid$SPEI_timescale <- factor(prediction_grid$SPEI_timescale, levels = levels(data$SPEI_timescale))
prediction_grid$coast_region <- factor(prediction_grid$coast_region, levels = levels(data$coast_region))
prediction_grid$water_class <- factor(prediction_grid$water_class, levels = levels(data$water_class))
prediction_grid$month_f <- factor(prediction_grid$month_f, levels = levels(data$month_f))
prediction_grid$site_name <- factor(prediction_grid$site_name, levels = levels(data$site_name))

prediction_grid$spei_coast <- interaction(prediction_grid$SPEI_timescale,
                                          prediction_grid$coast_region,
                                          sep = "__", drop = TRUE)
prediction_grid$spei_coast <- factor(prediction_grid$spei_coast,
                                     levels = levels(data$spei_coast))

# ---- predict ----
cat("Predicting...\\n")
pred <- predict(smooth_coast_model, newdata = prediction_grid,
                type = "link", se.fit = TRUE, exclude = "s(site_name)")
prediction_grid$predicted_WUE_T <- as.numeric(pred$fit)
prediction_grid$predicted_se <- as.numeric(pred$se.fit)
prediction_grid$predicted_lower <- prediction_grid$predicted_WUE_T - 1.96 * prediction_grid$predicted_se
prediction_grid$predicted_upper <- prediction_grid$predicted_WUE_T + 1.96 * prediction_grid$predicted_se

# ---- near‑normal baseline (per group) ----
cat("Computing baseline...\\n")
baseline <- prediction_grid[
    prediction_grid$SPEI_value >= -1 & prediction_grid$SPEI_value <= 1,
] |>
    group_by(SPEI_timescale, coast_region, water_class, month_f) |>
    summarise(predicted_near_normal_WUE_T = mean(predicted_WUE_T), .groups = "drop")

# ---- percentage changes ----
prediction_impact <- left_join(prediction_grid, baseline,
                               by = c("SPEI_timescale", "coast_region", "water_class", "month_f"))
prediction_impact$predicted_change <- prediction_impact$predicted_WUE_T -
                                       prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_pct_change <- 100 * prediction_impact$predicted_change /
                                           prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_lower_pct <- 100 * (prediction_impact$predicted_lower -
                                                prediction_impact$predicted_near_normal_WUE_T) /
                                          prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_upper_pct <- 100 * (prediction_impact$predicted_upper -
                                                prediction_impact$predicted_near_normal_WUE_T) /
                                          prediction_impact$predicted_near_normal_WUE_T

# ---- 1. Full prediction file (month‑specific) ----
write_table(prediction_impact,
            "Q3_reviewer_coast_threshold_predictions_all_ecosystems_all_months.csv")

# ---- 2. Month‑averaged prediction curves ----
curve_avg <- prediction_impact |>
    group_by(coast_region, water_class, SPEI_timescale, SPEI_value) |>
    summarise(
        mean_pct = mean(predicted_pct_change, na.rm = TRUE),
        lower_pct = mean(predicted_lower_pct, na.rm = TRUE),
        upper_pct = mean(predicted_upper_pct, na.rm = TRUE),
        n_months = n_distinct(month_f),
        .groups = "drop"
    )
write_table(curve_avg,
            "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv")

# ---- Helper: find SPEI threshold ----
find_spei_threshold <- function(group_data, pct_column, threshold_pct,
                                anomaly_side, impact_direction) {{
    pct_values <- group_data[[pct_column]]
    if (anomaly_side == "dry") {{
        side_filter <- group_data$SPEI_value < -1
    }} else {{
        side_filter <- group_data$SPEI_value > 1
    }}
    if (impact_direction == "decrease") {{
        direction_filter <- pct_values <= -threshold_pct
    }} else {{
        direction_filter <- pct_values >= threshold_pct
    }}
    threshold_values <- group_data$SPEI_value[side_filter & direction_filter]
    if (length(threshold_values) == 0) return(NA_real_)
    if (anomaly_side == "dry") {{
        return(max(threshold_values, na.rm = TRUE))
    }} else {{
        return(min(threshold_values, na.rm = TRUE))
    }}
}}

# ---- 3. Month‑specific threshold markers ----
cat("Computing month‑specific thresholds...\\n")
threshold_levels <- c(5, 10, 15, 20, 25, 30, 35, 40, 50, 75)

markers_month <- do.call(
    rbind,
    lapply(
        split(prediction_impact,
              list(prediction_impact$SPEI_timescale,
                   prediction_impact$coast_region,
                   prediction_impact$water_class,
                   prediction_impact$month_f), drop = TRUE),
        function(group) {{
            do.call(
                rbind,
                lapply(threshold_levels, function(thr) {{
                    # Mean curve crossings
                    dry_dec <- group$SPEI_value[group$SPEI_value < -1 & group$predicted_pct_change <= -thr]
                    dry_inc <- group$SPEI_value[group$SPEI_value < -1 & group$predicted_pct_change >= thr]
                    wet_dec <- group$SPEI_value[group$SPEI_value > 1 & group$predicted_pct_change <= -thr]
                    wet_inc <- group$SPEI_value[group$SPEI_value > 1 & group$predicted_pct_change >= thr]

                    data.frame(
                        SPEI_timescale = group$SPEI_timescale[1],
                        coast_region = group$coast_region[1],
                        water_class = group$water_class[1],
                        month_f = group$month_f[1],
                        threshold_pct = thr,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_dec) > 0, max(dry_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_inc) > 0, max(dry_inc, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_dec) > 0, min(wet_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_inc) > 0, min(wet_inc, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(-thr, thr, -thr, thr)
                    )
                }})
            )
        }}
    )
)
# Remove NA thresholds
markers_month <- markers_month[!is.na(markers_month$SPEI_threshold), ]

# Compute lower/upper uncertainty thresholds
markers_month$SPEI_threshold_lower <- NA_real_
markers_month$SPEI_threshold_upper <- NA_real_

for (i in seq_len(nrow(markers_month))) {{
    row <- markers_month[i, ]
    sub <- prediction_impact[
        prediction_impact$SPEI_timescale == row$SPEI_timescale &
        prediction_impact$coast_region == row$coast_region &
        prediction_impact$water_class == row$water_class &
        prediction_impact$month_f == row$month_f,
    ]
    lower <- find_spei_threshold(sub, "predicted_lower_pct",
                                 row$threshold_pct, row$anomaly_side, row$impact_direction)
    upper <- find_spei_threshold(sub, "predicted_upper_pct",
                                 row$threshold_pct, row$anomaly_side, row$impact_direction)
    range_vals <- c(row$SPEI_threshold, lower, upper)
    range_vals <- range_vals[!is.na(range_vals)]
    markers_month$SPEI_threshold_lower[i] <- min(range_vals)
    markers_month$SPEI_threshold_upper[i] <- max(range_vals)
}}

markers_month$threshold_pct <- factor(markers_month$threshold_pct,
                                      levels = threshold_levels,
                                      labels = paste0(threshold_levels, "%"))

write_table(markers_month,
            "Q3_reviewer_threshold_markers_month_specific_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv")

# ---- 4. Month‑averaged threshold markers (from averaged curve) ----
cat("Computing month‑averaged thresholds...\\n")
markers_avg <- do.call(
    rbind,
    lapply(
        split(curve_avg,
              list(curve_avg$SPEI_timescale,
                   curve_avg$coast_region,
                   curve_avg$water_class), drop = TRUE),
        function(group) {{
            do.call(
                rbind,
                lapply(threshold_levels, function(thr) {{
                    dry_dec <- group$SPEI_value[group$SPEI_value < -1 & group$mean_pct <= -thr]
                    dry_inc <- group$SPEI_value[group$SPEI_value < -1 & group$mean_pct >= thr]
                    wet_dec <- group$SPEI_value[group$SPEI_value > 1 & group$mean_pct <= -thr]
                    wet_inc <- group$SPEI_value[group$SPEI_value > 1 & group$mean_pct >= thr]

                    data.frame(
                        SPEI_timescale = group$SPEI_timescale[1],
                        coast_region = group$coast_region[1],
                        water_class = group$water_class[1],
                        threshold_pct = thr,
                        impact_direction = c("decrease", "increase", "decrease", "increase"),
                        anomaly_side = c("dry", "dry", "wet", "wet"),
                        SPEI_threshold = c(
                            ifelse(length(dry_dec) > 0, max(dry_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(dry_inc) > 0, max(dry_inc, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_dec) > 0, min(wet_dec, na.rm = TRUE), NA_real_),
                            ifelse(length(wet_inc) > 0, min(wet_inc, na.rm = TRUE), NA_real_)
                        ),
                        pct_change_threshold = c(-thr, thr, -thr, thr)
                    )
                }})
            )
        }}
    )
)
markers_avg <- markers_avg[!is.na(markers_avg$SPEI_threshold), ]

# Compute lower/upper from averaged curve uncertainty
markers_avg$SPEI_threshold_lower <- NA_real_
markers_avg$SPEI_threshold_upper <- NA_real_

for (i in seq_len(nrow(markers_avg))) {{
    row <- markers_avg[i, ]
    sub <- curve_avg[
        curve_avg$SPEI_timescale == row$SPEI_timescale &
        curve_avg$coast_region == row$coast_region &
        curve_avg$water_class == row$water_class,
    ]
    lower <- find_spei_threshold(sub, "lower_pct",
                                 row$threshold_pct, row$anomaly_side, row$impact_direction)
    upper <- find_spei_threshold(sub, "upper_pct",
                                 row$threshold_pct, row$anomaly_side, row$impact_direction)
    range_vals <- c(row$SPEI_threshold, lower, upper)
    range_vals <- range_vals[!is.na(range_vals)]
    markers_avg$SPEI_threshold_lower[i] <- min(range_vals)
    markers_avg$SPEI_threshold_upper[i] <- max(range_vals)
}}

markers_avg$threshold_pct <- factor(markers_avg$threshold_pct,
                                    levels = threshold_levels,
                                    labels = paste0(threshold_levels, "%"))

write_table(markers_avg,
            "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv")

# ---- console summary ----
cat("\\n=== Summary of month‑averaged threshold crossings ===\\n")
summary_avg <- markers_avg |>
    group_by(coast_region, water_class, SPEI_timescale, threshold_pct) |>
    summarise(n_crossings = n(), .groups = "drop")
print(summary_avg)

cat("\\n=== All reviewer outputs written to:", output_dir, "\\n")
'''

# ============================================================================
# SAVE R SCRIPT AND RUN
# ============================================================================
r_script_file = os.path.join(temp_dir, "run_reviewer.R")
with open(r_script_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print(f"\nR script saved to: {r_script_file}")

# Add R to PATH if needed
r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if r_path not in os.environ['PATH']:
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + r_path

print("\nRunning R script...")
result = subprocess.run(
    ["Rscript", r_script_file],
    capture_output=True,
    text=True
)

if result.returncode != 0:
    print("R script had errors:")
    print(result.stderr)
    raise RuntimeError("R execution failed")
else:
    print("R script completed successfully.")
    if result.stdout:
        print(result.stdout)

# ============================================================================
# VERIFY OUTPUTS (only the four new files)
# ============================================================================
print("\n" + "="*60)
print("Verifying reviewer outputs")
print("="*60)

expected_files = [
    "Q3_reviewer_coast_threshold_predictions_all_ecosystems_all_months.csv",
    "Q3_reviewer_coast_threshold_prediction_curves_month_averaged_all_ecosystems.csv",
    "Q3_reviewer_threshold_markers_month_specific_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv",
    "Q3_reviewer_threshold_markers_month_averaged_5_10_15_20_25_30_35_40_50_75_all_ecosystems.csv",
]

missing = []
for f in expected_files:
    if not os.path.exists(os.path.join(output_dir, f)):
        missing.append(f)

if missing:
    raise FileNotFoundError("Missing output files:\n" + "\n".join(missing))
else:
    print(f"All {len(expected_files)} reviewer outputs successfully created.")
    print("\nFiles saved to:")
    for f in expected_files:
        print(f"  {os.path.join(output_dir, f)}")

print("\n" + "="*60)
print("REVIEWER SCRIPT COMPLETE")
print("="*60)