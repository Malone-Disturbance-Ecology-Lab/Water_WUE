# Q3 WUE_T sensitivity to SPEI gradients
#
# Goal:
#   Determine how sensitive transpiration-based water-use efficiency (WUE_T)
#   is to changes in SPEI across short- and long-term moisture timescales.
#
# Response:
#   WUE_T = WUE_tra
#
# SPEI timescales:
#   SPEI_1, SPEI_3, SPEI_6, SPEI_12, SPEI_24, SPEI_36, SPEI_48
#
# Main interpretation:
#   A negative SPEI slope means WUE_T declines as moisture availability
#   increases, or rises as drought intensifies. A positive slope means WUE_T
#   increases as SPEI becomes wetter.
#
# Inputs:
#   Water_WUE/data/WUE_CUE_monthly_merged_indices_clean.csv
#
# Outputs:
#   Water_WUE/Malone_Workflow/results/Q3_WUE_T_SPEI_sensitivity/

suppressPackageStartupMessages({
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required. Install it with install.packages('mgcv').")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required. Install it with install.packages('scales').")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork').")
  }
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }

  candidates <- c(
    "Q3.WUE_Sensitivity_SPEI.R",
    file.path("Malone_Workflow", "Q3.WUE_Sensitivity_SPEI.R"),
    file.path("Water_WUE", "Malone_Workflow", "Q3.WUE_Sensitivity_SPEI.R")
  )
  matches <- candidates[file.exists(candidates)]
  if (length(matches) > 0) {
    return(normalizePath(matches[1], mustWork = TRUE))
  }

  stop("Could not locate Q3.WUE_Sensitivity_SPEI.R. Run from Water_WUE or Water_WUE/Malone_Workflow, or call with Rscript path/to/Q3.WUE_Sensitivity_SPEI.R.")
}

write_table <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}

capture_to_file <- function(expr, filename) {
  sink(file.path(output_dir, filename))
  on.exit(sink(), add = TRUE)
  force(expr)
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

theme_wue <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

classify_coast_region <- function(lat, long) {
  ifelse(
    lat > 50, "AK Coast",
    ifelse(
      long > -100, "Atlantic Coast",
      ifelse(long < -120, "Pacific Coast", "Gulf Coast")
    )
  )
}

script_dir <- dirname(get_script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

input_file <- file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv")
output_dir <- file.path(script_dir, "results", "Q3_WUE_T_SPEI_sensitivity")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
monthly[] <- lapply(monthly, function(x) {
  if (is.character(x)) trimws(x) else x
})

spei_cols <- c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48")
required_cols <- c("site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra", spei_cols)

missing_cols <- setdiff(required_cols, names(monthly))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ecosystem_classes <- c("Upland", "Freshwater", "Saline")

model_base <- monthly[complete.cases(monthly[, required_cols]), ]
model_base <- model_base[
    model_base$site_name != "" &
    model_base$water_class != "" &
    model_base$water_class %in% ecosystem_classes &
    is.finite(model_base$lat) &
    is.finite(model_base$long) &
    is.finite(model_base$WUE_tra),
]

coast_region_levels <- c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast")
model_base$coast_region <- factor(
  classify_coast_region(model_base$lat, model_base$long),
  levels = coast_region_levels
)
model_base$site_name <- factor(model_base$site_name)
model_base$water_class <- factor(model_base$water_class, levels = ecosystem_classes)
model_base$month_f <- factor(model_base$month, levels = 1:12)
model_base$WUE_T <- model_base$WUE_tra

spei_long <- reshape(
  model_base,
  varying = spei_cols,
  v.names = "SPEI_value",
  timevar = "SPEI_timescale",
  times = spei_cols,
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)
spei_long <- spei_long[is.finite(spei_long$SPEI_value) & is.finite(spei_long$WUE_T), ]
spei_long$SPEI_timescale <- factor(spei_long$SPEI_timescale, levels = spei_cols)
spei_long$site_name <- factor(spei_long$site_name)
spei_long$water_class <- factor(spei_long$water_class, levels = ecosystem_classes)
spei_long$coast_region <- factor(spei_long$coast_region, levels = coast_region_levels)
spei_long$month_f <- factor(spei_long$month, levels = 1:12)
spei_long$SPEI_z <- as.numeric(scale(spei_long$SPEI_value))
spei_long$spei_ecosystem <- interaction(
  spei_long$SPEI_timescale,
  spei_long$water_class,
  sep = "__",
  drop = TRUE
)
spei_long$spei_coast <- interaction(
  spei_long$SPEI_timescale,
  spei_long$coast_region,
  sep = "__",
  drop = TRUE
)

write_table(model_base, "Q3_WUE_T_SPEI_model_base_wide.csv")
write_table(spei_long, "Q3_WUE_T_SPEI_model_data_long.csv")

coverage <- spei_long |>
  dplyr::group_by(SPEI_timescale, water_class) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_SPEI = mean(SPEI_value),
    min_SPEI = min(SPEI_value),
    max_SPEI = max(SPEI_value),
    mean_WUE_T = mean(WUE_T),
    .groups = "drop"
  )
write_table(coverage, "Q3_WUE_T_SPEI_coverage.csv")

site_slopes <- spei_long |>
  dplyr::group_by(site_name, water_class, coast_region, SPEI_timescale) |>
  dplyr::filter(dplyr::n_distinct(SPEI_value) >= 4) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    SPEI_slope = stats::coef(stats::lm(WUE_T ~ SPEI_value))[2],
    SPEI_cor = stats::cor(WUE_T, SPEI_value),
    mean_WUE_T = mean(WUE_T),
    .groups = "drop"
  )
write_table(site_slopes, "Q3_WUE_T_SPEI_site_level_slopes.csv")

site_near_normal_baselines <- spei_long |>
  dplyr::filter(SPEI_value >= -1, SPEI_value <= 1) |>
  dplyr::group_by(site_name, water_class, coast_region, SPEI_timescale, month_f) |>
  dplyr::summarise(
    baseline_WUE_T = mean(WUE_T),
    baseline_sd_WUE_T = stats::sd(WUE_T),
    baseline_n = dplyr::n(),
    .groups = "drop"
  )

site_response_data <- dplyr::left_join(
  spei_long,
  site_near_normal_baselines,
  by = c("site_name", "water_class", "coast_region", "SPEI_timescale", "month_f")
)
site_response_data <- site_response_data[
  is.finite(site_response_data$baseline_WUE_T) &
    site_response_data$baseline_n >= 3,
]
site_response_data$WUE_T_change <- site_response_data$WUE_T - site_response_data$baseline_WUE_T
site_response_data$WUE_T_pct_change <- 100 * site_response_data$WUE_T_change / site_response_data$baseline_WUE_T
site_response_data$impact_class_10pct <- factor(
  ifelse(
    site_response_data$WUE_T_pct_change <= -10, "WUE decrease",
    ifelse(site_response_data$WUE_T_pct_change >= 10, "WUE increase", "No meaningful change")
  ),
  levels = c("WUE decrease", "No meaningful change", "WUE increase")
)
write_table(site_near_normal_baselines, "Q3_WUE_T_SPEI_site_near_normal_baselines.csv")
write_table(site_response_data, "Q3_WUE_T_SPEI_site_response_from_near_normal.csv")

site_sensitivity_rank <- site_response_data |>
  dplyr::group_by(site_name, water_class, coast_region, SPEI_timescale) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    baseline_WUE_T = dplyr::first(baseline_WUE_T),
    mean_abs_pct_change = mean(abs(WUE_T_pct_change)),
    max_abs_pct_change = max(abs(WUE_T_pct_change)),
    pct_months_decrease_10 = mean(impact_class_10pct == "WUE decrease") * 100,
    pct_months_increase_10 = mean(impact_class_10pct == "WUE increase") * 100,
    pct_months_impacted_10 = mean(impact_class_10pct != "No meaningful change") * 100,
    dry_threshold_observed = ifelse(
      any(SPEI_value < -1 & WUE_T_pct_change <= -10),
      max(SPEI_value[SPEI_value < -1 & WUE_T_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    wet_threshold_observed = ifelse(
      any(SPEI_value > 1 & WUE_T_pct_change <= -10),
      min(SPEI_value[SPEI_value > 1 & WUE_T_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(SPEI_timescale, dplyr::desc(mean_abs_pct_change))
write_table(site_sensitivity_rank, "Q3_WUE_T_SPEI_site_sensitivity_rank.csv")

ecosystem_slope_summary <- site_slopes |>
  dplyr::group_by(SPEI_timescale, water_class) |>
  dplyr::summarise(
    n_sites = dplyr::n(),
    mean_slope = mean(SPEI_slope),
    median_slope = median(SPEI_slope),
    sd_slope = sd(SPEI_slope),
    se_slope = sd_slope / sqrt(n_sites),
    ci95_slope = stats::qt(0.975, df = n_sites - 1) * se_slope,
    pct_negative = mean(SPEI_slope < 0) * 100,
    .groups = "drop"
  )
write_table(ecosystem_slope_summary, "Q3_WUE_T_SPEI_ecosystem_slope_summary.csv")

coast_slope_summary <- site_slopes |>
  dplyr::group_by(SPEI_timescale, coast_region) |>
  dplyr::summarise(
    n_sites = dplyr::n(),
    mean_slope = mean(SPEI_slope),
    median_slope = median(SPEI_slope),
    sd_slope = sd(SPEI_slope),
    se_slope = sd_slope / sqrt(n_sites),
    ci95_slope = stats::qt(0.975, df = n_sites - 1) * se_slope,
    pct_negative = mean(SPEI_slope < 0) * 100,
    .groups = "drop"
  )
write_table(coast_slope_summary, "Q3_WUE_T_SPEI_coast_region_slope_summary.csv")

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

# The GAM is the primary Q3 model because Q3 asks where SPEI gradients imply
# ecologically meaningful WUE_T change, not whether linear class contrasts are
# significant. Nonlinear SPEI responses are estimated separately for each
# SPEI-timescale-by-ecosystem combination.
smooth_model <- mgcv::gam(
  WUE_T ~
    SPEI_timescale * water_class +
    month_f +
    s(SPEI_value, by = spei_ecosystem, k = 6) +
    s(site_name, bs = "re"),
  data = spei_long,
  method = "ML",
  select = TRUE
)

smooth_coast_model <- mgcv::gam(
  WUE_T ~
    SPEI_timescale * coast_region +
    water_class +
    month_f +
    s(SPEI_value, by = spei_coast, k = 6) +
    s(site_name, bs = "re"),
  data = spei_long,
  method = "ML",
  select = TRUE
)

saveRDS(smooth_model, file.path(output_dir, "Q3_WUE_T_SPEI_smooth_gam_model.rds"))
saveRDS(smooth_coast_model, file.path(output_dir, "Q3_WUE_T_SPEI_coast_threshold_gam_model.rds"))

model_comparison <- data.frame(
  model = c("ecosystem_smooth_gam", "coast_threshold_gam"),
  AIC = c(AIC(smooth_model), AIC(smooth_coast_model)),
  BIC = c(BIC(smooth_model), BIC(smooth_coast_model)),
  logLik = c(
    as.numeric(stats::logLik(smooth_model)),
    as.numeric(stats::logLik(smooth_coast_model))
  ),
  deviance_explained = c(
    summary(smooth_model)$dev.expl,
    summary(smooth_coast_model)$dev.expl
  ),
  adjusted_r_squared = c(
    summary(smooth_model)$r.sq,
    summary(smooth_coast_model)$r.sq
  )
)
write_table(model_comparison, "Q3_WUE_T_SPEI_model_comparison.csv")

smooth_table <- as.data.frame(summary(smooth_model)$s.table)
smooth_table$smooth_term <- rownames(smooth_table)
rownames(smooth_table) <- NULL
smooth_table <- smooth_table[, c("smooth_term", setdiff(names(smooth_table), "smooth_term"))]
write_table(smooth_table, "Q3_WUE_T_SPEI_smooth_terms.csv")

parametric_table <- as.data.frame(summary(smooth_model)$p.table)
parametric_table$term <- rownames(parametric_table)
rownames(parametric_table) <- NULL
parametric_table <- parametric_table[, c("term", setdiff(names(parametric_table), "term"))]
write_table(parametric_table, "Q3_WUE_T_SPEI_smooth_parametric_terms.csv")

coast_smooth_table <- as.data.frame(summary(smooth_coast_model)$s.table)
coast_smooth_table$smooth_term <- rownames(coast_smooth_table)
rownames(coast_smooth_table) <- NULL
coast_smooth_table <- coast_smooth_table[, c("smooth_term", setdiff(names(coast_smooth_table), "smooth_term"))]
write_table(coast_smooth_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")

coast_parametric_table <- as.data.frame(summary(smooth_coast_model)$p.table)
coast_parametric_table$term <- rownames(coast_parametric_table)
rownames(coast_parametric_table) <- NULL
coast_parametric_table <- coast_parametric_table[, c("term", setdiff(names(coast_parametric_table), "term"))]
write_table(coast_parametric_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv")

# ---------------------------------------------------------------------------
# Predictions and plots
# ---------------------------------------------------------------------------

ecosystem_colors <- c(
  "Upland" = "#800080",
  "Freshwater" = "#0000FF",
  "Saline" = "#FFA500"
)

coast_colors <- c(
  "Atlantic Coast" = "#008B8B",
  "Pacific Coast" = "#C44E52",
  "Gulf Coast" = "#8C564B",
  "AK Coast" = "#4D4D4D"
)

spei_sequence <- seq(
  stats::quantile(spei_long$SPEI_value, 0.02),
  stats::quantile(spei_long$SPEI_value, 0.98),
  length.out = 120
)
prediction_grid <- expand.grid(
  SPEI_timescale = levels(spei_long$SPEI_timescale),
  water_class = levels(spei_long$water_class),
  SPEI_value = spei_sequence,
  month_f = factor("7", levels = levels(spei_long$month_f)),
  site_name = levels(spei_long$site_name)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
prediction_grid$SPEI_timescale <- factor(prediction_grid$SPEI_timescale, levels = levels(spei_long$SPEI_timescale))
prediction_grid$water_class <- factor(prediction_grid$water_class, levels = levels(spei_long$water_class))
prediction_grid$site_name <- factor(prediction_grid$site_name, levels = levels(spei_long$site_name))
prediction_grid$spei_ecosystem <- interaction(
  prediction_grid$SPEI_timescale,
  prediction_grid$water_class,
  sep = "__",
  drop = TRUE
)
prediction_grid$spei_ecosystem <- factor(prediction_grid$spei_ecosystem, levels = levels(spei_long$spei_ecosystem))

smooth_pred <- predict(
  smooth_model,
  newdata = prediction_grid,
  type = "link",
  se.fit = TRUE,
  exclude = "s(site_name)"
)
prediction_grid$predicted_WUE_T <- as.numeric(smooth_pred$fit)
prediction_grid$predicted_se <- as.numeric(smooth_pred$se.fit)
prediction_grid$predicted_lower <- prediction_grid$predicted_WUE_T - 1.96 * prediction_grid$predicted_se
prediction_grid$predicted_upper <- prediction_grid$predicted_WUE_T + 1.96 * prediction_grid$predicted_se
write_table(prediction_grid, "Q3_WUE_T_SPEI_smooth_predictions.csv")

near_normal_prediction_baseline <- prediction_grid[
  prediction_grid$SPEI_value >= -1 & prediction_grid$SPEI_value <= 1,
] |>
  dplyr::group_by(SPEI_timescale, water_class) |>
  dplyr::summarise(
    predicted_near_normal_WUE_T = mean(predicted_WUE_T),
    .groups = "drop"
  )

prediction_impact <- dplyr::left_join(
  prediction_grid,
  near_normal_prediction_baseline,
  by = c("SPEI_timescale", "water_class")
)
prediction_impact$predicted_change <- prediction_impact$predicted_WUE_T -
  prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_pct_change <- 100 * prediction_impact$predicted_change /
  prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_lower_pct <- 100 *
  (prediction_impact$predicted_lower - prediction_impact$predicted_near_normal_WUE_T) /
  prediction_impact$predicted_near_normal_WUE_T
prediction_impact$predicted_upper_pct <- 100 *
  (prediction_impact$predicted_upper - prediction_impact$predicted_near_normal_WUE_T) /
  prediction_impact$predicted_near_normal_WUE_T
prediction_impact$impact_class_10pct <- factor(
  ifelse(
    prediction_impact$predicted_pct_change <= -10, "WUE decrease",
    ifelse(prediction_impact$predicted_pct_change >= 10, "WUE increase", "No meaningful change")
  ),
  levels = c("WUE decrease", "No meaningful change", "WUE increase")
)

impact_thresholds <- prediction_impact |>
  dplyr::group_by(SPEI_timescale, water_class) |>
  dplyr::summarise(
    near_normal_WUE_T = dplyr::first(predicted_near_normal_WUE_T),
    dry_decrease_threshold = ifelse(
      any(SPEI_value < -1 & predicted_pct_change <= -10),
      max(SPEI_value[SPEI_value < -1 & predicted_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    dry_increase_threshold = ifelse(
      any(SPEI_value < -1 & predicted_pct_change >= 10),
      max(SPEI_value[SPEI_value < -1 & predicted_pct_change >= 10], na.rm = TRUE),
      NA_real_
    ),
    wet_decrease_threshold = ifelse(
      any(SPEI_value > 1 & predicted_pct_change <= -10),
      min(SPEI_value[SPEI_value > 1 & predicted_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    wet_increase_threshold = ifelse(
      any(SPEI_value > 1 & predicted_pct_change >= 10),
      min(SPEI_value[SPEI_value > 1 & predicted_pct_change >= 10], na.rm = TRUE),
      NA_real_
    ),
    min_predicted_pct_change = min(predicted_pct_change),
    max_predicted_pct_change = max(predicted_pct_change),
    .groups = "drop"
  )

write_table(near_normal_prediction_baseline, "Q3_WUE_T_SPEI_predicted_near_normal_baselines.csv")
write_table(prediction_impact, "Q3_WUE_T_SPEI_predicted_impact_classes.csv")
write_table(impact_thresholds, "Q3_WUE_T_SPEI_predicted_impact_thresholds_10pct.csv")

coast_prediction_grid <- expand.grid(
  SPEI_timescale = levels(spei_long$SPEI_timescale),
  coast_region = levels(spei_long$coast_region),
  SPEI_value = spei_sequence,
  water_class = factor("Upland", levels = levels(spei_long$water_class)),
  month_f = factor("7", levels = levels(spei_long$month_f)),
  site_name = levels(spei_long$site_name)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
coast_prediction_grid$SPEI_timescale <- factor(
  coast_prediction_grid$SPEI_timescale,
  levels = levels(spei_long$SPEI_timescale)
)
coast_prediction_grid$coast_region <- factor(
  coast_prediction_grid$coast_region,
  levels = levels(spei_long$coast_region)
)
coast_prediction_grid$water_class <- factor(
  coast_prediction_grid$water_class,
  levels = levels(spei_long$water_class)
)
coast_prediction_grid$site_name <- factor(
  coast_prediction_grid$site_name,
  levels = levels(spei_long$site_name)
)
coast_prediction_grid$spei_coast <- interaction(
  coast_prediction_grid$SPEI_timescale,
  coast_prediction_grid$coast_region,
  sep = "__",
  drop = TRUE
)
coast_prediction_grid$spei_coast <- factor(coast_prediction_grid$spei_coast, levels = levels(spei_long$spei_coast))

coast_smooth_pred <- predict(
  smooth_coast_model,
  newdata = coast_prediction_grid,
  type = "link",
  se.fit = TRUE,
  exclude = "s(site_name)"
)
coast_prediction_grid$predicted_WUE_T <- as.numeric(coast_smooth_pred$fit)
coast_prediction_grid$predicted_se <- as.numeric(coast_smooth_pred$se.fit)
coast_prediction_grid$predicted_lower <- coast_prediction_grid$predicted_WUE_T - 1.96 * coast_prediction_grid$predicted_se
coast_prediction_grid$predicted_upper <- coast_prediction_grid$predicted_WUE_T + 1.96 * coast_prediction_grid$predicted_se
write_table(coast_prediction_grid, "Q3_WUE_T_SPEI_coast_threshold_GAM_predictions_upland.csv")

coast_near_normal_prediction_baseline <- coast_prediction_grid[
  coast_prediction_grid$SPEI_value >= -1 & coast_prediction_grid$SPEI_value <= 1,
] |>
  dplyr::group_by(SPEI_timescale, coast_region) |>
  dplyr::summarise(
    predicted_near_normal_WUE_T = mean(predicted_WUE_T),
    .groups = "drop"
  )

coast_prediction_impact <- dplyr::left_join(
  coast_prediction_grid,
  coast_near_normal_prediction_baseline,
  by = c("SPEI_timescale", "coast_region")
)
coast_prediction_impact$predicted_change <- coast_prediction_impact$predicted_WUE_T -
  coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_pct_change <- 100 * coast_prediction_impact$predicted_change /
  coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_lower_pct <- 100 *
  (coast_prediction_impact$predicted_lower - coast_prediction_impact$predicted_near_normal_WUE_T) /
  coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$predicted_upper_pct <- 100 *
  (coast_prediction_impact$predicted_upper - coast_prediction_impact$predicted_near_normal_WUE_T) /
  coast_prediction_impact$predicted_near_normal_WUE_T
coast_prediction_impact$impact_class_10pct <- factor(
  ifelse(
    coast_prediction_impact$predicted_pct_change <= -10, "WUE decrease",
    ifelse(coast_prediction_impact$predicted_pct_change >= 10, "WUE increase", "No meaningful change")
  ),
  levels = c("WUE decrease", "No meaningful change", "WUE increase")
)

coast_impact_thresholds <- coast_prediction_impact |>
  dplyr::group_by(SPEI_timescale, coast_region) |>
  dplyr::summarise(
    near_normal_WUE_T = dplyr::first(predicted_near_normal_WUE_T),
    dry_decrease_threshold = ifelse(
      any(SPEI_value < -1 & predicted_pct_change <= -10),
      max(SPEI_value[SPEI_value < -1 & predicted_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    dry_increase_threshold = ifelse(
      any(SPEI_value < -1 & predicted_pct_change >= 10),
      max(SPEI_value[SPEI_value < -1 & predicted_pct_change >= 10], na.rm = TRUE),
      NA_real_
    ),
    wet_decrease_threshold = ifelse(
      any(SPEI_value > 1 & predicted_pct_change <= -10),
      min(SPEI_value[SPEI_value > 1 & predicted_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    wet_increase_threshold = ifelse(
      any(SPEI_value > 1 & predicted_pct_change >= 10),
      min(SPEI_value[SPEI_value > 1 & predicted_pct_change >= 10], na.rm = TRUE),
      NA_real_
    ),
    min_predicted_pct_change = min(predicted_pct_change),
    max_predicted_pct_change = max(predicted_pct_change),
    .groups = "drop"
  )

write_table(coast_near_normal_prediction_baseline, "Q3_WUE_T_SPEI_coast_threshold_GAM_near_normal_baselines_upland.csv")
write_table(coast_prediction_impact, "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv")
write_table(coast_impact_thresholds, "Q3_WUE_T_SPEI_coast_threshold_GAM_impact_thresholds_10pct_upland.csv")

find_spei_threshold <- function(group_data, pct_column, threshold_pct, anomaly_side, impact_direction) {
  pct_values <- group_data[[pct_column]]
  side_filter <- if (anomaly_side == "dry") group_data$SPEI_value < -1 else group_data$SPEI_value > 1
  direction_filter <- if (impact_direction == "decrease") {
    pct_values <= -threshold_pct
  } else {
    pct_values >= threshold_pct
  }
  threshold_values <- group_data$SPEI_value[side_filter & direction_filter]
  if (length(threshold_values) == 0) return(NA_real_)
  if (anomaly_side == "dry") {
    max(threshold_values, na.rm = TRUE)
  } else {
    min(threshold_values, na.rm = TRUE)
  }
}

coast_impact_threshold_markers <- do.call(
  rbind,
  lapply(
    split(coast_prediction_impact, list(coast_prediction_impact$SPEI_timescale, coast_prediction_impact$coast_region), drop = TRUE),
    function(group_data) {
      do.call(
        rbind,
        lapply(c(5, 10, 20), function(threshold_pct) {
          dry_decrease <- group_data$SPEI_value[
            group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -threshold_pct
          ]
          dry_increase <- group_data$SPEI_value[
            group_data$SPEI_value < -1 & group_data$predicted_pct_change >= threshold_pct
          ]
          wet_decrease <- group_data$SPEI_value[
            group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -threshold_pct
          ]
          wet_increase <- group_data$SPEI_value[
            group_data$SPEI_value > 1 & group_data$predicted_pct_change >= threshold_pct
          ]
          data.frame(
            SPEI_timescale = group_data$SPEI_timescale[1],
            coast_region = group_data$coast_region[1],
            water_class = group_data$water_class[1],
            threshold_pct = threshold_pct,
            impact_direction = c("decrease", "increase", "decrease", "increase"),
            anomaly_side = c("dry", "dry", "wet", "wet"),
            SPEI_threshold = c(
              ifelse(length(dry_decrease) > 0, max(dry_decrease, na.rm = TRUE), NA_real_),
              ifelse(length(dry_increase) > 0, max(dry_increase, na.rm = TRUE), NA_real_),
              ifelse(length(wet_decrease) > 0, min(wet_decrease, na.rm = TRUE), NA_real_),
              ifelse(length(wet_increase) > 0, min(wet_increase, na.rm = TRUE), NA_real_)
            ),
            pct_change_threshold = c(
              -threshold_pct,
              threshold_pct,
              -threshold_pct,
              threshold_pct
            )
          )
        })
      )
    }
  )
)
coast_impact_threshold_markers <- coast_impact_threshold_markers[
  !is.na(coast_impact_threshold_markers$SPEI_threshold),
]
coast_impact_threshold_markers$SPEI_threshold_lower <- NA_real_
coast_impact_threshold_markers$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(coast_impact_threshold_markers))) {
  marker_row <- coast_impact_threshold_markers[i, ]
  marker_group <- coast_prediction_impact[
    coast_prediction_impact$SPEI_timescale == marker_row$SPEI_timescale &
      coast_prediction_impact$coast_region == marker_row$coast_region,
  ]
  lower_threshold <- find_spei_threshold(
    marker_group,
    "predicted_lower_pct",
    marker_row$threshold_pct,
    marker_row$anomaly_side,
    marker_row$impact_direction
  )
  upper_threshold <- find_spei_threshold(
    marker_group,
    "predicted_upper_pct",
    marker_row$threshold_pct,
    marker_row$anomaly_side,
    marker_row$impact_direction
  )
  threshold_range <- c(marker_row$SPEI_threshold, lower_threshold, upper_threshold)
  threshold_range <- threshold_range[!is.na(threshold_range)]
  coast_impact_threshold_markers$SPEI_threshold_lower[i] <- min(threshold_range)
  coast_impact_threshold_markers$SPEI_threshold_upper[i] <- max(threshold_range)
}
coast_impact_threshold_markers$threshold_pct <- factor(
  coast_impact_threshold_markers$threshold_pct,
  levels = c(5, 10, 20),
  labels = c("5%", "10%", "20%")
)
coast_impact_threshold_markers$SPEI_timescale <- factor(
  coast_impact_threshold_markers$SPEI_timescale,
  levels = levels(spei_long$SPEI_timescale)
)
coast_impact_threshold_markers$coast_region <- factor(
  coast_impact_threshold_markers$coast_region,
  levels = levels(spei_long$coast_region)
)

coast_threshold_summary_table <- coast_impact_threshold_markers |>
  dplyr::arrange(SPEI_timescale, coast_region, anomaly_side, impact_direction, threshold_pct)

write_table(coast_impact_threshold_markers, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv")
write_table(coast_threshold_summary_table, "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv")

gam_impact_threshold_markers <- do.call(
  rbind,
  lapply(
    split(prediction_impact, list(prediction_impact$SPEI_timescale, prediction_impact$water_class), drop = TRUE),
    function(group_data) {
      do.call(
        rbind,
        lapply(c(5, 10, 20), function(threshold_pct) {
          dry_decrease <- group_data$SPEI_value[
            group_data$SPEI_value < -1 & group_data$predicted_pct_change <= -threshold_pct
          ]
          dry_increase <- group_data$SPEI_value[
            group_data$SPEI_value < -1 & group_data$predicted_pct_change >= threshold_pct
          ]
          wet_decrease <- group_data$SPEI_value[
            group_data$SPEI_value > 1 & group_data$predicted_pct_change <= -threshold_pct
          ]
          wet_increase <- group_data$SPEI_value[
            group_data$SPEI_value > 1 & group_data$predicted_pct_change >= threshold_pct
          ]
          data.frame(
            SPEI_timescale = group_data$SPEI_timescale[1],
            water_class = group_data$water_class[1],
            threshold_pct = threshold_pct,
            impact_direction = c("decrease", "increase", "decrease", "increase"),
            anomaly_side = c("dry", "dry", "wet", "wet"),
            SPEI_threshold = c(
              ifelse(length(dry_decrease) > 0, max(dry_decrease, na.rm = TRUE), NA_real_),
              ifelse(length(dry_increase) > 0, max(dry_increase, na.rm = TRUE), NA_real_),
              ifelse(length(wet_decrease) > 0, min(wet_decrease, na.rm = TRUE), NA_real_),
              ifelse(length(wet_increase) > 0, min(wet_increase, na.rm = TRUE), NA_real_)
            ),
            pct_change_threshold = c(
              -threshold_pct,
              threshold_pct,
              -threshold_pct,
              threshold_pct
            )
          )
        })
      )
    }
  )
)
gam_impact_threshold_markers <- gam_impact_threshold_markers[
  !is.na(gam_impact_threshold_markers$SPEI_threshold),
]
gam_impact_threshold_markers$SPEI_threshold_lower <- NA_real_
gam_impact_threshold_markers$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(gam_impact_threshold_markers))) {
  marker_row <- gam_impact_threshold_markers[i, ]
  marker_group <- prediction_impact[
    prediction_impact$SPEI_timescale == marker_row$SPEI_timescale &
      prediction_impact$water_class == marker_row$water_class,
  ]
  lower_threshold <- find_spei_threshold(
    marker_group,
    "predicted_lower_pct",
    marker_row$threshold_pct,
    marker_row$anomaly_side,
    marker_row$impact_direction
  )
  upper_threshold <- find_spei_threshold(
    marker_group,
    "predicted_upper_pct",
    marker_row$threshold_pct,
    marker_row$anomaly_side,
    marker_row$impact_direction
  )
  threshold_range <- c(marker_row$SPEI_threshold, lower_threshold, upper_threshold)
  threshold_range <- threshold_range[!is.na(threshold_range)]
  gam_impact_threshold_markers$SPEI_threshold_lower[i] <- min(threshold_range)
  gam_impact_threshold_markers$SPEI_threshold_upper[i] <- max(threshold_range)
}
gam_impact_threshold_markers$threshold_pct <- factor(
  gam_impact_threshold_markers$threshold_pct,
  levels = c(5, 10, 20),
  labels = c("5%", "10%", "20%")
)
gam_impact_threshold_markers$SPEI_timescale <- factor(
  gam_impact_threshold_markers$SPEI_timescale,
  levels = levels(spei_long$SPEI_timescale)
)
gam_impact_threshold_markers$water_class <- factor(
  gam_impact_threshold_markers$water_class,
  levels = levels(spei_long$water_class)
)

gam_threshold_summary_table <- gam_impact_threshold_markers |>
  dplyr::arrange(SPEI_timescale, water_class, anomaly_side, impact_direction, threshold_pct)

write_table(gam_impact_threshold_markers, "Q3_WUE_T_SPEI_GAM_impact_threshold_markers_5_10_20pct.csv")
write_table(gam_threshold_summary_table, "Q3_WUE_T_SPEI_GAM_threshold_summary_5_10_20pct.csv")

plot_site_slopes <- ggplot2::ggplot(
  site_slopes,
  ggplot2::aes(x = SPEI_timescale, y = SPEI_slope, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, color = "gray35") +
  ggplot2::geom_boxplot(outlier.alpha = 0.15) +
  ggplot2::labs(
    title = expression("Site-Level " * WUE[T] * " Sensitivity Slopes"),
    x = "SPEI timescale",
    y = expression("Site-level slope of " * WUE[T] * " vs SPEI"),
    fill = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

save_plot(plot_site_slopes, "Q3_plot_02_site_level_sensitivity_slopes.png", width = 11, height = 7)

plot_smooth <- ggplot2::ggplot(
  prediction_grid,
  ggplot2::aes(x = SPEI_value, y = predicted_WUE_T, color = water_class, fill = water_class)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = predicted_lower, ymax = predicted_upper),
    alpha = 0.14,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::facet_wrap(~ SPEI_timescale, scales = "free_x", ncol = 2) +
  ggplot2::labs(
    title = expression("GAM " * WUE[T] * " Response Across SPEI Gradients"),
    subtitle = "Dashed lines mark SPEI anomaly thresholds at -1 and 1",
    x = "SPEI",
    y = expression("Predicted " * WUE[T]),
    color = "Ecosystem class",
    fill = "Ecosystem class"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(plot_smooth, "Q3_plot_03_smooth_WUE_T_response_by_SPEI.png", width = 12, height = 12)

selected_timescales <- c("SPEI_1", "SPEI_3", "SPEI_48")
plot_selected_smooth <- ggplot2::ggplot(
  prediction_impact[prediction_impact$SPEI_timescale %in% selected_timescales, ],
  ggplot2::aes(x = SPEI_value, y = predicted_pct_change, color = water_class, fill = water_class)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45") +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray35") +
  ggplot2::geom_hline(yintercept = c(-10, 10), linetype = "dotted", color = "gray35") +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = 100 * (predicted_lower - predicted_near_normal_WUE_T) / predicted_near_normal_WUE_T,
      ymax = 100 * (predicted_upper - predicted_near_normal_WUE_T) / predicted_near_normal_WUE_T
    ),
    alpha = 0.14,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_segment(
    data = gam_impact_threshold_markers[
      gam_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold_lower,
      xend = SPEI_threshold_upper,
      y = pct_change_threshold,
      yend = pct_change_threshold,
      color = water_class
    ),
    linewidth = 1.05,
    alpha = 0.72,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = gam_impact_threshold_markers[
      gam_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      xend = SPEI_threshold,
      y = 0,
      yend = pct_change_threshold,
      color = water_class,
      linetype = threshold_pct
    ),
    linewidth = 0.65,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = gam_impact_threshold_markers[
      gam_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      y = pct_change_threshold,
      color = water_class,
      shape = threshold_pct
    ),
    size = 3.1,
    stroke = 1.1,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~ SPEI_timescale, scales = "free_x", nrow = 1) +
  ggplot2::labs(
    title = expression("GAM-Predicted " * WUE[T] * " Change from Near-Normal Conditions"),
    subtitle = "Markers show SPEI values where predicted change reaches 5%, 10%, and 20%; horizontal bars show prediction-interval threshold ranges",
    x = "SPEI",
    y = expression("Predicted " * WUE[T] * " change from near-normal (%)"),
    color = "Ecosystem class",
    fill = "Ecosystem class",
    shape = "Change threshold",
    linetype = "Change threshold"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  ggplot2::scale_shape_manual(values = c("5%" = 21, "10%" = 24, "20%" = 22)) +
  theme_wue()

save_plot(plot_selected_smooth, "Q3_plot_06_WUE_T_selected_SPEI_timescales.png", width = 14, height = 6)

coast_panel_labels <- data.frame(
  coast_region = factor(coast_region_levels, levels = coast_region_levels),
  SPEI_timescale = factor("SPEI_1", levels = levels(spei_long$SPEI_timescale)),
  coast_label = c("Atlantic", "Pacific", "Gulf", "AK")
)

plot_coast_threshold <- ggplot2::ggplot(
  coast_prediction_impact[coast_prediction_impact$SPEI_timescale %in% selected_timescales, ],
  ggplot2::aes(x = SPEI_value, y = predicted_pct_change, color = coast_region)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45") +
  ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "gray35") +
  ggplot2::geom_hline(yintercept = c(-20, -10, -5, 5, 10, 20), linetype = "dotted", color = "gray35") +
  ggplot2::geom_label(
    data = coast_panel_labels,
    ggplot2::aes(label = coast_label),
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.18,
    size = 5.2,
    linewidth = 0,
    fill = "white",
    alpha = 0.86,
    color = "gray20",
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = predicted_lower_pct, ymax = predicted_upper_pct, fill = coast_region),
    alpha = 0.10,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_segment(
    data = coast_impact_threshold_markers[
      coast_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold_lower,
      xend = SPEI_threshold_upper,
      y = pct_change_threshold,
      yend = pct_change_threshold,
      color = coast_region
    ),
    linewidth = 1.15,
    alpha = 0.72,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = coast_impact_threshold_markers[
      coast_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      xend = SPEI_threshold,
      y = 0,
      yend = pct_change_threshold,
      color = coast_region,
      linetype = threshold_pct
    ),
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = coast_impact_threshold_markers[
      coast_impact_threshold_markers$SPEI_timescale %in% selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      y = pct_change_threshold,
      color = coast_region,
      shape = threshold_pct
    ),
    size = 3.6,
    stroke = 1.1,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_grid(
    coast_region ~ SPEI_timescale,
    scales = "free"
  ) +
  ggplot2::labs(
    title = expression(bold("Coast-Specific GAM Thresholds for Upland " * WUE[T] * " Response")),
    subtitle = NULL,
    x = "SPEI",
    y = expression("Predicted upland " * WUE[T] * " change from near-normal (%)"),
    color = "Coast region",
    fill = "Coast region",
    shape = "Change threshold",
    linetype = "Change threshold"
  ) +
  ggplot2::scale_color_manual(values = coast_colors) +
  ggplot2::scale_fill_manual(values = coast_colors) +
  ggplot2::scale_shape_manual(values = c("5%" = 21, "10%" = 24, "20%" = 22)) +
  ggplot2::guides(color = "none", fill = "none") +
  theme_wue() +
  ggplot2::theme(
    strip.text.y = ggplot2::element_blank(),
    strip.background.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5, 12, 5, 5)
  )

save_plot(plot_coast_threshold, "Q3_plot_10_coast_specific_GAM_SPEI_thresholds_upland.png", width = 12, height = 10)

top_sensitive_sites <- site_sensitivity_rank |>
  dplyr::filter(SPEI_timescale %in% selected_timescales) |>
  dplyr::group_by(SPEI_timescale) |>
  dplyr::slice_max(mean_abs_pct_change, n = 12, with_ties = FALSE) |>
  dplyr::ungroup()
top_sensitive_sites$site_name <- factor(
  top_sensitive_sites$site_name,
  levels = unique(top_sensitive_sites$site_name[order(top_sensitive_sites$mean_abs_pct_change)])
)

plot_sensitive_sites <- ggplot2::ggplot(
  top_sensitive_sites,
  ggplot2::aes(x = mean_abs_pct_change, y = site_name, fill = coast_region)
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::facet_wrap(~ SPEI_timescale, scales = "free_y") +
  ggplot2::labs(
    title = "Most Sensitive Sites",
    subtitle = expression("Ranked by mean absolute " * WUE[T] * " percent change from each site's near-normal baseline"),
    x = expression("Mean absolute " * WUE[T] * " change (%)"),
    y = NULL,
    fill = "Coast region"
  ) +
  ggplot2::scale_fill_manual(values = coast_colors) +
  theme_wue()

save_plot(plot_sensitive_sites, "Q3_plot_07_most_sensitive_sites.png", width = 13, height = 8)

all_system_sensitivity <- site_sensitivity_rank |>
  dplyr::filter(SPEI_timescale %in% selected_timescales, n_months >= 6)
site_sensitivity_order <- all_system_sensitivity |>
  dplyr::group_by(site_name) |>
  dplyr::summarise(
    max_mean_abs_pct_change = max(mean_abs_pct_change, na.rm = TRUE),
    .groups = "drop"
  )
all_system_sensitivity$site_name <- factor(
  all_system_sensitivity$site_name,
  levels = site_sensitivity_order$site_name[order(site_sensitivity_order$max_mean_abs_pct_change)]
)

plot_all_systems <- ggplot2::ggplot(
  all_system_sensitivity,
  ggplot2::aes(x = mean_abs_pct_change, y = site_name, color = water_class, shape = SPEI_timescale)
) +
  ggplot2::geom_vline(xintercept = c(5, 10, 20), linetype = "dotted", color = "gray65") +
  ggplot2::geom_point(
    size = 2.8,
    alpha = 0.88,
    position = ggplot2::position_jitter(height = 0.12, width = 0)
  ) +
  ggplot2::labs(
    title = "Sensitivity of All Systems",
    subtitle = expression("Each point is one site-by-SPEI timescale; dotted lines mark 5%, 10%, and 20% mean absolute " * WUE[T] * " change"),
    x = expression("Mean absolute " * WUE[T] * " change from near-normal (%)"),
    y = NULL,
    color = "Ecosystem class",
    shape = "SPEI timescale"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  theme_wue(base_size = 11) +
  ggplot2::theme(
    legend.position = "right",
    axis.text.y = ggplot2::element_text(size = 6),
    panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.25)
  )

save_plot(plot_all_systems, "Q3_plot_09_all_system_sensitivity_heatmap.png", width = 13, height = 13)

coastline_sensitivity_summary <- all_system_sensitivity |>
  dplyr::group_by(coast_region, SPEI_timescale) |>
  dplyr::summarise(
    n_sites = dplyr::n_distinct(site_name),
    median_abs_pct_change = median(mean_abs_pct_change, na.rm = TRUE),
    q25_abs_pct_change = stats::quantile(mean_abs_pct_change, 0.25, na.rm = TRUE),
    q75_abs_pct_change = stats::quantile(mean_abs_pct_change, 0.75, na.rm = TRUE),
    mean_abs_pct_change = mean(mean_abs_pct_change, na.rm = TRUE),
    .groups = "drop"
  )
write_table(coastline_sensitivity_summary, "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv")

plot_coastline_sensitivity <- ggplot2::ggplot(
  all_system_sensitivity,
  ggplot2::aes(x = mean_abs_pct_change, y = coast_region)
) +
  ggplot2::geom_vline(xintercept = c(5, 10, 20), linetype = "dotted", color = "gray65") +
  ggplot2::geom_boxplot(
    ggplot2::aes(fill = coast_region),
    width = 0.58,
    alpha = 0.18,
    outlier.shape = NA,
    color = "gray35"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = water_class),
    size = 3.2,
    alpha = 0.82,
    position = ggplot2::position_jitter(width = 0, height = 0.12)
  ) +
  ggplot2::facet_wrap(~ SPEI_timescale, nrow = 1) +
  ggplot2::labs(
    title = "Site Sensitivity Grouped by Coastline",
    subtitle = expression("Each point is one site; boxes summarize mean absolute " * WUE[T] * " change by coastline"),
    x = expression("Mean absolute " * WUE[T] * " change from near-normal (%)"),
    y = NULL,
    fill = "Coast region",
    color = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = coast_colors) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.y = ggplot2::element_text(size = 11)
  )

save_plot(plot_coastline_sensitivity, "Q3_plot_11_coastline_grouped_site_sensitivity.png", width = 14, height = 6)

plot_coast_threshold_composite <- plot_coast_threshold +
  ggplot2::labs(
    subtitle = NULL,
    y = expression("Upland " * WUE[T] * " change (%)")
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    axis.title = ggplot2::element_text(size = 17),
    axis.text = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 15, face = "bold"),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    legend.key.size = ggplot2::unit(0.85, "lines"),
    plot.margin = ggplot2::margin(5, 28, 5, 5)
  )

plot_coastline_sensitivity_composite <- plot_coastline_sensitivity +
  ggplot2::labs(subtitle = NULL) +
  ggplot2::guides(fill = "none") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 20, face = "bold"),
    axis.title = ggplot2::element_text(size = 17),
    axis.text = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 15, face = "bold"),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    legend.key.size = ggplot2::unit(0.85, "lines")
  )

composite <- plot_coast_threshold_composite / plot_coastline_sensitivity_composite +
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 22),
    plot.tag.position = c(0.005, 0.98),
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_plot(composite, "Q3_plot_04_WUE_T_SPEI_sensitivity_multipanel.png", width = 10, height = 13.5)

capture_to_file({
  cat("Q3 WUE_T sensitivity to SPEI gradients\n")
  cat("Input:", input_file, "\n")
  cat("Output:", output_dir, "\n\n")
  cat("Rows in source monthly data:", nrow(monthly), "\n")
  cat("Rows after complete-case and ecosystem filters:", nrow(model_base), "\n")
  cat("Rows in SPEI long data:", nrow(spei_long), "\n")
  cat("Sites:", length(unique(spei_long$site_name)), "\n")
  cat("Years:", min(spei_long$Year), "-", max(spei_long$Year), "\n\n")

  cat("Primary model formula:\n")
  print(formula(smooth_model))

  cat("\nGAM model information:\n")
  print(model_comparison)

  cat("\nSmooth GAM summary:\n")
  print(summary(smooth_model))

  cat("\nCoast-threshold GAM summary:\n")
  print(summary(smooth_coast_model))

  cat("\nGAM-predicted ecological-impact threshold markers:\n")
  print(gam_threshold_summary_table)

  cat("\nCoast-specific GAM-predicted upland ecological-impact threshold markers:\n")
  print(coast_threshold_summary_table)
}, "Q3_WUE_T_SPEI_model_summary.txt")

message("Done. Q3 WUE_T sensitivity outputs written to: ", output_dir)
