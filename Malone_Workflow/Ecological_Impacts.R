# Ecological_Impacts.R
#
# Goal:
#   Build an Ecological Drought Index (EDI) that follows the final Q3 logic:
#   ecological impact is the modelled magnitude and direction of WUE_T change
#   from near-normal SPEI conditions, with the primary inference focused on
#   upland systems and coast-specific GAM thresholds.
#
# Primary Q3-aligned approach:
#   - Read the final Q3 coast-threshold GAM prediction and threshold outputs.
#   - Use those outputs as the source of truth for population-level upland
#     WUE_T response by coast.
#   - Classify response magnitude using 5%, 10%, and 20% thresholds and retain
#     response direction (increase vs decrease).
#
# Secondary/legacy diagnostics:
#   The older logistic EDI, P(WUE_T below baseline | SPEI-3), is retained as a
#   diagnostic output only. It should not be framed as the primary ecological
#   impact result because WUE_T decreases and increases can both be meaningful
#   and because Q3 showed the strongest signal is nonlinear and coast-specific.
#
# Inputs:
#   Water_WUE/data/WUE_CUE_monthly_merged_indices_clean.csv
#   Malone_Workflow/results/Q2_WUE_performance/Q2_site_resistance.csv
#   Malone_Workflow/results/Q2_WUE_performance/Q2_site_recovery.csv
#   Malone_Workflow/results/Q3_WUE_T_SPEI_sensitivity/
#     Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv
#     Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv
#
# Outputs:
#   Malone_Workflow/results/Ecological_Impacts/
#     EDI_logistic_model_summary.txt
#     EDI_Q3_aligned_model_summary.txt
#     EDI_Q3_aligned_coast_thresholds_5_10_20pct.csv
#     EDI_Q3_aligned_prediction_scores_upland.csv
#     EDI_Q3_aligned_monthly_scores.csv
#     EDI_plot_08_Q3_aligned_composite.png
#     EDI_logistic_coefficients.csv
#     EDI_threshold_summary.csv
#     EDI_monthly_site_scores.csv
#     EDI_validation_by_class.csv
#     EDI_plot_01_logistic_curve.png
#     EDI_plot_02_monthly_resistance_calibration.png
#     EDI_plot_03_spei_edi_by_ecosystem.png
#     EDI_plot_04_edi_timeseries_facet.png
#     EDI_plot_05_severity_class_proportions.png
#     EDI_plot_06_threshold_validation.png
#     EDI_plot_07_composite.png

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "tidyr", "ggplot2", "patchwork", "scales", "mgcv")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install with install.packages('%s').", pkg, pkg))
    }
  }
})

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  candidates <- c(
    "Ecological_Impacts.R",
    file.path("Malone_Workflow", "Ecological_Impacts.R"),
    file.path("Water_WUE", "Malone_Workflow", "Ecological_Impacts.R")
  )
  matches <- candidates[file.exists(candidates)]
  if (length(matches) > 0) {
    return(normalizePath(matches[1], mustWork = TRUE))
  }
  stop("Could not locate Ecological_Impacts.R. Run from Water_WUE or Water_WUE/Malone_Workflow, or call with Rscript path/to/Ecological_Impacts.R.")
}

write_table <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}

capture_to_file <- function(expr, filename) {
  sink(file.path(output_dir, filename))
  on.exit(sink(), add = TRUE)
  force(expr)
}

save_plot <- function(plot, filename, width = 11, height = 7) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = 300,
    bg       = "white"
  )
}

theme_wue <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom"
    )
}

logistic <- function(x, b0, b1) 1 / (1 + exp(-(b0 + b1 * x)))

# SPEI value at which P(impact) = p
spei_at_p <- function(p, b0, b1) (log(p / (1 - p)) - b0) / b1

classify_response_magnitude <- function(x) {
  out <- ifelse(
    x < 0.5, "Stable (<0.5%)",
    ifelse(
      x < 1, "Very mild (0.5-1%)",
      ifelse(
        x < 2.5, "Mild (1-2.5%)",
        ifelse(x < 5, "Moderate (2.5-5%)", "Strong (>=5%)")
      )
    )
  )
  factor(
    out,
    levels = c(
      "Stable (<0.5%)",
      "Very mild (0.5-1%)",
      "Mild (1-2.5%)",
      "Moderate (2.5-5%)",
      "Strong (>=5%)"
    )
  )
}

classify_response_direction <- function(x, stable_threshold = 0.5) {
  out <- ifelse(
    abs(x) < stable_threshold, "Stable",
    ifelse(x < 0, "WUE_T decrease", "WUE_T increase")
  )
  factor(out, levels = c("WUE_T decrease", "Stable", "WUE_T increase"))
}

classify_q3_response_magnitude <- function(x) {
  out <- ifelse(
    x < 5, "No meaningful change (<5%)",
    ifelse(
      x < 10, "Watch (5-10%)",
      ifelse(x < 20, "Stress (10-20%)", "Impact (>=20%)")
    )
  )
  factor(
    out,
    levels = c(
      "No meaningful change (<5%)",
      "Watch (5-10%)",
      "Stress (10-20%)",
      "Impact (>=20%)"
    )
  )
}

classify_q3_response_direction <- function(x, stable_threshold = 5) {
  out <- ifelse(
    abs(x) < stable_threshold, "No meaningful change",
    ifelse(x < 0, "WUE_T decrease", "WUE_T increase")
  )
  factor(out, levels = c("WUE_T decrease", "No meaningful change", "WUE_T increase"))
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

find_first_threshold <- function(group_data, threshold, side, direction = c("any", "increase", "decrease")) {
  direction <- match.arg(direction)
  side_filter <- if (side == "dry") group_data$SPEI_3 < -1 else group_data$SPEI_3 > 1
  response_filter <- switch(
    direction,
    any = abs(group_data$predicted_pct_change) >= threshold,
    increase = group_data$predicted_pct_change >= threshold,
    decrease = group_data$predicted_pct_change <= -threshold
  )
  threshold_values <- group_data$SPEI_3[side_filter & response_filter]
  if (length(threshold_values) == 0) return(NA_real_)
  if (side == "dry") max(threshold_values, na.rm = TRUE) else min(threshold_values, na.rm = TRUE)
}

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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

script_dir  <- dirname(get_script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
q2_dir      <- file.path(script_dir, "results", "Q2_WUE_performance")
q3_dir      <- file.path(script_dir, "results", "Q3_WUE_T_SPEI_sensitivity")
output_dir  <- file.path(script_dir, "results", "Ecological_Impacts")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

ecosystem_classes <- c("Upland", "Freshwater", "Saline")
required_cols     <- c("site_name", "Year", "month", "water_class",
                       "WUE_tra", "SPEI_3")

monthly <- read.csv(
  file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv"),
  stringsAsFactors = FALSE
)
monthly[] <- lapply(monthly, function(x) if (is.character(x)) trimws(x) else x)

resistance <- read.csv(file.path(q2_dir, "Q2_site_resistance.csv"),
                       stringsAsFactors = FALSE)
recovery   <- read.csv(file.path(q2_dir, "Q2_site_recovery.csv"),
                       stringsAsFactors = FALSE)

base_df <- monthly[
  complete.cases(monthly[, required_cols]) &
    monthly$site_name    != "" &
    monthly$water_class  %in% ecosystem_classes &
    is.finite(monthly$WUE_tra) &
    is.finite(monthly$SPEI_3),
]
base_df$WUE_T       <- base_df$WUE_tra
base_df$water_class <- factor(base_df$water_class, levels = ecosystem_classes)
base_df$month_f     <- factor(base_df$month, levels = 1:12)
base_df             <- base_df[order(base_df$site_name, base_df$Year, base_df$month), ]

message(sprintf("Monthly data: %d rows, %d sites",
                nrow(base_df), length(unique(base_df$site_name))))

# ---------------------------------------------------------------------------
# Step 1 — Compute per-site near-normal WUE_T baseline from SPEI_3
#
# Near-normal: -1 <= SPEI_3 <= 1 (same convention as Q2).
# Require >= 3 near-normal observations per site-month to establish a reliable
# seasonal baseline.
# ---------------------------------------------------------------------------

baselines <- base_df |>
  dplyr::filter(SPEI_3 >= -1, SPEI_3 <= 1) |>
  dplyr::group_by(site_name, month_f) |>
  dplyr::filter(dplyr::n() >= 3) |>
  dplyr::summarise(
    baseline_WUE_T    = mean(WUE_T),
    baseline_sd_WUE_T = sd(WUE_T),
    n_baseline_months = dplyr::n(),
    .groups = "drop"
  )

message(sprintf(
  "Near-normal site-month baselines computed for %d site-month combinations across %d sites",
  nrow(baselines), length(unique(baselines$site_name))
))

# Attach baseline to every row
base_df <- dplyr::left_join(base_df, baselines, by = c("site_name", "month_f"))

# Monthly resistance: WUE_T relative to that site's near-normal mean
# Only meaningful when a baseline exists; filter to those sites for modelling
base_df$monthly_resistance <- base_df$WUE_T / base_df$baseline_WUE_T
base_df$WUE_T_pct_change <- 100 * (base_df$WUE_T - base_df$baseline_WUE_T) / base_df$baseline_WUE_T

# ---------------------------------------------------------------------------
# Step 2 — Logistic calibration fitted to observed resistance
#
# For each site-month where SPEI_3 < 0 (any moisture deficit) and a baseline
# exists, define:
#   ecological_impact = 1  if monthly_resistance < 1  (WUE_T below baseline)
#                     = 0  if monthly_resistance >= 1
#
# Fit: P(ecological_impact = 1) ~ SPEI_3
#   — Pooled model (all ecosystems)
#   — Ecosystem-stratified models
#
# Using SPEI_3 < 0 rather than < -1 as the modelling domain gives the logistic
# more gradient to fit across; the -1 threshold emerges naturally from the
# fitted curve.
# ---------------------------------------------------------------------------

model_df <- base_df[
  !is.na(base_df$monthly_resistance) &
  is.finite(base_df$monthly_resistance) &
  base_df$SPEI_3 < 0,    # moisture-deficit months only
]
model_df$y <- as.integer(model_df$monthly_resistance < 1)

message(sprintf(
  "Logistic calibration data: %d site-months (SPEI_3 < 0), %d impacted (resistance < 1)",
  nrow(model_df), sum(model_df$y)
))

# Pooled model
pooled_model <- glm(
  y ~ SPEI_3,
  data   = model_df,
  family = binomial(link = "logit")
)

# Ecosystem-stratified models
eco_models <- setNames(
  lapply(ecosystem_classes, function(ec) {
    sub_df <- model_df[model_df$water_class == ec, ]
    if (sum(sub_df$y) < 5 | sum(1 - sub_df$y) < 5) {
      message(sprintf("  Skipping %s: too few events (%d impacted)", ec, sum(sub_df$y)))
      return(NULL)
    }
    tryCatch(
      glm(y ~ SPEI_3, data = sub_df, family = binomial(link = "logit")),
      error = function(e) { message("  Model failed for ", ec, ": ", e$message); NULL }
    )
  }),
  ecosystem_classes
)

# ---------------------------------------------------------------------------
# Step 3 — Apply EDI to all site-months
# ---------------------------------------------------------------------------

# Use ecosystem-specific model where available, pooled otherwise
base_df$EDI <- predict(pooled_model, newdata = base_df, type = "response")

for (ec in ecosystem_classes) {
  m <- eco_models[[ec]]
  if (!is.null(m)) {
    idx <- base_df$water_class == ec
    base_df$EDI[idx] <- predict(m, newdata = base_df[idx, ], type = "response")
  }
}

edi_breaks  <- c(0, 0.30, 0.50, 0.70, 1.0)
edi_labels  <- c("None", "Watch", "Stress", "Impact")
base_df$EDI_class <- cut(base_df$EDI,
                         breaks         = edi_breaks,
                         labels         = edi_labels,
                         include.lowest = TRUE)

# ---------------------------------------------------------------------------
# Step 4 — Response-based EDI: expected WUE_T change magnitude and direction
#
# The logistic EDI above asks whether WUE_T falls below baseline. The response
# EDI keeps the magnitude and direction of the WUE_T response:
#   EDI_response_magnitude = abs(predicted % change from site-month baseline)
#   EDI_response_direction = decrease, stable, or increase
# ---------------------------------------------------------------------------

response_model_df <- base_df[
  is.finite(base_df$WUE_T_pct_change) &
    is.finite(base_df$SPEI_3) &
    !is.na(base_df$baseline_WUE_T) &
    base_df$n_baseline_months >= 3,
]
response_model_df$site_name <- factor(response_model_df$site_name)
response_model_df$water_class <- factor(response_model_df$water_class, levels = ecosystem_classes)
response_model_df$month_f <- factor(response_model_df$month, levels = 1:12)

response_model <- mgcv::gam(
  WUE_T_pct_change ~
    water_class + month_f +
    s(SPEI_3, by = water_class, k = 6) +
    s(site_name, bs = "re"),
  data = response_model_df,
  method = "REML",
  select = TRUE
)

base_df$response_EDI_pct_change <- NA_real_
base_df$response_EDI_pct_change_lower <- NA_real_
base_df$response_EDI_pct_change_upper <- NA_real_
predictable_rows <- is.finite(base_df$SPEI_3) &
  !is.na(base_df$baseline_WUE_T) &
  base_df$n_baseline_months >= 3

response_pred <- predict(
  response_model,
  newdata = base_df[predictable_rows, ],
  se.fit = TRUE,
  exclude = "s(site_name)"
)
base_df$response_EDI_pct_change[predictable_rows] <- as.numeric(response_pred$fit)
base_df$response_EDI_pct_change_lower[predictable_rows] <-
  as.numeric(response_pred$fit - 1.96 * response_pred$se.fit)
base_df$response_EDI_pct_change_upper[predictable_rows] <-
  as.numeric(response_pred$fit + 1.96 * response_pred$se.fit)
base_df$response_EDI_magnitude <- abs(base_df$response_EDI_pct_change)
base_df$response_EDI_class <- classify_response_magnitude(base_df$response_EDI_magnitude)
base_df$response_EDI_direction <- classify_response_direction(base_df$response_EDI_pct_change)

spei_seq_response <- seq(
  stats::quantile(response_model_df$SPEI_3, 0.02),
  stats::quantile(response_model_df$SPEI_3, 0.98),
  length.out = 250
)
response_prediction_grid <- expand.grid(
  water_class = ecosystem_classes,
  SPEI_3 = spei_seq_response,
  month_f = factor("7", levels = levels(response_model_df$month_f)),
  site_name = levels(response_model_df$site_name)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
response_prediction_grid$water_class <- factor(
  response_prediction_grid$water_class,
  levels = ecosystem_classes
)
response_prediction_grid$site_name <- factor(
  response_prediction_grid$site_name,
  levels = levels(response_model_df$site_name)
)

response_grid_pred <- predict(
  response_model,
  newdata = response_prediction_grid,
  se.fit = TRUE,
  exclude = "s(site_name)"
)
response_prediction_grid$predicted_pct_change <- as.numeric(response_grid_pred$fit)
response_prediction_grid$predicted_pct_change_lower <-
  as.numeric(response_grid_pred$fit - 1.96 * response_grid_pred$se.fit)
response_prediction_grid$predicted_pct_change_upper <-
  as.numeric(response_grid_pred$fit + 1.96 * response_grid_pred$se.fit)
response_prediction_grid$predicted_magnitude <- abs(response_prediction_grid$predicted_pct_change)
response_prediction_grid$response_class <- classify_response_magnitude(
  response_prediction_grid$predicted_magnitude
)
response_prediction_grid$response_direction <- classify_response_direction(
  response_prediction_grid$predicted_pct_change
)

response_thresholds <- do.call(
  rbind,
  lapply(split(response_prediction_grid, response_prediction_grid$water_class), function(group_data) {
    do.call(
      rbind,
      lapply(c(0.5, 1, 2.5, 5), function(threshold) {
        stable_values <- group_data$SPEI_3[group_data$predicted_magnitude < threshold]
        data.frame(
          water_class = group_data$water_class[1],
          response_threshold_pct = threshold,
          stable_SPEI_min = ifelse(length(stable_values) > 0, min(stable_values), NA_real_),
          stable_SPEI_max = ifelse(length(stable_values) > 0, max(stable_values), NA_real_),
          dry_any_impact_threshold = find_first_threshold(group_data, threshold, "dry", "any"),
          dry_decrease_threshold = find_first_threshold(group_data, threshold, "dry", "decrease"),
          dry_increase_threshold = find_first_threshold(group_data, threshold, "dry", "increase"),
          wet_any_impact_threshold = find_first_threshold(group_data, threshold, "wet", "any"),
          wet_decrease_threshold = find_first_threshold(group_data, threshold, "wet", "decrease"),
          wet_increase_threshold = find_first_threshold(group_data, threshold, "wet", "increase"),
          min_predicted_pct_change = min(group_data$predicted_pct_change),
          max_predicted_pct_change = max(group_data$predicted_pct_change),
          max_predicted_magnitude = max(group_data$predicted_magnitude)
        )
      })
    )
  })
)

response_smooth_table <- as.data.frame(summary(response_model)$s.table)
response_smooth_table$smooth_term <- rownames(response_smooth_table)
rownames(response_smooth_table) <- NULL
response_smooth_table <- response_smooth_table[, c("smooth_term", setdiff(names(response_smooth_table), "smooth_term"))]

response_parametric_table <- as.data.frame(summary(response_model)$p.table)
response_parametric_table$term <- rownames(response_parametric_table)
rownames(response_parametric_table) <- NULL
response_parametric_table <- response_parametric_table[, c("term", setdiff(names(response_parametric_table), "term"))]

# ---------------------------------------------------------------------------
# Step 4b — Q3-aligned ecological impact model
#
# This is the primary EDI revision. It follows the final Q3 logic directly:
# ecological impact is the predicted magnitude and direction of WUE_T change
# from near-normal SPEI conditions. The manuscript-facing mapped/threshold
# result is upland-only and coast-specific.
# ---------------------------------------------------------------------------

spei_cols <- c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48")
q3_required_cols <- c("site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra", spei_cols)
q3_missing_cols <- setdiff(q3_required_cols, names(monthly))
if (length(q3_missing_cols) > 0) {
  stop("Missing required Q3-aligned EDI columns: ", paste(q3_missing_cols, collapse = ", "))
}

coast_region_levels <- c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast")
q3_base <- monthly[complete.cases(monthly[, q3_required_cols]), ]
q3_base <- q3_base[
  q3_base$site_name != "" &
    q3_base$water_class != "" &
    q3_base$water_class %in% ecosystem_classes &
    is.finite(q3_base$lat) &
    is.finite(q3_base$long) &
    is.finite(q3_base$WUE_tra),
]
q3_base$coast_region <- factor(
  classify_coast_region(q3_base$lat, q3_base$long),
  levels = coast_region_levels
)
q3_base$site_name <- factor(q3_base$site_name)
q3_base$water_class <- factor(q3_base$water_class, levels = ecosystem_classes)
q3_base$month_f <- factor(q3_base$month, levels = 1:12)
q3_base$WUE_T <- q3_base$WUE_tra

q3_spei_long <- reshape(
  q3_base,
  varying = spei_cols,
  v.names = "SPEI_value",
  timevar = "SPEI_timescale",
  times = spei_cols,
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)
q3_spei_long <- q3_spei_long[is.finite(q3_spei_long$SPEI_value) & is.finite(q3_spei_long$WUE_T), ]
q3_spei_long$SPEI_timescale <- factor(q3_spei_long$SPEI_timescale, levels = spei_cols)
q3_spei_long$site_name <- factor(q3_spei_long$site_name)
q3_spei_long$water_class <- factor(q3_spei_long$water_class, levels = ecosystem_classes)
q3_spei_long$coast_region <- factor(q3_spei_long$coast_region, levels = coast_region_levels)
q3_spei_long$month_f <- factor(q3_spei_long$month, levels = 1:12)
q3_spei_long$spei_coast <- interaction(
  q3_spei_long$SPEI_timescale,
  q3_spei_long$coast_region,
  sep = "__",
  drop = TRUE
)

q3_coast_model <- mgcv::gam(
  WUE_T ~
    SPEI_timescale * coast_region +
    water_class +
    month_f +
    s(SPEI_value, by = spei_coast, k = 6) +
    s(site_name, bs = "re"),
  data = q3_spei_long,
  method = "ML",
  select = TRUE
)
saveRDS(q3_coast_model, file.path(output_dir, "EDI_Q3_aligned_coast_threshold_gam_model.rds"))

q3_model_info <- data.frame(
  model = "Q3_aligned_coast_threshold_gam",
  response = "WUE_T",
  primary_ecosystem = "Upland",
  AIC = AIC(q3_coast_model),
  BIC = BIC(q3_coast_model),
  logLik = as.numeric(stats::logLik(q3_coast_model)),
  deviance_explained = summary(q3_coast_model)$dev.expl,
  adjusted_r_squared = summary(q3_coast_model)$r.sq,
  n_rows_long = nrow(q3_spei_long),
  n_sites = length(unique(q3_spei_long$site_name))
)

q3_coast_smooth_table <- as.data.frame(summary(q3_coast_model)$s.table)
q3_coast_smooth_table$smooth_term <- rownames(q3_coast_smooth_table)
rownames(q3_coast_smooth_table) <- NULL
q3_coast_smooth_table <- q3_coast_smooth_table[, c("smooth_term", setdiff(names(q3_coast_smooth_table), "smooth_term"))]

q3_coast_parametric_table <- as.data.frame(summary(q3_coast_model)$p.table)
q3_coast_parametric_table$term <- rownames(q3_coast_parametric_table)
rownames(q3_coast_parametric_table) <- NULL
q3_coast_parametric_table <- q3_coast_parametric_table[, c("term", setdiff(names(q3_coast_parametric_table), "term"))]

q3_spei_sequence <- seq(
  stats::quantile(q3_spei_long$SPEI_value, 0.02),
  stats::quantile(q3_spei_long$SPEI_value, 0.98),
  length.out = 120
)
q3_coast_prediction_grid <- expand.grid(
  SPEI_timescale = levels(q3_spei_long$SPEI_timescale),
  coast_region = levels(q3_spei_long$coast_region),
  SPEI_value = q3_spei_sequence,
  water_class = factor("Upland", levels = levels(q3_spei_long$water_class)),
  month_f = factor("7", levels = levels(q3_spei_long$month_f)),
  site_name = levels(q3_spei_long$site_name)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
q3_coast_prediction_grid$SPEI_timescale <- factor(
  q3_coast_prediction_grid$SPEI_timescale,
  levels = levels(q3_spei_long$SPEI_timescale)
)
q3_coast_prediction_grid$coast_region <- factor(
  q3_coast_prediction_grid$coast_region,
  levels = levels(q3_spei_long$coast_region)
)
q3_coast_prediction_grid$water_class <- factor(
  q3_coast_prediction_grid$water_class,
  levels = levels(q3_spei_long$water_class)
)
q3_coast_prediction_grid$site_name <- factor(
  q3_coast_prediction_grid$site_name,
  levels = levels(q3_spei_long$site_name)
)
q3_coast_prediction_grid$spei_coast <- interaction(
  q3_coast_prediction_grid$SPEI_timescale,
  q3_coast_prediction_grid$coast_region,
  sep = "__",
  drop = TRUE
)
q3_coast_prediction_grid$spei_coast <- factor(
  q3_coast_prediction_grid$spei_coast,
  levels = levels(q3_spei_long$spei_coast)
)

q3_coast_pred <- predict(
  q3_coast_model,
  newdata = q3_coast_prediction_grid,
  type = "link",
  se.fit = TRUE,
  exclude = "s(site_name)"
)
q3_coast_prediction_grid$predicted_WUE_T <- as.numeric(q3_coast_pred$fit)
q3_coast_prediction_grid$predicted_se <- as.numeric(q3_coast_pred$se.fit)
q3_coast_prediction_grid$predicted_lower <- q3_coast_prediction_grid$predicted_WUE_T -
  1.96 * q3_coast_prediction_grid$predicted_se
q3_coast_prediction_grid$predicted_upper <- q3_coast_prediction_grid$predicted_WUE_T +
  1.96 * q3_coast_prediction_grid$predicted_se

q3_coast_near_normal <- q3_coast_prediction_grid[
  q3_coast_prediction_grid$SPEI_value >= -1 & q3_coast_prediction_grid$SPEI_value <= 1,
] |>
  dplyr::group_by(SPEI_timescale, coast_region) |>
  dplyr::summarise(
    predicted_near_normal_WUE_T = mean(predicted_WUE_T),
    .groups = "drop"
  )

q3_coast_prediction_impact <- dplyr::left_join(
  q3_coast_prediction_grid,
  q3_coast_near_normal,
  by = c("SPEI_timescale", "coast_region")
)
q3_coast_prediction_impact$predicted_change <- q3_coast_prediction_impact$predicted_WUE_T -
  q3_coast_prediction_impact$predicted_near_normal_WUE_T
q3_coast_prediction_impact$predicted_pct_change <- 100 * q3_coast_prediction_impact$predicted_change /
  q3_coast_prediction_impact$predicted_near_normal_WUE_T
q3_coast_prediction_impact$predicted_lower_pct <- 100 *
  (q3_coast_prediction_impact$predicted_lower - q3_coast_prediction_impact$predicted_near_normal_WUE_T) /
  q3_coast_prediction_impact$predicted_near_normal_WUE_T
q3_coast_prediction_impact$predicted_upper_pct <- 100 *
  (q3_coast_prediction_impact$predicted_upper - q3_coast_prediction_impact$predicted_near_normal_WUE_T) /
  q3_coast_prediction_impact$predicted_near_normal_WUE_T
q3_coast_prediction_impact$EDI_Q3_magnitude <- abs(q3_coast_prediction_impact$predicted_pct_change)
q3_coast_prediction_impact$EDI_Q3_class <- classify_q3_response_magnitude(
  q3_coast_prediction_impact$EDI_Q3_magnitude
)
q3_coast_prediction_impact$EDI_Q3_direction <- classify_q3_response_direction(
  q3_coast_prediction_impact$predicted_pct_change
)
q3_coast_prediction_impact$prediction_interval_support <- factor(
  ifelse(
    q3_coast_prediction_impact$predicted_upper_pct <= -5, "PI supports decrease",
    ifelse(q3_coast_prediction_impact$predicted_lower_pct >= 5, "PI supports increase", "PI overlaps no-change")
  ),
  levels = c("PI supports decrease", "PI overlaps no-change", "PI supports increase")
)

q3_coast_threshold_markers <- do.call(
  rbind,
  lapply(
    split(q3_coast_prediction_impact, list(q3_coast_prediction_impact$SPEI_timescale, q3_coast_prediction_impact$coast_region), drop = TRUE),
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
            pct_change_threshold = c(-threshold_pct, threshold_pct, -threshold_pct, threshold_pct)
          )
        })
      )
    }
  )
)
q3_coast_threshold_markers <- q3_coast_threshold_markers[!is.na(q3_coast_threshold_markers$SPEI_threshold), ]
q3_coast_threshold_markers$SPEI_threshold_lower <- NA_real_
q3_coast_threshold_markers$SPEI_threshold_upper <- NA_real_
for (i in seq_len(nrow(q3_coast_threshold_markers))) {
  marker_row <- q3_coast_threshold_markers[i, ]
  marker_group <- q3_coast_prediction_impact[
    q3_coast_prediction_impact$SPEI_timescale == marker_row$SPEI_timescale &
      q3_coast_prediction_impact$coast_region == marker_row$coast_region,
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
  q3_coast_threshold_markers$SPEI_threshold_lower[i] <- min(threshold_range)
  q3_coast_threshold_markers$SPEI_threshold_upper[i] <- max(threshold_range)
}
q3_coast_threshold_markers$threshold_pct <- factor(
  q3_coast_threshold_markers$threshold_pct,
  levels = c(5, 10, 20),
  labels = c("5%", "10%", "20%")
)
q3_coast_threshold_markers$SPEI_timescale <- factor(
  q3_coast_threshold_markers$SPEI_timescale,
  levels = levels(q3_spei_long$SPEI_timescale)
)
q3_coast_threshold_markers$coast_region <- factor(
  q3_coast_threshold_markers$coast_region,
  levels = levels(q3_spei_long$coast_region)
)
q3_coast_threshold_summary <- q3_coast_threshold_markers |>
  dplyr::arrange(SPEI_timescale, coast_region, anomaly_side, impact_direction, threshold_pct)

q3_monthly_scores <- q3_base[q3_base$water_class == "Upland", ]
q3_monthly_scores$SPEI_timescale <- factor("SPEI_3", levels = levels(q3_spei_long$SPEI_timescale))
q3_monthly_scores$SPEI_value <- q3_monthly_scores$SPEI_3
q3_monthly_scores$spei_coast <- interaction(
  q3_monthly_scores$SPEI_timescale,
  q3_monthly_scores$coast_region,
  sep = "__",
  drop = TRUE
)
q3_monthly_scores$spei_coast <- factor(q3_monthly_scores$spei_coast, levels = levels(q3_spei_long$spei_coast))
q3_monthly_scores$site_name <- factor(q3_monthly_scores$site_name, levels = levels(q3_spei_long$site_name))
q3_monthly_scores$month_f <- factor(q3_monthly_scores$month, levels = levels(q3_spei_long$month_f))
q3_monthly_scores$water_class <- factor(q3_monthly_scores$water_class, levels = levels(q3_spei_long$water_class))
q3_monthly_pred <- predict(
  q3_coast_model,
  newdata = q3_monthly_scores,
  type = "link",
  se.fit = TRUE,
  exclude = "s(site_name)"
)
q3_monthly_scores$predicted_WUE_T <- as.numeric(q3_monthly_pred$fit)
q3_monthly_scores$predicted_se <- as.numeric(q3_monthly_pred$se.fit)
q3_monthly_scores$predicted_lower <- q3_monthly_scores$predicted_WUE_T - 1.96 * q3_monthly_scores$predicted_se
q3_monthly_scores$predicted_upper <- q3_monthly_scores$predicted_WUE_T + 1.96 * q3_monthly_scores$predicted_se
q3_monthly_scores <- dplyr::left_join(
  q3_monthly_scores,
  q3_coast_near_normal[q3_coast_near_normal$SPEI_timescale == "SPEI_3", ],
  by = c("SPEI_timescale", "coast_region")
)
q3_monthly_scores$EDI_Q3_pct_change <- 100 *
  (q3_monthly_scores$predicted_WUE_T - q3_monthly_scores$predicted_near_normal_WUE_T) /
  q3_monthly_scores$predicted_near_normal_WUE_T
q3_monthly_scores$EDI_Q3_pct_change_lower <- 100 *
  (q3_monthly_scores$predicted_lower - q3_monthly_scores$predicted_near_normal_WUE_T) /
  q3_monthly_scores$predicted_near_normal_WUE_T
q3_monthly_scores$EDI_Q3_pct_change_upper <- 100 *
  (q3_monthly_scores$predicted_upper - q3_monthly_scores$predicted_near_normal_WUE_T) /
  q3_monthly_scores$predicted_near_normal_WUE_T
q3_monthly_scores$EDI_Q3_magnitude <- abs(q3_monthly_scores$EDI_Q3_pct_change)
q3_monthly_scores$EDI_Q3_class <- classify_q3_response_magnitude(q3_monthly_scores$EDI_Q3_magnitude)
q3_monthly_scores$EDI_Q3_direction <- classify_q3_response_direction(q3_monthly_scores$EDI_Q3_pct_change)
q3_monthly_scores$prediction_interval_support <- factor(
  ifelse(
    q3_monthly_scores$EDI_Q3_pct_change_upper <= -5, "PI supports decrease",
    ifelse(q3_monthly_scores$EDI_Q3_pct_change_lower >= 5, "PI supports increase", "PI overlaps no-change")
  ),
  levels = c("PI supports decrease", "PI overlaps no-change", "PI supports increase")
)

q3_monthly_out <- q3_monthly_scores[, c(
  "site_name", "Year", "month", "water_class", "coast_region", "lat", "long",
  "WUE_T", "SPEI_3", "predicted_near_normal_WUE_T", "predicted_WUE_T",
  "EDI_Q3_pct_change", "EDI_Q3_pct_change_lower", "EDI_Q3_pct_change_upper",
  "EDI_Q3_magnitude", "EDI_Q3_direction", "EDI_Q3_class", "prediction_interval_support"
)]

q3_monthly_regional_summary <- q3_monthly_out |>
  dplyr::group_by(coast_region, EDI_Q3_class, EDI_Q3_direction, prediction_interval_support) |>
  dplyr::summarise(
    n_site_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_SPEI_3 = mean(SPEI_3, na.rm = TRUE),
    mean_pct_change = mean(EDI_Q3_pct_change, na.rm = TRUE),
    median_abs_pct_change = median(EDI_Q3_magnitude, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(coast_region) |>
  dplyr::mutate(prop_site_months = n_site_months / sum(n_site_months) * 100) |>
  dplyr::ungroup()

q3_regional_validation <- q3_monthly_out |>
  dplyr::group_by(coast_region) |>
  dplyr::summarise(
    n_site_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_abs_Q3_EDI_pct_change = mean(EDI_Q3_magnitude, na.rm = TRUE),
    median_abs_Q3_EDI_pct_change = median(EDI_Q3_magnitude, na.rm = TRUE),
    pct_months_5pct_or_more = mean(EDI_Q3_magnitude >= 5, na.rm = TRUE) * 100,
    pct_months_10pct_or_more = mean(EDI_Q3_magnitude >= 10, na.rm = TRUE) * 100,
    pct_months_PI_supported = mean(prediction_interval_support != "PI overlaps no-change", na.rm = TRUE) * 100,
    .groups = "drop"
  )
q3_site_sensitivity_file <- file.path(
  script_dir,
  "results",
  "Q3_WUE_T_SPEI_sensitivity",
  "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv"
)
if (file.exists(q3_site_sensitivity_file)) {
  q3_site_sensitivity_summary <- read.csv(q3_site_sensitivity_file, stringsAsFactors = FALSE)
  q3_site_sensitivity_summary <- q3_site_sensitivity_summary[
    q3_site_sensitivity_summary$SPEI_timescale == "SPEI_3",
  ]
  q3_regional_validation <- dplyr::left_join(
    q3_regional_validation,
    q3_site_sensitivity_summary[, c(
      "coast_region", "n_sites", "median_abs_pct_change",
      "q25_abs_pct_change", "q75_abs_pct_change", "mean_abs_pct_change"
    )],
    by = "coast_region",
    suffix = c("_EDI_monthly", "_Q3_site")
  )
}

q3_prediction_source_file <- file.path(
  q3_dir,
  "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv"
)
q3_threshold_source_file <- file.path(
  q3_dir,
  "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_summary_5_10_20pct_upland.csv"
)
q3_smooth_source_file <- file.path(
  q3_dir,
  "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv"
)
q3_parametric_source_file <- file.path(
  q3_dir,
  "Q3_WUE_T_SPEI_coast_threshold_GAM_parametric_terms.csv"
)
q3_model_comparison_source_file <- file.path(
  q3_dir,
  "Q3_WUE_T_SPEI_model_comparison.csv"
)

missing_q3_source_files <- c(
  q3_prediction_source_file,
  q3_threshold_source_file
)[!file.exists(c(q3_prediction_source_file, q3_threshold_source_file))]
if (length(missing_q3_source_files) > 0) {
  stop(
    "Missing final Q3 source output(s): ",
    paste(missing_q3_source_files, collapse = ", "),
    ". Run Q3.WUE_Sensitivity_SPEI.R before Ecological_Impacts.R."
  )
}

# The manuscript-facing EDI response surface now uses the final Q3 outputs
# directly. The model fit above is retained for site-level diagnostics, but the
# spatial workflow should consume these Q3-derived prediction and threshold
# tables rather than an independently refit response surface.
q3_coast_prediction_impact <- read.csv(q3_prediction_source_file, stringsAsFactors = FALSE)
q3_coast_prediction_impact$SPEI_timescale <- factor(
  q3_coast_prediction_impact$SPEI_timescale,
  levels = spei_cols
)
q3_coast_prediction_impact$coast_region <- factor(
  q3_coast_prediction_impact$coast_region,
  levels = coast_region_levels
)
q3_coast_prediction_impact$water_class <- factor(
  q3_coast_prediction_impact$water_class,
  levels = ecosystem_classes
)
q3_coast_prediction_impact$EDI_Q3_magnitude <- abs(q3_coast_prediction_impact$predicted_pct_change)
q3_coast_prediction_impact$EDI_Q3_class <- classify_q3_response_magnitude(
  q3_coast_prediction_impact$EDI_Q3_magnitude
)
q3_coast_prediction_impact$EDI_Q3_direction <- classify_q3_response_direction(
  q3_coast_prediction_impact$predicted_pct_change
)
q3_coast_prediction_impact$prediction_interval_support <- factor(
  ifelse(
    q3_coast_prediction_impact$predicted_upper_pct <= -5, "PI supports decrease",
    ifelse(q3_coast_prediction_impact$predicted_lower_pct >= 5, "PI supports increase", "PI overlaps no-change")
  ),
  levels = c("PI supports decrease", "PI overlaps no-change", "PI supports increase")
)

q3_coast_threshold_summary <- read.csv(q3_threshold_source_file, stringsAsFactors = FALSE)
q3_coast_threshold_summary$SPEI_timescale <- factor(
  q3_coast_threshold_summary$SPEI_timescale,
  levels = spei_cols
)
q3_coast_threshold_summary$coast_region <- factor(
  q3_coast_threshold_summary$coast_region,
  levels = coast_region_levels
)
q3_coast_threshold_summary$water_class <- factor(
  q3_coast_threshold_summary$water_class,
  levels = ecosystem_classes
)
q3_coast_threshold_summary <- q3_coast_threshold_summary |>
  dplyr::arrange(SPEI_timescale, coast_region, anomaly_side, impact_direction, threshold_pct)

if (file.exists(q3_smooth_source_file)) {
  q3_coast_smooth_table <- read.csv(q3_smooth_source_file, stringsAsFactors = FALSE)
}
if (file.exists(q3_parametric_source_file)) {
  q3_coast_parametric_table <- read.csv(q3_parametric_source_file, stringsAsFactors = FALSE)
}
if (file.exists(q3_model_comparison_source_file)) {
  q3_model_comparison <- read.csv(q3_model_comparison_source_file, stringsAsFactors = FALSE)
  q3_model_info <- q3_model_comparison[
    q3_model_comparison$model == "coast_threshold_gam",
  ]
}

write_table(q3_model_info, "EDI_Q3_aligned_model_info.csv")
write_table(q3_coast_smooth_table, "EDI_Q3_aligned_coast_GAM_smooth_terms.csv")
write_table(q3_coast_parametric_table, "EDI_Q3_aligned_coast_GAM_parametric_terms.csv")
write_table(q3_coast_prediction_impact, "EDI_Q3_aligned_prediction_scores_upland.csv")
write_table(q3_coast_threshold_summary, "EDI_Q3_aligned_coast_thresholds_5_10_20pct.csv")
write_table(q3_monthly_out, "EDI_Q3_aligned_monthly_scores.csv")
write_table(q3_monthly_regional_summary, "EDI_Q3_aligned_monthly_regional_summary.csv")
write_table(q3_regional_validation, "EDI_Q3_aligned_regional_validation.csv")

# ---------------------------------------------------------------------------
# Step 5 — SPEI-3 thresholds corresponding to logistic EDI class boundaries
# ---------------------------------------------------------------------------

b0 <- coef(pooled_model)[["(Intercept)"]]
b1 <- coef(pooled_model)[["SPEI_3"]]

thresholds_pooled <- data.frame(
  model        = "Pooled",
  water_class  = "All",
  EDI_0.30_SPEI3 = spei_at_p(0.30, b0, b1),
  EDI_0.50_SPEI3 = spei_at_p(0.50, b0, b1),
  EDI_0.70_SPEI3 = spei_at_p(0.70, b0, b1)
)

thresholds_eco <- do.call(rbind, lapply(ecosystem_classes, function(ec) {
  m <- eco_models[[ec]]
  if (is.null(m)) return(NULL)
  cc <- coef(m)
  data.frame(
    model          = "Ecosystem",
    water_class    = ec,
    EDI_0.30_SPEI3 = spei_at_p(0.30, cc[1], cc[2]),
    EDI_0.50_SPEI3 = spei_at_p(0.50, cc[1], cc[2]),
    EDI_0.70_SPEI3 = spei_at_p(0.70, cc[1], cc[2])
  )
}))

thresholds <- rbind(thresholds_pooled, thresholds_eco)

# ---------------------------------------------------------------------------
# Step 5 — Validate: do EDI classes align with Q2 resistance & recovery?
# ---------------------------------------------------------------------------

site_metrics <- dplyr::left_join(
  resistance[, c("site_name", "resistance")],
  recovery[, c("site_name", "mean_recovery", "pct_full_recovery")],
  by = "site_name"
)

# For drought months only, join site-level Q2 metrics and summarise by EDI class
validation_df <- base_df |>
  dplyr::filter(SPEI_3 < -1) |>
  dplyr::left_join(site_metrics, by = "site_name") |>
  dplyr::filter(!is.na(resistance)) |>
  dplyr::group_by(EDI_class, water_class) |>
  dplyr::summarise(
    n_months              = dplyr::n(),
    n_sites               = dplyr::n_distinct(site_name),
    mean_SPEI_3           = mean(SPEI_3),
    mean_EDI              = mean(EDI),
    mean_WUE_T            = mean(WUE_T),
    mean_monthly_resist   = mean(monthly_resistance, na.rm = TRUE),
    mean_site_resistance  = mean(resistance,    na.rm = TRUE),
    mean_recovery         = mean(mean_recovery, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# Coefficient table and output
# ---------------------------------------------------------------------------

coef_table <- rbind(
  data.frame(
    water_class  = "All (pooled)",
    intercept    = b0,
    slope_SPEI3  = b1,
    AIC          = AIC(pooled_model),
    n_obs        = nrow(model_df),
    n_impacted   = sum(model_df$y)
  ),
  do.call(rbind, lapply(ecosystem_classes, function(ec) {
    m <- eco_models[[ec]]
    if (is.null(m)) return(NULL)
    cc <- coef(m)
    sub <- model_df[model_df$water_class == ec, ]
    data.frame(
      water_class  = ec,
      intercept    = cc[1],
      slope_SPEI3  = cc[2],
      AIC          = AIC(m),
      n_obs        = nrow(sub),
      n_impacted   = sum(sub$y)
    )
  }))
)

monthly_out <- base_df[, c("site_name", "Year", "month", "water_class",
                            "WUE_T", "SPEI_3", "baseline_WUE_T",
                            "monthly_resistance", "WUE_T_pct_change",
                            "EDI", "EDI_class",
                            "response_EDI_pct_change",
                            "response_EDI_magnitude",
                            "response_EDI_direction",
                            "response_EDI_class")]

write_table(coef_table,   "EDI_logistic_coefficients.csv")
write_table(thresholds,   "EDI_threshold_summary.csv")
write_table(monthly_out,  "EDI_monthly_site_scores.csv")
write_table(validation_df,"EDI_validation_by_class.csv")
write_table(response_prediction_grid, "EDI_response_prediction_curves.csv")
write_table(response_thresholds, "EDI_response_threshold_summary.csv")
write_table(response_smooth_table, "EDI_response_gam_smooth_terms.csv")
write_table(response_parametric_table, "EDI_response_gam_parametric_terms.csv")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

ecosystem_colors <- c(Upland = "#800080", Freshwater = "#0000FF", Saline = "#FFA500")
coast_colors <- c(
  "Atlantic Coast" = "#008B8B",
  "Pacific Coast" = "#C44E52",
  "Gulf Coast" = "#8C564B",
  "AK Coast" = "#4D4D4D"
)
q3_class_colors <- c(
  "No meaningful change (<5%)" = "#F2F2F2",
  "Watch (5-10%)" = "#BFD3E6",
  "Stress (10-20%)" = "#F4A582",
  "Impact (>=20%)" = "#B2182B"
)
q3_pi_colors <- c(
  "PI supports decrease" = "#2166AC",
  "PI overlaps no-change" = "#D9D9D9",
  "PI supports increase" = "#B2182B"
)
severity_colors  <- c(None   = "#2166AC", Watch  = "#92C5DE",
                      Stress = "#F4A582", Impact = "#D6604D")
line_colors      <- c("All (pooled)" = "black",
                      Upland = "#800080", Freshwater = "#0000FF", Saline = "#FFA500")

spei_seq <- seq(-4, 2, length.out = 300)

# -- Plot 8: Q3-aligned EDI composite ---------------------------------------
q3_selected_timescales <- c("SPEI_1", "SPEI_3", "SPEI_48")
q3_coast_panel_labels <- data.frame(
  coast_region = factor(coast_region_levels, levels = coast_region_levels),
  SPEI_timescale = factor("SPEI_1", levels = levels(q3_spei_long$SPEI_timescale)),
  coast_label = c("Atlantic", "Pacific", "Gulf", "AK")
)

p_q3_threshold <- ggplot2::ggplot(
  q3_coast_prediction_impact[q3_coast_prediction_impact$SPEI_timescale %in% q3_selected_timescales, ],
  ggplot2::aes(x = SPEI_value, y = predicted_pct_change, color = coast_region)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45") +
  ggplot2::geom_hline(yintercept = 0, color = "gray35") +
  ggplot2::geom_hline(yintercept = c(-20, -10, -5, 5, 10, 20),
                      linetype = "dotted", color = "gray45", linewidth = 0.28) +
  ggplot2::geom_label(
    data = q3_coast_panel_labels,
    ggplot2::aes(label = coast_label),
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.18,
    size = 4.5,
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
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_segment(
    data = q3_coast_threshold_markers[
      q3_coast_threshold_markers$SPEI_timescale %in% q3_selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold_lower,
      xend = SPEI_threshold_upper,
      y = pct_change_threshold,
      yend = pct_change_threshold,
      color = coast_region
    ),
    linewidth = 1.05,
    alpha = 0.72,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = q3_coast_threshold_markers[
      q3_coast_threshold_markers$SPEI_timescale %in% q3_selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      xend = SPEI_threshold,
      y = 0,
      yend = pct_change_threshold,
      color = coast_region,
      linetype = threshold_pct
    ),
    linewidth = 0.72,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = q3_coast_threshold_markers[
      q3_coast_threshold_markers$SPEI_timescale %in% q3_selected_timescales,
    ],
    ggplot2::aes(
      x = SPEI_threshold,
      y = pct_change_threshold,
      color = coast_region,
      shape = threshold_pct
    ),
    size = 3.1,
    stroke = 1.0,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_grid(coast_region ~ SPEI_timescale, scales = "free") +
  ggplot2::scale_color_manual(values = coast_colors) +
  ggplot2::scale_fill_manual(values = coast_colors) +
  ggplot2::scale_shape_manual(values = c("5%" = 21, "10%" = 24, "20%" = 22)) +
  ggplot2::guides(color = "none", fill = "none") +
  ggplot2::labs(
    title = expression(bold("Q3-Aligned Coast-Specific Upland " * WUE[T] * " Thresholds")),
    x = "SPEI",
    y = expression("Upland " * WUE[T] * " change (%)"),
    shape = "Change threshold",
    linetype = "Change threshold"
  ) +
  theme_wue(base_size = 13) +
  ggplot2::theme(
    strip.text.y = ggplot2::element_blank(),
    strip.background.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5, 18, 5, 5)
  )

q3_class_props <- q3_monthly_out |>
  dplyr::group_by(coast_region, EDI_Q3_class) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(coast_region) |>
  dplyr::mutate(prop = n / sum(n) * 100) |>
  dplyr::ungroup()

p_q3_class_props <- ggplot2::ggplot(
  q3_class_props,
  ggplot2::aes(x = coast_region, y = prop, fill = EDI_Q3_class)
) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = q3_class_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  ggplot2::labs(
    title = "Monthly SPEI-3 Upland EDI Classes",
    x = NULL,
    y = "% of upland site-months",
    fill = "Q3 EDI class"
  ) +
  theme_wue(base_size = 13) +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.size = grid::unit(0.6, "lines")
  )

q3_pi_props <- q3_monthly_out |>
  dplyr::group_by(coast_region, prediction_interval_support) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(coast_region) |>
  dplyr::mutate(prop = n / sum(n) * 100) |>
  dplyr::ungroup()

p_q3_pi_props <- ggplot2::ggplot(
  q3_pi_props,
  ggplot2::aes(x = coast_region, y = prop, fill = prediction_interval_support)
) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = q3_pi_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  ggplot2::labs(
    title = "Prediction-Interval Support",
    x = NULL,
    y = "% of upland site-months",
    fill = "95% interval"
  ) +
  theme_wue(base_size = 13) +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.size = grid::unit(0.6, "lines")
  )

q3_composite <- p_q3_threshold / (p_q3_class_props | p_q3_pi_props) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = expression("Q3-Aligned Ecological Drought Index for Upland " * WUE[T]),
    subtitle = "Primary EDI is predicted WUE_T change from near-normal SPEI; older logistic EDI outputs are retained as diagnostics.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 10))
    )
  ) &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 18)
  )

save_plot(q3_composite, "EDI_plot_08_Q3_aligned_composite.png", width = 11, height = 15)

# Prediction curves for each model
pred_curves <- do.call(rbind, c(
  list(data.frame(SPEI_3 = spei_seq,
                  EDI    = logistic(spei_seq, b0, b1),
                  model  = "All (pooled)")),
  lapply(ecosystem_classes, function(ec) {
    m <- eco_models[[ec]]
    if (is.null(m)) return(NULL)
    cc <- coef(m)
    data.frame(SPEI_3 = spei_seq,
               EDI    = logistic(spei_seq, cc[1], cc[2]),
               model  = ec)
  })
))
pred_curves$model <- factor(pred_curves$model,
                             levels = c("All (pooled)", ecosystem_classes))

threshold_lines <- data.frame(y = c(0.30, 0.50, 0.70),
                               label = c("Watch", "Stress", "Impact"))

# -- Plot 0: response-based EDI curve ---------------------------------------
p_response <- ggplot2::ggplot(
  response_prediction_grid,
  ggplot2::aes(x = SPEI_3, y = predicted_pct_change, color = water_class, fill = water_class)
) +
  ggplot2::annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf,
                    fill = "gray85", alpha = 0.25) +
  ggplot2::geom_hline(yintercept = 0, color = "gray25", linewidth = 0.55) +
  ggplot2::geom_hline(yintercept = c(-5, -2.5, 2.5, 5),
                      linetype = "dotted", color = "gray70", linewidth = 0.28) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45", linewidth = 0.35) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = predicted_pct_change_lower, ymax = predicted_pct_change_upper),
    alpha = 0.10,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::facet_wrap(~ water_class, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  ggplot2::labs(
    title = "Response-Based EDI: Expected WUE_T Change Across SPEI-3",
    subtitle = "Gray band marks near-normal SPEI-3. Dotted lines mark +/-2.5% and +/-5% response thresholds.",
    x = "SPEI-3",
    y = "Predicted WUE_T change (%)",
    color = "Ecosystem class",
    fill = "Ecosystem class"
  ) +
  theme_wue(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 6)),
    plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 12)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 8)),
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 8)),
    strip.text = ggplot2::element_text(size = 12, face = "bold"),
    panel.spacing.y = grid::unit(1.0, "lines")
  )

save_plot(p_response, "EDI_plot_00_response_curve.png", width = 10, height = 12)

# -- Plot 0b: response class proportions ------------------------------------
response_props <- base_df |>
  dplyr::filter(!is.na(response_EDI_class)) |>
  dplyr::group_by(water_class, response_EDI_class, response_EDI_direction) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  dplyr::group_by(water_class) |>
  dplyr::mutate(prop = n / sum(n) * 100) |>
  dplyr::ungroup()

p_response_props <- ggplot2::ggplot(
  response_props,
  ggplot2::aes(x = water_class, y = prop, fill = response_EDI_class)
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  ggplot2::labs(
    title = "Response-Based EDI Class Proportions",
    subtitle = "Classes are based on absolute expected WUE_T change from near-normal",
    x = "Ecosystem class",
    y = "% of site-months",
    fill = "Response class"
  ) +
  theme_wue() +
  ggplot2::theme(legend.position = "right")

save_plot(p_response_props, "EDI_plot_00b_response_class_proportions.png", width = 10, height = 7)

# -- Plot 1: logistic curve --------------------------------------------------
p_logistic <- ggplot2::ggplot(pred_curves,
    ggplot2::aes(x = SPEI_3, y = EDI, color = model, linetype = model)) +
  ggplot2::annotate("rect", xmin=-4, xmax=2, ymin=0,    ymax=0.30,
                    fill=severity_colors["None"],   alpha=0.10) +
  ggplot2::annotate("rect", xmin=-4, xmax=2, ymin=0.30, ymax=0.50,
                    fill=severity_colors["Watch"],  alpha=0.10) +
  ggplot2::annotate("rect", xmin=-4, xmax=2, ymin=0.50, ymax=0.70,
                    fill=severity_colors["Stress"], alpha=0.10) +
  ggplot2::annotate("rect", xmin=-4, xmax=2, ymin=0.70, ymax=1.00,
                    fill=severity_colors["Impact"], alpha=0.10) +
  ggplot2::geom_hline(data=threshold_lines, ggplot2::aes(yintercept=y),
                      linetype="dashed", color="gray40", linewidth=0.45) +
  ggplot2::geom_vline(xintercept=-1, linetype="dotted", color="gray40") +
  ggplot2::geom_line(linewidth=1.0) +
  ggplot2::annotate("text", x=-3.8, y=c(0.15,0.40,0.60,0.85),
                    label=c("None","Watch","Stress","Impact"),
                    hjust=0, size=3.5, color="gray30") +
  ggplot2::scale_color_manual(values=line_colors) +
  ggplot2::scale_linetype_manual(
    values=c("All (pooled)"="solid", Upland="longdash",
             Freshwater="dashed", Saline="dotdash")) +
  ggplot2::scale_x_continuous(breaks=seq(-4,2,1)) +
  ggplot2::labs(
    title    = "Ecological Drought Index (EDI): Logistic Calibration to SPEI-3",
    subtitle = paste0(
      "EDI = P(WUE_T < near-normal baseline | SPEI-3)\n",
      "Fitted by logistic regression on monthly resistance data; ",
      "ecosystem-stratified where sufficient data exist.\n",
      "Vertical dotted line: SPEI-3 = -1 (meteorological drought onset)"
    ),
    x        = "SPEI-3",
    y        = "Ecological Drought Index (EDI)",
    color    = "Model",
    linetype = "Model"
  ) +
  theme_wue()

save_plot(p_logistic, "EDI_plot_01_logistic_curve.png", width=10, height=7)

# -- Plot 2: calibration scatter — monthly resistance vs SPEI_3 ------------
cal_df <- base_df[
  !is.na(base_df$monthly_resistance) &
  is.finite(base_df$monthly_resistance) &
  base_df$SPEI_3 < 0 &
  base_df$monthly_resistance < 5,  # trim extreme outliers for display
]

p_calib <- ggplot2::ggplot(cal_df,
    ggplot2::aes(x=SPEI_3, y=monthly_resistance, color=water_class)) +
  ggplot2::geom_hline(yintercept=1, linetype="dashed", color="gray35") +
  ggplot2::geom_vline(xintercept=-1, linetype="dotted", color="gray45") +
  ggplot2::geom_point(alpha=0.20, size=1.0) +
  ggplot2::geom_smooth(method="loess", se=TRUE, span=0.5, linewidth=1.0) +
  ggplot2::scale_color_manual(values=ecosystem_colors) +
  ggplot2::facet_wrap(~water_class, nrow=1) +
  ggplot2::labs(
    title    = "Monthly Resistance vs SPEI-3: Calibration Data",
    subtitle = "Each point = one site-month with SPEI_3 < 0; dashed line = resistance of 1 (no impairment)\nLogistic model is fitted to P(resistance < 1) across this domain",
    x        = "SPEI-3",
    y        = "Monthly resistance (WUE_T / near-normal WUE_T)",
    color    = "Ecosystem class"
  ) +
  theme_wue() +
  ggplot2::theme(legend.position="none")

save_plot(p_calib, "EDI_plot_02_monthly_resistance_calibration.png", width=13, height=6)

# -- Plot 3: EDI vs SPEI-3 all site-months ----------------------------------
p_scatter <- ggplot2::ggplot(
    base_df[is.finite(base_df$SPEI_3), ],
    ggplot2::aes(x=SPEI_3, y=EDI, color=water_class)) +
  ggplot2::geom_point(alpha=0.20, size=0.9) +
  ggplot2::geom_smooth(method="loess", se=FALSE, linewidth=1.1, span=0.4) +
  ggplot2::geom_hline(yintercept=c(0.30,0.50,0.70),
                      linetype="dashed", color="gray45", linewidth=0.4) +
  ggplot2::geom_vline(xintercept=-1, linetype="dotted", color="gray45") +
  ggplot2::scale_color_manual(values=ecosystem_colors) +
  ggplot2::facet_wrap(~water_class, nrow=1) +
  ggplot2::labs(
    title    = "EDI vs SPEI-3: All Site-Months",
    subtitle = "Each point is one site-month; lines are loess smooths; horizontal dashes mark severity thresholds",
    x        = "SPEI-3",
    y        = "EDI",
    color    = "Ecosystem class"
  ) +
  theme_wue() +
  ggplot2::theme(legend.position="none")

save_plot(p_scatter, "EDI_plot_03_spei_edi_by_ecosystem.png", width=13, height=6)

# -- Plot 4: EDI time series for selected sites -----------------------------
featured_sites <- unique(c(
  resistance$site_name[order(resistance$resistance)][1:4],
  resistance$site_name[order(-resistance$resistance)][1:4]
))
featured_sites <- featured_sites[featured_sites %in% base_df$site_name]
featured_sites <- unique(featured_sites)[1:min(12, length(featured_sites))]

ts_df       <- base_df[base_df$site_name %in% featured_sites, ]
ts_df$date  <- as.Date(paste(ts_df$Year, ts_df$month, "15", sep="-"))
ts_df$EDI_class <- factor(ts_df$EDI_class, levels=edi_labels)

p_ts <- ggplot2::ggplot(ts_df,
    ggplot2::aes(x=date, y=EDI, fill=EDI_class)) +
  ggplot2::geom_col(width=25) +
  ggplot2::geom_hline(yintercept=c(0.30,0.50,0.70),
                      linetype="dashed", color="gray40", linewidth=0.3) +
  ggplot2::scale_fill_manual(values=severity_colors) +
  ggplot2::scale_x_date(date_breaks="3 years", date_labels="%Y") +
  ggplot2::facet_wrap(~site_name, scales="free_x", ncol=3) +
  ggplot2::labs(
    title    = "EDI Time Series: Sites with Lowest and Highest Resistance",
    subtitle = "Bar height = EDI; colour = severity class (None / Watch / Stress / Impact)",
    x        = NULL, y = "EDI", fill = "Severity"
  ) +
  theme_wue(base_size=11) +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1))

save_plot(p_ts, "EDI_plot_04_edi_timeseries_facet.png", width=14, height=10)

# -- Plot 5: severity class proportions --------------------------------------
severity_props <- base_df |>
  dplyr::group_by(water_class, EDI_class) |>
  dplyr::summarise(n=dplyr::n(), .groups="drop") |>
  dplyr::group_by(water_class) |>
  dplyr::mutate(prop=n/sum(n)*100) |>
  dplyr::ungroup()
severity_props$EDI_class <- factor(severity_props$EDI_class, levels=edi_labels)

p_props <- ggplot2::ggplot(
    severity_props,
    ggplot2::aes(x=water_class, y=prop, fill=EDI_class)) +
  ggplot2::geom_col(position="stack", width=0.6) +
  ggplot2::scale_fill_manual(values=severity_colors) +
  ggplot2::scale_y_continuous(labels=scales::percent_format(scale=1)) +
  ggplot2::labs(
    title    = "EDI Severity Class Proportions by Ecosystem",
    subtitle = "Proportion of all site-months in each EDI severity class",
    x        = "Ecosystem class", y = "% of site-months", fill = "EDI class"
  ) +
  theme_wue() +
  ggplot2::theme(legend.position="right")

save_plot(p_props, "EDI_plot_05_severity_class_proportions.png", width=9, height=7)

# -- Plot 6: validation — EDI class vs Q2 resistance & recovery -------------
val_long <- validation_df |>
  dplyr::select(EDI_class, water_class,
                mean_monthly_resist, mean_site_resistance, mean_recovery) |>
  tidyr::pivot_longer(
    cols      = c(mean_monthly_resist, mean_site_resistance, mean_recovery),
    names_to  = "metric",
    values_to = "value"
  )
val_long$metric <- factor(
  val_long$metric,
  levels = c("mean_monthly_resist", "mean_site_resistance", "mean_recovery"),
  labels = c("Monthly resistance\n(WUE_T / baseline, this month)",
             "Site resistance\n(mean drought / near-normal, Q2)",
             "Site recovery\n(post / pre-drought WUE_T, Q2)")
)
val_long$EDI_class <- factor(val_long$EDI_class, levels=edi_labels)

p_valid <- ggplot2::ggplot(
    val_long,
    ggplot2::aes(x=EDI_class, y=value, color=water_class, group=water_class)) +
  ggplot2::geom_hline(yintercept=1, linetype="dashed", color="gray40") +
  ggplot2::geom_line(linewidth=0.8) +
  ggplot2::geom_point(size=3.0) +
  ggplot2::facet_wrap(~metric, scales="free_y", nrow=1) +
  ggplot2::scale_color_manual(values=ecosystem_colors) +
  ggplot2::labs(
    title    = "EDI Validation Against Q2 Performance Metrics",
    subtitle = paste0(
      "Mean resistance and recovery within each EDI severity class (drought months, SPEI-3 < -1 only)\n",
      "Dashed line = ratio of 1 (no impairment). If EDI is calibrated correctly, values should",
      " decline with increasing class."
    ),
    x        = "EDI class",
    y        = "Mean performance ratio",
    color    = "Ecosystem class"
  ) +
  theme_wue()

save_plot(p_valid, "EDI_plot_06_threshold_validation.png", width=14, height=6)

# -- Plot 7: composite -------------------------------------------------------
composite <- (p_response | p_response_props) /
              (p_logistic | p_props) /
              p_calib /
              p_valid +
  patchwork::plot_annotation(
    title = "Ecological Drought Index (EDI): SPEI-3 and WUE_T Response",
    subtitle = "Response EDI shows expected WUE_T change; logistic EDI is retained as secondary impairment probability.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 18)),
      plot.margin = ggplot2::margin(t = 10, r = 8, b = 10, l = 8)
    )
  ) +
  patchwork::plot_layout(heights=c(1.35, 1.05, 1.0, 1.0))

save_plot(composite, "EDI_plot_07_composite.png", width=16, height=34)

# ---------------------------------------------------------------------------
# Summary log
# ---------------------------------------------------------------------------

capture_to_file({
  cat("Ecological_Impacts.R — Q3-aligned EDI run summary\n")
  cat("Input:", file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv"), "\n")
  cat("Output:", output_dir, "\n\n")

  cat("=== Q3-aligned purpose ===\n")
  cat("Primary ecological impact score is predicted upland WUE_T percent change from near-normal SPEI.\n")
  cat("Near-normal SPEI is defined as -1 <= SPEI <= 1 within each coastline and SPEI timescale.\n")
  cat("The model matches the final Q3 coast-threshold GAM logic.\n\n")

  cat("=== Model formula ===\n")
  print(formula(q3_coast_model))

  cat("\n=== Model information ===\n")
  print(q3_model_info)

  cat("\n=== Coast-threshold GAM summary ===\n")
  print(summary(q3_coast_model))

  cat("\n=== Smooth terms ===\n")
  print(q3_coast_smooth_table)

  cat("\n=== Upland coast-specific 5/10/20% SPEI thresholds ===\n")
  print(q3_coast_threshold_summary)

  cat("\n=== SPEI-3 monthly upland EDI class counts by coastline ===\n")
  print(table(q3_monthly_out$coast_region, q3_monthly_out$EDI_Q3_class, useNA = "ifany"))

  cat("\n=== SPEI-3 monthly upland EDI direction counts by coastline ===\n")
  print(table(q3_monthly_out$coast_region, q3_monthly_out$EDI_Q3_direction, useNA = "ifany"))

  cat("\n=== SPEI-3 prediction-interval support by coastline ===\n")
  print(table(q3_monthly_out$coast_region, q3_monthly_out$prediction_interval_support, useNA = "ifany"))

  cat("\n=== Regional validation against Q3 site-sensitivity summaries when available ===\n")
  print(q3_regional_validation)
}, "EDI_Q3_aligned_model_summary.txt")

capture_to_file({
  cat("Ecological_Impacts.R — run summary\n")
  cat("SPEI timescale used: SPEI_3\n")
  cat("Output: ", output_dir, "\n\n")

  cat("=== Calibration data ===\n")
  cat("Total site-months with SPEI_3 < 0 and valid baseline:", nrow(model_df), "\n")
  cat("Of which resistance < 1 (impacted):                  ", sum(model_df$y), "\n")
  cat("Proportion impacted:", round(mean(model_df$y)*100, 1), "%\n\n")

  cat("=== Response-based EDI model ===\n")
  cat("Response: WUE_T percent change from site-month near-normal baseline\n")
  cat("Rows with valid response baseline:", nrow(response_model_df), "\n")
  print(summary(response_model))

  cat("\n=== Response-based EDI thresholds ===\n")
  print(response_thresholds)

  cat("\n=== Response class counts across all site-months ===\n")
  print(table(base_df$water_class, base_df$response_EDI_class, useNA = "ifany"))
  cat("\n=== Response direction counts across all site-months ===\n")
  print(table(base_df$water_class, base_df$response_EDI_direction, useNA = "ifany"))

  cat("=== Pooled logistic model ===\n")
  print(summary(pooled_model))

  cat("\n=== Ecosystem-stratified models ===\n")
  for (ec in ecosystem_classes) {
    m <- eco_models[[ec]]
    if (is.null(m)) {
      cat(ec, ": model not fitted (insufficient events)\n")
    } else {
      cat("\n---", ec, "---\n")
      print(summary(m))
    }
  }

  cat("\n=== SPEI-3 thresholds at EDI class boundaries ===\n")
  print(thresholds)

  cat("\n=== EDI validation against Q2 metrics (drought months, SPEI_3 < -1) ===\n")
  print(validation_df)

  cat("\n=== Severity class counts across all site-months ===\n")
  print(table(base_df$water_class, base_df$EDI_class))
}, "EDI_logistic_model_summary.txt")

message("Done. Outputs written to: ", output_dir)
