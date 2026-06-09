# Q4 WUE decline thresholds across SPEI gradients
#
# Goal:
#   Explore at what SPEI values WUE is likely to decline relative to
#   near-normal long-term moisture conditions.
#
# Decline definition:
#   Site-specific percent change from near-normal SPEI conditions.
#   Near-normal is -1 <= SPEI <= 1.
#   Declines are evaluated at 0.5%, 1%, 1.5%, 2.5%, 5%, 10%, and 20% thresholds.
#   Main model focuses on SPEI_48 and excludes coastline.
#
# WUE metrics:
#   WUE_ET = WUE
#   WUE_T  = WUE_tra
#
# Inputs:
#   Water_WUE/data/WUE_CUE_monthly_merged_indices_clean.csv
#
# Outputs:
#   Water_WUE/data/results/Q4_WUE_decline_SPEI/

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
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  }
  normalizePath("Water_WUE/function/Q4.Decline.R", mustWork = FALSE)
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

find_probability_threshold <- function(group_data, probability_column, side, probability_cutoff = 0.5) {
  side_filter <- if (side == "dry") group_data$SPEI_value < -1 else group_data$SPEI_value > 1
  threshold_values <- group_data$SPEI_value[
    side_filter & group_data[[probability_column]] >= probability_cutoff
  ]
  if (length(threshold_values) == 0) return(NA_real_)
  if (side == "dry") max(threshold_values, na.rm = TRUE) else min(threshold_values, na.rm = TRUE)
}

find_decline_threshold <- function(group_data, pct_column, side, decline_threshold) {
  side_filter <- if (side == "dry") group_data$SPEI_value < -1 else group_data$SPEI_value > 1
  threshold_values <- group_data$SPEI_value[
    side_filter & group_data[[pct_column]] <= -decline_threshold
  ]
  if (length(threshold_values) == 0) return(NA_real_)
  if (side == "dry") max(threshold_values, na.rm = TRUE) else min(threshold_values, na.rm = TRUE)
}

script_dir <- dirname(get_script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

input_file <- file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv")
output_dir <- file.path(project_dir, "data", "results", "Q4_WUE_decline_SPEI")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
monthly[] <- lapply(monthly, function(x) {
  if (is.character(x)) trimws(x) else x
})

spei_cols <- c("SPEI_48")
required_cols <- c(
  "site_name", "Year", "month", "water_class", "lat", "long",
  "WUE", "WUE_tra", spei_cols
)

missing_cols <- setdiff(required_cols, names(monthly))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ecosystem_classes <- c("Upland", "Freshwater", "Saline")
coast_region_levels <- c("Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast")

model_base <- monthly[complete.cases(monthly[, required_cols]), ]
model_base <- model_base[
  model_base$site_name != "" &
    model_base$water_class != "" &
    model_base$water_class %in% ecosystem_classes &
    is.finite(model_base$lat) &
    is.finite(model_base$long) &
    is.finite(model_base$WUE) &
    is.finite(model_base$WUE_tra),
]

model_base$coast_region <- factor(
  classify_coast_region(model_base$lat, model_base$long),
  levels = coast_region_levels
)
model_base$site_name <- factor(model_base$site_name)
model_base$water_class <- factor(model_base$water_class, levels = ecosystem_classes)
model_base$month_f <- factor(model_base$month, levels = 1:12)
model_base$WUE_ET <- model_base$WUE
model_base$WUE_T <- model_base$WUE_tra

metric_long <- reshape(
  model_base,
  varying = c("WUE_ET", "WUE_T"),
  v.names = "WUE_value",
  timevar = "metric",
  times = c("WUE_ET", "WUE_T"),
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)
metric_long$metric <- factor(metric_long$metric, levels = c("WUE_ET", "WUE_T"))

spei_long <- do.call(
  rbind,
  lapply(spei_cols, function(spei_col) {
    out <- metric_long
    out$SPEI_timescale <- spei_col
    out$SPEI_value <- out[[spei_col]]
    out
  })
)

spei_long <- spei_long[is.finite(spei_long$SPEI_value) & is.finite(spei_long$WUE_value), ]
spei_long$SPEI_timescale <- factor(spei_long$SPEI_timescale, levels = spei_cols)
spei_long$metric <- factor(spei_long$metric, levels = c("WUE_ET", "WUE_T"))
spei_long$site_name <- factor(spei_long$site_name)
spei_long$water_class <- factor(spei_long$water_class, levels = ecosystem_classes)
spei_long$coast_region <- factor(spei_long$coast_region, levels = coast_region_levels)
spei_long$month_f <- factor(spei_long$month, levels = 1:12)
spei_long$model_group <- spei_long$metric

site_baselines <- spei_long |>
  dplyr::filter(SPEI_value >= -1, SPEI_value <= 1) |>
  dplyr::group_by(site_name, metric, SPEI_timescale, month_f) |>
  dplyr::summarise(
    baseline_WUE = mean(WUE_value),
    baseline_n = dplyr::n(),
    .groups = "drop"
  )

analysis_data <- dplyr::left_join(
  spei_long,
  site_baselines,
  by = c("site_name", "metric", "SPEI_timescale", "month_f")
)
analysis_data <- analysis_data[
  is.finite(analysis_data$baseline_WUE) &
    analysis_data$baseline_WUE > 0 &
    analysis_data$baseline_n >= 3,
]
analysis_data$WUE_change <- analysis_data$WUE_value - analysis_data$baseline_WUE
analysis_data$WUE_pct_change <- 100 * analysis_data$WUE_change / analysis_data$baseline_WUE
analysis_data$decline_2_5 <- as.integer(analysis_data$WUE_pct_change <= -2.5)
analysis_data$decline_5 <- as.integer(analysis_data$WUE_pct_change <= -5)
analysis_data$decline_10 <- as.integer(analysis_data$WUE_pct_change <= -10)
analysis_data$decline_20 <- as.integer(analysis_data$WUE_pct_change <= -20)

write_table(model_base, "Q4_WUE_decline_model_base_wide.csv")
write_table(analysis_data, "Q4_WUE_decline_model_data_long.csv")
write_table(site_baselines, "Q4_WUE_decline_site_near_normal_baselines.csv")

coverage <- analysis_data |>
  dplyr::group_by(metric, SPEI_timescale, water_class, coast_region) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_pct_change = mean(WUE_pct_change),
    pct_decline_2_5 = mean(decline_2_5) * 100,
    pct_decline_5 = mean(decline_5) * 100,
    pct_decline_10 = mean(decline_10) * 100,
    pct_decline_20 = mean(decline_20) * 100,
    min_SPEI = min(SPEI_value),
    max_SPEI = max(SPEI_value),
    .groups = "drop"
  )
write_table(coverage, "Q4_WUE_decline_coverage.csv")

site_decline_thresholds <- analysis_data |>
  dplyr::group_by(site_name, water_class, coast_region, metric, SPEI_timescale) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    baseline_WUE = dplyr::first(baseline_WUE),
    mean_pct_change = mean(WUE_pct_change),
    min_pct_change = min(WUE_pct_change),
    pct_months_decline_2_5 = mean(decline_2_5) * 100,
    pct_months_decline_5 = mean(decline_5) * 100,
    pct_months_decline_10 = mean(decline_10) * 100,
    pct_months_decline_20 = mean(decline_20) * 100,
    dry_decline_2_5_threshold = ifelse(
      any(SPEI_value < -1 & WUE_pct_change <= -2.5),
      max(SPEI_value[SPEI_value < -1 & WUE_pct_change <= -2.5], na.rm = TRUE),
      NA_real_
    ),
    dry_decline_5_threshold = ifelse(
      any(SPEI_value < -1 & WUE_pct_change <= -5),
      max(SPEI_value[SPEI_value < -1 & WUE_pct_change <= -5], na.rm = TRUE),
      NA_real_
    ),
    dry_decline_10_threshold = ifelse(
      any(SPEI_value < -1 & WUE_pct_change <= -10),
      max(SPEI_value[SPEI_value < -1 & WUE_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    dry_decline_20_threshold = ifelse(
      any(SPEI_value < -1 & WUE_pct_change <= -20),
      max(SPEI_value[SPEI_value < -1 & WUE_pct_change <= -20], na.rm = TRUE),
      NA_real_
    ),
    wet_decline_5_threshold = ifelse(
      any(SPEI_value > 1 & WUE_pct_change <= -5),
      min(SPEI_value[SPEI_value > 1 & WUE_pct_change <= -5], na.rm = TRUE),
      NA_real_
    ),
    wet_decline_2_5_threshold = ifelse(
      any(SPEI_value > 1 & WUE_pct_change <= -2.5),
      min(SPEI_value[SPEI_value > 1 & WUE_pct_change <= -2.5], na.rm = TRUE),
      NA_real_
    ),
    wet_decline_10_threshold = ifelse(
      any(SPEI_value > 1 & WUE_pct_change <= -10),
      min(SPEI_value[SPEI_value > 1 & WUE_pct_change <= -10], na.rm = TRUE),
      NA_real_
    ),
    wet_decline_20_threshold = ifelse(
      any(SPEI_value > 1 & WUE_pct_change <= -20),
      min(SPEI_value[SPEI_value > 1 & WUE_pct_change <= -20], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )
write_table(site_decline_thresholds, "Q4_WUE_decline_site_thresholds_5_10_20pct.csv")

change_model <- mgcv::gam(
  WUE_pct_change ~
    metric + water_class + month_f +
    s(SPEI_value, by = model_group, k = 6) +
    s(site_name, bs = "re"),
  data = analysis_data,
  method = "REML",
  select = TRUE
)
saveRDS(change_model, file.path(output_dir, "Q4_WUE_pct_change_SPEI48_gam.rds"))

smooth_table <- as.data.frame(summary(change_model)$s.table)
smooth_table$smooth_term <- rownames(smooth_table)
rownames(smooth_table) <- NULL
smooth_table <- smooth_table[, c("smooth_term", setdiff(names(smooth_table), "smooth_term"))]
write_table(smooth_table, "Q4_WUE_pct_change_SPEI48_gam_smooth_terms.csv")

parametric_table <- as.data.frame(summary(change_model)$p.table)
parametric_table$term <- rownames(parametric_table)
rownames(parametric_table) <- NULL
parametric_table <- parametric_table[, c("term", setdiff(names(parametric_table), "term"))]
write_table(parametric_table, "Q4_WUE_pct_change_SPEI48_gam_parametric_terms.csv")

spei_sequence <- seq(
  stats::quantile(analysis_data$SPEI_value, 0.02),
  stats::quantile(analysis_data$SPEI_value, 0.98),
  length.out = 180
)
prediction_grid <- expand.grid(
  metric = levels(analysis_data$metric),
  water_class = levels(analysis_data$water_class),
  SPEI_value = spei_sequence,
  month_f = factor("7", levels = levels(analysis_data$month_f)),
  site_name = levels(analysis_data$site_name)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
prediction_grid$metric <- factor(prediction_grid$metric, levels = levels(analysis_data$metric))
prediction_grid$SPEI_timescale <- factor("SPEI_48", levels = levels(analysis_data$SPEI_timescale))
prediction_grid$water_class <- factor(prediction_grid$water_class, levels = levels(analysis_data$water_class))
prediction_grid$month_f <- factor(prediction_grid$month_f, levels = levels(analysis_data$month_f))
prediction_grid$site_name <- factor(prediction_grid$site_name, levels = levels(analysis_data$site_name))
prediction_grid$model_group <- prediction_grid$metric

pred_change <- predict(
  change_model,
  newdata = prediction_grid,
  se.fit = TRUE,
  exclude = "s(site_name)"
)
prediction_grid$predicted_pct_change <- as.numeric(pred_change$fit)
prediction_grid$predicted_pct_change_lower <- prediction_grid$predicted_pct_change - 1.96 * as.numeric(pred_change$se.fit)
prediction_grid$predicted_pct_change_upper <- prediction_grid$predicted_pct_change + 1.96 * as.numeric(pred_change$se.fit)
write_table(prediction_grid, "Q4_WUE_pct_change_SPEI48_predictions.csv")

decline_thresholds <- do.call(
  rbind,
  lapply(
    split(prediction_grid, list(prediction_grid$metric, prediction_grid$water_class), drop = TRUE),
    function(group_data) {
      do.call(
        rbind,
        lapply(c(0.5, 1, 1.5, 2.5, 5, 10, 20), function(threshold) {
          data.frame(
            metric = group_data$metric[1],
            SPEI_timescale = group_data$SPEI_timescale[1],
            water_class = group_data$water_class[1],
            decline_threshold_pct = threshold,
            dry_SPEI_threshold = find_decline_threshold(group_data, "predicted_pct_change", "dry", threshold),
            dry_SPEI_threshold_lower = find_decline_threshold(group_data, "predicted_pct_change_upper", "dry", threshold),
            dry_SPEI_threshold_upper = find_decline_threshold(group_data, "predicted_pct_change_lower", "dry", threshold),
            wet_SPEI_threshold = find_decline_threshold(group_data, "predicted_pct_change", "wet", threshold),
            wet_SPEI_threshold_lower = find_decline_threshold(group_data, "predicted_pct_change_upper", "wet", threshold),
            wet_SPEI_threshold_upper = find_decline_threshold(group_data, "predicted_pct_change_lower", "wet", threshold),
            min_predicted_pct_change = min(group_data$predicted_pct_change),
            max_predicted_pct_change = max(group_data$predicted_pct_change)
          )
        })
      )
    }
  )
)
write_table(decline_thresholds, "Q4_WUE_pct_change_SPEI48_decline_thresholds.csv")

ecosystem_colors <- c(
  "Upland" = "#800080",
  "Freshwater" = "#0000FF",
  "Saline" = "#FFA500"
)

selected_timescales <- c("SPEI_48")
plot_prediction_data <- prediction_grid[
  prediction_grid$metric == "WUE_T",
]

plot_pct_change <- ggplot2::ggplot(
  plot_prediction_data,
  ggplot2::aes(x = SPEI_value, y = predicted_pct_change, color = water_class, fill = water_class)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray45") +
  ggplot2::geom_hline(yintercept = 0, color = "gray35") +
  ggplot2::geom_hline(yintercept = c(-0.5, -1, -1.5, -2.5, -5, -10, -20), linetype = "dotted", color = "gray45") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = predicted_pct_change_lower, ymax = predicted_pct_change_upper),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  ggplot2::labs(
    title = "Predicted WUE_T Change Across SPEI-48",
    subtitle = "Continuous response model; dotted lines mark 0.5%, 1%, 1.5%, 2.5%, 5%, 10%, and 20% declines",
    x = "SPEI-48",
    y = "Predicted WUE_T change from near-normal (%)",
    color = "Ecosystem class",
    fill = "Ecosystem class"
  ) +
  theme_wue()

save_plot(plot_pct_change, "Q4_plot_01_WUE_T_pct_change_SPEI48.png", width = 10, height = 7)

observed_threshold_plot_data <- site_decline_thresholds |>
  dplyr::filter(metric == "WUE_T", SPEI_timescale %in% selected_timescales, n_months >= 6) |>
  dplyr::mutate(
    dry_decline_2_5_detected = !is.na(dry_decline_2_5_threshold),
    wet_decline_2_5_detected = !is.na(wet_decline_2_5_threshold),
    any_decline_2_5_detected = dry_decline_2_5_detected | wet_decline_2_5_detected
  )

site_order <- observed_threshold_plot_data |>
  dplyr::group_by(site_name) |>
  dplyr::summarise(
    max_pct_decline = max(pct_months_decline_2_5, na.rm = TRUE),
    .groups = "drop"
  )
observed_threshold_plot_data$site_name <- factor(
  observed_threshold_plot_data$site_name,
  levels = site_order$site_name[order(site_order$max_pct_decline)]
)

plot_site_decline <- ggplot2::ggplot(
  observed_threshold_plot_data,
  ggplot2::aes(x = pct_months_decline_2_5, y = site_name, color = water_class, shape = SPEI_timescale)
) +
  ggplot2::geom_vline(xintercept = c(25, 50, 75), linetype = "dotted", color = "gray65") +
  ggplot2::geom_point(
    size = 2.8,
    alpha = 0.88,
    position = ggplot2::position_jitter(height = 0.12, width = 0)
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::labs(
    title = "Site-Level Frequency of WUE_T Decline",
    subtitle = "Percent of months at each site with WUE_T at least 2.5% below near-normal baseline",
    x = "Months with >=2.5% WUE_T decline (%)",
    y = NULL,
    color = "Ecosystem class",
    shape = "SPEI timescale"
  ) +
  theme_wue(base_size = 11) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6),
    panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.25)
  )

save_plot(plot_site_decline, "Q4_plot_02_site_level_WUE_T_decline_frequency.png", width = 13, height = 13)

threshold_heatmap_data <- decline_thresholds |>
  dplyr::filter(metric == "WUE_T", decline_threshold_pct == 2.5)

plot_thresholds <- ggplot2::ggplot(
  threshold_heatmap_data,
  ggplot2::aes(x = water_class, y = "SPEI_48", fill = dry_SPEI_threshold)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(is.na(dry_SPEI_threshold), "none", round(dry_SPEI_threshold, 2))),
    size = 4
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#6A3D9A",
    mid = "white",
    high = "#F6AD55",
    midpoint = -1,
    na.value = "gray90"
  ) +
  ggplot2::labs(
    title = "Dry SPEI-48 Threshold for Predicted WUE_T Decline",
    subtitle = "Cells show dry-side SPEI where predicted WUE_T change crosses -2.5%; 'none' means no model threshold detected",
    x = "Ecosystem class",
    y = NULL,
    fill = "Dry SPEI\nthreshold"
  ) +
  theme_wue()

save_plot(plot_thresholds, "Q4_plot_03_dry_SPEI48_threshold_for_WUE_T_decline.png", width = 8, height = 4)

composite <- plot_pct_change / plot_thresholds / plot_site_decline +
  patchwork::plot_annotation(
    title = "Q4 SPEI-48 Thresholds for WUE Decline",
    subtitle = "Continuous percent-change model with declines defined relative to each site's near-normal SPEI baseline."
  )

save_plot(composite, "Q4_plot_04_WUE_decline_SPEI_multipanel.png", width = 14, height = 24)

capture_to_file({
  cat("Q4 WUE decline thresholds across SPEI gradients\n")
  cat("Input:", input_file, "\n")
  cat("Output:", output_dir, "\n\n")
  cat("Rows in source monthly data:", nrow(monthly), "\n")
  cat("Rows after complete-case and ecosystem filters:", nrow(model_base), "\n")
  cat("Rows in SPEI decline analysis data:", nrow(analysis_data), "\n")
  cat("Sites:", length(unique(analysis_data$site_name)), "\n")
  cat("Years:", min(analysis_data$Year), "-", max(analysis_data$Year), "\n\n")

  cat("Continuous WUE percent-change model:\n")
  print(summary(change_model))

  cat("\nSelected WUE_T decline thresholds from continuous SPEI-48 model:\n")
  print(
    decline_thresholds[
      decline_thresholds$metric == "WUE_T",
      c("SPEI_timescale", "water_class", "decline_threshold_pct", "dry_SPEI_threshold", "wet_SPEI_threshold", "min_predicted_pct_change", "max_predicted_pct_change")
    ]
  )
}, "Q4_WUE_decline_model_summary.txt")

cat("Done. Q4 WUE decline outputs written to:", output_dir, "\n")
