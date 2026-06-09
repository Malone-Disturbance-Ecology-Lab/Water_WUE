# Q2 moisture-anomaly analysis for WUE_ET and WUE_T
#
# Goal:
#   Determine whether short-term (SPEI-6) and long-term (SPEI-48)
#   moisture anomalies alter the magnitude and direction of WUE responses.
#
# Definitions:
#   Dry anomaly: SPEI < -1
#   Near-normal: -1 <= SPEI <= 1
#   Wet anomaly: SPEI > 1
#
# WUE metrics:
#   WUE_ET = WUE
#   WUE_T  = WUE_tra
#
# Inputs:
#   Water_WUE/data/WUE_CUE_monthly_merged_indices_clean.csv
#
# Outputs:
#   Water_WUE/data/results/Q2_anomalies/

suppressPackageStartupMessages({
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required. Install it with install.packages('lme4').")
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
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required. Install it with install.packages('emmeans').")
  }
  if (!requireNamespace("multcomp", quietly = TRUE)) {
    stop("Package 'multcomp' is required. Install it with install.packages('multcomp').")
  }
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  }
  normalizePath("Water_WUE/function/Q2_Anomolies.R", mustWork = FALSE)
}

write_table <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE)
}

capture_to_file <- function(expr, filename) {
  sink(file.path(output_dir, filename))
  on.exit(sink(), add = TRUE)
  force(expr)
}

safe_drop1 <- function(model) {
  out <- tryCatch(
    drop1(model, test = "Chisq"),
    error = function(e) {
      message("drop1 failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(out)) {
    return(NULL)
  }

  out_df <- as.data.frame(out)
  out_df$term <- rownames(out_df)
  rownames(out_df) <- NULL
  out_df[, c("term", setdiff(names(out_df), "term"))]
}

coef_table <- function(model) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df[, c("term", setdiff(names(coef_df), "term"))]
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

classify_spei <- function(x) {
  out <- ifelse(x < -1, "Dry", ifelse(x > 1, "Wet", "Near-normal"))
  factor(out, levels = c("Dry", "Near-normal", "Wet"))
}

p_value_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}

script_dir <- dirname(get_script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

input_file <- file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv")
output_dir <- file.path(project_dir, "data", "results", "Q2_anomalies")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
monthly[] <- lapply(monthly, function(x) {
  if (is.character(x)) trimws(x) else x
})

required_cols <- c(
  "site_name", "Year", "month", "water_class",
  "WUE", "WUE_tra", "Trans_ratio", "SPEI_6", "SPEI_48"
)

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
    is.finite(model_base$WUE) &
    is.finite(model_base$WUE_tra) &
    is.finite(model_base$Trans_ratio) &
    model_base$Trans_ratio >= 0 &
    model_base$Trans_ratio <= 1 &
    is.finite(model_base$SPEI_6) &
    is.finite(model_base$SPEI_48),
]

model_base$site_name <- factor(model_base$site_name)
model_base$water_class <- factor(model_base$water_class, levels = ecosystem_classes)
model_base$month_f <- factor(model_base$month, levels = 1:12)
model_base$Trans_ratio_z <- as.numeric(scale(model_base$Trans_ratio))
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

spei_long <- rbind(
  transform(metric_long, timescale = "SPEI_6", SPEI_value = SPEI_6),
  transform(metric_long, timescale = "SPEI_48", SPEI_value = SPEI_48)
)
spei_long$timescale <- factor(spei_long$timescale, levels = c("SPEI_6", "SPEI_48"))
spei_long$anomaly_class <- classify_spei(spei_long$SPEI_value)
spei_long$water_class <- factor(spei_long$water_class, levels = ecosystem_classes)
spei_long$site_name <- factor(spei_long$site_name)
spei_long$metric <- factor(spei_long$metric, levels = c("WUE_ET", "WUE_T"))
spei_long$month_f <- factor(spei_long$month, levels = 1:12)

baseline <- spei_long |>
  dplyr::filter(anomaly_class == "Near-normal") |>
  dplyr::group_by(site_name, timescale, metric, month_f) |>
  dplyr::summarise(
    baseline_mean = mean(WUE_value),
    baseline_median = median(WUE_value),
    baseline_n = dplyr::n(),
    .groups = "drop"
  )

analysis_data <- dplyr::left_join(
  spei_long,
  baseline,
  by = c("site_name", "timescale", "metric", "month_f")
)
analysis_data <- analysis_data[
  is.finite(analysis_data$baseline_mean) &
    analysis_data$baseline_n >= 3,
]
analysis_data$WUE_response <- analysis_data$WUE_value - analysis_data$baseline_mean
analysis_data$WUE_response_magnitude <- abs(analysis_data$WUE_response)
analysis_data$response_direction <- factor(
  ifelse(analysis_data$WUE_response < 0, "Decrease", "Increase_or_no_decrease"),
  levels = c("Increase_or_no_decrease", "Decrease")
)
analysis_data$decrease_binary <- as.integer(analysis_data$response_direction == "Decrease")

anomaly_data <- analysis_data[analysis_data$anomaly_class != "Near-normal", ]
anomaly_data$anomaly_class <- factor(as.character(anomaly_data$anomaly_class), levels = c("Dry", "Wet"))

write_table(model_base, "Q2_anomalies_model_base_wide.csv")
write_table(analysis_data, "Q2_anomalies_model_data_long.csv")
write_table(baseline, "Q2_anomalies_near_normal_baselines.csv")

coverage_by_class <- analysis_data |>
  dplyr::group_by(timescale, anomaly_class, metric, water_class) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_spei = mean(SPEI_value),
    mean_wue = mean(WUE_value),
    mean_response = mean(WUE_response),
    mean_abs_response = mean(WUE_response_magnitude),
    pct_decrease = mean(decrease_binary) * 100,
    .groups = "drop"
  )
write_table(coverage_by_class, "Q2_anomalies_summary_by_class_metric_ecosystem.csv")

summary_by_class <- analysis_data |>
  dplyr::group_by(timescale, anomaly_class, metric) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_wue = mean(WUE_value),
    median_wue = median(WUE_value),
    sd_wue = sd(WUE_value),
    se_wue = sd_wue / sqrt(n_months),
    ci95_wue = stats::qt(0.975, df = n_months - 1) * se_wue,
    mean_response = mean(WUE_response),
    se_response = sd(WUE_response) / sqrt(n_months),
    ci95_response = stats::qt(0.975, df = n_months - 1) * se_response,
    mean_abs_response = mean(WUE_response_magnitude),
    pct_decrease = mean(decrease_binary) * 100,
    .groups = "drop"
  )
write_table(summary_by_class, "Q2_anomalies_summary_by_class_metric.csv")

# ---------------------------------------------------------------------------
# Mixed models
# ---------------------------------------------------------------------------

wue_formula <- WUE_value ~ metric * timescale * anomaly_class * water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

response_formula <- WUE_response ~ metric * timescale * anomaly_class * water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

magnitude_formula <- WUE_response_magnitude ~ metric * timescale * anomaly_class * water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

direction_formula <- decrease_binary ~ metric * timescale * anomaly_class * water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

wue_model <- lme4::lmer(wue_formula, data = analysis_data, REML = FALSE)
response_model <- lme4::lmer(response_formula, data = analysis_data, REML = FALSE)
magnitude_model <- lme4::lmer(magnitude_formula, data = anomaly_data, REML = FALSE)
direction_model <- lme4::glmer(
  direction_formula,
  data = anomaly_data,
  family = stats::binomial(),
  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

saveRDS(wue_model, file.path(output_dir, "Q2_anomalies_wue_value_model.rds"))
saveRDS(response_model, file.path(output_dir, "Q2_anomalies_signed_response_model.rds"))
saveRDS(magnitude_model, file.path(output_dir, "Q2_anomalies_response_magnitude_model.rds"))
saveRDS(direction_model, file.path(output_dir, "Q2_anomalies_response_direction_model.rds"))

write_table(coef_table(wue_model), "Q2_anomalies_wue_value_fixed_effects.csv")
write_table(coef_table(response_model), "Q2_anomalies_signed_response_fixed_effects.csv")
write_table(coef_table(magnitude_model), "Q2_anomalies_response_magnitude_fixed_effects.csv")
write_table(coef_table(direction_model), "Q2_anomalies_response_direction_fixed_effects.csv")

wue_drop1 <- safe_drop1(wue_model)
response_drop1 <- safe_drop1(response_model)
magnitude_drop1 <- safe_drop1(magnitude_model)
direction_drop1 <- safe_drop1(direction_model)

if (!is.null(wue_drop1)) write_table(wue_drop1, "Q2_anomalies_wue_value_drop1_lrt.csv")
if (!is.null(response_drop1)) write_table(response_drop1, "Q2_anomalies_signed_response_drop1_lrt.csv")
if (!is.null(magnitude_drop1)) write_table(magnitude_drop1, "Q2_anomalies_response_magnitude_drop1_lrt.csv")
if (!is.null(direction_drop1)) write_table(direction_drop1, "Q2_anomalies_response_direction_drop1_lrt.csv")

model_results <- rbind(
  data.frame(
    model = "WUE value",
    AIC = AIC(wue_model),
    BIC = BIC(wue_model),
    logLik = as.numeric(stats::logLik(wue_model)),
    response = "WUE_value"
  ),
  data.frame(
    model = "Signed response",
    AIC = AIC(response_model),
    BIC = BIC(response_model),
    logLik = as.numeric(stats::logLik(response_model)),
    response = "WUE_response"
  ),
  data.frame(
    model = "Response magnitude",
    AIC = AIC(magnitude_model),
    BIC = BIC(magnitude_model),
    logLik = as.numeric(stats::logLik(magnitude_model)),
    response = "abs(WUE_response)"
  ),
  data.frame(
    model = "Direction",
    AIC = AIC(direction_model),
    BIC = BIC(direction_model),
    logLik = as.numeric(stats::logLik(direction_model)),
    response = "decrease_binary"
  )
)
write_table(model_results, "Q2_anomalies_model_results_table.csv")

wue_anomaly_emmeans <- emmeans::emmeans(
  wue_model,
  specs = ~ anomaly_class | timescale * metric
)
wue_anomaly_letters <- as.data.frame(multcomp::cld(
  wue_anomaly_emmeans,
  adjust = "tukey",
  Letters = letters
))
wue_anomaly_letters$.group <- gsub("\\s+", "", wue_anomaly_letters$.group)
write_table(wue_anomaly_letters, "Q2_anomalies_WUE_anomaly_class_letters.csv")

wue_anomaly_pairs <- as.data.frame(emmeans::contrast(
  wue_anomaly_emmeans,
  method = "pairwise",
  adjust = "tukey"
))
write_table(wue_anomaly_pairs, "Q2_anomalies_WUE_anomaly_class_pairwise_contrasts.csv")

response_anomaly_emmeans <- emmeans::emmeans(
  response_model,
  specs = ~ anomaly_class | timescale * metric
)
response_anomaly_letters <- as.data.frame(multcomp::cld(
  response_anomaly_emmeans,
  adjust = "tukey",
  Letters = letters
))
response_anomaly_letters$.group <- gsub("\\s+", "", response_anomaly_letters$.group)
write_table(response_anomaly_letters, "Q2_anomalies_signed_response_anomaly_class_letters.csv")

response_anomaly_pairs <- as.data.frame(emmeans::contrast(
  response_anomaly_emmeans,
  method = "pairwise",
  adjust = "tukey"
))
write_table(response_anomaly_pairs, "Q2_anomalies_signed_response_anomaly_class_pairwise_contrasts.csv")

analysis_data$wue_fitted <- fitted(wue_model)
analysis_data$wue_residual <- residuals(wue_model)
write_table(
  analysis_data[, c(
    "site_name", "Year", "month", "water_class", "timescale", "anomaly_class",
    "SPEI_value", "metric", "WUE_value", "baseline_mean", "WUE_response",
    "WUE_response_magnitude", "response_direction", "wue_fitted", "wue_residual"
  )],
  "Q2_anomalies_fitted_residuals.csv"
)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

ecosystem_colors <- c(
  "Upland" = "#800080",
  "Freshwater" = "#0000FF",
  "Saline" = "#FFA500"
)

anomaly_colors <- c(
  "Dry" = "#8C510A",
  "Near-normal" = "gray70",
  "Wet" = "#2B6CB0"
)

wue_letter_plot <- wue_anomaly_letters
wue_upper_col <- intersect(c("upper.CL", "asymp.UCL"), names(wue_letter_plot))[1]
if (is.na(wue_upper_col)) {
  wue_letter_plot$label_y <- wue_letter_plot$emmean
} else {
  wue_letter_plot$label_y <- wue_letter_plot[[wue_upper_col]]
}
wue_letter_plot <- wue_letter_plot |>
  dplyr::group_by(timescale, metric) |>
  dplyr::mutate(
    label_y = label_y + 0.08 * diff(range(label_y, na.rm = TRUE))
  ) |>
  dplyr::ungroup()
wue_letter_plot$anomaly_class <- factor(
  wue_letter_plot$anomaly_class,
  levels = levels(analysis_data$anomaly_class)
)
wue_letter_plot$metric <- factor(wue_letter_plot$metric, levels = levels(analysis_data$metric))

plot_wue_by_anomaly <- ggplot2::ggplot(
  summary_by_class,
  ggplot2::aes(x = anomaly_class, y = mean_wue, color = metric, group = metric)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "gray40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_wue - ci95_wue, ymax = mean_wue + ci95_wue),
    width = 0.12,
    linewidth = 0.7,
    position = ggplot2::position_dodge(width = 0.35)
  ) +
  ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.35)) +
  ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.35), linewidth = 0.8) +
  ggplot2::geom_text(
    data = wue_letter_plot,
    ggplot2::aes(x = anomaly_class, y = label_y, label = .group, color = metric, group = metric),
    position = ggplot2::position_dodge(width = 0.35),
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ timescale, nrow = 1, scales = "free_y") +
  ggplot2::labs(
    title = "WUE_ET and WUE_T Across Moisture Anomaly Classes",
    subtitle = "Letters show Tukey-adjusted model-based differences among anomaly classes within each metric and SPEI timescale",
    x = "Moisture anomaly class",
    y = "Mean WUE",
    color = "Metric"
  ) +
  ggplot2::scale_color_manual(values = c("WUE_ET" = "gray20", "WUE_T" = "#6A3D9A")) +
  theme_wue()

save_plot(plot_wue_by_anomaly, "Q2_plot_01_WUE_by_anomaly_class.png", width = 11, height = 6)

plot_signed_response <- ggplot2::ggplot(
  analysis_data,
  ggplot2::aes(x = anomaly_class, y = WUE_response, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, color = "gray35") +
  ggplot2::geom_boxplot(outlier.alpha = 0.08, width = 0.72) +
  ggplot2::facet_grid(metric ~ timescale, scales = "free_y") +
  ggplot2::labs(
    title = "Signed WUE Responses Relative to Near-Normal Baselines",
    subtitle = "Negative values indicate lower WUE than the near-normal mean for the same metric, timescale, and ecosystem class",
    x = "Moisture anomaly class",
    y = "WUE response",
    fill = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(plot_signed_response, "Q2_plot_02_signed_WUE_response_by_anomaly.png", width = 12, height = 8)

plot_magnitude <- ggplot2::ggplot(
  anomaly_data,
  ggplot2::aes(x = anomaly_class, y = WUE_response_magnitude, fill = water_class)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.08, width = 0.72) +
  ggplot2::facet_grid(metric ~ timescale, scales = "free_y") +
  ggplot2::labs(
    title = "Magnitude of WUE Responses During Dry and Wet Anomalies",
    x = "Moisture anomaly class",
    y = "|WUE response|",
    fill = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(plot_magnitude, "Q2_plot_03_WUE_response_magnitude_by_anomaly.png", width = 12, height = 8)

direction_summary <- anomaly_data |>
  dplyr::group_by(timescale, anomaly_class, metric, water_class) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    pct_decrease = mean(decrease_binary) * 100,
    .groups = "drop"
  )
write_table(direction_summary, "Q2_anomalies_direction_summary.csv")

plot_direction <- ggplot2::ggplot(
  direction_summary,
  ggplot2::aes(x = anomaly_class, y = pct_decrease, fill = water_class)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
  ggplot2::facet_grid(metric ~ timescale, scales = "free_y") +
  ggplot2::labs(
    title = "Probability of WUE Decrease During Moisture Anomalies",
    x = "Moisture anomaly class",
    y = "Months with WUE decrease (%)",
    fill = "Ecosystem class"
  ) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(scale = 1), limits = c(0, 100)) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(plot_direction, "Q2_plot_04_WUE_response_direction_by_anomaly.png", width = 12, height = 8)

plot_spei_gradient <- ggplot2::ggplot(
  analysis_data,
  ggplot2::aes(x = SPEI_value, y = WUE_response, color = water_class)
) +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray35") +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "gray35") +
  ggplot2::geom_point(alpha = 0.18, size = 0.8) +
  ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, linewidth = 0.9) +
  ggplot2::facet_grid(metric ~ timescale, scales = "free") +
  ggplot2::labs(
    title = "WUE Response Along Short- and Long-Term SPEI Gradients",
    subtitle = "Dashed vertical lines mark dry (< -1) and wet (> 1) anomaly thresholds",
    x = "SPEI value",
    y = "WUE response",
    color = "Ecosystem class"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(plot_spei_gradient, "Q2_plot_05_WUE_response_along_SPEI_gradient.png", width = 12, height = 8)

direction_heat <- direction_summary
direction_heat$label <- paste0(round(direction_heat$pct_decrease, 0), "%")

plot_direction_heat <- ggplot2::ggplot(
  direction_heat,
  ggplot2::aes(x = water_class, y = anomaly_class, fill = pct_decrease)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 4) +
  ggplot2::facet_grid(metric ~ timescale) +
  ggplot2::labs(
    title = "Directional WUE Response Summary",
    x = "Ecosystem class",
    y = "Moisture anomaly",
    fill = "Decrease (%)"
  ) +
  ggplot2::scale_fill_gradient(low = "white", high = "#2B6CB0", limits = c(0, 100)) +
  theme_wue() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

plot_wue_by_anomaly_panel <- plot_wue_by_anomaly +
  ggplot2::labs(title = "A. WUE_ET and WUE_T across moisture anomaly classes")

plot_direction_heat_panel <- plot_direction_heat +
  ggplot2::labs(title = "B. Directional WUE response summary")

composite <- plot_wue_by_anomaly_panel / plot_direction_heat_panel +
  patchwork::plot_annotation(
    title = "Q2 Moisture Anomalies and WUE Responses",
    subtitle = "Dry anomalies are SPEI < -1 and wet anomalies are SPEI > 1; models control for T:ET ratio, month, and site-level random effects."
  )

save_plot(composite, "Q2_plot_06_anomaly_model_results_multipanel.png", width = 13, height = 10)

capture_to_file({
  cat("Q2 moisture-anomaly analysis for WUE_ET and WUE_T\n")
  cat("Input:", input_file, "\n")
  cat("Output:", output_dir, "\n\n")
  cat("Rows in source monthly data:", nrow(monthly), "\n")
  cat("Rows after complete-case and ecosystem filters:", nrow(model_base), "\n")
  cat("Rows in metric x timescale long data:", nrow(analysis_data), "\n")
  cat("Rows in dry/wet anomaly subset:", nrow(anomaly_data), "\n")
  cat("Sites:", length(unique(analysis_data$site_name)), "\n")
  cat("Years:", min(analysis_data$Year), "-", max(analysis_data$Year), "\n\n")

  cat("Anomaly definitions:\n")
  cat("  Dry: SPEI < -1\n")
  cat("  Near-normal: -1 <= SPEI <= 1\n")
  cat("  Wet: SPEI > 1\n\n")

  cat("Coverage by timescale/anomaly/metric/ecosystem:\n")
  print(coverage_by_class)

  cat("\nWUE value model:\n")
  print(wue_formula)
  print(summary(wue_model))
  cat("\nWUE value drop1 LRT:\n")
  print(wue_drop1)

  cat("\nSigned response model:\n")
  print(response_formula)
  print(summary(response_model))
  cat("\nSigned response drop1 LRT:\n")
  print(response_drop1)

  cat("\nResponse magnitude model, dry/wet anomalies only:\n")
  print(magnitude_formula)
  print(summary(magnitude_model))
  cat("\nResponse magnitude drop1 LRT:\n")
  print(magnitude_drop1)

  cat("\nResponse direction model, dry/wet anomalies only:\n")
  print(direction_formula)
  print(summary(direction_model))
  cat("\nResponse direction drop1 LRT:\n")
  print(direction_drop1)
}, "Q2_anomalies_model_summary.txt")

message("Done. Q2 anomaly outputs written to: ", output_dir)
