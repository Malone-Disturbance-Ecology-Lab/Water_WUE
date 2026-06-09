# Simplified Q1 WUE metric-difference models
#
# Goal:
#   Ask what the Q1 results look like when the model includes metric
#   difference type, T:ET ratio, ecosystem type, and month/site random effects.
#
# Models:
#   Full-data mixed model:
#     WUE_difference ~ difference_type * water_class * Trans_ratio_z +
#       (1 | site_name) + (1 | month_f)
#
#   Near-normal GAM:
#     WUE_difference ~ difference_type * water_class +
#       s(Trans_ratio, by = difference_type:water_class) +
#       s(site_name, bs = "re") + s(month_f, bs = "re")
#
# Outputs:
#   Malone_Workflow/results/Q1_WUE_metric_difference_simple/

suppressPackageStartupMessages({
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required. Install it with install.packages('lme4').")
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required. Install it with install.packages('mgcv').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required. Install it with install.packages('scales').")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install it with install.packages('patchwork').")
  }
})

sig <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*",
          ifelse(p < 0.1, ".", "")
        )
      )
    )
  )
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

theme_wue <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

theme_wue_manuscript <- function() {
  ggplot2::theme_minimal(base_size = 27) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 25),
      plot.title = ggplot2::element_text(face = "bold", size = 33),
      plot.subtitle = ggplot2::element_text(size = 25),
      axis.title = ggplot2::element_text(size = 27),
      axis.text = ggplot2::element_text(size = 23),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 26),
      legend.text = ggplot2::element_text(size = 24)
    )
}

project_dir <- normalizePath(file.path(dirname(getwd()), "WUE"), mustWork = FALSE)
if (!dir.exists(file.path(project_dir, "Water_WUE"))) {
  project_dir <- normalizePath(getwd(), mustWork = TRUE)
}

input_file <- file.path(project_dir, "Water_WUE", "data", "WUE_CUE_monthly_merged_indices_clean.csv")
output_dir <- file.path(
  project_dir,
  "Water_WUE",
  "Malone_Workflow",
  "results",
  "Q1_WUE_metric_difference_simple"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
monthly[] <- lapply(monthly, function(x) if (is.character(x)) trimws(x) else x)

required_cols <- c(
  "site_name", "Year", "month", "water_class", "Trans_ratio",
  "WUE", "WUE_tra", "WUE_eva", "SPEI_1"
)
missing_cols <- setdiff(required_cols, names(monthly))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ecosystem_classes <- c("Upland", "Freshwater", "Saline")
ecosystem_colors <- c(Upland = "#800080", Freshwater = "#0000FF", Saline = "#FFA500")

model_base <- monthly[complete.cases(monthly[, required_cols]), ]
model_base <- model_base[
  model_base$site_name != "" &
    model_base$water_class != "" &
    model_base$water_class %in% ecosystem_classes &
    is.finite(model_base$Trans_ratio) &
    model_base$Trans_ratio >= 0 &
    model_base$Trans_ratio <= 1 &
    is.finite(model_base$WUE) &
    is.finite(model_base$WUE_tra) &
    is.finite(model_base$WUE_eva),
]

model_base$WUE_ET_minus_WUE_T <- model_base$WUE - model_base$WUE_tra
model_base$WUE_ET_minus_WUE_E <- model_base$WUE - model_base$WUE_eva
model_base$WUE_E_minus_WUE_T <- model_base$WUE_eva - model_base$WUE_tra

difference_cols <- c(
  "WUE_ET_minus_WUE_T",
  "WUE_ET_minus_WUE_E",
  "WUE_E_minus_WUE_T"
)

diff_data <- reshape(
  model_base,
  varying = difference_cols,
  v.names = "WUE_difference",
  timevar = "difference_type",
  times = difference_cols,
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)
diff_data <- diff_data[
  complete.cases(diff_data[, c("site_name", "Year", "month", "Trans_ratio", "WUE_difference")]),
]
diff_data <- diff_data[is.finite(diff_data$WUE_difference), ]
diff_data$site_name <- factor(diff_data$site_name)
diff_data$month_f <- factor(diff_data$month, levels = 1:12)
diff_data$water_class <- factor(diff_data$water_class, levels = ecosystem_classes)
diff_data$difference_type <- factor(diff_data$difference_type, levels = difference_cols)
diff_data$Trans_ratio_z <- as.numeric(scale(diff_data$Trans_ratio))
diff_data$diff_ecosystem <- interaction(
  diff_data$difference_type,
  diff_data$water_class,
  sep = "__",
  drop = TRUE
)

diff_labels <- c(
  WUE_ET_minus_WUE_T = "WUE_ET - WUE_T",
  WUE_ET_minus_WUE_E = "WUE_ET - WUE_E",
  WUE_E_minus_WUE_T = "WUE_E - WUE_T"
)
diff_label_levels <- unname(diff_labels[difference_cols])
diff_data$difference_label <- factor(
  diff_labels[as.character(diff_data$difference_type)],
  levels = diff_label_levels
)

# ---------------------------------------------------------------------------
# Full-data simple mixed model
# ---------------------------------------------------------------------------

simple_mixed <- lme4::lmer(
  WUE_difference ~ difference_type * water_class * Trans_ratio_z +
    (1 | site_name) + (1 | month_f),
  data = diff_data,
  REML = FALSE
)

simple_mixed_additive <- lme4::lmer(
  WUE_difference ~ difference_type + water_class + Trans_ratio_z +
    (1 | site_name) + (1 | month_f),
  data = diff_data,
  REML = FALSE
)

saveRDS(simple_mixed, file.path(output_dir, "Q1_simple_mixed_difference_TET_model.rds"))

diff_data$simple_mixed_fitted <- fitted(simple_mixed)
diff_data$simple_mixed_residual <- residuals(simple_mixed)
write.csv(
  diff_data[, c(
    "site_name", "Year", "month", "water_class", "difference_type",
    "Trans_ratio", "Trans_ratio_z", "WUE_difference",
    "simple_mixed_fitted", "simple_mixed_residual"
  )],
  file.path(output_dir, "Q1_simple_mixed_fitted_residuals.csv"),
  row.names = FALSE
)

mixed_fixed <- as.data.frame(summary(simple_mixed)$coefficients)
mixed_fixed$term <- rownames(mixed_fixed)
rownames(mixed_fixed) <- NULL
mixed_fixed <- mixed_fixed[, c("term", setdiff(names(mixed_fixed), "term"))]
names(mixed_fixed) <- gsub(" ", "_", names(mixed_fixed), fixed = TRUE)
mixed_fixed$p_value_normal_approx <- 2 * stats::pnorm(-abs(mixed_fixed$t_value))
mixed_fixed$signif <- sig(mixed_fixed$p_value_normal_approx)
write.csv(
  mixed_fixed,
  file.path(output_dir, "Q1_simple_mixed_fixed_effects_with_approx_p.csv"),
  row.names = FALSE
)

mixed_lrt <- as.data.frame(drop1(simple_mixed, test = "Chisq"))
mixed_lrt$term <- rownames(mixed_lrt)
rownames(mixed_lrt) <- NULL
mixed_lrt <- mixed_lrt[, c("term", setdiff(names(mixed_lrt), "term"))]
names(mixed_lrt) <- gsub(" ", "_", names(mixed_lrt), fixed = TRUE)
p_col <- grep("Pr", names(mixed_lrt), value = TRUE)[1]
mixed_lrt$signif <- sig(mixed_lrt[[p_col]])
write.csv(
  mixed_lrt,
  file.path(output_dir, "Q1_simple_mixed_likelihood_ratio_tests.csv"),
  row.names = FALSE
)

mixed_comparison <- as.data.frame(anova(simple_mixed_additive, simple_mixed))
p_col <- grep("Pr", names(mixed_comparison), value = TRUE)[1]
mixed_comparison$signif <- sig(mixed_comparison[[p_col]])
write.csv(
  mixed_comparison,
  file.path(output_dir, "Q1_simple_mixed_model_comparison_lrt.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Near-normal simple GAM
# ---------------------------------------------------------------------------

gam_data <- diff_data[
  is.finite(diff_data$SPEI_1) &
    diff_data$SPEI_1 >= -1 &
    diff_data$SPEI_1 <= 1,
]

simple_linear_gam <- mgcv::gam(
  WUE_difference ~ difference_type * water_class * Trans_ratio_z +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = gam_data,
  method = "ML"
)

simple_smooth_gam <- mgcv::gam(
  WUE_difference ~ difference_type * water_class +
    s(Trans_ratio, by = diff_ecosystem, k = 6) +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = gam_data,
  method = "ML",
  select = TRUE
)

saveRDS(simple_smooth_gam, file.path(output_dir, "Q1_simple_smooth_gam_difference_TET_model.rds"))

gam_data$simple_gam_fitted <- fitted(simple_smooth_gam)
gam_data$simple_gam_residual <- residuals(simple_smooth_gam)
write.csv(
  gam_data[, c(
    "site_name", "Year", "month", "water_class", "difference_type",
    "Trans_ratio", "Trans_ratio_z", "SPEI_1", "WUE_difference",
    "simple_gam_fitted", "simple_gam_residual"
  )],
  file.path(output_dir, "Q1_simple_gam_fitted_residuals.csv"),
  row.names = FALSE
)

gam_comparison <- data.frame(
  model = c("simple_linear_gam", "simple_smooth_gam"),
  AIC = c(AIC(simple_linear_gam), AIC(simple_smooth_gam)),
  BIC = c(BIC(simple_linear_gam), BIC(simple_smooth_gam)),
  deviance_explained = c(summary(simple_linear_gam)$dev.expl, summary(simple_smooth_gam)$dev.expl),
  adjusted_r_squared = c(summary(simple_linear_gam)$r.sq, summary(simple_smooth_gam)$r.sq),
  scale = c(summary(simple_linear_gam)$scale, summary(simple_smooth_gam)$scale)
)
gam_comparison$delta_AIC <- gam_comparison$AIC - min(gam_comparison$AIC)
write.csv(
  gam_comparison,
  file.path(output_dir, "Q1_simple_gam_model_comparison.csv"),
  row.names = FALSE
)

gam_param <- as.data.frame(summary(simple_smooth_gam)$p.table)
gam_param$term <- rownames(gam_param)
rownames(gam_param) <- NULL
gam_param <- gam_param[, c("term", setdiff(names(gam_param), "term"))]
names(gam_param) <- gsub(" ", "_", names(gam_param), fixed = TRUE)
p_col <- grep("Pr", names(gam_param), value = TRUE)[1]
names(gam_param)[names(gam_param) == p_col] <- "p_value"
gam_param$signif <- sig(gam_param$p_value)
write.csv(
  gam_param,
  file.path(output_dir, "Q1_simple_gam_parametric_terms.csv"),
  row.names = FALSE
)

gam_smooth <- as.data.frame(summary(simple_smooth_gam)$s.table)
gam_smooth$smooth_term <- rownames(gam_smooth)
rownames(gam_smooth) <- NULL
gam_smooth <- gam_smooth[, c("smooth_term", setdiff(names(gam_smooth), "smooth_term"))]
names(gam_smooth) <- gsub("p-value", "p_value", names(gam_smooth), fixed = TRUE)
gam_smooth$signif <- sig(gam_smooth$p_value)
write.csv(
  gam_smooth,
  file.path(output_dir, "Q1_simple_gam_smooth_terms.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Predictions and plots
# ---------------------------------------------------------------------------

tet_sequence <- seq(
  stats::quantile(diff_data$Trans_ratio, 0.05),
  stats::quantile(diff_data$Trans_ratio, 0.95),
  length.out = 100
)
tet_mean <- mean(diff_data$Trans_ratio)
tet_sd <- stats::sd(diff_data$Trans_ratio)

mixed_grid <- expand.grid(
  difference_type = levels(diff_data$difference_type),
  water_class = levels(diff_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(diff_data$site_name)[1],
  month_f = levels(diff_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
mixed_grid$difference_type <- factor(mixed_grid$difference_type, levels = levels(diff_data$difference_type))
mixed_grid$water_class <- factor(mixed_grid$water_class, levels = levels(diff_data$water_class))
mixed_grid$site_name <- factor(mixed_grid$site_name, levels = levels(diff_data$site_name))
mixed_grid$month_f <- factor(mixed_grid$month_f, levels = levels(diff_data$month_f))
mixed_grid$Trans_ratio_z <- (mixed_grid$Trans_ratio - tet_mean) / tet_sd
mixed_grid$difference_label <- factor(
  diff_labels[as.character(mixed_grid$difference_type)],
  levels = diff_label_levels
)
mixed_grid$ecosystem_label <- as.character(mixed_grid$water_class)
mixed_grid$prediction <- predict(simple_mixed, newdata = mixed_grid, re.form = NA)

fixed_terms <- stats::delete.response(stats::terms(suppressWarnings(lme4::nobars(formula(simple_mixed)))))
fixed_model_matrix <- stats::model.matrix(fixed_terms, mixed_grid)
fixed_vcov <- as.matrix(stats::vcov(simple_mixed))
fixed_model_matrix <- fixed_model_matrix[, colnames(fixed_vcov), drop = FALSE]
mixed_grid$se <- sqrt(diag(fixed_model_matrix %*% fixed_vcov %*% t(fixed_model_matrix)))
mixed_grid$lower <- mixed_grid$prediction - 1.96 * mixed_grid$se
mixed_grid$upper <- mixed_grid$prediction + 1.96 * mixed_grid$se
write.csv(
  mixed_grid,
  file.path(output_dir, "Q1_simple_mixed_predictions_TET.csv"),
  row.names = FALSE
)

gam_grid <- expand.grid(
  difference_type = levels(gam_data$difference_type),
  water_class = levels(gam_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(gam_data$site_name)[1],
  month_f = levels(gam_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
gam_grid$difference_type <- factor(gam_grid$difference_type, levels = levels(gam_data$difference_type))
gam_grid$water_class <- factor(gam_grid$water_class, levels = levels(gam_data$water_class))
gam_grid$site_name <- factor(gam_grid$site_name, levels = levels(gam_data$site_name))
gam_grid$month_f <- factor(gam_grid$month_f, levels = levels(gam_data$month_f))
gam_grid$Trans_ratio_z <- (gam_grid$Trans_ratio - mean(gam_data$Trans_ratio)) / stats::sd(gam_data$Trans_ratio)
gam_grid$difference_label <- factor(
  diff_labels[as.character(gam_grid$difference_type)],
  levels = diff_label_levels
)
gam_grid$ecosystem_label <- as.character(gam_grid$water_class)
gam_grid$diff_ecosystem <- interaction(
  gam_grid$difference_type,
  gam_grid$water_class,
  sep = "__",
  drop = TRUE
)
gam_grid$diff_ecosystem <- factor(gam_grid$diff_ecosystem, levels = levels(gam_data$diff_ecosystem))

gam_pred <- predict(
  simple_smooth_gam,
  newdata = gam_grid,
  type = "link",
  se.fit = TRUE,
  exclude = c("s(site_name)", "s(month_f)")
)
gam_grid$prediction <- as.numeric(gam_pred$fit)
gam_grid$se <- as.numeric(gam_pred$se.fit)
gam_grid$lower <- gam_grid$prediction - 1.96 * gam_grid$se
gam_grid$upper <- gam_grid$prediction + 1.96 * gam_grid$se
write.csv(
  gam_grid,
  file.path(output_dir, "Q1_simple_gam_predictions_TET.csv"),
  row.names = FALSE
)

summary_by_difference <- diff_data |>
  dplyr::group_by(difference_type, difference_label, water_class) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sites = dplyr::n_distinct(site_name),
    mean_difference = mean(WUE_difference),
    median_difference = median(WUE_difference),
    sd_difference = sd(WUE_difference),
    se_difference = sd_difference / sqrt(n_months),
    ci95_difference = stats::qt(0.975, df = n_months - 1) * se_difference,
    .groups = "drop"
  )
write.csv(
  summary_by_difference,
  file.path(output_dir, "Q1_simple_summary_by_difference_type.csv"),
  row.names = FALSE
)

plot_summary <- ggplot2::ggplot(
  summary_by_difference,
  ggplot2::aes(x = water_class, y = mean_difference, color = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_difference - ci95_difference, ymax = mean_difference + ci95_difference),
    width = 0.15,
    linewidth = 0.8
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~ difference_label, scales = "free_y") +
  ggplot2::labs(
    title = "Q1 Model Inputs: Mean WUE Metric Differences by Ecosystem",
    subtitle = "Points are means and bars are 95% CIs",
    x = "Ecosystem",
    y = "Mean WUE metric difference",
    color = "Ecosystem"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors, guide = "none") +
  theme_wue() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
save_plot(plot_summary, "Q1_simple_plot_01_mean_difference_type.png", width = 10, height = 5)

plot_mixed <- ggplot2::ggplot(
  mixed_grid,
  ggplot2::aes(x = Trans_ratio, y = prediction, color = water_class, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(~ difference_label, scales = "free_y") +
  ggplot2::labs(
    title = "Mixed Model: T:ET Response by Ecosystem",
    subtitle = "Fixed-effect predictions; site and month random effects excluded",
    x = "T:ET ratio (Trans_ratio)",
    y = "Predicted WUE metric difference",
    color = "Ecosystem",
    fill = "Ecosystem"
  ) +
  ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()
save_plot(plot_mixed, "Q1_simple_plot_02_mixed_TET_predictions.png", width = 10, height = 6)

plot_gam <- ggplot2::ggplot(
  gam_grid,
  ggplot2::aes(x = Trans_ratio, y = prediction, color = water_class, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::facet_wrap(~ difference_label, scales = "free_y") +
  ggplot2::labs(
    title = "Smooth GAM: T:ET Response by Ecosystem",
    subtitle = "Near-normal SPEI_1 subset; site and month random effects excluded",
    x = "T:ET ratio (Trans_ratio)",
    y = "Predicted WUE metric difference",
    color = "Ecosystem",
    fill = "Ecosystem"
  ) +
  ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue()
save_plot(plot_gam, "Q1_simple_plot_03_gam_TET_predictions.png", width = 10, height = 6)

model_fit_long <- data.frame(
  model = factor(
    c("Mixed model", "Linear GAM", "Smooth GAM"),
    levels = c("Mixed model", "Linear GAM", "Smooth GAM")
  ),
  deviance_explained = c(
    NA_real_,
    gam_comparison$deviance_explained[gam_comparison$model == "simple_linear_gam"] * 100,
    gam_comparison$deviance_explained[gam_comparison$model == "simple_smooth_gam"] * 100
  ),
  AIC = c(AIC(simple_mixed), gam_comparison$AIC)
)

plot_model_comparison <- ggplot2::ggplot(
  model_fit_long[!is.na(model_fit_long$deviance_explained), ],
  ggplot2::aes(x = model, y = deviance_explained, fill = model)
) +
  ggplot2::geom_col(width = 0.65, color = "gray30", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(round(deviance_explained, 1), "%")),
    vjust = -0.4,
    size = 6
  ) +
  ggplot2::labs(
    title = "Simple Q1 GAM Model Comparison",
    subtitle = paste0("Smooth GAM improves AIC by ", round(gam_comparison$delta_AIC[gam_comparison$model == "simple_linear_gam"], 1)),
    x = NULL,
    y = "Deviance explained (%)"
  ) +
  ggplot2::scale_fill_manual(values = c("Linear GAM" = "#8DA0CB", "Smooth GAM" = "#66C2A5"), guide = "none") +
  ggplot2::ylim(0, max(model_fit_long$deviance_explained, na.rm = TRUE) + 10) +
  theme_wue()
save_plot(plot_model_comparison, "Q1_simple_plot_04_gam_model_comparison.png", width = 7, height = 5)

plot_mixed_residuals <- ggplot2::ggplot(
  diff_data,
  ggplot2::aes(x = simple_mixed_fitted, y = simple_mixed_residual, color = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  ggplot2::geom_point(alpha = 0.25, size = 1) +
  ggplot2::facet_wrap(~ difference_label, scales = "free_x") +
  ggplot2::labs(
    title = "Simple Mixed Model Residuals",
    x = "Fitted difference",
    y = "Residual",
    color = "Difference type"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors, guide = "none") +
  theme_wue()
save_plot(plot_mixed_residuals, "Q1_simple_plot_05_mixed_residuals.png", width = 10, height = 6)

plot_gam_residuals <- ggplot2::ggplot(
  gam_data,
  ggplot2::aes(x = simple_gam_fitted, y = simple_gam_residual, color = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.4) +
  ggplot2::geom_point(alpha = 0.25, size = 1) +
  ggplot2::facet_wrap(~ difference_label, scales = "free_x") +
  ggplot2::labs(
    title = "Simple Smooth GAM Residuals",
    x = "Fitted difference",
    y = "Residual",
    color = "Difference type"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors, guide = "none") +
  theme_wue()
save_plot(plot_gam_residuals, "Q1_simple_plot_06_gam_residuals.png", width = 10, height = 6)

panel_summary <- plot_summary +
  ggplot2::labs(tag = "A", title = "Mean Metric Differences by Ecosystem", subtitle = NULL) +
  theme_wue_manuscript() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 23))

panel_model_comparison <- plot_model_comparison +
  ggplot2::labs(tag = "B", title = "GAM Model Comparison", subtitle = "Smooth GAM improves AIC by 638.7") +
  theme_wue_manuscript()

panel_mixed <- plot_mixed +
  ggplot2::labs(tag = "C", title = "Mixed Model T:ET Response", subtitle = NULL) +
  theme_wue_manuscript()

panel_gam <- plot_gam +
  ggplot2::labs(tag = "D", title = "Smooth GAM T:ET Response", subtitle = NULL) +
  theme_wue_manuscript()

manuscript_panel <- (
  panel_summary
) / (
  panel_model_comparison
) / (
  panel_mixed
) / (
  panel_gam
) +
  patchwork::plot_layout(heights = c(0.9, 0.8, 1.25, 1.25), guides = "collect") +
  patchwork::plot_annotation(
    title = "Q1 WUE Metric Difference Models",
    subtitle = "Difference type, T:ET ratio, ecosystem fixed effects, and site/month random effects"
  ) &
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 42),
    plot.subtitle = ggplot2::element_text(size = 32),
    plot.tag = ggplot2::element_text(face = "bold", size = 38),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )
save_plot(manuscript_panel, "Q1_simple_plot_07_manuscript_multipanel.png", width = 24, height = 25)

readme <- file.path(output_dir, "Q1_simple_outputs_README.md")
con <- file(readme, "w")
writeLines(c(
  "# Q1 difference type + T:ET + ecosystem outputs",
  "",
  "This is the retained streamlined Q1 analysis. It includes ecosystem class as a fixed effect and treats site and month as random effects.",
  "",
  "## Core models",
  "- Mixed model: `WUE_difference ~ difference_type * water_class * Trans_ratio_z + (1 | site_name) + (1 | month_f)`",
  "- Near-normal smooth GAM: `WUE_difference ~ difference_type * water_class + s(Trans_ratio, by = difference_type:water_class) + s(site_name, bs = \"re\") + s(month_f, bs = \"re\")`",
  "- Ecosystem colors match the Malone_Workflow figures: Upland purple, Freshwater blue, Saline orange.",
  "",
  "## Key tables",
  "- `Q1_simple_summary_by_difference_type.csv`",
  "- `Q1_simple_mixed_fixed_effects_with_approx_p.csv`",
  "- `Q1_simple_mixed_likelihood_ratio_tests.csv`",
  "- `Q1_simple_mixed_model_comparison_lrt.csv`",
  "- `Q1_simple_gam_parametric_terms.csv`",
  "- `Q1_simple_gam_smooth_terms.csv`",
  "- `Q1_simple_gam_model_comparison.csv`",
  "- `Q1_simple_mixed_predictions_TET.csv`",
  "- `Q1_simple_gam_predictions_TET.csv`",
  "- `Q1_simple_mixed_fitted_residuals.csv`",
  "- `Q1_simple_gam_fitted_residuals.csv`",
  "",
  "## Figures",
  "- `Q1_simple_plot_01_mean_difference_type.png`",
  "- `Q1_simple_plot_02_mixed_TET_predictions.png`",
  "- `Q1_simple_plot_03_gam_TET_predictions.png`",
  "- `Q1_simple_plot_04_gam_model_comparison.png`",
  "- `Q1_simple_plot_05_mixed_residuals.png`",
  "- `Q1_simple_plot_06_gam_residuals.png`",
  "- `Q1_simple_plot_07_manuscript_multipanel.png`"
), con)
close(con)

sink(file.path(output_dir, "Q1_simple_model_summary.txt"))
cat("Q1 difference type + T:ET + ecosystem models\n")
cat("Input:", input_file, "\n")
cat("Output:", output_dir, "\n\n")
cat("Full-data rows:", nrow(diff_data), "\n")
cat("Near-normal GAM rows:", nrow(gam_data), "\n")
cat("Sites:", length(unique(diff_data$site_name)), "\n\n")
cat("Ecosystem classes:", paste(levels(diff_data$water_class), collapse = ", "), "\n\n")
cat("Mixed model formula:\n")
print(formula(simple_mixed))
cat("\nMixed model summary:\n")
print(summary(simple_mixed))
cat("\nMixed model drop1 LRT:\n")
print(mixed_lrt)
cat("\nMixed additive vs interaction comparison:\n")
print(mixed_comparison)
cat("\nSimple GAM formula:\n")
print(formula(simple_smooth_gam))
cat("\nSimple GAM summary:\n")
print(summary(simple_smooth_gam))
cat("\nGAM comparison:\n")
print(gam_comparison)
sink()

message("Done. Simple Q1 outputs written to: ", output_dir)
