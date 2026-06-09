# Q2_WUE_performance.R
#
# Goal:
#   Characterise site-level WUE_T performance across five complementary
#   dimensions that together describe how each ecosystem navigates drought:
#
#   1. Stability index  = mean(WUE_T) / var(WUE_T)          (mean-to-variance ratio; higher = more stable)
#   2. Plasticity_slope = slope of WUE_T ~ SPEI              (per SPEI timescale; site-level OLS)
#   3. Plasticity_range = max(WUE_T) / min(WUE_T)            (dynamic range ratio)
#   4. Resistance       = mean WUE_T during drought / mean WUE_T during near-normal conditions
#   5. Recovery         = mean WUE_T in 12-months post-drought / mean WUE_T in 12-months pre-drought
#
# Drought is defined as SPEI_3 < -1.  Near-normal is -1 <= SPEI_3 <= 1.
# A drought event is a run of >= 2 consecutive months with SPEI_3 < -1.
# Pre- and post-drought windows are the 12 months immediately before/after
# each event; sites need at least 3 usable months in each window.
#
# Response variable:
#   WUE_T = WUE_tra
#
# SPEI timescales used for plasticity slopes:
#   SPEI_1, SPEI_3, SPEI_6, SPEI_12, SPEI_24, SPEI_36, SPEI_48
#
# Input:
#   Water_WUE/data/WUE_CUE_monthly_merged_indices_clean.csv
#
# Outputs:
#   Water_WUE/data/results/Q2_WUE_performance/
#     Q2_site_stability.csv
#     Q2_site_plasticity_slope.csv
#     Q2_site_plasticity_range.csv
#     Q2_site_resistance.csv
#     Q2_site_recovery.csv
#     Q2_site_performance_summary.csv
#     Q2_plot_01_stability.png
#     Q2_plot_02_plasticity_slope.png
#     Q2_plot_03_plasticity_range.png
#     Q2_plot_04_resistance.png
#     Q2_plot_05_recovery.png
#     Q2_plot_06_composite.png
#     Q2_WUE_performance_summary.txt

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "tidyr", "ggplot2", "patchwork", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install with install.packages('%s').", pkg, pkg))
    }
  }
})

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  }
  normalizePath("Malone_Workflow/Q2_WUE_performance.R", mustWork = FALSE)
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

ecosystem_metric_test <- function(df, metric_col, metric_label) {
  test_df <- df[is.finite(df[[metric_col]]) & !is.na(df$water_class), c("water_class", metric_col)]
  names(test_df)[2] <- "value"
  test_df$water_class <- factor(test_df$water_class, levels = ecosystem_classes)
  test_df <- test_df[!is.na(test_df$water_class), ]

  group_counts <- as.data.frame(table(test_df$water_class), stringsAsFactors = FALSE)
  names(group_counts) <- c("water_class", "n")
  get_n <- function(level) {
    out <- group_counts$n[group_counts$water_class == level]
    ifelse(length(out) == 0, 0L, as.integer(out))
  }

  result <- data.frame(
    metric = metric_label,
    response_variable = metric_col,
    test = "Kruskal-Wallis",
    n_sites = nrow(test_df),
    n_upland = get_n("Upland"),
    n_freshwater = get_n("Freshwater"),
    n_saline = get_n("Saline"),
    statistic = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    p_adjust_method = "BH/FDR for pairwise Wilcoxon",
    pairwise_upland_vs_freshwater_p_adj = NA_real_,
    pairwise_upland_vs_saline_p_adj = NA_real_,
    pairwise_freshwater_vs_saline_p_adj = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(test_df) >= 3 && length(unique(test_df$water_class)) >= 2) {
    kw <- stats::kruskal.test(value ~ water_class, data = test_df)
    result$statistic <- as.numeric(kw$statistic)
    result$df <- as.numeric(kw$parameter)
    result$p_value <- as.numeric(kw$p.value)

    pw <- stats::pairwise.wilcox.test(
      x = test_df$value,
      g = test_df$water_class,
      p.adjust.method = "BH",
      exact = FALSE
    )
    pw_mat <- pw$p.value
    get_pair_p <- function(row_name, col_name) {
      if (row_name %in% rownames(pw_mat) && col_name %in% colnames(pw_mat)) {
        return(as.numeric(pw_mat[row_name, col_name]))
      }
      if (col_name %in% rownames(pw_mat) && row_name %in% colnames(pw_mat)) {
        return(as.numeric(pw_mat[col_name, row_name]))
      }
      NA_real_
    }
    result$pairwise_upland_vs_freshwater_p_adj <- get_pair_p("Freshwater", "Upland")
    result$pairwise_upland_vs_saline_p_adj <- get_pair_p("Saline", "Upland")
    result$pairwise_freshwater_vs_saline_p_adj <- get_pair_p("Saline", "Freshwater")
  }

  result$signif <- sig(result$p_value)
  result
}

theme_wue <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "bottom"
    )
}

# Identify runs of consecutive months meeting a condition.
# Returns a data frame with: site_name, event_id, start_idx, end_idx,
# start_year, start_month, end_year, end_month, duration_months.
identify_drought_events <- function(df, spei_col = "SPEI_3", threshold = -1, min_duration = 2) {
  df <- df[order(df$site_name, df$Year, df$month), ]
  events <- do.call(rbind, lapply(split(df, df$site_name), function(site_df) {
    site_df <- site_df[order(site_df$Year, site_df$month), ]
    n <- nrow(site_df)
    in_drought <- site_df[[spei_col]] < threshold & is.finite(site_df[[spei_col]])
    # find run-length encoding
    rle_out <- rle(in_drought)
    ends     <- cumsum(rle_out$lengths)
    starts   <- c(1L, ends[-length(ends)] + 1L)
    event_rows <- data.frame(
      value  = rle_out$values,
      start  = starts,
      end    = ends,
      length = rle_out$lengths
    )
    drought_events <- event_rows[event_rows$value & event_rows$length >= min_duration, ]
    if (nrow(drought_events) == 0) return(NULL)
    data.frame(
      site_name     = site_df$site_name[1],
      event_id      = seq_len(nrow(drought_events)),
      start_idx     = drought_events$start,
      end_idx       = drought_events$end,
      start_year    = site_df$Year[drought_events$start],
      start_month   = site_df$month[drought_events$start],
      end_year      = site_df$Year[drought_events$end],
      end_month     = site_df$month[drought_events$end],
      duration_months = drought_events$length
    )
  }))
  if (is.null(events)) return(data.frame())
  rownames(events) <- NULL
  events
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

script_dir  <- dirname(get_script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

input_file  <- file.path(project_dir, "data", "WUE_CUE_monthly_merged_indices_clean.csv")
output_dir  <- file.path(script_dir, "results", "Q2_WUE_performance")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load and clean data
# ---------------------------------------------------------------------------

spei_cols    <- c("SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48")
single_spei_col <- "SPEI_3"
single_spei_label <- "SPEI-3"
required_cols <- c("site_name", "Year", "month", "water_class", "lat", "long",
                   "WUE_tra", "Trans_ratio", spei_cols)

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
monthly[] <- lapply(monthly, function(x) if (is.character(x)) trimws(x) else x)

missing_cols <- setdiff(required_cols, names(monthly))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ecosystem_classes <- c("Upland", "Freshwater", "Saline")

base_df <- monthly[
  complete.cases(monthly[, c("site_name", "Year", "month", "water_class",
                              "lat", "long", "WUE_tra", "Trans_ratio")]) &
  monthly$site_name  != "" &
  monthly$water_class != "" &
  monthly$water_class %in% ecosystem_classes &
  is.finite(monthly$lat) &
  is.finite(monthly$long) &
  is.finite(monthly$WUE_tra) &
  is.finite(monthly$Trans_ratio) &
  monthly$Trans_ratio >= 0 &
  monthly$Trans_ratio <= 1,
]

base_df$WUE_T       <- base_df$WUE_tra
base_df$site_name   <- as.character(base_df$site_name)
base_df$water_class <- factor(base_df$water_class, levels = ecosystem_classes)

# Ordered time index within each site (used for recovery windowing)
base_df <- base_df[order(base_df$site_name, base_df$Year, base_df$month), ]
base_df$row_idx <- ave(
  seq_len(nrow(base_df)),
  base_df$site_name,
  FUN = seq_along
)

message(sprintf(
  "Data loaded: %d rows, %d sites, years %d-%d",
  nrow(base_df),
  length(unique(base_df$site_name)),
  min(base_df$Year), max(base_df$Year)
))

# ---------------------------------------------------------------------------
# 1. Stability index  =  mean(WUE_T) / var(WUE_T)
# ---------------------------------------------------------------------------

stability <- base_df |>
  dplyr::group_by(site_name, water_class) |>
  dplyr::filter(dplyr::n() >= 6) |>
  dplyr::summarise(
    n_months      = dplyr::n(),
    mean_TET      = mean(Trans_ratio),
    mean_WUE_T    = mean(WUE_T),
    var_WUE_T     = var(WUE_T),
    sd_WUE_T      = sd(WUE_T),
    cv_WUE_T      = sd_WUE_T / mean_WUE_T,          # coefficient of variation
    stability     = mean_WUE_T / var_WUE_T,           # signal-to-variance ratio
    .groups = "drop"
  )

write_table(stability, "Q2_site_stability.csv")
message("Stability: ", nrow(stability), " sites")

# ---------------------------------------------------------------------------
# 2. Plasticity — slope  =  OLS slope of WUE_T ~ SPEI per site per timescale
# ---------------------------------------------------------------------------

plasticity_slope <- do.call(rbind, lapply(spei_cols, function(sc) {
  sub_df <- base_df[is.finite(base_df[[sc]]), c("site_name", "water_class", "WUE_T", sc)]
  names(sub_df)[4] <- "SPEI_value"
  sub_df |>
    dplyr::group_by(site_name, water_class) |>
    dplyr::filter(dplyr::n_distinct(SPEI_value) >= 4, dplyr::n() >= 6) |>
    dplyr::summarise(
      SPEI_timescale = sc,
      n_months       = dplyr::n(),
      slope          = stats::coef(stats::lm(WUE_T ~ SPEI_value))[["SPEI_value"]],
      r_squared      = summary(stats::lm(WUE_T ~ SPEI_value))$r.squared,
      correlation    = stats::cor(WUE_T, SPEI_value),
      .groups = "drop"
    )
}))

plasticity_slope$SPEI_timescale <- factor(plasticity_slope$SPEI_timescale, levels = spei_cols)
write_table(plasticity_slope, "Q2_site_plasticity_slope.csv")
message("Plasticity slope: ", nrow(plasticity_slope), " site × timescale rows")

# ---------------------------------------------------------------------------
# 3. Plasticity — range  =  max(WUE_T) / min(WUE_T) per site
# ---------------------------------------------------------------------------

plasticity_range <- base_df |>
  dplyr::group_by(site_name, water_class) |>
  dplyr::filter(dplyr::n() >= 6, min(WUE_T) > 0) |>
  dplyr::summarise(
    n_months        = dplyr::n(),
    WUE_T_max       = max(WUE_T),
    WUE_T_min       = min(WUE_T),
    plasticity_range = WUE_T_max / WUE_T_min,
    WUE_T_p95       = stats::quantile(WUE_T, 0.95),
    WUE_T_p05       = stats::quantile(WUE_T, 0.05),
    plasticity_p95_p05 = WUE_T_p95 / WUE_T_p05,  # robust version using 5th–95th percentile
    .groups = "drop"
  )

write_table(plasticity_range, "Q2_site_plasticity_range.csv")
message("Plasticity range: ", nrow(plasticity_range), " sites")

# ---------------------------------------------------------------------------
# 4. Resistance  =  mean WUE_T(drought) / mean WUE_T(near-normal)
#    Drought threshold: SPEI_3 < -1
#    Near-normal:       -1 <= SPEI_3 <= 1
# ---------------------------------------------------------------------------

resistance_df <- base_df[is.finite(base_df[[single_spei_col]]), ]
resistance_df$drought_class <- ifelse(
  resistance_df[[single_spei_col]] < -1,   "drought",
  ifelse(resistance_df[[single_spei_col]] <= 1, "near_normal", "wet")
)

# site-level means by drought class
resistance_wide <- resistance_df |>
  dplyr::filter(drought_class %in% c("drought", "near_normal")) |>
  dplyr::group_by(site_name, water_class, drought_class) |>
  dplyr::summarise(
    n_months    = dplyr::n(),
    mean_WUE_T  = mean(WUE_T),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    names_from   = drought_class,
    values_from  = c(mean_WUE_T, n_months),
    names_sep    = "_"
  )

# require at least 3 months in each class
resistance <- resistance_wide[
  !is.na(resistance_wide$mean_WUE_T_drought) &
  !is.na(resistance_wide$mean_WUE_T_near_normal) &
  resistance_wide$n_months_drought     >= 3 &
  resistance_wide$n_months_near_normal >= 3 &
  resistance_wide$mean_WUE_T_near_normal > 0,
]
resistance$resistance <- resistance$mean_WUE_T_drought / resistance$mean_WUE_T_near_normal

write_table(resistance, "Q2_site_resistance.csv")
message("Resistance: ", nrow(resistance), " sites")

# ---------------------------------------------------------------------------
# 5. Recovery  =  mean WUE_T (12 months post-drought) / mean WUE_T (12 months pre-drought)
#    Uses SPEI_3 < -1, minimum event duration 2 months.
# ---------------------------------------------------------------------------

events <- identify_drought_events(
  base_df[is.finite(base_df[[single_spei_col]]), ],
  spei_col = single_spei_col,
  threshold = -1,
  min_duration = 2
)
message("Drought events identified: ", nrow(events))

window_months <- 12L   # pre- and post-drought window length
min_window_n  <- 3L    # minimum months with valid WUE_T needed in each window

recovery_list <- lapply(seq_len(nrow(events)), function(i) {
  ev       <- events[i, ]
  site_df  <- base_df[base_df$site_name == ev$site_name, ]
  site_df  <- site_df[order(site_df$Year, site_df$month), ]
  n        <- nrow(site_df)

  pre_end   <- ev$start_idx - 1L
  pre_start <- max(1L, pre_end - window_months + 1L)
  post_start <- ev$end_idx + 1L
  post_end   <- min(n, post_start + window_months - 1L)

  # check windows exist
  if (pre_end < pre_start || post_start > n) return(NULL)

  pre_wue  <- site_df$WUE_T[pre_start:pre_end]
  post_wue <- site_df$WUE_T[post_start:post_end]
  pre_wue  <- pre_wue[is.finite(pre_wue)]
  post_wue <- post_wue[is.finite(post_wue)]

  if (length(pre_wue)  < min_window_n) return(NULL)
  if (length(post_wue) < min_window_n) return(NULL)

  mean_pre  <- mean(pre_wue)
  mean_post <- mean(post_wue)
  if (mean_pre <= 0) return(NULL)

  data.frame(
    site_name         = ev$site_name,
    water_class       = site_df$water_class[1],
    event_id          = ev$event_id,
    start_year        = ev$start_year,
    start_month       = ev$start_month,
    end_year          = ev$end_year,
    end_month         = ev$end_month,
    duration_months   = ev$duration_months,
    n_pre             = length(pre_wue),
    n_post            = length(post_wue),
    mean_WUE_T_pre    = mean_pre,
    mean_WUE_T_post   = mean_post,
    recovery          = mean_post / mean_pre
  )
})

recovery_events <- do.call(rbind, Filter(Negate(is.null), recovery_list))
if (!is.null(recovery_events)) rownames(recovery_events) <- NULL

# Site-level summary: mean recovery across all events at each site
recovery_site <- recovery_events |>
  dplyr::group_by(site_name, water_class) |>
  dplyr::summarise(
    n_events          = dplyr::n(),
    mean_recovery     = mean(recovery),
    median_recovery   = median(recovery),
    sd_recovery       = sd(recovery),
    min_recovery      = min(recovery),
    max_recovery      = max(recovery),
    pct_full_recovery = mean(recovery >= 1.0) * 100,  # % events where post >= pre
    .groups = "drop"
  )

write_table(recovery_events, "Q2_site_recovery_events.csv")
write_table(recovery_site,   "Q2_site_recovery.csv")
message("Recovery: ", nrow(recovery_events), " events across ",
        nrow(recovery_site), " sites")

# ---------------------------------------------------------------------------
# Summary table: one row per site, all metrics combined
# ---------------------------------------------------------------------------

# Use SPEI_3 slope as the representative plasticity slope when one SPEI is needed.
plasticity_slope_3 <- plasticity_slope |>
  dplyr::filter(SPEI_timescale == single_spei_col) |>
  dplyr::select(site_name, plasticity_slope_SPEI3 = slope,
                plasticity_r2_SPEI3 = r_squared)

summary_table <- stability |>
  dplyr::select(site_name, water_class, n_months, mean_TET, mean_WUE_T, var_WUE_T, cv_WUE_T, stability) |>
  dplyr::left_join(
    plasticity_slope_3,
    by = "site_name"
  ) |>
  dplyr::left_join(
    plasticity_range |> dplyr::select(site_name, plasticity_range, plasticity_p95_p05,
                                       WUE_T_max, WUE_T_min),
    by = "site_name"
  ) |>
  dplyr::left_join(
    resistance |> dplyr::select(site_name, resistance,
                                 mean_WUE_T_drought, mean_WUE_T_near_normal),
    by = "site_name"
  ) |>
  dplyr::left_join(
    recovery_site |> dplyr::select(site_name, n_events, mean_recovery,
                                    pct_full_recovery),
    by = "site_name"
  )

write_table(summary_table, "Q2_site_performance_summary.csv")
message("Summary table: ", nrow(summary_table), " sites")

# ---------------------------------------------------------------------------
# Ecosystem-class tests for site-level performance metrics
# ---------------------------------------------------------------------------

ecosystem_tests <- dplyr::bind_rows(
  ecosystem_metric_test(summary_table, "stability", "Stability index"),
  ecosystem_metric_test(summary_table, "plasticity_slope_SPEI3", "Plasticity slope (SPEI-3)"),
  ecosystem_metric_test(summary_table, "plasticity_p95_p05", "Plasticity range (95th/5th)"),
  ecosystem_metric_test(summary_table, "resistance", "Drought resistance"),
  ecosystem_metric_test(summary_table, "mean_recovery", "Drought recovery")
)

write_table(ecosystem_tests, "Q2_ecosystem_class_statistical_tests.csv")
message("Ecosystem-class tests: ", nrow(ecosystem_tests), " metrics")

# ---------------------------------------------------------------------------
# Colour palette (consistent with rest of workflow)
# ---------------------------------------------------------------------------

ecosystem_colors <- c(
  "Upland"     = "#800080",
  "Freshwater" = "#0000FF",
  "Saline"     = "#FFA500"
)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# -- 1. Stability across sites (dot + range) ---------------------------------
p_stability <- ggplot2::ggplot(
  stability,
  ggplot2::aes(x = mean_TET, y = stability, color = water_class)
) +
  ggplot2::geom_point(alpha = 0.75, size = 2.2) +
  ggplot2::scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1)) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title   = expression(WUE[T]~"Stability by Site"),
    subtitle = expression("Stability index = mean / variance; higher values indicate less variable "~WUE[T]~" relative to its mean"),
    x       = "Mean T:ET ratio",
    y       = "Stability index: mean / variance (log scale)",
    color   = "Ecosystem class"
  ) +
  ggplot2::scale_color_manual(values = ecosystem_colors) +
  theme_wue()

save_plot(p_stability, "Q2_plot_01_stability.png")

# -- 2. Plasticity slope across SPEI timescales --------------------------------
p_plasticity_slope <- ggplot2::ggplot(
  plasticity_slope,
  ggplot2::aes(x = SPEI_timescale, y = slope, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, color = "gray35") +
  ggplot2::geom_boxplot(outlier.alpha = 0.15, linewidth = 0.5) +
  ggplot2::labs(
    title    = expression(WUE[T]~"Plasticity: Slope of "~WUE[T]~" ~ SPEI"),
    subtitle = expression("Site-level OLS slopes; negative slope = "~WUE[T]~" rises with drought"),
    x        = "SPEI timescale",
    y        = "Slope (g C mm⁻¹ per SPEI unit)",
    fill     = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

save_plot(p_plasticity_slope, "Q2_plot_02_plasticity_slope.png")

# -- 3. Plasticity range (WUE_max / WUE_min) -----------------------------------
p_plasticity_range <- ggplot2::ggplot(
  plasticity_range,
  ggplot2::aes(x = water_class, y = plasticity_p95_p05, fill = water_class)
) +
  ggplot2::geom_violin(alpha = 0.5, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.2, linewidth = 0.5) +
  ggplot2::labs(
    title    = expression(WUE[T]~"Plasticity: Dynamic Range"),
    subtitle = expression("Ratio of 95th to 5th percentile "~WUE[T]~" per site (robust version of max/min)"),
    x        = "Ecosystem class",
    y        = expression(WUE[T]~"95th / 5th percentile"),
    fill     = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(legend.position = "none")

save_plot(p_plasticity_range, "Q2_plot_03_plasticity_range.png", width = 8, height = 7)

# -- 4. Resistance -------------------------------------------------------------
p_resistance <- ggplot2::ggplot(
  resistance,
  ggplot2::aes(x = water_class, y = resistance, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.5, linetype = "dashed", color = "gray35") +
  ggplot2::geom_violin(alpha = 0.5, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.2, linewidth = 0.5) +
  ggplot2::labs(
    title    = expression(WUE[T]~"Resistance to Drought"),
    subtitle = paste0("Ratio of mean transpiration-based WUE during drought (", single_spei_label, " < -1) to near-normal conditions\nValues > 1 indicate an increase under drought"),
    x        = "Ecosystem class",
    y        = expression("Resistance (drought "~WUE[T]~" / near-normal "~WUE[T]~")"),
    fill     = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(legend.position = "none")

save_plot(p_resistance, "Q2_plot_04_resistance.png", width = 8, height = 7)

# -- 5. Recovery ---------------------------------------------------------------
p_recovery <- ggplot2::ggplot(
  recovery_site,
  ggplot2::aes(x = water_class, y = mean_recovery, fill = water_class)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.5, linetype = "dashed", color = "gray35") +
  ggplot2::geom_violin(alpha = 0.5, color = NA) +
  ggplot2::geom_boxplot(width = 0.18, outlier.alpha = 0.2, linewidth = 0.5) +
  ggplot2::labs(
    title    = expression(WUE[T]~"Recovery Post-Drought"),
    subtitle = paste0("Mean post- to pre-drought transpiration-based WUE ratio per site across all ", single_spei_label, " drought events\nValues > 1 indicate full or overshoot recovery"),
    x        = "Ecosystem class",
    y        = expression("Recovery (post-drought "~WUE[T]~" / pre-drought "~WUE[T]~")"),
    fill     = "Ecosystem class"
  ) +
  ggplot2::scale_fill_manual(values = ecosystem_colors) +
  theme_wue() +
  ggplot2::theme(legend.position = "none")

save_plot(p_recovery, "Q2_plot_05_recovery.png", width = 8, height = 7)

# -- 6. Composite: all five metrics side by side (summary table heatmap) ------
# Rank sites by stability and show all five normalised scores

metrics_long <- summary_table |>
  dplyr::filter(
    is.finite(stability),
    is.finite(plasticity_slope_SPEI3),
    is.finite(plasticity_range),
    is.finite(resistance),
    is.finite(mean_recovery)
  ) |>
  dplyr::mutate(
    stability_z         = as.numeric(scale(stability)),
    plasticity_slope_z  = as.numeric(scale(plasticity_slope_SPEI3)),
    plasticity_range_z  = as.numeric(scale(plasticity_range)),
    resistance_z        = as.numeric(scale(resistance)),
    recovery_z          = as.numeric(scale(mean_recovery))
  ) |>
  dplyr::select(site_name, water_class,
                stability_z, plasticity_slope_z, plasticity_range_z,
                resistance_z, recovery_z) |>
  tidyr::pivot_longer(
    cols      = ends_with("_z"),
    names_to  = "metric",
    values_to = "z_score"
  )

metrics_long$metric <- factor(
  metrics_long$metric,
  levels = c("stability_z", "plasticity_slope_z", "plasticity_range_z",
             "resistance_z", "recovery_z"),
  labels = c("Stability index\n(mean/variance)",
             "Plasticity\n(slope~SPEI-3)",
             "Plasticity\n(max/min)",
             "Resistance\n(drought/normal)",
             "Recovery\n(post/pre)")
)

site_order <- metrics_long |>
  dplyr::filter(metric == "Stability index\n(mean/variance)") |>
  dplyr::arrange(z_score) |>
  dplyr::pull(site_name)


metrics_long$site_name <- factor(metrics_long$site_name, levels = site_order)

p_heatmap <- ggplot2::ggplot(
  metrics_long,
  ggplot2::aes(x = metric, y = site_name, fill = z_score)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::facet_grid(water_class ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_gradient2(
    low      = "#2C7BB6",
    mid      = "white",
    high     = "#D7191C",
    midpoint = 0,
    name     = "Z-score"
  ) +
  ggplot2::labs(
    title    = expression(WUE[T]~"Performance Profile: All Sites"),
    subtitle = "Z-scores within each metric; red = above average, blue = below average",
    x        = NULL,
    y        = NULL
  ) +
  theme_wue(base_size = 10) +
  ggplot2::theme(
    axis.text.y       = ggplot2::element_text(size = 6),
    strip.text.y      = ggplot2::element_text(angle = 0),
    legend.position   = "right",
    panel.grid        = ggplot2::element_blank()
  )

save_plot(p_heatmap, "Q2_plot_06_performance_heatmap.png", width = 10, height = 14)

# ---------------------------------------------------------------------------
# Composite panel (plots 1–5)
# ---------------------------------------------------------------------------

theme_q2_composite <- function() {
  theme_wue(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      plot.subtitle = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 15),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13)
    )
}

p_stability_compact <- p_stability +
  ggplot2::labs(tag = "A", title = "Stability", subtitle = NULL, y = "Stability index", color = "Ecosystem") +
  theme_q2_composite() +
  ggplot2::theme(legend.position = "bottom")

p_plasticity_slope_compact <- p_plasticity_slope +
  ggplot2::labs(tag = "B", title = "Plasticity slope", subtitle = NULL, y = "Slope", fill = "Ecosystem") +
  theme_q2_composite() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 12)) +
  ggplot2::guides(fill = "none")

p_plasticity_range_compact <- p_plasticity_range +
  ggplot2::labs(tag = "C", title = "Plasticity range", subtitle = NULL, y = expression(WUE[T]~"95th / 5th"), fill = "Ecosystem") +
  theme_q2_composite() +
  ggplot2::guides(fill = "none")

p_resistance_compact <- p_resistance +
  ggplot2::labs(tag = "D", title = "Drought Resistance", subtitle = NULL, y = "Resistance", fill = "Ecosystem") +
  theme_q2_composite() +
  ggplot2::guides(fill = "none")

p_recovery_compact <- p_recovery +
  ggplot2::labs(tag = "E", title = "Drought Recovery", subtitle = NULL, y = "Recovery", fill = "Ecosystem") +
  theme_q2_composite() +
  ggplot2::guides(fill = "none")

composite <- (p_stability_compact | p_plasticity_slope_compact) /
             (p_plasticity_range_compact | p_resistance_compact | p_recovery_compact) +
  patchwork::plot_layout(guides = "collect", heights = c(1.05, 1.0)) +
  patchwork::plot_annotation(
    title    = expression("Q2 "~WUE[T]~" Performance"),
    subtitle = paste0("Site-level metrics using transpiration-based WUE and ", single_spei_label, " drought classification")
  ) &
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold", size = 19),
    plot.subtitle = ggplot2::element_text(size = 14),
    plot.tag = ggplot2::element_text(face = "bold", size = 18)
  )

save_plot(composite, "Q2_plot_07_composite.png", width = 12, height = 10)

# ---------------------------------------------------------------------------
# Summary text log
# ---------------------------------------------------------------------------

capture_to_file({
  cat("Q2_WUE_performance.R — run summary\n")
  cat("Input:  ", input_file, "\n")
  cat("Output: ", output_dir, "\n\n")
  cat("Single-SPEI drought metric:", single_spei_label, "\n\n")

  cat("--- Data ---\n")
  cat("Rows in cleaned monthly data:", nrow(base_df), "\n")
  cat("Sites:", length(unique(base_df$site_name)), "\n")
  cat("Years:", min(base_df$Year), "-", max(base_df$Year), "\n\n")

  cat("--- 1. Stability ---\n")
  cat("Sites with stability estimates:", nrow(stability), "\n")
  print(summary(stability[, c("mean_WUE_T", "var_WUE_T", "cv_WUE_T", "stability")]))
  cat("\nBy ecosystem class:\n")
  print(
    stability |>
      dplyr::group_by(water_class) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_stability = mean(stability),
        median_stability = median(stability),
        .groups = "drop"
      )
  )

  cat("\n--- 2. Plasticity: slope ---\n")
  cat("Site × timescale rows:", nrow(plasticity_slope), "\n")
  print(
    plasticity_slope |>
      dplyr::group_by(SPEI_timescale, water_class) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_slope   = mean(slope),
        pct_negative = mean(slope < 0) * 100,
        .groups = "drop"
      )
  )

  cat("\n--- 3. Plasticity: range ---\n")
  cat("Sites with range estimates:", nrow(plasticity_range), "\n")
  print(summary(plasticity_range[, c("WUE_T_max", "WUE_T_min",
                                      "plasticity_range", "plasticity_p95_p05")]))

  cat("\n--- 4. Resistance ---\n")
  cat("Sites with resistance estimates:", nrow(resistance), "\n")
  print(summary(resistance[, c("mean_WUE_T_drought", "mean_WUE_T_near_normal", "resistance")]))
  cat("\nBy ecosystem class:\n")
  print(
    resistance |>
      dplyr::group_by(water_class) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_resistance   = mean(resistance),
        pct_above_1       = mean(resistance > 1) * 100,
        .groups = "drop"
      )
  )

  cat("\n--- 5. Recovery ---\n")
  cat("Drought events:", nrow(recovery_events), "\n")
  cat("Sites with recovery estimates:", nrow(recovery_site), "\n")
  print(summary(recovery_events[, c("duration_months", "mean_WUE_T_pre",
                                     "mean_WUE_T_post", "recovery")]))
  cat("\nBy ecosystem class:\n")
  print(
    recovery_site |>
      dplyr::group_by(water_class) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_recovery       = mean(mean_recovery),
        pct_full_recovery   = mean(pct_full_recovery),
        .groups = "drop"
      )
  )

  cat("\n--- Summary table coverage ---\n")
  cat("Sites in summary table:", nrow(summary_table), "\n")
  cat("Sites with all five metrics:",
      sum(complete.cases(summary_table[, c("stability", "plasticity_slope_SPEI3",
                                           "plasticity_range", "resistance",
                                           "mean_recovery")])), "\n")

  cat("\n--- Ecosystem-class statistical tests ---\n")
  print(ecosystem_tests)
}, "Q2_WUE_performance_summary.txt")

message("Done. All outputs written to: ", output_dir)
