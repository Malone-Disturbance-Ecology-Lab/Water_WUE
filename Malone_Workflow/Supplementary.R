# Supplementary.R
#
# Supplemental figures for the Malone WUE workflow.
#
# Current output:
#   results/EDI_spatial/Supplementary_plot_01_mean_SPEI3_timeseries_by_coast.png
#   results/EDI_spatial/Supplementary_plot_02_cumulative_SPEI3_by_coast.png
#
# This figure summarizes the raster-based monthly mean SPEI-3 signal for each
# coastal region used in the spatial EDI workflow.

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "ggplot2", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required. Install with install.packages('%s').", pkg, pkg))
    }
  }
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  candidates <- c(
    "Supplementary.R",
    file.path("Malone_Workflow", "Supplementary.R"),
    file.path("Water_WUE", "Malone_Workflow", "Supplementary.R")
  )
  matches <- candidates[file.exists(candidates)]
  if (length(matches) > 0) {
    return(normalizePath(matches[1], mustWork = TRUE))
  }
  stop("Could not locate Supplementary.R. Run from Water_WUE or Water_WUE/Malone_Workflow, or call with Rscript path/to/Supplementary.R.")
}

save_plot <- function(plot, filename, width = 7.2, height = 7.5) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

theme_supplement <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      legend.position = "bottom"
    )
}

script_dir <- dirname(get_script_path())
output_dir <- file.path(script_dir, "results", "EDI_spatial")
input_file <- file.path(output_dir, "EDI_response_v2_upland_monthly_coast_summary.csv")

if (!file.exists(input_file)) {
  fallback_file <- file.path(output_dir, "EDI_logistic_monthly_coast_summary.csv")
  if (!file.exists(fallback_file)) {
    stop(
      "Could not find monthly coast summary. Expected one of:\n",
      input_file, "\n",
      fallback_file, "\nRun function/24-EDI_spatial.py first."
    )
  }
  input_file <- fallback_file
}

monthly <- read.csv(input_file, stringsAsFactors = FALSE)
required_cols <- c("year", "month", "coast_region", "mean_SPEI3")
missing_cols <- setdiff(required_cols, names(monthly))
if (length(missing_cols) > 0) {
  stop("Missing required columns in monthly summary: ", paste(missing_cols, collapse = ", "))
}

coast_levels <- c("AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast")
coast_labels <- c(
  "AK Coast" = "Alaska Coast",
  "Pacific Coast" = "Pacific Coast",
  "Gulf Coast" = "Gulf Coast",
  "Atlantic Coast" = "Atlantic Coast"
)
coast_colors <- c(
  "AK Coast" = "#4D4D4D",
  "Pacific Coast" = "#C44E52",
  "Gulf Coast" = "#8C564B",
  "Atlantic Coast" = "#008B8B"
)

monthly$date <- as.Date(sprintf("%04d-%02d-15", monthly$year, monthly$month))
monthly$coast_region <- factor(monthly$coast_region, levels = coast_levels)
monthly <- monthly |>
  dplyr::filter(!is.na(coast_region), is.finite(mean_SPEI3)) |>
  dplyr::arrange(coast_region, date)
date_min <- min(monthly$date)
date_max <- max(monthly$date)

annual <- monthly |>
  dplyr::group_by(coast_region, year) |>
  dplyr::summarise(
    annual_mean_SPEI3 = mean(mean_SPEI3, na.rm = TRUE),
    .groups = "drop"
  )

cumulative <- monthly |>
  dplyr::group_by(coast_region) |>
  dplyr::arrange(date, .by_group = TRUE) |>
  dplyr::mutate(cumulative_SPEI3 = cumsum(mean_SPEI3)) |>
  dplyr::ungroup()

annual_cumulative <- annual |>
  dplyr::group_by(coast_region) |>
  dplyr::arrange(year, .by_group = TRUE) |>
  dplyr::mutate(cumulative_annual_SPEI3 = cumsum(annual_mean_SPEI3)) |>
  dplyr::ungroup()

write.csv(
  annual,
  file.path(output_dir, "Supplementary_table_01_annual_mean_SPEI3_by_coast.csv"),
  row.names = FALSE
)

write.csv(
  cumulative,
  file.path(output_dir, "Supplementary_table_02_monthly_cumulative_SPEI3_by_coast.csv"),
  row.names = FALSE
)

write.csv(
  annual_cumulative,
  file.path(output_dir, "Supplementary_table_03_annual_cumulative_SPEI3_by_coast.csv"),
  row.names = FALSE
)

spei_plot <- ggplot2::ggplot(monthly, ggplot2::aes(x = date, y = mean_SPEI3, color = coast_region)) +
  ggplot2::annotate("rect", xmin = date_min, xmax = date_max, ymin = -1, ymax = 1,
                    fill = "gray85", alpha = 0.30) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.45) +
  ggplot2::geom_hline(yintercept = c(-1, 1), color = "gray50", linetype = "dashed", linewidth = 0.40) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.85) +
  ggplot2::geom_line(
    data = annual,
    ggplot2::aes(
      x = as.Date(sprintf("%04d-07-01", year)),
      y = annual_mean_SPEI3,
      color = coast_region
    ),
    linewidth = 1.25,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~ coast_region, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(coast_labels)) +
  ggplot2::scale_color_manual(values = coast_colors, guide = "none") +
  ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  ggplot2::labs(
    title = "Monthly Mean SPEI-3 by Coastal Region",
    subtitle = "Thin lines show monthly raster means; thick lines show annual means. Gray band marks near-normal SPEI-3 (-1 to 1).",
    x = NULL,
    y = "Mean SPEI-3"
  ) +
  theme_supplement(base_size = 12) +
  ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    panel.spacing.y = grid::unit(0.8, "lines")
  )

save_plot(
  spei_plot,
  file.path(output_dir, "Supplementary_plot_01_mean_SPEI3_timeseries_by_coast.png"),
  width = 7.2,
  height = 8.2
)

cumulative_plot <- ggplot2::ggplot(
  cumulative,
  ggplot2::aes(x = date, y = cumulative_SPEI3, color = coast_region)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray35", linewidth = 0.55) +
  ggplot2::geom_line(linewidth = 0.70, alpha = 0.88) +
  ggplot2::geom_line(
    data = annual_cumulative,
    ggplot2::aes(
      x = as.Date(sprintf("%04d-07-01", year)),
      y = cumulative_annual_SPEI3,
      color = coast_region
    ),
    linewidth = 1.45,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~ coast_region, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(coast_labels)) +
  ggplot2::scale_color_manual(values = coast_colors, guide = "none") +
  ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
  ggplot2::labs(
    title = "Cumulative Mean SPEI-3 by Coastal Region",
    subtitle = "Thin lines show cumulative monthly raster means; thick lines show cumulative annual means.",
    x = NULL,
    y = "Cumulative mean SPEI-3"
  ) +
  theme_supplement(base_size = 12) +
  ggplot2::theme(
    plot.subtitle = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    panel.spacing.y = grid::unit(0.8, "lines")
  )

save_plot(
  cumulative_plot,
  file.path(output_dir, "Supplementary_plot_02_cumulative_SPEI3_by_coast.png"),
  width = 7.2,
  height = 8.2
)

message("Supplementary outputs written to: ", output_dir)
