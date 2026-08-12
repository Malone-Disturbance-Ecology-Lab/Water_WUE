# =============================================================================
# READ-ONLY DIAGNOSTIC:
# SPEI-3 vs SPEI-48 WUE_ET / WUE_T anomaly comparison
#
# Same logic as the original Q2 anomaly workflow, except:
# original = SPEI-6 vs SPEI-48
# here     = SPEI-3 vs SPEI-48
#
# NOTHING IS SAVED.
# =============================================================================

rm(list = ls())
options(width = 180)

suppressPackageStartupMessages({
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required.")
  }
})

library(lme4)
library(dplyr)
library(emmeans)

# =============================================================================
# INPUT
# =============================================================================

input_file <- paste0(
  "//corellia.environment.yale.edu/MaloneLab/Research/WUE_CUE/",
  "data_products/WUE_CUE_monthly_merged_indices_clean.csv"
)

ecosystem_classes <- c("Upland", "Freshwater", "Saline")

# =============================================================================
# HELPERS
# =============================================================================

classify_spei <- function(x) {
  out <- ifelse(
    x < -1,
    "Dry",
    ifelse(x > 1, "Wet", "Near-normal")
  )
  
  factor(
    out,
    levels = c("Dry", "Near-normal", "Wet")
  )
}

coef_table <- function(model) {
  out <- as.data.frame(summary(model)$coefficients)
  out$term <- rownames(out)
  rownames(out) <- NULL
  
  out$p_approx <- 2 * pnorm(-abs(out[, "t value"]))
  
  out[, c(
    "term",
    "Estimate",
    "Std. Error",
    "t value",
    "p_approx"
  )]
}

# =============================================================================
# LOAD + CLEAN — SAME GENERAL LOGIC AS ORIGINAL ANOMALY WORKFLOW
# =============================================================================

monthly <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

monthly[] <- lapply(
  monthly,
  function(x) {
    if (is.character(x)) trimws(x) else x
  }
)

required_cols <- c(
  "site_name",
  "Year",
  "month",
  "water_class",
  "WUE",
  "WUE_tra",
  "Trans_ratio",
  "SPEI_3",
  "SPEI_48"
)

missing_cols <- setdiff(
  required_cols,
  names(monthly)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

model_base <- monthly[
  complete.cases(
    monthly[, required_cols]
  ),
]

model_base <- model_base[
  model_base$site_name != "" &
    model_base$water_class != "" &
    model_base$water_class %in% ecosystem_classes &
    is.finite(model_base$WUE) &
    is.finite(model_base$WUE_tra) &
    is.finite(model_base$Trans_ratio) &
    model_base$Trans_ratio >= 0 &
    model_base$Trans_ratio <= 1 &
    is.finite(model_base$SPEI_3) &
    is.finite(model_base$SPEI_48),
]

model_base$site_name <- factor(
  model_base$site_name
)

model_base$water_class <- factor(
  model_base$water_class,
  levels = ecosystem_classes
)

model_base$month_f <- factor(
  model_base$month,
  levels = 1:12
)

model_base$Trans_ratio_z <- as.numeric(
  scale(model_base$Trans_ratio)
)

model_base$WUE_ET <- model_base$WUE
model_base$WUE_T <- model_base$WUE_tra

cat("\n============================================================\n")
cat("DATA CHECK\n")
cat("============================================================\n")

cat("Rows :", nrow(model_base), "\n")
cat(
  "Sites:",
  length(unique(model_base$site_name)),
  "\n"
)

# =============================================================================
# WUE METRICS TO LONG FORMAT
# =============================================================================

metric_long <- reshape(
  model_base,
  varying = c("WUE_ET", "WUE_T"),
  v.names = "WUE_value",
  timevar = "metric",
  times = c("WUE_ET", "WUE_T"),
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)

metric_long$metric <- factor(
  metric_long$metric,
  levels = c("WUE_ET", "WUE_T")
)

# =============================================================================
# SPEI-3 AND SPEI-48 TO LONG FORMAT
# =============================================================================

spei_long <- rbind(
  
  transform(
    metric_long,
    timescale = "SPEI_3",
    SPEI_value = SPEI_3
  ),
  
  transform(
    metric_long,
    timescale = "SPEI_48",
    SPEI_value = SPEI_48
  )
)

# IMPORTANT:
# SPEI-3 is the reference timescale.
spei_long$timescale <- factor(
  spei_long$timescale,
  levels = c("SPEI_3", "SPEI_48")
)

spei_long$anomaly_class <- classify_spei(
  spei_long$SPEI_value
)

spei_long$water_class <- factor(
  spei_long$water_class,
  levels = ecosystem_classes
)

spei_long$site_name <- factor(
  spei_long$site_name
)

spei_long$metric <- factor(
  spei_long$metric,
  levels = c("WUE_ET", "WUE_T")
)

spei_long$month_f <- factor(
  spei_long$month,
  levels = 1:12
)

# =============================================================================
# MONTH-SPECIFIC NEAR-NORMAL BASELINES
#
# Exactly the important logic from the original anomaly workflow:
# site × timescale × metric × calendar month
# Require at least 3 near-normal values.
# =============================================================================

baseline <- spei_long |>
  filter(
    anomaly_class == "Near-normal"
  ) |>
  group_by(
    site_name,
    timescale,
    metric,
    month_f
  ) |>
  summarise(
    baseline_mean = mean(WUE_value),
    baseline_median = median(WUE_value),
    baseline_n = n(),
    .groups = "drop"
  )

analysis_data <- left_join(
  spei_long,
  baseline,
  by = c(
    "site_name",
    "timescale",
    "metric",
    "month_f"
  )
)

analysis_data <- analysis_data[
  is.finite(analysis_data$baseline_mean) &
    analysis_data$baseline_n >= 3,
]

analysis_data$WUE_response <-
  analysis_data$WUE_value -
  analysis_data$baseline_mean

analysis_data$WUE_response_magnitude <-
  abs(analysis_data$WUE_response)

analysis_data$response_direction <- factor(
  ifelse(
    analysis_data$WUE_response < 0,
    "Decrease",
    "Increase_or_no_decrease"
  ),
  levels = c(
    "Increase_or_no_decrease",
    "Decrease"
  )
)

analysis_data$decrease_binary <- as.integer(
  analysis_data$response_direction == "Decrease"
)

anomaly_data <- analysis_data[
  analysis_data$anomaly_class != "Near-normal",
]

anomaly_data$anomaly_class <- factor(
  as.character(
    anomaly_data$anomaly_class
  ),
  levels = c("Dry", "Wet")
)

cat("\nAnalysis rows:", nrow(analysis_data), "\n")
cat(
  "Analysis sites:",
  length(unique(analysis_data$site_name)),
  "\n"
)

# =============================================================================
# 1. DESCRIPTIVE DRY RESPONSES
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("1. DRY RESPONSE: SPEI-3 VS SPEI-48\n")
cat("============================================================\n")

dry_summary <- analysis_data |>
  filter(
    anomaly_class == "Dry"
  ) |>
  group_by(
    timescale,
    metric
  ) |>
  summarise(
    n_months = n(),
    n_sites = n_distinct(site_name),
    mean_WUE = mean(WUE_value),
    mean_response = mean(WUE_response),
    median_response = median(WUE_response),
    mean_abs_response =
      mean(WUE_response_magnitude),
    pct_decrease =
      mean(decrease_binary) * 100,
    .groups = "drop"
  )

print(dry_summary)

# =============================================================================
# 2. FIT EXACT SAME TYPE OF SIGNED-RESPONSE MODEL
# =============================================================================

response_formula <-
  WUE_response ~
  metric *
  timescale *
  anomaly_class *
  water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

response_model <- lmer(
  response_formula,
  data = analysis_data,
  REML = FALSE
)

# =============================================================================
# 3. RESPONSE-MAGNITUDE MODEL
# =============================================================================

magnitude_formula <-
  WUE_response_magnitude ~
  metric *
  timescale *
  anomaly_class *
  water_class +
  Trans_ratio_z +
  month_f +
  (1 | site_name)

magnitude_model <- lmer(
  magnitude_formula,
  data = anomaly_data,
  REML = FALSE
)

# =============================================================================
# 4. PRINT RELEVANT MODEL COEFFICIENTS
#
# Because references are:
# metric        = WUE_ET
# timescale     = SPEI_3
# anomaly class = Dry
# ecosystem     = Upland
#
# metricWUE_T:timescaleSPEI_48 therefore directly asks:
# does WUE_T-vs-WUE_ET divergence differ at SPEI-48 relative to SPEI-3
# for the reference dry/upland condition?
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("2. SIGNED-RESPONSE MODEL: KEY SPEI-48 COEFFICIENTS\n")
cat("============================================================\n")

response_coefs <- coef_table(
  response_model
)

keep_signed <- grepl(
  "metricWUE_T:timescaleSPEI_48",
  response_coefs$term,
  fixed = TRUE
)

print(
  response_coefs[
    keep_signed,
  ],
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("3. RESPONSE-MAGNITUDE MODEL: KEY SPEI-48 COEFFICIENTS\n")
cat("============================================================\n")

magnitude_coefs <- coef_table(
  magnitude_model
)

keep_mag <- grepl(
  "metricWUE_T:timescaleSPEI_48",
  magnitude_coefs$term,
  fixed = TRUE
)

print(
  magnitude_coefs[
    keep_mag,
  ],
  row.names = FALSE
)

# =============================================================================
# 5. DIRECT MARGINAL COMPARISON:
#
# First estimate WUE_T - WUE_ET at Dry SPEI-3 and Dry SPEI-48,
# averaging over ecosystem classes.
#
# Then test whether those two metric divergences differ.
#
# This is the most useful test for your Discussion sentence.
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("4. MARGINAL WUE_T - WUE_ET DIVERGENCE DURING DRY CONDITIONS\n")
cat("============================================================\n")

emm_metric <- emmeans(
  response_model,
  ~ metric | timescale * anomaly_class
)

metric_difference <- contrast(
  emm_metric,
  method = "revpairwise",
  adjust = "none"
)

metric_difference_df <- as.data.frame(
  metric_difference
)

dry_metric_difference <- metric_difference_df[
  metric_difference_df$anomaly_class == "Dry",
]

print(
  dry_metric_difference,
  row.names = FALSE
)

# =============================================================================
# 6. DIRECT TEST:
# Is dry WUE_T-WUE_ET divergence different at SPEI-48 vs SPEI-3?
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("5. DIRECT TEST: SPEI-48 VS SPEI-3 METRIC DIVERGENCE\n")
cat("============================================================\n")

# Create EMMs for metric × timescale specifically under Dry conditions
emm_dry <- emmeans(
  response_model,
  ~ metric * timescale,
  at = list(
    anomaly_class = "Dry"
  )
)

# Coefficients correspond to:
# 1 = WUE_ET SPEI_3
# 2 = WUE_T  SPEI_3
# 3 = WUE_ET SPEI_48
# 4 = WUE_T  SPEI_48
#
# We test:
# (WUE_T - WUE_ET at SPEI48)
# minus
# (WUE_T - WUE_ET at SPEI3)
#
# = [-1, +1, +1, -1] depending on ordering check below.

cat("\nEMMEANS ordering:\n")
print(
  as.data.frame(emm_dry)[
    ,
    c(
      "metric",
      "timescale",
      "emmean",
      "SE"
    )
  ],
  row.names = FALSE
)

# Use named interaction contrast safely rather than assuming interpretation.
contrast_test <- contrast(
  emm_dry,
  interaction = c(
    "revpairwise",
    "revpairwise"
  ),
  adjust = "none"
)

cat("\nDifference-in-differences:\n")
print(
  as.data.frame(contrast_test),
  row.names = FALSE
)

# =============================================================================
# 7. OPTIONAL: SAME TEST WITHIN EACH ECOSYSTEM
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("6. DRY METRIC DIVERGENCE BY ECOSYSTEM\n")
cat("============================================================\n")

emm_by_ecosystem <- emmeans(
  response_model,
  ~ metric * timescale | water_class,
  at = list(
    anomaly_class = "Dry"
  )
)

ecosystem_interactions <- contrast(
  emm_by_ecosystem,
  interaction = c(
    "revpairwise",
    "revpairwise"
  ),
  adjust = "none"
)

print(
  as.data.frame(ecosystem_interactions),
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("DIAGNOSTIC COMPLETE — NOTHING WAS SAVED\n")
cat("============================================================\n")