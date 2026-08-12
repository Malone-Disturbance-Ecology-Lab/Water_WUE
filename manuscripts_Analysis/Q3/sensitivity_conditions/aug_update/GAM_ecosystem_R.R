# Q3_generate_ecosystem_GAM_predictions.R
# Generate ecosystem GAM predictions across all months and full SPEI range.
# Does NOT refit the model – uses the saved GAM.

library(mgcv)
library(dplyr)

# ---- Paths ----
base_dir <- "M:/Research/WUE_CUE/WUE_manuscript_version6/Q3"

model_rds <- file.path(
  base_dir,
  "Q3_WUE_T_SPEI_sensitivity_outputs",
  "Q3_WUE_T_SPEI_smooth_gam_model.rds"
)

data_csv <- file.path(
  base_dir,
  "Q3_WUE_T_SPEI_sensitivity_outputs",
  "Q3_WUE_T_SPEI_model_data_long.csv"
)

output_csv <- file.path(
  base_dir,
  "Q3_august_update",
  "Q3_ecosystem_GAM_predictions_all_ecosystems_all_months.csv"
)

# ---- Load model and data ----
smooth_model <- readRDS(model_rds)
data <- read.csv(data_csv)

# ---- Recreate factors exactly as in original Q3 workflow ----
data$site_name <- factor(data$site_name)

data$water_class <- factor(
  data$water_class,
  levels = c("Upland", "Freshwater", "Saline")
)

data$SPEI_timescale <- factor(
  data$SPEI_timescale,
  levels = c("SPEI_1", "SPEI_3", "SPEI_6",
             "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48")
)

data$month_f <- factor(data$month_f)

data$spei_ecosystem <- interaction(
  data$SPEI_timescale,
  data$water_class,
  sep = "__",
  drop = TRUE
)

# ---- Prediction grid ----
months <- levels(data$month_f)   # all months present

spei_lower <- min(data$SPEI_value, na.rm = TRUE)
spei_upper <- max(data$SPEI_value, na.rm = TRUE)
spei_sequence <- seq(spei_lower, spei_upper, length.out = 200)

grid <- expand.grid(
  SPEI_timescale = levels(data$SPEI_timescale),
  water_class    = levels(data$water_class),
  month_f        = months,
  SPEI_value     = spei_sequence,
  site_name      = levels(data$site_name)[1],
  stringsAsFactors = FALSE
)

# ---- Match factors to model-data levels ----
grid$SPEI_timescale <- factor(grid$SPEI_timescale,
                              levels = levels(data$SPEI_timescale))
grid$water_class <- factor(grid$water_class,
                           levels = levels(data$water_class))
grid$month_f <- factor(grid$month_f,
                       levels = levels(data$month_f))
grid$site_name <- factor(grid$site_name,
                         levels = levels(data$site_name))

# Recreate spei_ecosystem in the grid
grid$spei_ecosystem <- interaction(
  grid$SPEI_timescale,
  grid$water_class,
  sep = "__",
  drop = TRUE
)
grid$spei_ecosystem <- factor(grid$spei_ecosystem,
                              levels = levels(data$spei_ecosystem))

grid <- grid[order(grid$SPEI_timescale, grid$water_class,
                   grid$month_f, grid$SPEI_value), ]

# ---- Predict (exclude site random effect) ----
pred <- predict(
  smooth_model,
  newdata = grid,
  type = "link",
  se.fit = TRUE,
  exclude = "s(site_name)"
)

grid$predicted_WUE_T <- pred$fit
grid$predicted_se    <- pred$se.fit
grid$predicted_lower <- pred$fit - 1.96 * pred$se.fit
grid$predicted_upper <- pred$fit + 1.96 * pred$se.fit

# ---- Baseline (SPEI in [-1, 1]) ----
baseline <- grid %>%
  filter(SPEI_value >= -1 & SPEI_value <= 1) %>%
  group_by(water_class, SPEI_timescale, month_f) %>%
  summarise(
    predicted_near_normal_WUE_T = mean(predicted_WUE_T, na.rm = TRUE),
    .groups = "drop"
  )

grid <- grid %>%
  left_join(baseline, by = c("water_class", "SPEI_timescale", "month_f"))

# ---- Percent change ----
grid <- grid %>%
  mutate(
    predicted_pct_change = (predicted_WUE_T - predicted_near_normal_WUE_T) /
      predicted_near_normal_WUE_T * 100,
    predicted_lower_pct = (predicted_lower - predicted_near_normal_WUE_T) /
      predicted_near_normal_WUE_T * 100,
    predicted_upper_pct = (predicted_upper - predicted_near_normal_WUE_T) /
      predicted_near_normal_WUE_T * 100
  )

# ---- Select columns ----
out <- grid %>%
  select(
    water_class,
    SPEI_timescale,
    month_f,
    SPEI_value,
    predicted_WUE_T,
    predicted_lower,
    predicted_upper,
    predicted_near_normal_WUE_T,
    predicted_pct_change,
    predicted_lower_pct,
    predicted_upper_pct
  )

# ---- Save ----
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output_csv, row.names = FALSE)

cat("Ecosystem GAM predictions saved to:\n", output_csv, "\n")