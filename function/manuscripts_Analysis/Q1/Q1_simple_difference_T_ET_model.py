"""
COMPLETE Q1 ANALYSIS - Malone Workflow Replication (Python Version)
All Python outputs saved to separate directory structure for safe comparison

DIRECTORY STRUCTURE:
  Q1_mixed_model_outputs/   - CSV files, RDS files, README, model summary
  Q1_mixed_model_figures/   - PNG figure files  
  Q1_mixed_model_temp/      - Temporary R handoff files (cleaned at start)

WORKFLOW:
1. Python: Clean temp dir, data cleaning, filtering, duplicate checks, saves model_base
2. R: Reads wide data, reshapes to long, runs lme4/mgcv, outputs CSVs + RDS + summary
3. Python: Verifies R outputs, reads R output CSVs, creates figures
4. Python: Saves README and final output summary

OUTPUTS: M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\
"""

import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS - PYTHON OUTPUTS IN SEPARATE DIRECTORY (SAFE)
# ============================================================================

# Input data (read from Malone location)
monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"



# Python outputs - SEPARATE from Malone with clean organization
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1"

# Three separate subdirectories for organization
output_dir = os.path.join(base_output_dir, "Q1_mixed_model_outputs")     # CSV, RDS, README, summary
figure_dir = os.path.join(base_output_dir, "Q1_mixed_model_figures")     # PNG figures
temp_dir = os.path.join(base_output_dir, "Q1_mixed_model_temp")          # Temporary R files

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir]:
    os.makedirs(dir_path, exist_ok=True)

print("="*60)
print("Q1 WUE Metric-Difference Models - Python Replication")
print("="*60)
print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"(Separate from Malone outputs - safe for comparison)")

# ============================================================================
# CLEAN TEMP DIRECTORY AT START (BEFORE CREATING NEW TEMP FILES)
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning temp directory at start")
print("="*60)

# Clean temp directory BEFORE creating new temp files
for f in os.listdir(temp_dir):
    path = os.path.join(temp_dir, f)
    if os.path.isfile(path):
        os.remove(path)
        print(f"  Removed old temp file: {f}")

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def save_plot(fig, filename, width=10, height=7):
    """Save plot to figure directory"""
    fig.set_size_inches(width, height)
    fig.savefig(os.path.join(figure_dir, filename), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def theme_wue_manuscript():
    """Malone-style manuscript theme for Python plots"""
    plt.rcParams.update({
        'font.size': 27,
        'axes.labelsize': 27,
        'axes.titlesize': 33,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 23,
        'ytick.labelsize': 23,
        'legend.fontsize': 24,
        'legend.title_fontsize': 26,
        'figure.titlesize': 42,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
    })
    sns.set_style("whitegrid")
    sns.set_context("talk")

def theme_wue():
    """Malone-style base theme for Python plots"""
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.title_fontsize': 10,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.3,
    })
    sns.set_style("whitegrid")

# ============================================================================
# STEP 1: PYTHON DATA CLEANING AND FILTERING (WIDE FORMAT)
# ============================================================================

print("\n" + "="*60)
print("STEP 1: Python Data Cleaning (WIDE format)")
print("="*60)

print("\nLoading data...")
monthly = pd.read_csv(monthly_data_path)
print(f"  Loaded {len(monthly)} rows")

# Safer string trimming
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].where(
        monthly[col].isna(),
        monthly[col].astype(str).str.strip()
    )
    monthly[col] = monthly[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

# Required columns
required_cols = [
    "site_name", "Year", "month", "water_class", "Trans_ratio",
    "WUE", "WUE_tra", "WUE_eva", "SPEI_1"
]

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

ecosystem_classes = ["Upland", "Freshwater", "Saline"]

# Filter - exact Malone
model_base = monthly.dropna(subset=required_cols).copy()

model_base = model_base[
    (model_base['site_name'] != "") &
    (model_base['water_class'].isin(ecosystem_classes)) &
    np.isfinite(model_base['Trans_ratio']) &
    (model_base['Trans_ratio'] >= 0) &
    (model_base['Trans_ratio'] <= 1) &
    np.isfinite(model_base['WUE']) &
    np.isfinite(model_base['WUE_tra']) &
    np.isfinite(model_base['WUE_eva'])
].copy()

print(f"  After filtering: {len(model_base)} rows")

# ============================================================================
# CHECK FOR DUPLICATES (IMPORTANT FOR R reshape())
# ============================================================================

dup_count = model_base.duplicated(subset=["site_name", "Year", "month"]).sum()
print(f"\nDuplicate site-Year-month rows: {dup_count}")
if dup_count > 0:
    print("  WARNING: R reshape() with idvar=c('site_name','Year','month') may behave unexpectedly")
    print("  Consider investigating duplicates before proceeding")
else:
    print("  No duplicates found - safe for R reshape()")

# ============================================================================
# STEP 2: SAVE WIDE DATA FOR R (IN TEMP FOLDER)
# ============================================================================

temp_model_base_file = os.path.join(temp_dir, "temp_model_base_for_R.csv")
model_base.to_csv(temp_model_base_file, index=False)
print(f"\nWIDE data saved to: {temp_model_base_file}")
print(f"  Shape: {model_base.shape} (rows, columns)")

# Record original input path for summary
original_input_path = monthly_data_path.replace("\\", "/")

# ============================================================================
# STEP 3: CREATE R SCRIPT (R DOES THE RESHAPE AND MODELING)
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Creating R Script for Modeling")
print("="*60)

# Use forward slashes for R
output_dir_r = output_dir.replace("\\", "/")
temp_model_base_file_r = temp_model_base_file.replace("\\", "/")
figure_dir_r = figure_dir.replace("\\", "/")  # Not used in R but kept for consistency

r_script_content = f'''
# Simplified Q1 WUE metric-difference models
# R handles: reshaping to long, lme4 models, mgcv GAMs, predictions
#
# Models:
#   Mixed: WUE_difference ~ difference_type * water_class * Trans_ratio_z +
#           (1 | site_name) + (1 | month_f)
#   GAM:   WUE_difference ~ difference_type * water_class +
#           s(Trans_ratio, by = difference_type:water_class) +
#           s(site_name, bs = "re") + s(month_f, bs = "re")

suppressPackageStartupMessages({{
  if (!requireNamespace("lme4", quietly = TRUE)) {{
    stop("Package 'lme4' is required.")
  }}
  if (!requireNamespace("mgcv", quietly = TRUE)) {{
    stop("Package 'mgcv' is required.")
  }}
  if (!requireNamespace("dplyr", quietly = TRUE)) {{
    stop("Package 'dplyr' is required.")
  }}
}})

# EXACT MALONE SIG FUNCTION
sig <- function(p) {{
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
}}

# Load WIDE data (Python cleaned, NOT reshaped yet)
model_base <- read.csv("{temp_model_base_file_r}", stringsAsFactors = FALSE)

# Define ecosystem classes
ecosystem_classes <- c("Upland", "Freshwater", "Saline")
ecosystem_colors <- c(Upland = "#800080", Freshwater = "#0000FF", Saline = "#FFA500")

# Difference labels
diff_labels <- c(
  WUE_ET_minus_WUE_T = "WUE_ET - WUE_T",
  WUE_ET_minus_WUE_E = "WUE_ET - WUE_E",
  WUE_E_minus_WUE_T = "WUE_E - WUE_T"
)
difference_cols <- c("WUE_ET_minus_WUE_T", "WUE_ET_minus_WUE_E", "WUE_E_minus_WUE_T")
diff_label_levels <- unname(diff_labels[difference_cols])

# ---------------------------------------------------------------------------
# RESHAPE TO LONG FORMAT - R DOES THIS ONCE
# ---------------------------------------------------------------------------

model_base$WUE_ET_minus_WUE_T <- model_base$WUE - model_base$WUE_tra
model_base$WUE_ET_minus_WUE_E <- model_base$WUE - model_base$WUE_eva
model_base$WUE_E_minus_WUE_T <- model_base$WUE_eva - model_base$WUE_tra

diff_data <- reshape(
  model_base,
  varying = difference_cols,
  v.names = "WUE_difference",
  timevar = "difference_type",
  times = difference_cols,
  idvar = c("site_name", "Year", "month"),
  direction = "long"
)

# Clean up and set factors
diff_data <- diff_data[complete.cases(diff_data[, c("site_name", "Year", "month", "Trans_ratio", "WUE_difference")]), ]
diff_data <- diff_data[is.finite(diff_data$WUE_difference), ]
diff_data$site_name <- factor(diff_data$site_name)
diff_data$month_f <- factor(diff_data$month, levels = 1:12)
diff_data$water_class <- factor(diff_data$water_class, levels = ecosystem_classes)
diff_data$difference_type <- factor(diff_data$difference_type, levels = difference_cols)
diff_data$Trans_ratio_z <- as.numeric(scale(diff_data$Trans_ratio))
diff_data$diff_ecosystem <- interaction(diff_data$difference_type, diff_data$water_class, sep = "__", drop = TRUE)
diff_data$difference_label <- factor(diff_labels[as.character(diff_data$difference_type)], levels = diff_label_levels)

cat("\nR reshape complete:")
cat("\n  Long data rows:", nrow(diff_data))
cat("\n  Sites:", length(unique(diff_data$site_name)))
cat("\n  Months:", length(unique(diff_data$month_f)), "\n")

# ---------------------------------------------------------------------------
# MIXED MODEL - EXACT MALONE
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

# Save RDS
saveRDS(simple_mixed, file.path("{output_dir_r}", "Q1_simple_mixed_difference_TET_model.rds"))

# Fitted values
diff_data$simple_mixed_fitted <- fitted(simple_mixed)
diff_data$simple_mixed_residual <- residuals(simple_mixed)
write.csv(
  diff_data[, c("site_name", "Year", "month", "water_class", "difference_type",
                "Trans_ratio", "Trans_ratio_z", "WUE_difference",
                "simple_mixed_fitted", "simple_mixed_residual")],
  file.path("{output_dir_r}", "Q1_simple_mixed_fitted_residuals.csv"),
  row.names = FALSE
)

# Fixed effects
mixed_fixed <- as.data.frame(summary(simple_mixed)$coefficients)
mixed_fixed$term <- rownames(mixed_fixed)
rownames(mixed_fixed) <- NULL
mixed_fixed <- mixed_fixed[, c("term", setdiff(names(mixed_fixed), "term"))]
names(mixed_fixed) <- gsub(" ", "_", names(mixed_fixed), fixed = TRUE)
mixed_fixed$p_value_normal_approx <- 2 * stats::pnorm(-abs(mixed_fixed$t_value))
mixed_fixed$signif <- sig(mixed_fixed$p_value_normal_approx)
write.csv(
  mixed_fixed,
  file.path("{output_dir_r}", "Q1_simple_mixed_fixed_effects_with_approx_p.csv"),
  row.names = FALSE
)

# LRT
mixed_lrt <- as.data.frame(drop1(simple_mixed, test = "Chisq"))
mixed_lrt$term <- rownames(mixed_lrt)
rownames(mixed_lrt) <- NULL
mixed_lrt <- mixed_lrt[, c("term", setdiff(names(mixed_lrt), "term"))]
names(mixed_lrt) <- gsub(" ", "_", names(mixed_lrt), fixed = TRUE)
p_col <- grep("Pr", names(mixed_lrt), value = TRUE)[1]
mixed_lrt$signif <- sig(mixed_lrt[[p_col]])
write.csv(
  mixed_lrt,
  file.path("{output_dir_r}", "Q1_simple_mixed_likelihood_ratio_tests.csv"),
  row.names = FALSE
)

# Model comparison
mixed_comparison <- as.data.frame(anova(simple_mixed_additive, simple_mixed))
p_col <- grep("Pr", names(mixed_comparison), value = TRUE)[1]
mixed_comparison$signif <- sig(mixed_comparison[[p_col]])
write.csv(
  mixed_comparison,
  file.path("{output_dir_r}", "Q1_simple_mixed_model_comparison_lrt.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# GAM MODELS - EXACT MALONE
# ---------------------------------------------------------------------------

gam_data <- diff_data[is.finite(diff_data$SPEI_1) & diff_data$SPEI_1 >= -1 & diff_data$SPEI_1 <= 1, ]

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

saveRDS(simple_smooth_gam, file.path("{output_dir_r}", "Q1_simple_smooth_gam_difference_TET_model.rds"))

# Fitted values
gam_data$simple_gam_fitted <- fitted(simple_smooth_gam)
gam_data$simple_gam_residual <- residuals(simple_smooth_gam)
write.csv(
  gam_data[, c("site_name", "Year", "month", "water_class", "difference_type",
               "Trans_ratio", "Trans_ratio_z", "SPEI_1", "WUE_difference",
               "simple_gam_fitted", "simple_gam_residual")],
  file.path("{output_dir_r}", "Q1_simple_gam_fitted_residuals.csv"),
  row.names = FALSE
)

# GAM comparison
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
  file.path("{output_dir_r}", "Q1_simple_gam_model_comparison.csv"),
  row.names = FALSE
)

# Parametric terms
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
  file.path("{output_dir_r}", "Q1_simple_gam_parametric_terms.csv"),
  row.names = FALSE
)

# Smooth terms
gam_smooth <- as.data.frame(summary(simple_smooth_gam)$s.table)
gam_smooth$smooth_term <- rownames(gam_smooth)
rownames(gam_smooth) <- NULL
gam_smooth <- gam_smooth[, c("smooth_term", setdiff(names(gam_smooth), "smooth_term"))]
names(gam_smooth) <- gsub("p-value", "p_value", names(gam_smooth), fixed = TRUE)
gam_smooth$signif <- sig(gam_smooth$p_value)
write.csv(
  gam_smooth,
  file.path("{output_dir_r}", "Q1_simple_gam_smooth_terms.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# PREDICTIONS - EXACT MALONE
# ---------------------------------------------------------------------------

tet_sequence <- seq(quantile(diff_data$Trans_ratio, 0.05), quantile(diff_data$Trans_ratio, 0.95), length.out = 100)
tet_mean <- mean(diff_data$Trans_ratio)
tet_sd <- sd(diff_data$Trans_ratio)

# Mixed predictions
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
mixed_grid$difference_label <- factor(diff_labels[as.character(mixed_grid$difference_type)], levels = diff_label_levels)
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
  file.path("{output_dir_r}", "Q1_simple_mixed_predictions_TET.csv"),
  row.names = FALSE
)

# GAM predictions
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
gam_grid$Trans_ratio_z <- (gam_grid$Trans_ratio - mean(gam_data$Trans_ratio)) / sd(gam_data$Trans_ratio)
gam_grid$difference_label <- factor(diff_labels[as.character(gam_grid$difference_type)], levels = diff_label_levels)
gam_grid$ecosystem_label <- as.character(gam_grid$water_class)
gam_grid$diff_ecosystem <- interaction(gam_grid$difference_type, gam_grid$water_class, sep = "__", drop = TRUE)
gam_grid$diff_ecosystem <- factor(gam_grid$diff_ecosystem, levels = levels(gam_data$diff_ecosystem))

gam_pred <- predict(simple_smooth_gam, newdata = gam_grid, type = "link", se.fit = TRUE, 
                    exclude = c("s(site_name)", "s(month_f)"))
gam_grid$prediction <- as.numeric(gam_pred$fit)
gam_grid$se <- as.numeric(gam_pred$se.fit)
gam_grid$lower <- gam_grid$prediction - 1.96 * gam_grid$se
gam_grid$upper <- gam_grid$prediction + 1.96 * gam_grid$se
write.csv(
  gam_grid,
  file.path("{output_dir_r}", "Q1_simple_gam_predictions_TET.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SUMMARY STATISTICS
# ---------------------------------------------------------------------------

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
  file.path("{output_dir_r}", "Q1_simple_summary_by_difference_type.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# MODEL SUMMARY - FULL MALONE-STYLE
# ---------------------------------------------------------------------------

sink(file.path("{output_dir_r}", "Q1_simple_model_summary.txt"))
cat("Q1 difference type + T:ET + ecosystem models\n")
cat("Original input:", "{original_input_path}", "\n")
cat("R model input:", "{temp_model_base_file_r}", "\n")
cat("Output:", "{output_dir_r}", "\n\n")
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

cat("\nR analysis complete! All outputs saved to:", "{output_dir_r}", "\n")
'''

# Save R script in temp folder
r_script_file = os.path.join(temp_dir, "run_models.R")
with open(r_script_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print(f"  R script saved to: {r_script_file}")

# Verify temp files exist
print("\nVerifying temp files exist before running R:")
print(f"  temp_model_base_for_R.csv exists: {os.path.exists(temp_model_base_file)}")
print(f"  run_models.R exists: {os.path.exists(r_script_file)}")

# ============================================================================
# STEP 4: CLEAN OLD OUTPUTS AND FIGURES BEFORE RUNNING R
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Cleaning old outputs and figures (Python directories only)")
print("="*60)

# Clean R outputs
expected_r_outputs = [
    "Q1_simple_mixed_difference_TET_model.rds",
    "Q1_simple_smooth_gam_difference_TET_model.rds",
    "Q1_simple_mixed_fitted_residuals.csv",
    "Q1_simple_mixed_fixed_effects_with_approx_p.csv",
    "Q1_simple_mixed_likelihood_ratio_tests.csv",
    "Q1_simple_mixed_model_comparison_lrt.csv",
    "Q1_simple_gam_fitted_residuals.csv",
    "Q1_simple_gam_model_comparison.csv",
    "Q1_simple_gam_parametric_terms.csv",
    "Q1_simple_gam_smooth_terms.csv",
    "Q1_simple_mixed_predictions_TET.csv",
    "Q1_simple_gam_predictions_TET.csv",
    "Q1_simple_summary_by_difference_type.csv",
    "Q1_simple_model_summary.txt",
]

for f in expected_r_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

# Clean figure files
expected_figures = [
    "Q1_simple_plot_01_mean_difference_type.png",
    "Q1_simple_plot_02_mixed_TET_predictions.png",
    "Q1_simple_plot_03_gam_TET_predictions.png",
    "Q1_simple_plot_04_gam_model_comparison.png",
    "Q1_simple_plot_05_mixed_residuals.png",
    "Q1_simple_plot_06_gam_residuals.png",
    "Q1_simple_plot_07_manuscript_multipanel.png",
]

for f in expected_figures:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old figure: {f}")

# ============================================================================
# STEP 5: RUN R SCRIPT
# ============================================================================

print("\n" + "="*60)
print("STEP 4: Running R Modeling")
print("="*60)

# Add R to PATH
r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if r_path not in os.environ['PATH']:
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + r_path

try:
    result = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)
    print("R is available")
except Exception as e:
    print(f"R not found! Error: {e}")
    exit()

# Verify temp files still exist before running R
print("\nVerifying temp files still exist before R execution:")
print(f"  temp_model_base_for_R.csv exists: {os.path.exists(temp_model_base_file)}")
print(f"  run_models.R exists: {os.path.exists(r_script_file)}")

if not os.path.exists(temp_model_base_file):
    raise FileNotFoundError(f"Temp data file missing: {temp_model_base_file}")
if not os.path.exists(r_script_file):
    raise FileNotFoundError(f"R script file missing: {r_script_file}")

result = subprocess.run(["Rscript", r_script_file], capture_output=True, text=True)

if result.returncode != 0:
    print("R script failed!")
    print("\nSTDERR:")
    print(result.stderr)
    raise RuntimeError("Stopping because R analysis failed.")
else:
    print("R analysis completed successfully!")
    if result.stdout:
        print("\nR output:")
        print(result.stdout)

# ============================================================================
# STEP 6: VERIFY R OUTPUTS
# ============================================================================

print("\n" + "="*60)
print("STEP 5: Verifying R Outputs")
print("="*60)

missing_outputs = [
    f for f in expected_r_outputs
    if not os.path.exists(os.path.join(output_dir, f))
]

if missing_outputs:
    raise FileNotFoundError(
        "R finished, but these expected outputs are missing:\n" +
        "\n".join(missing_outputs)
    )
else:
    print(f"All {len(expected_r_outputs)} expected R outputs were created.")

# ============================================================================
# STEP 7: PYTHON READS R OUTPUT AND CREATES FIGURES
# ============================================================================

print("\n" + "="*60)
print("STEP 6: Creating Python Figures from R Outputs")
print("="*60)

# Load R outputs
print("\nLoading R output CSVs...")
summary_df = pd.read_csv(os.path.join(output_dir, "Q1_simple_summary_by_difference_type.csv"))
mixed_pred = pd.read_csv(os.path.join(output_dir, "Q1_simple_mixed_predictions_TET.csv"))
gam_pred = pd.read_csv(os.path.join(output_dir, "Q1_simple_gam_predictions_TET.csv"))
mixed_fitted = pd.read_csv(os.path.join(output_dir, "Q1_simple_mixed_fitted_residuals.csv"))
gam_fitted = pd.read_csv(os.path.join(output_dir, "Q1_simple_gam_fitted_residuals.csv"))
gam_comp = pd.read_csv(os.path.join(output_dir, "Q1_simple_gam_model_comparison.csv"))

# Define ecosystem colors
ecosystem_colors = {"Upland": "#800080", "Freshwater": "#0000FF", "Saline": "#FFA500"}
ecosystem_classes = ["Upland", "Freshwater", "Saline"]
diff_label_levels = ["WUE_ET - WUE_T", "WUE_ET - WUE_E", "WUE_E - WUE_T"]

print(f"  Summary: {len(summary_df)} rows")
print(f"  Mixed predictions: {len(mixed_pred)} rows")
print(f"  GAM predictions: {len(gam_pred)} rows")

# Apply theme
theme_wue_manuscript()

# ============================================================================
# FIGURE 1: Summary plot
# ============================================================================

print("\nCreating Plot 1: Summary...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes[idx]
    subset = summary_df[summary_df['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['water_class'] == ecosystem]
        if len(eco_subset) > 0:
            ax.errorbar(ecosystem, eco_subset['mean_difference'].values[0],
                       yerr=eco_subset['ci95_difference'].values[0],
                       fmt='o', color=ecosystem_colors[ecosystem], capsize=8,
                       markersize=12, elinewidth=3)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=16, fontweight='bold')
    ax.set_xlabel('Ecosystem', fontsize=14)
    if idx == 0: ax.set_ylabel('Mean WUE metric difference', fontsize=14)
    ax.tick_params(axis='x', rotation=30, labelsize=12)
    ax.grid(True, alpha=0.3)
plt.suptitle('Q1 Model Inputs: Mean WUE Metric Differences by Ecosystem', fontsize=18, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_01_mean_difference_type.png", width=10, height=5)

# ============================================================================
# FIGURE 2: Mixed predictions
# ============================================================================

print("Creating Plot 2: Mixed predictions...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes[idx]
    subset = mixed_pred[mixed_pred['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['ecosystem_label'] == ecosystem]
        if len(eco_subset) > 0:
            ax.fill_between(eco_subset['Trans_ratio'], eco_subset['lower'], eco_subset['upper'],
                           color=ecosystem_colors[ecosystem], alpha=0.15)
            ax.plot(eco_subset['Trans_ratio'], eco_subset['prediction'],
                   color=ecosystem_colors[ecosystem], linewidth=2.5, label=ecosystem)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=16, fontweight='bold')
    ax.set_xlabel('T:ET ratio (Trans_ratio)', fontsize=14)
    if idx == 0: ax.set_ylabel('Predicted WUE metric difference', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.3)
    if idx == 2: ax.legend(loc='best', fontsize=12, frameon=True)
plt.suptitle('Mixed Model: T:ET Response by Ecosystem', fontsize=18, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_02_mixed_TET_predictions.png", width=10, height=6)

# ============================================================================
# FIGURE 3: GAM predictions
# ============================================================================

print("Creating Plot 3: GAM predictions...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes[idx]
    subset = gam_pred[gam_pred['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['ecosystem_label'] == ecosystem]
        if len(eco_subset) > 0:
            ax.fill_between(eco_subset['Trans_ratio'], eco_subset['lower'], eco_subset['upper'],
                           color=ecosystem_colors[ecosystem], alpha=0.15)
            ax.plot(eco_subset['Trans_ratio'], eco_subset['prediction'],
                   color=ecosystem_colors[ecosystem], linewidth=2.5, label=ecosystem)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=16, fontweight='bold')
    ax.set_xlabel('T:ET ratio (Trans_ratio)', fontsize=14)
    if idx == 0: ax.set_ylabel('Predicted WUE metric difference', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.3)
    if idx == 2: ax.legend(loc='best', fontsize=12, frameon=True)
plt.suptitle('Smooth GAM: T:ET Response by Ecosystem', fontsize=18, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_03_gam_TET_predictions.png", width=10, height=6)

# ============================================================================
# FIGURE 4: Model comparison
# ============================================================================

print("Creating Plot 4: Model comparison...")
fig, ax = plt.subplots(figsize=(8, 6))
models = ['Linear GAM', 'Smooth GAM']
deviance = [gam_comp.loc[0, 'deviance_explained'] * 100, 
            gam_comp.loc[1, 'deviance_explained'] * 100]
colors = ['#8DA0CB', '#66C2A5']
bars = ax.bar(models, deviance, color=colors, edgecolor='#4D4D4D', linewidth=1.5)
for bar, val in zip(bars, deviance):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
            f'{val:.1f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')
ax.set_ylabel('Deviance explained (%)', fontsize=14)
ax.set_title('Simple Q1 GAM Model Comparison', fontsize=16, fontweight='bold')
ax.set_ylim(0, 95)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3, axis='y')
delta_aic = gam_comp.loc[1, 'delta_AIC'] - gam_comp.loc[0, 'delta_AIC']
ax.text(0.5, -0.15, f'Smooth GAM improves AIC by {abs(delta_aic):.1f}',
        transform=ax.transAxes, ha='center', va='center', fontsize=12)
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_04_gam_model_comparison.png", width=8, height=6)

# ============================================================================
# FIGURE 5: Mixed residuals
# ============================================================================

print("Creating Plot 5: Mixed residuals...")
mixed_fitted['difference_label'] = mixed_fitted['difference_type'].map(
    {'WUE_ET_minus_WUE_T': 'WUE_ET - WUE_T',
     'WUE_ET_minus_WUE_E': 'WUE_ET - WUE_E',
     'WUE_E_minus_WUE_T': 'WUE_E - WUE_T'}
)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes[idx]
    subset = mixed_fitted[mixed_fitted['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['water_class'] == ecosystem]
        if len(eco_subset) > 0:
            ax.scatter(eco_subset['simple_mixed_fitted'], eco_subset['simple_mixed_residual'],
                      color=ecosystem_colors[ecosystem], alpha=0.3, s=30)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=16, fontweight='bold')
    ax.set_xlabel('Fitted difference', fontsize=14)
    if idx == 0: ax.set_ylabel('Residual', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.3)
plt.suptitle('Simple Mixed Model Residuals', fontsize=18, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_05_mixed_residuals.png", width=10, height=6)

# ============================================================================
# FIGURE 6: GAM residuals
# ============================================================================

print("Creating Plot 6: GAM residuals...")
gam_fitted['difference_label'] = gam_fitted['difference_type'].map(
    {'WUE_ET_minus_WUE_T': 'WUE_ET - WUE_T',
     'WUE_ET_minus_WUE_E': 'WUE_ET - WUE_E',
     'WUE_E_minus_WUE_T': 'WUE_E - WUE_T'}
)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes[idx]
    subset = gam_fitted[gam_fitted['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['water_class'] == ecosystem]
        if len(eco_subset) > 0:
            ax.scatter(eco_subset['simple_gam_fitted'], eco_subset['simple_gam_residual'],
                      color=ecosystem_colors[ecosystem], alpha=0.3, s=30)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=16, fontweight='bold')
    ax.set_xlabel('Fitted difference', fontsize=14)
    if idx == 0: ax.set_ylabel('Residual', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.grid(True, alpha=0.3)
plt.suptitle('Simple Smooth GAM Residuals', fontsize=18, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q1_simple_plot_06_gam_residuals.png", width=10, height=6)

# ============================================================================
# FIGURE 7: Manuscript multi-panel
# ============================================================================

print("Creating Plot 7: Manuscript multi-panel...")

fig = plt.figure(figsize=(24, 25))
outer = gridspec.GridSpec(4, 1, height_ratios=[0.9, 0.8, 1.25, 1.25],
                          hspace=0.35, top=0.95, bottom=0.05)

# Row A: Summary
gs_a = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0], wspace=0.25)
axes_a = [fig.add_subplot(gs_a[0, i]) for i in range(3)]

# Row B: Model comparison
gs_b = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[1])
ax_b = fig.add_subplot(gs_b[0, 0])

# Row C: Mixed predictions
gs_c = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2], wspace=0.25)
axes_c = [fig.add_subplot(gs_c[0, i]) for i in range(3)]

# Row D: GAM predictions
gs_d = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[3], wspace=0.25)
axes_d = [fig.add_subplot(gs_d[0, i]) for i in range(3)]

# Row A
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes_a[idx]
    subset = summary_df[summary_df['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['water_class'] == ecosystem]
        if len(eco_subset) > 0:
            ax.errorbar(ecosystem, eco_subset['mean_difference'].values[0],
                       yerr=eco_subset['ci95_difference'].values[0],
                       fmt='o', color=ecosystem_colors[ecosystem], capsize=5,
                       markersize=10, elinewidth=2.5)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=22, fontweight='bold')
    ax.set_xlabel('Ecosystem', fontsize=18)
    if idx == 0: ax.set_ylabel('Mean WUE metric difference', fontsize=18)
    ax.tick_params(axis='x', rotation=25, labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid(True, alpha=0.3)
axes_a[0].text(-0.15, 1.05, 'A', transform=axes_a[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

# Row B
models = ['Linear GAM', 'Smooth GAM']
deviance = [gam_comp.loc[0, 'deviance_explained'] * 100, 
            gam_comp.loc[1, 'deviance_explained'] * 100]
bars = ax_b.bar(models, deviance, color=['#8DA0CB', '#66C2A5'], 
                edgecolor='#4D4D4D', linewidth=2)
for bar, val in zip(bars, deviance):
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f'{val:.1f}%', ha='center', va='bottom', fontsize=20, fontweight='bold')
ax_b.set_ylabel('Deviance explained (%)', fontsize=22)
ax_b.set_title('GAM Model Comparison', fontsize=25, fontweight='bold')
ax_b.set_ylim(0, 95)
ax_b.tick_params(labelsize=18)
ax_b.grid(True, alpha=0.3, axis='y')
delta_aic = gam_comp.loc[1, 'delta_AIC'] - gam_comp.loc[0, 'delta_AIC']
ax_b.text(0.5, -0.12, f'Smooth GAM improves AIC by {abs(delta_aic):.1f}',
         transform=ax_b.transAxes, ha='center', va='center', fontsize=18)
ax_b.text(-0.08, 1.05, 'B', transform=ax_b.transAxes, 
          fontsize=38, fontweight='bold', va='top')

# Row C
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes_c[idx]
    subset = mixed_pred[mixed_pred['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['ecosystem_label'] == ecosystem]
        if len(eco_subset) > 0:
            ax.fill_between(eco_subset['Trans_ratio'], eco_subset['lower'], eco_subset['upper'],
                           color=ecosystem_colors[ecosystem], alpha=0.15)
            ax.plot(eco_subset['Trans_ratio'], eco_subset['prediction'],
                   color=ecosystem_colors[ecosystem], linewidth=2.5, label=ecosystem)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=22, fontweight='bold')
    ax.set_xlabel('T:ET ratio (Trans_ratio)', fontsize=18)
    if idx == 0: ax.set_ylabel('Predicted WUE metric difference', fontsize=18)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    if idx == 2: ax.legend(loc='best', fontsize=16, frameon=True)
axes_c[0].text(-0.15, 1.05, 'C', transform=axes_c[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

# Row D
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes_d[idx]
    subset = gam_pred[gam_pred['difference_label'] == diff_type]
    for ecosystem in ecosystem_classes:
        eco_subset = subset[subset['ecosystem_label'] == ecosystem]
        if len(eco_subset) > 0:
            ax.fill_between(eco_subset['Trans_ratio'], eco_subset['lower'], eco_subset['upper'],
                           color=ecosystem_colors[ecosystem], alpha=0.15)
            ax.plot(eco_subset['Trans_ratio'], eco_subset['prediction'],
                   color=ecosystem_colors[ecosystem], linewidth=2.5, label=ecosystem)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.set_title(diff_type, fontsize=22, fontweight='bold')
    ax.set_xlabel('T:ET ratio (Trans_ratio)', fontsize=18)
    if idx == 0: ax.set_ylabel('Predicted WUE metric difference', fontsize=18)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    if idx == 2: ax.legend(loc='best', fontsize=16, frameon=True)
axes_d[0].text(-0.15, 1.05, 'D', transform=axes_d[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

fig.suptitle('Q1 WUE Metric Difference Models', fontsize=42, fontweight='bold', y=0.98)
save_plot(fig, "Q1_simple_plot_07_manuscript_multipanel.png", width=24, height=25)

# ============================================================================
# STEP 8: VERIFY FIGURES WERE CREATED
# ============================================================================

print("\n" + "="*60)
print("STEP 7: Verifying Figures")
print("="*60)

missing_figures = [
    f for f in expected_figures
    if not os.path.exists(os.path.join(figure_dir, f))
]

if missing_figures:
    raise FileNotFoundError(
        "These expected figures are missing:\n" + "\n".join(missing_figures)
    )
else:
    print(f"All {len(expected_figures)} expected figures were created.")

# ============================================================================
# STEP 9: CREATE README
# ============================================================================

print("\n" + "="*60)
print("STEP 8: Creating README")
print("="*60)

readme_content = """# Q1 difference type + T:ET + ecosystem outputs (Python Replication)

This is a Python-replicated version of the Malone Q1 analysis. 
The models are run in R (lme4, mgcv) for exact statistical replication.
Figures are rendered in Python using matplotlib/seaborn.

## Directory Structure
- `Q1_mixed_model_outputs/` - CSV files, RDS files, README, model summary
- `Q1_mixed_model_figures/` - PNG figure files
- `Q1_mixed_model_temp/` - Temporary R handoff files (cleaned at start)

## Workflow
- Python: Data cleaning, filtering, duplicate checks, orchestrating R
- R: Reshaping to long format, lme4 mixed model, mgcv GAM, predictions
- Python: Reading R outputs, creating figures, generating README

## Core models
- Mixed model: `WUE_difference ~ difference_type * water_class * Trans_ratio_z + (1 | site_name) + (1 | month_f)`
- Near-normal smooth GAM: `WUE_difference ~ difference_type * water_class + s(Trans_ratio, by = difference_type:water_class) + s(site_name, bs = "re") + s(month_f, bs = "re")`

## Key tables (in outputs folder)
- `Q1_simple_summary_by_difference_type.csv`
- `Q1_simple_mixed_fixed_effects_with_approx_p.csv`
- `Q1_simple_mixed_likelihood_ratio_tests.csv`
- `Q1_simple_mixed_model_comparison_lrt.csv`
- `Q1_simple_gam_parametric_terms.csv`
- `Q1_simple_gam_smooth_terms.csv`
- `Q1_simple_gam_model_comparison.csv`
- `Q1_simple_mixed_predictions_TET.csv`
- `Q1_simple_gam_predictions_TET.csv`
- `Q1_simple_mixed_fitted_residuals.csv`
- `Q1_simple_gam_fitted_residuals.csv`
- `Q1_simple_model_summary.txt`

## RDS files (in outputs folder)
- `Q1_simple_mixed_difference_TET_model.rds`
- `Q1_simple_smooth_gam_difference_TET_model.rds`

## Figures (in figures folder)
Python-rendered equivalent figures using R-generated model outputs.
Uses matplotlib/seaborn with Malone-style themes and color schemes.

- `Q1_simple_plot_01_mean_difference_type.png`
- `Q1_simple_plot_02_mixed_TET_predictions.png`
- `Q1_simple_plot_03_gam_TET_predictions.png`
- `Q1_simple_plot_04_gam_model_comparison.png`
- `Q1_simple_plot_05_mixed_residuals.png`
- `Q1_simple_plot_06_gam_residuals.png`
- `Q1_simple_plot_07_manuscript_multipanel.png`

## Note on figures
These are Python-rendered equivalent figures using R-generated model outputs. 
They use matplotlib/seaborn with Malone-style themes and color schemes.
"""

with open(os.path.join(output_dir, "Q1_simple_outputs_README.md"), 'w') as f:
    f.write(readme_content)

print("  Created: Q1_simple_outputs_README.md")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("COMPLETE! Malone Workflow Replicated in Python")
print("="*60)

print(f"\nPython outputs saved to:")
print(f"  Outputs (CSV/RDS/README): {output_dir}")
print(f"  Figures (PNG):           {figure_dir}")
print(f"  Temp files:              {temp_dir}")

# List output files
output_files = sorted([f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))])
print(f"\nOutput files ({len(output_files)} files):")
for f in output_files:
    size = os.path.getsize(os.path.join(output_dir, f))
    print(f"  {f} ({size:,} bytes)")

# List figure files
figure_files = sorted([f for f in os.listdir(figure_dir) if os.path.isfile(os.path.join(figure_dir, f))])
print(f"\nFigure files ({len(figure_files)} files):")
for f in figure_files:
    size = os.path.getsize(os.path.join(figure_dir, f))
    print(f"  {f} ({size:,} bytes)")

print("\n" + "="*60)
print("✓ Workflow complete - Python orchestrates, R models, Python plots")
print("  - No double reshaping (Python saves WIDE, R reshapes ONCE)")
print("  - Same model outputs as Malone")
print("  - Python figures with same filenames as Malone")
print("  - Clean directory structure (outputs/figures/temp separate)")
print("  - Temp cleaned at start (not after creating files)")
print("  - Outputs in separate directory (safe from Malone)")
print("="*60)