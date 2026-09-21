"""
COMPLETE Q1 FINAL ANALYSIS - Near-Normal SPEI-1 Workflow with Mixed Model Panel B
All outputs saved to Q1_updated_results directory structure

FINAL DECISION: All panels use near-normal hydroclimatic conditions (SPEI_1 >= -1 & SPEI_1 <= 1)
Main models: Mixed Model (Panel B) and Smooth GAM (Panel C)
Linear GAM retained only for Panel D model comparison

PANEL STRUCTURE:
  Panel A: Near-normal mean WUE metric differences by ecosystem
  Panel B: Near-normal mixed model T:ET response by ecosystem
  Panel C: Near-normal smooth GAM T:ET response by ecosystem
  Panel D: Linear GAM vs Smooth GAM model comparison

DIRECTORY STRUCTURE:
  Q1_updated_results/outputs/   - CSV files, README, model summary
  Q1_updated_results/figures/   - PNG figure files  
  Q1_updated_results/temp/      - Temporary R handoff files (cleaned at start)

OUTPUTS: M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\
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
# PATHS - UPDATED OUTPUT DIRECTORY (DO NOT OVERWRITE OLD)
# ============================================================================

# Input data (read from Malone location)
monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# Updated outputs - SEPARATE from old Q1 outputs
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"

# Three subdirectories for organization
output_dir = os.path.join(base_output_dir, "outputs")     # CSV, README, summary
figure_dir = os.path.join(base_output_dir, "figures")     # PNG figures
temp_dir = os.path.join(base_output_dir, "temp")          # Temporary R files

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir]:
    os.makedirs(dir_path, exist_ok=True)

print("="*60)
print("Q1 FINAL - Near-Normal SPEI-1 Workflow (Mixed Model Panel B)")
print("="*60)
print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"(Separate from old Q1 outputs - safe for comparison)")

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

# Required columns (including SPEI_1 for filtering)
required_cols = [
    "site_name", "Year", "month", "water_class", "Trans_ratio",
    "WUE", "WUE_tra", "WUE_eva", "SPEI_1"
]

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

ecosystem_classes = ["Upland", "Freshwater", "Saline"]

# Filter - exact Malone (keep SPEI_1 for later filtering in R)
model_base = monthly.dropna(subset=required_cols).copy()

model_base = model_base[
    (model_base['site_name'] != "") &
    (model_base['water_class'].isin(ecosystem_classes)) &
    np.isfinite(model_base['Trans_ratio']) &
    (model_base['Trans_ratio'] >= 0) &
    (model_base['Trans_ratio'] <= 1) &
    np.isfinite(model_base['WUE']) &
    np.isfinite(model_base['WUE_tra']) &
    np.isfinite(model_base['WUE_eva']) &
    np.isfinite(model_base['SPEI_1'])
].copy()

print(f"  After filtering: {len(model_base)} rows")

# ============================================================================
# CHECK FOR DUPLICATES (IMPORTANT FOR R reshape())
# ============================================================================

dup_count = model_base.duplicated(subset=["site_name", "Year", "month"]).sum()
print(f"\nDuplicate site-Year-month rows: {dup_count}")
if dup_count > 0:
    print("  WARNING: R reshape() with idvar=c('site_name','Year','month') may behave unexpectedly")
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
# STEP 3: CREATE R SCRIPT - FINAL WORKFLOW WITH MIXED MODEL PANEL B
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Creating R Script for Final Workflow (Mixed Model Panel B)")
print("="*60)

# Use forward slashes for R
output_dir_r = output_dir.replace("\\", "/")
temp_model_base_file_r = temp_model_base_file.replace("\\", "/")

r_script_content = f'''
# Q1 FINAL - Near-Normal SPEI-1 Workflow with Mixed Model Panel B
# R handles: reshaping to long, near-normal filtering, mixed model, GAM models, predictions
#
# PANEL B (near-normal SPEI-1 only):
#   Mixed Model: WUE_difference ~ difference_type * water_class * Trans_ratio_z_NN +
#                (1 | site_name) + (1 | month_f)
#
# PANEL C (near-normal SPEI-1 only):
#   Smooth GAM: WUE_difference ~ difference_type * water_class +
#               s(Trans_ratio, by = diff_ecosystem, k = 6) +
#               s(site_name, bs = "re") + s(month_f, bs = "re")
#
# PANEL D (near-normal SPEI-1 only, model comparison):
#   Linear GAM: WUE_difference ~ difference_type * water_class * Trans_ratio +
#               s(site_name, bs = "re") + s(month_f, bs = "re")

suppressPackageStartupMessages({{
  if (!requireNamespace("mgcv", quietly = TRUE)) {{
    stop("Package 'mgcv' is required.")
  }}
  if (!requireNamespace("lme4", quietly = TRUE)) {{
    stop("Package 'lme4' is required.")
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
diff_data <- diff_data[complete.cases(diff_data[, c("site_name", "Year", "month", "Trans_ratio", "WUE_difference", "SPEI_1")]), ]
diff_data <- diff_data[is.finite(diff_data$WUE_difference), ]
diff_data$site_name <- factor(diff_data$site_name)
diff_data$month_f <- factor(diff_data$month, levels = 1:12)
diff_data$water_class <- factor(diff_data$water_class, levels = ecosystem_classes)
diff_data$difference_type <- factor(diff_data$difference_type, levels = difference_cols)
diff_data$diff_ecosystem <- interaction(diff_data$difference_type, diff_data$water_class, sep = "__", drop = TRUE)
diff_data$difference_label <- factor(diff_labels[as.character(diff_data$difference_type)], levels = diff_label_levels)

cat("\nR reshape complete:")
cat("\n  Long data rows:", nrow(diff_data))
cat("\n  Sites:", length(unique(diff_data$site_name)))
cat("\n  Months:", length(unique(diff_data$month_f)), "\n")

# ---------------------------------------------------------------------------
# FILTER TO NEAR-NORMAL SPEI-1 (FINAL Q1 DECISION)
# ---------------------------------------------------------------------------

q1_data <- diff_data[
  is.finite(diff_data$SPEI_1) &
  diff_data$SPEI_1 >= -1 &
  diff_data$SPEI_1 <= 1,
]

# Drop unused factor levels to avoid prediction issues
q1_data <- droplevels(q1_data)
q1_data$diff_ecosystem <- interaction(
  q1_data$difference_type,
  q1_data$water_class,
  sep = "__",
  drop = TRUE
)

# For mixed model only: standardized Trans_ratio
q1_data$Trans_ratio_z_NN <- as.numeric(scale(q1_data$Trans_ratio))

cat("\nNear-normal SPEI-1 subset (SPEI_1 >= -1 & SPEI_1 <= 1):")
cat("\n  Rows:", nrow(q1_data))
cat("\n  Sites:", length(unique(q1_data$site_name)))
cat("\n  Months:", length(unique(q1_data$month_f)), "\n")

# =============================================================================
# PANEL A: SUMMARY STATISTICS (near-normal SPEI-1 only)
# =============================================================================

summary_by_difference <- q1_data |>
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
  file.path("{output_dir_r}", "Q1_final_near_normal_summary_by_difference_type.csv"),
  row.names = FALSE
)

# =============================================================================
# PANEL B: NEAR-NORMAL MIXED MODEL (Main Panel B)
# =============================================================================

cat("\n----------------------------------------------------------------")
cat("\nPANEL B: Fitting near-normal mixed model...")

# Fit mixed model with interaction
q1_nn_mixed <- lme4::lmer(
  WUE_difference ~ difference_type * water_class * Trans_ratio_z_NN +
    (1 | site_name) + (1 | month_f),
  data = q1_data,
  REML = FALSE
)

# Additive model for interaction check
q1_nn_mixed_additive <- lme4::lmer(
  WUE_difference ~ difference_type + water_class + Trans_ratio_z_NN +
    (1 | site_name) + (1 | month_f),
  data = q1_data,
  REML = FALSE
)

# Save RDS
saveRDS(q1_nn_mixed, file.path("{output_dir_r}", "Q1_final_near_normal_mixed_difference_TET_model.rds"))

# Fixed effects
mixed_fixed <- as.data.frame(summary(q1_nn_mixed)$coefficients)
mixed_fixed$term <- rownames(mixed_fixed)
rownames(mixed_fixed) <- NULL
mixed_fixed <- mixed_fixed[, c("term", setdiff(names(mixed_fixed), "term"))]
names(mixed_fixed) <- gsub(" ", "_", names(mixed_fixed), fixed = TRUE)
mixed_fixed$p_value_normal_approx <- 2 * stats::pnorm(-abs(mixed_fixed$t_value))
mixed_fixed$signif <- sig(mixed_fixed$p_value_normal_approx)
write.csv(
  mixed_fixed,
  file.path("{output_dir_r}", "Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv"),
  row.names = FALSE
)

# Drop1 LRT
mixed_lrt <- as.data.frame(drop1(q1_nn_mixed, test = "Chisq"))
mixed_lrt$term <- rownames(mixed_lrt)
rownames(mixed_lrt) <- NULL
mixed_lrt <- mixed_lrt[, c("term", setdiff(names(mixed_lrt), "term"))]
names(mixed_lrt) <- gsub(" ", "_", names(mixed_lrt), fixed = TRUE)
p_col <- grep("Pr", names(mixed_lrt), value = TRUE)[1]
mixed_lrt$signif <- sig(mixed_lrt[[p_col]])
write.csv(
  mixed_lrt,
  file.path("{output_dir_r}", "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv"),
  row.names = FALSE
)

# Additive vs interaction comparison
mixed_comparison <- as.data.frame(anova(q1_nn_mixed_additive, q1_nn_mixed))
p_col <- grep("Pr", names(mixed_comparison), value = TRUE)[1]
mixed_comparison$signif <- sig(mixed_comparison[[p_col]])
mixed_comparison$model <- rownames(mixed_comparison)
mixed_comparison <- mixed_comparison[, c("model", setdiff(names(mixed_comparison), "model"))]
write.csv(
  mixed_comparison,
  file.path("{output_dir_r}", "Q1_final_near_normal_mixed_model_comparison_lrt.csv"),
  row.names = FALSE
)

# Fitted values
q1_data$mixed_fitted <- fitted(q1_nn_mixed)
q1_data$mixed_residual <- residuals(q1_nn_mixed)
write.csv(
  q1_data[, c("site_name", "Year", "month", "water_class", "difference_type",
              "Trans_ratio", "Trans_ratio_z_NN", "SPEI_1", "WUE_difference",
              "mixed_fitted", "mixed_residual")],
  file.path("{output_dir_r}", "Q1_final_near_normal_mixed_fitted_residuals.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# MIXED MODEL PREDICTIONS (Panel B)
# ---------------------------------------------------------------------------

tet_sequence <- seq(
  quantile(q1_data$Trans_ratio, 0.05),
  quantile(q1_data$Trans_ratio, 0.95),
  length.out = 100
)

tet_mean_NN <- mean(q1_data$Trans_ratio)
tet_sd_NN <- sd(q1_data$Trans_ratio)

mixed_grid <- expand.grid(
  difference_type = levels(q1_data$difference_type),
  water_class = levels(q1_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(q1_data$site_name)[1],
  month_f = levels(q1_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
mixed_grid$difference_type <- factor(mixed_grid$difference_type, levels = levels(q1_data$difference_type))
mixed_grid$water_class <- factor(mixed_grid$water_class, levels = levels(q1_data$water_class))
mixed_grid$site_name <- factor(mixed_grid$site_name, levels = levels(q1_data$site_name))
mixed_grid$month_f <- factor(mixed_grid$month_f, levels = levels(q1_data$month_f))
mixed_grid$Trans_ratio_z_NN <- (mixed_grid$Trans_ratio - tet_mean_NN) / tet_sd_NN
mixed_grid$difference_label <- factor(diff_labels[as.character(mixed_grid$difference_type)], levels = diff_label_levels)
mixed_grid$ecosystem_label <- as.character(mixed_grid$water_class)

# Fixed-effect predictions
mixed_grid$prediction <- predict(q1_nn_mixed, newdata = mixed_grid, re.form = NA)

# SE using fixed-effect covariance matrix
fixed_terms <- stats::delete.response(stats::terms(lme4::nobars(formula(q1_nn_mixed))))
fixed_model_matrix <- stats::model.matrix(fixed_terms, mixed_grid)
fixed_vcov <- as.matrix(stats::vcov(q1_nn_mixed))
fixed_model_matrix <- fixed_model_matrix[, colnames(fixed_vcov), drop = FALSE]
mixed_grid$se <- sqrt(diag(fixed_model_matrix %*% fixed_vcov %*% t(fixed_model_matrix)))
mixed_grid$lower <- mixed_grid$prediction - 1.96 * mixed_grid$se
mixed_grid$upper <- mixed_grid$prediction + 1.96 * mixed_grid$se

write.csv(
  mixed_grid,
  file.path("{output_dir_r}", "Q1_final_near_normal_mixed_predictions_TET.csv"),
  row.names = FALSE
)

cat(" Mixed model complete!\n")

# =============================================================================
# PANEL C: SMOOTH GAM (Main Panel C)
# =============================================================================

cat("\nPANEL C: Fitting smooth GAM...")

q1_smooth_gam <- mgcv::gam(
  WUE_difference ~ difference_type * water_class +
    s(Trans_ratio, by = diff_ecosystem, k = 6) +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = q1_data,
  method = "ML",
  select = TRUE
)

# Save smooth-term table for edf and p-values
smooth_terms <- as.data.frame(summary(q1_smooth_gam)$s.table)
smooth_terms$term <- rownames(smooth_terms)
rownames(smooth_terms) <- NULL
smooth_terms <- smooth_terms[, c("term", setdiff(names(smooth_terms), "term"))]

# Add significance symbols for smooth-term p-values
p_col_smooth <- grep("p-value", names(smooth_terms), value = TRUE)[1]
if (!is.na(p_col_smooth)) {{
  smooth_terms$signif <- sig(smooth_terms[[p_col_smooth]])
}}

write.csv(
  smooth_terms,
  file.path("{output_dir_r}", "Q1_final_smooth_gam_smooth_terms.csv"),
  row.names = FALSE
)

# SMOOTH GAM PREDICTIONS (Panel C)
smooth_grid <- expand.grid(
  difference_type = levels(q1_data$difference_type),
  water_class = levels(q1_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(q1_data$site_name)[1],
  month_f = levels(q1_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
smooth_grid$difference_type <- factor(smooth_grid$difference_type, levels = levels(q1_data$difference_type))
smooth_grid$water_class <- factor(smooth_grid$water_class, levels = levels(q1_data$water_class))
smooth_grid$site_name <- factor(smooth_grid$site_name, levels = levels(q1_data$site_name))
smooth_grid$month_f <- factor(smooth_grid$month_f, levels = levels(q1_data$month_f))
smooth_grid$difference_label <- factor(diff_labels[as.character(smooth_grid$difference_type)], levels = diff_label_levels)
smooth_grid$ecosystem_label <- as.character(smooth_grid$water_class)
smooth_grid$diff_ecosystem <- interaction(smooth_grid$difference_type, smooth_grid$water_class, sep = "__", drop = TRUE)
smooth_grid$diff_ecosystem <- factor(smooth_grid$diff_ecosystem, levels = levels(q1_data$diff_ecosystem))

smooth_pred <- predict(q1_smooth_gam, newdata = smooth_grid, type = "link", se.fit = TRUE, 
                        exclude = c("s(site_name)", "s(month_f)"))
smooth_grid$prediction <- as.numeric(smooth_pred$fit)
smooth_grid$se <- as.numeric(smooth_pred$se.fit)
smooth_grid$lower <- smooth_grid$prediction - 1.96 * smooth_grid$se
smooth_grid$upper <- smooth_grid$prediction + 1.96 * smooth_grid$se

write.csv(
  smooth_grid,
  file.path("{output_dir_r}", "Q1_final_smooth_gam_predictions_TET.csv"),
  row.names = FALSE
)

cat(" Smooth GAM complete!\n")

# =============================================================================
# LINEAR GAM (For Panel D model comparison only - NOT plotted as main panel)
# =============================================================================

cat("\nPANEL D: Fitting linear GAM for model comparison...")

q1_linear_gam <- mgcv::gam(
  WUE_difference ~ difference_type * water_class * Trans_ratio +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = q1_data,
  method = "ML"
)

# LINEAR GAM PREDICTIONS (Diagnostic only - not used in main figure)
linear_grid <- expand.grid(
  difference_type = levels(q1_data$difference_type),
  water_class = levels(q1_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(q1_data$site_name)[1],
  month_f = levels(q1_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
linear_grid$difference_type <- factor(linear_grid$difference_type, levels = levels(q1_data$difference_type))
linear_grid$water_class <- factor(linear_grid$water_class, levels = levels(q1_data$water_class))
linear_grid$site_name <- factor(linear_grid$site_name, levels = levels(q1_data$site_name))
linear_grid$month_f <- factor(linear_grid$month_f, levels = levels(q1_data$month_f))
linear_grid$difference_label <- factor(diff_labels[as.character(linear_grid$difference_type)], levels = diff_label_levels)
linear_grid$ecosystem_label <- as.character(linear_grid$water_class)

linear_pred <- predict(q1_linear_gam, newdata = linear_grid, type = "link", se.fit = TRUE, 
                       exclude = c("s(site_name)", "s(month_f)"))
linear_grid$prediction <- as.numeric(linear_pred$fit)
linear_grid$se <- as.numeric(linear_pred$se.fit)
linear_grid$lower <- linear_grid$prediction - 1.96 * linear_grid$se
linear_grid$upper <- linear_grid$prediction + 1.96 * linear_grid$se

write.csv(
  linear_grid,
  file.path("{output_dir_r}", "Q1_final_linear_gam_predictions_TET.csv"),
  row.names = FALSE
)

# =============================================================================
# PANEL D: MODEL COMPARISON (Linear GAM vs Smooth GAM)
# =============================================================================

gam_comparison <- data.frame(
  model = c("q1_linear_gam", "q1_smooth_gam"),
  AIC = c(AIC(q1_linear_gam), AIC(q1_smooth_gam)),
  BIC = c(BIC(q1_linear_gam), BIC(q1_smooth_gam)),
  deviance_explained = c(summary(q1_linear_gam)$dev.expl, summary(q1_smooth_gam)$dev.expl),
  adjusted_r_squared = c(summary(q1_linear_gam)$r.sq, summary(q1_smooth_gam)$r.sq),
  scale = c(summary(q1_linear_gam)$scale, summary(q1_smooth_gam)$scale)
)
gam_comparison$delta_AIC <- gam_comparison$AIC - min(gam_comparison$AIC)
write.csv(
  gam_comparison,
  file.path("{output_dir_r}", "Q1_final_gam_model_comparison.csv"),
  row.names = FALSE
)

cat(" Model comparison complete!\n")

# =============================================================================
# MODEL SUMMARY
# =============================================================================

sink(file.path("{output_dir_r}", "Q1_final_model_summary.txt"))
cat("Q1 FINAL - Near-Normal SPEI-1 Workflow with Mixed Model Panel B\n")
cat("===============================================================\n")
cat("Original input:", "{original_input_path}", "\n")
cat("R model input:", "{temp_model_base_file_r}", "\n")
cat("Output:", "{output_dir_r}", "\n\n")
cat("Full data rows (after cleaning):", nrow(diff_data), "\n")
cat("Near-normal SPEI-1 rows (q1_data):", nrow(q1_data), "\n")
cat("Sites in q1_data:", length(unique(q1_data$site_name)), "\n")
cat("Months in q1_data:", length(unique(q1_data$month_f)), "\n\n")
cat("Ecosystem classes:", paste(levels(q1_data$water_class), collapse = ", "), "\n\n")
cat("Difference types:", paste(levels(q1_data$difference_type), collapse = ", "), "\n\n")
cat("----------------------------------------------------------------\n")
cat("PANEL A: Summary Statistics\n")
cat("  Near-normal mean WUE metric differences by ecosystem\n")
cat("  Saved to: Q1_final_near_normal_summary_by_difference_type.csv\n\n")
cat("----------------------------------------------------------------\n")
cat("PANEL B: Near-Normal Mixed Model\n")
cat("Formula:\n")
print(formula(q1_nn_mixed))
cat("\nMixed Model Summary:\n")
print(summary(q1_nn_mixed))
cat("\nMixed Drop1 LRT:\n")
print(mixed_lrt)
cat("\nMixed Additive vs Interaction Comparison:\n")
print(mixed_comparison)
cat("\n----------------------------------------------------------------\n")
cat("PANEL C: Smooth GAM\n")
cat("Formula:\n")
print(formula(q1_smooth_gam))
cat("\nSmooth GAM Summary:\n")
print(summary(q1_smooth_gam))
cat("\n----------------------------------------------------------------\n")
cat("PANEL D: Linear GAM vs Smooth GAM Comparison\n")
cat("Linear GAM formula (for comparison only):\n")
print(formula(q1_linear_gam))
cat("\nLinear GAM Summary:\n")
print(summary(q1_linear_gam))
cat("\nModel Comparison:\n")
print(gam_comparison)
cat("\nSmooth GAM has lower AIC than linear GAM:",
    ifelse(gam_comparison$AIC[2] < gam_comparison$AIC[1], "YES", "NO"), "\n")
sink()

cat("\nR analysis complete! All outputs saved to:", "{output_dir_r}", "\n")
'''

# Save R script in temp folder
r_script_file = os.path.join(temp_dir, "run_final_models.R")
with open(r_script_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print(f"  R script saved to: {r_script_file}")

# Verify temp files exist
print("\nVerifying temp files exist before running R:")
print(f"  temp_model_base_for_R.csv exists: {os.path.exists(temp_model_base_file)}")
print(f"  run_final_models.R exists: {os.path.exists(r_script_file)}")

# ============================================================================
# STEP 4: CLEAN OLD OUTPUTS AND FIGURES BEFORE RUNNING R
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Cleaning old outputs and figures")
print("="*60)

# Clean R outputs (main outputs)
expected_r_outputs = [
    # Panel A
    "Q1_final_near_normal_summary_by_difference_type.csv",
    # Panel B
    "Q1_final_near_normal_mixed_difference_TET_model.rds",
    "Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv",
    "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv",
    "Q1_final_near_normal_mixed_model_comparison_lrt.csv",
    "Q1_final_near_normal_mixed_fitted_residuals.csv",
    "Q1_final_near_normal_mixed_predictions_TET.csv",
    # Panel C
    "Q1_final_smooth_gam_predictions_TET.csv",
    "Q1_final_smooth_gam_smooth_terms.csv",
    # Panel D
    "Q1_final_linear_gam_predictions_TET.csv",
    "Q1_final_gam_model_comparison.csv",
    # Summary
    "Q1_final_model_summary.txt",
]

for f in expected_r_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

# Clean figure files
expected_figures = [
    "Q1_final_plot_07_manuscript_multipanel.png",
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
print(f"  run_final_models.R exists: {os.path.exists(r_script_file)}")

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
# STEP 7: PYTHON READS R OUTPUT AND CREATES FINAL FIGURE
# ============================================================================

print("\n" + "="*60)
print("STEP 6: Creating Final Python Figure from R Outputs")
print("="*60)

# Load R outputs
print("\nLoading R output CSVs...")
summary_df = pd.read_csv(os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv"))
mixed_pred = pd.read_csv(os.path.join(output_dir, "Q1_final_near_normal_mixed_predictions_TET.csv"))
smooth_pred = pd.read_csv(os.path.join(output_dir, "Q1_final_smooth_gam_predictions_TET.csv"))
smooth_terms = pd.read_csv(os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv"))
gam_comp = pd.read_csv(os.path.join(output_dir, "Q1_final_gam_model_comparison.csv"))

# Define ecosystem colors
ecosystem_colors = {"Upland": "#800080", "Freshwater": "#0000FF", "Saline": "#FFA500"}
ecosystem_classes = ["Upland", "Freshwater", "Saline"]
diff_label_levels = ["WUE_ET - WUE_T", "WUE_ET - WUE_E", "WUE_E - WUE_T"]

print(f"  Summary: {len(summary_df)} rows")
print(f"  Mixed predictions: {len(mixed_pred)} rows")
print(f"  Smooth GAM predictions: {len(smooth_pred)} rows")
print(f"  Smooth GAM terms: {len(smooth_terms)} rows")

# Apply theme
theme_wue_manuscript()

# ============================================================================
# FINAL MANUSCRIPT MULTI-PANEL FIGURE
# Panel A: Mean WUE Metric Differences under Near-Normal SPEI-1
# Panel B: Near-Normal Mixed Model T:ET Response
# Panel C: Near-Normal Smooth GAM T:ET Response
# Panel D: Linear vs Smooth GAM Model Comparison
# ============================================================================

print("\nCreating Final Manuscript Multi-Panel Figure...")

fig = plt.figure(figsize=(24, 25))
outer = gridspec.GridSpec(4, 1, height_ratios=[0.9, 1.25, 1.25, 0.8],
                          hspace=0.35, top=0.95, bottom=0.05)

# Row A: Summary (Panel A)
gs_a = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0], wspace=0.25)
axes_a = [fig.add_subplot(gs_a[0, i]) for i in range(3)]

# Row B: Near-normal mixed model predictions (Panel B)
gs_b = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1], wspace=0.25)
axes_b = [fig.add_subplot(gs_b[0, i]) for i in range(3)]

# Row C: Smooth GAM predictions (Panel C)
gs_c = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[2], wspace=0.25)
axes_c = [fig.add_subplot(gs_c[0, i]) for i in range(3)]

# Row D: Model comparison (Panel D)
gs_d = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[3])
ax_d = fig.add_subplot(gs_d[0, 0])

# ===== Panel A: Mean WUE Metric Differences =====
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
    if idx == 0: 
        ax.set_ylabel('Mean WUE metric difference', fontsize=18)
    ax.tick_params(axis='x', rotation=25, labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.grid(True, alpha=0.3)
axes_a[0].text(-0.15, 1.05, 'A', transform=axes_a[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

# ===== Panel B: Near-Normal Mixed Model =====
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes_b[idx]
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
    if idx == 0: 
        ax.set_ylabel('Predicted WUE metric difference', fontsize=18)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    if idx == 2: 
        ax.legend(loc='best', fontsize=16, frameon=True)
axes_b[0].text(-0.15, 1.05, 'B', transform=axes_b[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

# ===== Panel C: Smooth GAM =====
for idx, diff_type in enumerate(diff_label_levels):
    ax = axes_c[idx]
    subset = smooth_pred[smooth_pred['difference_label'] == diff_type]
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
    if idx == 0: 
        ax.set_ylabel('Predicted WUE metric difference', fontsize=18)
    ax.tick_params(labelsize=16)
    ax.grid(True, alpha=0.3)
    if idx == 2: 
        ax.legend(loc='best', fontsize=16, frameon=True)
axes_c[0].text(-0.15, 1.05, 'C', transform=axes_c[0].transAxes, 
               fontsize=38, fontweight='bold', va='top')

# ===== Panel D: Model Comparison =====
models = ['Linear GAM', 'Smooth GAM']
deviance = [gam_comp.loc[0, 'deviance_explained'] * 100, 
            gam_comp.loc[1, 'deviance_explained'] * 100]
bars = ax_d.bar(models, deviance, color=['#8DA0CB', '#66C2A5'], 
                edgecolor='#4D4D4D', linewidth=2)
for bar, val in zip(bars, deviance):
    ax_d.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f'{val:.1f}%', ha='center', va='bottom', fontsize=20, fontweight='bold')
ax_d.set_ylabel('Deviance explained (%)', fontsize=22)
ax_d.set_title('Linear vs Smooth GAM Model Comparison', fontsize=25, fontweight='bold')
ax_d.set_ylim(0, 95)
ax_d.tick_params(labelsize=18)
ax_d.grid(True, alpha=0.3, axis='y')

# Calculate AIC improvement conditionally
aic_linear = gam_comp.loc[gam_comp['model'] == 'q1_linear_gam', 'AIC'].values[0]
aic_smooth = gam_comp.loc[gam_comp['model'] == 'q1_smooth_gam', 'AIC'].values[0]
aic_improvement = aic_linear - aic_smooth

if aic_improvement > 0:
    aic_text = f"Smooth GAM improves AIC by {aic_improvement:.1f}"
else:
    aic_text = f"Smooth GAM does not improve AIC (ΔAIC = {aic_improvement:.1f})"

ax_d.text(0.5, -0.12, aic_text,
          transform=ax_d.transAxes, ha='center', va='center', fontsize=18)

ax_d.text(-0.08, 1.05, 'D', transform=ax_d.transAxes, 
          fontsize=38, fontweight='bold', va='top')

fig.suptitle('Q1 WUE Metric Difference Models (Near-Normal SPEI-1)', 
             fontsize=42, fontweight='bold', y=0.98)
save_plot(fig, "Q1_final_plot_07_manuscript_multipanel.png", width=24, height=25)

print("  Main figure saved: Q1_final_plot_07_manuscript_multipanel.png")

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

readme_content = """# Q1 Final Results - Near-Normal SPEI-1 Workflow (Mixed Model Panel B)

This is the FINAL Q1 analysis for the WUE manuscript. 

## Key Decision
**All four main panels use near-normal hydroclimatic conditions only:**
- SPEI_1 >= -1 & SPEI_1 <= 1

## Panel Structure
- **Panel A**: Near-normal mean WUE metric differences by ecosystem
- **Panel B**: Near-normal mixed model T:ET response by ecosystem
- **Panel C**: Near-normal smooth GAM T:ET response by ecosystem
- **Panel D**: Linear vs Smooth GAM model comparison

## Directory Structure
- `outputs/` - CSV files, model summary
- `figures/` - PNG figure files
- `temp/` - Temporary R handoff files (cleaned at start)

## Workflow
- Python: Data cleaning, filtering, orchestrating R
- R: Reshaping to long format, near-normal filtering, mixed model, GAM models, predictions
- Python: Reading R outputs, creating final figure, generating README

## MAIN Q1 Results (near-normal SPEI-1 only)

### Panel A: Summary Statistics
- Mean WUE metric differences by ecosystem under near-normal conditions

### Panel B: Near-Normal Mixed Model (Main Panel B)
`WUE_difference ~ difference_type * water_class * Trans_ratio_z_NN + (1 | site_name) + (1 | month_f)`

- Uses standardized Trans_ratio_z_NN internally
- Plotted against raw Trans_ratio on x-axis
- Fixed-effect predictions with site and month random effects excluded

### Panel C: Smooth GAM (Main Panel C)
`WUE_difference ~ difference_type * water_class + s(Trans_ratio, by = diff_ecosystem, k = 6) + s(site_name, bs = "re") + s(month_f, bs = "re")`

- Uses raw Trans_ratio
- Near-normal SPEI-1 subset only

### Panel D: Model Comparison
- Linear GAM vs Smooth GAM
- Deviance explained and AIC comparison
- Linear GAM retained for comparison only (not shown as response-curve panel)

## Key Tables (in outputs folder)
### Panel A
- `Q1_final_near_normal_summary_by_difference_type.csv`

### Panel B (Mixed Model)
- `Q1_final_near_normal_mixed_difference_TET_model.rds`
- `Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv`
- `Q1_final_near_normal_mixed_likelihood_ratio_tests.csv`
- `Q1_final_near_normal_mixed_model_comparison_lrt.csv`
- `Q1_final_near_normal_mixed_fitted_residuals.csv`
- `Q1_final_near_normal_mixed_predictions_TET.csv`

### Panel C (Smooth GAM)
- `Q1_final_smooth_gam_predictions_TET.csv`
- `Q1_final_smooth_gam_smooth_terms.csv`

### Panel D (Model Comparison)
- `Q1_final_linear_gam_predictions_TET.csv` (diagnostic only)
- `Q1_final_gam_model_comparison.csv`

### Summary
- `Q1_final_model_summary.txt`

## Figure
- `Q1_final_plot_07_manuscript_multipanel.png`
- 4-panel figure with panels A-D as described above

## Notes
- Near-normal SPEI-1 filter applied before all analyses
- Mixed model uses standardized Trans_ratio_z_NN internally
- GAMs use raw Trans_ratio
- Linear GAM retained only for Panel D model comparison
- Panel B and C both use raw Trans_ratio on x-axis for consistency
"""

with open(os.path.join(output_dir, "README.md"), 'w') as f:
    f.write(readme_content)

print("  Created: README.md")

# ============================================================================
# FINAL VALIDATION CHECKS
# ============================================================================

print("\n" + "="*60)
print("FINAL VALIDATION CHECKS")
print("="*60)

# Check original data
print(f"\n1. Original model_base rows: {len(model_base)}")

# Load R summary to get diff_data and q1_data counts
summary_file = os.path.join(output_dir, "Q1_final_model_summary.txt")
if os.path.exists(summary_file):
    with open(summary_file, 'r') as f:
        summary_text = f.read()
    
    # Extract counts from summary
    import re
    diff_data_match = re.search(r"Full data rows \(after cleaning\):\s*(\d+)", summary_text)
    q1_data_match = re.search(r"Near-normal SPEI-1 rows \(q1_data\):\s*(\d+)", summary_text)
    sites_match = re.search(r"Sites in q1_data:\s*(\d+)", summary_text)
    
    if diff_data_match:
        print(f"2. diff_data rows: {diff_data_match.group(1)}")
    if q1_data_match:
        print(f"3. q1_data rows (near-normal): {q1_data_match.group(1)}")
    if sites_match:
        print(f"4. Sites in q1_data: {sites_match.group(1)}")

# Check main final CSVs exist
print("\n5. Checking final CSV outputs:")
main_csvs = [
    "Q1_final_near_normal_summary_by_difference_type.csv",
    "Q1_final_near_normal_mixed_predictions_TET.csv",
    "Q1_final_smooth_gam_predictions_TET.csv",
    "Q1_final_smooth_gam_smooth_terms.csv",
    "Q1_final_gam_model_comparison.csv",
    "Q1_final_linear_gam_predictions_TET.csv",
    "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv",
    "Q1_final_near_normal_mixed_model_comparison_lrt.csv",
]
for csv_file in main_csvs:
    exists = os.path.exists(os.path.join(output_dir, csv_file))
    print(f"   {csv_file}: {'✓' if exists else '✗'}")

# Check main figure exists
print(f"\n6. Q1_final_plot_07_manuscript_multipanel.png: {'✓' if os.path.exists(os.path.join(figure_dir, 'Q1_final_plot_07_manuscript_multipanel.png')) else '✗'}")

# Print model comparison
print("\n7. Model Comparison:")
print(gam_comp.to_string(index=False))

# Check if smooth GAM has lower AIC
smooth_lower = gam_comp.loc[1, 'AIC'] < gam_comp.loc[0, 'AIC']
print(f"\n8. Smooth GAM has lower AIC than linear GAM: {'✓' if smooth_lower else '✗'}")

# Load and print mixed LRT if exists
lrt_file = os.path.join(output_dir, "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv")
if os.path.exists(lrt_file):
    print("\n9. Mixed Model Drop1 LRT:")
    lrt_df = pd.read_csv(lrt_file)
    print(lrt_df.to_string(index=False))

# Load and print mixed comparison if exists
comp_file = os.path.join(output_dir, "Q1_final_near_normal_mixed_model_comparison_lrt.csv")
if os.path.exists(comp_file):
    print("\n10. Mixed Additive vs Interaction Comparison:")
    comp_df = pd.read_csv(comp_file)
    print(comp_df.to_string(index=False))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("COMPLETE! Final Q1 Workflow Finished")
print("="*60)

print(f"\nOutputs saved to:")
print(f"  Outputs (CSV/README): {output_dir}")
print(f"  Figures (PNG):        {figure_dir}")
print(f"  Temp files:           {temp_dir}")

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
print("✓ Final Q1 workflow complete")
print("  - Panel order: A=Summary, B=Near-normal mixed model, C=Smooth GAM, D=Linear vs Smooth GAM comparison")
print("  - Mixed model uses standardized Trans_ratio_z_NN internally")
print("  - Linear and smooth GAMs use raw Trans_ratio")
print("  - All main panels use near-normal SPEI-1 (>= -1 & <= 1)")
print("  - Linear GAM retained only for Panel D model comparison")
print("  - All outputs in Q1_updated_results (separate from old outputs)")
print("="*60)