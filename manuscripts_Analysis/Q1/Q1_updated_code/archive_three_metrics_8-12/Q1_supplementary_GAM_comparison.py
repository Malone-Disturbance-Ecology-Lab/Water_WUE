# =============================================================================
# Q1 STANDALONE GAM REDUCED-MODEL COMPARISON
#
# Purpose:
#   Add reduced nonlinear GAM comparisons to replace the mixed-model interaction
#   evidence with GAM-only evidence.
#
# What this script saves:
#   1. Model A vs Model B comparison
#   2. Model B vs Model C comparison  <-- most important
#   3. Model A vs Model C comparison
#   4. AIC/BIC/deviance explained/adjusted R2 for Models A, B, and C
#   5. Reduced-model predictions for possible Supplementary figure
#
# Important:
#   - Your original workflow already has Model C full-GAM results.
#   - This script refits Model C only because the original workflow did not save
#     the full GAM object as an RDS.
#   - This script does NOT overwrite original Q1 outputs.
#
# New outputs:
#   M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs\GAM_comparison
# =============================================================================

import os
import subprocess
import numpy as np
import pandas as pd

# =============================================================================
# PATHS
# =============================================================================

monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

original_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs"

gam_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs\GAM_comparison"
os.makedirs(gam_output_dir, exist_ok=True)

temp_model_base_file = os.path.join(gam_output_dir, "temp_model_base_for_GAM_comparison.csv")
r_script_file = os.path.join(gam_output_dir, "run_GAM_reduced_BC_comparison.R")

print("=" * 90)
print("Q1 standalone GAM reduced-model comparison")
print("=" * 90)
print(f"Output folder: {gam_output_dir}")

# =============================================================================
# STEP 1: LOAD AND CLEAN DATA SAME WAY AS ORIGINAL Q1 WORKFLOW
# =============================================================================

print("\nLoading monthly data...")
monthly = pd.read_csv(monthly_data_path)
print(f"Loaded rows: {len(monthly)}")

for col in monthly.select_dtypes(include="object").columns:
    monthly[col] = monthly[col].where(
        monthly[col].isna(),
        monthly[col].astype(str).str.strip()
    )
    monthly[col] = monthly[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

required_cols = [
    "site_name", "Year", "month", "water_class", "Trans_ratio",
    "WUE", "WUE_tra", "WUE_eva", "SPEI_1"
]

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

ecosystem_classes = ["Upland", "Freshwater", "Saline"]

model_base = monthly.dropna(subset=required_cols).copy()

model_base = model_base[
    (model_base["site_name"] != "") &
    (model_base["water_class"].isin(ecosystem_classes)) &
    np.isfinite(model_base["Trans_ratio"]) &
    (model_base["Trans_ratio"] >= 0) &
    (model_base["Trans_ratio"] <= 1) &
    np.isfinite(model_base["WUE"]) &
    np.isfinite(model_base["WUE_tra"]) &
    np.isfinite(model_base["WUE_eva"]) &
    np.isfinite(model_base["SPEI_1"])
].copy()

print(f"Rows after cleaning: {len(model_base)}")

dup_count = model_base.duplicated(subset=["site_name", "Year", "month"]).sum()
print(f"Duplicate site-Year-month rows: {dup_count}")

model_base.to_csv(temp_model_base_file, index=False)
print(f"Saved temporary clean input for R: {temp_model_base_file}")

# =============================================================================
# STEP 2: CREATE R SCRIPT
# =============================================================================

gam_output_dir_r = gam_output_dir.replace("\\", "/")
original_output_dir_r = original_output_dir.replace("\\", "/")
temp_model_base_file_r = temp_model_base_file.replace("\\", "/")

r_script_content = r'''
# =============================================================================
# Q1 STANDALONE GAM REDUCED-MODEL COMPARISON
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required.")
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

output_dir <- "__GAM_OUTPUT_DIR__"
original_output_dir <- "__ORIGINAL_OUTPUT_DIR__"

model_base <- read.csv("__TEMP_MODEL_BASE_FILE__", stringsAsFactors = FALSE)

ecosystem_classes <- c("Upland", "Freshwater", "Saline")

diff_labels <- c(
  WUE_ET_minus_WUE_T = "WUE_ET - WUE_T",
  WUE_ET_minus_WUE_E = "WUE_ET - WUE_E",
  WUE_E_minus_WUE_T = "WUE_E - WUE_T"
)

difference_cols <- c(
  "WUE_ET_minus_WUE_T",
  "WUE_ET_minus_WUE_E",
  "WUE_E_minus_WUE_T"
)

# =============================================================================
# CREATE SIGNED WUE DIFFERENCES
# =============================================================================

model_base$WUE_ET_minus_WUE_T <- model_base$WUE - model_base$WUE_tra
model_base$WUE_ET_minus_WUE_E <- model_base$WUE - model_base$WUE_eva
model_base$WUE_E_minus_WUE_T <- model_base$WUE_eva - model_base$WUE_tra

# =============================================================================
# RESHAPE TO LONG FORMAT
# =============================================================================

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
  complete.cases(diff_data[, c("site_name", "Year", "month",
                               "Trans_ratio", "WUE_difference", "SPEI_1")]),
]

diff_data <- diff_data[is.finite(diff_data$WUE_difference), ]

diff_data$site_name <- factor(diff_data$site_name)
diff_data$month_f <- factor(diff_data$month, levels = 1:12)
diff_data$water_class <- factor(diff_data$water_class, levels = ecosystem_classes)
diff_data$difference_type <- factor(diff_data$difference_type, levels = difference_cols)
diff_data$difference_label <- factor(
  diff_labels[as.character(diff_data$difference_type)],
  levels = unname(diff_labels[difference_cols])
)

# =============================================================================
# NEAR-NORMAL SPEI-1 FILTER
# =============================================================================

q1_data <- diff_data[
  is.finite(diff_data$SPEI_1) &
    diff_data$SPEI_1 >= -1 &
    diff_data$SPEI_1 <= 1,
]

q1_data <- droplevels(q1_data)

q1_data$diff_ecosystem <- interaction(
  q1_data$difference_type,
  q1_data$water_class,
  sep = "__",
  drop = TRUE
)

cat("\nQ1 near-normal data:")
cat("\n  Rows:", nrow(q1_data))
cat("\n  Site-month observations:", nrow(q1_data) / 3)
cat("\n  Sites:", length(unique(q1_data$site_name)))
cat("\n  Months:", length(unique(q1_data$month_f)))
cat("\n")

# =============================================================================
# MODEL A: ADDITIVE / COMMON-SMOOTH NONLINEAR GAM
# =============================================================================
# This is the simple additive GAM:
# WUE comparison + ecosystem class + one common nonlinear T:ET smooth.

q1_gam_A_additive_common_smooth <- mgcv::gam(
  WUE_difference ~ difference_type + water_class +
    s(Trans_ratio, k = 6) +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = q1_data,
  method = "ML",
  select = TRUE
)

# =============================================================================
# MODEL B: MEAN-INTERACTION / COMMON-SMOOTH NONLINEAR GAM
# =============================================================================
# This model allows WUE comparison × ecosystem differences in mean WUE
# difference, but still uses one common nonlinear T:ET smooth for all groups.

q1_gam_B_mean_interaction_common_smooth <- mgcv::gam(
  WUE_difference ~ difference_type * water_class +
    s(Trans_ratio, k = 6) +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = q1_data,
  method = "ML",
  select = TRUE
)

# =============================================================================
# MODEL C: FULL SEPARATE-SMOOTH NONLINEAR GAM
# =============================================================================
# This is the same model structure as your original full nonlinear GAM.
# It is refit here only because model comparison requires fitted model objects.
# We do NOT save duplicate full-GAM smooth-term or prediction outputs here.

q1_gam_C_full_separate_smooths <- mgcv::gam(
  WUE_difference ~ difference_type * water_class +
    s(Trans_ratio, by = diff_ecosystem, k = 6) +
    s(site_name, bs = "re") +
    s(month_f, bs = "re"),
  data = q1_data,
  method = "ML",
  select = TRUE
)

# =============================================================================
# MODEL PERFORMANCE TABLE
# =============================================================================

get_stats <- function(model, model_name, description) {
  sm <- summary(model)
  data.frame(
    model = model_name,
    description = description,
    AIC = AIC(model),
    BIC = BIC(model),
    deviance_explained = sm$dev.expl,
    adjusted_r_squared = sm$r.sq,
    scale = sm$scale,
    n = nobs(model)
  )
}

model_performance <- rbind(
  get_stats(
    q1_gam_A_additive_common_smooth,
    "A_additive_common_smooth",
    "difference_type + water_class + one common nonlinear TET smooth"
  ),
  get_stats(
    q1_gam_B_mean_interaction_common_smooth,
    "B_mean_interaction_common_smooth",
    "difference_type * water_class + one common nonlinear TET smooth"
  ),
  get_stats(
    q1_gam_C_full_separate_smooths,
    "C_full_separate_smooths",
    "difference_type * water_class + separate nonlinear TET smooths by WUE comparison and ecosystem"
  )
)

model_performance$delta_AIC <- model_performance$AIC - min(model_performance$AIC)

write.csv(
  model_performance,
  file.path(output_dir, "Q1_GAM_A_B_C_model_performance.csv"),
  row.names = FALSE
)

# =============================================================================
# KEY MODEL COMPARISONS
# =============================================================================

extract_comparison <- function(model_small, model_big, comparison_label, question) {
  raw <- as.data.frame(anova(model_small, model_big, test = "Chisq"))
  
  p_col <- grep("Pr", names(raw), value = TRUE)[1]
  p_val <- NA
  if (!is.na(p_col)) {
    p_val <- raw[[p_col]][nrow(raw)]
  }
  
  data.frame(
    comparison = comparison_label,
    question = question,
    smaller_model_AIC = AIC(model_small),
    larger_model_AIC = AIC(model_big),
    delta_AIC_smaller_minus_larger = AIC(model_small) - AIC(model_big),
    smaller_deviance_explained = summary(model_small)$dev.expl,
    larger_deviance_explained = summary(model_big)$dev.expl,
    deviance_explained_gain = summary(model_big)$dev.expl - summary(model_small)$dev.expl,
    p_value = p_val,
    signif = sig(p_val)
  )
}

key_evidence <- rbind(
  extract_comparison(
    q1_gam_A_additive_common_smooth,
    q1_gam_B_mean_interaction_common_smooth,
    "A_vs_B",
    "Does adding WUE comparison × ecosystem mean interaction improve fit?"
  ),
  extract_comparison(
    q1_gam_B_mean_interaction_common_smooth,
    q1_gam_C_full_separate_smooths,
    "B_vs_C",
    "Do separate nonlinear TET smooths by WUE comparison × ecosystem improve fit?"
  ),
  extract_comparison(
    q1_gam_A_additive_common_smooth,
    q1_gam_C_full_separate_smooths,
    "A_vs_C",
    "Does the full separate-smooth GAM improve fit compared with the additive GAM?"
  )
)

write.csv(
  key_evidence,
  file.path(output_dir, "Q1_GAM_key_evidence_for_replacing_mixed_model.csv"),
  row.names = FALSE
)

# Save raw anova outputs too
save_raw_anova <- function(model_small, model_big, filename) {
  out <- as.data.frame(anova(model_small, model_big, test = "Chisq"))
  out$model_row <- rownames(out)
  rownames(out) <- NULL
  out <- out[, c("model_row", setdiff(names(out), "model_row"))]
  write.csv(out, file.path(output_dir, filename), row.names = FALSE)
}

save_raw_anova(
  q1_gam_A_additive_common_smooth,
  q1_gam_B_mean_interaction_common_smooth,
  "Q1_GAM_raw_LRT_A_vs_B.csv"
)

save_raw_anova(
  q1_gam_B_mean_interaction_common_smooth,
  q1_gam_C_full_separate_smooths,
  "Q1_GAM_raw_LRT_B_vs_C.csv"
)

save_raw_anova(
  q1_gam_A_additive_common_smooth,
  q1_gam_C_full_separate_smooths,
  "Q1_GAM_raw_LRT_A_vs_C.csv"
)

# =============================================================================
# CHECK THAT REFIT FULL GAM MATCHES ORIGINAL FULL GAM AIC
# =============================================================================

original_gam_comparison_file <- file.path(original_output_dir, "Q1_final_gam_model_comparison.csv")

if (file.exists(original_gam_comparison_file)) {
  original_comp <- read.csv(original_gam_comparison_file, stringsAsFactors = FALSE)
  original_smooth <- original_comp[original_comp$model == "q1_smooth_gam", ]
  
  if (nrow(original_smooth) == 1) {
    full_refit_check <- data.frame(
      original_full_GAM_AIC = original_smooth$AIC,
      refit_full_GAM_AIC = AIC(q1_gam_C_full_separate_smooths),
      AIC_difference_refit_minus_original = AIC(q1_gam_C_full_separate_smooths) - original_smooth$AIC,
      original_full_GAM_deviance_explained = original_smooth$deviance_explained,
      refit_full_GAM_deviance_explained = summary(q1_gam_C_full_separate_smooths)$dev.expl,
      deviance_explained_difference_refit_minus_original =
        summary(q1_gam_C_full_separate_smooths)$dev.expl - original_smooth$deviance_explained
    )
    
    write.csv(
      full_refit_check,
      file.path(output_dir, "Q1_GAM_full_refit_check_against_original.csv"),
      row.names = FALSE
    )
  }
}

# =============================================================================
# REDUCED-MODEL PREDICTIONS FOR POSSIBLE SUPPLEMENTARY FIGURE
# =============================================================================
# Original full-GAM predictions already exist:
#   Q1_final_smooth_gam_predictions_TET.csv
#
# Here we save only Model A and Model B predictions.

tet_sequence <- seq(
  quantile(q1_data$Trans_ratio, 0.05),
  quantile(q1_data$Trans_ratio, 0.95),
  length.out = 100
)

pred_grid <- expand.grid(
  difference_type = levels(q1_data$difference_type),
  water_class = levels(q1_data$water_class),
  Trans_ratio = tet_sequence,
  site_name = levels(q1_data$site_name)[1],
  month_f = levels(q1_data$month_f)[1],
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pred_grid$difference_type <- factor(pred_grid$difference_type, levels = levels(q1_data$difference_type))
pred_grid$water_class <- factor(pred_grid$water_class, levels = levels(q1_data$water_class))
pred_grid$site_name <- factor(pred_grid$site_name, levels = levels(q1_data$site_name))
pred_grid$month_f <- factor(pred_grid$month_f, levels = levels(q1_data$month_f))
pred_grid$difference_label <- factor(
  diff_labels[as.character(pred_grid$difference_type)],
  levels = unname(diff_labels[difference_cols])
)
pred_grid$ecosystem_label <- as.character(pred_grid$water_class)

predict_reduced_model <- function(model, grid, model_name) {
  pr <- predict(
    model,
    newdata = grid,
    type = "link",
    se.fit = TRUE,
    exclude = c("s(site_name)", "s(month_f)")
  )
  
  out <- grid
  out$model <- model_name
  out$prediction <- as.numeric(pr$fit)
  out$se <- as.numeric(pr$se.fit)
  out$lower <- out$prediction - 1.96 * out$se
  out$upper <- out$prediction + 1.96 * out$se
  out
}

pred_A <- predict_reduced_model(
  q1_gam_A_additive_common_smooth,
  pred_grid,
  "A_additive_common_smooth"
)

pred_B <- predict_reduced_model(
  q1_gam_B_mean_interaction_common_smooth,
  pred_grid,
  "B_mean_interaction_common_smooth"
)

reduced_predictions <- rbind(pred_A, pred_B)

write.csv(
  reduced_predictions,
  file.path(output_dir, "Q1_GAM_A_B_reduced_predictions_TET_for_supplement.csv"),
  row.names = FALSE
)

# =============================================================================
# TEXT SUMMARY
# =============================================================================

sink(file.path(output_dir, "Q1_GAM_A_B_C_comparison_summary.txt"))

cat("Q1 GAM A/B/C MODEL COMPARISON\n")
cat("=============================\n\n")

cat("Purpose:\n")
cat("This script adds GAM-only evidence to replace the mixed-model interaction tests.\n\n")

cat("Important note:\n")
cat("Model C is the same full separate-smooth GAM used in the original Q1 workflow.\n")
cat("It is refit here only because model comparison requires fitted model objects.\n")
cat("The original full-GAM smooth terms and prediction outputs are not recreated here.\n\n")

cat("Data:\n")
cat("Near-normal SPEI-1 rows:", nrow(q1_data), "\n")
cat("Site-month observations:", nrow(q1_data) / 3, "\n")
cat("Sites:", length(unique(q1_data$site_name)), "\n")
cat("Months:", length(unique(q1_data$month_f)), "\n\n")

cat("Model A: additive/common-smooth GAM\n")
print(formula(q1_gam_A_additive_common_smooth))
cat("\n\n")

cat("Model B: mean-interaction/common-smooth GAM\n")
print(formula(q1_gam_B_mean_interaction_common_smooth))
cat("\n\n")

cat("Model C: full separate-smooth GAM\n")
print(formula(q1_gam_C_full_separate_smooths))
cat("\n\n")

cat("Model performance:\n")
print(model_performance)
cat("\n\n")

cat("Key evidence:\n")
print(key_evidence)
cat("\n\n")

cat("Interpretation:\n")
cat("The most important comparison is B_vs_C.\n")
cat("B_vs_C tests whether separate nonlinear T:ET smooths for each WUE comparison × ecosystem class\n")
cat("improve fit compared with one common nonlinear T:ET smooth after allowing WUE comparison × ecosystem\n")
cat("mean differences. This is the GAM replacement for the mixed-model three-way interaction idea.\n\n")

sink()

cat("\nAll new GAM comparison outputs saved to:", output_dir, "\n")
'''.replace("__GAM_OUTPUT_DIR__", gam_output_dir_r).replace("__ORIGINAL_OUTPUT_DIR__", original_output_dir_r).replace("__TEMP_MODEL_BASE_FILE__", temp_model_base_file_r)

with open(r_script_file, "w", encoding="utf-8") as f:
    f.write(r_script_content)

print(f"Saved R script: {r_script_file}")

# =============================================================================
# STEP 3: RUN R SCRIPT
# =============================================================================

r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if os.path.exists(r_path) and r_path not in os.environ["PATH"]:
    os.environ["PATH"] = os.environ["PATH"] + os.pathsep + r_path

print("\nChecking R availability...")
result = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)

if result.returncode != 0:
    print(result.stderr)
    raise RuntimeError("Rscript was not found. Check your R installation/path.")

print("Rscript found.")

print("\nRunning GAM A/B/C comparison R script...")
result = subprocess.run(["Rscript", r_script_file], capture_output=True, text=True)

if result.returncode != 0:
    print("\nR script failed.")
    print("\nSTDOUT:")
    print(result.stdout)
    print("\nSTDERR:")
    print(result.stderr)
    raise RuntimeError("Stopping because GAM comparison failed.")

print("R GAM comparison completed successfully.")
if result.stdout:
    print(result.stdout)

# =============================================================================
# STEP 4: VERIFY OUTPUTS
# =============================================================================

expected_outputs = [
    "Q1_GAM_A_B_C_model_performance.csv",
    "Q1_GAM_key_evidence_for_replacing_mixed_model.csv",
    "Q1_GAM_raw_LRT_A_vs_B.csv",
    "Q1_GAM_raw_LRT_B_vs_C.csv",
    "Q1_GAM_raw_LRT_A_vs_C.csv",
    "Q1_GAM_A_B_reduced_predictions_TET_for_supplement.csv",
    "Q1_GAM_A_B_C_comparison_summary.txt",
]

print("\nChecking outputs:")
for f in expected_outputs:
    path = os.path.join(gam_output_dir, f)
    print(f"  {f}: {'FOUND' if os.path.exists(path) else 'MISSING'}")

print("\nMost important file:")
print("  Q1_GAM_key_evidence_for_replacing_mixed_model.csv")
print("\nMost important row:")
print("  B_vs_C")
print("\nB_vs_C answers:")
print("  Do separate nonlinear T:ET smooths by WUE comparison × ecosystem improve fit?")
print("\nDONE.")
print(f"All new outputs saved in: {gam_output_dir}")