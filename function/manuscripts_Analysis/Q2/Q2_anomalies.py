"""
Q2_anomalies.py - COMPLETE FINAL VERSION (UPDATED)
Full replication of Malone's Q2 Anomaly workflow
EXACT OUTPUT NAMES - Uppercase WUE where Malone uses it

FIXES:
1. Added ci95_response to summary_by_class (matches Malone)
2. EXACT column ordering for summary_by_class (matches Malone)
3. Keep all original columns in model_data_long (matches Malone)
4. R coefficient tables now match Malone column names exactly
5. Added RDS files to expected outputs
6. Separate R outputs from Python outputs for verification
7. Clean directory structure

INPUT: data_products copy (verified identical to Water_WUE/data)
OUTPUTS: Separate from Malone with clean directory structure
"""

import os
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("Q2: Moisture Anomaly Analysis for WUE_ET and WUE_T")
print("EXACT MALONE REPLICATION (UPDATED)")
print("="*60)

# ============================================================================
# PATHS - SEPARATE FROM MALONE
# ============================================================================

# Input data - data_products copy (verified identical to Water_WUE/data)
input_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

# Python outputs - SEPARATE from Malone with clean organization
base_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2_anomalies"

# Three separate subdirectories
output_dir = os.path.join(base_output_dir, "Q2_anomalies_outputs")     # CSV, RDS, summary
figure_dir = os.path.join(base_output_dir, "Q2_anomalies_figures")     # PNG figures
temp_dir = os.path.join(base_output_dir, "Q2_anomalies_temp")          # Temporary R files

# Create all directories
for dir_path in [output_dir, figure_dir, temp_dir]:
    os.makedirs(dir_path, exist_ok=True)

print(f"\nOutput directory: {output_dir}")
print(f"Figure directory: {figure_dir}")
print(f"Temp directory:   {temp_dir}")
print(f"(Separate from Malone outputs - safe for comparison)")

# ============================================================================
# CLEAN OLD OUTPUTS AT START
# ============================================================================

print("\n" + "="*60)
print("STEP 0: Cleaning old outputs (Python directories only)")
print("="*60)

# R outputs (created by R script)
r_outputs = [
    "Q2_anomalies_model_base_wide.csv",
    "Q2_anomalies_model_data_long.csv",
    "Q2_anomalies_near_normal_baselines.csv",
    "Q2_anomalies_summary_by_class_metric_ecosystem.csv",
    "Q2_anomalies_summary_by_class_metric.csv",
    "Q2_anomalies_direction_summary.csv",
    "Q2_anomalies_wue_value_model.rds",
    "Q2_anomalies_signed_response_model.rds",
    "Q2_anomalies_response_magnitude_model.rds",
    "Q2_anomalies_response_direction_model.rds",
    "Q2_anomalies_wue_value_fixed_effects.csv",
    "Q2_anomalies_signed_response_fixed_effects.csv",
    "Q2_anomalies_response_magnitude_fixed_effects.csv",
    "Q2_anomalies_response_direction_fixed_effects.csv",
    "Q2_anomalies_wue_value_drop1_lrt.csv",
    "Q2_anomalies_signed_response_drop1_lrt.csv",
    "Q2_anomalies_response_magnitude_drop1_lrt.csv",
    "Q2_anomalies_response_direction_drop1_lrt.csv",
    "Q2_anomalies_model_results_table.csv",
    "Q2_anomalies_WUE_anomaly_class_letters.csv",
    "Q2_anomalies_WUE_anomaly_class_pairwise_contrasts.csv",
    "Q2_anomalies_signed_response_anomaly_class_letters.csv",
    "Q2_anomalies_signed_response_anomaly_class_pairwise_contrasts.csv",
    "Q2_anomalies_fitted_residuals.csv",
]

# Python outputs (created by Python after R)
python_outputs = [
    "Q2_anomalies_model_summary.txt",
]

all_expected_outputs = r_outputs + python_outputs

expected_figures = [
    "Q2_plot_01_WUE_by_anomaly_class.png",
    "Q2_plot_02_signed_WUE_response_by_anomaly.png",
    "Q2_plot_03_WUE_response_magnitude_by_anomaly.png",
    "Q2_plot_04_WUE_response_direction_by_anomaly.png",
    "Q2_plot_05_WUE_response_along_SPEI_gradient.png",
    "Q2_plot_06_anomaly_model_results_multipanel.png",
]

# Clean output files
for f in all_expected_outputs:
    path = os.path.join(output_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old output: {f}")

# Clean figure files
for f in expected_figures:
    path = os.path.join(figure_dir, f)
    if os.path.exists(path):
        os.remove(path)
        print(f"  Removed old figure: {f}")

# Clean temp directory
for f in os.listdir(temp_dir):
    path = os.path.join(temp_dir, f)
    if os.path.isfile(path):
        os.remove(path)
        print(f"  Removed old temp file: {f}")

# ============================================================================
# CONSTANTS
# ============================================================================

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
ECOSYSTEM_COLORS = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}
ANOMALY_COLORS = {
    "Dry": "#8C510A",
    "Near-normal": "gray70",
    "Wet": "#2B6CB0"
}
METRIC_COLORS = {
    "WUE_ET": "#333333",
    "WUE_T": "#6A3D9A"
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def classify_spei(x):
    """Classify SPEI values into Dry, Near-normal, or Wet"""
    conditions = [x < -1, (x >= -1) & (x <= 1), x > 1]
    choices = ["Dry", "Near-normal", "Wet"]
    return np.select(conditions, choices, default="Near-normal")

def write_table(df, filename):
    """Save dataframe to CSV"""
    df.to_csv(os.path.join(output_dir, filename), index=False)

def save_plot(fig, filename, width=10, height=7, dpi=300):
    """Save figure to figure_dir"""
    fig.set_size_inches(width, height)
    fig.savefig(os.path.join(figure_dir, filename), dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def theme_wue(base_size=13):
    """Malone-style theme"""
    sns.set_style("whitegrid")
    plt.rcParams.update({
        'font.size': base_size,
        'axes.labelsize': base_size,
        'axes.titlesize': base_size + 2,
        'xtick.labelsize': base_size - 1,
        'ytick.labelsize': base_size - 1,
        'legend.fontsize': base_size - 1,
        'figure.titlesize': base_size + 4,
    })

# ============================================================================
# DATA LOADING AND PREPARATION
# ============================================================================

print("\n" + "="*60)
print("STEP 1: Loading and preparing data")
print("="*60)

monthly = pd.read_csv(input_file)

# Trim whitespace
for col in monthly.select_dtypes(include='object').columns:
    monthly[col] = monthly[col].astype(str).str.strip()
    monthly[col] = monthly[col].replace('nan', np.nan)

required_cols = [
    "site_name", "Year", "month", "water_class",
    "WUE", "WUE_tra", "Trans_ratio", "SPEI_6", "SPEI_48"
]

missing_cols = [col for col in required_cols if col not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

# Filter data - EXACT MALONE
model_base = monthly.dropna(subset=required_cols).copy()
model_base = model_base[
    model_base['site_name'].notna() &
    (model_base['site_name'] != "") &
    model_base['water_class'].notna() &
    (model_base['water_class'] != "") &
    model_base['water_class'].isin(ECOSYSTEM_CLASSES) &
    np.isfinite(model_base['WUE']) &
    np.isfinite(model_base['WUE_tra']) &
    np.isfinite(model_base['Trans_ratio']) &
    (model_base['Trans_ratio'] >= 0) &
    (model_base['Trans_ratio'] <= 1) &
    np.isfinite(model_base['SPEI_6']) &
    np.isfinite(model_base['SPEI_48'])
].copy()

model_base['site_name'] = model_base['site_name'].astype('category')
model_base['water_class'] = pd.Categorical(model_base['water_class'], categories=ECOSYSTEM_CLASSES)
model_base['month_f'] = pd.Categorical(model_base['month'], categories=range(1, 13))
model_base['Trans_ratio_z'] = (model_base['Trans_ratio'] - model_base['Trans_ratio'].mean()) / model_base['Trans_ratio'].std()
model_base['WUE_ET'] = model_base['WUE']
model_base['WUE_T'] = model_base['WUE_tra']

print(f"  Data loaded: {len(model_base)} rows")
print(f"  Sites: {model_base['site_name'].nunique()}")
print(f"  Years: {model_base['Year'].min()} - {model_base['Year'].max()}")

# ============================================================================
# RESHAPE TO LONG FORMAT - KEEP ALL ORIGINAL COLUMNS
# ============================================================================

print("\n" + "="*60)
print("STEP 2: Reshaping to long format (keeping all columns)")
print("="*60)

# Keep ALL columns except the two WUE metric columns
id_vars = [c for c in model_base.columns if c not in ['WUE_ET', 'WUE_T']]

# First reshape: WUE_ET and WUE_T to long format
metric_long = pd.melt(
    model_base,
    id_vars=id_vars,
    value_vars=['WUE_ET', 'WUE_T'],
    var_name='metric',
    value_name='WUE_value'
)
metric_long['metric'] = pd.Categorical(metric_long['metric'], categories=['WUE_ET', 'WUE_T'])

# Second reshape: SPEI_6 and SPEI_48 to long format
spei_long_list = []

# SPEI_6
spei_6 = metric_long.copy()
spei_6['timescale'] = 'SPEI_6'
spei_6['SPEI_value'] = spei_6['SPEI_6']
spei_long_list.append(spei_6)

# SPEI_48
spei_48 = metric_long.copy()
spei_48['timescale'] = 'SPEI_48'
spei_48['SPEI_value'] = spei_48['SPEI_48']
spei_long_list.append(spei_48)

spei_long = pd.concat(spei_long_list, ignore_index=True)
spei_long['timescale'] = pd.Categorical(spei_long['timescale'], categories=['SPEI_6', 'SPEI_48'])
spei_long['anomaly_class'] = spei_long['SPEI_value'].apply(classify_spei)
spei_long['anomaly_class'] = pd.Categorical(spei_long['anomaly_class'], categories=['Dry', 'Near-normal', 'Wet'])
spei_long['water_class'] = pd.Categorical(spei_long['water_class'], categories=ECOSYSTEM_CLASSES)
spei_long['metric'] = pd.Categorical(spei_long['metric'], categories=['WUE_ET', 'WUE_T'])
spei_long['month_f'] = pd.Categorical(spei_long['month'], categories=range(1, 13))

print(f"  Long data: {len(spei_long)} rows")

# ============================================================================
# BASELINE CALCULATION - EXACT MALONE LOGIC
# ============================================================================

print("\n" + "="*60)
print("STEP 3: Calculating near-normal baselines")
print("="*60)

# Step 1: Calculate baseline for ALL Near-normal months (before filtering)
baseline_full = spei_long[spei_long['anomaly_class'] == 'Near-normal'].groupby(
    ['site_name', 'timescale', 'metric', 'month_f'], observed=True
).agg(
    baseline_mean=('WUE_value', 'mean'),
    baseline_median=('WUE_value', 'median'),
    baseline_n=('WUE_value', 'size')
).reset_index()

# Save full baseline FIRST (matches Malone's write_table before filtering)
write_table(baseline_full, "Q2_anomalies_near_normal_baselines.csv")

# Step 2: LEFT JOIN baseline to all data (matches Malone's left_join)
analysis_data = spei_long.merge(
    baseline_full,
    on=['site_name', 'timescale', 'metric', 'month_f'],
    how='left'
)

# Step 3: Filter to rows with valid baseline and baseline_n >= 3 (matches Malone)
analysis_data = analysis_data[
    np.isfinite(analysis_data['baseline_mean']) &
    (analysis_data['baseline_n'] >= 3)
].copy()

# Calculate response variables
analysis_data['WUE_response'] = analysis_data['WUE_value'] - analysis_data['baseline_mean']
analysis_data['WUE_response_magnitude'] = np.abs(analysis_data['WUE_response'])
analysis_data['response_direction'] = np.where(
    analysis_data['WUE_response'] < 0, 'Decrease', 'Increase_or_no_decrease'
)
analysis_data['decrease_binary'] = (analysis_data['response_direction'] == 'Decrease').astype(int)

# Anomaly data (exclude Near-normal)
anomaly_data = analysis_data[analysis_data['anomaly_class'] != 'Near-normal'].copy()
anomaly_data['anomaly_class'] = pd.Categorical(anomaly_data['anomaly_class'], categories=['Dry', 'Wet'])

# Save intermediate tables
write_table(model_base, "Q2_anomalies_model_base_wide.csv")
write_table(analysis_data, "Q2_anomalies_model_data_long.csv")

print(f"  Analysis data: {len(analysis_data)} rows")
print(f"  Anomaly data: {len(anomaly_data)} rows")

# ============================================================================
# SUMMARY STATISTICS - EXACT MALONE COLUMN ORDER
# ============================================================================

print("\n" + "="*60)
print("STEP 4: Calculating summary statistics")
print("="*60)

coverage_by_class = analysis_data.groupby(
    ['timescale', 'anomaly_class', 'metric', 'water_class'], observed=True
).agg(
    n_months=('WUE_value', 'size'),
    n_sites=('site_name', 'nunique'),
    mean_spei=('SPEI_value', 'mean'),
    mean_wue=('WUE_value', 'mean'),
    mean_response=('WUE_response', 'mean'),
    mean_abs_response=('WUE_response_magnitude', 'mean'),
    pct_decrease=('decrease_binary', lambda x: x.mean() * 100)
).reset_index()

# Calculate summary_by_class with ALL Malone columns
summary_by_class = analysis_data.groupby(
    ['timescale', 'anomaly_class', 'metric'], observed=True
).agg(
    n_months=('WUE_value', 'size'),
    n_sites=('site_name', 'nunique'),
    mean_wue=('WUE_value', 'mean'),
    median_wue=('WUE_value', 'median'),
    sd_wue=('WUE_value', 'std'),
    mean_response=('WUE_response', 'mean'),
    sd_response=('WUE_response', 'std'),  # Needed for se_response calculation
    mean_abs_response=('WUE_response_magnitude', 'mean'),
    pct_decrease=('decrease_binary', lambda x: x.mean() * 100)
).reset_index()

# Calculate SE and CI
summary_by_class['se_wue'] = summary_by_class['sd_wue'] / np.sqrt(summary_by_class['n_months'])
summary_by_class['ci95_wue'] = stats.t.ppf(0.975, summary_by_class['n_months'] - 1) * summary_by_class['se_wue']
summary_by_class['se_response'] = summary_by_class['sd_response'] / np.sqrt(summary_by_class['n_months'])
summary_by_class['ci95_response'] = stats.t.ppf(0.975, summary_by_class['n_months'] - 1) * summary_by_class['se_response']

# REORDER COLUMNS TO MATCH MALONE EXACTLY
# Malone columns: timescale, anomaly_class, metric, n_months, n_sites,
# mean_wue, median_wue, sd_wue, se_wue, ci95_wue,
# mean_response, se_response, ci95_response, mean_abs_response, pct_decrease
summary_by_class = summary_by_class[[
    "timescale", "anomaly_class", "metric",
    "n_months", "n_sites",
    "mean_wue", "median_wue", "sd_wue", "se_wue", "ci95_wue",
    "mean_response", "se_response", "ci95_response",
    "mean_abs_response", "pct_decrease"
]]

write_table(coverage_by_class, "Q2_anomalies_summary_by_class_metric_ecosystem.csv")
write_table(summary_by_class, "Q2_anomalies_summary_by_class_metric.csv")

direction_summary = anomaly_data.groupby(
    ['timescale', 'anomaly_class', 'metric', 'water_class'], observed=True
).agg(
    n_months=('WUE_value', 'size'),
    pct_decrease=('decrease_binary', lambda x: x.mean() * 100)
).reset_index()
write_table(direction_summary, "Q2_anomalies_direction_summary.csv")

# ============================================================================
# SAVE DATA FOR R - IN TEMP FOLDER
# ============================================================================

print("\n" + "="*60)
print("STEP 5: Saving data for R models")
print("="*60)

r_data_file = os.path.join(temp_dir, "temp_r_data.csv")
analysis_data.to_csv(r_data_file, index=False)
print(f"  Data saved to: {r_data_file}")

# ============================================================================
# CREATE R SCRIPT - WITH EXACT MALONE COLUMN NAMES (NO RENAMING)
# ============================================================================

print("\n" + "="*60)
print("STEP 6: Creating R script for mixed models")
print("="*60)

output_dir_r = output_dir.replace("\\", "/")
temp_dir_r = temp_dir.replace("\\", "/")
r_data_file_r = r_data_file.replace("\\", "/")

r_script_content = f'''
# ============================================================================
# Q2 ANOMALY MODELS - EXACT MALONE REPLICATION
# Called from Python
# EXACT FILENAMES MATCHING MALONE
# NO COLUMN RENAMING - Keep R default names
# ============================================================================

sep_line <- function() {{
    cat(paste(rep("=", 60), collapse=""))
    cat("\n")
}}

# Load required packages
library(lme4)
library(emmeans)
library(multcomp)
library(dplyr)

cat("\nAll required packages loaded successfully\n")

# Read data from Python
data <- read.csv("{r_data_file_r}")

# Convert to factors (matching Malone's R)
data$site_name <- as.factor(data$site_name)
data$water_class <- factor(data$water_class, levels = c("Upland", "Freshwater", "Saline"))
data$month_f <- as.factor(data$month_f)
data$metric <- factor(data$metric, levels = c("WUE_ET", "WUE_T"))
data$timescale <- factor(data$timescale, levels = c("SPEI_6", "SPEI_48"))
data$anomaly_class <- factor(data$anomaly_class, levels = c("Dry", "Near-normal", "Wet"))

# Anomaly data (Dry and Wet only)
anomaly_data <- data[data$anomaly_class != "Near-normal", ]
anomaly_data$anomaly_class <- factor(anomaly_data$anomaly_class, levels = c("Dry", "Wet"))

# ============================================================================
# Helper functions - KEEP R DEFAULT COLUMN NAMES (MATCH MALONE)
# ============================================================================

coef_table_lmer <- function(model) {{
    coef_df <- as.data.frame(summary(model)$coefficients)
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL
    # Keep R default column names: Estimate, Std. Error, df, t value, Pr(>|t|)
    coef_df <- coef_df[, c("term", setdiff(names(coef_df), "term"))]
    return(coef_df)
}}

coef_table_glmer <- function(model) {{
    coef_df <- as.data.frame(summary(model)$coefficients)
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL
    # Keep R default column names: Estimate, Std. Error, z value, Pr(>|z|)
    coef_df <- coef_df[, c("term", setdiff(names(coef_df), "term"))]
    return(coef_df)
}}

safe_drop1 <- function(model) {{
    out <- tryCatch(
        drop1(model, test = "Chisq"),
        error = function(e) {{
            message("drop1 failed: ", conditionMessage(e))
            NULL
        }}
    )
    if (is.null(out)) return(NULL)
    out_df <- as.data.frame(out)
    out_df$term <- rownames(out_df)
    rownames(out_df) <- NULL
    out_df <- out_df[, c("term", setdiff(names(out_df), "term"))]
    return(out_df)
}}

write_emmeans_letters <- function(emmeans_obj, filename) {{
    letters_df <- as.data.frame(cld(emmeans_obj, adjust = "tukey", Letters = letters))
    letters_df$.group <- gsub(" ", "", letters_df$.group)
    write.csv(letters_df, file.path("{output_dir_r}", filename), row.names = FALSE)
}}

write_pairwise_contrasts <- function(contrasts_obj, filename) {{
    pairs_df <- as.data.frame(contrasts_obj)
    write.csv(pairs_df, file.path("{output_dir_r}", filename), row.names = FALSE)
}}

cat("\n")
sep_line()
cat("PART 1: MIXED MODELS")
cat("\n")
sep_line()
cat("\n")

# ----------------------------------------------------------------------------
# Model 1: WUE_value (lmer)
# ----------------------------------------------------------------------------
cat("\n--- Model 1: WUE_value ---\n")

wue_model <- lmer(
    WUE_value ~ metric * timescale * anomaly_class * water_class +
        Trans_ratio_z + month_f + (1 | site_name),
    data = data,
    REML = FALSE
)

saveRDS(wue_model, file.path("{output_dir_r}", "Q2_anomalies_wue_value_model.rds"))

wue_fixed <- coef_table_lmer(wue_model)
write.csv(wue_fixed, file.path("{output_dir_r}", "Q2_anomalies_wue_value_fixed_effects.csv"), row.names = FALSE)

wue_drop1 <- safe_drop1(wue_model)
if (!is.null(wue_drop1)) {{
    write.csv(wue_drop1, file.path("{output_dir_r}", "Q2_anomalies_wue_value_drop1_lrt.csv"), row.names = FALSE)
}}

# EXACT MALONE UPPERCASE WUE FILENAME
wue_emmeans <- emmeans(wue_model, specs = ~ anomaly_class | timescale * metric)
write_emmeans_letters(wue_emmeans, "Q2_anomalies_WUE_anomaly_class_letters.csv")

# EXACT MALONE UPPERCASE WUE FILENAME
wue_pairs <- contrast(wue_emmeans, method = "pairwise", adjust = "tukey")
write_pairwise_contrasts(wue_pairs, "Q2_anomalies_WUE_anomaly_class_pairwise_contrasts.csv")

# Fitted values and residuals
data$wue_fitted <- fitted(wue_model)
data$wue_residual <- residuals(wue_model)
write.csv(data[, c("site_name", "Year", "month", "water_class", "timescale", "anomaly_class",
                   "SPEI_value", "metric", "WUE_value", "baseline_mean", "WUE_response",
                   "WUE_response_magnitude", "response_direction", "wue_fitted", "wue_residual")],
          file.path("{output_dir_r}", "Q2_anomalies_fitted_residuals.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# Model 2: Signed Response (lmer)
# ----------------------------------------------------------------------------
cat("\n--- Model 2: Signed Response ---\n")

response_model <- lmer(
    WUE_response ~ metric * timescale * anomaly_class * water_class +
        Trans_ratio_z + month_f + (1 | site_name),
    data = data,
    REML = FALSE
)

saveRDS(response_model, file.path("{output_dir_r}", "Q2_anomalies_signed_response_model.rds"))

response_fixed <- coef_table_lmer(response_model)
write.csv(response_fixed, file.path("{output_dir_r}", "Q2_anomalies_signed_response_fixed_effects.csv"), row.names = FALSE)

response_drop1 <- safe_drop1(response_model)
if (!is.null(response_drop1)) {{
    write.csv(response_drop1, file.path("{output_dir_r}", "Q2_anomalies_signed_response_drop1_lrt.csv"), row.names = FALSE)
}}

response_emmeans <- emmeans(response_model, specs = ~ anomaly_class | timescale * metric)
write_emmeans_letters(response_emmeans, "Q2_anomalies_signed_response_anomaly_class_letters.csv")

response_pairs <- contrast(response_emmeans, method = "pairwise", adjust = "tukey")
write_pairwise_contrasts(response_pairs, "Q2_anomalies_signed_response_anomaly_class_pairwise_contrasts.csv")

# ----------------------------------------------------------------------------
# Model 3: Response Magnitude (lmer)
# ----------------------------------------------------------------------------
cat("\n--- Model 3: Response Magnitude ---\n")

magnitude_model <- lmer(
    WUE_response_magnitude ~ metric * timescale * anomaly_class * water_class +
        Trans_ratio_z + month_f + (1 | site_name),
    data = anomaly_data,
    REML = FALSE
)

saveRDS(magnitude_model, file.path("{output_dir_r}", "Q2_anomalies_response_magnitude_model.rds"))

magnitude_fixed <- coef_table_lmer(magnitude_model)
write.csv(magnitude_fixed, file.path("{output_dir_r}", "Q2_anomalies_response_magnitude_fixed_effects.csv"), row.names = FALSE)

magnitude_drop1 <- safe_drop1(magnitude_model)
if (!is.null(magnitude_drop1)) {{
    write.csv(magnitude_drop1, file.path("{output_dir_r}", "Q2_anomalies_response_magnitude_drop1_lrt.csv"), row.names = FALSE)
}}

# ----------------------------------------------------------------------------
# Model 4: Direction (Binomial GLMM - glmer)
# ----------------------------------------------------------------------------
cat("\n--- Model 4: Direction (Binomial) ---\n")

direction_model <- glmer(
    decrease_binary ~ metric * timescale * anomaly_class * water_class +
        Trans_ratio_z + month_f + (1 | site_name),
    data = anomaly_data,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

saveRDS(direction_model, file.path("{output_dir_r}", "Q2_anomalies_response_direction_model.rds"))

direction_fixed <- coef_table_glmer(direction_model)
write.csv(direction_fixed, file.path("{output_dir_r}", "Q2_anomalies_response_direction_fixed_effects.csv"), row.names = FALSE)

direction_drop1 <- safe_drop1(direction_model)
if (!is.null(direction_drop1)) {{
    write.csv(direction_drop1, file.path("{output_dir_r}", "Q2_anomalies_response_direction_drop1_lrt.csv"), row.names = FALSE)
}}

# ----------------------------------------------------------------------------
# Model summary table
# ----------------------------------------------------------------------------
model_results <- rbind(
    data.frame(
        model = "WUE value",
        AIC = AIC(wue_model),
        BIC = BIC(wue_model),
        logLik = as.numeric(logLik(wue_model)),
        response = "WUE_value"
    ),
    data.frame(
        model = "Signed response",
        AIC = AIC(response_model),
        BIC = BIC(response_model),
        logLik = as.numeric(logLik(response_model)),
        response = "WUE_response"
    ),
    data.frame(
        model = "Response magnitude",
        AIC = AIC(magnitude_model),
        BIC = BIC(magnitude_model),
        logLik = as.numeric(logLik(magnitude_model)),
        response = "abs(WUE_response)"
    ),
    data.frame(
        model = "Direction",
        AIC = AIC(direction_model),
        BIC = BIC(direction_model),
        logLik = as.numeric(logLik(direction_model)),
        response = "decrease_binary"
    )
)
write.csv(model_results, file.path("{output_dir_r}", "Q2_anomalies_model_results_table.csv"), row.names = FALSE)

cat("\n")
sep_line()
cat("R MODELS COMPLETE")
cat("\n")
sep_line()
cat("\n")
'''

# Save R script in temp folder
r_script_file = os.path.join(temp_dir, "run_models.R")
with open(r_script_file, 'w', encoding='utf-8') as f:
    f.write(r_script_content)

print(f"  R script saved to: {r_script_file}")

# ============================================================================
# VERIFY TEMP FILES EXIST
# ============================================================================

print("\nVerifying temp files exist before running R:")
print(f"  temp_r_data.csv exists: {os.path.exists(r_data_file)}")
print(f"  run_models.R exists: {os.path.exists(r_script_file)}")

if not os.path.exists(r_data_file):
    raise FileNotFoundError(f"Temp data file missing: {r_data_file}")
if not os.path.exists(r_script_file):
    raise FileNotFoundError(f"R script file missing: {r_script_file}")

# ============================================================================
# RUN R SCRIPT
# ============================================================================

print("\n" + "="*60)
print("STEP 7: Running R models")
print("="*60)

# Add R to PATH
r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
if r_path not in os.environ['PATH']:
    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + r_path

try:
    result = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)
    print("  R is available")
except Exception as e:
    print(f"  R not found! Error: {e}")
    exit()

# Run R script
print("\n  Running R models...")
result = subprocess.run(["Rscript", r_script_file], capture_output=True, text=True)

if result.returncode != 0:
    print("  R script had errors:")
    print(result.stderr)
    raise RuntimeError("R analysis failed.")
else:
    print("  R models completed successfully!")
    if result.stdout:
        print(result.stdout)

# ============================================================================
# VERIFY R OUTPUTS - HARD FAILURE IF MISSING (ONLY R OUTPUTS)
# ============================================================================

print("\n" + "="*60)
print("STEP 8: Verifying R outputs")
print("="*60)

missing_r_outputs = [
    f for f in r_outputs
    if not os.path.exists(os.path.join(output_dir, f))
]

if missing_r_outputs:
    raise FileNotFoundError(
        "Missing expected R output files:\n" + "\n".join(missing_r_outputs)
    )
else:
    print(f"All {len(r_outputs)} expected R outputs created.")

# ============================================================================
# CREATE FIGURES - EXACT MALONE NAMES (6 FIGURES ONLY)
# ============================================================================

print("\n" + "="*60)
print("STEP 9: Creating figures (Malone-style)")
print("="*60)

theme_wue()

# Figure 1: WUE by Anomaly Class
print("  Creating Figure 1...")
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

for idx, timescale in enumerate(['SPEI_6', 'SPEI_48']):
    ax = axes[idx]
    subset = summary_by_class[summary_by_class['timescale'] == timescale]
    
    for metric in ['WUE_ET', 'WUE_T']:
        metric_subset = subset[subset['metric'] == metric].copy()
        anomaly_order = ['Dry', 'Near-normal', 'Wet']
        metric_subset = metric_subset.set_index('anomaly_class').reindex(anomaly_order).reset_index()
        
        ax.errorbar(
            metric_subset['anomaly_class'],
            metric_subset['mean_wue'],
            yerr=metric_subset['ci95_wue'],
            fmt='o-',
            color=METRIC_COLORS[metric],
            label=metric,
            capsize=5,
            markersize=8,
            linewidth=2
        )
    
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    ax.set_title(f'SPEI-{timescale.split("_")[1]}')
    ax.set_xlabel('Moisture anomaly class')
    if idx == 0:
        ax.set_ylabel('Mean WUE')
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.suptitle('WUE_ET and WUE_T Across Moisture Anomaly Classes', fontsize=14, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_01_WUE_by_anomaly_class.png", width=12, height=6)

# Figure 2: Signed Response
print("  Creating Figure 2...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx_t, timescale in enumerate(['SPEI_6', 'SPEI_48']):
    for idx_m, metric in enumerate(['WUE_ET', 'WUE_T']):
        ax = axes[idx_m, idx_t]
        subset = analysis_data[
            (analysis_data['timescale'] == timescale) &
            (analysis_data['metric'] == metric)
        ]
        
        data_by_class = []
        positions = []
        colors = []
        anomaly_order = ['Dry', 'Near-normal', 'Wet']
        
        for i, anomaly in enumerate(anomaly_order):
            for j, eco in enumerate(ECOSYSTEM_CLASSES):
                data = subset[
                    (subset['anomaly_class'] == anomaly) &
                    (subset['water_class'] == eco)
                ]['WUE_response'].dropna().values
                if len(data) > 0:
                    data_by_class.append(data)
                    positions.append(i * 3 + j)
                    colors.append(ECOSYSTEM_COLORS[eco])
        
        if data_by_class:
            bp = ax.boxplot(data_by_class, positions=positions, widths=0.6, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.5)
        
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
        ax.set_xticks([0.5, 3.5, 6.5])
        ax.set_xticklabels(['Dry', 'Near-normal', 'Wet'])
        ax.set_title(f'{metric} - {timescale}')
        ax.set_xlabel('Moisture anomaly class')
        if idx_m == 0:
            ax.set_ylabel('WUE response')
        ax.grid(True, alpha=0.3)

plt.suptitle('Signed WUE Responses Relative to Near-Normal Baselines', fontsize=14, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_02_signed_WUE_response_by_anomaly.png", width=12, height=10)

# Figure 3: Response Magnitude
print("  Creating Figure 3...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx_t, timescale in enumerate(['SPEI_6', 'SPEI_48']):
    for idx_m, metric in enumerate(['WUE_ET', 'WUE_T']):
        ax = axes[idx_m, idx_t]
        subset = anomaly_data[
            (anomaly_data['timescale'] == timescale) &
            (anomaly_data['metric'] == metric)
        ]
        
        data_by_class = []
        positions = []
        colors = []
        anomaly_order = ['Dry', 'Wet']
        
        for i, anomaly in enumerate(anomaly_order):
            for j, eco in enumerate(ECOSYSTEM_CLASSES):
                data = subset[
                    (subset['anomaly_class'] == anomaly) &
                    (subset['water_class'] == eco)
                ]['WUE_response_magnitude'].dropna().values
                if len(data) > 0:
                    data_by_class.append(data)
                    positions.append(i * 3 + j)
                    colors.append(ECOSYSTEM_COLORS[eco])
        
        if data_by_class:
            bp = ax.boxplot(data_by_class, positions=positions, widths=0.6, patch_artist=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.5)
        
        ax.set_xticks([0.5, 3.5])
        ax.set_xticklabels(['Dry', 'Wet'])
        ax.set_title(f'{metric} - {timescale}')
        ax.set_xlabel('Moisture anomaly class')
        if idx_m == 0:
            ax.set_ylabel('|WUE response|')
        ax.grid(True, alpha=0.3)

plt.suptitle('Magnitude of WUE Responses During Dry and Wet Anomalies', fontsize=14, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_03_WUE_response_magnitude_by_anomaly.png", width=12, height=10)

# Figure 4: Direction Summary
print("  Creating Figure 4...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx_t, timescale in enumerate(['SPEI_6', 'SPEI_48']):
    for idx_m, metric in enumerate(['WUE_ET', 'WUE_T']):
        ax = axes[idx_m, idx_t]
        subset = direction_summary[
            (direction_summary['timescale'] == timescale) &
            (direction_summary['metric'] == metric)
        ]
        
        anomalies = subset['anomaly_class'].unique()
        x = np.arange(len(anomalies))
        width = 0.25
        
        for i, eco in enumerate(ECOSYSTEM_CLASSES):
            eco_subset = subset[subset['water_class'] == eco]
            if len(eco_subset) > 0:
                values = []
                for anomaly in anomalies:
                    val = eco_subset[eco_subset['anomaly_class'] == anomaly]['pct_decrease'].values
                    values.append(val[0] if len(val) > 0 else 0)
                ax.bar(
                    x + (i - 1) * width,
                    values,
                    width,
                    label=eco,
                    color=ECOSYSTEM_COLORS[eco],
                    alpha=0.7
                )
        
        ax.set_xticks(x)
        ax.set_xticklabels(anomalies)
        ax.set_ylim(0, 100)
        if idx_m == 0:
            ax.set_ylabel('Months with WUE decrease (%)')
        ax.set_title(f'{metric} - {timescale}')
        if idx_m == 0 and idx_t == 0:
            ax.legend()
        ax.grid(True, alpha=0.3)

plt.suptitle('Probability of WUE Decrease During Moisture Anomalies', fontsize=14, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_04_WUE_response_direction_by_anomaly.png", width=12, height=10)

# Figure 5: SPEI Gradient
print("  Creating Figure 5...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for idx_t, timescale in enumerate(['SPEI_6', 'SPEI_48']):
    for idx_m, metric in enumerate(['WUE_ET', 'WUE_T']):
        ax = axes[idx_m, idx_t]
        subset = analysis_data[
            (analysis_data['timescale'] == timescale) &
            (analysis_data['metric'] == metric)
        ]
        
        for eco in ECOSYSTEM_CLASSES:
            eco_subset = subset[subset['water_class'] == eco]
            if len(eco_subset) > 0:
                ax.scatter(
                    eco_subset['SPEI_value'],
                    eco_subset['WUE_response'],
                    color=ECOSYSTEM_COLORS[eco],
                    alpha=0.15,
                    s=10
                )
                try:
                    from statsmodels.nonparametric.smoothers_lowess import lowess
                    if len(eco_subset) > 10:
                        smoothed = lowess(eco_subset['WUE_response'], eco_subset['SPEI_value'], frac=0.3)
                        ax.plot(smoothed[:, 0], smoothed[:, 1], color=ECOSYSTEM_COLORS[eco], linewidth=2, label=eco)
                except:
                    pass
        
        ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.5)
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.5)
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
        ax.set_xlabel('SPEI value')
        if idx_m == 0:
            ax.set_ylabel('WUE response')
        ax.set_title(f'{metric} - {timescale}')
        if idx_m == 0 and idx_t == 0:
            ax.legend()
        ax.grid(True, alpha=0.3)

plt.suptitle('WUE Response Along Short- and Long-Term SPEI Gradients', fontsize=14, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_05_WUE_response_along_SPEI_gradient.png", width=12, height=10)

# Figure 6: Composite Multipanel - EXACT MALONE
print("  Creating Figure 6: Composite multipanel...")
fig = plt.figure(figsize=(13, 10))

# Panel A: WUE by Anomaly Class
ax_a = fig.add_subplot(2, 1, 1)
for timescale in ['SPEI_6', 'SPEI_48']:
    subset = summary_by_class[summary_by_class['timescale'] == timescale]
    for metric in ['WUE_ET', 'WUE_T']:
        metric_subset = subset[subset['metric'] == metric].copy()
        anomaly_order = ['Dry', 'Near-normal', 'Wet']
        metric_subset = metric_subset.set_index('anomaly_class').reindex(anomaly_order).reset_index()
        ax_a.errorbar(
            metric_subset['anomaly_class'],
            metric_subset['mean_wue'],
            yerr=metric_subset['ci95_wue'],
            fmt='o-',
            color=METRIC_COLORS[metric],
            label=f'{metric} ({timescale})',
            capsize=5,
            markersize=8,
            linewidth=2
        )
ax_a.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
ax_a.set_title('A. WUE_ET and WUE_T across moisture anomaly classes')
ax_a.set_xlabel('Moisture anomaly class')
ax_a.set_ylabel('Mean WUE')
ax_a.legend(loc='upper right')
ax_a.grid(True, alpha=0.3)

# Panel B: Direction Heatmap
ax_b = fig.add_subplot(2, 1, 2)
heat_data = direction_summary.pivot_table(
    index=['timescale', 'metric', 'anomaly_class'],
    columns='water_class',
    values='pct_decrease'
).fillna(0)

im = ax_b.imshow(heat_data.values, cmap='Blues', aspect='auto', vmin=0, vmax=100)
ax_b.set_xticks(range(len(heat_data.columns)))
ax_b.set_xticklabels(heat_data.columns, rotation=25, ha='right')
ax_b.set_yticks(range(len(heat_data.index)))
ax_b.set_yticklabels([f'{i[0]}_{i[1]}_{i[2]}' for i in heat_data.index], fontsize=8)
ax_b.set_title('B. Directional WUE response summary')
ax_b.set_ylabel('Timescale_Metric_Anomaly')
ax_b.set_xlabel('Ecosystem class')

for i in range(len(heat_data.index)):
    for j in range(len(heat_data.columns)):
        text = f"{heat_data.values[i, j]:.0f}%"
        ax_b.text(j, i, text, ha='center', va='center', color='white' if heat_data.values[i, j] > 50 else 'black')

plt.colorbar(im, ax=ax_b, label='Decrease (%)')
plt.suptitle('Q2 Moisture Anomalies and WUE Responses', fontsize=16, fontweight='bold')
plt.tight_layout()
save_plot(fig, "Q2_plot_06_anomaly_model_results_multipanel.png", width=13, height=10)

# ============================================================================
# VERIFY FIGURES - HARD FAILURE IF MISSING
# ============================================================================

print("\n" + "="*60)
print("STEP 10: Verifying figures")
print("="*60)

missing_figures = [
    f for f in expected_figures
    if not os.path.exists(os.path.join(figure_dir, f))
]

if missing_figures:
    raise FileNotFoundError(
        "Missing expected figure files:\n" + "\n".join(missing_figures)
    )
else:
    print(f"All {len(expected_figures)} expected figures created.")

# ============================================================================
# SUMMARY TEXT - EXACT MALONE (CREATED BY PYTHON)
# ============================================================================

print("\n" + "="*60)
print("STEP 11: Creating summary text")
print("="*60)

with open(os.path.join(output_dir, "Q2_anomalies_model_summary.txt"), 'w') as f:
    f.write("Q2 moisture-anomaly analysis for WUE_ET and WUE_T\n")
    f.write(f"Input: {input_file}\n")
    f.write(f"Output: {output_dir}\n\n")
    f.write(f"Rows in source monthly data: {len(monthly)}\n")
    f.write(f"Rows after complete-case and ecosystem filters: {len(model_base)}\n")
    f.write(f"Rows in metric x timescale long data: {len(analysis_data)}\n")
    f.write(f"Rows in dry/wet anomaly subset: {len(anomaly_data)}\n")
    f.write(f"Sites: {analysis_data['site_name'].nunique()}\n")
    f.write(f"Years: {analysis_data['Year'].min()} - {analysis_data['Year'].max()}\n\n")
    
    f.write("Anomaly definitions:\n")
    f.write("  Dry: SPEI < -1\n")
    f.write("  Near-normal: -1 <= SPEI <= 1\n")
    f.write("  Wet: SPEI > 1\n\n")
    
    f.write("Coverage by timescale/anomaly/metric/ecosystem:\n")
    f.write(coverage_by_class.to_string())
    f.write("\n\n")
    
    f.write("Model results summary:\n")
    model_results_file = os.path.join(output_dir, "Q2_anomalies_model_results_table.csv")
    if os.path.exists(model_results_file):
        model_results = pd.read_csv(model_results_file)
        f.write(model_results.to_string())
        f.write("\n")

print("  Saved: Q2_anomalies_model_summary.txt")

# ============================================================================
# FINAL VERIFICATION - ALL OUTPUTS INCLUDING SUMMARY
# ============================================================================

print("\n" + "="*60)
print("STEP 12: Final verification of all outputs")
print("="*60)

missing_all = [
    f for f in all_expected_outputs
    if not os.path.exists(os.path.join(output_dir, f))
]

if missing_all:
    print(f"WARNING: Some outputs are still missing:")
    for f in missing_all:
        print(f"  - {f}")
else:
    print(f"All {len(all_expected_outputs)} expected outputs created successfully.")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*60)
print("Q2 Anomaly Analysis COMPLETE!")
print("="*60)

print(f"\nAll outputs saved to:")
print(f"  Outputs (CSV/RDS/summary): {output_dir}")
print(f"  Figures (PNG):             {figure_dir}")
print(f"  Temp files:                {temp_dir}")

print("\nFiles generated:")
print("  CSV Files (Python):")
print("  - Q2_anomalies_model_base_wide.csv")
print("  - Q2_anomalies_model_data_long.csv")
print("  - Q2_anomalies_near_normal_baselines.csv")
print("  - Q2_anomalies_summary_by_class_metric.csv")
print("  - Q2_anomalies_summary_by_class_metric_ecosystem.csv")
print("  - Q2_anomalies_direction_summary.csv")
print("  - Q2_anomalies_model_summary.txt")

print("\n  CSV Files (R models):")
print("  - Q2_anomalies_wue_value_fixed_effects.csv")
print("  - Q2_anomalies_signed_response_fixed_effects.csv")
print("  - Q2_anomalies_response_magnitude_fixed_effects.csv")
print("  - Q2_anomalies_response_direction_fixed_effects.csv")
print("  - Q2_anomalies_wue_value_drop1_lrt.csv")
print("  - Q2_anomalies_signed_response_drop1_lrt.csv")
print("  - Q2_anomalies_response_magnitude_drop1_lrt.csv")
print("  - Q2_anomalies_response_direction_drop1_lrt.csv")
print("  - Q2_anomalies_WUE_anomaly_class_letters.csv")
print("  - Q2_anomalies_WUE_anomaly_class_pairwise_contrasts.csv")
print("  - Q2_anomalies_signed_response_anomaly_class_letters.csv")
print("  - Q2_anomalies_signed_response_anomaly_class_pairwise_contrasts.csv")
print("  - Q2_anomalies_fitted_residuals.csv")
print("  - Q2_anomalies_model_results_table.csv")

print("\n  Figures (6 files - Malone-style):")
print("  - Q2_plot_01_WUE_by_anomaly_class.png")
print("  - Q2_plot_02_signed_WUE_response_by_anomaly.png")
print("  - Q2_plot_03_WUE_response_magnitude_by_anomaly.png")
print("  - Q2_plot_04_WUE_response_direction_by_anomaly.png")
print("  - Q2_plot_05_WUE_response_along_SPEI_gradient.png")
print("  - Q2_plot_06_anomaly_model_results_multipanel.png")

print("\n" + "="*60)
print("✓ Q2 Anomaly workflow complete - EXACT MALONE REPLICATION")
print("  - Same output file names as Malone (uppercase WUE)")
print("  - Same figure file names as Malone (6 figures only)")
print("  - Separate directories (safe from Malone)")
print("  - Cleaned old outputs before run")
print("  - Verified all expected outputs created (hard failure if missing)")
print("  - Using verified data_products input (identical to Water_WUE/data)")
print("  - EXACT column order for summary_by_class (matches Malone)")
print("  - Keep all original columns in model_data_long")
print("  - R coefficient tables use Malone's default column names")
print("="*60)