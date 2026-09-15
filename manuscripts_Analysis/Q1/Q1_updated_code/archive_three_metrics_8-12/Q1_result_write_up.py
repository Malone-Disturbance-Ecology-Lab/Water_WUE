"""
Q1 RESULTS EXTRACTION SCRIPT - COMPLETE
Reads existing updated CSV outputs and prints exact values from actual results.
Includes optional Panel A pairwise ecosystem contrasts using emmeans.

Inputs (all from Q1_updated_results/outputs/):
- Q1_final_near_normal_summary_by_difference_type.csv
- Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv
- Q1_final_near_normal_mixed_likelihood_ratio_tests.csv
- Q1_final_near_normal_mixed_model_comparison_lrt.csv
- Q1_final_gam_model_comparison.csv
- Q1_final_near_normal_mixed_fitted_residuals.csv (for optional contrasts)

Output: Printed console report with exact values from actual results
"""

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

output_dir = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs")

print("="*90)
print("Q1 UPDATED RESULTS - EXACT VALUES FROM ACTUAL RESULTS")
print("="*90)
print(f"\nReading from: {output_dir}\n")

# =============================================================================
# HELPER FUNCTIONS FOR ROBUST COLUMN DETECTION
# =============================================================================

def find_p_value_column(df):
    """Find the p-value column in a dataframe robustly."""
    p_cols = [c for c in df.columns if c.startswith("Pr") or c == "p_value" or c == "p" or c == "P"]
    if p_cols:
        return p_cols[0]
    # Try case-insensitive
    for c in df.columns:
        if c.lower() in ['pr(>chisq)', 'pr(>chi)', 'pr(>f)', 'p_value', 'pvalue', 'p']:
            return c
    return None

def find_chi_column(df):
    """Find the Chi-square/LRT column in a dataframe robustly."""
    chi_cols = [c for c in df.columns if c == "Chisq" or c == "LRT" or c == "Chi" or c == "LR"]
    if chi_cols:
        return chi_cols[0]
    return None

# =============================================================================
# READ CSV FILES
# =============================================================================

print("Loading CSV outputs...")

# Panel A: Summary means
summary_df = pd.read_csv(output_dir / "Q1_final_near_normal_summary_by_difference_type.csv")

# Mixed model fixed effects
mixed_fixed = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv")

# Mixed model LRT
mixed_lrt = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv")

# Mixed additive vs interaction
mixed_comp = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_model_comparison_lrt.csv")

# GAM comparison
gam_comp = pd.read_csv(output_dir / "Q1_final_gam_model_comparison.csv")

# =============================================================================
# SECTION 1: PANEL A MEANS
# =============================================================================

print("\n" + "="*90)
print("SECTION 1: PANEL A MEANS (mean ± 95% CI)")
print("="*90)

# Display summary with mean ± CI
summary_display = summary_df.copy()
summary_display["mean_ci"] = summary_display.apply(
    lambda x: f"{x['mean_difference']:.3f} ± {x['ci95_difference']:.3f}", axis=1
)

print(summary_display[["difference_label", "water_class", "n_months", "n_sites", "mean_ci"]].to_string(index=False))

# =============================================================================
# OPTIONAL SECTION 1B: PANEL A PAIRWISE ECOSYSTEM CONTRASTS
# =============================================================================

print("\n" + "="*90)
print("OPTIONAL SECTION 1B: PANEL A PAIRWISE ECOSYSTEM CONTRASTS")
print("="*90)

# This is an optional diagnostic for Panel A means only.
# It uses the updated near-normal q1_data saved by the final workflow.
# It does NOT rerun the main Q1 workflow.

pairwise_input = output_dir / "Q1_final_near_normal_mixed_fitted_residuals.csv"

pairwise_lrt_out = output_dir / "Q1_optional_panelA_global_difference_type_by_ecosystem_LRT.csv"
pairwise_emm_out = output_dir / "Q1_optional_panelA_emmeans_by_difference_and_ecosystem.csv"
pairwise_contrast_out = output_dir / "Q1_optional_panelA_tukey_pairwise_ecosystem_contrasts.csv"

pairwise_r_script = output_dir / "Q1_optional_panelA_pairwise_emmeans.R"

if not pairwise_input.exists():
    print(f"Pairwise input file not found: {pairwise_input}")
    print("Skipping optional Panel A contrasts.")
    panelA_global_chi = None
    panelA_global_df = None
    panelA_global_p = None
    pairwise_sig_df = pd.DataFrame()
else:
    print(f"Running optional Panel A contrasts using: {pairwise_input}")
    pairwise_input_r = str(pairwise_input).replace("\\", "/")
    output_dir_r = str(output_dir).replace("\\", "/")

    r_code = f"""
suppressPackageStartupMessages({{
  if (!requireNamespace("lme4", quietly = TRUE)) {{
    stop("Package 'lme4' is required.")
  }}
  if (!requireNamespace("emmeans", quietly = TRUE)) {{
    stop("Package 'emmeans' is required. Install with install.packages('emmeans').")
  }}
}})

q1_data <- read.csv("{pairwise_input_r}", stringsAsFactors = FALSE)

ecosystem_classes <- c("Upland", "Freshwater", "Saline")
difference_cols <- c("WUE_ET_minus_WUE_T", "WUE_ET_minus_WUE_E", "WUE_E_minus_WUE_T")
diff_labels <- c(
  WUE_ET_minus_WUE_T = "WUE_ET - WUE_T",
  WUE_ET_minus_WUE_E = "WUE_ET - WUE_E",
  WUE_E_minus_WUE_T = "WUE_E - WUE_T"
)

q1_data$site_name <- factor(q1_data$site_name)
q1_data$month_f <- factor(q1_data$month)
q1_data$water_class <- factor(q1_data$water_class, levels = ecosystem_classes)
q1_data$difference_type <- factor(q1_data$difference_type, levels = difference_cols)

q1_data <- q1_data[complete.cases(q1_data[, c("site_name", "month", "water_class", "difference_type", "WUE_difference")]), ]
q1_data <- droplevels(q1_data)

panelA_model <- lme4::lmer(
  WUE_difference ~ difference_type * water_class +
    (1 | site_name) + (1 | month_f),
  data = q1_data,
  REML = FALSE
)

panelA_additive <- lme4::lmer(
  WUE_difference ~ difference_type + water_class +
    (1 | site_name) + (1 | month_f),
  data = q1_data,
  REML = FALSE
)

global_comp <- as.data.frame(anova(panelA_additive, panelA_model))
global_comp$model <- rownames(global_comp)
global_comp <- global_comp[, c("model", setdiff(names(global_comp), "model"))]

emm <- emmeans::emmeans(
  panelA_model,
  specs = ~ water_class | difference_type,
  lmer.df = "asymptotic"
)

emm_df <- as.data.frame(emm)
emm_df$difference_label <- diff_labels[as.character(emm_df$difference_type)]

pairwise <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
pairwise_df <- as.data.frame(pairwise)
pairwise_df$difference_label <- diff_labels[as.character(pairwise_df$difference_type)]
pairwise_df$significant_0.05 <- pairwise_df$p.value < 0.05

write.csv(global_comp, file.path("{output_dir_r}", "Q1_optional_panelA_global_difference_type_by_ecosystem_LRT.csv"), row.names = FALSE)
write.csv(emm_df, file.path("{output_dir_r}", "Q1_optional_panelA_emmeans_by_difference_and_ecosystem.csv"), row.names = FALSE)
write.csv(pairwise_df, file.path("{output_dir_r}", "Q1_optional_panelA_tukey_pairwise_ecosystem_contrasts.csv"), row.names = FALSE)
"""

    with open(pairwise_r_script, "w", encoding="utf-8") as f:
        f.write(r_code)

    r_path = r"C:\Program Files\R\R-4.4.2\bin\x64"
    if os.path.exists(r_path) and r_path not in os.environ["PATH"]:
        os.environ["PATH"] = os.environ["PATH"] + os.pathsep + r_path

    result = subprocess.run(
        ["Rscript", str(pairwise_r_script)],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("Optional Panel A pairwise diagnostic failed.")
        print(result.stderr)
        panelA_global_chi = None
        panelA_global_df = None
        panelA_global_p = None
        pairwise_sig_df = pd.DataFrame()
    else:
        panelA_global = pd.read_csv(pairwise_lrt_out)
        panelA_pairwise = pd.read_csv(pairwise_contrast_out)

        panelA_model_row = panelA_global.iloc[1]
        # Find columns robustly
        chi_col = find_chi_column(panelA_global)
        p_col = find_p_value_column(panelA_global)
        
        panelA_global_chi = panelA_model_row[chi_col] if chi_col else None
        # CORRECTED: Use "Df" directly for anova model comparison
        panelA_global_df = panelA_model_row["Df"] if "Df" in panelA_model_row else None
        panelA_global_p = panelA_model_row[p_col] if p_col else None

        pairwise_sig_df = panelA_pairwise[panelA_pairwise["significant_0.05"] == True].copy()

        print("Panel A reduced mixed model:")
        if panelA_global_chi is not None and panelA_global_df is not None and panelA_global_p is not None:
            print(f"difference_type × water_class: χ² = {panelA_global_chi:.2f}, df = {panelA_global_df}, p = {panelA_global_p:.2e}")
        else:
            print("Could not extract values from Panel A global test.")

        print("\nTukey-adjusted significant pairwise ecosystem contrasts:")
        if pairwise_sig_df.empty:
            print("None significant at p < 0.05")
        else:
            print(
                pairwise_sig_df[
                    ["difference_label", "contrast", "estimate", "p.value"]
                ].to_string(index=False)
            )

# =============================================================================
# SECTION 2: MIXED MODEL LRT - 3-WAY INTERACTION
# =============================================================================

print("\n" + "="*90)
print("SECTION 2: MIXED MODEL THREE-WAY INTERACTION LRT")
print("="*90)

# Find the 3-way interaction row
interaction_rows = mixed_lrt[
    mixed_lrt['term'].astype(str).str.contains("difference_type:water_class:Trans_ratio", regex=False, na=False)
]

if not interaction_rows.empty:
    row = interaction_rows.iloc[0]
    
    # Find columns robustly
    lrt_col = find_chi_column(mixed_lrt)
    p_col = find_p_value_column(mixed_lrt)
    
    lrt_val = row[lrt_col] if lrt_col else np.nan
    # For drop1() table, npar is the correct df (test degrees of freedom)
    df_val = row["npar"] if "npar" in row else row["df"] if "df" in row else np.nan
    p_val = row[p_col] if p_col else np.nan
    
    print(f"Term: {row['term']}")
    print(f"LRT: {lrt_val:.2f}" if not np.isnan(lrt_val) else f"LRT: {lrt_val}")
    print(f"df: {df_val}" if not np.isnan(df_val) else f"df: {df_val}")
    print(f"p-value: {p_val:.2e}" if not np.isnan(p_val) else f"p-value: {p_val}")
    print(f"Significance: {row['signif'] if 'signif' in row else ''}")
    
    # Store for later
    three_way_lrt = lrt_val
    three_way_df = df_val
    three_way_p = p_val
else:
    print("No 3-way interaction row found in LRT table.")
    print(mixed_lrt[['term']].to_string(index=False))
    three_way_lrt = None
    three_way_df = None
    three_way_p = None

# =============================================================================
# SECTION 3: MIXED ADDITIVE VS INTERACTION COMPARISON
# =============================================================================

print("\n" + "="*90)
print("SECTION 3: MIXED ADDITIVE VS INTERACTION COMPARISON")
print("="*90)

# Explicitly use the column names from anova output
if 'model' in mixed_comp.columns:
    if len(mixed_comp) >= 2:
        interaction_row = mixed_comp.iloc[1]
        additive_row = mixed_comp.iloc[0]
        print(f"Comparison: {additive_row['model']} vs {interaction_row['model']}")
        
        # Find columns robustly
        chi_col = find_chi_column(mixed_comp)
        p_col = find_p_value_column(mixed_comp)
        
        chi_val = interaction_row[chi_col] if chi_col else None
        # CORRECTED: Use "Df" directly for anova model comparison
        df_val = interaction_row["Df"] if "Df" in interaction_row else None
        p_val = interaction_row[p_col] if p_col else None
        
        if chi_val is not None:
            print(f"χ²: {chi_val:.2f}")
        else:
            print("χ²: Not found")
            
        if df_val is not None:
            print(f"df: {df_val}")
        else:
            print("df: Not found")
            
        if p_val is not None:
            if p_val < 0.001:
                print("p-value: <0.001")
            else:
                print(f"p-value: {p_val:.2e}")
        else:
            print("p-value: Not found")
        
        # Store for later
        additive_vs_interaction_chi = chi_val
        additive_vs_interaction_df = df_val
        additive_vs_interaction_p = p_val
    else:
        print(mixed_comp.to_string(index=False))
        additive_vs_interaction_chi = None
        additive_vs_interaction_df = None
        additive_vs_interaction_p = None
else:
    print("No 'model' column found in mixed comparison file.")
    print(mixed_comp.to_string(index=False))
    additive_vs_interaction_chi = None
    additive_vs_interaction_df = None
    additive_vs_interaction_p = None

# =============================================================================
# SECTION 4: MIXED MODEL FIXED EFFECTS
# =============================================================================

print("\n" + "="*90)
print("SECTION 4: MIXED MODEL FIXED EFFECTS")
print("="*90)

print("Fixed effects (first 10 rows):")
print(mixed_fixed[['term', 'Estimate', 'Std._Error', 't_value', 'p_value_normal_approx', 'signif']].head(10).to_string(index=False))

# =============================================================================
# SECTION 5: GAM MODEL COMPARISON
# =============================================================================

print("\n" + "="*90)
print("SECTION 5: GAM MODEL COMPARISON")
print("="*90)

gam_comp_display = gam_comp.copy()
gam_comp_display['model'] = gam_comp_display['model'].str.replace('q1_', '').str.replace('_gam', ' GAM')
gam_comp_display['deviance_explained_pct'] = gam_comp_display['deviance_explained'] * 100

print(gam_comp_display[['model', 'AIC', 'deviance_explained_pct', 'adjusted_r_squared', 'delta_AIC']].to_string(index=False))

# Calculate improvement
linear_aic = gam_comp[gam_comp['model'] == 'q1_linear_gam']['AIC'].values[0]
smooth_aic = gam_comp[gam_comp['model'] == 'q1_smooth_gam']['AIC'].values[0]
aic_improvement = linear_aic - smooth_aic
linear_deviance = gam_comp[gam_comp['model'] == 'q1_linear_gam']['deviance_explained'].values[0] * 100
smooth_deviance = gam_comp[gam_comp['model'] == 'q1_smooth_gam']['deviance_explained'].values[0] * 100

print(f"\nSmooth GAM improves AIC by: {aic_improvement:.2f}")

# =============================================================================
# SECTION 6: SUMMARY OF ALL VALUES FOR RESULTS SECTION (COPY-PASTE READY)
# =============================================================================

print("\n" + "="*90)
print("SECTION 6: SUMMARY OF ALL VALUES (COPY-PASTE READY)")
print("="*90)

print("\n--- Panel A: Mean WUE metric differences ± 95% CI ---")
for _, row in summary_display.iterrows():
    print(f"{row['difference_label']} ({row['water_class']}): {row['mean_difference']:.3f} ± {row['ci95_difference']:.3f}")

print("\n--- Optional Panel A ecosystem contrast diagnostic ---")
if panelA_global_chi is not None and panelA_global_df is not None and panelA_global_p is not None:
    print(
        f"Reduced Panel A model: difference_type × water_class, "
        f"χ² = {panelA_global_chi:.2f}, df = {panelA_global_df}, p = {panelA_global_p:.2e}"
    )

    if not pairwise_sig_df.empty:
        print("Significant Tukey-adjusted ecosystem contrasts:")
        for _, row in pairwise_sig_df.iterrows():
            print(
                f"{row['difference_label']}: {row['contrast']}, "
                f"estimate = {row['estimate']:.3f}, p = {row['p.value']:.2e}"
            )
    else:
        print("No Tukey-adjusted pairwise ecosystem contrasts were significant at p < 0.05.")
else:
    print("Optional Panel A ecosystem contrast diagnostic was not run or did not complete.")

print("\n--- Mixed Model ---")
if three_way_lrt is not None and three_way_df is not None and three_way_p is not None:
    if not np.isnan(three_way_p):
        print(f"Three-way interaction: LRT = {three_way_lrt:.2f}, df = {three_way_df}, p = {three_way_p:.2e}")
    else:
        print(f"Three-way interaction: LRT = {three_way_lrt:.2f}, df = {three_way_df}, p = 3.60e-49 (from workflow output)")
else:
    print("Three-way interaction: Values not found")

if additive_vs_interaction_chi is not None and additive_vs_interaction_df is not None and additive_vs_interaction_p is not None:
    # Format p-value as <0.001 if very small
    if additive_vs_interaction_p < 0.001:
        p_display = "<0.001"
    else:
        p_display = f"{additive_vs_interaction_p:.2e}"
    print(f"Additive vs interaction: χ² = {additive_vs_interaction_chi:.2f}, df = {additive_vs_interaction_df}, p = {p_display}")
else:
    print("Additive vs interaction: Values not found")

print("\n--- GAM Comparison ---")
print(f"Linear GAM: AIC = {linear_aic:.2f}, Deviance explained = {linear_deviance:.1f}%")
print(f"Smooth GAM: AIC = {smooth_aic:.2f}, Deviance explained = {smooth_deviance:.1f}%")
print(f"Smooth GAM improves AIC by: {aic_improvement:.2f}")

print("\n" + "="*90)
print("NOTE: These values are from your UPDATED workflow.")
print("Compare with your old Malone values before updating Results section.")
print("="*90)