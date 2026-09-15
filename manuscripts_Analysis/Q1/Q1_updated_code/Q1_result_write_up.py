"""
Q1 RESULTS EXTRACTION SCRIPT - COMPLETE (UPDATED for single-difference workflow)
Reads existing updated CSV outputs and prints exact values from actual results.

UPDATED: Only WUE_ET - WUE_T is modeled/reported.
  - Removed all three-difference-type logic.
  - Mixed model focus: water_class * Trans_ratio_z_NN interaction.
  - Added extraction of 3 ecosystem-specific GAM smooth terms (edf/F/p).
  - Added overall smooth GAM stats (adjusted R², deviance, AIC).
  - Added predicted WUE_ET - WUE_T values at the lower end of the modeled T:ET
    gradient (5th percentile = minimum of the GAM prediction range).
  - Removed the optional Panel A difference_type x ecosystem emmeans block.

Inputs (all from Q1_updated_results/outputs/):
- Q1_final_near_normal_summary_by_difference_type.csv
- Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv
- Q1_final_near_normal_mixed_likelihood_ratio_tests.csv
- Q1_final_near_normal_mixed_model_comparison_lrt.csv
- Q1_final_gam_model_comparison.csv
- Q1_final_smooth_gam_smooth_terms.csv
- Q1_final_smooth_gam_predictions_TET.csv
- Q1_final_near_normal_mixed_predictions_TET.csv

Output: Printed console report with exact values from actual results
"""

import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

output_dir = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs")

# Low T:ET point used in Results narrative.
# None  -> use the minimum of the GAM prediction range, which corresponds to
#          the 5th-percentile T:ET (consistent with previous Results reporting).
# number -> use that specific T:ET value.
LOW_TET_POINT = None

print("="*90)
print("Q1 UPDATED RESULTS - EXACT VALUES FROM ACTUAL RESULTS (WUE_ET - WUE_T ONLY)")
print("="*90)
print(f"\nReading from: {output_dir}\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def find_p_value_column(df):
    """Find the p-value column robustly."""
    for c in df.columns:
        if c.lower() in ['pr(>chisq)', 'pr(>chi)', 'pr(>f)', 'p_value', 'pvalue', 'p']:
            return c
    for c in df.columns:
        if c.startswith("Pr") or c.startswith("p_"):
            return c
    return None

def find_chi_column(df):
    """Find the Chi-square / LRT column robustly."""
    for c in df.columns:
        if c in ("Chisq", "LRT", "Chi", "LR"):
            return c
    for c in df.columns:
        if c.lower() in ("chisq", "lrt", "chi", "lr"):
            return c
    return None

def find_stat_column(df):
    """Find a t-value or z-value column."""
    for c in df.columns:
        if c.lower() in ('t_value', 't value', 'z_value', 'z value'):
            return c
    return None

def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    if p < 0.001:
        return f"{p:.2e} (<0.001)"
    return f"{p:.4f}"

# =============================================================================
# READ CSV FILES
# =============================================================================

print("Loading CSV outputs...")

summary_df   = pd.read_csv(output_dir / "Q1_final_near_normal_summary_by_difference_type.csv")
mixed_fixed  = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_fixed_effects_with_approx_p.csv")
mixed_lrt    = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv")
mixed_comp   = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_model_comparison_lrt.csv")
gam_comp     = pd.read_csv(output_dir / "Q1_final_gam_model_comparison.csv")
smooth_terms = pd.read_csv(output_dir / "Q1_final_smooth_gam_smooth_terms.csv")
smooth_pred  = pd.read_csv(output_dir / "Q1_final_smooth_gam_predictions_TET.csv")
mixed_pred   = pd.read_csv(output_dir / "Q1_final_near_normal_mixed_predictions_TET.csv")

# =============================================================================
# SECTION 1: PANEL A MEANS (single difference, three ecosystems)
# =============================================================================

print("\n" + "="*90)
print("SECTION 1: PANEL A MEANS (mean ± 95% CI) - WUE_ET - WUE_T by ecosystem")
print("="*90)

summary_display = summary_df.copy()
summary_display["mean_ci"] = summary_display.apply(
    lambda x: f"{x['mean_difference']:.3f} ± {x['ci95_difference']:.3f}", axis=1
)

print(summary_display[["difference_label", "water_class", "n_months", "n_sites", "mean_ci"]].to_string(index=False))

# =============================================================================
# SECTION 2: MIXED MODEL - water_class x T:ET INTERACTION LRT
# =============================================================================

print("\n" + "="*90)
print("SECTION 2: MIXED MODEL - water_class × Trans_ratio_z_NN INTERACTION LRT")
print("="*90)

# Look for the interaction row in the LRT table. drop1() output uses the model
# formula term names, so the interaction row contains both "water_class" and
# "Trans_ratio" (and a colon).
interaction_rows = mixed_lrt[
    mixed_lrt['term'].astype(str).str.contains("water_class", na=False) &
    mixed_lrt['term'].astype(str).str.contains("Trans_ratio", na=False) &
    mixed_lrt['term'].astype(str).str.contains(":", na=False, regex=False)
]

lrt_col = find_chi_column(mixed_lrt)
p_col   = find_p_value_column(mixed_lrt)

if not interaction_rows.empty:
    row = interaction_rows.iloc[0]
    lrt_val = row[lrt_col] if lrt_col else np.nan
    df_val  = row["npar"] if "npar" in row else (row["Df"] if "Df" in row else np.nan)
    p_val   = row[p_col] if p_col else np.nan

    print(f"Term: {row['term']}")
    print(f"LRT (χ²): {lrt_val:.2f}" if not np.isnan(lrt_val) else f"LRT: {lrt_val}")
    print(f"df: {df_val}")
    print(f"p-value: {fmt_p(p_val)}")
    print(f"Significance: {row['signif'] if 'signif' in row else ''}")

    interaction_lrt = lrt_val
    interaction_df  = df_val
    interaction_p   = p_val
else:
    print("No water_class × Trans_ratio interaction row found.")
    print("Available terms in mixed_lrt:")
    print(mixed_lrt[['term']].to_string(index=False))
    interaction_lrt = None
    interaction_df  = None
    interaction_p   = None

# =============================================================================
# SECTION 3: MIXED ADDITIVE VS INTERACTION COMPARISON
# =============================================================================

print("\n" + "="*90)
print("SECTION 3: MIXED ADDITIVE VS INTERACTION COMPARISON")
print("="*90)

if 'model' in mixed_comp.columns and len(mixed_comp) >= 2:
    interaction_row = mixed_comp.iloc[1]
    additive_row    = mixed_comp.iloc[0]
    print(f"Comparison: {additive_row['model']} vs {interaction_row['model']}")

    chi_col = find_chi_column(mixed_comp)
    p_col_c = find_p_value_column(mixed_comp)

    chi_val = interaction_row[chi_col] if chi_col else None
    df_val  = interaction_row["Df"] if "Df" in interaction_row else None
    p_val   = interaction_row[p_col_c] if p_col_c else None

    print(f"χ²: {chi_val:.2f}" if chi_val is not None else "χ²: Not found")
    print(f"df: {df_val}" if df_val is not None else "df: Not found")
    print(f"p-value: {fmt_p(p_val)}" if p_val is not None else "p-value: Not found")

    additive_vs_interaction_chi = chi_val
    additive_vs_interaction_df  = df_val
    additive_vs_interaction_p   = p_val
else:
    print(mixed_comp.to_string(index=False))
    additive_vs_interaction_chi = None
    additive_vs_interaction_df  = None
    additive_vs_interaction_p   = None

# =============================================================================
# SECTION 4: MIXED MODEL FIXED EFFECTS
# =============================================================================

print("\n" + "="*90)
print("SECTION 4: MIXED MODEL FIXED EFFECTS")
print("="*90)

est_col = next((c for c in mixed_fixed.columns if c.lower() in ("estimate",)), "Estimate")
se_col  = next((c for c in mixed_fixed.columns if c.lower() in ("std._error", "std_error", "stderr", "se")), None)
t_col   = find_stat_column(mixed_fixed)
p_col_f = "p_value_normal_approx" if "p_value_normal_approx" in mixed_fixed.columns else find_p_value_column(mixed_fixed)

show_cols = ['term', est_col]
if se_col:  show_cols.append(se_col)
if t_col:   show_cols.append(t_col)
if p_col_f: show_cols.append(p_col_f)
if 'signif' in mixed_fixed.columns: show_cols.append('signif')

print(mixed_fixed[show_cols].to_string(index=False))

# =============================================================================
# SECTION 5: GAM MODEL COMPARISON (linear vs smooth)
# =============================================================================

print("\n" + "="*90)
print("SECTION 5: GAM MODEL COMPARISON (linear vs smooth)")
print("="*90)

gam_comp_display = gam_comp.copy()
gam_comp_display['model'] = (
    gam_comp_display['model'].str.replace('q1_', '').str.replace('_gam', ' GAM')
)
gam_comp_display['deviance_explained_pct'] = gam_comp_display['deviance_explained'] * 100

print(gam_comp_display[['model', 'AIC', 'deviance_explained_pct',
                        'adjusted_r_squared', 'delta_AIC']].to_string(index=False))

linear_aic = gam_comp[gam_comp['model'] == 'q1_linear_gam']['AIC'].values[0]
smooth_aic = gam_comp[gam_comp['model'] == 'q1_smooth_gam']['AIC'].values[0]
aic_improvement = linear_aic - smooth_aic
linear_deviance = gam_comp[gam_comp['model'] == 'q1_linear_gam']['deviance_explained'].values[0] * 100
smooth_deviance = gam_comp[gam_comp['model'] == 'q1_smooth_gam']['deviance_explained'].values[0] * 100
smooth_r2 = gam_comp[gam_comp['model'] == 'q1_smooth_gam']['adjusted_r_squared'].values[0]

print(f"\nSmooth GAM improves AIC by: {aic_improvement:.2f}")

# =============================================================================
# SECTION 6: SMOOTH GAM TERMS (3 ecosystem-specific T:ET smooths)
# =============================================================================

print("\n" + "="*90)
print("SECTION 6: SMOOTH GAM TERMS (ecosystem-specific T:ET smooths)")
print("="*90)

edf_col = next((c for c in smooth_terms.columns if c.lower() == 'edf'), None)
f_col   = next((c for c in smooth_terms.columns if c.strip().lower() in ('f', 'f statistic', 'f_statistic')), None)
p_col_s = next((c for c in smooth_terms.columns if 'p' in c.lower() and 'value' in c.lower()), None)
if p_col_s is None:
    p_col_s = next((c for c in smooth_terms.columns if c.lower().startswith('p')), None)

ecosystem_smooths = smooth_terms[
    smooth_terms['term'].astype(str).str.contains("diff_ecosystem", na=False) &
    smooth_terms['term'].astype(str).str.contains("Trans_ratio", na=False)
].copy()

print(f"Number of ecosystem-specific T:ET smooth terms: {len(ecosystem_smooths)}")
print()

def extract_ecosystem(term):
    # term looks like: s(Trans_ratio):diff_ecosystemWUE_ET_minus_WUE_T__Upland
    m = re.search(r"diff_ecosystem.+?__(.+)$", str(term))
    return m.group(1) if m else "?"

print(f"{'Ecosystem':<15}{'edf':>8}{'F':>12}{'p-value':>16}{'signif':>10}")
print("-" * 61)
for _, row in ecosystem_smooths.iterrows():
    eco = extract_ecosystem(row['term'])
    edf = row[edf_col] if edf_col else np.nan
    fv  = row[f_col]   if f_col   else np.nan
    pv  = row[p_col_s] if p_col_s else np.nan
    sig = row['signif'] if 'signif' in row else ''
    print(f"{eco:<15}{edf:>8.3f}{fv:>12.3f}{fmt_p(pv):>16}{sig:>10}")

other_smooths = smooth_terms[
    smooth_terms['term'].astype(str).str.contains("site_name|month_f", na=False)
]
if not other_smooths.empty:
    print("\nRandom-effect smooths:")
    for _, row in other_smooths.iterrows():
        edf = row[edf_col] if edf_col else np.nan
        pv  = row[p_col_s] if p_col_s else np.nan
        print(f"  {row['term']}: edf = {edf:.3f}, p = {fmt_p(pv)}")

# =============================================================================
# SECTION 7: PREDICTED WUE_ET - WUE_T AT LOW T:ET
# =============================================================================

print("\n" + "="*90)
print("SECTION 7: PREDICTED WUE_ET - WUE_T AT LOW T:ET")
print("="*90)

smooth_pred_sorted = smooth_pred.copy()

if LOW_TET_POINT is None:
    target_tet = float(smooth_pred_sorted['Trans_ratio'].min())
    print(f"LOW_TET_POINT not set; using minimum of GAM prediction range "
          f"(5th-percentile T:ET): {target_tet:.4f}")
else:
    target_tet = float(LOW_TET_POINT)
    print(f"Target low T:ET point: {target_tet:.4f}")

print(f"\nSmooth GAM predictions at T:ET = {target_tet:.4f}:")
for eco in ["Upland", "Freshwater", "Saline"]:
    sub = smooth_pred_sorted[smooth_pred_sorted['ecosystem_label'] == eco]
    if sub.empty:
        continue
    idx = (sub['Trans_ratio'] - target_tet).abs().idxmin()
    row = sub.loc[idx]
    print(f"  {eco:<12} T:ET={row['Trans_ratio']:.4f}  "
          f"pred = {row['prediction']:.3f}  "
          f"95% CI [{row['lower']:.3f}, {row['upper']:.3f}]")

print(f"\nMixed model predictions at T:ET = {target_tet:.4f}:")
for eco in ["Upland", "Freshwater", "Saline"]:
    sub = mixed_pred[mixed_pred['ecosystem_label'] == eco]
    if sub.empty:
        continue
    idx = (sub['Trans_ratio'] - target_tet).abs().idxmin()
    row = sub.loc[idx]
    print(f"  {eco:<12} T:ET={row['Trans_ratio']:.4f}  "
          f"pred = {row['prediction']:.3f}  "
          f"95% CI [{row['lower']:.3f}, {row['upper']:.3f}]")

# =============================================================================
# SECTION 8: SUMMARY OF ALL VALUES (COPY-PASTE READY)
# =============================================================================

print("\n" + "="*90)
print("SECTION 8: SUMMARY OF ALL VALUES (COPY-PASTE READY)")
print("="*90)

print("\n--- Panel A: Mean WUE_ET - WUE_T ± 95% CI by ecosystem ---")
for _, row in summary_display.iterrows():
    print(f"{row['water_class']}: {row['mean_difference']:.3f} ± {row['ci95_difference']:.3f} "
          f"(n = {row['n_months']} site-months, {row['n_sites']} sites)")

print("\n--- Mixed Model: water_class × T:ET interaction ---")
if interaction_lrt is not None and interaction_df is not None and interaction_p is not None:
    print(f"water_class × Trans_ratio: χ² = {interaction_lrt:.2f}, "
          f"df = {interaction_df}, p = {interaction_p:.2e}")
else:
    print("water_class × Trans_ratio: Values not found")

print("\n--- Additive vs interaction model comparison ---")
if additive_vs_interaction_chi is not None:
    print(f"χ² = {additive_vs_interaction_chi:.2f}, "
          f"df = {additive_vs_interaction_df}, "
          f"p = {fmt_p(additive_vs_interaction_p)}")
else:
    print("Values not found")

print("\n--- GAM Comparison ---")
print(f"Linear GAM: AIC = {linear_aic:.2f}, Deviance explained = {linear_deviance:.1f}%")
print(f"Smooth GAM: AIC = {smooth_aic:.2f}, Deviance explained = {smooth_deviance:.1f}%, "
      f"adjusted R² = {smooth_r2:.3f}")
print(f"Smooth GAM improves AIC by: {aic_improvement:.2f}")

print("\n--- Smooth GAM ecosystem-specific T:ET smooths ---")
for _, row in ecosystem_smooths.iterrows():
    eco = extract_ecosystem(row['term'])
    edf = row[edf_col] if edf_col else np.nan
    fv  = row[f_col]   if f_col   else np.nan
    pv  = row[p_col_s] if p_col_s else np.nan
    print(f"{eco}: edf = {edf:.3f}, F = {fv:.3f}, p = {fmt_p(pv)}")

print("\n" + "="*90)
print("NOTE: These values are from your UPDATED workflow (WUE_ET - WUE_T only).")
print("The WUE_T/T:ET coast and ecosystem sentences in §3.2 come from the separate")
print("WUE_T/T:ET analysis and are NOT produced by this script.")
print("="*90)