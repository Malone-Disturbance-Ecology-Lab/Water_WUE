# -*- coding: utf-8 -*-
"""
WUE_TRA LINEAR ANALYSIS - LOGISTIC REGRESSION (LINEAR ONLY)
Percentile-based classification (10th-90th percentile)
SPEI-48 (48-month drought)
Model: Decrease vs Increase+NoChange (NoChange kept as baseline)
For Nature paper publication

UPDATED: Weighted GLM for inference (class-balanced, same as sklearn)
FIXED: CV now standardizes within folds to avoid data leakage
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import LeaveOneGroupOut
import statsmodels.api as sm
import warnings
import os

warnings.filterwarnings('ignore')

# ==============================================================================
# OUTPUT DIRECTORY
# ==============================================================================

OUTPUT_DIR = r"M:\Research\WUE_CUE\data_products\results\linear_models"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*80)
print("WUE_TRA LOGISTIC REGRESSION - SPEI-48")
print("Percentile-based classification (10th-90th percentile)")
print("Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
print("Class-balanced weighting (GLM with freq_weights)")
print("="*80)

# ==============================================================================
# LOAD DATA
# ==============================================================================

data_path = r"M:\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
df = pd.read_csv(data_path)
print(f"\n✅ Loaded: {df.shape[0]} rows, {df.shape[1]} columns")
print(f"   Unique sites: {df['site_name'].nunique()}")

# Focus on WUE_tra, SPEI-48
metric = 'WUE_tra'
class_col = f"{metric}_minus_NNmed_SPEI_48_Class"
spei_col = 'SPEI_48'

# ==============================================================================
# MODEL: Decrease = 1, Increase+NoChange = 0 (Keep ALL rows)
# ==============================================================================

print("\n" + "="*60)
print("MODEL: Decrease vs Increase+NoChange (NoChange kept as baseline)")
print("="*60)

# Keep ALL rows (Decrease, Increase, NoChange)
mask = df[spei_col].notna()
df_model = df[mask].copy()
df_model['target'] = (df_model[class_col] == 'Decrease').astype(int)

# Class distribution
n_decrease = df_model['target'].sum()
n_increase = (df_model[class_col] == 'Increase').sum()
n_nochange = (df_model[class_col] == 'NoChange').sum()
n_total = len(df_model)
pct_decrease = n_decrease / n_total * 100

print(f"\nClass distribution (ALL classes kept - NoChange as baseline):")
print(f"   Decrease: {n_decrease} ({pct_decrease:.1f}%) → target = 1")
print(f"   Increase: {n_increase} ({n_increase/n_total*100:.1f}%) → target = 0")
print(f"   NoChange: {n_nochange} ({n_nochange/n_total*100:.1f}%) → target = 0")
print(f"   Total: {n_total} rows")

# Features (RAW values for CV, scaled for inference)
X_raw = df_model[[spei_col]].values.reshape(-1, 1)
y = df_model['target'].values
sites = df_model['site_name'].values

# ==============================================================================
# STANDARDIZE SPEI-48 (for inference only)
# ==============================================================================

print("\n" + "="*60)
print("STANDARDIZING PREDICTOR (for inference)")
print("="*60)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)

print(f"SPEI-48 raw range: [{np.min(X_raw):.2f}, {np.max(X_raw):.2f}]")
print(f"SPEI-48 mean: {np.mean(X_raw):.2f}, std: {np.std(X_raw):.2f}")
print(f"SPEI-48 standardized: mean={np.mean(X_scaled):.2e}, std={np.std(X_scaled):.2f}")

# ==============================================================================
# COMPUTE CLASS WEIGHTS (same as sklearn class_weight='balanced')
# ==============================================================================

print("\n" + "="*60)
print("COMPUTING CLASS WEIGHTS")
print("="*60)

n_pos = n_decrease
n_neg = n_total - n_decrease

# sklearn's balanced weight formula: n_samples / (n_classes * n_class_samples)
weight_pos = n_total / (2 * n_pos)
weight_neg = n_total / (2 * n_neg)

print(f"  Weight for Decrease (class=1): {weight_pos:.4f}")
print(f"  Weight for Not Decrease (class=0): {weight_neg:.4f}")

# Create weights array
sample_weights = np.where(y == 1, weight_pos, weight_neg)
print(f"  Total weight sum: {np.sum(sample_weights):.2f}")

# ==============================================================================
# WEIGHTED STATSMODELS GLM (for inference, coefficients, p-values, CI)
# ==============================================================================

print("\n" + "="*60)
print("WEIGHTED STATSMODELS GLM (class-balanced inference)")
print("="*60)

# Prepare data with constant term
X_sm = sm.add_constant(X_scaled)

# Fit weighted GLM (binomial family with freq_weights)
glm_model = sm.GLM(y, X_sm, family=sm.families.Binomial(), freq_weights=sample_weights)
glm_result = glm_model.fit()

# Extract results
coef_standardized = glm_result.params[1]
intercept = glm_result.params[0]
se = glm_result.bse[1]
p_value = glm_result.pvalues[1]
ci_lower, ci_upper = glm_result.conf_int()[1]

# Compute raw equivalent coefficient (for 1-unit SPEI change)
coef_raw_equivalent = coef_standardized / scaler.scale_[0]

# Compute odds ratios
odds_ratio_1sd_decrease = np.exp(-coef_standardized)
odds_ratio_1unit_decrease = np.exp(-coef_raw_equivalent)

print(f"\n📊 Weighted GLM results (class-balanced):")
print(f"   Coefficient (standardized, per 1 SD increase in SPEI-48): {coef_standardized:+.4f}")
print(f"   Coefficient (raw equivalent, per 1 unit SPEI increase): {coef_raw_equivalent:+.4f}")
print(f"   Intercept: {intercept:.4f}")
print(f"   SE: {se:.4f}")
print(f"   p-value: {p_value:.4e}")
print(f"   95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
print(f"   Odds ratio (per 1 SD decrease in SPEI-48): {odds_ratio_1sd_decrease:.3f}")
print(f"   Odds ratio (per 1 unit SPEI decrease): {odds_ratio_1unit_decrease:.3f}")

# ==============================================================================
# SITE-BLOCKED CROSS-VALIDATION (FIXED - NO DATA LEAKAGE)
# ==============================================================================

print("\n" + "="*60)
print("SITE-BLOCKED CROSS-VALIDATION (sklearn - model evaluation only)")
print("CV FIX: SPEI-48 standardization is now fit within each training fold")
print("to avoid data leakage.")
print("="*60)

# Use RAW SPEI values for CV (not pre-standardized)
X_raw_cv = X_raw
y_cv = y
sites_cv = sites

logo = LeaveOneGroupOut()
auc_scores = []
pr_auc_scores = []
valid_folds = 0
skipped_folds = 0

for train_idx, test_idx in logo.split(X_raw_cv, y_cv, groups=sites_cv):
    # Split raw data
    X_train_raw = X_raw_cv[train_idx]
    X_test_raw = X_raw_cv[test_idx]
    y_train = y_cv[train_idx]
    y_test = y_cv[test_idx]
    
    # Skip fold if test set has only one class
    if len(np.unique(y_test)) < 2:
        skipped_folds += 1
        continue
    
    # Standardize within fold (fit only on training data)
    scaler_fold = StandardScaler()
    X_train_scaled = scaler_fold.fit_transform(X_train_raw)
    X_test_scaled = scaler_fold.transform(X_test_raw)
    
    # Fit sklearn LogisticRegression with class_weight='balanced'
    model_fold = LogisticRegression(
        class_weight='balanced',
        solver='liblinear',
        random_state=42,
        max_iter=1000
    )
    model_fold.fit(X_train_scaled, y_train)
    
    # Predict probabilities on held-out test set
    y_pred = model_fold.predict_proba(X_test_scaled)[:, 1]
    
    # Skip if predictions are constant
    if len(np.unique(y_pred)) < 2:
        skipped_folds += 1
        continue
    
    # Calculate metrics
    auc = roc_auc_score(y_test, y_pred)
    pr_auc = average_precision_score(y_test, y_pred)
    auc_scores.append(auc)
    pr_auc_scores.append(pr_auc)
    valid_folds += 1

auc_mean = np.mean(auc_scores) if auc_scores else 0
auc_std = np.std(auc_scores) if auc_scores else 0
pr_auc_mean = np.mean(pr_auc_scores) if pr_auc_scores else 0
baseline_pr = y.mean()

print(f"\nCross-validation results (site-blocked, no data leakage):")
print(f"  Valid folds: {valid_folds}")
print(f"  Skipped folds: {skipped_folds}")
print(f"  ROC-AUC = {auc_mean:.3f} ± {auc_std:.3f}")
print(f"  PR-AUC = {pr_auc_mean:.3f} (baseline = {baseline_pr:.3f})")
print(f"  Improvement over baseline = {(pr_auc_mean - baseline_pr) / baseline_pr * 100:.0f}%")

# ==============================================================================
# PROBABILITY TABLE (from weighted GLM)
# ==============================================================================

print("\n" + "="*60)
print("PROBABILITY TABLE (from weighted GLM)")
print("="*60)

spei_values = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]
probs = []

for sv in spei_values:
    sv_scaled = (sv - scaler.mean_[0]) / scaler.scale_[0]
    X_pred = np.array([1, sv_scaled]).reshape(1, -1)
    prob = glm_result.predict(X_pred)[0]
    probs.append(prob)

print("\nSPEI-48\tP(Decrease|SPEI-48)")
print("-" * 35)
for sv, p in zip(spei_values, probs):
    print(f"{sv:.1f}\t\t{p:.3f}")

# ==============================================================================
# GENERATE PREDICTION CURVE (from weighted GLM)
# ==============================================================================

spei_range = np.linspace(-3.5, 3.5, 200)
spei_range_scaled = (spei_range - scaler.mean_[0]) / scaler.scale_[0]
X_range = np.column_stack([np.ones(len(spei_range_scaled)), spei_range_scaled])
y_pred_proba = glm_result.predict(X_range)

# ==============================================================================
# CREATE FIGURE (using weighted GLM predictions)
# ==============================================================================

print("\n" + "="*60)
print("CREATING FIGURE")
print("="*60)

fig, ax = plt.subplots(figsize=(8, 7))

# Plot S-curve from weighted GLM
ax.plot(spei_range, y_pred_proba, color='#2ca02c', linewidth=2.5, alpha=0.9, label='Logistic fit (class-balanced)')

# Plot raw predictions with jitter (using sklearn for display only)
model_sklearn = LogisticRegression(class_weight='balanced', solver='liblinear', random_state=42, max_iter=1000)
model_sklearn.fit(X_scaled, y)
y_raw_pred = model_sklearn.predict_proba(X_scaled)[:, 1]
jitter_x = np.random.normal(0, 0.03, len(X_raw))
jitter_y = np.random.normal(0, 0.01, len(y_raw_pred))
ax.scatter(X_raw.flatten() + jitter_x, y_raw_pred + jitter_y, 
           alpha=0.15, color='#2ca02c', s=10, edgecolor='none')

# Visual elements
ax.axvspan(-0.5, 0.5, alpha=0.1, color='gray', label='Near-normal (NN) region')
ax.axhline(0.5, color='black', linestyle='--', alpha=0.4, linewidth=1)

# Labels
ax.set_xlabel('SPEI-48 (Drought Index)', fontsize=13, fontweight='bold')
ax.set_ylabel('Probability of WUE$_T$ Decrease', fontsize=13, fontweight='bold')
ax.set_title('Transpirational WUE (WUE$_T$): P(Decrease) vs SPEI-48\n(NoChange kept as baseline, class-balanced)', 
             fontsize=14, fontweight='bold')

ax.set_ylim(-0.05, 1.05)
ax.set_xlim(-3.2, 3.2)
ax.grid(True, alpha=0.2, linestyle='--')
ax.tick_params(axis='both', labelsize=11)

# Annotation box (use weighted GLM p-value)
p_display = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
stats_text = f"n = {n_total:,} (all observations)\n"
stats_text += f"Decrease: {pct_decrease:.1f}%\n"
stats_text += f"Increase+NoChange: {100-pct_decrease:.1f}%\n"
stats_text += f"ROC-AUC = {auc_mean:.3f}\n"
stats_text += f"PR-AUC = {pr_auc_mean:.3f}\n"
stats_text += f"β = {coef_standardized:+.4f} ({p_display})"

ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', 
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

ax.legend(loc='lower right', fontsize=10)

plt.tight_layout()

# Save figure (SAME NAMES - no downstream impact)
png_path = os.path.join(OUTPUT_DIR, "WUE_tra_Logistic_SPEI48_10_90_Decrease_vs_Baseline.png")
pdf_path = os.path.join(OUTPUT_DIR, "WUE_tra_Logistic_SPEI48_10_90_Decrease_vs_Baseline.pdf")
fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"✅ Saved: {png_path}")
print(f"✅ Saved: {pdf_path}")

# ==============================================================================
# SAVE SUMMARY CSV (from weighted GLM - SAME FILENAME, NO CHANGE)
# ==============================================================================

print("\n" + "="*60)
print("SAVING SUMMARY (from weighted GLM - class-balanced)")
print("="*60)

summary_path = os.path.join(OUTPUT_DIR, "WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv")

# SPEI range columns
spei_min_training = float(X_raw.min())
spei_max_training = float(X_raw.max())
spei_mean_training = float(X_raw.mean())
spei_std_training = float(X_raw.std())

summary_data = {
    # Existing columns (kept for downstream compatibility)
    'Metric': metric,
    'SPEI_Timescale': 48,
    'Percentile_Range': '10-90',
    'Model_Type': 'Decrease_vs_Increase_NoChange',
    'Coefficient': coef_standardized,
    'Intercept': intercept,
    'Odds_Ratio_1unit_decrease': odds_ratio_1unit_decrease,
    'ROC_AUC': auc_mean,
    'ROC_AUC_Std': auc_std,
    'PR_AUC': pr_auc_mean,
    'Baseline_PR': baseline_pr,
    'N_Total': n_total,
    'N_Decrease': n_decrease,
    'N_Increase': n_increase,
    'N_NoChange': n_nochange,
    'Pct_Decrease': pct_decrease,
    
    # Inference columns (from weighted GLM - NOW CORRECT!)
    'Coefficient_SE': se,
    'Coefficient_P_value': p_value,
    'Coefficient_CI_lower': ci_lower,
    'Coefficient_CI_upper': ci_upper,
    'Inference_Method': 'weighted statsmodels GLM (class-balanced)',
    
    # SPEI range for spatial model
    'spei_min': spei_min_training,
    'spei_max': spei_max_training,
    'spei_mean': spei_mean_training,
    'spei_std': spei_std_training,
    
    # Additional clarity columns
    'Coefficient_Standardized': coef_standardized,
    'Coefficient_RawEquivalent': coef_raw_equivalent,
    'Odds_Ratio_1SD_decrease': odds_ratio_1sd_decrease,
    'Weight_Decrease': weight_pos,
    'Weight_NotDecrease': weight_neg,
    'Model_Family': 'weighted GLM binomial',
    'Class_Weighting': 'balanced',
    'CV_Method': 'site-blocked LeaveOneGroupOut (no leakage: scaling within folds)',
    'Valid_Folds': valid_folds,
    'Skipped_Folds': skipped_folds
}

summary_df = pd.DataFrame([summary_data])
summary_df.to_csv(summary_path, index=False)
print(f"✅ Saved summary: {summary_path}")
print(f"   Weight_Decrease={weight_pos:.4f}, Weight_NotDecrease={weight_neg:.4f}")
print(f"   Inference now matches class-balanced model")
print(f"   CV method: scaling within folds (no data leakage)")

# ==============================================================================
# SAVE PROBABILITY TABLE (from weighted GLM - SAME FILENAME, NO CHANGE)
# ==============================================================================

print("\n" + "="*60)
print("SAVING PROBABILITY TABLE (from weighted GLM)")
print("="*60)

prob_df = pd.DataFrame({
    'SPEI': spei_values,
    'P_Decrease_SPEI48': probs
})
prob_path = os.path.join(OUTPUT_DIR, "WUE_tra_Probability_Table_10_90.csv")
prob_df.to_csv(prob_path, index=False)
print(f"✅ Saved probability table: {prob_path}")

# ==============================================================================
# FINAL RESULTS SUMMARY
# ==============================================================================

print("\n" + "="*60)
print("STEP 2 COMPLETE - WEIGHTED GLM (CLASS-BALANCED)")
print("="*60)
print(f"\n📊 Final results for {metric} (10-90 percentile, NoChange kept as baseline):")
print(f"   Class weighting: Decrease={weight_pos:.4f}, Not Decrease={weight_neg:.4f}")
print(f"   Coefficient (standardized): {coef_standardized:+.4f}")
print(f"   Coefficient (raw equivalent): {coef_raw_equivalent:+.4f}")
print(f"   SE: {se:.4f}")
print(f"   95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
print(f"   p-value: {p_value:.2e}")
print(f"   Odds ratio (per 1 SD decrease): {odds_ratio_1sd_decrease:.3f}")
print(f"   Odds ratio (per 1 unit decrease): {odds_ratio_1unit_decrease:.3f}")
print(f"   ROC-AUC = {auc_mean:.3f} ± {auc_std:.3f}")
print(f"   PR-AUC = {pr_auc_mean:.3f} (baseline = {baseline_pr:.3f})")
print(f"   N (total) = {n_total}")

print(f"\n   SPEI training range: [{spei_min_training:.3f}, {spei_max_training:.3f}]")
print(f"   SPEI training mean ± std: {spei_mean_training:.3f} ± {spei_std_training:.3f}")
print(f"\n   Model: weighted GLM (class-balanced, freq_weights)")
print(f"   Inference: standard SE (not cluster-robust)")
print(f"   ✓ p-value and CI are from the same class-balanced model")
print(f"   ✓ CV: site-blocked with scaling within folds (no data leakage)")

plt.show()