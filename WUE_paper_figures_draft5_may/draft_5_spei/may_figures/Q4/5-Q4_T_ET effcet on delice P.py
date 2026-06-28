# -*- coding: utf-8 -*-
"""
Q4 ANALYSIS: WUE_T DECLINE ~ SPEI-48 × T:ET RATIO INTERACTION
Test whether drought effect on WUE_T decline depends on ecosystem T:ET ratio

UPDATED: Now uses weighted statsmodels GLM with class-balanced weights
(matches Step 2 point-model approach)
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, average_precision_score
import warnings
import os

warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DATA_PATH = r"M:\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
OUTPUT_DIR = r"M:\Research\WUE_CUE\data_products\results\linear_models"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================

print("="*80)
print("Q4 ANALYSIS: WUE_T DECLINE ~ SPEI-48 × T:ET RATIO INTERACTION")
print("Weighted GLM (class-balanced) - Matches Step 2 workflow")
print("="*80)

df = pd.read_csv(DATA_PATH)
print(f"\n✅ Loaded final dataset: {df.shape[0]} rows, {df.shape[1]} columns")
print(f"   Unique sites: {df['site_name'].nunique()}")

# ==============================================================================
# STEP 2: CHECK REQUIRED COLUMNS
# ==============================================================================

print("\n" + "="*60)
print("STEP 2: CHECKING REQUIRED COLUMNS")
print("="*60)

# Define columns
METRIC = 'WUE_tra'
SPEI_COL = 'SPEI_48'
CLASS_COL = f"{METRIC}_minus_NNmed_{SPEI_COL}_Class"
TET_COL = 'Trans_ratio'

# Check if all columns exist
required_cols = [SPEI_COL, CLASS_COL, TET_COL, 'site_name']
missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    print(f"❌ Missing columns: {missing_cols}")
    print("\nAvailable columns in final dataset:")
    for i, col in enumerate(df.columns, 1):
        print(f"   {i:3d}. {col}")
    raise ValueError("Cannot proceed without required columns")

print(f"✅ All required columns found:")
print(f"   SPEI-48: {SPEI_COL}")
print(f"   WUE_T class: {CLASS_COL}")
print(f"   T:ET ratio: {TET_COL}")

# ==============================================================================
# STEP 3: PREPARE DATA FOR MODEL
# ==============================================================================

print("\n" + "="*60)
print("STEP 3: PREPARING DATA")
print("="*60)

# Filter to valid rows (SPEI_48, Trans_ratio, and class label not missing)
mask = (df[SPEI_COL].notna() & 
        df[TET_COL].notna() & 
        df[CLASS_COL].notna())

df_model = df[mask].copy()

# Keep only valid class labels (Decrease, Increase, NoChange)
valid_classes = ['Decrease', 'Increase', 'NoChange']
df_model = df_model[df_model[CLASS_COL].isin(valid_classes)].copy()

# Create binary target: Decrease = 1, Not decrease (Increase/NoChange) = 0
df_model['target'] = (df_model[CLASS_COL] == 'Decrease').astype(int)

print(f"\nObservations after filtering: {len(df_model):,}")
print(f"Unique sites: {df_model['site_name'].nunique()}")

# Class distribution
n_decrease = df_model['target'].sum()
n_not_decrease = len(df_model) - n_decrease
n_total = len(df_model)
pct_decrease = 100 * n_decrease / n_total

print(f"\nClass distribution:")
print(f"  Decrease (1): {n_decrease:,} ({pct_decrease:.1f}%)")
print(f"  Not decrease (0): {n_not_decrease:,} ({100-pct_decrease:.1f}%)")

# T:ET statistics
tet_values = df_model[TET_COL].values
tet_min, tet_max = tet_values.min(), tet_values.max()
tet_mean, tet_median = tet_values.mean(), np.median(tet_values)

print(f"\nT:ET ratio (Trans_ratio) statistics:")
print(f"  Range: [{tet_min:.4f}, {tet_max:.4f}]")
print(f"  Mean: {tet_mean:.4f}")
print(f"  Median: {tet_median:.4f}")

# SPEI-48 statistics
spei_values = df_model[SPEI_COL].values
spei_min, spei_max = spei_values.min(), spei_values.max()
spei_mean, spei_median = spei_values.mean(), np.median(spei_values)
spei_std = spei_values.std()

print(f"\nSPEI-48 statistics:")
print(f"  Range: [{spei_min:.2f}, {spei_max:.2f}]")
print(f"  Mean: {spei_mean:.2f}")
print(f"  Std: {spei_std:.2f}")

# ==============================================================================
# STEP 4: STANDARDIZE PREDICTORS
# ==============================================================================

print("\n" + "="*60)
print("STEP 4: STANDARDIZING PREDICTORS")
print("="*60)

# Standardize SPEI-48 and T:ET
scaler_spei = StandardScaler()
scaler_tet = StandardScaler()

SPEI_z = scaler_spei.fit_transform(df_model[[SPEI_COL]]).flatten()
TET_z = scaler_tet.fit_transform(df_model[[TET_COL]]).flatten()

# Create interaction term
SPEI_TET_interaction = SPEI_z * TET_z

print(f"SPEI-48 standardized: mean={SPEI_z.mean():.2e}, std={SPEI_z.std():.2f}")
print(f"T:ET standardized: mean={TET_z.mean():.2e}, std={TET_z.std():.2f}")
print(f"Interaction created: SPEI_z × TET_z")

# ==============================================================================
# STEP 5: COMPUTE CLASS WEIGHTS (sklearn balanced formula)
# ==============================================================================

print("\n" + "="*60)
print("STEP 5: COMPUTING CLASS WEIGHTS")
print("="*60)

# sklearn's balanced weight formula: n_samples / (n_classes * n_class_samples)
weight_pos = n_total / (2 * n_decrease)
weight_neg = n_total / (2 * n_not_decrease)

print(f"  Weight for Decrease (class=1): {weight_pos:.4f}")
print(f"  Weight for Not Decrease (class=0): {weight_neg:.4f}")

# Create weights array
sample_weights = np.where(df_model['target'].values == 1, weight_pos, weight_neg)
print(f"  Total weight sum: {np.sum(sample_weights):.2f}")

# ==============================================================================
# STEP 6: FIT WEIGHTED STATSMODELS GLM (BINOMIAL WITH FREQ_WEIGHTS)
# ==============================================================================

print("\n" + "="*60)
print("STEP 6: FIT WEIGHTED GLM (class-balanced inference)")
print("="*60)

# Create design matrix
X = np.column_stack([SPEI_z, TET_z, SPEI_TET_interaction])
X = sm.add_constant(X)

y = df_model['target'].values

# Fit weighted GLM (binomial family with freq_weights)
glm_model = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=sample_weights)
glm_result = glm_model.fit()

print("\n✅ Weighted GLM fitted successfully")
print(f"  Log-Likelihood: {glm_result.llf:.2f}")
print(f"  Deviance: {glm_result.deviance:.2f}")

# ==============================================================================
# STEP 7: EXTRACT COEFFICIENTS WITH 95% CI AND ODDS RATIOS
# ==============================================================================

print("\n" + "="*60)
print("STEP 7: COEFFICIENT TABLE (weighted GLM)")
print("="*60)

# Get coefficients, standard errors, p-values, and confidence intervals
params = glm_result.params
conf_int = glm_result.conf_int()
p_values = glm_result.pvalues
se = glm_result.bse

# Calculate odds ratios (note: for SPEI effect, OR per 1 SD decrease is exp(-coef))
# But we report standard OR = exp(coef) for interpretation
odds_ratios = np.exp(params)

# Create clean table
var_names = ['const', 'SPEI48_z', 'TET_z', 'SPEI48_z:TET_z']
coef_table = []

for i, var in enumerate(var_names):
    coef_table.append({
        'Variable': var,
        'Coefficient': params[i],
        'SE': se[i],
        'CI_lower': conf_int[i, 0],
        'CI_upper': conf_int[i, 1],
        'p_value': p_values[i],
        'Odds_Ratio': odds_ratios[i]
    })

coef_df = pd.DataFrame(coef_table)
print("\n", coef_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

# ==============================================================================
# STEP 8: INTERPRETATION (UPDATED - CORRECT DIRECTION)
# ==============================================================================

print("\n" + "="*60)
print("STEP 8: INTERPRETATION")
print("="*60)

# SPEI-48 effect
spei_p = p_values[1]
spei_coef = params[1]
if spei_p < 0.05:
    if spei_coef < 0:
        print("\n✅ SPEI-48 effect: SIGNIFICANT NEGATIVE")
        print("   → Lower SPEI (drought) INCREASES probability of WUE_T decline")
    else:
        print("\n✅ SPEI-48 effect: SIGNIFICANT POSITIVE")
        print("   → Higher SPEI (wetter) INCREASES probability of WUE_T decline")
else:
    print("\n⚠️ SPEI-48 effect: NOT SIGNIFICANT")
    print("   → No evidence that SPEI-48 affects WUE_T decline probability")

# T:ET main effect
tet_p = p_values[2]
if tet_p < 0.05:
    print(f"\n✅ T:ET main effect: SIGNIFICANT (p = {tet_p:.4f})")
else:
    print(f"\n⚠️ T:ET main effect: NOT SIGNIFICANT (p = {tet_p:.4f})")

# Interaction effect (UPDATED CORRECT INTERPRETATION)
interaction_p = p_values[3]
interaction_coef = params[3]
interaction_se = se[3]

if interaction_p < 0.05:
    print(f"\n✅ SPEI-48 × T:ET INTERACTION: SIGNIFICANT (p = {interaction_p:.4f})")
    print("   → The drought effect DIFFERS across the T:ET gradient")
    if interaction_coef < 0:
        print("   → NEGATIVE INTERACTION: SPEI-48 slope becomes MORE NEGATIVE at higher T:ET")
        print("   → Because lower SPEI represents drier conditions, this suggests")
        print("     STRONGER drought-associated WUE_T decline risk at higher T:ET")
    else:
        print("   → POSITIVE INTERACTION: SPEI-48 slope becomes LESS NEGATIVE at higher T:ET")
        print("   → Because lower SPEI represents drier conditions, this suggests")
        print("     WEAKER drought-associated WUE_T decline risk at higher T:ET")
else:
    print(f"\n⚠️ SPEI-48 × T:ET INTERACTION: NOT SIGNIFICANT (p = {interaction_p:.4f})")
    print("   → No evidence that drought effect depends on T:ET ratio")

# ==============================================================================
# STEP 9: PREDICTED PROBABILITIES (class-balanced fitted probabilities)
# ==============================================================================

print("\n" + "="*60)
print("STEP 9: PREDICTED PROBABILITIES BY T:ET GROUP")
print("Note: These are class-balanced fitted probabilities, not raw prevalence")
print("="*60)

# Define low and high T:ET based on percentiles
tet_25th = np.percentile(df_model[TET_COL], 25)
tet_75th = np.percentile(df_model[TET_COL], 75)

print(f"\nT:ET thresholds:")
print(f"  Low T:ET (25th percentile): {tet_25th:.4f}")
print(f"  High T:ET (75th percentile): {tet_75th:.4f}")

# Standardize these values
tet_low_z = (tet_25th - scaler_tet.mean_[0]) / scaler_tet.scale_[0]
tet_high_z = (tet_75th - scaler_tet.mean_[0]) / scaler_tet.scale_[0]

# SPEI values to predict
spei_values_pred = [1.0, 0.0, -1.0, -2.0, -3.0]

# Calculate predicted probabilities from weighted GLM
prob_table = []
for spei in spei_values_pred:
    spei_z = (spei - scaler_spei.mean_[0]) / scaler_spei.scale_[0]
    
    # Low T:ET
    interaction_low = spei_z * tet_low_z
    X_pred_low = np.array([1, spei_z, tet_low_z, interaction_low])
    log_odds_low = np.dot(X_pred_low, params)
    prob_low = 1 / (1 + np.exp(-log_odds_low))
    
    # High T:ET
    interaction_high = spei_z * tet_high_z
    X_pred_high = np.array([1, spei_z, tet_high_z, interaction_high])
    log_odds_high = np.dot(X_pred_high, params)
    prob_high = 1 / (1 + np.exp(-log_odds_high))
    
    prob_table.append({
        'SPEI-48': spei,
        'Class-balanced P(Decrease) | Low T:ET': prob_low,
        'Class-balanced P(Decrease) | High T:ET': prob_high,
        'Difference (Low - High)': prob_low - prob_high
    })

prob_df = pd.DataFrame(prob_table)
print("\n", prob_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

# ==============================================================================
# STEP 10: SITE-BLOCKED CROSS-VALIDATION (NO DATA LEAKAGE)
# ==============================================================================

print("\n" + "="*60)
print("STEP 10: SITE-BLOCKED CROSS-VALIDATION")
print("CV FIX: Standardization and SPEI × T:ET interaction are now computed")
print("within each training fold to avoid data leakage.")
print("Using sklearn LogisticRegression(class_weight='balanced') for evaluation only")
print("="*60)

# Use RAW predictors (not pre-standardized)
X_raw = df_model[[SPEI_COL, TET_COL]].values
y_raw = df_model['target'].values
sites = df_model['site_name'].values

logo = LeaveOneGroupOut()
auc_scores = []
pr_auc_scores = []
valid_folds = 0
skipped_folds = 0

for train_idx, test_idx in logo.split(X_raw, y_raw, groups=sites):
    # Split raw data
    X_train_raw = X_raw[train_idx]
    X_test_raw = X_raw[test_idx]
    y_train = y_raw[train_idx]
    y_test = y_raw[test_idx]
    
    # Skip fold if test set has only one class
    if len(np.unique(y_test)) < 2:
        skipped_folds += 1
        continue
    
    # Standardize within fold (fit only on training data)
    scaler_fold = StandardScaler()
    X_train_scaled = scaler_fold.fit_transform(X_train_raw)
    X_test_scaled = scaler_fold.transform(X_test_raw)
    
    # Create interaction term AFTER scaling
    train_interaction = X_train_scaled[:, 0] * X_train_scaled[:, 1]
    test_interaction = X_test_scaled[:, 0] * X_test_scaled[:, 1]
    
    # Build design matrices
    X_train_model = np.column_stack([
        X_train_scaled[:, 0],  # SPEI_z
        X_train_scaled[:, 1],  # TET_z
        train_interaction       # SPEI_z × TET_z
    ])
    
    X_test_model = np.column_stack([
        X_test_scaled[:, 0],
        X_test_scaled[:, 1],
        test_interaction
    ])
    
    # Fit sklearn LogisticRegression with class_weight='balanced'
    model_fold = LogisticRegression(
        class_weight='balanced',
        solver='liblinear',
        random_state=42,
        max_iter=1000
    )
    model_fold.fit(X_train_model, y_train)
    
    # Predict probabilities on held-out test set
    y_pred = model_fold.predict_proba(X_test_model)[:, 1]
    
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

# Calculate summary statistics
auc_mean = np.mean(auc_scores) if auc_scores else 0
auc_std = np.std(auc_scores) if auc_scores else 0
pr_auc_mean = np.mean(pr_auc_scores) if pr_auc_scores else 0
baseline_pr = y_raw.mean()

print(f"\nCross-validation results (site-blocked, no data leakage):")
print(f"  Valid folds: {valid_folds}")
print(f"  Skipped folds: {skipped_folds}")
print(f"  ROC-AUC = {auc_mean:.3f} ± {auc_std:.3f}")
print(f"  PR-AUC = {pr_auc_mean:.3f}")
print(f"  Baseline PR (Decrease prevalence) = {baseline_pr:.3f}")
print(f"  Improvement over baseline = {(pr_auc_mean - baseline_pr) / baseline_pr * 100:.0f}%")

# ==============================================================================
# STEP 11: SAVE OUTPUTS (SAME FILENAMES - NO DOWNSTREAM IMPACT)
# ==============================================================================

print("\n" + "="*60)
print("STEP 11: SAVING OUTPUTS")
print("="*60)

# Add metadata to coefficient table
coef_df['Model_Family'] = 'weighted GLM binomial'
coef_df['Class_Weighting'] = 'balanced'
coef_df['Weight_Decrease'] = weight_pos
coef_df['Weight_NotDecrease'] = weight_neg
coef_df['Inference_Method'] = 'weighted statsmodels GLM binomial (class-balanced)'
coef_df['SPEI_mean'] = spei_mean
coef_df['SPEI_std'] = spei_std
coef_df['TET_mean'] = tet_mean
coef_df['TET_std'] = tet_values.std()

# Add cross-validation and sample metadata to first row only
coef_df.loc[0, 'N_Total'] = n_total
coef_df.loc[0, 'N_Decrease'] = n_decrease
coef_df.loc[0, 'N_NotDecrease'] = n_not_decrease
coef_df.loc[0, 'ROC_AUC'] = auc_mean
coef_df.loc[0, 'ROC_AUC_Std'] = auc_std
coef_df.loc[0, 'PR_AUC'] = pr_auc_mean
coef_df.loc[0, 'Baseline_PR'] = baseline_pr
coef_df.loc[0, 'Valid_Folds'] = valid_folds
coef_df.loc[0, 'Skipped_Folds'] = skipped_folds

# Save coefficient table
coef_output_path = os.path.join(OUTPUT_DIR, "WUE_tra_SPEI48_TET_interaction_model_summary.csv")
coef_df.to_csv(coef_output_path, index=False)
print(f"✅ Saved coefficient table: {coef_output_path}")

# Save probability table
prob_output_path = os.path.join(OUTPUT_DIR, "WUE_tra_SPEI48_TET_probability_table.csv")
prob_df.to_csv(prob_output_path, index=False)
print(f"✅ Saved probability table: {prob_output_path}")

# Save interpretation text file
interpretation_path = os.path.join(OUTPUT_DIR, "WUE_tra_SPEI48_TET_interaction_interpretation.txt")

with open(interpretation_path, 'w', encoding='utf-8') as f:
    f.write("="*80 + "\n")
    f.write("WUE_T DECLINE ~ SPEI-48 × T:ET RATIO INTERACTION ANALYSIS\n")
    f.write("Weighted GLM (class-balanced) - Matches Step 2 workflow\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"DATA SUMMARY:\n")
    f.write(f"  Observations: {len(df_model):,}\n")
    f.write(f"  Unique sites: {df_model['site_name'].nunique()}\n")
    f.write(f"  Decrease events: {n_decrease} ({pct_decrease:.1f}%)\n")
    f.write(f"  Class weights: Decrease={weight_pos:.4f}, NotDecrease={weight_neg:.4f}\n")
    f.write(f"  T:ET column used: {TET_COL}\n")
    f.write(f"  T:ET range: [{tet_min:.4f}, {tet_max:.4f}]\n")
    f.write(f"  SPEI-48 range: [{spei_min:.2f}, {spei_max:.2f}]\n\n")
    
    f.write("MODEL COEFFICIENTS (weighted GLM, class-balanced):\n")
    f.write("-"*60 + "\n")
    f.write(coef_df[['Variable', 'Coefficient', 'SE', 'CI_lower', 'CI_upper', 'p_value', 'Odds_Ratio']].to_string(index=False, float_format=lambda x: f"{x:.4f}") + "\n\n")
    
    f.write("INFERENCE METHOD:\n")
    f.write("  weighted statsmodels GLM binomial (class-balanced)\n")
    f.write("  SE: standard SE (not cluster-robust)\n\n")
    
    f.write("CROSS-VALIDATION PERFORMANCE (site-blocked):\n")
    f.write(f"  ROC-AUC = {auc_mean:.3f} ± {auc_std:.3f}\n")
    f.write(f"  PR-AUC = {pr_auc_mean:.3f} (baseline = {baseline_pr:.3f})\n")
    f.write(f"  Valid folds: {valid_folds}, Skipped: {skipped_folds}\n\n")
    
    f.write("INTERPRETATION:\n")
    f.write("-"*60 + "\n")
    if spei_p < 0.05:
        if spei_coef < 0:
            f.write("✓ SPEI-48 effect: SIGNIFICANT NEGATIVE\n")
            f.write("  → Lower SPEI (drought) INCREASES probability of WUE_T decline\n\n")
        else:
            f.write("✓ SPEI-48 effect: SIGNIFICANT POSITIVE\n")
            f.write("  → Higher SPEI (wetter) INCREASES probability of WUE_T decline\n\n")
    else:
        f.write("✗ SPEI-48 effect: NOT SIGNIFICANT\n\n")
    
    if interaction_p < 0.05:
        f.write("✓ SPEI-48 × T:ET INTERACTION: SIGNIFICANT\n")
        f.write(f"  p = {interaction_p:.4f}, coefficient = {interaction_coef:.4f}\n")
        f.write("  → The drought effect DIFFERS across the T:ET gradient\n")
        if interaction_coef < 0:
            f.write("  → NEGATIVE INTERACTION: SPEI-48 slope becomes MORE NEGATIVE at higher T:ET\n")
            f.write("  → Because lower SPEI represents drier conditions, this suggests\n")
            f.write("    STRONGER drought-associated WUE_T decline risk at higher T:ET\n\n")
        else:
            f.write("  → POSITIVE INTERACTION: SPEI-48 slope becomes LESS NEGATIVE at higher T:ET\n")
            f.write("  → Because lower SPEI represents drier conditions, this suggests\n")
            f.write("    WEAKER drought-associated WUE_T decline risk at higher T:ET\n\n")
    else:
        f.write("✗ SPEI-48 × T:ET INTERACTION: NOT SIGNIFICANT\n")
        f.write(f"  p = {interaction_p:.4f}\n")
        f.write("  → No evidence that drought effect depends on T:ET ratio\n\n")
    
    f.write("PREDICTED PROBABILITIES (class-balanced fitted probabilities):\n")
    f.write("-"*60 + "\n")
    f.write("Note: These are class-balanced fitted probabilities, not raw prevalence.\n")
    f.write("The 0.50 threshold is a class-balanced decision boundary.\n\n")
    for _, row in prob_df.iterrows():
        f.write(f"  SPEI-48 = {row['SPEI-48']:.1f}:\n")
        f.write(f"    Low T:ET (25th %ile): P = {row['Class-balanced P(Decrease) | Low T:ET']:.3f}\n")
        f.write(f"    High T:ET (75th %ile): P = {row['Class-balanced P(Decrease) | High T:ET']:.3f}\n")
        f.write(f"    Difference (Low - High): {row['Difference (Low - High)']:+.3f}\n\n")

print(f"✅ Saved interpretation text: {interpretation_path}")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

print("\n" + "="*80)
print("Q4 ANALYSIS COMPLETE - Weighted GLM (class-balanced)")
print("="*80)
print(f"\n📁 Output files saved to: {OUTPUT_DIR}")
print("   • WUE_tra_SPEI48_TET_interaction_model_summary.csv")
print("   • WUE_tra_SPEI48_TET_probability_table.csv")
print("   • WUE_tra_SPEI48_TET_interaction_interpretation.txt")
print("\n✅ Model uses weighted GLM with class-balanced weights (matches Step 2)")
print("✅ Inference: standard SE from weighted GLM (not cluster-robust)")
print("✅ CV uses sklearn LogisticRegression(class_weight='balanced')")