# -*- coding: utf-8 -*-
"""
CHUNK 1: PERCENTILE THRESHOLD SENSITIVITY ANALYSIS - DATA PROCESSING (FIXED)
=============================================================================
Now saves Intercept to summary table for correct probability curves.
=============================================================================
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import LeaveOneGroupOut
import statsmodels.api as sm
import warnings
import os

warnings.filterwarnings('ignore')

# =============================================================================
# PARAMETERS
# =============================================================================

THRESHOLDS = [
    (5, 95, "5_95"),
    (10, 90, "10_90"),
    (15, 85, "15_85"),
    (20, 80, "20_80"),
    (25, 75, "25_75")
]

METRIC = 'WUE_tra'
SPEI_COL = 'SPEI_48'
SPEI_TIMESCALE = 48
ANOMALY_COL = f"{METRIC}_minus_NNmed_SPEI_{SPEI_TIMESCALE}"

INPUT_DATA_PATH = r"M:\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
OUTPUT_BASE_DIR = r"M:\Research\WUE_CUE\data_products\results\sensitivity_analysis_percentiles"

print("="*80)
print("CHUNK 1: PERCENTILE THRESHOLD SENSITIVITY ANALYSIS")
print("="*80)

os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
for lower, upper, label in THRESHOLDS:
    subdir = os.path.join(OUTPUT_BASE_DIR, f"{lower}_{upper}")
    os.makedirs(subdir, exist_ok=True)

# =============================================================================
# LOAD DATA
# =============================================================================

print("\nSTEP 1: LOADING DATA")
df = pd.read_csv(INPUT_DATA_PATH)
print(f"  Loaded: {df.shape[0]} rows, {df.shape[1]} columns")

mask = df[SPEI_COL].notna() & df[ANOMALY_COL].notna()
base_df = df[mask].copy()
print(f"  Valid rows: {len(base_df)}")

# =============================================================================
# FUNCTIONS
# =============================================================================

def classify_anomalies(df, anomaly_col, lower_pct, upper_pct):
    df_copy = df.copy()
    class_col = f"{anomaly_col}_Class_{lower_pct}_{upper_pct}"
    df_copy[class_col] = np.nan
    
    for site in df_copy['site_name'].unique():
        site_mask = df_copy['site_name'] == site
        site_anomalies = df_copy.loc[site_mask, anomaly_col].dropna().values
        
        if len(site_anomalies) < 3:
            df_copy.loc[site_mask, class_col] = "NoChange"
            continue
        
        lower_thresh = np.percentile(site_anomalies, lower_pct)
        upper_thresh = np.percentile(site_anomalies, upper_pct)
        
        for idx in df_copy.loc[site_mask].index:
            anomaly_val = df_copy.loc[idx, anomaly_col]
            if pd.isna(anomaly_val):
                df_copy.loc[idx, class_col] = np.nan
            elif anomaly_val < lower_thresh:
                df_copy.loc[idx, class_col] = "Decrease"
            elif anomaly_val > upper_thresh:
                df_copy.loc[idx, class_col] = "Increase"
            else:
                df_copy.loc[idx, class_col] = "NoChange"
    
    return df_copy, class_col

def fit_weighted_logistic(df, class_col, spei_col):
    model_df = df[df[spei_col].notna() & df[class_col].notna()].copy()
    model_df['target'] = (model_df[class_col] == 'Decrease').astype(int)
    model_df = model_df.dropna(subset=['target'])
    
    X_raw = model_df[[spei_col]].values.reshape(-1, 1)
    y = model_df['target'].values
    sites = model_df['site_name'].values
    
    n_total = len(y)
    n_decrease = y.sum()
    n_increase = (model_df[class_col] == 'Increase').sum()
    n_nochange = (model_df[class_col] == 'NoChange').sum()
    pct_decrease = n_decrease / n_total * 100
    
    n_pos = n_decrease
    n_neg = n_total - n_decrease
    weight_pos = n_total / (2 * n_pos) if n_pos > 0 else 1
    weight_neg = n_total / (2 * n_neg) if n_neg > 0 else 1
    sample_weights = np.where(y == 1, weight_pos, weight_neg)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_raw)
    spei_mean = float(X_raw.mean())
    spei_std = float(X_raw.std())
    
    # Weighted GLM
    X_sm = sm.add_constant(X_scaled)
    glm_model = sm.GLM(y, X_sm, family=sm.families.Binomial(), freq_weights=sample_weights)
    glm_result = glm_model.fit()
    
    coef_std = glm_result.params[1]
    intercept = glm_result.params[0]  # THIS IS THE CORRECT FITTED INTERCEPT
    se = glm_result.bse[1]
    p_value = glm_result.pvalues[1]
    ci_lower, ci_upper = glm_result.conf_int()[1]
    coef_raw = coef_std / scaler.scale_[0]
    odds_ratio_1sd = np.exp(-coef_std)
    odds_ratio_1unit = np.exp(-coef_raw)
    
    # Cross-validation
    logo = LeaveOneGroupOut()
    auc_scores = []
    pr_auc_scores = []
    valid_folds = 0
    
    for train_idx, test_idx in logo.split(X_raw, y, groups=sites):
        X_train_raw, X_test_raw = X_raw[train_idx], X_raw[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        if len(np.unique(y_test)) < 2:
            continue
        
        scaler_fold = StandardScaler()
        X_train_scaled = scaler_fold.fit_transform(X_train_raw)
        X_test_scaled = scaler_fold.transform(X_test_raw)
        
        model_fold = LogisticRegression(class_weight='balanced', solver='liblinear', 
                                       random_state=42, max_iter=1000)
        model_fold.fit(X_train_scaled, y_train)
        y_pred = model_fold.predict_proba(X_test_scaled)[:, 1]
        
        if len(np.unique(y_pred)) < 2:
            continue
        
        auc_scores.append(roc_auc_score(y_test, y_pred))
        pr_auc_scores.append(average_precision_score(y_test, y_pred))
        valid_folds += 1
    
    auc_mean = np.mean(auc_scores) if auc_scores else np.nan
    auc_std = np.std(auc_scores) if auc_scores else np.nan
    pr_auc_mean = np.mean(pr_auc_scores) if pr_auc_scores else np.nan
    baseline_pr = y.mean()
    
    return {
        'n_total': n_total, 'n_decrease': n_decrease, 'n_increase': n_increase,
        'n_nochange': n_nochange, 'pct_decrease': pct_decrease,
        'coef_standardized': coef_std, 'coef_raw': coef_raw, 'intercept': intercept,  # SAVED!
        'se': se, 'p_value': p_value, 'ci_lower': ci_lower, 'ci_upper': ci_upper,
        'odds_ratio_1sd': odds_ratio_1sd, 'odds_ratio_1unit': odds_ratio_1unit,
        'roc_auc': auc_mean, 'roc_auc_std': auc_std, 'pr_auc': pr_auc_mean,
        'baseline_pr': baseline_pr, 'valid_folds': valid_folds,
        'spei_mean': spei_mean, 'spei_std': spei_std,
        'weight_decrease': weight_pos, 'weight_not_decrease': weight_neg,
        'scaler': scaler, 'glm_result': glm_result, 'X_raw': X_raw, 'y': y
    }

# =============================================================================
# RUN ANALYSIS
# =============================================================================

print("\nSTEP 2: RUNNING ANALYSIS")
all_results = {}
threshold_info = {}

for lower_pct, upper_pct, label in THRESHOLDS:
    print(f"\n  Processing: {lower_pct}/{upper_pct}...")
    
    df_classified, class_col = classify_anomalies(base_df, ANOMALY_COL, lower_pct, upper_pct)
    mask_valid = df_classified[class_col].notna()
    df_model = df_classified[mask_valid].copy()
    results = fit_weighted_logistic(df_model, class_col, SPEI_COL)
    
    all_results[label] = results
    threshold_info[label] = (lower_pct, upper_pct)
    
    print(f"    n={results['n_total']}, Decrease={results['pct_decrease']:.1f}%, "
          f"β={results['coef_standardized']:+.3f}, intercept={results['intercept']:.3f}, "
          f"p={results['p_value']:.2e}, AUC={results['roc_auc']:.3f}")

# =============================================================================
# SAVE RESULTS (NOW INCLUDES INTERCEPT)
# =============================================================================

print("\nSTEP 3: SAVING RESULTS")

# Create master summary table with Intercept
summary_data = []
for label, results in all_results.items():
    lower_pct, upper_pct = threshold_info[label]
    summary_data.append({
        'Threshold': f"{lower_pct}/{upper_pct}",
        'Lower_Percentile': lower_pct,
        'Upper_Percentile': upper_pct,
        'N_Total': results['n_total'],
        'N_Decrease': results['n_decrease'],
        'N_Increase': results['n_increase'],
        'N_NoChange': results['n_nochange'],
        'Pct_Decrease': round(results['pct_decrease'], 1),
        'Coefficient_Standardized': results['coef_standardized'],
        'Intercept': results['intercept'],  # ADDED!
        'Coefficient_SE': results['se'],
        'CI_Lower': results['ci_lower'],
        'CI_Upper': results['ci_upper'],
        'P_Value': results['p_value'],
        'ROC_AUC': results['roc_auc'],
        'ROC_AUC_Std': results['roc_auc_std'],
        'PR_AUC': results['pr_auc'],
        'Baseline_PR': results['baseline_pr'],
        'Odds_Ratio_1SD': results['odds_ratio_1sd'],
        'Odds_Ratio_1unit': results['odds_ratio_1unit'],
        'SPEI_Mean': results['spei_mean'],
        'SPEI_Std': results['spei_std']
    })

master_df = pd.DataFrame(summary_data)
master_path = os.path.join(OUTPUT_BASE_DIR, "Sensitivity_Summary_AllThresholds.csv")
master_df.to_csv(master_path, index=False)
print(f"  Saved: {master_path}")
print(f"  Columns now include: Intercept")

print("\n✅ CHUNK 1 COMPLETE - Results saved with Intercept")
print(f"   Results stored in: {OUTPUT_BASE_DIR}")