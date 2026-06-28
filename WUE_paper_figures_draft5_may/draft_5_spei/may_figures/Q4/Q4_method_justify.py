# -*- coding: utf-8 -*-
"""
Created on Sun May  3 19:26:30 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
DIAGNOSTICS FOR PERCENTILE-BASED CLASSIFICATION
To determine if 10th-90th percentile is appropriate or needs adjustment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Output directory for diagnostic figures
DIAG_DIR = r"M:\Research\WUE_CUE\data_products\results\diagnostics"
os.makedirs(DIAG_DIR, exist_ok=True)

# Load the output from Step 1
data_path = r"M:\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
df = pd.read_csv(data_path)
print(f"✅ Loaded: {df.shape}")

# Focus on WUE_tra, SPEI-48
metric = 'WUE_tra'
class_col = f"{metric}_minus_NNmed_SPEI_48_Class"
anomaly_col = f"{metric}_minus_NNmed_SPEI_48"

print("\n" + "="*80)
print("DIAGNOSTIC 1: PER-SITE CLASS BALANCE")
print("="*80)

# =============================================================================
# DIAGNOSTIC 1: Per-site class balance (MOST IMPORTANT)
# =============================================================================

site_summary = []

for site in df['site_name'].unique():
    site_df = df[df['site_name'] == site]
    
    n_total = len(site_df)
    n_decrease = (site_df[class_col] == 'Decrease').sum()
    n_increase = (site_df[class_col] == 'Increase').sum()
    n_nochange = (site_df[class_col] == 'NoChange').sum()
    
    pct_decrease = n_decrease / n_total * 100 if n_total > 0 else 0
    pct_increase = n_increase / n_total * 100 if n_total > 0 else 0
    pct_nochange = n_nochange / n_total * 100 if n_total > 0 else 0
    
    site_summary.append({
        'site': site,
        'n_total': n_total,
        'n_decrease': n_decrease,
        'n_increase': n_increase,
        'n_nochange': n_nochange,
        'pct_decrease': pct_decrease,
        'pct_increase': pct_increase,
        'pct_nochange': pct_nochange
    })

site_df_summary = pd.DataFrame(site_summary)

print("\n📊 PER-SITE CLASS DISTRIBUTION SUMMARY:")
print(f"  Total sites: {len(site_df_summary)}")
print(f"  Sites with any Decrease: {(site_df_summary['n_decrease'] > 0).sum()}")
print(f"  Sites with any Increase: {(site_df_summary['n_increase'] > 0).sum()}")
print(f"  Sites with any NoChange: {(site_df_summary['n_nochange'] > 0).sum()}")

print("\n📊 Distribution of % Decrease across sites:")
print(f"  Mean % Decrease: {site_df_summary['pct_decrease'].mean():.1f}%")
print(f"  Median % Decrease: {site_df_summary['pct_decrease'].median():.1f}%")
print(f"  Min % Decrease: {site_df_summary['pct_decrease'].min():.1f}%")
print(f"  Max % Decrease: {site_df_summary['pct_decrease'].max():.1f}%")

# Flag potential issues
sites_no_extremes = site_df_summary[(site_df_summary['n_decrease'] == 0) & (site_df_summary['n_increase'] == 0)]
if len(sites_no_extremes) > 0:
    print(f"\n⚠️ WARNING: {len(sites_no_extremes)} sites have 0% Decrease AND 0% Increase (all NoChange)")
    print(f"   These sites: {sites_no_extremes['site'].tolist()}")
else:
    print("\n✅ GOOD: All sites have at least some Decrease or Increase observations")

print("\n📊 Sample sites (first 10):")
print(site_df_summary[['site', 'n_total', 'pct_decrease', 'pct_increase', 'pct_nochange']].head(10).to_string(index=False))

# =============================================================================
# DIAGNOSTIC 2: Anomaly distribution vs thresholds (for representative sites)
# =============================================================================
print("\n" + "="*80)
print("DIAGNOSTIC 2: ANOMALY DISTRIBUTION VS THRESHOLDS")
print("="*80)

# Get anomaly data
anomaly_data = df[[anomaly_col, 'site_name']].dropna()

# Select representative sites (one high, one medium, one low variability)
site_variability = anomaly_data.groupby('site_name')[anomaly_col].std().sort_values()
sites_to_plot = [
    site_variability.index[-1],  # highest variability
    site_variability.index[len(site_variability)//2],  # median variability
    site_variability.index[0]    # lowest variability
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, site in enumerate(sites_to_plot):
    ax = axes[idx]
    site_anomalies = df[df['site_name'] == site][anomaly_col].dropna()
    
    # Calculate percentiles
    p10 = np.percentile(site_anomalies, 10)
    p90 = np.percentile(site_anomalies, 90)
    
    # Plot histogram
    ax.hist(site_anomalies, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axvline(p10, color='red', linestyle='--', linewidth=2, label=f'10th percentile ({p10:.2f})')
    ax.axvline(p90, color='green', linestyle='--', linewidth=2, label=f'90th percentile ({p90:.2f})')
    ax.axvline(0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    ax.set_xlabel('Anomaly (WUE_T - NN median)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Site: {site}\n10th-90th percentile range = [{p10:.2f}, {p90:.2f}]')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.suptitle('Diagnostic 2: Anomaly Distribution with 10th-90th Percentile Thresholds', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(DIAG_DIR, 'Diagnostic2_Anomaly_Distribution.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"✅ Saved: Diagnostic2_Anomaly_Distribution.png")

# =============================================================================
# DIAGNOSTIC 3: Sensitivity to threshold choice
# =============================================================================
print("\n" + "="*80)
print("DIAGNOSTIC 3: SENSITIVITY TO THRESHOLD CHOICE")
print("="*80)

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneGroupOut

# Define thresholds to test
thresholds = [
    (10, 90, '10-90 (current)'),
    (15, 85, '15-85 (stricter)'),
    (20, 80, '20-80 (more balanced)'),
    (25, 75, '25-75 (even more balanced)')
]

results = []

for lower, upper, label in thresholds:
    # Reclassify based on new thresholds
    classifications = []
    for site in df['site_name'].unique():
        site_mask = df['site_name'] == site
        site_anomalies = df.loc[site_mask, anomaly_col].dropna()
        
        if len(site_anomalies) < 3:
            for idx in df.loc[site_mask].index:
                classifications.append('NoChange')
            continue
        
        ci_lower = np.percentile(site_anomalies, lower)
        ci_upper = np.percentile(site_anomalies, upper)
        
        for idx in df.loc[site_mask].index:
            anomaly = df.loc[idx, anomaly_col]
            if pd.isna(anomaly):
                classifications.append(np.nan)
            elif anomaly < ci_lower:
                classifications.append('Decrease')
            elif anomaly > ci_upper:
                classifications.append('Increase')
            else:
                classifications.append('NoChange')
    
    df_temp = df.copy()
    df_temp['temp_class'] = classifications
    
    # Model: Decrease vs Increase ONLY (drop NoChange)
    mask = (df_temp['temp_class'].isin(['Decrease', 'Increase'])) & (df_temp['SPEI_48'].notna())
    df_model = df_temp[mask].copy()
    
    if len(df_model) < 50:
        print(f"\n{label}: Too few samples ({len(df_model)}), skipping...")
        continue
    
    df_model['target'] = (df_model['temp_class'] == 'Decrease').astype(int)
    
    X = df_model[['SPEI_48']].values
    y = df_model['target'].values
    sites = df_model['site_name'].values
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Fit full model
    model = LogisticRegression(class_weight='balanced', solver='liblinear', random_state=42)
    model.fit(X_scaled, y)
    coef = model.coef_[0][0]
    
    # Cross-validation AUC
    logo = LeaveOneGroupOut()
    auc_scores = []
    for train_idx, test_idx in logo.split(X, y, groups=sites):
        if len(np.unique(y[test_idx])) < 2:
            continue
        scaler_fold = StandardScaler()
        X_train = scaler_fold.fit_transform(X[train_idx])
        X_test = scaler_fold.transform(X[test_idx])
        m_fold = LogisticRegression(class_weight='balanced', solver='liblinear', random_state=42)
        m_fold.fit(X_train, y[train_idx])
        y_pred = m_fold.predict_proba(X_test)[:, 1]
        if len(np.unique(y_pred)) > 1:
            auc_scores.append(roc_auc_score(y[test_idx], y_pred))
    
    auc_mean = np.mean(auc_scores) if auc_scores else 0
    
    # Predict at SPEI = -3
    prob_at_minus3 = model.predict_proba(scaler.transform([[-3]]))[0, 1]
    
    results.append({
        'Threshold': label,
        'N_Decrease_Increase': len(df_model),
        'N_NoChange': len(df_temp) - len(df_model),
        'Pct_NoChange': (len(df_temp) - len(df_model)) / len(df_temp) * 100,
        'Coefficient': coef,
        'AUC': auc_mean,
        'P(Decrease)_SPEI_-3': prob_at_minus3
    })

# Display results
print("\n📊 SENSITIVITY ANALYSIS RESULTS:")
print("-" * 80)
print(f"{'Threshold':<25} {'N (D+I)':<12} {'NoChange %':<12} {'Coef':<10} {'AUC':<8} {'P(Dec) at -3':<12}")
print("-" * 80)
for r in results:
    print(f"{r['Threshold']:<25} {r['N_Decrease_Increase']:<12} {r['Pct_NoChange']:.1f}%{'':<8} {r['Coefficient']:+.4f}   {r['AUC']:.3f}     {r['P(Decrease)_SPEI_-3']:.3f}")

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(DIAG_DIR, 'Diagnostic3_Threshold_Sensitivity.csv'), index=False)
print(f"\n✅ Saved: Diagnostic3_Threshold_Sensitivity.csv")

# =============================================================================
# FINAL RECOMMENDATION
# =============================================================================
print("\n" + "="*80)
print("FINAL RECOMMENDATION")
print("="*80)

# Check if coefficient stays negative across thresholds
all_negative = all(r['Coefficient'] < 0 for r in results)
if all_negative:
    print("✅ Coefficient remains NEGATIVE across all thresholds → robust pattern")
else:
    print("⚠️ Coefficient sign changes → need careful interpretation")

# Check if AUC is stable
auc_values = [r['AUC'] for r in results]
auc_range = max(auc_values) - min(auc_values)
print(f"   AUC range: {auc_range:.3f} (from {min(auc_values):.3f} to {max(auc_values):.3f})")

if auc_range < 0.05:
    print("✅ AUC is stable across thresholds → robust")
else:
    print("⚠️ AUC varies considerably → choose threshold carefully")

print("\n📌 RECOMMENDATION:")
print("   If coefficient stays negative and AUC is stable → keep current 10-90 threshold")
print("   If 20-80 gives higher AUC and more balanced classes → consider switching")
print("   If results are robust, report all thresholds in supplementary")

print("\n" + "="*80)
print("DIAGNOSTICS COMPLETE!")
print("="*80)