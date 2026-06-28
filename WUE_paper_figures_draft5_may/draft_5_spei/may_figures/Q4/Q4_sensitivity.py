
### show in sensitivity how diffrent model with diffrent split perfromce 


# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:30:09 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
PUBLICATION FIGURES & RESULTS REPORTING - WUE_TRA LOGISTIC REGRESSION
TWO PANELS WITH IMPROVED STRUCTURE:
Panel A: Cross-validated predictions (main result)
Panel B: Threshold sensitivity / Model robustness (supports main result)
Model: Decrease vs Increase+NoChange (NoChange kept as baseline)
Percentile-based classification (10th-90th percentile)
SPEI-48 (48-month drought)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
import warnings
import os

warnings.filterwarnings('ignore')

# ==============================================================================
# FIGURE PARAMETERS
# ==============================================================================

FIGURE_WIDTH = 16
FIGURE_HEIGHT = 7

LABEL_FONT_SIZE = 20
TITLE_FONT_SIZE = 22
TICK_FONT_SIZE = 18
ANNOTATION_FONT_SIZE = 14
LEGEND_FONT_SIZE = 16

LINE_COLOR = 'black'
POINT_COLOR = '#606060'
CV_POINT_COLOR = '#d62728'
LINE_WIDTH = 3.0
SCATTER_ALPHA = 0.15
SCATTER_SIZE = 15
CV_SCATTER_ALPHA = 0.25
CV_SCATTER_SIZE = 15
PANEL_SPACING = 0.35

# ==============================================================================
# THRESHOLD SENSITIVITY DATA (from your diagnostic)
# ==============================================================================

# Data from your diagnostic output (10-90, 15-85, 20-80, 25-75, etc.)
threshold_sensitivity = {
    'Percentile Range': ['2.5-97.5', '5-95', '10-90', '15-85', '20-80', '25-75', '0-100'],
    'Coverage': [95, 90, 80, 70, 60, 50, 100],
    'NoChange %': [91.5, 86.9, 77.5, 67.9, 58.6, 49.0, 1.6],
    'N_D_I': [136, 210, 362, 516, 666, 820, 1583],
    'Coefficient': [-0.077, -0.017, -0.104, -0.079, -0.081, -0.066, -0.074],
    'AUC': [0.530, 0.470, 0.567, 0.556, 0.560, 0.553, 0.553]
}

df_sensitivity = pd.DataFrame(threshold_sensitivity)

# ==============================================================================
# LOAD DATA
# ==============================================================================

def load_data():
    summary_path = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv"
    data_path = r"M:\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
    
    try:
        summary_df = pd.read_csv(summary_path)
        df = pd.read_csv(data_path)
        print("✅ Loaded data from Step 2")
        return summary_df, df
    except Exception as e:
        print(f"❌ Could not load: {e}")
        return None, None

# ==============================================================================
# PREPARE DATA
# ==============================================================================

def prepare_data(df):
    metric = 'WUE_tra'
    class_col = f"{metric}_minus_NNmed_SPEI_48_Class"
    spei_col = 'SPEI_48'
    
    mask = df[spei_col].notna()
    filtered_df = df[mask].copy()
    filtered_df['target'] = (filtered_df[class_col] == 'Decrease').astype(int)
    filtered_df = filtered_df.dropna(subset=['target'])
    
    X = filtered_df[[spei_col]].values.flatten()
    y = filtered_df['target'].values
    sites = filtered_df['site_name'].values
    
    return X, y, sites, filtered_df

# ==============================================================================
# GET SITE-BLOCKED CV PREDICTIONS (held-out predictions - MAIN RESULT)
# ==============================================================================

def get_cv_predictions(X, y, sites):
    X_2d = X.reshape(-1, 1)
    logo = LeaveOneGroupOut()
    cv_predictions = np.full(len(y), np.nan)
    valid_count = 0
    
    for train_idx, test_idx in logo.split(X_2d, y, groups=sites):
        X_train, X_test = X_2d[train_idx], X_2d[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        if len(np.unique(y_test)) < 2:
            continue
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        model = LogisticRegression(class_weight='balanced', solver='liblinear', random_state=42, max_iter=1000)
        model.fit(X_train_scaled, y_train)
        
        cv_predictions[test_idx] = model.predict_proba(X_test_scaled)[:, 1]
        valid_count += len(test_idx)
    
    return cv_predictions, valid_count

# ==============================================================================
# CREATE TWO-PANEL FIGURE (IMPROVED)
# ==============================================================================

def create_two_panel_figure():
    print("\n" + "="*80)
    print("CREATING TWO-PANEL FIGURE")
    print("Panel A: Cross-validated predictions (main result)")
    print("Panel B: Threshold sensitivity / Model robustness")
    print("="*80)
    
    summary_df, df = load_data()
    if summary_df is None or df is None:
        return None
    
    # Prepare data
    X, y, sites, filtered_df = prepare_data(df)
    
    # Extract model parameters from Step 2
    coef_step2 = summary_df['Coefficient'].values[0]
    intercept_step2 = summary_df['Intercept'].values[0]
    roc_auc = summary_df['ROC_AUC'].values[0]
    n_total = summary_df['N_Total'].values[0]
    pct_decrease = summary_df['Pct_Decrease'].values[0]
    
    # ==========================================================================
    # GET CV PREDICTIONS (MAIN RESULT)
    # ==========================================================================
    print("Computing site-blocked CV predictions (main result)...")
    cv_predictions, n_cv = get_cv_predictions(X, y, sites)
    
    # Remove NaN from CV predictions for plotting
    valid_cv_mask = ~np.isnan(cv_predictions)
    cv_X = X[valid_cv_mask]
    cv_preds = cv_predictions[valid_cv_mask]
    
    print(f"  Panel A (CV-held-out predictions) n = {n_cv}")
    print(f"  Panel B (Threshold sensitivity from diagnostic)")
    
    # Create figure with two panels
    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig.subplots_adjust(wspace=PANEL_SPACING)
    
    # Common SPEI range for S-curve
    spei_range = np.linspace(-3.5, 3.5, 200)
    log_odds = intercept_step2 + coef_step2 * spei_range
    probabilities = 1 / (1 + np.exp(-log_odds))
    
    # Unique SPEI values for points
    unique_spei = np.unique(X)
    log_odds_unique = intercept_step2 + coef_step2 * unique_spei
    probs_unique = 1 / (1 + np.exp(-log_odds_unique))
    sort_idx = np.argsort(unique_spei)
    unique_spei_sorted = unique_spei[sort_idx]
    probs_unique_sorted = probs_unique[sort_idx]
    
    # ==========================================================================
    # PANEL A: CROSS-VALIDATED PREDICTIONS (MAIN RESULT)
    # ==========================================================================
    ax1 = axes[0]
    
    # Plot S-curve (reference line)
    ax1.plot(spei_range, probabilities, color=LINE_COLOR, linewidth=LINE_WIDTH, label='Logistic fit')
    
    # Plot unique SPEI points
    ax1.scatter(unique_spei_sorted, probs_unique_sorted,
                facecolor=POINT_COLOR, edgecolor=POINT_COLOR,
                s=80, alpha=0.8, marker='o', zorder=2, label='Unique SPEI values')
    
    # Plot CV predictions (jittered points from held-out sites)
    jitter_x_cv = np.random.normal(0, 0.03, len(cv_X))
    jitter_y_cv = np.random.normal(0, 0.01, len(cv_preds))
    ax1.scatter(cv_X + jitter_x_cv, cv_preds + jitter_y_cv,
                alpha=CV_SCATTER_ALPHA, color=CV_POINT_COLOR, s=CV_SCATTER_SIZE, 
                edgecolor='none', rasterized=True, label='CV predictions (held-out sites)')
    
    # Reference lines
    ax1.axvspan(-0.5, 0.5, alpha=0.1, color='gray', label='Near-normal (NN) region')
    ax1.axhline(0.5, color='black', linestyle='--', alpha=0.4, linewidth=1.5)
    
    # Labels
    ax1.set_xlabel('SPEI-48', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax1.set_ylabel('P (WUE$_T$ decrease)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax1.set_title('(a) Cross-validated predictions (main result)', fontsize=TITLE_FONT_SIZE, fontweight='bold', loc='center', pad=10)
    
    ax1.set_ylim(0.25, 0.65)
    ax1.set_xlim(-3.2, 3.2)
    ax1.grid(True, alpha=0.2, linestyle='--')
    ax1.tick_params(axis='both', labelsize=TICK_FONT_SIZE, width=1.5, length=8)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.5)
    
    # Legend
    ax1.legend(loc='lower right', fontsize=LEGEND_FONT_SIZE, framealpha=0.9)
    
    # Annotation for Panel A (CV results)
    stats_text_a = f"n = {n_cv:,} (CV held-out)\n"
    stats_text_a += f"Decrease (ground truth): {pct_decrease:.1f}%\n"
    stats_text_a += f"CV ROC-AUC = {roc_auc:.3f}\n"
    stats_text_a += f"β = {coef_step2:+.4f}"
    
    ax1.text(0.02, 0.98, stats_text_a, transform=ax1.transAxes, fontsize=ANNOTATION_FONT_SIZE,
             verticalalignment='top', 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='gray', linewidth=1))
    
    # ==========================================================================
    # PANEL B: THRESHOLD SENSITIVITY (MODEL ROBUSTNESS)
    # ==========================================================================
    ax2 = axes[1]
    
    # Find the index of 10-90 (our selected threshold)
    selected_idx = df_sensitivity['Percentile Range'] == '10-90'
    
    # Plot 1: AUC by percentile range (bar chart)
    x_pos = np.arange(len(df_sensitivity))
    colors = ['#d62728' if r == '10-90' else '#1f77b4' for r in df_sensitivity['Percentile Range']]
    bars = ax2.bar(x_pos, df_sensitivity['AUC'], color=colors, alpha=0.7, edgecolor='black')
    ax2.axhline(0.5, color='red', linestyle='--', linewidth=2, label='Random (AUC=0.5)')
    ax2.set_xlabel('Percentile range', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax2.set_ylabel('ROC-AUC', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax2.set_title('(b) Threshold sensitivity (model robustness)', fontsize=TITLE_FONT_SIZE, fontweight='bold', loc='center', pad=10)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(df_sensitivity['Percentile Range'], rotation=45, ha='right', fontsize=TICK_FONT_SIZE-4)
    ax2.set_ylim(0.4, 0.65)
    ax2.tick_params(axis='y', labelsize=TICK_FONT_SIZE)
    ax2.grid(True, alpha=0.2, axis='y')
    
    # Add value labels on bars
    for bar, auc in zip(bars, df_sensitivity['AUC']):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005, 
                f'{auc:.1f}', ha='center', va='bottom', fontsize=ANNOTATION_FONT_SIZE)
    
    # Highlight the selected threshold (10-90)
    ax2.text(selected_idx.idxmax(), df_sensitivity.loc[selected_idx, 'AUC'].values[0] + 0.02, 
            'SELECTED', ha='center', va='bottom', fontsize=ANNOTATION_FONT_SIZE, 
            color='red', fontweight='bold')
    
    # Legend for panel B
    ax2.legend(loc='lower right', fontsize=LEGEND_FONT_SIZE)
    
    # Annotation for Panel B
    stats_text_b = f"Coefficient remains NEGATIVE across ALL thresholds\n"
    stats_text_b += f"Selected range: 10-90 (AUC={df_sensitivity.loc[selected_idx, 'AUC'].values[0]:.3f})"
    
    ax2.text(0.02, 0.98, stats_text_b, transform=ax2.transAxes, fontsize=ANNOTATION_FONT_SIZE,
             verticalalignment='top', 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='gray', linewidth=1))
    
    # ==========================================================================
    # SAVE FIGURE (SAME FILE NAMES - NO DOWNSTREAM IMPACT)
    # ==========================================================================
    output_dir = r"M:\Research\WUE_CUE\data_products\results\linear_models"
    
    png_path = os.path.join(output_dir, "Publication_Figure_WUE_tra_Training_Testing.png")
    pdf_path = os.path.join(output_dir, "Publication_Figure_WUE_tra_Training_Testing.pdf")
    
    plt.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"\n✅ Saved figure: {png_path}")
    print(f"✅ Saved figure: {pdf_path}")
    
    plt.show()
    
    return fig

# ==============================================================================
# RESULTS REPORTING (CONSOLE ONLY)
# ==============================================================================

def report_results():
    print("\n" + "="*80)
    print("RESULTS REPORTING - WUE_TRA LOGISTIC REGRESSION")
    print("="*80)
    
    summary_path = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv"
    prob_path = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Probability_Table_10_90.csv"
    
    try:
        summary_df = pd.read_csv(summary_path)
        prob_df = pd.read_csv(prob_path)
        print("✅ Loaded model results from Step 2")
    except Exception as e:
        print(f"❌ Could not load results: {e}")
        return None, None
    
    coef = summary_df['Coefficient'].values[0]
    roc_auc = summary_df['ROC_AUC'].values[0]
    pr_auc = summary_df['PR_AUC'].values[0]
    baseline_pr = summary_df['Baseline_PR'].values[0]
    n_total = summary_df['N_Total'].values[0]
    n_decrease = summary_df['N_Decrease'].values[0]
    pct_decrease = summary_df['Pct_Decrease'].values[0]
    
    se = 1 / np.sqrt(n_total)
    z_score = coef / se
    p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
    coef_ci_lower = coef - 1.96 * se
    coef_ci_upper = coef + 1.96 * se
    
    significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
    
    print("\n" + "="*60)
    print("MODEL STATISTICS")
    print("="*60)
    print(f"\nSample size: n = {n_total:,} observations")
    print(f"Class distribution: Decrease = {n_decrease} ({pct_decrease:.1f}%), Not Decrease = {n_total - n_decrease} ({100-pct_decrease:.1f}%)")
    
    print(f"\nCoefficient (SPEI-48): β₁ = {coef:.4f}")
    print(f"95% CI: [{coef_ci_lower:.4f}, {coef_ci_upper:.4f}]")
    print(f"p-value: {p_value:.2e} {significance}")
    print(f"Odds ratio (per 1-unit decrease in SPEI): {np.exp(-coef):.3f}")
    
    print(f"\nModel performance (site-blocked CV):")
    print(f"  ROC-AUC = {roc_auc:.3f} ± {summary_df['ROC_AUC_Std'].values[0]:.3f}")
    print(f"  PR-AUC = {pr_auc:.3f} (baseline = {baseline_pr:.3f})")
    print(f"  Improvement over random = {(pr_auc - baseline_pr) / baseline_pr * 100:.0f}%")
    
    print("\n" + "="*60)
    print("PROBABILITY TABLE")
    print("="*60)
    print("\nSPEI-48\tP(WUE$_T$ decrease)")
    print("-" * 35)
    for _, row in prob_df.iterrows():
        print(f"{row['SPEI']:.1f}\t\t{row['P_Decrease_SPEI48']:.3f}")
    
    return summary_df, prob_df

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

def publication_workflow():
    print("\n" + "="*80)
    print("PUBLICATION WORKFLOW - WUE_TRA LOGISTIC REGRESSION")
    print("Percentile-based classification (10th-90th percentile)")
    print("Model: Decrease vs Increase+NoChange (NoChange kept as baseline)")
    print("SPEI-48 (48-month drought)")
    print("="*80)
    
    print("\n📝 GENERATING RESULTS FOR MANUSCRIPT (CONSOLE OUTPUT)...")
    print("-" * 60)
    summary_df, prob_df = report_results()
    
    print("\n\n🎨 CREATING TWO-PANEL FIGURE...")
    print("   Panel A: Cross-validated predictions (main result)")
    print("   Panel B: Threshold sensitivity (model robustness)")
    print("-" * 60)
    fig = create_two_panel_figure()
    
    print("\n" + "="*80)
    print("WORKFLOW COMPLETE")
    print("="*80)
    print("\n📁 Output files (SAME NAMES - No downstream impact):")
    print("   • Publication_Figure_WUE_tra_Training_Testing.png")
    print("   • Publication_Figure_WUE_tra_Training_Testing.pdf")
    print("   • WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv (from Step 2)")
    print("   • WUE_tra_Probability_Table_10_90.csv (from Step 2)")
    
    return summary_df, prob_df, fig

# ==============================================================================
# EXECUTE
# ==============================================================================

if __name__ == "__main__":
    publication_workflow()

# To run manually:
# summary_df, prob_df, fig = publication_workflow()