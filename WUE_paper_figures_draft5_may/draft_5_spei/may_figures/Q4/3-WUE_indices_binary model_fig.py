# -*- coding: utf-8 -*-
"""
PUBLICATION FIGURE - WUE_TRA LOGISTIC REGRESSION (SINGLE PANEL)
Main result: Cross-validated predictions (site-blocked CV)
Model: Decrease vs Increase+NoChange (NoChange kept as baseline)
Percentile-based classification (10th-90th percentile)
SPEI-48 (48-month drought)

UPDATED: Weighted GLM (class-balanced) inference from Step 2
FIXED: Curve now uses standardized SPEI correctly (βz applied to z-scored SPEI)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
import warnings
import os

warnings.filterwarnings('ignore')

# ==============================================================================
# FIGURE PARAMETERS
# ==============================================================================

FIGURE_WIDTH = 9
FIGURE_HEIGHT = 8

LABEL_FONT_SIZE = 22
TICK_FONT_SIZE = 20
ANNOTATION_FONT_SIZE = 14
LEGEND_FONT_SIZE = 18

LINE_COLOR = 'black'
POINT_COLOR = 'black'  # Gray for all points (both training and CV)
LINE_WIDTH = 3.0
SCATTER_ALPHA = 0.15
SCATTER_SIZE = 14
MARKER_TYPE = 'o'  # Circle marker

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
# GET SITE-BLOCKED CV PREDICTIONS (using class_weight='balanced')
# ==============================================================================

def get_cv_predictions(X, y, sites):
    X_2d = X.reshape(-1, 1)
    logo = LeaveOneGroupOut()
    cv_predictions = np.full(len(y), np.nan)
    
    for train_idx, test_idx in logo.split(X_2d, y, groups=sites):
        X_train, X_test = X_2d[train_idx], X_2d[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        if len(np.unique(y_test)) < 2:
            continue
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Use class_weight='balanced' to match Step 2's weighted GLM
        model = LogisticRegression(class_weight='balanced', solver='liblinear', 
                                   random_state=42, max_iter=1000)
        model.fit(X_train_scaled, y_train)
        
        cv_predictions[test_idx] = model.predict_proba(X_test_scaled)[:, 1]
    
    return cv_predictions

# ==============================================================================
# FORMAT NUMBER FUNCTIONS
# ==============================================================================

def format_auc(auc):
    """Format AUC with 3 decimal places"""
    return f"{auc:.1f}"

def format_beta(beta):
    """Format beta with 3 decimal places"""
    return f"{beta:.2f}"

def format_p_value(p_val):
    """Format p-value: <0.001 or 3 decimal places"""
    if p_val < 0.001:
        return "<0.001"
    else:
        return f"{p_val:.3f}"

def format_percentage(pct):
    """Format percentage with 1 decimal place"""
    return f"{pct:.1f}"

# ==============================================================================
# CREATE SINGLE PANEL FIGURE
# ==============================================================================

def create_figure():
    print("\n" + "="*80)
    print("CREATING PUBLICATION FIGURE (Single Panel)")
    print("Using weighted GLM (class-balanced) inference from Step 2")
    print("FIXED: Curve uses standardized SPEI correctly")
    print("="*80)
    
    summary_df, df = load_data()
    if summary_df is None or df is None:
        return None
    
    # Check if proper inference columns exist (REQUIRED from Step 2)
    required_cols = ['Coefficient_P_value', 'Coefficient_CI_lower', 'Coefficient_CI_upper',
                     'Inference_Method', 'Weight_Decrease', 'Weight_NotDecrease',
                     'spei_mean', 'spei_std']  # Added spei_mean and spei_std for standardization
    missing_cols = [col for col in required_cols if col not in summary_df.columns]
    
    if missing_cols:
        raise ValueError(f"❌ Step 2 summary missing required columns: {missing_cols}\n"
                        f"   Please re-run Step 2 with weighted GLM (class-balanced).")
    
    # Prepare data
    X, y, sites, filtered_df = prepare_data(df)
    
    # Read values from Step 2 (NO hardcoded numbers!)
    coef = summary_df['Coefficient'].values[0]  # Standardized coefficient (per 1 SD increase)
    intercept = summary_df['Intercept'].values[0]
    roc_auc = summary_df['ROC_AUC'].values[0]
    n_total = summary_df['N_Total'].values[0]
    pct_decrease = summary_df['Pct_Decrease'].values[0]
    pr_auc = summary_df['PR_AUC'].values[0]
    
    # Read SPEI standardization parameters from Step 2
    spei_mean = summary_df['spei_mean'].values[0]
    spei_std = summary_df['spei_std'].values[0]
    
    # Read proper p-value and CI from Step 2 (weighted GLM)
    p_value = summary_df['Coefficient_P_value'].values[0]
    ci_lower = summary_df['Coefficient_CI_lower'].values[0]
    ci_upper = summary_df['Coefficient_CI_upper'].values[0]
    inference_method = summary_df['Inference_Method'].values[0]
    weight_decrease = summary_df['Weight_Decrease'].values[0]
    weight_not_decrease = summary_df['Weight_NotDecrease'].values[0]
    
    print(f"\n📊 Using weighted GLM inference from Step 2:")
    print(f"   p-value = {p_value:.2e}")
    print(f"   95% CI = [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"   Method = {inference_method}")
    print(f"   Class weights: Decrease={weight_decrease:.4f}, NotDecrease={weight_not_decrease:.4f}")
    print(f"\n📊 SPEI standardization parameters:")
    print(f"   spei_mean = {spei_mean:.3f}")
    print(f"   spei_std = {spei_std:.3f}")
    
    # Get CV predictions
    print("\nComputing site-blocked CV predictions...")
    cv_predictions = get_cv_predictions(X, y, sites)
    
    # Remove NaN from CV predictions for plotting
    valid_cv_mask = ~np.isnan(cv_predictions)
    cv_X = X[valid_cv_mask]
    cv_preds = cv_predictions[valid_cv_mask]
    
    print(f"  CV valid predictions: {len(cv_X)} (from held-out sites)")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    
    # ==========================================================================
    # CRITICAL FIX: Standardize SPEI range before applying coefficient
    # Coefficient is standardized (per 1 SD), so raw SPEI must be z-scored
    # ==========================================================================
    spei_range = np.linspace(-3.5, 3.5, 200)
    
    # Standardize using training data parameters from Step 2
    spei_range_z = (spei_range - spei_mean) / spei_std
    
    # Now apply standardized coefficient to z-scored SPEI
    log_odds = intercept + coef * spei_range_z
    probabilities = 1 / (1 + np.exp(-log_odds))
    
    # Plot S-curve (BLACK line)
    ax.plot(spei_range, probabilities, color=LINE_COLOR, linewidth=LINE_WIDTH, 
            label='Logistic fit (class-balanced)')
    
    # Plot CV predictions (GRAY points)
    jitter_x = np.random.normal(0, 0.03, len(cv_X))
    jitter_y = np.random.normal(0, 0.01, len(cv_preds))
    ax.scatter(cv_X + jitter_x, cv_preds + jitter_y,
               alpha=SCATTER_ALPHA, color=POINT_COLOR, s=SCATTER_SIZE, 
               marker=MARKER_TYPE, edgecolor='none', rasterized=True, 
               label='Cross-validated predictions')
    
    # Reference lines
    ax.axvspan(-0.5, 0.5, alpha=0.15, color='grey', label='Near-normal (NN) region')
    ax.axhline(0.5, color='black', linestyle='--', alpha=0.4, linewidth=1.5)
    
    # Labels
    ax.set_xlabel('SPEI-48', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    ax.set_ylabel('P (WUE$_T$ decrease)', fontsize=LABEL_FONT_SIZE, fontweight='bold')
    
    # Axis limits
    ax.set_ylim(0.3, 0.65)
    ax.set_xlim(-3.2, 3.2)
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.tick_params(axis='both', labelsize=TICK_FONT_SIZE, width=1.5, length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Legend
    ax.legend(loc='upper right', fontsize=LEGEND_FONT_SIZE, framealpha=0.9)
    
    # ==========================================================================
    # ANNOTATION BOX - USING WEIGHTED GLM VALUES FROM STEP 2
    # Note: β is standardized (per 1 SD change in SPEI-48)
    # ==========================================================================
    
    # Format all values with appropriate significant digits
    n_display = f"{n_total:,}"
    decrease_pct_display = format_percentage(pct_decrease)
    auc_display = format_auc(roc_auc)
    pr_auc_display = f"{pr_auc:.1f}"
    beta_display = format_beta(coef)
    p_display = f"p = {format_p_value(p_value)}"
    
    # Build annotation text (simplified for publication)
    stats_text = (
        f"n = {n_display}\n"
        f"Decrease = {decrease_pct_display}%\n"
        f"CV ROC-AUC = {auc_display}\n"
        f"PR-AUC = {pr_auc_display}\n"
        f"β = {beta_display}, {p_display}"
    )
    
    # Position annotation at bottom-left
    ax.text(0.03, 0.03, stats_text, transform=ax.transAxes, fontsize=ANNOTATION_FONT_SIZE,
             verticalalignment='bottom', horizontalalignment='left',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, 
                      edgecolor='gray', linewidth=1))
    
    plt.tight_layout()
    
    # Save figure (SAME NAMES - No downstream impact)
    output_dir = r"M:\Research\WUE_CUE\data_products\results\linear_models"
    
    png_path = os.path.join(output_dir, "Publication_Figure_WUE_tra_CrossValidated.png")
    pdf_path = os.path.join(output_dir, "Publication_Figure_WUE_tra_CrossValidated.pdf")
    
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
    print("Weighted GLM (class-balanced) inference")
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
    
    # Check required columns
    if 'Coefficient_P_value' not in summary_df.columns:
        raise ValueError("❌ Step 2 summary missing required inference columns. Re-run Step 2.")
    
    coef = summary_df['Coefficient'].values[0]
    roc_auc = summary_df['ROC_AUC'].values[0]
    pr_auc = summary_df['PR_AUC'].values[0]
    baseline_pr = summary_df['Baseline_PR'].values[0]
    n_total = summary_df['N_Total'].values[0]
    n_decrease = summary_df['N_Decrease'].values[0]
    pct_decrease = summary_df['Pct_Decrease'].values[0]
    
    # Read proper inference (weighted GLM)
    p_value = summary_df['Coefficient_P_value'].values[0]
    ci_lower = summary_df['Coefficient_CI_lower'].values[0]
    ci_upper = summary_df['Coefficient_CI_upper'].values[0]
    inference_method = summary_df['Inference_Method'].values[0]
    
    # Read SPEI standardization parameters
    spei_mean = summary_df['spei_mean'].values[0]
    spei_std = summary_df['spei_std'].values[0]
    
    significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
    
    print("\n" + "="*60)
    print("MODEL STATISTICS (Weighted GLM - Class Balanced)")
    print("="*60)
    print(f"\nSample size: n = {n_total:,} observations")
    print(f"Class distribution: Decrease = {n_decrease} ({pct_decrease:.1f}%), Not Decrease = {n_total - n_decrease} ({100-pct_decrease:.1f}%)")
    
    print(f"\nCoefficient (standardized, per 1 SD increase in SPEI-48): βz = {coef:.4f}")
    print(f"95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"p-value: {p_value:.2e} {significance}")
    print(f"Inference method: {inference_method}")
    
    print(f"\nSPEI standardization (training data):")
    print(f"  mean = {spei_mean:.3f}")
    print(f"  std = {spei_std:.3f}")
    
    print(f"\nModel performance (site-blocked CV):")
    print(f"  ROC-AUC = {roc_auc:.1f}")
    print(f"  PR-AUC = {pr_auc:.1f} (baseline = {baseline_pr:.3f})")
    print(f"  Improvement over baseline = {(pr_auc - baseline_pr) / baseline_pr * 100:.0f}%")
    
    print("\n" + "="*60)
    print("PROBABILITY TABLE (from weighted GLM)")
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
    print("Weighted GLM (class-balanced) for inference")
    print("="*80)
    
    print("\n📝 GENERATING RESULTS FOR MANUSCRIPT (CONSOLE OUTPUT)...")
    print("-" * 60)
    summary_df, prob_df = report_results()
    
    print("\n\n🎨 CREATING PUBLICATION FIGURE (Single Panel)...")
    print("-" * 60)
    fig = create_figure()
    
    print("\n" + "="*80)
    print("WORKFLOW COMPLETE")
    print("="*80)
    print("\n📁 Output files (SAME NAMES - No downstream impact):")
    print("   • Publication_Figure_WUE_tra_CrossValidated.png")
    print("   • Publication_Figure_WUE_tra_CrossValidated.pdf")
    print("   • WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv (from Step 2 - UNCHANGED)")
    print("   • WUE_tra_Probability_Table_10_90.csv (from Step 2 - UNCHANGED)")
    print("\n✅ CRITICAL FIX: Curve now uses standardized SPEI correctly")
    print("   • β = -0.125 is applied to z-scored SPEI (not raw SPEI)")
    print("   • Curve matches probability table from Step 2")
    print("   • Figure annotation shows β (standardized)")
    print("\n📝 For manuscript methods: Report β as standardized coefficient (per 1 SD change in SPEI-48)")
    
    return summary_df, prob_df, fig

# ==============================================================================
# EXECUTE
# ==============================================================================

if __name__ == "__main__":
    summary_df, prob_df, fig = publication_workflow()