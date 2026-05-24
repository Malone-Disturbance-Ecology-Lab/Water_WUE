"""
================================================================================
PUBLICATION FIGURES & RESULTS REPORTING - WUE LINEAR LOGISTIC REGRESSION
================================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import os

warnings.filterwarnings('ignore')

# ==============================================================================
# FIGURE PARAMETERS - ONLY LEGEND SIZE REDUCED
# ==============================================================================

FIGURE_WIDTH = 72  # Keep as is
FIGURE_HEIGHT = 32  # Keep as is

# Font sizes (ONLY ANNOTATION REDUCED by 20%)
YLABEL_FONT_SIZE = 106  # Keep as is
XLABEL_FONT_SIZE = 106  # Keep as is
XTICK_FONT_SIZE = 100   # Keep as is
YTICK_FONT_SIZE = 100   # Keep as is
LEGEND_FONT_SIZE = 88  # Keep as is
ANNOTATION_FONT_SIZE = 68  # REDUCED by 20% from 85 (85 * 0.8 = 68)

# Visual parameters (unchanged)
SITE_POINT_COLOR = '#606060'
SITE_LINE_COLOR = '#303030'
MIXED_LINE_COLOR = 'black'
AXIS_LINE_COLOR = '#000000'
AXIS_LINE_WIDTH = 6.0

SITE_LINE_ALPHA = 0.8
SITE_LINE_WIDTH = 6.0
MIXED_LINE_WIDTH = 9.0
SCATTER_ALPHA = 0.4
SCATTER_SIZE = 90

# Panel spacing
PANEL_SPACING = 0.35

# WUE labels - with subscripts (no dash)
WUE_LABELS = {
    "WUE": "WUE$_{ET}$",
    "WUE_eva": "WUE$_E$",  # Subscript E
    "WUE_tra": "WUE$_T$"   # Subscript T
}

# ==============================================================================
# PART 1: RESULTS REPORTING FUNCTIONS
# ==============================================================================

def report_detailed_results():
    """Print detailed statistical results for manuscript."""
    
    print("\n" + "="*80)
    print("RESULTS REPORTING - WUE LINEAR LOGISTIC REGRESSION")
    print("="*80)
    
    summary_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\linear_models\WUE_Linear_Model_Analysis_Summary.csv"
    
    try:
        summary_df = pd.read_csv(summary_path)
        print("📊 LOADED PREVIOUS MODEL RESULTS:")
        print("-" * 80)
        print(summary_df.to_string(index=False))
        print()
    except Exception as e:
        print(f"❌ Could not load summary: {e}")
        return None
    
    metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    
    print("\n📊 STATISTICAL SIGNIFICANCE FOR ALL METRICS:")
    print("-" * 80)
    
    model_stats = {}
    for metric in metrics:
        metric_data = summary_df[summary_df['Metric'] == metric]
        if len(metric_data) == 0:
            continue
            
        coef = metric_data['Linear_coef'].values[0]
        intercept = metric_data['Intercept'].values[0]
        n_obs = metric_data['N_raw_observations'].values[0]
        auc = metric_data['AUC'].values[0]
        auc_ci = metric_data['AUC_CI'].values[0]
        
        # Calculate p-value and confidence interval
        se = 1 / np.sqrt(n_obs)
        z_score = coef / se
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        
        # 95% CI for coefficient
        coef_ci_lower = coef - 1.96 * se
        coef_ci_upper = coef + 1.96 * se
        
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
        
        model_stats[metric] = {
            'coefficient': coef,
            'intercept': intercept,
            'p_value': p_value,
            'n_obs': n_obs,
            'auc': auc,
            'auc_ci': auc_ci,
            'coef_ci_lower': coef_ci_lower,
            'coef_ci_upper': coef_ci_upper,
            'significance': significance,
            'se': se
        }
        
        print(f"\n{WUE_LABELS[metric]}:")
        print(f"  Model: P(decrease) = 1 / (1 + exp(-({intercept:.2f} + {coef:.2f} × SPEI)))")
        print(f"  β₁ (coefficient): {coef:.2f} (95% CI: [{coef_ci_lower:.2f}, {coef_ci_upper:.2f}])")
        print(f"  p-value: {p_value:.6f} ({p_value:.3e}) {significance}")
        print(f"  AUC: {auc:.2f} ± {auc_ci:.2f}")
        print(f"  N: {n_obs:,}")
        print(f"  SE: {se:.2f}")
    
    return model_stats, summary_df

# ==============================================================================
# PART 2: PUBLICATION-READY FIGURES
# ==============================================================================

def create_publication_figures():
    """Create publication-ready figures with proper formatting."""
    
    print("\n" + "="*80)
    print("CREATING PUBLICATION-READY FIGURES")
    print("="*80)
    print("Figure Specifications:")
    print(f"- Size: {FIGURE_WIDTH} × {FIGURE_HEIGHT} inches")
    print(f"- Font sizes: Labels={XLABEL_FONT_SIZE}pt, Ticks={XTICK_FONT_SIZE}pt")
    print(f"- Annotation font: {ANNOTATION_FONT_SIZE}pt (20% smaller)")
    print(f"- WUE labels: {WUE_LABELS}")
    print("="*80)
    
    # Load data
    summary_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\linear_models\WUE_Linear_Model_Analysis_Summary.csv"
    
    try:
        summary_df = pd.read_csv(summary_path)
        print("✅ Loaded model summary")
    except Exception as e:
        print(f"❌ Could not load summary: {e}")
        return None
    
    data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
    
    try:
        df = pd.read_csv(data_path)
        print(f"✅ Loaded dataset: {df.shape}")
    except Exception as e:
        print(f"❌ Could not load data: {e}")
        return None
    
    # Enable LaTeX text rendering for proper beta symbol
    plt.rcParams.update({
        "text.usetex": False,  # Keep False for compatibility
        "mathtext.default": "regular",  # Use regular math text
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"]
    })
    
    # Set clean style
    plt.style.use('default')
    sns.set_style("whitegrid", {
        'grid.linestyle': ':', 
        'grid.alpha': 0.2,
        'axes.linewidth': AXIS_LINE_WIDTH,
        'axes.edgecolor': AXIS_LINE_COLOR
    })
    
    # Apply font sizes
    plt.rcParams.update({
        'font.size': XTICK_FONT_SIZE,
        'axes.labelsize': XLABEL_FONT_SIZE,
        'axes.titlesize': LEGEND_FONT_SIZE,
        'xtick.labelsize': XTICK_FONT_SIZE,
        'ytick.labelsize': YTICK_FONT_SIZE,
        'figure.titlesize': LEGEND_FONT_SIZE + 8,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })
    
    metrics = ['WUE', 'WUE_eva', 'WUE_tra']
    
    # Calculate model statistics
    print("\n📊 CALCULATING MODEL STATISTICS FOR FIGURE ANNOTATIONS:")
    print("-" * 80)
    model_stats = {}
    for metric in metrics:
        metric_data = summary_df[summary_df['Metric'] == metric]
        if len(metric_data) > 0:
            coef = metric_data['Linear_coef'].values[0]
            n_obs = metric_data['N_raw_observations'].values[0]
            se = 1 / np.sqrt(n_obs)
            z_score = coef / se
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
            
            # 95% CI for coefficient
            coef_ci_lower = coef - 1.96 * se
            coef_ci_upper = coef + 1.96 * se
            
            significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""
            
            model_stats[metric] = {
                'coefficient': coef,
                'p_value': p_value,
                'n_obs': n_obs,
                'coef_ci_lower': coef_ci_lower,
                'coef_ci_upper': coef_ci_upper,
                'significance': significance
            }
            
            print(f"{WUE_LABELS[metric]}:")
            print(f"  β₁ = {coef:.3f}, 95% CI = [{coef_ci_lower:.2f}, {coef_ci_upper:.2f}]")
            print(f"  p = {p_value:.2e} {significance}")
    
    # ==========================================================================
    # FIGURE 1: ALL OBSERVATIONS
    # ==========================================================================
    
    print("\n📊 Creating Figure 1: All observations...")
    
    fig1, axes1 = plt.subplots(1, 3, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig1.subplots_adjust(wspace=PANEL_SPACING)
    
    for idx, metric in enumerate(metrics):
        ax = axes1[idx]
        
        metric_data = summary_df[summary_df['Metric'] == metric]
        if len(metric_data) == 0:
            ax.text(0.5, 0.5, f"No data\nfor {metric}", 
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=LEGEND_FONT_SIZE)
            continue
        
        intercept = metric_data['Intercept'].values[0]
        coef = metric_data['Linear_coef'].values[0]
        
        column_map = {
            'WUE': 'WUE_minus_NNmed_SPEI_48_Class',
            'WUE_eva': 'WUE_eva_minus_NNmed_SPEI_48_Class',
            'WUE_tra': 'WUE_tra_minus_NNmed_SPEI_48_Class'
        }
        
        class_col = column_map[metric]
        spei_col = 'SPEI_48'
        
        mask = (df[class_col].isin(['Increase', 'Decrease'])) & (df[spei_col].notna())
        filtered_df = df[mask].copy()
        filtered_df['target'] = (filtered_df[class_col] == 'Decrease').astype(int)
        
        X_raw = filtered_df[spei_col].values
        n_samples = len(X_raw)
        
        # Create smooth curve
        spei_range = np.linspace(X_raw.min() - 0.5, X_raw.max() + 0.5, 200)
        log_odds = intercept + coef * spei_range
        probabilities = 1 / (1 + np.exp(-log_odds))
        
        # Plot BLACK fitted curve
        ax.plot(spei_range, probabilities, 
                color=MIXED_LINE_COLOR, 
                linewidth=MIXED_LINE_WIDTH,
                label=WUE_LABELS[metric])
        
        # Plot GREY jittered points
        raw_log_odds = intercept + coef * X_raw
        raw_probs = 1 / (1 + np.exp(-raw_log_odds))
        
        jitter_x = np.random.normal(0, 0.015, len(X_raw))
        jitter_y = np.random.normal(0, 0.008, len(raw_probs))
        
        ax.scatter(X_raw + jitter_x, raw_probs + jitter_y,
                  color=SITE_POINT_COLOR,
                  alpha=SCATTER_ALPHA,
                  s=SCATTER_SIZE,
                  edgecolor='none',
                  rasterized=True)
        
        # Add reference lines
        ax.axvspan(-0.5, 0.5, alpha=0.08, color='gray', zorder=0)
        ax.axhline(0.5, color='black', linestyle='-', alpha=0.3, linewidth=3.0)
        
        # Axis labels with HUGE font
        ax.set_xlabel('SPEI', fontsize=XLABEL_FONT_SIZE, fontweight='bold')
        if idx == 0:
            ax.set_ylabel('P(WUE decrease)', fontsize=YLABEL_FONT_SIZE, fontweight='bold')
        
        # WUE label as title (top center) - with HUGE font and subscript
        ax.set_title(WUE_LABELS[metric], fontsize=LEGEND_FONT_SIZE, pad=30, fontweight='bold')        
# Add panel label (a), (b), (c) in top left for Figure 2
        panel_labels = ['(a)', '(b)', '(c)']
        ax.text(0.5, 1.15, panel_labels[idx], 
                transform=ax.transAxes,
                fontsize=ANNOTATION_FONT_SIZE,  # Use annotation font size
                fontweight='bold',
                va='top',
                ha='left')


        
        # Annotations in top left - Multi-line with statistics
        if metric in model_stats:
            stats_info = model_stats[metric]
            
            # Format p-value in scientific notation
            p_val = stats_info['p_value']
            if p_val < 0.001:
                p_text = f"p < 0.001"
            else:
                p_text = f"p = {p_val:.2e}"
            
            # Create annotation text with FIXED beta symbol using Unicode
            annotation_text = f"n = {stats_info['n_obs']:,}\n"
         # Unicode beta with subscript 1
            annotation_text += r"$\beta_1$ = " + f"{stats_info['coefficient']:.3f}\n"
            annotation_text += f"95% CI = [{stats_info['coef_ci_lower']:.3f}, {stats_info['coef_ci_upper']:.3f}]\n"
            annotation_text += f"{p_text} {stats_info['significance']}"
            
            ax.text(0.04, 0.96, annotation_text,
                    transform=ax.transAxes,
                    fontsize=ANNOTATION_FONT_SIZE,  # USING REDUCED SIZE
                    verticalalignment='top',
                    linespacing=1.4,
                    bbox=dict(boxstyle='round,pad=1.0', facecolor='white', alpha=0.95,
                             edgecolor='black', linewidth=3.0))
            
            
        # Set limits
        ax.set_ylim([-0.05, 1.05])
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_xlim([X_raw.min() - 0.3, X_raw.max() + 0.3])
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(5))  # Set exactly 5 ticks
        
        # Customize ticks with HUGE font
        ax.tick_params(axis='both', which='major', 
                      labelsize=XTICK_FONT_SIZE,
                      width=AXIS_LINE_WIDTH * 0.8,
                      length=AXIS_LINE_WIDTH * 4)
        
        # Set axis line width
        for spine in ax.spines.values():
            spine.set_linewidth(AXIS_LINE_WIDTH)
    
    plt.tight_layout()
    
    # Save Figure 1 as BOTH PDF and PNG
    output_dir = os.path.dirname(summary_path)
    
    fig1_path_pdf = os.path.join(output_dir, "Publication_Figure1_WUE_Probability_All_Observations.pdf")
    fig1_path_png = os.path.join(output_dir, "Publication_Figure1_WUE_Probability_All_Observations.png")
    
    fig1.savefig(fig1_path_pdf, dpi=300, bbox_inches='tight')
    fig1.savefig(fig1_path_png, dpi=300, bbox_inches='tight')
    print(f"✅ Saved Figure 1 PDF: {os.path.basename(fig1_path_pdf)}")
    print(f"✅ Saved Figure 1 PNG: {os.path.basename(fig1_path_png)}")
    plt.show()
    
    # ==========================================================================
    # FIGURE 2: UNIQUE SPEI VALUES
    # ==========================================================================
    
    print("\n📊 Creating Figure 2: Unique SPEI values...")
    
    fig2, axes2 = plt.subplots(1, 3, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig2.subplots_adjust(wspace=PANEL_SPACING)
    
    for idx, metric in enumerate(metrics):
        ax = axes2[idx]
        
        metric_data = summary_df[summary_df['Metric'] == metric]
        if len(metric_data) == 0:
            ax.text(0.5, 0.5, f"No data\nfor {metric}", 
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=LEGEND_FONT_SIZE)
            continue
        
        intercept = metric_data['Intercept'].values[0]
        coef = metric_data['Linear_coef'].values[0]
        
        column_map = {
            'WUE': 'WUE_minus_NNmed_SPEI_48_Class',
            'WUE_eva': 'WUE_eva_minus_NNmed_SPEI_48_Class',
            'WUE_tra': 'WUE_tra_minus_NNmed_SPEI_48_Class'
        }
        
        class_col = column_map[metric]
        spei_col = 'SPEI_48'
        
        mask = (df[class_col].isin(['Increase', 'Decrease'])) & (df[spei_col].notna())
        filtered_df = df[mask].copy()
        filtered_df['target'] = (filtered_df[class_col] == 'Decrease').astype(int)
        
        # Get unique SPEI values
        unique_spei = np.unique(filtered_df[spei_col].values)
        
        # Calculate probabilities
        log_odds_unique = intercept + coef * unique_spei
        probs_unique = 1 / (1 + np.exp(-log_odds_unique))
        
        # Plot BLACK fitted curve
        sort_idx = np.argsort(unique_spei)
        unique_spei_sorted = unique_spei[sort_idx]
        probs_unique_sorted = probs_unique[sort_idx]
        
        ax.plot(unique_spei_sorted, probs_unique_sorted, 
                color=MIXED_LINE_COLOR, 
                linewidth=MIXED_LINE_WIDTH,
                label=WUE_LABELS[metric])
        
        # Plot GREY points
        ax.scatter(unique_spei, probs_unique,
                  facecolor=SITE_POINT_COLOR,
                  edgecolor=SITE_POINT_COLOR,
                  linewidth=3.0,
                  s=SCATTER_SIZE * 3.0,  # Much larger points
                  marker='o',
                  alpha=0.7,
                  zorder=1)
        
        # Add reference lines
        ax.axvspan(-0.5, 0.5, alpha=0.08, color='gray', zorder=0)
        ax.axhline(0.5, color='black', linestyle='-', alpha=0.3, linewidth=3.0)
        
        # Axis labels with HUGE font
        ax.set_xlabel('SPEI', fontsize=XLABEL_FONT_SIZE, fontweight='bold')
        if idx == 0:
            ax.set_ylabel('P(WUE decrease)', fontsize=YLABEL_FONT_SIZE, fontweight='bold')
        
        # WUE label as title - with HUGE font and subscript
        ax.set_title(WUE_LABELS[metric], fontsize=LEGEND_FONT_SIZE, pad=30, fontweight='bold')
        
        # Add panel label (a), (b), (c) in top left
        panel_labels = ['(a)', '(b)', '(c)']
        ax.text(-0.2, 0.98, panel_labels[idx], #0.5, 1.15
                transform=ax.transAxes,
                fontsize=ANNOTATION_FONT_SIZE+10,  # Use annotation font size
                fontweight='bold',
                va='bottom',
                ha='left')

        
        # Annotations in top left - Multi-line with statistics
        if metric in model_stats:
            stats_info = model_stats[metric]
            
            # Format p-value in scientific notation
            p_val = stats_info['p_value']
            if p_val < 0.001:
                p_text = f"p < 0.001"
            else:
                p_text = f"p = {p_val:.2e}"
            
            # Create annotation text with FIXED beta symbol using Unicode
            annotation_text = f"n = {len(unique_spei)}\n"
            annotation_text += r"$\beta_1$ = " + f"{stats_info['coefficient']:.2f}\n" # Unicode beta with subscript 1
            annotation_text += f"95% CI = [{stats_info['coef_ci_lower']:.2f}, {stats_info['coef_ci_upper']:.2f}]\n"
            annotation_text += f"{p_text} {stats_info['significance']}"
            
            ax.text(0.04, 0.96, annotation_text,
                    transform=ax.transAxes,
                    fontsize=ANNOTATION_FONT_SIZE,  # USING REDUCED SIZE
                    verticalalignment='top',
                    linespacing=1.4,
                    bbox=dict(boxstyle='round,pad=1.0', facecolor='white', alpha=0.95,
                             edgecolor='black', linewidth=3.0))
        
        
        # Set limits
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlim([unique_spei.min() - 0.3, unique_spei.max() + 0.3])
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))  # Set exactly 5 ticks
        
        # Customize ticks with HUGE font
        ax.tick_params(axis='both', which='major', 
                      labelsize=XTICK_FONT_SIZE,
                      width=AXIS_LINE_WIDTH * 0.8,
                      length=AXIS_LINE_WIDTH * 4)
        
        # Set axis line width
        for spine in ax.spines.values():
            spine.set_linewidth(AXIS_LINE_WIDTH)
    
    plt.tight_layout()
    
    # Save Figure 2 as BOTH PDF and PNG
    fig2_path_pdf = os.path.join(output_dir, "Publication_Figure2_WUE_Probability_Unique_SPEI.pdf")
    fig2_path_png = os.path.join(output_dir, "Publication_Figure2_WUE_Probability_Unique_SPEI.png")
    
    fig2.savefig(fig2_path_pdf, dpi=300, bbox_inches='tight')
    fig2.savefig(fig2_path_png, dpi=300, bbox_inches='tight')
    print(f"✅ Saved Figure 2 PDF: {os.path.basename(fig2_path_pdf)}")
    print(f"✅ Saved Figure 2 PNG: {os.path.basename(fig2_path_png)}")
    plt.show()
    
    print("\n" + "="*80)
    print("PUBLICATION FIGURES COMPLETE")
    print("="*80)
    print(f"\n📏 Figure dimensions: {FIGURE_WIDTH} × {FIGURE_HEIGHT} inches")
    print(f"📁 Generated files:")
    print(f"  • Publication_Figure1_WUE_Probability_All_Observations.pdf")
    print(f"  • Publication_Figure1_WUE_Probability_All_Observations.png")
    print(f"  • Publication_Figure2_WUE_Probability_Unique_SPEI.pdf")
    print(f"  • Publication_Figure2_WUE_Probability_Unique_SPEI.png")
    
    return fig1, fig2

# ==============================================================================
# MAIN WORKFLOW
# ==============================================================================

def publication_workflow():
    """Complete publication workflow."""
    
    print("\n" + "="*80)
    print("PUBLICATION WORKFLOW - WUE LINEAR LOGISTIC REGRESSION")
    print("="*80)
    
    print("\n📝 PART 1: GENERATING RESULTS FOR MANUSCRIPT...")
    print("-" * 80)
    model_stats, summary_df = report_detailed_results()
    
    print("\n\n🎨 PART 2: CREATING PUBLICATION-READY FIGURES...")
    print("-" * 80)
    figures = create_publication_figures()
    
    print("\n" + "="*80)
    print("WORKFLOW COMPLETE")
    print("="*80)
    
    return model_stats, figures

# ==============================================================================
# EXECUTE
# ==============================================================================

if __name__ == "__main__":
    publication_workflow()

# To run manually:
# model_stats, figures = publication_workflow()