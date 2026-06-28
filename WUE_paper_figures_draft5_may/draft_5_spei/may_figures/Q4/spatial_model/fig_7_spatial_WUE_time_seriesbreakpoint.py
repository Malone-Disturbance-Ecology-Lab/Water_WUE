"""
FIGURE 2 — DATA PROCESSING & FIGURE CREATION (Complete)
Processes CSV data, computes yearly aggregates, and creates the figure.
UPDATED: WUE_T only, consistent with upstream workflow
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input file - from Concise Spatial Analyzer output
INPUT_CSV = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\tables_breakpoint\monthly_metric_long.csv"

# Output files for data
OUTPUT_STATS = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\tables_breakpoint\Figure2_yearly_statistics_WUE_T.csv"
OUTPUT_DATA = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\tables_breakpoint\Figure2_yearly_data_WUE_T.csv"

# Output files for figure
OUTPUT_PNG = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure2_yearly_timeseries_WUE_T.png"
OUTPUT_PDF = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure2_yearly_timeseries_WUE_T.pdf"

# Analysis parameters
METRIC_FILTER = "WUE_tra"  # Changed from "WUE" to "WUE_tra"
BREAKPOINT_YEAR = 2014
BREAKPOINT_MONTH = 8  # August
BREAKPOINT_X = BREAKPOINT_YEAR + (BREAKPOINT_MONTH - 1) / 12.0  # ≈2014.58

# Plot parameters
PLOT_PARAMS = {
    'figsize': (22, 18),
    'line_width': 4.0,
    'marker_size': 150,
    'marker_edge_width': 2.0,
    'frac_less_color': 'black',
    'spei_color': 'black',
    'breakpoint_color': 'red',
    'breakpoint_linewidth': 2.5,
    'breakpoint_linestyle': '--',
    'breakpoint_alpha': 0.8,
    'title_fontsize': 24,
    'axis_label_fontsize': 28,
    'tick_label_fontsize': 26,
    'annotation_fontsize': 18,
    'panel_label_fontsize': 22,
    'annotation_x': 0.03,
    'annotation_y': 0.95,
    'show_statistics': True,
    'show_post_trend': True,
    'x_tick_frequency': 2,
}

# ============================================================================
# PART 1: DATA PROCESSING
# ============================================================================

def load_and_process_data():
    """Load CSV data and process for yearly analysis."""
    print("📊 Loading data...")
    df = pd.read_csv(INPUT_CSV)
    
    # Filter for WUE_tra metric
    df_wue = df[df['metric'] == METRIC_FILTER].copy()
    print(f"  Found {len(df_wue)} rows for metric='{METRIC_FILTER}'")
    
    # Determine time column and create year field
    if 'date' in df_wue.columns:
        df_wue['date'] = pd.to_datetime(df_wue['date'])
        df_wue['year'] = df_wue['date'].dt.year
        print("  Using 'date' column for year extraction")
    elif 'yearmonth_int' in df_wue.columns:
        df_wue['year'] = df_wue['yearmonth_int'] // 100
        print("  Using 'yearmonth_int' column for year extraction")
    else:
        raise ValueError("No time column found. Need 'date' or 'yearmonth_int'")
    
    # Check required columns
    required_cols = ['frac_less', 'spei_mean']
    missing_cols = [col for col in required_cols if col not in df_wue.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    return df_wue

def compute_yearly_aggregates(df):
    """Compute yearly aggregates for plotting."""
    print("📈 Computing yearly aggregates...")
    
    yearly_data = df.groupby('year').agg({
        'frac_less': 'mean',
        'spei_mean': 'mean'
    }).reset_index()
    
    yearly_data = yearly_data.rename(columns={
        'frac_less': 'frac_less_year',
        'spei_mean': 'spei48_year'
    })
    
    yearly_data = yearly_data.sort_values('year')
    
    print(f"  Yearly data: {len(yearly_data)} years ({yearly_data['year'].min()} to {yearly_data['year'].max()})")
    
    return yearly_data

def compute_statistics(yearly_data):
    """Compute pre/post statistics and trends."""
    print("📊 Computing statistics...")
    
    yearly_data['period'] = np.where(yearly_data['year'] <= BREAKPOINT_YEAR, 'PRE', 'POST')
    
    pre_data = yearly_data[yearly_data['period'] == 'PRE']
    post_data = yearly_data[yearly_data['period'] == 'POST']
    
    stats_dict = {}
    
    def compute_var_stats(var_name, pre_series, post_series, post_only_data=None):
        mean_pre = pre_series.mean()
        mean_post = post_series.mean()
        delta_mean = mean_post - mean_pre
        
        if len(pre_series) > 0 and len(post_series) > 0:
            stat, p_value = stats.ranksums(pre_series, post_series)
        else:
            stat, p_value = np.nan, np.nan
        
        slope = np.nan
        p_slope = np.nan
        
        if post_only_data is not None and len(post_only_data) > 1:
            years = post_only_data['year'].values
            values = post_only_data[var_name].values
            
            if len(years) > 1:
                slope, intercept, r_value, p_slope, std_err = stats.linregress(years, values)
        
        return {
            'mean_pre': mean_pre,
            'mean_post': mean_post,
            'delta_mean': delta_mean,
            'p_value': p_value,
            'slope_per_year': slope,
            'p_slope': p_slope
        }
    
    stats_dict['frac_less'] = compute_var_stats(
        'frac_less_year', 
        pre_data['frac_less_year'], 
        post_data['frac_less_year'],
        post_data
    )
    
    stats_dict['spei48'] = compute_var_stats(
        'spei48_year', 
        pre_data['spei48_year'], 
        post_data['spei48_year']
    )
    
    summary_rows = []
    for var_name, var_stats in stats_dict.items():
        row = {
            'variable': var_name,
            'mean_pre': var_stats['mean_pre'],
            'mean_post': var_stats['mean_post'],
            'delta_mean': var_stats['delta_mean'],
            'p_value': var_stats['p_value'],
            'slope_per_year': var_stats['slope_per_year'],
            'p_slope': var_stats['p_slope']
        }
        summary_rows.append(row)
    
    stats_df = pd.DataFrame(summary_rows)
    
    return stats_df

def save_processed_data(yearly_data, stats_df):
    """Save processed data and statistics."""
    print("💾 Saving processed data...")
    
    os.makedirs(os.path.dirname(OUTPUT_DATA), exist_ok=True)
    
    yearly_data.to_csv(OUTPUT_DATA, index=False)
    print(f"  Yearly data saved: {OUTPUT_DATA}")
    
    stats_df.to_csv(OUTPUT_STATS, index=False)
    print(f"  Statistics saved: {OUTPUT_STATS}")
    
    return yearly_data, stats_df

# ============================================================================
# PART 2: FORMATTING FUNCTIONS FOR FIGURE
# ============================================================================

def format_p_value(p_val):
    """Format p-value for display."""
    if pd.isna(p_val):
        return "p = NA"
    elif p_val < 0.001:
        return "p < 0.001"
    elif p_val < 0.01:
        return f"p = {p_val:.3f}"
    elif p_val < 0.05:
        return f"p = {p_val:.3f}"
    else:
        return f"p = {p_val:.3f}"

def format_value(val, decimals=2):
    """Format value with specified decimals."""
    if pd.isna(val):
        return "NA"
    return f"{val:.{decimals}f}"

def add_breakpoint_line(ax, params):
    """Add breakpoint vertical line and label."""
    ax.axvline(x=BREAKPOINT_X, 
               color=params['breakpoint_color'],
               linewidth=params['breakpoint_linewidth'],
               linestyle=params['breakpoint_linestyle'],
               alpha=params['breakpoint_alpha'],
               zorder=1)
    
    label_y = ax.get_ylim()[1] * 0.95
    ax.text(BREAKPOINT_X + 0.1, label_y, 'Break point\n(Aug 2014)',
            color=params['breakpoint_color'],
            fontsize=params['annotation_fontsize'],
            fontweight='bold',
            ha='left',
            va='top',
            bbox=dict(boxstyle='round,pad=0.3',
                     facecolor='white',
                     alpha=0.9,
                     edgecolor=params['breakpoint_color'],
                     linewidth=1.5))

def add_statistics_annotation(ax, var_stats, params, is_spei=False):
    """Add statistics annotation to panel."""
    if not params['show_statistics']:
        return
    
    delta = var_stats['delta_mean']
    p_val = var_stats['p_value']
    
    delta_str = f"Δ = {format_value(delta, 2)}"
    p_str = format_p_value(p_val)
    
    if is_spei:
        text = f"{delta_str}\n{p_str}"
    else:
        slope = var_stats['slope_per_year']
        p_slope = var_stats['p_slope']
        
        if params['show_post_trend'] and not pd.isna(slope) and not pd.isna(p_slope):
            slope_str = f"β = {format_value(slope, 3)} yr⁻¹\n{format_p_value(p_slope)}"
            text = f"{delta_str}\n{p_str}\n{slope_str}"
        else:
            text = f"{delta_str}\n{p_str}"
    
    ax.text(params['annotation_x'], params['annotation_y'], text,
           transform=ax.transAxes,
           fontsize=params['annotation_fontsize'],
           fontweight='bold',
           verticalalignment='top',
           bbox=dict(boxstyle='round,pad=0.3',
                    facecolor='white',
                    alpha=0.9,
                    edgecolor='gray',
                    linewidth=1.0))

# ============================================================================
# PART 3: FIGURE CREATION
# ============================================================================

def create_figure2_from_data(yearly_data, stats_dict):
    """Create figure from already loaded data."""
    print("\n🎨 Creating figure...")
    params = PLOT_PARAMS.copy()
    
    fig, axes = plt.subplots(2, 1, figsize=params['figsize'], 
                           sharex=True, gridspec_kw={'hspace': 0.08})
    
    # ===== PANEL A: Frequency of WUE_T Decline =====
    ax_a = axes[0]
    
    ax_a.plot(yearly_data['year'], yearly_data['frac_less_year'],
             color=params['frac_less_color'],
             linewidth=params['line_width'],
             marker='o',
             markersize=params['marker_size']**0.5,
             markeredgecolor='white',
             markeredgewidth=params['marker_edge_width'],
             zorder=2,
             label='Frequency of WUE$_T$ Decline')
    
    add_breakpoint_line(ax_a, params)
    add_statistics_annotation(ax_a, stats_dict['frac_less'], params, is_spei=False)
    
    ax_a.set_ylabel('Frequency of\nWUE$_T$ Decline',
                   fontsize=params['axis_label_fontsize'],
                   fontweight='bold')
    ax_a.tick_params(axis='y', labelsize=params['tick_label_fontsize'])
    ax_a.tick_params(axis='x', labelbottom=False)
    ax_a.grid(True, alpha=0.3, linestyle='--')
    
    ax_a.text(0.99, 0.98, '(a)',
             transform=ax_a.transAxes,
             fontsize=params['panel_label_fontsize'],
             fontweight='bold',
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white',
                      alpha=0.95,
                      edgecolor='black',
                      linewidth=2.0))
    
    # ===== PANEL B: Annual mean SPEI =====
    ax_b = axes[1]
    
    ax_b.plot(yearly_data['year'], yearly_data['spei48_year'],
             color=params['spei_color'],
             linewidth=params['line_width'],
             marker='s', linestyle=':',
             markersize=params['marker_size']**0.5,
             markeredgecolor='white',
             markeredgewidth=params['marker_edge_width'],
             zorder=2,
             label='Annual mean SPEI')
    
    add_breakpoint_line(ax_b, params)
    add_statistics_annotation(ax_b, stats_dict['spei48'], params, is_spei=True)
    
    ax_b.set_ylabel('Annual mean SPEI',
                   fontsize=params['axis_label_fontsize'],
                   fontweight='bold')
    ax_b.set_xlabel('Year',
                   fontsize=params['axis_label_fontsize'],
                   fontweight='bold')
    ax_b.tick_params(axis='both', labelsize=params['tick_label_fontsize'])
    ax_b.grid(True, alpha=0.3, linestyle='--')
    
    ax_b.text(0.99, 0.98, '(b)',
             transform=ax_b.transAxes,
             fontsize=params['panel_label_fontsize'],
             fontweight='bold',
             verticalalignment='top',
             horizontalalignment='right',
             bbox=dict(boxstyle='round,pad=0.3',
                      facecolor='white',
                      alpha=0.95,
                      edgecolor='black',
                      linewidth=2.0))
    
    # Set x-axis ticks
    years = yearly_data['year'].values
    step = params['x_tick_frequency']
    tick_years = years[::step]
    
    if years[0] not in tick_years:
        tick_years = np.insert(tick_years, 0, years[0])
    
    if len(tick_years) > 0:
        if years[-1] - tick_years[-1] >= 1.5:
            tick_years = np.append(tick_years, years[-1])
    
    ax_b.set_xticks(tick_years)
    ax_b.set_xlim(years.min() - 0.5, years.max() + 1.25)
    plt.setp(ax_b.get_xticklabels(), rotation=0, ha='center')
    
    # Add figure title
    fig.suptitle('Temporal Trends in WUE$_T$ and Hydroclimate (2000-2025)',
                fontsize=params['title_fontsize'],
                fontweight='bold',
                y=0.94)
    
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    # Save figure
    print("\n💾 Saving figure...")
    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)
    
    fig.savefig(OUTPUT_PNG, dpi=400, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, dpi=400, bbox_inches='tight', facecolor='white')
    
    print(f"✅ PNG saved: {OUTPUT_PNG}")
    print(f"✅ PDF saved: {OUTPUT_PDF}")
    
    # Print final statistics
    print("\n📊 FINAL STATISTICS (WUE$_T$):")
    frac_stats = stats_dict['frac_less']
    spei_stats = stats_dict['spei48']
    
    print(f"  Panel A (Frequency of WUE$_T$ Decline):")
    print(f"    PRE mean: {format_value(frac_stats['mean_pre'], 2)}")
    print(f"    POST mean: {format_value(frac_stats['mean_post'], 2)}")
    print(f"    Δ (POST-PRE): {format_value(frac_stats['delta_mean'], 2)}")
    print(f"    p-value: {format_p_value(frac_stats['p_value'])}")
    
    if params['show_post_trend'] and not pd.isna(frac_stats['slope_per_year']):
        print(f"    Post-2014 slope: {format_value(frac_stats['slope_per_year'], 3)} per year ({format_p_value(frac_stats['p_slope'])})")
    
    print(f"\n  Panel B (Annual mean SPEI):")
    print(f"    PRE mean: {format_value(spei_stats['mean_pre'], 2)}")
    print(f"    POST mean: {format_value(spei_stats['mean_post'], 2)}")
    print(f"    Δ (POST-PRE): {format_value(spei_stats['delta_mean'], 2)}")
    print(f"    p-value: {format_p_value(spei_stats['p_value'])}")
    
    plt.show()
    
    return fig, axes

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def run_complete_workflow():
    """Run the complete Figure 2 workflow (data processing + figure creation)."""
    print("=" * 80)
    print("FIGURE 2: COMPLETE WORKFLOW (WUE$_T$)")
    print("=" * 80)
    
    try:
        # Step 1: Load and process data
        print("\n" + "="*60)
        print("STEP 1: PROCESSING DATA")
        print("="*60)
        df_wue = load_and_process_data()
        yearly_data = compute_yearly_aggregates(df_wue)
        stats_df = compute_statistics(yearly_data)
        
        # Step 2: Save processed data
        yearly_data, stats_df = save_processed_data(yearly_data, stats_df)
        
        # Convert stats_df to dictionary for easy access
        stats_dict = {}
        for _, row in stats_df.iterrows():
            stats_dict[row['variable']] = {
                'mean_pre': row['mean_pre'],
                'mean_post': row['mean_post'],
                'delta_mean': row['delta_mean'],
                'p_value': row['p_value'],
                'slope_per_year': row['slope_per_year'],
                'p_slope': row['p_slope']
            }
        
        print("\n" + "="*60)
        print("STEP 2: CREATING FIGURE")
        print("="*60)
        
        # Step 3: Create figure
        fig, axes = create_figure2_from_data(yearly_data, stats_dict)
        
        print("\n" + "="*80)
        print("✅ FIGURE 2 COMPLETE!")
        print("="*80)
        
        return yearly_data, stats_dict, fig, axes
        
    except Exception as e:
        print(f"❌ ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None, None

# ============================================================================
# RUN EVERYTHING
# ============================================================================

if __name__ == "__main__":
    yearly_data, stats_dict, fig, axes = run_complete_workflow()