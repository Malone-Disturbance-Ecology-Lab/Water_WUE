# -*- coding: utf-8 -*-
"""
Created on Thu May 14 10:25:53 2026

@author: ammar
"""

"""
MERGED FIGURE: Breakpoint Analysis with SPEI-48 Gradient Lines (FINAL CORRECTED)
===============================================================================
- Left axis: p(WUE_T Decline) frequency
- PRE SPEI gradient: LEFT of breakpoint line
- POST SPEI gradient: RIGHT of breakpoint line
- Colormap: FireBrick → Light Pink → Light Blue → Deep Navy
- Min/Max SPEI labels at ends of each gradient line
- Single global colorbar at bottom
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
matplotlib.rcParams['font.weight'] = 'bold'

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
WUE_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\regional_analysis\all_regions_monthly_frac_less.csv"
SPEI_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\regional_analysis\all_regions_monthly_spei48.csv"

# Output
OUTPUT_PNG = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure_Merged_Breakpoint_SPEI.png"
OUTPUT_PDF = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\Figure_Merged_Breakpoint_SPEI.pdf"

# Region settings
REGION_COLORS = {
    "Pacific": "#E67E22",
    "Gulf": "#F1C40F",
    "Southeast Atlantic": "#9B59B6",
    "Atlantic North": "#3498DB",
    "Alaska": "#1ABC9C"
}

REGION_ORDER = ["Pacific", "Gulf", "Southeast Atlantic", "Atlantic North", "Alaska"]

# CORRECT SPEI COLORMAP (from your sample)
SPEI_CMAP = LinearSegmentedColormap.from_list(
    "enhanced_gradient",
    [
        "#B22222",  # FireBrick - very dry (< -2)
        "#DC143C",  # Crimson - dry
        "#F08080",  # LightCoral - moderately dry
        "#FADADD",  # Very light pink - near zero negative
        "#D4E6F1",  # Very light blue - near zero positive
        "#5DADE2",  # Light blue - moderately wet
        "#2874A6",  # Medium blue - wet
        "#1A5276",  # Dark blue - very wet (> 2)
        "#0B3B60"   # Deep Navy - extreme wet
    ]
)

# ============================================================================
# PETTITT TEST ON SPEI
# ============================================================================

def pettitt_test_spei(spei_series):
    """Pettitt test for change point detection on SPEI time series."""
    n = len(spei_series)
    if n < 10:
        return None, None
    
    valid_idx = np.isfinite(spei_series)
    if np.sum(valid_idx) < 10:
        return None, None
    
    series_clean = spei_series[valid_idx]
    indices_clean = np.where(valid_idx)[0]
    n_clean = len(series_clean)
    
    U = np.zeros(n_clean)
    for t in range(1, n_clean):
        for i in range(t):
            for j in range(t, n_clean):
                if series_clean[i] < series_clean[j]:
                    U[t] += 1
                elif series_clean[i] > series_clean[j]:
                    U[t] -= 1
    
    K = np.max(np.abs(U))
    tau = np.argmax(np.abs(U))
    p_value = 2 * np.exp(-6 * K**2 / (n_clean**3 + n_clean**2))
    breakpoint_idx = indices_clean[tau] if tau < len(indices_clean) else None
    
    return breakpoint_idx, p_value

# ============================================================================
# DATA LOADING AND PROCESSING
# ============================================================================

def load_and_process_data():
    """Load WUE and SPEI data."""
    df_wue = pd.read_csv(WUE_FILE)
    df_wue['date'] = pd.to_datetime(df_wue['date'])
    
    df_spei = pd.read_csv(SPEI_FILE)
    df_spei['date'] = pd.to_datetime(df_spei['date'])
    
    print("="*80)
    print("MERGED FIGURE: Breakpoint Analysis + SPEI-48 Gradient Lines")
    print("="*80)
    print(f"\n📊 WUE data: {len(df_wue)} records")
    print(f"📊 SPEI data: {len(df_spei)} records")
    print("\n⚠️ Breakpoint calculated from SPEI-48 (not WUE)")
    
    return df_wue, df_spei

def analyze_region(df_wue, df_spei, region_name, color):
    """Analyze one region: SPEI breakpoint + WUE statistics."""
    
    print(f"\n📈 Analyzing {region_name}...")
    
    # Get SPEI data first (for breakpoint calculation)
    df_spei_region = df_spei[df_spei['region'] == region_name].sort_values('date').reset_index(drop=True)
    spei_values = df_spei_region['mean_SPEI48'].values
    spei_dates = pd.to_datetime(df_spei_region['date'].values)
    
    # Calculate breakpoint from SPEI
    breakpoint_idx, p_value = pettitt_test_spei(spei_values)
    
    if breakpoint_idx is not None:
        breakpoint_date = spei_dates[breakpoint_idx]
        breakpoint_year = breakpoint_date.year
        breakpoint_str = breakpoint_date.strftime('%b %Y')
        
        # Split SPEI into PRE and POST
        pre_spei = spei_values[:breakpoint_idx]
        post_spei = spei_values[breakpoint_idx:]
        pre_spei_mean = np.mean(pre_spei)
        post_spei_mean = np.mean(post_spei)
        spei_delta = post_spei_mean - pre_spei_mean
        
        # SPEI min/max for gradient lines
        spei_min_pre = np.min(pre_spei)
        spei_max_pre = np.max(pre_spei)
        spei_min_post = np.min(post_spei)
        spei_max_post = np.max(post_spei)
        
        # Get WUE data (aligned with same breakpoint)
        df_wue_region = df_wue[df_wue['region'] == region_name].sort_values('date').reset_index(drop=True)
        wue_dates = pd.to_datetime(df_wue_region['date'].values)
        wue_values = df_wue_region['frac_less'].values
        
        # Find WUE values at same time points
        pre_wue = wue_values[wue_dates < breakpoint_date]
        post_wue = wue_values[wue_dates >= breakpoint_date]
        pre_wue_mean = np.mean(pre_wue) if len(pre_wue) > 0 else np.nan
        post_wue_mean = np.mean(post_wue) if len(post_wue) > 0 else np.nan
        wue_delta = post_wue_mean - pre_wue_mean
        wue_pct = (wue_delta / pre_wue_mean * 100) if pre_wue_mean != 0 else np.nan
    else:
        breakpoint_date = None
        breakpoint_year = None
        breakpoint_str = "None"
        p_value = np.nan
        pre_spei = []
        post_spei = []
        pre_spei_mean = np.nan
        post_spei_mean = np.nan
        spei_delta = np.nan
        spei_min_pre = np.nan
        spei_max_pre = np.nan
        spei_min_post = np.nan
        spei_max_post = np.nan
        pre_wue_mean = np.nan
        post_wue_mean = np.nan
        wue_delta = np.nan
        wue_pct = np.nan
        wue_values = []
        wue_dates = []
    
    print(f"  Breakpoint (from SPEI): {breakpoint_str} (p={p_value:.4f})")
    print(f"  SPEI: {pre_spei_mean:.3f} → {post_spei_mean:.3f} (Δ={spei_delta:+.3f})")
    print(f"  p(WUE_T Decline): {pre_wue_mean:.3f} → {post_wue_mean:.3f} (Δ={wue_delta:+.3f})")
    
    # Convert WUE to yearly for plotting
    if len(wue_values) > 0:
        df_yearly = pd.DataFrame({'date': wue_dates, 'frac_less': wue_values})
        df_yearly['year'] = df_yearly['date'].dt.year
        yearly = df_yearly.groupby('year')['frac_less'].mean().reset_index()
        years = yearly['year'].values
        yearly_wue = yearly['frac_less'].values
    else:
        years = np.array([])
        yearly_wue = np.array([])
    
    return {
        'region': region_name,
        'color': color,
        'wue_dates': wue_dates,
        'wue_values': wue_values,
        'spei_dates': spei_dates,
        'spei_values': spei_values,
        'breakpoint_idx': breakpoint_idx,
        'breakpoint_date': breakpoint_date,
        'breakpoint_year': breakpoint_year,
        'breakpoint_str': breakpoint_str,
        'p_value': p_value,
        'pre_wue_mean': pre_wue_mean,
        'post_wue_mean': post_wue_mean,
        'wue_delta': wue_delta,
        'wue_pct': wue_pct,
        'pre_spei_mean': pre_spei_mean,
        'post_spei_mean': post_spei_mean,
        'spei_delta': spei_delta,
        'spei_min_pre': spei_min_pre,
        'spei_max_pre': spei_max_pre,
        'spei_min_post': spei_min_post,
        'spei_max_post': spei_max_post,
        'years': years,
        'yearly_wue': yearly_wue
    }

# ============================================================================
# CREATE MERGED FIGURE
# ============================================================================

def create_merged_figure(results):
    """Create figure with WUE_T (left axis) and PRE/POST SPEI gradient lines."""
    
    n_regions = len(results)
    fig, axes = plt.subplots(n_regions, 1, figsize=(16, 3 * n_regions), sharex=False)
    
    if n_regions == 1:
        axes = [axes]
    
    panel_labels = ['a', 'b', 'c', 'd', 'e']
    
    # Global SPEI range for color normalization
    all_spei = []
    for r in results:
        all_spei.extend(r['spei_values'])
    global_spei_min, global_spei_max = np.min(all_spei), np.max(all_spei)
    spei_norm = Normalize(vmin=global_spei_min, vmax=global_spei_max)
    
    for i, result in enumerate(results):
        ax = axes[i]
        region = result['region']
        color = REGION_COLORS[region]
        
        # ========================================
        # LEFT AXIS: p(WUE_T Decline) Frequency
        # ========================================
        if len(result['years']) > 0:
            ax.plot(result['years'], result['yearly_wue'], 
                    color='black', linewidth=2.5, marker='o', markersize=7, 
                    alpha=0.85, zorder=2)
        
        # Breakpoint vertical line (from SPEI)
        if result['breakpoint_year'] is not None:
            ax.axvline(x=result['breakpoint_year'], color='red', linestyle='--', 
                      linewidth=2.5, alpha=0.9, zorder=3)
        
        # PRE mean line (p(WUE_T Decline))
        if not np.isnan(result['pre_wue_mean']) and len(result['years']) > 0:
            start_year = result['years'][0]
            bp_year = result['breakpoint_year']
            if bp_year is not None:
                ax.hlines(y=result['pre_wue_mean'], xmin=start_year, xmax=bp_year,
                         color='black', linestyle=':', linewidth=2, alpha=0.7)
        
        # POST mean line (p(WUE_T Decline))
        if not np.isnan(result['post_wue_mean']) and len(result['years']) > 0:
            end_year = result['years'][-1]
            bp_year = result['breakpoint_year']
            if bp_year is not None:
                ax.hlines(y=result['post_wue_mean'], xmin=bp_year, xmax=end_year,
                         color=color, linestyle=':', linewidth=2, alpha=0.7)
        
        # Left axis labels
        ax.set_ylabel('p(WUE$_T$ Decline)', fontsize=13, fontweight='bold', color='#2C3E50')
        ax.tick_params(axis='y', labelcolor='#2C3E50', labelsize=11)
        ax.set_ylim(0, 1)
        
        # ========================================
        # SPEI GRADIENT LINES (Side by side with breakpoint in between)
        # ========================================
        if result['breakpoint_year'] is not None:
            bp_year = result['breakpoint_year']
            y_min, y_max = ax.get_ylim()
            
            # Position for gradient lines (near bottom, above x-axis)
            y_position = 0.75  # 8% up from bottom
            spei_line_y = y_min + y_position * (y_max - y_min)
            
            # Line height/thickness
            line_height = 0.1 * (y_max - y_min)
            
            # PRE gradient line (LEFT of breakpoint)
            if not np.isnan(result['spei_min_pre']):
                x_start_pre = result['years'][0]
                x_end_pre = bp_year - 0.5  # Stop before breakpoint
                x_center_pre = (x_start_pre + x_end_pre) / 2
                line_width_pre = (x_end_pre - x_start_pre) * 0.8
                x_start_draw_pre = x_center_pre - line_width_pre/2
                
                # Draw PRE gradient line
                spei_range_pre = result['spei_max_pre'] - result['spei_min_pre']
                if spei_range_pre == 0:
                    spei_range_pre = 0.001
                
                n_segments = 200
                for seg in range(n_segments):
                    t = seg / n_segments
                    spei_value = result['spei_min_pre'] + t * spei_range_pre
                    line_color = SPEI_CMAP(spei_norm(spei_value))
                    
                    x_seg_start = x_start_draw_pre + t * line_width_pre
                    x_seg_end = x_start_draw_pre + (seg + 1) / n_segments * line_width_pre
                    
                    ax.hlines(y=spei_line_y, xmin=x_seg_start, xmax=x_seg_end,
                             color=line_color, linewidth=line_height*40, zorder=10, alpha=0.95)
                
                # PRE labels
                pre_min_color = SPEI_CMAP(spei_norm(result['spei_min_pre']))
                pre_max_color = SPEI_CMAP(spei_norm(result['spei_max_pre']))
                
                # Min label (left end)
                ax.text(x_start_draw_pre - 0.3, spei_line_y, f'{result["spei_min_pre"]:.1f}', 
                       ha='right', va='center', fontsize=15, fontweight='bold',
                       color=pre_min_color, zorder=11)
                
                # Max label (right end, near breakpoint)
                ax.text(x_start_draw_pre + line_width_pre + 0.1, spei_line_y, f'{result["spei_max_pre"]:.1f}', 
                       ha='left', va='center', fontsize=15, fontweight='bold',
                       color=pre_max_color, zorder=11)
                
                # PRE label text
                ax.text(x_center_pre, spei_line_y - 0.02, 'PRE', 
                       ha='center', va='top', fontsize=9, fontweight='bold',
                       color='gray', style='italic', zorder=11)
            
            # POST gradient line (RIGHT of breakpoint)
            if not np.isnan(result['spei_min_post']):
                x_start_post = bp_year + 0.5  # Start after breakpoint
                x_end_post = result['years'][-1]
                x_center_post = (x_start_post + x_end_post) / 2
                line_width_post = (x_end_post - x_start_post) * 0.8
                x_start_draw_post = x_center_post - line_width_post/2
                
                # Draw POST gradient line
                spei_range_post = result['spei_max_post'] - result['spei_min_post']
                if spei_range_post == 0:
                    spei_range_post = 0.001
                
                n_segments = 200
                for seg in range(n_segments):
                    t = seg / n_segments
                    spei_value = result['spei_min_post'] + t * spei_range_post
                    line_color = SPEI_CMAP(spei_norm(spei_value))
                    
                    x_seg_start = x_start_draw_post + t * line_width_post
                    x_seg_end = x_start_draw_post + (seg + 1) / n_segments * line_width_post
                    
                    ax.hlines(y=spei_line_y, xmin=x_seg_start, xmax=x_seg_end,
                             color=line_color, linewidth=line_height*40, zorder=10, alpha=0.95)
                
                # POST labels
                post_min_color = SPEI_CMAP(spei_norm(result['spei_min_post']))
                post_max_color = SPEI_CMAP(spei_norm(result['spei_max_post']))
                
                # Min label (left end, near breakpoint)
                ax.text(x_start_draw_post - 0.1, spei_line_y, f'{result["spei_min_post"]:.1f}', 
                       ha='right', va='center', fontsize=15, fontweight='bold',
                       color=post_min_color, zorder=11)
                
                # Max label (right end)
                ax.text(x_start_draw_post + line_width_post + 0.3, spei_line_y, f'{result["spei_max_post"]:.1f}', 
                       ha='left', va='center', fontsize=15, fontweight='bold',
                       color=post_max_color, zorder=11)
                
                # POST label text
                ax.text(x_center_post, spei_line_y - 0.02, 'POST', 
                       ha='center', va='top', fontsize=15, fontweight='bold',
                       color='gray', style='italic', zorder=11)
        
        # ========================================
        # CLEAN TEXT BOX (essential info only)
        # ========================================
        p_text = f"p = {result['p_value']:.3f}" if not np.isnan(result['p_value']) else "p = N/A"
        if not np.isnan(result['p_value']) and result['p_value'] < 0.001:
            p_text = "p < 0.001"
        
        stat_text = (
            f"Breakpoint: {result['breakpoint_str']}\n"
            f"{p_text}\n"
            f"Δp(WUE$_T$): {result['wue_delta']:+.2f} ({result['wue_pct']:+.0f}%)\n"
            f"ΔSPEI: {result['spei_delta']:+.2f}"
        )
        
        ax.text(0.02, 0.02, stat_text, transform=ax.transAxes,
               fontsize=13, fontweight='bold', verticalalignment='bottom',
               bbox=dict(boxstyle='round,pad=0.6', facecolor='white', alpha=0.85,
                        edgecolor=color, linewidth=1.5), zorder=20)
        
        # Panel label
        ax.text(0.02, 0.96, f'({panel_labels[i]})', transform=ax.transAxes,
               fontsize=14, fontweight='bold', verticalalignment='top',
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85,
                        edgecolor='black', linewidth=1), zorder=20)
        
        # Region title
        ax.set_title(f'{region}', fontsize=16, fontweight='bold', color=color, pad=12)
        
        # Grid and formatting
        ax.tick_params(axis='both', labelsize=10, width=1.5, length=5)
        ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        ax.set_xlim(1999.5, 2026)
        
        # Simple legend for WUE lines (only on first panel)
        if i == 0:
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color='black', linestyle=':', linewidth=2, label='p(WUE$_T$) PRE mean'),
                Line2D([0], [0], color=color, linestyle=':', linewidth=2, label='p(WUE$_T$) POST mean'),
                Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='Breakpoint (SPEI)')
            ]
            ax.legend(handles=legend_elements, loc='lower right', fontsize=13,
                     framealpha=0.9, edgecolor='gray', handlelength=1.5)
    
    # Format x-axis for bottom panel only
    #axes[-1].set_xlabel('Year', fontsize=13, fontweight='bold')
    
    # Adjust layout
    plt.subplots_adjust(hspace=0.35, left=0.1, right=0.92, top=0.95, bottom=0.12)
    
    return fig

# ============================================================================
# ADD GLOBAL SPEI COLORBAR AT BOTTOM
# ============================================================================

def add_global_colorbar(fig, results):
    """Add a single global colorbar for SPEI at the bottom of the figure."""
    
    # Collect all SPEI values
    all_spei = []
    for r in results:
        all_spei.extend(r['spei_values'])
    
    # Create colorbar axes (closer to x-axis)
    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.018])  # Moved closer to bottom
    
    # Create normalization and colormap
    norm = Normalize(vmin=np.min(all_spei), vmax=np.max(all_spei))
    sm = ScalarMappable(cmap=SPEI_CMAP, norm=norm)
    sm.set_array([])
    
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('SPEI-48', fontsize=20, fontweight='bold', labelpad=5)
    cbar.ax.tick_params(labelsize=15)
    
    # Set reasonable tick locations
    tick_vals = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]
    cbar.set_ticks([t for t in tick_vals if norm.vmin <= t <= norm.vmax])
    
    return fig

# ============================================================================
# MAIN WORKFLOW
# ============================================================================

def main():
    # Load data
    df_wue, df_spei = load_and_process_data()
    
    # Analyze each region (breakpoint from SPEI)
    results = []
    for region in REGION_ORDER:
        if region in REGION_COLORS:
            result = analyze_region(df_wue, df_spei, region, REGION_COLORS[region])
            results.append(result)
    
    # Create merged figure
    print("\n🎨 Creating merged figure...")
    fig = create_merged_figure(results)
    fig = add_global_colorbar(fig, results)
    
    # Save figure
    print("\n💾 Saving outputs...")
    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)
    
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches='tight', facecolor='white')
    
    print(f"✅ PNG saved: {OUTPUT_PNG}")
    print(f"✅ PDF saved: {OUTPUT_PDF}")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY TABLE (Breakpoint from SPEI-48)")
    print("="*80)
    print(f"{'Region':<20} {'Breakpoint':<12} {'Δp(WUE_T)':<12} {'ΔSPEI':<10} {'PRE SPEI Range':<18} {'POST SPEI Range'}")
    print("-" * 85)
    for r in results:
        pre_range = f"[{r['spei_min_pre']:.2f}, {r['spei_max_pre']:.2f}]" if not np.isnan(r['spei_min_pre']) else "N/A"
        post_range = f"[{r['spei_min_post']:.2f}, {r['spei_max_post']:.2f}]" if not np.isnan(r['spei_min_post']) else "N/A"
        print(f"{r['region']:<20} {r['breakpoint_str']:<12} {r['wue_delta']:+.2f}        {r['spei_delta']:+.2f}       {pre_range:<18} {post_range}")
    
    print("\n" + "="*80)
    print("✅ MERGED FIGURE COMPLETE!")
    print("="*80)
    print("\n📊 FIGURE FEATURES:")
    print("  • Breakpoint calculated from SPEI-48 (not WUE)")
    print("  • Colormap: FireBrick → Light Pink → Light Blue → Deep Navy")
    print("  • PRE gradient line: LEFT of breakpoint (with min/max labels)")
    print("  • POST gradient line: RIGHT of breakpoint (with min/max labels)")
    print("  • Left axis: p(WUE_T Decline) frequency")
    print("  • Text box: Δp(WUE_T) not ΔWUE_T")
    print("  • Single global colorbar at bottom")
    
    plt.show()
    
    return results, fig

if __name__ == "__main__":
    results, fig = main()