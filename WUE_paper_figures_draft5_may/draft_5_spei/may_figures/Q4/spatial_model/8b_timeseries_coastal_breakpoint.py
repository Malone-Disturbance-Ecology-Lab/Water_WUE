"""
SUPPLEMENTARY FIGURE: BREAKPOINT ANALYSIS BY COASTAL REGION
Uses REAL regional data extracted from NetCDF files.
Runs Pettitt breakpoint test independently for each region.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')
# Set global font properties for better readability
import matplotlib
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Helvetica']
matplotlib.rcParams['font.weight'] = 'bold'

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input: Regional data from Step 1
INPUT_FILE = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\regional_analysis\all_regions_monthly_frac_less.csv"

# Output files
OUTPUT_PNG = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\SuppFigure_Coastal_Regions_Breakpoint.png"
OUTPUT_PDF = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\spatial_SPEI\logistic_model\pre_figure_analysis\figures\SuppFigure_Coastal_Regions_Breakpoint.pdf"

# Colors matching Figure 3
REGION_COLORS = {
    "Pacific": "#E67E22",
    "Gulf": "#F1C40F",
    "Southeast Atlantic": "#9B59B6",
    "Atlantic North": "#3498DB",
    "Alaska": "#1ABC9C"
}

# Region order - Alaska at the END
REGION_ORDER = ["Pacific", "Gulf", "Southeast Atlantic", "Atlantic North", "Alaska"]

# Total months (from your data) - will be shown once in figure
TOTAL_MONTHS = 182

# ============================================================================
# PETTITT BREAKPOINT TEST
# ============================================================================

def pettitt_test(series):
    """Pettitt test for change point detection."""
    n = len(series)
    if n < 10:
        return None, None
    
    valid_idx = np.isfinite(series)
    if np.sum(valid_idx) < 10:
        return None, None
    
    series_clean = series[valid_idx]
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
# FORMAT P-VALUE FUNCTION
# ============================================================================

def format_p_value(p_val):
    """Format p-value: <0.001 or 2 decimal places."""
    if p_val is None or np.isnan(p_val):
        return "p = N/A"
    elif p_val < 0.001:
        return "p < 0.001"
    else:
        return f"p = {p_val:.3f}"

# ============================================================================
# BREAKPOINT ANALYSIS FOR A REGION
# ============================================================================

def analyze_region_breakpoint(df_region, region_name, color):
    """Perform breakpoint analysis on real regional data."""
    print(f"\n📈 Analyzing {region_name}...")
    
    df_sorted = df_region.sort_values('date').reset_index(drop=True)
    
    # Convert to pandas datetime
    dates = pd.to_datetime(df_sorted['date'].values)
    values = df_sorted['frac_less'].values
    
    breakpoint_idx, p_value = pettitt_test(values)
    
    if breakpoint_idx is not None:
        breakpoint_date = dates[breakpoint_idx]
        breakpoint_year = breakpoint_date.year
        breakpoint_str = breakpoint_date.strftime('%b %Y')
        
        pre_values = values[:breakpoint_idx]
        post_values = values[breakpoint_idx+1:]
        
        pre_mean = np.mean(pre_values)
        post_mean = np.mean(post_values)
        delta = post_mean - pre_mean
        pct_change = (delta / pre_mean * 100) if pre_mean != 0 else np.nan
    else:
        breakpoint_idx = None
        breakpoint_str = "None"
        breakpoint_year = None
        pre_mean = np.nan
        post_mean = np.nan
        delta = np.nan
        pct_change = np.nan
    
    print(f"  Breakpoint: {breakpoint_str}")
    print(f"  PRE: {pre_mean:.2f} → POST: {post_mean:.2f}")
    print(f"  Δ: {delta:+.2f}")
    
    return {
        'region': region_name,
        'color': color,
        'dates': dates,
        'values': values,
        'breakpoint_idx': breakpoint_idx,
        'breakpoint_str': breakpoint_str,
        'breakpoint_year': breakpoint_year,
        'p_value': p_value,
        'pre_mean': pre_mean,
        'post_mean': post_mean,
        'delta': delta,
        'pct_change': pct_change,
        'n_months': len(values)
    }

# ============================================================================
# CONVERT TO YEARLY DATA
# ============================================================================

def aggregate_to_yearly(dates, values):
    """Convert monthly data to yearly means."""
    df_temp = pd.DataFrame({'date': dates, 'value': values})
    df_temp['year'] = df_temp['date'].dt.year
    yearly = df_temp.groupby('year')['value'].mean().reset_index()
    return yearly['year'].values, yearly['value'].values

# ============================================================================
# FIGURE CREATION
# ============================================================================

def create_supplementary_figure(all_results):
    """Create multi-panel figure with yearly aggregated data."""
    n_regions = len(all_results)
    
    # Reduced height per panel with more space between panels
    fig, axes = plt.subplots(n_regions, 1, figsize=(14, 3.5 * n_regions), sharex=False)
    
    if n_regions == 1:
        axes = [axes]
    
    panel_labels = ['a', 'b', 'c', 'd', 'e']
    
    for i, result in enumerate(all_results):
        ax = axes[i]
        
        region = result['region']
        color = result['color']
        dates = result['dates']
        values = result['values']
        breakpoint_idx = result['breakpoint_idx']
        breakpoint_year = result['breakpoint_year']
        breakpoint_str = result['breakpoint_str']
        
        # Convert to yearly averages (remove monthly noise)
        years, yearly_values = aggregate_to_yearly(dates, values)
        
        # Plot yearly data only
        ax.plot(years, yearly_values, color=color, linewidth=2.5, 
               marker='o', markersize=8, alpha=0.85)
        
        # Add breakpoint vertical line
        if breakpoint_idx is not None and breakpoint_year is not None:
            ax.axvline(x=breakpoint_year, color='red', linestyle='--', 
                      linewidth=2.5, alpha=0.9)
        
        # PRE mean line - from START year to breakpoint year (touching the breakpoint)
        if not np.isnan(result['pre_mean']) and result['pre_mean'] > 0 and breakpoint_year is not None:
            # Get start year (first year in data)
            start_year = min(years)
            # Draw from start_year to breakpoint_year (touches the line)
            ax.hlines(y=result['pre_mean'], xmin=start_year, xmax=breakpoint_year,
                     color='black', linestyle=':', linewidth=2.5, alpha=0.8)
        
        # POST mean line - from breakpoint year to END year (touching the breakpoint)
        if not np.isnan(result['post_mean']) and result['post_mean'] > 0 and breakpoint_year is not None:
            # Get end year (last year in data)
            end_year = max(years)
            # Draw from breakpoint_year to end_year (touches the line)
            ax.hlines(y=result['post_mean'], xmin=breakpoint_year, xmax=end_year,
                     color=color, linestyle=':', linewidth=2.5, alpha=0.8)
        
        # Formatting - INCREASED FONT SIZES
        ax.set_ylabel('Frequency of\nWUE$_T$ Decline', fontsize=15, fontweight='bold')
        ax.tick_params(axis='both', labelsize=13, width=1.5, length=6)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        
        # Region title
        ax.set_title(f'{region}', fontsize=18, fontweight='bold', color=color, pad=15)
        
        # Statistics text box at BOTTOM-LEFT
        p_text = format_p_value(result['p_value'])
        
        stat_text = (
            f"Breakpoint: {result['breakpoint_str']}\n"
            f"{p_text}\n"
            f"PRE: {result['pre_mean']:.2f},  POST: {result['post_mean']:.2f}\n"
            f"Δ: {result['delta']:+.2f} ({result['pct_change']:+.1f}%)"
        )
        
        ax.text(0.02, 0.02, stat_text, transform=ax.transAxes,
               fontsize=12, fontweight='bold', verticalalignment='bottom', 
               fontfamily='sans-serif',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, 
                        edgecolor=color, linewidth=1.5))
        
        # Panel label at TOP-LEFT
        ax.text(0.02, 0.96, f'({panel_labels[i]})', transform=ax.transAxes,
               fontsize=16, fontweight='bold', verticalalignment='top', 
               horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, 
                        edgecolor='black', linewidth=1.5))
        
        # Legend at BOTTOM-RIGHT showing PRE (black) and POST (coast color)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='black', linestyle=':', linewidth=2.5, label='PRE mean'),
            Line2D([0], [0], color=color, linestyle=':', linewidth=2.5, label='POST mean')
        ]
        
        ax.legend(handles=legend_elements, loc='lower right', fontsize=11, 
                 framealpha=0.9, edgecolor='gray', handlelength=1.5,
                 prop={'weight': 'bold', 'family': 'sans-serif'})
    
    # Format x-axis for bottom panel only
    axes[-1].set_xlabel('Year', fontsize=15, fontweight='bold')
    
    # Set x-axis limits for all panels
    for ax in axes:
        ax.set_xlim(1999.5, 2026)
        ax.set_ylim(0, 1)
    
    # Add n=182 note at bottom
    fig.text(0.5, 0.01, f'n = {TOTAL_MONTHS} months per region (Mar 2000 - Sep 2025)', 
            ha='center', fontsize=11, fontstyle='italic')
    
    # Adjust layout
    plt.subplots_adjust(hspace=0.45, left=0.08, right=0.92, top=0.96, bottom=0.06)
    
    return fig
    """Create multi-panel figure with yearly aggregated data."""
    n_regions = len(all_results)
    
    # Reduced height per panel with more space between panels
    fig, axes = plt.subplots(n_regions, 1, figsize=(14, 3.5 * n_regions), sharex=False)
    
    if n_regions == 1:
        axes = [axes]
    
    panel_labels = ['a', 'b', 'c', 'd', 'e']
    
    for i, result in enumerate(all_results):
        ax = axes[i]
        
        region = result['region']
        color = result['color']
        dates = result['dates']
        values = result['values']
        breakpoint_idx = result['breakpoint_idx']
        breakpoint_year = result['breakpoint_year']
        breakpoint_str = result['breakpoint_str']
        
        # Convert to yearly averages (remove monthly noise)
        years, yearly_values = aggregate_to_yearly(dates, values)
        
        # Plot yearly data only
        ax.plot(years, yearly_values, color=color, linewidth=2.5, 
               marker='o', markersize=8, alpha=0.85)
        
        # Add breakpoint vertical line
        if breakpoint_idx is not None and breakpoint_year is not None:
            ax.axvline(x=breakpoint_year, color='red', linestyle='--', 
                      linewidth=2.5, alpha=0.9)
        
        # PRE mean line - ONLY for years BEFORE breakpoint
        if not np.isnan(result['pre_mean']) and result['pre_mean'] > 0 and breakpoint_year is not None:
            # Get pre-breakpoint years (years <= breakpoint_year - 1)
            pre_years = [y for y in years if y < breakpoint_year]
            if pre_years:
                # Draw horizontal line only over pre-breakpoint period
                ax.hlines(y=result['pre_mean'], xmin=min(pre_years), xmax=max(pre_years),
                         color='black', linestyle=':', linewidth=2.5, alpha=0.8)
        
        # POST mean line - ONLY for years AFTER breakpoint
        if not np.isnan(result['post_mean']) and result['post_mean'] > 0 and breakpoint_year is not None:
            # Get post-breakpoint years (years > breakpoint_year)
            post_years = [y for y in years if y > breakpoint_year]
            if post_years:
                # Draw horizontal line only over post-breakpoint period
                ax.hlines(y=result['post_mean'], xmin=min(post_years), xmax=max(post_years),
                         color=color, linestyle=':', linewidth=2.5, alpha=0.8)
        
        # Formatting - INCREASED FONT SIZES
        ax.set_ylabel('Frequency of\nWUE$_T$ Decline', fontsize=15, fontweight='bold')
        ax.tick_params(axis='both', labelsize=13, width=1.5, length=6)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        
        # Region title
        ax.set_title(f'{region}', fontsize=18, fontweight='bold', color=color, pad=15)
        
        # Statistics text box at BOTTOM-LEFT
        p_text = format_p_value(result['p_value'])
        
        stat_text = (
            f"Breakpoint: {result['breakpoint_str']}\n"
            f"{p_text}\n"
            f"PRE: {result['pre_mean']:.2f}  POST: {result['post_mean']:.2f}\n"
            f"Δ: {result['delta']:+.2f} ({result['pct_change']:+.1f}%)"
        )
        
        ax.text(0.02, 0.02, stat_text, transform=ax.transAxes,
               fontsize=12, fontweight='bold', verticalalignment='bottom', 
               fontfamily='sans-serif',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, 
                        edgecolor=color, linewidth=1.5))
        
        # Panel label at TOP-LEFT
        ax.text(0.02, 0.96, f'({panel_labels[i]})', transform=ax.transAxes,
               fontsize=16, fontweight='bold', verticalalignment='top', 
               horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, 
                        edgecolor='black', linewidth=1.5))
        
        # Legend at BOTTOM-RIGHT showing PRE (black) and POST (coast color)
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='black', linestyle=':', linewidth=2.5, label='PRE mean'),
            Line2D([0], [0], color=color, linestyle=':', linewidth=2.5, label='POST mean')
        ]
        
        ax.legend(handles=legend_elements, loc='lower right', fontsize=11, 
                 framealpha=0.9, edgecolor='gray', handlelength=1.5,
                 prop={'weight': 'bold', 'family': 'sans-serif'})
    
    # Format x-axis for bottom panel only
    axes[-1].set_xlabel('Year', fontsize=15, fontweight='bold')
    
    # Set x-axis limits for all panels
    for ax in axes:
        ax.set_xlim(1999.5, 2026)
        ax.set_ylim(0, 1)
    
    # Add n=182 note at bottom
    fig.text(0.5, 0.01, f'n = {TOTAL_MONTHS} months per region (Mar 2000 - Sep 2025)', 
            ha='center', fontsize=11, fontstyle='italic')
    
    # Adjust layout
    plt.subplots_adjust(hspace=0.45, left=0.08, right=0.92, top=0.96, bottom=0.06)
    
    return fig
    """Create multi-panel figure with yearly aggregated data."""
    n_regions = len(all_results)
    
    # Reduced height per panel with more space between panels
    fig, axes = plt.subplots(n_regions, 1, figsize=(14, 3.5 * n_regions), sharex=False)
    
    if n_regions == 1:
        axes = [axes]
    
    panel_labels = ['a', 'b', 'c', 'd', 'e']
    
    for i, result in enumerate(all_results):
        ax = axes[i]
        
        region = result['region']
        color = result['color']
        dates = result['dates']
        values = result['values']
        breakpoint_idx = result['breakpoint_idx']
        breakpoint_year = result['breakpoint_year']
        
        # Convert to yearly averages (remove monthly noise)
        years, yearly_values = aggregate_to_yearly(dates, values)
        
        # Plot yearly data only (no monthly dots/lines, no MA)
        ax.plot(years, yearly_values, color=color, linewidth=2.5, 
               marker='o', markersize=8, alpha=0.85)
        
        # Add breakpoint vertical line
        if breakpoint_idx is not None and breakpoint_year is not None:
            ax.axvline(x=breakpoint_year, color='red', linestyle='--', 
                      linewidth=2.5, alpha=0.9)
        
        # PRE mean line - BLACK (always black for all panels)
        if not np.isnan(result['pre_mean']) and result['pre_mean'] > 0:
            ax.axhline(y=result['pre_mean'], color='black', linestyle=':', 
                      linewidth=2.5, alpha=0.8)
        
        # POST mean line - SAME AS COAST COLOR (dotted)
        if not np.isnan(result['post_mean']) and result['post_mean'] > 0:
            ax.axhline(y=result['post_mean'], color=color, linestyle=':', 
                      linewidth=2.5, alpha=0.8)
        
        # Formatting - INCREASED FONT SIZES
        ax.set_ylabel('Frequency of\nWUE$_T$ Decline', fontsize=15, fontweight='bold')
        ax.tick_params(axis='both', labelsize=13, width=1.5, length=6)
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.8)
        
        # Region title
        ax.set_title(f'{region}', fontsize=18, fontweight='bold', color=color, pad=15)
        
        # Statistics text box at BOTTOM-LEFT - bold and Arial
        p_text = format_p_value(result['p_value'])
        
        stat_text = (
            f"Breakpoint: {result['breakpoint_str']}\n"
            f"{p_text}\n"
            f"PRE: {result['pre_mean']:.2f}  POST: {result['post_mean']:.2f}\n"
            f"Δ: {result['delta']:+.2f} ({result['pct_change']:+.1f}%)"
        )
        
        ax.text(0.02, 0.02, stat_text, transform=ax.transAxes,
               fontsize=12, fontweight='bold', verticalalignment='bottom', 
               fontfamily='sans-serif',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, 
                        edgecolor=color, linewidth=1.5))
        
        # Panel label at TOP-LEFT
        ax.text(0.02, 0.96, f'({panel_labels[i]})', transform=ax.transAxes,
               fontsize=16, fontweight='bold', verticalalignment='top', 
               horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, 
                        edgecolor='black', linewidth=1.5))
        
        # Add a small legend in the BOTTOM-RIGHT of each panel showing POST line color
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color=color, linestyle=':', linewidth=2.5, label=f'POST mean')
        ]
        
        ax.legend(handles=legend_elements, loc='lower right', fontsize=10, 
                 framealpha=0.9, edgecolor=color, handlelength=1.5,
                 prop={'weight': 'bold', 'family': 'sans-serif'})
    
    # Format x-axis for bottom panel only
    axes[-1].set_xlabel('Year', fontsize=15, fontweight='bold')
    
    # Set x-axis limits for all panels
    for ax in axes:
        ax.set_xlim(1999.5, 2026)
        ax.set_ylim(0, 1)
    
    # Add n=182 note at bottom
    fig.text(0.5, 0.01, f'n = {TOTAL_MONTHS} months per region (Mar 2000 - Sep 2025)', 
            ha='center', fontsize=11, fontstyle='italic')
    
    # Adjust layout
    plt.subplots_adjust(hspace=0.45, left=0.08, right=0.92, top=0.96, bottom=0.06)
    
    return fig
# ============================================================================
# LOAD DATA AND RUN ANALYSIS
# ============================================================================

def main():
    print("="*80)
    print("SUPPLEMENTARY FIGURE: COASTAL REGION BREAKPOINT ANALYSIS")
    print("="*80)
    
    # Load regional data
    df = pd.read_csv(INPUT_FILE)
    df['date'] = pd.to_datetime(df['date'])
    
    print(f"\n📊 Loaded {len(df)} records from {INPUT_FILE}")
    print(f"   Regions: {df['region'].unique().tolist()}")
    print(f"   Total months per region: {TOTAL_MONTHS}")
    
    # Analyze each region in specified order (Alaska last)
    all_results = []
    for region in REGION_ORDER:
        if region not in REGION_COLORS:
            continue
        df_region = df[df['region'] == region].copy()
        if len(df_region) > 0:
            result = analyze_region_breakpoint(df_region, region, REGION_COLORS[region])
            all_results.append(result)
    
    # Create figure
    print("\n🎨 Creating supplementary figure...")
    fig = create_supplementary_figure(all_results)
    
    # Save figure (high resolution)
    print("\n💾 Saving outputs...")
    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)
    
    fig.savefig(OUTPUT_PNG, dpi=1200, bbox_inches='tight', facecolor='white')
    fig.savefig(OUTPUT_PDF, dpi=1200, bbox_inches='tight', facecolor='white')
    
    print(f"✅ PNG saved: {OUTPUT_PNG}")
    print(f"✅ PDF saved: {OUTPUT_PDF}")
    
    # Print summary for manuscript
    print("\n" + "="*80)
    print("BREAKPOINT RESULTS SUMMARY")
    print("="*80)
    print(f"{'Region':<20} {'Breakpoint':<15} {'PRE':<8} {'POST':<8} {'Δ':<10} {'p-value'}")
    print("-" * 70)
    for r in all_results:
        p_str = format_p_value(r['p_value'])
        print(f"{r['region']:<20} {r['breakpoint_str']:<15} {r['pre_mean']:.2f}   {r['post_mean']:.2f}   {r['delta']:+.2f}     {p_str}")
    
    plt.show()
    print("\n✅ Figure complete!")
    
    return all_results, fig

if __name__ == "__main__":
    results, fig = main()