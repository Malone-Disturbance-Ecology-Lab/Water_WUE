# -*- coding: utf-8 -*-
"""
create_two_panel_resistance_figure.py

Reads:
- pixel_level_spatial_resistance_analogue_by_coast.csv (Panel a)
- EDI_response_v2_upland_monthly_coast_summary.csv (Panel b)

Calculates monthly resistance analogue = 1 + dry_mean_WUE_T_pct_change / 100
Filters months where dry_n_pixels > 0.

Creates a two-panel figure:
- a) Pixel-level mean dry-month resistance analogue (violin + jittered points + boxplot overlay)
- b) Coast-month dry-exposure resistance analogue (violin + jittered points + boxplot overlay)

Coasts are ordered by ascending mean resistance from the pixel data.
Saves PNG, PDF, derived monthly CSV, and summary CSV in talib_publication.
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
INPUT_DIR = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")
OUTPUT_DIR = INPUT_DIR / "talib_publication"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files (still in spatial_performance_metric)
PIXEL_CSV = INPUT_DIR / "spatial_performance_metric" / "pixel_level_spatial_resistance_analogue_by_coast.csv"
MONTHLY_CSV = INPUT_DIR / "EDI_response_v2_upland_monthly_coast_summary.csv"

# Output files (now in talib_publication)
FIGURE_PNG = OUTPUT_DIR / "spatial_resistance_analogue_two_panel.png"
FIGURE_PDF = OUTPUT_DIR / "spatial_resistance_analogue_two_panel.pdf"
DERIVED_CSV = OUTPUT_DIR / "monthly_spatial_resistance_analogue_by_coast.csv"
SUMMARY_CSV = OUTPUT_DIR / "spatial_resistance_analogue_two_panel_summary.csv"

# Original coast names and colors (will be reordered later)
ORIGINAL_COAST_ORDER = ["AK Coast", "Pacific Coast", "Gulf Coast", "Atlantic Coast"]
CUSTOM_COLORS = {
    "Atlantic Coast": "#2E8B57",   # SeaGreen
    "Pacific Coast":  "#DC143C",   # Crimson
    "Gulf Coast":     "#00CED1",   # DarkTurquoise
    "AK Coast":       "#8B4513"    # SaddleBrown
}
# Mapping for tick labels (remove "Coast" and shorten AK to Alaska)
TICK_LABELS = {
    "AK Coast": "Alaska",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf",
    "Atlantic Coast": "Atlantic"
}

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------
def read_pixel_data(file_path):
    """Read pixel-level CSV and return DataFrame with required columns."""
    df = pd.read_csv(file_path)
    res_cols = [c for c in df.columns if 'resistance' in c.lower()]
    if not res_cols:
        raise ValueError("No resistance column found in pixel CSV.")
    res_col = res_cols[0]
    print(f"Pixel CSV: using resistance column '{res_col}'")
    df['coast_region'] = df['coast_region'].str.strip()
    df = df[df['coast_region'].isin(ORIGINAL_COAST_ORDER)].copy()
    df = df.dropna(subset=['coast_region', res_col])
    df.rename(columns={res_col: 'resistance'}, inplace=True)
    # Keep as categorical but we'll reorder later based on mean
    return df

def read_monthly_data(file_path):
    """Read monthly summary CSV and compute monthly resistance."""
    df = pd.read_csv(file_path)
    required = ['year', 'month', 'coast_region', 'dry_n_pixels', 'dry_mean_WUE_T_pct_change']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in monthly CSV.")
    df = df[df['dry_n_pixels'] > 0].copy()
    df['monthly_resistance'] = 1 + df['dry_mean_WUE_T_pct_change'] / 100.0
    df['coast_region'] = df['coast_region'].str.strip()
    df = df[df['coast_region'].isin(ORIGINAL_COAST_ORDER)].copy()
    df = df.dropna(subset=['coast_region', 'monthly_resistance'])
    return df

def sample_points(data, n_max=500):
    """Subsample up to n_max points per coast for jitter plotting."""
    sampled = data.groupby('coast_region', group_keys=False).apply(
        lambda x: x.sample(min(len(x), n_max), random_state=42) if len(x) > 0 else x
    )
    return sampled

def add_violin_with_jitter(ax, data, y_col, coast_order, colors, ylabel, panel_label,
                           sample_n=500):
    """
    Draw violin plot, jittered points, and an overlaid boxplot on axis ax.
    Violins are made wider (width=0.85) for better visibility.
    """
    # Prepare data per coast in the given order
    groups = [data[data['coast_region'] == coast][y_col].dropna().values
              for coast in coast_order]
    positions = np.arange(1, len(coast_order) + 1)

    # Violin plot with wider violins
    vp = ax.violinplot(groups, positions=positions, showmeans=False,
                       showmedians=True, showextrema=False, widths=0.85)

    # Set median line properties
    vp['cmedians'].set_color('black')
    vp['cmedians'].set_linewidth(2)

    # Color violins
    for i, coast in enumerate(coast_order):
        vp['bodies'][i].set_facecolor(colors[coast])
        vp['bodies'][i].set_alpha(0.6)
        vp['bodies'][i].set_edgecolor('none')

    # Overlay a narrow boxplot to show quartiles and range
    box_width = 0.18
    for i, (coast, vals) in enumerate(zip(coast_order, groups)):
        if len(vals) == 0:
            continue
        bp = ax.boxplot(vals, positions=[positions[i]], widths=box_width,
                        showfliers=False, patch_artist=False,
                        boxprops=dict(linewidth=1.2, color=colors[coast]),
                        whiskerprops=dict(linewidth=1.2, color=colors[coast]),
                        capprops=dict(linewidth=1.2, color=colors[coast]),
                        medianprops=dict(linewidth=2, color='black'))

    # Jittered points (subsample)
    sampled = sample_points(data, n_max=sample_n)
    for i, coast in enumerate(coast_order):
        vals = sampled[sampled['coast_region'] == coast][y_col].values
        if len(vals) > 0:
            x_jitter = np.random.normal(i + 1, 0.04, size=len(vals))
            ax.scatter(x_jitter, vals, s=8, alpha=0.25,
                       color=colors[coast], rasterized=True)

    # Horizontal reference line: only y=1.0 (solid black)
    ax.axhline(y=1.0, color='black', linestyle='-', linewidth=1.2, alpha=0.7)

    # Axis labels
    ax.set_ylabel(ylabel, fontsize=18)          # increased y-label size
    ax.set_xlabel("Region", fontsize=18)        # add x-axis label with increased size

    # Set tick labels with custom mapping (remove "Coast")
    tick_labels = [TICK_LABELS[coast] for coast in coast_order]
    ax.set_xticks(positions)
    ax.set_xticklabels(tick_labels, fontsize=18)  # increased tick font
    ax.tick_params(axis='y', labelsize=18)        # increased y-tick font

    # Format y-axis ticks to 2 decimal places
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

    # Add panel label a) or b) in top-left corner, larger font
    ax.text(0.02, 0.95, panel_label, transform=ax.transAxes,
            fontsize=18, fontweight='bold', va='top', ha='left')

    # Add n= annotations, moved very low and font size increased
    n_vals = [len(g) for g in groups]
    ymin, ymax = ax.get_ylim()
    ypos = ymin + 0.001 * (ymax - ymin)   # very close to bottom axis
    for i, n in enumerate(n_vals):
        ax.text(positions[i], ypos, f'n={n}', ha='center', va='bottom',
                fontsize=16, color='black')   # increased to 16

    # Clean spines and add light grid
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.3, linewidth=0.5)
    ax.set_axisbelow(True)

    return vp

# ----------------------------------------------------------------------
# Main script
# ----------------------------------------------------------------------
def main():
    print("\n" + "=" * 80)
    print("TWO-PANEL SPATIAL RESISTANCE ANALOGUE FIGURE (VIOLIN + JITTER + BOXPLOT)")
    print("=" * 80)

    # 1. Read data
    print("\nReading pixel-level CSV...")
    df_pixel = read_pixel_data(PIXEL_CSV)
    print(f"  Total valid pixels: {len(df_pixel)}")

    print("\nReading monthly summary CSV...")
    df_monthly = read_monthly_data(MONTHLY_CSV)
    print(f"  Total valid dry coast-months: {len(df_monthly)}")

    # 2. Determine coast order by ascending mean resistance from pixel data
    pixel_means = df_pixel.groupby('coast_region')['resistance'].mean().sort_values()
    coast_order = pixel_means.index.tolist()  # ascending order
    print("\nCoasts ordered by mean resistance (ascending):")
    for coast in coast_order:
        print(f"  {coast}: mean = {pixel_means[coast]:.4f}")

    # 3. Save derived monthly CSV
    df_monthly_out = df_monthly[['year', 'month', 'coast_region',
                                 'dry_n_pixels', 'dry_mean_WUE_T_pct_change',
                                 'monthly_resistance']].copy()
    df_monthly_out.to_csv(DERIVED_CSV, index=False)
    print(f"\nDerived monthly CSV saved: {DERIVED_CSV}")

    # 4. Create summary statistics for both panels (using the new order)
    summary_rows = []
    for panel, data, ycol in [('Pixel', df_pixel, 'resistance'),
                              ('Monthly', df_monthly, 'monthly_resistance')]:
        for coast in coast_order:  # use ordered list
            vals = data[data['coast_region'] == coast][ycol].dropna().values
            if len(vals) == 0:
                continue
            summary_rows.append({
                'panel': panel,
                'coast_region': coast,
                'n': len(vals),
                'mean': np.mean(vals),
                'median': np.median(vals),
                'sd': np.std(vals, ddof=1),
                'q25': np.percentile(vals, 25),
                'q75': np.percentile(vals, 75),
                'iqr': np.percentile(vals, 75) - np.percentile(vals, 25),
                'min': np.min(vals),
                'max': np.max(vals)
            })
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_CSV, index=False)
    print(f"Summary CSV saved: {SUMMARY_CSV}")

    # Print per-coast stats to console
    print("\n" + "-" * 80)
    print("PER-COAST STATISTICS (ordered by pixel mean resistance)")
    print("-" * 80)
    for panel in ['Pixel', 'Monthly']:
        print(f"\n{panel} panel:")
        sub = df_summary[df_summary['panel'] == panel]
        for coast in coast_order:
            row = sub[sub['coast_region'] == coast]
            if row.empty:
                continue
            r = row.iloc[0]
            print(f"  {coast}: n={r['n']}, mean={r['mean']:.4f}, median={r['median']:.4f}, "
                  f"IQR={r['iqr']:.4f}")

    # 5. Create figure – larger size and bigger violins
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), dpi=100)
    fig.patch.set_facecolor('white')

    # Panel a: Pixel-level
    add_violin_with_jitter(ax1, df_pixel, 'resistance', coast_order,
                           CUSTOM_COLORS,
                           ylabel='Resistance analogue from grid-cell means',
                           panel_label='a)',
                           sample_n=500)

    # Panel b: Monthly
    add_violin_with_jitter(ax2, df_monthly, 'monthly_resistance', coast_order,
                           CUSTOM_COLORS,
                           ylabel='Resistance analogue from monthly means',
                           panel_label='b)',
                           sample_n=500)

    # Adjust y-limits to show reference lines clearly
    for ax in [ax1, ax2]:
        ymin, ymax = ax.get_ylim()
        y_range = ymax - ymin
        ax.set_ylim(ymin - 0.05 * y_range, ymax + 0.05 * y_range)

    plt.tight_layout()

    # Save
    fig.savefig(FIGURE_PNG, dpi=600, bbox_inches='tight')
    fig.savefig(FIGURE_PDF, bbox_inches='tight')
    print(f"\nFigure PNG saved: {FIGURE_PNG}")
    print(f"Figure PDF saved: {FIGURE_PDF}")

    # 6. Display
    plt.show(block=True)
    print("Figure displayed. If not visible, check Spyder Plots pane / graphics backend.")

    print("\n" + "=" * 80)
    print("SCRIPT COMPLETED SUCCESSFULLY")
    print("=" * 80)

if __name__ == "__main__":
    main()