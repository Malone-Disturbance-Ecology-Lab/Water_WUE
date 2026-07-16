# -*- coding: utf-8 -*-
"""
Panels g–h: Violin (g) and summary (h) with background points, p‑values placed right.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
fig_output = os.path.join(figure_dir, "Q1_panels_g_h_violin_summary.png")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
wue_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_WUE_T_monthly_by_coast_waterclass.csv"))
wue_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_WUE_T_site_level_by_coast_waterclass.csv"))

tet_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_TET_monthly_by_coast_waterclass.csv"))
tet_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_TET_site_level_by_coast_waterclass.csv"))

# ------------------------------------------------------------------------------
# Cap WUE_T at 99th percentile PER ECOSYSTEM CLASS (for violin display only)
# ------------------------------------------------------------------------------
def cap_by_group(df, group_col, value_col, percentile=99):
    capped = df.copy()
    capped['WUE_T_capped'] = np.nan
    for group in df[group_col].unique():
        mask = df[group_col] == group
        vals = df.loc[mask, value_col].dropna()
        if len(vals) > 0:
            p = np.percentile(vals, percentile)
            capped.loc[mask, 'WUE_T_capped'] = np.minimum(df.loc[mask, value_col], p)
    return capped

wue_monthly = cap_by_group(wue_monthly, 'water_class', 'WUE_T', 99)
print("\nGroup-specific 99th percentiles for WUE_T (violin cap):")
for group in wue_monthly['water_class'].unique():
    mask = wue_monthly['water_class'] == group
    vals = wue_monthly.loc[mask, 'WUE_T'].dropna()
    if len(vals) > 0:
        print(f"  {group}: {np.percentile(vals, 99):.3f}")

# ------------------------------------------------------------------------------
# Colours and short names
# ------------------------------------------------------------------------------
ecosystem_colors = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}
coast_colors = {
    'Atlantic Coast': '#2E8B57',
    'Gulf Coast':     '#00CED1',
    'Pacific Coast':  '#DC143C',
    'AK Coast':       '#8B4513'
}
coast_short = {
    'Atlantic Coast': 'Atlantic',
    'Pacific Coast': 'Pacific',
    'Gulf Coast': 'Gulf',
    'AK Coast': 'AK'
}

# ------------------------------------------------------------------------------
# Order groups
# ------------------------------------------------------------------------------
eco_order = ["Upland", "Freshwater", "Saline"]  # fixed for ecosystem

# Coast order by mean T:ET ascending (AK, Pacific, Gulf, Atlantic)
coast_means_order = tet_site.groupby('coast_region')['mean_TET'].mean().sort_values()
coast_order = coast_means_order.index.tolist()
print("\nCoast order (low to high mean T:ET):", coast_order)

# Convert to categorical for consistent ordering
wue_monthly['water_class'] = pd.Categorical(wue_monthly['water_class'], categories=eco_order, ordered=True)
wue_site['water_class'] = pd.Categorical(wue_site['water_class'], categories=eco_order, ordered=True)
tet_monthly['coast_region'] = pd.Categorical(tet_monthly['coast_region'], categories=coast_order, ordered=True)
tet_site['coast_region'] = pd.Categorical(tet_site['coast_region'], categories=coast_order, ordered=True)

# ------------------------------------------------------------------------------
# Compute Kruskal‑Wallis p‑values (site‑level means – uncapped)
# ------------------------------------------------------------------------------
eco_groups = [wue_site[wue_site['water_class'] == e]['mean_WUE_T'].dropna().values for e in eco_order]
h_eco, p_eco = kruskal(*eco_groups)
p_text_eco = f"p = {p_eco:.4f}" if p_eco >= 0.001 else "p < 0.001"

coast_groups = [tet_site[tet_site['coast_region'] == c]['mean_TET'].dropna().values for c in coast_order]
h_coast, p_coast = kruskal(*coast_groups)
p_text_coast = f"p = {p_coast:.3f}" if p_coast >= 0.001 else "p < 0.001"

# ------------------------------------------------------------------------------
# Group statistics for panel h (mean, CI)
# ------------------------------------------------------------------------------
coast_stats = tet_site.groupby('coast_region')['mean_TET'].agg(['mean', 'std', 'count']).reindex(coast_order)
coast_means = coast_stats['mean']
coast_sem = coast_stats['std'] / np.sqrt(coast_stats['count'])
coast_ci = 1.96 * coast_sem

# For panel g median and mean labels (from uncapped data)
eco_medians = wue_monthly.groupby('water_class')['WUE_T'].median().reindex(eco_order)
eco_means_site = wue_site.groupby('water_class')['mean_WUE_T'].mean().reindex(eco_order)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
fig, (ax_g, ax_h) = plt.subplots(1, 2, figsize=(24, 6))

plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

# Style both axes
for ax in [ax_g, ax_h]:
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ------------------------------------------------------------------------------
# Panel g: Ecosystem – violin plot (capped per group) with background points
# ------------------------------------------------------------------------------
positions = np.arange(len(eco_order))
data_violin = [wue_monthly[wue_monthly['water_class'] == e]['WUE_T_capped'].dropna().values for e in eco_order]

# Violin plot (capped for visual)
violin_parts = ax_g.violinplot(data_violin, positions=positions, showmeans=False, showmedians=False, widths=0.8)
for i, pc in enumerate(violin_parts['bodies']):
    pc.set_facecolor(ecosystem_colors[eco_order[i]])
    pc.set_edgecolor(ecosystem_colors[eco_order[i]])
    pc.set_alpha(0.7)

# Add jittered background points (not dim) – use capped values for consistency
for i, eco in enumerate(eco_order):
    vals = wue_monthly[wue_monthly['water_class'] == eco]['WUE_T_capped'].dropna().values
    x_jitter = np.random.normal(i, 0.06, size=len(vals))
    ax_g.scatter(x_jitter, vals, alpha=0.5, s=20, color=ecosystem_colors[eco],
                 edgecolors='none')  # no edge to keep clean

# Set x‑axis limits
ax_g.set_xlim(-0.5, len(eco_order) - 0.5)

# Adjust y‑axis top to give room for mean labels
max_data = max([np.max(vals) for vals in data_violin]) if data_violin else 1
max_label = eco_means_site.max() if not eco_means_site.empty else 1
ylim_top = max(max_data, max_label) * 1.15
ax_g.set_ylim(0, ylim_top)

# Median lines (black) – using uncapped monthly medians
for i, eco in enumerate(eco_order):
    med = eco_medians[eco]
    ax_g.plot([i - 0.15, i + 0.15], [med, med],
              color='black', linewidth=2.5, solid_capstyle='butt')

# Mean labels (uncapped site‑level means)
ymin_g, ymax_g = ax_g.get_ylim()
for i, eco in enumerate(eco_order):
    mean_val = eco_means_site[eco]
    ax_g.text(i, ymax_g - 0.02 * (ymax_g - ymin_g), f"{mean_val:.2f}",
              ha='center', va='top', fontsize=24, color='black', fontweight='bold')

# p‑value box – top right
ax_g.text(0.45, 0.80, p_text_eco, transform=ax_g.transAxes,
          fontsize=24, ha='right', va='top',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# x‑axis labels and colours
ax_g.set_xticks(positions)
ax_g.set_xticklabels(eco_order)
for tick, eco in zip(ax_g.get_xticklabels(), eco_order):
    tick.set_color(ecosystem_colors[eco])

ax_g.set_ylabel('WUE$_T$ (g C kg$^{-1}$ H$_2$O$^{-1}$)', fontsize=28, labelpad=10)
ax_g.set_xlabel('')
ax_g.set_title('')
ax_g.text(-0.15, 1.02, 'g)', transform=ax_g.transAxes,
          fontsize=34, fontweight='bold', va='bottom', ha='left')

# ------------------------------------------------------------------------------
# Panel h: Coast – mean ± 95% CI with background points (not dim)
# ------------------------------------------------------------------------------
x_pos_h = np.arange(len(coast_order))

# Background site‑level points (not dim)
for i, coast in enumerate(coast_order):
    vals = tet_site[tet_site['coast_region'] == coast]['mean_TET'].dropna().values
    x_jitter = np.random.normal(i, 0.06, size=len(vals))
    ax_h.scatter(x_jitter, vals, alpha=0.5, s=30, color=coast_colors[coast],
                 edgecolors='none')

# Mean ± CI
for i, coast in enumerate(coast_order):
    ax_h.errorbar(i, coast_means[coast], yerr=coast_ci[coast],
                  fmt='o', color=coast_colors[coast],
                  markersize=16, capsize=10, elinewidth=3,
                  markeredgecolor='black', markeredgewidth=1.2)

# Mean labels above
ymin_h, ymax_h = 0, 1.05  # fixed for T:ET
ax_h.set_ylim(ymin_h, ymax_h)
for i, coast in enumerate(coast_order):
    ax_h.text(i, ymax_h - 0.1*(ymax_h-ymin_h), f"{coast_means[coast]:.2f}",
              ha='center', va='bottom', fontsize=24, color='black', fontweight='bold')

ax_h.set_xticks(x_pos_h)
ax_h.set_xticklabels([coast_short[c] for c in coast_order])
for tick, coast in zip(ax_h.get_xticklabels(), coast_order):
    tick.set_color(coast_colors[coast])

ax_h.axhline(0.5, color='gray', linestyle='--', linewidth=2, alpha=0.6)

ax_h.set_ylabel('T:ET ratio', fontsize=28, labelpad=10)
ax_h.set_xlabel('')
ax_h.set_title('')
ax_h.text(-0.15, 1.02, 'h)', transform=ax_h.transAxes,
          fontsize=34, fontweight='bold', va='bottom', ha='left')

# p‑value box – bottom right (as originally)
ax_h.text(0.97, 0.05, p_text_coast, transform=ax_h.transAxes,
          fontsize=24, ha='right', va='bottom',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# ------------------------------------------------------------------------------
# Adjust spacing and save
# ------------------------------------------------------------------------------
plt.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.15, wspace=0.3)
fig.savefig(fig_output, dpi=600, bbox_inches='tight', pad_inches=0.5, facecolor='white')
print(f"Figure saved to: {fig_output}")

# ------------------------------------------------------------------------------
# Statistics (console only)
# ------------------------------------------------------------------------------
print("\n" + "=" * 80)
print("STATISTICS")
print("=" * 80)
print(f"Ecosystem (WUE_T) Kruskal-Wallis: H = {h_eco:.4f}, p = {p_eco:.8f}")
print(f"Coast (T:ET) Kruskal-Wallis: H = {h_coast:.4f}, p = {p_coast:.8f}")