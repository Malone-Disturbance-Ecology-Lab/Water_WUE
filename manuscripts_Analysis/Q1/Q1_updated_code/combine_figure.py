# -*- coding: utf-8 -*-
"""
Combined Q1 panels: A (a–c), B (d–f), and G/H (g–h)
Updated g/h: violin for ecosystem (capped at 99th percentile) and mean±CI for coast.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kruskal
import matplotlib.patches as patches

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
fig_output = os.path.join(figure_dir, "Q1_combined_panels_A_B_G_H.png")

# ------------------------------------------------------------------------------
# Load data for all panels
# ------------------------------------------------------------------------------
# Panel A
summary_file = os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv")
df_summary = pd.read_csv(summary_file)

# Panel B
smooth_file = os.path.join(output_dir, "Q1_final_smooth_gam_predictions_TET.csv")
gam_comp_file = os.path.join(output_dir, "Q1_final_gam_model_comparison.csv")
edf_file = os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv")
df_smooth = pd.read_csv(smooth_file)
gam_comp = pd.read_csv(gam_comp_file)
edf_df = pd.read_csv(edf_file)

# Panels g/h
wue_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_WUE_T_monthly_by_coast_waterclass.csv"))
wue_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_WUE_T_site_level_by_coast_waterclass.csv"))
tet_monthly = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                       "Q1_NN_TET_monthly_by_coast_waterclass.csv"))
tet_site = pd.read_csv(os.path.join(base_dir, "outputs", "T_ET_ratio_coast",
                                    "Q1_NN_TET_site_level_by_coast_waterclass.csv"))

# ------------------------------------------------------------------------------
# Common settings
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500",
}
coast_colors = {
    'Atlantic Coast': '#A50F15',   # dark red
    'Pacific Coast':  '#0072B2',   # blue
    'Gulf Coast':     '#4D4D4D',   # dark charcoal
    'AK Coast':       '#009E73'    # bluish green
}
coast_short = {
    'Atlantic Coast': 'Atlantic',
    'Pacific Coast': 'Pacific',
    'Gulf Coast': 'Gulf',
    'AK Coast': 'Alaska'
}

diff_label_mapping = {
    "WUE_ET - WUE_T": r"WUE$_{ET}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_ET - WUE_E": r"WUE$_{ET}$ – WUE$_{E}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_E - WUE_T": r"WUE$_{E}$ – WUE$_{T}$ (g C kg$^{-1}$ H$_2$O$^{-1}$)",
}
diff_labels = list(diff_label_mapping.keys())

# ------------------------------------------------------------------------------
# Data preparation for Panel A
# ------------------------------------------------------------------------------
df_summary = df_summary[df_summary["difference_label"].isin(diff_labels)]
df_summary["difference_label"] = pd.Categorical(df_summary["difference_label"],
                                                categories=diff_labels, ordered=True)
df_summary["water_class"] = pd.Categorical(df_summary["water_class"],
                                           categories=ecosystem_order, ordered=True)

# ------------------------------------------------------------------------------
# Data preparation for Panel B
# ------------------------------------------------------------------------------
df_smooth = df_smooth[df_smooth["difference_label"].isin(diff_labels)]
df_smooth["difference_label"] = pd.Categorical(df_smooth["difference_label"],
                                               categories=diff_labels, ordered=True)
df_smooth["ecosystem_label"] = pd.Categorical(df_smooth["ecosystem_label"],
                                              categories=ecosystem_order, ordered=True)

# Extract R² and p-value for Panel B
smooth_row = gam_comp[gam_comp["model"] == "q1_smooth_gam"].iloc[0]
smooth_r2 = smooth_row["adjusted_r_squared"]
edf_df_sub = edf_df[edf_df["term"].str.contains(r"s\(Trans_ratio\)", regex=True)].copy()
p_col = [c for c in edf_df_sub.columns if "p" in c.lower()][0]
max_tet_p = edf_df_sub[p_col].max()
smooth_p_text = "p < 0.001" if max_tet_p < 0.001 else f"p = {max_tet_p:.3f}"
model_text = f"Adjusted R² = {smooth_r2:.2f}; {smooth_p_text}"

# Parse edf values for Panel B
def parse_smooth_term(term):
    clean = term.replace("s(Trans_ratio):diff_ecosystem", "")
    diff_type, ecosystem = clean.split("__")
    diff_map = {
        "WUE_ET_minus_WUE_T": "WUE_ET - WUE_T",
        "WUE_ET_minus_WUE_E": "WUE_ET - WUE_E",
        "WUE_E_minus_WUE_T": "WUE_E - WUE_T",
    }
    return pd.Series({
        "difference_label": diff_map[diff_type],
        "ecosystem_label": ecosystem
    })
edf_parsed = edf_df_sub["term"].apply(parse_smooth_term)
edf_df_sub = pd.concat([edf_df_sub, edf_parsed], axis=1)
edf_lookup = {(row["difference_label"], row["ecosystem_label"]): row["edf"] for _, row in edf_df_sub.iterrows()}

# ------------------------------------------------------------------------------
# Data preparation for Panels g/h (capping for violin, ordering for coast)
# ------------------------------------------------------------------------------
# Cap WUE_T at 99th percentile per ecosystem for violin display
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

# Ecosystem order: fixed (Upland, Freshwater, Saline)
eco_order_gh = ["Upland", "Freshwater", "Saline"]

# Coast order: low to high mean T:ET (AK, Pacific, Gulf, Atlantic)
coast_means_order = tet_site.groupby('coast_region')['mean_TET'].mean().sort_values()
coast_order = coast_means_order.index.tolist()

# Convert to categorical for consistent ordering
wue_monthly['water_class'] = pd.Categorical(wue_monthly['water_class'], categories=eco_order_gh, ordered=True)
wue_site['water_class'] = pd.Categorical(wue_site['water_class'], categories=eco_order_gh, ordered=True)
tet_monthly['coast_region'] = pd.Categorical(tet_monthly['coast_region'], categories=coast_order, ordered=True)
tet_site['coast_region'] = pd.Categorical(tet_site['coast_region'], categories=coast_order, ordered=True)

# Statistics for g/h
eco_groups = [wue_site[wue_site['water_class'] == e]['mean_WUE_T'].dropna().values for e in eco_order_gh]
h_eco, p_eco = kruskal(*eco_groups)
p_text_eco = f"p = {p_eco:.4f}" if p_eco >= 0.001 else "p < 0.001"

coast_groups = [tet_site[tet_site['coast_region'] == c]['mean_TET'].dropna().values for c in coast_order]
h_coast, p_coast = kruskal(*coast_groups)
p_text_coast = f"p = {p_coast:.3f}" if p_coast >= 0.001 else "p < 0.001"

# Group statistics for plotting
# For panel g: medians and means (uncapped)
eco_medians = wue_monthly.groupby('water_class')['WUE_T'].median().reindex(eco_order_gh)
eco_means = wue_site.groupby('water_class')['mean_WUE_T'].mean().reindex(eco_order_gh)

# For panel h: means and CI
coast_stats = tet_site.groupby('coast_region')['mean_TET'].agg(['mean', 'std', 'count']).reindex(coast_order)
coast_means = coast_stats['mean']
coast_sem = coast_stats['std'] / np.sqrt(coast_stats['count'])
coast_ci = 1.96 * coast_sem

# ------------------------------------------------------------------------------
# Font settings (applied before figure creation)
# ------------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

# ------------------------------------------------------------------------------
# Create outer figure with 3 rows, 1 column
# Height ratios: Panel A = 6, Panel B = 7, Panel g/h = 6
# Total height = 6+7+6 + gap (0.5*2) = 20 inches
# ------------------------------------------------------------------------------
fig_height = 6 + 7 + 6 + 0.5 * 2   # 20 inches
fig = plt.figure(figsize=(24, fig_height))

outer = fig.add_gridspec(
    3, 1,
    height_ratios=[6, 7, 6],
    hspace=0.30,
    left=0.12,
    right=0.95,
    top=0.96,
    bottom=0.06
)

# Row 1: 3 columns (panels a–c)
gs_top = outer[0].subgridspec(1, 3, wspace=0.30)
# Row 2: 3 columns (panels d–f)
gs_mid = outer[1].subgridspec(1, 3, wspace=0.35)
# Row 3: 2 columns (panels g–h)
gs_bot = outer[2].subgridspec(1, 2, wspace=0.30)

# Create axes
axes_a = [fig.add_subplot(gs_top[0, i]) for i in range(3)]
axes_b = [fig.add_subplot(gs_mid[0, i]) for i in range(3)]
ax_g = fig.add_subplot(gs_bot[0, 0])
ax_h = fig.add_subplot(gs_bot[0, 1])

# ------------------------------------------------------------------------------
# Panel A (a, b, c) – unchanged
# ------------------------------------------------------------------------------
for i, diff_label in enumerate(diff_labels):
    ax = axes_a[i]
    subset = df_summary[df_summary["difference_label"] == diff_label]

    x_positions = np.arange(len(ecosystem_order))
    means, cis = [], []
    for eco in ecosystem_order:
        row = subset[subset["water_class"] == eco]
        if len(row) > 0:
            means.append(row["mean_difference"].values[0])
            cis.append(row["ci95_difference"].values[0])
        else:
            means.append(np.nan)
            cis.append(np.nan)

    for j, eco in enumerate(ecosystem_order):
        if not np.isnan(means[j]):
            ax.errorbar(x_positions[j], means[j], yerr=cis[j], fmt='o',
                        color=ecosystem_colours[eco], capsize=8, elinewidth=3,
                        markersize=12, markeredgecolor='none', capthick=3)

    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)
    ax.set_xticks(x_positions)
    tick_labels = ax.set_xticklabels(ecosystem_order, fontsize=26)
    for tick, eco in zip(tick_labels, ecosystem_order):
        tick.set_color(ecosystem_colours[eco])
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=24, labelpad=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.text(-0.15, 1.05, f"{chr(97+i)})", transform=ax.transAxes,
            fontsize=34, fontweight='bold', va='bottom', ha='left')
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)

# ------------------------------------------------------------------------------
# Panel B (d, e, f) – unchanged
# ------------------------------------------------------------------------------
for i, diff_label in enumerate(diff_labels):
    ax = axes_b[i]
    subset = df_smooth[df_smooth["difference_label"] == diff_label]

    for eco in ecosystem_order:
        eco_sub = subset[subset["ecosystem_label"] == eco]
        if eco_sub.empty:
            continue
        eco_sub = eco_sub.sort_values("Trans_ratio")
        x = eco_sub["Trans_ratio"].values
        y = eco_sub["prediction"].values
        lower = eco_sub["lower"].values
        upper = eco_sub["upper"].values
        color = ecosystem_colours[eco]

        ax.fill_between(x, lower, upper, color=color, alpha=0.15)
        ax.plot(x, y, color=color, linewidth=3)

    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)
    ax.set_xlabel("T:ET ratio", fontsize=28)
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=28, labelpad=10)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.text(-0.15, 1.02, f"{chr(100 + i)})", transform=ax.transAxes,
            fontsize=34, fontweight='bold', va='bottom', ha='left')

    # edf box
    lines = []
    for eco in ecosystem_order:
        edf_val = edf_lookup.get((diff_label, eco), np.nan)
        if not np.isnan(edf_val):
            lines.append((f"{eco}: edf = {edf_val:.2f}", ecosystem_colours[eco]))
        else:
            lines.append((f"{eco}: edf = NA", ecosystem_colours[eco]))
    x0, y0 = 0.90, 0.05
    line_height = 0.06
    rect_width = 0.70
    rect_height = len(lines) * line_height + 0.04
    rect = patches.Rectangle((x0 - rect_width, y0), rect_width, rect_height,
                             transform=ax.transAxes, facecolor='white',
                             edgecolor='black', linewidth=1.2, alpha=0.85, clip_on=False)
    ax.add_patch(rect)
    for j, (line, color) in enumerate(lines):
        ax.text(x0 - 0.02, y0 + rect_height - j * line_height - 0.02, line,
                transform=ax.transAxes, ha='right', va='top', fontsize=23,
                color=color, clip_on=False)

    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)

# Shared R² label for row 2 – place above panels d–f
mid_pos = outer[1].get_position(fig)
label_x = (mid_pos.x0 + mid_pos.x1) / 2
label_y = mid_pos.y1 + 0.01   # slight offset above
fig.text(label_x, label_y, model_text, ha="center", va="bottom", fontsize=24,
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                   edgecolor="black", linewidth=1.5, alpha=0.9))

# ------------------------------------------------------------------------------
# Panel g: Ecosystem – violin plot (capped per group) with background points
# ------------------------------------------------------------------------------
positions = np.arange(len(eco_order_gh))
data_violin = [wue_monthly[wue_monthly['water_class'] == e]['WUE_T_capped'].dropna().values for e in eco_order_gh]

# Violin plot (capped for visual)
violin_parts = ax_g.violinplot(data_violin, positions=positions, showmeans=False, showmedians=False, widths=0.8)
for i, pc in enumerate(violin_parts['bodies']):
    pc.set_facecolor(ecosystem_colours[eco_order_gh[i]])
    pc.set_edgecolor(ecosystem_colours[eco_order_gh[i]])
    pc.set_alpha(0.7)

# Add jittered background points (not dim) – use capped values for consistency
for i, eco in enumerate(eco_order_gh):
    vals = wue_monthly[wue_monthly['water_class'] == eco]['WUE_T_capped'].dropna().values
    x_jitter = np.random.normal(i, 0.06, size=len(vals))
    ax_g.scatter(x_jitter, vals, alpha=0.5, s=20, color=ecosystem_colours[eco],
                 edgecolors='none')  # no edge to keep clean

# Set x‑axis limits
ax_g.set_xlim(-0.5, len(eco_order_gh) - 0.5)

# Adjust y‑axis top to give room for mean labels
max_data = max([np.max(vals) for vals in data_violin]) if data_violin else 1
max_label = eco_means.max() if not eco_means.empty else 1
ylim_top = max(max_data, max_label) * 1.15
ax_g.set_ylim(0, ylim_top)
ax_g.locator_params(axis='y', nbins=3)
# Median lines (black) – using uncapped monthly medians
for i, eco in enumerate(eco_order_gh):
    med = eco_medians[eco]
    ax_g.plot([i - 0.15, i + 0.15], [med, med],
              color='black', linewidth=2.5, solid_capstyle='butt')

# Mean labels (uncapped site‑level means)
ymin_g, ymax_g = ax_g.get_ylim()
for i, eco in enumerate(eco_order_gh):
    mean_val = eco_means[eco]
    ax_g.text(i, ymax_g - 0.02 * (ymax_g - ymin_g), f"{mean_val:.2f}",
              ha='center', va='top', fontsize=24, color='black', fontweight='bold')

# p‑value box – inside panel, near top-left
ax_g.text(0.45, 0.80, p_text_eco, transform=ax_g.transAxes,
          fontsize=24, ha='right', va='top',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

# x‑axis labels and colours
ax_g.set_xticks(positions)
ax_g.set_xticklabels(eco_order_gh)
for tick, eco in zip(ax_g.get_xticklabels(), eco_order_gh):
    tick.set_color(ecosystem_colours[eco])

ax_g.set_ylabel('WUE$_T$ (g C kg$^{-1}$ H$_2$O$^{-1}$)', fontsize=28, labelpad=10)
ax_g.set_xlabel('')
ax_g.set_title('')
ax_g.text(-0.15, 1.02, 'g)', transform=ax_g.transAxes,
          fontsize=34, fontweight='bold', va='bottom', ha='left')
for spine in ax_g.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)
ax_g.tick_params(axis='both', colors='black', width=2, length=8)
ax_g.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ------------------------------------------------------------------------------
# Panel h: Coast – mean ± 95% CI with background site‑level points
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

# Mean labels – placed inside panel near top
#ymin_h, ymax_h = 0, 1.05  # fixed for T:ET
ymin_h, ymax_h = 0.2, 0.8  # fixed for T:ET
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

# p‑value box – bottom right
ax_h.text(0.97, 0.05, p_text_coast, transform=ax_h.transAxes,
          fontsize=24, ha='right', va='bottom',
          bbox=dict(boxstyle="round,pad=0.4", facecolor='white', edgecolor='black', linewidth=1.2, alpha=0.9))

for spine in ax_h.spines.values():
    spine.set_color('black')
    spine.set_linewidth(2)
ax_h.tick_params(axis='both', colors='black', width=2, length=8)
ax_h.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ------------------------------------------------------------------------------
# Save figure (without bbox_inches='tight' to preserve dimensions)
# ------------------------------------------------------------------------------
fig.savefig(fig_output, dpi=600, facecolor='white')
print(f"Combined figure saved to: {fig_output}")