import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from matplotlib.ticker import MaxNLocator

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------
data_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship\Q2_performance_by_coast_site_level_data.csv"
summary_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship\Q2_T_ET_performance_relationship_summary.csv"
coast_tests_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship\Q2_performance_by_coast_tests.csv"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
df = pd.read_csv(data_file)
df_summary = pd.read_csv(summary_file)
df_tests = pd.read_csv(coast_tests_file)

# ------------------------------------------------------------------------------
# Extract p-values and rho
# ------------------------------------------------------------------------------
stab_row = df_summary[df_summary["metric"] == "stability"].iloc[0]
rho = stab_row["spearman_rho"]
p_stab = stab_row["spearman_p_FDR"]
p_stab_text = "p < 0.001" if p_stab < 0.001 else f"p = {p_stab:.3f}"

p_stab_coast = df_tests[(df_tests["metric"] == "stability")]["kw_p_all_groups_FDR"].iloc[0]
p_stab_coast_text = "p < 0.001" if p_stab_coast < 0.001 else f"p = {p_stab_coast:.3f}"

p_slope = df_tests[(df_tests["metric"] == "plasticity_slope_SPEI3")]["kw_p_all_groups_FDR"].iloc[0]
p_slope_text = "p < 0.001" if p_slope < 0.001 else f"p = {p_slope:.3f}"

p_range = df_tests[(df_tests["metric"] == "plasticity_p95_p05")]["kw_p_all_groups_FDR"].iloc[0]
p_range_text = "p < 0.001" if p_range < 0.001 else f"p = {p_range:.3f}"

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500",
}

coast_order = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
# FINAL COAST COLOURS – consistent with all other figures
coast_colours = {
    "Atlantic Coast": "#A50F15",   # dark red
    "Pacific Coast":  "#0072B2",   # blue
    "Gulf Coast":     "#4D4D4D",   # dark charcoal
    "AK Coast":       "#009E73"    # bluish green
}
coast_labels = {
    "Atlantic Coast": "Atlantic",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf",
    "AK Coast": "Alaska"
}

df["coast_region"] = pd.Categorical(df["coast_region"], categories=coast_order, ordered=True)
df["water_class"] = pd.Categorical(df["water_class"], categories=ecosystem_order, ordered=True)

np.random.seed(42)

# ------------------------------------------------------------------------------
# Font sizes and rcParams – must be set BEFORE figure creation
# ------------------------------------------------------------------------------
TICK_FS = 28            # numerical tick labels
COAST_TICK_FS = 30      # coast name labels in panels c and d
SITE_TICK_FS = 26       # "Sites (N = ...)" labels in panels e and f

plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 32,
    "xtick.labelsize": TICK_FS,
    "ytick.labelsize": TICK_FS,
})

# ------------------------------------------------------------------------------
# Figure setup – standard 2x3 grid
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(26, 17))

for row in axes:
    for ax in row:
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(2)
        ax.tick_params(
            axis='both',
            which='major',
            colors='black',
            labelsize=TICK_FS,
            width=2,
            length=8
        )
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ==============================================================================
# Panel a – Stability vs T:ET by ecosystem
# ==============================================================================
ax_a = axes[0, 0]
for eco in ecosystem_order:
    sub = df[df["water_class"] == eco]
    ax_a.scatter(sub["mean_TET"], sub["stability"],
                 color=ecosystem_colours[eco], s=120, alpha=1.0,
                 label=eco, edgecolors='none')
    sub2 = sub.dropna(subset=["mean_TET", "stability"])
    if len(sub2) >= 3:
        slope, intercept, r, p, se = linregress(sub2["mean_TET"], sub2["stability"])
        x_vals = np.linspace(sub2["mean_TET"].min(), sub2["mean_TET"].max(), 50)
        ax_a.plot(x_vals, slope * x_vals + intercept,
                  color=ecosystem_colours[eco], linewidth=1.5, linestyle='--')

overall = df.dropna(subset=["mean_TET", "stability"])
if len(overall) >= 3:
    slope_all, intercept_all, _, _, _ = linregress(overall["mean_TET"], overall["stability"])
    x_all = np.linspace(overall["mean_TET"].min(), overall["mean_TET"].max(), 50)
    ax_a.plot(x_all, slope_all * x_all + intercept_all,
              color='black', linewidth=4, linestyle='-',
              label=f'Overall ({p_stab_text})')

ax_a.set_xlabel("Mean T:ET ratio")
ax_a.set_ylabel("Stability")
legend_a = ax_a.legend(loc='best', fontsize=24, labelspacing=0.2)
for text in legend_a.get_texts():
    if text.get_text() in ecosystem_order:
        text.set_color(ecosystem_colours[text.get_text()])
    else:
        text.set_color('black')
ax_a.text(-0.20, 1.06, 'a)', transform=ax_a.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_a.xaxis.set_major_locator(MaxNLocator(4))
ax_a.yaxis.set_major_locator(MaxNLocator(5))
ax_a.tick_params(axis='both', labelsize=TICK_FS)

# ==============================================================================
# Panel b – Stability vs T:ET by coast
# ==============================================================================
ax_b = axes[0, 1]
for coast in coast_order:
    sub = df[df["coast_region"] == coast]
    ax_b.scatter(sub["mean_TET"], sub["stability"],
                 color=coast_colours[coast], s=120, alpha=1.0,
                 label=coast_labels[coast], edgecolors='none')
    sub2 = sub.dropna(subset=["mean_TET", "stability"])
    if len(sub2) >= 3:
        slope, intercept, r, p, se = linregress(sub2["mean_TET"], sub2["stability"])
        x_vals = np.linspace(sub2["mean_TET"].min(), sub2["mean_TET"].max(), 50)
        ax_b.plot(x_vals, slope * x_vals + intercept,
                  color=coast_colours[coast], linewidth=2.5, linestyle='--')
ax_b.set_xlabel("Mean T:ET ratio")
ax_b.set_ylabel("Stability")
legend_b = ax_b.legend(loc='best', fontsize=24, labelspacing=0.2)
for text, coast in zip(legend_b.get_texts(), coast_order):
    text.set_color(coast_colours[coast])
ax_b.text(-0.20, 1.06, 'b)', transform=ax_b.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_b.text(0.02, 0.02, p_stab_coast_text, transform=ax_b.transAxes,
          fontsize=28, fontweight='bold', va='bottom', ha='left')
ax_b.xaxis.set_major_locator(MaxNLocator(4))
ax_b.yaxis.set_major_locator(MaxNLocator(5))
ax_b.tick_params(axis='both', labelsize=TICK_FS)

# ==============================================================================
# Panel c – SPEI‑3 plasticity slope by coast
# ==============================================================================
ax_c = axes[0, 2]
positions = np.arange(1, len(coast_order) + 1)

for i, coast in enumerate(coast_order):
    vals = df[df["coast_region"] == coast]["plasticity_slope_SPEI3"].dropna().values
    coast_p2 = np.percentile(vals, 2)
    capped_vals = np.maximum(vals, coast_p2)
    x_jitter = np.random.normal(i + 1, 0.06, size=len(vals))
    ax_c.scatter(x_jitter, capped_vals, alpha=0.5, s=140,
                 color=coast_colours[coast], edgecolors='none')
    bp = ax_c.boxplot(vals, positions=[i+1], widths=0.35, patch_artist=True,
                      showmeans=False, showfliers=False,
                      boxprops=dict(linewidth=2, color=coast_colours[coast],
                                    facecolor=coast_colours[coast], alpha=0.3),
                      whiskerprops=dict(linewidth=2, color=coast_colours[coast]),
                      capprops=dict(linewidth=2, color=coast_colours[coast]),
                      medianprops=dict(linewidth=3, color=coast_colours[coast]))

ax_c.axhline(0, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_c.set_xticks(positions)
ax_c.set_xticklabels([coast_labels[c] for c in coast_order], rotation=45, ha='right')
for tick, coast in zip(ax_c.get_xticklabels(), coast_order):
    tick.set_color(coast_colours[coast])
    tick.set_fontsize(COAST_TICK_FS)
    tick.set_fontweight('normal')
ax_c.set_ylabel("SPEI-3 plasticity slope")
ax_c.text(-0.20, 1.06, 'c)', transform=ax_c.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_c.text(0.98, 0.02, p_slope_text, transform=ax_c.transAxes,
          fontsize=28, fontweight='bold', va='bottom', ha='right')
ax_c.yaxis.set_major_locator(MaxNLocator(5))
ax_c.tick_params(axis='both', labelsize=TICK_FS)

# ==============================================================================
# Panel d – Plasticity range by coast
# ==============================================================================
ax_d = axes[1, 0]

violin_parts = ax_d.violinplot(
    [df[df["coast_region"] == coast]["plasticity_p95_p05"].dropna().values
     for coast in coast_order],
    positions=positions,
    showmeans=False,
    showmedians=False,
    showextrema=False,
    widths=0.7
)
for i, pc in enumerate(violin_parts['bodies']):
    color = coast_colours[coast_order[i]]
    pc.set_facecolor(color)
    pc.set_alpha(0.15)
    pc.set_edgecolor('none')

for i, coast in enumerate(coast_order):
    vals = df[df["coast_region"] == coast]["plasticity_p95_p05"].dropna().values
    x_jitter = np.random.normal(i + 1, 0.06, size=len(vals))
    ax_d.scatter(x_jitter, vals, alpha=0.6, s=140,
                 color=coast_colours[coast], edgecolors='none')

for i, coast in enumerate(coast_order):
    vals = df[df["coast_region"] == coast]["plasticity_p95_p05"].dropna().values
    bp = ax_d.boxplot(vals, positions=[i+1], widths=0.25, patch_artist=True,
                      showmeans=False, showfliers=False,
                      boxprops=dict(linewidth=2, color=coast_colours[coast],
                                    facecolor=coast_colours[coast], alpha=0.25),
                      whiskerprops=dict(linewidth=2, color=coast_colours[coast]),
                      capprops=dict(linewidth=2, color=coast_colours[coast]))
    for med in bp['medians']:
        med.set_color(coast_colours[coast])
        med.set_linewidth(3)

ax_d.set_xticks(positions)
ax_d.set_xticklabels([coast_labels[c] for c in coast_order], rotation=45, ha='right')
for tick, coast in zip(ax_d.get_xticklabels(), coast_order):
    tick.set_color(coast_colours[coast])
    tick.set_fontsize(COAST_TICK_FS)
    tick.set_fontweight('normal')

ax_d.set_xlim(0.5, len(coast_order) + 0.5)
ax_d.set_ylabel("Plasticity range")
ax_d.text(-0.20, 1.06, 'd)', transform=ax_d.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_d.text(0.02, 0.95, p_range_text, transform=ax_d.transAxes,
          fontsize=28, fontweight='bold', va='top', ha='left')
ax_d.yaxis.set_major_locator(MaxNLocator(5))
ax_d.tick_params(axis='both', labelsize=TICK_FS)

# ==============================================================================
# Panel e – Drought resistance
# ==============================================================================
ax_e = axes[1, 1]

x_jitter_e = np.random.normal(1, 0.08, size=len(df))

for coast in coast_order:
    mask = df["coast_region"] == coast
    vals = df.loc[mask, "resistance"].dropna().values
    x_vals = x_jitter_e[mask][~pd.isna(df.loc[mask, "resistance"])]
    ax_e.scatter(x_vals, vals, alpha=0.6, s=140,
                 color=coast_colours[coast], edgecolors='none', label=coast_labels[coast])

vals_e_all = df["resistance"].dropna().values
bp_e = ax_e.boxplot(vals_e_all, positions=[1], widths=0.4, patch_artist=True,
                    showmeans=False, showfliers=False,
                    boxprops=dict(linewidth=2, color='gray', facecolor='gray', alpha=0.25),
                    whiskerprops=dict(linewidth=2, color='gray'),
                    capprops=dict(linewidth=2, color='gray'),
                    medianprops=dict(linewidth=3, color='black'))
ax_e.axhline(1, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_e.set_xticks([1])
ax_e.set_xticklabels(['Sites\n(N = 39)'], fontsize=SITE_TICK_FS)
ax_e.set_xlabel("")
ax_e.set_ylabel("Drought resistance")
ax_e.text(-0.20, 1.06, 'e)', transform=ax_e.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_e.yaxis.set_major_locator(MaxNLocator(5))
ax_e.tick_params(axis='x', labelsize=SITE_TICK_FS, colors='black')
ax_e.tick_params(axis='y', labelsize=TICK_FS)

# ==============================================================================
# Panel f – Drought recovery
# ==============================================================================
ax_f = axes[1, 2]

x_jitter_f = np.random.normal(1, 0.08, size=len(df))

for coast in coast_order:
    mask = df["coast_region"] == coast
    vals = df.loc[mask, "mean_recovery"].dropna().values
    x_vals = x_jitter_f[mask][~pd.isna(df.loc[mask, "mean_recovery"])]
    ax_f.scatter(x_vals, vals, alpha=0.6, s=140,
                 color=coast_colours[coast], edgecolors='none', label=coast_labels[coast])

vals_f_all = df["mean_recovery"].dropna().values
bp_f = ax_f.boxplot(vals_f_all, positions=[1], widths=0.4, patch_artist=True,
                    showmeans=False, showfliers=False,
                    boxprops=dict(linewidth=2, color='gray', facecolor='gray', alpha=0.25),
                    whiskerprops=dict(linewidth=2, color='gray'),
                    capprops=dict(linewidth=2, color='gray'),
                    medianprops=dict(linewidth=3, color='black'))
ax_f.axhline(1, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_f.set_xticks([1])
ax_f.set_xticklabels(['Sites\n(N = 27)'], fontsize=SITE_TICK_FS)
ax_f.set_xlabel("")
ax_f.set_ylabel("Drought recovery")
ax_f.text(-0.20, 1.06, 'f)', transform=ax_f.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_f.yaxis.set_major_locator(MaxNLocator(5))
ax_f.tick_params(axis='x', labelsize=SITE_TICK_FS, colors='black')
ax_f.tick_params(axis='y', labelsize=TICK_FS)

# ------------------------------------------------------------------------------
# Final tick-label size and colour enforcement
# ------------------------------------------------------------------------------
for ax in [ax_a, ax_b, ax_c, ax_d, ax_e, ax_f]:
    ax.tick_params(axis='y', which='major', labelsize=TICK_FS, colors='black')

for ax in [ax_a, ax_b]:
    ax.tick_params(axis='x', which='major', labelsize=TICK_FS, colors='black')

# Coloured coastline labels in panels c and d
for ax in [ax_c, ax_d]:
    for tick, coast in zip(ax.get_xticklabels(), coast_order):
        tick.set_color(coast_colours[coast])
        tick.set_fontsize(COAST_TICK_FS)

# Site labels in panels e and f
ax_e.tick_params(axis='x', labelsize=SITE_TICK_FS, colors='black')
ax_f.tick_params(axis='x', labelsize=SITE_TICK_FS, colors='black')

# ------------------------------------------------------------------------------
# Adjust spacing and save
# ------------------------------------------------------------------------------
plt.subplots_adjust(left=0.06, right=0.95, top=0.92, bottom=0.13,
                    wspace=0.28, hspace=0.50)

save_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_figures"
os.makedirs(save_dir, exist_ok=True)
fig_output = os.path.join(save_dir, "Q2_performance_figure.png")
fig.savefig(fig_output, dpi=600, bbox_inches='tight', facecolor='white')
print(f"Figure saved to: {fig_output}")

plt.show()