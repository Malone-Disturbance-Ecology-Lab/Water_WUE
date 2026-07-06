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
coast_colours = {
    "Atlantic Coast": "#2E8B57",
    "Pacific Coast": "#DC143C",
    "Gulf Coast":     "#00CED1",
    "AK Coast":       "#8B4513"
}
coast_labels = {
    "Atlantic Coast": "Atlantic",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf",
    "AK Coast": "AK"
}

df["coast_region"] = pd.Categorical(df["coast_region"], categories=coast_order, ordered=True)
df["water_class"] = pd.Categorical(df["water_class"], categories=ecosystem_order, ordered=True)

np.random.seed(42)

# ------------------------------------------------------------------------------
# Figure setup
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(24, 17))

plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 32,
    "xtick.labelsize": 35,
    "ytick.labelsize": 35,
})

for row in axes:
    for ax in row:
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(2)
        ax.tick_params(axis='both', colors='black', width=2, length=8)
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

# ==============================================================================
# Panel a – Stability vs T:ET by ecosystem
# ==============================================================================
ax_a = axes[0, 0]
for eco in ecosystem_order:
    sub = df[df["water_class"] == eco]
    ax_a.scatter(sub["mean_TET"], sub["stability"],
                 color=ecosystem_colours[eco], s=80, alpha=1.0,
                 label=eco, edgecolors='none')
    sub2 = sub.dropna(subset=["mean_TET", "stability"])
    if len(sub2) >= 3:
        slope, intercept, r, p, se = linregress(sub2["mean_TET"], sub2["stability"])
        x_vals = np.linspace(sub2["mean_TET"].min(), sub2["mean_TET"].max(), 50)
        ax_a.plot(x_vals, slope * x_vals + intercept,
                  color=ecosystem_colours[eco], linewidth=2.5, linestyle='--')
ax_a.set_xlabel("Mean T:ET ratio")
ax_a.set_ylabel("Stability")
legend_a = ax_a.legend(loc='best', fontsize=28, labelspacing=0.2)
for text, eco in zip(legend_a.get_texts(), ecosystem_order):
    text.set_color(ecosystem_colours[eco])
ax_a.text(-0.20, 1.06, 'a)', transform=ax_a.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_a.text(0.5, 1.05, f'ρ = {rho:.2f}; {p_stab_text}', transform=ax_a.transAxes,
          fontsize=30, fontweight='bold', va='bottom', ha='center')
ax_a.xaxis.set_major_locator(MaxNLocator(5))
ax_a.yaxis.set_major_locator(MaxNLocator(5))

# ==============================================================================
# Panel b – Stability vs T:ET by coast
# ==============================================================================
ax_b = axes[0, 1]
for coast in coast_order:
    sub = df[df["coast_region"] == coast]
    ax_b.scatter(sub["mean_TET"], sub["stability"],
                 color=coast_colours[coast], s=80, alpha=1.0,
                 label=coast_labels[coast], edgecolors='none')
    sub2 = sub.dropna(subset=["mean_TET", "stability"])
    if len(sub2) >= 3:
        slope, intercept, r, p, se = linregress(sub2["mean_TET"], sub2["stability"])
        x_vals = np.linspace(sub2["mean_TET"].min(), sub2["mean_TET"].max(), 50)
        ax_b.plot(x_vals, slope * x_vals + intercept,
                  color=coast_colours[coast], linewidth=2.5, linestyle='--')
ax_b.set_xlabel("Mean T:ET ratio")
ax_b.set_ylabel("Stability")
legend_b = ax_b.legend(loc='best', fontsize=28, labelspacing=0.2)
for text, coast in zip(legend_b.get_texts(), coast_order):
    text.set_color(coast_colours[coast])
ax_b.text(-0.20, 1.06, 'b)', transform=ax_b.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_b.text(0.5, 1.05, f'{p_stab_coast_text}', transform=ax_b.transAxes,
          fontsize=36, fontweight='bold', va='bottom', ha='center')
ax_b.xaxis.set_major_locator(MaxNLocator(5))
ax_b.yaxis.set_major_locator(MaxNLocator(5))

# ==============================================================================
# Panel c – SPEI‑3 plasticity slope by coast
# ==============================================================================
ax_c = axes[0, 2]
positions = np.arange(1, len(coast_order) + 1)

for i, coast in enumerate(coast_order):
    vals = df[df["coast_region"] == coast]["plasticity_slope_SPEI3"].dropna().values
    x_jitter = np.random.normal(i + 1, 0.06, size=len(vals))
    ax_c.scatter(x_jitter, vals, alpha=0.5, s=50,
                 color=coast_colours[coast], edgecolors='none')

for i, coast in enumerate(coast_order):
    vals = df[df["coast_region"] == coast]["plasticity_slope_SPEI3"].dropna().values
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
ax_c.set_ylabel("SPEI-3 plasticity slope")
ax_c.text(-0.20, 1.06, 'c)', transform=ax_c.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_c.text(0.5, 1.05, f'{p_slope_text}', transform=ax_c.transAxes,
          fontsize=36, fontweight='bold', va='bottom', ha='center')
ax_c.yaxis.set_major_locator(MaxNLocator(5))
ax_c.tick_params(axis='x', labelsize=28)

# ==============================================================================
# Panel d – Plasticity range by coast – fixed: removed x-axis locator
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
    ax_d.scatter(x_jitter, vals, alpha=0.6, s=50,
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

ax_d.set_xlim(0.5, len(coast_order) + 0.5)   # ensure labels are inside
ax_d.set_ylabel("Plasticity range (95th / 5th)")
ax_d.text(-0.20, 1.06, 'd)', transform=ax_d.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_d.text(0.5, 1.05, f'{p_range_text}', transform=ax_d.transAxes,
          fontsize=36, fontweight='bold', va='bottom', ha='center')

# Keep y-axis ticks, but DO NOT set x-axis locator
ax_d.yaxis.set_major_locator(MaxNLocator(5))
# ax_d.xaxis.set_major_locator(MaxNLocator(5))   # REMOVED – this was causing label loss
ax_d.tick_params(axis='x', labelsize=28)

# ==============================================================================
# Panel e – Drought resistance (overall)
# ==============================================================================
ax_e = axes[1, 1]
vals_e = df["resistance"].dropna().values
x_jitter_e = np.random.normal(1, 0.08, size=len(vals_e))
ax_e.scatter(x_jitter_e, vals_e, alpha=0.6, s=50,
             color='gray', edgecolors='none')

bp_e = ax_e.boxplot(vals_e, positions=[1], widths=0.4, patch_artist=True,
                    showmeans=False, showfliers=False,
                    boxprops=dict(linewidth=2, color='gray', facecolor='gray', alpha=0.25),
                    whiskerprops=dict(linewidth=2, color='gray'),
                    capprops=dict(linewidth=2, color='gray'),
                    medianprops=dict(linewidth=3, color='black'))
ax_e.axhline(1, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_e.set_xticks([1])
ax_e.set_xticklabels(['Sites\n(N = 39)'])
ax_e.set_xlabel("")
ax_e.set_ylabel("Drought resistance")
ax_e.text(-0.20, 1.06, 'e)', transform=ax_e.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_e.yaxis.set_major_locator(MaxNLocator(5))

# ==============================================================================
# Panel f – Drought recovery (overall)
# ==============================================================================
ax_f = axes[1, 2]
vals_f = df["mean_recovery"].dropna().values
x_jitter_f = np.random.normal(1, 0.08, size=len(vals_f))
ax_f.scatter(x_jitter_f, vals_f, alpha=0.6, s=50,
             color='gray', edgecolors='none')

bp_f = ax_f.boxplot(vals_f, positions=[1], widths=0.4, patch_artist=True,
                    showmeans=False, showfliers=False,
                    boxprops=dict(linewidth=2, color='gray', facecolor='gray', alpha=0.25),
                    whiskerprops=dict(linewidth=2, color='gray'),
                    capprops=dict(linewidth=2, color='gray'),
                    medianprops=dict(linewidth=3, color='black'))
ax_f.axhline(1, color='gray', linestyle='--', linewidth=2, alpha=0.6)
ax_f.set_xticks([1])
ax_f.set_xticklabels(['Sites\n(N = 27)'])
ax_f.set_xlabel("")
ax_f.set_ylabel("Drought recovery")
ax_f.text(-0.20, 1.06, 'f)', transform=ax_f.transAxes,
          fontsize=40, fontweight='bold', va='bottom', ha='left')
ax_f.yaxis.set_major_locator(MaxNLocator(5))

# ------------------------------------------------------------------------------
# Adjust spacing
# ------------------------------------------------------------------------------
plt.subplots_adjust(left=0.08, right=0.95, top=0.90, bottom=0.25,
                    wspace=0.25, hspace=0.45)

save_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_figures"
os.makedirs(save_dir, exist_ok=True)
fig_output = os.path.join(save_dir, "Q2_performance_figure.png")
fig.savefig(fig_output, dpi=600, bbox_inches='tight', pad_inches=0.5, facecolor='white')
print(f"Figure saved to: {fig_output}")

plt.show()