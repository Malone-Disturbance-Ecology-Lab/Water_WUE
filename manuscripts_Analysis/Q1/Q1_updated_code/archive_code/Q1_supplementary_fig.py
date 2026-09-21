# -*- coding: utf-8 -*-
"""
Supplementary Figure: Q1 GAM comparison – 5‑panel version.
Top: Linear GAM predictions (3 subpanels labeled a, b, c).
Bottom left (d): Two ΔAIC bars.
Bottom right (e): Adjusted R² bar plot (gray, 2 decimals).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
gam_comp_dir = os.path.join(output_dir, "GAM_comparison")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
fig_output = os.path.join(figure_dir, "Supplementary_Q1_GAM_model_comparison.png")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
linear_pred = pd.read_csv(os.path.join(output_dir, "Q1_final_linear_gam_predictions_TET.csv"))

perf_linear = pd.read_csv(os.path.join(output_dir, "Q1_final_gam_model_comparison.csv"))
linear_row = perf_linear[perf_linear["model"] == "q1_linear_gam"].iloc[0]

perf_shared = pd.read_csv(os.path.join(gam_comp_dir, "Q1_GAM_A_B_C_model_performance.csv"))
shared_row = perf_shared[perf_shared["model"] == "B_mean_interaction_common_smooth"].iloc[0]

separate_row = perf_linear[perf_linear["model"] == "q1_smooth_gam"].iloc[0]

aic = {
    "linear": linear_row["AIC"],
    "shared": shared_row["AIC"],
    "separate": separate_row["AIC"]
}
r2 = {
    "linear": linear_row["adjusted_r_squared"],
    "shared": shared_row["adjusted_r_squared"],
    "separate": separate_row["adjusted_r_squared"]
}

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500"
}
diff_label_mapping = {
    "WUE_ET - WUE_T": r"WUE$_{ET}$ – WUE$_{T}$",
    "WUE_ET - WUE_E": r"WUE$_{ET}$ – WUE$_{E}$",
    "WUE_E - WUE_T": r"WUE$_{E}$ – WUE$_{T}$",
}
diff_labels = list(diff_label_mapping.keys())

# ------------------------------------------------------------------------------
# Panel d: two ΔAIC bars (short labels + note)
# ------------------------------------------------------------------------------
def plot_panel_d(ax):
    """
    Panel d: ΔAIC improvement toward the 9-curve nonlinear GAM.
    Short x‑labels and a clarifying note.
    """
    aic_linear = aic["linear"]
    aic_shared = aic["shared"]
    aic_separate = aic["separate"]

    delta_aic_lin_sep = aic_linear - aic_separate
    delta_aic_sha_sep = aic_shared - aic_separate

    x_pos = np.arange(2)
    delta_vals = [delta_aic_lin_sep, delta_aic_sha_sep]

    labels = [
        "Linear GAM\nvs Non-linear GAM (9-curve)",
        "1-curve\nvs 9-curve GAM"
    ]

    bars = ax.bar(
        x_pos,
        delta_vals,
        width=0.55,
        color=["#1f77b4", "#ff7f0e"],
        edgecolor="black",
        linewidth=1.2
    )

    ymax = max(delta_vals) * 1.18
    ax.set_ylim(0, ymax)

    for xi, yi in zip(x_pos, delta_vals):
        ax.text(
            xi,
            yi + ymax * 0.02,
            f"ΔAIC = {yi:.1f}",
            ha="center",
            va="bottom",
            fontsize=22,
            fontweight="bold"
        )

    ax.set_ylabel("ΔAIC improvement", fontsize=26, color="black")
    ax.tick_params(axis="y", labelsize=24, colors="black")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=18)

    ax.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="y")

    # Panel label
    ax.text(
        -0.12,
        1.04,
        "d)",
        transform=ax.transAxes,
        fontsize=32,
        fontweight="bold",
        va="bottom",
        ha="left"
    )

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.2)

# ------------------------------------------------------------------------------
# Panel e: adjusted R² bar plot (all bars gray, 2 decimals)
# ------------------------------------------------------------------------------
def plot_panel_e(ax):
    """
    Panel e: adjusted R² for the three models (all bars gray).
    R² values shown with 2 decimal places.
    """
    r2_linear = r2["linear"]
    r2_shared = r2["shared"]
    r2_separate = r2["separate"]

    x_pos = np.arange(3)
    r2_vals = [r2_linear, r2_shared, r2_separate]

    labels = [
        "Linear\nGAM",
        "Nonlinear GAM\n1 curve",
        "Nonlinear GAM\n9 curves"
    ]

    ax.bar(
        x_pos,
        r2_vals,
        width=0.55,
        color=["#A9A9A9", "#A9A9A9", "#A9A9A9"],
        edgecolor="black",
        linewidth=1.2
    )

    for xi, yi in zip(x_pos, r2_vals):
        ax.text(
            xi,
            yi + 0.015,
            f"{yi:.2f}",
            ha="center",
            va="bottom",
            fontsize=22,
            fontweight="bold"
        )

    ax.set_ylabel("Adjusted R²", fontsize=26, color="black")
    ax.tick_params(axis="y", labelsize=24, colors="black")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=18)

    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3, linestyle="-", linewidth=0.5, axis="y")

    # Panel label
    ax.text(
        -0.12,
        1.04,
        "e)",
        transform=ax.transAxes,
        fontsize=32,
        fontweight="bold",
        va="bottom",
        ha="left"
    )

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.2)

# ------------------------------------------------------------------------------
# Figure setup
# ------------------------------------------------------------------------------
fig = plt.figure(figsize=(20, 12))

outer = gridspec.GridSpec(
    2, 1,
    height_ratios=[1.15, 1.0],
    hspace=0.45,
    left=0.10,
    right=0.95,
    top=0.89,
    bottom=0.10
)

# Top row: 3 columns
gs_top = gridspec.GridSpecFromSubplotSpec(
    1, 3,
    subplot_spec=outer[0],
    wspace=0.30
)

# Bottom row: 2 columns
gs_bottom = gridspec.GridSpecFromSubplotSpec(
    1, 2,
    subplot_spec=outer[1],
    wspace=0.40
)

top_axes = [fig.add_subplot(gs_top[0, i]) for i in range(3)]

# ------------------------------------------------------------------------------
# Top row title (centered) – moved slightly higher
# ------------------------------------------------------------------------------
fig.text(
    0.5,
    0.975,                     # <-- moved higher (was 0.965)
    "Linear GAM predictions",
    fontsize=28,
    fontweight="bold",
    va="top",
    ha="center"
)

# ------------------------------------------------------------------------------
# Top row: Linear GAM predictions (with outside labels a, b, c)
# ------------------------------------------------------------------------------
for i, diff in enumerate(diff_labels):
    ax = top_axes[i]
    sub = linear_pred[linear_pred["difference_label"] == diff]
    for eco in ecosystem_order:
        eco_sub = sub[sub["ecosystem_label"] == eco].sort_values("Trans_ratio")
        if eco_sub.empty:
            continue
        x = eco_sub["Trans_ratio"].values
        y = eco_sub["prediction"].values
        lower = eco_sub["lower"].values
        upper = eco_sub["upper"].values
        color = ecosystem_colours[eco]
        ax.fill_between(x, lower, upper, color=color, alpha=0.15)
        ax.plot(x, y, color=color, linewidth=2.5, label=eco if i == 2 else "")
    ax.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.8)
    ax.set_xlabel("T:ET ratio", fontsize=26)
    if i == 0:
        ax.set_ylabel("Predicted difference", fontsize=26)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Outside panel label (a, b, c)
    ax.text(
        -0.12,
        1.04,
        f"{chr(97+i)})",
        transform=ax.transAxes,
        fontsize=32,
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False
    )
    
    ax.tick_params(axis='both', labelsize=24)

# Single legend
handles = [plt.Line2D([0],[0], color=ecosystem_colours[eco], lw=2.5, label=eco) for eco in ecosystem_order]
leg = top_axes[-1].legend(handles=handles, loc='upper left', fontsize=22, frameon=True, edgecolor='black')
leg.get_frame().set_linewidth(1)

# ------------------------------------------------------------------------------
# Bottom row: d) and e)
# ------------------------------------------------------------------------------
ax_d = fig.add_subplot(gs_bottom[0, 0])
plot_panel_d(ax_d)

ax_e = fig.add_subplot(gs_bottom[0, 1])
plot_panel_e(ax_e)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
fig.savefig(fig_output, dpi=600, bbox_inches='tight', facecolor='white')
print(f"Supplementary figure saved to: {fig_output}")