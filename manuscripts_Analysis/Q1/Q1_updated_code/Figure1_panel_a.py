import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)

summary_file = os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv")
fig_output = os.path.join(figure_dir, "Q1_panel_A_summary.png")

# ------------------------------------------------------------------------------
# Load summary data
# ------------------------------------------------------------------------------
df = pd.read_csv(summary_file)

required = ["difference_label", "water_class", "mean_difference", "ci95_difference"]
if not all(col in df.columns for col in required):
    raise ValueError("Missing required columns in summary CSV")

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
ecosystem_order = ["Upland", "Freshwater", "Saline"]
ecosystem_colours = {
    "Upland": "#800080",
    "Freshwater": "#0000FF",
    "Saline": "#FFA500",
}

# Two‑line y‑axis labels (difference + unit)
diff_label_mapping = {
    "WUE_ET - WUE_T": r"WUE$_{ET}$ – WUE$_{T}$\n(g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_ET - WUE_E": r"WUE$_{ET}$ – WUE$_{E}$\n(g C kg$^{-1}$ H$_2$O$^{-1}$)",
    "WUE_E - WUE_T": r"WUE$_{E}$ – WUE$_{T}$\n(g C kg$^{-1}$ H$_2$O$^{-1}$)",
}
diff_labels = list(diff_label_mapping.keys())

df = df[df["difference_label"].isin(diff_labels)]
df["difference_label"] = pd.Categorical(df["difference_label"], categories=diff_labels, ordered=True)
df["water_class"] = pd.Categorical(df["water_class"], categories=ecosystem_order, ordered=True)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(24, 6), sharey=False)

# Set font sizes
plt.rcParams.update({
    "font.size": 28,
    "axes.labelsize": 28,
    "xtick.labelsize": 26,
    "ytick.labelsize": 26,
})

# Make axes spines and tick lines black
for ax in axes:
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(2)
    ax.tick_params(axis='both', colors='black', width=2, length=8)

# Loop over each difference type
for i, diff_label in enumerate(diff_labels):
    ax = axes[i]
    subset = df[df["difference_label"] == diff_label]

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

    # Plot error bars – no black outlines
    for j, eco in enumerate(ecosystem_order):
        if not np.isnan(means[j]):
            ax.errorbar(
                x_positions[j], means[j],
                yerr=cis[j],
                fmt='o',
                color=ecosystem_colours[eco],
                capsize=8,
                elinewidth=3,
                markersize=12,
                markeredgecolor='none',
                capthick=3,
            )

    # Red dashed zero line
    ax.axhline(0, color='red', linestyle='--', linewidth=3, alpha=0.8)

    # X‑axis: set tick labels with ecosystem colours
    ax.set_xticks(x_positions)
    tick_labels = ax.set_xticklabels(ecosystem_order, fontsize=26)
    for tick, eco in zip(tick_labels, ecosystem_order):
        tick.set_color(ecosystem_colours[eco])   # colour each label to match ecosystem

    ax.set_xlabel("")

    # Y‑axis label
    ax.set_ylabel(diff_label_mapping[diff_label], fontsize=24, labelpad=10)

    ax.set_title("")
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Panel label (a), (b), (c)
    ax.text(-0.15, 1.02, f"{chr(97+i)})", transform=ax.transAxes,
            fontsize=34, fontweight='bold', va='bottom', ha='left')

# Adjust margins with more space between panels
plt.subplots_adjust(left=0.12, right=0.92, top=0.95, bottom=0.15, wspace=0.3)

# Save high‑resolution figure
fig.savefig(fig_output, dpi=600, bbox_inches='tight', pad_inches=0.5, facecolor='white')
print(f"Panel A figure saved to: {fig_output}")