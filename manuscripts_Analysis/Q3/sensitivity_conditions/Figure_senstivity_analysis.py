# Q3_sensitivity_four_panel_figure.py
# -----------------------------------------------------------------------------
# Figure script for Q3 reference-condition sensitivity analysis
# Produces only the combined 2x2 four-panel figure.
# Uses saved CSVs from the diagnostic; does NOT refit the GAM.
# Output: one PNG file, no separate panel images.
# -----------------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
results_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs\sensitivity_conditions_results"
output_file = os.path.join(results_dir, "Q3_sensitivity_four_panel.png")

# Ecosystem colors (original Q3 palette)
ecosystem_palette = {
    "Upland": "#800080",      # purple
    "Freshwater": "#0000FF",  # blue
    "Saline": "#FFA500"       # orange
}
water_class_order = ["Upland", "Freshwater", "Saline"]
coast_order_long = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]
coast_order_short = ["Atlantic", "Pacific", "Gulf", "AK"]
coast_mapping = dict(zip(coast_order_long, coast_order_short))
spei_order = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]

# Figure size and fonts
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 8
sns.set_style("whitegrid")

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
pred_curves = pd.read_csv(os.path.join(results_dir, "prediction_curves_all_reference_conditions.csv"))
resp_mag = pd.read_csv(os.path.join(results_dir, "response_magnitude_all_reference_conditions.csv"))
thresh_all = pd.read_csv(os.path.join(results_dir, "threshold_summary_all_reference_conditions.csv"))

# Convert month_f to numeric and set categorical order
for df in [pred_curves, resp_mag, thresh_all]:
    df["month_f"] = df["month_f"].astype(int)
    df["coast_region"] = pd.Categorical(df["coast_region"], categories=coast_order_long, ordered=True)
    df["water_class"] = pd.Categorical(df["water_class"], categories=water_class_order, ordered=True)

# -----------------------------------------------------------------------------
# CREATE FIGURE WITH 2x2 SUBPLOTS
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 12), constrained_layout=True)
ax_a, ax_b, ax_c, ax_d = axes.flatten()

# -----------------------------------------------------------------------------
# PANEL A – SPEI-3 sensitivity ranking
# -----------------------------------------------------------------------------
spei3_mag = resp_mag[resp_mag["SPEI_timescale"] == "SPEI_3"].copy()
spei3_summary = spei3_mag.groupby(["coast_region", "water_class"], observed=True).agg(
    mean_peak=("max_abs_pct_change", "mean"),
    sd_peak=("max_abs_pct_change", "std")
).reset_index()
# Short labels – convert to string to avoid TypeError
spei3_summary["coast_short"] = spei3_summary["coast_region"].map(coast_mapping)
spei3_summary["label"] = (
    spei3_summary["coast_short"].astype(str) + "\n" + spei3_summary["water_class"].astype(str)
)
spei3_summary = spei3_summary.sort_values("mean_peak", ascending=True)

y_pos = np.arange(len(spei3_summary))
colors = [ecosystem_palette[w] for w in spei3_summary["water_class"]]
ax_a.barh(y_pos, spei3_summary["mean_peak"], xerr=spei3_summary["sd_peak"],
          color=colors, edgecolor="black", capsize=3,
          error_kw={"elinewidth": 1.5, "ecolor": "black"})
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(spei3_summary["label"], fontsize=8)
ax_a.set_xlabel("Average monthly peak |% change|", fontsize=9)
ax_a.set_title("A) SPEI-3 sensitivity ranking", loc="left", fontsize=10)
handles = [Patch(facecolor=ecosystem_palette[w], edgecolor="black", label=w)
           for w in water_class_order]
ax_a.legend(handles=handles, loc="lower right", fontsize=8)

# -----------------------------------------------------------------------------
# PANEL B – Gulf Coast SPEI-3 signed response curves
# -----------------------------------------------------------------------------
gulf_spei3 = pred_curves[
    (pred_curves["coast_region"] == "Gulf Coast") &
    (pred_curves["SPEI_timescale"] == "SPEI_3")
].copy()
gulf_agg = gulf_spei3.groupby(["water_class", "SPEI_value"], observed=True).agg(
    mean_pct=("predicted_pct_change", "mean"),
    sd_pct=("predicted_pct_change", "std")
).reset_index()

for wc in water_class_order:
    subset = gulf_agg[gulf_agg["water_class"] == wc]
    ax_b.plot(subset["SPEI_value"], subset["mean_pct"],
              color=ecosystem_palette[wc], label=wc, linewidth=2)
    ax_b.fill_between(subset["SPEI_value"],
                      subset["mean_pct"] - subset["sd_pct"],
                      subset["mean_pct"] + subset["sd_pct"],
                      color=ecosystem_palette[wc], alpha=0.2)

ax_b.axhline(0, color="black", linestyle="--", linewidth=1)
ax_b.axvline(-1, color="grey", linestyle=":", linewidth=1)
ax_b.axvline(1, color="grey", linestyle=":", linewidth=1)
for yval in [-20, -10, -5, 5, 10, 20]:
    ax_b.axhline(yval, color="lightgrey", linestyle="-", linewidth=0.5, alpha=0.5)
ax_b.set_xlabel("SPEI", fontsize=9)
ax_b.set_ylabel("Predicted WUE_T change from near-normal (%)", fontsize=9)
ax_b.set_title("B) Gulf Coast SPEI-3 response across months", loc="left", fontsize=10)
ax_b.legend(title="Reference ecosystem", fontsize=8)

# -----------------------------------------------------------------------------
# PANEL C – Full-grid magnitude heatmap
# -----------------------------------------------------------------------------
mag_heat = resp_mag.groupby(["coast_region", "water_class", "SPEI_timescale"], observed=True).agg(
    mean_peak=("max_abs_pct_change", "mean")
).reset_index()
mag_heat["coast_short"] = mag_heat["coast_region"].map(coast_mapping)
mag_heat["row_label"] = (
    mag_heat["coast_short"].astype(str) + "\n" + mag_heat["water_class"].astype(str)
)
pivot_c = mag_heat.pivot(index="row_label", columns="SPEI_timescale", values="mean_peak")
pivot_c = pivot_c.reindex(columns=spei_order)
row_order = [coast_mapping[c] + "\n" + wc for c in coast_order_long for wc in water_class_order]
pivot_c = pivot_c.reindex(row_order)

sns.heatmap(pivot_c, annot=False, cmap="viridis",
            cbar_kws={"label": "Average monthly peak |% change|"},
            ax=ax_c, linewidths=0.5, linecolor="white")
max_val = np.nanmax(pivot_c.values)
for i, row in enumerate(pivot_c.index):
    for j, col in enumerate(pivot_c.columns):
        val = pivot_c.loc[row, col]
        if pd.notna(val):
            text_color = "white" if val < 0.55 * max_val else "black"
            label = "<0.1" if val < 0.1 else f"{val:.1f}"
            ax_c.text(j + 0.5, i + 0.5, label, ha="center", va="center",
                      color=text_color, fontsize=7)
ax_c.set_title("C) Full-grid sensitivity magnitude", loc="left", fontsize=10)
ax_c.set_xlabel("SPEI timescale", fontsize=9)
ax_c.set_ylabel("Coast / ecosystem", fontsize=9)
ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=45, ha="right", fontsize=7)
ax_c.set_yticklabels(ax_c.get_yticklabels(), fontsize=7)

# -----------------------------------------------------------------------------
# PANEL D – 10% dry-decrease threshold location (no stability filter)
# -----------------------------------------------------------------------------
thresh_10 = thresh_all[
    (thresh_all["threshold_pct"] == 10) &
    (thresh_all["anomaly_side"] == "dry") &
    (thresh_all["impact_direction"] == "decrease")
].copy()
thresh_summary = thresh_10.groupby(
    ["coast_region", "water_class", "SPEI_timescale"], observed=True
).agg(
    n_months_with=("SPEI_threshold", lambda x: x.notna().sum()),
    mean_threshold=("SPEI_threshold", "mean")
).reset_index()
thresh_summary["coast_short"] = thresh_summary["coast_region"].map(coast_mapping)
thresh_summary["row_label"] = (
    thresh_summary["coast_short"].astype(str) + "\n" + thresh_summary["water_class"].astype(str)
)

# Remove the >=6 filter: plot if detected in at least 1 month
thresh_summary["mean_threshold_plot"] = np.where(
    thresh_summary["n_months_with"] >= 1,
    thresh_summary["mean_threshold"],
    np.nan
)

pivot_d = thresh_summary.pivot(index="row_label", columns="SPEI_timescale",
                               values="mean_threshold_plot")
pivot_d = pivot_d.reindex(columns=spei_order).reindex(row_order)

sns.heatmap(pivot_d, annot=False, cmap="coolwarm",
            vmin=-3, vmax=0, cbar_kws={"label": "Mean SPEI threshold"},
            ax=ax_d, linewidths=0.5, linecolor="white")
for i, row in enumerate(pivot_d.index):
    for j, col in enumerate(pivot_d.columns):
        val = pivot_d.loc[row, col]
        if pd.notna(val):
            text_color = "white" if val <= -1.8 else "black"
            ax_d.text(j + 0.5, i + 0.5, f"{val:.2f}", ha="center", va="center",
                      color=text_color, fontsize=7)
ax_d.set_title("D) 10% dry-decrease SPEI threshold", loc="left", fontsize=10)
ax_d.set_xlabel("SPEI timescale", fontsize=9)
ax_d.set_ylabel("Coast / ecosystem", fontsize=9)
ax_d.set_xticklabels(ax_d.get_xticklabels(), rotation=45, ha="right", fontsize=7)
ax_d.set_yticklabels(ax_d.get_yticklabels(), fontsize=7)

# -----------------------------------------------------------------------------
# SAVE THE FINAL COMBINED FIGURE (ONLY THIS FILE)
# -----------------------------------------------------------------------------
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Four-panel figure saved to: {output_file}")
print("Input files used:")
print("  - prediction_curves_all_reference_conditions.csv")
print("  - response_magnitude_all_reference_conditions.csv")
print("  - threshold_summary_all_reference_conditions.csv")
print("Panel D now shows the mean 10% dry-decrease SPEI threshold across months where a threshold was detected. Blank cells indicate no detection in any month.")
print("Figure created successfully.")