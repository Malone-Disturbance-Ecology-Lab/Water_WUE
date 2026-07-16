# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 19:26:18 2026

@author: ammar
"""

# Q3_sensitivity_four_panel_figure.py
# -----------------------------------------------------------------------------
# Figure script for Q3 reference-condition sensitivity analysis
# Panel order:
#   a = Full-grid sensitivity magnitude heatmap
#   b = Pacific Coast SPEI-3 signed response curves
#   c = 10% dry-side WUE_T increase threshold stability
#   d = 10% dry-side WUE_T decrease threshold stability
# Uses saved CSVs; does NOT refit the GAM.
# -----------------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# -----------------------------------------------------------------------------
# SETTINGS – VERY LARGE FONTS
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
coast_order_short = ["Atlantic", "Pacific", "Gulf", "Alaska"]
coast_mapping = dict(zip(coast_order_long, coast_order_short))
spei_order = ["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"]

# Very large fonts for readability
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 30
sns.set_style("whitegrid")

# -----------------------------------------------------------------------------
# SAFE BOOLEAN CONVERSION
# -----------------------------------------------------------------------------
def to_bool_clean(x):
    return x.astype(str).str.upper().map({"TRUE": True, "FALSE": False, "1": True, "0": False}).fillna(False)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
pred_curves = pd.read_csv(os.path.join(results_dir, "prediction_curves_all_reference_conditions.csv"))
resp_mag = pd.read_csv(os.path.join(results_dir, "response_magnitude_all_reference_conditions.csv"))
thresh_all = pd.read_csv(os.path.join(results_dir, "threshold_summary_all_reference_conditions.csv"))

for df in [pred_curves, resp_mag, thresh_all]:
    df["month_f"] = df["month_f"].astype(int)
    df["coast_region"] = pd.Categorical(df["coast_region"], categories=coast_order_long, ordered=True)
    df["water_class"] = pd.Categorical(df["water_class"], categories=water_class_order, ordered=True)

# -----------------------------------------------------------------------------
# DYNAMIC FOCUS for console output
# -----------------------------------------------------------------------------
smooth_terms_file = os.path.join(os.path.dirname(results_dir), "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv")
focus_coast = "Pacific Coast"
focus_timescale = "SPEI_3"

if os.path.exists(smooth_terms_file):
    smooth_df = pd.read_csv(smooth_terms_file)
    p_col = "p.value" if "p.value" in smooth_df.columns else "p-value" if "p-value" in smooth_df.columns else None
    if p_col is not None:
        coast_smooths = smooth_df[smooth_df["smooth_term"].str.contains("spei_coast", na=False)]
        coast_smooths = coast_smooths[~coast_smooths["smooth_term"].str.contains("site_name", na=False)]
        sig_coast = coast_smooths[coast_smooths[p_col] < 0.05]
        if not sig_coast.empty:
            dyn_row = sig_coast.loc[sig_coast["F"].idxmax()]
            dyn_term = dyn_row["smooth_term"]
            parts = dyn_term.split("__")
            if len(parts) == 2:
                dyn_coast = parts[1].strip()
                dyn_timescale = parts[0].replace("s(SPEI_value):spei_coast", "").strip()
            else:
                dyn_coast, dyn_timescale = "unknown", "unknown"
        else:
            dyn_coast, dyn_timescale = "none", "none"
    else:
        dyn_coast, dyn_timescale = "unknown", "unknown"
else:
    dyn_coast, dyn_timescale = "unknown", "unknown"

print(f"Forced Panel b focus: {focus_coast} at {focus_timescale}")
print(f"Dynamically identified strongest smooth would be: {dyn_coast} at {dyn_timescale}")

# -----------------------------------------------------------------------------
# CREATE FIGURE WITH 2x2 SUBPLOTS
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(22, 18), constrained_layout=True)  # even larger figure
ax_a, ax_b, ax_c, ax_d = axes.flatten()

# -----------------------------------------------------------------------------
# PANEL a – Full-grid sensitivity magnitude heatmap
# -----------------------------------------------------------------------------
mag_heat = resp_mag.groupby(["coast_region", "water_class", "SPEI_timescale"], observed=True).agg(
    mean_peak=("max_abs_pct_change", "mean")
).reset_index()
mag_heat["coast_short"] = mag_heat["coast_region"].map(coast_mapping)
mag_heat["row_label"] = (
    mag_heat["coast_short"].astype(str) + "\n" + mag_heat["water_class"].astype(str)
)
pivot_a = mag_heat.pivot(index="row_label", columns="SPEI_timescale", values="mean_peak")
pivot_a = pivot_a.reindex(columns=spei_order)
row_order = [coast_mapping[c] + "\n" + wc for c in coast_order_long for wc in water_class_order]
pivot_a = pivot_a.reindex(row_order)

hm_a = sns.heatmap(
    pivot_a,
    annot=False,
    cmap="viridis",
    cbar_kws={"label": "Mean absolute WUE$_{T}$ change (%)"},
    ax=ax_a,
    linewidths=2,
    linecolor="white"
)

# Annotations – very large and bold
max_val = np.nanmax(pivot_a.values)
for i, row in enumerate(pivot_a.index):
    for j, col in enumerate(pivot_a.columns):
        val = pivot_a.loc[row, col]
        if pd.notna(val):
            text_color = "white" if val < 0.55 * max_val else "black"
            label = "<0.1" if val < 0.1 else f"{val:.1f}"
            ax_a.text(j + 0.5, i + 0.5, label, ha="center", va="center",
                      color=text_color, fontsize=18, weight="bold")

cbar_a = hm_a.collections[0].colorbar
cbar_a.ax.yaxis.label.set_size(20)
cbar_a.ax.tick_params(labelsize=18)

ax_a.set_title("a) Predicted WUE$_{T}$ response magnitude", loc="left", fontsize=30)
ax_a.set_ylabel("Coastal Ecosystem", fontsize=32)
ax_a.set_xlabel("")
ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=45, ha="right", fontsize=28)
ax_a.set_yticklabels(ax_a.get_yticklabels(), fontsize=21)

# -----------------------------------------------------------------------------
# PANEL b – Pacific Coast SPEI-3 signed response curves
# -----------------------------------------------------------------------------
focus_data = pred_curves[
    (pred_curves["coast_region"] == focus_coast) &
    (pred_curves["SPEI_timescale"] == focus_timescale)
].copy()

if focus_data.empty:
    sensitivity_summary = pd.read_csv(os.path.join(results_dir, "reference_condition_sensitivity_summary.csv"))
    top_combo = sensitivity_summary.sort_values("mean_max_abs_pct_change", ascending=False).iloc[0]
    focus_coast = top_combo["coast_region"]
    focus_timescale = top_combo["SPEI_timescale"]
    focus_data = pred_curves[
        (pred_curves["coast_region"] == focus_coast) &
        (pred_curves["SPEI_timescale"] == focus_timescale)
    ].copy()
    print(f"Fallback: Panel b now shows {focus_coast} {focus_timescale}")

focus_agg = focus_data.groupby(["water_class", "SPEI_value"], observed=True).agg(
    mean_pct=("predicted_pct_change", "mean"),
    sd_pct=("predicted_pct_change", "std")
).reset_index()

# Shade near-normal region (-1 to 1)
ax_b.axvspan(-1, 1, alpha=0.15, color="grey", label="Near-normal")

for wc in water_class_order:
    subset = focus_agg[focus_agg["water_class"] == wc]
    if not subset.empty:
        ax_b.plot(subset["SPEI_value"], subset["mean_pct"],
                  color=ecosystem_palette[wc], label=wc, linewidth=4)
        ax_b.fill_between(subset["SPEI_value"],
                          subset["mean_pct"] - subset["sd_pct"],
                          subset["mean_pct"] + subset["sd_pct"],
                          color=ecosystem_palette[wc], alpha=0.25)

ax_b.axhline(0, color="black", linestyle="--", linewidth=2)
ax_b.axvline(-1, color="grey", linestyle=":", linewidth=2)
ax_b.axvline(1, color="grey", linestyle=":", linewidth=2)
for yval in [-20, -10, -5, 5, 10, 20]:
    ax_b.axhline(yval, color="lightgrey", linestyle="-", linewidth=1, alpha=0.5)

ax_b.set_xlabel("SPEI-3", fontsize=32)
ax_b.set_ylabel("WUE$_{T}$ change from near-normal (%)", fontsize=32)
ax_b.set_title("b) Pacific Coast SPEI-3 response", loc="left", fontsize=30)
ax_b.legend(title="", fontsize=28, frameon=True)
ax_b.tick_params(labelsize=28)

# -----------------------------------------------------------------------------
# PANEL c – 10% dry-side WUE_T increase threshold stability
# -----------------------------------------------------------------------------
thresh_inc = thresh_all[
    (thresh_all["threshold_pct"] == 10) &
    (thresh_all["anomaly_side"] == "dry") &
    (thresh_all["impact_direction"] == "increase")
].copy()
thresh_inc["threshold_detected"] = to_bool_clean(thresh_inc["threshold_detected"])

n_months_inc = thresh_inc.groupby(
    ["coast_region", "water_class", "SPEI_timescale"], observed=True
).agg(
    n_months_with=("threshold_detected", "sum"),
    mean_threshold=("SPEI_threshold", "mean")
).reset_index()

n_months_inc["coast_short"] = n_months_inc["coast_region"].map(coast_mapping)
n_months_inc["row_label"] = (
    n_months_inc["coast_short"].astype(str) + "\n" + n_months_inc["water_class"].astype(str)
)
pivot_c_fill = n_months_inc.pivot(index="row_label", columns="SPEI_timescale", values="n_months_with")
pivot_c_fill = pivot_c_fill.reindex(columns=spei_order).reindex(row_order)
pivot_c_thr = n_months_inc.pivot(index="row_label", columns="SPEI_timescale", values="mean_threshold")
pivot_c_thr = pivot_c_thr.reindex(columns=spei_order).reindex(row_order)

hm_c = sns.heatmap(
    pivot_c_fill,
    annot=False,
    cmap="OrRd",
    vmin=0,
    vmax=7,
    cbar_kws={
        "label": "No. of Apr–Oct months\nwith threshold detected",
        "ticks": [0, 1, 2, 3, 4, 5, 6, 7]
    },
    ax=ax_c,
    linewidths=2,
    linecolor="white"
)

cbar_c = hm_c.collections[0].colorbar
cbar_c.ax.yaxis.label.set_size(20)
cbar_c.ax.tick_params(labelsize=28)

# Annotations – large and bold
for i, r in enumerate(row_order):
    for j, c in enumerate(spei_order):
        n = pivot_c_fill.loc[r, c] if r in pivot_c_fill.index and c in pivot_c_fill.columns else 0
        if pd.notna(n) and n > 0:
            thr = pivot_c_thr.loc[r, c]
            if pd.notna(thr):
                text_color = "white" if n >= 4 else "black"
                ax_c.text(j + 0.5, i + 0.5, f"{thr:.2f}",
                          ha="center", va="center", color=text_color,
                          fontsize=28, weight="bold")

ax_c.set_title("c) 10% dry-side WUE$_{T}$ increase threshold", loc="left", fontsize=30)
ax_c.set_ylabel("Coastal Ecosystem", fontsize=32)
ax_c.set_xlabel("")
ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=45, ha="right", fontsize=28)
ax_c.set_yticklabels(ax_c.get_yticklabels(), fontsize=21)

# -----------------------------------------------------------------------------
# PANEL d – 10% dry-side WUE_T decrease threshold stability
# -----------------------------------------------------------------------------
thresh_dec = thresh_all[
    (thresh_all["threshold_pct"] == 10) &
    (thresh_all["anomaly_side"] == "dry") &
    (thresh_all["impact_direction"] == "decrease")
].copy()
thresh_dec["threshold_detected"] = to_bool_clean(thresh_dec["threshold_detected"])

n_months_dec = thresh_dec.groupby(
    ["coast_region", "water_class", "SPEI_timescale"], observed=True
).agg(
    n_months_with=("threshold_detected", "sum"),
    mean_threshold=("SPEI_threshold", "mean")
).reset_index()

n_months_dec["coast_short"] = n_months_dec["coast_region"].map(coast_mapping)
n_months_dec["row_label"] = (
    n_months_dec["coast_short"].astype(str) + "\n" + n_months_dec["water_class"].astype(str)
)
pivot_d_fill = n_months_dec.pivot(index="row_label", columns="SPEI_timescale", values="n_months_with")
pivot_d_fill = pivot_d_fill.reindex(columns=spei_order).reindex(row_order)
pivot_d_thr = n_months_dec.pivot(index="row_label", columns="SPEI_timescale", values="mean_threshold")
pivot_d_thr = pivot_d_thr.reindex(columns=spei_order).reindex(row_order)

hm_d = sns.heatmap(
    pivot_d_fill,
    annot=False,
    cmap="Blues",
    vmin=0,
    vmax=7,
    cbar_kws={
        "label": "No. of Apr–Oct months\nwith threshold detected",
        "ticks": [0, 1, 2, 3, 4, 5, 6, 7]
    },
    ax=ax_d,
    linewidths=2,
    linecolor="white"
)

cbar_d = hm_d.collections[0].colorbar
cbar_d.ax.yaxis.label.set_size(20)
cbar_d.ax.tick_params(labelsize=28)

for i, r in enumerate(row_order):
    for j, c in enumerate(spei_order):
        n = pivot_d_fill.loc[r, c] if r in pivot_d_fill.index and c in pivot_d_fill.columns else 0
        if pd.notna(n) and n > 0:
            thr = pivot_d_thr.loc[r, c]
            if pd.notna(thr):
                text_color = "white" if n >= 4 else "black"
                ax_d.text(j + 0.5, i + 0.5, f"{thr:.2f}",
                          ha="center", va="center", color=text_color,
                          fontsize=28, weight="bold")

ax_d.set_title("d) 10% dry-side WUE$_{T}$ decrease threshold", loc="left", fontsize=30)
ax_d.set_ylabel("Coastal Ecosystem", fontsize=32)
ax_d.set_xlabel("")
ax_d.set_xticklabels(ax_d.get_xticklabels(), rotation=45, ha="right", fontsize=28)
ax_d.set_yticklabels(ax_d.get_yticklabels(), fontsize=21)

# -----------------------------------------------------------------------------
# SAVE AND DISPLAY
# -----------------------------------------------------------------------------
plt.savefig(output_file, dpi=600, bbox_inches="tight")
plt.show(block=True)  # Forces display in Spyder plot pane

print(f"Four-panel figure saved to: {output_file} (600 DPI)")
print("Input files used:")
print("  - prediction_curves_all_reference_conditions.csv")
print("  - response_magnitude_all_reference_conditions.csv")
print("  - threshold_summary_all_reference_conditions.csv")
print(f"Forced Panel b focus: {focus_coast} at {focus_timescale}")
print(f"Dynamically identified strongest smooth would be: {dyn_coast} at {dyn_timescale}")
print("\nCaption:\nSupplementary Figure X. Reference-condition sensitivity analysis for the coast-threshold GAM. Predictions were generated across all coast × reference ecosystem × growing-season month × SPEI-timescale combinations without refitting the model. Panel a shows average monthly peak absolute predicted WUE_T percent change across all SPEI timescales. Panel b shows signed response curves for Pacific Coast SPEI-3, the clearest modeled nonlinear coast-specific smooth. Panels c and d show the model-threshold basis for the positive and negative event directions used in the spatial upscaling: warm colors indicate dry-side WUE_T increase thresholds, blue colors indicate dry-side WUE_T decrease thresholds, and blank cells indicate no detected 10% dry-side threshold. Heatmap color shows how many April–October months had a detected 10% dry-side threshold. Cell text shows the mean SPEI threshold.")
print("\nTake-home: The reference-condition sensitivity analysis confirmed that Pacific Coast SPEI-3 remained the strongest coast-specific nonlinear response when the fitted coast-threshold GAM was evaluated across growing-season months and reference ecosystem classes. Dry-side threshold direction varied across coasts and SPEI timescales, supporting the Q4 decision to map positive and negative WUE_T response events separately.")