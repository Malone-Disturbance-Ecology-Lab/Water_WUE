# -*- coding: utf-8 -*-
"""
Q2: mean T:ET relationship with WUE_T performance metrics

Saves only:
  1) Q2_T_ET_performance_relationship_summary.csv
  2) Q2_T_ET_performance_relationship_long_for_plot.csv

No TXT files.
No figures.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

input_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q2_site_performance_summary.csv"

output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship"
os.makedirs(output_dir, exist_ok=True)

summary_csv = os.path.join(
    output_dir,
    "Q2_T_ET_performance_relationship_summary.csv"
)

plot_data_csv = os.path.join(
    output_dir,
    "Q2_T_ET_performance_relationship_long_for_plot.csv"
)

# Optional: remove old TXT file from earlier version if it exists
old_txt = os.path.join(
    output_dir,
    "Q2_T_ET_performance_relationship_summary.txt"
)
if os.path.exists(old_txt):
    os.remove(old_txt)
    print(f"Removed old TXT file: {old_txt}")

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

x_col = "mean_TET"

metrics_all = [
    ("stability", "Stability index", "A", "No SPEI used"),
    ("plasticity_slope_SPEI3", "Plasticity slope (SPEI-3)", "B", "Representative site-level plasticity slope using SPEI-3"),
    ("plasticity_p95_p05", "Plasticity range (95th/5th)", "C", "No SPEI used"),
    ("resistance", "Drought resistance", "D", "Drought defined using SPEI-3"),
    ("mean_recovery", "Drought recovery", "E", "Drought events defined using SPEI-3"),
]

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)

    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = np.empty(n, dtype=float)

    running_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        running_min = min(running_min, value)
        adjusted[i] = running_min

    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def sig_label(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "."
    return "ns"


def interpretation_text(rho, p_fdr):
    if pd.isna(rho) or pd.isna(p_fdr):
        return "not tested"

    if p_fdr < 0.05:
        if rho > 0:
            return "significant positive relationship"
        elif rho < 0:
            return "significant negative relationship"
        else:
            return "significant but near-zero relationship"

    if p_fdr < 0.10:
        if rho > 0:
            return "marginal positive relationship"
        elif rho < 0:
            return "marginal negative relationship"
        else:
            return "marginal but near-zero relationship"

    return "not significant"

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

df = pd.read_csv(input_file)

print("=" * 90)
print("Q2: mean T:ET relationship with WUE_T performance metrics")
print("=" * 90)
print(f"Input file: {input_file}")
print(f"Rows: {len(df)}")
print(f"Sites: {df['site_name'].nunique()}")

# -----------------------------------------------------------------------------
# Save long-format site-level plotting data
# -----------------------------------------------------------------------------

plot_rows = []

for metric, label, panel_label, notes in metrics_all:
    if metric not in df.columns:
        continue

    sub = df[
        ["site_name", "water_class", x_col, metric]
    ].copy()

    sub = sub.rename(columns={metric: "metric_value"})
    sub["metric"] = metric
    sub["label"] = label
    sub["panel_label_original_Q2"] = panel_label
    sub["notes"] = notes

    sub = sub.replace([np.inf, -np.inf], np.nan)
    sub = sub.dropna(subset=[x_col, "metric_value"])

    plot_rows.append(sub)

plot_df = pd.concat(plot_rows, ignore_index=True)
plot_df.to_csv(plot_data_csv, index=False)

# -----------------------------------------------------------------------------
# Spearman tests for all five performance metrics
# -----------------------------------------------------------------------------

results = []

for metric, label, panel_label, notes in metrics_all:
    if metric not in df.columns:
        results.append({
            "metric": metric,
            "label": label,
            "panel_label_original_Q2": panel_label,
            "notes": notes,
            "n_sites": 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "mean_TET_min": np.nan,
            "mean_TET_max": np.nan,
            "metric_min": np.nan,
            "metric_max": np.nan
        })
        continue

    sub = df[[x_col, metric, "site_name", "water_class"]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan)
    sub = sub.dropna(subset=[x_col, metric])

    n = len(sub)

    if n < 3:
        results.append({
            "metric": metric,
            "label": label,
            "panel_label_original_Q2": panel_label,
            "notes": notes,
            "n_sites": n,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "mean_TET_min": sub[x_col].min() if n > 0 else np.nan,
            "mean_TET_max": sub[x_col].max() if n > 0 else np.nan,
            "metric_min": sub[metric].min() if n > 0 else np.nan,
            "metric_max": sub[metric].max() if n > 0 else np.nan
        })
        continue

    rho, p = spearmanr(sub[x_col], sub[metric])

    results.append({
        "metric": metric,
        "label": label,
        "panel_label_original_Q2": panel_label,
        "notes": notes,
        "n_sites": n,
        "spearman_rho": rho,
        "spearman_p": p,
        "mean_TET_min": sub[x_col].min(),
        "mean_TET_max": sub[x_col].max(),
        "metric_min": sub[metric].min(),
        "metric_max": sub[metric].max()
    })

results_df = pd.DataFrame(results)

# FDR correction across the five primary Spearman tests
results_df["spearman_p_FDR"] = bh_adjust(results_df["spearman_p"].values)
results_df["sig_FDR"] = results_df["spearman_p_FDR"].apply(sig_label)

results_df["interpretation"] = results_df.apply(
    lambda row: interpretation_text(row["spearman_rho"], row["spearman_p_FDR"]),
    axis=1
)

results_df.to_csv(summary_csv, index=False)

# -----------------------------------------------------------------------------
# Console output
# -----------------------------------------------------------------------------

print("\n" + "=" * 90)
print("Summary results")
print("=" * 90)

print(
    results_df[
        [
            "metric",
            "label",
            "n_sites",
            "spearman_rho",
            "spearman_p",
            "spearman_p_FDR",
            "sig_FDR",
            "interpretation",
            "notes"
        ]
    ].to_string(index=False, float_format=lambda x: f"{x:.4g}")
)

print("\nSaved CSV files:")
print(f"- {summary_csv}")
print(f"- {plot_data_csv}")

print("\nDone. Only two CSV files were saved. No TXT files or figures were created.")