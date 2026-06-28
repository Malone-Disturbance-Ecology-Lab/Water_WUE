# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:40:17 2026

@author: ammar
"""

from pathlib import Path
import numpy as np
import pandas as pd

# =============================================================================
# Q1 DIAGNOSTIC: OLD MALONE WORKFLOW VS UPDATED FINAL WORKFLOW
#
# Updated final structure:
#   A = near-normal summary
#   B = near-normal mixed model
#   C = near-normal smooth GAM
#   D = linear vs smooth GAM comparison
# =============================================================================

# ---------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------

old_dir = Path(r"M:\Research\WUE_CUE\Water_WUE\Malone_Workflow\results\Q1_WUE_metric_difference_simple")

new_dir = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs")

diagnostic_out = new_dir / "Q1_diagnostic_old_malone_vs_updated_final"
diagnostic_out.mkdir(exist_ok=True)

# ---------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------

def read_csv_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path)


def sign_label(x):
    if pd.isna(x):
        return "NA"
    if x > 0:
        return "positive"
    if x < 0:
        return "negative"
    return "zero"


def clean_model_name(x):
    x = str(x).lower()
    if "linear" in x:
        return "Linear GAM"
    if "smooth" in x:
        return "Smooth GAM"
    return str(x)


def get_p_col(df: pd.DataFrame):
    candidates = [
        c for c in df.columns
        if "Pr" in c or "p_value" in c.lower() or c.lower() in ["p", "p.value"]
    ]
    return candidates[0] if candidates else None


def extract_interaction_lrt(lrt_df: pd.DataFrame):
    """
    Extracts the three-way mixed-model interaction row.
    Works for old Trans_ratio_z and updated Trans_ratio_z_NN.
    """
    if "term" not in lrt_df.columns:
        return None

    rows = lrt_df[
        lrt_df["term"].astype(str).str.contains(
            "difference_type:water_class:Trans_ratio",
            regex=False,
            na=False
        )
    ].copy()

    if rows.empty:
        return None

    p_col = get_p_col(rows)
    out = rows.iloc[0].to_dict()
    out["_p_col"] = p_col
    return out


def compare_prediction_curves(old_pred, new_pred, old_name, new_name, out_name):
    """
    Compares prediction curves by interpolating old and new predictions
    onto the same overlapping T:ET grid for each difference x ecosystem.
    """
    rows = []

    for df in [old_pred, new_pred]:
        if "ecosystem_label" not in df.columns and "water_class" in df.columns:
            df["ecosystem_label"] = df["water_class"]

    required = ["difference_label", "ecosystem_label", "Trans_ratio", "prediction"]

    for df_name, df in [(old_name, old_pred), (new_name, new_pred)]:
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"{df_name} missing columns: {missing}")

    groups = sorted(
        set(zip(old_pred["difference_label"], old_pred["ecosystem_label"]))
        .intersection(set(zip(new_pred["difference_label"], new_pred["ecosystem_label"])))
    )

    for diff_label, eco in groups:
        old_sub = old_pred[
            (old_pred["difference_label"] == diff_label) &
            (old_pred["ecosystem_label"] == eco)
        ].sort_values("Trans_ratio")

        new_sub = new_pred[
            (new_pred["difference_label"] == diff_label) &
            (new_pred["ecosystem_label"] == eco)
        ].sort_values("Trans_ratio")

        if old_sub.empty or new_sub.empty:
            continue

        xmin = max(old_sub["Trans_ratio"].min(), new_sub["Trans_ratio"].min())
        xmax = min(old_sub["Trans_ratio"].max(), new_sub["Trans_ratio"].max())

        if xmin >= xmax:
            continue

        xgrid = np.linspace(xmin, xmax, 200)

        old_y = np.interp(xgrid, old_sub["Trans_ratio"], old_sub["prediction"])
        new_y = np.interp(xgrid, new_sub["Trans_ratio"], new_sub["prediction"])

        delta = new_y - old_y

        low_x = np.quantile(xgrid, 0.05)
        mid_x = np.quantile(xgrid, 0.50)
        high_x = np.quantile(xgrid, 0.95)

        old_low = np.interp(low_x, xgrid, old_y)
        old_mid = np.interp(mid_x, xgrid, old_y)
        old_high = np.interp(high_x, xgrid, old_y)

        new_low = np.interp(low_x, xgrid, new_y)
        new_mid = np.interp(mid_x, xgrid, new_y)
        new_high = np.interp(high_x, xgrid, new_y)

        old_slope_direction = sign_label(old_high - old_low)
        new_slope_direction = sign_label(new_high - new_low)

        rows.append({
            "comparison": f"{old_name} vs {new_name}",
            "difference_label": diff_label,
            "ecosystem": eco,
            "tet_min_overlap": xmin,
            "tet_max_overlap": xmax,
            "old_low": old_low,
            "new_low": new_low,
            "delta_low": new_low - old_low,
            "old_mid": old_mid,
            "new_mid": new_mid,
            "delta_mid": new_mid - old_mid,
            "old_high": old_high,
            "new_high": new_high,
            "delta_high": new_high - old_high,
            "old_slope_direction": old_slope_direction,
            "new_slope_direction": new_slope_direction,
            "slope_direction_changed": old_slope_direction != new_slope_direction,
            "mean_abs_delta": np.mean(np.abs(delta)),
            "max_abs_delta": np.max(np.abs(delta)),
            "same_sign_percent": np.mean(np.sign(old_y) == np.sign(new_y)) * 100,
            "old_range": old_y.max() - old_y.min(),
            "new_range": new_y.max() - new_y.min(),
            "range_delta": (new_y.max() - new_y.min()) - (old_y.max() - old_y.min()),
        })

    out = pd.DataFrame(rows)
    out.to_csv(diagnostic_out / out_name, index=False)
    return out


# ---------------------------------------------------------------------
# LOAD OLD MALONE OUTPUTS
# ---------------------------------------------------------------------

old_summary = read_csv_required(old_dir / "Q1_simple_summary_by_difference_type.csv")
old_gam_comp = read_csv_required(old_dir / "Q1_simple_gam_model_comparison.csv")
old_mixed_lrt = read_csv_required(old_dir / "Q1_simple_mixed_likelihood_ratio_tests.csv")
old_mixed_comp = read_csv_required(old_dir / "Q1_simple_mixed_model_comparison_lrt.csv")
old_mixed_pred = read_csv_required(old_dir / "Q1_simple_mixed_predictions_TET.csv")
old_smooth_pred = read_csv_required(old_dir / "Q1_simple_gam_predictions_TET.csv")

# ---------------------------------------------------------------------
# LOAD UPDATED FINAL OUTPUTS
# ---------------------------------------------------------------------

new_summary = read_csv_required(new_dir / "Q1_final_near_normal_summary_by_difference_type.csv")
new_gam_comp = read_csv_required(new_dir / "Q1_final_gam_model_comparison.csv")
new_mixed_lrt = read_csv_required(new_dir / "Q1_final_near_normal_mixed_likelihood_ratio_tests.csv")
new_mixed_comp = read_csv_required(new_dir / "Q1_final_near_normal_mixed_model_comparison_lrt.csv")
new_mixed_pred = read_csv_required(new_dir / "Q1_final_near_normal_mixed_predictions_TET.csv")
new_smooth_pred = read_csv_required(new_dir / "Q1_final_smooth_gam_predictions_TET.csv")

# ---------------------------------------------------------------------
# 1. PANEL A SUMMARY COMPARISON
# ---------------------------------------------------------------------

panel_a = old_summary.merge(
    new_summary,
    on=["difference_type", "difference_label", "water_class"],
    suffixes=("_old", "_new")
)

panel_a["n_months_change"] = panel_a["n_months_new"] - panel_a["n_months_old"]
panel_a["mean_change"] = panel_a["mean_difference_new"] - panel_a["mean_difference_old"]
panel_a["abs_mean_change"] = panel_a["mean_change"].abs()
panel_a["old_sign"] = panel_a["mean_difference_old"].apply(sign_label)
panel_a["new_sign"] = panel_a["mean_difference_new"].apply(sign_label)
panel_a["direction_changed"] = panel_a["old_sign"] != panel_a["new_sign"]
panel_a["ci95_change"] = panel_a["ci95_difference_new"] - panel_a["ci95_difference_old"]

panel_a_out = panel_a[[
    "difference_label", "water_class",
    "n_months_old", "n_months_new", "n_months_change",
    "mean_difference_old", "mean_difference_new",
    "mean_change", "abs_mean_change",
    "old_sign", "new_sign", "direction_changed",
    "ci95_difference_old", "ci95_difference_new", "ci95_change"
]].copy()

panel_a_out.to_csv(
    diagnostic_out / "Q1_diagnostic_panel_A_old_full_vs_updated_near_normal.csv",
    index=False
)

# ---------------------------------------------------------------------
# 2. GAM MODEL COMPARISON
# ---------------------------------------------------------------------

old_gam_comp["model_clean"] = old_gam_comp["model"].apply(clean_model_name)
new_gam_comp["model_clean"] = new_gam_comp["model"].apply(clean_model_name)

model_comp = old_gam_comp.merge(
    new_gam_comp,
    on="model_clean",
    suffixes=("_old", "_new")
)

model_comp["AIC_change"] = model_comp["AIC_new"] - model_comp["AIC_old"]
model_comp["BIC_change"] = model_comp["BIC_new"] - model_comp["BIC_old"]
model_comp["deviance_explained_change"] = (
    model_comp["deviance_explained_new"] -
    model_comp["deviance_explained_old"]
)
model_comp["delta_AIC_change"] = model_comp["delta_AIC_new"] - model_comp["delta_AIC_old"]

model_comp_out = model_comp[[
    "model_clean",
    "AIC_old", "AIC_new", "AIC_change",
    "BIC_old", "BIC_new", "BIC_change",
    "deviance_explained_old", "deviance_explained_new", "deviance_explained_change",
    "delta_AIC_old", "delta_AIC_new", "delta_AIC_change"
]].copy()

model_comp_out.to_csv(
    diagnostic_out / "Q1_diagnostic_g