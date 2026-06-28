# -*- coding: utf-8 -*-
"""
Standalone manuscript-ready Q4 reporting script.

Uses outputs from updated Step 2:
- Weighted GLM, class-balanced
- Standardized SPEI-48 coefficient
- Probability table from weighted GLM

This script reports:
- manuscript-ready values
- corrected interpretation text
- class-balanced probability table
"""

import pandas as pd
import numpy as np
import os

# =============================================================================
# FILE PATHS
# =============================================================================

summary_path = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Linear_Summary_10_90_Decrease_vs_Baseline.csv"
prob_path = r"M:\Research\WUE_CUE\data_products\results\linear_models\WUE_tra_Probability_Table_10_90.csv"
output_dir = r"M:\Research\WUE_CUE\data_products\results\linear_models"

# =============================================================================
# LOAD FILES
# =============================================================================

summary_df = pd.read_csv(summary_path)
prob_df = pd.read_csv(prob_path)

print("✅ Loaded summary and probability table")
print(f"Summary shape: {summary_df.shape}")
print(f"Probability table shape: {prob_df.shape}")

# =============================================================================
# BASIC COLUMN CHECKS
# =============================================================================

required_summary_cols = [
    "N_Total",
    "N_Decrease",
    "N_Increase",
    "N_NoChange",
    "Pct_Decrease",
    "Coefficient",
    "Coefficient_CI_lower",
    "Coefficient_CI_upper",
    "Coefficient_P_value",
    "Odds_Ratio_1unit_decrease",
    "ROC_AUC",
    "PR_AUC",
    "Baseline_PR",
    "Inference_Method",
    "spei_mean",
    "spei_std"
]

required_prob_cols = [
    "SPEI",
    "P_Decrease_SPEI48"
]

missing_summary_cols = [col for col in required_summary_cols if col not in summary_df.columns]
missing_prob_cols = [col for col in required_prob_cols if col not in prob_df.columns]

if missing_summary_cols:
    raise ValueError(f"❌ Missing required summary columns: {missing_summary_cols}")

if missing_prob_cols:
    raise ValueError(f"❌ Missing required probability table columns: {missing_prob_cols}")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_probability_at_spei(prob_df, spei_value, prob_col="P_Decrease_SPEI48"):
    """
    Return class-balanced fitted probability at the closest SPEI value
    in the probability table.
    """
    idx = (prob_df["SPEI"] - spei_value).abs().idxmin()
    return float(prob_df.loc[idx, prob_col])


def print_manuscript_values(summary_df, prob_df):
    """
    Print manuscript-ready values and corrected manuscript sentence.
    """

    s = summary_df.iloc[0]

    n_total = int(s["N_Total"])
    n_decrease = int(s["N_Decrease"])
    n_increase = int(s["N_Increase"])
    n_nochange = int(s["N_NoChange"])
    pct_decrease = float(s["Pct_Decrease"])

    # This is the standardized coefficient from Step 2
    coef_z = float(s["Coefficient"])
    ci_lower = float(s["Coefficient_CI_lower"])
    ci_upper = float(s["Coefficient_CI_upper"])
    p_value = float(s["Coefficient_P_value"])

    odds_ratio_1unit = float(s["Odds_Ratio_1unit_decrease"])
    odds_pct = (odds_ratio_1unit - 1) * 100

    roc_auc = float(s["ROC_AUC"])
    pr_auc = float(s["PR_AUC"])
    baseline_pr = float(s["Baseline_PR"])

    inference_method = str(s["Inference_Method"])
    spei_mean = float(s["spei_mean"])
    spei_std = float(s["spei_std"])

    # Class-balanced fitted probabilities from weighted GLM table
    p_dry = get_probability_at_spei(prob_df, -3.0)
    p_wet = get_probability_at_spei(prob_df, 1.0)
    p_zero = get_probability_at_spei(prob_df, 0.0)
    p_nn_dry = get_probability_at_spei(prob_df, -0.5)

    abs_change = p_dry - p_wet

    print("\n" + "=" * 70)
    print("MANUSCRIPT-READY Q4 VALUES")
    print("=" * 70)

    print(f"N total: {n_total:,}")
    print(f"N decrease: {n_decrease:,} ({pct_decrease:.1f}%)")
    print(f"N increase: {n_increase:,}")
    print(f"N no change: {n_nochange:,}")

    print(f"\nModel type: weighted GLM, class-balanced")
    print(f"Inference method: {inference_method}")
    print(f"SPEI-48 standardization: mean = {spei_mean:.3f}, SD = {spei_std:.3f}")

    print(f"\nStandardized beta, βz: {coef_z:.3f}")
    print(f"95% CI for βz: [{ci_lower:.3f}, {ci_upper:.3f}]")
    print(f"p-value: {p_value:.3f}")

    print(f"\nOdds ratio per 1-unit SPEI decrease: {odds_ratio_1unit:.3f}")
    print(f"Percent increase in odds per 1-unit SPEI decrease: {odds_pct:.1f}%")

    print(f"\nROC-AUC: {roc_auc:.2f}")
    print(f"PR-AUC: {pr_auc:.2f}")
    print(f"Baseline PR: {baseline_pr:.2f}")

    print(f"\nClass-balanced fitted P decrease at SPEI-48 = 1.0: {p_wet:.2f}")
    print(f"Class-balanced fitted P decrease at SPEI-48 = -3.0: {p_dry:.2f}")
    print(f"Absolute fitted-probability increase from 1.0 to -3.0: {abs_change * 100:.1f} percentage points")
    print(f"Class-balanced fitted P decrease at SPEI-48 = 0.0: {p_zero:.2f}")
    print(f"Class-balanced fitted P decrease at SPEI-48 = -0.5: {p_nn_dry:.2f}")

    print("\n" + "=" * 70)
    print("MANUSCRIPT SENTENCE")
    print("=" * 70)

    print(
        f"The final class-balanced logistic model included {n_total:,} observations, "
        f"with WUE_T decreases representing {pct_decrease:.1f}% of the modeled observations. "
        f"SPEI-48 was negatively associated with the probability of WUE_T decrease "
        f"(standardized βz = {coef_z:.3f} per 1 SD increase in SPEI-48, "
        f"95% CI: {ci_lower:.3f} to {ci_upper:.3f}, p = {p_value:.3f}). "
        f"Each one-unit decrease in SPEI-48 increased the odds of WUE_T decline by "
        f"approximately {odds_pct:.1f}%. The class-balanced fitted probability increased "
        f"from {p_wet:.2f} at SPEI-48 = 1.0 to {p_dry:.2f} at SPEI-48 = -3.0, "
        f"an absolute increase of {abs_change * 100:.1f} percentage points. "
        f"Site-blocked cross-validation yielded ROC-AUC = {roc_auc:.2f} and "
        f"PR-AUC = {pr_auc:.2f}, compared with a baseline decline prevalence of "
        f"{baseline_pr:.2f}."
    )


def create_probability_table_for_manuscript(prob_df, output_dir):
    """
    Create a small manuscript/supplement table showing class-balanced fitted
    WUE_T decline probabilities at selected SPEI-48 values.
    """

    selected_spei = [1.0, 0.0, -1.0, -2.0, -3.0]

    interpretation = {
        1.0: "Wetter multi-year conditions",
        0.0: "Near-neutral multi-year moisture",
        -1.0: "Moderate long-term drying",
        -2.0: "Strong long-term drying",
        -3.0: "Severe long-term drought"
    }

    table_rows = []

    for spei_value in selected_spei:
        prob = get_probability_at_spei(prob_df, spei_value)

        table_rows.append({
            "SPEI-48": spei_value,
            "Class-balanced fitted P(WUE_T decrease)": round(prob, 2),
            "Class-balanced fitted P(WUE_T decrease, %)": round(prob * 100, 1),
            "Interpretation": interpretation[spei_value]
        })

    manuscript_table = pd.DataFrame(table_rows)

    output_path = os.path.join(
        output_dir,
        "WUE_tra_SPEI48_probability_table_for_manuscript.csv"
    )

    manuscript_table.to_csv(output_path, index=False)

    print("\n" + "=" * 70)
    print("MANUSCRIPT PROBABILITY TABLE")
    print("=" * 70)
    print(manuscript_table.to_string(index=False))
    print(f"\n✅ Saved table: {output_path}")

    return manuscript_table


# =============================================================================
# RUN
# =============================================================================

print_manuscript_values(summary_df, prob_df)

manuscript_table = create_probability_table_for_manuscript(
    prob_df=prob_df,
    output_dir=output_dir
)