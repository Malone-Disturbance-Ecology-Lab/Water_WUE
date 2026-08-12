import os
import pandas as pd
import numpy as np

# ============================================================
# Q1 DIAGNOSTIC:
# How strongly does WUE_ET - WUE_T divergence change with T:ET?
# Uses EXISTING final smooth-GAM predictions.
# SAVES NOTHING.
# ============================================================

base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")

pred_file = os.path.join(
    output_dir,
    "Q1_final_smooth_gam_predictions_TET.csv"
)

terms_file = os.path.join(
    output_dir,
    "Q1_final_smooth_gam_smooth_terms.csv"
)

pred = pd.read_csv(pred_file)
terms = pd.read_csv(terms_file)

# ------------------------------------------------------------
# Keep only the biologically important comparison:
# WUE_ET - WUE_T
# ------------------------------------------------------------
x = pred[
    pred["difference_label"] == "WUE_ET - WUE_T"
].copy()

ecosystems = ["Upland", "Freshwater", "Saline"]

# Find GAM p-value column
p_cols = [
    c for c in terms.columns
    if "p-value" in c.lower() or c.lower() in ["p", "p_value", "pvalue"]
]

p_col = p_cols[0] if len(p_cols) > 0 else None

results = []

print("\n" + "=" * 80)
print("WUE_ET - WUE_T DIVERGENCE ALONG THE T:ET GRADIENT")
print("Final smooth GAM; central 90% of observed T:ET")
print("=" * 80)

for eco in ecosystems:

    d = x[x["ecosystem_label"] == eco].sort_values("Trans_ratio").copy()

    if d.empty:
        continue

    # Low and high ends of the modeled T:ET range
    low = d.iloc[0]
    high = d.iloc[-1]

    tet_low = low["Trans_ratio"]
    tet_high = high["Trans_ratio"]

    pred_low = low["prediction"]
    pred_high = high["prediction"]

    # Signed endpoint change
    signed_change = pred_high - pred_low

    # Magnitude of divergence from zero
    abs_low = abs(pred_low)
    abs_high = abs(pred_high)
    magnitude_change = abs_high - abs_low

    # Average fitted change per +0.10 T:ET
    change_per_01 = (
        signed_change / (tet_high - tet_low)
    ) * 0.10

    magnitude_change_per_01 = (
        magnitude_change / (tet_high - tet_low)
    ) * 0.10

    # Maximum absolute WUE_ET - WUE_T divergence anywhere
    # across the modeled central 90% T:ET range
    idx_max = d["prediction"].abs().idxmax()
    max_row = d.loc[idx_max]

    max_abs_divergence = abs(max_row["prediction"])
    tet_at_max = max_row["Trans_ratio"]
    signed_at_max = max_row["prediction"]

    # Minimum absolute divergence across the curve
    min_abs_divergence = d["prediction"].abs().min()

    # Total variation in absolute divergence across T:ET
    divergence_range = max_abs_divergence - min_abs_divergence

    # Get p-value for this WUE_ET-WUE_T × ecosystem smooth
    p_val = np.nan

    if p_col is not None:
        target = terms[
            terms["term"].astype(str).str.contains(
                "WUE_ET_minus_WUE_T", regex=False
            )
            &
            terms["term"].astype(str).str.contains(
                eco, regex=False
            )
        ]

        if len(target) > 0:
            p_val = target.iloc[0][p_col]

    results.append({
        "Ecosystem": eco,
        "TET_low": tet_low,
        "TET_high": tet_high,
        "WUEdiff_low": pred_low,
        "WUEdiff_high": pred_high,
        "signed_change": signed_change,
        "magnitude_change": magnitude_change,
        "change_per_0.1_TET": change_per_01,
        "magnitude_change_per_0.1_TET": magnitude_change_per_01,
        "max_abs_divergence": max_abs_divergence,
        "TET_at_max_divergence": tet_at_max,
        "signed_value_at_max": signed_at_max,
        "divergence_range": divergence_range,
        "p_value": p_val
    })

    print(f"\n{eco}")
    print("-" * 50)

    print(
        f"T:ET modeled range: "
        f"{tet_low:.3f} to {tet_high:.3f}"
    )

    print(
        f"WUE_ET - WUE_T at low T:ET: "
        f"{pred_low:.3f}"
    )

    print(
        f"WUE_ET - WUE_T at high T:ET: "
        f"{pred_high:.3f}"
    )

    print(
        f"Signed change across T:ET gradient: "
        f"{signed_change:+.3f}"
    )

    print(
        f"Change in absolute divergence: "
        f"{magnitude_change:+.3f}"
    )

    print(
        f"Average fitted signed change per +0.10 T:ET: "
        f"{change_per_01:+.3f}"
    )

    print(
        f"Average change in divergence magnitude per +0.10 T:ET: "
        f"{magnitude_change_per_01:+.3f}"
    )

    print(
        f"Maximum |WUE_ET - WUE_T| divergence: "
        f"{max_abs_divergence:.3f}"
    )

    print(
        f"T:ET where maximum divergence occurs: "
        f"{tet_at_max:.3f}"
    )

    if not np.isnan(p_val):
        if p_val < 0.001:
            print("GAM T:ET relationship: p < 0.001")
        else:
            print(f"GAM T:ET relationship: p = {p_val:.4f}")


# ------------------------------------------------------------
# Compare ecosystems
# ------------------------------------------------------------

res = pd.DataFrame(results)

print("\n" + "=" * 80)
print("MAIN COMPARISONS")
print("=" * 80)

largest_absolute = res.loc[
    res["max_abs_divergence"].idxmax()
]

largest_gradient_change = res.loc[
    res["divergence_range"].idxmax()
]

largest_endpoint_change = res.loc[
    res["magnitude_change"].abs().idxmax()
]

print(
    "\nGreatest absolute WUE_ET - WUE_T divergence:"
)
print(
    f"{largest_absolute['Ecosystem']} "
    f"(|difference| = {largest_absolute['max_abs_divergence']:.3f}, "
    f"at T:ET = {largest_absolute['TET_at_max_divergence']:.3f})"
)

print(
    "\nGreatest variation in divergence across the T:ET gradient:"
)
print(
    f"{largest_gradient_change['Ecosystem']} "
    f"(range in |difference| = "
    f"{largest_gradient_change['divergence_range']:.3f})"
)

print(
    "\nLargest low-to-high change in divergence magnitude:"
)
print(
    f"{largest_endpoint_change['Ecosystem']} "
    f"(change = "
    f"{largest_endpoint_change['magnitude_change']:+.3f})"
)

print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)

cols = [
    "Ecosystem",
    "WUEdiff_low",
    "WUEdiff_high",
    "magnitude_change",
    "magnitude_change_per_0.1_TET",
    "max_abs_divergence",
    "TET_at_max_divergence",
    "p_value"
]

print(
    res[cols].round(3).to_string(index=False)
)