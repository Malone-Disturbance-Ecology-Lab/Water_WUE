# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:11:31 2026

@author: ammar
"""

import os
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================

OUT = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"

smooth_file = os.path.join(
    OUT,
    "Q3_WUE_T_SPEI_smooth_terms.csv"
)

pred_file = os.path.join(
    OUT,
    "Q3_WUE_T_SPEI_smooth_predictions.csv"
)

slope_file = os.path.join(
    OUT,
    "Q3_WUE_T_SPEI_ecosystem_slope_summary.csv"
)

# ============================================================
# LOAD
# ============================================================

smooth = pd.read_csv(smooth_file)
pred = pd.read_csv(pred_file)
slopes = pd.read_csv(slope_file)

print("=" * 100)
print("Q3 ECOSYSTEM-SPECIFIC WUE_T RESPONSE DIAGNOSTIC")
print("=" * 100)

print("\nSmooth-term columns:")
print(list(smooth.columns))

# ============================================================
# PARSE SMOOTH TERMS
# ============================================================

def parse_term(term):
    """
    Example:
    s(SPEI_value):spei_ecosystemSPEI_3__Upland
    """
    text = str(term)

    if "spei_ecosystem" not in text:
        return pd.Series([np.nan, np.nan])

    part = text.split("spei_ecosystem", 1)[1]

    if "__" not in part:
        return pd.Series([np.nan, np.nan])

    timescale, ecosystem = part.split("__", 1)

    return pd.Series([timescale, ecosystem])


smooth[["timescale", "ecosystem"]] = smooth["smooth_term"].apply(
    parse_term
)

# Robust column detection
edf_col = next(
    c for c in smooth.columns
    if c.lower() == "edf"
)

f_col = next(
    c for c in smooth.columns
    if c.lower() == "f"
)

p_col = next(
    c for c in smooth.columns
    if "p" in c.lower()
)

smooth["significant"] = smooth[p_col] < 0.05

timescale_order = [
    "SPEI_1", "SPEI_3", "SPEI_6",
    "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"
]

ecosystem_order = ["Upland", "Freshwater", "Saline"]

# ============================================================
# 1. EXACT GAM SMOOTH RESULTS
# ============================================================

print("\n" + "=" * 100)
print("1. ECOSYSTEM-SPECIFIC GAM SMOOTH TERMS")
print("=" * 100)

for eco in ecosystem_order:

    print(f"\n{eco}")
    print("-" * 80)

    sub = smooth[
        smooth["ecosystem"] == eco
    ].copy()

    sub["timescale"] = pd.Categorical(
        sub["timescale"],
        categories=timescale_order,
        ordered=True
    )

    sub = sub.sort_values("timescale")

    print(
        sub[
            [
                "timescale",
                edf_col,
                f_col,
                p_col,
                "significant"
            ]
        ]
        .round(4)
        .to_string(index=False)
    )

# ============================================================
# 2. STRONGEST SMOOTH FOR EACH ECOSYSTEM
# ============================================================

print("\n" + "=" * 100)
print("2. STRONGEST GAM SMOOTH BY ECOSYSTEM")
print("=" * 100)

for eco in ecosystem_order:

    sub = smooth[
        smooth["ecosystem"] == eco
    ].copy()

    if sub.empty:
        continue

    row = sub.loc[sub[f_col].idxmax()]

    print(
        f"{eco}: {row['timescale']} | "
        f"edf={row[edf_col]:.2f}, "
        f"F={row[f_col]:.2f}, "
        f"p={row[p_col]:.4g}"
    )

# ============================================================
# 3. NUMBER OF SIGNIFICANT TIMESCALES
# ============================================================

print("\n" + "=" * 100)
print("3. NUMBER OF STATISTICALLY DETECTABLE TIMESCALES")
print("=" * 100)

for eco in ecosystem_order:

    sub = smooth[
        smooth["ecosystem"] == eco
    ]

    sig = sub[sub["significant"]]

    print(
        f"{eco}: {len(sig)} of 7 significant | "
        f"{', '.join(sig['timescale'].astype(str).tolist()) if len(sig) else 'none'}"
    )

# ============================================================
# 4. FITTED DIRECTION FROM SAVED GAM PREDICTIONS
#
# Compare fitted WUE_T at the dry end with near-normal fitted WUE_T.
# This is descriptive only; significance comes from smooth-term table.
# ============================================================

print("\n" + "=" * 100)
print("4. FITTED DRY-SIDE DIRECTION FROM SAVED GAM PREDICTIONS")
print("=" * 100)

required_pred = [
    "SPEI_timescale",
    "water_class",
    "SPEI_value",
    "predicted_WUE_T"
]

missing = [
    c for c in required_pred
    if c not in pred.columns
]

if missing:
    print("Prediction columns differ from expected.")
    print("Available columns:")
    print(list(pred.columns))
else:

    for eco in ecosystem_order:

        print(f"\n{eco}")
        print("-" * 80)

        rows = []

        for ts in timescale_order:

            x = pred[
                (pred["water_class"] == eco) &
                (pred["SPEI_timescale"] == ts)
            ].copy()

            if x.empty:
                continue

            nn = x[
                (x["SPEI_value"] >= -1) &
                (x["SPEI_value"] <= 1)
            ]

            dry = x[
                x["SPEI_value"] < -1
            ]

            if nn.empty or dry.empty:
                continue

            baseline = nn["predicted_WUE_T"].mean()

            # Most negative SPEI prediction available
            dry_row = dry.loc[
                dry["SPEI_value"].idxmin()
            ]

            dry_pred = dry_row["predicted_WUE_T"]

            pct_change = (
                100 * (dry_pred - baseline) / baseline
                if baseline != 0 else np.nan
            )

            if pct_change > 0:
                direction = "WUE_T higher under dry conditions"
            elif pct_change < 0:
                direction = "WUE_T lower under dry conditions"
            else:
                direction = "no change"

            rows.append({
                "timescale": ts,
                "dry_end_SPEI": dry_row["SPEI_value"],
                "dry_end_pct_change": pct_change,
                "direction": direction
            })

        if rows:
            print(
                pd.DataFrame(rows)
                .round(2)
                .to_string(index=False)
            )

# ============================================================
# 5. OLD LINEAR SITE-LEVEL SLOPE SUMMARY
#
# Useful only to identify where old beta-like interpretation came from.
# Do NOT treat this as the final GAM coefficient.
# ============================================================

print("\n" + "=" * 100)
print("5. EXISTING LINEAR SITE-LEVEL SLOPE SUMMARY")
print("=" * 100)

print(
    slopes[
        [
            "SPEI_timescale",
            "water_class",
            "n_sites",
            "mean_slope",
            "median_slope",
            "pct_negative"
        ]
    ]
    .round(4)
    .to_string(index=False)
)

print("\nDONE. No files were written.")