# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:56:22 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
US-EKH DESCRIPTIVE DROUGHT DIAGNOSTIC
Internal reviewer check only.

READ ONLY:
- Does not save files.
- Does not change the manuscript Q3 workflow.
- Uses the 15 US-EKH observations that entered Q3 model_base.

Purpose:
1. Examine WUE_T against each SPEI timescale.
2. Compare WUE_T during dry (SPEI <= -1) vs non-dry conditions.
3. Calculate Spearman correlations.
4. Compare the same calendar months among years.
5. Inspect whether 2022 persistent drought coincided with relatively low WUE_T.

IMPORTANT:
This is descriptive because US-EKH only has 2022-2024 data.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ============================================================
# SETTINGS
# ============================================================

Q3_DIR = (
    r"M:\Research\WUE_CUE\WUE_manuscript_version6"
    r"\Q3\Q3_WUE_T_SPEI_sensitivity_outputs"
)

INPUT_FILE = os.path.join(
    Q3_DIR,
    "Q3_WUE_T_SPEI_model_base_wide.csv"
)

SITE = "US-EKH"

SPEI_COLS = [
    "SPEI_1",
    "SPEI_3",
    "SPEI_6",
    "SPEI_12",
    "SPEI_24",
    "SPEI_36",
    "SPEI_48"
]

DRY_THRESHOLD = -1.0

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT_FILE)

site = df[
    df["site_name"].astype(str).str.strip() == SITE
].copy()

if site.empty:
    raise ValueError(f"{SITE} not found in {INPUT_FILE}")

site = site.sort_values(["Year", "month"])

print("=" * 105)
print("US-EKH DESCRIPTIVE DROUGHT DIAGNOSTIC")
print("Elkhorn Slough Hester Marsh")
print("=" * 105)

print(f"\nUsable WUE_T observations: {len(site)}")
print(
    f"Years: {int(site['Year'].min())}–"
    f"{int(site['Year'].max())}"
)
print(
    "Ecosystem:",
    site["water_class"].dropna().unique().tolist()
)

# ============================================================
# 1. COMPLETE CHRONOLOGICAL TABLE
# ============================================================

print("\n" + "=" * 105)
print("1. CHRONOLOGICAL WUE_T AND SPEI")
print("=" * 105)

cols = [
    "Year", "month", "WUE_T",
    "SPEI_1", "SPEI_3", "SPEI_6",
    "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"
]

print(
    site[cols]
    .round(3)
    .to_string(index=False)
)

# ============================================================
# 2. SPEARMAN ASSOCIATION: WUE_T VS EACH SPEI TIMESCALE
# ============================================================

print("\n" + "=" * 105)
print("2. SPEARMAN ASSOCIATION BETWEEN WUE_T AND SPEI")
print("=" * 105)

cor_rows = []

for sp in SPEI_COLS:

    sub = site[["WUE_T", sp]].dropna()

    if len(sub) >= 4 and sub[sp].nunique() >= 4:
        rho, p = spearmanr(
            sub[sp],
            sub["WUE_T"]
        )
    else:
        rho, p = np.nan, np.nan

    cor_rows.append({
        "SPEI_timescale": sp,
        "n": len(sub),
        "rho": rho,
        "p_value": p
    })

cor_table = pd.DataFrame(cor_rows)

print(
    cor_table
    .round(3)
    .to_string(index=False)
)

print("""
Interpretation:
  rho > 0 : lower SPEI (drier) tends to coincide with lower WUE_T.
  rho < 0 : lower SPEI (drier) tends to coincide with higher WUE_T.

Do not interpret p-values strongly with only 15 observations.
""")

# ============================================================
# 3. DRY VS NON-DRY DESCRIPTIVE COMPARISON
# ============================================================

print("\n" + "=" * 105)
print("3. WUE_T DURING DRY VS NON-DRY CONDITIONS")
print("=" * 105)

comparison_rows = []

for sp in SPEI_COLS:

    sub = site[
        ["WUE_T", sp]
    ].dropna().copy()

    dry = sub[
        sub[sp] <= DRY_THRESHOLD
    ]["WUE_T"]

    nondry = sub[
        sub[sp] > DRY_THRESHOLD
    ]["WUE_T"]

    near_normal = sub[
        (sub[sp] >= -1) &
        (sub[sp] <= 1)
    ]["WUE_T"]

    comparison_rows.append({
        "SPEI_timescale": sp,

        "n_dry": len(dry),
        "mean_WUE_dry": (
            dry.mean() if len(dry) else np.nan
        ),
        "median_WUE_dry": (
            dry.median() if len(dry) else np.nan
        ),

        "n_non_dry": len(nondry),
        "mean_WUE_non_dry": (
            nondry.mean() if len(nondry) else np.nan
        ),

        "n_near_normal": len(near_normal),
        "mean_WUE_near_normal": (
            near_normal.mean()
            if len(near_normal) else np.nan
        ),

        "dry_minus_non_dry": (
            dry.mean() - nondry.mean()
            if len(dry) and len(nondry)
            else np.nan
        )
    })

comparison = pd.DataFrame(comparison_rows)

print(
    comparison
    .round(3)
    .to_string(index=False)
)

# ============================================================
# 4. SAME-MONTH COMPARISON AMONG YEARS
# ============================================================

print("\n" + "=" * 105)
print("4. SAME CALENDAR MONTH COMPARISON AMONG YEARS")
print("=" * 105)

pivot_wue = site.pivot_table(
    index="month",
    columns="Year",
    values="WUE_T",
    aggfunc="first"
)

print("\nObserved WUE_T:")
print(
    pivot_wue
    .round(3)
    .to_string()
)

# ============================================================
# 5. DIRECT 2022 VS 2024 SAME-MONTH COMPARISON
# ============================================================

print("\n" + "=" * 105)
print("5. 2022 DROUGHT VS 2024 SAME-MONTH WUE_T")
print("=" * 105)

y22 = site[
    site["Year"] == 2022
][
    ["month", "WUE_T"]
].rename(
    columns={"WUE_T": "WUE_T_2022"}
)

y24 = site[
    site["Year"] == 2024
][
    ["month", "WUE_T"]
].rename(
    columns={"WUE_T": "WUE_T_2024"}
)

paired = y22.merge(
    y24,
    on="month",
    how="inner"
)

paired["difference_2022_minus_2024"] = (
    paired["WUE_T_2022"] -
    paired["WUE_T_2024"]
)

paired["pct_lower_2022_vs_2024"] = (
    100 *
    (
        paired["WUE_T_2024"] -
        paired["WUE_T_2022"]
    ) /
    paired["WUE_T_2024"]
)

print(
    paired
    .round(2)
    .to_string(index=False)
)

if len(paired) > 0:

    print(
        "\nMean WUE_T for common months:"
    )
    print(
        f"  2022 = "
        f"{paired['WUE_T_2022'].mean():.3f}"
    )
    print(
        f"  2024 = "
        f"{paired['WUE_T_2024'].mean():.3f}"
    )

    pct = (
        100 *
        (
            paired["WUE_T_2024"].mean() -
            paired["WUE_T_2022"].mean()
        ) /
        paired["WUE_T_2024"].mean()
    )

    print(
        f"  2022 was {pct:.1f}% lower "
        f"than 2024 across common months."
    )

# ============================================================
# 6. JULY-AUGUST DETAIL
# ============================================================

print("\n" + "=" * 105)
print("6. JULY-AUGUST DETAIL")
print("=" * 105)

summer = site[
    site["month"].isin([7, 8])
][cols].copy()

print(
    summer
    .round(3)
    .to_string(index=False)
)

# ============================================================
# 7. LONG-TERM DROUGHT OBSERVATIONS
# ============================================================

print("\n" + "=" * 105)
print("7. MONTHS WITH PERSISTENT DROUGHT")
print("Criterion: SPEI-24 <= -1 OR SPEI-48 <= -1")
print("=" * 105)

persistent = site[
    (site["SPEI_24"] <= -1) |
    (site["SPEI_48"] <= -1)
][
    [
        "Year", "month", "WUE_T",
        "SPEI_3", "SPEI_12",
        "SPEI_24", "SPEI_36", "SPEI_48"
    ]
].copy()

print(
    persistent
    .round(3)
    .to_string(index=False)
)

# ============================================================
# 8. INTERNAL CAUTION
# ============================================================

print("\n" + "=" * 105)
print("8. INTERNAL INTERPRETATION")
print("=" * 105)

print("""
This analysis is descriptive only.

Potentially supportive pattern:
  - WUE_T is lower during persistent drought observations;
  - WUE_T increases in comparable calendar months after
    long-timescale SPEI recovers;
  - positive WUE_T-SPEI correlations are strongest at
    SPEI-12 to SPEI-48.

Do NOT claim direct multi-year drought causation because:
  - only 15 WUE_T observations are available;
  - only July-August WUE_T is available during the severe 2022 drought;
  - 2022 February-June WUE_T is missing;
  - no independent salinity or vegetation-cover measurement
    is included in this diagnostic.

Use the NERR / Elkhorn Slough vegetation literature as the
direct ecological evidence. This site diagnostic can only tell
us whether the flux-tower record is qualitatively consistent
with that evidence.
""")

print("=" * 105)
print("DONE. No files were written.")
print("=" * 105)