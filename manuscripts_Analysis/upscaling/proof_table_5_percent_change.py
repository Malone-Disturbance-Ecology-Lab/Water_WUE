# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 12:21:05 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
reviewer_5pct_WUE_Table1_validation.py

Purpose
-------
Validate the ±5% WUE_T ecological-response cutoff against the
site-level performance metrics reported in Table 1.

The analysis:
1. Uses observed site-level SPEI-3 responses from Q3.
2. Retains dry observations: SPEI-3 <= -1.
3. Calculates mean dry-period WUE_T percent change for each site.
4. Classifies sites as:
      <5% mean dry response magnitude
      >=5% mean dry response magnitude
5. Links these groups to the Q2 Table 1 metrics:
      stability
      linear plasticity slope (SPEI-3)
      plasticity range (95th/5th)
      drought resistance
      drought recovery
6. Reports median [IQR] and Mann-Whitney tests.
7. Produces both all-ecosystem and Upland-only comparisons.

No models are refit.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

# ============================================================
# PATHS
# ============================================================

BASE = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6")

Q2_FILE = (
    BASE
    / "Q2"
    / "Q2_WUE_performance_outputs"
    / "Q2_site_performance_summary.csv"
)

Q3_FILE = (
    BASE
    / "Q3"
    / "Q3_WUE_T_SPEI_sensitivity_outputs"
    / "Q3_WUE_T_SPEI_site_response_from_near_normal.csv"
)

OUTPUT_DIR = BASE / "upscaling" / "aug_figure"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TABLE_FILE = OUTPUT_DIR / "reviewer_5pct_WUE_Table1_validation.csv"
SITE_FILE = OUTPUT_DIR / "reviewer_5pct_WUE_site_classification.csv"

# ============================================================
# SETTINGS
# ============================================================

THRESHOLD = 5.0

METRICS = {
    "stability": "Stability",
    "plasticity_slope_SPEI3": "Plasticity slope (SPEI-3)",
    "plasticity_p95_p05": "Plasticity range (95th/5th)",
    "resistance": "Drought resistance",
    "mean_recovery": "Drought recovery",
}

# ============================================================
# HELPERS
# ============================================================

def bh_fdr(pvalues):
    """
    Benjamini-Hochberg FDR adjustment.
    Keeps the script independent of statsmodels.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full(len(pvalues), np.nan)

    valid = np.isfinite(pvalues)
    p = pvalues[valid]

    if len(p) == 0:
        return adjusted

    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)

    adj_ranked = ranked * n / np.arange(1, n + 1)

    # Enforce monotonicity from largest p downward
    adj_ranked = np.minimum.accumulate(adj_ranked[::-1])[::-1]
    adj_ranked = np.minimum(adj_ranked, 1.0)

    temp = np.empty(n)
    temp[order] = adj_ranked

    adjusted[np.where(valid)[0]] = temp

    return adjusted


def median_iqr(series):
    x = pd.to_numeric(series, errors="coerce").dropna()

    if len(x) == 0:
        return np.nan, np.nan, np.nan, 0

    return (
        float(x.median()),
        float(x.quantile(0.25)),
        float(x.quantile(0.75)),
        int(len(x)),
    )


def format_median_iqr(median, q25, q75):
    if not np.isfinite(median):
        return ""

    return f"{median:.3f} [{q25:.3f}, {q75:.3f}]"


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 80)
print("5% WUE_T THRESHOLD VALIDATION AGAINST TABLE 1 METRICS")
print("=" * 80)

if not Q2_FILE.exists():
    raise FileNotFoundError(f"Missing Q2 file:\n{Q2_FILE}")

if not Q3_FILE.exists():
    raise FileNotFoundError(f"Missing Q3 file:\n{Q3_FILE}")

q2 = pd.read_csv(Q2_FILE)
q3 = pd.read_csv(Q3_FILE)

print(f"\nQ2 site-performance rows: {len(q2)}")
print(f"Q3 site-response rows: {len(q3)}")

# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

q2_required = [
    "site_name",
    "water_class",
    "stability",
    "plasticity_slope_SPEI3",
    "plasticity_p95_p05",
    "resistance",
    "mean_recovery",
]

q3_required = [
    "site_name",
    "water_class",
    "coast_region",
    "SPEI_timescale",
    "SPEI_value",
    "WUE_T_pct_change",
]

missing_q2 = [c for c in q2_required if c not in q2.columns]
missing_q3 = [c for c in q3_required if c not in q3.columns]

if missing_q2:
    raise ValueError(f"Missing Q2 columns: {missing_q2}")

if missing_q3:
    raise ValueError(f"Missing Q3 columns: {missing_q3}")

# ============================================================
# SITE-LEVEL SPEI-3 DRY RESPONSE
# ============================================================

dry = q3[
    (q3["SPEI_timescale"] == "SPEI_3")
    & (q3["SPEI_value"] <= -1)
    & np.isfinite(q3["WUE_T_pct_change"])
].copy()

site_response = (
    dry
    .groupby(
        ["site_name", "water_class", "coast_region"],
        as_index=False
    )
    .agg(
        n_dry_months=("WUE_T_pct_change", "size"),
        mean_dry_WUE_T_pct_change=("WUE_T_pct_change", "mean"),
        median_dry_WUE_T_pct_change=("WUE_T_pct_change", "median"),
        mean_abs_dry_WUE_T_pct_change=(
            "WUE_T_pct_change",
            lambda x: np.mean(np.abs(x))
        ),
        max_abs_dry_WUE_T_pct_change=(
            "WUE_T_pct_change",
            lambda x: np.max(np.abs(x))
        ),
        pct_dry_months_crossing_5pct=(
            "WUE_T_pct_change",
            lambda x: np.mean(np.abs(x) >= THRESHOLD) * 100
        ),
    )
)

# Site-level response magnitude used for reviewer validation
site_response["abs_mean_dry_WUE_T_pct_change"] = (
    site_response["mean_dry_WUE_T_pct_change"].abs()
)

site_response["response_group"] = np.where(
    site_response["abs_mean_dry_WUE_T_pct_change"] >= THRESHOLD,
    "≥5% mean dry WUE_T change",
    "<5% mean dry WUE_T change",
)

# ============================================================
# MERGE WITH TABLE 1 METRICS
# ============================================================

q2_keep = q2[q2_required].copy()

merged = site_response.merge(
    q2_keep,
    on=["site_name", "water_class"],
    how="left",
    validate="one_to_one",
)

merged.to_csv(SITE_FILE, index=False)

print("\nSite classification:")
print(
    merged["response_group"]
    .value_counts()
    .to_string()
)

print("\nUpland only:")
print(
    merged.loc[merged["water_class"] == "Upland", "response_group"]
    .value_counts()
    .to_string()
)

# ============================================================
# CREATE VALIDATION TABLE
# ============================================================

results = []

scopes = {
    "All ecosystems": merged,
    "Upland only": merged[merged["water_class"] == "Upland"].copy(),
}

for scope_name, dat in scopes.items():

    scope_rows = []

    for metric_col, metric_label in METRICS.items():

        below = pd.to_numeric(
            dat.loc[
                dat["response_group"] == "<5% mean dry WUE_T change",
                metric_col
            ],
            errors="coerce",
        ).dropna()

        above = pd.to_numeric(
            dat.loc[
                dat["response_group"] == "≥5% mean dry WUE_T change",
                metric_col
            ],
            errors="coerce",
        ).dropna()

        med_b, q25_b, q75_b, n_b = median_iqr(below)
        med_a, q25_a, q75_a, n_a = median_iqr(above)

        if len(below) > 0 and len(above) > 0:
            test = mannwhitneyu(
                below,
                above,
                alternative="two-sided"
            )
            p_value = float(test.pvalue)
        else:
            p_value = np.nan

        scope_rows.append({
            "Scope": scope_name,
            "Metric": metric_label,

            "<5% n": n_b,
            "<5% median": med_b,
            "<5% Q25": q25_b,
            "<5% Q75": q75_b,
            "<5% median [IQR]": format_median_iqr(
                med_b, q25_b, q75_b
            ),

            "≥5% n": n_a,
            "≥5% median": med_a,
            "≥5% Q25": q25_a,
            "≥5% Q75": q75_a,
            "≥5% median [IQR]": format_median_iqr(
                med_a, q25_a, q75_a
            ),

            "Mann-Whitney p": p_value,
        })

    scope_df = pd.DataFrame(scope_rows)

    scope_df["BH-FDR adjusted p"] = bh_fdr(
        scope_df["Mann-Whitney p"].values
    )

    results.append(scope_df)

validation_table = pd.concat(results, ignore_index=True)

# ============================================================
# SAVE
# ============================================================

validation_table.to_csv(TABLE_FILE, index=False)

print("\n" + "=" * 80)
print("VALIDATION TABLE")
print("=" * 80)

display_cols = [
    "Scope",
    "Metric",
    "<5% n",
    "<5% median [IQR]",
    "≥5% n",
    "≥5% median [IQR]",
    "Mann-Whitney p",
    "BH-FDR adjusted p",
]

print(
    validation_table[display_cols]
    .to_string(index=False)
)

print("\nSaved:")
print(TABLE_FILE)
print(SITE_FILE)

print("\nDONE.")