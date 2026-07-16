# =============================================================================
# Q2 diagnostic:
# Coastal-region differences in WUE_T performance metrics
#
# Saves CSV files only:
#   1) site-level Q2 performance data with coast_region
#   2) long-format data for future coast-based plotting
#   3) coast × metric summary table
#   4) Kruskal-Wallis coast-region test summary
#   5) pairwise coast-region tests
#
# No TXT files.
# No figures.
#
# Important caution:
#   Resistance and recovery are tested, but Gulf and AK have very few valid sites
#   for those metrics. These coast-region comparisons should not be treated as
#   fair/strong coast comparisons unless sample sizes are adequate.
# =============================================================================

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

q2_summary_file = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q2_site_performance_summary.csv"

monthly_file = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs\Q_performace metric_T_ET_relationship"
os.makedirs(output_dir, exist_ok=True)

# -----------------------------------------------------------------------------
# Output CSV files
# -----------------------------------------------------------------------------

site_data_csv = os.path.join(
    output_dir,
    "Q2_performance_by_coast_site_level_data.csv"
)

long_plot_csv = os.path.join(
    output_dir,
    "Q2_performance_by_coast_long_for_plot.csv"
)

coast_summary_csv = os.path.join(
    output_dir,
    "Q2_performance_by_coast_group_summary.csv"
)

coast_tests_csv = os.path.join(
    output_dir,
    "Q2_performance_by_coast_tests.csv"
)

coast_pairwise_csv = os.path.join(
    output_dir,
    "Q2_performance_by_coast_pairwise_tests.csv"
)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]

SMALL_N_WARNING = 5

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

# --------------------------------------------------------------------------
# Proximity-based coast classifier (drop-in replacement)
# --------------------------------------------------------------------------
EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7), (35.4, -120.9),
        (36.6, -121.9), (37.8, -122.5), (39.0, -123.7), (41.0, -124.2),
        (43.5, -124.2), (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8), (29.5, -93.0),
        (29.2, -91.5), (29.3, -90.0), (29.5, -88.8), (30.2, -87.7),
        (30.1, -86.2), (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8), (25.3, -81.1),
        (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1), (28.4, -80.6),
        (29.0, -80.9), (29.9, -81.3), (30.4, -81.4), (31.2, -81.3),
        (32.0, -80.8), (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3), (40.5, -73.9),
        (41.3, -72.0), (42.4, -70.8), (43.5, -70.2)
    ],
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
    """Minimum distance from point to segment AB (equirectangular projection)."""
    import math
    lat0 = math.radians(lat)
    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y])
    p = np.array([0.0, 0.0])
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)
    ab = b - a
    denom = np.dot(ab, ab)
    if denom == 0:
        return float(np.linalg.norm(p - a))
    t = max(0.0, min(1.0, np.dot(p - a, ab) / denom))
    return float(np.linalg.norm(p - (a + t * ab)))

def _distance_to_polyline(lat, lon, polyline):
    """Minimum distance to a coastline polyline."""
    return min(
        _seg_dist(lat, lon, *polyline[i], *polyline[i+1])
        for i in range(len(polyline) - 1)
    )

def classify_coast_region(lat, long):
    """
    Proximity-based coast classifier (drop-in replacement).
    Returns the coast name (string).
    """
    if pd.isna(lat) or pd.isna(long):
        return "Other/Check"
    if lat > 50:
        return "AK Coast"
    
    lat, lon = float(lat), float(long)
    
    dists = {
        name: _distance_to_polyline(lat, lon, poly)
        for name, poly in COAST_POLYLINES.items()
    }
    
    return min(dists, key=dists.get)


def bh_adjust_ignore_nan(pvals):
    """Benjamini-Hochberg FDR correction while preserving NaNs."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)

    finite_mask = np.isfinite(pvals)
    finite_p = pvals[finite_mask]

    if len(finite_p) == 0:
        return out

    n = len(finite_p)
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    adjusted = np.empty(n, dtype=float)

    running_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        running_min = min(running_min, value)
        adjusted[i] = running_min

    adjusted_original = np.empty(n, dtype=float)
    adjusted_original[order] = np.minimum(adjusted, 1.0)
    out[finite_mask] = adjusted_original

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


def run_kw_test(df, group_col, value_col, min_group_n=2):
    """Run Kruskal-Wallis using groups with at least min_group_n valid sites."""
    test_df = df[[group_col, value_col]].replace([np.inf, -np.inf], np.nan).dropna()

    group_counts = test_df.groupby(group_col, observed=True)[value_col].size()
    valid_groups = [
        group for group, n in group_counts.items()
        if n >= min_group_n
    ]

    if len(valid_groups) < 2:
        return np.nan, np.nan, group_counts.to_dict(), valid_groups

    groups = [
        test_df.loc[test_df[group_col] == group, value_col].values
        for group in valid_groups
    ]

    stat, p = kruskal(*groups)
    return stat, p, group_counts.to_dict(), valid_groups


def run_pairwise_tests(df, group_col, value_col, metric, label, min_group_n=2):
    """Pairwise Mann-Whitney tests by coast."""
    test_df = df[[group_col, value_col]].replace([np.inf, -np.inf], np.nan).dropna()

    group_counts = test_df.groupby(group_col, observed=True)[value_col].size()
    valid_groups = [
        group for group, n in group_counts.items()
        if n >= min_group_n
    ]

    rows = []
    pvals = []

    for i, g1 in enumerate(valid_groups):
        for g2 in valid_groups[i + 1:]:
            x = test_df.loc[test_df[group_col] == g1, value_col].values
            y = test_df.loc[test_df[group_col] == g2, value_col].values

            if len(x) >= min_group_n and len(y) >= min_group_n:
                stat, p = mannwhitneyu(x, y, alternative="two-sided")

                rows.append({
                    "metric": metric,
                    "label": label,
                    "comparison": f"{g1} vs {g2}",
                    "group_1": g1,
                    "group_2": g2,
                    "n_1": len(x),
                    "n_2": len(y),
                    "mean_1": np.mean(x),
                    "mean_2": np.mean(y),
                    "median_1": np.median(x),
                    "median_2": np.median(y),
                    "difference_mean_1_minus_2": np.mean(x) - np.mean(y),
                    "p_value": p
                })

                pvals.append(p)

    if not rows:
        return pd.DataFrame()

    pair_df = pd.DataFrame(rows)
    pair_df["p_value_FDR_within_metric"] = bh_adjust_ignore_nan(pair_df["p_value"].values)
    pair_df["sig_FDR_within_metric"] = pair_df["p_value_FDR_within_metric"].apply(sig_label)

    return pair_df


def caution_text(metric, counts_all, counts_large):
    """Add caution about sample size, especially for resistance/recovery."""
    small_groups = [
        coast for coast in COAST_REGION_LEVELS
        if counts_all.get(coast, 0) > 0 and counts_all.get(coast, 0) < SMALL_N_WARNING
    ]

    if metric in ["resistance", "mean_recovery"]:
        return (
            "Caution: coast-region comparison is weak for this event/drought metric "
            "because some coasts have very few valid sites; do not treat as a fair "
            "coast comparison when Gulf/AK sample sizes are very small."
        )

    if small_groups:
        return (
            "Caution: at least one coast has fewer than 5 valid sites; interpret "
            "coast-region comparison carefully."
        )

    return "Sample sizes are acceptable for exploratory coast-region comparison."

# -----------------------------------------------------------------------------
# Load Q2 summary
# -----------------------------------------------------------------------------

q2 = pd.read_csv(q2_summary_file)

# -----------------------------------------------------------------------------
# Load monthly data only to assign coast_region
# -----------------------------------------------------------------------------

monthly = pd.read_csv(monthly_file)

for col in monthly.select_dtypes(include="object").columns:
    monthly[col] = monthly[col].where(
        monthly[col].isna(),
        monthly[col].astype(str).str.strip()
    )
    monthly[col] = monthly[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

site_lookup_base = monthly[
    monthly["site_name"].notna() &
    (monthly["site_name"] != "") &
    monthly["water_class"].notna() &
    monthly["water_class"].isin(ECOSYSTEM_CLASSES) &
    np.isfinite(monthly["lat"]) &
    np.isfinite(monthly["long"])
].copy()

site_lookup = (
    site_lookup_base
    .sort_values(["site_name", "Year", "month"])
    .groupby("site_name", as_index=False)
    .agg(
        lat=("lat", "first"),
        long=("long", "first")
    )
)

site_lookup["coast_region"] = site_lookup.apply(
    lambda row: classify_coast_region(row["lat"], row["long"]),
    axis=1
)

site_lookup["coast_region"] = pd.Categorical(
    site_lookup["coast_region"],
    categories=COAST_REGION_LEVELS,
    ordered=True
)

q2 = q2.merge(
    site_lookup[["site_name", "lat", "long", "coast_region"]],
    on="site_name",
    how="left"
)

q2["coast_region"] = pd.Categorical(
    q2["coast_region"],
    categories=COAST_REGION_LEVELS,
    ordered=True
)

# -----------------------------------------------------------------------------
# Save site-level data with coast_region
# -----------------------------------------------------------------------------

site_cols = [
    "site_name", "water_class", "coast_region", "lat", "long",
    "n_months", "mean_TET",
    "stability",
    "plasticity_slope_SPEI3",
    "plasticity_p95_p05",
    "resistance",
    "mean_recovery"
]

available_site_cols = [c for c in site_cols if c in q2.columns]
q2_site_data = q2[available_site_cols].copy()
q2_site_data.to_csv(site_data_csv, index=False)

# -----------------------------------------------------------------------------
# Create and save long-format data for future coast figure
# -----------------------------------------------------------------------------

long_rows = []

for metric, label, panel_label, notes in metrics_all:
    if metric not in q2.columns:
        continue

    sub = q2[
        ["site_name", "water_class", "coast_region", "lat", "long", "mean_TET", metric]
    ].copy()

    sub = sub.rename(columns={metric: "metric_value"})
    sub["metric"] = metric
    sub["label"] = label
    sub["panel_label_original_Q2"] = panel_label
    sub["notes"] = notes

    sub = sub.replace([np.inf, -np.inf], np.nan)
    sub = sub.dropna(subset=["coast_region", "metric_value"])

    long_rows.append(sub)

q2_long = pd.concat(long_rows, ignore_index=True)
q2_long.to_csv(long_plot_csv, index=False)

# -----------------------------------------------------------------------------
# Coast × metric summary and tests
# -----------------------------------------------------------------------------

summary_rows = []
test_rows = []
pairwise_all = []

for metric, label, panel_label, notes in metrics_all:
    if metric not in q2.columns:
        continue

    test_df = q2[
        ["site_name", "water_class", "coast_region", metric]
    ].copy()

    test_df = test_df.replace([np.inf, -np.inf], np.nan)
    test_df = test_df.dropna(subset=["coast_region", metric])

    # -------------------------------------------------------------------------
    # Coast summary for this metric
    # -------------------------------------------------------------------------

    metric_summary = (
        test_df
        .groupby("coast_region", observed=True)
        .agg(
            n_sites=("site_name", "nunique"),
            mean_value=(metric, "mean"),
            median_value=(metric, "median"),
            sd_value=(metric, "std"),
            min_value=(metric, "min"),
            max_value=(metric, "max")
        )
        .reindex(COAST_REGION_LEVELS)
        .reset_index()
    )

    metric_summary.insert(0, "metric", metric)
    metric_summary.insert(1, "label", label)
    metric_summary.insert(2, "panel_label_original_Q2", panel_label)
    metric_summary.insert(3, "notes", notes)

    summary_rows.append(metric_summary)

    # -------------------------------------------------------------------------
    # Kruskal-Wallis tests
    # -------------------------------------------------------------------------

    stat_all, p_all, counts_all, valid_groups_all = run_kw_test(
        df=test_df,
        group_col="coast_region",
        value_col=metric,
        min_group_n=2
    )

    stat_large, p_large, counts_large, valid_groups_large = run_kw_test(
        df=test_df,
        group_col="coast_region",
        value_col=metric,
        min_group_n=SMALL_N_WARNING
    )

    min_valid_n = min(
        [counts_all.get(coast, 0) for coast in valid_groups_all],
        default=0
    )

    test_rows.append({
        "metric": metric,
        "label": label,
        "panel_label_original_Q2": panel_label,
        "notes": notes,
        "n_sites_total": len(test_df),
        "n_atlantic": counts_all.get("Atlantic Coast", 0),
        "n_pacific": counts_all.get("Pacific Coast", 0),
        "n_gulf": counts_all.get("Gulf Coast", 0),
        "n_ak": counts_all.get("AK Coast", 0),
        "min_valid_coast_n": min_valid_n,
        "groups_all_n_ge_2": ", ".join(valid_groups_all),
        "kw_stat_all_groups": stat_all,
        "kw_p_all_groups": p_all,
        "groups_sensitivity_n_ge_5": ", ".join(valid_groups_large),
        "kw_stat_groups_n_ge_5": stat_large,
        "kw_p_groups_n_ge_5": p_large,
        "sample_size_caution": caution_text(metric, counts_all, counts_large)
    })

    # -------------------------------------------------------------------------
    # Pairwise tests
    # -------------------------------------------------------------------------

    pair_df = run_pairwise_tests(
        df=test_df,
        group_col="coast_region",
        value_col=metric,
        metric=metric,
        label=label,
        min_group_n=2
    )

    if not pair_df.empty:
        pair_df["sample_size_caution"] = caution_text(metric, counts_all, counts_large)
        pairwise_all.append(pair_df)

# Combine summaries
coast_summary = pd.concat(summary_rows, ignore_index=True)
coast_summary.to_csv(coast_summary_csv, index=False)

coast_tests = pd.DataFrame(test_rows)

# FDR correction across the five main coast-region tests
coast_tests["kw_p_all_groups_FDR"] = bh_adjust_ignore_nan(
    coast_tests["kw_p_all_groups"].values
)
coast_tests["sig_all_groups_FDR"] = coast_tests["kw_p_all_groups_FDR"].apply(sig_label)

coast_tests["kw_p_groups_n_ge_5_FDR"] = bh_adjust_ignore_nan(
    coast_tests["kw_p_groups_n_ge_5"].values
)
coast_tests["sig_groups_n_ge_5_FDR"] = coast_tests["kw_p_groups_n_ge_5_FDR"].apply(sig_label)

coast_tests.to_csv(coast_tests_csv, index=False)

if pairwise_all:
    coast_pairwise = pd.concat(pairwise_all, ignore_index=True)
else:
    coast_pairwise = pd.DataFrame()

coast_pairwise.to_csv(coast_pairwise_csv, index=False)

# -----------------------------------------------------------------------------
# Console output
# -----------------------------------------------------------------------------

print("=" * 100)
print("Q2 coast-region performance diagnostic saved")
print("=" * 100)

print(f"Q2 summary file: {q2_summary_file}")
print(f"Monthly file for coast lookup: {monthly_file}")
print(f"Rows in Q2 summary: {len(q2)}")
print(f"Sites in Q2 summary: {q2['site_name'].nunique()}")
print(f"Sites missing coast_region: {q2['coast_region'].isna().sum()}")

print("\nSites by coast region:")
print(
    q2.groupby("coast_region", observed=True)["site_name"]
    .nunique()
    .reindex(COAST_REGION_LEVELS)
    .to_string()
)

print("\nCoast-region test summary:")
print(
    coast_tests[
        [
            "label",
            "n_sites_total",
            "n_atlantic",
            "n_pacific",
            "n_gulf",
            "n_ak",
            "kw_p_all_groups",
            "kw_p_all_groups_FDR",
            "sig_all_groups_FDR",
            "kw_p_groups_n_ge_5",
            "kw_p_groups_n_ge_5_FDR",
            "sig_groups_n_ge_5_FDR",
            "sample_size_caution"
        ]
    ].to_string(index=False, float_format=lambda x: f"{x:.4g}")
)

print("\nSaved CSV files:")
print(f"- {site_data_csv}")
print(f"- {long_plot_csv}")
print(f"- {coast_summary_csv}")
print(f"- {coast_tests_csv}")
print(f"- {coast_pairwise_csv}")

print("\nDone. No TXT files or figures were created.")