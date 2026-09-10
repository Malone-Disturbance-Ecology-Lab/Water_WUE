import pandas as pd
import numpy as np
from datetime import datetime, timedelta

# ============================================================
# PATHS
# ============================================================

FINAL_MONTHLY = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
)

Q3_MODEL_DATA = (
    r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
    r"\Q3_WUE_T_SPEI_sensitivity_outputs"
    r"\Q3_WUE_T_SPEI_model_data_long.csv"
)

APRIL_START_DOY = 91      # April 1
OCTOBER_END_DOY = 304     # October 31

# ============================================================
# LOAD FINAL Q3 SITE LIST
# ============================================================

q3 = pd.read_csv(Q3_MODEL_DATA)

if "site_name" not in q3.columns:
    raise ValueError("Q3 model dataset does not contain 'site_name'.")

final_sites = set(
    q3["site_name"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)

print("=" * 90)
print("APRIL-OCTOBER VS GPP-DERIVED GROWING-SEASON DIAGNOSTIC")
print("FINAL Q3 ANALYTICAL SITES ONLY")
print("=" * 90)

print(f"Sites in Q3 model dataset: {len(final_sites)}")

# ============================================================
# LOAD MONTHLY DATA AND FILTER TO FINAL Q3 SITES
# ============================================================

df = pd.read_csv(FINAL_MONTHLY)

required = ["site_name", "Year", "avg_sos", "avg_eos"]
missing = [c for c in required if c not in df.columns]

if missing:
    raise ValueError(f"Missing required columns in monthly dataset: {missing}")

df["site_name"] = (
    df["site_name"]
    .astype(str)
    .str.strip()
)

print(f"Rows in monthly dataset before filtering: {len(df)}")
print(f"Sites in monthly dataset before filtering: {df['site_name'].nunique()}")

df = df[df["site_name"].isin(final_sites)].copy()

print(f"Rows after Q3 site filter: {len(df)}")
print(f"Sites after Q3 site filter: {df['site_name'].nunique()}")

missing_q3_sites = sorted(
    final_sites - set(df["site_name"].unique())
)

if missing_q3_sites:
    print("\nWARNING: Q3 sites not found in monthly dataset:")
    print(", ".join(missing_q3_sites))

# ============================================================
# CHECK OPTIONAL COAST COLUMN
# ============================================================

possible_coast_cols = [
    "coast_region",
    "coast_final",
    "coast",
    "region"
]

coast_col = None

for c in possible_coast_cols:
    if c in df.columns:
        coast_col = c
        break

if coast_col is not None:
    print(f"Using coastal-region column: {coast_col}")
else:
    print("No coastal-region column found in monthly dataset.")

# ============================================================
# ONE SOS/EOS VALUE PER SITE-YEAR
# ============================================================

cols = ["site_name", "Year", "avg_sos", "avg_eos"]

if coast_col is not None:
    cols.append(coast_col)

tmp = df[cols].copy()

tmp["avg_sos"] = pd.to_numeric(
    tmp["avg_sos"],
    errors="coerce"
)

tmp["avg_eos"] = pd.to_numeric(
    tmp["avg_eos"],
    errors="coerce"
)

tmp = tmp.dropna(
    subset=["avg_sos", "avg_eos"]
)

# Remove impossible values
tmp = tmp[
    (tmp["avg_sos"] >= 1) &
    (tmp["avg_sos"] <= 366) &
    (tmp["avg_eos"] >= 1) &
    (tmp["avg_eos"] <= 366) &
    (tmp["avg_eos"] >= tmp["avg_sos"])
].copy()

agg_dict = {
    "avg_sos": "median",
    "avg_eos": "median"
}

if coast_col is not None:
    agg_dict[coast_col] = "first"

pheno = (
    tmp
    .groupby(["site_name", "Year"], as_index=False)
    .agg(agg_dict)
)

print(f"Sites with valid SOS/EOS: {pheno['site_name'].nunique()}")
print(f"Valid analytical site-years: {len(pheno)}")

missing_pheno_sites = sorted(
    final_sites - set(pheno["site_name"].unique())
)

if missing_pheno_sites:
    print("\nQ3 sites without valid SOS/EOS:")
    print(", ".join(missing_pheno_sites))

# ============================================================
# DOY TO MONTH
# ============================================================

def doy_to_month(doy):
    ref = datetime(2021, 1, 1)
    date = ref + timedelta(
        days=int(round(doy)) - 1
    )
    return date.month

pheno["sos_month"] = (
    pheno["avg_sos"]
    .apply(doy_to_month)
)

pheno["eos_month"] = (
    pheno["avg_eos"]
    .apply(doy_to_month)
)

# ============================================================
# GROWING-SEASON LENGTH
# ============================================================

pheno["gs_length_days"] = (
    pheno["avg_eos"] -
    pheno["avg_sos"] +
    1
)

# ============================================================
# APRIL-OCTOBER OVERLAP
# ============================================================

overlap_start = np.maximum(
    pheno["avg_sos"],
    APRIL_START_DOY
)

overlap_end = np.minimum(
    pheno["avg_eos"],
    OCTOBER_END_DOY
)

pheno["apr_oct_overlap_days"] = np.maximum(
    0,
    overlap_end - overlap_start + 1
)

pheno["fraction_gs_in_apr_oct"] = (
    pheno["apr_oct_overlap_days"] /
    pheno["gs_length_days"]
)

pheno["any_overlap"] = (
    pheno["apr_oct_overlap_days"] > 0
)

pheno["fully_within_apr_oct"] = (
    (pheno["avg_sos"] >= APRIL_START_DOY) &
    (pheno["avg_eos"] <= OCTOBER_END_DOY)
)

# ============================================================
# OVERALL SUMMARY
# ============================================================

print("\n" + "=" * 90)
print("OVERALL PHENOLOGY")
print("=" * 90)

median_sos = pheno["avg_sos"].median()
median_eos = pheno["avg_eos"].median()
median_length = pheno["gs_length_days"].median()

median_overlap = (
    100 *
    pheno["fraction_gs_in_apr_oct"].median()
)

mean_overlap = (
    100 *
    pheno["fraction_gs_in_apr_oct"].mean()
)

any_overlap_pct = (
    100 *
    pheno["any_overlap"].mean()
)

fully_within_pct = (
    100 *
    pheno["fully_within_apr_oct"].mean()
)

weighted_overlap = (
    pheno["apr_oct_overlap_days"].sum() /
    pheno["gs_length_days"].sum()
)

print(
    f"Median SOS: DOY {median_sos:.0f} "
    f"(approximately month "
    f"{int(round(pheno['sos_month'].median()))})"
)

print(
    f"Median EOS: DOY {median_eos:.0f} "
    f"(approximately month "
    f"{int(round(pheno['eos_month'].median()))})"
)

print(
    f"Median growing-season length: "
    f"{median_length:.0f} days"
)

print(
    f"Median fraction of GPP-defined growing-season days "
    f"within April-October: "
    f"{median_overlap:.1f}%"
)

print(
    f"Mean fraction of GPP-defined growing-season days "
    f"within April-October: "
    f"{mean_overlap:.1f}%"
)

print(
    f"Site-years with any April-October overlap: "
    f"{any_overlap_pct:.1f}%"
)

print(
    f"Site-years entirely contained within April-October: "
    f"{fully_within_pct:.1f}%"
)

print(
    f"Fraction of all GPP-defined active-season days "
    f"falling within April-October: "
    f"{100 * weighted_overlap:.1f}%"
)

# ============================================================
# SOS MONTH DISTRIBUTION
# ============================================================

print("\n" + "=" * 90)
print("SOS MONTH DISTRIBUTION")
print("=" * 90)

print(
    pheno["sos_month"]
    .value_counts()
    .sort_index()
    .rename_axis("Month")
    .to_string()
)

# ============================================================
# EOS MONTH DISTRIBUTION
# ============================================================

print("\n" + "=" * 90)
print("EOS MONTH DISTRIBUTION")
print("=" * 90)

print(
    pheno["eos_month"]
    .value_counts()
    .sort_index()
    .rename_axis("Month")
    .to_string()
)

# ============================================================
# SITE-LEVEL SUMMARY
# ============================================================

site_summary = (
    pheno
    .groupby("site_name", as_index=False)
    .agg(
        n_years=("Year", "nunique"),
        median_sos=("avg_sos", "median"),
        median_eos=("avg_eos", "median"),
        median_length=("gs_length_days", "median"),
        mean_apr_oct_fraction=(
            "fraction_gs_in_apr_oct",
            "mean"
        )
    )
)

site_summary["mean_apr_oct_pct"] = (
    100 *
    site_summary["mean_apr_oct_fraction"]
)

# ============================================================
# LOWEST APRIL-OCTOBER COVERAGE
# ============================================================

print("\n" + "=" * 90)
print("15 SITES WITH LOWEST APRIL-OCTOBER COVERAGE")
print("=" * 90)

lowest = (
    site_summary
    .sort_values("mean_apr_oct_pct")
    .head(15)
)

for _, r in lowest.iterrows():
    print(
        f"{r['site_name']:8s} | "
        f"years={int(r['n_years']):2d} | "
        f"median SOS={r['median_sos']:.0f} | "
        f"median EOS={r['median_eos']:.0f} | "
        f"Apr-Oct coverage="
        f"{r['mean_apr_oct_pct']:.1f}%"
    )

# ============================================================
# LONG ACTIVE-SEASON SITES
# ============================================================

print("\n" + "=" * 90)
print("SITES WITH LONG GPP-DEFINED ACTIVE SEASONS")
print("=" * 90)

long_sites = (
    site_summary[
        site_summary["median_length"] >= 250
    ]
    .sort_values(
        "median_length",
        ascending=False
    )
)

if long_sites.empty:
    print(
        "No sites with median "
        "GPP-defined active season >= 250 days."
    )
else:
    for _, r in long_sites.iterrows():
        print(
            f"{r['site_name']:8s} | "
            f"median season="
            f"{r['median_length']:.0f} d | "
            f"median SOS="
            f"{r['median_sos']:.0f} | "
            f"median EOS="
            f"{r['median_eos']:.0f} | "
            f"Apr-Oct coverage="
            f"{r['mean_apr_oct_pct']:.1f}%"
        )

# ============================================================
# COASTAL REGION SUMMARY
# ============================================================

if coast_col is not None:

    coast_summary = (
        pheno
        .groupby(coast_col, as_index=False)
        .agg(
            n_sites=("site_name", "nunique"),
            n_site_years=("Year", "count"),
            median_sos=("avg_sos", "median"),
            median_eos=("avg_eos", "median"),
            median_length=("gs_length_days", "median"),
            median_apr_oct_fraction=(
                "fraction_gs_in_apr_oct",
                "median"
            ),
            mean_apr_oct_fraction=(
                "fraction_gs_in_apr_oct",
                "mean"
            ),
            total_overlap_days=(
                "apr_oct_overlap_days",
                "sum"
            ),
            total_gs_days=(
                "gs_length_days",
                "sum"
            )
        )
    )

    coast_summary["weighted_apr_oct_fraction"] = (
        coast_summary["total_overlap_days"] /
        coast_summary["total_gs_days"]
    )

    print("\n" + "=" * 90)
    print("SUMMARY BY COASTAL REGION")
    print("=" * 90)

    for _, r in coast_summary.iterrows():
        print(
            f"{str(r[coast_col]):16s} | "
            f"sites={int(r['n_sites']):2d} | "
            f"site-years="
            f"{int(r['n_site_years']):3d} | "
            f"median SOS="
            f"{r['median_sos']:.0f} | "
            f"median EOS="
            f"{r['median_eos']:.0f} | "
            f"median season="
            f"{r['median_length']:.0f} d | "
            f"median coverage="
            f"{100*r['median_apr_oct_fraction']:.1f}% | "
            f"weighted coverage="
            f"{100*r['weighted_apr_oct_fraction']:.1f}%"
        )

else:
    print(
        "\nNOTE: No coastal-region column was found, "
        "so no regional summary was calculated."
    )

# ============================================================
# REVIEWER-USEFUL SUMMARY
# ============================================================

print("\n" + "=" * 90)
print("REVIEWER-USEFUL SUMMARY")
print("=" * 90)

print(
    f"Across {pheno['site_name'].nunique()} Q3 analytical sites "
    f"and {len(pheno)} analytical site-years, "
    f"April-October contained a median of "
    f"{median_overlap:.1f}% "
    f"of each site-year's GPP-defined active-season days."
)

print(
    f"Across all site-years combined, "
    f"{100 * weighted_overlap:.1f}% "
    f"of GPP-defined active-season days "
    f"fell within April-October."
)

print("=" * 90)