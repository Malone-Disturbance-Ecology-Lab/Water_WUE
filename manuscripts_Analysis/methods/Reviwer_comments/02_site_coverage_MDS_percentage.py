# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:54:35 2026

@author: ammar
"""

import os
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================
MONTHLY_DATA = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
)

BEFORE_MDS_DIR = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\ameri_data\reddy_gaps\blended_gaps2"
)

AFTER_MDS_DIR = r"M:\Research\WUE_CUE\ameri_data\ameri_fill"


# ============================================================
# HELPERS
# ============================================================
def find_time_col(df):
    candidates = [
        "TIMESTAMP",
        "DateTime",
        "datetime",
        "timestamp",
        "TIMESTAMP_START"
    ]

    for c in candidates:
        if c in df.columns:
            return c

    return None


def longest_true_run(df, mask_col):
    """
    Longest continuous 30-minute run where mask_col is True.
    """
    if df.empty:
        return None

    d = df.sort_values("time").reset_index(drop=True).copy()
    mask = d[mask_col].fillna(False).astype(bool)

    if not mask.any():
        return None

    # Start new group when:
    # 1. True/False status changes, or
    # 2. timestamps are not exactly 30 min apart
    time_break = d["time"].diff().ne(pd.Timedelta(minutes=30))
    groups = ((mask != mask.shift()) | time_break).cumsum()

    runs = []

    for _, g in d[mask].groupby(groups[mask]):
        runs.append({
            "start": g["time"].iloc[0],
            "end": g["time"].iloc[-1],
            "n_hh": len(g),
            "days": len(g) / 48
        })

    if not runs:
        return None

    return max(runs, key=lambda x: x["n_hh"])


# ============================================================
# 1. LOAD FINAL ANALYTICAL DATASET
#    EXACT SAME FILTERS AS YOUR STUDY-AREA SCRIPT
# ============================================================
df = pd.read_csv(MONTHLY_DATA)

for col in df.select_dtypes(include="object").columns:
    df[col] = df[col].where(
        df[col].isna(),
        df[col].astype(str).str.strip()
    )

    df[col] = df[col].replace({
        "": np.nan,
        "nan": np.nan,
        "NaN": np.nan
    })


required = [
    "site_name",
    "Year",
    "month",
    "water_class",
    "lat",
    "long",
    "Trans_ratio",
    "WUE_tra"
]

study = df.dropna(subset=required).copy()

study = study[
    study["water_class"].isin(
        ["Upland", "Freshwater", "Saline"]
    )
].copy()

study = study[
    np.isfinite(study["lat"]) &
    np.isfinite(study["long"]) &
    np.isfinite(study["Trans_ratio"]) &
    (study["Trans_ratio"] >= 0) &
    (study["Trans_ratio"] <= 1) &
    np.isfinite(study["WUE_tra"])
].copy()

study["Year"] = study["Year"].astype(int)
study["month"] = study["month"].astype(int)


# ============================================================
# 2. IDENTIFY COAST / ALASKA
# ============================================================
coast_col = None

for candidate in [
    "coast",
    "Coast",
    "coastline",
    "Coastline",
    "region",
    "Region",
    "coastal_region"
]:
    if candidate in study.columns:
        coast_col = candidate
        break


site_meta = (
    study.groupby("site_name", as_index=False)
    .agg(
        lat=("lat", "first"),
        long=("long", "first"),
        ecosystem=("water_class", "first")
    )
)

if coast_col is not None:

    coast_meta = (
        study.groupby("site_name")[coast_col]
        .first()
        .reset_index()
        .rename(columns={coast_col: "coast"})
    )

    site_meta = site_meta.merge(
        coast_meta,
        on="site_name",
        how="left"
    )

    site_meta["is_alaska"] = (
        site_meta["coast"]
        .astype(str)
        .str.strip()
        .str.lower()
        .eq("alaska")
    )

    print(f"Using '{coast_col}' to identify Alaska sites.")

else:

    # Safe geographic fallback:
    # all contiguous U.S. sites are south of ~50°N.
    site_meta["coast"] = np.where(
        site_meta["lat"] > 50,
        "Alaska",
        "Other"
    )

    site_meta["is_alaska"] = site_meta["lat"] > 50

    print(
        "No coast/region column found; "
        "using latitude > 50°N to identify Alaska."
    )


# ============================================================
# 3. EXACT FINAL SITE-MONTHS USED
# ============================================================
used_months = (
    study[
        ["site_name", "Year", "month"]
    ]
    .drop_duplicates()
    .copy()
)

study_sites = sorted(
    used_months["site_name"].unique()
)

print("\n" + "=" * 110)
print("FINAL ANALYSIS DATA USED FOR GAP-FILL DIAGNOSTIC")
print("=" * 110)

print(f"Final sites:       {len(study_sites)}")
print(f"Final site-months: {len(used_months):,}")
print(
    f"Study years:       "
    f"{used_months['Year'].min()}–{used_months['Year'].max()}"
)

alaska_sites = sorted(
    site_meta.loc[
        site_meta["is_alaska"],
        "site_name"
    ].unique()
)

print(f"Alaska sites:      {len(alaska_sites)}")
print(", ".join(alaska_sites))


# ============================================================
# 4. COMPARE BEFORE-MDS VS AFTER-MDS
# ============================================================
site_year_results = []
site_month_results = []

mapping = {
    "NEE": "NEE_f",
    "LE": "LE_f"
}


for i, site in enumerate(study_sites, start=1):

    before_file = os.path.join(
        BEFORE_MDS_DIR,
        f"gaps_blend_{site}.csv"
    )

    after_file = os.path.join(
        AFTER_MDS_DIR,
        f"{site}_fill.csv"
    )

    if not os.path.exists(before_file):
        print(f"{site}: BEFORE-MDS file missing")
        continue

    if not os.path.exists(after_file):
        print(f"{site}: AFTER-MDS file missing")
        continue

    before = pd.read_csv(before_file)
    after = pd.read_csv(after_file)

    before_time = find_time_col(before)
    after_time = find_time_col(after)

    if before_time is None or after_time is None:
        print(f"{site}: timestamp column missing")
        continue

    before["time"] = pd.to_datetime(
        before[before_time],
        errors="coerce"
    )

    after["time"] = pd.to_datetime(
        after[after_time],
        errors="coerce"
    )

    before = (
        before
        .dropna(subset=["time"])
        .sort_values("time")
        .drop_duplicates("time")
    )

    after = (
        after
        .dropna(subset=["time"])
        .sort_values("time")
        .drop_duplicates("time")
    )


    # --------------------------------------------------------
    # Site-months that ACTUALLY survived final analysis
    # --------------------------------------------------------
    site_used = used_months[
        used_months["site_name"] == site
    ][["Year", "month"]].copy()

    if site_used.empty:
        continue


    for raw_var, filled_var in mapping.items():

        if raw_var not in before.columns:
            print(f"{site}: {raw_var} missing before MDS")
            continue

        if filled_var not in after.columns:
            print(f"{site}: {filled_var} missing after MDS")
            continue


        merged = before[
            ["time", raw_var]
        ].merge(
            after[["time", filled_var]],
            on="time",
            how="inner"
        )

        merged = merged.sort_values(
            "time"
        ).reset_index(drop=True)

        merged["Year"] = merged["time"].dt.year
        merged["month"] = merged["time"].dt.month


        # ----------------------------------------------------
        # KEEP ONLY EXACT MONTHS THAT ENTER FINAL ANALYSIS
        # ----------------------------------------------------
        merged = merged.merge(
            site_used.assign(used_final=True),
            on=["Year", "month"],
            how="inner"
        )

        if merged.empty:
            continue


        # ----------------------------------------------------
        # DEFINE MDS FILLING
        #
        # Raw missing in blended_gaps2
        # but available after REddyProc = supplied by MDS
        # ----------------------------------------------------
        merged["mds_filled"] = (
            merged[raw_var].isna() &
            merged[filled_var].notna()
        )

        merged["original_available"] = (
            merged[raw_var].notna()
        )

        merged["final_available"] = (
            merged[filled_var].notna()
        )


        # ====================================================
        # SITE-MONTH LEVEL
        # ====================================================
        for (year, month), m in merged.groupby(
            ["Year", "month"]
        ):

            n_total = len(m)

            n_original = int(
                m["original_available"].sum()
            )

            n_mds = int(
                m["mds_filled"].sum()
            )

            n_final = int(
                m["final_available"].sum()
            )

            longest = longest_true_run(
                m,
                "mds_filled"
            )

            site_month_results.append({
                "site": site,
                "Year": int(year),
                "month": int(month),
                "var": raw_var,

                "n_halfhours": n_total,

                "n_original": n_original,
                "n_mds_filled": n_mds,
                "n_final_available": n_final,

                "pct_original": (
                    100 * n_original / n_total
                    if n_total else np.nan
                ),

                "pct_mds": (
                    100 * n_mds / n_total
                    if n_total else np.nan
                ),

                "pct_final_available": (
                    100 * n_final / n_total
                    if n_total else np.nan
                ),

                "longest_gap_days": (
                    longest["days"]
                    if longest is not None
                    else 0
                )
            })


        # ====================================================
        # SITE-YEAR LEVEL
        # ONLY MONTHS ACTUALLY USED
        # ====================================================
        for year, y in merged.groupby("Year"):

            n_total = len(y)

            n_original = int(
                y["original_available"].sum()
            )

            n_mds = int(
                y["mds_filled"].sum()
            )

            n_final = int(
                y["final_available"].sum()
            )

            longest = longest_true_run(
                y,
                "mds_filled"
            )

            site_year_results.append({
                "site": site,
                "Year": int(year),
                "var": raw_var,

                "n_months_used": int(
                    y["month"].nunique()
                ),

                "n_halfhours": n_total,

                "n_original": n_original,
                "n_mds_filled": n_mds,
                "n_final_available": n_final,

                "pct_original": (
                    100 * n_original / n_total
                    if n_total else np.nan
                ),

                "pct_mds": (
                    100 * n_mds / n_total
                    if n_total else np.nan
                ),

                "pct_final_available": (
                    100 * n_final / n_total
                    if n_total else np.nan
                ),

                "longest_gap_days": (
                    longest["days"]
                    if longest is not None
                    else 0
                ),

                "gap_start": (
                    longest["start"]
                    if longest is not None
                    else pd.NaT
                ),

                "gap_end": (
                    longest["end"]
                    if longest is not None
                    else pd.NaT
                )
            })


# ============================================================
# 5. CREATE RESULT TABLES
# ============================================================
sy = pd.DataFrame(site_year_results)
sm = pd.DataFrame(site_month_results)

sy = sy.merge(
    site_meta[
        [
            "site_name",
            "lat",
            "long",
            "ecosystem",
            "coast",
            "is_alaska"
        ]
    ].rename(
        columns={"site_name": "site"}
    ),
    on="site",
    how="left"
)

sm = sm.merge(
    site_meta[
        [
            "site_name",
            "lat",
            "long",
            "ecosystem",
            "coast",
            "is_alaska"
        ]
    ].rename(
        columns={"site_name": "site"}
    ),
    on="site",
    how="left"
)


# ============================================================
# 6. OVERALL MDS CONTRIBUTION
# ============================================================
print("\n" + "=" * 110)
print("OVERALL MDS CONTRIBUTION — FINAL ANALYTICAL MONTHS ONLY")
print("=" * 110)

for var in ["NEE", "LE"]:

    x = sy[sy["var"] == var]

    total_hh = x["n_halfhours"].sum()
    total_mds = x["n_mds_filled"].sum()
    total_original = x["n_original"].sum()

    print(f"\n{var}")

    print(
        f"  Total half-hours represented: "
        f"{total_hh:,}"
    )

    print(
        f"  Original available: "
        f"{100 * total_original / total_hh:.2f}%"
    )

    print(
        f"  Supplied by MDS: "
        f"{100 * total_mds / total_hh:.2f}%"
    )

    print(
        f"  Median MDS % across site-years: "
        f"{x['pct_mds'].median():.2f}%"
    )

    print(
        f"  Mean MDS % across site-years: "
        f"{x['pct_mds'].mean():.2f}%"
    )

    print(
        f"  75th percentile MDS %: "
        f"{x['pct_mds'].quantile(0.75):.2f}%"
    )

    print(
        f"  90th percentile MDS %: "
        f"{x['pct_mds'].quantile(0.90):.2f}%"
    )


# ============================================================
# 7. GAP LENGTHS IN ACTUALLY USED SITE-YEARS
# ============================================================
print("\n" + "=" * 110)
print("LONGEST MDS GAPS IN ACTUALLY USED SITE-YEARS")
print("=" * 110)

for var in ["NEE", "LE"]:

    x = sy[sy["var"] == var].copy()

    print(f"\n{var}")

    print(
        f"  Site-years evaluated: "
        f"{len(x)}"
    )

    for threshold in [7, 14, 30]:

        n = (
            x["longest_gap_days"] > threshold
        ).sum()

        pct = (
            100 * n / len(x)
            if len(x)
            else np.nan
        )

        print(
            f"  Site-years with >{threshold}-day gap: "
            f"{n} ({pct:.1f}%)"
        )

    if not x.empty:

        idx = x["longest_gap_days"].idxmax()
        row = x.loc[idx]

        print(
            f"  Maximum gap: "
            f"{row['longest_gap_days']:.2f} days"
        )

        print(
            f"    Site:  {row['site']}"
        )

        print(
            f"    Year:  {int(row['Year'])}"
        )

        print(
            f"    Start: {row['gap_start']}"
        )

        print(
            f"    End:   {row['gap_end']}"
        )


# ============================================================
# 8. TOP SITE-YEARS BY PERCENT MDS
# ============================================================
print("\n" + "=" * 110)
print("TOP 20 SITE-YEARS BY MDS CONTRIBUTION")
print("=" * 110)

for var in ["NEE", "LE"]:

    print(f"\n{var}")

    cols = [
        "site",
        "Year",
        "ecosystem",
        "coast",
        "n_months_used",
        "pct_original",
        "pct_mds",
        "longest_gap_days"
    ]

    print(
        sy[
            sy["var"] == var
        ]
        .sort_values(
            "pct_mds",
            ascending=False
        )
        .head(20)[cols]
        .round({
            "pct_original": 2,
            "pct_mds": 2,
            "longest_gap_days": 2
        })
        .to_string(index=False)
    )


# ============================================================
# 9. MONTH-LEVEL HEAVY GAP FILLING
# ============================================================
print("\n" + "=" * 110)
print("FINAL SITE-MONTHS WITH LARGE MDS CONTRIBUTION")
print("=" * 110)

for var in ["NEE", "LE"]:

    x = sm[sm["var"] == var]

    print(f"\n{var}")

    for threshold in [10, 25, 50, 75]:

        n = (
            x["pct_mds"] > threshold
        ).sum()

        pct = (
            100 * n / len(x)
            if len(x)
            else np.nan
        )

        print(
            f"  Months with >{threshold}% MDS: "
            f"{n} / {len(x)} "
            f"({pct:.1f}%)"
        )


# ============================================================
# 10. ALASKA-SPECIFIC RESULTS
# ============================================================
print("\n" + "=" * 110)
print("ALASKA — MDS GAP-FILLING DIAGNOSTIC")
print("=" * 110)

ak = sy[sy["is_alaska"] == True].copy()

print(
    f"\nAlaska sites represented: "
    f"{ak['site'].nunique()}"
)

print(
    ", ".join(
        sorted(ak["site"].unique())
    )
)


for var in ["NEE", "LE"]:

    x = ak[ak["var"] == var]

    if x.empty:
        continue

    total_hh = x["n_halfhours"].sum()
    total_mds = x["n_mds_filled"].sum()

    print(f"\n{var}")

    print(
        f"  Site-years: "
        f"{len(x)}"
    )

    print(
        f"  Overall MDS contribution: "
        f"{100 * total_mds / total_hh:.2f}%"
    )

    print(
        f"  Median site-year MDS contribution: "
        f"{x['pct_mds'].median():.2f}%"
    )

    print(
        f"  Mean site-year MDS contribution: "
        f"{x['pct_mds'].mean():.2f}%"
    )

    print(
        f"  Maximum site-year MDS contribution: "
        f"{x['pct_mds'].max():.2f}%"
    )

    print(
        f"  Maximum continuous MDS gap: "
        f"{x['longest_gap_days'].max():.2f} days"
    )

    for threshold in [7, 14, 30]:

        n = (
            x["longest_gap_days"] > threshold
        ).sum()

        pct = (
            100 * n / len(x)
            if len(x)
            else np.nan
        )

        print(
            f"  Site-years with >{threshold}-day gap: "
            f"{n} ({pct:.1f}%)"
        )


# ============================================================
# 11. ALASKA SITE-YEAR DETAIL
# ============================================================
print("\n" + "=" * 110)
print("ALASKA SITE-YEAR DETAIL")
print("=" * 110)

ak_cols = [
    "site",
    "Year",
    "var",
    "ecosystem",
    "n_months_used",
    "pct_original",
    "pct_mds",
    "longest_gap_days",
    "gap_start",
    "gap_end"
]

print(
    ak[
        ak_cols
    ]
    .sort_values(
        ["site", "Year", "var"]
    )
    .round({
        "pct_original": 2,
        "pct_mds": 2,
        "longest_gap_days": 2
    })
    .to_string(index=False)
)


# ============================================================
# 12. OVERALL MOST EXTREME RETAINED GAPS
# ============================================================
print("\n" + "=" * 110)
print("TOP 30 LONGEST MDS GAPS THAT ACTUALLY ENTER THE ANALYSIS")
print("=" * 110)

cols = [
    "site",
    "Year",
    "var",
    "ecosystem",
    "coast",
    "n_months_used",
    "pct_mds",
    "longest_gap_days",
    "gap_start",
    "gap_end"
]

print(
    sy.sort_values(
        "longest_gap_days",
        ascending=False
    )
    .head(30)[cols]
    .round({
        "pct_mds": 2,
        "longest_gap_days": 2
    })
    .to_string(index=False)
)


print("\n" + "=" * 110)
print("DONE")
print("=" * 110)
print("Nothing was saved.")