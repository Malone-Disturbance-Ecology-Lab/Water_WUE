# ============================================================
# 13. CREATE SUPPLEMENTARY MDS GAP-DURATION TABLE
# ============================================================

OUTPUT_DIR = r"M:\Research\WUE_CUE\WUE_manuscript_version6\methods"

os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_CSV = os.path.join(
    OUTPUT_DIR,
    "Table_S_MDS_gap_duration_by_site.csv"
)

# ------------------------------------------------------------
# FULL COAST ASSIGNMENT
# Same classifier used in the coverage figure
# ------------------------------------------------------------
study = study.copy()

study["coast_region"] = study.apply(
    lambda row: assign_coast_region(
        row["lat"],
        row["long"]
    ),
    axis=1
)

coast_name_map = {
    "AK Coast": "Alaska",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf",
    "Atlantic Coast": "Atlantic"
}

study["coast_final"] = (
    study["coast_region"]
    .map(coast_name_map)
)

# ------------------------------------------------------------
# FINAL ANALYTICAL COVERAGE ONLY
# ------------------------------------------------------------
site_coverage = (
    study.groupby("site_name")
    .agg(
        Coast=("coast_final", "first"),
        Years_Used=("Year", "nunique"),
        Site_Months_Used=("Year", "size")
    )
    .reset_index()
    .rename(columns={"site_name": "site"})
)

# ------------------------------------------------------------
# SUMMARIZE CONTINUOUS MDS GAP DURATION BY SITE
# ------------------------------------------------------------
def summarize_gap_duration(sy, var):

    x = sy[
        sy["var"] == var
    ].copy()

    summary = (
        x.groupby("site")
        .agg(
            Years_GT14=(
                "longest_gap_days",
                lambda z: int((z > 14).sum())
            ),

            Years_GT30=(
                "longest_gap_days",
                lambda z: int((z > 30).sum())
            ),

            Max_Gap_Days=(
                "longest_gap_days",
                "max"
            )
        )
        .reset_index()
    )

    summary = summary.rename(
        columns={
            "Years_GT14": f"{var}_Years_GT14",
            "Years_GT30": f"{var}_Years_GT30",
            "Max_Gap_Days": f"{var}_Max_Gap_Days"
        }
    )

    return summary


# ------------------------------------------------------------
# NEE AND LE
# ------------------------------------------------------------
nee_summary = summarize_gap_duration(
    sy,
    "NEE"
)

le_summary = summarize_gap_duration(
    sy,
    "LE"
)

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------
table_s = (
    site_coverage
    .merge(
        nee_summary,
        on="site",
        how="left"
    )
    .merge(
        le_summary,
        on="site",
        how="left"
    )
)

# ------------------------------------------------------------
# ROUND GAP LENGTHS
# ------------------------------------------------------------
table_s["NEE_Max_Gap_Days"] = (
    table_s["NEE_Max_Gap_Days"]
    .round(1)
)

table_s["LE_Max_Gap_Days"] = (
    table_s["LE_Max_Gap_Days"]
    .round(1)
)

# ------------------------------------------------------------
# KEEP ONLY COLUMNS NEEDED FOR SUPPLEMENT
# ------------------------------------------------------------
table_s = table_s[
    [
        "site",
        "Coast",
        "Years_Used",
        "Site_Months_Used",

        "NEE_Years_GT14",
        "NEE_Years_GT30",
        "NEE_Max_Gap_Days",

        "LE_Years_GT14",
        "LE_Years_GT30",
        "LE_Max_Gap_Days"
    ]
]

# ------------------------------------------------------------
# MANUSCRIPT-READY COLUMN NAMES
# ------------------------------------------------------------
table_s = table_s.rename(
    columns={
        "site":
            "AmeriFlux Site",

        "Years_Used":
            "Years Used",

        "Site_Months_Used":
            "Site-Months Used",

        "NEE_Years_GT14":
            "NEE Years >14-d Gap",

        "NEE_Years_GT30":
            "NEE Years >30-d Gap",

        "NEE_Max_Gap_Days":
            "Max NEE Gap (days)",

        "LE_Years_GT14":
            "LE Years >14-d Gap",

        "LE_Years_GT30":
            "LE Years >30-d Gap",

        "LE_Max_Gap_Days":
            "Max LE Gap (days)"
    }
)

# ------------------------------------------------------------
# SORT BY COAST
# ------------------------------------------------------------
coast_order = {
    "Alaska": 0,
    "Pacific": 1,
    "Gulf": 2,
    "Atlantic": 3
}

table_s["coast_order"] = (
    table_s["Coast"]
    .map(coast_order)
)

table_s = (
    table_s
    .sort_values(
        [
            "coast_order",
            "AmeriFlux Site"
        ]
    )
    .drop(
        columns="coast_order"
    )
    .reset_index(drop=True)
)

# ------------------------------------------------------------
# CHECK
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("SUPPLEMENTARY MDS GAP-DURATION TABLE")
print("=" * 110)

print(f"\nTotal sites: {len(table_s)}")

print("\nSites per coast:")
print(table_s["Coast"].value_counts())

print("\n")
print(table_s.to_string(index=False))

# ------------------------------------------------------------
# SAVE CSV
# ------------------------------------------------------------
table_s.to_csv(
    OUTPUT_CSV,
    index=False
)

print(f"\nSaved to:\n{OUTPUT_CSV}")