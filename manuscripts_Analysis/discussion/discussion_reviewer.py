# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# ============================================================
# DIAGNOSTIC: PACIFIC / CALIFORNIA / DELTA REPRESENTATION
# Nothing is saved
# ============================================================

import pandas as pd
import numpy as np
import math

MONTHLY_DATA = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\WUE_CUE_monthly_merged_indices_clean.csv"
)

# ------------------------------------------------------------
# SAME COAST CLASSIFIER USED IN THE STUDY
# ------------------------------------------------------------
EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "Pacific Coast": [
        (32.6, -117.2), (33.6, -118.2), (34.4, -119.7),
        (35.4, -120.9), (36.6, -121.9), (37.8, -122.5),
        (39.0, -123.7), (41.0, -124.2), (43.5, -124.2),
        (45.5, -123.9), (47.0, -124.1), (48.8, -124.7)
    ],
    "Gulf Coast": [
        (25.8, -97.2), (28.0, -96.8), (29.3, -94.8),
        (29.5, -93.0), (29.2, -91.5), (29.3, -90.0),
        (29.5, -88.8), (30.2, -87.7), (30.1, -86.2),
        (29.9, -85.3), (29.7, -84.3), (28.8, -83.0),
        (27.8, -82.8), (26.6, -82.2), (25.9, -81.8),
        (25.3, -81.1), (25.0, -80.8)
    ],
    "Atlantic Coast": [
        (25.1, -80.3), (26.0, -80.1), (27.0, -80.1),
        (28.4, -80.6), (29.0, -80.9), (29.9, -81.3),
        (30.4, -81.4), (31.2, -81.3), (32.0, -80.8),
        (33.0, -79.5), (34.5, -77.8), (35.7, -75.6),
        (36.8, -75.9), (38.5, -75.0), (39.5, -74.3),
        (40.5, -73.9), (41.3, -72.0), (42.4, -70.8),
        (43.5, -70.2)
    ]
}

def _seg_dist(lat, lon, a_lat, a_lon, b_lat, b_lon):
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
    return min(
        _seg_dist(lat, lon, *polyline[i], *polyline[i + 1])
        for i in range(len(polyline) - 1)
    )

def assign_coast_region(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return np.nan

    lat = float(lat)
    lon = float(lon)

    if lat > 50:
        return "AK Coast"

    dists = {
        coast: _distance_to_polyline(lat, lon, polyline)
        for coast, polyline in COAST_POLYLINES.items()
    }

    return min(dists, key=dists.get)


# ------------------------------------------------------------
# LOAD AND APPLY EXACT STUDY-AREA FILTERS
# ------------------------------------------------------------
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

study["coast_region"] = study.apply(
    lambda row: assign_coast_region(
        row["lat"],
        row["long"]
    ),
    axis=1
)

# ------------------------------------------------------------
# PACIFIC COAST ONLY
# ------------------------------------------------------------
pac = study[
    study["coast_region"] == "Pacific Coast"
].copy()

# California border is 42 N.
pac["CA_status"] = np.where(
    pac["lat"] < 42,
    "California",
    "Oregon/Washington"
)

# ------------------------------------------------------------
# DELTA SITE DEFINITIONS
# ------------------------------------------------------------
# Sites clearly located in the Sacramento-San Joaquin Delta
DELTA_STRICT = {
    "US-Dmg",
    "US-Myb",
    "US-Snd",
    "US-Sne",
    "US-Tw1",
    "US-Tw2",
    "US-Tw3",
    "US-Tw4",
    "US-Tw5",
    "US-Twt",
}

# US-CGG is near the Suisun/Bay-Delta region.
# Print results both ways rather than assuming how it should
# be classified.
DELTA_EXPANDED = DELTA_STRICT | {"US-CGG"}

pac["Delta_strict"] = np.where(
    pac["site_name"].isin(DELTA_STRICT),
    "Delta",
    "Non-Delta"
)

pac["Delta_expanded"] = np.where(
    pac["site_name"].isin(DELTA_EXPANDED),
    "Delta/Bay-Delta",
    "Non-Delta"
)

# ------------------------------------------------------------
# SITE-LEVEL SUMMARY
# ------------------------------------------------------------
site_summary = (
    pac.groupby("site_name")
    .agg(
        latitude=("lat", "first"),
        longitude=("long", "first"),
        ecosystem=("water_class", "first"),
        first_year=("Year", "min"),
        last_year=("Year", "max"),
        n_years=("Year", "nunique"),
        n_site_months=("Year", "size"),
        CA_status=("CA_status", "first"),
        Delta_strict=("Delta_strict", "first"),
        Delta_expanded=("Delta_expanded", "first")
    )
    .reset_index()
)

year_lists = (
    pac.groupby("site_name")["Year"]
    .apply(
        lambda x: ", ".join(
            map(str, sorted(x.astype(int).unique()))
        )
    )
    .reset_index(name="years_present")
)

site_summary = site_summary.merge(
    year_lists,
    on="site_name",
    how="left"
)

site_summary["pct_of_Pacific_months"] = (
    100 *
    site_summary["n_site_months"] /
    len(pac)
)

ca_total = (
    site_summary.loc[
        site_summary["CA_status"] == "California",
        "n_site_months"
    ].sum()
)

site_summary["pct_of_CA_months"] = np.where(
    site_summary["CA_status"] == "California",
    100 * site_summary["n_site_months"] / ca_total,
    np.nan
)

print("\n" + "=" * 110)
print("PACIFIC COAST SITES IN FINAL ANALYTICAL DATASET")
print("=" * 110)

print(
    site_summary.sort_values(
        ["CA_status", "Delta_strict", "site_name"]
    ).to_string(
        index=False,
        float_format=lambda x: f"{x:.1f}"
    )
)


# ------------------------------------------------------------
# OVERALL PACIFIC / CALIFORNIA CONTRIBUTION
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("PACIFIC COAST REPRESENTATION")
print("=" * 110)

print(f"Pacific sites: {pac['site_name'].nunique()}")
print(f"Pacific site-months: {len(pac)}")

ca = pac[pac["CA_status"] == "California"].copy()

print(f"\nCalifornia Pacific sites: {ca['site_name'].nunique()}")
print(f"California Pacific site-months: {len(ca)}")
print(
    f"California contribution to all Pacific site-months: "
    f"{100 * len(ca) / len(pac):.1f}%"
)


# ------------------------------------------------------------
# STRICT DELTA DEFINITION
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("STRICT DELTA DEFINITION")
print("=" * 110)

strict_summary = (
    ca.groupby("Delta_strict")
    .agg(
        n_sites=("site_name", "nunique"),
        n_site_months=("site_name", "size"),
        first_year=("Year", "min"),
        last_year=("Year", "max")
    )
)

strict_summary["pct_CA_site_months"] = (
    100 *
    strict_summary["n_site_months"] /
    len(ca)
)

strict_summary["pct_Pacific_site_months"] = (
    100 *
    strict_summary["n_site_months"] /
    len(pac)
)

print(strict_summary.to_string())


# ------------------------------------------------------------
# EXPANDED DELTA / BAY-DELTA DEFINITION
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("EXPANDED DELTA / BAY-DELTA DEFINITION")
print("Includes US-CGG")
print("=" * 110)

expanded_summary = (
    ca.groupby("Delta_expanded")
    .agg(
        n_sites=("site_name", "nunique"),
        n_site_months=("site_name", "size"),
        first_year=("Year", "min"),
        last_year=("Year", "max")
    )
)

expanded_summary["pct_CA_site_months"] = (
    100 *
    expanded_summary["n_site_months"] /
    len(ca)
)

expanded_summary["pct_Pacific_site_months"] = (
    100 *
    expanded_summary["n_site_months"] /
    len(pac)
)

print(expanded_summary.to_string())


# ------------------------------------------------------------
# ECOSYSTEM COMPOSITION: DELTA VS NON-DELTA
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("ECOSYSTEM COMPOSITION OF CALIFORNIA SITE-MONTHS")
print("=" * 110)

eco_summary = pd.crosstab(
    ca["Delta_strict"],
    ca["water_class"],
    margins=True
)

print(eco_summary.to_string())


# ------------------------------------------------------------
# YEAR-BY-YEAR CONTRIBUTION
# ------------------------------------------------------------
print("\n" + "=" * 110)
print("YEAR-BY-YEAR CALIFORNIA CONTRIBUTION")
print("=" * 110)

year_summary = (
    ca.groupby(["Year", "Delta_strict"])
    .agg(
        n_sites=("site_name", "nunique"),
        n_site_months=("site_name", "size")
    )
    .reset_index()
)

print(
    year_summary.to_string(index=False)
)


# ------------------------------------------------------------
# OPTIONAL: SPEI-3 ELIGIBLE OBSERVATIONS
# This is closer to the drought analyses if SPEI_3 exists.
# ------------------------------------------------------------
spei_candidates = [
    "SPEI_3",
    "SPEI3",
    "spei_3"
]

spei_col = next(
    (c for c in spei_candidates if c in ca.columns),
    None
)

if spei_col is not None:

    ca_spei = ca[
        np.isfinite(
            pd.to_numeric(
                ca[spei_col],
                errors="coerce"
            )
        )
    ].copy()

    print("\n" + "=" * 110)
    print(
        f"CALIFORNIA REPRESENTATION AMONG "
        f"{spei_col}-ELIGIBLE OBSERVATIONS"
    )
    print("=" * 110)

    spei_summary = (
        ca_spei.groupby("Delta_strict")
        .agg(
            n_sites=("site_name", "nunique"),
            n_site_months=("site_name", "size")
        )
    )

    spei_summary["pct_site_months"] = (
        100 *
        spei_summary["n_site_months"] /
        len(ca_spei)
    )

    print(spei_summary.to_string())

else:
    print(
        "\nNo SPEI-3 column detected under expected names; "
        "base analytical representation was still calculated."
    )


# ------------------------------------------------------------
# SIMPLE FINAL DIAGNOSTIC
# ------------------------------------------------------------
delta_months = (
    ca["site_name"].isin(DELTA_STRICT)
).sum()

non_delta_months = len(ca) - delta_months

print("\n" + "=" * 110)
print("BOTTOM LINE")
print("=" * 110)

print(
    f"Strict Delta sites contribute "
    f"{delta_months:,} of {len(ca):,} California site-months "
    f"({100 * delta_months / len(ca):.1f}%)."
)

print(
    f"Non-Delta California sites contribute "
    f"{non_delta_months:,} of {len(ca):,} California site-months "
    f"({100 * non_delta_months / len(ca):.1f}%)."
)

print("\nNothing was saved.")