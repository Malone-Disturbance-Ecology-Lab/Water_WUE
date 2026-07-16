# -*- coding: utf-8 -*-
"""
Diagnostic plot: OLD vs NEW coast classification

Purpose:
- Compare old coast classification (fixed longitude cutoffs) vs updated 
  site-coordinate nearest-coastline classification.
- Uses SAME Q3/Q4 monthly-score filter:
    site_name, Year, month, water_class, lat, long, WUE_tra,
    SPEI_1, SPEI_3, SPEI_6, SPEI_12, SPEI_24, SPEI_36, SPEI_48
- Keeps only Upland, Freshwater, Saline.
- Uses only finite lat/long and finite WUE_tra.
- Expected plotted sites: 64.
- Plot only actual site markers.
- No changed-site rings.
- No site-name labels.
- Big legend with N sites per coast.
- Saves nothing.
"""

from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# =============================================================================
# INPUT
# =============================================================================

DATA_FILE = Path(
    r"M:\Research\WUE_CUE\Water_WUE\data\WUE_CUE_monthly_merged_indices_clean.csv"
)

print("\n" + "=" * 90)
print("DIAGNOSTIC MAP: OLD VS UPDATED COAST CLASSIFICATION")
print("FILTER: SAME Q3/Q4 monthly-score filter")
print("=" * 90)
print(f"Input file: {DATA_FILE}")

if not DATA_FILE.exists():
    raise FileNotFoundError(f"Could not find file:\n{DATA_FILE}")

df = pd.read_csv(DATA_FILE)


# =============================================================================
# Q3/Q4 FILTER SETTINGS
# =============================================================================

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]

SPEI_COLS = [
    "SPEI_1", "SPEI_3", "SPEI_6",
    "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"
]

required_cols = [
    "site_name", "Year", "month", "water_class", "lat", "long", "WUE_tra"
] + SPEI_COLS

missing = [c for c in required_cols if c not in df.columns]

if missing:
    raise ValueError(f"Missing required columns: {missing}")


# =============================================================================
# CLEAN TEXT AND NUMERIC COLUMNS
# =============================================================================

for col in df.select_dtypes(include=["object"]).columns:
    df[col] = df[col].where(
        df[col].isna(),
        df[col].astype(str).str.strip()
    )
    df[col] = df[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

numeric_cols = ["Year", "month", "lat", "long", "WUE_tra"] + SPEI_COLS

for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.replace([np.inf, -np.inf], np.nan)


# =============================================================================
# APPLY SAME Q3/Q4 FILTERS STEP BY STEP
# =============================================================================

df0 = df[required_cols].copy()

print("\nRaw data:")
print("Rows:", len(df))
print("Unique sites before Q3/Q4 filter:", df["site_name"].nunique())

print("\n" + "=" * 90)
print("FILTER STEP COUNTS")
print("=" * 90)

print("Start rows:", len(df0), "sites:", df0["site_name"].nunique())

df1 = df0.dropna(subset=required_cols).copy()
print("After dropna required cols:", len(df1), "sites:", df1["site_name"].nunique())

df2 = df1[df1["site_name"] != ""].copy()
print("After non-empty site_name:", len(df2), "sites:", df2["site_name"].nunique())

df3 = df2[df2["water_class"].isin(ECOSYSTEM_CLASSES)].copy()
print("After water_class in Upland/Freshwater/Saline:", len(df3), "sites:", df3["site_name"].nunique())

df4 = df3[np.isfinite(df3["lat"]) & np.isfinite(df3["long"])].copy()
print("After finite lat/long:", len(df4), "sites:", df4["site_name"].nunique())

df5 = df4[np.isfinite(df4["WUE_tra"])].copy()
print("After finite WUE_tra:", len(df5), "sites:", df5["site_name"].nunique())


# =============================================================================
# SITE TABLE: ONE ROW PER SITE AFTER SAME Q3/Q4 FILTERS
# =============================================================================

site_df = (
    df5
    .groupby("site_name", as_index=False)
    .agg(
        lat=("lat", "mean"),
        long=("long", "mean"),
        water_class=("water_class", "first"),
        n_valid_months=("WUE_tra", "count")
    )
)

print(f"\nSites plotted after SAME Q3/Q4 monthly-score filters: {len(site_df)}")


# =============================================================================
# OLD CLASSIFIER (fixed longitude cutoffs)
# =============================================================================

def classify_old(lat, long):
    """
    Old coast classifier using fixed longitude cutoffs.
    """
    if pd.isna(lat) or pd.isna(long):
        return "Other/Check"

    lat = float(lat)
    lon = float(long)

    if lat > 50:
        return "AK Coast"
    elif lon > -100:
        return "Atlantic Coast"
    elif lon < -120:
        return "Pacific Coast"
    else:
        return "Gulf Coast"


# =============================================================================
# UPDATED CLASSIFIER (site-coordinate nearest-coastline)
# =============================================================================

EARTH_RADIUS_KM = 6371.0088

COAST_POLYLINES = {
    "AK Coast": [
        (60.0, -165.0), (63.0, -165.0), (66.0, -163.5), (69.0, -166.0),
        (71.0, -160.0), (71.5, -156.0), (70.5, -150.0), (69.5, -145.0)
    ],

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
    """
    Minimum distance from point to coastline segment using local projection.
    """
    lat0 = math.radians(lat)

    def proj(lat2, lon2):
        x = EARTH_RADIUS_KM * math.radians(lon2 - lon) * math.cos(lat0)
        y = EARTH_RADIUS_KM * math.radians(lat2 - lat)
        return np.array([x, y], dtype=float)

    p = np.array([0.0, 0.0], dtype=float)
    a = proj(a_lat, a_lon)
    b = proj(b_lat, b_lon)

    ab = b - a
    denom = np.dot(ab, ab)

    if denom == 0:
        return float(np.linalg.norm(p - a))

    t = np.dot(p - a, ab) / denom
    t = max(0.0, min(1.0, t))

    closest = a + t * ab
    return float(np.linalg.norm(p - closest))


def _distance_to_polyline(lat, lon, polyline):
    """
    Minimum distance from site point to coast polyline.
    """
    distances = []

    for i in range(len(polyline) - 1):
        d = _seg_dist(lat, lon, *polyline[i], *polyline[i + 1])
        distances.append(d)

    return min(distances) if distances else float("inf")


def classify_updated(lat, long):
    """
    Updated coast classifier: site-coordinate nearest-coastline.
    """
    if pd.isna(lat) or pd.isna(long):
        return "Other/Check"

    lat = float(lat)
    lon = float(long)

    if lat >= 50:
        return "AK Coast"

    dists = {
        name: _distance_to_polyline(lat, lon, polyline)
        for name, polyline in COAST_POLYLINES.items()
        if name != "AK Coast"
    }

    return min(dists, key=dists.get)


# =============================================================================
# APPLY BOTH METHODS
# =============================================================================

site_df["old_coast"] = site_df.apply(
    lambda row: classify_old(row["lat"], row["long"]),
    axis=1
)

site_df["updated_coast"] = site_df.apply(
    lambda row: classify_updated(row["lat"], row["long"]),
    axis=1
)


# =============================================================================
# PRINT SUMMARY
# =============================================================================

print("\nOLD classification (fixed longitude cutoffs) counts after SAME Q3/Q4 filter:")
print(site_df["old_coast"].value_counts())

print("\nUPDATED classification (nearest-coastline) counts after SAME Q3/Q4 filter:")
print(site_df["updated_coast"].value_counts())

print("\nOld -> Updated crosstab after SAME Q3/Q4 filter:")
print(pd.crosstab(site_df["old_coast"], site_df["updated_coast"], margins=True))

print("\nAny unassigned / Other/Check sites in UPDATED method?")
other_updated = site_df[site_df["updated_coast"] == "Other/Check"]

if other_updated.empty:
    print("No. All 64 filtered sites were assigned to a coast.")
else:
    print(other_updated[["site_name", "lat", "long", "updated_coast"]].to_string(index=False))


# =============================================================================
# PLOT SETTINGS
# =============================================================================

COAST_COLORS = {
    "Atlantic Coast": "#0000FF",   # blue
    "Pacific Coast": "#FF0000",    # red
    "Gulf Coast": "#00AA00",       # green
    "AK Coast": "#8000FF",         # purple
}

COAST_ORDER = [
    "Atlantic Coast",
    "Pacific Coast",
    "Gulf Coast",
    "AK Coast"
]

X_LIM = (-170, -65)
Y_LIM = (23, 73)


def coast_counts_for_label(data, coast_col):
    """
    Count unique sites in each coast.
    """
    return (
        data.groupby(coast_col)["site_name"]
        .nunique()
        .reindex(COAST_ORDER)
        .fillna(0)
        .astype(int)
        .to_dict()
    )


def plot_coast_polylines(ax):
    """
    Draw approximate coast polylines used by the updated classifier.
    """
    for coast, polyline in COAST_POLYLINES.items():
        if coast not in COAST_COLORS:
            continue

        lats = [p[0] for p in polyline]
        lons = [p[1] for p in polyline]

        ax.plot(
            lons,
            lats,
            color=COAST_COLORS[coast],
            linewidth=2.5,
            alpha=0.55,
            zorder=1
        )


def plot_one_panel(ax, data, coast_col, title):
    """
    Plot one panel using only actual site markers.
    No changed-site rings.
    No site labels.
    """
    plot_coast_polylines(ax)

    counts = coast_counts_for_label(data, coast_col)

    for coast in COAST_ORDER:
        sub = data[data[coast_col] == coast].copy()

        if sub.empty:
            continue

        color = COAST_COLORS[coast]

        ax.scatter(
            sub["long"],
            sub["lat"],
            s=105,
            facecolor=color,
            edgecolor=color,
            linewidth=1.0,
            alpha=0.95,
            zorder=3
        )

    ax.set_title(
        f"{title}\n"
        f"Atlantic={counts['Atlantic Coast']}, "
        f"Pacific={counts['Pacific Coast']}, "
        f"Gulf={counts['Gulf Coast']}, "
        f"AK={counts['AK Coast']}",
        fontsize=14,
        fontweight="bold"
    )

    ax.set_xlabel("Longitude", fontsize=13)
    ax.set_ylabel("Latitude", fontsize=13)
    ax.set_xlim(X_LIM)
    ax.set_ylim(Y_LIM)
    ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="box")


# =============================================================================
# SIDE-BY-SIDE FIGURE
# =============================================================================

fig, axes = plt.subplots(
    1,
    2,
    figsize=(18, 8.5),
    sharex=True,
    sharey=True
)

plot_one_panel(
    axes[0],
    site_df,
    "old_coast",
    "OLD workflow: fixed longitude cutoffs"
)

plot_one_panel(
    axes[1],
    site_df,
    "updated_coast",
    "UPDATED workflow: site-coordinate nearest-coastline classification"
)


# =============================================================================
# BIG COAST-ONLY LEGEND
# =============================================================================

updated_counts = coast_counts_for_label(site_df, "updated_coast")

legend_handles = [
    Line2D(
        [0], [0],
        marker="o",
        linestyle="None",
        label=f"{coast} (N={updated_counts[coast]})",
        markerfacecolor=COAST_COLORS[coast],
        markeredgecolor=COAST_COLORS[coast],
        markersize=18
    )
    for coast in COAST_ORDER
]

fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=4,
    frameon=True,
    fontsize=16,
    title="Coast regions, updated workflow",
    title_fontsize=18,
    borderpad=1.3,
    labelspacing=1.1,
    handletextpad=0.9
)

fig.suptitle(
    "Coast Classification Comparison for WUE Sites",
    fontsize=19,
    fontweight="bold"
)

plt.tight_layout(rect=[0, 0.14, 1, 0.91])
plt.show()