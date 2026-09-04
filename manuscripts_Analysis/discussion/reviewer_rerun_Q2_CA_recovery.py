# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 13:27:23 2026

@author: ammar
"""

import os
import pandas as pd
import numpy as np

print("=" * 90)
print("Q2 PACIFIC POST-DROUGHT RECOVERY DIAGNOSTIC")
print("=" * 90)

# --------------------------------------------------------------------------
# Existing Q2 outputs
# --------------------------------------------------------------------------

base_dir = (
    r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2"
    r"\Q2_WUE_performance_outputs"
)

coast_dir = os.path.join(
    base_dir,
    "Q_performace metric_T_ET_relationship"
)

site_coast_file = os.path.join(
    coast_dir,
    "Q2_performance_by_coast_site_level_data.csv"
)

recovery_events_file = os.path.join(
    base_dir,
    "Q2_site_recovery_events.csv"
)

coast_tests_file = os.path.join(
    coast_dir,
    "Q2_performance_by_coast_tests.csv"
)

# --------------------------------------------------------------------------
# Load existing outputs
# --------------------------------------------------------------------------

sites = pd.read_csv(site_coast_file)
events = pd.read_csv(recovery_events_file)
coast_tests = pd.read_csv(coast_tests_file)

# --------------------------------------------------------------------------
# Site-level recovery
# --------------------------------------------------------------------------

recovery_sites = sites[
    sites["mean_recovery"].notna()
].copy()

print("\nSITE-LEVEL RECOVERY BY COAST")
print("-" * 90)

site_summary = (
    recovery_sites
    .groupby("coast_region", observed=True)
    .agg(
        n_sites=("site_name", "nunique"),
        mean_recovery=("mean_recovery", "mean"),
        median_recovery=("mean_recovery", "median"),
        min_recovery=("mean_recovery", "min"),
        max_recovery=("mean_recovery", "max")
    )
    .reset_index()
)

print(site_summary.to_string(index=False))

# --------------------------------------------------------------------------
# Pacific site-level recovery
# --------------------------------------------------------------------------

pac_sites = recovery_sites[
    recovery_sites["coast_region"] == "Pacific Coast"
].copy()

print("\n" + "=" * 90)
print("PACIFIC SITE-LEVEL RECOVERY")
print("=" * 90)

print(f"Recovery sites: {pac_sites['site_name'].nunique()}")

if len(pac_sites) > 0:

    vals = pac_sites["mean_recovery"].dropna()

    print(f"Mean recovery:   {vals.mean():.3f}")
    print(f"Median recovery: {vals.median():.3f}")
    print(f"Min:             {vals.min():.3f}")
    print(f"Max:             {vals.max():.3f}")

    q25 = vals.quantile(0.25)
    q75 = vals.quantile(0.75)

    print(f"IQR:             {q25:.3f} – {q75:.3f}")

    print(
        f"Sites recovery >= 1.0: "
        f"{(vals >= 1.0).sum()}/{len(vals)} "
        f"({100 * (vals >= 1.0).mean():.1f}%)"
    )

    print(
        f"Sites within 0.9–1.1: "
        f"{((vals >= 0.9) & (vals <= 1.1)).sum()}/{len(vals)} "
        f"({100 * ((vals >= 0.9) & (vals <= 1.1)).mean():.1f}%)"
    )

    print("\nPacific sites:")
    print(
        pac_sites[
            [
                "site_name",
                "water_class",
                "mean_recovery"
            ]
        ]
        .sort_values("mean_recovery")
        .to_string(index=False)
    )

# --------------------------------------------------------------------------
# Event-level recovery
# Add coast from existing site-level coast classification
# --------------------------------------------------------------------------

coast_lookup = (
    sites[
        ["site_name", "coast_region"]
    ]
    .drop_duplicates("site_name")
)

events = events.merge(
    coast_lookup,
    on="site_name",
    how="left"
)

pac_events = events[
    events["coast_region"] == "Pacific Coast"
].copy()

print("\n" + "=" * 90)
print("PACIFIC EVENT-LEVEL RECOVERY")
print("=" * 90)

print(f"Drought events: {len(pac_events)}")
print(f"Sites represented: {pac_events['site_name'].nunique()}")

if len(pac_events) > 0:

    vals = pac_events["recovery"].dropna()

    print(f"Mean recovery:   {vals.mean():.3f}")
    print(f"Median recovery: {vals.median():.3f}")
    print(f"Min:             {vals.min():.3f}")
    print(f"Max:             {vals.max():.3f}")

    print(
        f"Events recovery >= 1.0: "
        f"{(vals >= 1.0).sum()}/{len(vals)} "
        f"({100 * (vals >= 1.0).mean():.1f}%)"
    )

    print(
        f"Events within 0.9–1.1: "
        f"{((vals >= 0.9) & (vals <= 1.1)).sum()}/{len(vals)} "
        f"({100 * ((vals >= 0.9) & (vals <= 1.1)).mean():.1f}%)"
    )

# --------------------------------------------------------------------------
# Existing coast statistical test for recovery
# --------------------------------------------------------------------------

print("\n" + "=" * 90)
print("EXISTING COAST-REGION RECOVERY TEST")
print("=" * 90)

test = coast_tests[
    coast_tests["metric"] == "mean_recovery"
]

if len(test) > 0:

    row = test.iloc[0]

    print(f"Atlantic n: {row['n_atlantic']}")
    print(f"Pacific n:  {row['n_pacific']}")
    print(f"Gulf n:     {row['n_gulf']}")
    print(f"Alaska n:   {row['n_ak']}")

    print(
        f"Kruskal-Wallis raw p: "
        f"{row['kw_p_all_groups']:.4f}"
    )

    print(
        f"Kruskal-Wallis FDR p: "
        f"{row['kw_p_all_groups_FDR']:.4f}"
    )

print("\nNothing was saved.")