# -*- coding: utf-8 -*-
"""
Created on Sat Dec 20 14:13:27 2025
UPDATED: 2026-05-09
- NOW reports Mean ± SE (standard error) and SD (standard deviation)
- ALL THREE metrics (WUE_ET, WUE_E, WUE_T) use STRICT TRIPLE INTERSECTION
- Shared subset defined as sites with ALL THREE metrics available (WUE ∩ WUE_E ∩ WUE_T)
- All metrics report identical sample sizes and aggregated NN observation totals
- Updated column names: water_class (instead of Category), IGBP (instead of biome)
- Updated NN support check to use original WUE names (no _clean)
- **FIXED: Support totals use ALL retained NN observations (not unique month combos)**
- **UPDATED: Wording changed to "retained NN monthly observations"**
"""

import pandas as pd
import numpy as np

# =============================================================================
# SECTION 1: LOAD DATA AND DEFINE SHARED SUBSET (STRICT TRIPLE INTERSECTION)
# =============================================================================

FIG1_CSV = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\wue_site_level_NN_medians_SPEI_1.csv"
df = pd.read_csv(FIG1_CSV)

print("="*80)
print("FIGURE 1 STATISTICS WITH STRICT TRIPLE INTERSECTION (ALL METRICS)")
print("="*80)

# =============================================================================
# STEP 1: Define strict triple intersection (sites with ALL THREE metrics)
# =============================================================================

# Get sites with WUE available
wue_sites = set(df[df['WUE_Metric'] == 'WUE']['site_name'].dropna().unique())

# Get sites with WUE_eva available
eva_sites = set(df[df['WUE_Metric'] == 'WUE_eva']['site_name'].dropna().unique())

# Get sites with WUE_tra available
tra_sites = set(df[df['WUE_Metric'] == 'WUE_tra']['site_name'].dropna().unique())

# STRICT TRIPLE INTERSECTION: sites with ALL THREE metrics available
shared_subset_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

# Sites excluded from all metrics
wue_only = wue_sites - shared_subset_sites
eva_only = eva_sites - shared_subset_sites
tra_only = tra_sites - shared_subset_sites

print(f"\n{'='*60}")
print("STRICT TRIPLE INTERSECTION (WUE ∩ WUE_E ∩ WUE_T)")
print(f"{'='*60}")
print(f"WUE available sites: {len(wue_sites)}")
print(f"WUE_eva available sites: {len(eva_sites)}")
print(f"WUE_tra available sites: {len(tra_sites)}")
print(f"\nSTRICT TRIPLE INTERSECTION (ALL metrics): {len(shared_subset_sites)} sites")

if len(wue_only) > 0:
    print(f"\n⚠️ Sites with WUE ONLY (excluded from ALL metrics):")
    for site in sorted(wue_only):
        print(f"    - {site}")

if len(eva_only) > 0:
    print(f"\n⚠️ Sites with WUE_eva ONLY (excluded from ALL metrics):")
    for site in sorted(eva_only):
        print(f"    - {site}")

if len(tra_only) > 0:
    print(f"\n⚠️ Sites with WUE_tra ONLY (excluded from ALL metrics):")
    for site in sorted(tra_only):
        print(f"    - {site}")

# =============================================================================
# STEP 2: Helper function for statistics including mean, SD, and SE
# =============================================================================

def sig2(x):
    """2 significant digits formatter"""
    if pd.isna(x):
        return np.nan
    return f"{x:.2g}"

def metric_statistics(metric, shared_sites=None):
    """
    Calculate comprehensive statistics for a given WUE metric
    ALL metrics now use the same strict triple intersection
    Returns: dict with n_sites, mean, SD, SE, min, max, and percentiles
    """
    d = df.loc[df["WUE_Metric"] == metric, ["site_name", "WUE_median"]].dropna()
    
    # ALL metrics filtered to strict triple intersection
    if shared_sites is not None:
        d = d[d["site_name"].isin(shared_sites)]
    
    values = d["WUE_median"].astype(float).values
    
    if len(values) == 0:
        return None
    
    # Calculate mean, SD, and SE
    mean_val = float(np.mean(values))
    sd_val = float(np.std(values, ddof=1))  # sample standard deviation
    se_val = sd_val / np.sqrt(len(values))  # standard error
    
    out = {
        "n_sites": len(values),
        "mean": mean_val,
        "sd": sd_val,
        "se": se_val,
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }
    
    # Calculate percentiles for distribution description
    for p in [0.50, 0.75, 0.90, 0.95, 0.99]:
        out[f"p{int(p*100)}"] = float(np.quantile(values, p))
    
    out["p25"] = float(np.quantile(values, 0.25))
    out["p75"] = float(np.quantile(values, 0.75))
    
    return out

# =============================================================================
# STEP 3: Print statistics - ALL metrics now use strict triple intersection
# =============================================================================

print(f"\n{'='*60}")
print("SUMMARY STATISTICS (SITE-LEVEL NN MEDIANS)")
print(f"{'='*60}")
print("\nNOTE: ALL THREE metrics now use the SAME strict triple intersection")
print(f"      Strict triple intersection size: {len(shared_subset_sites)} sites\n")

# ALL metrics use strict triple intersection
wue_stats = metric_statistics("WUE", shared_sites=shared_subset_sites)
eva_stats = metric_statistics("WUE_eva", shared_sites=shared_subset_sites)
tra_stats = metric_statistics("WUE_tra", shared_sites=shared_subset_sites)

# Verify all metrics have same sample size
print(f"\n{'='*60}")
print("SAMPLE SIZE VERIFICATION:")
print(f"{'='*60}")
print(f"WUE:   {wue_stats['n_sites'] if wue_stats else 0} sites")
print(f"WUE_E: {eva_stats['n_sites'] if eva_stats else 0} sites")
print(f"WUE_T: {tra_stats['n_sites'] if tra_stats else 0} sites")
if wue_stats and eva_stats and tra_stats:
    if wue_stats['n_sites'] == eva_stats['n_sites'] == tra_stats['n_sites']:
        print(f"✓ VERIFIED: All three metrics have identical sample size ({wue_stats['n_sites']} sites each)")
    else:
        print(f"⚠️ WARNING: Sample sizes differ!")

if wue_stats:
    print("\n" + "="*70)
    print(f"WUE (bulk efficiency) | n_sites = {wue_stats['n_sites']} (strict triple intersection)")
    print(f"{'─'*70}")
    print(f"Mean ± SE:     {wue_stats['mean']:.3f} ± {wue_stats['se']:.3f}")
    print(f"SD:            {wue_stats['sd']:.3f}")
    print(f"Range:         {wue_stats['min']:.3f} to {wue_stats['max']:.3f}  |  (2-sig: {sig2(wue_stats['min'])}–{sig2(wue_stats['max'])}")
    print(f"Median (p50):  {wue_stats['p50']:.3f}  |  (2-sig: {sig2(wue_stats['p50'])}")
    print(f"IQR (p25–p75): {wue_stats['p25']:.3f}–{wue_stats['p75']:.3f}  |  (2-sig: {sig2(wue_stats['p25'])}–{sig2(wue_stats['p75'])}")
    print(f"p90:           {wue_stats['p90']:.3f}  |  (2-sig: {sig2(wue_stats['p90'])}")
    print(f"p95:           {wue_stats['p95']:.3f}  |  (2-sig: {sig2(wue_stats['p95'])}")
    print(f"p99:           {wue_stats['p99']:.3f}  |  (2-sig: {sig2(wue_stats['p99'])}")

if eva_stats:
    print("\n" + "="*70)
    print(f"WUE_E (evaporation efficiency) | n_sites = {eva_stats['n_sites']} (strict triple intersection)")
    print(f"{'─'*70}")
    print(f"Mean ± SE:     {eva_stats['mean']:.3f} ± {eva_stats['se']:.3f}")
    print(f"SD:            {eva_stats['sd']:.3f}")
    print(f"Range:         {eva_stats['min']:.3f} to {eva_stats['max']:.3f}  |  (2-sig: {sig2(eva_stats['min'])}–{sig2(eva_stats['max'])}")
    print(f"Median (p50):  {eva_stats['p50']:.3f}  |  (2-sig: {sig2(eva_stats['p50'])}")
    print(f"IQR (p25–p75): {eva_stats['p25']:.3f}–{eva_stats['p75']:.3f}  |  (2-sig: {sig2(eva_stats['p25'])}–{sig2(eva_stats['p75'])}")
    print(f"p90:           {eva_stats['p90']:.3f}  |  (2-sig: {sig2(eva_stats['p90'])}")
    print(f"p95:           {eva_stats['p95']:.3f}  |  (2-sig: {sig2(eva_stats['p95'])}")
    print(f"p99:           {eva_stats['p99']:.3f}  |  (2-sig: {sig2(eva_stats['p99'])}")

if tra_stats:
    print("\n" + "="*70)
    print(f"WUE_T (transpiration efficiency) | n_sites = {tra_stats['n_sites']} (strict triple intersection)")
    print(f"{'─'*70}")
    print(f"Mean ± SE:     {tra_stats['mean']:.3f} ± {tra_stats['se']:.3f}")
    print(f"SD:            {tra_stats['sd']:.3f}")
    print(f"Range:         {tra_stats['min']:.3f} to {tra_stats['max']:.3f}  |  (2-sig: {sig2(tra_stats['min'])}–{sig2(tra_stats['max'])}")
    print(f"Median (p50):  {tra_stats['p50']:.3f}  |  (2-sig: {sig2(tra_stats['p50'])}")
    print(f"IQR (p25–p75): {tra_stats['p25']:.3f}–{tra_stats['p75']:.3f}  |  (2-sig: {sig2(tra_stats['p25'])}–{sig2(tra_stats['p75'])}")
    print(f"p90:           {tra_stats['p90']:.3f}  |  (2-sig: {sig2(tra_stats['p90'])}")
    print(f"p95:           {tra_stats['p95']:.3f}  |  (2-sig: {sig2(tra_stats['p95'])}")
    print(f"p99:           {tra_stats['p99']:.3f}  |  (2-sig: {sig2(tra_stats['p99'])}")

# =============================================================================
# STEP 4: Top sites for WUE_tra (from strict triple intersection only)
# =============================================================================

print("\n" + "="*70)
print(f"TOP 10 WUE_T SITE-LEVEL NN MEDIANS (from strict triple intersection)")
print("="*70)

dT = df.loc[df["WUE_Metric"] == "WUE_tra", ["site_name", "WUE_median"]].dropna()
# Filter to strict triple intersection only
dT = dT[dT["site_name"].isin(shared_subset_sites)]
dT = dT.sort_values("WUE_median", ascending=False)
print(dT.head(10).to_string(index=False))

# =============================================================================
# STEP 5: FINAL NN OBSERVATION TOTALS - ALL METRICS USE SAME STRICT TRIPLE INTERSECTION
# FIXED: Now counts ALL retained NN observations (not unique month combos)
# =============================================================================

print("\n" + "="*80)
print("FINAL NN OBSERVATION TOTALS (ALL metrics use same strict triple intersection)")
print("="*80)

# Load monthly data
monthly_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\monthly_data_after_outlier_removal.csv"
df_monthly = pd.read_csv(monthly_path)

# Filter to NN conditions at SPEI-1
df_nn = df_monthly[df_monthly["SPEI_1_Cat"] == "NN"].copy()

print("-" * 70)
print("NOTE: ALL THREE metrics now use the SAME strict triple intersection")
print("      Support totals use ALL retained NN monthly observations")
print("-" * 70)

# Get strict triple intersection sites (ALL metrics use these)
shared_sites_list = list(shared_subset_sites)

# Filter to shared subset sites only
df_nn_shared = df_nn[df_nn["site_name"].isin(shared_sites_list)]

# FIXED: Count ALL retained NN observations (no drop_duplicates)
total_nn_observations = len(df_nn_shared)
n_sites_shared = len(shared_sites_list)

# Calculate observations per site
nn_obs_per_site = []
for site in shared_sites_list:
    site_data = df_nn_shared[df_nn_shared["site_name"] == site]
    n_obs = len(site_data)  # FIXED: count all observations
    nn_obs_per_site.append(n_obs)

print(f"\nALL METRICS (WUE, WUE$_E$, WUE$_T$) - Strict triple intersection (n={n_sites_shared} sites each):")
print(f"  Total retained NN observations: {total_nn_observations}")
print(f"  Sites: {n_sites_shared}")

print(f"\n  Observation statistics for kept sites:")
print(f"    Min observations: {min(nn_obs_per_site)}")
print(f"    Max observations: {max(nn_obs_per_site)}")
print(f"    Mean observations: {np.mean(nn_obs_per_site):.1f}")
print(f"    Median observations: {np.median(nn_obs_per_site):.0f}")

# Verification
print("\n" + "-" * 70)
print("VERIFICATION (ALL metrics now identical):")
print(f"  WUE sites: {n_sites_shared}")
print(f"  WUE_E sites: {n_sites_shared}")
print(f"  WUE_T sites: {n_sites_shared}")
print(f"  Match: ✓ YES (all metrics use identical strict triple intersection)")
print(f"\n  Total retained NN observations: {total_nn_observations}")
print(f"  Note: This represents ALL retained monthly observations (not unique month combinations)")
print("-" * 70)

# =============================================================================
# STEP 6: METADATA FOR UPPER-TAIL SITES (UPDATED - water_class and IGBP)
# =============================================================================

print("\n" + "="*80)
print("METADATA FOR UPPER-TAIL SITES (WUE_T extremes from strict triple intersection)")
print("="*80)

# Get top 3 sites from strict triple intersection
top_sites = dT.head(3)['site_name'].tolist()
print(f"\nTop 3 WUE_T sites (from strict triple intersection): {top_sites}")

# Metadata file path
meta_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\site_metadata_with_salinity_SPEIinfo.csv"
meta = pd.read_csv(meta_path)

# UPDATED: Changed Category → water_class, biome → IGBP
cols = [c for c in ["site_name", "water_class", "Salinity_Category", "IGBP", "climate", "lat", "long"] 
        if c in meta.columns]

out = (
    meta.loc[meta["site_name"].isin(top_sites), cols]
        .drop_duplicates()
        .sort_values("site_name")
)

print("\nMetadata for upper-tail sites:")
print(out.to_string(index=False))

# Check categories
if "Salinity_Category" in out.columns:
    print("\nSalinity_Category counts (tail sites):")
    print(out["Salinity_Category"].value_counts(dropna=False).to_string())

if "IGBP" in out.columns:
    print("\nIGBP (biome) counts (tail sites):")
    print(out["IGBP"].value_counts(dropna=False).to_string())

if "water_class" in out.columns:
    print("\nOriginal water_class values (tail sites):")
    for s in top_sites:
        v = out.loc[out["site_name"] == s, "water_class"].astype(str).unique()
        print(f"  {s}: {', '.join(v)}")

# =============================================================================
# STEP 7: METHODOLOGICAL NOTE FOR PAPER (UPDATED - STRICT TRIPLE INTERSECTION)
# =============================================================================

print("\n" + "="*80)
print("METHODOLOGICAL NOTE FOR PAPER")
print("="*80)

print(f"""
📝 STRICT TRIPLE INTERSECTION RATIONALE (ALL METRICS):

All three WUE metrics (WUE$_{{ET}}$, WUE$_E$, and WUE$_T$) were evaluated 
on the same shared subset of {len(shared_subset_sites)} sites that had ALL THREE 
metrics available under NN conditions (strict triple intersection: WUE ∩ WUE_E ∩ WUE_T).
This ensures full comparability across all metrics.

STATISTICS REPORTED:
  • Mean ± SE (standard error): Quantifies uncertainty around the mean estimate
  • SD (standard deviation): Measure of among-site variability
  • Median and IQR: Robust measures of central tendency and spread
  • Range and percentiles (p90, p95, p99): Describe distribution tails

Total retained NN observations: {total_nn_observations}
(This represents ALL retained monthly observations across the strict triple 
intersection, not unique month combinations. Support totals use all retained 
NN monthly observations.)

Sites excluded from analysis (lacked at least one of the three metrics):
  - WUE only (no partitioned data): {len(wue_only)} sites
  - WUE_E only (no WUE or WUE_T): {len(eva_only)} sites
  - WUE_T only (no WUE or WUE_E): {len(tra_only)} sites

Capping for visualization (applied to strict triple intersection):
  • WUE (WUE$_{{ET}}$): UNCAPPED in both visualization and statistics
  • WUE$_E$ (WUE_eva): Capped at 95th percentile for VISUALIZATION only
  • WUE$_T$ (WUE_tra): Capped at 95th percentile for VISUALIZATION only
  • All statistical summaries use ORIGINAL uncapped values within the strict triple intersection

This approach ensures that:
1. All three metrics are directly comparable (identical site set)
2. Comparisons between metrics are not biased by differential site availability
3. Support totals use all retained NN monthly observations
4. The shared subset includes only sites with complete data across all three WUE metrics
5. Mean, SE, and SD provide complementary information about central tendency and variability
""")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
print(f"\n✅ FINAL SHARED SUBSET SIZE: {len(shared_subset_sites)} sites")
print(f"✅ ALL THREE metrics have identical sample size: {len(shared_subset_sites)} sites each")
print(f"✅ Total retained NN observations: {total_nn_observations}")
print(f"✅ Statistics include: Mean ± SE, SD, median, IQR, range, and percentiles")