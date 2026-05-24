# =============================================================================
# DIAGNOSTIC: Do positive vs negative WUE_T-SPEI48 slopes differ in SPEI-48 exposure?
# Updated full version with merge-column fix
# =============================================================================

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results")

SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI48_slopes.csv"
MONTHLY_FILE = BASE_DIR / "monthly_data_after_outlier_removal.csv"
NN_MEDIANS_FILE = BASE_DIR / "wue_site_level_NN_medians_SPEI_1.csv"

OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "diagnostic_SPEI48_distribution_by_slope_sign"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 90)
print("DIAGNOSTIC: SPEI-48 DISTRIBUTION BY WUE_T SENSITIVITY SIGN")
print("=" * 90)

# =============================================================================
# STEP 1: LOAD SPEI-48 SLOPES
# =============================================================================

df_slopes = pd.read_csv(SLOPES_FILE)

required_cols = ["site_name", "slope_theilsen"]
missing = [c for c in required_cols if c not in df_slopes.columns]
if missing:
    raise ValueError(f"Missing required columns in slopes file: {missing}")

print("\n[STEP 1] Loaded SPEI-48 slopes")
print(f"  Sites in slopes file: {len(df_slopes)}")
print(f"  Columns: {df_slopes.columns.tolist()}")

# =============================================================================
# STEP 2: LOAD MONTHLY DATA AND APPLY SAME FILTERS
# =============================================================================

df_monthly = pd.read_csv(MONTHLY_FILE)
nn_data = pd.read_csv(NN_MEDIANS_FILE)

# Strict triple intersection sites
wue_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE"]["site_name"].dropna().unique())
eva_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE_eva"]["site_name"].dropna().unique())
tra_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE_tra"]["site_name"].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

df_monthly = df_monthly[df_monthly["site_name"].isin(shared_sites)].copy()

# Strict month filter: same valid WUE metrics
strict_mask = (
    df_monthly["WUE"].notna() &
    df_monthly["WUE_eva"].notna() &
    df_monthly["WUE_tra"].notna()
)
df_monthly = df_monthly[strict_mask].copy()

# Trans_ratio filter
if "Trans_ratio" in df_monthly.columns:
    before = len(df_monthly)
    df_monthly = df_monthly[df_monthly["Trans_ratio"] > 0].copy()
    print(f"\n[STEP 2] Trans_ratio > 0 filter removed {before - len(df_monthly)} rows")

# SPEI-48 and WUE_T valid rows
df_monthly = df_monthly.dropna(subset=["SPEI_48", "WUE_tra"]).copy()

print("\n[STEP 2] Monthly data after filters")
print(f"  Rows: {len(df_monthly)}")
print(f"  Sites: {df_monthly['site_name'].nunique()}")

# =============================================================================
# STEP 3: COMPUTE SITE-LEVEL SPEI-48 DISTRIBUTION METRICS
# =============================================================================

def fraction_between(x, low, high):
    return ((x >= low) & (x <= high)).mean()

site_rows = []

for site, g in df_monthly.groupby("site_name"):
    x = g["SPEI_48"].dropna()
    if len(x) >= 5:
        q25 = x.quantile(0.25)
        q75 = x.quantile(0.75)

        site_rows.append({
            "site_name": site,
            "n_months": len(x),
            "spei48_min": x.min(),
            "spei48_q25": q25,
            "spei48_median": x.median(),
            "spei48_q75": q75,
            "spei48_max": x.max(),
            "spei48_range": x.max() - x.min(),
            "spei48_iqr": q75 - q25,
            "spei48_mean": x.mean(),
            "frac_dry_neg": (x < 0).mean(),
            "frac_wet_pos": (x > 0).mean(),
            "frac_extreme_dry": (x <= -1).mean(),
            "frac_extreme_wet": (x >= 1).mean(),
            "frac_near_normal": fraction_between(x, -0.5, 0.5),
        })

df_spei = pd.DataFrame(site_rows)

print("\n[STEP 3] Site-level SPEI-48 distribution metrics computed")
print(f"  Sites: {len(df_spei)}")
print(f"  Columns: {df_spei.columns.tolist()}")

# =============================================================================
# STEP 4: MERGE WITH SLOPES AND FIX DUPLICATE COLUMN NAMES
# =============================================================================

df = df_slopes.merge(df_spei, on="site_name", how="inner", suffixes=("_slopesfile", "_fresh"))

# Prefer freshly computed monthly SPEI-48 metrics whenever duplicate columns exist
for col in ["spei48_min", "spei48_max", "spei48_range"]:
    fresh_col = f"{col}_fresh"
    old_col = f"{col}_slopesfile"

    if fresh_col in df.columns:
        df[col] = df[fresh_col]
    elif old_col in df.columns:
        df[col] = df[old_col]
    else:
        raise ValueError(f"Could not find either {fresh_col} or {old_col}")

# Median may exist in both files
if "spei48_median_fresh" in df.columns:
    df["spei48_median"] = df["spei48_median_fresh"]
elif "spei48_median_slopesfile" in df.columns:
    df["spei48_median"] = df["spei48_median_slopesfile"]
else:
    raise ValueError("Could not find SPEI-48 median column after merge")

# n_months may exist in both files
if "n_months_fresh" in df.columns:
    df["n_months"] = df["n_months_fresh"]
elif "n_months_slopesfile" in df.columns:
    df["n_months"] = df["n_months_slopesfile"]

# Define slope sign
df["slope_sign"] = np.where(
    df["slope_theilsen"] > 0, "Positive",
    np.where(df["slope_theilsen"] < 0, "Negative", "Zero")
)

# Safety check
needed_cols = [
    "spei48_min",
    "spei48_max",
    "spei48_median",
    "spei48_range",
    "spei48_q25",
    "spei48_q75",
    "spei48_iqr",
    "spei48_mean",
    "frac_dry_neg",
    "frac_wet_pos",
    "frac_extreme_dry",
    "frac_extreme_wet",
    "frac_near_normal",
    "n_months",
    "slope_theilsen",
    "slope_sign"
]

missing_cols = [c for c in needed_cols if c not in df.columns]
if missing_cols:
    raise ValueError(f"Missing columns after merge fix: {missing_cols}")

print("\n[STEP 4] Merged data")
print(f"  Sites merged: {len(df)}")
print(f"  Positive slopes: {(df['slope_theilsen'] > 0).sum()}")
print(f"  Negative slopes: {(df['slope_theilsen'] < 0).sum()}")
print(f"  Zero slopes: {(df['slope_theilsen'] == 0).sum()}")

print("\nColumns after merge/fix:")
print(df.columns.tolist())

# =============================================================================
# STEP 5: COMPARE SPEI-48 DISTRIBUTIONS BETWEEN POSITIVE AND NEGATIVE SLOPE SITES
# =============================================================================

df_sign = df[df["slope_sign"].isin(["Positive", "Negative"])].copy()

comparison_cols = [
    "spei48_min",
    "spei48_q25",
    "spei48_median",
    "spei48_q75",
    "spei48_max",
    "spei48_range",
    "spei48_iqr",
    "spei48_mean",
    "frac_dry_neg",
    "frac_wet_pos",
    "frac_extreme_dry",
    "frac_extreme_wet",
    "frac_near_normal",
    "n_months"
]

results = []

positive = df_sign[df_sign["slope_sign"] == "Positive"]
negative = df_sign[df_sign["slope_sign"] == "Negative"]

print("\n" + "=" * 90)
print("TEST 1: POSITIVE VS NEGATIVE SLOPE SITES — SPEI-48 DISTRIBUTION")
print("=" * 90)
print(f"Positive-slope sites: n = {len(positive)}")
print(f"Negative-slope sites: n = {len(negative)}")

for col in comparison_cols:
    pos_vals = positive[col].dropna()
    neg_vals = negative[col].dropna()

    if len(pos_vals) > 0 and len(neg_vals) > 0:
        u, p = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")

        results.append({
            "metric": col,
            "positive_median": pos_vals.median(),
            "negative_median": neg_vals.median(),
            "positive_mean": pos_vals.mean(),
            "negative_mean": neg_vals.mean(),
            "U": u,
            "p_value": p
        })

        sig = " *" if p < 0.05 else ""
        print(f"\n{col}:")
        print(f"  Positive median = {pos_vals.median():.3f}")
        print(f"  Negative median = {neg_vals.median():.3f}")
        print(f"  Mann-Whitney U = {u:.1f}, p = {p:.4f}{sig}")

results_df = pd.DataFrame(results)
results_csv = OUTPUT_DIR / "positive_vs_negative_spei48_distribution_tests.csv"
results_df.to_csv(results_csv, index=False)

print(f"\nSaved comparison table: {results_csv}")

# =============================================================================
# STEP 6: OVERLAP CHECK — COMMON SPEI-48 RANGE SUBSETS
# =============================================================================

print("\n" + "=" * 90)
print("TEST 2: OVERLAP SUBSET CHECK")
print("=" * 90)

# Subset A: includes moderate dry and wet exposure
overlap_subset = df[
    (df["spei48_min"] <= -0.5) &
    (df["spei48_max"] >= 0.5)
].copy()

print("\nSubset A: Sites with both dry and wet SPEI-48 exposure")
print("Criteria: spei48_min <= -0.5 and spei48_max >= 0.5")
print(f"  Sites retained: {len(overlap_subset)} / {len(df)}")
print(f"  Positive slopes: {(overlap_subset['slope_theilsen'] > 0).sum()}")
print(f"  Negative slopes: {(overlap_subset['slope_theilsen'] < 0).sum()}")
if len(overlap_subset) > 0:
    print(f"  Median slope: {overlap_subset['slope_theilsen'].median():.3f}")
    print(f"  Slope range: {overlap_subset['slope_theilsen'].min():.3f} to {overlap_subset['slope_theilsen'].max():.3f}")
    print(f"  SPEI-48 range: {overlap_subset['spei48_min'].min():.2f} to {overlap_subset['spei48_max'].max():.2f}")

# Subset B: includes stronger dry and wet exposure
strict_overlap_subset = df[
    (df["spei48_min"] <= -1.0) &
    (df["spei48_max"] >= 1.0)
].copy()

print("\nSubset B: Sites with both stronger dry and wet SPEI-48 exposure")
print("Criteria: spei48_min <= -1.0 and spei48_max >= 1.0")
print(f"  Sites retained: {len(strict_overlap_subset)} / {len(df)}")
print(f"  Positive slopes: {(strict_overlap_subset['slope_theilsen'] > 0).sum()}")
print(f"  Negative slopes: {(strict_overlap_subset['slope_theilsen'] < 0).sum()}")
if len(strict_overlap_subset) > 0:
    print(f"  Median slope: {strict_overlap_subset['slope_theilsen'].median():.3f}")
    print(f"  Slope range: {strict_overlap_subset['slope_theilsen'].min():.3f} to {strict_overlap_subset['slope_theilsen'].max():.3f}")
    print(f"  SPEI-48 range: {strict_overlap_subset['spei48_min'].min():.2f} to {strict_overlap_subset['spei48_max'].max():.2f}")

# =============================================================================
# STEP 7: SIMILAR SPEI-48 RANGE SUBSET
# =============================================================================

print("\n" + "=" * 90)
print("TEST 3: SIMILAR SPEI-48 RANGE SUBSET")
print("=" * 90)

range_q25 = df["spei48_range"].quantile(0.25)
range_q75 = df["spei48_range"].quantile(0.75)

similar_range_subset = df[
    (df["spei48_range"] >= range_q25) &
    (df["spei48_range"] <= range_q75)
].copy()

print("\nMiddle-IQR SPEI-48 range subset")
print(f"Criteria: spei48_range between Q25={range_q25:.3f} and Q75={range_q75:.3f}")
print(f"  Sites retained: {len(similar_range_subset)} / {len(df)}")
print(f"  Positive slopes: {(similar_range_subset['slope_theilsen'] > 0).sum()}")
print(f"  Negative slopes: {(similar_range_subset['slope_theilsen'] < 0).sum()}")
if len(similar_range_subset) > 0:
    print(f"  Median slope: {similar_range_subset['slope_theilsen'].median():.3f}")
    print(f"  Slope range: {similar_range_subset['slope_theilsen'].min():.3f} to {similar_range_subset['slope_theilsen'].max():.3f}")
    print(f"  SPEI-48 range subset: {similar_range_subset['spei48_range'].min():.3f} to {similar_range_subset['spei48_range'].max():.3f}")

# =============================================================================
# STEP 8: LOGISTIC CHECK — CAN SPEI DISTRIBUTION METRICS PREDICT SLOPE SIGN?
# =============================================================================

print("\n" + "=" * 90)
print("TEST 4: CAN SPEI-48 DISTRIBUTION METRICS PREDICT SLOPE SIGN?")
print("=" * 90)

feature_cols = [
    "spei48_min",
    "spei48_max",
    "spei48_median",
    "spei48_range",
    "spei48_iqr",
    "frac_dry_neg",
    "frac_extreme_dry",
    "frac_extreme_wet"
]

df_model = df_sign.dropna(subset=feature_cols + ["slope_theilsen"]).copy()
y = (df_model["slope_theilsen"] > 0).astype(int).values
X = df_model[feature_cols].values

if len(np.unique(y)) == 2 and len(df_model) >= 10:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = LogisticRegression(max_iter=1000, class_weight="balanced", random_state=42)
    clf.fit(X_scaled, y)

    pred_prob = clf.predict_proba(X_scaled)[:, 1]
    auc = roc_auc_score(y, pred_prob)

    print("\nLogistic check using SPEI-48 distribution summaries")
    print(f"  N = {len(df_model)} sites")
    print(f"  AUC = {auc:.3f}")
    print("  Note: this is in-sample and diagnostic only, not a predictive model.")

    coef_table = pd.DataFrame({
        "feature": feature_cols,
        "coefficient": clf.coef_[0]
    }).sort_values("coefficient", ascending=False)

    print("\nCoefficients:")
    print(coef_table.to_string(index=False))

    coef_csv = OUTPUT_DIR / "logistic_spei_distribution_predicts_slope_sign.csv"
    coef_table.to_csv(coef_csv, index=False)
    print(f"\nSaved coefficient table: {coef_csv}")
else:
    print("Not enough data or only one slope sign present for logistic check.")

# =============================================================================
# STEP 9: SAVE FULL SITE-LEVEL OUTPUT
# =============================================================================

output_site_csv = OUTPUT_DIR / "site_level_spei48_distribution_by_slope_sign.csv"
df.to_csv(output_site_csv, index=False)

print("\nSaved site-level output:")
print(output_site_csv)

# =============================================================================
# FINAL TAKE-HOME TEMPLATE
# =============================================================================

print("\n" + "=" * 90)
print("FINAL TAKE-HOME TEMPLATE")
print("=" * 90)

print("""
Interpretation guide:

1. If SPEI min/max/median/range/IQR do not differ strongly between positive and negative slope sites:
   Positive and negative sensitivities occurred across comparable SPEI-48 exposure distributions.

2. If both signs remain in the overlap subset:
   Sensitivity direction was not simply determined by whether sites experienced dry/wet SPEI-48 conditions.

3. If the middle-IQR SPEI range subset still has both signs:
   Opposite sensitivity directions also occurred among sites with similar SPEI-48 range lengths.

4. If the logistic AUC is low/moderate:
   Simple SPEI-48 distribution summaries did not strongly classify slope direction.
""")

print("\n" + "=" * 90)
print("DIAGNOSTIC COMPLETE")
print("=" * 90)