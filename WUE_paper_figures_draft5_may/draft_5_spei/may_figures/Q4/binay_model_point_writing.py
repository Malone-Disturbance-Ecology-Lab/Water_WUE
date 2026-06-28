# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 10:58:46 2025

@author: ammar
"""

import pandas as pd

# ------------------------------------------------------------------
# LOAD DATA USED FOR MODEL TRAINING
# ------------------------------------------------------------------
data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\final_dataset_with_spei_anomalies_and_classes.csv"
df = pd.read_csv(data_path)

metrics = {
    "WUE": "WUE_minus_NNmed_SPEI_48_Class",
    "WUE_eva": "WUE_eva_minus_NNmed_SPEI_48_Class",
    "WUE_tra": "WUE_tra_minus_NNmed_SPEI_48_Class"
}

spei_col = "SPEI_48"
site_col = "site_name"

print("\n" + "="*80)
print("DIAGNOSTIC SUMMARY: SPEI-48 PROBABILITY MODEL (FIGURE 1 / FIGURE 2)")
print("="*80)

for metric, class_col in metrics.items():
    
    # Filter exactly as done in the model
    mask = (
        df[class_col].isin(["Increase", "Decrease"]) &
        df[spei_col].notna()
    )
    sub = df.loc[mask].copy()
    
    n_obs = len(sub)
    n_sites = sub[site_col].nunique()
    n_unique_spei = sub[spei_col].nunique()
    
    n_dec = (sub[class_col] == "Decrease").sum()
    n_inc = (sub[class_col] == "Increase").sum()
    
    print(f"\nMETRIC: {metric}")
    print("-"*60)
    print(f"Observations used (n):      {n_obs}")
    print(f"Unique sites used:          {n_sites}")
    print(f"Unique SPEI-48 values:      {n_unique_spei}")
    print(f"Decrease class (n):         {n_dec}")
    print(f"Increase class (n):         {n_inc}")
    print(f"Decrease proportion (%):   {100*n_dec/n_obs:.1f}")
    
    # Site-level check (optional but useful)
    sites_per_class = (
        sub.groupby(site_col)[class_col]
        .value_counts()
        .unstack(fill_value=0)
    )
    
    sites_with_both = ((sites_per_class > 0).all(axis=1)).sum()
    print(f"Sites with both classes:    {sites_with_both} / {n_sites}")

print("\n" + "="*80)
print("END DIAGNOSTIC")
print("="*80)

###################################################################################################

# -*- coding: utf-8 -*-
"""
DIAGNOSTIC: Predicted P(Decrease) at selected SPEI-48 values (unique-SPEI model)
- Reads your saved model summary (intercept + slope) and computes:
  1) P(Decrease | SPEI) at SPEI = -1, -1.5, -2, -2.5, -3
  2) Stepwise % increase in probability moving between those SPEI levels

Output:
- Prints a clean, Results-ready block for each metric
- Saves a CSV table next to the summary file
"""

import os
import numpy as np
import pandas as pd

# =============================================================================
# USER CONFIG (edit if needed)
# =============================================================================
SUMMARY_CSV = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results\linear_models\WUE_Linear_Model_Analysis_Summary.csv"

# SPEI points you want to report in Results
SPEI_POINTS = [-1.0, -1.5, -2.0, -2.5, -3.0]

# Metrics to report (must match what's in your summary)
METRICS = ["WUE", "WUE_eva", "WUE_tra"]

# =============================================================================
# HELPERS
# =============================================================================
def logistic(z):
    return 1.0 / (1.0 + np.exp(-z))

def find_col(df, candidates):
    """Return the first column name in df that matches any candidate (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None

def get_params_for_metric(dfsum, metric):
    """
    Robustly extract intercept and slope for a given metric from your summary CSV.
    Handles common naming variants.
    """
    # Try to identify a "metric" column
    metric_col = find_col(dfsum, ["metric", "Metric", "WUE_metric", "wue_metric", "name", "Name"])
    if metric_col is None:
        raise ValueError("Could not find a metric/name column in the summary CSV.")

    # Filter row for metric
    sub = dfsum[dfsum[metric_col].astype(str).str.strip() == metric].copy()
    if sub.empty:
        # Sometimes labels are WUE$_E$ etc. Try partial match:
        sub = dfsum[dfsum[metric_col].astype(str).str.contains(metric, case=False, regex=False)].copy()
    if sub.empty:
        raise ValueError(f"Metric '{metric}' not found in summary CSV. Available: {dfsum[metric_col].unique()}")

    row = sub.iloc[0]

    # Candidate columns for intercept and slope
    intercept_col = find_col(dfsum, ["intercept", "Intercept", "beta0", "β0", "b0", "const", "Constant"])
    slope_col     = find_col(dfsum, ["coef", "Coef", "slope", "Slope", "beta1", "β1", "b1", "Linear coef", "linear_coef"])

    if intercept_col is None or slope_col is None:
        # If summary CSV only has β1 (slope), we can't compute absolute probabilities without β0.
        missing = []
        if intercept_col is None: missing.append("intercept (β0)")
        if slope_col is None: missing.append("slope (β1)")
        raise ValueError(
            "Missing required parameter(s) in summary CSV: " + ", ".join(missing) +
            ".\nYour summary must include BOTH β0 and β1 to compute probabilities."
        )

    b0 = float(row[intercept_col])
    b1 = float(row[slope_col])

    return b0, b1, metric_col, intercept_col, slope_col

# =============================================================================
# MAIN
# =============================================================================
dfsum = pd.read_csv(SUMMARY_CSV)

all_rows = []
print("="*90)
print("DIAGNOSTIC: P(Decrease | SPEI-48) at selected SPEI values (Results-ready)")
print("="*90)

for metric in METRICS:
    b0, b1, metric_col, intercept_col, slope_col = get_params_for_metric(dfsum, metric)

    xs = np.array(SPEI_POINTS, dtype=float)
    ps = logistic(b0 + b1 * xs)

    # Stepwise changes
    p_prev = ps[:-1]
    p_next = ps[1:]
    abs_pp_change = (p_next - p_prev) * 100.0                # percentage-points
    rel_pct_change = ((p_next - p_prev) / p_prev) * 100.0    # relative % increase

    # Print Results-ready block
    print(f"\nMETRIC: {metric}")
    print("-"*90)
    print(f"Logit parameters: β0 = {b0:+.4f}, β1 = {b1:+.4f}")
    print("\nPredicted probability of 'Decrease' at selected SPEI-48 values:")
    for x, p in zip(xs, ps):
        print(f"  SPEI = {x:>4.1f}  ->  P(Decrease) = {p:0.3f}  ({p*100:0.1f}%)")

    print("\nStepwise change in P(Decrease) as SPEI becomes more negative:")
    for i in range(len(xs)-1):
        print(
            f"  {xs[i]:>4.1f} → {xs[i+1]:>4.1f}: "
            f"+{abs_pp_change[i]:0.2f} pp  "
            f"({rel_pct_change[i]:0.1f}% relative increase)"
        )

    # Store tidy output rows for CSV
    for x, p in zip(xs, ps):
        all_rows.append({
            "metric": metric,
            "SPEI": x,
            "p_decrease": p,
            "p_decrease_percent": p*100.0,
            "beta0_intercept": b0,
            "beta1_slope": b1
        })

    # also store stepwise rows
    for i in range(len(xs)-1):
        all_rows.append({
            "metric": metric,
            "SPEI": f"{xs[i]:.1f}->{xs[i+1]:.1f}",
            "p_decrease": np.nan,
            "p_decrease_percent": np.nan,
            "beta0_intercept": b0,
            "beta1_slope": b1,
            "abs_change_pp": abs_pp_change[i],
            "rel_change_percent": rel_pct_change[i]
        })

# Save a CSV next to the summary
out_csv = os.path.join(os.path.dirname(SUMMARY_CSV), "Diagnostic_PDecrease_SelectedSPEI.csv")
pd.DataFrame(all_rows).to_csv(out_csv, index=False)

print("\n" + "="*90)
print("Saved table:", out_csv)
print("="*90)

##################################################################################################





