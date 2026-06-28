"""
VALUE-BASED DIAGNOSTIC: Compare Python workflow outputs against Malone reference
Uses the SAME file names for both sides (Q1_simple_* files)
"""

import os
import numpy as np
import pandas as pd

print("\n" + "="*80)
print("VALUE-BASED DIAGNOSTIC: PYTHON OUTPUTS VS MALONE REFERENCE")
print("="*80)

# ---------------------------------------------------------------------------
# PATHS - EDIT THESE IF NEEDED
# ---------------------------------------------------------------------------

# Malone reference outputs (the original Malone workflow)
malone_dir = r"M:\Research\WUE_CUE\Water_WUE\Malone_Workflow\results\Q1_WUE_metric_difference_simple"

# Python workflow outputs (your new workflow)
python_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_mixed_model_outputs"

# Diagnostic output location
diagnostic_out = os.path.join(python_dir, "Q1_value_comparison_against_malone.csv")

print(f"\nMalone reference directory: {malone_dir}")
print(f"Python output directory:   {python_dir}")

# ---------------------------------------------------------------------------
# FILE MAPPING: SAME NAMES ON BOTH SIDES (Q1_simple_* files)
# ---------------------------------------------------------------------------

file_map = {
    "mixed_fixed": {
        "malone_file": "Q1_simple_mixed_fixed_effects_with_approx_p.csv",
        "python_file": "Q1_simple_mixed_fixed_effects_with_approx_p.csv",
        "key_cols": ["term"]
    },
    "mixed_lrt": {
        "malone_file": "Q1_simple_mixed_likelihood_ratio_tests.csv",
        "python_file": "Q1_simple_mixed_likelihood_ratio_tests.csv",
        "key_cols": ["term"]
    },
    "mixed_comparison": {
        "malone_file": "Q1_simple_mixed_model_comparison_lrt.csv",
        "python_file": "Q1_simple_mixed_model_comparison_lrt.csv",
        "key_cols": None  # Compare row-by-row
    },
    "gam_comparison": {
        "malone_file": "Q1_simple_gam_model_comparison.csv",
        "python_file": "Q1_simple_gam_model_comparison.csv",
        "key_cols": ["model"]
    },
    "gam_param": {
        "malone_file": "Q1_simple_gam_parametric_terms.csv",
        "python_file": "Q1_simple_gam_parametric_terms.csv",
        "key_cols": ["term"]
    },
    "gam_smooth": {
        "malone_file": "Q1_simple_gam_smooth_terms.csv",
        "python_file": "Q1_simple_gam_smooth_terms.csv",
        "key_cols": ["smooth_term"]
    },
    "mixed_predictions": {
        "malone_file": "Q1_simple_mixed_predictions_TET.csv",
        "python_file": "Q1_simple_mixed_predictions_TET.csv",
        "key_cols": ["difference_type", "water_class", "Trans_ratio"]
    },
    "gam_predictions": {
        "malone_file": "Q1_simple_gam_predictions_TET.csv",
        "python_file": "Q1_simple_gam_predictions_TET.csv",
        "key_cols": ["difference_type", "water_class", "Trans_ratio"]
    },
    "mixed_fitted": {
        "malone_file": "Q1_simple_mixed_fitted_residuals.csv",
        "python_file": "Q1_simple_mixed_fitted_residuals.csv",
        "key_cols": ["site_name", "Year", "month", "difference_type"]
    },
    "gam_fitted": {
        "malone_file": "Q1_simple_gam_fitted_residuals.csv",
        "python_file": "Q1_simple_gam_fitted_residuals.csv",
        "key_cols": ["site_name", "Year", "month", "difference_type"]
    },
    "summary": {
        "malone_file": "Q1_simple_summary_by_difference_type.csv",
        "python_file": "Q1_simple_summary_by_difference_type.csv",
        "key_cols": ["difference_type", "difference_label", "water_class"]
    },
}

# ---------------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------------

def clean_colnames(df):
    """Normalize column names for better comparison."""
    df = df.copy()
    df.columns = (
        df.columns
        .str.strip()
        .str.replace(" ", "_", regex=False)
        .str.replace("Pr(Chi)", "Pr_Chi", regex=False)
        .str.replace("Pr(>Chisq)", "Pr_Chisq", regex=False)
        .str.replace("p-value", "p_value", regex=False)
        .str.replace("p.value", "p_value", regex=False)
    )
    return df

def add_row_id(df):
    df = df.copy()
    df["_row_id"] = np.arange(len(df))
    return df

def numeric_columns_common(a, b, key_cols):
    """Find numeric columns that exist in both dataframes."""
    common = sorted(set(a.columns).intersection(set(b.columns)))
    key_cols = [] if key_cols is None else key_cols
    candidates = [c for c in common if c not in key_cols]
    numeric = []
    for c in candidates:
        a_num = pd.to_numeric(a[c], errors="coerce")
        b_num = pd.to_numeric(b[c], errors="coerce")
        if a_num.notna().any() or b_num.notna().any():
            numeric.append(c)
    return numeric

def compare_table(label, malone_file, python_file, key_cols=None, atol=1e-5, rtol=1e-5):
    """Compare a single table between Malone and Python outputs."""
    rows = []

    malone_path = os.path.join(malone_dir, malone_file)
    python_path = os.path.join(python_dir, python_file)

    # Check file existence
    if not os.path.exists(malone_path):
        rows.append({
            "table": label,
            "check": "FILE_EXISTS",
            "status": "FAIL",
            "message": f"Missing Malone file: {malone_file}",
            "n_compared": np.nan,
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
        return rows

    if not os.path.exists(python_path):
        rows.append({
            "table": label,
            "check": "FILE_EXISTS",
            "status": "FAIL",
            "message": f"Missing Python file: {python_file}",
            "n_compared": np.nan,
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
        return rows

    # Load and clean
    a = clean_colnames(pd.read_csv(malone_path))
    b = clean_colnames(pd.read_csv(python_path))

    # Check row counts
    if len(a) != len(b):
        rows.append({
            "table": label,
            "check": "ROW_COUNT",
            "status": "FAIL",
            "message": f"Malone rows={len(a)}; Python rows={len(b)}",
            "n_compared": min(len(a), len(b)),
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
    else:
        rows.append({
            "table": label,
            "check": "ROW_COUNT",
            "status": "PASS",
            "message": f"Both have {len(a)} rows",
            "n_compared": len(a),
            "max_abs_diff": 0.0,
            "max_rel_diff": 0.0,
        })

    # Check column sets
    missing_in_python = sorted(set(a.columns) - set(b.columns))
    extra_in_python = sorted(set(b.columns) - set(a.columns))

    if missing_in_python or extra_in_python:
        rows.append({
            "table": label,
            "check": "COLUMN_SET",
            "status": "FAIL",
            "message": f"Missing in Python={missing_in_python}; Extra in Python={extra_in_python}",
            "n_compared": np.nan,
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
    else:
        rows.append({
            "table": label,
            "check": "COLUMN_SET",
            "status": "PASS",
            "message": f"Column sets match ({len(a.columns)} columns)",
            "n_compared": len(a.columns),
            "max_abs_diff": 0.0,
            "max_rel_diff": 0.0,
        })

    # If no key columns, use row-by-row comparison
    if key_cols is None:
        a = add_row_id(a)
        b = add_row_id(b)
        key_cols = ["_row_id"]
    else:
        key_cols = [k for k in key_cols if k in a.columns and k in b.columns]

    if not key_cols:
        rows.append({
            "table": label,
            "check": "KEY_MATCH",
            "status": "FAIL",
            "message": "No valid common key columns found",
            "n_compared": np.nan,
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
        return rows

    # For prediction grids, round Trans_ratio to avoid precision issues
    if "Trans_ratio" in key_cols:
        a["Trans_ratio"] = pd.to_numeric(a["Trans_ratio"], errors="coerce").round(10)
        b["Trans_ratio"] = pd.to_numeric(b["Trans_ratio"], errors="coerce").round(10)

    # Merge on key columns
    merged = a.merge(
        b,
        on=key_cols,
        how="outer",
        suffixes=("_malone", "_python"),
        indicator=True
    )

    only_malone = (merged["_merge"] == "left_only").sum()
    only_python = (merged["_merge"] == "right_only").sum()
    both = (merged["_merge"] == "both").sum()

    if only_malone or only_python:
        rows.append({
            "table": label,
            "check": "KEY_MATCH",
            "status": "FAIL",
            "message": f"Rows matched={both}; Only Malone={only_malone}; Only Python={only_python}",
            "n_compared": both,
            "max_abs_diff": np.nan,
            "max_rel_diff": np.nan,
        })
    else:
        rows.append({
            "table": label,
            "check": "KEY_MATCH",
            "status": "PASS",
            "message": f"All rows matched: {both}",
            "n_compared": both,
            "max_abs_diff": 0.0,
            "max_rel_diff": 0.0,
        })

    # Compare numeric values
    common_numeric = numeric_columns_common(a, b, key_cols)

    for col in common_numeric:
        col_a = f"{col}_malone"
        col_b = f"{col}_python"

        if col_a not in merged.columns or col_b not in merged.columns:
            continue

        x = pd.to_numeric(merged.loc[merged["_merge"] == "both", col_a], errors="coerce")
        y = pd.to_numeric(merged.loc[merged["_merge"] == "both", col_b], errors="coerce")

        valid = x.notna() & y.notna()
        if valid.sum() == 0:
            continue

        diff = (x[valid] - y[valid]).abs()
        rel = diff / np.maximum(np.abs(x[valid]), 1e-12)

        max_abs = diff.max()
        max_rel = rel.max()
        passed = np.allclose(x[valid], y[valid], atol=atol, rtol=rtol)

        rows.append({
            "table": label,
            "check": f"VALUES_{col}",
            "status": "PASS" if passed else "FAIL",
            "message": "values match" if passed else f"values differ (max_abs={max_abs:.2e}, max_rel={max_rel:.2e})",
            "n_compared": int(valid.sum()),
            "max_abs_diff": float(max_abs),
            "max_rel_diff": float(max_rel),
        })

    return rows

# ---------------------------------------------------------------------------
# RUN ALL COMPARISONS
# ---------------------------------------------------------------------------

all_rows = []

for label, config in file_map.items():
    print(f"\nComparing {label}...")
    all_rows.extend(
        compare_table(
            label=label,
            malone_file=config["malone_file"],
            python_file=config["python_file"],
            key_cols=config["key_cols"],
            atol=1e-5,
            rtol=1e-5,
        )
    )

# Save diagnostic results
diagnostics = pd.DataFrame(all_rows)
diagnostics.to_csv(diagnostic_out, index=False)

# ---------------------------------------------------------------------------
# PRINT SUMMARY
# ---------------------------------------------------------------------------

print("\n" + "="*80)
print("VALUE COMPARISON SUMMARY")
print("="*80)

# Summary by table and status
summary = (
    diagnostics
    .groupby(["table", "status"])
    .size()
    .reset_index(name="n_checks")
)
print("\nSummary by table and status:")
print(summary.to_string(index=False))

# Overall status
failed = diagnostics[diagnostics["status"] == "FAIL"]
passed = diagnostics[diagnostics["status"] == "PASS"]

print(f"\nTotal checks: {len(diagnostics)}")
print(f"  Passed: {len(passed)}")
print(f"  Failed: {len(failed)}")

if len(failed) == 0:
    print("\n✓ ALL CHECKS PASSED! Python outputs match Malone reference.")
else:
    print(f"\n✗ {len(failed)} checks failed:")

# Show detailed failures
if len(failed) > 0:
    print("\nFailed checks (detailed):")
    print("="*80)
    
    # Group failures by table
    for table in failed["table"].unique():
        table_failures = failed[failed["table"] == table]
        print(f"\n{table}:")
        for _, row in table_failures.iterrows():
            print(f"  - {row['check']}: {row['message']}")
            
        # Show max diffs for value failures
        value_failures = table_failures[table_failures["check"].str.startswith("VALUES_")]
        if len(value_failures) > 0:
            print(f"    Max abs diff: {value_failures['max_abs_diff'].max():.2e}")
            print(f"    Max rel diff: {value_failures['max_rel_diff'].max():.2e}")

print(f"\nDiagnostic saved to: {diagnostic_out}")

# ---------------------------------------------------------------------------
# OPTIONAL: CHECK FOR STALE FILES IN OUTPUT FOLDER
# ---------------------------------------------------------------------------

print("\n" + "="*80)
print("CLEANUP CHECK: Stale files in Python output folder")
print("="*80)

expected_files = {
    "Q1_simple_mixed_difference_TET_model.rds",
    "Q1_simple_smooth_gam_difference_TET_model.rds",
    "Q1_simple_mixed_fitted_residuals.csv",
    "Q1_simple_mixed_fixed_effects_with_approx_p.csv",
    "Q1_simple_mixed_likelihood_ratio_tests.csv",
    "Q1_simple_mixed_model_comparison_lrt.csv",
    "Q1_simple_gam_fitted_residuals.csv",
    "Q1_simple_gam_model_comparison.csv",
    "Q1_simple_gam_parametric_terms.csv",
    "Q1_simple_gam_smooth_terms.csv",
    "Q1_simple_mixed_predictions_TET.csv",
    "Q1_simple_gam_predictions_TET.csv",
    "Q1_simple_summary_by_difference_type.csv",
    "Q1_simple_model_summary.txt",
    "Q1_simple_outputs_README.md",
    "Q1_value_comparison_against_malone.csv",  # This diagnostic output
}

all_files = set(os.listdir(python_dir))
stale_files = all_files - expected_files

# Exclude temp files if they somehow ended up here
stale_files = {f for f in stale_files if not f.startswith("temp_") and not f.startswith("run_")}

if stale_files:
    print("\n⚠ WARNING: Stale files found in output folder (not from current workflow):")
    for f in sorted(stale_files):
        print(f"  - {f}")
    print("\nConsider removing these files to avoid confusion.")
else:
    print("\n✓ No stale files found. Output folder is clean.")

print("\n" + "="*80)
print("DIAGNOSTIC COMPLETE")
print("="*80)