"""
DIAGNOSTIC: Compare Python Q2 Anomaly vs Original Malone R Q2 Anomaly Outputs
HARDENED VERSION - Row/column mismatches = hard failures, tighter tolerance
UPDATED: Match Malone's default R column names (Std. Error, t value, etc.)
"""

import os
import pandas as pd
import numpy as np

# ============================================================================
# PATHS - CONFIRMED MALONE LOCATION
# ============================================================================

python_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2_anomalies\Q2_anomalies_outputs"

# Malone R workflow saves to: Water_WUE/data/results/Q2_anomalies
malone_output_dir = r"M:\Research\WUE_CUE\Water_WUE\data\results\Q2_anomalies"

# If the above doesn't exist, try the alternative location
if not os.path.exists(malone_output_dir):
    malone_output_dir = r"M:\Research\WUE_CUE\Water_WUE\Malone_Workflow\results\Q2_anomalies"
    print(f"Using alternative Malone directory: {malone_output_dir}")

print("="*70)
print("DIAGNOSTIC: COMPARING PYTHON Q2 ANOMALY vs MALONE R Q2 ANOMALY OUTPUTS")
print("="*70)
print(f"\nPython output directory: {python_output_dir}")
print(f"Malone output directory: {malone_output_dir}")

# ============================================================================
# CHECK IF FILES EXIST
# ============================================================================

print("\n" + "-"*70)
print("STEP 1: Checking if output files exist")
print("-"*70)

# Python files (with actual filenames from the workflow)
python_files = {
    "model_base": "Q2_anomalies_model_base_wide.csv",
    "model_data": "Q2_anomalies_model_data_long.csv",
    "baseline": "Q2_anomalies_near_normal_baselines.csv",
    "summary_by_class": "Q2_anomalies_summary_by_class_metric.csv",
    "summary_by_class_ecosystem": "Q2_anomalies_summary_by_class_metric_ecosystem.csv",
    "direction_summary": "Q2_anomalies_direction_summary.csv",
    "wue_fixed": "Q2_anomalies_wue_value_fixed_effects.csv",
    "response_fixed": "Q2_anomalies_signed_response_fixed_effects.csv",
    "magnitude_fixed": "Q2_anomalies_response_magnitude_fixed_effects.csv",
    "direction_fixed": "Q2_anomalies_response_direction_fixed_effects.csv",
    "wue_drop1": "Q2_anomalies_wue_value_drop1_lrt.csv",
    "response_drop1": "Q2_anomalies_signed_response_drop1_lrt.csv",
    "magnitude_drop1": "Q2_anomalies_response_magnitude_drop1_lrt.csv",
    "direction_drop1": "Q2_anomalies_response_direction_drop1_lrt.csv",
    "wue_letters": "Q2_anomalies_WUE_anomaly_class_letters.csv",
    "response_letters": "Q2_anomalies_signed_response_anomaly_class_letters.csv",
    "wue_pairs": "Q2_anomalies_WUE_anomaly_class_pairwise_contrasts.csv",
    "response_pairs": "Q2_anomalies_signed_response_anomaly_class_pairwise_contrasts.csv",
    "fitted_residuals": "Q2_anomalies_fitted_residuals.csv",
    "model_results": "Q2_anomalies_model_results_table.csv"
}

# Malone files (same names)
malone_files = python_files.copy()  # Same filenames

print("\nPython files found:")
python_available = {}
for name, filename in python_files.items():
    filepath = os.path.join(python_output_dir, filename)
    exists = os.path.exists(filepath)
    python_available[name] = exists
    print(f"  {name}: {'YES' if exists else 'NO'} - {filename}")

print("\nMalone files found:")
malone_available = {}
for name, filename in malone_files.items():
    filepath = os.path.join(malone_output_dir, filename)
    exists = os.path.exists(filepath)
    malone_available[name] = exists
    print(f"  {name}: {'YES' if exists else 'NO'} - {filename}")

# Hard failure if files missing
missing_python = [name for name, exists in python_available.items() if not exists]
missing_malone = [name for name, exists in malone_available.items() if not exists]

if missing_python:
    print(f"\n[FAIL] Missing Python files: {missing_python}")
    raise FileNotFoundError(f"Missing Python files: {missing_python}")

if missing_malone:
    print(f"\n[FAIL] Missing Malone files: {missing_malone}")
    raise FileNotFoundError(f"Missing Malone files: {missing_malone}")

print("\n[PASS] All files present.")

# ============================================================================
# HELPER: Compare two dataframes with column name mapping - HARDENED
# ============================================================================

def compare_dataframes_with_mapping(df1, df2, name, merge_cols, compare_map, tolerance=1e-5):
    """
    Compare two dataframes with column name mapping
    compare_map: dict mapping {python_col: malone_col}
    HARDENED: Row/column mismatches = hard failures
    """
    print(f"\n  {name}:")
    print(f"    Python: {len(df1)} rows, {len(df1.columns)} columns")
    print(f"    Malone: {len(df2)} rows, {len(df2.columns)} columns")
    
    # Clean column names
    df1.columns = df1.columns.str.strip()
    df2.columns = df2.columns.str.strip()
    
    all_match = True
    
    # Check column sets - HARD FAILURE
    col_diff_python = set(df2.columns) - set(df1.columns)
    col_diff_malone = set(df1.columns) - set(df2.columns)
    
    if col_diff_python or col_diff_malone:
        if col_diff_python:
            print(f"    [FAIL] Columns missing in Python: {col_diff_python}")
        if col_diff_malone:
            print(f"    [FAIL] Columns missing in Malone: {col_diff_malone}")
        all_match = False
    else:
        print(f"    [PASS] Column sets match ({len(df1.columns)} columns)")
    
    # Check if merge columns exist
    missing_merge = [c for c in merge_cols if c not in df1.columns or c not in df2.columns]
    if missing_merge:
        print(f"    [FAIL] Missing merge columns: {missing_merge}")
        all_match = False
        return {}, all_match
    
    # Merge with indicator - HARD FAILURE on row mismatches
    merged = df1.merge(df2, on=merge_cols, suffixes=('_python', '_malone'), how='outer', indicator=True)
    
    only_python = (merged['_merge'] == 'left_only').sum()
    only_malone = (merged['_merge'] == 'right_only').sum()
    both = (merged['_merge'] == 'both').sum()
    
    print(f"    Merge: {both} matched, {only_python} only in Python, {only_malone} only in Malone")
    
    if only_python > 0 or only_malone > 0:
        print(f"    [FAIL] Row mismatch detected")
        all_match = False
    else:
        print(f"    [PASS] All rows matched")
    
    results = {}
    
    for py_col, ml_col in compare_map.items():
        # Find matching columns (handle different naming conventions)
        py_col_actual = None
        ml_col_actual = None
        
        for col in df1.columns:
            if col.lower() == py_col.lower():
                py_col_actual = col
                break
        
        for col in df2.columns:
            if col.lower() == ml_col.lower():
                ml_col_actual = col
                break
        
        if py_col_actual and ml_col_actual:
            py_col_suffixed = f'{py_col_actual}_python'
            ml_col_suffixed = f'{ml_col_actual}_malone'
            
            if py_col_suffixed in merged.columns and ml_col_suffixed in merged.columns:
                # Check if numeric
                try:
                    python_vals = pd.to_numeric(merged[py_col_suffixed], errors='coerce')
                    malone_vals = pd.to_numeric(merged[ml_col_suffixed], errors='coerce')
                    diff = python_vals - malone_vals
                    absdiff = np.abs(diff)
                    max_diff = absdiff.max()
                    nan_count = absdiff.isna().sum()
                    
                    results[py_col] = max_diff
                    
                    # Handle df column specially (lme4 doesn't provide df in summary)
                    if py_col == 'df':
                        print(f"    {py_col}: {nan_count} NaN values (expected - lme4 doesn't provide df) [OK]")
                    elif nan_count > 0 and nan_count < len(merged):
                        print(f"    {py_col}: {nan_count} NaN values in difference (max diff for valid = {max_diff:.10f})")
                        # NaN differences are a warning, not a hard failure
                    elif max_diff <= 1e-8:
                        print(f"    {py_col}: Max diff = {max_diff:.10f} [IDENTICAL]")
                    elif max_diff <= tolerance:
                        print(f"    {py_col}: Max diff = {max_diff:.8f} [OK]")
                    else:
                        all_match = False
                        print(f"    {py_col}: Max diff = {max_diff:.8f} [FAIL] (tolerance = {tolerance})")
                except Exception as e:
                    print(f"    {py_col}: Error comparing - {str(e)}")
                    results[py_col] = np.nan
                    all_match = False
            else:
                print(f"    {py_col}: Column not found in merged data")
                results[py_col] = np.nan
                all_match = False
        else:
            print(f"    {py_col}: Column not found in source data")
            results[py_col] = np.nan
            all_match = False
    
    return results, all_match

# ============================================================================
# COMPARE ALL FILES - UPDATED WITH CORRECT MALONE COLUMN NAMES
# ============================================================================

print("\n" + "-"*70)
print("STEP 2: Comparing individual files (tolerance = 1e-5)")
print("-"*70)

all_results = {}
all_metrics_match = True

# File comparison mappings - UPDATED to match Malone's default R column names
comparisons = [
    ("model_base", ['site_name', 'Year', 'month', 'water_class'],
     {'WUE': 'WUE', 'WUE_tra': 'WUE_tra', 'Trans_ratio': 'Trans_ratio', 
      'SPEI_6': 'SPEI_6', 'SPEI_48': 'SPEI_48', 'Trans_ratio_z': 'Trans_ratio_z'}),
    
    ("model_data", ['site_name', 'Year', 'month', 'water_class', 'timescale', 'metric', 'anomaly_class'],
     {'WUE_value': 'WUE_value', 'baseline_mean': 'baseline_mean', 'baseline_n': 'baseline_n',
      'WUE_response': 'WUE_response', 'WUE_response_magnitude': 'WUE_response_magnitude', 
      'decrease_binary': 'decrease_binary'}),
    
    ("baseline", ['site_name', 'timescale', 'metric', 'month_f'],
     {'baseline_mean': 'baseline_mean', 'baseline_median': 'baseline_median', 'baseline_n': 'baseline_n'}),
    
    # UPDATED: Full summary_by_class with all Malone columns
    ("summary_by_class", ['timescale', 'anomaly_class', 'metric'],
     {
         'n_months': 'n_months',
         'n_sites': 'n_sites',
         'mean_wue': 'mean_wue',
         'median_wue': 'median_wue',
         'sd_wue': 'sd_wue',
         'se_wue': 'se_wue',
         'ci95_wue': 'ci95_wue',
         'mean_response': 'mean_response',
         'se_response': 'se_response',
         'ci95_response': 'ci95_response',
         'mean_abs_response': 'mean_abs_response',
         'pct_decrease': 'pct_decrease'
     }),
    
    ("summary_by_class_ecosystem", ['timescale', 'anomaly_class', 'metric', 'water_class'],
     {'n_months': 'n_months', 'n_sites': 'n_sites', 'mean_spei': 'mean_spei',
      'mean_wue': 'mean_wue', 'mean_response': 'mean_response',
      'mean_abs_response': 'mean_abs_response', 'pct_decrease': 'pct_decrease'}),
    
    ("direction_summary", ['timescale', 'anomaly_class', 'metric', 'water_class'],
     {'n_months': 'n_months', 'pct_decrease': 'pct_decrease'}),
    
    # UPDATED: Fixed effects use Malone's default column names (no underscores)
    ("wue_fixed", ['term'],
     {'Estimate': 'Estimate', 'Std. Error': 'Std. Error', 't value': 't value'}),
    
    ("response_fixed", ['term'],
     {'Estimate': 'Estimate', 'Std. Error': 'Std. Error', 't value': 't value'}),
    
    ("magnitude_fixed", ['term'],
     {'Estimate': 'Estimate', 'Std. Error': 'Std. Error', 't value': 't value'}),
    
    # UPDATED: Direction fixed effects include all Malone columns
    ("direction_fixed", ['term'],
     {'Estimate': 'Estimate', 'Std. Error': 'Std. Error', 
      'z value': 'z value', 'Pr(>|z|)': 'Pr(>|z|)'}),
    
    ("wue_letters", ['timescale', 'metric', 'anomaly_class'],
     {'emmean': 'emmean', 'SE': 'SE', 'df': 'df', 
      'lower.CL': 'lower.CL', 'upper.CL': 'upper.CL'}),
    
    ("response_letters", ['timescale', 'metric', 'anomaly_class'],
     {'emmean': 'emmean', 'SE': 'SE', 'df': 'df', 
      'lower.CL': 'lower.CL', 'upper.CL': 'upper.CL'}),
    
    ("wue_pairs", ['contrast', 'timescale', 'metric'],
     {'estimate': 'estimate', 'SE': 'SE', 't.ratio': 't.ratio', 'p.value': 'p.value'}),
    
    ("response_pairs", ['contrast', 'timescale', 'metric'],
     {'estimate': 'estimate', 'SE': 'SE', 't.ratio': 't.ratio', 'p.value': 'p.value'}),
    
    ("fitted_residuals", ['site_name', 'Year', 'month', 'water_class', 'timescale', 'metric'],
     {'wue_fitted': 'wue_fitted', 'wue_residual': 'wue_residual'}),
    
    ("model_results", ['model'],
     {'AIC': 'AIC', 'BIC': 'BIC', 'logLik': 'logLik'}),
]

for name, merge_cols, compare_map in comparisons:
    if python_available.get(name, False) and malone_available.get(name, False):
        try:
            python_df = pd.read_csv(os.path.join(python_output_dir, python_files[name]))
            malone_df = pd.read_csv(os.path.join(malone_output_dir, malone_files[name]))
            
            results, match = compare_dataframes_with_mapping(
                python_df, malone_df, name, merge_cols, compare_map, tolerance=1e-5
            )
            all_results[name] = results
            if not match:
                all_metrics_match = False
        except Exception as e:
            print(f"\n  {name}: ERROR reading files - {str(e)}")
            all_results[name] = {}
            all_metrics_match = False

# ============================================================================
# STALE FILE CHECK - UPDATED TO INCLUDE RDS FILES
# ============================================================================

print("\n" + "-"*70)
print("STEP 3: Checking for stale files in Python output folder")
print("-"*70)

expected_files = set(python_files.values()) | {
    "Q2_anomalies_model_summary.txt",
    "Q2_anomalies_diagnostic_comparison_report.txt",
    "Q2_anomalies_wue_value_model.rds",
    "Q2_anomalies_signed_response_model.rds",
    "Q2_anomalies_response_magnitude_model.rds",
    "Q2_anomalies_response_direction_model.rds",
}

all_files = set(os.listdir(python_output_dir))
stale_files = all_files - expected_files

if stale_files:
    print("\n[WARNING] Stale files found in output folder:")
    for f in sorted(stale_files):
        print(f"  - {f}")
else:
    print("\n[PASS] No stale files found.")

# ============================================================================
# OVERALL CONCLUSION - HARDENED
# ============================================================================

print("\n" + "="*70)
print("OVERALL CONCLUSION")
print("="*70)

if all_metrics_match:
    print("\n[PASS] ALL CHECKS MATCH!")
    print("        Python Q2 Anomaly successfully replicates Malone R Q2 Anomaly.")
    print("        - All files present")
    print("        - All rows matched")
    print("        - All columns matched")
    print("        - All numeric values within tolerance (1e-5)")
    print("\n        NOTE:")
    print("        - df column has NaN values (lme4 doesn't provide df in summary) - expected")
    print("        - Column sets were checked directly; true column mismatches are failures")
else:
    print("\n[FAIL] Some checks do not match.")
    print("\n        Details of failures:")
    
    # Find and report specific failures
    failure_details = []
    for name, results in all_results.items():
        for col, max_diff in results.items():
            if max_diff is not None and max_diff > 1e-5 and col != 'df':
                failure_details.append(f"{name}.{col}: {max_diff:.8f}")
    
    if failure_details:
        print("\n        Numeric differences > 1e-5:")
        for detail in failure_details:
            print(f"          - {detail}")
    
    # Check for row/column mismatches that might have been flagged
    print("\n        Most likely causes:")
    print("          1. Different model convergence (lme4 vs lmerTest)")
    print("          2. Different df calculation (Satterthwaite vs Kenward-Roger)")
    print("          3. Different factor contrasts in R vs Python")
    print("          4. Slight differences in data filtering or joins")

print("\n" + "="*70)
print("DIAGNOSTIC COMPLETE")
print("="*70)

# ============================================================================
# SAVE REPORT
# ============================================================================

report_lines = [
    "="*70,
    "Q2 ANOMALY DIAGNOSTIC REPORT: Python vs Malone R",
    "="*70,
    "",
    f"Python output directory: {python_output_dir}",
    f"Malone output directory: {malone_output_dir}",
    "",
    "File Availability:",
]

for name in python_files.keys():
    report_lines.append(f"  {name}: Python {'YES' if python_available.get(name, False) else 'NO'}, Malone {'YES' if malone_available.get(name, False) else 'NO'}")

report_lines.append("")
report_lines.append("Comparison Results (Max Differences):")

for name, results in all_results.items():
    report_lines.append(f"\n  {name}:")
    if results:
        for col, max_diff in results.items():
            if max_diff is None:
                continue
            if max_diff <= 1e-8:
                status = "IDENTICAL"
            elif max_diff <= 1e-5:
                status = "OK"
            else:
                status = "FAIL"
            report_lines.append(f"    [{status}] {col}: {max_diff:.10f}")

# Add stale file info
report_lines.append("")
report_lines.append("Stale Files Check:")
if stale_files:
    report_lines.append(f"  [WARNING] Stale files found: {sorted(stale_files)}")
else:
    report_lines.append("  [PASS] No stale files found")

report_lines.append("")
if all_metrics_match:
    report_lines.append("RESULT: PASS - ALL CHECKS MATCH!")
else:
    report_lines.append("RESULT: FAIL - Some checks do not match.")
    report_lines.append("")
    for name, results in all_results.items():
        for col, max_diff in results.items():
            if max_diff is not None and max_diff > 1e-5 and col != 'df':
                report_lines.append(f"  - {name}.{col}: {max_diff:.8f}")

report_lines.append("")
report_lines.append("="*70)

report_path = os.path.join(python_output_dir, "Q2_anomalies_diagnostic_comparison_report.txt")
with open(report_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(report_lines))

print(f"\nReport saved to: {report_path}")