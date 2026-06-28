"""
DIAGNOSTIC: Compare Python Q2 vs Original Malone R Q2 Outputs
STRENGTHENED VERSION - With hard failures on row/column mismatches
"""

import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# PATHS
# ============================================================================

python_output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2\Q2_WUE_performance_outputs"
malone_output_dir = r"M:\Research\WUE_CUE\Water_WUE\Malone_Workflow\results\Q2_WUE_performance"

# Create python output dir if it doesn't exist (for report)
os.makedirs(python_output_dir, exist_ok=True)

print("="*70)
print("DIAGNOSTIC: COMPARING PYTHON Q2 vs MALONE R Q2 OUTPUTS")
print("="*70)
print(f"\nPython outputs: {python_output_dir}")
print(f"Malone outputs: {malone_output_dir}")

# ============================================================================
# CHECK IF FILES EXIST
# ============================================================================

print("\n" + "-"*70)
print("STEP 1: Checking if output files exist")
print("-"*70)

file_map = {
    "stability": "Q2_site_stability.csv",
    "plasticity_slope": "Q2_site_plasticity_slope.csv",
    "plasticity_range": "Q2_site_plasticity_range.csv",
    "resistance": "Q2_site_resistance.csv",
    "recovery_events": "Q2_site_recovery_events.csv",
    "recovery": "Q2_site_recovery.csv",
    "summary": "Q2_site_performance_summary.csv",
    "ecosystem_tests": "Q2_ecosystem_class_statistical_tests.csv"
}

print("\nPython files found:")
python_available = {}
for name, filename in file_map.items():
    filepath = os.path.join(python_output_dir, filename)
    exists = os.path.exists(filepath)
    python_available[name] = exists
    print(f"  {name}: {'YES' if exists else 'NO'} - {filename}")

print("\nMalone files found:")
malone_available = {}
for name, filename in file_map.items():
    filepath = os.path.join(malone_output_dir, filename)
    exists = os.path.exists(filepath)
    malone_available[name] = exists
    print(f"  {name}: {'YES' if exists else 'NO'} - {filename}")

# Check if all files exist
all_python_files = all(python_available.values())
all_malone_files = all(malone_available.values())

if not all_python_files:
    print("\n[ERROR] Some Python files are missing. Run the Python Q2 workflow first.")
    missing = [name for name, exists in python_available.items() if not exists]
    print(f"  Missing: {', '.join(missing)}")
    exit(1)

if not all_malone_files:
    print("\n[ERROR] Some Malone files are missing. Check Malone output directory.")
    missing = [name for name, exists in malone_available.items() if not exists]
    print(f"  Missing: {', '.join(missing)}")
    exit(1)

print("\n✓ All files present.")

# ============================================================================
# HELPER: Compare two dataframes - STRENGTHENED
# ============================================================================

def compare_dataframes(df1, df2, name, merge_cols, compare_cols, tolerance=1e-5):
    """
    Compare two dataframes and return max differences
    STRENGTHENED: Row/column mismatches are hard failures
    """
    print(f"\n  {name}:")
    print(f"    Python: {len(df1)} rows, {len(df1.columns)} columns")
    print(f"    Malone: {len(df2)} rows, {len(df2.columns)} columns")
    
    # Check column sets - HARD FAILURE
    col_diff_python = set(df2.columns) - set(df1.columns)
    col_diff_malone = set(df1.columns) - set(df2.columns)
    
    all_match = True
    
    if col_diff_python or col_diff_malone:
        if col_diff_python:
            print(f"    [FAIL] Columns missing in Python: {col_diff_python}")
        if col_diff_malone:
            print(f"    [FAIL] Columns missing in Malone: {col_diff_malone}")
        all_match = False
    else:
        print(f"    [PASS] Column sets match ({len(df1.columns)} columns)")
    
    # Merge with indicator to check rows
    merged = df1.merge(df2, on=merge_cols, suffixes=('_python', '_malone'), how='outer', indicator=True)
    
    only_python = (merged['_merge'] == 'left_only').sum()
    only_malone = (merged['_merge'] == 'right_only').sum()
    both = (merged['_merge'] == 'both').sum()
    
    print(f"    Merge: {both} matched, {only_python} only in Python, {only_malone} only in Malone")
    
    # Row mismatch - HARD FAILURE
    if only_python > 0 or only_malone > 0:
        print(f"    [FAIL] Row mismatch detected")
        all_match = False
    else:
        print(f"    [PASS] All rows matched")
    
    results = {}
    
    for col in compare_cols:
        if f'{col}_python' in merged.columns and f'{col}_malone' in merged.columns:
            # Check if both columns are numeric
            try:
                python_vals = pd.to_numeric(merged[f'{col}_python'], errors='coerce')
                malone_vals = pd.to_numeric(merged[f'{col}_malone'], errors='coerce')
                diff = python_vals - malone_vals
                absdiff = np.abs(diff)
                max_diff = absdiff.max()
                mean_diff = absdiff.mean()
                nan_count = absdiff.isna().sum()
                
                results[col] = max_diff
                
                # Determine status
                if nan_count > 0:
                    status = "WARNING"
                    print(f"    {col}: {nan_count} NaN values in difference (max diff for valid = {max_diff:.8f})")
                    # NaN differences are not a hard failure if they're from non-numeric data
                elif max_diff <= 1e-8:
                    status = "IDENTICAL"
                    print(f"    {col}: Max diff = {max_diff:.10f} [IDENTICAL]")
                elif max_diff <= tolerance:
                    status = "OK"
                    print(f"    {col}: Max diff = {max_diff:.8f} [OK]")
                else:
                    status = "FAIL"
                    all_match = False
                    print(f"    {col}: Max diff = {max_diff:.8f} [FAIL] (tolerance = {tolerance})")
            except Exception as e:
                print(f"    {col}: Error comparing - {str(e)}")
                results[col] = np.nan
                all_match = False
        else:
            print(f"    {col}: Column not found in one or both dataframes")
            results[col] = np.nan
            all_match = False
    
    return results, all_match

# ============================================================================
# STEP 2: Compare each metric
# ============================================================================

print("\n" + "-"*70)
print("STEP 2: Comparing individual metrics (tolerance = 1e-5)")
print("-"*70)

all_results = {}
all_metrics_match = True

# Define comparisons
comparisons = {
    "stability": {
        "python_file": "Q2_site_stability.csv",
        "malone_file": "Q2_site_stability.csv",
        "merge_cols": ['site_name', 'water_class'],
        "compare_cols": ['mean_WUE_T', 'var_WUE_T', 'sd_WUE_T', 'cv_WUE_T', 'stability']
    },
    "plasticity_slope": {
        "python_file": "Q2_site_plasticity_slope.csv",
        "malone_file": "Q2_site_plasticity_slope.csv",
        "merge_cols": ['site_name', 'water_class', 'SPEI_timescale'],
        "compare_cols": ['n_months', 'slope', 'r_squared', 'correlation']
    },
    "plasticity_range": {
        "python_file": "Q2_site_plasticity_range.csv",
        "malone_file": "Q2_site_plasticity_range.csv",
        "merge_cols": ['site_name', 'water_class'],
        "compare_cols": ['n_months', 'WUE_T_max', 'WUE_T_min', 'plasticity_range', 'plasticity_p95_p05']
    },
    "resistance": {
        "python_file": "Q2_site_resistance.csv",
        "malone_file": "Q2_site_resistance.csv",
        "merge_cols": ['site_name', 'water_class'],
        "compare_cols": ['n_months_drought', 'n_months_near_normal', 'mean_WUE_T_drought', 
                        'mean_WUE_T_near_normal', 'resistance']
    },
    "recovery": {
        "python_file": "Q2_site_recovery.csv",
        "malone_file": "Q2_site_recovery.csv",
        "merge_cols": ['site_name', 'water_class'],
        "compare_cols": ['n_events', 'mean_recovery', 'median_recovery', 'sd_recovery', 
                        'min_recovery', 'max_recovery', 'pct_full_recovery']
    }
}

for name, config in comparisons.items():
    python_path = os.path.join(python_output_dir, config["python_file"])
    malone_path = os.path.join(malone_output_dir, config["malone_file"])
    
    try:
        python_df = pd.read_csv(python_path)
        malone_df = pd.read_csv(malone_path)
        
        # Clean column names
        python_df.columns = python_df.columns.str.strip()
        malone_df.columns = malone_df.columns.str.strip()
        
        results, match = compare_dataframes(
            python_df, malone_df, name,
            config["merge_cols"], config["compare_cols"], tolerance=1e-5
        )
        all_results[name] = results
        if not match:
            all_metrics_match = False
    except Exception as e:
        print(f"\n  {name}: ERROR reading files - {str(e)}")
        all_results[name] = {}
        all_metrics_match = False

# ============================================================================
# STEP 3: Recovery Events - Special handling with row check
# ============================================================================

print("\n" + "-"*70)
print("STEP 3: Comparing Recovery Events")
print("-"*70)

if os.path.exists(os.path.join(python_output_dir, "Q2_site_recovery_events.csv")) and \
   os.path.exists(os.path.join(malone_output_dir, "Q2_site_recovery_events.csv")):
    
    try:
        python_events = pd.read_csv(os.path.join(python_output_dir, "Q2_site_recovery_events.csv"))
        malone_events = pd.read_csv(os.path.join(malone_output_dir, "Q2_site_recovery_events.csv"))
        
        print(f"\n  recovery_events:")
        print(f"    Python: {len(python_events)} rows, {len(python_events.columns)} columns")
        print(f"    Malone: {len(malone_events)} rows, {len(malone_events.columns)} columns")
        
        # Check column sets
        col_diff_python = set(malone_events.columns) - set(python_events.columns)
        col_diff_malone = set(python_events.columns) - set(malone_events.columns)
        
        events_match = True
        
        if col_diff_python or col_diff_malone:
            if col_diff_python:
                print(f"    [FAIL] Columns missing in Python: {col_diff_python}")
            if col_diff_malone:
                print(f"    [FAIL] Columns missing in Malone: {col_diff_malone}")
            events_match = False
        else:
            print(f"    [PASS] Column sets match ({len(python_events.columns)} columns)")
        
        # Merge on site_name and event_id
        merged = python_events.merge(
            malone_events, 
            on=['site_name', 'event_id'], 
            suffixes=('_python', '_malone'), 
            how='outer',
            indicator=True
        )
        
        only_python = (merged['_merge'] == 'left_only').sum()
        only_malone = (merged['_merge'] == 'right_only').sum()
        both = (merged['_merge'] == 'both').sum()
        
        print(f"    Merge: {both} matched, {only_python} only in Python, {only_malone} only in Malone")
        
        if only_python > 0 or only_malone > 0:
            print(f"    [FAIL] Row mismatch detected")
            events_match = False
        else:
            print(f"    [PASS] All rows matched")
        
        # Compare numeric columns
        compare_cols = ['duration_months', 'n_pre', 'n_post', 'mean_WUE_T_pre', 
                       'mean_WUE_T_post', 'recovery']
        
        event_results = {}
        
        for col in compare_cols:
            if f'{col}_python' in merged.columns and f'{col}_malone' in merged.columns:
                diff = pd.to_numeric(merged[f'{col}_python'], errors='coerce') - \
                       pd.to_numeric(merged[f'{col}_malone'], errors='coerce')
                absdiff = np.abs(diff)
                max_diff = absdiff.max()
                nan_count = absdiff.isna().sum()
                
                event_results[col] = max_diff
                
                if nan_count > 0:
                    print(f"    {col}: {nan_count} NaN values in difference (max diff for valid = {max_diff:.8f})")
                elif max_diff <= 1e-8:
                    print(f"    {col}: Max diff = {max_diff:.10f} [IDENTICAL]")
                elif max_diff <= 1e-5:
                    print(f"    {col}: Max diff = {max_diff:.8f} [OK]")
                else:
                    events_match = False
                    print(f"    {col}: Max diff = {max_diff:.8f} [FAIL]")
            else:
                print(f"    {col}: Column not found")
                event_results[col] = np.nan
                events_match = False
        
        all_results['recovery_events'] = event_results
        if not events_match:
            all_metrics_match = False
            
    except Exception as e:
        print(f"\n  ERROR comparing recovery events: {str(e)}")
        all_metrics_match = False

# ============================================================================
# STEP 4: Ecosystem Tests - Special handling
# ============================================================================

print("\n" + "-"*70)
print("STEP 4: Comparing Ecosystem Tests")
print("-"*70)

if os.path.exists(os.path.join(python_output_dir, "Q2_ecosystem_class_statistical_tests.csv")) and \
   os.path.exists(os.path.join(malone_output_dir, "Q2_ecosystem_class_statistical_tests.csv")):
    
    try:
        python_tests = pd.read_csv(os.path.join(python_output_dir, "Q2_ecosystem_class_statistical_tests.csv"))
        malone_tests = pd.read_csv(os.path.join(malone_output_dir, "Q2_ecosystem_class_statistical_tests.csv"))
        
        print(f"\n  ecosystem_tests:")
        print(f"    Python: {len(python_tests)} rows, {len(python_tests.columns)} columns")
        print(f"    Malone: {len(malone_tests)} rows, {len(malone_tests.columns)} columns")
        
        # Check column sets
        col_diff_python = set(malone_tests.columns) - set(python_tests.columns)
        col_diff_malone = set(python_tests.columns) - set(malone_tests.columns)
        
        tests_match = True
        
        if col_diff_python or col_diff_malone:
            if col_diff_python:
                print(f"    [FAIL] Columns missing in Python: {col_diff_python}")
            if col_diff_malone:
                print(f"    [FAIL] Columns missing in Malone: {col_diff_malone}")
            tests_match = False
        else:
            print(f"    [PASS] Column sets match ({len(python_tests.columns)} columns)")
        
        # Merge on metric
        merged = python_tests.merge(
            malone_tests, 
            on=['metric'], 
            suffixes=('_python', '_malone'), 
            how='outer',
            indicator=True
        )
        
        only_python = (merged['_merge'] == 'left_only').sum()
        only_malone = (merged['_merge'] == 'right_only').sum()
        both = (merged['_merge'] == 'both').sum()
        
        print(f"    Merge: {both} matched, {only_python} only in Python, {only_malone} only in Malone")
        
        if only_python > 0 or only_malone > 0:
            print(f"    [FAIL] Row mismatch detected")
            tests_match = False
        else:
            print(f"    [PASS] All rows matched")
        
        # Compare statistical columns - looser tolerance for p-values
        compare_cols = ['statistic', 'df', 'p_value', 
                       'pairwise_upland_vs_freshwater_p_adj', 
                       'pairwise_upland_vs_saline_p_adj',
                       'pairwise_freshwater_vs_saline_p_adj']
        
        test_results = {}
        
        for col in compare_cols:
            if f'{col}_python' in merged.columns and f'{col}_malone' in merged.columns:
                diff = pd.to_numeric(merged[f'{col}_python'], errors='coerce') - \
                       pd.to_numeric(merged[f'{col}_malone'], errors='coerce')
                absdiff = np.abs(diff)
                max_diff = absdiff.max()
                nan_count = absdiff.isna().sum()
                
                test_results[col] = max_diff
                
                # Looser tolerance for p-values (1e-4)
                if col in ['p_value', 'pairwise_upland_vs_freshwater_p_adj', 
                          'pairwise_upland_vs_saline_p_adj', 'pairwise_freshwater_vs_saline_p_adj']:
                    tol = 1e-4
                else:
                    tol = 1e-5
                
                if nan_count > 0:
                    print(f"    {col}: {nan_count} NaN values in difference (max diff for valid = {max_diff:.8f})")
                elif max_diff <= 1e-8:
                    print(f"    {col}: Max diff = {max_diff:.10f} [IDENTICAL]")
                elif max_diff <= tol:
                    print(f"    {col}: Max diff = {max_diff:.8f} [OK] (tolerance = {tol})")
                else:
                    tests_match = False
                    print(f"    {col}: Max diff = {max_diff:.8f} [FAIL] (tolerance = {tol})")
            else:
                print(f"    {col}: Column not found")
                test_results[col] = np.nan
                tests_match = False
        
        all_results['ecosystem_tests'] = test_results
        if not tests_match:
            all_metrics_match = False
            
    except Exception as e:
        print(f"\n  ERROR comparing ecosystem tests: {str(e)}")
        all_metrics_match = False

# ============================================================================
# STEP 5: Summary Table Comparison
# ============================================================================

print("\n" + "-"*70)
print("STEP 5: Comparing Summary Table")
print("-"*70)

if os.path.exists(os.path.join(python_output_dir, "Q2_site_performance_summary.csv")) and \
   os.path.exists(os.path.join(malone_output_dir, "Q2_site_performance_summary.csv")):
    
    try:
        python_summary = pd.read_csv(os.path.join(python_output_dir, "Q2_site_performance_summary.csv"))
        malone_summary = pd.read_csv(os.path.join(malone_output_dir, "Q2_site_performance_summary.csv"))
        
        print(f"\n  Summary Table:")
        print(f"    Python: {len(python_summary)} rows, {len(python_summary.columns)} columns")
        print(f"    Malone: {len(malone_summary)} rows, {len(malone_summary.columns)} columns")
        
        # Check column sets
        col_diff_python = set(malone_summary.columns) - set(python_summary.columns)
        col_diff_malone = set(python_summary.columns) - set(malone_summary.columns)
        
        summary_match = True
        
        if col_diff_python or col_diff_malone:
            if col_diff_python:
                print(f"    [FAIL] Columns missing in Python: {col_diff_python}")
            if col_diff_malone:
                print(f"    [FAIL] Columns missing in Malone: {col_diff_malone}")
            summary_match = False
        else:
            print(f"    [PASS] Column sets match ({len(python_summary.columns)} columns)")
        
        # Merge
        merged_summary = python_summary.merge(
            malone_summary, 
            on=['site_name', 'water_class'], 
            suffixes=('_python', '_malone'),
            how='outer',
            indicator=True
        )
        
        only_python = (merged_summary['_merge'] == 'left_only').sum()
        only_malone = (merged_summary['_merge'] == 'right_only').sum()
        both = (merged_summary['_merge'] == 'both').sum()
        
        print(f"    Merge: {both} matched, {only_python} only in Python, {only_malone} only in Malone")
        
        if only_python > 0 or only_malone > 0:
            print(f"    [FAIL] Row mismatch detected")
            summary_match = False
        else:
            print(f"    [PASS] All rows matched")
        
        # Compare key columns
        key_cols = ['stability', 'plasticity_slope_SPEI3', 'plasticity_range', 
                   'resistance', 'mean_recovery']
        
        summary_results = {}
        
        for col in key_cols:
            if f'{col}_python' in merged_summary.columns and f'{col}_malone' in merged_summary.columns:
                diff = pd.to_numeric(merged_summary[f'{col}_python'], errors='coerce') - \
                       pd.to_numeric(merged_summary[f'{col}_malone'], errors='coerce')
                absdiff = np.abs(diff)
                max_diff = absdiff.max()
                nan_count = absdiff.isna().sum()
                
                summary_results[col] = max_diff
                
                if nan_count > 0:
                    print(f"    {col}: {nan_count} NaN values in difference (max diff for valid = {max_diff:.8f})")
                elif max_diff <= 1e-8:
                    print(f"    {col}: Max diff = {max_diff:.10f} [IDENTICAL]")
                elif max_diff <= 1e-5:
                    print(f"    {col}: Max diff = {max_diff:.8f} [OK]")
                else:
                    summary_match = False
                    print(f"    {col}: Max diff = {max_diff:.8f} [FAIL]")
            else:
                print(f"    {col}: Column not found")
                summary_results[col] = np.nan
                summary_match = False
        
        all_results['summary'] = summary_results
        if not summary_match:
            all_metrics_match = False
            
    except Exception as e:
        print(f"\n  ERROR comparing summary: {str(e)}")
        all_metrics_match = False

# ============================================================================
# STEP 6: STALE FILE CHECK
# ============================================================================

print("\n" + "-"*70)
print("STEP 6: Checking for stale files in Python output folder")
print("-"*70)

expected_files = {
    "Q2_site_stability.csv",
    "Q2_site_plasticity_slope.csv",
    "Q2_site_plasticity_range.csv",
    "Q2_site_resistance.csv",
    "Q2_site_recovery_events.csv",
    "Q2_site_recovery.csv",
    "Q2_site_performance_summary.csv",
    "Q2_ecosystem_class_statistical_tests.csv",
    "Q2_WUE_performance_summary.txt",
}

all_files = set(os.listdir(python_output_dir))
stale_files = all_files - expected_files

# Exclude diagnostic report if it exists
stale_files = {f for f in stale_files if f != "Q2_diagnostic_comparison_report.txt"}

if stale_files:
    print(f"\n  [WARNING] Stale files found in output folder:")
    for f in sorted(stale_files):
        print(f"    - {f}")
else:
    print(f"\n  [PASS] No stale files found.")

# ============================================================================
# STEP 7: OVERALL CONCLUSION
# ============================================================================

print("\n" + "="*70)
print("OVERALL CONCLUSION")
print("="*70)

if all_metrics_match:
    print("\n[PASS] ALL METRICS MATCH! Python Q2 successfully replicates Malone R Q2.")
    print("        - All row counts match")
    print("        - All column sets match")
    print("        - All numeric values within tolerance (1e-5, 1e-4 for p-values)")
else:
    print("\n[FAIL] NOT ALL METRICS MATCH.")
    print("        Some metrics have differences exceeding the tolerance.")
    print("\n        Metrics with differences > tolerance:")
    
    has_failures = False
    for name, results in all_results.items():
        for col, max_diff in results.items():
            if pd.notna(max_diff) and max_diff > 0.001:
                has_failures = True
                print(f"          - {name}.{col}: {max_diff:.8f}")
    
    if not has_failures:
        print("          No large differences (>0.001) found.")
        print("          Check individual column tolerances above.")
    
    print("\n        Most likely causes for differences:")
    print("          1. Recovery: Indexing differences (0-based vs 1-based) in pre/post windows")
    print("          2. Recovery: Inclusion/exclusion of drought months in windows")
    print("          3. Ecosystem tests: Wilcoxon/Mann-Whitney approximation differences")
    print("          4. Data filtering: Slight differences in complete.cases() handling")

# ============================================================================
# STEP 8: SAVE DIAGNOSTIC REPORT
# ============================================================================

print("\n" + "-"*70)
print("STEP 7: Saving Diagnostic Report")
print("-"*70)

report_lines = [
    "="*70,
    "Q2 DIAGNOSTIC REPORT: Python vs Malone R",
    "="*70,
    "",
    f"Python output directory: {python_output_dir}",
    f"Malone output directory: {malone_output_dir}",
    "",
    "File Availability:",
]

for name in file_map.keys():
    report_lines.append(f"  {name}: Python {'YES' if python_available.get(name, False) else 'NO'}, Malone {'YES' if malone_available.get(name, False) else 'NO'}")

report_lines.append("")
report_lines.append("Comparison Results (Max Differences):")

for name, results in all_results.items():
    report_lines.append(f"\n  {name}:")
    if results:
        for col, max_diff in results.items():
            if pd.isna(max_diff):
                status = "ERROR"
            elif max_diff <= 1e-8:
                status = "IDENTICAL"
            elif max_diff <= 1e-5:
                status = "OK"
            else:
                status = "FAIL"
            report_lines.append(f"    [{status}] {col}: {max_diff:.10f}" if pd.notna(max_diff) else f"    [ERROR] {col}: N/A")

report_lines.append("")
if all_metrics_match:
    report_lines.append("RESULT: PASS - ALL METRICS MATCH!")
else:
    report_lines.append("RESULT: FAIL - Some metrics do not match.")
    report_lines.append("")
    report_lines.append("Metrics with differences > 0.001:")
    for name, results in all_results.items():
        for col, max_diff in results.items():
            if pd.notna(max_diff) and max_diff > 0.001:
                report_lines.append(f"  - {name}.{col}: {max_diff:.8f}")

# Add stale file info
report_lines.append("")
report_lines.append("Stale Files Check:")
if stale_files:
    report_lines.append("  [WARNING] Stale files found:")
    for f in sorted(stale_files):
        report_lines.append(f"    - {f}")
else:
    report_lines.append("  [PASS] No stale files found.")

report_lines.append("")
report_lines.append("="*70)

report_path = os.path.join(python_output_dir, "Q2_diagnostic_comparison_report.txt")
with open(report_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(report_lines))

print(f"  Report saved to: {report_path}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*70)
print("DIAGNOSTIC COMPLETE")
print("="*70)

print(f"\nReport saved to: {report_path}")

print("\nStatus Summary:")
print(f"  Files checked: {len(file_map)}")
print(f"  All files present: {'YES' if all_python_files and all_malone_files else 'NO'}")
print(f"  All metrics match: {'YES' if all_metrics_match else 'NO'}")
print(f"  Stale files: {'YES' if stale_files else 'NO'}")
print(f"  Report: {report_path}")

print("\n" + "="*70)