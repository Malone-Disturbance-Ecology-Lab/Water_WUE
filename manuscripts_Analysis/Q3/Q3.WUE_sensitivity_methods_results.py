#!/usr/bin/env python3
"""
Q3_RESULT_CHECKER.py
Read-only Q3 result-checking script.
Reads the updated CSVs and prints model interpretation to the console.
Does NOT write any files. Does NOT overwrite Malone's official summary.
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy.stats import kruskal
from datetime import datetime

# ----------------------------------------------------------------------------
# PATHS (adjust if necessary)
# ----------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
output_dir = os.path.join(base_dir, "Q3_WUE_T_SPEI_sensitivity_outputs")

# List of required input files (only these are read)
required_files = [
    "Q3_WUE_T_SPEI_model_base_wide.csv",
    "Q3_WUE_T_SPEI_model_data_long.csv",
    "Q3_WUE_T_SPEI_model_comparison.csv",
    "Q3_WUE_T_SPEI_smooth_terms.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_smooth_terms.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_threshold_markers_5_10_20pct_upland.csv",
    "Q3_WUE_T_SPEI_coast_threshold_GAM_predicted_impact_classes_upland.csv",
    "Q3_WUE_T_SPEI_coastline_grouped_site_sensitivity_summary.csv",
    "Q3_WUE_T_SPEI_site_sensitivity_rank.csv",
]

# Check for missing files
missing = []
for fname in required_files:
    if not os.path.exists(os.path.join(output_dir, fname)):
        missing.append(fname)
if missing:
    print("ERROR: The following required files are missing:")
    for f in missing:
        print(f"  - {f}")
    print("Please ensure the Q3 workflow has been run successfully.")
    sys.exit(1)

print("="*80)
print("Q3 RESULT CHECKER (read-only) - FROM UPDATED CSVs")
print(f"Checked: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------
def get_required_col(df, possible_names):
    """Return the first column name from possible_names that exists in df."""
    for name in possible_names:
        if name in df.columns:
            return name
    raise ValueError(
        f"Missing column. Tried {possible_names}. "
        f"Available columns: {list(df.columns)}"
    )

def parse_smooth_term(term, prefix):
    """
    Extract timescale and group (ecosystem or coast) from a smooth term.
    Example: "s(SPEI_value):spei_ecosystemSPEI_3__Upland"
    Returns (timescale, group) or (nan, nan).
    """
    if prefix in term:
        parts = term.split(prefix)[-1]
        if "__" in parts:
            timescale, group = parts.split("__", 1)
            return timescale, group
    return np.nan, np.nan

# ----------------------------------------------------------------------------
# LOAD REQUIRED CSVs (all exist)
# ----------------------------------------------------------------------------
print("\nLoading CSVs...")
model_base = pd.read_csv(os.path.join(output_dir, required_files[0]))
model_long = pd.read_csv(os.path.join(output_dir, required_files[1]))
model_comp = pd.read_csv(os.path.join(output_dir, required_files[2]))
eco_smooth = pd.read_csv(os.path.join(output_dir, required_files[3]))
coast_smooth = pd.read_csv(os.path.join(output_dir, required_files[4]))
coast_threshold_markers = pd.read_csv(os.path.join(output_dir, required_files[5]))
# coast_pred_impact is not used directly in this script, but we keep it loaded if needed.
coast_pred_impact = pd.read_csv(os.path.join(output_dir, required_files[6]))
coastline_summary = pd.read_csv(os.path.join(output_dir, required_files[7]))
site_sensitivity = pd.read_csv(os.path.join(output_dir, required_files[8]))

print("All CSVs loaded successfully.\n")

# ----------------------------------------------------------------------------
# 1. Dataset size and coverage
# ----------------------------------------------------------------------------
n_obs = len(model_base)
n_sites = model_base['site_name'].nunique()
year_min = model_base['Year'].min()
year_max = model_base['Year'].max()
n_years = year_max - year_min + 1
n_long = len(model_long)
n_timescales = model_long['SPEI_timescale'].nunique()

# ----------------------------------------------------------------------------
# 2. Ecosystem distribution (observations and sites)
# ----------------------------------------------------------------------------
eco_dist = model_base.groupby('water_class').agg(
    obs=('WUE_T', 'count'),
    sites=('site_name', 'nunique')
).reset_index()
eco_dist_dict = {row['water_class']: {'obs': row['obs'], 'sites': row['sites']}
                 for _, row in eco_dist.iterrows()}

# ----------------------------------------------------------------------------
# 3. Coast distribution (observations and sites)
# ----------------------------------------------------------------------------
coast_dist = model_base.groupby('coast_region').agg(
    obs=('WUE_T', 'count'),
    sites=('site_name', 'nunique')
).reset_index()
coast_dist_dict = {row['coast_region']: {'obs': row['obs'], 'sites': row['sites']}
                   for _, row in coast_dist.iterrows()}

# ----------------------------------------------------------------------------
# 4. Model comparison
# ----------------------------------------------------------------------------
eco_row = model_comp[model_comp['model'] == 'ecosystem_smooth_gam']
coast_row = model_comp[model_comp['model'] == 'coast_threshold_gam']
if len(eco_row) == 0 or len(coast_row) == 0:
    raise ValueError("Model comparison CSV missing expected model names.")
eco_comp = eco_row.iloc[0]
coast_comp = coast_row.iloc[0]

aic_eco = eco_comp['AIC']
aic_coast = coast_comp['AIC']
bic_eco = eco_comp['BIC']
bic_coast = coast_comp['BIC']
loglik_eco = eco_comp['logLik']
loglik_coast = coast_comp['logLik']
dev_eco = eco_comp['deviance_explained']
dev_coast = coast_comp['deviance_explained']
r2_eco = eco_comp['adjusted_r_squared']
r2_coast = coast_comp['adjusted_r_squared']

delta_aic = aic_coast - aic_eco   # positive if coast is worse (higher AIC)
if delta_aic < -2:
    model_message = "The coast-threshold GAM was AIC-favored over the ecosystem GAM."
elif delta_aic > 2:
    model_message = "The ecosystem GAM was AIC-favored over the coast-threshold GAM."
else:
    model_message = "The two GAMs had similar AIC support."

# ----------------------------------------------------------------------------
# 5. Ecosystem smooth terms: parse, significance, strongest, Upland check
# ----------------------------------------------------------------------------
p_col_eco = get_required_col(eco_smooth, ["p-value", "p.value", "p"])
f_col_eco = get_required_col(eco_smooth, ["F"])

eco_smooth[["timescale", "ecosystem"]] = eco_smooth["smooth_term"].apply(
    lambda x: pd.Series(parse_smooth_term(x, "spei_ecosystem"))
)
eco_smooth["significant"] = eco_smooth[p_col_eco] < 0.05
eco_smooth_sig = eco_smooth[eco_smooth["significant"]].copy()

# Strongest ecosystem smooth
idx_eco_max = eco_smooth[f_col_eco].idxmax()
strongest_eco = eco_smooth.loc[idx_eco_max]

# Upland significance across all seven timescales
all_timescales = set(["SPEI_1", "SPEI_3", "SPEI_6", "SPEI_12", "SPEI_24", "SPEI_36", "SPEI_48"])
upland_sig_timescales = set(
    eco_smooth.loc[
        (eco_smooth["ecosystem"] == "Upland") & (eco_smooth["significant"]),
        "timescale"
    ].dropna()
)
all_upland_sig = upland_sig_timescales == all_timescales

# ----------------------------------------------------------------------------
# 6. Coast-threshold smooth terms: parse, significance, strongest, Gulf status
# ----------------------------------------------------------------------------
p_col_coast = get_required_col(coast_smooth, ["p-value", "p.value", "p"])
f_col_coast = get_required_col(coast_smooth, ["F"])

coast_smooth[["timescale", "coast"]] = coast_smooth["smooth_term"].apply(
    lambda x: pd.Series(parse_smooth_term(x, "spei_coast"))
)
coast_smooth["significant"] = coast_smooth[p_col_coast] < 0.05
coast_smooth_sig = coast_smooth[coast_smooth["significant"]].copy()

# Strongest coast smooth
idx_coast_max = coast_smooth[f_col_coast].idxmax()
strongest_coast = coast_smooth.loc[idx_coast_max]

gulf_sig = coast_smooth_sig[coast_smooth_sig["coast"] == "Gulf Coast"]
if len(gulf_sig) == 0:
    gulf_status = "nonsignificant (no smooth terms with p<0.05)"
else:
    gulf_status = f"significant for {gulf_sig['timescale'].nunique()} timescales"

# ----------------------------------------------------------------------------
# 7. Threshold markers (5%, 10%, 20%)
# ----------------------------------------------------------------------------
threshold_markers = coast_threshold_markers.copy()
if 'threshold_pct' in threshold_markers.columns:
    if threshold_markers['threshold_pct'].dtype in ['int64', 'float64']:
        threshold_markers['threshold_pct'] = threshold_markers['threshold_pct'].astype(int).astype(str) + '%'

# ----------------------------------------------------------------------------
# 8. Observed coastline sensitivity medians (site-level)
# ----------------------------------------------------------------------------
selected_timescales = ["SPEI_1", "SPEI_3", "SPEI_48"]
sens_filtered = site_sensitivity[
    (site_sensitivity['SPEI_timescale'].isin(selected_timescales)) &
    (site_sensitivity['n_months'] >= 6)
]
coast_sens_medians = sens_filtered.groupby(['coast_region', 'SPEI_timescale'])['mean_abs_pct_change'].median().unstack()
# Reindex to ensure all three timescales appear, even if missing
coast_sens_medians = coast_sens_medians.reindex(columns=selected_timescales)

# Ranking by mean median across the three timescales
if not coast_sens_medians.empty:
    coast_rank = coast_sens_medians.mean(axis=1).sort_values(ascending=False)
    highest_coast = coast_rank.index[0]
    second_coast = coast_rank.index[1] if len(coast_rank) > 1 else None
    lowest_coast = coast_rank.index[-1]
else:
    highest_coast = second_coast = lowest_coast = None

# ----------------------------------------------------------------------------
# 9. Kruskal-Wallis tests for coast differences (for each timescale)
# ----------------------------------------------------------------------------
kw_results = {}
for ts in selected_timescales:
    data_groups = []
    coast_names = []
    for coast in coast_sens_medians.index:
        vals = sens_filtered[
            (sens_filtered['coast_region'] == coast) &
            (sens_filtered['SPEI_timescale'] == ts)
        ]['mean_abs_pct_change'].dropna().values
        if len(vals) > 0:
            data_groups.append(vals)
            coast_names.append(coast)
    if len(data_groups) >= 2:
        h, p = kruskal(*data_groups)
        kw_results[ts] = {'H': h, 'p': p, 'n_groups': len(data_groups)}
    else:
        kw_results[ts] = {'H': np.nan, 'p': np.nan, 'n_groups': len(data_groups)}

# ----------------------------------------------------------------------------
# 10. Warning/Check for Gulf markers if it is not the strongest coast smooth
# ----------------------------------------------------------------------------
if strongest_coast['coast'] != 'Gulf Coast':
    gulf_markers = threshold_markers[threshold_markers['coast_region'] == 'Gulf Coast']
    if len(gulf_markers) > 0:
        print("\nCHECK: Gulf Coast still has threshold markers, but it is not the strongest coast smooth.")
        print("         Verify whether these are real updated main-model threshold markers.")
    else:
        print("\nNo Gulf markers found; consistency okay.")
else:
    print("\nStrongest coast smooth is Gulf Coast, so markers are consistent.")

# ----------------------------------------------------------------------------
# 11. Assemble report lines (to be printed, not saved)
# ----------------------------------------------------------------------------
report_lines = []
report_lines.append("="*80)
report_lines.append("Q3 RESULT CHECK - SUMMARY FROM UPDATED CSVs")
report_lines.append("="*80)

# Dataset
report_lines.append("\nDATASET:")
report_lines.append(f"  • {n_obs:,} site-month observations")
report_lines.append(f"  • {n_sites} sites")
report_lines.append(f"  • {year_min} to {year_max} ({n_years} years)")
report_lines.append(f"  • {n_long:,} long-format (site-month-timescale) rows")
report_lines.append(f"  • {n_timescales} SPEI timescales")

report_lines.append("\nEcosystem distribution:")
for eco, stats in eco_dist_dict.items():
    report_lines.append(f"  • {eco}: {stats['obs']:,} obs, {stats['sites']} sites")

report_lines.append("\nCoast distribution:")
for coast, stats in coast_dist_dict.items():
    report_lines.append(f"  • {coast}: {stats['obs']:,} obs, {stats['sites']} sites")

# Model comparison
report_lines.append("\nMODEL COMPARISON:")
report_lines.append(f"  Ecosystem GAM: AIC={aic_eco:.2f}, Dev={dev_eco*100:.1f}%, R²={r2_eco:.3f}")
report_lines.append(f"  Coast GAM:     AIC={aic_coast:.2f}, Dev={dev_coast*100:.1f}%, R²={r2_coast:.3f}")
report_lines.append(f"  Delta AIC (Coast - Eco) = {delta_aic:.2f}  -> {model_message}")

# Ecosystem smooths
report_lines.append("\nECOSYSTEM SMOOTH TERMS (p<0.05):")
report_lines.append(f"  • {len(eco_smooth_sig)} significant out of {len(eco_smooth)} total.")
report_lines.append(f"  • Strongest: {strongest_eco['smooth_term']} (F={strongest_eco[f_col_eco]:.2f}, p={strongest_eco[p_col_eco]:.3g})")
if all_upland_sig:
    report_lines.append("  • Upland: significant across all seven SPEI timescales.")
else:
    missing = all_timescales - upland_sig_timescales
    report_lines.append(f"  • Upland: significant for {len(upland_sig_timescales)} of 7 timescales; missing: {', '.join(sorted(missing))}.")

# Coast smooths
report_lines.append("\nCOAST-THRESHOLD SMOOTH TERMS (p<0.05):")
report_lines.append(f"  • {len(coast_smooth_sig)} significant out of {len(coast_smooth)} total.")
if not pd.isna(strongest_coast['coast']):
    report_lines.append(
        f"  • Strongest: {strongest_coast['coast']} {strongest_coast['timescale']} "
        f"(F={strongest_coast[f_col_coast]:.2f}, p={strongest_coast[p_col_coast]:.3g})"
    )
else:
    report_lines.append("  • No coast-specific smooth term parsed.")
report_lines.append(f"  • Gulf Coast: {gulf_status}")

# Threshold markers summary
report_lines.append("\nTHRESHOLD MARKERS (5%, 10%, 20%):")
if not threshold_markers.empty:
    # Count per coast
    for coast in coast_dist_dict.keys():
        sub = threshold_markers[threshold_markers['coast_region'] == coast]
        if not sub.empty:
            report_lines.append(f"  • {coast}: {len(sub)} markers")
            # Optionally list a few details
            for _, row in sub.iterrows():
                report_lines.append(
                    f"      {row['SPEI_timescale']} {row['anomaly_side']} {row['impact_direction']} at {row['threshold_pct']}: "
                    f"SPEI {row['SPEI_threshold']:.2f} (range {row['SPEI_threshold_lower']:.2f}-{row['SPEI_threshold_upper']:.2f})"
                )
else:
    report_lines.append("  No threshold markers found.")

# Observed sensitivity medians
report_lines.append("\nOBSERVED SITE-LEVEL SENSITIVITY (median of mean abs % change, n_months>=6):")
if not coast_sens_medians.empty:
    report_lines.append(f"{'Coast':<15} {'SPEI-1':>10} {'SPEI-3':>10} {'SPEI-48':>10}")
    report_lines.append("-"*45)
    for coast in coast_sens_medians.index:
        vals = coast_sens_medians.loc[coast]
        report_lines.append(
            f"{coast:<15} {vals['SPEI_1']:>10.1f} {vals['SPEI_3']:>10.1f} {vals['SPEI_48']:>10.1f}"
        )
    # Ranking
    if highest_coast is not None:
        report_lines.append(f"\n  Ranking (mean of medians):")
        for i, (coast, mean_val) in enumerate(coast_rank.items(), 1):
            report_lines.append(f"    {i}. {coast}: {mean_val:.1f}%")
else:
    report_lines.append("  Insufficient data for medians.")

# Kruskal-Wallis
report_lines.append("\nKRUSKAL-WALLIS TESTS (coast differences):")
for ts, res in kw_results.items():
    if not np.isnan(res['p']):
        report_lines.append(f"  • {ts}: H={res['H']:.2f}, p={res['p']:.3g} (n_groups={res['n_groups']})")
    else:
        report_lines.append(f"  • {ts}: insufficient data for test.")

# ----------------------------------------------------------------------------
# 12. Manuscript interpretation (conditional, not hardcoded)
# ----------------------------------------------------------------------------
report_lines.append("\nMANUSCRIPT-READY INTERPRETATION:")
# Model selection
if "ecosystem GAM was AIC-favored" in model_message:
    report_lines.append("  • The ecosystem GAM was AIC-favored, indicating ecosystem class captured more structure.")
elif "coast-threshold GAM was AIC-favored" in model_message:
    report_lines.append("  • The coast-threshold GAM was AIC-favored, suggesting coast-specific responses improve fit.")
else:
    report_lines.append("  • Models had similar AIC support; coast-specific effects not strongly distinguished.")

# Upland
if all_upland_sig:
    report_lines.append("  • Upland smooths significant across all SPEI timescales, indicating broad moisture sensitivity.")
else:
    report_lines.append(f"  • Upland smooths significant for {len(upland_sig_timescales)} of 7 timescales.")

# Strongest coast
if not pd.isna(strongest_coast['coast']):
    report_lines.append(f"  • Strongest coast-specific smooth: {strongest_coast['coast']} {strongest_coast['timescale']}.")
else:
    report_lines.append("  • No dominant coast smooth.")

# Gulf status
if "nonsignificant" in gulf_status:
    report_lines.append("  • Gulf Coast smooths were largely nonsignificant after updated classification.")
else:
    report_lines.append(f"  • Gulf Coast smooths significant for {gulf_sig['timescale'].nunique()} timescales.")

# Observed sensitivity
if highest_coast is not None:
    report_lines.append(f"  • Observed sensitivity highest in {highest_coast} (mean median {coast_rank.iloc[0]:.1f}%).")
    if second_coast is not None:
        report_lines.append(f"  • Second highest: {second_coast} ({coast_rank.iloc[1]:.1f}%).")

# Kruskal-Wallis significance
any_sig = any(res['p'] < 0.05 for res in kw_results.values() if not np.isnan(res['p']))
if any_sig:
    report_lines.append("  • Kruskal-Wallis tests indicated significant coast differences (p<0.05) for at least one timescale.")
else:
    report_lines.append("  • Kruskal-Wallis tests did not detect significant coast differences (p≥0.05).")

# Final note
report_lines.append("\n  • Coast-specific SPEI responses were evaluated as a secondary framework; ecosystem-level Upland remained stable.")
report_lines.append("  • Use threshold markers with uncertainty; rerun sensitivity diagnostics if classification changes.")

# ----------------------------------------------------------------------------
# 13. Print report (no file writing)
# ----------------------------------------------------------------------------
print("\n" + '\n'.join(report_lines))
print("\n" + "="*80)
print("Q3 result check complete. No files were written.")
print("Official model summary remains: Q3_WUE_T_SPEI_model_summary.txt")
print("="*80)