import pandas as pd
import numpy as np
import os
from scipy.stats import kruskal

# ------------------------------------------------------------------------------
# PATHS
# ------------------------------------------------------------------------------
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q2"
input_dir = os.path.join(base_dir, "Q2_WUE_performance_outputs", "Q_performace metric_T_ET_relationship")
output_dir = os.path.join(base_dir, "Q2_WUE_performance_figures")
os.makedirs(output_dir, exist_ok=True)

site_data_file = os.path.join(input_dir, "Q2_performance_by_coast_site_level_data.csv")
corr_summary_file = os.path.join(input_dir, "Q2_T_ET_performance_relationship_summary.csv")
coast_tests_file = os.path.join(input_dir, "Q2_performance_by_coast_tests.csv")
table_output = os.path.join(output_dir, "Table_S2_Q2_summary.csv")

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------
df_sites = pd.read_csv(site_data_file)
df_corr = pd.read_csv(corr_summary_file)
df_tests = pd.read_csv(coast_tests_file)

# ------------------------------------------------------------------------------
# METRICS
# ------------------------------------------------------------------------------
metrics = [
    "stability",
    "plasticity_slope_SPEI3",
    "plasticity_p95_p05",
    "resistance",
    "mean_recovery"
]
display_names = {
    "stability": "Stability",
    "plasticity_slope_SPEI3": "Plasticity slope (SPEI-3)",
    "plasticity_p95_p05": "Plasticity range",
    "resistance": "Drought resistance (SPEI-3)",
    "mean_recovery": "Post-drought recovery"
}

# ------------------------------------------------------------------------------
# HELPER FOR P‑VALUE FORMATTING
# ------------------------------------------------------------------------------
def fmt_p_text(p):
    """Format p-value: '< 0.001' if very small, otherwise '= 0.xxx'."""
    if np.isnan(p):
        return "--"
    if p < 0.001:
        return "< 0.001"
    return f"= {p:.3f}"

# ------------------------------------------------------------------------------
# BUILD TABLE ROWS
# ------------------------------------------------------------------------------
rows = []
for metric in metrics:
    vals = df_sites[metric].dropna()
    n = len(vals)
    med = vals.median() if n > 0 else np.nan
    minv = vals.min() if n > 0 else np.nan
    maxv = vals.max() if n > 0 else np.nan

    # Spearman from correlation summary
    corr_row = df_corr[df_corr["metric"] == metric]
    if not corr_row.empty:
        rho = corr_row["spearman_rho"].values[0]
        p_corr = corr_row["spearman_p_FDR"].values[0]
    else:
        rho, p_corr = np.nan, np.nan

    # Coast Kruskal-Wallis p (FDR-adjusted) from coast tests file
    test_row = df_tests[df_tests["metric"] == metric]
    if not test_row.empty:
        p_coast = test_row["kw_p_all_groups_FDR"].values[0]
    else:
        p_coast = np.nan

    # Ecosystem Kruskal-Wallis p (raw, not FDR-adjusted across metrics)
    eco_groups = []
    for eco in ["Upland", "Freshwater", "Saline"]:
        eco_vals = df_sites[df_sites["water_class"] == eco][metric].dropna()
        if len(eco_vals) > 0:
            eco_groups.append(eco_vals)
    if len(eco_groups) >= 2:
        h_stat, p_eco = kruskal(*eco_groups)
    else:
        p_eco = np.nan

    rows.append({
        "metric_display": display_names[metric],
        "n": n,
        "median": med,
        "min": minv,
        "max": maxv,
        "rho": rho,
        "p_corr": p_corr,
        "p_eco": p_eco,
        "p_coast": p_coast
    })

df_table = pd.DataFrame(rows)

# ------------------------------------------------------------------------------
# FORMAT COLUMNS
# ------------------------------------------------------------------------------
df_table["Overall median (min-max)"] = df_table.apply(
    lambda r: f"{r['median']:.2f} ({r['min']:.2f}-{r['max']:.2f})"
    if not np.isnan(r['median']) else "--", axis=1
)

df_table["Association with mean T:ET, ρ (FDR p)"] = df_table.apply(
    lambda r: f"ρ = {r['rho']:.2f}, p {fmt_p_text(r['p_corr'])}"
    if not np.isnan(r['rho']) else "--", axis=1
)

df_table["Ecosystem class p-value"] = df_table["p_eco"].apply(
    lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}" if not np.isnan(p) else "--"
)
df_table["Coastal region p-value (FDR)"] = df_table["p_coast"].apply(
    lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}" if not np.isnan(p) else "--"
)

# ------------------------------------------------------------------------------
# FINAL SELECTION AND SORTING
# ------------------------------------------------------------------------------
final_columns = [
    "metric_display",
    "n",
    "Overall median (min-max)",
    "Association with mean T:ET, ρ (FDR p)",
    "Ecosystem class p-value",
    "Coastal region p-value (FDR)"
]
df_final = df_table[final_columns]
df_final.rename(columns={"metric_display": "Metric"}, inplace=True)

# Sort by the order in the metrics list
df_final["Metric"] = pd.Categorical(df_final["Metric"], categories=[display_names[m] for m in metrics], ordered=True)
df_final = df_final.sort_values("Metric").reset_index(drop=True)

# ------------------------------------------------------------------------------
# SAVE WITH UTF-8-SIG (BOM) TO AVOID ENCODING ISSUES IN EXCEL
# ------------------------------------------------------------------------------
df_final.to_csv(table_output, index=False, encoding='utf-8-sig')
print(f"✅ Q2 summary table saved to: {table_output}")

# Print to console for quick inspection
print("\n" + "="*80)
print("Q2 SUMMARY TABLE")
print("="*80)
print(df_final.to_string(index=False))