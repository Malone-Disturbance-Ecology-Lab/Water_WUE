# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:03:09 2026

@author: ammar
"""

# =============================================================================
# REVIEW SAVED RESULTS
# Quickly identify:
#   - Significant vs nonsignificant results
#   - Which tests matter
#   - Which findings may be worth discussion/figures
# =============================================================================

import pandas as pd

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\results\SPEI_analysis_results"
)

KRUSKAL_FILE = BASE_DIR + r"\ClimateBiome_Q2_Kruskal.csv"
PAIRWISE_FILE = BASE_DIR + r"\ClimateBiome_Q2_Pairwise.csv"
WILCOXON_FILE = BASE_DIR + r"\ClimateBiome_Q2_Wilcoxon.csv"

# =============================================================================
# LOAD FILES
# =============================================================================

kruskal_df = pd.read_csv(KRUSKAL_FILE)

wilcoxon_df = pd.read_csv(WILCOXON_FILE)

# ---------------------------------------------------------------------
# PAIRWISE FILE MAY BE EMPTY
# ---------------------------------------------------------------------

try:

    pairwise_df = pd.read_csv(PAIRWISE_FILE)

except pd.errors.EmptyDataError:

    print("\nPairwise file is empty.")

    print(
        "No significant omnibus Kruskal tests "
        "triggered pairwise comparisons."
    )

    pairwise_df = pd.DataFrame()

print("="*100)
print("CLIMATE–BIOME Q2 RESULT SUMMARY")
print("="*100)

# =============================================================================
# 1. KRUSKAL RESULTS
# =============================================================================

print("\n" + "="*100)
print("1. KRUSKAL–WALLIS RESULTS")
print("="*100)

kruskal_df['Significant'] = kruskal_df['p_value'] < 0.05

for _, row in kruskal_df.iterrows():

    sig = "SIGNIFICANT" if row['Significant'] else "ns"

    print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

    print(f"  H = {row['H_statistic']:.3f}")

    print(f"  p = {row['p_value']:.4f} --> {sig}")

# =============================================================================
# KRUSKAL SUMMARY
# =============================================================================

n_sig_kw = kruskal_df['Significant'].sum()
n_total_kw = len(kruskal_df)

print("\n" + "-"*60)
print("KRUSKAL SUMMARY")
print("-"*60)

print(f"Significant omnibus tests: {n_sig_kw}/{n_total_kw}")

if n_sig_kw == 0:

    print(
        "\n→ No evidence that climate–biome groups "
        "differ systematically in drought/wet response magnitude."
    )

# =============================================================================
# 2. PAIRWISE RESULTS
# =============================================================================

print("\n" + "="*100)
print("2. PAIRWISE RESULTS")
print("="*100)

if len(pairwise_df) == 0:

    print("\nNo pairwise tests were triggered.")

    print(
        "This occurred because no omnibus Kruskal–Wallis "
        "tests were significant."
    )

else:

    sig_pair = pairwise_df[
        pairwise_df['p_corrected'] < 0.05
    ].copy()

    if len(sig_pair) == 0:

        print("\nNo significant pairwise comparisons.")

    else:

        print("\nSIGNIFICANT PAIRWISE RESULTS:")

        for _, row in sig_pair.iterrows():

            print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

            print(f"  {row['Comparison']}")

            print(
                f"  Corrected p = {row['p_corrected']:.4f}"
            )

            print(
                f"  Median comparison = "
                f"{row['Median_1']:+.1f}% vs "
                f"{row['Median_2']:+.1f}%"
            )

# =============================================================================
# 3. WILCOXON RESULTS
# =============================================================================

print("\n" + "="*100)
print("3. WITHIN-GROUP WILCOXON RESULTS")
print("="*100)

wilcoxon_df['Significant'] = wilcoxon_df['p_value'] < 0.05

# ---------------------------------------------------------------------
# SIGNIFICANT RESULTS
# ---------------------------------------------------------------------

sig_wilcox = wilcoxon_df[
    wilcoxon_df['Significant']
].copy()

if len(sig_wilcox) == 0:

    print("\nNo significant within-group directional responses.")

else:

    print("\nSIGNIFICANT RESULTS:")

    for _, row in sig_wilcox.iterrows():

        print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

        print(f"  Group = {row['Final_Group']}")

        print(
            f"  Median % change = "
            f"{row['Median_Percent_Change']:+.1f}%"
        )

        print(
            f"  p = {row['p_value']:.4f}"
        )

        print(
            f"  n = {row['N_sites']}"
        )

# =============================================================================
# NEAR-SIGNIFICANT RESULTS
# =============================================================================

near_sig = wilcoxon_df[
    (wilcoxon_df['p_value'] >= 0.05) &
    (wilcoxon_df['p_value'] < 0.10)
].copy()

print("\n" + "-"*60)
print("NEAR-SIGNIFICANT RESULTS (0.05–0.10)")
print("-"*60)

if len(near_sig) == 0:

    print("\nNo near-significant results.")

else:

    for _, row in near_sig.iterrows():

        print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

        print(f"  Group = {row['Final_Group']}")

        print(
            f"  Median % change = "
            f"{row['Median_Percent_Change']:+.1f}%"
        )

        print(
            f"  p = {row['p_value']:.4f}"
        )

# =============================================================================
# STRONGEST RESPONSES
# =============================================================================

print("\n" + "="*100)
print("4. STRONGEST MAGNITUDE RESPONSES")
print("="*100)

top_declines = wilcoxon_df.sort_values(
    'Median_Percent_Change'
).head(5)

print("\nLARGEST NEGATIVE RESPONSES:")

for _, row in top_declines.iterrows():

    print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

    print(f"  {row['Final_Group']}")

    print(
        f"  Median % change = "
        f"{row['Median_Percent_Change']:+.1f}%"
    )

    print(
        f"  p = {row['p_value']:.4f}"
    )

top_increases = wilcoxon_df.sort_values(
    'Median_Percent_Change',
    ascending=False
).head(5)

print("\n" + "-"*60)
print("LARGEST POSITIVE RESPONSES")
print("-"*60)

for _, row in top_increases.iterrows():

    print(f"\n{row['Condition']} | {row['SPEI_Timescale']} | {row['Metric']}")

    print(f"  {row['Final_Group']}")

    print(
        f"  Median % change = "
        f"{row['Median_Percent_Change']:+.1f}%"
    )

    print(
        f"  p = {row['p_value']:.4f}"
    )

# =============================================================================
# FINAL INTERPRETATION SUMMARY
# =============================================================================

print("\n" + "="*100)
print("FINAL INTERPRETATION SUMMARY")
print("="*100)

print("\n1. Omnibus climate–biome separation:")

print(f"   Significant tests = {n_sig_kw}/{n_total_kw}")

if n_sig_kw == 0:

    print(
        "   → No strong evidence that climate–biome groups "
        "differ systematically in drought/wet response magnitude."
    )

print("\n2. Pairwise climate–biome differences:")

if len(pairwise_df) == 0:

    print(
        "   → No pairwise group differences detected."
    )

print("\n3. Within-group directional responses:")

if len(sig_wilcox) > 0:

    print(
        "   → Some individual climate–biome systems showed "
        "significant directional WUE responses."
    )

    print(
        "   → These responses were stronger for WUE_T "
        "than WUE_ET."
    )

print("\n4. Overall implication:")

print(
    "   → Baseline ecological organization exists, but "
    "drought-response trajectories remain highly heterogeneous "
    "across climate–biome systems."
)

print("\n" + "="*100)
print("SUMMARY COMPLETE")
print("="*100)