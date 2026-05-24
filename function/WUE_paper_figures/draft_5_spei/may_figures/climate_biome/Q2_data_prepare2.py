# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:53:24 2026

@author: ammar
"""

# =============================================================================
# STEP 2: CLIMATE–BIOME Q2 INFERENTIAL TESTING
# Uses ONLY:
#   - Dry(all)
#   - Wet(all)
# For:
#   - SPEI_6
#   - SPEI_48
#
# Same logic as original Q2 workflow
# =============================================================================

import pandas as pd
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

# =============================================================================
# LOAD DATASET
# =============================================================================

FILE = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\results\SPEI_analysis_results"
    r"\Q2_ClimateBiome_Dataset.csv"
)

OUTPUT_DIR = (
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE"
    r"\data_products\results\SPEI_analysis_results"
)

df = pd.read_csv(FILE)

print("="*100)
print("STEP 2: CLIMATE–BIOME Q2 INFERENTIAL TESTING")
print("="*100)

# =============================================================================
# SETTINGS
# =============================================================================

metrics = ['WUE', 'WUE_tra']

metric_display = {
    'WUE': 'WUE_ET',
    'WUE_tra': 'WUE_T'
}

conditions_keep = [
    ('PASS C', 'Dry (all)'),
    ('PASS C', 'Wet (all)')
]

timescales = ['SPEI_6', 'SPEI_48']

# =============================================================================
# STORAGE LISTS
# =============================================================================

all_site_level_results = []

all_descriptive_results = []

all_kruskal_results = []

all_pairwise_results = []

all_wilcoxon_results = []

# =============================================================================
# MAIN LOOP
# =============================================================================

for pass_name, condition in conditions_keep:

    for timescale in timescales:

        print("\n" + "="*100)
        print(f"{condition} | {timescale}")
        print("="*100)

        subset = df[
            (df['Pass'] == pass_name) &
            (df['Condition'] == condition) &
            (df['SPEI_Timescale'] == timescale)
        ].copy()

        # =========================================================================
        # STRICT SHARED-SITE INTERSECTION
        # =========================================================================

        wue_sites = set(
            subset[
                subset['WUE_Metric'] == 'WUE'
            ]['Site'].unique()
        )

        tra_sites = set(
            subset[
                subset['WUE_Metric'] == 'WUE_tra'
            ]['Site'].unique()
        )

        shared_sites = wue_sites.intersection(tra_sites)

        subset = subset[
            subset['Site'].isin(shared_sites)
        ].copy()

        print(f"\nShared sites retained: {len(shared_sites)}")

        # =========================================================================
        # METRIC LOOP
        # =========================================================================

        for metric in metrics:

            print("\n" + "-"*80)
            print(f"{metric_display[metric]}")
            print("-"*80)

            metric_df = subset[
                subset['WUE_Metric'] == metric
            ].copy()

            # ---------------------------------------------------------------------
            # COLLAPSE TO ONE ROW PER SITE
            # ---------------------------------------------------------------------

            site_level = (
                metric_df
                .groupby(['Site', 'Final_Group'])
                .agg({
                    'Median_%_Change': 'median',
                    'Direction': 'first'
                })
                .reset_index()
            )

            # ---------------------------------------------------------------------
            # ADD METADATA
            # ---------------------------------------------------------------------

            site_level['Condition'] = condition
            site_level['SPEI_Timescale'] = timescale
            site_level['Metric'] = metric_display[metric]

            all_site_level_results.append(site_level.copy())

            # ---------------------------------------------------------------------
            # REMOVE GROUPS WITH n < 3
            # ---------------------------------------------------------------------

            valid_groups = (
                site_level['Final_Group']
                .value_counts()
            )

            valid_groups = valid_groups[
                valid_groups >= 3
            ].index.tolist()

            site_level = site_level[
                site_level['Final_Group'].isin(valid_groups)
            ].copy()

            print("\nGroups retained (n >= 3):")

            for grp in valid_groups:

                n_grp = len(
                    site_level[
                        site_level['Final_Group'] == grp
                    ]
                )

                print(f"  {grp}: {n_grp}")

            # ---------------------------------------------------------------------
            # DESCRIPTIVE STATISTICS
            # ---------------------------------------------------------------------

            print("\n" + "-"*60)
            print("DESCRIPTIVE STATISTICS")
            print("-"*60)

            groups_for_kw = []
            group_names = []

            for grp in valid_groups:

                vals = site_level[
                    site_level['Final_Group'] == grp
                ]['Median_%_Change'].dropna().values

                groups_for_kw.append(vals)
                group_names.append(grp)

                median = np.median(vals)
                q25 = np.percentile(vals, 25)
                q75 = np.percentile(vals, 75)

                mean = np.mean(vals)
                sd = np.std(vals, ddof=1)
                se = sd / np.sqrt(len(vals))

                dir_counts = (
                    site_level[
                        site_level['Final_Group'] == grp
                    ]['Direction']
                    .value_counts()
                )

                inc = dir_counts.get('Increase', 0)
                dec = dir_counts.get('Decrease', 0)
                no_change = dir_counts.get('No change', 0)

                print(f"\n{grp}")

                print(f"  n = {len(vals)}")

                print(
                    f"  Median % change = {median:+.1f}%"
                )

                print(
                    f"  IQR = {q25:+.1f}% to {q75:+.1f}%"
                )

                print(
                    f"  Mean ± SE = {mean:+.1f}% ± {se:.1f}%"
                )

                print(
                    f"  Direction counts: "
                    f"Increase={inc}, "
                    f"Decrease={dec}, "
                    f"No change={no_change}"
                )

                # SAVE DESCRIPTIVE RESULTS

                all_descriptive_results.append({

                    'Condition': condition,
                    'SPEI_Timescale': timescale,
                    'Metric': metric_display[metric],
                    'Final_Group': grp,
                    'N_sites': len(vals),
                    'Median_Percent_Change': median,
                    'Q25': q25,
                    'Q75': q75,
                    'Mean': mean,
                    'SE': se,
                    'Increase_Count': inc,
                    'Decrease_Count': dec,
                    'NoChange_Count': no_change
                })

            # ---------------------------------------------------------------------
            # KRUSKAL-WALLIS
            # ---------------------------------------------------------------------

            if len(groups_for_kw) >= 2:

                h_stat, p_kw = kruskal(*groups_for_kw)

                print("\n" + "-"*60)
                print("KRUSKAL-WALLIS TEST")
                print("-"*60)

                print(f"H = {h_stat:.3f}")
                print(f"p = {p_kw:.6f}")

                all_kruskal_results.append({

                    'Condition': condition,
                    'SPEI_Timescale': timescale,
                    'Metric': metric_display[metric],
                    'H_statistic': h_stat,
                    'p_value': p_kw,
                    'N_groups': len(valid_groups)
                })

                # -------------------------------------------------------------
                # PAIRWISE TESTS
                # -------------------------------------------------------------

                if p_kw < 0.05:

                    print("\n✓ Significant group differences detected")

                    pairwise_results = []

                    for i in range(len(group_names)):
                        for j in range(i+1, len(group_names)):

                            g1 = group_names[i]
                            g2 = group_names[j]

                            vals1 = groups_for_kw[i]
                            vals2 = groups_for_kw[j]

                            u_stat, p_pair = mannwhitneyu(
                                vals1,
                                vals2,
                                alternative='two-sided'
                            )

                            pairwise_results.append({
                                'comparison': f"{g1} vs {g2}",
                                'u_stat': u_stat,
                                'p_raw': p_pair,
                                'median1': np.median(vals1),
                                'median2': np.median(vals2)
                            })

                    raw_p = [x['p_raw'] for x in pairwise_results]

                    reject, p_corr, _, _ = multipletests(
                        raw_p,
                        method='bonferroni'
                    )

                    print("\n" + "-"*60)
                    print("PAIRWISE MANN-WHITNEY TESTS")
                    print("-"*60)

                    for i, res in enumerate(pairwise_results):

                        sig = (
                            '***' if p_corr[i] < 0.001 else
                            '**' if p_corr[i] < 0.01 else
                            '*' if p_corr[i] < 0.05 else
                            'ns'
                        )

                        print(f"\n{res['comparison']}")

                        print(
                            f"  U = {res['u_stat']:.1f}"
                        )

                        print(
                            f"  p_raw = {res['p_raw']:.4f}"
                        )

                        print(
                            f"  p_corr = {p_corr[i]:.4f} ({sig})"
                        )

                        all_pairwise_results.append({

                            'Condition': condition,
                            'SPEI_Timescale': timescale,
                            'Metric': metric_display[metric],
                            'Comparison': res['comparison'],
                            'U_statistic': res['u_stat'],
                            'p_raw': res['p_raw'],
                            'p_corrected': p_corr[i],
                            'Significance': sig,
                            'Median_1': res['median1'],
                            'Median_2': res['median2']
                        })

                else:

                    print("\nNo significant differences among groups")

            # ---------------------------------------------------------------------
            # WITHIN-GROUP WILCOXON TESTS
            # ---------------------------------------------------------------------

            print("\n" + "-"*60)
            print("WITHIN-GROUP WILCOXON TESTS")
            print("-"*60)

            is_dry = 'Dry' in condition

            alternative = 'less' if is_dry else 'greater'

            for grp in valid_groups:

                vals = site_level[
                    site_level['Final_Group'] == grp
                ]['Median_%_Change'].dropna().values

                if len(vals) < 3:
                    continue

                try:

                    if np.all(vals == 0):

                        print(f"\n{grp}")
                        print("  All values are zero")
                        continue

                    result = wilcoxon(
                        vals,
                        alternative=alternative,
                        zero_method='wilcox'
                    )

                    w_stat = result.statistic
                    p_val = result.pvalue

                    sig = (
                        '***' if p_val < 0.001 else
                        '**' if p_val < 0.01 else
                        '*' if p_val < 0.05 else
                        'ns'
                    )

                    median_val = np.median(vals)

                    print(f"\n{grp}")

                    print(
                        f"  Median = {median_val:+.1f}%"
                    )

                    print(
                        f"  W = {w_stat:.1f}"
                    )

                    print(
                        f"  p = {p_val:.4f} ({sig})"
                    )

                    all_wilcoxon_results.append({

                        'Condition': condition,
                        'SPEI_Timescale': timescale,
                        'Metric': metric_display[metric],
                        'Final_Group': grp,
                        'Median_Percent_Change': median_val,
                        'W_statistic': w_stat,
                        'p_value': p_val,
                        'Significance': sig,
                        'N_sites': len(vals)
                    })

                except Exception as e:

                    print(f"\n{grp}")
                    print(f"  Wilcoxon unavailable: {e}")

# =============================================================================
# SAVE RESULTS
# =============================================================================

site_df = pd.concat(
    all_site_level_results,
    ignore_index=True
)

site_df.to_csv(
    OUTPUT_DIR + r"\ClimateBiome_Q2_SiteLevel.csv",
    index=False
)

desc_df = pd.DataFrame(all_descriptive_results)

desc_df.to_csv(
    OUTPUT_DIR + r"\ClimateBiome_Q2_DescriptiveStats.csv",
    index=False
)

kruskal_df = pd.DataFrame(all_kruskal_results)

kruskal_df.to_csv(
    OUTPUT_DIR + r"\ClimateBiome_Q2_Kruskal.csv",
    index=False
)

pairwise_df = pd.DataFrame(all_pairwise_results)

pairwise_df.to_csv(
    OUTPUT_DIR + r"\ClimateBiome_Q2_Pairwise.csv",
    index=False
)

wilcoxon_df = pd.DataFrame(all_wilcoxon_results)

wilcoxon_df.to_csv(
    OUTPUT_DIR + r"\ClimateBiome_Q2_Wilcoxon.csv",
    index=False
)

print("\n" + "="*100)
print("RESULT FILES SAVED")
print("="*100)

print("\nSaved CSV files:")

print("  ClimateBiome_Q2_SiteLevel.csv")
print("  ClimateBiome_Q2_DescriptiveStats.csv")
print("  ClimateBiome_Q2_Kruskal.csv")
print("  ClimateBiome_Q2_Pairwise.csv")
print("  ClimateBiome_Q2_Wilcoxon.csv")

print("\n" + "="*100)
print("STEP 2 COMPLETE")
print("="*100)