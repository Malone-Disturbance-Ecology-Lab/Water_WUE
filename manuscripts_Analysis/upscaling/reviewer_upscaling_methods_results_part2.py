
#please make sure to run the first part of the script which is upscaling_methods_results



# ============================================================
# 1B. UNCERTAINTY AROUND MEAN SPEI-3 BY COAST
# ============================================================
print_section("1B. UNCERTAINTY AROUND MEAN SPEI-3 BY COAST")

from scipy import stats

BOOT_N = 10000
RNG = np.random.default_rng(42)

uncertainty_rows = []

for coast in COAST_ORDER:

    sub = (
        monthly.loc[
            monthly["coast_region"] == coast,
            ["year", "month", "mean_SPEI3"]
        ]
        .dropna()
        .copy()
    )

    x = sub["mean_SPEI3"].to_numpy()

    n_months = len(x)
    n_years = sub["year"].nunique()

    # --------------------------------------------------------
    # Basic descriptive uncertainty
    # --------------------------------------------------------
    mean_val = np.mean(x)
    sd_val = np.std(x, ddof=1)
    se_val = sd_val / np.sqrt(n_months)

    # Conventional 95% CI treating months as independent
    tcrit = stats.t.ppf(0.975, df=n_months - 1)
    naive_ci_low = mean_val - tcrit * se_val
    naive_ci_high = mean_val + tcrit * se_val

    # --------------------------------------------------------
    # YEAR-BLOCK BOOTSTRAP
    # Resample whole years rather than individual months.
    # This keeps April-October observations from the same year
    # together and is preferable for the long-term mean.
    # --------------------------------------------------------
    years = np.sort(sub["year"].unique())

    boot_means = np.empty(BOOT_N)

    for b in range(BOOT_N):
        sampled_years = RNG.choice(
            years,
            size=len(years),
            replace=True
        )

        sampled_values = []

        for yr in sampled_years:
            vals = sub.loc[
                sub["year"] == yr,
                "mean_SPEI3"
            ].to_numpy()

            sampled_values.extend(vals)

        boot_means[b] = np.mean(sampled_values)

    block_ci_low, block_ci_high = np.percentile(
        boot_means,
        [2.5, 97.5]
    )

    uncertainty_rows.append({
        "coast_region": coast,
        "n_years": n_years,
        "n_months": n_months,
        "mean_SPEI3": mean_val,
        "SD_monthly": sd_val,
        "SE_monthly": se_val,
        "naive_95CI_low": naive_ci_low,
        "naive_95CI_high": naive_ci_high,
        "year_block_95CI_low": block_ci_low,
        "year_block_95CI_high": block_ci_high,
    })


uncertainty = pd.DataFrame(uncertainty_rows)
uncertainty = order_coasts(uncertainty)

print(
    uncertainty.round(3).to_string(index=False)
)


# ============================================================
# SENTENCE-READY OUTPUT
# ============================================================
print("\nSentence-ready values:\n")

for _, r in uncertainty.iterrows():

    coast = r["coast_region"]
    label = DISPLAY.get(coast, coast)

    print(
        f"{label}: mean SPEI-3 = "
        f"{r['mean_SPEI3']:.2f} "
        f"(95% year-block bootstrap CI: "
        f"{r['year_block_95CI_low']:.2f} to "
        f"{r['year_block_95CI_high']:.2f}); "
        f"monthly SD = {r['SD_monthly']:.2f}; "
        f"monthly SE = {r['SE_monthly']:.2f}."
    )