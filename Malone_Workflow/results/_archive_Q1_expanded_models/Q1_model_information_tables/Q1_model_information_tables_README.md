# Q1 model information tables

Generated from existing Q1 mixed-model and GAM outputs.

## Mixed model
- `Q1_mixed_primary_fixed_effects_with_approx_p.csv`: fixed effects from primary `lmer`; p-values are normal approximations from t values.
- `Q1_mixed_primary_likelihood_ratio_tests.csv`: term-level likelihood-ratio tests from `drop1(primary_model)`.
- `Q1_mixed_model_comparison_lrt.csv`: additive vs T:ET ecosystem vs primary model comparison.
- `Q1_mixed_single_contrast_fixed_effects_with_approx_p.csv`: fixed effects for each single metric contrast; p-values are normal approximations.
- `Q1_mixed_single_contrast_likelihood_ratio_tests.csv`: LRTs for each single metric contrast.

## GAM
- `Q1_gam_parametric_terms.csv`: smooth GAM parametric terms with p-values.
- `Q1_gam_smooth_terms.csv`: smooth terms with edf, F, and p-values.
- `Q1_gam_model_comparison.csv`: linear vs smooth GAM fit statistics.
- `Q1_gam_ecosystem_pairwise_contrasts.csv`: ecosystem pairwise contrasts at representative T:ET levels.
