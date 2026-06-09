# Q1 difference type + T:ET + ecosystem outputs

This is the retained streamlined Q1 analysis. It includes ecosystem class as a fixed effect and treats site and month as random effects.

## Core models
- Mixed model: `WUE_difference ~ difference_type * water_class * Trans_ratio_z + (1 | site_name) + (1 | month_f)`
- Near-normal smooth GAM: `WUE_difference ~ difference_type * water_class + s(Trans_ratio, by = difference_type:water_class) + s(site_name, bs = "re") + s(month_f, bs = "re")`
- Ecosystem colors match the Malone_Workflow figures: Upland purple, Freshwater blue, Saline orange.

## Key tables
- `Q1_simple_summary_by_difference_type.csv`
- `Q1_simple_mixed_fixed_effects_with_approx_p.csv`
- `Q1_simple_mixed_likelihood_ratio_tests.csv`
- `Q1_simple_mixed_model_comparison_lrt.csv`
- `Q1_simple_gam_parametric_terms.csv`
- `Q1_simple_gam_smooth_terms.csv`
- `Q1_simple_gam_model_comparison.csv`
- `Q1_simple_mixed_predictions_TET.csv`
- `Q1_simple_gam_predictions_TET.csv`
- `Q1_simple_mixed_fitted_residuals.csv`
- `Q1_simple_gam_fitted_residuals.csv`

## Figures
- `Q1_simple_plot_01_mean_difference_type.png`
- `Q1_simple_plot_02_mixed_TET_predictions.png`
- `Q1_simple_plot_03_gam_TET_predictions.png`
- `Q1_simple_plot_04_gam_model_comparison.png`
- `Q1_simple_plot_05_mixed_residuals.png`
- `Q1_simple_plot_06_gam_residuals.png`
- `Q1_simple_plot_07_manuscript_multipanel.png`
