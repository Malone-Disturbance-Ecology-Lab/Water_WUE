# Q1 simplified GAM workflow
#
# Retained Q1 GAM:
#   WUE_difference ~ difference_type +
#     s(Trans_ratio, by = difference_type) +
#     s(site_name, bs = "re")
#
# The full simplified workflow, tables, figures, and companion mixed model live
# in Q1_simple_difference_TET_models.R. This wrapper keeps the original GAM
# entry point while ensuring only the simplified difference type + T:ET analysis
# runs.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  "Water_WUE/Malone_Workflow/WUE_metric_test_Q1_GAM.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "Q1_simple_difference_TET_models.R"))
