# Q1 simplified WUE metric-difference workflow
#
# Retained Q1 model:
#   WUE_difference ~ difference_type * Trans_ratio_z + (1 | site_name)
#
# The full simplified workflow, tables, figures, and GAM companion model live in
# Q1_simple_difference_TET_models.R. This wrapper keeps the original Q1 entry
# point while ensuring only the simplified difference type + T:ET analysis runs.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  "Water_WUE/Malone_Workflow/WUE_metric_test_Q1.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "Q1_simple_difference_TET_models.R"))
