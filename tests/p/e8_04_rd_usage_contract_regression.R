pkg <- "/Users/cxy/Desktop/lwdid_r/lwdid-r"

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing expected file: %s", path), call. = FALSE)
  }
}

assert_contains <- function(text, needle, label) {
  if (!grepl(needle, text, fixed = TRUE)) {
    stop(sprintf("Missing contract '%s': %s", label, needle), call. = FALSE)
  }
}

assert_absent <- function(text, needle, label) {
  if (grepl(needle, text, fixed = TRUE)) {
    stop(sprintf("Unexpected contract '%s': %s", label, needle), call. = FALSE)
  }
}

devtools::document(pkg = pkg, quiet = TRUE)

lwdid_rd_path <- file.path(pkg, "man", "lwdid.Rd")
assert_file_exists(lwdid_rd_path)
lwdid_rd <- read_text(lwdid_rd_path)
for (arg in c("event_time_range", "df_strategy")) {
  assert_contains(lwdid_rd, sprintf("\\item{%s}{", arg), sprintf("lwdid %s docs", arg))
}

new_result_rd_path <- file.path(pkg, "man", "new_lwdid_result.Rd")
assert_file_exists(new_result_rd_path)
new_result_rd <- read_text(new_result_rd_path)
for (arg in c(
  "n_units", "n_periods", "n_cohorts",
  "att_cohort_agg", "se_cohort_agg", "t_stat_cohort_agg",
  "pvalue_cohort_agg", "ci_cohort_agg"
)) {
  assert_contains(
    new_result_rd,
    sprintf("\\item{%s}{", arg),
    sprintf("new_lwdid_result %s docs", arg)
  )
}

estimate_common_timing_rd_path <- file.path(pkg, "man", "dot-estimate_common_timing.Rd")
assert_file_exists(estimate_common_timing_rd_path)
estimate_common_timing_rd <- read_text(estimate_common_timing_rd_path)
for (arg in c("vce", "cluster_var", "alpha")) {
  assert_absent(
    estimate_common_timing_rd,
    sprintf("\\item{%s}{", arg),
    sprintf(".estimate_common_timing stray %s docs", arg)
  )
}

to_latex_comparison_rd_path <- file.path(pkg, "man", "to_latex_comparison.Rd")
assert_file_exists(to_latex_comparison_rd_path)
to_latex_comparison_rd <- read_text(to_latex_comparison_rd_path)
assert_contains(
  to_latex_comparison_rd,
  "\\name{to_latex_comparison}",
  "to_latex_comparison man page"
)

conditions_rd_path <- file.path(pkg, "man", "lwdid-conditions.Rd")
assert_file_exists(conditions_rd_path)
conditions_rd <- read_text(conditions_rd_path)
for (arg in c(
  "param", "value", "allowed", "call",
  "n", "n_treated", "n_control",
  "min_obs", "max_obs", "n_incomplete_units", "pct_unbalanced",
  "gaps", "ivar", "tvar",
  "column", "available",
  "reps_completed", "error_detail",
  "plot_type", "detail",
  "gvar", "invalid_values",
  "aggregate_type", "method", "vce_type", "treat_var", "n_units",
  "required", "rolling", "cohort", "excluded",
  "unit", "missing_quarters", "pre_quarters", "post_quarters",
  "aggregate", "control_group", "violation",
  "min_cell_size", "max_observed_size", "n_cells", "period",
  "n_extreme", "ps_range", "max_weight",
  "source_function", "condition_number",
  "action_taken", "model", "iterations", "max_iterations", "tolerance",
  "original", "switched_to"
)) {
  assert_contains(
    conditions_rd,
    sprintf("\\item{%s}{", arg),
    sprintf("lwdid-conditions %s docs", arg)
  )
}

cat("PASS: E8-04 Rd usage/documentation contracts satisfied.\n")
