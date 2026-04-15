suppressPackageStartupMessages({
  library(devtools)
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This audit must be run via Rscript --file.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
repo_root <- dirname(dirname(dirname(dirname(script_path))))
package_dir <- file.path(repo_root, "lwdid-r")
helper_path <- file.path(
  package_dir,
  "tests",
  "testthat",
  "helper-epic7-dgp.R"
)
json_path <- file.path(
  dirname(script_path),
  "20260324-qa-parity-e8-03-fatal-inheritance-regression.json"
)
md_path <- file.path(
  dirname(script_path),
  "20260324-qa-parity-e8-03-fatal-inheritance-regression.md"
)

source(helper_path, local = TRUE)
load_all(package_dir, quiet = TRUE)

tolerance <- 1e-10

collect_warnings <- function(expr) {
  warnings <- list()
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, list(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    result = result,
    warnings = warnings
  )
}

dt <- generate_staggered_panel(
  n_per_cohort = 8L,
  cohorts = c(4L, 6L, 8L),
  n_never_treated = 6L,
  T_total = 10L,
  tau_base = 2.0,
  tau_dynamic = 0.5,
  seed = 42L
)

strict_wrapper <- suppressWarnings(suppressMessages(
  lwdid_sensitivity(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    type = "no_anticipation",
    control_group = "not_yet_treated",
    aggregate = "none",
    verbose = FALSE
  )
))

direct_nyt <- suppressWarnings(suppressMessages(
  lwdid(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ra",
    control_group = "not_yet_treated",
    aggregate = "none",
    verbose = "quiet"
  )
))

direct_nt <- suppressWarnings(suppressMessages(
  lwdid(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ra",
    control_group = "never_treated",
    aggregate = "none",
    verbose = "quiet"
  )
))

aggregation_capture <- collect_warnings(
  lwdid_sensitivity(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    type = "transformation",
    control_group = "not_yet_treated",
    aggregate = "overall",
    verbose = FALSE
  )
)

aggregation_wrapper <- aggregation_capture$result
aggregation_warning_messages <- vapply(
  aggregation_capture$warnings,
  conditionMessage,
  character(1)
)
aggregation_warning_classes <- lapply(aggregation_capture$warnings, class)
switch_warning_detected <- any(vapply(
  aggregation_capture$warnings,
  function(w) inherits(w, "lwdid_control_group_switch"),
  logical(1)
))

direct_overall_demean <- suppressWarnings(suppressMessages(
  lwdid(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ra",
    control_group = "never_treated",
    aggregate = "overall",
    verbose = "quiet"
  )
))

direct_overall_detrend <- suppressWarnings(suppressMessages(
  lwdid(
    data = dt,
    y = "y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "detrend",
    estimator = "ra",
    control_group = "never_treated",
    aggregate = "overall",
    verbose = "quiet"
  )
))

strict_wrapper_att <- strict_wrapper$baseline_estimate$att
aggregation_comp <- aggregation_wrapper$transformation_comparison

strict_matches_nyt <- isTRUE(all.equal(
  strict_wrapper_att,
  direct_nyt$att,
  tolerance = tolerance
))
strict_differs_from_nt <- abs(strict_wrapper_att - direct_nt$att) > 1e-6

aggregation_matches_nt <- all(
  isTRUE(all.equal(
    aggregation_comp$demean_att,
    direct_overall_demean$att,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    aggregation_comp$demean_se,
    direct_overall_demean$se_att,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    aggregation_comp$detrend_att,
    direct_overall_detrend$att,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    aggregation_comp$detrend_se,
    direct_overall_detrend$se_att,
    tolerance = tolerance
  ))
)

exact_status <- if (strict_matches_nyt &&
                    strict_differs_from_nt &&
                    switch_warning_detected &&
                    aggregation_matches_nt) {
  "passed"
} else {
  "failed"
}

numeric_status <- if (strict_matches_nyt &&
                      strict_differs_from_nt &&
                      aggregation_matches_nt) {
  "passed"
} else {
  "failed"
}

payload <- list(
  story = "story-E8-03",
  task = "8.5",
  role = "qa-parity",
  generated_at = format(
    Sys.time(),
    tz = "Asia/Macau",
    format = "%Y-%m-%d %H:%M:%S CST"
  ),
  audit_seed = 42L,
  tolerance = tolerance,
  strict_mask_case = list(
    wrapper_baseline_att = unname(strict_wrapper_att),
    direct_not_yet_treated_att = unname(direct_nyt$att),
    direct_never_treated_att = unname(direct_nt$att),
    abs_diff_vs_direct_nyt = unname(abs(strict_wrapper_att - direct_nyt$att)),
    abs_diff_vs_direct_nt = unname(abs(strict_wrapper_att - direct_nt$att)),
    wrapper_matches_direct_nyt = strict_matches_nyt,
    wrapper_differs_from_direct_nt = strict_differs_from_nt
  ),
  aggregation_case = list(
    switch_warning_detected = switch_warning_detected,
    warning_messages = unname(aggregation_warning_messages),
    warning_classes = aggregation_warning_classes,
    wrapper_demean_att = unname(aggregation_comp$demean_att),
    wrapper_demean_se = unname(aggregation_comp$demean_se),
    wrapper_detrend_att = unname(aggregation_comp$detrend_att),
    wrapper_detrend_se = unname(aggregation_comp$detrend_se),
    direct_never_treated_demean_att = unname(direct_overall_demean$att),
    direct_never_treated_demean_se = unname(direct_overall_demean$se_att),
    direct_never_treated_detrend_att = unname(direct_overall_detrend$att),
    direct_never_treated_detrend_se = unname(direct_overall_detrend$se_att),
    wrapper_matches_direct_nt = aggregation_matches_nt
  ),
  exact_status = exact_status,
  numeric_status = numeric_status,
  blocker_boundary = "narrowed",
  remaining_gap = "FATAL-003 story-level RI regression"
)

writeLines(
  toJSON(
    payload,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = 16
  ),
  json_path,
  useBytes = TRUE
)

md_lines <- c(
  "# story-E8-03 FATAL inheritance regression",
  "",
  "## 时间",
  "",
  paste0("- `", payload$generated_at, "`"),
  "",
  "## 范围",
  "",
  "- 使用 `generate_staggered_panel(seed = 42)` 对 `story-E8-03` comprehensive sensitivity wrapper 做 Layer 5 FATAL inheritance regression。",
  "- 本轮覆盖 `FATAL-001` 的 `not_yet_treated` strict-mask 数值路径，以及 `FATAL-002` / `FATAL-004` 的 aggregation auto-switch + NT weighted-average 路径。",
  "- `FATAL-003` 仍未形成 story-level `ri=TRUE` executable regression，本证据只把 blocker 缩窄到该项。",
  "",
  "## 关键结果",
  "",
  paste0(
    "- strict-mask baseline：wrapper = `",
    format(payload$strict_mask_case$wrapper_baseline_att, digits = 15),
    "`；direct NYT = `",
    format(payload$strict_mask_case$direct_not_yet_treated_att, digits = 15),
    "`；direct NT = `",
    format(payload$strict_mask_case$direct_never_treated_att, digits = 15),
    "`。"
  ),
  paste0(
    "- strict-mask 判定：`wrapper_matches_direct_nyt = ",
    tolower(as.character(payload$strict_mask_case$wrapper_matches_direct_nyt)),
    "`；`wrapper_differs_from_direct_nt = ",
    tolower(as.character(payload$strict_mask_case$wrapper_differs_from_direct_nt)),
    "`。"
  ),
  paste0(
    "- aggregation 判定：`switch_warning_detected = ",
    tolower(as.character(payload$aggregation_case$switch_warning_detected)),
    "`；`wrapper_matches_direct_nt = ",
    tolower(as.character(payload$aggregation_case$wrapper_matches_direct_nt)),
    "`。"
  ),
  paste0(
    "- overall/demean：wrapper ATT/SE = `",
    format(payload$aggregation_case$wrapper_demean_att, digits = 15),
    "` / `",
    format(payload$aggregation_case$wrapper_demean_se, digits = 15),
    "`。"
  ),
  paste0(
    "- overall/detrend：wrapper ATT/SE = `",
    format(payload$aggregation_case$wrapper_detrend_att, digits = 15),
    "` / `",
    format(payload$aggregation_case$wrapper_detrend_se, digits = 15),
    "`。"
  ),
  "",
  "## 产出文件",
  "",
  paste0("- JSON oracle: `", json_path, "`"),
  paste0("- 当前说明: `", md_path, "`"),
  "",
  "## 裁决",
  "",
  paste0("- `exact_status = ", payload$exact_status, "`"),
  paste0("- `numeric_status = ", payload$numeric_status, "`"),
  paste0("- `remaining_gap = ", payload$remaining_gap, "`")
)

writeLines(md_lines, md_path, useBytes = TRUE)

cat(json_path, "\n")
cat(md_path, "\n")
