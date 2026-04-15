suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This audit must be run via Rscript --file.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
artifact_dir <- dirname(script_path)
helper_locator <- file.path(
  dirname(dirname(dirname(script_path))),
  "helper-package-source-paths.R"
)
source(helper_locator, local = TRUE)
package_dir <- resolve_package_source_root()
helper_path <- resolve_package_source_file(
  "tests",
  "testthat",
  "helper-visualization-export.R"
)
legacy_probe_path <- file.path(
  artifact_dir,
  "e8_04_package_hardening_regression.R"
)
json_path <- file.path(
  artifact_dir,
  "20260324-qa-parity-e8-04-layer5-release-regression.json"
)
md_path <- file.path(
  artifact_dir,
  "20260324-qa-parity-e8-04-layer5-release-regression.md"
)

run_rscript <- function(expr_or_path, use_expression = FALSE) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  args <- if (isTRUE(use_expression)) c("-e", shQuote(expr_or_path)) else expr_or_path
  output <- system2(rscript_bin, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  list(
    status = as.integer(status),
    output = output
  )
}

legacy_probe <- run_rscript(legacy_probe_path)

visualization_expr <- sprintf(
  "devtools::test(pkg=%s, filter=%s, reporter=%s)",
  dQuote(package_dir),
  dQuote("visualization-export"),
  dQuote("summary")
)
visualization_probe <- run_rscript(visualization_expr, use_expression = TRUE)

extract_suite_section <- function(lines, start_pattern, stop_patterns) {
  start_idx <- grep(start_pattern, lines)
  if (length(start_idx) == 0L) {
    return(character())
  }
  start_idx <- start_idx[[1L]]
  remaining <- lines[(start_idx + 1L):length(lines)]
  if (length(remaining) == 0L) {
    return(character())
  }

  stop_idx <- integer()
  for (pattern in stop_patterns) {
    hits <- grep(pattern, remaining)
    if (length(hits) > 0L) {
      stop_idx <- c(stop_idx, hits[[1L]] - 1L)
    }
  }

  if (length(stop_idx) == 0L) {
    remaining
  } else {
    remaining[seq_len(max(0L, min(stop_idx)))]
  }
}

suite_warning_lines <- extract_suite_section(
  visualization_probe$output,
  "^== Warnings ",
  c("^== Failed ", "^== DONE ")
)
suite_warning_lines <- suite_warning_lines[nzchar(suite_warning_lines)]
suite_failure_lines <- grep("^-- ", visualization_probe$output, value = TRUE)
suite_failed <- length(suite_failure_lines) > 0L

source(helper_path, local = TRUE)
devtools::load_all(package_dir, quiet = TRUE)

obj_common <- .mock_common_timing_result()
csv_path <- tempfile(fileext = ".csv")
on.exit(unlink(csv_path), add = TRUE)

by_cohort_error <- tryCatch(
  {
    to_csv(obj_common, file = csv_path, what = "by_cohort")
    NA_character_
  },
  error = function(e) conditionMessage(e)
)

expected_by_cohort_literal <- "No cohort-specific results (Staggered mode only)."
by_cohort_literal_match <- !is.na(by_cohort_error) &&
  grepl(expected_by_cohort_literal, by_cohort_error, fixed = TRUE)
by_cohort_contains_u_escape <- !is.na(by_cohort_error) &&
  grepl("<U\\+[0-9A-F]{4}>", by_cohort_error)

obj_ri <- .mock_result_with_ri()
ri_summary_lines <- capture.output(print(summary(obj_ri)))
ri_summary_text <- paste(ri_summary_lines, collapse = "\n")
ri_summary_header_present <- grepl(
  "Randomization Inference",
  ri_summary_text,
  fixed = TRUE
)
ri_summary_valid_token_present <- grepl(
  "valid=999",
  ri_summary_text,
  fixed = TRUE
)
ri_summary_valid_fraction_present <- grepl(
  "valid=999/999",
  ri_summary_text,
  fixed = TRUE
)
ri_detail_line <- grep("valid=", ri_summary_lines, value = TRUE)
if (length(ri_detail_line) == 0L) {
  ri_detail_line <- NA_character_
} else {
  ri_detail_line <- ri_detail_line[[1L]]
}

failures <- character()
if (legacy_probe$status != 0L) {
  failures <- c(
    failures,
    sprintf(
      "[legacy-hardening] %s",
      paste(tail(legacy_probe$output, 5L), collapse = " | ")
    )
  )
}
if (suite_failed) {
  failures <- c(
    failures,
    sprintf(
      "[visualization-suite] %s",
      paste(suite_failure_lines, collapse = " | ")
    )
  )
}
if (!by_cohort_literal_match) {
  failures <- c(
    failures,
    sprintf(
      "[by-cohort-export] expected literal '%s' not found in error message: %s",
      expected_by_cohort_literal,
      by_cohort_error
    )
  )
}
if (!ri_summary_valid_fraction_present) {
  failures <- c(
    failures,
    sprintf(
      "[ri-summary] expected 'valid=999/999' but observed: %s",
      ri_detail_line
    )
  )
}

payload <- list(
  story = "story-E8-04",
  layer = "layer_5",
  role = "qa-parity",
  generated_at = format(
    Sys.time(),
    tz = "Asia/Macau",
    format = "%Y-%m-%d %H:%M:%S CST"
  ),
  legacy_hardening_probe = list(
    script = basename(legacy_probe_path),
    status = if (legacy_probe$status == 0L) "passed" else "failed",
    output_tail = tail(legacy_probe$output, 10L)
  ),
  visualization_export_suite = list(
    status = if (suite_failed) "failed" else "passed",
    failure_labels = unname(suite_failure_lines),
    warning_lines = unname(suite_warning_lines)
  ),
  direct_probes = list(
    by_cohort_export = list(
      expected_literal = expected_by_cohort_literal,
      actual_message = by_cohort_error,
      literal_match = by_cohort_literal_match,
      contains_u_escape = by_cohort_contains_u_escape
    ),
    ri_summary = list(
      header_present = ri_summary_header_present,
      valid_token_present = ri_summary_valid_token_present,
      valid_fraction_present = ri_summary_valid_fraction_present,
      detail_line = ri_detail_line
    )
  ),
  exact_status = if (length(failures) == 0L) "passed" else "failed",
  numeric_status = "not-applicable",
  blocker_boundary = if (length(failures) == 0L) {
    "targeted-release-contracts-cleared-full-check-pending"
  } else if (legacy_probe$status == 0L) {
    "legacy-ac17-parse-cleared-live-blocker-shifted-to-export-summary-contracts"
  } else {
    "legacy-ac17-parse-still-live"
  },
  remaining_gap = if (length(failures) == 0L) {
    "none"
  } else {
    "visualization-export suite and Layer 5 release contracts still fail"
  }
)

writeLines(
  toJSON(
    payload,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = 16
  ),
  json_path
)

md_lines <- c(
  "# story-E8-04 Layer 5 release regression",
  "",
  "## 时间",
  sprintf("- `%s`", payload$generated_at),
  "",
  "## 范围",
  sprintf(
    "- 重跑 `%s`，确认旧的 `Rd/nobs/parse` hardening probe 是否仍失败。",
    basename(legacy_probe_path)
  ),
  "- 运行 `devtools::test(filter = \"visualization-export\")` 作为当前 release blocker suite。",
  "- 直接复现 `to_csv(..., what = \"by_cohort\")` 与 RI summary print contract。",
  "",
  "## 关键结果",
  sprintf(
    "- legacy hardening probe: `%s`。",
    payload$legacy_hardening_probe$status
  ),
  sprintf("- visualization-export suite: `%s`。", payload$visualization_export_suite$status),
  sprintf(
    "- suite failure labels: `%s`。",
    if (length(suite_failure_lines) == 0L) "none" else paste(suite_failure_lines, collapse = "`；`")
  ),
  sprintf(
    "- by_cohort export literal match: `%s`；actual message: `%s`。",
    if (isTRUE(by_cohort_literal_match)) "true" else "false",
    by_cohort_error
  ),
  sprintf(
    "- RI summary header present: `%s`；valid fraction present: `%s`；detail line: `%s`。",
    if (isTRUE(ri_summary_header_present)) "true" else "false",
    if (isTRUE(ri_summary_valid_fraction_present)) "true" else "false",
    ri_detail_line
  ),
  "",
  "## 产出文件",
  sprintf("- JSON: `%s`", json_path),
  sprintf("- 说明: `%s`", md_path),
  "",
  "## 裁决",
  if (legacy_probe$status == 0L) {
    "- 旧的 `aggregate_to_overall.Rd` / `compute_inference.Rd` / `test-sensitivity-pre-period.R` / namespace `nobs` blocker 已退出 live blocker 集。"
  } else {
    "- 旧的 AC-17 parse / namespace blocker 仍未完全退出 live blocker 集。"
  },
  if (length(failures) == 0L) {
    "- 当前 targeted Layer 5 release probe 已转绿：`visualization-export` suite、`by_cohort` export literal contract 与 RI summary `valid=999/999` contract 当前均通过。"
  } else {
    "- 当前 Layer 5 live blocker 已收敛为 `visualization-export` suite 两条 contract："
  },
  if (length(failures) == 0L) {
    "- 后续若 `story-E8-04` 仍未 closure-ready，剩余 blocker 应由 fresh `devtools::check(--no-manual)` 的 test / Rd / documentation hygiene 结果单独维护，不得继续把 `TC-10.6.18` / `TC-10.6.34` 写成 live blocker。"
  } else {
    "  `TC-10.6.18` 的 by-cohort export error message 在当前会话下未被 literal regex 稳定捕获；"
  },
  if (length(failures) == 0L) {
    "- 该证据只覆盖 targeted release blocker probe；完整 `devtools::check(--no-manual)` 仍需单独作为 closure-readiness 证据维护。"
  } else {
    "  `TC-10.6.34` 的 RI summary 仍未打印 `valid=999/999`。"
  },
  if (length(failures) > 0L) {
    "- 该证据只覆盖 targeted release blocker probe；完整 `devtools::check(--no-manual)` 仍需单独作为 closure-readiness 证据维护。"
  }
)

writeLines(md_lines, md_path)

if (length(failures) > 0L) {
  stop(
    paste(
      c("Layer 5 release regression failed:", failures),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("Layer 5 release regression passed.")
