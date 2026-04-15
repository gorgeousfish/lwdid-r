suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This probe must be run via Rscript --file.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
artifact_dir <- dirname(script_path)
helper_locator <- file.path(
  dirname(dirname(dirname(script_path))),
  "helper-package-source-paths.R"
)
source(helper_locator, local = TRUE)
package_dir <- resolve_package_source_root()

json_path <- file.path(
  artifact_dir,
  "20260324-qa-parity-e8-04-layer5-interactive-bucket.json"
)
md_path <- file.path(
  artifact_dir,
  "20260324-qa-parity-e8-04-layer5-interactive-bucket.md"
)

run_rscript <- function(expr_or_path, use_expression = FALSE) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  args <- if (isTRUE(use_expression)) {
    c("-e", shQuote(expr_or_path, type = "sh"))
  } else {
    expr_or_path
  }
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

extract_section <- function(lines, start_pattern, stop_patterns) {
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

parse_test_probe <- function(probe) {
  failure_lines <- grep("^-- ", probe$output, value = TRUE)
  warning_lines <- extract_section(
    probe$output,
    "^== Warnings ",
    c("^== Failed ", "^== DONE ")
  )
  warning_lines <- warning_lines[nzchar(warning_lines)]

  list(
    status_code = probe$status,
    status = if (probe$status == 0L && length(failure_lines) == 0L) {
      "passed"
    } else {
      "failed"
    },
    failure_labels = unname(failure_lines),
    warning_lines = unname(warning_lines),
    output_tail = tail(probe$output, 12L)
  )
}

make_test_expr <- function(filter) {
  sprintf(
    "devtools::test(pkg=%s, filter=%s, reporter=%s)",
    dQuote(package_dir),
    dQuote(filter),
    dQuote("summary")
  )
}

release_probe <- run_rscript(file.path(artifact_dir, "e8_04_layer5_release_regression.R"))
transform_probe <- parse_test_probe(
  run_rscript(
    make_test_expr("transform-detrend|utils|validate"),
    use_expression = TRUE
  )
)
vce_probe <- parse_test_probe(
  run_rscript(
    make_test_expr("vce-return-values|vce-integration"),
    use_expression = TRUE
  )
)
visualization_probe <- parse_test_probe(
  run_rscript(
    make_test_expr("visualization-export"),
    use_expression = TRUE
  )
)

partial_check_expr <- sprintf(
  paste0(
    "devtools::check(",
    "pkg=%s, ",
    "args=c(%s, %s, %s), ",
    "error_on=%s, ",
    "document=FALSE",
    ")"
  ),
  dQuote(package_dir),
  dQuote("--no-manual"),
  dQuote("--no-tests"),
  dQuote("--no-examples"),
  dQuote("never")
)
partial_check_raw <- run_rscript(partial_check_expr, use_expression = TRUE)
partial_check_warning_lines <- grep(" WARNING$", partial_check_raw$output, value = TRUE)
partial_check_note_lines <- grep(" NOTE$", partial_check_raw$output, value = TRUE)
partial_check_error_lines <- grep(" ERROR$", partial_check_raw$output, value = TRUE)

partial_check_probe <- list(
  status = if (length(partial_check_error_lines) == 0L) {
    "passed-no-errors"
  } else {
    "failed"
  },
  warning_lines = unname(partial_check_warning_lines),
  note_lines = unname(partial_check_note_lines),
  error_lines = unname(partial_check_error_lines),
  output_tail = tail(partial_check_raw$output, 20L)
)

failures <- character()
if (release_probe$status != 0L) {
  failures <- c(
    failures,
    sprintf(
      "[release-regression] %s",
      paste(tail(release_probe$output, 6L), collapse = " | ")
    )
  )
}

for (probe_name in c("transform_probe", "vce_probe", "visualization_probe")) {
  probe <- get(probe_name, inherits = FALSE)
  if (!identical(probe$status, "passed")) {
    detail <- c(
      probe$failure_labels,
      if (length(probe$failure_labels) == 0L) {
        paste0("status=", probe$status_code)
      } else {
        character()
      },
      probe$output_tail
    )
    detail <- detail[nzchar(detail)]
    failures <- c(
      failures,
      sprintf(
        "[%s] %s",
        sub("_probe$", "", probe_name),
        paste(detail, collapse = " | ")
      )
    )
  }
}

if (!identical(partial_check_probe$status, "passed-no-errors")) {
  failures <- c(
    failures,
    sprintf(
      "[partial-check] %s",
      paste(partial_check_probe$error_lines, collapse = " | ")
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
  release_regression = list(
    status = if (release_probe$status == 0L) "passed" else "failed",
    output_tail = tail(release_probe$output, 12L)
  ),
  interactive_tests = list(
    transform_utils_validate = transform_probe,
    vce = vce_probe,
    visualization_export = visualization_probe
  ),
  static_hardening_partial_check = partial_check_probe,
  exact_status = if (length(failures) == 0L) "passed" else "failed",
  numeric_status = "not-applicable",
  blocker_boundary = if (length(failures) == 0L) {
    "interactive-test-bucket-cleared-static-doc-bucket-persists"
  } else {
    "interactive-test-bucket-still-open"
  },
  remaining_gap = if (length(failures) == 0L) {
    "full-check-time-tests-docs-still-need-fresh-completed-audit"
  } else {
    "interactive-targeted-bucket-not-fully-cleared"
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
  "# story-E8-04 Layer 5 interactive bucket probe",
  "",
  "## 时间",
  sprintf("- `%s`", payload$generated_at),
  "",
  "## 范围",
  "- 重跑 `e8_04_layer5_release_regression.R`，确认 targeted release probe 继续通过。",
  "- fresh rerun `transform-detrend|utils|validate`，验证 18:24 旧 tests bucket 中的三条交互式入口是否仍失败。",
  "- fresh rerun `vce-return-values|vce-integration` 与 `visualization-export`，确认旧的 stale release/VCE bucket 未回流。",
  "- 运行 `devtools::check(--no-manual, --no-tests, --no-examples, document = FALSE)`，只冻结 static hardening 的 WARNING/NOTE 桶，不替代 full check-time 结论。",
  "",
  "## 关键结果",
  sprintf(
    "- targeted release regression: `%s`。",
    payload$release_regression$status
  ),
  sprintf(
    "- `transform-detrend|utils|validate`: `%s`。",
    payload$interactive_tests$transform_utils_validate$status
  ),
  sprintf(
    "- `vce-return-values|vce-integration`: `%s`。",
    payload$interactive_tests$vce$status
  ),
  sprintf(
    "- `visualization-export`: `%s`。",
    payload$interactive_tests$visualization_export$status
  ),
  sprintf(
    "- static hardening partial check: `%s`。",
    payload$static_hardening_partial_check$status
  ),
  sprintf(
    "- partial check WARNING lines: `%s`。",
    if (length(partial_check_warning_lines) == 0L) {
      "none"
    } else {
      paste(partial_check_warning_lines, collapse = "`；`")
    }
  ),
  sprintf(
    "- partial check NOTE lines: `%s`。",
    if (length(partial_check_note_lines) == 0L) {
      "none"
    } else {
      paste(partial_check_note_lines, collapse = "`；`")
    }
  ),
  "",
  "## 裁决",
  if (length(failures) == 0L) {
    "- 交互式 Layer 5 test bucket 当前为绿色：`transform-detrend|utils|validate`、`vce-*` 与 `visualization-export` fresh rerun 均未复现 live failure。"
  } else {
    "- 交互式 Layer 5 test bucket 仍存在 failure，不能把旧 blocker 移出当前叙述。"
  },
  "- static hardening partial check 仍保留 `R code for possible problems` 与 `Rd/codoc/missing documentation/usage` WARNING/NOTE 桶。",
  "- 该 probe 不覆盖 fresh completed `devtools::check(--no-manual)` 的 check-time tests 结论；是否还能在 `R CMD check` 环境中复现 tests bucket，仍需单独完整审计。",
  "",
  "## 产出文件",
  sprintf("- JSON: `%s`", json_path),
  sprintf("- 说明: `%s`", md_path)
)

writeLines(md_lines, md_path)

if (length(failures) > 0L) {
  stop(
    paste(
      c("Layer 5 interactive bucket probe failed:", failures),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("Layer 5 interactive bucket probe passed.")
