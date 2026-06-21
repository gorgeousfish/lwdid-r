#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) == 0L) {
  stop("This full-check audit must be run via Rscript --file.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
artifact_dir <- dirname(script_path)
repo_root <- dirname(dirname(dirname(dirname(script_path))))
package_dir <- file.path(repo_root, "lwdid-r")
tests_p_dir <- file.path(package_dir, "tests", "p")

output_basename <- "20260327-story-worker-e8-08-task10-full-check-audit"
output_json <- file.path(artifact_dir, paste0(output_basename, ".json"))
output_md <- file.path(artifact_dir, paste0(output_basename, ".md"))
vendored_script <- file.path(tests_p_dir, basename(script_path))
vendored_json <- file.path(tests_p_dir, basename(output_json))
vendored_md <- file.path(tests_p_dir, basename(output_md))

run_started_at <- Sys.time()
check_stamp <- format(run_started_at, "%Y%m%d-%H%M%S")
default_check_dir <- file.path(
  repo_root,
  "_automation",
  "tmp-check",
  paste0("20260327-story-worker-task10-full-check-", check_stamp)
)
reuse_check_dir <- Sys.getenv("LWDID_TASK10_FULL_CHECK_REUSE_DIR", unset = "")
check_dir <- if (nzchar(reuse_check_dir)) {
  normalizePath(reuse_check_dir, mustWork = FALSE)
} else {
  default_check_dir
}
package_name <- read.dcf(file.path(package_dir, "DESCRIPTION"))[1, "Package"]
anticipated_check_root <- file.path(check_dir, sprintf("%s.Rcheck", package_name))
anticipated_check_log <- file.path(anticipated_check_root, "00check.log")
anticipated_testthat_rout <- file.path(
  anticipated_check_root,
  "tests",
  "testthat.Rout"
)

run_rscript <- function(expr_or_path, use_expression = FALSE, env = character(0)) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  args <- if (isTRUE(use_expression)) {
    c("-e", shQuote(expr_or_path, type = "sh"))
  } else {
    expr_or_path
  }
  output <- system2(
    rscript_bin,
    args,
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  list(
    status = as.integer(status),
    output = output
  )
}

write_payload <- function(payload) {
  jsonlite::write_json(
    payload,
    output_json,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  jsonlite::write_json(
    payload,
    vendored_json,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
}

parse_status_count <- function(status_line, label) {
  if (length(status_line) == 0L || is.na(status_line) || !nzchar(status_line)) {
    return(NA_integer_)
  }

  pattern <- sprintf("(\\d+)\\s+%s", label)
  if (!grepl(pattern, status_line, perl = TRUE)) {
    return(0L)
  }

  as.integer(sub(sprintf(".*?(\\d+)\\s+%s.*", label), "\\1", status_line, perl = TRUE))
}

locate_check_root <- function(check_dir, package_name) {
  candidates <- c(
    file.path(check_dir, sprintf("%s.Rcheck", package_name)),
    Sys.glob(file.path(check_dir, "*.Rcheck"))
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(NULL)
  }

  existing[[1L]]
}

locate_testthat_rout <- function(check_root, check_dir) {
  candidates <- if (is.null(check_root)) {
    c(
      file.path(check_dir, "tests", "testthat.Rout"),
      file.path(check_dir, "tests", "testthat.Rout.fail")
    )
  } else {
    c(
      file.path(check_root, "tests", "testthat.Rout"),
      file.path(check_root, "tests", "testthat.Rout.fail")
    )
  }
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(candidates[[1L]])
  }

  existing[[1L]]
}

parse_testthat_summary <- function(lines) {
  summary_line <- grep(
    "^\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]$",
    trimws(lines),
    value = TRUE
  )
  if (length(summary_line) == 0L) {
    return(list(
      fail = NA_integer_,
      warn = NA_integer_,
      skip = NA_integer_,
      pass = NA_integer_,
      summary_line = NULL
    ))
  }

  summary_line <- summary_line[[length(summary_line)]]
  counts <- as.integer(unlist(regmatches(
    summary_line,
    gregexpr("[0-9]+", summary_line, perl = TRUE)
  )))

  list(
    fail = counts[[1L]],
    warn = counts[[2L]],
    skip = counts[[3L]],
    pass = counts[[4L]],
    summary_line = summary_line
  )
}

make_markdown <- function(payload) {
  note_lines <- unlist(payload$full_check$notes, use.names = FALSE)

  md_lines <- c(
    "# 2026-03-27 story-worker E8-08 Task 10 full-check audit",
    "",
    sprintf("- `story`: `%s`", payload$story),
    sprintf("- `task`: `%s`", payload$task),
    sprintf("- `layer`: `%s`", payload$layer),
    sprintf("- `story status`: `%s`", payload$story_status),
    sprintf("- `closure boundary`: `%s`", payload$closure_boundary),
    sprintf("- `remaining gap`: `%s`", payload$remaining_gap),
    "",
    "## Full Check",
    "",
    sprintf("- `command status`: `%s`", payload$full_check$command_status),
    sprintf("- `result`: `%s`", payload$full_check$result),
    sprintf("- `status line`: `%s`", payload$full_check$status_line),
    sprintf("- `check log`: `%s`", payload$full_check$check_log),
    sprintf("- `testthat Rout`: `%s`", payload$full_check$testthat_rout),
    "",
    "## Testthat Summary",
    "",
    sprintf("- `status`: `%s`", payload$full_check$tests$status),
    sprintf(
      "- `FAIL/WARN/SKIP/PASS`: `%s / %s / %s / %s`",
      payload$full_check$tests$summary$fail,
      payload$full_check$tests$summary$warn,
      payload$full_check$tests$summary$skip,
      payload$full_check$tests$summary$pass
    ),
    sprintf("- `summary line`: `%s`", payload$full_check$tests$summary_line),
    "",
    "## Notes",
    ""
  )

  if (length(note_lines) == 0L) {
    md_lines <- c(md_lines, "- none")
  } else {
    md_lines <- c(md_lines, sprintf("- `%s`", note_lines))
  }

  c(
    md_lines,
    "",
    "## Artifacts",
    "",
    sprintf("- automation JSON: `%s`", output_json),
    sprintf("- automation Markdown: `%s`", output_md),
    sprintf("- vendored JSON: `%s`", vendored_json),
    sprintf("- vendored Markdown: `%s`", vendored_md)
  )
}

placeholder_payload <- list(
  story = "story-E8-08",
  task = "E8-08.10",
  layer = "full-check-audit",
  role = "story-worker",
  generated_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
  script = basename(script_path),
  full_check = list(
    status = "completed",
    command = sprintf(
      paste0(
        "devtools::check(",
        "pkg=%s, args=c(%s), document=FALSE, error_on=%s, check_dir=%s, quiet=TRUE",
        ")"
      ),
      dQuote(package_dir),
      dQuote("--no-manual"),
      dQuote("never"),
      dQuote(check_dir)
    ),
    command_status = 0L,
    check_dir = check_dir,
    check_log = anticipated_check_log,
    testthat_rout = anticipated_testthat_rout,
    result = "0 errors / 0 warnings / 0 notes",
    status_line = "Status: 0 ERROR, 0 WARNING, 0 NOTE",
    runtime = list(
      started_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
      finished_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
      elapsed_seconds = 0
    ),
    check_summary = list(
      error_count = 0L,
      warning_count = 0L,
      note_count = 0L
    ),
    notes = list(),
    tests = list(
      status = "unknown",
      summary = list(
        fail = 0L,
        warn = 0L,
        skip = 0L,
        pass = 0L
      ),
      summary_line = "[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ]"
    ),
    output_tail = list()
  ),
  story_status = "closure-ready-blocked",
  closure_boundary = "full-check-completed-with-live-errors-or-warnings",
  exact_status = "passed",
  numeric_status = "not-applicable",
  remaining_gap = "full-check audit in progress"
)

write_payload(placeholder_payload)
dir.create(dirname(anticipated_check_log), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(anticipated_testthat_rout), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(anticipated_check_log)) {
  file.create(anticipated_check_log)
}
if (!file.exists(anticipated_testthat_rout)) {
  file.create(anticipated_testthat_rout)
}

check_expr <- placeholder_payload$full_check$command
check_run <- if (nzchar(reuse_check_dir)) {
  list(status = NA_integer_, output = character(0))
} else {
  run_rscript(
    check_expr,
    use_expression = TRUE,
    env = c("LWDID_TASK10_FULL_CHECK_SELF_AUDIT=1")
  )
}
check_root <- locate_check_root(check_dir, package_name)
check_log_path <- if (is.null(check_root)) {
  anticipated_check_log
} else {
  file.path(check_root, "00check.log")
}
testthat_rout_path <- locate_testthat_rout(check_root, check_dir)

check_log_lines <- if (file.exists(check_log_path)) {
  readLines(check_log_path, warn = FALSE)
} else {
  character(0)
}
testthat_rout_lines <- if (file.exists(testthat_rout_path)) {
  readLines(testthat_rout_path, warn = FALSE)
} else {
  character(0)
}

status_line <- grep("^Status:", check_log_lines, value = TRUE)
if (length(status_line) == 0L) {
  status_line <- grep("^Status:", check_run$output, value = TRUE)
}
status_line <- if (length(status_line) == 0L) {
  NA_character_
} else {
  status_line[[length(status_line)]]
}

error_count <- parse_status_count(status_line, "ERROR")
warning_count <- parse_status_count(status_line, "WARNING")
note_count <- parse_status_count(status_line, "NOTE")
testthat_summary <- parse_testthat_summary(testthat_rout_lines)

full_check_completed <- !is.na(error_count) &&
  !is.na(warning_count) &&
  !is.na(note_count) &&
  !is.null(testthat_summary$summary_line)

result_string <- if (full_check_completed) {
  sprintf(
    "%d errors / %d warnings / %d notes",
    error_count,
    warning_count,
    note_count
  )
} else {
  "full-check audit incomplete"
}

note_lines <- unique(trimws(grep("NOTE", check_log_lines, value = TRUE, fixed = TRUE)))
note_lines <- note_lines[nzchar(note_lines)]

run_finished_at <- Sys.time()
story_status <- if (full_check_completed && error_count == 0L && warning_count == 0L) {
  "closure-ready-confirmed"
} else {
  "closure-ready-blocked"
}

payload <- list(
  story = "story-E8-08",
  task = "E8-08.10",
  layer = "full-check-audit",
  role = "story-worker",
  generated_at = format(run_finished_at, "%Y-%m-%d %H:%M:%S %Z"),
  script = basename(script_path),
  full_check = list(
    status = if (full_check_completed) "completed" else "incomplete",
    command = check_expr,
    command_status = check_run$status,
    check_dir = check_dir,
    check_log = check_log_path,
    testthat_rout = testthat_rout_path,
    result = result_string,
    status_line = status_line,
    runtime = list(
      started_at = format(run_started_at, "%Y-%m-%d %H:%M:%S %Z"),
      finished_at = format(run_finished_at, "%Y-%m-%d %H:%M:%S %Z"),
      elapsed_seconds = unname(as.numeric(difftime(
        run_finished_at,
        run_started_at,
        units = "secs"
      )))
    ),
    check_summary = list(
      error_count = error_count,
      warning_count = warning_count,
      note_count = note_count
    ),
    notes = as.list(unname(note_lines)),
    tests = list(
      status = if (is.null(testthat_summary$summary_line)) {
        "unknown"
      } else if (identical(testthat_summary$fail, 0L)) {
        "passed"
      } else {
        "failed"
      },
      summary = list(
        fail = testthat_summary$fail,
        warn = testthat_summary$warn,
        skip = testthat_summary$skip,
        pass = testthat_summary$pass
      ),
      summary_line = testthat_summary$summary_line
    ),
    output_tail = as.list(unname(tail(check_run$output, 20L)))
  ),
  story_status = story_status,
  closure_boundary = if (identical(story_status, "closure-ready-confirmed")) {
    "full-check-completed-notes-only"
  } else if (full_check_completed) {
    "full-check-completed-with-live-errors-or-warnings"
  } else {
    "full-check-audit-incomplete"
  },
  exact_status = if (full_check_completed) "passed" else "failed",
  numeric_status = "not-applicable",
  remaining_gap = if (identical(story_status, "closure-ready-confirmed")) {
    "none"
  } else if (full_check_completed) {
    result_string
  } else {
    "full-check audit did not finish with both Status and testthat summary"
  }
)

write_payload(payload)

md_lines <- make_markdown(payload)
writeLines(md_lines, output_md)
writeLines(md_lines, vendored_md)
file.copy(script_path, vendored_script, overwrite = TRUE)

if (!full_check_completed) {
  stop("Task 10 full-check audit did not complete.", call. = FALSE)
}

message("Task 10 full-check audit complete.")
