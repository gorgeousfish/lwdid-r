test_that("E8-08 Task 10 closure boundary wrapper freezes full-check pending state", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_parity_fixture_path(
    "20260327-story-worker-e8-08-task10-closure-boundary.json"
  )

  if (!file.exists(oracle_path)) {
    fail("Expected Task E8-08.10 closure boundary wrapper to exist.")
    return(invisible(NULL))
  }

  oracle <- jsonlite::read_json(oracle_path, simplifyVector = FALSE)

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.10")
  expect_identical(oracle$layer, "closure-boundary")
  expect_identical(oracle$exact_status, "passed")
  expect_true(isTRUE(oracle$targeted_release_contracts_cleared))
  expect_identical(oracle$unified_release_probe$task_status, "probe-green")
  expect_identical(oracle$interactive_bucket$exact_status, "passed")
  expect_identical(
    oracle$interactive_bucket$remaining_gap,
    "full-check-time-tests-docs-still-need-fresh-completed-audit"
  )
  expect_identical(
    oracle$closure_boundary,
    "full-check-closure-readiness-pending"
  )
  expect_identical(
    oracle$remaining_gap,
    "full-check-time-tests-docs-still-need-fresh-completed-audit"
  )
})

test_that("E8-08 Task 10 full-check audit wrapper freezes completed audit contract", {
  skip_if_not_installed("jsonlite")
  skip_if(
    identical(Sys.getenv("LWDID_TASK10_FULL_CHECK_SELF_AUDIT"), "1"),
    "Task E8-08.10 full-check wrapper is generated after the self-audit check completes."
  )

  oracle_path <- resolve_parity_fixture_path(
    "20260327-story-worker-e8-08-task10-full-check-audit.json"
  )

  if (!file.exists(oracle_path)) {
    fail("Expected Task E8-08.10 full-check audit wrapper to exist.")
    return(invisible(NULL))
  }

  oracle <- jsonlite::read_json(oracle_path, simplifyVector = FALSE)

  expect_identical(oracle$story, "story-E8-08")
  expect_identical(oracle$task, "E8-08.10")
  expect_identical(oracle$layer, "full-check-audit")

  if (identical(Sys.getenv("LWDID_TASK10_FULL_CHECK_SELF_AUDIT", ""), "1")) {
    expect_true(length(oracle$remaining_gap) == 1L)
    return(invisible(NULL))
  }

  expect_identical(oracle$full_check$status, "completed")
  expect_true(nzchar(oracle$full_check$result))
  expect_true(grepl("^[0-9]+ errors / [0-9]+ warnings / [0-9]+ notes$", oracle$full_check$result))
  expect_false(is.null(oracle$full_check$command_status))
  expect_true(nzchar(oracle$full_check$status_line))
  expect_true(file.exists(oracle$full_check$check_log))
  expect_true(file.exists(oracle$full_check$testthat_rout))
  expect_false(is.null(oracle$full_check$tests$summary_line))
  expect_true(length(oracle$remaining_gap) == 1L)
  expect_true(
    oracle$story_status %in% c(
      "closure-ready-confirmed",
      "closure-ready-blocked"
    )
  )
})

test_that("E8-08 Task 10 unified suite bootstrap avoids undeclared pkgload dependency", {
  bootstrap_path <- resolve_package_source_file(
    "tests",
    "testthat",
    "test-sensitivity-diagnostics.R"
  )

  bootstrap_lines <- readLines(bootstrap_path, warn = FALSE)

  expect_false(any(grepl("pkgload", bootstrap_lines, fixed = TRUE)))
  expect_true(any(grepl('requireNamespace\\("devtools"', bootstrap_lines)))
  expect_true(any(grepl('getExportedValue\\("devtools", "load_all"\\)', bootstrap_lines)))
})
