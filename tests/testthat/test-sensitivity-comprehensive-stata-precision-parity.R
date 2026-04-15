library(testthat)

resolve_smoking_stata_precision_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

test_that("Stata precision oracle resolver prefers vendored package fixtures", {
  raw_path <- tryCatch(
    testthat::test_path(
      "..", "p", "20260324-qa-parity-e8-03-smoking-stata-precision.json"
    ),
    error = function(...) NULL
  )
  skip_if(
    is.null(raw_path) || !file.exists(raw_path),
    "Vendored parity fixture not available (excluded from build)"
  )
  fixture_path <- normalizePath(raw_path, mustWork = TRUE)

  expect_true(file.exists(fixture_path))
  expect_identical(
    resolve_smoking_stata_precision_path(
      "20260324-qa-parity-e8-03-smoking-stata-precision.json"
    ),
    fixture_path
  )
})

test_that("E8-03 smoking comprehensive transformation stays on Stata precision oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_smoking_stata_precision_path(
    "20260324-qa-parity-e8-03-smoking-stata-precision.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing Stata precision oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  data("smoking", package = "lwdid", envir = environment())

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_identical(oracle$comparison$status, "passed")
  expect_identical(oracle$comparison$numeric_status, "passed")
  expect_equal(
    sort(unique(warnings_seen)),
    sort(unique(unlist(oracle$r_current$warnings, use.names = FALSE)))
  )

  transform_ref <- oracle$stata_reference$transformation
  tolerances <- oracle$tolerances

  expect_equal(
    result$transformation_comparison$demean_att,
    transform_ref$demean_att,
    tolerance = tolerances$att
  )
  expect_equal(
    result$transformation_comparison$demean_se,
    transform_ref$demean_se,
    tolerance = tolerances$se
  )
  expect_equal(
    result$transformation_comparison$detrend_att,
    transform_ref$detrend_att,
    tolerance = tolerances$att
  )
  expect_equal(
    result$transformation_comparison$detrend_se,
    transform_ref$detrend_se,
    tolerance = tolerances$se
  )
  expect_equal(
    result$transformation_comparison$difference,
    transform_ref$difference,
    tolerance = tolerances$att
  )
  expect_equal(
    result$transformation_comparison$rel_diff,
    transform_ref$rel_diff,
    tolerance = tolerances$att
  )

  current_numeric <- unname(unlist(result$transformation_comparison, use.names = FALSE))
  expect_true(all(is.finite(as.numeric(current_numeric))))
  expect_null(result$estimator_comparison)
})

test_that("E8-03 smoking comprehensive pre-period and no-anticipation stay on Stata precision oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_smoking_stata_precision_path(
    "20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing Stata sensitivity-branches oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  data("smoking", package = "lwdid", envir = environment())

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_identical(oracle$comparison$numeric_status, "passed")
  expect_equal(
    sort(unique(warnings_seen)),
    sort(unique(unlist(oracle$r_current$warnings, use.names = FALSE)))
  )
  expect_true(
    any(vapply(
      oracle$comparison$blocked_specs,
      function(item) {
        identical(item$branch, "pre_period") &&
          identical(as.integer(item$index), 1L) &&
          identical(as.integer(item$stata_rc), 2001L)
      },
      logical(1)
    ))
  )

  tolerances <- oracle$tolerances
  pre_ref <- Filter(
    function(item) identical(as.integer(item$stata_rc), 0L),
    oracle$stata_reference$pre_period$specs
  )
  pre_current <- result$pre_period_result$specifications

  pre_current_by_n <- stats::setNames(
    pre_current,
    vapply(pre_current, function(item) item$n_pre_periods, integer(1))
  )
  expect_length(pre_ref, length(pre_current) - 1L)
  for (ref in pre_ref) {
    current <- pre_current_by_n[[as.character(ref$n_pre_periods)]]
    expect_true(isTRUE(current$converged))
    expect_equal(current$n_pre_periods, ref$n_pre_periods)
    expect_equal(current$att, ref$att, tolerance = tolerances$att)
    expect_equal(current$se, ref$se, tolerance = tolerances$se)
  }

  noa_ref <- Filter(
    function(item) identical(as.integer(item$stata_rc), 0L),
    oracle$stata_reference$no_anticipation$estimates
  )
  noa_current <- result$no_anticipation_result$estimates

  expect_length(noa_current, length(noa_ref))
  for (idx in seq_along(noa_ref)) {
    expect_true(isTRUE(noa_current[[idx]]$converged))
    expect_equal(
      noa_current[[idx]]$excluded_periods,
      noa_ref[[idx]]$excluded_periods
    )
    expect_equal(
      noa_current[[idx]]$n_pre_periods_used,
      noa_ref[[idx]]$n_pre_periods_used
    )
    expect_equal(noa_current[[idx]]$att, noa_ref[[idx]]$att, tolerance = tolerances$att)
    expect_equal(noa_current[[idx]]$se, noa_ref[[idx]]$se, tolerance = tolerances$se)
  }

  expect_identical(
    result$no_anticipation_result$anticipation_detected,
    oracle$r_current$no_anticipation$anticipation_detected
  )
  expect_equal(
    result$no_anticipation_result$recommended_exclusion,
    oracle$r_current$no_anticipation$recommended_exclusion
  )
})

test_that("E8-03 smoking n_pre=1 Stata boundary is explicitly waived in the oracle contract", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_smoking_stata_precision_path(
    "20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing Stata sensitivity-branches oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)

  decisions <- oracle$comparison$boundary_decisions
  expect_length(decisions, 1L)

  decision <- decisions[[1L]]
  expect_identical(decision$branch, "pre_period")
  expect_identical(as.integer(decision$index), 1L)
  expect_identical(as.integer(decision$stata_rc), 2001L)
  expect_identical(decision$decision, "waive")
  expect_identical(decision$reason, "stata-implementation-boundary")
  expect_true(isTRUE(decision$allows_story_task_closure))
  expect_true("papers-and-math" %in% unlist(decision$truth_basis, use.names = FALSE))
  expect_true("r-python-agreement" %in% unlist(decision$truth_basis, use.names = FALSE))
  expect_true(any(grepl(
    "20260324-theory-parity-e8-03-pre-period-npre1-manual-oracle\\.json$",
    unlist(decision$evidence, use.names = FALSE)
  )))
})

test_that("E8-03 smoking n_pre=1 manual oracle is executable in the QA contract", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_smoking_stata_precision_path(
    "20260324-qa-parity-e8-03-smoking-stata-sensitivity-branches.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing Stata sensitivity-branches oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  manual_oracle <- oracle$manual_oracle$pre_period_n_pre1

  expect_true(
    !is.null(manual_oracle),
    info = "QA oracle should embed the n_pre=1 manual oracle contract"
  )
  if (is.null(manual_oracle)) {
    return(invisible(NULL))
  }

  data("smoking", package = "lwdid", envir = environment())
  result <- suppressWarnings(
    lwdid_sensitivity(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      type = "all",
      verbose = FALSE
    )
  )

  current <- Filter(
    function(item) identical(as.integer(item$n_pre_periods), 1L),
    result$pre_period_result$specifications
  )

  expect_length(current, 1L)
  current <- current[[1L]]

  expect_equal(current$att, manual_oracle$att, tolerance = oracle$tolerances$att)
  expect_equal(current$se, manual_oracle$se, tolerance = oracle$tolerances$se)
  expect_equal(
    current$ci_lower,
    manual_oracle$ci_lower,
    tolerance = oracle$tolerances$se
  )
  expect_equal(
    current$ci_upper,
    manual_oracle$ci_upper,
    tolerance = oracle$tolerances$se
  )
  expect_equal(
    current$pvalue,
    manual_oracle$pvalue,
    tolerance = oracle$tolerances$att
  )
})
