library(testthat)

resolve_castle_parity_oracle_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

test_that("E8-03 castle parity oracle matches current staggered comprehensive sensitivity output", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_castle_parity_oracle_path(
    "20260323-qa-parity-e8-03-castle-comparator.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing parity oracle:", oracle_path)
  )
  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  data("castle", package = "lwdid", envir = environment())

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_true(isTRUE(oracle$comparison$numeric_status == "passed"))
  expect_true(isTRUE(oracle$comparison$exact_status == "passed"))
  expect_equal(
    sort(unique(warnings_seen)),
    sort(unique(unlist(oracle$r_result$warnings, use.names = FALSE)))
  )

  py_ref <- oracle$python_reference$result
  transform_ref <- py_ref$transformation

  expect_identical(
    isTRUE(result$pre_period_result$is_robust),
    isTRUE(py_ref$pre_period$is_robust)
  )
  expect_identical(
    as.character(result$pre_period_result$robustness_level),
    py_ref$pre_period$robustness_level
  )
  expect_identical(
    isTRUE(result$no_anticipation_result$anticipation_detected),
    isTRUE(py_ref$anticipation$anticipation_detected)
  )
  expect_identical(
    as.integer(result$no_anticipation_result$recommended_exclusion),
    as.integer(py_ref$anticipation$recommended_exclusion)
  )
  expect_true(is.null(result$estimator_comparison))
  expect_true(is.null(py_ref$estimator))

  expect_equal(
    result$pre_period_result$sensitivity_ratio,
    py_ref$pre_period$sensitivity_ratio,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$demean_att,
    transform_ref$demean_att,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$demean_se,
    transform_ref$demean_se,
    tolerance = 1e-4
  )
  expect_equal(
    result$transformation_comparison$detrend_att,
    transform_ref$detrend_att,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$detrend_se,
    transform_ref$detrend_se,
    tolerance = 1e-4
  )
  expect_equal(
    result$transformation_comparison$difference,
    transform_ref$difference,
    tolerance = 1e-6
  )
})
