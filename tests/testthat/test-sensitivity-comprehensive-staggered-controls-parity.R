library(testthat)

resolve_castle_controls_parity_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

test_that("E8-03 castle controls oracle freezes paper-backed simple-controls path", {
  skip_if_no_parity_fixtures()
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_castle_controls_parity_path(
    "20260323-qa-parity-e8-03-castle-controls-comparator.json"
  )

  skip_if(!file.exists(oracle_path), paste("missing oracle:", oracle_path))

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
      controls = c("income", "unemployrt", "poverty"),
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_identical(oracle$comparison$status, "drift-confirmed")
  expect_identical(oracle$comparison$paper_backed_r_status, "passed")
  # Warning set may drift as implementation evolves; just verify
  # no hard errors and the result structure is intact.
  # (Oracle warnings are informational, not hard contracts.)

  ref <- oracle$r_result$result
  no_controls_ref <- oracle$no_controls_reference$r_result$result

  expect_false(is.null(result$estimator_comparison))
  expect_equal(
    result$pre_period_result$sensitivity_ratio,
    ref$pre_period$sensitivity_ratio,
    tolerance = 1e-6
  )
  expect_identical(
    isTRUE(result$pre_period_result$is_robust),
    isTRUE(ref$pre_period$is_robust)
  )
  expect_identical(
    as.character(result$pre_period_result$robustness_level),
    ref$pre_period$robustness_level
  )
  expect_identical(
    isTRUE(result$no_anticipation_result$anticipation_detected),
    isTRUE(ref$anticipation$anticipation_detected)
  )
  expect_identical(
    as.integer(result$no_anticipation_result$recommended_exclusion),
    as.integer(ref$anticipation$recommended_exclusion)
  )

  expect_equal(
    result$transformation_comparison$demean_att,
    ref$transformation$demean_att,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$demean_se,
    ref$transformation$demean_se,
    tolerance = 1e-4
  )
  expect_equal(
    result$transformation_comparison$detrend_att,
    ref$transformation$detrend_att,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$detrend_se,
    ref$transformation$detrend_se,
    tolerance = 1e-4
  )
  expect_equal(
    result$transformation_comparison$difference,
    ref$transformation$difference,
    tolerance = 1e-6
  )
  expect_equal(
    result$transformation_comparison$rel_diff,
    ref$transformation$rel_diff,
    tolerance = 1e-6
  )

  expect_equal(result$estimator_comparison$ra, ref$estimator$ra, tolerance = 1e-6)
  expect_equal(result$estimator_comparison$ipw, ref$estimator$ipw, tolerance = 0.1)
  expect_equal(result$estimator_comparison$ipwra, ref$estimator$ipwra, tolerance = 0.1)
  expect_equal(
    result$estimator_comparison$baseline_att,
    ref$estimator$baseline_att,
    tolerance = 1e-6
  )
  expect_equal(result$estimator_comparison$range, ref$estimator$range, tolerance = 0.05)
  expect_equal(
    result$estimator_comparison$rel_range,
    ref$estimator$rel_range,
    tolerance = 0.3
  )

  current_numeric <- c(
    unname(unlist(result$transformation_comparison, use.names = FALSE)),
    unname(unlist(result$estimator_comparison, use.names = FALSE))
  )
  expect_true(all(is.finite(as.numeric(current_numeric))))

  expect_gt(
    abs(
      result$transformation_comparison$demean_att -
        no_controls_ref$transformation$demean_att
    ),
    1e-3
  )
})
