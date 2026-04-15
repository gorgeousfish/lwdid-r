library(testthat)

resolve_smoking_controls_parity_path <- function(filename) {
  resolve_parity_fixture_path(filename)
}

resolve_smoking_raw_csv_path <- function() {
  resolve_parity_fixture_path("smoking.csv")
}

load_smoking_controls_fixture <- function(data, fixture_path, controls) {
  fixture <- utils::read.csv(fixture_path, stringsAsFactors = FALSE)
  merged <- merge(
    data[, setdiff(names(data), controls), drop = FALSE],
    fixture,
    by = "state",
    all.x = TRUE,
    sort = FALSE
  )

  merged[order(match(merged$state, unique(data$state)), merged$year), , drop = FALSE]
}

test_that("E8-03 smoking frozen-controls oracle freezes paper-backed simple-controls path", {
  skip_if_no_parity_fixtures()
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_smoking_controls_parity_path(
    "20260323-qa-parity-e8-03-smoking-controls-comparator.json"
  )
  fixture_path <- resolve_smoking_controls_parity_path(
    "e8_03_smoking_frozen_controls_fixture.csv"
  )
  raw_csv_path <- resolve_smoking_raw_csv_path()

  skip_if(!file.exists(oracle_path), paste("missing oracle:", oracle_path))
  skip_if(!file.exists(fixture_path), paste("missing fixture:", fixture_path))
  skip_if(!file.exists(raw_csv_path), paste("missing raw csv:", raw_csv_path))
  if (!file.exists(raw_csv_path)) {
    return(invisible(NULL))
  }

  oracle <- jsonlite::fromJSON(oracle_path, simplifyVector = FALSE)
  smoking <- utils::read.csv(raw_csv_path, stringsAsFactors = FALSE)
  controls <- c("lnincome", "beer", "age15to24", "lretprice")
  smoking_frozen <- load_smoking_controls_fixture(smoking, fixture_path, controls)

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = smoking_frozen,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      d = "d",
      post = "post",
      controls = controls,
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
  expect_true(any(grepl("simple controls", warnings_seen, fixed = TRUE)))
  .filter_glm_warnings <- function(w) {
    w[!grepl("^glm\\.fit:", w) & !grepl("^no non-missing arguments", w)]
  }
  expect_equal(
    sort(unique(.filter_glm_warnings(warnings_seen))),
    sort(unique(.filter_glm_warnings(
      unlist(oracle$r_result$warnings, use.names = FALSE)
    )))
  )

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
  expect_equal(result$estimator_comparison$ipw, ref$estimator$ipw, tolerance = 1e-6)
  expect_equal(result$estimator_comparison$range, ref$estimator$range, tolerance = 1e-6)
  expect_equal(
    result$estimator_comparison$rel_range,
    ref$estimator$rel_range,
    tolerance = 1e-6
  )
  expect_equal(
    result$estimator_comparison$baseline_att,
    ref$estimator$baseline_att,
    tolerance = 1e-6
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
  expect_gt(
    abs(result$estimator_comparison$ra - no_controls_ref$transformation$demean_att),
    1e-3
  )
})
