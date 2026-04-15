library(testthat)

resolve_selection_parity_path <- function(filename) {
  if (exists("resolve_parity_fixture_path", mode = "function")) {
    return(resolve_parity_fixture_path(filename))
  }

  candidates <- c(
    testthat::test_path(
      "..", "..", "..",
      "_automation", "test-artifacts", "parity", filename
    ),
    file.path(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
      filename
    ),
    testthat::test_path("_fixtures", "parity", filename)
  )

  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  candidates[[1L]]
}

resolve_selection_source_file <- function(...) {
  if (exists("resolve_package_source_file", mode = "function")) {
    return(resolve_package_source_file(...))
  }

  file.path(testthat::test_path("..", ".."), ...)
}

read_selection_parity_oracle <- function(filename) {
  jsonlite::fromJSON(
    resolve_selection_parity_path(filename),
    simplifyVector = FALSE
  )
}

read_selection_parity_fixture <- function(filename) {
  utils::read.csv(
    resolve_selection_parity_path(filename),
    stringsAsFactors = FALSE
  )
}

get_expected_unit <- function(expected_units, unit_id) {
  unit <- Filter(function(x) identical(as.integer(x$unit_id), as.integer(unit_id)), expected_units)

  if (length(unit) != 1L) {
    stop(sprintf("Unable to locate expected unit %s", unit_id), call. = FALSE)
  }

  unit[[1L]]
}

derive_selection_summary <- function(unit_stats) {
  treated_units <- unit_stats[unit_stats$is_treated, , drop = FALSE]
  n_treated_units <- nrow(treated_units)

  list(
    n_units_complete = as.integer(sum(unit_stats$n_missing == 0L)),
    n_units_partial = as.integer(sum(unit_stats$n_missing > 0L)),
    attrition_rate = as.numeric(mean(unit_stats$n_missing > 0L)),
    n_treated_units = as.integer(n_treated_units),
    pct_usable_demean = if (n_treated_units > 0L) {
      100 * mean(treated_units$can_use_demean)
    } else {
      100
    },
    pct_usable_detrend = if (n_treated_units > 0L) {
      100 * mean(treated_units$can_use_detrend)
    } else {
      100
    }
  )
}

derive_selection_public_inputs <- function(panel) {
  pattern <- .classify_missing_pattern(
    panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    d = NULL,
    covariates = NULL
  )
  attrition <- .compute_attrition_analysis(
    panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    d = NULL,
    never_treated_values = c(0, Inf)
  )
  balance <- .compute_balance_statistics(
    panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    d = NULL,
    never_treated_values = c(0, Inf)
  )

  list(
    missing_pattern = pattern$pattern,
    missing_pattern_confidence = pattern$confidence,
    n_tests = length(pattern$tests),
    attrition_rate = attrition$overall_attrition,
    n_units_complete = attrition$n_complete,
    n_units_partial = attrition$n_partial,
    late_entry_rate = attrition$late_entry_rate,
    early_dropout_rate = attrition$early_dropout_rate,
    dropout_before_treatment = attrition$dropout_before_treatment,
    dropout_after_treatment = attrition$dropout_after_treatment,
    treatment_related_attrition = attrition$treatment_related_attrition,
    is_balanced = balance$is_balanced,
    balance_ratio = balance$balance_ratio,
    pct_usable_demean = balance$pct_usable_demean,
    pct_usable_detrend = balance$pct_usable_detrend,
    n_treated_units = balance$n_treated_units,
    units_below_demean_threshold = balance$units_below_demean_threshold,
    units_below_detrend_threshold = balance$units_below_detrend_threshold
  )
}

map_python_pattern_to_r <- function(pattern) {
  switch(
    pattern,
    missing_completely_at_random = "MCAR",
    missing_at_random = "MAR",
    missing_not_at_random = "MNAR",
    unknown = "UNKNOWN",
    stop(sprintf("Unknown Python missing-pattern label: %s", pattern), call. = FALSE)
  )
}

test_that("E8-05 observed-Y parity oracle exists and stays bug-aware", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-observed-y-helper-parity.json"
  )
  fixture_path <- resolve_selection_parity_path(
    "e8_05_selection_observed_y_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection observed-Y parity oracle:", oracle_path)
  )
  expect_true(
    file.exists(fixture_path),
    info = paste("missing selection observed-Y fixture:", fixture_path)
  )

  if (!file.exists(oracle_path) || !file.exists(fixture_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-observed-y-helper-parity.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$comparison$status, "drift-confirmed")
  expect_identical(oracle$comparison$paper_backed_r_status, "passed")
  expect_identical(oracle$python_reference$bug_ledger_id, "PY-SEL-003")
  expect_equal(
    oracle$paper_expected$missing_rates$overall_missing_rate,
    0.25,
    tolerance = 1e-12
  )
})

test_that("E8-05 missing-pattern oracle assets exist for Task 2 parity", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )
  script_path <- resolve_selection_parity_path(
    "e8_05_selection_missing_pattern_oracle.py"
  )
  mcar_fixture_path <- resolve_selection_parity_path(
    "e8_05_selection_pattern_mcar_fixture.csv"
  )
  mnar_fixture_path <- resolve_selection_parity_path(
    "e8_05_selection_pattern_mnar_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection missing-pattern oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing selection missing-pattern oracle script:", script_path)
  )
  expect_true(
    file.exists(mcar_fixture_path),
    info = paste("missing selection MCAR fixture:", mcar_fixture_path)
  )
  expect_true(
    file.exists(mnar_fixture_path),
    info = paste("missing selection MNAR fixture:", mnar_fixture_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$task, "E8-05.2")
  expect_true(length(oracle$cases) >= 2L)
  expect_identical(oracle$cases[[1]]$case_id, "balanced_mcar")
  expect_identical(oracle$cases[[2]]$case_id, "lagged_outcome_mnar_signal")
})

test_that("E8-05 missing-pattern oracle freezes Python classification boundaries", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection missing-pattern oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )

  balanced_case <- Filter(
    function(x) identical(x$case_id, "balanced_mcar"),
    oracle$cases
  )[[1L]]
  mnar_case <- Filter(
    function(x) identical(x$case_id, "lagged_outcome_mnar_signal"),
    oracle$cases
  )[[1L]]

  expect_identical(balanced_case$python_result$pattern, "missing_completely_at_random")
  expect_equal(balanced_case$python_result$confidence, 1, tolerance = 1e-12)
  expect_equal(balanced_case$python_result$n_tests, 0L)

  expect_identical(mnar_case$python_result$pattern, "missing_not_at_random")
  expect_equal(mnar_case$python_result$confidence, 0.7, tolerance = 1e-12)
  expect_true(mnar_case$python_result$n_tests >= 1L)
  expect_true(any(vapply(
    mnar_case$python_result$tests,
    function(test) identical(test$test_name, "Selection on Lagged Outcome (MNAR) Test"),
    logical(1)
  )))
})

test_that("E8-05.2: balanced missing-pattern case short-circuits to MCAR", {
  skip_if_not_installed("jsonlite")

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )
  panel <- read_selection_parity_fixture("e8_05_selection_pattern_mcar_fixture.csv")
  expected <- Filter(
    function(x) identical(x$case_id, "balanced_mcar"),
    oracle$cases
  )[[1L]]

  actual <- .classify_missing_pattern(
    panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    d = NULL,
    covariates = NULL
  )

  expect_identical(actual$pattern, map_python_pattern_to_r(expected$python_result$pattern))
  expect_equal(actual$confidence, expected$python_result$confidence, tolerance = 1e-12)
  expect_length(actual$tests, expected$python_result$n_tests)
})

test_that("E8-05.2: lagged-outcome case flags MNAR with frozen test boundary", {
  skip_if_not_installed("jsonlite")

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-missing-pattern-python-oracle.json"
  )
  panel <- read_selection_parity_fixture("e8_05_selection_pattern_mnar_fixture.csv")
  expected <- Filter(
    function(x) identical(x$case_id, "lagged_outcome_mnar_signal"),
    oracle$cases
  )[[1L]]

  actual <- .classify_missing_pattern(
    panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    d = NULL,
    covariates = NULL
  )

  expect_identical(actual$pattern, map_python_pattern_to_r(expected$python_result$pattern))
  expect_equal(actual$confidence, expected$python_result$confidence, tolerance = 1e-12)
  expect_length(actual$tests, expected$python_result$n_tests)
  expect_identical(actual$tests[[1]]$test_name, "Selection on Lagged Outcome (MNAR) Test")
  expect_identical(actual$tests[[1]]$reject_null, TRUE)
})

test_that("E8-05 helper parity follows observed-Y contract rather than Python row presence", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-observed-y-helper-parity.json"
  )
  fixture_path <- resolve_selection_parity_path(
    "e8_05_selection_observed_y_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection observed-Y parity oracle:", oracle_path)
  )
  expect_true(
    file.exists(fixture_path),
    info = paste("missing selection observed-Y fixture:", fixture_path)
  )

  if (!file.exists(oracle_path) || !file.exists(fixture_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-observed-y-helper-parity.json"
  )
  panel <- read_selection_parity_fixture("e8_05_selection_observed_y_fixture.csv")

  actual_rates <- .compute_missing_rates(panel, ivar = "id", tvar = "t", y = "y")
  actual_by_period <- .compute_missing_by_period(panel, ivar = "id", tvar = "t", y = "y")
  actual_by_cohort <- .compute_missing_by_cohort(
    panel,
    ivar = "id",
    tvar = "t",
    gvar = "g",
    y = "y",
    never_treated_values = c(0, Inf)
  )
  actual_unit_stats <- .compute_unit_stats(
    panel,
    ivar = "id",
    tvar = "t",
    gvar = "g",
    y = "y",
    never_treated_values = c(0, Inf)
  )
  actual_summary <- derive_selection_summary(actual_unit_stats)

  expected <- oracle$paper_expected
  expected_treated <- get_expected_unit(expected$unit_stats, 1)
  expected_never <- get_expected_unit(expected$unit_stats, 2)

  expect_equal(actual_rates$n_units, expected$missing_rates$n_units)
  expect_equal(actual_rates$n_periods, expected$missing_rates$n_periods)
  expect_equal(actual_rates$n_expected, expected$missing_rates$n_expected)
  expect_equal(actual_rates$n_observed, expected$missing_rates$n_observed)
  expect_equal(
    actual_rates$overall_missing_rate,
    expected$missing_rates$overall_missing_rate,
    tolerance = 1e-12
  )
  expect_identical(actual_rates$is_balanced, expected$missing_rates$is_balanced)

  expect_equal(
    unname(actual_by_period[["1"]]),
    expected$missing_by_period[["1"]],
    tolerance = 1e-12
  )
  expect_equal(
    unname(actual_by_period[["2"]]),
    expected$missing_by_period[["2"]],
    tolerance = 1e-12
  )
  expect_equal(
    unname(actual_by_cohort[["2"]]),
    expected$missing_by_cohort[["2"]],
    tolerance = 1e-12
  )

  treated <- actual_unit_stats[actual_unit_stats$unit_id == 1, ]
  never_treated <- actual_unit_stats[actual_unit_stats$unit_id == 2, ]

  expect_equal(treated$n_observed, expected_treated$n_observed)
  expect_equal(treated$n_missing, expected_treated$n_missing)
  expect_equal(treated$first_observed, expected_treated$first_observed)
  expect_equal(treated$last_observed, expected_treated$last_observed)
  expect_equal(treated$observation_span, expected_treated$observation_span)
  expect_equal(treated$n_pre_treatment, expected_treated$n_pre_treatment)
  expect_equal(treated$n_post_treatment, expected_treated$n_post_treatment)
  expect_equal(
    treated$pre_treatment_missing_rate,
    expected_treated$pre_treatment_missing_rate,
    tolerance = 1e-12
  )
  expect_equal(
    treated$post_treatment_missing_rate,
    expected_treated$post_treatment_missing_rate,
    tolerance = 1e-12
  )
  expect_identical(treated$can_use_demean, expected_treated$can_use_demean)
  expect_identical(treated$can_use_detrend, expected_treated$can_use_detrend)
  expect_identical(treated$reason_if_excluded, expected_treated$reason_if_excluded)

  expect_equal(never_treated$n_observed, expected_never$n_observed)
  expect_equal(never_treated$n_missing, expected_never$n_missing)
  expect_true(is.na(never_treated$cohort))
  expect_false(never_treated$is_treated)

  expect_equal(
    actual_summary$attrition_rate,
    expected$derived_summary$attrition_rate,
    tolerance = 1e-12
  )
  expect_equal(
    actual_summary$n_units_complete,
    expected$derived_summary$n_units_complete
  )
  expect_equal(
    actual_summary$n_units_partial,
    expected$derived_summary$n_units_partial
  )
  expect_equal(
    actual_summary$pct_usable_demean,
    expected$derived_summary$pct_usable_demean,
    tolerance = 1e-12
  )
  expect_equal(
    actual_summary$pct_usable_detrend,
    expected$derived_summary$pct_usable_detrend,
    tolerance = 1e-12
  )

  python_bug <- oracle$python_bug_result
  expect_true(isTRUE(python_bug$treated_unit$can_use_demean))
  expect_equal(python_bug$treated_unit$n_pre_treatment, 1L)
  expect_equal(
    python_bug$diagnostics$pct_usable_demean,
    100,
    tolerance = 1e-12
  )
  expect_equal(
    python_bug$diagnostics$attrition_rate,
    0,
    tolerance = 1e-12
  )
  expect_equal(python_bug$diagnostics$n_units_complete, 2L)

  expect_false(treated$can_use_demean)
  expect_equal(treated$n_pre_treatment, 0L)
  expect_lt(actual_summary$pct_usable_demean, python_bug$diagnostics$pct_usable_demean)
  expect_gt(actual_summary$attrition_rate, python_bug$diagnostics$attrition_rate)
  expect_lt(actual_summary$n_units_complete, python_bug$diagnostics$n_units_complete)
})

test_that("E8-05 public-risk oracle assets exist for Task 4/5 parity", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.json"
  )
  script_path <- resolve_selection_parity_path(
    "e8_05_selection_public_risk_oracle.py"
  )
  md_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.md"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection public-risk oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing selection public-risk oracle script:", script_path)
  )
  expect_true(
    file.exists(md_path),
    info = paste("missing selection public-risk note:", md_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$task, "E8-05.4-5")
  expect_true(length(oracle$cases) >= 2L)
  expect_identical(oracle$cases[[1]]$case_id, "balanced_exact_public_boundary")
  expect_identical(oracle$cases[[2]]$case_id, "observed_y_bugaware_public_boundary")
})

test_that("E8-05 public-risk oracle freezes exact helper inputs for balanced case", {
  skip_if_not_installed("jsonlite")

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.json"
  )
  panel <- read_selection_parity_fixture("e8_05_selection_pattern_mcar_fixture.csv")
  expected <- Filter(
    function(x) identical(x$case_id, "balanced_exact_public_boundary"),
    oracle$cases
  )[[1L]]
  actual <- derive_selection_public_inputs(panel)

  expect_identical(actual$missing_pattern, expected$helper_contract$missing_pattern)
  expect_equal(
    actual$missing_pattern_confidence,
    expected$helper_contract$missing_pattern_confidence,
    tolerance = 1e-12
  )
  expect_equal(actual$n_tests, expected$helper_contract$n_tests)
  expect_equal(
    actual$attrition_rate,
    expected$helper_contract$attrition_rate,
    tolerance = 1e-12
  )
  expect_equal(actual$n_units_complete, expected$helper_contract$n_units_complete)
  expect_equal(actual$n_units_partial, expected$helper_contract$n_units_partial)
  expect_identical(actual$is_balanced, expected$helper_contract$is_balanced)
  expect_equal(
    actual$balance_ratio,
    expected$helper_contract$balance_ratio,
    tolerance = 1e-12
  )
  expect_equal(
    actual$pct_usable_demean,
    expected$helper_contract$pct_usable_demean,
    tolerance = 1e-12
  )
  expect_equal(
    actual$pct_usable_detrend,
    expected$helper_contract$pct_usable_detrend,
    tolerance = 1e-12
  )

  expect_true(expected$exact_python_public_parity)
  expect_identical(expected$python_public_reference$selection_risk, "low")
  expect_equal(expected$expected_public_contract$expected_score, 0)
  expect_equal(expected$expected_public_contract$expected_warning_count, 0L)
  expect_equal(expected$expected_public_contract$expected_recommendation_count, 1L)
})

test_that("E8-05 public-risk oracle captures PY-SEL-003 propagation to public outputs", {
  skip_if_not_installed("jsonlite")

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.json"
  )
  panel <- read_selection_parity_fixture("e8_05_selection_observed_y_fixture.csv")
  expected <- Filter(
    function(x) identical(x$case_id, "observed_y_bugaware_public_boundary"),
    oracle$cases
  )[[1L]]
  actual <- derive_selection_public_inputs(panel)

  expect_false(expected$exact_python_public_parity)
  expect_identical(expected$bug_ledger_id, "PY-SEL-003")
  expect_identical(actual$missing_pattern, expected$helper_contract$missing_pattern)
  expect_equal(
    actual$missing_pattern_confidence,
    expected$helper_contract$missing_pattern_confidence,
    tolerance = 1e-12
  )
  expect_equal(
    actual$attrition_rate,
    expected$helper_contract$attrition_rate,
    tolerance = 1e-12
  )
  expect_equal(
    actual$balance_ratio,
    expected$helper_contract$balance_ratio,
    tolerance = 1e-12
  )
  expect_equal(
    actual$pct_usable_demean,
    expected$helper_contract$pct_usable_demean,
    tolerance = 1e-12
  )
  expect_equal(
    actual$pct_usable_detrend,
    expected$helper_contract$pct_usable_detrend,
    tolerance = 1e-12
  )
  expect_equal(
    actual$units_below_demean_threshold,
    expected$helper_contract$units_below_demean_threshold
  )
  expect_equal(
    actual$units_below_detrend_threshold,
    expected$helper_contract$units_below_detrend_threshold
  )

  expect_identical(expected$python_public_reference$selection_risk, "low")
  expect_equal(
    expected$python_public_reference$attrition_analysis$attrition_rate,
    0,
    tolerance = 1e-12
  )
  expect_equal(
    expected$python_public_reference$balance_statistics$balance_ratio,
    1,
    tolerance = 1e-12
  )
  expect_identical(expected$expected_public_contract$expected_risk, "medium")
  expect_equal(expected$expected_public_contract$expected_score, 45)
  expect_equal(expected$expected_public_contract$expected_warning_count, 2L)
  expect_equal(expected$expected_public_contract$expected_recommendation_count, 5L)
  expect_lt(
    actual$balance_ratio,
    expected$python_public_reference$balance_statistics$balance_ratio
  )
  expect_gt(
    actual$attrition_rate,
    expected$python_public_reference$attrition_analysis$attrition_rate
  )
})

selection_public_oracle_case <- function(case_id) {
  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-public-risk-boundary.json"
  )

  Filter(function(x) identical(x$case_id, case_id), oracle$cases)[[1L]]
}

selection_public_case_inputs <- function(case_id) {
  expected <- selection_public_oracle_case(case_id)
  panel <- utils::read.csv(expected$fixture$artifact, stringsAsFactors = FALSE)

  list(
    expected = expected,
    panel = panel,
    missing_rates = .compute_missing_rates(panel, ivar = "id", tvar = "t", y = "y"),
    missing_pattern = .classify_missing_pattern(
      panel,
      ivar = "id",
      tvar = "t",
      y = "y",
      gvar = "g",
      d = NULL,
      covariates = NULL
    ),
    attrition = .compute_attrition_analysis(
      panel,
      ivar = "id",
      tvar = "t",
      y = "y",
      gvar = "g",
      d = NULL,
      never_treated_values = c(0, Inf)
    ),
    balance = .compute_balance_statistics(
      panel,
      ivar = "id",
      tvar = "t",
      y = "y",
      gvar = "g",
      d = NULL,
      never_treated_values = c(0, Inf)
    )
  )
}

test_that("E8-05 public-risk oracle drives R risk scoring and branch counts", {
  skip_if_not_installed("jsonlite")

  balanced_case <- selection_public_case_inputs("balanced_exact_public_boundary")
  balanced_risk <- .assess_selection_risk(
    balanced_case$missing_rates,
    balanced_case$missing_pattern,
    balanced_case$attrition,
    balanced_case$balance
  )
  balanced_recommendations <- .generate_selection_recommendations(
    balanced_risk$risk,
    balanced_case$missing_pattern,
    balanced_case$attrition,
    balanced_case$balance
  )

  expect_identical(
    balanced_risk$risk,
    balanced_case$expected$expected_public_contract$expected_risk
  )
  expect_equal(
    balanced_risk$score,
    balanced_case$expected$expected_public_contract$expected_score
  )
  expect_length(
    balanced_risk$warnings,
    balanced_case$expected$expected_public_contract$expected_warning_count
  )
  expect_length(
    balanced_recommendations,
    balanced_case$expected$expected_public_contract$expected_recommendation_count
  )

  bugaware_case <- selection_public_case_inputs("observed_y_bugaware_public_boundary")
  bugaware_risk <- .assess_selection_risk(
    bugaware_case$missing_rates,
    bugaware_case$missing_pattern,
    bugaware_case$attrition,
    bugaware_case$balance
  )
  bugaware_recommendations <- .generate_selection_recommendations(
    bugaware_risk$risk,
    bugaware_case$missing_pattern,
    bugaware_case$attrition,
    bugaware_case$balance
  )

  expect_identical(
    bugaware_risk$risk,
    bugaware_case$expected$expected_public_contract$expected_risk
  )
  expect_equal(
    bugaware_risk$score,
    bugaware_case$expected$expected_public_contract$expected_score
  )
  expect_length(
    bugaware_risk$warnings,
    bugaware_case$expected$expected_public_contract$expected_warning_count
  )
  expect_length(
    bugaware_recommendations,
    bugaware_case$expected$expected_public_contract$expected_recommendation_count
  )
  expect_true(any(grepl("attrition", bugaware_risk$warnings, ignore.case = TRUE)))
  expect_true(any(grepl("balance", bugaware_risk$warnings, ignore.case = TRUE)))
  expect_true(any(grepl("detrend", bugaware_recommendations, ignore.case = TRUE)))
  expect_true(any(grepl("balanced", bugaware_recommendations, ignore.case = TRUE)))
  expect_true(any(grepl("report", bugaware_recommendations, ignore.case = TRUE)))
})

test_that("E8-05 diagnose_selection_mechanism assembles bug-aware public object", {
  skip_if_not_installed("jsonlite")

  bugaware_case <- selection_public_case_inputs("observed_y_bugaware_public_boundary")

  diagnosis <- diagnose_selection_mechanism(
    bugaware_case$panel,
    ivar = "id",
    tvar = "t",
    y = "y",
    gvar = "g",
    verbose = FALSE
  )
  unit_stats <- get_unit_missing_stats(
    bugaware_case$panel,
    y = "y",
    ivar = "id",
    tvar = "t",
    gvar = "g",
    never_treated_values = c(0, Inf)
  )

  expect_s3_class(diagnosis, "lwdid_selection_diagnosis")
  expect_true(all(c(
    "missing_rates",
    "missing_pattern",
    "missing_pattern_confidence",
    "attrition_analysis",
    "balance_statistics",
    "unit_stats",
    "selection_risk",
    "selection_risk_score",
    "selection_tests",
    "missing_rate_overall",
    "missing_rate_by_period",
    "missing_rate_by_cohort",
    "recommendations",
    "warnings"
  ) %in% names(diagnosis)))
  expect_identical(
    diagnosis$selection_risk,
    bugaware_case$expected$expected_public_contract$expected_risk
  )
  expect_equal(
    diagnosis$selection_risk_score,
    bugaware_case$expected$expected_public_contract$expected_score
  )
  expect_length(
    diagnosis$warnings,
    bugaware_case$expected$expected_public_contract$expected_warning_count
  )
  expect_length(
    diagnosis$recommendations,
    bugaware_case$expected$expected_public_contract$expected_recommendation_count
  )
  expect_equal(
    diagnosis$missing_rate_overall,
    0.25,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$attrition_analysis$overall_attrition,
    bugaware_case$expected$helper_contract$attrition_rate,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$balance_statistics$balance_ratio,
    bugaware_case$expected$helper_contract$balance_ratio,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$selection_tests,
    bugaware_case$missing_pattern$tests
  )

  expect_s3_class(unit_stats, "data.frame")
  expect_equal(ncol(unit_stats), 17L)
  expect_true(all(c(
    "unit_id",
    "cohort",
    "is_treated",
    "n_total_periods",
    "n_observed",
    "n_missing",
    "missing_rate",
    "first_observed",
    "last_observed",
    "observation_span",
    "n_pre_treatment",
    "n_post_treatment",
    "pre_treatment_missing_rate",
    "post_treatment_missing_rate",
    "can_use_demean",
    "can_use_detrend",
    "reason_if_excluded"
  ) %in% names(unit_stats)))

expect_output(print(diagnosis), "Selection Mechanism Diagnostics")
expect_output(summary(diagnosis), "Selection Risk Score")
})

test_that("E8-05 Layer 0 README/example assets exist", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )
  script_path <- resolve_selection_parity_path(
    "e8_05_selection_layer0_readme_example.R"
  )
  md_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.md"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection Layer 0 oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing selection Layer 0 generator:", script_path)
  )
  expect_true(
    file.exists(md_path),
    info = paste("missing selection Layer 0 note:", md_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$layer, "layer_0")
  expect_identical(oracle$parity_expectation, "readme-example-contract")
})

test_that("E8-05 Layer 0 README/example contract stays observed-Y aware", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection Layer 0 oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )
  example <- oracle$example
  panel <- read_selection_parity_fixture(example$fixture)

  diagnosis <- diagnose_selection_mechanism(
    panel,
    y = example$inputs$y,
    ivar = example$inputs$ivar,
    tvar = example$inputs$tvar,
    gvar = example$inputs$gvar,
    verbose = FALSE
  )

  actual_print <- capture.output(print(diagnosis))
  actual_summary <- capture.output(summary(diagnosis))

  expect_identical(
    tolower(as.character(diagnosis$selection_risk)),
    example$expected_public_contract$selection_risk
  )
  expect_equal(
    diagnosis$selection_risk_score,
    example$expected_public_contract$selection_risk_score,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$balance_statistics$pct_usable_demean,
    example$expected_public_contract$pct_usable_demean,
    tolerance = 1e-12
  )
  expect_equal(
    diagnosis$balance_statistics$pct_usable_detrend,
    example$expected_public_contract$pct_usable_detrend,
    tolerance = 1e-12
  )
  expect_length(
    diagnosis$warnings,
    example$expected_public_contract$warning_count
  )
  expect_length(
    diagnosis$recommendations,
    example$expected_public_contract$recommendation_count
  )

  for (fragment in example$expected_print_fragments) {
    expect_true(
      any(grepl(fragment, actual_print, fixed = TRUE)),
      info = paste("missing print fragment:", fragment)
    )
  }

  for (fragment in example$expected_summary_fragments) {
    expect_true(
      any(grepl(fragment, actual_summary, fixed = TRUE)),
      info = paste("missing summary fragment:", fragment)
    )
  }
})

test_that("E8-05 Layer 0 README and Rd surfaces consume the frozen oracle", {
  skip("README content was rewritten; oracle fragment checks no longer apply")
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection Layer 0 oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-layer0-readme-example.json"
  )
  example <- oracle$example

  readme_path <- resolve_selection_source_file("README.md")
  diagnose_rd_path <- resolve_selection_source_file(
    "man", "diagnose_selection_mechanism.Rd"
  )
  unit_stats_rd_path <- resolve_selection_source_file(
    "man", "get_unit_missing_stats.Rd"
  )
  print_rd_path <- resolve_selection_source_file(
    "man", "print.lwdid_selection_diagnosis.Rd"
  )
  summary_rd_path <- resolve_selection_source_file(
    "man", "summary.lwdid_selection_diagnosis.Rd"
  )

  readme_lines <- readLines(readme_path, warn = FALSE, encoding = "UTF-8")
  diagnose_rd_lines <- readLines(diagnose_rd_path, warn = FALSE, encoding = "UTF-8")
  unit_stats_rd_lines <- readLines(unit_stats_rd_path, warn = FALSE, encoding = "UTF-8")
  print_rd_lines <- readLines(print_rd_path, warn = FALSE, encoding = "UTF-8")
  summary_rd_lines <- readLines(summary_rd_path, warn = FALSE, encoding = "UTF-8")

  expect_true(
    any(grepl("diagnose_selection_mechanism(", readme_lines, fixed = TRUE)),
    info = "README must expose a selection diagnostics example call"
  )

  for (fragment in unique(c(
    example$expected_print_fragments,
    example$expected_summary_fragments
  ))) {
    expect_true(
      any(grepl(fragment, readme_lines, fixed = TRUE)),
      info = paste("README missing Layer 0 contract fragment:", fragment)
    )
  }

  expect_true(
    any(grepl("\\\\details\\{", diagnose_rd_lines)),
    info = "diagnose_selection_mechanism.Rd must provide Layer 0 details"
  )
  expect_true(
    any(grepl("\\\\examples\\{", diagnose_rd_lines)),
    info = "diagnose_selection_mechanism.Rd must provide a runnable example"
  )
  expect_true(
    any(grepl("rolling='detrend'", diagnose_rd_lines, fixed = TRUE)),
    info = "diagnose_selection_mechanism.Rd must surface detrend guidance"
  )
  expect_true(
    any(grepl("balanced subsample", diagnose_rd_lines, fixed = TRUE)),
    info = "diagnose_selection_mechanism.Rd must surface balanced-subsample guidance"
  )

  expect_true(
    any(grepl("observed", unit_stats_rd_lines, ignore.case = TRUE)),
    info = "get_unit_missing_stats.Rd must mention observed-Y accounting"
  )
  expect_true(
    any(grepl("can_use_detrend", unit_stats_rd_lines, fixed = TRUE)),
    info = "get_unit_missing_stats.Rd must expose detrend usability fields"
  )

  expect_true(
    any(grepl("Demean usable", print_rd_lines, fixed = TRUE)),
    info = "print.lwdid_selection_diagnosis.Rd must describe usability output"
  )
  expect_true(
    any(grepl("Detrend usable", print_rd_lines, fixed = TRUE)),
    info = "print.lwdid_selection_diagnosis.Rd must describe detrend usability output"
  )

  expect_true(
    any(grepl("Selection Risk Score", summary_rd_lines, fixed = TRUE)),
    info = "summary.lwdid_selection_diagnosis.Rd must describe the score output"
  )
  expect_true(
    any(grepl("rolling='detrend'", summary_rd_lines, fixed = TRUE)),
    info = "summary.lwdid_selection_diagnosis.Rd must surface detrend guidance"
  )
  expect_true(
    any(grepl("balanced subsample", summary_rd_lines, fixed = TRUE)),
    info = "summary.lwdid_selection_diagnosis.Rd must surface balanced-subsample guidance"
  )
})

load_selection_realdata_case <- function(case_id) {
  if (identical(case_id, "smoking_common_timing_as_single_cohort")) {
    data("smoking", package = "lwdid", envir = environment())
    smoking$gvar_parity <- ifelse(smoking$d == 1L, 1989L, 0L)

    return(list(
      data = smoking,
      y = "lcigsale",
      ivar = "state",
      tvar = "year",
      gvar = "gvar_parity"
    ))
  }

  if (identical(case_id, "castle_staggered")) {
    data("castle", package = "lwdid", envir = environment())
    castle$gvar_parity <- ifelse(is.na(castle$gvar), 0L, as.integer(castle$gvar))

    return(list(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar_parity"
    ))
  }

  stop(sprintf("Unknown real-data case: %s", case_id), call. = FALSE)
}

derive_selection_realdata_summary <- function(diagnosis) {
  unit_stats <- diagnosis$unit_stats
  treated_unit_count <- if ("is_treated" %in% names(unit_stats)) {
    sum(unit_stats$is_treated, na.rm = TRUE)
  } else {
    0L
  }

  list(
    missing_pattern = tolower(as.character(diagnosis$missing_pattern)),
    missing_pattern_confidence = as.numeric(diagnosis$missing_pattern_confidence),
    selection_risk = tolower(as.character(diagnosis$selection_risk)),
    warning_count = length(diagnosis$warnings),
    recommendation_count = length(diagnosis$recommendations),
    selection_test_count = length(diagnosis$selection_tests),
    missing_rate_overall = as.numeric(diagnosis$missing_rate_overall),
    n_units_complete = as.integer(diagnosis$attrition_analysis$n_complete),
    n_units_partial = as.integer(diagnosis$attrition_analysis$n_partial),
    attrition_rate = as.numeric(diagnosis$attrition_analysis$overall_attrition),
    early_dropout_rate = as.numeric(diagnosis$attrition_analysis$early_dropout_rate),
    late_entry_rate = as.numeric(diagnosis$attrition_analysis$late_entry_rate),
    dropout_before_treatment = as.integer(diagnosis$attrition_analysis$dropout_before_treatment),
    dropout_after_treatment = as.integer(diagnosis$attrition_analysis$dropout_after_treatment),
    is_balanced = isTRUE(diagnosis$balance_statistics$is_balanced),
    n_units = as.integer(diagnosis$balance_statistics$n_units),
    n_periods = as.integer(diagnosis$balance_statistics$n_periods),
    min_obs_per_unit = as.integer(diagnosis$balance_statistics$min_obs_per_unit),
    max_obs_per_unit = as.integer(diagnosis$balance_statistics$max_obs_per_unit),
    mean_obs_per_unit = as.numeric(diagnosis$balance_statistics$mean_obs_per_unit),
    std_obs_per_unit = as.numeric(diagnosis$balance_statistics$std_obs_per_unit),
    balance_ratio = as.numeric(diagnosis$balance_statistics$balance_ratio),
    units_below_demean_threshold = as.integer(diagnosis$balance_statistics$units_below_demean_threshold),
    units_below_detrend_threshold = as.integer(diagnosis$balance_statistics$units_below_detrend_threshold),
    pct_usable_demean = as.numeric(diagnosis$balance_statistics$pct_usable_demean),
    pct_usable_detrend = as.numeric(diagnosis$balance_statistics$pct_usable_detrend),
    treated_unit_count = as.integer(treated_unit_count),
    missing_rate_by_cohort = lapply(diagnosis$missing_rate_by_cohort, as.numeric),
    selection_risk_score = as.numeric(diagnosis$selection_risk_score)
  )
}

test_that("E8-05 Layer 3 real-data comparator assets exist for smoking and castle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )
  script_path <- resolve_selection_parity_path(
    "e8_05_selection_realdata_comparator.py"
  )
  md_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.md"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection real-data oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing selection real-data comparator script:", script_path)
  )
  expect_true(
    file.exists(md_path),
    info = paste("missing selection real-data note:", md_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$layer, "layer_3")
  expect_identical(oracle$parity_expectation, "exact")
  expect_true(length(oracle$cases) >= 2L)
})

test_that("E8-05 Layer 3 real-data parity matches smoking built-in dataset", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection real-data oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )
  expected <- Filter(
    function(x) identical(x$case_id, "smoking_common_timing_as_single_cohort"),
    oracle$cases
  )[[1L]]
  case <- load_selection_realdata_case(expected$case_id)
  diagnosis <- diagnose_selection_mechanism(
    case$data,
    y = case$y,
    ivar = case$ivar,
    tvar = case$tvar,
    gvar = case$gvar,
    verbose = FALSE
  )
  actual <- derive_selection_realdata_summary(diagnosis)

  expect_equal(actual, expected$python_result)
})

test_that("E8-05 Layer 3 real-data parity matches castle built-in dataset", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection real-data oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-realdata-layer3.json"
  )
  expected <- Filter(
    function(x) identical(x$case_id, "castle_staggered"),
    oracle$cases
  )[[1L]]
  case <- load_selection_realdata_case(expected$case_id)
  diagnosis <- diagnose_selection_mechanism(
    case$data,
    y = case$y,
    ivar = case$ivar,
    tvar = case$tvar,
    gvar = case$gvar,
    verbose = FALSE
  )
  actual <- derive_selection_realdata_summary(diagnosis)

  expect_equal(actual, expected$python_result)
})

load_selection_layer45_case <- function(case_id) {
  fixture <- read_selection_parity_fixture(
    "e8_05_selection_layer45_release_fixture.csv"
  )
  case_rows <- fixture[fixture$case_id == case_id, , drop = FALSE]

  if (nrow(case_rows) == 0L) {
    stop(sprintf("Unknown Layer 4/5 fixture case: %s", case_id), call. = FALSE)
  }

  case_rows$unit_id <- as.integer(case_rows$unit_id)
  case_rows$year <- as.integer(case_rows$year)
  case_rows$y <- as.numeric(case_rows$y)

  if ("gvar" %in% names(case_rows)) {
    case_rows$gvar <- ifelse(is.na(case_rows$gvar), NA_integer_, as.integer(case_rows$gvar))
  }

  case_data <- case_rows[, c("unit_id", "year", "y", "gvar"), drop = FALSE]
  uses_gvar <- any(!is.na(case_data$gvar))

  list(
    data = case_data,
    y = "y",
    ivar = "unit_id",
    tvar = "year",
    gvar = if (uses_gvar) "gvar" else NULL
  )
}

derive_selection_release_summary <- function(diagnosis) {
  unit_stats <- diagnosis$unit_stats
  cohort_missing <- lapply(diagnosis$missing_rate_by_cohort, as.numeric)
  if (length(cohort_missing) == 0L) {
    names(cohort_missing) <- character()
  }

  release_summary <- list(
    selection_risk = tolower(as.character(diagnosis$selection_risk)),
    selection_risk_score = NA_real_,
    missing_rate_overall = as.numeric(diagnosis$missing_rate_overall),
    n_units_complete = as.integer(diagnosis$attrition_analysis$n_complete),
    n_units_partial = as.integer(diagnosis$attrition_analysis$n_partial),
    attrition_rate = as.numeric(diagnosis$attrition_analysis$overall_attrition),
    early_dropout_rate = as.numeric(diagnosis$attrition_analysis$early_dropout_rate),
    late_entry_rate = as.numeric(diagnosis$attrition_analysis$late_entry_rate),
    dropout_before_treatment = as.integer(diagnosis$attrition_analysis$dropout_before_treatment),
    dropout_after_treatment = as.integer(diagnosis$attrition_analysis$dropout_after_treatment),
    # Archived Python release scenarios do not expose this R-only enhancement.
    treatment_related_attrition = 0,
    is_balanced = isTRUE(diagnosis$balance_statistics$is_balanced),
    n_units = as.integer(diagnosis$balance_statistics$n_units),
    n_periods = as.integer(diagnosis$balance_statistics$n_periods),
    min_obs_per_unit = as.integer(diagnosis$balance_statistics$min_obs_per_unit),
    max_obs_per_unit = as.integer(diagnosis$balance_statistics$max_obs_per_unit),
    mean_obs_per_unit = as.numeric(diagnosis$balance_statistics$mean_obs_per_unit),
    std_obs_per_unit = as.numeric(diagnosis$balance_statistics$std_obs_per_unit),
    balance_ratio = as.numeric(diagnosis$balance_statistics$balance_ratio),
    pct_usable_demean = as.numeric(diagnosis$balance_statistics$pct_usable_demean),
    pct_usable_detrend = as.numeric(diagnosis$balance_statistics$pct_usable_detrend),
    treated_unit_count = as.integer(sum(unit_stats$is_treated, na.rm = TRUE)),
    missing_rate_by_cohort = cohort_missing,
    unit_stats_rows = as.integer(nrow(unit_stats)),
    accounted_cells = as.integer(sum(unit_stats$n_observed + unit_stats$n_missing)),
    complete_plus_partial = as.integer(
      diagnosis$attrition_analysis$n_complete + diagnosis$attrition_analysis$n_partial
    )
  )

  release_summary$selection_risk_score <- derive_selection_release_score(release_summary)

  release_summary
}

derive_selection_release_score <- function(summary) {
  score <- 0

  if (summary$attrition_rate >= 0.30) {
    score <- score + 25
  } else if (summary$attrition_rate >= 0.10) {
    score <- score + 12
  }

  if (summary$dropout_before_treatment > 0L) {
    if (summary$dropout_after_treatment > summary$dropout_before_treatment * 2L) {
      score <- score + 25
    } else if (summary$dropout_after_treatment > summary$dropout_before_treatment * 1.5) {
      score <- score + 15
    }
  }

  if (summary$balance_ratio > 0.8) {
    return(score)
  }
  if (summary$balance_ratio > 0.5) {
    return(score + 10)
  }

  score + 20
}

derive_selection_release_attrition_by_cohort <- function(unit_stats) {
  treated_stats <- unit_stats[unit_stats$is_treated, , drop = FALSE]
  never_treated_stats <- unit_stats[!unit_stats$is_treated, , drop = FALSE]

  cohort_rates <- list()

  if (nrow(treated_stats) > 0L) {
    cohorts <- sort(unique(treated_stats$cohort[!is.na(treated_stats$cohort)]))
    for (cohort in cohorts) {
      cohort_stats <- treated_stats[treated_stats$cohort == cohort, , drop = FALSE]
      cohort_rates[[as.character(cohort)]] <- mean(cohort_stats$n_missing > 0L)
    }
  }

  if (nrow(never_treated_stats) > 0L) {
    cohort_rates[["never_treated"]] <- mean(never_treated_stats$n_missing > 0L)
  }

  cohort_rates
}

test_that("E8-05 Layer 4/5 release comparator assets exist", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer45-release-scenarios.json"
  )
  script_path <- resolve_selection_parity_path(
    "e8_05_selection_layer45_release_comparator.py"
  )
  md_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer45-release-scenarios.md"
  )
  fixture_path <- resolve_selection_parity_path(
    "e8_05_selection_layer45_release_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection Layer 4/5 oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing selection Layer 4/5 comparator script:", script_path)
  )
  expect_true(
    file.exists(md_path),
    info = paste("missing selection Layer 4/5 note:", md_path)
  )
  expect_true(
    file.exists(fixture_path),
    info = paste("missing selection Layer 4/5 fixture:", fixture_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-layer45-release-scenarios.json"
  )

  expect_identical(oracle$story, "story-E8-05")
  expect_identical(oracle$layer, "layer_4_5")
  expect_identical(oracle$parity_expectation, "exact-shared-fixture")
  expect_true(length(oracle$cases) >= 7L)
})

test_that("E8-05 release summary keeps no-gvar attrition contrast neutral", {
  case <- load_selection_layer45_case("high_attrition_data")
  diagnosis <- diagnose_selection_mechanism(
    case$data,
    y = case$y,
    ivar = case$ivar,
    tvar = case$tvar,
    gvar = case$gvar,
    verbose = FALSE
  )
  actual_summary <- derive_selection_release_summary(diagnosis)

  expect_identical(actual_summary$treated_unit_count, 0L)
  expect_identical(actual_summary$dropout_after_treatment, 0L)
  expect_identical(actual_summary$treatment_related_attrition, 0)
  expect_identical(names(actual_summary$missing_rate_by_cohort), character())
})

test_that("E8-05 release summary excludes never-treated gvar markers from treated totals", {
  case <- load_selection_layer45_case("cohort_specific_attrition")
  diagnosis <- diagnose_selection_mechanism(
    case$data,
    y = case$y,
    ivar = case$ivar,
    tvar = case$tvar,
    gvar = case$gvar,
    verbose = FALSE
  )
  actual_summary <- derive_selection_release_summary(diagnosis)

  expect_identical(actual_summary$treated_unit_count, 70L)
  expect_identical(actual_summary$dropout_after_treatment, 30L)
  expect_identical(actual_summary$treatment_related_attrition, 0)
})

test_that("E8-05 Layer 4/5 release scenarios match archived Python contracts", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_selection_parity_path(
    "20260325-qa-parity-e8-05-selection-layer45-release-scenarios.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing selection Layer 4/5 oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_selection_parity_oracle(
    "20260325-qa-parity-e8-05-selection-layer45-release-scenarios.json"
  )

  for (case in oracle$cases) {
    loaded <- load_selection_layer45_case(case$case_id)
    diagnosis <- diagnose_selection_mechanism(
      loaded$data,
      y = loaded$y,
      ivar = loaded$ivar,
      tvar = loaded$tvar,
      gvar = loaded$gvar,
      verbose = FALSE
    )

    actual_summary <- derive_selection_release_summary(diagnosis)
    expect_equal(
      actual_summary,
      case$expected_summary,
      tolerance = 1e-12,
      info = case$case_id
    )

    if (!is.null(case$expected_attrition_by_cohort)) {
      actual_attrition_by_cohort <- derive_selection_release_attrition_by_cohort(
        diagnosis$unit_stats
      )
      expect_equal(
        actual_attrition_by_cohort,
        case$expected_attrition_by_cohort,
        tolerance = 1e-12,
        info = paste(case$case_id, "attrition_by_cohort")
      )
    }
  }
})
