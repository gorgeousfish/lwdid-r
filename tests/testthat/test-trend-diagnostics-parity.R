library(testthat)

resolve_trend_parity_path <- function(filename) {
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
    )
  )

  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  candidates[[1L]]
}

read_trend_parity_oracle <- function(filename) {
  jsonlite::fromJSON(
    resolve_trend_parity_path(filename),
    simplifyVector = FALSE
  )
}

read_trend_parity_fixture <- function(filename) {
  utils::read.csv(
    resolve_trend_parity_path(filename),
    stringsAsFactors = FALSE
  )
}

test_that("E8-06 helper parity assets exist and record the common-timing bug boundary", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_helper_oracle.py"
  )
  staggered_fixture_path <- resolve_trend_parity_path(
    "e8_06_trend_helper_staggered_fixture.csv"
  )
  common_timing_fixture_path <- resolve_trend_parity_path(
    "e8_06_trend_common_timing_fixture.csv"
  )
  short_fixture_path <- resolve_trend_parity_path(
    "e8_06_trend_short_panel_fixture.csv"
  )
  seasonal_fixture_path <- resolve_trend_parity_path(
    "e8_06_trend_seasonal_panel_fixture.csv"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend helper parity oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend helper oracle script:", script_path)
  )
  expect_true(
    file.exists(staggered_fixture_path),
    info = paste("missing trend helper staggered fixture:", staggered_fixture_path)
  )
  expect_true(
    file.exists(common_timing_fixture_path),
    info = paste("missing trend helper common-timing fixture:", common_timing_fixture_path)
  )
  expect_true(
    file.exists(short_fixture_path),
    info = paste("missing trend helper short fixture:", short_fixture_path)
  )
  expect_true(
    file.exists(seasonal_fixture_path),
    info = paste("missing trend helper seasonal fixture:", seasonal_fixture_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )

  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.1")
  expect_identical(
    oracle$cases$common_timing_without_gvar$python_public_api$recommend_transformation$status,
    "bug-confirmed"
  )
  expect_identical(
    oracle$cases$common_timing_without_gvar$python_public_api$recommend_transformation$bug_ledger_id,
    "PY-TRD-004"
  )
})

test_that("E8-06 helper parity matches archived Python helper outputs on staggered fixture", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend helper parity oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  case <- oracle$cases$staggered_helper_baseline
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  expect_equal(
    .get_valid_cohorts(
      fixture,
      gvar = "first_treat",
      ivar = "unit",
      never_treated_values = c(0, Inf)
    ),
    as.integer(unlist(case$python_helper_results$valid_cohorts, use.names = FALSE))
  )

  expect_equal(
    .compute_pre_period_range(
      fixture,
      ivar = "unit",
      tvar = "time",
      gvar = "first_treat",
      never_treated_values = c(0, Inf)
    ),
    c(
      as.integer(case$python_helper_results$pre_period_range$min),
      as.integer(case$python_helper_results$pre_period_range$max)
    )
  )

  expect_identical(
    .check_panel_balance(fixture, ivar = "unit", tvar = "time"),
    isTRUE(case$python_helper_results$is_balanced)
  )
})

test_that("E8-06 common-timing helper path follows paper-backed range instead of Python crash", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend helper parity oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  case <- oracle$cases$common_timing_without_gvar
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  expect_equal(
    .compute_pre_period_range(
      fixture,
      ivar = "unit",
      tvar = "time",
      gvar = NULL,
      never_treated_values = c(0, Inf)
    ),
    c(
      as.integer(case$paper_backed_expected$pre_period_range$min),
      as.integer(case$paper_backed_expected$pre_period_range$max)
    )
  )

  expect_identical(
    case$python_public_api$test_parallel_trends$status,
    "passed"
  )
  expect_identical(
    case$python_public_api$diagnose_heterogeneous_trends$status,
    "passed"
  )
  expect_identical(
    case$python_public_api$recommend_transformation$status,
    "bug-confirmed"
  )
  expect_match(
    case$python_public_api$recommend_transformation$error_message,
    "_dummy",
    fixed = TRUE
  )
})

test_that("E8-06 seasonal helper parity matches Python helper booleans", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend helper parity oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-helper-common-timing.json"
  )

  short_case <- oracle$cases$short_panel_no_seasonality
  short_fixture <- read_trend_parity_fixture(short_case$fixture_csv)
  seasonal_case <- oracle$cases$quarterly_pattern_seasonality
  seasonal_fixture <- read_trend_parity_fixture(seasonal_case$fixture_csv)

  expect_identical(
    .detect_seasonal_patterns(
      short_fixture,
      y = "Y",
      ivar = "unit",
      tvar = "time",
      threshold = short_case$threshold
    ),
    isTRUE(short_case$python_helper_results$has_seasonal)
  )

  expect_identical(
    .detect_seasonal_patterns(
      seasonal_fixture,
      y = "Y",
      ivar = "unit",
      tvar = "time",
      threshold = seasonal_case$threshold
    ),
    isTRUE(seasonal_case$python_helper_results$has_seasonal)
  )
})

test_that("E8-06 public common-timing comparator asset exists and records the Layer 2 match", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-common-timing.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend public comparator oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-common-timing.json"
  )
  case <- oracle$cases$common_timing_placebo_public_api

  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.2")
  expect_identical(case$comparison$status, "matched")
  expect_length(case$comparison$missing_event_times_in_r, 0L)
  expect_null(case$comparison$current_blocker)
  expect_identical(case$python_result$recommendation, "demean")
  expect_identical(case$r_current_behavior$recommendation, "demean")
})

test_that("E8-06 public common-timing comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-common-timing.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend public comparator oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-common-timing.json"
  )
  case <- oracle$cases$common_timing_placebo_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_test_parallel_trends(
    fixture,
    y = "Y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  simplify_pre_trend <- function(estimates) {
    lapply(estimates, function(estimate) {
      list(
        event_time = as.integer(estimate$event_time),
        att = as.numeric(estimate$att),
        se = as.numeric(estimate$se),
        pvalue = as.numeric(estimate$pvalue),
        df = as.integer(estimate$df)
      )
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(result$reject_null, isTRUE(expected$reject_null))
  expect_equal(unname(as.integer(result$joint_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(length(result$pre_trend_estimates), as.integer(expected$n_pre))
  expect_equal(
    simplify_pre_trend(result$pre_trend_estimates),
    expected$pre_trend_estimates,
    tolerance = 1e-10
  )
})

test_that("E8-06 public staggered fallback comparator asset records the Layer 2 match", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_staggered_fallback_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend staggered fallback oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend staggered fallback comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  )
  case <- oracle$cases$staggered_simple_did_fallback_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing trend staggered fallback fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.2")
  expect_identical(case$comparison$status, "matched")
  expect_equal(as.integer(case$python_result$n_pre), 0L)
  expect_equal(as.integer(case$r_current_behavior$n_pre), 0L)
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 public staggered fallback comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend staggered fallback oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-fallback.json"
  )
  case <- oracle$cases$staggered_simple_did_fallback_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_test_parallel_trends(
    fixture,
    y = "Y",
    ivar = "firm",
    tvar = "year",
    gvar = "first_treat",
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  simplify_pre_trend <- function(estimates) {
    lapply(estimates, function(estimate) {
      list(
        event_time = as.integer(estimate$event_time),
        att = as.numeric(estimate$att),
        se = as.numeric(estimate$se),
        pvalue = as.numeric(estimate$pvalue),
        df = as.integer(estimate$df),
        n_treated = as.integer(estimate$n_treated),
        n_control = as.integer(estimate$n_control)
      )
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(result$reject_null, isTRUE(expected$reject_null))
  expect_equal(unname(as.integer(result$joint_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(length(result$pre_trend_estimates), as.integer(expected$n_pre))
  expect_equal(
    unname(as.character(result$warnings)),
    unname(as.character(unlist(expected$warnings, use.names = FALSE)))
  )
  expect_equal(
    simplify_pre_trend(result$pre_trend_estimates),
    expected$pre_trend_estimates,
    tolerance = 1e-10
  )
})

test_that("E8-06 public staggered module-unavailable comparator asset records the Layer 2 match", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-module-unavailable.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_staggered_module_unavailable_comparator.py"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend module-unavailable oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend module-unavailable comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-module-unavailable.json"
  )
  case <- oracle$cases$staggered_module_unavailable_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing trend module-unavailable fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.2")
  expect_identical(case$comparison$status, "matched")
  expect_equal(as.integer(case$python_result$n_pre), 1L)
  expect_equal(as.integer(case$r_current_behavior$n_pre), 1L)
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 public staggered module-unavailable comparator reproduces current R mocked-unavailable snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-module-unavailable.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend module-unavailable oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-module-unavailable.json"
  )
  case <- oracle$cases$staggered_module_unavailable_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)
  unavailable_condition <- structure(
    list(message = "forced unavailable", call = NULL),
    class = c("lwdid_staggered_unavailable", "error", "condition")
  )

  result <- testthat::with_mocked_bindings(
    lwdid_test_parallel_trends(
      fixture,
      y = "Y",
      ivar = "firm",
      tvar = "year",
      gvar = "first_treat",
      method = "placebo",
      rolling = "demean",
      alpha = 0.05,
      verbose = FALSE
    ),
    estimate_pre_treatment_staggered = function(...) {
      stop(unavailable_condition)
    },
    .package = "lwdid"
  )

  expected <- case$r_current_behavior

  simplify_pre_trend <- function(estimates) {
    lapply(estimates, function(estimate) {
      list(
        event_time = as.integer(estimate$event_time),
        att = as.numeric(estimate$att),
        se = as.numeric(estimate$se),
        pvalue = as.numeric(estimate$pvalue),
        df = as.integer(estimate$df),
        n_treated = as.integer(estimate$n_treated),
        n_control = as.integer(estimate$n_control)
      )
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(result$reject_null, isTRUE(expected$reject_null))
  expect_equal(result$joint_f_stat, as.numeric(expected$joint_f_stat), tolerance = 1e-10)
  expect_equal(result$joint_pvalue, as.numeric(expected$joint_pvalue), tolerance = 1e-10)
  expect_equal(unname(as.integer(result$joint_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(length(result$pre_trend_estimates), as.integer(expected$n_pre))
  expect_equal(
    unname(as.character(result$warnings)),
    unname(as.character(unlist(expected$warnings, use.names = FALSE)))
  )
  expect_equal(
    simplify_pre_trend(result$pre_trend_estimates),
    expected$pre_trend_estimates,
    tolerance = 1e-10
  )
})

test_that("E8-06 public staggered estimable comparator asset records the Layer 2 match", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-estimable.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_staggered_estimable_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend staggered estimable oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend staggered estimable comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-estimable.json"
  )
  case <- oracle$cases$staggered_estimable_placebo_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing trend staggered estimable fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.2")
  expect_identical(case$comparison$status, "matched")
  expect_equal(as.integer(case$python_result$n_pre), 12L)
  expect_equal(as.integer(case$r_current_behavior$n_pre), 12L)
  expect_equal(
    as.integer(unlist(case$comparison$python_joint_df)),
    c(12L, 148L)
  )
  expect_equal(
    as.integer(unlist(case$comparison$r_joint_df)),
    c(12L, 148L)
  )
  expect_length(case$comparison$missing_pre_estimates_in_r, 0L)
  expect_length(case$comparison$extra_pre_estimates_in_r, 0L)
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 public staggered estimable comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-staggered-estimable.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend staggered estimable oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-staggered-estimable.json"
  )
  case <- oracle$cases$staggered_estimable_placebo_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_test_parallel_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  simplify_pre_trend <- function(estimates) {
    template_names <- names(expected$pre_trend_estimates[[1L]])

    lapply(estimates, function(estimate) {
      fields <- list(
        cohort = as.integer(estimate$cohort),
        event_time = as.integer(estimate$event_time),
        att = as.numeric(estimate$att),
        se = as.numeric(estimate$se),
        pvalue = as.numeric(estimate$pvalue),
        df = as.integer(estimate$df),
        n_treated = as.integer(estimate$n_treated),
        n_control = as.integer(estimate$n_control)
      )

      fields[template_names]
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(result$reject_null, isTRUE(expected$reject_null))
  expect_equal(unname(as.integer(result$joint_df)), as.integer(unlist(expected$joint_df)))
  expect_equal(length(result$pre_trend_estimates), as.integer(expected$n_pre))
  expect_equal(
    unname(as.character(result$warnings)),
    unname(as.character(unlist(expected$warnings, use.names = FALSE)))
  )
  expect_equal(
    simplify_pre_trend(result$pre_trend_estimates),
    expected$pre_trend_estimates,
    tolerance = 1e-10
  )
})

test_that("E8-06 public diagnostics comparator asset records heterogeneity and recommendation matches", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_diagnostics_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend diagnostics oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend diagnostics comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  heterogeneity_case <- oracle$cases$heterogeneous_trends_public_api
  recommendation_case <- oracle$cases$recommend_transformation_public_api
  fixture_path <- resolve_trend_parity_path(heterogeneity_case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing trend diagnostics fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(heterogeneity_case$comparison$status, "matched")
  expect_identical(recommendation_case$comparison$status, "matched")
  expect_identical(heterogeneity_case$python_result$recommendation, "detrend")
  expect_identical(recommendation_case$python_result$recommended_method, "detrend")
  expect_null(heterogeneity_case$comparison$current_blocker)
  expect_null(recommendation_case$comparison$current_blocker)
})

test_that("E8-06 public heterogeneity comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend diagnostics oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  case <- oracle$cases$heterogeneous_trends_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  simplify_trend <- function(trends) {
    lapply(trends, function(trend) {
      list(
        cohort = as.integer(trend$cohort),
        slope = as.numeric(trend$slope),
        slope_se = as.numeric(trend$slope_se),
        slope_pvalue = as.numeric(trend$slope_pvalue),
        n_units = as.integer(trend$n_units),
        n_pre_periods = as.integer(trend$n_pre_periods)
      )
    })
  }

  simplify_control <- function(trend) {
    if (is.null(trend)) {
      return(NULL)
    }

    list(
      cohort = as.integer(trend$cohort),
      slope = as.numeric(trend$slope),
      slope_se = as.numeric(trend$slope_se),
      slope_pvalue = as.numeric(trend$slope_pvalue),
      n_units = as.integer(trend$n_units),
      n_pre_periods = as.integer(trend$n_pre_periods)
    )
  }

  simplify_differences <- function(differences) {
    lapply(differences, function(difference) {
      list(
        cohort_1 = as.integer(difference$cohort_1),
        cohort_2 = as.integer(difference$cohort_2),
        slope_diff = as.numeric(difference$slope_diff),
        slope_diff_se = as.numeric(difference$slope_diff_se),
        pvalue = as.numeric(difference$pvalue),
        df = as.integer(difference$df)
      )
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(
    isTRUE(result$has_heterogeneous_trends),
    isTRUE(expected$has_heterogeneous_trends)
  )
  expect_equal(
    unname(as.numeric(result$recommendation_confidence)),
    as.numeric(expected$recommendation_confidence),
    tolerance = 1e-10
  )
  expect_equal(
    list(
      f_stat = unname(as.numeric(result$trend_heterogeneity_test$f_stat)),
      pvalue = unname(as.numeric(result$trend_heterogeneity_test$pvalue)),
      df_num = as.integer(result$trend_heterogeneity_test$df_num),
      df_den = as.integer(result$trend_heterogeneity_test$df_den),
      reject_null = isTRUE(result$trend_heterogeneity_test$reject_null)
    ),
    expected$trend_heterogeneity_test,
    tolerance = 1e-10
  )
  expect_equal(
    simplify_trend(result$trend_by_cohort),
    expected$trend_by_cohort,
    tolerance = 1e-10
  )
  expect_equal(
    simplify_control(result$control_group_trend),
    expected$control_group_trend,
    tolerance = 1e-10
  )
  expect_equal(
    simplify_differences(result$trend_differences),
    expected$trend_differences,
    tolerance = 1e-10
  )
})

test_that("E8-06 sparse heterogeneity comparator asset records the <10 obs default branch", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_sparse_heterogeneity_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing sparse heterogeneity oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing sparse heterogeneity comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
  )
  case <- oracle$cases$heterogeneous_trends_sparse_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing sparse heterogeneity fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.3")
  expect_identical(case$comparison$status, "matched")
  expect_equal(
    case$python_result$trend_heterogeneity_test,
    list(
      f_stat = 0,
      pvalue = 1,
      df_num = 0,
      df_den = 0,
      reject_null = FALSE
    )
  )
  expect_equal(
    case$r_current_behavior$trend_heterogeneity_test,
    list(
      f_stat = 0,
      pvalue = 1,
      df_num = 0,
      df_den = 0,
      reject_null = FALSE
    )
  )
  expect_length(case$python_result$trend_differences, 1L)
  expect_length(case$r_current_behavior$trend_differences, 1L)
  expect_equal(case$python_result$trend_differences[[1L]]$df, 1L)
  expect_equal(case$r_current_behavior$trend_differences[[1L]]$df, 1L)
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 sparse heterogeneity comparator is materialized in workspace parity artifacts", {
  expected_oracle <- normalizePath(
    file.path(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
      "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
    ),
    mustWork = FALSE
  )
  expected_script <- normalizePath(
    file.path(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
      "e8_06_trend_public_sparse_heterogeneity_comparator.py"
    ),
    mustWork = FALSE
  )
  expected_fixture <- normalizePath(
    file.path(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
      "e8_06_trend_public_sparse_heterogeneity_fixture.csv"
    ),
    mustWork = FALSE
  )

  expect_true(
    file.exists(expected_oracle),
    info = paste("missing workspace sparse heterogeneity oracle:", expected_oracle)
  )
  expect_true(
    file.exists(expected_script),
    info = paste("missing workspace sparse heterogeneity comparator:", expected_script)
  )
  expect_true(
    file.exists(expected_fixture),
    info = paste("missing workspace sparse heterogeneity fixture:", expected_fixture)
  )
  expect_identical(
    normalizePath(
      resolve_trend_parity_path(
        "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
      ),
      mustWork = FALSE
    ),
    expected_oracle
  )
})

test_that("E8-06 sparse heterogeneity comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing sparse heterogeneity oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-sparse-heterogeneity.json"
  )
  case <- oracle$cases$heterogeneous_trends_sparse_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    include_control_group = TRUE,
    alpha = 0.05,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  simplify_trend <- function(trends) {
    lapply(trends, function(trend) {
      list(
        cohort = as.integer(trend$cohort),
        slope = as.numeric(trend$slope),
        slope_se = as.numeric(trend$slope_se),
        slope_pvalue = as.numeric(trend$slope_pvalue),
        n_units = as.integer(trend$n_units),
        n_pre_periods = as.integer(trend$n_pre_periods)
      )
    })
  }

  simplify_differences <- function(differences) {
    lapply(differences, function(difference) {
      list(
        cohort_1 = as.integer(difference$cohort_1),
        cohort_2 = as.integer(difference$cohort_2),
        slope_diff = as.numeric(difference$slope_diff),
        slope_diff_se = as.numeric(difference$slope_diff_se),
        pvalue = as.numeric(difference$pvalue),
        df = as.integer(difference$df)
      )
    })
  }

  expect_identical(result$recommendation, expected$recommendation)
  expect_identical(
    isTRUE(result$has_heterogeneous_trends),
    isTRUE(expected$has_heterogeneous_trends)
  )
  expect_equal(
    list(
      f_stat = unname(as.numeric(result$trend_heterogeneity_test$f_stat)),
      pvalue = unname(as.numeric(result$trend_heterogeneity_test$pvalue)),
      df_num = as.integer(result$trend_heterogeneity_test$df_num),
      df_den = as.integer(result$trend_heterogeneity_test$df_den),
      reject_null = isTRUE(result$trend_heterogeneity_test$reject_null)
    ),
    expected$trend_heterogeneity_test,
    tolerance = 1e-10
  )
  expect_equal(
    simplify_trend(result$trend_by_cohort),
    expected$trend_by_cohort,
    tolerance = 1e-10
  )
  expect_null(result$control_group_trend)
  expect_equal(
    simplify_differences(result$trend_differences),
    expected$trend_differences,
    tolerance = 1e-10
  )
})

test_that("E8-06 Welch df comparator asset records non-floor numeric coverage", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-story-worker-e8-06-trend-public-welch-df.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_welch_df_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing Welch df oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing Welch df comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-story-worker-e8-06-trend-public-welch-df.json"
  )
  case <- oracle$cases$heterogeneous_trends_welch_df_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing Welch df fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.3")
  expect_identical(case$comparison$status, "matched")
  expect_true(all(unlist(case$comparison$df_matches, use.names = FALSE)))
  expect_equal(
    as.integer(unlist(case$comparison$welch_df_contract, use.names = FALSE)),
    c(71L, 60L, 52L, 87L, 63L, 76L)
  )
  expect_true(all(as.integer(unlist(case$comparison$welch_df_contract, use.names = FALSE)) > 1L))
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 Welch df comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-story-worker-e8-06-trend-public-welch-df.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing Welch df oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-story-worker-e8-06-trend-public-welch-df.json"
  )
  case <- oracle$cases$heterogeneous_trends_welch_df_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_diagnose_heterogeneous_trends(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    alpha = 0.05,
    verbose = FALSE
  )

  expected_differences <- case$r_current_behavior$trend_differences
  actual_differences <- lapply(result$trend_differences, function(difference) {
    list(
      cohort_1 = as.integer(difference$cohort_1),
      cohort_2 = as.integer(difference$cohort_2),
      slope_diff = as.numeric(difference$slope_diff),
      slope_diff_se = as.numeric(difference$slope_diff_se),
      df = as.integer(difference$df)
    )
  })

  expect_equal(actual_differences, expected_differences, tolerance = 1e-10)
  expect_equal(
    vapply(result$trend_differences, `[[`, integer(1), "df"),
    c(71L, 60L, 52L, 87L, 63L, 76L)
  )
})

test_that("E8-06 public recommendation comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend diagnostics oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-diagnostics.json"
  )
  case <- oracle$cases$recommend_transformation_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_recommend_transformation(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  expect_identical(result$recommended_method, expected$recommended_method)
  expect_equal(
    unname(as.numeric(result$confidence)),
    as.numeric(expected$confidence),
    tolerance = 1e-10
  )
  expect_identical(result$confidence_level, expected$confidence_level)
  expect_equal(result$n_pre_periods_min, as.integer(expected$n_pre_periods_min))
  expect_equal(result$n_pre_periods_max, as.integer(expected$n_pre_periods_max))
  expect_identical(
    isTRUE(result$has_seasonal_pattern),
    isTRUE(expected$has_seasonal_pattern)
  )
  expect_identical(
    isTRUE(result$is_balanced_panel),
    isTRUE(expected$is_balanced_panel)
  )
  expect_identical(result$alternative_method, expected$alternative_method)
  expect_identical(result$alternative_reason, expected$alternative_reason)
  expect_equal(unname(as.character(result$warnings)), unname(as.character(expected$warnings)))
  expect_equal(unname(as.character(result$reasons)), unname(as.character(expected$reasons)))
  expect_equal(
    unname(as.numeric(result$scores)),
    as.numeric(unlist(expected$scores, use.names = FALSE)),
    tolerance = 1e-10
  )
})

test_that("E8-06 public seasonal recommendation comparator asset records the Layer 2 match", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-seasonal.json"
  )
  script_path <- resolve_trend_parity_path(
    "e8_06_trend_public_seasonal_comparator.py"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend seasonal oracle:", oracle_path)
  )
  expect_true(
    file.exists(script_path),
    info = paste("missing trend seasonal comparator:", script_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-seasonal.json"
  )
  case <- oracle$cases$recommend_transformation_seasonal_public_api
  fixture_path <- resolve_trend_parity_path(case$fixture_csv)

  expect_true(
    file.exists(fixture_path),
    info = paste("missing trend seasonal fixture:", fixture_path)
  )
  expect_identical(oracle$story, "story-E8-06")
  expect_identical(oracle$task, "E8-06.4")
  expect_identical(case$comparison$status, "matched")
  expect_identical(case$python_result$recommended_method, "demeanq")
  expect_identical(case$r_current_behavior$recommended_method, "demeanq")
  expect_identical(case$r_current_behavior$alternative_method, "detrendq")
  expect_null(case$comparison$current_blocker)
})

test_that("E8-06 public seasonal recommendation comparator matches current R snapshot", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_trend_parity_path(
    "20260325-qa-parity-e8-06-trend-public-seasonal.json"
  )
  expect_true(
    file.exists(oracle_path),
    info = paste("missing trend seasonal oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_trend_parity_oracle(
    "20260325-qa-parity-e8-06-trend-public-seasonal.json"
  )
  case <- oracle$cases$recommend_transformation_seasonal_public_api
  fixture <- read_trend_parity_fixture(case$fixture_csv)

  result <- lwdid_recommend_transformation(
    fixture,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    run_all_diagnostics = FALSE,
    verbose = FALSE
  )

  expected <- case$r_current_behavior

  expect_identical(result$recommended_method, expected$recommended_method)
  expect_equal(
    unname(as.numeric(result$confidence)),
    as.numeric(expected$confidence),
    tolerance = 1e-10
  )
  expect_identical(result$confidence_level, expected$confidence_level)
  expect_equal(result$n_pre_periods_min, as.integer(expected$n_pre_periods_min))
  expect_equal(result$n_pre_periods_max, as.integer(expected$n_pre_periods_max))
  expect_identical(
    isTRUE(result$has_seasonal_pattern),
    isTRUE(expected$has_seasonal_pattern)
  )
  expect_identical(
    isTRUE(result$is_balanced_panel),
    isTRUE(expected$is_balanced_panel)
  )
  expect_identical(result$alternative_method, expected$alternative_method)
  expect_identical(result$alternative_reason, expected$alternative_reason)
  expect_equal(unname(as.character(result$warnings)), unname(as.character(expected$warnings)))
  expect_equal(unname(as.character(result$reasons)), unname(as.character(expected$reasons)))
  expect_equal(
    unname(as.numeric(result$scores)),
    as.numeric(unlist(expected$scores, use.names = FALSE)),
    tolerance = 1e-10
  )
})
