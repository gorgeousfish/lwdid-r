library(testthat)

make_trend_helper_panel <- function() {
  data.frame(
    firm = rep(1:6, each = 3),
    year = rep(1:3, times = 6),
    cohort = rep(c(5, 3, NA, 0, Inf, 99), each = 3),
    y = seq_len(18),
    stringsAsFactors = FALSE
  )
}

make_trend_common_timing_panel <- function() {
  data.frame(
    firm = rep(1:4, each = 6),
    year = rep(1:6, times = 4),
    y = rep(c(0, 1, 2, 3, 4, 5), times = 4),
    stringsAsFactors = FALSE
  )
}

make_trend_balanced_panel <- function() {
  data.frame(
    firm = rep(1:3, each = 4),
    year = rep(1:4, times = 3),
    y = seq_len(12),
    stringsAsFactors = FALSE
  )
}

make_trend_short_panel <- function() {
  data.frame(
    firm = rep(1:3, each = 6),
    year = rep(1:6, times = 3),
    y = rep(c(0, 1, 0, -1, 0, 1), times = 3),
    stringsAsFactors = FALSE
  )
}

make_trend_seasonal_panel <- function() {
  pattern <- c(0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1)

  data.frame(
    firm = rep(1:4, each = length(pattern)),
    year = rep(seq_along(pattern), times = 4),
    y = rep(pattern, times = 4),
    stringsAsFactors = FALSE
  )
}

make_parallel_trends_common_panel <- function(delta = 0) {
  panel <- expand.grid(
    firm = 1:12,
    year = 1:8,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel <- panel[order(panel$firm, panel$year), ]
  panel$eps <- ((panel$firm * 7 + panel$year * 3) %% 11 - 5) / 10
  panel$y <- with(
    panel,
    2 + firm / 5 + 0.4 * year + eps +
      ifelse(firm <= 6 & year < 4, delta * (year - 1), 0) +
      ifelse(firm <= 6 & year >= 4, 1.0, 0)
  )
  panel
}

make_parallel_trends_staggered_fallback_panel <- function() {
  data.frame(
    firm = rep(1:4, each = 4),
    year = rep(1:4, times = 4),
    cohort = rep(c(4, 4, Inf, Inf), each = 4),
    y = c(
      2, 4, NA, 7,
      4, 6, NA, 9,
      1, 2, NA, 4,
      3, 4, NA, 6
    ),
    stringsAsFactors = FALSE
  )
}

read_parallel_trends_staggered_estimable_fixture <- function() {
  fixture_path <- if (exists("resolve_parity_fixture_path", mode = "function")) {
    resolve_parity_fixture_path("e8_06_trend_public_staggered_estimable_fixture.csv")
  } else {
    candidates <- c(
      testthat::test_path(
        "..", "..", "..",
        "_automation", "test-artifacts", "parity",
        "e8_06_trend_public_staggered_estimable_fixture.csv"
      ),
      file.path(
        "/Users/cxy/Desktop/lwdid_r",
        "_automation", "test-artifacts", "parity",
        "e8_06_trend_public_staggered_estimable_fixture.csv"
      )
    )
    candidates <- unique(normalizePath(candidates, mustWork = FALSE))
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) existing[[1L]] else candidates[[1L]]
  }

  utils::read.csv(
    fixture_path,
    stringsAsFactors = FALSE
  )
}

make_heterogeneous_trend_panel <- function() {
  panel <- expand.grid(
    firm = 1:15,
    year = 1:6,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  panel <- panel[order(panel$firm, panel$year), ]

  cohort_map <- c(rep(5, 5), rep(6, 5), rep(Inf, 5))
  base_map <- rep(seq(10, by = 1, length.out = 5), times = 3)
  slope_map <- c(rep(0.5, 5), rep(4, 5), rep(1, 5))

  panel$cohort <- cohort_map[panel$firm]
  panel$y <- base_map[panel$firm] + slope_map[panel$firm] * panel$year
  panel
}

make_singular_trend_panel <- function() {
  data.frame(
    firm = rep(1:2, each = 3),
    year = rep(1:3, times = 2),
    y = c(1, 2, 3, 2, 3, 4),
    x_year = rep(1:3, times = 2),
    stringsAsFactors = FALSE
  )
}

test_that("E8-06.1: .validate_trend_test_inputs accepts valid inputs", {
  panel <- make_trend_helper_panel()

  expect_invisible(
    .validate_trend_test_inputs(
      panel,
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      method = "placebo"
    )
  )

  expect_invisible(
    .validate_trend_test_inputs(
      make_trend_common_timing_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = NULL,
      method = "joint"
    )
  )
})

test_that("E8-06.1: .validate_trend_test_inputs rejects missing columns and invalid methods", {
  panel <- make_trend_helper_panel()

  expect_error(
    .validate_trend_test_inputs(
      panel,
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "missing_gvar",
      method = "placebo"
    ),
    class = "lwdid_missing_column"
  )

  expect_error(
    .validate_trend_test_inputs(
      panel,
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      method = "unsupported"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("E8-06.1: .get_valid_cohorts excludes never-treated markers and sorts cohorts", {
  panel <- make_trend_helper_panel()

  expect_equal(
    .get_valid_cohorts(
      panel,
      gvar = "cohort",
      ivar = "firm",
      never_treated_values = c(0, 99, Inf)
    ),
    c(3L, 5L)
  )
})

test_that("E8-06.1: .compute_pre_period_range returns zeros when no valid cohorts", {
  panel <- data.frame(
    firm = rep(1:3, each = 4),
    year = rep(1:4, times = 3),
    cohort = rep(c(0, Inf, NA), each = 4),
    y = seq_len(12),
    stringsAsFactors = FALSE
  )

  expect_equal(
    .compute_pre_period_range(
      panel,
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      never_treated_values = c(0, Inf)
    ),
    c(0L, 0L)
  )
})

test_that("E8-06.1: .compute_pre_period_range infers common-timing path without dummy columns", {
  panel <- make_trend_common_timing_panel()

  expect_equal(
    .compute_pre_period_range(
      panel,
      ivar = "firm",
      tvar = "year",
      gvar = NULL,
      never_treated_values = c(0, Inf)
    ),
    c(2L, 2L)
  )
})

test_that("E8-06.1: .check_panel_balance distinguishes balanced and unbalanced panels", {
  balanced <- make_trend_balanced_panel()
  unbalanced <- balanced[-1, ]

  expect_true(.check_panel_balance(balanced, ivar = "firm", tvar = "year"))
  expect_false(.check_panel_balance(unbalanced, ivar = "firm", tvar = "year"))
})

test_that("E8-06.1: .detect_seasonal_patterns respects minimum length and threshold", {
  expect_false(
    .detect_seasonal_patterns(
      make_trend_short_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      threshold = 0.1
    )
  )

  expect_true(
    .detect_seasonal_patterns(
      make_trend_seasonal_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      threshold = 0.1
    )
  )
})

test_that("E8-06.2: .compute_joint_f_test matches the Python Wald-over-k contract", {
  estimates <- list(
    list(att = 0.5, se = 0.1, df = 50L, pvalue = 0.001),
    list(att = -0.3, se = 0.2, df = 40L, pvalue = 0.14),
    list(att = 0.8, se = 0.15, df = 60L, pvalue = 0.002),
    list(att = NA_real_, se = 0.3, df = 20L, pvalue = NA_real_),
    list(att = 0.2, se = 0, df = 30L, pvalue = 0.5)
  )

  result <- .compute_joint_f_test(estimates)

  expect_equal(result$f_stat, 18.564814814814813, tolerance = 1e-10)
  expect_equal(
    result$pvalue,
    stats::pf(18.564814814814813, df1 = 3, df2 = 40, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(result$df, c(3L, 40L))
})

test_that("E8-06.2: .compute_joint_f_test falls back to denominator df 100", {
  estimates <- list(
    list(att = 0.5, se = 0.1, df = 0L, pvalue = 0.001),
    list(att = -0.3, se = 0.2, df = -5L, pvalue = 0.14)
  )

  result <- .compute_joint_f_test(estimates)

  expect_equal(result$df, c(2L, 100L))
})

test_that("E8-06.2: .estimate_placebo_att uses Python-style cell means and variances", {
  placebo_data <- data.frame(
    firm = c(1, 1, 2, 2, 3, 3, 4, 4),
    year = c(1, 2, 1, 2, 1, 2, 1, 2),
    y = c(2, 4, 4, 6, 1, 2, 3, 4),
    post_placebo = c(0, 1, 0, 1, 0, 1, 0, 1),
    treat_cohort = c(1, 1, 1, 1, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )

  result <- .estimate_placebo_att(
    placebo_data,
    y = "y",
    ivar = "firm",
    tvar = "year"
  )

  expect_equal(result$att, 1.0, tolerance = 1e-12)
  expect_equal(result$se, 2.0, tolerance = 1e-12)
  expect_equal(result$df, 4L)
  expect_equal(result$n_treated, 2L)
  expect_equal(result$n_control, 2L)
})

test_that("E8-06.2: lwdid_test_parallel_trends infers common timing and recommends demean", {
  result <- lwdid_test_parallel_trends(
    make_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_parallel_trends")
  expect_equal(result$cohorts, 4L)
  expect_match(
    paste(result$warnings, collapse = " "),
    "Assuming common timing with treatment at period 4",
    fixed = TRUE
  )
  expect_equal(result$gvar, "_gvar_dummy")
  expect_equal(result$recommendation, "demean")
  expect_false(result$reject_null)
  expect_true(length(result$pre_trend_estimates) >= 1L)
  expect_false(any(vapply(
    result$pre_trend_estimates,
    function(estimate) identical(estimate$event_time, -1L),
    logical(1)
  )))
})

test_that("E8-06.2: common-timing placebo keeps all estimable demean pre-periods", {
  result <- lwdid_test_parallel_trends(
    make_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_equal(
    vapply(result$pre_trend_estimates, `[[`, integer(1), "event_time"),
    c(-2L, -3L)
  )
  expect_equal(result$joint_df, c(2L, 16L))
  expect_equal(result$joint_f_stat, 0.23947750362844697, tolerance = 1e-12)
  expect_equal(result$joint_pvalue, 0.7898097985657896, tolerance = 1e-12)
})

test_that("E8-06.2: lwdid_test_parallel_trends recommends detrend when placebo effects reject PT", {
  result <- lwdid_test_parallel_trends(
    make_parallel_trends_common_panel(delta = 1),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_equal(result$recommendation, "detrend")
  expect_true(result$reject_null)
  expect_lt(result$joint_pvalue, 0.05)
})

test_that("E8-06.2: lwdid_test_parallel_trends keeps no-estimate contract when rolling placebo estimates are unavailable", {
  result <- lwdid_test_parallel_trends(
    make_parallel_trends_staggered_fallback_panel(),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = "cohort",
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_length(result$pre_trend_estimates, 0L)
  expect_match(
    paste(result$warnings, collapse = " "),
    "No valid pre-treatment estimates computed. Cannot perform joint test.",
    fixed = TRUE
  )
  expect_equal(result$joint_df, c(0L, 0L))
  expect_true(is.nan(result$joint_f_stat))
  expect_equal(result$joint_pvalue, 1.0, tolerance = 1e-12)
  expect_false(result$reject_null)
  expect_equal(result$recommendation, "demean")
})

test_that("E8-06.2: lwdid_test_parallel_trends falls back to simple DiD when the staggered helper is unavailable", {
  sentinel_estimates <- list(list(
    cohort = 4L,
    period = 2L,
    event_time = -2L,
    att = 0.25,
    se = 0.1,
    pvalue = 0.04,
    df = 3L,
    is_anchor = FALSE,
    n_treated = 2L,
    n_control = 2L
  ))

  unavailable_condition <- structure(
    list(message = "forced unavailable", call = NULL),
    class = c("lwdid_staggered_unavailable", "error", "condition")
  )

  result <- testthat::with_mocked_bindings(
    lwdid_test_parallel_trends(
      make_parallel_trends_staggered_fallback_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      method = "placebo",
      rolling = "demean",
      alpha = 0.05,
      verbose = FALSE
    ),
    estimate_pre_treatment_staggered = function(...) {
      stop(unavailable_condition)
    },
    .estimate_placebo_with_simple_did = function(...) {
      sentinel_estimates
    },
    .package = "lwdid"
  )

  expect_length(result$pre_trend_estimates, 1L)
  expect_equal(result$pre_trend_estimates[[1L]]$event_time, -2L)
  expect_equal(result$pre_trend_estimates[[1L]]$att, 0.25)
  expect_true(any(grepl(
    "Staggered module not available. Using simple 2x2 DiD",
    result$warnings,
    fixed = TRUE
  )))
})

test_that("E8-06.2: lwdid_test_parallel_trends falls back when staggered internals are unavailable", {
  sentinel_estimates <- list(list(
    cohort = 4L,
    period = 2L,
    event_time = -2L,
    att = 0.4,
    se = 0.2,
    pvalue = 0.1,
    df = 4L,
    is_anchor = FALSE,
    n_treated = 2L,
    n_control = 2L
  ))

  result <- testthat::with_mocked_bindings(
    lwdid_test_parallel_trends(
      make_parallel_trends_staggered_fallback_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      method = "placebo",
      rolling = "demean",
      alpha = 0.05,
      verbose = FALSE
    ),
    estimate_pre_treatment_staggered = function(...) {
      stop(simpleError("could not find function \"transform_staggered_demean_pre\""))
    },
    .estimate_placebo_with_simple_did = function(...) {
      sentinel_estimates
    },
    .package = "lwdid"
  )

  expect_length(result$pre_trend_estimates, 1L)
  expect_equal(result$pre_trend_estimates[[1L]]$att, 0.4)
  expect_true(any(grepl(
    "Staggered module not available. Using simple 2x2 DiD",
    result$warnings,
    fixed = TRUE
  )))
})

test_that("E8-06.2: module-unavailable fallback matches the Python simple 2x2 contract", {
  unavailable_condition <- structure(
    list(message = "forced unavailable", call = NULL),
    class = c("lwdid_staggered_unavailable", "error", "condition")
  )

  result <- testthat::with_mocked_bindings(
    lwdid_test_parallel_trends(
      make_parallel_trends_staggered_fallback_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
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

  expect_length(result$pre_trend_estimates, 1L)
  expect_equal(result$pre_trend_estimates[[1L]]$event_time, -2L)
  expect_equal(result$pre_trend_estimates[[1L]]$att, 1.0, tolerance = 1e-12)
  expect_equal(result$pre_trend_estimates[[1L]]$se, 2.0, tolerance = 1e-12)
  expect_equal(result$pre_trend_estimates[[1L]]$df, 4L)
  expect_equal(result$pre_trend_estimates[[1L]]$n_treated, 2L)
  expect_equal(result$pre_trend_estimates[[1L]]$n_control, 2L)
  expect_equal(result$joint_f_stat, 0.25, tolerance = 1e-12)
  expect_equal(result$joint_pvalue, 0.6433299631818633, tolerance = 1e-12)
  expect_equal(result$joint_df, c(1L, 4L))
  expect_false(result$reject_null)
  expect_equal(result$recommendation, "demean")
})

test_that("E8-06.2: staggered placebo keeps the earliest estimable pre-period for each cohort", {
  fixture <- read_parallel_trends_staggered_estimable_fixture()

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

  observed_pairs <- vapply(
    result$pre_trend_estimates,
    function(estimate) paste(estimate$cohort, estimate$event_time, sep = ":"),
    character(1)
  )

  expect_equal(length(result$pre_trend_estimates), 12L)
  expect_equal(result$joint_df, c(12L, 148L))
  expect_true(all(c("4:-3", "6:-5", "8:-7") %in% observed_pairs))
})

test_that("E8-06.3: .estimate_cohort_trend matches pooled OLS slope and short-sample boundary", {
  linear_panel <- data.frame(
    firm = rep(1:3, each = 3),
    year = rep(1:3, times = 3),
    y = c(1, 3, 5, 2, 4, 6, 3, 5, 7),
    stringsAsFactors = FALSE
  )

  result <- .estimate_cohort_trend(
    linear_panel,
    y = "y",
    ivar = "firm",
    tvar = "year",
    controls = NULL
  )

  expect_equal(result$slope, 2, tolerance = 1e-12)
  expect_gt(result$slope_se, 0)
  expect_lt(result$slope_pvalue, 0.01)
  expect_equal(result$n_units, 3L)
  expect_equal(result$n_pre_periods, 3L)

  short_result <- .estimate_cohort_trend(
    linear_panel[1:2, ],
    y = "y",
    ivar = "firm",
    tvar = "year",
    controls = NULL
  )

  expect_true(is.nan(short_result$slope))
  expect_true(is.nan(short_result$slope_se))
  expect_equal(short_result$slope_pvalue, 1)
})

test_that("E8-06.3: .estimate_cohort_trend returns NaN defaults for singular trend designs", {
  singular_result <- .estimate_cohort_trend(
    make_singular_trend_panel(),
    y = "y",
    ivar = "firm",
    tvar = "year",
    controls = "x_year"
  )

  expect_true(is.nan(singular_result$intercept))
  expect_true(is.nan(singular_result$intercept_se))
  expect_true(is.nan(singular_result$slope))
  expect_true(is.nan(singular_result$slope_se))
  expect_equal(singular_result$slope_pvalue, 1)
  expect_equal(singular_result$n_units, 2L)
  expect_equal(singular_result$n_pre_periods, 3L)
  expect_true(is.nan(singular_result$r_squared))
  expect_true(is.nan(singular_result$residual_std))
})

test_that("E8-06.3: lwdid_diagnose_heterogeneous_trends returns public heterogeneity diagnostics", {
  result <- lwdid_diagnose_heterogeneous_trends(
    make_heterogeneous_trend_panel(),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = "cohort",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_heterogeneous_trends")
  expect_false("warnings" %in% names(result))
  expect_true(result$has_heterogeneous_trends)
  expect_equal(result$recommendation, "detrend")
  expect_gt(result$recommendation_confidence, 0.5)
  expect_equal(vapply(result$trend_by_cohort, `[[`, integer(1), "cohort"), c(5L, 6L))
  expect_equal(result$control_group_trend$cohort, 0L)
  expect_identical(result$trend_heterogeneity_test$reject_null, TRUE)

  pair_45 <- Filter(function(x) {
    identical(x$cohort_1, 5L) && identical(x$cohort_2, 6L)
  }, result$trend_differences)
  expect_length(pair_45, 1L)
  expect_true(all(c("slope_1", "slope_2", "slope_diff_se", "df") %in% names(pair_45[[1L]])))
  expect_lt(pair_45[[1L]]$pvalue, 0.05)
})

test_that("E8-06.3: sparse heterogeneity panel keeps the <10 obs default F-test branch", {
  sparse_panel <- data.frame(
    unit = c(rep(1L, 5), rep(2L, 5)),
    time = rep(1:5, times = 2),
    Y = c(1.0, 1.7, 2.9, 5.1, 5.4, 2.0, 2.1, 3.0, 4.4, 5.8),
    first_treat = c(rep(4L, 5), rep(5L, 5)),
    stringsAsFactors = FALSE
  )

  result <- lwdid_diagnose_heterogeneous_trends(
    sparse_panel,
    y = "Y",
    ivar = "unit",
    tvar = "time",
    gvar = "first_treat",
    include_control_group = TRUE,
    alpha = 0.05,
    verbose = FALSE
  )

  expect_false(result$has_heterogeneous_trends)
  expect_equal(
    result$trend_heterogeneity_test,
    list(
      f_stat = 0,
      pvalue = 1,
      df_num = 0L,
      df_den = 0L,
      reject_null = FALSE
    )
  )
  expect_equal(result$recommendation, "demean")
  expect_equal(result$recommendation_confidence, 1.0, tolerance = 1e-12)
  expect_null(result$control_group_trend)
  expect_length(result$trend_by_cohort, 2L)
  expect_length(result$trend_differences, 1L)
  expect_equal(result$trend_differences[[1L]]$df, 1L)
})

test_that("E8-06.3: .compute_pairwise_trend_differences preserves the Welch df floor", {
  trend_1 <- list(
    cohort = 4L,
    slope = 0.95,
    slope_se = 0.14433756729740643,
    n_units = 1L
  )
  trend_2 <- list(
    cohort = 5L,
    slope = 0.81,
    slope_se = 0.20663978319771836,
    n_units = 1L
  )

  result <- .compute_pairwise_trend_differences(
    trend_by_cohort = list(trend_1, trend_2),
    control_group_trend = NULL,
    alpha = 0.05
  )

  slope_diff_se <- sqrt(trend_1$slope_se^2 + trend_2$slope_se^2)
  t_stat <- (trend_1$slope - trend_2$slope) / slope_diff_se
  df_raw <- (trend_1$slope_se^2 + trend_2$slope_se^2)^2 / (
    trend_1$slope_se^4 / max(trend_1$n_units - 1L, 1L) +
      trend_2$slope_se^4 / max(trend_2$n_units - 1L, 1L)
  )

  expect_length(result, 1L)
  expect_lt(df_raw, 2)
  expect_equal(result[[1L]]$df, 1L)
  expect_equal(result[[1L]]$slope_diff_se, slope_diff_se, tolerance = 1e-12)
  expect_equal(
    result[[1L]]$pvalue,
    2 * stats::pt(-abs(t_stat), df = 1L),
    tolerance = 1e-12
  )
})

test_that("E8-06.4: lwdid_recommend_transformation handles common timing without a dummy column", {
  short_common_panel <- data.frame(
    firm = rep(1:4, each = 3),
    year = rep(1:3, times = 4),
    y = c(1, 2, 3, 2, 3, 4, 1, 1, 2, 2, 2, 3),
    stringsAsFactors = FALSE
  )

  result <- lwdid_recommend_transformation(
    short_common_panel,
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_transformation_recommendation")
  expect_equal(result$recommended_method, "demean")
  expect_equal(result$n_pre_periods_min, 1L)
  expect_equal(result$n_pre_periods_max, 1L)
  expect_equal(result$confidence, 1.0, tolerance = 1e-12)
  expect_equal(unname(result$scores[c("demean", "detrend")]), c(1, 0))
  expect_null(result$alternative_method)
  expect_false(isTRUE(result$has_seasonal_pattern))
  expect_true(any(grepl("Detrending requires", result$warnings, fixed = TRUE)))
})

test_that("E8-06.4: equal transformation scores trigger low-confidence warning", {
  result <- .resolve_transformation_recommendation(
    score_demean = 0.5,
    score_detrend = 0.5,
    detrend_feasible = TRUE,
    has_seasonal = FALSE,
    warnings_list = character(0),
    reasons = "Diagnostics are split across methods."
  )

  expect_identical(result$recommended_method, "demean")
  expect_equal(result$confidence, 0.5, tolerance = 1e-12)
  expect_identical(result$confidence_level, "Low")
  expect_identical(result$alternative_method, "detrend")
  expect_equal(
    unname(result$scores[c("demean", "detrend")]),
    c(0.5, 0.5),
    tolerance = 1e-12
  )
  expect_true(any(grepl("Low confidence in recommendation.", result$warnings, fixed = TRUE)))
})

test_that("E8-06.4: summary methods preserve trend diagnostic objects", {
  parallel_result <- lwdid_test_parallel_trends(
    make_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "joint",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )
  heterogeneity_result <- lwdid_diagnose_heterogeneous_trends(
    make_heterogeneous_trend_panel(),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = "cohort",
    alpha = 0.05,
    verbose = FALSE
  )
  recommendation_result <- lwdid_recommend_transformation(
    data.frame(
      firm = rep(1:4, each = 3),
      year = rep(1:3, times = 4),
      y = c(1, 2, 3, 2, 3, 4, 1, 1, 2, 2, 2, 3),
      stringsAsFactors = FALSE
    ),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    run_all_diagnostics = TRUE,
    verbose = FALSE
  )

  expect_identical(summary(parallel_result), parallel_result)
  expect_identical(summary(heterogeneity_result), heterogeneity_result)
  expect_identical(summary(recommendation_result), recommendation_result)
})

test_that("E8-06.4: plot surface returns ggplot objects for trend diagnostics", {
  skip_if_not_installed("ggplot2")

  parallel_result <- lwdid_test_parallel_trends(
    make_parallel_trends_common_panel(delta = 0),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = NULL,
    method = "placebo",
    rolling = "demean",
    alpha = 0.05,
    verbose = FALSE
  )
  heterogeneity_result <- lwdid_diagnose_heterogeneous_trends(
    make_heterogeneous_trend_panel(),
    y = "y",
    ivar = "firm",
    tvar = "year",
    gvar = "cohort",
    alpha = 0.05,
    verbose = FALSE
  )

  expect_s3_class(plot(parallel_result), "ggplot")
  expect_s3_class(plot(heterogeneity_result), "ggplot")
  expect_s3_class(
    plot_cohort_trends(
      make_heterogeneous_trend_panel(),
      y = "y",
      ivar = "firm",
      tvar = "year",
      gvar = "cohort",
      alpha = 0.05
    ),
    "ggplot"
  )
})
