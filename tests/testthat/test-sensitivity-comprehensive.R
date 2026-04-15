make_ct_sensitivity_comprehensive <- function(
    n_units = 40L,
    n_periods = 10L,
    n_treated = 20L,
    treatment_start = 2006L,
    seed = 812L
) {
  set.seed(seed)

  years <- seq.int(treatment_start - (n_periods - 5L), treatment_start + 4L)
  df <- data.frame(
    id = rep(seq_len(n_units), each = n_periods),
    year = rep(years, times = n_units)
  )

  df$d <- ifelse(df$id <= n_treated, 1, 0)
  df$post <- ifelse(df$year >= treatment_start, 1, 0)

  unit_fe <- rnorm(n_units, sd = 0.6)
  x1_unit <- rnorm(n_units)
  x2_unit <- runif(n_units, min = -0.5, max = 0.5)

  df$x1 <- x1_unit[df$id]
  df$x2 <- x2_unit[df$id]

  time_index <- df$year - min(df$year)
  df$y <- unit_fe[df$id] +
    0.15 * time_index +
    0.4 * df$x1 -
    0.25 * df$x2 +
    1.5 * df$d * df$post +
    rnorm(nrow(df), sd = 0.15)

  df
}

make_ct_sensitivity_short_pre <- function(seed = 411L) {
  set.seed(seed)

  df <- data.frame(
    id = rep(seq_len(20L), each = 3L),
    year = rep(1:3, times = 20L)
  )

  df$d <- ifelse(df$id <= 10L, 1, 0)
  df$post <- ifelse(df$year >= 2L, 1, 0)

  unit_fe <- rnorm(20L, sd = 0.4)
  df$y <- unit_fe[df$id] +
    0.2 * df$year +
    1.0 * df$d * df$post +
    rnorm(nrow(df), sd = 0.05)

  df
}

test_that("E8-03 regression: public wrapper keeps no-anticipation tuning args", {
  sensitivity_args <- names(formals(lwdid_sensitivity))

  expect_true("max_anticipation" %in% sensitivity_args)
  expect_true("detection_threshold" %in% sensitivity_args)
})

test_that("E8-03 regression: comprehensive sensitivity returns all four blocks", {
  df <- make_ct_sensitivity_comprehensive()

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x2"),
    type = "all",
    max_anticipation = 2L,
    detection_threshold = 0.90,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_setequal(
    names(result),
    c(
      "pre_period_result",
      "no_anticipation_result",
      "transformation_comparison",
      "estimator_comparison",
      "overall_assessment",
      "recommendations"
    )
  )
  expect_false(is.null(result$pre_period_result))
  expect_false(is.null(result$no_anticipation_result))
  expect_false(is.null(result$transformation_comparison))
  expect_false(is.null(result$estimator_comparison))
  expect_equal(result$no_anticipation_result$max_anticipation_tested, 2L)
  expect_equal(result$no_anticipation_result$detection_threshold, 0.90)
})

test_that("E8-03 dispatch: wrapper returns direct pre-period and no-anticipation results", {
  df <- make_ct_sensitivity_comprehensive()

  pre_period_result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    type = "pre_period",
    verbose = FALSE
  )

  no_anticipation_result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    type = "no_anticipation",
    max_anticipation = 2L,
    detection_threshold = 0.90,
    verbose = FALSE
  )

  expect_s3_class(pre_period_result, "lwdid_sensitivity")
  expect_identical(pre_period_result$type, "pre_period")
  expect_s3_class(no_anticipation_result, "lwdid_sensitivity")
  expect_identical(no_anticipation_result$type, "no_anticipation")
  expect_equal(no_anticipation_result$max_anticipation_tested, 2L)
  expect_equal(no_anticipation_result$detection_threshold, 0.90)
})

test_that("E8-03 comprehensive omits estimator block without controls", {
  df <- make_ct_sensitivity_comprehensive()

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = NULL,
    type = "all",
    max_anticipation = 2L,
    detection_threshold = 0.90,
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_null(result$estimator_comparison)
  expect_null(result$overall_robust)
  expect_identical(
    result$overall_assessment,
    "Estimates are robust across multiple robustness checks"
  )
  expect_identical(
    result$recommendations,
    "No major robustness issues found"
  )
})

test_that("E8-03 constructor falls back to default comprehensive messages", {
  result <- lwdid:::.new_lwdid_sensitivity_comprehensive(
    overall_assessment = NULL,
    recommendations = NULL
  )

  expect_identical(
    result$overall_assessment,
    "No major robustness issues found"
  )
  expect_identical(
    result$recommendations,
    "No major robustness issues found"
  )
})

test_that("E8-03 edge: transformation comparison skips detrend when only one pre-period exists", {
  df <- make_ct_sensitivity_short_pre()

  expect_warning(
    result <- lwdid_sensitivity(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      type = "transformation",
      verbose = FALSE
    ),
    class = "lwdid_data"
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_null(result$transformation_comparison)
})

test_that("E8-03 edge: comprehensive wrapper is silent when verbose is FALSE", {
  df <- make_ct_sensitivity_comprehensive()

  expect_silent(
    lwdid_sensitivity(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      max_anticipation = 2L,
      detection_threshold = 0.90,
      verbose = FALSE
    )
  )
})

test_that("E8-03 edge: comprehensive collapses to empty object when all analyses fail", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      stop("forced failure")
    },
    .package = "lwdid"
  )

  expect_warning(
    result <- lwdid_sensitivity(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      max_anticipation = 2L,
      detection_threshold = 0.90,
      verbose = FALSE
    ),
    class = "lwdid_data"
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_null(result$pre_period_result)
  expect_null(result$no_anticipation_result)
  expect_null(result$transformation_comparison)
  expect_null(result$estimator_comparison)
})

test_that("TC-8.3.4: plot() reports the ggplot2 dependency clearly", {
  df <- make_ct_sensitivity_comprehensive()

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    type = "pre_period",
    verbose = FALSE
  )

  err <- tryCatch(
    testthat::with_mocked_bindings(
      plot(result),
      requireNamespace = function(...) FALSE,
      .package = "base"
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(
    conditionMessage(err),
    "ggplot2.*required.*plotting"
  )
})

test_that("E8-03 regression: staggered control_group reaches every all-mode sub-analysis", {
  df <- data.frame(
    id = rep(seq_len(6L), each = 5L),
    year = rep(seq_len(5L), times = 6L),
    y = 0
  )
  df$gvar <- rep(c(4L, 4L, 5L, 5L, Inf, Inf), each = 5L)

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      args <- list(...)
      control_group <- if (!is.null(args$control_group)) {
        args$control_group
      } else {
        "default"
      }
      att <- if (identical(control_group, "never_treated")) 11 else 7

      list(
        att = att,
        se_att = 0.1,
        t_stat = att / 0.1,
        pvalue = 0.01,
        ci_lower = att - 0.2,
        ci_upper = att + 0.2,
        n_treated = 2L,
        n_control = 2L,
        df_inference = 8L
      )
    },
    .package = "lwdid"
  )

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    gvar = "gvar",
    type = "all",
    control_group = "never_treated",
    verbose = FALSE
  )

  expect_equal(result$pre_period_result$baseline_spec$att, 11)
  expect_equal(result$no_anticipation_result$baseline_estimate$att, 11)
  expect_equal(result$transformation_comparison$demean_att, 11)
  expect_equal(result$transformation_comparison$detrend_att, 11)
})

test_that("TC-8.3.10: single-step failure warns and preserves remaining comprehensive blocks", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    sensitivity_no_anticipation = function(...) {
      stop("anticipation failed")
    },
    .package = "lwdid"
  )

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    lwdid_sensitivity(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      type = "all",
      verbose = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_false(is.null(result$pre_period_result))
  expect_null(result$no_anticipation_result)
  expect_false(is.null(result$transformation_comparison))
  expect_false(is.null(result$estimator_comparison))
  expect_true(any(grepl("anticipation failed", warnings_seen, fixed = TRUE)))
})

test_that("TC-8.3.14: single-estimator comparison keeps baseline ATT while skipping range", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      args <- list(...)
      est <- if (is.null(args$estimator)) "ra" else args$estimator
      if (!identical(est, "ra")) {
        stop(sprintf("forced %s failure", est))
      }

      list(
        att = 1.23,
        se_att = 0.11,
        t_stat = 11.18,
        pvalue = 0.01,
        ci_lower = 1.0,
        ci_upper = 1.46,
        n_treated = 10L,
        n_control = 10L,
        df_inference = 18L
      )
    },
    .package = "lwdid"
  )

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x2"),
    type = "estimator",
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_equal(result$estimator_comparison$ra, 1.23)
  expect_null(result$estimator_comparison$range)
  expect_null(result$estimator_comparison$rel_range)
  expect_equal(result$estimator_comparison$baseline_att, 1.23)
})

test_that("TC-8.3 zero-baseline guard skips transformation sensitivity issue", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      args <- list(...)
      rolling <- if (is.null(args$rolling)) "demean" else args$rolling
      att <- if (identical(rolling, "detrend")) 0.3 else 1e-12

      list(
        att = att,
        se_att = 0.05,
        t_stat = att / 0.05,
        pvalue = 0.5,
        ci_lower = att - 0.1,
        ci_upper = att + 0.1,
        n_treated = 10L,
        n_control = 10L,
        df_inference = 18L
      )
    },
    .package = "lwdid"
  )

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    type = "transformation",
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_true(is.na(result$transformation_comparison$rel_diff))
  expect_identical(
    result$overall_assessment,
    "Estimates are robust across multiple robustness checks"
  )
  expect_identical(
    result$recommendations,
    "No major robustness issues found"
  )
})

test_that("Task 8.6: estimator comparison drops non-finite ATT values", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      args <- list(...)
      est <- if (is.null(args$estimator)) "ra" else args$estimator

      att <- switch(
        est,
        ra = 1.2,
        ipw = Inf,
        ipwra = 1.6,
        stop("unexpected estimator")
      )

      list(
        att = att,
        se_att = 0.1,
        t_stat = att / 0.1,
        pvalue = 0.01,
        ci_lower = att - 0.2,
        ci_upper = att + 0.2,
        n_treated = 10L,
        n_control = 10L,
        df_inference = 18L
      )
    },
    .package = "lwdid"
  )

  result <- lwdid_sensitivity(
    data = df,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x2"),
    type = "estimator",
    verbose = FALSE
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_named(
    result$estimator_comparison,
    c("ra", "ra_se", "ipwra", "ipwra_se", "range", "rel_range", "baseline_att")
  )
  expect_equal(result$estimator_comparison$ra, 1.2)
  expect_equal(result$estimator_comparison$ipwra, 1.6)
  expect_false("ipw" %in% names(result$estimator_comparison))
  expect_true(is.finite(result$estimator_comparison$range))
  expect_true(is.finite(result$estimator_comparison$rel_range))
})

test_that("Task 8.6: transformation comparison rejects non-finite ATT or SE", {
  df <- make_ct_sensitivity_comprehensive()

  testthat::local_mocked_bindings(
    lwdid = function(...) {
      args <- list(...)
      rolling <- if (is.null(args$rolling)) "demean" else args$rolling

      if (identical(rolling, "demean")) {
        return(list(
          att = 1.1,
          se_att = 0.15,
          t_stat = 7.33,
          pvalue = 0.01,
          ci_lower = 0.8,
          ci_upper = 1.4,
          n_treated = 10L,
          n_control = 10L,
          df_inference = 18L
        ))
      }

      list(
        att = 1.4,
        se_att = Inf,
        t_stat = 14,
        pvalue = 0.01,
        ci_lower = 1.2,
        ci_upper = 1.6,
        n_treated = 10L,
        n_control = 10L,
        df_inference = 18L
      )
    },
    .package = "lwdid"
  )

  expect_warning(
    result <- lwdid_sensitivity(
      data = df,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      type = "transformation",
      verbose = FALSE
    ),
    class = "lwdid_data"
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_null(result$transformation_comparison)
  expect_identical(
    result$overall_assessment,
    "Estimates are robust across multiple robustness checks"
  )
})

test_that("Task 8.1: transformation threshold arithmetic stays source-backed", {
  comparison <- lwdid:::.summarize_transformation_comparison(
    demean_estimate = list(att = 2.0, se_att = 0.2),
    detrend_estimate = list(att = 2.8, se_att = 0.3)
  )

  expect_equal(comparison$difference, 0.8)
  expect_equal(comparison$rel_diff, 0.4)

  assessment <- lwdid:::.compute_overall_assessment(
    pre_period_result = NULL,
    no_anticipation_result = NULL,
    transformation_comparison = comparison,
    estimator_comparison = list(
      ra = 1.5,
      ipw = 2.3,
      ipwra = 1.8,
      range = 0.8,
      rel_range = 0.8 / 1.5,
      baseline_att = 1.5
    )
  )

  expect_identical(
    assessment$overall_assessment,
    "Multiple issues: transformation sensitivity, estimator sensitivity, interpret with caution"
  )
  expect_true(any(grepl("demean", assessment$recommendations, fixed = TRUE)))
  expect_true(any(grepl("inter-estimator difference", assessment$recommendations, fixed = TRUE)))
})

test_that("Task 8.2: zero-baseline guard boundary is explicit in helper arithmetic", {
  near_zero <- lwdid:::.summarize_transformation_comparison(
    demean_estimate = list(att = 1e-11, se_att = 0.2),
    detrend_estimate = list(att = 2e-11, se_att = 0.3)
  )
  expect_equal(near_zero$difference, 1e-11)
  expect_true(is.na(near_zero$rel_diff))

  flagged <- lwdid:::.summarize_transformation_comparison(
    demean_estimate = list(att = 1e-9, se_att = 0.2),
    detrend_estimate = list(att = 2e-9, se_att = 0.3)
  )
  expect_equal(flagged$difference, 1e-9)
  expect_equal(flagged$rel_diff, 1)

  zero_assessment <- lwdid:::.compute_overall_assessment(
    pre_period_result = NULL,
    no_anticipation_result = NULL,
    transformation_comparison = near_zero,
    estimator_comparison = NULL
  )
  expect_identical(
    zero_assessment$overall_assessment,
    "Estimates are robust across multiple robustness checks"
  )

  flagged_assessment <- lwdid:::.compute_overall_assessment(
    pre_period_result = NULL,
    no_anticipation_result = NULL,
    transformation_comparison = flagged,
    estimator_comparison = NULL
  )
  expect_identical(
    flagged_assessment$overall_assessment,
    "Caution: transformation sensitivity detected, see recommendations"
  )
  expect_true(any(grepl("detrend", flagged_assessment$recommendations, fixed = TRUE)))
})
