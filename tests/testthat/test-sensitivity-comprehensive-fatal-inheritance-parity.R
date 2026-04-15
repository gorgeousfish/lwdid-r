fatal_inheritance_oracle_path <- resolve_parity_fixture_path(
  "20260324-qa-parity-e8-03-fatal-inheritance-regression.json"
)

test_that("Task 8.5 parity oracle exists for fatal inheritance regression", {
  expect_true(
    file.exists(fatal_inheritance_oracle_path),
    info = "Run the qa-parity fatal inheritance audit to refresh the oracle."
  )
})

test_that("Task 8.5 parity: common-timing transformation wrapper uses c/N RI p-values", {
  dt <- generate_ct_panel(
    N = 60L,
    T_total = 8L,
    tpost1 = 5L,
    treat_frac = 0.5,
    with_controls = TRUE,
    seed = 101L
  )
  dt$post <- as.integer(dt$time >= 5L)

  original_lwdid <- get("lwdid", asNamespace("lwdid"))
  captured <- list()

  result <- testthat::with_mocked_bindings(
    lwdid_sensitivity(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      d = "d",
      post = "post",
      type = "transformation",
      controls = c("x1", "x2"),
      ri = TRUE,
      rireps = 30L,
      ri_method = "permutation",
      seed = 123L,
      verbose = FALSE
    ),
    lwdid = function(...) {
      out <- original_lwdid(...)
      captured[[length(captured) + 1L]] <<- out
      out
    },
    .package = "lwdid"
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_false(is.null(result$transformation_comparison))
  expect_length(captured, 2L)

  for (fit in captured) {
    expect_false(is.null(fit$ri_distribution))
    expect_gt(length(fit$ri_distribution), 0L)

    expected_pvalue <- sum(abs(fit$ri_distribution) >= abs(fit$att)) /
      length(fit$ri_distribution)

    expect_equal(fit$ri_pvalue, expected_pvalue, tolerance = 0)
  }
})

test_that("Task 8.5 parity: staggered transformation wrapper uses c/N RI p-values", {
  dt <- generate_staggered_panel(
    n_per_cohort = 12L,
    cohorts = c(4L, 6L, 8L),
    n_never_treated = 12L,
    T_total = 10L,
    with_controls = TRUE,
    seed = 202L
  )

  original_lwdid <- get("lwdid", asNamespace("lwdid"))
  captured <- list()

  result <- testthat::with_mocked_bindings(
    suppressWarnings(
      lwdid_sensitivity(
        data = dt,
        y = "y",
        ivar = "id",
        tvar = "time",
        gvar = "gvar",
        type = "transformation",
        controls = "x1",
        aggregate = "overall",
        ri = TRUE,
        rireps = 60L,
        ri_method = "permutation",
        seed = 321L,
        verbose = FALSE
      )
    ),
    lwdid = function(...) {
      out <- original_lwdid(...)
      captured[[length(captured) + 1L]] <<- out
      out
    },
    .package = "lwdid"
  )

  expect_s3_class(result, "lwdid_sensitivity_comprehensive")
  expect_false(is.null(result$transformation_comparison))
  expect_length(captured, 2L)

  for (fit in captured) {
    expect_false(is.null(fit$ri_distribution))
    expect_gt(length(fit$ri_distribution), 0L)

    expected_pvalue <- sum(abs(fit$ri_distribution) >= abs(fit$att)) /
      length(fit$ri_distribution)

    expect_equal(fit$ri_pvalue, expected_pvalue, tolerance = 0)
  }
})

test_that("Task 8.5 parity: comprehensive wrapper preserves staggered FATAL inheritance evidence", {
  skip_if_not(file.exists(fatal_inheritance_oracle_path))
  skip_if_not_installed("jsonlite")

  oracle <- jsonlite::fromJSON(fatal_inheritance_oracle_path, simplifyVector = TRUE)

  dt <- generate_staggered_panel(
    n_per_cohort = 8L,
    cohorts = c(4L, 6L, 8L),
    n_never_treated = 6L,
    T_total = 10L,
    tau_base = 2.0,
    tau_dynamic = 0.5,
    seed = 42L
  )

  no_anticipation_result <- suppressWarnings(suppressMessages(
    lwdid_sensitivity(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      type = "no_anticipation",
      control_group = "not_yet_treated",
      aggregate = "none",
      verbose = FALSE
    )
  ))

  direct_nyt <- suppressWarnings(suppressMessages(
    lwdid(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ra",
      control_group = "not_yet_treated",
      aggregate = "none",
      verbose = "quiet"
    )
  ))

  direct_nt <- suppressWarnings(suppressMessages(
    lwdid(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ra",
      control_group = "never_treated",
      aggregate = "none",
      verbose = "quiet"
    )
  ))

  seen_warnings <- list()
  transformation_result <- withCallingHandlers(
    lwdid_sensitivity(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      type = "transformation",
      control_group = "not_yet_treated",
      aggregate = "overall",
      verbose = FALSE
    ),
    warning = function(w) {
      seen_warnings <<- c(seen_warnings, list(w))
      invokeRestart("muffleWarning")
    }
  )

  direct_overall_demean <- suppressWarnings(suppressMessages(
    lwdid(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ra",
      control_group = "never_treated",
      aggregate = "overall",
      verbose = "quiet"
    )
  ))

  direct_overall_detrend <- suppressWarnings(suppressMessages(
    lwdid(
      data = dt,
      y = "y",
      ivar = "id",
      tvar = "time",
      gvar = "gvar",
      rolling = "detrend",
      estimator = "ra",
      control_group = "never_treated",
      aggregate = "overall",
      verbose = "quiet"
    )
  ))

  expect_identical(oracle$exact_status, "passed")
  expect_identical(oracle$numeric_status, "passed")
  expect_identical(
    oracle$remaining_gap,
    "FATAL-003 story-level RI regression"
  )

  expect_equal(
    no_anticipation_result$baseline_estimate$att,
    oracle$strict_mask_case$wrapper_baseline_att,
    tolerance = 1e-10
  )
  expect_equal(
    direct_nyt$att,
    oracle$strict_mask_case$direct_not_yet_treated_att,
    tolerance = 1e-10
  )
  expect_equal(
    direct_nt$att,
    oracle$strict_mask_case$direct_never_treated_att,
    tolerance = 1e-10
  )
  expect_equal(
    no_anticipation_result$baseline_estimate$att,
    direct_nyt$att,
    tolerance = 1e-10
  )
  expect_gt(
    abs(no_anticipation_result$baseline_estimate$att - direct_nt$att),
    1e-6
  )

  expect_false(is.null(transformation_result$transformation_comparison))
  expect_true(any(vapply(seen_warnings, function(w) {
    inherits(w, "lwdid_control_group_switch")
  }, logical(1))))

  expect_equal(
    transformation_result$transformation_comparison$demean_att,
    oracle$aggregation_case$wrapper_demean_att,
    tolerance = 1e-10
  )
  expect_equal(
    transformation_result$transformation_comparison$demean_se,
    oracle$aggregation_case$wrapper_demean_se,
    tolerance = 1e-10
  )
  expect_equal(
    transformation_result$transformation_comparison$detrend_att,
    oracle$aggregation_case$wrapper_detrend_att,
    tolerance = 1e-10
  )
  expect_equal(
    transformation_result$transformation_comparison$detrend_se,
    oracle$aggregation_case$wrapper_detrend_se,
    tolerance = 1e-10
  )
  expect_equal(
    direct_overall_demean$att,
    oracle$aggregation_case$direct_never_treated_demean_att,
    tolerance = 1e-10
  )
  expect_equal(
    direct_overall_demean$se_att,
    oracle$aggregation_case$direct_never_treated_demean_se,
    tolerance = 1e-10
  )
  expect_equal(
    direct_overall_detrend$att,
    oracle$aggregation_case$direct_never_treated_detrend_att,
    tolerance = 1e-10
  )
  expect_equal(
    direct_overall_detrend$se_att,
    oracle$aggregation_case$direct_never_treated_detrend_se,
    tolerance = 1e-10
  )
  expect_equal(
    transformation_result$transformation_comparison$demean_att,
    direct_overall_demean$att,
    tolerance = 1e-10
  )
  expect_equal(
    transformation_result$transformation_comparison$detrend_att,
    direct_overall_detrend$att,
    tolerance = 1e-10
  )
})
