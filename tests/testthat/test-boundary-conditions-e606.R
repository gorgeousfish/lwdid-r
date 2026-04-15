library(testthat)

if (!exists("lwdid", mode = "function")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run test-boundary-conditions-e606.R directly.")
  }

  devtools::load_all(
    "/Users/cxy/Desktop/lwdid_r/lwdid-r",
    export_all = FALSE,
    quiet = TRUE
  )
}

quiet_e606_boundary <- function(expr) {
  withCallingHandlers(
    expr,
    lwdid_small_sample = function(w) invokeRestart("muffleWarning"),
    lwdid_overlap = function(w) invokeRestart("muffleWarning"),
    lwdid_data = function(w) invokeRestart("muffleWarning"),
    lwdid_numerical = function(w) invokeRestart("muffleWarning"),
    lwdid_convergence = function(w) invokeRestart("muffleWarning"),
    warning = function(w) invokeRestart("muffleWarning")
  )
}

e606_estimate_propensity_score <- function(...) {
  lwdid:::estimate_propensity_score(...)
}

e606_estimate_outcome_model <- function(...) {
  lwdid:::estimate_outcome_model(...)
}

e606_estimate_ipw <- function(...) {
  lwdid:::estimate_ipw(...)
}

e606_estimate_ipwra <- function(...) {
  lwdid:::estimate_ipwra(...)
}

e606_estimate_psm <- function(...) {
  lwdid::estimate_psm(...)
}

build_e606_minimal_sample_data <- function() {
  data.frame(
    Y = c(0.4, 0.9, 1.1, 1.9, 2.3),
    D = c(0L, 0L, 0L, 1L, 1L),
    X1 = c(-1.2, -0.4, 0.2, 0.9, 1.4)
  )
}

build_e606_drop_trim_failure_data <- function() {
  data.frame(
    Y = c(0.1, 0.3, 0.2, 0.4, 2.2, 2.4, 2.5, 2.7),
    D = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    X1 = c(-8, -7, -6, -5, 5, 6, 7, 8)
  )
}

build_e606_complete_separation_data <- function() {
  data.frame(
    Y = c(0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4),
    D = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    X1 = c(-8, -7, -6, -5, 5, 6, 7, 8)
  )
}

build_e606_highdimensional_data <- function() {
  set.seed(60645)
  n <- 160L
  x_mat <- replicate(10L, stats::rnorm(n))
  colnames(x_mat) <- paste0("X", seq_len(ncol(x_mat)))
  linpred <- 0.3 * x_mat[, 1L] - 0.2 * x_mat[, 2L] + 0.15 * x_mat[, 3L]
  d <- stats::rbinom(n, 1L, stats::plogis(linpred))

  while (sum(d) < 30L || sum(d == 0L) < 30L) {
    d <- stats::rbinom(n, 1L, stats::plogis(linpred))
  }

  data.frame(
    Y = 0.5 + x_mat[, 1L] - 0.5 * x_mat[, 2L] + 0.8 * d +
      stats::rnorm(n, sd = 0.5),
    D = d,
    x_mat,
    check.names = FALSE
  )
}

build_e606_same_ps_data <- function() {
  data.frame(
    Y = c(0.2, 1.2, 0.4, 1.4, 0.6, 1.6, 0.8, 1.8, 1.0, 2.0),
    D = rep(c(0L, 1L), times = 5L),
    X1 = rep(c(-2, -1, 0, 1, 2), each = 2L)
  )
}

build_e606_singular_wls_data <- function() {
  data.frame(
    Y = c(0.2, 0.5, 0.7, 1.0, 1.3, 1.5, 1.8, 2.1),
    D = c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    X1 = c(-2, -1, 0, 1, -2, -1, 0, 1),
    X2 = c(-2, -1, 0, 1, -2, -1, 0, 1)
  )
}

build_e606_zero_se_data <- function() {
  set.seed(60651)
  n <- 40L
  data.frame(
    Y = rep(0, n),
    D = rep(c(0L, 1L), each = n / 2L),
    X1 = stats::rnorm(n)
  )
}

build_e606_bootstrap_drop_data <- function() {
  set.seed(60653)
  n <- 120L
  x1 <- stats::rnorm(n, sd = 1.4)
  x2 <- stats::rnorm(n, sd = 0.8)
  d <- stats::rbinom(n, 1L, stats::plogis(1.1 * x1 - 0.3 * x2))

  while (sum(d) < 25L || sum(d == 0L) < 25L) {
    d <- stats::rbinom(n, 1L, stats::plogis(1.1 * x1 - 0.3 * x2))
  }

  data.frame(
    Y = 0.5 + 0.7 * x1 - 0.4 * x2 + 1.2 * d + stats::rnorm(n, sd = 0.5),
    D = d,
    X1 = x1,
    X2 = x2
  )
}

prepare_e606_ipwra_bootstrap_input <- function(data, trim_method) {
  all_vars <- c("Y", "D", "X1", "X2")
  data_clean <- data[stats::complete.cases(data[, all_vars, drop = FALSE]), , drop = FALSE]

  ps_result <- e606_estimate_propensity_score(
    data_clean,
    "D",
    c("X1", "X2"),
    trim_threshold = 0.1,
    trim_method = trim_method
  )

  if (identical(trim_method, "clip")) {
    return(data_clean)
  }

  valid_mask <- !is.na(ps_result$propensity_scores) &
    ps_result$propensity_scores >= 0.1 &
    ps_result$propensity_scores <= 0.9

  data_clean[valid_mask, , drop = FALSE]
}

manual_e606_ipwra_bootstrap <- function(data, trim_method, seed, n_bootstrap) {
  bootstrap_input <- prepare_e606_ipwra_bootstrap_input(data, trim_method = trim_method)
  n <- nrow(bootstrap_input)
  set.seed(seed)

  att_boot <- numeric(n_bootstrap)
  for (b in seq_len(n_bootstrap)) {
    idx <- sample.int(n, n, replace = TRUE)
    att_boot[[b]] <- lwdid:::.ipwra_point_estimate(
      data = bootstrap_input[idx, , drop = FALSE],
      y = "Y",
      d = "D",
      controls = c("X1", "X2"),
      propensity_controls = c("X1", "X2"),
      trim_threshold = 0.1,
      trim_method = trim_method
    )
  }

  att_valid <- att_boot[!is.na(att_boot)]
  list(
    n_valid = length(att_valid),
    se = stats::sd(att_valid),
    ci = as.numeric(stats::quantile(att_valid, probs = c(0.025, 0.975)))
  )
}

test_that("TC-6.6.43: complete-separation propensity scores emit standardized boundary warnings", {
  df <- build_e606_complete_separation_data()
  boundary_warnings <- list()
  raw_glm_warnings <- character(0)

  ps_result <- withCallingHandlers(
    e606_estimate_propensity_score(
      df,
      "D",
      "X1",
      trim_threshold = 0.01,
      trim_method = "clip"
    ),
    lwdid_numerical = function(w) {
      if (identical(w$detail, "ps_boundary")) {
        boundary_warnings[[length(boundary_warnings) + 1L]] <<- w
      }
      invokeRestart("muffleWarning")
    },
    simpleWarning = function(w) {
      raw_glm_warnings <<- c(raw_glm_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(boundary_warnings, 1L)
  expect_match(conditionMessage(boundary_warnings[[1L]]), "boundary convergence")
  expect_true(all(ps_result$propensity_scores >= 0.01 & ps_result$propensity_scores <= 0.99))
  expect_true(any(ps_result$propensity_scores %in% c(0.01, 0.99)))
  expect_false(
    any(grepl("fitted probabilities numerically 0 or 1 occurred", raw_glm_warnings, fixed = TRUE))
  )
})

test_that("TC-6.6.42: minimum treated/control counts still return finite estimates", {
  df <- build_e606_minimal_sample_data()

  ipw <- quiet_e606_boundary(
    e606_estimate_ipw(df, "Y", "D", "X1")
  )
  ipwra <- quiet_e606_boundary(
    e606_estimate_ipwra(df, "Y", "D", "X1")
  )
  psm <- quiet_e606_boundary(
    e606_estimate_psm(df, "Y", "D", "X1", n_neighbors = 1L, with_replacement = TRUE)
  )

  for (result in list(ipw, ipwra, psm)) {
    expect_true(is.finite(result$att))
    expect_true(is.finite(result$se))
    expect_true(is.finite(result$ci_lower))
    expect_true(is.finite(result$ci_upper))
  }
})

test_that("TC-6.6.44: single-covariate inputs keep the PS and outcome design matrices scalar-valued", {
  df <- build_e606_minimal_sample_data()

  ps_result <- quiet_e606_boundary(
    e606_estimate_propensity_score(df, "D", "X1")
  )
  ipw <- quiet_e606_boundary(
    e606_estimate_ipw(df, "Y", "D", "X1")
  )
  ipwra <- quiet_e606_boundary(
    e606_estimate_ipwra(
      df,
      "Y",
      "D",
      controls = "X1",
      propensity_controls = "X1"
    )
  )
  psm <- quiet_e606_boundary(
    e606_estimate_psm(
      df,
      "Y",
      "D",
      propensity_controls = "X1",
      n_neighbors = 1L,
      with_replacement = TRUE
    )
  )

  expect_identical(ps_result$removed_cols, character(0))
  expect_identical(names(ps_result$coefficients), c("_intercept", "X1"))
  expect_equal(dim(stats::vcov(ps_result$glm_object)), c(2L, 2L))
  expect_identical(colnames(stats::model.matrix(ps_result$glm_object)), c("(Intercept)", "X1"))
  expect_identical(names(ipw$propensity_model_coef), c("_intercept", "X1"))
  expect_identical(names(ipwra$propensity_model_coef), c("_intercept", "X1"))
  expect_identical(names(ipwra$outcome_model_coef), c("_intercept", "X1"))
  expect_true(is.finite(psm$att))
  expect_true(is.finite(psm$se))
})

test_that("TC-6.6.46: shared control-model inputs reject NA controls with the offending column name", {
  df <- build_e606_minimal_sample_data()
  df$X1[[2]] <- NA_real_

  expect_error(
    e606_estimate_propensity_score(df, "D", "X1"),
    class = "lwdid_invalid_param",
    regexp = "X1"
  )
  expect_error(
    e606_estimate_outcome_model(df, "Y", "D", "X1"),
    class = "lwdid_invalid_param",
    regexp = "X1"
  )
  expect_error(
    e606_estimate_ipw(df, "Y", "D", "X1"),
    class = "lwdid_invalid_param",
    regexp = "X1"
  )
})

test_that("TC-6.6.48: an absolute near-zero caliper errors when every treated unit is unmatched", {
  df <- data.frame(
    Y = c(0.2, 0.4, 0.3, 2.0, 2.2),
    D = c(0L, 0L, 0L, 1L, 1L),
    X1 = c(-4, -3, -2, 2.5, 3)
  )

  expect_error(
    quiet_e606_boundary(
      e606_estimate_psm(
        df,
        "Y",
        "D",
        "X1",
        caliper = 1e-10,
        caliper_scale = "absolute",
        with_replacement = TRUE
      )
    ),
    class = "lwdid_estimation_failed",
    regexp = "All treated units failed to match"
  )
})

test_that("TC-6.6.49: drop trimming errors when it removes the estimable sample", {
  df <- build_e606_drop_trim_failure_data()

  err <- tryCatch(
    quiet_e606_boundary(
      e606_estimate_ipw(
        df,
        "Y",
        "D",
        "X1",
        trim_threshold = 0.49,
        trim_method = "drop"
      )
    ),
    error = identity
  )

  expect_s3_class(err, "error")
  expect_true(
    inherits(err, "lwdid_no_treated") || inherits(err, "lwdid_no_control"),
    info = paste("unexpected condition classes:", paste(class(err), collapse = ", "))
  )
  expect_match(conditionMessage(err), "after trimming")
})

test_that("TC-6.6.45: high-dimensional controls keep all advanced estimators numerically stable", {
  df <- build_e606_highdimensional_data()
  controls <- paste0("X", 1:10)

  ipw <- quiet_e606_boundary(
    e606_estimate_ipw(df, "Y", "D", controls)
  )
  ipwra <- quiet_e606_boundary(
    e606_estimate_ipwra(
      df,
      "Y",
      "D",
      controls = controls,
      propensity_controls = controls
    )
  )
  psm <- quiet_e606_boundary(
    e606_estimate_psm(
      df,
      "Y",
      "D",
      propensity_controls = controls,
      n_neighbors = 1L,
      with_replacement = TRUE
    )
  )

  for (result in list(ipw, ipwra, psm)) {
    expect_true(is.finite(result$att))
    expect_true(is.finite(result$se))
    expect_true(is.finite(result$ci_lower))
    expect_true(is.finite(result$ci_upper))
  }

  expect_equal(ipwra$df_resid, nrow(df) - 12L)
  expect_identical(names(ipwra$outcome_model_coef), c("_intercept", controls))
})

test_that("TC-6.6.47: constant propensity scores do not spuriously trigger overlap warnings", {
  df <- build_e606_same_ps_data()
  overlap_warnings <- list()

  ps_result <- e606_estimate_propensity_score(df, "D", "X1")
  ipw <- withCallingHandlers(
    e606_estimate_ipw(df, "Y", "D", "X1"),
    lwdid_overlap = function(w) {
      overlap_warnings[[length(overlap_warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(unique(round(ps_result$propensity_scores, 12L)), 0.5, tolerance = 1e-12)
  expect_length(overlap_warnings, 0L)
  expect_true(is.finite(ipw$att))
  expect_true(is.finite(ipw$se))
  expect_false(is.nan(ipw$se))
})

test_that("TC-6.6.51: zero standard errors yield NaN t statistics instead of Inf", {
  df <- build_e606_zero_se_data()

  ipw <- quiet_e606_boundary(
    e606_estimate_ipw(df, "Y", "D", "X1")
  )
  ipwra <- quiet_e606_boundary(
    e606_estimate_ipwra(df, "Y", "D", "X1")
  )
  psm <- quiet_e606_boundary(
    e606_estimate_psm(df, "Y", "D", "X1")
  )

  for (result in list(ipw, ipwra, psm)) {
    expect_equal(result$att, 0, tolerance = 1e-12)
    expect_equal(result$se, 0, tolerance = 1e-12)
    expect_true(is.nan(result$t_stat))
    expect_true(is.nan(result$pvalue))
  }
})

test_that("TC-6.6.50: singular WLS errors cleanly while singular PS Hessians fall back via ginv", {
  df <- build_e606_singular_wls_data()
  weights <- ifelse(df$D == 1L, 1.0, 0.5)

  expect_error(
    e606_estimate_outcome_model(
      df,
      "Y",
      "D",
      controls = c("X1", "X2"),
      sample_weights = weights
    ),
    class = "lwdid_estimation_failed",
    regexp = "WLS design matrix singular"
  )

  ipwra <- quiet_e606_boundary(
    e606_estimate_ipwra(
      df,
      "Y",
      "D",
      controls = "X1",
      propensity_controls = c("X1", "X2")
    )
  )

  expect_true(is.finite(ipwra$att))
  expect_true(is.finite(ipwra$se))
  expect_true(is.finite(ipwra$ci_lower))
  expect_true(is.finite(ipwra$ci_upper))
})

test_that("TC-6.6.52: mocked degenerate propensity outputs trigger the IPW non-positive weight-sum guard", {
  df <- build_e606_minimal_sample_data()
  mocked_ps <- quiet_e606_boundary(
    e606_estimate_propensity_score(df, "D", "X1")
  )
  mocked_ps$propensity_scores[df$D == 0L] <- 0

  err <- tryCatch(
    testthat::with_mocked_bindings(
      e606_estimate_ipw(df, "Y", "D", "X1"),
      estimate_propensity_score = function(...) mocked_ps,
      .package = "lwdid"
    ),
    error = identity
  )

  expect_s3_class(err, "error")
  expect_true(
    inherits(err, "lwdid_estimation_failed"),
    info = paste("unexpected condition classes:", paste(class(err), collapse = ", "))
  )
  expect_match(conditionMessage(err), "weight sum is non-positive")
})

test_that("TC-6.6.53: IPWRA bootstrap threads trim_method='drop' through bootstrap point estimates", {
  df <- build_e606_bootstrap_drop_data()

  analytical_drop <- quiet_e606_boundary(
    e606_estimate_ipwra(
      df,
      "Y",
      "D",
      controls = c("X1", "X2"),
      propensity_controls = c("X1", "X2"),
      trim_threshold = 0.1,
      trim_method = "drop"
    )
  )

  expect_gt(analytical_drop$n_trimmed, 0L)

  bootstrap_drop <- quiet_e606_boundary(
    e606_estimate_ipwra(
      df,
      "Y",
      "D",
      controls = c("X1", "X2"),
      propensity_controls = c("X1", "X2"),
      trim_threshold = 0.1,
      trim_method = "drop",
      vce = "bootstrap",
      boot_reps = 40L,
      seed = 1353L
    )
  )

  manual_drop <- manual_e606_ipwra_bootstrap(
    data = df,
    trim_method = "drop",
    seed = 1353L,
    n_bootstrap = 40L
  )
  manual_clip <- manual_e606_ipwra_bootstrap(
    data = df,
    trim_method = "clip",
    seed = 1353L,
    n_bootstrap = 40L
  )

  expect_gte(manual_drop$n_valid, 10L)
  expect_true(is.finite(bootstrap_drop$se))
  expect_equal(bootstrap_drop$se, manual_drop$se, tolerance = 1e-12)
  expect_equal(bootstrap_drop$ci_lower, manual_drop$ci[[1L]], tolerance = 1e-12)
  expect_equal(bootstrap_drop$ci_upper, manual_drop$ci[[2L]], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(manual_drop$se, manual_clip$se)))
})
