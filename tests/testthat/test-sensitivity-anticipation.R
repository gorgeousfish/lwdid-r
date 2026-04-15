# ===========================================================================
# test-sensitivity-anticipation.R
# Unit tests for sensitivity_no_anticipation() (Story E8-02)
# ===========================================================================

# === Group 1: .detect_anticipation_effects() unit tests (TC-01~08) ===

# TC-01: No anticipation — stable ATTs don't trigger detection
test_that("TC-01: no anticipation with stable ATTs", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.0, converged = TRUE),
    list(excluded_periods = 1L, att = 1.02, converged = TRUE),
    list(excluded_periods = 2L, att = 0.98, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "none_detected")
  expect_equal(result$recommended_exclusion, 0L)
})

# TC-02: Method 1 coefficient change detection
test_that("TC-02: coefficient change method detects anticipation", {
  estimates <- list(
    list(excluded_periods = 0L, att = 0.5, converged = TRUE),
    list(excluded_periods = 1L, att = 0.8, converged = TRUE),
    list(excluded_periods = 2L, att = 0.9, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_true(result$detected)
  expect_equal(result$method, "coefficient_change")
  expect_equal(result$recommended_exclusion, 1L)
})

# TC-03: Direction check prevents false positive
test_that("TC-03: ATT magnitude decreasing does not trigger", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.0, converged = TRUE),
    list(excluded_periods = 1L, att = 0.5, converged = TRUE),
    list(excluded_periods = 2L, att = 0.3, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "none_detected")
})

# TC-04: Method 2 trend break detection
test_that("TC-04: trend break method detects anticipation", {
  estimates <- list(
    list(excluded_periods = 0L, att = 0.5, converged = TRUE),
    list(excluded_periods = 1L, att = 0.55, converged = TRUE),
    list(excluded_periods = 2L, att = 0.9, converged = TRUE),
    list(excluded_periods = 3L, att = 0.95, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.95)
  expect_true(result$detected)
  expect_equal(result$method, "trend_break")
  expect_equal(result$recommended_exclusion, 2L)
})

# TC-04: Method 2 trend break detection
test_that("detect_anticipation: trend break method detects", {
  estimates <- list(
    list(excluded_periods = 0L, att = 0.5, converged = TRUE),
    list(excluded_periods = 1L, att = 0.55, converged = TRUE),
    list(excluded_periods = 2L, att = 0.9, converged = TRUE),
    list(excluded_periods = 3L, att = 0.95, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  # High threshold so method 1 doesn't trigger
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.95)
  expect_true(result$detected)
  expect_equal(result$method, "trend_break")
  expect_equal(result$recommended_exclusion, 2L)
})

# TC-05: Insufficient valid estimates
test_that("detect_anticipation: insufficient data with < 2 valid", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.0, converged = TRUE),
    list(excluded_periods = 1L, att = NA_real_, converged = FALSE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "insufficient_data")
})

# TC-06: Baseline failed
test_that("detect_anticipation: baseline failed returns baseline_failed", {
  estimates <- list(
    list(excluded_periods = 0L, att = NA_real_, converged = FALSE),
    list(excluded_periods = 1L, att = 1.0, converged = TRUE),
    list(excluded_periods = 2L, att = 1.2, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "baseline_failed")
})

# TC-07: Baseline near zero skips method 1, method 2 still runs
test_that("detect_anticipation: baseline near zero skips method 1, method 2 runs", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1e-15, converged = TRUE),
    list(excluded_periods = 1L, att = 0.3, converged = TRUE),
    list(excluded_periods = 2L, att = 0.8, converged = TRUE),
    list(excluded_periods = 3L, att = 0.85, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_true(result$detected)
  expect_equal(result$method, "trend_break")
})

# TC-08: Baseline near zero and no trend break
test_that("detect_anticipation: baseline near zero no trend returns none_detected", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1e-15, converged = TRUE),
    list(excluded_periods = 1L, att = 0.5, converged = TRUE),
    list(excluded_periods = 2L, att = 0.3, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "none_detected")
})


# TC-05: Insufficient valid estimates
test_that("TC-05: insufficient valid estimates returns insufficient_data", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1.0, converged = TRUE),
    list(excluded_periods = 1L, att = NA_real_, converged = FALSE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "insufficient_data")
})

# TC-06: Baseline failed
test_that("TC-06: baseline failed returns baseline_failed", {
  estimates <- list(
    list(excluded_periods = 0L, att = NA_real_, converged = FALSE),
    list(excluded_periods = 1L, att = 1.0, converged = TRUE),
    list(excluded_periods = 2L, att = 1.2, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "baseline_failed")
})

# TC-07: Baseline near zero skips method 1, method 2 still runs
test_that("TC-07: baseline near zero skips method 1, method 2 detects", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1e-15, converged = TRUE),
    list(excluded_periods = 1L, att = 0.3, converged = TRUE),
    list(excluded_periods = 2L, att = 0.8, converged = TRUE),
    list(excluded_periods = 3L, att = 0.85, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_true(result$detected)
  expect_equal(result$method, "trend_break")
})

# TC-08: Baseline near zero and no trend break
test_that("TC-08: baseline near zero no trend break returns none_detected", {
  estimates <- list(
    list(excluded_periods = 0L, att = 1e-15, converged = TRUE),
    list(excluded_periods = 1L, att = 0.5, converged = TRUE),
    list(excluded_periods = 2L, att = 0.3, converged = TRUE)
  )
  baseline <- estimates[[1L]]
  result <- lwdid:::.detect_anticipation_effects(estimates, baseline, 0.10)
  expect_false(result$detected)
  expect_equal(result$method, "none_detected")
})


# ============================================================================
# Group 2: .filter_excluding_periods() unit tests (TC-09 to TC-12a)
# ============================================================================

# TC-09: CT time recoding after exclusion
test_that("filter_excluding_periods: CT recodes time to consecutive integers", {
  ct_data <- data.frame(
    id = rep(1:2, each = 6),
    time = rep(1:6, 2),
    post = rep(c(0, 0, 0, 0, 1, 1), 2),
    y = rnorm(12)
  )
  # exclude=1: remove last pre-period (time=4)
  filtered <- lwdid:::.filter_excluding_periods(ct_data, "id", "time", NULL, "post", 1L)
  remaining_times <- sort(unique(filtered$time))
  expect_equal(remaining_times, 1:5)  # consecutive integers

  # exclude=2: remove last 2 pre-periods (time=3,4)
  filtered2 <- lwdid:::.filter_excluding_periods(ct_data, "id", "time", NULL, "post", 2L)
  remaining_times2 <- sort(unique(filtered2$time))
  expect_equal(remaining_times2, 1:4)  # consecutive integers
})

# TC-10: CT overflow warning
test_that("filter_excluding_periods: CT overflow warns and returns original", {
  ct_data <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, 2),
    post = rep(c(0, 0, 1, 1), 2),
    y = rnorm(8)
  )
  expect_warning(
    filtered <- lwdid:::.filter_excluding_periods(ct_data, "id", "time", NULL, "post", 2L),
    class = "lwdid_data"
  )
  expect_equal(nrow(filtered), nrow(ct_data))
})

# TC-11: Staggered excludes cohort-specific periods
test_that("filter_excluding_periods: Staggered excludes cohort periods", {
  stag_data <- data.frame(
    id = rep(1:4, each = 7),
    time = rep(1:7, 4),
    gvar = rep(c(5, 5, Inf, Inf), each = 7),
    y = rnorm(28)
  )
  # exclude=2: cohort g=5 excludes times {3,4}
  filtered <- lwdid:::.filter_excluding_periods(stag_data, "id", "time", "gvar", NULL, 2L)
  treated_times <- sort(unique(filtered$time[filtered$gvar == 5]))
  expect_false(3 %in% treated_times)
  expect_false(4 %in% treated_times)
  expect_true(all(c(1, 2, 5, 6, 7) %in% treated_times))
})

# TC-12: Never-treated (gvar=Inf) unaffected
test_that("filter_excluding_periods: gvar=Inf unaffected by exclusion", {
  stag_data <- data.frame(
    id = rep(1:4, each = 7),
    time = rep(1:7, 4),
    gvar = rep(c(5, 5, Inf, Inf), each = 7),
    y = rnorm(28)
  )
  filtered <- lwdid:::.filter_excluding_periods(stag_data, "id", "time", "gvar", NULL, 2L)
  never_times <- sort(unique(filtered$time[is.infinite(filtered$gvar)]))
  expect_equal(never_times, 1:7)  # all periods retained
})

# TC-12a: Never-treated (gvar=NA) unaffected — critical R-specific test
test_that("filter_excluding_periods: gvar=NA unaffected by exclusion", {
  stag_data <- data.frame(
    id = rep(1:4, each = 7),
    time = rep(1:7, 4),
    gvar = rep(c(5, 5, NA, NA), each = 7),
    y = rnorm(28)
  )
  filtered <- lwdid:::.filter_excluding_periods(stag_data, "id", "time", "gvar", NULL, 2L)

  # gvar=NA units retain all periods
  na_rows <- filtered[is.na(filtered$gvar), ]
  na_times <- sort(unique(na_rows$time))
  expect_equal(na_times, 1:7)
  expect_equal(nrow(na_rows), 14L)

  # Treated units correctly excluded times {3,4}
  treated_rows <- filtered[!is.na(filtered$gvar) & filtered$gvar == 5, ]
  treated_times <- sort(unique(treated_rows$time))
  expect_false(3 %in% treated_times)
  expect_false(4 %in% treated_times)
  expect_true(all(c(1, 2, 5, 6, 7) %in% treated_times))
})

# Short-circuit: exclude=0 returns original data
test_that("filter_excluding_periods: exclude=0 returns original data", {
  ct_data <- data.frame(
    id = rep(1:2, each = 4),
    time = rep(1:4, 2),
    post = rep(c(0, 0, 1, 1), 2),
    y = rnorm(8)
  )
  filtered <- lwdid:::.filter_excluding_periods(ct_data, "id", "time", NULL, "post", 0L)
  expect_equal(nrow(filtered), nrow(ct_data))
  expect_equal(filtered$time, ct_data$time)
})


# ============================================================================
# Group 3: .get_max_pre_periods() unit tests (TC-13, TC-14)
# ============================================================================

# TC-13: Staggered returns min cohort pre-periods
test_that("TC-13: get_max_pre_periods staggered min cohort", {
  stag_data <- data.frame(
    id = rep(1:3, each = 8),
    time = rep(1:8, 3),
    gvar = rep(c(5, 7, Inf), each = 8)
  )
  result <- lwdid:::.get_max_pre_periods(
    stag_data, "id", "time", "gvar", NULL
  )
  # min(sum(1:8 < 5), sum(1:8 < 7)) = min(4, 6) = 4
  expect_equal(result, 4L)
})

# TC-14: CT returns unique pre-treatment period count
test_that("TC-14: get_max_pre_periods CT pre-period count", {
  ct_data <- data.frame(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    post = rep(c(0, 0, 0, 0, 1, 1), 3)
  )
  result <- lwdid:::.get_max_pre_periods(
    ct_data, "id", "time", NULL, "post"
  )
  expect_equal(result, 4L)
})


# ============================================================================
# Group 4: End-to-end tests (TC-15 to TC-18)
# ============================================================================

# TC-15: CT no-anticipation DGP end-to-end
test_that("TC-15: CT no-anticipation DGP end-to-end", {
  set.seed(42)
  n_units <- 40
  n_periods <- 10
  treat_time <- 6
  n_treated <- 20
  alpha_i <- rnorm(n_units)
  delta_t <- seq(0, 0.9, length.out = n_periods)

  test_data <- expand.grid(id = 1:n_units, time = 1:n_periods)
  test_data$treat <- as.integer(test_data$id <= n_treated)
  test_data$post <- as.integer(test_data$time >= treat_time)
  test_data$y <- alpha_i[test_data$id] +
    delta_t[test_data$time] +
    2.0 * test_data$treat * test_data$post +
    rnorm(nrow(test_data), sd = 0.5)

  result <- lwdid::sensitivity_no_anticipation(
    data = test_data, y = "y",
    ivar = "id", tvar = "time",
    d = "treat", post = "post",
    rolling = "demean",
    max_anticipation = 3L, verbose = FALSE
  )
  expect_s3_class(result, "lwdid_sensitivity")
  expect_equal(result$type, "no_anticipation")
  expect_false(result$anticipation_detected)

  # Numerical: converged ATTs should be near 2.0
  converged <- Filter(
    function(e) isTRUE(e$converged), result$estimates
  )
  for (est in converged) {
    expect_true(abs(est$att - 2.0) < 1.0)
  }
})

# TC-16: CT with anticipation DGP end-to-end
test_that("TC-16: CT with anticipation DGP end-to-end", {
  set.seed(42)
  n_units <- 40
  n_periods <- 10
  treat_time <- 6
  n_treated <- 20
  alpha_i <- rnorm(n_units)
  delta_t <- seq(0, 0.9, length.out = n_periods)

  test_data <- expand.grid(id = 1:n_units, time = 1:n_periods)
  test_data$treat <- as.integer(test_data$id <= n_treated)
  test_data$post <- as.integer(test_data$time >= treat_time)
  anticipation <- ifelse(
    test_data$treat == 1 &
      test_data$time == (treat_time - 1),
    1.5, 0
  )
  test_data$y <- alpha_i[test_data$id] +
    delta_t[test_data$time] +
    3.0 * test_data$treat * test_data$post +
    anticipation +
    rnorm(nrow(test_data), sd = 0.3)

  result <- lwdid::sensitivity_no_anticipation(
    data = test_data, y = "y",
    ivar = "id", tvar = "time",
    d = "treat", post = "post",
    rolling = "demean",
    max_anticipation = 3L, verbose = FALSE
  )
  expect_s3_class(result, "lwdid_sensitivity")
  expect_true(result$anticipation_detected)
  expect_gte(result$recommended_exclusion, 1L)
})

# TC-17: Staggered end-to-end (no anticipation)
test_that("TC-17: Staggered no-anticipation end-to-end", {
  set.seed(42)
  n_per_group <- 10
  n_periods <- 10
  ids <- 1:(3 * n_per_group)
  stag_data <- expand.grid(id = ids, time = 1:n_periods)
  stag_data$gvar <- rep(
    c(rep(5, n_per_group),
      rep(7, n_per_group),
      rep(Inf, n_per_group)),
    each = n_periods
  )
  alpha_i <- rnorm(length(ids))
  treated_post <- as.numeric(
    stag_data$gvar != Inf &
      stag_data$time >= stag_data$gvar
  )
  stag_data$y <- alpha_i[stag_data$id] +
    0.1 * stag_data$time +
    2.0 * treated_post +
    rnorm(nrow(stag_data), sd = 0.5)

  result <- lwdid::sensitivity_no_anticipation(
    data = stag_data, y = "y",
    ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    verbose = FALSE
  )
  expect_s3_class(result, "lwdid_sensitivity")
  expect_equal(result$type, "no_anticipation")
  expect_false(result$anticipation_detected)
})

# TC-18: S3 methods (print/summary) work
test_that("TC-18: S3 print and summary produce output", {
  set.seed(42)
  n_units <- 20
  n_periods <- 8
  test_data <- expand.grid(
    id = 1:n_units, time = 1:n_periods
  )
  test_data$treat <- as.integer(test_data$id <= 10)
  test_data$post <- as.integer(test_data$time >= 5)
  test_data$y <- rnorm(nrow(test_data)) +
    2.0 * test_data$treat * test_data$post

  result <- lwdid::sensitivity_no_anticipation(
    data = test_data, y = "y",
    ivar = "id", tvar = "time",
    d = "treat", post = "post",
    rolling = "demean", verbose = FALSE
  )
  expect_output(print(result))
  expect_output(summary(result))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- plot(result)
    expect_s3_class(p, "ggplot")
  }
})


# ============================================================================
# Group 5: Numerical reasonableness tests (TC-19 to TC-22)
# ============================================================================

test_that("TC-19~22: numerical reasonableness checks", {
  set.seed(42)
  n_units <- 40
  n_periods <- 10
  treat_time <- 6
  n_treated <- 20
  alpha_i <- rnorm(n_units)
  delta_t <- seq(0, 0.9, length.out = n_periods)

  test_data <- expand.grid(id = 1:n_units, time = 1:n_periods)
  test_data$treat <- as.integer(test_data$id <= n_treated)
  test_data$post <- as.integer(test_data$time >= treat_time)
  test_data$y <- alpha_i[test_data$id] +
    delta_t[test_data$time] +
    2.0 * test_data$treat * test_data$post +
    rnorm(nrow(test_data), sd = 0.5)

  result <- lwdid::sensitivity_no_anticipation(
    data = test_data, y = "y",
    ivar = "id", tvar = "time",
    d = "treat", post = "post",
    rolling = "demean",
    max_anticipation = 3L, verbose = FALSE
  )

  converged <- Filter(
    function(e) isTRUE(e$converged), result$estimates
  )
  expect_true(length(converged) >= 1L)

  for (est in converged) {
    k <- est$excluded_periods
    # TC-19: ATT is finite
    expect_true(is.finite(est$att),
      info = sprintf("exclude=%d: ATT not finite", k))
    # TC-20: SE > 0
    expect_true(est$se > 0,
      info = sprintf("exclude=%d: SE not positive", k))
    # TC-21: CI contains ATT
    expect_true(est$ci_lower <= est$att,
      info = sprintf("exclude=%d: ci_lower > att", k))
    expect_true(est$ci_upper >= est$att,
      info = sprintf("exclude=%d: ci_upper < att", k))
    # TC-22: p-value in [0, 1]
    expect_true(est$pvalue >= 0,
      info = sprintf("exclude=%d: pvalue < 0", k))
    expect_true(est$pvalue <= 1,
      info = sprintf("exclude=%d: pvalue > 1", k))
  }
})
