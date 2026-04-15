# ============================================================================
# test-estimate-ra.R — E2-06 Layers 2-6: RA Estimator & Integration Tests
#
# Test Groups:
#   Layer 2 (T2-01 to T2-07): RA estimator correctness
#   Layer 3 (T3-01, T3-02): Stata numerical consistency
#   Layer 4 (T4-01, T4-04 to T4-06): Estimator edge cases
#   Layer 5 (T5-01 to T5-07): Three-tier fallback strategy
#   Layer 6 (T6-01 to T6-08): Period-specific effects
# ============================================================================

# Helper: capture warnings while preserving result
capture_with_warnings <- function(expr) {
  warnings_caught <- list()
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings_caught)
}

# Helper: check if any warning has a given class
has_warning_class <- function(warnings_list, cls) {
  any(vapply(warnings_list, function(w) inherits(w, cls), logical(1)))
}

# ============================================================================
# Layer 2: RA Estimator Correctness Tests
# ============================================================================

# --- T2-01: No controls simple difference ---
test_that("T2-01: no controls ATT equals simple difference of means", {
  y_trans <- c(5, 3, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  est <- estimate_ra_common(y_trans, d)

  expect_equal(est$att, 2.5, tolerance = 1e-12)
  expect_equal(est$df, 2L)
  expect_equal(est$controls_tier, "none")
  expect_equal(est$n, 4L)
  expect_equal(est$n_treated, 2L)
  expect_equal(est$n_control, 2L)
  expect_equal(est$K, 0L)
  expect_true(is.finite(est$att))
  expect_true(est$se > 0)
  expect_equal(est$t_stat, est$att / est$se, tolerance = 1e-12)
})

# --- T2-02: Tier 1 full interaction model ---
test_that("T2-02: Tier 1 with N1=10, N0=10, K=2 uses full_interaction", {
  set.seed(42)
  n1 <- 10L; n0 <- 10L; K <- 2L
  n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(n, sd = 0.5)
  est <- estimate_ra_common(y_trans, d, x)

  expect_equal(est$controls_tier, "full_interaction")
  expect_equal(est$df, 14L)
  expect_equal(est$K, K)
  expect_equal(est$n, n)
  expect_true(abs(est$att - 2) < 1.5)
  expect_equal(ncol(est$X_design), 6L)
})

# --- T2-03: Tier 2 simple controls model ---
test_that("T2-03: Tier 2 with N1=2, N0=10, K=2 uses simple controls", {
  set.seed(42)
  n1 <- 2L; n0 <- 10L; K <- 2L
  n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(n, sd = 0.5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_true(has_warning_class(cw$warnings, "lwdid_data"))
  expect_equal(est$controls_tier, "simple")
  expect_equal(est$df, 8L)
  expect_equal(ncol(est$X_design), 4L)
})

# --- T2-04: BUG-009 centering in Tier 1 ---
test_that("T2-04: Tier 1 centers at treated-group mean (BUG-009)", {
  set.seed(42)
  d <- c(rep(1L, 3L), rep(0L, 3L))
  x <- matrix(c(10, 20, 30, 1, 2, 3), ncol = 1)
  y_trans <- c(15, 25, 35, 5, 6, 7)
  est <- estimate_ra_common(y_trans, d, x)

  expect_equal(est$controls_tier, "full_interaction")
  x_c_treated <- est$X_design[d == 1, 3]
  expect_equal(mean(x_c_treated), 0, tolerance = 1e-12)
  simple_diff <- mean(y_trans[d == 1]) - mean(y_trans[d == 0])
  expect_true(abs(est$att - simple_diff) > 0.01)
})

# --- T2-05: Tier 3 degradation ---
test_that("T2-05: Tier 3 drops controls when N <= K+2", {
  d <- c(1L, 0L, 0L)
  x <- matrix(c(1, 2, 3), ncol = 1)
  y_trans <- c(10, 3, 5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_equal(est$controls_tier, "none")
  expect_equal(est$K, 0L)
  expect_equal(est$df, 1L)
})

# --- T2-06: Rank deficient design matrix ---
test_that("T2-06: constant control variable throws lwdid_singular_design", {
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  x <- matrix(rep(5, 6), ncol = 1)
  y_trans <- c(10, 11, 12, 3, 4, 5)

  expect_error(
    estimate_ra_common(y_trans, d, x),
    class = "lwdid_singular_design"
  )
})

# --- T2-07: Degenerate SE (zero residuals) ---
test_that("T2-07: degenerate SE sets inference to NA with warning", {
  y_trans <- c(5, 5, 2, 2)
  d <- c(1L, 1L, 0L, 0L)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d))
  est <- cw$result

  expect_equal(est$att, 3.0, tolerance = 1e-12)
  expect_true(is.na(est$se))
  expect_true(is.na(est$t_stat))
  expect_true(is.na(est$pvalue))
  expect_true(is.na(est$ci_lower))
  expect_true(is.na(est$ci_upper))
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))
})

# ============================================================================
# Layer 3: Stata Numerical Consistency Tests
# ============================================================================

# --- T3-01: Demean end-to-end with smoking data ---
test_that("T3-01: demean ATT/SE/df/t-stat/N match Stata benchmarks", {
  data(smoking)
  result <- lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
                  d = "d", post = "post", rolling = "demean")

  expect_equal(result$att, -0.4221746, tolerance = 1e-6)
  expect_equal(result$se_att, 0.1207995, tolerance = 1e-4)
  expect_equal(result$df_resid, 37L)
  expect_equal(result$t_stat, -3.4948, tolerance = 1e-3)
  expect_equal(result$nobs, 39L)
  expect_true(result$att < 0)
  expect_true(result$se_att > 0)
  expect_true(result$pvalue < 0.05)
})

# --- T3-02: Detrend end-to-end with smoking data ---
test_that("T3-02: detrend ATT/SE/df/t-stat/N match Stata benchmarks", {
  data(smoking)
  cw <- capture_with_warnings(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "detrend")
  )
  result <- cw$result

  expect_equal(result$att, -0.2269887, tolerance = 1e-6)
  expect_equal(result$se_att, 0.0940689, tolerance = 1e-4)
  expect_equal(result$df_resid, 37L)
  expect_equal(result$t_stat, -2.4130, tolerance = 1e-3)
  expect_equal(result$nobs, 39L)
  expect_true(result$att < 0)
  expect_true(abs(result$att) < abs(-0.4221746))
  expect_true(result$pvalue < 0.05)
})

# ============================================================================
# Layer 4: Estimator Edge Cases
# ============================================================================

# --- T4-01: Extreme small sample N1=1, N0=2 ---
test_that("T4-01: extreme small sample N1=1, N0=2 produces valid inference", {
  y_trans <- c(10, 3, 5)
  d <- c(1L, 0L, 0L)
  est <- estimate_ra_common(y_trans, d)

  expect_equal(est$n, 3L)
  expect_equal(est$n_treated, 1L)
  expect_equal(est$n_control, 2L)
  expect_equal(est$df, 1L)
  expect_equal(est$controls_tier, "none")
  expect_equal(est$att, 6, tolerance = 1e-12)
  expect_true(is.finite(est$se))
  expect_true(est$se > 0)
  expect_true(is.finite(est$t_stat))
})

# --- T4-04: N1=0 throws lwdid_insufficient_data ---
test_that("T4-04: N1=0 throws lwdid_insufficient_data error", {
  expect_error(
    estimate_ra_common(c(1, 2, 3), c(0L, 0L, 0L)),
    class = "lwdid_insufficient_data"
  )
})

# --- T4-05: N0=0 throws lwdid_insufficient_data ---
test_that("T4-05: N0=0 throws lwdid_insufficient_data error", {
  expect_error(
    estimate_ra_common(c(1, 2, 3), c(1L, 1L, 1L)),
    class = "lwdid_insufficient_data"
  )
})

# --- T4-06: N=2 (N1=1, N0=1) throws lwdid_insufficient_data ---
test_that("T4-06: N=2 with df=0 throws lwdid_insufficient_data error", {
  expect_error(
    estimate_ra_common(c(5, 3), c(1L, 0L)),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Layer 5: Three-Tier Fallback Strategy Tests
# ============================================================================

# --- T5-01: Tier 1 full interaction (N1=5, N0=8, K=3) ---
test_that("T5-01: Tier 1 with N1=5, N0=8, K=3 uses full_interaction", {
  set.seed(101)
  n1 <- 5L; n0 <- 8L; K <- 3L; n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(n, sd = 0.3)
  est <- estimate_ra_common(y_trans, d, x)

  expect_equal(est$controls_tier, "full_interaction")
  expect_equal(est$K, K)
  expect_equal(est$df, 5L)
  expect_equal(ncol(est$X_design), 8L)
  expect_equal(est$n, n)
})

# --- T5-02: Tier 2 simple controls (N1=3, N0=15, K=3) ---
test_that("T5-02: Tier 2 with N1=3, N0=15, K=3 uses simple controls", {
  set.seed(102)
  n1 <- 3L; n0 <- 15L; K <- 3L; n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + rnorm(n, sd = 0.5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_true(has_warning_class(cw$warnings, "lwdid_data"))
  expect_equal(est$controls_tier, "simple")
  expect_equal(est$K, K)
  expect_equal(est$df, 13L)
  expect_equal(ncol(est$X_design), 5L)
})

# --- T5-03: Tier 2 boundary on N0 side (N1=10, N0=2, K=2) ---
test_that("T5-03: Tier 2 with N1=10, N0=2, K=2 (N0<=K+1) uses simple", {
  set.seed(103)
  n1 <- 10L; n0 <- 2L; K <- 2L; n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + rnorm(n, sd = 0.5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_equal(est$controls_tier, "simple")
  expect_equal(est$df, 8L)
})

# --- T5-04: Tier 3 drops all controls (N1=1, N0=2, K=2) ---
test_that("T5-04: Tier 3 with N=3, K=2 (N<=K+2) drops all controls", {
  d <- c(1L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  y_trans <- c(10, 3, 5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_equal(est$controls_tier, "none")
  expect_equal(est$K, 0L)
  expect_equal(est$df, 1L)
  expect_equal(ncol(est$X_design), 2L)
})

# --- T5-05: Boundary — N1=K+1 does NOT satisfy Tier 1 ---
test_that("T5-05: N1=K+1=4 fails Tier 1 strict inequality, degrades to Tier 2", {
  set.seed(105)
  n1 <- 4L; n0 <- 10L; K <- 3L; n <- n1 + n0
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  y_trans <- 1 + 2 * d + rnorm(n, sd = 0.5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_equal(est$controls_tier, "simple")
  expect_true(has_warning_class(cw$warnings, "lwdid_data"))
})

# --- T5-06: Boundary — N=K+2 does NOT satisfy Tier 2 ---
test_that("T5-06: N=K+2=4 fails Tier 2 strict inequality, degrades to Tier 3", {
  d <- c(1L, 1L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2)
  y_trans <- c(10, 12, 3, 5)

  cw <- capture_with_warnings(estimate_ra_common(y_trans, d, x))
  est <- cw$result

  expect_equal(est$controls_tier, "none")
  expect_equal(est$K, 0L)
  expect_equal(est$df, 2L)
})

# --- T5-07: No controls does not trigger degradation ---
test_that("T5-07: no controls does not trigger any degradation warning", {
  y_trans <- c(5, 3, 1, 2)
  d <- c(1L, 1L, 0L, 0L)

  expect_silent(estimate_ra_common(y_trans, d))
  est <- estimate_ra_common(y_trans, d)
  expect_equal(est$controls_tier, "none")
  expect_equal(est$df, 2L)
})

# ============================================================================
# Layer 6: Period-Specific Effects Tests
# ============================================================================

# --- T6-01: tvar parameter correctly passed ---
test_that("T6-01: period effects use correct tvar and return expected periods", {
  set.seed(601)
  dt <- data.table::data.table(
    id   = rep(1:3, each = 5),
    year = rep(2001:2005, 3),
    d    = rep(c(1L, 0L, 0L), each = 5)
  )
  dt[, y_trans := ifelse(d == 1L, 5 + rnorm(.N, sd = 0.5),
                         2 + rnorm(.N, sd = 0.5))]

  result <- estimate_period_effects(
    dt, y_trans_col = "y_trans", d_col = "d",
    tvar = "year", x = NULL,
    periods = c(2003L, 2004L, 2005L), alpha = 0.05
  )

  expect_equal(nrow(result), 3L)
  expect_equal(result$period, c(2003L, 2004L, 2005L))
  expect_true(all(!is.na(result$att)))
  expect_true(all(result$n_obs == 3L))
  expect_true(all(result$n_treated == 1L))
  expect_true(all(result$n_control == 2L))
})

# --- T6-02: controls_tier column returned correctly ---
test_that("T6-02: period effects return correct controls_tier", {
  set.seed(602)
  dt_no_x <- data.table::data.table(
    id   = rep(1:4, each = 3),
    t    = rep(3:5, 4),
    d    = rep(c(1L, 1L, 0L, 0L), each = 3),
    y_trans = c(5, 6, 7, 4, 5, 6, 1, 2, 3, 2, 3, 4)
  )
  res_no_x <- estimate_period_effects(
    dt_no_x, "y_trans", "d", "t", x = NULL,
    periods = 3:5, alpha = 0.05
  )
  expect_true(all(res_no_x$controls_tier == "none"))

  dt_x <- data.table::data.table(
    id   = rep(1:6, each = 3),
    t    = rep(3:5, 6),
    d    = rep(c(1L, 1L, 1L, 0L, 0L, 0L), each = 3),
    y_trans = rnorm(18, mean = 3),
    x1   = rnorm(18)
  )
  res_x <- estimate_period_effects(
    dt_x, "y_trans", "d", "t", x = "x1",
    periods = 3:5, alpha = 0.05
  )
  expect_true(all(res_x$controls_tier == "full_interaction"))
})

# --- T6-03: Complete 13-column set ---
test_that("T6-03: period effects return all 13 required columns", {
  dt <- data.table::data.table(
    id   = rep(1:4, each = 3),
    t    = rep(3:5, 4),
    d    = rep(c(1L, 1L, 0L, 0L), each = 3),
    y_trans = c(5, 6, 7, 4, 5, 6, 1, 2, 3, 2, 3, 4)
  )
  result <- estimate_period_effects(
    dt, "y_trans", "d", "t", x = NULL,
    periods = 3:5, alpha = 0.05
  )

  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated", "n_control",
                     "df", "vce_type", "n_clusters", "controls_tier")
  expect_equal(names(result), expected_cols)
  expect_equal(ncol(result), 15L)
})

# --- T6-04: NA units filtered in period estimation ---
test_that("T6-04: units with NA y_trans are filtered per period", {
  dt <- data.table::data.table(
    id   = rep(1:4, each = 3),
    t    = rep(3:5, 4),
    d    = rep(c(1L, 0L, 0L, 0L), each = 3),
    y_trans = c(5, 6, 7, 1, 2, 3, 2, 3, 4, NA, NA, NA)
  )
  result <- estimate_period_effects(
    dt, "y_trans", "d", "t", x = NULL,
    periods = 3:5, alpha = 0.05
  )

  expect_true(all(result$n_obs == 3L))
  expect_true(all(result$n_treated == 1L))
  expect_true(all(result$n_control == 2L))
  expect_true(all(!is.na(result$att)))
})

# --- T6-05: Unbalanced panel per-period degradation ---
test_that("T6-05: per-period tier degradation with small N1", {
  set.seed(605)
  n1 <- 2L; n0 <- 8L
  dt <- data.table::data.table(
    id   = rep(1:(n1 + n0), each = 3),
    t    = rep(3:5, n1 + n0),
    d    = rep(c(rep(1L, n1), rep(0L, n0)), each = 3),
    y_trans = rnorm((n1 + n0) * 3, mean = 3),
    x1   = rnorm((n1 + n0) * 3),
    x2   = rnorm((n1 + n0) * 3)
  )

  cw <- capture_with_warnings(
    estimate_period_effects(
      dt, "y_trans", "d", "t", x = c("x1", "x2"),
      periods = 3:5, alpha = 0.05
    )
  )
  result <- cw$result
  expect_true(all(result$controls_tier == "simple"))
})

# --- T6-06: tryCatch fallback when period has N1=0 ---
test_that("T6-06: period with N1=0 returns NA row with warning", {
  # Period 3: N1=1 (normal), Period 4: N1=0 (treated unit missing)
  dt <- data.table::data.table(
    id   = c(1, 1, 2, 2, 3, 3),
    t    = c(3, 4, 3, 4, 3, 4),
    d    = c(1L, NA, 0L, 0L, 0L, 0L),
    y_trans = c(5, NA, 1, 2, 3, 4)
  )
  # Remove the treated unit from period 4 entirely
  dt2 <- data.table::data.table(
    id   = c(1, 2, 2, 3, 3),
    t    = c(3, 3, 4, 3, 4),
    d    = c(1L, 0L, 0L, 0L, 0L),
    y_trans = c(5, 1, 2, 3, 4)
  )

  cw <- capture_with_warnings(
    estimate_period_effects(
      dt2, "y_trans", "d", "t", x = NULL,
      periods = c(3L, 4L), alpha = 0.05
    )
  )
  result <- cw$result

  # Period 3: normal (N1=1, N0=2)
  expect_true(!is.na(result$att[1]))
  expect_equal(result$n_obs[1], 3L)

  # Period 4: N1=0 -> error caught, NA row
  expect_true(is.na(result$att[2]))
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))
})

# --- T6-07: Degenerate SE in period estimation ---
test_that("T6-07: degenerate period SE sets all periods to NA", {
  # Treated all 5, control all 2 -> perfect fit, se=0
  dt <- data.table::data.table(
    id   = rep(1:4, each = 3),
    t    = rep(3:5, 4),
    d    = rep(c(1L, 1L, 0L, 0L), each = 3),
    y_trans = rep(c(5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2))
  )

  cw <- capture_with_warnings(
    estimate_period_effects(
      dt, "y_trans", "d", "t", x = NULL,
      periods = 3:5, alpha = 0.05
    )
  )
  result <- cw$result

  # All periods should have degenerate results (ATT=NA due to se=0)
  expect_true(all(is.na(result$att)))
  expect_true(all(is.na(result$se)))
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))
})

# --- T6-08: Period ATT mean equals summary ATT (balanced, no controls) ---
test_that("T6-08: mean of period ATTs equals summary ATT for balanced panel", {
  set.seed(608)
  # 6 units x 5 periods (2 treated + 4 control), g=3
  n_units <- 6L
  n_periods <- 5L
  dt <- data.table::data.table(
    id = rep(1:n_units, each = n_periods),
    t  = rep(1:n_periods, n_units),
    d  = rep(c(1L, 1L, 0L, 0L, 0L, 0L), each = n_periods)
  )
  dt[, y := ifelse(d == 1L & t >= 3, 5 + rnorm(.N, sd = 0.3),
                   2 + rnorm(.N, sd = 0.3))]

  # Transform (demean)
  transformed <- transform_demean(dt, "y", "id", "t", g = 3L)

  # Period effects (post periods: t=3,4,5)
  period_result <- estimate_period_effects(
    transformed, "y_trans", "d", "t", x = NULL,
    periods = 3:5, alpha = 0.05
  )

  # Summary ATT: mean of post-period y_trans per unit, then RA
  summary_dt <- transformed[t >= 3, .(y_trans_summary = mean(y_trans)), by = "id"]
  summary_dt <- merge(summary_dt, unique(dt[, c("id", "d")]), by = "id")
  summary_est <- estimate_ra_common(summary_dt$y_trans_summary, summary_dt$d)

  # Mean of period ATTs should equal summary ATT
  mean_period_att <- mean(period_result$att)
  expect_equal(mean_period_att, summary_est$att, tolerance = 1e-12)
})
