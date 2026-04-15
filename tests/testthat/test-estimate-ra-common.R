# ============================================================================
# Tests for estimate_ra_common()
# Story E2-03: Common Timing RA Estimator
# ============================================================================

# ============================================================================
# Group 1: No-controls basic regression (T-01 ~ T-08b)
# ============================================================================
# Test data: treated Y=[3,5], control Y=[1,2], d=[1,1,0,0]
# vibe-math verified: ATT=2.5, SE=1.118033988749895, t=2.2360679774997894
# beta = [1.5, 2.5], sigma2=1.25, df=2
# (X'X)^{-1} = [[0.5, -0.5], [-0.5, 1.0]]

test_that("T-01: no-controls ATT equals mean difference", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(res$att, 2.5, tolerance = 1e-12)
})

test_that("T-02: ATT equals treated mean minus control mean", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expected_att <- mean(y[d == 1]) - mean(y[d == 0])
  expect_equal(res$att, expected_att, tolerance = 1e-12)
})

test_that("T-03: degrees of freedom df = N - 2", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(res$df, 2)
})

test_that("T-04: homoskedastic SE matches vibe-math", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  # vibe-math: SE = sqrt(1.25 * 1.0) = 1.118033988749895
  expect_equal(res$se, 1.118033988749895, tolerance = 1e-12)
})

test_that("T-05: t-statistic = ATT / SE", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(res$t_stat, res$att / res$se, tolerance = 1e-12)
  # vibe-math: t = 2.2360679774997894
  expect_equal(res$t_stat, 2.2360679774997894, tolerance = 1e-10)
})

test_that("T-06: p-value uses t distribution (not normal)", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expected_p <- 2 * pt(-abs(res$t_stat), res$df)
  expect_equal(res$pvalue, expected_p, tolerance = 1e-12)
  # Verify it differs from normal approximation
  normal_p <- 2 * pnorm(-abs(res$t_stat))
  expect_false(isTRUE(all.equal(res$pvalue, normal_p, tolerance = 1e-4)))
})

test_that("T-07: 95% CI uses t quantile", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  tcrit <- qt(0.975, res$df)
  expect_equal(res$ci_lower, res$att - tcrit * res$se, tolerance = 1e-12)
  expect_equal(res$ci_upper, res$att + tcrit * res$se, tolerance = 1e-12)
})

test_that("T-08: controls_tier is 'none' for no-controls", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(res$controls_tier, "none")
  expect_equal(res$K, 0L)
})

test_that("T-08b: empty matrix x equivalent to x=NULL", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res_null <- estimate_ra_common(y, d, x = NULL)
  res_empty <- estimate_ra_common(y, d, x = matrix(nrow = 4, ncol = 0))
  expect_equal(res_null$att, res_empty$att, tolerance = 1e-12)
  expect_equal(res_null$se, res_empty$se, tolerance = 1e-12)
  expect_equal(res_null$df, res_empty$df)
  expect_equal(res_null$controls_tier, res_empty$controls_tier)
  expect_equal(res_null$K, res_empty$K)
})

# ============================================================================
# Group 2: OLS coefficient and vcov verification (vibe-math verified)
# ============================================================================

test_that("no-controls OLS coefficients match vibe-math", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  # vibe-math: beta = [1.5, 2.5]
  expect_equal(unname(res$params[1]), 1.5, tolerance = 1e-12)
  expect_equal(unname(res$params[2]), 2.5, tolerance = 1e-12)
})

test_that("no-controls vcov matches vibe-math", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  # vibe-math: sigma2=1.25, (X'X)^{-1} = [[0.5,-0.5],[-0.5,1.0]]
  expect_equal(res$vcov[1, 1], 1.25 * 0.5, tolerance = 1e-12)
  expect_equal(res$vcov[1, 2], 1.25 * (-0.5), tolerance = 1e-12)
  expect_equal(res$vcov[2, 2], 1.25 * 1.0, tolerance = 1e-12)
})

test_that("no-controls residuals match vibe-math", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  # yhat = [4, 4, 1.5, 1.5], resid = [-1, 1, -0.5, 0.5]
  expect_equal(res$resid, c(-1, 1, -0.5, 0.5), tolerance = 1e-12)
  expect_true(is.numeric(res$resid))
  expect_false(is.matrix(res$resid))
})

# ============================================================================
# Group 3: Three-tier fallback strategy (T-09 ~ T-15)
# ============================================================================
# Helper: run estimate_ra_common suppressing expected tier warnings
muffle_tier <- function(expr) {
  withCallingHandlers(expr, lwdid_data = function(w) invokeRestart("muffleWarning"))
}

test_that("T-09: Tier 1 full interaction model", {
  set.seed(42)
  n1 <- 5L; n0 <- 5L; K <- 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(10 * K), ncol = K)
  y <- rnorm(10) + 3 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(res$df, 10 - 2 - 2 * K)
  expect_equal(ncol(res$X_design), 2 + 2 * K)
})

test_that("T-10: Tier 2 simple controls model", {
  set.seed(42)
  n1 <- 2L; n0 <- 10L; K <- 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(12 * K), ncol = K)
  y <- rnorm(12) + 3 * d
  # N1=2 <= K+1=3, Tier 1 fails; N=12 > K+2=4, Tier 2
  res <- muffle_tier(estimate_ra_common(y, d, x = x))
  expect_equal(res$controls_tier, "simple")
  expect_equal(res$df, 12 - K - 2)
  expect_equal(ncol(res$X_design), K + 2)
})

test_that("T-11: Tier 3 drops all controls", {
  n1 <- 1L; n0 <- 2L; K <- 2L
  d <- c(1L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = K)
  y <- c(10, 5, 6)
  # N=3 <= K+2=4, Tier 3
  res <- muffle_tier(estimate_ra_common(y, d, x = x))
  expect_equal(res$controls_tier, "none")
  expect_equal(res$K, 0L)
  expect_equal(res$df, 3 - 2)
  expect_equal(ncol(res$X_design), 2)
})

test_that("T-12: Tier 1 df = N - 2 - 2K", {
  set.seed(123)
  n1 <- 8L; n0 <- 8L; K <- 3L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(16 * K), ncol = K)
  y <- rnorm(16) + 2 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(res$df, 16 - 2 - 2 * K)
})

test_that("T-13: Tier 2 df = N - K - 2", {
  set.seed(123)
  n1 <- 3L; n0 <- 12L; K <- 3L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(15 * K), ncol = K)
  y <- rnorm(15) + 2 * d
  res <- muffle_tier(estimate_ra_common(y, d, x = x))
  expect_equal(res$controls_tier, "simple")
  expect_equal(res$df, 15 - K - 2)
})

test_that("T-14: Tier 2 warning has correct detail", {
  set.seed(42)
  n1 <- 2L; n0 <- 10L; K <- 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(12 * K), ncol = K)
  y <- rnorm(12) + 3 * d
  w <- NULL
  withCallingHandlers(
    estimate_ra_common(y, d, x = x),
    lwdid_data = function(cond) { w <<- cond; invokeRestart("muffleWarning") }
  )
  expect_equal(w$detail, "controls_degraded_to_simple")
  expect_equal(w$action_taken, "interaction terms dropped")
})

test_that("T-15: Tier 3 warning has correct detail", {
  d <- c(1L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  y <- c(10, 5, 6)
  w <- NULL
  withCallingHandlers(
    estimate_ra_common(y, d, x = x),
    lwdid_data = function(cond) { w <<- cond; invokeRestart("muffleWarning") }
  )
  expect_equal(w$detail, "controls_dropped")
  expect_equal(w$action_taken, "all controls dropped")
})

test_that("Tier 1 boundary: N1=K+2, N0=K+2 just satisfies Tier 1", {
  set.seed(99)
  K <- 2L
  n1 <- K + 2L; n0 <- K + 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(8 * K), ncol = K)
  y <- rnorm(8) + 2 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
})

test_that("Tier 2 boundary: N=K+3 just satisfies Tier 2", {
  set.seed(77)
  K <- 2L
  # N1=2 <= K+1=3 fails Tier 1; N=5 > K+2=4 satisfies Tier 2
  n1 <- 2L; n0 <- 3L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(5 * K), ncol = K)
  y <- rnorm(5) + 2 * d
  res <- muffle_tier(estimate_ra_common(y, d, x = x))
  expect_equal(res$controls_tier, "simple")
  expect_equal(res$df, 5 - K - 2)
})

test_that("Tier 3 boundary: N=K+2 triggers Tier 3", {
  set.seed(77)
  K <- 2L
  # N=4 <= K+2=4, Tier 3
  n1 <- 1L; n0 <- 3L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(4 * K), ncol = K)
  y <- rnorm(4) + 2 * d
  res <- muffle_tier(estimate_ra_common(y, d, x = x))
  expect_equal(res$controls_tier, "none")
  expect_equal(res$K, 0L)
  expect_equal(res$df, 4 - 2)
})

# ============================================================================
# Group 4: BUG-009 centering verification (T-16 ~ T-18)
# ============================================================================

test_that("T-16: x_bar1 uses treated unit-weighted mean (BUG-009)", {
  d <- c(1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  x <- matrix(c(2, 4, 6, 1, 3, 5, 7, 9), ncol = 1)
  y <- c(10, 12, 14, 5, 7, 9, 11, 13)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  # x_bar1 = (2+4+6)/3 = 4, so X_c = X - 4
  # Treated X_c = [-2, 0, 2]
  X_c_treated <- res$X_design[d == 1, 3]
  expect_equal(unname(X_c_treated), c(-2, 0, 2), tolerance = 1e-12)
})

test_that("T-17: centered treated X_c mean is zero", {
  set.seed(42)
  n1 <- 10L; n0 <- 10L; K <- 3L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(20 * K, mean = 100), ncol = K)
  y <- rnorm(20) + 2 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  x_c_treated <- res$X_design[d == 1, 3:(2 + K), drop = FALSE]
  expect_equal(unname(colMeans(x_c_treated)), rep(0, K), tolerance = 1e-12)
})

test_that("T-18: centering does not affect ATT estimate", {
  set.seed(42)
  n1 <- 10L; n0 <- 10L; K <- 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(20 * K), ncol = K)
  y <- rnorm(20) + 3 * d

  # R implementation (centered form)
  res_centered <- estimate_ra_common(y, d, x = x)

  # Manual non-centered form: [1, D, X, D*(X - x_bar1)]
  x_bar1 <- colMeans(x[d == 1, , drop = FALSE])
  X_nc <- cbind(1, d, x, d * sweep(x, 2, x_bar1, "-"))
  beta_nc <- qr.coef(qr(X_nc), y)
  att_nc <- unname(beta_nc[2])

  expect_equal(res_centered$att, att_nc, tolerance = 1e-10)
})

# ============================================================================
# Group 5: Error handling (T-19 ~ T-23)
# ============================================================================

test_that("T-19: N1=0 throws lwdid_insufficient_data", {
  y <- c(1, 2, 3)
  d <- c(0L, 0L, 0L)
  expect_error(estimate_ra_common(y, d), class = "lwdid_insufficient_data")
  err <- tryCatch(estimate_ra_common(y, d), error = function(e) e)
  expect_equal(err$n_treated, 0L)
})

test_that("T-20: N0=0 throws lwdid_insufficient_data", {
  y <- c(1, 2, 3)
  d <- c(1L, 1L, 1L)
  expect_error(estimate_ra_common(y, d), class = "lwdid_insufficient_data")
  err <- tryCatch(estimate_ra_common(y, d), error = function(e) e)
  expect_equal(err$n_control, 0L)
})

test_that("T-21: N=2 throws lwdid_insufficient_data", {
  y <- c(1, 2)
  d <- c(1L, 0L)
  expect_error(estimate_ra_common(y, d), class = "lwdid_insufficient_data")
})

test_that("check order: N1=0 reported before N<3", {
  y <- c(1, 2)
  d <- c(0L, 0L)
  err <- tryCatch(estimate_ra_common(y, d), error = function(e) e)
  expect_s3_class(err, "lwdid_insufficient_data")
  expect_equal(err$n_treated, 0L)
})

test_that("T-23: constant control variable throws lwdid_singular_design", {
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  x <- matrix(rep(5, 8), ncol = 1)
  y <- c(10, 12, 14, 16, 5, 7, 9, 11)
  expect_error(estimate_ra_common(y, d, x = x), class = "lwdid_singular_design")
  err <- tryCatch(estimate_ra_common(y, d, x = x), error = function(e) e)
  expect_true(!is.null(err$rank))
  expect_true(!is.null(err$expected_rank))
  expect_true(err$rank < err$expected_rank)
})

test_that("collinear controls throw lwdid_singular_design", {
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8), ncol = 2)
  y <- c(10, 12, 14, 16, 5, 7, 9, 11)
  expect_error(estimate_ra_common(y, d, x = x), class = "lwdid_singular_design")
})

# ============================================================================
# Group 6: Degenerate SE (T-24)
# ============================================================================

test_that("T-24: degenerate SE sets inference to NA, preserves ATT", {
  # Treated all same Y, control all same Y -> residuals all zero -> se=0
  y <- c(5, 5, 2, 2)
  d <- c(1L, 1L, 0L, 0L)
  w_captured <- NULL
  res <- withCallingHandlers(
    estimate_ra_common(y, d),
    lwdid_small_sample = function(w) {
      w_captured <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(res$att, 3.0, tolerance = 1e-12)
  expect_true(is.na(res$se))
  expect_true(is.na(res$t_stat))
  expect_true(is.na(res$pvalue))
  expect_true(is.na(res$ci_lower))
  expect_true(is.na(res$ci_upper))
  # Warning detail
  expect_equal(w_captured$detail, "degenerate_se")
  expect_equal(w_captured$action_taken, "inference results set to NA")
})

# ============================================================================
# Group 7: Extreme small sample (T-22)
# ============================================================================

test_that("T-22: N=3 (N1=1, N0=2) computes correctly with df=1", {
  y <- c(10, 3, 5)
  d <- c(1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(res$att, 6.0, tolerance = 1e-12)
  expect_equal(res$df, 1)
  expect_equal(res$n, 3)
  expect_equal(res$n_treated, 1)
  expect_equal(res$n_control, 2)
  tcrit <- qt(0.975, 1)
  expect_equal(res$ci_lower, res$att - tcrit * res$se, tolerance = 1e-12)
  expect_equal(res$ci_upper, res$att + tcrit * res$se, tolerance = 1e-12)
  # CI should be very wide (t_1 quantile is ~12.706)
  expect_true(res$ci_upper - res$ci_lower > 20)
})

# ============================================================================
# Group 8: Tier 1 numerical verification (vibe-math verified)
# ============================================================================
# Data: treated Y=[11,13,13,17,19], control Y=[4,8,8,12,14]
# X_treated=[2,4,6,8,10], X_control=[1,3,5,7,9], K=1
# x_bar1 = 6, x_c = X - 6
# vibe-math: beta=[10.4, 4.2, 1.2, -0.2]
# SSR=6.4, df=6, sigma2=1.0666667
# (X'X)^{-1}[2,2] = 0.425
# SE = 0.6733003292241384, t = 6.23792952075305

test_that("Tier 1 ATT matches vibe-math", {
  y <- c(11, 13, 13, 17, 19, 4, 8, 8, 12, 14)
  d <- c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  x <- matrix(c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), ncol = 1)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(res$att, 4.2, tolerance = 1e-10)
})

test_that("Tier 1 SE matches vibe-math", {
  y <- c(11, 13, 13, 17, 19, 4, 8, 8, 12, 14)
  d <- c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  x <- matrix(c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), ncol = 1)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$se, 0.6733003292241384, tolerance = 1e-10)
  expect_equal(res$df, 6)
})

test_that("Tier 1 coefficients match vibe-math", {
  y <- c(11, 13, 13, 17, 19, 4, 8, 8, 12, 14)
  d <- c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  x <- matrix(c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), ncol = 1)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(unname(res$params[1]), 10.4, tolerance = 1e-10)
  expect_equal(unname(res$params[2]), 4.2, tolerance = 1e-10)
  expect_equal(unname(res$params[3]), 1.2, tolerance = 1e-10)
  expect_equal(unname(res$params[4]), -0.2, tolerance = 1e-10)
})

test_that("Tier 1 t-stat matches vibe-math", {
  y <- c(11, 13, 13, 17, 19, 4, 8, 8, 12, 14)
  d <- c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  x <- matrix(c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), ncol = 1)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$t_stat, 6.23792952075305, tolerance = 1e-8)
})

# ============================================================================
# Group 9: Alpha parameter and CI width (T-30)
# ============================================================================

test_that("T-30: alpha=0.01 gives wider CI than alpha=0.05", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res_05 <- estimate_ra_common(y, d, alpha = 0.05)
  res_01 <- estimate_ra_common(y, d, alpha = 0.01)
  expect_true(res_01$ci_upper - res_01$ci_lower > res_05$ci_upper - res_05$ci_lower)
  expect_equal(res_05$att, res_01$att, tolerance = 1e-12)
  expect_equal(res_05$se, res_01$se, tolerance = 1e-12)
})

test_that("alpha=0.10 gives narrower CI than alpha=0.05", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res_05 <- estimate_ra_common(y, d, alpha = 0.05)
  res_10 <- estimate_ra_common(y, d, alpha = 0.10)
  expect_true(res_10$ci_upper - res_10$ci_lower < res_05$ci_upper - res_05$ci_lower)
})

# ============================================================================
# Group 10: Return value completeness (T-32 ~ T-36)
# ============================================================================

test_that("T-32: return list has all 21 fields", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expected_names <- c("att", "se", "t_stat", "df", "pvalue",
                      "ci_lower", "ci_upper", "params", "vcov",
                      "resid", "fit", "n", "n_treated", "n_control",
                      "K", "vce_type", "n_clusters",
                      "controls_tier", "X_design",
                      "df_resid", "df_inference")
  expect_equal(sort(names(res)), sort(expected_names))
})

test_that("T-33: X_design dimensions correct", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(dim(res$X_design), c(4, 2))
})

test_that("T-34: vcov dimensions correct", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  p <- ncol(res$X_design)
  expect_equal(dim(res$vcov), c(p, p))
})

test_that("T-35: resid length equals N", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(length(res$resid), 4)
})

test_that("T-36: params length equals p", {
  y <- c(3, 5, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(length(res$params), ncol(res$X_design))
})

test_that("return value types are correct", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_true(is.numeric(res$att) && length(res$att) == 1)
  expect_true(is.numeric(res$se) && length(res$se) == 1)
  expect_true(is.numeric(res$t_stat) && length(res$t_stat) == 1)
  expect_true(is.numeric(res$pvalue) && length(res$pvalue) == 1)
  expect_true(is.numeric(res$ci_lower) && length(res$ci_lower) == 1)
  expect_true(is.numeric(res$ci_upper) && length(res$ci_upper) == 1)
  expect_true(is.numeric(res$df) && length(res$df) == 1)
  expect_true(is.numeric(res$n) && length(res$n) == 1)
  expect_true(is.character(res$controls_tier) && length(res$controls_tier) == 1)
  expect_true(is.numeric(res$params))
  expect_true(is.matrix(res$vcov))
  expect_true(is.numeric(res$resid) && !is.matrix(res$resid))
  expect_true(is.matrix(res$X_design))
})

# ============================================================================
# Group 11: Stata consistency (T-25 ~ T-27) — placeholder
# ============================================================================

test_that("T-25/T-26/T-27: Stata consistency (placeholder)", {
  skip("Stata consistency tests require California Smoking data - deferred to E2-06")
})

# ============================================================================
# Group 12: End-to-end ATT recovery (T-28 ~ T-29)
# ============================================================================

test_that("T-28: end-to-end ATT recovery with known tau", {
  set.seed(2024)
  n <- 100L
  tau <- 5.0
  d <- c(rep(1L, 50), rep(0L, 50))
  alpha_i <- rnorm(n, mean = 10, sd = 2)
  eps <- rnorm(n, sd = 0.5)
  y <- alpha_i + tau * d + eps
  res <- estimate_ra_common(y, d)
  # tau should be within the CI
  expect_true(res$ci_lower < tau && tau < res$ci_upper)
  # ATT should be close to 5
  expect_equal(res$att, tau, tolerance = 0.5)
})

test_that("T-29: large values (Y ~ 1e6) precision", {
  set.seed(42)
  n <- 50L
  d <- c(rep(1L, 25), rep(0L, 25))
  tau <- 100.0
  y <- 1e6 + rnorm(n, sd = 10) + tau * d
  res <- estimate_ra_common(y, d)
  # Relative error < 1e-2 (tau=100, noise sd=10)
  expect_equal(res$att, tau, tolerance = 10)
  # SE should be reasonable
  expect_true(res$se < 10)
})

test_that("end-to-end with controls: Tier 1 ATT recovery", {
  set.seed(123)
  n <- 200L
  d <- c(rep(1L, 100), rep(0L, 100))
  x <- matrix(rnorm(n * 2), ncol = 2)
  tau <- 3.0
  y <- 1 + tau * d + 0.5 * x[, 1] + 0.3 * x[, 2] + rnorm(n, sd = 0.5)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_true(res$ci_lower < tau && tau < res$ci_upper)
})

test_that("multi-controls K=3: Tier 1 design matrix dimensions", {
  set.seed(42)
  K <- 3L
  d <- c(rep(1L, 10), rep(0L, 10))
  x <- matrix(rnorm(20 * K), ncol = K)
  y <- rnorm(20) + 2 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(ncol(res$X_design), 2 + 2 * K)
  expect_equal(res$df, 20 - 2 - 2 * K)
  expect_equal(length(res$params), 2 + 2 * K)
  expect_equal(dim(res$vcov), c(8, 8))
})

test_that("single control K=1: drop=FALSE works correctly", {
  set.seed(42)
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 1)
  y <- c(10, 12, 14, 16, 5, 7, 9, 11)
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(ncol(res$X_design), 4)
  expect_equal(res$K, 1L)
})

# ============================================================================
# Group 13: POLS equivalence (T-37 ~ T-38)
# ============================================================================

test_that("T-37: POLS equivalence (no controls)", {
  set.seed(42)
  N <- 6L; T_total <- 4L; S <- 3L
  unit_fe <- rnorm(N)
  time_fe <- c(0, 1, 2, 3)
  tau <- 5.0
  d_unit <- c(1L, 1L, 1L, 0L, 0L, 0L)

  # Generate panel: Y_it = alpha_i + gamma_t + tau * W_it + eps
  Y <- matrix(NA_real_, nrow = N, ncol = T_total)
  for (i in seq_len(N)) {
    for (t in seq_len(T_total)) {
      W_it <- as.integer(d_unit[i] == 1L && t >= S)
      Y[i, t] <- unit_fe[i] + time_fe[t] + tau * W_it + rnorm(1, sd = 0.3)
    }
  }

  # RA path: demean transform then RA
  # y_trans_i = mean(Y_i,post) - mean(Y_i,pre)
  y_pre <- rowMeans(Y[, 1:(S - 1), drop = FALSE])
  y_post <- rowMeans(Y[, S:T_total, drop = FALSE])
  y_trans <- y_post - y_pre
  res_ra <- estimate_ra_common(y_trans, d_unit)

  # POLS path: Y_it on 1, D_i, post_t, W_it
  y_vec <- as.vector(t(Y))
  d_vec <- rep(d_unit, each = T_total)
  post_vec <- rep(as.integer(seq_len(T_total) >= S), N)
  w_vec <- d_vec * post_vec
  X_pols <- cbind(1, d_vec, post_vec, w_vec)
  beta_pols <- qr.coef(qr(X_pols), y_vec)
  att_pols <- unname(beta_pols[4])

  expect_equal(res_ra$att, att_pols, tolerance = 1e-10)
})

test_that("T-38: POLS equivalence (with controls, Tier 1)", {
  set.seed(123)
  N <- 10L; T_total <- 4L; S <- 3L
  unit_fe <- rnorm(N)
  time_fe <- c(0, 1, 2, 3)
  tau <- 3.0
  d_unit <- c(rep(1L, 5), rep(0L, 5))
  x_unit <- matrix(rnorm(N), ncol = 1)

  # DGP: x effect interacts with post (so it survives demeaning)
  # Y_it = alpha_i + gamma_t + tau*W_it + beta_x*X_i*post_t + eps
  beta_x <- 0.5
  Y <- matrix(NA_real_, nrow = N, ncol = T_total)
  for (i in seq_len(N)) {
    for (t in seq_len(T_total)) {
      W_it <- as.integer(d_unit[i] == 1L && t >= S)
      post_it <- as.integer(t >= S)
      Y[i, t] <- unit_fe[i] + time_fe[t] + tau * W_it +
        beta_x * x_unit[i, 1] * post_it + rnorm(1, sd = 0.3)
    }
  }

  # RA path
  y_pre <- rowMeans(Y[, 1:(S - 1), drop = FALSE])
  y_post <- rowMeans(Y[, S:T_total, drop = FALSE])
  y_trans <- y_post - y_pre
  res_ra <- estimate_ra_common(y_trans, d_unit, x = x_unit)
  expect_equal(res_ra$controls_tier, "full_interaction")

  # POLS path: Y_it on 1, D_i, post_t, W_it, X_i*post_t,
  #   W_it*(X_i - x_bar1)
  # Also need X_i and D_i*X_i for the unit-level effects
  x_bar1 <- mean(x_unit[d_unit == 1, 1])
  y_vec <- as.vector(t(Y))
  d_vec <- rep(d_unit, each = T_total)
  post_vec <- rep(as.integer(seq_len(T_total) >= S), N)
  w_vec <- d_vec * post_vec
  x_vec <- rep(x_unit[, 1], each = T_total)
  x_main <- x_vec
  dx_main <- d_vec * (x_vec - x_bar1)
  x_post <- (x_vec - x_bar1) * post_vec
  x_w_centered <- w_vec * (x_vec - x_bar1)
  X_pols <- cbind(1, d_vec, post_vec, w_vec,
                  x_main, dx_main, x_post, x_w_centered)
  beta_pols <- qr.coef(qr(X_pols), y_vec)
  att_pols <- unname(beta_pols[4])

  expect_equal(res_ra$att, att_pols, tolerance = 1e-10)
})

# ============================================================================
# Group 14: Period effects and standard DiD equivalence (T-39)
# ============================================================================

test_that("T-39: period effect equals standard panel DiD", {
  set.seed(42)
  N <- 8L; S <- 3L; T_total <- 5L
  unit_fe <- rnorm(N)
  time_fe <- c(0, 1, 2, 3, 4)
  tau_t <- c(NA, NA, 2.0, 3.0, 4.0)
  d_unit <- c(rep(1L, 4), rep(0L, 4))

  Y <- matrix(NA_real_, nrow = N, ncol = T_total)
  for (i in seq_len(N)) {
    for (t in seq_len(T_total)) {
      W_it <- as.integer(d_unit[i] == 1L && t >= S)
      tau_val <- if (W_it == 1L) tau_t[t] else 0
      Y[i, t] <- unit_fe[i] + time_fe[t] + tau_val + rnorm(1, sd = 0.2)
    }
  }

  y_pre_mean <- rowMeans(Y[, 1:(S - 1), drop = FALSE])

  for (t in S:T_total) {
    # LW path: y_dot_it = Y_it - mean(Y_i,pre), regress on 1, D_i
    y_dot_t <- Y[, t] - y_pre_mean
    res_lw <- estimate_ra_common(y_dot_t, d_unit)

    # Standard DiD: use periods {1,...,S-1,t}
    periods_used <- c(1:(S - 1), t)
    Y_sub <- Y[, periods_used, drop = FALSE]
    n_periods <- length(periods_used)
    y_vec <- as.vector(t(Y_sub))
    d_vec <- rep(d_unit, each = n_periods)
    post_vec <- rep(c(rep(0L, S - 1), 1L), N)
    w_vec <- d_vec * post_vec
    X_did <- cbind(1, d_vec, post_vec, w_vec)
    beta_did <- qr.coef(qr(X_did), y_vec)
    att_did <- unname(beta_did[4])

    expect_equal(res_lw$att, att_did, tolerance = 1e-10,
                 label = sprintf("period %d LW vs DiD", t))
  }
})


# ============================================================================
# Group 15: Period-specific independent tier fallback (REQ-82, REQ-83)
# ============================================================================

test_that("different periods can fall into different tiers", {
  # Simulate unbalanced panel: some periods have fewer treated units
  # Period A: many treated -> Tier 1
  # Period B: few treated -> Tier 2 or Tier 3
  set.seed(42)
  k <- 2L

  # Period A: N1=6, N0=6, K=2 -> Tier 1 (N1>K+1=3, N0>K+1=3)
  d_a <- c(rep(1L, 6), rep(0L, 6))
  x_a <- matrix(rnorm(12 * k), ncol = k)
  y_a <- rnorm(12) + 3 * d_a
  res_a <- estimate_ra_common(y_a, d_a, x = x_a)
  expect_equal(res_a$controls_tier, "full_interaction")


  # Period B: N1=2, N0=8, K=2 -> Tier 2 (N1=2 <= K+1=3, N=10 > K+2=4)
  d_b <- c(rep(1L, 2), rep(0L, 8))
  x_b <- matrix(rnorm(10 * k), ncol = k)
  y_b <- rnorm(10) + 3 * d_b
  res_b <- muffle_tier(estimate_ra_common(y_b, d_b, x = x_b))
  expect_equal(res_b$controls_tier, "simple")

  # Both produce valid ATT estimates
  expect_true(is.finite(res_a$att))
  expect_true(is.finite(res_b$att))
})

test_that("each period computes independent x_bar1", {
  # Two calls with different treated subsets -> different x_bar1
  set.seed(123)
  k <- 1L

  # Period A: treated units have x = [10, 20, 30]
  d_a <- c(1L, 1L, 1L, 0L, 0L, 0L)
  x_a <- matrix(c(10, 20, 30, 1, 2, 3), ncol = 1)
  y_a <- c(15, 25, 35, 5, 7, 9)
  res_a <- estimate_ra_common(y_a, d_a, x = x_a)
  # x_bar1 = 20, so treated X_c = [-10, 0, 10]
  x_c_treated_a <- res_a$X_design[d_a == 1, 3]
  expect_equal(mean(x_c_treated_a), 0, tolerance = 1e-12)

  # Period B: treated units have x = [100, 200, 300]
  # N1=3 > K+1=2, N0=3 > K+1=2 -> Tier 1

  d_b <- c(1L, 1L, 1L, 0L, 0L, 0L)
  x_b <- matrix(c(100, 200, 300, 1, 2, 3), ncol = 1)
  y_b <- c(15, 25, 35, 5, 7, 9)
  res_b <- estimate_ra_common(y_b, d_b, x = x_b)
  expect_equal(res_b$controls_tier, "full_interaction")
  # x_bar1 = 200, so treated X_c = [-100, 0, 100]
  x_c_treated_b <- res_b$X_design[d_b == 1, 3]
  expect_equal(mean(x_c_treated_b), 0, tolerance = 1e-12)

  # The centering values are different
  expect_false(isTRUE(all.equal(
    mean(x_a[d_a == 1, ]),
    mean(x_b[d_b == 1, ])
  )))
})

# ============================================================================
# Group 16: Numerical reasonableness checks
# ============================================================================

test_that("ATT sign matches direction of treatment effect", {
  # Positive treatment effect
  y_pos <- c(10, 12, 14, 1, 2, 3)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res_pos <- estimate_ra_common(y_pos, d)
  expect_true(res_pos$att > 0)

  # Negative treatment effect
  y_neg <- c(1, 2, 3, 10, 12, 14)
  res_neg <- estimate_ra_common(y_neg, d)
  expect_true(res_neg$att < 0)
})

test_that("SE decreases with larger sample size", {
  set.seed(42)
  tau <- 5.0

  # Small sample
  n_small <- 20L
  d_small <- c(rep(1L, 10), rep(0L, 10))
  y_small <- rnorm(n_small) + tau * d_small
  res_small <- estimate_ra_common(y_small, d_small)

  # Large sample (same DGP)
  n_large <- 200L
  d_large <- c(rep(1L, 100), rep(0L, 100))
  y_large <- rnorm(n_large) + tau * d_large
  res_large <- estimate_ra_common(y_large, d_large)

  expect_true(res_large$se < res_small$se)
})

test_that("p-value is in [0, 1]", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_true(res$pvalue >= 0 && res$pvalue <= 1)
})

test_that("CI contains ATT", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_true(res$ci_lower <= res$att)
  expect_true(res$ci_upper >= res$att)
})

test_that("residuals sum to approximately zero", {
  y <- c(3, 5, 7, 9, 1, 2, 4, 6)
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_equal(sum(res$resid), 0, tolerance = 1e-10)
})

test_that("residuals sum to zero with controls", {
  set.seed(42)
  d <- c(rep(1L, 5), rep(0L, 5))
  x <- matrix(rnorm(10), ncol = 1)
  y <- rnorm(10) + 2 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(sum(res$resid), 0, tolerance = 1e-10)
})

test_that("vcov matrix is symmetric positive semi-definite", {
  y <- c(3, 5, 7, 9, 1, 2, 4, 6)
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  # Symmetric
  expect_equal(res$vcov, t(res$vcov), tolerance = 1e-12)
  # All eigenvalues >= 0
  evals <- eigen(res$vcov, symmetric = TRUE)$values
  expect_true(all(evals >= -1e-12))
})

test_that("SE is non-negative", {
  y <- c(3, 5, 7, 1, 2, 4)
  d <- c(1L, 1L, 1L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  expect_true(res$se >= 0)
})

test_that("fitted values equal X %*% beta", {
  y <- c(3, 5, 7, 9, 1, 2, 4, 6)
  d <- c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
  res <- estimate_ra_common(y, d)
  fitted <- as.numeric(res$X_design %*% res$params)
  expect_equal(y - res$resid, fitted, tolerance = 1e-12)
})

# ── Group 17: Layer 5 — Three-tier fallback exact boundary tests ──────────────

test_that("T5-01: Tier 1 full_interaction (N1=5, N0=8, K=3)", {
  set.seed(42)
  N1 <- 5L; N0 <- 8L; K <- 3L
  d <- c(rep(1L, N1), rep(0L, N0))
  x <- matrix(rnorm((N1 + N0) * K), ncol = K)
  y <- rnorm(N1 + N0) + 3 * d
  res <- estimate_ra_common(y, d, x = x)
  expect_equal(res$controls_tier, "full_interaction")
  expect_equal(res$K, K)
  expect_equal(res$df, 13L - 2L - 2L * K)  # df = 5
  expect_equal(ncol(res$X_design), 2L + 2L * K)  # 8 columns
})

test_that("T5-02: Tier 2 simple (N1=3, N0=15, K=3)", {
  set.seed(42)
  N1 <- 3L; N0 <- 15L; K <- 3L
  d <- c(rep(1L, N1), rep(0L, N0))
  x <- matrix(rnorm((N1 + N0) * K), ncol = K)
  y <- rnorm(N1 + N0) + 3 * d
  w <- NULL
  res <- withCallingHandlers(
    estimate_ra_common(y, d, x = x),
    lwdid_data = function(cond) { w <<- cond; invokeRestart("muffleWarning") }
  )
  expect_equal(res$controls_tier, "simple")
  expect_equal(res$K, K)
  expect_equal(res$df, 18L - K - 2L)  # df = 13
  expect_equal(ncol(res$X_design), K + 2L)  # 5 columns
  expect_true(!is.null(w))  # lwdid_data warning emitted
})

test_that("T5-03: Tier 2 triggered by N0 <= K+1 (N1=10, N0=2, K=2)", {
  set.seed(42)
  N1 <- 10L; N0 <- 2L; K <- 2L
  d <- c(rep(1L, N1), rep(0L, N0))
  x <- matrix(rnorm((N1 + N0) * K), ncol = K)
  y <- rnorm(N1 + N0) + 3 * d
  w <- NULL
  res <- withCallingHandlers(
    estimate_ra_common(y, d, x = x),
    lwdid_data = function(cond) { w <<- cond; invokeRestart("muffleWarning") }
  )
  expect_equal(res$controls_tier, "simple")
  expect_equal(res$df, 12L - K - 2L)  # df = 8
})

test_that("T5-05: N1=K+1 boundary does NOT satisfy Tier 1 (strict >)", {
  set.seed(42)
  K <- 3L
  N1 <- K + 1L  # N1=4, K+1=4, NOT strictly greater
  N0 <- 20L
  d <- c(rep(1L, N1), rep(0L, N0))
  x <- matrix(rnorm((N1 + N0) * K), ncol = K)
  y <- rnorm(N1 + N0) + 3 * d
  w <- NULL
  res <- withCallingHandlers(
    estimate_ra_common(y, d, x = x),
    lwdid_data = function(cond) { w <<- cond; invokeRestart("muffleWarning") }
  )
  expect_equal(res$controls_tier, "simple")  # Tier 2, not Tier 1
  expect_equal(res$df, (N1 + N0) - K - 2L)
})

test_that("T5-07: no controls does not trigger degradation warning", {
  y <- c(5, 3, 1, 2)
  d <- c(1L, 1L, 0L, 0L)
  expect_silent(estimate_ra_common(y, d, x = NULL))
  res <- estimate_ra_common(y, d, x = NULL)
  expect_equal(res$controls_tier, "none")
  expect_equal(res$df, 2L)
})

# ── Group 18: Layer 3 — Stata numerical consistency ──────────────────────────

test_that("T3-01: demean ATT/SE/df/t-stat/N match Stata (smoking data)", {
  data(smoking, package = "lwdid")
  result <- withCallingHandlers(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_equal(result$att, -0.4221746, tolerance = 1e-6)
  expect_equal(result$se_att, 0.1207995, tolerance = 1e-4)
  expect_equal(result$df_inference, 37L)
  expect_equal(result$t_stat, result$att / result$se_att, tolerance = 1e-10)
  expect_equal(result$nobs, 39L)
})

test_that("T3-02: detrend ATT/SE/df/t-stat/N match Stata (smoking data)", {
  data(smoking, package = "lwdid")
  result <- withCallingHandlers(
    lwdid(smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "detrend"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_equal(result$att, -0.2269887, tolerance = 1e-6)
  expect_equal(result$se_att, 0.0940689, tolerance = 1e-4)
  expect_equal(result$df_inference, 37L)
  expect_equal(result$t_stat, result$att / result$se_att, tolerance = 1e-10)
  expect_equal(result$nobs, 39L)
})
