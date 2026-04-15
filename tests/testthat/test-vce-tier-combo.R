# ============================================================================
# Tests for Three-tier Degradation + VCE Combination (Task E3-05.5)
# Story E3-05: VCE Integration Tests for estimate_ra_common()
#
# Comprehensive tests verifying correct interaction between the
# three-tier control degradation system and VCE types:
#   T5b-01: Tier 1 + vce=NULL
#   T5b-02: Tier 2 + vce="hc3"
#   T5b-03: Tier 3 + vce=NULL
#   T5b-04: Tier 2 + vce="cluster" (df = G-1)
#   T5b-05: Tier 1 + vce="cluster" (df = G-1)
#   T5b-06: Degradation warnings contain actual N_1/N_0/N values
#   T5b-07: fit$df.residual automatically reflects correct tier
# ============================================================================

library(testthat)
library(sandwich)

# ============================================================================
# T5b-01: Tier 1 + vce=NULL: df = N - 2 - 2K, controls_tier = "full_interaction"
# ============================================================================
test_that("T5b-01: Tier 1 + vce=NULL gives full_interaction, df = N-2-2K", {
  set.seed(100)
  N <- 30L
  K <- 2L
  n1 <- 15L
  n0 <- 15L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(N, 0, 0.5)

  result <- estimate_ra_common(y, d, x, vce = NULL)

  # Tier verification
  expect_equal(result$controls_tier, "full_interaction")

  # df = N - 2 - 2K = 30 - 2 - 4 = 24
  expected_df <- N - 2L - 2L * K
  expect_equal(result$df, expected_df)

  # df matches lm() internal df.residual
  expect_equal(result$df, result$fit$df.residual)

  # Numerical reasonableness
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$att)
  expect_true(result$ci_upper > result$att)

  # ATT should be near true value of 2
  expect_true(abs(result$att - 2) < 2)

  # VCE type
  expect_equal(result$vce_type, "homoskedastic")
  expect_null(result$n_clusters)
})

# ============================================================================
# T5b-02: Tier 2 + vce="hc3": df = N - K - 2, controls_tier = "simple"
# ============================================================================
test_that("T5b-02: Tier 2 + vce=hc3 gives simple, df = N-K-2", {
  set.seed(200)
  N <- 10L
  K <- 2L
  n1 <- 2L
  n0 <- 8L
  # n1=2 <= K+1=3, so NOT Tier 1; N=10 > K+2=4, so Tier 2
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, 0, 0.5)

  # Should emit degradation warning
  expect_warning(
    result <- estimate_ra_common(y, d, x, vce = "hc3"),
    regexp = "degraded to simple"
  )

  # Tier verification
  expect_equal(result$controls_tier, "simple")

  # df = N - K - 2 = 10 - 2 - 2 = 6
  expected_df <- N - K - 2L
  expect_equal(result$df, expected_df)

  # df matches lm() internal df.residual
  expect_equal(result$df, result$fit$df.residual)

  # VCE type is HC3
  expect_equal(result$vce_type, "HC3")

  # Verify SE matches sandwich::vcovHC manually
  d_idx <- which(names(coef(result$fit)) == "D")
  manual_vcov <- sandwich::vcovHC(result$fit, type = "HC3")
  manual_se <- sqrt(manual_vcov[d_idx, d_idx])
  expect_equal(result$se, manual_se, tolerance = 1e-12)

  # Numerical reasonableness
  expect_true(result$se > 0)
  expect_null(result$n_clusters)
})

# ============================================================================
# T5b-03: Tier 3 + vce=NULL: df = N - 2, controls_tier = "none"
# ============================================================================
test_that("T5b-03: Tier 3 + vce=NULL gives none, df = N-2", {
  set.seed(300)
  N <- 5L
  K <- 3L
  n1 <- 2L
  n0 <- 3L
  # N=5 <= K+2=5, triggers Tier 3 (controls dropped)
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, 0, 0.5)

  # Should emit controls_dropped warning
  expect_warning(
    result <- estimate_ra_common(y, d, x, vce = NULL),
    regexp = "controls dropped"
  )

  # Tier verification
  expect_equal(result$controls_tier, "none")

  # df = N - 2 = 5 - 2 = 3
  expected_df <- N - 2L
  expect_equal(result$df, expected_df)

  # df matches lm() internal df.residual
  expect_equal(result$df, result$fit$df.residual)

  # K should be 0 (controls dropped)
  expect_equal(result$K, 0L)

  # VCE type
  expect_equal(result$vce_type, "homoskedastic")

  # Numerical reasonableness
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$att)
  expect_true(result$ci_upper > result$att)
})

# ============================================================================
# T5b-04: Tier 2 + vce="cluster": df = G - 1 (NOT N - K - 2)
# ============================================================================
test_that("T5b-04: Tier 2 + cluster VCE gives df = G-1, not N-K-2", {
  set.seed(400)
  N <- 20L
  K <- 2L
  n1 <- 2L
  n0 <- 18L
  G <- 10L
  # n1=2 <= K+1=3, so Tier 2; N=20 > K+2=4
  cluster <- rep(1:G, each = 2L)
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, 0, 0.5)

  # Suppress both tier degradation and small-cluster warnings
  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "cluster", cluster = cluster)
  )

  # Tier verification
  expect_equal(result$controls_tier, "simple")

  # df = G - 1 = 9 (NOT N - K - 2 = 16)
  expect_equal(result$df, G - 1L)
  expect_false(result$df == N - K - 2L)

  # Cluster info
  expect_equal(result$n_clusters, G)
  expect_equal(result$vce_type, "cluster")

  # Numerical reasonableness
  expect_true(result$se > 0)
})

# ============================================================================
# T5b-05: Tier 1 + vce="cluster": df = G - 1 (NOT N - 2 - 2K)
# ============================================================================
test_that("T5b-05: Tier 1 + cluster VCE gives df = G-1, not N-2-2K", {
  set.seed(500)
  N <- 40L
  K <- 2L
  G <- 20L
  # Each cluster has 1 treated, 1 control -> n1=20, n0=20
  # n1=20 > K+1=3 and n0=20 > K+1=3 -> Tier 1
  cluster <- rep(1:G, each = 2L)
  d <- rep(c(1L, 0L), G)
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(N, 0, 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "cluster", cluster = cluster)
  )

  # Tier verification
  expect_equal(result$controls_tier, "full_interaction")

  # df = G - 1 = 19 (NOT N - 2 - 2K = 34)
  expect_equal(result$df, G - 1L)
  expect_false(result$df == N - 2L - 2L * K)

  # Cluster info
  expect_equal(result$n_clusters, G)
  expect_equal(result$vce_type, "cluster")

  # Numerical reasonableness
  expect_true(result$se > 0)
  expect_true(result$ci_lower < result$att)
  expect_true(result$ci_upper > result$att)
})

# ============================================================================
# T5b-06: Degradation warnings contain actual N_1/N_0/N values
# ============================================================================
test_that("T5b-06a: Tier 2 degradation warning contains N_1 and N_0", {
  set.seed(600)
  N <- 10L
  K <- 2L
  n1 <- 2L
  n0 <- 8L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, 0, 0.5)

  # Capture the warning
  captured <- NULL
  withCallingHandlers(
    estimate_ra_common(y, d, x, vce = NULL),
    lwdid_data = function(w) {
      captured <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_false(is.null(captured))
  msg <- conditionMessage(captured)
  # Warning should contain actual sample sizes
  expect_true(grepl("N_1=2", msg, fixed = TRUE))
  expect_true(grepl("N_0=8", msg, fixed = TRUE))
})

test_that("T5b-06b: Tier 3 degradation warning contains actual N", {
  set.seed(601)
  N <- 4L
  K <- 3L
  n1 <- 2L
  n0 <- 2L
  d <- c(rep(1L, n1), rep(0L, n0))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, 0, 0.5)

  # Capture the warning
  captured <- NULL
  withCallingHandlers(
    estimate_ra_common(y, d, x, vce = NULL),
    lwdid_data = function(w) {
      captured <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_false(is.null(captured))
  msg <- conditionMessage(captured)
  # Warning should contain actual N
  expect_true(grepl("N=4", msg, fixed = TRUE))
})

# ============================================================================
# T5b-07: fit$df.residual automatically reflects correct tier
# ============================================================================
test_that("T5b-07a: Tier 1 fit$df.residual = N - 2 - 2K", {
  set.seed(701)
  N <- 30L
  K <- 2L
  d <- c(rep(1L, 15), rep(0L, 15))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, 0, 0.5)

  result <- estimate_ra_common(y, d, x, vce = NULL)

  expect_equal(result$controls_tier, "full_interaction")
  expect_equal(result$fit$df.residual, N - 2L - 2L * K)
  expect_equal(result$df, result$fit$df.residual)
})

test_that("T5b-07b: Tier 2 fit$df.residual = N - K - 2", {
  set.seed(702)
  N <- 10L
  K <- 2L
  d <- c(rep(1L, 2), rep(0L, 8))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, 0, 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = NULL)
  )

  expect_equal(result$controls_tier, "simple")
  expect_equal(result$fit$df.residual, N - K - 2L)
  expect_equal(result$df, result$fit$df.residual)
})

test_that("T5b-07c: Tier 3 fit$df.residual = N - 2", {
  set.seed(703)
  N <- 5L
  K <- 3L
  d <- c(rep(1L, 2), rep(0L, 3))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, 0, 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = NULL)
  )

  expect_equal(result$controls_tier, "none")
  expect_equal(result$fit$df.residual, N - 2L)
  expect_equal(result$df, result$fit$df.residual)
})
