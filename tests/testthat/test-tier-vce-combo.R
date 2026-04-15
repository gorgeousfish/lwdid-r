# ============================================================================
# Tests for Three-Tier Degradation + VCE Combination (Task E3-05.5)
# Story E3-05: VCE Integration Tests
#
# These tests verify that the three-tier model degradation strategy works
# correctly in combination with different VCE types. Each test verifies
# the correct tier selection, df calculation, and VCE behavior.
#
# Three-Tier Strategy:
#   Tier 1 (full_interaction): N1 > K+1 AND N0 > K+1
#     Design: [1, D, Xc, D*Xc], df = N - 2 - 2K
#   Tier 2 (simple): N > K+2 but Tier 1 conditions not met
#     Design: [1, D, X], df = N - K - 2
#   Tier 3 (none): N <= K+2
#     Design: [1, D], df = N - 2
#
# Test Blocks:
#   TV-01: Tier 1 + vce=NULL → df = N - 2 - 2K
#   TV-02: Tier 2 + vce="hc3" → df = N - K - 2 (NOT N - 2 - 2K)
#   TV-03: Tier 3 + vce=NULL → df = N - 2
#   TV-04: Tier 2 + vce="cluster" → df = G - 1 (NOT N - K - 2)
#   TV-05: Tier 1 + vce="cluster" → df = G - 1 (NOT N - 2 - 2K)
#   TV-06: Tier 2 degradation warning with actual N1/N0 values
#   TV-07: Tier 3 degradation warning with actual N value
#   TV-08: fit$df.residual matches tier automatically
#   TV-09: Tier 1 + HC0/HC1/HC2/HC3/HC4 all match sandwich
#   TV-10: Tier 2 + HC3 SE matches sandwich for the simple model
# ============================================================================

library(testthat)
library(sandwich)

# ============================================================================
# TV-01: Tier 1 + vce=NULL: df = N - 2 - 2K, controls_tier = "full_interaction"
# Data: N=20, K=2, N1=10, N0=10
# Tier 1 conditions: N1=10 > K+1=3 AND N0=10 > K+1=3 → full_interaction
# df = 20 - 2 - 2*2 = 14
# ============================================================================
test_that("TV-01: Tier 1 + vce=NULL: df = N - 2 - 2K", {
  set.seed(42)
  N <- 20L; K <- 2L
  d <- c(rep(1L, 10), rep(0L, 10))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(N, sd = 0.5)

  result <- estimate_ra_common(y, d, x, vce = NULL)

  # Tier selection
  expect_identical(result$controls_tier, "full_interaction")

  # df = N - 2 - 2K = 20 - 2 - 4 = 14
  expect_equal(result$df, 14L)

  # VCE type
  expect_identical(result$vce_type, "homoskedastic")

  # n_clusters must be NULL for non-cluster VCE
  expect_null(result$n_clusters)

  # Numerical reasonableness
  expect_true(result$se > 0)
  expect_true(result$ci_lower <= result$att)
  expect_true(result$ci_upper >= result$att)
})

# ============================================================================
# TV-02: Tier 2 + vce="hc3": df = N - K - 2 (NOT N - 2 - 2K)
# Data: N=8, K=2, N1=2, N0=6
# Tier 2 conditions: N1=2 NOT > K+1=3, but N=8 > K+2=4 → simple
# df = 8 - 2 - 2 = 4 (NOT 8 - 2 - 4 = 2)
# Verify SE matches sandwich::vcovHC(fit, type="HC3") for the simple model
# ============================================================================
test_that("TV-02: Tier 2 + vce='hc3': df = N - K - 2", {
  set.seed(42)
  N <- 8L; K <- 2L
  d <- c(rep(1L, 2), rep(0L, 6))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, sd = 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "hc3")
  )

  # Tier selection
  expect_identical(result$controls_tier, "simple")

  # df = N - K - 2 = 8 - 2 - 2 = 4 (NOT N - 2 - 2K = 8 - 2 - 4 = 2)
  expect_equal(result$df, 4L)

  # VCE type
  expect_identical(result$vce_type, "HC3")

  # SE must match sandwich::vcovHC on the simple model [1, D, X1, X2]
  fit_ref <- lm(y ~ d + x[, 1] + x[, 2])
  se_hc3_ref <- sqrt(sandwich::vcovHC(fit_ref, type = "HC3")[2, 2])
  expect_equal(result$se, se_hc3_ref, tolerance = 1e-12)
})

# ============================================================================
# TV-03: Tier 3 + vce=NULL: df = N - 2, controls_tier = "none"
# Data: N=4, K=3, N1=2, N0=2
# Tier 3 conditions: N=4 <= K+2=5 → none (all controls dropped)
# df = 4 - 2 = 2
# ============================================================================
test_that("TV-03: Tier 3 + vce=NULL: df = N - 2", {
  set.seed(42)
  N <- 4L; K <- 3L
  d <- c(rep(1L, 2), rep(0L, 2))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, sd = 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = NULL)
  )

  # Tier selection
  expect_identical(result$controls_tier, "none")

  # df = N - 2 = 4 - 2 = 2
  expect_equal(result$df, 2L)

  # VCE type
  expect_identical(result$vce_type, "homoskedastic")

  # K should be reset to 0 for Tier 3
  expect_equal(result$K, 0L)

  # Numerical reasonableness
  expect_true(result$se > 0)
})


# ============================================================================
# TV-04: Tier 2 + vce="cluster": df = G - 1 (NOT N - K - 2)
# Data: N=24, K=2, N1=3, N0=21, G=8 clusters
# Tier 2 conditions: N1=3 NOT > K+1=3, but N=24 > K+2=4 → simple
# df = G - 1 = 7 (NOT N - K - 2 = 20)
# ============================================================================
test_that("TV-04: Tier 2 + vce='cluster': df = G - 1", {
  set.seed(42)
  N <- 24L; K <- 2L; G <- 8L
  cluster <- rep(1:G, each = N %/% G)  # 3 obs per cluster
  d <- c(rep(1L, 3), rep(0L, 21))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, sd = 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "cluster", cluster = cluster)
  )

  # Tier selection
  expect_identical(result$controls_tier, "simple")

  # df = G - 1 = 7 (NOT N - K - 2 = 20)
  expect_equal(result$df, G - 1L)
  expect_false(result$df == N - K - 2L,
               label = "df must be G-1, not N-K-2")

  # VCE type and cluster count
  expect_identical(result$vce_type, "cluster")
  expect_identical(result$n_clusters, G)

  # SE must match sandwich::vcovCL on the simple model
  fit_ref <- lm(y ~ d + x[, 1] + x[, 2])
  se_cl_ref <- sqrt(sandwich::vcovCL(fit_ref, cluster = cluster,
                                      type = "HC1")[2, 2])
  expect_equal(result$se, se_cl_ref, tolerance = 1e-12)
})

# ============================================================================
# TV-05: Tier 1 + vce="cluster": df = G - 1 (NOT N - 2 - 2K)
# Data: N=60, K=2, N1=30, N0=30, G=20 clusters
# Tier 1 conditions: N1=30 > K+1=3 AND N0=30 > K+1=3 → full_interaction
# df = G - 1 = 19 (NOT N - 2 - 2K = 54)
# ============================================================================
test_that("TV-05: Tier 1 + vce='cluster': df = G - 1", {
  set.seed(42)
  N <- 60L; K <- 2L; G <- 20L
  cluster <- rep(1:G, each = N %/% G)  # 3 obs per cluster
  d <- c(rep(1L, 30), rep(0L, 30))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(N, sd = 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "cluster", cluster = cluster)
  )

  # Tier selection
  expect_identical(result$controls_tier, "full_interaction")

  # df = G - 1 = 19 (NOT N - 2 - 2K = 54)
  expect_equal(result$df, G - 1L)
  expect_false(result$df == N - 2L - 2L * K,
               label = "df must be G-1, not N-2-2K")

  # VCE type and cluster count
  expect_identical(result$vce_type, "cluster")
  expect_identical(result$n_clusters, G)

  # SE must match sandwich::vcovCL on the full interaction model
  # Reconstruct the Tier 1 model: y ~ D + Xc1 + Xc2 + D*Xc1 + D*Xc2
  x_bar1 <- colMeans(x[d == 1, , drop = FALSE])
  x_c <- sweep(x, 2, x_bar1, "-")
  fit_ref <- lm(y ~ d + x_c[, 1] + x_c[, 2] +
                  I(d * x_c[, 1]) + I(d * x_c[, 2]))
  se_cl_ref <- sqrt(sandwich::vcovCL(fit_ref, cluster = cluster,
                                      type = "HC1")[2, 2])
  expect_equal(result$se, se_cl_ref, tolerance = 1e-12)
})

# ============================================================================
# TV-06: Tier 2 degradation warning with actual N1/N0 values
# Verify lwdid_data warning is emitted with detail = "controls_degraded_to_simple"
# ============================================================================
test_that("TV-06: Tier 2 degradation warning with actual N1/N0 values", {
  set.seed(42)
  N <- 8L; K <- 2L
  d <- c(rep(1L, 2), rep(0L, 6))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, sd = 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, x, vce = NULL),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # Filter for lwdid_data warnings
  data_warnings <- Filter(
    function(w) inherits(w, "lwdid_data"),
    warnings_caught
  )
  expect_true(length(data_warnings) >= 1L,
              info = "Expected at least one lwdid_data warning for Tier 2 degradation")

  # Check detail field
  tier2_warnings <- Filter(
    function(w) identical(w$detail, "controls_degraded_to_simple"),
    data_warnings
  )
  expect_true(length(tier2_warnings) >= 1L,
              info = "Expected detail = 'controls_degraded_to_simple'")

  # Check that the warning message contains actual N1/N0 values
  msg <- conditionMessage(tier2_warnings[[1]])
  expect_true(grepl("N_1=2", msg),
              info = "Warning should contain actual N_1 value")
  expect_true(grepl("N_0=6", msg),
              info = "Warning should contain actual N_0 value")
})

# ============================================================================
# TV-07: Tier 3 degradation warning with actual N value
# Verify lwdid_data warning is emitted with detail = "controls_dropped"
# ============================================================================
test_that("TV-07: Tier 3 degradation warning with actual N value", {
  set.seed(42)
  N <- 4L; K <- 3L
  d <- c(rep(1L, 2), rep(0L, 2))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + rnorm(N, sd = 0.5)

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_ra_common(y, d, x, vce = NULL),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # Filter for lwdid_data warnings
  data_warnings <- Filter(
    function(w) inherits(w, "lwdid_data"),
    warnings_caught
  )
  expect_true(length(data_warnings) >= 1L,
              info = "Expected at least one lwdid_data warning for Tier 3 degradation")

  # Check detail field
  tier3_warnings <- Filter(
    function(w) identical(w$detail, "controls_dropped"),
    data_warnings
  )
  expect_true(length(tier3_warnings) >= 1L,
              info = "Expected detail = 'controls_dropped'")

  # Check that the warning message contains actual N value
  msg <- conditionMessage(tier3_warnings[[1]])
  expect_true(grepl("N=4", msg),
              info = "Warning should contain actual N value")
})


# ============================================================================
# TV-08: fit$df.residual matches tier automatically
# Verify that lm()'s df.residual automatically reflects the correct model
# tier without manual calculation.
#   Tier 1: fit$df.residual = N - 2 - 2K
#   Tier 2: fit$df.residual = N - K - 2
#   Tier 3: fit$df.residual = N - 2
# ============================================================================
test_that("TV-08: fit$df.residual matches tier automatically", {
  # --- Tier 1: N=20, K=2, N1=10, N0=10 ---
  set.seed(42)
  N1 <- 20L; K1 <- 2L
  d1 <- c(rep(1L, 10), rep(0L, 10))
  x1 <- matrix(rnorm(N1 * K1), N1, K1)
  y1 <- 1 + 2 * d1 + 0.5 * x1[, 1] - 0.3 * x1[, 2] + rnorm(N1, sd = 0.5)

  r1 <- estimate_ra_common(y1, d1, x1, vce = NULL)
  expect_identical(r1$controls_tier, "full_interaction")
  # Tier 1: df.residual = N - 2 - 2K = 20 - 2 - 4 = 14
  expect_equal(r1$fit$df.residual, N1 - 2L - 2L * K1)
  expect_equal(r1$df, r1$fit$df.residual)

  # --- Tier 2: N=8, K=2, N1=2, N0=6 ---
  set.seed(42)
  N2 <- 8L; K2 <- 2L
  d2 <- c(rep(1L, 2), rep(0L, 6))
  x2 <- matrix(rnorm(N2 * K2), N2, K2)
  y2 <- 1 + 2 * d2 + 0.5 * x2[, 1] + rnorm(N2, sd = 0.5)

  r2 <- suppressWarnings(estimate_ra_common(y2, d2, x2, vce = NULL))
  expect_identical(r2$controls_tier, "simple")
  # Tier 2: df.residual = N - K - 2 = 8 - 2 - 2 = 4
  expect_equal(r2$fit$df.residual, N2 - K2 - 2L)
  expect_equal(r2$df, r2$fit$df.residual)

  # --- Tier 3: N=4, K=3, N1=2, N0=2 ---
  set.seed(42)
  N3 <- 4L; K3 <- 3L
  d3 <- c(rep(1L, 2), rep(0L, 2))
  x3 <- matrix(rnorm(N3 * K3), N3, K3)
  y3 <- 1 + 2 * d3 + rnorm(N3, sd = 0.5)

  r3 <- suppressWarnings(estimate_ra_common(y3, d3, x3, vce = NULL))
  expect_identical(r3$controls_tier, "none")
  # Tier 3: df.residual = N - 2 = 4 - 2 = 2
  expect_equal(r3$fit$df.residual, N3 - 2L)
  expect_equal(r3$df, r3$fit$df.residual)
})

# ============================================================================
# TV-09: Tier 1 + HC0/HC1/HC2/HC3/HC4: all HC types work, SE matches sandwich
# Quick check: for each HC type, SE matches sandwich::vcovHC(fit, type=...)
# exactly on the full interaction model.
# ============================================================================
test_that("TV-09: Tier 1 + all HC types match sandwich", {
  set.seed(42)
  N <- 20L; K <- 2L
  d <- c(rep(1L, 10), rep(0L, 10))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(N, sd = 0.5)

  # Reconstruct the Tier 1 reference model
  x_bar1 <- colMeans(x[d == 1, , drop = FALSE])
  x_c <- sweep(x, 2, x_bar1, "-")
  fit_ref <- lm(y ~ d + x_c[, 1] + x_c[, 2] +
                  I(d * x_c[, 1]) + I(d * x_c[, 2]))

  hc_types <- c("hc0", "hc1", "hc2", "hc3", "hc4")

  for (hc in hc_types) {
    result <- estimate_ra_common(y, d, x, vce = hc)

    # Verify Tier 1
    expect_identical(result$controls_tier, "full_interaction",
                     info = paste("Tier check for", hc))

    # SE must match sandwich::vcovHC on the reference model
    se_ref <- sqrt(sandwich::vcovHC(fit_ref, type = toupper(hc))[2, 2])
    expect_equal(result$se, se_ref, tolerance = 1e-12,
                 info = paste("SE match for", hc))

    # vce_type must be uppercase
    expect_identical(result$vce_type, toupper(hc),
                     info = paste("vce_type for", hc))
  }
})

# ============================================================================
# TV-10: Tier 2 + HC3: SE matches sandwich for the simple model
# Verify the SE is computed on the correct (degraded) model, not the
# interaction model. The simple model is [1, D, X1, X2], not
# [1, D, Xc1, Xc2, D*Xc1, D*Xc2].
# ============================================================================
test_that("TV-10: Tier 2 + HC3 SE matches sandwich for simple model", {
  set.seed(42)
  N <- 8L; K <- 2L
  d <- c(rep(1L, 2), rep(0L, 6))
  x <- matrix(rnorm(N * K), N, K)
  y <- 1 + 2 * d + 0.5 * x[, 1] + rnorm(N, sd = 0.5)

  result <- suppressWarnings(
    estimate_ra_common(y, d, x, vce = "hc3")
  )

  # Verify Tier 2
  expect_identical(result$controls_tier, "simple")

  # Reference: simple model [1, D, X1, X2] (NOT interaction model)
  fit_simple <- lm(y ~ d + x[, 1] + x[, 2])
  se_hc3_simple <- sqrt(sandwich::vcovHC(fit_simple, type = "HC3")[2, 2])

  # SE must match the simple model's HC3 SE
  expect_equal(result$se, se_hc3_simple, tolerance = 1e-12)

  # Verify the full VCE matrix values match the simple model's VCE
  # (ignore dimnames since estimate_ra_common uses D/X1/X2 vs d/x[,1]/x[,2])
  vcov_ref <- sandwich::vcovHC(fit_simple, type = "HC3")
  expect_equal(unname(result$vcov), unname(vcov_ref), tolerance = 1e-12)

  # Verify df is for the simple model (N - K - 2), not interaction (N - 2 - 2K)
  expect_equal(result$df, N - K - 2L)  # 8 - 2 - 2 = 4
  expect_false(result$df == N - 2L - 2L * K)  # NOT 8 - 2 - 4 = 2
})
