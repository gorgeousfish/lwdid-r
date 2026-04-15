# ============================================================================
# test-vce-unit.R — E3-06 VCE Unit Test Suite
#
# Layers covered:
#   L1:  sandwich package consistency (HC0-HC4, cluster)
#   L4:  boundary/error handling
#   L5:  HC hierarchy (HC1 > HC0)
#   L6:  t-distribution verification
#   L8:  near-zero SE warning
#   L9:  cluster imbalance warning
#   L10: graded cluster warnings
#   L11: compute_inference parameter validation
#   L12: HC input validation
#   L13: HC2/HC3 numerical properties
#   L14: cluster variable type handling
#   L16: robust alias and vce_type fields
#
# Test IDs: L1-xx, L4-xx, L5-xx, L6-xx, L8-xx, L9-xx, L10-xx,
#           L11-xx, L12-xx, L13-xx, L14-xx, L16-xx
# ============================================================================

# --- Shared test data generators ---

#' Generate simple regression data for VCE testing
#' @param n Total sample size
#' @param n_treated Number of treated units
#' @param k Number of control variables (0 for no controls)
#' @param seed Random seed
#' @return List with y, d, x, cluster components
generate_vce_test_data <- function(n = 30, n_treated = 15, k = 2,
                                   seed = 42) {
  set.seed(seed)
  d <- c(rep(1L, n_treated), rep(0L, n - n_treated))
  x <- if (k > 0) matrix(rnorm(n * k), n, k) else NULL
  y <- 1 + 2 * d + (if (k > 0) rowSums(x) else 0) + rnorm(n)
  list(y = y, d = d, x = x)
}

#' Generate clustered data for VCE testing
#' @param g Number of clusters
#' @param obs_per_cluster Observations per cluster (scalar or vector)
#' @param seed Random seed
#' @return List with y, d, x, cluster components
generate_clustered_data <- function(g = 10, obs_per_cluster = 6,
                                    seed = 42) {
  set.seed(seed)
  if (length(obs_per_cluster) == 1) {
    obs_per_cluster <- rep(obs_per_cluster, g)
  }
  n <- sum(obs_per_cluster)
  cluster <- rep(seq_len(g), times = obs_per_cluster)
  # Cluster-level random effect
  cluster_effect <- rnorm(g, sd = 2)[cluster]
  d <- c(rep(1L, n %/% 2), rep(0L, n - n %/% 2))
  y <- 1 + 2 * d + cluster_effect + rnorm(n)
  list(y = y, d = d, x = NULL, cluster = cluster, n = n, g = g)
}

# ============================================================================
# Layer 1: Sandwich Package Consistency (HC0-HC4 + Cluster)
# ============================================================================

# L1-01: HC0-HC4 loop test against sandwich::vcovHC
test_that("L1-01: HC0-HC4 VCE matches sandwich::vcovHC exactly", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  for (hc_type in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- estimate_ra_common(dat$y, dat$d, dat$x, vce = hc_type)

    # Compare with sandwich directly
    sandwich_vcov <- sandwich::vcovHC(result$fit, type = toupper(hc_type))

    expect_equal(
      result$vcov, sandwich_vcov,
      tolerance = 1e-12,
      label = sprintf("vcov for %s", hc_type)
    )

    # Verify SE is positive and finite
    expect_true(result$se > 0, label = sprintf("SE > 0 for %s", hc_type))
    expect_true(is.finite(result$se), label = sprintf("SE finite for %s", hc_type))

    # Verify df = fit$df.residual (not cluster df)
    expect_equal(result$df, result$fit$df.residual,
                 label = sprintf("df for %s", hc_type))
  }
})

# L1-02: Cluster SE matches sandwich::vcovCL exactly
test_that("L1-02: cluster VCE matches sandwich::vcovCL exactly", {
  dat <- generate_clustered_data(g = 10, obs_per_cluster = 6, seed = 42)

  result <- suppressWarnings(
    estimate_ra_common(dat$y, dat$d, dat$x, vce = "cluster", cluster = dat$cluster)
  )

  # Compare with sandwich directly
  sandwich_vcov <- sandwich::vcovCL(result$fit, cluster = dat$cluster, type = "HC1")

  expect_equal(result$vcov, sandwich_vcov, tolerance = 1e-12)

  # Verify cluster df = G - 1
  expect_equal(result$df, dat$g - 1L)

  # Verify n_clusters = G
  expect_equal(result$n_clusters, dat$g)

  # Verify SE is positive and finite
  expect_true(result$se > 0)
  expect_true(is.finite(result$se))
})

# L1-03: Homoskedastic VCE matches vcov(fit) exactly
test_that("L1-03: homoskedastic VCE matches vcov(fit) exactly", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  result <- estimate_ra_common(dat$y, dat$d, dat$x, vce = NULL)

  # Compare with base R vcov()
  base_vcov <- vcov(result$fit)

  expect_equal(result$vcov, base_vcov, tolerance = 1e-15)
  expect_equal(result$vce_type, "homoskedastic")
  expect_null(result$n_clusters)
})

# ============================================================================
# Layer 4: Boundary Cases and Error Handling
# ============================================================================

# L4-01: Invalid VCE type throws lwdid_invalid_vce
test_that("L4-01: invalid VCE type throws lwdid_invalid_vce", {
  fit <- lm(y ~ x, data = data.frame(y = 1:5, x = 1:5))
  expect_error(compute_vce(fit, vce = "invalid"), class = "lwdid_invalid_vce")
  expect_error(compute_vce(fit, vce = "hc5"), class = "lwdid_invalid_vce")
  expect_error(compute_vce(fit, vce = ""), class = "lwdid_invalid_vce")
})

# L4-02: cluster VCE without cluster argument throws lwdid_invalid_parameter
test_that("L4-02: cluster VCE without cluster throws lwdid_invalid_parameter", {
  fit <- lm(y ~ x, data = data.frame(y = 1:5, x = 1:5))
  expect_error(compute_vce(fit, vce = "cluster"),
               class = "lwdid_invalid_parameter")
})

# L4-03: G=1 throws lwdid_insufficient_data
test_that("L4-03: G=1 throws lwdid_insufficient_data", {
  fit <- lm(y ~ x, data = data.frame(y = 1:5, x = 1:5))
  expect_error(compute_vce(fit, vce = "cluster", cluster = rep(1, 5)),
               class = "lwdid_insufficient_data")
})

# L4-04: G=2 boundary works correctly
test_that("L4-04: G=2 boundary computes correctly", {
  set.seed(42)
  n <- 10
  d <- c(rep(1, 5), rep(0, 5))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:2, each = 5)
  fit <- lm(y ~ d)
  result <- suppressWarnings(
    compute_vce(fit, vce = "cluster", cluster = cluster_id)
  )
  expect_equal(result$df, 1L)  # G - 1 = 1
  expect_equal(result$n_clusters, 2L)
  expect_true(is.matrix(result$vcov))
  expect_true(all(is.finite(diag(result$vcov))))
})

# L4-05: G=10 triggers informational warning (10<=G<20), not strong warning
test_that("L4-05: G=10 triggers informational warning", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:10, each = 3)  # G=10
  fit <- lm(y ~ d)
  # Should trigger "建议使用Wild Cluster Bootstrap" (informational)
  # NOT "高度不可靠" (strong warning for G<10)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster_id),
    class = "lwdid_small_sample"
  )
})

# L4-06: Extreme small sample N=3, df=1
test_that("L4-06: extreme small sample N=3 df=1 works", {
  y <- c(5, 2, 3)
  d <- c(1, 0, 0)
  result <- estimate_ra_common(y, d, vce = NULL)
  expect_equal(result$df, 1L)  # N - 2 = 1
  expect_true(!is.na(result$pvalue))
  expect_true(!is.na(result$ci_lower))
  expect_true(result$ci_upper > result$ci_lower)
})

# ============================================================================
# Layer 5: HC hierarchy (HC1 > HC0 strict, general ordering)
# ============================================================================

test_that("L5-01: HC1 > HC0 strict scalar relationship", {
  set.seed(42)
  n <- 15
  d <- c(rep(1, 5), rep(0, 10))
  x <- rnorm(n)
  y <- 1 + 2 * d + 3 * x + rnorm(n, sd = abs(x) + 0.1)
  fit <- lm(y ~ d + x + d:x)
  p <- length(coef(fit))

  se_hc0 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC0")))["d"]
  se_hc1 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC1")))["d"]

  # HC1 = sqrt(N/(N-p)) * HC0 — exact scalar relationship
  expect_equal(se_hc1, se_hc0 * sqrt(n / (n - p)),
               tolerance = 1e-12)
  expect_true(se_hc1 > se_hc0)
})

test_that("L5-02: small sample HC ordering HC3 >= HC2 >= HC1 >= HC0", {
  set.seed(42)
  n <- 15
  d <- c(rep(1, 5), rep(0, 10))
  x <- rnorm(n)
  y <- 1 + 2 * d + 3 * x + rnorm(n, sd = abs(x) + 0.1)
  fit <- lm(y ~ d + x + d:x)

  se_hc0 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC0")))["d"]
  se_hc1 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC1")))["d"]
  se_hc2 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC2")))["d"]
  se_hc3 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC3")))["d"]

  expect_true(se_hc1 > se_hc0)
  expect_true(se_hc2 >= se_hc1 - 1e-10)
  expect_true(se_hc3 >= se_hc2 - 1e-10)
})

# ============================================================================
# Layer 6: t-distribution verification
# ============================================================================

test_that("L6-01: p-value uses t-distribution not normal", {
  y <- c(5, 3, 1, 2, 4)
  d <- c(1, 1, 0, 0, 0)
  result <- estimate_ra_common(y, d, vce = NULL)

  # Manual t-distribution p-value
  p_t <- 2 * pt(-abs(result$t_stat), df = result$df)
  # Normal distribution p-value (should differ)
  p_z <- 2 * pnorm(-abs(result$t_stat))

  expect_equal(result$pvalue, p_t, tolerance = 1e-12)
  # Small df → t-distribution p-value > normal p-value
  if (result$df < 30) {
    expect_true(result$pvalue > p_z)
  }
})

# ============================================================================
# Layer 8: Near-Zero SE Warning
# ============================================================================

# L8-01: Perfect fit data triggers lwdid_numerical warning
test_that("L8-01: near-zero SE triggers lwdid_numerical warning", {
  # Perfect fit: y = 1 + 2*D (no noise)
  # compute_vce checks for coefficient named "D" (uppercase)
  y <- c(3, 3, 3, 1, 1, 1)
  D <- c(1, 1, 1, 0, 0, 0)
  fit <- lm(y ~ D)
  expect_warning(
    compute_vce(fit, vce = NULL),
    class = "lwdid_numerical"
  )
})

# ============================================================================
# Layer 9: Cluster Size Imbalance Warning
# ============================================================================

# L9-01: Highly unbalanced clusters trigger warning (G>=20 to avoid small-sample confusion)
test_that("L9-01: highly unbalanced clusters trigger warning", {
  set.seed(42)
  # G=25 (>=20, no small-sample warning), but extreme imbalance
  n <- 73
  d <- c(rep(1, 36), rep(0, 37))
  y <- 1 + 2 * d + rnorm(n)
  # 1 cluster with 49 obs, 24 clusters with 1 obs each → CV >> 1.0
  cluster_id <- c(rep(1, 49), 2:25)
  fit <- lm(y ~ d)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster_id),
    class = "lwdid_small_sample"
  )
})

# L9-02: Balanced clusters do NOT trigger imbalance warning
test_that("L9-02: balanced clusters no imbalance warning", {
  set.seed(42)
  n <- 60
  d <- c(rep(1, 30), rep(0, 30))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:20, each = 3)  # balanced: 3 per cluster, CV=0
  fit <- lm(y ~ d)
  # G=20: no small-sample warning; balanced: no imbalance warning
  result <- compute_vce(fit, vce = "cluster", cluster = cluster_id)
  expect_equal(result$n_clusters, 20L)
  expect_true(is.matrix(result$vcov))
})

# ============================================================================
# Layer 10: Graded Cluster Warnings
# ============================================================================

# L10-01: G<10 triggers strong warning
test_that("L10-01: G<10 triggers strong warning", {
  set.seed(42)
  n <- 15
  d <- c(rep(1, 7), rep(0, 8))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:5, each = 3)  # G=5 < 10
  fit <- lm(y ~ d)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster_id),
    regexp = "highly unreliable|高度不可靠",
    class = "lwdid_small_sample"
  )
})

# L10-02: 10<=G<20 triggers informational warning
test_that("L10-02: 10<=G<20 triggers informational warning", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:15, each = 2)  # G=15
  fit <- lm(y ~ d)
  expect_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster_id),
    regexp = "Wild Cluster Bootstrap",
    class = "lwdid_small_sample"
  )
})

# L10-03: G>=20 no small-sample warning
test_that("L10-03: G>=20 no small-sample warning", {
  set.seed(42)
  n <- 60
  d <- c(rep(1, 30), rep(0, 30))
  y <- 1 + 2 * d + rnorm(n)
  cluster_id <- rep(1:20, each = 3)  # G=20
  fit <- lm(y ~ d)
  # Should NOT trigger any lwdid_small_sample warning
  expect_no_warning(
    compute_vce(fit, vce = "cluster", cluster = cluster_id),
    class = "lwdid_small_sample"
  )
})


# ============================================================================
# Layer 5: HC hierarchy verification (via estimate_ra_common)
# ============================================================================

test_that("L5-01b: HC1 > HC0 strict relationship (scalar multiple) via estimate_ra_common", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  result_hc0 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc0")
  result_hc1 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc1")

  N <- length(dat$y)
  p <- length(coef(result_hc0$fit))

  # HC1 SE = HC0 SE * sqrt(N/(N-p)) — exact scalar relationship

  expect_equal(result_hc1$se, result_hc0$se * sqrt(N / (N - p)),
               tolerance = 1e-12)

  # HC1 SE > HC0 SE strictly
  expect_true(result_hc1$se > result_hc0$se)
})

test_that("L5-02b: General HC hierarchy HC3 >= HC2 >= HC1 >= HC0 via estimate_ra_common", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  se_hc0 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc0")$se
  se_hc1 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc1")$se
  se_hc2 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc2")$se
  se_hc3 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc3")$se

  # Strict ordering with tolerance for >=
  expect_true(se_hc3 >= se_hc2 - 1e-10)
  expect_true(se_hc2 >= se_hc1 - 1e-10)
  expect_true(se_hc1 >= se_hc0 - 1e-10)
})

# ============================================================================
# Layer 6: t-distribution verification (via estimate_ra_common)
# ============================================================================

test_that("L6-01b: p-value uses t-distribution (not normal) via estimate_ra_common", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  result <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc3")

  # Manual t-distribution p-value
  p_t <- 2 * pt(-abs(result$t_stat), result$df)
  # Normal distribution p-value
  p_norm <- 2 * pnorm(-abs(result$t_stat))

  # result$pvalue must match t-distribution exactly
  expect_equal(result$pvalue, p_t, tolerance = 1e-15)

  # result$pvalue should NOT match normal distribution
  expect_true(abs(result$pvalue - p_norm) > 1e-10)
})

test_that("L6-02b: Small df -> t-distribution p-value > normal p-value", {
  y <- c(5, 2, 3, 1)
  d <- c(1, 1, 0, 0)

  result <- estimate_ra_common(y, d, x = NULL, vce = NULL)

  # df should be N - p = 4 - 2 = 2
  expect_equal(result$df, 2L)

  # Manual p-values
  p_t <- 2 * pt(-abs(result$t_stat), result$df)
  p_norm <- 2 * pnorm(-abs(result$t_stat))

  # t-distribution has heavier tails -> larger p-values
  expect_true(p_t > p_norm)
})

test_that("L6-03b: CI width uses t-distribution quantiles", {
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)

  result <- estimate_ra_common(dat$y, dat$d, dat$x, vce = "hc3")

  # CI width should be 2 * qt(0.975, df) * se
  expected_width <- 2 * qt(0.975, result$df) * result$se
  actual_width <- result$ci_upper - result$ci_lower

  expect_equal(actual_width, expected_width, tolerance = 1e-12)
})


# ============================================================================
# Layer 11: compute_inference() Parameter Validation
# ============================================================================

# L11-01: alpha parameter validation
test_that("L11-01: alpha parameter validation", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  # alpha <= 0
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, alpha = 0),
               class = "lwdid_invalid_parameter")
  # alpha >= 1
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, alpha = 1),
               class = "lwdid_invalid_parameter")
  # alpha = NA
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, alpha = NA),
               class = "lwdid_invalid_parameter")
  # alpha non-numeric
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, alpha = "0.05"),
               class = "lwdid_invalid_parameter")
})

# L11-02: att NA/NaN throws lwdid_numerical
test_that("L11-02: att NA/NaN throws lwdid_numerical", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_error(compute_inference(att = NA_real_,
                                 vcov_mat = vcov_mat, df = 10),
               class = "lwdid_numerical")
  expect_error(compute_inference(att = NaN,
                                 vcov_mat = vcov_mat, df = 10),
               class = "lwdid_numerical")
})

# L11-03: coef_index out of bounds
test_that("L11-03: coef_index out of bounds throws error", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, coef_index = 0L),
               class = "lwdid_invalid_parameter")
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, coef_index = 3L),
               class = "lwdid_invalid_parameter")
})

# L11-04: VCE diagonal NA throws lwdid_numerical
test_that("L11-04: VCE diagonal NA throws lwdid_numerical", {
  vcov_mat <- matrix(c(NA, 0, 0, 1), nrow = 2)
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, coef_index = 1L),
               class = "lwdid_numerical")
})

# L11-05: VCE diagonal negative throws lwdid_numerical
test_that("L11-05: VCE diagonal negative throws lwdid_numerical", {
  vcov_mat <- matrix(c(-1e-15, 0, 0, 1), nrow = 2)
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 10, coef_index = 1L),
               class = "lwdid_numerical")
})

# L11-06: df <= 0 or df = NA throws lwdid_insufficient_data
test_that("L11-06: df <= 0 or NA throws lwdid_insufficient_data", {
  vcov_mat <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = 0),
               class = "lwdid_insufficient_data")
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = -1),
               class = "lwdid_insufficient_data")
  # Epic change #103: df=NA needs is.na() check before df<=0L
  expect_error(compute_inference(att = 1.0, vcov_mat = vcov_mat,
                                 df = NA),
               class = "lwdid_insufficient_data")
})

# L11-07: SE=0 warns but doesn't stop; att!=0 -> t_stat=Inf, pvalue=0
test_that("L11-07: SE=0 warns but computes; att!=0 gives Inf t_stat", {
  vcov_mat <- matrix(c(1, 0, 0, 0), nrow = 2)
  expect_warning(
    result <- compute_inference(att = 2.0, vcov_mat = vcov_mat,
                                df = 10, coef_index = 2L),
    class = "lwdid_numerical"
  )
  expect_equal(result$se, 0)
  expect_true(is.infinite(result$t_stat))
  expect_equal(result$pvalue, 0)
})

# L11-08: att=0 and se>0 -> t_stat=0, pvalue=1.0, CI symmetric
test_that("L11-08: att=0 se>0 gives t_stat=0 pvalue=1", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  result <- compute_inference(att = 0, vcov_mat = vcov_mat,
                              df = 10, coef_index = 2L)
  expect_equal(result$t_stat, 0)
  expect_equal(result$pvalue, 1.0, tolerance = 1e-12)
  # CI symmetric around 0
  expect_equal(result$ci_lower, -result$ci_upper,
               tolerance = 1e-12)
})

# L11-09: alpha=0.10 uses qt(0.95, df), 90% CI narrower than 95%
test_that("L11-09: alpha=0.10 gives narrower CI than alpha=0.05", {
  vcov_mat <- matrix(c(1, 0, 0, 0.25), nrow = 2)
  result_90 <- compute_inference(
    att = 1.0, vcov_mat = vcov_mat, df = 20,
    alpha = 0.10, coef_index = 2L
  )
  result_95 <- compute_inference(
    att = 1.0, vcov_mat = vcov_mat, df = 20,
    alpha = 0.05, coef_index = 2L
  )
  width_90 <- result_90$ci_upper - result_90$ci_lower
  width_95 <- result_95$ci_upper - result_95$ci_lower
  expect_true(width_90 < width_95)
  # Verify 90% CI uses qt(0.95, 20)
  se <- sqrt(0.25)
  t_crit_90 <- qt(0.95, 20)
  expect_equal(result_90$ci_lower, 1.0 - t_crit_90 * se,
               tolerance = 1e-12)
  expect_equal(result_90$ci_upper, 1.0 + t_crit_90 * se,
               tolerance = 1e-12)
})

# ============================================================================
# Layer 12: HC Type Input Validation and Warnings
# ============================================================================

# L12-01: Invalid HC type throws lwdid_invalid_parameter
test_that("L12-01: invalid HC type throws lwdid_invalid_parameter", {
  fit <- lm(y ~ x, data = data.frame(y = 1:10, x = 1:10))
  expect_error(compute_hc_vce(fit, type = "HC5"),
               class = "lwdid_invalid_parameter")
  expect_error(compute_hc_vce(fit, type = "invalid"),
               class = "lwdid_invalid_parameter")
})

# L12-02: HC1/HC3 small sample (N_treated<2) warns; HC0 does not
test_that("L12-02: HC1/HC3 small sample warns; HC0 does not", {
  set.seed(42)
  # N_treated=1, N_control=9
  # Variable name must be uppercase "D" because compute_hc_vce()
  # checks "D" %in% names(model.frame(fit)) (epic change #107)
  D <- c(1, rep(0, 9))
  y <- c(5, rnorm(9, mean = 2))
  fit <- lm(y ~ D)
  # HC1 should warn
  expect_warning(compute_hc_vce(fit, type = "HC1"),
                 class = "lwdid_small_sample")
  # HC3 should warn
  expect_warning(compute_hc_vce(fit, type = "HC3"),
                 class = "lwdid_small_sample")
  # HC0 should NOT warn about small sample
  expect_no_warning(
    result_hc0 <- compute_hc_vce(fit, type = "HC0"),
    class = "lwdid_small_sample"
  )
  expect_true(!is.null(result_hc0$vcov))
})

# L12-03: HC2-HC4 high leverage warns lwdid_numerical
test_that("L12-03: HC2-HC4 high leverage warns lwdid_numerical", {
  set.seed(42)
  n <- 10
  # x[1] far from others -> h_11 near 1
  x <- c(100, rnorm(n - 1))
  d <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n, sd = 0.1)
  fit <- lm(y ~ d + x)
  h <- hatvalues(fit)
  expect_true(max(h) > 0.99)
  expect_warning(compute_hc_vce(fit, type = "HC2"),
                 class = "lwdid_numerical")
  expect_warning(compute_hc_vce(fit, type = "HC3"),
                 class = "lwdid_numerical")
  expect_warning(compute_hc_vce(fit, type = "HC4"),
                 class = "lwdid_numerical")
})

# L12-04: HC4 d_i values: high leverage > 3.5, normal < 1.5
test_that("L12-04: HC4 d_i high leverage > 3.5, normal < 1.5", {
  set.seed(42)
  n <- 20
  x <- c(200, rnorm(n - 1))
  d <- c(rep(1, 10), rep(0, 10))
  y <- 1 + 2 * d + 0.3 * x + rnorm(n, sd = 0.5)
  fit <- lm(y ~ d + x)
  h <- hatvalues(fit)
  p <- length(coef(fit))
  h_bar <- p / n
  d_i <- pmin(4, h / h_bar)
  # High leverage observation d_i should be near 4
  expect_true(d_i[1] > 3.5)
  # Normal observations d_i should be near or below 1
  expect_true(mean(d_i[-1]) < 1.5)
})


# =============================================================================
# L13: HC2 unbiasedness and HC3-jackknife relationship
# =============================================================================

test_that("L13-01: HC2 unbiased under homoskedastic data", {
  set.seed(42)
  n <- 50
  d <- c(rep(1, 25), rep(0, 25))
  y <- 1 + 2 * d + rnorm(n, sd = 1)  # homoskedastic
  fit <- lm(y ~ d)
  se_ols <- sqrt(diag(vcov(fit)))["d"]
  se_hc2 <- sqrt(diag(sandwich::vcovHC(fit, type = "HC2")))["d"]
  # Under homoskedasticity, HC2 should be close to OLS SE
  # because E[e_i^2/(1-h_ii)] = sigma^2
  relative_diff <- abs(se_hc2 - se_ols) / se_ols
  expect_true(relative_diff < 0.20,
              label = sprintf("HC2 vs OLS SE relative diff=%.4f should be <0.20", relative_diff))
})

test_that("L13-02: HC3-jackknife exact relationship V_jack = (N-1)/N * V_HC3", {
  set.seed(42)
  n <- 20
  d <- c(rep(1, 10), rep(0, 10))
  x <- rnorm(n)
  y <- 1 + 2 * d + 0.5 * x + rnorm(n)
  fit <- lm(y ~ d + x)
  # HC3 variance matrix
  V_hc3 <- sandwich::vcovHC(fit, type = "HC3")
  # Manual delete-one jackknife variance
  beta_full <- coef(fit)
  beta_jack <- matrix(NA, nrow = n, ncol = length(beta_full))
  for (i in seq_len(n)) {
    fit_i <- lm(y ~ d + x, subset = -i)
    beta_jack[i, ] <- coef(fit_i)
  }
  # Jackknife variance: V_jack = (N-1)/N * sum((beta_{(-i)} - beta_full)^2)
  # CRITICAL: Must use beta_full as center point (NOT mean(beta_{(-i)}))
  # This is because HC3 derivation uses Sherman-Morrison formula:
  # beta_{(-i)} - beta_full = -(X'X)^{-1} x_i e_i / (1 - h_ii)
  V_jack <- ((n - 1) / n) *
    crossprod(sweep(beta_jack, 2, beta_full, "-"))
  # Verify V_jack ≈ (N-1)/N * V_HC3
  V_expected <- ((n - 1) / n) * V_hc3
  # Strip dimnames for comparison (crossprod has no names,
  # sandwich::vcovHC has coefficient names)
  expect_equal(
    unname(V_jack), unname(V_expected),
    tolerance = 1e-10,
    label = "jackknife variance = (N-1)/N * HC3 variance"
  )
})

# =============================================================================
# L14: Cluster variable validation
# =============================================================================

test_that("L14-01: cluster vector length mismatch raises error", {
  set.seed(42)
  fit <- lm(y ~ d, data = data.frame(y = rnorm(10), d = rep(0:1, 5)))
  # Length mismatch: fit has 10 obs, cluster has 8
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = rep(1:4, each = 2)),
    class = "lwdid_invalid_parameter"
  )
  # Also test compute_cluster_vce() directly
  expect_error(
    compute_cluster_vce(fit, cluster = rep(1:2, each = 3)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("L14-02: cluster vector with NA raises error", {
  set.seed(42)
  fit <- lm(y ~ d, data = data.frame(y = rnorm(10), d = rep(0:1, 5)))
  cluster_with_na <- c(1, 1, 2, 2, NA, 3, 3, 4, 4, 5)
  expect_error(
    compute_cluster_vce(fit, cluster = cluster_with_na),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    compute_vce(fit, vce = "cluster", cluster = cluster_with_na),
    class = "lwdid_invalid_parameter"
  )
})

test_that("L14-03: character cluster variable works correctly", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  cluster_char <- rep(c("CA", "NY", "TX", "FL", "WA"), each = 6)
  fit <- lm(y ~ d)
  result <- suppressWarnings(
    compute_vce(fit, vce = "cluster", cluster = cluster_char)
  )
  expect_equal(result$n_clusters, 5L)
  expect_equal(result$df, 4L)  # G - 1 = 4
  expect_true(is.matrix(result$vcov))
})

test_that("L14-04: numeric cluster variable works correctly", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  cluster_num <- rep(1:5, each = 6)
  fit <- lm(y ~ d)
  result <- suppressWarnings(
    compute_vce(fit, vce = "cluster", cluster = cluster_num)
  )
  expect_equal(result$n_clusters, 5L)
  expect_equal(result$df, 4L)
})

test_that("L14-05: character and numeric cluster produce identical VCE", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  fit <- lm(y ~ d)
  cluster_char <- rep(c("A", "B", "C", "D", "E"), each = 6)
  cluster_num <- rep(1:5, each = 6)
  result_char <- suppressWarnings(
    compute_vce(fit, vce = "cluster", cluster = cluster_char)
  )
  result_num <- suppressWarnings(
    compute_vce(fit, vce = "cluster", cluster = cluster_num)
  )
  expect_equal(result_char$vcov, result_num$vcov, tolerance = 1e-12)
  expect_equal(result_char$df, result_num$df)
})

# =============================================================================
# L16: robust alias and vce_type field validation
# =============================================================================

test_that("L16-01: vce='robust' equivalent to vce='hc1'", {
  set.seed(42)
  n <- 20
  d <- c(rep(1, 10), rep(0, 10))
  y <- 1 + 2 * d + rnorm(n)
  fit <- lm(y ~ d)
  result_robust <- compute_vce(fit, vce = "robust")
  result_hc1 <- compute_vce(fit, vce = "hc1")
  expect_equal(result_robust$vcov, result_hc1$vcov, tolerance = 1e-12)
})

test_that("L16-02: homoskedastic vce_type and n_clusters=NULL", {
  set.seed(42)
  fit <- lm(y ~ d, data = data.frame(y = rnorm(10), d = rep(0:1, 5)))
  result <- compute_vce(fit, vce = NULL)
  expect_equal(result$vce_type, "homoskedastic")
  expect_null(result$n_clusters)
})

test_that("L16-03: HC types have n_clusters=NULL", {
  set.seed(42)
  fit <- lm(y ~ d, data = data.frame(y = rnorm(10), d = rep(0:1, 5)))
  for (hc in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- compute_vce(fit, vce = hc)
    expect_null(result$n_clusters, label = paste("n_clusters for", hc))
  }
})

test_that("L16-04: vcov matrix dimensions match parameter count p", {
  set.seed(42)
  n <- 20; K <- 2
  d <- c(rep(1, 10), rep(0, 10))
  x <- matrix(rnorm(n * K), ncol = K)
  y <- 1 + 2 * d + rnorm(n)
  # Tier 1: p = 2 + 2K = 6 (intercept + D + 2*Xc + 2*D:Xc)
  result <- estimate_ra_common(y, d, x, vce = "hc3")
  p <- 2L + 2L * K
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
})

test_that("L16-05: vce parameter case-insensitive (BC-38/BC-39)", {
  set.seed(42)
  n <- 20
  d <- c(rep(1, 10), rep(0, 10))
  y <- 1 + 2 * d + rnorm(n)
  fit <- lm(y ~ d)
  # HC3 case variants all dispatch to HC3
  result_lower <- compute_vce(fit, vce = "hc3")
  result_upper <- compute_vce(fit, vce = "HC3")
  result_mixed <- compute_vce(fit, vce = "Hc3")
  expect_equal(result_lower$vcov, result_upper$vcov, tolerance = 1e-12)
  expect_equal(result_lower$vcov, result_mixed$vcov, tolerance = 1e-12)
  # robust case variants all dispatch to HC1
  result_robust_lower <- compute_vce(fit, vce = "robust")
  result_robust_upper <- compute_vce(fit, vce = "ROBUST")
  result_robust_mixed <- compute_vce(fit, vce = "Robust")
  result_hc1 <- compute_vce(fit, vce = "hc1")
  expect_equal(result_robust_lower$vcov, result_hc1$vcov, tolerance = 1e-12)
  expect_equal(result_robust_upper$vcov, result_hc1$vcov, tolerance = 1e-12)
  expect_equal(result_robust_mixed$vcov, result_hc1$vcov, tolerance = 1e-12)
})

test_that("L16-06: homoskedastic VCE matches vcov(fit) exactly", {
  set.seed(42)
  n <- 20
  d <- c(rep(1, 10), rep(0, 10))
  y <- 1 + 2 * d + rnorm(n)
  fit <- lm(y ~ d)
  result <- compute_vce(fit, vce = NULL)
  # Homoskedastic VCE should exactly match lm()'s built-in vcov()
  expect_equal(result$vcov, vcov(fit), tolerance = 1e-15)
})


# ============================================================================
# Layer 16: robust alias and vce_type fields
# ============================================================================

test_that("L16-01: robust and hc1 produce identical VCE", {
  set.seed(42)
  n <- 30
  d <- c(rep(1, 15), rep(0, 15))
  y <- 1 + 2 * d + rnorm(n)
  fit <- lm(y ~ d)
  res_robust <- compute_vce(fit, vce = "robust")
  res_hc1 <- compute_vce(fit, vce = "hc1")
  expect_equal(res_robust$vcov, res_hc1$vcov, tolerance = 1e-12)
})

test_that("L16-02: homoskedastic vce_type='homoskedastic' n_clusters=NULL", {
  set.seed(42)
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 0, seed = 42)
  result <- estimate_ra_common(dat$y, dat$d, vce = NULL)
  expect_identical(result$vce_type, "homoskedastic")
  expect_null(result$n_clusters)
})

test_that("L16-03: HC types have n_clusters = NULL", {
  set.seed(42)
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 0, seed = 42)
  for (hc in c("hc0", "hc1", "hc2", "hc3", "hc4")) {
    result <- estimate_ra_common(dat$y, dat$d, vce = hc)
    expect_null(result$n_clusters, label = paste("n_clusters for", hc))
  }
})

test_that("L16-04: vcov dimension matches parameter count", {
  set.seed(42)
  # Tier 3: p=2
  result_t3 <- estimate_ra_common(c(1,2,3,4), c(1,1,0,0), vce = NULL)
  expect_equal(nrow(result_t3$vcov), 2L)
  expect_equal(ncol(result_t3$vcov), 2L)
  # Tier 1: p=2+2K
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)
  result_t1 <- estimate_ra_common(dat$y, dat$d, dat$x, vce = NULL)
  expect_equal(nrow(result_t1$vcov), 2L + 2L * 2L)  # 6
  expect_equal(ncol(result_t1$vcov), 2L + 2L * 2L)
})

test_that("L16-05: vce parameter case insensitive", {
  set.seed(42)
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 0, seed = 42)
  # HC3 variants
  res_lower <- estimate_ra_common(dat$y, dat$d, vce = "hc3")
  res_upper <- estimate_ra_common(dat$y, dat$d, vce = "HC3")
  res_mixed <- estimate_ra_common(dat$y, dat$d, vce = "Hc3")
  expect_equal(res_lower$vcov, res_upper$vcov, tolerance = 1e-12)
  expect_equal(res_lower$vcov, res_mixed$vcov, tolerance = 1e-12)
  # robust variants
  res_r1 <- estimate_ra_common(dat$y, dat$d, vce = "robust")
  res_r2 <- estimate_ra_common(dat$y, dat$d, vce = "ROBUST")
  res_r3 <- estimate_ra_common(dat$y, dat$d, vce = "Robust")
  expect_equal(res_r1$vcov, res_r2$vcov, tolerance = 1e-12)
  expect_equal(res_r1$vcov, res_r3$vcov, tolerance = 1e-12)
})

test_that("L16-06: homoskedastic VCE = vcov(fit) exactly", {
  set.seed(42)
  dat <- generate_vce_test_data(n = 30, n_treated = 15, k = 2, seed = 42)
  result <- estimate_ra_common(dat$y, dat$d, dat$x, vce = NULL)
  expect_equal(result$vcov, vcov(result$fit), tolerance = 1e-15)
})
