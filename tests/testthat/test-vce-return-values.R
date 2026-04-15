# test-vce-return-values.R - E3-05.8 Return Value Completeness Tests
# ============================================================================
# Tests for return value completeness and structural correctness of:
#   - estimate_ra_common(): 21-field named list with lm object
#   - estimate_period_effects(): 15-column data.frame per period
#
# Covers:
#   T5e-01: estimate_ra_common() returns complete lm object
#   T5e-02: lm object contains sandwich-required components
#   T5e-03: estimate_period_effects() returns 15 columns with correct names
#   T5e-04: vcov matrix dimensions match parameter count p by tier
#   T5e-05: Non-cluster VCE in period results: n_clusters is NA_integer_
#   T5e-06: vce_type field correctly reflects VCE type
#   T5e-07: Non-cluster VCE n_clusters is NULL in estimate_ra_common return
# ============================================================================
library(testthat)
library(data.table)

# ============================================================================
# Helper: generate simple test data for estimate_ra_common()
# ============================================================================
make_ra_data <- function(n = 100, n_treated = 30, K = 0, seed = 5800) {
  set.seed(seed)
  d <- c(rep(1L, n_treated), rep(0L, n - n_treated))
  y <- rnorm(n) + 2.0 * d
  if (K > 0) {
    x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  } else {
    x <- NULL
  }
  list(y = y, d = d, x = x)
}

# ============================================================================
# Helper: generate panel data for estimate_period_effects()
# ============================================================================
make_period_data <- function(n_per_period = 60, n_treated_per = 20,
                             periods = 1:3, K = 0, n_clusters = NULL,
                             seed = 5801) {
  set.seed(seed)
  dt_list <- lapply(periods, function(r) {
    n <- n_per_period
    d <- c(rep(1L, n_treated_per), rep(0L, n - n_treated_per))
    y <- rnorm(n) + 1.5 * d
    row <- data.table(tindex = r, y_trans = y, d = d)
    if (K > 0) {
      for (k in seq_len(K)) {
        set(row, j = paste0("x", k), value = rnorm(n))
      }
    }
    if (!is.null(n_clusters)) {
      row[, cluster := rep(seq_len(n_clusters), length.out = n)]
    }
    row
  })
  rbindlist(dt_list)
}


# ============================================================================
# T5e-01: estimate_ra_common() returns complete lm object
# ============================================================================

test_that("T5e-01a: fit element is an lm object (no controls)", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_s3_class(result$fit, "lm")
})

test_that("T5e-01b: fit element is an lm object (with controls, Tier 1)", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  expect_s3_class(result$fit, "lm")
  expect_equal(result$controls_tier, "full_interaction")
})

test_that("T5e-01c: fit element is an lm object (HC3 VCE)", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_s3_class(result$fit, "lm")
})

test_that("T5e-01d: fit element is an lm object (cluster VCE)", {
  set.seed(5801)
  n <- 60
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL, vce = "cluster", cluster = cl)
  expect_s3_class(result$fit, "lm")
})

test_that("T5e-01e: return value is a named list with all 21 fields", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")

  expected_names <- c("att", "se", "t_stat", "df", "pvalue",
                      "ci_lower", "ci_upper",
                      "params", "vcov", "resid", "fit",
                      "n", "n_treated", "n_control", "K",
                      "vce_type", "n_clusters", "controls_tier", "X_design",
                      "df_resid", "df_inference")
  expect_true(is.list(result))
  expect_equal(length(expected_names), 21L)
  for (nm in expected_names) {
    expect_true(nm %in% names(result),
                info = sprintf("Missing field: %s", nm))
  }
  # No extra fields
  expect_equal(length(result), 21L)
})


# ============================================================================
# T5e-02: lm object contains sandwich-required components
# ============================================================================
# sandwich::vcovHC and sandwich::vcovCL require specific components from
# the lm object: model (model frame), qr (QR decomposition), residuals,
# and df.residual. Without these, sandwich VCE computation will fail.

test_that("T5e-02a: lm object has model component (model frame)", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$model),
              info = "lm object must contain $model for sandwich::bread()")
  expect_true(is.data.frame(fit$model))
  expect_true("y" %in% names(fit$model))
})

test_that("T5e-02b: lm object has qr component (QR decomposition)", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$qr),
              info = "lm object must contain $qr for sandwich::bread()")
  expect_s3_class(fit$qr, "qr")
})

test_that("T5e-02c: lm object has residuals component", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$residuals),
              info = "lm object must contain $residuals for sandwich::estfun()")
  expect_true(is.numeric(fit$residuals))
  expect_equal(length(fit$residuals), 80L)
})

test_that("T5e-02d: lm object has df.residual component", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$df.residual),
              info = "lm object must contain $df.residual for VCE df computation")
  expect_true(is.numeric(fit$df.residual))
  # Tier 3 (no controls): p = 2 (intercept + D), df = N - p = 80 - 2 = 78
  expect_equal(fit$df.residual, 78L)
})

test_that("T5e-02e: lm object has all four sandwich components with controls", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 3)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")
  fit <- result$fit
  # Tier 1 full_interaction: p = 1 + 1 + K + K = 2 + 2*3 = 8
  expect_true(!is.null(fit$model))
  expect_true(!is.null(fit$qr))
  expect_true(!is.null(fit$residuals))
  expect_true(!is.null(fit$df.residual))
  expect_equal(fit$df.residual, 100L - 8L)  # N - p = 92
  expect_equal(length(fit$residuals), 100L)
})

test_that("T5e-02f: lm object coefficients are consistent with params field", {
  dat <- make_ra_data(n = 100, n_treated = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  fit_coefs <- coef(result$fit)
  expect_equal(fit_coefs, result$params)
  # D coefficient should match ATT
  expect_equal(unname(fit_coefs[["D"]]), result$att)
})

test_that("T5e-02g: lm residuals are consistent with resid field", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(as.numeric(residuals(result$fit)), result$resid)
})


# ============================================================================
# T5e-03: estimate_period_effects() returns 15 columns with correct names
# ============================================================================

test_that("T5e-03a: period results have exactly 15 columns with correct names", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:3, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL, cluster_var = NULL)
  )

  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated",
                     "n_control", "df", "vce_type", "n_clusters",
                     "controls_tier")
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 3L)
  for (col in expected_cols) {
    expect_true(col %in% names(result),
                info = sprintf("Missing column: %s", col))
  }
})

test_that("T5e-03b: period results column order matches specification", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "hc3")
  )

  expected_order <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                      "ci_lower", "ci_upper", "n_obs", "n_treated",
                      "n_control", "df", "vce_type", "n_clusters",
                      "controls_tier")
  expect_equal(names(result), expected_order)
})

test_that("T5e-03c: period results with controls still have 15 columns", {
  dt <- make_period_data(n_per_period = 80, n_treated_per = 25,
                         periods = 1:2, K = 2)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"),
                            periods = 1:2, vce = "hc3")
  )
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 2L)
})

test_that("T5e-03d: period results with cluster VCE still have 15 columns", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:3, K = 0, n_clusters = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 3L)
})

test_that("T5e-03e: period results column types are correct", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0, n_clusters = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )
  # Numeric columns
  expect_true(is.numeric(result$tindex))
  expect_true(is.numeric(result$period))
  expect_true(is.numeric(result$att))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$t_stat))
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  # Integer columns
  expect_true(is.integer(result$n_obs))
  expect_true(is.integer(result$n_treated))
  expect_true(is.integer(result$n_control))
  # Character columns
  expect_true(is.character(result$vce_type))
  expect_true(is.character(result$controls_tier))
})

test_that("T5e-03f: tindex and period columns match requested periods", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = c(2L, 5L, 8L), K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = c(2L, 5L, 8L), vce = NULL)
  )
  expect_equal(result$tindex, c(2, 5, 8))
  expect_equal(result$period, c(2, 5, 8))
})

test_that("T5e-03g: numerical reasonableness of period results", {
  set.seed(5803)
  dt <- make_period_data(n_per_period = 200, n_treated_per = 60,
                         periods = 1:2, K = 0, seed = 5803)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  for (i in seq_len(nrow(result))) {
    row <- result[i, ]
    # ATT should be in a reasonable range (true effect is 1.5)
    expect_true(abs(row$att - 1.5) < 1.0,
                info = sprintf("Period %d: ATT=%f far from true 1.5", row$period, row$att))
    # SE should be positive and reasonable
    expect_true(row$se > 0 && row$se < 2.0,
                info = sprintf("Period %d: SE=%f unreasonable", row$period, row$se))
    # CI should bracket ATT
    expect_true(row$ci_lower < row$att && row$att < row$ci_upper)
    # p-value in [0, 1]
    expect_true(row$pvalue >= 0 && row$pvalue <= 1)
    # n_obs = n_treated + n_control
    expect_equal(row$n_obs, row$n_treated + row$n_control)
    # df should be positive
    expect_true(row$df > 0)
  }
})


# ============================================================================
# T5e-04: vcov matrix dimensions match parameter count p by tier
# ============================================================================
# Three-tier degradation:
#   Tier 1 (full_interaction): p = 2 + 2K  (intercept + D + K centered + K interactions)
#   Tier 2 (simple):           p = K + 2    (intercept + D + K controls)
#   Tier 3 (none):             p = 2        (intercept + D)

test_that("T5e-04a: Tier 3 (no controls) vcov is 2x2", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(result$controls_tier, "none")
  p <- 2L
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  # Verify symmetry
  expect_equal(result$vcov, t(result$vcov))
  # Verify diagonal is non-negative (variances)
  expect_true(all(diag(result$vcov) >= 0))
})

test_that("T5e-04b: Tier 1 (full_interaction) vcov is (2+2K)x(2+2K)", {
  K <- 3L
  dat <- make_ra_data(n = 100, n_treated = 30, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K  # = 8
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  expect_equal(result$K, K)
  # Verify symmetry
  expect_true(max(abs(result$vcov - t(result$vcov))) < 1e-12)
  # Verify diagonal is non-negative
  expect_true(all(diag(result$vcov) >= 0))
})

test_that("T5e-04c: Tier 1 vcov with K=1 is 4x4", {
  K <- 1L
  dat <- make_ra_data(n = 80, n_treated = 25, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K  # = 4
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
})

test_that("T5e-04d: Tier 2 (simple) vcov is (K+2)x(K+2)", {
  # Force Tier 2: need n1 <= K+1 or n0 <= K+1, but n > K+2
  # With K=3, need n1 <= 4 or n0 <= 4, and n > 5
  # Use n=20, n_treated=4 (n1=4 <= K+1=4), n0=16 > K+1=4
  K <- 3L
  set.seed(5804)
  n <- 20L
  n1 <- 4L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  result <- suppressWarnings(
    estimate_ra_common(y, d, x = x, vce = NULL)
  )
  expect_equal(result$controls_tier, "simple")
  p <- K + 2L  # = 5
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
})

test_that("T5e-04e: Tier 3 forced by insufficient sample with controls", {
  # Force Tier 3: need n <= K+2
  # With K=5, need n <= 7
  K <- 5L
  set.seed(5805)
  n <- 7L
  n1 <- 3L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  result <- suppressWarnings(
    estimate_ra_common(y, d, x = x, vce = NULL)
  )
  expect_equal(result$controls_tier, "none")
  p <- 2L
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  # K should be reset to 0 when controls are dropped
  expect_equal(result$K, 0L)
})

test_that("T5e-04f: vcov dimensions match params length across all tiers", {
  # Tier 3
  dat0 <- make_ra_data(n = 80, n_treated = 25, K = 0)
  r0 <- estimate_ra_common(dat0$y, dat0$d, x = NULL, vce = NULL)
  expect_equal(nrow(r0$vcov), length(r0$params))
  expect_equal(ncol(r0$vcov), length(r0$params))

  # Tier 1
  dat2 <- make_ra_data(n = 100, n_treated = 30, K = 2)
  r2 <- estimate_ra_common(dat2$y, dat2$d, x = dat2$x, vce = "hc3")
  expect_equal(nrow(r2$vcov), length(r2$params))
  expect_equal(ncol(r2$vcov), length(r2$params))
})

test_that("T5e-04g: vcov dimensions match with cluster VCE", {
  K <- 2L
  set.seed(5806)
  n <- 60L
  n1 <- 20L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = x, vce = "cluster", cluster = cl)
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K  # = 6
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  expect_equal(length(result$params), p)
})

test_that("T5e-04h: X_design dimensions match vcov dimensions", {
  K <- 2L
  n <- 100L
  dat <- make_ra_data(n = n, n_treated = 30, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  p <- ncol(result$vcov)
  expect_equal(ncol(result$X_design), p)
  expect_equal(nrow(result$X_design), n)
})


# ============================================================================
# T5e-05: Non-cluster VCE in period results: n_clusters is NA_integer_
# ============================================================================
# In data.frame columns, NULL cannot be stored per-row. When VCE is not
# "cluster", the period results must use NA_integer_ (not NULL) for
# n_clusters to maintain consistent column types across rows.

test_that("T5e-05a: homoskedastic VCE period results have n_clusters = NA_integer_", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:3, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL)
  )
  expect_true(is.integer(result$n_clusters),
              info = "n_clusters column must be integer type")
  expect_true(all(is.na(result$n_clusters)),
              info = "n_clusters must be NA for non-cluster VCE")
  # Verify it's specifically NA_integer_, not NA_real_ or NA
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_,
                     info = sprintf("Row %d: n_clusters should be NA_integer_", i))
  }
})

test_that("T5e-05b: HC3 VCE period results have n_clusters = NA_integer_", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "hc3")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_)
  }
})

test_that("T5e-05c: robust VCE period results have n_clusters = NA_integer_", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "robust")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
})

test_that("T5e-05d: cluster VCE period results have integer n_clusters > 0", {
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0, n_clusters = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(!is.na(result$n_clusters)))
  expect_true(all(result$n_clusters == 10L))
})

test_that("T5e-05e: mixed success/failure periods maintain n_clusters type consistency", {
  # Period 1: normal data (should succeed)
  # Period 2: only treated units (should fail → NA row)
  set.seed(5805)
  dt1 <- data.table(tindex = 1L, y_trans = rnorm(40),
                     d = c(rep(1L, 15), rep(0L, 25)))
  dt2 <- data.table(tindex = 2L, y_trans = rnorm(10),
                     d = rep(1L, 10))  # all treated → will fail
  dt <- rbindlist(list(dt1, dt2))

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  # Column type must be consistent even when some rows are NA
  expect_true(is.integer(result$n_clusters))
  expect_equal(nrow(result), 2L)
  # Period 1: NA_integer_ (non-cluster VCE)
  expect_identical(result$n_clusters[1], NA_integer_)
  # Period 2: NA_integer_ (failed estimation)
  expect_identical(result$n_clusters[2], NA_integer_)
})


# ============================================================================
# T5e-06: vce_type field correctly reflects VCE type
# ============================================================================
# Mapping:
#   vce=NULL    → vce_type = "homoskedastic"
#   vce="hc3"   → vce_type = "HC3"
#   vce="robust" → vce_type = "HC1" (robust is alias for hc1)
#   vce="cluster"→ vce_type = "cluster"

test_that("T5e-06a: vce=NULL → vce_type='homoskedastic' in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(result$vce_type, "homoskedastic")
})

test_that("T5e-06b: vce='hc3' → vce_type='HC3' in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_equal(result$vce_type, "HC3")
})

test_that("T5e-06c: vce='robust' → vce_type='HC1' in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_equal(result$vce_type, "HC1")
})

test_that("T5e-06d: vce='cluster' → vce_type='cluster' in estimate_ra_common", {
  set.seed(5806)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL, vce = "cluster", cluster = cl)
  expect_equal(result$vce_type, "cluster")
})

test_that("T5e-06e: vce='hc0' → vce_type='HC0'", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc0")
  expect_equal(result$vce_type, "HC0")
})

test_that("T5e-06f: vce='hc1' → vce_type='HC1'", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc1")
  expect_equal(result$vce_type, "HC1")
})

test_that("T5e-06g: vce='hc2' → vce_type='HC2'", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc2")
  expect_equal(result$vce_type, "HC2")
})

test_that("T5e-06h: vce='hc4' → vce_type='HC4'", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc4")
  expect_equal(result$vce_type, "HC4")
})

test_that("T5e-06i: vce_type in period results matches VCE type", {
  # Homoskedastic
  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1:2, K = 0)
  r1 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  expect_true(all(r1$vce_type == "homoskedastic"))

  # HC3
  r2 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "hc3")
  )
  expect_true(all(r2$vce_type == "HC3"))

  # Robust → HC1
  r3 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "robust")
  )
  expect_true(all(r3$vce_type == "HC1"))

  # Cluster
  dt_cl <- make_period_data(n_per_period = 60, n_treated_per = 20,
                            periods = 1:2, K = 0, n_clusters = 10)
  r4 <- suppressWarnings(
    estimate_period_effects(dt_cl, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(all(r4$vce_type == "cluster"))
})

test_that("T5e-06j: case insensitivity - vce='HC3' and 'Hc3' both work", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  r1 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "HC3")
  r2 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "Hc3")
  expect_equal(r1$vce_type, "HC3")
  expect_equal(r2$vce_type, "HC3")
})


# ============================================================================
# T5e-07: Non-cluster VCE n_clusters is NULL in estimate_ra_common return
# ============================================================================
# In the estimate_ra_common() return value (a named list, not a data.frame),
# n_clusters should be NULL (not NA_integer_) when VCE is not cluster.
# This is distinct from period results where NA_integer_ is required for
# data.frame column type consistency.

test_that("T5e-07a: vce=NULL → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_null(result$n_clusters)
})

test_that("T5e-07b: vce='hc3' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(result$n_clusters)
})

test_that("T5e-07c: vce='robust' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_null(result$n_clusters)
})

test_that("T5e-07d: vce='hc0' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc0")
  expect_null(result$n_clusters)
})

test_that("T5e-07e: vce='hc1' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc1")
  expect_null(result$n_clusters)
})

test_that("T5e-07f: vce='hc2' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc2")
  expect_null(result$n_clusters)
})

test_that("T5e-07g: vce='hc4' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc4")
  expect_null(result$n_clusters)
})

test_that("T5e-07h: vce='cluster' → n_clusters is integer > 0 in estimate_ra_common", {
  set.seed(5807)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL, vce = "cluster", cluster = cl)
  expect_false(is.null(result$n_clusters))
  expect_true(is.integer(result$n_clusters) || is.numeric(result$n_clusters))
  expect_equal(result$n_clusters, 10L)
})

test_that("T5e-07i: NULL vs NA_integer_ distinction between list and data.frame", {
  # This test verifies the design contract:
  # estimate_ra_common() returns NULL for n_clusters (list context)
  # estimate_period_effects() returns NA_integer_ for n_clusters (data.frame context)
  dat <- make_ra_data(n = 80, n_treated = 25, K = 0)
  ra_result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(ra_result$n_clusters)

  dt <- make_period_data(n_per_period = 60, n_treated_per = 20,
                         periods = 1L, K = 0)
  pe_result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1L, vce = "hc3")
  )
  expect_identical(pe_result$n_clusters[1], NA_integer_)
})
