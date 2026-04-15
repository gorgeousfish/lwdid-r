# test-return-completeness.R - E3-05.8 返回值完整性测试
# ============================================================================
# Return Value Completeness Tests for:
#   - estimate_ra_common():      Named list with lm object + VCE results
#   - estimate_period_effects(): 15-column data.frame per period
#
# Test blocks:
#   T8-01: estimate_ra_common() returns complete lm object
#   T8-02: estimate_period_effects() returns 15-column data.frame
#   T8-03: vcov matrix dimensions match model tier
#   T8-04: n_clusters handling (NULL vs NA_integer_)
#   T8-05: vce_type field correctness
#   T8-06: lm object has sandwich-compatible components
# ============================================================================
library(testthat)
library(data.table)

# ============================================================================
# Helpers
# ============================================================================

#' Generate simple cross-sectional data for estimate_ra_common()
make_cs_data <- function(n = 100, n1 = 30, K = 0, seed = 8000) {
  set.seed(seed)
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.5 * d
  x <- if (K > 0) matrix(rnorm(n * K), nrow = n, ncol = K) else NULL
  list(y = y, d = d, x = x)
}

#' Generate panel data for estimate_period_effects()
make_panel_data <- function(n_per = 60, n1_per = 20, periods = 1:3,
                            K = 0, n_cl = NULL, seed = 8001) {
  set.seed(seed)
  dt_list <- lapply(periods, function(r) {
    d <- c(rep(1L, n1_per), rep(0L, n_per - n1_per))
    y <- rnorm(n_per) + 1.8 * d
    row <- data.table(tindex = r, y_trans = y, d = d)
    if (K > 0) {
      for (k in seq_len(K)) {
        set(row, j = paste0("x", k), value = rnorm(n_per))
      }
    }
    if (!is.null(n_cl)) {
      row[, cluster := rep(seq_len(n_cl), length.out = n_per)]
    }
    row
  })
  rbindlist(dt_list)
}


# ============================================================================
# T8-01: estimate_ra_common() returns complete lm object
# ============================================================================

test_that("T8-01a: fit is lm class — no controls (Tier 3)", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_s3_class(result$fit, "lm")
  expect_equal(result$controls_tier, "none")
})

test_that("T8-01b: fit is lm class — with controls (Tier 1)", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  expect_s3_class(result$fit, "lm")
  expect_equal(result$controls_tier, "full_interaction")
})

test_that("T8-01c: fit is lm class — HC3 VCE", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_s3_class(result$fit, "lm")
})

test_that("T8-01d: fit is lm class — cluster VCE", {
  set.seed(8010)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL,
                               vce = "cluster", cluster = cl)
  expect_s3_class(result$fit, "lm")
})

test_that("T8-01e: lm fit has model, qr, residuals, df.residual", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$model))
  expect_true(!is.null(fit$qr))
  expect_true(!is.null(fit$residuals))
  expect_true(!is.null(fit$df.residual))
})

test_that("T8-01f: return list has all expected fields", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")

  expected <- c(
    "att", "se", "t_stat", "df", "pvalue", "ci_lower", "ci_upper",
    "params", "vcov", "resid", "fit",
    "n", "n_treated", "n_control", "K",
    "vce_type", "n_clusters", "controls_tier", "X_design",
    "df_resid", "df_inference"
  )
  expect_true(is.list(result))
  for (nm in expected) {
    expect_true(nm %in% names(result),
                info = sprintf("Missing field: %s", nm))
  }
  expect_equal(length(result), length(expected))
})


# ============================================================================
# T8-02: estimate_period_effects() returns 15-column data.frame
# ============================================================================

test_that("T8-02a: 15 columns with correct names — no controls", {
  dt <- make_panel_data(n_per = 60, n1_per = 20, periods = 1:3, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:3, vce = NULL)
  )
  expected_cols <- c(
    "tindex", "period", "att", "se", "t_stat", "pvalue",
    "ci_lower", "ci_upper", "n_obs", "n_treated",
    "n_control", "df", "vce_type", "n_clusters", "controls_tier"
  )
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 3L)
  for (col in expected_cols) {
    expect_true(col %in% names(result),
                info = sprintf("Missing column: %s", col))
  }
})

test_that("T8-02b: column order matches specification", {
  dt <- make_panel_data(n_per = 60, n1_per = 20, periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "hc3")
  )
  expected_order <- c(
    "tindex", "period", "att", "se", "t_stat", "pvalue",
    "ci_lower", "ci_upper", "n_obs", "n_treated",
    "n_control", "df", "vce_type", "n_clusters", "controls_tier"
  )
  expect_equal(names(result), expected_order)
})

test_that("T8-02c: 15 columns with controls (Tier 1)", {
  dt <- make_panel_data(n_per = 80, n1_per = 25,
                        periods = 1:2, K = 2)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"),
                            periods = 1:2, vce = "hc3")
  )
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 2L)
})

test_that("T8-02d: 15 columns with cluster VCE", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:3, K = 0, n_cl = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:3,
                            vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_equal(ncol(result), 15L)
  expect_equal(nrow(result), 3L)
})

test_that("T8-02e: column types are correct", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0, n_cl = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2,
                            vce = "cluster",
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

test_that("T8-02f: tindex and period match requested periods", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = c(3L, 7L, 11L), K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL,
                            periods = c(3L, 7L, 11L), vce = NULL)
  )
  expect_equal(result$tindex, c(3, 7, 11))
  expect_equal(result$period, c(3, 7, 11))
})

test_that("T8-02g: numerical reasonableness of period ATT estimates", {
  dt <- make_panel_data(n_per = 200, n1_per = 60,
                        periods = 1:2, K = 0, seed = 8020)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = NULL)
  )
  for (i in seq_len(nrow(result))) {
    row <- result[i, ]
    # True effect is 1.8; ATT should be within 1.0 of it
    expect_true(abs(row$att - 1.8) < 1.0,
                info = sprintf("Period %d: ATT=%.3f far from 1.8",
                               row$period, row$att))
    # SE positive and reasonable
    expect_true(row$se > 0 && row$se < 2.0,
                info = sprintf("Period %d: SE=%.3f unreasonable",
                               row$period, row$se))
    # CI brackets ATT
    expect_true(row$ci_lower < row$att)
    expect_true(row$att < row$ci_upper)
    # p-value in [0, 1]
    expect_true(row$pvalue >= 0 && row$pvalue <= 1)
    # n_obs = n_treated + n_control
    expect_equal(row$n_obs, row$n_treated + row$n_control)
    # df positive
    expect_true(row$df > 0)
  }
})


# ============================================================================
# T8-03: vcov matrix dimensions match model tier
# ============================================================================
# Three-tier degradation:
#   Tier 1 (full_interaction): p = 2 + 2K
#   Tier 2 (simple):           p = K + 2
#   Tier 3 (none):             p = 2

test_that("T8-03a: Tier 3 (no controls) — vcov is 2x2", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(result$controls_tier, "none")
  expect_equal(nrow(result$vcov), 2L)
  expect_equal(ncol(result$vcov), 2L)
  # Symmetry
  expect_equal(result$vcov, t(result$vcov))
  # Non-negative diagonal (variances)
  expect_true(all(diag(result$vcov) >= 0))
})

test_that("T8-03b: Tier 1 (full_interaction) K=3 — vcov is 8x8", {
  K <- 3L
  dat <- make_cs_data(n = 100, n1 = 30, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  expect_equal(result$K, K)
  # Near-symmetry (numerical tolerance)
  expect_true(max(abs(result$vcov - t(result$vcov))) < 1e-12)
  expect_true(all(diag(result$vcov) >= 0))
})

test_that("T8-03c: Tier 1 K=1 — vcov is 4x4", {
  K <- 1L
  dat <- make_cs_data(n = 80, n1 = 25, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
})

test_that("T8-03d: Tier 2 (simple) — vcov is (K+2)x(K+2)", {

  # Force Tier 2: n1 <= K+1 but n > K+2
  # K=3 → need n1 <= 4, n > 5
  K <- 3L
  set.seed(8030)
  n <- 20L
  n1 <- 4L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  result <- suppressWarnings(
    estimate_ra_common(y, d, x = x, vce = NULL)
  )
  expect_equal(result$controls_tier, "simple")
  p <- K + 2L
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
})

test_that("T8-03e: Tier 3 forced by insufficient sample with controls", {
  # Force Tier 3: n <= K+2
  # K=5 → need n <= 7
  K <- 5L
  set.seed(8031)
  n <- 7L
  n1 <- 3L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  result <- suppressWarnings(
    estimate_ra_common(y, d, x = x, vce = NULL)
  )
  expect_equal(result$controls_tier, "none")
  expect_equal(nrow(result$vcov), 2L)
  expect_equal(ncol(result$vcov), 2L)
  # K reset to 0 when controls dropped
  expect_equal(result$K, 0L)
})

test_that("T8-03f: vcov dims match params length across all tiers", {
  # Tier 3
  dat0 <- make_cs_data(n = 80, n1 = 25, K = 0)
  r0 <- estimate_ra_common(dat0$y, dat0$d, x = NULL, vce = NULL)
  expect_equal(nrow(r0$vcov), length(r0$params))
  expect_equal(ncol(r0$vcov), length(r0$params))

  # Tier 1
  dat2 <- make_cs_data(n = 100, n1 = 30, K = 2)
  r2 <- estimate_ra_common(dat2$y, dat2$d, x = dat2$x, vce = "hc3")
  expect_equal(nrow(r2$vcov), length(r2$params))
  expect_equal(ncol(r2$vcov), length(r2$params))
})

test_that("T8-03g: vcov dims match with cluster VCE", {
  K <- 2L
  set.seed(8032)
  n <- 60L
  n1 <- 20L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = x,
                               vce = "cluster", cluster = cl)
  expect_equal(result$controls_tier, "full_interaction")
  p <- 2L + 2L * K
  expect_equal(nrow(result$vcov), p)
  expect_equal(ncol(result$vcov), p)
  expect_equal(length(result$params), p)
})

test_that("T8-03h: X_design dims match vcov dims", {
  K <- 2L
  n <- 100L
  dat <- make_cs_data(n = n, n1 = 30, K = K)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  p <- ncol(result$vcov)
  expect_equal(ncol(result$X_design), p)
  expect_equal(nrow(result$X_design), n)
})


# ============================================================================
# T8-04: n_clusters handling — NULL vs NA_integer_
# ============================================================================
# Design contract:
#   estimate_ra_common() → n_clusters is NULL for non-cluster VCE (list)
#   estimate_period_effects() → n_clusters is NA_integer_ for non-cluster
#     VCE (data.frame column type consistency)

test_that("T8-04a: vce=NULL → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_null(result$n_clusters)
})

test_that("T8-04b: vce='hc3' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(result$n_clusters)
})

test_that("T8-04c: vce='robust' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_null(result$n_clusters)
})

test_that("T8-04d: vce='cluster' → n_clusters is integer > 0", {
  set.seed(8040)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL,
                               vce = "cluster", cluster = cl)
  expect_false(is.null(result$n_clusters))
  expect_true(is.integer(result$n_clusters) ||
              is.numeric(result$n_clusters))
  expect_equal(result$n_clusters, 10L)
})

test_that("T8-04e: homoskedastic period results → n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:3, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:3, vce = NULL)
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_)
  }
})

test_that("T8-04f: HC3 period results → n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "hc3")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_)
  }
})

test_that("T8-04g: robust period results → n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "robust")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
})

test_that("T8-04h: cluster period results → integer n_clusters > 0", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0, n_cl = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2,
                            vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(!is.na(result$n_clusters)))
  expect_true(all(result$n_clusters == 10L))
})

test_that("T8-04i: NULL vs NA_integer_ distinction — list vs data.frame", {
  # estimate_ra_common returns NULL (list context)
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  ra <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(ra$n_clusters)

  # estimate_period_effects returns NA_integer_ (data.frame context)
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1L, K = 0)
  pe <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1L, vce = "hc3")
  )
  expect_identical(pe$n_clusters[1], NA_integer_)
})

test_that("T8-04j: mixed success/failure periods maintain n_clusters type", {
  # Period 1: normal data (succeeds)
  # Period 2: all treated (fails → NA row)
  set.seed(8041)
  dt1 <- data.table(tindex = 1L, y_trans = rnorm(40),
                    d = c(rep(1L, 15), rep(0L, 25)))
  dt2 <- data.table(tindex = 2L, y_trans = rnorm(10),
                    d = rep(1L, 10))
  dt <- rbindlist(list(dt1, dt2))

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = NULL)
  )
  expect_true(is.integer(result$n_clusters))
  expect_equal(nrow(result), 2L)
  expect_identical(result$n_clusters[1], NA_integer_)
  expect_identical(result$n_clusters[2], NA_integer_)
})


# ============================================================================
# T8-05: vce_type field correctness
# ============================================================================
# Mapping (from compute_vce in inference.R):
#   vce=NULL     → "homoskedastic"
#   vce="hc0"    → "HC0"   (toupper)
#   vce="hc1"    → "HC1"   (toupper)
#   vce="hc2"    → "HC2"   (toupper)
#   vce="hc3"    → "HC3"   (toupper)
#   vce="hc4"    → "HC4"   (toupper)
#   vce="robust" → "HC1"   (alias for hc1, then toupper)
#   vce="cluster"→ "cluster"

test_that("T8-05a: vce=NULL → 'homoskedastic'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(result$vce_type, "homoskedastic")
})

test_that("T8-05b: vce='hc3' → 'HC3'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_equal(result$vce_type, "HC3")
})

test_that("T8-05c: vce='robust' → 'HC1'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_equal(result$vce_type, "HC1")
})

test_that("T8-05d: vce='cluster' → 'cluster'", {
  set.seed(8050)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL,
                               vce = "cluster", cluster = cl)
  expect_equal(result$vce_type, "cluster")
})

test_that("T8-05e: vce='hc0' → 'HC0'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc0")
  expect_equal(result$vce_type, "HC0")
})

test_that("T8-05f: vce='hc1' → 'HC1'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc1")
  expect_equal(result$vce_type, "HC1")
})

test_that("T8-05g: vce='hc2' → 'HC2'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc2")
  expect_equal(result$vce_type, "HC2")
})

test_that("T8-05h: vce='hc4' → 'HC4'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc4")
  expect_equal(result$vce_type, "HC4")
})

test_that("T8-05i: case insensitivity — 'HC3' and 'Hc3' both → 'HC3'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  r1 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "HC3")
  r2 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "Hc3")
  expect_equal(r1$vce_type, "HC3")
  expect_equal(r2$vce_type, "HC3")
})

test_that("T8-05j: vce_type in period results matches VCE type", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)

  # Homoskedastic
  r1 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = NULL)
  )
  expect_true(all(r1$vce_type == "homoskedastic"))

  # HC3
  r2 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "hc3")
  )
  expect_true(all(r2$vce_type == "HC3"))

  # Robust → HC1
  r3 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "robust")
  )
  expect_true(all(r3$vce_type == "HC1"))

  # Cluster
  dt_cl <- make_panel_data(n_per = 60, n1_per = 20,
                           periods = 1:2, K = 0, n_cl = 10)
  r4 <- suppressWarnings(
    estimate_period_effects(dt_cl, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2,
                            vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(all(r4$vce_type == "cluster"))
})


# ============================================================================
# T8-04: n_clusters handling (NULL vs NA_integer_)
# ============================================================================
# Design contract:
#   estimate_ra_common() returns NULL for n_clusters (list context)
#   estimate_period_effects() returns NA_integer_ for n_clusters (data.frame)

test_that("T8-04a: vce=NULL → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_null(result$n_clusters)
})

test_that("T8-04b: vce='hc3' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(result$n_clusters)
})

test_that("T8-04c: vce='robust' → n_clusters is NULL in estimate_ra_common", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_null(result$n_clusters)
})

test_that("T8-04d: vce='cluster' → n_clusters is integer > 0", {
  set.seed(8040)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL,
                               vce = "cluster", cluster = cl)
  expect_false(is.null(result$n_clusters))
  expect_equal(result$n_clusters, 10L)
})

test_that("T8-04e: homoskedastic period results have n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:3, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:3, vce = NULL)
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_)
  }
})

test_that("T8-04f: HC3 period results have n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "hc3")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
  for (i in seq_len(nrow(result))) {
    expect_identical(result$n_clusters[i], NA_integer_)
  }
})

test_that("T8-04g: robust period results have n_clusters = NA_integer_", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "robust")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(is.na(result$n_clusters)))
})

test_that("T8-04h: cluster period results have integer n_clusters > 0", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0, n_cl = 10)
  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2,
                            vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(is.integer(result$n_clusters))
  expect_true(all(!is.na(result$n_clusters)))
  expect_true(all(result$n_clusters == 10L))
})

test_that("T8-04i: NULL vs NA_integer_ distinction between list and df", {
  # estimate_ra_common() → NULL (list context)
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  ra_result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_null(ra_result$n_clusters)

  # estimate_period_effects() → NA_integer_ (data.frame context)
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1L, K = 0)
  pe_result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1L, vce = "hc3")
  )
  expect_identical(pe_result$n_clusters[1], NA_integer_)
})

test_that("T8-04j: mixed success/failure periods maintain n_clusters type", {
  # Period 1: normal, Period 2: all treated → fails
  set.seed(8041)
  dt1 <- data.table(tindex = 1L, y_trans = rnorm(40),
                    d = c(rep(1L, 15), rep(0L, 25)))
  dt2 <- data.table(tindex = 2L, y_trans = rnorm(10),
                    d = rep(1L, 10))
  dt <- rbindlist(list(dt1, dt2))

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = NULL)
  )
  expect_true(is.integer(result$n_clusters))
  expect_equal(nrow(result), 2L)
  expect_identical(result$n_clusters[1], NA_integer_)
  expect_identical(result$n_clusters[2], NA_integer_)
})


# ============================================================================
# T8-05: vce_type field correctness
# ============================================================================
# Mapping (from compute_vce):
#   vce=NULL     → "homoskedastic"
#   vce="hc0"    → "HC0"
#   vce="hc1"    → "HC1"
#   vce="hc2"    → "HC2"
#   vce="hc3"    → "HC3"
#   vce="hc4"    → "HC4"
#   vce="robust" → "HC1" (alias for hc1)
#   vce="cluster"→ "cluster"

test_that("T8-05a: vce=NULL → 'homoskedastic'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(result$vce_type, "homoskedastic")
})

test_that("T8-05b: vce='hc3' → 'HC3'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc3")
  expect_equal(result$vce_type, "HC3")
})

test_that("T8-05c: vce='robust' → 'HC1'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "robust")
  expect_equal(result$vce_type, "HC1")
})

test_that("T8-05d: vce='cluster' → 'cluster'", {
  set.seed(8050)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL,
                               vce = "cluster", cluster = cl)
  expect_equal(result$vce_type, "cluster")
})

test_that("T8-05e: vce='hc0' → 'HC0'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc0")
  expect_equal(result$vce_type, "HC0")
})

test_that("T8-05f: vce='hc1' → 'HC1'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc1")
  expect_equal(result$vce_type, "HC1")
})

test_that("T8-05g: vce='hc2' → 'HC2'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc2")
  expect_equal(result$vce_type, "HC2")
})

test_that("T8-05h: vce='hc4' → 'HC4'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "hc4")
  expect_equal(result$vce_type, "HC4")
})

test_that("T8-05i: case insensitivity — 'HC3' and 'Hc3' both → 'HC3'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  r1 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "HC3")
  r2 <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = "Hc3")
  expect_equal(r1$vce_type, "HC3")
  expect_equal(r2$vce_type, "HC3")
})

test_that("T8-05j: vce_type in period results matches VCE type", {
  dt <- make_panel_data(n_per = 60, n1_per = 20,
                        periods = 1:2, K = 0)

  # Homoskedastic
  r1 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = NULL)
  )
  expect_true(all(r1$vce_type == "homoskedastic"))

  # HC3
  r2 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "hc3")
  )
  expect_true(all(r2$vce_type == "HC3"))

  # Robust → HC1
  r3 <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2, vce = "robust")
  )
  expect_true(all(r3$vce_type == "HC1"))

  # Cluster
  dt_cl <- make_panel_data(n_per = 60, n1_per = 20,
                           periods = 1:2, K = 0, n_cl = 10)
  r4 <- suppressWarnings(
    estimate_period_effects(dt_cl, "y_trans", "d", "tindex",
                            x = NULL, periods = 1:2,
                            vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_true(all(r4$vce_type == "cluster"))
})


# ============================================================================
# T8-06: lm object has sandwich-compatible components
# ============================================================================

test_that("T8-06a: model frame is a data.frame with 'y' column", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(is.data.frame(fit$model))
  expect_true("y" %in% names(fit$model))
})

test_that("T8-06b: qr component is of class 'qr'", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_s3_class(result$fit$qr, "qr")
})

test_that("T8-06c: residuals length matches sample size", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_true(is.numeric(result$fit$residuals))
  expect_equal(length(result$fit$residuals), 80L)
})

test_that("T8-06d: df.residual = N - p for Tier 3", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  # Tier 3: p = 2, df = 80 - 2 = 78
  expect_equal(result$fit$df.residual, 78L)
})

test_that("T8-06e: all four sandwich components present with controls", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 3)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")
  fit <- result$fit
  # Tier 1: p = 2 + 2*3 = 8
  expect_true(!is.null(fit$model))
  expect_true(!is.null(fit$qr))
  expect_true(!is.null(fit$residuals))
  expect_true(!is.null(fit$df.residual))
  expect_equal(fit$df.residual, 100L - 8L)
  expect_equal(length(fit$residuals), 100L)
})

test_that("T8-06f: lm coefficients match params field", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  expect_equal(coef(result$fit), result$params)
  expect_equal(unname(coef(result$fit)[["D"]]), result$att)
})

test_that("T8-06g: lm residuals match resid field", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(as.numeric(residuals(result$fit)), result$resid)
})


# ============================================================================
# T8-06: lm object has sandwich-compatible components
# ============================================================================
# sandwich::vcovHC and sandwich::vcovCL require:
#   - model:       model frame (data.frame used in fitting)
#   - qr:          QR decomposition of design matrix
#   - residuals:   OLS residuals vector
#   - df.residual: residual degrees of freedom (N - p)

test_that("T8-06a: model component is a data.frame with 'y' column", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$model))
  expect_true(is.data.frame(fit$model))
  expect_true("y" %in% names(fit$model))
})

test_that("T8-06b: qr component is a qr object", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$qr))
  expect_s3_class(fit$qr, "qr")
})

test_that("T8-06c: residuals are numeric with correct length", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$residuals))
  expect_true(is.numeric(fit$residuals))
  expect_equal(length(fit$residuals), 80L)
})

test_that("T8-06d: df.residual is correct for Tier 3 (p=2)", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  fit <- result$fit
  expect_true(!is.null(fit$df.residual))
  expect_true(is.numeric(fit$df.residual))
  # Tier 3: p = 2, df = N - p = 80 - 2 = 78
  expect_equal(fit$df.residual, 78L)
})

test_that("T8-06e: all four components present with Tier 1 controls", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 3)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = "hc3")
  fit <- result$fit
  # Tier 1: p = 2 + 2*3 = 8
  expect_true(!is.null(fit$model))
  expect_true(!is.null(fit$qr))
  expect_true(!is.null(fit$residuals))
  expect_true(!is.null(fit$df.residual))
  expect_equal(fit$df.residual, 100L - 8L)
  expect_equal(length(fit$residuals), 100L)
})

test_that("T8-06f: lm coefficients consistent with params field", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  fit_coefs <- coef(result$fit)
  expect_equal(fit_coefs, result$params)
  # D coefficient = ATT
  expect_equal(unname(fit_coefs[["D"]]), result$att)
})

test_that("T8-06g: lm residuals consistent with resid field", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  expect_equal(as.numeric(residuals(result$fit)), result$resid)
})

test_that("T8-06h: sandwich::vcovHC works on returned lm object", {
  dat <- make_cs_data(n = 80, n1 = 25, K = 0)
  result <- estimate_ra_common(dat$y, dat$d, x = NULL, vce = NULL)
  # If the lm object is sandwich-compatible, vcovHC should not error
  hc_vcov <- sandwich::vcovHC(result$fit, type = "HC3")
  expect_true(is.matrix(hc_vcov))
  expect_equal(nrow(hc_vcov), 2L)
  expect_equal(ncol(hc_vcov), 2L)
})

test_that("T8-06i: sandwich::vcovHC works with Tier 1 controls", {
  dat <- make_cs_data(n = 100, n1 = 30, K = 2)
  result <- estimate_ra_common(dat$y, dat$d, x = dat$x, vce = NULL)
  hc_vcov <- sandwich::vcovHC(result$fit, type = "HC3")
  p <- 2L + 2L * 2L  # = 6
  expect_true(is.matrix(hc_vcov))
  expect_equal(nrow(hc_vcov), p)
  expect_equal(ncol(hc_vcov), p)
})

test_that("T8-06j: sandwich::vcovCL works on returned lm object", {
  set.seed(8060)
  n <- 60L
  d <- c(rep(1L, 20), rep(0L, 40))
  y <- rnorm(n) + 2.0 * d
  cl <- rep(1:10, each = 6)
  result <- estimate_ra_common(y, d, x = NULL, vce = NULL)
  # vcovCL should work on the lm object
  cl_vcov <- sandwich::vcovCL(result$fit, cluster = cl, type = "HC1")
  expect_true(is.matrix(cl_vcov))
  expect_equal(nrow(cl_vcov), 2L)
  expect_equal(ncol(cl_vcov), 2L)
})

test_that("T8-06k: df.residual correct for Tier 2", {
  # Force Tier 2: K=3, n1=4 (n1 <= K+1=4), n=20 > K+2=5
  K <- 3L
  set.seed(8061)
  n <- 20L
  n1 <- 4L
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- rnorm(n) + 2.0 * d
  x <- matrix(rnorm(n * K), nrow = n, ncol = K)
  result <- suppressWarnings(
    estimate_ra_common(y, d, x = x, vce = NULL)
  )
  expect_equal(result$controls_tier, "simple")
  # Tier 2: p = K + 2 = 5, df = 20 - 5 = 15
  expect_equal(result$fit$df.residual, n - (K + 2L))
})
