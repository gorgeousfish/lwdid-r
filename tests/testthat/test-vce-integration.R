# test-vce-integration.R - E3-05.4 VCE Type Integration Tests
library(testthat)
library(sandwich)

capture_with_warnings <- function(expr) {
  warnings_caught <- list()
  result <- withCallingHandlers(expr, warning = function(w) {
    warnings_caught[[length(warnings_caught) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(result = result, warnings = warnings_caught)
}

has_warning_class <- function(wlist, cls) {
  any(vapply(wlist, function(w) inherits(w, cls), logical(1)))
}

make_simple <- function(seed = 42L, n = 50L) {
  set.seed(seed)
  n1 <- as.integer(n * 0.4)
  d <- c(rep(1L, n1), rep(0L, n - n1))
  y <- 1 + 2 * d + rnorm(n, sd = 0.5)
  list(y = y, d = d, n = n)
}

make_ctrl <- function(seed = 123L, n = 40L, k = 2L) {
  set.seed(seed)
  x <- matrix(rnorm(n * k), n, k)
  d <- c(rep(1L, 20L), rep(0L, n - 20L))
  y <- 1 + 2 * d + 0.5 * x[, 1] - 0.3 * x[, 2] + rnorm(n, sd = 0.3)
  list(y = y, d = d, x = x, n = n, k = k)
}

make_clust <- function(seed = 42L, g = 20L, per = 3L) {
  set.seed(seed)
  n <- g * per
  cl <- rep(seq_len(g), each = per)
  d <- rep(c(1L, 0L, 0L), g)
  y <- 1 + 2 * d + rnorm(n, sd = 0.5)
  list(y = y, d = d, cl = cl, n = n, g = g)
}

make_perfect <- function(n = 20L) {
  n1 <- as.integer(n / 2)
  d <- c(rep(1L, n1), rep(0L, n - n1))
  list(y = 1 + 2 * d, d = d, n = n)
}


# Spec 1: vce=NULL homoskedastic SE
test_that("T5-01a: vce=NULL SE matches lm() vcov", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = NULL)
  fit <- lm(s$y ~ s$d)
  expect_equal(res$se, sqrt(vcov(fit)[2, 2]), tolerance = 1e-12)
})
test_that("T5-01b: vce=NULL vce_type is homoskedastic", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = NULL)
  expect_identical(res$vce_type, "homoskedastic")
})
test_that("T5-01c: vce=NULL n_clusters is NULL", {
  s <- make_simple()
  expect_null(estimate_ra_common(s$y, s$d, vce = NULL)$n_clusters)
})
test_that("T5-01d: vce=NULL inference fields finite", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = NULL)
  expect_true(is.finite(res$att))
  expect_true(is.finite(res$se))
  expect_true(is.finite(res$t_stat))
  expect_true(is.finite(res$pvalue))
  expect_true(is.finite(res$ci_lower))
  expect_true(is.finite(res$ci_upper))
})
test_that("T5-01e: vce=NULL SE > 0 and CI contains ATT", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = NULL)
  expect_true(res$se > 0)
  expect_true(res$ci_lower <= res$att)
  expect_true(res$ci_upper >= res$att)
})


# Spec 2: vce="hc3" HC3 SE correct
test_that("T5-02a: vce=hc3 SE matches sandwich HC3", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = "hc3")
  fit <- lm(s$y ~ s$d)
  expect_equal(res$se,
    sqrt(sandwich::vcovHC(fit, type = "HC3")[2, 2]),
    tolerance = 1e-12)
})
test_that("T5-02b: vce=hc3 vce_type is HC3", {
  s <- make_simple()
  expect_identical(
    estimate_ra_common(s$y, s$d, vce = "hc3")$vce_type, "HC3")
})
test_that("T5-02c: HC3 SE differs from homoskedastic", {
  s <- make_simple()
  r3 <- estimate_ra_common(s$y, s$d, vce = "hc3")
  rh <- estimate_ra_common(s$y, s$d, vce = NULL)
  expect_false(isTRUE(all.equal(r3$se, rh$se, tolerance = 1e-10)))
})
test_that("T5-02d: HC3 full vcov matches sandwich", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = "hc3")
  fit <- lm(y ~ D, data = data.frame(y = s$y, D = s$d))
  expect_equal(res$vcov,
    sandwich::vcovHC(fit, type = "HC3"), tolerance = 1e-12)
})


# Spec 3: vce="robust" -> HC1 SE, vce_type="HC1"
test_that("T5-03a: vce=robust SE matches sandwich HC1", {
  s <- make_simple()
  res <- estimate_ra_common(s$y, s$d, vce = "robust")
  fit <- lm(s$y ~ s$d)
  expect_equal(res$se,
    sqrt(sandwich::vcovHC(fit, type = "HC1")[2, 2]),
    tolerance = 1e-12)
})
test_that("T5-03b: vce=robust vce_type is HC1 not robust", {
  s <- make_simple()
  expect_identical(
    estimate_ra_common(s$y, s$d, vce = "robust")$vce_type, "HC1")
})

# Spec 4: vce="cluster" cluster SE, n_clusters correct
test_that("T5-04a: cluster SE matches sandwich vcovCL", {
  s <- make_clust()
  res <- suppressWarnings(
    estimate_ra_common(s$y, s$d, vce = "cluster", cluster = s$cl))
  fit <- lm(s$y ~ s$d)
  se_ref <- sqrt(sandwich::vcovCL(
    fit, cluster = s$cl, type = "HC1")[2, 2])
  expect_equal(res$se, se_ref, tolerance = 1e-12)
})
test_that("T5-04b: cluster n_clusters correct", {
  s <- make_clust()
  res <- suppressWarnings(
    estimate_ra_common(s$y, s$d, vce = "cluster", cluster = s$cl))
  expect_identical(res$n_clusters, s$g)
})
test_that("T5-04c: cluster df = G - 1", {
  s <- make_clust()
  res <- suppressWarnings(
    estimate_ra_common(s$y, s$d, vce = "cluster", cluster = s$cl))
  expect_equal(res$df, s$g - 1L)
})
test_that("T5-04d: cluster vce_type is cluster", {
  s <- make_clust()
  res <- suppressWarnings(
    estimate_ra_common(s$y, s$d, vce = "cluster", cluster = s$cl))
  expect_identical(res$vce_type, "cluster")
})
test_that("T5-04e: cluster SE > 0", {
  s <- make_clust()
  res <- suppressWarnings(
    estimate_ra_common(s$y, s$d, vce = "cluster", cluster = s$cl))
  expect_true(res$se > 0)
})
