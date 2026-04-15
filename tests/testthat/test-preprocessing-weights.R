# test-preprocessing-weights.R — Tests for Story 11.2 computation functions

# ============================================================
# .kahan_sum() tests
# ============================================================

test_that("Kahan sum correct for simple input", {
  expect_equal(.kahan_sum(c(1, 2, 3)), 6.0)
  expect_equal(.kahan_sum(c(0.1, 0.2, 0.3)), 0.6, tolerance = 1e-15)
})

test_that("Kahan sum returns 0 for empty vector", {
  expect_equal(.kahan_sum(numeric(0)), 0.0)
})

test_that("Kahan sum handles alternating large values", {
  # Note: c(1e16, 1, -1e16) loses precision in individual addition
  # (1e16 + 1 == 1e16 in double). Use a case where Kahan helps.
  x <- c(1, 1e-16, -1e-16)
  expect_equal(.kahan_sum(x), 1.0)
})

test_that("Kahan sum large-scale precision", {
  n <- 1e5
  delta <- 0.5
  x <- rep(c(1e15 + delta, -(1e15)), n / 2)
  analytical_true <- (n / 2) * delta

  kahan_result <- .kahan_sum(x)
  naive_result <- sum(x)

  kahan_err <- abs(kahan_result - analytical_true)
  naive_err <- abs(naive_result - analytical_true)

  expect_equal(kahan_result, analytical_true, tolerance = 1e-6)
  expect_true(kahan_err <= naive_err + .Machine$double.eps)
})

# ============================================================
# .normalize_weights() tests
# ============================================================

test_that("weights normalize to sum 1", {
  w <- c(1, 2, 3, 4)
  nw <- .normalize_weights(w)
  expect_equal(sum(nw), 1.0, tolerance = 1e-15)
  expect_equal(nw, c(0.1, 0.2, 0.3, 0.4))
})

test_that("all-zero weights return equal weights", {
  w <- c(0, 0, 0, 0)
  nw <- .normalize_weights(w)
  expect_equal(nw, rep(0.25, 4))
  expect_equal(sum(nw), 1.0)
})

test_that("single element normalizes to 1", {
  nw <- .normalize_weights(c(5.0))
  expect_equal(nw, 1.0)
})

test_that("near-zero weight sum returns equal weights", {
  w <- c(1e-16, 1e-16, 1e-16)
  nw <- .normalize_weights(w)
  expect_equal(nw, rep(1 / 3, 3), tolerance = 1e-15)
})

test_that("empty vector returns empty vector", {
  nw <- .normalize_weights(numeric(0))
  expect_identical(nw, numeric(0))
})

test_that("large-scale normalization precision", {
  set.seed(123)
  n <- 1e5
  raw_w <- runif(n, 0.01, 100)
  nw <- .normalize_weights(raw_w)
  expect_equal(.kahan_sum(nw), 1.0, tolerance = 1e-12)
  expect_true(all(nw > 0))
})

# ============================================================
# WEIGHT_SUM_TOLERANCE constant test
# ============================================================

test_that("WEIGHT_SUM_TOLERANCE is 1e-9", {
  expect_equal(WEIGHT_SUM_TOLERANCE, 1e-9)
})

# ============================================================
# .compute_cell_weighted_mean() tests
# ============================================================

test_that("equal weights return arithmetic mean", {
  outcomes <- c(10, 20, 30)
  result <- .compute_cell_weighted_mean(outcomes)
  expect_equal(result$mean, 20.0)
  expect_equal(result$n_obs, 3L)
  expect_null(result$variance)
  expect_null(result$ess)
})

test_that("survey weights return correct weighted mean", {
  outcomes <- c(10, 20, 30)
  weights <- c(1, 2, 1)
  result <- .compute_cell_weighted_mean(outcomes, weights)
  expect_equal(result$mean, 20.0, tolerance = 1e-12)
})

test_that("unequal weights produce correct weighted mean", {
  outcomes <- c(100, 200)
  weights <- c(3, 1)
  result <- .compute_cell_weighted_mean(outcomes, weights)
  expect_equal(result$mean, 125.0, tolerance = 1e-12)
})

test_that("all-NA input returns NaN and n_obs=0", {
  outcomes <- c(NA, NA, NA)
  result <- .compute_cell_weighted_mean(outcomes)
  expect_true(is.nan(result$mean))
  expect_equal(result$n_obs, 0L)
})

test_that("partial NA correctly excluded", {
  outcomes <- c(10, NA, 30)
  result <- .compute_cell_weighted_mean(outcomes)
  expect_equal(result$mean, 20.0)
  expect_equal(result$n_obs, 2L)
})

test_that("partial NA with weights correctly aligned", {
  outcomes <- c(10, NA, 30)
  weights <- c(1, 100, 3)
  result <- .compute_cell_weighted_mean(outcomes, weights)
  expect_equal(result$mean, 25.0, tolerance = 1e-12)
  expect_equal(result$n_obs, 2L)
})

test_that("compute_variance=TRUE returns weighted variance", {
  outcomes <- c(10, 20, 30)
  result <- .compute_cell_weighted_mean(outcomes,
                                         compute_variance = TRUE)
  expect_equal(result$variance, 200 / 3, tolerance = 1e-10)
})

test_that("n_obs==1 returns NULL variance", {
  outcomes <- c(42)
  result <- .compute_cell_weighted_mean(outcomes,
                                         compute_variance = TRUE)
  expect_null(result$variance)
  expect_equal(result$n_obs, 1L)
})

test_that("single obs with survey weight correct", {
  outcomes <- c(42)
  weights <- c(5.0)
  result <- .compute_cell_weighted_mean(outcomes, weights,
                                         compute_variance = TRUE)
  expect_equal(result$mean, 42.0)
  expect_null(result$variance)
  expect_equal(result$ess, 1.0, tolerance = 1e-10)
  expect_equal(result$n_obs, 1L)
})

test_that("survey weights return ESS", {
  outcomes <- c(10, 20, 30)
  weights <- c(1, 2, 1)
  result <- .compute_cell_weighted_mean(outcomes, weights)
  expect_equal(result$ess, 16 / 6, tolerance = 1e-10)
  expect_false(is.null(result$ess))
})

test_that("equal weights return NULL ess", {
  outcomes <- c(10, 20, 30)
  result <- .compute_cell_weighted_mean(outcomes)
  expect_null(result$ess)
})

test_that("outcomes/weights length mismatch throws error", {
  outcomes <- c(10, 20, 30)
  weights <- c(1, 2)
  expect_error(
    .compute_cell_weighted_mean(outcomes, weights),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================
# .compute_effective_sample_size() tests
# ============================================================

test_that("equal weights ESS equals n", {
  w <- rep(1.0, 10)
  expect_equal(.compute_effective_sample_size(w), 10.0,
               tolerance = 1e-10)
})

test_that("unequal weights ESS less than n", {
  w <- c(1, 1, 1, 100)
  ess <- .compute_effective_sample_size(w)
  expect_true(ess < 4)
  expect_true(ess > 1)
  expect_equal(ess, 10609 / 10003, tolerance = 1e-10)
})

test_that("all-zero weights ESS is 0", {
  w <- c(0, 0, 0)
  expect_equal(.compute_effective_sample_size(w), 0.0)
})

test_that("single weight ESS is 1", {
  expect_equal(.compute_effective_sample_size(c(5.0)), 1.0,
               tolerance = 1e-15)
  expect_equal(.compute_effective_sample_size(c(0.001)), 1.0,
               tolerance = 1e-15)
})

test_that("empty vector ESS is 0", {
  expect_equal(.compute_effective_sample_size(numeric(0)), 0.0)
})

test_that("equal weights ESS exactly n (large scale)", {
  w <- rep(3.7, 100)
  expect_equal(.compute_effective_sample_size(w), 100.0,
               tolerance = 1e-10)
})

# ============================================================
# .compute_control_weighted_mean() tests
# ============================================================

test_that("control equal weight mean correct", {
  values <- c(5, 10, 15)
  expect_equal(.compute_control_weighted_mean(values), 10.0)
})

test_that("control weighted mean correct", {
  values <- c(5, 10, 15)
  weights <- c(1, 2, 1)
  result <- .compute_control_weighted_mean(values, weights)
  expect_equal(result, 10.0, tolerance = 1e-12)
})

test_that("control all-NA returns NaN", {
  values <- c(NA, NA)
  expect_true(is.nan(.compute_control_weighted_mean(values)))
})

test_that("control partial NA with weights correctly aligned", {
  values <- c(5, NA, 15)
  weights <- c(1, 100, 3)
  result <- .compute_control_weighted_mean(values, weights)
  expect_equal(result, 12.5, tolerance = 1e-12)
})

test_that("control values/weights length mismatch throws error", {
  values <- c(5, 10, 15)
  weights <- c(1, 2)
  expect_error(
    .compute_control_weighted_mean(values, weights),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================
# Numerical precision end-to-end tests
# ============================================================

test_that("large-scale weighted mean precision", {
  set.seed(42)
  n <- 100000L
  outcomes <- rnorm(n, mean = 1000, sd = 100)
  weights <- runif(n, 0.1, 10)

  result <- .compute_cell_weighted_mean(outcomes, weights)

  norm_w <- weights / sum(weights)
  ref_mean <- sum(norm_w * outcomes)

  expect_equal(result$mean, ref_mean, tolerance = 1e-10)
})

test_that("large-scale Kahan vs naive sum precision", {
  set.seed(99)
  n <- 100000L
  raw_w <- runif(n, 0.01, 100)
  norm_w <- raw_w / sum(raw_w)
  outcomes <- rnorm(n, mean = 50000, sd = 10000)
  products <- norm_w * outcomes

  kahan_result <- .kahan_sum(products)
  naive_result <- sum(products)

  rel_diff <- abs(kahan_result - naive_result) / abs(naive_result)
  expect_true(rel_diff < 1e-10)
})
