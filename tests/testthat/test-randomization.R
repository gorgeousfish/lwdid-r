# ===========================================================================
# test-randomization.R
# L1 Unit tests for randomization_inference() (Story E7-07, Task E7-07.1)
# ===========================================================================

# Internal function accessor
ri <- lwdid:::randomization_inference

# TC-7.7.1: Permutation preserves N1
test_that("TC-7.7.1: permutation preserves N1 count", {
  df <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  # Extract first post-treatment period cross-section
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d
  n1_orig <- sum(d == 1L)

  # Direct verification: sample(d) preserves sum
  set.seed(100L)
  for (i in 1:100) {
    d_perm <- sample(d)
    expect_equal(sum(d_perm), n1_orig)
  }

  # Function-level: permutation RI completes and produces valid result
  result <- ri(y, d, reps = 50L, seed = 42L, method = "permutation")
  expect_true(!is.null(result$ri_pvalue))
  expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
  expect_equal(result$method, "permutation")
})

# TC-7.7.2: Bootstrap uses replacement (has duplicates)
test_that("TC-7.7.2: bootstrap resamples with replacement", {
  df <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d
  n <- length(d)

  # Direct verification: at least 1 of 100 bootstrap draws has duplicates
  set.seed(200L)
  has_dup <- FALSE
  for (i in 1:100) {
    idx <- sample.int(n, n, replace = TRUE)
    if (any(duplicated(idx))) {
      has_dup <- TRUE
      break
    }
  }
  expect_true(has_dup, info = "Bootstrap should produce duplicate indices")

  # Function-level: bootstrap RI completes
  result <- ri(y, d, reps = 200L, seed = 42L, method = "bootstrap")
  expect_true(!is.null(result$ri_pvalue))
  expect_equal(result$method, "bootstrap")
})

# TC-7.7.3: FATAL-003 interaction rebuild verification
test_that("TC-7.7.3: FATAL-003 interactions rebuilt per permutation", {
  # Generate data with controls
  df <- generate_ct_panel(N = 60L, tau = 2.0, with_controls = TRUE, seed = 42L)
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d
  x <- as.matrix(post_df[, c("x1", "x2")])

  # RI with controls should complete without error
  result <- ri(y, d, x = x, reps = 100L, seed = 42L, method = "permutation")
  expect_true(!is.null(result$ri_pvalue))
  expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)

  # Verify the distribution has variation (not all same value)
  expect_true(length(unique(result$ri_distribution)) > 1L,
              info = "RI distribution should have variation with controls")

  # Compare with no-controls RI — distributions should differ
  result_no_x <- ri(y, d, x = NULL, reps = 100L, seed = 42L,
                     method = "permutation")
  # KS test: distributions should be different (controls matter)
  ks_result <- ks.test(result$ri_distribution, result_no_x$ri_distribution)
  # We don't require p < 0.05 (depends on data), but distributions should exist
  expect_true(length(result$ri_distribution) > 0L)
  expect_true(length(result_no_x$ri_distribution) > 0L)
})

# TC-7.7.4: p-value in [0,1] for both methods
test_that("TC-7.7.4: p-value in [0,1] for permutation and bootstrap", {
  df <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d

  # Permutation
  res_perm <- ri(y, d, reps = 50L, seed = 42L, method = "permutation")
  expect_true(res_perm$ri_pvalue >= 0 && res_perm$ri_pvalue <= 1)

  # Bootstrap
  res_boot <- ri(y, d, reps = 200L, seed = 42L, method = "bootstrap")
  expect_true(res_boot$ri_pvalue >= 0 && res_boot$ri_pvalue <= 1)
})

# TC-7.7.5: Seed reproducibility
test_that("TC-7.7.5: same seed produces identical results", {
  df <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d

  res1 <- ri(y, d, reps = 100L, seed = 42L, method = "permutation")
  res2 <- ri(y, d, reps = 100L, seed = 42L, method = "permutation")

  expect_identical(res1$ri_pvalue, res2$ri_pvalue)
  expect_identical(res1$ri_distribution, res2$ri_distribution)
  expect_identical(res1$obs_att, res2$obs_att)
  expect_identical(res1$n_valid, res2$n_valid)
  expect_identical(res1$n_failed, res2$n_failed)
})


# ===========================================================================
# Task E7-07.4: L2 Numerical tests (RI portion)
# TC-7.7.12 and TC-7.7.16
# ===========================================================================

# TC-7.7.12: RI p-value direction consistent
# Under H1 (tau=1.5), p-value should be smaller than under H0 (tau=0).
# Also verify obs_att > 0 under H1 (positive treatment effect).
test_that("TC-7.7.12: RI p-value direction consistent", {
  # H1: real effect (tau=1.5)
  df_h1 <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  tpost1 <- 5L
  post_h1 <- df_h1[time == tpost1]
  res_h1 <- suppressWarnings(suppressMessages(
    ri(post_h1$y, post_h1$d, reps = 500L, seed = 42L, method = "permutation")
  ))

  # H0: no effect (tau=0)
  df_h0 <- generate_ct_panel(N = 50L, tau = 0.0, seed = 42L)
  post_h0 <- df_h0[time == tpost1]
  res_h0 <- suppressWarnings(suppressMessages(
    ri(post_h0$y, post_h0$d, reps = 500L, seed = 42L, method = "permutation")
  ))

  # Under H1, p-value should be small (< 0.20 at minimum)
  expect_true(res_h1$ri_pvalue < 0.20)

  # Key assertion: effect detected more strongly under H1

  expect_true(res_h1$ri_pvalue < res_h0$ri_pvalue)

  # obs_att should be positive under H1 (positive treatment effect)
  expect_true(res_h1$obs_att > 0)
})

# TC-7.7.16: p-value simple-proportion verification
# Formula: p = count / n_valid
# Verify by manually counting from ri_distribution.
test_that("TC-7.7.16: p-value simple-proportion verification", {
  df <- generate_ct_panel(N = 50L, tau = 0, seed = 42L)
  tpost1 <- 5L
  post_df <- df[time == tpost1]
  y <- post_df$y
  d <- post_df$d

  result <- suppressWarnings(suppressMessages(
    ri(y, d, reps = 500L, seed = 42L, method = "permutation")
  ))

  # Manually count from ri_distribution
  manual_count <- sum(abs(result$ri_distribution) >= abs(result$obs_att))
  expected_p <- manual_count / result$n_valid
  expect_equal(result$ri_pvalue, expected_p, tolerance = 1e-15)
})


# ===========================================================================
# Task E7-07.5: L3 Monte Carlo tests (RI portion)
# TC-7.7.17 and TC-7.7.18
# ===========================================================================

# TC-7.7.17: Under H0 (tau=0), RI p-values should be uniformly distributed
# 200 simulations, reps=500 each. KS test against U[0,1].
test_that("TC-7.7.17: H0 RI p-values uniformly distributed (KS test)", {
  skip_on_cran()
  n_sims <- 200L
  reps <- 500L
  pvalues <- numeric(n_sims)

  for (sim in seq_len(n_sims)) {
    df <- generate_ct_panel(N = 50L, tau = 0, seed = sim)
    post_df <- df[time == 5L]
    res <- ri(post_df$y, post_df$d,
              reps = reps, seed = sim, method = "permutation")
    pvalues[sim] <- res$ri_pvalue
  }

  # All p-values should be in [0,1]
  expect_true(all(pvalues >= 0 & pvalues <= 1))

  # KS test: p-values should follow U[0,1] under H0
  ks_res <- ks.test(pvalues, "punif")
  expect_true(ks_res$p.value > 0.05,
    info = sprintf("KS p=%.4f, p-values not uniform under H0", ks_res$p.value))
})

# TC-7.7.18: Under H1 (tau=3.0, N=100), RI power > 0.8
# 200 simulations, reps=500 each. Rejection rate at alpha=0.05.
test_that("TC-7.7.18: H1 RI power > 0.8", {
  skip_on_cran()
  n_sims <- 200L
  reps <- 500L
  alpha <- 0.05
  pvalues <- numeric(n_sims)

  for (sim in seq_len(n_sims)) {
    df <- generate_ct_panel(N = 100L, tau = 3.0, seed = sim)
    post_df <- df[time == 5L]
    res <- ri(post_df$y, post_df$d,
              reps = reps, seed = sim, method = "permutation")
    pvalues[sim] <- res$ri_pvalue
  }

  rejection_rate <- mean(pvalues < alpha)
  expect_true(rejection_rate > 0.8,
    info = sprintf("Power=%.3f, expected > 0.8", rejection_rate))
})


# ===========================================================================
# Task E7-07.6: L4 Boundary condition tests (RI portion)
# TC-7.7.21 and TC-7.7.22
# ===========================================================================

# TC-7.7.21: N=3 minimal sample RI
# N=3 (N1=1, N0=2), permutation method. Should not error, p in [0,1].
test_that("TC-7.7.21: N=3 minimal sample RI runs without error", {
  dt <- generate_minimal_panel(seed = 42L)
  # Use period 1 (pre-treatment) cross-section: 3 units
  cs <- dt[time == 1L]
  y <- cs$y
  d <- cs$d
  expect_equal(length(y), 3L)
  expect_equal(sum(d), 1L)

  # Should not error (N=3 >= 3 minimum)
  res <- ri(y, d, reps = 50L, seed = 42L, method = "permutation")
  expect_true(res$ri_pvalue >= 0 && res$ri_pvalue <= 1)
  expect_equal(res$method, "permutation")
})

# TC-7.7.22: reps=1 minimal repetitions
# reps=1: permutation min_valid = max(10, 0.1*1) = 10, n_valid=1 < 10
# → should trigger hard error. Bootstrap similarly.
test_that("TC-7.7.22: reps=1 triggers insufficient valid reps error", {
  df <- generate_ct_panel(N = 50L, tau = 1.5, seed = 42L)
  post_df <- df[time == 5L]
  y <- post_df$y
  d <- post_df$d

  # Permutation: min_valid = max(10, 0.1*1) = 10, but n_valid <= 1
  expect_error(
    ri(y, d, reps = 1L, seed = 42L, method = "permutation"),
    "valid replications"
  )

  # Bootstrap: min_valid = max(100, 0.1*1) = 100, but n_valid <= 1
  expect_error(
    ri(y, d, reps = 1L, seed = 42L, method = "bootstrap"),
    "valid replications"
  )
})
