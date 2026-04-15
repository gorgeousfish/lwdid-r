# test-psm-unit-e606.R — PSM Unit Tests (E6-06.2, Part 3)
# TC-6.6.19 to TC-6.6.23

if (!exists(".lwdid_env") || is.null(.lwdid_env$warning_registry)) {
  .lwdid_env$warning_registry <- new_warning_registry()
}

# --------------------------------------------------------------------------
# TC-6.6.19: 1:1 and 1:3 matching (match_counts verification)
# --------------------------------------------------------------------------
test_that("TC-6.6.19: 1:1 match_counts are all 1; 1:3 are all min(3, n_ctrl)", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1 + X1 + D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1)

  # 1:1 matching
  res1 <- estimate_psm(df, "Y", "D", "X1",
                        n_neighbors = 1L, with_replacement = TRUE)
  n_treat <- sum(D == 1)
  # All matched treated units should have match_counts == 1
  matched_mask <- res1$match_counts > 0L
  expect_true(all(res1$match_counts[matched_mask] == 1L))
  # No units should be dropped (no caliper)
  expect_equal(res1$n_dropped, 0L)

  # 1:3 matching
  res3 <- estimate_psm(df, "Y", "D", "X1",
                        n_neighbors = 3L, with_replacement = TRUE)
  n_ctrl <- sum(D == 0)
  expected_k <- min(3L, n_ctrl)
  matched_mask3 <- res3$match_counts > 0L
  expect_true(all(res3$match_counts[matched_mask3] == expected_k))
  expect_equal(res3$n_dropped, 0L)
})

# --------------------------------------------------------------------------
# TC-6.6.20: Narrow caliper drops treated units (n_dropped > 0)
# --------------------------------------------------------------------------
test_that("TC-6.6.20: Narrow caliper drops treated units", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n, sd = 2)
  D <- rbinom(n, 1, plogis(1.5 * X1))
  Y <- 1 + X1 + D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1)

  .lwdid_env$warning_registry <- new_warning_registry()

  # Very narrow caliper (absolute scale) should drop some treated units
  res <- suppressWarnings(
    estimate_psm(df, "Y", "D", "X1",
                 caliper = 0.001, caliper_scale = "absolute",
                 with_replacement = TRUE)
  )

  expect_gt(res$n_dropped, 0L)

  # Dropped units should have match_counts == 0
  dropped_idx <- which(res$match_counts == 0L)
  expect_equal(length(dropped_idx), res$n_dropped)

  # Dropped units should have empty matched_control_ids
  for (idx in dropped_idx) {
    expect_length(res$matched_control_ids[[idx]], 0L)
  }

  # ATT should only be based on matched units
  n_valid <- sum(res$match_counts > 0L)
  expect_equal(res$match_success_rate, n_valid / res$n_treated)
})

# --------------------------------------------------------------------------
# TC-6.6.21: with_replacement=TRUE allows same control matched multiple times
# --------------------------------------------------------------------------
test_that("TC-6.6.21: with_replacement=TRUE allows duplicate control matches", {
  set.seed(42)
  n <- 100
  X1 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1 + X1 + D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1)

  res <- estimate_psm(df, "Y", "D", "X1",
                       n_neighbors = 1L, with_replacement = TRUE)

  # Collect all matched control indices
  all_ctrl_ids <- unlist(res$matched_control_ids)
  all_ctrl_ids <- all_ctrl_ids[!is.na(all_ctrl_ids)]

  # With replacement: some controls should appear more than once
  # (especially with moderate n and 1:1 matching)
  freq_table <- table(all_ctrl_ids)
  has_duplicates <- any(freq_table > 1)

  # With n=100 and replacement, duplicates are very likely
  expect_true(has_duplicates)
})

# --------------------------------------------------------------------------
# TC-6.6.22: with_replacement=FALSE ensures no duplicate control indices
# --------------------------------------------------------------------------
test_that("TC-6.6.22: with_replacement=FALSE has no duplicate controls", {
  set.seed(42)
  n <- 200
  X1 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1))
  Y <- 1 + X1 + D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1)

  res <- estimate_psm(df, "Y", "D", "X1",
                       n_neighbors = 1L, with_replacement = FALSE)

  # Collect all matched control indices
  all_ctrl_ids <- unlist(res$matched_control_ids)
  all_ctrl_ids <- all_ctrl_ids[!is.na(all_ctrl_ids)]

  # Without replacement: no control index should appear more than once
  expect_equal(length(all_ctrl_ids), length(unique(all_ctrl_ids)))
})


# --------------------------------------------------------------------------
# TC-6.6.23: SE = sqrt(var(individual_effects)/N_valid) manual verification
# --------------------------------------------------------------------------
test_that("TC-6.6.23: SE matches manual sqrt(var(effects)/N_valid)", {
  set.seed(42)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + X1 + 2 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  res <- estimate_psm(df, "Y", "D", c("X1", "X2"),
                       n_neighbors = 1L, with_replacement = TRUE)

  # Manually reconstruct individual effects
  treat_idx <- which(D == 1)
  ctrl_idx <- which(D == 0)
  Y_treat <- Y[treat_idx]
  Y_ctrl <- Y[ctrl_idx]

  individual_effects <- numeric(0)
  for (i in seq_along(res$matched_control_ids)) {
    matches <- res$matched_control_ids[[i]]
    if (length(matches) > 0L) {
      effect_i <- Y_treat[i] - mean(Y_ctrl[matches])
      individual_effects <- c(individual_effects, effect_i)
    }
  }

  n_valid <- length(individual_effects)
  expect_gt(n_valid, 1L)

  # Manual SE: sqrt(var(effects) / N_valid), var() uses ddof=1
  se_manual <- sqrt(var(individual_effects) / n_valid)

  # Manual ATT
  att_manual <- mean(individual_effects)

  # Compare with estimate_psm output

  expect_equal(res$att, att_manual, tolerance = 1e-12)
  expect_equal(res$se, se_manual, tolerance = 1e-12)

  # Verify SE is positive and finite
  expect_true(is.finite(res$se))
  expect_gt(res$se, 0)

  # Verify CI uses normal distribution (qnorm, not qt)
  z_crit <- qnorm(0.975)
  ci_lower_manual <- att_manual - z_crit * se_manual
  ci_upper_manual <- att_manual + z_crit * se_manual
  expect_equal(res$ci_lower, ci_lower_manual, tolerance = 1e-12)
  expect_equal(res$ci_upper, ci_upper_manual, tolerance = 1e-12)
})
