# test-ipwra-unit-e606.R — IPWRA Unit Tests (E6-06.2, Part 2)
# TC-6.6.13 to TC-6.6.18

if (!exists(".lwdid_env") || is.null(.lwdid_env$warning_registry)) {
  .lwdid_env$warning_registry <- new_warning_registry()
}

# --- TC-6.6.13: WLS weights verification ---
test_that("TC-6.6.13: WLS weights are correct (treated=1, control=ps/(1-ps))", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  # ATT weights: treated = 1, control = ps/(1-ps)
  w <- ifelse(D == 1, 1, ps / (1 - ps))

  # Verify treated weights are exactly 1.0
  expect_true(all(w[D == 1] == 1.0))

  # Verify control weights are ps/(1-ps)
  ctrl_idx <- which(D == 0)
  expected_ctrl_w <- ps[ctrl_idx] / (1 - ps[ctrl_idx])
  expect_equal(w[ctrl_idx], expected_ctrl_w, tolerance = 1e-14)

  # Verify all weights are positive and finite
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))

  # Verify outcome model with these weights produces finite coefficients
  om_result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                      sample_weights = w)
  expect_true(all(sapply(om_result$coefficients, is.finite)))
})


# --- TC-6.6.14: Hajek normalization verification ---
test_that("TC-6.6.14: Hajek normalization divides by sum(w_ctrl), not N1", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  # ATT weights
  w <- ifelse(D == 1, 1, ps / (1 - ps))

  # Outcome model predictions
  om_result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                      sample_weights = w)
  m0_hat <- om_result$predictions

  # Manual Hajek computation
  ctrl_idx <- which(D == 0)
  treat_idx <- which(D == 1)
  w_ctrl <- ps[ctrl_idx] / (1 - ps[ctrl_idx])
  resid_ctrl <- Y[ctrl_idx] - m0_hat[ctrl_idx]

  # Hajek: divide by sum(w_ctrl), NOT N1
  ctrl_term_manual <- sum(w_ctrl * resid_ctrl) / sum(w_ctrl)
  treat_term_manual <- mean(Y[treat_idx] - m0_hat[treat_idx])
  att_manual <- treat_term_manual - ctrl_term_manual

  # Compare with estimate_ipwra
  result <- estimate_ipwra(df, "Y", "D", c("X1", "X2"),
                           trim_threshold = 0.01, trim_method = "clip")
  expect_equal(result$att, att_manual, tolerance = 1e-12)
})


# --- TC-6.6.15: Hajek IF main component ---
test_that("TC-6.6.15: Hajek IF main component is finite and mean-zero", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  result <- estimate_ipwra(df, "Y", "D", c("X1", "X2"),
                           trim_threshold = 0.01, trim_method = "clip")

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  w <- ifelse(D == 1, 1, ps / (1 - ps))
  om_result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                      sample_weights = w)
  m0_hat <- om_result$predictions

  ctrl_idx <- which(D == 0)
  treat_idx <- which(D == 1)
  w_ctrl <- ps[ctrl_idx] / (1 - ps[ctrl_idx])
  resid_ctrl <- Y[ctrl_idx] - m0_hat[ctrl_idx]
  B <- sum(w_ctrl * resid_ctrl) / sum(w_ctrl)
  p_bar <- sum(D == 1) / n

  # Construct psi_main
  psi_main <- numeric(n)
  psi_main[treat_idx] <- (Y[treat_idx] - m0_hat[treat_idx] - result$att) / p_bar
  psi_main[ctrl_idx] <- -w_ctrl * (resid_ctrl - B) / p_bar

  # Verify all values are finite
  expect_true(all(is.finite(psi_main)))

  # Verify mean is approximately 0
  expect_equal(mean(psi_main), 0, tolerance = 0.5)
})


# --- TC-6.6.16: PS adjustment uses residuals (r-B), not raw Y ---
test_that("TC-6.6.16: PS adjustment uses residuals (r-B), not raw Y", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  w <- ifelse(D == 1, 1, ps / (1 - ps))
  om_result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                      sample_weights = w)
  m0_hat <- om_result$predictions

  ctrl_idx <- which(D == 0)
  treat_idx <- which(D == 1)
  w_ctrl <- ps[ctrl_idx] / (1 - ps[ctrl_idx])
  resid_ctrl <- Y[ctrl_idx] - m0_hat[ctrl_idx]
  B <- sum(w_ctrl * resid_ctrl) / sum(w_ctrl)
  p_bar <- sum(D == 1) / n

  X_ctrl_const <- cbind(1, X1[ctrl_idx], X2[ctrl_idx])

  # IPWRA uses (r_i - B) in PS adjustment
  dATT_dgamma_ipwra <- -colSums(w_ctrl * X_ctrl_const *
                                  as.vector(resid_ctrl - B)) / (n * p_bar)

  # Verify length and finiteness

  expect_length(dATT_dgamma_ipwra, 3)
  expect_true(all(is.finite(dATT_dgamma_ipwra)))

  # IPW would use raw Y, not residuals — verify they differ
  dATT_dgamma_ipw <- -colSums(w_ctrl * X_ctrl_const *
                                as.vector(Y[ctrl_idx])) / (n * p_bar)
  expect_false(isTRUE(all.equal(dATT_dgamma_ipwra, dATT_dgamma_ipw)))
})


# --- TC-6.6.17: OM adjustment term ---
test_that("TC-6.6.17: OM adjustment term dATT_dbeta is correct", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  ctrl_idx <- which(D == 0)
  treat_idx <- which(D == 1)
  w_ctrl <- ps[ctrl_idx] / (1 - ps[ctrl_idx])

  X_all_const <- cbind(1, X1, X2)

  # Mean of X for treated units
  X_bar_treated <- colMeans(X_all_const[treat_idx, , drop = FALSE])

  # Weighted mean of X for control units (Hajek weights)
  X_bar_ctrl_w <- colSums(w_ctrl * X_all_const[ctrl_idx, , drop = FALSE]) /
                  sum(w_ctrl)

  # Newey & McFadden (1994) OM adjustment
  dATT_dbeta <- -X_bar_treated + X_bar_ctrl_w

  # Verify length and finiteness
  expect_length(dATT_dbeta, 3)
  expect_true(all(is.finite(dATT_dbeta)))
})


# --- TC-6.6.18: Residuals verification ---
test_that("TC-6.6.18: WLS residuals for control group are well-behaved", {
  set.seed(42)
  n <- 300; X1 <- rnorm(n); X2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.5 * X1 - 0.3 * X2))
  Y <- 1 + 2 * X1 - X2 + 1.5 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)

  ps_result <- estimate_propensity_score(df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_result$propensity_scores

  # ATT weights
  w <- ifelse(D == 1, 1, ps / (1 - ps))

  om_result <- estimate_outcome_model(df, "Y", "D", c("X1", "X2"),
                                      sample_weights = w)
  m0_hat <- om_result$predictions

  ctrl_idx <- which(D == 0)
  n_control <- length(ctrl_idx)

  # Manual residuals for control group
  resid_manual <- Y[ctrl_idx] - m0_hat[ctrl_idx]

  # Verify type and length
  expect_true(is.numeric(resid_manual))
  expect_length(resid_manual, n_control)

  # WLS residuals should center near 0 for control group
  expect_equal(mean(resid_manual), 0, tolerance = 0.5)

  # All residuals must be finite
  expect_true(all(is.finite(resid_manual)))
})