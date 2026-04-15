# test-numerical-consistency-e606.R — L2 Numerical Consistency Tests
# TC-6.6.24 to TC-6.6.32
#
# Verifies that IPW, IPWRA, and PSM estimators produce results that
# exactly match hand-computed values using the same algorithms.

# ============================================================================
# Warning registry initialization
# ============================================================================
if (exists(".lwdid_env")) {
  .lwdid_env$warning_registry <- new_warning_registry()
}

# ============================================================================
# Helper: shared data setup
# ============================================================================
setup_l2_data <- function() {
  set.seed(123)
  N <- 20
  X <- seq(-2, 2, length.out = N)
  D <- c(rep(0, 10), rep(1, 10))
  Y0 <- 1.0 + 0.5 * X + rnorm(N, 0, 0.5)
  Y1 <- Y0 + 1.0
  Y <- ifelse(D == 1, Y1, Y0)
  data.frame(Y = Y, D = D, X = X)
}

# ============================================================================
# TC-6.6.24: IPW ATT hand calculation
# ============================================================================
test_that("TC-6.6.24: IPW ATT matches hand-computed Hajek estimator", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D; X <- data_l2$X

  # Step 1: Estimate propensity scores
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  # Step 2: Manual Hajek ATT
  treat_mask <- D == 1
  ctrl_mask <- D == 0
  w_ctrl <- ps[ctrl_mask] / (1 - ps[ctrl_mask])
  att_hand <- mean(Y[treat_mask]) - sum(w_ctrl * Y[ctrl_mask]) / sum(w_ctrl)

  # Step 3: Compare with estimate_ipw
  ipw_result <- estimate_ipw(data_l2, "Y", "D", "X")
  expect_equal(ipw_result$att, att_hand, tolerance = 1e-10,
               label = "IPW ATT matches hand-computed Hajek estimator")
})

# ============================================================================
# TC-6.6.25: IPW SE hand calculation (full M-estimation IF)
# ============================================================================
test_that("TC-6.6.25: IPW SE matches hand-computed M-estimation IF", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D; X <- data_l2$X

  # PS estimation
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  treat_mask <- D == 1
  ctrl_mask <- D == 0
  n1 <- sum(treat_mask)
  n0 <- sum(ctrl_mask)
  p_bar <- n1 / N
  D_float <- as.numeric(D)

  # IPW weights and ATT
  w_ctrl <- ps[ctrl_mask] / (1 - ps[ctrl_mask])
  att <- mean(Y[treat_mask]) - sum(w_ctrl * Y[ctrl_mask]) / sum(w_ctrl)

  # Component 1: Horvitz-Thompson IF
  psi_ht <- numeric(N)
  psi_ht[treat_mask] <- (Y[treat_mask] - att) / p_bar
  psi_ht[ctrl_mask] <- -w_ctrl * Y[ctrl_mask] / p_bar

  # Component 2: PS estimation uncertainty adjustment
  X_const <- cbind(1, X)  # N x 2 design matrix with intercept

  # Logit score: S_{gamma,i} = (D_i - ps_i) * X_i
  S_gamma <- (D_float - ps) * X_const  # N x 2

  # Logit Hessian: H_gamma = -(1/N) X' diag(ps*(1-ps)) X
  W_ps <- ps * (1 - ps)
  H_gamma <- -crossprod(X_const * W_ps, X_const) / N

  # dATT/dgamma: uses raw Y (not residuals) for IPW
  X_ctrl_const <- X_const[ctrl_mask, , drop = FALSE]
  Y_ctrl <- Y[ctrl_mask]
  dATT_dgamma <- -colSums(w_ctrl * X_ctrl_const * Y_ctrl) / (N * p_bar)

  # Adjustment: S_gamma %*% t(solve(H_gamma)) %*% dATT_dgamma
  H_gamma_inv <- solve(H_gamma)
  adjustment <- as.numeric((S_gamma %*% t(H_gamma_inv)) %*% dATT_dgamma)

  # Complete IF
  psi <- psi_ht - adjustment

  # SE = sqrt(var(psi)/N) with ddof=1
  se_hand <- sqrt(var(psi) / N)

  # Compare
  ipw_result <- estimate_ipw(data_l2, "Y", "D", "X")
  expect_equal(ipw_result$se, se_hand, tolerance = 1e-10,
               label = "IPW SE matches hand-computed M-estimation IF")
})

# ============================================================================
# TC-6.6.26: IPWRA ATT hand calculation
# ============================================================================
test_that("TC-6.6.26: IPWRA ATT matches hand-computed doubly-robust estimator", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D; X <- data_l2$X

  # PS estimation
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  treat_mask <- D == 1
  ctrl_mask <- D == 0

  # ATT weights: treated=1, control=ps/(1-ps)
  att_weights <- ifelse(D == 1, 1.0, ps / (1 - ps))

  # WLS outcome model
  om_result <- estimate_outcome_model(data_l2, "Y", "D", "X",
                                       sample_weights = att_weights)
  m0_hat <- om_result$predictions

  # IPWRA ATT: Hajek normalization
  w_ctrl <- ps[ctrl_mask] / (1 - ps[ctrl_mask])
  w_ctrl_sum <- sum(w_ctrl)
  treat_term <- mean(Y[treat_mask] - m0_hat[treat_mask])
  resid_ctrl <- Y[ctrl_mask] - m0_hat[ctrl_mask]
  ctrl_term <- sum(w_ctrl * resid_ctrl) / w_ctrl_sum
  att_hand <- treat_term - ctrl_term

  # Compare
  ipwra_result <- estimate_ipwra(data_l2, "Y", "D", "X")
  expect_equal(ipwra_result$att, att_hand, tolerance = 1e-10,
               label = "IPWRA ATT matches hand-computed DR estimator")
})

# ============================================================================
# TC-6.6.27: IPWRA SE hand calculation (3-component IF)
# ============================================================================
test_that("TC-6.6.27: IPWRA SE matches hand-computed 3-component IF", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D; X <- data_l2$X

  # PS estimation
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  treat_mask <- D == 1
  ctrl_mask <- D == 0
  n1 <- sum(treat_mask)
  p_bar <- n1 / N
  D_float <- as.numeric(D)

  # ATT weights and outcome model
  att_weights <- ifelse(D == 1, 1.0, ps / (1 - ps))
  om_result <- estimate_outcome_model(data_l2, "Y", "D", "X",
                                       sample_weights = att_weights)
  m0_hat <- om_result$predictions

  # IPW weights and ATT
  w_ctrl <- ps[ctrl_mask] / (1 - ps[ctrl_mask])
  w_ctrl_sum <- sum(w_ctrl)
  treat_term <- mean(Y[treat_mask] - m0_hat[treat_mask])
  resid_ctrl <- Y[ctrl_mask] - m0_hat[ctrl_mask]
  ctrl_term <- sum(w_ctrl * resid_ctrl) / w_ctrl_sum
  att <- treat_term - ctrl_term
  B <- ctrl_term  # alias for clarity

  # ================================================================
  # Component 1: Hajek IF (main)
  # ================================================================
  psi_main <- numeric(N)
  psi_main[treat_mask] <- (Y[treat_mask] - m0_hat[treat_mask] - att) / p_bar
  psi_main[ctrl_mask] <- -w_ctrl * (resid_ctrl - B) / w_ctrl_sum * N

  # ================================================================
  # Component 2: PS estimation uncertainty adjustment
  # ================================================================
  X_const <- cbind(1, X)  # N x 2

  # Logit score
  S_gamma <- (D_float - ps) * X_const

  # Logit Hessian
  W_ps <- ps * (1 - ps)
  H_gamma <- -crossprod(X_const * W_ps, X_const) / N
  H_gamma_inv <- solve(H_gamma)

  # dATT/dgamma: uses residuals (r_i - B), NOT raw Y
  X_ctrl_const <- X_const[ctrl_mask, , drop = FALSE]
  r_minus_B <- resid_ctrl - B
  dATT_dgamma <- -colSums(w_ctrl * X_ctrl_const * r_minus_B) / w_ctrl_sum

  # PS adjustment
  ps_adj <- as.numeric((S_gamma %*% H_gamma_inv) %*% dATT_dgamma)

  # ================================================================
  # Component 3: Outcome model estimation uncertainty adjustment
  # ================================================================
  X_out_const <- cbind(1, X)  # N x 2 (same controls for OM)
  X_out_ctrl_const <- X_out_const[ctrl_mask, , drop = FALSE]
  X_out_treat_const <- X_out_const[treat_mask, , drop = FALSE]

  # WLS score: w_i * resid_i * X_i for controls, 0 for treated
  S_beta <- matrix(0, nrow = N, ncol = ncol(X_out_const))
  S_beta[ctrl_mask, ] <- w_ctrl * resid_ctrl * X_out_ctrl_const

  # WLS Hessian: -(1/N) X_ctrl' diag(w_ctrl) X_ctrl
  H_beta <- -crossprod(X_out_ctrl_const * w_ctrl, X_out_ctrl_const) / N
  H_beta_inv <- solve(H_beta)

  # dATT/dbeta = -X_bar_treated + X_bar_ctrl_weighted
  X_bar_treated <- colMeans(X_out_treat_const)
  X_bar_ctrl_w <- colSums(w_ctrl * X_out_ctrl_const) / w_ctrl_sum
  dATT_dbeta <- -X_bar_treated + X_bar_ctrl_w

  # OM adjustment
  om_adj <- as.numeric((S_beta %*% H_beta_inv) %*% dATT_dbeta)

  # ================================================================
  # Combine full IF and compute SE
  # ================================================================
  psi_full <- psi_main - ps_adj - om_adj
  se_hand <- sqrt(var(as.numeric(psi_full)) / N)

  # Compare
  ipwra_result <- estimate_ipwra(data_l2, "Y", "D", "X")
  expect_equal(ipwra_result$se, se_hand, tolerance = 1e-8,
               label = "IPWRA SE matches hand-computed 3-component IF")
})

# ============================================================================
# TC-6.6.28: PSM ATT hand calculation
# ============================================================================
test_that("TC-6.6.28: PSM ATT matches hand-computed nearest-neighbor matching", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D

  # PS estimation
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  treat_mask <- D == 1
  ctrl_mask <- D == 0
  ps_treat <- ps[treat_mask]
  ps_ctrl <- ps[ctrl_mask]
  Y_treat <- Y[treat_mask]
  Y_ctrl <- Y[ctrl_mask]

  # 1:1 nearest neighbor matching (no caliper)
  n_treat <- sum(treat_mask)
  matched_ctrl_Y <- numeric(n_treat)
  for (i in seq_len(n_treat)) {
    dists <- abs(ps_treat[i] - ps_ctrl)
    best_idx <- which.min(dists)
    matched_ctrl_Y[i] <- Y_ctrl[best_idx]
  }

  individual_effects <- Y_treat - matched_ctrl_Y
  att_hand <- mean(individual_effects)

  # Compare
  psm_result <- estimate_psm(data_l2, "Y", "D", "X")
  expect_equal(psm_result$att, att_hand, tolerance = 1e-10,
               label = "PSM ATT matches hand-computed NN matching")
})

# ============================================================================
# TC-6.6.29: PSM SE hand calculation
# ============================================================================
test_that("TC-6.6.29: PSM SE matches hand-computed matched-pair SE", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()
  N <- nrow(data_l2)
  Y <- data_l2$Y; D <- data_l2$D

  # PS estimation
  ps_result <- estimate_propensity_score(data_l2, "D", "X", 0.01, "clip")
  ps <- ps_result$propensity_scores

  treat_mask <- D == 1
  ctrl_mask <- D == 0
  ps_treat <- ps[treat_mask]
  ps_ctrl <- ps[ctrl_mask]
  Y_treat <- Y[treat_mask]
  Y_ctrl <- Y[ctrl_mask]

  # 1:1 nearest neighbor matching
  n_treat <- sum(treat_mask)
  matched_ctrl_Y <- numeric(n_treat)
  for (i in seq_len(n_treat)) {
    dists <- abs(ps_treat[i] - ps_ctrl)
    best_idx <- which.min(dists)
    matched_ctrl_Y[i] <- Y_ctrl[best_idx]
  }

  individual_effects <- Y_treat - matched_ctrl_Y
  n_valid <- length(individual_effects)
  se_hand <- sqrt(var(individual_effects) / n_valid)

  # Compare
  psm_result <- estimate_psm(data_l2, "Y", "D", "X")
  expect_equal(psm_result$se, se_hand, tolerance = 1e-10,
               label = "PSM SE matches hand-computed matched-pair SE")
})

# ============================================================================
# TC-6.6.30: Normal inference verification (qnorm, not qt)
# ============================================================================
test_that("TC-6.6.30: All estimators use normal inference (qnorm CI, pnorm pvalue)", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()

  # Run all three estimators
  ipw_result <- estimate_ipw(data_l2, "Y", "D", "X")
  ipwra_result <- estimate_ipwra(data_l2, "Y", "D", "X")
  psm_result <- estimate_psm(data_l2, "Y", "D", "X")

  z_crit <- qnorm(0.975)

  # --- IPW ---
  ipw_ci_lower <- ipw_result$att - z_crit * ipw_result$se
  ipw_ci_upper <- ipw_result$att + z_crit * ipw_result$se
  ipw_pvalue <- 2 * (1 - pnorm(abs(ipw_result$att / ipw_result$se)))

  expect_equal(ipw_result$ci_lower, ipw_ci_lower, tolerance = 1e-10,
               label = "IPW CI lower uses qnorm")
  expect_equal(ipw_result$ci_upper, ipw_ci_upper, tolerance = 1e-10,
               label = "IPW CI upper uses qnorm")
  expect_equal(ipw_result$pvalue, ipw_pvalue, tolerance = 1e-10,
               label = "IPW pvalue uses pnorm")

  # --- IPWRA ---
  ipwra_ci_lower <- ipwra_result$att - z_crit * ipwra_result$se
  ipwra_ci_upper <- ipwra_result$att + z_crit * ipwra_result$se
  ipwra_pvalue <- 2 * (1 - pnorm(abs(ipwra_result$att / ipwra_result$se)))

  expect_equal(ipwra_result$ci_lower, ipwra_ci_lower, tolerance = 1e-10,
               label = "IPWRA CI lower uses qnorm")
  expect_equal(ipwra_result$ci_upper, ipwra_ci_upper, tolerance = 1e-10,
               label = "IPWRA CI upper uses qnorm")
  expect_equal(ipwra_result$pvalue, ipwra_pvalue, tolerance = 1e-10,
               label = "IPWRA pvalue uses pnorm")

  # --- PSM ---
  psm_ci_lower <- psm_result$att - z_crit * psm_result$se
  psm_ci_upper <- psm_result$att + z_crit * psm_result$se
  psm_pvalue <- 2 * (1 - pnorm(abs(psm_result$att / psm_result$se)))

  expect_equal(psm_result$ci_lower, psm_ci_lower, tolerance = 1e-10,
               label = "PSM CI lower uses qnorm")
  expect_equal(psm_result$ci_upper, psm_ci_upper, tolerance = 1e-10,
               label = "PSM CI upper uses qnorm")
  expect_equal(psm_result$pvalue, psm_pvalue, tolerance = 1e-10,
               label = "PSM pvalue uses pnorm")
})

# ============================================================================
# TC-6.6.31: IPW df_resid formula verification
# ============================================================================
test_that("TC-6.6.31: IPW df_resid equals n1 + n0 - 2", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()

  ipw_result <- estimate_ipw(data_l2, "Y", "D", "X")

  n1 <- ipw_result$n_treated
  n0 <- ipw_result$n_control
  expected_df <- as.integer(n1 + n0 - 2L)

  expect_identical(ipw_result$df_resid, expected_df,
                   label = "IPW df_resid == n1 + n0 - 2")
})

# ============================================================================
# TC-6.6.32: IPWRA df_resid formula verification
# ============================================================================
test_that("TC-6.6.32: IPWRA df_resid equals N - (2 + ncol(X_out))", {
  if (exists(".lwdid_env")) .lwdid_env$warning_registry <- new_warning_registry()
  data_l2 <- setup_l2_data()

  ipwra_result <- estimate_ipwra(data_l2, "Y", "D", "X")

  N <- nrow(data_l2)
  # X_out = as.matrix(data[, controls, drop = FALSE]) — controls without intercept
  # For controls = "X", ncol(X_out) = 1
  ncol_X_out <- 1L  # single covariate "X"
  expected_df <- as.integer(N - (2L + ncol_X_out))

  expect_identical(ipwra_result$df_resid, expected_df,
                   label = "IPWRA df_resid == N - (2 + ncol(X_out))")
})
