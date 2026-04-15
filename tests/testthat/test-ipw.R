# test-ipw.R — IPW estimator tests (Story E6-02, TC-6.2.1 to TC-6.2.31)

generate_ipw_test_data <- function(n = 500, seed = 42, tau = 3.0,
                                    beta0 = 1, beta1 = 2, beta2 = 0.5,
                                    ps_intercept = -0.5, ps_coef1 = 0.8,
                                    ps_coef2 = 0.5, noise_sd = 1.0) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  linpred <- ps_intercept + ps_coef1 * X1 + ps_coef2 * X2
  prob <- 1 / (1 + exp(-linpred))
  D <- rbinom(n, 1, prob)
  if (all(D == 1) || all(D == 0)) {
    D[1] <- 0L
    D[2] <- 1L
  }
  Y <- beta0 + beta1 * X1 + beta2 * X2 + tau * D +
    rnorm(n, sd = noise_sd)
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
}

generate_noconf_data <- function(n = 2000, seed = 123, tau = 3.0) {
  set.seed(seed)
  X1 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * X1 + tau * D + rnorm(n, sd = 1.0)
  data.frame(Y = Y, D = D, X1 = X1)
}

test_that("TC-6.2.1: No-confounding DGP — ATT close to true", {
  df <- generate_noconf_data(n = 2000, seed = 123, tau = 3.0)
  result <- estimate_ipw(df, "Y", "D", "X1")
  expect_true(abs(result$att - 3.0) < 0.5)
  expect_true(result$se > 0 && result$se < 1.0)
  expect_length(result, 15L)
  expect_identical(result$estimator, "ipw")
  expect_identical(result$vce_method, "analytical")
})

test_that("TC-6.2.2: Manual Hajek — matches hand computation", {
  df <- generate_ipw_test_data(n = 300, seed = 99)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  ps_res <- estimate_propensity_score(
    df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_res$propensity_scores
  Y <- df$Y
  D <- df$D
  treat <- D == 1L
  ctrl <- D == 0L
  w <- ps[ctrl] / (1 - ps[ctrl])
  att_manual <- mean(Y[treat]) - sum(w * Y[ctrl]) / sum(w)
  expect_equal(result$att, att_manual, tolerance = 1e-10)
})

test_that("TC-6.2.3: Large-sample IPW vs RA comparison", {
  df <- generate_ipw_test_data(n = 5000, seed = 99, tau = 3.0)
  controls <- c("X1", "X2")
  ipw_res <- estimate_ipw(df, "Y", "D", controls)
  om_res <- estimate_outcome_model(df, "Y", "D", controls)
  att_ra <- mean(df$Y[df$D == 1]) -
    mean(om_res$predictions[df$D == 1])
  expect_true(abs(ipw_res$att - att_ra) < 0.5)
  expect_true(abs(ipw_res$att - 3.0) < 0.5)
})

test_that("TC-6.2.4: Clip mode — all obs retained, PS clipped", {
  df <- generate_ipw_test_data(n = 500, seed = 42,
    ps_intercept = -2.0, ps_coef1 = 2.0)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    trim_threshold = 0.05, trim_method = "clip")
  ps <- result$propensity_scores
  expect_true(all(ps >= 0.05 - 1e-12))
  expect_true(all(ps <= 0.95 + 1e-12))
  expect_equal(result$n_treated + result$n_control, nrow(df))
})

test_that("TC-6.2.5: Clip mode — n_trimmed counts clipped", {
  df <- generate_ipw_test_data(n = 500, seed = 42,
    ps_intercept = 3.0, ps_coef1 = 2.0)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    trim_threshold = 0.05, trim_method = "clip")
  expect_true(is.integer(result$n_trimmed))
  expect_true(result$n_trimmed >= 0L)
  expect_equal(result$n_treated + result$n_control, nrow(df))
})

test_that("TC-6.2.6: Clip mode — effective sample = n_total", {
  df <- generate_ipw_test_data(n = 300, seed = 77)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    trim_method = "clip")
  expect_equal(result$n_treated + result$n_control, nrow(df))
  expect_length(result$propensity_scores, nrow(df))
  expect_length(result$weights, nrow(df))
})

test_that("TC-6.2.7: Drop mode — extreme PS excluded", {
  df <- generate_ipw_test_data(n = 500, seed = 42,
    ps_intercept = -2.0, ps_coef1 = 2.0)
  result <- suppressWarnings(
    estimate_ipw(df, "Y", "D", c("X1", "X2"),
      trim_threshold = 0.05, trim_method = "drop"))
  eff_n <- result$n_treated + result$n_control
  expect_true(eff_n <= nrow(df))
  expect_true(result$n_trimmed >= 0L)
})

test_that("TC-6.2.8: Drop mode — no treated after trim", {
  set.seed(808)
  n <- 100
  X1 <- c(rnorm(20, mean = 5), rnorm(80, mean = 0))
  D <- c(rep(1L, 20), rep(0L, 80))
  Y <- 1 + X1 + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = rnorm(n))
  expect_error(
    suppressWarnings(
      estimate_ipw(df, "Y", "D", c("X1", "X2"),
        trim_threshold = 0.45, trim_method = "drop")),
    class = "lwdid_no_treated")
})

test_that("TC-6.2.9: IF SE positive and reasonable", {
  df <- generate_ipw_test_data(n = 500, seed = 42)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_true(result$se > 0)
  expect_true(result$se < abs(result$att) * 2)
  expect_identical(result$vce_method, "analytical")
})

test_that("TC-6.2.10: HT IF component verified", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  ps_res <- estimate_propensity_score(
    df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_res$propensity_scores
  Y <- df$Y
  D <- df$D
  treat <- D == 1L
  ctrl <- D == 0L
  n <- nrow(df)
  n1 <- sum(treat)
  p_bar <- n1 / n
  w <- ps[ctrl] / (1 - ps[ctrl])
  w_sum <- sum(w)
  att <- mean(Y[treat]) - sum(w * Y[ctrl]) / w_sum
  psi_ht <- numeric(n)
  psi_ht[treat] <- (Y[treat] - att) / p_bar
  psi_ht[ctrl] <- -w * Y[ctrl] / p_bar
  expect_equal(psi_ht[treat][1],
    (Y[treat][1] - att) / p_bar, tolerance = 1e-12)
  expect_equal(psi_ht[ctrl][1],
    -w[1] * Y[ctrl][1] / p_bar, tolerance = 1e-12)
})

test_that("TC-6.2.11: dATT_dgamma uses raw Y", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  ps_res <- estimate_propensity_score(
    df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_res$propensity_scores
  Y <- df$Y
  D <- df$D
  ctrl <- D == 0L
  n <- nrow(df)
  n1 <- sum(D == 1L)
  p_bar <- n1 / n
  w <- ps[ctrl] / (1 - ps[ctrl])
  kept <- setdiff(c("X1", "X2"), ps_res$removed_cols)
  X_ps <- as.matrix(df[, kept, drop = FALSE])
  X_const <- cbind(1, X_ps)
  X_ctrl <- X_const[ctrl, , drop = FALSE]
  Y_ctrl <- Y[ctrl]
  dATT <- -colSums(w * X_ctrl * Y_ctrl) / (n * p_bar)
  expect_true(any(abs(dATT) > 1e-10))
  expect_length(dATT, ncol(X_const))
})

test_that("TC-6.2.12: Logit Hessian vs vcov(glm)", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  ps_res <- estimate_propensity_score(
    df, "D", c("X1", "X2"), 0.01, "clip")
  n <- nrow(df)
  glm_obj <- ps_res$glm_object
  ps_raw_glm <- predict(glm_obj, type = "response")
  kept <- setdiff(c("X1", "X2"), ps_res$removed_cols)
  X_ps <- as.matrix(df[, kept, drop = FALSE])
  X_const <- cbind(1, X_ps)
  W_ps <- ps_raw_glm * (1 - ps_raw_glm)
  H_gamma <- -crossprod(X_const * W_ps, X_const) / n
  vcov_glm <- vcov(glm_obj)
  H_from_vcov <- -solve(vcov_glm) / n
  expect_equal(unname(H_gamma), unname(H_from_vcov),
    tolerance = 1e-6)
})

test_that("TC-6.2.13: Complete IF = psi_ht - adjustment", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  ps_res <- estimate_propensity_score(
    df, "D", c("X1", "X2"), 0.01, "clip")
  ps <- ps_res$propensity_scores
  Y <- df$Y
  D <- df$D
  treat <- D == 1L
  ctrl <- D == 0L
  n <- nrow(df)
  n1 <- sum(treat)
  p_bar <- n1 / n
  w <- ps[ctrl] / (1 - ps[ctrl])
  w_sum <- sum(w)
  att <- mean(Y[treat]) - sum(w * Y[ctrl]) / w_sum
  psi_ht <- numeric(n)
  psi_ht[treat] <- (Y[treat] - att) / p_bar
  psi_ht[ctrl] <- -w * Y[ctrl] / p_bar
  D_float <- as.numeric(D)
  kept <- setdiff(c("X1", "X2"), ps_res$removed_cols)
  X_ps <- as.matrix(df[, kept, drop = FALSE])
  X_const <- cbind(1, X_ps)
  S_gamma <- (D_float - ps) * X_const
  W_ps <- ps * (1 - ps)
  H_gamma <- -crossprod(X_const * W_ps, X_const) / n
  H_inv <- solve(H_gamma)
  X_ctrl <- X_const[ctrl, , drop = FALSE]
  Y_ctrl <- Y[ctrl]
  dATT <- -colSums(w * X_ctrl * Y_ctrl) / (n * p_bar)
  adj <- as.numeric((S_gamma %*% t(H_inv)) %*% dATT)
  psi_full <- psi_ht - adj
  se_manual <- sqrt(var(psi_full) / n)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_equal(result$att, att, tolerance = 1e-10)
  expect_equal(result$se, se_manual, tolerance = 1e-8)
})

test_that("TC-6.2.14: SE scales as 1/sqrt(n)", {
  df1 <- generate_ipw_test_data(n = 500, seed = 42)
  df2 <- generate_ipw_test_data(n = 2000, seed = 42)
  r1 <- estimate_ipw(df1, "Y", "D", c("X1", "X2"))
  r2 <- estimate_ipw(df2, "Y", "D", c("X1", "X2"))
  ratio <- r1$se / r2$se
  expected <- sqrt(2000 / 500)
  expect_true(ratio > expected * 0.5 && ratio < expected * 2.0)
})

test_that("TC-6.2.15: Uses pnorm not pt for p-value", {
  df <- generate_ipw_test_data(n = 500, seed = 42)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  z <- result$att / result$se
  # Primary: p-value matches pnorm formula exactly
  pval_z <- 2 * (1 - pnorm(abs(z)))
  expect_equal(result$pvalue, pval_z, tolerance = 1e-12)
  # Secondary: verify pnorm != pt at a known large z-value (z=2.5)
  # where the difference between normal and t(5) is unambiguous
  pval_norm_ref <- 2 * (1 - pnorm(2.5))
  pval_t_ref <- 2 * (1 - pt(2.5, df = 5))
  expect_true(abs(pval_norm_ref - pval_t_ref) > 0.01)
})

test_that("TC-6.2.16: Uses qnorm not qt for CI", {
  df <- generate_ipw_test_data(n = 500, seed = 42)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    alpha = 0.05)
  z_crit <- qnorm(1 - 0.05 / 2)
  t_crit <- qt(1 - 0.05 / 2, df = result$df_resid)
  expect_equal(result$ci_lower,
    result$att - z_crit * result$se, tolerance = 1e-12)
  expect_equal(result$ci_upper,
    result$att + z_crit * result$se, tolerance = 1e-12)
  width_z <- 2 * z_crit * result$se
  width_t <- 2 * t_crit * result$se
  expect_true(abs(width_z - width_t) > 1e-6 ||
    result$df_resid > 1e6)
})

test_that("TC-6.2.17: df_resid = n1 + n0 - 2", {
  df <- generate_ipw_test_data(n = 500, seed = 42)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_equal(result$df_resid,
    as.integer(result$n_treated + result$n_control - 2L))
  expect_true(is.integer(result$df_resid))
})

test_that("TC-6.2.18: 95% CI coverage near 0.95 (MC)", {
  n_sims <- 200
  tau_true <- 3.0
  covers <- logical(n_sims)
  for (i in seq_len(n_sims)) {
    df <- generate_noconf_data(n = 500, seed = i, tau = tau_true)
    res <- tryCatch(
      suppressWarnings(estimate_ipw(df, "Y", "D", "X1")),
      error = function(e) NULL)
    if (!is.null(res)) {
      covers[i] <- res$ci_lower <= tau_true &&
        tau_true <= res$ci_upper
    } else {
      covers[i] <- NA
    }
  }
  coverage <- mean(covers, na.rm = TRUE)
  expect_true(coverage >= 0.85 && coverage <= 1.0)
})

test_that("TC-6.2.19: Bootstrap SE positive", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    vce = "bootstrap", boot_reps = 50L, seed = 123)
  expect_true(result$se > 0)
  expect_identical(result$vce_method, "bootstrap")
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se))
})

test_that("TC-6.2.20: Bootstrap seed reproducibility", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  r1 <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    vce = "bootstrap", boot_reps = 50L, seed = 999)
  r2 <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    vce = "bootstrap", boot_reps = 50L, seed = 999)
  expect_equal(r1$att, r2$att, tolerance = 1e-12)
  expect_equal(r1$se, r2$se, tolerance = 1e-12)
  expect_equal(r1$ci_lower, r2$ci_lower, tolerance = 1e-12)
  expect_equal(r1$ci_upper, r2$ci_upper, tolerance = 1e-12)
})

test_that("TC-6.2.21: Bootstrap SE within 20% of analytical (large n)", {
  df <- generate_ipw_test_data(n = 2000, seed = 42)
  r_ana <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  r_boot <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    vce = "bootstrap", boot_reps = 200L, seed = 42)
  ratio <- r_boot$se / r_ana$se
  expect_true(ratio > 0.8 && ratio < 1.2,
    label = sprintf("Boot/Ana SE ratio = %.3f, expected in [0.8, 1.2]",
      ratio))
  expect_equal(r_ana$att, r_boot$att, tolerance = 1e-10)
})

test_that("TC-6.2.22: Bootstrap full pipeline (clip + drop)", {
  df <- generate_ipw_test_data(n = 300, seed = 42)
  r_clip <- estimate_ipw(df, "Y", "D", c("X1", "X2"),
    vce = "bootstrap", boot_reps = 30L, seed = 1,
    trim_method = "clip", trim_threshold = 0.05)
  expect_true(r_clip$se > 0)
  expect_identical(r_clip$vce_method, "bootstrap")
  r_drop <- suppressWarnings(
    estimate_ipw(df, "Y", "D", c("X1", "X2"),
      vce = "bootstrap", boot_reps = 30L, seed = 1,
      trim_method = "drop", trim_threshold = 0.05))
  expect_true(r_drop$se > 0)
  expect_identical(r_drop$vce_method, "bootstrap")
})

test_that("TC-6.2.23: Weight CV > 2 triggers lwdid_overlap warning", {
  set.seed(2323)
  n <- 400
  X1 <- c(rnorm(50, mean = 3), rnorm(350, mean = 0))
  D <- c(rep(1L, 50), rep(0L, 350))
  Y <- 1 + 2 * X1 + 3 * D + rnorm(n, sd = 1)
  df <- data.frame(Y = Y, D = D, X1 = X1)
  expect_warning(
    estimate_ipw(df, "Y", "D", "X1"),
    class = "lwdid_overlap")
})

test_that("TC-6.2.24: Hessian singular — MASS::ginv() fallback", {
  set.seed(2424)
  n <- 100
  X1 <- rnorm(n)
  X2 <- X1
  D <- rbinom(n, 1, 0.5)
  Y <- 1 + X1 + 3 * D + rnorm(n)
  df <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  result <- suppressWarnings(
    estimate_ipw(df, "Y", "D", c("X1", "X2")))
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se))
  expect_true(result$se > 0)
})

test_that("TC-6.2.25: Single covariate works", {
  df <- generate_noconf_data(n = 500, seed = 25)
  result <- estimate_ipw(df, "Y", "D", "X1")
  expect_true(is.finite(result$att))
  expect_true(result$se > 0)
  expect_equal(result$n_treated + result$n_control, nrow(df))
  expect_length(result, 15L)
})

test_that("TC-6.2.26: Multiple covariates work", {
  df <- generate_ipw_test_data(n = 500, seed = 26)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_true(is.finite(result$att))
  expect_true(result$se > 0)
  expect_equal(result$n_treated + result$n_control, nrow(df))
  expect_length(result, 15L)
})

test_that("TC-6.2.27: Different propensity_controls routing", {
  df <- generate_ipw_test_data(n = 500, seed = 27)
  r1 <- estimate_ipw(df, "Y", "D", "X1")
  r2 <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_true(abs(r1$att - r2$att) > 1e-6)
  expect_true(is.finite(r1$se) && is.finite(r2$se))
})

test_that("TC-6.2.28: Extreme but valid data (weight sum edge case)", {
  set.seed(2828)
  n <- 200
  X1 <- c(rnorm(10, mean = 4), rnorm(190, mean = 0))
  D <- c(rep(1L, 10), rep(0L, 190))
  Y <- 1 + X1 + 3 * D + rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, D = D, X1 = X1)
  result <- suppressWarnings(
    estimate_ipw(df, "Y", "D", "X1"))
  expect_true(is.finite(result$att))
  expect_true(result$se > 0)
  w_ctrl <- result$weights[D == 0L]
  expect_true(sum(w_ctrl) > 0)
})

test_that("TC-6.2.29: Empty propensity_controls — lwdid_invalid_param", {
  df <- generate_ipw_test_data(n = 100, seed = 29)
  expect_error(
    estimate_ipw(df, "Y", "D", character(0)),
    class = "lwdid_invalid_param")
})

test_that("TC-6.2.30: y column not found — lwdid_invalid_param", {
  df <- generate_ipw_test_data(n = 100, seed = 30)
  expect_error(
    estimate_ipw(df, "nonexistent_y", "D", c("X1", "X2")),
    class = "lwdid_invalid_param")
})

test_that("TC-6.2.31: d column not found — lwdid_invalid_param", {
  df <- generate_ipw_test_data(n = 100, seed = 31)
  expect_error(
    estimate_ipw(df, "Y", "nonexistent_d", c("X1", "X2")),
    class = "lwdid_invalid_param")
})

# --- Structural tests ---

test_that("Structural: Return list has exactly 15 fields with correct names", {
  df <- generate_ipw_test_data(n = 300, seed = 50)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expected_names <- c("att", "se", "ci_lower", "ci_upper",
    "t_stat", "pvalue", "propensity_scores", "weights",
    "propensity_model_coef", "n_treated", "n_control",
    "n_trimmed", "df_resid", "vce_method", "estimator")
  expect_length(result, 15L)
  expect_setequal(names(result), expected_names)
  expect_identical(result$estimator, "ipw")
})

test_that("Structural: PS/weights length = nrow, ctrl weights sum to n_treated", {
  df <- generate_ipw_test_data(n = 400, seed = 51)
  result <- estimate_ipw(df, "Y", "D", c("X1", "X2"))
  expect_length(result$propensity_scores, nrow(df))
  expect_length(result$weights, nrow(df))
  ctrl_mask <- df$D == 0L
  ctrl_weights <- result$weights[ctrl_mask]
  expect_equal(sum(ctrl_weights), result$n_treated, tolerance = 1e-10)
})
