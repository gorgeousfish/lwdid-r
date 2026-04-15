# ============================================================================
# estimate_ipwra.R — IPWRA Doubly-Robust ATT Estimator
# ============================================================================
# Implements Hajek IPWRA-ATT with 3-component IF standard errors.
#
# Paper references:
#   - lw2025 Procedure 3.1 Step 2 / Procedure 4.1 Step 3
#   - Wooldridge (2007): IPWRA doubly-robust estimator
#   - Wooldridge (2025b §19.4): WLS for IPWRA ATT
#   - Cattaneo (2010): semiparametric influence function
#   - Newey & McFadden (1994): M-estimation framework
#   - Hajek (1971): normalized weights
# ============================================================================

#' Estimate IPWRA-ATT (Doubly-Robust)
#'
#' @param data data.frame, cross-sectional data (one row per unit)
#' @param y character, outcome variable name
#' @param d character, treatment indicator name (1=treated, 0=control)
#' @param controls character vector, outcome model control variables
#' @param propensity_controls character vector or NULL, PS model controls
#' @param trim_threshold numeric, PS trimming threshold (default 0.01)
#' @param trim_method character, "clip" (default) or "drop"
#' @param vce NULL (analytical 3-component IF SE) or "bootstrap"
#' @param alpha numeric, significance level (default 0.05)
#' @param boot_reps integer, bootstrap replications (default 200)
#' @param seed integer or NULL, random seed for bootstrap
#'
#' @return list with 15 fields (see design doc)
#' @keywords internal
estimate_ipwra <- function(data, y, d, controls,
                           propensity_controls = NULL,
                           trim_threshold = 0.01,
                           trim_method = "clip",
                           vce = NULL, alpha = 0.05,
                           boot_reps = 200L, seed = NULL) {

  # ---- Step 1: Parameter defaults ----
  if (is.null(propensity_controls)) {
    propensity_controls <- controls
  }

  # ---- Step 2: Input validation ----
  if (length(controls) == 0L) {
    stop_lwdid("IPWRA requires at least one control variable.",
               class = "lwdid_invalid_param", param = "controls")
  }
  vce_method <- if (is.null(vce)) "analytical" else match.arg(vce, "bootstrap")
  trim_method <- match.arg(trim_method, c("clip", "drop"))

  # ---- Step 3: Missing value handling ----
  all_vars <- unique(c(y, d, controls, propensity_controls))
  missing_cols <- setdiff(all_vars, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(sprintf("Column(s) not found in data: %s",
                       paste(missing_cols, collapse = ", ")),
               class = "lwdid_invalid_param", param = "data")
  }
  complete_mask <- stats::complete.cases(data[, all_vars, drop = FALSE])
  data_clean <- data[complete_mask, , drop = FALSE]
  n_excluded <- sum(!complete_mask)

  # ---- Step 4: Sample size checks ----
  D_all <- as.integer(data_clean[[d]])
  n_treated_raw <- sum(D_all == 1L)
  n_control_raw <- sum(D_all == 0L)

  if (n_treated_raw < 2L) {
    stop_lwdid(sprintf("Insufficient treated units (%d). Need >= 2.", n_treated_raw),
               class = "lwdid_insufficient_data")
  }
  if (n_control_raw < 2L) {
    stop_lwdid(sprintf("Insufficient control units (%d). Need >= 2.", n_control_raw),
               class = "lwdid_insufficient_data")
  }
  if (n_treated_raw < 5L) {
    warn_lwdid(sprintf("Small treated sample (%d). IPWRA estimates may be unstable.",
                       n_treated_raw),
               class = "lwdid_small_sample")
  }
  if (n_control_raw < 10L) {
    warn_lwdid(sprintf("Small control sample (%d). Outcome model may be unstable.",
                       n_control_raw),
               class = "lwdid_small_sample")
  }

  # ---- Step 5: Estimate propensity scores ----
  ps_result <- estimate_propensity_score(data_clean, d, propensity_controls,
                                          trim_threshold, trim_method)
  ps_raw <- ps_result$propensity_scores

  # ---- Step 6: Trimming strategy dispatch ----
  n_total <- nrow(data_clean)
  if (trim_method == "clip") {
    valid_mask <- rep(TRUE, n_total)
    n_trimmed <- 0L
  } else {
    # drop mode
    valid_mask <- !is.na(ps_raw) & ps_raw >= trim_threshold &
                  ps_raw <= (1 - trim_threshold)
    n_trimmed <- sum(!valid_mask)
    pct_trimmed <- 100 * n_trimmed / n_total
    if (pct_trimmed > 10) {
      warn_lwdid(
        sprintf("IPWRA drop trimming excluded %.1f%% of observations (%d/%d).",
                pct_trimmed, n_trimmed, n_total),
        class = "lwdid_data", detail = "ipwra_high_trimming")
    }
  }

  # Apply mask
  data_valid <- data_clean[valid_mask, , drop = FALSE]
  pscores <- ps_raw[valid_mask]
  Y <- as.numeric(data_valid[[y]])
  D <- as.integer(data_valid[[d]])
  n <- length(Y)

  treat_mask <- D == 1L
  control_mask <- D == 0L
  n1 <- sum(treat_mask)
  n0 <- sum(control_mask)

  if (n1 == 0L) {
    stop_lwdid("No treated units after trimming.", class = "lwdid_no_treated")
  }
  if (n0 == 0L) {
    stop_lwdid("No control units after trimming.", class = "lwdid_no_control")
  }

  # ---- Step 7: ATT weights ----
  att_weights <- ifelse(D == 1L, 1.0, pscores / (1 - pscores))

  # ---- Step 8: WLS outcome model ----
  om_result <- estimate_outcome_model(data_valid, y, d, controls,
                                       sample_weights = att_weights)
  m0_hat <- om_result$predictions

  # ---- Step 9: IPWRA-ATT (Hajek normalization) ----
  w_ctrl <- pscores[control_mask] / (1 - pscores[control_mask])
  w_ctrl_sum <- sum(w_ctrl)
  if (w_ctrl_sum <= 0) {
    stop_lwdid("IPW weight sum is non-positive. PS model may be misspecified.",
               class = "lwdid_estimation_failed")
  }

  treat_term <- mean(Y[treat_mask] - m0_hat[treat_mask])
  resid_ctrl <- Y[control_mask] - m0_hat[control_mask]
  ctrl_term <- sum(w_ctrl * resid_ctrl) / w_ctrl_sum
  att <- treat_term - ctrl_term

  # ---- Step 10: Weight CV diagnostic ----
  if (length(w_ctrl) > 1L && mean(w_ctrl) > 0) {
    weights_cv <- sd(w_ctrl) / mean(w_ctrl)
    if (weights_cv > 2.0) {
      warn_lwdid(
        sprintf("IPWRA weight CV = %.2f > 2.0. Possible overlap violation.",
                weights_cv),
        class = "lwdid_overlap", detail = "ipwra_extreme_weights")
    }
  }

  # ---- Step 11: Extreme PS diagnostic ----
  pct_low <- mean(pscores < 0.05)
  pct_high <- mean(pscores > 0.95)
  if (pct_low > 0.1 || pct_high > 0.1) {
    warn_lwdid(
      sprintf("%.1f%% of PS < 0.05 and %.1f%% > 0.95. Overlap assumption may be violated.",
              100 * pct_low, 100 * pct_high),
      class = "lwdid_overlap", detail = "ipwra_extreme_ps")
  }

  # ---- Step 12: SE computation dispatch ----
  if (vce_method == "analytical") {
    se_result <- .compute_ipwra_se_analytical(
      data = data_valid, y = y, d = d, controls = controls,
      att = att, pscores = pscores, m0_hat = m0_hat,
      weights = w_ctrl, alpha = alpha,
      propensity_controls = propensity_controls,
      ps_result = ps_result
    )
    se <- se_result$se
    ci_lower <- se_result$ci_lower
    ci_upper <- se_result$ci_upper
  } else {
    boot_result <- .compute_ipwra_se_bootstrap(
      data = data_valid, y = y, d = d, controls = controls,
      propensity_controls = propensity_controls,
      trim_threshold = trim_threshold,
      trim_method = trim_method,
      n_bootstrap = boot_reps, seed = seed, alpha = alpha
    )
    se <- boot_result$se
    ci_lower <- boot_result$ci_lower
    ci_upper <- boot_result$ci_upper
  }

  # ---- Step 13: Normal distribution inference ----
  if (se > 0) {
    t_stat <- att / se
    pvalue <- 2 * (1 - stats::pnorm(abs(t_stat)))
  } else {
    t_stat <- NaN
    pvalue <- NaN
  }

  # ---- Step 14: df_resid (interface compatibility) ----
  X_out <- as.matrix(data_valid[, controls, drop = FALSE])
  df_resid <- as.integer(n - (2L + ncol(X_out)))

  # ---- Step 15: Return value assembly (15 fields) ----
  # Full-sample weight vector
  weights_all <- numeric(n)
  weights_all[which(treat_mask)] <- 1.0
  weights_all[which(control_mask)] <- w_ctrl

  list(
    att                   = att,
    se                    = se,
    ci_lower              = ci_lower,
    ci_upper              = ci_upper,
    t_stat                = t_stat,
    pvalue                = pvalue,
    propensity_scores     = pscores,
    weights               = weights_all,
    outcome_model_coef    = om_result$coefficients,
    propensity_model_coef = ps_result$coefficients,
    n_treated             = as.integer(n1),
    n_control             = as.integer(n0),
    n_trimmed             = as.integer(n_trimmed),
    df_resid              = df_resid,
    estimator             = "ipwra"
  )
}


# ============================================================================
# .compute_ipwra_se_analytical() — 3-Component IF Standard Error
# ============================================================================
# Implements Cattaneo (2010) semiparametric IF + Newey & McFadden (1994)
# M-estimation first-order Taylor expansion.
#
# Three components:
#   1. Hajek IF (IPWRA-ATT ratio estimator linearization)
#   2. PS estimation uncertainty adjustment
#   3. Outcome model estimation uncertainty adjustment
# ============================================================================

#' @keywords internal
.compute_ipwra_se_analytical <- function(data, y, d, controls,
                                         att, pscores, m0_hat, weights,
                                         alpha = 0.05,
                                         propensity_controls = NULL,
                                         ps_result = NULL) {

  # ---- Step 0: Basic quantities ----
  D <- as.numeric(data[[d]])
  Y <- as.numeric(data[[y]])
  n <- length(D)

  treat_mask <- D == 1
  control_mask <- D == 0
  p_bar <- sum(treat_mask) / n

  # w_ctrl = pscores[ctrl] / (1 - pscores[ctrl])
  w_ctrl <- weights  # passed in as w_ctrl from caller
  w_ctrl_sum <- sum(w_ctrl)

  # ================================================================
  # Component 1: Hajek Influence Function
  # ================================================================
  resid_ctrl <- Y[control_mask] - m0_hat[control_mask]
  ctrl_term <- sum(w_ctrl * resid_ctrl) / w_ctrl_sum  # B

  psi <- numeric(n)
  psi[treat_mask] <- (Y[treat_mask] - m0_hat[treat_mask] - att) / p_bar
  psi[control_mask] <- -w_ctrl * (resid_ctrl - ctrl_term) / w_ctrl_sum * n

  # ================================================================
  # Component 2: PS Estimation Uncertainty Adjustment
  # ================================================================
  if (is.null(propensity_controls)) propensity_controls <- controls

  # Handle removed constant columns from PS estimation
  kept_ps_controls <- propensity_controls
  if (!is.null(ps_result) && length(ps_result$removed_cols) > 0L) {
    kept_ps_controls <- setdiff(propensity_controls, ps_result$removed_cols)
  }

  X_ps <- as.matrix(data[, kept_ps_controls, drop = FALSE])
  X_ps_const <- cbind(1, X_ps)  # N x (k_ps+1)

  # Logit score: S_{gamma,i} = (D_i - p_hat_i) * X_i
  S_gamma <- (D - pscores) * X_ps_const  # N x (k_ps+1)

  # Logit Hessian: H_gamma = -(1/N) X' diag(p*(1-p)) X
  W_ps <- pscores * (1 - pscores)
  H_gamma <- -crossprod(X_ps_const * W_ps, X_ps_const) / n

  # Hessian inverse with MASS::ginv() fallback
  H_gamma_inv <- tryCatch(
    solve(H_gamma),
    error = function(e) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(H_gamma)
      } else {
        stop_lwdid("PS Hessian singular and MASS package unavailable.",
                   class = "lwdid_estimation_failed")
      }
    }
  )

  # dATT/dgamma: uses residuals (r_i - B), NOT raw Y_i
  r_minus_B <- resid_ctrl - ctrl_term
  X_ps_ctrl_const <- X_ps_const[control_mask, , drop = FALSE]
  dATT_dgamma <- -colSums(w_ctrl * X_ps_ctrl_const * r_minus_B) / w_ctrl_sum

  # PS adjustment
  ps_adjustment <- as.numeric(
    (S_gamma %*% H_gamma_inv) %*% dATT_dgamma
  )

  # ================================================================
  # Component 3: Outcome Model Estimation Uncertainty Adjustment
  # ================================================================
  X_out <- as.matrix(data[, controls, drop = FALSE])
  X_out_const <- cbind(1, X_out)  # N x (k_out+1)
  X_ctrl_const <- X_out_const[control_mask, , drop = FALSE]

  # WLS score (only control group non-zero)
  S_beta <- matrix(0, nrow = n, ncol = ncol(X_out_const))
  S_beta[control_mask, ] <- w_ctrl * resid_ctrl * X_ctrl_const

  # WLS Hessian
  H_beta <- -crossprod(X_ctrl_const * w_ctrl, X_ctrl_const) / n

  # Hessian inverse with MASS::ginv() fallback
  H_beta_inv <- tryCatch(
    solve(H_beta),
    error = function(e) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(H_beta)
      } else {
        stop_lwdid("OM Hessian singular and MASS package unavailable.",
                   class = "lwdid_estimation_failed")
      }
    }
  )

  # dATT/dbeta = -X_bar_treated + X_bar_ctrl_weighted
  X_bar_treated <- colMeans(X_out_const[treat_mask, , drop = FALSE])
  X_bar_ctrl_w <- colSums(w_ctrl * X_ctrl_const) / w_ctrl_sum
  dATT_dbeta <- -X_bar_treated + X_bar_ctrl_w

  # OM adjustment
  om_adjustment <- as.numeric(
    (S_beta %*% H_beta_inv) %*% dATT_dbeta
  )

  # ================================================================
  # Combine full IF
  # ================================================================
  psi_full <- psi - ps_adjustment - om_adjustment

  # Variance: var(psi_full)/n, var() uses ddof=1
  se <- sqrt(stats::var(as.numeric(psi_full)) / n)

  # Confidence interval (normal distribution)
  z_crit <- stats::qnorm(1 - alpha / 2)
  ci_lower <- att - z_crit * se
  ci_upper <- att + z_crit * se

  list(se = se, ci_lower = ci_lower, ci_upper = ci_upper)
}


# ============================================================================
# .compute_ipwra_se_bootstrap() — Bootstrap Standard Error
# ============================================================================

#' @keywords internal
.compute_ipwra_se_bootstrap <- function(data, y, d, controls,
                                        propensity_controls,
                                        trim_threshold,
                                        trim_method = "clip",
                                        n_bootstrap, seed, alpha) {

  n <- nrow(data)
  if (!is.null(seed)) set.seed(seed)

  att_boot <- numeric(n_bootstrap)
  n_failed <- 0L

  for (b in seq_len(n_bootstrap)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_data <- data[idx, , drop = FALSE]

    # Check minimum group sizes
    D_boot <- as.integer(boot_data[[d]])
    if (sum(D_boot == 1L) < 2L || sum(D_boot == 0L) < 2L) {
      att_boot[b] <- NA_real_
      n_failed <- n_failed + 1L
      next
    }

    att_boot[b] <- .ipwra_point_estimate(
      data = boot_data, y = y, d = d, controls = controls,
      propensity_controls = propensity_controls,
      trim_threshold = trim_threshold,
      trim_method = trim_method
    )
    if (is.na(att_boot[b])) n_failed <- n_failed + 1L
  }

  # Success rate check
  n_valid <- sum(!is.na(att_boot))
  success_rate <- n_valid / n_bootstrap

  if (success_rate < 0.5) {
    warn_lwdid(
      sprintf("IPWRA bootstrap success rate = %.1f%% (< 50%%). SE may be unreliable.",
              100 * success_rate),
      class = "lwdid_convergence", detail = "ipwra_bootstrap_low_success")
  }
  if (n_valid < 10L) {
    stop_lwdid(
      sprintf("IPWRA bootstrap: only %d valid replicates (need >= 10).", n_valid),
      class = "lwdid_estimation_failed")
  }

  att_valid <- att_boot[!is.na(att_boot)]
  se <- sd(att_valid)
  ci <- as.numeric(stats::quantile(att_valid,
                                    probs = c(alpha / 2, 1 - alpha / 2)))

  list(se = se, ci_lower = ci[1], ci_upper = ci[2])
}


# ============================================================================
# .ipwra_point_estimate() — Bootstrap Helper (Point Estimate Only)
# ============================================================================

#' @keywords internal
.ipwra_point_estimate <- function(data, y, d, controls,
                                  propensity_controls,
                                  trim_threshold,
                                  trim_method = "clip") {

  Y <- as.numeric(data[[y]])
  D <- as.integer(data[[d]])

  # Step 1: PS estimation
  fml <- stats::as.formula(
    paste(d, "~", paste(propensity_controls, collapse = " + "))
  )
  glm_fit <- tryCatch(
    suppressWarnings(
      stats::glm(fml, data = data, family = stats::binomial(link = "logit"))
    ),
    error = function(e) return(NULL)
  )
  if (is.null(glm_fit)) return(NA_real_)

  ps <- as.numeric(stats::predict(glm_fit, type = "response"))

  # Step 2: Trimming
  if (trim_method == "clip") {
    ps <- pmax(trim_threshold, pmin(1 - trim_threshold, ps))
    valid <- rep(TRUE, length(ps))
  } else {
    valid <- ps >= trim_threshold & ps <= (1 - trim_threshold)
  }
  Y <- Y[valid]; D <- D[valid]; ps <- ps[valid]
  data_sub <- data[valid, , drop = FALSE]

  treat <- D == 1L; ctrl <- D == 0L
  if (sum(treat) == 0L || sum(ctrl) == 0L) return(NA_real_)

  # Step 3: ATT weights
  att_weights <- ifelse(D == 1L, 1.0, ps / (1 - ps))

  # Step 4: WLS outcome model
  om_result <- tryCatch(
    estimate_outcome_model(data_sub, y, d, controls,
                           sample_weights = att_weights),
    error = function(e) return(NULL)
  )
  if (is.null(om_result)) return(NA_real_)

  m0_hat <- om_result$predictions

  # Step 5: Hajek ATT
  w_ctrl <- ps[ctrl] / (1 - ps[ctrl])
  w_ctrl_sum <- sum(w_ctrl)
  if (w_ctrl_sum <= 0) return(NA_real_)

  treat_term <- mean(Y[treat] - m0_hat[treat])
  resid_ctrl <- Y[ctrl] - m0_hat[ctrl]
  ctrl_term <- sum(w_ctrl * resid_ctrl) / w_ctrl_sum

  treat_term - ctrl_term
}
