# ============================================================================
# estimate_ipw.R — IPW ATT Estimator
# ============================================================================
# Implements Hajek IPW-ATT with M-estimation IF standard errors.
#
# Paper references:
#   - lw2025 Procedure 3.1 Step 2 / Procedure 4.1 Step 3
#   - Lunceford & Davidian (2004): Hajek IPW-ATT + M-estimation IF SE
#   - Wooldridge (2010) Ch.21: IPW asymptotic theory
#   - Newey & McFadden (1994): M-estimation framework
# ============================================================================

#' Estimate IPW-ATT (Inverse Probability Weighting)
#'
#' @param data data.frame, cross-sectional data (one row per unit)
#' @param y character, outcome variable name
#' @param d character, treatment indicator name (1=treated, 0=control)
#' @param propensity_controls character vector, PS model covariates
#' @param trim_threshold numeric, PS trimming threshold (default 0.01)
#' @param trim_method character, "clip" (default) or "drop"
#' @param vce NULL (analytical IF SE) or "bootstrap"
#' @param alpha numeric, significance level (default 0.05)
#' @param boot_reps integer, bootstrap replications (default 200)
#' @param seed integer or NULL, random seed for bootstrap
#'
#' @return list with 15 fields (see design doc)
#' @keywords internal
estimate_ipw <- function(data, y, d, propensity_controls,
                         trim_threshold = 0.01,
                         trim_method = c("clip", "drop"),
                         vce = NULL, alpha = 0.05,
                         boot_reps = 200L, seed = NULL) {

  # ---- Step 1: Input validation ----
  if (!y %in% names(data)) {
    stop_lwdid(sprintf("Column '%s' not found in data.", y),
               class = "lwdid_invalid_param", param = "y")
  }
  if (!d %in% names(data)) {
    stop_lwdid(sprintf("Column '%s' not found in data.", d),
               class = "lwdid_invalid_param", param = "d")
  }
  if (length(propensity_controls) == 0L) {
    stop_lwdid("IPW requires at least one propensity control variable.",
               class = "lwdid_invalid_param", param = "propensity_controls")
  }

  # ---- Step 2: Data extraction ----
  Y_all <- as.numeric(data[[y]])
  D_all <- as.integer(data[[d]])
  n_total <- length(Y_all)

  # ---- Step 3: Trim method dispatch ----
  trim_method <- match.arg(trim_method, c("clip", "drop"))

  # ---- Step 4: Estimate propensity scores ----
  ps_result <- estimate_propensity_score(data, d, controls = propensity_controls,
                                          trim_threshold = trim_threshold,
                                          trim_method = trim_method)
  ps_raw <- ps_result$propensity_scores

  # ---- Step 5: Build valid mask ----
  if (trim_method == "clip") {
    valid_mask <- rep(TRUE, n_total)
    n_trimmed <- ps_result$n_trimmed
  } else {
    # drop mode
    valid_mask <- !is.na(ps_raw) & ps_raw >= trim_threshold & ps_raw <= (1 - trim_threshold)
    n_trimmed <- sum(!valid_mask)
    pct_trimmed <- 100 * n_trimmed / n_total
    if (pct_trimmed > 10) {
      warn_lwdid(
        sprintf("IPW drop trimming excluded %.1f%% of observations (%d/%d).",
                pct_trimmed, n_trimmed, n_total),
        class = "lwdid_data", detail = "ipw_high_trimming")
    }
  }

  # ---- Step 6: Apply mask and check ----
  Y <- Y_all[valid_mask]
  D <- D_all[valid_mask]
  ps <- ps_raw[valid_mask]
  n <- length(Y)

  treated_mask <- D == 1L
  control_mask <- D == 0L
  n1 <- sum(treated_mask)
  n0 <- sum(control_mask)

  if (n1 == 0L) {
    stop_lwdid("No treated units after trimming.",
               class = "lwdid_no_treated")
  }
  if (n0 == 0L) {
    stop_lwdid("No control units after trimming.",
               class = "lwdid_no_control")
  }

  # ---- Step 7: IPW weights ----
  w_ctrl <- ps[control_mask] / (1 - ps[control_mask])
  w_sum <- sum(w_ctrl)
  if (w_sum <= 0) {
    stop_lwdid("IPW weight sum is non-positive. PS model may be misspecified.",
               class = "lwdid_estimation_failed")
  }

  # Normalized weights (control group sums to n1)
  weights_valid <- numeric(n)
  weights_valid[which(treated_mask)] <- 1.0
  weights_valid[which(control_mask)] <- w_ctrl * n1 / w_sum

  # ---- Step 8: Hajek ATT ----
  att <- mean(Y[treated_mask]) - sum(w_ctrl * Y[control_mask]) / w_sum

  # ---- Step 9: SE computation ----
  vce_method <- if (is.null(vce)) "analytical" else match.arg(vce, "bootstrap")

  if (vce_method == "analytical") {
    # ---- Analytical IF SE (Lunceford & Davidian 2004 Theorem 1) ----
    D_float <- as.numeric(D)
    p_bar <- n1 / n

    # Component 1: Horvitz-Thompson IF
    psi_ht <- numeric(n)
    psi_ht[treated_mask] <- (Y[treated_mask] - att) / p_bar
    psi_ht[control_mask] <- -w_ctrl * Y[control_mask] / p_bar

    # Component 2: PS estimation uncertainty adjustment
    X_ps <- as.matrix(data[valid_mask, propensity_controls, drop = FALSE])
    # Handle removed constant columns from PS estimation
    kept_controls <- setdiff(propensity_controls, ps_result$removed_cols)
    if (length(kept_controls) < length(propensity_controls)) {
      X_ps <- as.matrix(data[valid_mask, kept_controls, drop = FALSE])
    }
    X_const <- cbind(1, X_ps)

    # Logit score: S_{gamma,i} = (D_i - p_hat_i) * X_i
    S_gamma <- (D_float - ps) * X_const

    # Logit Hessian: H_gamma = -(1/N) X' diag(p*(1-p)) X
    W_ps <- ps * (1 - ps)
    H_gamma <- -crossprod(X_const * W_ps, X_const) / n

    # Hessian inverse with MASS::ginv() fallback
    H_gamma_inv <- tryCatch(
      solve(H_gamma),
      error = function(e) {
        if (requireNamespace("MASS", quietly = TRUE)) {
          MASS::ginv(H_gamma)
        } else {
          stop_lwdid("Logit Hessian singular and MASS package unavailable.",
                     class = "lwdid_estimation_failed")
        }
      }
    )

    # ATT sensitivity to PS coefficients: dATT/dgamma
    X_ctrl_const <- X_const[control_mask, , drop = FALSE]
    Y_ctrl <- Y[control_mask]
    dATT_dgamma <- -colSums(w_ctrl * X_ctrl_const * Y_ctrl) / (n * p_bar)

    # Adjustment term
    adjustment <- as.numeric((S_gamma %*% t(H_gamma_inv)) %*% dATT_dgamma)

    # Complete IF
    psi <- psi_ht - adjustment

    # Variance: var(psi)/n with N-1 degrees of freedom
    var_att <- var(psi) / n
    se <- sqrt(var_att)

  } else {
    # ---- Bootstrap SE ----
    if (!is.null(seed)) set.seed(seed)
    boot_atts <- numeric(boot_reps)
    for (b in seq_len(boot_reps)) {
      idx <- sample.int(n_total, n_total, replace = TRUE)
      boot_atts[b] <- .ipw_point_estimate(
        data = data[idx, , drop = FALSE],
        y = y, d = d,
        propensity_controls = propensity_controls,
        trim_threshold = trim_threshold,
        trim_method = trim_method
      )
    }
    se <- sd(boot_atts, na.rm = TRUE)
  }

  # ---- Step 10: Normal distribution inference ----
  if (se > 0) {
    t_stat <- att / se
    pvalue <- 2 * (1 - pnorm(abs(t_stat)))
  } else {
    t_stat <- NaN
    pvalue <- NaN
  }
  z_crit <- qnorm(1 - alpha / 2)
  ci_lower <- att - z_crit * se
  ci_upper <- att + z_crit * se

  # ---- Step 11: Weight CV diagnostic ----
  w_ctrl_norm <- weights_valid[which(control_mask)]
  if (mean(w_ctrl_norm) > 0) {
    weights_cv <- sd(w_ctrl_norm) / mean(w_ctrl_norm)
    if (weights_cv > 2.0) {
      warn_lwdid(
        sprintf("IPW weight CV = %.2f > 2.0. Possible overlap violation.", weights_cv),
        class = "lwdid_overlap", detail = "ipw_extreme_weights")
    }
  }

  # ---- Step 12: Degrees of freedom (interface compatibility) ----
  df_resid <- max(1L, n1 + n0 - 2L)

  # ---- Step 13: Full-sample weight vector ----
  weights_all <- numeric(n_total)
  weights_all[which(valid_mask)] <- weights_valid

  # ---- Step 14: Return 15-field list ----
  list(
    att                   = att,
    se                    = se,
    ci_lower              = ci_lower,
    ci_upper              = ci_upper,
    t_stat                = t_stat,
    pvalue                = pvalue,
    propensity_scores     = ps_raw,
    weights               = weights_all,
    propensity_model_coef = ps_result$coefficients,
    n_treated             = as.integer(n1),
    n_control             = as.integer(n0),
    n_trimmed             = as.integer(n_trimmed),
    df_resid              = as.integer(df_resid),
    vce_method            = vce_method,
    estimator             = "ipw"
  )
}


#' Bootstrap helper: IPW point estimate only
#'
#' @description Minimal IPW-ATT computation for bootstrap iterations.
#'   Executes full pipeline: glm -> predict -> trim -> weights -> Hajek ATT.
#'   Returns NA_real_ on any failure.
#'
#' @keywords internal
.ipw_point_estimate <- function(data, y, d, propensity_controls,
                                trim_threshold, trim_method = "clip") {
  Y <- as.numeric(data[[y]])
  D <- as.integer(data[[d]])

  # PS estimation (tryCatch: return NA on failure)
  fml <- as.formula(paste(d, "~", paste(propensity_controls, collapse = " + ")))
  glm_fit <- tryCatch(
    glm(fml, data = data, family = binomial(link = "logit")),
    error = function(e) return(NULL),
    warning = function(w) {
      # Suppress glm warnings (convergence etc.) in bootstrap
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(glm_fit)) return(NA_real_)

  ps <- predict(glm_fit, type = "response")

  # Trimming
  if (trim_method == "clip") {
    ps <- pmax(trim_threshold, pmin(1 - trim_threshold, ps))
    valid <- rep(TRUE, length(ps))
  } else {
    valid <- ps >= trim_threshold & ps <= (1 - trim_threshold)
  }
  Y <- Y[valid]; D <- D[valid]; ps <- ps[valid]

  # Hajek ATT
  treat <- D == 1L; ctrl <- D == 0L
  if (sum(treat) == 0L || sum(ctrl) == 0L) return(NA_real_)
  w <- ps[ctrl] / (1 - ps[ctrl])
  w_sum <- sum(w)
  if (w_sum <= 0) return(NA_real_)

  mean(Y[treat]) - sum(w * Y[ctrl]) / w_sum
}
