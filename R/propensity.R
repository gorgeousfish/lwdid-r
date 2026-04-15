# ============================================================================
# propensity.R — Propensity Score & Outcome Model Shared Infrastructure
# ============================================================================
# Implements estimate_propensity_score() and estimate_outcome_model()
# for IPW (Story 6.2), IPWRA (Story 6.3), and PSM (Story 6.4).
#
# Paper references:
#   - lw2025 Procedure 3.1 Step 2: apply standard TE methods
#   - lw2025 Procedure 4.1 Step 3: apply to (g,r) subsamples
#   - Wooldridge (2025b §19.4): WLS for IPWRA ATT
#   - Wooldridge (2007): IPWRA doubly-robust estimator
# ============================================================================

#' Estimate propensity scores via logistic regression
#'
#' @description Uses unpenalized logistic regression (MLE via IRLS) to estimate
#'   P(D=1|X). Returns propensity scores, coefficients, and the glm object
#'   for downstream IF-based SE computation.
#'
#' @param data data.frame or data.table containing treatment indicator and controls
#' @param d character, treatment indicator variable name (1=treated, 0=control)
#' @param controls character vector, control variable names (length >= 1)
#' @param trim_threshold numeric, trimming threshold in (0, 0.5), default 0.01
#' @param trim_method character, "clip" (default) or "drop"
#'
#' @return list with 8 elements:
#'   \describe{
#'     \item{propensity_scores}{numeric vector (clip: clipped PS; drop: PS with NA)}
#'     \item{coefficients}{named list with "_intercept" key}
#'     \item{converged}{logical, whether glm converged}
#'     \item{n_trimmed}{integer, number of trimmed observations}
#'     \item{trimmed_mask}{logical vector (meaningful only for drop mode)}
#'     \item{trim_method}{character, actual trim method used}
#'     \item{glm_object}{glm object for downstream vcov() calls}
#'     \item{removed_cols}{character vector of removed constant columns}
#'   }
#'
#' @references
#'   lw2025 Procedure 3.1 Step 2;
#'   MLE for \eqn{P(D=1 \mid X)} under a Bernoulli log-likelihood.
#'
#' @keywords internal
estimate_propensity_score <- function(data, d, controls,
                                      trim_threshold = 0.01,
                                      trim_method = c("clip", "drop")) {

  # --- Step 0: Parameter validation ---
  trim_method <- match.arg(trim_method)

  if (!is.character(d) || length(d) != 1L) {
    stop_lwdid("'d' must be a single character string.",
               class = "lwdid_invalid_param",
               param = "d", value = d)
  }
  if (!is.character(controls) || length(controls) < 1L) {
    stop_lwdid("'controls' must be a character vector of length >= 1.",
               class = "lwdid_invalid_param",
               param = "controls", value = controls)
  }
  if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
      trim_threshold <= 0 || trim_threshold >= 0.5) {
    stop_lwdid(
      sprintf("'trim_threshold' must be in (0, 0.5), got %s.", trim_threshold),
      class = "lwdid_invalid_param",
      param = "trim_threshold", value = trim_threshold)
  }

  # --- Step 1: Column existence ---
  if (!d %in% names(data)) {
    stop_lwdid(sprintf("Column '%s' not found in data.", d),
               class = "lwdid_invalid_param", param = "d")
  }
  missing_ctrl <- setdiff(controls, names(data))
  if (length(missing_ctrl) > 0L) {
    stop_lwdid(
      sprintf("Control columns not found in data: %s",
              paste(missing_ctrl, collapse = ", ")),
      class = "lwdid_invalid_param", param = "controls")
  }

  # --- Step 2: NA checks ---
  D <- data[[d]]
  if (anyNA(D)) {
    stop_lwdid(
      sprintf("Treatment indicator '%s' contains NA values.", d),
      class = "lwdid_invalid_param", param = d)
  }
  na_cols <- vapply(controls, function(v) anyNA(data[[v]]), logical(1))
  if (any(na_cols)) {
    stop_lwdid(
      sprintf("Control variable(s) contain NA: %s",
              paste(controls[na_cols], collapse = ", ")),
      class = "lwdid_invalid_param", param = "controls")
  }

  # --- Step 3: Binary check ---
  D_num <- as.numeric(D)
  if (!all(D_num %in% c(0, 1))) {
    stop_lwdid(
      sprintf("Treatment indicator '%s' must be binary (0/1).", d),
      class = "lwdid_invalid_param", param = d)
  }

  # --- Step 4: Group existence ---
  n_treated <- sum(D_num == 1)
  n_control <- sum(D_num == 0)
  if (n_treated == 0L) {
    stop_lwdid("No treated units (all D=0), cannot estimate propensity score.",
               class = "lwdid_insufficient_data")
  }
  if (n_control == 0L) {
    stop_lwdid("No control units (all D=1), cannot estimate propensity score.",
               class = "lwdid_insufficient_data")
  }

  # --- Step 5: Constant column detection ---
  X <- as.matrix(data[, controls, drop = FALSE])
  col_sd <- apply(X, 2, sd, na.rm = TRUE)
  constant_mask <- col_sd == 0 | is.na(col_sd)
  removed_cols <- controls[constant_mask]

  if (length(removed_cols) > 0L) {
    warn_lwdid(
      sprintf("Removed %d constant column(s) from PS model: %s",
              length(removed_cols), paste(removed_cols, collapse = ", ")),
      class = "lwdid_data")
    controls <- setdiff(controls, removed_cols)
  }

  if (length(controls) == 0L) {
    stop_lwdid(
      "No non-constant control variables remaining after removal.",
      class = "lwdid_insufficient_data")
  }

  # --- Step 6: Fit glm ---
  fml <- as.formula(paste(d, "~", paste(controls, collapse = " + ")))
  glm_nonconvergence_warning <- FALSE
  glm_boundary_warning <- FALSE
  glm_fit <- withCallingHandlers(
    tryCatch(
      stats::glm(fml, data = data, family = stats::binomial(link = "logit")),
      error = function(e) {
        stop_lwdid(
          sprintf("glm fitting failed: %s", conditionMessage(e)),
          class = "lwdid_estimation_failed")
      }
    ),
    warning = function(w) {
      warning_msg <- conditionMessage(w)

      if (grepl("did not converge", warning_msg, fixed = TRUE)) {
        glm_nonconvergence_warning <<- TRUE
        invokeRestart("muffleWarning")
      }

      if (grepl("fitted probabilities numerically 0 or 1 occurred",
                warning_msg, fixed = TRUE)) {
        glm_boundary_warning <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )

  # --- Step 7: Convergence checks ---
  converged <- glm_fit$converged
  if (!converged || glm_nonconvergence_warning) {
    warn_lwdid("glm did not converge. PS coefficients may be unreliable.",
               class = "lwdid_numerical", detail = "ps_not_converged")
  }
  if (isTRUE(glm_fit$boundary) || glm_boundary_warning) {
    warn_lwdid(
      "glm boundary convergence: possible complete/quasi-complete separation.",
      class = "lwdid_numerical", detail = "ps_boundary")
  }

  # --- Step 8: Extract PS ---
  ps_raw <- as.numeric(stats::predict(glm_fit, type = "response"))

  # --- Step 9: Trimming ---
  if (trim_method == "clip") {
    ps_trimmed <- pmax(trim_threshold, pmin(1 - trim_threshold, ps_raw))
    n_trimmed <- sum(ps_raw != ps_trimmed)
    trimmed_mask <- rep(FALSE, length(ps_raw))
  } else {
    # drop mode
    trimmed_mask <- ps_raw < trim_threshold | ps_raw > (1 - trim_threshold)
    ps_trimmed <- ps_raw
    ps_trimmed[trimmed_mask] <- NA_real_
    n_trimmed <- sum(trimmed_mask)
  }

  # --- Step 10: Extract coefficients ---
  coef_vec <- stats::coef(glm_fit)
  coef_list <- as.list(coef_vec)
  intercept_idx <- which(names(coef_list) == "(Intercept)")
  if (length(intercept_idx) > 0L) {
    names(coef_list)[intercept_idx] <- "_intercept"
  }

  # --- Step 11: Return ---
  list(
    propensity_scores = ps_trimmed,
    coefficients      = coef_list,
    converged         = converged,
    n_trimmed         = as.integer(n_trimmed),
    trimmed_mask      = trimmed_mask,
    trim_method       = trim_method,
    glm_object        = glm_fit,
    removed_cols      = removed_cols
  )
}


#' Estimate outcome model (control group conditional expectation)
#'
#' @description Fits a linear regression for the control-group conditional mean
#'   \eqn{E[Y \mid X, D = 0]},
#'   generates predictions for all units. Supports OLS (default) and WLS
#'   (when sample_weights provided).
#'
#' @param data data.frame or data.table
#' @param y character, outcome variable name
#' @param d character, treatment indicator variable name
#' @param controls character vector, control variable names
#' @param sample_weights numeric vector or NULL, sample weights (length = nrow(data)).
#'   For IPWRA ATT: caller passes w_i = p_hat_i / (1 - p_hat_i).
#'   Internally extracts control group subset for WLS fitting.
#'
#' @return list with 2 elements:
#'   \describe{
#'     \item{predictions}{numeric vector, m0_hat for all units (length = nrow(data))}
#'     \item{coefficients}{named list with "_intercept" key}
#'   }
#'
#' @references
#'   Wooldridge (2025b Section 19.4): WLS for IPWRA ATT;
#'   OLS: beta_hat = (X0'X0)^{-1} X0'Y0;
#'   WLS: beta_hat = (X0'WX0)^{-1} X0'WY0
#'
#' @keywords internal
estimate_outcome_model <- function(data, y, d, controls,
                                   sample_weights = NULL) {

  # --- Step 1: Column existence ---
  for (col_name in c(y, d)) {
    if (!col_name %in% names(data)) {
      stop_lwdid(sprintf("Column '%s' not found in data.", col_name),
                 class = "lwdid_invalid_param", param = col_name)
    }
  }
  missing_ctrl <- setdiff(controls, names(data))
  if (length(missing_ctrl) > 0L) {
    stop_lwdid(
      sprintf("Control columns not found in data: %s",
              paste(missing_ctrl, collapse = ", ")),
      class = "lwdid_invalid_param", param = "controls")
  }

  # --- Step 2: NA checks ---
  for (col_name in c(y, d)) {
    if (anyNA(data[[col_name]])) {
      stop_lwdid(
        sprintf("Column '%s' contains NA values.", col_name),
        class = "lwdid_invalid_param", param = col_name)
    }
  }
  na_cols <- vapply(controls, function(v) anyNA(data[[v]]), logical(1))
  if (any(na_cols)) {
    stop_lwdid(
      sprintf("Control variable(s) contain NA: %s",
              paste(controls[na_cols], collapse = ", ")),
      class = "lwdid_invalid_param", param = "controls")
  }

  # --- Step 3: Extract data ---
  D <- as.numeric(data[[d]])
  Y <- as.numeric(data[[y]])
  X <- as.matrix(data[, controls, drop = FALSE])
  ctrl_mask <- D == 0
  X_ctrl <- X[ctrl_mask, , drop = FALSE]
  Y_ctrl <- Y[ctrl_mask]
  X_ctrl_const <- cbind(1, X_ctrl)

  # --- Step 4: Fit model ---
  if (!is.null(sample_weights)) {
    # --- Step 4a: Weight validation ---
    if (length(sample_weights) != nrow(data)) {
      stop_lwdid(
        sprintf("sample_weights length (%d) != nrow(data) (%d).",
                length(sample_weights), nrow(data)),
        class = "lwdid_invalid_param", param = "sample_weights")
    }
    if (any(!is.finite(sample_weights))) {
      stop_lwdid("sample_weights contains Inf or NaN values.",
                 class = "lwdid_invalid_param", param = "sample_weights")
    }
    w_ctrl <- sample_weights[ctrl_mask]
    if (mean(w_ctrl) <= 0) {
      stop_lwdid(
        sprintf("Control group weight mean <= 0 (got %.6e).", mean(w_ctrl)),
        class = "lwdid_invalid_param", param = "sample_weights")
    }

    # --- WLS mode ---
    w_ctrl <- pmax(w_ctrl, 1e-10)  # weight floor
    sqrt_w <- sqrt(w_ctrl)
    X_weighted <- X_ctrl_const * sqrt_w  # R recycling: row-wise multiply
    Y_weighted <- Y_ctrl * sqrt_w
    XtWX <- crossprod(X_weighted)
    XtWY <- crossprod(X_weighted, Y_weighted)
    beta <- tryCatch(
      solve(XtWX, XtWY),
      error = function(e) {
        stop_lwdid(
          sprintf("WLS design matrix singular: %s", conditionMessage(e)),
          class = "lwdid_estimation_failed")
      }
    )
  } else {
    # --- OLS mode ---
    XtX <- crossprod(X_ctrl_const)
    XtY <- crossprod(X_ctrl_const, Y_ctrl)
    beta <- tryCatch(
      solve(XtX, XtY),
      error = function(e) {
        stop_lwdid(
          sprintf("OLS design matrix singular: %s", conditionMessage(e)),
          class = "lwdid_estimation_failed")
      }
    )
  }

  # --- Step 5: Full-sample prediction ---
  X_all_const <- cbind(1, X)
  m0_hat <- as.numeric(X_all_const %*% beta)

  # --- Step 6: Extract coefficients ---
  coef_list <- list()
  coef_list[["_intercept"]] <- as.numeric(beta[1])
  for (j in seq_along(controls)) {
    coef_list[[controls[j]]] <- as.numeric(beta[j + 1L])
  }

  # --- Step 7: Return ---
  list(
    predictions  = m0_hat,
    coefficients = coef_list
  )
}
