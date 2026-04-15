# ============================================================================
# lwdid: Estimation Functions
# ============================================================================
# This file contains estimation functions for the lwdid package:
#   - estimate_ra_common():      RA estimator (lm-based, full VCE support)
#   - estimate_ra_common_fast(): RA estimator (QR-based, ATT only)
#   - estimate_period_effects(): Period-specific ATT estimation
#
# Reference: Lee & Wooldridge (2025, 2026) [lw2026]
# ============================================================================

#' @title Common Timing RA Estimator with VCE Support
#'
#' @description
#' Regression Adjustment estimator for Common Timing designs
#' (lw2026 equations 2.18--2.19) using \code{lm()} internally for
#' sandwich-based VCE compatibility. Implements a three-tier fallback
#' strategy for controls: (1) full interaction, (2) simple controls,
#' (3) no controls, depending on sample size relative to the number
#' of covariates.
#'
#' @param y_trans Numeric vector, transformed outcome (length N).
#' @param d Integer vector, treatment indicator (length N).
#' @param x Numeric matrix or NULL, controls (N x K).
#' @param vce Character or NULL. VCE type: \code{NULL}
#'   (homoskedastic, default), \code{"hc0"}--\code{"hc4"}
#'   (heteroskedasticity-robust), \code{"robust"} (alias for
#'   \code{"hc1"}), \code{"cluster"} (cluster-robust).
#'   Case-insensitive.
#' @param cluster Numeric/character vector or NULL. Cluster
#'   identifiers (length N), required when \code{vce = "cluster"}.
#' @param alpha Numeric in (0,1). Significance level for confidence
#'   intervals. Default 0.05.
#'
#' @return Named list with elements:
#'   \describe{
#'     \item{att}{Numeric, estimated ATT.}
#'     \item{se}{Numeric, standard error of ATT.}
#'     \item{t_stat}{Numeric, t-statistic.}
#'     \item{df}{Integer, degrees of freedom.}
#'     \item{pvalue}{Numeric, two-sided p-value.}
#'     \item{ci_lower}{Numeric, lower confidence bound.}
#'     \item{ci_upper}{Numeric, upper confidence bound.}
#'     \item{params}{Named numeric, all regression coefficients.}
#'     \item{vcov}{Matrix, variance-covariance matrix.}
#'     \item{resid}{Numeric vector, OLS residuals.}
#'     \item{fit}{lm object, fitted model.}
#'     \item{n}{Integer, total sample size.}
#'     \item{n_treated}{Integer, number of treated units.}
#'     \item{n_control}{Integer, number of control units.}
#'     \item{K}{Integer, number of control variables used.}
#'     \item{vce_type}{Character, VCE type actually applied.}
#'     \item{n_clusters}{Integer or NULL, number of clusters.}
#'     \item{controls_tier}{Character, tier selected by fallback.}
#'     \item{X_design}{Matrix, design matrix used in estimation.}
#'     \item{df_resid}{Integer, OLS residual degrees of freedom (n - k).}
#'     \item{df_inference}{Integer, inference degrees of freedom (same as
#'       df_resid for homoskedastic/HC; G_cluster - 1 for cluster).}
#'   }
#'
#' @references Lee, S. and Wooldridge, J.M. (2025, 2026).
#'   Equations 2.18--2.19.
#'
#' @keywords internal
estimate_ra_common <- function(y_trans, d, x = NULL,
                               vce = NULL, cluster = NULL,
                               alpha = 0.05) {
  # Step 0: Sample size validation
  n <- length(y_trans)
  n1 <- sum(d == 1)
  n0 <- sum(d == 0)
  if (n1 == 0L) {
    stop_lwdid(
      "No treated units in regression sample (N_1=0), cannot estimate ATT",
      class = "lwdid_insufficient_data",
      n = n, n_treated = 0L, n_control = n0)
  }
  if (n0 == 0L) {
    stop_lwdid(
      "No control units in regression sample (N_0=0), cannot estimate ATT",
      class = "lwdid_insufficient_data",
      n = n, n_treated = n1, n_control = 0L)
  }
  if (n < 3L) {
    stop_lwdid(
      sprintf("Insufficient sample size (N=%d, need >=3), df too low for inference", n),
      class = "lwdid_insufficient_data",
      n = n, n_treated = n1, n_control = n0)
  }

  # Step 1: Three-tier fallback - build data.frame + lm() + X_design
  controls_tier <- "none"
  K <- 0L

  if (is.null(x) || ncol(x) == 0) {
    reg_df <- data.frame(y = y_trans, D = d)
    X_design <- cbind(1, d)
    fit <- lm(y ~ ., data = reg_df)
    K <- 0L
    controls_tier <- "none"
  } else {
    K <- ncol(x)
    if (n1 > K + 1L && n0 > K + 1L) {
      # Tier 1: Full interaction (lw2026 eq. 2.18 post-text)
      x_bar1 <- colMeans(x[d == 1, , drop = FALSE])
      x_c <- sweep(x, 2, x_bar1, "-")
      xc_names <- paste0("Xc", seq_len(K))
      dxc_names <- paste0("DXc", seq_len(K))
      reg_df <- data.frame(y = y_trans, D = d)
      reg_df[xc_names] <- as.data.frame(x_c)
      reg_df[dxc_names] <- as.data.frame(d * x_c)
      X_design <- cbind(1, d, x_c, d * x_c)
      fit <- lm(y ~ ., data = reg_df)
      controls_tier <- "full_interaction"
    } else if (n > K + 2L) {
      # Tier 2: Simple controls (lw2026 eq. 2.18)
      x_names <- paste0("X", seq_len(K))
      reg_df <- data.frame(y = y_trans, D = d)
      reg_df[x_names] <- as.data.frame(x)
      X_design <- cbind(1, d, x)
      fit <- lm(y ~ ., data = reg_df)
      controls_tier <- "simple"
      warn_lwdid(
        sprintf(
          paste0("Insufficient sample for full interaction model ",
                 "(need N_1>%d and N_0>%d, got N_1=%d, N_0=%d), ",
                 "degraded to simple controls (lw2026 eq. 2.18)"),
          K + 1L, K + 1L, n1, n0),
        class = "lwdid_data",
        detail = "controls_degraded_to_simple",
        action_taken = "interaction terms dropped")
    } else {
      # Tier 3: Drop all controls
      warn_lwdid(
        sprintf(
          paste0("Insufficient sample for any controls model ",
                 "(need N>%d, got N=%d), all controls dropped"),
          K + 2L, n),
        class = "lwdid_data",
        detail = "controls_dropped",
        action_taken = "all controls dropped")
      reg_df <- data.frame(y = y_trans, D = d)
      X_design <- cbind(1, d)
      fit <- lm(y ~ ., data = reg_df)
      K <- 0L
      controls_tier <- "none"
    }
  }

  # Step 1b: Rank-deficiency check
  fit_coefs <- coef(fit)
  if (any(is.na(fit_coefs))) {
    fit_rank <- fit$rank
    expected_rank <- length(fit_coefs)
    stop_lwdid(
      sprintf("Design matrix rank-deficient (rank=%d, need=%d)",
              fit_rank, expected_rank),
      class = "lwdid_singular_design",
      rank = fit_rank, expected_rank = expected_rank)
  }

  # Step 2: Extract ATT and D position
  att <- unname(fit_coefs[["D"]])
  d_idx <- which(names(fit_coefs) == "D")

  # Step 3: VCE
  vce_result <- compute_vce(fit, vce = vce, cluster = cluster)

  # Step 4: Inference (with degenerate SE check)
  var_d <- vce_result$vcov[d_idx, d_idx]
  se_d <- sqrt(var_d)

  if (!is.finite(se_d) || se_d < 1e-10) {
    warn_lwdid(
      paste0("Standard error is effectively zero (SE < 1e-10). ",
             "Inference results set to NA. Possible causes: ",
             "perfect fit, zero residual variance, or degenerate data."),
      class = "lwdid_small_sample",
      detail = "degenerate_se",
      action_taken = "inference results set to NA")
    inf <- list(att = att, se = NA_real_, t_stat = NA_real_,
                df = vce_result$df, pvalue = NA_real_,
                ci_lower = NA_real_, ci_upper = NA_real_)
  } else {
    inf <- compute_inference(
      att = att, vcov_mat = vce_result$vcov,
      df = vce_result$df, alpha = alpha, coef_index = d_idx)
  }

  # Step 5: Return
  # Cross-Epic dependency: Story E5-02 requires separate df_resid and
  # df_inference fields for cohort effect aggregation.
  #   df_resid: OLS residual df (n - k), independent of VCE type.
  #   df_inference: inference df for t-distribution.
  #     Homoskedastic/HC: same as df_resid
  #     Cluster: G_cluster - 1
  # The original `df` field (from compute_inference) is retained for
  # backward compatibility.
  c(inf, list(
    params        = fit_coefs,
    vcov          = vce_result$vcov,
    resid         = as.numeric(residuals(fit)),
    fit           = fit,
    n             = n,
    n_treated     = n1,
    n_control     = n0,
    K             = K,
    vce_type      = vce_result$vce_type,
    n_clusters    = vce_result$n_clusters,
    controls_tier = controls_tier,
    X_design      = X_design,
    df_resid      = fit$df.residual,
    df_inference  = vce_result$df))
}


# ============================================================================
# Fast QR Path (ATT only, for RI/Bootstrap inner loops)
# ============================================================================

#' @title Common Timing RA Estimator -- Fast QR Path
#'
#' @description
#' QR-based fast path for ATT-only estimation. Designed for Epic 7
#' (Randomization Inference) and Epic 9 (Wild Cluster Bootstrap)
#' inner loops where \code{lm()} overhead is prohibitive.
#' Numerically equivalent to
#' \code{estimate_ra_common(vce = NULL)} for ATT (< 1e-12).
#'
#' @param y_trans Numeric vector, transformed outcome (length N).
#' @param d Integer vector, treatment indicator (length N).
#' @param x Numeric matrix or NULL, controls. Default NULL.
#' @return Named list: att, n, n_treated, n_control, K,
#'   controls_tier.
#' @keywords internal
estimate_ra_common_fast <- function(y_trans, d, x = NULL) {
  n <- length(y_trans)
  n1 <- sum(d == 1)
  n0 <- sum(d == 0)
  if (n1 == 0L) {
    stop_lwdid(
      "No treated units in regression sample (N_1=0), cannot estimate ATT",
      class = "lwdid_insufficient_data",
      n = n, n_treated = 0L, n_control = n0)
  }
  if (n0 == 0L) {
    stop_lwdid(
      "No control units in regression sample (N_0=0), cannot estimate ATT",
      class = "lwdid_insufficient_data",
      n = n, n_treated = n1, n_control = 0L)
  }
  if (n < 3L) {
    stop_lwdid(
      sprintf("Insufficient sample size (N=%d, need >=3)", n),
      class = "lwdid_insufficient_data",
      n = n, n_treated = n1, n_control = n0)
  }

  controls_tier <- "none"
  K <- 0L
  if (is.null(x) || ncol(x) == 0) {
    X_design <- cbind(1, d)
    K <- 0L
    controls_tier <- "none"
  } else {
    K <- ncol(x)
    if (n1 > K + 1L && n0 > K + 1L) {
      x_bar1 <- colMeans(x[d == 1, , drop = FALSE])
      x_c <- sweep(x, 2, x_bar1, "-")
      X_design <- cbind(1, d, x_c, d * x_c)
      controls_tier <- "full_interaction"
    } else if (n > K + 2L) {
      X_design <- cbind(1, d, x)
      controls_tier <- "simple"
    } else {
      X_design <- cbind(1, d)
      K <- 0L
      controls_tier <- "none"
    }
  }

  qr_fit <- qr(X_design)
  if (qr_fit$rank < ncol(X_design)) {
    stop_lwdid(
      sprintf("Design matrix rank-deficient (rank=%d, need=%d)",
              qr_fit$rank, ncol(X_design)),
      class = "lwdid_singular_design",
      rank = qr_fit$rank, expected_rank = ncol(X_design))
  }
  beta <- qr.coef(qr_fit, y_trans)
  att <- unname(beta[2])

  list(att = att, n = n, n_treated = n1, n_control = n0,
       K = K, controls_tier = controls_tier)
}


# ============================================================================
# Period-specific effect estimation (Story E2-05)
# ============================================================================

#' @title Construct all-NA period row (15 columns)
#' @param r Numeric, period value.
#' @return data.frame with 1 row and 15 columns.
#' @keywords internal
.make_na_period_row <- function(r) {
  data.frame(
    tindex = r, period = r,
    att = NA_real_, se = NA_real_,
    t_stat = NA_real_, pvalue = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    n_obs = 0L, n_treated = 0L, n_control = 0L,
    df = NA_integer_,
    vce_type = NA_character_, n_clusters = NA_integer_,
    controls_tier = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' @title Construct degenerate period row (15 columns)
#' @param r Numeric, period value.
#' @param est List, return value from estimate_ra_common().
#' @return data.frame with 1 row and 15 columns.
#' @keywords internal
.make_degenerate_period_row <- function(r, est) {
  data.frame(
    tindex = r, period = r,
    att = NA_real_, se = NA_real_,
    t_stat = NA_real_, pvalue = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    n_obs = est$n, n_treated = est$n_treated, n_control = est$n_control,
    df = NA_integer_,
    vce_type = est$vce_type %||% NA_character_,
    n_clusters = if (is.null(est$n_clusters)) NA_integer_ else est$n_clusters,
    controls_tier = est$controls_tier,
    stringsAsFactors = FALSE
  )
}

#' @title Period-specific effect estimation with VCE support
#'
#' @description
#' For each post-treatment period r, independently calls
#' \code{\link{estimate_ra_common}} to estimate the period-specific
#' ATT (lw2026 equation 2.20, Procedure 2.1).
#'
#' @param dt_transformed data.table, panel data after transformation.
#' @param y_trans_col Character, name of transformed outcome column.
#' @param d_col Character, name of treatment indicator column.
#' @param tvar Character, name of time variable column.
#' @param x Character vector or NULL, control variable column names.
#' @param periods Numeric vector, post-treatment period values.
#' @param vce Character or NULL. VCE type: \code{NULL}
#'   (homoskedastic, default), \code{"hc0"}--\code{"hc4"}
#'   (heteroskedasticity-robust), \code{"robust"} (alias for
#'   \code{"hc1"}), \code{"cluster"} (cluster-robust).
#'   Case-insensitive.
#' @param cluster_var Character or NULL. Column name of cluster
#'   variable in \code{dt_transformed}. Required when
#'   \code{vce = "cluster"}. Note: unlike
#'   \code{\link{estimate_ra_common}}'s \code{cluster} parameter
#'   (which takes a vector of values), this parameter takes a
#'   column name.
#' @param alpha Numeric in (0,1). Significance level for confidence
#'   intervals. Default 0.05.
#' @param include_pretreatment Logical, stub for Epic 8.
#'   Default FALSE.
#'
#' @return A \code{data.frame} with one row per period and 15
#'   columns: \code{tindex}, \code{period}, \code{att}, \code{se},
#'   \code{t_stat}, \code{pvalue}, \code{ci_lower},
#'   \code{ci_upper}, \code{n_obs}, \code{n_treated},
#'   \code{n_control}, \code{df}, \code{vce_type},
#'   \code{n_clusters}, \code{controls_tier}.
#'
#' @references Lee, S. and Wooldridge, J.M. (2025, 2026).
#'   Equation 2.20, Procedure 2.1.
#'
#' @keywords internal
estimate_period_effects <- function(dt_transformed, y_trans_col, d_col,
                                    tvar, x = NULL, periods,
                                    vce = NULL, cluster_var = NULL,
                                    alpha = 0.05,
                                    include_pretreatment = FALSE) {

  results <- lapply(periods, function(r) {
    # Step 1: Extract cross-section for period r
    sub <- dt_transformed[get(tvar) == r]
    y_r <- sub[[y_trans_col]]
    d_r <- sub[[d_col]]
    x_r <- if (!is.null(x)) as.matrix(sub[, x, with = FALSE]) else NULL
    cl_r <- if (!is.null(cluster_var)) sub[[cluster_var]] else NULL

    # Step 2: NA filtering (qr() does not handle NA)
    valid <- !is.na(y_r) & !is.na(d_r)
    if (!is.null(x_r)) {
      valid <- valid & complete.cases(x_r)
    }
    if (!is.null(cl_r)) {
      valid <- valid & !is.na(cl_r)
    }

    # Step 3: No valid observations -> all-NA row
    if (sum(valid) == 0L) {
      return(.make_na_period_row(r))
    }

    y_r <- y_r[valid]
    d_r <- d_r[valid]
    if (!is.null(x_r)) x_r <- x_r[valid, , drop = FALSE]
    if (!is.null(cl_r)) cl_r <- cl_r[valid]

    # Step 4: Call RA estimator (three-tier fallback automatic)
    tryCatch({
      est <- estimate_ra_common(y_r, d_r, x_r,
                                vce = vce, cluster = cl_r,
                                alpha = alpha)

      # Step 5: Degenerate SE check
      if (is.na(est$se) || !is.finite(est$se) || est$se == 0) {
        warn_lwdid(
          sprintf(
            "Period %s regression produced degenerate result (se=%s), set to NA",
            r, est$se),
          class = "lwdid_small_sample",
          detail = "degenerate_period_regression")
        return(.make_degenerate_period_row(r, est))
      }

      # Step 6: Normal result (15 columns)
      data.frame(
        tindex = r, period = r,
        att = est$att, se = est$se,
        t_stat = est$t_stat, pvalue = est$pvalue,
        ci_lower = est$ci_lower, ci_upper = est$ci_upper,
        n_obs = est$n, n_treated = est$n_treated, n_control = est$n_control,
        df = est$df,
        vce_type = est$vce_type %||% NA_character_,
        n_clusters = if (is.null(est$n_clusters)) NA_integer_ else est$n_clusters,
        controls_tier = est$controls_tier,
        stringsAsFactors = FALSE
      )
    },
    lwdid_error = function(e) {
      # Layer 1: lwdid errors (insufficient data, singular design, etc.)
      warn_lwdid(
        sprintf("Period %s estimation failed: %s", r, conditionMessage(e)),
        class = "lwdid_small_sample",
        detail = "period_regression_failed")
      data.frame(
        tindex = r, period = r,
        att = NA_real_, se = NA_real_,
        t_stat = NA_real_, pvalue = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        n_obs = length(y_r), n_treated = sum(d_r == 1L),
        n_control = sum(d_r == 0L),
        df = NA_integer_,
        vce_type = NA_character_, n_clusters = NA_integer_,
        controls_tier = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      # Layer 2: Generic errors (sandwich internal, numerical singularity)
      warn_lwdid(
        sprintf("Period %s unexpected error: %s", r, conditionMessage(e)),
        class = "lwdid_numerical",
        detail = "period_regression_failed")
      data.frame(
        tindex = r, period = r,
        att = NA_real_, se = NA_real_,
        t_stat = NA_real_, pvalue = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_,
        n_obs = length(y_r), n_treated = sum(d_r == 1L),
        n_control = sum(d_r == 0L),
        df = NA_integer_,
        vce_type = NA_character_, n_clusters = NA_integer_,
        controls_tier = NA_character_,
        stringsAsFactors = FALSE
      )
    })
  })

  do.call(rbind, results)
}
