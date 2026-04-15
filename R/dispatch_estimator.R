# ============================================================================
# dispatch_estimator.R — Unified Estimator Routing
# ============================================================================
# Paper: lw2025 Procedure 3.1 Step 2 / Procedure 4.1 Step 3
# ============================================================================
#' Dispatch to Estimator Function
#'
#' @description
#' Unified routing function that dispatches estimation to one of four
#' treatment effect estimators: Regression Adjustment (RA), Inverse
#' Probability Weighting (IPW), Inverse Probability Weighted Regression
#' Adjustment (IPWRA), or Propensity Score Matching (PSM).
#'
#' Implements the estimator selection mechanism described in
#' Lee & Wooldridge (2025) Procedure 3.1 Step 2: "apply standard TE
#' methods -- such as linear RA, IPW, IPWRA, matching -- to the
#' cross section". Called by \code{\link{lwdid}} for both Common Timing
#' and Staggered adoption designs (Procedure 4.1 Step 3).
#'
#' @param data data.frame. The (transformed) cross-sectional data for
#'   a single estimation unit.
#' @param y character(1). Name of the outcome variable column.
#' @param d character(1). Name of the treatment indicator column.
#' @param controls character vector or NULL. Names of outcome model
#'   control variables. Required for IPWRA; required for IPW/PSM if
#'   \code{ps_controls} is NULL.
#' @param ps_controls character vector or NULL. Names of propensity
#'   score model controls. If NULL, \code{controls} is used.
#' @param estimator character(1). One of \code{"ra"}, \code{"ipw"},
#'   \code{"ipwra"}, \code{"psm"}. Default \code{"ra"}.
#' @param vce character(1) or NULL. Variance-covariance estimator type
#'   passed to RA estimator.
#' @param cluster_var character(1) or NULL. Cluster variable for RA.
#' @param alpha numeric(1). Significance level in (0, 1). Default 0.05.
#' @param trim_threshold numeric(1). Propensity score trimming threshold
#'   in (0, 0.5). Default 0.01. Used by IPW, IPWRA, PSM.
#' @param trim_method character(1). One of \code{"clip"}, \code{"drop"}.
#'   Default \code{"clip"}.
#' @param n_neighbors integer(1). Number of PSM nearest neighbors.
#'   Default 1L. Only used when \code{estimator = "psm"}.
#' @param caliper numeric(1) or NULL. PSM caliper distance. Must be
#'   positive and finite if specified.
#' @param caliper_scale character(1). One of \code{"sd"},
#'   \code{"absolute"}. Default \code{"sd"}.
#' @param with_replacement logical(1). Whether PSM uses replacement.
#'   Default TRUE.
#' @param match_order character(1). PSM match order. One of
#'   \code{"data"}, \code{"random"}, \code{"largest"},
#'   \code{"smallest"}. Default \code{"data"}.
#' @param se_method character(1) or NULL. Standard error method.
#'   For IPW/IPWRA: \code{"analytical"} or \code{"bootstrap"}.
#'   For PSM: \code{"abadie_imbens"} or \code{"bootstrap"}.
#'   NULL uses estimator default. Ignored for RA.
#' @param boot_reps integer(1). Number of bootstrap replications.
#'   Default 200L.
#' @param seed integer(1) or NULL. Random seed for bootstrap/PSM.
#' @param return_diagnostics logical(1). Whether to return propensity
#'   score diagnostics. Default FALSE.
#'
#' @return A list with estimator-specific results, plus:
#'   \describe{
#'     \item{estimator}{character(1). The estimator used.}
#'     \item{inference_dist}{character(1). \code{"t"} for RA
#'       (Lee & Wooldridge 2026, Equation 2.10 exact t-inference),
#'       \code{"normal"} for IPW/IPWRA/PSM (asymptotic normality).}
#'   }
#'
#' @references
#' Lee, S. and Wooldridge, J.M. (2025). "Difference-in-Differences
#'   with a Continuous Treatment."
#'
#' Lee, S. and Wooldridge, J.M. (2026). "Simple Difference-in-Differences
#'   Estimation in Fixed-T Panels."
#'
#' Lunceford, J.K. and Davidian, M. (2004). "Stratification and Weighting
#'   via the Propensity Score in Estimation of Causal Treatment Effects."
#'   \emph{Statistics in Medicine}, 23(19), 2937-2960.
#'
#' Cattaneo, M.D. (2010). "Efficient Semiparametric Estimation of
#'   Multi-Valued Treatment Effects under Ignorability."
#'   \emph{Journal of Econometrics}, 155(2), 138-154.
#'
#' Abadie, A. and Imbens, G.W. (2006). "Large Sample Properties of
#'   Matching Estimators for Average Treatment Effects."
#'   \emph{Econometrica}, 74(1), 235-267.
#'
#' @seealso \code{\link{estimate_ra_common}}, \code{\link{estimate_ipw}},
#'   \code{\link{estimate_ipwra}}, \code{\link{estimate_psm}},
#'   \code{\link{lwdid}}
#'
#' @keywords internal
dispatch_estimator <- function(data, y, d, controls = NULL,
    ps_controls = NULL, estimator = c("ra", "ipw", "ipwra", "psm"),
    vce = NULL, cluster_var = NULL, alpha = 0.05,
    trim_threshold = 0.01, trim_method = c("clip", "drop"),
    n_neighbors = 1L, caliper = NULL,
    caliper_scale = c("sd", "absolute"), with_replacement = TRUE,
    match_order = c("data", "random", "largest", "smallest"),
    se_method = NULL, boot_reps = 200L, seed = NULL,
    return_diagnostics = FALSE) {
  estimator <- match.arg(estimator)
  trim_method <- match.arg(trim_method)
  caliper_scale <- match.arg(caliper_scale)
  match_order <- match.arg(match_order)
  if (!is.null(controls) && length(controls) == 0L) controls <- NULL
  if (!is.null(ps_controls) && length(ps_controls) == 0L) ps_controls <- NULL
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(sprintf("alpha must be in (0, 1), got: %s", deparse(alpha)), class = "lwdid_invalid_parameter")
  }
  if (estimator == "ipwra" && is.null(controls)) {
    stop_lwdid("estimator='ipwra' requires 'controls' for the outcome model, but controls is NULL.", class = "lwdid_invalid_parameter")
  }
  if (estimator %in% c("ipw", "psm") && is.null(controls) && is.null(ps_controls)) {
    stop_lwdid(sprintf("estimator='%s' requires 'controls' or 'ps_controls', but both are NULL.", estimator), class = "lwdid_invalid_parameter")
  }
  if (estimator == "ra") {
    if (!is.null(ps_controls)) message("estimator='ra': ps_controls ignored.")
    if (!identical(trim_threshold, 0.01)) message("estimator='ra': trim_threshold ignored.")
    if (!is.null(se_method)) message("estimator='ra': se_method ignored.")
    if (isTRUE(return_diagnostics)) message("estimator='ra': return_diagnostics ignored.")
  }
  if (estimator %in% c("ipw", "ipwra", "psm")) {
    if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L || is.na(trim_threshold) || trim_threshold <= 0 || trim_threshold >= 0.5) {
      stop_lwdid(sprintf("trim_threshold must be in (0, 0.5), got: %s", deparse(trim_threshold)), class = "lwdid_invalid_parameter")
    }
  }
  if (estimator == "psm") {
    if (!is.numeric(n_neighbors) || length(n_neighbors) != 1L || is.na(n_neighbors) || n_neighbors < 1 || n_neighbors != floor(n_neighbors)) {
      stop_lwdid(sprintf("n_neighbors must be a positive integer (>= 1), got: %s", deparse(n_neighbors)), class = "lwdid_invalid_parameter")
    }
    n_neighbors <- as.integer(n_neighbors)
    if (!is.null(caliper) && (!is.numeric(caliper) || length(caliper) != 1L || !is.finite(caliper) || caliper <= 0)) {
      stop_lwdid(sprintf("caliper must be a positive finite number or NULL, got: %s", deparse(caliper)), class = "lwdid_invalid_parameter")
    }
    if (!is.logical(with_replacement) || length(with_replacement) != 1L || is.na(with_replacement)) {
      stop_lwdid("with_replacement must be TRUE or FALSE.", class = "lwdid_invalid_parameter")
    }
  }
  if (!is.null(se_method) && estimator != "ra") {
    valid_se <- switch(estimator, ipw = c("analytical", "bootstrap"), ipwra = c("analytical", "bootstrap"), psm = c("abadie_imbens", "bootstrap"))
    if (!se_method %in% valid_se) {
      stop_lwdid(sprintf("Invalid se_method '%s' for estimator='%s'. Valid: %s", se_method, estimator, paste(valid_se, collapse = ", ")), class = "lwdid_invalid_parameter")
    }
  }
  ps_vars <- if (!is.null(ps_controls)) ps_controls else controls
  result <- switch(estimator,
    ra = {
      x_mat <- if (!is.null(controls)) as.matrix(data[, controls, drop = FALSE]) else NULL
      estimate_ra_common(y_trans = as.numeric(data[[y]]), d = as.integer(data[[d]]), x = x_mat, vce = vce, cluster = cluster_var, alpha = alpha)
    },
    ipw = {
      ipw_vce <- if (!is.null(se_method) && se_method == "bootstrap") "bootstrap" else NULL
      estimate_ipw(data = data, y = y, d = d, propensity_controls = ps_vars, trim_threshold = trim_threshold, trim_method = trim_method, vce = ipw_vce, alpha = alpha, boot_reps = boot_reps, seed = seed)
    },
    ipwra = {
      ipwra_vce <- if (!is.null(se_method) && se_method == "bootstrap") "bootstrap" else NULL
      estimate_ipwra(data = data, y = y, d = d, controls = controls, propensity_controls = ps_vars, trim_threshold = trim_threshold, trim_method = trim_method, vce = ipwra_vce, alpha = alpha, boot_reps = boot_reps, seed = seed)
    },
    psm = {
      psm_se <- if (is.null(se_method)) "abadie_imbens" else se_method
      estimate_psm(data = data, y = y, d = d, propensity_controls = ps_vars, n_neighbors = n_neighbors, with_replacement = with_replacement, caliper = caliper, caliper_scale = caliper_scale, trim_threshold = trim_threshold, se_method = psm_se, n_bootstrap = boot_reps, seed = seed, alpha = alpha, match_order = match_order)
    }
  )
  result$estimator <- estimator
  result$inference_dist <- if (estimator == "ra") "t" else "normal"
  result
}
