# ============================================================================
# lwdid: Variance-Covariance Estimation (VCE) Inference Module
# ============================================================================
# This file implements the VCE dispatcher and robust variance estimation
# for the lwdid package:
#   - compute_vce():         Main VCE dispatcher (homoskedastic / HC / cluster)
#   - compute_hc_vce():      HC0-HC4 heteroskedasticity-robust VCE
#   - compute_cluster_vce(): Cluster-robust VCE (Liang-Zeger / sandwich::vcovCL)
#
# Reference: Lee & Wooldridge (2025, 2026) [lw2026]
#   - Section 2.4: Variance estimation under homoskedasticity
#   - Section 3.3: Robust and cluster-robust inference
#
# Dependencies: sandwich (>= 3.0-0) for vcovHC / vcovCL
# ============================================================================

#' Compute variance-covariance matrix for lwdid estimates
#'
#' Main dispatcher that routes to the appropriate VCE estimator based on
#' the \code{vce} argument. Supports homoskedastic (OLS), HC0--HC4
#' heteroskedasticity-robust, and cluster-robust variance estimation.
#'
#' @param fit A fitted model object (e.g., from \code{lm}) or an lwdid
#'   internal result list containing \code{vcov} and \code{df.residual}
#'   (or \code{df}) fields.
#' @param vce Character string or NULL specifying the VCE type.
#'   \describe{
#'     \item{\code{NULL}}{Homoskedastic OLS variance (default).
#'       Uses \code{vcov(fit)} and \code{fit$df.residual}.}
#'     \item{\code{"hc0"}, \code{"hc1"}, \code{"hc2"}, \code{"hc3"},
#'       \code{"hc4"}}{Heteroskedasticity-consistent (HC) robust
#'       variance. Delegated to \code{compute_hc_vce()}.}
#'     \item{\code{"robust"}}{Alias for \code{"hc1"}.}
#'     \item{\code{"cluster"}}{Cluster-robust variance. Requires
#'       \code{cluster} argument. Delegated to
#'       \code{compute_cluster_vce()}.}
#'   }
#' @param cluster A vector of cluster identifiers (same length as
#'   observations), or NULL. Required when \code{vce = "cluster"}.
#' @param verbose Character string controlling verbosity. Currently
#'   unused; reserved for future diagnostic output. Default
#'   \code{"default"}.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{vcov}{Numeric matrix. The estimated variance-covariance
#'     matrix of the coefficient estimates.}
#'   \item{df}{Integer. Degrees of freedom for inference.
#'     For homoskedastic and HC: \code{fit$df.residual}.
#'     For cluster: number of clusters minus one (from
#'     \code{compute_cluster_vce()}).}
#'   \item{vce_type}{Character string. The canonical VCE type label
#'     (e.g., \code{"homoskedastic"}, \code{"HC1"}, \code{"cluster"}).}
#'   \item{n_clusters}{Integer or NULL. Number of clusters (cluster
#'     path only); NULL for homoskedastic and HC paths.}
#' }
#'
#' @section VCE Type Resolution:
#' The \code{vce} argument is case-insensitive. After
#' \code{tolower()}, the following mappings apply:
#' \itemize{
#'   \item \code{NULL} \eqn{\to} homoskedastic
#'   \item \code{"robust"} \eqn{\to} \code{"hc1"}
#'   \item \code{"hc0"} .. \code{"hc4"} \eqn{\to} HC robust
#'   \item \code{"cluster"} \eqn{\to} cluster-robust
#' }
#'
#' @section Error Conditions:
#' \describe{
#'   \item{\code{lwdid_invalid_vce}}{Thrown when \code{vce} is not a
#'     recognized type.}
#'   \item{\code{lwdid_invalid_parameter}}{Thrown when
#'     \code{vce = "cluster"} but \code{cluster} is NULL.}
#' }
#'
#' @keywords internal
compute_vce <- function(fit, vce = NULL, cluster = NULL,
                        verbose = "default") {

  # ------------------------------------------------------------------
  # Valid VCE types (after tolower)
  # ------------------------------------------------------------------
  hc_types <- c("hc0", "hc1", "hc2", "hc3", "hc4")
  all_valid <- c(hc_types, "robust", "cluster")

  # ------------------------------------------------------------------
  # NULL path: homoskedastic OLS variance
  # ------------------------------------------------------------------
  if (is.null(vce)) {
    vcov_mat <- vcov(fit)
    df_val <- fit$df.residual
    result <- list(
      vcov       = vcov_mat,
      df         = df_val,
      vce_type   = "homoskedastic",
      n_clusters = NULL
    )
  } else {

    # ------------------------------------------------------------------
    # Normalize vce to lowercase
    # ------------------------------------------------------------------
    vce <- tolower(vce)

    # ------------------------------------------------------------------
    # Validate vce type
    # ------------------------------------------------------------------
    if (!vce %in% all_valid) {
      stop_lwdid(
        sprintf(
          "Invalid VCE type '%s'. Must be one of: %s",
          vce, paste(c("NULL", all_valid), collapse = ", ")
        ),
        class = c("lwdid_invalid_vce", "lwdid_invalid_parameter"),
        vce_type = vce,
        allowed = c("NULL", all_valid)
      )
    }

    # ------------------------------------------------------------------
    # "robust" is an alias for "hc1"
    # ------------------------------------------------------------------
    if (vce == "robust") {
      vce <- "hc1"
    }

    # ------------------------------------------------------------------
    # HC path: delegate to compute_hc_vce()
    # ------------------------------------------------------------------
    if (vce %in% hc_types) {
      hc_result <- compute_hc_vce(fit, type = vce)
      result <- list(
        vcov       = hc_result$vcov,
        df         = fit$df.residual,
        vce_type   = toupper(vce),
        n_clusters = NULL
      )

    # ------------------------------------------------------------------
    # Cluster path: validate cluster, delegate to compute_cluster_vce()
    # ------------------------------------------------------------------
    } else if (vce == "cluster") {
      if (is.null(cluster)) {
        stop_lwdid(
          paste0(
            "VCE type 'cluster' requires a non-NULL 'cluster' argument ",
            "specifying cluster identifiers."
          ),
          class = "lwdid_invalid_parameter",
          param = "cluster",
          value = "NULL",
          allowed = "non-NULL vector of cluster identifiers"
        )
      }

      cl_result <- compute_cluster_vce(fit, cluster, type = "HC1")
      result <- list(
        vcov       = cl_result$vcov,
        df         = cl_result$df,
        vce_type   = "cluster",
        n_clusters = cl_result$n_clusters
      )
    }
  }

  # ------------------------------------------------------------------
  # Near-zero SE check for treatment effect "D"
  # Runs after all VCE paths have computed their result, before return.
  # If the model includes a "D" coefficient, check whether its SE is
  # suspiciously small (< 1e-10), which may indicate perfect fit or
  # numerical issues.
  # ------------------------------------------------------------------
  d_idx <- which(names(coef(fit)) == "D")
  if (length(d_idx) == 1L) {
    se_d <- sqrt(diag(result$vcov))[d_idx]
    if (is.finite(se_d) && se_d < 1e-10) {
      warn_lwdid(
        "Near-zero standard error for treatment effect D (SE < 1e-10). This may indicate perfect fit or numerical issues.",
        class = "lwdid_numerical",
        detail = "Near-zero standard error for treatment effect D",
        source_function = "compute_vce"
      )
    }
  }

  return(result)
}


# ============================================================================
# HC Robust VCE (Task E3-01.3)
# ============================================================================

#' Compute HC heteroskedasticity-robust variance-covariance matrix
#'
#' @title HC Robust Variance-Covariance Estimation
#'
#' @description
#' Computes heteroskedasticity-consistent (HC) robust variance-covariance
#' matrices for OLS coefficient estimates using the sandwich package.
#' Supports HC0 through HC4 estimators:
#' \itemize{
#'   \item \strong{HC0}: White (1980) heteroskedasticity-consistent
#'     estimator. No degrees-of-freedom correction.
#'   \item \strong{HC1}: White (1980) with small-sample correction
#'     \eqn{n / (n - k)}. Equivalent to Stata's \code{robust} option.
#'   \item \strong{HC2}: MacKinnon & White (1985). Divides squared
#'     residuals by \eqn{(1 - h_{ii})}, where \eqn{h_{ii}} are hat
#'     matrix diagonal elements.
#'   \item \strong{HC3}: Davidson & MacKinnon (1993). Divides squared
#'     residuals by \eqn{(1 - h_{ii})^2}. Jackknife-like correction
#'     with better finite-sample properties.
#'   \item \strong{HC4}: Cribari-Neto (2004). Adjusts leverage
#'     exponent adaptively based on \eqn{h_{ii}} relative to
#'     \eqn{k/n}, providing further robustness to influential
#'     observations.
#' }
#'
#' @param fit A fitted model object (typically from \code{lm}).
#' @param type Character string specifying the HC type. One of
#'   \code{"hc0"}, \code{"hc1"}, \code{"hc2"}, \code{"hc3"}, or
#'   \code{"hc4"} (case-insensitive). Default \code{"HC3"}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{vcov}{Numeric matrix. The HC robust variance-covariance
#'     matrix of the coefficient estimates.}
#'   \item{se}{Numeric vector. HC robust standard errors
#'     (\code{sqrt(diag(vcov))}).}
#' }
#'
#' @references
#' White, H. (1980). A heteroskedasticity-consistent covariance
#'   matrix estimator and a direct test for heteroskedasticity.
#'   \emph{Econometrica}, 48(4), 817--838.
#'
#' MacKinnon, J. G. and White, H. (1985). Some
#'   heteroskedasticity-consistent covariance matrix estimators with
#'   improved finite sample properties. \emph{Journal of
#'   Econometrics}, 29(3), 305--325.
#'
#' Davidson, R. and MacKinnon, J. G. (1993). \emph{Estimation and
#'   Inference in Econometrics}. Oxford University Press.
#'
#' Cribari-Neto, F. (2004). Asymptotic inference under
#'   heteroskedasticity of unknown form. \emph{Computational
#'   Statistics & Data Analysis}, 45(2), 215--233.
#'
#' @importFrom sandwich vcovHC
#' @keywords internal
compute_hc_vce <- function(fit, type = "HC3") {

  # ------------------------------------------------------------------
  # Input validation: normalize and check type
  # ------------------------------------------------------------------
  type <- tolower(type)
  valid_types <- c("hc0", "hc1", "hc2", "hc3", "hc4")

  if (!type %in% valid_types) {
    stop_lwdid(
      sprintf(
        "Invalid HC type '%s'. Must be one of: hc0, hc1, hc2, hc3, hc4",
        type
      ),
      class = "lwdid_invalid_parameter",
      param = "type",
      value = type,
      allowed = valid_types
    )
  }

  # ------------------------------------------------------------------
  # Small sample warning (HC1, HC3 only)
  # Check treated/control counts via "D" column in model frame.
  # ------------------------------------------------------------------
  if (type %in% c("hc1", "hc3")) {
    mf <- model.frame(fit)
    if ("D" %in% names(mf)) {
      d_vec <- mf[["D"]]
      n_treated <- sum(d_vec == 1L)
      n_control <- sum(d_vec == 0L)
      if (n_treated < 2L || n_control < 2L) {
        warn_lwdid(
          sprintf(
            paste0(
              "Small sample for %s: N_treated=%d, ",
              "N_control=%d. HC standard errors ",
              "may be unreliable."
            ),
            toupper(type), n_treated, n_control
          ),
          class = "lwdid_small_sample",
          n = nobs(fit),
          n_treated = n_treated,
          n_control = n_control
        )
      }
    }
  }

  # ------------------------------------------------------------------
  # High leverage warning (HC2, HC3, HC4 only)
  # ------------------------------------------------------------------
  if (type %in% c("hc2", "hc3", "hc4")) {
    h <- hatvalues(fit)
    max_h <- max(h)
    if (max_h >= 1 - .Machine$double.eps) {
      warn_lwdid(
        sprintf(
          paste0(
            "Hat matrix diagonal contains h_ii = 1 (or near 1). ",
            "%s divides by (1-h_ii)^%s which is zero, producing ",
            "infinite/NA standard errors. This typically occurs ",
            "with a single treated unit (N_1=1). Use vce=NULL ",
            "(homoskedastic) or vce='hc1' instead."
          ),
          toupper(type),
          if (type == "hc3") "2" else if (type == "hc2") "1" else "delta"
        ),
        class = "lwdid_numerical",
        detail = sprintf(
          "Hat value h_ii = %.6f (effectively 1)", max_h
        ),
        source_function = "compute_hc_vce"
      )
    } else if (max_h > 0.99) {
      warn_lwdid(
        sprintf(
          paste0(
            "High leverage point detected ",
            "(max h_ii = %.6f > 0.99). ",
            "%s standard errors may be unreliable."
          ),
          max_h, toupper(type)
        ),
        class = "lwdid_numerical",
        detail = sprintf(
          "High leverage point (max h_ii = %.6f)", max_h
        ),
        source_function = "compute_hc_vce"
      )
    }
  }

  # ------------------------------------------------------------------
  # Compute VCE via sandwich::vcovHC (expects uppercase type)
  # ------------------------------------------------------------------
  vcov_mat <- sandwich::vcovHC(fit, type = toupper(type))

  # ------------------------------------------------------------------
  # Return vcov matrix and standard errors
  # ------------------------------------------------------------------
  list(
    vcov = vcov_mat,
    se   = sqrt(diag(vcov_mat))
  )
}


# ============================================================================
# Cluster-Robust VCE (Task E3-01.4)
# ============================================================================

#' Compute cluster-robust variance-covariance matrix
#'
#' @title Cluster-Robust Variance-Covariance Estimation
#'
#' @description
#' Computes cluster-robust (Liang--Zeger) variance-covariance matrices
#' for OLS coefficient estimates using \code{sandwich::vcovCL}.
#'
#' Cluster-robust standard errors account for arbitrary within-cluster
#' correlation of the error terms, which is essential when observations
#' are grouped (e.g., by individual, firm, or region) and independence
#' holds only across clusters, not within them.
#'
#' The estimator is based on the seminal work of Liang & Zeger (1986)
#' and Arellano (1987), with finite-sample corrections discussed in
#' Cameron, Gelbach & Miller (2008) and Cameron & Miller (2015). The
#' implementation delegates to \code{sandwich::vcovCL}, which supports
#' HC0--HC4 small-sample correction types applied at the cluster level.
#'
#' With the default \code{type = "HC1"}, the small-sample correction
#' factor applied is:
#' \deqn{\frac{G}{G - 1} \cdot \frac{N - 1}{N - p}}{G/(G-1) * (N-1)/(N-p)}
#' where \eqn{G} is the number of clusters, \eqn{N} is the total
#' number of observations, and \eqn{p} is the number of estimated
#' parameters. This matches Stata's default \code{vce(cluster)}
#' behaviour; see Zeileis, K\"oll & Graham (2020, \emph{Journal of
#' Statistical Software}, 95(1)) for a detailed comparison of R and
#' Stata cluster-robust implementations.
#'
#' The \code{cluster} variable may be numeric, character, or factor.
#' \code{sandwich::vcovCL} internally converts it to a factor for
#' grouping, so any type that can be meaningfully coerced to factor
#' is accepted.
#'
#' @param fit A fitted model object (typically from \code{lm}).
#' @param cluster A vector of cluster identifiers with the same length
#'   as the number of observations in \code{fit}. Must not contain NA
#'   values. Can be numeric, character, or factor.
#' @param type Character string specifying the small-sample correction
#'   type passed to \code{sandwich::vcovCL}. Default \code{"HC1"},
#'   which applies the \eqn{G / (G - 1)} degrees-of-freedom correction
#'   (where \eqn{G} is the number of clusters), analogous to Stata's
#'   default cluster-robust standard errors.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{vcov}{Numeric matrix. The cluster-robust variance-covariance
#'     matrix of the coefficient estimates.}
#'   \item{df}{Integer. Degrees of freedom for inference, equal to
#'     \eqn{G - 1} (number of clusters minus one).}
#'   \item{n_clusters}{Integer. The number of unique clusters \eqn{G}.}
#' }
#'
#' @references
#' Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis
#'   using generalized linear models. \emph{Biometrika}, 73(1), 13--22.
#'
#' Arellano, M. (1987). Computing robust standard errors for
#'   within-groups estimators. \emph{Oxford Bulletin of Economics and
#'   Statistics}, 49(4), 431--434.
#'
#' Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2008).
#'   Bootstrap-based improvements for inference with clustered errors.
#'   \emph{Review of Economics and Statistics}, 90(3), 414--427.
#'
#' Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
#'   cluster-robust inference. \emph{Journal of Human Resources},
#'   50(2), 317--372.
#'
#' Zeileis, A., K\"oll, S. and Graham, N. (2020). Various versatile
#'   variances: An object-oriented implementation of clustered
#'   covariances in R. \emph{Journal of Statistical Software},
#'   95(1), 1--36.
#'
#' @importFrom sandwich vcovCL
#' @keywords internal
compute_cluster_vce <- function(fit, cluster, type = "HC1") {

  # ------------------------------------------------------------------
  # 1. NA validation (before G calculation)
  # ------------------------------------------------------------------
  if (anyNA(cluster)) {
    n_na <- sum(is.na(cluster))
    stop_lwdid(
      sprintf("Cluster variable contains %d NA value%s. All observations must have valid cluster identifiers.",
              n_na, if (n_na == 1L) "" else "s"),
      class = "lwdid_invalid_parameter",
      param = "cluster", value = sprintf("contains %d NA", n_na),
      allowed = "non-NA cluster identifiers"
    )
  }

  # ------------------------------------------------------------------
  # 2. Length validation
  # ------------------------------------------------------------------
  if (length(cluster) != nobs(fit)) {
    stop_lwdid(
      sprintf("Cluster variable length (%d) does not match number of observations (%d).",
              length(cluster), nobs(fit)),
      class = "lwdid_invalid_parameter",
      param = "cluster", value = length(cluster),
      allowed = sprintf("length %d", nobs(fit))
    )
  }

  # ------------------------------------------------------------------
  # 3. Cluster count: require G >= 2
  # Cameron, Gelbach & Miller (2008) and Cameron & Miller (2015)
  # emphasize that cluster-robust inference requires a sufficient
  # number of clusters G for asymptotic validity.
  # ------------------------------------------------------------------
  G <- length(unique(cluster))

  if (G < 2L) {
    stop_lwdid(
      sprintf("Insufficient number of clusters (G=%d). At least 2 clusters required for cluster-robust VCE.", G),
      class = "lwdid_insufficient_data",
      n = nobs(fit), n_treated = NA_integer_, n_control = NA_integer_
    )
  }

  # ------------------------------------------------------------------
  # 4. Graded cluster warnings
  # Cameron, Gelbach & Miller (2008, Section II) show that cluster-
  # robust SEs can severely over-reject with few clusters. Cameron &
  # Miller (2015, Section V) recommend G >= 20 as a minimum and
  # suggest Wild Cluster Bootstrap when G is small.
  # ------------------------------------------------------------------
  if (G < 10L) {
    # Strong warning: very few clusters
    warn_lwdid(
      sprintf("Very few clusters (G=%d). Cluster-robust standard errors are highly unreliable. Consider using Wild Cluster Bootstrap.", G),
      class = "lwdid_small_sample",
      n = nobs(fit), n_treated = NA_integer_, n_control = NA_integer_
    )
  } else if (G < 20L) {
    # Informational warning: few clusters
    warn_lwdid(
      sprintf("Few clusters (G=%d). Consider using Wild Cluster Bootstrap for more reliable inference.", G),
      class = "lwdid_small_sample",
      n = nobs(fit), n_treated = NA_integer_, n_control = NA_integer_
    )
  }
  # G >= 20: no warning

  # ------------------------------------------------------------------
  # 5. CV imbalance check
  # Cameron & Miller (2015, Section VII) discuss how unbalanced
  # cluster sizes can distort cluster-robust inference.
  # ------------------------------------------------------------------
  sizes <- as.integer(table(cluster))
  cv <- sd(sizes) / mean(sizes)

  if (cv > 1.0) {
    warn_lwdid(
      sprintf("Highly unbalanced cluster sizes (CV=%.2f > 1.0). Cluster-robust inference may be unreliable.", cv),
      class = "lwdid_small_sample",
      n = nobs(fit), n_treated = NA_integer_, n_control = NA_integer_
    )
  }

  # ------------------------------------------------------------------
  # 6. Compute cluster-robust VCE via sandwich::vcovCL
  # Liang & Zeger (1986) / Arellano (1987) sandwich estimator with
  # small-sample correction per Cameron, Gelbach & Miller (2008).
  # sandwich::vcovCL with type="HC1" applies the G/(G-1)*(N-1)/(N-p)
  # correction, matching Stata's vce(cluster) default; see Zeileis,
  # Köll & Graham (2020, JSS 95(1)) for implementation details.
  # ------------------------------------------------------------------
  vcov_mat <- sandwich::vcovCL(fit, cluster = cluster, type = type)

  # ------------------------------------------------------------------
  # 7. Return
  # ------------------------------------------------------------------
  list(
    vcov       = vcov_mat,
    df         = G - 1L,
    n_clusters = G
  )
}


# ============================================================================
# Unified t-Inference (Task E3-04.1)
# ============================================================================

#' @title Unified t-Inference for ATT
#' @description Compute complete t-distribution inference results from an
#'   ATT point estimate, variance-covariance matrix, and degrees of freedom.
#' @details Uses the t-distribution, rather than the normal distribution,
#'   for p-values and confidence intervals, following lw2026 Equation 2.10.
#'   This is a pure computation helper with no dependency on the sandwich
#'   package.
#'
#' @param att Numeric scalar. The ATT point estimate.
#' @param vcov_mat Numeric matrix of coefficient variances and covariances.
#' @param df Integer degrees of freedom for the t-distribution.
#'   Use residual degrees of freedom for homoskedastic or HC inference,
#'   and `G - 1` for cluster inference.
#' @param alpha Numeric scalar in (0, 1). Significance level for
#'   confidence intervals. Default 0.05 for a 95 percent interval.
#' @param coef_index Integer index of the ATT coefficient on the VCE
#'   diagonal. Default 2L because the ATT term follows the intercept in
#'   the current regression layout.
#'
#' @return Named list with elements `att`, `se`, `t_stat`, `df`, `pvalue`,
#'   `ci_lower`, and `ci_upper`.
#'
#' @section Validation Order:
#' Parameters are validated in this order: alpha, att, coef_index, df,
#' var_att, SE. This ensures the most informative error is thrown first.
#'
#' @section Distribution Choice:
#' Uses \code{pt()}/\code{qt()} (t-distribution), NOT
#' \code{pnorm()}/\code{qnorm()} (normal distribution). Under the
#' normality assumption (lw2026 Eq. 2.9/2.19), the t-statistic follows
#' an exact t-distribution. For small samples, the t-distribution has
#' heavier tails than the normal, yielding wider confidence intervals
#' and larger p-values.
#'
#' @references
#' Lee, S. and Wooldridge, J. M. (2025, 2026). A simple approach to
#'   difference-in-differences estimation. Equation 2.10.
#'
#' @keywords internal
compute_inference <- function(att, vcov_mat, df, alpha = 0.05,
                              coef_index = 2L) {

  # --- 1. Alpha validation (checked first) ---
  # alpha must be in (0, 1) open interval:
  # alpha <= 0 or >= 1 causes qt() to return Inf/NaN
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(
      sprintf(
        "alpha must be a numeric value in (0, 1). Got: %s",
        deparse(alpha)
      ),
      class = "lwdid_invalid_parameter",
      param = "alpha",
      value = alpha
    )
  }

  # --- 2. ATT validation ---
  # NA/NaN ATT (e.g., from rank-deficient lm()) must be caught
  # before t_stat = att / se to prevent silent NA propagation
  if (is.na(att) || is.nan(att)) {
    stop_lwdid(
      paste0(
        "ATT estimate is NA/NaN. Cannot perform inference. ",
        "Possible causes: (1) rank-deficient design matrix ",
        "(D coefficient not estimable), (2) data issues causing ",
        "regression failure."
      ),
      class = "lwdid_numerical",
      detail = "ATT is NA/NaN",
      source_function = "compute_inference"
    )
  }

  # --- 3. coef_index bounds check ---
  if (coef_index < 1L || coef_index > nrow(vcov_mat)) {
    stop_lwdid(
      sprintf(
        "coef_index=%d is out of VCE matrix bounds (1 to %d).",
        coef_index, nrow(vcov_mat)
      ),
      class = "lwdid_invalid_parameter",
      param = "coef_index",
      value = coef_index,
      allowed = sprintf("1 to %d", nrow(vcov_mat))
    )
  }

  # --- 4. df validation ---
  # df = NA or df <= 0 means t-distribution is undefined
  if (is.na(df) || df <= 0L) {
    stop_lwdid(
      sprintf(
        paste0(
          "Invalid degrees of freedom (df=%s). ",
          "Cannot perform t-inference. ",
          "Possible causes: number of parameters (p) >= ",
          "number of observations (N), or insufficient ",
          "clusters (G)."
        ),
        deparse(df)
      ),
      class = "lwdid_insufficient_data",
      n = NA_integer_,
      n_treated = NA_integer_,
      n_control = NA_integer_
    )
  }

  # --- 5. Variance extraction and validation ---
  var_att <- vcov_mat[coef_index, coef_index]

  # NA/NaN: numerical issues or perfect collinearity
  if (is.na(var_att) || is.nan(var_att)) {
    stop_lwdid(
      sprintf(
        paste0(
          "VCE matrix diagonal element is NA/NaN at position ",
          "[%d,%d]. Cannot compute standard error. ",
          "Possible causes: (1) perfect collinearity, ",
          "(2) numerical overflow, (3) singular design matrix."
        ),
        coef_index, coef_index
      ),
      class = "lwdid_numerical",
      detail = sprintf(
        "VCE diagonal NA/NaN at [%d,%d]", coef_index, coef_index
      ),
      source_function = "compute_inference"
    )
  }

  # Negative variance: numerical precision issue
  if (var_att < 0) {
    stop_lwdid(
      sprintf(
        paste0(
          "VCE matrix diagonal element is negative ",
          "(Var=%.2e at position [%d,%d]). ",
          "Cannot compute standard error (sqrt of negative). ",
          "This is typically caused by numerical precision issues."
        ),
        var_att, coef_index, coef_index
      ),
      class = "lwdid_numerical",
      detail = sprintf(
        "Negative variance %.2e at [%d,%d]",
        var_att, coef_index, coef_index
      ),
      source_function = "compute_inference"
    )
  }

  # --- 6. SE computation ---
  se <- sqrt(var_att)

  # SE = 0: perfect fit or degenerate VCE
  # Warn but don't stop — R handles Inf/NaN from division gracefully:
  # - att != 0: t_stat = Inf, pvalue = 0, CI = [att, att]
  # - att == 0: t_stat = NaN, pvalue = NaN, CI = [0, 0]
  if (se == 0) {
    warn_lwdid(
      paste0(
        "Standard error is zero (SE=0). Inference results ",
        "(t-statistic, p-value, CI) are meaningless. ",
        "Possible causes: perfect fit, degenerate VCE matrix, ",
        "or data issues."
      ),
      class = "lwdid_numerical",
      detail = "SE is zero",
      source_function = "compute_inference"
    )
  }

  # --- 7. t-inference (lw2026 Eq. 2.10) ---
  # t-statistic
  t_stat <- att / se

  # Two-sided p-value using t-distribution (NOT normal)
  # Equivalent to Python: 2 * stats.t.cdf(-abs(t_stat), df)
  pvalue <- 2 * pt(-abs(t_stat), df)

  # Confidence interval using t-distribution quantiles (NOT normal)
  # Equivalent to Python: att +/- stats.t.ppf(1 - alpha/2, df) * se
  t_crit <- qt(1 - alpha / 2, df)
  ci_lower <- att - t_crit * se
  ci_upper <- att + t_crit * se

  # --- 8. Return ---
  list(
    att      = att,
    se       = se,
    t_stat   = t_stat,
    df       = df,
    pvalue   = pvalue,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
}
