#' @title lwdid Internal Utility Functions and Constants
#' @description Cross-module shared utility functions and numerical constants.
#' @name lwdid-utils
#' @keywords internal
NULL

# ============================================================================
# Numerical Constants
# ============================================================================

#' @keywords internal
LWDID_NT_ZERO_TOLERANCE <- 1e-10

#' @keywords internal
LWDID_COHORT_FLOAT_TOLERANCE <- 1e-9

#' @keywords internal
LWDID_WEIGHT_SUM_TOLERANCE <- 1e-6

#' @keywords internal
LWDID_NUMERICAL_TOLERANCE <- 1e-12

#' @keywords internal
LWDID_VARIANCE_THRESHOLD <- 1e-10

#' @keywords internal
LWDID_PS_TRIM_DEFAULT <- 0.01

#' @keywords internal
LWDID_STATA_ATT_TOLERANCE <- 1e-6

#' @keywords internal
LWDID_PYTHON_ATT_TOLERANCE <- 1e-5

# ============================================================================
# Valid Value Sets
# ============================================================================

#' @keywords internal
LWDID_VALID_ESTIMATORS <- c("ra", "ipw", "ipwra", "psm")

#' @keywords internal
LWDID_VALID_TRANSFORMATIONS <- c("demean", "detrend", "demeanq", "detrendq")

#' @keywords internal
LWDID_VALID_BASE_PERIODS <- c("universal", "varying")

#' @keywords internal
LWDID_VALID_BOOT_TYPES <- c("weighted", "multiplier", "bayesian")

#' @keywords internal
LWDID_VALID_AGGREGATE_LEVELS <- c("none", "cohort", "overall", "event_time")

# ============================================================================
# NULL Coalescing Operator
# ============================================================================

#' NULL coalescing operator
#'
#' If left operand is NULL, return right operand; otherwise return left.
#' Equivalent to Python's `a if a is not None else b`.
#' R 4.4.0+ includes this in base, but this package targets R >= 4.0.0,
#' so we define it explicitly. Package namespace definition does not
#' conflict with base.
#'
#' @param a Any R object.
#' @param b Any R object, fallback when a is NULL.
#' @return a if non-NULL, otherwise b.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================================
# NT Unit Identification Functions
# ============================================================================

#' Check if gvar values indicate never-treated units
#'
#' Vectorized identification of never-treated (NT) units from cohort
#' variable values. Supports three NT markers: NA/NaN, 0 (including
#' near-zero values within tolerance), and +Inf. Negative infinity
#' (-Inf) is treated as an invalid value and triggers an error.
#'
#' @param g Numeric scalar or vector of gvar (cohort variable) values.
#' @return Logical vector of same length as \code{g}. TRUE indicates
#'   a never-treated unit.
#' @export
#' @examples
#' is_never_treated(c(NA, 0, Inf, 2005, 2010))
#' # [1]  TRUE  TRUE  TRUE FALSE FALSE
#'
#' is_never_treated(c(1e-11, -1e-11, 1e-9))
#' # [1]  TRUE  TRUE FALSE
is_never_treated <- function(g) {
  if (length(g) == 0L) return(logical(0))
  # Check for -Inf first (error condition)
  if (any(!is.na(g) & is.infinite(g) & g < 0)) {
    stop_lwdid(
      "Negative infinity (-Inf) is not a valid gvar value.",
      class = "lwdid_invalid_staggered_data",
      gvar = "gvar", invalid_values = "-Inf",
      detail = "Use NA, 0, or Inf to indicate never-treated units."
    )
  }
  is_na <- is.na(g)
  is_pos_inf <- !is_na & is.infinite(g) & g > 0
  is_near_zero <- !is_na & abs(g) < LWDID_NT_ZERO_TOLERANCE
  is_na | is_pos_inf | is_near_zero
}

#' Extract sorted unique treatment cohorts from panel data
#'
#' @param data data.frame or data.table, panel data.
#' @param gvar Character, name of the cohort variable column.
#' @param ivar Character, name of the unit identifier column.
#' @param never_treated_values Numeric vector of explicit NT marker values.
#' @return Sorted integer vector of valid cohort values.
#' @keywords internal
get_cohorts <- function(data, gvar, ivar, never_treated_values = c(0, Inf)) {
  unit_gvar_dt <- get_unit_level_gvar(data.table::as.data.table(data), gvar, ivar)
  unit_gvar <- unit_gvar_dt[[gvar]]
  nt_mask <- is_never_treated(unit_gvar)
  g <- unit_gvar[!nt_mask]
  # Also exclude custom never_treated_values
  g <- g[!(g %in% never_treated_values)]
  if (length(g) == 0L) return(integer(0))
  sort(as.integer(unique(g)))
}

#' Check if panel data contains never-treated units
#'
#' @param data data.frame or data.table, panel data.
#' @param gvar Character, name of the cohort variable column.
#' @param ivar Character, name of the unit identifier column.
#' @return Single logical value.
#' @keywords internal
has_never_treated <- function(data, gvar, ivar) {
  unit_gvar_dt <- get_unit_level_gvar(data.table::as.data.table(data), gvar, ivar)
  any(is_never_treated(unit_gvar_dt[[gvar]]))
}

#' Get boolean mask for never-treated units
#'
#' @param unit_gvar Numeric vector of unit-level gvar values.
#' @return Logical vector (TRUE = never-treated).
#' @keywords internal
get_never_treated_mask <- function(unit_gvar) {
  is_never_treated(unit_gvar)
}

#' Get boolean mask for units belonging to a specific cohort
#'
#' Uses floating-point tolerance for comparison.
#'
#' @param unit_gvar Numeric vector of unit-level gvar values.
#' @param g Scalar integer, target cohort value.
#' @return Logical vector (TRUE = belongs to cohort g).
#' @keywords internal
get_cohort_mask <- function(unit_gvar, g) {
  result <- abs(unit_gvar - g) < LWDID_COHORT_FLOAT_TOLERANCE
  result[is.na(result)] <- FALSE
  result
}

# ============================================================================
# Period and Subset Extraction Functions
# ============================================================================

#' Get valid post-treatment periods for a cohort
#'
#' @param cohort Scalar integer, treatment cohort (first treatment period g).
#' @param T_max Scalar integer, maximum period in the panel.
#' @return Integer vector `{g, g+1, ..., T_max}`, or `integer(0)` if
#'   `cohort > T_max`.
#' @keywords internal
get_valid_periods_for_cohort <- function(cohort, T_max) {
  cohort <- as.integer(cohort)
  T_max <- as.integer(T_max)
  if (cohort > T_max) return(integer(0))
  seq.int(cohort, T_max)
}

#' Extract cross-section for a single period
#'
#' @param dt data.table, panel data.
#' @param tvar Character, name of the time variable column.
#' @param period Scalar numeric, target period.
#' @return data.table subset with rows where `tvar == period`.
#' @keywords internal
extract_cross_section <- function(dt, tvar, period) {
  dt[dt[[tvar]] == period, ]
}

# ============================================================================
# Numerical Utility Functions
# ============================================================================

#' Safe mean that returns NA_real_ for empty or all-NA inputs
#'
#' @param x Numeric vector.
#' @param na.rm Logical, whether to remove NAs (default TRUE).
#' @return Single numeric value.
#' @keywords internal
safe_mean <- function(x, na.rm = TRUE) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

#' Center control variables by treated group mean
#'
#' Computes X_centered = X - X_bar_treated (equation 3.3).
#' Centering is applied to ALL units, not just controls.
#'
#' @param X Numeric matrix (N x K), control variables.
#' @param d Numeric vector (length N), treatment indicator (0/1).
#' @return Named list with `$X_centered` (N x K matrix) and
#'   `$X_mean_treated` (length K vector).
#' @keywords internal
center_by_treated <- function(X, d) {
  if (!is.matrix(X)) X <- as.matrix(X)
  treated_mask <- d == 1
  n_treated <- sum(treated_mask)
  if (n_treated == 0L) {
    stop_lwdid(
      "No treated units found for centering.",
      class = c("lwdid_no_treated", "lwdid_insufficient_data"),
      treat_var = "d", n_units = length(d)
    )
  }
  X_mean_treated <- colMeans(X[treated_mask, , drop = FALSE])
  X_centered <- sweep(X, 2, X_mean_treated, "-")
  list(X_centered = X_centered, X_mean_treated = X_mean_treated)
}

#' Numerically stable OLS via QR decomposition
#'
#' @param X Numeric matrix (N x K), design matrix (including intercept).
#' @param y Numeric vector (length N), dependent variable.
#' @param w Numeric vector (length N) or NULL, observation weights.
#' @return Named list with `$coefficients`, `$residuals`, `$fitted.values`,
#'   `$rank`, `$df.residual`, `$qr`, `$is_rank_deficient`, `$weights`.
#' @importFrom stats .lm.fit lm.wfit
#' @keywords internal
stable_ols <- function(X, y, w = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (!is.null(w)) {
    # Validate weights
    if (length(w) != n) {
      stop_lwdid(
        sprintf("Weight vector length (%d) != number of observations (%d).",
                length(w), n),
        class = "lwdid_invalid_parameter",
        param = "w", value = length(w), allowed = n
      )
    }
    if (any(is.na(w))) {
      stop_lwdid(
        "Weights contain NA values.",
        class = "lwdid_invalid_parameter",
        param = "w", value = "contains NA", allowed = "no NA"
      )
    }
    if (any(w <= 0)) {
      stop_lwdid(
        "Weights must be strictly positive.",
        class = "lwdid_invalid_parameter",
        param = "w", value = "contains non-positive", allowed = "all > 0"
      )
    }
    fit <- lm.wfit(X, y, w)
  } else {
    fit <- .lm.fit(X, y)
  }

  coefficients <- fit$coefficients
  if (is.null(names(coefficients)) && !is.null(colnames(X))) {
    names(coefficients) <- colnames(X)
  }

  rank_val <- if (!is.null(fit$rank)) fit$rank else p

  # CRITICAL: Use fit$residuals directly — never compute via X %*% coefficients

  # because coefficients contain NA when rank-deficient.
  # Both .lm.fit() and lm.wfit() return correct unweighted residuals.
  residuals_val <- fit$residuals
  fitted_vals <- y - residuals_val

  list(
    coefficients    = coefficients,
    residuals       = residuals_val,
    fitted.values   = fitted_vals,
    rank            = rank_val,
    df.residual     = as.integer(n - rank_val),
    qr              = fit$qr,
    is_rank_deficient = rank_val < p,
    weights         = w
  )
}

#' Check matrix rank via QR decomposition
#'
#' @param X Numeric matrix, design matrix.
#' @param tol Numeric scalar or NULL, rank determination tolerance.
#'   NULL uses R's default (1e-7).
#' @return Named list with `$rank`, `$ncol`, `$is_full_rank`,
#'   `$rank_deficiency`, `$problematic_columns`.
#' @keywords internal
check_rank <- function(X, tol = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  qr_args <- list(x = X)
  if (!is.null(tol)) qr_args$tol <- tol
  qr_obj <- do.call(qr, qr_args)
  r <- qr_obj$rank
  p <- ncol(X)
  problematic <- if (r < p) qr_obj$pivot[(r + 1):p] else integer(0)
  list(
    rank = r,
    ncol = p,
    is_full_rank = (r == p),
    rank_deficiency = p - r,
    problematic_columns = problematic
  )
}


# ============================================================================
# Memory Management Utilities (Epic 13, Story 13.4)
# ============================================================================

#' Estimate memory usage of a data object
#'
#' Reports the approximate memory footprint of \code{x} in bytes,
#' kilobytes, and megabytes. Useful for deciding whether to enable
#' parallel execution (where each worker gets a copy).
#'
#' @param x Any R object.
#' @return Named list with \code{bytes}, \code{kb}, \code{mb}.
#' @keywords internal
check_memory <- function(x) {
  sz <- as.numeric(utils::object.size(x))
  list(
    bytes = sz,
    kb    = sz / 1024,
    mb    = sz / 1024^2
  )
}

#' Determine optimal chunk size for processing large data
#'
#' Computes chunk size based on available memory and data size.
#' Used for chunked processing when data exceeds memory budget.
#'
#' @param n_total Integer. Total number of items to process.
#' @param item_size_bytes Numeric. Approximate size of each item in bytes.
#' @param memory_budget_mb Numeric. Memory budget in megabytes (default 500).
#' @param min_chunk Integer. Minimum chunk size (default 10).
#' @return Integer chunk size.
#' @keywords internal
auto_chunk <- function(n_total, item_size_bytes,
                       memory_budget_mb = 500,
                       min_chunk = 10L) {
  budget_bytes <- memory_budget_mb * 1024^2
  chunk <- max(min_chunk, as.integer(budget_bytes / max(item_size_bytes, 1)))
  min(chunk, n_total)
}
