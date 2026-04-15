# ============================================================================
# simulate.R — Simulated panel data generation
#
# Generates simulated panel data for testing and demonstration.
# Supports both Common Timing and Staggered modes.
#
# DGP: y_it = alpha_i + lambda_t + delta_i * t + tau * D_it + X_it * beta + eps_it
#   alpha_i ~ N(0, unit_fe_sd^2)       unit fixed effects
#   lambda_t = time_trend_coef * t      time trend
#   delta_i ~ N(0, trend^2)            unit-specific trend (when trend != NULL)
#   eps_it ~ N(0, error_sd^2)          idiosyncratic error
#   D_it = treatment indicator
#   tau = true_att (scalar or function)
# ============================================================================

#' Generate Simulated Panel Data for DiD Analysis
#'
#' Creates simulated panel data with known treatment effects for testing
#' and demonstration. Supports both common timing and staggered adoption
#' designs.
#'
#' The data generating process (DGP) is:
#' \deqn{y_{it} = \alpha_i + \lambda_t + \delta_i \cdot t +
#'   \tau \cdot D_{it} + X_{it} \beta + \varepsilon_{it}}
#'
#' @param n_units Integer, total number of units (default 20).
#' @param n_periods Integer, total number of time periods (default 10).
#'   Years are labeled 2001 through \code{2000 + n_periods}.
#' @param n_treated Integer or NULL. Number of treated units for common
#'   timing mode. Must be > 0 and < n_units.
#' @param treatment_period Integer or NULL. Period index (1-based) when
#'   treatment begins for common timing mode. The corresponding year is
#'   \code{2000 + treatment_period}. Must be >= 2 and <= n_periods.
#' @param cohorts Named list or NULL. For staggered mode, maps cohort
#'   year (as character) to number of units. E.g.
#'   \code{list("2005" = 20, "2007" = 30)}.
#' @param n_never_treated Integer or NULL. Number of never-treated units
#'   in staggered mode (default NULL). When NULL, computed as
#'   \code{n_units - sum(cohorts)}. When provided, must satisfy
#'   \code{sum(cohorts) + n_never_treated == n_units}.
#' @param true_att Numeric scalar or function. If scalar, constant ATT.
#'   If function, must accept \code{(i, t, g)} and return the
#'   unit-time-cohort specific effect.
#' @param unit_fe_sd Numeric, SD of unit fixed effects (default 1.0).
#' @param time_trend_coef Numeric, coefficient for linear time trend
#'   (default 0.5).
#' @param error_sd Numeric, SD of idiosyncratic errors (default 0.5).
#' @param n_controls Integer, number of control variables to generate
#'   (default 0). Variables named x1, x2, etc.
#' @param control_effect Numeric, coefficient for each control variable
#'   (default 0.3).
#' @param trend Numeric or NULL. If non-NULL, SD of unit-specific
#'   linear trends \eqn{\delta_i \sim N(0, trend^2)} (default NULL).
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param balanced Logical. If TRUE (default), generate balanced panel.
#'   If FALSE, randomly drop observations (each with 0.8 probability
#'   of retention, minimum 2 periods per unit).
#'
#' @return A data.frame. For common timing mode: columns id, year, y,
#'   d, post (and x1, x2, ... if n_controls > 0). For staggered mode:
#'   columns id, year, y, gvar (and x1, x2, ... if n_controls > 0).
#'   NT units have \code{gvar = NA}.
#'
#' @export
#' @examples
#' # Common timing
#' ct <- simulate_panel_data(
#'   n_units = 50, n_periods = 8, n_treated = 20,
#'   treatment_period = 5, true_att = 2.0, seed = 42
#' )
#' head(ct)
#'
#' # Staggered adoption
#' stag <- simulate_panel_data(
#'   n_units = 60, n_periods = 10,
#'   cohorts = list("2005" = 20, "2007" = 20),
#'   n_never_treated = 20, true_att = 1.5, seed = 42
#' )
#' head(stag)
simulate_panel_data <- function(
    n_units = 20L,
    n_periods = 10L,
    n_treated = NULL,
    treatment_period = NULL,
    cohorts = NULL,
    n_never_treated = NULL,
    true_att = 2.0,
    unit_fe_sd = 1.0,
    time_trend_coef = 0.5,
    error_sd = 0.5,
    n_controls = 0L,
    control_effect = 0.3,
    trend = NULL,
    seed = NULL,
    balanced = TRUE
) {
  # --- Mode dispatch ---
  has_ct <- !is.null(n_treated) && !is.null(treatment_period)
  has_stag <- !is.null(cohorts)

  if (!has_ct && !has_stag) {
    stop("Must specify either (n_treated + treatment_period) for ",
         "common timing or cohorts for staggered mode.")
  }
  if (has_ct && has_stag) {
    stop("Cannot specify both (n_treated + treatment_period) and ",
         "cohorts. Choose one mode.")
  }

  if (!is.null(seed)) set.seed(seed)

  n_units <- as.integer(n_units)
  n_periods <- as.integer(n_periods)
  years <- seq.int(2001L, 2000L + n_periods)

  # --- Common DGP components ---
  # Unit fixed effects
  alpha <- rnorm(n_units, mean = 0, sd = unit_fe_sd)
  # Time trend
  lambda <- time_trend_coef * seq_len(n_periods)
  # Unit-specific trends (optional)
  delta <- if (!is.null(trend)) rnorm(n_units, mean = 0, sd = trend) else NULL
  # Errors
  eps <- matrix(rnorm(n_units * n_periods, mean = 0, sd = error_sd),
                nrow = n_units, ncol = n_periods)

  # Control variables (time-invariant: same value across all periods per unit)
  X_list <- list()
  if (n_controls > 0L) {
    for (k in seq_len(n_controls)) {
      x_unit <- rnorm(n_units)
      X_list[[paste0("x", k)]] <- matrix(
        rep(x_unit, each = 1L), nrow = n_units, ncol = n_periods
      )
    }
  }

  if (has_ct) {
    result <- .simulate_common_timing(
      n_units, n_periods, years, n_treated, treatment_period,
      true_att, alpha, lambda, delta, eps, X_list, control_effect
    )
  } else {
    result <- .simulate_staggered(
      n_units, n_periods, years, cohorts, n_never_treated,
      true_att, alpha, lambda, delta, eps, X_list, control_effect
    )
  }

  # --- Unbalanced panel (optional) ---
  if (!balanced) {
    result <- .make_unbalanced(result, n_units, n_periods)
  }

  result
}


# ============================================================================
# Internal helpers
# ============================================================================

#' @keywords internal
.simulate_common_timing <- function(
    n_units, n_periods, years, n_treated, treatment_period,
    true_att, alpha, lambda, delta, eps, X_list, control_effect
) {
  n_treated <- as.integer(n_treated)
  treatment_period <- as.integer(treatment_period)

  stopifnot(n_treated > 0L, n_treated < n_units)
  stopifnot(treatment_period >= 2L, treatment_period <= n_periods)

  # Treatment assignment: first n_treated units are treated
  d_unit <- c(rep(1L, n_treated), rep(0L, n_units - n_treated))
  # Post indicator: periods >= treatment_period
  # treatment_period is 1-based index, corresponding year = 2000 + treatment_period
  post_period <- rep(0L, n_periods)
  post_period[treatment_period:n_periods] <- 1L

  # Build panel
  id_vec <- integer(n_units * n_periods)
  year_vec <- integer(n_units * n_periods)
  y_vec <- numeric(n_units * n_periods)
  d_vec <- integer(n_units * n_periods)
  post_vec <- integer(n_units * n_periods)

  idx <- 0L
  for (i in seq_len(n_units)) {
    for (t in seq_len(n_periods)) {
      idx <- idx + 1L
      id_vec[idx] <- i
      year_vec[idx] <- years[t]
      d_vec[idx] <- d_unit[i]
      post_vec[idx] <- post_period[t]

      # D_it = d * post
      D_it <- d_unit[i] * post_period[t]

      # Treatment effect
      if (D_it == 1L) {
        if (is.function(true_att)) {
          gvar_val <- 2000L + treatment_period
          tau <- true_att(i, years[t], gvar_val)
        } else {
          tau <- true_att
        }
      } else {
        tau <- 0
      }

      # y = alpha_i + lambda_t + delta_i*t + tau*D_it + X*beta + eps
      y_val <- alpha[i] + lambda[t] + tau + eps[i, t]
      if (!is.null(delta)) {
        y_val <- y_val + delta[i] * t
      }
      # Control variables contribution
      for (nm in names(X_list)) {
        y_val <- y_val + control_effect * X_list[[nm]][i, t]
      }

      y_vec[idx] <- y_val
    }
  }

  df <- data.frame(
    id = id_vec, year = year_vec, y = y_vec,
    d = d_vec, post = post_vec
  )

  # Add control variable columns
  if (length(X_list) > 0L) {
    for (nm in names(X_list)) {
      x_col <- numeric(n_units * n_periods)
      idx2 <- 0L
      for (i in seq_len(n_units)) {
        for (t in seq_len(n_periods)) {
          idx2 <- idx2 + 1L
          x_col[idx2] <- X_list[[nm]][i, t]
        }
      }
      df[[nm]] <- x_col
    }
  }

  df
}

#' @keywords internal
.simulate_staggered <- function(
    n_units, n_periods, years, cohorts, n_never_treated,
    true_att, alpha, lambda, delta, eps, X_list, control_effect
) {
  # --- Cohort key auto-detection (design.md §6.8) ---
  # If all keys (as integers) <= n_periods, treat as period indices
  # and convert to years: year = 2000 + period_index.
  # If any key > n_periods, treat as actual years.
  raw_keys <- as.integer(names(cohorts))
  if (all(raw_keys <= n_periods)) {
    # Period indices → convert to years
    cohort_years <- 2000L + raw_keys
  } else {
    # Already actual years
    cohort_years <- raw_keys
  }
  cohort_sizes <- as.integer(unlist(cohorts))

  total_treated <- sum(cohort_sizes)

  # --- n_never_treated handling ---
  if (is.null(n_never_treated)) {
    n_never_treated <- n_units - total_treated
  } else {
    n_never_treated <- as.integer(n_never_treated)
  }
  stopifnot(total_treated + n_never_treated == n_units)

  # Assign gvar to each unit
  gvar_unit <- integer(n_units)
  unit_idx <- 1L
  for (k in seq_along(cohort_years)) {
    for (j in seq_len(cohort_sizes[k])) {
      gvar_unit[unit_idx] <- cohort_years[k]
      unit_idx <- unit_idx + 1L
    }
  }
  # Remaining units are never-treated (NA)
  if (n_never_treated > 0L) {
    gvar_unit[unit_idx:n_units] <- NA_integer_
  }

  # Build panel
  n_obs <- n_units * n_periods
  id_vec <- integer(n_obs)
  year_vec <- integer(n_obs)
  y_vec <- numeric(n_obs)
  gvar_vec <- integer(n_obs)

  idx <- 0L
  for (i in seq_len(n_units)) {
    g <- gvar_unit[i]
    for (t in seq_len(n_periods)) {
      idx <- idx + 1L
      id_vec[idx] <- i
      year_vec[idx] <- years[t]
      gvar_vec[idx] <- g

      # D_it: treated if gvar is not NA and year >= gvar
      D_it <- if (!is.na(g) && years[t] >= g) 1L else 0L

      # Treatment effect
      if (D_it == 1L) {
        if (is.function(true_att)) {
          tau <- true_att(i, years[t], g)
        } else {
          tau <- true_att
        }
      } else {
        tau <- 0
      }

      y_val <- alpha[i] + lambda[t] + tau + eps[i, t]
      if (!is.null(delta)) {
        y_val <- y_val + delta[i] * t
      }
      for (nm in names(X_list)) {
        y_val <- y_val + control_effect * X_list[[nm]][i, t]
      }

      y_vec[idx] <- y_val
    }
  }

  # gvar_vec has NA for NT units — need to preserve NA_integer_
  gvar_vec[is.na(gvar_vec)] <- NA_integer_

  df <- data.frame(
    id = id_vec, year = year_vec, y = y_vec, gvar = gvar_vec
  )

  # Add control variable columns
  if (length(X_list) > 0L) {
    for (nm in names(X_list)) {
      x_col <- numeric(n_obs)
      idx2 <- 0L
      for (i in seq_len(n_units)) {
        for (t in seq_len(n_periods)) {
          idx2 <- idx2 + 1L
          x_col[idx2] <- X_list[[nm]][i, t]
        }
      }
      df[[nm]] <- x_col
    }
  }

  df
}

#' @keywords internal
.make_unbalanced <- function(df, n_units, n_periods) {
  # Each observation retained with probability 0.8
  keep <- runif(nrow(df)) < 0.8

  # Ensure each unit has at least 2 periods
  id_col <- if ("id" %in% names(df)) "id" else "id"
  for (uid in unique(df[[id_col]])) {
    unit_rows <- which(df[[id_col]] == uid)
    unit_keep <- keep[unit_rows]
    if (sum(unit_keep) < 2L) {
      # Force keep at least 2 random rows
      forced <- sample(unit_rows, min(2L, length(unit_rows)))
      keep[forced] <- TRUE
    }
  }

  df[keep, , drop = FALSE]
}
