# transform.R
# Rolling transformation methods for lwdid package.
# Implements four transformation methods from Lee & Wooldridge (2025, 2026):
#   - demean:   Unit-specific demeaning (Procedure 2.1)
#   - detrend:  Unit-specific detrending (Procedure 3.1)
#   - demeanq:  Seasonal demeaning (common-timing path)
#   - detrendq: Seasonal detrending (common-timing path)

#' Resolve seasonal transformation config for transform_common()
#'
#' @param dt data.table containing panel data.
#' @param post Optional post indicator column name.
#' @param season_var Optional seasonal column name.
#' @param quarter Optional backward-compatible seasonal alias.
#' @param Q Optional seasonal frequency.
#' @return Named list with resolved `post`, `season_var`, and `Q`.
#' @keywords internal
.resolve_transform_common_seasonal_config <- function(
  dt,
  post = NULL,
  season_var = NULL,
  quarter = NULL,
  Q = NULL,
  rolling = NULL
) {
  seasonal_config <- attr(dt, "lwdid_seasonal_config", exact = TRUE)

  effective_season_var <- if (!is.null(season_var) || !is.null(quarter)) {
    .resolve_season_var(season_var = season_var, quarter = quarter)
  } else {
    seasonal_config$season_var %||% NULL
  }

  if (is.null(effective_season_var)) {
    example_method <- rolling %||% "demeanq"
    stop_lwdid(
      message = paste(
        "Rolling seasonal transformations require 'season_var' (or",
        "backward-compatible 'quarter') metadata.",
        sprintf(
          "Example: lwdid(..., rolling='%s', season_var='quarter', Q=4)",
          example_method
        )
      ),
      class = "lwdid_invalid_parameter",
      param = "season_var",
      value = NULL,
      allowed = "existing seasonal column name"
    )
  }

  list(
    post = post %||% seasonal_config$post %||% "post_",
    season_var = effective_season_var,
    Q = as.integer(Q %||% seasonal_config$Q %||% 4L)
  )
}

#' Build a seasonal design matrix using observed pre-period seasons
#'
#' @param season_values Seasonal values for the rows being encoded.
#' @param observed_seasons Sorted vector of seasons observed in valid pre data.
#' @return Numeric design matrix with intercept and seasonal dummies.
#' @keywords internal
.build_seasonal_design_matrix <- function(season_values, observed_seasons) {
  n_obs <- length(season_values)
  missing_season_rows <- is.na(season_values)
  design <- matrix(1, nrow = n_obs, ncol = 1L)
  colnames(design) <- "(Intercept)"

  if (length(observed_seasons) <= 1L) {
    if (any(missing_season_rows)) {
      design[missing_season_rows, ] <- NA_real_
    }
    return(design)
  }

  dummy_matrix <- vapply(
    observed_seasons[-1L],
    function(season_code) {
      as.numeric(!is.na(season_values) & season_values == season_code)
    },
    numeric(n_obs)
  )

  if (is.null(dim(dummy_matrix))) {
    dummy_matrix <- matrix(dummy_matrix, ncol = 1L)
  }

  colnames(dummy_matrix) <- paste0("season_", observed_seasons[-1L])
  design <- cbind(design, dummy_matrix)
  if (any(missing_season_rows)) {
    design[missing_season_rows, ] <- NA_real_
  }
  design
}

#' Seasonal demean fit for a single unit
#'
#' @param unit_data data.table/data.frame for one unit.
#' @param y Outcome column name.
#' @param season_var Seasonal column name.
#' @param post Effective post indicator column name.
#' @param Q Seasonal frequency (kept for API parity).
#' @param unit_label Optional unit label used in warnings.
#' @return List with `yhat`, `ydot`, `n_pre`, and fitted seasonal `coefficients`.
#' @keywords internal
.demeanq_unit <- function(unit_data, y, season_var, post, Q = 4L,
                          unit_label = NULL) {
  n_obs <- nrow(unit_data)
  post_values <- unit_data[[post]]
  unit_pre <- unit_data[post_values == 0L, ]
  valid_pre_mask <- !is.na(unit_pre[[y]]) & !is.na(unit_pre[[season_var]])
  valid_pre <- unit_pre[valid_pre_mask, ]
  n_valid <- nrow(valid_pre)
  observed_seasons <- sort(unique(valid_pre[[season_var]]))
  n_params <- length(observed_seasons)
  resolved_unit_label <- if (is.null(unit_label)) "Current unit" else as.character(unit_label)

  if (n_valid <= n_params) {
    warn_lwdid(
      message = sprintf(
        paste(
          "Unit %s has %d valid pre-treatment observation(s) for seasonal",
          "demeaning after NA filtering; seasonal transformed outcomes are",
          "set to NA."
        ),
        resolved_unit_label,
        n_valid
      ),
      class = "lwdid_small_sample",
      detail = "seasonal_demean_insufficient_valid_pre",
      action_taken = "y_trans set to NA"
    )
    return(list(
      yhat = rep(NA_real_, n_obs),
      ydot = rep(NA_real_, n_obs),
      n_pre = as.integer(n_valid)
    ))
  }

  X_pre <- .build_seasonal_design_matrix(
    season_values = valid_pre[[season_var]],
    observed_seasons = observed_seasons
  )
  fit <- tryCatch(
    stats::lm.fit(x = X_pre, y = valid_pre[[y]]),
    error = function(cnd) cnd
  )

  if (inherits(fit, "error")) {
    warn_lwdid(
      message = sprintf(
        paste(
          "Unit %s failed seasonal demean estimation before coefficients could",
          "be checked; seasonal transformed outcomes are set to NA."
        ),
        resolved_unit_label
      ),
      class = "lwdid_data",
      detail = "seasonal_demean_fit_failure",
      action_taken = "y_trans set to NA"
    )
    return(list(
      yhat = rep(NA_real_, n_obs),
      ydot = rep(NA_real_, n_obs),
      n_pre = as.integer(n_valid)
    ))
  }

  coefficients <- fit$coefficients

  if (any(!is.finite(coefficients))) {
    warn_lwdid(
      message = sprintf(
        paste(
          "Unit %s produced non-finite seasonal demean coefficients;",
          "seasonal transformed outcomes are set to NA."
        ),
        resolved_unit_label
      ),
      class = "lwdid_data",
      detail = "seasonal_demean_invalid_coefficients",
      action_taken = "y_trans set to NA"
    )
    return(list(
      yhat = rep(NA_real_, n_obs),
      ydot = rep(NA_real_, n_obs),
      n_pre = as.integer(n_valid)
    ))
  }

  X_all <- .build_seasonal_design_matrix(
    season_values = unit_data[[season_var]],
    observed_seasons = observed_seasons
  )
  fitted_values <- as.vector(X_all %*% coefficients)

  list(
    yhat = fitted_values,
    ydot = unit_data[[y]] - fitted_values,
    n_pre = as.integer(n_valid),
    coefficients = stats::setNames(as.numeric(coefficients), colnames(X_pre))
  )
}

#' Build a centered linear-trend seasonal design matrix
#'
#' @param time_values Time index values.
#' @param season_values Seasonal values for the rows being encoded.
#' @param observed_seasons Sorted vector of seasons observed in valid pre data.
#' @param t_bar_pre Mean pre-treatment time index used for centering.
#' @return Numeric design matrix with intercept, centered time, and seasonal dummies.
#' @keywords internal
.build_seasonal_trend_design_matrix <- function(
  time_values,
  season_values,
  observed_seasons,
  t_bar_pre
) {
  seasonal_design <- .build_seasonal_design_matrix(
    season_values = season_values,
    observed_seasons = observed_seasons
  )

  trend_column <- matrix(time_values - t_bar_pre, ncol = 1L)
  colnames(trend_column) <- "t_centered"

  cbind(
    seasonal_design[, 1L, drop = FALSE],
    trend_column,
    seasonal_design[, -1L, drop = FALSE]
  )
}

#' Attach common-timing post-treatment summaries to seasonal transforms
#'
#' @param data data.table containing transformed outcomes.
#' @param ivar Unit identifier column name.
#' @param tindex Time index column name.
#' @param post Original post indicator column name.
#' @param y_col Transformed outcome column name.
#' @return data.table with `ydot_postavg` and `firstpost`.
#' @keywords internal
.attach_common_timing_post_summaries <- function(
  data,
  ivar,
  tindex,
  post,
  y_col = "y_trans"
) {
  ivar_col <- ivar
  tindex_col <- tindex
  post_col <- post
  y_col_name <- y_col
  out <- data.table::copy(data.table::as.data.table(data))
  out[["ydot_postavg"]] <- NA_real_
  out[["firstpost"]] <- FALSE

  post_rows <- !is.na(out[[post_col]]) & out[[post_col]] == 1L
  if (!any(post_rows)) {
    return(out)
  }

  post_summary <- out[post_rows, .(
    ydot_postavg = mean(get(y_col_name), na.rm = TRUE),
    tpost1 = min(get(tindex_col), na.rm = TRUE)
  ), by = ivar_col]

  if (nrow(post_summary) == 0L) {
    return(out)
  }

  out[post_summary, on = ivar_col, ydot_postavg := i.ydot_postavg]
  tpost1_values <- post_summary$tpost1[match(out[[ivar_col]], post_summary[[ivar_col]])]
  out[["firstpost"]] <- !is.na(out[["ydot_postavg"]]) &
    !is.na(tpost1_values) &
    out[[tindex_col]] == tpost1_values

  out
}

#' Seasonal demean transform for multiple units
#'
#' @param data data.table containing panel data.
#' @param y Outcome column name.
#' @param ivar Unit identifier column name.
#' @param tindex Time index column name.
#' @param season_var Seasonal column name.
#' @param post Post indicator column name.
#' @param Q Seasonal frequency (kept for API parity).
#' @param exclude_pre_periods Number of trailing pre-periods excluded from fitting.
#' @return data.table with `seasonal_fit`, `y_trans`, and `n_pre`.
#' @keywords internal
.demeanq_transform <- function(data, y, ivar, tindex, season_var, post,
                               Q = 4L, exclude_pre_periods = 0L) {
  out <- data.table::copy(data.table::as.data.table(data))

  out[[".seasonal_post"]] <- .compute_effective_seasonal_post(
    data = out,
    ivar = ivar,
    tindex = tindex,
    post = post,
    exclude_pre_periods = as.integer(exclude_pre_periods)
  )
  out[["y_trans"]] <- NA_real_
  out[["seasonal_fit"]] <- NA_real_
  out[["n_pre"]] <- 0L

  for (uid in unique(out[[ivar]])) {
    unit_idx <- which(out[[ivar]] == uid)
    unit_result <- .demeanq_unit(
      unit_data = out[unit_idx],
      y = y,
      season_var = season_var,
      post = ".seasonal_post",
      Q = Q,
      unit_label = uid
    )

    seasonal_fit <- unit_result$yhat
    y_trans <- unit_result$ydot
    if (all(is.na(seasonal_fit))) {
      seasonal_fit <- rep(NaN, length(unit_idx))
    }
    if (all(is.na(y_trans))) {
      y_trans <- rep(NaN, length(unit_idx))
    }

    data.table::set(out, i = unit_idx, j = "seasonal_fit", value = seasonal_fit)
    data.table::set(out, i = unit_idx, j = "y_trans", value = y_trans)
    data.table::set(
      out,
      i = unit_idx,
      j = "n_pre",
      value = rep.int(unit_result$n_pre, length(unit_idx))
    )
  }

  out[, ".seasonal_post" := NULL]
  .attach_common_timing_post_summaries(
    data = out,
    ivar = ivar,
    tindex = tindex,
    post = post,
    y_col = "y_trans"
  )
}

#' Seasonal demean transform for common-timing designs
#'
#' @param dt data.table containing panel data.
#' @param y character, outcome variable name.
#' @param ivar character, unit identifier name.
#' @param tvar character, time index column name.
#' @param g integer, treatment onset period (unused, kept for API parity).
#' @param season_var character, seasonal column name.
#' @param Q integer, number of seasonal buckets.
#' @param post character, post indicator column name.
#' @param exclude_pre_periods integer, trailing pre-periods excluded from fitting.
#' @return data.table with seasonal demean transform columns.
#' @keywords internal
transform_demeanq <- function(dt, y, ivar, tvar, g,
                              season_var, Q = 4L, post = "post_",
                              exclude_pre_periods = 0L) {
  out <- data.table::copy(dt)

  .validate_seasonal_inputs(
    data = out,
    y = y,
    season_var = season_var,
    Q = Q,
    ivar = ivar,
    tindex = tvar,
    post = post,
    method = "demeanq",
    exclude_pre_periods = exclude_pre_periods,
    min_global_pre_periods = 1L
  )

  .demeanq_transform(
    data = out,
    y = y,
    ivar = ivar,
    tindex = tvar,
    season_var = season_var,
    post = post,
    Q = Q,
    exclude_pre_periods = exclude_pre_periods
  )
}

#' Seasonal detrend transform for common-timing designs
#'
#' @param dt data.table containing panel data.
#' @param y character, outcome variable name.
#' @param ivar character, unit identifier name.
#' @param tvar character, time index column name.
#' @param g integer, treatment onset period (unused, kept for API parity).
#' @param season_var character, seasonal column name.
#' @param Q integer, number of seasonal buckets.
#' @param post character, post indicator column name.
#' @param exclude_pre_periods integer, trailing pre-periods excluded from fitting.
#' @return data.table with seasonal detrend transform columns.
#' @keywords internal
transform_detrendq <- function(dt, y, ivar, tvar, g,
                               season_var, Q = 4L, post = "post_",
                               exclude_pre_periods = 0L) {
  out <- data.table::copy(dt)

  .validate_seasonal_inputs(
    data = out,
    y = y,
    season_var = season_var,
    Q = Q,
    ivar = ivar,
    tindex = tvar,
    post = post,
    method = "detrendq",
    exclude_pre_periods = exclude_pre_periods,
    min_global_pre_periods = 2L
  )

  out[[".seasonal_post"]] <- .compute_effective_seasonal_post(
    data = out,
    ivar = ivar,
    tindex = tvar,
    post = post,
    exclude_pre_periods = as.integer(exclude_pre_periods)
  )
  out[["y_trans"]] <- NA_real_
  out[["seasonal_fit"]] <- NA_real_
  out[["n_pre"]] <- 0L
  out[["intercept_c"]] <- NA_real_
  out[["slope"]] <- NA_real_
  out[["t_bar_pre"]] <- NA_real_

  for (uid in unique(out[[ivar]])) {
    unit_idx <- which(out[[ivar]] == uid)
    unit_data <- out[unit_idx]
    unit_pre <- unit_data[unit_data[[".seasonal_post"]] == 0L]

    t_bar_pre <- mean(unit_pre[[tvar]], na.rm = TRUE)
    t_variance <- stats::var(unit_pre[[tvar]], na.rm = TRUE)

    if (!is.finite(t_bar_pre) ||
        is.na(t_variance) ||
        t_variance < .Machine$double.eps * 100) {
      warn_lwdid(
        message = sprintf(
          paste(
            "Unit %s has degenerate pre-treatment time variation for seasonal",
            "detrending; seasonal transformed outcomes are set to NA."
          ),
          as.character(uid)
        ),
        class = "lwdid_data",
        detail = "seasonal_detrend_degenerate_time_variance",
        action_taken = "y_trans set to NA"
      )
      next
    }

    valid_pre_mask <- !is.na(unit_pre[[y]]) &
      !is.na(unit_pre[[tvar]]) &
      !is.na(unit_pre[[season_var]])
    valid_pre <- unit_pre[valid_pre_mask, ]
    n_valid <- nrow(valid_pre)
    observed_seasons <- sort(unique(valid_pre[[season_var]]))
    n_params <- length(observed_seasons) + 1L

    if (n_valid <= n_params) {
      warn_lwdid(
        message = sprintf(
          paste(
            "Unit %s has %d valid pre-treatment observation(s) for seasonal",
            "detrending after NA filtering; seasonal transformed outcomes are",
            "set to NA."
          ),
          as.character(uid),
          n_valid
        ),
        class = "lwdid_small_sample",
        detail = "seasonal_detrend_insufficient_valid_pre",
        action_taken = "y_trans set to NA"
      )
      next
    }

    X_pre <- .build_seasonal_trend_design_matrix(
      time_values = valid_pre[[tvar]],
      season_values = valid_pre[[season_var]],
      observed_seasons = observed_seasons,
      t_bar_pre = t_bar_pre
    )
    fit <- stats::lm.fit(x = X_pre, y = valid_pre[[y]])
    coefficients <- fit$coefficients

    if (any(!is.finite(coefficients))) {
      warn_lwdid(
        message = sprintf(
          paste(
            "Unit %s produced non-finite seasonal detrend coefficients;",
            "seasonal transformed outcomes are set to NA."
          ),
          as.character(uid)
        ),
        class = "lwdid_data",
        detail = "seasonal_detrend_invalid_coefficients",
        action_taken = "y_trans set to NA"
      )
      next
    }

    X_all <- .build_seasonal_trend_design_matrix(
      time_values = unit_data[[tvar]],
      season_values = unit_data[[season_var]],
      observed_seasons = observed_seasons,
      t_bar_pre = t_bar_pre
    )
    fitted_values <- as.vector(X_all %*% coefficients)

    data.table::set(out, i = unit_idx, j = "seasonal_fit", value = fitted_values)
    data.table::set(out, i = unit_idx, j = "y_trans", value = unit_data[[y]] - fitted_values)
    data.table::set(out, i = unit_idx, j = "n_pre", value = rep.int(as.integer(n_valid), length(unit_idx)))
    data.table::set(out, i = unit_idx, j = "intercept_c", value = rep.int(as.numeric(coefficients[[1L]]), length(unit_idx)))
    data.table::set(out, i = unit_idx, j = "slope", value = rep.int(as.numeric(coefficients[[2L]]), length(unit_idx)))
    data.table::set(out, i = unit_idx, j = "t_bar_pre", value = rep.int(as.numeric(t_bar_pre), length(unit_idx)))
  }

  out[, ".seasonal_post" := NULL]
  .attach_common_timing_post_summaries(
    data = out,
    ivar = ivar,
    tindex = tvar,
    post = post,
    y_col = "y_trans"
  )
}

#' Unit-specific demeaning transformation
#'
#' Implements the unit-specific demeaning transformation from
#' Lee & Wooldridge (2025, 2026) Procedure 2.1 (equations 2.11-2.12).
#' When \code{exclude_pre_periods > 0}, the pre-treatment window is
#' shortened from the end following equation 2.22 (anticipation handling).
#'
#' The transformation is computed for ALL periods (including pre-treatment),
#' not just post-treatment. Each unit's outcome is demeaned by its own
#' pre-treatment mean: \eqn{\tilde{Y}_{it} = Y_{it} - \bar{Y}_i^{pre}}.
#'
#' @param dt data.table containing panel data.
#' @param y character, name of the outcome variable column.
#' @param ivar character, name of the unit identifier column.
#' @param tvar character, name of the time variable column.
#' @param g integer, treatment time S (shared by all units under Common Timing).
#' @param exclude_pre_periods non-negative integer, number of pre-treatment
#'   periods to exclude from the end of the pre-treatment window
#'   (anticipation handling, lw2026 equation 2.22). Default is 0.
#'
#' @return A \code{data.table} with the following columns added:
#'   \describe{
#'     \item{y_trans}{Transformed outcome (Y minus pre-treatment mean), for all periods.}
#'     \item{pre_mean}{Unit-specific pre-treatment mean of the outcome.}
#'     \item{n_pre}{Integer count of valid (non-NA) pre-treatment observations per unit.}
#'   }
#'
#' @keywords internal
transform_demean <- function(dt, y, ivar, tvar, g, exclude_pre_periods = 0L) {

  # --- Step 1: Pre-treatment cutoff (calendar-based) ---
  # Exclude from the END of pre-treatment: g-1 is last pre-period,

  # then subtract exclude_pre_periods for anticipation (eq 2.22)
  pre_end <- g - 1L - exclude_pre_periods

  # --- Step 2: Panel minimum time ---
  t_min <- min(dt[[tvar]], na.rm = TRUE)

  # --- Step 3: Panel-level validation ---
  # If the pre-treatment window is empty (pre_end < t_min), no unit

  # can have valid pre-period data — this is a hard error.
  if (pre_end < t_min) {
    stop_lwdid(
      paste0(
        "No valid pre-treatment periods available for demeaning. ",
        "Pre-treatment window ends at t = ", pre_end,
        " but panel starts at t = ", t_min,
        " (g = ", g, ", exclude_pre_periods = ", exclude_pre_periods, ")."
      ),
      class = "lwdid_insufficient_pre_periods",
      n_pre = 0L,
      required = 1L,
      rolling = "demean"
    )
  }

  # --- Step 4: Per-unit pre-treatment means (vectorized via data.table) ---
  # Work on a copy to avoid modifying the input
  out <- data.table::copy(dt)

  # Subset to pre-treatment window: t_min <= t <= pre_end
  pre_data <- out[out[[tvar]] >= t_min & out[[tvar]] <= pre_end]

  # Compute per-unit mean and count of non-NA observations
  # Use .SD[[1L]] since .SDcols restricts to exactly the y column
  pre_means <- pre_data[,
    .(
      pre_mean = mean(.SD[[1L]], na.rm = TRUE),
      n_pre    = sum(!is.na(.SD[[1L]]))
    ),
    by = ivar,
    .SDcols = y
  ]

  # --- Step 5: Merge pre_means back to main data ---
  out <- merge(out, pre_means, by = ivar, all.x = TRUE)

  # --- Step 6: Handle units with no valid pre-period data ---
  # Three cases produce missing pre_mean:
  #   (A) Unit present in pre-window but all Y are NA → NaN from mean(), n_pre = 0
  #   (B) Unit absent from pre-window entirely → NA after merge
  #   (C) n_pre == 0 from either (A) or (B)
  # Normalize all to: pre_mean = NA_real_, n_pre = 0L

  has_no_pre <- is.na(out[["pre_mean"]]) | is.nan(out[["pre_mean"]]) |
                is.na(out[["n_pre"]]) | out[["n_pre"]] == 0L

  if (any(has_no_pre)) {
    # Identify affected units for the warning message
    bad_units <- unique(out[[ivar]][has_no_pre])

    warn_lwdid(
      paste0(
        length(bad_units), " unit(s) have no valid pre-treatment observations. ",
        "Their transformed outcome (y_trans) is set to NA."
      ),
      class = "lwdid_data",
      detail = "units_no_pre_periods",
      action_taken = "y_trans set to NA"
    )

    # Set pre_mean to NA_real_ and n_pre to 0L for affected rows
    data.table::set(out, i = which(has_no_pre), j = "pre_mean", value = NA_real_)
    data.table::set(out, i = which(has_no_pre), j = "n_pre",    value = 0L)
  }

  # --- Step 7: Compute transformed outcome ---
  # y_trans = y - pre_mean (NA propagates naturally for units without pre-data)
  data.table::set(out, j = "y_trans", value = out[[y]] - out[["pre_mean"]])

  out
}



#' @title Unit-specific detrending transformation
#' @description Implements lw2026 Procedure 3.1 detrending transformation
#'   (equations 3.1-3.2). For each unit, estimates a linear time trend from
#'   pre-treatment data using QR decomposition OLS, then removes the predicted
#'   trend from all periods. Time variable is centered at the pre-treatment
#'   mean for numerical stability.
#'
#'   Detrending is a generalization of demeaning: when the estimated slope
#'   B_hat = 0, detrending reduces to demeaning. Units with insufficient
#'   pre-treatment observations (n_valid < 2) are gracefully degraded to
#'   demeaning with a warning.
#'
#'   The exclude_pre_periods parameter corresponds to lw2026 Section 2.4
#'   equation 2.22 for anticipation handling.
#' @param dt data.table containing panel data
#' @param y character, name of outcome variable column
#' @param ivar character, name of unit identifier column
#' @param tvar character, name of time variable column
#' @param g integer, treatment time S (shared by all units under
#'   Common Timing)
#' @param exclude_pre_periods non-negative integer, number of
#'   pre-treatment periods to exclude (anticipation handling,
#'   corresponds to lw2026 equation 2.22, default 0)
#' @return data.table with added columns: y_trans (detrended outcome,
#'   all periods), n_pre (number of valid pre-treatment observations),
#'   slope (estimated slope B_hat, 0 for degraded units),
#'   intercept_c (centered intercept A*), t_bar_pre (pre-treatment
#'   time mean)
#' @keywords internal
transform_detrend <- function(dt, y, ivar, tvar, g,
                              exclude_pre_periods = 0L) {
  # Pre-treatment cutoff (lw2026 equation 2.22)
  pre_end <- g - 1L - exclude_pre_periods
  t_min <- min(dt[[tvar]], na.rm = TRUE)

  # Panel-level validation: detrend requires at least 2 pre-treatment
  # time points (need intercept + slope = 2 parameters)
  n_pre_periods_available <- length(
    unique(dt[[tvar]][dt[[tvar]] >= t_min & dt[[tvar]] <= pre_end])
  )
  if (n_pre_periods_available < 2L) {
    stop_lwdid(
      sprintf(
        paste("Detrend transform requires at least 2 pre-treatment",
              "periods, but only %d available.",
              "Consider using rolling='demean'"),
        n_pre_periods_available
      ),
      class = "lwdid_insufficient_pre_periods",
      n_pre = n_pre_periods_available,
      required = 2L,
      rolling = "detrend"
    )
  }

  # Informational message when pre-periods < 5 (lw2026 recommendation)
  if (n_pre_periods_available < 5L) {
    message(
      sprintf(
        paste("Note: only %d pre-treatment time points available.",
              "With few pre-periods, rolling='demean' may provide",
              "more stable estimates than detrend."),
        n_pre_periods_available
      )
    )
  }

  # Work on a copy to avoid modifying the input
  out <- data.table::copy(dt)

  # Per-unit linear trend fitting (vectorized via data.table by)
  # Each branch returns a unified 6-field list to avoid
  # data.table by inconsistency errors
  unit_trends <- out[get(tvar) >= t_min & get(tvar) <= pre_end, {
    t_vals <- .SD[[2L]]
    y_vals <- .SD[[1L]]
    valid <- !is.na(y_vals)
    t_valid <- t_vals[valid]
    y_valid <- y_vals[valid]
    n_valid <- length(y_valid)

    if (n_valid < 2L) {
      # Degradation condition 1: insufficient valid pre-period obs
      # Degrade to demean (slope = 0)
      pre_mean_fallback <- if (n_valid >= 1L) mean(y_valid) else NA_real_
      list(
        intercept_c = pre_mean_fallback,
        slope = 0,
        t_bar_pre = if (n_valid >= 1L) mean(t_valid) else NA_real_,
        n_pre = n_valid,
        degraded = TRUE,
        exact_fit = FALSE
      )
    } else {
      # Time centering for numerical stability (equation 3.1)
      t_bar <- mean(t_valid)
      t_centered <- t_valid - t_bar

      # Degradation condition 2: degenerate time variance
      # All valid pre-period time values identical -> cannot estimate slope
      t_var <- var(t_valid)
      if (is.na(t_var) || t_var < .Machine$double.eps * 100) {
        pre_mean_fallback <- mean(y_valid)
        list(
          intercept_c = pre_mean_fallback,
          slope = 0,
          t_bar_pre = t_bar,
          n_pre = n_valid,
          degraded = TRUE,
          exact_fit = FALSE
        )
      } else {
        # Mark exact fit when n_valid == 2 (df = 0)
        exact_fit_flag <- (n_valid == 2L)

        # QR decomposition OLS: Y = A* + B*(t - t_bar) (equation 3.1)
        X_pre <- cbind(1, t_centered)
        qr_fit <- qr(X_pre)
        coefs <- qr.coef(qr_fit, y_valid)

        # Degradation condition 3: NaN/Inf coefficients
        if (any(is.na(coefs)) || any(is.infinite(coefs))) {
          pre_mean_fallback <- mean(y_valid)
          list(
            intercept_c = pre_mean_fallback,
            slope = 0,
            t_bar_pre = t_bar,
            n_pre = n_valid,
            degraded = TRUE,
            exact_fit = FALSE
          )
        } else {
          list(
            intercept_c = coefs[1],  # A* = centered intercept
            slope = coefs[2],         # B = slope
            t_bar_pre = t_bar,
            n_pre = n_valid,
            degraded = FALSE,
            exact_fit = exact_fit_flag
          )
        }
      }
    }
  }, by = ivar, .SDcols = c(y, tvar)]

  # Report degraded units (warning issued outside by-loop)
  n_degraded <- sum(unit_trends$degraded)
  if (n_degraded > 0L) {
    warn_lwdid(
      sprintf(
        paste("%d unit(s) degraded to demean transform (possible",
              "causes: <2 valid pre-period obs, degenerate time",
              "variance, or ill-conditioned regression coefficients)"),
        n_degraded
      ),
      class = "lwdid_data",
      detail = "detrend_degraded_to_demean",
      action_taken = sprintf("%d units degraded", n_degraded)
    )
  }

  # Report exact fit units (n_valid == 2, df = 0)
  n_exact_fit <- sum(unit_trends$exact_fit, na.rm = TRUE)
  if (n_exact_fit > 0L) {
    warn_lwdid(
      sprintf(
        paste("%d unit(s) have exactly 2 valid pre-period obs;",
              "detrend fit has zero residual degrees of freedom",
              "(df=0), variance may be underestimated"),
        n_exact_fit
      ),
      class = "lwdid_small_sample",
      detail = "detrend_exact_fit",
      action_taken = sprintf(
        "%d units with exact fit (df=0)", n_exact_fit
      )
    )
  }

  # Merge trend parameters to main data
  out <- merge(out, unit_trends, by = ivar, all.x = TRUE)

  # Handle units with no pre-period rows at all (absent from
  # unit_trends after merge -> n_pre is NA)
  no_pre_mask <- is.na(out[["n_pre"]])
  no_pre_unit_ids <- unique(out[[ivar]][no_pre_mask])
  if (length(no_pre_unit_ids) > 0L) {
    data.table::set(out, i = which(no_pre_mask), j = "n_pre", value = 0L)
    data.table::set(out, i = which(no_pre_mask), j = "slope", value = 0)
    warn_lwdid(
      sprintf(
        paste("%d unit(s) have no pre-treatment observation rows;",
              "their transformed outcomes will be NA"),
        length(no_pre_unit_ids)
      ),
      class = "lwdid_data",
      detail = "units_no_pre_periods",
      action_taken = "y_trans set to NA"
    )
  }

  # Detrend: y_trans = Y - (A* + B*(t - t_bar_pre)) (equation 3.2)
  data.table::set(out, j = "y_trans",
    value = out[[y]] - (out[["intercept_c"]] + out[["slope"]] * (out[[tvar]] - out[["t_bar_pre"]])))

  # Clean up internal columns (keep diagnostic columns)
  cols_to_remove <- intersect(
    c("degraded", "exact_fit"), names(out)
  )
  if (length(cols_to_remove) > 0L) {
    out[, (cols_to_remove) := NULL]
  }

  out
}


#' Dispatch to the appropriate rolling transformation method
#'
#' Routes to the specific transformation function based on the
#' \code{rolling} parameter. Supports \code{"demean"}
#' (Procedure 2.1), \code{"detrend"} (Procedure 3.1), and the
#' common-timing \code{"demeanq"} seasonal extension.
#'
#' @param dt data.table containing panel data.
#' @param y character, name of the outcome variable column.
#' @param ivar character, name of the unit identifier column.
#' @param tvar character, name of the time variable column.
#' @param g integer, treatment time S.
#' @param rolling character, transformation method. One of
#'   \code{"demean"}, \code{"detrend"}, \code{"demeanq"}, \code{"detrendq"}.
#'   Default is \code{"demean"}.
#' @param exclude_pre_periods non-negative integer, number of pre-treatment
#'   periods to exclude (default 0).
#'
#' @return A \code{data.table} with transformation results (columns depend
#'   on the specific method).
#' @keywords internal
transform_common <- function(dt, y, ivar, tvar, g,
                             rolling = "demean",
                             exclude_pre_periods = 0L,
                             post = NULL,
                             season_var = NULL,
                             Q = NULL,
                             quarter = NULL) {
  switch(rolling,
    "demean"   = transform_demean(dt, y, ivar, tvar, g, exclude_pre_periods),
    "detrend"  = transform_detrend(dt, y, ivar, tvar, g, exclude_pre_periods),
    "demeanq"  = {
      seasonal_config <- .resolve_transform_common_seasonal_config(
        dt = dt,
        post = post,
        season_var = season_var,
        quarter = quarter,
        Q = Q,
        rolling = rolling
      )

      transform_demeanq(
        dt = dt,
        y = y,
        ivar = ivar,
        tvar = tvar,
        g = g,
        season_var = seasonal_config$season_var,
        Q = seasonal_config$Q,
        post = seasonal_config$post,
        exclude_pre_periods = exclude_pre_periods
      )
    },
    "detrendq" = {
      seasonal_config <- .resolve_transform_common_seasonal_config(
        dt = dt,
        post = post,
        season_var = season_var,
        quarter = quarter,
        Q = Q,
        rolling = rolling
      )

      transform_detrendq(
        dt = dt,
        y = y,
        ivar = ivar,
        tvar = tvar,
        g = g,
        season_var = seasonal_config$season_var,
        Q = seasonal_config$Q,
        post = seasonal_config$post,
        exclude_pre_periods = exclude_pre_periods
      )
    },
    stop_lwdid(
      sprintf("Invalid rolling method: '%s'", rolling),
      class = "lwdid_invalid_rolling", method = rolling)
  )
}

#' @title Precompute cohort-specific transformation statistics
#' @description Precomputes per-unit transformation statistics for each
#'   treatment cohort in staggered DiD designs. This is the first stage
#'   of the lazy transformation strategy (SA section 3.3).
#'
#'   CRITICAL: Each cohort g uses its OWN pre-treatment period
#'   \eqn{T_{pre}(g) = \{T_{min}, \ldots, g - 1 - k\}}. Different cohorts have
#'   different
#'   pre-treatment windows. NEVER share pre-treatment periods across
#'   cohorts.
#'
#'   For demean (lw2025 eq 4.11, lw2026 eq 7.3-7.4):
#'     Computes pre_mean = mean(Y) over T_pre(g) for each unit.
#'
#'   For detrend (lw2026 eq 7.5-7.6):
#'     Computes (pre_mean, slope, t_bar_pre) via centered OLS on
#'     T_pre(g) for each unit. Uses QR decomposition for numerical
#'     stability with time centering at pre-treatment mean.
#'
#' @param dt data.table containing panel data
#' @param y character, name of outcome variable column
#' @param ivar character, name of unit identifier column
#' @param tvar character, name of time variable column
#' @param cohorts integer vector, treatment cohort values (sorted,
#'   from get_cohorts())
#' @param rolling character, transformation method: "demean" or
#'   "detrend"
#' @param exclude_pre_periods non-negative integer, number of
#'   pre-treatment periods to exclude (lw2026 eq 2.22, default 0)
#' @return named list of data.tables, one per cohort (NULL entries
#'   removed for skipped cohorts). Each data.table contains:
#'   \describe{
#'     \item{ivar column}{unit identifier}
#'     \item{pre_mean}{pre-treatment mean (demean) or centered
#'       intercept A* (detrend); numerically equal}
#'     \item{n_pre}{number of non-NA pre-period observations}
#'     \item{slope}{(detrend only) estimated slope B_hat, 0 for
#'       degraded units}
#'     \item{t_bar_pre}{(detrend only) pre-treatment time mean;
#'       0 for cohort-level degraded (see design doc section 2.1.5)}
#'     \item{degraded}{(detrend only) logical, TRUE if unit or
#'       cohort was degraded to demean}
#'   }
#' @keywords internal
precompute_transforms <- function(dt, y, ivar, tvar, cohorts, rolling,
                                  exclude_pre_periods = 0L) {
  # Avoid data.table scoping conflict: if the column name equals the

  # parameter name (e.g., y="y"), get(y) inside [.data.table resolves
  # to the column (numeric) instead of the parameter (string).
  # Use distinct local aliases that cannot collide with column names.
  y_var <- y
  tvar_var <- tvar
  t_min <- min(dt[[tvar]], na.rm = TRUE)

  pre_stats <- lapply(cohorts, function(g) {
    pre_end <- g - 1L - exclude_pre_periods

    # Validate: need sufficient pre-periods
    if (pre_end < t_min) {
      warn_lwdid(
        sprintf("Cohort g=%s: pre-period range [%d, %d] is empty (exclude_pre_periods=%d), skipping cohort",
                g, t_min, pre_end, exclude_pre_periods),
        class = "lwdid_data",
        detail = "cohort_insufficient_pre_periods"
      )
      return(NULL)
    }

    if (rolling == "demean") {
      # Compute per-unit pre-period mean for cohort g
      # Cohort-specific pre-treatment: T_pre(g) = {t_min, ..., g-1-k}
      stats <- dt[get(tvar_var) >= t_min & get(tvar_var) <= pre_end,
                  .(pre_mean = mean(get(y_var), na.rm = TRUE),
                    n_pre = sum(!is.na(get(y_var)))),
                  by = ivar]

      # Handle all-NA pre-period units: NaN -> NA_real_
      stats[is.nan(pre_mean), pre_mean := NA_real_]

    } else if (rolling == "detrend") {
      # Validate: detrend needs >= 2 pre-period time points at cohort level
      n_pre_periods <- length(unique(
        dt[[tvar]][dt[[tvar]] >= t_min & dt[[tvar]] <= pre_end]
      ))
      if (n_pre_periods < 2L) {
        # Cohort-level degradation to demean
        warn_lwdid(
          sprintf("Cohort g=%s: only %d pre-period time points, detrend requires >= 2. Degrading to demean",
                  g, n_pre_periods),
          class = "lwdid_data",
          detail = "staggered_cohort_detrend_degraded"
        )
        stats <- dt[get(tvar_var) >= t_min & get(tvar_var) <= pre_end,
                    .(pre_mean = mean(get(y_var), na.rm = TRUE),
                      n_pre = sum(!is.na(get(y_var)))),
                    by = ivar]
        stats[is.nan(pre_mean), pre_mean := NA_real_]
        # t_bar_pre MUST be finite (use 0), NOT NA_real_!
        # In R: 0 * (r - NA) = NA (not 0), which would make entire
        # transform NA instead of y - pre_mean. See design doc §2.1.5.
        stats[, c("slope", "t_bar_pre", "degraded") := .(0, 0, TRUE)]
        return(stats)
      }

      # Per-unit QR decomposition OLS with 3 degradation conditions
      # Cohort-specific pre-treatment: T_pre(g) = {t_min, ..., g-1-k}
      stats <- dt[get(tvar_var) >= t_min & get(tvar_var) <= pre_end, {
        t_vals <- get(tvar_var)
        y_vals <- get(y_var)
        valid <- !is.na(y_vals)
        t_valid <- t_vals[valid]
        y_valid <- y_vals[valid]
        n_valid <- length(y_valid)

        if (n_valid < 2L) {
          # Unit-level degradation condition 1: insufficient valid obs
          pre_mean_fb <- if (n_valid >= 1L) mean(y_valid) else NA_real_
          list(pre_mean = pre_mean_fb, slope = 0,
               t_bar_pre = if (n_valid >= 1L) mean(t_valid) else NA_real_,
               n_pre = n_valid, degraded = TRUE)
        } else {
          # Time centering for numerical stability
          t_bar <- mean(t_valid)
          t_centered <- t_valid - t_bar

          # Unit-level degradation condition 2: degenerate time variance
          if (var(t_valid) < .Machine$double.eps * 100) {
            list(pre_mean = mean(y_valid), slope = 0,
                 t_bar_pre = t_bar,
                 n_pre = n_valid, degraded = TRUE)
          } else {
            X_pre <- cbind(1, t_centered)
            qr_fit <- qr(X_pre)
            coefs <- qr.coef(qr_fit, y_valid)

            # Unit-level degradation condition 3: ill-conditioned QR
            if (any(!is.finite(coefs))) {
              list(pre_mean = mean(y_valid), slope = 0,
                   t_bar_pre = t_bar,
                   n_pre = n_valid, degraded = TRUE)
            } else {
              list(pre_mean = coefs[1],
                   slope = coefs[2],
                   t_bar_pre = t_bar,
                   n_pre = n_valid, degraded = FALSE)
            }
          }
        }
      }, by = ivar]

    } else {
      stop_lwdid(sprintf("Invalid rolling method: '%s'", rolling),
                 class = "lwdid_invalid_parameter")
    }

    stats
  })

  names(pre_stats) <- as.character(cohorts)

  # Remove NULL entries (cohorts skipped due to insufficient pre-periods)
  pre_stats <- pre_stats[!vapply(pre_stats, is.null, logical(1))]

  pre_stats
}


#' @title Apply precomputed transformation to a subsample
#' @description Applies cohort-specific transformation using
#'   precomputed statistics from precompute_transforms(). This is
#'   the second stage of the lazy transformation strategy (SA
#'   section 3.3), called within the (g,r) effect estimation loop.
#'
#'   The function does NOT write back to the main data.table.
#'   It returns a temporary numeric vector for immediate use in
#'   the estimation function.
#'
#'   For demean (lw2026 eq 7.3-7.4):
#'     y_trans = Y_ir - pre_mean_ig
#'
#'   For detrend (lw2026 eq 7.5-7.6):
#'     y_trans = Y_ir - (pre_mean_ig + slope_ig * (r - t_bar_pre_ig))
#'     Mathematically equivalent to Y - (A + B*r) (see design doc
#'     section 2.1.3 equivalence proof).
#'
#'   Uses match() for unit mapping; unmatched units (not in pre_stat)
#'   receive NA in the transformed output (unbalanced panel support).
#'
#' @param y_vals numeric vector, outcome values for current subsample
#' @param unit_ids vector, unit identifiers for current subsample
#'   (same length as y_vals)
#' @param pre_stat data.table, precomputed statistics for current
#'   cohort g (one element from precompute_transforms() output)
#' @param rolling character, transformation method: "demean" or
#'   "detrend"
#' @param r integer or NULL, current estimation period. Required
#'   (non-NULL) for detrend; may be NULL for demean.
#' @param ivar character, name of unit identifier column in pre_stat.
#'   Defaults to \code{names(pre_stat)[1]} (first column).
#' @return numeric vector of transformed Y values, same length as
#'   y_vals. Units not found in pre_stat receive NA.
#' @keywords internal
apply_precomputed_transform <- function(y_vals, unit_ids, pre_stat, rolling,
                                        r = NULL,
                                        ivar = names(pre_stat)[1]) {
  # Match units via match() — lightweight, no merge overhead
  idx <- match(unit_ids, pre_stat[[ivar]])

  if (rolling == "demean") {
    y_trans <- y_vals - pre_stat$pre_mean[idx]
  } else if (rolling == "detrend") {
    # Defensive check: detrend requires current period r
    if (is.null(r)) {
      stop_lwdid(
        "apply_precomputed_transform(): detrend requires parameter r (current period), but r is NULL",
        class = "lwdid_invalid_parameter"
      )
    }
    # lw2026 eq 7.6: Y_trans = Y - (A* + B * (r - t_bar_pre))
    # Mathematically equivalent to Y - (A + B*r) — see equivalence proof
    y_trans <- y_vals - (pre_stat$pre_mean[idx] +
                         pre_stat$slope[idx] * (r - pre_stat$t_bar_pre[idx]))
  } else {
    stop_lwdid(sprintf("Invalid rolling method: '%s'", rolling),
               class = "lwdid_invalid_parameter")
  }

  y_trans
}
