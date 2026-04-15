# ============================================================================
# trend_diagnostics.R -- Trend diagnostics helpers and public entry points
#
# Foundations for trend diagnostics in Lee-Wooldridge DiD applications.
# The current implementation covers the helper layer, placebo-based parallel
# trends testing, and the first heterogeneous-trends diagnostic slice.
# ============================================================================

#' @title Trend Diagnostics Internals
#' @description Internal helpers for trend diagnostics.
#' @name trend_diagnostics
#' @keywords internal
NULL


#' Validate inputs for trend diagnostics
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param method character(1). One of visual/regression/placebo/joint.
#' @return NULL invisibly. Raises errors on failure.
#' @keywords internal
.validate_trend_test_inputs <- function(data, y, ivar, tvar, gvar, method) {
  valid_methods <- c("visual", "regression", "placebo", "joint")

  if (!is.data.frame(data)) {
    stop_lwdid(
      "data must be a data.frame",
      class = "lwdid_invalid_parameter",
      param = "data",
      value = class(data)[1L] %||% "unknown",
      allowed = "data.frame"
    )
  }

  required_columns <- c(y, ivar, tvar)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  if (!is.character(method) || length(method) != 1L ||
      !tolower(method) %in% valid_methods) {
    stop(lwdid_invalid_parameter_error(
      param = "method",
      value = method,
      allowed = valid_methods
    ))
  }

  invisible(NULL)
}

.validate_never_treated_values <- function(never_treated_values) {
  if (!is.atomic(never_treated_values) ||
      !is.numeric(never_treated_values) ||
      length(never_treated_values) < 1L) {
    invalid_value <- if (is.atomic(never_treated_values) &&
        length(never_treated_values) > 0L) {
      paste(never_treated_values, collapse = ", ")
    } else {
      class(never_treated_values)[[1L]] %||% "unknown"
    }

    stop(lwdid_invalid_parameter_error(
      param = "never_treated_values",
      value = invalid_value,
      allowed = "numeric vector"
    ))
  }

  never_treated_values
}


#' Extract sorted valid treatment cohorts
#'
#' @param data data.frame or data.table.
#' @param gvar character(1). Cohort variable.
#' @param ivar character(1). Unit identifier.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return Integer vector of sorted valid cohorts.
#' @keywords internal
.get_valid_cohorts <- function(
    data, gvar, ivar, never_treated_values = c(0, Inf)
) {
  if (is.null(gvar)) {
    return(integer(0))
  }

  never_treated_values <- .validate_never_treated_values(never_treated_values)

  missing_columns <- setdiff(c(ivar, gvar), names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  get_cohorts(data, gvar, ivar, never_treated_values = never_treated_values)
}


#' Compute the pre-treatment period range
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return integer(2). c(min_pre_periods, max_pre_periods).
#' @keywords internal
.compute_pre_period_range <- function(
    data, ivar, tvar, gvar, never_treated_values = c(0, Inf)
) {
  required_columns <- c(ivar, tvar)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  if (nrow(data) == 0L) {
    return(c(0L, 0L))
  }

  t_min <- suppressWarnings(min(data[[tvar]], na.rm = TRUE))
  t_max <- suppressWarnings(max(data[[tvar]], na.rm = TRUE))

  if (!is.finite(t_min) || !is.finite(t_max)) {
    return(c(0L, 0L))
  }

  if (is.null(gvar)) {
    treatment_period <- as.integer((t_min + t_max) %/% 2)
    n_pre <- max(as.integer(treatment_period - t_min), 0L)
    return(c(n_pre, n_pre))
  }

  cohorts <- .get_valid_cohorts(
    data,
    gvar = gvar,
    ivar = ivar,
    never_treated_values = never_treated_values
  )

  if (length(cohorts) == 0L) {
    return(c(0L, 0L))
  }

  pre_periods <- as.integer(cohorts - t_min)
  pre_periods <- pre_periods[pre_periods > 0L]

  if (length(pre_periods) == 0L) {
    return(c(0L, 0L))
  }

  c(min(pre_periods), max(pre_periods))
}


#' Check whether a panel is balanced
#'
#' @param data data.frame or data.table.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @return Logical scalar.
#' @keywords internal
.check_panel_balance <- function(data, ivar, tvar) {
  missing_columns <- setdiff(c(ivar, tvar), names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  dt <- data.table::as.data.table(data)
  n_units <- data.table::uniqueN(dt[[ivar]])
  n_periods <- data.table::uniqueN(dt[[tvar]])

  identical(as.integer(nrow(dt)), as.integer(n_units * n_periods))
}


#' Detect seasonal patterns in the outcome series
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier (unused, kept for API parity).
#' @param tvar character(1). Time variable.
#' @param threshold Numeric autocorrelation threshold.
#' @return Logical scalar.
#' @keywords internal
.detect_seasonal_patterns <- function(
    data, y, ivar, tvar, threshold = 0.1
) {
  missing_columns <- setdiff(c(y, ivar, tvar), names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  tryCatch({
    dt <- data.table::as.data.table(data)
    time_means <- dt[
      ,
      .(mean_y = mean(.SD[[1L]], na.rm = TRUE)),
      by = c(tvar),
      .SDcols = y
    ]
    data.table::setorderv(time_means, tvar)

    values <- time_means$mean_y
    values <- values[is.finite(values)]

    if (length(values) < 8L) {
      return(FALSE)
    }

    mean_value <- mean(values)
    variance_value <- mean((values - mean_value)^2)

    if (!is.finite(variance_value) || variance_value <= 0) {
      return(FALSE)
    }

    n_values <- length(values)
    for (lag in c(4L, 12L)) {
      if (n_values > lag * 2L) {
        lead_values <- values[(lag + 1L):n_values]
        lagged_values <- values[seq_len(n_values - lag)]
        autocorr <- mean(
          (lagged_values - mean_value) * (lead_values - mean_value)
        ) / variance_value

        if (is.finite(autocorr) && abs(autocorr) > threshold) {
          return(TRUE)
        }
      }
    }

    FALSE
  }, error = function(e) {
    FALSE
  })
}


#' Compute the joint F test for pre-treatment estimates
#'
#' @param pre_trend_estimates List of pre-trend estimate records.
#' @return Named list with f_stat, pvalue, and df.
#' @keywords internal
.compute_joint_f_test <- function(pre_trend_estimates) {
  if (length(pre_trend_estimates) == 0L) {
    return(list(f_stat = 0.0, pvalue = 1.0, df = c(0L, 0L)))
  }

  valid_estimates <- Filter(function(estimate) {
    att <- suppressWarnings(as.numeric(estimate$att %||% NA_real_))
    se <- suppressWarnings(as.numeric(estimate$se %||% NA_real_))

    is.finite(att) && is.finite(se) && se > 0
  }, pre_trend_estimates)

  if (length(valid_estimates) == 0L) {
    return(list(f_stat = 0.0, pvalue = 1.0, df = c(0L, 0L)))
  }

  atts <- vapply(valid_estimates, function(estimate) {
    as.numeric(estimate$att)
  }, numeric(1))
  ses <- vapply(valid_estimates, function(estimate) {
    as.numeric(estimate$se)
  }, numeric(1))
  dfs <- vapply(valid_estimates, function(estimate) {
    df_value <- estimate$df %||% estimate$df_inference %||% NA_integer_
    as.integer(df_value)
  }, integer(1))

  k <- length(valid_estimates)
  wald_stat <- sum((atts^2) / (ses^2))
  df_candidates <- dfs[is.finite(dfs) & dfs > 0L]
  df_den <- if (length(df_candidates) > 0L) {
    min(df_candidates)
  } else {
    100L
  }
  f_stat <- wald_stat / k
  pvalue <- stats::pf(f_stat, df1 = k, df2 = df_den, lower.tail = FALSE)

  list(
    f_stat = f_stat,
    pvalue = pvalue,
    df = c(as.integer(k), as.integer(df_den))
  )
}


#' Create a placebo two-by-two dataset
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1). Cohort variable.
#' @param cohort Numeric cohort identifier.
#' @param placebo_period Numeric placebo period.
#' @return data.table containing placebo indicators.
#' @keywords internal
.create_placebo_dataset <- function(
    data,
    y,
    ivar,
    tvar,
    gvar,
    cohort,
    placebo_period
) {
  dt <- data.table::as.data.table(data.table::copy(data))

  treat_mask <- !is_never_treated(dt[[gvar]]) &
    abs(dt[[gvar]] - cohort) < 1e-8
  control_mask <- is_never_treated(dt[[gvar]]) | dt[[gvar]] > placebo_period
  keep_mask <- dt[[tvar]] <= placebo_period & (treat_mask | control_mask)

  placebo_data <- dt[keep_mask]
  placebo_data[, post_placebo := as.integer(.SD[[1L]] == placebo_period), .SDcols = tvar]
  placebo_data[, treat_cohort := as.integer(
    !is_never_treated(.SD[[1L]]) & abs(.SD[[1L]] - cohort) < 1e-8
  ), .SDcols = gvar]

  placebo_data
}


#' Estimate a placebo ATT using simple two-by-two DiD
#'
#' @param placebo_data data.frame or data.table with placebo indicators.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param controls Optional control names (unused in simple fallback).
#' @param estimator Estimator label (unused in simple fallback).
#' @param n_bootstrap Integer bootstrap count placeholder.
#' @return Named list with att, se, n_treated, n_control, and df.
#' @keywords internal
.estimate_placebo_att <- function(
    placebo_data,
    y,
    ivar,
    tvar,
    controls = NULL,
    estimator = "ra",
    n_bootstrap = 0L
) {
  dt <- data.table::as.data.table(data.table::copy(placebo_data))

  if (nrow(dt) == 0L) {
    return(list(
      att = NA_real_,
      se = NA_real_,
      n_treated = 0L,
      n_control = 0L,
      df = 0L
    ))
  }

  treat_mask <- dt$treat_cohort == 1L
  control_mask <- dt$treat_cohort == 0L
  post_mask <- dt$post_placebo == 1L
  pre_mask <- dt$post_placebo == 0L

  n_treated <- as.integer(sum(treat_mask, na.rm = TRUE) %/% 2L)
  n_control <- as.integer(sum(control_mask, na.rm = TRUE) %/% 2L)

  if (n_treated < 1L || n_control < 1L) {
    return(list(
      att = NA_real_,
      se = NA_real_,
      n_treated = n_treated,
      n_control = n_control,
      df = 0L
    ))
  }

  y_treat_post <- suppressWarnings(mean(dt[[y]][treat_mask & post_mask], na.rm = TRUE))
  y_treat_pre <- suppressWarnings(mean(dt[[y]][treat_mask & pre_mask], na.rm = TRUE))
  y_control_post <- suppressWarnings(mean(dt[[y]][control_mask & post_mask], na.rm = TRUE))
  y_control_pre <- suppressWarnings(mean(dt[[y]][control_mask & pre_mask], na.rm = TRUE))

  att <- (y_treat_post - y_treat_pre) - (y_control_post - y_control_pre)
  n_total <- nrow(dt)

  if (n_total > 4L) {
    var_treat_post <- stats::var(dt[[y]][treat_mask & post_mask], na.rm = TRUE)
    var_treat_pre <- stats::var(dt[[y]][treat_mask & pre_mask], na.rm = TRUE)
    var_control_post <- stats::var(dt[[y]][control_mask & post_mask], na.rm = TRUE)
    var_control_pre <- stats::var(dt[[y]][control_mask & pre_mask], na.rm = TRUE)

    n_tp <- sum(treat_mask & post_mask, na.rm = TRUE)
    n_tr <- sum(treat_mask & pre_mask, na.rm = TRUE)
    n_cp <- sum(control_mask & post_mask, na.rm = TRUE)
    n_cr <- sum(control_mask & pre_mask, na.rm = TRUE)

    se_sq <- 0
    if (n_tp > 0L && is.finite(var_treat_post)) {
      se_sq <- se_sq + var_treat_post / n_tp
    }
    if (n_tr > 0L && is.finite(var_treat_pre)) {
      se_sq <- se_sq + var_treat_pre / n_tr
    }
    if (n_cp > 0L && is.finite(var_control_post)) {
      se_sq <- se_sq + var_control_post / n_cp
    }
    if (n_cr > 0L && is.finite(var_control_pre)) {
      se_sq <- se_sq + var_control_pre / n_cr
    }

    se <- if (se_sq > 0) sqrt(se_sq) else NA_real_
    df <- as.integer(n_total - 4L)
  } else {
    se <- NA_real_
    df <- 0L
  }

  list(
    att = att,
    se = se,
    n_treated = as.integer(n_treated),
    n_control = as.integer(n_control),
    df = df
  )
}


#' Estimate placebo ATTs via simple two-by-two DiD fallback
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1). Cohort variable.
#' @param controls Optional control names.
#' @param estimator Estimator label.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param alpha Numeric significance level.
#' @param warnings_list Character vector of warnings.
#' @param T_min Minimum period in the sample.
#' @param cohorts Integer vector of cohorts.
#' @param n_bootstrap Integer bootstrap count placeholder.
#' @return List of placebo estimate records.
#' @keywords internal
.estimate_placebo_with_simple_did <- function(
    data,
    y,
    ivar,
    tvar,
    gvar,
    controls = NULL,
    estimator = "ra",
    never_treated_values = c(0, Inf),
    alpha = 0.05,
    warnings_list = character(0),
    T_min,
    cohorts,
    n_bootstrap = 0L
) {
  dt <- data.table::as.data.table(data.table::copy(data))
  estimates <- list()

  if (length(cohorts) == 0L) {
    return(estimates)
  }

  for (cohort in cohorts) {
    n_pre <- as.integer(cohort - T_min)
    if (n_pre < 2L) {
      next
    }

    placebo_periods <- seq.int(as.integer(T_min + 1L), as.integer(cohort - 2L))
    for (placebo_period in placebo_periods) {
      placebo_data <- .create_placebo_dataset(
        data = dt,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        cohort = cohort,
        placebo_period = placebo_period
      )

      att_result <- .estimate_placebo_att(
        placebo_data = placebo_data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        controls = controls,
        estimator = estimator,
        n_bootstrap = n_bootstrap
      )

      if (!is.finite(att_result$att) ||
          !is.finite(att_result$se) ||
          att_result$se <= 0) {
        next
      }

      t_stat <- att_result$att / att_result$se
      df <- max(as.integer(att_result$df), 1L)
      pvalue <- 2 * stats::pt(-abs(t_stat), df = df)
      t_crit <- stats::qt(1 - alpha / 2, df = df)

      estimates[[length(estimates) + 1L]] <- list(
        cohort = as.integer(cohort),
        period = as.integer(placebo_period),
        event_time = as.integer(placebo_period - cohort),
        att = as.numeric(att_result$att),
        se = as.numeric(att_result$se),
        pvalue = as.numeric(pvalue),
        df = as.integer(att_result$df),
        is_anchor = FALSE,
        n_treated = as.integer(att_result$n_treated),
        n_control = as.integer(att_result$n_control),
        fallback_source = "simple_did",
        t_stat = as.numeric(t_stat),
        ci_lower = as.numeric(att_result$att - t_crit * att_result$se),
        ci_upper = as.numeric(att_result$att + t_crit * att_result$se)
      )
    }
  }

  estimates
}


#' Coerce placebo estimates into a tabular container
#'
#' @param pre_effects data.frame/data.table or list of estimate records.
#' @return data.table with one row per estimate.
#' @keywords internal
.coerce_placebo_effects <- function(pre_effects) {
  if (is.null(pre_effects)) {
    return(NULL)
  }

  if (inherits(pre_effects, "data.table") || inherits(pre_effects, "data.frame")) {
    return(data.table::as.data.table(data.table::copy(pre_effects)))
  }

  if (is.list(pre_effects)) {
    if (length(pre_effects) == 0L) {
      return(data.table::data.table())
    }

    return(data.table::rbindlist(pre_effects, fill = TRUE))
  }

  stop("Unsupported placebo estimate container.", call. = FALSE)
}


#' Identify errors indicating the staggered placebo helper is unavailable
#'
#' @param error condition object.
#' @return Logical scalar.
#' @keywords internal
.is_staggered_module_unavailable_error <- function(error) {
  if (inherits(error, "lwdid_staggered_unavailable")) {
    return(TRUE)
  }

  error_message <- conditionMessage(error)
  unavailable_patterns <- c(
    "could not find function \"estimate_pre_treatment_staggered\"",
    "could not find function \"transform_staggered_demean_pre\"",
    "could not find function \"transform_staggered_detrend_pre\"",
    "attempt to apply non-function",
    "Staggered module not available"
  )

  any(vapply(
    unavailable_patterns,
    function(pattern) grepl(pattern, error_message, fixed = TRUE),
    logical(1)
  ))
}


#' Estimate placebo effects with simple-DiD fallback when staggered helper is unavailable
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1). Cohort variable.
#' @param controls Optional control names.
#' @param estimator Estimator label.
#' @param rolling character(1). Rolling transform choice.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param alpha Numeric significance level.
#' @param warnings_list Character vector of warnings.
#' @param T_min Minimum period in the sample.
#' @param cohorts Integer vector of cohorts.
#' @param n_bootstrap Integer bootstrap count placeholder.
#' @return List with `pre_effects` and `warnings_list`.
#' @keywords internal
.estimate_placebo_with_staggered_fallback <- function(
    data,
    y,
    ivar,
    tvar,
    gvar,
    controls = NULL,
    estimator = "ra",
    rolling = "demean",
    never_treated_values = c(0, Inf),
    alpha = 0.05,
    warnings_list = character(0),
    T_min,
    cohorts,
    n_bootstrap = 0L
) {
  tryCatch(
    list(
      pre_effects = .coerce_placebo_effects(estimate_pre_treatment_staggered(
        data = data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        rolling = rolling,
        controls = controls,
        estimator = estimator,
        control_group = "not_yet_treated",
        alpha = alpha,
        include_earliest_pre_period = TRUE
      )),
      warnings_list = warnings_list
    ),
    error = function(error) {
      if (!.is_staggered_module_unavailable_error(error)) {
        stop(error)
      }

      fallback_warning <- paste(
        "Staggered module not available. Using simple 2x2 DiD for placebo test.",
        "This may produce conservative standard errors. For proper implementation,",
        "ensure the staggered module is importable."
      )

      list(
        pre_effects = .coerce_placebo_effects(.estimate_placebo_with_simple_did(
          data = data,
          y = y,
          ivar = ivar,
          tvar = tvar,
          gvar = gvar,
          controls = controls,
          estimator = estimator,
          never_treated_values = never_treated_values,
          alpha = alpha,
          warnings_list = warnings_list,
          T_min = T_min,
          cohorts = cohorts,
          n_bootstrap = n_bootstrap
        )),
        warnings_list = c(warnings_list, fallback_warning)
      )
    }
  )
}


#' Return a default cohort-trend estimate record
#'
#' @param n_units Integer unit count.
#' @param n_pre_periods Integer pre-period count.
#' @return Named list with NaN trend statistics.
#' @keywords internal
.default_cohort_trend <- function(n_units = 0L, n_pre_periods = 0L) {
  list(
    intercept = NaN,
    intercept_se = NaN,
    slope = NaN,
    slope_se = NaN,
    slope_pvalue = 1.0,
    n_units = as.integer(n_units),
    n_pre_periods = as.integer(n_pre_periods),
    r_squared = NaN,
    residual_std = NaN
  )
}


#' Build an optional control matrix for trend regressions
#'
#' @param data data.frame or data.table.
#' @param controls Optional character vector of control names.
#' @return Numeric matrix or NULL.
#' @keywords internal
.build_trend_control_matrix <- function(data, controls = NULL) {
  if (is.null(controls) || length(controls) == 0L) {
    return(NULL)
  }

  control_df <- as.data.frame(data[, controls, with = FALSE])
  if (ncol(control_df) == 0L) {
    return(NULL)
  }

  stats::model.matrix(~ . - 1, data = control_df)
}


#' Estimate a pooled cohort-specific linear pre-trend
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param controls Optional character vector of control names.
#' @return Named list with pooled OLS trend statistics.
#' @keywords internal
.estimate_cohort_trend <- function(
    data,
    y,
    ivar,
    tvar,
    controls = NULL
) {
  required_columns <- c(y, ivar, tvar, controls %||% character(0))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  n_units <- data.table::uniqueN(dt[[ivar]])
  n_pre_periods <- data.table::uniqueN(dt[[tvar]])

  if (nrow(dt) < 3L) {
    return(.default_cohort_trend(n_units = n_units, n_pre_periods = n_pre_periods))
  }

  y_vals <- as.numeric(dt[[y]])
  t_vals <- as.numeric(dt[[tvar]])
  controls_matrix <- .build_trend_control_matrix(dt, controls)
  X <- cbind(`(Intercept)` = 1, t_centered = t_vals - mean(t_vals, na.rm = TRUE))
  if (!is.null(controls_matrix)) {
    X <- cbind(X, controls_matrix)
  }

  keep_mask <- stats::complete.cases(cbind(y_vals, X))
  y_vals <- y_vals[keep_mask]
  X <- X[keep_mask, , drop = FALSE]
  t_vals <- t_vals[keep_mask]

  if (length(y_vals) < 3L || nrow(X) < 3L) {
    return(.default_cohort_trend(n_units = n_units, n_pre_periods = n_pre_periods))
  }

  fit <- tryCatch(
    stats::lm.fit(x = X, y = y_vals),
    error = function(...) NULL
  )
  if (is.null(fit) || fit$rank < ncol(X)) {
    return(.default_cohort_trend(n_units = n_units, n_pre_periods = n_pre_periods))
  }

  df_resid <- nrow(X) - ncol(X)
  if (df_resid <= 0L) {
    return(.default_cohort_trend(n_units = n_units, n_pre_periods = n_pre_periods))
  }

  XtX_inv <- tryCatch(
    solve(crossprod(X)),
    error = function(...) NULL
  )
  if (is.null(XtX_inv)) {
    return(.default_cohort_trend(n_units = n_units, n_pre_periods = n_pre_periods))
  }

  residual_ss <- sum(fit$residuals^2)
  sigma2 <- residual_ss / df_resid
  vcov_beta <- sigma2 * XtX_inv
  coefficients <- as.numeric(fit$coefficients)
  slope_se <- sqrt(vcov_beta[2L, 2L])
  slope_t <- coefficients[2L] / slope_se
  slope_pvalue <- 2 * stats::pt(-abs(slope_t), df = df_resid)
  ss_tot <- sum((y_vals - mean(y_vals))^2)

  list(
    intercept = coefficients[1L],
    intercept_se = sqrt(vcov_beta[1L, 1L]),
    slope = coefficients[2L],
    slope_se = slope_se,
    slope_pvalue = slope_pvalue,
    n_units = as.integer(n_units),
    n_pre_periods = as.integer(data.table::uniqueN(t_vals)),
    r_squared = if (ss_tot > 0) 1 - residual_ss / ss_tot else 0.0,
    residual_std = sqrt(sigma2)
  )
}


#' Return a default heterogeneous-trend test result
#'
#' @return Named list with zeroed F-test statistics.
#' @keywords internal
.default_trend_heterogeneity_test <- function() {
  list(
    f_stat = 0.0,
    pvalue = 1.0,
    df_num = 0L,
    df_den = 0L,
    reject_null = FALSE
  )
}


#' Test whether cohort-specific linear trends are homogeneous
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1). Cohort variable.
#' @param controls Optional character vector of control names.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param alpha Numeric significance level.
#' @return Named list with F-test statistics.
#' @keywords internal
.test_trend_heterogeneity <- function(
    data,
    y,
    ivar,
    tvar,
    gvar,
    controls = NULL,
    never_treated_values = c(0, Inf),
    alpha = 0.05
) {
  required_columns <- c(y, ivar, tvar, gvar, controls %||% character(0))
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  custom_never <- setdiff(never_treated_values, c(0, Inf))
  if (length(custom_never) > 0L) {
    custom_mask <- dt[[gvar]] %in% custom_never
    custom_mask[is.na(custom_mask)] <- FALSE
    if (any(custom_mask)) {
      data.table::set(dt, i = which(custom_mask), j = gvar, value = Inf)
    }
  }

  cohorts <- .get_valid_cohorts(
    dt,
    gvar = gvar,
    ivar = ivar,
    never_treated_values = never_treated_values
  )
  if (length(cohorts) < 2L) {
    return(.default_trend_heterogeneity_test())
  }

  g_values <- dt[[gvar]]
  pre_mask <- is_never_treated(g_values) | dt[[tvar]] < g_values
  pre_data <- dt[pre_mask]
  if (nrow(pre_data) < 10L) {
    return(.default_trend_heterogeneity_test())
  }

  y_vals <- as.numeric(pre_data[[y]])
  t_centered <- as.numeric(pre_data[[tvar]] - mean(pre_data[[tvar]], na.rm = TRUE))
  controls_matrix <- .build_trend_control_matrix(pre_data, controls)
  x_restricted <- cbind(`(Intercept)` = 1, t_centered = t_centered)
  if (!is.null(controls_matrix)) {
    x_restricted <- cbind(x_restricted, controls_matrix)
  }

  x_full <- x_restricted
  for (cohort in cohorts[-1L]) {
    cohort_indicator <- as.numeric(abs(pre_data[[gvar]] - cohort) < LWDID_COHORT_FLOAT_TOLERANCE)
    cohort_indicator[is.na(cohort_indicator)] <- 0
    interaction_column <- matrix(
      cohort_indicator * t_centered,
      ncol = 1L,
      dimnames = list(NULL, sprintf("cohort_%s_trend", cohort))
    )
    x_full <- cbind(x_full, interaction_column)
  }

  keep_mask <- stats::complete.cases(cbind(y_vals, x_full))
  y_vals <- y_vals[keep_mask]
  x_restricted <- x_restricted[keep_mask, , drop = FALSE]
  x_full <- x_full[keep_mask, , drop = FALSE]

  if (length(y_vals) < 10L) {
    return(.default_trend_heterogeneity_test())
  }

  fit_restricted <- tryCatch(
    stats::lm.fit(x = x_restricted, y = y_vals),
    error = function(...) NULL
  )
  fit_full <- tryCatch(
    stats::lm.fit(x = x_full, y = y_vals),
    error = function(...) NULL
  )
  if (is.null(fit_restricted) ||
      is.null(fit_full) ||
      fit_restricted$rank < ncol(x_restricted) ||
      fit_full$rank < ncol(x_full)) {
    return(.default_trend_heterogeneity_test())
  }

  df_num <- length(cohorts) - 1L
  df_den <- length(y_vals) - ncol(x_full)
  if (df_num <= 0L || df_den <= 0L) {
    return(.default_trend_heterogeneity_test())
  }

  ssr_restricted <- sum(fit_restricted$residuals^2)
  ssr_full <- sum(fit_full$residuals^2)
  if (!is.finite(ssr_restricted) || !is.finite(ssr_full) || ssr_full <= 0) {
    return(.default_trend_heterogeneity_test())
  }

  numerator <- max(ssr_restricted - ssr_full, 0) / df_num
  denominator <- ssr_full / df_den
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(.default_trend_heterogeneity_test())
  }

  f_stat <- numerator / denominator
  pvalue <- stats::pf(f_stat, df1 = df_num, df2 = df_den, lower.tail = FALSE)
  if (!is.finite(f_stat) || !is.finite(pvalue)) {
    return(.default_trend_heterogeneity_test())
  }

  list(
    f_stat = f_stat,
    pvalue = pvalue,
    df_num = as.integer(df_num),
    df_den = as.integer(df_den),
    reject_null = isTRUE(pvalue < alpha)
  )
}


#' Compute pairwise cohort-trend differences
#'
#' @param trend_by_cohort List of cohort-trend estimate records.
#' @param control_group_trend Optional control-group trend record.
#' @param alpha Numeric significance level.
#' @return List of pairwise trend-difference records.
#' @keywords internal
.compute_pairwise_trend_differences <- function(
    trend_by_cohort,
    control_group_trend = NULL,
    alpha = 0.05
) {
  all_trends <- trend_by_cohort
  if (!is.null(control_group_trend)) {
    all_trends <- c(all_trends, list(control_group_trend))
  }

  if (length(all_trends) < 2L) {
    return(list())
  }

  differences <- list()

  for (i in seq_len(length(all_trends) - 1L)) {
    for (j in seq.int(i + 1L, length(all_trends))) {
      trend_1 <- all_trends[[i]]
      trend_2 <- all_trends[[j]]
      slope_1 <- as.numeric(trend_1$slope %||% NaN)
      slope_2 <- as.numeric(trend_2$slope %||% NaN)
      slope_se_1 <- as.numeric(trend_1$slope_se %||% NaN)
      slope_se_2 <- as.numeric(trend_2$slope_se %||% NaN)

      if (!is.finite(slope_1) ||
          !is.finite(slope_2) ||
          !is.finite(slope_se_1) ||
          !is.finite(slope_se_2) ||
          slope_se_1 <= 0 ||
          slope_se_2 <= 0) {
        next
      }

      slope_diff <- slope_1 - slope_2
      slope_diff_se <- sqrt(slope_se_1^2 + slope_se_2^2)

      if (slope_diff_se > 0) {
        t_stat <- slope_diff / slope_diff_se
        df_num <- (slope_se_1^2 + slope_se_2^2)^2
        df_den <- (slope_se_1^4 / max(as.integer(trend_1$n_units %||% 0L) - 1L, 1L)) +
          (slope_se_2^4 / max(as.integer(trend_2$n_units %||% 0L) - 1L, 1L))
        df <- if (df_den > 0) max(as.integer(floor(df_num / df_den)), 1L) else 1L
        pvalue <- 2 * stats::pt(-abs(t_stat), df = df)
      } else {
        t_stat <- 0.0
        pvalue <- 1.0
        df <- 1L
      }

      differences[[length(differences) + 1L]] <- list(
        cohort_1 = as.integer(trend_1$cohort),
        cohort_2 = as.integer(trend_2$cohort),
        slope_1 = slope_1,
        slope_2 = slope_2,
        slope_diff = slope_diff,
        slope_diff_se = slope_diff_se,
        t_stat = t_stat,
        pvalue = pvalue,
        df = as.integer(df)
      )
    }
  }

  differences
}


#' Parallel trends diagnostics via placebo pre-trend estimates
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param controls Optional character vector of controls.
#' @param method character(1). Currently supports placebo/joint.
#' @param estimator character(1). Estimator type.
#' @param alpha Numeric significance level.
#' @param n_bootstrap Integer bootstrap count placeholder.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param rolling character(1). Rolling transform choice.
#' @param verbose Logical scalar.
#' @return `lwdid_parallel_trends` object.
#' @export
lwdid_test_parallel_trends <- function(
    data,
    y,
    ivar,
    tvar,
    gvar = NULL,
    controls = NULL,
    method = "placebo",
    estimator = "ra",
    alpha = 0.05,
    n_bootstrap = 0L,
    never_treated_values = c(0, Inf),
    rolling = "demean",
    verbose = TRUE
) {
  .validate_trend_test_inputs(
    data = data,
    y = y,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    method = method
  )
  never_treated_values <- .validate_never_treated_values(never_treated_values)

  method <- tolower(method)
  if (!method %in% c("placebo", "joint")) {
    stop(lwdid_invalid_parameter_error(
      param = "method",
      value = method,
      allowed = c("placebo", "joint")
    ))
  }

  dt <- data.table::as.data.table(data.table::copy(data))
  warnings_list <- character(0)
  t_min <- suppressWarnings(min(dt[[tvar]], na.rm = TRUE))

  if (!is.null(gvar)) {
    custom_never <- setdiff(never_treated_values, c(0, Inf))
    if (length(custom_never) > 0L) {
      never_mask <- dt[[gvar]] %in% custom_never
      never_mask[is.na(never_mask)] <- FALSE
      if (any(never_mask)) {
        data.table::set(dt, i = which(never_mask), j = gvar, value = Inf)
      }
    }
  }

  if (is.null(gvar)) {
    unique_units <- unique(dt[[ivar]])
    n_units <- length(unique_units)
    if (n_units < 2L) {
      stop_lwdid(
        "Common-timing fallback requires at least two units.",
        class = "lwdid_insufficient_data",
        n = n_units,
        n_treated = 0L,
        n_control = n_units
      )
    }

    t_max <- suppressWarnings(max(dt[[tvar]], na.rm = TRUE))
    treatment_period <- as.integer((t_min + t_max) %/% 2)
    cohorts <- as.integer(treatment_period)
    treated_units <- unique_units[seq_len(n_units %/% 2L)]

    warnings_list <- c(
      warnings_list,
      sprintf(
        "No gvar specified. Assuming common timing with treatment at period %d.",
        treatment_period
      )
    )

    dt[, `_gvar_dummy` := ifelse(
      .SD[[1L]] %in% treated_units,
      treatment_period,
      Inf
    ), .SDcols = ivar]
    dt[, `.trend_treat_dummy` := as.integer(.SD[[1L]] %in% treated_units),
      .SDcols = ivar]

    staggered_result <- .estimate_placebo_with_staggered_fallback(
      data = dt,
      y = y,
      ivar = ivar,
      tvar = tvar,
      gvar = "_gvar_dummy",
      controls = controls,
      rolling = rolling,
      estimator = estimator,
      never_treated_values = never_treated_values,
      alpha = alpha,
      warnings_list = warnings_list,
      T_min = t_min,
      cohorts = cohorts,
      n_bootstrap = n_bootstrap
    )
    pre_effects <- staggered_result$pre_effects
    warnings_list <- staggered_result$warnings_list
    gvar_used <- "_gvar_dummy"
  } else {
    cohorts <- .get_valid_cohorts(
      dt,
      gvar = gvar,
      ivar = ivar,
      never_treated_values = never_treated_values
    )

    if (length(cohorts) == 0L) {
      stop_lwdid(
        "No valid treatment cohorts found in data.",
        class = "lwdid_invalid_staggered_data",
        gvar = gvar,
        invalid_values = unique(dt[[gvar]]),
        detail = "No valid treatment cohorts remain after excluding never-treated markers."
      )
    }

    staggered_result <- .estimate_placebo_with_staggered_fallback(
      data = dt,
      y = y,
      ivar = ivar,
      tvar = tvar,
      gvar = gvar,
      estimator = estimator,
      controls = controls,
      rolling = rolling,
      never_treated_values = never_treated_values,
      alpha = alpha,
      warnings_list = warnings_list,
      T_min = t_min,
      cohorts = cohorts,
      n_bootstrap = n_bootstrap
    )
    pre_effects <- staggered_result$pre_effects
    warnings_list <- staggered_result$warnings_list
    gvar_used <- gvar
  }

  unit_level_gvar <- if (!is.null(gvar_used) && gvar_used %in% names(dt)) {
    get_unit_level_gvar(dt, gvar_used, ivar)
  } else {
    NULL
  }

  pre_trend_estimates <- if (!is.null(pre_effects) && nrow(pre_effects) > 0L) {
    lapply(seq_len(nrow(pre_effects)), function(index) {
      row <- pre_effects[index, , drop = FALSE]
      period_value <- as.integer(row$period[[1L]])
      public_control_count <- as.integer(row$n_control[[1L]])
      simple_fallback <- "fallback_source" %in% names(row) &&
        identical(row$fallback_source[[1L]], "simple_did")
      public_df <- if ("df" %in% names(row)) as.integer(row$df[[1L]]) else NA_integer_

      if (!simple_fallback && !is.null(unit_level_gvar)) {
        cohort_values <- unit_level_gvar[[gvar_used]]
        control_mask <- is_never_treated(cohort_values) | cohort_values > period_value
        public_control_count <- as.integer(sum(control_mask, na.rm = TRUE))
        public_df <- as.integer(max(row$n_treated[[1L]] + public_control_count - 2L, 0L))
      }

      list(
        cohort = row$cohort[[1L]],
        period = row$period[[1L]],
        event_time = row$event_time[[1L]],
        att = row$att[[1L]],
        se = row$se[[1L]],
        pvalue = row$pvalue[[1L]],
        df = public_df,
        is_anchor = isTRUE(row$is_anchor[[1L]]),
        n_treated = row$n_treated[[1L]],
        n_control = public_control_count
      )
    })
  } else {
    list()
  }
  pre_trend_estimates <- Filter(function(estimate) {
    se_value <- suppressWarnings(as.numeric(estimate$se %||% NA_real_))

    !isTRUE(estimate$is_anchor) &&
      is.finite(se_value) &&
      se_value > 0
  }, pre_trend_estimates)

  joint_test <- if (length(pre_trend_estimates) > 0L) {
    .compute_joint_f_test(pre_trend_estimates)
  } else {
    list(f_stat = NaN, pvalue = 1.0, df = c(0L, 0L))
  }
  significant_estimates <- Filter(function(estimate) {
    pvalue <- suppressWarnings(as.numeric(estimate$pvalue %||% NA_real_))

    is.finite(pvalue) && pvalue < alpha
  }, pre_trend_estimates)

  if (length(pre_trend_estimates) == 0L) {
    warnings_list <- c(
      warnings_list,
      "No valid pre-treatment estimates computed. Cannot perform joint test."
    )
  }

  reject_null <- is.finite(joint_test$pvalue) && joint_test$pvalue < alpha
  n_significant <- length(significant_estimates)
  n_total <- length(pre_trend_estimates)

  if (reject_null || (n_total > 0L && n_significant > n_total * 0.2)) {
    recommendation <- "detrend"
    recommendation_reason <- sprintf(
      paste(
        "Parallel trends assumption appears violated:",
        "joint F-test p=%.4f, %d/%d pre-treatment estimates significant.",
        "Detrending removes unit-specific linear trends under Assumption CHT."
      ),
      joint_test$pvalue,
      n_significant,
      n_total
    )
  } else {
    recommendation <- "demean"
    recommendation_reason <- sprintf(
      paste(
        "Parallel trends assumption appears to hold:",
        "joint F-test p=%.4f, %d/%d pre-treatment estimates significant.",
        "Demeaning is more efficient when PT holds."
      ),
      joint_test$pvalue,
      n_significant,
      n_total
    )
  }

  result <- structure(
    list(
      method = method,
      reject_null = reject_null,
      pvalue = joint_test$pvalue,
      test_statistic = joint_test$f_stat,
      pre_trend_estimates = pre_trend_estimates,
      joint_f_stat = joint_test$f_stat,
      joint_pvalue = joint_test$pvalue,
      joint_df = joint_test$df,
      recommendation = recommendation,
      recommendation_reason = recommendation_reason,
      warnings = warnings_list,
      figure = NULL,
      cohorts = as.integer(cohorts),
      gvar = gvar_used,
      rolling = rolling,
      alpha = alpha,
      n_bootstrap = as.integer(n_bootstrap)
    ),
    class = c("lwdid_parallel_trends", "list")
  )

  if (isTRUE(verbose)) {
    message(result$recommendation_reason)
  }

  result
}


#' Diagnose heterogeneous cohort trends
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param controls Optional character vector of controls.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param include_control_group Logical scalar.
#' @param alpha Numeric significance level.
#' @param verbose Logical scalar.
#' @return `lwdid_heterogeneous_trends` object.
#' @export
lwdid_diagnose_heterogeneous_trends <- function(
    data,
    y,
    ivar,
    tvar,
    gvar = NULL,
    controls = NULL,
    never_treated_values = c(0, Inf),
    include_control_group = TRUE,
    alpha = 0.05,
    verbose = TRUE
) {
  required_columns <- c(y, ivar, tvar, controls %||% character(0))
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }
  never_treated_values <- .validate_never_treated_values(never_treated_values)

  dt <- data.table::as.data.table(data.table::copy(data))
  if (!is.null(gvar)) {
    custom_never <- setdiff(never_treated_values, c(0, Inf))
    if (length(custom_never) > 0L) {
      custom_mask <- dt[[gvar]] %in% custom_never
      custom_mask[is.na(custom_mask)] <- FALSE
      if (any(custom_mask)) {
        data.table::set(dt, i = which(custom_mask), j = gvar, value = Inf)
      }
    }
  }

  if (is.null(gvar)) {
    t_min <- suppressWarnings(min(dt[[tvar]], na.rm = TRUE))
    t_max <- suppressWarnings(max(dt[[tvar]], na.rm = TRUE))
    cohorts <- as.integer((t_min + t_max) %/% 2)
  } else {
    cohorts <- .get_valid_cohorts(
      dt,
      gvar = gvar,
      ivar = ivar,
      never_treated_values = never_treated_values
    )
  }

  trend_by_cohort <- list()

  for (cohort_value in cohorts) {
    cohort_data <- if (is.null(gvar)) {
      dt
    } else {
      dt[abs(get(gvar) - cohort_value) < LWDID_COHORT_FLOAT_TOLERANCE]
    }
    pre_data <- cohort_data[cohort_data[[tvar]] < cohort_value]
    if (nrow(pre_data) < 3L) {
      next
    }

    trend_estimate <- .estimate_cohort_trend(
      pre_data,
      y = y,
      ivar = ivar,
      tvar = tvar,
      controls = controls
    )
    trend_estimate$cohort <- as.integer(cohort_value)
    trend_by_cohort[[length(trend_by_cohort) + 1L]] <- trend_estimate
  }

  control_group_trend <- NULL
  if (isTRUE(include_control_group) && !is.null(gvar)) {
    control_mask <- is_never_treated(dt[[gvar]])
    control_data <- dt[control_mask]
    if (nrow(control_data) >= 3L) {
      control_group_trend <- .estimate_cohort_trend(
        control_data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        controls = controls
      )
      control_group_trend$cohort <- 0L
    }
  }

  trend_heterogeneity_test <- if (!is.null(gvar) && length(trend_by_cohort) >= 2L) {
    .test_trend_heterogeneity(
      dt,
      y = y,
      ivar = ivar,
      tvar = tvar,
      gvar = gvar,
      controls = controls,
      never_treated_values = never_treated_values,
      alpha = alpha
    )
  } else {
    .default_trend_heterogeneity_test()
  }

  trend_differences <- if (length(trend_by_cohort) >= 2L) {
    .compute_pairwise_trend_differences(
      trend_by_cohort = trend_by_cohort,
      control_group_trend = control_group_trend,
      alpha = alpha
    )
  } else {
    list()
  }

  has_heterogeneous_trends <- isTRUE(trend_heterogeneity_test$reject_null)
  if (has_heterogeneous_trends) {
    recommendation <- "detrend"
    recommendation_confidence <- min(0.95, 1 - trend_heterogeneity_test$pvalue)
    recommendation_reason <- sprintf(
      paste(
        "Significant trend heterogeneity detected (F-test p=%.4f).",
        "Detrending removes cohort-specific linear trends under Assumption CHT."
      ),
      trend_heterogeneity_test$pvalue
    )
  } else {
    recommendation <- "demean"
    recommendation_confidence <- trend_heterogeneity_test$pvalue
    recommendation_reason <- sprintf(
      paste(
        "No significant trend heterogeneity detected (F-test p=%.4f).",
        "Demeaning is more efficient when trends are parallel."
      ),
      trend_heterogeneity_test$pvalue
    )
  }

  result <- structure(
    list(
      trend_by_cohort = trend_by_cohort,
      trend_heterogeneity_test = trend_heterogeneity_test,
      trend_differences = trend_differences,
      control_group_trend = control_group_trend,
      has_heterogeneous_trends = has_heterogeneous_trends,
      recommendation = recommendation,
      recommendation_confidence = recommendation_confidence,
      recommendation_reason = recommendation_reason,
      figure = NULL,
      cohorts = as.integer(cohorts),
      gvar = gvar,
      alpha = alpha,
      include_control_group = isTRUE(include_control_group)
    ),
    class = c("lwdid_heterogeneous_trends", "list")
  )

  if (isTRUE(verbose)) {
    message(result$recommendation_reason)
  }

  result
}


#' Resolve a transformation recommendation from normalized evidence scores
#'
#' @param score_demean Numeric score for demeaning.
#' @param score_detrend Numeric score for detrending.
#' @param detrend_feasible Logical scalar.
#' @param has_seasonal Logical scalar.
#' @param warnings_list Character vector of accumulated warnings.
#' @param reasons Character vector of accumulated reasons.
#' @return Named list with recommendation fields.
#' @keywords internal
.resolve_transformation_recommendation <- function(
    score_demean,
    score_detrend,
    detrend_feasible,
    has_seasonal,
    warnings_list = character(0),
    reasons = character(0)
) {
  total_score <- score_demean + score_detrend
  if (total_score > 0) {
    score_demean <- score_demean / total_score
    score_detrend <- score_detrend / total_score
  }

  if (!detrend_feasible) {
    recommended_method <- "demean"
    confidence <- 1.0
    reasons <- c(
      "Detrending not feasible (insufficient pre-treatment periods).",
      reasons
    )
  } else if (has_seasonal) {
    if (score_detrend > score_demean) {
      recommended_method <- "detrendq"
      confidence <- score_detrend
    } else {
      recommended_method <- "demeanq"
      confidence <- score_demean
    }
    reasons <- c(reasons, "Seasonal patterns detected - using seasonal variant.")
  } else if (score_detrend > score_demean) {
    recommended_method <- "detrend"
    confidence <- score_detrend
  } else {
    recommended_method <- "demean"
    confidence <- score_demean
  }

  confidence_level <- if (confidence > 0.8) {
    "High"
  } else if (confidence > 0.5) {
    "Medium"
  } else {
    warnings_list <- c(
      warnings_list,
      paste(
        "Low confidence in recommendation.",
        "Consider running sensitivity analysis comparing demean and detrend."
      )
    )
    "Low"
  }

  if (recommended_method %in% c("demean", "demeanq") && detrend_feasible) {
    alternative_method <- if (has_seasonal) "detrendq" else "detrend"
    alternative_reason <- "Use if parallel trends assumption is questionable."
  } else if (recommended_method %in% c("detrend", "detrendq")) {
    alternative_method <- if (has_seasonal) "demeanq" else "demean"
    alternative_reason <- paste(
      "Use if parallel trends assumption is believed to hold",
      "(more efficient)."
    )
  } else {
    alternative_method <- NULL
    alternative_reason <- NULL
  }

  if (has_seasonal && detrend_feasible) {
    scores <- c(demeanq = score_demean, detrendq = score_detrend)
  } else {
    scores <- c(demean = score_demean, detrend = score_detrend)
  }

  list(
    recommended_method = recommended_method,
    confidence = as.numeric(confidence),
    confidence_level = confidence_level,
    reasons = reasons,
    alternative_method = alternative_method,
    alternative_reason = alternative_reason,
    warnings = warnings_list,
    scores = scores
  )
}


#' Recommend a transformation for trend-robust DiD estimation
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param controls Optional character vector of controls.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param run_all_diagnostics Logical scalar. Whether to run PT and HT diagnostics.
#' @param verbose Logical scalar.
#' @return `lwdid_transformation_recommendation` object.
#' @export
lwdid_recommend_transformation <- function(
    data,
    y,
    ivar,
    tvar,
    gvar = NULL,
    controls = NULL,
    never_treated_values = c(0, Inf),
    run_all_diagnostics = TRUE,
    verbose = TRUE
) {
  required_columns <- c(y, ivar, tvar, controls %||% character(0))
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(lwdid_missing_column_error(missing_columns[[1L]], names(data)))
  }
  never_treated_values <- .validate_never_treated_values(never_treated_values)

  warnings_list <- character(0)
  reasons <- character(0)

  pre_period_range <- .compute_pre_period_range(
    data,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    never_treated_values = never_treated_values
  )
  n_pre_min <- as.integer(pre_period_range[[1L]])
  n_pre_max <- as.integer(pre_period_range[[2L]])

  if (n_pre_min < 1L) {
    stop("At least one pre-treatment period required for any transformation.")
  }

  detrend_feasible <- n_pre_min >= 2L
  if (!detrend_feasible) {
    warnings_list <- c(
      warnings_list,
      sprintf(
        paste(
          "Detrending requires >=2 pre-treatment periods.",
          "Minimum found: %d. Only demeaning is feasible."
        ),
        n_pre_min
      )
    )
  }

  is_balanced <- .check_panel_balance(data, ivar = ivar, tvar = tvar)
  if (!is_balanced) {
    reasons <- c(
      reasons,
      "Unbalanced panel detected. Detrending is more robust to selection."
    )
  }

  has_seasonal <- .detect_seasonal_patterns(
    data,
    y = y,
    ivar = ivar,
    tvar = tvar
  )

  pt_test_result <- NULL
  ht_diag_result <- NULL

  if (isTRUE(run_all_diagnostics) && !is.null(gvar)) {
    pt_test_result <- tryCatch(
      lwdid_test_parallel_trends(
        data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        controls = controls,
        method = "placebo",
        alpha = 0.05,
        verbose = FALSE,
        never_treated_values = never_treated_values
      ),
      error = function(e) {
        warnings_list <<- c(
          warnings_list,
          sprintf("Parallel trends test failed: %s", conditionMessage(e))
        )
        NULL
      }
    )

    ht_diag_result <- tryCatch(
      lwdid_diagnose_heterogeneous_trends(
        data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        controls = controls,
        never_treated_values = never_treated_values,
        include_control_group = TRUE,
        alpha = 0.05,
        verbose = FALSE
      ),
      error = function(e) {
        warnings_list <<- c(
          warnings_list,
          sprintf(
            "Heterogeneous trends diagnosis failed: %s",
            conditionMessage(e)
          )
        )
        NULL
      }
    )
  }

  score_demean <- 0.5
  score_detrend <- 0.5

  if (!is.null(pt_test_result)) {
    if (isTRUE(pt_test_result$reject_null)) {
      score_detrend <- score_detrend + 0.4
      reasons <- c(
        reasons,
        sprintf(
          "Parallel trends test rejected (p=%.4f).",
          pt_test_result$joint_pvalue %||% pt_test_result$pvalue
        )
      )
    } else {
      score_demean <- score_demean + 0.4
      reasons <- c(
        reasons,
        sprintf(
          "Parallel trends test not rejected (p=%.4f).",
          pt_test_result$joint_pvalue %||% pt_test_result$pvalue
        )
      )
    }
  }

  if (!is.null(ht_diag_result)) {
    ht_pvalue <- ht_diag_result$trend_heterogeneity_test$pvalue %||% 1.0
    if (isTRUE(ht_diag_result$has_heterogeneous_trends)) {
      score_detrend <- score_detrend + 0.4
      reasons <- c(
        reasons,
        sprintf("Heterogeneous trends detected (F-test p=%.4f).", ht_pvalue)
      )
    } else {
      score_demean <- score_demean + 0.4
      reasons <- c(
        reasons,
        sprintf("No heterogeneous trends detected (F-test p=%.4f).", ht_pvalue)
      )
    }
  }

  if (!is_balanced) {
    score_detrend <- score_detrend + 0.1
  } else {
    score_demean <- score_demean + 0.1
  }

  if (detrend_feasible) {
    score_demean <- score_demean + 0.05
    score_detrend <- score_detrend + 0.05
  } else {
    score_demean <- score_demean + 0.1
    score_detrend <- 0
  }

  recommendation_decision <- .resolve_transformation_recommendation(
    score_demean = score_demean,
    score_detrend = score_detrend,
    detrend_feasible = detrend_feasible,
    has_seasonal = has_seasonal,
    warnings_list = warnings_list,
    reasons = reasons
  )

  result <- structure(
    list(
      recommended_method = recommendation_decision$recommended_method,
      confidence = recommendation_decision$confidence,
      confidence_level = recommendation_decision$confidence_level,
      reasons = recommendation_decision$reasons,
      parallel_trends_test = pt_test_result,
      heterogeneous_trends_diag = ht_diag_result,
      n_pre_periods_min = n_pre_min,
      n_pre_periods_max = n_pre_max,
      has_seasonal_pattern = has_seasonal,
      is_balanced_panel = is_balanced,
      alternative_method = recommendation_decision$alternative_method,
      alternative_reason = recommendation_decision$alternative_reason,
      warnings = recommendation_decision$warnings,
      scores = recommendation_decision$scores
    ),
    class = c("lwdid_transformation_recommendation", "list")
  )

  if (isTRUE(verbose)) {
    message(
      sprintf(
        "Recommended %s (confidence: %.2f, level: %s).",
        result$recommended_method,
        result$confidence,
        result$confidence_level
      )
    )
  }

  result
}

#' @keywords internal
.format_trend_numeric <- function(value, digits = 4L) {
  if (is.null(value) || length(value) == 0L) {
    return("NA")
  }

  numeric_value <- suppressWarnings(as.numeric(value[[1L]]))
  if (!is.finite(numeric_value)) {
    return("NA")
  }

  sprintf("%.*f", digits, numeric_value)
}

#' @keywords internal
.emit_parallel_trends_summary <- function(x, detailed = FALSE) {
  joint_df <- as.integer(x$joint_df %||% c(0L, 0L))
  joint_label <- if (length(joint_df) == 2L && all(joint_df >= 0L)) {
    sprintf(
      "F(%d, %d) = %s, p = %s",
      joint_df[[1L]],
      joint_df[[2L]],
      .format_trend_numeric(x$joint_f_stat),
      .format_pvalue(x$joint_pvalue, 4L)
    )
  } else {
    "Joint F-test unavailable"
  }

  cat("Parallel Trends Diagnostics\n")
  cat(sprintf("  Recommendation: %s\n", x$recommendation %||% "NA"))
  cat(sprintf("  Joint test: %s\n", joint_label))
  cat(sprintf(
    "  Pre-treatment estimates: %d\n",
    length(x$pre_trend_estimates %||% list())
  ))

  if (isTRUE(detailed) && length(x$pre_trend_estimates %||% list()) > 0L) {
    event_times <- vapply(
      x$pre_trend_estimates,
      function(estimate) as.integer(estimate$event_time %||% NA_integer_),
      integer(1)
    )
    cat(sprintf(
      "  Event times: %s\n",
      paste(event_times, collapse = ", ")
    ))
  }

  if (length(x$warnings %||% character(0)) > 0L) {
    cat("  Warnings:\n")
    for (warning_message in x$warnings) {
      cat(sprintf("    - %s\n", warning_message))
    }
  }

  invisible(x)
}

#' @keywords internal
.emit_heterogeneous_trends_summary <- function(x, detailed = FALSE) {
  test_result <- x$trend_heterogeneity_test %||% list()

  cat("Heterogeneous Trend Diagnostics\n")
  cat(sprintf("  Recommendation: %s\n", x$recommendation %||% "NA"))
  cat(sprintf(
    "  Heterogeneity test: F(%d, %d) = %s, p = %s\n",
    as.integer(test_result$df_num %||% 0L),
    as.integer(test_result$df_den %||% 0L),
    .format_trend_numeric(test_result$f_stat),
    .format_pvalue(test_result$pvalue, 4L)
  ))
  cat(sprintf(
    "  Cohort trends estimated: %d\n",
    length(x$trend_by_cohort %||% list())
  ))
  cat(sprintf(
    "  Pairwise differences: %d\n",
    length(x$trend_differences %||% list())
  ))

  if (isTRUE(detailed) && length(x$trend_by_cohort %||% list()) > 0L) {
    slopes <- vapply(
      x$trend_by_cohort,
      function(trend) {
        sprintf(
          "%s=%s",
          as.integer(trend$cohort %||% NA_integer_),
          .format_trend_numeric(trend$slope)
        )
      },
      character(1)
    )
    cat(sprintf("  Cohort slopes: %s\n", paste(slopes, collapse = "; ")))
  }

  invisible(x)
}

#' @keywords internal
.emit_transformation_recommendation_summary <- function(x, detailed = FALSE) {
  cat("Transformation Recommendation\n")
  cat(sprintf(
    "  Recommended method: %s\n",
    x$recommended_method %||% "NA"
  ))
  cat(sprintf(
    "  Confidence: %s (%s)\n",
    .format_trend_numeric(x$confidence),
    x$confidence_level %||% "NA"
  ))
  cat(sprintf(
    "  Pre-treatment periods: min=%d, max=%d\n",
    as.integer(x$n_pre_periods_min %||% 0L),
    as.integer(x$n_pre_periods_max %||% 0L)
  ))

  if (isTRUE(detailed) && !is.null(x$scores)) {
    score_labels <- paste(
      sprintf("%s=%s", names(x$scores), vapply(x$scores, .format_trend_numeric, character(1))),
      collapse = ", "
    )
    cat(sprintf("  Scores: %s\n", score_labels))
  }

  if (!is.null(x$alternative_method)) {
    cat(sprintf("  Alternative: %s\n", x$alternative_method))
  }

  if (length(x$warnings %||% character(0)) > 0L) {
    cat("  Warnings:\n")
    for (warning_message in x$warnings) {
      cat(sprintf("    - %s\n", warning_message))
    }
  }

  invisible(x)
}

#' Print method for parallel trends diagnostics
#'
#' @param x An object of class `lwdid_parallel_trends`.
#' @param ... Additional arguments (ignored).
#' @return `x` invisibly.
#' @export
print.lwdid_parallel_trends <- function(x, ...) {
  .emit_parallel_trends_summary(x, detailed = FALSE)
}

#' Summary method for parallel trends diagnostics
#'
#' @param object An object of class `lwdid_parallel_trends`.
#' @param ... Additional arguments (ignored).
#' @return `object` invisibly.
#' @export
summary.lwdid_parallel_trends <- function(object, ...) {
  .emit_parallel_trends_summary(object, detailed = TRUE)
}

#' Print method for heterogeneous trend diagnostics
#'
#' @param x An object of class `lwdid_heterogeneous_trends`.
#' @param ... Additional arguments (ignored).
#' @return `x` invisibly.
#' @export
print.lwdid_heterogeneous_trends <- function(x, ...) {
  .emit_heterogeneous_trends_summary(x, detailed = FALSE)
}

#' Summary method for heterogeneous trend diagnostics
#'
#' @param object An object of class `lwdid_heterogeneous_trends`.
#' @param ... Additional arguments (ignored).
#' @return `object` invisibly.
#' @export
summary.lwdid_heterogeneous_trends <- function(object, ...) {
  .emit_heterogeneous_trends_summary(object, detailed = TRUE)
}

#' Print method for transformation recommendations
#'
#' @param x An object of class `lwdid_transformation_recommendation`.
#' @param ... Additional arguments (ignored).
#' @return `x` invisibly.
#' @export
print.lwdid_transformation_recommendation <- function(x, ...) {
  .emit_transformation_recommendation_summary(x, detailed = FALSE)
}

#' Summary method for transformation recommendations
#'
#' @param object An object of class `lwdid_transformation_recommendation`.
#' @param ... Additional arguments (ignored).
#' @return `object` invisibly.
#' @export
summary.lwdid_transformation_recommendation <- function(object, ...) {
  .emit_transformation_recommendation_summary(object, detailed = TRUE)
}
