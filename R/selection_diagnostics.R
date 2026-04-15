# ============================================================================
# selection_diagnostics.R — Internal selection diagnostics primitives
#
# Foundations for selection diagnostics in Lee-Wooldridge DiD applications.
# This initial slice fixes the observed-Y availability contract required by
# story E8-05 for missing-rate helpers and unit-level usability summaries.
# ============================================================================

#' @title Selection Diagnostics Internals
#' @description Internal helpers for selection diagnostics.
#' @name selection_diagnostics
#' @keywords internal
NULL


#' Selection-specific never-treated helper
#'
#' @param g Numeric scalar or vector of cohort values.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return Logical vector indicating never-treated status.
#' @keywords internal
.is_never_treated <- function(g, never_treated_values = c(0, Inf)) {
  if (length(g) == 0L) {
    return(logical(0))
  }

  is_missing <- is.na(g)
  is_infinite_value <- !is_missing & is.infinite(g)
  is_explicit_marker <- !is_missing & (g %in% never_treated_values)

  is_missing | is_infinite_value | is_explicit_marker
}


#' Validate core selection inputs
#'
#' @param data data.frame or data.table.
#' @param columns Character vector of required column names.
#' @return NULL (invisible). Raises errors on failure.
#' @keywords internal
.validate_selection_inputs <- function(data, columns) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }
  .validate_data_not_empty(data)

  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Column '%s' not found in data", missing_columns[[1L]]),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Build the full unit-period grid used by observed-Y diagnostics
#'
#' @param dt data.table panel.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @return data.table full grid with merged outcome column.
#' @keywords internal
.build_selection_grid <- function(dt, ivar, tvar, y) {
  units <- sort(unique(dt[[ivar]]))
  periods <- sort(unique(dt[[tvar]]))

  full_grid <- data.table::CJ(unit = units, period = periods)
  data.table::setnames(full_grid, c(ivar, tvar))

  merge(
    full_grid,
    dt[, c(ivar, tvar, y), with = FALSE],
    by = c(ivar, tvar),
    all.x = TRUE,
    sort = FALSE
  )
}


#' Compute missing-rate summaries from observed outcomes
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @return Named list of missing-rate diagnostics.
#' @keywords internal
.compute_missing_rates <- function(data, ivar, tvar, y) {
  .validate_selection_inputs(data, c(ivar, tvar, y))

  dt <- data.table::as.data.table(data)
  merged <- .build_selection_grid(dt, ivar, tvar, y)

  n_expected <- as.integer(nrow(merged))
  n_observed <- as.integer(sum(!is.na(merged[[y]])))

  list(
    n_units = as.integer(data.table::uniqueN(merged[[ivar]])),
    n_periods = as.integer(data.table::uniqueN(merged[[tvar]])),
    n_expected = n_expected,
    n_observed = n_observed,
    overall_missing_rate = if (n_expected > 0L) {
      as.numeric(mean(is.na(merged[[y]])))
    } else {
      0
    },
    is_balanced = identical(n_observed, n_expected)
  )
}


#' Compute missing rates by period from observed outcomes
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @return Named list keyed by period.
#' @keywords internal
.compute_missing_by_period <- function(data, ivar, tvar, y) {
  .validate_selection_inputs(data, c(ivar, tvar, y))

  dt <- data.table::as.data.table(data)
  merged <- .build_selection_grid(dt, ivar, tvar, y)

  period_rates <- merged[
    ,
    .(missing_rate = mean(is.na(.SD[[1L]]))),
    by = c(tvar),
    .SDcols = y
  ]

  stats::setNames(
    as.list(as.numeric(period_rates$missing_rate)),
    as.character(period_rates[[tvar]])
  )
}


#' Compute missing rates by treated cohort from observed outcomes
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param gvar Character cohort column.
#' @param y Character outcome column.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return Named list keyed by cohort.
#' @keywords internal
.compute_missing_by_cohort <- function(
    data, ivar, tvar, gvar, y, never_treated_values = c(0, Inf)
) {
  .validate_selection_inputs(data, c(ivar, tvar, gvar, y))

  dt <- data.table::as.data.table(data)
  merged <- .build_selection_grid(dt, ivar, tvar, y)
  unit_gvar <- get_unit_level_gvar(dt, gvar, ivar)

  treated_cohorts <- sort(unique(unit_gvar[[gvar]][
    !.is_never_treated(unit_gvar[[gvar]], never_treated_values)
  ]))

  if (length(treated_cohorts) == 0L) {
    return(list())
  }

  cohort_rates <- lapply(treated_cohorts, function(cohort) {
    cohort_mask <- unit_gvar[[gvar]] == cohort
    cohort_units <- unit_gvar[cohort_mask][[ivar]]
    cohort_data <- merged[merged[[ivar]] %in% cohort_units]
    as.numeric(mean(is.na(cohort_data[[y]])))
  })

  stats::setNames(cohort_rates, as.character(treated_cohorts))
}


#' Classify missing-data mechanism from observed-Y availability
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @param gvar Unused cohort column placeholder for interface consistency.
#' @param d Unused common-timing placeholder for interface consistency.
#' @param covariates Character vector of optional MAR controls.
#' @return Named list with pattern, confidence, tests, and description.
#' @keywords internal
.classify_missing_pattern <- function(
    data,
    ivar,
    tvar,
    y,
    gvar = NULL,
    d = NULL,
    covariates = NULL
) {
  required_columns <- c(ivar, tvar, y)
  available_covariates <- character(0)

  if (!is.null(covariates)) {
    available_covariates <- covariates[covariates %in% names(data)]
    required_columns <- unique(c(required_columns, available_covariates))
  }

  .validate_selection_inputs(data, required_columns)

  dt <- data.table::as.data.table(data)
  merged <- .build_selection_grid(dt, ivar, tvar, y)
  merged[, `_missing` := as.integer(is.na(.SD[[1L]])), .SDcols = y]

  if (sum(merged[["_missing"]]) == 0L) {
    return(list(
      pattern = "MCAR",
      confidence = 1,
      tests = list(),
      description = "No missing outcomes detected; panel is classified as MCAR."
    ))
  }

  if (sum(!is.na(merged[[y]])) == 0L) {
    return(list(
      pattern = "UNKNOWN",
      confidence = 0.3,
      tests = list(),
      description = paste(
        "All outcomes are missing;",
        "the selection mechanism cannot be classified."
      )
    ))
  }

  tests <- list()

  obs_per_unit <- merged[
    ,
    .(n_missing = sum(.SD[[1L]])),
    by = c(ivar),
    .SDcols = "_missing"
  ]
  complete_units <- obs_per_unit[n_missing == 0L][[ivar]]
  incomplete_units <- obs_per_unit[n_missing > 0L][[ivar]]

  if (length(complete_units) >= 5L && length(incomplete_units) >= 5L) {
    observed_y_mask <- !is.na(dt[[y]])
    complete_y <- dt[[y]][dt[[ivar]] %in% complete_units & observed_y_mask]
    incomplete_y <- dt[[y]][dt[[ivar]] %in% incomplete_units & observed_y_mask]

    if (length(complete_y) > 1L && length(incomplete_y) > 1L) {
      t_result <- tryCatch(
        stats::t.test(complete_y, incomplete_y),
        error = function(e) NULL
      )

      if (!is.null(t_result)) {
        pvalue <- as.numeric(t_result$p.value)
        reject_null <- isTRUE(is.finite(pvalue) && pvalue < 0.05)
        tests[[length(tests) + 1L]] <- list(
          test_name = "Simplified Little's MCAR Test",
          statistic = as.numeric(unname(t_result$statistic)),
          pvalue = pvalue,
          reject_null = reject_null,
          interpretation = if (reject_null) {
            paste(
              "Reject MCAR: significant difference in outcomes between",
              "complete and incomplete units"
            )
          } else {
            "Cannot reject MCAR: no significant difference detected"
          }
        )
      }
    }
  }

  if (length(available_covariates) > 0L) {
    unit_controls <- dt[!duplicated(dt[[ivar]]), c(ivar, available_covariates), with = FALSE]
    if (nrow(unit_controls) > 0L) {
      complete_rows <- stats::complete.cases(unit_controls[, available_covariates, with = FALSE])
      unit_controls <- unit_controls[complete_rows]
    }

    unit_missing <- merged[
      ,
      .(unit_missing_rate = mean(.SD[[1L]])),
      by = c(ivar),
      .SDcols = "_missing"
    ]
    test_data <- merge(unit_missing, unit_controls, by = ivar)

    if (nrow(test_data) > length(available_covariates) + 2L) {
      formula <- stats::as.formula(
        paste("unit_missing_rate ~", paste(available_covariates, collapse = " + "))
      )
      fit <- tryCatch(
        stats::lm(formula, data = test_data),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        f_stat <- summary(fit)$fstatistic

        if (!is.null(f_stat) && all(is.finite(f_stat))) {
          pvalue <- as.numeric(stats::pf(f_stat[[1L]], f_stat[[2L]], f_stat[[3L]], lower.tail = FALSE))
          reject_null <- isTRUE(is.finite(pvalue) && pvalue < 0.05)
          tests[[length(tests) + 1L]] <- list(
            test_name = "Selection on Observables (MAR) Test",
            statistic = as.numeric(unname(f_stat[[1L]])),
            pvalue = pvalue,
            reject_null = reject_null,
            interpretation = if (reject_null) {
              "Selection depends on observed controls (MAR)"
            } else {
              "No evidence of selection on observed controls"
            }
          )
        }
      }
    }
  }

  data_sorted <- data.table::copy(dt)
  data.table::setorderv(data_sorted, c(ivar, tvar))
  data_sorted[
    ,
    `_y_lag` := data.table::shift(.SD[[1L]], 1L, type = "lag"),
    by = ivar,
    .SDcols = y
  ]

  lag_test <- merge(
    merged[, c(ivar, tvar, "_missing"), with = FALSE],
    data_sorted[, c(ivar, tvar, "_y_lag"), with = FALSE],
    by = c(ivar, tvar),
    all.x = TRUE,
    sort = FALSE
  )
  lag_test <- lag_test[!is.na(`_y_lag`)]

  if (nrow(lag_test) > 10L) {
    cor_result <- tryCatch(
      stats::cor.test(lag_test[["_missing"]], lag_test[["_y_lag"]], method = "pearson"),
      error = function(e) NULL
    )

    if (!is.null(cor_result)) {
      pvalue <- as.numeric(cor_result$p.value)
      reject_null <- isTRUE(is.finite(pvalue) && pvalue < 0.05)
      tests[[length(tests) + 1L]] <- list(
        test_name = "Selection on Lagged Outcome (MNAR) Test",
        statistic = as.numeric(unname(cor_result$estimate)),
        pvalue = pvalue,
        reject_null = reject_null,
        interpretation = if (reject_null) {
          paste(
            "WARNING: Missingness correlates with lagged outcomes.",
            "This suggests potential MNAR and may violate the selection mechanism assumption."
          )
        } else {
          "No evidence of selection on lagged outcomes"
        }
      )
    }
  }

  mcar_rejected <- any(vapply(
    tests,
    function(test) identical(test$test_name, "Simplified Little's MCAR Test") &&
      isTRUE(test$reject_null),
    logical(1)
  ))
  mar_detected <- any(vapply(
    tests,
    function(test) identical(test$test_name, "Selection on Observables (MAR) Test") &&
      isTRUE(test$reject_null),
    logical(1)
  ))
  mnar_detected <- any(vapply(
    tests,
    function(test) identical(test$test_name, "Selection on Lagged Outcome (MNAR) Test") &&
      isTRUE(test$reject_null),
    logical(1)
  ))

  if (mnar_detected) {
    pattern <- "MNAR"
    confidence <- 0.7
    description <- "Missingness correlates with lagged outcomes and is classified as MNAR."
  } else if (mar_detected) {
    pattern <- "MAR"
    confidence <- 0.8
    description <- "Missingness appears related to observed covariates and is classified as MAR."
  } else if (!mcar_rejected) {
    pattern <- "MCAR"
    confidence <- if (length(tests) > 0L) 0.9 else 0.5
    description <- "No diagnostic test rejected MCAR."
  } else {
    pattern <- "UNKNOWN"
    confidence <- 0.3
    description <- "Available tests do not support a stable missingness classification."
  }

  list(
    pattern = pattern,
    confidence = confidence,
    tests = tests,
    description = description
  )
}


#' Compute observed-Y unit statistics
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param gvar Character cohort column or NULL.
#' @param y Character outcome column.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return data.table with one row per unit.
#' @keywords internal
.compute_unit_stats <- function(
    data, ivar, tvar, gvar, y, never_treated_values = c(0, Inf), d = NULL
) {
  required_columns <- c(ivar, tvar, y)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  if (!is.null(d)) {
    required_columns <- c(required_columns, d)
  }
  .validate_selection_inputs(data, required_columns)

  dt <- data.table::as.data.table(data)
  merged <- .build_selection_grid(dt, ivar, tvar, y)
  periods <- sort(unique(merged[[tvar]]))
  n_total_periods <- as.integer(length(periods))

  if (!is.null(gvar)) {
    unit_gvar <- get_unit_level_gvar(dt, gvar, ivar)
  } else {
    unit_gvar <- NULL
  }
  if (!is.null(d)) {
    unit_treatment <- dt[
      ,
      .(
        .is_treated = any(!is.na(.SD[[1L]]) & .SD[[1L]] == 1L),
        .treatment_time = if (any(!is.na(.SD[[1L]]) & .SD[[1L]] == 1L)) {
          as.numeric(min(.SD[[2L]][!is.na(.SD[[1L]]) & .SD[[1L]] == 1L]))
        } else {
          NA_real_
        }
      ),
      by = ivar,
      .SDcols = c(d, tvar)
    ]
  } else {
    unit_treatment <- NULL
  }

  unit_rows <- lapply(sort(unique(merged[[ivar]])), function(unit_id) {
    current_unit_id <- unit_id
    unit_mask <- merged[[ivar]] == current_unit_id
    unit_grid <- merged[unit_mask]
    observed_mask <- !is.na(unit_grid[[y]])
    observed_grid <- unit_grid[observed_mask]

    n_observed <- as.integer(nrow(observed_grid))
    n_missing <- as.integer(n_total_periods - n_observed)

    if (n_observed > 0L) {
      first_observed <- as.integer(min(observed_grid[[tvar]]))
      last_observed <- as.integer(max(observed_grid[[tvar]]))
      observation_span <- as.integer(last_observed - first_observed + 1L)
    } else {
      first_observed <- NA_integer_
      last_observed <- NA_integer_
      observation_span <- NA_integer_
    }

    row <- list(
      unit_id = current_unit_id,
      cohort = NA_integer_,
      is_treated = FALSE,
      n_total_periods = n_total_periods,
      n_observed = n_observed,
      n_missing = n_missing,
      missing_rate = if (n_total_periods > 0L) n_missing / n_total_periods else 0,
      first_observed = first_observed,
      last_observed = last_observed,
      observation_span = observation_span,
      n_pre_treatment = NA_integer_,
      n_post_treatment = NA_integer_,
      pre_treatment_missing_rate = NA_real_,
      post_treatment_missing_rate = NA_real_,
      can_use_demean = TRUE,
      can_use_detrend = TRUE,
      reason_if_excluded = NA_character_
    )

    if (!is.null(unit_gvar)) {
      g_value <- unit_gvar[[gvar]][match(current_unit_id, unit_gvar[[ivar]])]

      if (!.is_never_treated(g_value, never_treated_values)) {
        pre_periods <- periods[periods < g_value]
        post_periods <- periods[periods >= g_value]
        observed_periods <- observed_grid[[tvar]]
        n_pre_treatment <- as.integer(sum(observed_periods < g_value))
        n_post_treatment <- as.integer(sum(observed_periods >= g_value))

        row$cohort <- as.integer(g_value)
        row$is_treated <- TRUE
        row$n_pre_treatment <- n_pre_treatment
        row$n_post_treatment <- n_post_treatment
        row$pre_treatment_missing_rate <- if (length(pre_periods) > 0L) {
          1 - n_pre_treatment / length(pre_periods)
        } else {
          NA_real_
        }
        row$post_treatment_missing_rate <- if (length(post_periods) > 0L) {
          1 - n_post_treatment / length(post_periods)
        } else {
          NA_real_
        }

        if (n_pre_treatment < 1L) {
          row$can_use_demean <- FALSE
          row$reason_if_excluded <- "No pre-treatment observations"
        }
        if (n_pre_treatment < 2L) {
          row$can_use_detrend <- FALSE
          if (is.na(row$reason_if_excluded)) {
            row$reason_if_excluded <- "Fewer than 2 pre-treatment observations"
          }
        }
      }
    } else if (!is.null(unit_treatment)) {
      treatment_row <- unit_treatment[match(current_unit_id, unit_treatment[[ivar]])]
      is_treated <- isTRUE(treatment_row$.is_treated[[1L]])
      treatment_time <- treatment_row$.treatment_time[[1L]]

      if (is_treated && is.finite(treatment_time)) {
        pre_periods <- periods[periods < treatment_time]
        post_periods <- periods[periods >= treatment_time]
        observed_periods <- observed_grid[[tvar]]
        n_pre_treatment <- as.integer(sum(observed_periods < treatment_time))
        n_post_treatment <- as.integer(sum(observed_periods >= treatment_time))

        row$cohort <- as.integer(treatment_time)
        row$is_treated <- TRUE
        row$n_pre_treatment <- n_pre_treatment
        row$n_post_treatment <- n_post_treatment
        row$pre_treatment_missing_rate <- if (length(pre_periods) > 0L) {
          1 - n_pre_treatment / length(pre_periods)
        } else {
          NA_real_
        }
        row$post_treatment_missing_rate <- if (length(post_periods) > 0L) {
          1 - n_post_treatment / length(post_periods)
        } else {
          NA_real_
        }

        if (n_pre_treatment < 1L) {
          row$can_use_demean <- FALSE
          row$reason_if_excluded <- "No pre-treatment observations"
        }
        if (n_pre_treatment < 2L) {
          row$can_use_detrend <- FALSE
          if (is.na(row$reason_if_excluded)) {
            row$reason_if_excluded <- "Fewer than 2 pre-treatment observations"
          }
        }
      }
    }

    row
  })

  data.table::rbindlist(unit_rows)
}


#' Compute observed-Y attrition diagnostics
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @param gvar Character cohort column or NULL.
#' @param d Unused common-timing convenience parameter.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return Named list of attrition diagnostics.
#' @keywords internal
.compute_attrition_analysis <- function(
    data,
    ivar,
    tvar,
    y,
    gvar = NULL,
    d = NULL,
    never_treated_values = c(0, Inf)
) {
  required_columns <- c(ivar, tvar, y)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  .validate_selection_inputs(data, required_columns)

  dt <- data.table::as.data.table(data)
  unit_stats <- .compute_unit_stats(
    dt,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    y = y,
    never_treated_values = never_treated_values,
    d = d
  )
  all_periods <- sort(unique(dt[[tvar]]))
  t_min <- min(all_periods)
  t_max <- max(all_periods)
  n_units <- nrow(unit_stats)
  n_complete <- as.integer(sum(unit_stats$n_missing == 0L))
  n_partial <- as.integer(n_units - n_complete)

  treated_stats <- unit_stats[unit_stats$is_treated, , drop = FALSE]
  control_stats <- unit_stats[!unit_stats$is_treated, , drop = FALSE]

  dropout_before_treatment <- 0L
  dropout_after_treatment <- 0L

  if (nrow(treated_stats) > 0L) {
    observed_last <- !is.na(treated_stats$last_observed)
    dropout_before_treatment <- as.integer(sum(
      observed_last & treated_stats$last_observed < treated_stats$cohort
    ))
    dropout_after_treatment <- as.integer(sum(
      observed_last &
        treated_stats$last_observed >= treated_stats$cohort &
        treated_stats$last_observed < t_max
    ))
  }

  treated_attrition <- if (nrow(treated_stats) > 0L) {
    mean(treated_stats$n_missing > 0L)
  } else {
    0
  }
  control_attrition <- if (nrow(control_stats) > 0L) {
    mean(control_stats$n_missing > 0L)
  } else {
    0
  }

  attrition_by_period <- .compute_missing_by_period(
    dt,
    ivar = ivar,
    tvar = tvar,
    y = y
  )

  attrition_by_cohort <- list()
  if (!is.null(gvar) && nrow(treated_stats) > 0L) {
    cohort_rates <- tapply(
      treated_stats$n_missing > 0L,
      treated_stats$cohort,
      mean
    )
    attrition_by_cohort <- as.list(as.numeric(cohort_rates))
    names(attrition_by_cohort) <- names(cohort_rates)
  }

  list(
    n_complete = n_complete,
    n_partial = n_partial,
    overall_attrition = if (n_units > 0L) n_partial / n_units else 0,
    attrition_by_cohort = attrition_by_cohort,
    attrition_by_period = attrition_by_period,
    early_dropout_rate = if (n_units > 0L) {
      mean(!is.na(unit_stats$last_observed) & unit_stats$last_observed < t_max)
    } else {
      0
    },
    late_entry_rate = if (n_units > 0L) {
      mean(!is.na(unit_stats$first_observed) & unit_stats$first_observed > t_min)
    } else {
      0
    },
    dropout_before_treatment = dropout_before_treatment,
    dropout_after_treatment = dropout_after_treatment,
    treatment_related_attrition = if (
      nrow(treated_stats) > 0L &&
      nrow(control_stats) > 0L
    ) {
      as.numeric(treated_attrition - control_attrition)
    } else {
      0
    }
  )
}


#' Compute observed-Y balance diagnostics
#'
#' @param data data.frame or data.table.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @param gvar Character cohort column or NULL.
#' @param d Unused common-timing convenience parameter.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @return Named list of balance diagnostics.
#' @keywords internal
.compute_balance_statistics <- function(
    data,
    ivar,
    tvar,
    y,
    gvar = NULL,
    d = NULL,
    never_treated_values = c(0, Inf)
) {
  required_columns <- c(ivar, tvar, y)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  .validate_selection_inputs(data, required_columns)

  unit_stats <- .compute_unit_stats(
    data,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    y = y,
    never_treated_values = never_treated_values,
    d = d
  )
  obs_per_unit <- unit_stats$n_observed
  n_units <- length(obs_per_unit)
  n_periods <- if (n_units > 0L) as.integer(unit_stats$n_total_periods[[1L]]) else 0L
  min_obs_per_unit <- if (n_units > 0L) as.integer(min(obs_per_unit)) else 0L
  max_obs_per_unit <- if (n_units > 0L) as.integer(max(obs_per_unit)) else 0L
  mean_obs_per_unit <- if (n_units > 0L) mean(obs_per_unit) else 0
  std_obs_per_unit <- if (n_units > 1L) stats::sd(obs_per_unit) else 0
  treated_stats <- unit_stats[unit_stats$is_treated, , drop = FALSE]
  n_treated_units <- as.integer(nrow(treated_stats))
  units_below_demean_threshold <- if (n_treated_units > 0L) {
    as.integer(sum(!treated_stats$can_use_demean))
  } else {
    0L
  }
  units_below_detrend_threshold <- if (n_treated_units > 0L) {
    as.integer(sum(!treated_stats$can_use_detrend))
  } else {
    0L
  }

  list(
    is_balanced = data.table::uniqueN(obs_per_unit) <= 1L,
    n_units = as.integer(n_units),
    n_periods = n_periods,
    min_obs_per_unit = min_obs_per_unit,
    max_obs_per_unit = max_obs_per_unit,
    mean_obs_per_unit = as.numeric(mean_obs_per_unit),
    std_obs_per_unit = as.numeric(std_obs_per_unit),
    balance_ratio = if (max_obs_per_unit > 0L) {
      min_obs_per_unit / max_obs_per_unit
    } else {
      0
    },
    n_treated_units = n_treated_units,
    units_below_demean_threshold = units_below_demean_threshold,
    units_below_detrend_threshold = units_below_detrend_threshold,
    pct_usable_demean = if (n_treated_units > 0L) {
      100 * (1 - units_below_demean_threshold / n_treated_units)
    } else {
      100
    },
    pct_usable_detrend = if (n_treated_units > 0L) {
      100 * (1 - units_below_detrend_threshold / n_treated_units)
    } else {
      100
    }
  )
}


#' Assess selection-bias risk from observed-Y diagnostics
#'
#' @param missing_rates Named list returned by \code{.compute_missing_rates()}.
#' @param missing_pattern Named list returned by
#'   \code{.classify_missing_pattern()} or a scalar pattern label.
#' @param attrition Named list returned by \code{.compute_attrition_analysis()}.
#' @param balance Named list returned by \code{.compute_balance_statistics()}.
#' @return Named list with risk label, numeric score, warning messages, and
#'   factor contributions.
#' @keywords internal
.assess_selection_risk <- function(
    missing_rates,
    missing_pattern,
    attrition,
    balance
) {
  pattern <- if (is.list(missing_pattern)) {
    missing_pattern$pattern
  } else {
    missing_pattern
  }
  pattern <- toupper(as.character(pattern)[[1L]])

  missing_rate_overall <- as.numeric(missing_rates$overall_missing_rate %||% NA_real_)
  panel_is_balanced <- isTRUE(balance$is_balanced) &&
    isTRUE(missing_rates$is_balanced) &&
    is.finite(missing_rate_overall) &&
    missing_rate_overall == 0

  zero_factors <- list(
    missing_pattern = 0,
    attrition = 0,
    differential_attrition = 0,
    balance = 0
  )

  if (panel_is_balanced) {
    return(list(
      risk = "low",
      score = 0,
      warnings = character(0),
      factors = zero_factors
    ))
  }

  warnings <- character(0)
  factors <- zero_factors

  if (identical(pattern, "MAR")) {
    factors$missing_pattern <- 15
  } else if (identical(pattern, "MNAR")) {
    factors$missing_pattern <- 30
    warnings <- c(
      warnings,
      paste(
        "Missing data pattern suggests selection on unobservables.",
        "This may violate the selection mechanism assumption."
      )
    )
  } else if (identical(pattern, "UNKNOWN")) {
    factors$missing_pattern <- 10
  }

  attrition_rate <- as.numeric(attrition$overall_attrition %||% 0)
  if (attrition_rate >= 0.30) {
    factors$attrition <- 25
    warnings <- c(
      warnings,
      sprintf(
        "High attrition rate (%.1f%%). Consider using detrending which is more robust to selection on trends.",
        100 * attrition_rate
      )
    )
  } else if (attrition_rate >= 0.10) {
    factors$attrition <- 12
  }

  dropout_before <- as.numeric(attrition$dropout_before_treatment %||% 0)
  dropout_after <- as.numeric(attrition$dropout_after_treatment %||% 0)
  if (dropout_before > 0 && dropout_after > dropout_before * 2) {
    factors$differential_attrition <- 25
    warnings <- c(
      warnings,
      sprintf(
        paste(
          "Significantly more dropout after treatment (%d) than before (%d).",
          "This may indicate selection related to treatment effects."
        ),
        as.integer(dropout_after),
        as.integer(dropout_before)
      )
    )
  } else if (dropout_before > 0 && dropout_after > dropout_before * 1.5) {
    factors$differential_attrition <- 15
  }

  balance_ratio <- as.numeric(balance$balance_ratio %||% 0)
  if (balance_ratio > 0.8) {
    factors$balance <- 0
  } else if (balance_ratio > 0.5) {
    factors$balance <- 10
  } else {
    factors$balance <- 20
    warnings <- c(
      warnings,
      sprintf(
        "Low balance ratio (%.1f%%). Some units have much fewer observations than others.",
        100 * balance_ratio
      )
    )
  }

  score <- sum(unlist(factors), na.rm = TRUE)
  risk <- if (score < 25) {
    "low"
  } else if (score < 50) {
    "medium"
  } else {
    "high"
  }

  list(
    risk = risk,
    score = as.numeric(score),
    warnings = warnings,
    factors = factors
  )
}


#' Generate actionable recommendations from selection diagnostics
#'
#' @param risk Character scalar or list returned by \code{.assess_selection_risk()}.
#' @param pattern Named list returned by \code{.classify_missing_pattern()} or
#'   a scalar pattern label.
#' @param attrition Named list returned by \code{.compute_attrition_analysis()}.
#' @param balance Named list returned by \code{.compute_balance_statistics()}.
#' @return Character vector of recommendations.
#' @keywords internal
.generate_selection_recommendations <- function(
    risk,
    pattern,
    attrition,
    balance
) {
  risk_label <- if (is.list(risk)) risk$risk else risk
  risk_label <- toupper(as.character(risk_label)[[1L]])

  recommendations <- switch(
    risk_label,
    LOW = c(
      "Selection risk is low. Proceed with estimation. The selection mechanism assumption appears reasonable."
    ),
    MEDIUM = c(
      "Moderate selection risk detected. Consider the following:",
      "1. Use rolling='detrend' for additional robustness to selection on trends",
      "2. Compare results with a balanced subsample as a sensitivity check",
      "3. Report both demean and detrend results for transparency"
    ),
    HIGH = c(
      "High selection risk detected. Strongly recommend:",
      "1. Use rolling='detrend' method for greater robustness to selection on trends",
      "2. Conduct sensitivity analysis with a balanced subsample",
      "3. Report diagnostics and discuss potential selection bias",
      "4. Consider alternative identification strategies if possible"
    ),
    c("Selection diagnostics are inconclusive. Review the missing-data pattern before proceeding.")
  )

  pct_usable_detrend <- as.numeric(balance$pct_usable_detrend %||% 100)
  if (is.finite(pct_usable_detrend) && pct_usable_detrend < 90) {
    recommendations <- c(
      recommendations,
      sprintf(
        paste(
          "Note: Only %.1f%% of treated units have sufficient pre-treatment periods for detrending.",
          "Consider using demean if detrending excludes too many units."
        ),
        pct_usable_detrend
      )
    )
  }

  recommendations
}


#' Emit a concise or detailed selection-diagnostics summary
#'
#' @param x An object of class \code{lwdid_selection_diagnosis}.
#' @param detailed Logical; whether to print the detailed report.
#' @return \code{x} invisibly.
#' @keywords internal
.emit_selection_diagnosis_summary <- function(x, detailed = FALSE) {
  sep <- strrep("=", 60)

  cat(sep, "\n")
  cat("Selection Mechanism Diagnostics\n")
  cat(sep, "\n")
  cat(sprintf("Pattern: %s (confidence %.0f%%)\n", x$missing_pattern, 100 * x$missing_pattern_confidence))
  cat(sprintf("Overall Missing Rate: %.1f%%\n", 100 * x$missing_rate_overall))
  risk_display <- paste0(
    toupper(substr(x$selection_risk, 1L, 1L)),
    substr(x$selection_risk, 2L, nchar(x$selection_risk))
  )
  cat(sprintf("Selection Risk: %s\n", risk_display))
  cat(sprintf("Attrition Rate: %.1f%%\n", 100 * x$attrition_analysis$overall_attrition))
  cat(sprintf("Demean usable: %.1f%% of treated units\n", x$balance_statistics$pct_usable_demean))
  cat(sprintf("Detrend usable: %.1f%% of treated units\n", x$balance_statistics$pct_usable_detrend))

  if (isTRUE(detailed)) {
    cat(sprintf("Selection Risk Score: %.0f\n", x$selection_risk_score))
    cat("\nFactor Contributions:\n")
    cat(sprintf("  Missing pattern: %.0f\n", x$selection_risk_factors$missing_pattern))
    cat(sprintf("  Attrition: %.0f\n", x$selection_risk_factors$attrition))
    cat(sprintf(
      "  Differential attrition: %.0f\n",
      x$selection_risk_factors$differential_attrition
    ))
    cat(sprintf("  Balance: %.0f\n", x$selection_risk_factors$balance))
  }

  if (length(x$selection_tests) > 0L) {
    cat("\nStatistical Tests:\n")
    for (test in x$selection_tests) {
      status <- if (isTRUE(test$reject_null)) "REJECT" else "FAIL TO REJECT"
      cat(sprintf(
        "  - %s: statistic=%.4f, p-value=%.4f (%s)\n",
        test$test_name,
        test$statistic,
        test$pvalue,
        status
      ))
    }
  }

  if (length(x$warnings) > 0L) {
    cat("\nWarnings:\n")
    for (warning in x$warnings) {
      cat(sprintf("  - %s\n", warning))
    }
  }

  if (length(x$recommendations) > 0L) {
    cat("\nRecommendations:\n")
    for (recommendation in x$recommendations) {
      cat(sprintf("  - %s\n", recommendation))
    }
  }

  cat(sep, "\n")
  invisible(x)
}


#' Diagnose missing-data selection mechanisms in panel data
#'
#' @param data data.frame or data.table in long format.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param y Character outcome column.
#' @param gvar Character cohort column or \code{NULL}.
#' @param d Character treatment indicator or \code{NULL}.
#' @param covariates Optional character vector of MAR-test controls.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @param verbose Logical; whether to print a concise diagnosis summary.
#' @details
#' This function diagnoses whether missing outcomes are compatible with the
#' Lee-Wooldridge selection assumptions for unbalanced panels. It is a
#' diagnostic aid, not a repair for shock-related selection. Observed-outcome
#' availability drives the usability checks for demeaning and detrending, and
#' the public report can surface guidance such as
#' \code{Use rolling='detrend' for additional robustness to selection on trends}
#' and \code{Compare results with a balanced subsample as a sensitivity check}
#' when treated units lack enough observed pre-treatment outcomes.
#' @examples
#' example_panel <- data.frame(
#'   id = c(1, 1, 2, 2),
#'   t = c(1, 2, 1, 2),
#'   g = c(2, 2, 0, 0),
#'   y = c(NA, 5, 3, 4)
#' )
#'
#' diag <- diagnose_selection_mechanism(
#'   example_panel,
#'   y = "y",
#'   ivar = "id",
#'   tvar = "t",
#'   gvar = "g",
#'   verbose = FALSE
#' )
#'
#' print(diag)
#' summary(diag)
#' @return An object of class \code{lwdid_selection_diagnosis}.
#' @export
diagnose_selection_mechanism <- function(
    data,
    ivar,
    tvar,
    y,
    gvar = NULL,
    d = NULL,
    covariates = NULL,
    never_treated_values = c(0, Inf),
    verbose = TRUE
) {
  required_columns <- c(ivar, tvar, y)
  if (!is.null(gvar)) {
    required_columns <- c(required_columns, gvar)
  }
  if (!is.null(d)) {
    required_columns <- c(required_columns, d)
  }
  if (!is.null(covariates)) {
    required_columns <- c(required_columns, covariates)
  }
  .validate_selection_inputs(data, unique(required_columns))

  missing_rates <- .compute_missing_rates(data, ivar = ivar, tvar = tvar, y = y)
  missing_pattern <- .classify_missing_pattern(
    data,
    ivar = ivar,
    tvar = tvar,
    y = y,
    gvar = gvar,
    d = d,
    covariates = covariates
  )
  attrition <- .compute_attrition_analysis(
    data,
    ivar = ivar,
    tvar = tvar,
    y = y,
    gvar = gvar,
    d = d,
    never_treated_values = never_treated_values
  )
  balance <- .compute_balance_statistics(
    data,
    ivar = ivar,
    tvar = tvar,
    y = y,
    gvar = gvar,
    d = d,
    never_treated_values = never_treated_values
  )
  unit_stats <- as.data.frame(
    .compute_unit_stats(
      data,
      ivar = ivar,
      tvar = tvar,
      gvar = gvar,
      y = y,
      never_treated_values = never_treated_values,
      d = d
    ),
    stringsAsFactors = FALSE
  )
  risk <- .assess_selection_risk(
    missing_rates,
    missing_pattern,
    attrition,
    balance
  )
  recommendations <- .generate_selection_recommendations(
    risk$risk,
    missing_pattern,
    attrition,
    balance
  )
  missing_by_period <- .compute_missing_by_period(data, ivar = ivar, tvar = tvar, y = y)
  missing_by_cohort <- if (!is.null(gvar)) {
    .compute_missing_by_cohort(
      data,
      ivar = ivar,
      tvar = tvar,
      gvar = gvar,
      y = y,
      never_treated_values = never_treated_values
    )
  } else {
    list()
  }

  result <- structure(
    list(
      missing_rates = missing_rates,
      missing_pattern = missing_pattern$pattern,
      missing_pattern_confidence = as.numeric(missing_pattern$confidence),
      missing_pattern_description = missing_pattern$description,
      attrition_analysis = attrition,
      balance_statistics = balance,
      unit_stats = unit_stats,
      selection_risk = risk$risk,
      selection_risk_score = as.numeric(risk$score),
      selection_risk_factors = risk$factors,
      selection_tests = missing_pattern$tests,
      missing_rate_overall = as.numeric(missing_rates$overall_missing_rate),
      missing_rate_by_period = missing_by_period,
      missing_rate_by_cohort = missing_by_cohort,
      recommendations = recommendations,
      warnings = risk$warnings,
      attrition_rate = as.numeric(attrition$overall_attrition),
      attrition_by_period = attrition$attrition_by_period
    ),
    class = "lwdid_selection_diagnosis"
  )

  if (isTRUE(verbose)) {
    .emit_selection_diagnosis_summary(result, detailed = FALSE)
  }

  result
}


#' Get per-unit observed-Y missing-data statistics
#'
#' @param data data.frame or data.table in long format.
#' @param y Character outcome column.
#' @param ivar Character unit identifier column.
#' @param tvar Character time column.
#' @param gvar Character cohort column or \code{NULL}.
#' @param never_treated_values Numeric vector of explicit never-treated markers.
#' @details
#' The returned fields are computed from observed outcomes rather than row
#' presence. In particular, \code{n_pre_treatment}, \code{can_use_demean}, and
#' \code{can_use_detrend} reflect whether a treated unit has enough observed
#' pre-treatment outcomes to support the corresponding transformation.
#' @examples
#' example_panel <- data.frame(
#'   id = c(1, 1, 2, 2),
#'   t = c(1, 2, 1, 2),
#'   g = c(2, 2, 0, 0),
#'   y = c(NA, 5, 3, 4)
#' )
#'
#' get_unit_missing_stats(
#'   example_panel,
#'   y = "y",
#'   ivar = "id",
#'   tvar = "t",
#'   gvar = "g"
#' )
#' @return A data.frame with one row per unit and 17 diagnostic fields.
#' @export
get_unit_missing_stats <- function(
    data,
    y,
    ivar,
    tvar,
    gvar = NULL,
    never_treated_values = c(0, Inf)
) {
  as.data.frame(
    .compute_unit_stats(
      data,
      ivar = ivar,
      tvar = tvar,
      gvar = gvar,
      y = y,
      never_treated_values = never_treated_values
    ),
    stringsAsFactors = FALSE
  )
}


#' Print method for selection diagnosis objects
#'
#' @param x An object of class \code{lwdid_selection_diagnosis}.
#' @param ... Additional arguments (ignored).
#' @details
#' The printed report is a concise public summary. It surfaces the main risk
#' label alongside observed-outcome usability lines such as
#' \code{Demean usable: 0.0\% of treated units} and
#' \code{Detrend usable: 0.0\% of treated units}. When the observed
#' pre-treatment support is weak, the printed recommendations can include
#' \code{Use rolling='detrend' for additional robustness to selection on trends}
#' and \code{Compare results with a balanced subsample as a sensitivity check}.
#' @return \code{x} invisibly.
#' @export
print.lwdid_selection_diagnosis <- function(x, ...) {
  .emit_selection_diagnosis_summary(x, detailed = FALSE)
}


#' Summary method for selection diagnosis objects
#'
#' @param object An object of class \code{lwdid_selection_diagnosis}.
#' @param ... Additional arguments (ignored).
#' @details
#' The detailed summary adds the score breakdown on top of the concise report,
#' including lines such as \code{Selection Risk Score: 45},
#' \code{Attrition: 25}, and \code{Balance: 20}. It also preserves the
#' observed-outcome guidance around
#' \code{Use rolling='detrend' for additional robustness to selection on trends}
#' and \code{Compare results with a balanced subsample as a sensitivity check}.
#' @return \code{object} invisibly.
#' @export
summary.lwdid_selection_diagnosis <- function(object, ...) {
  .emit_selection_diagnosis_summary(object, detailed = TRUE)
}
