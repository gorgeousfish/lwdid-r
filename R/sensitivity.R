# ============================================================================
# sensitivity.R — Pre-period robustness sensitivity analysis for lwdid
#
# Internal helper functions for assessing ATT estimate robustness under
# varying pre-treatment period configurations.
#
# Sensitivity ratio thresholds:
#   highly_robust:      SR < 10%
#   moderately_robust:  10% <= SR < 25%
#   sensitive:          25% <= SR < 50%
#   highly_sensitive:   SR >= 50%
#
# Dependencies: conditions.R (stop_lwdid)
# ============================================================================


# ---- E8-01.1: Input Validation ---------------------------------------------

#' Validate inputs for sensitivity analysis
#'
#' Checks: (1) required columns exist, (2) rolling method valid,
#' (3) design mode consistency (gvar or d+post).
#'
#' @param data data.frame or data.table.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method.
#' @return NULL (invisible). Raises errors on failure.
#' @keywords internal
.validate_sensitivity_inputs <- function(
    data, y, ivar, tvar, gvar, d, post, rolling
) {
  required_cols <- c(y, ivar, tvar)
  if (!is.null(gvar)) required_cols <- c(required_cols, gvar)
  if (!is.null(d)) required_cols <- c(required_cols, d)
  if (!is.null(post)) required_cols <- c(required_cols, post)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      message = sprintf("Missing required columns: %s",
                        paste(missing_cols, collapse = ", ")),
      class = "lwdid_missing_column",
      column = missing_cols, available = names(data))
  }
  valid_rolling <- c("demean", "detrend", "demeanq", "detrendq")
  if (!is.character(rolling) || length(rolling) != 1L ||
      !tolower(rolling) %in% valid_rolling) {
    stop_lwdid(
      message = sprintf("rolling must be one of {%s}, got '%s'",
                        paste(valid_rolling, collapse = ", "),
                        as.character(rolling)),
      class = c("lwdid_invalid_rolling",
                "lwdid_invalid_parameter"),
      method = rolling, allowed = valid_rolling)
  }
  is_staggered <- !is.null(gvar)
  is_common <- !is.null(d) && !is.null(post)
  if (!is_staggered && !is_common) {
    stop_lwdid(
      message = paste0("Must specify either gvar (staggered) ",
                       "or both d and post (common timing)"),
      class = "lwdid_invalid_parameter",
      param = "mode", value = "incomplete",
      allowed = "gvar OR (d + post)")
  }
  invisible(NULL)
}


# ---- E8-01.2.1: Sensitivity Ratio ------------------------------------------

#' Compute sensitivity ratio
#'
#' SR = (max(ATT) - min(ATT)) / |ATT_baseline|
#'
#' @param atts numeric vector. ATT estimates.
#' @param baseline_att numeric(1). Baseline ATT.
#' @return numeric(1). Sensitivity ratio.
#' @keywords internal
.compute_sensitivity_ratio <- function(atts, baseline_att) {
  if (length(atts) == 0L) return(0)
  att_range <- max(atts) - min(atts)
  if (abs(baseline_att) > 1e-10) {
    return(att_range / abs(baseline_att))
  } else {
    if (att_range > 1e-10) return(Inf)
    return(0)
  }
}


# ---- E8-01.2.2: Robustness Level -------------------------------------------

#' Determine robustness level from sensitivity ratio
#'
#' @param sensitivity_ratio numeric(1).
#' @return character(1). Robustness level string.
#' @keywords internal
.determine_robustness_level <- function(sensitivity_ratio) {
  if (sensitivity_ratio < 0.10) return("highly_robust")
  if (sensitivity_ratio < 0.25) return("moderately_robust")
  if (sensitivity_ratio < 0.50) return("sensitive")
  return("highly_sensitive")
}


# ---- E8-01.2.3: Robustness Recommendations ---------------------------------

#' Generate robustness recommendations
#'
#' Four mutually exclusive main branches, independent
#' method-specific checks, and monotonic trend detection.
#'
#' @param specifications list of spec results (each with
#'   \code{$att}, \code{$converged}, \code{$n_pre_periods}).
#' @param baseline_spec list. Baseline spec result.
#' @param sensitivity_ratio numeric(1).
#' @param is_robust logical(1).
#' @param all_same_sign logical(1).
#' @param all_significant logical(1).
#' @param rolling character(1). Transformation method.
#' @return list(main_rec, detailed_recommendations,
#'   result_warnings)
#' @keywords internal
.generate_robustness_recommendations <- function(
    specifications, baseline_spec, sensitivity_ratio,
    is_robust, all_same_sign, all_significant, rolling
) {
  recs <- character(0)
  warns <- character(0)
  # Branch 1: robust + same sign + all significant
  if (is_robust && all_same_sign && all_significant) {
    main_rec <- paste0(
      "Results are robust to pre-treatment period ",
      "selection. The ATT estimate is stable across ",
      "specifications.")
  # Branch 2: robust + same sign (not all significant)
  } else if (is_robust && all_same_sign) {
    main_rec <- paste0(
      "Results are moderately robust. Sign is ",
      "consistent but significance varies across ",
      "specifications.")
    recs <- c(recs,
      paste0("Consider reporting the range of estimates ",
             "for transparency."))
  # Branch 3: sign changes
  } else if (!all_same_sign) {
    main_rec <- paste0(
      "CAUTION: Results are sensitive to ",
      "pre-treatment period selection. ",
      "Sign changes detected across specifications.")
    warns <- c(warns,
      paste0("Sign change detected - interpret results ",
             "with caution."))
    recs <- c(recs,
      paste0("Investigate why estimates change sign ",
             "with different pre-periods."),
      paste0("Consider using detrend method if trends ",
             "may be heterogeneous."))
  # Branch 4: else (moderate sensitivity)
  } else {
    main_rec <- sprintf(
      paste0("Results show moderate sensitivity ",
             "(ratio = %.1f%%). ",
             "Consider additional robustness checks."),
      sensitivity_ratio * 100)
  }
  # Method-specific (independent, always run)
  rolling_lower <- tolower(rolling)
  if (rolling_lower == "demean" && sensitivity_ratio > 0.25) {
    recs <- c(recs,
      paste0("High sensitivity with demean suggests ",
             "potential heterogeneous trends. ",
             "Consider using rolling='detrend' instead."))
  }
  if (rolling_lower == "detrend" && sensitivity_ratio > 0.50) {
    recs <- c(recs,
      paste0("High sensitivity even with detrend ",
             "suggests potential model misspecification ",
             "or data quality issues."))
  }
  # Monotonic trend detection (>=3 converged specs)
  conv <- Filter(function(s) isTRUE(s$converged), specifications)
  if (length(conv) >= 3L) {
    npre <- vapply(conv, function(s) s$n_pre_periods, numeric(1))
    conv <- conv[order(npre)]
    atts <- vapply(conv, function(s) s$att, numeric(1))
    diffs <- diff(atts)
    if (length(diffs) > 0L &&
        (all(diffs > 0) || all(diffs < 0))) {
      warns <- c(warns,
        paste0("ATT estimates show monotonic trend with ",
               "pre-period count. This may indicate ",
               "time-varying confounding."))
    }
  }
  list(main_rec = main_rec,
       detailed_recommendations = recs,
       result_warnings = warns)
}


# ---- E8-01.2.4: Auto-detect Pre-period Range -------------------------------

#' Auto-detect valid pre-treatment period range
#'
#' @param data data.frame. Panel data.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method.
#' @return integer(2). c(min_pre, max_pre).
#' @keywords internal
.auto_detect_pre_period_range <- function(
    data, ivar, tvar, gvar, d, post, rolling
) {
  min_req <- c(demean = 1L, detrend = 2L,
               demeanq = 1L, detrendq = 2L)
  rl <- tolower(rolling)
  min_pre <- if (rl %in% names(min_req)) min_req[[rl]] else 1L
  if (!is.null(gvar)) {
    cohorts <- unique(data[[gvar]])
    cohorts <- cohorts[!is.na(cohorts)]
    cohorts <- cohorts[is.finite(cohorts) & cohorts > 0]
    if (length(cohorts) == 0L) {
      return(as.integer(c(min_pre, min_pre)))
    }
    min_time <- min(data[[tvar]], na.rm = TRUE)
    max_pre_by_cohort <- as.integer(cohorts - min_time)
    max_pre <- min(max_pre_by_cohort)
  } else {
    pre_data <- data[data[[post]] == 0, , drop = FALSE]
    max_pre <- length(unique(pre_data[[tvar]]))
  }
  max_pre <- max(max_pre, min_pre)
  as.integer(c(min_pre, max_pre))
}


# ---- E8-01.3.1+3.2: Filter to N Pre-periods --------------------------------

#' Filter data to use specified number of pre-treatment periods
#'
#' For staggered designs, filters cohort-by-cohort.
#' For common timing, filters the entire dataset.
#'
#' @param data data.frame. Panel data.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param n_pre_periods integer(1). Number of pre-treatment periods.
#' @param exclude_periods integer(1). Periods to exclude before treatment.
#' @return data.frame. Filtered data.
#' @keywords internal
.filter_to_n_pre_periods <- function(
    data, ivar, tvar, gvar, d, post,
    n_pre_periods, exclude_periods
) {
  if (!is.null(gvar)) {
    # ---- Staggered mode ----
    filtered_list <- list()
    gvar_vals <- data[[gvar]]

    # Identify valid cohorts (non-NA, finite, > 0)
    cohorts <- unique(gvar_vals)
    cohorts <- cohorts[!is.na(cohorts)]
    cohorts <- cohorts[is.finite(cohorts) & cohorts > 0]

    for (g in cohorts) {
      cohort_data <- data[which(gvar_vals == g), , drop = FALSE]
      pre_end <- g - 1L - exclude_periods
      pre_start <- pre_end - n_pre_periods + 1L

      # Keep [pre_start, pre_end] union [g, Inf)
      tvec <- cohort_data[[tvar]]
      time_mask <- (tvec >= pre_start & tvec <= pre_end) |
        (tvec >= g)
      filtered_list[[length(filtered_list) + 1L]] <-
        cohort_data[time_mask, , drop = FALSE]
    }

    # Never-treated units (gvar is NA, 0, or Inf)
    nt_mask <- is_never_treated(gvar_vals)
    if (any(nt_mask)) {
      never_data <- data[nt_mask, , drop = FALSE]
      if (length(cohorts) > 0L) {
        min_cohort <- min(cohorts)
        pre_end_nt <- min_cohort - 1L - exclude_periods
        pre_start_nt <- pre_end_nt - n_pre_periods + 1L
        nt_time_mask <- never_data[[tvar]] >= pre_start_nt
        filtered_list[[length(filtered_list) + 1L]] <-
          never_data[nt_time_mask, , drop = FALSE]
      } else {
        filtered_list[[length(filtered_list) + 1L]] <- never_data
      }
    }

    if (length(filtered_list) > 0L) {
      return(do.call(rbind, filtered_list))
    }
    # Empty data.frame preserving column structure
    return(data[0L, , drop = FALSE])

  } else {
    # ---- Common Timing mode ----
    pre_data <- data[data[[post]] == 0, , drop = FALSE]
    post_data <- data[data[[post]] != 0, , drop = FALSE]

    pre_times <- sort(unique(pre_data[[tvar]]))

    # Exclude last exclude_periods pre-treatment periods
    if (exclude_periods > 0L) {
      if (length(pre_times) > exclude_periods) {
        n_keep <- length(pre_times) - exclude_periods
        pre_times <- pre_times[seq_len(n_keep)]
      } else {
        # Edge guard: exclude_periods >= available pre-periods
        warn_lwdid(
          message = sprintf(
            paste0("Cannot exclude %d periods: only %d ",
                   "pre-treatment periods available. ",
                   "Using all available periods."),
            exclude_periods, length(pre_times)),
          class = "lwdid_data")
        # Keep all pre_times (don't exclude any)
      }
    }

    # Select last n_pre_periods
    if (length(pre_times) <= n_pre_periods) {
      selected_pre_times <- pre_times
    } else {
      selected_pre_times <- tail(pre_times, n_pre_periods)
    }

    filtered_pre <- pre_data[
      pre_data[[tvar]] %in% selected_pre_times, , drop = FALSE]
    rbind(filtered_pre, post_data)
  }
}


# ---- E8-01.3.3: Run Single Specification ------------------------------------

#' Run estimation for a single specification
#'
#' Filters data to the specified number of pre-treatment periods
#' and runs lwdid estimation. Failures are caught and returned
#' as converged=FALSE.
#'
#' @param data data.frame. Panel data.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param d character(1) or NULL. Treatment indicator.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method.
#' @param estimator character(1). Estimation method.
#' @param controls character vector or NULL. Control variables.
#' @param vce character(1) or NULL. Variance estimator.
#' @param cluster_var character(1) or NULL. Cluster variable.
#' @param n_pre_periods integer(1). Pre-treatment periods to use.
#' @param exclude_periods integer(1). Periods to exclude.
#' @param alpha numeric(1). Significance level.
#' @param spec_id integer(1). Specification identifier.
#' @param extra_args list. Additional arguments forwarded to \code{lwdid()}.
#' @return list. SpecificationResult.
#' @keywords internal
.run_single_specification <- function(
    data, y, ivar, tvar, gvar, d, post,
    rolling, estimator, controls, vce, cluster_var,
    n_pre_periods, exclude_periods, alpha, spec_id,
    extra_args = list()
) {
  tryCatch({
    # 1. Filter data
    filtered <- .filter_to_n_pre_periods(
      data = data, ivar = ivar, tvar = tvar,
      gvar = gvar, d = d, post = post,
      n_pre_periods = n_pre_periods,
      exclude_periods = exclude_periods)

    if (nrow(filtered) == 0L) {
      stop("No data remaining after filtering")
    }

    # 2. Determine start/end periods
    if (!is.null(gvar)) {
      filt_gvar <- filtered[[gvar]]
      cohorts_f <- unique(filt_gvar[!is.na(filt_gvar)])
      cohorts_f <- cohorts_f[is.finite(cohorts_f) & cohorts_f > 0]
      if (length(cohorts_f) > 0L) {
        min_cohort <- min(cohorts_f)
        end_period <- as.integer(min_cohort - 1L - exclude_periods)
        start_period <- end_period - n_pre_periods + 1L
      } else {
        start_period <- 0L
        end_period <- 0L
      }
    } else {
      pre_mask <- filtered[[post]] == 0
      pre_times <- filtered[[tvar]][pre_mask]
      if (length(pre_times) > 0L) {
        start_period <- as.integer(min(pre_times))
        end_period <- as.integer(max(pre_times))
      } else {
        start_period <- 0L
        end_period <- 0L
      }
    }

    # 3. Run lwdid estimation
    result <- do.call(
      lwdid,
      c(
        list(
          data = filtered, y = y, ivar = ivar, tvar = tvar,
          gvar = gvar, d = d, post = post,
          rolling = rolling, estimator = estimator,
          controls = controls, vce = vce,
          cluster_var = cluster_var, alpha = alpha,
          verbose = "quiet"
        ),
        .filter_sensitivity_args(extra_args, lwdid)
      )
    )

    # 4. Return SpecificationResult
    list(
      specification_id = as.integer(spec_id),
      n_pre_periods = as.integer(n_pre_periods),
      start_period = as.integer(start_period),
      end_period = as.integer(end_period),
      excluded_periods = as.integer(exclude_periods),
      att = result$att,
      se = result$se_att,
      t_stat = result$t_stat,
      pvalue = result$pvalue,
      ci_lower = result$ci_lower,
      ci_upper = result$ci_upper,
      n_treated = result$n_treated,
      n_control = result$n_control,
      df = result$df_inference,
      converged = TRUE,
      spec_warnings = character(0))

  }, error = function(e) {
    warning(sprintf(
      "Specification %d (n_pre=%d) failed: %s",
      spec_id, n_pre_periods, conditionMessage(e)))
    list(
      specification_id = as.integer(spec_id),
      n_pre_periods = as.integer(n_pre_periods),
      start_period = 0L,
      end_period = 0L,
      excluded_periods = as.integer(exclude_periods),
      att = NA_real_,
      se = NA_real_,
      t_stat = NA_real_,
      pvalue = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      n_treated = 0L,
      n_control = 0L,
      df = 0L,
      converged = FALSE,
      spec_warnings = conditionMessage(e))
  })
}


# ---- E8-01.4: Main Function ------------------------------------------------

#' Pre-treatment Period Robustness Analysis
#'
#' Assess robustness of ATT estimates to pre-treatment period selection.
#' Tests how ATT estimates vary when using different numbers of
#' pre-treatment periods for the unit-specific transformation.
#'
#' @param data data.frame or data.table. Panel data in long format.
#' @param y character(1). Outcome variable column name.
#' @param ivar character(1). Unit identifier column name.
#' @param tvar character(1). Time variable column name.
#' @param gvar character(1) or NULL. Cohort variable (staggered).
#' @param d character(1) or NULL. Treatment indicator (common timing).
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method. Default "demean".
#' @param estimator character(1). Estimation method. Default "ra".
#' @param controls character vector or NULL. Control variables.
#' @param vce character(1) or NULL. Variance estimator type.
#' @param cluster_var character(1) or NULL. Cluster variable.
#' @param pre_period_range integer(2) or NULL. Range of pre-treatment
#'   periods to test c(min, max). NULL for auto-detection.
#' @param step integer(1). Step size. Default 1L.
#' @param exclude_periods_before_treatment integer(1). Periods to
#'   exclude before treatment. Default 0L.
#' @param robustness_threshold numeric(1). SR threshold for robustness.
#'   Default 0.25.
#' @param alpha numeric(1). Significance level. Default 0.05.
#' @param verbose logical(1). Print progress. Default TRUE.
#' @param ... Additional arguments forwarded to underlying \code{lwdid()} calls.
#' @return An S3 object of class \code{lwdid_sensitivity}.
#' @export
lwdid_sensitivity_pre_period <- function(
    data, y, ivar, tvar,
    gvar = NULL, d = NULL, post = NULL,
    rolling = "demean", estimator = "ra",
    controls = NULL, vce = NULL, cluster_var = NULL,
    pre_period_range = NULL, step = 1L,
    exclude_periods_before_treatment = 0L,
    robustness_threshold = 0.25, alpha = 0.05,
    verbose = TRUE,
    ...
) {
  extra_args <- list(...)

  # Step 1: Input validation
  .validate_sensitivity_inputs(
    data, y, ivar, tvar, gvar, d, post, rolling)

  # Step 2: Determine pre-period range
  if (is.null(pre_period_range)) {
    pre_period_range <- .auto_detect_pre_period_range(
      data, ivar, tvar, gvar, d, post, rolling)
  } else {
    if (pre_period_range[1L] > pre_period_range[2L]) {
      stop_lwdid(
        message = sprintf(
          paste0("Invalid pre_period_range: min (%d) > max (%d). ",
                 "Ensure min_pre <= max_pre."),
          pre_period_range[1L], pre_period_range[2L]),
        class = "lwdid_invalid_parameter",
        param = "pre_period_range",
        value = pre_period_range,
        allowed = "c(min_pre, max_pre) where min_pre <= max_pre")
    }
  }

  min_pre <- pre_period_range[1L]
  max_pre <- pre_period_range[2L]

  if (verbose) {
    message("Pre-treatment period robustness analysis")
    message(sprintf("Testing pre-period range: %d to %d", min_pre, max_pre))
    message(strrep("-", 50L))
  }

  # Step 3: Generate specification list
  n_pre_values <- seq(min_pre, max_pre, by = step)

  if (length(n_pre_values) < 2L) {
    warning(sprintf(
      paste0("Only %d specification(s) possible. ",
             "Consider expanding pre_period_range or reducing step."),
      length(n_pre_values)))
  }

  # Step 4: Run estimations
  specifications <- vector("list", length(n_pre_values))
  for (i in seq_along(n_pre_values)) {
    if (verbose) {
      message(sprintf(
        "Running specification %d/%d: n_pre=%d",
        i, length(n_pre_values), n_pre_values[i]))
    }
    specifications[[i]] <- .run_single_specification(
      data = data, y = y, ivar = ivar, tvar = tvar,
      gvar = gvar, d = d, post = post,
      rolling = rolling, estimator = estimator,
      controls = controls, vce = vce,
      cluster_var = cluster_var,
      n_pre_periods = n_pre_values[i],
      exclude_periods = exclude_periods_before_treatment,
      alpha = alpha, spec_id = i,
      extra_args = extra_args)
  }

  # Step 5: Identify baseline (max n_pre converged spec)
  converged_specs <- Filter(
    function(s) isTRUE(s$converged), specifications)

  if (length(converged_specs) == 0L) {
    stop_lwdid(
      message = "All specifications failed to converge",
      class = "lwdid_convergence")
  }

  npre_vals <- vapply(
    converged_specs, function(s) s$n_pre_periods, integer(1))
  baseline_spec <- converged_specs[[which.max(npre_vals)]]

  # Step 6: Compute sensitivity metrics
  atts <- vapply(
    converged_specs, function(s) s$att, numeric(1))
  att_range <- c(min(atts), max(atts))
  att_mean <- mean(atts)
  att_std <- if (length(atts) > 1L) sd(atts) else 0.0

  sensitivity_ratio <- .compute_sensitivity_ratio(
    atts, baseline_spec$att)

  # Step 7: Assess robustness
  robustness_level <- .determine_robustness_level(sensitivity_ratio)
  is_robust <- sensitivity_ratio < robustness_threshold

  # Step 8: Sign and significance stability
  baseline_sign <- sign(baseline_spec$att)
  all_same_sign <- all(sign(atts) == baseline_sign)
  pvals <- vapply(
    converged_specs, function(s) s$pvalue, numeric(1))
  n_significant <- sum(pvals < alpha)
  all_significant <- n_significant == length(converged_specs)
  n_sign_changes <- sum(sign(atts) != baseline_sign)

  # Step 9: Generate recommendations
  rec <- .generate_robustness_recommendations(
    specifications = specifications,
    baseline_spec = baseline_spec,
    sensitivity_ratio = sensitivity_ratio,
    is_robust = is_robust,
    all_same_sign = all_same_sign,
    all_significant = all_significant,
    rolling = rolling)

  # Step 10: Construct lwdid_sensitivity object
  result <- structure(list(
    type = "pre_period",
    specifications = specifications,
    baseline_spec = baseline_spec,
    att_range = att_range,
    att_mean = att_mean,
    att_std = att_std,
    sensitivity_ratio = sensitivity_ratio,
    robustness_level = robustness_level,
    is_robust = is_robust,
    robustness_threshold = robustness_threshold,
    all_same_sign = all_same_sign,
    all_significant = all_significant,
    n_significant = as.integer(n_significant),
    n_sign_changes = as.integer(n_sign_changes),
    rolling_method = rolling,
    estimator = estimator,
    n_specifications = length(specifications),
    pre_period_range_tested = as.integer(c(min_pre, max_pre)),
    recommendation = rec$main_rec,
    detailed_recommendations = rec$detailed_recommendations,
    result_warnings = rec$result_warnings
  ), class = c("lwdid_sensitivity", "list"))

  if (verbose) {
    message("")
    print(result)
  }

  result
}


# ---- E8-01.5: S3 Methods ---------------------------------------------------

#' Print method for lwdid_sensitivity
#'
#' @param x An lwdid_sensitivity object.
#' @param ... Additional arguments (ignored).
#' @return x invisibly.
#' @export
print.lwdid_sensitivity <- function(x, ...) {
  if (x$type == "no_anticipation") {
    cat("\n=== No-Anticipation Sensitivity Analysis ===\n")
    cat(sprintf("Rolling method: %s | Estimator: %s\n", x$rolling_method, x$estimator))
    cat(sprintf("Anticipation detected: %s\n", ifelse(x$anticipation_detected, "YES", "NO")))
    cat(sprintf("Detection method: %s\n", x$detection_method))
    cat(sprintf("Recommended exclusion: %d periods\n", x$recommended_exclusion))
    cat(sprintf("Estimates: %d total, %d converged\n", x$n_estimates, x$n_converged))
    cat(sprintf("\nRecommendation: %s\n", x$recommendation))
    return(invisible(x))
  }
  cat("Pre-treatment Period Robustness Analysis\n")
  cat(sprintf("  Transformation: %s | Estimator: %s\n",
              x$rolling_method, x$estimator))
  cat(sprintf("  Sensitivity Ratio: %.1f%% (%s)\n",
              x$sensitivity_ratio * 100,
              gsub("_", " ", x$robustness_level)))
  cat(sprintf("  Baseline ATT: %.4f (n_pre=%d)\n",
              x$baseline_spec$att,
              x$baseline_spec$n_pre_periods))
  cat(sprintf("  Recommendation: %s\n", x$recommendation))
  invisible(x)
}

#' Summary method for lwdid_sensitivity
#'
#' Detailed report matching Python PrePeriodRobustnessResult.summary().
#'
#' @param object An lwdid_sensitivity object.
#' @param ... Additional arguments (ignored).
#' @return object invisibly.
#' @export
summary.lwdid_sensitivity <- function(object, ...) {
  if (object$type == "no_anticipation") {
    cat("\n", strrep("=", 60), "\n")
    cat("  No-Anticipation Sensitivity Analysis Report\n")
    cat(strrep("=", 60), "\n\n")
    cat("Configuration:\n")
    cat(sprintf("  Rolling method: %s\n", object$rolling_method))
    cat(sprintf("  Estimator: %s\n", object$estimator))
    cat(sprintf("  Detection threshold: %.2f\n", object$detection_threshold))
    cat(sprintf("  Max anticipation tested: %d\n\n", object$max_anticipation_tested))
    cat("Estimation Results:\n")
    cat(sprintf("  %-10s %-12s %-12s %-12s %-6s\n", "Excluded", "ATT", "SE", "p-value", "Sig"))
    cat(strrep("-", 55), "\n")
    for (est in object$estimates) {
      if (isTRUE(est$converged)) {
        sig <- if (est$pvalue < 0.05) "***" else if (est$pvalue < 0.10) "*" else ""
        cat(sprintf("  %-10d %-12.4f %-12.4f %-12.4f %-6s\n",
                    est$excluded_periods, est$att, est$se, est$pvalue, sig))
      } else {
        cat(sprintf("  %-10d %-12s %-12s %-12s %-6s\n",
                    est$excluded_periods, "FAILED", "-", "-", ""))
      }
    }
    cat(sprintf("\nDetection Result: %s\n", object$detection_method))
    cat(sprintf("Anticipation Detected: %s\n", ifelse(object$anticipation_detected, "YES", "NO")))
    if (object$anticipation_detected) {
      cat(sprintf("Recommended Exclusion: %d periods\n", object$recommended_exclusion))
    }
    cat(sprintf("\nRecommendation: %s\n", object$recommendation))
    if (length(object$result_warnings) > 0L) {
      cat("\nWarnings:\n")
      for (w in object$result_warnings) cat(sprintf("  - %s\n", w))
    }
    return(invisible(object))
  }
  sep <- strrep("=", 60)
  subsep <- strrep("-", 60)

  cat(sep, "\n")
  cat("Pre-treatment Period Robustness Analysis\n")
  cat(sep, "\n\n")

  # Configuration
  cat("Configuration:\n")
  cat(sprintf("  Transformation method: %s\n", object$rolling_method))
  cat(sprintf("  Estimator: %s\n", object$estimator))
  cat(sprintf("  Pre-period range tested: %d to %d\n",
              object$pre_period_range_tested[1L],
              object$pre_period_range_tested[2L]))
  cat(sprintf("  Number of specifications: %d\n",
              object$n_specifications))
  cat("\n")

  # Baseline estimate
  cat(subsep, "\n")
  cat("Baseline Estimate:\n")
  bs <- object$baseline_spec
  cat(sprintf("  ATT: %.6f\n", bs$att))
  cat(sprintf("  SE:  %.6f\n", bs$se))
  cat(sprintf("  t-stat: %.4f\n", bs$t_stat))
  cat(sprintf("  p-value: %.4f\n", bs$pvalue))
  cat(sprintf("  95%% CI: [%.6f, %.6f]\n", bs$ci_lower, bs$ci_upper))
  cat(sprintf("  n_pre_periods: %d\n", bs$n_pre_periods))
  cat("\n")

  # Sensitivity metrics
  cat(subsep, "\n")
  cat("Sensitivity Metrics:\n")
  cat(sprintf("  ATT range: [%.6f, %.6f]\n",
              object$att_range[1L], object$att_range[2L]))
  cat(sprintf("  ATT mean: %.6f\n", object$att_mean))
  cat(sprintf("  ATT std: %.6f\n", object$att_std))
  cat(sprintf("  Sensitivity ratio: %.4f (%.1f%%)\n",
              object$sensitivity_ratio,
              object$sensitivity_ratio * 100))
  cat("\n")

  # Robustness assessment
  cat(subsep, "\n")
  cat("Robustness Assessment:\n")
  cat(sprintf("  Level: %s\n",
              gsub("_", " ", object$robustness_level)))
  cat(sprintf("  Is robust (SR < %.0f%%): %s\n",
              object$robustness_threshold * 100,
              if (object$is_robust) "YES" else "NO"))
  cat(sprintf("  All same sign: %s\n",
              if (object$all_same_sign) "YES" else "NO"))
  cat(sprintf("  All significant: %s (%d/%d)\n",
              if (object$all_significant) "YES" else "NO",
              object$n_significant,
              length(Filter(
                function(s) isTRUE(s$converged),
                object$specifications))))
  cat(sprintf("  Sign changes: %d\n", object$n_sign_changes))
  cat("\n")

  # Specification details table
  cat(subsep, "\n")
  cat("Specification Details:\n")
  cat(sprintf("  %-8s %-12s %-12s %-10s %-5s %-5s\n",
              "n_pre", "ATT", "SE", "p-value", "Sig", "Base"))
  cat(sprintf("  %s\n", strrep("-", 52)))
  for (s in object$specifications) {
    if (isTRUE(s$converged)) {
      sig_mark <- if (!is.na(s$pvalue) && s$pvalue < 0.05) "*" else ""
      base_mark <- if (s$n_pre_periods ==
                       object$baseline_spec$n_pre_periods) "<-" else ""
      cat(sprintf("  %-8d %-12.6f %-12.6f %-10.4f %-5s %-5s\n",
                  s$n_pre_periods, s$att, s$se,
                  s$pvalue, sig_mark, base_mark))
    } else {
      cat(sprintf("  %-8d %-12s %-12s %-10s %-5s %-5s\n",
                  s$n_pre_periods, "FAILED", "-", "-", "", ""))
    }
  }
  cat("\n")

  # Recommendations
  cat(subsep, "\n")
  cat("Recommendation:\n")
  cat(sprintf("  %s\n", object$recommendation))
  if (length(object$detailed_recommendations) > 0L) {
    cat("\nDetailed Recommendations:\n")
    for (r in object$detailed_recommendations) {
      cat(sprintf("  - %s\n", r))
    }
  }
  if (length(object$result_warnings) > 0L) {
    cat("\nWarnings:\n")
    for (w in object$result_warnings) {
      cat(sprintf("  ! %s\n", w))
    }
  }
  cat(sep, "\n")

  invisible(object)
}

#' Plot method for lwdid_sensitivity
#'
#' Specification curve plot with CI error bars, baseline reference,
#' and robustness band. Requires ggplot2.
#'
#' @param x An lwdid_sensitivity object.
#' @param type Optional plot type override.
#' @param ci_level Confidence level used to interpret significance colors.
#' @param show_threshold Logical; whether to show the fixed +/-25% robustness band.
#' @param show_ci Logical; whether to show confidence interval error bars.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.lwdid_sensitivity <- function(x, type = NULL, ci_level = 0.95,
                                   show_threshold = TRUE, show_ci = TRUE, ...) {
  stopifnot(inherits(x, "lwdid_sensitivity"))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' required for plotting. Install via: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  if (is.null(type)) {
    type <- x$type
  }

  if (identical(type, "no_anticipation")) {
    converged <- Filter(function(e) isTRUE(e$converged), x$estimates)
    if (length(converged) == 0L) {
      warning("All estimates failed to converge, cannot plot valid results", call. = FALSE)
      return(
        ggplot2::ggplot() +
          ggplot2::labs(
            title = "No-Anticipation Sensitivity",
            subtitle = "No converged estimates"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.subtitle = ggplot2::element_text(color = "red")
          )
      )
    }

    plot_df <- data.frame(
      excluded = vapply(converged, function(e) e$excluded_periods, integer(1)),
      att = vapply(converged, function(e) e$att, numeric(1)),
      ci_lower = vapply(converged, function(e) e$ci_lower, numeric(1)),
      ci_upper = vapply(converged, function(e) e$ci_upper, numeric(1))
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = excluded, y = att)) +
      ggplot2::geom_hline(
        yintercept = x$baseline_estimate$att,
        linetype = "dashed",
        color = "gray50"
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "gray80")

    if (show_ci) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        width = 0.2,
        color = "gray40"
      )
    }

    p <- p +
      ggplot2::geom_point(size = 3, color = "#2166AC") +
      ggplot2::geom_line(color = "#2166AC", alpha = 0.5) +
      ggplot2::scale_x_continuous(breaks = plot_df$excluded)

    if (isTRUE(x$anticipation_detected) && x$recommended_exclusion > 0L) {
      rec_row <- plot_df[plot_df$excluded == x$recommended_exclusion, , drop = FALSE]
      if (nrow(rec_row) > 0L) {
        p <- p + ggplot2::geom_point(
          data = rec_row,
          color = "red",
          size = 5,
          shape = 1,
          stroke = 1.5
        )
      }
    }

    return(
      p +
        ggplot2::labs(
          title = sprintf(
            "No-Anticipation Sensitivity (%s)",
            ifelse(x$anticipation_detected, "Detected", "Not Detected")
          ),
          x = "Excluded Pre-treatment Periods",
          y = "ATT Estimate"
        ) +
        ggplot2::theme_minimal()
    )
  }

  if (!identical(type, "pre_period")) {
    stop(sprintf("Unknown sensitivity analysis type: '%s'", type), call. = FALSE)
  }

  alpha <- 1 - ci_level
  conv <- Filter(function(s) isTRUE(s$converged), x$specifications)

  if (length(conv) == 0L) {
    warning("All specifications failed to converge, cannot plot valid results", call. = FALSE)
    return(
      ggplot2::ggplot() +
        ggplot2::labs(
          title = "Pre-period Robustness",
          subtitle = "No converged specifications"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(color = "red")
        )
    )
  }

  plot_df <- data.frame(
    n_pre = vapply(conv, function(s) s$n_pre_periods, integer(1)),
    att = vapply(conv, function(s) s$att, numeric(1)),
    ci_lower = vapply(conv, function(s) s$ci_lower, numeric(1)),
    ci_upper = vapply(conv, function(s) s$ci_upper, numeric(1)),
    significant = vapply(conv, function(s) {
      !is.na(s$pvalue) && s$pvalue < alpha
    }, logical(1))
  )

  baseline_att <- x$baseline_spec$att
  band_half <- 0.25 * abs(baseline_att)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = n_pre, y = att))

  if (show_threshold) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = baseline_att - band_half,
      ymax = baseline_att + band_half,
      fill = "gray90",
      alpha = 0.5
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, color = "gray80", linewidth = 0.5) +
    ggplot2::geom_hline(
      yintercept = baseline_att,
      color = "red",
      linetype = "dashed",
      linewidth = 0.7
    )

  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      color = "gray40"
    )
  }

  p +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 3) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = "Significance"
    ) +
    ggplot2::scale_x_continuous(breaks = plot_df$n_pre) +
    ggplot2::labs(
      title = sprintf(
        "Pre-period Robustness (SR=%.1f%%, %s)",
        x$sensitivity_ratio * 100,
        gsub("_", " ", x$robustness_level)
      ),
      x = "Number of Pre-treatment Periods",
      y = "ATT Estimate"
    ) +
    ggplot2::theme_minimal()
}


# ---- E8-02.1: Get Maximum Pre-treatment Periods ----------------------------

#' Get maximum pre-treatment periods available
#'
#' Staggered: returns minimum across all cohorts of sum(periods < g).
#' Common Timing: returns count of unique pre-treatment periods.
#'
#' @details
#'   Paper reference: lw2025 Section 4.4
#'   Python reference: lwdid-py sensitivity.py _get_max_pre_periods()
#'   R enhancement: uses sum(periods < g) instead of Python's int(c - min_time),
#'   more robust for non-consecutive time indices.
#'
#' @param data data.frame.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @return integer(1). Maximum pre-treatment periods.
#' @keywords internal
.get_max_pre_periods <- function(data, ivar, tvar, gvar, post) {
  periods <- sort(unique(data[[tvar]]))

  if (!is.null(gvar)) {
    # Staggered: actual unique period count per cohort [R enhancement]
    cohorts <- unique(data[[gvar]])
    cohorts <- cohorts[!is.na(cohorts) & is.finite(cohorts) & cohorts > 0]

    if (length(cohorts) == 0L) return(0L)

    max_pre_by_cohort <- vapply(
      cohorts,
      function(g) sum(periods < g),
      integer(1)
    )
    return(min(max_pre_by_cohort))

  } else if (!is.null(post)) {
    # Common Timing: unique pre-treatment period count
    return(length(unique(data[[tvar]][data[[post]] == 0])))

  } else {
    # Fallback
    return(length(periods) %/% 2L)
  }
}


# ---- E8-02.2: Filter Excluding Periods -------------------------------------

#' Filter data by excluding periods before treatment
#'
#' Staggered: excludes cohort-specific periods from \eqn{g-k} through
#'   \eqn{g-1}.
#' Common Timing: excludes last k pre-treatment periods and recodes time.
#'
#' @details
#'   Paper reference: lw2025 Section 4.4, lw2026 equation 2.22
#'   Python reference: lwdid-py sensitivity.py _filter_excluding_periods()
#'   CRITICAL: Staggered mode must use
#'   \code{!is.na(data[[gvar]])} explicitly
#'   because R's NA == g returns NA (not FALSE), which would incorrectly
#'   exclude never-treated units with gvar=NA.
#'
#' @param data data.frame.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param exclude_periods integer(1). Number of periods to exclude.
#' @return data.frame. Filtered data.
#' @keywords internal
.filter_excluding_periods <- function(data, ivar, tvar, gvar, post,
                                       exclude_periods) {
  # Short-circuit: no exclusion needed
  if (exclude_periods == 0L) return(data)

  if (!is.null(gvar)) {
    # === Staggered mode ===
    cohorts <- unique(data[[gvar]])
    cohorts <- cohorts[!is.na(cohorts) & is.finite(cohorts) & cohorts > 0]

    exclude_mask <- rep(FALSE, nrow(data))
    for (g in cohorts) {
      excluded_times <- seq(as.integer(g - exclude_periods), as.integer(g - 1L))
      # CRITICAL: must check !is.na(gvar) explicitly
      # R: NA == g returns NA (not FALSE), causing NA propagation
      # Without this check, gvar=NA never-treated units get incorrectly excluded
      exclude_mask <- exclude_mask |
        (!is.na(data[[gvar]]) & data[[gvar]] == g & data[[tvar]] %in% excluded_times)
    }
    return(data[!exclude_mask, , drop = FALSE])

  } else if (!is.null(post)) {
    # === Common Timing mode ===
    pre_data <- data[data[[post]] == 0, , drop = FALSE]
    post_data <- data[data[[post]] != 0, , drop = FALSE]
    pre_times <- sort(unique(pre_data[[tvar]]))

    # Overflow check
    if (exclude_periods >= length(pre_times)) {
      warn_lwdid(
        message = sprintf(
          "Cannot exclude %d periods: only %d pre-treatment periods available. Returning original data.",
          exclude_periods, length(pre_times)),
        class = "lwdid_data"
      )
      return(data)
    }

    # Exclude last k pre-treatment periods
    excluded_times <- tail(pre_times, exclude_periods)
    filtered_pre <- pre_data[!pre_data[[tvar]] %in% excluded_times, , drop = FALSE]
    result <- rbind(filtered_pre, post_data)

    # Recode time to consecutive integers (avoid lwdid time discontinuity error)
    remaining_times <- sort(unique(result[[tvar]]))
    time_mapping <- stats::setNames(
      seq_along(remaining_times),
      as.character(remaining_times)
    )
    result[[tvar]] <- as.integer(time_mapping[as.character(result[[tvar]])])

    return(result)
  }

  data  # Fallback
}


# ---- E8-02.3: Detect Anticipation Effects -----------------------------------

#' Detect anticipation effects using dual-method sequential detection
#'
#' Method 1 (Coefficient Change): detects when ATT changes by more than
#' threshold AND magnitude increases after excluding periods.
#' Method 2 (Trend Break): detects monotonically increasing |ATT| pattern
#' with growth rate declining >= 50% at stabilization point.
#'
#' @details
#'   Python reference: lwdid-py sensitivity.py _detect_anticipation_effects()
#'   R enhancements:
#'   - Additional converged check in valid estimate filtering
#'   - New "baseline_failed" explicit status
#'   - Trend break returns excluded_periods value (not list index)
#'
#' @param estimates list. Estimation results per exclusion level.
#' @param baseline list. Baseline estimate (k=0).
#' @param detection_threshold numeric(1). Detection threshold.
#' @return list(detected, recommended_exclusion, method).
#' @keywords internal
.detect_anticipation_effects <- function(estimates, baseline, detection_threshold) {
  # Pre-check: filter valid estimates
  valid_estimates <- Filter(
    function(e) isTRUE(e$converged) && !is.na(e$att),
    estimates
  )
  if (length(valid_estimates) < 2L) {
    return(list(detected = FALSE, recommended_exclusion = 0L,
                method = "insufficient_data"))
  }

  # Pre-check: baseline validity [R enhancement]
  if (!isTRUE(baseline$converged) || is.na(baseline$att)) {
    return(list(detected = FALSE, recommended_exclusion = 0L,
                method = "baseline_failed"))
  }

  baseline_att <- baseline$att

  # === Method 1: Coefficient Change ===
  for (est in valid_estimates[-1L]) {
    if (abs(baseline_att) > 1e-10) {
      relative_change <- abs(est$att - baseline_att) / abs(baseline_att)
      if (relative_change > detection_threshold) {
        # Direction check: ATT magnitude should increase after excluding anticipation periods
        if (abs(est$att) > abs(baseline_att)) {
          return(list(
            detected = TRUE,
            recommended_exclusion = est$excluded_periods,
            method = "coefficient_change"
          ))
        }
      }
    }
  }

  # === Method 2: Trend Break ===
  if (length(valid_estimates) >= 3L) {
    abs_atts <- vapply(valid_estimates, function(e) abs(e$att), numeric(1))

    # Check monotonic non-decreasing
    is_monotone <- all(diff(abs_atts) >= 0)

    if (is_monotone) {
      # Find stabilization point: growth rate drops >= 50%
      for (i in seq(2L, length(abs_atts) - 1L)) {
        current_increase <- abs_atts[i] - abs_atts[i - 1L]
        next_increase <- abs_atts[i + 1L] - abs_atts[i]
        if (current_increase > 0 && next_increase < current_increase * 0.5) {
          # [R fix] Return excluded_periods value, not list index
          return(list(
            detected = TRUE,
            recommended_exclusion = valid_estimates[[i]]$excluded_periods,
            method = "trend_break"
          ))
        }
      }
    }
  }

  list(detected = FALSE, recommended_exclusion = 0L,
       method = "none_detected")
}


# ---- E8-02.4: Main No-Anticipation Sensitivity Function ---------------------

#' No-Anticipation Sensitivity Analysis
#'
#' Assess sensitivity of ATT estimates to potential anticipation effects
#' by systematically excluding pre-treatment periods and detecting
#' anticipation patterns using dual-method sequential detection.
#'
#' @details
#'   Paper reference:
#'   - lw2025 Section 4.4: NA assumption violation handling
#'   - lw2026 equation 2.14: NA assumption definition
#'   - lw2026 equation 2.22: truncated pre-period mean
#'   - lw2025 equation 4.4: Conditional NA assumption (CNAS)
#'
#'   Python reference: lwdid-py sensitivity.py sensitivity_no_anticipation()
#'
#' @param data data.frame or data.table. Panel data in long format.
#' @param y character(1). Outcome variable column name.
#' @param ivar character(1). Unit identifier column name.
#' @param tvar character(1). Time variable column name.
#' @param gvar character(1) or NULL. Cohort variable (staggered).
#' @param d character(1) or NULL. Treatment indicator (common timing).
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method. Default "demean".
#' @param estimator character(1). Estimation method. Default "ra".
#' @param controls character vector or NULL. Control variables.
#' @param vce character(1) or NULL. Variance estimator type.
#' @param cluster_var character(1) or NULL. Cluster variable.
#' @param max_anticipation integer(1). Maximum exclusion periods. Default 3L.
#' @param detection_threshold numeric(1). Detection threshold. Default 0.10.
#' @param alpha numeric(1). Significance level. Default 0.05.
#' @param verbose logical(1). Print progress. Default TRUE.
#' @param ... Additional arguments forwarded to underlying \code{lwdid()} calls.
#' @return An S3 object of class \code{lwdid_sensitivity} with type="no_anticipation".
#' @export
#'
#' @examples
#' \dontrun{
#' # Common Timing
#' result <- sensitivity_no_anticipation(
#'   data = panel_data, y = "outcome", ivar = "unit", tvar = "time",
#'   d = "treat", post = "post", rolling = "detrend"
#' )
#' print(result)
#' summary(result)
#' plot(result)
#' }
sensitivity_no_anticipation <- function(
  data, y, ivar, tvar,
  gvar = NULL, d = NULL, post = NULL,
  rolling = "demean", estimator = "ra",
  controls = NULL, vce = NULL, cluster_var = NULL,
  max_anticipation = 3L, detection_threshold = 0.10,
  alpha = 0.05, verbose = TRUE,
  ...
) {
  lwdid_extra_args <- .filter_sensitivity_args(list(...), lwdid)

  # Step 1: Input validation (reuse E8-01 validator)
  .validate_sensitivity_inputs(data, y, ivar, tvar, gvar, d, post, rolling)

  # Step 2: Determine maximum feasible exclusion
  min_required <- if (tolower(rolling) %in% c("detrend", "detrendq")) 2L else 1L
  max_available <- .get_max_pre_periods(data, ivar, tvar, gvar, post)
  max_feasible <- max(0L, max_available - min_required)
  max_anticipation <- min(as.integer(max_anticipation), max_feasible)

  if (max_anticipation < 1L) {
    warn_lwdid(
      message = sprintf(
        "Insufficient pre-treatment periods for anticipation analysis. Need at least %d, have %d.",
        min_required + 1L, max_available),
      class = "lwdid_data"
    )
  }

  if (verbose) {
    message("No-anticipation sensitivity analysis")
    message(sprintf("Testing exclusion range: 0 to %d", max_anticipation))
    message(strrep("-", 50))
  }

  # Step 3: Estimate at each exclusion level
  estimates <- vector("list", max_anticipation + 1L)
  result_warnings <- character(0)

  for (exclude in 0L:max_anticipation) {
    if (verbose) {
      message(sprintf("Testing exclusion = %d periods...", exclude))
    }

    tryCatch({
      filtered <- .filter_excluding_periods(
        data, ivar, tvar, gvar, post, exclude
      )

      result <- do.call(
        lwdid,
        c(
          list(
            data = filtered, y = y, ivar = ivar, tvar = tvar,
            gvar = gvar, d = d, post = post,
            rolling = rolling, estimator = estimator,
            controls = controls, vce = vce, cluster_var = cluster_var,
            alpha = alpha, verbose = "quiet"
          ),
          lwdid_extra_args
        )
      )

      n_pre_used <- max_available - exclude

      estimates[[exclude + 1L]] <- list(
        excluded_periods = exclude,
        att = result$att,
        se = result$se_att,
        t_stat = result$t_stat,
        pvalue = result$pvalue,
        ci_lower = result$ci_lower,
        ci_upper = result$ci_upper,
        n_pre_periods_used = n_pre_used,
        converged = TRUE
      )
    }, error = function(e) {
      result_warnings <<- c(
        result_warnings,
        sprintf("Exclusion of %d periods failed: %s", exclude, conditionMessage(e))
      )
      estimates[[exclude + 1L]] <<- list(
        excluded_periods = exclude,
        att = NA_real_, se = NA_real_, t_stat = NA_real_,
        pvalue = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        n_pre_periods_used = 0L,
        converged = FALSE
      )
    })
  }

  # Step 4: Identify baseline
  baseline <- estimates[[1L]]

  # Step 5: Detect anticipation effects
  detection <- .detect_anticipation_effects(
    estimates, baseline, detection_threshold
  )

  # Step 6: Generate recommendation
  if (detection$detected) {
    recommendation <- sprintf(
      "Anticipation effects detected (method: %s). Recommend excluding %d pre-treatment periods. Use lwdid(..., exclude_pre_periods = %d).",
      detection$method, detection$recommended_exclusion,
      detection$recommended_exclusion
    )
  } else {
    recommendation <- sprintf(
      "No significant anticipation effects detected (method: %s). No-anticipation assumption is reasonable.",
      detection$method
    )
  }

  # Step 7: Summary statistics
  n_converged <- sum(vapply(estimates, function(e) isTRUE(e$converged), logical(1)))

  # Step 8: Construct S3 object
  structure(list(
    type = "no_anticipation",
    estimates = estimates,
    baseline_estimate = baseline,
    anticipation_detected = detection$detected,
    recommended_exclusion = detection$recommended_exclusion,
    detection_method = detection$method,
    recommendation = recommendation,
    result_warnings = result_warnings,
    max_anticipation_tested = max_anticipation,
    detection_threshold = detection_threshold,
    rolling_method = rolling,
    estimator = estimator,
    n_estimates = length(estimates),
    n_converged = n_converged
  ), class = c("lwdid_sensitivity", "list"))
}


# ============================================================================
# E8-03: Comprehensive Sensitivity Analysis
# ============================================================================

#' Filter extra arguments for a specific function
#' @keywords internal
.filter_sensitivity_args <- function(extra_args, fun) {
  if (length(extra_args) == 0L) {
    return(list())
  }

  formal_names <- names(formals(fun))
  if (is.null(formal_names) || "..." %in% formal_names) {
    return(extra_args)
  }

  extra_args[intersect(names(extra_args), formal_names)]
}

#' Create a comprehensive sensitivity result object
#'
#' @param pre_period_result Optional pre-period sensitivity result.
#' @param no_anticipation_result Optional no-anticipation sensitivity result.
#' @param transformation_comparison Optional transformation-comparison result.
#' @param estimator_comparison Optional estimator-comparison result.
#' @param overall_assessment Optional overall assessment summary. When
#'   \code{NULL}, the default no-issue summary is used.
#' @param recommendations Optional recommendation text or vector. When
#'   \code{NULL}, the default no-issue message is used.
#' @keywords internal
.new_lwdid_sensitivity_comprehensive <- function(
  pre_period_result = NULL,
  no_anticipation_result = NULL,
  transformation_comparison = NULL,
  estimator_comparison = NULL,
  overall_assessment = NULL,
  recommendations = NULL
) {
  default_message <- "No major robustness issues found"
  if (is.null(overall_assessment)) {
    overall_assessment <- default_message
  }
  if (is.null(recommendations)) {
    recommendations <- default_message
  }
  structure(
    list(
      pre_period_result = pre_period_result,
      no_anticipation_result = no_anticipation_result,
      transformation_comparison = transformation_comparison,
      estimator_comparison = estimator_comparison,
      overall_assessment = overall_assessment,
      recommendations = recommendations
    ),
    class = c("lwdid_sensitivity_comprehensive", "list")
  )
}

#' Resolve a comprehensive result component with legacy fallbacks
#' @keywords internal
.resolve_comprehensive_component <- function(x, primary, fallback = NULL) {
  value <- x[[primary]]
  if (!is.null(value)) {
    return(value)
  }

  if (!is.null(fallback)) {
    return(x[[fallback]])
  }

  NULL
}

#' Determine whether a no-anticipation result has usable estimates
#' @keywords internal
.has_valid_no_anticipation_estimates <- function(result) {
  if (is.null(result) ||
      !inherits(result, "lwdid_sensitivity") ||
      !identical(result$type, "no_anticipation")) {
    return(FALSE)
  }

  estimates <- result$estimates
  if (!is.list(estimates) || length(estimates) == 0L) {
    return(FALSE)
  }

  any(vapply(estimates, function(e) {
    isTRUE(e$converged) &&
      is.numeric(e$att) &&
      length(e$att) == 1L &&
      !is.na(e$att)
  }, logical(1)))
}

.is_finite_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.extract_finite_lwdid_estimate <- function(
  result,
  label,
  require_se = FALSE,
  warn = TRUE
) {
  att <- result$att
  se_att <- result$se_att

  if (!.is_finite_scalar(att) ||
      (isTRUE(require_se) && !.is_finite_scalar(se_att))) {
    if (isTRUE(warn)) {
      warn_lwdid(
        message = sprintf("%s returned non-finite ATT/SE, skipping this comparison.", label),
        class = "lwdid_data"
      )
    }
    return(NULL)
  }

  list(att = att, se_att = se_att)
}

#' Summarize transformation comparison arithmetic
#' @keywords internal
.summarize_transformation_comparison <- function(
  demean_estimate,
  detrend_estimate
) {
  difference <- abs(demean_estimate$att - detrend_estimate$att)
  rel_diff <- NA_real_
  if (abs(demean_estimate$att) > 1e-10) {
    rel_diff <- difference / abs(demean_estimate$att)
  }

  list(
    demean_att = demean_estimate$att,
    demean_se = demean_estimate$se_att,
    detrend_att = detrend_estimate$att,
    detrend_se = detrend_estimate$se_att,
    difference = difference,
    rel_diff = rel_diff
  )
}

#' Summarize estimator comparison arithmetic
#' @keywords internal
.summarize_estimator_comparison <- function(results) {
  if (length(results) == 0L) {
    return(NULL)
  }

  att_names <- intersect(c("ra", "ipw", "ipwra"), names(results))
  valid_atts <- unlist(results[att_names], use.names = FALSE)
  baseline_att <- if (!is.null(results$ra) && !is.na(results$ra)) {
    results$ra
  } else {
    valid_atts[1L]
  }

  if (length(valid_atts) < 2L) {
    results$range <- NULL
    results$rel_range <- NULL
    results$baseline_att <- baseline_att
    return(results)
  }

  results$range <- max(valid_atts) - min(valid_atts)
  results$rel_range <- NA_real_
  if (abs(baseline_att) > 1e-10) {
    results$rel_range <- results$range / abs(baseline_att)
  }
  results$baseline_att <- baseline_att

  results
}

#' Run transformation comparison (demean vs detrend)
#' @keywords internal
.run_transformation_comparison <- function(
  data, y, ivar, tvar, gvar, d, post,
  estimator, controls, vce, cluster_var,
  alpha, extra_args = list()
) {
  max_pre <- .get_max_pre_periods(data, ivar, tvar, gvar, post)
  if (max_pre < 2L) {
    warn_lwdid(
      message = sprintf(
        "Insufficient pre-treatment periods for detrend comparison (need >=2, have %d)",
        max_pre
      ),
      class = "lwdid_data"
    )
    return(NULL)
  }

  lwdid_extra_args <- .filter_sensitivity_args(extra_args, lwdid)

  demean_result <- do.call(
    lwdid,
    c(
      list(
        data = data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        d = d,
        post = post,
        rolling = "demean",
        estimator = estimator,
        controls = controls,
        vce = vce,
        cluster_var = cluster_var,
        alpha = alpha,
        verbose = "quiet"
      ),
      lwdid_extra_args
    )
  )

  detrend_result <- do.call(
    lwdid,
    c(
      list(
        data = data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        d = d,
        post = post,
        rolling = "detrend",
        estimator = estimator,
        controls = controls,
        vce = vce,
        cluster_var = cluster_var,
        alpha = alpha,
        verbose = "quiet"
      ),
      lwdid_extra_args
    )
  )

  demean_estimate <- .extract_finite_lwdid_estimate(
    demean_result,
    label = "demean result",
    require_se = TRUE
  )
  detrend_estimate <- .extract_finite_lwdid_estimate(
    detrend_result,
    label = "detrend result",
    require_se = TRUE
  )

  if (is.null(demean_estimate) || is.null(detrend_estimate)) {
    return(NULL)
  }

  .summarize_transformation_comparison(
    demean_estimate = demean_estimate,
    detrend_estimate = detrend_estimate
  )
}

#' Run estimator comparison (RA, IPW, IPWRA)
#' @keywords internal
.run_estimator_comparison <- function(
  data, y, ivar, tvar, gvar, d, post,
  rolling, controls, vce, cluster_var,
  alpha, extra_args = list()
) {
  if (is.null(controls)) {
    return(NULL)
  }

  lwdid_extra_args <- .filter_sensitivity_args(extra_args, lwdid)
  results <- list()

  for (est in c("ra", "ipw", "ipwra")) {
    tryCatch({
      est_result <- do.call(
        lwdid,
        c(
          list(
            data = data,
            y = y,
            ivar = ivar,
            tvar = tvar,
            gvar = gvar,
            d = d,
            post = post,
            rolling = rolling,
            estimator = est,
            controls = controls,
            vce = vce,
            cluster_var = cluster_var,
            alpha = alpha,
            verbose = "quiet"
          ),
          lwdid_extra_args
        )
      )
      estimate <- .extract_finite_lwdid_estimate(
        est_result,
        label = sprintf("%s estimator result", est),
        require_se = TRUE,
        warn = FALSE
      )
      if (!is.null(estimate)) {
        results[[est]] <- estimate$att
        results[[paste0(est, "_se")]] <- estimate$se_att
      }
    }, error = function(e) {
      NULL
    })
  }

  if (length(results) == 0L) {
    return(NULL)
  }

  .summarize_estimator_comparison(results)
}

#' Compute the overall assessment for comprehensive sensitivity analysis
#' @keywords internal
.compute_overall_assessment <- function(
  pre_period_result,
  no_anticipation_result,
  transformation_comparison,
  estimator_comparison
) {
  issues <- character(0)
  recommendations <- character(0)

  if (!is.null(pre_period_result) && !isTRUE(pre_period_result$is_robust)) {
    issues <- c(issues, "pre-period sensitivity")
    recommendations <- c(
      recommendations,
      sprintf(
        "Pre-period sensitivity detected (SR=%.1f%%), consider using detrend or checking data quality",
        pre_period_result$sensitivity_ratio * 100
      )
    )
  }

  if (!is.null(no_anticipation_result) &&
      isTRUE(no_anticipation_result$anticipation_detected)) {
    issues <- c(issues, "anticipation effects")
    recommendations <- c(
      recommendations,
      sprintf(
        "Anticipation effects detected, recommend excluding %d pre-treatment periods",
        no_anticipation_result$recommended_exclusion
      )
    )
  }

  if (!is.null(transformation_comparison) &&
      !is.na(transformation_comparison$rel_diff) &&
      transformation_comparison$rel_diff > 0.25) {
    issues <- c(issues, "transformation sensitivity")
    recommendations <- c(
      recommendations,
      sprintf(
        "Significant difference between demean and detrend (%.1f%%), possible heterogeneous trends",
        transformation_comparison$rel_diff * 100
      )
    )
  }

  if (!is.null(estimator_comparison) &&
      !is.null(estimator_comparison$rel_range) &&
      !is.na(estimator_comparison$rel_range) &&
      estimator_comparison$rel_range > 0.25) {
    issues <- c(issues, "estimator sensitivity")
    recommendations <- c(
      recommendations,
      sprintf(
        "Significant inter-estimator difference (%.1f%%), check model specification",
        estimator_comparison$rel_range * 100
      )
    )
  }

  overall_assessment <- switch(
    as.character(length(issues)),
    "0" = "Estimates are robust across multiple robustness checks",
    "1" = sprintf("Caution: %s detected, see recommendations", issues[1L]),
    sprintf(
      "Multiple issues: %s, interpret with caution",
      paste(issues, collapse = ", ")
    )
  )

  if (length(recommendations) == 0L) {
    recommendations <- "No major robustness issues found"
  }

  list(
    overall_assessment = overall_assessment,
    recommendations = recommendations
  )
}

#' Comprehensive Sensitivity Analysis
#'
#' Unified entry point for sensitivity analysis of LWDID estimates.
#'
#' @param data data.frame. Panel data in long format.
#' @param y character(1). Outcome variable.
#' @param ivar character(1). Unit identifier.
#' @param tvar character(1). Time variable.
#' @param gvar character(1) or NULL. Cohort variable (staggered).
#' @param d character(1) or NULL. Treatment indicator.
#' @param post character(1) or NULL. Post-treatment indicator.
#' @param rolling character(1). Transformation method. Default "demean".
#' @param estimator character(1). Estimation method. Default "ra".
#' @param controls character vector or NULL. Control variables.
#' @param vce character(1) or NULL. Variance estimator type.
#' @param cluster_var character(1) or NULL. Cluster variable.
#' @param type character(1). Analysis type: "all", "pre_period",
#'   "no_anticipation", "transformation", "estimator".
#' @param max_anticipation integer(1). Maximum anticipation window. Default 3L.
#' @param detection_threshold numeric(1). Anticipation detection threshold.
#' @param alpha numeric(1). Significance level. Default 0.05.
#' @param verbose logical(1). Print progress. Default TRUE.
#' @param ... Additional arguments routed to the relevant lower-level routines.
#' @return An S3 object of class \code{lwdid_sensitivity} or
#'   \code{lwdid_sensitivity_comprehensive}.
#' @export
lwdid_sensitivity <- function(
  data, y, ivar, tvar,
  gvar = NULL, d = NULL, post = NULL,
  rolling = "demean", estimator = "ra",
  controls = NULL, vce = NULL, cluster_var = NULL,
  type = c("all", "pre_period", "no_anticipation",
           "transformation", "estimator"),
  max_anticipation = 3L,
  detection_threshold = 0.10,
  alpha = 0.05,
  verbose = TRUE,
  ...
) {
  type <- match.arg(type)
  extra_args <- list(...)

  if (type == "pre_period") {
    return(do.call(
      lwdid_sensitivity_pre_period,
      c(
        list(
          data = data,
          y = y,
          ivar = ivar,
          tvar = tvar,
          gvar = gvar,
          d = d,
          post = post,
          rolling = rolling,
          estimator = estimator,
          controls = controls,
          vce = vce,
          cluster_var = cluster_var,
          alpha = alpha,
          verbose = verbose
        ),
        .filter_sensitivity_args(extra_args, lwdid_sensitivity_pre_period)
      )
    ))
  }

  if (type == "no_anticipation") {
    return(do.call(
      sensitivity_no_anticipation,
      c(
        list(
          data = data,
          y = y,
          ivar = ivar,
          tvar = tvar,
          gvar = gvar,
          d = d,
          post = post,
          rolling = rolling,
          estimator = estimator,
          controls = controls,
          vce = vce,
          cluster_var = cluster_var,
          max_anticipation = max_anticipation,
          detection_threshold = detection_threshold,
          alpha = alpha,
          verbose = verbose
        ),
        .filter_sensitivity_args(extra_args, sensitivity_no_anticipation)
      )
    ))
  }

  .validate_sensitivity_inputs(data, y, ivar, tvar, gvar, d, post, rolling)

  pre_period_result <- NULL
  no_anticipation_result <- NULL
  transformation_comparison <- NULL
  estimator_comparison <- NULL

  if (type == "all") {
    tryCatch({
      pre_period_result <- do.call(
        lwdid_sensitivity_pre_period,
        c(
          list(
            data = data,
            y = y,
            ivar = ivar,
            tvar = tvar,
            gvar = gvar,
            d = d,
            post = post,
            rolling = rolling,
            estimator = estimator,
            controls = controls,
            vce = vce,
            cluster_var = cluster_var,
            alpha = alpha,
            verbose = FALSE
          ),
          .filter_sensitivity_args(extra_args, lwdid_sensitivity_pre_period)
        )
      )
    }, error = function(e) {
      warn_lwdid(
        message = sprintf("Pre-period robustness analysis failed: %s", conditionMessage(e)),
        class = "lwdid_data"
      )
    })

    tryCatch({
      no_anticipation_result <- do.call(
        sensitivity_no_anticipation,
        c(
          list(
            data = data,
            y = y,
            ivar = ivar,
            tvar = tvar,
            gvar = gvar,
            d = d,
            post = post,
            rolling = rolling,
            estimator = estimator,
            controls = controls,
            vce = vce,
            cluster_var = cluster_var,
            max_anticipation = max_anticipation,
            detection_threshold = detection_threshold,
            alpha = alpha,
            verbose = FALSE
          ),
          .filter_sensitivity_args(extra_args, sensitivity_no_anticipation)
        )
      )
      if (!.has_valid_no_anticipation_estimates(no_anticipation_result)) {
        no_anticipation_result <- NULL
      }
    }, error = function(e) {
      warn_lwdid(
        message = sprintf("No-anticipation analysis failed: %s", conditionMessage(e)),
        class = "lwdid_data"
      )
    })
  }

  if (type %in% c("all", "transformation")) {
    tryCatch({
      transformation_comparison <- .run_transformation_comparison(
        data = data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        d = d,
        post = post,
        estimator = estimator,
        controls = controls,
        vce = vce,
        cluster_var = cluster_var,
        alpha = alpha,
        extra_args = extra_args
      )
    }, error = function(e) {
      warn_lwdid(
        message = sprintf("Transformation comparison failed: %s", conditionMessage(e)),
        class = "lwdid_data"
      )
    })
  }

  if (type %in% c("all", "estimator")) {
    tryCatch({
      estimator_comparison <- .run_estimator_comparison(
        data = data,
        y = y,
        ivar = ivar,
        tvar = tvar,
        gvar = gvar,
        d = d,
        post = post,
        rolling = rolling,
        controls = controls,
        vce = vce,
        cluster_var = cluster_var,
        alpha = alpha,
        extra_args = extra_args
      )
    }, error = function(e) {
      warn_lwdid(
        message = sprintf("Estimator comparison failed: %s", conditionMessage(e)),
        class = "lwdid_data"
      )
    })
  }

  assessment <- .compute_overall_assessment(
    pre_period_result = pre_period_result,
    no_anticipation_result = no_anticipation_result,
    transformation_comparison = transformation_comparison,
    estimator_comparison = estimator_comparison
  )

  result <- .new_lwdid_sensitivity_comprehensive(
    pre_period_result = pre_period_result,
    no_anticipation_result = no_anticipation_result,
    transformation_comparison = transformation_comparison,
    estimator_comparison = estimator_comparison,
    overall_assessment = assessment$overall_assessment,
    recommendations = assessment$recommendations
  )

  if (isTRUE(verbose) && type == "all") {
    print(result)
  }

  result
}

#' Print method for lwdid_sensitivity_comprehensive
#'
#' @param x An lwdid_sensitivity_comprehensive object.
#' @param ... Additional arguments (ignored).
#' @return x invisibly.
#' @export
print.lwdid_sensitivity_comprehensive <- function(x, ...) {
  pre_period_result <- .resolve_comprehensive_component(
    x,
    "pre_period_result",
    "pre_period"
  )
  no_anticipation_result <- .resolve_comprehensive_component(
    x,
    "no_anticipation_result",
    "no_anticipation"
  )

  cat("\n=== Comprehensive Sensitivity Analysis ===\n\n")
  cat(sprintf("Overall: %s\n\n", x$overall_assessment))

  if (!is.null(pre_period_result)) {
    cat(sprintf(
      "  Pre-period: SR=%.1f%% (%s)\n",
      pre_period_result$sensitivity_ratio * 100,
      gsub("_", " ", pre_period_result$robustness_level)
    ))
  } else {
    cat("  Pre-period: not available\n")
  }

  if (!is.null(no_anticipation_result)) {
    cat(sprintf(
      "  Anticipation: %s\n",
      ifelse(no_anticipation_result$anticipation_detected, "detected", "not detected")
    ))
  } else {
    cat("  Anticipation: not available\n")
  }

  if (!is.null(x$transformation_comparison)) {
    tc <- x$transformation_comparison
    cat(sprintf("  Transformation: diff=%.4f", tc$difference))
    if (!is.na(tc$rel_diff)) {
      cat(sprintf(" (rel=%.1f%%)", tc$rel_diff * 100))
    }
    cat("\n")
  } else {
    cat("  Transformation: not available\n")
  }

  if (!is.null(x$estimator_comparison)) {
    ec <- x$estimator_comparison
    if (!is.null(ec$range)) {
      cat(sprintf("  Estimator: range=%.4f", ec$range))
      if (!is.null(ec$rel_range) && !is.na(ec$rel_range)) {
        cat(sprintf(" (rel=%.1f%%)", ec$rel_range * 100))
      }
      cat("\n")
    } else {
      cat("  Estimator: insufficient results\n")
    }
  } else {
    cat("  Estimator: not available\n")
  }

  cat("\nRecommendations:\n")
  for (rec in x$recommendations) {
    cat(sprintf("  - %s\n", rec))
  }

  invisible(x)
}

#' Summary method for lwdid_sensitivity_comprehensive
#'
#' @param object An lwdid_sensitivity_comprehensive object.
#' @param ... Additional arguments (ignored).
#' @return object invisibly.
#' @export
summary.lwdid_sensitivity_comprehensive <- function(object, ...) {
  pre_period_result <- .resolve_comprehensive_component(
    object,
    "pre_period_result",
    "pre_period"
  )
  no_anticipation_result <- .resolve_comprehensive_component(
    object,
    "no_anticipation_result",
    "no_anticipation"
  )

  sep <- strrep("=", 60)
  subsep <- strrep("-", 60)

  cat(sep, "\n")
  cat("Comprehensive Sensitivity Analysis Report\n")
  cat(sep, "\n\n")

  cat(subsep, "\n")
  cat("1. Pre-period Robustness\n")
  if (!is.null(pre_period_result)) {
    cat(sprintf(
      "   SR: %.4f (%.1f%%)\n",
      pre_period_result$sensitivity_ratio,
      pre_period_result$sensitivity_ratio * 100
    ))
    cat(sprintf(
      "   Level: %s\n",
      gsub("_", " ", pre_period_result$robustness_level)
    ))
    cat(sprintf(
      "   Baseline ATT: %.6f (n_pre=%d)\n",
      pre_period_result$baseline_spec$att,
      pre_period_result$baseline_spec$n_pre_periods
    ))
  } else {
    cat("   Not available\n")
  }
  cat("\n")

  cat(subsep, "\n")
  cat("2. Anticipation Effects\n")
  if (!is.null(no_anticipation_result)) {
    cat(sprintf(
      "   Detected: %s\n",
      ifelse(no_anticipation_result$anticipation_detected, "YES", "NO")
    ))
    cat(sprintf("   Method: %s\n", no_anticipation_result$detection_method))
    if (isTRUE(no_anticipation_result$anticipation_detected)) {
      cat(sprintf(
        "   Recommended exclusion: %d periods\n",
        no_anticipation_result$recommended_exclusion
      ))
    }
  } else {
    cat("   Not available\n")
  }
  cat("\n")

  cat(subsep, "\n")
  cat("3. Transformation Comparison\n")
  if (!is.null(object$transformation_comparison)) {
    tc <- object$transformation_comparison
    cat(sprintf("   Demean ATT:  %.6f\n", tc$demean_att))
    cat(sprintf("   Detrend ATT: %.6f\n", tc$detrend_att))
    cat(sprintf("   Difference:  %.6f\n", tc$difference))
    if (!is.na(tc$rel_diff)) {
      cat(sprintf("   Relative:    %.1f%%\n", tc$rel_diff * 100))
    }
  } else {
    cat("   Not available\n")
  }
  cat("\n")

  cat(subsep, "\n")
  cat("4. Estimator Comparison\n")
  if (!is.null(object$estimator_comparison)) {
    ec <- object$estimator_comparison
    for (est in c("ra", "ipw", "ipwra")) {
      if (!is.null(ec[[est]]) && !is.na(ec[[est]])) {
        cat(sprintf("   %s ATT: %.6f\n", toupper(est), ec[[est]]))
      }
    }
    if (!is.null(ec$range)) {
      cat(sprintf("   Range: %.6f\n", ec$range))
      if (!is.null(ec$rel_range) && !is.na(ec$rel_range)) {
        cat(sprintf("   Relative: %.1f%%\n", ec$rel_range * 100))
      }
    }
  } else {
    cat("   Not available (no controls)\n")
  }
  cat("\n")

  cat(sep, "\n")
  cat(sprintf("Overall Assessment: %s\n\n", object$overall_assessment))
  cat("Recommendations:\n")
  for (rec in object$recommendations) {
    cat(sprintf("  - %s\n", rec))
  }
  cat(sep, "\n")

  invisible(object)
}

#' Plot method for lwdid_sensitivity_comprehensive
#'
#' @param x An lwdid_sensitivity_comprehensive object.
#' @param show_ci Logical; whether to propagate CI drawing to child plots.
#' @param show_threshold Logical; whether to show the pre-period robustness band.
#' @param ... Additional arguments passed to sub-plots.
#' @return A ggplot object, a combined grob, or \code{NULL} invisibly.
#' @export
plot.lwdid_sensitivity_comprehensive <- function(x, show_ci = TRUE,
                                                 show_threshold = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting. Install via: install.packages('ggplot2')", call. = FALSE)
  }

  pre_period_result <- .resolve_comprehensive_component(
    x,
    "pre_period_result",
    "pre_period"
  )
  no_anticipation_result <- .resolve_comprehensive_component(
    x,
    "no_anticipation_result",
    "no_anticipation"
  )

  plots <- list()

  if (!is.null(pre_period_result)) {
    plots[[length(plots) + 1L]] <- plot(
      pre_period_result,
      show_ci = show_ci,
      show_threshold = show_threshold,
      ...
    )
  }

  if (!is.null(no_anticipation_result)) {
    plots[[length(plots) + 1L]] <- plot(
      no_anticipation_result,
      show_ci = show_ci,
      ...
    )
  }

  if (length(plots) == 0L) {
    warning("\u65e0\u53ef\u7ed8\u5236\u7ed3\u679c (no results available to plot)")
    return(invisible(NULL))
  }

  if (length(plots) == 1L) {
    return(plots[[1L]])
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plots, ncol = length(plots)))
  }

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(grobs = plots, ncol = length(plots)))
  }

  warning("Install 'gridExtra' for a combined layout. Returning a plot list.")
  plots
}
