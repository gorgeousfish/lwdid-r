# aggregate.R
# Aggregation utilities for staggered DiD.
# Implements cohort aggregated variable computation (lw2026 eq 7.9),
# unit-level gvar extraction, and NT mask identification for
# cohort/overall/event_time aggregation (lw2026 Section 7).

#' @title Extract unit-level cohort variable
#' @description Extracts each unit's cohort value from panel data,
#'   deduplicated to unit level.
#'
#'   Corresponds to Python's \code{_get_unit_level_gvar()}: Python uses
#'   \code{groupby(ivar)[gvar].first()}, which implicitly assumes \code{gvar} is
#'   time-invariant but does not validate this. R uses unique() +
#'   uniqueness validation: if the same unit has different \code{gvar}
#'   values across periods (data error), an lwdid_invalid_data
#'   error is raised immediately. This is a stricter data integrity
#'   guarantee unique to the R implementation.
#'
#' @param dt data.table
#' @param gvar character, cohort variable name
#' @param ivar character, unit identifier variable name
#' @return data.table with columns \code{ivar} and \code{gvar} (one row per
#'   unit)
#' @keywords internal
get_unit_level_gvar <- function(dt, gvar, ivar) {
  if (!data.table::is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  }
  unit_gvar <- unique(dt[, c(ivar, gvar), with = FALSE])
  dup_units <- unit_gvar[, .N, by = ivar][N > 1L]
  if (nrow(dup_units) > 0L) {
    stop_lwdid(
      sprintf(
        paste0("Unit(s) %s have non-unique cohort variable ",
               "(different gvar values across periods)"),
        paste(head(dup_units[[ivar]], 3), collapse = ", ")),
      class = "lwdid_invalid_data"
    )
  }
  unit_gvar
}

#' @title Compute cohort-level time-averaged transformed outcome
#' @description Implements lw2026 equation 7.9: for all units, compute
#'   the time-averaged transformed Y for cohort g over post-treatment
#'   periods \{g, g+1, ..., T_max\}.
#'
#'   For each unit i, computes:
#'     Y_bar_ig = (1 / n_valid_i) * sum_\{r in valid\} dot_Y_\{irg\}
#'   where dot_Y_irg is computed on-the-fly via
#'   apply_precomputed_transform(), and n_valid_i is the number of
#'   post-treatment periods with non-missing data for unit i.
#'
#' @param dt data.table, complete panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param g integer, cohort (first treatment period)
#' @param T_max integer, maximum time period
#' @param pre_stat data.table, pre-computed statistics for cohort g
#' @param rolling character, transform method ("demean" or "detrend")
#' @return data.table with columns ivar, Y_bar_ig, n_periods;
#'   NULL if no post-treatment periods
#' @keywords internal
compute_cohort_aggregated_variable <- function(dt, y, ivar, tvar,
                                                g, T_max,
                                                pre_stat, rolling) {
  if (g > T_max) {
    warn_lwdid(
      sprintf("Cohort g=%d: no post-treatment periods (g=%d > T_max=%d)",
              g, g, T_max),
      class = "lwdid_data",
      detail = "cohort_no_post_periods"
    )
    return(NULL)
  }
  post_periods <- seq.int(g, T_max)
  n_post <- length(post_periods)
  unit_ids <- pre_stat[[ivar]]
  n_units <- length(unit_ids)
  y_sum <- rep(0, n_units)
  n_valid <- rep(0L, n_units)
  for (r in post_periods) {
    period_data <- dt[get(tvar) == r]
    if (nrow(period_data) == 0L) next
    idx_in_pre <- match(period_data[[ivar]], unit_ids)
    y_trans <- apply_precomputed_transform(
      y_vals = period_data[[y]],
      unit_ids = period_data[[ivar]],
      pre_stat = pre_stat,
      rolling = rolling,
      r = r,
      ivar = ivar
    )
    for (j in seq_along(y_trans)) {
      pos <- idx_in_pre[j]
      if (!is.na(pos) && !is.na(y_trans[j])) {
        y_sum[pos] <- y_sum[pos] + y_trans[j]
        n_valid[pos] <- n_valid[pos] + 1L
      }
    }
  }
  Y_bar_ig <- ifelse(n_valid > 0L, y_sum / n_valid, NA_real_)
  result <- data.table::data.table(
    unit_id = unit_ids,
    Y_bar_ig = Y_bar_ig,
    n_periods = n_valid
  )
  data.table::setnames(result, "unit_id", ivar)
  n_missing <- sum(is.na(Y_bar_ig))
  if (n_missing > 0L) {
    warn_lwdid(
      sprintf(
        paste0("Cohort g=%d: %d/%d units have NA aggregated ",
               "variable (missing post-treatment data)"),
        g, n_missing, n_units),
      class = "lwdid_data",
      detail = "cohort_agg_missing_units"
    )
  }
  result
}

#' @title Cohort effect aggregation
#' @description Implements lw2026 equations 7.9-7.10: estimate
#'   time-averaged ATT for each cohort g.
#'
#'   Corresponds to Python aggregate_to_cohort():
#'   1. For each cohort g, compute aggregated variable Y_bar_ig
#'   2. Restrict sample to cohort g + NT (D_ig + D_i_inf = 1)
#'   3. Run OLS regression Y_bar_ig ~ 1 + D_ig (+ controls)
#'   4. Extract D_ig coefficient as tau_g
#'
#' @param dt data.table, complete panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param gvar character, cohort variable name
#' @param cohorts integer vector, all cohorts to estimate
#' @param T_max integer, maximum time period
#' @param pre_stats named list, precomputed transform statistics
#' @param rolling character, transform method
#' @param vce character or NULL, VCE type
#' @param cluster_var character or NULL, cluster variable name
#' @param alpha numeric, significance level (default 0.05)
#' @param controls character vector or NULL, control variable names
#' @return list of cohort effect results, each with 13 fields
#' @keywords internal
aggregate_to_cohort <- function(dt, y, ivar, tvar, gvar,
                                 cohorts, T_max,
                                 pre_stats, rolling,
                                 vce = NULL, cluster_var = NULL,
                                 alpha = 0.05, controls = NULL) {
  unit_gvar <- get_unit_level_gvar(dt, gvar, ivar)
  nt_mask <- is_never_treated(unit_gvar[[gvar]])
  if (sum(nt_mask) == 0L) {
    stop_lwdid(
      paste0("Cohort effect aggregation requires NT units, ",
             "but none found in data."),
      class = "lwdid_no_never_treated"
    )
  }
  nt_units <- unit_gvar[[ivar]][nt_mask]
  results <- list()
  for (g in cohorts) {
    g_char <- as.character(g)
    if (!(g_char %in% names(pre_stats))) {
      warn_lwdid(
        sprintf("Cohort g=%d: no precomputed statistics, skipping", g),
        class = "lwdid_data", detail = "cohort_agg_skipped"
      )
      next
    }
    # Step 1: Compute aggregated variable (equation 7.9)
    Y_bar_ig <- compute_cohort_aggregated_variable(
      dt, y, ivar, tvar, g, T_max, pre_stats[[g_char]], rolling
    )
    if (is.null(Y_bar_ig)) next
    # Step 2: Restrict sample (equation 7.10: D_ig + D_i_inf = 1)
    # FATAL-004: control group is NT-only
    cohort_units <- unit_gvar[[ivar]][
      unit_gvar[[gvar]] == g & !is_never_treated(unit_gvar[[gvar]])
    ]
    sample_units <- c(cohort_units, nt_units)
    reg_data <- Y_bar_ig[get(ivar) %in% sample_units]
    reg_data[, D_ig := as.integer(get(ivar) %in% cohort_units)]
    # Add control variables (unit-level, time-invariant)
    K <- 0L
    if (!is.null(controls) && length(controls) > 0L) {
      unit_controls <- unique(
        dt[, c(ivar, controls), with = FALSE]
      )
      unit_controls <- unit_controls[
        !duplicated(unit_controls[[ivar]])
      ]
      reg_data <- merge(
        reg_data, unit_controls, by = ivar, all.x = TRUE
      )
      K <- length(controls)
    }
    # Drop NA rows
    dropna_cols <- c("Y_bar_ig", controls)
    complete_mask <- complete.cases(
      reg_data[, dropna_cols, with = FALSE]
    )
    reg_data <- reg_data[complete_mask]
    # Sample size checks
    n_treat <- sum(reg_data$D_ig == 1L)
    n_ctrl <- sum(reg_data$D_ig == 0L)
    n_total <- nrow(reg_data)
    if (n_total < 2L || n_treat < 1L || n_ctrl < 1L) {
      warn_lwdid(
        sprintf(
          paste0("Cohort g=%d: insufficient sample ",
                 "(total=%d, treat=%d, ctrl=%d), skipping"),
          g, n_total, n_treat, n_ctrl),
        class = "lwdid_small_sample"
      )
      next
    }
    if (n_total == 2L) {
      warn_lwdid(
        sprintf(
          paste0("Cohort g=%d: only 2 units, ",
                 "standard errors may be unreliable"), g),
        class = "lwdid_small_sample"
      )
    }
    # Step 3: Regression estimation (equation 7.10)
    cluster_vals <- NULL
    if (!is.null(vce) && vce == "cluster") {
      if (is.null(cluster_var)) {
        stop_lwdid(
          "vce='cluster' requires cluster_var to be specified",
          class = "lwdid_invalid_parameter"
        )
      }
      unit_cluster <- unique(
        dt[, c(ivar, cluster_var), with = FALSE]
      )
      unit_cluster <- unit_cluster[
        !duplicated(unit_cluster[[ivar]])
      ]
      match_idx <- match(
        reg_data[[ivar]], unit_cluster[[ivar]]
      )
      cluster_vals <- unit_cluster[[cluster_var]][match_idx]
      n_clusters <- length(
        unique(cluster_vals[!is.na(cluster_vals)])
      )
      if (n_clusters < 20L) {
        warn_lwdid(
          sprintf(
            paste0("Cohort g=%d: few clusters (%d), ",
                   "cluster SE may be unreliable"),
            g, n_clusters),
          class = "lwdid_small_sample"
        )
      }
    }
    # 2-tier degradation for controls
    x_mat <- NULL
    if (K > 0L && n_treat > K + 1L && n_ctrl > K + 1L) {
      x_mat <- as.matrix(reg_data[, controls, with = FALSE])
    } else if (K > 0L) {
      warn_lwdid(
        sprintf(
          paste0("Cohort g=%d: controls degraded ",
                 "(N_treat=%d or N_ctrl=%d insufficient, ",
                 "need >%d)"),
          g, n_treat, n_ctrl, K + 1L),
        class = "lwdid_data",
        detail = "cohort_agg_controls_degraded"
      )
      K <- 0L
    }
    # Call RA estimator
    est <- tryCatch(
      estimate_ra_common(
        y_trans = reg_data$Y_bar_ig,
        d = reg_data$D_ig,
        x = x_mat,
        vce = vce,
        cluster = cluster_vals,
        alpha = alpha
      ),
      error = function(e) {
        warn_lwdid(
          sprintf("Cohort g=%d: regression failed - %s",
                  g, conditionMessage(e)),
          class = "lwdid_numerical"
        )
        NULL
      }
    )
    if (is.null(est)) next
    # Compute n_periods
    treat_n_periods <- Y_bar_ig[
      get(ivar) %in% cohort_units
    ]$n_periods
    n_periods <- if (length(treat_n_periods) > 0L &&
                     any(!is.na(treat_n_periods))) {
      max(treat_n_periods, na.rm = TRUE)
    } else {
      length(seq.int(g, T_max))
    }
    # Verify df fields exist
    if (is.null(est$df_resid) ||
        is.null(est$df_inference)) {
      warn_lwdid(
        sprintf(
          paste0("Cohort g=%d: estimate_ra_common() missing ",
                 "df_resid or df_inference, skipping"), g),
        class = "lwdid_numerical"
      )
      next
    }
    # Construct result (13 fields)
    results[[length(results) + 1L]] <- list(
      cohort = as.integer(g),
      att = est$att,
      se = est$se,
      ci_lower = est$ci_lower,
      ci_upper = est$ci_upper,
      t_stat = est$t_stat,
      pvalue = est$pvalue,
      n_periods = n_periods,
      n_units = n_treat,
      n_control = n_ctrl,
      df_resid = est$df_resid,
      df_inference = est$df_inference,
      K = est$K
    )
  }  # end for loop
  # Sort by cohort
  if (length(results) > 0L) {
    cohort_order <- order(
      vapply(results, `[[`, integer(1), "cohort")
    )
    results <- results[cohort_order]
  }
  # Report success/failure counts
  n_requested <- length(cohorts)
  n_success <- length(results)
  if (n_success == 0L) {
    warn_lwdid(
      sprintf(
        paste0("All %d cohort effect estimates failed.\n",
               "Possible causes:\n",
               "  1. Missing transform columns\n",
               "  2. Insufficient sample size\n",
               "  3. Aggregated variable computation failed"),
        n_requested),
      class = "lwdid_numerical"
    )
  } else if (n_success < n_requested) {
    failed_cohorts <- setdiff(
      cohorts,
      vapply(results, `[[`, integer(1), "cohort")
    )
    warn_lwdid(
      sprintf(
        paste0("Partial cohort effect failure: ",
               "%d/%d succeeded. Failed cohorts: %s"),
        n_success, n_requested,
        paste(sort(failed_cohorts), collapse = ", ")),
      class = "lwdid_numerical"
    )
  }
  results
}

#' @title Construct aggregated outcome for overall effect estimation
#' @description Implements lw2026 equation 7.18: construct the aggregated
#'   outcome variable Y_bar_i by treatment status.
#'
#'   Corresponds to Python construct_aggregated_outcome():
#'   - Treated units: use own cohort's time-averaged transform
#'   - NT units: use weighted average of all cohort transforms (FATAL-002 core)
#'   - Missing data: renormalize weights to available cohorts
#'
#'   Key difference from Python: returns list(result, diagnostics) instead
#'   of pd.Series. The diagnostics field provides richer data quality
#'   feedback (n_cohorts_requested, n_cohorts_succeeded, cohort_errors,
#'   treated_nan_count, nt_normalized_count, nt_excluded_count, yield_rate).
#'
#' @param dt data.table, complete panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param gvar character, cohort variable name
#' @param cohorts integer vector, valid cohorts
#' @param weights named numeric vector, cohort weights (omega_g = N_g / N_treat)
#' @param T_max integer, maximum time period
#' @param pre_stats named list, pre-computed transform statistics
#' @param rolling character, transform method ("demean" or "detrend")
#' @return list with result (data.table: ivar, Y_bar, D_ever) and
#'   diagnostics (n_cohorts_requested, n_cohorts_succeeded, cohort_errors,
#'   treated_nan_count, nt_normalized_count, nt_excluded_count, yield_rate)
#' @keywords internal
construct_aggregated_outcome <- function(dt, y, ivar, tvar, gvar,
                                          cohorts, weights, T_max,
                                          pre_stats, rolling) {
  # --- Input validation ---

  # Validate cohorts non-empty
  if (length(cohorts) == 0L) {
    stop_lwdid(
      paste0("cohorts cannot be empty.\n",
             "At least one treatment cohort is required.\n\n",
             "Possible causes:\n",
             "  1. No treated units found in data\n",
             "  2. All cohorts filtered out due to data issues\n",
             "  3. Incorrect gvar column specification"),
      class = "lwdid_invalid_input"
    )
  }

  # Validate weights keys match cohorts
  cohort_chars <- as.character(cohorts)
  weights_keys <- names(weights)
  if (!setequal(weights_keys, cohort_chars)) {
    missing_in_w <- setdiff(cohort_chars, weights_keys)
    extra_in_w <- setdiff(weights_keys, cohort_chars)
    stop_lwdid(
      sprintf(
        paste0("weights keys must match cohorts.\n",
               "  cohorts: %s\n  weights keys: %s\n",
               "  Missing in weights: %s\n  Extra in weights: %s"),
        paste(sort(cohorts), collapse = ", "),
        paste(sort(weights_keys), collapse = ", "),
        if (length(missing_in_w) > 0L)
          paste(sort(missing_in_w), collapse = ", ") else "None",
        if (length(extra_in_w) > 0L)
          paste(sort(extra_in_w), collapse = ", ") else "None"),
      class = "lwdid_invalid_input"
    )
  }

  # Validate weight sum (LWDID_WEIGHT_SUM_TOLERANCE = 1e-6, defined in constants.R;
  # counterpart: Python WEIGHT_SUM_TOLERANCE = 1e-6)
  weights_sum <- sum(weights)
  if (abs(weights_sum - 1.0) > LWDID_WEIGHT_SUM_TOLERANCE) {
    warn_lwdid(
      sprintf(
        paste0("Cohort weights sum to %.10f, expected 1.0. ",
               "This may indicate incorrect weight calculation."),
        weights_sum),
      class = "lwdid_numerical"
    )
  }

  # --- Get unit-level gvar and NT mask ---
  unit_gvar <- get_unit_level_gvar(dt, gvar, ivar)
  nt_mask <- is_never_treated(unit_gvar[[gvar]])

  # --- Per-cohort precomputation with tryCatch ---
  cohort_Y_bars <- list()
  cohort_errors <- list()

  for (g in cohorts) {
    g_char <- as.character(g)
    if (!(g_char %in% names(pre_stats))) {
      cohort_errors[[g_char]] <- "No pre-computed statistics"
      next
    }
    Y_bar_ig <- tryCatch(
      compute_cohort_aggregated_variable(
        dt, y, ivar, tvar, g, T_max, pre_stats[[g_char]], rolling
      ),
      error = function(e) {
        cohort_errors[[g_char]] <<- sprintf(
          "Aggregation failed: %s", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(Y_bar_ig)) {
      cohort_Y_bars[[g_char]] <- Y_bar_ig
    } else if (!(g_char %in% names(cohort_errors))) {
      cohort_errors[[g_char]] <- sprintf(
        "NULL result (cohort g=%s may have no post-treatment periods)",
        g_char)
    }
  }

  if (length(cohort_Y_bars) == 0L) {
    error_details <- paste(
      vapply(names(cohort_errors), function(g) {
        sprintf("  - cohort %s: %s", g, cohort_errors[[g]])
      }, character(1)),
      collapse = "\n"
    )
    stop_lwdid(
      sprintf(
        "All cohort aggregated variable computations failed.\n%s",
        error_details),
      class = "lwdid_insufficient_data"
    )
  }

  # --- Build wide table (each column = one cohort's Y_bar_ig) ---
  all_units <- unit_gvar[[ivar]]
  n_all <- length(all_units)
  succeeded_chars <- names(cohort_Y_bars)

  wide_dt <- data.table::data.table(uid = all_units)
  for (g_char in succeeded_chars) {
    Y_bar_ig_dt <- cohort_Y_bars[[g_char]]
    col_name <- paste0("Ybar_", g_char)
    match_idx <- match(all_units, Y_bar_ig_dt[[ivar]])
    wide_dt[, (col_name) := Y_bar_ig_dt$Y_bar_ig[match_idx]]
  }

  # --- Treated units: own cohort's transform (eq 7.18 LHS) ---
  Y_bar <- rep(NA_real_, n_all)
  D_ever <- as.integer(!nt_mask)
  treat_mask <- !nt_mask
  non_integer_units <- list()

  for (idx in which(treat_mask)) {
    g_i <- unit_gvar[[gvar]][idx]
    if (is.na(g_i) || !is.finite(g_i)) next
    g_rounded <- round(g_i)
    # Float tolerance matching for cohort assignment
    # (LWDID_COHORT_FLOAT_TOLERANCE = 1e-9, defined in constants.R;
    # counterpart: Python COHORT_FLOAT_TOLERANCE = 1e-9)
    if (abs(g_i - g_rounded) > LWDID_COHORT_FLOAT_TOLERANCE) {
      non_integer_units[[length(non_integer_units) + 1L]] <- list(
        unit = all_units[idx], gvar_val = g_i
      )
      next
    }
    g_char <- as.character(as.integer(g_rounded))
    col_name <- paste0("Ybar_", g_char)
    if (col_name %in% names(wide_dt)) {
      Y_bar[idx] <- wide_dt[[col_name]][idx]
    }
  }

  # Batch warning for non-integer gvar units
  if (length(non_integer_units) > 0L) {
    n_non_int <- length(non_integer_units)
    examples <- vapply(
      head(non_integer_units, 3L),
      function(x) sprintf("unit %s (gvar=%.4f)", x$unit, x$gvar_val),
      character(1)
    )
    example_str <- paste(examples, collapse = ", ")
    if (n_non_int > 3L) {
      example_str <- paste0(
        example_str,
        sprintf(" ... and %d more", n_non_int - 3L))
    }
    warn_lwdid(
      sprintf(
        paste0("%d units have non-integer gvar values, ",
               "skipped in aggregation.\nExamples: %s"),
        n_non_int, example_str),
      class = "lwdid_data",
      detail = "non_integer_gvar"
    )
  }

  # --- NT units: weighted average (FATAL-002, eq 7.18 RHS) ---
  # FATAL-002 protection: NT units must use the weighted average of ALL
  # cohort transforms, not a single cohort's transform. Using a single
  # cohort would bias the counterfactual baseline. Per-unit weight
  # renormalization handles partial missing data across cohorts.
  # Matrix multiply approach with per-unit weight renormalization
  nt_indices <- which(nt_mask)
  nt_normalized_count <- 0L
  nt_excluded_count <- 0L

  if (length(nt_indices) > 0L) {
    ybar_cols <- paste0("Ybar_", succeeded_chars)
    w_vec <- weights[succeeded_chars]
    nt_mat <- as.matrix(wide_dt[nt_indices, ybar_cols, with = FALSE])
    # Handle single-cohort case: ensure matrix shape
    if (!is.matrix(nt_mat)) {
      nt_mat <- matrix(nt_mat, ncol = length(succeeded_chars))
    }
    valid_mask <- !is.na(nt_mat)
    nt_mat_zero <- nt_mat
    nt_mat_zero[is.na(nt_mat_zero)] <- 0
    weighted_sums <- nt_mat_zero %*% w_vec
    valid_weight_sums <- valid_mask %*% w_vec

    for (j in seq_along(nt_indices)) {
      idx <- nt_indices[j]
      vws <- valid_weight_sums[j]
      if (vws > LWDID_WEIGHT_SUM_TOLERANCE) {
        Y_bar[idx] <- weighted_sums[j] / vws
        if (abs(vws - 1.0) > LWDID_WEIGHT_SUM_TOLERANCE) {
          nt_normalized_count <- nt_normalized_count + 1L
        }
      } else {
        nt_excluded_count <- nt_excluded_count + 1L
      }
    }
  }

  if (nt_normalized_count > 0L) {
    warn_lwdid(
      sprintf(
        paste0("%d NT units had weights renormalized due to ",
               "partial cohort data missing."),
        nt_normalized_count),
      class = "lwdid_data",
      detail = "nt_weight_renormalized"
    )
  }
  if (nt_excluded_count > 0L) {
    warn_lwdid(
      sprintf(
        paste0("%d NT units excluded due to all cohort data ",
               "missing (dropped from regression)."),
        nt_excluded_count),
      class = "lwdid_data",
      detail = "nt_excluded_all_missing"
    )
  }

  # --- Diagnostics ---
  treated_nan_count <- sum(is.na(Y_bar[treat_mask]))
  n_valid <- sum(!is.na(Y_bar))
  yield_rate <- if (n_all > 0L) n_valid / n_all else 0

  if (n_valid == 0L) {
    stop_lwdid(
      sprintf(
        paste0("All aggregated outcomes are NA.\n",
               "Successful cohorts: %s\nFailed cohorts: %s"),
        paste(succeeded_chars, collapse = ", "),
        paste(names(cohort_errors), collapse = ", ")),
      class = "lwdid_insufficient_data"
    )
  }

  if (treated_nan_count > 0L) {
    warn_lwdid(
      sprintf(
        paste0("%d treated units have NA aggregated outcome, ",
               "will be excluded from overall effect regression."),
        treated_nan_count),
      class = "lwdid_data",
      detail = "treated_nan_ybar"
    )
  }

  # --- Construct result ---
  result_dt <- data.table::data.table(
    uid_tmp = all_units,
    Y_bar = Y_bar,
    D_ever = D_ever
  )
  data.table::setnames(result_dt, "uid_tmp", ivar)

  diagnostics <- list(
    n_cohorts_requested = length(cohorts),
    n_cohorts_succeeded = length(cohort_Y_bars),
    cohort_errors = cohort_errors,
    treated_nan_count = treated_nan_count,
    nt_normalized_count = nt_normalized_count,
    nt_excluded_count = nt_excluded_count,
    yield_rate = yield_rate
  )

  list(result = result_dt, diagnostics = diagnostics)
}



#' @title Overall weighted effect estimation
#' @description Implements lw2026 equations 7.18-7.19 to estimate the overall
#'   weighted ATT tau_omega via OLS regression. The workflow computes cohort
#'   weights omega_g = N_g / N_treat (eq. 7.12), constructs the aggregated
#'   outcome with `construct_aggregated_outcome()` (eq. 7.18, including
#'   FATAL-002 protection), runs `estimate_ra_common()` on `Y_bar ~ D_ever`
#'   with optional controls (eq. 7.19), and extracts tau_omega as the
#'   coefficient on `D_ever`.
#'
#'   This function corresponds to Python `aggregate_to_overall()`, but the R
#'   implementation keeps both `cohort_weights` and `effective_weights`, uses an
#'   explicit three-tier controls degradation path through
#'   `estimate_ra_common()`, and stores `df_resid`, `df_inference`, `n_total`,
#'   `K`, `controls_tier`, `controls_used`, and `diagnostics` as explicit result
#'   fields.
#'
#' @param dt data.table, complete panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param gvar character, cohort variable name
#' @param cohorts integer vector, all cohorts to estimate
#' @param T_max integer, maximum time period
#' @param pre_stats named list, precomputed transform statistics
#' @param rolling character, transform method ("demean" or "detrend")
#' @param vce character or NULL, VCE type ("homoskedastic", "HC1", "cluster")
#' @param cluster_var character or NULL, cluster variable name
#'   (required when vce="cluster")
#' @param alpha numeric, significance level (default 0.05)
#' @param controls character vector or NULL, control variable names
#'   (optional extension beyond paper eq 7.19)
#' @return A named list containing the overall ATT estimate (`att`), its
#'   standard error (`se`), confidence interval bounds (`ci_lower`,
#'   `ci_upper`), test statistics (`t_stat`, `pvalue`), original and effective
#'   cohort weights (`cohort_weights`, `effective_weights`), sample-size and
#'   degrees-of-freedom diagnostics (`n_treated`, `n_control`, `n_total`,
#'   `df_resid`, `df_inference`, `K`), controls metadata (`controls_tier`,
#'   `controls_used`), and the aggregated-outcome diagnostics returned by
#'   `construct_aggregated_outcome()`.
#' @keywords internal
aggregate_to_overall <- function(dt, y, ivar, tvar, gvar,
                                  cohorts, T_max,
                                  pre_stats, rolling,
                                  vce = NULL, cluster_var = NULL,
                                  alpha = 0.05, controls = NULL) {
  # --- Data preparation ---
  unit_gvar <- get_unit_level_gvar(dt, gvar, ivar)
  nt_mask <- is_never_treated(unit_gvar[[gvar]])

  if (sum(nt_mask) == 0L) {
    stop_lwdid(
      paste0("Overall effect estimation requires never-treated (NT) units, ",
             "but none found in data.\n",
             "Different cohorts use different pre-treatment transforms; ",
             "only NT units provide consistent counterfactual baselines."),
      class = "lwdid_no_never_treated"
    )
  }

  # --- Step 1: Compute cohort weights omega_g = N_g / N_treat (eq 7.12) ---
  # Counterpart: Python cohort_sizes with get_cohort_mask using
  # COHORT_FLOAT_TOLERANCE = 1e-9
  cohort_sizes <- vapply(cohorts, function(g) {
    gvar_vals <- unit_gvar[[gvar]]
    is_nt <- is_never_treated(gvar_vals)
    sum(abs(gvar_vals - g) < LWDID_COHORT_FLOAT_TOLERANCE & !is_nt,
        na.rm = TRUE)
  }, integer(1))
  names(cohort_sizes) <- as.character(cohorts)

  # Exclude zero-size cohorts (counterpart: Python zero_cohorts filtering)
  zero_cohorts <- cohorts[cohort_sizes == 0L]
  if (length(zero_cohorts) > 0L) {
    warn_lwdid(
      sprintf("Zero-size cohorts excluded: %s. May indicate data filtering issues.",
              paste(sort(zero_cohorts), collapse = ", ")),
      class = "lwdid_data"
    )
  }

  valid_cohorts <- cohorts[cohort_sizes > 0L]
  if (length(valid_cohorts) == 0L) {
    stop_lwdid(
      "No valid cohorts after filtering (all cohort sizes are 0). Check data integrity.",
      class = "lwdid_insufficient_data"
    )
  }
  valid_sizes <- cohort_sizes[as.character(valid_cohorts)]
  N_treat <- sum(as.numeric(valid_sizes))

  if (N_treat == 0) {
    stop_lwdid("No treated units", class = "lwdid_insufficient_data")
  }

  # Computation weights (used for aggregated outcome construction)
  computation_weights <- valid_sizes / N_treat
  names(computation_weights) <- as.character(valid_cohorts)

  # Weight sum validation (counterpart: Python WEIGHT_SUM_TOLERANCE)
  w_sum <- sum(computation_weights)
  if (abs(w_sum - 1.0) > LWDID_WEIGHT_SUM_TOLERANCE) {
    warn_lwdid(
      sprintf("Cohort weight sum is %.10f, expected 1.0", w_sum),
      class = "lwdid_numerical"
    )
    computation_weights <- computation_weights / w_sum
  }

  # --- Step 2: Construct aggregated outcome (eq 7.18, FATAL-002 core) ---
  agg_result <- construct_aggregated_outcome(
    dt, y, ivar, tvar, gvar,
    weights = computation_weights,
    cohorts = valid_cohorts,
    T_max = T_max,
    pre_stats = pre_stats,
    rolling = rolling
  )
  agg_data <- agg_result$result
  agg_diagnostics <- agg_result$diagnostics

  # --- Step 3: Construct regression data (eq 7.19) ---
  reg_data <- data.frame(
    Y_bar = agg_data$Y_bar,
    D_ever = agg_data$D_ever,
    stringsAsFactors = FALSE
  )
  all_units <- agg_data[[ivar]]

  # Controls extension (optional, not part of eq 7.19)
  # 3-tier degradation delegated to estimate_ra_common()
  # [R enhancement] Difference from Story E5-02: E5-02 pre-filters at
  # caller level (2-tier: Tier 1 -> Tier 3); E5-03 passes controls
  # unconditionally to estimate_ra_common() for internal tier selection.
  K <- 0L
  controls_used <- FALSE
  if (!is.null(controls) && length(controls) > 0L) {
    warn_lwdid(
      paste0("Overall effect regression includes controls. ",
             "Note: eq 7.19 regression equivalence (7.13-7.17 derivation) ",
             "assumes no controls.\n",
             "With controls, tau_omega is a conditional weighted effect, ",
             "not strictly equal to sum(omega_g * tau_g^ctrl)."),
      class = "lwdid_data", detail = "overall_controls_extension"
    )
    unit_controls <- unique(dt[, c(ivar, controls), with = FALSE])
    unit_controls <- unit_controls[!duplicated(unit_controls[[ivar]])]
    for (ctrl in controls) {
      match_idx <- match(all_units, unit_controls[[ivar]])
      reg_data[[ctrl]] <- unit_controls[[ctrl]][match_idx]
    }
    K <- length(controls)
    controls_used <- TRUE
  }

  # Drop NA rows (Y_bar NA + control NA)
  dropna_cols <- "Y_bar"
  if (controls_used) dropna_cols <- c(dropna_cols, controls)
  complete_mask <- complete.cases(reg_data[, dropna_cols, drop = FALSE])
  reg_data <- reg_data[complete_mask, , drop = FALSE]
  units_in_reg <- all_units[complete_mask]

  n_treated <- sum(reg_data$D_ever == 1L)
  n_control <- sum(reg_data$D_ever == 0L)
  n_total <- nrow(reg_data)

  # --- Effective weights check ---
  # [R enhancement] Preserve both computation and effective weights
  # (Python updates cohort_weights to effective_weights when deviation > 1%)
  effective_gvar <- unit_gvar[match(units_in_reg, unit_gvar[[ivar]]), ]
  effective_sizes <- vapply(valid_cohorts, function(g) {
    gvar_vals <- effective_gvar[[gvar]]
    is_nt <- is_never_treated(gvar_vals)
    sum(abs(gvar_vals - g) < LWDID_COHORT_FLOAT_TOLERANCE & !is_nt,
        na.rm = TRUE)
  }, integer(1))
  names(effective_sizes) <- as.character(valid_cohorts)
  N_treat_eff <- sum(as.numeric(effective_sizes))

  effective_weights <- computation_weights  # default: same as computation
  if (N_treat_eff > 0) {
    effective_weights <- effective_sizes / N_treat_eff
    names(effective_weights) <- as.character(valid_cohorts)
    max_diff <- max(abs(
      computation_weights - effective_weights[names(computation_weights)]
    ))
    if (max_diff > 0.01) {
      warn_lwdid(
        sprintf(
          paste0("Post-dropna cohort weight deviation > 1%% ",
                 "(max deviation: %.3f).\n",
                 "May be caused by differential NA rates across cohorts.\n",
                 "Computation weights: %s\n",
                 "Effective weights: %s"),
          max_diff,
          paste(names(computation_weights), "=",
                round(computation_weights, 4), collapse = ", "),
          paste(names(effective_weights), "=",
                round(effective_weights, 4), collapse = ", ")),
        class = "lwdid_numerical"
      )
    }
  }

  # Sample size checks
  if (n_total < 2L || n_treated < 1L || n_control < 1L) {
    stop_lwdid(
      sprintf(
        paste0("Overall effect regression sample size insufficient: ",
               "total=%d, treated=%d, control=%d"),
        n_total, n_treated, n_control),
      class = "lwdid_insufficient_data"
    )
  }

  if (n_total == 2L) {
    warn_lwdid(
      "Sample size is only 2 units, standard errors may be unreliable",
      class = "lwdid_small_sample"
    )
  }

  # --- VCE validation and cluster variable extraction ---
  # Counterpart: Python raises ValueError when vce=='cluster' and
  # cluster_var is None
  if (!is.null(vce) && vce == "cluster" && is.null(cluster_var)) {
    stop_lwdid(
      paste0("vce='cluster' requires cluster_var parameter.\n",
             "Please specify the clustering variable name."),
      class = "lwdid_invalid_input"
    )
  }

  cluster_vals <- NULL
  if (!is.null(cluster_var) && !is.null(vce) && vce == "cluster") {
    unit_cluster <- unique(dt[, c(ivar, cluster_var), with = FALSE])
    unit_cluster <- unit_cluster[!duplicated(unit_cluster[[ivar]])]
    cluster_match <- match(units_in_reg, unit_cluster[[ivar]])
    cluster_vals <- unit_cluster[[cluster_var]][cluster_match]

    n_missing_cluster <- sum(is.na(cluster_vals))
    if (n_missing_cluster > 0L) {
      warn_lwdid(
        sprintf("%d units missing cluster variable values, may cause cluster SE computation errors",
                n_missing_cluster),
        class = "lwdid_data"
      )
    }

    n_clusters <- length(unique(cluster_vals[!is.na(cluster_vals)]))
    if (n_clusters < 20L) {
      warn_lwdid(
        sprintf(
          paste0("Number of clusters (%d) is small, cluster standard ",
                 "errors may be unreliable (recommend >= 20)"),
          n_clusters),
        class = "lwdid_small_sample"
      )
    }
  }

  # Prepare control variable matrix
  # 3-tier degradation (Tier 1 full interaction / Tier 2 simple controls /
  # Tier 3 no controls) delegated to estimate_ra_common().
  x_mat <- NULL
  if (K > 0L) {
    x_mat <- as.matrix(reg_data[, controls, drop = FALSE])
  }

  # --- Step 4: Regression estimation (eq 7.19) ---
  est <- estimate_ra_common(
    y_trans = reg_data$Y_bar,
    d = reg_data$D_ever,
    x = x_mat,
    vce = vce,
    cluster = cluster_vals,
    alpha = alpha
  )

  # Validate estimate_ra_common() returned df_resid and df_inference
  # (counterpart: Python aggregate_to_overall() explicit validation)
  if (is.null(est$df_resid) || is.null(est$df_inference)) {
    stop_lwdid(
      sprintf(
        paste0("Overall effect regression result missing df_resid or ",
               "df_inference field.\nAvailable fields: %s"),
        paste(names(est), collapse = ", ")),
      class = "lwdid_numerical"
    )
  }

  # --- Return 17-field result list ---
  list(
    att = est$att,
    se = est$se,
    ci_lower = est$ci_lower,
    ci_upper = est$ci_upper,
    t_stat = est$t_stat,
    pvalue = est$pvalue,
    # Computation weights: passed to construct_aggregated_outcome
    # (based on original cohort sizes before dropna)
    cohort_weights = as.list(computation_weights),
    # [R enhancement] Effective weights: post-dropna treated unit
    # cohort proportions (for diagnostics).
    # Python OverallEffect has no effective_weights field; Python
    # updates cohort_weights to effective_weights when deviation > 1%.
    effective_weights = as.list(effective_weights),
    n_treated = n_treated,
    n_control = n_control,
    # [R enhancement] n_total: Python OverallEffect has no this field
    n_total = n_total,
    # df_resid: OLS residual df (n - k), independent of VCE type
    # df_inference: inference df (homoskedastic/HC: n-k; cluster: G-1)
    df_resid = est$df_resid,
    df_inference = est$df_inference,
    K = est$K,
    # [R enhancement] controls_tier/controls_used
    controls_tier = if (!is.null(est$controls_tier)) est$controls_tier else "none",
    controls_used = controls_used,
    # [R enhancement] diagnostics from construct_aggregated_outcome
    diagnostics = agg_diagnostics
  )
}

# ============================================================================
# Event-Time Aggregation Helper Functions
# ============================================================================

#' @title Compute event-time weights for WATT aggregation
#' @description Compute normalized weights for available cohorts at a given
#'   event time. Weights are proportional to cohort sizes (N_g).
#'   Counterpart: Python compute_event_time_weights().
#'
#' @param cohort_sizes named integer vector, cohort sizes keyed by cohort
#'   (as character). E.g., c("3" = 30, "5" = 20).
#' @param available_cohorts character vector, cohorts available at this
#'   event time.
#' @return named numeric vector, normalized weights summing to 1.0.
#' @keywords internal
compute_event_time_weights <- function(cohort_sizes, available_cohorts) {
  # Empty available_cohorts check
  if (length(available_cohorts) == 0L) {
    stop_lwdid(
      "No available cohorts for event-time weight computation.",
      class = "lwdid_invalid_input"
    )
  }

  # Extract sizes for available cohorts (missing keys default to 0)
  sizes <- vapply(available_cohorts, function(g) {
    g_char <- as.character(g)
    if (g_char %in% names(cohort_sizes)) {
      as.integer(cohort_sizes[[g_char]])
    } else {
      0L
    }
  }, integer(1))
  names(sizes) <- as.character(available_cohorts)

  total_size <- sum(as.numeric(sizes))
  if (total_size <= 0) {
    stop_lwdid(
      sprintf("Total cohort size is %d for available cohorts [%s]. Cannot compute weights.",
              as.integer(total_size), paste(available_cohorts, collapse = ", ")),
      class = "lwdid_insufficient_data"
    )
  }

  # Normalize
  weights <- sizes / total_size
  names(weights) <- as.character(available_cohorts)
  weights
}


#' @title Validate that weights sum to 1.0 within tolerance
#' @description Check whether a weight vector sums to 1.0 within
#'   LWDID_WEIGHT_SUM_TOLERANCE. Returns validation result.
#'   Counterpart: Python validate_weight_sum().
#'
#' @param weights numeric vector, weights to validate.
#' @param tolerance numeric, absolute tolerance (default LWDID_WEIGHT_SUM_TOLERANCE).
#' @return list with is_valid (logical) and weight_sum (numeric).
#' @keywords internal
validate_weight_sum <- function(weights, tolerance = LWDID_WEIGHT_SUM_TOLERANCE) {
  weight_sum <- sum(weights)
  is_valid <- abs(weight_sum - 1.0) <= tolerance
  list(is_valid = is_valid, weight_sum = weight_sum)
}


#' @title Select degrees of freedom for event-time inference
#' @description Select inference degrees of freedom from per-cohort df values
#'   using one of three strategies: conservative (min), weighted (weighted
#'   average), or fallback (n_cohorts - 1).
#'   Counterpart: Python select_degrees_of_freedom().
#'
#' @param cohort_dfs numeric vector, per-cohort df_inference values (may contain NA).
#' @param weights numeric vector, cohort weights (same length as cohort_dfs).
#' @param strategy character, one of "conservative", "weighted", "fallback".
#' @param n_cohorts integer, number of cohorts (used for fallback).
#' @return integer, selected degrees of freedom (at least 1).
#' @keywords internal
select_degrees_of_freedom <- function(cohort_dfs, weights, strategy, n_cohorts) {
  # Validate strategy
  valid_strategies <- c("conservative", "weighted", "fallback")
  if (!strategy %in% valid_strategies) {
    stop_lwdid(
      sprintf("Invalid df_strategy '%s'. Must be one of: %s",
              strategy, paste(valid_strategies, collapse = ", ")),
      class = "lwdid_invalid_input"
    )
  }

  # Filter valid dfs (non-NA, positive)
  valid_mask <- !is.na(cohort_dfs) & cohort_dfs > 0
  valid_dfs <- cohort_dfs[valid_mask]
  valid_weights <- weights[valid_mask]

  # No valid dfs -> auto fallback
  if (length(valid_dfs) == 0L) {
    return(max(1L, as.integer(n_cohorts - 1L)))
  }

  df <- switch(strategy,
    "conservative" = max(1L, as.integer(min(valid_dfs))),
    "weighted" = {
      # Renormalize weights to valid cohorts
      w_sum <- sum(valid_weights)
      if (w_sum > 0) {
        w_norm <- valid_weights / w_sum
        max(1L, as.integer(round(sum(w_norm * valid_dfs))))
      } else {
        max(1L, as.integer(n_cohorts - 1L))
      }
    },
    "fallback" = max(1L, as.integer(n_cohorts - 1L))
  )

  df
}


# ============================================================================
# Event-Time Aggregation Main Function
# ============================================================================

#' @title Aggregate cohort-period effects to event-time WATT
#' @description Aggregate (g,r) treatment effect estimates to event-time
#'   weighted average treatment effects WATT(e). Implements lw2026
#'   Section 7 event-time aggregation.
#'
#'   Counterpart: Python aggregate_to_event_time().
#'
#'   Key differences from Python:
#'   - Supports both data.frame and list input for cohort_time_effects
#'   - R-side extension: pre-validates \code{alpha}, \code{df_strategy}, and
#'     \code{event_time_range} before the main loop
#'   - R-side extension: detects and deduplicates \code{(cohort, period)} pairs
#'   - R-side extension: uses \code{qnorm(1 - alpha / 2)} when
#'     \code{df <= 0} instead of a hard-coded 1.96
#'   - R-side extension: verbose summary uses \code{message()} rather than
#'     \code{warning()}
#'
#' @param cohort_time_effects data.frame or list of lists. Required
#'   columns/fields: cohort (integer), period (integer), att (numeric),
#'   se (numeric). Optional: df_inference (integer).
#' @param cohort_sizes named integer vector or named list, cohort sizes
#'   (N_g). Names are cohort values as character.
#' @param alpha numeric in (0,1), significance level. Default 0.05.
#' @param event_time_range integer vector of length 2 or NULL,
#'   \code{c(min_e, max_e)} filter range. Default NULL (no filtering).
#' @param df_strategy character, one of "conservative" (default), "weighted",
#'   "fallback".
#' @param verbose logical, whether to output diagnostic messages. Default FALSE.
#' @return list of lists, each with: event_time, att, se, ci_lower, ci_upper,
#'   t_stat, pvalue, df_inference, n_cohorts, cohort_contributions,
#'   weight_sum, alpha.
#' @keywords internal
aggregate_to_event_time <- function(cohort_time_effects,
                                    cohort_sizes,
                                    alpha = 0.05,
                                    event_time_range = NULL,
                                    df_strategy = "conservative",
                                    verbose = FALSE) {

  # ========================================================================
  # Phase 1: Input validation
  # ========================================================================

  # --- Pre-validate scalar parameters (fail fast, R enhancement) ---
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(
      sprintf("alpha must be a single numeric value in (0, 1), got: %s",
              deparse(alpha)),
      class = "lwdid_invalid_input"
    )
  }

  valid_df_strategies <- c("conservative", "weighted", "fallback")
  if (!is.character(df_strategy) || length(df_strategy) != 1L ||
      !(df_strategy %in% valid_df_strategies)) {
    stop_lwdid(
      sprintf("df_strategy must be one of: %s. Got: '%s'",
              paste(valid_df_strategies, collapse = ", "),
              as.character(df_strategy)),
      class = "lwdid_invalid_input"
    )
  }

  if (!is.null(event_time_range)) {
    if (!is.numeric(event_time_range) || length(event_time_range) != 2L) {
      stop_lwdid(
        sprintf(
          "event_time_range must be NULL or numeric vector of length 2, got length %d",
          length(event_time_range)),
        class = "lwdid_invalid_input"
      )
    }
    if (any(is.na(event_time_range))) {
      stop_lwdid(
        "event_time_range must not contain NA values.",
        class = "lwdid_invalid_input"
      )
    }
  }

  # --- Convert list input to data.frame ---
  if (is.list(cohort_time_effects) && !is.data.frame(cohort_time_effects)) {
    # List of lists: each element should have cohort, period, att, se
    if (length(cohort_time_effects) == 0L) {
      stop_lwdid(
        "cohort_time_effects is an empty list.",
        class = "lwdid_invalid_input"
      )
    }
    cohort_time_effects <- do.call(rbind, lapply(cohort_time_effects, function(x) {
      data.frame(
        cohort = as.integer(x$cohort),
        period = as.integer(x$period),
        att = as.numeric(x$att),
        se = as.numeric(x$se),
        df_inference = if (!is.null(x$df_inference)) as.integer(x$df_inference)
                       else NA_integer_,
        stringsAsFactors = FALSE
      )
    }))
  }

  if (!is.data.frame(cohort_time_effects)) {
    stop_lwdid(
      "cohort_time_effects must be a data.frame or list of lists.",
      class = "lwdid_invalid_input"
    )
  }

  # --- Validate required columns ---
  required_cols <- c("cohort", "period", "att", "se")
  missing_cols <- setdiff(required_cols, names(cohort_time_effects))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      sprintf("cohort_time_effects missing required columns: %s",
              paste(missing_cols, collapse = ", ")),
      class = "lwdid_invalid_input"
    )
  }

  # Handle empty data.frame
  if (nrow(cohort_time_effects) == 0L) {
    return(list())
  }

  # Ensure df_inference column exists
  if (!"df_inference" %in% names(cohort_time_effects)) {
    cohort_time_effects$df_inference <- NA_integer_
  }

  # --- Validate cohort_sizes ---
  # Convert named list to named numeric vector
  if (is.list(cohort_sizes) && !is.null(names(cohort_sizes))) {
    cohort_sizes <- vapply(cohort_sizes, function(x) as.numeric(x), numeric(1))
  }
  if (!is.numeric(cohort_sizes) && !is.integer(cohort_sizes)) {
    stop_lwdid(
      "cohort_sizes must be a named numeric/integer vector or named list.",
      class = "lwdid_invalid_input"
    )
  }
  if (length(cohort_sizes) == 0L) {
    stop_lwdid(
      "cohort_sizes must be non-empty.",
      class = "lwdid_invalid_input"
    )
  }
  if (is.null(names(cohort_sizes)) ||
      any(is.na(names(cohort_sizes))) ||
      any(names(cohort_sizes) == "")) {
    stop_lwdid(
      "cohort_sizes must have valid (non-NA, non-empty) names.",
      class = "lwdid_invalid_input"
    )
  }
  if (any(is.na(cohort_sizes))) {
    stop_lwdid(
      "cohort_sizes must not contain NA values.",
      class = "lwdid_invalid_input"
    )
  }
  if (any(cohort_sizes <= 0)) {
    stop_lwdid(
      sprintf("cohort_sizes must all be positive. Found non-positive values: %s",
              paste(cohort_sizes[cohort_sizes <= 0], collapse = ", ")),
      class = "lwdid_invalid_input"
    )
  }
  # Ensure numeric type for arithmetic
  cohort_sizes_vec <- as.numeric(cohort_sizes)
  names(cohort_sizes_vec) <- names(cohort_sizes)

  # ========================================================================
  # Phase 2: Data preparation
  # ========================================================================

  # --- Compute event_time ---
  cohort_time_effects$event_time <- cohort_time_effects$period -
                                    cohort_time_effects$cohort

  # --- Filter by event_time_range if specified ---
  if (!is.null(event_time_range)) {
    min_e <- event_time_range[1L]
    max_e <- event_time_range[2L]
    keep_mask <- cohort_time_effects$event_time >= min_e &
                 cohort_time_effects$event_time <= max_e
    n_before <- nrow(cohort_time_effects)
    cohort_time_effects <- cohort_time_effects[keep_mask, , drop = FALSE]
    n_after <- nrow(cohort_time_effects)
    if (verbose && n_before != n_after) {
      message(sprintf(
        "[aggregate_to_event_time] Filtered event_time to [%d, %d]: %d -> %d rows",
        min_e, max_e, n_before, n_after))
    }
    if (n_after == 0L) {
      if (verbose) {
        message("[aggregate_to_event_time] No rows remain after event_time_range filter.")
      }
      return(list())
    }
  }

  # --- Detect and deduplicate (cohort, period) pairs (R enhancement) ---
  dup_key <- paste(cohort_time_effects$cohort,
                   cohort_time_effects$period, sep = "_")
  dup_mask <- duplicated(dup_key)
  if (any(dup_mask)) {
    n_dups <- sum(dup_mask)
    dup_pairs <- unique(dup_key[dup_mask])
    warn_lwdid(
      sprintf(
        paste0("Detected %d duplicate (cohort, period) rows. ",
               "Keeping first occurrence. Duplicated pairs: %s"),
        n_dups,
        paste(head(dup_pairs, 5L), collapse = ", ")),
      class = "lwdid_data"
    )
    cohort_time_effects <- cohort_time_effects[!dup_mask, , drop = FALSE]
  }

  # ========================================================================
  # Phase 3: Per-event-time aggregation loop
  # ========================================================================

  unique_event_times <- sort(unique(cohort_time_effects$event_time))
  results <- vector("list", length(unique_event_times))
  n_nan_results <- 0L

  for (i in seq_along(unique_event_times)) {
    e <- unique_event_times[i]

    # --- 9a: Extract rows for this event time ---
    rows_e <- cohort_time_effects[cohort_time_effects$event_time == e, ,
                                  drop = FALSE]

    # --- 9b: Filter out rows with NA/NaN/Inf att or se ---
    valid_att <- is.finite(rows_e$att)
    valid_se <- is.finite(rows_e$se)
    valid_rows <- valid_att & valid_se
    n_invalid <- sum(!valid_rows)

    if (n_invalid > 0L && verbose) {
      message(sprintf(
        "[aggregate_to_event_time] event_time=%d: %d/%d rows have invalid att/se, excluded.",
        e, n_invalid, nrow(rows_e)))
    }

    rows_valid <- rows_e[valid_rows, , drop = FALSE]

    # --- 9c: If all invalid -> store NaN result ---
    if (nrow(rows_valid) == 0L) {
      n_nan_results <- n_nan_results + 1L
      results[[i]] <- list(
        event_time = as.integer(e),
        att = NaN,
        se = NaN,
        ci_lower = NaN,
        ci_upper = NaN,
        t_stat = NaN,
        pvalue = NaN,
        df_inference = NA_integer_,
        n_cohorts = 0L,
        cohort_contributions = list(),
        weight_sum = NaN,
        alpha = alpha
      )
      next
    }

    # --- 9d: Get available cohorts (as character) ---
    available_cohorts <- as.character(rows_valid$cohort)
    n_cohorts <- length(available_cohorts)

    # --- 9e: Compute event-time weights ---
    weights <- tryCatch(
      compute_event_time_weights(cohort_sizes_vec, available_cohorts),
      error = function(err) {
        if (verbose) {
          message(sprintf(
            "[aggregate_to_event_time] event_time=%d: weight computation failed: %s",
            e, conditionMessage(err)))
        }
        NULL
      }
    )

    if (is.null(weights)) {
      n_nan_results <- n_nan_results + 1L
      results[[i]] <- list(
        event_time = as.integer(e),
        att = NaN,
        se = NaN,
        ci_lower = NaN,
        ci_upper = NaN,
        t_stat = NaN,
        pvalue = NaN,
        df_inference = NA_integer_,
        n_cohorts = n_cohorts,
        cohort_contributions = list(),
        weight_sum = NaN,
        alpha = alpha
      )
      next
    }

    # --- 9f: Validate weight sum ---
    wv <- validate_weight_sum(weights)
    if (!wv$is_valid && verbose) {
      message(sprintf(
        "[aggregate_to_event_time] event_time=%d: weight sum = %.10f (expected 1.0)",
        e, wv$weight_sum))
    }

    # --- 9g: Compute WATT(e) = sum(w_g * att_g) ---
    att_vec <- rows_valid$att
    se_vec <- rows_valid$se
    w_vec <- weights[available_cohorts]
    watt <- sum(w_vec * att_vec)

    # --- 9h: Compute SE(WATT(e)) = sqrt(sum(w_g^2 * se_g^2)) ---
    se_watt <- sqrt(sum(w_vec^2 * se_vec^2))

    # --- 9i: Collect per-cohort df_inference values ---
    cohort_dfs <- rows_valid$df_inference

    # --- 9j: Select degrees of freedom ---
    df <- select_degrees_of_freedom(cohort_dfs, w_vec, df_strategy, n_cohorts)

    # Verbose: check if all df are NA (R enhancement)
    if (verbose && all(is.na(cohort_dfs))) {
      message(sprintf(
        "[aggregate_to_event_time] event_time=%d: all cohort df_inference are NA, using fallback df=%d",
        e, df))
    }

    # --- 9k: Compute t_stat (handle se_watt == 0) ---
    if (se_watt == 0) {
      if (watt > 0) {
        t_stat <- Inf
      } else if (watt < 0) {
        t_stat <- -Inf
      } else {
        t_stat <- NaN
      }
    } else {
      t_stat <- watt / se_watt
    }

    # --- 9l: Compute p-value ---
    if (is.finite(t_stat)) {
      pvalue <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
    } else if (is.infinite(t_stat)) {
      pvalue <- 0
    } else {
      # NaN t_stat
      pvalue <- NaN
    }

    # --- 9m-o: Compute CI ---
    if (df > 0) {
      t_crit <- stats::qt(1 - alpha / 2, df = df)
    } else {
      # R improvement: use qnorm instead of hardcoded 1.96
      t_crit <- stats::qnorm(1 - alpha / 2)
    }
    ci_lower <- watt - t_crit * se_watt
    ci_upper <- watt + t_crit * se_watt

    # --- 9p: Build cohort_contributions ---
    cohort_contributions <- vector("list", n_cohorts)
    for (j in seq_len(n_cohorts)) {
      g_char <- available_cohorts[j]
      cohort_contributions[[j]] <- list(
        cohort = as.integer(rows_valid$cohort[j]),
        weight = as.numeric(w_vec[g_char]),
        att = as.numeric(att_vec[j]),
        se = as.numeric(se_vec[j])
      )
    }

    # --- 9q: Store result ---
    results[[i]] <- list(
      event_time = as.integer(e),
      att = watt,
      se = se_watt,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      t_stat = t_stat,
      pvalue = pvalue,
      df_inference = as.integer(df),
      n_cohorts = as.integer(n_cohorts),
      cohort_contributions = cohort_contributions,
      weight_sum = wv$weight_sum,
      alpha = alpha
    )
  }  # end for loop over event times

  # ========================================================================
  # Phase 4: Sort and summary
  # ========================================================================

  # Results are already sorted since we iterated over sorted unique_event_times

  # --- Verbose summary via message() ---
  if (verbose) {
    n_total <- length(results)
    n_valid <- n_total - n_nan_results
    et_range <- if (n_total > 0L) {
      sprintf("[%d, %d]",
              results[[1L]]$event_time,
              results[[n_total]]$event_time)
    } else {
      "N/A"
    }
    message(sprintf(
      paste0("[aggregate_to_event_time] Summary: %d event times, ",
             "%d valid, %d NaN. Range: %s. df_strategy='%s', alpha=%.3f"),
      n_total, n_valid, n_nan_results, et_range, df_strategy, alpha))
  }

  results
}
