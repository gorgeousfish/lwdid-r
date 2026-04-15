# pretreatment.R
# Pre-treatment effect estimation for Common Timing DiD.
# Paper: lw2025 Section 5 (equations 5.1-5.7, Procedure 5.1)
# Symmetric transforms use future pre-treatment periods {t+1,...,S-1}
# as reference, not all pre-treatment periods {1,...,S-1}.

#' @title Pre-treatment effect estimation (Common Timing)
#' @description Estimates ATT for each pre-treatment period using
#'   symmetric transforms. Under parallel trends, pre-treatment
#'   ATTs should be approximately zero (placebo test).
#'   Receives raw panel data and performs symmetric transforms
#'   internally for each pre-treatment period.
#' @param data data.table, raw panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param d character, treatment indicator variable name (0/1)
#' @param tvar character, time variable name
#' @param tpost1 integer, first post-treatment period (treatment
#'   time S)
#' @param rolling character, transform method ("demean"/"detrend")
#' @param controls character vector or NULL, control variable names
#' @param estimator character, estimator type
#'   ("ra"/"ipw"/"ipwra"/"psm")
#' @param vce character or NULL, VCE type
#' @param cluster_var character or NULL, cluster variable name
#' @param alpha numeric, significance level (default 0.05)
#' @param exclude_pre_periods integer, exclude last k pre-treatment
#'   periods (default 0, PRD section 8e)
#' @return data.frame with 14 columns, or NULL if no pre-treatment
#'   periods
#' @keywords internal
estimate_pre_treatment_common <- function(
    data, y, ivar, d, tvar, tpost1,
    rolling = "demean",
    controls = NULL,
    estimator = "ra",
    vce = NULL,
    cluster_var = NULL,
    alpha = 0.05,
    exclude_pre_periods = 0L) {

  # ===== Step 0: Input validation =====
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  if (nrow(data) == 0L) {
    stop_lwdid("data is empty", class = "lwdid_invalid_param")
  }
  required_cols <- c(y, ivar, d, tvar)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      sprintf("Missing required columns: %s",
              paste(missing_cols, collapse = ", ")),
      class = "lwdid_invalid_param"
    )
  }
  if (!rolling %in% c("demean", "detrend")) {
    stop_lwdid(
      sprintf("rolling must be 'demean' or 'detrend', got '%s'",
              rolling),
      class = "lwdid_invalid_param"
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(
      sprintf("alpha must be in (0, 1), got: %s",
              deparse(alpha)),
      class = "lwdid_invalid_param"
    )
  }

  # ===== Step 1: Get pre-treatment periods =====
  all_periods <- sort(unique(data[[tvar]]))
  effective_tpost1 <- tpost1
  if (exclude_pre_periods > 0L) {
    effective_tpost1 <- tpost1 - exclude_pre_periods
  }
  pre_periods <- all_periods[all_periods < effective_tpost1]

  if (length(pre_periods) == 0L) {
    warn_lwdid(
      "No pre-treatment periods available",
      class = "lwdid_insufficient_pre_periods"
    )
    return(NULL)
  }

  # ===== Step 2: Determine anchor and estimable periods =====
  # Anchor period: last pre-treatment period before effective_tpost1
  # (event_time = -1 relative to effective_tpost1). ATT = 0 by
  # construction. When exclude_pre_periods > 0, anchor shifts with
  # the effective treatment boundary.
  anchor_period <- effective_tpost1 - 1L
  # [R design] Exclude T_min (first pre-period): demean
  # produces low variance at T_min
  estimable_periods <- pre_periods[-1]  # drop T_min
  estimable_periods <- estimable_periods[
    estimable_periods != anchor_period
  ]

  # ===== Step 3: Build anchor result =====
  # lw2025 Section 5: anchor (t=S-1) has ATT=0 by symmetric
  # transform definition (no reference periods).
  # Matches Python estimation_pre.py L296-L320.
  anchor_mask <- data[[tvar]] == anchor_period
  n_treat_anchor <- sum(data[[d]][anchor_mask] == 1L)
  n_ctrl_anchor <- sum(data[[d]][anchor_mask] == 0L)
  anchor_result <- data.frame(
    cohort = NA_integer_,
    period = as.integer(anchor_period),
    event_time = as.integer(anchor_period - tpost1),
    att = 0.0,
    se = 0.0,
    ci_lower = 0.0,
    ci_upper = 0.0,
    t_stat = NaN,
    pvalue = NaN,
    n_treated = as.integer(n_treat_anchor),
    n_control = as.integer(n_ctrl_anchor),
    is_anchor = TRUE,
    rolling_window_size = 0L,
    df_inference = as.integer(
      max(n_treat_anchor + n_ctrl_anchor - 2L, 0L)
    ),
    stringsAsFactors = FALSE
  )

  # ===== Step 4: Loop over estimable periods =====
  results <- vector("list", length(estimable_periods))
  for (k in seq_along(estimable_periods)) {
    t_k <- estimable_periods[k]
    event_time_k <- as.integer(t_k - tpost1)

    # --- 4a: Reference periods {t_k+1, ..., S-1} ---
    # Uses original tpost1 (not effective_tpost1) for ref periods
    ref_periods <- all_periods[
      all_periods > t_k & all_periods < tpost1
    ]
    if (length(ref_periods) == 0L) next

    # --- 4b: Extract reference + current period data ---
    ref_data <- data[data[[tvar]] %in% ref_periods]
    t_data <- data.table::copy(data[data[[tvar]] == t_k])

    # --- 4c: Symmetric transform (by ivar) ---
    # lw2025 eq 2.12 symmetric version (demean) /
    # Procedure 5.1 eq 5.6 (detrend)
    # Capture column names for data.table NSE scoping
    y_col <- y
    t_col <- tvar
    if (rolling == "demean") {
      # demean: y_trans = Y_it - mean(Y_i, {t+1,...,S-1})
      ref_means <- ref_data[,
        .(ref_val = mean(.SD[[y_col]], na.rm = TRUE)),
        by = c(ivar), .SDcols = y_col
      ]
      t_data <- merge(t_data, ref_means, by = ivar,
                       all.x = TRUE)
      t_data[, y_trans_pre := .SD[[y_col]] - ref_val,
             .SDcols = y_col]
    } else {
      # detrend: y_trans = Y_it - (alpha_hat + beta_hat * t)
      # OLS on ref periods per unit
      if (length(ref_periods) < 2L) {
        # [R design] Degrade to demean when < 2 ref periods
        ref_means <- ref_data[,
          .(ref_val = mean(.SD[[y_col]], na.rm = TRUE)),
          by = c(ivar), .SDcols = y_col
        ]
        t_data <- merge(t_data, ref_means, by = ivar,
                         all.x = TRUE)
        t_data[, y_trans_pre := .SD[[y_col]] - ref_val,
               .SDcols = y_col]
      } else {
        # Full detrend: fit Y_iq = alpha + beta*q per unit
        sd_cols <- c(y_col, t_col)
        trend_fits <- ref_data[, {
          y_val <- .SD[[y_col]]
          t_val <- .SD[[t_col]]
          fit <- lm.fit(cbind(1, t_val), y_val)
          list(
            alpha_hat = fit$coefficients[1],
            beta_hat = fit$coefficients[2]
          )
        }, by = c(ivar), .SDcols = sd_cols]
        t_data <- merge(t_data, trend_fits, by = ivar,
                         all.x = TRUE)
        t_data[,
          y_trans_pre := .SD[[y_col]] -
            (alpha_hat + beta_hat * .SD[[t_col]]),
          .SDcols = sd_cols
        ]
      }
    }

    # Remove units with NA transform (missing in ref periods)
    t_data <- t_data[!is.na(y_trans_pre)]

    # --- 4d: Sample size check ---
    d_k <- as.integer(t_data[[d]])
    n_treat_k <- sum(d_k == 1L)
    n_ctrl_k <- sum(d_k == 0L)
    n_total_k <- n_treat_k + n_ctrl_k

    if (n_total_k < 3L || n_treat_k < 1L || n_ctrl_k < 1L) {
      results[[k]] <- data.frame(
        cohort = NA_integer_,
        period = as.integer(t_k),
        event_time = event_time_k,
        att = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        t_stat = NA_real_,
        pvalue = NA_real_,
        n_treated = as.integer(n_treat_k),
        n_control = as.integer(n_ctrl_k),
        is_anchor = FALSE,
        rolling_window_size = NA_integer_,
        df_inference = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # --- 4e: Estimate pre-treatment ATT ---
    # Route through dispatch_estimator() for ra/ipw/ipwra/psm
    # Build estimation data.frame for dispatch_estimator
    est_result <- tryCatch({
      est_df <- data.frame(
        .y_outcome = t_data$y_trans_pre,
        .d_treat = d_k,
        stringsAsFactors = FALSE
      )
      ctrl_names <- NULL
      if (!is.null(controls) && length(controls) > 0L) {
        ctrl_df <- as.data.frame(
          t_data[, controls, with = FALSE]
        )
        est_df <- cbind(est_df, ctrl_df)
        ctrl_names <- controls
      }
      dispatch_estimator(
        data = est_df,
        y = ".y_outcome",
        d = ".d_treat",
        controls = ctrl_names,
        estimator = estimator,
        vce = vce,
        cluster_var = cluster_var,
        alpha = alpha
      )
    }, error = function(e) {
      warn_lwdid(
        sprintf("Pre-period %d estimator '%s' failed: %s",
                t_k, estimator, conditionMessage(e)),
        class = "lwdid_data"
      )
      NULL
    })

    if (is.null(est_result)) {
      results[[k]] <- data.frame(
        cohort = NA_integer_,
        period = as.integer(t_k),
        event_time = event_time_k,
        att = NA_real_,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        t_stat = NA_real_,
        pvalue = NA_real_,
        n_treated = as.integer(n_treat_k),
        n_control = as.integer(n_ctrl_k),
        is_anchor = FALSE,
        rolling_window_size = NA_integer_,
        df_inference = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Extract df_inference with fallback chain
    df_inf <- if (!is.null(est_result$df_inference)) {
      est_result$df_inference
    } else if (!is.null(est_result$df)) {
      est_result$df
    } else {
      NA_integer_
    }

    rws <- as.integer(tpost1 - 1L - t_k)  # |{t+1,...,S-1}|

    results[[k]] <- data.frame(
      cohort = NA_integer_,
      period = as.integer(t_k),
      event_time = event_time_k,
      att = est_result$att,
      se = est_result$se,
      ci_lower = est_result$ci_lower,
      ci_upper = est_result$ci_upper,
      t_stat = est_result$t_stat,
      pvalue = est_result$pvalue,
      n_treated = as.integer(n_treat_k),
      n_control = as.integer(n_ctrl_k),
      is_anchor = FALSE,
      rolling_window_size = rws,
      df_inference = as.integer(df_inf),
      stringsAsFactors = FALSE
    )
  }

  # ===== Step 5: Combine and return =====
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) > 0L) {
    all_results <- do.call(rbind, c(list(anchor_result), results))
  } else {
    all_results <- anchor_result
  }
  # Sort by event_time descending (anchor at top)
  all_results <- all_results[order(-all_results$event_time), ]
  rownames(all_results) <- NULL
  all_results
}


# ============================================================================
# Staggered Pre-treatment Effect Estimation
# Paper: lw2025 Section 5, Procedure 5.1 (rolling window adaptation)
# ============================================================================

#' @title Pre-treatment effect estimation (Staggered)
#' @description Estimates ATT for each pre-treatment period per cohort
#'   using symmetric transforms. Under parallel trends, pre-treatment
#'   ATTs should be approximately zero (placebo test).
#'   Implements lw2025 Procedure 5.1 (staggered entry) adapted for
#'   pre-treatment periods with rolling window reference sets.
#' @param data data.table, raw panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param gvar character, cohort (treatment timing) variable name
#' @param rolling character, transform method ("demean"/"detrend")
#' @param estimator character, estimator type ("ra"/"ipw"/"ipwra"/"psm")
#' @param controls character vector or NULL, control variable names
#' @param control_group character, "not_yet_treated" or "never_treated"
#' @param vce character or NULL, VCE type
#' @param cluster_var character or NULL, cluster variable name
#' @param alpha numeric, significance level (default 0.05)
#' @return data.frame with 14 columns, or NULL if no estimable periods
#' @keywords internal
estimate_pre_treatment_staggered <- function(
    data, y, ivar, tvar, gvar,
    rolling = "demean",
    estimator = "ra",
    controls = NULL,
    control_group = "not_yet_treated",
    vce = NULL,
    cluster_var = NULL,
    alpha = 0.05,
    include_earliest_pre_period = FALSE) {

  # ===== Step 0: Input validation =====
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  if (nrow(data) == 0L) {
    stop_lwdid("data is empty", class = "lwdid_invalid_param")
  }
  required_cols <- c(y, ivar, tvar, gvar)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop_lwdid(
      sprintf("Missing required columns: %s",
              paste(missing_cols, collapse = ", ")),
      class = "lwdid_invalid_param"
    )
  }
  if (!rolling %in% c("demean", "detrend")) {
    stop_lwdid(
      sprintf("rolling must be 'demean' or 'detrend', got '%s'",
              rolling),
      class = "lwdid_invalid_param"
    )
  }
  if (!control_group %in% c("not_yet_treated", "never_treated")) {
    stop_lwdid(
      sprintf("control_group must be 'not_yet_treated' or 'never_treated', got '%s'",
              control_group),
      class = "lwdid_invalid_param"
    )
  }
  if (!is.null(vce) && vce == "cluster" && is.null(cluster_var)) {
    stop_lwdid(
      "cluster_var must be specified when vce='cluster'",
      class = "lwdid_invalid_param"
    )
  }

  # ===== Step 1: Get cohorts and periods =====
  cohorts <- get_cohorts(data, gvar, ivar)
  all_periods <- sort(unique(data[[tvar]]))
  t_min <- min(all_periods)

  # Column name captures for data.table NSE
  y_col <- y
  t_col <- tvar
  g_col <- gvar

  # ===== Step 2: Main loop over cohorts =====
  results_list <- list()
  idx <- 0L

  for (g in cohorts) {
    # --- 2a: Estimable pre-periods (exclude T_min and anchor g-1) ---
    # Default package path keeps T_min excluded for the legacy pre-treatment
    # interface. Trend diagnostics can opt in to the full Python-style
    # placebo enumeration and retain the earliest estimable pre-period.
    if (isTRUE(include_earliest_pre_period)) {
      pre_periods_g <- all_periods[
        all_periods >= t_min & all_periods <= (g - 2L)
      ]
    } else {
      pre_periods_g <- all_periods[
        all_periods > t_min & all_periods <= (g - 2L)
      ]
    }

    # --- 2b: Anchor (t = g-1) ---
    anchor_period <- g - 1L
    anchor_mask <- data[[tvar]] == anchor_period
    g_vals_anchor <- data[[gvar]][anchor_mask]

    n_treat_anchor <- sum(
      !is_never_treated(g_vals_anchor) &
        abs(g_vals_anchor - g) < 1e-8,
      na.rm = TRUE
    )

    # Control group for anchor
    if (control_group == "not_yet_treated") {
      # G_i > g-1 (strict) AND G_i != g, OR never-treated
      ctrl_anchor <- sum(
        (g_vals_anchor > anchor_period &
           abs(g_vals_anchor - g) >= 1e-8 &
           !is_never_treated(g_vals_anchor)) |
          is_never_treated(g_vals_anchor),
        na.rm = TRUE
      )
    } else {
      # never_treated only
      ctrl_anchor <- sum(is_never_treated(g_vals_anchor),
                         na.rm = TRUE)
    }

    idx <- idx + 1L
    results_list[[idx]] <- data.frame(
      cohort = as.integer(g),
      period = as.integer(anchor_period),
      event_time = -1L,
      att = 0.0,
      se = 0.0,
      ci_lower = 0.0,
      ci_upper = 0.0,
      t_stat = NaN,
      pvalue = NaN,
      n_treated = as.integer(n_treat_anchor),
      n_control = as.integer(ctrl_anchor),
      is_anchor = TRUE,
      rolling_window_size = 0L,
      df_inference = as.integer(
        max(n_treat_anchor + ctrl_anchor - 2L, 0L)
      ),
      stringsAsFactors = FALSE
    )

    # --- 2c: Loop over estimable pre-periods ---
    if (length(pre_periods_g) == 0L) next

    for (t_pre in pre_periods_g) {
      event_time_k <- as.integer(t_pre - g)

      # Reference periods: {t_pre+1, ..., g-1}
      ref_periods <- all_periods[
        all_periods > t_pre & all_periods < g
      ]
      if (length(ref_periods) == 0L) next

      # Extract reference and current period data
      ref_data <- data[data[[tvar]] %in% ref_periods]
      t_data <- data.table::copy(data[data[[tvar]] == t_pre])

      # --- Symmetric transform (by ivar) ---
      # lw2025 Procedure 5.1 rolling window adaptation
      if (rolling == "demean") {
        # demean: y_trans = Y_it - mean(Y_i, ref_periods)
        ref_means <- ref_data[,
          .(ref_val = mean(.SD[[y_col]], na.rm = TRUE)),
          by = c(ivar), .SDcols = y_col
        ]
        t_data <- merge(t_data, ref_means, by = ivar,
                         all.x = TRUE)
        t_data[, y_trans_pre := .SD[[y_col]] - ref_val,
               .SDcols = y_col]
      } else {
        # detrend: y_trans = Y_it - (alpha_hat + beta_hat * t)
        if (length(ref_periods) < 2L) {
          # Degrade to demean when < 2 ref periods
          ref_means <- ref_data[,
            .(ref_val = mean(.SD[[y_col]], na.rm = TRUE)),
            by = c(ivar), .SDcols = y_col
          ]
          t_data <- merge(t_data, ref_means, by = ivar,
                           all.x = TRUE)
          t_data[, y_trans_pre := .SD[[y_col]] - ref_val,
                 .SDcols = y_col]
        } else {
          # Full detrend: OLS Y_iq = alpha + beta*q per unit
          sd_cols <- c(y_col, t_col)
          trend_fits <- ref_data[, {
            yv <- .SD[[y_col]]
            tv <- .SD[[t_col]]
            fit <- lm.fit(cbind(1, tv), yv)
            list(
              alpha_hat = fit$coefficients[1],
              beta_hat = fit$coefficients[2]
            )
          }, by = c(ivar), .SDcols = sd_cols]
          t_data <- merge(t_data, trend_fits, by = ivar,
                           all.x = TRUE)
          t_data[,
            y_trans_pre := .SD[[y_col]] -
              (alpha_hat + beta_hat * .SD[[t_col]]),
            .SDcols = sd_cols
          ]
        }
      }

      # Remove units with NA transform
      t_data <- t_data[!is.na(y_trans_pre)]

      # --- Control group selection (FATAL-001) ---
      g_vals_k <- t_data[[gvar]]
      nt_mask <- is_never_treated(g_vals_k)

      if (control_group == "not_yet_treated") {
        # Keep: treated (G_i==g) OR not-yet-treated (G_i>t) OR NT
        keep_mask <- (!nt_mask & abs(g_vals_k - g) < 1e-8) |
          (!nt_mask & g_vals_k > t_pre) |
          nt_mask
      } else {
        # Keep: treated (G_i==g) OR never-treated
        keep_mask <- (!nt_mask & abs(g_vals_k - g) < 1e-8) |
          nt_mask
      }
      t_data <- t_data[keep_mask]

      if (nrow(t_data) == 0L) next

      # Treatment indicator
      g_vals_sub <- t_data[[gvar]]
      d_pre <- as.integer(
        !is_never_treated(g_vals_sub) &
          abs(g_vals_sub - g) < 1e-8
      )

      # --- Sample size check ---
      n_treat_k <- sum(d_pre == 1L)
      n_ctrl_k <- sum(d_pre == 0L)
      n_total_k <- n_treat_k + n_ctrl_k

      if (n_total_k < 3L || n_treat_k < 1L || n_ctrl_k < 1L) {
        idx <- idx + 1L
        results_list[[idx]] <- data.frame(
          cohort = as.integer(g),
          period = as.integer(t_pre),
          event_time = event_time_k,
          att = NA_real_,
          se = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          t_stat = NA_real_,
          pvalue = NA_real_,
          n_treated = as.integer(n_treat_k),
          n_control = as.integer(n_ctrl_k),
          is_anchor = FALSE,
          rolling_window_size = NA_integer_,
          df_inference = NA_integer_,
          stringsAsFactors = FALSE
        )
        next
      }

      # --- Estimate pre-treatment ATT ---
      est_result <- tryCatch({
        est_df <- data.frame(
          .y_outcome = t_data$y_trans_pre,
          .d_treat = d_pre,
          stringsAsFactors = FALSE
        )
        ctrl_names <- NULL
        if (!is.null(controls) && length(controls) > 0L) {
          ctrl_df <- as.data.frame(
            t_data[, controls, with = FALSE]
          )
          est_df <- cbind(est_df, ctrl_df)
          ctrl_names <- controls
        }
        dispatch_estimator(
          data = est_df,
          y = ".y_outcome",
          d = ".d_treat",
          controls = ctrl_names,
          estimator = estimator,
          vce = vce,
          cluster_var = cluster_var,
          alpha = alpha
        )
      }, error = function(e) {
        warn_lwdid(
          sprintf("Cohort %d pre-period %d estimator '%s' failed: %s",
                  g, t_pre, estimator, conditionMessage(e)),
          class = "lwdid_data"
        )
        NULL
      })

      if (is.null(est_result)) {
        idx <- idx + 1L
        results_list[[idx]] <- data.frame(
          cohort = as.integer(g),
          period = as.integer(t_pre),
          event_time = event_time_k,
          att = NA_real_,
          se = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          t_stat = NA_real_,
          pvalue = NA_real_,
          n_treated = as.integer(n_treat_k),
          n_control = as.integer(n_ctrl_k),
          is_anchor = FALSE,
          rolling_window_size = NA_integer_,
          df_inference = NA_integer_,
          stringsAsFactors = FALSE
        )
        next
      }

      # Extract df_inference with fallback
      df_inf <- if (!is.null(est_result$df_inference)) {
        est_result$df_inference
      } else if (!is.null(est_result$df)) {
        est_result$df
      } else {
        NA_integer_
      }

      rws <- as.integer(g - t_pre - 1L)

      idx <- idx + 1L
      results_list[[idx]] <- data.frame(
        cohort = as.integer(g),
        period = as.integer(t_pre),
        event_time = event_time_k,
        att = est_result$att,
        se = est_result$se,
        ci_lower = est_result$ci_lower,
        ci_upper = est_result$ci_upper,
        t_stat = est_result$t_stat,
        pvalue = est_result$pvalue,
        n_treated = as.integer(n_treat_k),
        n_control = as.integer(n_ctrl_k),
        is_anchor = FALSE,
        rolling_window_size = rws,
        df_inference = as.integer(df_inf),
        stringsAsFactors = FALSE
      )
    }
  }

  # ===== Step 3: Combine and sort =====
  if (length(results_list) == 0L) return(NULL)
  result_df <- do.call(rbind, results_list)
  # Sort: cohort ascending, event_time descending
  result_df <- result_df[order(result_df$cohort,
                                -result_df$event_time), ]
  rownames(result_df) <- NULL
  result_df
}


# ============================================================================
# Parallel Trends Joint Test
# Paper: LW2025 Equation (2.3) / (2.15) Parallel Trends Hypothesis
# PRD §3b Joint Test; PRD §8d Pre-treatment Parameters
# Python: lwdid-py_v0.2.3/src/lwdid/staggered/parallel_trends.py
# ============================================================================

#' @title Parallel trends joint test
#' @description Tests whether all pre-treatment ATTs are jointly zero
#'   (LW2025 equation 2.3 parallel trends hypothesis). Supports F-test
#'   and Wald chi-squared test. Receives pre-treatment effect estimates
#'   from Story E7-03/E7-04 and performs joint null hypothesis test.
#' @param pre_treatment_effects data.frame, pre-treatment effect
#'   estimates (from estimate_pre_treatment_common/staggered). Must
#'   contain columns: att, se, t_stat, pvalue, is_anchor, n_treated,
#'   n_control, event_time, cohort, period.
#' @param alpha numeric, significance level (default 0.05)
#' @param test_type character, test type: "f" (F-test, default) or
#'   "wald" (Wald chi-squared test)
#' @param min_pre_periods integer, minimum pre-treatment periods
#'   required for joint test (default 2L)
#' @return list with 14 fields:
#'   joint_stat, joint_pvalue, joint_df1, joint_df2, reject_null,
#'   n_pre_periods, excluded_periods, n_significant_05,
#'   n_significant_10, max_abs_pre_att, individual_tests, alpha,
#'   test_type, interpretation
#' @keywords internal
test_parallel_trends_joint <- function(
    pre_treatment_effects,
    alpha = 0.05,
    test_type = c("f", "wald"),
    min_pre_periods = 2L) {

  # ===== Step 0: Input validation =====
  if (nrow(pre_treatment_effects) == 0L) {
    stop_lwdid(
      "pre_treatment_effects is empty",
      class = "lwdid_invalid_param"
    )
  }
  test_type <- match.arg(test_type)

  if (!is.numeric(min_pre_periods) || length(min_pre_periods) != 1L ||
      is.na(min_pre_periods) || min_pre_periods < 1L ||
      min_pre_periods != as.integer(min_pre_periods)) {
    stop_lwdid(
      sprintf("min_pre_periods must be a positive integer, got: %s",
              deparse(min_pre_periods)),
      class = "lwdid_invalid_param"
    )
  }
  min_pre_periods <- as.integer(min_pre_periods)

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop_lwdid(
      sprintf("alpha must be in (0, 1), got: %s", deparse(alpha)),
      class = "lwdid_invalid_param"
    )
  }

  # ===== Step 1: Exclusion filtering =====
  # Exclude: anchor rows, NA att, NA/non-positive se
  n_rows <- nrow(pre_treatment_effects)
  keep <- rep(TRUE, n_rows)
  for (i in seq_len(n_rows)) {
    row_i <- pre_treatment_effects[i, ]
    if (isTRUE(row_i$is_anchor)) {
      keep[i] <- FALSE
    } else if (is.na(row_i$att)) {
      keep[i] <- FALSE
    } else if (is.na(row_i$se) || row_i$se <= 0) {
      keep[i] <- FALSE
    }
  }
  excluded_periods <- pre_treatment_effects$event_time[!keep]
  valid <- pre_treatment_effects[keep, , drop = FALSE]
  K <- nrow(valid)

  # ===== Step 2: Empty individual_tests column names =====
  indiv_colnames <- c("event_time", "cohort", "period", "att",
                       "se", "t_stat", "pvalue", "significant",
                       "n_treated", "n_control")

  # ===== Step 3: K < min_pre_periods warning =====
  if (K < min_pre_periods) {
    warn_lwdid(
      sprintf(
        "Only %d valid pre-treatment period(s), need >= %d for reliable joint test",
        K, min_pre_periods
      ),
      class = "lwdid_insufficient_pre_periods"
    )
  }

  # ===== Step 4: K=0 early return =====
  if (K == 0L) {
    warn_lwdid(
      "No valid pre-treatment periods for joint test (all excluded)",
      class = "lwdid_insufficient_pre_periods"
    )
    empty_indiv <- data.frame(
      event_time = integer(0),
      cohort = integer(0),
      period = integer(0),
      att = numeric(0),
      se = numeric(0),
      t_stat = numeric(0),
      pvalue = numeric(0),
      significant = logical(0),
      n_treated = integer(0),
      n_control = integer(0),
      stringsAsFactors = FALSE
    )
    return(list(
      joint_stat = NA_real_,
      joint_pvalue = NA_real_,
      joint_df1 = 0L,
      joint_df2 = 0L,
      reject_null = FALSE,
      n_pre_periods = 0L,
      excluded_periods = excluded_periods,
      n_significant_05 = 0L,
      n_significant_10 = 0L,
      max_abs_pre_att = NA_real_,
      individual_tests = empty_indiv,
      alpha = alpha,
      test_type = test_type,
      interpretation = "No valid pre-treatment periods available for joint test."
    ))
  }

  # ===== Step 5: Extract t-statistics =====
  t_stats <- valid$t_stat
  sum_t2 <- sum(t_stats^2)

  # ===== Step 6: Joint test statistic =====
  if (test_type == "f") {
    # F-test: F = mean(t²), df2 from average sample sizes
    joint_stat <- sum_t2 / K
    df1 <- as.integer(K)
    avg_n <- mean(valid$n_treated) + mean(valid$n_control)
    df2 <- max(as.integer(avg_n - 2), 1L)
    joint_pvalue <- stats::pf(joint_stat, df1, df2,
                               lower.tail = FALSE)
  } else {
    # Wald chi-squared: W = sum(t²)
    joint_stat <- sum_t2
    df1 <- as.integer(K)
    df2 <- 0L
    joint_pvalue <- stats::pchisq(joint_stat, df = K,
                                   lower.tail = FALSE)
  }

  # ===== Step 7: Summary statistics =====
  n_significant_05 <- sum(valid$pvalue < 0.05, na.rm = TRUE)
  n_significant_10 <- sum(valid$pvalue < 0.10, na.rm = TRUE)
  max_abs_pre_att <- max(abs(valid$att), na.rm = TRUE)

  # ===== Step 8: reject_null =====
  reject_null <- if (!is.na(joint_pvalue)) {
    joint_pvalue < alpha
  } else {
    FALSE
  }

  # ===== Step 9: Interpretation =====
  if (reject_null) {
    interpretation <- sprintf(
      "Joint test rejects H0 (p=%.4f < alpha=%.2f): evidence against parallel trends.",
      joint_pvalue, alpha
    )
  } else {
    interpretation <- sprintf(
      "Joint test does not reject H0 (p=%.4f >= alpha=%.2f): consistent with parallel trends.",
      joint_pvalue, alpha
    )
  }

  # ===== Step 10: individual_tests data.frame =====
  individual_tests <- data.frame(
    event_time = valid$event_time,
    cohort = valid$cohort,
    period = valid$period,
    att = valid$att,
    se = valid$se,
    t_stat = valid$t_stat,
    pvalue = valid$pvalue,
    significant = valid$pvalue < alpha,
    n_treated = valid$n_treated,
    n_control = valid$n_control,
    stringsAsFactors = FALSE
  )
  rownames(individual_tests) <- NULL

  # ===== Step 11: Return =====
  list(
    joint_stat = joint_stat,
    joint_pvalue = joint_pvalue,
    joint_df1 = df1,
    joint_df2 = df2,
    reject_null = reject_null,
    n_pre_periods = as.integer(K),
    excluded_periods = excluded_periods,
    n_significant_05 = as.integer(n_significant_05),
    n_significant_10 = as.integer(n_significant_10),
    max_abs_pre_att = max_abs_pre_att,
    individual_tests = individual_tests,
    alpha = alpha,
    test_type = test_type,
    interpretation = interpretation
  )
}
