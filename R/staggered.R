# ============================================================================
# lwdid: Staggered DiD Functions
# ============================================================================
# This file contains functions for the Staggered DiD estimation mode:
#   - prepare_staggered_controls(): Control variable extraction for (g,r) subsample
#   - estimate_staggered_effects(): (g,r) effect estimation loop (to be added)
#
# Reference: Lee & Wooldridge (2025, 2026) [lw2025, lw2026]
# ============================================================================

#' @title Prepare control variable matrix for staggered (g,r) subsample
#' @description Extracts raw control variable matrix from a (g,r) subsample
#'   and handles anomalous columns (all-NA, constant columns).
#'
#'   This function does NOT perform:
#'   \itemize{
#'     \item Sample size checks (delegated to \code{estimate_ra_common()})
#'     \item Centering (delegated to \code{estimate_ra_common()})
#'     \item Partial NA row filtering (delegated to caller via
#'       \code{complete.cases()})
#'   }
#'
#'   The returned matrix contains raw (uncentered) control variable values
#'   corresponding to the time-invariant covariates \eqn{\mathbf{X}_i} in
#'   lw2026 equation 2.17. \code{estimate_ra_common()} receives this matrix
#'   along with treatment indicator D, and internally handles three-tier
#'   degradation: Tier 1 (full interaction, lw2026 Section 2.2), Tier 2
#'   (simple controls, lw2026 eq 2.18), Tier 3 (no controls).
#'
#' @param sub data.table, (g,r) subsample (period r cross-section with
#'   treated and control units).
#' @param controls character vector, control variable names corresponding
#'   to \eqn{\mathbf{X}_i} columns in \code{sub}.
#' @param ivar character, unit identifier variable name (reserved for
#'   future diagnostics).
#' @return Matrix or NULL. Returns NULL when no controls specified, all
#'   columns are all-NA, or all columns are constant. Returned matrix
#'   preserves column names from data.table column names. Dimensions
#'   are \eqn{N_{gr} \times K'} where \eqn{K'} is the number of valid
#'   control variables after removing all-NA and constant columns.
#' @keywords internal
prepare_staggered_controls <- function(sub, controls, ivar) {
  # ── Step 0: Empty controls check ──
  K <- length(controls)
  if (K == 0L) return(NULL)

  # ── Step 1: Defensive type check ──
  # validate_inputs() should have already verified this; this is a
  # safety net for direct internal calls.
  non_numeric <- vapply(sub[, controls, with = FALSE],
                        function(col) !is.numeric(col), logical(1))
  if (any(non_numeric)) {
    stop_lwdid(
      sprintf("Control variable(s) %s are not numeric",
              paste(controls[non_numeric], collapse = ", ")),
      class = "lwdid_input"
    )
  }

  # ── Step 2: Extract raw control variable matrix ──
  # as.matrix() preserves data.table column names; downstream code
  # can access variable names via colnames(x).
  x <- as.matrix(sub[, controls, with = FALSE])

  # ── Step 3: Detect and remove all-NA columns ──
  # This typically occurs when a control variable only has values in
  # certain time periods, and the current (g,r) subsample period
  # happens to have no values for that variable.
  col_all_na <- apply(x, 2, function(col) all(is.na(col)))
  if (any(col_all_na)) {
    na_vars <- controls[col_all_na]
    warn_lwdid(
      sprintf("Control variable(s) %s are all NA in current (g,r) subsample, removed",
              paste(na_vars, collapse = ", ")),
      class = "lwdid_data",
      detail = "staggered_controls_all_na"
    )
    x <- x[, !col_all_na, drop = FALSE]
    if (ncol(x) == 0L) return(NULL)
  }

  # ── Step 4: Detect and remove constant columns (variance == 0) ──
  # Constant columns are collinear with the intercept (Tier 2) or
  # with the treatment indicator D after centering (Tier 1), causing
  # rank deficiency in QR decomposition. While qr.coef() handles
  # rank-deficient matrices gracefully (returns NA for redundant
  # coefficients), removing constant columns improves ATT estimate
  # stability and avoids confusing downstream diagnostics.
  col_var <- apply(x, 2, function(col) {
    vals <- col[!is.na(col)]
    if (length(vals) < 2L) return(0)
    var(vals)
  })
  col_const <- col_var == 0
  if (any(col_const)) {
    const_vars <- colnames(x)[col_const]
    warn_lwdid(
      sprintf("Control variable(s) %s are constant (variance=0) in current (g,r) subsample, removed",
              paste(const_vars, collapse = ", ")),
      class = "lwdid_data",
      detail = "staggered_controls_constant"
    )
    x <- x[, !col_const, drop = FALSE]
    if (ncol(x) == 0L) return(NULL)
  }

  # ── Step 5: Return raw matrix ──
  x
}


#' @title Estimate staggered (g,r) treatment effects
#' @description Estimates ATT for each valid (cohort, period) pair in
#'   staggered DiD designs. Implements lw2025 Procedure 4.1 and
#'   lw2026 equation 7.8.
#'
#'   The estimation proceeds in 8 steps for each (g,r) pair:
#'   \enumerate{
#'     \item Extract period r cross-section (period_data)
#'     \item Select control group on period_data (get_valid_controls)
#'     \item Build treatment mask and estimation sample (sub)
#'     \item Construct treatment indicator (d_sub) + min sample checks
#'     \item Apply precomputed transform (apply_precomputed_transform)
#'     \item Extract controls + two-step NA filter + min_obs check
#'     \item Run RA estimator (estimate_ra_common)
#'     \item Store results (17-column data.frame row)
#'   }
#'
#'   Corresponds to Python \code{estimate_cohort_time_effects()}.
#'
#' @param dt data.table, complete panel data
#' @param y character, outcome variable name
#' @param ivar character, unit identifier variable name
#' @param tvar character, time variable name
#' @param gvar character, cohort variable name
#' @param rolling character, transformation method
#'   (\code{"demean"} or \code{"detrend"})
#' @param control_group character, resolved control group strategy
#'   (not \code{"auto"})
#' @param controls character vector or NULL, control variable names
#' @param vce character or NULL, VCE method
#' @param cluster_var character or NULL, cluster variable name
#' @param alpha numeric, significance level (default 0.05)
#' @param pre_stats named list from \code{precompute_transforms()}
#' @param min_obs integer, minimum total observations (default 3L)
#' @param min_treated integer, minimum treated units (default 1L)
#' @param min_control integer, minimum control units (default 1L)
#' @return data.frame with 17 columns, one row per (g,r) effect
#' @keywords internal
estimate_staggered_effects <- function(
    dt, y, ivar, tvar, gvar,
    rolling, control_group, controls,
    vce, cluster_var, alpha,
    pre_stats,
    min_obs = 3L,
    min_treated = 1L,
    min_control = 1L,
    estimator = "ra",
    ps_controls = NULL,
    trim_threshold = 0.01,
    trim_method = "clip",
    n_neighbors = 1L,
    caliper = NULL,
    caliper_scale = "sd",
    with_replacement = TRUE,
    match_order = "data",
    se_method = NULL,
    boot_reps = 200L,
    seed = NULL,
    return_diagnostics = FALSE,
    parallel = FALSE,
    n_cores = NULL
) {
  # Ensure results are ordered by (cohort, period)
  cohorts <- sort(as.integer(names(pre_stats)))
  T_max <- max(dt[[tvar]], na.rm = TRUE)

  # Pre-compute all (g, r) pairs for potential parallelization
  gr_pairs <- vector("list", length = 0L)
  for (g in cohorts) {
    periods <- seq.int(g, T_max)
    for (r in periods) {
      gr_pairs[[length(gr_pairs) + 1L]] <- list(g = g, r = r)
    }
  }

  # Closure for processing a single (g, r) pair
  # Returns list(type = "result"|"skip", value = ...)
  .estimate_one_gr <- function(pair) {
    g <- pair$g
    r <- pair$r

    # ── Step 1: Extract period r cross-section ──
    period_data <- dt[get(tvar) == r]

    if (nrow(period_data) == 0L) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "no_data_in_period")))
    }

    # ── Step 2: Select control group on period_data ──
    control_mask <- tryCatch(
      get_valid_controls(period_data, gvar, g, r,
                         control_group),
      lwdid_no_control = function(e) NULL,
      error = function(e) NULL
    )
    if (is.null(control_mask)) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "no_control")))
    }

    # ── Step 3: Build treatment mask + estimation sample ──
    treat_mask <- !is_never_treated(period_data[[gvar]]) &
                  period_data[[gvar]] == g
    sample_mask <- treat_mask | control_mask
    sub <- period_data[sample_mask]

    if (nrow(sub) == 0L) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "empty_sample")))
    }

    # ── Step 4: Treatment indicator + min sample checks ──
    d_sub <- as.integer(
      !is_never_treated(sub[[gvar]]) & sub[[gvar]] == g
    )
    n_treat <- sum(d_sub == 1L)
    n_ctrl <- sum(d_sub == 0L)

    if (n_treat < min_treated) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "insufficient_treated")))
    }
    if (n_ctrl < min_control) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "insufficient_control")))
    }

    # ── Step 5: On-the-fly transformation ──
    y_trans <- apply_precomputed_transform(
      y_vals = sub[[y]],
      unit_ids = sub[[ivar]],
      pre_stat = pre_stats[[as.character(g)]],
      rolling = rolling,
      r = r,
      ivar = ivar
    )
    valid_trans <- !is.na(y_trans)

    # ── Step 6: Controls + two-step filter + min_obs ──
    x_sub <- NULL
    if (!is.null(controls) && length(controls) > 0L) {
      x_sub <- prepare_staggered_controls(
        sub, controls, ivar
      )
    }

    valid_obs <- valid_trans
    if (!is.null(x_sub)) {
      valid_obs <- valid_obs & complete.cases(x_sub)
    }

    if (sum(valid_obs) < min_obs) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "insufficient_obs")))
    }

    if (any(!valid_obs)) {
      y_trans <- y_trans[valid_obs]
      d_sub <- d_sub[valid_obs]
      sub <- sub[valid_obs]
      if (!is.null(x_sub)) {
        x_sub <- x_sub[valid_obs, , drop = FALSE]
      }
    }

    cluster_vals <- NULL
    if (identical(vce, "cluster") &&
        !is.null(cluster_var)) {
      cluster_vals <- sub[[cluster_var]]
    }

    # ── Step 7: RA estimation ──
    est <- tryCatch({
      if (estimator == "ra") {
        ra_res <- estimate_ra_common(
          y_trans, d_sub, x = x_sub,
          vce = vce,
          cluster = cluster_vals,
          alpha = alpha
        )
        ra_res$estimator <- "ra"
        ra_res$inference_dist <- "t"
        ra_res
      } else {
        est_df <- data.frame(
          .y_outcome = y_trans,
          .d_treat = d_sub
        )
        if (!is.null(x_sub)) {
          est_df <- cbind(est_df, as.data.frame(x_sub))
        }
        dispatch_estimator(
          data = est_df,
          y = ".y_outcome",
          d = ".d_treat",
          controls = controls,
          ps_controls = ps_controls,
          estimator = estimator,
          vce = vce,
          cluster_var = NULL,
          alpha = alpha,
          trim_threshold = trim_threshold,
          trim_method = trim_method,
          n_neighbors = n_neighbors,
          caliper = caliper,
          caliper_scale = caliper_scale,
          with_replacement = with_replacement,
          match_order = match_order,
          se_method = se_method,
          boot_reps = boot_reps,
          seed = seed,
          return_diagnostics = return_diagnostics
        )
      }
    },
      error = function(e) {
        NULL
      }
    )
    if (is.null(est)) {
      return(list(type = "skip",
                  value = list(g = g, r = r,
                               reason = "estimation_error")))
    }

    # ── Step 8: Store results ──
    df_inf <- if (!is.null(est$df_inference)) {
      est$df_inference
    } else if (!is.null(est$df)) {
      est$df
    } else {
      NA_real_
    }
    vce_actual <- if (is.null(vce)) "homoskedastic" else vce

    est_df <- if (!is.null(est$df)) est$df else {
      if (!is.null(est$df_resid)) est$df_resid else NA_real_
    }
    est_n <- if (!is.null(est$n)) {
      est$n
    } else if (!is.null(est$n_treated) && !is.null(est$n_control)) {
      est$n_treated + est$n_control
    } else {
      NA_integer_
    }
    est_K <- if (!is.null(est$K)) est$K else NA_integer_
    est_tier <- if (!is.null(est$controls_tier)) {
      est$controls_tier
    } else {
      NA_character_
    }

    list(type = "result", value = data.frame(
      cohort = g,
      period = r,
      event_time = r - g,
      att = est$att,
      se = est$se,
      ci_lower = est$ci_lower,
      ci_upper = est$ci_upper,
      t_stat = est$t_stat,
      pvalue = est$pvalue,
      df = est_df,
      df_inference = df_inf,
      n = est_n,
      n_treated = est$n_treated,
      n_control = est$n_control,
      K = est_K,
      controls_tier = est_tier,
      vce_type = vce_actual,
      stringsAsFactors = FALSE
    ))
  }

  # Execute: use run_parallel for unified parallel/sequential dispatch
  all_outcomes <- run_parallel(
    X = gr_pairs,
    FUN = .estimate_one_gr,
    parallel = parallel,
    n_cores = n_cores,
    fail_threshold = 1.0,
    future.seed = NULL,
    future.scheduling = 2.0,
    task_name = "Staggered (g,r) effect"
  )

  # Separate results and skipped
  results <- list()
  skipped <- list()
  for (out in all_outcomes) {
    if (identical(out$type, "result")) {
      results[[length(results) + 1L]] <- out$value
    } else {
      skipped[[length(skipped) + 1L]] <- out$value
    }
  }

  # ── Report skipped (g,r) pairs ──
  if (length(skipped) > 0L) {
    n_skipped <- length(skipped)
    n_total_pairs <- length(gr_pairs)
    reasons <- vapply(
      skipped, function(x) x$reason, character(1)
    )
    reason_counts <- table(reasons)
    reason_summary <- paste(
      sprintf("%s: %d", names(reason_counts),
              as.integer(reason_counts)),
      collapse = ", "
    )
    warn_lwdid(
      sprintf(
        "%d/%d (g,r) pairs skipped (%s)",
        n_skipped, n_total_pairs, reason_summary
      ),
      class = "lwdid_data",
      detail = "gr_pairs_skipped",
      action_taken = sprintf(
        "%d/%d pairs skipped",
        n_skipped, n_total_pairs
      ),
      skipped_pairs = skipped
    )
  }

  if (length(results) == 0L) {
    stop_lwdid(
      "No valid (g,r) effect estimates produced",
      class = "lwdid_insufficient_data"
    )
  }

  do.call(rbind, results)
}
