#' @title lwdid Result Class
#' @description S3 class for lwdid estimation results. Mirrors the Python
#'   `LWDIDResults` class structure with all core, staggered, RI, and
#'   diagnostic attributes.
#' @name lwdid_result
#' @family lwdid-results
NULL

# -- Module-level constants -------------------------------------------------
.VCE_DISPLAY_MAP <- c(
  "ols"       = "OLS (Homoskedastic)",
  "robust"    = "HC1 (Heteroskedasticity-robust)",
  "hc0"       = "HC0 (White)",
  "hc1"       = "HC1 (Heteroskedasticity-robust)",
  "hc2"       = "HC2 (Bell-McCaffrey)",
  "hc3"       = "HC3 (Small-sample adjusted)",
  "hc4"       = "HC4 (Cribari-Neto)",
  "bootstrap" = "Bootstrap"
)

.VALID_AGGREGATE_TYPES <- c("none", "cohort", "overall", "event_time")

.ESTIMATOR_DISPLAY_MAP <- c(
  "ra"    = "Regression Adjustment (RA)",
  "ipw"   = "Inverse Probability Weighting (IPW)",
  "ipwra" = "IPW with Regression Adjustment (IPWRA)",
  "psm"   = "Propensity Score Matching (PSM)"
)

.resolve_top_level_inference_dist <- function(estimator = NULL) {
  if (!is.character(estimator) || length(estimator) != 1L || is.na(estimator)) {
    return(NA_character_)
  }

  estimator <- tolower(estimator)
  if (identical(estimator, "ra")) {
    return("t")
  }
  if (estimator %in% c("ipw", "ipwra", "psm")) {
    return("normal")
  }

  NA_character_
}

#' Construct a new lwdid_result object
#'
#' @param att Numeric, ATT point estimate.
#' @param se_att Numeric, standard error of ATT.
#' @param t_stat Numeric, t-statistic.
#' @param pvalue Numeric, two-sided p-value.
#' @param ci_lower Numeric, lower CI bound.
#' @param ci_upper Numeric, upper CI bound.
#' @param df_resid Integer, residual degrees of freedom.
#' @param df_inference Integer, inference degrees of freedom.
#' @param nobs Integer, number of observations.
#' @param n_treated Integer, number of treated units.
#' @param n_control Integer, number of control units.
#' @param K Integer, last pre-treatment period index.
#' @param tpost1 Integer, first post-treatment period index.
#' @param depvar Character, dependent variable name.
#' @param rolling Character, transformation method used.
#' @param vce_type Character or NULL, VCE type used.
#' @param cluster_var Character or NULL, cluster variable name.
#' @param n_clusters Integer or NULL, number of clusters.
#' @param estimator Character or NULL, estimator used.
#' @param method Character, "common_timing" or "staggered".
#' @param alpha Numeric, significance level.
#' @param is_staggered Logical, whether staggered mode.
#' @param controls_used Logical, whether controls were used.
#' @param controls Character vector, control variable names.
#' @param include_pretreatment Logical.
#' @param control_group Character or NULL, user-specified control group.
#' @param control_group_used Character or NULL, actual control group used.
#' @param params Numeric vector or NULL, full coefficient vector.
#' @param bse Numeric vector or NULL, coefficient standard errors.
#' @param vcov_matrix Matrix or NULL, variance-covariance matrix.
#' @param resid Numeric vector or NULL, residuals.
#' @param data data.table or NULL, transformed regression data.
#' @param att_by_period data.frame or NULL, period-specific ATT.
#' @param att_pre_treatment data.frame or NULL, pre-treatment ATT.
#' @param parallel_trends_test List or NULL, joint F-test results.
#' @param aggregate Character or NULL, aggregation level (staggered).
#' @param cohorts Integer vector or NULL, treatment cohorts.
#' @param cohort_sizes Named integer or NULL, units per cohort.
#' @param n_never_treated Integer or NULL, number of never-treated units.
#' @param n_units Integer or NULL, number of units contributing to staggered
#'   aggregation output.
#' @param n_periods Integer or NULL, number of periods contributing to staggered
#'   aggregation output.
#' @param n_cohorts Integer or NULL, number of contributing cohorts.
#' @param att_by_cohort data.frame or NULL, cohort-level ATT.
#' @param att_by_cohort_time data.frame or NULL, (g,r)-specific ATT.
#' @param att_overall Numeric or NULL, overall weighted ATT.
#' @param se_overall Numeric or NULL, overall ATT SE.
#' @param ci_overall_lower Numeric or NULL.
#' @param ci_overall_upper Numeric or NULL.
#' @param t_stat_overall Numeric or NULL.
#' @param pvalue_overall Numeric or NULL.
#' @param cohort_effects List or NULL.
#' @param event_time_effects data.frame or NULL.
#' @param att_cohort_agg Numeric or NULL, cohort-aggregated ATT.
#' @param se_cohort_agg Numeric or NULL, standard error for the
#'   cohort-aggregated ATT.
#' @param t_stat_cohort_agg Numeric or NULL, t-statistic for the
#'   cohort-aggregated ATT.
#' @param pvalue_cohort_agg Numeric or NULL, p-value for the
#'   cohort-aggregated ATT.
#' @param ci_cohort_agg Numeric vector of length 2 or NULL, confidence interval
#'   for the cohort-aggregated ATT.
#' @param cohort_weights Named numeric or NULL.
#' @param ri_pvalue Numeric or NULL, RI p-value.
#' @param ri_seed Integer or NULL.
#' @param rireps Integer or NULL.
#' @param ri_method Character or NULL.
#' @param ri_valid Integer or NULL.
#' @param ri_failed Integer or NULL.
#' @param ri_error Character or NULL.
#' @param ri_target Character or NULL.
#' @param ri_distribution Numeric vector or NULL.
#' @param diagnostics List or NULL.
#' @param warning_diagnostics List, warning diagnostics from registry.
#' @param propensity_scores Numeric vector or NULL.
#' @param matched_data data.frame or NULL.
#' @param n_matched Integer or NULL.
#' @param match_rate Numeric or NULL.
#' @param weights_cv Numeric or NULL.
#' @param warnings_log List, raw warnings log.
#' @param call Call object or NULL.
#' @param lwdid_version Character or NULL.
#' @param ivar Character, unit identifier variable name.
#' @param tvar Character (length 1 or 2), time variable name(s).
#' @param is_quarterly Logical, whether quarterly data.
#' @return An object of class `lwdid_result`.
#' @export
#' @family lwdid-results
new_lwdid_result <- function(
    att = NA_real_,
    se_att = NA_real_,
    t_stat = NA_real_,
    pvalue = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    df_resid = NA_integer_,
    df_inference = NA_integer_,
    nobs = NA_integer_,
    n_treated = NA_integer_,
    n_control = NA_integer_,
    K = NA_integer_,
    tpost1 = NA_integer_,
    depvar = NA_character_,
    rolling = NA_character_,
    vce_type = NULL,
    cluster_var = NULL,
    n_clusters = NULL,
    estimator = NULL,
    method = if (is_staggered) "staggered" else "common_timing",
    alpha = 0.05,
    is_staggered = FALSE,
    controls_used = FALSE,
    controls = character(0),
    include_pretreatment = FALSE,
    control_group = NULL,
    control_group_used = NULL,
    params = numeric(0),
    bse = numeric(0),
    vcov_matrix = NULL,
    resid = numeric(0),
    data = NULL,
    att_by_period = NULL,
    att_pre_treatment = NULL,
    parallel_trends_test = NULL,
    aggregate = NULL,
    cohorts = NULL,
    cohort_sizes = NULL,
    n_never_treated = NULL,
    n_units = NULL,
    n_periods = NULL,
    n_cohorts = NULL,
    att_by_cohort = NULL,
    att_by_cohort_time = NULL,
    att_overall = NULL,
    se_overall = NULL,
    ci_overall_lower = NULL,
    ci_overall_upper = NULL,
    t_stat_overall = NULL,
    pvalue_overall = NULL,
    cohort_effects = NULL,
    event_time_effects = NULL,
    # Cohort aggregation statistics (E5-05.3)
    att_cohort_agg = NULL,
    se_cohort_agg = NULL,
    t_stat_cohort_agg = NULL,
    pvalue_cohort_agg = NULL,
    ci_cohort_agg = NULL,
    cohort_weights = NULL,
    ri_pvalue = NULL,
    ri_seed = NULL,
    rireps = NULL,
    ri_method = NULL,
    ri_valid = NULL,
    ri_failed = NULL,
    ri_error = NULL,
    ri_target = NULL,
    ri_distribution = NULL,
    diagnostics = NULL,
    warning_diagnostics = list(),
    propensity_scores = NULL,
    matched_data = NULL,
    n_matched = NULL,
    match_rate = NULL,
    weights_cv = NULL,
    warnings_log = list(),
    call = NULL,
    lwdid_version = tryCatch(
      as.character(utils::packageVersion("lwdid")),
      error = function(e) "dev"
    ),
    ivar = NA_character_,
    tvar = NA_character_,
    is_quarterly = FALSE
) {
  obj <- list(
    att = att, se_att = se_att, t_stat = t_stat, pvalue = pvalue,
    ci_lower = ci_lower, ci_upper = ci_upper,
    df_resid = df_resid, df_inference = df_inference,
    nobs = nobs, n_treated = n_treated, n_control = n_control,
    K = K, tpost1 = tpost1, depvar = depvar, rolling = rolling,
    vce_type = vce_type, cluster_var = cluster_var, n_clusters = n_clusters,
    cmd = "lwdid", estimator = estimator, method = method, alpha = alpha,
    is_staggered = is_staggered, controls_used = controls_used,
    controls = controls, include_pretreatment = include_pretreatment,
    control_group = control_group, control_group_used = control_group_used,
    params = params, bse = bse, vcov_matrix = vcov_matrix, resid = resid,
    data = data,
    att_by_period = att_by_period, att_pre_treatment = att_pre_treatment,
    parallel_trends_test = parallel_trends_test,
    aggregate = aggregate, cohorts = cohorts, cohort_sizes = cohort_sizes,
    n_never_treated = n_never_treated,
    att_by_cohort = att_by_cohort, att_by_cohort_time = att_by_cohort_time,
    att_overall = att_overall, se_overall = se_overall,
    ci_overall_lower = ci_overall_lower, ci_overall_upper = ci_overall_upper,
    t_stat_overall = t_stat_overall, pvalue_overall = pvalue_overall,
    cohort_effects = cohort_effects, event_time_effects = event_time_effects,
    att_cohort_agg = att_cohort_agg, se_cohort_agg = se_cohort_agg,
    t_stat_cohort_agg = t_stat_cohort_agg,
    pvalue_cohort_agg = pvalue_cohort_agg,
    ci_cohort_agg = ci_cohort_agg,
    n_units = n_units, n_periods = n_periods, n_cohorts = n_cohorts,
    cohort_weights = cohort_weights,
    ri_pvalue = ri_pvalue, ri_seed = ri_seed, rireps = rireps,
    ri_method = ri_method, ri_valid = ri_valid, ri_failed = ri_failed,
    ri_error = ri_error, ri_target = ri_target,
    ri_distribution = ri_distribution,
    diagnostics = diagnostics, warning_diagnostics = warning_diagnostics,
    propensity_scores = propensity_scores, matched_data = matched_data,
    n_matched = n_matched, match_rate = match_rate, weights_cv = weights_cv,
    warnings_log = warnings_log,
    call = call, lwdid_version = lwdid_version,
    ivar = ivar, tvar = tvar, is_quarterly = is_quarterly
  )

  top_level_inference_dist <- .resolve_top_level_inference_dist(estimator)
  if (is.na(top_level_inference_dist) &&
      is.numeric(df_inference) && length(df_inference) == 1L &&
      is.finite(df_inference) && df_inference > 0L) {
    top_level_inference_dist <- "t"
  }
  obj$inference_dist <- top_level_inference_dist

  # --- Compute top-level inference (E5-05.3) ---
  if (is.numeric(att) && length(att) == 1L && is.finite(att) &&
      is.numeric(se_att) && length(se_att) == 1L && is.finite(se_att) &&
      se_att > 0) {
    has_t_stat <- is.numeric(obj$t_stat) && length(obj$t_stat) == 1L &&
      is.finite(obj$t_stat)
    has_pvalue <- is.numeric(obj$pvalue) && length(obj$pvalue) == 1L &&
      is.finite(obj$pvalue)
    has_ci_lower <- is.numeric(obj$ci_lower) && length(obj$ci_lower) == 1L &&
      is.finite(obj$ci_lower)
    has_ci_upper <- is.numeric(obj$ci_upper) && length(obj$ci_upper) == 1L &&
      is.finite(obj$ci_upper)

    if (!has_t_stat) {
      obj$t_stat <- att / se_att
    }

    if (identical(top_level_inference_dist, "normal")) {
      z_crit <- stats::qnorm(1 - alpha / 2)
      if (!has_pvalue) {
        obj$pvalue <- 2 * (1 - stats::pnorm(abs(obj$t_stat)))
      }
      if (!has_ci_lower) {
        obj$ci_lower <- att - z_crit * se_att
      }
      if (!has_ci_upper) {
        obj$ci_upper <- att + z_crit * se_att
      }
    } else if (is.numeric(df_inference) && length(df_inference) == 1L &&
               is.finite(df_inference) && df_inference > 0L) {
      if (!has_pvalue) {
        obj$pvalue <- 2 * stats::pt(
          abs(obj$t_stat), df = df_inference, lower.tail = FALSE
        )
      }
      t_crit <- stats::qt(1 - alpha / 2, df = df_inference)
      if (!has_ci_lower) {
        obj$ci_lower <- att - t_crit * se_att
      }
      if (!has_ci_upper) {
        obj$ci_upper <- att + t_crit * se_att
      }
    }
  }

  # --- Compute count fields (E5-05.3) ---
  obj$n_gr_effects <- if (!is.null(att_by_cohort_time) &&
                          is.data.frame(att_by_cohort_time)) {
    nrow(att_by_cohort_time)
  } else { 0L }
  obj$n_cohort_effects <- if (!is.null(cohort_effects) &&
                              is.list(cohort_effects)) {
    length(cohort_effects)
  } else { 0L }
  obj$n_event_time_effects <- if (!is.null(event_time_effects) &&
                                  is.list(event_time_effects)) {
    length(event_time_effects)
  } else { 0L }

  # --- Split ci_cohort_agg into lower/upper (E5-05.3) ---
  if (!is.null(ci_cohort_agg) && length(ci_cohort_agg) == 2L) {
    obj$ci_lower_cohort_agg <- ci_cohort_agg[1]
    obj$ci_upper_cohort_agg <- ci_cohort_agg[2]
  }

  # --- Default control_group_used (E5-05.3) ---
  if (is.null(obj$control_group_used) && !is.null(obj$control_group)) {
    obj$control_group_used <- obj$control_group
  }

  class(obj) <- "lwdid_result"
  obj
}

#' Validate internal consistency of an lwdid_result object
#'
#' @param x An `lwdid_result` object.
#' @return The object invisibly if valid; throws error otherwise.
#' @keywords internal
validate_lwdid_result <- function(x) {
  tol <- LWDID_NUMERICAL_TOLERANCE

  # 1. Type check: x must inherit "lwdid_result"
  if (!inherits(x, "lwdid_result")) {
    stop_lwdid("Object does not inherit from 'lwdid_result'.",
               class = "lwdid_invalid_parameter",
               param = "x", value = class(x), allowed = "lwdid_result")
  }

  # 2. se_att >= 0 (skip if NA)
  if (!is.na(x$se_att) && x$se_att < 0) {
    stop_lwdid("se_att must be non-negative.",
               class = "lwdid_invalid_parameter",
               param = "se_att", value = x$se_att, allowed = ">= 0")
  }

  # 3. t_stat ~= att / se_att within tolerance (skip if any NA)
  if (!is.na(x$att) && !is.na(x$se_att) && !is.na(x$t_stat)) {
    if (x$se_att != 0) {
      expected_t <- x$att / x$se_att
      if (abs(x$t_stat - expected_t) > tol) {
        stop_lwdid(
          sprintf("t_stat (%.6e) != att/se_att (%.6e)", x$t_stat, expected_t),
          class = "lwdid_invalid_parameter",
          param = "t_stat", value = x$t_stat, allowed = expected_t
        )
      }
    }
    # se_att == 0 allows Inf/NaN t_stat; no check needed
  }

  # 4. df_inference > 0 (skip if NA)
  if (!is.na(x$df_inference) && x$df_inference <= 0) {
    stop_lwdid("df_inference must be positive.",
               class = "lwdid_invalid_parameter",
               param = "df_inference", value = x$df_inference, allowed = "> 0")
  }

  # 5. df_resid > 0 (skip if NA)
  if (!is.na(x$df_resid) && x$df_resid <= 0) {
    stop_lwdid("df_resid must be positive.",
               class = "lwdid_invalid_parameter",
               param = "df_resid", value = x$df_resid, allowed = "> 0")
  }

  # 6. pvalue in [0, 1] (skip if NA)
  if (!is.na(x$pvalue) && (x$pvalue < 0 || x$pvalue > 1)) {
    stop_lwdid("pvalue must be in [0, 1].",
               class = "lwdid_invalid_parameter",
               param = "pvalue", value = x$pvalue, allowed = "[0, 1]")
  }

  # 7. ci_lower <= att <= ci_upper within tolerance (skip if any NA)
  if (!is.na(x$ci_lower) && !is.na(x$att) && !is.na(x$ci_upper)) {
    if (x$ci_lower - tol > x$att || x$att > x$ci_upper + tol) {
      stop_lwdid(
        sprintf("CI [%.6e, %.6e] does not contain att %.6e",
                x$ci_lower, x$ci_upper, x$att),
        class = "lwdid_invalid_parameter",
        param = "ci", value = c(x$ci_lower, x$ci_upper), allowed = "contains att"
      )
    }
  }

  # 8. n_treated >= 0 (skip if NA)
  if (!is.na(x$n_treated) && x$n_treated < 0) {
    stop_lwdid("n_treated must be non-negative.",
               class = "lwdid_invalid_parameter",
               param = "n_treated", value = x$n_treated, allowed = ">= 0")
  }

  # 9. n_control >= 0 (skip if NA)
  if (!is.na(x$n_control) && x$n_control < 0) {
    stop_lwdid("n_control must be non-negative.",
               class = "lwdid_invalid_parameter",
               param = "n_control", value = x$n_control, allowed = ">= 0")
  }

  # 10. nobs >= n_treated + n_control - only for Common Timing (skip if any NA)
  if (!isTRUE(x$is_staggered)) {
    if (!is.na(x$nobs) && !is.na(x$n_treated) && !is.na(x$n_control)) {
      if (x$nobs < x$n_treated + x$n_control) {
        stop_lwdid(
          sprintf("nobs (%d) < n_treated + n_control (%d + %d = %d).",
                  x$nobs, x$n_treated, x$n_control,
                  x$n_treated + x$n_control),
          class = "lwdid_invalid_parameter",
          param = "nobs", value = x$nobs,
          allowed = sprintf(">= %d", x$n_treated + x$n_control)
        )
      }
    }
  }

  # 11. cluster VCE: n_clusters >= 2 (only when vce_type == "cluster")
  if (!is.null(x$vce_type) && identical(x$vce_type, "cluster")) {
    if (!is.null(x$n_clusters) && !is.na(x$n_clusters) && x$n_clusters < 2) {
      stop_lwdid("Cluster VCE requires n_clusters >= 2.",
                 class = "lwdid_invalid_parameter",
                 param = "n_clusters", value = x$n_clusters, allowed = ">= 2")
    }
  }

  # 12. alpha in (0, 1) (skip if NA)
  if (!is.na(x$alpha) && (x$alpha <= 0 || x$alpha >= 1)) {
    stop_lwdid("alpha must be in (0, 1).",
               class = "lwdid_invalid_parameter",
               param = "alpha", value = x$alpha, allowed = "(0, 1)")
  }

  # 13. is_staggered must be TRUE or FALSE (not NA, not NULL)
  if (is.null(x$is_staggered) || is.na(x$is_staggered)) {
    stop_lwdid("is_staggered must be TRUE or FALSE (not NA or NULL).",
               class = "lwdid_invalid_parameter",
               param = "is_staggered", value = as.character(x$is_staggered),
               allowed = "TRUE or FALSE")
  }

  # 14. Staggered consistency: if is_staggered==TRUE, aggregate must be non-NULL
  if (isTRUE(x$is_staggered)) {
    if (is.null(x$aggregate)) {
      stop_lwdid("Staggered result must have non-NULL aggregate.",
                 class = "lwdid_invalid_parameter",
                 param = "aggregate", value = "NULL",
                 allowed = "none/cohort/overall")
    }
  }

  # 15. CT consistency: if is_staggered==FALSE, att_overall must be NULL
  if (isFALSE(x$is_staggered)) {
    if (!is.null(x$att_overall)) {
      stop_lwdid("Common Timing result must not have att_overall.",
                 class = "lwdid_invalid_parameter",
                 param = "att_overall", value = x$att_overall,
                 allowed = "NULL")
    }
  }

  # 16. df_inference <= df_resid (skip if any NA)
  #     non-cluster: stop_lwdid error; cluster: warn_lwdid warning only
  if (!is.na(x$df_inference) && !is.na(x$df_resid)) {
    if (x$df_inference > x$df_resid) {
      is_cluster <- !is.null(x$vce_type) && identical(x$vce_type, "cluster")
      if (is_cluster) {
        warn_lwdid(
          sprintf("df_inference (%s) > df_resid (%s) under cluster VCE.",
                  x$df_inference, x$df_resid),
          class = "lwdid_data",
          detail = "df_inference exceeds df_resid",
          action_taken = "proceeding with cluster-robust inference"
        )
      } else {
        stop_lwdid(
          sprintf("df_inference (%s) must be <= df_resid (%s).",
                  x$df_inference, x$df_resid),
          class = "lwdid_invalid_parameter",
          param = "df_inference", value = x$df_inference,
          allowed = sprintf("<= %s", x$df_resid)
        )
      }
    }
  }

  # 17. K >= 0 (skip if NA)
  if (!is.na(x$K) && x$K < 0) {
    stop_lwdid("K must be non-negative.",
               class = "lwdid_invalid_parameter",
               param = "K", value = x$K, allowed = ">= 0")
  }

  # 18. tpost1 > K (skip if both non-NA)
  if (!is.na(x$tpost1) && !is.na(x$K)) {
    if (x$tpost1 <= x$K) {
      stop_lwdid(
        sprintf("tpost1 (%s) must be > K (%s).", x$tpost1, x$K),
        class = "lwdid_invalid_parameter",
        param = "tpost1", value = x$tpost1,
        allowed = sprintf("> %s", x$K)
      )
    }
  }

  invisible(x)
}

# -- extract_effects --------------------------------------------------------

#' Extract effects from lwdid result as a tidy data.frame
#'
#' @param result An `lwdid_result` object.
#' @param type Character or NULL. One of "gr", "cohort", "overall",
#'   "event_time". If NULL, auto-inferred from `result$aggregate`.
#' @return A data.frame of extracted effects.
#' @export
#' @family lwdid-results
extract_effects <- function(result, type = NULL) {
  # --- Input validation ---
  if (!inherits(result, "lwdid_result")) {
    stop_lwdid(
      "result must be a lwdid_result object.",
      class = "lwdid_invalid_input",
      param = "result", value = class(result),
      allowed = "lwdid_result"
    )
  }

  valid_types <- c("gr", "cohort", "overall", "event_time")

  # --- Auto-infer type from aggregate ---
  if (is.null(type)) {
    agg <- result$aggregate
    if (is.null(agg) || identical(agg, "none")) {
      type <- "gr"
    } else {
      type <- agg
    }
  }

  if (!type %in% valid_types) {
    stop_lwdid(
      sprintf("Invalid type '%s'. Must be one of: %s",
              type, paste(valid_types, collapse = ", ")),
      class = "lwdid_invalid_input",
      param = "type", value = type,
      allowed = paste(valid_types, collapse = ", ")
    )
  }

  # --- Dispatch by type ---
  switch(type,
    "gr" = {
      # Return cohort_time_effects (alias for att_by_cohort_time)
      eff <- result$cohort_time_effects
      if (is.null(eff) && !is.null(result$att_by_cohort_time)) {
        eff <- result$att_by_cohort_time
      }
      if (is.null(eff) || (is.data.frame(eff) && nrow(eff) == 0L)) {
        warn_lwdid(
          "No (g,r) effects available to extract.",
          class = "lwdid_data",
          detail = "cohort_time_effects is NULL or empty",
          action_taken = "returning empty data.frame"
        )
        return(data.frame())
      }
      if (is.data.frame(eff)) return(eff)
      # If list of lists, rbind
      do.call(rbind, lapply(eff, as.data.frame, stringsAsFactors = FALSE))
    },

    "cohort" = {
      eff <- result$cohort_effects
      if (is.null(eff) || length(eff) == 0L) {
        warn_lwdid(
          "No cohort effects available to extract.",
          class = "lwdid_data",
          detail = "cohort_effects is NULL or empty",
          action_taken = "returning empty data.frame"
        )
        return(data.frame())
      }
      # rbind list of cohort effect lists -> data.frame
      rows <- lapply(eff, function(e) {
        data.frame(
          cohort       = as.integer(e$cohort),
          att          = as.numeric(e$att),
          se           = as.numeric(e$se),
          ci_lower     = as.numeric(e$ci_lower),
          ci_upper     = as.numeric(e$ci_upper),
          t_stat       = as.numeric(e$t_stat),
          pvalue       = as.numeric(e$pvalue),
          n_periods    = as.integer(e$n_periods),
          n_units      = as.integer(e$n_units),
          # Defensive extraction for n_control (design doc spec)
          n_control    = if (!is.null(e$n_control)) {
            as.integer(e$n_control)
          } else {
            NA_integer_
          },
          df_resid     = as.integer(e$df_resid),
          df_inference = as.integer(e$df_inference),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    },

    "overall" = {
      # Check for overall_effect list first, then individual fields
      oe <- result$overall_effect
      if (!is.null(oe)) {
        data.frame(
          att          = as.numeric(oe$att),
          se           = as.numeric(oe$se),
          ci_lower     = as.numeric(oe$ci_lower),
          ci_upper     = as.numeric(oe$ci_upper),
          t_stat       = as.numeric(oe$t_stat),
          pvalue       = as.numeric(oe$pvalue),
          n_treated    = as.integer(oe$n_treated %||% NA_integer_),
          n_control    = as.integer(oe$n_control %||% NA_integer_),
          df_resid     = as.integer(oe$df_resid %||% NA_integer_),
          df_inference = as.integer(oe$df_inference %||% NA_integer_),
          stringsAsFactors = FALSE
        )
      } else if (!is.null(result$att_overall)) {
        # Fallback to individual fields
        data.frame(
          att          = as.numeric(result$att_overall),
          se           = as.numeric(result$se_overall %||% NA_real_),
          ci_lower     = as.numeric(result$ci_overall_lower %||% NA_real_),
          ci_upper     = as.numeric(result$ci_overall_upper %||% NA_real_),
          t_stat       = as.numeric(result$t_stat_overall %||% NA_real_),
          pvalue       = as.numeric(result$pvalue_overall %||% NA_real_),
          n_treated    = NA_integer_,
          n_control    = NA_integer_,
          df_resid     = NA_integer_,
          df_inference = NA_integer_,
          stringsAsFactors = FALSE
        )
      } else {
        warn_lwdid(
          "No overall effect available to extract.",
          class = "lwdid_data",
          detail = "overall_effect is NULL",
          action_taken = "returning empty data.frame"
        )
        return(data.frame())
      }
    },

    "event_time" = {
      eff <- result$event_time_effects
      if (is.null(eff) || length(eff) == 0L) {
        warn_lwdid(
          "No event-time effects available to extract.",
          class = "lwdid_data",
          detail = "event_time_effects is NULL or empty",
          action_taken = "returning empty data.frame"
        )
        return(data.frame())
      }
      rows <- lapply(eff, function(e) {
        data.frame(
          event_time   = as.integer(e$event_time),
          att          = as.numeric(e$att),
          se           = as.numeric(e$se),
          ci_lower     = as.numeric(e$ci_lower),
          ci_upper     = as.numeric(e$ci_upper),
          t_stat       = as.numeric(e$t_stat),
          pvalue       = as.numeric(e$pvalue),
          df_inference = as.integer(e$df_inference),
          n_cohorts    = as.integer(e$n_cohorts),
          weight_sum   = as.numeric(e$weight_sum),
          stringsAsFactors = FALSE
        )
      })
      do.call(rbind, rows)
    }
  )
}

# -- Helper functions (shared across S3 methods) ----------------------------

#' @keywords internal
.vce_description <- function(vce_type, cluster_var = NULL) {
  if (is.null(vce_type) || length(vce_type) == 0L) {
    return("OLS (Homoskedastic)")
  }
  mapping <- c(
    "ols"       = "OLS (Homoskedastic)",
    "robust"    = "HC1 (Heteroskedasticity-robust)",
    "hc0"       = "HC0 (White)",
    "hc1"       = "HC1 (Heteroskedasticity-robust)",
    "hc2"       = "HC2 (Bell-McCaffrey)",
    "hc3"       = "HC3 (Small-sample adjusted)",
    "hc4"       = "HC4 (Cribari-Neto)",
    "bootstrap" = "Bootstrap",
    "wild_cluster_bootstrap" = "Wild Cluster Bootstrap"
  )
  key <- tolower(vce_type)
  if (key == "cluster") {
    return(sprintf("Cluster-robust (clustered by %s)",
                   if (!is.null(cluster_var)) cluster_var else "?"))
  }
  if (key %in% names(mapping)) return(unname(mapping[key]))
  vce_type
}

#' @keywords internal
.format_pvalue <- function(p, digits = 4L) {
  if (is.null(p) || length(p) == 0L) return("NA")
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.*f", digits, p)
}

#' @keywords internal
.ri_total_permutations <- function(valid = NULL, failed = NULL,
                                   reps = NULL, distribution = NULL) {
  if (!is.null(distribution)) {
    total <- length(distribution)
    if (total > 0L) {
      return(as.integer(total))
    }
  }

  if (!is.null(reps) && !is.na(reps) && reps > 0L) {
    return(as.integer(reps))
  }

  if (!is.null(valid) && !is.null(failed) &&
      !is.na(valid) && !is.na(failed)) {
    total <- as.integer(valid) + as.integer(failed)
    if (total > 0L) {
      return(total)
    }
  }

  NULL
}

#' @keywords internal
.format_ri_valid_count <- function(valid = NULL, total = NULL) {
  if (is.null(valid) || is.na(valid)) {
    return("?")
  }

  valid_int <- as.integer(valid)
  if (!is.null(total) && !is.na(total) && total > 0L) {
    return(sprintf("%d/%d", valid_int, as.integer(total)))
  }

  as.character(valid_int)
}

#' @keywords internal
.signif_stars <- function(p) {
  if (is.na(p)) return(" ")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.10)  return(".")
  " "
}

#' @keywords internal
.build_diagnostics_summary <- function(diagnostics) {
  summary_list <- list()
  if (!is.null(diagnostics$clustering)) {
    d <- diagnostics$clustering
    summary_list$clustering <- sprintf(
      "n_clusters=%d, effective=%.1f, reliability=%s",
      d$n_clusters, d$effective_clusters, d$reliability_level
    )
  }
  if (!is.null(diagnostics$selection)) {
    d <- diagnostics$selection
    summary_list$selection <- sprintf(
      "attrition=%.1f%%, risk=%s",
      d$attrition_rate * 100, d$selection_risk
    )
  }
  if (!is.null(diagnostics$parallel_trends)) {
    d <- diagnostics$parallel_trends
    # f_test may be nested (lwdid_trends class) or flat
    f_stat <- d$f_stat %||% d$f_test$f_stat
    f_pval <- d$f_pvalue %||% d$f_test$f_pvalue
    if (!is.null(f_stat) && !is.null(f_pval)) {
      summary_list$trends <- sprintf(
        "F=%.3f, p=%s",
        f_stat, .format_pvalue(f_pval)
      )
    }
  }
  summary_list
}

# -- print.lwdid_result -----------------------------------------------------

#' @title Print lwdid result
#' @description Compact output of ATT estimate, SE, CI, and sample info.
#'   Staggered mode shows meta information, top-level ATT summary,
#'   and aggregate-specific effect summaries.
#' @param x lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param ... ignored
#' @return invisible(x)
#' @export
print.lwdid_result <- function(x, digits = 4L, ...) {
  stopifnot(inherits(x, "lwdid_result"))
  has_wcb <- identical(tolower(x$vce_type %||% ""), "wild_cluster_bootstrap") &&
    !is.null(x$wcb_details)

  cat("\nLocal Wald DID Estimation\n")
  cat("Lee-Wooldridge DiD Estimation\n")
  cat(sprintf("Method: %s | Estimator: %s | Rolling: %s\n",
              x$method, x$estimator, x$rolling))
  vce_desc <- .vce_description(x$vce_type, x$cluster_var)
  cat(sprintf("VCE: %s | Dep.var: %s\n", vce_desc, x$depvar))
  if (!is.null(x$vce_type) && identical(tolower(x$vce_type), "cluster") &&
      !is.null(x$n_clusters) && !is.na(x$n_clusters)) {
    cat(sprintf("Number of clusters: %d\n", as.integer(x$n_clusters)))
  }
  cat(paste(rep("-", 50), collapse = ""), "\n")

  if (isTRUE(x$is_staggered)) {
    # --- Meta information (E5-05.5 step 1) ---
    agg_label <- x$aggregate %||% "none"
    cat(sprintf("Staggered DID | Aggregate: %s\n", agg_label))

    # Control group with auto-switch hint
    if (!is.null(x$control_group_used)) {
      cg_str <- x$control_group_used
      if (!is.null(x$control_group) &&
          !identical(x$control_group_used, x$control_group)) {
        cg_str <- sprintf("%s [auto-switched from %s]",
                          cg_str, x$control_group)
      }
      cat(sprintf("Control group: %s\n", cg_str))
    }

    # Counts: units, periods, cohorts, NT
    n_u <- x$n_units %||% NA_integer_
    n_p <- x$n_periods %||% NA_integer_
    n_c <- x$n_cohorts %||% {
      if (!is.null(x$cohorts)) length(x$cohorts) else NA_integer_
    }
    n_nt <- x$n_never_treated %||% NA_integer_
    count_parts <- character(0)
    if (!is.na(n_u))  count_parts <- c(count_parts, sprintf("Units: %d", n_u))
    if (!is.na(n_p))  count_parts <- c(count_parts, sprintf("Periods: %d", n_p))
    if (!is.na(n_c))  count_parts <- c(count_parts, sprintf("Cohorts: %d", n_c))
    if (!is.na(n_nt)) count_parts <- c(count_parts, sprintf("NT: %d", n_nt))
    if (length(count_parts) > 0L) {
      cat(paste(count_parts, collapse = " | "), "\n")
    }
    if (!is.na(n_nt)) {
      cat(sprintf("N Never-treated: %d\n", n_nt))
    }

    # (g,r) effects count
    n_gr <- x$n_gr_effects %||% {
      if (!is.null(x$att_by_cohort_time) &&
          is.data.frame(x$att_by_cohort_time)) {
        nrow(x$att_by_cohort_time)
      } else {
        0L
      }
    }
    if (n_gr > 0L) cat(sprintf("(g,r) effects: %d\n", n_gr))

    cat(paste(rep("-", 50), collapse = ""), "\n")

    # --- Top-level ATT/SE summary (E5-05.5 step 2) ---
    has_att <- !is.null(x$att) && !is.na(x$att)
    if (has_att) {
      # Use tau_omega label when overall ATT exists while keeping source ASCII-safe.
      has_overall <- !is.null(x$overall_effect) || !is.null(x$att_overall)
      att_label <- if (has_overall) "ATT(tau_omega)" else "ATT"
      cat(sprintf("%s = %.*f\n", att_label, digits, x$att))
      if (!is.null(x$se_att) && !is.na(x$se_att)) {
        se_label <- if (has_wcb) "WCB SE" else "SE"
        cat(sprintf("%-8s = %.*f\n", se_label, digits, x$se_att))
      }
      if (!is.null(x$t_stat) && !is.na(x$t_stat)) {
        t_label <- if (has_wcb) "t-stat(original)" else "t-stat"
        cat(sprintf("%s = %.*f\n", t_label, digits, x$t_stat))
      }
      if (!is.null(x$pvalue) && !is.na(x$pvalue)) {
        p_label <- if (has_wcb) "WCB p-value" else "p-value"
        cat(sprintf("%s = %s %s\n",
                    p_label,
                    .format_pvalue(x$pvalue, digits),
                    .signif_stars(x$pvalue)))
      }
      if (!has_wcb && !is.null(x$df_inference) && !is.na(x$df_inference)) {
        cat(sprintf("df      = %d\n", as.integer(x$df_inference)))
      }
      if (!is.null(x$ci_lower) && !is.na(x$ci_lower) &&
          !is.null(x$ci_upper) && !is.na(x$ci_upper)) {
        ci_pct <- as.integer((1 - (x$alpha %||% 0.05)) * 100)
        ci_label <- if (has_wcb) "WCB CI" else "CI"
        cat(sprintf("%s [%d%%] = [%.*f, %.*f]\n",
                    ci_label, ci_pct, digits, x$ci_lower, digits, x$ci_upper))
      }
      if (has_wcb) {
        cat(sprintf(
          "Weight type: %s | Requested bootstrap: %d | Actual bootstrap: %d\n",
          x$wcb_details$weight_type,
          as.integer(x$wcb_details$requested_n_bootstrap),
          as.integer(x$wcb_details$actual_n_bootstrap)
        ))
        cat(sprintf(
          "Wild Cluster Bootstrap | Clusters: %d | Method: %s | Restricted model: %s | Full enumeration: %s\n",
          as.integer(x$wcb_details$n_clusters),
          x$wcb_details$method %||% "native",
          x$wcb_details$restricted_model %||% "with_controls",
          if (isTRUE(x$wcb_details$full_enumeration)) "TRUE" else "FALSE"
        ))
      }
    } else {
      cat("ATT     = NA (no valid top-level estimate)\n")
      cat(sprintf("Aggregate: %s\n", x$aggregate %||% "none"))
    }

    # --- Aggregate-specific summaries (E5-05.5 step 3) ---
    if (identical(agg_label, "cohort")) {
      # Cohort-aggregated ATT
      if (!is.null(x$att_cohort_agg) && !is.na(x$att_cohort_agg)) {
        cat(paste(rep("-", 50), collapse = ""), "\n")
        cat("Cohort-Aggregated Effect:\n")
        cat(sprintf("  ATT(cohort) = %.*f", digits, x$att_cohort_agg))
        if (!is.null(x$se_cohort_agg) && !is.na(x$se_cohort_agg)) {
          cat(sprintf("  SE = %.*f", digits, x$se_cohort_agg))
        }
        if (!is.null(x$pvalue_cohort_agg) && !is.na(x$pvalue_cohort_agg)) {
          cat(sprintf("  p = %s %s",
                      .format_pvalue(x$pvalue_cohort_agg, digits),
                      .signif_stars(x$pvalue_cohort_agg)))
        }
        cat("\n")
      }
      # Cohort effects table
      if (!is.null(x$cohort_effects) && length(x$cohort_effects) > 0L) {
        cat(paste(rep("-", 50), collapse = ""), "\n")
        cat("Cohort Effects:\n")
        tbl <- do.call(rbind, lapply(x$cohort_effects, function(e) {
          data.frame(
            cohort = as.integer(e$cohort),
            att    = round(as.numeric(e$att), digits),
            se     = round(as.numeric(e$se), digits),
            pvalue = as.numeric(e$pvalue),
            stringsAsFactors = FALSE
          )
        }))
        tbl$sig <- vapply(tbl$pvalue, .signif_stars, character(1))
        tbl$pvalue <- vapply(tbl$pvalue, .format_pvalue, character(1),
                             digits = digits)
        print(tbl, row.names = FALSE, right = FALSE)
      }

    } else if (identical(agg_label, "overall")) {
      cat(paste(rep("-", 50), collapse = ""), "\n")
      cat("Overall Effect:\n")
      oe <- x$overall_effect
      if (!is.null(oe)) {
        cat(sprintf("  ATT(overall) = %.*f\n", digits, oe$att))
        if (!is.null(oe$se) && !is.na(oe$se)) {
          cat(sprintf("  SE           = %.*f\n", digits, oe$se))
        }
        if (!is.null(oe$t_stat) && !is.na(oe$t_stat)) {
          cat(sprintf("  t-stat       = %.*f\n", digits, oe$t_stat))
        }
        if (!is.null(oe$pvalue) && !is.na(oe$pvalue)) {
          cat(sprintf("  p-value      = %s %s\n",
                      .format_pvalue(oe$pvalue, digits),
                      .signif_stars(oe$pvalue)))
        }
        if (!is.null(oe$ci_lower) && !is.na(oe$ci_lower) &&
            !is.null(oe$ci_upper) && !is.na(oe$ci_upper)) {
          ci_pct <- as.integer((1 - (x$alpha %||% 0.05)) * 100)
          cat(sprintf("  CI [%d%%]     = [%.*f, %.*f]\n",
                      ci_pct, digits, oe$ci_lower, digits, oe$ci_upper))
        }
      } else if (!is.null(x$att_overall) && !is.na(x$att_overall)) {
        cat(sprintf("  ATT(overall) = %.*f\n", digits, x$att_overall))
        if (!is.null(x$se_overall) && !is.na(x$se_overall)) {
          cat(sprintf("  SE           = %.*f\n", digits, x$se_overall))
        }
      }
      # Hint about cohort effects
      if (!is.null(x$cohort_effects) && length(x$cohort_effects) > 0L) {
        cat("  Use extract_effects(x, \"cohort\") for cohort-level details.\n")
      }

    } else if (identical(agg_label, "event_time")) {
      if (!is.null(x$event_time_effects) &&
          length(x$event_time_effects) > 0L) {
        cat(paste(rep("-", 50), collapse = ""), "\n")
        cat("Event-Time Effects:\n")
        tbl <- do.call(rbind, lapply(x$event_time_effects, function(e) {
          data.frame(
            event_time = as.integer(e$event_time),
            att        = round(as.numeric(e$att), digits),
            se         = round(as.numeric(e$se), digits),
            pvalue     = as.numeric(e$pvalue),
            n_cohorts  = as.integer(e$n_cohorts),
            stringsAsFactors = FALSE
          )
        }))
        tbl$sig <- vapply(tbl$pvalue, .signif_stars, character(1))
        tbl$pvalue <- vapply(tbl$pvalue, .format_pvalue, character(1),
                             digits = digits)
        print(tbl, row.names = FALSE, right = FALSE)
      }
    }
    # none: no additional aggregate summary needed

    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(sprintf("Pre-treatment: %s\n",
                if (isTRUE(x$include_pretreatment)) "TRUE" else "FALSE"))

  } else {
    # --- Common Timing path ---
    has_wcb <- identical(tolower(x$vce_type %||% ""), "wild_cluster_bootstrap") &&
      !is.null(x$wcb_details)
    cat(sprintf("ATT     = %.*f\n", digits, x$att))
    se_label <- if (has_wcb) "WCB SE" else "SE"
    cat(sprintf("%-8s = %.*f\n", se_label, digits, x$se_att))
    if (!is.null(x$t_stat) && !is.na(x$t_stat)) {
      t_label <- if (has_wcb) "t-stat(original)" else "t-stat"
      cat(sprintf("%s = %.*f\n", t_label, digits, x$t_stat))
    }
    p_label <- if (has_wcb) "WCB p-value" else "p-value"
    cat(sprintf("%s = %s\n", p_label, .format_pvalue(x$pvalue, digits)))
    if (!has_wcb) {
      cat(sprintf("df      = %d\n", x$df_inference))
    }
    ci_pct <- as.integer((1 - x$alpha) * 100)
    ci_label <- if (has_wcb) "WCB CI" else "CI"
    cat(sprintf("%s [%d%%] = [%.*f, %.*f]\n",
                ci_label, ci_pct, digits, x$ci_lower, digits, x$ci_upper))
    if (has_wcb) {
      cat(sprintf(
        "Weight type: %s | Requested bootstrap: %d | Actual bootstrap: %d\n",
        x$wcb_details$weight_type,
        as.integer(x$wcb_details$requested_n_bootstrap),
        as.integer(x$wcb_details$actual_n_bootstrap)
      ))
      cat(sprintf(
        "Wild Cluster Bootstrap | Clusters: %d | Method: %s | Restricted model: %s | Full enumeration: %s\n",
        as.integer(x$wcb_details$n_clusters),
        x$wcb_details$method %||% "native",
        x$wcb_details$restricted_model %||% "with_controls",
        if (isTRUE(x$wcb_details$full_enumeration)) "TRUE" else "FALSE"
      ))
    }

    if (!is.null(x$att_by_period) &&
        is.data.frame(x$att_by_period) &&
        nrow(x$att_by_period) > 0L) {
      cat("\nPeriod-by-period effects available in results$att_by_period\n")
      n_period_rows <- nrow(x$att_by_period)
      print(utils::head(x$att_by_period, 5L), digits = digits, row.names = FALSE)
      if (n_period_rows > 5L) {
        cat(sprintf("  ... (%d more periods)\n", n_period_rows - 5L))
      }
    }
  }
  cat(paste(rep("-", 50), collapse = ""), "\n")

  # Sample size line
  nobs_val <- x$nobs %||% NA_integer_
  nt_val <- x$n_treated %||% NA_integer_
  nc_val <- x$n_control %||% NA_integer_
  if (!is.na(nobs_val)) {
    cat(sprintf("N = %d", nobs_val))
    if (!is.na(nt_val)) cat(sprintf(" | N_treated = %d", nt_val))
    if (!is.na(nc_val)) cat(sprintf(" | N_control = %d", nc_val))
    if (!is.null(x$K) && !is.na(x$K)) {
      cat(sprintf(" | K = %d", x$K))
    }
    if (!is.null(x$tpost1) && !is.na(x$tpost1)) {
      cat(sprintf(" | tpost1 = %s", x$tpost1))
    }
    cat("\n")
  }

  # RI section (E7-06.4.1, E7-06.4.2)
  if (!is.null(x$ri_pvalue) || !is.null(x$ri_error)) {
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("Randomization Inference:\n")
    if (!is.null(x$ri_pvalue)) {
      cat(sprintf("  RI p-value = %s\n", .format_pvalue(x$ri_pvalue, digits)))
    }
    details <- character(0)
    ri_valid <- x$ri_valid %||% x$ri_n_valid
    ri_failed <- x$ri_failed %||% x$ri_n_failed
    ri_total <- .ri_total_permutations(
      valid = ri_valid,
      failed = ri_failed,
      reps = x$rireps,
      distribution = x$ri_distribution
    )
    if (!is.null(x$ri_method)) {
      details <- c(details, sprintf("method=%s", x$ri_method))
    }
    if (!is.null(x$rireps)) {
      details <- c(details, sprintf("reps=%d", as.integer(x$rireps)))
    }
    if (!is.null(x$ri_seed)) {
      details <- c(details, sprintf("seed=%d", as.integer(x$ri_seed)))
    }
    if (!is.null(ri_valid)) {
      details <- c(
        details,
        sprintf("valid=%s", .format_ri_valid_count(ri_valid, ri_total))
      )
    }
    if (!is.null(ri_failed)) {
      details <- c(details, sprintf("failed=%d", as.integer(ri_failed)))
    }
    if (length(details) > 0L) {
      cat(sprintf("  %s\n", paste(details, collapse = " | ")))
    }
    if (!is.null(x$ri_target)) {
      cat(sprintf("  Target: %s\n", x$ri_target))
    }
    # E7-06.4.2: RI error display
    if (!is.null(x$ri_error)) {
      cat(sprintf("  RI Error: %s\n", x$ri_error))
    }
  }

  # Pre-treatment effects summary (E7-06.4.3)
  if (isTRUE(x$include_pretreatment) && !is.null(x$att_pre_treatment) &&
      is.data.frame(x$att_pre_treatment) && nrow(x$att_pre_treatment) > 0L) {
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("Pre-treatment Effects:\n")
    pre_df <- x$att_pre_treatment
    n_pre <- nrow(pre_df)
    n_sig_05 <- sum(pre_df$pvalue < 0.05, na.rm = TRUE)
    max_abs <- max(abs(pre_df$att), na.rm = TRUE)
    cat(sprintf("  Periods: %d | Significant (p<0.05): %d | Max|ATT|: %.*f\n",
                n_pre, n_sig_05, digits, max_abs))
  }

  # Parallel trends test (E7-06.4.4)
  if (!is.null(x$parallel_trends_test)) {
    pt <- x$parallel_trends_test
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("Parallel Trends Test (Joint Test):\n")
    stat_val <- pt$joint_stat %||% pt$f_stat
    p_val <- pt$joint_pvalue %||% pt$f_pvalue
    reject <- pt$reject_null %||% (p_val < (x$alpha %||% 0.05))
    if (!is.null(stat_val) && !is.null(p_val)) {
      cat(sprintf("  Statistic = %.*f | p-value = %s\n",
                  digits, stat_val, .format_pvalue(p_val, digits)))
      if (isTRUE(reject)) {
        cat("  => Reject H0: Evidence against parallel trends\n")
      } else {
        cat("  => Fail to reject H0: No evidence against parallel trends\n")
      }
    }
  }

  # exclude_pre_periods (E7-06.4.5)
  if (!is.null(x$exclude_pre_periods) && x$exclude_pre_periods > 0L) {
    cat(sprintf("Excluded pre-periods: %d\n", x$exclude_pre_periods))
  }

  if (length(x$warnings_log) > 0L) {
    cat(sprintf("\n[%d warning(s) issued]\n", length(x$warnings_log)))
  }

  if (isTRUE(x$is_staggered)) {
    cat("Hint: Use coef(x, type='all') or results$att_by_cohort_time for (g,r)-level effects; ")
    cat("use plot_event_study(x) or plot(x) for event study.\n")
  }

  invisible(x)
}

# -- summary.lwdid_result ---------------------------------------------------

#' @title lwdid result summary
#' @description Generate a full regression summary including period effects,
#'   cohort effects, parallel trends test, and diagnostics.
#' @param object lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param ... ignored
#' @return summary.lwdid_result object (invisibly)
#' @export
summary.lwdid_result <- function(object, digits = 4L, ...) {
  stopifnot(inherits(object, "lwdid_result"))

  out <- list()
  out$call       <- object$metadata$call
  out$method     <- object$method
  out$estimator  <- object$estimator
  out$rolling    <- object$rolling
  out$vce_type   <- object$vce_type
  out$depvar     <- object$depvar
  out$nobs       <- object$nobs
  out$n_treated  <- object$n_treated
  out$n_control  <- object$n_control
  out$df_inference <- object$df_inference
  out$alpha      <- object$alpha
  out$K          <- object$K
  out$tpost1     <- object$tpost1
  out$params     <- object$params
  out$bse        <- object$bse
  out$wcb_details <- object$wcb_details

  # VCE human-readable description
  out$vce_description <- .vce_description(object$vce_type, object$cluster_var)

  # Overall ATT coefficient table
  out$coefficients <- data.frame(
    Estimate  = object$att,
    Std.Error = object$se_att,
    t.value   = object$t_stat,
    Pr.t      = object$pvalue,
    CI.lower  = object$ci_lower,
    CI.upper  = object$ci_upper,
    row.names = "ATT"
  )

  if (length(object$params) > 0L) {
    params <- object$params
    ses <- object$bse
    if (length(ses) != length(params)) {
      ses <- rep(NA_real_, length(params))
      names(ses) <- names(params)
    }
    t_vals <- rep(NA_real_, length(params))
    valid_se <- is.finite(ses) & ses > 0
    t_vals[valid_se] <- params[valid_se] / ses[valid_se]
    p_vals <- rep(NA_real_, length(params))
    if (!is.null(object$df_inference) && !is.na(object$df_inference)) {
      p_vals[is.finite(t_vals)] <- 2 * stats::pt(
        -abs(t_vals[is.finite(t_vals)]),
        df = object$df_inference
      )
    }
    out$regression_coefficients <- data.frame(
      Estimate = as.numeric(params),
      Std.Error = as.numeric(ses),
      t.value = as.numeric(t_vals),
      Pr.t = as.numeric(p_vals),
      row.names = names(params),
      stringsAsFactors = FALSE
    )
  }

  # Period-specific effects
  if (!is.null(object$att_by_period)) {
    out$period_effects <- object$att_by_period
    out$n_period_effects <- nrow(object$att_by_period)
  }

  # Staggered-specific
  if (isTRUE(object$is_staggered)) {
    out$is_staggered <- TRUE
    out$cohort_effects     <- object$att_by_cohort
    out$event_time_effects <- object$event_time_effects
    out$cohort_weights     <- object$cohort_weights
    out$cohort_list        <- if (!is.null(object$att_by_cohort)) {
      unique(object$att_by_cohort$cohort)
    } else if (!is.null(object$cohorts)) {
      object$cohorts
    } else {
      integer(0)
    }
    out$control_group_used <- object$control_group_used
    out$control_group      <- object$control_group
    out$control_group_auto_switched <- !is.null(object$control_group) &&
      !is.null(object$control_group_used) &&
      !identical(object$control_group_used, object$control_group)
    out$has_overall_effect <- !is.null(object$att) && !is.na(object$att)
    out$n_never_treated    <- object$n_never_treated
    out$aggregate          <- object$aggregate %||% "none"
    if (!is.null(object$cohort_sample_sizes)) {
      out$cohort_sample_sizes <- object$cohort_sample_sizes
    }
    # Cohort effects with n_units and n_periods (AC-30)
    if (!is.null(object$cohort_effects) && is.list(object$cohort_effects)) {
      ce_names <- names(object$cohort_effects)
      ce_rows <- Map(function(e, nm) {
        cohort_val <- e$cohort
        if (is.null(cohort_val) && !is.null(nm) && nzchar(nm)) {
          cohort_val <- suppressWarnings(as.integer(nm))
        }
        data.frame(
          cohort   = as.integer(cohort_val %||% NA_integer_),
          att      = as.numeric(e$att),
          se       = as.numeric(e$se),
          n_units  = as.integer(e$n_units %||% NA_integer_),
          n_periods = as.integer(e$n_periods %||% NA_integer_),
          stringsAsFactors = FALSE
        )
      }, object$cohort_effects, ce_names %||% rep("", length(object$cohort_effects)))
      out$cohort_effects_detail <- do.call(rbind, ce_rows)
    }
  } else {
    out$is_staggered <- FALSE
  }

  # Parallel trends test (with reject conclusion)
  if (!is.null(object$parallel_trends_test)) {
    pt <- object$parallel_trends_test
    # Normalize field names: test_parallel_trends_joint returns
    # joint_stat/joint_pvalue/joint_df1/joint_df2/reject_null
    out$parallel_trends <- list(
      f_stat  = pt$joint_stat  %||% pt$f_stat,
      f_pvalue = pt$joint_pvalue %||% pt$f_pvalue,
      df1     = pt$joint_df1   %||% pt$df1 %||% 0L,
      df2     = pt$joint_df2   %||% pt$df2 %||% 0L,
      reject  = {
        p_val <- pt$joint_pvalue %||% pt$f_pvalue
        pt$reject_null %||%
          ((!is.null(p_val) && !is.na(p_val)) &&
             p_val < (object$alpha %||% 0.05))
      }
    )
  }

  # RI details
  if (!is.null(object$ri_pvalue) || !is.null(object$ri_error)) {
    ri_valid <- object$ri_valid %||% object$ri_n_valid
    ri_failed <- object$ri_failed %||% object$ri_n_failed
    out$ri_pvalue         <- object$ri_pvalue
    out$ri_n_permutations <- if (!is.null(object$ri_distribution)) {
      length(object$ri_distribution)
    } else {
      object$rireps %||% 0L
    }
    out$ri_method         <- object$ri_method
    out$ri_seed           <- object$ri_seed
    out$ri_n_valid        <- ri_valid
    out$ri_n_failed       <- ri_failed
    out$ri_error          <- object$ri_error
    out$ri_target         <- object$ri_target
    out$rireps            <- object$rireps
    # RI distribution summary (E7-06.4.6)
    if (!is.null(object$ri_distribution) &&
        length(object$ri_distribution) > 0L) {
      ri_dist <- object$ri_distribution[is.finite(object$ri_distribution)]
      if (length(ri_dist) > 0L) {
        out$ri_dist_summary <- list(
          mean   = mean(ri_dist),
          median = stats::median(ri_dist),
          sd     = stats::sd(ri_dist)
        )
      }
    }
  }

  # Exclude pre-periods
  if (!is.null(object$exclude_pre_periods) &&
      object$exclude_pre_periods > 0L) {
    out$exclude_pre_periods <- object$exclude_pre_periods
  }

  # Diagnostics summary
  if (!is.null(object$diagnostics)) {
    out$diagnostics_summary <- .build_diagnostics_summary(object$diagnostics)
  }

  # Pre-treatment dynamics
  if (!is.null(object$att_pre_treatment)) {
    out$att_pre_treatment <- object$att_pre_treatment
    out$include_pretreatment <- isTRUE(object$include_pretreatment)
  }

  class(out) <- "summary.lwdid_result"
  out
}

# -- print.summary.lwdid_result ---------------------------------------------

#' @title Print lwdid summary
#' @description Full regression summary output with printCoefmat, cohort
#'   weights with N, pre-treatment dynamics with anchor marking, parallel
#'   trends reject/fail-to-reject conclusions, and RI details.
#' @param x summary.lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param signif.stars logical, show significance stars (default TRUE)
#' @param ... ignored
#' @return invisible(x)
#' @export
print.summary.lwdid_result <- function(x, digits = 4L,
                                        signif.stars = TRUE, ...) {
  cat("\nLocal Wald DID Estimation Summary\n")
  cat("Lee-Wooldridge DiD Estimation\n")
  if (!is.null(x$call)) cat("Call:", deparse(x$call), "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Model info
  cat(sprintf("Method: %s | Estimator: %s | Rolling: %s\n",
              x$method, x$estimator, x$rolling))
  cat(sprintf("VCE: %s\n", x$vce_description))
  if (!is.null(x$wcb_details)) {
    cat(sprintf("Dep.var: %s | df: WCB\n", x$depvar))
  } else {
    cat(sprintf("Dep.var: %s | df: %d\n", x$depvar, x$df_inference))
  }
  if (!is.null(x$K)) {
    cat(sprintf("K: %d", x$K))
    if (!is.null(x$tpost1)) cat(sprintf(" | tpost1: %s", x$tpost1))
    cat("\n")
  }
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # Staggered header
  if (isTRUE(x$is_staggered)) {
    cat(sprintf("\nStaggered DID | Cohorts: %s\n",
                paste(x$cohort_list, collapse = ", ")))
    if (!is.null(x$control_group_used)) {
      cat(sprintf("Control group: %s", x$control_group_used))
      if (isTRUE(x$control_group_auto_switched) && !is.null(x$control_group)) {
        cat(sprintf(" [auto-switched from %s]", x$control_group))
      }
      cat("\n")
    }
    if (!is.null(x$n_never_treated)) {
      cat(sprintf("Number of Never Treated Units: %d\n", x$n_never_treated))
    }
    if (!is.null(x$aggregate)) {
      cat(sprintf("Aggregation: %s\n", x$aggregate))
    }
  }

  # Overall ATT
  if (isTRUE(x$is_staggered) && isTRUE(x$has_overall_effect)) {
    cat("\nOverall Weighted Effect (tau_omega):\n")
  } else {
    cat("\nOverall ATT:\n")
  }
  coef_mat <- as.matrix(x$coefficients[, 1:4])
  colnames(coef_mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  stats::printCoefmat(coef_mat, digits = digits,
                      signif.stars = signif.stars,
                      P.values = TRUE, has.Pvalue = TRUE,
                      na.print = "")
  ci_pct <- as.integer((1 - x$alpha) * 100)
  cat(sprintf("  CI [%d%%]: [%.*f, %.*f]\n",
              ci_pct,
              digits, x$coefficients$CI.lower,
              digits, x$coefficients$CI.upper))

  if (!is.null(x$wcb_details)) {
    cat("\nWild Cluster Bootstrap:\n")
    cat(sprintf("  Weight type: %s\n", x$wcb_details$weight_type))
    cat(sprintf(
      "  Requested bootstrap: %d | Actual bootstrap: %d\n",
      as.integer(x$wcb_details$requested_n_bootstrap),
      as.integer(x$wcb_details$actual_n_bootstrap)
    ))
    cat(sprintf(
      "  Clusters: %d | Method: %s | Restricted model: %s\n",
      as.integer(x$wcb_details$n_clusters),
      x$wcb_details$method %||% "native",
      x$wcb_details$restricted_model %||% "with_controls"
    ))
    cat(sprintf(
      "  Full enumeration: %s\n",
      if (isTRUE(x$wcb_details$full_enumeration)) "TRUE" else "FALSE"
    ))
  }

  if (!is.null(x$regression_coefficients)) {
    cat("\nRegression Coefficients:\n")
    reg_mat <- as.matrix(x$regression_coefficients[, 1:4])
    colnames(reg_mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    stats::printCoefmat(reg_mat, digits = digits,
                        signif.stars = signif.stars,
                        P.values = TRUE, has.Pvalue = TRUE,
                        na.print = "")
  }

  # Cohort weights (with N)
  if (!is.null(x$cohort_weights)) {
    cat("\n--- Cohort Weights ---\n")
    if (!is.null(x$cohort_sample_sizes)) {
      for (nm in names(x$cohort_weights)) {
        cat(sprintf("  %s: w=%.3f, N=%d\n",
                    nm, x$cohort_weights[nm],
                    x$cohort_sample_sizes[nm]))
      }
    } else {
      cat(paste(sprintf("  %s: w=%.3f",
                        names(x$cohort_weights),
                        x$cohort_weights), collapse = "\n"))
      cat("\n")
    }
  }

  # Cohort-specific effects
  if (!is.null(x$cohort_effects)) {
    cat("\n--- Cohort Effects ---\n")
    print(x$cohort_effects, digits = digits, row.names = FALSE)
  }

  # Event-time effects
  if (!is.null(x$event_time_effects)) {
    cat("\n--- Event-Time Effects ---\n")
    print(x$event_time_effects, digits = digits, row.names = FALSE)
  }

  # Period-specific effects (truncate >5 rows)
  if (!is.null(x$period_effects)) {
    cat("\n--- Period-Specific Effects ---\n")
    print(x$period_effects, digits = digits, row.names = FALSE)
    if (nrow(x$period_effects) > 5L) {
      cat(sprintf("  ... (%d more rows)\n", nrow(x$period_effects) - 5L))
    }
  }

  # Pre-treatment dynamics (with anchor point marking)
  if (!is.null(x$att_pre_treatment) && isTRUE(x$include_pretreatment)) {
    cat("\n--- Pre-treatment Dynamics ---\n")
    pre_df <- x$att_pre_treatment
    anchor_idx <- which(pre_df$event_time == -1)
    if (length(anchor_idx) > 0L) {
      pre_display <- pre_df
      pre_display$note <- ""
      pre_display$note[anchor_idx] <- "<- anchor (by construction)"
      print(pre_display, digits = digits, row.names = FALSE)
    } else {
      print(pre_df, digits = digits, row.names = FALSE)
      cat("  e=-1: ATT=0 by construction (anchor point)\n")
    }
  }

  # Parallel trends test (with reject conclusion)
  if (!is.null(x$parallel_trends)) {
    pt <- x$parallel_trends
    cat("\n--- Parallel Trends Test ---\n")
    f_stat_val <- pt$f_stat %||% NA_real_
    f_pval_val <- pt$f_pvalue %||% NA_real_
    df1_val <- pt$df1 %||% 0L
    df2_val <- pt$df2 %||% 0L
    if (!is.na(f_stat_val) && !is.na(f_pval_val)) {
      cat(sprintf("Joint F-test: F(%d,%d) = %.3f, p = %s\n",
                  df1_val, df2_val,
                  f_stat_val,
                  .format_pvalue(f_pval_val, digits)))
    }
    alpha_val <- x$alpha %||% 0.05
    if (isTRUE(pt$reject)) {
      cat(sprintf("  => Reject H0 at alpha=%.2f: ", alpha_val))
      cat("Evidence against parallel trends\n")
    } else {
      cat(sprintf("  => Fail to reject H0 at alpha=%.2f: ", alpha_val))
      cat("No evidence against parallel trends\n")
    }
  }

  # RI details
  if (!is.null(x$ri_pvalue) || !is.null(x$ri_error)) {
    cat(sprintf("\nRandomization Inference (RI):\n"))
    if (!is.null(x$ri_pvalue)) {
      cat(sprintf("  RI p-value = %s\n", .format_pvalue(x$ri_pvalue, digits)))
    }
    details <- c()
    ri_total <- .ri_total_permutations(
      valid = x$ri_n_valid,
      failed = x$ri_n_failed,
      reps = x$ri_n_permutations
    )
    if (!is.null(x$ri_method)) {
      details <- c(details, sprintf("method=%s", x$ri_method))
    }
    if (!is.null(x$rireps)) {
      details <- c(details, sprintf("reps=%d", x$rireps))
    }
    if (!is.null(x$ri_seed)) {
      details <- c(details, sprintf("seed=%d", x$ri_seed))
    }
    if (!is.null(x$ri_n_valid)) {
      details <- c(
        details,
        sprintf("valid=%s", .format_ri_valid_count(x$ri_n_valid, ri_total))
      )
    }
    if (!is.null(x$ri_n_failed)) {
      details <- c(details, sprintf("failed=%d", x$ri_n_failed))
    }
    if (length(details) > 0L) {
      cat(sprintf("  %s\n", paste(details, collapse = " | ")))
    }
    if (!is.null(x$ri_target)) {
      cat(sprintf("  Target: %s\n", x$ri_target))
    }
    # RI distribution summary (E7-06.4.6)
    if (!is.null(x$ri_dist_summary)) {
      ds <- x$ri_dist_summary
      cat(sprintf("  Distribution: mean=%.*f, median=%.*f, sd=%.*f\n",
                  digits, ds$mean, digits, ds$median, digits, ds$sd))
    }
    # RI error
    if (!is.null(x$ri_error)) {
      cat(sprintf("  RI Error: %s\n", x$ri_error))
    }
  }

  # Diagnostics summary
  if (!is.null(x$diagnostics_summary) &&
      length(x$diagnostics_summary) > 0L) {
    cat("\n--- Diagnostics ---\n")
    for (nm in names(x$diagnostics_summary)) {
      cat(sprintf("  %s: %s\n", nm, x$diagnostics_summary[[nm]]))
    }
  }

  # Exclude pre-periods
  if (!is.null(x$exclude_pre_periods) && x$exclude_pre_periods > 0L) {
    cat(sprintf("\nExcluded pre-periods: %d\n", x$exclude_pre_periods))
  }

  # Sample info
  cat(paste(rep("-", 60), collapse = ""), "\n")
  cat(sprintf("N = %d | N_treated = %d | N_control = %d",
              x$nobs, x$n_treated, x$n_control))
  if (!is.null(x$K)) cat(sprintf(" | K = %d", x$K))
  if (!is.null(x$tpost1)) cat(sprintf(" | tpost1 = %s", x$tpost1))
  cat("\n")

  # Staggered usage hints
  if (isTRUE(x$is_staggered)) {
    cat("\nHint: Use coef(x, type='all') for (g,r)-level effects, plot(x) for event study.\n")
  }

  invisible(x)
}

# -- coef.lwdid_result ------------------------------------------------------

#' @title Extract ATT coefficients
#' @description Extract ATT coefficients from lwdid_result object.
#' @param object lwdid_result object
#' @param type character, granularity: "overall" (default), "by_period",
#'   "by_cohort", "all"
#' @param ... ignored
#' @return named numeric vector
#' @export
coef.lwdid_result <- function(object, type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "all"))

  switch(type,
    "overall" = {
      out <- object$att
      names(out) <- "ATT"
      out
    },
    "by_period" = {
      if (is.null(object$att_by_period))
        stop("No period-specific results available. Use aggregate='none' or include att_by_period data.",
             call. = FALSE)
      out <- object$att_by_period$att
      names(out) <- as.character(object$att_by_period$period)
      out
    },
    "by_cohort" = {
      if (is.null(object$att_by_cohort))
        stop("No cohort-specific results (Staggered mode only).",
             call. = FALSE)
      out <- object$att_by_cohort$att
      names(out) <- as.character(object$att_by_cohort$cohort)
      out
    },
    "all" = {
      if (!is.null(object$att_by_cohort_time)) {
        out <- object$att_by_cohort_time$att
        names(out) <- sprintf("g%s.r%s",
          object$att_by_cohort_time$cohort,
          object$att_by_cohort_time$period)
      } else if (!is.null(object$att_by_period)) {
        out <- object$att_by_period$att
        names(out) <- as.character(object$att_by_period$period)
      } else {
        out <- object$att
        names(out) <- "ATT"
      }
      out
    }
  )
}

# -- confint.lwdid_result ---------------------------------------------------

#' @title Confidence intervals
#' @description Compute confidence intervals for ATT estimates.
#' @param object lwdid_result object
#' @param parm character vector or NULL, select specific parameters
#' @param level numeric, confidence level (default 0.95)
#' @param type character, granularity (matches coef): "overall", "by_period",
#'   "by_cohort", "all"
#' @param ... ignored
#' @return matrix with lower and upper bounds
#' @export
confint.lwdid_result <- function(object, parm = NULL, level = 0.95,
                                  type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  if (level <= 0 || level >= 1)
    stop("level must be in (0, 1)", call. = FALSE)
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "all"))

  info <- switch(type,
    "overall" = {
      list(est = object$att, se = object$se_att, rn = "ATT")
    },
    "by_period" = {
      if (is.null(object$att_by_period))
        stop("No period-specific results for confint.", call. = FALSE)
      list(est = object$att_by_period$att,
           se  = object$att_by_period$se,
           rn  = as.character(object$att_by_period$period))
    },
    "by_cohort" = {
      if (is.null(object$att_by_cohort))
        stop("No cohort-specific results (Staggered mode only).",
             call. = FALSE)
      list(est = object$att_by_cohort$att,
           se  = object$att_by_cohort$se,
           rn  = as.character(object$att_by_cohort$cohort))
    },
    "all" = {
      if (!is.null(object$att_by_cohort_time)) {
        list(est = object$att_by_cohort_time$att,
             se  = object$att_by_cohort_time$se,
             rn  = sprintf("g%s.r%s",
                           object$att_by_cohort_time$cohort,
                           object$att_by_cohort_time$period))
      } else if (!is.null(object$att_by_period)) {
        list(est = object$att_by_period$att,
             se  = object$att_by_period$se,
             rn  = as.character(object$att_by_period$period))
      } else {
        list(est = object$att, se = object$se_att, rn = "ATT")
      }
    }
  )

  alpha <- 1 - level
  df <- object$df_inference
  t_crit <- stats::qt(1 - alpha / 2, df = df)

  ci_lower <- info$est - t_crit * info$se
  ci_upper <- info$est + t_crit * info$se

  ci <- cbind(ci_lower, ci_upper)
  pct <- sprintf("%.1f %%", c(alpha / 2, 1 - alpha / 2) * 100)
  colnames(ci) <- pct
  rownames(ci) <- info$rn

  if (!is.null(parm)) {
    if (is.character(parm)) {
      ci <- ci[parm, , drop = FALSE]
    } else if (is.numeric(parm)) {
      ci <- ci[parm, , drop = FALSE]
    }
  }

  ci
}

# -- vcov.lwdid_result ------------------------------------------------------

#' @title Variance-covariance matrix
#' @description Extract variance-covariance matrix for ATT estimates.
#' @param object lwdid_result object
#' @param type character, granularity: "overall" (1x1), "by_period",
#'   "by_cohort", "full" (complete OLS parameter VCE)
#' @param ... ignored
#' @return variance-covariance matrix
#' @export
vcov.lwdid_result <- function(object, type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "full"))

  switch(type,
    "overall" = {
      V <- matrix(object$se_att^2, nrow = 1L, ncol = 1L)
      rownames(V) <- colnames(V) <- "ATT"
      V
    },
    "by_period" = {
      if (is.null(object$att_by_period))
        stop("No period-specific results for vcov.", call. = FALSE)
      if (!is.null(object$vcov_att_periods)) {
        V <- object$vcov_att_periods
      } else {
        ses <- object$att_by_period$se
        V <- diag(ses^2, nrow = length(ses))
      }
      pnames <- as.character(object$att_by_period$period)
      rownames(V) <- colnames(V) <- pnames
      V
    },
    "by_cohort" = {
      if (is.null(object$att_by_cohort))
        stop("No cohort-specific results for vcov.", call. = FALSE)
      if (!is.null(object$vcov_att_cohorts)) {
        V <- object$vcov_att_cohorts
      } else {
        ses <- object$att_by_cohort$se
        V <- diag(ses^2, nrow = length(ses))
      }
      cnames <- as.character(object$att_by_cohort$cohort)
      rownames(V) <- colnames(V) <- cnames
      V
    },
    "full" = {
      if (is.null(object$vcov_full))
        stop("Full OLS VCE matrix not available. ",
             "Set store_vcov_full=TRUE in lwdid() call.",
             call. = FALSE)
      object$vcov_full
    }
  )
}

# -- nobs.lwdid_result ------------------------------------------------------

#' @title Number of observations
#' @description Extract the number of observations used in the lwdid model.
#'   Ensures compatibility with modelsummary and broom::glance().
#' @param object lwdid_result object
#' @param ... ignored
#' @return integer, number of observations
#' @importFrom stats nobs
#' @exportS3Method stats::nobs
nobs.lwdid_result <- function(object, ...) {
  stopifnot(inherits(object, "lwdid_result"))
  object$nobs
}


# plot_event_study generic - S3 method in plot_event_study.R
#' @title Event Study Plot
#' @description Generic function for event study visualization.
#' @param x An object to plot.
#' @param ... Additional arguments passed to methods.
#' @return A ggplot object or data.
#' @export
plot_event_study <- function(x, ...) UseMethod("plot_event_study")

# plot.lwdid_result() - Moved to plot.R (Story 10.5)
