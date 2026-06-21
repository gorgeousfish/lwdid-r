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

.lwdid_effect_inference_stats <- function(att, se, df, alpha,
                                          inference_dist = NULL) {
  att <- as.numeric(att)
  se <- as.numeric(se)
  n <- length(att)

  if (length(se) == 1L && n != 1L) {
    se <- rep(se, n)
  }
  if (length(se) != n) {
    stop("Effect standard errors must match the number of ATT estimates.",
         call. = FALSE)
  }

  if (is.null(df)) {
    df <- rep(NA_real_, n)
  } else {
    df <- as.numeric(df)
    if (length(df) == 1L && n != 1L) {
      df <- rep(df, n)
    }
    if (length(df) != n) {
      stop("Inference degrees of freedom must match the number of ATT estimates.",
           call. = FALSE)
    }
  }

  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).", call. = FALSE)
  }

  inference_dist <- inference_dist %||% NA_character_
  use_normal <- identical(inference_dist, "normal")

  t_stat <- rep(NA_real_, n)
  valid_se <- is.finite(se) & se > 0
  t_stat[valid_se] <- att[valid_se] / se[valid_se]

  zero_se <- is.finite(se) & se == 0 & is.finite(att)
  t_stat[zero_se & att > 0] <- Inf
  t_stat[zero_se & att < 0] <- -Inf
  t_stat[zero_se & att == 0] <- NaN

  pvalue <- rep(NA_real_, n)
  finite_t <- is.finite(t_stat)
  if (use_normal) {
    pvalue[finite_t] <- 2 * stats::pnorm(
      abs(t_stat[finite_t]), lower.tail = FALSE
    )
  } else {
    valid_df <- finite_t & is.finite(df) & df > 0
    pvalue[valid_df] <- 2 * stats::pt(
      abs(t_stat[valid_df]), df = df[valid_df], lower.tail = FALSE
    )
  }
  pvalue[is.infinite(t_stat)] <- 0
  pvalue[is.nan(t_stat)] <- NaN

  if (use_normal) {
    crit <- rep(stats::qnorm(1 - alpha / 2), n)
  } else {
    crit <- rep(NA_real_, n)
    valid_df <- is.finite(df) & df > 0
    crit[valid_df] <- stats::qt(1 - alpha / 2, df = df[valid_df])
  }

  data.frame(
    se = se,
    t_stat = t_stat,
    pvalue = pvalue,
    ci_lower = att - crit * se,
    ci_upper = att + crit * se,
    stringsAsFactors = FALSE
  )
}

.lwdid_result_inference_dist <- function(object) {
  inference_dist <- object$inference_dist %||%
    .resolve_top_level_inference_dist(object$estimator)
  if (!is.character(inference_dist) || length(inference_dist) != 1L ||
      is.na(inference_dist) || !(inference_dist %in% c("t", "normal"))) {
    inference_dist <- "t"
  }
  inference_dist
}

.lwdid_effect_rows_with_supplied_vcov <- function(df, object, type) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    return(df)
  }
  if (!"att" %in% names(df)) {
    stop("Effect rows must contain an att column.", call. = FALSE)
  }
  if (identical(type, "event_time")) {
    df <- .event_time_rows_with_default_vcov_metadata(df)
  }

  supplied <- switch(type,
    "by_period" = !is.null(object$vcov_att_periods),
    "by_cohort" = !is.null(object$vcov_att_cohorts),
    "event_time" = !is.null(object$vcov_att_event_time),
    FALSE
  )
  if (!isTRUE(supplied)) {
    return(df)
  }

  V <- switch(type,
    "by_period" = vcov(object, type = "by_period"),
    "by_cohort" = vcov(object, type = "by_cohort"),
    "event_time" = .validate_lwdid_effect_vcov(
      object$vcov_att_event_time,
      nrow(df),
      "vcov_att_event_time",
      "event-time effects"
    ),
    stop("Unsupported effect type for supplied VCE alignment.", call. = FALSE)
  )

  alpha <- object$alpha %||% 0.05
  stats <- .lwdid_effect_inference_stats(
    att = df$att,
    se = sqrt(diag(V)),
    df = df$df_inference %||% object$df_inference,
    alpha = alpha,
    inference_dist = .lwdid_result_inference_dist(object)
  )

  df$se <- stats$se
  df$t_stat <- stats$t_stat
  df$pvalue <- stats$pvalue
  df$ci_lower <- stats$ci_lower
  df$ci_upper <- stats$ci_upper
  if (identical(type, "event_time")) {
    covariance_assumption <- attr(V, "covariance_assumption", exact = TRUE)
    if (is.null(covariance_assumption)) {
      covariance_assumption <- "provided_joint_event_time_covariance"
    }
    se_aggregation <- attr(V, "se_aggregation", exact = TRUE)
    if (is.null(se_aggregation)) {
      se_aggregation <- "provided_event_time_vcov_diagonal"
    }
    df$se_aggregation <- as.character(se_aggregation)[1L]
    df$covariance_assumption <- as.character(covariance_assumption)[1L]
  }
  df
}

.event_time_se_aggregation <- function(value = NULL) {
  value <- value %||% "diagonal_weighted_cohort_se"
  value <- as.character(value)[1L]
  if (is.na(value) || !nzchar(value)) {
    value <- "diagonal_weighted_cohort_se"
  }
  value
}

.event_time_covariance_assumption <- function(value = NULL) {
  value <- value %||% "zero_cross_cohort_covariance"
  value <- as.character(value)[1L]
  if (is.na(value) || !nzchar(value)) {
    value <- "zero_cross_cohort_covariance"
  }
  value
}

.event_time_rows_with_default_vcov_metadata <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0L) {
    return(df)
  }

  if (!"se_aggregation" %in% names(df)) {
    df$se_aggregation <- "diagonal_weighted_cohort_se"
  } else {
    se_aggregation <- as.character(df$se_aggregation)
    missing_se <- is.na(se_aggregation) | !nzchar(se_aggregation)
    se_aggregation[missing_se] <- "diagonal_weighted_cohort_se"
    df$se_aggregation <- se_aggregation
  }

  if (!"covariance_assumption" %in% names(df)) {
    df$covariance_assumption <- "zero_cross_cohort_covariance"
  } else {
    covariance_assumption <- as.character(df$covariance_assumption)
    missing_cov <- is.na(covariance_assumption) | !nzchar(covariance_assumption)
    covariance_assumption[missing_cov] <- "zero_cross_cohort_covariance"
    df$covariance_assumption <- covariance_assumption
  }

  df
}

.event_time_metadata_values <- function(value, default) {
  value <- unique(as.character(value))
  value <- value[!is.na(value) & nzchar(value)]
  if (length(value) == 0L) {
    default
  } else {
    value
  }
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
#' @param event_time_omissions List of omitted event-time WATT diagnostics.
#' @param event_time_omission_summary data.frame with `reason` and `n` columns
#'   summarizing omitted event-time WATT rows.
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
#' @param skipped_pairs List of skipped staggered `(g,r)` pair diagnostics.
#' @param skipped_summary data.frame with `reason` and `n` columns summarizing
#'   skipped staggered `(g,r)` pairs.
#' @param ri_pvalue Numeric or NULL, RI p-value.
#' @param ri_seed Integer or NULL.
#' @param rireps Integer or NULL.
#' @param ri_method Character or NULL.
#' @param ri_valid Integer or NULL.
#' @param ri_failed Integer or NULL.
#' @param ri_error Character or NULL.
#' @param ri_target Character or NULL.
#' @param ri_observed_stat Numeric or NULL, observed statistic used by RI.
#' @param ri_estimator Character or NULL, estimator used by RI refits.
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
    event_time_omissions = list(),
    event_time_omission_summary = data.frame(
      reason = character(0),
      n = integer(0),
      stringsAsFactors = FALSE
    ),
    # Cohort aggregation statistics (E5-05.3)
    att_cohort_agg = NULL,
    se_cohort_agg = NULL,
    t_stat_cohort_agg = NULL,
    pvalue_cohort_agg = NULL,
    ci_cohort_agg = NULL,
    cohort_weights = NULL,
    skipped_pairs = list(),
    skipped_summary = data.frame(
      reason = character(0),
      n = integer(0),
      stringsAsFactors = FALSE
    ),
    ri_pvalue = NULL,
    ri_seed = NULL,
    rireps = NULL,
    ri_method = NULL,
    ri_valid = NULL,
    ri_failed = NULL,
    ri_error = NULL,
    ri_target = NULL,
    ri_observed_stat = NULL,
    ri_estimator = NULL,
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
    cohort_sample_sizes = cohort_sizes,
    n_never_treated = n_never_treated,
    att_by_cohort = att_by_cohort, att_by_cohort_time = att_by_cohort_time,
    att_overall = att_overall, se_overall = se_overall,
    ci_overall_lower = ci_overall_lower, ci_overall_upper = ci_overall_upper,
    t_stat_overall = t_stat_overall, pvalue_overall = pvalue_overall,
    cohort_effects = cohort_effects, event_time_effects = event_time_effects,
    event_time_omissions = event_time_omissions,
    event_time_omission_summary = event_time_omission_summary,
    att_cohort_agg = att_cohort_agg, se_cohort_agg = se_cohort_agg,
    t_stat_cohort_agg = t_stat_cohort_agg,
    pvalue_cohort_agg = pvalue_cohort_agg,
    ci_cohort_agg = ci_cohort_agg,
    n_units = n_units, n_periods = n_periods, n_cohorts = n_cohorts,
    cohort_weights = cohort_weights,
    skipped_pairs = skipped_pairs,
    skipped_summary = skipped_summary,
    ri_pvalue = ri_pvalue, ri_seed = ri_seed, rireps = rireps,
    ri_method = ri_method, ri_valid = ri_valid, ri_failed = ri_failed,
    ri_error = ri_error, ri_target = ri_target,
    ri_observed_stat = ri_observed_stat,
    ri_estimator = ri_estimator,
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

  inferred_design_counts <- .infer_lwdid_design_counts(
    data = data,
    ivar = ivar,
    tvar = tvar,
    cohorts = cohorts,
    att_by_cohort = att_by_cohort,
    cohort_effects = cohort_effects
  )
  if (is.null(obj$n_units)) {
    obj$n_units <- inferred_design_counts$n_units
  }
  if (is.null(obj$n_periods)) {
    obj$n_periods <- inferred_design_counts$n_periods
  }
  if (is.null(obj$n_cohorts)) {
    obj$n_cohorts <- inferred_design_counts$n_cohorts
  }

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
                                  is.data.frame(event_time_effects)) {
    nrow(event_time_effects)
  } else if (!is.null(event_time_effects) &&
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
#'   "event_time", or "event_time_contributions". If NULL, auto-inferred from
#'   `result$aggregate`.
#' @return A data.frame of extracted effects. Event-time extractions include
#'   `se_aggregation` and `covariance_assumption` when the fitted object
#'   carries those metadata. For non-RA staggered fits estimated with
#'   propensity diagnostics, event-time rows also carry overlap summaries such
#'   as `max_weight_cv` and `weighted_weight_cv`, plus event-level support
#'   counts such as `min_n_treated` and `min_n_control` when available.
#'   If `vcov_att_event_time` is supplied, event-time rows use its diagonal
#'   for uncertainty and report its covariance metadata instead of the default
#'   diagonal cohort-SE boundary.
#'   Event-time contribution extractions carry the same event-level metadata
#'   and any available cell-level propensity diagnostics on each cohort
#'   contribution row, including bootstrap replicate counts and success rates
#'   when the contributing non-RA cells were estimated with bootstrap standard
#'   errors.
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

  valid_types <- c(
    "gr", "cohort", "overall", "event_time",
    "event_time_contributions"
  )

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
        return(data.frame(
          event_time = integer(0),
          att = numeric(0),
          se = numeric(0),
          ci_lower = numeric(0),
          ci_upper = numeric(0),
          t_stat = numeric(0),
          pvalue = numeric(0),
          df_inference = integer(0),
          n_cohorts = integer(0),
          weight_sum = numeric(0),
          inference_dist = character(0),
          se_aggregation = character(0),
          covariance_assumption = character(0)
        ))
      }
      optional_event_fields <- c(
        "max_weight_cv", "weighted_weight_cv", "min_ps", "max_ps",
        "total_trimmed", "min_n_treated", "max_n_treated",
        "min_n_control", "max_n_control", "min_bootstrap_success_rate",
        "min_bootstrap_reps_valid", "max_bootstrap_reps_failed"
      )
      if (is.data.frame(eff)) {
        return(.lwdid_effect_rows_with_supplied_vcov(
          eff,
          result,
          "event_time"
        ))
      }
      integer_event_fields <- c(
        "total_trimmed", "min_n_treated", "max_n_treated",
        "min_n_control", "max_n_control", "min_bootstrap_reps_valid",
        "max_bootstrap_reps_failed"
      )
      present_event_fields <- optional_event_fields[
        vapply(optional_event_fields, function(field) {
          any(vapply(eff, function(e) !is.null(e[[field]]), logical(1L)))
        }, logical(1L))
      ]
      rows <- lapply(eff, function(e) {
        row <- data.frame(
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
          inference_dist = as.character(
            e$inference_dist %||% result$inference_dist %||% NA_character_
          ),
          se_aggregation = .event_time_se_aggregation(e$se_aggregation),
          covariance_assumption = .event_time_covariance_assumption(
            e$covariance_assumption
          ),
          stringsAsFactors = FALSE
        )
        for (field in present_event_fields) {
          value <- e[[field]]
          if (is.null(value) || length(value) == 0L) {
            row[[field]] <- if (field %in% integer_event_fields) {
              NA_integer_
            } else {
              NA_real_
            }
          } else if (field %in% integer_event_fields) {
            row[[field]] <- as.integer(value[1L])
          } else {
            row[[field]] <- as.numeric(value[1L])
          }
        }
        row
      })
      out <- do.call(rbind, rows)
      .lwdid_effect_rows_with_supplied_vcov(out, result, "event_time")
    },

    "event_time_contributions" = {
      eff <- result$event_time_effects
      if (is.null(eff) || length(eff) == 0L) {
        warn_lwdid(
          "No event-time effects available to extract contributions from.",
          class = "lwdid_data",
          detail = "event_time_effects is NULL or empty",
          action_taken = "returning empty data.frame"
        )
        return(data.frame(
          event_time = integer(0),
          event_att = numeric(0),
          event_se = numeric(0),
          n_cohorts = integer(0),
          weight_sum = numeric(0),
          cohort = integer(0),
          weight = numeric(0),
          contribution_att = numeric(0),
          contribution_se = numeric(0),
          se_aggregation = character(0),
          covariance_assumption = character(0)
        ))
      }
      if (is.data.frame(eff)) {
        stop_lwdid(
          paste0(
            "Event-time contribution extraction requires list-based ",
            "event_time_effects with cohort_contributions."
          ),
          class = "lwdid_invalid_input",
          param = "event_time_effects",
          value = "data.frame",
          allowed = "list elements containing cohort_contributions"
        )
      }

      optional_contribution_fields <- c(
        "n", "n_treated", "n_control", "ps_min", "ps_max", "ps_mean",
        "n_trimmed", "weights_cv", "n_matched", "match_rate", "n_dropped",
        "bootstrap_reps_requested", "bootstrap_reps_valid",
        "bootstrap_reps_failed", "bootstrap_success_rate"
      )
      integer_contribution_fields <- c(
        "n", "n_treated", "n_control", "n_trimmed", "n_matched", "n_dropped",
        "bootstrap_reps_requested", "bootstrap_reps_valid",
        "bootstrap_reps_failed"
      )
      present_contribution_fields <- optional_contribution_fields[
        vapply(optional_contribution_fields, function(field) {
          any(vapply(eff, function(e) {
            contributions <- e$cohort_contributions
            if (is.null(contributions) || length(contributions) == 0L) {
              return(FALSE)
            }
            any(vapply(
              contributions,
              function(cn) !is.null(cn[[field]]),
              logical(1L)
            ))
          }, logical(1L)))
        }, logical(1L))
      ]

      rows <- lapply(eff, function(e) {
        contributions <- e$cohort_contributions
        if (is.null(contributions) || length(contributions) == 0L) {
          return(NULL)
        }
        do.call(rbind, lapply(contributions, function(cn) {
          row <- data.frame(
            event_time = as.integer(e$event_time),
            event_att = as.numeric(e$att),
            event_se = as.numeric(e$se),
            n_cohorts = as.integer(e$n_cohorts),
            weight_sum = as.numeric(e$weight_sum),
            cohort = as.integer(cn$cohort),
            weight = as.numeric(cn$weight),
            contribution_att = as.numeric(cn$att),
            contribution_se = as.numeric(cn$se),
            se_aggregation = .event_time_se_aggregation(e$se_aggregation),
            covariance_assumption = .event_time_covariance_assumption(
              e$covariance_assumption
            ),
            stringsAsFactors = FALSE
          )
          for (field in present_contribution_fields) {
            value <- cn[[field]]
            if (is.null(value) || length(value) == 0L) {
              row[[field]] <- if (field %in% integer_contribution_fields) {
                NA_integer_
              } else {
                NA_real_
              }
            } else if (field %in% integer_contribution_fields) {
              row[[field]] <- as.integer(value[1L])
            } else {
              row[[field]] <- as.numeric(value[1L])
            }
          }
          row
        }))
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows) == 0L) {
        warn_lwdid(
          "No event-time cohort contributions are available to extract.",
          class = "lwdid_data",
          detail = "cohort_contributions is empty for every event time",
          action_taken = "returning empty data.frame"
        )
        return(data.frame())
      }
      out <- do.call(rbind, rows)
      if (!is.null(result$vcov_att_event_time)) {
        V <- .validate_lwdid_effect_vcov(
          result$vcov_att_event_time,
          length(eff),
          "vcov_att_event_time",
          "event-time effects"
        )
        covariance_assumption <- attr(V, "covariance_assumption", exact = TRUE)
        if (is.null(covariance_assumption)) {
          covariance_assumption <- "provided_joint_event_time_covariance"
        }
        se_aggregation <- attr(V, "se_aggregation", exact = TRUE)
        if (is.null(se_aggregation)) {
          se_aggregation <- "provided_event_time_vcov_diagonal"
        }
        event_se <- sqrt(diag(V))
        event_names <- vapply(
          eff,
          function(e) as.character(as.integer(e$event_time)),
          character(1L)
        )
        out$event_se <- unname(event_se[match(
          as.character(out$event_time), event_names
        )])
        out$se_aggregation <- as.character(se_aggregation)[1L]
        out$covariance_assumption <- as.character(covariance_assumption)[1L]
      }
      out
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
    "analytical" = "Analytical (asymptotic normal)",
    "abadie_imbens" = "Abadie-Imbens matching SE",
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
  if (!is.null(diagnostics$propensity)) {
    d <- diagnostics$propensity
    if (is.data.frame(d)) {
      ps_min <- suppressWarnings(min(d$ps_min, na.rm = TRUE))
      ps_max <- suppressWarnings(max(d$ps_max, na.rm = TRUE))
      max_cv <- suppressWarnings(max(d$weights_cv, na.rm = TRUE))
      cv_part <- if (is.finite(max_cv)) {
        sprintf(", max_weight_cv=%.2f", max_cv)
      } else {
        ""
      }
      summary_list$propensity <- sprintf(
        "cells=%d, ps_range=[%.3f, %.3f]%s",
        nrow(d), ps_min, ps_max, cv_part
      )
    } else if (is.list(d)) {
      cv_part <- if (!is.null(d$weights_cv) && is.finite(d$weights_cv)) {
        sprintf(", weight_cv=%.2f", d$weights_cv)
      } else {
        ""
      }
      summary_list$propensity <- sprintf(
        "n=%d, ps_range=[%.3f, %.3f]%s",
        as.integer(d$n %||% NA_integer_),
        as.numeric(d$ps_min %||% NA_real_),
        as.numeric(d$ps_max %||% NA_real_),
        cv_part
      )
    }
  }
  summary_list
}

#' @keywords internal
.event_time_extreme <- function(df, value_col, fun) {
  if (!value_col %in% names(df) || !"event_time" %in% names(df)) {
    return(list(value = NA_real_, event_time = NA_integer_))
  }
  values <- suppressWarnings(as.numeric(df[[value_col]]))
  finite <- is.finite(values)
  if (!any(finite)) {
    return(list(value = NA_real_, event_time = NA_integer_))
  }
  finite_values <- values[finite]
  target <- if (identical(fun, "min")) {
    min(finite_values)
  } else {
    max(finite_values)
  }
  idx <- which(finite & values == target)[1L]
  list(
    value = target,
    event_time = as.integer(df$event_time[[idx]])
  )
}

#' @keywords internal
.summarize_event_time_support_metadata <- function(result) {
  eff <- result$event_time_effects
  if (is.null(eff) || length(eff) == 0L) {
    return(NULL)
  }
  event_df <- if (is.data.frame(eff)) {
    eff
  } else {
    extract_effects(result, type = "event_time")
  }
  if (!is.data.frame(event_df) || nrow(event_df) == 0L) {
    return(NULL)
  }

  max_cv <- .event_time_extreme(event_df, "max_weight_cv", "max")
  max_weighted_cv <- .event_time_extreme(
    event_df, "weighted_weight_cv", "max"
  )
  min_treated <- .event_time_extreme(event_df, "min_n_treated", "min")
  min_control <- .event_time_extreme(event_df, "min_n_control", "min")

  has_any_metadata <- any(vapply(
    list(max_cv, max_weighted_cv, min_treated, min_control),
    function(x) is.finite(x$value),
    logical(1L)
  ))
  if (!has_any_metadata) {
    return(NULL)
  }

  cue_parts <- character(0L)
  if (is.finite(max_cv$value) && max_cv$value > 2) {
    cue_parts <- c(cue_parts, sprintf("max weight CV=%.2f", max_cv$value))
  }
  if (is.finite(min_treated$value) && min_treated$value <= 5) {
    cue_parts <- c(
      cue_parts,
      sprintf("min treated cell N=%d", as.integer(min_treated$value))
    )
  }

  list(
    n_event_time_rows = as.integer(nrow(event_df)),
    max_weight_cv = as.numeric(max_cv$value),
    max_weight_cv_event_time = as.integer(max_cv$event_time),
    max_weighted_weight_cv = as.numeric(max_weighted_cv$value),
    max_weighted_weight_cv_event_time =
      as.integer(max_weighted_cv$event_time),
    min_n_treated = as.integer(min_treated$value),
    min_n_treated_event_time = as.integer(min_treated$event_time),
    min_n_control = as.integer(min_control$value),
    min_n_control_event_time = as.integer(min_control$event_time),
    cue = if (length(cue_parts) > 0L) {
      paste(cue_parts, collapse = ", ")
    } else {
      NA_character_
    }
  )
}

#' @keywords internal
.summarize_event_time_bootstrap_metadata <- function(result) {
  eff <- result$event_time_effects
  if (is.null(eff) || length(eff) == 0L) {
    return(NULL)
  }
  event_df <- if (is.data.frame(eff)) {
    eff
  } else {
    extract_effects(result, type = "event_time")
  }
  if (!is.data.frame(event_df) || nrow(event_df) == 0L) {
    return(NULL)
  }

  min_success <- .event_time_extreme(
    event_df, "min_bootstrap_success_rate", "min"
  )
  min_valid <- .event_time_extreme(
    event_df, "min_bootstrap_reps_valid", "min"
  )
  max_failed <- .event_time_extreme(
    event_df, "max_bootstrap_reps_failed", "max"
  )

  has_any_metadata <- any(vapply(
    list(min_success, min_valid, max_failed),
    function(x) is.finite(x$value),
    logical(1L)
  ))
  if (!has_any_metadata) {
    return(NULL)
  }

  cue_parts <- character(0L)
  if (is.finite(min_success$value)) {
    cue_parts <- c(
      cue_parts,
      sprintf("min bootstrap success=%.2f", min_success$value)
    )
  }
  if (is.finite(min_valid$value)) {
    cue_parts <- c(
      cue_parts,
      sprintf("min valid reps=%d", as.integer(min_valid$value))
    )
  }
  if (is.finite(max_failed$value)) {
    cue_parts <- c(
      cue_parts,
      sprintf("max failed reps=%d", as.integer(max_failed$value))
    )
  }

  list(
    n_event_time_rows = as.integer(nrow(event_df)),
    min_bootstrap_success_rate = as.numeric(min_success$value),
    min_bootstrap_success_rate_event_time =
      as.integer(min_success$event_time),
    min_bootstrap_reps_valid = as.integer(min_valid$value),
    min_bootstrap_reps_valid_event_time = as.integer(min_valid$event_time),
    max_bootstrap_reps_failed = as.integer(max_failed$value),
    max_bootstrap_reps_failed_event_time =
      as.integer(max_failed$event_time),
    cue = if (length(cue_parts) > 0L) {
      paste(cue_parts, collapse = ", ")
    } else {
      NA_character_
    }
  )
}

#' @keywords internal
.empty_skipped_summary <- function() {
  data.frame(
    reason = character(0),
    n = integer(0),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.validate_skipped_diagnostics <- function(skipped_summary,
                                          skipped_pairs = list()) {
  if (is.null(skipped_pairs)) {
    skipped_pairs <- list()
  }
  if (!is.list(skipped_pairs)) {
    stop("skipped_pairs must be a list.", call. = FALSE)
  }

  if (is.null(skipped_summary)) {
    skipped_summary <- .empty_skipped_summary()
  }
  if (!is.data.frame(skipped_summary)) {
    stop("skipped_summary must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("reason", "n")
  if (!all(required_cols %in% names(skipped_summary))) {
    stop("skipped_summary must contain reason and n columns.", call. = FALSE)
  }

  out <- skipped_summary[, required_cols, drop = FALSE]
  out$reason <- as.character(out$reason)
  if (any(is.na(out$reason) | !nzchar(out$reason))) {
    stop("skipped_summary$reason must be non-missing and non-empty.",
         call. = FALSE)
  }
  n_raw <- out$n
  if (!is.numeric(n_raw)) {
    stop("skipped_summary$n must contain non-negative integer counts.",
         call. = FALSE)
  }
  n_numeric <- as.numeric(n_raw)
  if (any(!is.finite(n_numeric) | n_numeric < 0 | n_numeric != floor(n_numeric))) {
    stop("skipped_summary$n must contain non-negative integer counts.",
         call. = FALSE)
  }
  out$n <- as.integer(n_numeric)
  if (anyDuplicated(out$reason)) {
    stop("skipped_summary must have one row per skip reason.", call. = FALSE)
  }

  total_skipped <- sum(out$n)
  if (length(skipped_pairs) > 0L && nrow(out) == 0L) {
    stop("skipped_summary is required when skipped_pairs is non-empty.",
         call. = FALSE)
  }
  if (length(skipped_pairs) > 0L && nrow(out) > 0L &&
      total_skipped != length(skipped_pairs)) {
    stop("skipped_summary counts must match skipped_pairs length.",
         call. = FALSE)
  }

  list(
    skipped_pairs = skipped_pairs,
    skipped_summary = out,
    n_skipped_pairs = as.integer(total_skipped)
  )
}

#' @keywords internal
.validate_event_time_omissions <- function(omission_summary,
                                           omissions = list()) {
  if (is.null(omissions)) {
    omissions <- list()
  }
  if (!is.list(omissions)) {
    stop("event_time_omissions must be a list.", call. = FALSE)
  }

  if (is.null(omission_summary)) {
    omission_summary <- .empty_skipped_summary()
  }
  if (!is.data.frame(omission_summary)) {
    stop("event_time_omission_summary must be a data.frame.", call. = FALSE)
  }
  required_cols <- c("reason", "n")
  if (!all(required_cols %in% names(omission_summary))) {
    stop("event_time_omission_summary must contain reason and n columns.",
         call. = FALSE)
  }

  out <- omission_summary[, required_cols, drop = FALSE]
  out$reason <- as.character(out$reason)
  if (any(is.na(out$reason) | !nzchar(out$reason))) {
    stop("event_time_omission_summary$reason must be non-missing and non-empty.",
         call. = FALSE)
  }
  n_raw <- out$n
  if (!is.numeric(n_raw)) {
    stop("event_time_omission_summary$n must contain non-negative integer counts.",
         call. = FALSE)
  }
  n_numeric <- as.numeric(n_raw)
  if (any(!is.finite(n_numeric) | n_numeric < 0 | n_numeric != floor(n_numeric))) {
    stop("event_time_omission_summary$n must contain non-negative integer counts.",
         call. = FALSE)
  }
  out$n <- as.integer(n_numeric)
  if (anyDuplicated(out$reason)) {
    stop("event_time_omission_summary must have one row per omission reason.",
         call. = FALSE)
  }

  total_omitted <- sum(out$n)
  if (length(omissions) > 0L && nrow(out) == 0L) {
    stop(
      "event_time_omission_summary is required when event_time_omissions is non-empty.",
      call. = FALSE
    )
  }
  if (length(omissions) > 0L && nrow(out) > 0L &&
      total_omitted != length(omissions)) {
    stop(
      "event_time_omission_summary counts must match event_time_omissions length.",
      call. = FALSE
    )
  }

  list(
    event_time_omissions = omissions,
    event_time_omission_summary = out,
    n_omitted_event_times = as.integer(total_omitted)
  )
}

# -- Formatting helpers for Stata-quality output ----------------------------

#' Format a numeric value with fixed width
#' @keywords internal
.fmt_num <- function(x, width = 10L, digits = 4L) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(formatC("", width = width))
  formatC(x, format = "f", digits = digits, width = width)
}

#' Format a p-value with consistent display
#' @keywords internal
.fmt_pval <- function(p, width = 8L) {
  if (is.null(p) || length(p) == 0L || is.na(p)) return(formatC("", width = width))
  if (p < 0.001) return(formatC("<0.001", width = width))
  formatC(p, format = "f", digits = 4L, width = width)
}

#' Format significance stars (fixed width 4)
#' @keywords internal
.fmt_stars <- function(p) {
  if (is.null(p) || is.na(p)) return("    ")
  if (p < 0.001) return(" ***")
  if (p < 0.01)  return("  **")
  if (p < 0.05)  return("   *")
  if (p < 0.10)  return("   .")
  "    "
}

#' Print a formatted effects table (period / cohort / event-time)
#'
#' @param df data.frame with at least 'att', 'se', 'pvalue' columns.
#' @param label_col character, name of the leftmost label column (e.g. "Period").
#' @param label_name character, column name in df for labels.
#' @param extra_cols named list of list(name, col, width, fmt) for extra columns.
#' @param digits integer, decimal places.
#' @param max_rows integer, max rows to show before truncation.
#' @keywords internal
.print_effects_table <- function(df, label_col = "Period", label_name = "period",
                                  extra_cols = NULL, digits = 4L,
                                  max_rows = 20L) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(invisible(NULL))

  n_rows <- nrow(df)
  show_df <- if (n_rows > max_rows) df[seq_len(max_rows), , drop = FALSE] else df

  # Build header
  hdr <- sprintf("  %8s %10s %10s %10s %8s     %10s %10s",
                 label_col, "ATT", "Std. Err.", "t-stat", "p-value",
                 "CI lower", "CI upper")
  rule <- paste0("  ", strrep("-", nchar(hdr) - 2L))

  cat(rule, "\n")
  cat(hdr, "\n")
  cat(rule, "\n")

  for (i in seq_len(nrow(show_df))) {
    row <- show_df[i, ]
    lbl <- if (label_name %in% names(row)) as.character(row[[label_name]]) else ""
    att_val <- if ("att" %in% names(row)) row$att else NA_real_
    se_val <- if ("se" %in% names(row)) row$se else NA_real_
    t_val <- if ("t_stat" %in% names(row)) row$t_stat else NA_real_
    p_val <- if ("pvalue" %in% names(row)) row$pvalue else NA_real_
    ci_lo <- if ("ci_lower" %in% names(row)) row$ci_lower else NA_real_
    ci_hi <- if ("ci_upper" %in% names(row)) row$ci_upper else NA_real_
    if (identical(label_name, "event_time") &&
        suppressWarnings(as.integer(row[[label_name]])) == -1L &&
        isTRUE(all.equal(as.numeric(att_val), 0, tolerance = 1e-12)) &&
        isTRUE(all.equal(as.numeric(se_val), 0, tolerance = 1e-12))) {
      lbl <- paste0(lbl, " (anchor)")
    }

    cat(sprintf("  %8s %s %s %s %s%s %s %s\n",
                lbl,
                .fmt_num(att_val, 10L, digits),
                .fmt_num(se_val, 10L, digits),
                .fmt_num(t_val, 10L, digits),
                .fmt_pval(p_val),
                .fmt_stars(p_val),
                .fmt_num(ci_lo, 10L, digits),
                .fmt_num(ci_hi, 10L, digits)))
  }

  if (n_rows > max_rows) {
    hidden <- n_rows - max_rows
    more_label <- if (identical(label_col, "Period")) {
      "periods"
    } else if (identical(label_col, "Cohort")) {
      "cohorts"
    } else {
      "rows"
    }
    cat(sprintf("  ... (%d more %s)\n", hidden, more_label))
  }
  cat(rule, "\n")
}

#' Print two-column key-value header
#' @keywords internal
.print_header_kv <- function(pairs) {
  # pairs: list of c(key1, val1, key2, val2) vectors
  for (p in pairs) {
    if (length(p) == 4L) {
      cat(sprintf("  %-13s %-24s %-14s %s\n", p[1], p[2], p[3], p[4]))
    } else if (length(p) == 2L) {
      cat(sprintf("  %-13s %s\n", p[1], p[2]))
    }
  }
}

.RULE_DOUBLE <- strrep("=", 78L)
.RULE_SINGLE <- strrep("-", 78L)

# -- print.lwdid_result -----------------------------------------------------

#' @title Print lwdid result
#' @description Compact output with Stata-quality formatting. Delegates to
#'   \code{summary()} for consistent display.
#' @param x lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param ... ignored
#' @return invisible(x)
#' @export
print.lwdid_result <- function(x, digits = 4L, ...) {
  stopifnot(inherits(x, "lwdid_result"))
  s <- summary(x, digits = digits)
  print(s, digits = digits, signif.stars = TRUE, compact = TRUE)
  invisible(x)
}

# -- summary.lwdid_result ---------------------------------------------------

#' @title lwdid result summary
#' @description Generate a full regression summary including period effects,
#'   cohort effects, skipped staggered `(g,r)` diagnostics, omitted event-time
#'   WATT diagnostics, event-time overlap/support extremes when available,
#'   event-time bootstrap replicate diagnostics when available,
#'   parallel trends test, and diagnostics.
#' @param object lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param ... ignored
#' @return summary.lwdid_result object (invisibly)
#' @export
summary.lwdid_result <- function(object, digits = 4L, ...) {
  stopifnot(inherits(object, "lwdid_result"))

  out <- list()
  out$call       <- if (!is.null(object$call)) {
    object$call
  } else if (!is.null(object$metadata) && !is.null(object$metadata$call)) {
    object$metadata$call
  } else {
    NULL
  }
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
  out$n_units <- object$n_units
  out$n_periods <- object$n_periods
  out$n_cohorts <- object$n_cohorts
  out$n_clusters <- object$n_clusters
  out$cohort_sample_sizes <- if (!is.null(object$cohort_sample_sizes)) {
    object$cohort_sample_sizes
  } else {
    object$cohort_sizes
  }
  out$warnings_count <- if (is.null(object$warnings_log)) {
    0L
  } else if (is.data.frame(object$warnings_log)) {
    nrow(object$warnings_log)
  } else {
    length(object$warnings_log)
  }

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
    out$period_effects <- .lwdid_effect_rows_with_supplied_vcov(
      object$att_by_period,
      object,
      "by_period"
    )
    out$n_period_effects <- nrow(object$att_by_period)
  }

  # Staggered-specific
  if (isTRUE(object$is_staggered)) {
    out$is_staggered <- TRUE
    out$cohort_effects <- if (!is.null(object$att_by_cohort)) {
      .lwdid_effect_rows_with_supplied_vcov(
        object$att_by_cohort,
        object,
        "by_cohort"
      )
    } else {
      NULL
    }
    out$event_time_effects <- if (!is.null(object$event_time_effects)) {
      if (is.data.frame(object$event_time_effects)) {
        .lwdid_effect_rows_with_supplied_vcov(
          object$event_time_effects,
          object,
          "event_time"
        )
      } else if (length(object$event_time_effects) == 0L) {
        data.frame()
      } else {
        extract_effects(object, type = "event_time")
      }
    } else {
      NULL
    }
    out$cohort_weights     <- object$cohort_weights
    out$effective_weights  <- object$effective_weights
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
    skipped_diagnostics <- .validate_skipped_diagnostics(
      object$skipped_summary,
      object$skipped_pairs
    )
    event_time_omissions <- .validate_event_time_omissions(
      object$event_time_omission_summary,
      object$event_time_omissions
    )
    out$skipped_pairs <- skipped_diagnostics$skipped_pairs
    out$skipped_summary <- skipped_diagnostics$skipped_summary
    out$n_skipped_pairs <- skipped_diagnostics$n_skipped_pairs
    out$event_time_omissions <- event_time_omissions$event_time_omissions
    out$event_time_omission_summary <-
      event_time_omissions$event_time_omission_summary
    out$n_omitted_event_times <- event_time_omissions$n_omitted_event_times
    out$event_time_support_summary <- .summarize_event_time_support_metadata(
      object
    )
    out$event_time_bootstrap_summary <-
      .summarize_event_time_bootstrap_metadata(object)
    out$has_overall_effect <- identical(object$aggregate, "overall") &&
      !is.null(object$att) && !is.na(object$att)
    out$n_never_treated    <- object$n_never_treated
    out$aggregate          <- object$aggregate %||% "none"
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
      ce_detail <- do.call(rbind, ce_rows)
      if (!is.null(out$cohort_effects) &&
          is.data.frame(out$cohort_effects) &&
          "cohort" %in% names(out$cohort_effects)) {
        cohort_match <- match(
          as.character(ce_detail$cohort),
          as.character(out$cohort_effects$cohort)
        )
        matched <- !is.na(cohort_match)
        aligned_cols <- intersect(
          c("att", "se", "t_stat", "pvalue", "ci_lower", "ci_upper"),
          names(out$cohort_effects)
        )
        for (nm in aligned_cols) {
          if (!nm %in% names(ce_detail)) {
            ce_detail[[nm]] <- NA_real_
          }
          ce_detail[[nm]][matched] <- out$cohort_effects[[nm]][cohort_match[matched]]
        }
      }
      out$cohort_effects_detail <- ce_detail
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
    out$ri_observed_stat  <- object$ri_observed_stat
    out$ri_estimator      <- object$ri_estimator
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
#' @description Full regression summary with Stata-quality formatting,
#'   including concise skipped staggered `(g,r)` and omitted event-time WATT
#'   diagnostics when present.
#' @param x summary.lwdid_result object
#' @param digits integer, decimal places (default 4)
#' @param signif.stars logical, show significance stars (default TRUE)
#' @param compact logical, print shortened diagnostic and effect tables
#'   for inline display.
#' @param ... ignored
#' @return invisible(x)
#' @export
print.summary.lwdid_result <- function(x, digits = 4L,
                                        signif.stars = TRUE,
                                        compact = FALSE, ...) {
  W <- 78L
  RULE2 <- strrep("=", W)
  RULE1 <- strrep("-", W)

  # ── Title ──────────────────────────────────────────────────────────────────
  cat("\nLee-Wooldridge DiD Estimation\n")
  cat(RULE2, "\n")

  # ── Header: two-column key-value ───────────────────────────────────────────
  est_label <- x$estimator %||% "ra"
  vce_label <- x$vce_description %||% "homoskedastic"
  df_label <- if (!is.null(x$wcb_details)) "WCB" else {
    if (!is.null(x$df_inference) && !is.na(x$df_inference))
      as.character(x$df_inference)
    else "NA"
  }

  .print_header_kv(list(
    c("Method:", x$method %||% "common_timing",
      "Dep. var:", x$depvar %||% ""),
    c("Estimator:", est_label,
      "VCE:", vce_label),
    c("Transform:", x$rolling %||% "demean",
      "df:", df_label)
  ))

  # Common-timing specific header fields
  if (!isTRUE(x$is_staggered)) {
    k_str <- if (!is.null(x$K)) sprintf("K = %d", x$K) else ""
    tp_str <- if (!is.null(x$tpost1)) sprintf("tpost1 = %s", x$tpost1) else ""
    if (nzchar(k_str) || nzchar(tp_str)) {
      .print_header_kv(list(c("Pre-periods:", k_str, "Post start:", tp_str)))
    }
  }

  # Staggered specific header fields
  if (isTRUE(x$is_staggered)) {
    agg_label <- x$aggregate %||% "none"
    .print_header_kv(list(
      c("Design:", "Staggered DID", "Aggregate:", agg_label),
      c("Aggregation:", agg_label, "", "")
    ))
    sample_parts <- c()
    if (!is.null(x$n_units) && !is.na(x$n_units)) {
      sample_parts <- c(sample_parts, sprintf("Units: %d", x$n_units))
    }
    if (!is.null(x$n_periods) && !is.na(x$n_periods)) {
      sample_parts <- c(sample_parts, sprintf("Periods: %d", x$n_periods))
    }
    if (!is.null(x$n_cohorts) && !is.na(x$n_cohorts)) {
      sample_parts <- c(sample_parts, sprintf("Cohorts: %d", x$n_cohorts))
    }
    if (!is.null(x$n_never_treated) && !is.na(x$n_never_treated)) {
      sample_parts <- c(sample_parts,
                        sprintf("NT: %d", x$n_never_treated),
                        sprintf("N Never-treated: %d", x$n_never_treated))
    }
    if (length(sample_parts) > 0L) {
      cat(sprintf("  %-13s %s\n", "Sample:", paste(sample_parts, collapse = " | ")))
    }
    if (!is.null(x$cohort_list) && length(x$cohort_list) > 0L) {
      cat(sprintf("  %-13s %s\n", "Cohorts:",
                  paste(x$cohort_list, collapse = ", ")))
    }
    if (!is.null(x$control_group_used)) {
      cg_str <- x$control_group_used
      if (isTRUE(x$control_group_auto_switched) && !is.null(x$control_group)) {
        cg_str <- sprintf("%s [auto-switched from %s]", cg_str, x$control_group)
      }
      nt_str <- if (!is.null(x$n_never_treated))
        sprintf("N Never-treated: %d", x$n_never_treated) else ""
      .print_header_kv(list(c("Control:", cg_str, nt_str, "")))
    }
    if (!is.null(x$skipped_summary) &&
        is.data.frame(x$skipped_summary) &&
        nrow(x$skipped_summary) > 0L) {
      reason_summary <- paste(
        sprintf("%s: %d", x$skipped_summary$reason,
                as.integer(x$skipped_summary$n)),
        collapse = ", "
      )
      cat(sprintf(
        "  %-15s %d (%s)\n",
        "Skipped (g,r):",
        as.integer(x$n_skipped_pairs %||% sum(x$skipped_summary$n)),
        reason_summary
      ))
    }
    if (!is.null(x$event_time_omission_summary) &&
        is.data.frame(x$event_time_omission_summary) &&
        nrow(x$event_time_omission_summary) > 0L) {
      omission_summary <- paste(
        sprintf("%s: %d", x$event_time_omission_summary$reason,
                as.integer(x$event_time_omission_summary$n)),
        collapse = ", "
      )
      cat(sprintf(
        "  %-15s %d (%s)\n",
        "Omitted WATT(e):",
        as.integer(
          x$n_omitted_event_times %||% sum(x$event_time_omission_summary$n)
        ),
        omission_summary
      ))
    }
  }

  if (!is.null(x$n_clusters) && !is.na(x$n_clusters)) {
    .print_header_kv(list(c("Clusters:", sprintf("%d", as.integer(x$n_clusters)))))
  }

  cat(RULE1, "\n")

  # ── Overall ATT ────────────────────────────────────────────────────────────
  att_val <- x$coefficients$Estimate
  se_val <- x$coefficients$Std.Error
  t_val <- x$coefficients$t.value
  p_val <- x$coefficients$Pr.t
  ci_lo <- x$coefficients$CI.lower
  ci_hi <- x$coefficients$CI.upper
  has_se <- !is.null(se_val) && !is.na(se_val)

  if (isTRUE(x$is_staggered) && isTRUE(x$has_overall_effect)) {
    cat("\nOverall ATT (tau_omega):\n")
  } else if (isTRUE(x$is_staggered) &&
             identical(x$aggregate, "event_time")) {
    cat("\nCohort-period ATT summary:\n")
  } else {
    cat("\nOverall ATT:\n")
  }

  has_event_time_rows <- isTRUE(x$is_staggered) &&
    identical(x$aggregate, "event_time") &&
    !is.null(x$event_time_effects) &&
    (
      (is.data.frame(x$event_time_effects) && nrow(x$event_time_effects) > 0L) ||
        (!is.data.frame(x$event_time_effects) && length(x$event_time_effects) > 0L)
    )

  # Header line for ATT table
  hdr <- sprintf("              %10s %10s %10s %8s     %10s %10s",
                 "Estimate", "Std. Err.", "t value", "Pr(>|t|)",
                 "CI lower", "CI upper")
  cat(hdr, "\n")

  # ATT row
  cat(sprintf("  ATT     %s %s %s %s%s %s %s\n",
              .fmt_num(att_val, 10L, digits),
              .fmt_num(se_val, 10L, digits),
              .fmt_num(t_val, 10L, digits),
              .fmt_pval(p_val),
              if (has_se) .fmt_stars(p_val) else "    ",
              .fmt_num(ci_lo, 10L, digits),
              .fmt_num(ci_hi, 10L, digits)))

  if (signif.stars && has_se) {
    cat("  ---\n")
    cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  if (isTRUE(x$is_staggered) && identical(x$aggregate, "event_time")) {
    if (has_event_time_rows) {
      cat("  Note: this scalar summarizes valid (g,r) effects; event-time WATT(e) rows are below.\n")
      cat("  Note: WATT(e) standard errors use diagonal cohort-SE aggregation; cross-cohort covariance is not modeled.\n")
    } else {
      cat("  Note: no finite event-time WATT(e) rows are available for display.\n")
    }
  }

  # ── Regression Coefficients (common timing) ───────────────────────────────
  if (!is.null(x$regression_coefficients) &&
      is.data.frame(x$regression_coefficients) &&
      nrow(x$regression_coefficients) > 0L) {
    cat("\nRegression Coefficients:\n")
    rc <- x$regression_coefficients
    hdr_rc <- sprintf("  %-16s %10s %10s %10s %8s",
                      "Term", "Estimate", "Std. Err.", "t value", "Pr(>|t|)")
    rule_rc <- paste0("  ", strrep("-", nchar(hdr_rc) - 2L))
    cat(rule_rc, "\n")
    cat(hdr_rc, "\n")
    cat(rule_rc, "\n")
    term_names <- rownames(rc)
    for (i in seq_len(nrow(rc))) {
      cat(sprintf("  %-16s %s %s %s %s%s\n",
                  term_names[[i]],
                  .fmt_num(rc$Estimate[[i]], 10L, digits),
                  .fmt_num(rc$Std.Error[[i]], 10L, digits),
                  .fmt_num(rc$t.value[[i]], 10L, digits),
                  .fmt_pval(rc$Pr.t[[i]]),
                  .fmt_stars(rc$Pr.t[[i]])))
    }
    cat(rule_rc, "\n")
  }

  # ── WCB details ────────────────────────────────────────────────────────────
  if (!is.null(x$wcb_details)) {
    cat("\n  Wild Cluster Bootstrap:\n")
    cat(sprintf("    Weight type: %s\n", x$wcb_details$weight_type))
    cat(sprintf("    Clusters: %d\n", as.integer(x$wcb_details$n_clusters)))
    cat(sprintf("    Actual bootstrap draws: %d\n",
                as.integer(x$wcb_details$actual_n_bootstrap)))
    cat(sprintf("    Requested bootstrap draws: %d\n",
                as.integer(x$wcb_details$requested_n_bootstrap)))
    if (!is.null(x$wcb_details$restricted_model)) {
      cat(sprintf("    Restricted model: %s\n", x$wcb_details$restricted_model))
    }
  }

  # ── Cohort Effects (staggered) ─────────────────────────────────────────────
  if (isTRUE(x$is_staggered) && !is.null(x$cohort_effects_detail) &&
      is.data.frame(x$cohort_effects_detail) && nrow(x$cohort_effects_detail) > 0L) {
    cat("\nCohort Effects:\n")
    # Build cohort effects with CI from the full cohort_effects list
    ce <- x$cohort_effects_detail
    ce_full <- NULL
    if (!is.null(x$cohort_effects) && is.data.frame(x$cohort_effects)) {
      ce_full <- x$cohort_effects
    }
    # Construct display df with CI
    ce_df <- data.frame(
      cohort = ce$cohort,
      att = ce$att,
      se = ce$se,
      stringsAsFactors = FALSE
    )
    if ("pvalue" %in% names(ce)) {
      ce_df$pvalue <- ce$pvalue
    }
    if ("ci_lower" %in% names(ce)) {
      ce_df$ci_lower <- ce$ci_lower
    }
    if ("ci_upper" %in% names(ce)) {
      ce_df$ci_upper <- ce$ci_upper
    }
    # Add pvalue, ci_lower, ci_upper from full cohort_effects if available
    if (!"pvalue" %in% names(ce_df) &&
        !is.null(ce_full) && "pvalue" %in% names(ce_full)) {
      ce_df$pvalue <- ce_full$pvalue
    } else if (!"pvalue" %in% names(ce_df)) {
      ce_df$pvalue <- NA_real_
    }
    if (!"ci_lower" %in% names(ce_df) &&
        !is.null(ce_full) && "ci_lower" %in% names(ce_full)) {
      ce_df$ci_lower <- ce_full$ci_lower
    } else if (!"ci_lower" %in% names(ce_df)) {
      ce_df$ci_lower <- NA_real_
    }
    if (!"ci_upper" %in% names(ce_df) &&
        !is.null(ce_full) && "ci_upper" %in% names(ce_full)) {
      ce_df$ci_upper <- ce_full$ci_upper
    } else if (!"ci_upper" %in% names(ce_df)) {
      ce_df$ci_upper <- NA_real_
    }
    ce_df$t_stat <- if ("t_stat" %in% names(ce)) {
      ce$t_stat
    } else {
      ifelse(ce_df$se > 0, ce_df$att / ce_df$se, NA_real_)
    }

    # Print header
    hdr_ce <- sprintf("  %8s %7s %7s %10s %10s %10s %8s     %10s %10s",
                      "Cohort", "N_units", "Periods",
                      "ATT", "Std. Err.", "t-stat", "p-value",
                      "CI lower", "CI upper")
    rule_ce <- paste0("  ", strrep("-", nchar(hdr_ce) - 2L))
    cat(rule_ce, "\n")
    cat(hdr_ce, "\n")
    cat(rule_ce, "\n")

    for (i in seq_len(nrow(ce_df))) {
      r <- ce_df[i, ]
      nu <- if (!is.null(ce$n_units)) as.integer(ce$n_units[i]) else NA_integer_
      np <- if (!is.null(ce$n_periods)) as.integer(ce$n_periods[i]) else NA_integer_
      cat(sprintf("  %8d %7s %7s %s %s %s %s%s %s %s\n",
                  r$cohort,
                  if (!is.na(nu)) formatC(nu, width = 7L) else formatC("", width = 7L),
                  if (!is.na(np)) formatC(np, width = 7L) else formatC("", width = 7L),
                  .fmt_num(r$att, 10L, digits),
                  .fmt_num(r$se, 10L, digits),
                  .fmt_num(r$t_stat, 10L, digits),
                  .fmt_pval(r$pvalue),
                  .fmt_stars(r$pvalue),
                  .fmt_num(r$ci_lower, 10L, digits),
                  .fmt_num(r$ci_upper, 10L, digits)))
    }
    cat(rule_ce, "\n")
  }

  if (isTRUE(x$is_staggered) && !is.null(x$cohort_weights) &&
      length(x$cohort_weights) > 0L) {
    cat("\nCohort Weights:\n")
    has_effective_weight_deviation <- FALSE
    for (cohort_id in names(x$cohort_weights)) {
      weight_val <- x$cohort_weights[[cohort_id]]
      effective_val <- NULL
      if (!is.null(x$effective_weights) &&
          cohort_id %in% names(x$effective_weights)) {
        effective_val <- x$effective_weights[[cohort_id]]
      }
      n_val <- NA_integer_
      if (!is.null(x$cohort_sample_sizes) && cohort_id %in% names(x$cohort_sample_sizes)) {
        n_val <- as.integer(x$cohort_sample_sizes[[cohort_id]])
      }
      line <- sprintf("  %s  weight=%s",
                      cohort_id,
                      .fmt_num(weight_val, 8L, digits))
      if (!is.na(n_val)) {
        line <- sprintf("%s  N=%d", line, n_val)
      }
      if (!is.null(effective_val) &&
          is.finite(as.numeric(weight_val)) &&
          is.finite(as.numeric(effective_val)) &&
          abs(as.numeric(weight_val) - as.numeric(effective_val)) > 1e-12) {
        line <- sprintf(
          "%s  effective=%s",
          line,
          .fmt_num(effective_val, 8L, digits)
        )
        has_effective_weight_deviation <- TRUE
      }
      cat(line, "\n")
    }
    if (has_effective_weight_deviation) {
      cat("  Note: effective weights reflect the post-dropna regression sample.\n")
    }
  }

  # ── Event-Time Effects (staggered) ─────────────────────────────────────────
  if (!is.null(x$event_time_effects) && length(x$event_time_effects) > 0L) {
    cat("\nEvent-Time Effects:\n")
    # Convert to data.frame
    if (is.data.frame(x$event_time_effects)) {
      et_df <- x$event_time_effects
    } else {
      et_df <- do.call(rbind, lapply(x$event_time_effects, function(e) {
        data.frame(
          event_time = as.integer(e$event_time),
          att = as.numeric(e$att),
          se = as.numeric(e$se),
          t_stat = as.numeric(e$t_stat %||% (e$att / e$se)),
          pvalue = as.numeric(e$pvalue),
          ci_lower = as.numeric(e$ci_lower),
          ci_upper = as.numeric(e$ci_upper),
          stringsAsFactors = FALSE
        )
      }))
    }
    .print_effects_table(et_df, label_col = "Rel.time",
                         label_name = "event_time", digits = digits)
    if (isTRUE(x$is_staggered) && identical(x$aggregate, "event_time")) {
      cat("  WATT(e) standard errors use diagonal cohort-SE aggregation; cross-cohort covariance is not modeled.\n")
      has_weight_cv <- "max_weight_cv" %in% names(et_df)
      has_treated_support <- "min_n_treated" %in% names(et_df)
      max_weight_cv <- if (has_weight_cv) {
        suppressWarnings(max(et_df$max_weight_cv, na.rm = TRUE))
      } else {
        NA_real_
      }
      min_n_treated <- if (has_treated_support) {
        suppressWarnings(min(et_df$min_n_treated, na.rm = TRUE))
      } else {
        NA_integer_
      }
      show_overlap_cue <- is.finite(max_weight_cv) && max_weight_cv > 2
      show_support_cue <- is.finite(min_n_treated) && min_n_treated <= 5
      if (show_overlap_cue || show_support_cue) {
        cue_parts <- character(0)
        if (show_overlap_cue) {
          cue_parts <- c(
            cue_parts,
            sprintf("max weight CV=%.2f", max_weight_cv)
          )
        }
        if (show_support_cue) {
          cue_parts <- c(
            cue_parts,
            sprintf("min treated cell N=%d", as.integer(min_n_treated))
          )
        }
        cat(sprintf(
          "  Diagnostics: event-time rows retain overlap/support metadata (%s).\n",
          paste(cue_parts, collapse = ", ")
        ))
      }
      if (!is.null(x$event_time_bootstrap_summary) &&
          is.character(x$event_time_bootstrap_summary$cue) &&
          !is.na(x$event_time_bootstrap_summary$cue)) {
        cat(sprintf(
          "  Diagnostics: event-time bootstrap replicate metadata retained (%s).\n",
          x$event_time_bootstrap_summary$cue
        ))
      }
    }
  }

  # ── Period-Specific Effects (common timing) ────────────────────────────────
  if (!is.null(x$period_effects) && is.data.frame(x$period_effects) &&
      nrow(x$period_effects) > 0L) {
    cat("\nPeriod-Specific Effects:\n")
    .print_effects_table(x$period_effects, label_col = "Period",
                         label_name = "period", digits = digits,
                         max_rows = if (isTRUE(compact)) 5L else Inf)
    cat("  Period-by-period effects are available in results$att_by_period.\n")
  }

  # ── Pre-treatment Dynamics ─────────────────────────────────────────────────
  if (!is.null(x$att_pre_treatment) && isTRUE(x$include_pretreatment)) {
    cat("\nPre-treatment Dynamics:\n")
    .print_effects_table(x$att_pre_treatment, label_col = "e(t)",
                         label_name = "event_time", digits = digits)
  }

  # ── Parallel Trends Test ───────────────────────────────────────────────────
  if (!is.null(x$parallel_trends)) {
    pt <- x$parallel_trends
    cat("\nParallel Trends Test:\n")
    f_stat_val <- pt$f_stat %||% NA_real_
    f_pval_val <- pt$f_pvalue %||% NA_real_
    df1_val <- pt$df1 %||% 0L
    df2_val <- pt$df2 %||% 0L
    if (!is.na(f_stat_val) && !is.na(f_pval_val)) {
      cat(sprintf("  Joint F-test: F(%d,%d) = %.4f, p = %s\n",
                  df1_val, df2_val, f_stat_val, .fmt_pval(f_pval_val, 1L)))
    }
    alpha_val <- x$alpha %||% 0.05
    if (isTRUE(pt$reject)) {
      cat(sprintf("  => Reject H0 at alpha=%.2f: evidence against parallel trends\n",
                  alpha_val))
    } else {
      cat(sprintf("  => Fail to reject H0 at alpha=%.2f: no evidence against PT\n",
                  alpha_val))
    }
  }

  # ── Randomization Inference ────────────────────────────────────────────────
  if (!is.null(x$ri_pvalue) || !is.null(x$ri_error)) {
    cat("\nRandomization Inference (RI):\n")
    if (!is.null(x$ri_pvalue)) {
      cat(sprintf("  RI p-value = %s%s\n",
                  .fmt_pval(x$ri_pvalue, 1L), .fmt_stars(x$ri_pvalue)))
    }
    details <- c()
    ri_total <- .ri_total_permutations(
      valid = x$ri_n_valid, failed = x$ri_n_failed,
      reps = x$ri_n_permutations
    )
    if (!is.null(x$ri_method))  details <- c(details, sprintf("method=%s", x$ri_method))
    if (!is.null(x$ri_estimator)) details <- c(details, sprintf("estimator=%s", x$ri_estimator))
    if (!is.null(x$ri_target))  details <- c(details, sprintf("target=%s", x$ri_target))
    if (!is.null(x$rireps))     details <- c(details, sprintf("reps=%d", x$rireps))
    if (!is.null(x$ri_seed))    details <- c(details, sprintf("seed=%d", x$ri_seed))
    if (!is.null(x$ri_n_valid)) {
      details <- c(details,
                   sprintf("valid=%s", .format_ri_valid_count(x$ri_n_valid, ri_total)))
    }
    if (!is.null(x$ri_n_failed)) details <- c(details, sprintf("failed=%d", x$ri_n_failed))
    if (length(details) > 0L) cat(sprintf("  %s\n", paste(details, collapse = " | ")))
    if (!is.null(x$ri_observed_stat) && length(x$ri_observed_stat) > 0L &&
        is.finite(x$ri_observed_stat[1L])) {
      cat(sprintf("  Observed statistic = %.*f\n",
                  digits, x$ri_observed_stat[1L]))
    }
    if (!is.null(x$ri_dist_summary)) {
      ds <- x$ri_dist_summary
      cat(sprintf("  Distribution: mean=%.*f, median=%.*f, sd=%.*f\n",
                  digits, ds$mean, digits, ds$median, digits, ds$sd))
    }
    if (!is.null(x$ri_error)) cat(sprintf("  RI Error: %s\n", x$ri_error))
  }

  # ── Diagnostics ────────────────────────────────────────────────────────────
  if (!is.null(x$diagnostics_summary) && length(x$diagnostics_summary) > 0L) {
    cat("\nDiagnostics:\n")
    for (nm in names(x$diagnostics_summary)) {
      cat(sprintf("  %s: %s\n", nm, x$diagnostics_summary[[nm]]))
    }
  }

  if (!is.null(x$exclude_pre_periods) && x$exclude_pre_periods > 0L) {
    cat(sprintf("  Excluded pre-periods: %d\n", x$exclude_pre_periods))
  }

  if (!is.null(x$warnings_count) && is.finite(x$warnings_count) &&
      x$warnings_count > 0L) {
    cat(sprintf("  Warnings: %d\n", as.integer(x$warnings_count)))
  }

  # ── Footer ─────────────────────────────────────────────────────────────────
  cat(RULE2, "\n")
  parts <- c(sprintf("N = %d", x$nobs),
             sprintf("N_treated = %d", x$n_treated),
             sprintf("N_control = %d", x$n_control))
  if (!isTRUE(x$is_staggered)) {
    if (!is.null(x$K)) parts <- c(parts, sprintf("K = %d", x$K))
    if (!is.null(x$tpost1)) parts <- c(parts, sprintf("tpost1 = %s", x$tpost1))
  }
  cat(sprintf("  %s\n", paste(parts, collapse = "    ")))

  if (isTRUE(x$is_staggered)) {
    if (identical(x$aggregate, "none")) {
      cat("  Hint: inspect results$att_by_cohort_time for (g,r)-level effects; use plot_event_study(x) for event-study plots.\n")
    } else if (identical(x$aggregate, "overall")) {
      cat("  Hint: use extract_effects(x, type='overall') for the Overall Effect and extract_effects(x, type='cohort') for cohort details.\n")
    } else {
      cat("  Hint: coef(x, type='all') for (g,r)-level effects, plot(x) for event study\n")
    }
  }

  invisible(x)
}

# -- coef.lwdid_result ------------------------------------------------------

#' @title Extract ATT coefficients
#' @description Extract ATT coefficients from lwdid_result object.
#' @param object lwdid_result object
#' @param type character, granularity: "overall" (default), "by_period",
#'   "by_cohort", "event_time", "all"
#' @param ... ignored
#' @return named numeric vector
#' @export
coef.lwdid_result <- function(object, type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "event_time", "all"))

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
    "event_time" = {
      et <- extract_effects(object, type = "event_time")
      if (nrow(et) == 0L)
        stop("No event-time results available.", call. = FALSE)
      out <- et$att
      names(out) <- paste0("e", et$event_time)
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
#'   "by_cohort", "event_time", "all"
#' @param ... ignored
#' @return matrix with lower and upper bounds
#' @export
confint.lwdid_result <- function(object, parm = NULL, level = 0.95,
                                  type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  if (level <= 0 || level >= 1)
    stop("level must be in (0, 1)", call. = FALSE)
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "event_time", "all"))

  info <- switch(type,
    "overall" = {
      list(est = object$att, se = object$se_att, rn = "ATT",
           df = object$df_inference)
    },
    "by_period" = {
      if (is.null(object$att_by_period))
        stop("No period-specific results for confint.", call. = FALSE)
      period_vcov <- vcov(object, type = "by_period")
      list(est = object$att_by_period$att,
           se  = sqrt(diag(period_vcov)),
           rn  = as.character(object$att_by_period$period),
           df = object$df_inference)
    },
    "by_cohort" = {
      if (is.null(object$att_by_cohort))
        stop("No cohort-specific results (Staggered mode only).",
             call. = FALSE)
      cohort_vcov <- vcov(object, type = "by_cohort")
      list(est = object$att_by_cohort$att,
           se  = sqrt(diag(cohort_vcov)),
           rn  = as.character(object$att_by_cohort$cohort),
           df = object$df_inference)
    },
    "event_time" = {
      et <- extract_effects(object, type = "event_time")
      if (nrow(et) == 0L)
        stop("No event-time results for confint.", call. = FALSE)
      event_vcov <- vcov(object, type = "event_time")
      list(est = et$att,
           se = sqrt(diag(event_vcov)),
           rn = paste0("e", et$event_time),
           df = et$df_inference)
    },
    "all" = {
      if (!is.null(object$att_by_cohort_time)) {
        list(est = object$att_by_cohort_time$att,
             se  = object$att_by_cohort_time$se,
             rn  = sprintf("g%s.r%s",
                           object$att_by_cohort_time$cohort,
                           object$att_by_cohort_time$period),
             df = object$df_inference)
      } else if (!is.null(object$att_by_period)) {
        period_vcov <- vcov(object, type = "by_period")
        list(est = object$att_by_period$att,
             se  = sqrt(diag(period_vcov)),
             rn  = as.character(object$att_by_period$period),
             df = object$df_inference)
      } else {
        list(est = object$att, se = object$se_att, rn = "ATT",
             df = object$df_inference)
      }
    }
  )

  alpha <- 1 - level
  df <- info$df
  if (length(df) == 1L) {
    df <- rep(df, length(info$est))
  }
  inference_dist <- object$inference_dist %||%
    .resolve_top_level_inference_dist(object$estimator)
  if (identical(inference_dist, "normal")) {
    t_crit <- stats::qnorm(1 - alpha / 2)
  } else {
    t_crit <- stats::qt(1 - alpha / 2, df = df)
  }

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

.validate_lwdid_effect_vcov <- function(V, expected_n, label, effects_label) {
  if (!is.matrix(V) || !is.numeric(V)) {
    stop(sprintf("%s must be a numeric matrix.", label), call. = FALSE)
  }

  expected_dim <- c(expected_n, expected_n)
  if (!identical(dim(V), expected_dim)) {
    stop(sprintf(
      "%s must be a %dx%d matrix matching %s.",
      label, expected_dim[1L], expected_dim[2L], effects_label
    ), call. = FALSE)
  }

  v_diag <- diag(V)
  if (any(!is.finite(v_diag) | v_diag < 0)) {
    stop(
      sprintf("%s diagonal variances must be finite and non-negative.", label),
      call. = FALSE
    )
  }

  if (any(!is.finite(V))) {
    stop(sprintf("%s covariance entries must be finite.", label), call. = FALSE)
  }

  if (!isTRUE(all.equal(V, t(V), tolerance = 1e-10, check.attributes = FALSE))) {
    stop(sprintf("%s must be symmetric.", label), call. = FALSE)
  }

  V
}

# -- vcov.lwdid_result ------------------------------------------------------

#' @title Variance-covariance matrix
#' @description Extract variance-covariance matrix for ATT estimates.
#' @param object lwdid_result object
#' @param type character, granularity: "overall" (1x1), "by_period",
#'   "by_cohort", "event_time", "full" (complete OLS parameter VCE)
#' @param ... ignored
#' @return variance-covariance matrix. For `type = "event_time"`, this is the
#'   diagonal matrix formed from marginal WATT(e) standard errors unless a
#'   fitted object carries a numeric, symmetric `vcov_att_event_time` matrix
#'   with dimensions matching the event-time rows. Attributes record the
#'   standard-error aggregation rule and covariance boundary; supplied joint
#'   matrices keep their own `covariance_assumption` attribute when present.
#' @export
vcov.lwdid_result <- function(object, type = "overall", ...) {
  stopifnot(inherits(object, "lwdid_result"))
  type <- match.arg(type, c("overall", "by_period", "by_cohort", "event_time", "full"))

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
        V <- .validate_lwdid_effect_vcov(
          object$vcov_att_periods,
          nrow(object$att_by_period),
          "vcov_att_periods",
          "period-specific effects"
        )
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
        V <- .validate_lwdid_effect_vcov(
          object$vcov_att_cohorts,
          nrow(object$att_by_cohort),
          "vcov_att_cohorts",
          "cohort-specific effects"
        )
      } else {
        ses <- object$att_by_cohort$se
        V <- diag(ses^2, nrow = length(ses))
      }
      cnames <- as.character(object$att_by_cohort$cohort)
      rownames(V) <- colnames(V) <- cnames
      V
    },
    "event_time" = {
      et <- extract_effects(object, type = "event_time")
      if (nrow(et) == 0L)
        stop("No event-time results for vcov.", call. = FALSE)
      enames <- paste0("e", et$event_time)
      if (!is.null(object$vcov_att_event_time)) {
        V <- .validate_lwdid_effect_vcov(
          object$vcov_att_event_time,
          nrow(et),
          "vcov_att_event_time",
          "event-time effects"
        )
        covariance_assumption <- attr(V, "covariance_assumption", exact = TRUE)
        if (is.null(covariance_assumption)) {
          covariance_assumption <- "provided_joint_event_time_covariance"
        }
        se_aggregation <- attr(V, "se_aggregation", exact = TRUE)
        if (is.null(se_aggregation)) {
          se_aggregation <- "provided_event_time_vcov_diagonal"
        }
      } else {
        V <- diag(et$se^2, nrow = nrow(et))
        covariance_assumption <- .event_time_metadata_values(
          et$covariance_assumption,
          "zero_cross_cohort_covariance"
        )
        se_aggregation <- .event_time_metadata_values(
          et$se_aggregation,
          "diagonal_weighted_cohort_se"
        )
      }
      rownames(V) <- colnames(V) <- enames
      attr(V, "se_aggregation") <- se_aggregation
      attr(V, "covariance_assumption") <- covariance_assumption
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
