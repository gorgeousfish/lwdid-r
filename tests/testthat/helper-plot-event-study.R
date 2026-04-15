# ============================================================================
# helper-plot-event-study.R — Mock data factories & fixtures for Event Study
# plot tests.
#
# All helpers build lwdid_result objects directly via
#   structure(list(...), class = "lwdid_result")
# to avoid any dependency on the package constructor during testing.
# ============================================================================

# ── Mock factories ──────────────────────────────────────────────────────────

#' Create a mock lwdid_result for Staggered mode
#'
#' @param cohort_data  data.frame with columns: cohort, period, att, se,
#'   df_inference.
#' @param cohort_sizes Named numeric vector, e.g. c("2005" = 100, "2006" = 300).
#' @param att_pre_treatment Optional data.frame (same schema as cohort_data).
#' @param include_pretreatment Logical.
#' @return An object of class "lwdid_result".
create_mock_staggered_result <- function(cohort_data,
                                         cohort_sizes,
                                         att_pre_treatment = NULL,
                                         include_pretreatment = FALSE) {
  cohort_weights <- cohort_sizes / sum(cohort_sizes)

  obj <- list(
    # Core estimates
    att              = mean(cohort_data$att),
    se_att           = 0.1,
    t_stat           = NA_real_,
    pvalue           = NA_real_,
    ci_lower         = NA_real_,
    ci_upper         = NA_real_,
    df_resid         = NA_integer_,
    df_inference     = min(cohort_data$df_inference, na.rm = TRUE),
    nobs             = as.integer(sum(cohort_sizes)),
    n_treated        = NA_integer_,
    n_control        = NA_integer_,
    K                = NA_integer_,
    tpost1           = NA_integer_,
    depvar           = NA_character_,
    rolling          = NA_character_,
    vce_type         = NULL,
    cluster_var      = NULL,
    n_clusters       = NULL,
    cmd              = "lwdid",
    estimator        = NULL,
    method           = "staggered",
    alpha            = 0.05,
    is_staggered     = TRUE,
    controls_used    = FALSE,
    controls         = character(0),
    include_pretreatment = include_pretreatment,
    control_group    = NULL,
    control_group_used = NULL,
    params           = numeric(0),
    bse              = numeric(0),
    vcov_matrix      = NULL,
    resid            = numeric(0),
    data             = NULL,
    # Period / cohort-time tables
    att_by_period    = NULL,
    att_pre_treatment = att_pre_treatment,
    parallel_trends_test = NULL,
    aggregate        = "none",
    cohorts          = unique(cohort_data$cohort),
    cohort_sizes     = cohort_sizes,
    n_never_treated  = NULL,
    att_by_cohort    = NULL,
    att_by_cohort_time = cohort_data,
    att_overall      = NULL,
    se_overall       = NULL,
    ci_overall_lower = NULL,
    ci_overall_upper = NULL,
    t_stat_overall   = NULL,
    pvalue_overall   = NULL,
    cohort_effects   = NULL,
    event_time_effects = NULL,
    cohort_weights   = cohort_weights,
    # RI fields
    ri_pvalue        = NULL,
    ri_seed          = NULL,
    rireps           = NULL,
    ri_method        = NULL,
    ri_valid         = NULL,
    ri_failed        = NULL,
    ri_error         = NULL,
    ri_target        = NULL,
    ri_distribution  = NULL,
    # Diagnostics
    diagnostics      = NULL,
    warning_diagnostics = list(),
    propensity_scores = NULL,
    matched_data     = NULL,
    n_matched        = NULL,
    match_rate       = NULL,
    weights_cv       = NULL,
    warnings_log     = list(),
    call             = NULL,
    lwdid_version    = "dev",
    ivar             = NA_character_,
    tvar             = NA_character_,
    is_quarterly     = FALSE
  )
  structure(obj, class = "lwdid_result")
}


#' Create a mock lwdid_result for Common Timing mode
#'
#' @param period_data  data.frame with columns: period, att, se (optionally
#'   event_time, df_inference).
#' @param df_inference Global degrees of freedom (used when period_data lacks
#'   its own df_inference column).
#' @param treatment_time Integer, the treatment time point.
#' @param att_pre_treatment Optional data.frame.
#' @param include_pretreatment Logical.
#' @return An object of class "lwdid_result".
create_mock_common_timing_result <- function(period_data,
                                             df_inference = 50L,
                                             treatment_time = 2000L,
                                             att_pre_treatment = NULL,
                                             include_pretreatment = FALSE) {
  obj <- list(
    # Core estimates
    att              = period_data$att[1L],
    se_att           = 0.1,
    t_stat           = NA_real_,
    pvalue           = NA_real_,
    ci_lower         = NA_real_,
    ci_upper         = NA_real_,
    df_resid         = NA_integer_,
    df_inference     = as.integer(df_inference),
    nobs             = 100L,
    n_treated        = NA_integer_,
    n_control        = NA_integer_,
    K                = NA_integer_,
    tpost1           = NA_integer_,
    depvar           = NA_character_,
    rolling          = NA_character_,
    vce_type         = NULL,
    cluster_var      = NULL,
    n_clusters       = NULL,
    cmd              = "lwdid",
    estimator        = NULL,
    method           = "common_timing",
    alpha            = 0.05,
    is_staggered     = FALSE,
    controls_used    = FALSE,
    controls         = character(0),
    include_pretreatment = include_pretreatment,
    control_group    = NULL,
    control_group_used = NULL,
    params           = numeric(0),
    bse              = numeric(0),
    vcov_matrix      = NULL,
    resid            = numeric(0),
    data             = NULL,
    # Period table
    att_by_period    = period_data,
    att_pre_treatment = att_pre_treatment,
    parallel_trends_test = NULL,
    treatment_time   = as.integer(treatment_time),
    # Staggered fields (all NULL for CT)
    aggregate        = NULL,
    cohorts          = NULL,
    cohort_sizes     = NULL,
    n_never_treated  = NULL,
    att_by_cohort    = NULL,
    att_by_cohort_time = NULL,
    att_overall      = NULL,
    se_overall       = NULL,
    ci_overall_lower = NULL,
    ci_overall_upper = NULL,
    t_stat_overall   = NULL,
    pvalue_overall   = NULL,
    cohort_effects   = NULL,
    event_time_effects = NULL,
    cohort_weights   = NULL,
    # RI fields
    ri_pvalue        = NULL,
    ri_seed          = NULL,
    rireps           = NULL,
    ri_method        = NULL,
    ri_valid         = NULL,
    ri_failed        = NULL,
    ri_error         = NULL,
    ri_target        = NULL,
    ri_distribution  = NULL,
    # Diagnostics
    diagnostics      = NULL,
    warning_diagnostics = list(),
    propensity_scores = NULL,
    matched_data     = NULL,
    n_matched        = NULL,
    match_rate       = NULL,
    weights_cv       = NULL,
    warnings_log     = list(),
    call             = NULL,
    lwdid_version    = "dev",
    ivar             = NA_character_,
    tvar             = NA_character_,
    is_quarterly     = FALSE
  )
  structure(obj, class = "lwdid_result")
}


# ── Pre-defined fixtures ────────────────────────────────────────────────────
# Each fixture is a *function* so every test gets a fresh copy.

#' Standard multi-cohort staggered result (2 cohorts: 2005, 2007)
fixture_staggered_result <- function() {
  cohort_data <- data.frame(
    cohort       = as.integer(rep(c(2005L, 2007L), each = 5L)),
    period       = as.integer(c(
      2003L, 2004L, 2005L, 2006L, 2007L,
      2005L, 2006L, 2007L, 2008L, 2009L
    )),
    att          = c(-0.1, -0.05, 0.5, 0.8, 1.0,
                     -0.2, -0.1,  0.6, 0.9, 1.2),
    se           = c( 0.2,  0.15, 0.3, 0.25, 0.2,
                      0.25, 0.2,  0.35, 0.3, 0.25),
    df_inference = rep(c(10L, 20L), each = 5L)
  )
  cohort_sizes <- c("2005" = 100, "2007" = 300)

  create_mock_staggered_result(cohort_data, cohort_sizes)
}


#' Unbalanced cohort set (3 cohorts with different period counts)
fixture_staggered_result_unbalanced <- function() {
  # Cohort 2004: event_times -1, 0, 1  (3 periods)
  c2004 <- data.frame(
    cohort       = 2004L,
    period       = c(2003L, 2004L, 2005L),
    att          = c(-0.08, 0.45, 0.70),
    se           = c( 0.18, 0.28, 0.22),
    df_inference = 12L
  )
  # Cohort 2006: event_times -2, -1, 0, 1, 2  (5 periods)
  c2006 <- data.frame(
    cohort       = 2006L,
    period       = c(2004L, 2005L, 2006L, 2007L, 2008L),
    att          = c(-0.15, -0.07, 0.55, 0.85, 1.10),
    se           = c( 0.22,  0.17, 0.32, 0.27, 0.23),
    df_inference = 18L
  )
  # Cohort 2008: event_times -2, -1, 0, 1, 2, 3  (6 periods)
  c2008 <- data.frame(
    cohort       = 2008L,
    period       = c(2006L, 2007L, 2008L, 2009L, 2010L, 2011L),
    att          = c(-0.12, -0.06, 0.50, 0.78, 1.05, 1.25),
    se           = c( 0.20,  0.16, 0.30, 0.26, 0.21, 0.19),
    df_inference = 25L
  )
  cohort_data  <- rbind(c2004, c2006, c2008)
  cohort_sizes <- c("2004" = 50, "2006" = 200, "2008" = 150)

  create_mock_staggered_result(cohort_data, cohort_sizes)
}


#' Staggered result with pre-treatment data
fixture_staggered_result_with_pre <- function() {
  # Same post-treatment data as fixture_staggered_result()
  cohort_data <- data.frame(
    cohort       = as.integer(rep(c(2005L, 2007L), each = 5L)),
    period       = as.integer(c(
      2003L, 2004L, 2005L, 2006L, 2007L,
      2005L, 2006L, 2007L, 2008L, 2009L
    )),
    att          = c(-0.1, -0.05, 0.5, 0.8, 1.0,
                     -0.2, -0.1,  0.6, 0.9, 1.2),
    se           = c( 0.2,  0.15, 0.3, 0.25, 0.2,
                      0.25, 0.2,  0.35, 0.3, 0.25),
    df_inference = rep(c(10L, 20L), each = 5L)
  )
  cohort_sizes <- c("2005" = 100, "2007" = 300)

  # Pre-treatment periods: event_times -4, -3 for both cohorts
  att_pre <- data.frame(
    cohort       = as.integer(rep(c(2005L, 2007L), each = 2L)),
    period       = as.integer(c(2001L, 2002L, 2003L, 2004L)),
    att          = c( 0.02, -0.01,  0.03, -0.02),
    se           = c( 0.10,  0.09,  0.12,  0.11),
    df_inference = rep(c(10L, 20L), each = 2L)
  )

  create_mock_staggered_result(
    cohort_data,
    cohort_sizes,
    att_pre_treatment    = att_pre,
    include_pretreatment = TRUE
  )
}


#' Single-cohort staggered result (cohort 2005 only)
fixture_staggered_result_single_cohort <- function() {
  cohort_data <- data.frame(
    cohort       = 2005L,
    period       = c(2003L, 2004L, 2005L, 2006L, 2007L),
    att          = c(-0.05, -0.02, 0.4, 0.7, 0.9),
    se           = c( 0.15,  0.10, 0.25, 0.20, 0.15),
    df_inference = 15L
  )
  cohort_sizes <- c("2005" = 200)

  create_mock_staggered_result(cohort_data, cohort_sizes)
}


#' All-zero ATT staggered result
fixture_staggered_result_zero_att <- function() {
  cohort_data <- data.frame(
    cohort       = as.integer(rep(c(2005L, 2007L), each = 5L)),
    period       = as.integer(c(
      2003L, 2004L, 2005L, 2006L, 2007L,
      2005L, 2006L, 2007L, 2008L, 2009L
    )),
    att          = rep(0, 10L),
    se           = rep(0.1, 10L),
    df_inference = 10L
  )
  cohort_sizes <- c("2005" = 100, "2007" = 100)

  create_mock_staggered_result(cohort_data, cohort_sizes)
}


#' Staggered result containing zero-SE rows
fixture_staggered_result_zero_se <- function() {
  cohort_data <- data.frame(
    cohort       = as.integer(rep(c(2005L, 2007L), each = 2L)),
    period       = as.integer(c(2005L, 2006L, 2007L, 2008L)),
    att          = c(0.5, 0.8, 0.6, 0.9),
    se           = c(0.0, 0.3, 0.2, 0.0),
    df_inference = 10L
  )
  cohort_sizes <- c("2005" = 100, "2007" = 100)

  create_mock_staggered_result(cohort_data, cohort_sizes)
}


#' Common Timing result
fixture_common_timing_result <- function() {
  period_data <- data.frame(
    period = c(1998L, 1999L, 2000L, 2001L, 2002L),
    att    = c(-0.1, -0.05, 0.5, 0.8, 1.0),
    se     = c( 0.2,  0.15, 0.3, 0.25, 0.2)
  )

  create_mock_common_timing_result(
    period_data,
    df_inference   = 50L,
    treatment_time = 2000L
  )
}
