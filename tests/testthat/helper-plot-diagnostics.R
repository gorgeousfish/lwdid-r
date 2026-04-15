# helper-plot-diagnostics.R
# Mock object factories for E10-02 diagnostic plot tests.
# Each factory produces an S3 object matching the class structure
# expected by the corresponding plot functions.

create_mock_sensitivity_pre_period <- function(
    specifications = NULL,
    baseline_att = 2.0,
    robustness_threshold = 0.10,
    sensitivity_ratio = 0.85,
    robustness_level = "Robust"
) {
  if (is.null(specifications)) {
    specifications <- list(
      list(n_pre_periods = 2L, att = 1.8, se = 0.3, pvalue = 0.001,
           ci_lower = 1.2, ci_upper = 2.4, converged = TRUE),
      list(n_pre_periods = 3L, att = 2.1, se = 0.25, pvalue = 0.0001,
           ci_lower = 1.6, ci_upper = 2.6, converged = TRUE),
      list(n_pre_periods = 4L, att = 1.9, se = 0.35, pvalue = 0.002,
           ci_lower = 1.2, ci_upper = 2.6, converged = TRUE),
      list(n_pre_periods = 5L, att = 2.3, se = 0.4, pvalue = 0.001,
           ci_lower = 1.5, ci_upper = 3.1, converged = TRUE)
    )
  }
  obj <- list(
    type = "pre_period",
    specifications = specifications,
    baseline_spec = list(att = baseline_att),
    robustness_threshold = robustness_threshold,
    sensitivity_ratio = sensitivity_ratio,
    robustness_level = robustness_level
  )
  class(obj) <- "lwdid_sensitivity"
  obj
}

create_mock_sensitivity_no_anticipation <- function(
    estimates = NULL,
    baseline_att = 2.0,
    anticipation_detected = TRUE,
    recommended_exclusion = 2L
) {
  if (is.null(estimates)) {
    estimates <- list(
      list(excluded_periods = 0L, att = 2.0, se = 0.3,
           ci_lower = 1.4, ci_upper = 2.6, converged = TRUE),
      list(excluded_periods = 1L, att = 2.2, se = 0.28,
           ci_lower = 1.64, ci_upper = 2.76, converged = TRUE),
      list(excluded_periods = 2L, att = 1.9, se = 0.32,
           ci_lower = 1.26, ci_upper = 2.54, converged = TRUE),
      list(excluded_periods = 3L, att = 2.1, se = 0.35,
           ci_lower = 1.4, ci_upper = 2.8, converged = TRUE)
    )
  }
  obj <- list(
    type = "no_anticipation",
    estimates = estimates,
    baseline_estimate = list(att = baseline_att),
    anticipation_detected = anticipation_detected,
    recommended_exclusion = recommended_exclusion
  )
  class(obj) <- "lwdid_sensitivity"
  obj
}

create_mock_sensitivity_comprehensive <- function() {
  obj <- list(
    pre_period = create_mock_sensitivity_pre_period(),
    no_anticipation = create_mock_sensitivity_no_anticipation()
  )
  class(obj) <- "lwdid_sensitivity_comprehensive"
  obj
}

create_mock_heterogeneous_trends <- function(
    trend_by_cohort = NULL,
    control_group_trend = NULL,
    has_heterogeneous_trends = TRUE,
    trend_heterogeneity_test = NULL
) {
  if (is.null(trend_by_cohort)) {
    trend_by_cohort <- list(
      list(cohort = 2004L, intercept = 10.0, slope = 0.5,
           slope_se = 0.1, slope_pvalue = 0.001,
           n_units = 50L, n_pre_periods = 4L,
           r_squared = 0.85, residual_std = 0.5),
      list(cohort = 2006L, intercept = 12.0, slope = 0.3,
           slope_se = 0.15, slope_pvalue = 0.08,
           n_units = 40L, n_pre_periods = 6L,
           r_squared = 0.72, residual_std = 0.6)
    )
  }
  if (is.null(control_group_trend)) {
    control_group_trend <- list(
      intercept = 8.0, slope = 0.2, slope_se = 0.05,
      slope_pvalue = 0.001, n_units = 100L
    )
  }
  if (is.null(trend_heterogeneity_test)) {
    trend_heterogeneity_test <- list(
      f_stat = 3.45, f_pvalue = 0.035,
      df1 = 1L, df2 = 88L
    )
  }
  obj <- list(
    trend_by_cohort = trend_by_cohort,
    control_group_trend = control_group_trend,
    has_heterogeneous_trends = has_heterogeneous_trends,
    trend_heterogeneity_test = trend_heterogeneity_test
  )
  class(obj) <- "lwdid_heterogeneous_trends"
  obj
}

create_mock_transformation_recommendation <- function(
    scores = NULL,
    recommended_method = "demeanq",
    confidence_level = "High",
    sub_scores = NULL,
    method_labels = NULL
) {
  if (is.null(scores)) {
    scores <- c(demean = 0.85, detrend = 0.72,
                demeanq = 0.90, detrendq = 0.65)
  }
  obj <- list(
    scores = scores,
    recommended_method = recommended_method,
    confidence_level = confidence_level,
    sub_scores = sub_scores,
    method_labels = method_labels
  )
  class(obj) <- "lwdid_transformation_recommendation"
  obj
}
