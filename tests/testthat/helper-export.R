# helper-export.R
# Mock object factories for E10-03 export tests.

# TC-10.3.0a: Basic common timing lwdid_result mock
.mock_common_timing_result <- function(
    att = 1.5, se_att = 0.3, t_stat = 5.0, pvalue = 0.001,
    ci_lower = 0.9, ci_upper = 2.1, alpha = 0.05,
    nobs = 500L, n_treated = 200L, n_control = 300L,
    estimator = "DR", vce_type = "HC1", rolling = TRUE,
    ri_pvalue = NULL, ri_seed = NULL, rireps = NULL,
    diagnostics = NULL, att_by_period = NULL,
    is_staggered = FALSE, att_by_cohort = NULL,
    cohort_weights = NULL, att_by_cohort_time = NULL
) {
  obj <- list(
    att = att, se_att = se_att, t_stat = t_stat, pvalue = pvalue,
    ci_lower = ci_lower, ci_upper = ci_upper, alpha = alpha,
    nobs = nobs, n_treated = n_treated, n_control = n_control,
    estimator = estimator, vce_type = vce_type, rolling = rolling,
    ri_pvalue = ri_pvalue, ri_seed = ri_seed, rireps = rireps,
    diagnostics = diagnostics, att_by_period = att_by_period,
    is_staggered = is_staggered, att_by_cohort = att_by_cohort,
    cohort_weights = cohort_weights, att_by_cohort_time = att_by_cohort_time
  )
  class(obj) <- "lwdid_result"
  obj
}

# TC-10.3.0b: Staggered lwdid_result mock
.mock_staggered_result <- function() {
  att_by_cohort <- data.frame(
    cohort = c(2004L, 2006L, 2008L),
    att = c(1.2, 1.8, 1.5),
    se = c(0.25, 0.30, 0.28),
    n = c(100L, 150L, 120L),
    weight = c(0.30, 0.40, 0.30),
    pvalue = c(0.001, 0.0001, 0.002),
    stringsAsFactors = FALSE
  )
  att_by_cohort_time <- data.frame(
    cohort = rep(c(2004L, 2006L, 2008L), each = 3),
    time = rep(1:3, 3),
    att = rnorm(9, 1.5, 0.3),
    se = rep(0.3, 9),
    stringsAsFactors = FALSE
  )
  .mock_common_timing_result(
    is_staggered = TRUE,
    att_by_cohort = att_by_cohort,
    cohort_weights = c("2004" = 0.30, "2006" = 0.40, "2008" = 0.30),
    att_by_cohort_time = att_by_cohort_time
  )
}

# Freeze export-specific factories so later helper files with the same
# symbol names cannot overwrite the objects used by test-export.R.
.mock_export_common_timing_result <- local({
  factory <- .mock_common_timing_result
  function(...) factory(...)
})

.mock_export_staggered_result <- local({
  factory_ct <- .mock_common_timing_result
  function() {
    att_by_cohort <- data.frame(
      cohort = c(2004L, 2006L, 2008L),
      att = c(1.2, 1.8, 1.5),
      se = c(0.25, 0.30, 0.28),
      n = c(100L, 150L, 120L),
      weight = c(0.30, 0.40, 0.30),
      pvalue = c(0.001, 0.0001, 0.002),
      stringsAsFactors = FALSE
    )
    att_by_cohort_time <- data.frame(
      cohort = rep(c(2004L, 2006L, 2008L), each = 3),
      time = rep(1:3, 3),
      att = stats::rnorm(9, 1.5, 0.3),
      se = rep(0.3, 9),
      stringsAsFactors = FALSE
    )
    factory_ct(
      is_staggered = TRUE,
      att_by_cohort = att_by_cohort,
      cohort_weights = c("2004" = 0.30, "2006" = 0.40, "2008" = 0.30),
      att_by_cohort_time = att_by_cohort_time
    )
  }
})
