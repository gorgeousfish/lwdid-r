# ============================================================================
# helper-visualization-export.R — Mock object factories for Epic 10
# comprehensive visualization/export/S3 method tests (Story 10.6).
#
# 7 factory functions: 4 lwdid_result variants + 3 diagnostic objects.
# All use fixed deterministic values for reproducibility.
# ============================================================================

# ── lwdid_result mock factories ─────────────────────────────────────────────

#' Common Timing lwdid_result mock
#' CI values: 0.5 ± qt(0.975, 48) * 0.12 = 0.5 ± 2.0106 * 0.12
.mock_common_timing_result <- function() {
  result <- list(
    att = 0.5,
    se_att = 0.12,
    t_stat = 4.167,
    pvalue = 0.0003,
    ci_lower = 0.259,
    ci_upper = 0.741,
    df_inference = 48L,
    nobs = 500L,
    n_treated = 50L,
    n_control = 450L,
    vce_type = "HC1",
    estimator = "RA",
    method = "common_timing",
    rolling = "demean",
    depvar = "y",
    K = 4L,
    tpost1 = 5L,
    is_staggered = FALSE,
    alpha = 0.05,
    control_group_used = "not_yet_treated",
    include_pretreatment = TRUE,
    att_by_period = data.frame(
      period = 5:8,
      att = c(0.45, 0.50, 0.55, 0.48),
      se = c(0.13, 0.12, 0.14, 0.11),
      t_stat = c(3.46, 4.17, 3.93, 4.36),
      pvalue = c(0.001, 0.0003, 0.0005, 0.0001),
      ci_lower = c(0.189, 0.259, 0.269, 0.259),
      ci_upper = c(0.711, 0.741, 0.831, 0.701)
    ),
    att_pre_treatment = data.frame(
      event_time = c(-4, -3, -2),
      att = c(0.02, -0.01, 0.03),
      se = c(0.10, 0.11, 0.09),
      ci_lower = c(-0.181, -0.231, -0.151),
      ci_upper = c(0.221, 0.211, 0.211),
      pvalue = c(0.84, 0.93, 0.74)
    ),
    parallel_trends_test = list(
      f_stat = 0.45, f_pvalue = 0.72, df1 = 3L, df2 = 45L
    ),
    ri_pvalue = NULL,
    ri_distribution = NULL,
    diagnostics = NULL,
    sensitivity = NULL,
    warnings_log = list(),
    metadata = list(
      call = quote(lwdid(y ~ d, data = df)),
      version = "0.1.0",
      parameters = list(),
      plot_data = NULL
    ),
    # Fields needed by S3 methods
    cluster_var = NULL,
    n_clusters = NULL,
    params = numeric(0),
    bse = numeric(0),
    vcov_matrix = NULL,
    att_by_cohort = NULL,
    att_by_cohort_time = NULL,
    cohort_weights = NULL,
    cohort_sample_sizes = NULL,
    event_time_effects = NULL,
    ri_seed = NULL,
    ri_method = NULL,
    ri_n_valid = NULL,
    rireps = NULL,
    vcov_att_periods = NULL,
    vcov_att_cohorts = NULL,
    vcov_full = NULL,
    control_group_auto_switched = FALSE
  )
  class(result) <- "lwdid_result"
  result
}

#' Staggered lwdid_result mock
.mock_staggered_result <- function() {
  result <- .mock_common_timing_result()
  result$method <- "staggered"
  result$is_staggered <- TRUE
  result$control_group_auto_switched <- FALSE
  result$att_by_cohort <- data.frame(
    cohort = c(5, 7),
    att = c(0.48, 0.53),
    se = c(0.14, 0.11),
    ci_lower = c(0.199, 0.309),
    ci_upper = c(0.761, 0.751),
    n = c(30L, 20L)
  )
  result$att_by_cohort_time <- data.frame(
    cohort = c(5, 5, 5, 7, 7),
    period = c(5, 6, 7, 7, 8),
    event_time = c(0, 1, 2, 0, 1),
    att = c(0.42, 0.50, 0.52, 0.55, 0.51),
    se = c(0.15, 0.13, 0.16, 0.12, 0.14),
    ci_lower = c(0.113, 0.234, 0.192, 0.298, 0.216),
    ci_upper = c(0.727, 0.766, 0.848, 0.802, 0.804),
    n_treated = c(30L, 30L, 30L, 20L, 20L),
    n_control = c(470L, 470L, 450L, 450L, 450L),
    df_inference = c(28L, 28L, 28L, 18L, 18L)
  )
  result$event_time_effects <- data.frame(
    event_time = 0:2,
    att = c(0.47, 0.505, 0.52),
    se = c(0.10, 0.09, 0.16),
    ci_lower = c(0.269, 0.324, 0.198),
    ci_upper = c(0.671, 0.686, 0.842),
    pvalue = c(0.0001, 0.00005, 0.001)
  )
  # Staggered mode: pre-treatment data should be cohort-specific
  # Remove inherited common-timing att_pre_treatment (no cohort column)
  result$att_pre_treatment <- data.frame(
    cohort = c(5, 5, 5, 7, 7),
    period = c(2, 3, 4, 4, 5),
    event_time = c(-3, -2, -1, -3, -2),
    att = c(0.02, -0.01, 0.0, 0.03, 0.0),
    se = c(0.10, 0.11, 0.0, 0.09, 0.0),
    ci_lower = c(-0.181, -0.231, 0.0, -0.151, 0.0),
    ci_upper = c(0.221, 0.211, 0.0, 0.211, 0.0),
    pvalue = c(0.84, 0.93, NA, 0.74, NA),
    df_inference = c(28L, 28L, 28L, 18L, 18L),
    data_source = rep("pre_treatment", 5)
  )
  result$cohort_weights <- c("5" = 0.6, "7" = 0.4)
  result$cohort_df <- c("5" = 28L, "7" = 18L)
  result$cohort_sample_sizes <- c("5" = 30L, "7" = 20L)
  result$cohort_sizes <- c("5" = 30L, "7" = 20L)
  result
}

#' lwdid_result with diagnostics
.mock_result_with_diagnostics <- function() {
  result <- .mock_common_timing_result()
  result$diagnostics <- list(
    clustering = structure(list(
      n_clusters = 25L,
      cluster_sizes = setNames(
        c(18L, 22L, 15L, 25L, 20L, 19L, 21L, 17L, 23L, 16L,
          24L, 20L, 18L, 22L, 19L, 21L, 17L, 23L, 20L, 16L,
          24L, 18L, 22L, 19L, 21L),
        paste0("c", 1:25)),
      balance_ratio = 3.2,
      icc = 0.15,
      effective_clusters = 22.5,
      reliability_level = "High"
    ), class = "lwdid_clustering_diagnosis"),
    selection = structure(list(
      attrition_rate = 0.05,
      selection_risk = "Low",
      attrition_by_period = data.frame(
        period = 1:8,
        attrition_rate = c(0.01, 0.02, 0.03, 0.04,
                           0.05, 0.06, 0.07, 0.08)
      )
    ), class = "lwdid_selection_diagnosis"),
    parallel_trends = structure(list(
      pre_treatment_coefficients = data.frame(
        period = 1:4,
        coefficient = c(0.02, -0.01, 0.03, -0.02),
        se = c(0.10, 0.11, 0.09, 0.12),
        ci_lower = c(-0.18, -0.23, -0.15, -0.26),
        ci_upper = c(0.22, 0.21, 0.21, 0.22)
      ),
      f_test = list(
        f_stat = 0.45, f_pvalue = 0.72, df1 = 3L, df2 = 45L
      ),
      group_trends = NULL
    ), class = "lwdid_trends")
  )
  result
}

#' lwdid_result with RI results
.mock_result_with_ri <- function() {
  result <- .mock_common_timing_result()
  result$ri_pvalue <- 0.004
  result$ri_distribution <- seq(-3, 3, length.out = 999)
  result$ri_seed <- 42L
  result$ri_method <- "permutation"
  result$ri_n_valid <- 999L
  result$rireps <- 999L
  result
}

# ── Diagnostic object mock factories ────────────────────────────────────────

#' Parallel trends diagnostic mock
.mock_trends_obj <- function() {
  structure(list(
    pre_treatment_coefficients = data.frame(
      period = 1:4,
      coefficient = c(0.02, -0.01, 0.03, -0.02),
      se = c(0.10, 0.11, 0.09, 0.12),
      ci_lower = c(-0.18, -0.23, -0.15, -0.26),
      ci_upper = c(0.22, 0.21, 0.21, 0.22)
    ),
    f_test = list(
      f_stat = 0.45, f_pvalue = 0.72, df1 = 3L, df2 = 45L
    ),
    group_trends = NULL
  ), class = "lwdid_trends")
}

#' Clustering diagnostic mock (fixed sizes, no randomness)
.mock_clustering_obj <- function() {
  fixed_sizes <- c(18L, 22L, 15L, 25L, 20L, 19L, 21L, 17L,
                   23L, 16L, 24L, 20L, 18L, 22L, 19L, 21L,
                   17L, 23L, 20L, 16L, 24L, 18L, 22L, 19L,
                   21L)
  structure(list(
    n_clusters = 25L,
    cluster_sizes = setNames(fixed_sizes, paste0("c", 1:25)),
    balance_ratio = 3.2,
    icc = 0.15,
    effective_clusters = 22.5,
    reliability_level = "High"
  ), class = "lwdid_clustering_diagnosis")
}

#' Selection diagnostic mock
.mock_selection_obj <- function() {
  structure(list(
    attrition_rate = 0.05,
    selection_risk = "Low",
    attrition_by_period = data.frame(
      period = 1:8,
      attrition_rate = c(0.01, 0.02, 0.03, 0.04,
                         0.05, 0.06, 0.07, 0.08)
    )
  ), class = "lwdid_selection_diagnosis")
}
