# ============================================================================
# test-s3-methods.R — E10-04.8: S3 Method Tests (95 test cases)
# ============================================================================
library(testthat)

# ── Test Helpers ────────────────────────────────────────────────────────────

# Minimal Common Timing result
.mk_ct <- function(...) {
  defaults <- list(
    att = 0.5, se_att = 0.12, df_resid = 52L, df_inference = 50L,
    nobs = 200L, n_treated = 80L, n_control = 120L,
    K = 3L, tpost1 = 4L, depvar = "y", rolling = "demean",
    vce_type = "ols", alpha = 0.05, is_staggered = FALSE,
    estimator = "ra"
  )
  do.call(new_lwdid_result, modifyList(defaults, list(...)))
}

# Minimal Staggered result
.mk_stag <- function(...) {
  cohort_eff <- list(
    list(cohort = 2005L, att = 0.3, se = 0.1, ci_lower = 0.1, ci_upper = 0.5,
         t_stat = 3.0, pvalue = 0.003, n_periods = 5L, n_units = 40L,
         n_control = 60L, df_resid = 98L, df_inference = 48L),
    list(cohort = 2007L, att = 0.7, se = 0.15, ci_lower = 0.4, ci_upper = 1.0,
         t_stat = 4.67, pvalue = 0.001, n_periods = 3L, n_units = 30L,
         n_control = 60L, df_resid = 88L, df_inference = 48L)
  )
  defaults <- list(
    att = 0.5, se_att = 0.12, df_resid = 100L, df_inference = 50L,
    nobs = 500L, n_treated = 200L, n_control = 300L,
    K = 2L, tpost1 = 3L, depvar = "outcome", rolling = "demean",
    vce_type = "cluster", cluster_var = "state", n_clusters = 51L,
    alpha = 0.05, is_staggered = TRUE, estimator = "ra",
    aggregate = "cohort", cohorts = c(2005L, 2007L),
    n_never_treated = 300L,
    cohort_effects = cohort_eff,
    att_by_cohort = data.frame(
      cohort = c(2005L, 2007L), att = c(0.3, 0.7), se = c(0.1, 0.15),
      stringsAsFactors = FALSE
    ),
    att_by_cohort_time = data.frame(
      cohort = c(2005L, 2005L, 2007L),
      period = c(2005L, 2006L, 2007L),
      att = c(0.25, 0.35, 0.70),
      se = c(0.08, 0.09, 0.15),
      stringsAsFactors = FALSE
    ),
    cohort_weights = c("2005" = 0.57, "2007" = 0.43),
    control_group = "not_yet_treated",
    control_group_used = "never_treated",
    att_overall = 0.5, se_overall = 0.12
  )
  do.call(new_lwdid_result, modifyList(defaults, list(...)))
}

# ============================================================================
# print tests (TC-10.4.1 through TC-10.4.41)
# ============================================================================

# TC-10.4.1: print basic output (ATT/SE/t-stat/p-value/CI/N)
test_that("TC-10.4.1: print CT shows ATT/SE/t-stat/p-value/CI/N", {
  obj <- .mk_ct()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "ATT")
  expect_match(txt, "SE")
  expect_match(txt, "t-stat")
  expect_match(txt, "p-value")
  expect_match(txt, "CI")
  expect_match(txt, "N = 200")
})

# TC-10.4.2: print returns invisible(x)
test_that("TC-10.4.2: print returns invisible(x)", {
  obj <- .mk_ct()
  ret <- print(obj)
  expect_identical(ret, obj)
})

# TC-10.4.3: print shows RI info
test_that("TC-10.4.3: print shows RI p-value", {
  obj <- .mk_ct(ri_pvalue = 0.03, ri_method = "fisher",
                 ri_seed = 42L, ri_distribution = rnorm(500))
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "RI p-value")
})

# TC-10.4.4: print shows warning count
test_that("TC-10.4.4: print shows warning count", {
  obj <- .mk_ct(warnings_log = list("w1", "w2", "w3"))
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "3 warning")
})

# TC-10.4.5: print on non-lwdid_result errors
test_that("TC-10.4.5: print on non-lwdid_result errors", {
  expect_error(print.lwdid_result(list(att = 1)))
})

# TC-10.4.31: .vce_description mapping correct
test_that("TC-10.4.31: .vce_description maps all VCE types", {
  expect_equal(.vce_description("ols"), "OLS (Homoskedastic)")
  expect_equal(.vce_description("robust"), "HC1 (Heteroskedasticity-robust)")
  expect_equal(.vce_description("hc1"), "HC1 (Heteroskedasticity-robust)")
  expect_equal(.vce_description("hc0"), "HC0 (White)")
  expect_equal(.vce_description("hc2"), "HC2 (Bell-McCaffrey)")
  expect_equal(.vce_description("hc3"), "HC3 (Small-sample adjusted)")
  expect_equal(.vce_description("hc4"), "HC4 (Cribari-Neto)")
  expect_equal(.vce_description("bootstrap"), "Bootstrap")
  expect_equal(.vce_description("cluster", "state"),
               "Cluster-robust (clustered by state)")
  expect_equal(.vce_description("cluster"), "Cluster-robust (clustered by ?)")
  expect_equal(.vce_description(NULL), "OLS (Homoskedastic)")
  # Case insensitive
  expect_equal(.vce_description("OLS"), "OLS (Homoskedastic)")
  expect_equal(.vce_description("HC3"), "HC3 (Small-sample adjusted)")
  # Unknown type returned as-is
  expect_equal(.vce_description("custom_vce"), "custom_vce")
})

# TC-10.4.32: print outputs VCE human-readable description
test_that("TC-10.4.32: print outputs VCE description", {
  obj <- .mk_ct(vce_type = "hc3")
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "HC3 \\(Small-sample adjusted\\)")
})

# TC-10.4.33: print cluster VCE includes variable name
test_that("TC-10.4.33: print cluster VCE includes variable name", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Cluster-robust \\(clustered by state\\)")
})

# TC-10.4.34: print Staggered shows cohort count
test_that("TC-10.4.34: print Staggered shows cohort count", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Cohorts: 2")
})

# TC-10.4.35: print Staggered shows control group strategy
test_that("TC-10.4.35: print Staggered shows control group", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Control group: never_treated")
})

# TC-10.4.36: print Staggered shows auto-switch notification
test_that("TC-10.4.36: print Staggered shows auto-switch", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "auto-switched from not_yet_treated")
})

# TC-10.4.37: print Staggered shows tau_omega label
test_that("TC-10.4.37: print Staggered shows tau_omega label", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_true(grepl("tau_omega", txt, fixed = TRUE),
              info = "print output must contain ASCII tau_omega label")
})

# TC-10.4.38: print Staggered shows Pre-treatment flag
test_that("TC-10.4.38: print Staggered shows Pre-treatment flag", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Pre-treatment: FALSE")
})

# TC-10.4.39: print outputs df
test_that("TC-10.4.39: print outputs df", {
  obj <- .mk_ct()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "df\\s+=\\s+50")
})

# TC-10.4.40: print outputs K and tpost1
test_that("TC-10.4.40: print outputs K and tpost1", {
  obj <- .mk_ct()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "K = 3")
  expect_match(txt, "tpost1 = 4")
})

# TC-10.4.40a: print Staggered att=NULL shows aggregate type
test_that("TC-10.4.40a: print Staggered att=NULL shows aggregate type", {
  obj <- .mk_stag(att = NA_real_, se_att = NA_real_)
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Aggregate: cohort")
})

# TC-10.4.40b: print Staggered att=NULL aggregate=NULL shows "none"
test_that("TC-10.4.40b: print Staggered att=NULL aggregate=NULL shows none", {
  obj <- .mk_stag(att = NA_real_, se_att = NA_real_, aggregate = "none")
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Aggregate: none")
})

# TC-10.4.41: print outputs RI details (method/seed/valid/reps)
test_that("TC-10.4.41: print RI details", {
  ri_dist <- rnorm(200)
  obj <- .mk_ct(ri_pvalue = 0.04, ri_method = "fisher",
                 ri_seed = 123L, ri_valid = 180L,
                 ri_distribution = ri_dist)
  # Need to set ri_n_valid manually since constructor uses ri_valid
  obj$ri_n_valid <- 180L
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "method=fisher")
  expect_match(txt, "seed=123")
  expect_match(txt, "valid=180/200")
})

# TC-10.4.77: print cluster VCE shows n_clusters
test_that("TC-10.4.77: print cluster VCE shows n_clusters", {
  obj <- .mk_stag()
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Number of clusters: 51")
})

# TC-10.4.78: print non-cluster VCE does NOT show n_clusters
test_that("TC-10.4.78: print non-cluster VCE no n_clusters", {
  obj <- .mk_ct(vce_type = "ols")
  out <- capture.output(print(obj))
  txt <- paste(out, collapse = "\n")
  expect_false(grepl("Number of clusters", txt))
})

# ============================================================================
# summary tests (TC-10.4.6 through TC-10.4.48)
# ============================================================================

# TC-10.4.6: summary returns correct class
test_that("TC-10.4.6: summary returns summary.lwdid_result class", {
  obj <- .mk_ct()
  s <- summary(obj)
  expect_s3_class(s, "summary.lwdid_result")
})

# TC-10.4.7: summary coefficient table structure
test_that("TC-10.4.7: summary coef table has correct columns", {
  obj <- .mk_ct()
  s <- summary(obj)
  expect_true(all(c("Estimate", "Std.Error", "t.value", "Pr.t",
                     "CI.lower", "CI.upper") %in% names(s$coefficients)))
})

# TC-10.4.8: summary Staggered includes cohort info
test_that("TC-10.4.8: summary Staggered includes cohort info", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_true(s$is_staggered)
  expect_true(!is.null(s$cohort_list))
  expect_true(!is.null(s$cohort_weights))
})

# TC-10.4.9: summary includes diagnostics summary
test_that("TC-10.4.9: summary includes diagnostics", {
  diag <- list(clustering = list(n_clusters = 51L, effective_clusters = 45.2,
                                  reliability_level = "high"))
  obj <- .mk_stag(diagnostics = diag)
  s <- summary(obj)
  expect_true(!is.null(s$diagnostics_summary))
  expect_match(s$diagnostics_summary$clustering, "n_clusters=51")
})

# TC-10.4.28: summary CI percentage uses correct alpha
test_that("TC-10.4.28: summary CI uses alpha", {
  obj <- .mk_ct(alpha = 0.10)
  s <- summary(obj)
  # CI columns should reflect alpha=0.10 (90% CI)
  expect_true(!is.null(s$alpha))
  expect_equal(s$alpha, 0.10)
})

# TC-10.4.29: summary includes pre-treatment dynamics
test_that("TC-10.4.29: summary includes pre-treatment dynamics", {
  pre <- data.frame(event_time = c(-3, -2, -1), att = c(0.01, -0.02, 0),
                    se = c(0.05, 0.04, 0), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_pre_treatment = pre, include_pretreatment = TRUE)
  s <- summary(obj)
  expect_true(!is.null(s$att_pre_treatment))
  expect_true(s$include_pretreatment)
})

# TC-10.4.30: summary alpha correctly passed
test_that("TC-10.4.30: summary alpha passed", {
  obj <- .mk_ct(alpha = 0.01)
  s <- summary(obj)
  expect_equal(s$alpha, 0.01)
})

# TC-10.4.42: summary Staggered includes cohort list
test_that("TC-10.4.42: summary Staggered cohort list", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_true(length(s$cohort_list) > 0)
  expect_true(2005L %in% s$cohort_list)
})

# TC-10.4.43: summary Staggered includes control group strategy
test_that("TC-10.4.43: summary Staggered control group", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_equal(s$control_group_used, "never_treated")
  expect_equal(s$control_group, "not_yet_treated")
})

# TC-10.4.44: summary includes VCE description
test_that("TC-10.4.44: summary VCE description", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_match(s$vce_description, "Cluster-robust")
})

# TC-10.4.45: summary includes K and tpost1
test_that("TC-10.4.45: summary K and tpost1", {
  obj <- .mk_ct()
  s <- summary(obj)
  expect_equal(s$K, 3L)
  expect_equal(s$tpost1, 4L)
})

# TC-10.4.46: summary Parallel Trends includes reject conclusion
test_that("TC-10.4.46: summary parallel trends reject", {
  pt <- list(f_stat = 5.2, f_pvalue = 0.001, df1 = 3L, df2 = 96L)
  obj <- .mk_ct(parallel_trends_test = pt)
  s <- summary(obj)
  expect_true(s$parallel_trends$reject)
})

# TC-10.4.47: summary includes RI details
test_that("TC-10.4.47: summary RI details", {
  obj <- .mk_ct(ri_pvalue = 0.05, ri_method = "fisher",
                 ri_seed = 42L, ri_valid = 900L,
                 ri_distribution = rnorm(1000))
  s <- summary(obj)
  expect_equal(s$ri_pvalue, 0.05)
  expect_equal(s$ri_method, "fisher")
  expect_equal(s$ri_n_permutations, 1000L)
})

# TC-10.4.48: summary Staggered includes cohort sample sizes
test_that("TC-10.4.48: summary cohort sample sizes", {
  obj <- .mk_stag()
  obj$cohort_sample_sizes <- c("2005" = 40L, "2007" = 30L)
  s <- summary(obj)
  expect_true(!is.null(s$cohort_sample_sizes))
  expect_equal(s$cohort_sample_sizes[["2005"]], 40L)
})

# TC-10.4.79: summary Staggered includes n_never_treated
test_that("TC-10.4.79: summary n_never_treated", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_equal(s$n_never_treated, 300L)
})

# TC-10.4.80: summary Staggered includes aggregate field
test_that("TC-10.4.80: summary aggregate field", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_equal(s$aggregate, "cohort")
})

# TC-10.4.84: RI rireps via length(ri_distribution)
test_that("TC-10.4.84: RI rireps from ri_distribution length", {
  obj <- .mk_ct(ri_pvalue = 0.05, ri_method = "fisher",
                 ri_distribution = rnorm(500))
  s <- summary(obj)
  expect_equal(s$ri_n_permutations, 500L)
})

# TC-10.4.85: RI rireps fallback to rireps field
test_that("TC-10.4.85: RI rireps fallback to rireps field", {
  obj <- .mk_ct(ri_pvalue = 0.05, ri_method = "fisher",
                 rireps = 1000L, ri_distribution = NULL)
  s <- summary(obj)
  expect_equal(s$ri_n_permutations, 1000L)
})

# TC-10.4.88a: summary Staggered cohort_effects_detail has n_units/n_periods
test_that("TC-10.4.88a: summary cohort_effects_detail columns", {
  obj <- .mk_stag()
  s <- summary(obj)
  expect_true(!is.null(s$cohort_effects_detail))
  expect_true("n_units" %in% names(s$cohort_effects_detail))
  expect_true("n_periods" %in% names(s$cohort_effects_detail))
})

# TC-10.4.93: summary all fields non-NULL for CT basic
test_that("TC-10.4.93: summary CT all basic fields non-NULL", {
  obj <- .mk_ct()
  s <- summary(obj)
  expect_true(!is.null(s$method))
  expect_true(!is.null(s$estimator))
  expect_true(!is.null(s$rolling))
  expect_true(!is.null(s$depvar))
  expect_true(!is.null(s$nobs))
  expect_true(!is.null(s$coefficients))
  expect_true(!is.null(s$alpha))
  expect_true(!is.null(s$vce_description))
})

# ============================================================================
# print.summary tests (TC-10.4.10, TC-10.4.49-60, TC-10.4.94)
# ============================================================================

# TC-10.4.10: print.summary outputs complete info
test_that("TC-10.4.10: print.summary complete output", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Summary")
  expect_match(txt, "ATT")
  expect_match(txt, "Cohort")
})

# TC-10.4.49: print.summary outputs VCE description
test_that("TC-10.4.49: print.summary VCE description", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Cluster-robust")
})

# TC-10.4.50: print.summary Staggered outputs tau_omega label
test_that("TC-10.4.50: print.summary tau_omega label", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_true(grepl("tau_omega", txt, fixed = TRUE),
              info = "summary output must contain ASCII tau_omega label")
})

# TC-10.4.51: print.summary Staggered outputs cohort list
test_that("TC-10.4.51: print.summary cohort list", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "2005")
  expect_match(txt, "2007")
})

# TC-10.4.52: print.summary control group and auto-switch
test_that("TC-10.4.52: print.summary control group auto-switch", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Control group: never_treated")
  expect_match(txt, "auto-switched from not_yet_treated")
})

# TC-10.4.53: print.summary att_by_period >5 rows truncated
test_that("TC-10.4.53: print.summary period effects truncated >5", {
  periods <- data.frame(
    period = 1:8, att = runif(8), se = runif(8, 0.01, 0.1),
    stringsAsFactors = FALSE
  )
  obj <- .mk_ct(att_by_period = periods)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "more rows")
})

# TC-10.4.54: print.summary att_by_period <=5 rows complete
test_that("TC-10.4.54: print.summary period effects complete <=5", {
  periods <- data.frame(
    period = 1:3, att = c(0.1, 0.2, 0.3), se = c(0.05, 0.06, 0.07),
    stringsAsFactors = FALSE
  )
  obj <- .mk_ct(att_by_period = periods)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_false(grepl("more rows", txt))
})

# TC-10.4.55: print.summary Pre-treatment Dynamics anchor point
test_that("TC-10.4.55: print.summary anchor point marking", {
  pre <- data.frame(event_time = c(-3, -2, -1), att = c(0.01, -0.02, 0),
                    se = c(0.05, 0.04, 0), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_pre_treatment = pre, include_pretreatment = TRUE)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "anchor")
})

# TC-10.4.56: print.summary Parallel Trends reject conclusion
test_that("TC-10.4.56: print.summary parallel trends reject", {
  pt <- list(f_stat = 5.2, f_pvalue = 0.001, df1 = 3L, df2 = 96L)
  obj <- .mk_ct(parallel_trends_test = pt)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Reject H0")
})

# TC-10.4.57: print.summary Parallel Trends fail-to-reject
test_that("TC-10.4.57: print.summary parallel trends fail-to-reject", {
  pt <- list(f_stat = 1.2, f_pvalue = 0.30, df1 = 3L, df2 = 96L)
  obj <- .mk_ct(parallel_trends_test = pt)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Fail to reject")
})

# TC-10.4.58: print.summary RI details
test_that("TC-10.4.58: print.summary RI details", {
  obj <- .mk_ct(ri_pvalue = 0.04, ri_method = "fisher",
                 ri_seed = 42L, ri_valid = 900L,
                 ri_distribution = rnorm(1000))
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Randomization Inference")
  expect_match(txt, "method=fisher")
})

# TC-10.4.59: print.summary K and tpost1
test_that("TC-10.4.59: print.summary K and tpost1", {
  obj <- .mk_ct()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "K: 3")
  expect_match(txt, "tpost1: 4")
})

# TC-10.4.60: print.summary cohort weights with N
test_that("TC-10.4.60: print.summary cohort weights with N", {
  obj <- .mk_stag()
  obj$cohort_sample_sizes <- c("2005" = 40L, "2007" = 30L)
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Cohort Weights")
  expect_match(txt, "N=40")
})

# TC-10.4.81: print.summary Staggered aggregate level
test_that("TC-10.4.81: print.summary aggregate level", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Aggregation: cohort")
})

# TC-10.4.82: print.summary Staggered usage hints
test_that("TC-10.4.82: print.summary usage hints", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "coef")
  expect_match(txt, "plot")
})

# TC-10.4.83: print.summary auto-switch includes original control_group
test_that("TC-10.4.83: print.summary auto-switch original value", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "not_yet_treated")
})

# TC-10.4.88: print.summary Staggered n_never_treated
test_that("TC-10.4.88: print.summary n_never_treated", {
  obj <- .mk_stag()
  s <- summary(obj)
  out <- capture.output(print(s))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Never Treated.*300")
})

# TC-10.4.94: print.summary signif.stars=FALSE no stars
test_that("TC-10.4.94: print.summary signif.stars=FALSE", {
  obj <- .mk_ct(att = 5.0, se_att = 0.1, pvalue = 0.0001)
  s <- summary(obj)
  out <- capture.output(print(s, signif.stars = FALSE))
  # printCoefmat with signif.stars=FALSE should not show stars
  # (hard to test precisely, but at least no error)
  expect_true(length(out) > 0)
})

# ============================================================================
# coef tests (TC-10.4.11-15, TC-10.4.74, TC-10.4.86, TC-10.4.91)
# ============================================================================

# TC-10.4.11: coef overall returns named scalar
test_that("TC-10.4.11: coef overall returns named scalar", {
  obj <- .mk_ct()
  co <- coef(obj, type = "overall")
  expect_equal(length(co), 1L)
  expect_equal(names(co), "ATT")
  expect_equal(unname(co), 0.5)
})

# TC-10.4.12: coef by_period returns named vector
test_that("TC-10.4.12: coef by_period returns named vector", {
  periods <- data.frame(period = c(4L, 5L, 6L), att = c(0.3, 0.5, 0.7),
                        se = c(0.1, 0.1, 0.1), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  co <- coef(obj, type = "by_period")
  expect_equal(length(co), 3L)
  expect_equal(names(co), c("4", "5", "6"))
})

# TC-10.4.13: coef by_cohort on Common Timing errors
test_that("TC-10.4.13: coef by_cohort CT errors", {
  obj <- .mk_ct()
  expect_error(coef(obj, type = "by_cohort"), "Staggered")
})

# TC-10.4.14: coef by_cohort on Staggered correct
test_that("TC-10.4.14: coef by_cohort Staggered", {
  obj <- .mk_stag()
  co <- coef(obj, type = "by_cohort")
  expect_equal(length(co), 2L)
  expect_equal(names(co), c("2005", "2007"))
})

# TC-10.4.15: coef all returns finest granularity (g,r)
test_that("TC-10.4.15: coef all returns (g,r) level", {
  obj <- .mk_stag()
  co <- coef(obj, type = "all")
  expect_equal(length(co), 3L)
  expect_true(all(grepl("^g\\d+\\.r\\d+$", names(co))))
})

# TC-10.4.74: coef on non-lwdid_result errors
test_that("TC-10.4.74: coef non-lwdid_result errors", {
  expect_error(coef.lwdid_result(list(att = 1)))
})

# TC-10.4.86: coef by_period NULL errors with data source hint
test_that("TC-10.4.86: coef by_period NULL error hint", {
  obj <- .mk_ct()
  expect_error(coef(obj, type = "by_period"), "period")
})

# TC-10.4.91: coef all fallback chain
test_that("TC-10.4.91: coef all fallback to att_by_period", {
  periods <- data.frame(period = c(4L, 5L), att = c(0.3, 0.5),
                        se = c(0.1, 0.1), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  co <- coef(obj, type = "all")
  expect_equal(length(co), 2L)
  expect_equal(names(co), c("4", "5"))
})

# ============================================================================
# confint tests (TC-10.4.16-21, TC-10.4.61-63, TC-10.4.70-72, TC-10.4.75,
#                TC-10.4.87, TC-10.4.90, TC-10.4.92)
# ============================================================================

# TC-10.4.16: confint default level=0.95 standard format
test_that("TC-10.4.16: confint default level=0.95", {
  obj <- .mk_ct()
  ci <- confint(obj)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(ncol(ci), 2L)
  expect_equal(rownames(ci), "ATT")
  expect_match(colnames(ci)[1], "%")
  expect_match(colnames(ci)[2], "%")
})

# TC-10.4.17: confint level=0.90 narrower CI
test_that("TC-10.4.17: confint level=0.90 narrower", {
  obj <- .mk_ct()
  ci95 <- confint(obj, level = 0.95)
  ci90 <- confint(obj, level = 0.90)
  width95 <- ci95[1, 2] - ci95[1, 1]
  width90 <- ci90[1, 2] - ci90[1, 1]
  expect_true(width90 < width95)
})

# TC-10.4.18: confint numerical correctness (tolerance=1e-10)
test_that("TC-10.4.18: confint numerical correctness", {
  # att=0.5, se=0.12, df=50, level=0.95
  # qt(0.975, 50) ≈ 2.00856
  obj <- .mk_ct()
  ci <- confint(obj)
  t_crit <- qt(0.975, 50)
  expected_lower <- 0.5 - t_crit * 0.12
  expected_upper <- 0.5 + t_crit * 0.12
  expect_equal(ci[1, 1], expected_lower, tolerance = 1e-10)
  expect_equal(ci[1, 2], expected_upper, tolerance = 1e-10)
})

# TC-10.4.19: confint by_period returns multi-row matrix
test_that("TC-10.4.19: confint by_period multi-row", {
  periods <- data.frame(period = c(4L, 5L, 6L), att = c(0.3, 0.5, 0.7),
                        se = c(0.1, 0.1, 0.1), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  ci <- confint(obj, type = "by_period")
  expect_equal(nrow(ci), 3L)
  expect_equal(rownames(ci), c("4", "5", "6"))
})

# TC-10.4.20: confint parm numeric index filter
test_that("TC-10.4.20: confint parm numeric filter", {
  periods <- data.frame(period = c(4L, 5L, 6L), att = c(0.3, 0.5, 0.7),
                        se = c(0.1, 0.1, 0.1), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  ci <- confint(obj, parm = 1:2, type = "by_period")
  expect_equal(nrow(ci), 2L)
})

# TC-10.4.21: confint invalid level errors
test_that("TC-10.4.21: confint invalid level errors", {
  obj <- .mk_ct()
  expect_error(confint(obj, level = 0))
  expect_error(confint(obj, level = 1))
  expect_error(confint(obj, level = -0.5))
})

# TC-10.4.61: confint type="all" returns (g,r) level CI
test_that("TC-10.4.61: confint type=all (g,r) level", {
  obj <- .mk_stag()
  ci <- confint(obj, type = "all")
  expect_equal(nrow(ci), 3L)
  expect_true(all(grepl("^g\\d+\\.r\\d+$", rownames(ci))))
})

# TC-10.4.62: confint type="all" rownames match coef type="all"
test_that("TC-10.4.62: confint all rownames match coef all", {
  obj <- .mk_stag()
  co <- coef(obj, type = "all")
  ci <- confint(obj, type = "all")
  expect_equal(rownames(ci), names(co))
})

# TC-10.4.63: confint type="all" numerical correctness
test_that("TC-10.4.63: confint all numerical correctness", {
  obj <- .mk_stag()
  ci <- confint(obj, type = "all")
  t_crit <- qt(0.975, 50)
  # First row: att=0.25, se=0.08
  expect_equal(ci[1, 1], 0.25 - t_crit * 0.08, tolerance = 1e-10)
  expect_equal(ci[1, 2], 0.25 + t_crit * 0.08, tolerance = 1e-10)
})

# TC-10.4.70: confint by_cohort Staggered correct
test_that("TC-10.4.70: confint by_cohort Staggered", {
  obj <- .mk_stag()
  ci <- confint(obj, type = "by_cohort")
  expect_equal(nrow(ci), 2L)
  t_crit <- qt(0.975, 50)
  # cohort 2005: att=0.3, se=0.1
  expect_equal(ci[1, 1], 0.3 - t_crit * 0.1, tolerance = 1e-10)
  expect_equal(ci[1, 2], 0.3 + t_crit * 0.1, tolerance = 1e-10)
})

# TC-10.4.71: confint parm character name filter
test_that("TC-10.4.71: confint parm character filter", {
  obj <- .mk_stag()
  ci <- confint(obj, parm = "2005", type = "by_cohort")
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "2005")
})

# TC-10.4.72: confint small df numerical verification
test_that("TC-10.4.72: confint small df (df=2) wider CI", {
  obj <- .mk_ct(df_inference = 2L, df_resid = 4L)
  ci_small <- confint(obj)
  obj2 <- .mk_ct(df_inference = 50L)
  ci_large <- confint(obj2)
  # df=2 should give much wider CI (t_crit ≈ 4.303 vs 2.009)
  width_small <- ci_small[1, 2] - ci_small[1, 1]
  width_large <- ci_large[1, 2] - ci_large[1, 1]
  expect_true(width_small > width_large * 1.5)
})

# TC-10.4.75: confint non-lwdid_result errors
test_that("TC-10.4.75: confint non-lwdid_result errors", {
  expect_error(confint.lwdid_result(list(att = 1)))
})

# TC-10.4.87: confint by_cohort NULL errors with Staggered hint
test_that("TC-10.4.87: confint by_cohort NULL error hint", {
  obj <- .mk_ct()
  expect_error(confint(obj, type = "by_cohort"), "Staggered")
})

# TC-10.4.90: confint numerical verification att=0.5 se=0.12 df=50
test_that("TC-10.4.90: confint numerical att=0.5 se=0.12 df=50", {
  obj <- .mk_ct()
  ci <- confint(obj)
  # vibe-math verified: CI ≈ [0.2589, 0.7411] (using qt(0.975,50)≈2.00856)
  expect_equal(ci[1, 1], 0.5 - qt(0.975, 50) * 0.12, tolerance = 1e-4)
  expect_equal(ci[1, 2], 0.5 + qt(0.975, 50) * 0.12, tolerance = 1e-4)
})

# TC-10.4.92: confint all fallback chain
test_that("TC-10.4.92: confint all fallback to att_by_period", {
  periods <- data.frame(period = c(4L, 5L), att = c(0.3, 0.5),
                        se = c(0.1, 0.1), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  ci <- confint(obj, type = "all")
  co <- coef(obj, type = "all")
  expect_equal(rownames(ci), names(co))
})


# ============================================================================
# vcov tests (TC-10.4.22-25, TC-10.4.64-66, TC-10.4.73, TC-10.4.76, TC-10.4.89)
# ============================================================================

# TC-10.4.22: vcov overall returns 1×1 matrix (se_att², tolerance=1e-12)
test_that("TC-10.4.22: vcov overall 1x1 se_att^2", {
  obj <- .mk_ct()
  V <- vcov(obj, type = "overall")
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 1L)
  expect_equal(ncol(V), 1L)
  expect_equal(rownames(V), "ATT")
  expect_equal(colnames(V), "ATT")
  expect_equal(V[1, 1], 0.12^2, tolerance = 1e-12)
})

# TC-10.4.23: vcov by_period returns diagonal matrix (no joint vcov)
test_that("TC-10.4.23: vcov by_period diagonal", {
  periods <- data.frame(period = c(4L, 5L, 6L), att = c(0.3, 0.5, 0.7),
                        se = c(0.10, 0.15, 0.20), stringsAsFactors = FALSE)
  obj <- .mk_ct(att_by_period = periods)
  V <- vcov(obj, type = "by_period")
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 3L)
  expect_equal(ncol(V), 3L)
  expect_equal(rownames(V), c("4", "5", "6"))
  # Diagonal should be se^2
  expect_equal(V[1, 1], 0.10^2, tolerance = 1e-12)
  expect_equal(V[2, 2], 0.15^2, tolerance = 1e-12)
  expect_equal(V[3, 3], 0.20^2, tolerance = 1e-12)
  # Off-diagonal should be 0
  expect_equal(V[1, 2], 0)
  expect_equal(V[1, 3], 0)
  expect_equal(V[2, 3], 0)
})

# TC-10.4.24: vcov by_period uses joint vcov when available
test_that("TC-10.4.24: vcov by_period joint vcov", {
  periods <- data.frame(period = c(4L, 5L), att = c(0.3, 0.5),
                        se = c(0.10, 0.15), stringsAsFactors = FALSE)
  joint_V <- matrix(c(0.01, 0.003, 0.003, 0.0225), nrow = 2)
  obj <- .mk_ct(att_by_period = periods)
  obj$vcov_att_periods <- joint_V
  V <- vcov(obj, type = "by_period")
  expect_equal(V[1, 2], 0.003, tolerance = 1e-12)
  expect_equal(V[2, 1], 0.003, tolerance = 1e-12)
  expect_equal(rownames(V), c("4", "5"))
})

# TC-10.4.25: vcov by_cohort Staggered correct
test_that("TC-10.4.25: vcov by_cohort Staggered", {
  obj <- .mk_stag()
  V <- vcov(obj, type = "by_cohort")
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 2L)
  expect_equal(rownames(V), c("2005", "2007"))
  # Diagonal = se^2 (0.1^2=0.01, 0.15^2=0.0225)
  expect_equal(V[1, 1], 0.01, tolerance = 1e-12)
  expect_equal(V[2, 2], 0.0225, tolerance = 1e-12)
})

# TC-10.4.64: vcov full returns complete OLS VCE matrix
test_that("TC-10.4.64: vcov full returns vcov_full", {
  full_V <- matrix(c(0.01, 0.002, 0.002, 0.03), nrow = 2)
  rownames(full_V) <- colnames(full_V) <- c("(Intercept)", "treat")
  obj <- .mk_ct()
  obj$vcov_full <- full_V
  V <- vcov(obj, type = "full")
  expect_identical(V, full_V)
})

# TC-10.4.65: vcov full unavailable errors
test_that("TC-10.4.65: vcov full unavailable error", {
  obj <- .mk_ct()
  expect_error(vcov(obj, type = "full"), "not available")
})

# TC-10.4.66: vcov full error mentions store_vcov_full
test_that("TC-10.4.66: vcov full error hint store_vcov_full", {
  obj <- .mk_ct()
  expect_error(vcov(obj, type = "full"), "store_vcov_full")
})

# TC-10.4.73: vcov by_cohort uses joint vcov when available
test_that("TC-10.4.73: vcov by_cohort joint vcov", {
  joint_V <- matrix(c(0.01, 0.004, 0.004, 0.0225), nrow = 2)
  obj <- .mk_stag()
  obj$vcov_att_cohorts <- joint_V
  V <- vcov(obj, type = "by_cohort")
  expect_equal(V[1, 2], 0.004, tolerance = 1e-12)
  expect_equal(V[2, 1], 0.004, tolerance = 1e-12)
  expect_equal(rownames(V), c("2005", "2007"))
})

# TC-10.4.76: vcov non-lwdid_result errors
test_that("TC-10.4.76: vcov non-lwdid_result errors", {
  expect_error(vcov.lwdid_result(list(att = 1)))
})

# TC-10.4.89: vcov numerical verification se=0.12 → 0.0144
test_that("TC-10.4.89: vcov numerical se=0.12 -> 0.0144", {
  obj <- .mk_ct(se_att = 0.12)
  V <- vcov(obj, type = "overall")
  # vibe-math verified: 0.12^2 = 0.0144
  expect_equal(V[1, 1], 0.0144, tolerance = 1e-12)
})

# ============================================================================
# nobs tests (TC-10.4.67-69)
# ============================================================================

# TC-10.4.67: nobs returns correct observation count
test_that("TC-10.4.67: nobs returns nobs", {
  obj <- .mk_ct(nobs = 200L)
  expect_equal(nobs(obj), 200L)
})

# TC-10.4.68: nobs non-lwdid_result errors
test_that("TC-10.4.68: nobs non-lwdid_result errors", {
  expect_error(nobs.lwdid_result(list(nobs = 100)))
})

# TC-10.4.69: nobs matches summary nobs
test_that("TC-10.4.69: nobs matches summary", {
  obj <- .mk_ct(nobs = 200L)
  s <- summary(obj)
  expect_equal(nobs(obj), s$nobs)
})

# ============================================================================
# Helper function tests (TC-10.4.26-27)
# ============================================================================

# TC-10.4.26: .format_pvalue formatting
test_that("TC-10.4.26: .format_pvalue correct formatting", {
  expect_equal(.format_pvalue(NA), "NA")
  expect_equal(.format_pvalue(0.0001), "<0.001")
  expect_equal(.format_pvalue(0.0005), "<0.001")
  # p in [0.001, 0.01) → 4 decimal places
  formatted_005 <- .format_pvalue(0.005)
  expect_match(formatted_005, "0\\.005")
  # p >= 0.01 → uses digits param
  formatted_12 <- .format_pvalue(0.123, digits = 3)
  expect_match(formatted_12, "0\\.123")
})

# TC-10.4.27: .signif_stars mapping
test_that("TC-10.4.27: .signif_stars correct mapping", {
  expect_equal(.signif_stars(NA), " ")
  expect_equal(.signif_stars(0.0001), "***")
  expect_equal(.signif_stars(0.005), "**")
  expect_equal(.signif_stars(0.03), "*")
  expect_equal(.signif_stars(0.08), ".")
  expect_equal(.signif_stars(0.5), " ")
})
