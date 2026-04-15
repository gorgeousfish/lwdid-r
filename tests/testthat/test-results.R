# ============================================================================
# test-results.R — Comprehensive tests for the lwdid_result S3 class
# ============================================================================

# ── Helpers ─────────────────────────────────────────────────────────────────

# Build a minimal consistent CT object for tests
.make_ct_obj <- function(...) {
  defaults <- list(
    att = 2.5, se_att = 0.8, t_stat = 2.5 / 0.8,
    pvalue = 0.002, ci_lower = 0.9, ci_upper = 4.1,
    df_resid = 100L, df_inference = 98L,
    nobs = 500L, n_treated = 200L, n_control = 300L,
    K = 3L, tpost1 = 4L,
    depvar = "y", rolling = "level",
    vce_type = "ols", alpha = 0.05, is_staggered = FALSE
  )
  args <- modifyList(defaults, list(...))
  do.call(new_lwdid_result, args)
}

# Build a minimal consistent Staggered object for tests
.make_stag_obj <- function(...) {
  defaults <- list(
    att = 1.0, se_att = 0.5, t_stat = 2.0,
    pvalue = 0.05, ci_lower = 0.0, ci_upper = 2.0,
    df_resid = 200L, df_inference = 50L,
    nobs = 1000L, n_treated = 400L, n_control = 600L,
    K = 2L, tpost1 = 3L,
    depvar = "outcome", rolling = "level",
    vce_type = "cluster", cluster_var = "state", n_clusters = 51L,
    alpha = 0.05, is_staggered = TRUE,
    aggregate = "cohort",
    cohorts = c(2005L, 2007L),
    n_never_treated = 600L
  )
  args <- modifyList(defaults, list(...))
  do.call(new_lwdid_result, args)
}

# ============================================================================
# Group 1: Constructor tests
# ============================================================================

test_that("default constructor returns lwdid_result class", {
  obj <- new_lwdid_result()
  expect_s3_class(obj, "lwdid_result")
})

test_that("cmd is always 'lwdid'", {
  obj <- new_lwdid_result()
  expect_identical(obj$cmd, "lwdid")
})

test_that("method derived: FALSE -> common_timing, TRUE -> staggered", {
  obj_ct <- new_lwdid_result(is_staggered = FALSE)
  expect_identical(obj_ct$method, "common_timing")
  obj_st <- new_lwdid_result(is_staggered = TRUE)
  expect_identical(obj_st$method, "staggered")
})

test_that("default numeric fields are NA_real_", {
  obj <- new_lwdid_result()
  for (f in c("att", "se_att", "t_stat", "pvalue", "ci_lower", "ci_upper")) {
    expect_true(is.na(obj[[f]]) && is.double(obj[[f]]),
                info = sprintf("%s should be NA_real_", f))
  }
})

test_that("default integer fields are NA_integer_", {
  obj <- new_lwdid_result()
  for (f in c("df_resid", "df_inference", "nobs", "n_treated",
              "n_control", "K", "tpost1")) {
    expect_true(is.na(obj[[f]]), info = sprintf("%s should be NA", f))
  }
})

test_that("Staggered-specific attributes NULL when is_staggered=FALSE", {
  obj <- new_lwdid_result(is_staggered = FALSE)
  stag_attrs <- c("aggregate", "cohorts", "cohort_sizes", "n_never_treated",
                  "att_by_cohort", "att_by_cohort_time", "att_overall",
                  "se_overall", "ci_overall_lower", "ci_overall_upper",
                  "t_stat_overall", "pvalue_overall", "cohort_effects",
                  "event_time_effects", "cohort_weights")
  for (a in stag_attrs) {
    expect_null(obj[[a]], info = sprintf("%s should be NULL", a))
  }
})

test_that("warning_diagnostics and warnings_log default to empty list", {
  obj <- new_lwdid_result()
  expect_identical(obj$warning_diagnostics, list())
  expect_identical(obj$warnings_log, list())
})

test_that("params/bse/resid default numeric(0); controls default character(0)", {
  obj <- new_lwdid_result()
  expect_identical(obj$params, numeric(0))
  expect_identical(obj$bse, numeric(0))
  expect_identical(obj$resid, numeric(0))
  expect_identical(obj$controls, character(0))
})

test_that("all attributes accessible via $ accessor", {
  obj <- new_lwdid_result(att = 1.5, se_att = 0.3, is_staggered = FALSE)
  expect_equal(obj$att, 1.5)
  expect_equal(obj$se_att, 0.3)
  expect_false(obj$is_staggered)
  expect_equal(obj$alpha, 0.05)
})

test_that("full Common Timing construction", {
  att_val <- 2.5; se_val <- 0.8
  t_val <- att_val / se_val
  df_i <- 100L
  pval <- 2 * pt(-abs(t_val), df = df_i)
  tc <- qt(0.975, df_i)

  obj <- new_lwdid_result(
    att = att_val, se_att = se_val, t_stat = t_val,
    pvalue = pval, ci_lower = att_val - tc * se_val,
    ci_upper = att_val + tc * se_val,
    df_resid = 100L, df_inference = df_i,
    nobs = 500L, n_treated = 200L, n_control = 300L,
    K = 3L, tpost1 = 4L, depvar = "y", rolling = "level",
    vce_type = "robust", alpha = 0.05, is_staggered = FALSE,
    ivar = "id", tvar = "year"
  )
  expect_s3_class(obj, "lwdid_result")
  expect_identical(obj$method, "common_timing")
  expect_equal(obj$att, att_val)
  expect_null(obj$att_overall)
})

test_that("full Staggered construction", {
  obj <- new_lwdid_result(
    att = 1.0, se_att = 0.5, t_stat = 2.0,
    pvalue = 0.05, ci_lower = 0.0, ci_upper = 2.0,
    df_resid = 200L, df_inference = 50L,
    nobs = 1000L, n_treated = 400L, n_control = 600L,
    K = 2L, tpost1 = 3L, depvar = "outcome", rolling = "level",
    vce_type = "cluster", cluster_var = "state", n_clusters = 51L,
    alpha = 0.05, is_staggered = TRUE,
    aggregate = "cohort",
    cohorts = c(2005L, 2007L, 2010L),
    cohort_sizes = c("2005" = 100L, "2007" = 150L, "2010" = 150L),
    n_never_treated = 600L, att_overall = 1.2, se_overall = 0.4,
    ivar = "state_id", tvar = "year"
  )
  expect_s3_class(obj, "lwdid_result")
  expect_identical(obj$method, "staggered")
  expect_true(obj$is_staggered)
  expect_equal(obj$aggregate, "cohort")
  expect_equal(obj$att_overall, 1.2)
})

# ============================================================================
# Group 2: Validation tests (18 checks)
# ============================================================================

test_that("validate: default (NA) object passes", {
  expect_silent(validate_lwdid_result(new_lwdid_result()))
})

test_that("validate: non-lwdid_result fails (check 1)", {
  expect_error(validate_lwdid_result(list(att = 1)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: negative se_att fails (check 2)", {
  expect_error(validate_lwdid_result(new_lwdid_result(se_att = -0.5)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: t_stat mismatch fails (check 3)", {
  expect_error(
    validate_lwdid_result(new_lwdid_result(att = 2.0, se_att = 1.0, t_stat = 999)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: se_att=0 with t_stat=Inf passes (check 3 edge)", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 2.0, se_att = 0.0, t_stat = Inf)))
})

test_that("validate: se_att=0 with t_stat=NaN passes (check 3 edge)", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 0.0, se_att = 0.0, t_stat = NaN)))
})

test_that("validate: df_inference <= 0 fails (check 4)", {
  expect_error(validate_lwdid_result(new_lwdid_result(df_inference = 0L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: df_resid <= 0 fails (check 5)", {
  expect_error(validate_lwdid_result(new_lwdid_result(df_resid = -1L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: pvalue < 0 fails (check 6)", {
  expect_error(validate_lwdid_result(new_lwdid_result(pvalue = -0.01)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: pvalue > 1 fails (check 6)", {
  expect_error(validate_lwdid_result(new_lwdid_result(pvalue = 1.01)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: ci_lower > att fails (check 7)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(att = 1.0, ci_lower = 2.0, ci_upper = 3.0)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: att > ci_upper fails (check 7)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(att = 5.0, ci_lower = 1.0, ci_upper = 3.0)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: n_treated < 0 fails (check 8)", {
  expect_error(validate_lwdid_result(new_lwdid_result(n_treated = -1L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: n_control < 0 fails (check 9)", {
  expect_error(validate_lwdid_result(new_lwdid_result(n_control = -1L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: nobs < n_treated+n_control fails for CT (check 10)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(nobs = 10L, n_treated = 6L, n_control = 6L,
                     is_staggered = FALSE)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: nobs check skipped for Staggered (check 10)", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(nobs = 10L, n_treated = 6L, n_control = 6L,
                     is_staggered = TRUE, aggregate = "cohort")))
})

test_that("validate: cluster n_clusters < 2 fails (check 11)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(vce_type = "cluster", n_clusters = 1L)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: alpha <= 0 fails (check 12)", {
  expect_error(validate_lwdid_result(new_lwdid_result(alpha = 0)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: alpha >= 1 fails (check 12)", {
  expect_error(validate_lwdid_result(new_lwdid_result(alpha = 1)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: is_staggered=NA fails (check 13)", {
  obj <- new_lwdid_result()
  obj$is_staggered <- NA
  expect_error(validate_lwdid_result(obj), class = "lwdid_invalid_parameter")
})

test_that("validate: staggered but aggregate=NULL fails (check 14)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(is_staggered = TRUE, aggregate = NULL)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: CT but att_overall non-NULL fails (check 15)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(is_staggered = FALSE, att_overall = 1.5)),
    class = "lwdid_invalid_parameter")
})

test_that("validate: df_inference > df_resid non-cluster fails (check 16)", {
  expect_error(validate_lwdid_result(
    new_lwdid_result(df_inference = 100L, df_resid = 50L, vce_type = "robust")),
    class = "lwdid_invalid_parameter")
})

test_that("validate: df_inference > df_resid cluster warns only (check 16)", {
  expect_warning(validate_lwdid_result(
    new_lwdid_result(df_inference = 100L, df_resid = 50L,
                     vce_type = "cluster", cluster_var = "st", n_clusters = 51L)),
    class = "lwdid_data")
})

test_that("validate: K < 0 fails (check 17)", {
  expect_error(validate_lwdid_result(new_lwdid_result(K = -1L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: tpost1 <= K fails (check 18)", {
  expect_error(validate_lwdid_result(new_lwdid_result(tpost1 = 3L, K = 3L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: tpost1 < K fails (check 18)", {
  expect_error(validate_lwdid_result(new_lwdid_result(tpost1 = 2L, K = 3L)),
               class = "lwdid_invalid_parameter")
})

test_that("validate: consistent CT object passes all checks", {
  att_val <- 2.5; se_val <- 0.8; t_val <- att_val / se_val; df_i <- 100L
  pval <- 2 * pt(-abs(t_val), df = df_i)
  tc <- qt(0.975, df_i)
  obj <- new_lwdid_result(
    att = att_val, se_att = se_val, t_stat = t_val, pvalue = pval,
    ci_lower = att_val - tc * se_val, ci_upper = att_val + tc * se_val,
    df_resid = 100L, df_inference = df_i,
    nobs = 500L, n_treated = 200L, n_control = 300L,
    K = 3L, tpost1 = 4L, alpha = 0.05, is_staggered = FALSE)
  expect_silent(validate_lwdid_result(obj))
})

test_that("validate: consistent Staggered object passes all checks", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(
      att = 1.0, se_att = 0.5, t_stat = 2.0,
      pvalue = 0.05, ci_lower = 0.0, ci_upper = 2.0,
      df_resid = 200L, df_inference = 50L,
      nobs = 1000L, n_treated = 400L, n_control = 600L,
      K = 2L, tpost1 = 3L, alpha = 0.05,
      is_staggered = TRUE, aggregate = "cohort",
      vce_type = "cluster", cluster_var = "st", n_clusters = 51L)))
})

# ============================================================================
# Group 3: print format tests
# ============================================================================

test_that("print CT: contains title and key values", {
  out <- capture.output(print(.make_ct_obj()))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Lee-Wooldridge DiD Estimation", combined))
  expect_true(grepl("2.5000", combined))   # ATT
  expect_true(grepl("0.8000", combined))   # SE
  expect_true(grepl("K = 3", combined))
  expect_true(grepl("tpost1 = 4", combined))
})

test_that("print CT: att_by_period shows hint text", {
  abp <- data.frame(period = 4:6, att = c(1, 1.5, 2), se = c(.3, .4, .5))
  out <- capture.output(print(.make_ct_obj(att_by_period = abp)))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Period-by-period", combined))
  expect_true(grepl("results\\$att_by_period", combined))
})

test_that("print CT: att_by_period > 5 rows shows 'more periods'", {
  abp <- data.frame(period = 4:10, att = seq(1, 7), se = rep(0.3, 7))
  out <- capture.output(print(.make_ct_obj(att_by_period = abp)))
  expect_true(any(grepl("more periods", out)))
})

test_that("print CT: RI shows 'Randomization Inference'", {
  obj <- .make_ct_obj(ri_pvalue = 0.03, ri_method = "fisher",
                      ri_seed = 42L, ri_valid = 999L, rireps = 1000L)
  out <- capture.output(print(obj))
  expect_true(any(grepl("Randomization Inference", out)))
})

test_that("print Staggered: contains 'Staggered' and 'N Never-treated'", {
  out <- capture.output(print(.make_stag_obj()))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Staggered", combined))
  expect_true(grepl("N Never-treated", combined))
})

test_that("print Staggered: aggregate='none' shows hint, no (g,r) data", {
  out <- capture.output(print(.make_stag_obj(aggregate = "none")))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("results\\$att_by_cohort_time", combined))
  expect_true(grepl("plot_event_study", combined))
})

test_that("print Staggered: control_group auto-switch shown", {
  obj <- .make_stag_obj(control_group = "not_yet_treated",
                        control_group_used = "never_treated")
  out <- capture.output(print(obj))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("never_treated", combined))
  expect_true(grepl("not_yet_treated", combined))
})

test_that("print: include_pretreatment shows 'Parallel Trends Test'", {
  pt <- list(f_stat = 1.5, pvalue = 0.22, reject_null = FALSE, alpha = 0.05)
  obj <- .make_ct_obj(include_pretreatment = TRUE, parallel_trends_test = pt)
  out <- capture.output(print(obj))
  expect_true(any(grepl("Parallel Trends Test", out)))
})

test_that("print: VCE ols -> 'OLS (Homoskedastic)'", {
  out <- capture.output(print(.make_ct_obj(vce_type = "ols")))
  expect_true(any(grepl("OLS \\(Homoskedastic\\)", out)))
})

test_that("print: VCE cluster -> 'Cluster-robust' with var name", {
  obj <- .make_ct_obj(vce_type = "cluster", cluster_var = "firm", n_clusters = 30L)
  out <- capture.output(print(obj))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Cluster-robust", combined))
  expect_true(grepl("firm", combined))
})

# ============================================================================
# Group 4: coef tests
# ============================================================================

test_that("coef: no params returns c(ATT = att_value)", {
  result <- coef(new_lwdid_result(att = 3.14))
  expect_equal(result, c(ATT = 3.14))
  expect_equal(names(result), "ATT")
})

test_that("constructor stores full params vector unchanged", {
  p <- c(intercept = 1.0, treat = 2.0, post = 0.5, did = 3.0)
  obj <- new_lwdid_result(att = 3.0, params = p)
  expect_identical(obj$params, p)
  expect_equal(coef(obj), c(ATT = 3.0))
})

# ============================================================================
# Group 5: confint tests (numerical correctness)
# ============================================================================

test_that("confint: default level=0.95 matches stored CI", {
  att <- 2.0; se <- 0.5; df_i <- 30L
  tc <- qt(0.975, df_i)
  ci_lo <- att - tc * se; ci_hi <- att + tc * se
  obj <- new_lwdid_result(
    att = att, se_att = se, t_stat = att / se,
    pvalue = 2 * pt(-abs(att / se), df_i),
    ci_lower = ci_lo, ci_upper = ci_hi,
    df_resid = 100L, df_inference = df_i,
    nobs = 200L, n_treated = 80L, n_control = 100L,
    K = 3L, tpost1 = 4L, is_staggered = FALSE)
  ci <- confint(obj, level = 0.95)
  expect_equal(ci[1, 1], ci_lo, tolerance = 1e-10)
  expect_equal(ci[1, 2], ci_hi, tolerance = 1e-10)
})

test_that("confint: level=0.90 narrower, level=0.99 wider than 0.95", {
  obj <- new_lwdid_result(att = 1.0, se_att = 0.5, df_inference = 30L,
                          df_resid = 100L, is_staggered = FALSE)
  ci_90 <- confint(obj, level = 0.90)
  ci_95 <- confint(obj, level = 0.95)
  ci_99 <- confint(obj, level = 0.99)
  # 90% narrower than 95%
  expect_true(ci_90[1, 1] > ci_95[1, 1])
  expect_true(ci_90[1, 2] < ci_95[1, 2])
  # 99% wider than 95%
  expect_true(ci_99[1, 1] < ci_95[1, 1])
  expect_true(ci_99[1, 2] > ci_95[1, 2])
})

test_that("confint: att=1.0, se=0.5, df=30 numerical correctness", {
  # Pre-computed: qt(0.975, 30) = 2.042272456301238
  # ci_lo = 1.0 - 2.042272456301238 * 0.5 = -0.021136228150619
  # ci_hi = 1.0 + 2.042272456301238 * 0.5 =  2.021136228150619
  obj <- new_lwdid_result(att = 1.0, se_att = 0.5, df_inference = 30L,
                          df_resid = 100L, is_staggered = FALSE)
  ci <- confint(obj, level = 0.95)
  tc <- qt(0.975, 30)
  expect_equal(ci[1, 1], 1.0 - tc * 0.5, tolerance = 1e-10)
  expect_equal(ci[1, 2], 1.0 + tc * 0.5, tolerance = 1e-10)
  # Numerical reasonableness: lower bound slightly negative, upper > 1.5
  expect_true(ci[1, 1] < 0)
  expect_true(ci[1, 2] > 1.5)
  expect_true(ci[1, 2] - ci[1, 1] > 2.0)  # width > 2*se
})

test_that("confint: return matrix dimnames correct", {
  obj <- new_lwdid_result(att = 1.0, se_att = 0.5, df_inference = 30L,
                          df_resid = 100L, is_staggered = FALSE)
  ci95 <- confint(obj, level = 0.95)
  expect_true(is.matrix(ci95))
  expect_equal(dim(ci95), c(1L, 2L))
  expect_equal(rownames(ci95), "ATT")
  expect_equal(colnames(ci95), c("2.5 %", "97.5 %"))

  ci90 <- confint(obj, level = 0.90)
  expect_equal(colnames(ci90), c("5.0 %", "95.0 %"))
})

test_that("confint: NA att/se/df returns NA matrix", {
  ci <- confint(new_lwdid_result(), level = 0.95)
  expect_true(is.matrix(ci))
  expect_true(is.na(ci[1, 1]))
  expect_true(is.na(ci[1, 2]))
})

test_that("confint: width ordering 90% < 95% < 99%", {
  obj <- new_lwdid_result(att = 1.0, se_att = 0.5, df_inference = 30L,
                          df_resid = 100L, is_staggered = FALSE)
  w90 <- diff(confint(obj, level = 0.90)[1, ])
  w95 <- diff(confint(obj, level = 0.95)[1, ])
  w99 <- diff(confint(obj, level = 0.99)[1, ])
  expect_true(w90 < w95)
  expect_true(w95 < w99)
  # All widths positive and reasonable for se=0.5
  expect_true(w90 > 1.0)
  expect_true(w99 < 4.0)
})

# ============================================================================
# Group 6: vcov tests
# ============================================================================

test_that("constructor stores vcov_matrix attribute", {
  m <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
  expect_identical(new_lwdid_result(vcov_matrix = m)$vcov_matrix, m)
})

test_that("constructor leaves vcov_matrix NULL when not set", {
  expect_null(new_lwdid_result()$vcov_matrix)
})

test_that("stored vcov_matrix dimensions can match params length", {
  p <- c(b0 = 1.0, b1 = 2.0, b2 = 3.0)
  m <- matrix(runif(9), nrow = 3, ncol = 3)
  obj <- new_lwdid_result(params = p, vcov_matrix = m)
  expect_equal(nrow(obj$vcov_matrix), length(p))
  expect_equal(ncol(obj$vcov_matrix), length(p))
})

# ============================================================================
# Group 7: summary tests
# ============================================================================

test_that("summary CT: returns object invisibly, prints output", {
  result <- summary(.make_ct_obj())
  expect_s3_class(result, "summary.lwdid_result")
  out <- capture.output(print(result))
  expect_true(any(grepl("Local Wald DID Estimation Summary", out)))
})

test_that("summary CT: shows coefficient table when params present", {
  p <- c(intercept = 1.0, treat = 2.0, post = 0.5, did = 3.0)
  bse <- c(intercept = 0.1, treat = 0.3, post = 0.2, did = 0.8)
  out <- capture.output(summary(.make_ct_obj(params = p, bse = bse)))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Regression Coefficients", combined))
  expect_true(grepl("intercept", combined))
  expect_true(grepl("did", combined))
})

test_that("summary CT: att_by_period shows ALL rows (not truncated)", {
  abp <- data.frame(period = 4:12, att = seq(1, 9), se = rep(0.3, 9))
  out <- capture.output(summary(.make_ct_obj(att_by_period = abp)))
  combined <- paste(out, collapse = "\n")
  # summary should NOT contain "more periods" (unlike print)
  expect_false(grepl("more periods", combined))
})

test_that("summary Staggered: contains cohort info", {
  obj <- .make_stag_obj(
    att_overall = 1.5, se_overall = 0.3,
    t_stat_overall = 5.0, pvalue_overall = 0.001,
    ci_overall_lower = 0.9, ci_overall_upper = 2.1,
    cohort_effects = list(
      "2005" = list(att = 1.2, se = 0.4),
      "2007" = list(att = 1.8, se = 0.5)))
  out <- capture.output(summary(obj))
  combined <- paste(out, collapse = "\n")
  expect_true(grepl("Staggered", combined))
  expect_true(grepl("Cohort Effects", combined))
})

# ============================================================================
# Group 8: Plot surface tests
# ============================================================================

test_that("plot.lwdid_result: empty result errors without event-study data", {
  expect_error(plot(new_lwdid_result()), "att_by_period")
})

test_that("plot.lwdid_result: signature exposes type dispatch", {
  args <- names(formals(plot.lwdid_result))
  expect_true("type" %in% args)
  expect_true("..." %in% args)
})

test_that("plot_event_study generic is defined", {
  expect_true(is.function(plot_event_study))
})

test_that("plot_event_study.lwdid_result: empty result errors without period data", {
  expect_error(plot_event_study(new_lwdid_result()), "att_by_period")
})

test_that("plot_event_study.lwdid_result: signature matches current surface", {
  args <- names(formals(plot_event_study.lwdid_result))
  for (p in c("ci_level", "reference_line", "show_pre_treatment",
              "facet_by_cohort", "color_significant", "ref_period",
              "aggregation", "df_strategy", "return_data", "point_size",
              "errorbar_width", "title", "theme")) {
    expect_true(p %in% args, info = sprintf("Missing param: %s", p))
  }
})

# ============================================================================
# Group 9: Attribute completeness tests
# ============================================================================

test_that("CT object: all staggered attrs NULL", {
  obj <- .make_ct_obj()
  for (a in c("aggregate", "cohorts", "cohort_sizes", "n_never_treated",
              "att_by_cohort", "att_by_cohort_time", "att_overall",
              "se_overall", "ci_overall_lower", "ci_overall_upper",
              "t_stat_overall", "pvalue_overall", "cohort_effects",
              "event_time_effects", "cohort_weights")) {
    expect_null(obj[[a]], info = sprintf("Expected NULL: %s", a))
  }
})

test_that("Staggered object: aggregate and cohorts non-NULL", {
  obj <- .make_stag_obj()
  expect_false(is.null(obj$aggregate))
  expect_false(is.null(obj$cohorts))
})

test_that("data attribute: set and accessible", {
  dt <- data.frame(x = 1:5, y = 6:10)
  obj <- new_lwdid_result(data = dt)
  expect_equal(nrow(obj$data), 5L)
})

test_that("ri_error/ri_target: default NULL, settable via constructor", {
  expect_null(new_lwdid_result()$ri_error)
  expect_null(new_lwdid_result()$ri_target)
  obj <- new_lwdid_result(ri_error = "err", ri_target = "att")
  expect_equal(obj$ri_error, "err")
  expect_equal(obj$ri_target, "att")
})

test_that("ivar, tvar, is_quarterly correctly stored", {
  obj <- new_lwdid_result(ivar = "uid", tvar = "year", is_quarterly = FALSE)
  expect_equal(obj$ivar, "uid")
  expect_equal(obj$tvar, "year")
  expect_false(obj$is_quarterly)
})

test_that("tvar length 1 = annual; length 2 = quarterly", {
  obj_a <- new_lwdid_result(tvar = "year", is_quarterly = FALSE)
  expect_equal(length(obj_a$tvar), 1L)
  obj_q <- new_lwdid_result(tvar = c("year", "qtr"), is_quarterly = TRUE)
  expect_equal(length(obj_q$tvar), 2L)
  expect_true(obj_q$is_quarterly)
})

# ============================================================================
# Group 10: Edge case / boundary tests
# ============================================================================

test_that("se_att=0, att=5 -> t_stat=Inf passes validation", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 5.0, se_att = 0, t_stat = Inf)))
})

test_that("se_att=0, att=0 -> t_stat=NaN passes validation", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 0, se_att = 0, t_stat = NaN)))
})

test_that("single-cohort Staggered object valid", {
  obj <- new_lwdid_result(
    att = 1.5, se_att = 0.3, t_stat = 5.0,
    pvalue = 0.001, ci_lower = 0.9, ci_upper = 2.1,
    df_resid = 200L, df_inference = 50L,
    nobs = 300L, n_treated = 100L, n_control = 200L,
    K = 3L, tpost1 = 4L,
    is_staggered = TRUE, aggregate = "cohort", cohorts = 2005L)
  expect_equal(length(obj$cohorts), 1L)
  expect_silent(validate_lwdid_result(obj))
})

test_that("extreme att values: 1e15 and 1e-15", {
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 1e15, se_att = 1e14, t_stat = 10.0)))
  expect_silent(validate_lwdid_result(
    new_lwdid_result(att = 1e-15, se_att = 1e-16, t_stat = 10.0)))
})

# ============================================================================
# Group 11: End-to-end numerical reasonableness
# ============================================================================

test_that("e2e: att=2.0, se=0.5, df=30 full numerical check", {
  att <- 2.0; se <- 0.5; df_inf <- 30L
  expected_t <- att / se  # 4.0
  expected_pval <- 2 * pt(-abs(expected_t), df = df_inf)
  crit <- qt(0.975, df_inf)
  expected_ci_lo <- att - crit * se
  expected_ci_hi <- att + crit * se

  # Verify intermediate values

  expect_equal(expected_t, 4.0, tolerance = 1e-12)
  expect_true(expected_pval > 0.0003 && expected_pval < 0.0004)

  obj <- new_lwdid_result(
    att = att, se_att = se, t_stat = expected_t, pvalue = expected_pval,
    ci_lower = expected_ci_lo, ci_upper = expected_ci_hi,
    df_resid = 100L, df_inference = df_inf,
    nobs = 200L, n_treated = 80L, n_control = 100L,
    K = 4L, tpost1 = 5L, is_staggered = FALSE)
  expect_silent(validate_lwdid_result(obj))

  ci <- confint(obj, level = 0.95)
  expect_equal(ci[1, 1], expected_ci_lo, tolerance = 1e-10)
  expect_equal(ci[1, 2], expected_ci_hi, tolerance = 1e-10)
  # CI bounds reasonable for t=4
  expect_true(ci[1, 1] > 0.5)
  expect_true(ci[1, 2] < 4.0)
})

test_that("e2e: confint width ordering 90% < 95% < 99%", {
  obj <- new_lwdid_result(att = 1.0, se_att = 0.5, df_inference = 30L,
                          df_resid = 100L, is_staggered = FALSE)
  w90 <- diff(confint(obj, level = 0.90)[1, ])
  w95 <- diff(confint(obj, level = 0.95)[1, ])
  w99 <- diff(confint(obj, level = 0.99)[1, ])
  expect_true(w90 < w95 && w95 < w99)
  # Reasonable range for se=0.5, df=30
  expect_true(w90 > 1.0 && w99 < 4.0)
})


# ── E5-05.7: new_lwdid_result() and extract_effects() tests ─────────────────

test_that("new_lwdid_result: returns lwdid_result class", {
  result <- new_lwdid_result(
    att_by_cohort_time = data.frame(
      cohort = c(3L, 5L), ref_period = c(2L, 4L),
      att = c(0.1, 0.2), se = c(0.01, 0.02),
      n_treated = c(3L, 2L)
    ),
    is_staggered = TRUE, aggregate = "none"
  )
  expect_s3_class(result, "lwdid_result")
  expect_true(result$is_staggered)
})

test_that("new_lwdid_result: top-level inference computed correctly", {
  result <- new_lwdid_result(
    att = 0.1, se_att = 0.02, df_inference = 50L,
    is_staggered = TRUE, aggregate = "cohort"
  )
  # t_stat = att / se_att = 0.1 / 0.02 = 5.0
  expect_equal(result$t_stat, 0.1 / 0.02, tolerance = 1e-10)
  expect_true(!is.na(result$pvalue))
  expect_true(result$pvalue > 0 && result$pvalue < 1)
  # pvalue = 2 * pt(5.0, df=50, lower.tail=FALSE)
  expected_p <- 2 * pt(abs(5.0), df = 50, lower.tail = FALSE)
  expect_equal(result$pvalue, expected_p, tolerance = 1e-10)
  # CI should contain att
  expect_true(!is.na(result$ci_lower))
  expect_true(!is.na(result$ci_upper))
  expect_true(result$ci_lower < result$att)
  expect_true(result$ci_upper > result$att)
  # Verify CI formula: att ± qt(0.975, 50) * se
  t_crit <- qt(0.975, df = 50)
  expect_equal(result$ci_lower, 0.1 - t_crit * 0.02, tolerance = 1e-10)
  expect_equal(result$ci_upper, 0.1 + t_crit * 0.02, tolerance = 1e-10)
})

test_that("new_lwdid_result: NA att gives NA inference", {
  result <- new_lwdid_result(
    att = NA_real_, se_att = NA_real_,
    is_staggered = TRUE, aggregate = "none"
  )
  expect_true(is.na(result$t_stat))
  expect_true(is.na(result$pvalue))
  expect_true(is.na(result$ci_lower))
  expect_true(is.na(result$ci_upper))
})

test_that("new_lwdid_result: cohort_agg fields stored correctly", {
  result <- new_lwdid_result(
    att_cohort_agg = 0.15, se_cohort_agg = 0.03,
    t_stat_cohort_agg = 5.0, pvalue_cohort_agg = 0.001,
    ci_cohort_agg = c(0.09, 0.21),
    is_staggered = TRUE, aggregate = "cohort"
  )
  expect_equal(result$att_cohort_agg, 0.15)
  expect_equal(result$se_cohort_agg, 0.03)
  expect_equal(result$t_stat_cohort_agg, 5.0)
  expect_equal(result$pvalue_cohort_agg, 0.001)
  # ci_cohort_agg split into lower/upper
  expect_equal(result$ci_lower_cohort_agg, 0.09)
  expect_equal(result$ci_upper_cohort_agg, 0.21)
})

test_that("new_lwdid_result: control_group_used defaults to control_group", {
  result <- new_lwdid_result(
    control_group = "not_yet_treated",
    is_staggered = TRUE, aggregate = "none"
  )
  expect_equal(result$control_group_used, "not_yet_treated")
})

test_that("extract_effects: non-lwdid_result raises error", {
  expect_error(
    extract_effects(list(att = 0.1)),
    class = "lwdid_invalid_input"
  )
})

test_that("extract_effects: type=NULL auto-infers from aggregate=none as gr", {
  gr_effects <- data.frame(att = c(0.1, 0.2), se = c(0.01, 0.02))
  result <- new_lwdid_result(
    att_by_cohort_time = gr_effects,
    is_staggered = TRUE, aggregate = "none"
  )
  df <- extract_effects(result)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2L)
})

test_that("extract_effects: type=gr returns cohort_time_effects", {
  gr_effects <- data.frame(att = c(0.1, 0.2), se = c(0.01, 0.02))
  result <- new_lwdid_result(
    att_by_cohort_time = gr_effects,
    is_staggered = TRUE, aggregate = "none"
  )
  df <- extract_effects(result, type = "gr")
  expect_equal(nrow(df), 2L)
  expect_true("att" %in% names(df))
})

test_that("extract_effects: type=cohort converts list to data.frame", {
  cohort_effects <- list(
    list(cohort = 3L, att = 0.1, se = 0.01, ci_lower = 0.08,
         ci_upper = 0.12, t_stat = 10.0, pvalue = 0.001,
         n_periods = 4L, n_units = 3L, n_control = 2L,
         df_resid = 10L, df_inference = 8L),
    list(cohort = 5L, att = 0.2, se = 0.02, ci_lower = 0.16,
         ci_upper = 0.24, t_stat = 10.0, pvalue = 0.001,
         n_periods = 2L, n_units = 2L, n_control = 2L,
         df_resid = 8L, df_inference = 6L)
  )
  result <- new_lwdid_result(
    cohort_effects = cohort_effects,
    is_staggered = TRUE, aggregate = "cohort"
  )
  df <- extract_effects(result, type = "cohort")
  expect_equal(nrow(df), 2L)
  expect_true("cohort" %in% names(df))
  expect_true("n_control" %in% names(df))
  expect_equal(df$cohort, c(3L, 5L))
  expect_equal(df$att, c(0.1, 0.2))
  # n_control defensively extracted
  expect_equal(df$n_control, c(2L, 2L))
})

test_that("extract_effects: type=overall returns single-row data.frame", {
  overall_effect <- list(
    att = 0.15, se = 0.02, ci_lower = 0.11, ci_upper = 0.19,
    t_stat = 7.5, pvalue = 0.0001,
    n_treated = 5L, n_control = 2L,
    df_resid = 10L, df_inference = 8L
  )
  result <- new_lwdid_result(
    is_staggered = TRUE, aggregate = "overall"
  )
  # Store overall_effect on the result object
  result$overall_effect <- overall_effect
  df <- extract_effects(result, type = "overall")
  expect_equal(nrow(df), 1L)
  expect_equal(df$att, 0.15)
  expect_equal(df$se, 0.02)
  expect_equal(df$n_treated, 5L)
})

test_that("extract_effects: type=event_time converts list to data.frame", {
  et_effects <- list(
    list(event_time = 0L, att = 0.1, se = 0.01, ci_lower = 0.08,
         ci_upper = 0.12, t_stat = 10.0, pvalue = 0.001,
         df_inference = 8L, n_cohorts = 2L, weight_sum = 5.0),
    list(event_time = 1L, att = 0.2, se = 0.02, ci_lower = 0.16,
         ci_upper = 0.24, t_stat = 10.0, pvalue = 0.001,
         df_inference = 8L, n_cohorts = 2L, weight_sum = 5.0)
  )
  result <- new_lwdid_result(
    event_time_effects = et_effects,
    is_staggered = TRUE, aggregate = "event_time"
  )
  df <- extract_effects(result, type = "event_time")
  expect_equal(nrow(df), 2L)
  expect_true("event_time" %in% names(df))
  expect_true("weight_sum" %in% names(df))
  expect_equal(df$event_time, c(0L, 1L))
  expect_equal(df$att, c(0.1, 0.2))
})

test_that("extract_effects: missing effects warns and returns empty df", {
  result <- new_lwdid_result(
    is_staggered = TRUE, aggregate = "none"
  )
  expect_warning(
    df <- extract_effects(result, type = "cohort"),
    class = "lwdid_data"
  )
  expect_equal(nrow(df), 0L)
})

test_that("extract_effects: invalid type raises error", {
  result <- new_lwdid_result(
    is_staggered = TRUE, aggregate = "none"
  )
  expect_error(
    extract_effects(result, type = "invalid"),
    class = "lwdid_invalid_input"
  )
})


# ── E5-05.8: print.lwdid_result() tests ─────────────────────────────────────

test_that("print.lwdid_result: staggered displays header and metadata", {
  result <- new_lwdid_result(
    att = 0.15, se_att = 0.02,
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated",
    control_group_used = "not_yet_treated",
    aggregate = "none",
    n_units = 7L, n_periods = 6L, n_cohorts = 2L, n_never_treated = 2L,
    df_inference = 50L, alpha = 0.05,
    is_staggered = TRUE, depvar = "y",
    method = "staggered", vce_type = "ols"
  )
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Local Wald DID Estimation", combined))
  expect_true(grepl("Staggered DID", combined))
  expect_true(grepl("Aggregate: none", combined))
  expect_true(grepl("Units: 7", combined))
  expect_true(grepl("Periods: 6", combined))
  expect_true(grepl("Cohorts: 2", combined))
  expect_true(grepl("NT: 2", combined))
})

test_that("print.lwdid_result: displays top-level ATT summary", {
  result <- new_lwdid_result(
    att = 0.092, se_att = 0.025,
    df_inference = 50L, alpha = 0.05,
    is_staggered = TRUE, aggregate = "none",
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("ATT", combined))
  expect_true(grepl("0.0920", combined))
  expect_true(grepl("SE", combined))
  expect_true(grepl("0.0250", combined))
  # Should have CI since att, se, df are all valid
  expect_true(grepl("CI", combined))
})

test_that("print.lwdid_result: returns invisible(x)", {
  result <- new_lwdid_result(
    is_staggered = TRUE, aggregate = "none",
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  ret <- withVisible(print(result))
  expect_false(ret$visible)
  expect_identical(ret$value, result)
})

test_that("print.lwdid_result: shows control group auto-switch info", {
  result <- new_lwdid_result(
    control_group = "not_yet_treated",
    control_group_used = "never_treated",
    is_staggered = TRUE, aggregate = "cohort",
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("never_treated", combined))
  expect_true(grepl("auto-switched from not_yet_treated", combined))
})

test_that("print.lwdid_result: cohort aggregate shows cohort table", {
  cohort_effects <- list(
    list(cohort = 3L, att = 0.085, se = 0.030,
         ci_lower = 0.02, ci_upper = 0.15,
         t_stat = 2.83, pvalue = 0.005,
         n_periods = 4L, n_units = 3L,
         df_resid = 10L, df_inference = 8L),
    list(cohort = 5L, att = 0.095, se = 0.028,
         ci_lower = 0.04, ci_upper = 0.15,
         t_stat = 3.39, pvalue = 0.0008,
         n_periods = 2L, n_units = 2L,
         df_resid = 8L, df_inference = 6L)
  )
  result <- new_lwdid_result(
    cohort_effects = cohort_effects,
    aggregate = "cohort",
    is_staggered = TRUE,
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Cohort Effects", combined))
})

test_that("print.lwdid_result: overall aggregate shows overall + cohort hint", {
  overall_effect <- list(
    att = 0.15, se = 0.02, ci_lower = 0.11, ci_upper = 0.19,
    t_stat = 7.5, pvalue = 0.0001
  )
  cohort_effects <- list(
    list(cohort = 3L, att = 0.1, se = 0.01, ci_lower = 0.08,
         ci_upper = 0.12, t_stat = 10.0, pvalue = 0.001,
         n_periods = 4L, n_units = 3L,
         df_resid = 10L, df_inference = 8L)
  )
  result <- new_lwdid_result(
    cohort_effects = cohort_effects,
    aggregate = "overall",
    is_staggered = TRUE,
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  result$overall_effect <- overall_effect
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Overall Effect", combined))
  expect_true(grepl("extract_effects", combined))
})

test_that("print.lwdid_result: event_time aggregate shows event time table", {
  et_effects <- list(
    list(event_time = 0L, att = 0.1, se = 0.01,
         ci_lower = 0.08, ci_upper = 0.12,
         t_stat = 10.0, pvalue = 0.001,
         df_inference = 8L, n_cohorts = 2L, weight_sum = 5.0)
  )
  result <- new_lwdid_result(
    event_time_effects = et_effects,
    aggregate = "event_time",
    is_staggered = TRUE,
    method = "staggered", estimator = "ra",
    rolling = "demean", depvar = "y", vce_type = "ols"
  )
  output <- capture.output(print(result))
  combined <- paste(output, collapse = "\n")
  expect_true(grepl("Event-Time Effects", combined))
})
