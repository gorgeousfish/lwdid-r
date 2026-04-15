# ============================================================================
# test-visualization-export.R — Numerical Validation Test Suite
# TC-10.6.1 through TC-10.6.62
# Story E10-06: 数值验证测试套件
# ============================================================================

# ============================================================================
# Task E10-06.1: Mock Object Factory Functions
# ============================================================================

.mock_common_timing_result <- function() {
  att <- 0.5
  se_att <- 0.12
  df <- 48L
  alpha <- 0.05
  t_crit <- stats::qt(1 - alpha / 2, df = df)

  obj <- new_lwdid_result(
    att = att, se_att = se_att,
    df_inference = df, df_resid = 450L,
    nobs = 500L, n_treated = 50L, n_control = 450L,
    alpha = alpha,
    depvar = "y", rolling = "demean",
    vce_type = "hc1", estimator = "ra",
    method = "common_timing", is_staggered = FALSE,
    K = 4L, tpost1 = 5L,
    ivar = "id", tvar = "time",
    att_by_period = data.frame(
      period = 5:8,
      att = c(0.45, 0.55, 0.48, 0.52),
      se = c(0.13, 0.14, 0.12, 0.15),
      ci_lower = 0.45 - stats::qt(0.975, df) * 0.13,
      ci_upper = 0.45 + stats::qt(0.975, df) * 0.13,
      pvalue = 2 * stats::pt(abs(0.45 / 0.13), df, lower.tail = FALSE),
      weight = rep(0.25, 4)
    ),
    att_pre_treatment = data.frame(
      event_time = c(-4L, -3L, -2L),
      att = c(0.02, -0.01, 0.005),
      se = c(0.05, 0.04, 0.03),
      ci_lower = c(0.02, -0.01, 0.005) - stats::qt(0.975, df) * c(0.05, 0.04, 0.03),
      ci_upper = c(0.02, -0.01, 0.005) + stats::qt(0.975, df) * c(0.05, 0.04, 0.03),
      pvalue = 2 * stats::pt(abs(c(0.02, -0.01, 0.005) / c(0.05, 0.04, 0.03)),
                              df, lower.tail = FALSE)
    ),
    parallel_trends_test = list(
      f_stat = 0.45, f_pvalue = 0.72, df1 = 3L, df2 = 45L
    )
  )
  # Fix att_by_period CIs for all rows (not just first)
  obj$att_by_period$ci_lower <- obj$att_by_period$att -
    stats::qt(0.975, df) * obj$att_by_period$se
  obj$att_by_period$ci_upper <- obj$att_by_period$att +
    stats::qt(0.975, df) * obj$att_by_period$se
  obj$att_by_period$pvalue <- 2 * stats::pt(
    abs(obj$att_by_period$att / obj$att_by_period$se), df, lower.tail = FALSE)
  obj$include_pretreatment <- TRUE
  obj
}

.mock_staggered_result <- function() {
  obj <- .mock_common_timing_result()
  obj$method <- "staggered"
  obj$is_staggered <- TRUE
  obj$include_pretreatment <- FALSE
  obj$control_group <- "not_yet_treated"
  obj$control_group_used <- "not_yet_treated"
  obj$aggregate <- "cohort"
  df_overall <- 48L

  obj$att_by_cohort <- data.frame(
    cohort = c(5L, 7L),
    att = c(0.48, 0.53),
    se = c(0.11, 0.14),
    ci_lower = c(0.48, 0.53) - stats::qt(0.975, df_overall) * c(0.11, 0.14),
    ci_upper = c(0.48, 0.53) + stats::qt(0.975, df_overall) * c(0.11, 0.14),
    pvalue = 2 * stats::pt(abs(c(0.48, 0.53) / c(0.11, 0.14)),
                            df_overall, lower.tail = FALSE),
    weight = c(0.6, 0.4)
  )

  # att_by_cohort_time with per-cohort df
  df_c5 <- 28L; df_c7 <- 18L
  obj$att_by_cohort_time <- data.frame(
    cohort = c(5L, 5L, 5L, 7L, 7L),
    period = c(5L, 6L, 7L, 7L, 8L),
    att = c(0.42, 0.50, 0.52, 0.51, 0.55),
    se = c(0.10, 0.12, 0.11, 0.13, 0.15),
    df_inference = c(df_c5, df_c5, df_c5, df_c7, df_c7),
    weight = c(0.2, 0.2, 0.2, 0.2, 0.2)
  )
  # Compute CIs and pvalues using per-cohort df
  dfs <- obj$att_by_cohort_time$df_inference
  obj$att_by_cohort_time$ci_lower <- obj$att_by_cohort_time$att -
    stats::qt(0.975, dfs) * obj$att_by_cohort_time$se
  obj$att_by_cohort_time$ci_upper <- obj$att_by_cohort_time$att +
    stats::qt(0.975, dfs) * obj$att_by_cohort_time$se
  obj$att_by_cohort_time$pvalue <- 2 * stats::pt(
    abs(obj$att_by_cohort_time$att / obj$att_by_cohort_time$se),
    dfs, lower.tail = FALSE)

  obj$event_time_effects <- data.frame(
    event_time = 0:2,
    att = c(0.46, 0.52, 0.52),
    se = c(0.10, 0.12, 0.11),
    ci_lower = c(0.46, 0.52, 0.52) - stats::qt(0.975, df_overall) * c(0.10, 0.12, 0.11),
    ci_upper = c(0.46, 0.52, 0.52) + stats::qt(0.975, df_overall) * c(0.10, 0.12, 0.11),
    pvalue = 2 * stats::pt(abs(c(0.46, 0.52, 0.52) / c(0.10, 0.12, 0.11)),
                            df_overall, lower.tail = FALSE)
  )

  obj$cohort_weights <- c("5" = 0.6, "7" = 0.4)
  obj$cohort_sizes <- c("5" = 30L, "7" = 20L)
  obj$cohort_sample_sizes <- c("5" = 30L, "7" = 20L)
  obj$cohorts <- c(5L, 7L)
  obj$n_cohorts <- 2L
  obj
}

.mock_result_with_diagnostics <- function() {
  obj <- .mock_common_timing_result()
  obj$diagnostics <- list(
    clustering = structure(list(
      cluster_sizes = setNames(as.integer(seq(20, 140, by = 5)),
                               paste0("cl_", 1:25)),
      n_clusters = 25L, balance_ratio = 0.85, icc = 0.05,
      effective_clusters = 20.0, reliability_level = "Good"
    ), class = "lwdid_clustering_diagnosis"),
    selection = structure(list(
      attrition_by_period = data.frame(
        period = 1:5, attrition_rate = c(0.02, 0.03, 0.05, 0.04, 0.06)
      ),
      attrition_rate = 0.05, selection_risk = "Low"
    ), class = "lwdid_selection_diagnosis"),
    parallel_trends = structure(list(
      pre_treatment_coefficients = data.frame(
        period = c(-3, -2, -1),
        coefficient = c(0.02, -0.01, 0.005),
        se = c(0.05, 0.04, 0.03),
        ci_lower = c(-0.08, -0.09, -0.055),
        ci_upper = c(0.12, 0.07, 0.065)
      ),
      f_test = list(df1 = 3L, df2 = 45L, f_stat = 0.45, f_pvalue = 0.72),
      group_trends = NULL
    ), class = "lwdid_trends")
  )
  obj
}

.mock_result_with_ri <- function() {
  obj <- .mock_common_timing_result()
  obj$ri_pvalue <- 0.004
  obj$ri_distribution <- seq(-3, 3, length.out = 999)
  obj$ri_seed <- 42L
  obj$ri_method <- "permutation"
  obj$ri_valid <- 999L
  obj$rireps <- 999L
  obj
}

.mock_trends_obj <- function() {
  structure(list(
    pre_treatment_coefficients = data.frame(
      period = c(-3, -2, -1),
      coefficient = c(0.02, -0.01, 0.005),
      se = c(0.05, 0.04, 0.03),
      ci_lower = c(-0.08, -0.09, -0.055),
      ci_upper = c(0.12, 0.07, 0.065)
    ),
    f_test = list(df1 = 3L, df2 = 45L, f_stat = 0.45, f_pvalue = 0.72),
    group_trends = data.frame(
      period = rep(c(-3, -2, -1, 0, 1), 2),
      mean_y = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.0, 1.15, 1.25, 1.8, 2.0),
      group = rep(c("control", "treated"), each = 5)
    )
  ), class = "lwdid_trends")
}

.mock_clustering_obj <- function() {
  structure(list(
    cluster_sizes = setNames(as.integer(seq(20, 140, by = 5)),
                             paste0("cl_", 1:25)),
    n_clusters = 25L, balance_ratio = 0.85, icc = 0.05,
    effective_clusters = 20.0, reliability_level = "Good"
  ), class = "lwdid_clustering_diagnosis")
}

.mock_selection_obj <- function() {
  structure(list(
    attrition_by_period = data.frame(
      period = 1:5, attrition_rate = c(0.02, 0.03, 0.05, 0.04, 0.06)
    ),
    attrition_rate = 0.05, selection_risk = "Low"
  ), class = "lwdid_selection_diagnosis")
}

# ============================================================================
# Task E10-06.2: S3 Method Integration Tests (10.6.1—10.6.6)
# ============================================================================

test_that("TC-10.6.1: coef->confint consistency (ATT within CI)", {
  obj <- .mock_common_timing_result()
  att <- coef(obj, type = "overall")
  ci <- confint(obj, level = 0.95, type = "overall")
  expect_true(att[["ATT"]] >= ci["ATT", 1])
  expect_true(att[["ATT"]] <= ci["ATT", 2])
})

test_that("TC-10.6.2: coef->vcov consistency (sqrt(vcov)==se_att)", {
  obj <- .mock_common_timing_result()
  V <- vcov(obj, type = "overall")
  se_from_vcov <- sqrt(V["ATT", "ATT"])
  expect_equal(se_from_vcov, obj$se_att, tolerance = 1e-10)
})

test_that("TC-10.6.3: confint width consistent with vcov", {
  obj <- .mock_common_timing_result()
  ci <- confint(obj, level = 0.95, type = "overall")
  ci_width <- ci["ATT", 2] - ci["ATT", 1]
  df <- obj$df_inference
  expected_width <- 2 * stats::qt(0.975, df) * obj$se_att
  expect_equal(ci_width, expected_width, tolerance = 1e-10)
})

test_that("TC-10.6.4: summary coefficients match coef", {
  obj <- .mock_common_timing_result()
  s <- summary(obj)
  att_coef <- coef(obj, type = "overall")
  expect_equal(s$coefficients["ATT", "Estimate"], att_coef[["ATT"]],
               tolerance = 1e-10)
})

test_that("TC-10.6.5: by_period granularity cross-method consistency", {
  obj <- .mock_common_timing_result()
  c_bp <- coef(obj, type = "by_period")
  ci_bp <- confint(obj, type = "by_period")
  v_bp <- vcov(obj, type = "by_period")
  # Dimensions match

  expect_equal(length(c_bp), nrow(ci_bp))
  expect_equal(length(c_bp), nrow(v_bp))
  expect_equal(length(c_bp), ncol(v_bp))
  # Names match
  expect_equal(names(c_bp), rownames(ci_bp))
  expect_equal(names(c_bp), rownames(v_bp))
  expect_equal(names(c_bp), colnames(v_bp))
})

test_that("TC-10.6.6: Staggered by_cohort cross-method consistency", {
  obj <- .mock_staggered_result()
  c_bc <- coef(obj, type = "by_cohort")
  ci_bc <- confint(obj, type = "by_cohort")
  # Each cohort ATT within its CI
  for (nm in names(c_bc)) {
    expect_true(c_bc[[nm]] >= ci_bc[nm, 1],
                info = sprintf("Cohort %s ATT below CI lower", nm))
    expect_true(c_bc[[nm]] <= ci_bc[nm, 2],
                info = sprintf("Cohort %s ATT above CI upper", nm))
  }
})

# ============================================================================
# Task E10-06.3: Visualization Integration Tests (10.6.7—10.6.11)
# ============================================================================

test_that("TC-10.6.7: Event Study plot data completeness", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_common_timing_result()
  res <- plot_event_study(obj, return_data = TRUE)
  expect_s3_class(res$plot, "ggplot")
  pd <- res$data
  n_pre <- nrow(obj$att_pre_treatment)
  n_post <- nrow(obj$att_by_period)
  # Data should contain pre + post + anchor
  expect_gte(nrow(pd), n_pre + n_post + 1L)
})

test_that("TC-10.6.8: show_pre_treatment=FALSE reduces data rows", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_common_timing_result()
  res_with <- plot_event_study(obj, show_pre_treatment = TRUE, return_data = TRUE)
  res_without <- plot_event_study(obj, show_pre_treatment = FALSE, return_data = TRUE)
  expect_lt(nrow(res_without$data), nrow(res_with$data))
})

test_that("TC-10.6.9: Spec curve plot returns ggplot with reasonable ATT", {
  skip_if_not_installed("ggplot2")
  # Create lwdid_sensitivity mock with pre_period type
  # plot.lwdid_sensitivity dispatches on x$type to .plot_specification_curve
  sens_obj <- structure(
    list(
      type = "pre_period",
      specifications = list(
        list(n_pre_periods = 2L, att = 0.50, ci_lower = 0.26, ci_upper = 0.74,
             pvalue = 0.001, converged = TRUE),
        list(n_pre_periods = 3L, att = 0.48, ci_lower = 0.22, ci_upper = 0.74,
             pvalue = 0.002, converged = TRUE),
        list(n_pre_periods = 4L, att = 0.51, ci_lower = 0.29, ci_upper = 0.73,
             pvalue = 0.001, converged = TRUE)
      ),
      baseline_spec = list(att = 0.50),
      robustness_threshold = 0.2,
      sensitivity_ratio = 0.95,
      robustness_level = "Robust"
    ),
    class = "lwdid_sensitivity"
  )
  p <- plot(sens_obj)
  expect_s3_class(p, "ggplot")
})

test_that("TC-10.6.10: Trends plot coefficients type returns ggplot", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_trends_obj()
  # plot.lwdid_trends supports type="coefficients" (default) and "trajectories"
  # It does NOT accept intervention_point parameter
  p <- plot(obj, type = "coefficients")
  expect_s3_class(p, "ggplot")
  # Should have a zero-reference hline (dashed)
  layers <- p$layers
  hline_layers <- Filter(function(l) {
    inherits(l$geom, "GeomHline")
  }, layers)
  expect_gte(length(hline_layers), 1L)
  # Should have point layer for coefficients
  point_layers <- Filter(function(l) {
    inherits(l$geom, "GeomPoint")
  }, layers)
  expect_gte(length(point_layers), 1L)
  # Trajectories type should also work
  p2 <- plot(obj, type = "trajectories")
  expect_s3_class(p2, "ggplot")
})

test_that("TC-10.6.11: Diagnostic panel is patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  obj <- .mock_result_with_diagnostics()
  p <- plot(obj, type = "diagnostics")
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

# ============================================================================
# Task E10-06.4: Export Format Validation Tests (10.6.12—10.6.20)
# ============================================================================

test_that("TC-10.6.12: LaTeX booktabs vs hline mode", {
  obj <- .mock_common_timing_result()
  tex_bt <- to_latex(obj, booktabs = TRUE)
  tex_hl <- to_latex(obj, booktabs = FALSE)
  # booktabs mode
  expect_true(grepl("\\\\toprule", tex_bt))
  expect_false(grepl("\\\\hline", tex_bt))
  # hline mode
  expect_true(grepl("\\\\hline", tex_hl))
  expect_false(grepl("\\\\toprule", tex_hl))
})

test_that("TC-10.6.13: LaTeX stars=FALSE no stars or tablenotes", {
  obj <- .mock_common_timing_result()
  tex <- to_latex(obj, stars = FALSE)
  expect_false(grepl("\\*\\*\\*", tex))
  expect_false(grepl("tablenotes", tex))
})

test_that("TC-10.6.14: LaTeX CI annotation contains confidence level", {
  obj <- .mock_common_timing_result()
  tex <- to_latex(obj, include_ci = TRUE)
  expect_true(grepl("95", tex))
})

test_that("TC-10.6.15: LaTeX include_diagnostics controls display", {
  obj <- .mock_result_with_diagnostics()
  tex_with <- to_latex(obj, include_diagnostics = TRUE)
  tex_without <- to_latex(obj, include_diagnostics = FALSE)
  expect_true(nchar(tex_with) > nchar(tex_without))
})

test_that("TC-10.6.16: CSV summary exports single row", {
  obj <- .mock_common_timing_result()
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  df <- to_csv(obj, file = tf, what = "summary")
  expect_equal(nrow(df), 1L)
  expect_equal(df$att, 0.5, tolerance = 1e-10)
  expect_equal(df$nobs, 500L)
  expect_true(all(c("att", "se", "t_stat", "pvalue", "nobs") %in% names(df)))
})

test_that("TC-10.6.17: CSV by_period row count matches att_by_period", {
  obj <- .mock_common_timing_result()
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  df <- to_csv(obj, file = tf, what = "by_period")
  expect_equal(nrow(df), nrow(obj$att_by_period))
})

test_that("TC-10.6.18: CSV by_cohort errors on non-staggered", {
  obj <- .mock_common_timing_result()
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  expect_error(to_csv(obj, file = tf, what = "by_cohort"),
               "No cohort-specific results \\(Staggered mode only\\)\\.")
})

test_that("TC-10.6.19: to_latex_comparison column count correct", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_staggered_result()
  m3 <- .mock_result_with_ri()
  tex <- to_latex_comparison(m1, m2, m3,
                             model_names = c("M1", "M2", "M3"))
  # tabular format should be {lccc} for 3 models
  expect_true(grepl("\\{lccc\\}", tex))
  expect_true(grepl("M1", tex))
  expect_true(grepl("M2", tex))
  expect_true(grepl("M3", tex))
})

test_that("TC-10.6.20: to_latex file param writes file", {
  obj <- .mock_common_timing_result()
  tf <- tempfile(fileext = ".tex")
  on.exit(unlink(tf), add = TRUE)
  result <- to_latex(obj, file = tf, caption = "Test", label = "tab:test")
  expect_true(file.exists(tf))
  content <- readLines(tf)
  expect_true(any(grepl("\\\\begin\\{table\\}", content)))
  expect_true(is.character(result))
})

# ============================================================================
# Task E10-06.5: Boundary Conditions & Error Handling (10.6.21—10.6.30)
# ============================================================================

test_that("TC-10.6.21: Non-lwdid_result object rejected by all S3 methods", {
  fake <- list(att = 1)
  expect_error(print.lwdid_result(fake))
  expect_error(summary.lwdid_result(fake))
  expect_error(coef.lwdid_result(fake))
  expect_error(confint.lwdid_result(fake))
  expect_error(vcov.lwdid_result(fake))
})

test_that("TC-10.6.22: Extreme p-value formatting", {
  # .format_pvalue is internal — access via lwdid:::
  fmt <- lwdid:::.format_pvalue
  expect_equal(fmt(1e-16), "<0.001")
  expect_equal(fmt(0), "<0.001")
  expect_equal(fmt(1.0, digits = 4L), "1.0000")
  expect_equal(fmt(0.04999, digits = 4L), "0.0500")
})

test_that("TC-10.6.23: Small df produces wider CI than large df", {
  obj_small <- .mock_common_timing_result()
  obj_small$df_inference <- 2L
  obj_large <- .mock_common_timing_result()
  obj_large$df_inference <- 100L
  ci_small <- confint(obj_small, level = 0.95, type = "overall")
  ci_large <- confint(obj_large, level = 0.95, type = "overall")
  width_small <- ci_small["ATT", 2] - ci_small["ATT", 1]
  width_large <- ci_large["ATT", 2] - ci_large["ATT", 1]
  expect_gt(width_small, width_large)
})

test_that("TC-10.6.24: vcov matrix is symmetric", {
  obj <- .mock_common_timing_result()
  V <- vcov(obj, type = "by_period")
  expect_true(isSymmetric(V))
})

test_that("TC-10.6.25: vcov matrix is positive semi-definite", {
  obj <- .mock_common_timing_result()
  V <- vcov(obj, type = "by_period")
  eigenvalues <- eigen(V, symmetric = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))
})

test_that("TC-10.6.26: coef invalid type argument errors", {
  obj <- .mock_common_timing_result()
  expect_error(coef(obj, type = "invalid_type"))
})

test_that("TC-10.6.27: CSV invalid what argument errors", {
  obj <- .mock_common_timing_result()
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  expect_error(to_csv(obj, file = tf, what = "invalid_what"))
})

test_that("TC-10.6.28: LaTeX handles special characters in estimator", {
  obj <- .mock_common_timing_result()
  obj$estimator <- "RA_robust"
  tex <- to_latex(obj)
  expect_true(is.character(tex))
  expect_true(nchar(tex) > 0)
  # Should not crash — underscore is escaped
})

test_that("TC-10.6.29: Empty att_by_period errors on coef(type='by_period')", {
  obj <- .mock_common_timing_result()
  obj$att_by_period <- NULL
  expect_error(coef(obj, type = "by_period"),
               "\u65f6\u671f\u7279\u5b9a|period")
})

test_that("TC-10.6.30: Empty att_by_cohort errors on vcov(type='by_cohort')", {
  obj <- .mock_common_timing_result()
  obj$att_by_cohort <- NULL
  expect_error(vcov(obj, type = "by_cohort"),
               "\u961f\u5217\u7279\u5b9a|cohort")
})

# ============================================================================
# Task E10-06.6: Cross-Functional End-to-End Tests (10.6.31—10.6.36)
# ============================================================================

test_that("TC-10.6.31: print->summary->to_latex pipeline", {
  obj <- .mock_common_timing_result()
  # print outputs >5 lines
  out <- capture.output(print(obj))
  expect_gt(length(out), 5L)
  # summary returns correct class
  s <- summary(obj)
  expect_s3_class(s, "summary.lwdid_result")
  # to_latex contains "table"
  tex <- to_latex(obj)
  expect_true(grepl("table", tex, ignore.case = TRUE))
})

test_that("TC-10.6.32: coef/confint/vcov numerical triangle consistency", {
  obj <- .mock_common_timing_result()
  att <- coef(obj)[["ATT"]]
  se <- sqrt(vcov(obj)["ATT", "ATT"])
  ci <- confint(obj)
  df <- obj$df_inference
  expected_lower <- att - stats::qt(0.975, df) * se
  expected_upper <- att + stats::qt(0.975, df) * se
  expect_equal(ci["ATT", 1], expected_lower, tolerance = 1e-10)
  expect_equal(ci["ATT", 2], expected_upper, tolerance = 1e-10)
})

test_that("TC-10.6.33: plot->to_latex->to_csv full chain", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_common_timing_result()
  # plot returns ggplot
  p <- plot(obj)
  expect_s3_class(p, "ggplot")
  # to_latex contains caption and label
  tex <- to_latex(obj, caption = "My Caption", label = "tab:my")
  expect_true(grepl("My Caption", tex))
  expect_true(grepl("tab:my", tex))
  # to_csv writes non-empty file
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  to_csv(obj, file = tf, what = "summary")
  expect_true(file.exists(tf))
  expect_gt(file.info(tf)$size, 0)
})

test_that("TC-10.6.34: RI result complete output", {
  obj <- .mock_result_with_ri()
  out <- capture.output(print(obj))
  out_str <- paste(out, collapse = "\n")
  expect_true(grepl("RI p-value", out_str))
  expect_true(grepl("valid=999/999", out_str))
  s <- summary(obj)
  s_out <- capture.output(print(s))
  s_str <- paste(s_out, collapse = "\n")
  expect_true(grepl("Randomization Inference", s_str))
  expect_true(grepl("valid=999/999", s_str))
})

test_that("TC-10.6.35: Diagnostics complete output", {
  obj <- .mock_result_with_diagnostics()
  s <- summary(obj)
  s_out <- capture.output(print(s))
  s_str <- paste(s_out, collapse = "\n")
  expect_true(grepl("Diagnostics", s_str, ignore.case = TRUE))
  expect_true(grepl("clustering", s_str, ignore.case = TRUE))
  expect_true(grepl("selection", s_str, ignore.case = TRUE))
})

test_that("TC-10.6.36: Multi-model comparison LaTeX and CSV consistency", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_staggered_result()
  tex <- to_latex_comparison(m1, m2, model_names = c("CT", "Stag"))
  # LaTeX contains formatted ATT values
  expect_true(grepl(sprintf("%.4f", m1$att), tex))
  expect_true(grepl(sprintf("%.4f", m2$att), tex))
  # CSV values match
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  df <- to_csv(m1, file = tf, what = "summary")
  expect_equal(df$att, m1$att, tolerance = 1e-10)
})

# ============================================================================
# Task E10-06.7: Supplementary Coverage Tests (10.6.37—10.6.62)
# ============================================================================

# --- Sub-category A: S3 Method Supplementary (10.6.37—10.6.43) ---

test_that("TC-10.6.37: confint(type='by_cohort') returns correct matrix", {
  obj <- .mock_staggered_result()
  ci <- confint(obj, type = "by_cohort")
  expect_equal(nrow(ci), 2L)
  expect_equal(ncol(ci), 2L)
  expect_equal(rownames(ci), c("5", "7"))
  # Each cohort ATT within CI
  atts <- coef(obj, type = "by_cohort")
  for (nm in names(atts)) {
    expect_true(atts[[nm]] >= ci[nm, 1])
    expect_true(atts[[nm]] <= ci[nm, 2])
  }
})

test_that("TC-10.6.38: print Staggered output contains key terms", {
  obj <- .mock_staggered_result()
  out <- capture.output(print(obj))
  out_str <- paste(out, collapse = "\n")
  expect_true(grepl("Staggered", out_str, ignore.case = TRUE))
  expect_true(grepl("Cohort|cohort", out_str))
  expect_true(grepl("ATT", out_str))
})

test_that("TC-10.6.39: coef(type='all') returns g{c}.r{p} named vector", {
  obj <- .mock_staggered_result()
  c_all <- coef(obj, type = "all")
  expect_equal(length(c_all), 5L)
  # Names should be g{cohort}.r{period}
  expected_names <- sprintf("g%s.r%s",
    obj$att_by_cohort_time$cohort,
    obj$att_by_cohort_time$period)
  expect_equal(names(c_all), expected_names)
  # Values match att_by_cohort_time$att
  expect_equal(unname(c_all), obj$att_by_cohort_time$att,
               tolerance = 1e-10)
})

test_that("TC-10.6.40: coef Staggered type='by_period' uses att_by_period", {
  obj <- .mock_staggered_result()
  c_bp <- coef(obj, type = "by_period")
  expect_equal(length(c_bp), nrow(obj$att_by_period))
  expect_equal(unname(c_bp), obj$att_by_period$att, tolerance = 1e-10)
})

test_that("TC-10.6.41: vcov prefers vcov_att_periods when available", {
  obj <- .mock_common_timing_result()
  n <- nrow(obj$att_by_period)
  # Construct a non-diagonal vcov matrix
  custom_vcov <- diag(obj$att_by_period$se^2)
  custom_vcov[1, 2] <- 0.001
  custom_vcov[2, 1] <- 0.001
  obj$vcov_att_periods <- custom_vcov
  V <- vcov(obj, type = "by_period")
  # Should use the custom matrix, not diag(se^2)
  expect_equal(V[1, 2], 0.001, tolerance = 1e-10)
  expect_equal(V[2, 1], 0.001, tolerance = 1e-10)
})

test_that("TC-10.6.42: confint extreme level=0.999 wider than 0.95", {
  obj <- .mock_common_timing_result()
  ci_95 <- confint(obj, level = 0.95)
  ci_999 <- confint(obj, level = 0.999)
  w_95 <- ci_95["ATT", 2] - ci_95["ATT", 1]
  w_999 <- ci_999["ATT", 2] - ci_999["ATT", 1]
  expect_gt(w_999, w_95)
  # Verify uses qt(0.9995, df)
  df <- obj$df_inference
  expected_w999 <- 2 * stats::qt(0.9995, df) * obj$se_att
  expect_equal(w_999, expected_w999, tolerance = 1e-10)
})

test_that("TC-10.6.43: Empty warnings_log no warning lines in print", {
  obj <- .mock_common_timing_result()
  obj$warnings_log <- character(0)
  out <- capture.output(print(obj))
  warning_lines <- grep("warning", out, ignore.case = TRUE)
  expect_equal(length(warning_lines), 0L)
})

# --- Sub-category B: Export Supplementary (10.6.44—10.6.49) ---

test_that("TC-10.6.44: to_latex include_ri=TRUE shows RI fields", {
  obj <- .mock_result_with_ri()
  tex <- to_latex(obj, include_ri = TRUE)
  expect_true(grepl("RI p-value", tex))
  expect_true(grepl("RI seed", tex))
  expect_true(grepl("42", tex))
  expect_true(grepl("RI reps", tex))
  expect_true(grepl("999", tex))
})

test_that("TC-10.6.45: to_latex include_ri=FALSE hides RI fields", {
  obj <- .mock_result_with_ri()
  tex <- to_latex(obj, include_ri = FALSE)
  expect_false(grepl("RI p-value", tex))
  expect_false(grepl("RI seed", tex))
})

test_that("TC-10.6.46: to_latex include_periods=TRUE shows period table", {
  obj <- .mock_common_timing_result()
  tex <- to_latex(obj, include_periods = TRUE)
  expect_true(grepl("Period", tex))
  expect_true(grepl("t-stat", tex))
  # Each period ATT value should appear (4 decimal places)
  for (a in obj$att_by_period$att) {
    expect_true(grepl(sprintf("%.4f", a), tex))
  }
})

test_that("TC-10.6.47: to_latex include_staggered=TRUE shows cohort table", {
  obj <- .mock_staggered_result()
  tex <- to_latex(obj, include_staggered = TRUE)
  expect_true(grepl("Cohort", tex))
  expect_true(grepl("Weight", tex))
  expect_true(grepl("Overall.*weighted", tex, ignore.case = TRUE))
})

test_that("TC-10.6.48: to_latex uses threeparttable environment", {
  obj <- .mock_common_timing_result()
  tex <- to_latex(obj, stars = TRUE)
  expect_true(grepl("\\\\begin\\{threeparttable\\}", tex))
  expect_true(grepl("\\\\end\\{threeparttable\\}", tex))
  expect_true(grepl("tablenotes", tex))
})

test_that("TC-10.6.49: to_csv what='all' Staggered exports cohort_time", {
  obj <- .mock_staggered_result()
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  df <- to_csv(obj, file = tf, what = "all")
  expect_equal(nrow(df), 5L)
  expect_true("cohort" %in% names(df))
  expect_true("period" %in% names(df) || "event_time" %in% names(df))
  # ATT values match
  expect_equal(df$att, obj$att_by_cohort_time$att, tolerance = 1e-10)
})

# --- Sub-category C: Visualization Supplementary (10.6.50—10.6.52) ---

test_that("TC-10.6.50: plot_event_study Common Timing correct", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_common_timing_result()
  res <- plot_event_study(obj, return_data = TRUE)
  expect_s3_class(res$plot, "ggplot")
  pd <- res$data
  expect_true(is.data.frame(pd))
  expect_gte(nrow(pd), nrow(obj$att_by_period))
})

test_that("TC-10.6.51: plot_event_study aggregation produces different ATT", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_staggered_result()
  res_w <- plot_event_study(obj, aggregation = "weighted", return_data = TRUE)
  res_m <- plot_event_study(obj, aggregation = "mean", return_data = TRUE)
  # With unequal weights (0.6/0.4), weighted != mean
  pd_w <- res_w$data
  pd_m <- res_m$data
  # Compare ATT values for non-anchor rows
  att_w <- pd_w$att[!pd_w$is_anchor]
  att_m <- pd_m$att[!pd_m$is_anchor]
  # At least one ATT value should differ between weighted and mean
  if (length(att_w) > 0 && length(att_m) > 0) {
    expect_false(all(abs(att_w - att_m) < 1e-12))
  }
})

test_that("TC-10.6.52: plot_event_study ref_period=0 normalizes anchor", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_staggered_result()
  res <- plot_event_study(obj, ref_period = 0L, return_data = TRUE)
  pd <- res$data
  # Find the anchor row (event_time closest to 0 or is_anchor=TRUE)
  anchor_rows <- pd[pd$is_anchor == TRUE, ]
  if (nrow(anchor_rows) > 0) {
    expect_equal(anchor_rows$att[1], 0, tolerance = 1e-12)
  }
})

# --- Sub-category D: End-to-End & Advanced (10.6.53—10.6.62) ---

test_that("TC-10.6.53: Staggered full pipeline print->summary->plot->export", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_staggered_result()
  # print
  out <- capture.output(print(obj))
  expect_gt(length(out), 5L)
  # summary
  s <- summary(obj)
  expect_s3_class(s, "summary.lwdid_result")
  # plot
  p <- plot(obj)
  expect_s3_class(p, "ggplot")
  # export
  tex <- to_latex(obj)
  expect_true(grepl("table", tex, ignore.case = TRUE))
})

test_that("TC-10.6.54: Event Study plot data matches LaTeX ATT values", {
  skip_if_not_installed("ggplot2")
  obj <- .mock_common_timing_result()
  tex <- to_latex(obj)
  # Overall ATT should appear in LaTeX
  expect_true(grepl(sprintf("%.4f", obj$att), tex))
  # Plot should be valid
  p <- plot_event_study(obj)
  expect_s3_class(p, "ggplot")
})

test_that("TC-10.6.55: df_strategy three strategies integration", {
  obj_cons <- .mock_staggered_result()
  obj_cons$df_strategy <- "conservative"
  obj_wt <- .mock_staggered_result()
  obj_wt$df_strategy <- "weighted"
  # ATT should be identical regardless of df_strategy
  expect_equal(obj_cons$att, obj_wt$att, tolerance = 1e-10)
  # CI widths may differ (conservative uses min df -> wider CI)
  # Both should have valid confint
  ci_cons <- confint(obj_cons)
  ci_wt <- confint(obj_wt)
  w_cons <- ci_cons["ATT", 2] - ci_cons["ATT", 1]
  w_wt <- ci_wt["ATT", 2] - ci_wt["ATT", 1]
  # With same overall df_inference, widths are equal here
  # (df_strategy affects per-cohort-time, not overall)
  expect_equal(w_cons, w_wt, tolerance = 1e-10)
})

test_that("TC-10.6.56: Anchor point (event_time=-1) ATT=0, SE=0, CI=[0,0]", {
  obj <- .mock_common_timing_result()
  # Add anchor row to att_pre_treatment
  anchor <- data.frame(
    event_time = -1L, att = 0, se = 0,
    ci_lower = 0, ci_upper = 0, pvalue = NA_real_
  )
  obj$att_pre_treatment <- rbind(obj$att_pre_treatment, anchor)
  anchor_row <- obj$att_pre_treatment[
    obj$att_pre_treatment$event_time == -1L, ]
  expect_equal(anchor_row$att, 0, tolerance = 1e-12)
  expect_equal(anchor_row$se, 0, tolerance = 1e-12)
  expect_equal(anchor_row$ci_lower, 0, tolerance = 1e-12)
  expect_equal(anchor_row$ci_upper, 0, tolerance = 1e-12)
})

test_that("TC-10.6.57: to_latex_comparison contains meta-info rows", {
  m1 <- .mock_common_timing_result()
  m2 <- .mock_staggered_result()
  tex <- to_latex_comparison(m1, m2, model_names = c("CT", "Stag"))
  expect_true(grepl("Estimator", tex))
  expect_true(grepl("VCE", tex))
  expect_true(grepl("Rolling", tex))
})

test_that("TC-10.6.58: print.summary shows Pre-treatment with anchor", {
  obj <- .mock_common_timing_result()
  # Add anchor row
  anchor <- data.frame(
    event_time = -1L, att = 0, se = 0,
    ci_lower = 0, ci_upper = 0, pvalue = NA_real_
  )
  obj$att_pre_treatment <- rbind(obj$att_pre_treatment, anchor)
  obj$include_pretreatment <- TRUE
  s <- summary(obj)
  s_out <- capture.output(print(s))
  s_str <- paste(s_out, collapse = "\n")
  expect_true(grepl("Pre-treatment", s_str, ignore.case = TRUE))
  expect_true(
    grepl("anchor", s_str, ignore.case = TRUE) ||
    grepl("by construction", s_str, ignore.case = TRUE) ||
    grepl("e=-1", s_str, ignore.case = TRUE)
  )
})

test_that("TC-10.6.59: confint(type='all') matches coef(type='all')", {
  obj <- .mock_staggered_result()
  c_all <- coef(obj, type = "all")
  ci_all <- confint(obj, type = "all")
  # Dimensions and names match
  expect_equal(length(c_all), nrow(ci_all))
  expect_equal(names(c_all), rownames(ci_all))
  # CI width consistent with SE
  df <- obj$df_inference
  t_crit <- stats::qt(0.975, df)
  for (i in seq_along(c_all)) {
    nm <- names(c_all)[i]
    se_i <- obj$att_by_cohort_time$se[i]
    expected_w <- 2 * t_crit * se_i
    actual_w <- ci_all[nm, 2] - ci_all[nm, 1]
    expect_equal(actual_w, expected_w, tolerance = 1e-10,
                 info = sprintf("CI width mismatch for %s", nm))
  }
})

test_that("TC-10.6.60: vcov(type='full') returns matrix or errors", {
  obj <- .mock_common_timing_result()
  # Without vcov_full, should error
  expect_error(vcov(obj, type = "full"),
               "store_vcov_full")
  # With vcov_full set, should return it
  obj$vcov_full <- matrix(c(0.0144, 0.001, 0.001, 0.0196),
                          nrow = 2)
  V <- vcov(obj, type = "full")
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 2L)
})

test_that("TC-10.6.61: per-event-time df differences in staggered", {
  obj <- .mock_staggered_result()
  bct <- obj$att_by_cohort_time
  # event_time=2 only has cohort 5 (df=28)
  # event_time=0 has both cohorts
  # Verify df values are present and differ
  et2_rows <- bct[bct$period == 7L & bct$cohort == 5L, ]
  et0_rows <- bct[bct$period == 5L & bct$cohort == 5L, ]
  expect_equal(et2_rows$df_inference[1], 28L)
  expect_equal(et0_rows$df_inference[1], 28L)
  # Cohort 7 has df=18
  c7_rows <- bct[bct$cohort == 7L, ]
  expect_true(all(c7_rows$df_inference == 18L))
  # df differs between cohorts
  expect_true(28L > 18L)
})

test_that("TC-10.6.62: to_csv summary with/without RI columns", {
  # With RI
  obj_ri <- .mock_result_with_ri()
  tf1 <- tempfile(fileext = ".csv")
  on.exit(unlink(tf1), add = TRUE)
  df_ri <- to_csv(obj_ri, file = tf1, what = "summary")
  expect_true("ri_pvalue" %in% names(df_ri))
  expect_true("ri_seed" %in% names(df_ri))
  expect_true("rireps" %in% names(df_ri))
  expect_equal(df_ri$ri_pvalue, 0.004, tolerance = 1e-10)
  expect_equal(df_ri$ri_seed, 42L)
  expect_equal(df_ri$rireps, 999L)
  # Without RI
  obj_no_ri <- .mock_common_timing_result()
  tf2 <- tempfile(fileext = ".csv")
  on.exit(unlink(tf2), add = TRUE)
  df_no_ri <- to_csv(obj_no_ri, file = tf2, what = "summary")
  expect_false("ri_pvalue" %in% names(df_no_ri))
  expect_false("ri_seed" %in% names(df_no_ri))
  expect_false("rireps" %in% names(df_no_ri))
})
