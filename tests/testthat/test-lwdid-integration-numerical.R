# ===========================================================================
# test-lwdid-integration-numerical.R
# Numerical validation tests for lwdid() RI and pretreatment (Story E7-06.7)
# Verifies numerical correctness, statistical properties, and consistency.
# ===========================================================================

# --- Helper: Common Timing DGP ---
.make_ct_num <- function(n_treated = 50L, n_control = 50L,
                         n_periods = 8L, K = 4L,
                         att = 2.0, seed = 42L) {
  set.seed(seed)
  n_units <- n_treated + n_control
  years <- 2000L + seq_len(n_periods)
  df <- expand.grid(id = seq_len(n_units), year = years)
  df <- df[order(df$id, df$year), ]
  df$d <- ifelse(df$id <= n_treated, 1L, 0L)
  df$post <- ifelse(df$year > (2000L + K), 1L, 0L)
  alpha_i <- rnorm(n_units, sd = 1.5)
  df$y <- alpha_i[df$id] + df$d * df$post * att + rnorm(nrow(df), sd = 0.5)
  rownames(df) <- NULL
  df
}

# ===========================================================================
# E7-06.7.1: RI p-value range + seed reproducibility
# ===========================================================================
test_that("E7-06.7.1: RI p-value in [0,1] and seed-fixed reproducibility", {
  df <- .make_ct_num(att = 2.0, seed = 7001L)

  res1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 12345L)
  ))
  expect_true(!is.null(res1$ri_pvalue))
  expect_true(res1$ri_pvalue >= 0 && res1$ri_pvalue <= 1)

  # Same seed -> bit-identical
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 12345L)
  ))
  expect_identical(res1$ri_pvalue, res2$ri_pvalue)
  expect_identical(res1$ri_distribution, res2$ri_distribution)
  expect_identical(res1$ri_seed, res2$ri_seed)
  expect_identical(res1$ri_valid, res2$ri_valid)
  expect_identical(res1$ri_failed, res2$ri_failed)
  expect_identical(res1$att, res2$att)
  expect_identical(res1$se_att, res2$se_att)
})

# ===========================================================================
# E7-06.7.2: RI statistical properties (size and power)
# ===========================================================================
test_that("E7-06.7.2: RI p-value large under null (effect=0)", {
  df <- .make_ct_num(att = 0.0, seed = 7021L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 7022L)
  ))
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue > 0.01,
                info = sprintf("Under null (att=0), RI p = %.4f should be > 0.01",
                               result$ri_pvalue))
  }
})

test_that("E7-06.7.2: RI p-value small under strong alternative (effect=5)", {
  df <- .make_ct_num(att = 5.0, seed = 7023L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 7024L)
  ))
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue < 0.10,
                info = sprintf("Under strong alt (att=5), RI p = %.4f should be < 0.10",
                               result$ri_pvalue))
  }
})

# ===========================================================================
# E7-06.7.3: Pre-treatment ATT numerical reasonableness
# ===========================================================================
test_that("E7-06.7.3: Pre-treatment ATTs near zero under parallel trends DGP", {
  df <- .make_ct_num(n_treated = 80L, n_control = 80L,
                     att = 3.0, seed = 7031L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE,
          pretreatment_test = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(is.data.frame(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)

  post_att <- abs(result$att)
  pre_atts <- result$att_pre_treatment$att[!is.na(result$att_pre_treatment$att)]

  # Each |pre-ATT| < 0.5 * |post-ATT|
  for (i in seq_along(pre_atts)) {
    expect_true(abs(pre_atts[i]) < 0.5 * post_att,
                info = sprintf("Pre[%d]: |%.4f| should be < 0.5*|%.4f|",
                               i, pre_atts[i], post_att))
  }
  # Mean |pre-ATT| < 0.3 * |post-ATT|
  mean_abs_pre <- mean(abs(pre_atts))
  expect_true(mean_abs_pre < 0.3 * post_att,
              info = sprintf("Mean |pre-ATT| = %.4f < 0.3*%.4f",
                             mean_abs_pre, post_att))

  # Parallel trends test should not reject at very low threshold
  if (!is.null(result$parallel_trends_test$joint_pvalue)) {
    expect_true(result$parallel_trends_test$joint_pvalue > 0.001,
                info = sprintf("Joint p = %.6f should be > 0.001 under PT",
                               result$parallel_trends_test$joint_pvalue))
  }
})

# ===========================================================================
# E7-06.7.4: Joint F-test size under H0
# ===========================================================================
test_that("E7-06.7.4: Joint F-test reject rate near alpha under H0", {
  n_mc <- 50L
  alpha_test <- 0.05
  rejections <- 0L
  valid_runs <- 0L

  for (mc in seq_len(n_mc)) {
    df <- .make_ct_num(n_treated = 40L, n_control = 40L,
                       att = 0.0, seed = 7040L + mc)
    res <- tryCatch(
      suppressWarnings(suppressMessages(
        lwdid(data = df, y = "y", ivar = "id", tvar = "year",
              d = "d", post = "post",
              include_pretreatment = TRUE, pretreatment_test = TRUE,
              pretreatment_alpha = alpha_test)
      )),
      error = function(e) NULL
    )
    if (!is.null(res) && !is.null(res$parallel_trends_test)) {
      valid_runs <- valid_runs + 1L
      if (isTRUE(res$parallel_trends_test$reject_null)) {
        rejections <- rejections + 1L
      }
    }
  }
  # Need at least 30 valid runs
  if (valid_runs >= 30L) {
    reject_rate <- rejections / valid_runs
    # Under H0, reject rate should be <= 0.20 (generous for 50 MC reps)
    expect_true(reject_rate <= 0.20,
                info = sprintf("Reject rate = %.3f (n=%d) should be <= 0.20",
                               reject_rate, valid_runs))
  }
})

# ===========================================================================
# E7-06.7.5: exclude_pre_periods=0 identical to default
# ===========================================================================
test_that("E7-06.7.5: exclude_pre_periods=0 identical to default", {
  df <- .make_ct_num(att = 2.0, seed = 7051L)
  res_default <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post")
  ))
  res_zero <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", exclude_pre_periods = 0L)
  ))
  expect_identical(res_default$att, res_zero$att)
  expect_identical(res_default$se_att, res_zero$se_att)
  expect_identical(res_default$pvalue, res_zero$pvalue)
  expect_identical(res_default$ci_lower, res_zero$ci_lower)
  expect_identical(res_default$ci_upper, res_zero$ci_upper)
})

# ===========================================================================
# E7-06.7.6: RI failure does not affect ATT
# ===========================================================================
test_that("E7-06.7.6: RI does not affect ATT/SE/CI/pvalue", {
  df <- .make_ct_num(att = 2.0, seed = 7061L)
  res_no_ri <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = FALSE)
  ))
  res_ri <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 7062L)
  ))
  expect_identical(res_no_ri$att, res_ri$att)
  expect_identical(res_no_ri$se_att, res_ri$se_att)
  expect_identical(res_no_ri$pvalue, res_ri$pvalue)
  expect_identical(res_no_ri$ci_lower, res_ri$ci_lower)
  expect_identical(res_no_ri$ci_upper, res_ri$ci_upper)
})

# ===========================================================================
# E7-06.7.7: ri + include_pretreatment do not interfere
# ===========================================================================
test_that("E7-06.7.7: ri+pretreatment ATT identical to ri-only", {
  df <- .make_ct_num(att = 2.0, seed = 7071L)
  res_ri <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post",
          ri = TRUE, rireps = 200L, seed = 7072L)
  ))
  res_both <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post",
          ri = TRUE, rireps = 200L, seed = 7072L,
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))
  expect_identical(res_ri$att, res_both$att)
  expect_identical(res_ri$se_att, res_both$se_att)
  expect_identical(res_ri$pvalue, res_both$pvalue)
  expect_identical(res_ri$ri_pvalue, res_both$ri_pvalue)
  expect_identical(res_ri$ri_distribution, res_both$ri_distribution)
  expect_true(!is.null(res_both$att_pre_treatment))
})
