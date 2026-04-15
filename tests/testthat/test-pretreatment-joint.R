# ============================================================================
# test-pretreatment-joint.R
# Unit tests for test_parallel_trends_joint()
# Story E7-05: Parallel Trends Joint Test
# Paper: LW2025 Equation (2.3) Parallel Trends Hypothesis
# ============================================================================

# ============================================================================
# Section 0: Test Helper — Build pre_treatment_effects data.frame
# ============================================================================

#' Build a pre_treatment_effects data.frame for testing
#' @param att numeric vector of ATT estimates
#' @param se numeric vector of standard errors
#' @param t_stat numeric vector of t-statistics (default att/se)
#' @param pvalue numeric vector of p-values (default 2*pt(-abs(t),df=48))
#' @param is_anchor logical vector
#' @param n_treated integer vector
#' @param n_control integer vector
#' @param event_time integer vector
#' @param cohort integer vector (NA for common timing)
#' @param period integer vector
#' @return data.frame with 10 required columns
make_pre_effects <- function(att, se,
                              t_stat = NULL,
                              pvalue = NULL,
                              is_anchor = NULL,
                              n_treated = NULL,
                              n_control = NULL,
                              event_time = NULL,
                              cohort = NULL,
                              period = NULL) {
  K <- length(att)
  if (is.null(t_stat)) t_stat <- att / se
  if (is.null(pvalue)) pvalue <- 2 * stats::pt(-abs(t_stat), df = 48)
  if (is.null(is_anchor)) is_anchor <- rep(FALSE, K)
  if (is.null(n_treated)) n_treated <- rep(25L, K)
  if (is.null(n_control)) n_control <- rep(25L, K)
  if (is.null(event_time)) event_time <- seq(-K, -1L)
  if (is.null(cohort)) cohort <- rep(NA_integer_, K)
  if (is.null(period)) period <- seq_len(K)
  data.frame(
    att = att, se = se, t_stat = t_stat, pvalue = pvalue,
    is_anchor = is_anchor, n_treated = n_treated,
    n_control = n_control, event_time = event_time,
    cohort = cohort, period = period,
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# Section 1: Basic Functionality (TC-7.5.1 to TC-7.5.5)
# ============================================================================

test_that("TC-7.5.1: parallel trends hold DGP — does not reject H0", {
  # Small ATTs near zero → F-test should not reject
  set.seed(42)
  att_vals <- rnorm(5, mean = 0, sd = 0.1)
  se_vals <- rep(1.0, 5)
  pe <- make_pre_effects(att = att_vals, se = se_vals,
                          n_treated = rep(25L, 5),
                          n_control = rep(25L, 5))
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_false(res$reject_null)
  expect_true(res$joint_pvalue > 0.05)
  expect_equal(res$n_pre_periods, 5L)
  expect_equal(res$test_type, "f")
  expect_equal(res$alpha, 0.05)
  expect_equal(length(res), 14L)
})

test_that("TC-7.5.2: parallel trends violated DGP — rejects H0", {
  # Large ATTs → F-test should reject
  att_vals <- c(3.0, 4.0, 5.0, 3.5, 4.5)
  se_vals <- rep(1.0, 5)
  pe <- make_pre_effects(att = att_vals, se = se_vals,
                          n_treated = rep(25L, 5),
                          n_control = rep(25L, 5))
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_true(res$reject_null)
  expect_true(res$joint_pvalue < 0.05)
})

test_that("TC-7.5.3: single pre-period (K=1) — F degenerates to t²", {
  att_val <- 3.0
  se_val <- 1.0
  t_val <- 3.0
  pe <- make_pre_effects(att = att_val, se = se_val, t_stat = t_val,
                          n_treated = 25L, n_control = 25L)
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f",
                                     min_pre_periods = 1L)
  # F = sum(t²)/K = 9/1 = 9 = t²

  expect_equal(res$joint_stat, t_val^2, tolerance = 1e-12)
  expect_equal(res$joint_df1, 1L)
  expect_equal(res$n_pre_periods, 1L)
})

test_that("TC-7.5.4: multiple pre-periods (K=5) — F = sum(t²)/5", {
  t_stats <- c(2.0, 1.5, 0.5, -1.0, 0.8)
  att_vals <- t_stats * 1.0  # se = 1
  se_vals <- rep(1.0, 5)
  pe <- make_pre_effects(att = att_vals, se = se_vals, t_stat = t_stats,
                          n_treated = rep(25L, 5),
                          n_control = rep(25L, 5))
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  # sum(t²) = 4 + 2.25 + 0.25 + 1 + 0.64 = 8.14
  expected_sum_t2 <- sum(t_stats^2)
  expected_F <- expected_sum_t2 / 5
  expect_equal(res$joint_stat, expected_F, tolerance = 1e-10)
  expect_equal(res$joint_df1, 5L)
})

test_that("TC-7.5.5: all pre-period ATT=0 exactly — F=0, p≈1", {
  att_vals <- rep(0.0, 3)
  se_vals <- rep(1.0, 3)
  t_stats <- rep(0.0, 3)
  pe <- make_pre_effects(att = att_vals, se = se_vals, t_stat = t_stats,
                          n_treated = rep(25L, 3),
                          n_control = rep(25L, 3))
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_equal(res$joint_stat, 0.0, tolerance = 1e-12)
  expect_equal(res$joint_pvalue, 1.0, tolerance = 1e-10)
  expect_false(res$reject_null)
})


# ============================================================================
# Section 2: Exclusion Logic (TC-7.5.6 to TC-7.5.7)
# ============================================================================

test_that("TC-7.5.6: NA ATT/SE rows excluded, excluded_periods correct", {
  att_vals <- c(0.5, NA_real_, 0.3, 0.2)
  se_vals <- c(0.1, 0.1, NA_real_, 0.1)
  t_stats <- c(5.0, NA_real_, NA_real_, 2.0)
  pvals <- c(0.001, NA_real_, NA_real_, 0.05)
  pe <- make_pre_effects(
    att = att_vals, se = se_vals, t_stat = t_stats, pvalue = pvals,
    event_time = c(-4L, -3L, -2L, -1L),
    is_anchor = rep(FALSE, 4)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  # Rows 2 (NA att) and 3 (NA se) excluded

  expect_equal(res$n_pre_periods, 2L)
  expect_true(-3L %in% res$excluded_periods)
  expect_true(-2L %in% res$excluded_periods)
  expect_equal(length(res$excluded_periods), 2L)
})

test_that("TC-7.5.7: anchor point (is_anchor=TRUE) excluded from test", {
  att_vals <- c(0.5, 0.0, 0.3)
  se_vals <- c(0.1, 0.0, 0.1)
  t_stats <- c(5.0, NaN, 3.0)
  pvals <- c(0.001, NaN, 0.01)
  pe <- make_pre_effects(
    att = att_vals, se = se_vals, t_stat = t_stats, pvalue = pvals,
    event_time = c(-3L, -1L, -2L),
    is_anchor = c(FALSE, TRUE, FALSE)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  # Anchor (event_time=-1) excluded
  expect_equal(res$n_pre_periods, 2L)
  expect_true(-1L %in% res$excluded_periods)
  # Only non-anchor rows in individual_tests
  expect_equal(nrow(res$individual_tests), 2L)
  expect_false(-1L %in% res$individual_tests$event_time)
})

# ============================================================================
# Section 3: Wald Test & Parameter Tests (TC-7.5.8 to TC-7.5.9)
# ============================================================================

test_that("TC-7.5.8: Wald test — stat = sum(t²), chi-squared p-value", {
  t_stats <- c(2.0, 1.5, 0.5)
  att_vals <- t_stats
  se_vals <- rep(1.0, 3)
  pe <- make_pre_effects(att = att_vals, se = se_vals, t_stat = t_stats,
                          n_treated = rep(25L, 3),
                          n_control = rep(25L, 3))
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "wald")
  # W = 4 + 2.25 + 0.25 = 6.5
  expect_equal(res$joint_stat, 6.5, tolerance = 1e-12)
  expect_equal(res$joint_df1, 3L)
  expect_equal(res$joint_df2, 0L)
  expect_equal(res$test_type, "wald")
  # p-value from chi-squared(3)
  expected_p <- stats::pchisq(6.5, df = 3, lower.tail = FALSE)
  expect_equal(res$joint_pvalue, expected_p, tolerance = 1e-12)
})

test_that("TC-7.5.9: min_pre_periods insufficient triggers warning", {
  att_vals <- c(0.5)
  se_vals <- c(0.1)
  pe <- make_pre_effects(att = att_vals, se = se_vals)
  expect_warning(
    test_parallel_trends_joint(pe, min_pre_periods = 2L),
    class = "lwdid_insufficient_pre_periods"
  )
})

# ============================================================================
# Section 4: Edge Cases & Robustness (TC-7.5.10 to TC-7.5.13)
# ============================================================================

test_that("TC-7.5.10: staggered mode (cohort = actual values) works", {
  att_vals <- c(0.3, 0.2, 0.4, 0.1)
  se_vals <- rep(0.1, 4)
  pe <- make_pre_effects(
    att = att_vals, se = se_vals,
    cohort = c(5L, 5L, 7L, 7L),
    event_time = c(-3L, -2L, -3L, -2L),
    n_treated = rep(10L, 4),
    n_control = rep(20L, 4)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_equal(res$n_pre_periods, 4L)
  expect_true(all(res$individual_tests$cohort %in% c(5L, 7L)))
})

test_that("TC-7.5.11: K=0 returns NA result with correct structure", {
  # All rows are anchors → K=0
  pe <- make_pre_effects(
    att = c(0.0, 0.0), se = c(0.0, 0.0),
    is_anchor = c(TRUE, TRUE),
    event_time = c(-2L, -1L)
  )
  expect_warning(
    res <- test_parallel_trends_joint(pe, alpha = 0.05),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_true(is.na(res$joint_stat))
  expect_true(is.na(res$joint_pvalue))
  expect_false(res$reject_null)
  expect_equal(res$n_pre_periods, 0L)
  # individual_tests: empty data.frame with 10 columns
  expect_equal(nrow(res$individual_tests), 0L)
  expect_equal(ncol(res$individual_tests), 10L)
  expected_cols <- c("event_time", "cohort", "period", "att", "se",
                      "t_stat", "pvalue", "significant",
                      "n_treated", "n_control")
  expect_equal(names(res$individual_tests), expected_cols)
  # Full 14 fields
  expect_equal(length(res), 14L)
})

test_that("TC-7.5.12: individual_tests has significant column, values correct", {
  t_stats <- c(3.0, 0.5, 2.5)
  att_vals <- t_stats
  se_vals <- rep(1.0, 3)
  # Manually set pvalues: first and third significant at 0.05
  pvals <- c(0.003, 0.62, 0.013)
  pe <- make_pre_effects(att = att_vals, se = se_vals, t_stat = t_stats,
                          pvalue = pvals)
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_true("significant" %in% names(res$individual_tests))
  expect_equal(res$individual_tests$significant,
               c(TRUE, FALSE, TRUE))
  expect_equal(res$n_significant_05, 2L)
})

test_that("TC-7.5.13: reject_null is FALSE when joint_pvalue is NA", {
  # This tests the NaN guard: if somehow pvalue is NA, reject_null = FALSE
  # We can't easily make pf() return NA, but we test the K=0 path
  pe <- make_pre_effects(
    att = c(NA_real_), se = c(NA_real_),
    t_stat = c(NA_real_), pvalue = c(NA_real_),
    is_anchor = c(FALSE)
  )
  expect_warning(
    res <- test_parallel_trends_joint(pe, alpha = 0.05),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_false(res$reject_null)
})


# ============================================================================
# Section 5: Supplementary Edge Tests (TC-7.5.14 to TC-7.5.17)
# ============================================================================

test_that("TC-7.5.14: empty input (nrow=0) triggers stop_lwdid()", {
  empty_pe <- data.frame(
    att = numeric(0), se = numeric(0), t_stat = numeric(0),
    pvalue = numeric(0), is_anchor = logical(0),
    n_treated = integer(0), n_control = integer(0),
    event_time = integer(0), cohort = integer(0),
    period = integer(0), stringsAsFactors = FALSE
  )
  expect_error(
    test_parallel_trends_joint(empty_pe),
    class = "lwdid_invalid_param"
  )
})

test_that("TC-7.5.15: SE all zero — all rows excluded", {
  pe <- make_pre_effects(
    att = c(0.5, 0.3, 0.2),
    se = c(0.0, 0.0, 0.0),
    t_stat = c(Inf, Inf, Inf),
    pvalue = c(0.0, 0.0, 0.0),
    event_time = c(-3L, -2L, -1L)
  )
  expect_warning(
    res <- test_parallel_trends_joint(pe, alpha = 0.05),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_equal(res$n_pre_periods, 0L)
  expect_true(is.na(res$joint_stat))
})

test_that("TC-7.5.16: df2 lower bound — avg_n very small, df2 = 1", {
  pe <- make_pre_effects(
    att = c(1.0, 2.0),
    se = c(0.5, 0.5),
    n_treated = c(1L, 1L),
    n_control = c(1L, 1L)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  # avg_n = mean(c(1,1)) + mean(c(1,1)) = 1 + 1 = 2
  # df2 = max(as.integer(2 - 2), 1) = max(0, 1) = 1
  expect_equal(res$joint_df2, 1L)
})

test_that("TC-7.5.17: common timing input (cohort = NA_integer_)", {
  pe <- make_pre_effects(
    att = c(0.1, 0.2, 0.3),
    se = c(0.05, 0.05, 0.05),
    cohort = rep(NA_integer_, 3)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  expect_equal(res$n_pre_periods, 3L)
  expect_true(all(is.na(res$individual_tests$cohort)))
})


# ============================================================================
# Section 6: Integration Verification (E7-05.4)
# ============================================================================

test_that("E7-05.4.1: return list structure matches E7-06 parallel_trends_test field", {
  pe <- make_pre_effects(
    att = c(0.1, 0.2, 0.3),
    se = c(0.05, 0.05, 0.05)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")
  # E7-06 expects these 14 fields in parallel_trends_test
  expected_fields <- c(
    "joint_stat", "joint_pvalue", "joint_df1", "joint_df2",
    "reject_null", "n_pre_periods", "excluded_periods",
    "n_significant_05", "n_significant_10", "max_abs_pre_att",
    "individual_tests", "alpha", "test_type", "interpretation"
  )
  expect_equal(sort(names(res)), sort(expected_fields))
  # Type checks for downstream compatibility

  expect_true(is.numeric(res$joint_stat))
  expect_true(is.numeric(res$joint_pvalue))
  expect_true(is.integer(res$joint_df1))
  expect_true(is.integer(res$joint_df2))
  expect_true(is.logical(res$reject_null))
  expect_true(is.integer(res$n_pre_periods))
  expect_true(is.integer(res$n_significant_05))
  expect_true(is.integer(res$n_significant_10))
  expect_true(is.numeric(res$max_abs_pre_att))
  expect_true(is.data.frame(res$individual_tests))
  expect_true(is.numeric(res$alpha))
  expect_true(is.character(res$test_type))
  expect_true(is.character(res$interpretation))
})

test_that("E7-05.4.2: Common Timing and Staggered inputs both work", {
  # Common Timing: cohort = NA_integer_
  pe_ct <- make_pre_effects(
    att = c(0.1, 0.2), se = c(0.05, 0.05),
    cohort = rep(NA_integer_, 2)
  )
  res_ct <- test_parallel_trends_joint(pe_ct, test_type = "f")
  expect_equal(res_ct$n_pre_periods, 2L)

  # Staggered: cohort = actual values
  pe_stag <- make_pre_effects(
    att = c(0.1, 0.2, 0.3), se = c(0.05, 0.05, 0.05),
    cohort = c(5L, 5L, 7L)
  )
  res_stag <- test_parallel_trends_joint(pe_stag, test_type = "f")
  expect_equal(res_stag$n_pre_periods, 3L)
  expect_true(all(res_stag$individual_tests$cohort %in% c(5L, 7L)))
})

test_that("E7-05.4.3: K=0 path does not cause downstream errors", {
  # All anchors → K=0
  pe <- make_pre_effects(
    att = c(0.0, 0.0), se = c(0.0, 0.0),
    is_anchor = c(TRUE, TRUE)
  )
  expect_warning(
    res <- test_parallel_trends_joint(pe),
    class = "lwdid_insufficient_pre_periods"
  )
  # Downstream code should be able to access all fields without error
  expect_true(is.na(res$joint_stat))
  expect_true(is.na(res$joint_pvalue))
  expect_false(res$reject_null)
  expect_equal(nrow(res$individual_tests), 0L)
  # Can safely check individual_tests columns
  expect_true("significant" %in% names(res$individual_tests))
  expect_true("event_time" %in% names(res$individual_tests))
  # Can safely compute on max_abs_pre_att
  expect_true(is.na(res$max_abs_pre_att))
})
