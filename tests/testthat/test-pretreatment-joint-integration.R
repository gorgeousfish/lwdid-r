# ============================================================================
# test-pretreatment-joint-integration.R
# Story E7-05 Task 4: Integration verification for E7-06 compatibility
#
# Verifies test_parallel_trends_joint() return structure is compatible
# with Story E7-06 lwdid() integration expectations:
#   - result$parallel_trends_test field structure
#   - Common Timing and Staggered input paths
#   - K=0 path does not cause downstream errors
# ============================================================================

# ---- Helper: make pre-treatment effects data.frame ----
make_integration_pre_effects <- function(
    n = 5L,
    cohort_vals = rep(NA_integer_, 5),
    att_vals = rnorm(5, 0, 0.1),
    se_vals = rep(0.5, 5),
    include_anchor = TRUE) {
  t_stats <- att_vals / se_vals
  pvals <- 2 * pt(abs(t_stats), df = 48, lower.tail = FALSE)
  pe <- data.frame(
    cohort = cohort_vals,
    period = seq_len(n),
    event_time = seq(-n, -1L),
    att = att_vals,
    se = se_vals,
    t_stat = t_stats,
    pvalue = pvals,
    ci_lower = att_vals - 1.96 * se_vals,
    ci_upper = att_vals + 1.96 * se_vals,
    n_treated = rep(25L, n),
    n_control = rep(25L, n),
    is_anchor = rep(FALSE, n),
    rolling_window_size = as.integer(seq(n - 1L, 0L)),
    df_inference = rep(48L, n),
    stringsAsFactors = FALSE
  )
  if (include_anchor) {
    anchor <- data.frame(
      cohort = cohort_vals[1],
      period = n + 1L,
      event_time = -1L,
      att = 0.0, se = 0.0,
      t_stat = NaN, pvalue = NaN,
      ci_lower = 0.0, ci_upper = 0.0,
      n_treated = 25L, n_control = 25L,
      is_anchor = TRUE,
      rolling_window_size = 0L,
      df_inference = 48L,
      stringsAsFactors = FALSE
    )
    pe <- rbind(pe, anchor)
  }
  pe
}

# ============================================================
# E7-05.4.1: Return list structure matches E7-06 expectations
# ============================================================
test_that("E7-05.4.1: return list has all 14 fields for E7-06", {
  pe <- make_integration_pre_effects()
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")

  # E7-06 expects parallel_trends_test to be a list with these fields
  expected_fields <- c(
    "joint_stat", "joint_pvalue", "joint_df1", "joint_df2",
    "reject_null", "n_pre_periods", "excluded_periods",
    "n_significant_05", "n_significant_10", "max_abs_pre_att",
    "individual_tests", "alpha", "test_type", "interpretation"
  )
  expect_true(is.list(res))
  expect_equal(length(res), 14L)
  for (field in expected_fields) {
    expect_true(field %in% names(res),
                info = sprintf("Missing field: %s", field))
  }

  # Type checks for E7-06 downstream usage
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

  # individual_tests has expected 10 columns
  indiv_cols <- c("event_time", "cohort", "period", "att",
                  "se", "t_stat", "pvalue", "significant",
                  "n_treated", "n_control")
  expect_equal(sort(names(res$individual_tests)), sort(indiv_cols))

  # joint_pvalue in [0, 1]
  expect_true(res$joint_pvalue >= 0 && res$joint_pvalue <= 1)

  # Wald test also returns same structure
  res_w <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "wald")
  expect_equal(length(res_w), 14L)
  for (field in expected_fields) {
    expect_true(field %in% names(res_w),
                info = sprintf("Wald missing field: %s", field))
  }
  expect_equal(res_w$test_type, "wald")
  expect_equal(res_w$joint_df2, 0L)  # Wald has df2=0
})

# ============================================================
# E7-05.4.2: Common Timing and Staggered inputs both work
# ============================================================
test_that("E7-05.4.2: Common Timing input (cohort=NA) works", {
  pe <- make_integration_pre_effects(
    n = 4L,
    cohort_vals = rep(NA_integer_, 4),
    att_vals = c(0.02, -0.03, 0.01, -0.01),
    se_vals = rep(0.5, 4)
  )
  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")

  expect_equal(res$n_pre_periods, 4L)
  expect_true(all(is.na(res$individual_tests$cohort)))
  expect_true(is.numeric(res$joint_stat))
  expect_true(is.numeric(res$joint_pvalue))
})

test_that("E7-05.4.2: Staggered input (cohort=actual) works", {
  # Two cohorts: g=4 (3 periods) and g=6 (2 periods)
  pe <- data.frame(
    cohort = c(4L, 4L, 4L, 6L, 6L,
               4L, 6L),
    period = c(1L, 2L, 3L, 3L, 4L,
               3L, 5L),
    event_time = c(-3L, -2L, -1L, -3L, -2L,
                   -1L, -1L),
    att = c(0.1, -0.05, 0.0, 0.08, -0.02,
            0.0, 0.0),
    se = c(0.3, 0.4, 0.0, 0.5, 0.6,
           0.0, 0.0),
    t_stat = c(0.333, -0.125, NaN, 0.16, -0.033,
               NaN, NaN),
    pvalue = c(0.74, 0.90, NaN, 0.87, 0.97,
               NaN, NaN),
    ci_lower = c(-0.49, -0.83, 0.0, -0.90, -1.22,
                 0.0, 0.0),
    ci_upper = c(0.69, 0.73, 0.0, 1.06, 1.18,
                 0.0, 0.0),
    n_treated = c(10L, 10L, 10L, 5L, 5L,
                  10L, 5L),
    n_control = c(15L, 15L, 15L, 20L, 20L,
                  15L, 20L),
    is_anchor = c(FALSE, FALSE, FALSE, FALSE, FALSE,
                  TRUE, TRUE),
    rolling_window_size = c(2L, 1L, 0L, 2L, 1L,
                            0L, 0L),
    df_inference = rep(23L, 7),
    stringsAsFactors = FALSE
  )

  res <- test_parallel_trends_joint(pe, alpha = 0.05, test_type = "f")

  # Anchors excluded, 4 valid periods remain
  expect_equal(res$n_pre_periods, 4L)
  expect_true(4L %in% res$individual_tests$cohort)
  expect_true(6L %in% res$individual_tests$cohort)
  expect_true(is.numeric(res$joint_stat))
  expect_false(res$reject_null)  # small t-stats → don't reject
})

# ============================================================
# E7-05.4.3: K=0 path does not cause E7-06 downstream errors
# ============================================================
test_that("E7-05.4.3: K=0 returns valid structure for E7-06", {
  # All rows are anchors → K=0
  pe <- data.frame(
    cohort = c(4L, 6L),
    period = c(3L, 5L),
    event_time = c(-1L, -1L),
    att = c(0.0, 0.0),
    se = c(0.0, 0.0),
    t_stat = c(NaN, NaN),
    pvalue = c(NaN, NaN),
    ci_lower = c(0.0, 0.0),
    ci_upper = c(0.0, 0.0),
    n_treated = c(10L, 5L),
    n_control = c(15L, 20L),
    is_anchor = c(TRUE, TRUE),
    rolling_window_size = c(0L, 0L),
    df_inference = c(23L, 23L),
    stringsAsFactors = FALSE
  )

  # K=0 emits two warnings (insufficient + no valid periods)
  expect_warning(
    test_parallel_trends_joint(pe, alpha = 0.05),
    class = "lwdid_insufficient_pre_periods"
  )
  res <- suppressWarnings(
    test_parallel_trends_joint(pe, alpha = 0.05)
  )

  # Still returns full 14-field list
  expect_true(is.list(res))
  expect_equal(length(res), 14L)

  # NA values for statistics
  expect_true(is.na(res$joint_stat))
  expect_true(is.na(res$joint_pvalue))
  expect_false(res$reject_null)
  expect_equal(res$n_pre_periods, 0L)

  # individual_tests is empty but has correct columns
  expect_true(is.data.frame(res$individual_tests))
  expect_equal(nrow(res$individual_tests), 0L)
  expect_equal(ncol(res$individual_tests), 10L)

  # E7-06 downstream: accessing fields should not error
  expect_true(is.character(res$interpretation))
  expect_true(is.integer(res$n_significant_05))
  expect_true(is.integer(res$n_significant_10))
  expect_true(is.na(res$max_abs_pre_att))

  # Simulate E7-06 downstream access patterns
  # These should all work without error
  if (!is.na(res$joint_pvalue)) {
    reject <- res$joint_pvalue < 0.05
  } else {
    reject <- FALSE
  }
  expect_false(reject)

  # Accessing individual_tests columns on empty df
  expect_equal(length(res$individual_tests$att), 0L)
  expect_equal(length(res$individual_tests$significant), 0L)
})
