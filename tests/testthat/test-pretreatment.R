# ===========================================================================
# test-pretreatment.R
# L1 Unit tests for pretreatment internal logic & Anchor properties
# (Story E7-07, Task E7-07.2)
# ===========================================================================

# Internal function accessors
.pre_ct <- lwdid:::estimate_pre_treatment_common
.pre_stag <- lwdid:::estimate_pre_treatment_staggered
.pre_joint <- lwdid:::test_parallel_trends_joint

# ===================================================================
# TC-7.7.6: Symmetric transform reference periods correct
# CT mode (S=5), period t=2: ref periods = {3,4}
# ===================================================================
test_that("TC-7.7.6: symmetric transform uses correct reference periods", {
  set.seed(42L)
  n_units <- 20L
  tpost1 <- 5L
  dt <- data.table::CJ(unit = seq_len(n_units), time = 1:8)
  dt[, d := as.integer(unit <= 8L)]
  alpha_vals <- seq(1, by = 0.5, length.out = n_units)
  dt[, y := alpha_vals[unit] + 0.5 * time]

  result <- suppressWarnings(suppressMessages(
    .pre_ct(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = tpost1, rolling = "demean"
    )
  ))

  row_t2 <- result[result$event_time == -3L &
                    !result$is_anchor, ]
  expect_equal(nrow(row_t2), 1L)

  # Manual: y_trans_i = Y_i2 - mean(Y_i3, Y_i4)
  # Y_it = alpha_i + 0.5*t (deterministic, no noise)
  # y_trans_i = (alpha_i+1) - (alpha_i+1.75) = -0.75
  # Same for all units => ATT = 0.0
  treated_ids <- which(seq_len(n_units) <= 8L)
  control_ids <- which(seq_len(n_units) > 8L)
  y_trans <- numeric(n_units)
  for (i in seq_len(n_units)) {
    yi2 <- alpha_vals[i] + 0.5 * 2
    yi3 <- alpha_vals[i] + 0.5 * 3
    yi4 <- alpha_vals[i] + 0.5 * 4
    y_trans[i] <- yi2 - mean(c(yi3, yi4))
  }
  manual_att <- mean(y_trans[treated_ids]) -
    mean(y_trans[control_ids])


  expect_equal(row_t2$att, manual_att, tolerance = 1e-10)
  expect_equal(row_t2$att, 0.0, tolerance = 1e-10)
  expect_equal(row_t2$rolling_window_size, 2L)
})

# ===================================================================
# TC-7.7.7: Anchor point ATT=0 and SE=0 (staggered)
# Cohort g=4, anchor period t=3 (event_time=-1)
# ===================================================================
test_that("TC-7.7.7: anchor point ATT=0 and SE=0 in staggered mode", {
  dt <- generate_staggered_panel(
    n_per_cohort = 5L, cohorts = c(4L, 6L, 8L),
    n_never_treated = 5L, T_total = 10L, seed = 42L
  )
  result <- suppressWarnings(suppressMessages(
    .pre_stag(
      dt, y = "y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean"
    )
  ))
  anchor_g4 <- result[result$cohort == 4L &
                       result$is_anchor, ]
  expect_equal(nrow(anchor_g4), 1L)
  expect_identical(anchor_g4$att, 0.0)
  expect_identical(anchor_g4$se, 0.0)
  expect_equal(anchor_g4$period, 3L)
  expect_equal(anchor_g4$event_time, -1L)
})


# ===================================================================
# TC-7.7.8: Control group based on t selection (FATAL-004)
# not_yet_treated: control = G_i > t (not G_i > ref period)
# ===================================================================
test_that("TC-7.7.8: control group uses G_i > t for not_yet_treated", {
  dt <- generate_staggered_panel(
    n_per_cohort = 5L, cohorts = c(4L, 6L, 8L),
    n_never_treated = 5L, T_total = 10L, seed = 42L
  )
  result <- suppressWarnings(suppressMessages(
    .pre_stag(
      dt, y = "y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated"
    )
  ))
  row_g4_t2 <- result[result$cohort == 4L &
                       result$period == 2L, ]
  expect_equal(nrow(row_g4_t2), 1L)
  expect_equal(row_g4_t2$n_treated, 5L)
  expect_equal(row_g4_t2$n_control, 15L)
})


# ===================================================================
# TC-7.7.9: event_time calculation correct
# CT: event_time = period - tpost1
# Staggered: event_time = period - cohort
# ===================================================================
test_that("TC-7.7.9: event_time computed correctly for CT and staggered", {
  dt_ct <- generate_ct_panel(
    N = 40L, T_total = 8L, tpost1 = 5L, seed = 42L
  )
  result_ct <- suppressWarnings(suppressMessages(
    .pre_ct(
      dt_ct, y = "y", ivar = "id", d = "d",
      tvar = "time", tpost1 = 5L
    )
  ))
  expect_true(all(result_ct$event_time < 0L))
  anchor_ct <- result_ct[result_ct$is_anchor, ]
  expect_equal(anchor_ct$event_time, -1L)
  row_t2 <- result_ct[result_ct$period == 2L, ]
  if (nrow(row_t2) > 0L) {
    expect_equal(row_t2$event_time, -3L)
  }

  dt_stag <- generate_staggered_panel(
    n_per_cohort = 5L, cohorts = c(4L, 6L, 8L),
    n_never_treated = 5L, T_total = 10L, seed = 42L
  )
  result_stag <- suppressWarnings(suppressMessages(
    .pre_stag(
      dt_stag, y = "y", ivar = "id",
      tvar = "time", gvar = "gvar"
    )
  ))
  row_g6_t3 <- result_stag[result_stag$cohort == 6L &
                            result_stag$period == 3L, ]
  if (nrow(row_g6_t3) > 0L) {
    expect_equal(row_g6_t3$event_time, -3L)
  }
  non_anchor <- result_stag[!result_stag$is_anchor, ]
  expect_true(all(non_anchor$event_time < -1L, na.rm = TRUE))
})

# ===================================================================
# TC-7.7.10: Anchor is_anchor=TRUE and ATT=0.0 exact (both modes)
# ===================================================================
test_that("TC-7.7.10: anchor is_anchor=TRUE and ATT=0.0 exact", {
  dt_ct <- generate_ct_panel(
    N = 40L, T_total = 8L, tpost1 = 5L, seed = 42L
  )
  result_ct <- suppressWarnings(suppressMessages(
    .pre_ct(
      dt_ct, y = "y", ivar = "id", d = "d",
      tvar = "time", tpost1 = 5L
    )
  ))
  anchor_ct <- result_ct[result_ct$is_anchor, ]
  expect_equal(nrow(anchor_ct), 1L)
  expect_true(anchor_ct$is_anchor)
  expect_identical(anchor_ct$att, 0.0)
  expect_identical(anchor_ct$se, 0.0)

  dt_stag <- generate_staggered_panel(
    n_per_cohort = 5L, cohorts = c(4L, 6L, 8L),
    n_never_treated = 5L, T_total = 10L, seed = 42L
  )
  result_stag <- suppressWarnings(suppressMessages(
    .pre_stag(
      dt_stag, y = "y", ivar = "id",
      tvar = "time", gvar = "gvar"
    )
  ))
  anchors_stag <- result_stag[result_stag$is_anchor, ]
  expect_equal(nrow(anchors_stag), 3L)
  for (i in seq_len(nrow(anchors_stag))) {
    expect_identical(anchors_stag$att[i], 0.0)
    expect_identical(anchors_stag$se[i], 0.0)
    expect_true(anchors_stag$is_anchor[i])
  }
})


# ===================================================================
# TC-7.7.11: Anchor excluded from joint test
# ===================================================================
test_that("TC-7.7.11: anchor period excluded from joint test", {
  dt_ct <- generate_ct_panel(
    N = 40L, T_total = 8L, tpost1 = 5L, seed = 42L
  )
  pre_effects <- suppressWarnings(suppressMessages(
    .pre_ct(
      dt_ct, y = "y", ivar = "id", d = "d",
      tvar = "time", tpost1 = 5L
    )
  ))
  joint_result <- suppressWarnings(suppressMessages(
    .pre_joint(pre_effects)
  ))
  expect_true(-1L %in% joint_result$excluded_periods)
  if (nrow(joint_result$individual_tests) > 0L) {
    expect_false(
      any(joint_result$individual_tests$event_time == -1L)
    )
  }
})


# ===========================================================================
# Task E7-07.4: L2 Numerical tests (pretreatment portion)
# TC-7.7.13, TC-7.7.14, TC-7.7.15
# ===========================================================================

# TC-7.7.13: Pre-treatment ATT near zero under parallel trends
# Generate staggered panel with NO treatment effect in pre-period (pt_violation=0).
# All non-anchor ATTs should be close to 0 (within reasonable tolerance given noise).
test_that("TC-7.7.13: pre-period ATT near zero under parallel trends", {
  dt <- generate_staggered_panel(
    n_per_cohort = 10L, cohorts = c(4L, 6L, 8L),
    n_never_treated = 10L, T_total = 10L,
    tau_base = 2.0, tau_dynamic = 0.5,
    pt_violation = 0, seed = 42L
  )

  result <- suppressWarnings(suppressMessages(
    .pre_stag(dt, y = "y", ivar = "id", tvar = "time",
              gvar = "gvar", rolling = "demean")
  ))

  # Under parallel trends (pt_violation=0), pre-period ATTs should be
  # close to 0 (only noise). With N=10 per cohort, |ATT| < 2.0 is
  # a generous bound given noise sd=0.5.
  non_anchor <- result[!result$is_anchor, ]
  valid_rows <- non_anchor[!is.na(non_anchor$att), ]

  expect_true(nrow(valid_rows) > 0L)
  for (i in seq_len(nrow(valid_rows))) {
    expect_true(abs(valid_rows$att[i]) < 2.0,
      info = sprintf("cohort=%d, period=%d, att=%.4f too large",
                     valid_rows$cohort[i], valid_rows$period[i],
                     valid_rows$att[i]))
  }

  # All ATTs should be finite
  expect_true(all(is.finite(valid_rows$att)))
})

# TC-7.7.14: Joint F-statistic manual calculation
# t_stats = c(1.5, -0.8, 2.1)
# F = mean(t^2) = (2.25 + 0.64 + 4.41) / 3 = 7.3 / 3 = 2.4333...
test_that("TC-7.7.14: joint F-statistic matches manual calculation", {
  mock_pre <- data.frame(
    cohort = rep(NA_integer_, 3),
    period = c(1L, 2L, 3L),
    event_time = c(-4L, -3L, -2L),
    att = c(0.3, -0.16, 0.42),
    se = c(0.2, 0.2, 0.2),
    t_stat = c(1.5, -0.8, 2.1),
    pvalue = c(0.13, 0.42, 0.04),
    n_treated = rep(20L, 3),
    n_control = rep(30L, 3),
    is_anchor = c(FALSE, FALSE, FALSE),
    rolling_window_size = c(2L, 1L, 0L),
    df_inference = rep(48L, 3),
    ci_lower = c(0.3 - 1.96 * 0.2, -0.16 - 1.96 * 0.2, 0.42 - 1.96 * 0.2),
    ci_upper = c(0.3 + 1.96 * 0.2, -0.16 + 1.96 * 0.2, 0.42 + 1.96 * 0.2),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(suppressMessages(
    .pre_joint(mock_pre, test_type = "f")
  ))

  # Manual: F = (1.5^2 + 0.8^2 + 2.1^2) / 3 = (2.25 + 0.64 + 4.41) / 3
  #         = 7.3 / 3 = 2.43333...
  expect_equal(res$joint_stat, 2.43333333333333, tolerance = 1e-10)
  expect_equal(res$joint_df1, 3L)
  expect_equal(res$n_pre_periods, 3L)

  # df2 = max(as.integer(mean(20,20,20) + mean(30,30,30) - 2), 1) = 48
  expect_equal(res$joint_df2, 48L)

  # p-value from F distribution
  expected_p <- stats::pf(7.3 / 3, 3, 48, lower.tail = FALSE)
  expect_equal(res$joint_pvalue, expected_p, tolerance = 1e-12)
})

# TC-7.7.15: Wald statistic manual calculation
# t_stats = c(1.5, -0.8, 2.1)
# W = sum(t^2) = 2.25 + 0.64 + 4.41 = 7.3
test_that("TC-7.7.15: Wald statistic matches manual calculation", {
  mock_pre <- data.frame(
    cohort = rep(NA_integer_, 3),
    period = c(1L, 2L, 3L),
    event_time = c(-4L, -3L, -2L),
    att = c(0.3, -0.16, 0.42),
    se = c(0.2, 0.2, 0.2),
    t_stat = c(1.5, -0.8, 2.1),
    pvalue = c(0.13, 0.42, 0.04),
    n_treated = rep(20L, 3),
    n_control = rep(30L, 3),
    is_anchor = c(FALSE, FALSE, FALSE),
    rolling_window_size = c(2L, 1L, 0L),
    df_inference = rep(48L, 3),
    ci_lower = c(0.3 - 1.96 * 0.2, -0.16 - 1.96 * 0.2, 0.42 - 1.96 * 0.2),
    ci_upper = c(0.3 + 1.96 * 0.2, -0.16 + 1.96 * 0.2, 0.42 + 1.96 * 0.2),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(suppressMessages(
    .pre_joint(mock_pre, test_type = "wald")
  ))

  # Manual: W = 1.5^2 + 0.8^2 + 2.1^2 = 2.25 + 0.64 + 4.41 = 7.3
  expect_equal(res$joint_stat, 7.3, tolerance = 1e-10)
  expect_equal(res$joint_df1, 3L)
  expect_equal(res$joint_df2, 0L)

  # p-value from chi-squared distribution
  expected_p <- stats::pchisq(7.3, df = 3, lower.tail = FALSE)
  expect_equal(res$joint_pvalue, expected_p, tolerance = 1e-12)
})


# ===========================================================================
# Task E7-07.5: L3 Monte Carlo tests (pretreatment portion)
# TC-7.7.19 and TC-7.7.20
# ===========================================================================

# TC-7.7.19: When parallel trends hold (pt_violation=0),
# joint F-test rejection rate should be approximately alpha=0.05.
# 200 simulations. Acceptable range: [0.02, 0.08] (alpha ± 3%).
test_that("TC-7.7.19: PT holds => rejection rate approx alpha", {
  skip_on_cran()
  n_sims <- 200L
  alpha <- 0.05
  rejections <- logical(n_sims)

  for (sim in seq_len(n_sims)) {
    dt <- generate_ct_panel(N = 100L, T_total = 8L, tpost1 = 5L,
                            tau = 2.0, pt_violation = 0, seed = sim)
    pre_eff <- suppressWarnings(suppressMessages(
      .pre_ct(dt, y = "y", ivar = "id", d = "d",
              tvar = "time", tpost1 = 5L, rolling = "demean")
    ))
    joint_res <- suppressWarnings(suppressMessages(
      .pre_joint(pre_eff, alpha = alpha, test_type = "f")
    ))
    rejections[sim] <- isTRUE(joint_res$reject_null)
  }

  rejection_rate <- mean(rejections)
  expect_true(rejection_rate >= 0.02 && rejection_rate <= 0.08,
    info = sprintf("Rejection rate=%.3f, expected in [0.02, 0.08]",
                   rejection_rate))
})

# TC-7.7.20: When parallel trends are violated (pt_violation=2.0),
# joint F-test power > 0.8.
# 200 simulations with strong PT violation.
test_that("TC-7.7.20: PT violated => power > 0.8", {
  skip_on_cran()
  n_sims <- 200L
  alpha <- 0.05
  rejections <- logical(n_sims)

  for (sim in seq_len(n_sims)) {
    dt <- generate_ct_panel(N = 100L, T_total = 8L, tpost1 = 5L,
                            tau = 2.0, pt_violation = 2.0, seed = sim)
    pre_eff <- suppressWarnings(suppressMessages(
      .pre_ct(dt, y = "y", ivar = "id", d = "d",
              tvar = "time", tpost1 = 5L, rolling = "demean")
    ))
    joint_res <- suppressWarnings(suppressMessages(
      .pre_joint(pre_eff, alpha = alpha, test_type = "f")
    ))
    rejections[sim] <- isTRUE(joint_res$reject_null)
  }

  rejection_rate <- mean(rejections)
  expect_true(rejection_rate > 0.8,
    info = sprintf("Power=%.3f, expected > 0.8", rejection_rate))
})


# ===========================================================================
# Task E7-07.6: L4 Boundary condition tests (pretreatment portion)
# TC-7.7.23, TC-7.7.24, TC-7.7.25, TC-7.7.26, TC-7.7.27,
# TC-7.7.40, TC-7.7.41, TC-7.7.51
# ===========================================================================

# TC-7.7.23: Single pre-period (K=1) — F-test degenerates to t²
test_that("TC-7.7.23: K=1 F-test degenerates to t_stat^2", {
  pe <- data.frame(
    event_time = c(-2L, -1L),
    cohort = rep(NA_integer_, 2),
    period = c(3L, 4L),
    att = c(0.5, 0.0),
    se = c(0.25, 0.0),
    t_stat = c(2.0, NA_real_),
    pvalue = c(0.05, NA_real_),
    is_anchor = c(FALSE, TRUE),
    n_treated = c(20L, 20L),
    n_control = c(20L, 20L),
    rolling_window_size = c(1L, 0L),
    df_inference = c(38L, NA_integer_),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(.pre_joint(pe, alpha = 0.05, test_type = "f",
                                      min_pre_periods = 1L))
  # K=1: F = t^2 / 1 = t^2 = 4.0
  expect_equal(res$joint_stat, 2.0^2, tolerance = 1e-12)
  expect_equal(res$n_pre_periods, 1L)
  expect_equal(res$joint_df1, 1L)
})

# TC-7.7.24: All pre-periods NA → K=0 return
test_that("TC-7.7.24: all pre-periods NA returns K=0 result", {
  pe <- data.frame(
    event_time = c(-3L, -2L, -1L),
    cohort = rep(NA_integer_, 3),
    period = c(2L, 3L, 4L),
    att = c(NA_real_, NA_real_, 0.0),
    se = c(NA_real_, NA_real_, 0.0),
    t_stat = c(NA_real_, NA_real_, NA_real_),
    pvalue = c(NA_real_, NA_real_, NA_real_),
    is_anchor = c(FALSE, FALSE, TRUE),
    n_treated = c(5L, 5L, 5L),
    n_control = c(5L, 5L, 5L),
    rolling_window_size = c(2L, 1L, 0L),
    df_inference = c(NA_integer_, NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(.pre_joint(pe, alpha = 0.05))
  expect_true(is.na(res$joint_stat))
  expect_equal(res$n_pre_periods, 0L)
  expect_false(res$reject_null)
})

# TC-7.7.25: exclude_pre_periods = S-2 (only 1 pre-period left)
test_that("TC-7.7.25: exclude_pre_periods=3 leaves minimal periods", {
  dt <- generate_ct_panel(N = 40L, T_total = 8L, tpost1 = 5L, seed = 42L)

  # exclude_pre_periods=3: effective_tpost1 = 5-3 = 2
  # pre_periods = {1} (only period 1 < 2)
  # anchor = period 1 (effective_tpost1 - 1 = 1)
  # estimable = {} (after removing anchor and T_min)
  # → should return NULL or very few rows
  result <- suppressWarnings(suppressMessages(
    .pre_ct(dt, y = "y", ivar = "id", d = "d",
            tvar = "time", tpost1 = 5L, rolling = "demean",
            exclude_pre_periods = 3L)
  ))

  # With exclude=3, very few or no estimable periods remain
  if (!is.null(result)) {
    # At most 1 non-anchor row + 1 anchor row
    non_anchor <- result[!result$is_anchor, ]
    expect_true(nrow(non_anchor) <= 1L)
  }
})

# TC-7.7.26: Extreme value outcome does not overflow
test_that("TC-7.7.26: extreme Y ~ 1e10 produces finite ATT", {
  dt <- generate_extreme_panel(N = 20L, T_total = 6L, seed = 42L)

  result <- suppressWarnings(suppressMessages(
    .pre_ct(dt, y = "y", ivar = "id", d = "d",
            tvar = "time", tpost1 = 4L, rolling = "demean")
  ))

  expect_false(is.null(result))
  valid_rows <- result[!result$is_anchor & !is.na(result$att), ]
  if (nrow(valid_rows) > 0L) {
    expect_true(all(is.finite(valid_rows$att)),
      info = "ATT should be finite even with Y ~ 1e10")
  }
})

# TC-7.7.27: Constant outcome Y=5.0 → pre-period ATT all 0.0
test_that("TC-7.7.27: constant Y=5.0 produces ATT=0 after demean", {
  dt <- generate_constant_panel(N = 20L, T_total = 6L, seed = 42L)

  result <- suppressWarnings(suppressMessages(
    .pre_ct(dt, y = "y", ivar = "id", d = "d",
            tvar = "time", tpost1 = 4L, rolling = "demean")
  ))

  expect_false(is.null(result))
  non_anchor <- result[!result$is_anchor & !is.na(result$att), ]
  for (i in seq_len(nrow(non_anchor))) {
    expect_equal(non_anchor$att[i], 0.0, tolerance = 1e-10,
      info = sprintf("period=%d ATT should be 0 for constant Y",
                     non_anchor$period[i]))
  }
})

# TC-7.7.40: K=0 empty data.frame has all 10 columns
test_that("TC-7.7.40: K=0 individual_tests has 10 columns", {
  pe <- data.frame(
    event_time = c(-1L),
    cohort = NA_integer_,
    period = 4L,
    att = 0.0,
    se = 0.0,
    t_stat = NA_real_,
    pvalue = NA_real_,
    is_anchor = TRUE,
    n_treated = 10L,
    n_control = 10L,
    rolling_window_size = 0L,
    df_inference = NA_integer_,
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(.pre_joint(pe, alpha = 0.05))
  expect_equal(res$n_pre_periods, 0L)

  indiv <- res$individual_tests
  expect_equal(nrow(indiv), 0L)
  expected_cols <- c("event_time", "cohort", "period", "att",
                     "se", "t_stat", "pvalue", "significant",
                     "n_treated", "n_control")
  expect_equal(sort(names(indiv)), sort(expected_cols))
})

# TC-7.7.41: min_obs < 3 returns NA for that period
test_that("TC-7.7.41: period with < 3 obs handled gracefully", {
  # Create tiny panel: 1 treated + 1 control, 4 periods, tpost1=3
  set.seed(42L)
  dt <- data.table::data.table(
    id = rep(1:2, each = 4L),
    time = rep(1:4, 2L),
    y = rnorm(8),
    d = rep(c(1L, 0L), each = 4L)
  )

  # With only 2 units total (1 treated + 1 control),
  # each period has n_treated + n_control = 2 < 3
  # Function should not crash — either returns NULL or data with NA ATTs
  result <- suppressWarnings(suppressMessages(
    .pre_ct(dt, y = "y", ivar = "id", d = "d",
            tvar = "time", tpost1 = 3L, rolling = "demean")
  ))

  # Key assertion: function completes without error
  # Result is either NULL (no estimable periods) or a data.frame
  expect_true(is.null(result) || is.data.frame(result))

  # If result exists, verify structure is valid
  if (!is.null(result) && nrow(result) > 0L) {
    expect_true("att" %in% names(result))
    expect_true("is_anchor" %in% names(result))
  }
})

# TC-7.7.51: K=0 return list field completeness
test_that("TC-7.7.51: K=0 return has all required fields", {
  pe <- data.frame(
    event_time = c(-1L),
    cohort = NA_integer_,
    period = 4L,
    att = 0.0,
    se = 0.0,
    t_stat = NA_real_,
    pvalue = NA_real_,
    is_anchor = TRUE,
    n_treated = 10L,
    n_control = 10L,
    rolling_window_size = 0L,
    df_inference = NA_integer_,
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(.pre_joint(pe, alpha = 0.05))

  # Verify all required fields exist
  expect_equal(res$n_significant_05, 0L)
  expect_equal(res$n_significant_10, 0L)
  expect_true(is.na(res$max_abs_pre_att))
  expect_false(is.null(res$interpretation))
  expect_true(is.character(res$interpretation))
  expect_true(nchar(res$interpretation) > 0L)
  expect_false(res$reject_null)
  expect_equal(res$n_pre_periods, 0L)
})
