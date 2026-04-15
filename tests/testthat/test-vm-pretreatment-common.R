# test-vm-pretreatment-common.R
# Vibe-math end-to-end numerical verification for E7-03 Tasks 4 & 5
# Tests VM-01 through VM-08 + boundary cases

# ============================================================
# Task 4.1: VM-01 demean symmetric transform
# ============================================================
test_that("VM-01: demean symmetric transform matches vibe-math (Task 4.1)", {
  # 3 units x 6 periods, S=4 (tpost1=4)
  # Unit 1 (treated): Y = 10,12,15,20,25,30
  # Unit 2 (control): Y = 8,10,13,16,19,22
  # Unit 3 (control): Y = 5,7,9,12,15,18
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # Pre-periods: 1,2,3. Anchor=3, T_min=1 excluded. Estimable: {2}
  # t=2, ref={3}: demean = Y_i2 - Y_i3
  # Unit1: 12-15=-3, Unit2: 10-13=-3, Unit3: 7-9=-2
  # ATT = -3 - (-3+-2)/2 = -3 - (-2.5) = -0.5
  expect_true(!is.null(result))
  expect_true(is.data.frame(result))

  anchor_row <- result[result$is_anchor == TRUE, ]
  expect_equal(nrow(anchor_row), 1L)
  expect_equal(anchor_row$att, 0.0)
  expect_equal(anchor_row$event_time, -1L)

  est_row <- result[result$event_time == -2L, ]
  expect_equal(nrow(est_row), 1L)
  expect_equal(est_row$att, -0.5, tolerance = 1e-10)
  expect_equal(est_row$rolling_window_size, 1L)
})


# ============================================================
# Task 4.2: VM-02 detrend symmetric transform
# ============================================================
test_that("VM-02: detrend symmetric transform matches vibe-math (Task 4.2)", {
  # Same panel but S=5 (tpost1=5) so t=2 has ref={3,4} (2 periods for OLS)
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 5L, rolling = "detrend", estimator = "ra", alpha = 0.05
  )

  # Pre-periods: 1,2,3,4. Anchor=4, T_min=1 excluded. Estimable: {2,3}
  # t=2, ref={3,4}: OLS Y_iq = alpha + beta*q
  #   Unit1: (3,15),(4,20) -> beta=5, alpha=0 -> pred(2)=10 -> detrend=12-10=2
  #   Unit2: (3,13),(4,16) -> beta=3, alpha=4 -> pred(2)=10 -> detrend=10-10=0
  #   Unit3: (3,9),(4,12) -> beta=3, alpha=0 -> pred(2)=6 -> detrend=7-6=1
  #   ATT = 2 - (0+1)/2 = 1.5

  expect_true(!is.null(result))
  est_t2 <- result[result$event_time == (2L - 5L), ]  # event_time = -3
  expect_equal(nrow(est_t2), 1L)
  expect_equal(est_t2$att, 1.5, tolerance = 1e-10)

  # t=3, ref={4}: only 1 ref period -> degrades to demean
  # Unit1: 15-20=-5, Unit2: 13-16=-3, Unit3: 9-12=-3
  # ATT = -5 - (-3+-3)/2 = -5 - (-3) = -2
  est_t3 <- result[result$event_time == (3L - 5L), ]  # event_time = -2
  expect_equal(nrow(est_t3), 1L)
  expect_equal(est_t3$att, -2.0, tolerance = 1e-10)
})

# ============================================================
# Task 4.3: VM-03 no-controls RA estimate
# ============================================================
test_that("VM-03: no-controls RA = simple diff of means (Task 4.3)", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # ATT at t=2 (demean, no controls) = mean(treated_trans) - mean(ctrl_trans)
  # = -3 - (-2.5) = -0.5
  est_row <- result[result$event_time == -2L, ]
  expect_equal(est_row$att, -0.5, tolerance = 1e-10)
})

# ============================================================
# Task 4.4: VM-04 rolling_window_size
# ============================================================
test_that("VM-04: rolling_window_size = S-1-t (Task 4.4)", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  # S=5 -> estimable: t=2 (rws=5-1-2=2), t=3 (rws=5-1-3=1)
  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 5L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  est_t2 <- result[result$event_time == (2L - 5L), ]
  est_t3 <- result[result$event_time == (3L - 5L), ]
  anchor <- result[result$is_anchor == TRUE, ]

  expect_equal(est_t2$rolling_window_size, 2L)  # |{3,4}| = 2
  expect_equal(est_t3$rolling_window_size, 1L)  # |{4}| = 1
  expect_equal(anchor$rolling_window_size, 0L)   # anchor = 0
})


# ============================================================
# Task 4.5: VM-05 symmetric vs standard transform difference
# ============================================================
test_that("VM-05: symmetric vs standard transform produce different ATT (Task 4.5)", {
  # With S=4, t=1 (if we could estimate it):
  # Symmetric ref={2,3}: Y_i1 - mean(Y_i2, Y_i3)
  #   Unit1: 10-(12+15)/2 = -3.5
  #   Unit2: 8-(10+13)/2 = -3.5
  #   Unit3: 5-(7+9)/2 = -3.0
  #   ATT_sym = -3.5 - (-3.5+-3.0)/2 = -0.25
  #
  # Standard ref={1,2,3}: Y_i1 - mean(Y_i1, Y_i2, Y_i3)
  #   Unit1: 10-(10+12+15)/3 = -2.333
  #   Unit2: 8-(8+10+13)/3 = -2.333
  #   Unit3: 5-(5+7+9)/3 = -2.0
  #   ATT_std = -2.333 - (-2.333+-2.0)/2 = -0.1667
  #
  # Difference = |-0.25 - (-0.1667)| = 0.0833 > 1e-6

  # We verify this by manually computing the symmetric transform
  # (the function uses symmetric, so we just verify the function output
  # and compare with the standard transform computed by hand)

  # Use S=5 so t=2 has ref={3,4} (symmetric) vs ref={2,3,4} (standard)
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 5L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # Symmetric at t=2, ref={3,4}: Y_i2 - mean(Y_i3, Y_i4)
  # Unit1: 12-(15+20)/2 = -5.5
  # Unit2: 10-(13+16)/2 = -4.5
  # Unit3: 7-(9+12)/2 = -3.5
  # ATT_sym = -5.5 - (-4.5+-3.5)/2 = -5.5 - (-4.0) = -1.5
  est_t2 <- result[result$event_time == (2L - 5L), ]
  att_symmetric <- est_t2$att

  # Standard at t=2, ref={2,3,4}: Y_i2 - mean(Y_i2, Y_i3, Y_i4)
  # Unit1: 12-(12+15+20)/3 = 12-15.667 = -3.667
  # Unit2: 10-(10+13+16)/3 = 10-13 = -3.0
  # Unit3: 7-(7+9+12)/3 = 7-9.333 = -2.333
  # ATT_std = -3.667 - (-3.0+-2.333)/2 = -3.667 - (-2.667) = -1.0
  att_standard <- -1.0  # hand-computed

  expect_true(abs(att_symmetric - att_standard) > 1e-6)
})

# ============================================================
# Task 4.6: VM-06 detrend degradation (ref=1 period)
# ============================================================
test_that("VM-06: detrend degrades to demean when ref=1 (Task 4.6)", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  # S=4: t=2 has ref={3} (1 period) -> detrend degrades to demean
  result_detrend <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "detrend", estimator = "ra", alpha = 0.05
  )
  result_demean <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # Both should produce identical ATT at t=2
  att_detrend <- result_detrend[result_detrend$event_time == -2L, "att"]
  att_demean <- result_demean[result_demean$event_time == -2L, "att"]
  expect_equal(att_detrend, att_demean, tolerance = 1e-10)
})

# ============================================================
# Task 4.7: VM-08 triple-difference interpretation
# ============================================================
test_that("VM-08: detrend with 1 ref = first-difference DiD (Task 4.7)", {
  # When detrend degrades to demean (ref=1 period), the pre-treatment
  # ATT equals a difference-in-first-differences:
  # ATT_t = (Y_1t - Y_1,t+1)_treat - (Y_0t - Y_0,t+1)_ctrl
  # This is the triple-difference interpretation (lw2025 eq 5.7)

  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  # S=4, t=2, ref={3}: detrend degrades to demean = Y_i2 - Y_i3
  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "detrend", estimator = "ra", alpha = 0.05
  )

  att_func <- result[result$event_time == -2L, "att"]

  # Manual triple-diff: (Y_12-Y_13)_treat - mean((Y_22-Y_23)_ctrl, (Y_32-Y_33)_ctrl)
  # = (12-15) - ((10-13)+(7-9))/2 = -3 - (-3+-2)/2 = -3-(-2.5) = -0.5
  att_manual <- (12 - 15) - ((10 - 13) + (7 - 9)) / 2
  expect_equal(att_func, att_manual, tolerance = 1e-10)
})


# ============================================================
# Task 5: Boundary cases and regression tests
# ============================================================

# Task 5.1: Only 1 pre-treatment period (T_min = S-1, that period is anchor)
test_that("Task 5.1: single pre-period returns anchor only (not NULL)", {
  # S=2, periods={1,2,...}. Pre-periods={1}. Anchor=1. T_min=1 excluded.
  # But anchor=1 should still be returned.
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    y = c(10, 20, 30, 8, 16, 24, 5, 12, 18),
    d = rep(c(1L, 0L, 0L), each = 3)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 2L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # Should return non-NULL with anchor row
  expect_true(!is.null(result))
  expect_true(is.data.frame(result))
  expect_true(any(result$is_anchor == TRUE))
  anchor <- result[result$is_anchor == TRUE, ]
  expect_equal(anchor$att, 0.0)
  expect_equal(anchor$event_time, -1L)
})

# Task 5.2: No pre-treatment periods -> NULL + warning
test_that("Task 5.2: no pre-treatment periods returns NULL with warning", {
  # All periods >= tpost1
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(5:7, 3),
    y = c(10, 20, 30, 8, 16, 24, 5, 12, 18),
    d = rep(c(1L, 0L, 0L), each = 3)
  )

  expect_warning(
    result <- estimate_pre_treatment_common(
      data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
      tpost1 = 5L, rolling = "demean", estimator = "ra", alpha = 0.05
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_null(result)
})

# Task 5.3: exclude_pre_periods parameter
test_that("Task 5.3: exclude_pre_periods reduces pre-period set", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 8),
    time = rep(1:8, 3),
    y = c(10, 12, 15, 18, 20, 25, 30, 35,
          8, 10, 13, 15, 16, 19, 22, 25,
          5, 7, 9, 11, 12, 15, 18, 21),
    d = rep(c(1L, 0L, 0L), each = 8)
  )

  # S=6, exclude_pre_periods=2 -> effective_tpost1=4
  # Pre-periods < 4: {1,2,3}. Anchor=3. T_min=1 excluded. Estimable={2}
  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 6L, rolling = "demean", estimator = "ra", alpha = 0.05,
    exclude_pre_periods = 2L
  )

  expect_true(!is.null(result))
  # Should only have anchor (event_time relative to tpost1=6) and t=2
  # event_times should be negative relative to tpost1=6
  expect_true(all(result$event_time < 0L))

  # Without exclusion: S=6, pre={1,2,3,4,5}, anchor=5, estimable={2,3,4}
  result_full <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 6L, rolling = "demean", estimator = "ra", alpha = 0.05,
    exclude_pre_periods = 0L
  )

  # Excluded version should have fewer rows
  expect_true(nrow(result) < nrow(result_full))
})

# Task 5.4: All pre-periods have insufficient sample -> all NA + anchor
test_that("Task 5.4: all periods insufficient sample returns NA rows + anchor", {
  # Only 2 units total (1 treated, 1 control) -> n_total=2 < 3
  dt <- data.table::data.table(
    id = rep(1:2, each = 6),
    time = rep(1:6, 2),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22),
    d = rep(c(1L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 4L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  # Should return data.frame (not NULL), anchor row exists
  expect_true(!is.null(result))
  expect_true(is.data.frame(result))
  expect_true(any(result$is_anchor == TRUE))

  # Non-anchor rows should have NA att (n_total=2 < 3)
  non_anchor <- result[result$is_anchor == FALSE, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})

# Task 5.5: cohort column is always NA_integer_
test_that("Task 5.5: cohort column is always NA_integer_", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 6),
    time = rep(1:6, 3),
    y = c(10, 12, 15, 20, 25, 30,
          8, 10, 13, 16, 19, 22,
          5, 7, 9, 12, 15, 18),
    d = rep(c(1L, 0L, 0L), each = 6)
  )

  result <- estimate_pre_treatment_common(
    data = dt, y = "y", ivar = "id", d = "d", tvar = "time",
    tpost1 = 5L, rolling = "demean", estimator = "ra", alpha = 0.05
  )

  expect_true("cohort" %in% names(result))
  expect_true(all(is.na(result$cohort)))
  expect_true(is.integer(result$cohort))
})
