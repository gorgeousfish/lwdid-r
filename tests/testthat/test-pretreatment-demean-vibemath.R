# ============================================================================
# test-pretreatment-demean-vibemath.R
# Task 4.1 (E7-03): Verify demean symmetric transformation against
# vibe-math manual computation.
#
# Panel: 3 units × 6 periods, S=4 (tpost1=4)
#   Unit 1 (treated):  Y = [3, 8, 11, 15, 19, 23]
#   Unit 2 (control):  Y = [2, 5, 10, 13, 18, 21]
#   Unit 3 (control):  Y = [4, 9, 12, 17, 20, 25]
#
# Pre-periods: {1, 2, 3}, anchor = 3 (S-1), T_min = 1 excluded
# Estimable periods: {2}
#
# Demean transform for period t:
#   ref_periods = {t+1, ..., S-1} = {t+1, ..., 3}
#   Y_dot_it = Y_it - mean(Y_i,ref_periods)
#
# vibe-math computed values:
#   t=1 (T_min, excluded by function):
#     ref = {2, 3}
#     Unit1: 3 - (8+11)/2 = 3 - 9.5 = -6.5
#     Unit2: 2 - (5+10)/2 = 2 - 7.5 = -5.5
#     Unit3: 4 - (9+12)/2 = 4 - 10.5 = -6.5
#     ATT_t1 = -6.5 - (-5.5 + -6.5)/2 = -6.5 - (-6.0) = -0.5
#
#   t=2 (estimable):
#     ref = {3}
#     Unit1: 8 - 11 = -3
#     Unit2: 5 - 10 = -5
#     Unit3: 9 - 12 = -3
#     ATT_t2 = -3 - (-5 + -3)/2 = -3 - (-4) = 1.0
# ============================================================================

test_that("Task 4.1: demean symmetric transform matches vibe-math manual computation", {
  # --- Construct 3-unit × 6-period panel ---
  dt <- data.table::data.table(
    unit = rep(1:3, each = 6L),
    time = rep(1:6, 3L),
    y = c(
      3, 8, 11, 15, 19, 23,   # Unit 1 (treated)
      2, 5, 10, 13, 18, 21,   # Unit 2 (control)
      4, 9, 12, 17, 20, 25    # Unit 3 (control)
    ),
    d = rep(c(1L, 0L, 0L), each = 6L)
  )

  # --- Run function with demean ---
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L, rolling = "demean"
  )

  # --- Verify structure ---
  expect_s3_class(res, "data.frame")
  # Should have 2 rows: anchor (t=3) + estimable (t=2)
  expect_equal(nrow(res), 2L)

  # --- Verify anchor row (t=3, event_time=-1) ---
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  expect_equal(anchor$period, 3L)
  expect_equal(anchor$event_time, -1L)
  expect_equal(anchor$att, 0.0)
  expect_equal(anchor$rolling_window_size, 0L)

  # --- Verify estimable period t=2 ---
  t2_row <- res[res$period == 2L, ]
  expect_equal(nrow(t2_row), 1L)
  expect_false(t2_row$is_anchor)
  expect_equal(t2_row$event_time, -2L)

  # vibe-math expected ATT at t=2: 1.0
  # Tolerance < 1e-10 as specified
  expect_equal(t2_row$att, 1.0, tolerance = 1e-10)

  # rolling_window_size = |{t+1,...,S-1}| = S-1-t = 4-1-2 = 1
  expect_equal(t2_row$rolling_window_size, 1L)

  # --- Verify t=1 is NOT in output (T_min excluded) ---
  t1_rows <- res[res$period == 1L, ]
  expect_equal(nrow(t1_rows), 0L)

  # --- Verify SE is positive and finite ---
  expect_true(t2_row$se > 0)
  expect_true(is.finite(t2_row$se))

  # --- Verify sample counts ---
  expect_equal(t2_row$n_treated, 1L)
  expect_equal(t2_row$n_control, 2L)
})


# ============================================================================
# Additional verification: 4-unit panel with 2 treated, 2 control
# to get a richer ATT computation
# ============================================================================

test_that("Task 4.1 extended: 4-unit demean matches vibe-math", {
  # 4 units × 6 periods, S=4
  # Unit 1 (treated): Y = [1, 5, 9, 14, 18, 22]
  # Unit 2 (treated): Y = [2, 7, 10, 16, 20, 24]
  # Unit 3 (control): Y = [3, 6, 11, 15, 19, 23]
  # Unit 4 (control): Y = [4, 8, 12, 17, 21, 25]
  #
  # t=2, ref={3}:
  #   Unit1: 5 - 9 = -4
  #   Unit2: 7 - 10 = -3
  #   Unit3: 6 - 11 = -5
  #   Unit4: 8 - 12 = -4
  #   ATT = mean_treat - mean_ctrl = (-4 + -3)/2 - (-5 + -4)/2
  #       = -3.5 - (-4.5) = 1.0

  dt <- data.table::data.table(
    unit = rep(1:4, each = 6L),
    time = rep(1:6, 4L),
    y = c(
      1, 5, 9, 14, 18, 22,    # Unit 1 (treated)
      2, 7, 10, 16, 20, 24,   # Unit 2 (treated)
      3, 6, 11, 15, 19, 23,   # Unit 3 (control)
      4, 8, 12, 17, 21, 25    # Unit 4 (control)
    ),
    d = rep(c(1L, 1L, 0L, 0L), each = 6L)
  )

  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L, rolling = "demean"
  )

  t2_row <- res[res$period == 2L, ]
  expect_equal(nrow(t2_row), 1L)

  # vibe-math expected ATT at t=2: 1.0
  expect_equal(t2_row$att, 1.0, tolerance = 1e-10)
  expect_equal(t2_row$n_treated, 2L)
  expect_equal(t2_row$n_control, 2L)
})
