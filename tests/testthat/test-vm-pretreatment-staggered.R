# test-vm-pretreatment-staggered.R
# Vibe-math cross-validation for E7-04 Task 3 (V-01 to V-08)

# Helper: create staggered panel
make_staggered_panel <- function() {
  # 5 units x 6 periods
  # U1(g=4): Y=10,20,25,30,35,40
  # U2(g=4): Y=8,15,18,22,26,30
  # U3(g=6): Y=12,22,28,35,42,50
  # U4(NT):  Y=5,10,12,15,18,21
  # U5(NT):  Y=3,8,9,11,14,17
  data.table::data.table(
    id = rep(1:5, each = 6),
    time = rep(1:6, 5),
    y = c(10, 20, 25, 30, 35, 40,
          8, 15, 18, 22, 26, 30,
          12, 22, 28, 35, 42, 50,
          5, 10, 12, 15, 18, 21,
          3, 8, 9, 11, 14, 17),
    g = rep(c(4L, 4L, 6L, Inf, Inf), each = 6)
  )
}

# V-01: demean transform hand-calculation
test_that("V-01: demean transform matches vibe-math (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # Cohort g=4: pre={1,2,3}, anchor=3, T_min=1 excluded, estimable={2}
  # t=2, ref={3}: demean = Y_i2 - Y_i3
  # U1: 20-25=-5, U2: 15-18=-3 (treated)
  # Controls at t=2 with not_yet_treated: G_i>2 and G_i!=4 -> U3(g=6>2,!=4)
  # + NT: U4, U5
  # U3: 22-28=-6, U4: 10-12=-2, U5: 8-9=-1
  # ATT = mean(-5,-3) - mean(-6,-2,-1) = -4 - (-3) = -1
  g4_t2 <- result[result$cohort == 4L & result$event_time == -2L, ]
  expect_equal(nrow(g4_t2), 1L)
  expect_equal(g4_t2$att, -1.0, tolerance = 1e-10)
  expect_equal(g4_t2$n_treated, 2L)
  expect_equal(g4_t2$n_control, 3L)
  expect_equal(g4_t2$rolling_window_size, 1L)
})

# V-02/V-03: detrend OLS coefficients and prediction
test_that("V-02/V-03: detrend OLS matches vibe-math (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "detrend", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # Cohort g=6: pre={1,2,3,4,5}, anchor=5, T_min=1 excluded
  # estimable={2,3,4}
  # t=3, ref={4,5}: OLS Y_iq = alpha + beta*q
  #   U3(g=6,treated): (4,35),(5,42) -> beta=7, alpha=35-28=7 -> pred(3)=7+21=28 -> dt=28-28=0
  #   U1(g=4,ctrl g>3): (4,30),(5,35) -> beta=5, alpha=30-20=10 -> pred(3)=10+15=25 -> dt=25-25=0
  #   U2(g=4,ctrl g>3): (4,22),(5,26) -> beta=4, alpha=22-16=6 -> pred(3)=6+12=18 -> dt=18-18=0
  #   U4(NT): (4,15),(5,18) -> beta=3, alpha=15-12=3 -> pred(3)=3+9=12 -> dt=12-12=0
  #   U5(NT): (4,11),(5,14) -> beta=3, alpha=11-12=-1 -> pred(3)=-1+9=8 -> dt=9-8=1
  # ATT = 0 - (0+0+0+1)/4 = -0.25
  g6_t3 <- result[result$cohort == 6L & result$event_time == (3L - 6L), ]
  expect_equal(nrow(g6_t3), 1L)
  expect_equal(g6_t3$att, -0.25, tolerance = 1e-10)
  expect_equal(g6_t3$rolling_window_size, 2L)  # |{4,5}| = 2
})

# V-04: rolling_window_size = g - t_pre - 1
test_that("V-04: rolling_window_size matches vibe-math (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # g=4, t=2: rws = 4-2-1 = 1 (ref={3})
  g4_t2 <- result[result$cohort == 4L & result$event_time == -2L, ]
  expect_equal(g4_t2$rolling_window_size, 1L)

  # g=6, t=2: rws = 6-2-1 = 3 (ref={3,4,5})
  g6_t2 <- result[result$cohort == 6L & result$event_time == (2L - 6L), ]
  expect_equal(g6_t2$rolling_window_size, 3L)

  # g=6, t=3: rws = 6-3-1 = 2 (ref={4,5})
  g6_t3 <- result[result$cohort == 6L & result$event_time == (3L - 6L), ]
  expect_equal(g6_t3$rolling_window_size, 2L)

  # g=6, t=4: rws = 6-4-1 = 1 (ref={5})
  g6_t4 <- result[result$cohort == 6L & result$event_time == (4L - 6L), ]
  expect_equal(g6_t4$rolling_window_size, 1L)

  # Anchors: rws = 0
  g4_anchor <- result[result$cohort == 4L & result$is_anchor == TRUE, ]
  g6_anchor <- result[result$cohort == 6L & result$is_anchor == TRUE, ]
  expect_equal(g4_anchor$rolling_window_size, 0L)
  expect_equal(g6_anchor$rolling_window_size, 0L)
})

# V-05: anchor df_inference = n_treat + n_ctrl - 2
test_that("V-05: anchor df_inference matches vibe-math (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # g=4 anchor (t=3): n_treat=2 (U1,U2), n_ctrl=3 (U3,U4,U5)
  # df = 2+3-2 = 3
  g4_anchor <- result[result$cohort == 4L & result$is_anchor == TRUE, ]
  expect_equal(g4_anchor$n_treated, 2L)
  expect_equal(g4_anchor$n_control, 3L)
  expect_equal(g4_anchor$df_inference, 3L)

  # g=6 anchor (t=5): n_treat=1 (U3), n_ctrl=2 (U4,U5)
  # G_i>5 and G_i!=6: no finite cohort qualifies (g=4 is NOT >5)
  # Only NT: U4,U5 -> n_ctrl=2
  # df = 1+2-2 = 1
  g6_anchor <- result[result$cohort == 6L & result$is_anchor == TRUE, ]
  expect_equal(g6_anchor$n_treated, 1L)
  expect_equal(g6_anchor$n_control, 2L)
  expect_equal(g6_anchor$df_inference, 1L)
})

# V-06: symmetric reference periods = {t+1,...,g-1}
test_that("V-06: symmetric ref periods verified (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # g=4, t=2: ref = {3} -> rws=1
  # g=6, t=2: ref = {3,4,5} -> rws=3
  # g=6, t=3: ref = {4,5} -> rws=2
  # g=6, t=4: ref = {5} -> rws=1
  # Verify via rolling_window_size = |ref| = g-t-1
  g4_t2 <- result[result$cohort == 4L & result$event_time == -2L, ]
  g6_t2 <- result[result$cohort == 6L & result$event_time == (2L - 6L), ]
  g6_t3 <- result[result$cohort == 6L & result$event_time == (3L - 6L), ]
  g6_t4 <- result[result$cohort == 6L & result$event_time == (4L - 6L), ]

  # Symmetric: ref uses FUTURE pre-treatment periods {t+1,...,g-1}
  # NOT all pre-treatment periods. Verify rws = g-t-1.
  expect_equal(g4_t2$rolling_window_size, as.integer(4L - 2L - 1L))
  expect_equal(g6_t2$rolling_window_size, as.integer(6L - 2L - 1L))
  expect_equal(g6_t3$rolling_window_size, as.integer(6L - 3L - 1L))
  expect_equal(g6_t4$rolling_window_size, as.integer(6L - 4L - 1L))

  # Also verify ATT at g=6 t=2 (demean, ref={3,4,5}):
  # U3(treated): Y_2 - mean(Y_3,Y_4,Y_5) = 22 - (28+35+42)/3 = 22-35 = -13
  # U1(ctrl): 20 - (25+30+35)/3 = 20-30 = -10
  # U2(ctrl): 15 - (18+22+26)/3 = 15-22 = -7
  # U4(ctrl): 10 - (12+15+18)/3 = 10-15 = -5
  # U5(ctrl): 8 - (9+11+14)/3 = 8-34/3 = -10/3
  # ATT = -13 - (-10-7-5-10/3)/4 = -13 - (-76/12) = -13+19/3 = -20/3
  expect_equal(g6_t2$att, -20 / 3, tolerance = 1e-10)
})

# V-07: detrend degrades to demean when |ref|=1
test_that("V-07: detrend degrades to demean when ref=1 (staggered)", {
  dt <- make_staggered_panel()

  result_detrend <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "detrend", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )
  result_demean <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # g=4, t=2: ref={3} (1 period) -> detrend degrades to demean
  g4_t2_dt <- result_detrend[result_detrend$cohort == 4L &
                               result_detrend$event_time == -2L, ]
  g4_t2_dm <- result_demean[result_demean$cohort == 4L &
                              result_demean$event_time == -2L, ]
  expect_equal(g4_t2_dt$att, g4_t2_dm$att, tolerance = 1e-10)
  # Both should equal -1.0 (verified in V-01)
  expect_equal(g4_t2_dt$att, -1.0, tolerance = 1e-10)

  # g=6, t=4: ref={5} (1 period) -> detrend degrades to demean
  # U3(treated): Y_4 - Y_5 = 35-42 = -7
  # Controls at t=4: G_i>4 and G_i!=6 -> none (g=4 NOT >4) + NT: U4,U5
  # U4: 15-18=-3, U5: 11-14=-3
  # ATT = -7 - (-3-3)/2 = -7-(-3) = -4
  g6_t4_dt <- result_detrend[result_detrend$cohort == 6L &
                               result_detrend$event_time == (4L - 6L), ]
  g6_t4_dm <- result_demean[result_demean$cohort == 6L &
                              result_demean$event_time == (4L - 6L), ]
  expect_equal(g6_t4_dt$att, g6_t4_dm$att, tolerance = 1e-10)
  expect_equal(g6_t4_dt$att, -4.0, tolerance = 1e-10)
})

# V-08: multi-cohort reference period independence
test_that("V-08: multi-cohort ref independence (staggered)", {
  dt <- make_staggered_panel()
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "g",
    rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated", alpha = 0.05
  )

  # Same period t=2 but different cohorts -> different ref sets
  # g=4, t=2: ref={3} -> ATT=-1.0 (V-01)
  # g=6, t=2: ref={3,4,5} -> ATT=-20/3 (V-06)
  g4_t2 <- result[result$cohort == 4L & result$event_time == -2L, ]
  g6_t2 <- result[result$cohort == 6L & result$event_time == (2L - 6L), ]

  expect_equal(g4_t2$att, -1.0, tolerance = 1e-10)
  expect_equal(g6_t2$att, -20 / 3, tolerance = 1e-10)

  # ATTs must differ because ref sets differ
  expect_true(abs(g4_t2$att - g6_t2$att) > 1e-6)

  # rolling_window_size also differs
  expect_equal(g4_t2$rolling_window_size, 1L)  # |{3}|
  expect_equal(g6_t2$rolling_window_size, 3L)  # |{3,4,5}|
})
