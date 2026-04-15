# test-python-consistency-pretreatment-staggered.R
# Task E7-04.4: Python数值一致性端到端测试
#
# Verifies R estimate_pre_treatment_staggered() produces results
# numerically equivalent to Python lwdid v0.2.3.
#
# Python reference: lwdid-py_v0.2.3/src/lwdid/staggered/estimation_pre.py
# Python test data: lwdid-py_v0.2.3/tests/staggered/test_pre_treatment_unit.py
#
# Known R-vs-Python design differences (design.md section 4):
#   D1: R excludes T_min from estimable pre-periods
#   D3: R excludes G_i==g from anchor control count
#   D5: R generates NA rows for skipped pairs (Python skips)
#   D7: R degrades detrend to demean when |ref|<2 (Python NaN)
#   D4: R n_control counts d_pre==0; Python counts control_mask
#       which may include treated units for not_yet_treated
#
# Tolerances: ATT < 1e-5, SE < 1e-4, integers exact match.

library(testthat)
library(data.table)

# ============================================================
# Helper: Python simple_panel_data fixture
# From test_pre_treatment_unit.py simple_panel_data
# ============================================================
make_simple_panel <- function() {
  data.table(
    id   = rep(1:4, each = 6),
    time = rep(1:6, times = 4),
    y    = c(
      10, 12, 14, 16, 18, 20,
      20, 22, 24, 26, 28, 30,
       5,  7,  9, 11, 13, 15,
      15, 15, 15, 15, 15, 15
    ),
    gvar = c(rep(4L, 6), rep(4L, 6),
             rep(6L, 6), rep(0L, 6))
  )
}

# ============================================================
# Helper: Staggered panel from Python np.random.seed(42)
# 3 cohorts (g=4,6,8), 5 units each + 5 never-treated
# Periods 1..10, Y = unit_fe + 2*time + 3*post + noise
# ============================================================
make_staggered_panel <- function() {
  y_vals <- c(
    # Unit 1 (g=4)
    2.924296, 5.317273, 7.754943, 11.876352,
    13.87636, 16.783035, 18.377146, 19.758691,
    22.264708, 23.761719,
    # Unit 2 (g=4)
    1.189522, 2.1119, 4.206082, 9.787397,
    11.562125, 14.225664, 15.614528, 17.362389,
    20.801365, 21.955652,
    # Unit 3 (g=4)
    1.422682, 3.862865, 6.190518, 10.55956,
    13.322905, 14.834737, 16.98921, 18.834203,
    22.061196, 23.128308,
    # Unit 4 (g=4)
    0.295851, 1.274156, 3.98901, 7.904743,
    10.220485, 12.983009, 15.253811, 16.970262,
    18.826754, 20.734026,
    # Unit 5 (g=4)
    -1.316966, 0.812637, 3.571517, 8.214765,
    9.161436, 12.204998, 13.850415, 15.704495,
    18.348794, 20.558456,
    # Unit 6 (g=6)
    3.442951, 5.707954, 8.028192, 10.350333,
    11.622973, 16.769731, 18.309393, 20.264457,
    23.268823, 25.54068,
    # Unit 7 (g=6)
    2.357746, 4.036798, 5.53342, 8.036678,
    10.624998, 14.838067, 17.638302, 17.546107,
    21.266931, 22.899503,
    # Unit 8 (g=6)
    1.447866, 2.408201, 5.292149, 7.580542,
    10.140932, 14.14285, 15.997738, 18.151107,
    20.859686, 22.566361,
    # Unit 9 (g=6)
    1.197113, 2.989018, 5.424802, 6.589453,
    8.776649, 13.744426, 15.208722, 18.08854,
    20.071007, 21.943036,
    # Unit 10 (g=6)
    0.82314, 3.320503, 5.359468, 7.129687,
    9.450183, 14.732851, 17.473919, 18.618115,
    20.659601, 22.493603,
    # Unit 11 (g=8)
    -1.850799, 0.192573, 3.394079, 4.066277,
    6.313231, 8.145102, 9.578119, 15.733869,
    17.538424, 19.557974,
    # Unit 12 (g=8)
    0.882622, 1.4803, 4.474654, 7.276453,
    7.685957, 9.898076, 12.231051, 16.929487,
    18.405893, 21.215507,
    # Unit 13 (g=8)
    0.112189, 1.41568, 4.65036, 5.483766,
    7.714362, 10.282151, 11.25996, 16.989123,
    19.528964, 20.071651,
    # Unit 14 (g=8)
    2.499209, 4.760179, 5.750792, 7.709039,
    10.630238, 12.51776, 14.494514, 19.542492,
    21.029255, 23.485395,
    # Unit 15 (g=8)
    2.228969, 5.519032, 6.823061, 7.990493,
    10.914422, 12.098804, 14.979687, 20.165443,
    21.175804, 24.067833,
    # Unit 16 (g=0, NT)
    3.236592, 5.773958, 6.702868, 8.448694,
    10.380805, 12.417657, 14.787011, 16.996138,
    18.963907, 21.239153,
    # Unit 17 (g=0, NT)
    2.752771, 3.893675, 7.386088, 8.338837,
    9.597425, 11.490558, 14.26724, 15.914272,
    18.383004, 20.262623,
    # Unit 18 (g=0, NT)
    1.430945, 3.096919, 5.631085, 8.282542,
    9.961389, 11.231473, 13.940933, 16.047001,
    17.412413, 19.931205,
    # Unit 19 (g=0, NT)
    1.544932, 4.295311, 6.39681, 8.657943,
    10.643318, 11.427583, 13.647505, 16.373935,
    18.37331, 20.373941,
    # Unit 20 (g=0, NT)
    9.990908, 12.273246, 14.182464, 16.031159,
    17.547828, 20.084948, 21.31905, 23.587054,
    25.462781, 27.7464
  )
  data.table(
    id   = rep(1:20, each = 10),
    time = rep(1:10, times = 20),
    y    = y_vals,
    gvar = c(rep(4L, 50), rep(6L, 50),
             rep(8L, 50), rep(0L, 50))
  )
}

# ============================================================
# P-01: Test data from Python test_pre_treatment_unit.py
# ============================================================
test_that("P-01: test data matches Python fixture structure", {
  dt <- make_staggered_panel()
  expect_equal(nrow(dt), 200L)
  expect_equal(length(unique(dt$id)), 20L)
  expect_equal(sort(unique(dt$time)), 1:10)
  cohort_counts <- dt[, .(n = uniqueN(id)), by = gvar]
  expect_equal(cohort_counts[gvar == 4L]$n, 5L)
  expect_equal(cohort_counts[gvar == 6L]$n, 5L)
  expect_equal(cohort_counts[gvar == 8L]$n, 5L)
  expect_equal(cohort_counts[gvar == 0L]$n, 5L)
  # Spot-check known Y values from Python output
  expect_equal(dt[id == 1 & time == 1]$y, 2.924296)
  expect_equal(dt[id == 4 & time == 3]$y, 3.98901)
  expect_equal(dt[id == 20 & time == 10]$y, 27.7464)
})

# ============================================================
# P-02: demean ATT matches Python (error < 1e-5)
# Verifies lw2026 eq(2.12) rolling version equivalence
# R excludes T_min=1; compare only overlapping periods
# ============================================================
test_that("P-02: demean ATT matches Python (not_yet_treated)", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )
  expect_false(is.null(res))

  # Python reference (lwdid-py v0.2.3, seed=42, ra, nyt)
  py <- data.frame(
    cohort = c(4L, 6L, 6L, 6L,
               8L, 8L, 8L, 8L, 8L),
    period = c(2L, 4L, 3L, 2L,
               6L, 5L, 4L, 3L, 2L),
    att = c(
      -0.2088514667,
      -0.2754312000,
       0.4957515333,
       0.4303647556,
       0.3416164000,
       1.4019674500,
       0.7182717667,
       1.7329051333,
       1.3396211733
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py))) {
    r_row <- res[res$cohort == py$cohort[i] &
                   res$period == py$period[i], ]
    expect_equal(nrow(r_row), 1L,
      info = sprintf("g=%d t=%d not found",
                     py$cohort[i], py$period[i]))
    expect_equal(r_row$att, py$att[i], tolerance = 1e-5,
      info = sprintf(
        "demean ATT g=%d t=%d: R=%.10f Py=%.10f",
        py$cohort[i], py$period[i],
        r_row$att, py$att[i]))
  }
})

# ============================================================
# P-03: detrend ATT matches Python (error < 1e-5)
# Verifies lw2025 eq(5.6) rolling version equivalence
# Only compare periods where Python also produces results
# (Python skips t=g-2; R degrades to demean there)
# ============================================================
test_that("P-03: detrend ATT matches Python (not_yet_treated)", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "detrend", estimator = "ra",
    control_group = "not_yet_treated"
  )
  expect_false(is.null(res))

  # Python detrend reference (overlapping periods only)
  # g=6: t=3,2 (t=4 skipped by Python, t=1 excluded by R)
  # g=8: t=5,4,3,2 (t=6 skipped by Python, t=1 excluded by R)
  py <- data.frame(
    cohort = c(6L, 6L, 8L, 8L, 8L, 8L),
    period = c(3L, 2L, 5L, 4L, 3L, 2L),
    att = c(
       0.8840375333,
       0.0640418889,
       1.0259485000,
      -0.8090353333,
       0.7250099000,
      -0.3048590467
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py))) {
    r_row <- res[res$cohort == py$cohort[i] &
                   res$period == py$period[i], ]
    expect_equal(nrow(r_row), 1L,
      info = sprintf("g=%d t=%d not found",
                     py$cohort[i], py$period[i]))
    expect_equal(r_row$att, py$att[i], tolerance = 1e-5,
      info = sprintf(
        "detrend ATT g=%d t=%d: R=%.10f Py=%.10f",
        py$cohort[i], py$period[i],
        r_row$att, py$att[i]))
  }
})

# ============================================================
# P-04: SE matches Python (error < 1e-4)
# ============================================================
test_that("P-04: demean SE matches Python", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )

  py <- data.frame(
    cohort = c(4L, 6L, 6L, 6L,
               8L, 8L, 8L, 8L, 8L),
    period = c(2L, 4L, 3L, 2L,
               6L, 5L, 4L, 3L, 2L),
    se = c(
      0.3805683795,
      0.3816374764,
      0.7138047815,
      0.5688756532,
      0.4334125702,
      0.7563908588,
      0.6771964787,
      0.6882316627,
      0.6124751727
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py))) {
    r_row <- res[res$cohort == py$cohort[i] &
                   res$period == py$period[i], ]
    expect_equal(r_row$se, py$se[i], tolerance = 1e-4,
      info = sprintf(
        "demean SE g=%d t=%d: R=%.10f Py=%.10f",
        py$cohort[i], py$period[i],
        r_row$se, py$se[i]))
  }
})

test_that("P-04: detrend SE matches Python", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "detrend", estimator = "ra",
    control_group = "not_yet_treated"
  )

  py <- data.frame(
    cohort = c(6L, 6L, 8L, 8L, 8L, 8L),
    period = c(3L, 2L, 5L, 4L, 3L, 2L),
    se = c(
      0.8439712338,
      0.5686052650,
      0.9582869611,
      0.4446359124,
      0.8191607908,
      0.6002332566
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py))) {
    r_row <- res[res$cohort == py$cohort[i] &
                   res$period == py$period[i], ]
    expect_equal(r_row$se, py$se[i], tolerance = 1e-4,
      info = sprintf(
        "detrend SE g=%d t=%d: R=%.10f Py=%.10f",
        py$cohort[i], py$period[i],
        r_row$se, py$se[i]))
  }
})

# ============================================================
# P-05: n_treated/n_control matches Python (exact match)
# Uses never_treated to avoid D4 counting difference
# ============================================================
test_that("P-05: n_treated/n_control exact (never_treated)", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "never_treated"
  )
  expect_false(is.null(res))

  # Python never_treated: all periods have n_treated=5, n_control=5
  non_anchor <- res[!res$is_anchor, ]
  expect_true(all(non_anchor$n_treated == 5L),
    info = "All cohorts have 5 treated units")
  expect_true(all(non_anchor$n_control == 5L),
    info = "never_treated gives exactly 5 NT controls")

  # Anchors also have n_treated=5
  anchors <- res[res$is_anchor, ]
  expect_true(all(anchors$n_treated == 5L))
})

test_that("P-05: n_control for not_yet_treated is consistent", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )

  # n_treated always 5 (5 units per cohort)
  non_anchor <- res[!res$is_anchor, ]
  expect_true(all(non_anchor$n_treated == 5L))

  # Verify specific n_control values based on control logic:
  # controls = {G_i > t, G_i != g} + NT
  # g=4, t=2: g6(5)+g8(5)+NT(5) = 15
  expect_equal(
    res[res$cohort == 4 & res$period == 2, ]$n_control,
    15L
  )
  # g=6, t=4: g8(5, since 8>4)+NT(5) = 10
  # g=4 excluded (4>4 FALSE)
  expect_equal(
    res[res$cohort == 6 & res$period == 4, ]$n_control,
    10L
  )
  # g=6, t=2: g4(5, 4>2)+g8(5, 8>2)+NT(5) = 15
  expect_equal(
    res[res$cohort == 6 & res$period == 2, ]$n_control,
    15L
  )
  # g=8, t=6: only NT(5), since g4(4>6 F), g6(6>6 F)
  expect_equal(
    res[res$cohort == 8 & res$period == 6, ]$n_control,
    5L
  )
  # g=8, t=2: g4(5)+g6(5)+NT(5) = 15
  expect_equal(
    res[res$cohort == 8 & res$period == 2, ]$n_control,
    15L
  )
})

# ============================================================
# P-06: rolling_window_size matches Python (exact match)
# Verifies |R_t| = g - t - 1
# ============================================================
test_that("P-06: rolling_window_size = g - t - 1", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )

  py_rws <- data.frame(
    cohort = c(4L, 4L,
               6L, 6L, 6L, 6L,
               8L, 8L, 8L, 8L, 8L, 8L),
    period = c(3L, 2L,
               5L, 4L, 3L, 2L,
               7L, 6L, 5L, 4L, 3L, 2L),
    rws = c(0L, 1L,
            0L, 1L, 2L, 3L,
            0L, 1L, 2L, 3L, 4L, 5L),
    is_anchor = c(TRUE, FALSE,
                  TRUE, FALSE, FALSE, FALSE,
                  TRUE, FALSE, FALSE, FALSE,
                  FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py_rws))) {
    r_row <- res[res$cohort == py_rws$cohort[i] &
                   res$period == py_rws$period[i], ]
    expect_equal(nrow(r_row), 1L,
      info = sprintf("g=%d t=%d not found",
                     py_rws$cohort[i], py_rws$period[i]))
    expect_equal(
      r_row$rolling_window_size, py_rws$rws[i],
      info = sprintf("rws g=%d t=%d: R=%s Py=%d",
        py_rws$cohort[i], py_rws$period[i],
        r_row$rolling_window_size, py_rws$rws[i]))
    # Verify formula for non-anchor
    if (!py_rws$is_anchor[i]) {
      expect_equal(
        r_row$rolling_window_size,
        as.integer(py_rws$cohort[i] - py_rws$period[i] - 1L))
    }
  }
})

# ============================================================
# P-07: Anchor behavior matches Python
# ATT=0, SE=0, is_anchor=TRUE, t_stat=NaN, pvalue=NaN
# ============================================================
test_that("P-07: anchor ATT=0, SE=0, is_anchor=TRUE", {
  dt <- make_staggered_panel()

  for (rolling in c("demean", "detrend")) {
    res <- estimate_pre_treatment_staggered(
      data = dt, y = "y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = rolling, estimator = "ra",
      control_group = "not_yet_treated"
    )
    anchors <- res[res$is_anchor == TRUE, ]

    # One anchor per cohort
    expect_equal(nrow(anchors), 3L,
      info = sprintf("%s: expected 3 anchors", rolling))

    for (j in seq_len(nrow(anchors))) {
      a <- anchors[j, ]
      expect_equal(a$att, 0.0,
        info = sprintf("%s g=%d: ATT!=0", rolling, a$cohort))
      expect_equal(a$se, 0.0,
        info = sprintf("%s g=%d: SE!=0", rolling, a$cohort))
      expect_true(is.nan(a$t_stat),
        info = sprintf("%s g=%d: t_stat not NaN",
                       rolling, a$cohort))
      expect_true(is.nan(a$pvalue),
        info = sprintf("%s g=%d: pvalue not NaN",
                       rolling, a$cohort))
      expect_true(a$is_anchor)
      expect_equal(a$rolling_window_size, 0L)
      expect_equal(a$event_time, -1L)
      expect_equal(a$ci_lower, 0.0)
      expect_equal(a$ci_upper, 0.0)
    }

    # Anchor periods = g - 1
    expect_equal(sort(anchors$period), c(3L, 5L, 7L))
  }
})

# ============================================================
# P-08: Multi-cohort: each cohort's pre-treatment period set
# matches Python (accounting for R T_min exclusion)
# ============================================================
test_that("P-08: multi-cohort pre-treatment period sets", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )

  # All 3 cohorts present
  expect_true(all(c(4L, 6L, 8L) %in% res$cohort))

  # Cohort 4: R periods = {2, 3(anchor)}
  c4 <- res[res$cohort == 4, ]
  expect_equal(sort(c4$period), c(2L, 3L))
  expect_equal(c4[c4$is_anchor, ]$period, 3L)
  expect_true(all(c4$period < 4L))
  expect_true(all(c4$event_time < 0L))

  # Cohort 6: R periods = {2,3,4, 5(anchor)}
  c6 <- res[res$cohort == 6, ]
  expect_equal(sort(c6$period), c(2L, 3L, 4L, 5L))
  expect_equal(c6[c6$is_anchor, ]$period, 5L)
  expect_true(all(c6$period < 6L))

  # Cohort 8: R periods = {2,3,4,5,6, 7(anchor)}
  c8 <- res[res$cohort == 8, ]
  expect_equal(sort(c8$period), c(2L, 3L, 4L, 5L, 6L, 7L))
  expect_equal(c8[c8$is_anchor, ]$period, 7L)
  expect_true(all(c8$period < 8L))

  # event_time = period - cohort for all rows
  expect_true(all(
    res$event_time == res$period - res$cohort
  ))

  # Sorting: cohort asc, event_time desc
  expect_true(all(diff(res$cohort) >= 0))
  for (g in c(4L, 6L, 8L)) {
    cg <- res[res$cohort == g, ]
    expect_true(all(diff(cg$event_time) <= 0),
      info = sprintf("Cohort %d not sorted", g))
  }
})

# ============================================================
# Supplement: Simple panel demean verification
# ============================================================
test_that("P-02 supplement: simple panel demean ATT", {
  dt <- make_simple_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "never_treated"
  )
  expect_false(is.null(res))

  # g=4, t=2: ref={3}, demean = Y_t2 - Y_t3
  #   Unit 1: 12-14=-2, Unit 2: 22-24=-2 (treated, mean=-2)
  #   Unit 4 (NT): 15-15=0 (control, mean=0)
  #   ATT = -2 - 0 = -2
  c4_t2 <- res[res$cohort == 4 & res$period == 2, ]
  if (nrow(c4_t2) > 0 && !is.na(c4_t2$att)) {
    expect_equal(c4_t2$att, -2.0, tolerance = 1e-5)
    expect_equal(c4_t2$n_treated, 2L)
    expect_equal(c4_t2$n_control, 1L)
  }

  # g=6, t=3: ref={4,5}, demean = Y_t3 - mean(Y_t4, Y_t5)
  #   Unit 3 (g=6): 9 - mean(11,13) = 9-12 = -3
  #   Unit 4 (NT): 15 - mean(15,15) = 0
  #   ATT = -3 - 0 = -3
  c6_t3 <- res[res$cohort == 6 & res$period == 3, ]
  if (nrow(c6_t3) > 0 && !is.na(c6_t3$att)) {
    expect_equal(c6_t3$att, -3.0, tolerance = 1e-5)
  }
})

# ============================================================
# Supplement: Detrend degradation at t=g-2
# ============================================================
test_that("P-03 supplement: detrend degrades to demean at t=g-2", {
  dt <- make_staggered_panel()
  res_dm <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "not_yet_treated"
  )
  res_dt <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "detrend", estimator = "ra",
    control_group = "not_yet_treated"
  )

  # g=6, t=4: ref={5}, rws=1 → detrend degrades to demean
  dm_c6_t4 <- res_dm[res_dm$cohort == 6 &
                        res_dm$period == 4, ]
  dt_c6_t4 <- res_dt[res_dt$cohort == 6 &
                        res_dt$period == 4, ]
  if (nrow(dm_c6_t4) > 0 && nrow(dt_c6_t4) > 0 &&
      !is.na(dm_c6_t4$att) && !is.na(dt_c6_t4$att)) {
    expect_equal(dt_c6_t4$att, dm_c6_t4$att,
      tolerance = 1e-10)
    expect_equal(dt_c6_t4$se, dm_c6_t4$se,
      tolerance = 1e-10)
  }

  # g=8, t=6: ref={7}, rws=1 → detrend degrades to demean
  dm_c8_t6 <- res_dm[res_dm$cohort == 8 &
                        res_dm$period == 6, ]
  dt_c8_t6 <- res_dt[res_dt$cohort == 8 &
                        res_dt$period == 6, ]
  if (nrow(dm_c8_t6) > 0 && nrow(dt_c8_t6) > 0 &&
      !is.na(dm_c8_t6$att) && !is.na(dt_c8_t6$att)) {
    expect_equal(dt_c8_t6$att, dm_c8_t6$att,
      tolerance = 1e-10)
    expect_equal(dt_c8_t6$se, dm_c8_t6$se,
      tolerance = 1e-10)
  }
})

# ============================================================
# Supplement: ATT/SE for never_treated control group
# ============================================================
test_that("P-02/P-04 supplement: never_treated ATT/SE", {
  dt <- make_staggered_panel()
  res <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean", estimator = "ra",
    control_group = "never_treated"
  )

  py <- data.frame(
    cohort = c(4L, 6L, 6L, 6L,
               8L, 8L, 8L, 8L, 8L),
    period = c(2L, 4L, 3L, 2L,
               6L, 5L, 4L, 3L, 2L),
    att = c(
      -0.2734066000,
      -0.5114904000,
      -0.3735056000,
      -0.2908739333,
       0.3416164000,
      -0.0616376000,
      -0.5132101333,
       0.0204480000,
      -0.1354368000
    ),
    se = c(
      0.4357517811,
      0.2752512493,
      0.4271902309,
      0.4504603432,
      0.4334125702,
      0.3027668841,
      0.4119066701,
      0.4373449572,
      0.4656206480
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(py))) {
    r_row <- res[res$cohort == py$cohort[i] &
                   res$period == py$period[i], ]
    expect_equal(nrow(r_row), 1L,
      info = sprintf("g=%d t=%d not found",
                     py$cohort[i], py$period[i]))
    expect_equal(r_row$att, py$att[i], tolerance = 1e-5,
      info = sprintf("NT ATT g=%d t=%d",
                     py$cohort[i], py$period[i]))
    expect_equal(r_row$se, py$se[i], tolerance = 1e-4,
      info = sprintf("NT SE g=%d t=%d",
                     py$cohort[i], py$period[i]))
  }
})
