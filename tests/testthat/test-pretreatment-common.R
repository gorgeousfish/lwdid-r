# ============================================================================
# test-pretreatment-common.R
# Unit tests for estimate_pre_treatment_common()
# Story E7-03: Pre-treatment Effect Estimation (Common Timing)
# Paper: lw2025 Section 5 (equations 5.1-5.7, Procedure 5.1)
# ============================================================================

# ============================================================================
# Section 0: Test DGP Helper (TC-7.3.1 / TC-7.3.2 setup)
# ============================================================================

#' Generate panel data for pre-treatment testing
#' @param n_treat Number of treated units
#' @param n_ctrl Number of control units
#' @param T_total Total time periods
#' @param tpost1 First post-treatment period (treatment time S)
#' @param pt_violation Logical, if TRUE add pre-treatment trend
#'   difference (delta_t != 0)
#' @param delta_slope Slope of pre-treatment violation trend
#' @param seed Random seed
#' @return data.table with columns: unit, time, y, d, x1
make_pretreat_dgp <- function(n_treat = 20L, n_ctrl = 20L,
                               T_total = 8L, tpost1 = 5L,
                               pt_violation = FALSE,
                               delta_slope = 0.5,
                               seed = 42L) {
  set.seed(seed)
  N <- n_treat + n_ctrl
  d_vec <- c(rep(1L, n_treat), rep(0L, n_ctrl))

  # Unit fixed effects
  alpha_i <- rnorm(N, mean = 0, sd = 1)
  # Common time trend
  gamma_t <- seq(0, by = 0.3, length.out = T_total)

  rows <- vector("list", N * T_total)
  idx <- 0L
  for (i in seq_len(N)) {
    for (t in seq_len(T_total)) {
      idx <- idx + 1L
      # Base outcome: unit FE + time trend + noise
      y_it <- alpha_i[i] + gamma_t[t] + rnorm(1, sd = 0.5)
      # Pre-treatment violation: treated units get extra trend
      if (pt_violation && d_vec[i] == 1L && t < tpost1) {
        y_it <- y_it + delta_slope * (t - tpost1)
      }
      # Post-treatment effect (not relevant for pre-treatment test
      # but included for completeness)
      if (d_vec[i] == 1L && t >= tpost1) {
        y_it <- y_it + 2.0
      }
      # Control variable (time-invariant)
      x1_i <- alpha_i[i] * 0.5 + rnorm(1, sd = 0.1)
      rows[[idx]] <- list(
        unit = i, time = t, y = y_it,
        d = d_vec[i], x1 = x1_i
      )
    }
  }
  data.table::rbindlist(rows)
}


# ============================================================================
# Section 1: Input Validation (TC-7.3.1 prerequisites)
# ============================================================================

test_that("empty data raises lwdid_invalid_param", {
  dt <- data.table::data.table(
    unit = integer(0), time = integer(0),
    y = numeric(0), d = integer(0)
  )
  expect_error(
    estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L
    ),
    class = "lwdid_invalid_param"
  )
})

test_that("missing columns raises lwdid_invalid_param", {
  dt <- data.table::data.table(unit = 1:3, time = 1:3, y = 1:3)
  expect_error(
    estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L
    ),
    class = "lwdid_invalid_param"
  )
})

test_that("invalid rolling raises lwdid_invalid_param", {
  dt <- make_pretreat_dgp()
  expect_error(
    estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L, rolling = "bad"
    ),
    class = "lwdid_invalid_param"
  )
})

test_that("invalid alpha raises lwdid_invalid_param", {
  dt <- make_pretreat_dgp()
  expect_error(
    estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L, alpha = 1.5
    ),
    class = "lwdid_invalid_param"
  )
  expect_error(
    estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L, alpha = 0
    ),
    class = "lwdid_invalid_param"
  )
})


# ============================================================================
# Section 2: TC-7.3.1 — PT holds, pre-treatment ATT ≈ 0 (AC-13)
# ============================================================================

test_that("TC-7.3.1: parallel trends hold => pre-ATT near zero", {
  dt <- make_pretreat_dgp(
    n_treat = 20L, n_ctrl = 20L, T_total = 8L,
    tpost1 = 5L, pt_violation = FALSE, seed = 123L
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, rolling = "demean"
  )
  expect_s3_class(res, "data.frame")
  # Non-anchor rows should have ATT near 0
  non_anchor <- res[!res$is_anchor, ]
  # With N=40 and no violation, ATTs should be small
  # (within ~2 SE of zero for most)
  expect_true(all(abs(non_anchor$att) < 2.0),
    info = "Pre-treatment ATTs should be near zero under PT"
  )
})

# ============================================================================
# Section 3: TC-7.3.2 — PT violated, pre-treatment ATT ≠ 0 (AC-14)
# ============================================================================

test_that("TC-7.3.2: parallel trends violated => pre-ATT nonzero", {
  dt <- make_pretreat_dgp(
    n_treat = 20L, n_ctrl = 20L, T_total = 8L,
    tpost1 = 5L, pt_violation = TRUE,
    delta_slope = 1.0, seed = 456L
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, rolling = "demean"
  )
  non_anchor <- res[!res$is_anchor, ]
  # At least some pre-treatment ATTs should be meaningfully
  # different from zero
  expect_true(any(abs(non_anchor$att) > 0.3),
    info = "Pre-treatment ATTs should detect PT violation"
  )
})


# ============================================================================
# Section 4: TC-7.3.3 — event_time correct (AC-11)
# ============================================================================

test_that("TC-7.3.3: event_time = period - tpost1 (negative)", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 100L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  # All event_times should be negative (pre-treatment)
  expect_true(all(res$event_time < 0))
  # event_time = period - tpost1
  expect_equal(res$event_time, res$period - 5L)
  # Anchor should have event_time = -1
  anchor <- res[res$is_anchor, ]
  expect_equal(anchor$event_time, -1L)
})

# ============================================================================
# Section 5: TC-7.3.11 — rolling_window_size correct (AC-09)
# ============================================================================

test_that("TC-7.3.11: rolling_window_size = S-1-t, anchor=0", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 200L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  anchor <- res[res$is_anchor, ]
  expect_equal(anchor$rolling_window_size, 0L)

  non_anchor <- res[!res$is_anchor, ]
  # rolling_window_size = |{t+1,...,S-1}| = S-1-t = tpost1-1-period
  expected_rws <- as.integer(5L - 1L - non_anchor$period)
  expect_equal(non_anchor$rolling_window_size, expected_rws)
})


# ============================================================================
# Section 6: TC-7.3.4 — Controls with RA estimator
# ============================================================================

test_that("TC-7.3.4: RA estimator with controls runs", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 300L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L,
    controls = "x1", estimator = "ra"
  )
  expect_s3_class(res, "data.frame")
  expect_true(ncol(res) == 14L)
  # Should have anchor + estimable periods
  expect_true(nrow(res) >= 1L)
  # ATT values should be finite for non-NA rows
  non_na <- res[!is.na(res$att), ]
  expect_true(all(is.finite(non_na$att) | non_na$is_anchor))
})

# ============================================================================
# Section 7: TC-7.3.5 — VCE="hc3" SE correct
# ============================================================================

test_that("TC-7.3.5: VCE=hc3 produces valid SE", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 400L)
  res_default <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  res_hc3 <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, vce = "hc3"
  )
  # HC3 SEs should be >= homoskedastic SEs (generally)
  non_anchor_def <- res_default[!res_default$is_anchor, ]
  non_anchor_hc3 <- res_hc3[!res_hc3$is_anchor, ]
  # Both should have same structure
  expect_equal(nrow(non_anchor_def), nrow(non_anchor_hc3))
  # HC3 SEs should be positive
  expect_true(all(non_anchor_hc3$se > 0, na.rm = TRUE))
})


# ============================================================================
# Section 8: TC-7.3.8 — IPWRA estimator routing
# ============================================================================

test_that("TC-7.3.8: estimator=ipwra routes correctly", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 500L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L,
    controls = "x1", estimator = "ipwra"
  )
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 14L)
  # Should produce results (may have some NA if estimation fails)
  expect_true(nrow(res) >= 1L)
})

# ============================================================================
# Section 9: TC-7.3.9 — PSM estimator routing
# ============================================================================

test_that("TC-7.3.9: estimator=psm routes correctly", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 600L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L,
    controls = "x1", estimator = "psm"
  )
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 14L)
  expect_true(nrow(res) >= 1L)
})


# ============================================================================
# Section 10: TC-7.3.6 — Minimal pre-periods (only 2)
# ============================================================================

test_that("TC-7.3.6: only 2 pre-periods => anchor only", {
  # With T_total=3, tpost1=3: periods 1,2 are pre-treatment
  # T_min=1 excluded, anchor=2, no estimable periods
  dt <- make_pretreat_dgp(
    T_total = 3L, tpost1 = 3L, seed = 700L
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 3L
  )
  expect_s3_class(res, "data.frame")
  # Should have only anchor row
  expect_equal(nrow(res), 1L)
  expect_true(res$is_anchor[1])
  expect_equal(res$att[1], 0.0)
})

# ============================================================================
# Section 11: TC-7.3.7 — Insufficient sample returns NA
# ============================================================================

test_that("TC-7.3.7: insufficient sample returns NA row", {
  # Create tiny dataset: 1 treated, 0 control at some period
  dt <- data.table::data.table(
    unit = c(1L, 1L, 1L, 1L, 1L),
    time = c(1L, 2L, 3L, 4L, 5L),
    y = c(1.0, 2.0, 3.0, 4.0, 5.0),
    d = c(1L, 1L, 1L, 1L, 1L)  # all treated, no control
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L
  )
  # Anchor should still exist
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  # Non-anchor rows should have NA att (insufficient control)
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})

# ============================================================================
# Section 12: TC-7.3.12 — n_total < 3 returns NA
# ============================================================================

test_that("TC-7.3.12: n_total < 3 returns NA row", {
  # 1 treated + 1 control = 2 < 3
  dt <- data.table::data.table(
    unit = rep(1:2, each = 5L),
    time = rep(1:5, 2L),
    y = rnorm(10),
    d = rep(c(1L, 0L), each = 5L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L
  )
  # Non-anchor rows should have NA (n_total=2 < 3)
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})


# ============================================================================
# Section 13: TC-7.3.10 — Estimator failure degrades to NA
# ============================================================================

test_that("TC-7.3.10: estimator failure degrades to NA row", {
  # Force estimator failure by requesting ipwra without controls
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 800L)
  # ipwra requires controls; passing NULL should cause failure
  # which gets caught by tryCatch and returns NA
  expect_warning(
    res <- estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 5L,
      estimator = "ipwra", controls = NULL
    )
  )
  # Should still return a data.frame (not error out)
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 14L)
  # Non-anchor rows should be NA (estimator failed)
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})


# ============================================================================
# Section 14: TC-7.3.13 — Demean symmetric transform verification
#   (AC-02/AC-15): uses {t+1,...,S-1} reference, NOT {1,...,S-1}
# ============================================================================

test_that("TC-7.3.13: demean uses symmetric {t+1,...,S-1} ref", {
  # Construct small panel: 3 units (2 treat, 1 ctrl), 6 periods
  # S=4 (tpost1=4), pre-periods: 1,2,3
  # T_min=1 excluded, anchor=3, estimable={2}
  # For period t=2, ref_periods = {3} (one period)
  dt <- data.table::data.table(
    unit = rep(1:3, each = 6L),
    time = rep(1:6, 3L),
    y = c(
      # Unit 1 (treated): 1,2,3,4,5,6
      1, 2, 3, 4, 5, 6,
      # Unit 2 (treated): 2,4,6,8,10,12
      2, 4, 6, 8, 10, 12,
      # Unit 3 (control): 1,3,5,7,9,11
      1, 3, 5, 7, 9, 11
    ),
    d = rep(c(1L, 1L, 0L), each = 6L)
  )

  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L, rolling = "demean"
  )

  # For t=2, ref_periods = {3}
  # Unit 1: y_trans = Y_{1,2} - Y_{1,3} = 2 - 3 = -1
  # Unit 2: y_trans = Y_{2,2} - Y_{2,3} = 4 - 6 = -2
  # Unit 3: y_trans = Y_{3,2} - Y_{3,3} = 3 - 5 = -2
  # ATT = mean_treat - mean_ctrl = (-1 + -2)/2 - (-2) = -1.5 + 2 = 0.5
  t2_row <- res[res$period == 2L, ]
  expect_equal(nrow(t2_row), 1L)
  expect_equal(t2_row$att, 0.5, tolerance = 1e-10)
})

# ============================================================================
# Section 15: TC-7.3.14 — Detrend symmetric transform verification
#   (AC-03/AC-16): OLS on {t+1,...,S-1}
# ============================================================================

test_that("TC-7.3.14: detrend OLS on symmetric ref periods", {
  # 3 units, 7 periods, S=5 (tpost1=5)
  # Pre-periods: 1,2,3,4. T_min=1 excluded, anchor=4
  # Estimable: {2, 3}
  # For t=2, ref_periods = {3, 4} (2 periods, can fit OLS)
  dt <- data.table::data.table(
    unit = rep(1:3, each = 7L),
    time = rep(1:7, 3L),
    y = c(
      # Unit 1 (treated): linear trend y = 1 + 0.5*t
      1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
      # Unit 2 (treated): y = 2 + 1*t
      3, 4, 5, 6, 7, 8, 9,
      # Unit 3 (control): y = 0 + 0.8*t
      0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 5.6
    ),
    d = rep(c(1L, 1L, 0L), each = 7L)
  )

  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, rolling = "detrend"
  )

  # For t=2, ref_periods = {3, 4}
  # Unit 1: OLS on (3,2.5),(4,3.0) => alpha=1, beta=0.5
  #   y_trans = Y_{1,2} - (1 + 0.5*2) = 2.0 - 2.0 = 0.0
  # Unit 2: OLS on (3,5),(4,6) => alpha=2, beta=1
  #   y_trans = Y_{2,2} - (2 + 1*2) = 4.0 - 4.0 = 0.0
  # Unit 3: OLS on (3,2.4),(4,3.2) => alpha=0, beta=0.8
  #   y_trans = Y_{3,2} - (0 + 0.8*2) = 1.6 - 1.6 = 0.0
  # ATT = 0 - 0 = 0 (perfect linear trends)
  t2_row <- res[res$period == 2L, ]
  expect_equal(nrow(t2_row), 1L)
  expect_equal(t2_row$att, 0.0, tolerance = 1e-10)

  # For t=3, ref_periods = {4} (only 1 period => degrade to demean)
  # Unit 1: y_trans = Y_{1,3} - Y_{1,4} = 2.5 - 3.0 = -0.5
  # Unit 2: y_trans = Y_{2,3} - Y_{2,4} = 5.0 - 6.0 = -1.0
  # Unit 3: y_trans = Y_{3,3} - Y_{3,4} = 2.4 - 3.2 = -0.8
  # ATT = (-0.5 + -1.0)/2 - (-0.8) = -0.75 + 0.8 = 0.05
  t3_row <- res[res$period == 3L, ]
  expect_equal(nrow(t3_row), 1L)
  expect_equal(t3_row$att, 0.05, tolerance = 1e-10)
})


# ============================================================================
# Section 16: TC-7.3.15 — Detrend degrades to demean when ref < 2
#   (AC-04)
# ============================================================================

test_that("TC-7.3.15: detrend degrades to demean with 1 ref period", {
  # 3 units, 4 periods, S=3 (tpost1=3)
  # Pre-periods: 1,2. T_min=1 excluded, anchor=2
  # No estimable periods (only anchor). But let's use S=4:
  # Pre-periods: 1,2,3. T_min=1 excluded, anchor=3, estimable={2}
  # For t=2, ref_periods = {3} (1 period => degrade to demean)
  dt <- data.table::data.table(
    unit = rep(1:3, each = 5L),
    time = rep(1:5, 3L),
    y = c(
      # Unit 1 (treated)
      10, 20, 30, 40, 50,
      # Unit 2 (treated)
      15, 25, 35, 45, 55,
      # Unit 3 (control)
      12, 22, 32, 42, 52
    ),
    d = rep(c(1L, 1L, 0L), each = 5L)
  )

  res_detrend <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L, rolling = "detrend"
  )
  res_demean <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L, rolling = "demean"
  )

  # For t=2, ref={3} (1 period): detrend should equal demean
  t2_detrend <- res_detrend[res_detrend$period == 2L, ]
  t2_demean <- res_demean[res_demean$period == 2L, ]
  expect_equal(t2_detrend$att, t2_demean$att, tolerance = 1e-10)
})

# ============================================================================
# Section 17: TC-7.3.16 — Function receives raw data, internal
#   transform matches manual calculation (AC-01)
# ============================================================================

test_that("TC-7.3.16: raw data input, transform matches manual", {
  # 4 units (2 treat, 2 ctrl), 6 periods, S=5
  # Pre: 1,2,3,4. T_min=1 excluded, anchor=4, estimable={2,3}
  dt <- data.table::data.table(
    unit = rep(1:4, each = 6L),
    time = rep(1:6, 4L),
    y = c(
      # Unit 1 (treated): 10,20,30,40,50,60
      10, 20, 30, 40, 50, 60,
      # Unit 2 (treated): 5,15,25,35,45,55
      5, 15, 25, 35, 45, 55,
      # Unit 3 (control): 8,18,28,38,48,58
      8, 18, 28, 38, 48, 58,
      # Unit 4 (control): 12,22,32,42,52,62
      12, 22, 32, 42, 52, 62
    ),
    d = rep(c(1L, 1L, 0L, 0L), each = 6L)
  )

  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, rolling = "demean"
  )

  # Manual calculation for t=2, ref_periods = {3, 4}
  # Unit 1: ref_mean = (30+40)/2 = 35, y_trans = 20 - 35 = -15
  # Unit 2: ref_mean = (25+35)/2 = 30, y_trans = 15 - 30 = -15
  # Unit 3: ref_mean = (28+38)/2 = 33, y_trans = 18 - 33 = -15
  # Unit 4: ref_mean = (32+42)/2 = 37, y_trans = 22 - 37 = -15
  # ATT = mean_treat - mean_ctrl = -15 - (-15) = 0
  t2_row <- res[res$period == 2L, ]
  expect_equal(t2_row$att, 0.0, tolerance = 1e-10)

  # Manual for t=3, ref_periods = {4}
  # Unit 1: y_trans = 30 - 40 = -10
  # Unit 2: y_trans = 25 - 35 = -10
  # Unit 3: y_trans = 28 - 38 = -10
  # Unit 4: y_trans = 32 - 42 = -10
  # ATT = -10 - (-10) = 0
  t3_row <- res[res$period == 3L, ]
  expect_equal(t3_row$att, 0.0, tolerance = 1e-10)
})


# ============================================================================
# Section 18: Return structure — 14 columns consistent
# ============================================================================

test_that("return structure has exactly 14 columns", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 900L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expected_cols <- c(
    "cohort", "period", "event_time", "att", "se",
    "ci_lower", "ci_upper", "t_stat", "pvalue",
    "n_treated", "n_control", "is_anchor",
    "rolling_window_size", "df_inference"
  )
  expect_equal(ncol(res), 14L)
  expect_equal(names(res), expected_cols)
})

test_that("cohort is always NA_integer_ (AC-10)", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 950L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expect_true(all(is.na(res$cohort)))
  expect_true(is.integer(res$cohort))
})

test_that("anchor row has correct values", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 960L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  expect_equal(anchor$att, 0.0)
  expect_equal(anchor$se, 0.0)
  expect_equal(anchor$ci_lower, 0.0)
  expect_equal(anchor$ci_upper, 0.0)
  expect_true(is.nan(anchor$t_stat))
  expect_true(is.nan(anchor$pvalue))
  expect_equal(anchor$rolling_window_size, 0L)
  expect_equal(anchor$event_time, -1L)
  expect_equal(anchor$period, 4L)  # tpost1 - 1
})

test_that("results sorted by event_time descending", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 970L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  # event_time should be in descending order
  expect_true(all(diff(res$event_time) <= 0))
  # First row should be anchor (event_time = -1)
  expect_true(res$is_anchor[1])
})


# ============================================================================
# Section 19: Edge Cases (Task 5)
# ============================================================================

# 5.1: Only 1 pre-treatment period (T_min = S-1, anchor only)
test_that("5.1: only 1 pre-period => anchor result (not NULL)", {
  # tpost1=2, periods={1,2,...}. Pre-periods={1}.
  # T_min=1 excluded => estimable empty, but anchor=1 (=tpost1-1)
  dt <- data.table::data.table(
    unit = rep(1:4, each = 3L),
    time = rep(1:3, 4L),
    y = rnorm(12),
    d = rep(c(1L, 1L, 0L, 0L), each = 3L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 2L
  )
  # Should return anchor (not NULL)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(res$is_anchor[1])
  expect_equal(res$att[1], 0.0)
})

# 5.2: No pre-treatment periods => NULL + warning
test_that("5.2: no pre-periods => NULL + warning", {
  # tpost1=1, all periods >= 1 => no pre-periods
  dt <- data.table::data.table(
    unit = rep(1:4, each = 3L),
    time = rep(1:3, 4L),
    y = rnorm(12),
    d = rep(c(1L, 1L, 0L, 0L), each = 3L)
  )
  expect_warning(
    res <- estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 1L
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_null(res)
})

# 5.3: exclude_pre_periods parameter
test_that("5.3: exclude_pre_periods reduces pre-period set", {
  dt <- make_pretreat_dgp(
    T_total = 8L, tpost1 = 5L, seed = 1100L
  )
  # Without exclusion: pre_periods = {1,2,3,4}
  res_full <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 0L
  )
  # With exclude_pre_periods=1: effective_S=4,
  # pre_periods = {1,2,3}, T_min=1 excluded, anchor=3,
  # estimable={2}
  res_excl <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 1L
  )
  # Excluded version should have fewer rows
  expect_true(nrow(res_excl) < nrow(res_full))
  # Excluded version anchor should be at period 3 (not 4)
  anchor_excl <- res_excl[res_excl$is_anchor, ]
  expect_equal(anchor_excl$period, 3L)
})

# 5.4: All pre-periods have insufficient sample
test_that("5.4: all pre-periods insufficient => all NA + anchor", {
  # 1 treated + 1 control = 2 < 3 at every period
  dt <- data.table::data.table(
    unit = rep(1:2, each = 6L),
    time = rep(1:6, 2L),
    y = rnorm(12),
    d = rep(c(1L, 0L), each = 6L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expect_s3_class(res, "data.frame")
  # Anchor should exist
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  expect_equal(anchor$att, 0.0)
  # Non-anchor rows should all be NA
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})

# 5.5: cohort column always NA_integer_
test_that("5.5: cohort always NA_integer_ for E7-05 compat", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 1200L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expect_true(all(is.na(res$cohort)))
  expect_true(is.integer(res$cohort))
  # Also check with detrend
  res_dt <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, rolling = "detrend"
  )
  expect_true(all(is.na(res_dt$cohort)))
  expect_true(is.integer(res_dt$cohort))
})


# ── Boundary case tests (Task 5) ─────────────────────────────────────────────

test_that("5.1: single pre-period (anchor only) returns non-NULL", {
  dt <- data.table::data.table(
    unit = rep(1:6, each = 2L),
    time = rep(1:2, 6L),
    y = rnorm(12),
    d = rep(c(1L, 1L, 1L, 0L, 0L, 0L), each = 2L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 2L
  )
  expect_false(is.null(res))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(res$is_anchor[1])
  expect_equal(res$att[1], 0.0)
  expect_equal(res$rolling_window_size[1], 0L)
})

test_that("5.2: no pre-treatment periods returns NULL with warning", {
  dt <- data.table::data.table(
    unit = rep(1:6, each = 3L),
    time = rep(1:3, 6L),
    y = rnorm(18),
    d = rep(c(1L, 1L, 1L, 0L, 0L, 0L), each = 3L)
  )
  expect_warning(
    res <- estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 1L
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_null(res)
})

test_that("5.3: exclude_pre_periods reduces pre-period set", {
  dt <- make_pretreat_dgp(
    n_treat = 10L, n_ctrl = 10L, T_total = 8L,
    tpost1 = 5L, seed = 901L
  )
  res_full <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 0L
  )
  res_excl <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 2L
  )
  # Excluded version should have fewer rows
  expect_true(nrow(res_excl) < nrow(res_full))
  # Both should be valid data.frames with 14 columns
  expect_equal(ncol(res_full), 14L)
  expect_equal(ncol(res_excl), 14L)
})

test_that("5.4: all periods insufficient sample returns NA df", {
  dt <- data.table::data.table(
    unit = rep(1:3, each = 5L),
    time = rep(1:5, 3L),
    y = rnorm(15),
    d = rep(1L, 15L)  # all treated, no control
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L
  )
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 14L)
  # Anchor should exist
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  # Non-anchor rows should have NA att
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})

test_that("5.5: cohort column is always NA_integer_", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 902L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expect_true(all(is.na(res$cohort)))
  expect_true(is.integer(res$cohort))
})


# ============================================================================
# Section 15: Boundary Cases (Task 5)
# ============================================================================

# --- 5.1: Single pre-period (anchor only) returns non-NULL (AC-12) ---
test_that("5.1: single pre-period (anchor only) returns non-NULL", {
  dt <- data.table::data.table(
    unit = rep(1:6, each = 2L),
    time = rep(1:2, 6L),
    y = rnorm(12),
    d = rep(c(1L, 1L, 1L, 0L, 0L, 0L), each = 2L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 2L
  )
  expect_false(is.null(res))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(res$is_anchor[1])
  expect_equal(res$att[1], 0.0)
  expect_equal(res$rolling_window_size[1], 0L)
})

# --- 5.2: No pre-treatment periods → NULL + warning (AC-19) ---
test_that("5.2: no pre-treatment periods returns NULL with warning", {
  dt <- data.table::data.table(
    unit = rep(1:6, each = 3L),
    time = rep(1:3, 6L),
    y = rnorm(18),
    d = rep(c(1L, 1L, 1L, 0L, 0L, 0L), each = 3L)
  )
  expect_warning(
    res <- estimate_pre_treatment_common(
      dt, y = "y", ivar = "unit", d = "d",
      tvar = "time", tpost1 = 1L
    ),
    class = "lwdid_insufficient_pre_periods"
  )
  expect_null(res)
})

# --- 5.3: exclude_pre_periods reduces pre-period set (REQ-16) ---
test_that("5.3: exclude_pre_periods reduces pre-period set", {
  dt <- make_pretreat_dgp(
    n_treat = 10L, n_ctrl = 10L, T_total = 8L,
    tpost1 = 5L, seed = 901L
  )
  res_full <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 0L
  )
  res_excl <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L, exclude_pre_periods = 2L
  )
  # Excluded version should have fewer rows
  expect_true(nrow(res_excl) < nrow(res_full))
  # Both should be valid data.frames with 14 columns
  expect_equal(ncol(res_full), 14L)
  expect_equal(ncol(res_excl), 14L)
})

# --- 5.4: All periods insufficient sample → NA df (AC-20) ---
test_that("5.4: all periods insufficient sample returns NA df", {
  dt <- data.table::data.table(
    unit = rep(1:3, each = 5L),
    time = rep(1:5, 3L),
    y = rnorm(15),
    d = rep(1L, 15L)
  )
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 4L
  )
  expect_s3_class(res, "data.frame")
  expect_equal(ncol(res), 14L)
  anchor <- res[res$is_anchor, ]
  expect_equal(nrow(anchor), 1L)
  non_anchor <- res[!res$is_anchor, ]
  if (nrow(non_anchor) > 0L) {
    expect_true(all(is.na(non_anchor$att)))
  }
})

# --- 5.5: cohort column is always NA_integer_ (AC-10) ---
test_that("5.5: cohort column is always NA_integer_", {
  dt <- make_pretreat_dgp(tpost1 = 5L, seed = 902L)
  res <- estimate_pre_treatment_common(
    dt, y = "y", ivar = "unit", d = "d",
    tvar = "time", tpost1 = 5L
  )
  expect_true(all(is.na(res$cohort)))
  expect_true(is.integer(res$cohort))
})
