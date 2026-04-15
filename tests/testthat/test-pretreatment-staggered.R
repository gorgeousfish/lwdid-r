# test-pretreatment-staggered.R
# Unit tests for estimate_pre_treatment_staggered()
# Story E7-04: Pre-treatment effect estimation (Staggered)
# Paper: lw2025 Section 5, Procedure 5.1

library(testthat)
library(data.table)

# ===== Test DGP helper =====
# Creates a staggered panel with known properties
# Units: N_per_cohort per cohort, plus N_nt never-treated
# Cohorts: g1, g2 (treatment times)
# Periods: 1..T_max
# Y = unit_fe + time_fe + treatment_effect * post
make_staggered_panel <- function(
    cohorts = c(5L, 7L),
    n_per_cohort = 10L,
    n_nt = 10L,
    T_max = 10L,
    unit_fe_sd = 2,
    time_trend = 1.5,
    treatment_effect = 0,
    seed = 42L) {
  set.seed(seed)
  n_cohorts <- length(cohorts)
  n_total <- n_cohorts * n_per_cohort + n_nt

  # Unit IDs and cohort assignments
  ids <- seq_len(n_total)
  gvar_vals <- c(
    rep(cohorts, each = n_per_cohort),
    rep(Inf, n_nt)
  )

  # Build panel
  dt <- CJ(id = ids, time = seq_len(T_max))
  dt[, gvar := gvar_vals[id]]
  dt[, d := as.integer(!is.na(gvar) & gvar != Inf &
                          gvar != 0 & time >= gvar)]

  # Y = unit_fe + time_trend * time + treatment_effect * d + noise
  unit_fes <- rnorm(n_total, sd = unit_fe_sd)
  dt[, y := unit_fes[id] + time_trend * time +
       treatment_effect * d + rnorm(.N, sd = 0.5)]
  dt
}

# ===== TC-7.4.1: Two cohorts, demean, ra =====
test_that("TC-7.4.1: two cohorts pre-period ATT correct (demean, ra)", {
  dt <- make_staggered_panel(
    cohorts = c(5L, 8L), n_per_cohort = 15L, n_nt = 15L,
    T_max = 10L, treatment_effect = 0, seed = 101
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("cohort", "period", "event_time", "att", "se",
    "is_anchor", "rolling_window_size") %in% names(result)))
  # Under PT (treatment_effect=0), pre-period ATTs should be near 0
  non_anchor <- result[!result$is_anchor, ]
  expect_true(all(abs(non_anchor$att) < 5),
    info = "Pre-period ATTs should be small under PT")
  # Both cohorts present
  expect_true(5L %in% result$cohort)
  expect_true(8L %in% result$cohort)
})


# ===== TC-7.4.2: Symmetric transform ref periods correct =====
test_that("TC-7.4.2: ref_periods = {t+1,...,g-1} verified", {
  # Small deterministic panel: 3 treated (g=5), 3 NT
  dt <- data.table(
    id = rep(1:6, each = 8),
    time = rep(1:8, 6),
    gvar = rep(c(5L, 5L, 5L, Inf, Inf, Inf), each = 8)
  )
  # Y = id*10 + time*3 (deterministic, no noise)
  dt[, y := id * 10 + time * 3]

  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  expect_false(is.null(result))
  # Cohort 5: anchor=4, T_min=1 excluded, estimable={2,3}
  # Period 2: ref={3,4}, rws=2
  # Period 3: ref={4}, rws=1
  r2 <- result[result$period == 2 & result$cohort == 5, ]
  r3 <- result[result$period == 3 & result$cohort == 5, ]
  expect_equal(r2$rolling_window_size, 2L)
  expect_equal(r3$rolling_window_size, 1L)
})

# ===== TC-7.4.3: Anchor properties =====
test_that("TC-7.4.3: anchor ATT=0, SE=0, is_anchor=TRUE", {
  dt <- make_staggered_panel(cohorts = c(5L, 8L), seed = 201)
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  anchors <- result[result$is_anchor == TRUE, ]
  expect_equal(nrow(anchors), 2L)  # one per cohort
  expect_true(all(anchors$att == 0))
  expect_true(all(anchors$se == 0))
  expect_true(all(is.nan(anchors$t_stat)))
  expect_true(all(is.nan(anchors$pvalue)))
  expect_true(all(anchors$rolling_window_size == 0L))
  # Anchor event_time = -1
  expect_true(all(anchors$event_time == -1L))
})

# ===== TC-7.4.4: Control group based on t (not r) =====
test_that("TC-7.4.4: n_control varies with t_pre (not_yet_treated)", {
  # Cohorts 4, 7: at t=2, cohort 4 units are not-yet-treated controls for cohort 7
  # At t=5, cohort 4 units are already treated → excluded from controls for cohort 7
  dt <- make_staggered_panel(
    cohorts = c(4L, 7L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 301
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "not_yet_treated"
  )
  # For cohort 7: earlier pre-periods should have more controls
  c7 <- result[result$cohort == 7 & !result$is_anchor, ]
  if (nrow(c7) >= 2) {
    # Earlier periods (smaller t) → more not-yet-treated → more controls
    earliest <- c7[which.min(c7$period), ]
    latest <- c7[which.max(c7$period), ]
    expect_true(earliest$n_control >= latest$n_control,
      info = "Earlier pre-periods should have >= controls")
  }
})

# ===== TC-7.4.5: never_treated control group =====
test_that("TC-7.4.5: never_treated uses only NT units as controls", {
  dt <- make_staggered_panel(
    cohorts = c(5L, 8L), n_per_cohort = 10L, n_nt = 15L,
    T_max = 10L, seed = 401
  )
  result_nt <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "never_treated"
  )
  # All non-anchor rows should have n_control = 15 (the NT units)
  non_anchor <- result_nt[!result_nt$is_anchor, ]
  expect_true(all(non_anchor$n_control == 15L),
    info = "never_treated should use exactly 15 NT units")
})

# ===== TC-7.4.6: FATAL-001 strict inequality G_i > t =====
test_that("TC-7.4.6: not_yet_treated uses G_i > t strict inequality", {
  # Cohort 5 at t=5: units with G_i=5 should NOT be controls (G_i > 5 fails)
  # But at t=4: units with G_i=5 have G_i > 4, so they ARE treated (d=1), not controls
  dt <- make_staggered_panel(
    cohorts = c(5L, 8L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 501
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "not_yet_treated"
  )
  # For cohort 8, at pre-period t=3: G_i=5 units have G_i > 3 → valid controls
  # For cohort 8, at pre-period t=5: G_i=5 units have G_i > 5 FALSE → excluded
  c8 <- result[result$cohort == 8 & !result$is_anchor, ]
  if (nrow(c8) >= 2) {
    t3_row <- c8[c8$period == 3, ]
    t5_row <- c8[c8$period == 5, ]
    if (nrow(t3_row) > 0 && nrow(t5_row) > 0) {
      # t=3 should have more controls than t=5
      expect_true(t3_row$n_control > t5_row$n_control)
    }
  }
})


# ===== TC-7.4.8: Only 1 pre-period (anchor only) =====
test_that("TC-7.4.8: cohort with 1 pre-period returns anchor only", {
  # Cohort g=2: pre-periods={1}, anchor=1, T_min=1 excluded → no estimable
  dt <- make_staggered_panel(
    cohorts = c(2L, 6L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 8L, seed = 601
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  # Cohort 2 should only have anchor row
  c2 <- result[result$cohort == 2, ]
  expect_equal(nrow(c2), 1L)
  expect_true(c2$is_anchor)
  expect_equal(c2$att, 0)
})

# ===== TC-7.4.9: ipwra estimator routes correctly =====
test_that("TC-7.4.9: estimator='ipwra' routes correctly", {
  dt <- make_staggered_panel(
    cohorts = c(5L), n_per_cohort = 20L, n_nt = 20L,
    T_max = 8L, seed = 701
  )
  # Add a control variable
  dt[, x1 := rnorm(.N)]
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ipwra", controls = "x1"
  )
  expect_false(is.null(result))
  non_anchor <- result[!result$is_anchor, ]
  # Should have some non-NA ATT estimates
  expect_true(any(!is.na(non_anchor$att)))
})

# ===== TC-7.4.10: Result sorting =====
test_that("TC-7.4.10: sorted by cohort asc, event_time desc", {
  dt <- make_staggered_panel(
    cohorts = c(4L, 7L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 801
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  # Cohort should be non-decreasing
  expect_true(all(diff(result$cohort) >= 0))
  # Within each cohort, event_time should be non-increasing
  for (g in unique(result$cohort)) {
    cg <- result[result$cohort == g, ]
    expect_true(all(diff(cg$event_time) <= 0),
      info = sprintf("Cohort %d event_time not descending", g))
  }
})

# ===== TC-7.4.11: rolling_window_size correct =====
test_that("TC-7.4.11: rolling_window_size = g-t-1, anchor=0", {
  dt <- make_staggered_panel(
    cohorts = c(6L), n_per_cohort = 15L, n_nt = 15L,
    T_max = 8L, seed = 901
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  for (i in seq_len(nrow(result))) {
    row <- result[i, ]
    if (row$is_anchor) {
      expect_equal(row$rolling_window_size, 0L)
    } else if (!is.na(row$rolling_window_size)) {
      expected_rws <- as.integer(row$cohort - row$period - 1L)
      expect_equal(row$rolling_window_size, expected_rws)
    }
  }
})

# ===== TC-7.4.12: Anchor n_control via control group logic =====
test_that("TC-7.4.12: anchor n_control computed via control group", {
  dt <- make_staggered_panel(
    cohorts = c(5L), n_per_cohort = 10L, n_nt = 12L,
    T_max = 8L, seed = 1001
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "never_treated"
  )
  anchor <- result[result$is_anchor, ]
  # With never_treated control, n_control should be 12
  expect_equal(anchor$n_control, 12L)
  expect_equal(anchor$n_treated, 10L)
})

# ===== TC-7.4.13: n_total < 3 returns NA row =====
test_that("TC-7.4.13: insufficient sample returns NA row", {
  # Create tiny panel: 1 treated, 1 control
  dt <- data.table(
    id = rep(1:2, each = 6),
    time = rep(1:6, 2),
    gvar = rep(c(5L, Inf), each = 6),
    y = rnorm(12)
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  # With only 2 units total (1 treat + 1 ctrl), n_total=2 < 3
  non_anchor <- result[!result$is_anchor, ]
  if (nrow(non_anchor) > 0) {
    expect_true(all(is.na(non_anchor$att)),
      info = "ATT should be NA when n_total < 3")
  }
})

# ===== TC-7.4.14: Estimator failure degrades to NA =====
test_that("TC-7.4.14: estimator failure produces NA row + warning", {
  dt <- make_staggered_panel(
    cohorts = c(5L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 8L, seed = 1101
  )
  # Use PSM with a non-existent control variable to force failure
  expect_warning(
    result <- estimate_pre_treatment_staggered(
      data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
      rolling = "demean", estimator = "psm", controls = "nonexistent_var"
    ),
    class = "lwdid_data"
  )
})

# ===== TC-7.4.15: NA row has NA_integer_ for rws and df =====
test_that("TC-7.4.15: NA rows have NA_integer_ rws and df_inference", {
  dt <- data.table(
    id = rep(1:2, each = 6),
    time = rep(1:6, 2),
    gvar = rep(c(5L, Inf), each = 6),
    y = rnorm(12)
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  na_rows <- result[is.na(result$att) & !result$is_anchor, ]
  if (nrow(na_rows) > 0) {
    expect_true(all(is.na(na_rows$rolling_window_size)))
    expect_true(all(is.na(na_rows$df_inference)))
  }
})


# ===== TC-7.4.16: Detrend differs from demean =====
test_that("TC-7.4.16: detrend uses linear trend, differs from demean", {
  dt <- make_staggered_panel(
    cohorts = c(6L), n_per_cohort = 15L, n_nt = 15L,
    T_max = 8L, seed = 1201
  )
  res_dm <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  res_dt <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "detrend", estimator = "ra"
  )
  # For periods with >= 2 ref periods, detrend should differ from demean
  # Cohort 6: period 2 has ref={3,4,5}, rws=3 → full detrend
  dm_t2 <- res_dm[res_dm$period == 2 & res_dm$cohort == 6, ]
  dt_t2 <- res_dt[res_dt$period == 2 & res_dt$cohort == 6, ]
  if (nrow(dm_t2) > 0 && nrow(dt_t2) > 0 &&
      !is.na(dm_t2$att) && !is.na(dt_t2$att)) {
    expect_false(isTRUE(all.equal(dm_t2$att, dt_t2$att, tolerance = 1e-6)),
      info = "Detrend ATT should differ from demean ATT")
  }
})

# ===== TC-7.4.17: Detrend degrades to demean with 1 ref =====
test_that("TC-7.4.17: detrend with 1 ref period degrades to demean", {
  dt <- make_staggered_panel(
    cohorts = c(6L), n_per_cohort = 15L, n_nt = 15L,
    T_max = 8L, seed = 1301
  )
  res_dm <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  res_dt <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "detrend", estimator = "ra"
  )
  # Period 4 for cohort 6: ref={5}, rws=1 → detrend degrades to demean
  dm_t4 <- res_dm[res_dm$period == 4 & res_dm$cohort == 6, ]
  dt_t4 <- res_dt[res_dt$period == 4 & res_dt$cohort == 6, ]
  if (nrow(dm_t4) > 0 && nrow(dt_t4) > 0 &&
      !is.na(dm_t4$att) && !is.na(dt_t4$att)) {
    expect_equal(dm_t4$att, dt_t4$att, tolerance = 1e-10,
      info = "Detrend with 1 ref should equal demean")
  }
})

# ===== TC-7.4.18: Detrend degradation numerical match =====
test_that("TC-7.4.18: detrend degradation matches demean exactly", {
  # Same as TC-7.4.17 but also check SE
  dt <- make_staggered_panel(
    cohorts = c(6L), n_per_cohort = 15L, n_nt = 15L,
    T_max = 8L, seed = 1301
  )
  res_dm <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  res_dt <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "detrend", estimator = "ra"
  )
  dm_t4 <- res_dm[res_dm$period == 4 & res_dm$cohort == 6, ]
  dt_t4 <- res_dt[res_dt$period == 4 & res_dt$cohort == 6, ]
  if (nrow(dm_t4) > 0 && nrow(dt_t4) > 0) {
    expect_equal(dm_t4$se, dt_t4$se, tolerance = 1e-10)
  }
})

# ===== TC-7.4.19: Missing required columns =====
test_that("TC-7.4.19: missing columns triggers stop_lwdid", {
  dt <- data.table(id = 1:5, time = 1:5, y = rnorm(5))
  expect_error(
    estimate_pre_treatment_staggered(
      data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar"
    ),
    class = "lwdid_invalid_param"
  )
})

# ===== TC-7.4.20: vce='cluster' without cluster_var =====
test_that("TC-7.4.20: vce='cluster' without cluster_var errors", {
  dt <- make_staggered_panel(cohorts = c(5L), seed = 1401)
  expect_error(
    estimate_pre_treatment_staggered(
      data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
      vce = "cluster", cluster_var = NULL
    ),
    class = "lwdid_invalid_param"
  )
})

# ===== TC-7.4.21: All cohorts no estimable periods → NULL =====
test_that("TC-7.4.21: no estimable periods returns result with anchors only", {
  # Cohort g=2: only pre-period is 1 (T_min), excluded → anchor only
  dt <- data.table(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    gvar = rep(c(rep(2L, 5), rep(Inf, 5)), each = 3),
    y = rnorm(30)
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  # Should have anchor only (or NULL if no anchor generated)
  if (!is.null(result)) {
    expect_true(all(result$is_anchor))
  }
})

# ===== TC-7.4.23: Three+ cohorts independent pre-period sets =====
test_that("TC-7.4.23: multiple cohorts have independent pre-periods", {
  dt <- make_staggered_panel(
    cohorts = c(4L, 6L, 8L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 1501
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  # Each cohort should have its own set of periods
  for (g in c(4L, 6L, 8L)) {
    cg <- result[result$cohort == g, ]
    expect_true(nrow(cg) > 0, info = sprintf("Cohort %d missing", g))
    # All periods should be < g
    expect_true(all(cg$period < g))
    # Anchor should be at g-1
    anchor <- cg[cg$is_anchor, ]
    expect_equal(anchor$period, g - 1L)
  }
})

# ===== TC-7.4.24: Different cohorts have dynamic control groups =====
test_that("TC-7.4.24: not_yet_treated controls vary by cohort and t", {
  dt <- make_staggered_panel(
    cohorts = c(4L, 7L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 1601
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "not_yet_treated"
  )
  # Cohort 4 at t=2: controls include cohort 7 (G=7>2) + NT
  # Cohort 7 at t=2: controls include cohort 4 (G=4>2) + NT
  c4_t2 <- result[result$cohort == 4 & result$period == 2, ]
  c7_t2 <- result[result$cohort == 7 & result$period == 2, ]
  if (nrow(c4_t2) > 0 && nrow(c7_t2) > 0) {
    # Both should have controls (cohort 4 has more: 10 c7 + 10 NT = 20)
    expect_true(c4_t2$n_control > 0)
    expect_true(c7_t2$n_control > 0)
  }
})

# ===== TC-7.4.25: Anchor df_inference =====
test_that("TC-7.4.25: anchor df_inference = n_treat+n_ctrl-2", {
  dt <- make_staggered_panel(
    cohorts = c(5L), n_per_cohort = 10L, n_nt = 12L,
    T_max = 8L, seed = 1701
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra", control_group = "never_treated"
  )
  anchor <- result[result$is_anchor, ]
  expected_df <- anchor$n_treated + anchor$n_control - 2L
  expect_equal(anchor$df_inference, expected_df)
})

# ===== TC-7.4.26: Anchor event_time always -1 =====
test_that("TC-7.4.26: anchor event_time is always -1", {
  dt <- make_staggered_panel(
    cohorts = c(3L, 5L, 8L), n_per_cohort = 10L, n_nt = 10L,
    T_max = 10L, seed = 1801
  )
  result <- estimate_pre_treatment_staggered(
    data = dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
    rolling = "demean", estimator = "ra"
  )
  anchors <- result[result$is_anchor, ]
  expect_true(all(anchors$event_time == -1L))
})
