# ============================================================================
# test-transform.R — E2-06 Layer 1: Transform Correctness Tests
#
# Test Groups:
#   Layer 1 (T1-01 to T1-06): Transform correctness
#   Layer 4 (T4-02, T4-03, T4-07): Transform-related edge cases
# ============================================================================

# Helper: capture warnings while preserving result
capture_with_warnings <- function(expr) {
  warnings_caught <- list()
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings_caught)
}

# Helper: check if any warning has a given class
has_warning_class <- function(warnings_list, cls) {
  any(vapply(warnings_list, function(w) inherits(w, cls), logical(1)))
}

# ============================================================================
# Layer 1: Transform Correctness Tests
# ============================================================================

# --- T1-01: Demean transform basic mathematical correctness ---
test_that("T1-01: demean transform produces correct y_trans values", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    t  = rep(1:5, 3),
    y  = c(1, 2, 3, 4, 5,
           2, 3, 4, 5, 6,
           3, 4, 5, 6, 7),
    d    = rep(c(1, 0, 0), each = 5),
    post = rep(c(0, 0, 1, 1, 1), 3)
  )

  result <- transform_demean(dt, "y", "id", "t", g = 3L)

  # Unit 1: pre_mean = (1+2)/2 = 1.5
  expect_equal(result[id == 1, pre_mean][1], 1.5, tolerance = 1e-15)
  expect_equal(result[id == 1 & t == 3, y_trans], 3 - 1.5, tolerance = 1e-15)
  expect_equal(result[id == 1 & t == 4, y_trans], 4 - 1.5, tolerance = 1e-15)
  expect_equal(result[id == 1 & t == 5, y_trans], 5 - 1.5, tolerance = 1e-15)

  # Unit 2: pre_mean = (2+3)/2 = 2.5
  expect_equal(result[id == 2, pre_mean][1], 2.5, tolerance = 1e-15)
  expect_equal(result[id == 2 & t == 3, y_trans], 4 - 2.5, tolerance = 1e-15)
  expect_equal(result[id == 2 & t == 4, y_trans], 5 - 2.5, tolerance = 1e-15)
  expect_equal(result[id == 2 & t == 5, y_trans], 6 - 2.5, tolerance = 1e-15)

  # Unit 3: pre_mean = (3+4)/2 = 3.5
  expect_equal(result[id == 3, pre_mean][1], 3.5, tolerance = 1e-15)
  expect_equal(result[id == 3 & t == 3, y_trans], 5 - 3.5, tolerance = 1e-15)
  expect_equal(result[id == 3 & t == 4, y_trans], 6 - 3.5, tolerance = 1e-15)
  expect_equal(result[id == 3 & t == 5, y_trans], 7 - 3.5, tolerance = 1e-15)

  # n_pre should be 2 for all units (t=1,2 are pre-treatment)
  expect_true(all(result[["n_pre"]] == 2L))

  # y_trans = y - pre_mean for ALL rows
  expect_equal(result[["y_trans"]], result[["y"]] - result[["pre_mean"]],
               tolerance = 1e-15)
})


# --- T1-02: Perfect linear trend detrend residuals are zero ---
test_that("T1-02: detrend with perfect linear trend Y=1+2t yields zero post-period residuals", {
  # Y = 1 + 2*t for t=1,...,6, g=4 (pre: t=1,2,3; post: t=4,5,6)
  dt <- data.table::data.table(
    id = rep(1:2, each = 6),
    t  = rep(1:6, 2),
    y  = rep(1 + 2 * (1:6), 2)  # Y = 3,5,7,9,11,13
  )

  result <- suppressMessages(
    transform_detrend(dt, "y", "id", "t", g = 4L)
  )

  # Pre-period: t=1,2,3 -> Y=3,5,7
  # t_bar_pre = (1+2+3)/3 = 2
  # OLS on centered time: Y = A* + B*(t - 2)
  # A* = mean(Y_pre) = 5, B = 2 (perfect fit)
  # Post-period predictions:
  #   t=4: 5 + 2*(4-2) = 9, residual = 9-9 = 0
  #   t=5: 5 + 2*(5-2) = 11, residual = 11-11 = 0
  #   t=6: 5 + 2*(6-2) = 13, residual = 13-13 = 0

  for (uid in unique(result[["id"]])) {
    post_resid <- result[id == uid & t >= 4, y_trans]
    expect_equal(post_resid, rep(0, 3), tolerance = 1e-10)
  }

  # Verify slope = 2 and intercept_c = 5
  expect_equal(result[id == 1, slope][1], 2, tolerance = 1e-10)
  expect_equal(result[id == 1, intercept_c][1], 5, tolerance = 1e-10)
  expect_equal(result[id == 1, t_bar_pre][1], 2, tolerance = 1e-10)
})

# --- T1-03: Individual unit with <2 pre-period obs degrades to demean ---
test_that("T1-03: unit with only 1 pre-period obs degrades to demean with warning", {
  # 3 units, g=4 (pre: t<=3). Unit 3 only has t=3 in pre-period (1 obs).
  dt <- data.table::data.table(
    id = c(rep(1, 6), rep(2, 6), 3, 3, 3, 3),
    t  = c(1:6, 1:6, 3, 4, 5, 6),
    y  = c(1, 2, 3, 4, 5, 6,
           2, 3, 4, 5, 6, 7,
           10, 11, 12, 13)
  )

  cw <- capture_with_warnings(
    suppressMessages(transform_detrend(dt, "y", "id", "t", g = 4L))
  )
  result <- cw$result

  # Should have lwdid_data warning about degradation
  expect_true(has_warning_class(cw$warnings, "lwdid_data"))

  # Unit 3 should be degraded: slope = 0
  u3 <- result[result[["id"]] == 3]
  expect_equal(u3[["slope"]][1], 0)

  # Unit 3 degraded to demean: y_trans = Y - pre_mean
  # pre_mean = Y at t=3 = 10
  u3_t4 <- result[result[["id"]] == 3 & result[["t"]] == 4]
  u3_t5 <- result[result[["id"]] == 3 & result[["t"]] == 5]
  u3_t6 <- result[result[["id"]] == 3 & result[["t"]] == 6]
  expect_equal(u3_t4[["y_trans"]], 11 - 10, tolerance = 1e-15)
  expect_equal(u3_t5[["y_trans"]], 12 - 10, tolerance = 1e-15)
  expect_equal(u3_t6[["y_trans"]], 13 - 10, tolerance = 1e-15)

  # Units 1 and 2 should NOT be degraded (3 pre-period obs each)
  u1 <- result[result[["id"]] == 1]
  u2 <- result[result[["id"]] == 2]
  expect_true(u1[["slope"]][1] != 0)
  expect_true(u2[["slope"]][1] != 0)
})

# --- T1-04: Exactly 2 pre-period obs triggers exact_fit warning ---
test_that("T1-04: exactly 2 pre-period obs triggers lwdid_small_sample warning", {
  # 1 unit x 4 periods, g=3 (pre: t=1,2; post: t=3,4)
  # Data chosen so A*=3, B=2, t_bar=1.5 (per task spec)
  dt <- data.table::data.table(
    id = rep(1, 4),
    t  = 1:4,
    y  = c(2, 4, 7, 9)
  )

  cw <- capture_with_warnings(
    suppressMessages(transform_detrend(dt, "y", "id", "t", g = 3L))
  )
  result <- cw$result

  # Should have lwdid_small_sample warning
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))

  # Pre: t=1,2 -> Y=2,4
  # t_bar = (1+2)/2 = 1.5
  # Centered: t*=(-0.5, 0.5), Y=(2,4)
  # A* = mean(2,4) = 3
  # B = sum(t* * Y) / sum(t*^2) = ((-0.5)*2+(0.5)*4)/0.5 = 1/0.5 = 2
  u1 <- result[result[["id"]] == 1]
  expect_equal(u1[["intercept_c"]][1], 3, tolerance = 1e-10)
  expect_equal(u1[["slope"]][1], 2, tolerance = 1e-10)
  expect_equal(u1[["t_bar_pre"]][1], 1.5, tolerance = 1e-10)

  # Post predictions:
  # t=3: predicted = 3 + 2*(3-1.5) = 6, y_trans = 7 - 6 = 1
  # t=4: predicted = 3 + 2*(4-1.5) = 8, y_trans = 9 - 8 = 1
  u1_t3 <- result[result[["id"]] == 1 & result[["t"]] == 3]
  u1_t4 <- result[result[["id"]] == 1 & result[["t"]] == 4]
  expect_equal(u1_t3[["y_trans"]], 1, tolerance = 1e-10)
  expect_equal(u1_t4[["y_trans"]], 1, tolerance = 1e-10)
})

# --- T1-05: Degenerate time variance degrades to demean ---
test_that("T1-05: degenerate time variance (all same t) degrades to demean", {
  # Unit 2 has 2 pre-period rows both at t=3 (var(t)=0)
  # g=4, pre is t<=3. Unit 1: t=1,2,3 (normal). Unit 2: t=3,3 (degenerate)
  dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 1, 2, 2, 2, 2),
    t  = c(1, 2, 3, 4, 5, 3, 3, 4, 5),
    y  = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  )

  cw <- capture_with_warnings(
    suppressMessages(transform_detrend(dt, "y", "id", "t", g = 4L))
  )
  result <- cw$result

  # Should have lwdid_data warning about degradation
  expect_true(has_warning_class(cw$warnings, "lwdid_data"))

  # Unit 2 should be degraded: slope = 0
  u2 <- result[result[["id"]] == 2]
  expect_equal(u2[["slope"]][1], 0)

  # Unit 2 degraded to demean: y_trans = Y - pre_mean
  # pre_mean = mean(6, 7) = 6.5
  pre_mean_u2 <- mean(c(6, 7))
  u2_t4 <- result[result[["id"]] == 2 & result[["t"]] == 4]
  expect_equal(
    u2_t4[["y_trans"]][1], 8 - pre_mean_u2, tolerance = 1e-15
  )
  u2_t5 <- result[result[["id"]] == 2 & result[["t"]] == 5]
  expect_equal(
    u2_t5[["y_trans"]][1], 9 - pre_mean_u2, tolerance = 1e-15
  )
})

# --- T1-06: Panel-level <2 pre-treatment time points -> error ---
test_that("T1-06: panel with <2 pre-treatment time points throws lwdid_insufficient_pre_periods", {
  # 2 units x 2 periods, g=2 (pre: only t=1, one time point)
  dt <- data.table::data.table(
    id = rep(1:2, each = 2),
    t  = rep(1:2, 2),
    y  = c(1, 2, 3, 4)
  )

  expect_error(
    transform_detrend(dt, "y", "id", "t", g = 2L),
    class = "lwdid_insufficient_pre_periods"
  )
})


# ============================================================================
# Layer 4: Transform-related Edge Cases
# ============================================================================

# --- T4-02: Unbalanced panel missing pre-period ---
test_that("T4-02: unbalanced panel unit missing t=1 uses only available pre-period", {
  # Unit 1 missing t=1, only has t=2 in pre-period (g=3, pre: t=1,2)
  # Unit 2 has both t=1,2
  dt <- data.table::data.table(
    id = c(1, 1, 1, 2, 2, 2, 2),
    t  = c(2, 3, 4, 1, 2, 3, 4),
    y  = c(5, 7, 9, 1, 2, 3, 4)
  )

  result <- transform_demean(dt, "y", "id", "t", g = 3L)

  # Unit 1: pre = {t=2, y=5}, pre_mean = 5
  u1 <- result[result[["id"]] == 1]
  expect_equal(u1[["pre_mean"]][1], 5, tolerance = 1e-15)
  expect_equal(u1[["n_pre"]][1], 1L)

  # Unit 1 post: y_trans = y - 5
  u1_t3 <- result[result[["id"]] == 1 & result[["t"]] == 3]
  expect_equal(u1_t3[["y_trans"]], 7 - 5, tolerance = 1e-15)

  # Unit 2: pre = {t=1,2, y=1,2}, pre_mean = 1.5
  u2 <- result[result[["id"]] == 2]
  expect_equal(u2[["pre_mean"]][1], 1.5, tolerance = 1e-15)
  expect_equal(u2[["n_pre"]][1], 2L)
})

# --- T4-03: exclude_pre_periods=1 excludes last pre-period ---
test_that("T4-03: exclude_pre_periods=1 excludes period g-1 from pre-treatment", {
  # g=4, pre normally t=1,2,3. With exclude_pre_periods=1, pre = t=1,2 (exclude t=3)
  dt <- data.table::data.table(
    id = rep(1:2, each = 5),
    t  = rep(1:5, 2),
    y  = c(1, 2, 3, 4, 5,
           2, 4, 6, 8, 10)
  )

  result <- transform_demean(dt, "y", "id", "t", g = 4L,
                             exclude_pre_periods = 1L)

  # Unit 1: pre = {t=1,2}, pre_mean = (1+2)/2 = 1.5
  u1 <- result[result[["id"]] == 1]
  expect_equal(u1[["pre_mean"]][1], 1.5, tolerance = 1e-15)
  expect_equal(u1[["n_pre"]][1], 2L)

  # Unit 2: pre = {t=1,2}, pre_mean = (2+4)/2 = 3
  u2 <- result[result[["id"]] == 2]
  expect_equal(u2[["pre_mean"]][1], 3, tolerance = 1e-15)
  expect_equal(u2[["n_pre"]][1], 2L)

  # Post-period y_trans
  u1_t4 <- result[result[["id"]] == 1 & result[["t"]] == 4]
  expect_equal(u1_t4[["y_trans"]], 4 - 1.5, tolerance = 1e-15)
})

# --- T4-07: NA filtering through transform to estimator ---
test_that("T4-07: unit with no pre-period obs gets NA y_trans, filtered in estimation", {
  # 4 units, unit 4 has no pre-period rows at all
  # g=3, pre: t=1,2
  dt <- data.table::data.table(
    id = c(rep(1, 4), rep(2, 4), rep(3, 4), 4, 4),
    t  = c(1:4, 1:4, 1:4, 3, 4),
    y  = c(1, 2, 3, 4,
           2, 3, 4, 5,
           3, 4, 5, 6,
           10, 11),
    d    = c(rep(1, 4), rep(0, 4), rep(0, 4), 0, 0),
    post = c(rep(c(0, 0, 1, 1), 3), 1, 1)
  )

  cw <- capture_with_warnings(
    transform_demean(dt, "y", "id", "t", g = 3L)
  )
  result <- cw$result

  # Unit 4 should have n_pre = 0 and y_trans = NA
  u4 <- result[result[["id"]] == 4]
  expect_equal(u4[["n_pre"]][1], 0L)
  expect_true(all(is.na(u4[["y_trans"]])))

  # Units 1-3 should have valid y_trans
  u1 <- result[result[["id"]] == 1]
  expect_true(all(!is.na(u1[u1[["t"]] >= 3, ][["y_trans"]])))

  # After summarization and NA filtering, estimation should work with n=3
  # Summarize: mean of post-period y_trans per unit
  post_dt <- result[result[["t"]] >= 3]
  summary_dt <- post_dt[, .(y_trans_summary = mean(y_trans)), by = "id"]
  summary_dt <- merge(summary_dt,
                      unique(dt[, c("id", "d")]),
                      by = "id")

  valid <- !is.na(summary_dt$y_trans_summary)
  expect_equal(sum(valid), 3L)  # unit 4 filtered out

  est <- suppressWarnings(
    estimate_ra_common(
      summary_dt[valid]$y_trans_summary,
      summary_dt[valid]$d
    )
  )
  expect_equal(est$n, 3L)
  expect_true(is.finite(est$att))
})
