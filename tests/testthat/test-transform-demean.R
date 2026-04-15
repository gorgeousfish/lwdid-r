# ============================================================================
# test-transform-demean.R — Tests for transform_demean()
#
# Test Groups:
#   1: Basic transform correctness
#   2: Hand-calculated value verification
#   3: exclude_pre_periods parameter
#   4: Time range handling
# ============================================================================

# Helper: create simple balanced panel
make_panel <- function(n_units = 3, n_periods = 5, g = 4) {
  dt <- data.table::CJ(id = seq_len(n_units), time = seq_len(n_periods))
  set.seed(42)
  dt[, y := rnorm(.N, mean = id * 10, sd = 1)]
  dt
}

# ==========================================================================
# Group 1 — Basic transform correctness
# ==========================================================================

test_that("1.1: balanced panel y_trans equals y minus pre_mean", {
  dt <- make_panel()
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(result[["y_trans"]], result[["y"]] - result[["pre_mean"]],
               tolerance = 1e-15)
})

test_that("1.2: pre_mean equals arithmetic mean of pre-period Y", {
  dt <- make_panel()
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  for (uid in unique(result[["id"]])) {
    pre_y <- dt[id == uid & time < 4L, y]
    expect_equal(result[id == uid, pre_mean][1], mean(pre_y),
                 tolerance = 1e-15)
  }
})

test_that("1.3: n_pre equals number of pre-period time points", {
  dt <- make_panel()
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_true(all(result[["n_pre"]] == 3L))
})

test_that("1.4: output contains y_trans, pre_mean, n_pre columns", {
  dt <- make_panel()
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_true(all(c("y_trans", "pre_mean", "n_pre") %in% names(result)))
})

test_that("transform_demean preserves row count", {
  dt <- make_panel()
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(nrow(res), nrow(dt))
})

test_that("transform_demean does not modify input dt", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    time = rep(1:3, 2),
    y = 1:6
  )
  dt_copy <- data.table::copy(dt)
  res <- transform_demean(dt, "y", "id", "time", g = 3L)
  expect_identical(dt, dt_copy)
})

# ============================================================================
# Group 2: Hand-calculated value verification (3 tests)
# ============================================================================
test_that("hand-calculated y_trans for known data (MCP-verified)", {
  # Verified: pre_mean = (10+12+11)/3 = 11
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res[id == 1, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-15)
  expect_equal(res[id == 2, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-15)
  expect_equal(res[id == 3, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-15)
})

test_that("each unit pre_mean computed independently", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 4),
    time = rep(1:4, 2),
    y = c(100, 200, 300, 400, 10, 20, 30, 40)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 3L)
  # Unit 1: pre_mean = (100+200)/2 = 150
  # Unit 2: pre_mean = (10+20)/2 = 15
  expect_equal(res[id == 1, pre_mean][1], 150,
               tolerance = 1e-15)
  expect_equal(res[id == 2, pre_mean][1], 15,
               tolerance = 1e-15)
})

test_that("pre-period y_trans sums to zero (MCP-verified)", {
  # Arithmetic sequence: y={3,7,11,15,19,23}, g=4
  # pre_mean=7, y_trans={-4,0,4,8,12,16}, pre_sum=0
  dt <- data.table::data.table(
    id = rep(1, 6), time = 1:6,
    y = c(3, 7, 11, 15, 19, 23)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res$pre_mean[1], 7, tolerance = 1e-15)
  expect_equal(res$y_trans,
               c(-4, 0, 4, 8, 12, 16), tolerance = 1e-15)
  # Mathematical property: pre-period y_trans sums to zero
  pre_sum <- sum(res[time < 4L, y_trans])
  expect_equal(pre_sum, 0, tolerance = 1e-14)
})

# ============================================================================
# Group 3: exclude_pre_periods parameter (4 tests)
# ============================================================================
test_that("exclude_pre_periods=0 uses all pre-treatment periods", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  res0 <- transform_demean(
    dt, "y", "id", "time", g = 4L,
    exclude_pre_periods = 0L
  )
  res_default <- transform_demean(
    dt, "y", "id", "time", g = 4L
  )
  expect_equal(res0$y_trans, res_default$y_trans,
               tolerance = 1e-15)
  expect_equal(res0$pre_mean, res_default$pre_mean,
               tolerance = 1e-15)
  expect_equal(res0$n_pre, res_default$n_pre)
})

test_that("exclude_pre_periods=1 excludes last pre-period", {
  # exclude 1 -> pre_end=2, uses t={1,2}
  # Unit 1: pre_mean = (10+12)/2 = 11, n_pre = 2
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  res <- transform_demean(
    dt, "y", "id", "time", g = 4L,
    exclude_pre_periods = 1L
  )
  expect_equal(res[id == 1, pre_mean][1], 11,
               tolerance = 1e-15)
  expect_equal(res[id == 1, n_pre][1], 2L)
})

test_that("exclude_pre_periods=2 excludes last 2 pre-periods", {
  # exclude 2 -> pre_end=1, uses t={1} only
  # Unit 1: pre_mean = 10, n_pre = 1
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  res <- transform_demean(
    dt, "y", "id", "time", g = 4L,
    exclude_pre_periods = 2L
  )
  expect_equal(res[id == 1, pre_mean][1], 10,
               tolerance = 1e-15)
  expect_equal(res[id == 1, n_pre][1], 1L)
})

test_that("exclude_pre_periods uses calendar cutoff not observation count", {
  # Non-contiguous time: t={2,4,6,8,10}, g=8, exclude=2
  # pre_end = 8-1-2 = 5, so uses t <= 5: {2, 4}
  # y at t=2 is 100, y at t=4 is 200 -> pre_mean = 150
  dt <- data.table::data.table(
    id = rep(1, 5),
    time = c(2, 4, 6, 8, 10),
    y = c(100, 200, 300, 400, 500)
  )
  res <- transform_demean(
    dt, "y", "id", "time", g = 8L,
    exclude_pre_periods = 2L
  )
  expect_equal(res$pre_mean[1], 150, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 2L)
})

# ============================================================================
# Group 4: Time range handling (2 tests)
# ============================================================================
test_that("time not starting from 1 works correctly", {
  # tvar = {5,6,7,8,9}, g=8
  # pre_end = 7, t_min = 5, uses t={5,6,7}
  dt <- data.table::data.table(
    id = rep(1, 5), time = 5:9,
    y = c(50, 60, 70, 80, 90)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 8L)
  # pre_mean = (50+60+70)/3 = 60
  expect_equal(res$pre_mean[1], 60, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 3L)
  expect_equal(res$y_trans,
               c(-10, 0, 10, 20, 30), tolerance = 1e-15)
})

test_that("non-contiguous time works correctly", {
  # tvar = {2,4,6,8,10}, g=8
  # pre_end = 7, t_min = 2, uses t={2,4,6}
  dt <- data.table::data.table(
    id = rep(1, 5),
    time = c(2, 4, 6, 8, 10),
    y = c(20, 40, 60, 80, 100)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 8L)
  # pre_mean = (20+40+60)/3 = 40
  expect_equal(res$pre_mean[1], 40, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 3L)
  expect_equal(res$y_trans,
               c(-20, 0, 20, 40, 60), tolerance = 1e-15)
})


# ============================================================================
# Group 5: Unbalanced panel (4 tests)
# ============================================================================
test_that("partial pre-period missing: correct n_pre and pre_mean", {
  # Unit 1: has t={1,2,3,4,5}, Unit 2: has t={3,4,5} only
  # g=4, pre_end=3
  # Unit 1: pre uses t={1,2,3}, pre_mean = (10+20+30)/3 = 20
  # Unit 2: pre uses t={3}, pre_mean = 300, n_pre = 1
  dt <- data.table::data.table(
    id = c(1,1,1,1,1, 2,2,2),
    time = c(1,2,3,4,5, 3,4,5),
    y = c(10,20,30,40,50, 300,400,500)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res[id == 1, pre_mean][1], 20, tolerance = 1e-15)
  expect_equal(res[id == 1, n_pre][1], 3L)
  expect_equal(res[id == 2, pre_mean][1], 300, tolerance = 1e-15)
  expect_equal(res[id == 2, n_pre][1], 1L)
})

test_that("all-NA pre-period: y_trans is NA, n_pre is 0", {
  # Unit 1: normal, Unit 2: pre-period Y all NA
  dt <- data.table::data.table(
    id = c(1,1,1, 2,2,2),
    time = c(1,2,3, 1,2,3),
    y = c(10,20,30, NA,NA,60)
  )
  res <- withCallingHandlers(
    transform_demean(dt, "y", "id", "time", g = 3L),
    lwdid_data = function(w) invokeRestart("muffleWarning")
  )
  # Unit 1: pre uses t={1,2}, pre_mean = (10+20)/2 = 15, n_pre = 2
  expect_equal(res[id == 1, pre_mean][1], 15, tolerance = 1e-15)
  expect_equal(res[id == 1, n_pre][1], 2L)
  expect_true(is.na(res[id == 2, pre_mean][1]))
  expect_equal(res[id == 2, n_pre][1], 0L)
  expect_true(all(is.na(res[id == 2, y_trans])))
})

test_that("no pre-period rows: y_trans is NA, n_pre is 0", {
  # Unit 1: has t={1,2,3}, Unit 2: has t={3,4} only
  # g=3, pre_end=2, Unit 2 has no rows with t <= 2
  dt <- data.table::data.table(
    id = c(1,1,1, 2,2),
    time = c(1,2,3, 3,4),
    y = c(10,20,30, 50,60)
  )
  res <- withCallingHandlers(
    transform_demean(dt, "y", "id", "time", g = 3L),
    lwdid_data = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(res[id == 1, pre_mean][1], 15, tolerance = 1e-15)
  expect_true(is.na(res[id == 2, pre_mean][1]))
  expect_equal(res[id == 2, n_pre][1], 0L)
  expect_true(all(is.na(res[id == 2, y_trans])))
})

test_that("mixed: normal units unaffected by units with no pre-data", {
  # Unit 1: normal (t=1..5), Unit 2: no pre rows (t=4,5 only)
  # g=4, pre_end=3
  dt <- data.table::data.table(
    id = c(1,1,1,1,1, 2,2),
    time = c(1,2,3,4,5, 4,5),
    y = c(100,200,300,400,500, 40,50)
  )
  res <- withCallingHandlers(
    transform_demean(dt, "y", "id", "time", g = 4L),
    lwdid_data = function(w) invokeRestart("muffleWarning")
  )
  # Unit 1: pre_mean = (100+200+300)/3 = 200 (MCP-verified)
  expect_equal(res[id == 1, pre_mean][1], 200, tolerance = 1e-15)
  expect_equal(res[id == 1, n_pre][1], 3L)
  expect_equal(res[id == 1, y_trans],
               c(-100, 0, 100, 200, 300), tolerance = 1e-15)
  # Unit 2: no pre data
  expect_true(is.na(res[id == 2, pre_mean][1]))
  expect_equal(res[id == 2, n_pre][1], 0L)
})

# ============================================================================
# Group 6: Warnings and errors (5 tests)
# ============================================================================
test_that("units with no pre-data trigger lwdid_data warning", {
  dt <- data.table::data.table(
    id = c(1,1,1, 2,2),
    time = c(1,2,3, 3,4),
    y = c(10,20,30, 50,60)
  )
  expect_warning(
    transform_demean(dt, "y", "id", "time", g = 3L),
    class = "lwdid_data"
  )
})

test_that("warning message contains affected unit count", {
  dt <- data.table::data.table(
    id = c(1,1,1, 2,2, 3,3),
    time = c(1,2,3, 3,4, 3,4),
    y = c(10,20,30, 50,60, 70,80)
  )
  w <- NULL
  withCallingHandlers(
    transform_demean(dt, "y", "id", "time", g = 3L),
    lwdid_data = function(cond) {
      w <<- cond
      invokeRestart("muffleWarning")
    }
  )
  expect_true(grepl("2 unit", w$message))
  expect_equal(w$detail, "units_no_pre_periods")
  expect_equal(w$action_taken, "y_trans set to NA")
})

test_that("exclude_pre_periods too large triggers error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 5),
    time = rep(1:5, 2),
    y = 1:10
  )
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 4L,
                     exclude_pre_periods = 3L),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("error condition has required fields", {
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5, y = 1:5
  )
  err <- tryCatch(
    transform_demean(dt, "y", "id", "time", g = 4L,
                     exclude_pre_periods = 3L),
    lwdid_insufficient_pre_periods = function(e) e
  )
  expect_equal(err$n_pre, 0L)
  expect_equal(err$required, 1L)
  expect_equal(err$rolling, "demean")
})

test_that("g <= t_min triggers error even with exclude=0", {
  # tvar starts at 3, g=3 → pre_end = 2 < t_min = 3
  dt <- data.table::data.table(
    id = rep(1, 3), time = 3:5, y = c(10, 20, 30)
  )
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 3L),
    class = "lwdid_insufficient_pre_periods"
  )
  # Also test g < t_min
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 2L),
    class = "lwdid_insufficient_pre_periods"
  )
})

# ============================================================================
# Group 7: Edge cases (6 tests)
# ============================================================================
test_that("single pre-period (g=2) works", {
  # g=2, only t=1 is pre-period
  dt <- data.table::data.table(
    id = rep(1, 3), time = 1:3,
    y = c(42, 100, 200)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 2L)
  # pre_mean = 42 (MCP-verified)
  expect_equal(res$pre_mean[1], 42, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 1L)
  expect_equal(res$y_trans, c(0, 58, 158), tolerance = 1e-15)
})

test_that("single unit works", {
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5,
    y = c(10, 20, 30, 40, 50)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res$pre_mean[1], 20, tolerance = 1e-15)
  expect_equal(nrow(res), 5L)
})

test_that("single post-period (g = last time) works", {
  # g=5, pre uses t={1,2,3,4}
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5,
    y = c(10, 20, 30, 40, 50)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 5L)
  # pre_mean = (10+20+30+40)/4 = 25 (MCP-verified)
  expect_equal(res$pre_mean[1], 25, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 4L)
})

test_that("exclude_pre_periods exactly eliminates all pre-periods errors", {
  # t={1,2,3,4,5}, g=4, exclude=3 → pre_end = 0 < t_min = 1
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5, y = 1:5
  )
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 4L,
                     exclude_pre_periods = 3L),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("g = t_min triggers error", {
  # tvar = {3,4,5}, g = 3 → pre_end = 2 < t_min = 3
  dt <- data.table::data.table(
    id = rep(1, 3), time = 3:5, y = c(10, 20, 30)
  )
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 3L),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("g < t_min triggers error", {
  # tvar = {3,4,5}, g = 2 → pre_end = 1 < t_min = 3
  dt <- data.table::data.table(
    id = rep(1, 3), time = 3:5, y = c(10, 20, 30)
  )
  expect_error(
    transform_demean(dt, "y", "id", "time", g = 2L),
    class = "lwdid_insufficient_pre_periods"
  )
})

# ============================================================================
# Group 8: NA handling (3 tests)
# ============================================================================
test_that("pre-period partial NA: pre_mean based on non-NA values", {
  # Unit 1: y at t={1,2,3} is {10, NA, 30}, g=4
  # pre_mean = (10+30)/2 = 20 (MCP-verified), n_pre = 2
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5,
    y = c(10, NA, 30, 40, 50)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res$pre_mean[1], 20, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 2L)
  # y_trans at t=1: 10-20 = -10
  expect_equal(res$y_trans[1], -10, tolerance = 1e-15)
  # y_trans at t=2: NA - 20 = NA
  expect_true(is.na(res$y_trans[2]))
})

test_that("post-period NA: y_trans is NA", {
  dt <- data.table::data.table(
    id = rep(1, 4), time = 1:4,
    y = c(10, 20, 30, NA)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 3L)
  # pre_mean = (10+20)/2 = 15
  expect_equal(res$pre_mean[1], 15, tolerance = 1e-15)
  # y_trans at t=4: NA - 15 = NA
  expect_true(is.na(res[time == 4, y_trans]))
  # y_trans at t=3: 30 - 15 = 15
  expect_equal(res[time == 3, y_trans], 15, tolerance = 1e-15)
})

test_that("NA in both pre and post periods handled independently", {
  dt <- data.table::data.table(
    id = rep(1, 5), time = 1:5,
    y = c(10, NA, 30, NA, 50)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  # pre: t={1,2,3}, non-NA: {10, 30}, pre_mean = 20, n_pre = 2
  expect_equal(res$pre_mean[1], 20, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 2L)
  # t=1: 10-20=-10, t=2: NA, t=3: 30-20=10, t=4: NA, t=5: 50-20=30
  expect_equal(res$y_trans[1], -10, tolerance = 1e-15)
  expect_true(is.na(res$y_trans[2]))
  expect_equal(res$y_trans[3], 10, tolerance = 1e-15)
  expect_true(is.na(res$y_trans[4]))
  expect_equal(res$y_trans[5], 30, tolerance = 1e-15)
})


# ============================================================================
# Group 9: Python numerical consistency (3 tests)
# ============================================================================
test_that("9.1: Python-consistent values for standard panel", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  # Unit 1: pre_mean = mean(10,12,11) = 11
  # y_trans = {-1, 1, 0, 4, 5}
  expect_equal(result[id == 1, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-10)
  # Unit 2: pre_mean = mean(20,22,21) = 21
  expect_equal(result[id == 2, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-10)
  # Unit 3: pre_mean = mean(30,32,31) = 31
  expect_equal(result[id == 3, y_trans],
               c(-1, 1, 0, 4, 5), tolerance = 1e-10)
})

test_that("9.2: all units y_trans consistent with expected", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  for (uid in 1:3) {
    pre_y <- dt[id == uid & time <= 3, y]
    pm <- mean(pre_y)
    expected_ytrans <- dt[id == uid, y] - pm
    expect_equal(result[id == uid, y_trans], expected_ytrans,
                 tolerance = 1e-10)
  }
})

test_that("9.3: all-NA pre-period semantically equivalent to Python NaN", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(1, 2, 3, 1, 2, 3),
    y = c(10, 20, 30, NA, NA, 60)
  )
  w <- NULL
  res <- withCallingHandlers(
    transform_demean(dt, "y", "id", "time", g = 3L),
    lwdid_data = function(cond) {
      w <<- cond
      invokeRestart("muffleWarning")
    }
  )
  # Python: DataWarning + ydot = NaN; R: lwdid_data warning + y_trans = NA
  expect_true(!is.null(w))
  expect_true(all(is.na(res[id == 2, y_trans])))
})


# ============================================================================
# Group 10: Mathematical property verification (4 tests)
# ============================================================================
test_that("10.1: pre-period demeaned values sum to zero", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  for (uid in 1:3) {
    pre_sum <- sum(result[id == uid & time < 4L, y_trans])
    expect_equal(pre_sum, 0, tolerance = 1e-14)
  }
})

test_that("10.2: post-period y_trans mean = post mean minus pre mean", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  for (uid in 1:3) {
    post_ytrans_mean <- mean(result[id == uid & time >= 4L, y_trans])
    post_y_mean <- mean(dt[id == uid & time >= 4L, y])
    pre_y_mean <- mean(dt[id == uid & time < 4L, y])
    expect_equal(post_ytrans_mean, post_y_mean - pre_y_mean,
                 tolerance = 1e-14)
  }
})

test_that("10.3: transform preserves post-period inter-unit ranking", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    time = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16,
          20, 22, 21, 25, 26,
          30, 32, 31, 35, 36)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L)
  # Post-period mean of y_trans for each unit
  post_means <- result[time >= 4L, .(m = mean(y_trans)), by = id]
  # Original post-period mean minus pre-period mean
  orig_diffs <- dt[, {
    pre_m <- mean(y[time < 4L])
    post_m <- mean(y[time >= 4L])
    .(d = post_m - pre_m)
  }, by = id]
  # Rankings should be identical
  expect_equal(order(post_means$m), order(orig_diffs$d))
})

test_that("10.4: large panel numerical stability", {
  set.seed(123)
  n_units <- 100L
  n_periods <- 20L
  dt <- data.table::CJ(id = seq_len(n_units),
                        time = seq_len(n_periods))
  dt[, y := rnorm(.N, mean = id * 100, sd = 10)]
  result <- transform_demean(dt, "y", "id", "time", g = 11L)
  # Pre-period y_trans sums to zero for each unit
  for (uid in sample(seq_len(n_units), 10)) {
    pre_sum <- sum(result[id == uid & time < 11L, y_trans])
    expect_equal(pre_sum, 0, tolerance = 1e-10)
  }
})


# ============================================================================
# Group 11: exclude_pre_periods numerical verification (3 tests)
# ============================================================================
test_that("11.1: exclude=1 pre_mean based on t_min..S-2", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 5),
    time = rep(1:5, 2),
    y = c(10, 20, 30, 40, 50,
          100, 200, 300, 400, 500)
  )
  result <- transform_demean(dt, "y", "id", "time", g = 4L,
                             exclude_pre_periods = 1L)
  # pre_end = 2, uses t={1,2}
  # Unit 1: pre_mean = (10+20)/2 = 15
  expect_equal(result[id == 1, pre_mean][1], 15, tolerance = 1e-10)
  expect_equal(result[id == 1, n_pre][1], 2L)
  # Unit 2: pre_mean = (100+200)/2 = 150
  expect_equal(result[id == 2, pre_mean][1], 150, tolerance = 1e-10)
  expect_equal(result[id == 2, n_pre][1], 2L)
})

test_that("11.2: calendar cutoff vs observed-time difference", {
  # Unit with gaps: observed at t={1,2,4}, g=5, k=2
  # R: pre_end = 5-1-2 = 2, uses t={1,2}, pre_mean = mean(10,20) = 15
  dt <- data.table::data.table(
    id = c(1, 1, 1),
    time = c(1L, 2L, 4L),
    y = c(10, 20, 40)
  )
  result <- transform_demean(dt, "y", "id", "time",
                             g = 5L, exclude_pre_periods = 2L)
  expect_equal(result$pre_mean[1], 15, tolerance = 1e-10)
  expect_equal(result$n_pre[1], 2L)
})

test_that("11.3: n_pre decreases monotonically with exclude", {
  dt <- data.table::data.table(
    id = rep(1, 6), time = 1:6,
    y = c(10, 20, 30, 40, 50, 60)
  )
  # g=5, max valid exclude = 3 (pre_end = 1 = t_min)
  n_pre_vals <- integer(4)
  for (k in 0:3) {
    res <- transform_demean(dt, "y", "id", "time", g = 5L,
                            exclude_pre_periods = as.integer(k))
    n_pre_vals[k + 1L] <- res$n_pre[1]
  }
  # n_pre should be 4, 3, 2, 1
  expect_equal(n_pre_vals, c(4L, 3L, 2L, 1L))
  # Monotonically decreasing
  expect_true(all(diff(n_pre_vals) < 0))
})



# ============================================================================
# Group 12: vibe-math MCP verified exact values (2 tests)
# ============================================================================
test_that("12.1: arithmetic panel MCP-verified values", {
  # 5 units x 6 periods, g=4, arithmetic sequences
  # MCP-verified: unit1 pre_mean=7, unit2=9, unit3=4, unit4=6, unit5=12
  dt <- data.table::data.table(
    id = rep(1:5, each = 6),
    time = rep(1:6, 5),
    y = c(3, 7, 11, 15, 19, 23,
          5, 9, 13, 17, 21, 25,
          1, 4, 7, 10, 13, 16,
          2, 6, 10, 14, 18, 22,
          8, 12, 16, 20, 24, 28)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  # MCP calculate("(3+7+11)/3") = 7
  expect_equal(res[id == 1, pre_mean][1], 7, tolerance = 1e-15)
  # MCP calculate("(5+9+13)/3") = 9
  expect_equal(res[id == 2, pre_mean][1], 9, tolerance = 1e-15)
  # MCP calculate("(1+4+7)/3") = 4
  expect_equal(res[id == 3, pre_mean][1], 4, tolerance = 1e-15)
  # MCP calculate("(2+6+10)/3") = 6
  expect_equal(res[id == 4, pre_mean][1], 6, tolerance = 1e-15)
  # MCP calculate("(8+12+16)/3") = 12
  expect_equal(res[id == 5, pre_mean][1], 12, tolerance = 1e-15)
  # Unit 1 y_trans: {3-7, 7-7, 11-7, 15-7, 19-7, 23-7}
  #               = {-4, 0, 4, 8, 12, 16}
  expect_equal(res[id == 1, y_trans],
               c(-4, 0, 4, 8, 12, 16), tolerance = 1e-15)
  # MCP calculate("-4 + 0 + 4") = 0 (pre-period sum)
  pre_sum <- sum(res[id == 1 & time < 4L, y_trans])
  expect_equal(pre_sum, 0, tolerance = 1e-15)
})

test_that("12.2: unbalanced panel MCP-verified values", {
  # Unit with y={10, NA, 30} in pre-period
  # MCP calculate("(10+30)/2") = 20
  dt <- data.table::data.table(
    id = rep(1, 5),
    time = 1:5,
    y = c(10, NA, 30, 40, 50)
  )
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  expect_equal(res$pre_mean[1], 20, tolerance = 1e-15)
  expect_equal(res$n_pre[1], 2L)
  # y_trans: {10-20, NA-20, 30-20, 40-20, 50-20}
  #        = {-10, NA, 10, 20, 30}
  expect_equal(res$y_trans[1], -10, tolerance = 1e-15)
  expect_true(is.na(res$y_trans[2]))
  expect_equal(res$y_trans[3], 10, tolerance = 1e-15)
  expect_equal(res$y_trans[4], 20, tolerance = 1e-15)
  expect_equal(res$y_trans[5], 30, tolerance = 1e-15)
})


# ============================================================================
# Group 13: Fixed-effect elimination property (2 tests)
# ============================================================================
test_that("13.1: alpha_i eliminated in noisy panel", {
  # Y_it(0) = alpha_i + lambda_t + eps_it
  # alpha = {100, 200, 300}, lambda = {1,2,3,4,5}, g=4
  # After demean: y_trans should not depend on alpha_i
  set.seed(999)
  n_units <- 3L
  n_periods <- 5L
  alpha <- c(100, 200, 300)
  lambda <- c(1, 2, 3, 4, 5)
  dt <- data.table::CJ(id = 1:n_units, time = 1:n_periods)
  eps <- rnorm(nrow(dt), sd = 0.001)
  dt[, y := alpha[id] + lambda[time] + eps]
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  # MCP: lambda_pre_mean = (1+2+3)/3 = 2
  # For post-period t=4: y_trans(0) ~ (lambda_4 - lambda_pre) + (eps - eps_pre)
  #                     ~ (4 - 2) + small = 2 + small
  # All units should have nearly identical y_trans at each post-period
  for (tt in 4:5) {
    vals <- res[time == tt, y_trans]
    # All units' y_trans at this time should be very close
    expect_true(max(vals) - min(vals) < 0.01)
  }
})

test_that("13.2: noiseless panel: all units identical y_trans", {
  # Y_it(0) = alpha_i + lambda_t, eps = 0
  # alpha = {100, 200, 300}, lambda = {1,2,3,4,5}, g=4
  # After demean: y_trans = lambda_t - lambda_pre_mean for all units
  alpha <- c(100, 200, 300)
  lambda <- c(1, 2, 3, 4, 5)
  dt <- data.table::CJ(id = 1:3, time = 1:5)
  dt[, y := alpha[id] + lambda[time]]
  res <- transform_demean(dt, "y", "id", "time", g = 4L)
  # MCP: lambda_pre_mean = (1+2+3)/3 = 2
  # y_trans = lambda_t - 2 = {-1, 0, 1, 2, 3} for ALL units
  # MCP: 4 - 2 = 2, 5 - 2 = 3
  expected <- c(-1, 0, 1, 2, 3)
  for (uid in 1:3) {
    expect_equal(res[id == uid, y_trans], expected,
                 tolerance = 1e-14)
  }
})

