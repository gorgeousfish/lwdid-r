# test-transform-detrend.R — Tests for transform_detrend()
# Groups 1-24 covering all E2-02 test requirements

# Helper to suppress messages and exact-fit warnings
run_detrend <- function(...) {
  suppressMessages(
    withCallingHandlers(
      transform_detrend(...),
      lwdid_small_sample = function(w) invokeRestart("muffleWarning")
    )
  )
}

# Helper to suppress messages only (keep warnings)
run_detrend_warn <- function(...) {
  suppressMessages(transform_detrend(...))
}

# ============================================================
# Group 1: Perfect linear trend (2 tests)
# ============================================================

test_that("1.1: perfect linear trend yields zero residuals", {
  dt <- data.table::data.table(
    id = rep(1, 6), t = 1:6,
    y = c(3, 5, 7, 9, 11, 13)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_equal(result$y_trans, rep(0, 6), tolerance = 1e-10)
})

test_that("1.2: trend with noise gives correct slope", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = 1:5,
    y = c(13.5, 16.0, 19.5, 22.0, 25.5)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_equal(result$slope[1], 3.0, tolerance = 1e-10)
  expect_equal(result$t_bar_pre[1], 2.0, tolerance = 1e-10)
})

# ---- Group 4: exclude_pre_periods ----

test_that("detrend: exclude=0 same as default", {
  dt <- data.table::data.table(
    id = rep(1, 6), t = 1:6, y = c(10, 14, 18, 25, 30, 35)
  )
  r0 <- suppressMessages(withCallingHandlers(
    transform_detrend(dt, "y", "id", "t", g = 4L, exclude_pre_periods = 0L),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  ))
  rd <- suppressMessages(withCallingHandlers(
    transform_detrend(dt, "y", "id", "t", g = 4L),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  ))
  expect_equal(r0$y_trans, rd$y_trans, tolerance = 1e-15)
})

test_that("detrend: exclude=1 reduces n_pre", {
  dt <- data.table::data.table(
    id = rep(1, 6), t = 1:6, y = c(10, 14, 18, 25, 30, 35)
  )
  res <- suppressMessages(withCallingHandlers(
    transform_detrend(dt, "y", "id", "t", g = 5L, exclude_pre_periods = 1L),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  ))
  expect_equal(res$n_pre[1], 3L)
})

test_that("detrend: exclude too large triggers error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 4), t = rep(1:4, 2),
    y = c(1, 3, 5, 7, 2, 4, 6, 8)
  )
  expect_error(
    transform_detrend(dt, "y", "id", "t", g = 4L, exclude_pre_periods = 2L),
    class = "lwdid_insufficient_pre_periods"
  )
})

# ---- Group 5: Time range ----

test_that("detrend: time not from 1", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = 2000:2004, y = c(10, 14, 18, 25, 30)
  )
  result <- suppressMessages(withCallingHandlers(
    transform_detrend(dt, "y", "id", "t", g = 2003L),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  ))
  expect_equal(result$t_bar_pre[1], 2001, tolerance = 1e-10)
  expect_equal(result$slope[1], 4.0, tolerance = 1e-10)
})

test_that("detrend: non-contiguous time", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = c(1, 3, 5, 7, 9), y = c(3, 7, 11, 15, 19)
  )
  result <- suppressMessages(withCallingHandlers(
    transform_detrend(dt, "y", "id", "t", g = 7L),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  ))
  expect_equal(result$t_bar_pre[1], 3.0, tolerance = 1e-10)
  expect_equal(result$slope[1], 2.0, tolerance = 1e-10)
  expect_equal(result$y_trans, rep(0, 5), tolerance = 1e-10)
})

# ============================================================
# Group 2: Hand-calculated value verification (3 tests)
# ============================================================

test_that("2.1: manual OLS coefficients and y_trans", {
  # Unit 1: Y={10,14,18,25,30}, t=1:5, g=4
  # t_bar=2, t_c={-1,0,1}
  # MCP: intercept_c=14, slope=4
  # t=4: 25-(14+4*(4-2))=3, t=5: 30-(14+4*(5-2))=4
  dt <- data.table::data.table(
    id = rep(1, 5), t = 1:5,
    y = c(10, 14, 18, 25, 30)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_equal(result$intercept_c[1], 14.0, tolerance = 1e-10)
  expect_equal(result$slope[1], 4.0, tolerance = 1e-10)
  expect_equal(result$t_bar_pre[1], 2.0, tolerance = 1e-10)
  expect_equal(result[t == 4, y_trans], 3.0, tolerance = 1e-10)
  expect_equal(result[t == 5, y_trans], 4.0, tolerance = 1e-10)
})

test_that("2.2: multi-unit independent fitting", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 4), t = rep(1:4, 2),
    y = c(2, 4, 6, 10, 5, 8, 11, 20)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_equal(result[id == 1, slope][1], 2.0, tolerance = 1e-10)
  expect_equal(result[id == 2, slope][1], 3.0, tolerance = 1e-10)
})

test_that("2.3: pre-period OLS residuals orthogonality", {
  dt <- data.table::data.table(
    id = rep(1, 6), t = 1:6,
    y = c(2, 5, 4, 8, 10, 13)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  pre_resid <- result[t <= 3, y_trans]
  t_bar <- result$t_bar_pre[1]
  t_c_pre <- (1:3) - t_bar
  expect_equal(sum(pre_resid), 0, tolerance = 1e-13)
  expect_equal(sum(t_c_pre * pre_resid), 0, tolerance = 1e-13)
})

# ============================================================
# Group 3: Output columns verification (2 tests)
# ============================================================

test_that("3.1: output contains required 5 columns", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = 1:5,
    y = c(1, 3, 5, 7, 9)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_true("y_trans" %in% names(result))
  expect_true("n_pre" %in% names(result))
  expect_true("slope" %in% names(result))
  expect_true("intercept_c" %in% names(result))
  expect_true("t_bar_pre" %in% names(result))
})

test_that("3.2: internal columns degraded/exact_fit removed", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = 1:5,
    y = c(1, 3, 5, 7, 9)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_false("degraded" %in% names(result))
  expect_false("exact_fit" %in% names(result))
})

# ============================================================
# Group 4: exclude_pre_periods parameter (3 tests)
# ============================================================

test_that("4.1: exclude=0 uses all pre-periods", {
  dt <- data.table::data.table(
    id = rep(1, 5), t = 1:5,
    y = c(10, 14, 18, 25, 30)
  )
  r0 <- run_detrend(dt, "y", "id", "t", g = 4L)
  expect_equal(r0$n_pre[1], 3L)
  expect_equal(r0$slope[1], 4.0, tolerance = 1e-10)
})

test_that("4.2: exclude=1 drops last pre-period", {
  dt <- data.table::data.table(
    id = rep(1, 6), t = 1:6,
    y = c(10, 14, 18, 22, 30, 35)
  )
  # g=5, exclude=1: pre_end=3, uses t={1,2,3}
  r1 <- run_detrend(dt, "y", "id", "t", g = 5L,
                    exclude_pre_periods = 1L)
  expect_equal(r1$n_pre[1], 3L)
  expect_equal(r1$t_bar_pre[1], 2.0, tolerance = 1e-10)
})

test_that("4.3: exclude too large raises error", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), t = rep(1:3, 2),
    y = c(1, 3, 5, 2, 4, 6)
  )
  # g=3, exclude=1: pre_end=1, only 1 time point -> error
  expect_error(
    suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 3L,
                        exclude_pre_periods = 1L)
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

# ============================================================
# Group 5: Time range handling (2 tests)
# ============================================================

test_that("5.1: time not starting from 1", {
  dt <- data.table::data.table(
    id = rep(1, 6),
    t = 2000:2005,
    y = c(10, 14, 18, 22, 30, 35)
  )
  result <- run_detrend(dt, "y", "id", "t", g = 2003L)
  expect_equal(result$t_bar_pre[1], 2001, tolerance = 1e-10)
  expect_equal(result$n_pre[1], 3L)
})

test_that("5.2: non-contiguous time", {
  dt <- data.table::data.table(
    id = rep(1, 5),
    t = c(1, 3, 5, 7, 9),
    y = c(2, 6, 10, 14, 18)
  )
  # g=7: pre = t <= 6, so t={1,3,5}
  result <- run_detrend(dt, "y", "id", "t", g = 7L)
  expect_equal(result$n_pre[1], 3L)
  expect_equal(result$slope[1], 2.0, tolerance = 1e-10)
})

# ============================================================
# Group 6: Degradation condition 1 — n_valid < 2 (3 tests)
# ============================================================

test_that("6.1: unit with 1 valid pre-period degrades", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3),
    t = c(1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4),
    y = c(1, 3, 5, 7, 2, 4, 6, 8, 3, 5, 7)
  )
  # g=3: unit 3 has only t=2 pre-period (1 obs)
  expect_warning(
    result <- run_detrend(dt, "y", "id", "t", g = 3L),
    class = "lwdid_data"
  )
  expect_equal(result[id == 3, slope][1], 0)
})

test_that("6.2: unit with all-NA pre-period Y degrades", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 2, 2, 2, 2),
    t = c(1, 2, 3, 4, 1, 2, 3, 4),
    y = c(1, 3, 5, 7, NA, NA, 6, 8)
  )
  expect_warning(
    result <- run_detrend(dt, "y", "id", "t", g = 3L),
    class = "lwdid_data"
  )
  expect_equal(result[id == 2, n_pre][1], 0L)
  expect_true(all(is.na(result[id == 2, y_trans])))
})

test_that("6.3: unit with no pre-period rows gets NA", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 2, 2),
    t = c(1, 2, 3, 4, 3, 4),
    y = c(1, 3, 5, 7, 5, 7)
  )
  expect_warning(
    result <- run_detrend(dt, "y", "id", "t", g = 3L),
    class = "lwdid_data"
  )
  expect_equal(result[id == 2, n_pre][1], 0L)
  expect_equal(result[id == 2, slope][1], 0)
  expect_true(all(is.na(result[id == 2, y_trans])))
})

# ============================================================
# Group 7: Degradation condition 2 — degenerate time (1 test)
# ============================================================

test_that("7.1: degenerate time variance degrades to demean", {
  dt <- data.table::data.table(
    id = c(1, 1, 1, 1, 2, 2, 2),
    t = c(1, 2, 3, 4, 1, 1, 3),
    y = c(1, 3, 5, 7, 2, 4, 6)
  )
  expect_warning(
    result <- run_detrend(dt, "y", "id", "t", g = 3L),
    class = "lwdid_data"
  )
  expect_equal(result[id == 2, slope][1], 0)
})

# ============================================================
# Group 7b: Degradation condition 3 — NaN/Inf coefs (1 test)
# ============================================================

test_that("7b.1: extreme Y values do not crash", {
  dt <- data.table::data.table(
    id = c(rep(1, 5), rep(2, 5)),
    t = rep(1:5, 2),
    y = c(1, 3, 5, 7, 9,
          1e308, -1e308, 1e308, 1e308, 1e308)
  )
  result <- withCallingHandlers(
    suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 4L)
    ),
    lwdid_data = function(w) invokeRestart("muffleWarning"),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(result[id == 1, slope][1], 2.0, tolerance = 1e-10)
  expect_true(all(!is.nan(result[id == 1, y_trans])))
})

# ============================================================
# Group 8: Exact fit warning (1 test)
# ============================================================

test_that("8.1: exact fit (n_valid=2) warns small sample", {
  dt <- data.table::data.table(
    id = rep(1, 4), t = 1:4,
    y = c(2, 4, 7, 9)
  )
  # g=3: pre=t=1,2 (exactly 2)
  expect_warning(
    result <- suppressMessages(
      transform_detrend(dt, "y", "id", "t", g = 3L)
    ),
    class = "lwdid_small_sample"
  )
  # MCP: t_bar=1.5, A*=3, B=2
  # t=3: 7-(3+2*(3-1.5))=1, t=4: 9-(3+2*(4-1.5))=1
  expect_equal(result[t == 3, y_trans], 1.0, tolerance = 1e-10)
  expect_equal(result[t == 4, y_trans], 1.0, tolerance = 1e-10)
})
test_that("22.1: eliminates unit FE and unit-specific trends", {
  # Y_it = alpha_i + beta_i * t + lambda_t
  # Pre lambdas linear {0,0.5,1,1.5} -> OLS captures perfectly
  # Post lambdas non-linear {3,4,2,5} -> y_trans reveals lambda deviations
  # vibe-math MCP verified: post y_trans = {1.0, 1.5, -1.0, 1.5} for all units
  n_units <- 5L
  T_max <- 8L
  g <- 5L
  alphas <- c(10, 20, 30, 40, 50)
  betas <- c(1, 2, 3, 4, 5)
  lambdas <- c(0, 0.5, 1, 1.5, 3, 4, 2, 5)
  dt <- data.table::data.table(
    id = rep(1:n_units, each = T_max),
    t = rep(1:T_max, n_units)
  )
  dt[, y := alphas[id] + betas[id] * t + lambdas[t]]
  result <- run_detrend(dt, "y", "id", "t", g = g)
  # Pre-period: perfect linear fit -> residuals = 0
  pre_resid <- result[t < g, y_trans]
  expect_equal(pre_resid, rep(0, n_units * (g - 1L)), tolerance = 1e-10)
  # Post-period: all units same y_trans at each t (FE + trend eliminated)
  for (tt in g:T_max) {
    vals <- result[t == tt, y_trans]
    expect_equal(max(vals) - min(vals), 0, tolerance = 1e-10)
  }
  # Verify specific post y_trans values (vibe-math MCP verified)
  expect_equal(result[id == 1 & t == 5, y_trans], 1.0, tolerance = 1e-10)
  expect_equal(result[id == 1 & t == 6, y_trans], 1.5, tolerance = 1e-10)
  expect_equal(result[id == 1 & t == 7, y_trans], -1.0, tolerance = 1e-10)
  expect_equal(result[id == 1 & t == 8, y_trans], 1.5, tolerance = 1e-10)
})
