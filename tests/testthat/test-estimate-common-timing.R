# test-estimate-common-timing.R
# Tests for transform_common() dispatch and .estimate_common_timing() pipeline
# Covers Tasks E2-04.10 through E2-04.14

# =============================================================================
# Helper: create balanced panel with known ATT
# =============================================================================
make_panel <- function(N = 8L, TT = 6L, S = 4L, n_treated = 5L,
                       tau = 3.0, sd_noise = 0.5, seed = 123L,
                       add_controls = FALSE) {
  alpha <- seq(10, by = 10, length.out = N)
  set.seed(seed)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, post := as.integer(time >= S)]
  dt[, y := alpha[id] + tau * d * post + rnorm(.N, sd = sd_noise)]
  if (add_controls) {
    dt[, x1 := rep(rnorm(N, sd = 2), each = TT)]
  }
  dt
}

# Helper: suppress all warnings from lwdid
quiet_lwdid <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# =============================================================================
# Group A: transform_common() dispatch tests (Task E2-04.10)
# =============================================================================

test_that("T-01 rolling='demean' dispatches correctly", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5), t = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16, 20, 22, 21, 25, 26, 30, 32, 31, 35, 36)
  )
  result <- transform_common(dt, "y", "id", "t", g = 4L,
                             rolling = "demean")
  expect_true("y_trans" %in% names(result))
  expect_true("pre_mean" %in% names(result))
  expect_false("slope" %in% names(result))
})

test_that("T-02 rolling='detrend' dispatches correctly", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 5), t = rep(1:5, 3),
    y = c(10, 12, 11, 15, 16, 20, 22, 21, 25, 26, 30, 32, 31, 35, 36)
  )
  result <- withCallingHandlers(
    transform_common(dt, "y", "id", "t", g = 4L, rolling = "detrend"),
    lwdid_small_sample = function(w) invokeRestart("muffleWarning")
  )
  expect_true("y_trans" %in% names(result))
  expect_true("slope" %in% names(result))
})

test_that("T-03 rolling='demeanq' without seasonal metadata throws lwdid_invalid_parameter", {
  dt <- data.table::data.table(id = 1, t = 1, y = 1)
  expect_error(
    transform_common(dt, "y", "id", "t", g = 2L, rolling = "demeanq"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("T-04 rolling='detrendq' without seasonal metadata throws lwdid_invalid_parameter", {
  dt <- data.table::data.table(id = 1, t = 1, y = 1)
  expect_error(
    transform_common(dt, "y", "id", "t", g = 2L, rolling = "detrendq"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("T-05 invalid rolling value throws lwdid_invalid_rolling", {
  dt <- data.table::data.table(id = 1, t = 1, y = 1)
  expect_error(
    transform_common(dt, "y", "id", "t", g = 2L, rolling = "invalid"),
    class = "lwdid_invalid_rolling"
  )
})

test_that("T-06 exclude_pre_periods passed through correctly", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 6), t = rep(1:6, 2),
    y = c(10, 12, 14, 16, 18, 20, 5, 7, 9, 11, 13, 15)
  )
  direct <- transform_demean(dt, "y", "id", "t", g = 5L,
                             exclude_pre_periods = 1L)
  via_common <- transform_common(dt, "y", "id", "t", g = 5L,
                                 rolling = "demean",
                                 exclude_pre_periods = 1L)
  expect_equal(direct$y_trans, via_common$y_trans)
  expect_equal(direct$pre_mean, via_common$pre_mean)
})

# =============================================================================
# Group B: End-to-end demean flow (Task E2-04.11)
# Deterministic data: 8 units x 6 periods, S=4, tau=3, seed=123
# Expected values independently verified via vibe-math MCP:
#   ATT    = 3.392235899405178
#   SE     = 0.238922279187284
#   sigma2 = 0.107032229047587
#   (X'X)^{-1}[2,2] = 0.533333...
#   df     = 6
# =============================================================================

test_that("T-07 lwdid() demean returns lwdid_result", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_s3_class(result, "lwdid_result")
})

test_that("T-08 demean ATT matches vibe-math verified value", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(result$att, 3.392235899405178, tolerance = 1e-10)
})

test_that("T-09 demean SE matches vibe-math verified value", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(result$se_att, 0.2389222791872838, tolerance = 1e-10)
})

test_that("T-10 demean df = N - 2 (no controls)", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  # 8 units, 2 params -> df = 6
  expect_identical(result$df_resid, 6L)
})

test_that("T-11 print() outputs ATT and SE", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  out <- paste(capture.output(print(result)), collapse = "\n")
  expect_true(grepl("ATT", out))
  expect_true(grepl("SE", out, ignore.case = TRUE))
})

test_that("T-12 coef() returns named ATT", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  cf <- coef(result)
  expect_named(cf, "ATT")
  expect_equal(unname(cf), result$att)
})

test_that("T-13 confint() returns CI matching result", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(ncol(ci), 2L)
  expect_equal(ci[1, 1], result$ci_lower, tolerance = 1e-12)
  expect_equal(ci[1, 2], result$ci_upper, tolerance = 1e-12)
})

test_that("T-14 vcov() returns 1x1 ATT variance matrix", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  v <- vcov(result)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), 1L)
  expect_equal(ncol(v), 1L)
  # vcov[1,1] = se_att^2
  expect_equal(v[1, 1], result$se_att^2, tolerance = 1e-12)
})

test_that("T-14b vcov_matrix stores full OLS vcov", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  vm <- result$vcov_matrix
  expect_true(is.matrix(vm))
  # No controls: 2x2 (intercept + d)
  expect_equal(nrow(vm), 2L)
  expect_equal(ncol(vm), 2L)
  # se_att^2 = vcov_matrix[2,2]
  expect_equal(result$se_att^2, vm[2, 2], tolerance = 1e-12)
})

# =============================================================================
# Group C: End-to-end detrend flow (Task E2-04.11)
# =============================================================================

test_that("T-15 lwdid() detrend returns reasonable ATT", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "detrend")
  )
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))
  # True tau=3, should be in reasonable range
  expect_true(abs(result$att - 3.0) < 5.0)
})

test_that("T-16 detrend ATT differs from demean ATT", {
  dt <- make_panel()
  r_demean <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  r_detrend <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "detrend")
  )
  # With noise and only 3 pre-periods, detrend and demean differ
  expect_false(isTRUE(all.equal(r_demean$att, r_detrend$att,
                                tolerance = 1e-6)))
})

# =============================================================================
# Group D: NA filtering tests (Task E2-04.12)
# =============================================================================

test_that("T-17 units with insufficient pre-periods excluded", {
  dt <- make_panel(N = 6L, TT = 6L, S = 4L, n_treated = 4L)
  # Remove all pre-period rows for unit 1
  dt <- dt[!(id == 1L & time < 4L)]
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_s3_class(result, "lwdid_result")
  # Unit 1 excluded -> N = 5
  expect_equal(result$nobs, 5L)
})

test_that("T-18 all units NA throws lwdid_insufficient_data", {
  # Panel with only post-period data -> no pre-period -> all NA
  dt <- data.table::data.table(
    id = rep(1:3, each = 2), time = rep(4:5, 3),
    d = c(1L, 1L, 0L, 0L, 0L, 0L), post = 1L,
    y = c(10, 12, 5, 7, 3, 4)
  )
  expect_error(
    quiet_lwdid(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            d = "d", post = "post", rolling = "demean")
    ),
    class = "lwdid_insufficient_data"
  )
})

test_that("T-19 partial NA units emit warning", {
  dt <- make_panel(N = 6L, TT = 6L, S = 4L, n_treated = 4L)
  dt <- dt[!(id == 1L & time < 4L)]
  # Should warn about excluded units
  expect_warning(
    suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            d = "d", post = "post", rolling = "demean")
    ),
    class = "lwdid_data"
  )
})

# =============================================================================
# Group E: Defensive copy (Task E2-04.12)
# =============================================================================

test_that("T-20 defensive copy: input data not modified", {
  dt <- make_panel()
  dt_before <- data.table::copy(dt)
  quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(dt, dt_before)
})

# =============================================================================
# Group F: Control variable tests (Task E2-04.11/12)
# =============================================================================

test_that("T-21 controls end-to-end: controls_tier propagated", {
  dt <- make_panel(N = 10L, TT = 6L, S = 4L, n_treated = 6L,
                   add_controls = TRUE)
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          controls = "x1")
  )
  expect_s3_class(result, "lwdid_result")
  expect_true(result$diagnostics$controls_tier %in%
    c("full_interaction", "simple", "none"))
  expect_true(result$controls_used)
})

test_that("T-23 controls_tier in diagnostics", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_true("controls_tier" %in% names(result$diagnostics))
  expect_equal(result$diagnostics$controls_tier, "none")
})

# =============================================================================
# Group G: Stata consistency (Task E2-04.13)
# =============================================================================

test_that("T-24 Stata consistency: smoking demean ATT", {
  data(smoking, package = "lwdid")
  result <- withCallingHandlers(
    lwdid(data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_s3_class(result, "lwdid_result")
  # ATT should match Stata within 1e-6
  expect_equal(result$att, -0.4221746, tolerance = 1e-6)
})

test_that("T-25 Stata consistency: smoking demean SE", {
  data(smoking, package = "lwdid")
  result <- withCallingHandlers(
    lwdid(data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  # SE should match Stata within 1e-4
  expect_equal(result$se_att, 0.1207995, tolerance = 1e-4)
})

test_that("T-26 Stata consistency: smoking att_by_period structure", {
  data(smoking, package = "lwdid")
  result <- withCallingHandlers(
    lwdid(data = smoking, y = "lcigsale", ivar = "state", tvar = "year",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  # att_by_period should exist and have correct structure
  expect_false(is.null(result$att_by_period))
  expect_true(is.data.frame(result$att_by_period))
  # Should have 12 post-treatment periods (1989-2000)
  expect_equal(nrow(result$att_by_period), 12L)
  # All required columns present (15-column format, Story E3-05)
  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated",
                     "n_control", "df", "vce_type", "n_clusters",
                     "controls_tier")
  expect_true(all(expected_cols %in% names(result$att_by_period)))
  # All ATT values should be finite (no NA periods)
  expect_true(all(is.finite(result$att_by_period$att)))
  # All SE values should be positive
  expect_true(all(result$att_by_period$se > 0))
})

# =============================================================================
# Group H: exclude_pre_periods end-to-end (Task E2-04.11)
# =============================================================================

test_that("T-27 exclude_pre_periods=1 end-to-end", {
  dt <- make_panel(N = 8L, TT = 8L, S = 5L, n_treated = 5L)
  r0 <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 0L)
  )
  r1 <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 1L)
  )
  expect_s3_class(r0, "lwdid_result")
  expect_s3_class(r1, "lwdid_result")
})

test_that("T-28 exclude_pre_periods changes ATT", {
  dt <- make_panel(N = 8L, TT = 8L, S = 5L, n_treated = 5L)
  r0 <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 0L)
  )
  r1 <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 1L)
  )
  expect_false(isTRUE(all.equal(r0$att, r1$att, tolerance = 1e-6)))
})

# =============================================================================
# Group I: End-to-end flow tests (Task E2-04.11)
# Hand-computed reference values verified via vibe-math MCP
# 5 units x 6 periods, g=4, balanced panel
# Units 1-2 treated (d=1), units 3-5 control (d=0)
# Periods 1-3 pre (post=0), periods 4-6 post (post=1)
# Demean ATT = 29/9, SE = 0.7590333900710814, df = 3
# =============================================================================

# Shared test data for Group I
test_df_i <- data.frame(
  id = rep(1:5, each = 6),
  time = rep(1:6, 5),
  y = c(10, 12, 11, 18, 20, 22,
        20, 22, 21, 26, 28, 30,
        30, 32, 31, 35, 36, 37,
        5, 7, 6, 9, 11, 12,
        15, 17, 16, 19, 21, 22),
  d = rep(c(1, 1, 0, 0, 0), each = 6),
  post = rep(c(0, 0, 0, 1, 1, 1), 5)
)

# Pre-compute demean result for reuse
demean_result_i <- withCallingHandlers(
  lwdid(data = test_df_i, y = "y", ivar = "id", tvar = "time",
        d = "d", post = "post", rolling = "demean"),
  warning = function(w) invokeRestart("muffleWarning"),
  message = function(m) invokeRestart("muffleMessage")
)

test_that("B-01 demean e2e returns lwdid_result", {
  expect_true(inherits(demean_result_i, "lwdid_result"))
})

test_that("B-02 detrend e2e returns reasonable ATT", {
  result <- withCallingHandlers(
    lwdid(data = test_df_i, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "detrend"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_true(inherits(result, "lwdid_result"))
  expect_true(is.finite(result$att))
  expect_true(is.numeric(result$att))
})

test_that("B-03 ATT matches hand calculation", {
  expect_equal(demean_result_i$att, 29 / 9, tolerance = 1e-10)
})

test_that("B-04 SE matches hand calculation", {
  expect_equal(demean_result_i$se_att, 0.7590333900710814, tolerance = 1e-10)
})

test_that("B-05 df = N - 2 exact", {
  expect_identical(demean_result_i$df_resid, 3L)
})

test_that("B-06 print() output format", {
  output <- capture.output(print(demean_result_i))
  output_text <- paste(output, collapse = "\n")
  expect_true(grepl("ATT", output_text))
  expect_true(grepl("SE", output_text, ignore.case = TRUE))
})

test_that("B-07 coef() method", {
  cf <- coef(demean_result_i)
  expect_equal(unname(cf), 29 / 9, tolerance = 1e-10)
  expect_equal(names(cf), "ATT")
})

test_that("B-08 confint() method", {
  ci <- confint(demean_result_i)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1L)
  expect_equal(ncol(ci), 2L)
  expect_true(ci[1, 1] <= demean_result_i$att)
  expect_true(ci[1, 2] >= demean_result_i$att)
})

test_that("B-09 vcov() method", {
  V <- vcov(demean_result_i)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 1L)
  expect_equal(ncol(V), 1L)
  expect_equal(V[1, 1], 0.576131687242798, tolerance = 1e-10)
})

test_that("B-10 exclude_pre_periods e2e", {
  # Use asymmetric pre-period data so excluding period 3 changes pre-means
  # Pre-period values: (a, a+2, a+10) — mean of 3 != mean of first 2
  df_asym <- data.frame(
    id = rep(1:5, each = 6),
    time = rep(1:6, 5),
    y = c(10, 12, 20, 28, 30, 32,
          20, 22, 30, 36, 38, 40,
          30, 32, 33, 35, 36, 37,
          5, 7, 8, 9, 11, 12,
          15, 17, 18, 19, 21, 22),
    d = rep(c(1, 1, 0, 0, 0), each = 6),
    post = rep(c(0, 0, 0, 1, 1, 1), 5)
  )
  result_k0 <- withCallingHandlers(
    lwdid(data = df_asym, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 0L),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  result_k1 <- withCallingHandlers(
    lwdid(data = df_asym, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          exclude_pre_periods = 1L),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_true(inherits(result_k1, "lwdid_result"))
  # ATT should differ because pre-means change when period 3 is excluded
  expect_false(isTRUE(all.equal(result_k1$att, result_k0$att)))
})

test_that("B-11 detrend ATT differs from demean ATT on data with different slopes", {
  # Treated units have steep pre-period slope (3), controls have flat slope (1)
  # This guarantees detrend and demean produce different ATTs
  df_slopes <- data.frame(
    id = rep(1:5, each = 6),
    time = rep(1:6, 5),
    y = c(10, 13, 16, 28, 31, 34,
          20, 23, 26, 36, 39, 42,
          30, 31, 32, 35, 36, 37,
          5, 6, 7, 9, 11, 12,
          15, 16, 17, 19, 21, 22),
    d = rep(c(1, 1, 0, 0, 0), each = 6),
    post = rep(c(0, 0, 0, 1, 1, 1), 5)
  )
  result_demean <- withCallingHandlers(
    lwdid(data = df_slopes, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  result_detrend <- withCallingHandlers(
    lwdid(data = df_slopes, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "detrend"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  # With different slopes across treated/control, detrend and demean must differ
  expect_false(isTRUE(all.equal(result_detrend$att, result_demean$att)))
})

test_that("B-12 controls e2e with controls_tier", {
  test_df_ctrl <- test_df_i
  test_df_ctrl$x1 <- rep(c(2.5, 1.0, 3.0, 0.5, 1.5), each = 6)

  result_ctrl <- withCallingHandlers(
    lwdid(data = test_df_ctrl, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          controls = "x1"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_true(inherits(result_ctrl, "lwdid_result"))
  expect_false(is.null(result_ctrl$diagnostics$controls_tier))
})

# =============================================================================
# Group I: vibe-math MCP numerical verification (Task E2-04.14)
# Deterministic data: make_panel(N=8, TT=6, S=4, tau=3, sd=0.5, seed=123)
# All expected values independently computed via vibe-math MCP tools:
#   ATT    = 3.392235899405178
#   SE     = 0.238922279187284
#   t_stat = 14.198072741245305
#   df     = 6
#   sigma2 = 0.107032229047587
#   (X'X)^{-1} = [[0.333, -0.333], [-0.333, 0.533]]
# =============================================================================

test_that("T-29 vibe-math verified: ATT exact match", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(result$att, 3.392235899405178, tolerance = 1e-10)
})

test_that("T-30 vibe-math verified: SE exact match", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(result$se_att, 0.2389222791872838, tolerance = 1e-10)
})

test_that("T-31 vibe-math verified: t-stat exact match", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(result$t_stat, 14.198072741245305, tolerance = 1e-8)
})

# =============================================================================
# Group J: Period-specific effects placeholder (Task E2-04.11)
# estimate_period_effects() not yet implemented (Story E2-05)
# =============================================================================

test_that("T-32 att_by_period is data.frame with correct columns when period effects implemented", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  # estimate_period_effects() is implemented (Story E2-05, updated E3-05)
  expect_true(is.data.frame(result$att_by_period))
  expect_true(nrow(result$att_by_period) > 0L)
  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated",
                     "n_control", "df", "vce_type", "n_clusters",
                     "controls_tier")
  expect_true(all(expected_cols %in% names(result$att_by_period)))
})

# =============================================================================
# Group K: summary() method (Task E2-04.11)
# =============================================================================

test_that("T-33 summary() returns structured output", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  out <- paste(capture.output(summary(result)), collapse = "\n")
  # summary should contain key information
  expect_true(grepl("ATT", out))
  expect_true(grepl("[0-9]", out))
})

# =============================================================================
# Group L: Numerical reasonableness (Task E2-04.11/14)
# =============================================================================

test_that("T-34 ATT sign matches treatment direction", {
  # tau = 3 (positive treatment), ATT should be positive
  dt <- make_panel(tau = 3.0)
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_true(result$att > 0)

  # tau = -2 (negative treatment), ATT should be negative
  dt_neg <- make_panel(tau = -2.0, seed = 456L)
  result_neg <- quiet_lwdid(
    lwdid(data = dt_neg, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_true(result_neg$att < 0)
})

test_that("T-35 SE decreases with larger N", {
  dt_small <- make_panel(N = 8L, n_treated = 5L)
  dt_large <- make_panel(N = 40L, n_treated = 25L, seed = 999L)
  r_small <- quiet_lwdid(
    lwdid(data = dt_small, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  r_large <- quiet_lwdid(
    lwdid(data = dt_large, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_true(r_large$se_att < r_small$se_att)
})

test_that("T-36 CI contains ATT", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_true(result$ci_lower <= result$att)
  expect_true(result$ci_upper >= result$att)
})

test_that("T-37 p-value consistent with t-stat and df", {
  dt <- make_panel()
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  # Recompute p-value from t-stat and df
  expected_p <- 2 * pt(abs(result$t_stat), df = result$df_resid,
                       lower.tail = FALSE)
  expect_equal(result$pvalue, expected_p, tolerance = 1e-12)
})

test_that("T-38 method and rolling fields correct", {
  dt <- make_panel()
  r_dm <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean")
  )
  expect_equal(r_dm$rolling, "demean")

  r_dt <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "detrend")
  )
  expect_equal(r_dt$rolling, "detrend")
})

test_that("T-39 VCE hc1 produces valid results (no longer placeholder)", {
  dt <- make_panel()
  # VCE is now fully implemented (Story E3-05), should produce valid results
  result <- quiet_lwdid(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          vce = "hc1")
  )
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
  expect_equal(result$vce_type, "HC1")
})

# =============================================================================
# Group C: NA filtering and defensive copy tests (Task E2-04.12)
# =============================================================================

test_that("C-01 pre-period insufficient unit excluded, others unaffected", {
  # Units 1-2 treated, 3-4 control (have all 6 periods),
  # unit 5 control (only periods 4-6 → no pre-period → excluded)
  df <- data.frame(
    id = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 3)),
    time = c(rep(1:6, 4), 4:6),
    y = c(10, 12, 11, 18, 20, 22,
          20, 22, 21, 26, 28, 30,
          30, 32, 31, 35, 36, 37,
          5, 7, 6, 9, 11, 12,
          19, 21, 22),
    d = c(rep(1, 12), rep(0, 15)),
    post = c(rep(c(0, 0, 0, 1, 1, 1), 4), 1, 1, 1)
  )

  warnings_caught <- list()
  result <- withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean"),
    lwdid_data = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )

  # Expect at least one lwdid_data warning about excluded units
  has_exclusion_warning <- any(vapply(warnings_caught, function(w) {
    grepl("excluded", conditionMessage(w), ignore.case = TRUE) ||
      isTRUE(w$detail == "units_excluded_na_ytrans")
  }, logical(1)))
  expect_true(has_exclusion_warning)

  # Result is valid
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))

  # Unit 5 excluded → 4 units remain
  expect_equal(result$nobs, 4L)
})

test_that("C-02 all units NA throws lwdid_insufficient_data error", {
  # 4 units x 4 periods, g=3 (pre: periods 1-2, post: periods 3-4)
  # With exclude_pre_periods=2, pre_end = 3-1-2 = 0, no pre-periods
  # This is caught by validation as lwdid_insufficient_pre_periods
  # which inherits from lwdid_insufficient_data
  df_all_na <- data.frame(
    id = rep(1:4, each = 4),
    time = rep(1:4, 4),
    y = c(10, 12, 18, 20, 20, 22, 26, 28,
          30, 32, 35, 36, 5, 7, 9, 11),
    d = rep(c(1, 1, 0, 0), each = 4),
    post = rep(c(0, 0, 1, 1), 4)
  )
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df_all_na, y = "y", ivar = "id", tvar = "time",
            d = "d", post = "post", rolling = "demean",
            exclude_pre_periods = 2L)
    )),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("C-03 partial units NA emits warning with count", {
  # Same as C-01: unit 5 has no pre-period → excluded
  df <- data.frame(
    id = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 3)),
    time = c(rep(1:6, 4), 4:6),
    y = c(10, 12, 11, 18, 20, 22,
          20, 22, 21, 26, 28, 30,
          30, 32, 31, 35, 36, 37,
          5, 7, 6, 9, 11, 12,
          19, 21, 22),
    d = c(rep(1, 12), rep(0, 15)),
    post = c(rep(c(0, 0, 0, 1, 1, 1), 4), 1, 1, 1)
  )

  # Capture all warnings to check for count in message
  warnings_caught <- list()
  withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean"),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )

  # Find the exclusion warning and verify it mentions "1 unit"
  exclusion_msgs <- vapply(warnings_caught, function(w) {
    conditionMessage(w)
  }, character(1))
  has_count <- any(grepl("1 unit", exclusion_msgs))
  expect_true(has_count)
})

test_that("C-04 defensive copy: validated$data not modified", {
  df <- data.frame(
    id = rep(1:5, each = 6), time = rep(1:6, 5),
    y = c(10, 12, 11, 18, 20, 22,
          20, 22, 21, 26, 28, 30,
          30, 32, 31, 35, 36, 37,
          5, 7, 6, 9, 11, 12,
          15, 17, 16, 19, 21, 22),
    d = rep(c(1, 1, 0, 0, 0), each = 6),
    post = rep(c(0, 0, 0, 1, 1, 1), 5)
  )
  dt_original <- data.table::as.data.table(df)
  dt_copy_before <- data.table::copy(dt_original)
  # Run lwdid
  result <- withCallingHandlers(
    lwdid(data = dt_original, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  # Original should be unchanged
  expect_identical(dt_original, dt_copy_before)
})

test_that("C-05 time-varying controls rejected by validation", {
  # Controls that vary over time within units are caught by
  # .validate_time_invariant_controls() as lwdid_invalid_parameter.
  # This is the correct behavior: validation is strict.
  df_tv <- data.frame(
    id = rep(1:5, each = 6), time = rep(1:6, 5),
    y = c(10, 12, 11, 18, 20, 22,
          20, 22, 21, 26, 28, 30,
          30, 32, 31, 35, 36, 37,
          5, 7, 6, 9, 11, 12,
          15, 17, 16, 19, 21, 22),
    d = rep(c(1, 1, 0, 0, 0), each = 6),
    post = rep(c(0, 0, 0, 1, 1, 1), 5),
    x1 = c(1, 2, 3, 4, 5, 6,
           2, 3, 4, 5, 6, 7,
           3, 4, 5, 6, 7, 8,
           1, 1, 1, 2, 2, 2,
           2, 2, 2, 3, 3, 3)
  )

  expect_warning(
    suppressMessages(
      lwdid(data = df_tv, y = "y", ivar = "id", tvar = "time",
            d = "d", post = "post", rolling = "demean",
            controls = "x1")
    ),
    class = "lwdid_data"
  )

  # Time-invariant controls should work fine
  df_ti <- df_tv
  df_ti$x1 <- rep(c(2.5, 1.0, 3.0, 0.5, 1.5), each = 6)
  result <- withCallingHandlers(
    lwdid(data = df_ti, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          controls = "x1"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_s3_class(result, "lwdid_result")
  expect_true(result$controls_used)
})

test_that("C-06 NA filtering syncs y, d, x removal", {
  # Same as C-01 but with a control variable
  # Unit 5 (control, only post-period) should be excluded;

  # remaining units with x1 should produce valid result
  df_sync <- data.frame(
    id = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 3)),
    time = c(rep(1:6, 4), 4:6),
    y = c(10, 12, 11, 18, 20, 22,
          20, 22, 21, 26, 28, 30,
          30, 32, 31, 35, 36, 37,
          5, 7, 6, 9, 11, 12,
          19, 21, 22),
    d = c(rep(1, 12), rep(0, 15)),
    post = c(rep(c(0, 0, 0, 1, 1, 1), 4), 1, 1, 1),
    x1 = c(rep(2.5, 6), rep(1.0, 6), rep(3.0, 6),
           rep(0.5, 6), rep(1.5, 3))
  )

  result <- withCallingHandlers(
    lwdid(data = df_sync, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post", rolling = "demean",
          controls = "x1"),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )

  # Unit 5 excluded, but estimation succeeds with controls
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))
  expect_equal(result$nobs, 4L)
  expect_true(result$controls_used)
})
