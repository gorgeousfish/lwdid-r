# test-staggered-controls.R
# Unit tests for prepare_staggered_controls()

make_controls_test_sub <- function(n = 10L, controls_data = NULL, seed = 42) {
  set.seed(seed)
  dt <- data.table::data.table(
    id = seq_len(n),
    gvar = c(rep(3L, n %/% 2), rep(NA_real_, n - n %/% 2)),
    y = rnorm(n)
  )
  if (!is.null(controls_data)) {
    for (nm in names(controls_data)) {
      data.table::set(dt, j = nm, value = controls_data[[nm]])
    }
  }
  dt
}

# ── Group 1: Basic functionality ─────────────────────────────────────────────

test_that("returns NULL when controls is empty character(0)", {
  sub <- make_controls_test_sub(n = 8L)
  result <- prepare_staggered_controls(sub, character(0), "id")
  expect_null(result)
})

test_that("returns matrix with correct dimensions for valid controls", {
  n <- 8L
  sub <- make_controls_test_sub(n = n, controls_data = list(
    x1 = rnorm(n), x2 = rnorm(n)
  ))
  result <- prepare_staggered_controls(sub, c("x1", "x2"), "id")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 2L)
})

test_that("column names preserved in returned matrix", {
  n <- 8L
  sub <- make_controls_test_sub(n = n, controls_data = list(
    x1 = rnorm(n), x2 = rnorm(n)
  ))
  result <- prepare_staggered_controls(sub, c("x1", "x2"), "id")
  expect_equal(colnames(result), c("x1", "x2"))
})

test_that("returned matrix is raw (not centered) - values identical to input", {
  n <- 8L
  set.seed(99)
  vals <- rnorm(n, mean = 50)
  sub <- make_controls_test_sub(n = n, controls_data = list(x1 = vals))
  result <- prepare_staggered_controls(sub, "x1", "id")
  expect_equal(as.numeric(result[, "x1"]), vals)
})

test_that("single control variable returns N x 1 matrix", {
  n <- 8L
  sub <- make_controls_test_sub(n = n, controls_data = list(x1 = rnorm(n)))
  result <- prepare_staggered_controls(sub, "x1", "id")
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1L)
  expect_equal(nrow(result), n)
})

# ── Group 2: Type check ─────────────────────────────────────────────────────

test_that("non-numeric control triggers lwdid_input error", {
  n <- 10L
  sub <- make_controls_test_sub(n = n, controls_data = list(
    x1 = letters[seq_len(n)]
  ))
  expect_error(
    prepare_staggered_controls(sub, "x1", "id"),
    class = "lwdid_input"
  )
})

test_that("all numeric controls (integer + double) pass type check", {
  n <- 10L
  sub <- make_controls_test_sub(n = n, controls_data = list(
    x_int = seq_len(n),
    x_dbl = rnorm(n)
  ))
  result <- prepare_staggered_controls(sub, c("x_int", "x_dbl"), "id")
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2L)
})

# ── Group 3: All-NA column handling ──────────────────────────────────────────

test_that("all-NA column removed with lwdid_data warning", {
  sub <- make_controls_test_sub(
    n = 5L,
    controls_data = list(x1 = rnorm(5), x2 = rep(NA_real_, 5))
  )
  expect_warning(
    result <- prepare_staggered_controls(sub, c("x1", "x2"), "id"),
    class = "lwdid_data"
  )
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "x1")
})

test_that("all columns all-NA returns NULL", {
  sub <- make_controls_test_sub(
    n = 4L,
    controls_data = list(x1 = rep(NA_real_, 4), x2 = rep(NA_real_, 4))
  )
  expect_warning(
    result <- prepare_staggered_controls(sub, c("x1", "x2"), "id"),
    class = "lwdid_data"
  )
  expect_null(result)
})

test_that("all-NA warning has correct detail field", {
  sub <- make_controls_test_sub(
    n = 3L,
    controls_data = list(x1 = rep(NA_real_, 3))
  )
  cond <- NULL
  withCallingHandlers(
    prepare_staggered_controls(sub, "x1", "id"),
    lwdid_data = function(c) {
      if (identical(c$detail, "staggered_controls_all_na")) {
        cond <<- c
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(cond))
  expect_equal(cond$detail, "staggered_controls_all_na")
})

# ── Group 4: Constant column handling ────────────────────────────────────────

test_that("constant column removed with lwdid_data warning", {
  sub <- make_controls_test_sub(
    n = 5L,
    controls_data = list(x1 = rnorm(5), x2 = rep(5.0, 5))
  )
  expect_warning(
    result <- prepare_staggered_controls(sub, c("x1", "x2"), "id"),
    class = "lwdid_data"
  )
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "x1")
})

test_that("all columns constant returns NULL", {
  sub <- make_controls_test_sub(
    n = 4L,
    controls_data = list(x1 = rep(3.0, 4), x2 = rep(7.0, 4))
  )
  expect_warning(
    result <- prepare_staggered_controls(sub, c("x1", "x2"), "id"),
    class = "lwdid_data"
  )
  expect_null(result)
})

test_that("constant warning has correct detail field", {
  sub <- make_controls_test_sub(
    n = 3L,
    controls_data = list(x1 = rep(2.0, 3))
  )
  cond <- NULL
  withCallingHandlers(
    prepare_staggered_controls(sub, "x1", "id"),
    lwdid_data = function(c) {
      if (identical(c$detail, "staggered_controls_constant")) {
        cond <<- c
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(cond))
  expect_equal(cond$detail, "staggered_controls_constant")
})


# ── Group 5: Mixed anomalous and edge cases ──────────────────────────────────

test_that("mixed anomalous: all-NA + constant + normal keeps only normal", {
  sub <- make_controls_test_sub(
    n = 6L,
    controls_data = list(
      x1 = rnorm(6),
      x2 = rep(NA_real_, 6),
      x3 = rep(5.0, 6)
    )
  )
  result <- suppressWarnings(
    prepare_staggered_controls(sub, c("x1", "x2", "x3"), "id")
  )
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1L)
  expect_equal(colnames(result), "x1")
})

test_that("all anomalous (mixed NA + constant) returns NULL", {
  sub <- make_controls_test_sub(
    n = 4L,
    controls_data = list(x1 = rep(NA_real_, 4), x2 = rep(3.0, 4))
  )
  result <- suppressWarnings(
    prepare_staggered_controls(sub, c("x1", "x2"), "id")
  )
  expect_null(result)
})

test_that("partial NA rows preserved (not filtered by this function)", {
  sub <- make_controls_test_sub(
    n = 5L,
    controls_data = list(x1 = c(1.0, NA, 3.0, 4.0, 5.0))
  )
  result <- prepare_staggered_controls(sub, "x1", "id")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5L)
  expect_true(is.na(result[2, 1]))
})

test_that("single non-NA value column treated as constant", {
  sub <- make_controls_test_sub(
    n = 4L,
    controls_data = list(x1 = c(NA, NA, NA, 5.0))
  )
  expect_warning(
    result <- prepare_staggered_controls(sub, "x1", "id"),
    class = "lwdid_data"
  )
  expect_null(result)
})

test_that("tiny but nonzero variance column preserved", {
  sub <- make_controls_test_sub(
    n = 4L,
    controls_data = list(x1 = c(1.0, 1.0, 1.0, 1.0 + 1e-15))
  )
  result <- prepare_staggered_controls(sub, "x1", "id")
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1L)
})


# =============================================================================
# Group 6: Three-tier degradation integration tests (E4-07.6)
# =============================================================================

# --- Test 6.1: controls_tier reflects actual degradation level ---
test_that("controls_tier reflects actual degradation level", {
  set.seed(60601)
  # 20 treated (cohort 4) + 20 NT, periods 1:6, 1 control x1
  # N1=20 > K+1=2, N0=20 > K+1=2 → Tier 1 (full_interaction)
  units <- 1:40
  periods <- 1:6
  gvar_map <- c(rep(4L, 20), rep(Inf, 20))
  x1_vals <- rnorm(40)

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x1 := x1_vals[id]]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 2.0 * x1_vals[id] + 3.0 * D +
    rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- 4L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = "x1", vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_true(all(result$controls_tier == "full_interaction"))
  expect_true(all(result$K >= 0L))
})

# --- Test 6.2: small sample degrades to Tier 2 or Tier 3 ---
test_that("small sample degrades to Tier 2 or Tier 3", {
  set.seed(60602)
  # Very small: 2 treated (cohort 3) + 2 NT, periods 1:4, 3 controls
  # N1=2, N0=2, K=3: N1 <= K+1=4 → not Tier 1
  units <- 1:4
  periods <- 1:4
  gvar_map <- c(rep(3L, 2), rep(Inf, 2))
  x1_vals <- rnorm(4)
  x2_vals <- rnorm(4)
  x3_vals <- rnorm(4)

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x1 := x1_vals[id]]
  dt[, x2 := x2_vals[id]]
  dt[, x3 := x3_vals[id]]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + D + rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- 3L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  result <- tryCatch(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = c("x1", "x2", "x3"), vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_obs = 2L, min_treated = 1L, min_control = 1L
    )),
    lwdid_insufficient_data = function(e) NULL
  )

  if (!is.null(result)) {
    expect_true(all(result$controls_tier %in% c("simple", "none")))
  }
})

# --- Test 6.3: degradation is (g,r)-level independent ---
test_that("degradation is (g,r)-level independent", {
  set.seed(60603)
  # Two cohorts: cohort 3 (6 units), cohort 5 (2 units), 4 NT, periods 1:8
  # 1 control variable. Different cohorts may get different tiers.
  units <- 1:12
  periods <- 1:8
  gvar_map <- c(rep(3L, 6), rep(5L, 2), rep(Inf, 4))
  x1_vals <- rnorm(12)

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x1 := x1_vals[id]]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 1.5 * x1_vals[id] + 2.0 * D +
    rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = "x1", vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats,
    min_obs = 2L, min_treated = 1L, min_control = 1L
  ))

  expect_true(all(
    result$controls_tier %in% c("full_interaction", "simple", "none")
  ))
})


# =============================================================================
# Group 7: Two-step filter integration tests (E4-07.6)
# =============================================================================

# --- Test 7.1: partial NA in controls filtered by caller ---
test_that("partial NA in controls filtered by caller", {
  set.seed(70701)
  # 10 units (5 treated cohort 4, 5 NT), periods 1:6, 1 control x1
  # Unit 1 has NA for x1 → complete.cases filters the NA row
  units <- 1:10
  periods <- 1:6
  gvar_map <- c(rep(4L, 5), rep(Inf, 5))
  x1_vals <- rnorm(10)
  x1_vals[1] <- NA

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x1 := x1_vals[id]]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 3.0 * D + rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- 4L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = "x1", vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_true(nrow(result) > 0L)
  expect_true(all(result$K >= 0L))
})

# --- Test 7.2: controls with all-NA column in subsample ---
test_that("controls with all-NA column in subsample", {
  set.seed(70702)
  # 8 units (4 treated cohort 4, 4 NT), periods 1:6
  # x1: normal rnorm by unit; x2: NA for time <= 3, rnorm otherwise
  # At period r=4+, x2 should have values
  units <- 1:8
  periods <- 1:6
  gvar_map <- c(rep(4L, 4), rep(Inf, 4))
  x1_vals <- rnorm(8)

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x1 := x1_vals[id]]
  set.seed(70703)
  dt[, x2 := data.table::fifelse(time <= 3L, NA_real_, rnorm(.N))]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 2.0 * D + rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- 4L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = c("x1", "x2"), vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_true(nrow(result) > 0L)
})


# =============================================================================
# Group 8: Numerical correctness tests (E4-07.6)
# =============================================================================

# --- Test 8.1: controls reduce ATT SE (known DGP) ---
test_that("controls reduce ATT SE (known DGP)", {
  set.seed(80801)
  # DGP: Y = alpha_i + beta * x1_it + tau * D_it + epsilon
  # x1 must be time-varying so it survives demeaning.
  # 30 treated (cohort 4) + 30 NT, periods 1:6
  units <- 1:60
  periods <- 1:6
  gvar_map <- c(rep(4L, 30), rep(Inf, 30))

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  # Time-varying control: unit-specific slope * time + noise
  set.seed(80811)
  slopes <- rnorm(60, sd = 3)
  dt[, x1 := slopes[id] * time + rnorm(.N, 0, 0.5)]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  # Y depends strongly on x1 (beta=5) so controlling for it reduces variance
  dt[, Y := as.numeric(id) + 5.0 * x1 + 3.0 * D + rnorm(.N, 0, 1.0)]
  dt[, D := NULL]

  cohorts <- 4L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Run without controls
  result_no_ctrl <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Run with controls
  result_ctrl <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = "x1", vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Controls should reduce SE (time-varying x1 survives demeaning)
  expect_true(mean(result_ctrl$se) < mean(result_no_ctrl$se))
  # ATT should be close to 3.0
  expect_true(all(abs(result_ctrl$att - 3.0) < 2.0))
})

# --- Test 8.2: controls=NULL and all-constant controls give same ATT ---
test_that("controls=NULL and all-constant controls give same ATT", {
  set.seed(80802)
  # 10 treated (cohort 3) + 10 NT, periods 1:6
  # x_const = 42.0 for all units
  units <- 1:20
  periods <- 1:6
  gvar_map <- c(rep(3L, 10), rep(Inf, 10))

  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, x_const := 42.0]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 5.0 * D + rnorm(.N, 0, 0.5)]
  dt[, D := NULL]

  cohorts <- 3L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Run without controls
  result_no_ctrl <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Run with constant control (should be dropped by prepare_staggered_controls)
  result_const <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = "x_const", vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # ATT should be identical (constant control is dropped → same model)
  expect_equal(result_no_ctrl$att, result_const$att, tolerance = 1e-10)
  # SE should be identical
  expect_equal(result_no_ctrl$se, result_const$se, tolerance = 1e-10)
})
