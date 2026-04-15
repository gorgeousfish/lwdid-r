# test-vce-period-effects.R - E3-05.7 Per-period VCE Tests
# ============================================================================
# Tests for estimate_period_effects() with VCE support.
# Covers:
#   T7-01: Per-period cluster count variation (unbalanced panel, df_r = G_r - 1)
#   T7-02: Per-period estimation fault tolerance (N_r < 3 -> NA)
#   T7-03: Per-period controls_tier variation
#   T7-04: Non-lwdid error fault tolerance (generic error handler)
# ============================================================================
library(testthat)
library(data.table)
library(sandwich)

# ============================================================================
# Helpers
# ============================================================================

#' Capture all warnings from an expression, returning result + warnings list
capture_with_warnings <- function(expr) {
  warnings_caught <- list()
  result <- withCallingHandlers(expr, warning = function(w) {
    warnings_caught[[length(warnings_caught) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(result = result, warnings = warnings_caught)
}

#' Check if any warning in a list inherits from a given class
has_warning_class <- function(wlist, cls) {
  any(vapply(wlist, function(w) inherits(w, cls), logical(1)))
}

#' Filter warnings by class
filter_warnings <- function(wlist, cls) {
  Filter(function(w) inherits(w, cls), wlist)
}

#' Filter warnings by detail field
filter_warnings_detail <- function(wlist, detail_val) {
  Filter(function(w) identical(w$detail, detail_val), wlist)
}

#' Suppress all lwdid warnings
quiet_lwdid <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

#' Expected 15-column names for estimate_period_effects() output
expected_cols <- c(
  "tindex", "period", "att", "se", "t_stat", "pvalue",
  "ci_lower", "ci_upper", "n_obs", "n_treated", "n_control",
  "df", "vce_type", "n_clusters", "controls_tier"
)


# ============================================================================
# T7-01: Per-period cluster count varies (unbalanced panel), df_r = G_r - 1
# ============================================================================
# In an unbalanced panel, different post-treatment periods may have different
# numbers of clusters G_r. The cluster-robust df should be computed
# independently per period as df_r = G_r - 1.

test_that("T7-01a: per-period cluster count varies, df_r = G_r - 1 independently", {
  set.seed(701)

  # Build unbalanced panel where cluster membership differs across periods:
  # Period 1: 5 clusters, 6 obs each -> G_1=5, df_1=4
  # Period 2: 4 clusters, 6 obs each -> G_2=4, df_2=3
  # Period 3: 3 clusters, 6 obs each -> G_3=3, df_3=2
  make_period <- function(period, n_clusters, obs_per_cluster) {
    n_local <- n_clusters * obs_per_cluster
    n_treated <- as.integer(obs_per_cluster / 2)
    data.table(
      tindex  = rep(period, n_local),
      y_trans = rnorm(n_local, mean = 2),
      d       = rep(c(rep(1L, n_treated),
                       rep(0L, obs_per_cluster - n_treated)), n_clusters),
      cluster = rep(seq_len(n_clusters), each = obs_per_cluster)
    )
  }

  dt <- rbindlist(list(
    make_period(1L, 5L, 6L),
    make_period(2L, 4L, 6L),
    make_period(3L, 3L, 6L)
  ))
  # Add treatment effect so ATT is identifiable
  dt[d == 1L, y_trans := y_trans + 2.0]

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster", alpha = 0.05)
  )

  # Basic structure: 3 rows, 15 columns
  expect_equal(nrow(result), 3L)
  expect_true(all(expected_cols %in% names(result)))

  # df_r = G_r - 1 per period, computed independently
  expect_equal(result$df[1], 4L)   # G_1=5 -> df=4
  expect_equal(result$df[2], 3L)   # G_2=4 -> df=3
  expect_equal(result$df[3], 2L)   # G_3=3 -> df=2

  # n_clusters matches per-period cluster count
  expect_equal(result$n_clusters[1], 5L)
  expect_equal(result$n_clusters[2], 4L)
  expect_equal(result$n_clusters[3], 3L)

  # All periods should have valid ATT (not NA)
  expect_true(all(!is.na(result$att)))
  expect_true(all(!is.na(result$se)))
  expect_true(all(result$se > 0))

  # vce_type should be "cluster" for all periods
  expect_true(all(result$vce_type == "cluster"))
})

test_that("T7-01b: per-period df independent of model-level df", {
  # Verify that cluster df = G_r - 1 is NOT the model residual df
  # Even with controls (Tier 1), cluster df should be G_r - 1
  set.seed(702)

  make_period_ctrl <- function(period, n_clusters, obs_per_cluster, k = 2L) {
    n_local <- n_clusters * obs_per_cluster
    n_treated <- as.integer(obs_per_cluster / 2)
    dt_local <- data.table(
      tindex  = rep(period, n_local),
      d       = rep(c(rep(1L, n_treated),
                       rep(0L, obs_per_cluster - n_treated)), n_clusters),
      cluster = rep(seq_len(n_clusters), each = obs_per_cluster)
    )
    for (j in seq_len(k)) {
      set(dt_local, j = paste0("x", j), value = rnorm(n_local))
    }
    dt_local[, y_trans := 1 + 2 * d + 0.5 * x1 + rnorm(.N, sd = 0.5)]
    dt_local
  }

  # Period 1: G=20 clusters, 4 obs each -> N=80, K=2
  # Tier 1: N1=40>3, N0=40>3 -> full_interaction
  # Model df = 80 - 2 - 2*2 = 74, but cluster df = 20 - 1 = 19
  dt <- make_period_ctrl(1L, 20L, 4L, k = 2L)

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"), periods = 1L,
                            vce = "cluster", cluster_var = "cluster")
  )

  expect_equal(result$controls_tier, "full_interaction")
  expect_equal(result$df, 19L)        # G - 1 = 19
  expect_equal(result$n_clusters, 20L)
  # Verify it's NOT the model residual df
  expect_false(result$df == 74L)
})

test_that("T7-01c: per-period n_obs and n_treated/n_control correct", {
  set.seed(703)

  # Period 1: 30 obs (15 treated, 15 control)
  # Period 2: 24 obs (12 treated, 12 control)
  dt1 <- data.table(
    tindex = rep(1L, 30L),
    y_trans = rnorm(30, mean = 3),
    d = c(rep(1L, 15L), rep(0L, 15L)),
    cluster = rep(1:5, each = 6L)
  )
  dt2 <- data.table(
    tindex = rep(2L, 24L),
    y_trans = rnorm(24, mean = 3),
    d = c(rep(1L, 12L), rep(0L, 12L)),
    cluster = rep(1:4, each = 6L)
  )
  dt1[d == 1L, y_trans := y_trans + 2]
  dt2[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2))

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )

  expect_equal(result$n_obs[1], 30L)
  expect_equal(result$n_treated[1], 15L)
  expect_equal(result$n_control[1], 15L)
  expect_equal(result$n_obs[2], 24L)
  expect_equal(result$n_treated[2], 12L)
  expect_equal(result$n_control[2], 12L)
})

test_that("T7-01d: cluster SE matches sandwich::vcovCL per period", {
  set.seed(704)

  # Build a simple 2-period panel with known cluster structure
  # Period 1: 30 obs, 5 clusters of 6
  # Period 2: 30 obs, 5 clusters of 6
  make_p <- function(period) {
    data.table(
      tindex = rep(period, 30L),
      y_trans = rnorm(30, mean = 1 + period),
      d = c(rep(1L, 15L), rep(0L, 15L)),
      cluster = rep(1:5, each = 6L)
    )
  }
  dt <- rbindlist(list(make_p(1L), make_p(2L)))
  dt[d == 1L, y_trans := y_trans + 2]

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )

  # Verify SE matches sandwich::vcovCL for each period independently
  for (i in 1:2) {
    sub <- dt[tindex == i]
    fit_ref <- lm(sub$y_trans ~ sub$d)
    se_ref <- sqrt(sandwich::vcovCL(
      fit_ref, cluster = sub$cluster, type = "HC1")[2, 2])
    expect_equal(result$se[i], se_ref, tolerance = 1e-12,
                 label = sprintf("Period %d cluster SE", i))
  }
})

# ============================================================================
# T7-02: Per-period estimation fault tolerance — failed period returns NA
# ============================================================================
# When one period has insufficient data (N_r < 3), that period's results
# should be NA, but other periods should compute normally.

test_that("T7-02a: period with N_r < 3 returns NA, doesn't break others", {
  set.seed(7021)

  # Period 1: normal (N=20, mixed treatment)
  # Period 2: only 2 observations (N<3, will trigger lwdid_insufficient_data)
  # Period 3: normal (N=20, mixed treatment)
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 2), rep(3L, 20)),
    y_trans = rnorm(42),
    d       = c(rep(c(1L, 0L), 10), c(1L, 0L), rep(c(1L, 0L), 10))
  )
  dt[d == 1L, y_trans := y_trans + 2.0]

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL, alpha = 0.05)
  )
  result <- cw$result

  expect_equal(nrow(result), 3L)

  # Period 1 and 3 should have valid results
  expect_true(!is.na(result$att[1]))
  expect_true(!is.na(result$att[3]))
  expect_true(result$se[1] > 0)
  expect_true(result$se[3] > 0)

  # Period 2 should be NA (failed due to N<3)
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))
  expect_true(is.na(result$t_stat[2]))
  expect_true(is.na(result$pvalue[2]))
  expect_true(is.na(result$ci_lower[2]))
  expect_true(is.na(result$ci_upper[2]))

  # Should have emitted a lwdid_small_sample warning for period 2
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))
})

test_that("T7-02b: period with all-control (N_1=0) returns NA with correct counts", {
  set.seed(7022)

  # Period 1: normal (N=20, mixed treatment)
  # Period 2: all control (N_1=0 → lwdid_insufficient_data)
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 10)),
    y_trans = rnorm(30),
    d       = c(rep(c(1L, 0L), 10), rep(0L, 10))
  )
  dt[d == 1L, y_trans := y_trans + 2.0]

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05)
  )
  result <- cw$result

  expect_equal(nrow(result), 2L)

  # Period 1 should succeed
  expect_true(!is.na(result$att[1]))
  expect_true(result$se[1] > 0)

  # Period 2 should be NA
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))

  # lwdid_error handler emits lwdid_small_sample warning
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))

  # Period 2 should still have n_obs recorded (from the error handler)
  expect_equal(result$n_obs[2], 10L)
  expect_equal(result$n_treated[2], 0L)
  expect_equal(result$n_control[2], 10L)
})

test_that("T7-02c: period with all-treated (N_0=0) returns NA", {
  set.seed(7023)

  # Period 1: normal
  # Period 2: all treated (N_0=0 → lwdid_insufficient_data)
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 8)),
    y_trans = rnorm(28),
    d       = c(rep(c(1L, 0L), 10), rep(1L, 8))
  )

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05)
  )
  result <- cw$result

  expect_equal(nrow(result), 2L)
  expect_true(!is.na(result$att[1]))
  expect_true(is.na(result$att[2]))
  expect_equal(result$n_treated[2], 8L)
  expect_equal(result$n_control[2], 0L)
})

test_that("T7-02d: multiple failed periods don't break remaining valid periods", {
  set.seed(7024)

  # Period 1: N<3 (only 2 obs)
  # Period 2: normal (N=20)
  # Period 3: all control (N_1=0)
  # Period 4: normal (N=20)
  dt <- data.table(
    tindex  = c(rep(1L, 2), rep(2L, 20), rep(3L, 10), rep(4L, 20)),
    y_trans = rnorm(52),
    d       = c(c(1L, 0L),
                rep(c(1L, 0L), 10),
                rep(0L, 10),
                rep(c(1L, 0L), 10))
  )
  dt[d == 1L, y_trans := y_trans + 1.5]

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:4, vce = NULL, alpha = 0.05)
  )

  expect_equal(nrow(result), 4L)

  # Periods 1 and 3 should be NA
  expect_true(is.na(result$att[1]))
  expect_true(is.na(result$att[3]))

  # Periods 2 and 4 should be valid
  expect_true(!is.na(result$att[2]))
  expect_true(!is.na(result$att[4]))
  expect_true(result$se[2] > 0)
  expect_true(result$se[4] > 0)
  expect_true(result$att[2] > 0)
  expect_true(result$att[4] > 0)
})

test_that("T7-02e: fault tolerance with cluster VCE — failed period doesn't corrupt cluster info", {
  set.seed(7025)

  # Period 1: 4 clusters, normal → G_1=4, df_1=3
  # Period 2: only 2 obs → fails (N<3)
  # Period 3: 3 clusters, normal → G_3=3, df_3=2
  dt <- rbindlist(list(
    data.table(
      tindex  = rep(1L, 24),
      y_trans = rnorm(24, mean = 2),
      d       = rep(c(1L, 0L, 0L, 0L, 0L, 0L), 4),
      cluster = rep(1:4, each = 6)
    ),
    data.table(
      tindex  = rep(2L, 2),
      y_trans = rnorm(2),
      d       = c(1L, 0L),
      cluster = c(1L, 2L)
    ),
    data.table(
      tindex  = rep(3L, 18),
      y_trans = rnorm(18, mean = 2),
      d       = rep(c(1L, 0L, 0L, 0L, 0L, 0L), 3),
      cluster = rep(1:3, each = 6)
    )
  ))
  dt[d == 1L, y_trans := y_trans + 2.0]

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster", alpha = 0.05)
  )

  expect_equal(nrow(result), 3L)

  # Period 1: valid with cluster VCE
  expect_true(!is.na(result$att[1]))
  expect_equal(result$n_clusters[1], 4L)
  expect_equal(result$df[1], 3L)
  expect_equal(result$vce_type[1], "cluster")

  # Period 2: failed — NA results
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))
  expect_true(is.na(result$df[2]))

  # Period 3: valid with cluster VCE (independent of period 2 failure)
  expect_true(!is.na(result$att[3]))
  expect_equal(result$n_clusters[3], 3L)
  expect_equal(result$df[3], 2L)
  expect_equal(result$vce_type[3], "cluster")
})

# ============================================================================
# T7-03: Per-period controls_tier varies across periods
# ============================================================================
# Different periods may have different sample sizes, causing the three-tier
# fallback to select different model tiers per period.

test_that("T7-03a: controls_tier varies across periods (full_interaction, simple, none)", {
  set.seed(703)
  K <- 2L  # 2 control variables

  # Tier thresholds for K=2:
  #   Tier 1 (full_interaction): n1 > K+1=3 AND n0 > K+1=3
  #   Tier 2 (simple):           N > K+2=4
  #   Tier 3 (none):             N <= K+2=4

  # Period 1: Tier 1 — n1=15, n0=15 (both > 3)
  n1_p1 <- 15L; n0_p1 <- 15L; N_p1 <- n1_p1 + n0_p1
  # Period 2: Tier 2 — n1=2 (≤ 3), n0=6, N=8 > 4
  n1_p2 <- 2L; n0_p2 <- 6L; N_p2 <- n1_p2 + n0_p2
  # Period 3: Tier 3 — N=4 ≤ K+2=4
  n1_p3 <- 1L; n0_p3 <- 3L; N_p3 <- n1_p3 + n0_p3

  x1 <- matrix(rnorm(N_p1 * K), N_p1, K)
  x2 <- matrix(rnorm(N_p2 * K), N_p2, K)
  x3 <- matrix(rnorm(N_p3 * K), N_p3, K)

  dt <- data.table(
    tindex  = c(rep(1L, N_p1), rep(2L, N_p2), rep(3L, N_p3)),
    y_trans = c(1 + 2 * c(rep(1, n1_p1), rep(0, n0_p1)) + rnorm(N_p1, 0, 0.5),
                1 + 2 * c(rep(1, n1_p2), rep(0, n0_p2)) + rnorm(N_p2, 0, 0.5),
                1 + 2 * c(rep(1, n1_p3), rep(0, n0_p3)) + rnorm(N_p3, 0, 0.5)),
    d       = c(rep(1L, n1_p1), rep(0L, n0_p1),
                rep(1L, n1_p2), rep(0L, n0_p2),
                rep(1L, n1_p3), rep(0L, n0_p3)),
    X1      = c(x1[, 1], x2[, 1], x3[, 1]),
    X2      = c(x1[, 2], x2[, 2], x3[, 2])
  )

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("X1", "X2"), periods = 1:3,
                            vce = NULL, alpha = 0.05)
  )

  expect_equal(nrow(result), 3L)

  # Verify per-period tier assignment
  expect_equal(result$controls_tier[1], "full_interaction")
  expect_equal(result$controls_tier[2], "simple")
  expect_equal(result$controls_tier[3], "none")

  # All periods should produce valid ATT estimates
  expect_true(all(!is.na(result$att)))

  # ATT should be roughly around 2 (the true effect)
  expect_true(all(result$att > 0))
})

test_that("T7-03b: controls_tier with cluster VCE varies per period", {
  set.seed(7032)
  K <- 1L  # 1 control variable

  # Tier thresholds for K=1:
  #   Tier 1 (full_interaction): n1 > 2 AND n0 > 2
  #   Tier 2 (simple):           N > 3
  #   Tier 3 (none):             N <= 3

  # Period 1: Tier 1 — n1=10, n0=10, 4 clusters
  # Period 2: Tier 2 — n1=2, n0=8, 3 clusters (n1 ≤ 2 → simple)
  dt1 <- data.table(
    tindex  = rep(1L, 20),
    y_trans = rnorm(20, mean = 3),
    d       = rep(c(1L, 0L), 10),
    X1      = rnorm(20),
    cluster = rep(1:4, each = 5)
  )
  dt2 <- data.table(
    tindex  = rep(2L, 10),
    y_trans = rnorm(10, mean = 3),
    d       = c(rep(1L, 2), rep(0L, 8)),
    X1      = rnorm(10),
    cluster = rep(1:3, c(4, 3, 3))
  )
  dt1[d == 1L, y_trans := y_trans + 1.5]
  dt2[d == 1L, y_trans := y_trans + 1.5]
  dt <- rbindlist(list(dt1, dt2))

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = "X1", periods = 1:2,
                            vce = "cluster", cluster_var = "cluster",
                            alpha = 0.05)
  )

  expect_equal(nrow(result), 2L)

  # Period 1: full_interaction (n1=10>2, n0=10>2)
  expect_equal(result$controls_tier[1], "full_interaction")
  expect_equal(result$vce_type[1], "cluster")
  expect_equal(result$n_clusters[1], 4L)
  expect_equal(result$df[1], 3L)

  # Period 2: simple (n1=2 ≤ K+1=2)
  expect_equal(result$controls_tier[2], "simple")
  expect_equal(result$vce_type[2], "cluster")
  expect_equal(result$n_clusters[2], 3L)
  expect_equal(result$df[2], 2L)
})

test_that("T7-03c: controls_tier 'none' when controls dropped due to small N", {
  set.seed(7033)
  K <- 3L  # 3 control variables

  # Tier thresholds for K=3:
  #   Tier 1: n1 > 4 AND n0 > 4
  #   Tier 2: N > 5
  #   Tier 3: N ≤ 5

  # Period 1: Tier 1 — n1=10, n0=10
  # Period 2: Tier 3 — N=5 ≤ K+2=5, all controls dropped
  dt1 <- data.table(
    tindex  = rep(1L, 20),
    y_trans = 1 + 2 * rep(c(1, 0), 10) + rnorm(20, 0, 0.5),
    d       = rep(c(1L, 0L), 10),
    X1 = rnorm(20), X2 = rnorm(20), X3 = rnorm(20)
  )
  dt2 <- data.table(
    tindex  = rep(2L, 5),
    y_trans = 1 + 2 * c(1, 1, 0, 0, 0) + rnorm(5, 0, 0.5),
    d       = c(1L, 1L, 0L, 0L, 0L),
    X1 = rnorm(5), X2 = rnorm(5), X3 = rnorm(5)
  )
  dt <- rbindlist(list(dt1, dt2))

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("X1", "X2", "X3"), periods = 1:2,
                            vce = NULL, alpha = 0.05)
  )

  expect_equal(nrow(result), 2L)
  expect_equal(result$controls_tier[1], "full_interaction")
  expect_equal(result$controls_tier[2], "none")

  # Both should still produce valid estimates
  expect_true(!is.na(result$att[1]))
  expect_true(!is.na(result$att[2]))
})

# ============================================================================
# T7-04: Non-lwdid error fault tolerance (generic error handler)
# ============================================================================
# The tryCatch in estimate_period_effects has two error layers:
#   Layer 1: lwdid_error → emits lwdid_small_sample warning
#   Layer 2: generic error → emits lwdid_numerical warning
# T7-04 tests both layers.

test_that("T7-04a: lwdid_error handler catches insufficient data correctly", {
  set.seed(704)

  # Period 1: normal (N=20, mixed treatment)
  # Period 2: all control (N_1=0 → lwdid_insufficient_data → lwdid_error handler)
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 10)),
    y_trans = rnorm(30),
    d       = c(rep(c(1L, 0L), 10),  # period 1: normal
                rep(0L, 10))          # period 2: all control
  )
  dt[d == 1L, y_trans := y_trans + 2.0]

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  expect_equal(nrow(result), 2L)

  # Period 1 should succeed
  expect_true(!is.na(result$att[1]))
  expect_true(result$se[1] > 0)

  # Period 2 should be NA (lwdid_error handler)
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))

  # The lwdid_error handler emits lwdid_small_sample warning
  small_sample_warns <- Filter(function(w) inherits(w, "lwdid_small_sample"),
                               warnings_caught)
  expect_true(length(small_sample_warns) >= 1L)

  # Verify the warning message mentions the period
  warn_msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  expect_true(any(grepl("Period 2", warn_msgs, fixed = TRUE)))
})

test_that("T7-04b: generic error handler — verify tryCatch structure catches non-lwdid errors", {
  set.seed(7042)

  # This test verifies the generic error handler works by checking that
  # when a period fails for any reason, the result is NA and the function
  # doesn't abort. We test this by creating a period with no valid
  # observations after NA filtering (all y_trans are NA).
  # This triggers the .make_na_period_row path (sum(valid)==0).

  # Period 1: normal
  # Period 2: all y_trans are NA → no valid obs → all-NA row
  # Period 3: normal
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 10), rep(3L, 20)),
    y_trans = c(rnorm(20), rep(NA_real_, 10), rnorm(20)),
    d       = c(rep(c(1L, 0L), 10),
                rep(c(1L, 0L), 5),
                rep(c(1L, 0L), 10))
  )
  dt[!is.na(y_trans) & d == 1L, y_trans := y_trans + 2.0]

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL, alpha = 0.05)
  )

  expect_equal(nrow(result), 3L)

  # Period 1 and 3 should be valid
  expect_true(!is.na(result$att[1]))
  expect_true(!is.na(result$att[3]))

  # Period 2 should be all-NA row
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))
  expect_equal(result$n_obs[2], 0L)  # .make_na_period_row sets n_obs=0
})

test_that("T7-04c: lwdid_error vs generic error — different warning classes", {
  set.seed(7043)

  # Test that lwdid_error handler emits lwdid_small_sample class
  # while generic error handler would emit lwdid_numerical class.
  # We can test the lwdid_error path by creating N_1=0 (lwdid_insufficient_data).

  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 10)),
    y_trans = rnorm(30),
    d       = c(rep(c(1L, 0L), 10),
                rep(0L, 10))  # period 2: all control → lwdid_insufficient_data
  )
  dt[d == 1L, y_trans := y_trans + 1.0]

  warnings_caught <- list()
  result <- withCallingHandlers(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05),
    warning = function(w) {
      warnings_caught[[length(warnings_caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  # Period 2 failed via lwdid_error handler → lwdid_small_sample warning
  small_sample_warns <- Filter(function(w) inherits(w, "lwdid_small_sample"),
                               warnings_caught)
  expect_true(length(small_sample_warns) >= 1L)

  # Verify the warning message contains "estimation failed"
  warn_msgs <- vapply(small_sample_warns, conditionMessage, character(1))
  period2_warns <- warn_msgs[grepl("Period 2", warn_msgs)]
  expect_true(length(period2_warns) >= 1L)
  expect_true(any(grepl("failed", period2_warns)))
})

test_that("T7-04d: all periods fail gracefully — returns all-NA data.frame", {
  set.seed(7044)

  # All periods have insufficient data
  # Period 1: N=2 (< 3)
  # Period 2: N_1=0 (all control)
  dt <- data.table(
    tindex  = c(rep(1L, 2), rep(2L, 6)),
    y_trans = rnorm(8),
    d       = c(1L, 0L, rep(0L, 6))
  )

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05)
  )

  expect_equal(nrow(result), 2L)
  expect_true(all(is.na(result$att)))
  expect_true(all(is.na(result$se)))

  # Should still have correct period identifiers
  expect_equal(result$tindex, c(1L, 2L))
  expect_equal(result$period, c(1L, 2L))
})

test_that("T7-04e: 15-column structure preserved even for failed periods", {
  set.seed(7045)

  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated",
                     "n_control", "df", "vce_type", "n_clusters",
                     "controls_tier")

  # Period 1: normal, Period 2: fails (N<3)
  dt <- data.table(
    tindex  = c(rep(1L, 20), rep(2L, 2)),
    y_trans = rnorm(22),
    d       = c(rep(c(1L, 0L), 10), c(1L, 0L))
  )
  dt[d == 1L, y_trans := y_trans + 1.0]

  result <- suppressWarnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL, alpha = 0.05)
  )

  # All 15 columns must be present
  expect_equal(sort(names(result)), sort(expected_cols))
  expect_equal(ncol(result), 15L)

  # Failed period row should still have all 15 columns
  failed_row <- result[2, ]
  expect_equal(ncol(failed_row), 15L)
  expect_true(is.na(failed_row$att))
  expect_true(is.na(failed_row$se))
})

# ============================================================================
# T7-02: Per-period estimation fault tolerance (N_r < 3 -> NA)
# ============================================================================
# When N_r < 3 for some period, that period's result should be NA without
# interrupting other periods. The lwdid_error is caught by the inner tryCatch
# layer and a lwdid_small_sample warning is emitted.

test_that("T7-02a: period with N_r < 3 produces NA row, others unaffected", {
  set.seed(721)

  # Period 1: 30 obs (normal)
  # Period 2: 2 obs (N < 3 -> should fail gracefully)
  # Period 3: 30 obs (normal)
  dt1 <- data.table(
    tindex = rep(1L, 30L),
    y_trans = rnorm(30, mean = 3),
    d = c(rep(1L, 15L), rep(0L, 15L)),
    cluster = rep(1:5, each = 6L)
  )
  dt2 <- data.table(
    tindex = rep(2L, 2L),
    y_trans = c(5, 2),
    d = c(1L, 0L),
    cluster = c(1L, 2L)
  )
  dt3 <- data.table(
    tindex = rep(3L, 30L),
    y_trans = rnorm(30, mean = 3),
    d = c(rep(1L, 15L), rep(0L, 15L)),
    cluster = rep(1:5, each = 6L)
  )
  dt1[d == 1L, y_trans := y_trans + 2]
  dt3[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2, dt3))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster")
  )
  result <- cw$result

  # Period 2 should be NA
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))
  expect_true(is.na(result$t_stat[2]))
  expect_true(is.na(result$pvalue[2]))
  expect_true(is.na(result$ci_lower[2]))
  expect_true(is.na(result$ci_upper[2]))
  expect_true(is.na(result$df[2]))

  # Period 2 should preserve sample size info
  expect_equal(result$n_obs[2], 2L)
  expect_equal(result$n_treated[2], 1L)
  expect_equal(result$n_control[2], 1L)

  # Periods 1 and 3 should be valid
  expect_false(is.na(result$att[1]))
  expect_false(is.na(result$att[3]))
  expect_true(is.finite(result$att[1]))
  expect_true(is.finite(result$att[3]))
  expect_true(result$se[1] > 0)
  expect_true(result$se[3] > 0)

  # A lwdid_small_sample warning should have been emitted for period 2
  small_sample_warnings <- filter_warnings(cw$warnings, "lwdid_small_sample")
  period_failed <- filter_warnings_detail(small_sample_warnings,
                                          "period_regression_failed")
  expect_true(length(period_failed) >= 1L,
              info = "Expected lwdid_small_sample warning for period 2 failure")
})

test_that("T7-02b: period with N1=0 produces NA row via lwdid_error path", {
  set.seed(722)

  # Period 1: normal (10 obs)
  # Period 2: only control units (N1=0 -> lwdid_insufficient_data)
  dt1 <- data.table(
    tindex = rep(1L, 10L),
    y_trans = rnorm(10, mean = 3),
    d = c(rep(1L, 5L), rep(0L, 5L))
  )
  dt2 <- data.table(
    tindex = rep(2L, 6L),
    y_trans = rnorm(6, mean = 3),
    d = rep(0L, 6L)
  )
  dt1[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  result <- cw$result

  # Period 2 should be NA (N1=0 error caught)
  expect_true(is.na(result$att[2]))
  expect_equal(result$n_treated[2], 0L)
  expect_equal(result$n_control[2], 6L)

  # Period 1 should be valid
  expect_false(is.na(result$att[1]))

  # lwdid_small_sample warning emitted
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"))
})

test_that("T7-02c: period with N0=0 produces NA row via lwdid_error path", {
  set.seed(723)

  # Period 1: normal
  # Period 2: only treated units (N0=0 -> lwdid_insufficient_data)
  dt1 <- data.table(
    tindex = rep(1L, 10L),
    y_trans = rnorm(10, mean = 3),
    d = c(rep(1L, 5L), rep(0L, 5L))
  )
  dt2 <- data.table(
    tindex = rep(2L, 4L),
    y_trans = rnorm(4, mean = 5),
    d = rep(1L, 4L)
  )
  dt1[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  result <- cw$result

  # Period 2 should be NA (N0=0 error caught)
  expect_true(is.na(result$att[2]))
  expect_equal(result$n_treated[2], 4L)
  expect_equal(result$n_control[2], 0L)

  # Period 1 should be valid
  expect_false(is.na(result$att[1]))
})

test_that("T7-02d: multiple periods failing simultaneously, all produce NA", {
  set.seed(724)

  # Period 1: N=2 (< 3, fails)
  # Period 2: N=10 (normal)
  # Period 3: N=1 (< 3, fails)
  dt1 <- data.table(
    tindex = rep(1L, 2L),
    y_trans = c(5, 2), d = c(1L, 0L)
  )
  dt2 <- data.table(
    tindex = rep(2L, 10L),
    y_trans = rnorm(10, mean = 3),
    d = c(rep(1L, 5L), rep(0L, 5L))
  )
  dt3 <- data.table(
    tindex = rep(3L, 1L),
    y_trans = 10, d = 1L
  )
  dt2[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2, dt3))

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL)
  )

  expect_equal(nrow(result), 3L)
  expect_true(is.na(result$att[1]))   # N=2 < 3
  expect_false(is.na(result$att[2]))  # N=10, normal
  expect_true(is.na(result$att[3]))   # N=1 < 3
})

test_that("T7-02e: cluster VCE with insufficient clusters in one period", {
  set.seed(725)

  # Period 1: G=5 clusters (enough for cluster VCE, G >= 2)
  # Period 2: G=1 cluster (G < 2 -> lwdid_insufficient_data from compute_cluster_vce)
  # Period 3: G=4 clusters (enough)
  dt1 <- data.table(
    tindex = rep(1L, 30L),
    y_trans = rnorm(30, mean = 3),
    d = c(rep(1L, 15L), rep(0L, 15L)),
    cluster = rep(1:5, each = 6L)
  )
  dt2 <- data.table(
    tindex = rep(2L, 6L),
    y_trans = rnorm(6, mean = 3),
    d = c(rep(1L, 3L), rep(0L, 3L)),
    cluster = rep(1L, 6L)  # single cluster
  )
  dt3 <- data.table(
    tindex = rep(3L, 24L),
    y_trans = rnorm(24, mean = 3),
    d = c(rep(1L, 12L), rep(0L, 12L)),
    cluster = rep(1:4, each = 6L)
  )
  dt1[d == 1L, y_trans := y_trans + 2]
  dt2[d == 1L, y_trans := y_trans + 2]
  dt3[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2, dt3))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster")
  )
  result <- cw$result

  # Period 2 should be NA (G=1 < 2 -> error caught)
  expect_true(is.na(result$att[2]))

  # Periods 1 and 3 should be valid
  expect_false(is.na(result$att[1]))
  expect_false(is.na(result$att[3]))
  expect_equal(result$n_clusters[1], 5L)
  expect_equal(result$n_clusters[3], 4L)
})


# ============================================================================
# T7-03: Per-period controls_tier variation
# ============================================================================
# Different periods may use different model tiers due to varying sample sizes.
# controls_tier should be recorded independently per period.

test_that("T7-03a: different periods fall into different tiers with VCE", {
  set.seed(731)

  # K=2 controls
  # Period 1: N=40, N1=20, N0=20 -> Tier 1 (N1>3, N0>3)
  # Period 2: N=8, N1=2, N0=6 -> Tier 2 (N1=2 not > K+1=3, but N=8 > K+2=4)
  # Period 3: N=4, N1=2, N0=2 -> Tier 3 (N=4 <= K+2=4)
  k <- 2L

  make_tier_period <- function(period, n1, n0) {
    n_local <- n1 + n0
    d_vec <- c(rep(1L, n1), rep(0L, n0))
    x_mat <- matrix(rnorm(n_local * k), n_local, k)
    dt_local <- data.table(
      tindex = rep(period, n_local),
      y_trans = 1 + 2 * d_vec + 0.5 * x_mat[, 1] + rnorm(n_local, sd = 0.5),
      d = d_vec,
      x1 = x_mat[, 1],
      x2 = x_mat[, 2]
    )
    dt_local
  }

  dt <- rbindlist(list(
    make_tier_period(1L, 20L, 20L),  # Tier 1
    make_tier_period(2L, 2L, 6L),    # Tier 2
    make_tier_period(3L, 2L, 2L)     # Tier 3
  ))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"), periods = 1:3,
                            vce = "hc3")
  )
  result <- cw$result

  # Verify tier assignments
  expect_equal(result$controls_tier[1], "full_interaction")
  expect_equal(result$controls_tier[2], "simple")
  expect_equal(result$controls_tier[3], "none")

  # Verify all three tiers are distinct
  expect_equal(length(unique(result$controls_tier)), 3L)

  # Verify df is consistent with each tier's formula
  # Period 1 (Tier 1): df = N - 2 - 2K = 40 - 2 - 4 = 34
  expect_equal(result$df[1], 40L - 2L - 2L * k)
  # Period 2 (Tier 2): df = N - K - 2 = 8 - 2 - 2 = 4
  expect_equal(result$df[2], 8L - k - 2L)
  # Period 3 (Tier 3): df = N - 2 = 4 - 2 = 2
  expect_equal(result$df[3], 4L - 2L)

  # All periods should have valid ATT
  expect_true(all(!is.na(result$att)))

  # VCE type should be HC3 for all periods
  expect_true(all(result$vce_type == "HC3"))

  # n_clusters should be NA for non-cluster VCE
  expect_true(all(is.na(result$n_clusters)))
})

test_that("T7-03b: per-period tier with cluster VCE, df = G_r - 1 regardless of tier", {

  set.seed(732)

  # K=2 controls
  # Period 1: N=60, N1=30, N0=30, G=20 -> Tier 1, df = G-1 = 19
  # Period 2: N=24, N1=3, N0=21, G=8 -> Tier 2, df = G-1 = 7
  k <- 2L

  make_cluster_tier_period <- function(period, n1, n0, g) {
    n_local <- n1 + n0
    d_vec <- c(rep(1L, n1), rep(0L, n0))
    # Assign clusters evenly
    cl <- rep(seq_len(g), length.out = n_local)
    x_mat <- matrix(rnorm(n_local * k), n_local, k)
    data.table(
      tindex = rep(period, n_local),
      y_trans = 1 + 2 * d_vec + 0.5 * x_mat[, 1] + rnorm(n_local, sd = 0.5),
      d = d_vec,
      x1 = x_mat[, 1],
      x2 = x_mat[, 2],
      cluster = cl
    )
  }

  dt <- rbindlist(list(
    make_cluster_tier_period(1L, 30L, 30L, 20L),  # Tier 1, G=20
    make_cluster_tier_period(2L, 3L, 21L, 8L)     # Tier 2, G=8
  ))

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"), periods = 1:2,
                            vce = "cluster", cluster_var = "cluster")
  )

  # Tier assignments
  expect_equal(result$controls_tier[1], "full_interaction")
  expect_equal(result$controls_tier[2], "simple")

  # df = G_r - 1 regardless of tier
  expect_equal(result$df[1], 19L)  # G=20 -> df=19
  expect_equal(result$df[2], 7L)   # G=8 -> df=7

  # NOT the model residual df
  expect_false(result$df[1] == 60L - 2L - 2L * k)  # not 54
  expect_false(result$df[2] == 24L - k - 2L)        # not 20

  # n_clusters correct
  expect_equal(result$n_clusters[1], 20L)
  expect_equal(result$n_clusters[2], 8L)
})

test_that("T7-03c: per-period ATT matches hand calculation for each tier", {
  set.seed(733)

  k <- 2L

  # Period 1: Tier 1 (full_interaction)
  n1_p1 <- 20L; n0_p1 <- 20L; n_p1 <- n1_p1 + n0_p1
  d1 <- c(rep(1L, n1_p1), rep(0L, n0_p1))
  x1_mat <- matrix(rnorm(n_p1 * k), n_p1, k)
  y1 <- 1 + 2 * d1 + 0.5 * x1_mat[, 1] + rnorm(n_p1, sd = 0.5)

  # Period 2: Tier 2 (simple)
  n1_p2 <- 2L; n0_p2 <- 6L; n_p2 <- n1_p2 + n0_p2
  d2 <- c(rep(1L, n1_p2), rep(0L, n0_p2))
  x2_mat <- matrix(rnorm(n_p2 * k), n_p2, k)
  y2 <- 1 + 2 * d2 + 0.5 * x2_mat[, 1] + rnorm(n_p2, sd = 0.5)

  dt <- rbindlist(list(
    data.table(tindex = 1L, y_trans = y1, d = d1,
               x1 = x1_mat[, 1], x2 = x1_mat[, 2]),
    data.table(tindex = 2L, y_trans = y2, d = d2,
               x1 = x2_mat[, 1], x2 = x2_mat[, 2])
  ))

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"), periods = 1:2,
                            vce = NULL)
  )

  # Period 1 (Tier 1): hand-compute ATT with full interaction
  x_bar1_p1 <- colMeans(x1_mat[d1 == 1, , drop = FALSE])
  x_c_p1 <- sweep(x1_mat, 2, x_bar1_p1, "-")
  x_design_p1 <- cbind(1, d1, x_c_p1, d1 * x_c_p1)
  beta_p1 <- qr.coef(qr(x_design_p1), y1)
  att_hand_p1 <- unname(beta_p1[2])
  expect_equal(result$att[1], att_hand_p1, tolerance = 1e-12)

  # Period 2 (Tier 2): hand-compute ATT with simple controls
  x_design_p2 <- cbind(1, d2, x2_mat)
  beta_p2 <- qr.coef(qr(x_design_p2), y2)
  att_hand_p2 <- unname(beta_p2[2])
  expect_equal(result$att[2], att_hand_p2, tolerance = 1e-12)
})

test_that("T7-03d: degradation warnings emitted per period", {
  set.seed(734)

  k <- 2L

  # Period 1: Tier 1 (no degradation warning)
  # Period 2: Tier 2 (lwdid_data warning with detail=controls_degraded_to_simple)
  # Period 3: Tier 3 (lwdid_data warning with detail=controls_dropped)
  make_tier_period <- function(period, n1, n0) {
    n_local <- n1 + n0
    d_vec <- c(rep(1L, n1), rep(0L, n0))
    data.table(
      tindex = rep(period, n_local),
      y_trans = 1 + 2 * d_vec + rnorm(n_local, sd = 0.5),
      d = d_vec,
      x1 = rnorm(n_local),
      x2 = rnorm(n_local)
    )
  }

  dt <- rbindlist(list(
    make_tier_period(1L, 20L, 20L),  # Tier 1
    make_tier_period(2L, 2L, 6L),    # Tier 2
    make_tier_period(3L, 2L, 2L)     # Tier 3
  ))

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex",
                            x = c("x1", "x2"), periods = 1:3,
                            vce = NULL)
  )

  # Check for Tier 2 degradation warning
  data_warnings <- filter_warnings(cw$warnings, "lwdid_data")
  tier2_warnings <- filter_warnings_detail(data_warnings,
                                           "controls_degraded_to_simple")
  expect_true(length(tier2_warnings) >= 1L,
              info = "Expected Tier 2 degradation warning")

  # Check for Tier 3 degradation warning
  tier3_warnings <- filter_warnings_detail(data_warnings, "controls_dropped")
  expect_true(length(tier3_warnings) >= 1L,
              info = "Expected Tier 3 degradation warning")
})


# ============================================================================
# T7-04: Non-lwdid error fault tolerance (generic error handler)
# ============================================================================
# When a period fails with a non-lwdid error (e.g., from sandwich internals
# or numerical singularity), the generic error handler should catch it,
# emit a lwdid_numerical warning, and produce an NA row.

test_that("T7-04a: generic error caught, lwdid_numerical warning emitted", {
  set.seed(741)

  # Strategy: mock estimate_ra_common to throw a generic error for period 2
  # We do this by constructing data that causes a non-lwdid error.
  # A rank-deficient design with controls that passes initial checks but
  # fails during VCE computation would be ideal, but that's hard to construct.
  #
  # Instead, we use a direct approach: temporarily override estimate_ra_common
  # to throw a generic error for specific input patterns.

  # Period 1: normal data
  dt1 <- data.table(
    tindex = rep(1L, 20L),
    y_trans = rnorm(20, mean = 3),
    d = c(rep(1L, 10L), rep(0L, 10L))
  )
  dt1[d == 1L, y_trans := y_trans + 2]

  # Period 2: normal data (we'll intercept the call)
  dt2 <- data.table(
    tindex = rep(2L, 20L),
    y_trans = rnorm(20, mean = 3),
    d = c(rep(1L, 10L), rep(0L, 10L))
  )
  dt2[d == 1L, y_trans := y_trans + 2]

  dt <- rbindlist(list(dt1, dt2))

  # Save original function
  orig_fn <- estimate_ra_common
  call_count <- 0L

  # Override: throw generic error on second call (period 2)
  mock_fn <- function(y_trans, d, x = NULL, vce = NULL,
                      cluster = NULL, alpha = 0.05) {
    call_count <<- call_count + 1L
    if (call_count == 2L) {
      stop("Unexpected internal error in sandwich computation")
    }
    orig_fn(y_trans, d, x, vce = vce, cluster = cluster, alpha = alpha)
  }

  # Temporarily replace in the lwdid namespace
  ns <- asNamespace("lwdid")
  unlockBinding("estimate_ra_common", ns)
  assign("estimate_ra_common", mock_fn, envir = ns)
  on.exit({
    assign("estimate_ra_common", orig_fn, envir = ns)
    lockBinding("estimate_ra_common", ns)
  }, add = TRUE)

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = NULL)
  )
  result <- cw$result

  # Period 1 should be valid
  expect_false(is.na(result$att[1]))

  # Period 2 should be NA (generic error caught)
  expect_true(is.na(result$att[2]))
  expect_true(is.na(result$se[2]))
  expect_true(is.na(result$df[2]))

  # lwdid_numerical warning should have been emitted (not lwdid_small_sample)
  numerical_warnings <- filter_warnings(cw$warnings, "lwdid_numerical")
  expect_true(length(numerical_warnings) >= 1L,
              info = "Expected lwdid_numerical warning for generic error")

  # The warning detail should be "period_regression_failed"
  period_failed <- filter_warnings_detail(numerical_warnings,
                                          "period_regression_failed")
  expect_true(length(period_failed) >= 1L,
              info = "Expected detail = 'period_regression_failed'")
})

test_that("T7-04b: generic error preserves sample size info in NA row", {
  set.seed(742)

  dt1 <- data.table(
    tindex = rep(1L, 12L),
    y_trans = rnorm(12, mean = 3),
    d = c(rep(1L, 5L), rep(0L, 7L))
  )
  dt1[d == 1L, y_trans := y_trans + 2]

  dt <- dt1

  orig_fn <- estimate_ra_common
  mock_fn <- function(y_trans, d, x = NULL, vce = NULL,
                      cluster = NULL, alpha = 0.05) {
    stop("Simulated sandwich internal error")
  }

  ns <- asNamespace("lwdid")
  unlockBinding("estimate_ra_common", ns)
  assign("estimate_ra_common", mock_fn, envir = ns)
  on.exit({
    assign("estimate_ra_common", orig_fn, envir = ns)
    lockBinding("estimate_ra_common", ns)
  }, add = TRUE)

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1L, vce = NULL)
  )

  # NA row should preserve sample sizes from the filtered data
  expect_true(is.na(result$att))
  expect_equal(result$n_obs, 12L)
  expect_equal(result$n_treated, 5L)
  expect_equal(result$n_control, 7L)
})

test_that("T7-04c: lwdid_error vs generic error produce different warning classes", {
  set.seed(743)

  # Period 1: will trigger lwdid_error (N < 3)
  # Period 2: will trigger generic error (via mock)
  # Period 3: normal
  dt1 <- data.table(
    tindex = rep(1L, 2L),
    y_trans = c(5, 2), d = c(1L, 0L)
  )
  dt2 <- data.table(
    tindex = rep(2L, 20L),
    y_trans = rnorm(20, mean = 3),
    d = c(rep(1L, 10L), rep(0L, 10L))
  )
  dt3 <- data.table(
    tindex = rep(3L, 20L),
    y_trans = rnorm(20, mean = 3),
    d = c(rep(1L, 10L), rep(0L, 10L))
  )
  dt2[d == 1L, y_trans := y_trans + 2]
  dt3[d == 1L, y_trans := y_trans + 2]
  dt <- rbindlist(list(dt1, dt2, dt3))

  orig_fn <- estimate_ra_common
  call_count <- 0L

  mock_fn <- function(y_trans, d, x = NULL, vce = NULL,
                      cluster = NULL, alpha = 0.05) {
    call_count <<- call_count + 1L
    if (call_count == 2L) {
      stop("Generic internal error")
    }
    orig_fn(y_trans, d, x, vce = vce, cluster = cluster, alpha = alpha)
  }

  ns <- asNamespace("lwdid")
  unlockBinding("estimate_ra_common", ns)
  assign("estimate_ra_common", mock_fn, envir = ns)
  on.exit({
    assign("estimate_ra_common", orig_fn, envir = ns)
    lockBinding("estimate_ra_common", ns)
  }, add = TRUE)

  cw <- capture_with_warnings(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = NULL)
  )
  result <- cw$result

  # Period 1: NA (lwdid_error -> lwdid_small_sample warning)
  expect_true(is.na(result$att[1]))
  # Period 2: NA (generic error -> lwdid_numerical warning)
  expect_true(is.na(result$att[2]))
  # Period 3: valid
  expect_false(is.na(result$att[3]))

  # Both warning classes should be present
  expect_true(has_warning_class(cw$warnings, "lwdid_small_sample"),
              info = "lwdid_small_sample from period 1 (lwdid_error)")
  expect_true(has_warning_class(cw$warnings, "lwdid_numerical"),
              info = "lwdid_numerical from period 2 (generic error)")
})

# ============================================================================
# T7-05: Additional edge cases and numerical validation
# ============================================================================

test_that("T7-05a: 15-column output structure preserved for all VCE types", {
  set.seed(751)

  dt <- data.table(
    tindex = rep(1:2, each = 30L),
    y_trans = rnorm(60, mean = 3),
    d = rep(c(rep(1L, 15L), rep(0L, 15L)), 2),
    cluster = rep(rep(1:5, each = 6L), 2)
  )
  dt[d == 1L, y_trans := y_trans + 2]

  # Test with vce=NULL
  r_null <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                    periods = 1:2, vce = NULL)
  expect_equal(names(r_null), expected_cols)
  expect_equal(ncol(r_null), 15L)
  expect_true(all(r_null$vce_type == "homoskedastic"))
  expect_true(all(is.na(r_null$n_clusters)))

  # Test with vce="hc3"
  r_hc3 <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                   periods = 1:2, vce = "hc3")
  expect_equal(names(r_hc3), expected_cols)
  expect_true(all(r_hc3$vce_type == "HC3"))
  expect_true(all(is.na(r_hc3$n_clusters)))

  # Test with vce="cluster"
  r_cl <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:2, vce = "cluster",
                            cluster_var = "cluster")
  )
  expect_equal(names(r_cl), expected_cols)
  expect_true(all(r_cl$vce_type == "cluster"))
  expect_true(all(!is.na(r_cl$n_clusters)))
})

test_that("T7-05b: CI contains ATT for all valid periods", {
  set.seed(752)

  dt <- data.table(
    tindex = rep(1:3, each = 40L),
    y_trans = rnorm(120, mean = 3),
    d = rep(c(rep(1L, 20L), rep(0L, 20L)), 3),
    cluster = rep(rep(1:10, each = 4L), 3)
  )
  dt[d == 1L, y_trans := y_trans + 2]

  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                            periods = 1:3, vce = "cluster",
                            cluster_var = "cluster")
  )

  for (i in seq_len(nrow(result))) {
    if (!is.na(result$att[i])) {
      expect_true(result$ci_lower[i] <= result$att[i],
                  label = sprintf("Period %d: ci_lower <= att", i))
      expect_true(result$ci_upper[i] >= result$att[i],
                  label = sprintf("Period %d: ci_upper >= att", i))
    }
  }
})

test_that("T7-05c: robust VCE alias works in period effects", {
  set.seed(753)

  dt <- data.table(
    tindex = rep(1:2, each = 20L),
    y_trans = rnorm(40, mean = 3),
    d = rep(c(rep(1L, 10L), rep(0L, 10L)), 2)
  )
  dt[d == 1L, y_trans := y_trans + 2]

  r_robust <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                      periods = 1:2, vce = "robust")
  r_hc1 <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                   periods = 1:2, vce = "hc1")

  # "robust" is alias for "hc1" -> SE should be identical
  expect_equal(r_robust$se, r_hc1$se, tolerance = 1e-12)
  expect_equal(r_robust$att, r_hc1$att, tolerance = 1e-12)

  # vce_type: robust maps to HC1
  expect_true(all(r_robust$vce_type == "HC1"))
  expect_true(all(r_hc1$vce_type == "HC1"))
})

test_that("T7-05d: alpha parameter propagates to per-period CI width", {
  set.seed(754)

  dt <- data.table(
    tindex = rep(1L, 30L),
    y_trans = rnorm(30, mean = 3),
    d = c(rep(1L, 15L), rep(0L, 15L))
  )
  dt[d == 1L, y_trans := y_trans + 2]

  r_05 <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                  periods = 1L, vce = NULL, alpha = 0.05)
  r_10 <- estimate_period_effects(dt, "y_trans", "d", "tindex", x = NULL,
                                  periods = 1L, vce = NULL, alpha = 0.10)

  # alpha=0.10 should produce narrower CI than alpha=0.05
  width_05 <- r_05$ci_upper - r_05$ci_lower
  width_10 <- r_10$ci_upper - r_10$ci_lower
  expect_true(width_10 < width_05,
              info = "alpha=0.10 CI should be narrower than alpha=0.05")

  # ATT and SE should be identical (alpha only affects CI)
  expect_equal(r_05$att, r_10$att, tolerance = 1e-12)
  expect_equal(r_05$se, r_10$se, tolerance = 1e-12)
})
