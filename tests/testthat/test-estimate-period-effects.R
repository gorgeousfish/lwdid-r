# test-estimate-period-effects.R
# Tests for estimate_period_effects() and helper functions
# Covers Tasks E2-05.8 through E2-05.15

# =============================================================================
# Helper: create balanced panel with known DGP
# =============================================================================
make_period_panel <- function(N = 10L, TT = 6L, S = 4L, n_treated = 5L,
                              tau_fn = function(t, S) t - S + 1L,
                              sd_noise = 0.5, seed = 42L,
                              add_controls = FALSE, n_controls = 1L) {
  alpha_i <- seq(10, by = 5, length.out = N)
  delta_t <- seq(0, by = 2, length.out = TT)
  set.seed(seed)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, post := as.integer(time >= S)]
  dt[, tau_t := ifelse(post == 1L, tau_fn(time, S), 0)]
  dt[, y := alpha_i[id] + delta_t[time] + tau_t * d + rnorm(.N, sd = sd_noise)]
  if (add_controls) {
    for (k in seq_len(n_controls)) {
      col_name <- paste0("x", k)
      dt[, (col_name) := rep(rnorm(N, sd = 2), each = TT)]
    }
  }
  dt
}

# Helper: transform and get period effects in one call
run_period_effects <- function(dt, S, TT, controls = NULL, alpha = 0.05,
                               include_pretreatment = FALSE) {
  dt_trans <- suppressWarnings(
    transform_demean(dt, "y", "id", "time", g = S)
  )
  estimate_period_effects(dt_trans, "y_trans", "d", "time",
                          x = controls, periods = S:TT,
                          alpha = alpha,
                          include_pretreatment = include_pretreatment)
}

# Helper: suppress lwdid warnings
quiet_lwdid <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

# =============================================================================
# Group A: Basic period estimation correctness (T-01 to T-08) — Task E2-05.8
# Data: 10 units x 6 periods, S=4, 5 treated + 5 control, seed=42
# =============================================================================

test_that("T-01 period ATTs match hand-computed mean differences", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  for (i in seq_len(nrow(result))) {
    r <- result$period[i]
    sub <- dt_trans[time == r]
    expected_att <- mean(sub[d == 1]$y_trans) - mean(sub[d == 0]$y_trans)
    expect_equal(result$att[i], expected_att, tolerance = 1e-12,
                 label = sprintf("ATT period %d", r))
  }
})

test_that("T-02 output data.frame has all 13 columns", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  expected_cols <- c("tindex", "period", "att", "se", "t_stat", "pvalue",
                     "ci_lower", "ci_upper", "n_obs", "n_treated", "n_control",
                     "df", "vce_type", "n_clusters", "controls_tier")
  expect_equal(names(result), expected_cols)
  expect_equal(ncol(result), 15L)
})

test_that("T-03 tindex and period columns match input periods", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  expect_equal(result$tindex, 4:6)
  expect_equal(result$period, 4:6)
})

test_that("T-04 period SEs match hand-computed OLS SE", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  for (i in seq_len(nrow(result))) {
    r <- result$period[i]
    sub <- dt_trans[time == r]
    y_r <- sub$y_trans; d_r <- sub$d
    X <- cbind(1, d_r)
    XtX_inv <- solve(t(X) %*% X)
    beta <- XtX_inv %*% t(X) %*% y_r
    resid <- y_r - X %*% beta
    sigma2 <- sum(resid^2) / (length(y_r) - 2L)
    se_r <- sqrt(sigma2 * XtX_inv[2, 2])
    expect_equal(result$se[i], se_r, tolerance = 1e-12,
                 label = sprintf("SE period %d", r))
  }
})

test_that("T-05 period df = N_t - 2 (no controls)", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  expect_equal(result$df, rep(8L, 3))
})

test_that("T-06 period CIs use t-distribution critical value", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  for (i in seq_len(nrow(result))) {
    tcrit <- qt(0.975, result$df[i])
    expect_equal(result$ci_lower[i], result$att[i] - tcrit * result$se[i],
                 tolerance = 1e-12)
    expect_equal(result$ci_upper[i], result$att[i] + tcrit * result$se[i],
                 tolerance = 1e-12)
  }
})

test_that("T-07 period t_stat = ATT / SE", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  for (i in seq_len(nrow(result))) {
    expect_equal(result$t_stat[i], result$att[i] / result$se[i],
                 tolerance = 1e-12)
  }
})

test_that("T-08 period pvalue = 2*pt(-|t|, df)", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  for (i in seq_len(nrow(result))) {
    expected_p <- 2 * pt(-abs(result$t_stat[i]), result$df[i])
    expect_equal(result$pvalue[i], expected_p, tolerance = 1e-10)
  }
})

# =============================================================================
# Group B: Summary ATT equivalence (T-09 to T-11) — Task E2-05.9
# =============================================================================

test_that("T-09 balanced panel no controls: mean(period ATTs) = summary ATT", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  period_result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                            x = NULL, periods = 4:6)
  mean_period_att <- mean(period_result$att)
  # Summary ATT: post-period mean regression
  post_dt <- dt_trans[time >= 4L]
  y_bar <- post_dt[, .(y_bar = mean(y_trans)), by = id]
  d_map <- unique(dt[, .(id, d)])
  y_bar <- merge(y_bar, d_map, by = "id")
  summary_att <- mean(y_bar[d == 1]$y_bar) - mean(y_bar[d == 0]$y_bar)
  expect_equal(mean_period_att, summary_att, tolerance = 1e-12)
})

test_that("T-10 unbalanced panel no controls: mean(period ATTs) != summary ATT", {
  dt <- make_period_panel(N = 10L, TT = 6L, S = 4L, n_treated = 5L, seed = 42L)
  # Remove treated unit 1 from period 5, treated unit 2 from period 6
  dt <- dt[!(id == 1L & time == 5L)]
  dt <- dt[!(id == 2L & time == 6L)]
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  period_result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                            x = NULL, periods = 4:6)
  mean_period_att <- mean(period_result$att)
  post_dt <- dt_trans[time >= 4L]
  y_bar <- post_dt[, .(y_bar = mean(y_trans)), by = id]
  d_map <- unique(dt[, .(id, d)])
  y_bar <- merge(y_bar, d_map, by = "id")
  summary_att <- mean(y_bar[d == 1]$y_bar) - mean(y_bar[d == 0]$y_bar)
  expect_false(isTRUE(all.equal(mean_period_att, summary_att, tolerance = 1e-8)))
})

test_that("T-11 balanced panel with controls: mean(period ATTs) != summary ATT", {
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 1L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  # Period ATTs with controls
  period_result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = "x1", periods = 4:6)
  )
  mean_period_att <- mean(period_result$att)
  # Summary ATT: post-period mean regression WITHOUT controls (simple mean diff)
  # vs period ATTs WITH controls — these should differ
  post_dt <- dt_trans[time >= 4L]
  y_bar <- post_dt[, .(y_bar = mean(y_trans)), by = id]
  d_map <- unique(dt[, .(id, d)])
  y_bar <- merge(y_bar, d_map, by = "id")
  summary_att_no_ctrl <- mean(y_bar[d == 1]$y_bar) - mean(y_bar[d == 0]$y_bar)
  # Period ATTs with controls should differ from summary ATT without controls
  expect_false(isTRUE(all.equal(mean_period_att, summary_att_no_ctrl, tolerance = 1e-8)))
})

# =============================================================================
# Group C: Controls and per-period fallback (T-12 to T-16) — Task E2-05.10
# =============================================================================

test_that("T-12 with controls Tier 1: each period regresses independently", {
  # 20 units, 2 controls, enough for Tier 1 (N1=10 > K+1=3, N0=10 > K+1=3)
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 2L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = c("x1", "x2"), periods = 4:6)
  )

  # All periods should have valid ATT (not NA)
  expect_true(all(!is.na(result$att)))
  # All periods should use full_interaction tier
  expect_equal(result$controls_tier, rep("full_interaction", 3))

  # Verify each period's ATT matches independent hand calculation
  # using the Tier 1 design matrix [1, D, X_c, D*X_c]
  K <- 2L
  for (i in seq_len(nrow(result))) {
    r <- result$period[i]
    sub <- dt_trans[time == r]
    y_r <- sub$y_trans
    d_r <- sub$d
    x_r <- as.matrix(sub[, c("x1", "x2"), with = FALSE])

    # Tier 1: center at treated group mean
    x_bar1 <- colMeans(x_r[d_r == 1, , drop = FALSE])
    x_c <- sweep(x_r, 2, x_bar1, "-")
    X_design <- cbind(1, d_r, x_c, d_r * x_c)

    qr_fit <- qr(X_design)
    beta <- qr.coef(qr_fit, y_r)
    att_hand <- unname(beta[2])

    expect_equal(result$att[i], att_hand, tolerance = 1e-12,
                 label = sprintf("Tier 1 ATT period %d", r))
  }

  # Verify control coefficients differ across periods (independence)
  # DGP noise ensures different periods have different OLS coefficients
  coefs_by_period <- lapply(4:6, function(r) {
    sub <- dt_trans[time == r]
    y_r <- sub$y_trans
    d_r <- sub$d
    x_r <- as.matrix(sub[, c("x1", "x2"), with = FALSE])
    x_bar1 <- colMeans(x_r[d_r == 1, , drop = FALSE])
    x_c <- sweep(x_r, 2, x_bar1, "-")
    X_design <- cbind(1, d_r, x_c, d_r * x_c)
    qr.coef(qr(X_design), y_r)
  })
  # x1 main-effect coefficient (position 3) should differ across periods
  gamma1_coefs <- vapply(coefs_by_period, function(b) b[3], numeric(1))
  expect_false(
    isTRUE(all.equal(gamma1_coefs, rep(gamma1_coefs[1], 3), tolerance = 1e-6)),
    label = "x1 coefficients should differ across periods"
  )
})

test_that("T-13 with controls Tier 1: df = N_t - 2 - 2K", {
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 2L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = c("x1", "x2"), periods = 4:6)
  )

  # K=2, N_t=20 per period (balanced), df = 20 - 2 - 2*2 = 14
  expect_equal(result$df, rep(14L, 3))
})

test_that("T-14 controls_tier values are valid strings", {
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 2L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = c("x1", "x2"), periods = 4:6)
  )

  valid_tiers <- c("full_interaction", "simple", "none", NA_character_)
  expect_true(all(result$controls_tier %in% valid_tiers))
})

test_that("T-15 unbalanced panel: different periods fall into different tiers", {
  # K=2 controls: Tier 1 needs N1>3 AND N0>3; Tier 2 needs N>4; Tier 3: N<=4
  # Start with 20 units (10 treated ids 1-10, 10 control ids 11-20)
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 2L)

  # Period 4: keep all 20 units -> Tier 1 (N1=10>3, N0=10>3)
  # Period 5: keep 3 treated + 4 control = 7 units -> Tier 2
  #   N1=3 not > K+1=3, so NOT Tier 1; N=7 > K+2=4, so Tier 2
  dt <- dt[!(id %in% 4:10 & time == 5L)]   # remove treated 4-10 from period 5
  dt <- dt[!(id %in% 15:20 & time == 5L)]  # remove control 15-20 from period 5

  # Period 6: keep 2 treated + 2 control = 4 units -> Tier 3
  #   N=4 <= K+2=4, so Tier 3
  dt <- dt[!(id %in% 3:10 & time == 6L)]   # remove treated 3-10 from period 6
  dt <- dt[!(id %in% 13:20 & time == 6L)]  # remove control 13-20 from period 6

  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))

  # Verify unit counts per period
  for (r in 4:6) {
    sub <- dt_trans[time == r]
    n1_r <- sum(sub$d == 1)
    n0_r <- sum(sub$d == 0)
    if (r == 4L) {
      expect_equal(n1_r, 10L, label = "period 4: N1=10")
      expect_equal(n0_r, 10L, label = "period 4: N0=10")
    } else if (r == 5L) {
      expect_equal(n1_r, 3L, label = "period 5: N1=3")
      expect_equal(n0_r, 4L, label = "period 5: N0=4")
    } else if (r == 6L) {
      expect_equal(n1_r, 2L, label = "period 6: N1=2")
      expect_equal(n0_r, 2L, label = "period 6: N0=2")
    }
  }

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = c("x1", "x2"), periods = 4:6)
  )

  # Period 4: Tier 1 (full_interaction)
  expect_equal(result$controls_tier[1], "full_interaction")
  # Period 5: Tier 2 (simple) — N1=3 not > K+1=3
  expect_equal(result$controls_tier[2], "simple")
  # Period 6: Tier 3 (none) — N=4 <= K+2=4
  expect_equal(result$controls_tier[3], "none")

  # Verify tiers are not all the same
  expect_true(length(unique(result$controls_tier)) == 3L,
              label = "all three tiers should be present")

  # Verify df is consistent with each tier's formula
  K <- 2L
  # Period 4 (Tier 1): df = N_t - 2 - 2K = 20 - 6 = 14
  expect_equal(result$df[1], 20L - 2L - 2L * K)
  # Period 5 (Tier 2): df = N_t - K - 2 = 7 - 2 - 2 = 3
  expect_equal(result$df[2], 7L - K - 2L)
  # Period 6 (Tier 3): df = N_t - 2 = 4 - 2 = 2
  expect_equal(result$df[3], 4L - 2L)

  # Verify ATT for each period matches hand calculation with appropriate tier
  for (i in seq_len(nrow(result))) {
    r <- result$period[i]
    sub <- dt_trans[time == r]
    y_r <- sub$y_trans
    d_r <- sub$d
    x_r <- as.matrix(sub[, c("x1", "x2"), with = FALSE])
    tier <- result$controls_tier[i]

    if (tier == "full_interaction") {
      x_bar1 <- colMeans(x_r[d_r == 1, , drop = FALSE])
      x_c <- sweep(x_r, 2, x_bar1, "-")
      X_design <- cbind(1, d_r, x_c, d_r * x_c)
    } else if (tier == "simple") {
      X_design <- cbind(1, d_r, x_r)
    } else {
      X_design <- cbind(1, d_r)
    }

    qr_fit <- qr(X_design)
    beta <- qr.coef(qr_fit, y_r)
    att_hand <- unname(beta[2])

    expect_equal(result$att[i], att_hand, tolerance = 1e-12,
                 label = sprintf("period %d ATT (tier=%s)", r, tier))
  }
})

test_that("T-16 per-period X_bar_1 computed independently", {
  # Unbalanced panel: different treated units present in different periods
  # -> per-period X_bar_1 should differ, producing different ATTs
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 1L)

  # Remove treated units 1-3 from period 5 only
  # Period 4: treated units 1-10, Period 5: treated units 4-10
  # Both periods have N1>2 and N0>2 -> Tier 1 (K=1, need N1>2 and N0>2)
  dt <- dt[!(id %in% 1:3 & time == 5L)]

  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))

  # Verify both periods qualify for Tier 1 (K=1, need N1>2 and N0>2)
  sub4 <- dt_trans[time == 4L]
  sub5 <- dt_trans[time == 5L]
  expect_true(sum(sub4$d == 1) > 2L && sum(sub4$d == 0) > 2L,
              label = "period 4 qualifies for Tier 1")
  expect_true(sum(sub5$d == 1) > 2L && sum(sub5$d == 0) > 2L,
              label = "period 5 qualifies for Tier 1")

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = "x1", periods = 4:5)
  )

  # Both periods should produce valid results with full_interaction
  expect_true(all(!is.na(result$att)))
  expect_equal(result$controls_tier, rep("full_interaction", 2))

  # Verify X_bar_1 differs between periods
  x4_bar1 <- mean(sub4[d == 1]$x1)
  x5_bar1 <- mean(sub5[d == 1]$x1)
  expect_false(
    isTRUE(all.equal(x4_bar1, x5_bar1, tolerance = 1e-6)),
    label = "X_bar_1 should differ between periods"
  )

  # Compute period 5 ATT with SHARED X_bar_1 (from period 4's treated units)
  y_5 <- sub5$y_trans
  d_5 <- sub5$d
  x_5 <- as.matrix(sub5[, "x1", with = FALSE])

  # Shared X_bar_1 (wrong: from period 4's treated)
  x_bar1_shared <- mean(sub4[d == 1]$x1)
  x_c_shared <- x_5 - x_bar1_shared
  X_shared <- cbind(1, d_5, x_c_shared, d_5 * x_c_shared)
  beta_shared <- qr.coef(qr(X_shared), y_5)
  att_shared <- unname(beta_shared[2])

  # Period-specific X_bar_1 (correct: from period 5's treated)
  x_bar1_period <- mean(sub5[d == 1]$x1)
  x_c_period <- x_5 - x_bar1_period
  X_period <- cbind(1, d_5, x_c_period, d_5 * x_c_period)
  beta_period <- qr.coef(qr(X_period), y_5)
  att_period <- unname(beta_period[2])

  # Actual result should match period-specific calculation
  expect_equal(result$att[2], att_period, tolerance = 1e-12,
               label = "period 5 ATT matches period-specific X_bar_1")

  # Actual result should NOT match shared X_bar_1 calculation
  expect_false(
    isTRUE(all.equal(result$att[2], att_shared, tolerance = 1e-6)),
    label = "period 5 ATT differs from shared X_bar_1"
  )
})

# =============================================================================
# Group D: NA filtering (T-17 to T-20) — Task E2-05.11
# =============================================================================

test_that("T-17 y_trans NA observations are excluded", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  dt_trans[id %in% c(1L, 6L) & time == 4L, y_trans := NA_real_]
  result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                     x = NULL, periods = 4:6)
  expect_equal(result$n_obs[result$period == 4L], 8L)
  expect_equal(result$n_obs[result$period == 5L], 10L)
  expect_equal(result$n_obs[result$period == 6L], 10L)
})

test_that("T-18 d NA observations are excluded", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  dt_trans[id == 3L & time == 5L, d := NA_integer_]
  result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                     x = NULL, periods = 4:6)
  expect_equal(result$n_obs[result$period == 5L], 9L)
  expect_equal(result$n_obs[result$period == 4L], 10L)
})

test_that("T-19 control variable NA observations are excluded", {
  dt <- make_period_panel(N = 20L, TT = 6L, S = 4L, n_treated = 10L,
                          seed = 42L, add_controls = TRUE, n_controls = 1L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  dt_trans[id %in% c(2L, 5L, 15L) & time == 6L, x1 := NA_real_]
  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = "x1", periods = 4:6)
  )
  expect_equal(result$n_obs[result$period == 6L], 17L)
  expect_equal(result$n_obs[result$period == 4L], 20L)
})

test_that("T-20 all observations NA in a period returns all-NA row", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  dt_trans[time == 5L, y_trans := NA_real_]
  result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                     x = NULL, periods = 4:6)
  row5 <- result[result$period == 5L, ]
  expect_true(is.na(row5$att))
  expect_true(is.na(row5$se))
  expect_equal(row5$n_obs, 0L)
  expect_equal(row5$n_treated, 0L)
  expect_equal(row5$n_control, 0L)
  expect_true(is.na(row5$df))
  expect_true(is.na(row5$controls_tier))
  expect_false(is.na(result$att[result$period == 4L]))
  expect_false(is.na(result$att[result$period == 6L]))
})

# =============================================================================
# Group E: Fault tolerance (T-21 to T-29) — Task E2-05.12
# =============================================================================

test_that("T-21 period with N1=0 triggers tryCatch + lwdid_small_sample warning", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  # Remove all treated units from period 5
  dt_trans <- dt_trans[!(d == 1L & time == 5L)]

  w <- NULL
  result <- withCallingHandlers(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = 4:6),
    lwdid_small_sample = function(cnd) {
      if (!is.null(cnd$detail) && cnd$detail == "period_regression_failed") {
        w <<- cnd
      }
      invokeRestart("muffleWarning")
    }
  )

  row5 <- result[result$period == 5L, ]
  expect_true(is.na(row5$att))
  expect_equal(row5$n_treated, 0L)
  expect_false(is.null(w))
  expect_equal(w$detail, "period_regression_failed")
})

test_that("T-22 period with N0=0 triggers tryCatch + lwdid_small_sample warning", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  # Remove all control units from period 6
  dt_trans <- dt_trans[!(d == 0L & time == 6L)]

  w <- NULL
  result <- withCallingHandlers(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = 4:6),
    lwdid_small_sample = function(cnd) {
      if (!is.null(cnd$detail) && cnd$detail == "period_regression_failed") {
        w <<- cnd
      }
      invokeRestart("muffleWarning")
    }
  )

  row6 <- result[result$period == 6L, ]
  expect_true(is.na(row6$att))
  expect_equal(row6$n_control, 0L)
  expect_false(is.null(w))
  expect_equal(w$detail, "period_regression_failed")
})

test_that("T-23 perfect separation triggers degenerate SE path", {
  # Perfect separation: all treated have same y_trans, all control have same y_trans
  # Residuals are exactly 0, SE is exactly 0 (< 1e-10 threshold)
  # estimate_ra_common detects degenerate SE, returns se=0 with inference NAs
  # Period-level degenerate check catches se==0, returns NA row
  dt <- data.table::data.table(
    id = 1:6, time = rep(4L, 6),
    y_trans = c(5, 5, 5, 2, 2, 2),
    d = c(1L, 1L, 1L, 0L, 0L, 0L)
  )

  w <- NULL
  result <- withCallingHandlers(
    estimate_period_effects(dt, "y_trans", "d", "time",
                            x = NULL, periods = 4L),
    lwdid_small_sample = function(cnd) {
      if (!is.null(cnd$detail) && cnd$detail == "degenerate_period_regression") {
        w <<- cnd
      }
      invokeRestart("muffleWarning")
    }
  )

  # Degenerate SE -> ATT set to NA at period level
  expect_true(is.na(result$att))
  expect_true(is.na(result$se))
  # Sample sizes preserved in degenerate row
  expect_equal(result$n_obs, 6L)
  expect_equal(result$n_treated, 3L)
  expect_equal(result$n_control, 3L)
  # Degenerate warning was emitted
  expect_false(is.null(w))
  expect_equal(w$detail, "degenerate_period_regression")
})

test_that("T-24 N=2 period triggers tryCatch (N<3 error)", {
  # N=2 triggers estimate_ra_common's N<3 check -> lwdid_error -> inner tryCatch
  dt <- data.table::data.table(
    id = 1:2, time = rep(4L, 2),
    y_trans = c(5, 2),
    d = c(1L, 0L)
  )

  w <- NULL
  result <- withCallingHandlers(
    estimate_period_effects(dt, "y_trans", "d", "time",
                            x = NULL, periods = 4L),
    lwdid_small_sample = function(cnd) {
      if (!is.null(cnd$detail) && cnd$detail == "period_regression_failed") {
        w <<- cnd
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(is.na(result$att))
  expect_true(is.na(result$se))
  # Sample sizes from filtered data (post-error)
  expect_equal(result$n_obs, 2L)
  expect_equal(result$n_treated, 1L)
  expect_equal(result$n_control, 1L)
  expect_false(is.null(w))
  expect_equal(w$detail, "period_regression_failed")
})

test_that("T-25 N=1 period triggers tryCatch (N<3 error)", {
  # Single observation -> N<3 error -> lwdid_error -> inner tryCatch
  dt <- data.table::data.table(
    id = 1L, time = 5L,
    y_trans = 10,
    d = 1L
  )

  w <- NULL
  result <- withCallingHandlers(
    estimate_period_effects(dt, "y_trans", "d", "time",
                            x = NULL, periods = 5L),
    lwdid_small_sample = function(cnd) {
      if (!is.null(cnd$detail) && cnd$detail == "period_regression_failed") {
        w <<- cnd
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_true(is.na(result$att))
  expect_equal(result$n_obs, 1L)
  expect_false(is.null(w))
  expect_equal(w$detail, "period_regression_failed")
})

test_that("T-26 individual period failure does not affect other periods", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  # Remove all treated from period 5 -> period 5 fails
  dt_trans <- dt_trans[!(d == 1L & time == 5L)]

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = 4:6)
  )

  expect_true(is.na(result$att[result$period == 5L]))
  expect_false(is.na(result$att[result$period == 4L]))
  expect_false(is.na(result$att[result$period == 6L]))
  expect_true(is.finite(result$att[result$period == 4L]))
  expect_true(is.finite(result$att[result$period == 6L]))
})

test_that("T-27 multiple periods failing simultaneously", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  # Remove all treated from periods 4 and 6
  dt_trans <- dt_trans[!(d == 1L & time %in% c(4L, 6L))]

  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = 4:6)
  )

  expect_true(is.na(result$att[result$period == 4L]))
  expect_true(is.na(result$att[result$period == 6L]))
  expect_false(is.na(result$att[result$period == 5L]))
  expect_equal(nrow(result), 3L)
})

test_that("T-28 include_pretreatment=FALSE (default) works normally", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L,
                                            include_pretreatment = FALSE))
  expect_equal(nrow(result), 3L)
  expect_true(all(!is.na(result$att)))
})

test_that("T-29 include_pretreatment=TRUE produces same results (stub)", {
  dt <- make_period_panel(seed = 42L)
  result_false <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L,
                                                  include_pretreatment = FALSE))
  result_true <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L,
                                                 include_pretreatment = TRUE))
  expect_equal(result_true$att, result_false$att, tolerance = 1e-15)
  expect_equal(result_true$se, result_false$se, tolerance = 1e-15)
  expect_equal(result_true$df, result_false$df)
  expect_equal(result_true$controls_tier, result_false$controls_tier)
})

# =============================================================================
# Group G: Edge cases and output format (T-30 to T-33, T-38 to T-40) — Task E2-05.13
# =============================================================================

test_that("T-30 single post period returns single-row data.frame", {
  dt <- make_period_panel(N = 10L, TT = 4L, S = 4L, n_treated = 5L, seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                     x = NULL, periods = 4L)
  expect_equal(nrow(result), 1L)
  expect_s3_class(result, "data.frame")
  expect_false(is.na(result$att))
})

test_that("T-31 minimal sample (1 treated + 1 control) triggers N<3 error", {
  dt <- data.table::data.table(
    id = c(1L, 2L), time = rep(4L, 2),
    y_trans = c(10, 3),
    d = c(1L, 0L)
  )
  result <- quiet_lwdid(
    estimate_period_effects(dt, "y_trans", "d", "time",
                            x = NULL, periods = 4L)
  )
  # N=2 < 3 -> error caught by tryCatch -> NA row
  expect_equal(result$n_obs, 2L)
  expect_equal(result$n_treated, 1L)
  expect_equal(result$n_control, 1L)
  expect_true(is.na(result$att))
})

test_that("T-32 many post periods (>20) returns correct number of rows", {
  N <- 30L; TT <- 25L; S <- 5L; n_treated <- 15L
  alpha_i <- seq(10, by = 2, length.out = N)
  delta_t <- seq(0, by = 1, length.out = TT)
  set.seed(99L)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, y := alpha_i[id] + delta_t[time] + 2 * d * as.integer(time >= S) +
       rnorm(.N, sd = 0.3)]
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = S))
  result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                     x = NULL, periods = S:TT)
  expect_equal(nrow(result), TT - S + 1L)  # 21 post periods
  expect_true(all(!is.na(result$att)))
})

test_that("T-33 all periods fail returns all-NA data.frame", {
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  dt_trans <- dt_trans[!(d == 1L & time >= 4L)]
  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = 4:6)
  )
  expect_equal(nrow(result), 3L)
  expect_true(all(is.na(result$att)))
  expect_true(all(is.na(result$se)))
})

test_that("T-38 output rows = length(periods)", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  expect_equal(nrow(result), length(4:6))
})

test_that("T-39 column types are correct", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  expect_true(is.integer(result$tindex))
  expect_true(is.integer(result$period))
  expect_true(is.integer(result$n_obs))
  expect_true(is.integer(result$n_treated))
  expect_true(is.integer(result$n_control))
  expect_true(is.integer(result$df))
  expect_true(is.numeric(result$att))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$t_stat))
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(is.character(result$vce_type))
  expect_true(is.integer(result$n_clusters))
  expect_true(is.character(result$controls_tier))
})

test_that("T-40 no anomalous rownames", {
  dt <- make_period_panel(seed = 42L)
  result <- quiet_lwdid(run_period_effects(dt, S = 4L, TT = 6L))
  rn <- rownames(result)
  expect_true(all(rn %in% as.character(seq_len(nrow(result)))))
})

# =============================================================================
# Group H: Numerical precision and Stata consistency (T-34 to T-37) — Task E2-05.14
# =============================================================================

test_that("T-34 California Smoking data period ATTs (Stata deferred)", {
  skip("Stata lwdid.ado period-level results not yet available")
})

test_that("T-35 California Smoking data period SEs (Stata deferred)", {
  skip("Stata lwdid.ado period-level results not yet available")
})

test_that("T-36 noiseless DGP: period ATTs trigger degenerate SE path", {
  # Noiseless DGP produces perfect fit (residuals = 0), SE = 0 < 1e-10
  # estimate_ra_common detects degenerate SE, period-level returns NA
  # This is correct behavior: with zero residual variance, inference is undefined
  N <- 10L; TT <- 6L; S <- 4L; n_treated <- 5L
  alpha_i <- seq(10, by = 5, length.out = N)
  delta_t <- seq(0, by = 2, length.out = TT)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, post := as.integer(time >= S)]
  dt[, tau_t := ifelse(post == 1L, time - S + 1L, 0L)]
  dt[, y := alpha_i[id] + delta_t[time] + tau_t * d]
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = S))
  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = S:TT)
  )
  # Degenerate SE -> ATT set to NA at period level
  expect_true(all(is.na(result$att)))
  expect_true(all(is.na(result$se)))
  # But sample sizes are preserved
  expect_true(all(result$n_obs == 10L))
  expect_true(all(result$n_treated == 5L))
  expect_true(all(result$n_control == 5L))

  # Verify the underlying ATT is correct via estimate_ra_common directly
  # (it returns att even when SE is degenerate)
  sub4 <- dt_trans[time == 4L]
  ra_res <- suppressWarnings(
    estimate_ra_common(sub4$y_trans, sub4$d)
  )
  expect_equal(ra_res$att, 1.0, tolerance = 1e-10)
  expect_true(is.na(ra_res$se))  # degenerate SE
})

test_that("T-37 noiseless DGP: zero effect triggers degenerate SE path", {
  # Zero treatment effect with noiseless DGP -> perfect fit -> SE = 0
  # Period-level returns NA for ATT (degenerate SE path)
  N <- 10L; TT <- 6L; S <- 4L; n_treated <- 5L
  alpha_i <- seq(10, by = 5, length.out = N)
  delta_t <- seq(0, by = 2, length.out = TT)
  dt <- data.table::CJ(id = seq_len(N), time = seq_len(TT))
  dt[, d := as.integer(id <= n_treated)]
  dt[, y := alpha_i[id] + delta_t[time]]
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = S))
  result <- quiet_lwdid(
    estimate_period_effects(dt_trans, "y_trans", "d", "time",
                            x = NULL, periods = S:TT)
  )
  # Degenerate SE -> ATT set to NA at period level
  expect_true(all(is.na(result$att)))
  expect_true(all(is.na(result$se)))

  # Verify the underlying ATT is correct via estimate_ra_common directly
  sub4 <- dt_trans[time == 4L]
  ra_res <- suppressWarnings(
    estimate_ra_common(sub4$y_trans, sub4$d)
  )
  expect_equal(ra_res$att, 0.0, tolerance = 1e-10)
  expect_true(is.na(ra_res$se))  # degenerate SE
})

# =============================================================================
# Group J: vibe-math MCP independent verification (T-41 to T-44) — Task E2-05.15
# All expected values independently computed via vibe-math MCP tools.
# =============================================================================

test_that("T-41 vibe-math verified: period ATT = mean(treated) - mean(control)", {
  # vibe-math: ATT = (3+5)/2 - (1+2)/2 = 4.0 - 1.5 = 2.5
  dt <- data.table::data.table(
    id = 1:4, time = rep(2L, 4),
    y_trans = c(3, 5, 1, 2),
    d = c(1L, 1L, 0L, 0L)
  )
  result <- estimate_period_effects(dt, "y_trans", "d", "time",
                                     x = NULL, periods = 2L)
  expect_equal(result$att, 2.5, tolerance = 1e-12)
})

test_that("T-42 vibe-math verified: period SE via OLS formula", {
  # vibe-math: (X'X)^{-1} = [[0.5,-0.5],[-0.5,1.0]]
  # SSR = 2.5, sigma2 = 2.5/2 = 1.25, SE = sqrt(1.25 * 1.0) = 1.118033988749895
  dt <- data.table::data.table(
    id = 1:4, time = rep(2L, 4),
    y_trans = c(3, 5, 1, 2),
    d = c(1L, 1L, 0L, 0L)
  )
  result <- estimate_period_effects(dt, "y_trans", "d", "time",
                                     x = NULL, periods = 2L)
  expect_equal(result$se, 1.118033988749895, tolerance = 1e-12)
  expect_equal(result$df, 2L)
})

test_that("T-43 vibe-math verified: period t_stat and pvalue", {
  # vibe-math: t = 2.5 / 1.118033988749895 = 2.23606797749979
  # R: pvalue = 2*pt(-2.23606797749979, 2) = 0.154845745271483
  dt <- data.table::data.table(
    id = 1:4, time = rep(2L, 4),
    y_trans = c(3, 5, 1, 2),
    d = c(1L, 1L, 0L, 0L)
  )
  result <- estimate_period_effects(dt, "y_trans", "d", "time",
                                     x = NULL, periods = 2L)
  expect_equal(result$t_stat, 2.23606797749979, tolerance = 1e-10)
  expect_equal(result$pvalue, 0.154845745271483, tolerance = 1e-10)
  tcrit <- qt(0.975, 2)
  expect_equal(result$ci_lower, 2.5 - tcrit * 1.118033988749895, tolerance = 1e-10)
  expect_equal(result$ci_upper, 2.5 + tcrit * 1.118033988749895, tolerance = 1e-10)
})

test_that("T-44 vibe-math verified: summary ATT equivalence (balanced, no controls)", {
  # vibe-math: mean of period ATTs = 1.9115989304272543
  dt <- make_period_panel(seed = 42L)
  dt_trans <- suppressWarnings(transform_demean(dt, "y", "id", "time", g = 4L))
  period_result <- estimate_period_effects(dt_trans, "y_trans", "d", "time",
                                            x = NULL, periods = 4:6)
  mean_period_att <- mean(period_result$att)
  post_dt <- dt_trans[time >= 4L]
  y_bar <- post_dt[, .(y_bar = mean(y_trans)), by = id]
  d_map <- unique(dt[, .(id, d)])
  y_bar <- merge(y_bar, d_map, by = "id")
  summary_att <- mean(y_bar[d == 1]$y_bar) - mean(y_bar[d == 0]$y_bar)
  expect_equal(mean_period_att, 1.9115989304272543, tolerance = 1e-10)
  expect_equal(summary_att, 1.9115989304272543, tolerance = 1e-10)
  expect_equal(mean_period_att, summary_att, tolerance = 1e-12)
})
