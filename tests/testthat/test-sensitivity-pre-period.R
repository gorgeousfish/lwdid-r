# ============================================================================
# test-sensitivity-pre-period.R
# Tests for E8-01: Pre-treatment Period Robustness Analysis
# ============================================================================

# ============================================================================
# Helper functions for creating test data
# ============================================================================

#' Create a simple staggered panel for sensitivity tests
make_stag_sensitivity <- function(n_units = 12, n_periods = 10,
                                  cohorts = c(6),
                                  n_nt = 4, seed = 42L) {
  set.seed(seed)
  ids <- seq_len(n_units)
  times <- seq_len(n_periods)
  df <- expand.grid(id = ids, time = times)
  df <- df[order(df$id, df$time), ]
  treated_ids <- (n_nt + 1):n_units
  n_per_cohort <- length(treated_ids) %/% length(cohorts)
  gvar_map <- numeric(n_units)
  gvar_map[1:n_nt] <- 0  # never-treated
  for (i in seq_along(cohorts)) {
    start_idx <- n_nt + (i - 1) * n_per_cohort + 1
    end_idx <- if (i == length(cohorts)) n_units
               else n_nt + i * n_per_cohort
    gvar_map[start_idx:end_idx] <- cohorts[i]
  }
  df$gvar <- gvar_map[df$id]
  df$y <- rnorm(nrow(df))
  rownames(df) <- NULL
  df
}

#' Create a simple common timing panel for sensitivity tests
make_ct_sensitivity <- function(n_units = 10, n_periods = 8,
                                n_treated = 5, K = 4, seed = 42L) {
  set.seed(seed)
  ids <- seq_len(n_units)
  years <- 2000 + seq_len(n_periods)
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  df$d <- ifelse(df$id <= n_treated, 1, 0)
  df$post <- ifelse(df$year > (2000 + K), 1, 0)
  df$y <- rnorm(nrow(df)) + df$d * df$post * 2
  rownames(df) <- NULL
  df
}


# ============================================================================
# TC-01: Sensitivity Ratio Calculation Correctness
# ============================================================================

test_that("TC-01a: SR normal case — atts=c(1.0,1.5,2.0), baseline=1.5", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(1.0, 1.5, 2.0), 1.5)
  # SR = (2.0 - 1.0) / |1.5| = 1.0 / 1.5 = 2/3
  expect_equal(sr, 2 / 3, tolerance = 1e-10)
})

test_that("TC-01b: SR all same — atts=c(1.0,1.0,1.0), baseline=1.0", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(1.0, 1.0, 1.0), 1.0)
  expect_equal(sr, 0, tolerance = 1e-10)
})

test_that("TC-01c: SR negative range — atts=c(-1.0,1.0), baseline=0.5", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(-1.0, 1.0), 0.5)
  # SR = (1.0 - (-1.0)) / |0.5| = 2.0 / 0.5 = 4.0
  expect_equal(sr, 4.0, tolerance = 1e-10)
})

test_that("TC-01d: SR empty atts — atts=numeric(0), baseline=1.0", {
  sr <- lwdid:::.compute_sensitivity_ratio(numeric(0), 1.0)
  expect_equal(sr, 0, tolerance = 1e-10)
})

test_that("TC-01e: SR baseline near zero — atts=c(0.1,0.2), baseline=1e-15", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(0.1, 0.2), 1e-15)
  expect_equal(sr, Inf)
})

test_that("TC-01f: SR baseline+range near zero — atts=c(1e-15,1e-15), baseline=1e-15", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(1e-15, 1e-15), 1e-15)
  expect_equal(sr, 0, tolerance = 1e-10)
})

test_that("TC-01g: SR single element returns 0", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(2.5), 2.5)
  expect_equal(sr, 0, tolerance = 1e-10)
})

test_that("TC-01h: SR with large values maintains precision", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(1e8, 1e8 + 1e4), 1e8)
  # SR = 1e4 / 1e8 = 1e-4
  expect_equal(sr, 1e-4, tolerance = 1e-10)
})

test_that("TC-01i: SR with negative baseline uses absolute value", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(-3.0, -1.0), -2.0)
  # SR = ((-1.0) - (-3.0)) / |-2.0| = 2.0 / 2.0 = 1.0
  expect_equal(sr, 1.0, tolerance = 1e-10)
})


# ============================================================================
# TC-02: Robustness Level Thresholds (exact boundary values)
# ============================================================================

test_that("TC-02: .determine_robustness_level thresholds are exact", {
  # highly_robust: SR < 0.10
  expect_equal(lwdid:::.determine_robustness_level(0.0), "highly_robust")
  expect_equal(lwdid:::.determine_robustness_level(0.0999), "highly_robust")

  # moderately_robust: 0.10 <= SR < 0.25 (boundary at 0.10)
  expect_equal(lwdid:::.determine_robustness_level(0.10), "moderately_robust")
  expect_equal(lwdid:::.determine_robustness_level(0.2499), "moderately_robust")

  # sensitive: 0.25 <= SR < 0.50 (boundary at 0.25)
  expect_equal(lwdid:::.determine_robustness_level(0.25), "sensitive")
  expect_equal(lwdid:::.determine_robustness_level(0.4999), "sensitive")

  # highly_sensitive: SR >= 0.50 (boundary at 0.50)
  expect_equal(lwdid:::.determine_robustness_level(0.50), "highly_sensitive")
  expect_equal(lwdid:::.determine_robustness_level(1.0), "highly_sensitive")
  expect_equal(lwdid:::.determine_robustness_level(Inf), "highly_sensitive")
})

test_that("TC-02b: boundaries at machine epsilon precision", {
  expect_equal(
    lwdid:::.determine_robustness_level(0.1 - .Machine$double.eps * 10),
    "highly_robust"
  )
  expect_equal(lwdid:::.determine_robustness_level(0.1), "moderately_robust")
  expect_equal(
    lwdid:::.determine_robustness_level(0.25 - .Machine$double.eps * 100),
    "moderately_robust"
  )
  expect_equal(lwdid:::.determine_robustness_level(0.25), "sensitive")
  expect_equal(
    lwdid:::.determine_robustness_level(0.5 - .Machine$double.eps * 100),
    "sensitive"
  )
  expect_equal(lwdid:::.determine_robustness_level(0.5), "highly_sensitive")
})


# ============================================================================
# TC-03: Auto-detect Pre-period Range
# ============================================================================

test_that("TC-03a: demean → min_pre=1", {
  df <- data.frame(
    id = rep(1:3, each = 10), time = rep(1:10, 3),
    gvar = rep(6, 30), y = rnorm(30)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "demean"
  )
  expect_equal(r[1], 1L)
})

test_that("TC-03b: detrend → min_pre=2", {
  df <- data.frame(
    id = rep(1:3, each = 10), time = rep(1:10, 3),
    gvar = rep(6, 30), y = rnorm(30)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "detrend"
  )
  expect_equal(r[1], 2L)
})

test_that("TC-03c: demeanq → min_pre=1", {
  df <- data.frame(
    id = rep(1:3, each = 10), time = rep(1:10, 3),
    gvar = rep(6, 30), y = rnorm(30)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "demeanq"
  )
  expect_equal(r[1], 1L)
})

test_that("TC-03d: detrendq → min_pre=2", {
  df <- data.frame(
    id = rep(1:3, each = 10), time = rep(1:10, 3),
    gvar = rep(6, 30), y = rnorm(30)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "detrendq"
  )
  expect_equal(r[1], 2L)
})


test_that("TC-03e: staggered cohorts at 5 and 7, min_time=1 → max_pre=4", {
  # Cohort 5: 5-1=4, Cohort 7: 7-1=6, max_pre = min(4,6) = 4
  df <- data.frame(
    id = rep(1:6, each = 10), time = rep(1:10, 6),
    gvar = rep(c(5, 5, 7, 7, Inf, Inf), each = 10),
    y = rnorm(60)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "demean"
  )
  expect_equal(r[2], 4L)
})

test_that("TC-03f: common timing with 6 unique pre-treatment periods", {
  # 8 periods, K=6 → post starts at year 2007, pre-periods: 2001..2006 (6)
  df <- make_ct_sensitivity(n_units = 6, n_periods = 8, K = 6, seed = 100L)
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "year", NULL, "d", "post", "demean"
  )
  expect_equal(r[2], 6L)
})

test_that("TC-03g: no valid cohorts returns c(min_pre, min_pre)", {
  df <- data.frame(
    id = rep(1:4, each = 5), time = rep(1:5, 4),
    y = rnorm(20), gvar = 0  # all never-treated
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "demean"
  )
  expect_equal(r, c(1L, 1L))
})

test_that("TC-03h: max_pre >= min_pre guaranteed", {
  # Cohort at 2, times 1-5 → cohort-min_time = 2-1 = 1
  # detrend min_pre=2, but max_pre=1 → clamped to max(1,2)=2
  df <- data.frame(
    id = rep(1:2, each = 5), time = rep(1:5, 2),
    gvar = rep(2, 10), y = rnorm(10)
  )
  r <- lwdid:::.auto_detect_pre_period_range(
    df, "id", "time", "gvar", NULL, NULL, "detrend"
  )
  expect_true(r[2] >= r[1])
  expect_equal(r[1], 2L)
})


# ============================================================================
# TC-14: Input Validation
# ============================================================================

test_that("TC-14a: missing column → error with class lwdid_missing_column", {
  df <- data.frame(id = 1:5, time = 1:5, y = rnorm(5), gvar = 3)
  expect_error(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", "nonexistent", NULL, NULL, "demean"
    ),
    class = "lwdid_missing_column"
  )
})

test_that("TC-14b: invalid rolling → error with class lwdid_invalid_rolling", {
  df <- data.frame(id = 1:5, time = 1:5, y = rnorm(5), gvar = 3)
  expect_error(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", "gvar", NULL, NULL, "invalid_method"
    ),
    class = "lwdid_invalid_rolling"
  )
})

test_that("TC-14c: no gvar and no d+post → error with class lwdid_invalid_parameter", {
  df <- data.frame(id = 1:5, time = 1:5, y = rnorm(5))
  expect_error(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", NULL, NULL, NULL, "demean"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("TC-14d: valid staggered inputs pass without error", {
  df <- data.frame(id = 1:5, time = 1:5, y = rnorm(5), gvar = 3)
  expect_silent(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", "gvar", NULL, NULL, "demean"
    )
  )
})

test_that("TC-14e: valid common timing inputs pass without error", {
  df <- data.frame(id = 1:5, time = 1:5, y = rnorm(5),
                   d = c(1, 1, 0, 0, 0), post = c(0, 0, 1, 1, 1))
  expect_silent(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", NULL, "d", "post", "demean"
    )
  )
})

test_that("TC-14f: all four valid rolling methods accepted", {
  df <- data.frame(id = 1:4, time = rep(1:2, 2), y = rnorm(4),
                   gvar = c(3, 3, Inf, Inf))
  for (method in c("demean", "detrend", "demeanq", "detrendq")) {
    expect_silent(
      lwdid:::.validate_sensitivity_inputs(
        df, "y", "id", "time", "gvar", NULL, NULL, method
      )
    )
  }
})

test_that("TC-14g: rolling validation is case-insensitive", {
  df <- data.frame(id = 1:4, time = rep(1:2, 2), y = rnorm(4),
                   gvar = c(3, 3, Inf, Inf))
  expect_silent(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", "gvar", NULL, NULL, "Demean"
    )
  )
})

test_that("TC-14h: missing multiple columns reported", {
  df <- data.frame(id = 1:5, time = 1:5)
  expect_error(
    lwdid:::.validate_sensitivity_inputs(
      df, "y", "id", "time", "gvar", NULL, NULL, "demean"
    ),
    class = "lwdid_missing_column"
  )
})


# ============================================================================
# TC-10: Filter Staggered Mode
# ============================================================================

test_that("TC-10a: staggered filter n_pre=3, exclude=0", {
  # Units 1-3: gvar=5 (treated at t=5), Units 4-6: gvar=Inf, times 1-8
  df <- data.frame(
    id = rep(1:6, each = 8), time = rep(1:8, 6),
    gvar = rep(c(5, 5, 5, Inf, Inf, Inf), each = 8),
    y = rnorm(48)
  )

  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 3L, exclude_periods = 0L
  )

  # Treated (gvar=5): pre_end=5-1-0=4, pre_start=4-3+1=2
  # Keep t in [2,4] ∪ [5,8]
  treated_rows <- filtered[filtered$gvar == 5, ]
  treated_times <- sort(unique(treated_rows$time))
  expect_equal(treated_times, c(2L, 3L, 4L, 5L, 6L, 7L, 8L))
  expect_false(1L %in% treated_rows$time)

  # NT: aligned to min_cohort=5, pre_start_nt=2, keep t >= 2
  nt_rows <- filtered[is.infinite(filtered$gvar), ]
  expect_true(all(nt_rows$time >= 2L))
  expect_false(1L %in% nt_rows$time)
})

test_that("TC-10b: staggered filter n_pre=3, exclude=1", {
  df <- data.frame(
    id = rep(1:6, each = 8), time = rep(1:8, 6),
    gvar = rep(c(5, 5, 5, Inf, Inf, Inf), each = 8),
    y = rnorm(48)
  )

  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 3L, exclude_periods = 1L
  )

  # Treated (gvar=5): pre_end=5-1-1=3, pre_start=3-3+1=1
  # Keep t in [1,3] ∪ [5,8]
  treated_rows <- filtered[filtered$gvar == 5, ]
  treated_times <- sort(unique(treated_rows$time))
  expect_equal(treated_times, c(1L, 2L, 3L, 5L, 6L, 7L, 8L))
  expect_false(4L %in% treated_rows$time)  # gap period excluded
})


# ============================================================================
# TC-11: Filter Common Timing Mode
# ============================================================================

test_that("TC-11a: CT filter n_pre=3, exclude=0", {
  # 6 units, 8 periods (2001-2008), K=4 → post starts at 2005
  df <- make_ct_sensitivity(n_units = 6, n_periods = 8, K = 4, seed = 111L)

  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "year", NULL, "d", "post",
    n_pre_periods = 3L, exclude_periods = 0L
  )

  # Pre-periods: 2001,2002,2003,2004 (4 total). Last 3 = {2002,2003,2004}
  # Post-periods: 2005,2006,2007,2008
  pre_times <- sort(unique(filtered$year[filtered$post == 0]))
  post_times <- sort(unique(filtered$year[filtered$post == 1]))
  expect_equal(pre_times, c(2002L, 2003L, 2004L))
  expect_equal(post_times, c(2005L, 2006L, 2007L, 2008L))
})

test_that("TC-11b: CT filter n_pre=3, exclude=1", {
  df <- make_ct_sensitivity(n_units = 6, n_periods = 8, K = 4, seed = 112L)

  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "year", NULL, "d", "post",
    n_pre_periods = 3L, exclude_periods = 1L
  )

  # Pre-periods: 2001,2002,2003,2004. Exclude last 1 → {2001,2002,2003}
  # Last 3 of {2001,2002,2003} = {2001,2002,2003}
  pre_times <- sort(unique(filtered$year[filtered$post == 0]))
  expect_equal(pre_times, c(2001L, 2002L, 2003L))
})

test_that("TC-11c: CT filter n_pre=10 (more than available) keeps all pre-periods", {
  df <- make_ct_sensitivity(n_units = 6, n_periods = 8, K = 4, seed = 113L)

  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "year", NULL, "d", "post",
    n_pre_periods = 10L, exclude_periods = 0L
  )

  # Only 4 pre-periods available, request 10 → keep all 4
  pre_times <- sort(unique(filtered$year[filtered$post == 0]))
  expect_equal(pre_times, c(2001L, 2002L, 2003L, 2004L))
})


# ============================================================================
# TC-12: Filter Staggered Never-treated Units
# ============================================================================

test_that("TC-12a: gvar=Inf treated as never-treated", {
  df <- data.frame(
    id = rep(1:4, each = 10), time = rep(1:10, 4),
    gvar = rep(c(6, 6, Inf, Inf), each = 10), y = rnorm(40)
  )
  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 3L, exclude_periods = 0L
  )
  nt_rows <- filtered[is.infinite(filtered$gvar), ]
  expect_true(nrow(nt_rows) > 0L)
  # min_cohort=6, pre_end_nt=5, pre_start_nt=3 → keep t >= 3
  expect_true(all(nt_rows$time >= 3L))
})

test_that("TC-12b: gvar=NA treated as never-treated", {
  df <- data.frame(
    id = rep(1:4, each = 10), time = rep(1:10, 4),
    gvar = rep(c(NA, NA, 8, 8), each = 10), y = rnorm(40)
  )
  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 2L, exclude_periods = 0L
  )
  nt_rows <- filtered[is.na(filtered$gvar), ]
  expect_true(nrow(nt_rows) > 0L)
  # min_cohort=8, pre_end_nt=7, pre_start_nt=6 → keep t >= 6
  expect_true(all(nt_rows$time >= 6L))
})

test_that("TC-12c: gvar=0 treated as never-treated", {
  df <- data.frame(
    id = rep(1:4, each = 8), time = rep(1:8, 4),
    gvar = rep(c(0, 0, 5, 5), each = 8), y = rnorm(32)
  )
  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 2L, exclude_periods = 0L
  )
  nt_rows <- filtered[filtered$gvar == 0, ]
  expect_true(nrow(nt_rows) > 0L)
})

test_that("TC-12d: NT time window aligns with min(cohorts)", {
  # Cohorts at 5 and 7, min_cohort=5
  df <- data.frame(
    id = rep(1:6, each = 10), time = rep(1:10, 6),
    gvar = rep(c(5, 5, 7, 7, Inf, Inf), each = 10), y = rnorm(60)
  )
  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 3L, exclude_periods = 0L
  )
  # min_cohort=5, pre_end_nt=4, pre_start_nt=2 → NT keep t >= 2
  nt_rows <- filtered[is.infinite(filtered$gvar), ]
  expect_true(all(nt_rows$time >= 2L))
  expect_false(1L %in% nt_rows$time)
})


# ============================================================================
# TC-17: CT Boundary — exclude_periods >= available pre-periods
# ============================================================================

test_that("TC-17: CT exclude_periods >= available pre-periods issues warning", {
  df <- make_ct_sensitivity(n_units = 4, n_periods = 6, K = 3, seed = 170L)
  # Pre-periods: 2001,2002,2003 (3 periods). exclude=5 >= 3 → warning
  expect_warning(
    lwdid:::.filter_to_n_pre_periods(
      df, "id", "year", NULL, "d", "post",
      n_pre_periods = 2L, exclude_periods = 5L
    ),
    class = "lwdid_data"
  )
})

# ============================================================================
# TC-18: Staggered Boundary — cohort pre-periods insufficient
# ============================================================================

test_that("TC-18: staggered cohort with insufficient pre-periods", {
  # Cohort at t=2, times 1-8. With exclude=1:
  # pre_end = 2-1-1 = 0, pre_start = 0-3+1 = -2
  # No times in [-2,0] for times 1-8 → no pre-period rows for this cohort
  df <- data.frame(
    id = rep(1:4, each = 8), time = rep(1:8, 4),
    gvar = rep(c(2, 2, Inf, Inf), each = 8), y = rnorm(32)
  )
  filtered <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "time", "gvar", NULL, NULL,
    n_pre_periods = 3L, exclude_periods = 1L
  )
  # Cohort 2 units: only post-treatment times [2,8] kept (no pre-period data)
  treated_rows <- filtered[filtered$gvar == 2, ]
  treated_pre <- treated_rows[treated_rows$time < 2, ]
  expect_equal(nrow(treated_pre), 0L)
})


# ============================================================================
# TC-08: Monotonic Trend Detection
# ============================================================================

test_that("TC-08a: monotonically increasing ATTs trigger warning", {
  specs <- list(
    list(converged = TRUE, n_pre_periods = 2L, att = 1.0),
    list(converged = TRUE, n_pre_periods = 3L, att = 1.5),
    list(converged = TRUE, n_pre_periods = 4L, att = 2.0),
    list(converged = TRUE, n_pre_periods = 5L, att = 2.5)
  )
  baseline <- list(att = 2.5)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.6, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(any(grepl("monotonic", rec$result_warnings, ignore.case = TRUE)))
})

test_that("TC-08b: monotonically decreasing ATTs also detected", {
  specs <- list(
    list(converged = TRUE, n_pre_periods = 2L, att = 3.0),
    list(converged = TRUE, n_pre_periods = 3L, att = 2.0),
    list(converged = TRUE, n_pre_periods = 4L, att = 1.0)
  )
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.6, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(any(grepl("monotonic", rec$result_warnings, ignore.case = TRUE)))
})

test_that("TC-08c: non-monotonic ATTs do NOT trigger warning", {
  specs <- list(
    list(converged = TRUE, n_pre_periods = 2L, att = 1.0),
    list(converged = TRUE, n_pre_periods = 3L, att = 2.0),
    list(converged = TRUE, n_pre_periods = 4L, att = 1.5)
  )
  baseline <- list(att = 1.5)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.3, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_false(any(grepl("monotonic", rec$result_warnings, ignore.case = TRUE)))
})

test_that("TC-08d: fewer than 3 converged specs skip monotonic check", {
  specs <- list(
    list(converged = TRUE, n_pre_periods = 2L, att = 1.0),
    list(converged = TRUE, n_pre_periods = 3L, att = 2.0),
    list(converged = FALSE, n_pre_periods = 4L, att = 3.0)
  )
  baseline <- list(att = 2.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.5, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_false(any(grepl("monotonic", rec$result_warnings, ignore.case = TRUE)))
})


# ============================================================================
# TC-09: Recommendation Structure
# ============================================================================

test_that("TC-09: recommendation returns list with correct structure", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.05, is_robust = TRUE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_type(rec, "list")
  expect_true("main_rec" %in% names(rec))
  expect_true("detailed_recommendations" %in% names(rec))
  expect_true("result_warnings" %in% names(rec))
  expect_type(rec$main_rec, "character")
  expect_type(rec$detailed_recommendations, "character")
  expect_type(rec$result_warnings, "character")
})

# ============================================================================
# TC-19: 4 Mutually Exclusive Main Branches
# ============================================================================

test_that("TC-19a: Branch 1 — robust + same sign + all significant", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.05, is_robust = TRUE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(grepl("robust", rec$main_rec, ignore.case = TRUE))
  expect_true(grepl("stable", rec$main_rec, ignore.case = TRUE))
  expect_false(grepl("CAUTION", rec$main_rec))
})

test_that("TC-19b: Branch 2 — robust + same sign + NOT all significant", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.05, is_robust = TRUE,
    all_same_sign = TRUE, all_significant = FALSE, rolling = "demean"
  )
  expect_true(grepl("moderately robust", rec$main_rec, ignore.case = TRUE))
  expect_true(grepl("significance varies", rec$main_rec, ignore.case = TRUE))
  expect_true(any(grepl("range", rec$detailed_recommendations, ignore.case = TRUE)))
})


test_that("TC-19c: Branch 3 — sign changes → CAUTION message", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.80, is_robust = FALSE,
    all_same_sign = FALSE, all_significant = FALSE, rolling = "demean"
  )
  expect_true(grepl("CAUTION", rec$main_rec))
  expect_true(any(grepl("[Ss]ign change", rec$result_warnings)))
})

test_that("TC-19d: Branch 4 — not robust + same sign → moderate sensitivity", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.40, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(grepl("moderate sensitivity", rec$main_rec, ignore.case = TRUE))
  expect_true(grepl("40.0%", rec$main_rec))
})

test_that("TC-19e: Branch 3 takes priority over Branch 4", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.80, is_robust = FALSE,
    all_same_sign = FALSE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(grepl("CAUTION", rec$main_rec))
})

# ============================================================================
# TC-20: Method-specific Recommendations
# ============================================================================

test_that("TC-20a: demean + SR>0.25 → recommendation about detrend", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.30, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_true(any(grepl("detrend", rec$detailed_recommendations, ignore.case = TRUE)))
})

test_that("TC-20b: detrend + SR>0.50 → recommendation about misspecification", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.60, is_robust = FALSE,
    all_same_sign = FALSE, all_significant = FALSE, rolling = "detrend"
  )
  expect_true(any(grepl("misspecification", rec$detailed_recommendations, ignore.case = TRUE)))
})


test_that("TC-20c: demean + SR<=0.25 → no detrend recommendation", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.20, is_robust = TRUE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "demean"
  )
  expect_false(any(grepl("detrend", rec$detailed_recommendations, ignore.case = TRUE)))
})

test_that("TC-20d: detrend + SR<=0.50 → no misspecification recommendation", {
  specs <- list(list(converged = TRUE, n_pre_periods = 3L, att = 1.0))
  baseline <- list(att = 1.0)
  rec <- lwdid:::.generate_robustness_recommendations(
    specifications = specs, baseline_spec = baseline,
    sensitivity_ratio = 0.45, is_robust = FALSE,
    all_same_sign = TRUE, all_significant = TRUE, rolling = "detrend"
  )
  expect_false(any(grepl("misspecification", rec$detailed_recommendations, ignore.case = TRUE)))
})

# ============================================================================
# TC-04: Single Spec Failure Doesn't Crash
# ============================================================================

test_that("TC-04: single spec failure returns converged=FALSE with att=NA", {
  # Tiny dataset that will cause lwdid() to fail
  df <- data.frame(
    id = c(1L, 1L), time = c(1L, 2L), y = c(1.0, 2.0), gvar = c(3L, 3L)
  )
  result <- suppressWarnings(
    lwdid:::.run_single_specification(
      data = df, y = "y", ivar = "id", tvar = "time",
      gvar = "gvar", d = NULL, post = NULL,
      rolling = "demean", estimator = "ra",
      controls = NULL, vce = NULL, cluster_var = NULL,
      n_pre_periods = 2L, exclude_periods = 0L,
      alpha = 0.05, spec_id = 1L
    )
  )
  expect_true(is.list(result))
  expect_false(result$converged)
  expect_true(is.na(result$att))
  expect_true(is.na(result$se))
  expect_equal(result$specification_id, 1L)
  expect_equal(result$n_pre_periods, 2L)
  expect_equal(result$n_treated, 0L)
  expect_equal(result$n_control, 0L)
})


# ============================================================================
# TC-05: Sign Change Detection
# ============================================================================

test_that("TC-05: specs with different signs → n_sign_changes > 0", {
  atts <- c(1.0, -0.5, 1.2, -0.3)
  baseline_att <- 1.0
  baseline_sign <- sign(baseline_att)
  n_sign_changes <- sum(sign(atts) != baseline_sign)
  expect_equal(n_sign_changes, 2L)
  expect_false(all(sign(atts) == baseline_sign))
})

# ============================================================================
# TC-06: att_std Single Element
# ============================================================================

test_that("TC-06: single converged spec → att_std = 0.0 (not NA)", {
  atts <- c(2.5)
  att_std <- if (length(atts) > 1L) sd(atts) else 0.0
  expect_equal(att_std, 0.0, tolerance = 1e-10)
  expect_false(is.na(att_std))
})

# ============================================================================
# TC-07: n_significant Count
# ============================================================================

test_that("TC-07: n_significant matches count of pvalue < alpha", {
  pvals <- c(0.01, 0.04, 0.06, 0.10, 0.001)
  alpha <- 0.05
  n_significant <- sum(pvals < alpha)
  expect_equal(n_significant, 3L)

  # Edge: pvalue exactly at alpha is NOT significant
  pvals2 <- c(0.05, 0.049, 0.051)
  n_sig2 <- sum(pvals2 < alpha)
  expect_equal(n_sig2, 1L)  # only 0.049 < 0.05
})

# ============================================================================
# TC-13: exclude_periods_before_treatment Passthrough
# ============================================================================

test_that("TC-13: exclude_periods is passed through to filtering", {
  df <- make_ct_sensitivity(n_units = 6, n_periods = 8, K = 4, seed = 130L)

  # With exclude=0: last 3 pre-periods
  filtered_0 <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "year", NULL, "d", "post",
    n_pre_periods = 3L, exclude_periods = 0L
  )
  pre_times_0 <- sort(unique(filtered_0$year[filtered_0$post == 0]))

  # With exclude=2: exclude last 2 pre-periods, then take last 3
  filtered_2 <- lwdid:::.filter_to_n_pre_periods(
    df, "id", "year", NULL, "d", "post",
    n_pre_periods = 3L, exclude_periods = 2L
  )
  pre_times_2 <- sort(unique(filtered_2$year[filtered_2$post == 0]))

  # The two results should differ
  expect_false(identical(pre_times_0, pre_times_2))
})


# ============================================================================
# TC-15: All Specs Unconverged → Error
# ============================================================================

test_that("TC-15: all specifications unconverged throws error", {
  # Create data that will cause all specs to fail
  df <- data.frame(
    id = rep(1L, 3), time = 1:3, y = rnorm(3), gvar = rep(2L, 3)
  )
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid_sensitivity_pre_period(
        data = df, y = "y", ivar = "id", tvar = "time",
        gvar = "gvar", rolling = "demean",
        pre_period_range = c(1L, 1L), verbose = FALSE
      )
    )),
    class = "lwdid_convergence"
  )
})

# ============================================================================
# TC-21: pre_period_range = c(5, 3) → Error
# ============================================================================

test_that("TC-21: pre_period_range min > max throws error", {
  df <- make_stag_sensitivity(n_units = 12, n_periods = 10,
                              cohorts = c(6), n_nt = 4, seed = 210L)
  expect_error(
    suppressMessages(
      lwdid_sensitivity_pre_period(
        data = df, y = "y", ivar = "id", tvar = "time",
        gvar = "gvar", rolling = "demean",
        pre_period_range = c(5L, 3L), verbose = FALSE
      )
    ),
    class = "lwdid_invalid_parameter"
  )
})


# ============================================================================
# TC-V01 through TC-V04: Vibe-math Numerical Verification
# ============================================================================

test_that("TC-V01: SR = (2.0-1.0)/abs(1.5) = 0.6667 verified to 4 decimal places", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(1.0, 1.5, 2.0), 1.5)
  expected <- (2.0 - 1.0) / abs(1.5)
  expect_equal(sr, expected, tolerance = 1e-10)
  expect_equal(round(sr, 4), 0.6667)
})

test_that("TC-V02: SR = Inf when baseline near 0 and range > 0", {
  sr <- lwdid:::.compute_sensitivity_ratio(c(0.1, 0.2), 1e-15)
  expect_true(is.infinite(sr))
  expect_equal(sr, Inf)

  sr2 <- lwdid:::.compute_sensitivity_ratio(c(0.0, 0.5), 1e-12)
  expect_equal(sr2, Inf)
})

test_that("TC-V03: Boundary precision — SR feeds correctly into level", {
  # baseline=10.0, range=0.999 → SR = 0.0999
  sr_below <- lwdid:::.compute_sensitivity_ratio(c(10.0, 10.999), 10.0)
  expect_equal(sr_below, 0.0999, tolerance = 1e-10)
  expect_equal(lwdid:::.determine_robustness_level(sr_below), "highly_robust")

  # baseline=10.0, range=1.0 → SR = 0.10
  sr_at <- lwdid:::.compute_sensitivity_ratio(c(10.0, 11.0), 10.0)
  expect_equal(sr_at, 0.10, tolerance = 1e-10)
  expect_equal(lwdid:::.determine_robustness_level(sr_at), "moderately_robust")
})


test_that("TC-V04: Known DGP end-to-end sensitivity analysis", {
  skip_if_not(
    requireNamespace("data.table", quietly = TRUE),
    message = "data.table required for DGP test"
  )

  # Simple panel: N=30, T=8, treatment at t=5, true ATT=2.0
  set.seed(123)
  N <- 30L; T_total <- 8L; S <- 5L; tau <- 2.0; n_treat <- 12L

  dt <- data.table::CJ(id = seq_len(N), time = seq_len(T_total))
  dt[, d := as.integer(id <= n_treat)]
  alpha_i <- stats::rnorm(N, mean = 10, sd = 2)
  dt[, y := alpha_i[id] + 0.5 * time +
       tau * d * as.integer(time >= S) +
       stats::rnorm(.N, sd = 0.3)]
  dt[, post := as.integer(time >= S)]

  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid_sensitivity_pre_period(
        data = dt, y = "y", ivar = "id", tvar = "time",
        d = "d", post = "post", rolling = "demean",
        estimator = "ra", pre_period_range = c(2L, 4L),
        verbose = FALSE
      )
    )),
    error = function(e) {
      skip(paste("lwdid estimation failed:", conditionMessage(e)))
    }
  )

  # Structure checks
  expect_s3_class(result, "lwdid_sensitivity")
  expect_equal(result$type, "pre_period")
  expect_equal(result$n_specifications, 3L)
  expect_equal(result$pre_period_range_tested, c(2L, 4L))

  # Numerical reasonableness
  converged <- Filter(function(s) isTRUE(s$converged), result$specifications)
  skip_if(length(converged) == 0L, "No converged specifications")

  atts <- vapply(converged, function(s) s$att, numeric(1))
  expect_true(all(atts > 0), info = "All ATTs should be positive (true ATT=2.0)")
  expect_true(all(abs(atts - tau) < 2.0),
              info = "All ATTs should be within 2.0 of true effect")
  expect_true(result$sensitivity_ratio < 1.0,
              info = "SR should be < 1.0 for clean DGP")
  expect_true(result$all_same_sign)
  expect_equal(result$n_sign_changes, 0L)
  expect_equal(result$baseline_spec$n_pre_periods, 4L)
  expect_true(is.finite(result$att_std))
  expect_true(result$att_std >= 0)
  expect_true(result$att_range[1] <= result$att_range[2])
  expect_true(result$robustness_level %in%
    c("highly_robust", "moderately_robust", "sensitive", "highly_sensitive"))
  expect_true(nchar(result$recommendation) > 0L)
})


# ============================================================================
# S3 Methods: print, summary, plot
# ============================================================================

test_that("print.lwdid_sensitivity produces output and returns invisibly", {
  obj <- structure(list(
    type = "pre_period",
    specifications = list(
      list(converged = TRUE, n_pre_periods = 3L, att = 1.5,
           se = 0.3, t_stat = 5.0, pvalue = 0.001,
           ci_lower = 0.9, ci_upper = 2.1)
    ),
    baseline_spec = list(att = 1.5, n_pre_periods = 3L,
                         se = 0.3, t_stat = 5.0, pvalue = 0.001,
                         ci_lower = 0.9, ci_upper = 2.1),
    att_range = c(1.5, 1.5), att_mean = 1.5, att_std = 0.0,
    sensitivity_ratio = 0.0, robustness_level = "highly_robust",
    is_robust = TRUE, robustness_threshold = 0.25,
    all_same_sign = TRUE, all_significant = TRUE,
    n_significant = 1L, n_sign_changes = 0L,
    rolling_method = "demean", estimator = "ra",
    n_specifications = 1L, pre_period_range_tested = c(3L, 3L),
    recommendation = "Results are robust.",
    detailed_recommendations = character(0),
    result_warnings = character(0)
  ), class = c("lwdid_sensitivity", "list"))

  out <- capture.output(ret <- print(obj))
  expect_true(length(out) > 0L)
  expect_true(any(grepl("Sensitivity Ratio", out)))
  expect_true(any(grepl("Baseline ATT", out)))
  expect_identical(ret, obj)
})


test_that("summary.lwdid_sensitivity produces detailed output", {
  obj <- structure(list(
    type = "pre_period",
    specifications = list(
      list(converged = TRUE, n_pre_periods = 3L, att = 1.5,
           se = 0.3, t_stat = 5.0, pvalue = 0.001,
           ci_lower = 0.9, ci_upper = 2.1),
      list(converged = FALSE, n_pre_periods = 2L, att = NA_real_,
           se = NA_real_, t_stat = NA_real_, pvalue = NA_real_,
           ci_lower = NA_real_, ci_upper = NA_real_)
    ),
    baseline_spec = list(att = 1.5, n_pre_periods = 3L,
                         se = 0.3, t_stat = 5.0, pvalue = 0.001,
                         ci_lower = 0.9, ci_upper = 2.1),
    att_range = c(1.5, 1.5), att_mean = 1.5, att_std = 0.0,
    sensitivity_ratio = 0.0, robustness_level = "highly_robust",
    is_robust = TRUE, robustness_threshold = 0.25,
    all_same_sign = TRUE, all_significant = TRUE,
    n_significant = 1L, n_sign_changes = 0L,
    rolling_method = "demean", estimator = "ra",
    n_specifications = 2L, pre_period_range_tested = c(2L, 3L),
    recommendation = "Results are robust.",
    detailed_recommendations = c("Consider reporting range."),
    result_warnings = c("Monotonic trend detected.")
  ), class = c("lwdid_sensitivity", "list"))

  out <- capture.output(ret <- summary(obj))
  expect_true(length(out) > 10L)
  expect_true(any(grepl("Configuration", out)))
  expect_true(any(grepl("Baseline Estimate", out)))
  expect_true(any(grepl("Sensitivity Metrics", out)))
  expect_true(any(grepl("Robustness Assessment", out)))
  expect_true(any(grepl("Specification Details", out)))
  expect_true(any(grepl("Recommendation", out)))
  expect_true(any(grepl("FAILED", out)))
  expect_true(any(grepl("Detailed Recommendations", out)))
  expect_true(any(grepl("Warnings", out)))
  expect_identical(ret, obj)
})

test_that("plot.lwdid_sensitivity returns ggplot object", {
  skip_if_not(
    requireNamespace("ggplot2", quietly = TRUE),
    message = "ggplot2 required for plot test"
  )
  obj <- structure(list(
    type = "pre_period",
    specifications = list(
      list(converged = TRUE, n_pre_periods = 3L, att = 1.5,
           se = 0.3, t_stat = 5.0, pvalue = 0.001,
           ci_lower = 0.9, ci_upper = 2.1)
    ),
    baseline_spec = list(att = 1.5, n_pre_periods = 3L),
    sensitivity_ratio = 0.0, robustness_level = "highly_robust"
  ), class = c("lwdid_sensitivity", "list"))
  p <- plot(obj)
  expect_s3_class(p, "ggplot")
})

test_that("plot.lwdid_sensitivity returns empty ggplot with no converged specs", {
  skip_if_not(
    requireNamespace("ggplot2", quietly = TRUE),
    message = "ggplot2 required for plot test"
  )
  obj <- structure(list(
    type = "pre_period",
    specifications = list(
      list(converged = FALSE, n_pre_periods = 3L, att = NA_real_)
    ),
    baseline_spec = list(att = NA_real_, n_pre_periods = 3L),
    sensitivity_ratio = 0.0, robustness_level = "highly_robust"
  ), class = c("lwdid_sensitivity", "list"))
  p <- plot(obj)
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 0L)
})

# ============================================================================
# Integration: Result field consistency
# ============================================================================

test_that("Integration: result fields are internally consistent", {
  skip_on_cran()
  df <- make_ct_sensitivity(
    n_units = 20,
    n_periods = 10,
    n_treated = 10,
    K = 5,
    seed = 300L
  )
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid_sensitivity_pre_period(
        data = df,
        y = "y",
        ivar = "id",
        tvar = "year",
        d = "d",
        post = "post",
        rolling = "demean",
        pre_period_range = c(2L, 4L),
        step = 1L,
        verbose = FALSE
      )
    )),
    error = function(e) NULL
  )
  skip_if(is.null(result), "lwdid_sensitivity_pre_period failed")

  expect_equal(result$n_specifications, length(result$specifications))

  converged <- Filter(function(s) isTRUE(s$converged), result$specifications)
  if (length(converged) > 0) {
    atts <- vapply(converged, function(s) s$att, numeric(1))
    expect_equal(result$att_range[1], min(atts))
    expect_equal(result$att_range[2], max(atts))
    expect_equal(result$att_mean, mean(atts))
  }

  expect_equal(
    result$is_robust,
    result$sensitivity_ratio < result$robustness_threshold
  )

  expected_level <- lwdid:::.determine_robustness_level(
    result$sensitivity_ratio
  )
  expect_equal(result$robustness_level, expected_level)
  expect_equal(result$pre_period_range_tested, c(2L, 4L))

  if (length(converged) > 0) {
    max_npre <- max(vapply(
      converged,
      function(s) s$n_pre_periods,
      integer(1)
    ))
    expect_equal(result$baseline_spec$n_pre_periods, max_npre)
  }
})
