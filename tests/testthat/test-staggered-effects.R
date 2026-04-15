# =============================================================================
# Tests for estimate_staggered_effects() - Core Functionality
# Verifies return structure, ordering, event_time, and (g,r) pairs
# =============================================================================

# --- Helper: build balanced staggered panel with known ATT ---
make_staggered_test_data <- function(seed = 12345, tau = 5.0, sigma = 0.01) {
  set.seed(seed)
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 5), rep(5L, 5), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, alpha_i := as.numeric(id)]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := alpha_i + tau * treated + rnorm(.N, 0, sigma)]
  dt[, c("alpha_i", "treated") := NULL]
  dt
}

# --- Test 1: returns 17-column data.frame ---
test_that("returns 17-column data.frame", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 17L)
  expected_cols <- c(
    "cohort", "period", "event_time", "att", "se",
    "ci_lower", "ci_upper", "t_stat", "pvalue",
    "df", "df_inference", "n", "n_treated",
    "n_control", "K", "controls_tier", "vce_type")
  expect_equal(names(result), expected_cols)
})

# --- Test 2: results ordered by (cohort, period) ---
test_that("results ordered by (cohort, period)", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_true(all(diff(result$cohort) >= 0))
  for (g in unique(result$cohort)) {
    sub <- result[result$cohort == g, ]
    expect_true(all(diff(sub$period) > 0))
  }
})

# --- Test 3: event_time equals period minus cohort ---
test_that("event_time equals period minus cohort", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_equal(result$event_time, result$period - result$cohort)
  expect_true(all(result$event_time >= 0))
})

# --- Test 4: all valid (g,r) pairs estimated ---
test_that("all valid (g,r) pairs estimated", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_equal(nrow(result), 8L)
  expect_equal(sum(result$cohort == 3), 5L)
  expect_equal(sum(result$cohort == 5), 3L)
  expect_equal(result$period[result$cohort == 3], 3:7)
  expect_equal(result$period[result$cohort == 5], 5:7)
})


# =============================================================================
# Group 2: Skip Reasons and Error Tolerance Tests (E4-04.4)
# =============================================================================

# --- Test 5: AET last cohort has no controls and is skipped ---
test_that("AET last cohort has no controls and is skipped", {
  set.seed(77701)
  # All-eventually-treated data: cohorts 3 and 5, NO never-treated units.
  # With control_group = "not_yet_treated":
  #   Cohort 3 at r=3,4: cohort 5 (G=5 > r) serves as NYT control -> OK
  #   Cohort 5 at r=5,6,7: no units with G > r exist -> no controls -> skip
  units <- 1:20
  periods <- 1:7
  # 10 units in cohort 3, 10 units in cohort 5 (no NT units)
  gvar_map <- c(rep(3L, 10), rep(5L, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Should produce results for cohort 3 (r=3,4) and warn about skipped pairs
  w_captured <- NULL
  result <- withCallingHandlers(
    estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats
    ),
    warning = function(w) {
      if (inherits(w, "lwdid_data")) {
        w_captured <<- w
      }
      invokeRestart("muffleWarning")
    }
  )

  # Cohort 3 should have results (at least r=3,4 where cohort 5 is NYT)
  expect_true(nrow(result) > 0L)
  expect_true(all(result$cohort == 3L) || any(result$cohort == 3L))
  # Cohort 5 should have NO results (no controls available)
  expect_equal(sum(result$cohort == 5L), 0L)
  # Warning should have been issued with gr_pairs_skipped detail
  expect_false(is.null(w_captured))
  expect_equal(w_captured$detail, "gr_pairs_skipped")
})


# --- Test 6: min_treated not met causes skip ---
test_that("min_treated not met causes skip", {
  set.seed(77702)
  # Cohort 3: only 1 unit (id=1), cohort 5: 5 units (id=2..6), NT: 10 units
  units <- 1:16
  periods <- 1:7
  gvar_map <- c(3L, rep(5L, 5), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # min_treated = 2 should skip cohort 3 (only 1 treated unit)
  w_captured <- NULL
  result <- withCallingHandlers(
    estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_treated = 2L
    ),
    warning = function(w) {
      if (inherits(w, "lwdid_data")) {
        w_captured <<- w
      }
      invokeRestart("muffleWarning")
    }
  )

  # Cohort 3 pairs should all be skipped (only 1 treated unit < min_treated=2)
  expect_equal(sum(result$cohort == 3L), 0L)
  # Cohort 5 should still have results (5 treated units >= 2)
  expect_true(sum(result$cohort == 5L) > 0L)
  # Warning issued about skipped pairs
  expect_false(is.null(w_captured))
  expect_equal(w_captured$detail, "gr_pairs_skipped")
})


# --- Test 7: min_obs not met causes lwdid_insufficient_data error ---
test_that("min_obs not met causes lwdid_insufficient_data error", {
  set.seed(77703)
  # Small dataset: 5 treated + 5 NT, 7 periods
  units <- 1:10
  periods <- 1:7
  gvar_map <- c(rep(3L, 5), rep(Inf, 5))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # min_obs = 100 is way above the 10 units per period -> all pairs skipped
  # When ALL pairs are skipped, function throws lwdid_insufficient_data
  expect_error(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_obs = 100L
    )),
    class = "lwdid_insufficient_data"
  )
})


# --- Test 8: OLS failure caught by tryCatch (estimation_error skip) ---
test_that("OLS failure caught by tryCatch does not crash loop", {
  set.seed(77704)
  # Strategy: create data where after transformation, all y_trans values

  # are identical (zero variance) for one cohort, causing OLS to fail
  # or produce degenerate results. We use a cohort with very few units
  # and engineer the outcome so that after demeaning, the transformed
  # values become constant (all units have identical Y trajectories).
  units <- 1:15
  periods <- 1:7
  # Cohort 3: 2 units with IDENTICAL outcomes (will cause issues after transform)
  # Cohort 5: 5 normal units, NT: 8 units
  gvar_map <- c(rep(3L, 2), rep(5L, 5), rep(Inf, 8))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  # Cohort 3 units: identical Y values (same alpha_i = 100 for both)
  # This means after demeaning, y_trans for treated units will be identical
  dt[, Y := data.table::fifelse(
    gvar == 3L,
    100.0 + 5.0 * treated,
    as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.5)
  )]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # The function should NOT crash even if cohort 3 estimation has issues.
  # It should either skip those pairs or produce results.
  # Cohort 5 should definitely produce valid results.
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)
  # Cohort 5 should have valid results regardless
  expect_true(sum(result$cohort == 5L) > 0L)
})


# --- Test 9: Skip summary warning contains correct detail and attributes ---
test_that("skip summary warning contains correct detail and skipped_pairs", {
  set.seed(77705)
  # Same AET setup as test 5: cohort 5 pairs will be skipped (no NYT controls)
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 10), rep(5L, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Capture the warning object to inspect its attributes
  w_captured <- NULL
  result <- withCallingHandlers(
    estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats
    ),
    warning = function(w) {
      if (inherits(w, "lwdid_data")) {
        w_captured <<- w
      }
      invokeRestart("muffleWarning")
    }
  )

  # Warning must exist
  expect_false(is.null(w_captured))
  # Check class hierarchy
  expect_true(inherits(w_captured, "lwdid_data"))
  expect_true(inherits(w_captured, "lwdid_warning"))
  # Check detail field
  expect_equal(w_captured$detail, "gr_pairs_skipped")
  # Check skipped_pairs attribute is a list of lists with g, r, reason
  sp <- w_captured$skipped_pairs
  expect_true(is.list(sp))
  expect_true(length(sp) > 0L)
  # Each skipped pair should have g, r, reason fields
  for (pair in sp) {
    expect_true("g" %in% names(pair))
    expect_true("r" %in% names(pair))
    expect_true("reason" %in% names(pair))
  }
  # At least some skipped pairs should be from cohort 5 with "no_control" reason
  skipped_g5 <- Filter(function(p) p$g == 5L, sp)
  expect_true(length(skipped_g5) > 0L)
  reasons_g5 <- vapply(skipped_g5, function(p) p$reason, character(1))
  expect_true("no_control" %in% reasons_g5)
  # Warning message should mention "skipped"
  expect_true(grepl("skipped", w_captured$message, ignore.case = TRUE))
})


# --- Test 10: No valid results throws lwdid_insufficient_data ---
test_that("no valid results throws lwdid_insufficient_data", {
  set.seed(77706)
  # Create data where ALL (g,r) pairs will be skipped.
  # Use a single cohort with very high min_obs so every pair fails.
  units <- 1:6
  periods <- 1:5
  gvar_map <- c(rep(3L, 3), rep(Inf, 3))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # min_obs = 1000 ensures all pairs are skipped -> no valid results
  expect_error(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_obs = 1000L
    )),
    class = "lwdid_insufficient_data"
  )

  # Also verify the error message mentions "No valid"
  err <- tryCatch(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_obs = 1000L
    )),
    lwdid_insufficient_data = function(e) e
  )
  expect_true(grepl("No valid", err$message))
})


# =============================================================================
# Group 3: NA Protection Tests (E4-04.5)
# Verifies that NT units produce d_sub = 0 (not NA), and that units
# with NA y_trans are filtered without crashing estimation.
# =============================================================================

# --- Test 11: NT units have d_sub = 0 (not NA) ---
test_that("NT units have d_sub = 0 via is_never_treated short-circuit", {
  set.seed(88801)
  # Panel: 5 treated (cohort 3), 10 NT (gvar = Inf), 5 NT (gvar = 0)
  # The NA protection pattern:
  #   !is_never_treated(gvar) & gvar == g
  # For NT units: !TRUE & (Inf == 3) -> FALSE & NA -> FALSE (not NA)
  # For gvar=0:   !TRUE & (0 == 3)   -> FALSE & FALSE -> FALSE
  # So NT units get d_sub = 0, not NA.
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 5), rep(Inf, 10), rep(0L, 5))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    !is_never_treated(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)
  # NT units (gvar=Inf and gvar=0) should be counted as controls
  # n_control > 0 confirms they were NOT dropped due to NA d_sub
  expect_true(all(result$n_control > 0L))
  # n_treated should be exactly 5 (the cohort 3 units)
  expect_true(all(result$n_treated == 5L))
  # ATT should be close to 5.0 (the true treatment effect)
  expect_true(all(abs(result$att - 5.0) < 0.5))
})


# --- Test 12: All-NA transform units are filtered, estimation continues ---
test_that("all-NA transform units are filtered, estimation continues", {
  set.seed(88802)
  # Unbalanced panel: units 11-12 appear only in periods 5-7 (no pre-period
  # data for cohort 3 whose pre-period is t < 3). After precompute_transforms,

  # these units won't have pre_stats entries. apply_precomputed_transform
  # uses match() which returns NA for unmatched units -> y_trans = NA.
  # The valid_trans filter in Step 5 should remove them.
  units_full <- 1:10   # full panel
  units_late <- 11:12  # late entrants (only periods 5-7)
  periods <- 1:7

  # Full panel units: 5 treated (cohort 3), 5 NT (Inf)
  gvar_full <- c(rep(3L, 5), rep(Inf, 5))
  dt_full <- data.table::CJ(id = units_full, time = periods)
  dt_full[, gvar := gvar_full[id]]

  # Late entrants: NT units appearing only in periods 5-7
  dt_late <- data.table::CJ(id = units_late, time = 5:7)
  dt_late[, gvar := Inf]

  dt <- data.table::rbindlist(list(dt_full, dt_late))
  dt[, treated := data.table::fifelse(
    !is_never_treated(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Units 11-12 should NOT be in pre_stats (no pre-period rows)
  expect_false(11L %in% pre_stats[["3"]]$id)
  expect_false(12L %in% pre_stats[["3"]]$id)

  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)
  # For periods 5-7, units 11-12 would be in the cross-section but
  # their y_trans is NA -> filtered out. So n should be less than
  # the total units present in those periods.
  late_periods <- result[result$period >= 5L, ]
  # Full panel has 10 units; with late entrants, cross-section has 12.
  # After NA filter, n should be 10 (late entrants removed).
  expect_true(all(late_periods$n <= 10L))
  # ATT should still be close to 5.0
  expect_true(all(abs(result$att - 5.0) < 0.5))
})


# --- Test 13: valid_trans filter reduces n count ---
test_that("valid_trans filter reduces n count for late entrants", {
  set.seed(88803)
  # Create a combined panel with late entrants that have no
  # pre-period data. The valid_trans filter should remove them,
  # so n in the result should reflect only units with valid
  # y_trans (i.e., those present in pre-period).
  units_full <- 1:10
  units_late <- 13:15  # 3 late entrants (periods 5-7 only)
  periods <- 1:7

  gvar_full <- c(rep(3L, 5), rep(Inf, 5))
  dt_full <- data.table::CJ(id = units_full, time = periods)
  dt_full[, gvar := gvar_full[id]]

  dt_late <- data.table::CJ(id = units_late, time = 5:7)
  dt_late[, gvar := Inf]

  dt_combined <- data.table::rbindlist(list(dt_full, dt_late))
  dt_combined[, treated := data.table::fifelse(
    !is_never_treated(gvar) & time >= gvar, 1L, 0L
  )]
  dt_combined[, Y := as.numeric(id) + 5.0 * treated +
    rnorm(.N, 0, 0.01)]
  dt_combined[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt_combined, "Y", "id", "time",
    cohorts = cohorts, rolling = "demean"
  )

  # Late entrants should NOT be in pre_stats
  expect_false(any(c(13L, 14L, 15L) %in% pre_stats[["3"]]$id))

  # Run estimation on combined panel
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt_combined, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # For periods 3-4: only full-panel units present (10 units)
  for (r in 3:4) {
    n_r <- result$n[result$period == r]
    expect_equal(n_r, 10L,
      info = sprintf("period %d: n=%d, expected 10", r, n_r))
  }

  # For periods 5-7: cross-section has 13 units (10 full + 3 late)
  # but after valid_trans filter, n should be 10 (late entrants
  # have NA y_trans and are removed)
  for (r in 5:7) {
    n_r <- result$n[result$period == r]
    # n should be 10, NOT 13 — late entrants filtered out
    expect_equal(n_r, 10L,
      info = sprintf("period %d: n=%d, expected 10", r, n_r))
  }
})


# =============================================================================
# Group 4: VCE Parameter Tests (E4-04.5)
# Verifies vce_type column reflects the vce argument correctly,
# and that NULL vce does not cause errors in cluster extraction.
# =============================================================================

# --- Test 14: vce=NULL produces "homoskedastic" vce_type ---
test_that("vce=NULL produces homoskedastic vce_type", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_true(all(result$vce_type == "homoskedastic"))
})


# --- Test 15: vce="hc1" produces "hc1" vce_type ---
test_that("vce='hc1' produces hc1 vce_type", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = "hc1",
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_true(all(result$vce_type == "hc1"))
})


# --- Test 16: identical() protects NULL vce in cluster extraction ---
test_that("identical() protects NULL vce in cluster extraction", {
  # The function uses identical(vce, "cluster") instead of
  # vce == "cluster" to prevent NULL == "cluster" returning
  # logical(0), which would cause if() to error.
  # This test explicitly verifies vce=NULL + cluster_var=NULL
  # does not error and produces valid results.
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # Should NOT error — identical(NULL, "cluster") returns FALSE
  expect_no_error(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats
    ))
  )

  # Also verify the result is valid
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)
  # No NA in critical columns
  expect_false(any(is.na(result$att)))
  expect_false(any(is.na(result$se)))
  expect_false(any(is.na(result$vce_type)))
})


# =============================================================================
# Group 5: Demean/Detrend numerical consistency tests
# =============================================================================

# --- Test 5.1: Known ATT=5 DGP with demean recovers ATT within 0.5 ---
test_that("demean: known ATT=5 DGP recovers ATT within 0.5", {
  set.seed(50501)
  units <- 1:20
  periods <- 1:7
  # 10 units in cohort 3, 10 NT units
  gvar_map <- c(rep(3L, 10), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, alpha_i := as.numeric(id)]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := alpha_i + 5.0 * D + rnorm(.N, 0, 0.01)]
  dt[, c("alpha_i", "D") := NULL]

  cohorts <- 3L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # 5 post-treatment periods: 3,4,5,6,7

  expect_equal(nrow(result), 5L)
  # All ATT estimates within 0.5 of true ATT=5.0
  expect_true(all(abs(result$att - 5.0) < 0.5))
  # All standard errors positive
  expect_true(all(result$se > 0))
  # Confidence intervals contain the point estimate
  expect_true(all(result$ci_lower < result$att))
  expect_true(all(result$att < result$ci_upper))
})


# --- Test 5.2: Known ATT=3 DGP with detrend recovers ATT within 0.5 ---
test_that("detrend: known ATT=3 DGP with unit trends recovers ATT within 0.5", {
  set.seed(50502)
  units <- 1:20
  periods <- 1:8
  # 10 units in cohort 4, 10 NT units
  gvar_map <- c(rep(4L, 10), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, alpha_i := as.numeric(id)]
  dt[, beta_i := 0.5 * as.numeric(id)]
  dt[, D := data.table::fifelse(is.finite(gvar) & time >= gvar, 1L, 0L)]
  dt[, Y := alpha_i + beta_i * time + 3.0 * D + rnorm(.N, 0, 0.01)]
  dt[, c("alpha_i", "beta_i", "D") := NULL]

  cohorts <- 4L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "detrend"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "detrend",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # 5 post-treatment periods: 4,5,6,7,8
  expect_equal(nrow(result), 5L)
  # All ATT estimates within 0.5 of true ATT=3.0
  expect_true(all(abs(result$att - 3.0) < 0.5))
  # All standard errors positive
  expect_true(all(result$se > 0))
})


# --- Test 5.3: Zero treatment effect DGP, estimates near zero ---
test_that("demean: zero treatment effect DGP produces estimates near zero", {
  set.seed(50503)
  units <- 1:20
  periods <- 1:7
  # 8 units in cohort 3, 12 NT units
  gvar_map <- c(rep(3L, 8), rep(Inf, 12))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, alpha_i := as.numeric(id)]
  dt[, Y := alpha_i + rnorm(.N, 0, 0.1)]
  dt[, c("alpha_i") := NULL]

  cohorts <- 3L
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # 5 post-treatment periods: 3,4,5,6,7
  expect_equal(nrow(result), 5L)
  # All ATT estimates near zero (true ATT = 0)
  expect_true(all(abs(result$att) < 1.0))
  # All standard errors positive
  expect_true(all(result$se > 0))
})


# =============================================================================
# Group 6: Multi-cohort and control group strategy tests
# =============================================================================

# --- Test 6.1: Multi-cohort independent estimation, each cohort has results ---
test_that("multi-cohort: each cohort produces independent results", {
  set.seed(10701)
  units <- 1:30
  periods <- 1:10
  # 8 units in cohort 3, 7 in cohort 5, 5 in cohort 7, 10 NT
  gvar_map <- c(rep(3L, 8), rep(5L, 7), rep(7L, 5), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L, 7L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Each cohort should have results
  expect_true(any(result$cohort == 3L))
  expect_true(any(result$cohort == 5L))
  expect_true(any(result$cohort == 7L))
  # Cohort 3: periods 3-10 = 8 pairs
  expect_equal(sum(result$cohort == 3L), 8L)
  # Cohort 5: periods 5-10 = 6 pairs
  expect_equal(sum(result$cohort == 5L), 6L)
  # Cohort 7: periods 7-10 = 4 pairs
  expect_equal(sum(result$cohort == 7L), 4L)
  # Total: 18 pairs
  expect_equal(nrow(result), 18L)
})

# --- Test 6.2: never_treated strategy: control group size constant ---
test_that("never_treated: control group size constant across periods", {
  set.seed(10702)
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 5), rep(5L, 5), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "never_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # With never_treated, control group = 10 NT units for ALL (g,r) pairs
  expect_true(all(result$n_control == 10L))
  # n_treated should match cohort sizes
  expect_true(all(result$n_treated[result$cohort == 3] == 5L))
  expect_true(all(result$n_treated[result$cohort == 5] == 5L))
})

# =============================================================================
# Group 7: min parameter configuration tests
# =============================================================================

# --- Test 7.1: min_obs custom value correctly filters ---
test_that("min_obs custom value correctly filters pairs", {
  set.seed(10703)
  # 3 treated + 5 NT = 8 units total per (g,r) pair
  units <- 1:8
  periods <- 1:7
  gvar_map <- c(rep(3L, 3), rep(Inf, 5))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # min_obs = 10: all pairs skipped (8 obs < 10)
  expect_error(
    suppressWarnings(estimate_staggered_effects(
      dt = dt, y = "Y", ivar = "id", tvar = "time",
      gvar = "gvar", rolling = "demean",
      control_group = "not_yet_treated",
      controls = NULL, vce = NULL,
      cluster_var = NULL, alpha = 0.05,
      pre_stats = pre_stats,
      min_obs = 10L
    )),
    class = "lwdid_insufficient_data"
  )

  # min_obs = 5: all pairs should succeed (8 obs >= 5)
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats,
    min_obs = 5L
  ))
  expect_equal(nrow(result), 5L)  # periods 3-7
})

# --- Test 7.2: min_treated=1, min_control=1 allows small samples ---
test_that("min_treated=1 min_control=1 allows small samples", {
  set.seed(10704)
  # 1 treated unit + 2 NT units = 3 total per (g,r) pair
  units <- 1:3
  periods <- 1:5
  gvar_map <- c(3L, Inf, Inf)
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  cohorts <- c(3L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )

  # With min_treated=1, min_control=1, min_obs=3: should work
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats,
    min_treated = 1L, min_control = 1L, min_obs = 3L
  ))

  expect_true(nrow(result) > 0L)
  expect_true(all(result$n_treated == 1L))
  expect_true(all(result$n_control == 2L))
  expect_true(all(result$n == 3L))
  # ATT should still be reasonable
  expect_true(all(abs(result$att - 5.0) < 1.0))
})

# =============================================================================
# Group 8: End-to-End Integration Tests (E4-04.8)
# =============================================================================

# --- Test 8.1: Demean end-to-end pipeline works ---
test_that("E2E demean: full pipeline with known ATT=5 recovers treatment effect", {
  dt <- make_staggered_test_data(seed = 20801, tau = 5.0, sigma = 0.01)
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Result is a valid data.frame with rows

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)

  # ATT within 0.5 of known tau = 5.0
  expect_true(all(abs(result$att - 5.0) < 0.5),
    info = paste("ATT values:", paste(round(result$att, 4), collapse = ", ")))

  # Standard errors positive
  expect_true(all(result$se > 0))

  # Confidence intervals bracket ATT
  expect_true(all(result$ci_lower < result$att))
  expect_true(all(result$ci_upper > result$att))

  # p-values in [0, 1]
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))

  # Sample sizes positive
  expect_true(all(result$n > 0L))
  expect_true(all(result$n_treated > 0L))
  expect_true(all(result$n_control > 0L))
})

# --- Test 8.2: Detrend end-to-end pipeline works ---
test_that("E2E detrend: unit-specific linear trends removed, ATT=3 recovered", {
  set.seed(20802)
  units <- 1:20
  periods <- 1:7
  # Cohorts 3 and 5, 10 NT units
  gvar_map <- c(rep(3L, 5), rep(5L, 5), rep(Inf, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  # DGP: Y = alpha_i + beta_i * time + tau * D_it + epsilon
  dt[, alpha_i := as.numeric(id)]
  dt[, beta_i := 0.5 * as.numeric(id)]
  dt[, D_it := data.table::fifelse(
    is.finite(gvar) & time >= gvar, 1L, 0L
  )]
  dt[, Y := alpha_i + beta_i * time + 3.0 * D_it + rnorm(.N, 0, 0.1)]
  dt[, c("alpha_i", "beta_i", "D_it") := NULL]

  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "detrend"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "detrend",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # Result is a valid data.frame with rows
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0L)

  # ATT within 1.0 of known tau = 3.0 (detrend has more variance)
  expect_true(all(abs(result$att - 3.0) < 1.0),
    info = paste("ATT values:", paste(round(result$att, 4), collapse = ", ")))

  # Standard errors positive
  expect_true(all(result$se > 0))

  # p-values in [0, 1]
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
})

# --- Test 8.3: Result field types correct and numerically reasonable ---
test_that("E2E field types: correct column types and numerical reasonableness", {
  dt <- make_staggered_test_data()
  cohorts <- c(3L, 5L)
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time", cohorts = cohorts, rolling = "demean"
  )
  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt, y = "Y", ivar = "id", tvar = "time",
    gvar = "gvar", rolling = "demean",
    control_group = "not_yet_treated",
    controls = NULL, vce = NULL,
    cluster_var = NULL, alpha = 0.05,
    pre_stats = pre_stats
  ))

  # --- Column type checks ---
  # cohort, period, event_time: integer or numeric whole numbers
  expect_true(is.numeric(result$cohort))
  expect_true(is.numeric(result$period))
  expect_true(is.integer(result$event_time) || is.numeric(result$event_time))
  if (is.numeric(result$event_time) && !is.integer(result$event_time)) {
    expect_true(all(result$event_time == as.integer(result$event_time)))
  }

  # att, se, ci_lower, ci_upper, t_stat, pvalue, df, df_inference: numeric
  expect_true(is.numeric(result$att))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(is.numeric(result$t_stat))
  expect_true(is.numeric(result$pvalue))
  expect_true(is.numeric(result$df))
  expect_true(is.numeric(result$df_inference))

  # n, n_treated, n_control, K: integer or numeric whole numbers
  expect_true(is.numeric(result$n))
  expect_true(is.numeric(result$n_treated))
  expect_true(is.numeric(result$n_control))
  expect_true(is.numeric(result$K))

  # controls_tier: character or integer
  expect_true(is.character(result$controls_tier) || is.integer(result$controls_tier))

  # vce_type: character
  expect_true(is.character(result$vce_type))

  # --- Numerical reasonableness checks ---
  # pvalue in [0, 1] for all rows
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))

  # se >= 0 for all rows
  expect_true(all(result$se >= 0))

  # ci_lower < ci_upper for all rows
  expect_true(all(result$ci_lower < result$ci_upper))

  # n == n_treated + n_control for all rows
  expect_true(all(result$n == result$n_treated + result$n_control))

  # event_time == period - cohort for all rows
  expect_true(all(result$event_time == result$period - result$cohort))

  # t_stat approximately equals att / se (within 1e-6 relative tolerance)
  expected_t <- result$att / result$se
  expect_true(all(abs(result$t_stat - expected_t) <= 1e-6 * abs(expected_t) + 1e-10),
    info = paste("t_stat mismatch: max rel diff =",
      max(abs(result$t_stat - expected_t) / (abs(expected_t) + 1e-10))))

  # df > 0 for all rows
  expect_true(all(result$df > 0))
})
