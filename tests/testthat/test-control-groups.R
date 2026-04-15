# test-control-groups.R
# Tests for control group selection utilities (control_groups.R)
# Note: is_never_treated() and get_cohorts() are tested in test-utils.R.
# This file tests check_never_treated() from control_groups.R,
# plus additional integration tests for the control group pipeline.

# ============================================================================
# Group 1: is_never_treated() — Five NT marker values
# (Adapted to match actual utils.R behavior)
# ============================================================================

test_that("is_never_treated identifies NA as NT", {
  expect_true(is_never_treated(NA))
  expect_true(is_never_treated(NA_real_))
  expect_true(is_never_treated(NA_integer_))
})

test_that("is_never_treated identifies NaN as NT", {
  expect_true(is_never_treated(NaN))
})

test_that("is_never_treated identifies 0 as NT", {
  expect_true(is_never_treated(0))
  expect_true(is_never_treated(0L))
  expect_true(is_never_treated(0.0))
})

test_that("is_never_treated identifies Inf as NT", {
  expect_true(is_never_treated(Inf))
})

test_that("is_never_treated rejects -Inf with error", {
  # utils.R implementation throws lwdid_invalid_staggered_data for -Inf
  expect_error(
    is_never_treated(-Inf),
    class = "lwdid_invalid_staggered_data"
  )
})

# ============================================================================
# Group 2: is_never_treated() — Normal cohort values
# ============================================================================

test_that("is_never_treated returns FALSE for positive integers", {
  expect_equal(is_never_treated(c(3, 5, 7)), c(FALSE, FALSE, FALSE))
})

test_that("is_never_treated returns FALSE for negative values", {
  # Negative values are unusual but not NT markers
  expect_equal(is_never_treated(c(-1, -2)), c(FALSE, FALSE))
})

test_that("is_never_treated returns FALSE for float cohort values", {
  expect_equal(is_never_treated(c(3.0, 5.0)), c(FALSE, FALSE))
})

# ============================================================================
# Group 3: is_never_treated() — Mixed vectors and edge cases
# ============================================================================

test_that("is_never_treated handles mixed vector correctly", {
  # Note: -Inf excluded since it throws error in utils.R implementation
  g <- c(NA, 0, 3, Inf, 5, NaN)
  expected <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  expect_equal(is_never_treated(g), expected)
})

test_that("is_never_treated handles empty vector", {
  expect_equal(is_never_treated(numeric(0)), logical(0))
})

test_that("is_never_treated handles all valid NT markers together", {
  # NA, NaN, 0, Inf (not -Inf which throws error)
  g <- c(NA, NaN, 0, Inf)
  expect_true(all(is_never_treated(g)))
})

test_that("is_never_treated handles year-based cohort values", {
  # BC-11: large numeric values used as year-based cohorts
  g <- c(2005, 2010, 2020)
  expect_equal(is_never_treated(g), c(FALSE, FALSE, FALSE))
})

# ============================================================================
# Helper: create staggered panel data
# ============================================================================

make_staggered_panel <- function(
  cohort_sizes = c(g3 = 3, g5 = 2),
  n_nt = 2,
  n_periods = 6
) {
  units <- list()
  uid <- 1L
  # Treated units
  for (i in seq_along(cohort_sizes)) {
    g_val <- as.integer(gsub("g", "", names(cohort_sizes)[i]))
    for (j in seq_len(cohort_sizes[i])) {
      units[[uid]] <- data.table::data.table(
        id = uid, time = seq_len(n_periods), gvar = g_val
      )
      uid <- uid + 1L
    }
  }
  # NT units (use NA as NT marker)
  for (j in seq_len(n_nt)) {
    units[[uid]] <- data.table::data.table(
      id = uid, time = seq_len(n_periods), gvar = NA_real_
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

# ============================================================================
# Group 4: get_cohorts() — Basic functionality
# ============================================================================

test_that("get_cohorts extracts sorted unique cohorts", {
  dt <- make_staggered_panel(c(g5 = 2, g3 = 3), n_nt = 1)
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, c(3L, 5L))
})

test_that("get_cohorts returns integer vector", {
  dt <- make_staggered_panel(c(g3 = 2), n_nt = 1)
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_true(is.integer(cohorts))
})

test_that("get_cohorts filters NT markers (NA, 0, Inf)", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    gvar = rep(c(3, NA, 0, Inf), each = 3)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, 3L)
})

test_that("get_cohorts handles single cohort", {
  dt <- make_staggered_panel(c(g4 = 5), n_nt = 3)
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, 4L)
  expect_equal(length(cohorts), 1L)
})

test_that("get_cohorts handles many cohorts", {
  dt <- data.table::data.table(
    id = rep(1:10, each = 5),
    time = rep(1:5, 10),
    gvar = rep(c(2, 3, 4, 5, 6, 7, 8, 9, NA, 0), each = 5)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, 2L:9L)
})

test_that("get_cohorts handles year-based cohort values", {
  dt <- data.table::data.table(
    id = rep(1:5, each = 10),
    time = rep(2001:2010, 5),
    gvar = rep(c(2005, 2007, 2009, NA, Inf), each = 10)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, c(2005L, 2007L, 2009L))
  expect_true(is.integer(cohorts))
})

# ============================================================================
# Group 5: get_cohorts() — Edge cases
# ============================================================================

test_that("get_cohorts returns integer(0) when all NT", {
  # utils.R version returns integer(0), not error
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    gvar = rep(c(NA, 0, Inf), each = 3)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, integer(0))
})

test_that("get_cohorts returns integer(0) when all NA", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    time = rep(1:3, 2),
    gvar = NA_real_
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, integer(0))
})

# ============================================================================
# Group 6: check_never_treated() — Basic functionality
# ============================================================================

test_that("check_never_treated returns correct structure", {
  dt <- make_staggered_panel(c(g3 = 3, g5 = 2), n_nt = 2)
  result <- check_never_treated(dt, "id", "gvar")
  expect_true(is.list(result))
  expect_named(result, c("has_nt", "n_nt", "n_treated", "cohort_sizes"))
})

test_that("check_never_treated counts correctly in mixed scenario", {
  dt <- make_staggered_panel(c(g3 = 3, g5 = 2), n_nt = 2)
  result <- check_never_treated(dt, "id", "gvar")
  expect_true(result$has_nt)
  expect_equal(result$n_nt, 2L)
  expect_equal(result$n_treated, 5L)
  expect_equal(as.integer(result$cohort_sizes[["3"]]), 3L)
  expect_equal(as.integer(result$cohort_sizes[["5"]]), 2L)
})

test_that("check_never_treated handles all-NT scenario", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    gvar = rep(c(NA, 0, Inf), each = 3)
  )
  result <- check_never_treated(dt, "id", "gvar")
  expect_true(result$has_nt)
  expect_equal(result$n_nt, 3L)
  expect_equal(result$n_treated, 0L)
  expect_equal(length(result$cohort_sizes), 0L)
})

test_that("check_never_treated handles no-NT scenario", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    gvar = rep(c(3, 5, 7), each = 3)
  )
  result <- check_never_treated(dt, "id", "gvar")
  expect_false(result$has_nt)
  expect_equal(result$n_nt, 0L)
  expect_equal(result$n_treated, 3L)
})

test_that("check_never_treated deduplicates at unit level", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 5),
    time = rep(1:5, 2),
    gvar = rep(c(3, NA), each = 5)
  )
  result <- check_never_treated(dt, "id", "gvar")
  expect_equal(result$n_nt, 1L)
  expect_equal(result$n_treated, 1L)
})

test_that("check_never_treated handles single unit scenario", {
  dt <- data.table::data.table(
    id = rep(1L, 4),
    time = 1:4,
    gvar = rep(5, 4)
  )
  result <- check_never_treated(dt, "id", "gvar")
  expect_false(result$has_nt)
  expect_equal(result$n_nt, 0L)
  expect_equal(result$n_treated, 1L)
  expect_equal(result$n_nt + result$n_treated, 1L)

  # Single NT unit
  dt_nt <- data.table::data.table(
    id = rep(1L, 4),
    time = 1:4,
    gvar = rep(NA_real_, 4)
  )
  result_nt <- check_never_treated(dt_nt, "id", "gvar")
  expect_true(result_nt$has_nt)
  expect_equal(result_nt$n_nt, 1L)
  expect_equal(result_nt$n_treated, 0L)
})

# ============================================================================
# Group 7: Numerical consistency with Python reference
# ============================================================================

test_that("get_cohorts matches Python sorted(valid_cohorts) output", {
  # Replicate Python test scenario: mixed cohorts with NT
  dt <- data.table::data.table(
    id = rep(1:6, each = 5),
    time = rep(1:5, 6),
    gvar = rep(c(4, 4, 6, 6, 0, NA), each = 5)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  # Python: sorted([4, 6]) = [4, 6]
  expect_equal(cohorts, c(4L, 6L))
})

test_that("is_never_treated matches Python identify_never_treated_units", {
  # Python default: never_treated_values = [0, np.inf], plus NaN
  # R additionally handles near-zero values
  g_vals <- c(3, 5, 0, Inf, NA, NaN)
  r_result <- is_never_treated(g_vals)
  # Python would return: [F, F, T, T, T, T] for [3, 5, 0, inf, nan, nan]
  expect_equal(r_result, c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
})

test_that("check_never_treated cohort_sizes match Python groupby counts", {
  # Replicate Python scenario
  dt <- data.table::data.table(
    id = rep(1:8, each = 4),
    time = rep(1:4, 8),
    gvar = rep(c(3, 3, 3, 5, 5, NA, NA, 0), each = 4)
  )
  result <- check_never_treated(dt, "id", "gvar")
  expect_equal(result$n_treated, 5L)  # 3 in cohort 3, 2 in cohort 5
  expect_equal(result$n_nt, 3L)       # 2 NA + 1 zero
  expect_equal(as.integer(result$cohort_sizes[["3"]]), 3L)
  expect_equal(as.integer(result$cohort_sizes[["5"]]), 2L)
})

# ============================================================================
# Group 8: End-to-end integration tests
# ============================================================================

test_that("end-to-end: get_cohorts + check_never_treated consistency", {
  dt <- make_staggered_panel(c(g3 = 4, g5 = 3, g7 = 2), n_nt = 5)
  cohorts <- get_cohorts(dt, "gvar", "id")
  nt_info <- check_never_treated(dt, "id", "gvar")

  # Cohorts from get_cohorts should match cohort_sizes keys
  expect_equal(cohorts, as.integer(sort(as.numeric(
    names(nt_info$cohort_sizes)))))

  # Total units = n_nt + n_treated
  n_total_units <- length(unique(dt$id))
  expect_equal(nt_info$n_nt + nt_info$n_treated, n_total_units)

  # cohort_sizes sum = n_treated
  expect_equal(sum(nt_info$cohort_sizes), nt_info$n_treated)
})

test_that("end-to-end: AET scenario (no NT units)", {
  # All-Eventually-Treated: no NT units (lw2025 Section 4.3)
  dt <- data.table::data.table(
    id = rep(1:6, each = 5),
    time = rep(1:5, 6),
    gvar = rep(c(3, 3, 4, 4, 5, 5), each = 5)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, c(3L, 4L, 5L))

  nt_info <- check_never_treated(dt, "id", "gvar")
  expect_false(nt_info$has_nt)
  expect_equal(nt_info$n_nt, 0L)
  expect_equal(nt_info$n_treated, 6L)

  # Consistency: cohort sets match
  expect_equal(cohorts, as.integer(sort(as.numeric(
    names(nt_info$cohort_sizes)))))
  # n_nt + n_treated = total units
  expect_equal(nt_info$n_nt + nt_info$n_treated,
               length(unique(dt$id)))
})

test_that("end-to-end: realistic staggered panel scenario", {
  # Simulate a realistic dataset: 100 units, 10 periods,
  # 3 cohorts (g=4,6,8), 20 NT units
  set.seed(42)
  n_units <- 100
  n_periods <- 10
  # Assign cohorts: 30 to g=4, 25 to g=6, 25 to g=8, 20 NT
  g_assign <- c(rep(4, 30), rep(6, 25), rep(8, 25), rep(NA, 20))
  dt <- data.table::CJ(id = seq_len(n_units),
                        time = seq_len(n_periods))
  dt[, gvar := g_assign[id]]

  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, c(4L, 6L, 8L))

  nt_info <- check_never_treated(dt, "id", "gvar")
  expect_true(nt_info$has_nt)
  expect_equal(nt_info$n_nt, 20L)
  expect_equal(nt_info$n_treated, 80L)
  expect_equal(as.integer(nt_info$cohort_sizes[["4"]]), 30L)
  expect_equal(as.integer(nt_info$cohort_sizes[["6"]]), 25L)
  expect_equal(as.integer(nt_info$cohort_sizes[["8"]]), 25L)
})

# ============================================================================
# Group 9: Independent numerical verification
# ============================================================================

test_that("get_cohorts output verified by independent calculation", {
  # Known data: gvar values = {3, 3, 5, 5, 7, NA, 0, Inf}
  # Expected: unique non-NT sorted = {3, 5, 7}
  dt <- data.table::data.table(
    id = rep(1:8, each = 2),
    time = rep(1:2, 8),
    gvar = rep(c(3, 3, 5, 5, 7, NA, 0, Inf), each = 2)
  )
  cohorts <- get_cohorts(dt, "gvar", "id")
  # Verified: sort(unique(setdiff({3,3,5,5,7,NA,0,Inf},
  #   {NA,NaN,0,Inf,-Inf}))) = {3, 5, 7}
  expect_equal(cohorts, c(3L, 5L, 7L))
})

test_that("check_never_treated counts verified by independent calculation", {
  dt <- data.table::data.table(
    id = rep(1:8, each = 2),
    time = rep(1:2, 8),
    gvar = rep(c(3, 3, 5, 5, 7, NA, 0, Inf), each = 2)
  )
  result <- check_never_treated(dt, "id", "gvar")
  # Verified:
  # NT units: id=6(NA), id=7(0), id=8(Inf) -> n_nt = 3
  # Treated: id=1,2(g=3), id=3,4(g=5), id=5(g=7) -> n_treated = 5
  # cohort_sizes: 3->2, 5->2, 7->1
  expect_equal(result$n_nt, 3L)
  expect_equal(result$n_treated, 5L)
  expect_equal(as.integer(result$cohort_sizes[["3"]]), 2L)
  expect_equal(as.integer(result$cohort_sizes[["5"]]), 2L)
  expect_equal(as.integer(result$cohort_sizes[["7"]]), 1L)
})

# ============================================================================
# Group 10: FATAL-001 Strict Inequality Tests (Story E4-02)
# ============================================================================

test_that("FATAL-001: G_i == r is excluded from not_yet_treated", {
  dt <- data.table::data.table(
    id = 1:4,
    time = rep(5L, 4),
    gvar = c(3, 5, 7, NA)
  )
  mask <- get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated")
  # gvar=3: focal cohort (g==3, and 3 < 5 so not NYT) -> FALSE

  # gvar=5: G_i == r, STRICT inequality excludes this -> FALSE
  # gvar=7: G_i > r -> TRUE
  # gvar=NA: never-treated -> TRUE
  expect_equal(mask, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

test_that("FATAL-001: G_i > r is included in not_yet_treated", {
  dt <- data.table::data.table(
    id = 1:4,
    time = rep(5L, 4),
    gvar = c(3, 6, 8, NA)
  )
  mask <- get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated")
  # gvar=3: focal cohort (g==3, and 3 < 5 so not NYT) -> FALSE
  # gvar=6: G_i > r -> TRUE
  # gvar=8: G_i > r -> TRUE
  # gvar=NA: never-treated -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 3L)
})

test_that("FATAL-001: strict vs non-strict inequality difference", {
  dt <- data.table::data.table(
    id = 1:8,
    time = rep(5L, 8),
    gvar = c(3, 4, 5, 6, 7, NA, 0, Inf)
  )
  mask <- get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated")
  # gvar=3: focal cohort -> FALSE
  # gvar=4: G_i=4 < 5, not NYT, not NT -> FALSE
  # gvar=5: G_i == r, STRICT excludes -> FALSE
  # gvar=6: G_i > r -> TRUE
  # gvar=7: G_i > r -> TRUE
  # gvar=NA: NT -> TRUE
  # gvar=0: NT -> TRUE
  # gvar=Inf: NT -> TRUE
  # With strict >:  2 NYT + 3 NT = 5 controls
  # With non-strict >=: 3 NYT + 3 NT = 6 controls (WRONG!)
  expect_equal(mask, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 5L)
  # Explicitly verify gvar=5 unit (index 3) is excluded
  expect_false(mask[3])
})

test_that("FATAL-001: r=g boundary excludes focal cohort at treatment onset", {
  dt <- data.table::data.table(
    id = 1:4,
    time = rep(5L, 4),
    gvar = c(5, 5, 7, NA)
  )
  mask <- get_valid_controls(dt, "gvar", 5, 5, "not_yet_treated")
  # gvar=5: focal cohort AND G_i == r -> FALSE
  # gvar=5: same -> FALSE
  # gvar=7: G_i > r -> TRUE
  # gvar=NA: NT -> TRUE
  expect_equal(mask, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

# ============================================================================
# Group 11: Three Strategy Unit Tests (Story E4-02)
# ============================================================================

# Helper: create single-period cross-section for control group tests
make_control_test_panel <- function(gvar_values) {
  data.table::data.table(
    id = seq_along(gvar_values),
    time = rep(1L, length(gvar_values)),
    gvar = gvar_values
  )
}

# --- not_yet_treated strategy (5 tests) ---

test_that("NYT-01: basic not_yet_treated functionality", {
  dt <- make_control_test_panel(c(3, 5, 7, 9, NA, 0))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  # gvar=3: focal, 3>4=F, not NT -> FALSE
  # gvar=5: 5>4=T -> TRUE
  # gvar=7: 7>4=T -> TRUE
  # gvar=9: 9>4=T -> TRUE
  # gvar=NA: NT -> TRUE
  # gvar=0: NT -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 5L)
})

test_that("NYT-02: control group shrinks as r increases", {
  dt <- make_control_test_panel(c(3, 5, 7, 9, NA))
  mask_r4 <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  mask_r6 <- get_valid_controls(dt, "gvar", 3, 6, "not_yet_treated")
  mask_r8 <- get_valid_controls(dt, "gvar", 3, 8, "not_yet_treated")
  # r=4: NYT={5,7,9} + NT={NA} = 4
  # r=6: NYT={7,9} + NT={NA} = 3
  # r=8: NYT={9} + NT={NA} = 2
  expect_equal(sum(mask_r4), 4L)
  expect_equal(sum(mask_r6), 3L)
  expect_equal(sum(mask_r8), 2L)
})

test_that("NYT-03: five NT marker values all included", {
  dt <- make_control_test_panel(c(3, NA, NaN, 0, Inf))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  # gvar=3: 3>4=F, not NT -> FALSE
  # NA, NaN, 0, Inf: all NT -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

test_that("NYT-04: NA safety in not_yet_treated", {
  dt <- make_control_test_panel(c(3, NA, 7))
  mask <- get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated")
  # gvar=3: 3>5=F, not NT -> FALSE
  # gvar=NA: NA>5=NA coerced to F, but is_never_treated(NA)=T -> TRUE
  # gvar=7: 7>5=T -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

test_that("NYT-05: focal cohort included when G_i > r (pre-treatment)", {
  dt <- make_control_test_panel(c(5, 7, NA))
  mask <- get_valid_controls(dt, "gvar", 5, 3, "not_yet_treated")
  # gvar=5: focal cohort g=5, but 5>3=T -> TRUE (pre-treatment period)
  # gvar=7: 7>3=T -> TRUE
  # gvar=NA: NT -> TRUE
  expect_equal(mask, c(TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 3L)
})

# --- never_treated strategy (3 tests) ---

test_that("NT-06: only NT units included in never_treated", {
  dt <- make_control_test_panel(c(3, 5, 7, NA, 0))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "never_treated")
  # Only NA and 0 are NT
  expect_equal(mask, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

test_that("NT-07: never_treated does not depend on r", {
  dt <- make_control_test_panel(c(3, 5, NA, Inf))
  mask_r1   <- get_valid_controls(dt, "gvar", 3, 1,   "never_treated")
  mask_r100 <- get_valid_controls(dt, "gvar", 3, 100, "never_treated")
  # Both should be identical: only NT units (NA, Inf)
  expect_equal(mask_r1,   c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(mask_r100, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(sum(mask_r1), 2L)
  expect_equal(sum(mask_r100), 2L)
})

test_that("NT-08: multiple NT markers all recognized in never_treated", {
  dt <- make_control_test_panel(c(3, NA, NaN, 0, Inf))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "never_treated")
  # gvar=3: not NT -> FALSE
  # NA, NaN, 0, Inf: all NT -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

# --- all_others strategy (4 tests) ---

test_that("AO-09: all_others includes already-treated units", {
  dt <- make_control_test_panel(c(3, 4, 5, NA))
  # Suppress the all_others bias warning
  suppressWarnings(
    mask <- get_valid_controls(dt, "gvar", 5, 5, "all_others")
  )
  # gvar=3: G!=g -> TRUE
  # gvar=4: G!=g -> TRUE
  # gvar=5: focal -> FALSE
  # gvar=NA: NT -> TRUE
  expect_equal(mask, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(sum(mask), 3L)
})

test_that("AO-10: all_others handles NT NA values correctly", {
  dt <- make_control_test_panel(c(5, NA, 0, Inf, 3))
  suppressWarnings(
    mask <- get_valid_controls(dt, "gvar", 5, 4, "all_others")
  )
  # gvar=5: focal -> FALSE
  # gvar=NA: NT -> TRUE
  # gvar=0: NT -> TRUE
  # gvar=Inf: NT -> TRUE
  # gvar=3: G!=g -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

test_that("AO-11: all_others excludes only focal cohort", {
  dt <- make_control_test_panel(c(3, 3, 5, 7, NA))
  suppressWarnings(
    mask <- get_valid_controls(dt, "gvar", 3, 5, "all_others")
  )
  # gvar=3: focal -> FALSE
  # gvar=3: focal -> FALSE
  # gvar=5: G!=g -> TRUE
  # gvar=7: G!=g -> TRUE
  # gvar=NA: NT -> TRUE
  expect_equal(mask, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 3L)
})

test_that("AO-12: all_others includes units with G < r (already treated)", {
  dt <- make_control_test_panel(c(5, 2, 3, NA))
  suppressWarnings(
    mask <- get_valid_controls(dt, "gvar", 5, 4, "all_others")
  )
  # gvar=5: focal -> FALSE
  # gvar=2: G<r, already treated, G!=g -> TRUE
  # gvar=3: G<r, already treated, G!=g -> TRUE
  # gvar=NA: NT -> TRUE
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 3L)
})

# ============================================================================
# Group 12: Auto Strategy & Auto-Switching Tests (Story E4-02)
# ============================================================================

# --- Auto strategy: 8 combinations (has_nt × aggregate) ---

test_that("AUTO-01: auto + has_nt + cohort -> never_treated", {
  result <- resolve_control_group("auto", "cohort", TRUE)
  expect_equal(result$resolved, "never_treated")
  expect_false(result$switched)
})

test_that("AUTO-02: auto + has_nt + overall -> never_treated", {
  result <- resolve_control_group("auto", "overall", TRUE)
  expect_equal(result$resolved, "never_treated")
  expect_false(result$switched)
})

test_that("AUTO-03: auto + has_nt + none -> not_yet_treated", {
  result <- resolve_control_group("auto", "none", TRUE)
  expect_equal(result$resolved, "not_yet_treated")
  expect_false(result$switched)
})

test_that("AUTO-04: auto + has_nt + event_time -> not_yet_treated", {
  result <- resolve_control_group("auto", "event_time", TRUE)
  expect_equal(result$resolved, "not_yet_treated")
  expect_false(result$switched)
})

test_that("AUTO-05: auto + AET + cohort -> error", {
  expect_error(
    resolve_control_group("auto", "cohort", FALSE),
    class = "lwdid_no_never_treated"
  )
})

test_that("AUTO-06: auto + AET + overall -> error", {
  expect_error(
    resolve_control_group("auto", "overall", FALSE),
    class = "lwdid_no_never_treated"
  )
})

test_that("AUTO-07: auto + AET + none -> not_yet_treated", {
  result <- resolve_control_group("auto", "none", FALSE)
  expect_equal(result$resolved, "not_yet_treated")
  expect_false(result$switched)
})

test_that("AUTO-08: auto + AET + event_time -> not_yet_treated", {
  result <- resolve_control_group("auto", "event_time", FALSE)
  expect_equal(result$resolved, "not_yet_treated")
  expect_false(result$switched)
})

# --- Auto-switching tests (6 tests) ---

test_that("SWITCH-09: nyt + cohort -> auto-switch to never_treated", {
  expect_message(
    result <- resolve_control_group("not_yet_treated", "cohort", TRUE),
    "Switching"
  )
  expect_equal(result$resolved, "never_treated")
  expect_true(result$switched)
})

test_that("SWITCH-10: nyt + overall -> auto-switch to never_treated", {
  expect_message(
    result <- resolve_control_group("not_yet_treated", "overall", TRUE),
    "Switching"
  )
  expect_equal(result$resolved, "never_treated")
  expect_true(result$switched)
})

test_that("SWITCH-11: ao + cohort -> auto-switch to never_treated", {
  expect_message(
    result <- resolve_control_group("all_others", "cohort", TRUE),
    "Switching"
  )
  expect_equal(result$resolved, "never_treated")
  expect_true(result$switched)
})

test_that("SWITCH-12: ao + overall -> auto-switch to never_treated", {
  expect_message(
    result <- resolve_control_group("all_others", "overall", TRUE),
    "Switching"
  )
  expect_equal(result$resolved, "never_treated")
  expect_true(result$switched)
})

test_that("SWITCH-13: nyt + none -> no switch", {
  result <- resolve_control_group("not_yet_treated", "none", TRUE)
  expect_equal(result$resolved, "not_yet_treated")
  expect_false(result$switched)
})

test_that("SWITCH-14: nt + cohort -> no switch (already correct)", {
  result <- resolve_control_group("never_treated", "cohort", TRUE)
  expect_equal(result$resolved, "never_treated")
  expect_false(result$switched)
})

# --- all_others warning tests (2 tests) ---

test_that("WARN-15: all_others + none -> lwdid_data warning", {
  expect_warning(
    result <- resolve_control_group("all_others", "none", TRUE),
    class = "lwdid_data"
  )
  expect_equal(result$resolved, "all_others")
  expect_false(result$switched)
})

test_that("WARN-16: all_others + event_time -> lwdid_data warning", {
  expect_warning(
    result <- resolve_control_group("all_others", "event_time", TRUE),
    class = "lwdid_data"
  )
  expect_equal(result$resolved, "all_others")
  expect_false(result$switched)
})

# ============================================================================
# Group 13: Error Handling & Boundary Condition Tests (Story E4-02)
# ============================================================================

test_that("ERR-01: empty nyt control group -> lwdid_no_control", {
  dt <- make_control_test_panel(c(3, 4))
  expect_error(
    get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated"),
    class = "lwdid_no_control"
  )
})

test_that("ERR-02: lwdid_no_control error has correct fields", {
  dt <- make_control_test_panel(c(3, 4))
  err <- tryCatch(
    get_valid_controls(dt, "gvar", 3, 5, "not_yet_treated"),
    lwdid_no_control = function(e) e
  )
  expect_equal(err$g, 3)
  expect_equal(err$r, 5)
  expect_equal(err$control_group, "not_yet_treated")
})

test_that("ERR-03: empty nt control group -> lwdid_no_control", {
  dt <- make_control_test_panel(c(3, 5, 7))
  expect_error(
    get_valid_controls(dt, "gvar", 3, 4, "never_treated"),
    class = "lwdid_no_control"
  )
})

test_that("BC-04: pre-treatment period r < g works", {
  dt <- make_control_test_panel(c(5, 7, NA))
  mask <- get_valid_controls(dt, "gvar", 5, 3, "not_yet_treated")
  expect_equal(mask, c(TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 3L)
})

test_that("BC-05: pre-treatment focal cohort included as control", {
  dt <- make_control_test_panel(c(5, 5, 7, NA))
  mask <- get_valid_controls(dt, "gvar", 5, 2, "not_yet_treated")
  # All units satisfy G_i > 2 or are NT
  expect_equal(mask, c(TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

test_that("BC-06: pre-treatment never_treated unaffected", {
  dt <- make_control_test_panel(c(5, NA, 0))
  mask <- get_valid_controls(dt, "gvar", 5, 2, "never_treated")
  expect_equal(mask, c(FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

test_that("OBS-07: mask length equals nrow(data)", {
  dt <- make_control_test_panel(c(3, 5, 7, NA, 0, Inf))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  expect_equal(length(mask), nrow(dt))
  expect_equal(length(mask), 6L)
})

test_that("OBS-08: mask can subset data directly", {
  dt <- make_control_test_panel(c(3, 5, 7, NA))
  mask <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  sub <- dt[mask]
  # gvar=5 (5>4), gvar=7 (7>4), gvar=NA (NT) → 3 rows
  expect_equal(nrow(sub), 3L)
  # Verify the subset contains expected gvar values
  expect_true(5 %in% sub$gvar)
  expect_true(7 %in% sub$gvar)
  expect_true(any(is.na(sub$gvar)))
})

test_that("OBS-09: observation-level mask consistent with unique logic", {
  # Panel with 2 rows per unit (multi-period), filter to single period
  dt_full <- data.table::data.table(
    id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    gvar = rep(c(3, 5, 7, NA), each = 2)
  )
  # Single-period cross-section
  dt_t1 <- dt_full[time == 1]
  mask <- get_valid_controls(dt_t1, "gvar", 3, 4, "not_yet_treated")
  expect_equal(length(mask), nrow(dt_t1))
  expect_equal(length(mask), 4L)
  # Same result as working with unique gvar values
  unique_g <- unique(dt_t1$gvar)
  unique_mask <- is_never_treated(unique_g) | (unique_g > 4 & !is.na(unique_g))
  unique_mask[is.na(unique_mask)] <- FALSE
  expect_equal(sum(mask), sum(unique_mask))
})

test_that("ERR-10: auto passed to get_valid_controls -> error", {
  dt <- make_control_test_panel(c(3, 5, NA))
  expect_error(
    get_valid_controls(dt, "gvar", 3, 4, "auto"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("ERR-11: invalid strategy string -> error", {
  dt <- make_control_test_panel(c(3, 5, NA))
  expect_error(
    get_valid_controls(dt, "gvar", 3, 4, "invalid_strategy"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("STRUCT-12: resolve returns list(resolved, switched)", {
  result <- resolve_control_group("auto", "none", TRUE)
  expect_true(is.list(result))
  expect_named(result, c("resolved", "switched"))
  expect_true(is.character(result$resolved))
  expect_true(is.logical(result$switched))
})

test_that("STRUCT-13: switched field is logical", {
  # Non-switched case
  r1 <- resolve_control_group("auto", "none", TRUE)
  expect_identical(r1$switched, FALSE)
  # Switched case
  expect_message(
    r2 <- resolve_control_group("not_yet_treated", "cohort", TRUE)
  )
  expect_identical(r2$switched, TRUE)
})

test_that("PURE-14: get_valid_controls does not modify input data", {
  dt <- make_control_test_panel(c(3, 5, 7, NA, 0))
  dt_copy <- data.table::copy(dt)
  mask <- get_valid_controls(dt, "gvar", 3, 4, "not_yet_treated")
  # Verify data unchanged
  expect_identical(dt, dt_copy)
})

# ============================================================================
# Group 14: Python Consistency & End-to-End Tests (Story E4-02)
# ============================================================================

test_that("PY-01: not_yet_treated matches Python reference", {
  dt <- make_control_test_panel(c(2004, 2006, 2008, NA, 0))
  mask <- get_valid_controls(dt, "gvar", 2004, 2005, "not_yet_treated")
  # Python: unit_gvar > period → {2006, 2008} + NT = 4
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

test_that("PY-02: never_treated matches Python NEVER_TREATED", {
  dt <- make_control_test_panel(c(2004, 2006, 2008, NA, 0))
  mask <- get_valid_controls(dt, "gvar", 2004, 2005, "never_treated")
  # Python: only NT units
  expect_equal(mask, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(sum(mask), 2L)
})

test_that("PY-03: Castle Law control group evolution", {
  dt <- make_control_test_panel(c(2004, 2006, 2008, NA))
  # Control group shrinks as r increases
  expect_equal(sum(get_valid_controls(dt, "gvar", 2004, 2004, "not_yet_treated")), 3L)
  expect_equal(sum(get_valid_controls(dt, "gvar", 2004, 2005, "not_yet_treated")), 3L)
  expect_equal(sum(get_valid_controls(dt, "gvar", 2004, 2006, "not_yet_treated")), 2L)
  expect_equal(sum(get_valid_controls(dt, "gvar", 2004, 2007, "not_yet_treated")), 2L)
  expect_equal(sum(get_valid_controls(dt, "gvar", 2004, 2008, "not_yet_treated")), 1L)
})

test_that("PY-04: all_others matches Python ALL_OTHERS", {
  dt <- make_control_test_panel(c(2004, 2006, 2008, NA, 0))
  suppressWarnings(
    mask <- get_valid_controls(dt, "gvar", 2006, 2005, "all_others")
  )
  # Python: unit_gvar != cohort → {2004, 2008} + NT = 4
  expect_equal(mask, c(TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(sum(mask), 4L)
})

test_that("E2E-05: resolve + get_valid_controls for cohort", {
  # Step 1: resolve
  result <- resolve_control_group("auto", "cohort", TRUE)
  expect_equal(result$resolved, "never_treated")

  # Step 2: apply resolved strategy
  dt <- make_control_test_panel(c(3, 5, 7, NA, 0, Inf))
  mask <- get_valid_controls(dt, "gvar", 3, 4, result$resolved)
  # never_treated: only NT = {NA, 0, Inf} = 3
  expect_equal(sum(mask), 3L)
  expect_equal(mask, c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that("E2E-06: AET + none -> not_yet_treated end-to-end", {
  # Step 1: resolve
  result <- resolve_control_group("auto", "none", FALSE)
  expect_equal(result$resolved, "not_yet_treated")

  # Step 2: apply on AET data (no NT units)
  dt <- make_control_test_panel(c(3, 5, 7, 9))
  mask <- get_valid_controls(dt, "gvar", 3, 4, result$resolved)
  # not_yet_treated: G>4 = {5, 7, 9} = 3 controls
  expect_equal(sum(mask), 3L)
  expect_equal(mask, c(FALSE, TRUE, TRUE, TRUE))
})

test_that("E2E-07: AET + cohort -> error at resolve stage", {
  # AET data cannot use cohort aggregation
  expect_error(
    resolve_control_group("auto", "cohort", FALSE),
    class = "lwdid_no_never_treated"
  )
  # The error prevents get_valid_controls from ever being called
})
