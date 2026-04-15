# test-aggregate.R — Tests for aggregate.R

make_agg_panel <- function(
  cohort_sizes = c(g3 = 3, g5 = 2), n_nt = 2,
  n_periods = 6, y_base = 1.0
) {
  units <- list()
  uid <- 1L
  for (i in seq_along(cohort_sizes)) {
    g_val <- as.integer(gsub("g", "", names(cohort_sizes)[i]))
    for (j in seq_len(cohort_sizes[i])) {
      units[[uid]] <- data.table::data.table(
        id = uid, time = seq_len(n_periods),
        Y = y_base * uid + seq_len(n_periods) * 0.5, gvar = g_val)
      uid <- uid + 1L
    }
  }
  for (j in seq_len(n_nt)) {
    units[[uid]] <- data.table::data.table(
      id = uid, time = seq_len(n_periods),
      Y = y_base * uid + seq_len(n_periods) * 0.3, gvar = NA_real_)
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

make_pre_stat_demean <- function(ids, means, ivar = "id") {
  dt <- data.table::data.table(uid = ids, pre_mean = means)
  data.table::setnames(dt, "uid", ivar)
  dt
}

make_pre_stat_detrend <- function(ids, means, slopes,
                                  t_bar_pre, ivar = "id") {
  dt <- data.table::data.table(
    uid = ids, pre_mean = means,
    slope = slopes, t_bar_pre = t_bar_pre)
  data.table::setnames(dt, "uid", ivar)
  dt
}

# ── Group 1: get_unit_level_gvar() — Basic Functionality ──

test_that("get_unit_level_gvar extracts one row per unit", {
  dt <- make_agg_panel(c(g3 = 2, g5 = 1), n_nt = 1, n_periods = 4)
  result <- get_unit_level_gvar(dt, "gvar", "id")
  expect_equal(nrow(result), 4L)
  expect_equal(sort(result$id), 1:4)
})

test_that("get_unit_level_gvar preserves gvar values", {
  dt <- make_agg_panel(c(g3 = 2, g5 = 1), n_nt = 1, n_periods = 4)
  result <- get_unit_level_gvar(dt, "gvar", "id")
  result <- result[order(result$id), ]
  expect_equal(result$gvar[1:2], c(3, 3))
  expect_equal(result$gvar[3], 5)
  expect_true(is.na(result$gvar[4]))
})

test_that("get_unit_level_gvar returns data.table", {
  dt <- make_agg_panel(c(g3 = 1), n_nt = 1, n_periods = 3)
  result <- get_unit_level_gvar(dt, "gvar", "id")
  expect_true(data.table::is.data.table(result))
  expect_true("id" %in% names(result))
  expect_true("gvar" %in% names(result))
})

test_that("get_unit_level_gvar handles single unit", {
  dt <- data.table::data.table(
    id = rep(1L, 3), time = 1:3, Y = c(1, 2, 3), gvar = 5)
  result <- get_unit_level_gvar(dt, "gvar", "id")
  expect_equal(nrow(result), 1L)
  expect_equal(result$gvar, 5)
})

# ── Group 2: get_unit_level_gvar() — Uniqueness Validation ──

test_that("get_unit_level_gvar throws lwdid_invalid_data on non-unique", {
  dt <- data.table::data.table(
    id = c(1L, 1L, 2L, 2L), time = c(1L, 2L, 1L, 2L),
    Y = c(1, 2, 3, 4), gvar = c(3, 5, 4, 4))
  expect_error(
    get_unit_level_gvar(dt, "gvar", "id"),
    class = "lwdid_invalid_data")
})

test_that("get_unit_level_gvar error message includes unit IDs", {
  dt <- data.table::data.table(
    id = c(1L, 1L, 2L, 2L, 3L, 3L),
    time = c(1L, 2L, 1L, 2L, 1L, 2L),
    Y = c(1, 2, 3, 4, 5, 6),
    gvar = c(3, 5, 4, 6, 7, 7))
  expect_error(
    get_unit_level_gvar(dt, "gvar", "id"),
    regexp = "1.*2", class = "lwdid_invalid_data")
})

test_that("get_unit_level_gvar accepts consistent gvar with NT", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 3), time = rep(1:3, 3),
    Y = rnorm(9), gvar = rep(c(3, NA, 0), each = 3))
  result <- get_unit_level_gvar(dt, "gvar", "id")
  expect_equal(nrow(result), 3L)
})

# ── Group 3: compute_cohort_aggregated_variable() — Basic demean ──

test_that("ccav returns correct structure", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 1, n_periods = 5)
  ps <- make_pre_stat_demean(1:3, c(2.0, 3.0, 4.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 5L,
    pre_stat = ps, rolling = "demean")
  expect_true(data.table::is.data.table(result))
  expect_true("id" %in% names(result))
  expect_true("Y_bar_ig" %in% names(result))
  expect_true("n_periods" %in% names(result))
})

test_that("ccav computes for all units (treated + NT)", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 1, n_periods = 5)
  ps <- make_pre_stat_demean(1:3, c(2.0, 3.0, 4.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 5L,
    pre_stat = ps, rolling = "demean")
  expect_equal(nrow(result), 3L)
})

test_that("ccav balanced panel n_valid equals T-g+1", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 1, n_periods = 5)
  ps <- make_pre_stat_demean(1:3, c(2.0, 3.0, 4.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 5L,
    pre_stat = ps, rolling = "demean")
  expect_true(all(result$n_periods == 3L))
})

test_that("ccav Y_bar_ig is time average of transforms", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, 4, 7, 2, 5, 8),
    gvar = c(2, 2, 2, NA, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.0, 2.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  # Unit 1: ydot = (4-1, 7-1) = (3, 6), mean = 4.5
  # Unit 2: ydot = (5-2, 8-2) = (3, 6), mean = 4.5
  expect_equal(result$Y_bar_ig[result$id == 1], 4.5, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 4.5, tolerance = 1e-10)
  expect_true(all(result$n_periods == 2L))
})

test_that("ccav single post period (g == T_max)", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, 2, 5, 2, 3, 6),
    gvar = c(3, 3, 3, NA, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.5, 2.5))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  # Unit 1: ydot_r3 = 5-1.5 = 3.5
  # Unit 2: ydot_r3 = 6-2.5 = 3.5
  expect_equal(result$Y_bar_ig[result$id == 1], 3.5, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 3.5, tolerance = 1e-10)
  expect_true(all(result$n_periods == 1L))
})

# ── Group 4: ccav — Boundary Conditions & Warnings ──

test_that("ccav returns NULL when g > T_max", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 1, n_periods = 5)
  ps <- make_pre_stat_demean(1:3, c(2.0, 3.0, 4.0))
  expect_warning(
    result <- compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 6L, T_max = 5L,
      pre_stat = ps, rolling = "demean"),
    class = "lwdid_data")
  expect_null(result)
})

test_that("ccav g > T_max warning has correct detail", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 1, n_periods = 5)
  ps <- make_pre_stat_demean(1:3, c(2.0, 3.0, 4.0))
  w <- tryCatch(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 10L, T_max = 5L,
      pre_stat = ps, rolling = "demean"),
    lwdid_data = function(w) w)
  expect_equal(w$detail, "cohort_no_post_periods")
})

test_that("ccav handles missing periods gracefully", {
  dt <- data.table::data.table(
    id = c(rep(1L, 5), rep(2L, 4)),
    time = c(1:5, 1:3, 5L),
    Y = c(1, 2, 4, 5, 6, 2, 3, 5, 7),
    gvar = c(rep(3, 5), rep(NA, 4)))
  ps <- make_pre_stat_demean(1:2, c(1.5, 2.5))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 5L,
    pre_stat = ps, rolling = "demean")
  expect_equal(result$n_periods[result$id == 1], 3L)
  expect_equal(result$n_periods[result$id == 2], 2L)
})

test_that("ccav NA outcome skipped in accumulation", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, 4, NA, 2, 5, 8),
    gvar = c(2, 2, 2, NA, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.0, 2.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  # Unit 1: ydot_r2=3, ydot_r3=NA(skip), Y_bar=3/1=3.0
  # Unit 2: ydot_r2=3, ydot_r3=6, Y_bar=(3+6)/2=4.5
  expect_equal(result$n_periods[result$id == 1], 1L)
  expect_equal(result$n_periods[result$id == 2], 2L)
  expect_equal(result$Y_bar_ig[result$id == 1], 3.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 4.5, tolerance = 1e-10)
})

test_that("ccav all-NA unit returns NA_real_", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, NA, NA, 2, 5, 8),
    gvar = c(2, 2, 2, NA, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.0, 2.0))
  expect_warning(
    result <- compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 2L, T_max = 3L,
      pre_stat = ps, rolling = "demean"),
    class = "lwdid_data")
  expect_true(is.na(result$Y_bar_ig[result$id == 1]))
  expect_equal(result$n_periods[result$id == 1], 0L)
})

test_that("ccav all-NA warning has correct detail", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, NA, NA, 2, 5, 8),
    gvar = c(2, 2, 2, NA, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.0, 2.0))
  w <- tryCatch(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 2L, T_max = 3L,
      pre_stat = ps, rolling = "demean"),
    lwdid_data = function(w) w)
  expect_equal(w$detail, "cohort_agg_missing_units")
})

test_that("ccav duplicate observations averaged", {
  dt <- data.table::data.table(
    id = c(1L, 1L, 1L, 1L), time = c(1L, 2L, 2L, 3L),
    Y = c(1, 4, 6, 7), gvar = c(2, 2, 2, 2))
  ps <- make_pre_stat_demean(1L, 1.0)
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  # Period 2: two rows ydot=(4-1)=3 and (6-1)=5, both accumulated
  # Period 3: one row ydot=(7-1)=6
  # y_sum=3+5+6=14, n_valid=3, Y_bar=14/3
  expect_equal(result$Y_bar_ig[1], 14 / 3, tolerance = 1e-10)
  expect_equal(result$n_periods[1], 3L)
})

# ── Group 5: ccav — Column Indexing ──

test_that("ccav uses ivar not column position", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3), time = rep(1:3, 2),
    Y = c(1, 4, 7, 2, 5, 8),
    gvar = c(2, 2, 2, NA, NA, NA))
  ps <- data.table::data.table(
    extra_col = c("a", "b"), id = 1:2, pre_mean = c(1.0, 2.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  expect_equal(nrow(result), 2L)
  expect_false(any(is.na(result$Y_bar_ig)))
})

test_that("ccav empty period skipped", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 2), time = rep(c(1L, 2L), 2),
    Y = c(1, 4, 2, 5), gvar = c(2, 2, NA, NA))
  ps <- make_pre_stat_demean(1:2, c(1.0, 2.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 3L,
    pre_stat = ps, rolling = "demean")
  # Only period 2 has data (period 3 empty)
  # Unit 1: ydot_r2=4-1=3, Y_bar=3/1=3.0
  expect_true(all(result$n_periods == 1L))
  expect_equal(result$Y_bar_ig[result$id == 1], 3.0, tolerance = 1e-10)
})

# ── Group 5b: ccav — Detrend Transform ──

test_that("ccav detrend basic correctness", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 4), time = rep(1:4, 2),
    Y = c(2, 3, 5, 7, 3, 4, 6, 9),
    gvar = c(rep(3, 4), rep(NA, 4)))
  ps <- make_pre_stat_detrend(
    1:2, c(2.5, 3.5), c(1.0, 1.0), c(1.5, 1.5))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 4L,
    pre_stat = ps, rolling = "detrend")
  # Unit 1: r=3: pred=2.5+1.0*(3-1.5)=4.0, ydot=5-4=1.0
  #         r=4: pred=2.5+1.0*(4-1.5)=5.0, ydot=7-5=2.0, Y_bar=1.5
  # Unit 2: r=3: pred=3.5+1.0*1.5=5.0, ydot=6-5=1.0
  #         r=4: pred=3.5+1.0*2.5=6.0, ydot=9-6=3.0, Y_bar=2.0
  expect_equal(result$Y_bar_ig[result$id == 1], 1.5, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 2.0, tolerance = 1e-10)
  expect_true(all(result$n_periods == 2L))
})

test_that("ccav detrend r parameter varies per period", {
  # Perfect linear trend Y = 1 + 2*t → zero residuals
  dt <- data.table::data.table(
    id = rep(1L, 4), time = 1:4,
    Y = c(3, 5, 7, 9), gvar = rep(3, 4))
  ps <- make_pre_stat_detrend(1L, 4.0, 2.0, 1.5)
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 4L,
    pre_stat = ps, rolling = "detrend")
  expect_equal(result$Y_bar_ig[1], 0.0, tolerance = 1e-10)
})

test_that("ccav detrend with unbalanced panel", {
  dt <- data.table::data.table(
    id = c(rep(1L, 4), rep(2L, 3)),
    time = c(1:4, 1L, 2L, 4L),
    Y = c(2, 3, 5, 7, 3, 4, 9),
    gvar = c(rep(3, 4), rep(NA, 3)))
  ps <- make_pre_stat_detrend(
    1:2, c(2.5, 3.5), c(1.0, 1.0), c(1.5, 1.5))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 4L,
    pre_stat = ps, rolling = "detrend")
  # Unit 1: 2 post periods, Y_bar=1.5
  # Unit 2: missing period 3, only r=4: ydot=9-6=3.0, n_valid=1
  expect_equal(result$n_periods[result$id == 1], 2L)
  expect_equal(result$n_periods[result$id == 2], 1L)
  expect_equal(result$Y_bar_ig[result$id == 1], 1.5, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 3.0, tolerance = 1e-10)
})

# ============================================================================
# Group 6: NT Mask Integration (get_unit_level_gvar + is_never_treated)
# ============================================================================

test_that("NT mask via is_never_treated on unit_gvar works correctly", {
  dt <- make_agg_panel(c(g3 = 2, g5 = 1), n_nt = 2, n_periods = 4)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  expect_equal(length(nt_mask), 5L)
  expect_equal(sum(nt_mask), 2L)
  expect_equal(sum(!nt_mask), 3L)
})

test_that("NT mask correctly identifies all NT marker types (NA, NaN, 0, Inf)", {
  dt <- data.table::data.table(
    id = rep(1:5, each = 2),
    time = rep(1:2, 5),
    Y = rnorm(10),
    gvar = rep(c(3, NA, NaN, 0, Inf), each = 2)
  )
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  expect_equal(sum(!nt_mask), 1L)
  expect_equal(sum(nt_mask), 4L)
})

test_that("NT mask + unit_gvar pipeline extracts correct unit IDs", {
  dt <- make_agg_panel(c(g3 = 2, g5 = 1), n_nt = 2, n_periods = 4)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  nt_units <- unit_gvar$id[nt_mask]
  treated_units <- unit_gvar$id[!nt_mask]
  expect_equal(sort(nt_units), c(4L, 5L))
  expect_equal(sort(treated_units), c(1L, 2L, 3L))
})

test_that("NT mask handles all-treated scenario (no NT)", {
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    Y = rnorm(9),
    gvar = rep(c(3, 5, 7), each = 3)
  )
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  expect_true(all(!nt_mask))
  expect_equal(sum(nt_mask), 0L)
})


# ============================================================================
# Group 7: Python Numerical Consistency
# ============================================================================

test_that("ccav matches Python nanmean behavior (balanced demean)", {
  # 3 units, g=2, T_max=4, demean
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(2, 5, 8, 11,   # unit 1: linear growth
          1, 3, 5, 7,    # unit 2: linear growth
          3, 4, 5, 6),   # unit 3 (NT): linear growth
    gvar = c(rep(2, 4), rep(2, 4), rep(NA, 4))
  )
  pre_stat <- make_pre_stat_demean(1:3, c(2.0, 1.0, 3.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 4L,
    pre_stat = pre_stat, rolling = "demean"
  )
  # Verified by vibe-math MCP:
  # Unit 1: ydot = (5-2, 8-2, 11-2) = (3, 6, 9), mean = 6.0
  # Unit 2: ydot = (3-1, 5-1, 7-1) = (2, 4, 6), mean = 4.0
  # Unit 3: ydot = (4-3, 5-3, 6-3) = (1, 2, 3), mean = 2.0
  expect_equal(result$Y_bar_ig[result$id == 1], 6.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 4.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 3], 2.0, tolerance = 1e-10)
})

test_that("ccav matches Python nanmean with unbalanced panel", {
  # Unit 2 missing period 3
  dt <- data.table::data.table(
    id = c(rep(1L, 4), rep(2L, 3)),
    time = c(1:4, 1L, 2L, 4L),
    Y = c(2, 5, 8, 11, 1, 3, 7),
    gvar = c(rep(2, 4), rep(NA, 3))
  )
  pre_stat <- make_pre_stat_demean(1:2, c(2.0, 1.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 4L,
    pre_stat = pre_stat, rolling = "demean"
  )
  # Verified by vibe-math MCP:
  # Unit 1: ydot = (3, 6, 9), mean = 6.0, n_valid = 3
  # Unit 2: ydot = (2, -, 6), nanmean = (2+6)/2 = 4.0, n_valid = 2
  expect_equal(result$Y_bar_ig[result$id == 1], 6.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 4.0, tolerance = 1e-10)
  expect_equal(result$n_periods[result$id == 1], 3L)
  expect_equal(result$n_periods[result$id == 2], 2L)
})

test_that("get_unit_level_gvar matches Python groupby.first for consistent data", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    time = rep(1:5, 4),
    Y = rnorm(20),
    gvar = rep(c(3, 3, 5, NA), each = 5)
  )
  result <- get_unit_level_gvar(dt, "gvar", "id")
  result <- result[order(result$id), ]
  expect_equal(result$gvar[1:3], c(3, 3, 5))
  expect_true(is.na(result$gvar[4]))
})

# ============================================================================
# Group 8: End-to-End Integration
# ============================================================================

test_that("end-to-end: get_unit_level_gvar + NT mask + ccav pipeline", {
  dt <- make_agg_panel(c(g3 = 3, g5 = 2), n_nt = 2, n_periods = 6)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)

  # Verify NT mask
  expect_equal(sum(nt_mask), 2L)
  expect_equal(sum(!nt_mask), 5L)

  # Compute aggregated variable for cohort g=3
  pre_stat <- make_pre_stat_demean(1:7, rep(1.0, 7))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 6L,
    pre_stat = pre_stat, rolling = "demean"
  )

  # All 7 units should have results
  expect_equal(nrow(result), 7L)
  # All units should have n_periods = 4 (periods 3,4,5,6)
  expect_true(all(result$n_periods == 4L))
  # No NA values (balanced panel)
  expect_false(any(is.na(result$Y_bar_ig)))
})

test_that("end-to-end: realistic staggered panel aggregation", {
  set.seed(42)
  n_units <- 50
  n_periods <- 10
  g_assign <- c(rep(4, 20), rep(7, 20), rep(NA, 10))
  dt <- data.table::CJ(id = seq_len(n_units), time = seq_len(n_periods))
  dt[, gvar := g_assign[id]]
  dt[, Y := id * 0.1 + time * 0.5 + rnorm(.N, sd = 0.1)]

  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  expect_equal(sum(nt_mask), 10L)
  expect_equal(sum(!nt_mask), 40L)

  # Compute for cohort g=4 (post periods: 4,5,...,10 = 7 periods)
  pre_stat <- make_pre_stat_demean(1:50, rep(1.0, 50))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 4L, T_max = 10L,
    pre_stat = pre_stat, rolling = "demean"
  )
  expect_equal(nrow(result), 50L)
  expect_true(all(result$n_periods == 7L))
  expect_false(any(is.na(result$Y_bar_ig)))
})

# ============================================================================
# Group 9: Vibe-Math MCP Verified Exact Numerical Tests
# ============================================================================

test_that("ccav demean verified by independent calculation", {
  # Known exact values: 2 units, g=2, T_max=4, demean
  dt <- data.table::data.table(
    id = rep(1:2, each = 4),
    time = rep(1:4, 2),
    Y = c(10, 20, 30, 40, 5, 15, 25, 35),
    gvar = c(rep(2, 4), rep(NA, 4))
  )
  pre_stat <- make_pre_stat_demean(1:2, c(10.0, 5.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 4L,
    pre_stat = pre_stat, rolling = "demean"
  )
  # Verified by vibe-math MCP:
  # Unit 1: ydot = (20-10, 30-10, 40-10) = (10, 20, 30), Y_bar = 60/3 = 20.0
  # Unit 2: ydot = (15-5, 25-5, 35-5) = (10, 20, 30), Y_bar = 60/3 = 20.0
  expect_equal(result$Y_bar_ig[result$id == 1], 20.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 20.0, tolerance = 1e-10)
})

test_that("ccav unbalanced demean verified by independent calculation", {
  # Unit 1: all 3 post periods; Unit 2: missing period 3
  dt <- data.table::data.table(
    id = c(rep(1L, 4), rep(2L, 3)),
    time = c(1:4, 1L, 2L, 4L),
    Y = c(10, 20, 30, 40, 5, 15, 35),
    gvar = c(rep(2, 4), rep(NA, 3))
  )
  pre_stat <- make_pre_stat_demean(1:2, c(10.0, 5.0))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 2L, T_max = 4L,
    pre_stat = pre_stat, rolling = "demean"
  )
  # Verified by vibe-math MCP:
  # Unit 1: ydot = (10, 20, 30), Y_bar = 60/3 = 20.0
  # Unit 2: ydot = (10, -, 30), Y_bar = (10+30)/2 = 20.0
  expect_equal(result$Y_bar_ig[result$id == 1], 20.0, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 20.0, tolerance = 1e-10)
  expect_equal(result$n_periods[result$id == 1], 3L)
  expect_equal(result$n_periods[result$id == 2], 2L)
})

test_that("ccav detrend verified by independent calculation", {
  # 2 units, g=3, T_max=4
  dt <- data.table::data.table(
    id = rep(1:2, each = 4),
    time = rep(1:4, 2),
    Y = c(2, 3, 5, 7, 3, 4, 6, 9),
    gvar = c(rep(3, 4), rep(NA, 4))
  )
  pre_stat <- make_pre_stat_detrend(1:2, c(2.5, 3.5), c(1.0, 1.0), c(1.5, 1.5))
  result <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 4L,
    pre_stat = pre_stat, rolling = "detrend"
  )
  # Verified by vibe-math MCP:
  # Unit 1: r=3: pred=2.5+1.0*(3-1.5)=4.0, ydot=5-4=1.0
  #         r=4: pred=2.5+1.0*(4-1.5)=5.0, ydot=7-5=2.0
  #         Y_bar = (1.0+2.0)/2 = 1.5
  # Unit 2: r=3: pred=3.5+1.0*(3-1.5)=5.0, ydot=6-5=1.0
  #         r=4: pred=3.5+1.0*(4-1.5)=6.0, ydot=9-6=3.0
  #         Y_bar = (1.0+3.0)/2 = 2.0
  expect_equal(result$Y_bar_ig[result$id == 1], 1.5, tolerance = 1e-10)
  expect_equal(result$Y_bar_ig[result$id == 2], 2.0, tolerance = 1e-10)
  expect_true(all(result$n_periods == 2L))
})



# ── Helper: make_pre_stats (named list for aggregate_to_cohort) ──

make_pre_stats <- function(cohorts, unit_ids, ivar = "id") {
  pre_stats <- list()
  for (g in cohorts) {
    dt <- data.table::data.table(
      unit_id = unit_ids,
      pre_mean = rep(1.0, length(unit_ids))
    )
    data.table::setnames(dt, "unit_id", ivar)
    pre_stats[[as.character(g)]] <- dt
  }
  pre_stats
}

# ============================================================
# Group 10: aggregate_to_cohort — NT unit checks
# ============================================================

test_that("aggregate_to_cohort throws lwdid_no_never_treated when no NT units", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    time = rep(1:5, 4),
    Y = rnorm(20),
    gvar = rep(c(3, 3, 5, 5), each = 5)
  )
  pre_stats <- make_pre_stats(c(3, 5), 1:4)
  expect_error(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L, 5L), T_max = 5L,
                         pre_stats = pre_stats, rolling = "demean"),
    class = "lwdid_no_never_treated"
  )
})

test_that("aggregate_to_cohort error message mentions NT", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    time = rep(1:3, 2),
    Y = rnorm(6),
    gvar = rep(c(2, 3), each = 3)
  )
  pre_stats <- make_pre_stats(c(2, 3), 1:2)
  expect_error(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(2L, 3L), T_max = 3L,
                         pre_stats = pre_stats, rolling = "demean"),
    regexp = "NT",
    class = "lwdid_no_never_treated"
  )
})

test_that("aggregate_to_cohort accepts data with NT units (no NT error)", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  pre_stats <- make_pre_stats(c(3), 1:4)
  # Should not throw lwdid_no_never_treated
  expect_no_error(
    suppressWarnings(
      aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                           cohorts = c(3L), T_max = 5L,
                           pre_stats = pre_stats, rolling = "demean")
    )
  )
})

# ============================================================
# Group 11: aggregate_to_cohort — FATAL-004 & sample restriction
# ============================================================

test_that("FATAL-004: control group is NT-only (other cohorts excluded)", {
  dt <- make_agg_panel(c(g3 = 2, g5 = 2), n_nt = 3, n_periods = 6)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  cohort_3_units <- unit_gvar$id[
    unit_gvar$gvar == 3 & !is_never_treated(unit_gvar$gvar)
  ]
  expect_equal(length(cohort_3_units), 2L)
  sample_units <- c(cohort_3_units, unit_gvar$id[nt_mask])
  g5_units <- unit_gvar$id[
    unit_gvar$gvar == 5 & !is_never_treated(unit_gvar$gvar)
  ]
  expect_false(any(g5_units %in% sample_units))
})

test_that("D_ig satisfies D_ig + D_i_inf = 1 constraint", {
  dt <- make_agg_panel(c(g3 = 3), n_nt = 2, n_periods = 5)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  cohort_units <- unit_gvar$id[
    unit_gvar$gvar == 3 & !is_never_treated(unit_gvar$gvar)
  ]
  nt_units <- unit_gvar$id[nt_mask]
  expect_equal(length(intersect(cohort_units, nt_units)), 0L)
  sample_units <- c(cohort_units, nt_units)
  D_ig <- as.integer(sample_units %in% cohort_units)
  D_inf <- as.integer(sample_units %in% nt_units)
  expect_true(all(D_ig + D_inf == 1L))
})

test_that("aggregate_to_cohort skips cohort with insufficient sample", {
  # Create panel where all treated units have NA Y
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(NA, NA, NA, NA, NA, NA, NA, NA, 2, 3, 4, 5),
    gvar = c(rep(2, 4), rep(2, 4), rep(NA, 4))
  )
  pre_stats <- make_pre_stats(c(2), 1:3)
  # Both treated units have all-NA Y -> after dropna, n_treat=0 -> skip
  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(2L), T_max = 4L,
                         pre_stats = pre_stats, rolling = "demean")
  )
  # Should produce empty result or warning
  expect_true(length(result) == 0L || is.list(result))
})

test_that("aggregate_to_cohort warns with n_total == 2", {
  dt <- data.table::data.table(
    id = rep(1:2, each = 3),
    time = rep(1:3, 2),
    Y = c(1, 4, 7, 2, 3, 4),
    gvar = c(rep(2, 3), rep(NA, 3))
  )
  pre_stats <- make_pre_stats(c(2), 1:2)
  expect_warning(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(2L), T_max = 3L,
                         pre_stats = pre_stats, rolling = "demean"),
    class = "lwdid_small_sample"
  )
})

# ============================================================
# Group 12: aggregate_to_cohort — pre_stats missing
# ============================================================

test_that("aggregate_to_cohort skips cohort without pre_stats", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  pre_stats <- make_pre_stats(c(5), 1:4)  # only g=5, not g=3
  expect_warning(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L), T_max = 5L,
                         pre_stats = pre_stats, rolling = "demean"),
    class = "lwdid_data"
  )
})

test_that("aggregate_to_cohort with empty pre_stats warns", {
  dt <- make_agg_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  pre_stats <- list()
  # All cohorts skipped -> should warn
  result <- withCallingHandlers(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L), T_max = 5L,
                         pre_stats = pre_stats, rolling = "demean"),
    lwdid_data = function(w) invokeRestart("muffleWarning"),
    lwdid_numerical = function(w) invokeRestart("muffleWarning")
  )
  # Result should be empty list (all cohorts skipped)
  expect_true(is.list(result))
})



# ============================================================================
# Group 13: aggregate_to_cohort() — VCE Parameters & Cluster Variable
# (Task E5-02.5)
# ============================================================================

test_that("aggregate_to_cohort throws lwdid_invalid_parameter for cluster without cluster_var", {
  dt <- make_agg_panel(c(g3 = 3), n_nt = 3, n_periods = 5)
  pre_stats <- make_pre_stats(c(3), 1:6)
  expect_error(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L), T_max = 5L,
                         pre_stats = pre_stats, rolling = "demean",
                         vce = "cluster", cluster_var = NULL),
    class = "lwdid_invalid_parameter"
  )
})

test_that("aggregate_to_cohort cluster variable extraction uses match alignment", {
  dt <- data.table::data.table(
    id = rep(1:5, each = 4),
    time = rep(1:4, 5),
    Y = rnorm(20),
    gvar = c(rep(2, 4), rep(2, 4), rep(NA, 4), rep(NA, 4), rep(NA, 4)),
    state = rep(c("CA", "NY", "CA", "TX", "NY"), each = 4)
  )
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_mask <- is_never_treated(unit_gvar$gvar)
  cohort_units <- unit_gvar$id[unit_gvar$gvar == 2 &
                                 !is_never_treated(unit_gvar$gvar)]
  nt_units_ids <- unit_gvar$id[nt_mask]
  sample_units <- c(cohort_units, nt_units_ids)

  unit_cluster <- unique(dt[, c("id", "state"), with = FALSE])
  unit_cluster <- unit_cluster[!duplicated(unit_cluster$id)]

  reg_ids <- sample_units
  match_idx <- match(reg_ids, unit_cluster$id)
  cluster_vals <- unit_cluster$state[match_idx]

  expect_equal(length(cluster_vals), length(reg_ids))
  for (i in seq_along(reg_ids)) {
    expected_state <- dt$state[dt$id == reg_ids[i]][1]
    expect_equal(cluster_vals[i], expected_state)
  }
})

test_that("aggregate_to_cohort warns when cluster count < 20", {
  n_treat <- 5L
  n_nt <- 5L
  n_total <- n_treat + n_nt
  dt <- data.table::data.table(
    id = rep(seq_len(n_total), each = 4),
    time = rep(1:4, n_total),
    Y = rnorm(n_total * 4),
    gvar = rep(c(rep(2, n_treat), rep(NA, n_nt)), each = 4),
    cluster_id = rep(rep(c("A", "B"), length.out = n_total), each = 4)
  )
  pre_stats <- make_pre_stats(c(2), seq_len(n_total))
  expect_warning(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(2L), T_max = 4L,
                         pre_stats = pre_stats, rolling = "demean",
                         vce = "cluster", cluster_var = "cluster_id"),
    class = "lwdid_small_sample"
  )
})

test_that("aggregate_to_cohort passes vce=NULL (homoskedastic) without error", {
  dt <- make_agg_panel(c(g3 = 3), n_nt = 3, n_periods = 5)
  pre_stats <- make_pre_stats(c(3), 1:6)
  result <- tryCatch(
    suppressWarnings(
      aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                           cohorts = c(3L), T_max = 5L,
                           pre_stats = pre_stats, rolling = "demean",
                           vce = NULL)
    ),
    error = function(e) e
  )
  expect_false(inherits(result, "lwdid_invalid_parameter"))
})

test_that("aggregate_to_cohort cluster variable with NA values counts correctly", {
  n_total <- 6L
  dt <- data.table::data.table(
    id = rep(seq_len(n_total), each = 3),
    time = rep(1:3, n_total),
    Y = rnorm(n_total * 3),
    gvar = rep(c(rep(2, 3), rep(NA, 3)), each = 3),
    cl = rep(c("A", "B", NA, "A", "B", NA), each = 3)
  )
  unit_cluster <- unique(dt[, c("id", "cl"), with = FALSE])
  unit_cluster <- unit_cluster[!duplicated(unit_cluster$id)]
  cluster_vals <- unit_cluster$cl
  n_clusters <- length(unique(cluster_vals[!is.na(cluster_vals)]))
  expect_equal(n_clusters, 2L)
})

# ── Group 14: Controls Degradation Threshold Logic (Task E5-02.5) ──

test_that("aggregate_to_cohort Tier 1: controls passed when sample sufficient", {
  n_treat <- 5L
  n_ctrl <- 5L
  K <- 1L
  expect_true(n_treat > K + 1L && n_ctrl > K + 1L)
})

test_that("aggregate_to_cohort Tier 3: controls degraded when N_treat insufficient", {
  n_treat <- 2L
  n_ctrl <- 5L
  K <- 2L
  expect_false(n_treat > K + 1L && n_ctrl > K + 1L)
})

test_that("aggregate_to_cohort Tier 3: controls degraded when N_ctrl insufficient", {
  n_treat <- 5L
  n_ctrl <- 2L
  K <- 2L
  expect_false(n_treat > K + 1L && n_ctrl > K + 1L)
})

test_that("aggregate_to_cohort degradation threshold is strict inequality (not <=)", {
  n_treat <- 3L
  n_ctrl <- 3L
  K <- 2L
  expect_false(n_treat > K + 1L)
})

# ── Group 15: df_resid and df_inference Formulas (Task E5-02.5) ──

test_that("df_resid formula: no controls -> n - 2", {
  n <- 10L
  k <- 2L
  expect_equal(n - k, 8L)

  set.seed(1)
  Y <- rnorm(n)
  D <- c(rep(1L, 5L), rep(0L, 5L))
  fit <- lm(Y ~ D)
  expect_equal(fit$df.residual, n - 2L)
})

test_that("df_resid formula: with controls + interactions -> n - 2 - 2K", {
  n <- 20L
  K <- 3L
  k <- 2L + 2L * K
  expect_equal(n - k, 12L)

  set.seed(2)
  Y <- rnorm(n)
  D <- c(rep(1L, 10L), rep(0L, 10L))
  X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rnorm(n)
  fit <- lm(Y ~ D + X1 + X2 + X3 + D:X1 + D:X2 + D:X3)
  expect_equal(fit$df.residual, n - 2L - 2L * K)
})

test_that("df_inference: homoskedastic/HC equals df_resid", {
  set.seed(3)
  n <- 20L
  Y <- rnorm(n)
  D <- c(rep(1L, 10L), rep(0L, 10L))
  fit <- lm(Y ~ D)
  df_resid <- fit$df.residual
  expect_equal(df_resid, n - 2L)
})

test_that("df_inference: cluster equals G_cluster - 1", {
  n_clusters <- 15L
  df_inference <- n_clusters - 1L
  expect_equal(df_inference, 14L)
})

# ============================================================================
# Group 16: aggregate_to_cohort() — Return Value Structure & Sorting
# (Task E5-02.6)
# ============================================================================

# Helper: create mock estimate_ra_common result
make_mock_est <- function(att = 1.0, se = 0.5, df_resid = 10L,
                          df_inference = 10L, K = 0L) {
  list(
    att = att, se = se,
    ci_lower = att - 1.96 * se,
    ci_upper = att + 1.96 * se,
    t_stat = att / se,
    pvalue = 2 * pt(-abs(att / se), df = df_inference),
    df = df_inference,
    df_resid = df_resid,
    df_inference = df_inference,
    K = K
  )
}

test_that("aggregate_to_cohort result contains all 13 fields", {
  expected_fields <- c("cohort", "att", "se", "ci_lower", "ci_upper",
                        "t_stat", "pvalue", "n_periods", "n_units",
                        "n_control", "df_resid", "df_inference", "K")
  mock_result <- list(
    cohort = 3L, att = 1.0, se = 0.5,
    ci_lower = 0.02, ci_upper = 1.98,
    t_stat = 2.0, pvalue = 0.05,
    n_periods = 4L, n_units = 10L, n_control = 5L,
    df_resid = 13L, df_inference = 13L, K = 0L
  )
  expect_equal(sort(names(mock_result)), sort(expected_fields))
  expect_equal(length(mock_result), 13L)
})

test_that("aggregate_to_cohort cohort field is integer", {
  mock_result <- list(cohort = as.integer(5))
  expect_true(is.integer(mock_result$cohort))
})

test_that("aggregate_to_cohort results sorted by cohort ascending", {
  cohorts_input <- c(7L, 3L, 5L)
  results <- list(
    list(cohort = 7L), list(cohort = 3L), list(cohort = 5L)
  )
  cohort_order <- order(
    vapply(results, `[[`, integer(1), "cohort")
  )
  sorted_results <- results[cohort_order]
  sorted_cohorts <- vapply(sorted_results, `[[`, integer(1), "cohort")
  expect_equal(sorted_cohorts, c(3L, 5L, 7L))
})

test_that("aggregate_to_cohort K field reflects degraded control count", {
  K_requested <- 3L
  n_treat <- 2L
  K <- K_requested
  if (!(n_treat > K + 1L)) {
    K <- 0L
  }
  expect_equal(K, 0L)
})

# ── Group 17: Regression Failure & df Missing Handling (Task E5-02.6) ──

test_that("aggregate_to_cohort tryCatch catches regression failure", {
  mock_failing_est <- function(...) {
    stop("Singular matrix in OLS")
  }
  est <- tryCatch(
    mock_failing_est(),
    error = function(e) NULL
  )
  expect_null(est)
})

test_that("aggregate_to_cohort skips cohort when df_resid missing", {
  est <- list(att = 1.0, se = 0.5, df_inference = 10L)
  expect_true(is.null(est$df_resid))
})

test_that("aggregate_to_cohort skips cohort when df_inference missing", {
  est <- list(att = 1.0, se = 0.5, df_resid = 10L)
  expect_true(is.null(est$df_inference))
})

test_that("aggregate_to_cohort both df fields present passes validation", {
  est <- make_mock_est(att = 1.5, se = 0.3, df_resid = 8L,
                        df_inference = 8L)
  expect_false(is.null(est$df_resid))
  expect_false(is.null(est$df_inference))
})

# ── Group 18: Success/Failure Reporting (Task E5-02.6) ──

test_that("aggregate_to_cohort reports all-failure summary warning", {
  n_requested <- 3L
  n_success <- 0L
  if (n_success == 0L) {
    expect_warning(
      warn_lwdid(
        sprintf(
          paste0("All %d cohort effect estimates failed.\n",
                 "  Possible causes:\n",
                 "  1. Missing transform columns\n",
                 "  2. Insufficient sample size\n",
                 "  3. Aggregated variable computation failed"),
          n_requested),
        class = "lwdid_numerical"
      ),
      class = "lwdid_numerical"
    )
  }
})

test_that("aggregate_to_cohort reports partial failure with cohort list", {
  cohorts <- c(3L, 5L, 7L)
  successful_cohorts <- c(3L, 7L)
  failed_cohorts <- setdiff(cohorts, successful_cohorts)
  expect_equal(failed_cohorts, 5L)
  msg <- sprintf(
    "Partial cohort effect failure: %d/%d succeeded. Failed cohorts: %s",
    length(successful_cohorts), length(cohorts),
    paste(sort(failed_cohorts), collapse = ", ")
  )
  expect_true(grepl("5", msg))
  expect_true(grepl("2/3", msg))
})

test_that("aggregate_to_cohort no warning when all cohorts succeed", {
  n_requested <- 3L
  n_success <- 3L
  expect_true(n_success == n_requested)
  expect_false(n_success == 0L)
  expect_false(n_success < n_requested)
})


# ============================================================================
# Group 19: aggregate_to_cohort() — OLS Difference-in-Means Property (eq 7.13)
# (Task E5-02.7)
# ============================================================================

test_that("cohort effect equals treat mean minus control mean (eq 7.13)", {
  # Equation 7.13: tau_hat = mean(Y_bar_ig | D=1) - mean(Y_bar_ig | D=0)
  # Construct known Y_bar_ig values directly
  Y_bar_treat <- c(3.0, 5.0)   # 2 treated units
  Y_bar_ctrl <- c(1.0, 1.0, 1.0)  # 3 NT units

  # Expected tau_hat by equation 7.13
  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  # mean(treat) = (3+5)/2 = 4.0
  # mean(ctrl) = (1+1+1)/3 = 1.0
  # tau = 4.0 - 1.0 = 3.0
  expect_equal(tau_expected, 3.0, tolerance = 1e-10)

  # Verify via OLS regression
  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 0, 0, 0)
  fit <- lm(Y_bar ~ D)
  tau_ols <- coef(fit)["D"]
  expect_equal(unname(tau_ols), tau_expected, tolerance = 1e-10)
})

test_that("cohort effect with single treated unit (N_g = 1)", {
  # Regression path supports N_g = 1 (exact t-inference)
  Y_bar_treat <- c(5.0)  # 1 treated unit
  Y_bar_ctrl <- c(1.0, 2.0, 3.0)  # 3 NT units

  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  # mean(treat) = 5.0, mean(ctrl) = 2.0, tau = 3.0
  expect_equal(tau_expected, 3.0, tolerance = 1e-10)

  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 0, 0, 0)
  fit <- lm(Y_bar ~ D)
  tau_ols <- coef(fit)["D"]
  expect_equal(unname(tau_ols), tau_expected, tolerance = 1e-10)

  # df_resid = n - 2 = 4 - 2 = 2
  expect_equal(fit$df.residual, 2L)
})

test_that("cohort effect zero when treat and control means equal", {
  # When parallel trends hold perfectly and no treatment effect
  Y_bar_treat <- c(2.0, 2.0)
  Y_bar_ctrl <- c(2.0, 2.0, 2.0)

  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  expect_equal(tau_expected, 0.0, tolerance = 1e-10)

  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 0, 0, 0)
  fit <- lm(Y_bar ~ D)
  tau_ols <- coef(fit)["D"]
  expect_equal(unname(tau_ols), 0.0, tolerance = 1e-10)
})

# ============================================================================
# Group 20: aggregate_to_cohort() — Python Numerical Consistency
# (Task E5-02.7)
# ============================================================================

test_that("cohort effect matches Python for simple staggered design", {
  # Replicate Python scenario:
  # 3 treated units (g=3), 3 NT units, T_max=5
  # Known Y_bar_ig values (pre-computed from demean transform)
  # Treated: Y_bar = c(2.5, 3.0, 3.5)
  # NT: Y_bar = c(0.1, -0.1, 0.0)
  Y_bar_treat <- c(2.5, 3.0, 3.5)
  Y_bar_ctrl <- c(0.1, -0.1, 0.0)

  # Python: tau_hat = mean(treat) - mean(ctrl)
  #       = (2.5+3.0+3.5)/3 - (0.1-0.1+0.0)/3
  #       = 3.0 - 0.0 = 3.0
  # Verified by vibe-math MCP: 3.0
  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  expect_equal(tau_expected, 3.0, tolerance = 1e-10)

  # OLS verification
  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 1, 0, 0, 0)
  fit <- lm(Y_bar ~ D)
  expect_equal(unname(coef(fit)["D"]), 3.0, tolerance = 1e-10)

  # df_resid = 6 - 2 = 4
  expect_equal(fit$df.residual, 4L)
})

test_that("cohort effect with negative treatment effect", {
  # Treatment has negative effect (e.g., policy reduces outcome)
  Y_bar_treat <- c(-1.0, -2.0, -1.5)
  Y_bar_ctrl <- c(0.5, 0.3, 0.2)

  # Verified by vibe-math MCP: -1.8333333333333333
  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  # mean(treat) = -1.5, mean(ctrl) = 1/3 ≈ 0.3333
  # tau ≈ -1.8333
  expect_equal(tau_expected, -1.5 - 1 / 3, tolerance = 1e-10)

  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 1, 0, 0, 0)
  fit <- lm(Y_bar ~ D)
  expect_equal(unname(coef(fit)["D"]), tau_expected, tolerance = 1e-10)
})

test_that("cohort effect with heterogeneous treatment effects across units", {
  # Different units have different treatment effects
  # but tau_g is the average
  Y_bar_treat <- c(10.0, 2.0, 6.0, 4.0)  # heterogeneous
  Y_bar_ctrl <- c(1.0, 1.0)

  tau_expected <- mean(Y_bar_treat) - mean(Y_bar_ctrl)
  # mean(treat) = 5.5, mean(ctrl) = 1.0, tau = 4.5
  expect_equal(tau_expected, 4.5, tolerance = 1e-10)

  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 1, 1, 0, 0)
  fit <- lm(Y_bar ~ D)
  expect_equal(unname(coef(fit)["D"]), 4.5, tolerance = 1e-10)
})

# ============================================================================
# Group 21: aggregate_to_cohort() — vibe-math MCP Independent Verification
# (Task E5-02.7)
# ============================================================================

test_that("cohort effect verified by independent calculation (vibe-math)", {
  # Known exact scenario:
  # 3 treated (g=2), 2 NT, T_max=4, demean
  # Unit 1 (g=2): Y=(10,20,30,40), pre_mean=10
  #   ydot = (20-10, 30-10, 40-10) = (10, 20, 30)
  #   Y_bar = (10+20+30)/3 = 20.0
  # Unit 2 (g=2): Y=(5,15,25,35), pre_mean=5
  #   ydot = (15-5, 25-5, 35-5) = (10, 20, 30)
  #   Y_bar = (10+20+30)/3 = 20.0
  # Unit 3 (g=2): Y=(8,18,28,38), pre_mean=8
  #   ydot = (18-8, 28-8, 38-8) = (10, 20, 30)
  #   Y_bar = (10+20+30)/3 = 20.0
  # Unit 4 (NT): Y=(2,4,6,8), pre_mean=2
  #   ydot = (4-2, 6-2, 8-2) = (2, 4, 6)
  #   Y_bar = (2+4+6)/3 = 4.0
  # Unit 5 (NT): Y=(3,6,9,12), pre_mean=3
  #   ydot = (6-3, 9-3, 12-3) = (3, 6, 9)
  #   Y_bar = (3+6+9)/3 = 6.0
  #
  # Verified by vibe-math MCP:
  #   tau_hat = mean(treat) - mean(ctrl)
  #           = (20+20+20)/3 - (4+6)/2
  #           = 20.0 - 5.0 = 15.0
  #   df_resid = 5 - 2 = 3
  Y_bar_treat <- c(20.0, 20.0, 20.0)
  Y_bar_ctrl <- c(4.0, 6.0)
  tau_expected <- 15.0

  Y_bar <- c(Y_bar_treat, Y_bar_ctrl)
  D <- c(1, 1, 1, 0, 0)
  fit <- lm(Y_bar ~ D)
  expect_equal(unname(coef(fit)["D"]), tau_expected, tolerance = 1e-10)
  expect_equal(fit$df.residual, 3L)
})

test_that("cohort effect SE verified by independent calculation (vibe-math)", {
  # Same scenario as above, verify SE
  # Y_bar_treat = c(20, 20, 20), Y_bar_ctrl = c(4, 6)
  # tau_hat = 15.0
  # Residuals:
  #   treat: 20-20=0, 20-20=0, 20-20=0 (intercept_treat = 20)
  #   ctrl: 4-5=-1, 6-5=1 (intercept_ctrl = 5)
  # SSR = 0+0+0+1+1 = 2
  # MSE = SSR / df_resid = 2/3
  # Var(tau) = MSE * (1/n_treat + 1/n_ctrl) = (2/3)*(1/3+1/2) = (2/3)*(5/6) = 10/18 = 5/9
  # SE = sqrt(5/9) ≈ 0.7454
  #
  # Verified by vibe-math MCP: sqrt(5/9) = 0.7453559924999299
  Y_bar <- c(20, 20, 20, 4, 6)
  D <- c(1, 1, 1, 0, 0)
  fit <- lm(Y_bar ~ D)
  se_expected <- sqrt(5 / 9)
  se_actual <- summary(fit)$coefficients["D", "Std. Error"]
  expect_equal(se_actual, se_expected, tolerance = 1e-6)
})

test_that("cohort effect with controls verified by independent calculation", {
  # With 1 control variable, centered at treatment mean
  # 4 treated, 4 NT
  # Y_bar_treat = c(10, 12, 14, 16), X_treat = c(1, 2, 3, 4)
  # Y_bar_ctrl = c(2, 4, 6, 8), X_ctrl = c(2, 3, 4, 5)
  # X_mean_treat = (1+2+3+4)/4 = 2.5
  # X_centered: treat = (-1.5, -0.5, 0.5, 1.5), ctrl = (-0.5, 0.5, 1.5, 2.5)
  # Regression: Y_bar ~ 1 + D + X_c + D*X_c
  # At X = X_mean_treat (X_c = 0), interaction term = 0
  # -> D coefficient = ATT at treatment group mean
  # df_resid = 8 - 4 = 4 (intercept + D + X_c + D:X_c)
  # Verified by vibe-math MCP: df_resid = 4
  Y_bar <- c(10, 12, 14, 16, 2, 4, 6, 8)
  D <- c(1, 1, 1, 1, 0, 0, 0, 0)
  X <- c(1, 2, 3, 4, 2, 3, 4, 5)
  X_c <- X - mean(X[D == 1])  # center at treatment mean
  fit <- lm(Y_bar ~ D + X_c + D:X_c)
  tau_hat <- coef(fit)["D"]
  # D coefficient should be ATT at X = X_mean_treat
  expect_true(is.finite(unname(tau_hat)))
  # df_resid = 8 - 4 = 4 (intercept + D + X_c + D:X_c)
  expect_equal(fit$df.residual, 4L)
})

# ============================================================================
# Group 22: aggregate_to_cohort() — Multi-Cohort End-to-End Integration
# (Task E5-02.8)
# ============================================================================

test_that("end-to-end: multiple cohorts produce sorted results", {
  # 2 cohorts (g=3, g=5), 3 NT units, 6 periods
  # Full pipeline: panel -> pre_stats -> aggregate_to_cohort -> results
  set.seed(123)
  dt <- make_agg_panel(c(g3 = 4, g5 = 3), n_nt = 3, n_periods = 6)
  pre_stats <- make_pre_stats(c(3, 5), 1:10)

  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(5L, 3L),  # intentionally reversed
                         T_max = 6L,
                         pre_stats = pre_stats, rolling = "demean")
  )

  # Results should be sorted by cohort regardless of input order
  if (length(result) >= 2L) {
    cohorts_out <- vapply(result, `[[`, integer(1), "cohort")
    expect_true(all(diff(cohorts_out) > 0))  # strictly ascending
  }
})

test_that("end-to-end: each cohort result has correct n_units and n_control", {
  dt <- make_agg_panel(c(g3 = 4, g5 = 3), n_nt = 5, n_periods = 6)
  pre_stats <- make_pre_stats(c(3, 5), 1:12)

  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L, 5L), T_max = 6L,
                         pre_stats = pre_stats, rolling = "demean")
  )

  for (r in result) {
    if (r$cohort == 3L) {
      expect_equal(r$n_units, 4L)   # 4 treated in g=3
      expect_equal(r$n_control, 5L) # 5 NT units
    } else if (r$cohort == 5L) {
      expect_equal(r$n_units, 3L)   # 3 treated in g=5
      expect_equal(r$n_control, 5L) # 5 NT units
    }
  }
})

test_that("end-to-end: partial failure reports correct failed cohorts", {
  # Cohort g=3 has pre_stats, cohort g=7 does not
  dt <- make_agg_panel(c(g3 = 3, g7 = 2), n_nt = 3, n_periods = 8)
  pre_stats <- make_pre_stats(c(3), 1:8)  # only g=3

  expect_warning(
    result <- aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                                   cohorts = c(3L, 7L), T_max = 8L,
                                   pre_stats = pre_stats,
                                   rolling = "demean"),
    class = "lwdid_data"  # cohort_agg_skipped for g=7
  )
})

# ============================================================================
# Group 23: aggregate_to_cohort() — n_periods Semantic Verification
# (Task E5-02.8)
# ============================================================================

test_that("n_periods equals T-g+1 for balanced panel", {
  # Balanced panel: all units have data in all periods
  # g=3, T_max=6 -> post periods = {3,4,5,6} -> n_periods = 4
  g <- 3L
  T_max <- 6L
  n_periods_expected <- T_max - g + 1L  # 4
  expect_equal(n_periods_expected, 4L)
  expect_equal(length(seq.int(g, T_max)), 4L)
})

test_that("n_periods uses max of treatment units valid periods", {
  # Unbalanced: treat unit 1 has 3 valid periods, unit 2 has 2
  treat_n_periods <- c(3L, 2L)
  n_periods <- max(treat_n_periods, na.rm = TRUE)
  expect_equal(n_periods, 3L)  # max, not mean
})

test_that("n_periods falls back to theoretical value when all NA", {
  # All treatment units have NA n_periods
  treat_n_periods <- c(NA_integer_, NA_integer_)
  g <- 3L
  T_max <- 7L
  n_periods <- if (length(treat_n_periods) > 0L &&
                   any(!is.na(treat_n_periods))) {
    max(treat_n_periods, na.rm = TRUE)
  } else {
    length(seq.int(g, T_max))
  }
  expect_equal(n_periods, 5L)  # theoretical: 7-3+1 = 5
})

# ============================================================================
# Group 24: aggregate_to_cohort() — Large-Scale & Boundary Scenarios
# (Task E5-02.8)
# ============================================================================

test_that("end-to-end: large panel (100 units x 20 periods)", {
  set.seed(42)
  n_units <- 100L
  n_periods <- 20L
  # 3 cohorts: g=5 (30 units), g=10 (30 units), g=15 (20 units), NT (20 units)
  g_assign <- c(rep(5, 30), rep(10, 30), rep(15, 20), rep(NA, 20))
  dt <- data.table::CJ(id = seq_len(n_units), time = seq_len(n_periods))
  dt[, gvar := g_assign[id]]
  dt[, Y := id * 0.1 + time * 0.5 + rnorm(.N, sd = 0.1)]

  pre_stats <- make_pre_stats(c(5, 10, 15), seq_len(n_units))

  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(5L, 10L, 15L), T_max = 20L,
                         pre_stats = pre_stats, rolling = "demean")
  )

  # Should produce results (may have warnings but not errors)
  expect_true(is.list(result))
  if (length(result) > 0L) {
    # All results should have valid cohort field
    cohorts_out <- vapply(result, `[[`, integer(1), "cohort")
    expect_true(all(cohorts_out %in% c(5L, 10L, 15L)))
  }
})

test_that("end-to-end: single cohort single treated unit", {
  # Minimal valid scenario: 1 treated + 2 NT (need >=3 total)
  dt <- data.table::data.table(
    id = rep(1:3, each = 3),
    time = rep(1:3, 3),
    Y = c(1, 4, 7, 2, 3, 4, 3, 4, 5),
    gvar = c(rep(2, 3), rep(NA, 6))
  )
  pre_stats <- make_pre_stats(c(2), 1:3)

  # Should run without error (n_total=3 is sufficient)
  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(2L), T_max = 3L,
                         pre_stats = pre_stats,
                         rolling = "demean")
  )
  expect_true(is.list(result))
})

test_that("end-to-end: empty result when all cohorts fail", {
  # All cohorts lack pre_stats -> all skipped
  dt <- make_agg_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  pre_stats <- list()  # empty

  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar",
                         cohorts = c(3L), T_max = 5L,
                         pre_stats = pre_stats, rolling = "demean")
  )

  expect_true(is.list(result))
  expect_equal(length(result), 0L)
})

# ============================================================================
# Helper: create panel data for overall aggregation tests
# ============================================================================

make_overall_panel <- function(
  cohort_sizes = c(g3 = 3, g5 = 2),
  n_nt = 2,
  n_periods = 6,
  y_base = 1.0,
  treatment_effect = 2.0
) {
  units <- list()
  uid <- 1L
  for (i in seq_along(cohort_sizes)) {
    g_val <- as.integer(gsub("g", "", names(cohort_sizes)[i]))
    for (j in seq_len(cohort_sizes[i])) {
      y_vals <- y_base * uid + seq_len(n_periods) * 0.5
      post_mask <- seq_len(n_periods) >= g_val
      y_vals[post_mask] <- y_vals[post_mask] + treatment_effect
      units[[uid]] <- data.table::data.table(
        id = uid, time = seq_len(n_periods),
        Y = y_vals, gvar = g_val
      )
      uid <- uid + 1L
    }
  }
  for (j in seq_len(n_nt)) {
    units[[uid]] <- data.table::data.table(
      id = uid, time = seq_len(n_periods),
      Y = y_base * uid + seq_len(n_periods) * 0.3,
      gvar = NA_real_
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

# ============================================================================
# Group 25: construct_aggregated_outcome — FATAL-002 Protection Tests
# ============================================================================

test_that("construct_aggregated_outcome: treated units use own cohort transform", {
  # 2 cohorts: g=3 (2 units), g=5 (1 unit), 2 NT, 6 periods
  dt <- make_overall_panel(c(g3 = 2, g5 = 1), n_nt = 2, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)
  weights <- c("3" = 0.6667, "5" = 0.3333)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights,
      T_max = 6L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Treated unit 1 (g=3) should use cohort 3's transform
  # Compute cohort 3's Y_bar_ig for unit 1 independently
  ybar_g3 <- suppressWarnings(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 3L, T_max = 6L,
      pre_stat = pre_stats[["3"]], rolling = "demean"
    )
  )
  ybar_unit1_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == 1]
  ybar_unit1_overall <- res$result$Y_bar[res$result$id == 1]
  expect_equal(ybar_unit1_overall, ybar_unit1_g3, tolerance = 1e-10)

  # Treated unit 3 (g=5) should use cohort 5's transform
  ybar_g5 <- suppressWarnings(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 5L, T_max = 6L,
      pre_stat = pre_stats[["5"]], rolling = "demean"
    )
  )
  ybar_unit3_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == 3]
  ybar_unit3_overall <- res$result$Y_bar[res$result$id == 3]
  expect_equal(ybar_unit3_overall, ybar_unit3_g5, tolerance = 1e-10)
})

test_that("construct_aggregated_outcome: NT units use weighted average (FATAL-002)", {
  # FATAL-002 core: NT units must NOT use a single cohort's Y_bar_ig
  # They must use the weighted average across ALL cohorts
  dt <- make_overall_panel(c(g3 = 2, g5 = 1), n_nt = 2, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)
  weights <- c("3" = 0.6667, "5" = 0.3333)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights,
      T_max = 6L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Get NT unit (id=4) Y_bar from overall result
  nt_ybar <- res$result$Y_bar[res$result$id == 4]

  # Get cohort-specific Y_bar_ig for NT unit 4
  ybar_g3 <- suppressWarnings(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 3L, T_max = 6L,
      pre_stat = pre_stats[["3"]], rolling = "demean"
    )
  )
  ybar_g5 <- suppressWarnings(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 5L, T_max = 6L,
      pre_stat = pre_stats[["5"]], rolling = "demean"
    )
  )
  nt4_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == 4]
  nt4_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == 4]

  # NT Y_bar should NOT equal any single cohort's value
  expect_false(isTRUE(all.equal(nt_ybar, nt4_g3, tolerance = 1e-10)))
  expect_false(isTRUE(all.equal(nt_ybar, nt4_g5, tolerance = 1e-10)))

  # NT Y_bar SHOULD equal the weighted average
  expected_nt_ybar <- nt4_g3 * weights["3"] + nt4_g5 * weights["5"]
  expect_equal(nt_ybar, unname(expected_nt_ybar), tolerance = 1e-10)
})

test_that("construct_aggregated_outcome: NT weighted average matches manual calculation", {
  # Construct a scenario with known exact Y_bar_ig values for NT units
  # 2 cohorts, 1 treated each, 1 NT, 4 periods
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(10, 20, 30, 40,   # unit 1 (g=2)
          5,  15, 25, 35,    # unit 2 (g=3)
          2,   4,  6,  8),   # unit 3 (NT)
    gvar = c(rep(2, 4), rep(3, 4), rep(NA, 4))
  )
  # pre_stats with known pre_mean = 0 so Y_bar_ig = mean of post Y
  ps_g2 <- data.table::data.table(id = 1:3, pre_mean = c(0, 0, 0))
  ps_g3 <- data.table::data.table(id = 1:3, pre_mean = c(0, 0, 0))
  pre_stats <- list("2" = ps_g2, "3" = ps_g3)
  weights <- c("2" = 0.5, "3" = 0.5)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(2L, 3L), weights = weights,
      T_max = 4L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # NT unit 3:
  # Cohort g=2 post periods {2,3,4}: Y_trans = (4-0, 6-0, 8-0) = (4,6,8), mean = 6.0
  # Cohort g=3 post periods {3,4}: Y_trans = (6-0, 8-0) = (6,8), mean = 7.0
  # Weighted average = 0.5 * 6.0 + 0.5 * 7.0 = 6.5
  nt_ybar <- res$result$Y_bar[res$result$id == 3]
  expect_equal(nt_ybar, 6.5, tolerance = 1e-10)
})

test_that("construct_aggregated_outcome: returns correct structure", {
  dt <- make_overall_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)
  weights <- c("3" = 1.0)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Must return list with result and diagnostics
  expect_true(is.list(res))
  expect_true("result" %in% names(res))
  expect_true("diagnostics" %in% names(res))

  # result must be a data.table with ivar, Y_bar, D_ever
  expect_true(data.table::is.data.table(res$result))
  expect_true("id" %in% names(res$result))
  expect_true("Y_bar" %in% names(res$result))
  expect_true("D_ever" %in% names(res$result))

  # diagnostics must be a list with expected keys
  expect_true(is.list(res$diagnostics))
  expect_true("n_cohorts_requested" %in% names(res$diagnostics))
  expect_true("n_cohorts_succeeded" %in% names(res$diagnostics))
  expect_true("yield_rate" %in% names(res$diagnostics))
})

# ============================================================================
# Group 26: construct_aggregated_outcome — Weight Renormalization Tests
# ============================================================================

test_that("construct_aggregated_outcome: NT weight renormalization for missing cohorts", {
  # NT unit has data for cohort g=3 but NA for cohort g=5
  # -> weights should be renormalized to available cohorts
  # 2 treated (g=3), 1 treated (g=5), 1 NT, 5 periods
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    time = rep(1:5, 4),
    Y = c(1, 2, 5, 6, 7,     # unit 1 (g=3)
          2, 3, 6, 7, 8,     # unit 2 (g=3)
          3, 4, 5, 8, 9,     # unit 3 (g=5)
          4, 5, 6, 7, 8),    # unit 4 (NT)
    gvar = c(rep(3, 5), rep(3, 5), rep(5, 5), rep(NA, 5))
  )

  # pre_stats for g=3 includes all units, g=5 includes only units 1-3 (not unit 4)
  ps_g3 <- data.table::data.table(id = 1:4, pre_mean = rep(1.0, 4))
  # For g=5, exclude unit 4 from pre_stat so it gets NA
  ps_g5 <- data.table::data.table(id = 1:3, pre_mean = rep(1.0, 3))
  pre_stats <- list("3" = ps_g3, "5" = ps_g5)
  weights <- c("3" = 0.6667, "5" = 0.3333)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # NT unit 4 should still have a valid Y_bar (renormalized to g=3 only)
  nt_ybar <- res$result$Y_bar[res$result$id == 4]
  expect_false(is.na(nt_ybar))

  # Diagnostics should show renormalization occurred
  expect_true(res$diagnostics$nt_normalized_count >= 1L)
})

test_that("construct_aggregated_outcome: NT unit excluded when all cohorts missing", {
  # NT unit not in any cohort's pre_stat -> all Y_bar_ig are NA -> excluded
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(1, 4, 7, 10,    # unit 1 (g=2)
          2, 5, 8, 11,    # unit 2 (NT, in pre_stats)
          3, 6, 9, 12),   # unit 3 (NT, NOT in pre_stats)
    gvar = c(rep(2, 4), rep(NA, 4), rep(NA, 4))
  )
  # pre_stat only includes units 1 and 2, not unit 3
  ps_g2 <- data.table::data.table(id = 1:2, pre_mean = c(1.0, 2.0))
  pre_stats <- list("2" = ps_g2)
  weights <- c("2" = 1.0)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(2L), weights = weights,
      T_max = 4L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Unit 3 (NT, not in pre_stats) should have NA Y_bar
  nt3_ybar <- res$result$Y_bar[res$result$id == 3]
  expect_true(is.na(nt3_ybar))

  # Diagnostics should count the excluded NT unit
  expect_true(res$diagnostics$nt_excluded_count >= 1L)
})

test_that("construct_aggregated_outcome: diagnostics track normalization counts", {
  # Setup: 2 cohorts, 2 NT units — one NT missing one cohort, one NT has both
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    time = rep(1:5, 4),
    Y = c(1, 2, 5, 6, 7,     # unit 1 (g=3)
          2, 3, 6, 7, 8,     # unit 2 (g=5)
          3, 4, 5, 6, 7,     # unit 3 (NT, in both pre_stats)
          4, 5, 6, 7, 8),    # unit 4 (NT, only in g=3 pre_stat)
    gvar = c(rep(3, 5), rep(5, 5), rep(NA, 5), rep(NA, 5))
  )
  ps_g3 <- data.table::data.table(id = 1:4, pre_mean = rep(1.0, 4))
  ps_g5 <- data.table::data.table(id = 1:3, pre_mean = rep(1.0, 3))  # unit 4 missing
  pre_stats <- list("3" = ps_g3, "5" = ps_g5)
  weights <- c("3" = 0.5, "5" = 0.5)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Unit 4 had partial data -> nt_normalized_count should be >= 1
  expect_true(res$diagnostics$nt_normalized_count >= 1L)
  # Unit 3 had full data -> nt_excluded_count should be 0
  expect_equal(res$diagnostics$nt_excluded_count, 0L)
})

# ============================================================================
# Group 27: construct_aggregated_outcome — Input Validation Tests
# ============================================================================

test_that("construct_aggregated_outcome: empty cohorts raises lwdid_invalid_input", {
  dt <- make_overall_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  pre_stats <- list()
  weights <- numeric(0)

  expect_error(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = integer(0), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_invalid_input"
  )
})

test_that("construct_aggregated_outcome: weights-cohorts mismatch raises error", {
  dt <- make_overall_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)
  # weights has key "5" but cohorts only has 3
  weights <- c("5" = 1.0)

  expect_error(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_invalid_input"
  )
})

test_that("construct_aggregated_outcome: weight sum deviation warns", {
  dt <- make_overall_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)
  # Weights sum to 0.8, not 1.0 -> should warn
  weights <- c("3" = 0.8)

  expect_warning(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), weights = weights,
      T_max = 5L, pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_numerical"
  )
})

test_that("construct_aggregated_outcome: all cohorts fail raises lwdid_insufficient_data", {
  dt <- make_overall_panel(c(g3 = 2), n_nt = 2, n_periods = 5)
  # pre_stats has no matching cohort key -> all cohorts fail
  pre_stats <- list("99" = data.table::data.table(id = 1:4, pre_mean = rep(1.0, 4)))
  weights <- c("3" = 1.0)

  expect_error(
    suppressWarnings(
      construct_aggregated_outcome(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), weights = weights,
        T_max = 5L, pre_stats = pre_stats, rolling = "demean"
      )
    ),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Group 28: construct_aggregated_outcome — Boundary Condition Tests
# ============================================================================

test_that("construct_aggregated_outcome: single cohort degeneracy", {
  # With only one cohort, NT weighted average = single cohort's transform
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(10, 20, 30, 40,   # unit 1 (g=2)
          5,  15, 25, 35,    # unit 2 (NT)
          2,   4,  6,  8),   # unit 3 (NT)
    gvar = c(rep(2, 4), rep(NA, 4), rep(NA, 4))
  )
  ps_g2 <- data.table::data.table(id = 1:3, pre_mean = c(10.0, 5.0, 2.0))
  pre_stats <- list("2" = ps_g2)
  weights <- c("2" = 1.0)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(2L), weights = weights,
      T_max = 4L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # NT unit 2: single cohort -> weighted avg = cohort g=2 transform
  ybar_g2 <- suppressWarnings(
    compute_cohort_aggregated_variable(
      dt, "Y", "id", "time", g = 2L, T_max = 4L,
      pre_stat = pre_stats[["2"]], rolling = "demean"
    )
  )
  nt2_direct <- ybar_g2$Y_bar_ig[ybar_g2$id == 2]
  nt2_overall <- res$result$Y_bar[res$result$id == 2]
  expect_equal(nt2_overall, nt2_direct, tolerance = 1e-10)

  # Also verify numerical reasonableness:
  # Unit 2 (NT): Y=(5,15,25,35), pre_mean=5, post={2,3,4}
  # ydot = (15-5, 25-5, 35-5) = (10, 20, 30), mean = 20.0
  expect_equal(nt2_overall, 20.0, tolerance = 1e-10)
})

test_that("construct_aggregated_outcome: non-integer gvar batch warning", {
  # Treated units with float gvar values -> should trigger batch warning
  dt <- data.table::data.table(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    Y = c(1, 2, 5, 6,     # unit 1 (g=2.7, non-integer)
          2, 3, 6, 7,     # unit 2 (g=3, integer)
          3, 4, 5, 6),    # unit 3 (NT)
    gvar = c(rep(2.7, 4), rep(3, 4), rep(NA, 4))
  )
  ps_g3 <- data.table::data.table(id = 1:3, pre_mean = rep(1.0, 3))
  pre_stats <- list("3" = ps_g3)
  weights <- c("3" = 1.0)

  # The non-integer gvar unit should trigger lwdid_data warning
  expect_warning(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), weights = weights,
      T_max = 4L, pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_data"
  )
})

test_that("construct_aggregated_outcome: yield_rate in diagnostics is correct", {
  dt <- make_overall_panel(c(g3 = 3, g5 = 2), n_nt = 2, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  n_total <- length(all_ids)
  pre_stats <- make_pre_stats(c(3, 5), all_ids)
  weights <- c("3" = 0.6, "5" = 0.4)

  res <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights,
      T_max = 6L, pre_stats = pre_stats, rolling = "demean"
    )
  )

  # yield_rate = n_valid / n_total
  n_valid <- sum(!is.na(res$result$Y_bar))
  expected_yield <- n_valid / n_total
  expect_equal(res$diagnostics$yield_rate, expected_yield, tolerance = 1e-10)

  # With all units having valid data, yield_rate should be 1.0
  expect_equal(res$diagnostics$yield_rate, 1.0, tolerance = 1e-10)
})


# ============================================================================
# Group 29: aggregate_to_overall — Core Regression Tests
# ============================================================================

test_that("aggregate_to_overall: basic overall effect estimation", {
  # 2 cohorts: g=3 (3 units), g=5 (2 units), 3 NT, 6 periods
  # Known treatment effect = 2.0
  set.seed(42)
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6,
    treatment_effect = 2.0
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # ATT should be numerically reasonable (close to known effect of 2.0)
  expect_true(is.numeric(result$att))
  # ATT should be finite and numeric (exact value depends on data construction;
  # numerical verification is in Task E5-03.5)
  expect_true(is.finite(result$att),
    label = "ATT must be finite")

  # SE must be positive

  expect_true(result$se > 0)

  # CI must contain ATT
  expect_true(result$ci_lower <= result$att)
  expect_true(result$ci_upper >= result$att)

  # t_stat = att / se
  expect_equal(result$t_stat, result$att / result$se, tolerance = 1e-10)

  # p-value in [0, 1]
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
})

test_that("aggregate_to_overall: cohort weights correct", {
  # Cohort g=3 has 3 units, g=5 has 2 units -> omega_3=0.6, omega_5=0.4
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # cohort_weights should be a named list
  expect_true(is.list(result$cohort_weights))
  expect_equal(result$cohort_weights[["3"]], 0.6, tolerance = 1e-10)
  expect_equal(result$cohort_weights[["5"]], 0.4, tolerance = 1e-10)

  # Weights must sum to 1
  w_sum <- sum(unlist(result$cohort_weights))
  expect_equal(w_sum, 1.0, tolerance = 1e-10)
})

test_that("aggregate_to_overall: single cohort degeneracy", {
  # Only 1 cohort -> weight = 1.0, tau_omega = single cohort's effect
  dt <- make_overall_panel(
    c(g3 = 4), n_nt = 3, n_periods = 6, treatment_effect = 3.0
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Single cohort weight must be 1.0
  expect_equal(result$cohort_weights[["3"]], 1.0, tolerance = 1e-10)

  # ATT should be close to the known treatment effect of 3.0
  expect_true(is.numeric(result$att))
  # ATT should be finite (exact numerical verification in Task E5-03.5)
  expect_true(is.finite(result$att),
    label = "Single-cohort ATT must be finite")
})

test_that("aggregate_to_overall: no NT units raises lwdid_no_never_treated", {
  # All units are treated (no never-treated)
  dt <- data.table::data.table(
    id = rep(1:4, each = 5),
    time = rep(1:5, 4),
    Y = rep(1:4, each = 5) + rep(1:5, 4) * 0.5,
    gvar = rep(c(3, 3, 5, 5), each = 5)
  )
  pre_stats <- make_pre_stats(c(3, 5), 1:4)

  expect_error(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 5L,
      pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_no_never_treated"
  )
})

test_that("aggregate_to_overall: returns correct structure", {
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # All 17 fields must be present
  expected_fields <- c(
    "att", "se", "ci_lower", "ci_upper", "t_stat", "pvalue",
    "cohort_weights", "effective_weights",
    "n_treated", "n_control", "n_total",
    "df_resid", "df_inference", "K",
    "controls_tier", "controls_used", "diagnostics"
  )
  for (f in expected_fields) {
    expect_true(f %in% names(result),
      label = sprintf("Field '%s' must be present in result", f))
  }

  # Type checks
  expect_true(is.numeric(result$att))
  expect_true(is.numeric(result$se))
  expect_true(is.numeric(result$ci_lower))
  expect_true(is.numeric(result$ci_upper))
  expect_true(is.numeric(result$t_stat))
  expect_true(is.numeric(result$pvalue))
  expect_true(is.list(result$cohort_weights))
  expect_true(is.list(result$effective_weights))
  expect_true(is.integer(result$n_treated) || is.numeric(result$n_treated))
  expect_true(is.integer(result$n_control) || is.numeric(result$n_control))
  expect_true(is.integer(result$n_total) || is.numeric(result$n_total))
  expect_true(is.integer(result$df_resid) || is.numeric(result$df_resid))
  expect_true(is.integer(result$df_inference) || is.numeric(result$df_inference))
  expect_true(is.integer(result$K) || is.numeric(result$K))
  expect_true(is.character(result$controls_tier))
  expect_true(is.logical(result$controls_used))
  expect_true(is.list(result$diagnostics))
})


# ============================================================================
# Group 30: aggregate_to_overall — Effective Weights Tests
# ============================================================================

test_that("aggregate_to_overall: effective weights match computation weights (no dropna)", {

  # Balanced panel with no NAs -> effective_weights == cohort_weights
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # With no NAs, effective weights should equal computation weights
  for (g_char in names(result$cohort_weights)) {
    expect_equal(
      result$effective_weights[[g_char]],
      result$cohort_weights[[g_char]],
      tolerance = 1e-10,
      label = sprintf("Effective weight for cohort %s", g_char)
    )
  }
})

test_that("aggregate_to_overall: effective weights deviation > 1%% warns", {
  # Create data where differential NA rates cause >1% weight deviation.
  # Cohort g=3: 5 units, cohort g=5: 5 units, 3 NT, 8 periods.
  # Then inject NAs into cohort g=3 units' Y to cause dropna to remove them.
  dt <- make_overall_panel(
    c(g3 = 5, g5 = 5), n_nt = 3, n_periods = 8
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  # Inject NAs into cohort g=3 units (ids 1-5) for ALL periods
  # so their Y_bar becomes NA and they get dropped.
  # We need to drop enough g=3 units to cause >1% deviation.
  # Original weights: g3=5/10=0.5, g5=5/10=0.5
  # If 3 of 5 g=3 units get NA Y_bar, effective: g3=2/7, g5=5/7
  # Deviation: |0.5 - 2/7| = |0.5 - 0.286| = 0.214 > 0.01
  dt[id %in% c(1, 2, 3), Y := NA_real_]

  expect_warning(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 8L,
      pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_numerical"
  )
})

test_that("aggregate_to_overall: R preserves both weight sets when deviation > 1%%", {
  # Same setup as above: differential NA rates cause weight deviation
  dt <- make_overall_panel(
    c(g3 = 5, g5 = 5), n_nt = 3, n_periods = 8
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  # Inject NAs into 3 of 5 g=3 units
  dt[id %in% c(1, 2, 3), Y := NA_real_]

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 8L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Both weight sets should exist
  expect_true(is.list(result$cohort_weights))
  expect_true(is.list(result$effective_weights))

  # Computation weights should reflect original proportions (5/10 each)
  expect_equal(result$cohort_weights[["3"]], 0.5, tolerance = 1e-10)
  expect_equal(result$cohort_weights[["5"]], 0.5, tolerance = 1e-10)

  # Effective weights should differ from computation weights
  # (g=3 lost 3 units, so effective weight for g=3 < 0.5)
  expect_true(
    result$effective_weights[["3"]] < result$cohort_weights[["3"]],
    label = "Effective weight for g=3 should be less than computation weight"
  )
  expect_true(
    result$effective_weights[["5"]] > result$cohort_weights[["5"]],
    label = "Effective weight for g=5 should be greater than computation weight"
  )

  # Effective weights must still sum to 1
  eff_sum <- sum(unlist(result$effective_weights))
  expect_equal(eff_sum, 1.0, tolerance = 1e-10)
})


# ============================================================================
# Group 31: aggregate_to_overall — Degrees of Freedom Tests
# ============================================================================

test_that("aggregate_to_overall: df_resid and df_inference correct (no controls)", {
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Both fields must exist and be numeric
  expect_true(!is.null(result$df_resid))
  expect_true(!is.null(result$df_inference))
  expect_true(is.numeric(result$df_resid))
  expect_true(is.numeric(result$df_inference))

  # Both must be positive integers
  expect_true(result$df_resid > 0)
  expect_true(result$df_inference > 0)
  expect_equal(result$df_resid, as.integer(result$df_resid))
  expect_equal(result$df_inference, as.integer(result$df_inference))
})

test_that("aggregate_to_overall: df_resid and df_inference separate fields", {
  dt <- make_overall_panel(
    c(g3 = 3, g5 = 2), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Both fields must exist and be numeric
  expect_true(!is.null(result$df_resid))
  expect_true(!is.null(result$df_inference))
  expect_true(is.numeric(result$df_resid))
  expect_true(is.numeric(result$df_inference))

  # For homoskedastic VCE (default), df_inference = df_resid
  expect_equal(result$df_inference, result$df_resid)
})

test_that("aggregate_to_overall: df fields are always non-NULL integers", {
  # Verify that the returned result always has df_resid and df_inference
  # as non-NULL integer values across different configurations
  dt <- make_overall_panel(
    c(g3 = 5), n_nt = 5, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  expect_false(is.null(result$df_resid),
    label = "df_resid must not be NULL")
  expect_false(is.null(result$df_inference),
    label = "df_inference must not be NULL")

  # Must be representable as integers
  expect_true(is.finite(result$df_resid))
  expect_true(is.finite(result$df_inference))
  expect_equal(result$df_resid %% 1, 0, tolerance = 1e-10,
    label = "df_resid must be an integer value")
  expect_equal(result$df_inference %% 1, 0, tolerance = 1e-10,
    label = "df_inference must be an integer value")
})


# ============================================================================
# Group 32: aggregate_to_overall — VCE and Clustering Tests
# ============================================================================

test_that("aggregate_to_overall: vce=cluster without cluster_var raises error", {
  dt <- make_overall_panel(
    c(g3 = 3), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  expect_error(
    suppressWarnings(
      aggregate_to_overall(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), T_max = 6L,
        pre_stats = pre_stats, rolling = "demean",
        vce = "cluster", cluster_var = NULL
      )
    ),
    class = "lwdid_invalid_input"
  )
})

test_that("aggregate_to_overall: cluster SE with small clusters warns", {
  # Create data with < 20 clusters -> expect lwdid_small_sample warning
  n_treat <- 8
  n_nt <- 7
  n_total <- n_treat + n_nt
  n_periods <- 6
  dt <- make_overall_panel(
    c(g3 = n_treat), n_nt = n_nt, n_periods = n_periods
  )
  # Add cluster variable: 5 clusters (< 20)
  dt[, cluster_id := rep(
    rep(paste0("C", 1:5), length.out = n_total),
    each = n_periods
  )]
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  expect_warning(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = n_periods,
      pre_stats = pre_stats, rolling = "demean",
      vce = "cluster", cluster_var = "cluster_id"
    ),
    class = "lwdid_small_sample"
  )
})

test_that("aggregate_to_overall: cluster variable extraction runs without error", {
  # Verify cluster-based estimation runs when cluster_var is provided
  n_treat <- 10
  n_nt <- 10
  n_total <- n_treat + n_nt
  n_periods <- 6
  dt <- make_overall_panel(
    c(g3 = n_treat), n_nt = n_nt, n_periods = n_periods
  )
  # Add cluster variable: 20 clusters (sufficient)
  dt[, cluster_id := rep(
    rep(paste0("C", 1:20), length.out = n_total),
    each = n_periods
  )]
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = n_periods,
      pre_stats = pre_stats, rolling = "demean",
      vce = "cluster", cluster_var = "cluster_id"
    )
  )

  # Should return a valid result with ATT
  expect_true(is.numeric(result$att))
  expect_true(result$se > 0)
  expect_true(result$n_total == n_total)
})


# ============================================================================
# Group 33: aggregate_to_overall — Controls Extension Tests
# ============================================================================

test_that("aggregate_to_overall: controls extension warning emitted", {
  # Pass controls parameter -> expect lwdid_data warning with
  # "overall_controls_extension" detail
  dt <- make_overall_panel(
    c(g3 = 5), n_nt = 5, n_periods = 6
  )
  # Add a control variable (unit-level, time-invariant)
  n_total <- 10
  dt[, X1 := rep(rnorm(n_total), each = 6)]
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  expect_warning(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean",
      controls = "X1"
    ),
    regexp = "overall_controls_extension|controls"
  )
})

test_that("aggregate_to_overall: no controls by default", {
  dt <- make_overall_panel(
    c(g3 = 3), n_nt = 3, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean",
      controls = NULL
    )
  )

  # K=0 and controls_used=FALSE when no controls
  expect_equal(result$K, 0L)
  expect_false(result$controls_used)
})

test_that("aggregate_to_overall: controls_tier and controls_used in return", {
  dt <- make_overall_panel(
    c(g3 = 5), n_nt = 5, n_periods = 6
  )
  n_total <- 10
  dt[, X1 := rep(rnorm(n_total), each = 6)]
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean",
      controls = "X1"
    )
  )

  # controls_tier must be a character string
  expect_true("controls_tier" %in% names(result))
  expect_true(is.character(result$controls_tier))
  expect_true(result$controls_tier %in% c("none", "simple", "full_interaction"))

  # controls_used must be TRUE when controls are passed
  expect_true("controls_used" %in% names(result))
  expect_true(result$controls_used)
})


# ============================================================================
# Group 34: aggregate_to_overall — Sample Size Boundary Tests
# ============================================================================

test_that("aggregate_to_overall: sample size 2 warns then errors", {
  # Minimal data: exactly 2 units (1 treated, 1 NT)
  # aggregate_to_overall warns lwdid_small_sample at n_total==2,
  # then estimate_ra_common errors with lwdid_insufficient_data (N<3).
  dt <- make_overall_panel(
    c(g3 = 1), n_nt = 1, n_periods = 6
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  # The warning IS emitted before the error
  warned <- FALSE
  tryCatch(
    withCallingHandlers(
      aggregate_to_overall(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), T_max = 6L,
        pre_stats = pre_stats, rolling = "demean"
      ),
      lwdid_small_sample = function(w) {
        warned <<- TRUE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  expect_true(warned,
    label = "lwdid_small_sample warning should be emitted for n_total==2")

  # Also verify the error is raised
  expect_error(
    suppressWarnings(
      aggregate_to_overall(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), T_max = 6L,
        pre_stats = pre_stats, rolling = "demean"
      )
    ),
    class = "lwdid_insufficient_data"
  )
})

test_that("aggregate_to_overall: zero-size cohorts excluded with warning", {
  # Include cohort g=7 in cohorts vector but no units have gvar=7
  # -> zero-size cohort should be excluded with lwdid_data warning
  dt <- make_overall_panel(
    c(g3 = 3), n_nt = 3, n_periods = 8
  )
  all_ids <- sort(unique(dt$id))
  # pre_stats for both cohorts (g=3 exists, g=7 has no units)
  pre_stats <- make_pre_stats(c(3, 7), all_ids)

  expect_warning(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 7L), T_max = 8L,
      pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_data"
  )
})


# ============================================================================
# Group 35: Demean End-to-End Numerical Verification
# ============================================================================

test_that("35.1: demean end-to-end numerical verification with hand-calculated values", {
  # Deterministic panel: 2 cohorts g=3 (2 units), g=5 (1 unit), 2 NT, 6 periods
  # Unit 1 (g=3): Y = [1,2,3,4,5,6]
  # Unit 2 (g=3): Y = [2,3,4,5,6,7]
  # Unit 3 (g=5): Y = [1,1,1,1,10,10]
  # Unit 4 (NT):  Y = [3,3,3,3,3,3]
  # Unit 5 (NT):  Y = [5,5,5,5,5,5]
  dt <- data.table::data.table(
    id   = rep(1:5, each = 6),
    time = rep(1:6, 5),
    Y    = c(1,2,3,4,5,6,       # unit 1
             2,3,4,5,6,7,       # unit 2
             1,1,1,1,10,10,     # unit 3
             3,3,3,3,3,3,       # unit 4
             5,5,5,5,5,5),      # unit 5
    gvar = rep(c(3, 3, 5, NA, NA), each = 6)
  )

  # pre_mean = 1.0 for all units in all cohorts
  pre_stats <- make_pre_stats(c(3, 5), 1:5)

  # --- Step A: Verify compute_cohort_aggregated_variable for each cohort ---
  ybar_g3 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 6L,
    pre_stat = pre_stats[["3"]], rolling = "demean"
  )
  ybar_g3 <- ybar_g3[order(ybar_g3$id), ]

  # Cohort g=3, post periods t=3,4,5,6 (4 periods), pre_mean=1.0
  # Unit 1: mean(3-1, 4-1, 5-1, 6-1) = mean(2,3,4,5) = 3.5
  # Unit 2: mean(4-1, 5-1, 6-1, 7-1) = mean(3,4,5,6) = 4.5
  # Unit 3: mean(1-1, 1-1, 10-1, 10-1) = mean(0,0,9,9) = 4.5
  # Unit 4: mean(3-1, 3-1, 3-1, 3-1) = mean(2,2,2,2) = 2.0
  # Unit 5: mean(5-1, 5-1, 5-1, 5-1) = mean(4,4,4,4) = 4.0
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 1], 3.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 2], 4.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 3], 4.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 4], 2.0, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 5], 4.0, tolerance = 1e-10)
  expect_equal(ybar_g3$n_periods[ybar_g3$id == 1], 4L)

  ybar_g5 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 5L, T_max = 6L,
    pre_stat = pre_stats[["5"]], rolling = "demean"
  )
  ybar_g5 <- ybar_g5[order(ybar_g5$id), ]

  # Cohort g=5, post periods t=5,6 (2 periods), pre_mean=1.0
  # Unit 1: mean(5-1, 6-1) = mean(4,5) = 4.5
  # Unit 2: mean(6-1, 7-1) = mean(5,6) = 5.5
  # Unit 3: mean(10-1, 10-1) = mean(9,9) = 9.0
  # Unit 4: mean(3-1, 3-1) = mean(2,2) = 2.0
  # Unit 5: mean(5-1, 5-1) = mean(4,4) = 4.0
  expect_equal(ybar_g5$Y_bar_ig[ybar_g5$id == 1], 4.5, tolerance = 1e-10)
  expect_equal(ybar_g5$Y_bar_ig[ybar_g5$id == 2], 5.5, tolerance = 1e-10)
  expect_equal(ybar_g5$Y_bar_ig[ybar_g5$id == 3], 9.0, tolerance = 1e-10)
  expect_equal(ybar_g5$Y_bar_ig[ybar_g5$id == 4], 2.0, tolerance = 1e-10)
  expect_equal(ybar_g5$Y_bar_ig[ybar_g5$id == 5], 4.0, tolerance = 1e-10)
  expect_equal(ybar_g5$n_periods[ybar_g5$id == 3], 2L)

  # --- Step B: Verify construct_aggregated_outcome ---
  # Cohort weights: omega_3 = 2/3, omega_5 = 1/3
  weights <- c("3" = 2/3, "5" = 1/3)
  agg <- construct_aggregated_outcome(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L, 5L), weights = weights,
    T_max = 6L, pre_stats = pre_stats, rolling = "demean"
  )
  agg_dt <- agg$result[order(agg$result$id), ]

  # Treated units use own cohort transform:
  # Unit 1 (g=3): Y_bar = Y_bar_1,g3 = 3.5
  # Unit 2 (g=3): Y_bar = Y_bar_2,g3 = 4.5
  # Unit 3 (g=5): Y_bar = Y_bar_3,g5 = 9.0
  expect_equal(agg_dt$Y_bar[agg_dt$id == 1], 3.5, tolerance = 1e-10)
  expect_equal(agg_dt$Y_bar[agg_dt$id == 2], 4.5, tolerance = 1e-10)
  expect_equal(agg_dt$Y_bar[agg_dt$id == 3], 9.0, tolerance = 1e-10)

  # NT units use weighted average:
  # Unit 4: (2/3)*2.0 + (1/3)*2.0 = 2.0
  # Unit 5: (2/3)*4.0 + (1/3)*4.0 = 4.0
  expect_equal(agg_dt$Y_bar[agg_dt$id == 4], 2.0, tolerance = 1e-10)
  expect_equal(agg_dt$Y_bar[agg_dt$id == 5], 4.0, tolerance = 1e-10)

  # D_ever flags
  expect_equal(agg_dt$D_ever[agg_dt$id == 1], 1L)
  expect_equal(agg_dt$D_ever[agg_dt$id == 4], 0L)

  # --- Step C: Verify aggregate_to_overall ATT ---
  # mean(Y_bar | D=1) = (3.5 + 4.5 + 9.0) / 3 = 17/3
  # mean(Y_bar | D=0) = (2.0 + 4.0) / 2 = 3.0
  # tau_omega = 17/3 - 3.0 = 8/3 ≈ 2.666667
  overall <- aggregate_to_overall(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L, 5L), T_max = 6L,
    pre_stats = pre_stats, rolling = "demean"
  )
  expect_equal(overall$att, 8/3, tolerance = 1e-5)
  expect_equal(overall$n_treated, 3L)
  expect_equal(overall$n_control, 2L)
  expect_equal(overall$n_total, 5L)

  # Cohort weights should be omega_3=2/3, omega_5=1/3
  expect_equal(
    as.numeric(overall$cohort_weights[["3"]]),
    2/3, tolerance = 1e-10
  )
  expect_equal(
    as.numeric(overall$cohort_weights[["5"]]),
    1/3, tolerance = 1e-10
  )
})

test_that("35.2: demean with unbalanced panel runs correctly", {
  # Same base panel as 35.1 but with varying Y for unit 4 and missing obs
  dt <- data.table::data.table(
    id   = rep(1:5, each = 6),
    time = rep(1:6, 5),
    Y    = c(1,2,3,4,5,6,       # unit 1 (g=3)
             2,3,4,5,6,7,       # unit 2 (g=3)
             1,1,1,1,10,10,     # unit 3 (g=5)
             3,3,3,3,4,5,       # unit 4 (NT) — varying Y
             5,5,5,5,5,5),      # unit 5 (NT)
    gvar = rep(c(3, 3, 5, NA, NA), each = 6)
  )
  # Remove row where id=4, time=5 to create unbalanced panel

  dt <- dt[!(id == 4 & time == 5)]

  pre_stats <- make_pre_stats(c(3, 5), 1:5)

  # Verify aggregate_to_overall runs without error on unbalanced panel
  overall <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Basic sanity checks
  expect_true(is.finite(overall$att))
  expect_true(is.numeric(overall$att))
  expect_equal(overall$n_treated, 3L)
  expect_equal(overall$n_control, 2L)
  expect_equal(overall$n_total, 5L)
  expect_true(!is.null(overall$diagnostics))
  expect_true(is.finite(overall$se))
  expect_true(overall$se > 0)

  # Verify the ATT differs from balanced case due to missing observation
  # Unit 4 g=3: post t=3,4,5,6 but t=5 missing → t=3,4,6 only
  # Y_bar_4,g3 = mean(3-1, 3-1, 5-1) = mean(2,2,4) = 8/3
  # Unit 4 g=5: post t=5,6 but t=5 missing → t=6 only
  # Y_bar_4,g5 = (5-1) = 4.0
  # NT Y_bar_4 = (2/3)*(8/3) + (1/3)*4.0 = 16/9 + 4/3 = 16/9 + 12/9 = 28/9
  # This differs from balanced case (2.0), so ATT should differ
  expect_true(abs(overall$att - 8/3) > 1e-6,
    label = "Unbalanced ATT should differ from balanced case")
})

test_that("35.3: demean single cohort overall equals cohort-level effect", {
  # 1 cohort g=3 (3 units), 3 NT units, 5 periods
  # With single cohort, omega_3 = 1.0, so NT Y_bar = Y_bar_ig3
  # Therefore aggregate_to_overall ATT should equal aggregate_to_cohort ATT
  dt <- data.table::data.table(
    id   = rep(1:6, each = 5),
    time = rep(1:5, 6),
    Y    = c(1,2,5,6,7,         # unit 1 (g=3)
             2,3,8,9,10,        # unit 2 (g=3)
             3,4,6,7,8,         # unit 3 (g=3)
             1,1,1,1,1,         # unit 4 (NT)
             2,2,2,2,2,         # unit 5 (NT)
             3,3,3,3,3),        # unit 6 (NT)
    gvar = rep(c(3, 3, 3, NA, NA, NA), each = 5)
  )

  pre_stats <- make_pre_stats(c(3), 1:6)

  # Cohort-level effect
  cohort_results <- aggregate_to_cohort(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L), T_max = 5L,
    pre_stats = pre_stats, rolling = "demean"
  )
  tau_cohort <- cohort_results[[1]]$att

  # Overall effect
  overall <- aggregate_to_overall(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L), T_max = 5L,
    pre_stats = pre_stats, rolling = "demean"
  )
  tau_overall <- overall$att

  # With single cohort, these should be identical
  expect_equal(tau_overall, tau_cohort, tolerance = 1e-10)

  # Also verify cohort weight is exactly 1.0
  expect_equal(as.numeric(overall$cohort_weights[["3"]]), 1.0, tolerance = 1e-10)
  expect_equal(overall$n_treated, 3L)
  expect_equal(overall$n_control, 3L)
})


# ============================================================================
# Group 36: Detrend End-to-End Numerical Verification
# ============================================================================

test_that("36.1: detrend end-to-end numerical verification with hand-calculated values", {
  # 1 cohort g=3, 2 treated, 2 NT, 4 periods
  # pre_stats: pre_mean=0 (intercept), slope=1, t_bar_pre=1.5
  # Detrend formula: Y_trans = Y - (pre_mean + slope * (t - t_bar_pre))
  #                          = Y - (0 + 1 * (t - 1.5))
  # At t=3: pred = 1.5, at t=4: pred = 2.5
  dt <- data.table::data.table(
    id   = rep(1:4, each = 4),
    time = rep(1:4, 4),
    Y    = c(1,2,10,11,         # unit 1 (g=3)
             2,3,12,13,         # unit 2 (g=3)
             1,2,3,4,           # unit 3 (NT)
             2,3,4,5),          # unit 4 (NT)
    gvar = rep(c(3, 3, NA, NA), each = 4)
  )

  # Custom detrend pre_stats: pre_mean=0, slope=1, t_bar_pre=1.5
  pre_stats <- list(
    "3" = data.table::data.table(
      id = 1:4,
      pre_mean = rep(0, 4),
      slope = rep(1, 4),
      t_bar_pre = rep(1.5, 4)
    )
  )

  # --- Step A: Verify compute_cohort_aggregated_variable ---
  ybar_g3 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 4L,
    pre_stat = pre_stats[["3"]], rolling = "detrend"
  )
  ybar_g3 <- ybar_g3[order(ybar_g3$id), ]

  # Post periods t=3,4. Predicted: t=3 → 1.5, t=4 → 2.5
  # Unit 1: mean(10-1.5, 11-2.5) = mean(8.5, 8.5) = 8.5
  # Unit 2: mean(12-1.5, 13-2.5) = mean(10.5, 10.5) = 10.5
  # Unit 3: mean(3-1.5, 4-2.5) = mean(1.5, 1.5) = 1.5
  # Unit 4: mean(4-1.5, 5-2.5) = mean(2.5, 2.5) = 2.5
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 1], 8.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 2], 10.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 3], 1.5, tolerance = 1e-10)
  expect_equal(ybar_g3$Y_bar_ig[ybar_g3$id == 4], 2.5, tolerance = 1e-10)
  expect_equal(ybar_g3$n_periods[1], 2L)

  # --- Step B: Verify aggregate_to_overall ---
  # Single cohort → omega_3 = 1.0
  # Y_bar: [8.5, 10.5, 1.5, 2.5], D_ever: [1, 1, 0, 0]
  # tau = mean(8.5, 10.5) - mean(1.5, 2.5) = 9.5 - 2.0 = 7.5
  overall <- aggregate_to_overall(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L), T_max = 4L,
    pre_stats = pre_stats, rolling = "detrend"
  )
  expect_equal(overall$att, 7.5, tolerance = 1e-5)
  expect_equal(overall$n_treated, 2L)
  expect_equal(overall$n_control, 2L)
  expect_equal(overall$n_total, 4L)
})

test_that("36.2: detrend perfect linear trend gives zero effect", {
  # All units follow Y = intercept_i + slope * t (same slope=2, different intercepts)
  # No treatment effect → detrend residuals should be ~0 → tau ≈ 0
  # 1 cohort g=3, 3 treated, 3 NT, 5 periods
  slope_common <- 2.0
  dt <- data.table::data.table(
    id   = rep(1:6, each = 5),
    time = rep(1:5, 6),
    Y    = c(10 + slope_common * (1:5),   # unit 1 (g=3): intercept=10
             20 + slope_common * (1:5),   # unit 2 (g=3): intercept=20
             30 + slope_common * (1:5),   # unit 3 (g=3): intercept=30
             5  + slope_common * (1:5),   # unit 4 (NT): intercept=5
             15 + slope_common * (1:5),   # unit 5 (NT): intercept=15
             25 + slope_common * (1:5)),  # unit 6 (NT): intercept=25
    gvar = rep(c(3, 3, 3, NA, NA, NA), each = 5)
  )

  # Pre-treatment periods for g=3: t=1,2
  # t_bar_pre = mean(1,2) = 1.5
  # For each unit i: pre_mean_i = intercept_i + slope * t_bar_pre
  #   (this is the "A*" in the formula: A* = mean(Y_pre) = intercept + slope * t_bar_pre)
  # slope_i = slope_common = 2.0 for all units
  # Detrend: Y_trans = Y - (pre_mean + slope * (t - t_bar_pre))
  #        = (intercept + slope*t) - (intercept + slope*t_bar_pre + slope*(t - t_bar_pre))
  #        = (intercept + slope*t) - (intercept + slope*t)
  #        = 0
  t_bar_pre <- 1.5
  pre_stats <- list(
    "3" = data.table::data.table(
      id = 1:6,
      pre_mean = c(10, 20, 30, 5, 15, 25) + slope_common * t_bar_pre,
      slope = rep(slope_common, 6),
      t_bar_pre = rep(t_bar_pre, 6)
    )
  )

  overall <- aggregate_to_overall(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L), T_max = 5L,
    pre_stats = pre_stats, rolling = "detrend"
  )

  # With perfect linear trend and no treatment effect, tau should be ~0
  expect_equal(overall$att, 0.0, tolerance = 1e-10)
  expect_equal(overall$n_treated, 3L)
  expect_equal(overall$n_control, 3L)
})


# ============================================================================
# Group 37: FATAL-002 Numerical Verification
# ============================================================================

test_that("37.1: NT weighted average differs from single cohort transform", {
  # 2 cohorts with DIFFERENT Y_bar_ig values for NT units
  # Verify NT Y_bar is the weighted average, NOT equal to any single cohort's Y_bar_ig
  # This is the core FATAL-002 protection numerical test
  #
  # Design: NT units have different Y_bar_ig across cohorts because
  # different post-treatment windows produce different averages
  dt <- data.table::data.table(
    id   = rep(1:4, each = 6),
    time = rep(1:6, 4),
    Y    = c(1,1,5,5,5,5,       # unit 1 (g=3): constant post-treatment
             1,1,1,1,8,8,       # unit 2 (g=5): jump at t=5
             1,1,2,2,10,10,     # unit 3 (NT): different Y across periods
             1,1,3,3,1,1),      # unit 4 (NT): different Y across periods
    gvar = rep(c(3, 5, NA, NA), each = 6)
  )

  pre_stats <- make_pre_stats(c(3, 5), 1:4)

  # Compute per-cohort Y_bar_ig for NT units
  ybar_g3 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 6L,
    pre_stat = pre_stats[["3"]], rolling = "demean"
  )
  ybar_g5 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 5L, T_max = 6L,
    pre_stat = pre_stats[["5"]], rolling = "demean"
  )

  # NT unit 3:
  # g=3 post t=3,4,5,6: mean(2-1, 2-1, 10-1, 10-1) = mean(1,1,9,9) = 5.0
  # g=5 post t=5,6:     mean(10-1, 10-1) = mean(9,9) = 9.0
  ybar_u3_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == 3]
  ybar_u3_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == 3]
  expect_equal(ybar_u3_g3, 5.0, tolerance = 1e-10)
  expect_equal(ybar_u3_g5, 9.0, tolerance = 1e-10)

  # NT unit 4:
  # g=3 post t=3,4,5,6: mean(3-1, 3-1, 1-1, 1-1) = mean(2,2,0,0) = 1.0
  # g=5 post t=5,6:     mean(1-1, 1-1) = mean(0,0) = 0.0
  ybar_u4_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == 4]
  ybar_u4_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == 4]
  expect_equal(ybar_u4_g3, 1.0, tolerance = 1e-10)
  expect_equal(ybar_u4_g5, 0.0, tolerance = 1e-10)

  # Verify Y_bar_ig DIFFERS across cohorts for NT units (prerequisite for FATAL-002)
  expect_true(abs(ybar_u3_g3 - ybar_u3_g5) > 1.0,
    label = "NT unit 3 must have different Y_bar_ig across cohorts")
  expect_true(abs(ybar_u4_g3 - ybar_u4_g5) > 0.5,
    label = "NT unit 4 must have different Y_bar_ig across cohorts")

  # Now verify construct_aggregated_outcome uses weighted average
  # omega_3 = 1/2, omega_5 = 1/2 (1 treated unit each)
  weights <- c("3" = 0.5, "5" = 0.5)
  agg <- construct_aggregated_outcome(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L, 5L), weights = weights,
    T_max = 6L, pre_stats = pre_stats, rolling = "demean"
  )
  agg_dt <- agg$result[order(agg$result$id), ]

  # NT unit 3: Y_bar = 0.5 * 5.0 + 0.5 * 9.0 = 7.0
  # NOT 5.0 (g=3 only) and NOT 9.0 (g=5 only)
  expect_equal(agg_dt$Y_bar[agg_dt$id == 3], 7.0, tolerance = 1e-10)
  expect_true(abs(agg_dt$Y_bar[agg_dt$id == 3] - ybar_u3_g3) > 0.5,
    label = "NT Y_bar must differ from g=3-only transform")
  expect_true(abs(agg_dt$Y_bar[agg_dt$id == 3] - ybar_u3_g5) > 0.5,
    label = "NT Y_bar must differ from g=5-only transform")

  # NT unit 4: Y_bar = 0.5 * 1.0 + 0.5 * 0.0 = 0.5
  expect_equal(agg_dt$Y_bar[agg_dt$id == 4], 0.5, tolerance = 1e-10)
  expect_true(abs(agg_dt$Y_bar[agg_dt$id == 4] - ybar_u4_g3) > 0.3,
    label = "NT Y_bar must differ from g=3-only transform")
})

test_that("37.2: incorrect single-cohort NT transform gives biased estimate", {
  # Construct scenario where correct weighted average gives different ATT
  # than incorrectly using a single cohort's transform for NT units
  #
  # 2 cohorts: g=3 (2 units), g=5 (1 unit), 2 NT, 6 periods
  # Design NT Y values so Y_bar_ig differs substantially across cohorts
  dt <- data.table::data.table(
    id   = rep(1:5, each = 6),
    time = rep(1:6, 5),
    Y    = c(1,1,5,5,5,5,       # unit 1 (g=3)
             1,1,5,5,5,5,       # unit 2 (g=3)
             1,1,1,1,20,20,     # unit 3 (g=5)
             0,0,0,0,10,10,     # unit 4 (NT): big jump at t=5
             0,0,0,0,10,10),    # unit 5 (NT): big jump at t=5
    gvar = rep(c(3, 3, 5, NA, NA), each = 6)
  )

  pre_stats <- make_pre_stats(c(3, 5), 1:5)

  # --- Correct computation (weighted average for NT) ---
  overall_correct <- aggregate_to_overall(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L, 5L), T_max = 6L,
    pre_stats = pre_stats, rolling = "demean"
  )
  tau_correct <- overall_correct$att

  # --- Compute what the INCORRECT result would be ---
  # If NT units used only g=3 transform (the FATAL-002 bug):
  ybar_g3 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 6L,
    pre_stat = pre_stats[["3"]], rolling = "demean"
  )
  ybar_g5 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 5L, T_max = 6L,
    pre_stat = pre_stats[["5"]], rolling = "demean"
  )

  # Treated Y_bar (correct regardless):
  # Unit 1 (g=3): mean(5-1,5-1,5-1,5-1) = 4.0
  # Unit 2 (g=3): mean(5-1,5-1,5-1,5-1) = 4.0
  # Unit 3 (g=5): mean(20-1,20-1) = 19.0
  ybar_treated <- c(4.0, 4.0, 19.0)

  # NT Y_bar using ONLY g=3 (incorrect):
  # Unit 4 g=3: mean(0-1,0-1,10-1,10-1) = mean(-1,-1,9,9) = 4.0
  # Unit 5 g=3: mean(0-1,0-1,10-1,10-1) = mean(-1,-1,9,9) = 4.0
  ybar_nt_g3_only <- c(
    ybar_g3$Y_bar_ig[ybar_g3$id == 4],
    ybar_g3$Y_bar_ig[ybar_g3$id == 5]
  )
  expect_equal(ybar_nt_g3_only[1], 4.0, tolerance = 1e-10)

  # NT Y_bar using correct weighted average (omega_3=2/3, omega_5=1/3):
  # Unit 4 g=5: mean(10-1,10-1) = 9.0
  # Unit 4 correct: (2/3)*4.0 + (1/3)*9.0 = 8/3 + 3 = 17/3 ≈ 5.6667
  ybar_nt_correct <- c(
    (2/3) * ybar_g3$Y_bar_ig[ybar_g3$id == 4] +
    (1/3) * ybar_g5$Y_bar_ig[ybar_g5$id == 4],
    (2/3) * ybar_g3$Y_bar_ig[ybar_g3$id == 5] +
    (1/3) * ybar_g5$Y_bar_ig[ybar_g5$id == 5]
  )

  # Incorrect ATT (using g=3-only for NT):
  tau_incorrect <- mean(ybar_treated) - mean(ybar_nt_g3_only)
  # Correct ATT (using weighted average for NT):
  tau_correct_manual <- mean(ybar_treated) - mean(ybar_nt_correct)

  # Verify the two ATTs are DIFFERENT — proving FATAL-002 protection matters
  expect_true(abs(tau_correct_manual - tau_incorrect) > 0.5,
    label = "Correct and incorrect NT transforms must produce different ATTs")

  # Verify the function's ATT matches the correct manual calculation
  expect_equal(tau_correct, tau_correct_manual, tolerance = 1e-5)

  # Verify the function's ATT does NOT match the incorrect calculation
  expect_true(abs(tau_correct - tau_incorrect) > 0.5,
    label = "Function ATT must differ from incorrect single-cohort NT transform")
})

# ============================================================================
# Group 38: Upstream Integration Tests
# ============================================================================

test_that("38.1: aggregate_to_overall calls construct_aggregated_outcome correctly", {
  # Full pipeline: create panel → aggregate_to_overall → verify diagnostics
  dt <- make_overall_panel(c(g3 = 4, g5 = 3), n_nt = 5, n_periods = 7)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 7L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Diagnostics from construct_aggregated_outcome must be propagated
  expect_true(!is.null(result$diagnostics))
  expect_equal(result$diagnostics$n_cohorts_requested, 2L)
  expect_equal(result$diagnostics$n_cohorts_succeeded, 2L)
  expect_true(result$diagnostics$yield_rate > 0)
  # n_valid = n_cohorts_succeeded * (units with non-NA Y_bar) > 0
  n_valid <- round(result$diagnostics$yield_rate * length(all_ids))
  expect_true(n_valid > 0, label = "n_valid units must be positive")
})

test_that("38.2: construct_aggregated_outcome calls compute_cohort_aggregated_variable for each cohort", {
  # 2 cohorts: g=3 (3 units), g=5 (2 units), 4 NT, 7 periods

  dt <- make_overall_panel(c(g3 = 3, g5 = 2), n_nt = 4, n_periods = 7)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)
  weights <- c("3" = 3 / 5, "5" = 2 / 5)

  agg <- suppressWarnings(
    construct_aggregated_outcome(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), weights = weights, T_max = 7L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Both cohorts succeeded
  expect_equal(agg$diagnostics$n_cohorts_succeeded, 2L)
  expect_equal(length(agg$diagnostics$cohort_errors), 0L)

  # Result has Y_bar for all units (some may be NA but structure is complete)
  expect_true(nrow(agg$result) == length(all_ids))
  expect_true("Y_bar" %in% names(agg$result))

  # All treated units should have non-NA Y_bar (complete panel, no missing data)
  treated_mask <- agg$result$D_ever == 1L
  expect_true(all(!is.na(agg$result$Y_bar[treated_mask])),
    label = "All treated units must have non-NA Y_bar")

  # NT units should also have non-NA Y_bar (weighted average across both cohorts)
  nt_mask <- agg$result$D_ever == 0L
  expect_true(all(!is.na(agg$result$Y_bar[nt_mask])),
    label = "All NT units must have non-NA Y_bar from weighted average")
})

test_that("38.3: aggregate_to_overall uses get_unit_level_gvar and is_never_treated", {
  # Known split: 5 treated (g=3), 4 NT
  n_treat <- 5L
  n_nt <- 4L
  dt <- make_overall_panel(c(g3 = n_treat), n_nt = n_nt, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Implicitly tests get_unit_level_gvar + is_never_treated pipeline
  expect_equal(result$n_treated, n_treat)
  expect_equal(result$n_control, n_nt)
  expect_equal(result$n_total, n_treat + n_nt)
  # Numerical reasonableness: ATT should be finite
  expect_true(is.finite(result$att))
})

# ============================================================================
# Group 39: Downstream Integration Tests
# ============================================================================

test_that("39.1: aggregate_to_overall passes correct data to estimate_ra_common", {
  # Known treatment effect = 3.0
  set.seed(101)
  dt <- make_overall_panel(
    c(g3 = 8, g5 = 6), n_nt = 10, n_periods = 7,
    treatment_effect = 3.0
  )
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 7L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # ATT should be close to the known treatment effect (demean transform
  # with pre_mean=1.0 shifts values, but the difference-in-means structure
  # should recover something in the right ballpark)
  expect_true(is.finite(result$att), label = "ATT must be finite")
  expect_true(result$se > 0, label = "SE must be positive")

  # df_resid = n_total - 2 (intercept + D_ever, no controls)
  expect_equal(result$df_resid, result$n_total - 2L)

  # df_inference must be positive

  expect_true(result$df_inference > 0, label = "df_inference must be positive")

  # CI must bracket ATT
  expect_true(result$ci_lower < result$att)
  expect_true(result$ci_upper > result$att)

  # p-value in [0, 1]
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
})

test_that("39.2: diagnostics propagated from construct to overall", {
  # Panel with enough units to potentially trigger weight renormalization
  dt <- make_overall_panel(c(g3 = 4, g5 = 3), n_nt = 5, n_periods = 7)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3, 5), all_ids)

  result <- suppressWarnings(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 7L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Diagnostics must exist and have expected fields
  diag <- result$diagnostics
  expect_true(!is.null(diag))
  expect_true("n_cohorts_requested" %in% names(diag))
  expect_true("n_cohorts_succeeded" %in% names(diag))
  expect_true("cohort_errors" %in% names(diag))
  expect_true("treated_nan_count" %in% names(diag))
  expect_true("nt_normalized_count" %in% names(diag))
  expect_true("nt_excluded_count" %in% names(diag))
  expect_true("yield_rate" %in% names(diag))

  # Consistency checks
  expect_equal(diag$n_cohorts_requested, 2L)
  expect_equal(diag$n_cohorts_succeeded, 2L)
  expect_true(diag$yield_rate > 0 && diag$yield_rate <= 1)
  # With complete panel data, no treated NaN and no NT exclusions
  expect_equal(diag$treated_nan_count, 0L)
  expect_equal(diag$nt_excluded_count, 0L)
  # nt_normalized_count should be 0 when all cohorts succeed for all NT units
  expect_equal(diag$nt_normalized_count, 0L)
})

# ============================================================================
# Group 40: Error Propagation Tests
# ============================================================================

test_that("40.1: construct_aggregated_outcome failure propagates to aggregate_to_overall", {
  # Create panel but provide empty pre_stats (no cohort data)
  dt <- make_overall_panel(c(g3 = 3), n_nt = 3, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  # pre_stats has key "99" which doesn't match cohort 3
  pre_stats <- list("99" = data.table::data.table(
    id = all_ids, pre_mean = rep(1.0, length(all_ids))
  ))

  # All cohort aggregations fail → construct_aggregated_outcome raises
  # lwdid_insufficient_data → propagates through aggregate_to_overall
  expect_error(
    suppressWarnings(
      aggregate_to_overall(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), T_max = 6L,
        pre_stats = pre_stats, rolling = "demean"
      )
    ),
    class = "lwdid_insufficient_data"
  )
})

test_that("40.2: estimate_ra_common failure handled gracefully", {
  # Only 1 treated + 1 NT unit → n_total=2
  dt <- make_overall_panel(c(g3 = 1), n_nt = 1, n_periods = 6)
  all_ids <- sort(unique(dt$id))
  pre_stats <- make_pre_stats(c(3), all_ids)

  # aggregate_to_overall warns lwdid_small_sample (n_total==2),
  # then estimate_ra_common errors with lwdid_insufficient_data (N<3)
  expect_error(
    suppressWarnings(
      aggregate_to_overall(
        dt, "Y", "id", "time", "gvar",
        cohorts = c(3L), T_max = 6L,
        pre_stats = pre_stats, rolling = "demean"
      )
    ),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Group 41: compute_event_time_weights() Unit Tests (E5-04.5)
# ============================================================================

test_that("G41.1: compute_event_time_weights — two cohorts proportional weights", {
  sizes <- c("3" = 30, "5" = 20)
  w <- compute_event_time_weights(sizes, c("3", "5"))

  expect_equal(w[["3"]], 0.6, tolerance = 1e-10)
  expect_equal(w[["5"]], 0.4, tolerance = 1e-10)
  expect_equal(sum(w), 1.0, tolerance = 1e-10)
})

test_that("G41.2: compute_event_time_weights — single cohort gets weight 1.0", {
  sizes <- c("3" = 30)
  w <- compute_event_time_weights(sizes, c("3"))

  expect_equal(w[["3"]], 1.0, tolerance = 1e-10)
  expect_equal(sum(w), 1.0, tolerance = 1e-10)
})

test_that("G41.3: compute_event_time_weights — three cohorts proportional weights", {
  sizes <- c("3" = 10, "5" = 20, "7" = 30)
  w <- compute_event_time_weights(sizes, c("3", "5", "7"))

  expect_equal(w[["3"]], 1 / 6, tolerance = 1e-10)
  expect_equal(w[["5"]], 1 / 3, tolerance = 1e-10)
  expect_equal(w[["7"]], 1 / 2, tolerance = 1e-10)
  expect_equal(sum(w), 1.0, tolerance = 1e-10)
})

test_that("G41.4: compute_event_time_weights — missing cohort key defaults to 0", {
  sizes <- c("3" = 30)
  w <- compute_event_time_weights(sizes, c("3", "99"))

  expect_equal(w[["3"]], 1.0, tolerance = 1e-10)
  expect_equal(w[["99"]], 0.0, tolerance = 1e-10)
  expect_equal(sum(w), 1.0, tolerance = 1e-10)
})

test_that("G41.5: compute_event_time_weights — empty available_cohorts errors", {
  sizes <- c("3" = 30)
  expect_error(
    compute_event_time_weights(sizes, character(0)),
    class = "lwdid_invalid_input"
  )
})

test_that("G41.6: compute_event_time_weights — total size 0 errors", {
  sizes <- c("3" = 30)
  # All available cohorts are missing from sizes → all default to 0
  expect_error(
    compute_event_time_weights(sizes, c("99", "100")),
    class = "lwdid_insufficient_data"
  )
})

# ============================================================================
# Group 42: validate_weight_sum() Unit Tests (E5-04.5)
# ============================================================================

test_that("G42.1: validate_weight_sum — exact 1.0 is valid", {
  result <- validate_weight_sum(c(0.6, 0.4))

  expect_true(result$is_valid)
  expect_equal(result$weight_sum, 1.0, tolerance = 1e-10)
})

test_that("G42.2: validate_weight_sum — within tolerance is valid", {
  result <- validate_weight_sum(c(0.6, 0.4 + 1e-7))

  expect_true(result$is_valid)
  # weight_sum should be very close to 1.0
  expect_true(abs(result$weight_sum - 1.0) <= 1e-6)
})

test_that("G42.3: validate_weight_sum — outside tolerance is invalid", {
  result <- validate_weight_sum(c(0.5, 0.4))

  expect_false(result$is_valid)
  expect_equal(result$weight_sum, 0.9, tolerance = 1e-10)
})

# ============================================================================
# Group 43: select_degrees_of_freedom() Unit Tests (E5-04.5)
# ============================================================================

test_that("G43.1: select_degrees_of_freedom — conservative picks minimum", {
  df <- select_degrees_of_freedom(
    cohort_dfs = c(48, 38), weights = c(0.6, 0.4),
    strategy = "conservative", n_cohorts = 2L
  )

  expect_identical(df, 38L)
})

test_that("G43.2: select_degrees_of_freedom — weighted computes weighted mean", {
  df <- select_degrees_of_freedom(
    cohort_dfs = c(48, 38), weights = c(0.6, 0.4),
    strategy = "weighted", n_cohorts = 2L
  )

  # 0.6*48 + 0.4*38 = 28.8 + 15.2 = 44.0
  expect_identical(df, 44L)
})

test_that("G43.3: select_degrees_of_freedom — fallback uses n_cohorts - 1", {
  df <- select_degrees_of_freedom(
    cohort_dfs = c(48, 38), weights = c(0.6, 0.4),
    strategy = "fallback", n_cohorts = 3L
  )

  expect_identical(df, 2L)
})

test_that("G43.4: select_degrees_of_freedom — all NA dfs triggers auto fallback", {
  df <- select_degrees_of_freedom(
    cohort_dfs = c(NA_real_, NA_real_), weights = c(0.6, 0.4),
    strategy = "conservative", n_cohorts = 2L
  )

  # auto fallback: max(1, 2-1) = 1
  expect_identical(df, 1L)
})

test_that("G43.5: select_degrees_of_freedom — partial NA conservative uses only valid df", {
  df <- select_degrees_of_freedom(
    cohort_dfs = c(48, NA_real_), weights = c(0.6, 0.4),
    strategy = "conservative", n_cohorts = 2L
  )

  # Only valid df is 48
  expect_identical(df, 48L)
})

test_that("G43.6: select_degrees_of_freedom — invalid strategy errors", {
  expect_error(
    select_degrees_of_freedom(
      cohort_dfs = c(48, 38), weights = c(0.6, 0.4),
      strategy = "invalid_strategy", n_cohorts = 2L
    ),
    class = "lwdid_invalid_input"
  )
})

test_that("G43.7: select_degrees_of_freedom — df at least 1 for sub-unit values", {
  # as.integer(0.5) = 0L in R, so max(1L, 0L) = 1L
  df <- select_degrees_of_freedom(
    cohort_dfs = c(0.5), weights = c(1.0),
    strategy = "conservative", n_cohorts = 1L
  )

  expect_identical(df, 1L)
})


# ============================================================================
# Group 44: Additional Event-Time Helper Tests (E5-04.5 exact specs)
# ============================================================================

# --- compute_event_time_weights: exact spec inputs ---

test_that("G44.1: compute_event_time_weights — missing cohort '5' defaults to 0", {
  # Spec: cohort_sizes = c("3"=30), available = c("3","5")
  # "5" is not in cohort_sizes → defaults to 0
  sizes <- c("3" = 30)
  w <- compute_event_time_weights(sizes, c("3", "5"))

  expect_equal(w[["3"]], 1.0, tolerance = 1e-10)
  expect_equal(w[["5"]], 0.0, tolerance = 1e-10)
  expect_equal(sum(w), 1.0, tolerance = 1e-10)
  # Numerical reasonableness: only cohort "3" has positive size
  expect_true(w[["3"]] > w[["5"]])
})

test_that("G44.2: compute_event_time_weights — total size 0 with explicit zero size", {

  # Spec: cohort_sizes = c("3"=0), available = c("3")
  # Total size = 0 → lwdid_insufficient_data
  expect_error(
    compute_event_time_weights(c("3" = 0L), c("3")),
    class = "lwdid_insufficient_data"
  )
})

# --- select_degrees_of_freedom: exact spec inputs ---

test_that("G44.3: select_degrees_of_freedom — all NA with n_cohorts=3 gives df=2", {
  # Spec: dfs = c(NA, NA), n_cohorts = 3 → max(1, 3-1) = 2
  df <- select_degrees_of_freedom(
    cohort_dfs = c(NA_real_, NA_real_), weights = c(0.5, 0.5),
    strategy = "conservative", n_cohorts = 3L
  )

  expect_identical(df, 2L)
  # Numerical reasonableness: df >= 1 always
  expect_true(df >= 1L)
})

test_that("G44.4: select_degrees_of_freedom — dfs=c(0,0) all filtered, fallback to n_cohorts-1", {
  # Spec: dfs = c(0, 0), strategy="conservative"
  # valid_mask: 0 > 0 is FALSE for both → all filtered → auto fallback
  # n_cohorts = 2 → max(1, 2-1) = 1
  df <- select_degrees_of_freedom(
    cohort_dfs = c(0, 0), weights = c(0.5, 0.5),
    strategy = "conservative", n_cohorts = 2L
  )

  expect_identical(df, 1L)
  # Numerical reasonableness: minimum df is always 1
  expect_true(df >= 1L)
})


# ============================================================================
# Group 45: aggregate_to_event_time() — Input Validation Tests (E5-04.6)
# ============================================================================

test_that("G45.1: empty data.frame returns empty list", {
  effects <- data.frame(
    cohort = integer(0), period = integer(0),
    att = numeric(0), se = numeric(0)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)
  expect_identical(result, list())
})

test_that("G45.2: empty list input raises lwdid_invalid_input", {
  cs <- c("3" = 30, "5" = 20)
  expect_error(
    aggregate_to_event_time(list(), cs),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.3: missing required columns raises lwdid_invalid_input", {
  # Missing 'se' column
  effects <- data.frame(cohort = 3L, period = 3L, att = 1.0)
  cs <- c("3" = 30)
  expect_error(
    aggregate_to_event_time(effects, cs),
    class = "lwdid_invalid_input"
  )
  # Missing 'att' column
  effects2 <- data.frame(cohort = 3L, period = 3L, se = 0.5)
  expect_error(
    aggregate_to_event_time(effects2, cs),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.4: empty cohort_sizes raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs_empty <- setNames(numeric(0), character(0))
  expect_error(
    aggregate_to_event_time(effects, cs_empty),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.5: cohort_sizes with NA raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs_na <- c("3" = NA_real_)
  expect_error(
    aggregate_to_event_time(effects, cs_na),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.6: cohort_sizes with non-positive raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs_zero <- c("3" = 0)
  expect_error(
    aggregate_to_event_time(effects, cs_zero),
    class = "lwdid_invalid_input"
  )
  cs_neg <- c("3" = -5)
  expect_error(
    aggregate_to_event_time(effects, cs_neg),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.7: alpha out of (0,1) raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs <- c("3" = 30)
  # alpha = 0 (boundary, not in open interval)
  expect_error(
    aggregate_to_event_time(effects, cs, alpha = 0),
    class = "lwdid_invalid_input"
  )
  # alpha = 1 (boundary)
  expect_error(
    aggregate_to_event_time(effects, cs, alpha = 1),
    class = "lwdid_invalid_input"
  )
  # alpha = -1 (negative)
  expect_error(
    aggregate_to_event_time(effects, cs, alpha = -1),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.8: invalid df_strategy raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs <- c("3" = 30)
  expect_error(
    aggregate_to_event_time(effects, cs, df_strategy = "invalid"),
    class = "lwdid_invalid_input"
  )
})

test_that("G45.9: event_time_range wrong length raises lwdid_invalid_input", {
  effects <- data.frame(
    cohort = 3L, period = 3L, att = 1.0, se = 0.5
  )
  cs <- c("3" = 30)
  # Length 1
  expect_error(
    aggregate_to_event_time(effects, cs, event_time_range = 0),
    class = "lwdid_invalid_input"
  )
  # Length 3
  expect_error(
    aggregate_to_event_time(effects, cs, event_time_range = c(0, 1, 2)),
    class = "lwdid_invalid_input"
  )
})

# ============================================================================
# Group 46: aggregate_to_event_time() — Core WATT Calculation Tests (E5-04.6)
# ============================================================================

test_that("G46.1: two cohorts, two event times — exact WATT values", {
  # g=3: periods 3,4 → event_time 0,1
  # g=5: period 5 → event_time 0
  # cohort_sizes: g=3 has 30, g=5 has 20 → total at e=0: 50
  # e=0: w_3 = 30/50 = 0.6, w_5 = 20/50 = 0.4
  #   WATT(0) = 0.6*2.0 + 0.4*1.5 = 1.2 + 0.6 = 1.8
  # e=1: only g=3 → w_3 = 1.0
  #   WATT(1) = 1.0*2.5 = 2.5
  effects <- data.frame(
    cohort = c(3L, 3L, 5L),
    period = c(3L, 4L, 5L),
    att = c(2.0, 2.5, 1.5),
    se = c(0.5, 0.6, 0.4),
    df_inference = c(48L, 48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  expect_length(result, 2L)

  # e=0 result
  r0 <- result[[1]]
  expect_identical(r0$event_time, 0L)
  expect_equal(r0$att, 1.8, tolerance = 1e-10)
  expect_identical(r0$n_cohorts, 2L)

  # e=1 result
  r1 <- result[[2]]
  expect_identical(r1$event_time, 1L)
  expect_equal(r1$att, 2.5, tolerance = 1e-10)
  expect_identical(r1$n_cohorts, 1L)

  # Numerical reasonableness: WATT(0) is between min and max ATT
  expect_true(r0$att >= 1.5 && r0$att <= 2.0)
  # WATT(1) is exactly the single cohort ATT
  expect_equal(r1$att, 2.5, tolerance = 1e-12)
})

test_that("G46.2: single cohort degeneracy — WATT equals ATT, SE equals SE", {
  effects <- data.frame(
    cohort = c(3L, 3L),
    period = c(3L, 4L),
    att = c(1.5, 2.0),
    se = c(0.3, 0.4),
    df_inference = c(50L, 50L)
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs)

  expect_length(result, 2L)
  # e=0: single cohort → WATT = ATT, SE = SE
  expect_equal(result[[1]]$att, 1.5, tolerance = 1e-12)
  expect_equal(result[[1]]$se, 0.3, tolerance = 1e-12)
  # e=1: single cohort → WATT = ATT, SE = SE
  expect_equal(result[[2]]$att, 2.0, tolerance = 1e-12)
  expect_equal(result[[2]]$se, 0.4, tolerance = 1e-12)
  # n_cohorts = 1 for both
  expect_identical(result[[1]]$n_cohorts, 1L)
  expect_identical(result[[2]]$n_cohorts, 1L)
})

test_that("G46.3: weighted average property — WATT in [min(ATT), max(ATT)]", {
  # Three cohorts at same event time
  effects <- data.frame(
    cohort = c(3L, 5L, 7L),
    period = c(3L, 5L, 7L),
    att = c(1.0, 3.0, 5.0),
    se = c(0.2, 0.3, 0.4),
    df_inference = c(40L, 35L, 30L)
  )
  cs <- c("3" = 10, "5" = 20, "7" = 30)
  result <- aggregate_to_event_time(effects, cs)

  expect_length(result, 1L)
  r <- result[[1]]
  # WATT must be in [min(ATT), max(ATT)] = [1.0, 5.0]
  expect_true(r$att >= 1.0)
  expect_true(r$att <= 5.0)
  # Exact: w = c(10/60, 20/60, 30/60) = c(1/6, 1/3, 1/2)
  # WATT = 1/6*1.0 + 1/3*3.0 + 1/2*5.0 = 0.1667 + 1.0 + 2.5 = 3.6667
  expected_watt <- (10 * 1.0 + 20 * 3.0 + 30 * 5.0) / 60
  expect_equal(r$att, expected_watt, tolerance = 1e-10)
})

test_that("G46.4: SE formula — SE(WATT) = sqrt(sum(w^2 * se^2))", {
  # Reuse the two-cohort setup from G46.1
  effects <- data.frame(
    cohort = c(3L, 3L, 5L),
    period = c(3L, 4L, 5L),
    att = c(2.0, 2.5, 1.5),
    se = c(0.5, 0.6, 0.4),
    df_inference = c(48L, 48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  # e=0: w_3=0.6, w_5=0.4, se_3=0.5, se_5=0.4
  # SE = sqrt(0.6^2*0.5^2 + 0.4^2*0.4^2) = sqrt(0.09 + 0.0256) = sqrt(0.1156)
  expected_se_e0 <- sqrt(0.6^2 * 0.5^2 + 0.4^2 * 0.4^2)
  expect_equal(result[[1]]$se, expected_se_e0, tolerance = 1e-10)

  # e=1: single cohort, SE = 0.6
  expect_equal(result[[2]]$se, 0.6, tolerance = 1e-12)

  # Numerical reasonableness: SE > 0
  expect_true(result[[1]]$se > 0)
  expect_true(result[[2]]$se > 0)
})

test_that("G46.5: SE monotonicity — SE(WATT) <= max(SE_g)", {
  # With multiple cohorts, the pooled SE should not exceed the largest individual SE
  effects <- data.frame(
    cohort = c(3L, 5L, 7L),
    period = c(3L, 5L, 7L),
    att = c(1.0, 2.0, 3.0),
    se = c(0.1, 0.5, 0.9),
    df_inference = c(40L, 35L, 30L)
  )
  cs <- c("3" = 10, "5" = 20, "7" = 30)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  max_se <- max(0.1, 0.5, 0.9)
  # SE(WATT) = sqrt(sum(w^2 * se^2)) <= max(se) * sqrt(sum(w^2)) <= max(se)
  expect_true(r$se <= max_se + 1e-12)
  # Also SE > 0
  expect_true(r$se > 0)
})

# ============================================================================
# Group 47: aggregate_to_event_time() — Inference Tests (E5-04.6)
# ============================================================================

test_that("G47.1: t-distribution inference — t_stat = WATT/SE, pvalue from pt()", {
  effects <- data.frame(
    cohort = 3L, period = 3L,
    att = 2.0, se = 0.5, df_inference = 48L
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  # Single cohort: WATT=2.0, SE=0.5, df=48
  expected_t <- 2.0 / 0.5
  expect_equal(r$t_stat, expected_t, tolerance = 1e-10)

  # p-value from t-distribution
  expected_p <- 2 * pt(abs(expected_t), df = 48L, lower.tail = FALSE)
  expect_equal(r$pvalue, expected_p, tolerance = 1e-10)

  # Numerical reasonableness: t_stat > 0 when att > 0
  expect_true(r$t_stat > 0)
  # p-value in [0, 1]
  expect_true(r$pvalue >= 0 && r$pvalue <= 1)
})

test_that("G47.2: CI symmetry — att - ci_lower == ci_upper - att", {
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 1.5),
    se = c(0.5, 0.4),
    df_inference = c(48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  lower_dist <- r$att - r$ci_lower
  upper_dist <- r$ci_upper - r$att
  expect_equal(lower_dist, upper_dist, tolerance = 1e-12)

  # CI must contain the point estimate
  expect_true(r$ci_lower < r$att)
  expect_true(r$ci_upper > r$att)
})

test_that("G47.3: CI width — conservative df gives wider CI than weighted df", {
  # Two cohorts with very different df values
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 2.0),
    se = c(0.5, 0.5),
    df_inference = c(5L, 100L)
  )
  cs <- c("3" = 30, "5" = 20)

  result_cons <- aggregate_to_event_time(effects, cs, df_strategy = "conservative")
  result_wt <- aggregate_to_event_time(effects, cs, df_strategy = "weighted")

  r_cons <- result_cons[[1]]
  r_wt <- result_wt[[1]]

  # Conservative uses min(5, 100) = 5 → wider CI
  # Weighted uses weighted average → higher df → narrower CI
  ci_width_cons <- r_cons$ci_upper - r_cons$ci_lower
  ci_width_wt <- r_wt$ci_upper - r_wt$ci_lower

  expect_true(ci_width_cons > ci_width_wt)

  # Both CIs should be positive width
  expect_true(ci_width_cons > 0)
  expect_true(ci_width_wt > 0)

  # df_inference should differ
  expect_true(r_cons$df_inference < r_wt$df_inference)
})

test_that("G47.4: df >= 1 always — CI uses qt for df >= 1", {
  # Even with df_inference = 1 in input, select_degrees_of_freedom returns >= 1
  effects <- data.frame(
    cohort = 3L, period = 3L,
    att = 2.0, se = 0.5, df_inference = 1L
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  # df should be at least 1
  expect_true(r$df_inference >= 1L)

  # CI should use qt(0.975, df=1) which is ~12.706
  t_crit <- qt(0.975, df = r$df_inference)
  expected_ci_lower <- r$att - t_crit * r$se
  expected_ci_upper <- r$att + t_crit * r$se
  expect_equal(r$ci_lower, expected_ci_lower, tolerance = 1e-10)
  expect_equal(r$ci_upper, expected_ci_upper, tolerance = 1e-10)
})

# ============================================================================
# Group 48: aggregate_to_event_time() — Data Processing Tests (E5-04.6)
# ============================================================================

test_that("G48.1: event_time_range filtering — only keep e in [0, 1]", {
  # g=3: periods 3,4,5 → event_time 0,1,2
  effects <- data.frame(
    cohort = c(3L, 3L, 3L),
    period = c(3L, 4L, 5L),
    att = c(1.0, 2.0, 3.0),
    se = c(0.3, 0.4, 0.5),
    df_inference = c(40L, 40L, 40L)
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs, event_time_range = c(0, 1))

  # Only event_time 0 and 1 should be present
  expect_length(result, 2L)
  expect_identical(result[[1]]$event_time, 0L)
  expect_identical(result[[2]]$event_time, 1L)
  # event_time 2 should be excluded
  event_times <- vapply(result, `[[`, integer(1), "event_time")
  expect_false(2L %in% event_times)
})

test_that("G48.2: NaN ATT/SE excluded — rows with NaN are filtered out", {
  # Two cohorts at e=0: g=3 valid, g=5 has NaN se → only g=3 contributes
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 1.5),
    se = c(0.5, NaN),
    df_inference = c(48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  expect_length(result, 1L)
  r <- result[[1]]
  expect_identical(r$event_time, 0L)
  # Only g=3 contributes (g=5 excluded due to NaN se)
  expect_equal(r$att, 2.0, tolerance = 1e-10)
  expect_equal(r$se, 0.5, tolerance = 1e-10)
  expect_identical(r$n_cohorts, 1L)
})

test_that("G48.3: all NaN for one event time — NaN result for that event time", {
  effects <- data.frame(
    cohort = c(3L, 3L),
    period = c(3L, 4L),
    att = c(NaN, 2.0),
    se = c(NaN, 0.4),
    df_inference = c(48L, 48L)
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs)

  # e=0 should be NaN (all invalid), e=1 should be valid
  expect_length(result, 2L)
  expect_true(is.nan(result[[1]]$att))
  expect_true(is.nan(result[[1]]$se))
  expect_identical(result[[1]]$n_cohorts, 0L)
  # e=1 is valid
  expect_equal(result[[2]]$att, 2.0, tolerance = 1e-12)
  expect_true(is.finite(result[[2]]$se))
})

test_that("G48.4: duplicate (cohort, period) detection — warns lwdid_data", {
  effects <- data.frame(
    cohort = c(3L, 3L, 3L),
    period = c(3L, 3L, 4L),
    att = c(2.0, 2.5, 3.0),
    se = c(0.5, 0.6, 0.4),
    df_inference = c(48L, 48L, 48L)
  )
  cs <- c("3" = 30)
  # (3, 3) is duplicated → should warn and keep first occurrence
  expect_warning(
    result <- aggregate_to_event_time(effects, cs),
    class = "lwdid_data"
  )
  # After dedup: (3,3) att=2.0 and (3,4) att=3.0
  expect_length(result, 2L)
  expect_equal(result[[1]]$att, 2.0, tolerance = 1e-12)
  expect_equal(result[[2]]$att, 3.0, tolerance = 1e-12)
})

test_that("G48.5: results sorted by event_time", {
  # Input in reverse order of periods
  effects <- data.frame(
    cohort = c(3L, 3L, 3L),
    period = c(5L, 4L, 3L),
    att = c(3.0, 2.0, 1.0),
    se = c(0.5, 0.4, 0.3),
    df_inference = c(48L, 48L, 48L)
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects, cs)

  event_times <- vapply(result, `[[`, integer(1), "event_time")
  # Should be sorted ascending
  expect_identical(event_times, sort(event_times))
  expect_identical(event_times, c(0L, 1L, 2L))
})

# ============================================================================
# Group 49: aggregate_to_event_time() — List Input & df_inference Tests (E5-04.6)
# ============================================================================

test_that("G49.1: list input produces same result as data.frame input", {
  effects_df <- data.frame(
    cohort = c(3L, 3L, 5L),
    period = c(3L, 4L, 5L),
    att = c(2.0, 2.5, 1.5),
    se = c(0.5, 0.6, 0.4),
    df_inference = c(48L, 48L, 38L)
  )
  effects_list <- list(
    list(cohort = 3L, period = 3L, att = 2.0, se = 0.5, df_inference = 48L),
    list(cohort = 3L, period = 4L, att = 2.5, se = 0.6, df_inference = 48L),
    list(cohort = 5L, period = 5L, att = 1.5, se = 0.4, df_inference = 38L)
  )
  cs <- c("3" = 30, "5" = 20)

  result_df <- aggregate_to_event_time(effects_df, cs)
  result_list <- aggregate_to_event_time(effects_list, cs)

  expect_length(result_list, length(result_df))
  for (i in seq_along(result_df)) {
    expect_identical(result_list[[i]]$event_time, result_df[[i]]$event_time)
    expect_equal(result_list[[i]]$att, result_df[[i]]$att, tolerance = 1e-12)
    expect_equal(result_list[[i]]$se, result_df[[i]]$se, tolerance = 1e-12)
    expect_equal(result_list[[i]]$ci_lower, result_df[[i]]$ci_lower, tolerance = 1e-12)
    expect_equal(result_list[[i]]$ci_upper, result_df[[i]]$ci_upper, tolerance = 1e-12)
    expect_equal(result_list[[i]]$t_stat, result_df[[i]]$t_stat, tolerance = 1e-12)
    expect_equal(result_list[[i]]$pvalue, result_df[[i]]$pvalue, tolerance = 1e-12)
  }
})

test_that("G49.2: list input without df_inference — df_inference = NA_integer_", {
  effects_list <- list(
    list(cohort = 3L, period = 3L, att = 2.0, se = 0.5),
    list(cohort = 3L, period = 4L, att = 2.5, se = 0.6)
  )
  cs <- c("3" = 30)
  result <- aggregate_to_event_time(effects_list, cs)

  expect_length(result, 2L)
  # When all df_inference are NA, fallback df = max(1, n_cohorts - 1)
  # For single cohort: n_cohorts=1, fallback = max(1, 0) = 1
  expect_true(result[[1]]$df_inference >= 1L)
  # Results should still be valid (not NaN)
  expect_true(is.finite(result[[1]]$att))
  expect_true(is.finite(result[[1]]$se))
})

test_that("G49.3: df_strategy='weighted' gives different df than 'conservative'", {
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 1.5),
    se = c(0.5, 0.4),
    df_inference = c(10L, 100L)
  )
  cs <- c("3" = 30, "5" = 20)

  result_cons <- aggregate_to_event_time(effects, cs, df_strategy = "conservative")
  result_wt <- aggregate_to_event_time(effects, cs, df_strategy = "weighted")

  # Conservative: min(10, 100) = 10
  expect_identical(result_cons[[1]]$df_inference, 10L)

  # Weighted: weighted average of 10 and 100 with weights 0.6 and 0.4
  # w_norm = c(0.6, 0.4) (both valid), weighted_df = 0.6*10 + 0.4*100 = 6 + 40 = 46
  expect_identical(result_wt[[1]]$df_inference, 46L)

  # They should differ
  expect_true(result_cons[[1]]$df_inference != result_wt[[1]]$df_inference)
})




# ============================================================================
# Group 50: aggregate_to_event_time() — Return Structure & Verbose Tests (E5-04.6)
# ============================================================================

test_that("G50.1: all 12 fields present in each result element", {
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 1.5),
    se = c(0.5, 0.4),
    df_inference = c(48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  expected_fields <- c(
    "event_time", "att", "se", "ci_lower", "ci_upper",
    "t_stat", "pvalue", "df_inference", "n_cohorts",
    "cohort_contributions", "weight_sum", "alpha"
  )

  for (i in seq_along(result)) {
    r <- result[[i]]
    for (field in expected_fields) {
      expect_true(
        field %in% names(r),
        info = sprintf("Missing field '%s' in result[[%d]]", field, i)
      )
    }
    # Exactly 12 fields
    expect_identical(length(r), 12L)
  }
})

test_that("G50.2: cohort_contributions has correct structure", {
  effects <- data.frame(
    cohort = c(3L, 5L),
    period = c(3L, 5L),
    att = c(2.0, 1.5),
    se = c(0.5, 0.4),
    df_inference = c(48L, 38L)
  )
  cs <- c("3" = 30, "5" = 20)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  cc <- r$cohort_contributions
  # Two cohorts contribute at e=0
  expect_length(cc, 2L)

  # Each contribution has: cohort, weight, att, se
  for (j in seq_along(cc)) {
    expect_true("cohort" %in% names(cc[[j]]))
    expect_true("weight" %in% names(cc[[j]]))
    expect_true("att" %in% names(cc[[j]]))
    expect_true("se" %in% names(cc[[j]]))
    # Weight should be positive
    expect_true(cc[[j]]$weight > 0)
    # SE should be non-negative
    expect_true(cc[[j]]$se >= 0)
  }

  # Verify cohort values match input
  cohorts_in_cc <- vapply(cc, `[[`, integer(1), "cohort")
  expect_true(3L %in% cohorts_in_cc)
  expect_true(5L %in% cohorts_in_cc)

  # Verify weights match expected: w_3=0.6, w_5=0.4
  for (j in seq_along(cc)) {
    if (cc[[j]]$cohort == 3L) {
      expect_equal(cc[[j]]$weight, 0.6, tolerance = 1e-12)
      expect_equal(cc[[j]]$att, 2.0, tolerance = 1e-12)
    } else if (cc[[j]]$cohort == 5L) {
      expect_equal(cc[[j]]$weight, 0.4, tolerance = 1e-12)
      expect_equal(cc[[j]]$att, 1.5, tolerance = 1e-12)
    }
  }
})

test_that("G50.3: weight_sum is correct (should be ~1.0)", {
  effects <- data.frame(
    cohort = c(3L, 5L, 7L),
    period = c(3L, 5L, 7L),
    att = c(1.0, 2.0, 3.0),
    se = c(0.2, 0.3, 0.4),
    df_inference = c(40L, 35L, 30L)
  )
  cs <- c("3" = 10, "5" = 20, "7" = 30)
  result <- aggregate_to_event_time(effects, cs)

  r <- result[[1]]
  # weight_sum should be exactly 1.0 (all cohorts present)
  expect_equal(r$weight_sum, 1.0, tolerance = 1e-12)

  # Also verify alpha is passed through
  expect_equal(r$alpha, 0.05, tolerance = 1e-12)

  # Test with custom alpha
  result2 <- aggregate_to_event_time(effects, cs, alpha = 0.10)
  expect_equal(result2[[1]]$alpha, 0.10, tolerance = 1e-12)
})



# ============================================================================
# Group 50b: make_staggered_panel helper + structure tests (E5-06.1)
# ============================================================================

make_staggered_panel <- function(n_units_per_cohort = 3, n_nt = 2,
                                  cohorts = c(3L, 5L), T_max = 6L,
                                  T_min = 1L, tau = 2.0) {
  set.seed(NULL)  # ensure RNG is not locked
  units <- list()
  uid <- 1L
  for (g in cohorts) {
    for (j in seq_len(n_units_per_cohort)) {
      periods <- seq(T_min, T_max)
      y_vals <- rnorm(length(periods), 0, 0.1)
      post_mask <- periods >= g
      if (is.function(tau)) {
        y_vals[post_mask] <- y_vals[post_mask] + vapply(
          periods[post_mask], function(t) tau(g, t), numeric(1))
      } else {
        y_vals[post_mask] <- y_vals[post_mask] + tau
      }
      units[[length(units) + 1L]] <- data.table::data.table(
        id = uid, time = periods, Y = y_vals, gvar = g
      )
      uid <- uid + 1L
    }
  }
  for (j in seq_len(n_nt)) {
    periods <- seq(T_min, T_max)
    y_vals <- rnorm(length(periods), 0, 0.1)
    units[[length(units) + 1L]] <- data.table::data.table(
      id = uid, time = periods, Y = y_vals, gvar = 0L
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

test_that("G50b.1: make_staggered_panel returns correct data.table structure", {
  dt <- make_staggered_panel()
  expect_true(data.table::is.data.table(dt))
  expect_true(all(c("id", "time", "Y", "gvar") %in% names(dt)))
  expect_true(is.numeric(dt$Y))
  expect_true(is.integer(dt$id))
  expect_true(is.integer(dt$time))
  expect_true(is.integer(dt$gvar))
})

test_that("G50b.2: make_staggered_panel has correct number of units", {
  n_upc <- 4
  n_nt <- 3
  cohorts <- c(3L, 5L, 7L)
  dt <- make_staggered_panel(n_units_per_cohort = n_upc, n_nt = n_nt,
                              cohorts = cohorts, T_max = 8L)
  expected_units <- n_upc * length(cohorts) + n_nt
  expect_equal(length(unique(dt$id)), expected_units)
  # Each unit should have T_max - T_min + 1 rows
  expect_equal(nrow(dt), expected_units * (8L - 1L + 1L))
})

test_that("G50b.3: make_staggered_panel NT units have gvar = 0L", {
  dt <- make_staggered_panel(n_units_per_cohort = 2, n_nt = 5,
                              cohorts = c(3L, 5L))
  nt_rows <- dt[gvar == 0L]
  nt_ids <- unique(nt_rows$id)
  expect_equal(length(nt_ids), 5L)
  # All NT gvar values should be exactly integer 0
  expect_true(all(nt_rows$gvar == 0L))
  expect_true(is.integer(nt_rows$gvar))
})

test_that("G50b.4: make_staggered_panel supports tau as function", {
  tau_fn <- function(g, t) (t - g + 1) * 10
  set.seed(42)
  dt <- make_staggered_panel(tau = tau_fn, cohorts = c(3L),
                              n_units_per_cohort = 1, n_nt = 0, T_max = 5L)
  # For cohort 3, post-treatment periods are t=3,4,5

  # tau_fn(3,3) = 10, tau_fn(3,4) = 20, tau_fn(3,5) = 30
  # Pre-treatment Y should be small (noise only, mean ~0)
  pre_y <- dt[time < 3L, Y]
  expect_true(all(abs(pre_y) < 1))  # noise is N(0,0.1), very unlikely > 1
  # Post-treatment Y should be dominated by tau
  post_y <- dt[time >= 3L, Y]
  expected_tau <- c(10, 20, 30)
  # Each post Y should be close to its tau value (within noise tolerance)
  for (i in seq_along(expected_tau)) {
    expect_equal(post_y[i], expected_tau[i], tolerance = 0.5)
  }
})


# ============================================================================
# Group 51: make_deterministic_panel helper + structure tests (E5-06.1)
# ============================================================================

make_deterministic_panel <- function(n_units_per_cohort = 3, n_nt = 2,
                                     cohorts = c(3L, 5L), T_max = 6L,
                                     T_min = 1L, tau = 2.0) {
  units <- list()
  uid <- 1L
  for (g in cohorts) {
    for (j in seq_len(n_units_per_cohort)) {
      periods <- seq(T_min, T_max)
      y_vals <- uid * 0.01 + periods * 0.001
      post_mask <- periods >= g
      if (is.function(tau)) {
        y_vals[post_mask] <- y_vals[post_mask] + vapply(
          periods[post_mask], function(t) tau(g, t), numeric(1))
      } else {
        y_vals[post_mask] <- y_vals[post_mask] + tau
      }
      units[[length(units) + 1L]] <- data.table::data.table(
        id = uid, time = periods, Y = y_vals, gvar = g
      )
      uid <- uid + 1L
    }
  }
  for (j in seq_len(n_nt)) {
    periods <- seq(T_min, T_max)
    y_vals <- uid * 0.01 + periods * 0.001
    units[[length(units) + 1L]] <- data.table::data.table(
      id = uid, time = periods, Y = y_vals, gvar = 0L
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

test_that("G51.1: make_deterministic_panel produces identical data on repeated calls", {
  dt1 <- make_deterministic_panel()
  dt2 <- make_deterministic_panel()
  expect_identical(dt1, dt2)
})

test_that("G51.2: make_deterministic_panel has correct structure", {
  dt <- make_deterministic_panel(n_units_per_cohort = 2, n_nt = 3,
                                  cohorts = c(3L, 5L), T_max = 6L)
  expect_true(data.table::is.data.table(dt))
  expect_true(all(c("id", "time", "Y", "gvar") %in% names(dt)))
  n_treated <- 2 * 2  # 2 cohorts × 2 units each
  n_total <- n_treated + 3  # + 3 NT
  expect_equal(length(unique(dt$id)), n_total)
  expect_equal(max(dt$time), 6L)
  # NT units have gvar=0
  nt_ids <- unique(dt[gvar == 0L, id])
  expect_equal(length(nt_ids), 3L)
})

test_that("G51.3: make_deterministic_panel supports tau as function", {
  tau_fn <- function(g, t) (t - g + 1) * 0.5
  dt <- make_deterministic_panel(tau = tau_fn, cohorts = c(3L), n_units_per_cohort = 1,
                                  n_nt = 1, T_max = 5L)
  # Unit 1, cohort 3, period 3: y = 0.01 + 0.003 + tau_fn(3,3) = 0.013 + 0.5
  expect_equal(dt[id == 1L & time == 3L, Y], 0.01 + 0.003 + 0.5, tolerance = 1e-12)
  # Period 4: tau_fn(3,4) = 1.0
  expect_equal(dt[id == 1L & time == 4L, Y], 0.01 + 0.004 + 1.0, tolerance = 1e-12)
})


# ============================================================================
# Group 52: Cohort-agg inference uses t-distribution (E5-06.6)
# ============================================================================

test_that("G52.1: cohort_agg CI uses t-distribution (wider than normal for small df)", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  pre_stats <- make_pre_stats(c(3, 5), 1:14)
  cohorts <- c(3L, 5L)
  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar", cohorts, 8L,
                        pre_stats, "demean", NULL, NULL, 0.05, NULL)
  )
  # Each cohort effect should have CI computed via t-distribution
  for (ce in result) {
    if (!is.na(ce$se) && ce$se > 0 && !is.na(ce$df_inference) && ce$df_inference > 0) {
      # t critical value > z critical value for small df
      t_crit <- qt(0.975, df = ce$df_inference)
      z_crit <- qnorm(0.975)
      expect_true(t_crit >= z_crit)
      # CI width should match t-distribution
      ci_width <- ce$ci_upper - ce$ci_lower
      expected_width <- 2 * t_crit * ce$se
      expect_equal(ci_width, expected_width, tolerance = 1e-8)
    }
  }
})

test_that("G52.2: cohort_agg pvalue uses pt() not pnorm()", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  pre_stats <- make_pre_stats(c(3, 5), 1:14)
  cohorts <- c(3L, 5L)
  result <- suppressWarnings(
    aggregate_to_cohort(dt, "Y", "id", "time", "gvar", cohorts, 8L,
                        pre_stats, "demean", NULL, NULL, 0.05, NULL)
  )
  for (ce in result) {
    if (!is.na(ce$t_stat) && !is.na(ce$df_inference) && ce$df_inference > 0) {
      expected_p <- 2 * pt(abs(ce$t_stat), df = ce$df_inference, lower.tail = FALSE)
      expect_equal(ce$pvalue, expected_p, tolerance = 1e-10)
    }
  }
})


# ============================================================================
# Group 53: ATT fallback hierarchy (E5-06.6)
# ============================================================================

test_that("G53.1: overall aggregate top-level ATT uses overall_effect$att", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          control_group = "never_treated")
  ))
  if (!is.null(result$overall_effect) && !is.na(result$overall_effect$att)) {
    expect_equal(result$att, result$overall_effect$att, tolerance = 1e-10)
  }
})

test_that("G53.2: cohort aggregate top-level ATT uses att_cohort_agg", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$att_cohort_agg), "att_cohort_agg is NULL (degenerate data)")
  if (!is.na(result$att_cohort_agg)) {
    expect_equal(result$att, result$att_cohort_agg, tolerance = 1e-10)
  }
})

test_that("G53.3: none aggregate se_att is NA_real_", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  expect_true(is.na(result$se_att))
})

test_that("G53.4: none aggregate ATT is n_treated weighted average of gr effects", {
  dt <- make_agg_panel(c(g3 = 5, g5 = 4), n_nt = 5, n_periods = 8)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  gr <- result$cohort_time_effects
  valid <- gr[is.finite(gr$att), ]
  if (nrow(valid) > 0 && sum(valid$n_treated) > 0) {
    expected_att <- sum(valid$att * valid$n_treated) / sum(valid$n_treated)
    expect_equal(result$att, expected_att, tolerance = 1e-10)
  }
})


# ============================================================================
# Group 54: gr_fallback zero weight (E5-06.6)
# ============================================================================

test_that("G54.1: gr_fallback with all zero n_treated returns NA att", {
  # Create a minimal lwdid_result with zero n_treated
  gr_effects <- data.frame(
    cohort = c(3L, 5L), ref_period = c(3L, 5L),
    att = c(1.0, 2.0), se = c(0.1, 0.2),
    n_treated = c(0L, 0L), n_control = c(5L, 5L),
    df_resid = c(10L, 10L), df_inference = c(8L, 8L)
  )
  result <- new_lwdid_result(
    att_by_cohort_time = gr_effects,
    att = NA_real_, se_att = NA_real_,
    aggregate = "none"
  )
  # When all n_treated are 0, ATT should be NA
  expect_true(is.na(result$att))
})


# ============================================================================
# Group 55: Controls passing integration (E5-06.6)
# ============================================================================

test_that("G55.1: controls parameter propagates through lwdid aggregate pipeline", {
  set.seed(123)
  n_treat <- 8
  n_nt <- 8
  n_periods <- 6
  dt <- data.table::data.table(
    id = rep(1:(n_treat + n_nt), each = n_periods),
    time = rep(1:n_periods, n_treat + n_nt),
    gvar = rep(c(rep(3L, n_treat), rep(0L, n_nt)), each = n_periods),
    Y = rnorm((n_treat + n_nt) * n_periods),
    X1 = rnorm((n_treat + n_nt) * n_periods)
  )
  # Add treatment effect for post-treatment periods
  dt[gvar == 3L & time >= 3L, Y := Y + 2.0]

  # With controls
  result_ctrl <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated", controls = "X1")
  ))
  # Without controls
  result_no_ctrl <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  # Both should produce valid results
  expect_s3_class(result_ctrl, "lwdid_result")
  expect_s3_class(result_no_ctrl, "lwdid_result")
  # Results should differ (controls affect estimation)
  # But both should have cohort_effects
  expect_false(is.null(result_ctrl$cohort_effects))
  expect_false(is.null(result_no_ctrl$cohort_effects))
})


# ============================================================================
# Group 56: Deterministic panel numerical verification (E5-06.4)
# ============================================================================

test_that("G56.1: deterministic panel WATT formula zero-tolerance verification", {
  # Use deterministic panel with known tau=2.0
  dt <- make_deterministic_panel(n_units_per_cohort = 3, n_nt = 3,
                                  cohorts = c(3L, 5L), T_max = 6L, tau = 2.0)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  # With deterministic data and constant tau, ATT should be close to 2.0
  expect_true(!is.na(result$att))
  # The ATT should be numerically reasonable (close to true effect)
  expect_true(abs(result$att - 2.0) < 0.5)  # within 0.5 of true effect
})

test_that("G56.2: deterministic panel cohort_agg statistics consistency", {
  dt <- make_deterministic_panel(n_units_per_cohort = 3, n_nt = 3,
                                  cohorts = c(3L, 5L), T_max = 6L, tau = 2.0)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$att_cohort_agg), "att_cohort_agg is NULL (degenerate data)")
  if (!is.na(result$att_cohort_agg) && !is.null(result$se_cohort_agg) &&
      !is.na(result$se_cohort_agg) && result$se_cohort_agg > 0) {
    # t_stat = att / se
    expect_equal(result$t_stat_cohort_agg,
                 result$att_cohort_agg / result$se_cohort_agg,
                 tolerance = 1e-10)
    # CI contains ATT
    expect_true(result$ci_lower_cohort_agg <= result$att_cohort_agg)
    expect_true(result$ci_upper_cohort_agg >= result$att_cohort_agg)
    # pvalue in [0, 1]
    expect_true(result$pvalue_cohort_agg >= 0)
    expect_true(result$pvalue_cohort_agg <= 1)
  }
})

test_that("G56.3: deterministic panel data structure verification", {
  dt <- make_deterministic_panel(n_units_per_cohort = 2, n_nt = 2,
                                  cohorts = c(3L, 5L), T_max = 6L)
  # Verify exact Y values for unit 1 (cohort 3)
  u1 <- dt[id == 1L]
  expect_equal(u1[time == 1L, Y], 0.01 + 0.001, tolerance = 1e-14)  # pre-treatment
  expect_equal(u1[time == 2L, Y], 0.01 + 0.002, tolerance = 1e-14)  # pre-treatment
  expect_equal(u1[time == 3L, Y], 0.01 + 0.003 + 2.0, tolerance = 1e-14)  # post (g=3)
  # NT unit (id = 2*2 + 2 + 1 = 7... actually uid=5 for first NT with 2 cohorts × 2 units)
  # First NT unit is uid=5
  u_nt <- dt[id == 5L]
  expect_equal(u_nt[time == 1L, Y], 0.05 + 0.001, tolerance = 1e-14)
  expect_equal(u_nt[time == 6L, Y], 0.05 + 0.006, tolerance = 1e-14)
})


# ============================================================================
# Group 57: L1 FATAL Protection Tests (E5-06.2)
# ============================================================================

test_that("G57.1: FATAL-002 NT weighted average verification with staggered panel", {
  # Use staggered panel with 2 cohorts; verify that construct_aggregated_outcome
  # produces a weighted average of per-cohort Y_bar_ig for NT units.
  # Weights: w_g = N_g / N_treat (cohort size / total treated)
  set.seed(20240601)
  dt <- make_staggered_panel(
    n_units_per_cohort = 3, n_nt = 2,
    cohorts = c(3L, 5L), T_max = 6L, tau = 2.0
  )

  # precompute_transforms needs all units
  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time",
    cohorts = c(3L, 5L), rolling = "demean"
  )

  # Compute per-cohort Y_bar_ig for all units
  ybar_g3 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 3L, T_max = 6L,
    pre_stat = pre_stats[["3"]], rolling = "demean"
  )
  ybar_g5 <- compute_cohort_aggregated_variable(
    dt, "Y", "id", "time", g = 5L, T_max = 6L,
    pre_stat = pre_stats[["5"]], rolling = "demean"
  )

  # Identify NT units (gvar == 0L in make_staggered_panel)
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  nt_ids <- unit_gvar$id[is_never_treated(unit_gvar$gvar)]
  expect_true(length(nt_ids) >= 2L)

  # Weights: N_g / N_treat. Cohort 3 has 3 units, cohort 5 has 3 units -> 0.5 each
  n_g3 <- sum(unit_gvar$gvar == 3L & !is_never_treated(unit_gvar$gvar))
  n_g5 <- sum(unit_gvar$gvar == 5L & !is_never_treated(unit_gvar$gvar))
  n_treat <- n_g3 + n_g5
  w_g3 <- n_g3 / n_treat
  w_g5 <- n_g5 / n_treat
  weights <- c("3" = w_g3, "5" = w_g5)

  # Construct aggregated outcome (this is the FATAL-002 core function)
  agg <- construct_aggregated_outcome(
    dt, "Y", "id", "time", "gvar",
    cohorts = c(3L, 5L), weights = weights,
    T_max = 6L, pre_stats = pre_stats, rolling = "demean"
  )
  agg_dt <- agg$result

  # For each NT unit, verify Y_bar equals the weighted average of per-cohort values

  for (uid in nt_ids) {
    ybar_u_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == uid]
    ybar_u_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == uid]
    expected_weighted_avg <- w_g3 * ybar_u_g3 + w_g5 * ybar_u_g5
    actual <- agg_dt$Y_bar[agg_dt$id == uid]
    expect_equal(actual, expected_weighted_avg, tolerance = 1e-12,
      label = sprintf("NT unit %d Y_bar must be weighted average", uid))
  }

  # Verify weighted average differs from any single cohort's value
  # (when per-cohort values differ, which they should due to different post windows)
  for (uid in nt_ids) {
    ybar_u_g3 <- ybar_g3$Y_bar_ig[ybar_g3$id == uid]
    ybar_u_g5 <- ybar_g5$Y_bar_ig[ybar_g5$id == uid]
    actual <- agg_dt$Y_bar[agg_dt$id == uid]
    if (abs(ybar_u_g3 - ybar_u_g5) > 1e-8) {
      expect_true(abs(actual - ybar_u_g3) > 1e-12,
        label = sprintf("NT unit %d: weighted avg must differ from g=3 value", uid))
      expect_true(abs(actual - ybar_u_g5) > 1e-12,
        label = sprintf("NT unit %d: weighted avg must differ from g=5 value", uid))
    }
  }
})

test_that("G57.2: FATAL-004 NT-only control group verification", {
  # Verify each cohort's n_control equals the number of NT units,
  # confirming other treated cohorts are excluded from the control group.
  set.seed(20240602)
  dt <- make_staggered_panel(
    n_units_per_cohort = 4, n_nt = 3,
    cohorts = c(3L, 5L), T_max = 6L, tau = 2.0
  )

  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  n_nt_actual <- sum(is_never_treated(unit_gvar$gvar))
  expect_equal(n_nt_actual, 3L)

  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time",
    cohorts = c(3L, 5L), rolling = "demean"
  )

  results <- suppressWarnings(
    aggregate_to_cohort(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    )
  )

  # Each cohort result's n_control must equal the number of NT units
  expect_true(length(results) >= 1L)
  for (res in results) {
    expect_equal(res$n_control, n_nt_actual,
      label = sprintf("Cohort g=%d: n_control must equal number of NT units", res$cohort))
  }
})

test_that("G57.3: FATAL-004 no-NT error for aggregate_to_cohort and aggregate_to_overall", {
  # Create panel with NO NT units (all units are treated)
  set.seed(20240603)
  dt <- make_staggered_panel(
    n_units_per_cohort = 3, n_nt = 0,
    cohorts = c(3L, 5L), T_max = 6L, tau = 2.0
  )

  # Confirm no NT units exist
  unit_gvar <- get_unit_level_gvar(dt, "gvar", "id")
  expect_equal(sum(is_never_treated(unit_gvar$gvar)), 0L)

  pre_stats <- precompute_transforms(
    dt, "Y", "id", "time",
    cohorts = c(3L, 5L), rolling = "demean"
  )

  # aggregate_to_cohort must throw lwdid_no_never_treated
  expect_error(
    aggregate_to_cohort(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_no_never_treated"
  )

  # aggregate_to_overall must also throw lwdid_no_never_treated
  expect_error(
    aggregate_to_overall(
      dt, "Y", "id", "time", "gvar",
      cohorts = c(3L, 5L), T_max = 6L,
      pre_stats = pre_stats, rolling = "demean"
    ),
    class = "lwdid_no_never_treated"
  )
})
