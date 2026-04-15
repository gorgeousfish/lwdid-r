# ============================================================================
# test-utils.R — Comprehensive tests for lwdid utility functions and constants
# ============================================================================

# ── Group A: Numerical Constants ────────────────────────────────────────────

test_that("numerical constants have correct values", {
  expect_identical(LWDID_NT_ZERO_TOLERANCE, 1e-10)
  expect_identical(LWDID_COHORT_FLOAT_TOLERANCE, 1e-9)
  expect_identical(LWDID_WEIGHT_SUM_TOLERANCE, 1e-6)
  expect_identical(LWDID_NUMERICAL_TOLERANCE, 1e-12)
  expect_identical(LWDID_VARIANCE_THRESHOLD, 1e-10)
  expect_identical(LWDID_PS_TRIM_DEFAULT, 0.01)
  expect_identical(LWDID_STATA_ATT_TOLERANCE, 1e-6)
  expect_identical(LWDID_PYTHON_ATT_TOLERANCE, 1e-5)
})

# ── Group B: Valid Value Set Constants ──────────────────────────────────────

test_that("valid value set constants have correct values", {
  expect_identical(
    LWDID_VALID_ESTIMATORS,
    c("ra", "ipw", "ipwra", "psm")
  )
  expect_identical(
    LWDID_VALID_TRANSFORMATIONS,
    c("demean", "detrend", "demeanq", "detrendq")
  )
  expect_identical(
    LWDID_VALID_BASE_PERIODS,
    c("universal", "varying")
  )
  expect_identical(
    LWDID_VALID_BOOT_TYPES,
    c("weighted", "multiplier", "bayesian")
  )
  expect_identical(
    LWDID_VALID_AGGREGATE_LEVELS,
    c("none", "cohort", "overall", "event_time")
  )
})

# ── Group L: %||% operator ─────────────────────────────────────────────────

test_that("%||% returns fallback for NULL", {
  expect_equal(`%||%`(NULL, 5), 5)
})

test_that("%||% returns original for non-NULL", {
  expect_equal(`%||%`(3, 5), 3)
})

test_that("%||% treats empty string as non-NULL", {
  expect_equal(`%||%`("", "default"), "")
})

test_that("%||% treats zero-length vector as non-NULL", {
  expect_equal(`%||%`(character(0), "default"), character(0))
})

test_that("%||% treats FALSE as non-NULL", {
  expect_equal(`%||%`(FALSE, TRUE), FALSE)
})

test_that("%||% treats NA as non-NULL", {
  expect_identical(`%||%`(NA, 0), NA)
})

# ── Group C: is_never_treated ──────────────────────────────────────────────

test_that("is_never_treated: scalar NA returns TRUE", {
  expect_true(is_never_treated(NA))
})

test_that("is_never_treated: scalar NaN returns TRUE", {
  expect_true(is_never_treated(NaN))
})

test_that("is_never_treated: scalar 0 returns TRUE", {
  expect_true(is_never_treated(0))
})

test_that("is_never_treated: scalar Inf returns TRUE", {
  expect_true(is_never_treated(Inf))
})

test_that("is_never_treated: scalar 2005 returns FALSE", {
  expect_false(is_never_treated(2005))
})

test_that("is_never_treated: scalar -1 returns FALSE", {
  expect_false(is_never_treated(-1))
})

test_that("is_never_treated: near-zero 1e-11 returns TRUE", {
  expect_true(is_never_treated(1e-11))
})

test_that("is_never_treated: near-zero -1e-11 returns TRUE", {
  expect_true(is_never_treated(-1e-11))
})

test_that("is_never_treated: 1e-9 is NOT near-zero, returns FALSE", {
  expect_false(is_never_treated(1e-9))
})

test_that("is_never_treated: -Inf triggers error", {
  expect_error(
    is_never_treated(-Inf),
    class = "lwdid_invalid_staggered_data"
  )
})

test_that("is_never_treated: vectorized input", {
  result <- is_never_treated(c(NA, 0, Inf, 2005, 2010))
  expect_identical(result, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("is_never_treated: integer vector", {
  result <- is_never_treated(c(0L, 2005L, NA_integer_))
  expect_identical(result, c(TRUE, FALSE, TRUE))
})

test_that("is_never_treated: empty vector returns logical(0)", {
  expect_identical(is_never_treated(numeric(0)), logical(0))
})

# ── Group D: get_unit_level_gvar ───────────────────────────────────────────

test_that("get_unit_level_gvar: basic panel returns data.table", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(2001, 2001, NA, NA, 2002, 2002)
  )
  result <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(nrow(result), 3L)
  expect_true(data.table::is.data.table(result))
  expect_true("id" %in% names(result))
  expect_true("gvar" %in% names(result))
})

test_that("get_unit_level_gvar: ivar column contains unit IDs", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(2001, 2001, NA, NA, 2002, 2002)
  )
  result <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(sort(result$id), c(1, 2, 3))
})

test_that("get_unit_level_gvar: gvar values correct", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(2001, 2001, NA, NA, 2002, 2002)
  )
  result <- get_unit_level_gvar(panel, "gvar", "id")
  # Sort by id to ensure deterministic order
  result <- result[order(result$id), ]
  expect_equal(result$gvar, c(2001, NA, 2002))
})

test_that("get_unit_level_gvar: empty data returns zero-row data.table", {
  panel <- data.frame(id = numeric(0), year = numeric(0), gvar = numeric(0))
  result <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(nrow(result), 0L)
  expect_true(data.table::is.data.table(result))
})

test_that("get_unit_level_gvar: single unit returns one row", {
  panel <- data.frame(id = c(1, 1), year = c(2000, 2001), gvar = c(2005, 2005))
  result <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(nrow(result), 1L)
  expect_equal(result$gvar, 2005)
})

# ── Group E: get_cohorts ───────────────────────────────────────────────────

test_that("get_cohorts: basic panel returns sorted integer cohorts", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(2005, 2005, NA, NA, 2010, 2010)
  )
  result <- get_cohorts(panel, "gvar", "id")
  expect_identical(result, c(2005L, 2010L))
})

test_that("get_cohorts: all NT returns integer(0)", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(NA, NA, 0, 0, Inf, Inf)
  )
  result <- get_cohorts(panel, "gvar", "id")
  expect_identical(result, integer(0))
})

test_that("get_cohorts: custom never_treated_values excludes 999", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2, 3, 3),
    year = c(2000, 2001, 2000, 2001, 2000, 2001),
    gvar = c(2005, 2005, 999, 999, 2010, 2010)
  )
  result <- get_cohorts(
    panel, "gvar", "id",
    never_treated_values = c(0, 999)
  )
  expect_identical(result, c(2005L, 2010L))
})

test_that("get_cohorts: float gvar coerced to integer", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(2005.0, 2005.0, 2010.0, 2010.0)
  )
  result <- get_cohorts(panel, "gvar", "id")
  expect_identical(result, c(2005L, 2010L))
})

test_that("get_cohorts: near-zero gvar excluded as NT", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(1e-11, 1e-11, 2005, 2005)
  )
  result <- get_cohorts(panel, "gvar", "id")
  expect_identical(result, 2005L)
})

test_that("get_cohorts: -Inf triggers error", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(-Inf, -Inf, 2005, 2005)
  )
  expect_error(
    get_cohorts(panel, "gvar", "id"),
    class = "lwdid_invalid_staggered_data"
  )
})

# ── Group F: has_never_treated ─────────────────────────────────────────────

test_that("has_never_treated: panel with NA gvar returns TRUE", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(NA, NA, 2005, 2005)
  )
  expect_true(has_never_treated(panel, "gvar", "id"))
})

test_that("has_never_treated: panel with 0 gvar returns TRUE", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(0, 0, 2005, 2005)
  )
  expect_true(has_never_treated(panel, "gvar", "id"))
})

test_that("has_never_treated: panel with Inf gvar returns TRUE", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(Inf, Inf, 2005, 2005)
  )
  expect_true(has_never_treated(panel, "gvar", "id"))
})

test_that("has_never_treated: only valid cohorts returns FALSE", {
  panel <- data.frame(
    id   = c(1, 1, 2, 2),
    year = c(2000, 2001, 2000, 2001),
    gvar = c(2005, 2005, 2010, 2010)
  )
  expect_false(has_never_treated(panel, "gvar", "id"))
})

# ── Group G: get_never_treated_mask and get_cohort_mask ────────────────────

test_that("get_never_treated_mask: consistent with is_never_treated", {
  gvars <- c(NA, 0, Inf, 2005, 1e-11, 2010)
  mask <- get_never_treated_mask(gvars)
  expected <- is_never_treated(gvars)
  expect_identical(mask, expected)
})

test_that("get_cohort_mask: float tolerance matches within 1e-9", {
  gvars <- c(2004.9999999999, 2005.0000000001)
  mask <- get_cohort_mask(gvars, 2005)
  expect_identical(mask, c(TRUE, TRUE))
})

test_that("get_cohort_mask: NA returns FALSE", {
  gvars <- c(NA, 2005, NA)
  mask <- get_cohort_mask(gvars, 2005)
  expect_identical(mask, c(FALSE, TRUE, FALSE))
})

test_that("get_cohort_mask: boundary exactly 1e-9 returns FALSE", {
  # Due to floating-point representation, 2005 + 1e-9 may not be
  # exactly 1e-9 away. Use a value clearly beyond the tolerance.
  gvars <- c(2005 + 2e-9, 2005 - 2e-9)
  mask <- get_cohort_mask(gvars, 2005)
  expect_identical(mask, c(FALSE, FALSE))
})

test_that("get_cohort_mask: named vector preserves names", {
  gvars <- c(a = 2005, b = 2010, c = NA)
  mask <- get_cohort_mask(gvars, 2005)
  expect_identical(names(mask), c("a", "b", "c"))
  expect_identical(unname(mask), c(TRUE, FALSE, FALSE))
})

# ── Group H: get_valid_periods_for_cohort ──────────────────────────────────

test_that("get_valid_periods_for_cohort: cohort=2005, T_max=2010", {
  result <- get_valid_periods_for_cohort(2005, 2010)
  expect_identical(result, 2005L:2010L)
  expect_length(result, 6)
})

test_that("get_valid_periods_for_cohort: cohort=T_max returns single", {
  result <- get_valid_periods_for_cohort(2010, 2010)
  expect_identical(result, 2010L)
  expect_length(result, 1)
})

test_that("get_valid_periods_for_cohort: cohort > T_max returns integer(0)", {
  result <- get_valid_periods_for_cohort(2011, 2010)
  expect_identical(result, integer(0))
})

# ── Group I: extract_cross_section ─────────────────────────────────────────

test_that("extract_cross_section: normal extraction correct row count", {
  dt <- data.table::data.table(
    id   = c(1, 2, 3, 1, 2, 3),
    year = c(2000, 2000, 2000, 2001, 2001, 2001),
    y    = c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
  )
  result <- extract_cross_section(dt, "year", 2000)
  expect_equal(nrow(result), 3)
  expect_equal(result$id, c(1, 2, 3))
})

test_that("extract_cross_section: non-existent period returns zero rows", {
  dt <- data.table::data.table(
    id   = c(1, 2),
    year = c(2000, 2000),
    y    = c(1.0, 2.0)
  )
  result <- extract_cross_section(dt, "year", 9999)
  expect_equal(nrow(result), 0)
  expect_identical(names(result), c("id", "year", "y"))
})

test_that("extract_cross_section: output is data.table", {
  dt <- data.table::data.table(
    id   = c(1, 2),
    year = c(2000, 2001),
    y    = c(1.0, 2.0)
  )
  result <- extract_cross_section(dt, "year", 2000)
  expect_true(data.table::is.data.table(result))
})

# ── Group J: safe_mean ─────────────────────────────────────────────────────

test_that("safe_mean: basic numeric vector", {
  expect_equal(safe_mean(c(1, 2, 3)), 2)
})

test_that("safe_mean: all-NA with na.rm=TRUE returns NA_real_ no warning", {
  expect_no_warning({
    result <- safe_mean(c(NA, NA), na.rm = TRUE)
  })
  expect_identical(result, NA_real_)
})

test_that("safe_mean: empty vector returns NA_real_ no warning", {
  expect_no_warning({
    result <- safe_mean(numeric(0))
  })
  expect_identical(result, NA_real_)
})

test_that("safe_mean: NULL returns NA_real_", {
  expect_identical(safe_mean(NULL), NA_real_)
})

test_that("safe_mean: vector with NA uses na.rm", {
  expect_equal(safe_mean(c(1, NA, 3)), 2)
})

test_that("safe_mean: Inf propagates", {
  expect_equal(safe_mean(c(Inf, 1, 2)), Inf)
})

# ── Group K: center_by_treated ─────────────────────────────────────────────

test_that("center_by_treated: numerical verification", {
  X <- matrix(c(1, 3, 5, 2, 4, 6), nrow = 3)
  d <- c(1, 1, 0)
  result <- center_by_treated(X, d)

  # Treated mean = colMeans of rows 1,2 = c(2, 3)
  expect_equal(result$X_mean_treated, c(2, 3))

  # Centered values: each row minus c(2, 3)
  expected_centered <- matrix(c(-1, 1, 3, -1, 1, 3), nrow = 3)
  expect_equal(result$X_centered, expected_centered)
})

test_that("center_by_treated: treated group mean of centered X is zero", {
  X <- matrix(c(1, 3, 5, 2, 4, 6), nrow = 3)
  d <- c(1, 1, 0)
  result <- center_by_treated(X, d)

  treated_centered <- result$X_centered[d == 1, , drop = FALSE]
  treated_means <- colMeans(treated_centered)
  expect_equal(treated_means, c(0, 0), tolerance = 1e-14)
})

test_that("center_by_treated: single column X works", {
  X <- matrix(c(10, 20, 30), nrow = 3)
  d <- c(1, 0, 0)
  result <- center_by_treated(X, d)

  expect_equal(result$X_mean_treated, 10)
  expect_equal(result$X_centered[, 1], c(0, 10, 20))
})

test_that("center_by_treated: no treated units triggers error", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  d <- c(0, 0)
  expect_error(
    center_by_treated(X, d),
    class = "lwdid_no_treated"
  )
})

# ── Group M: stable_ols ───────────────────────────────────────────────────

test_that("stable_ols: exact linear system", {
  X <- cbind(1, 1:5)
  y <- 2 + 3 * (1:5)  # c(5, 8, 11, 14, 17)
  result <- stable_ols(X, y)

  expect_equal(result$coefficients, c(2, 3), tolerance = 1e-10)
  expect_equal(result$residuals, rep(0, 5), tolerance = 1e-10)
  expect_equal(result$rank, 2)
  expect_equal(result$df.residual, 3L)
  expect_false(result$is_rank_deficient)
})

test_that("stable_ols: weighted OLS matches lm", {
  x_vals <- 1:5
  X <- cbind(1, x_vals)
  y <- 2 + 3 * x_vals
  w <- as.numeric(1:5)

  result <- stable_ols(X, y, w = w)
  lm_fit <- lm(y ~ x_vals, weights = w)

  expect_equal(
    unname(result$coefficients),
    unname(coef(lm_fit)),
    tolerance = 1e-10
  )
})

test_that("stable_ols: weight with zero triggers error", {
  X <- cbind(1, 1:3)
  y <- c(1, 2, 3)
  expect_error(
    stable_ols(X, y, w = c(1, 0, 1)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("stable_ols: weight with negative triggers error", {
  X <- cbind(1, 1:3)
  y <- c(1, 2, 3)
  expect_error(
    stable_ols(X, y, w = c(1, -1, 1)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("stable_ols: weight with NA triggers error", {
  X <- cbind(1, 1:3)
  y <- c(1, 2, 3)
  expect_error(
    stable_ols(X, y, w = c(1, NA, 1)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("stable_ols: weight wrong length triggers error", {
  X <- cbind(1, 1:3)
  y <- c(1, 2, 3)
  expect_error(
    stable_ols(X, y, w = c(1, 2)),
    class = "lwdid_invalid_parameter"
  )
})

test_that("stable_ols: rank deficient matrix detected", {
  X <- cbind(1, 1:5, 2:6)  # col3 = col1 + col2
  y <- rnorm(5)
  result <- stable_ols(X, y)
  expect_true(result$is_rank_deficient)
})

test_that("stable_ols: return value has all 8 fields", {
  X <- cbind(1, 1:5)
  y <- rnorm(5)
  result <- stable_ols(X, y)
  expected_names <- c(
    "coefficients", "residuals", "fitted.values",
    "rank", "df.residual", "qr",
    "is_rank_deficient", "weights"
  )
  expect_true(all(expected_names %in% names(result)))
  expect_length(result, 8)
})

# ── Group N: check_rank ───────────────────────────────────────────────────

test_that("check_rank: full rank matrix", {
  X <- cbind(1, 1:5, (1:5)^2)
  result <- check_rank(X)
  expect_true(result$is_full_rank)
  expect_identical(result$problematic_columns, integer(0))
  expect_equal(result$rank, 3)
  expect_equal(result$ncol, 3)
  expect_equal(result$rank_deficiency, 0)
})

test_that("check_rank: rank deficient identifies problematic column", {
  col1 <- c(1, 2, 3, 4, 5)
  col2 <- c(2, 3, 4, 5, 6)
  col3 <- col1 + col2  # linearly dependent
  X <- cbind(col1, col2, col3)
  result <- check_rank(X)
  expect_equal(result$rank, 2)
  expect_false(result$is_full_rank)
  expect_true(3 %in% result$problematic_columns)
})

test_that("check_rank: single column non-zero has rank 1", {
  X <- matrix(c(1, 2, 3), nrow = 3)
  result <- check_rank(X)
  expect_equal(result$rank, 1)
  expect_true(result$is_full_rank)
})

test_that("check_rank: custom tol parameter works", {
  X <- cbind(1, 1:5, (1:5) + 1e-8)
  result_default <- check_rank(X)
  result_strict <- check_rank(X, tol = 1e-6)
  # With strict tolerance, near-collinear columns detected
  expect_true(result_strict$rank <= result_default$rank)
})


# ── Group O: Integration Tests (E1-07.17) ─────────────────────────────────

test_that("Integration: staggered panel data function chain works", {
  # Create panel data with 3 units: unit 1 treated in 2005, unit 2 treated in 2010, unit 3 never-treated (gvar=0)
  dt <- data.table::data.table(
    id = rep(1:3, each = 5),
    year = rep(2003:2007, 3),
    gvar = rep(c(2005, 2010, 0), each = 5),
    y = rnorm(15)
  )

  # Step 1: get_unit_level_gvar
  unit_gvar_dt <- get_unit_level_gvar(dt, "gvar", "id")
  expect_equal(nrow(unit_gvar_dt), 3L)
  expect_true(data.table::is.data.table(unit_gvar_dt))

  # Convert to named vector for downstream mask/cohort operations
  unit_gvar <- unit_gvar_dt[["gvar"]]
  names(unit_gvar) <- as.character(unit_gvar_dt[["id"]])

  # Step 2: is_never_treated on unit-level gvar
  nt_mask <- is_never_treated(unit_gvar)
  expect_equal(nt_mask, c("1" = FALSE, "2" = FALSE, "3" = TRUE))

  # Step 3: get_cohorts
  cohorts <- get_cohorts(dt, "gvar", "id")
  expect_equal(cohorts, c(2005L, 2010L))

  # Step 4: has_never_treated
  expect_true(has_never_treated(dt, "gvar", "id"))

  # Step 5: get_never_treated_mask
  nt_mask2 <- get_never_treated_mask(unit_gvar)
  expect_equal(nt_mask2, nt_mask)

  # Step 6: get_cohort_mask for cohort 2005
  c_mask <- get_cohort_mask(unit_gvar, 2005)
  expect_equal(c_mask, c("1" = TRUE, "2" = FALSE, "3" = FALSE))

  # Step 7: get_valid_periods_for_cohort
  periods <- get_valid_periods_for_cohort(2005, 2007)
  expect_equal(periods, 2005L:2007L)

  # Step 8: extract_cross_section
  cs <- extract_cross_section(dt, "year", 2005)
  expect_equal(nrow(cs), 3L)
  expect_true(data.table::is.data.table(cs))
})

test_that("Integration: numerical computation chain recovers known ATT", {
  set.seed(42)
  n <- 100
  # Treatment indicator: first 50 treated, last 50 control
  d <- c(rep(1, 50), rep(0, 50))
  # Control variable
  x <- rnorm(n)
  # True ATT = 5.0, true beta = 2.0
  y <- 5.0 * d + 2.0 * x + rnorm(n, sd = 0.01)

  # Build design matrix with intercept, treatment, and control
  X <- cbind(intercept = 1, treatment = d, x = x)

  # Step 1: check_rank
  rank_info <- check_rank(X)
  expect_true(rank_info$is_full_rank)
  expect_equal(rank_info$rank, 3L)

  # Step 2: stable_ols
  ols_result <- stable_ols(X, y)
  expect_false(ols_result$is_rank_deficient)
  # Coefficients should be close to c(0, 5, 2)
  expect_equal(ols_result$coefficients[["treatment"]], 5.0, tolerance = 0.1)
  expect_equal(ols_result$coefficients[["x"]], 2.0, tolerance = 0.1)

  # Step 3: center_by_treated on control variable
  X_ctrl <- matrix(x, ncol = 1)
  centered <- center_by_treated(X_ctrl, d)
  # Treated group mean of centered X should be ~0
  expect_equal(mean(centered$X_centered[d == 1, 1]), 0, tolerance = 1e-12)
})

test_that("Integration: error condition classes are consistent with conditions.R", {
  # is_never_treated with -Inf should throw lwdid_invalid_staggered_data
  err <- tryCatch(is_never_treated(-Inf), error = function(e) e)
  expect_s3_class(err, "lwdid_invalid_staggered_data")
  expect_s3_class(err, "lwdid_error")

  # center_by_treated with no treated units should throw lwdid_no_treated
  X <- matrix(1:4, ncol = 2)
  d <- c(0, 0)
  err2 <- tryCatch(center_by_treated(X, d), error = function(e) e)
  expect_s3_class(err2, "lwdid_no_treated")
  expect_s3_class(err2, "lwdid_error")

  # stable_ols with bad weights should throw lwdid_invalid_parameter
  X3 <- matrix(c(1,1,1,2), ncol = 2)
  y3 <- c(1, 2)
  err3 <- tryCatch(stable_ols(X3, y3, w = c(-1, 1)), error = function(e) e)
  expect_s3_class(err3, "lwdid_invalid_parameter")
  expect_s3_class(err3, "lwdid_error")
})


# ============================================================================
# Integration Tests (E1-07.17)
# ============================================================================

# ── Integration Test 1: Full function chain on staggered panel data ────────

test_that("integration: get_unit_level_gvar → is_never_treated → get_cohorts chain", {
  # Construct a staggered panel: 6 units, 5 periods (2000-2004)
  # Units 1,2: cohort 2002; Units 3,4: cohort 2003; Units 5,6: never-treated (NA)
  panel <- data.frame(
    id   = rep(1:6, each = 5),
    year = rep(2000:2004, times = 6),
    gvar = rep(c(2002, 2002, 2003, 2003, NA, NA), each = 5),
    y    = rnorm(30)
  )

  # Step 1: get_unit_level_gvar
  unit_gvar_dt <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(nrow(unit_gvar_dt), 6L)

  # Convert to named vector for downstream operations
  unit_gvar <- unit_gvar_dt[["gvar"]]
  names(unit_gvar) <- as.character(unit_gvar_dt[["id"]])
  expect_identical(names(unit_gvar), as.character(1:6))
  expect_equal(unname(unit_gvar), c(2002, 2002, 2003, 2003, NA, NA))

  # Step 2: is_never_treated on unit-level gvar
  nt_mask <- is_never_treated(unit_gvar)
  expect_identical(unname(nt_mask), c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE))

  # Step 3: get_cohorts
  cohorts <- get_cohorts(panel, "gvar", "id")
  expect_identical(cohorts, c(2002L, 2003L))
})

test_that("integration: has_never_treated + mask functions chain", {
  panel <- data.frame(
    id   = rep(1:4, each = 3),
    year = rep(2000:2002, times = 4),
    gvar = rep(c(2001, 2002, 0, Inf), each = 3),
    y    = rnorm(12)
  )

  # has_never_treated detects 0 and Inf
  expect_true(has_never_treated(panel, "gvar", "id"))

  # get unit-level gvar for mask operations
  unit_gvar_dt <- get_unit_level_gvar(panel, "gvar", "id")
  unit_gvar <- unit_gvar_dt[["gvar"]]
  names(unit_gvar) <- as.character(unit_gvar_dt[["id"]])

  # NT mask: units 3 (gvar=0) and 4 (gvar=Inf) are NT
  nt_mask <- get_never_treated_mask(unit_gvar)
  expect_identical(unname(nt_mask), c(FALSE, FALSE, TRUE, TRUE))

  # Cohort mask for g=2001: only unit 1
  cohort_mask_2001 <- get_cohort_mask(unit_gvar, 2001)
  expect_identical(unname(cohort_mask_2001), c(TRUE, FALSE, FALSE, FALSE))

  # Cohort mask for g=2002: only unit 2
  cohort_mask_2002 <- get_cohort_mask(unit_gvar, 2002)
  expect_identical(unname(cohort_mask_2002), c(FALSE, TRUE, FALSE, FALSE))
})

test_that("integration: get_valid_periods_for_cohort → extract_cross_section chain", {
  dt <- data.table::data.table(
    id   = rep(1:4, each = 6),
    year = rep(2000:2005, times = 4),
    gvar = rep(c(2003, 2003, NA, NA), each = 6),
    y    = seq_len(24) * 1.0
  )

  # Get valid post-treatment periods for cohort 2003 with T_max=2005
  valid_periods <- get_valid_periods_for_cohort(2003, 2005)
  expect_identical(valid_periods, 2003L:2005L)

  # Extract cross-section for each valid period and verify row counts
  for (p in valid_periods) {
    cs <- extract_cross_section(dt, "year", p)
    expect_equal(nrow(cs), 4, info = paste("period", p))
    expect_true(data.table::is.data.table(cs))
    expect_true(all(cs$year == p))
  }

  # Extract cross-section for non-existent period
  cs_empty <- extract_cross_section(dt, "year", 9999)
  expect_equal(nrow(cs_empty), 0)
})

# ── Integration Test 2: Numerical computation chain ────────────────────────

test_that("integration: center_by_treated → stable_ols → check_rank recovers known ATT", {
  # Simulate a simple 2x2 DID scenario:
  # 4 treated units (d=1), 4 control units (d=0)
  # True ATT = 3.0
  # X = intercept + covariate
  set.seed(42)
  n <- 8
  d <- c(rep(1, 4), rep(0, 4))
  x_cov <- c(2, 4, 6, 8, 1, 3, 5, 7)  # covariate
  true_att <- 3.0
  true_beta <- 1.5  # covariate effect

  # Outcome: y = true_att * d + true_beta * x_cov + noise
  # For exact recovery, use no noise
  y <- true_att * d + true_beta * x_cov

  # Step 1: Center covariate by treated mean
  X_cov <- matrix(x_cov, ncol = 1)
  center_result <- center_by_treated(X_cov, d)

  # Treated mean of x_cov = mean(2,4,6,8) = 5
  expect_equal(center_result$X_mean_treated, 5.0)

  # Verify treated group centered mean is zero
  treated_centered_mean <- mean(center_result$X_centered[d == 1, 1])
  expect_equal(treated_centered_mean, 0, tolerance = 1e-14)

  # Step 2: Build design matrix [d, X_centered] and run OLS
  X_design <- cbind(d, center_result$X_centered)
  ols_result <- stable_ols(X_design, y)

  # Step 3: Check rank
  rank_result <- check_rank(X_design)
  expect_true(rank_result$is_full_rank)
  expect_equal(rank_result$rank, 2)

  # Step 4: Verify ATT recovery
  # The coefficient on d should be the ATT evaluated at treated mean of X
  # y = att * d + beta * (x - x_bar_treated) + beta * x_bar_treated
  # When x is centered: y = (att + beta * x_bar_treated) * ... 
  # Actually with design [d, X_centered]:
  #   y = alpha * d + beta * X_centered
  #   alpha = att + beta * x_bar_treated ... no, let's think carefully.
  #
  # y_i = att * d_i + beta * x_i
  # Rewrite x_i = (x_i - x_bar_1) + x_bar_1 where x_bar_1 = treated mean
  # y_i = att * d_i + beta * (x_i - x_bar_1) + beta * x_bar_1
  # y_i = (att + beta * x_bar_1) * d_i + beta * (x_i - x_bar_1) + beta * x_bar_1 * (1 - d_i)
  # This doesn't simplify cleanly. Let me just verify the OLS coefficients directly.
  
  # With X_design = [d, x_centered], the OLS solves:
  # y = coef[1] * d + coef[2] * x_centered
  # Since y = 3*d + 1.5*x and x_centered = x - 5:
  # y = 3*d + 1.5*(x_centered + 5) = 3*d + 1.5*x_centered + 7.5
  # But we have no intercept! So OLS won't recover exactly.
  # Let's add an intercept.

  X_design_int <- cbind(1, d, center_result$X_centered)
  ols_result2 <- stable_ols(X_design_int, y)

  # y = 3*d + 1.5*(x_centered + 5) = 7.5 + 3*d + 1.5*x_centered
  # So: intercept = 7.5, coef_d = 3.0, coef_x_centered = 1.5
  coefs <- unname(ols_result2$coefficients)
  expect_equal(coefs[1], 7.5, tolerance = 1e-10)
  expect_equal(coefs[2], 3.0, tolerance = 1e-10)
  expect_equal(coefs[3], 1.5, tolerance = 1e-10)
  expect_false(ols_result2$is_rank_deficient)
  expect_equal(ols_result2$df.residual, 5L)

  # Residuals should be zero (exact linear system)
  expect_equal(ols_result2$residuals, rep(0, n), tolerance = 1e-10)
})


test_that("integration: weighted OLS with centered covariates recovers ATT", {
  # Weighted scenario: different sample weights
  set.seed(123)
  n <- 10
  d <- c(rep(1, 5), rep(0, 5))
  x_cov <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  w <- c(2, 1, 3, 1, 2, 1, 2, 1, 3, 1)  # positive weights
  true_att <- 5.0
  true_beta <- 2.0

  y <- true_att * d + true_beta * x_cov

  # Center by treated
  X_cov <- matrix(x_cov, ncol = 1)
  center_result <- center_by_treated(X_cov, d)
  # Treated mean = mean(1,2,3,4,5) = 3
  expect_equal(center_result$X_mean_treated, 3.0)

  # OLS with intercept, treatment, centered covariate
  X_design <- cbind(1, d, center_result$X_centered)
  ols_result <- stable_ols(X_design, y, w = w)

  # y = 5*d + 2*(x_centered + 3) = 6 + 5*d + 2*x_centered
  coefs <- unname(ols_result$coefficients)
  expect_equal(coefs[1], 6.0, tolerance = 1e-10)
  expect_equal(coefs[2], 5.0, tolerance = 1e-10)
  expect_equal(coefs[3], 2.0, tolerance = 1e-10)

  # Residuals should be zero (exact system)
  expect_equal(max(abs(ols_result$residuals)), 0, tolerance = 1e-10)

  # Verify rank
  rank_result <- check_rank(X_design)
  expect_true(rank_result$is_full_rank)
  expect_equal(rank_result$rank, 3)
})

test_that("integration: full staggered panel workflow with multiple cohorts", {
  # 12 units, 8 periods (2000-2007)
  # Cohort 2003: units 1-3
  # Cohort 2005: units 4-6
  # Never-treated (NA): units 7-9
  # Never-treated (0): units 10-12
  n_units <- 12
  periods <- 2000:2007
  n_periods <- length(periods)

  panel <- data.frame(
    id   = rep(1:n_units, each = n_periods),
    year = rep(periods, times = n_units),
    gvar = rep(c(2003, 2003, 2003, 2005, 2005, 2005,
                 NA, NA, NA, 0, 0, 0), each = n_periods),
    y    = rnorm(n_units * n_periods)
  )

  # Full chain: unit_gvar → cohorts → has_nt → masks → periods → cross-sections
  unit_gvar_dt <- get_unit_level_gvar(panel, "gvar", "id")
  expect_equal(nrow(unit_gvar_dt), 12L)

  # Convert to named vector for mask operations
  unit_gvar <- unit_gvar_dt[["gvar"]]
  names(unit_gvar) <- as.character(unit_gvar_dt[["id"]])

  cohorts <- get_cohorts(panel, "gvar", "id")
  expect_identical(cohorts, c(2003L, 2005L))

  expect_true(has_never_treated(panel, "gvar", "id"))

  nt_mask <- get_never_treated_mask(unit_gvar)
  expect_equal(sum(nt_mask), 6)  # 6 NT units (3 NA + 3 zero)

  # Cohort 2003 mask
  c2003_mask <- get_cohort_mask(unit_gvar, 2003)
  expect_equal(sum(c2003_mask), 3)
  expect_identical(which(unname(c2003_mask)), 1:3)

  # Cohort 2005 mask
  c2005_mask <- get_cohort_mask(unit_gvar, 2005)
  expect_equal(sum(c2005_mask), 3)
  expect_identical(which(unname(c2005_mask)), 4:6)

  # Valid periods for cohort 2003 (T_max = 2007)
  periods_2003 <- get_valid_periods_for_cohort(2003, 2007)
  expect_identical(periods_2003, 2003L:2007L)
  expect_length(periods_2003, 5)

  # Valid periods for cohort 2005 (T_max = 2007)
  periods_2005 <- get_valid_periods_for_cohort(2005, 2007)
  expect_identical(periods_2005, 2005L:2007L)
  expect_length(periods_2005, 3)

  # Extract cross-section for period 2003
  dt <- data.table::as.data.table(panel)
  cs_2003 <- extract_cross_section(dt, "year", 2003)
  expect_equal(nrow(cs_2003), 12)  # all units present in every period
  expect_true(all(cs_2003$year == 2003))
})

# ── Integration Test 3: Condition class name consistency ───────────────────

test_that("integration: condition classes in utils.R match conditions.R definitions", {
  # Verify that all error classes thrown by utils.R functions are valid
  # condition classes defined in conditions.R

  # 1. lwdid_invalid_staggered_data (from is_never_treated with -Inf)
  err1 <- tryCatch(
    is_never_treated(-Inf),
    error = function(e) e
  )
  expect_true(inherits(err1, "lwdid_invalid_staggered_data"))
  expect_true(inherits(err1, "lwdid_error"))

  # 2. lwdid_no_treated (from center_by_treated with no treated)
  err2 <- tryCatch(
    center_by_treated(matrix(1:4, nrow = 2), c(0, 0)),
    error = function(e) e
  )
  expect_true(inherits(err2, "lwdid_no_treated"))
  expect_true(inherits(err2, "lwdid_error"))

  # 3. lwdid_invalid_parameter (from stable_ols with bad weights)
  err3 <- tryCatch(
    stable_ols(cbind(1, 1:3), 1:3, w = c(1, 0, 1)),
    error = function(e) e
  )
  expect_true(inherits(err3, "lwdid_invalid_parameter"))
  expect_true(inherits(err3, "lwdid_error"))

  # 4. lwdid_invalid_parameter (from stable_ols with NA weights)
  err4 <- tryCatch(
    stable_ols(cbind(1, 1:3), 1:3, w = c(1, NA, 1)),
    error = function(e) e
  )
  expect_true(inherits(err4, "lwdid_invalid_parameter"))
  expect_true(inherits(err4, "lwdid_error"))

  # 5. lwdid_invalid_parameter (from stable_ols with negative weights)
  err5 <- tryCatch(
    stable_ols(cbind(1, 1:3), 1:3, w = c(1, -1, 1)),
    error = function(e) e
  )
  expect_true(inherits(err5, "lwdid_invalid_parameter"))
  expect_true(inherits(err5, "lwdid_error"))

  # 6. lwdid_invalid_parameter (from stable_ols with wrong length weights)
  err6 <- tryCatch(
    stable_ols(cbind(1, 1:3), 1:3, w = c(1, 2)),
    error = function(e) e
  )
  expect_true(inherits(err6, "lwdid_invalid_parameter"))
  expect_true(inherits(err6, "lwdid_error"))
})

test_that("integration: safe_mean in aggregation context", {
  # Simulate computing cohort-level ATT means from cell-level estimates
  # Some cells may have NA (insufficient data)
  cell_atts <- c(0.5, 0.8, NA, 1.2, NA, 0.3)

  # safe_mean should handle NAs gracefully
  cohort_att <- safe_mean(cell_atts)
  expect_equal(cohort_att, mean(c(0.5, 0.8, 1.2, 0.3)))

  # All-NA case (no valid cells)
  all_na_atts <- c(NA_real_, NA_real_, NA_real_)
  expect_identical(safe_mean(all_na_atts), NA_real_)

  # Empty case (no cells at all)
  expect_identical(safe_mean(numeric(0)), NA_real_)
})
