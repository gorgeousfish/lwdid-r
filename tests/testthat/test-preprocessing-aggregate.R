# test-preprocessing-aggregate.R — Tests for Story 11.3 aggregate_to_panel()

# ============================================================
# Test data constructor
# ============================================================

make_test_data <- function() {
  set.seed(42)
  data.table::data.table(
    state = rep(c("CA", "TX", "NY"), each = 6),
    year = rep(rep(c(2000L, 2001L), each = 3), 3),
    income = c(50000, 55000, 60000,
               62000, 65000, 68000,
               45000, 48000, 52000,
               53000, 56000, 59000,
               70000, 72000, 75000,
               78000, 80000, 82000),
    weight = c(1.0, 1.2, 0.8,
               1.1, 0.9, 1.0,
               1.0, 1.5, 0.5,
               0.8, 1.2, 1.0,
               1.0, 1.0, 1.0,
               1.1, 0.9, 1.0),
    treat = rep(c(1L, 0L, 0L), each = 6),
    gvar = rep(c(2001L, 0L, 0L), each = 6),
    x1 = rnorm(18, 10, 2)
  )
}

# ============================================================
# Basic aggregation tests
# ============================================================

test_that("basic equal-weight aggregation correct", {
  dt <- make_test_data()
  result <- aggregate_to_panel(dt, "state", "year", "income")
  expect_s3_class(result, "lwdid_aggregation_result")
  expect_equal(result$n_cells, 6L)
  expect_equal(result$n_units, 3L)
  expect_equal(result$n_periods, 2L)
  expect_equal(result$n_original_obs, 18L)

  panel <- result$panel_data
  ca_2000 <- panel[state == "CA" & year == 2000L]
  expect_equal(ca_2000$income, 55000, tolerance = 1e-10)
})

test_that("survey weight aggregation correct", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight"
  )
  panel <- result$panel_data

  ca_2000 <- panel[state == "CA" & year == 2000L]
  w <- c(1.0, 1.2, 0.8)
  nw <- w / sum(w)
  expected <- sum(nw * c(50000, 55000, 60000))
  expect_equal(ca_2000$income, expected, tolerance = 1e-8)
})

test_that("min_cell_size filtering correct", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    income = c(50000, 55000, 60000, 45000)
  )
  result <- NULL
  expect_warning(
    result <- aggregate_to_panel(dt, "state", "year", "income",
                                  min_cell_size = 2L),
    class = "lwdid_data"
  )
  expect_equal(result$n_cells, 1L)
  expect_equal(result$n_excluded_cells, 1L)
})

test_that("all-excluded throws lwdid_insufficient_cell_size", {
  dt <- data.table::data.table(
    state = c("CA", "TX"),
    year = c(2000L, 2000L),
    income = c(50000, 45000)
  )
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = 5L),
    class = "lwdid_insufficient_cell_size"
  )
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = 5L),
    class = "lwdid_aggregation_error"
  )
})

test_that("fewer than 3 units emits lwdid_small_sample warning", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2001L, 2000L, 2001L),
    income = c(50000, 55000, 45000, 48000)
  )
  expect_warning(
    aggregate_to_panel(dt, "state", "year", "income"),
    class = "lwdid_small_sample"
  )
})

test_that("compute_variance=TRUE outputs variance column", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    compute_variance = TRUE
  )
  expect_true("income_var" %in% names(result$panel_data))
})

test_that("survey weights output ESS column", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight"
  )
  expect_true("_ess" %in% names(result$panel_data))
})

test_that("_n_obs column always present", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income"
  )
  expect_true("_n_obs" %in% names(result$panel_data))
  expect_true(all(result$panel_data[["_n_obs"]] == 3L))
})

test_that("controls correctly aggregated", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    controls = "x1"
  )
  expect_true("x1" %in% names(result$panel_data))
})

test_that("treatment_var correctly passed through", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    treatment_var = "treat"
  )
  expect_true("treat" %in% names(result$panel_data))
  panel <- result$panel_data
  expect_equal(panel[state == "CA"]$treat[1], 1L)
  expect_equal(panel[state == "TX"]$treat[1], 0L)
})

test_that("gvar correctly passed through", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    gvar = "gvar"
  )
  expect_true("gvar" %in% names(result$panel_data))
})

# ============================================================
# Parameter validation tests
# ============================================================

test_that("invalid frequency throws lwdid_invalid_parameter", {
  dt <- make_test_data()
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       frequency = "daily"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("frequency is case insensitive", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    frequency = "Annual"
  )
  expect_equal(result$frequency, "annual")

  result2 <- aggregate_to_panel(
    dt, "state", "year", "income",
    frequency = "QUARTERLY"
  )
  expect_equal(result2$frequency, "quarterly")
})

test_that("invalid min_cell_size throws lwdid_invalid_parameter", {
  dt <- make_test_data()
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = 0L),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = -1L),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = 2.5),
    class = "lwdid_invalid_parameter"
  )
})

test_that("min_cell_size accepts double integer values", {
  dt <- make_test_data()
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    min_cell_size = 2
  )
  expect_s3_class(result, "lwdid_aggregation_result")
  result2 <- aggregate_to_panel(
    dt, "state", "year", "income",
    min_cell_size = 2.0
  )
  expect_s3_class(result2, "lwdid_aggregation_result")
})

test_that("min_cell_size string throws lwdid_invalid_parameter", {
  dt <- make_test_data()
  expect_error(
    aggregate_to_panel(dt, "state", "year", "income",
                       min_cell_size = "2"),
    class = "lwdid_invalid_parameter"
  )
})

# ============================================================
# Defensive copy test
# ============================================================

test_that("original data not modified", {
  dt <- make_test_data()
  dt_copy <- data.table::copy(dt)
  result <- aggregate_to_panel(dt, "state", "year", "income")
  expect_identical(dt, dt_copy)
})

# ============================================================
# Multi-frequency test
# ============================================================

test_that("quarterly data correctly aggregated", {
  dt <- data.table::data.table(
    state = rep("CA", 12),
    year = rep(2000L, 12),
    quarter = rep(1:4, each = 3),
    income = rnorm(12, 50000, 5000)
  )
  result <- aggregate_to_panel(
    dt, "state", c("year", "quarter"), "income",
    frequency = "quarterly"
  )
  expect_equal(result$n_cells, 4L)
  expect_equal(result$n_periods, 4L)
})

# ============================================================
# All-NA outcome test
# ============================================================

test_that("all-NA outcome cell excluded", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    income = c(NA, NA, 45000, 48000)
  )
  result <- NULL
  expect_warning(
    result <- aggregate_to_panel(dt, "state", "year", "income"),
    class = "lwdid_data"
  )
  expect_equal(result$n_cells, 1L)
  expect_equal(result$n_excluded_cells, 1L)
})

# ============================================================
# Large-scale numerical precision test
# ============================================================

test_that("large-scale aggregation numerical precision", {
  set.seed(123)
  n_per_cell <- 1000L
  states <- c("CA", "TX", "NY", "FL", "IL")
  years <- 2000:2004
  dt <- data.table::CJ(
    state = states, year = years, obs = 1:n_per_cell
  )
  dt[, income := rnorm(.N, 50000, 10000)]
  dt[, weight := runif(.N, 0.5, 2.0)]
  dt[, obs := NULL]

  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight"
  )

  expect_equal(result$n_cells, 25L)
  expect_equal(result$n_units, 5L)
  expect_equal(result$n_periods, 5L)
  expect_equal(result$n_original_obs,
               as.integer(5 * 5 * n_per_cell))

  expect_equal(result$min_cell_size_stat, n_per_cell)
  expect_equal(result$max_cell_size, n_per_cell)
})

# ============================================================
# NA group key retention test
# ============================================================

test_that("NA group key retained (not silently dropped)", {
  dt <- data.table::data.table(
    state = c("CA", "CA", NA, NA, "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L, 2000L, 2000L),
    income = c(50000, 55000, 60000, 65000, 45000, 48000)
  )
  result <- aggregate_to_panel(dt, "state", "year", "income")
  expect_equal(result$n_cells, 3L)
})

# ============================================================
# All-NA weights test
# ============================================================

test_that("all-NA weights correctly throws error", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    income = c(50000, 55000, 45000, 48000),
    w = c(NA_real_, NA_real_, NA_real_, NA_real_)
  )
  expect_error(
    suppressWarnings(
      aggregate_to_panel(dt, "state", "year", "income",
                         weight_var = "w")
    ),
    class = "lwdid_insufficient_cell_size"
  )
})

# ============================================================
# Single observation cell test
# ============================================================

test_that("single obs cell correctly aggregated (min_cell_size=1)", {
  dt <- data.table::data.table(
    state = c("CA", "TX", "NY"),
    year = c(2000L, 2000L, 2000L),
    income = c(50000, 45000, 70000)
  )
  result <- aggregate_to_panel(dt, "state", "year", "income")
  expect_equal(result$n_cells, 3L)
  panel <- result$panel_data
  expect_equal(panel[state == "CA"]$income, 50000)
  expect_true(all(panel[["_n_obs"]] == 1L))
})

# ============================================================
# Control variable numerical verification
# ============================================================

test_that("control weighted mean numerically correct", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA"),
    year = c(2000L, 2000L, 2000L),
    income = c(50000, 55000, 60000),
    weight = c(1.0, 2.0, 1.0),
    age = c(30, 40, 50)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight", controls = "age"
  )
  panel <- result$panel_data
  expect_equal(panel$age[1], 40, tolerance = 1e-10)
})

# ============================================================
# Variance numerical verification
# ============================================================

test_that("compute_variance=TRUE weighted variance correct", {
  dt <- data.table::data.table(
    state = rep("CA", 3),
    year = rep(2000L, 3),
    income = c(50000, 55000, 60000),
    weight = c(1.0, 1.2, 0.8)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight", compute_variance = TRUE
  )
  panel <- result$panel_data
  expect_equal(panel$income_var[1], 14888888.888888888,
               tolerance = 1e-4)
})

test_that("compute_variance=TRUE equal weight variance correct", {
  dt <- data.table::data.table(
    state = rep("CA", 3),
    year = rep(2000L, 3),
    income = c(50000, 55000, 60000)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    compute_variance = TRUE
  )
  panel <- result$panel_data
  expect_equal(panel$income_var[1], 16666666.666666666,
               tolerance = 1e-4)
})

# ============================================================
# ESS numerical verification
# ============================================================

test_that("ESS numerically correct", {
  dt <- data.table::data.table(
    state = rep("CA", 3),
    year = rep(2000L, 3),
    income = c(50000, 55000, 60000),
    weight = c(1.0, 1.2, 0.8)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight"
  )
  panel <- result$panel_data
  expect_equal(panel[["_ess"]][1], 2.922077922077922,
               tolerance = 1e-10)
})

# ============================================================
# Output column order verification
# ============================================================

test_that("output panel column order correct", {
  dt <- data.table::data.table(
    state = rep("CA", 3),
    year = rep(2000L, 3),
    income = c(50000, 55000, 60000),
    weight = c(1.0, 1.2, 0.8),
    treat = rep(1L, 3),
    gvar = rep(2001L, 3),
    x1 = c(10, 20, 30)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight", controls = "x1",
    treatment_var = "treat", gvar = "gvar",
    compute_variance = TRUE
  )
  col_names <- names(result$panel_data)
  expected_order <- c("state", "year", "income", "_n_obs",
                      "income_var", "_ess", "x1", "treat", "gvar")
  expect_equal(col_names, expected_order)
})

# ============================================================
# Multiple controls test
# ============================================================

test_that("multiple controls all correctly aggregated", {
  dt <- data.table::data.table(
    state = rep("CA", 3),
    year = rep(2000L, 3),
    income = c(50000, 55000, 60000),
    weight = c(1.0, 2.0, 1.0),
    age = c(30, 40, 50),
    edu = c(12, 16, 20)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight", controls = c("age", "edu")
  )
  panel <- result$panel_data
  expect_true("age" %in% names(panel))
  expect_true("edu" %in% names(panel))
  expect_equal(panel$age[1], 40, tolerance = 1e-10)
  expect_equal(panel$edu[1], 16, tolerance = 1e-10)
})

# ============================================================
# Comprehensive end-to-end numerical verification
# ============================================================

test_that("comprehensive e2e: all fields numerically verified", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L, 2000L),
    income = c(50000, 55000, 60000, 45000, 48000),
    weight = c(1.0, 1.2, 0.8, 1.0, 1.0),
    treat = c(1L, 1L, 1L, 0L, 0L),
    gvar = c(2001L, 2001L, 2001L, 0L, 0L),
    age = c(30, 40, 50, 25, 35)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    weight_var = "weight", controls = "age",
    treatment_var = "treat", gvar = "gvar",
    compute_variance = TRUE
  )

  # Metadata
  expect_equal(result$n_original_obs, 5L)
  expect_equal(result$n_cells, 2L)
  expect_equal(result$n_units, 2L)
  expect_equal(result$n_periods, 1L)
  expect_equal(result$n_excluded_cells, 0L)
  expect_equal(result$frequency, "annual")
  expect_equal(result$unit_var, "state")
  expect_equal(result$time_var, "year")
  expect_equal(result$outcome_var, "income")
  expect_equal(result$weight_var, "weight")

  panel <- result$panel_data

  # CA-2000
  ca <- panel[state == "CA"]
  expect_equal(ca$income[1], 54666.66666666666, tolerance = 1e-6)
  expect_equal(ca[["_n_obs"]][1], 3L)
  expect_equal(ca$income_var[1], 14888888.888888888, tolerance = 1e-2)
  expect_equal(ca[["_ess"]][1], 2.922077922077922, tolerance = 1e-8)
  expect_equal(ca$age[1], 10 + 16 + 50 * 4 / 15, tolerance = 1e-8)
  expect_equal(ca$treat[1], 1L)
  expect_equal(ca$gvar[1], 2001L)

  # TX-2000
  tx <- panel[state == "TX"]
  expect_equal(tx$income[1], 46500, tolerance = 1e-10)
  expect_equal(tx[["_n_obs"]][1], 2L)
  expect_equal(tx$treat[1], 0L)
  expect_equal(tx$gvar[1], 0L)

  # Cell stats
  cs <- result$cell_stats
  expect_equal(nrow(cs), 2L)

  # Cell size stats
  expect_equal(result$min_cell_size_stat, 2L)
  expect_equal(result$max_cell_size, 3L)
  expect_equal(result$mean_cell_size, 2.5)
  expect_equal(result$median_cell_size, 2.5)
})

# ============================================================
# Single observation variance boundary test
# ============================================================

test_that("single obs compute_variance=TRUE returns NA variance", {
  dt <- data.table::data.table(
    state = c("CA"),
    year = c(2000L),
    income = c(50000)
  )
  result <- aggregate_to_panel(
    dt, "state", "year", "income",
    compute_variance = TRUE
  )
  panel <- result$panel_data
  expect_true(is.na(panel$income_var[1]))
})
