# test-preprocessing-cellstats.R — Tests for Story 11.5 CellStatistics

# ============================================================
# new_cell_statistics() tests
# ============================================================

test_that("new_cell_statistics returns named list with all fields", {
  cs <- new_cell_statistics(
    unit = "CA",
    period = 2000L,
    n_obs = 100L,
    outcome_mean = 50000.0,
    outcome_variance = 1000.0,
    effective_sample_size = 95.0,
    weight_type = "survey"
  )
  expect_type(cs, "list")
  expect_equal(cs$unit, "CA")
  expect_equal(cs$period, 2000L)
  expect_equal(cs$n_obs, 100L)
  expect_equal(cs$outcome_mean, 50000.0)
  expect_equal(cs$outcome_variance, 1000.0)
  expect_equal(cs$effective_sample_size, 95.0)
  expect_equal(cs$weight_type, "survey")
})

test_that("NULL variance converted to NA_real_", {
  cs <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 3L,
    outcome_mean = 50000.0
  )
  expect_true(is.na(cs$outcome_variance))
  expect_type(cs$outcome_variance, "double")
})

test_that("NULL ESS converted to NA_real_", {
  cs <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 3L,
    outcome_mean = 50000.0
  )
  expect_true(is.na(cs$effective_sample_size))
  expect_type(cs$effective_sample_size, "double")
})

test_that("default weight_type is 'equal'", {
  cs <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 3L,
    outcome_mean = 50000.0
  )
  expect_equal(cs$weight_type, "equal")
})

test_that("n_obs is integer type", {
  cs <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 5,
    outcome_mean = 50000.0
  )
  expect_type(cs$n_obs, "integer")
})

test_that("outcome_mean is double type", {
  cs <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 3L,
    outcome_mean = 50000L
  )
  expect_type(cs$outcome_mean, "double")
})

test_that("list period collapsed to string", {
  cs <- new_cell_statistics(
    unit = "CA",
    period = list(year = 2000L, quarter = 1L),
    n_obs = 3L,
    outcome_mean = 50000.0
  )
  expect_type(cs$period, "character")
})

test_that("rbindlist works with multiple CellStatistics", {
  cs1 <- new_cell_statistics(
    unit = "CA", period = 2000L, n_obs = 3L,
    outcome_mean = 50000.0, outcome_variance = 100.0,
    effective_sample_size = 2.5, weight_type = "survey"
  )
  cs2 <- new_cell_statistics(
    unit = "TX", period = 2000L, n_obs = 2L,
    outcome_mean = 45000.0
  )
  dt <- data.table::rbindlist(list(cs1, cs2))
  expect_equal(nrow(dt), 2L)
  expect_equal(dt$unit, c("CA", "TX"))
  expect_true(is.na(dt$outcome_variance[2]))
  expect_true(is.na(dt$effective_sample_size[2]))
})
