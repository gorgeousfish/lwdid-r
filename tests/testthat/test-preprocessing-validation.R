# test-preprocessing-validation.R — Tests for Story 11.1 validation functions

# ============================================================
# .validate_aggregation_inputs() tests
# ============================================================

test_that("non-data.frame input throws lwdid_invalid_parameter", {
  expect_error(
    .validate_aggregation_inputs(
      data = list(a = 1), unit_var = "state",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    .validate_aggregation_inputs(
      data = matrix(1:4, 2, 2), unit_var = "state",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("empty data throws lwdid_invalid_parameter", {
  empty_df <- data.frame(state = character(0),
                         year = integer(0),
                         income = numeric(0))
  expect_error(
    .validate_aggregation_inputs(
      data = empty_df, unit_var = "state",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("missing columns throw lwdid_missing_column", {
  df <- data.frame(state = "CA", year = 2000L, income = 50000)
  expect_error(
    .validate_aggregation_inputs(
      data = df, unit_var = "region",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_missing_column"
  )
  expect_error(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income",
      weight_var = "w"
    ),
    class = "lwdid_missing_column"
  )
  expect_error(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income",
      controls = c("age", "education")
    ),
    class = "lwdid_missing_column"
  )
})

test_that("non-numeric outcome throws lwdid_invalid_parameter", {
  df <- data.frame(state = "CA", year = 2000L, income = "high")
  expect_error(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_invalid_parameter"
  )
  df2 <- data.frame(state = "CA", year = 2000L,
                    income = factor("high"))
  expect_error(
    .validate_aggregation_inputs(
      data = df2, unit_var = "state",
      time_var = "year", outcome_var = "income"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("integer outcome is accepted", {
  df <- data.frame(state = "CA", year = 2000L, income = 50000L)
  expect_silent(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income"
    )
  )
})

test_that("logical outcome throws lwdid_invalid_parameter", {
  df <- data.frame(state = "CA", year = 2000L, treated = TRUE)
  expect_error(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "treated"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("valid input passes silently", {
  df <- data.frame(state = "CA", year = 2000L, income = 50000)
  expect_silent(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income"
    )
  )
})

test_that("data.table input is accepted", {
  dt <- data.table::data.table(state = "CA", year = 2000L,
                                income = 50000)
  expect_silent(
    .validate_aggregation_inputs(
      data = dt, unit_var = "state",
      time_var = "year", outcome_var = "income"
    )
  )
})

test_that("multi time_var handled correctly", {
  df <- data.frame(state = "CA", year = 2000L,
                   quarter = 1L, income = 50000)
  expect_silent(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = c("year", "quarter"),
      outcome_var = "income"
    )
  )
  df2 <- data.frame(state = "CA", year = 2000L, income = 50000)
  expect_error(
    .validate_aggregation_inputs(
      data = df2, unit_var = "state",
      time_var = c("year", "quarter"),
      outcome_var = "income"
    ),
    class = "lwdid_missing_column"
  )
})

test_that("multiple controls columns all checked", {
  df <- data.frame(state = "CA", year = 2000L,
                   income = 50000, age = 35, education = 16)
  expect_silent(
    .validate_aggregation_inputs(
      data = df, unit_var = "state",
      time_var = "year", outcome_var = "income",
      controls = c("age", "education")
    )
  )
})

# ============================================================
# .validate_weights() tests
# ============================================================

test_that("negative weights throw lwdid_invalid_parameter", {
  dt <- data.table::data.table(w = c(1.0, -0.5, 2.0))
  expect_error(
    .validate_weights(dt, "w"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("missing weights emit lwdid_data warning and remove rows", {
  dt <- data.table::data.table(
    x = 1:5, w = c(1.0, NA, 2.0, NA, 3.0)
  )
  result <- NULL
  expect_warning(
    result <- .validate_weights(dt, "w"),
    class = "lwdid_data"
  )
  expect_equal(nrow(result$data), 3L)
  expect_equal(result$n_missing, 2L)
})

test_that("all-valid weights produce no warning", {
  dt <- data.table::data.table(
    x = 1:3, w = c(1.0, 2.0, 3.0)
  )
  result <- expect_silent(.validate_weights(dt, "w"))
  expect_equal(nrow(result$data), 3L)
  expect_equal(result$n_missing, 0L)
})

test_that("zero weights are valid non-negative values", {
  dt <- data.table::data.table(
    x = 1:3, w = c(0.0, 0.0, 1.0)
  )
  result <- expect_silent(.validate_weights(dt, "w"))
  expect_equal(nrow(result$data), 3L)
  expect_equal(result$n_missing, 0L)
})

test_that("all-NA weights emit warning and return 0 rows", {
  dt <- data.table::data.table(
    x = 1:3, w = c(NA_real_, NA_real_, NA_real_)
  )
  result <- NULL
  expect_warning(
    result <- .validate_weights(dt, "w"),
    class = "lwdid_data"
  )
  expect_equal(nrow(result$data), 0L)
  expect_equal(result$n_missing, 3L)
})

test_that("mixed negative and NA: negative check takes priority", {
  dt <- data.table::data.table(
    x = 1:4, w = c(1.0, -0.5, NA_real_, 2.0)
  )
  expect_error(
    .validate_weights(dt, "w"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("no-missing weights return deep copy", {
  dt <- data.table::data.table(
    x = 1:3, w = c(1.0, 2.0, 3.0)
  )
  result <- expect_silent(.validate_weights(dt, "w"))
  result$data[, new_col := 999]
  expect_false("new_col" %in% names(dt))
})

test_that("Inf weight is valid non-negative value", {
  dt <- data.table::data.table(
    x = 1:3, w = c(1.0, Inf, 3.0)
  )
  result <- expect_silent(.validate_weights(dt, "w"))
  expect_equal(nrow(result$data), 3L)
  expect_equal(result$n_missing, 0L)
})

test_that("-Inf weight is caught as negative", {
  dt <- data.table::data.table(
    x = 1:3, w = c(1.0, -Inf, 3.0)
  )
  expect_error(
    .validate_weights(dt, "w"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("NaN weight treated as missing", {
  dt <- data.table::data.table(
    x = 1:3, w = c(1.0, NaN, 3.0)
  )
  result <- NULL
  expect_warning(
    result <- .validate_weights(dt, "w"),
    class = "lwdid_data"
  )
  expect_equal(nrow(result$data), 2L)
  expect_equal(result$n_missing, 1L)
})

# ============================================================
# .validate_treatment_consistency() tests
# ============================================================

test_that("consistent treatment within cell passes", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    treat = c(1, 1, 0, 0)
  )
  expect_silent(
    .validate_treatment_consistency(dt, "state", "year", "treat")
  )
})

test_that("inconsistent treatment throws lwdid_invalid_aggregation", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    treat = c(1, 0, 0, 0)
  )
  expect_error(
    .validate_treatment_consistency(dt, "state", "year", "treat"),
    class = "lwdid_invalid_aggregation"
  )
  expect_error(
    .validate_treatment_consistency(dt, "state", "year", "treat"),
    class = "lwdid_aggregation_error"
  )
})

test_that("single obs cell treated as consistent", {
  dt <- data.table::data.table(
    state = c("CA", "TX"),
    year = c(2000L, 2000L),
    treat = c(1, 0)
  )
  expect_silent(
    .validate_treatment_consistency(dt, "state", "year", "treat")
  )
})

test_that("multi time_var grouping detects inconsistency", {
  dt_ok <- data.table::data.table(
    state = c("CA", "CA", "CA", "CA"),
    year = c(2000L, 2000L, 2000L, 2000L),
    quarter = c(1L, 1L, 2L, 2L),
    treat = c(1, 1, 0, 0)
  )
  expect_silent(
    .validate_treatment_consistency(
      dt_ok, "state", c("year", "quarter"), "treat"
    )
  )
  dt_bad <- data.table::data.table(
    state = c("CA", "CA", "CA", "CA"),
    year = c(2000L, 2000L, 2000L, 2000L),
    quarter = c(1L, 1L, 2L, 2L),
    treat = c(1, 0, 0, 0)
  )
  expect_error(
    .validate_treatment_consistency(
      dt_bad, "state", c("year", "quarter"), "treat"
    ),
    class = "lwdid_invalid_aggregation"
  )
})

test_that("treatment with NA not misreported as inconsistent", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L, 2000L),
    treat = c(1, NA, 1, 0, 0)
  )
  expect_silent(
    .validate_treatment_consistency(dt, "state", "year", "treat")
  )
})

test_that("all-NA treatment within cell treated as consistent", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    treat = c(NA, NA, 0, 0)
  )
  expect_silent(
    .validate_treatment_consistency(dt, "state", "year", "treat")
  )
})

test_that("treatment with NA but actually inconsistent still errors", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA"),
    year = c(2000L, 2000L, 2000L),
    treat = c(1, NA, 0)
  )
  expect_error(
    .validate_treatment_consistency(dt, "state", "year", "treat"),
    class = "lwdid_invalid_aggregation"
  )
})

# ============================================================
# .validate_gvar_consistency() tests
# ============================================================

test_that("consistent gvar within unit passes", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2001L, 2000L, 2001L),
    gvar = c(2001L, 2001L, 0L, 0L)
  )
  expect_silent(
    .validate_gvar_consistency(dt, "state", "gvar")
  )
})

test_that("inconsistent gvar throws lwdid_invalid_aggregation", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2001L, 2000L, 2001L),
    gvar = c(2001L, 2002L, 0L, 0L)
  )
  expect_error(
    .validate_gvar_consistency(dt, "state", "gvar"),
    class = "lwdid_invalid_aggregation"
  )
  expect_error(
    .validate_gvar_consistency(dt, "state", "gvar"),
    class = "lwdid_aggregation_error"
  )
})

test_that("single obs unit treated as consistent", {
  dt <- data.table::data.table(
    state = c("CA", "TX"),
    year = c(2000L, 2000L),
    gvar = c(2001L, 0L)
  )
  expect_silent(
    .validate_gvar_consistency(dt, "state", "gvar")
  )
})

test_that("gvar with NA not misreported", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA", "TX", "TX"),
    year = c(2000L, 2001L, 2002L, 2000L, 2001L),
    gvar = c(2001L, NA, 2001L, 0L, 0L)
  )
  expect_silent(
    .validate_gvar_consistency(dt, "state", "gvar")
  )
})

test_that("all-NA gvar within unit treated as consistent", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2001L, 2000L, 2001L),
    gvar = c(NA_integer_, NA_integer_, 0L, 0L)
  )
  expect_silent(
    .validate_gvar_consistency(dt, "state", "gvar")
  )
})

test_that("gvar with NA but actually inconsistent still errors", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "CA"),
    year = c(2000L, 2001L, 2002L),
    gvar = c(2001L, NA, 2002L)
  )
  expect_error(
    .validate_gvar_consistency(dt, "state", "gvar"),
    class = "lwdid_invalid_aggregation"
  )
})

# ============================================================
# Error message context tests
# ============================================================

test_that("error messages contain context info", {
  df <- data.frame(state = "CA", year = 2000L, income = 50000)
  err <- tryCatch(
    .validate_aggregation_inputs(
      data = df, unit_var = "region",
      time_var = "year", outcome_var = "income"
    ),
    lwdid_missing_column = function(e) e
  )
  expect_true(grepl("region", err$message))
  expect_true(grepl("state", err$message))
  expect_equal(err$column, "region")
  expect_true("state" %in% err$available)
})

test_that("negative weight error includes count and examples", {
  dt <- data.table::data.table(w = c(1.0, -0.5, -2.0, 3.0))
  err <- tryCatch(
    .validate_weights(dt, "w"),
    lwdid_invalid_parameter = function(e) e
  )
  expect_true(grepl("2", err$message))
  expect_true(grepl("-0.5", err$message))
})

test_that("treatment inconsistency error includes cell count and examples", {
  dt <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX"),
    year = c(2000L, 2000L, 2000L, 2000L),
    treat = c(1, 0, 1, 0)
  )
  err <- tryCatch(
    .validate_treatment_consistency(dt, "state", "year", "treat"),
    lwdid_invalid_aggregation = function(e) e
  )
  expect_true(grepl("2 cell", err$message))
  expect_true(grepl("CA", err$message))
})
