# test-preprocessing-result.R — Tests for Story 11.4 AggregationResult S3 class

# ============================================================
# Helper: create test result
# ============================================================

make_test_result <- function() {
  panel <- data.table::data.table(
    state = c("CA", "CA", "TX", "TX", "NY", "NY"),
    year = c(2000L, 2001L, 2000L, 2001L, 2000L, 2001L),
    income = c(55000, 65000, 48333, 56000, 72333, 80000),
    `_n_obs` = rep(3L, 6)
  )
  cell_stats <- data.table::data.table(
    unit = rep(c("CA", "TX", "NY"), each = 2),
    period = rep(c(2000L, 2001L), 3),
    n_obs = rep(3L, 6),
    outcome_mean = panel$income,
    outcome_variance = rep(NA_real_, 6),
    effective_sample_size = rep(NA_real_, 6),
    weight_type = rep("equal", 6)
  )
  new_lwdid_aggregation_result(
    panel_data = panel,
    n_original_obs = 18L,
    n_cells = 6L,
    n_units = 3L,
    n_periods = 2L,
    cell_stats = cell_stats,
    min_cell_size_stat = 3L,
    max_cell_size = 3L,
    mean_cell_size = 3.0,
    median_cell_size = 3.0,
    unit_var = "state",
    time_var = "year",
    outcome_var = "income",
    weight_var = NULL,
    frequency = "annual",
    n_excluded_cells = 0L,
    excluded_cells_info = list()
  )
}

# ============================================================
# Constructor tests
# ============================================================

test_that("constructor returns correct S3 class", {
  result <- make_test_result()
  expect_s3_class(result, "lwdid_aggregation_result")
  expect_true(inherits(result, "lwdid_aggregation_result"))
})

test_that("all fields accessible via $", {
  result <- make_test_result()
  expect_equal(result$n_original_obs, 18L)
  expect_equal(result$n_cells, 6L)
  expect_equal(result$n_units, 3L)
  expect_equal(result$n_periods, 2L)
  expect_equal(result$unit_var, "state")
  expect_equal(result$time_var, "year")
  expect_equal(result$outcome_var, "income")
  expect_null(result$weight_var)
  expect_equal(result$frequency, "annual")
  expect_equal(result$n_excluded_cells, 0L)
})

test_that("integer fields are integer type", {
  result <- make_test_result()
  expect_type(result$n_original_obs, "integer")
  expect_type(result$n_cells, "integer")
  expect_type(result$n_units, "integer")
  expect_type(result$n_periods, "integer")
  expect_type(result$min_cell_size_stat, "integer")
  expect_type(result$max_cell_size, "integer")
})

test_that("double fields are double type", {
  result <- make_test_result()
  expect_type(result$mean_cell_size, "double")
  expect_type(result$median_cell_size, "double")
})

test_that("constructor rejects non-data.table panel_data", {
  expect_error(
    new_lwdid_aggregation_result(
      panel_data = data.frame(x = 1),
      n_original_obs = 1L, n_cells = 1L,
      n_units = 1L, n_periods = 1L,
      cell_stats = data.table::data.table(x = 1),
      min_cell_size_stat = 1L, max_cell_size = 1L,
      mean_cell_size = 1.0, median_cell_size = 1.0,
      unit_var = "x", time_var = "t",
      outcome_var = "y", weight_var = NULL,
      frequency = "annual"
    ),
    "panel_data must be a data.table"
  )
})

test_that("constructor rejects non-data.table cell_stats", {
  expect_error(
    new_lwdid_aggregation_result(
      panel_data = data.table::data.table(x = 1),
      n_original_obs = 1L, n_cells = 1L,
      n_units = 1L, n_periods = 1L,
      cell_stats = data.frame(x = 1),
      min_cell_size_stat = 1L, max_cell_size = 1L,
      mean_cell_size = 1.0, median_cell_size = 1.0,
      unit_var = "x", time_var = "t",
      outcome_var = "y", weight_var = NULL,
      frequency = "annual"
    ),
    "cell_stats must be a data.table"
  )
})

# ============================================================
# print method tests
# ============================================================

test_that("print outputs key info", {
  result <- make_test_result()
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Aggregation Summary")
  expect_match(output_str, "18")
  expect_match(output_str, "6")
  expect_match(output_str, "3")
  expect_match(output_str, "2")
  expect_match(output_str, "annual")
  expect_match(output_str, "equal weights")
})

test_that("print returns invisible(x)", {
  result <- make_test_result()
  ret <- withVisible(print(result))
  expect_false(ret$visible)
  expect_identical(ret$value, result)
})

test_that("print shows excluded cells count", {
  result <- make_test_result()
  result$n_excluded_cells <- 2L
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Excluded cells: 2")
})

test_that("print shows weight variable name", {
  result <- make_test_result()
  result$weight_var <- "survey_wt"
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "survey_wt")
})

test_that("print does not contain Cell Size Statistics", {
  result <- make_test_result()
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_false(grepl("Cell Size Statistics", output_str))
  expect_false(grepl("Configuration", output_str))
  expect_false(grepl("Minimum", output_str))
})

test_that("print hides Excluded cells when 0", {
  result <- make_test_result()
  expect_equal(result$n_excluded_cells, 0L)
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_false(grepl("Excluded cells", output_str))
})

test_that("print uses comma separator for large numbers", {
  result <- make_test_result()
  result$n_original_obs <- 10000L
  output <- capture.output(print(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "10,000")
})

# ============================================================
# summary method tests
# ============================================================

test_that("summary contains Cell Size Statistics", {
  result <- make_test_result()
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Cell Size Statistics")
  expect_match(output_str, "Minimum")
  expect_match(output_str, "Maximum")
  expect_match(output_str, "Mean")
  expect_match(output_str, "Median")
})

test_that("summary contains Configuration", {
  result <- make_test_result()
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Configuration")
  expect_match(output_str, "Unit variable: state")
  expect_match(output_str, "Time variable: year")
  expect_match(output_str, "Outcome variable: income")
})

test_that("summary shows multi time_var correctly", {
  result <- make_test_result()
  result$time_var <- c("year", "quarter")
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "year, quarter")
})

test_that("summary returns invisible(object)", {
  result <- make_test_result()
  ret <- withVisible(summary(result))
  expect_false(ret$visible)
  expect_identical(ret$value, result)
})

test_that("summary hides Excluded cells when 0", {
  result <- make_test_result()
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_false(grepl("Excluded cells", output_str))
})

test_that("summary shows excluded cells count", {
  result <- make_test_result()
  result$n_excluded_cells <- 5L
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Excluded cells: 5")
})

test_that("summary shows weight variable name", {
  result <- make_test_result()
  result$weight_var <- "survey_wt"
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "Weight variable: survey_wt")
})

test_that("summary shows None (equal weights) when no weight", {
  result <- make_test_result()
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "None \\(equal weights\\)")
})

test_that("summary uses comma separator for large numbers", {
  result <- make_test_result()
  result$n_original_obs <- 10000L
  result$min_cell_size_stat <- 1500L
  result$max_cell_size <- 25000L
  output <- capture.output(summary(result))
  output_str <- paste(output, collapse = "\n")
  expect_match(output_str, "10,000")
  expect_match(output_str, "Minimum: 1,500")
  expect_match(output_str, "Maximum: 25,000")
})

# ============================================================
# as.data.frame method tests
# ============================================================

test_that("as.data.frame returns data.frame", {
  result <- make_test_result()
  df <- as.data.frame(result)
  expect_s3_class(df, "data.frame")
  expect_false(data.table::is.data.table(df))
  expect_equal(nrow(df), 6L)
  expect_true("state" %in% names(df))
  expect_true("income" %in% names(df))
})

test_that("as.data.frame preserves data values", {
  result <- make_test_result()
  df <- as.data.frame(result)
  expect_equal(ncol(df), ncol(result$panel_data))
  expect_equal(df$state, c("CA", "CA", "TX", "TX", "NY", "NY"))
  expect_equal(df$year, c(2000L, 2001L, 2000L, 2001L, 2000L, 2001L))
})

test_that("as.data.frame supports row.names", {
  result <- make_test_result()
  rn <- paste0("row", seq_len(6))
  df <- as.data.frame(result, row.names = rn)
  expect_equal(rownames(df), rn)
})

# ============================================================
# to_csv method tests
# ============================================================

test_that("to_csv writes CSV with metadata", {
  result <- make_test_result()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = TRUE)

  lines <- readLines(tmp)
  expect_match(lines[1], "^# Aggregation Metadata")
  expect_match(lines[2], "^# Original observations")

  df <- read.csv(tmp, comment.char = "#")
  expect_equal(nrow(df), 6L)
})

test_that("to_csv writes CSV without metadata", {
  result <- make_test_result()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = FALSE)

  lines <- readLines(tmp)
  expect_false(grepl("^#", lines[1]))
})

test_that("to_csv returns invisible(x)", {
  result <- make_test_result()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  ret <- withVisible(to_csv(result, tmp))
  expect_false(ret$visible)
  expect_identical(ret$value, result)
})

test_that("to_csv shows multi time_var in metadata", {
  result <- make_test_result()
  result$time_var <- c("year", "quarter")
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = TRUE)

  lines <- readLines(tmp)
  time_line <- grep("^# Time variable:", lines, value = TRUE)
  expect_match(time_line, "year, quarter")
})

test_that("to_csv shows weight_var name in metadata", {
  result <- make_test_result()
  result$weight_var <- "survey_wt"
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = TRUE)

  lines <- readLines(tmp)
  wt_line <- grep("^# Weight variable:", lines, value = TRUE)
  expect_match(wt_line, "survey_wt")
})

test_that("to_csv shows None when weight_var is NULL", {
  result <- make_test_result()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = TRUE)

  lines <- readLines(tmp)
  wt_line <- grep("^# Weight variable:", lines, value = TRUE)
  expect_match(wt_line, "None")
})

test_that("to_csv metadata has correct number of lines", {
  result <- make_test_result()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  to_csv(result, tmp, include_metadata = TRUE)

  lines <- readLines(tmp)
  comment_lines <- grep("^#", lines)
  expect_equal(length(comment_lines), 11L)
})

# ============================================================
# to_dict method tests
# ============================================================

test_that("to_dict returns correct named list", {
  result <- make_test_result()
  d <- to_dict(result)
  expect_type(d, "list")
  expect_equal(d$n_original_obs, 18L)
  expect_equal(d$n_cells, 6L)
  expect_equal(d$n_units, 3L)
  expect_equal(d$n_periods, 2L)
  expect_equal(d$min_cell_size, 3L)
  expect_equal(d$max_cell_size, 3L)
  expect_equal(d$mean_cell_size, 3.0)
  expect_equal(d$median_cell_size, 3.0)
  expect_equal(d$unit_var, "state")
  expect_equal(d$time_var, "year")
  expect_equal(d$outcome_var, "income")
  expect_null(d$weight_var)
  expect_equal(d$frequency, "annual")
  expect_equal(d$n_excluded_cells, 0L)
})

test_that("to_dict excludes panel_data, cell_stats, excluded_cells_info", {
  result <- make_test_result()
  d <- to_dict(result)
  expect_null(d$panel_data)
  expect_null(d$cell_stats)
  expect_null(d$excluded_cells_info)
})

test_that("to_dict with weight_var", {
  result <- make_test_result()
  result$weight_var <- "survey_wt"
  d <- to_dict(result)
  expect_equal(d$weight_var, "survey_wt")
})

test_that("to_dict with vector time_var", {
  result <- make_test_result()
  result$time_var <- c("year", "quarter")
  d <- to_dict(result)
  expect_equal(d$time_var, c("year", "quarter"))
})

test_that("to_dict key names match Python", {
  result <- make_test_result()
  d <- to_dict(result)
  expected_keys <- c(
    "n_original_obs", "n_cells", "n_units", "n_periods",
    "min_cell_size", "max_cell_size", "mean_cell_size",
    "median_cell_size", "unit_var", "time_var",
    "outcome_var", "weight_var", "frequency",
    "n_excluded_cells"
  )
  expect_equal(sort(names(d)), sort(expected_keys))
  expect_equal(length(d), 14L)
})
