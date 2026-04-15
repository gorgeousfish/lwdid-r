make_e901_sparse_season_panel <- function() {
  data.frame(
    id = rep(1L, 5L),
    tindex = 1:5,
    y = c(10, 11, 12, 13, 14),
    quarter = c(1, 1, 2, 3, 4),
    post = c(0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_validate_inputs_panel <- function() {
  template <- make_e901_sparse_season_panel()
  names(template)[names(template) == "tindex"] <- "time"
  panel <- do.call(
    rbind,
    lapply(
      seq_len(3),
      function(uid) {
        unit_panel <- template
        unit_panel$id <- uid
        unit_panel$d <- if (uid < 3L) 1L else 0L
        unit_panel
      }
    )
  )
  panel
}

make_e901_transform_panel <- function() {
  seasonal_pattern <- c(10, 20, 30, 40, 10, 20, 30)
  unit_offsets <- c(0, 100, 200, 300)
  post_adjustments <- list(
    c(6, 4),
    c(4, 6),
    c(1, -1),
    c(-1, 1)
  )

  do.call(
    rbind,
    lapply(seq_along(unit_offsets), function(idx) {
      data.frame(
        id = idx,
        time = 1:7,
        quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L),
        d = as.integer(idx <= 2L),
        post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L),
        y = unit_offsets[[idx]] +
          seasonal_pattern +
          c(0, 0, 0, 0, 0, post_adjustments[[idx]])
      )
    })
  )
}

make_e901_detrendq_panel <- function() {
  seasonal_effects <- c(`1` = 0, `2` = 1, `3` = 0.5, `4` = 1.5)
  quarter_pattern <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L)
  post_pattern <- c(rep(0L, 6L), rep(1L, 4L))
  intercepts <- c(5, 15, 25, 35)
  slopes <- c(0.5, 0.7, 0.3, 0.4)

  do.call(
    rbind,
    lapply(seq_along(intercepts), function(idx) {
      effect <- if (idx <= 2L) 4 else 0
      data.frame(
        id = idx,
        time = 1:10,
        quarter = quarter_pattern,
        d = as.integer(idx <= 2L),
        post = post_pattern,
        y = intercepts[[idx]] +
          slopes[[idx]] * (1:10) +
          unname(seasonal_effects[as.character(quarter_pattern)]) +
          effect * post_pattern
      )
    })
  )
}

make_e901_detrendq_exclude_panel <- function() {
  seasonal_effects <- c(`1` = 0, `2` = 1, `3` = 0.5, `4` = 1.5)
  quarter_pattern <- c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L)
  post_pattern <- c(rep(0L, 7L), rep(1L, 4L))
  intercepts <- c(5, 15, 25, 35)
  slopes <- c(0.5, 0.7, 0.3, 0.4)

  do.call(
    rbind,
    lapply(seq_along(intercepts), function(idx) {
      effect <- if (idx <= 2L) 4 else 0
      data.frame(
        id = idx,
        time = 1:11,
        quarter = quarter_pattern,
        d = as.integer(idx <= 2L),
        post = post_pattern,
        y = intercepts[[idx]] +
          slopes[[idx]] * (1:11) +
          unname(seasonal_effects[as.character(quarter_pattern)]) +
          effect * post_pattern
      )
    })
  )
}

make_e901_transform_gate_fail_panel <- function() {
  data.frame(
    id = rep(1:3, each = 4),
    time = rep(1:4, times = 3),
    quarter = c(1L, 2L, 1L, 3L, 1L, 2L, 1L, 3L, 1L, 2L, 1L, 3L),
    d = c(rep(1L, 4), rep(1L, 4), rep(0L, 4)),
    post = rep(c(0L, 0L, 1L, 1L), times = 3),
    y = c(10, 20, 12, 33, 100, 200, 102, 203, 50, 60, 52, 63)
  )
}

make_e901_all_na_pre_tindex_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = c(NA, NA, NA, NA, NA, 6),
    y = c(10, 11, 12, 13, 14, 15),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_balanced_quarter_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = 1:6,
    y = c(10, 11, 12, 13, 14, 15),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_unit_without_pre_periods_panel <- function() {
  rbind(
    data.frame(
      id = 1L,
      tindex = 6L,
      y = 15,
      quarter = 2L,
      post = 1L
    ),
    data.frame(
      id = 2L,
      tindex = 1:6,
      y = 20:25,
      quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
      post = c(0L, 0L, 0L, 0L, 0L, 1L)
    )
  )
}

make_e901_short_pre_period_unit_panel <- function() {
  rbind(
    data.frame(
      id = 1L,
      tindex = c(5L, 6L),
      y = c(14, 15),
      quarter = c(1L, 2L),
      post = c(0L, 1L)
    ),
    data.frame(
      id = 2L,
      tindex = 1:6,
      y = 20:25,
      quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
      post = c(0L, 0L, 0L, 0L, 0L, 1L)
    )
  )
}

make_e901_na_pre_outcome_panel <- function() {
  data.frame(
    id = rep(1L, 5L),
    tindex = 1:5,
    y = c(10, 11, 12, NA, 13),
    quarter = c(1L, 2L, 3L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_na_season_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = 1:6,
    y = c(10, 11, 12, 13, 14, 15),
    quarter = c(1L, 2L, 3L, NA, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_detrendq_missing_tindex_panel <- function() {
  data.frame(
    id = rep(1L, 6L),
    tindex = c(1, 2, 3, NA, 5, 6),
    y = c(10, 12, 14, 16, 18, 20),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L)
  )
}

make_e901_detrendq_one_pre_period_panel <- function() {
  data.frame(
    id = rep(1L, 3L),
    tindex = c(1, 1, 2),
    y = c(10, 11, 12),
    quarter = c(1L, 2L, 1L),
    post = c(0L, 0L, 1L)
  )
}

make_e901_monthly_panel <- function() {
  data.frame(
    id = rep(1L, 14L),
    tindex = 1:14,
    y = 100 + seq_len(14),
    month = c(1:12, 1L, 2L),
    post = c(rep(0L, 13L), 1L)
  )
}

make_e901_weekly_panel <- function() {
  data.frame(
    id = rep(1L, 54L),
    tindex = 1:54,
    y = 200 + seq_len(54),
    week = c(1:52, 1L, 2L),
    post = c(rep(0L, 53L), 1L)
  )
}

read_e901_required_parity_json <- function(filename) {
  skip_if_not_installed("jsonlite")

  path <- resolve_parity_fixture_path(filename)
  expect_true(file.exists(path))
  if (!file.exists(path)) {
    return(NULL)
  }

  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

test_that("story-E9-01 seasonal validation suite lands all planned TC-9.1 labels", {
  lines <- readLines(testthat::test_path("test-seasonal-validation.R"), warn = FALSE)
  labels <- unlist(regmatches(
    lines,
    gregexpr("TC-9\\.1\\.[0-9]+", lines)
  ))
  labels <- sort(unique(as.integer(sub("^TC-9\\.1\\.", "", labels))))
  missing <- setdiff(1:36, labels)

  expect_identical(missing, integer(0))
})

test_that("TC-9.1.3/TC-9.1.9: unified seasonal gate warns on uncovered post season but accepts sparse pre seasons via n_unique_seasons", {
  data <- make_e901_sparse_season_panel()

  expect_warning(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(is.list(result))
  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 4L)
  expect_identical(result$n_units_valid, 1L)
  expect_identical(result$n_units_invalid, 0L)
  expect_identical(result$season_coverage$all_covered, FALSE)
  expect_identical(result$season_coverage$total_uncovered, 1L)
  expect_identical(
    result$season_coverage$units_with_gaps[["1"]],
    4
  )
})

test_that("TC-9.1.4: seasonal demeanq validation errors when a unit lacks df >= 1", {
  data <- data.frame(
    id = rep(1L, 3L),
    tindex = 1:3,
    y = c(10, 11, 12),
    quarter = c(1L, 2L, 3L),
    post = c(0L, 0L, 1L)
  )

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "requires at least 3 observation\\(s\\) so that df >= 1",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.1/TC-9.1.11: balanced quarterly data passes validation and returns the max pre-treatment K", {
  data <- make_e901_balanced_quarter_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 5L)
  expect_identical(result$n_units_valid, 1L)
  expect_identical(result$season_coverage$all_covered, TRUE)
  expect_identical(result$season_coverage$total_uncovered, 0L)
})

test_that("TC-9.1.2: seasonal values outside 1:Q raise a hard error", {
  data <- make_e901_balanced_quarter_panel()
  data$quarter[[2]] <- 5L

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "invalid values: 5|Expected integer values",
    class = "lwdid_invalid_parameter"
  )
})

test_that("TC-9.1.7: min_global_pre_periods rejects globally short pre windows", {
  data <- make_e901_balanced_quarter_panel()

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 6L
    ),
    regexp = "requires at least 6 pre-treatment period\\(s\\)|Found: 5",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.12: units with zero pre-period rows raise a hard error", {
  data <- make_e901_unit_without_pre_periods_panel()

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "Unit 1 has only 0 pre-period observation\\(s\\)",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.13: units below min_global_pre_periods raise a hard error", {
  data <- make_e901_short_pre_period_unit_panel()

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 2L
    ),
    regexp = "Unit 1 has only 1 pre-period observation\\(s\\)\\. rolling\\('demeanq'\\) requires at least 2",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.14: unit_pre_count uses all pre rows rather than only valid outcomes", {
  data <- make_e901_na_pre_outcome_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 4L)
})

test_that("TC-9.1.17: NA seasonal values are ignored by range validation", {
  data <- make_e901_na_season_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$season_coverage$all_covered, TRUE)
})

test_that("TC-9.1.19: zero-coded seasonal values ask users to recode to 1:Q", {
  data <- make_e901_balanced_quarter_panel()
  data$quarter[[1]] <- 0L

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "recode them to 1-4|invalid values: 0",
    class = "lwdid_invalid_parameter"
  )
})

test_that("TC-9.1.20: integer-like floating seasonal codes are accepted", {
  data <- make_e901_balanced_quarter_panel()
  data$quarter <- as.numeric(data$quarter)

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 5L)
  expect_true(result$season_coverage$all_covered)
})

test_that("TC-9.1.21: exclude_pre_periods beyond the pre window leaves the seasonal gate unchanged", {
  data <- make_e901_balanced_quarter_panel()
  baseline <- lwdid:::.validate_seasonal_inputs(
    data = data,
    y = "y",
    season_var = "quarter",
    Q = 4L,
    ivar = "id",
    tindex = "tindex",
    post = "post",
    method = "demeanq",
    exclude_pre_periods = 0L,
    min_global_pre_periods = 1L
  )

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 10L,
      min_global_pre_periods = 1L
    )
  )

  expect_identical(result$K, baseline$K)
  expect_identical(result$n_units_valid, baseline$n_units_valid)
  expect_identical(result$season_coverage, baseline$season_coverage)
})

test_that("TC-9.1.22: season coverage emits no warning when there are no post rows", {
  data <- make_e901_balanced_quarter_panel()
  data$post <- 0L

  expect_silent(
    result <- lwdid:::.check_season_coverage(
      data = data,
      season_var = "quarter",
      ivar = "id",
      post = "post",
      Q = 4L
    )
  )

  expect_true(result$all_covered)
  expect_identical(result$total_uncovered, 0L)
  expect_length(result$warnings, 0L)
})

test_that("TC-9.1.23: .check_season_coverage() returns structured gap details", {
  data <- make_e901_sparse_season_panel()

  expect_warning(
    result <- lwdid:::.check_season_coverage(
      data = data,
      season_var = "quarter",
      ivar = "id",
      post = "post",
      Q = 4L
    ),
    class = "lwdid_data"
  )

  expect_false(result$all_covered)
  expect_identical(result$total_uncovered, 1L)
  expect_identical(result$units_with_gaps[["1"]], 4)
  expect_length(result$warnings, 1L)
})

test_that("TC-9.1.10/TC-9.1.29: detrendq valid_mask excludes missing pre-treatment tindex rows", {
  data <- make_e901_detrendq_missing_tindex_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "detrendq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 2L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 5L)
})

test_that("TC-9.1.30: detrendq requires at least two distinct pre-treatment periods globally", {
  data <- make_e901_detrendq_one_pre_period_panel()

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "detrendq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 2L
    ),
    regexp = "requires at least 2 pre-treatment period\\(s\\)|Found: 1",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.31/9.1.32: validate_inputs consumes unified seasonal gate for common-timing demeanq", {
  data <- make_e901_validate_inputs_panel()

  expect_warning(
    result <- lwdid:::validate_inputs(
      data = data,
      y = "y",
      d = "d",
      post = "post",
      ivar = "id",
      tvar = "time",
      rolling = "demeanq",
      season_var = "quarter",
      Q = 4L
    )
  )

  expect_identical(result$validated_params$rolling, "demeanq")
  expect_identical(result$n_pre_periods, 4L)
  expect_identical(result$K, 4L)
})

test_that("TC-9.1.31: transform_common() executes demeanq after seasonal gate passes", {
  data <- make_e901_transform_panel()

  result <- suppressWarnings(
    lwdid:::transform_common(
      data.table::as.data.table(data),
      y = "y",
      ivar = "id",
      tvar = "time",
      g = 6L,
      rolling = "demeanq",
      post = "post",
      season_var = "quarter",
      Q = 4L
    ),
  )

  expect_equal(
    result[result$id == 1L, "y_trans"][[1]],
    c(0, 0, 0, 0, 0, 6, 4)
  )
  expect_equal(
    result[result$id == 3L, "y_trans"][[1]],
    c(0, 0, 0, 0, 0, 1, -1)
  )
})

test_that("TC-9.1.32: transform_common() blocks demeanq when seasonal gate fails", {
  data <- make_e901_transform_gate_fail_panel()

  expect_error(
    lwdid:::transform_common(
      data.table::as.data.table(data),
      y = "y",
      ivar = "id",
      tvar = "time",
      g = 3L,
      rolling = "demeanq",
      post = "post",
      season_var = "quarter",
      Q = 4L
    ),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.31: lwdid common-timing demeanq runs through transform_common() with season_var metadata", {
  data <- make_e901_transform_panel()

  result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "demeanq",
        season_var = "quarter",
        Q = 4L
      )
    )
  )

  expect_s3_class(result, "lwdid_result")
  expect_equal(result$att, 5, tolerance = 1e-10)
})

test_that("TC-9.1.28/TC-9.1.31: transform_common() executes detrendq after seasonal gate passes", {
  data <- make_e901_detrendq_panel()

  result <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expect_equal(
    result[result$id == 1L, "y_trans"][[1]],
    c(rep(0, 6), rep(4, 4)),
    tolerance = 1e-8
  )
  expect_equal(
    result[result$id == 3L, "y_trans"][[1]],
    rep(0, 10),
    tolerance = 1e-8
  )
})

test_that("TC-9.1.31: lwdid common-timing detrendq runs through transform_common() with season_var metadata", {
  data <- make_e901_detrendq_panel()

  result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "detrendq",
        season_var = "quarter",
        Q = 4L
      )
    )
  )

  expect_s3_class(result, "lwdid_result")
  expect_equal(result$att, 4, tolerance = 1e-10)
  expect_equal(result$K, 6L)
  expect_equal(result$tpost1, 7L)
})

test_that("TC-9.1.33: validated seasonal K propagates to the top-level common-timing result", {
  data <- make_e901_detrendq_panel()
  validated <- lwdid:::.validate_seasonal_inputs(
    data = data,
    y = "y",
    season_var = "quarter",
    Q = 4L,
    ivar = "id",
    tindex = "time",
    post = "post",
    method = "detrendq",
    exclude_pre_periods = 0L,
    min_global_pre_periods = 2L
  )

  result <- suppressWarnings(
    suppressMessages(
      lwdid(
        data = data,
        y = "y",
        ivar = "id",
        tvar = "time",
        d = "d",
        post = "post",
        rolling = "detrendq",
        season_var = "quarter",
        Q = 4L
      )
    )
  )

  expect_identical(result$K, validated$K)
  expect_identical(result$tpost1, validated$K + 1L)
})

test_that("story-E9-01 quarterly public result surface includes average row and quarter-formatted labels", {
  data_path <- normalizePath(
    testthat::test_path("../../../lwdid-py_v0.2.3/tests/data/smoking_quarterly.csv"),
    mustWork = FALSE
  )
  skip_if_not(file.exists(data_path),
              "External quarterly smoking data not available")
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)

  for (rolling in c("demeanq", "detrendq")) {
    result <- suppressWarnings(
      suppressMessages(
        lwdid(
          data = panel,
          y = "lcigsale",
          ivar = "state",
          tvar = c("year", "quarter"),
          d = "d",
          post = "post",
          rolling = rolling,
          vce = NULL
        )
      )
    )

    expect_identical(nrow(result$att_by_period), 49L, info = rolling)
    expect_identical(as.character(result$att_by_period$period[[1]]), "average", info = rolling)
    expect_identical(as.character(result$att_by_period$tindex[[1]]), "-", info = rolling)
    expect_identical(as.character(result$att_by_period$period[[2]]), "1989q1", info = rolling)
    expect_identical(
      as.character(result$att_by_period$period[[nrow(result$att_by_period)]]),
      "2000q4",
      info = rolling
    )
    expect_equal(result$att_by_period$att[[1]], result$att, tolerance = 1e-10, info = rolling)
  }
})

test_that("TC-9.1.8: transform_common() detrendq accepts quarter alias when season_var is omitted", {
  data <- make_e901_detrendq_panel()

  baseline <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )
  aliased <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 7L,
    rolling = "detrendq",
    post = "post",
    quarter = "quarter",
    Q = 4L
  )

  expect_equal(aliased$y_trans, baseline$y_trans, tolerance = 1e-10)
  expect_equal(aliased$n_pre, baseline$n_pre)
})

test_that("TC-9.1.6: transform_common() detrendq honors exclude_pre_periods when enough pre-periods remain", {
  data <- make_e901_detrendq_exclude_panel()

  result <- lwdid:::transform_common(
    data.table::as.data.table(data),
    y = "y",
    ivar = "id",
    tvar = "time",
    g = 8L,
    rolling = "detrendq",
    post = "post",
    season_var = "quarter",
    Q = 4L,
    exclude_pre_periods = 1L
  )

  expect_equal(
    result[result$id == 1L, "y_trans"][[1]],
    c(rep(0, 7), rep(4, 4)),
    tolerance = 1e-10
  )
  expect_equal(
    result[result$id == 3L, "y_trans"][[1]],
    rep(0, 11),
    tolerance = 1e-10
  )
  expect_true(all(result[result$id == 1L, "n_pre"][[1]] == 6L))
})

test_that("TC-9.1.18: all pre-treatment tindex values missing raises K-specific hard error", {
  data <- make_e901_all_na_pre_tindex_panel()

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "all pre-treatment tindex values are missing|valid time index values",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.5: .auto_detect_frequency() returns quarterly contract for repeated yearly quarters", {
  data <- data.frame(
    id = rep(1:3, each = 20),
    year = rep(rep(2001:2005, each = 4), 3),
    quarter = rep(rep(1:4, 5), 3),
    y = rnorm(60)
  )

  result <- lwdid:::.auto_detect_frequency(
    data = data,
    tvar = "year",
    ivar = "id"
  )

  expect_true(is.list(result))
  expect_identical(result$frequency, "quarterly")
  expect_identical(result$Q, 4L)
  expect_gt(result$confidence, 0.5)
})

test_that("TC-9.1.15: .auto_detect_frequency() returns NULL for ambiguous consecutive integer index", {
  data <- data.frame(
    id = rep(1:2, each = 10),
    t = rep(1:10, 2),
    y = rnorm(20)
  )

  expect_warning(
    result <- lwdid:::.auto_detect_frequency(
      data = data,
      tvar = "t",
      ivar = "id"
    ),
    class = "lwdid_data"
  )

  expect_null(result)
})

test_that("TC-9.1.24/TC-9.1.26: monthly seasonal validation and frequency detection are supported", {
  validation_data <- make_e901_monthly_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = validation_data,
      y = "y",
      season_var = "month",
      Q = 12L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 13L)

  detection_data <- data.frame(
    id = rep(1:2, each = 24),
    year = rep(rep(2001:2002, each = 12), 2),
    month = rep(rep(1:12, 2), 2),
    y = rnorm(48)
  )

  detection <- lwdid:::.auto_detect_frequency(
    data = detection_data,
    tvar = "year",
    ivar = "id"
  )

  expect_identical(detection$frequency, "monthly")
  expect_identical(detection$Q, 12L)
})

test_that("TC-9.1.25/TC-9.1.27: weekly seasonal validation and frequency detection are supported", {
  validation_data <- make_e901_weekly_panel()

  expect_silent(
    result <- lwdid:::.validate_seasonal_inputs(
      data = validation_data,
      y = "y",
      season_var = "week",
      Q = 52L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 53L)

  detection_data <- data.frame(
    id = rep(1:2, each = 104),
    year = rep(rep(2001:2002, each = 52), 2),
    week = rep(rep(1:52, 2), 2),
    y = rnorm(208)
  )

  detection <- lwdid:::.auto_detect_frequency(
    data = detection_data,
    tvar = "year",
    ivar = "id"
  )

  expect_identical(detection$frequency, "weekly")
  expect_identical(detection$Q, 52L)
})

test_that("TC-9.1.16: missing season_var column raises explicit error", {
  data <- make_e901_sparse_season_panel()
  names(data)[names(data) == "quarter"] <- "quarter_alt"

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "Required column 'quarter' not found",
    class = "lwdid_missing_column"
  )
})

test_that("TC-9.1.17: NA seasonal values are ignored during range validation", {
  data <- data.frame(
    id = rep(1L, 5L),
    tindex = 1:5,
    y = c(10, 11, 12, 13, 14),
    quarter = c(1L, NA, 2L, 3L, 1L),
    post = c(0L, 0L, 0L, 0L, 1L)
  )

  expect_no_warning(
    result <- lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    )
  )

  expect_true(isTRUE(result$is_valid))
  expect_identical(result$K, 4L)
  expect_identical(result$season_coverage$total_uncovered, 0L)
})

test_that("TC-9.1.19: zero-coded seasonal values prompt recoding hint", {
  data <- make_e901_sparse_season_panel()
  data$quarter[1] <- 0L

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "demeanq",
      exclude_pre_periods = 0L,
      min_global_pre_periods = 1L
    ),
    regexp = "recode them to 1-4",
    class = "lwdid_invalid_parameter"
  )
})

test_that("TC-9.1.23: .check_season_coverage() returns structured coverage summary", {
  data <- make_e901_sparse_season_panel()

  expect_warning(
    coverage <- lwdid:::.check_season_coverage(
      data = data,
      season_var = "quarter",
      ivar = "id",
      post = "post",
      Q = 4L
    ),
    class = "lwdid_data"
  )

  expect_true(is.list(coverage))
  expect_named(
    coverage,
    c("all_covered", "units_with_gaps", "total_uncovered", "warnings")
  )
  expect_identical(coverage$all_covered, FALSE)
  expect_identical(coverage$total_uncovered, 1L)
  expect_identical(coverage$units_with_gaps[["1"]], 4)
  expect_length(coverage$warnings, 1L)
})

test_that("TC-9.1.30: detrendq helper defaults min_global_pre_periods to two", {
  data <- data.frame(
    id = c(1L, 1L, 1L, 1L),
    tindex = c(1L, 1L, 1L, 2L),
    y = c(10, 11, 12, 13),
    quarter = c(1L, 1L, 1L, 1L),
    post = c(0L, 0L, 0L, 1L)
  )

  expect_error(
    lwdid:::.validate_seasonal_inputs(
      data = data,
      y = "y",
      season_var = "quarter",
      Q = 4L,
      ivar = "id",
      tindex = "tindex",
      post = "post",
      method = "detrendq",
      exclude_pre_periods = 0L
    ),
    regexp = "requires at least 2 pre-treatment period\\(s\\)",
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("TC-9.1.34: seasonal parameter requirements match demeanq seasonal-dummy counts", {
  contract <- lwdid:::.seasonal_parameter_requirements(
    method = "demeanq",
    n_unique_seasons = 4L
  )

  expect_identical(contract$method, "demeanq")
  expect_identical(contract$n_unique_seasons, 4L)
  expect_identical(contract$n_params, 4L)
  expect_identical(contract$min_required, 5L)
})

test_that("TC-9.1.35: seasonal parameter requirements add one trend term for detrendq", {
  contract <- lwdid:::.seasonal_parameter_requirements(
    method = "detrendq",
    n_unique_seasons = 4L
  )

  expect_identical(contract$method, "detrendq")
  expect_identical(contract$n_unique_seasons, 4L)
  expect_identical(contract$n_params, 5L)
  expect_identical(contract$min_required, 6L)
})

test_that("TC-9.1.36: seasonal parameter requirements enforce the df >= 1 boundary", {
  contract <- lwdid:::.seasonal_parameter_requirements(
    method = "demeanq",
    n_unique_seasons = 3L
  )

  expect_identical(contract$n_params, 3L)
  expect_identical(contract$min_required, 4L)
  expect_identical(contract$min_required - contract$n_params, 1L)
})

test_that("TC-9.1.34: vibe-math confirms demeanq seasonal parameter count equals observed pre-period seasons", {
  oracle <- read_e901_required_parity_json(
    "20260328-theory-parity-e9-01-task6-vibemath.json"
  )

  case <- oracle$cases$tc_9_1_34

  expect_true(isTRUE(case$matched))
  expect_identical(case$comparison$status, "matched")
  expect_equal(case$manual_expected, 4)
  expect_equal(case$vibemath_result, 4)
  expect_equal(case$variables$n_seasons, 4)
})

test_that("TC-9.1.35: vibe-math confirms detrendq adds exactly one trend parameter beyond demeanq", {
  oracle <- read_e901_required_parity_json(
    "20260328-theory-parity-e9-01-task6-vibemath.json"
  )

  case <- oracle$cases$tc_9_1_35

  expect_true(isTRUE(case$matched))
  expect_identical(case$comparison$status, "matched")
  expect_equal(case$manual_expected, 5)
  expect_equal(case$vibemath_result, 5)
  expect_equal(case$variables$n_seasons, 4)
})

test_that("TC-9.1.36: vibe-math confirms min_required = n_params + 1 is the df >= 1 gate", {
  oracle <- read_e901_required_parity_json(
    "20260328-theory-parity-e9-01-task6-vibemath.json"
  )

  case <- oracle$cases$tc_9_1_36
  checks <- case$checks

  expect_identical(case$comparison$status, "matched")
  expect_true(all(vapply(checks, function(check) isTRUE(check$matched), logical(1))))
  expect_equal(checks[[1]]$vibemath_result, 5)
  expect_equal(checks[[2]]$vibemath_result, 6)
  expect_equal(checks[[3]]$vibemath_result, 1)
  expect_equal(checks[[4]]$vibemath_result, 0)
})
