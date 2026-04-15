# test-demeanq.R
# Story E9-02 dedicated demeanq helper/package-surface regressions

resolve_demeanq_parity_path <- function(filename) {
  if (exists("resolve_parity_fixture_path", mode = "function")) {
    resolved <- resolve_parity_fixture_path(filename)
    if (file.exists(resolved)) {
      return(resolved)
    }
  }

  candidates <- c(
    testthat::test_path(
      "..", "..", "..",
      "_automation", "test-artifacts", "parity", filename
    ),
    file.path(
      "/Users/cxy/Desktop/lwdid_r/_automation/test-artifacts/parity",
      filename
    ),
    testthat::test_path("_fixtures", "parity", filename)
  )

  candidates <- unique(normalizePath(candidates, mustWork = FALSE))
  existing <- candidates[file.exists(candidates)]

  if (length(existing) > 0L) {
    return(existing[[1L]])
  }

  candidates[[1L]]
}

read_demeanq_parity_oracle <- function(filename) {
  jsonlite::fromJSON(
    resolve_demeanq_parity_path(filename),
    simplifyVector = FALSE
  )
}

read_demeanq_parity_fixture <- function(filename) {
  utils::read.csv(
    resolve_demeanq_parity_path(filename),
    stringsAsFactors = FALSE
  )
}

run_demeanq_parity_producer <- function(filename) {
  script_path <- resolve_demeanq_parity_path(filename)
  output <- system2("Rscript", script_path, stdout = TRUE, stderr = TRUE)
  exit_code <- attr(output, "status")
  if (is.null(exit_code)) {
    exit_code <- 0L
  }

  list(
    script = script_path,
    output = output,
    exit_code = as.integer(exit_code)
  )
}

test_that("TC-9.2.1, TC-9.2.3: .demeanq_unit removes a known seasonal mean pattern", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 9L),
    tindex = seq_len(9L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    y = c(10, 12, 13, 14, 10, 17, 18, 19, 15)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_equal(
    result$yhat,
    c(10, 12, 13, 14, 10, 12, 13, 14, 10),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot,
    c(rep.int(0, 5L), rep.int(5, 4L)),
    tolerance = 1e-12
  )
  expect_identical(result$n_pre, 5L)
  expect_equal(mean(result$ydot[unit_data$post == 0L]), 0, tolerance = 1e-12)
})

test_that("TC-9.2.2: .demeanq_unit matches the frozen Python helper oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_demeanq_parity_path("e9_02_layer2_python_helper_oracle.json")
  fixture_path <- resolve_demeanq_parity_path("e9_02_demeanq_helper_fixture.csv")

  expect_true(
    file.exists(oracle_path),
    info = paste("missing demeanq python helper oracle:", oracle_path)
  )
  expect_true(
    file.exists(fixture_path),
    info = paste("missing demeanq helper fixture:", fixture_path)
  )

  if (!file.exists(oracle_path) || !file.exists(fixture_path)) {
    return(invisible(NULL))
  }

  oracle <- read_demeanq_parity_oracle("e9_02_layer2_python_helper_oracle.json")
  fixture <- data.table::as.data.table(
    read_demeanq_parity_fixture("e9_02_demeanq_helper_fixture.csv")
  )

  expect_identical(oracle$story, "story-E9-02")
  expect_identical(oracle$layer, "layer-2")
  expect_identical(oracle$exact_status, "passed")
  expect_identical(oracle$numeric_status, "passed")
  expect_identical(
    oracle$blocker_boundary,
    "no-python-helper-oracle-blocker"
  )

  case <- oracle$cases$known_pattern
  result <- lwdid:::.demeanq_unit(
    unit_data = fixture,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_identical(case$fixture_rows, 9L)
  expect_identical(case$n_pre, 5L)
  expect_equal(
    result$yhat,
    unname(as.numeric(unlist(case$python_result$yhat, use.names = FALSE))),
    tolerance = 1e-10
  )
  expect_equal(
    result$ydot,
    unname(as.numeric(unlist(case$python_result$ydot, use.names = FALSE))),
    tolerance = 1e-10
  )
})

test_that("TC-9.2.22: .demeanq_unit exposes coefficients that match the vibe-math oracle", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_demeanq_parity_path("20260328-theory-parity-e9-02-task22-vibemath.json")
  expect_true(
    file.exists(oracle_path),
    info = paste("missing demeanq vibe-math oracle:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_demeanq_parity_oracle("20260328-theory-parity-e9-02-task22-vibemath.json")
  oracle_values <- stats::setNames(
    vapply(oracle$checks, function(item) as.numeric(item$manual_expected), numeric(1)),
    vapply(oracle$checks, function(item) item$id, character(1))
  )

  unit_data <- data.table::data.table(
    id = rep.int(1L, 9L),
    tindex = seq_len(9L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    y = c(10, 12, 13, 14, 10, 17, 18, 19, 15)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_identical(oracle$comparison$status, "matched")
  expect_identical(oracle$comparison$story_local_test_status, "oracle-ready-but-not-yet-wired")
  expect_equal(
    result$coefficients,
    c(
      "(Intercept)" = oracle_values[["TC-9.2.22-mu-hat"]],
      "season_2" = oracle_values[["TC-9.2.22-gamma-q2"]],
      "season_3" = oracle_values[["TC-9.2.22-gamma-q3"]],
      "season_4" = oracle_values[["TC-9.2.22-gamma-q4"]]
    ),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot[6L],
    oracle_values[["TC-9.2.22-post-residual-q2"]],
    tolerance = 1e-12
  )
  expect_equal(
    mean(result$ydot[unit_data$post == 1L]),
    oracle_values[["TC-9.2.22-postavg"]],
    tolerance = 1e-12
  )
})

test_that("TC-9.2.4, TC-9.2.7: .demeanq_unit extrapolates unseen post seasons from the intercept", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 7L),
    tindex = seq_len(7L),
    quarter = c(1L, 2L, 1L, 2L, 1L, 3L, 4L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L),
    y = c(10, 12, 10, 12, 10, 16, 17)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  post_idx <- which(unit_data$post == 1L)
  expect_equal(result$yhat[post_idx], c(10, 10), tolerance = 1e-12)
  expect_equal(result$ydot[post_idx], c(6, 7), tolerance = 1e-12)
  expect_false(any(is.na(result$ydot[post_idx])))
})

test_that("TC-9.2.8: .demeanq_unit supports monthly seasonality with Q = 12", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 16L),
    tindex = seq_len(16L),
    month = c(1:12, 1L, 2L, 3L, 4L),
    post = c(rep.int(0L, 13L), rep.int(1L, 3L)),
    y = c(10:21, 10, 16, 17, 18)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "month",
    post = "post",
    Q = 12L
  )

  expect_identical(result$n_pre, 13L)
  expect_equal(
    result$yhat[c(1L, 12L, 14L, 16L)],
    c(10, 21, 11, 13),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot[unit_data$post == 0L],
    rep.int(0, 13L),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot[unit_data$post == 1L],
    rep.int(5, 3L),
    tolerance = 1e-12
  )
})

test_that("TC-9.2.9: .demeanq_unit supports weekly seasonality with Q = 52", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 56L),
    tindex = seq_len(56L),
    week = c(1:52, 1L, 2L, 3L, 4L),
    post = c(rep.int(0L, 53L), rep.int(1L, 3L)),
    y = c(101:152, 101, 109, 110, 111)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "week",
    post = "post",
    Q = 52L
  )

  expect_identical(result$n_pre, 53L)
  expect_equal(
    result$yhat[c(1L, 52L, 54L, 56L)],
    c(101, 152, 102, 104),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot[unit_data$post == 0L],
    rep.int(0, 53L),
    tolerance = 1e-12
  )
  expect_equal(
    result$ydot[unit_data$post == 1L],
    rep.int(7, 3L),
    tolerance = 1e-12
  )
})

test_that("TC-9.2.5: .demeanq_unit preserves the current small-sample warning and NA surface", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 4L),
    tindex = seq_len(4L),
    quarter = c(1L, 2L, 3L, 4L),
    post = c(0L, 0L, 0L, 1L),
    y = c(10, 12, 13, 15)
  )

  expect_warning(
    result <- lwdid:::.demeanq_unit(
      unit_data = unit_data,
      y = "y",
      season_var = "quarter",
      post = "post",
      Q = 4L
    ),
    class = "lwdid_small_sample"
  )

  expect_true(all(is.na(result$yhat)))
  expect_true(all(is.na(result$ydot)))
  expect_identical(result$n_pre, 3L)
})

test_that("TC-9.2.6: .demeanq_unit degrades OLS failures into the NA warning surface", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 6L),
    tindex = seq_len(6L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 1L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L),
    y = c(10, Inf, 13, 14, 10, 15)
  )

  expect_warning(
    result <- lwdid:::.demeanq_unit(
      unit_data = unit_data,
      y = "y",
      season_var = "quarter",
      post = "post",
      Q = 4L
    ),
    class = "lwdid_data"
  )

  expect_true(all(is.na(result$yhat)))
  expect_true(all(is.na(result$ydot)))
  expect_identical(result$n_pre, 5L)
})

test_that("TC-9.2.21: .demeanq_unit returns NA for rows with missing season values", {
  unit_data <- data.table::data.table(
    id = rep.int(1L, 8L),
    tindex = seq_len(8L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L, NA, 4L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L),
    y = c(10, 12, 13, 14, 10, 17, 18, 19)
  )

  result <- lwdid:::.demeanq_unit(
    unit_data = unit_data,
    y = "y",
    season_var = "quarter",
    post = "post",
    Q = 4L
  )

  expect_equal(result$yhat[6L], 12, tolerance = 1e-12)
  expect_equal(result$ydot[6L], 5, tolerance = 1e-12)
  expect_true(is.na(result$yhat[7L]))
  expect_true(is.na(result$ydot[7L]))
  expect_equal(result$yhat[8L], 14, tolerance = 1e-12)
  expect_equal(result$ydot[8L], 5, tolerance = 1e-12)
})

test_that("TC-9.2.10, TC-9.2.11: .demeanq_transform keeps the seasonal surface and honors exclude_pre_periods", {
  panel <- data.table::data.table(
    id = rep(1:2, each = 8L),
    tindex = rep(seq_len(8L), 2L),
    quarter = rep(c(1L, 2L, 3L, 4L, 1L, 2L, 1L, 2L), 2L),
    post = rep(c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L), 2L),
    y = c(
      10, 12, 13, 14, 11, 15, 16, 18,
      20, 22, 23, 24, 21, 25, 26, 28
    )
  )

  baseline <- lwdid:::.demeanq_transform(
    data = panel,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L,
    exclude_pre_periods = 0L
  )
  excluded <- lwdid:::.demeanq_transform(
    data = panel,
    y = "y",
    ivar = "id",
    tindex = "tindex",
    season_var = "quarter",
    post = "post",
    Q = 4L,
    exclude_pre_periods = 1L
  )

  expect_true(all(c("seasonal_fit", "y_trans", "n_pre") %in% names(excluded)))
  expect_false(".seasonal_post" %in% names(excluded))
  expect_false(isTRUE(all.equal(baseline$y_trans, excluded$y_trans)))
  expect_identical(unique(excluded[id == 1L]$n_pre), 5L)
  expect_identical(unique(excluded[id == 2L]$n_pre), 5L)
})

test_that("TC-9.2.12, TC-9.2.16: transform_common exposes the demeanq seasonal_fit surface", {
  panel <- data.table::data.table(
    id = rep.int(1L, 9L),
    tindex = seq_len(9L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    y = c(10, 12, 13, 14, 10, 17, 18, 19, 15)
  )

  result <- lwdid:::transform_common(
    dt = panel,
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    rolling = "demeanq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )

  expect_true(all(c("seasonal_fit", "y_trans", "n_pre") %in% names(result)))
  expect_equal(
    result[id == 1L, seasonal_fit],
    c(10, 12, 13, 14, 10, 12, 13, 14, 10),
    tolerance = 1e-12
  )
  expect_equal(
    result[id == 1L, y_trans],
    c(rep.int(0, 5L), rep.int(5, 4L)),
    tolerance = 1e-12
  )
  expect_identical(unique(result[id == 1L]$n_pre), 5L)
})

test_that("TC-9.2.13: demeanq missing season_var error includes a concrete example", {
  panel <- data.table::data.table(
    id = rep.int(1L, 6L),
    tindex = seq_len(6L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L),
    y = c(10, 12, 13, 14, 10, 17)
  )

  expect_error(
    lwdid:::transform_common(
      dt = panel,
      y = "y",
      ivar = "id",
      tvar = "tindex",
      g = 6L,
      rolling = "demeanq",
      post = "post"
    ),
    regexp = paste0(
      "season_var'.*Example: lwdid\\(",
      "\\.\\.\\., rolling='demeanq', season_var='quarter', Q=4\\)"
    ),
    class = "lwdid_invalid_parameter"
  )
})

test_that("TC-9.2.14: demeanq quarter alias matches season_var", {
  panel <- data.table::data.table(
    id = rep.int(1L, 9L),
    tindex = seq_len(9L),
    quarter = c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L),
    post = c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L),
    y = c(10, 12, 13, 14, 10, 17, 18, 19, 15)
  )

  baseline <- lwdid:::transform_common(
    dt = panel,
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    rolling = "demeanq",
    post = "post",
    season_var = "quarter",
    Q = 4L
  )
  aliased <- lwdid:::transform_common(
    dt = panel,
    y = "y",
    ivar = "id",
    tvar = "tindex",
    g = 6L,
    rolling = "demeanq",
    post = "post",
    quarter = "quarter",
    Q = 4L
  )

  expect_equal(aliased$seasonal_fit, baseline$seasonal_fit, tolerance = 1e-12)
  expect_equal(aliased$y_trans, baseline$y_trans, tolerance = 1e-12)
  expect_identical(aliased$n_pre, baseline$n_pre)
})

test_that("TC-9.2.15: demeanq validator surface rejects out-of-range season codes", {
  skip_if_not(file.exists(resolve_demeanq_parity_path("e9_02_validator_surface_probe.R")),
              "Fixture files not available (excluded from built package)")
  producer <- run_demeanq_parity_producer("e9_02_validator_surface_probe.R")
  expect_equal(producer$exit_code, 0L, info = paste(producer$output, collapse = "\n"))

  payload <- read_demeanq_parity_oracle(
    "20260328-qa-parity-e9-02-validator-surface-probe.json"
  )
  case <- payload$cases$tc_9_2_15

  expect_identical(payload$blocker_boundary, "no-validator-surface-blocker")
  expect_identical(case$status, "error")
  expect_match(case$message, "invalid values: 0, 5", fixed = TRUE)
  expect_match(case$message, "recode them to 1-4", fixed = TRUE)
})

test_that("TC-9.2.19: demeanq validator surface rejects float non-integer season codes", {
  skip_if_not(file.exists(resolve_demeanq_parity_path("e9_02_validator_surface_probe.R")),
              "Fixture files not available (excluded from built package)")
  producer <- run_demeanq_parity_producer("e9_02_validator_surface_probe.R")
  expect_equal(producer$exit_code, 0L, info = paste(producer$output, collapse = "\n"))

  payload <- read_demeanq_parity_oracle(
    "20260328-qa-parity-e9-02-validator-surface-probe.json"
  )
  case <- payload$cases$tc_9_2_19

  expect_identical(case$status, "error")
  expect_match(case$message, "invalid values: 1.5", fixed = TRUE)
})

test_that("TC-9.2.20: demeanq validator surface accepts float integer season codes", {
  skip_if_not(file.exists(resolve_demeanq_parity_path("e9_02_validator_surface_probe.R")),
              "Fixture files not available (excluded from built package)")
  producer <- run_demeanq_parity_producer("e9_02_validator_surface_probe.R")
  expect_equal(producer$exit_code, 0L, info = paste(producer$output, collapse = "\n"))

  payload <- read_demeanq_parity_oracle(
    "20260328-qa-parity-e9-02-validator-surface-probe.json"
  )
  case <- payload$cases$tc_9_2_20

  expect_identical(case$status, "passed")
  expect_equal(
    as.numeric(unlist(case$normalized_values, use.names = FALSE)),
    c(1, 2, 4),
    tolerance = 0
  )
})

test_that("TC-9.2.17: top-level demeanq exclude_pre_periods changes the common-timing seasonal surface", {
  skip_if_not(file.exists(resolve_demeanq_parity_path("e9_02_top_level_surface_probe.R")),
              "Fixture files not available (excluded from built package)")
  producer <- run_demeanq_parity_producer("e9_02_top_level_surface_probe.R")
  expect_equal(producer$exit_code, 0L, info = paste(producer$output, collapse = "\n"))

  payload <- read_demeanq_parity_oracle(
    "20260328-story-worker-e9-02-top-level-surface-probe.json"
  )
  case <- payload$cases$tc_9_2_17

  expect_identical(payload$story, "story-E9-02")
  expect_identical(payload$role, "story-worker")
  expect_identical(payload$scenario, "top-level-surface-probe")
  expect_identical(
    payload$blocker_boundary,
    "no-top-level-contract-blocker"
  )

  expect_identical(case$transform_common$baseline$exclude_pre_periods, 0L)
  expect_identical(case$transform_common$excluded$exclude_pre_periods, 1L)
  expect_true(isTRUE(case$transform_common$changed$seasonal_fit))
  expect_true(isTRUE(case$transform_common$changed$y_trans))
  expect_true(isTRUE(case$transform_common$changed$n_pre))
  expect_equal(
    as.numeric(unlist(case$transform_common$baseline$unit_1$seasonal_fit, use.names = FALSE)),
    c(10.5, 13.5, 13, 14, 10.5, 13.5, 10.5, 13.5),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(unlist(case$transform_common$excluded$unit_1$seasonal_fit, use.names = FALSE)),
    c(10.5, 12, 13, 14, 10.5, 12, 10.5, 12),
    tolerance = 1e-12
  )
  expect_equal(
    as.numeric(unlist(case$transform_common$excluded$unit_1$y_trans, use.names = FALSE)),
    c(-0.5, 0, 0, 0, 0.5, 3, 5.5, 6),
    tolerance = 1e-12
  )
  expect_identical(case$transform_common$baseline$unit_1$n_pre, 6L)
  expect_identical(case$transform_common$excluded$unit_1$n_pre, 5L)

  expect_identical(case$lwdid$baseline$exclude_pre_periods, 0L)
  expect_identical(case$lwdid$excluded$exclude_pre_periods, 1L)
  expect_true(is.finite(case$lwdid$baseline$att))
  expect_true(is.finite(case$lwdid$excluded$att))
})

test_that("TC-9.2.18: top-level uncovered post seasons warn and pass while logging diagnostics", {
  skip_if_not(file.exists(resolve_demeanq_parity_path("e9_02_top_level_surface_probe.R")),
              "Fixture files not available (excluded from built package)")
  producer <- run_demeanq_parity_producer("e9_02_top_level_surface_probe.R")
  expect_equal(producer$exit_code, 0L, info = paste(producer$output, collapse = "\n"))

  payload <- read_demeanq_parity_oracle(
    "20260328-story-worker-e9-02-top-level-surface-probe.json"
  )
  case <- payload$cases$tc_9_2_18

  expect_identical(case$seasonal_gate$status, "warning+pass")
  expect_false(isTRUE(case$seasonal_gate$season_coverage$all_covered))
  expect_identical(case$seasonal_gate$season_coverage$total_uncovered, 3L)
  expect_length(case$seasonal_gate$warnings, 3L)
  expect_true(all(vapply(
    case$seasonal_gate$warnings,
    function(message) grepl("allows the seasonal gate to pass", message, fixed = TRUE),
    logical(1)
  )))

  expect_identical(case$lwdid$status, "warning+pass")
  expect_identical(case$lwdid$warning_diagnostics_count, 1L)
  expect_identical(case$lwdid$warnings_log_count, 3L)
  expect_identical(case$lwdid$warning_diagnostics[[1L]]$category, "lwdid_data")
  expect_identical(case$lwdid$warning_diagnostics[[1L]]$count, 3L)
  expect_true(all(vapply(
    case$lwdid$warnings_log,
    function(record) identical(record$category, "lwdid_data"),
    logical(1)
  )))
  expect_true(all(vapply(
    case$lwdid$warnings_log,
    function(record) grepl("allows the seasonal gate to pass", record$message, fixed = TRUE),
    logical(1)
  )))

  expect_identical(case$contextual_comparator$python_public_api, "hard-error")
  expect_identical(case$contextual_comparator$stata_public_api, "hard-error")
})
