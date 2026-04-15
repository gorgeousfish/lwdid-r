library(testthat)

if (!exists("lwdid", mode = "function")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run test-estimator-realdata-parity-e606.R directly.")
  }

  devtools::load_all(
    "/Users/cxy/Desktop/lwdid_r/lwdid-r",
    export_all = FALSE,
    quiet = TRUE
  )
}

read_e606_realdata_oracle <- function(filename) {
  jsonlite::fromJSON(
    resolve_parity_fixture_path(filename),
    simplifyVector = FALSE
  )
}

normalize_e606_group_time_rows <- function(df) {
  cols <- c("cohort", "period", "event_time", "att", "se", "pvalue")
  if (is.null(df) || !nrow(df)) {
    out <- as.data.frame(
      setNames(replicate(length(cols), numeric(0), simplify = FALSE), cols)
    )
    return(out)
  }

  out <- df[, cols, drop = FALSE]
  out <- out[order(out$cohort, out$period), , drop = FALSE]
  out$cohort <- as.integer(out$cohort)
  out$period <- as.integer(out$period)
  out$event_time <- as.integer(out$event_time)
  out$att <- as.numeric(out$att)
  out$se <- as.numeric(out$se)
  out$pvalue <- as.numeric(out$pvalue)
  rownames(out) <- NULL
  out
}

oracle_group_time_rows <- function(oracle, estimator) {
  rows <- oracle$estimators[[estimator]]$shared_group_time_rows
  if (is.null(rows) || length(rows) == 0L) {
    return(normalize_e606_group_time_rows(NULL))
  }

  if (is.data.frame(rows)) {
    return(normalize_e606_group_time_rows(rows))
  }

  build_row <- function(row) {
    data.frame(
      cohort = as.integer(row$cohort),
      period = as.integer(row$period),
      event_time = as.integer(row$event_time),
      att = if (is.null(row$att)) NA_real_ else as.numeric(row$att),
      se = if (is.null(row$se)) NA_real_ else as.numeric(row$se),
      pvalue = if (is.null(row$pvalue)) NA_real_ else as.numeric(row$pvalue),
      stringsAsFactors = FALSE
    )
  }

  normalize_e606_group_time_rows(
    do.call(rbind, lapply(rows, build_row))
  )
}

current_e606_castle_result <- function(estimator, aggregate = "none") {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)

  args <- list(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = estimator,
    aggregate = aggregate
  )

  if (estimator != "ra") {
    args$controls <- c("income", "unemployrt", "poverty")
  }
  if (identical(estimator, "psm")) {
    args$caliper <- 0.2
  }

  warnings_seen <- character(0)
  result <- withCallingHandlers(
    do.call(lwdid, args),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(
    result = result,
    warnings = unique(warnings_seen)
  )
}

expect_numeric_or_na_equal <- function(actual, expected, tolerance, info = NULL) {
  expect_equal(length(actual), length(expected), info = info)

  for (idx in seq_along(actual)) {
    if (is.na(actual[[idx]]) && is.na(expected[[idx]])) {
      next
    }

    expect_equal(
      actual[[idx]],
      expected[[idx]],
      tolerance = tolerance,
      info = info
    )
  }
}

test_that("E6-06 Layer 3 castle oracle exists and freezes the current parity boundary", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_parity_fixture_path(
    "e6_06_layer3_castle_realdata_comparator.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing E6-06 Layer 3 real-data comparator:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_e606_realdata_oracle("e6_06_layer3_castle_realdata_comparator.json")

  expect_identical(oracle$story, "story-E6-06")
  expect_identical(oracle$layer, "layer-3")
  expect_identical(oracle$scenario, "castle-staggered-shared-realdata")
  expect_identical(oracle$group_time_status, "partial-drift-confirmed")
  expect_identical(oracle$overall_aggregation_status, "drift-confirmed")
  expect_identical(
    oracle$blocker_boundary,
    "castle-shared-ipw-group-time-and-overall-aggregation-drift"
  )
})

test_that("E6-06 Layer 3 castle group-time rows stay aligned with shared comparator", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_parity_fixture_path(
    "e6_06_layer3_castle_realdata_comparator.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing E6-06 Layer 3 real-data comparator:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_e606_realdata_oracle("e6_06_layer3_castle_realdata_comparator.json")

  for (estimator in c("ipwra", "psm")) {
    current <- current_e606_castle_result(estimator = estimator, aggregate = "none")
    current_rows <- normalize_e606_group_time_rows(current$result$att_by_cohort_time)
    expected_rows <- oracle_group_time_rows(oracle, estimator)

    expect_equal(
      nrow(current_rows),
      nrow(expected_rows),
      info = paste("row-count drift for", estimator)
    )

    if (!nrow(expected_rows)) {
      next
    }

    expect_equal(
      current_rows$cohort,
      expected_rows$cohort,
      info = paste("cohort drift for", estimator)
    )
    expect_equal(
      current_rows$period,
      expected_rows$period,
      info = paste("period drift for", estimator)
    )
    expect_equal(
      current_rows$event_time,
      expected_rows$event_time,
      info = paste("event-time drift for", estimator)
    )
    expect_numeric_or_na_equal(
      current_rows$att,
      expected_rows$att,
      tolerance = 1e-5,
      info = paste("ATT drift for", estimator)
    )
    expect_numeric_or_na_equal(
      current_rows$se,
      expected_rows$se,
      tolerance = 1e-5,
      info = paste("SE drift for", estimator)
    )
    expect_numeric_or_na_equal(
      current_rows$pvalue,
      expected_rows$pvalue,
      tolerance = 1e-5,
      info = paste("p-value drift for", estimator)
    )
  }
})

test_that("E6-06 Layer 3 castle keeps the current IPW group-time drift visible", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_parity_fixture_path(
    "e6_06_layer3_castle_realdata_comparator.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing E6-06 Layer 3 real-data comparator:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_e606_realdata_oracle("e6_06_layer3_castle_realdata_comparator.json")
  current <- current_e606_castle_result(estimator = "ipw", aggregate = "none")
  current_rows <- normalize_e606_group_time_rows(current$result$att_by_cohort_time)

  expect_identical(oracle$estimators$ipw$group_time_status, "row-count-drift")
  expect_equal(
    nrow(current_rows),
    oracle$estimators$ipw$group_time_comparison$r_rows
  )
  expect_equal(
    length(oracle$estimators$ipw$shared_group_time_rows),
    oracle$estimators$ipw$group_time_comparison$python_rows
  )
  expect_gt(nrow(current_rows), 0L)
})

test_that("E6-06 Layer 3 castle keeps the current overall aggregation drift visible", {
  skip_if_not_installed("jsonlite")

  oracle_path <- resolve_parity_fixture_path(
    "e6_06_layer3_castle_realdata_comparator.json"
  )

  expect_true(
    file.exists(oracle_path),
    info = paste("missing E6-06 Layer 3 real-data comparator:", oracle_path)
  )

  if (!file.exists(oracle_path)) {
    return(invisible(NULL))
  }

  oracle <- read_e606_realdata_oracle("e6_06_layer3_castle_realdata_comparator.json")

  for (estimator in c("ipw", "ipwra", "psm")) {
    current <- current_e606_castle_result(estimator = estimator, aggregate = "overall")
    current_overall <- list(
      att = as.numeric(current$result$att_overall),
      se = as.numeric(current$result$se_overall),
      ci_lower = as.numeric(current$result$ci_overall_lower),
      ci_upper = as.numeric(current$result$ci_overall_upper),
      pvalue = as.numeric(current$result$pvalue_overall)
    )
    expected_r <- oracle$estimators[[estimator]]$r_overall
    expected_python <- oracle$estimators[[estimator]]$python_overall

    expect_identical(
      oracle$estimators[[estimator]]$overall_status,
      "drift-confirmed"
    )
    expect_equal(current_overall$att, expected_r$att, tolerance = 1e-10)
    expect_equal(current_overall$se, expected_r$se, tolerance = 1e-10)
    expect_equal(current_overall$ci_lower, expected_r$ci_lower, tolerance = 1e-10)
    expect_equal(current_overall$ci_upper, expected_r$ci_upper, tolerance = 1e-10)
    expect_equal(current_overall$pvalue, expected_r$pvalue, tolerance = 1e-10)
    expect_true(
      abs(current_overall$att - expected_python$att) >
        oracle$estimators[[estimator]]$overall_att_tolerance,
      info = paste("overall ATT drift unexpectedly closed for", estimator)
    )
  }
})
