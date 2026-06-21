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

test_that("E6-06 Layer 3 castle IPW accepts propensity-only controls", {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)
  ps_vars <- c("income", "unemployrt", "poverty")

  control_alias <- suppressWarnings(lwdid(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipw",
    aggregate = "none",
    controls = ps_vars
  ))
  ps_explicit <- suppressWarnings(lwdid(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipw",
    aggregate = "none",
    controls = NULL,
    ps_controls = ps_vars
  ))

  alias_rows <- normalize_e606_group_time_rows(control_alias$att_by_cohort_time)
  explicit_rows <- normalize_e606_group_time_rows(ps_explicit$att_by_cohort_time)

  expect_equal(nrow(explicit_rows), 20L)
  expect_equal(explicit_rows$cohort, alias_rows$cohort)
  expect_equal(explicit_rows$period, alias_rows$period)
  expect_equal(explicit_rows$event_time, alias_rows$event_time)
  expect_equal(explicit_rows$att, alias_rows$att, tolerance = 1e-12)
  expect_equal(explicit_rows$se, alias_rows$se, tolerance = 1e-12)
})

test_that("E6-06 Layer 3 castle IPW and PSM use explicit ps_controls only", {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)
  ps_vars <- c("income", "unemployrt")

  panel_with_unused_na <- panel
  panel_with_unused_na$poverty <- NA_real_

  for (estimator in c("ipw", "psm")) {
    base_args <- list(
      data = panel,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = estimator,
      aggregate = "none",
      controls = NULL,
      ps_controls = ps_vars
    )
    with_unused_args <- base_args
    with_unused_args$data <- panel_with_unused_na
    with_unused_args$controls <- "poverty"
    if (identical(estimator, "psm")) {
      base_args$caliper <- 0.2
      with_unused_args$caliper <- 0.2
    }

    base <- suppressWarnings(do.call(lwdid, base_args))
    with_unused <- suppressWarnings(do.call(lwdid, with_unused_args))

    base_rows <- normalize_e606_group_time_rows(base$att_by_cohort_time)
    unused_rows <- normalize_e606_group_time_rows(with_unused$att_by_cohort_time)

    expect_equal(nrow(unused_rows), nrow(base_rows), info = estimator)
    expect_equal(unused_rows$cohort, base_rows$cohort, info = estimator)
    expect_equal(unused_rows$period, base_rows$period, info = estimator)
    expect_equal(unused_rows$att, base_rows$att, tolerance = 1e-12, info = estimator)
    expect_equal(unused_rows$se, base_rows$se, tolerance = 1e-12, info = estimator)
  }
})

test_that("E6-06 Layer 3 castle pretreatment uses explicit ps_controls", {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)
  ps_vars <- c("income", "unemployrt", "poverty")

  result <- suppressWarnings(lwdid(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipw",
    controls = NULL,
    ps_controls = ps_vars,
    aggregate = "cohort",
    control_group = "never_treated",
    include_pretreatment = TRUE,
    pretreatment_test = FALSE
  ))

  pre <- result$att_pre_treatment
  expect_s3_class(result, "lwdid_result")
  expect_true(is.data.frame(pre))
  expect_gt(sum(!pre$is_anchor & is.finite(pre$att)), 0L)
  warning_text <- paste(unlist(result$warnings_log), collapse = "\n")
  expect_false(
    grepl("requires 'controls' or 'ps_controls'", warning_text, fixed = TRUE),
    info = "pre-treatment dispatch should receive explicit ps_controls"
  )
})

test_that("E6-06 Layer 3 castle RI reuses non-RA ps_controls", {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)
  panel$unused_all_na <- NA_real_

  result <- suppressWarnings(lwdid(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipw",
    controls = "unused_all_na",
    ps_controls = c("income", "unemployrt", "poverty"),
    aggregate = "none",
    control_group = "never_treated",
    ri = TRUE,
    ri_method = "permutation",
    rireps = 50L,
    seed = 99L
  ))

  expect_s3_class(result, "lwdid_result")
  expect_equal(result$ri_valid, 50L)
  expect_equal(result$ri_failed, 0L)
  expect_null(result$ri_error)
  expect_true(is.finite(result$ri_pvalue))
  expect_equal(result$ri_target, "cohort_time")
  expect_equal(result$ri_estimator, "ipw")

  valid_rows <- result$att_by_cohort_time[
    is.finite(result$att_by_cohort_time$att), ,
    drop = FALSE
  ]
  expect_gt(nrow(valid_rows), 0L)
  expect_equal(result$ri_observed_stat, valid_rows$att[1L], tolerance = 1e-12)
  expect_false(
    isTRUE(all.equal(result$ri_observed_stat, result$att, tolerance = 1e-12)),
    info = "staggered aggregate='none' RI must report the cohort-time target, not the scalar ATT summary"
  )
  expect_equal(
    result$ri_pvalue,
    sum(abs(result$ri_distribution) >= abs(result$ri_observed_stat)) /
      length(result$ri_distribution),
    tolerance = 0
  )

  s <- summary(result)
  expect_equal(s$ri_observed_stat, result$ri_observed_stat, tolerance = 1e-12)
  expect_equal(s$ri_estimator, "ipw")

  d <- to_dict(result)
  expect_equal(d$ri_observed_stat, result$ri_observed_stat, tolerance = 1e-12)
  expect_equal(d$ri_estimator, "ipw")
})

test_that("E6-06 Layer 3 castle IPWRA accepts distinct outcome and ps controls", {
  data_path <- "/Users/cxy/Desktop/lwdid_r/lwdid-py_v0.2.3/data/castle.csv"
  panel <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  panel$gvar <- ifelse(is.na(panel$effyear), 0, panel$effyear)

  result <- suppressWarnings(lwdid(
    data = panel,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ipwra",
    aggregate = "none",
    controls = "income",
    ps_controls = c("unemployrt", "poverty")
  ))

  rows <- normalize_e606_group_time_rows(result$att_by_cohort_time)
  expect_gt(nrow(rows), 0L)
  expect_true(all(is.finite(rows$att)))
  expect_true(all(is.finite(rows$se)))
})

test_that("E6-06 Layer 3 castle non-RA overall follows cohort-effect aggregation", {
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
    cohort_df <- current$result$att_by_cohort[
      is.finite(current$result$att_by_cohort$att) &
        is.finite(current$result$att_by_cohort$se), ,
      drop = FALSE
    ]
    weights <- unlist(current$result$cohort_weights)
    cohort_weights <- as.numeric(weights[as.character(cohort_df$cohort)])
    expected_att <- sum(cohort_weights * cohort_df$att)
    expected_se <- sqrt(sum(cohort_weights^2 * cohort_df$se^2))

    expect_identical(current$result$estimator, estimator)
    expect_equal(
      current$result$overall_effect$aggregation_method,
      "cohort_size_weighted_non_ra_cohort_effects"
    )
    expect_equal(current$result$att_overall, expected_att, tolerance = 1e-10)
    expect_equal(current$result$se_overall, expected_se, tolerance = 1e-10)
    expect_true(
      abs(current$result$att_overall - oracle$estimators[[estimator]]$python_overall$att) >
        oracle$estimators[[estimator]]$overall_att_tolerance,
      info = paste("overall ATT drift unexpectedly closed against the frozen Python oracle for", estimator)
    )
  }
})
