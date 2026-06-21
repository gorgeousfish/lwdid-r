# ============================================================================
# test-event-study-doc-contract.R
# Event-study documentation examples must match the public plotting contract.
# ============================================================================

test_that("event-study docs expose stricter anchor and ref_period examples", {
  doc_paths <- c(
    resolve_package_source_file("README.md"),
    resolve_package_source_file("vignettes", "getting-started.Rmd"),
    resolve_package_source_file("inst", "doc", "getting-started.Rmd")
  )

  for (path in doc_paths) {
    expect_true(file.exists(path), info = sprintf("Missing doc file: %s", path))
    text <- paste(readLines(path, warn = FALSE), collapse = "\n")

    expect_match(
      text,
      "plot\\(result_(event|es), show_pre_treatment = FALSE\\)",
      info = sprintf("%s must document hiding pre-treatment periods", path)
    )
    expect_match(
      text,
      "plot\\(result_(event|es), ref_period = 0\\)",
      info = sprintf("%s must document observed reference-period normalization", path)
    )
    expect_match(
      text,
      "does not relabel\\s+(that observed row|the observed\\s+reference row) as a visual anchor",
      info = sprintf("%s must document that ref_period is not a visual anchor", path)
    )
    expect_match(
      text,
      "printed scalar ATT summarizes valid\\s+cohort-period `\\(g, r\\)` effects",
      info = sprintf("%s must distinguish scalar ATT from WATT(e) rows", path)
    )
    expect_match(
      text,
      "event-study table and plot report WATT\\(e\\)\\s+rows by relative time",
      info = sprintf("%s must document event-time WATT(e) rows", path)
    )
    expect_match(
      text,
      "Plot data returned by `plot_event_study\\(result_(event|es), return_data = TRUE\\)`\\s+carries the same standard-error metadata",
      info = sprintf("%s must document event-study plot-data SE metadata", path)
    )
    expect_match(
      text,
      "plot\\(result_(event|es)\\)\\s*```\\s*For `aggregate = \"event_time\"`",
      info = sprintf("%s must close the event-study code block before prose", path)
    )
    expect_match(
      text,
      "`se_aggregation` and `covariance_assumption`",
      fixed = TRUE,
      info = sprintf("%s must name the event-study SE metadata columns", path)
    )
    expect_match(
      text,
      "extract_effects\\(result_(event|es), type = \"event_time_contributions\"\\)",
      info = sprintf("%s must document event-time contribution extraction", path)
    )
    expect_match(
      text,
      "what = \"event_time_contributions\"",
      info = sprintf("%s must document event-time contribution CSV export", path)
    )
    expect_match(
      text,
      "to_dict\\(result_(event|es)\\)\\$event_time_contributions",
      info = sprintf("%s must document event-time contribution dictionary export", path)
    )
  }
})

test_that("event-study docs examples run on Castle and preserve plot semantics", {
  skip_if_not_installed("ggplot2")

  data("castle", package = "lwdid")

  result_event <- lwdid(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = "ra",
    aggregate = "event_time"
  )

  p_default <- suppressWarnings(plot(result_event))
  p_hidden <- suppressWarnings(plot(result_event, show_pre_treatment = FALSE))
  p_ref0 <- suppressWarnings(plot(result_event, ref_period = 0))

  expect_s3_class(p_default, "ggplot")
  expect_s3_class(p_hidden, "ggplot")
  expect_s3_class(p_ref0, "ggplot")

  hidden_data <- suppressWarnings(
    plot_event_study(
      result_event,
      show_pre_treatment = FALSE,
      return_data = TRUE
    )$data
  )
  ref0_data <- suppressWarnings(
    plot_event_study(result_event, ref_period = 0, return_data = TRUE)$data
  )

  expect_true(all(hidden_data$event_time >= 0L))
  expect_false(any(hidden_data$is_anchor))
  ref0_row <- ref0_data[ref0_data$event_time == 0L, , drop = FALSE]
  expect_equal(nrow(ref0_row), 1L)
  expect_equal(ref0_row$att, 0, tolerance = 1e-12)
  expect_false(ref0_row$is_anchor)
  expect_true("se_aggregation" %in% names(hidden_data))
  expect_true("covariance_assumption" %in% names(hidden_data))
  expect_true(all(
    hidden_data$se_aggregation[!hidden_data$is_anchor] ==
      "diagonal_weighted_cohort_se"
  ))
  expect_true(all(
    hidden_data$covariance_assumption[!hidden_data$is_anchor] ==
      "zero_cross_cohort_covariance"
  ))
})

test_that("plot_event_study Rd default aggregation matches implementation", {
  rd_path <- resolve_package_source_file(
    "man", "plot_event_study.lwdid_result.Rd"
  )
  text <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")

  expect_match(
    text,
    'aggregation = c\\("weighted", "mean"\\)',
    info = "Rd usage must preserve the weighted default used by the implementation"
  )
})
