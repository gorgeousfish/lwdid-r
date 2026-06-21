test_that("Castle event-time aggregation renormalizes cohort-size weights", {
  data(castle, package = "lwdid")

  fit <- suppressWarnings(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ra",
      aggregate = "event_time",
      control_group = "never_treated",
      verbose = "quiet"
    )
  )

  expect_equal(fit$aggregate, "event_time")
  expect_false(is.null(fit$event_time_effects))

  cohort_sizes <- unlist(fit$cohort_sizes, use.names = TRUE)
  expected_sizes <- c(
    "2005" = 1,
    "2006" = 13,
    "2007" = 4,
    "2008" = 2,
    "2009" = 1
  )
  expect_equal(cohort_sizes[names(expected_sizes)], expected_sizes)

  event_times <- vapply(fit$event_time_effects, `[[`, integer(1), "event_time")
  expect_equal(event_times, 0:5)

  for (effect in fit$event_time_effects) {
    contributions <- effect$cohort_contributions
    cohorts <- vapply(contributions, `[[`, integer(1), "cohort")
    weights <- vapply(contributions, `[[`, numeric(1), "weight")
    available_sizes <- expected_sizes[as.character(cohorts)]
    expected_weights <- as.numeric(available_sizes / sum(available_sizes))

    expect_equal(effect$n_cohorts, length(contributions))
    expect_equal(effect$weight_sum, 1, tolerance = 1e-12)
    expect_equal(weights, expected_weights, tolerance = 1e-12)
    expect_equal(effect$se_aggregation, "diagonal_weighted_cohort_se")
    expect_equal(effect$covariance_assumption, "zero_cross_cohort_covariance")
    expect_equal(
      effect$se,
      sqrt(sum(weights^2 * vapply(contributions, `[[`, numeric(1), "se")^2)),
      tolerance = 1e-12
    )
  }

  last_effect <- fit$event_time_effects[[length(fit$event_time_effects)]]
  expect_equal(last_effect$event_time, 5L)
  expect_equal(last_effect$n_cohorts, 1L)
  expect_equal(last_effect$cohort_contributions[[1L]]$cohort, 2005L)
  expect_equal(last_effect$cohort_contributions[[1L]]$weight, 1)

  contribution_df <- extract_effects(fit, type = "event_time_contributions")
  expect_equal(nrow(contribution_df), 20L)
  expect_true("se_aggregation" %in% names(contribution_df))
  expect_true("covariance_assumption" %in% names(contribution_df))
  expect_equal(
    unique(contribution_df$se_aggregation),
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    unique(contribution_df$covariance_assumption),
    "zero_cross_cohort_covariance"
  )
  expect_equal(
    sort(unique(contribution_df$event_time)),
    0:5
  )
  for (event_time in unique(contribution_df$event_time)) {
    event_rows <- contribution_df[contribution_df$event_time == event_time, ]
    expect_equal(sum(event_rows$weight), 1, tolerance = 1e-12)
    expect_equal(
      unique(event_rows$n_cohorts),
      nrow(event_rows)
    )
  }
  expect_equal(
    contribution_df$contribution_att[contribution_df$event_time == 5L],
    contribution_df$event_att[contribution_df$event_time == 5L],
    tolerance = 1e-12
  )

  event_time_df <- extract_effects(fit, type = "event_time")
  expect_true("se_aggregation" %in% names(event_time_df))
  expect_true("covariance_assumption" %in% names(event_time_df))
  expect_equal(
    unique(event_time_df$se_aggregation),
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    unique(event_time_df$covariance_assumption),
    "zero_cross_cohort_covariance"
  )

  event_time_coef <- coef(fit, type = "event_time")
  event_time_ci <- confint(fit, type = "event_time")
  event_time_vcov <- vcov(fit, type = "event_time")
  event_time_names <- paste0("e", event_time_df$event_time)

  expect_equal(names(event_time_coef), event_time_names)
  expect_equal(unname(event_time_coef), event_time_df$att, tolerance = 1e-12)
  expect_equal(rownames(event_time_ci), event_time_names)
  expect_equal(rownames(event_time_vcov), event_time_names)
  expect_equal(colnames(event_time_vcov), event_time_names)
  expect_equal(
    unname(diag(event_time_vcov)),
    event_time_df$se^2,
    tolerance = 1e-12
  )
  expect_equal(
    attr(event_time_vcov, "se_aggregation"),
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    attr(event_time_vcov, "covariance_assumption"),
    "zero_cross_cohort_covariance"
  )

  contribution_file <- tempfile(fileext = ".csv")
  on.exit(unlink(contribution_file), add = TRUE)
  exported_contributions <- suppressMessages(
    to_csv(
      fit,
      file = contribution_file,
      what = "event_time_contributions"
    )
  )
  expect_equal(exported_contributions, contribution_df)
  expect_equal(nrow(utils::read.csv(contribution_file)), 20L)

  dict <- to_dict(fit)
  expect_true("event_time_contributions" %in% names(dict))
  expect_equal(dict$event_time_contributions, contribution_df)

  default_tex <- to_latex(fit)
  expect_false(grepl("Event time & WATT\\(e\\)", default_tex))

  event_time_tex <- to_latex(fit, include_event_time = TRUE)
  expect_true(grepl("Event time & WATT\\(e\\)", event_time_tex))
  expect_true(grepl(
    "p-value & CI & Cohorts & Min treated & Max weight CV",
    event_time_tex,
    fixed = TRUE
  ))
  expect_true(grepl("0 &", event_time_tex, fixed = TRUE))
  expect_true(grepl("5 &", event_time_tex, fixed = TRUE))
  expect_true(grepl("diagonal cohort-SE aggregation", event_time_tex, fixed = TRUE))
  expect_true(grepl("Cross-cohort covariance is not modeled", event_time_tex, fixed = TRUE))
})

test_that("Castle IPW event-time inference uses the normal scale consistently", {
  data(castle, package = "lwdid")
  controls <- c("income", "unemployrt", "poverty")

  fit <- suppressWarnings(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ipw",
      aggregate = "event_time",
      control_group = "never_treated",
      controls = controls,
      return_diagnostics = TRUE,
      verbose = "quiet"
    )
  )

  event_time_df <- extract_effects(fit, type = "event_time")
  expect_equal(unique(event_time_df$inference_dist), "normal")
  expect_true(all(is.na(event_time_df$df_inference)))

  z_crit <- stats::qnorm(1 - fit$alpha / 2)
  expected_p <- 2 * stats::pnorm(
    abs(event_time_df$att / event_time_df$se),
    lower.tail = FALSE
  )
  expect_equal(event_time_df$pvalue, expected_p, tolerance = 1e-12)
  expect_equal(
    event_time_df$ci_lower,
    event_time_df$att - z_crit * event_time_df$se,
    tolerance = 1e-12
  )
  expect_equal(
    event_time_df$ci_upper,
    event_time_df$att + z_crit * event_time_df$se,
    tolerance = 1e-12
  )

  plot_df <- plot_event_study(fit, return_data = TRUE)$data
  plot_post <- plot_df[
    !plot_df$is_anchor & plot_df$event_time >= 0L,
    , drop = FALSE
  ]
  plot_post <- plot_post[order(plot_post$event_time), , drop = FALSE]
  event_time_df <- event_time_df[order(event_time_df$event_time), , drop = FALSE]

  expect_equal(plot_post$event_time, event_time_df$event_time)
  expect_equal(plot_post$att, event_time_df$att, tolerance = 1e-12)
  expect_equal(plot_post$se, event_time_df$se, tolerance = 1e-12)
  expect_equal(plot_post$ci_lower, event_time_df$ci_lower, tolerance = 1e-12)
  expect_equal(plot_post$ci_upper, event_time_df$ci_upper, tolerance = 1e-12)

  contribution_df <- extract_effects(fit, type = "event_time_contributions")
  for (event_time in unique(contribution_df$event_time)) {
    event_rows <- contribution_df[
      contribution_df$event_time == event_time,
      , drop = FALSE
    ]
    event_row <- event_time_df[
      event_time_df$event_time == event_time,
      , drop = FALSE
    ]
    expect_equal(sum(event_rows$weight), 1, tolerance = 1e-12)
    expect_equal(
      event_row$att,
      sum(event_rows$weight * event_rows$contribution_att),
      tolerance = 1e-12
    )
    expect_equal(
      event_row$se,
      sqrt(sum(event_rows$weight^2 * event_rows$contribution_se^2)),
      tolerance = 1e-12
    )
    expect_equal(event_row$min_n_treated, min(event_rows$n_treated))
    expect_equal(event_row$max_n_treated, max(event_rows$n_treated))
    expect_equal(event_row$min_n_control, min(event_rows$n_control))
    expect_equal(event_row$max_n_control, max(event_rows$n_control))
  }
  expect_equal(min(event_time_df$min_n_treated), 1L)
  expect_true(all(event_time_df$min_n_control > 0L))

  event_time_tex <- to_latex(fit, include_event_time = TRUE)
  expect_true(grepl("Event time & WATT\\(e\\)", event_time_tex))
  expect_true(grepl(
    "p-value & CI & Cohorts & Min treated & Max weight CV",
    event_time_tex,
    fixed = TRUE
  ))
  expect_true(grepl("Min treated and max weight CV summarize contributing cells", event_time_tex, fixed = TRUE))

  contribution_file <- tempfile(fileext = ".csv")
  on.exit(unlink(contribution_file), add = TRUE)
  exported_contributions <- suppressMessages(
    to_csv(
      fit,
      file = contribution_file,
      what = "event_time_contributions"
    )
  )
  expect_equal(exported_contributions, contribution_df)
})

test_that("event-time aggregation preserves bootstrap replicate diagnostics", {
  cohort_time_effects <- data.frame(
    cohort = c(2005L, 2007L, 2005L, 2007L),
    period = c(2005L, 2007L, 2006L, 2008L),
    att = c(0.4, 0.8, 0.6, 1.0),
    se = c(0.20, 0.25, 0.30, 0.35),
    df_inference = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
    n_treated = c(8L, 12L, 8L, 12L),
    n_control = c(30L, 28L, 30L, 28L),
    bootstrap_reps_requested = c(40L, 40L, 40L, 40L),
    bootstrap_reps_valid = c(38L, 35L, 36L, 34L),
    bootstrap_reps_failed = c(2L, 5L, 4L, 6L),
    bootstrap_success_rate = c(0.95, 0.875, 0.90, 0.85)
  )

  effects <- aggregate_to_event_time(
    cohort_time_effects = cohort_time_effects,
    cohort_sizes = c("2005" = 80L, "2007" = 120L),
    inference_dist = "normal"
  )
  result <- new_lwdid_result(
    event_time_effects = effects,
    is_staggered = TRUE,
    aggregate = "event_time",
    estimator = "ipw"
  )

  event_rows <- extract_effects(result, type = "event_time")
  contribution_rows <- extract_effects(result, type = "event_time_contributions")

  expect_true(all(c(
    "min_bootstrap_success_rate", "min_bootstrap_reps_valid",
    "max_bootstrap_reps_failed"
  ) %in% names(event_rows)))
  expect_equal(event_rows$min_bootstrap_success_rate, c(0.875, 0.85))
  expect_equal(event_rows$min_bootstrap_reps_valid, c(35L, 34L))
  expect_equal(event_rows$max_bootstrap_reps_failed, c(5L, 6L))
  expect_true(all(c(
    "bootstrap_reps_requested", "bootstrap_reps_valid",
    "bootstrap_reps_failed", "bootstrap_success_rate"
  ) %in% names(contribution_rows)))
  expect_equal(
    contribution_rows$bootstrap_success_rate,
    contribution_rows$bootstrap_reps_valid /
      contribution_rows$bootstrap_reps_requested,
    tolerance = 1e-12
  )
})

test_that("Castle PSM event-time output omits non-estimable WATT rows", {
  data(castle, package = "lwdid")
  controls <- c("income", "unemployrt", "poverty")

  fit <- suppressWarnings(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "psm",
      aggregate = "event_time",
      control_group = "never_treated",
      controls = controls,
      seed = 123,
      verbose = "quiet"
    )
  )

  event_time_df <- extract_effects(fit, type = "event_time")
  expect_true(all(is.finite(event_time_df$att)))
  expect_true(all(is.finite(event_time_df$se)))
  expect_false(5L %in% event_time_df$event_time)

  contribution_df <- extract_effects(fit, type = "event_time_contributions")
  expect_equal(
    sort(unique(contribution_df$event_time)),
    sort(event_time_df$event_time)
  )

  plot_df <- plot_event_study(fit, return_data = TRUE)$data
  plot_post <- plot_df[
    !plot_df$is_anchor & plot_df$event_time >= 0L,
    , drop = FALSE
  ]
  expect_equal(
    sort(plot_post$event_time),
    sort(event_time_df$event_time)
  )
  expect_false(any(!is.finite(plot_post$att)))
  expect_false(any(!is.finite(plot_post$se)))
})

test_that("AET event-time print shows skipped diagnostics and WATT boundary", {
  set.seed(77709)
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 10), rep(5L, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 5.0 * treated + rnorm(.N, 0, 0.01)]
  dt[, treated := NULL]

  fit <- suppressWarnings(lwdid(
    data = dt,
    y = "Y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "demean",
    control_group = "not_yet_treated",
    aggregate = "event_time",
    verbose = "quiet"
  ))

  fit_summary <- summary(fit)
  event_time_df <- extract_effects(fit, type = "event_time")
  printed <- paste(capture.output(print(fit)), collapse = "\n")

  expect_true(length(fit$skipped_pairs) > 0L)
  expect_equal(fit_summary$n_skipped_pairs, length(fit$skipped_pairs))
  expect_identical(fit_summary$skipped_summary, fit$skipped_summary)
  expect_true("no_control" %in% fit_summary$skipped_summary$reason)
  expect_true(nrow(event_time_df) > 0L)
  expect_equal(
    unique(event_time_df$se_aggregation),
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    unique(event_time_df$covariance_assumption),
    "zero_cross_cohort_covariance"
  )
  expect_true(grepl("Skipped \\(g,r\\):", printed))
  expect_true(grepl("no_control:", printed, fixed = TRUE))
  expect_true(grepl("event-time WATT\\(e\\) rows are below", printed))
  expect_true(grepl("cross-cohort covariance is not modeled", printed))
  expect_false(grepl("no finite event-time WATT\\(e\\) rows", printed))

  dict <- to_dict(fit)
  expect_equal(dict$n_skipped_pairs, length(fit$skipped_pairs))
  expect_identical(dict$skipped_summary, fit$skipped_summary)

  summary_file <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_file), add = TRUE)
  summary_csv <- suppressMessages(
    to_csv(fit, file = summary_file, what = "summary")
  )
  expect_equal(summary_csv$n_skipped_pairs, length(fit$skipped_pairs))
  expect_true(grepl("no_control:", summary_csv$skipped_summary, fixed = TRUE))
})

test_that("AET event-time omitted WATT rows remain inspectable", {
  units <- 1:20
  periods <- 1:7
  gvar_map <- c(rep(3L, 10), rep(5L, 10))
  dt <- data.table::CJ(id = units, time = periods)
  dt[, gvar := gvar_map[id]]
  dt[, treated := data.table::fifelse(time >= gvar, 1L, 0L)]
  dt[, Y := as.numeric(id) + 5.0 * treated]
  dt[, treated := NULL]

  fit <- suppressWarnings(lwdid(
    data = dt,
    y = "Y",
    ivar = "id",
    tvar = "time",
    gvar = "gvar",
    rolling = "demean",
    control_group = "not_yet_treated",
    aggregate = "event_time",
    verbose = "quiet"
  ))

  expect_equal(length(fit$event_time_effects), 0L)
  expect_equal(length(fit$event_time_omissions), 2L)
  expect_s3_class(fit$event_time_omission_summary, "data.frame")
  expect_equal(
    fit$event_time_omission_summary,
    data.frame(
      reason = "non_finite_att_or_se",
      n = 2L,
      stringsAsFactors = FALSE
    )
  )
  expect_equal(
    vapply(fit$event_time_omissions, `[[`, integer(1L), "event_time"),
    0:1
  )
  expect_true(all(vapply(
    fit$event_time_omissions,
    function(x) identical(x$reason, "non_finite_att_or_se"),
    logical(1L)
  )))

  expect_silent(fit_summary <- summary(fit))
  expect_equal(fit_summary$n_omitted_event_times, 2L)
  expect_identical(
    fit_summary$event_time_omission_summary,
    fit$event_time_omission_summary
  )

  expect_silent(printed <- paste(capture.output(print(fit)), collapse = "\n"))
  expect_true(grepl("Omitted WATT\\(e\\):[[:space:]]+2", printed))
  expect_true(grepl("non_finite_att_or_se: 2", printed, fixed = TRUE))
  expect_true(grepl("no finite event-time WATT\\(e\\) rows are available", printed))
  expect_false(grepl("event-time WATT\\(e\\) rows are below", printed))

  expect_silent(dict <- to_dict(fit))
  expect_equal(dict$n_omitted_event_times, 2L)
  expect_identical(dict$event_time_omissions, fit$event_time_omissions)
  expect_identical(
    dict$event_time_omission_summary,
    fit$event_time_omission_summary
  )

  summary_file <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_file), add = TRUE)
  summary_csv <- suppressMessages(
    to_csv(fit, file = summary_file, what = "summary")
  )
  expect_equal(summary_csv$n_omitted_event_times, 2L)
  expect_equal(
    summary_csv$event_time_omission_summary,
    "non_finite_att_or_se: 2"
  )
  expect_equal(summary_csv$n_skipped_pairs, length(fit$skipped_pairs))
  expect_true(grepl("no_control:", summary_csv$skipped_summary, fixed = TRUE))
})

test_that("Castle event-time range is respected in aggregated and faceted plots", {
  data(castle, package = "lwdid")
  controls <- c("income", "unemployrt", "poverty")

  fit <- suppressWarnings(
    lwdid(
      data = castle,
      y = "lhomicide",
      ivar = "sid",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      estimator = "ipw",
      aggregate = "event_time",
      control_group = "never_treated",
      controls = controls,
      include_pretreatment = TRUE,
      event_time_range = c(0, 2),
      verbose = "quiet"
    )
  )

  event_time_df <- extract_effects(fit, type = "event_time")
  expect_equal(event_time_df$event_time, 0:2)

  weighted_plot <- plot_event_study(fit, return_data = TRUE)$data
  weighted_post <- weighted_plot[
    !weighted_plot$is_anchor & weighted_plot$event_time >= 0L,
    , drop = FALSE
  ]
  expect_equal(sort(unique(weighted_post$event_time)), 0:2)
  expect_true(any(weighted_plot$event_time < 0L))

  mean_plot <- suppressWarnings(
    plot_event_study(fit, aggregation = "mean", return_data = TRUE)$data
  )
  mean_post <- mean_plot[
    !mean_plot$is_anchor & mean_plot$event_time >= 0L,
    , drop = FALSE
  ]
  expect_equal(sort(unique(mean_post$event_time)), 0:2)
  expect_true(any(mean_plot$event_time < 0L))

  faceted_plot <- plot_event_study(
    fit,
    facet_by_cohort = TRUE,
    return_data = TRUE
  )$data
  faceted_post <- faceted_plot[
    !faceted_plot$is_anchor & faceted_plot$event_time >= 0L,
    , drop = FALSE
  ]
  expect_equal(sort(unique(faceted_post$event_time)), 0:2)
  expect_true(any(faceted_plot$event_time < 0L))
})
