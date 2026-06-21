# ============================================================================
# test-staggered-control-group-realdata.R
# Real-data control group contracts for non-RA staggered estimation.
# ============================================================================

expect_realdata_staggered_control_counts <- function(data, y, ivar, tvar, gvar,
                                                     controls, dataset_name) {
  dt <- data.table::as.data.table(data)

  fit_nyt <- suppressWarnings(lwdid(
    data = data,
    y = y,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    rolling = "demean",
    estimator = "ipw",
    aggregate = "none",
    control_group = "not_yet_treated",
    controls = controls,
    verbose = "quiet"
  ))
  fit_nt <- suppressWarnings(lwdid(
    data = data,
    y = y,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    rolling = "demean",
    estimator = "ipw",
    aggregate = "none",
    control_group = "never_treated",
    controls = controls,
    verbose = "quiet"
  ))

  cte_nyt <- fit_nyt$att_by_cohort_time
  cte_nt <- fit_nt$att_by_cohort_time
  key_cols <- c("cohort", "period")

  expected_counts <- lapply(seq_len(nrow(cte_nyt)), function(i) {
    g <- cte_nyt$cohort[i]
    r <- cte_nyt$period[i]
    period_data <- dt[get(tvar) == r]
    data.frame(
      cohort = g,
      period = r,
      n_treated_expected = sum(
        !is_never_treated(period_data[[gvar]]) & period_data[[gvar]] == g
      ),
      n_nyt_expected = sum(
        get_valid_controls(period_data, gvar, g, r, "not_yet_treated")
      ),
      n_nt_expected = sum(
        get_valid_controls(period_data, gvar, g, r, "never_treated")
      )
    )
  })
  expected_counts <- do.call(rbind, expected_counts)

  merged_nyt <- merge(
    cte_nyt[, c(key_cols, "n_treated", "n_control", "vce_type")],
    expected_counts,
    by = key_cols
  )
  merged_nt <- merge(
    cte_nt[, c(key_cols, "n_treated", "n_control", "vce_type")],
    expected_counts,
    by = key_cols
  )

  expect_equal(
    merged_nyt$n_treated,
    merged_nyt$n_treated_expected,
    info = sprintf("%s not_yet_treated treated counts", dataset_name)
  )
  expect_equal(
    merged_nyt$n_control,
    merged_nyt$n_nyt_expected,
    info = sprintf("%s not_yet_treated control counts", dataset_name)
  )
  expect_equal(
    merged_nt$n_control,
    merged_nt$n_nt_expected,
    info = sprintf("%s never_treated control counts", dataset_name)
  )
  expect_true(
    all(merged_nyt$n_control >= merged_nt$n_control),
    info = sprintf("%s NYT controls should include NT controls", dataset_name)
  )
  expect_true(
    any(merged_nyt$n_control > merged_nt$n_control),
    info = sprintf("%s should have real not-yet-treated extra controls", dataset_name)
  )
  expect_identical(unique(merged_nyt$vce_type), "analytical")
  expect_identical(unique(merged_nt$vce_type), "analytical")
  expect_identical(fit_nyt$vce_type, "analytical")
  expect_identical(fit_nt$vce_type, "analytical")
  expect_match(summary(fit_nyt)$vce_description, "Analytical")
  expect_match(summary(fit_nt)$vce_description, "Analytical")
}

test_that("Castle IPW control groups match direct period masks", {
  data("castle", package = "lwdid")

  expect_realdata_staggered_control_counts(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    dataset_name = "Castle"
  )
})

test_that("Walmart IPW control groups match direct period masks", {
  data("walmart", package = "lwdid")

  expect_realdata_staggered_control_counts(
    data = walmart,
    y = "log_retail_emp",
    ivar = "fips",
    tvar = "year",
    gvar = "g",
    controls = c("total_pop", "emp1964", "emp1977"),
    dataset_name = "Walmart"
  )
})

expect_staggered_propensity_diagnostics <- function(data, y, ivar, tvar, gvar,
                                                    controls, dataset_name) {
  fit <- suppressWarnings(lwdid(
    data = data,
    y = y,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    rolling = "demean",
    estimator = "ipw",
    aggregate = "none",
    control_group = "not_yet_treated",
    controls = controls,
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  propensity_diag <- get_diagnostics(fit, "propensity")
  expect_s3_class(fit, "lwdid_result")
  expect_s3_class(propensity_diag, "data.frame")
  expect_equal(
    nrow(propensity_diag),
    nrow(fit$att_by_cohort_time),
    info = sprintf("%s staggered propensity diagnostics row count", dataset_name)
  )
  expect_true(all(c("cohort", "period", "ps_min", "ps_max", "weights_cv") %in%
                    names(propensity_diag)))
  expect_true(
    all(is.finite(propensity_diag$ps_min)),
    info = sprintf("%s propensity minima finite", dataset_name)
  )
  expect_true(
    all(is.finite(propensity_diag$ps_max)),
    info = sprintf("%s propensity maxima finite", dataset_name)
  )
  expect_true(
    all(propensity_diag$ps_min >= 0 & propensity_diag$ps_max <= 1),
    info = sprintf("%s propensity ranges bounded", dataset_name)
  )
  expect_match(summary(fit)$diagnostics_summary$propensity, "cells=")
}

test_that("Castle IPW exposes staggered propensity diagnostics", {
  data("castle", package = "lwdid")

  expect_staggered_propensity_diagnostics(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    dataset_name = "Castle"
  )
})

test_that("Castle IPW retains real-data skipped-cell diagnostics", {
  data("castle", package = "lwdid")
  dt <- data.table::as.data.table(castle)
  cohorts <- sort(unique(dt$gvar[is.finite(dt$gvar)]))
  pre_stats <- precompute_transforms(
    dt = dt,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    cohorts = cohorts,
    rolling = "demean"
  )

  result <- suppressWarnings(estimate_staggered_effects(
    dt = dt,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    control_group = "not_yet_treated",
    controls = c("income", "unemployrt", "poverty"),
    vce = NULL,
    cluster_var = NULL,
    alpha = 0.05,
    pre_stats = pre_stats,
    estimator = "ipw",
    min_control = 30L,
    return_diagnostics = TRUE
  ))

  skipped_pairs <- attr(result, "skipped_pairs", exact = TRUE)
  skipped_summary <- attr(result, "skipped_summary", exact = TRUE)
  skipped_reasons <- vapply(
    skipped_pairs, function(x) x$reason, character(1)
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10L)
  expect_equal(nrow(result) + length(skipped_pairs), 20L)
  expect_equal(unique(skipped_summary$reason), "insufficient_control")
  expect_equal(skipped_summary$n, 10L)
  expect_true(all(skipped_reasons == "insufficient_control"))
  expect_true(all(result$n_control >= 30L))
  expect_true(all(vapply(skipped_pairs, function(x) {
    period_data <- dt[year == x$r]
    control_mask <- get_valid_controls(
      period_data, "gvar", x$g, x$r, "not_yet_treated"
    )
    sum(control_mask) < 30L
  }, logical(1))))
})

test_that("Walmart IPW exposes staggered propensity diagnostics", {
  data("walmart", package = "lwdid")

  expect_staggered_propensity_diagnostics(
    data = walmart,
    y = "log_retail_emp",
    ivar = "fips",
    tvar = "year",
    gvar = "g",
    controls = c("total_pop", "emp1964", "emp1977"),
    dataset_name = "Walmart"
  )
})

expect_realdata_event_time_export_metadata <- function(data, y, ivar, tvar, gvar,
                                                       controls, dataset_name) {
  fit <- suppressWarnings(lwdid(
    data = data,
    y = y,
    ivar = ivar,
    tvar = tvar,
    gvar = gvar,
    rolling = "demean",
    estimator = "ipw",
    aggregate = "event_time",
    control_group = "not_yet_treated",
    controls = controls,
    event_time_range = c(0, 2),
    return_diagnostics = TRUE,
    verbose = "quiet"
  ))

  event_time_df <- extract_effects(fit, type = "event_time")
  contribution_df <- extract_effects(fit, type = "event_time_contributions")
  plot_df <- plot_event_study(fit, return_data = TRUE)$data
  plot_post <- plot_df[
    !plot_df$is_anchor & plot_df$event_time >= 0L,
    , drop = FALSE
  ]
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  dict <- to_dict(fit)
  summary_file <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_file), add = TRUE)
  summary_csv <- suppressMessages(
    to_csv(fit, file = summary_file, what = "summary")
  )

  expect_true(
    nrow(event_time_df) > 0L,
    info = sprintf("%s event-time rows", dataset_name)
  )
  expect_true(
    all(is.finite(event_time_df$att)),
    info = sprintf("%s event-time ATT finite", dataset_name)
  )
  expect_true(
    all(is.finite(event_time_df$se)),
    info = sprintf("%s event-time SE finite", dataset_name)
  )
  expect_true(
    "event_time_contributions" %in% names(dict),
    info = sprintf("%s dictionary contribution export", dataset_name)
  )
  expect_true(
    all(c(
      "max_weight_cv", "weighted_weight_cv", "min_ps", "max_ps",
      "min_n_treated", "max_n_treated", "min_n_control", "max_n_control"
    ) %in%
          names(event_time_df)),
    info = sprintf("%s event-time support and overlap diagnostics", dataset_name)
  )
  expect_true(
    all(c("weights_cv", "ps_min", "ps_max", "n_treated", "n_control") %in%
          names(contribution_df)),
    info = sprintf("%s contribution overlap diagnostics", dataset_name)
  )
  expect_true(
    all(c("max_weight_cv", "weighted_weight_cv", "min_n_control") %in%
          names(plot_post)),
    info = sprintf("%s plot-data overlap diagnostics", dataset_name)
  )
  expect_true(
    all(c("weights_cv", "ps_min", "ps_max", "n_treated", "n_control") %in%
          names(dict$event_time_contributions)),
    info = sprintf("%s dictionary contribution overlap diagnostics", dataset_name)
  )
  expect_true(
    all(is.finite(event_time_df$max_weight_cv)),
    info = sprintf("%s event-time max weight CV finite", dataset_name)
  )
  expect_true(
    all(is.finite(contribution_df$weights_cv)),
    info = sprintf("%s contribution weight CV finite", dataset_name)
  )
  for (event_time in event_time_df$event_time) {
    event_row <- event_time_df[event_time_df$event_time == event_time, ]
    contribution_rows <- contribution_df[
      contribution_df$event_time == event_time,
      , drop = FALSE
    ]
    expect_equal(
      event_row$max_weight_cv,
      max(contribution_rows$weights_cv),
      tolerance = 1e-12,
      info = sprintf("%s event-time max CV e=%d", dataset_name, event_time)
    )
    expect_equal(
      event_row$weighted_weight_cv,
      sum(contribution_rows$weight * contribution_rows$weights_cv),
      tolerance = 1e-12,
      info = sprintf("%s event-time weighted CV e=%d", dataset_name, event_time)
    )
    plot_row <- plot_post[plot_post$event_time == event_time, , drop = FALSE]
    expect_equal(
      plot_row$max_weight_cv,
      event_row$max_weight_cv,
      tolerance = 1e-12,
      info = sprintf("%s plot-data max CV e=%d", dataset_name, event_time)
    )
    expect_equal(
      plot_row$weighted_weight_cv,
      event_row$weighted_weight_cv,
      tolerance = 1e-12,
      info = sprintf("%s plot-data weighted CV e=%d", dataset_name, event_time)
    )
    expect_equal(
      event_row$min_n_treated,
      min(contribution_rows$n_treated),
      info = sprintf("%s event-time min treated support e=%d", dataset_name, event_time)
    )
    expect_equal(
      event_row$max_n_treated,
      max(contribution_rows$n_treated),
      info = sprintf("%s event-time max treated support e=%d", dataset_name, event_time)
    )
    expect_equal(
      event_row$min_n_control,
      min(contribution_rows$n_control),
      info = sprintf("%s event-time min control support e=%d", dataset_name, event_time)
    )
    expect_equal(
      event_row$max_n_control,
      max(contribution_rows$n_control),
      info = sprintf("%s event-time max control support e=%d", dataset_name, event_time)
    )
  }
  if (identical(dataset_name, "Castle")) {
    expect_true(
      min(event_time_df$min_n_treated) == 1L,
      info = "Castle event-time rows expose sparse treated support"
    )
    expect_true(
      grepl("min treated cell N=1", printed, fixed = TRUE),
      info = "Castle print cue exposes sparse treated support"
    )
  }
  if (identical(dataset_name, "Walmart")) {
    expect_true(
      max(event_time_df$max_weight_cv) > 2,
      info = "Walmart event-time rows expose high IPW weight concentration"
    )
    expect_true(
      grepl("max weight CV=", printed, fixed = TRUE),
      info = "Walmart print cue exposes high IPW weight concentration"
    )
  }
  expect_equal(
    dict$n_skipped_pairs,
    length(fit$skipped_pairs),
    info = sprintf("%s dictionary skipped count", dataset_name)
  )
  expect_equal(
    summary_csv$n_skipped_pairs,
    length(fit$skipped_pairs),
    info = sprintf("%s summary CSV skipped count", dataset_name)
  )
  expect_equal(
    dict$n_omitted_event_times,
    length(fit$event_time_omissions),
    info = sprintf("%s dictionary omitted event-time count", dataset_name)
  )
  expect_equal(
    summary_csv$n_omitted_event_times,
    length(fit$event_time_omissions),
    info = sprintf("%s summary CSV omitted event-time count", dataset_name)
  )
  expect_equal(
    summary_csv$event_time_max_weight_cv,
    max(event_time_df$max_weight_cv),
    tolerance = 1e-12,
    info = sprintf("%s summary CSV max event-time weight CV", dataset_name)
  )
  expect_equal(
    summary_csv$event_time_min_n_treated,
    min(event_time_df$min_n_treated),
    info = sprintf("%s summary CSV minimum event-time treated support", dataset_name)
  )
}

test_that("Castle and Walmart IPW event-time exports preserve diagnostics", {
  data("castle", package = "lwdid")
  data("walmart", package = "lwdid")
  walmart_overlap_panel <- walmart[
    walmart$state_fips %in% c(1L, 5L, 12L, 13L, 28L, 45L) &
      walmart$year <= 1991L,
    , drop = FALSE
  ]

  expect_realdata_event_time_export_metadata(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    controls = c("income", "unemployrt", "poverty"),
    dataset_name = "Castle"
  )
  expect_realdata_event_time_export_metadata(
    data = walmart_overlap_panel,
    y = "log_retail_emp",
    ivar = "fips",
    tvar = "year",
    gvar = "g",
    controls = c("total_pop", "emp1964", "emp1977"),
    dataset_name = "Walmart"
  )
})
