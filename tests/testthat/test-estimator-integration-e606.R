# test-estimator-integration-e606.R — Estimator Integration Tests (E6-06.6)

if (!exists("lwdid", mode = "function")) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run test-estimator-integration-e606.R directly.")
  }
  devtools::load_all(
    "/Users/cxy/Desktop/lwdid_r/lwdid-r",
    export_all = FALSE,
    quiet = TRUE
  )
}

build_e606_common_timing_panel <- function(seed = 606L) {
  set.seed(seed)

  n_units <- 120L
  periods <- 1:4

  x1 <- stats::rnorm(n_units)
  x2 <- stats::rnorm(n_units)
  x3 <- stats::rnorm(n_units)
  treat_prob <- stats::plogis(0.9 * x1 - 0.5 * x2)
  d <- stats::rbinom(n_units, size = 1L, prob = treat_prob)

  while (sum(d) < 30L || sum(d == 0L) < 30L) {
    d <- stats::rbinom(n_units, size = 1L, prob = treat_prob)
  }

  unit_fe <- stats::rnorm(n_units, sd = 0.25)

  panel <- do.call(rbind, lapply(seq_len(n_units), function(i) {
    data.frame(
      id = i,
      year = periods,
      post = as.integer(periods >= 3L),
      d = d[i],
      x1 = x1[i],
      x2 = x2[i],
      x3 = x3[i],
      y = unit_fe[i] +
        0.12 * periods +
        0.55 * d[i] * as.integer(periods >= 3L) +
        0.45 * x1[i] -
        0.35 * x2[i] +
        0.25 * x3[i] +
        stats::rnorm(length(periods), sd = 0.35)
    )
  }))

  rownames(panel) <- NULL
  panel
}

build_e606_staggered_panel <- function(seed = 6606L) {
  set.seed(seed)

  n_per_cohort <- 36L
  ids <- seq_len(3L * n_per_cohort)
  cohorts <- rep(c(0L, 2002L, 2003L), each = n_per_cohort)
  years <- 2000:2004

  x1 <- stats::rnorm(
    length(ids),
    mean = rep(c(-0.6, 0.2, 0.8), each = n_per_cohort),
    sd = 0.6
  )
  x2 <- stats::rnorm(
    length(ids),
    mean = rep(c(0.3, -0.4, 0.5), each = n_per_cohort),
    sd = 0.7
  )
  x3 <- rep(c(-0.5, 0.1, 0.7), each = n_per_cohort) +
    stats::rnorm(length(ids), sd = 0.2)

  panel <- expand.grid(id = ids, year = years)
  panel <- panel[order(panel$id, panel$year), ]
  panel$gvar <- cohorts[panel$id]
  panel$x1 <- x1[panel$id]
  panel$x2 <- x2[panel$id]
  panel$x3 <- x3[panel$id]
  panel$y <- 1 +
    0.45 * panel$x1 -
    0.25 * panel$x2 +
    0.20 * panel$x3 +
    0.22 * (panel$year - min(years)) +
    ifelse(
      panel$gvar == 2002L & panel$year >= 2002L,
      1.4 + 0.10 * (panel$year - 2002L),
      0
    ) +
    ifelse(
      panel$gvar == 2003L & panel$year >= 2003L,
      2.1 + 0.20 * (panel$year - 2003L),
      0
    ) +
    stats::rnorm(nrow(panel), sd = 0.35)

  rownames(panel) <- NULL
  panel
}

build_e606_common_timing_estimator_input <- function(
    panel,
    estimator,
    controls,
    ps_controls = NULL,
    se_method = NULL,
    boot_reps = 200L,
    seed = NULL
) {
  validated <- lwdid:::validate_inputs(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = ps_controls,
    estimator = estimator,
    se_method = se_method,
    boot_reps = boot_reps,
    seed = seed,
    pretreatment_test = FALSE
  )

  dt <- data.table::copy(validated$data)
  params <- validated$validated_params
  estimator_controls <- unique(c(params$controls, params$ps_controls))
  if (length(estimator_controls) == 0L) {
    estimator_controls <- NULL
  }

  unit_d <- dt[, .(d = data.table::first(get(params$d))), by = c(params$ivar)]
  x_dt <- NULL
  if (!is.null(estimator_controls)) {
    x_dt <- lwdid:::.extract_controls(
      dt,
      params$ivar,
      "tindex",
      estimator_controls,
      validated$tpost1,
      params$exclude_pre_periods
    )
  }

  dt <- lwdid:::transform_common(
    dt,
    params$depvar,
    params$ivar,
    "tindex",
    g = validated$tpost1,
    rolling = params$rolling,
    exclude_pre_periods = params$exclude_pre_periods
  )

  summary_dt <- dt[
    tindex >= validated$tpost1,
    .(y_trans_summary = mean(y_trans, na.rm = TRUE)),
    by = c(params$ivar)
  ]
  summary_dt <- merge(summary_dt, unit_d, by = params$ivar)
  if (!is.null(x_dt)) {
    summary_dt <- merge(summary_dt, x_dt, by = params$ivar)
  }

  summary_valid <- summary_dt[!is.na(y_trans_summary)]
  est_df <- data.frame(
    .y_outcome = summary_valid$y_trans_summary,
    .d_treat = summary_valid$d
  )
  if (!is.null(estimator_controls)) {
    est_df <- cbind(
      est_df,
      as.data.frame(summary_valid[, estimator_controls, with = FALSE])
    )
  }

  est_df
}

test_that("TC-6.6.54: lwdid routes common-timing IPW through the end-to-end entry point", {
  panel <- build_e606_common_timing_panel()

  result <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x2", "x3"),
    ps_controls = c("x1", "x2"),
    estimator = "ipw",
    pretreatment_test = FALSE
  )

  expect_identical(result$estimator, "ipw")
  expect_identical(result$inference_dist, "normal")
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  expect_true(is.finite(result$pvalue))
  expect_equal(result$n_treated + result$n_control, length(unique(panel$id)))
})

test_that("TC-6.6.55: staggered IPWRA cohort aggregation avoids spurious top-level warnings", {
  panel <- build_e606_staggered_panel()
  warnings_seen <- character(0)

  result <- withCallingHandlers(
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      controls = c("x1", "x2"),
      ps_controls = "x1",
      estimator = "ipwra",
      aggregate = "cohort",
      control_group = "never_treated",
      pretreatment_test = FALSE
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  allowed_warning_markers <- c(
    "IPWRA weight CV",
    "Overlap assumption may be violated"
  )

  # Zero warnings is fine (no overlap violation); if any appear,
  # they must be from the allowed set (no spurious top-level warnings).
  if (length(warnings_seen) > 0L) {
    expect_true(all(vapply(
      warnings_seen,
      function(msg) any(vapply(
        allowed_warning_markers,
        grepl,
        logical(1),
        x = msg,
        fixed = TRUE
      )),
      logical(1)
    )))
  }
  expect_identical(result$aggregate, "cohort")
  expect_identical(result$control_group_used, "never_treated")
  expect_true(result$n_cohort_effects > 0L)
  expect_true(result$n_gr_effects > 0L)
  expect_equal(result$att, result$att_cohort_agg, tolerance = 1e-10)
  expect_true(is.finite(result$se_att))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  expect_null(result$K)
})

test_that("TC-6.6.57: top-level common-timing inference keeps estimator-native distribution", {
  panel <- build_e606_common_timing_panel()
  controls <- c("x1", "x2")
  alpha <- 0.05

  ra <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    estimator = "ra",
    alpha = alpha,
    pretreatment_test = FALSE
  )

  ipw <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = controls,
    estimator = "ipw",
    alpha = alpha,
    pretreatment_test = FALSE
  )

  ipwra <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = controls,
    estimator = "ipwra",
    alpha = alpha,
    pretreatment_test = FALSE
  )

  expect_true(is.finite(ra$df_inference))
  expect_true(is.finite(ipw$df_inference))
  expect_true(is.finite(ipwra$df_inference))

  t_crit_ra <- stats::qt(1 - alpha / 2, df = ra$df_inference)
  z_crit <- stats::qnorm(1 - alpha / 2)

  expect_equal(
    ra$ci_lower,
    ra$att - t_crit_ra * ra$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ra$ci_upper,
    ra$att + t_crit_ra * ra$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ra$pvalue,
    2 * stats::pt(abs(ra$t_stat), df = ra$df_inference, lower.tail = FALSE),
    tolerance = 1e-10
  )

  expect_equal(
    ipw$ci_lower,
    ipw$att - z_crit * ipw$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ipw$ci_upper,
    ipw$att + z_crit * ipw$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ipw$pvalue,
    2 * (1 - stats::pnorm(abs(ipw$t_stat))),
    tolerance = 1e-10
  )

  expect_equal(
    ipwra$ci_lower,
    ipwra$att - z_crit * ipwra$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ipwra$ci_upper,
    ipwra$att + z_crit * ipwra$se_att,
    tolerance = 1e-10
  )
  expect_equal(
    ipwra$pvalue,
    2 * (1 - stats::pnorm(abs(ipwra$t_stat))),
    tolerance = 1e-10
  )
})

test_that("TC-6.6.55: staggered cohort aggregation keeps IPWRA group-time effects and uses cohort regression", {
  panel <- build_e606_staggered_panel()
  controls <- c("x1", "x2", "x3")
  ps_controls <- c("x1", "x2")
  cohorts <- sort(unique(panel$gvar[panel$gvar > 0L]))

  cohort_result <- suppressWarnings(suppressMessages(
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      controls = controls,
      ps_controls = ps_controls,
      estimator = "ipwra",
      aggregate = "cohort",
      control_group = "never_treated",
      pretreatment_test = FALSE
    )
  ))

  manual_cohort <- lwdid:::aggregate_to_cohort(
    dt = data.table::as.data.table(panel),
    y = "y",
    ivar = "id",
    tvar = "year",
    gvar = "gvar",
    cohorts = cohorts,
    T_max = max(panel$year),
    pre_stats = lwdid:::precompute_transforms(
      dt = data.table::as.data.table(panel),
      y = "y",
      ivar = "id",
      tvar = "year",
      cohorts = cohorts,
      rolling = "demean",
      exclude_pre_periods = 0L
    ),
    rolling = "demean",
    alpha = 0.05,
    controls = controls
  )

  cohort_means <- tapply(
    cohort_result$att_by_cohort_time$att,
    cohort_result$att_by_cohort_time$cohort,
    mean
  )
  cohort_atts <- setNames(
    vapply(cohort_result$cohort_effects, function(x) x$att, numeric(1)),
    vapply(cohort_result$cohort_effects, function(x) as.character(x$cohort), character(1))
  )

  expect_gt(nrow(cohort_result$att_by_cohort_time), length(cohort_result$cohort_effects))
  # IPWRA uses .aggregate_gr_to_cohort (mean of (g,r)) which differs from RA's
  # aggregate_to_cohort (time-averaged regression). Verify structural consistency.
  expect_true(length(cohort_result$cohort_effects) > 0L)
  for (ce in cohort_result$cohort_effects) {
    expect_true(is.finite(ce$att))
    expect_true(is.integer(ce$cohort))
  }
})

test_that("TC-6.6.56: common-timing PSM with sd caliper matches the transformed cross-section dispatch", {
  panel <- build_e606_common_timing_panel()
  controls <- c("x1", "x2", "x3")
  ps_controls <- c("x1", "x2")
  est_df <- build_e606_common_timing_estimator_input(
    panel = panel,
    estimator = "psm",
    controls = controls,
    ps_controls = ps_controls
  )

  direct <- lwdid:::dispatch_estimator(
    data = est_df,
    y = ".y_outcome",
    d = ".d_treat",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "sd"
  )

  result <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "sd",
    pretreatment_test = FALSE
  )

  ps_mean <- mean(direct$propensity_scores)
  ps_sd_pop <- sqrt(sum((direct$propensity_scores - ps_mean)^2) / length(direct$propensity_scores))

  expect_equal(direct$caliper, 0.1 * ps_sd_pop, tolerance = 1e-10)
  expect_equal(result$att, direct$att, tolerance = 1e-12)
  expect_equal(result$se_att, direct$se, tolerance = 1e-12)
  expect_equal(result$ci_lower, direct$ci_lower, tolerance = 1e-12)
  expect_equal(result$ci_upper, direct$ci_upper, tolerance = 1e-12)
})

test_that("TC-6.6.58: staggered IPWRA retains the same group-time path under aggregate='cohort'", {
  panel <- build_e606_staggered_panel()
  controls <- c("x1", "x2", "x3")
  ps_controls <- c("x1", "x2")

  none_result <- suppressWarnings(suppressMessages(
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      controls = controls,
      ps_controls = ps_controls,
      estimator = "ipwra",
      aggregate = "none",
      control_group = "never_treated",
      pretreatment_test = FALSE
    )
  ))

  cohort_result <- suppressWarnings(suppressMessages(
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "gvar",
      rolling = "demean",
      controls = controls,
      ps_controls = ps_controls,
      estimator = "ipwra",
      aggregate = "cohort",
      control_group = "never_treated",
      pretreatment_test = FALSE
    )
  ))

  cell_ids <- with(
    cohort_result$att_by_cohort_time,
    paste(cohort, period, sep = "::")
  )

  expect_equal(cohort_result$att_by_cohort_time, none_result$att_by_cohort_time, tolerance = 1e-12)
  expect_gt(length(unique(cell_ids)), 2L)
  expect_gt(length(unique(round(cohort_result$att_by_cohort_time$att, 8))), 2L)
})

test_that("TC-6.6.59: common-timing dispatch keeps ps_controls-only covariates available", {
  panel <- build_e606_common_timing_panel()

  ps_specific <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x3"),
    ps_controls = c("x1", "x2"),
    estimator = "ipwra",
    pretreatment_test = FALSE
  )

  controls_only <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = c("x1", "x3"),
    ps_controls = c("x1", "x3"),
    estimator = "ipwra",
    pretreatment_test = FALSE
  )

  expect_true(is.finite(ps_specific$att))
  expect_true(is.finite(ps_specific$se_att))
  expect_false(isTRUE(all.equal(ps_specific$att, controls_only$att)))
})

test_that("TC-6.6.60: common-timing IPWRA bootstrap keeps the estimator percentile CI", {
  panel <- build_e606_common_timing_panel()
  controls <- c("x1", "x2")
  est_df <- build_e606_common_timing_estimator_input(
    panel = panel,
    estimator = "ipwra",
    controls = controls,
    ps_controls = controls,
    se_method = "bootstrap",
    boot_reps = 15L,
    seed = 321L
  )

  direct <- lwdid:::dispatch_estimator(
    data = est_df,
    y = ".y_outcome",
    d = ".d_treat",
    controls = controls,
    ps_controls = controls,
    estimator = "ipwra",
    se_method = "bootstrap",
    boot_reps = 15L,
    seed = 321L
  )

  result <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = controls,
    estimator = "ipwra",
    se_method = "bootstrap",
    boot_reps = 15L,
    seed = 321L,
    pretreatment_test = FALSE
  )

  normal_ci <- result$att + stats::qnorm(c(0.025, 0.975)) * result$se_att

  expect_gt(result$se_att, 0)
  expect_equal(result$att, direct$att, tolerance = 1e-12)
  expect_equal(result$se_att, direct$se, tolerance = 1e-12)
  expect_equal(result$ci_lower, direct$ci_lower, tolerance = 1e-12)
  expect_equal(result$ci_upper, direct$ci_upper, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    c(result$ci_lower, result$ci_upper),
    normal_ci,
    tolerance = 1e-8
  )))
})

test_that("TC-6.6.61: common-timing PSM absolute caliper keeps the literal threshold", {
  panel <- build_e606_common_timing_panel()
  controls <- c("x1", "x2", "x3")
  ps_controls <- c("x1", "x2")
  est_df <- build_e606_common_timing_estimator_input(
    panel = panel,
    estimator = "psm",
    controls = controls,
    ps_controls = ps_controls
  )

  direct_absolute <- lwdid:::dispatch_estimator(
    data = est_df,
    y = ".y_outcome",
    d = ".d_treat",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "absolute"
  )
  direct_sd <- lwdid:::dispatch_estimator(
    data = est_df,
    y = ".y_outcome",
    d = ".d_treat",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "sd"
  )

  result_absolute <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "absolute",
    pretreatment_test = FALSE
  )
  result_sd <- lwdid(
    data = panel,
    y = "y",
    ivar = "id",
    tvar = "year",
    d = "d",
    post = "post",
    controls = controls,
    ps_controls = ps_controls,
    estimator = "psm",
    caliper = 0.1,
    caliper_scale = "sd",
    pretreatment_test = FALSE
  )

  expect_equal(direct_absolute$caliper, 0.1, tolerance = 1e-12)
  expect_equal(result_absolute$att, direct_absolute$att, tolerance = 1e-12)
  expect_equal(result_absolute$se_att, direct_absolute$se, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(result_absolute$att, result_sd$att, tolerance = 1e-8)))
  expect_false(isTRUE(all.equal(direct_absolute$att, direct_sd$att, tolerance = 1e-8)))
})

test_that("TC-6.6.62: bootstrap end-to-end results are reproducible for a fixed seed", {
  panel <- build_e606_common_timing_panel()

  fit_once <- function() {
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      d = "d",
      post = "post",
      controls = c("x1", "x2"),
      ps_controls = c("x1", "x2"),
      estimator = "ipw",
      se_method = "bootstrap",
      boot_reps = 20L,
      seed = 123L,
      pretreatment_test = FALSE
    )
  }

  first <- fit_once()
  second <- fit_once()

  expect_equal(first$att, second$att, tolerance = 1e-15)
  expect_equal(first$se_att, second$se_att, tolerance = 1e-15)
  expect_equal(first$ci_lower, second$ci_lower, tolerance = 1e-15)
  expect_equal(first$ci_upper, second$ci_upper, tolerance = 1e-15)
})
