# ============================================================================
# test-non-ra-bootstrap-diagnostics.R
# Bootstrap and overlap diagnostics for non-RA estimators.
# ============================================================================

make_extreme_ipw_cross_section <- function() {
  set.seed(20260525)
  n_treated <- 20L
  n_control <- 80L
  x <- c(
    stats::rnorm(n_treated, mean = 3, sd = 0.15),
    stats::rnorm(n_control - 3L, mean = 0, sd = 0.20),
    stats::rnorm(3L, mean = 3, sd = 0.10)
  )
  d <- c(rep(1L, n_treated), rep(0L, n_control))
  y <- 0.2 * x + 0.5 * d + stats::rnorm(n_treated + n_control, sd = 0.1)
  data.frame(y = y, d = d, x = x)
}

make_extreme_ipw_staggered_panel <- function() {
  cs <- make_extreme_ipw_cross_section()
  years <- 2000:2004
  panel <- expand.grid(id = seq_len(nrow(cs)), year = years)
  panel <- panel[order(panel$id, panel$year), ]
  panel$x <- cs$x[panel$id]
  panel$g <- ifelse(cs$d[panel$id] == 1L, 2003, Inf)
  panel$y <- 0.2 * panel$x +
    0.1 * (panel$year - min(panel$year)) +
    0.5 * (panel$g == 2003 & panel$year >= 2003) +
    stats::rnorm(nrow(panel), sd = 0.1)
  rownames(panel) <- NULL
  panel
}

test_that("IPW bootstrap handles glm warnings without losing SE", {
  df <- make_extreme_ipw_cross_section()
  overlap_warnings <- character()

  result <- withCallingHandlers(
    estimate_ipw(
      data = df,
      y = "y",
      d = "d",
      propensity_controls = "x",
      vce = "bootstrap",
      boot_reps = 20L,
      seed = 20260525
    ),
    warning = function(w) {
      if (inherits(w, "lwdid_overlap")) {
        overlap_warnings <<- c(overlap_warnings, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(result$vce_method, "bootstrap")
  expect_true(is.finite(result$se))
  expect_true(result$se > 0)
  expect_equal(result$bootstrap_reps_requested, 20L)
  expect_true(result$bootstrap_reps_valid >= 10L)
  expect_equal(
    result$bootstrap_success_rate,
    result$bootstrap_reps_valid / result$bootstrap_reps_requested,
    tolerance = 1e-12
  )
  expect_true(result$weights_cv > 2)
  expect_true(any(grepl("IPW weight CV", overlap_warnings, fixed = TRUE)))
})

test_that("staggered IPW reports aggregated high weight concentration", {
  panel <- make_extreme_ipw_staggered_panel()
  overlap_warnings <- character()

  fit <- withCallingHandlers(
    lwdid(
      data = panel,
      y = "y",
      ivar = "id",
      tvar = "year",
      gvar = "g",
      rolling = "demean",
      estimator = "ipw",
      controls = "x",
      aggregate = "none",
      return_diagnostics = TRUE,
      verbose = "quiet"
    ),
    warning = function(w) {
      if (inherits(w, "lwdid_overlap")) {
        overlap_warnings <<- c(overlap_warnings, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )

  propensity_diag <- get_diagnostics(fit, "propensity")

  expect_s3_class(fit, "lwdid_result")
  expect_true(max(propensity_diag$weights_cv, na.rm = TRUE) > 2)
  expect_true(any(grepl("max weight CV", overlap_warnings, fixed = TRUE)))
  expect_match(summary(fit)$diagnostics_summary$propensity, "max_weight_cv")
})

fit_castle_non_ra <- function(estimator, se_method = NULL, boot_reps = 20L,
                              aggregate = "none", event_time_range = NULL) {
  data("castle", package = "lwdid")
  args <- list(
    data = castle,
    y = "lhomicide",
    ivar = "sid",
    tvar = "year",
    gvar = "gvar",
    rolling = "demean",
    estimator = estimator,
    controls = c("income", "unemployrt", "poverty"),
    aggregate = aggregate,
    control_group = "not_yet_treated",
    return_diagnostics = TRUE,
    verbose = "quiet"
  )
  if (!is.null(event_time_range)) {
    args$event_time_range <- event_time_range
  }
  if (!is.null(se_method)) {
    args$se_method <- se_method
    args$boot_reps <- boot_reps
    args$seed <- 20260525L
  }
  suppressWarnings(do.call(lwdid, args))
}

expect_event_time_bootstrap_boundary <- function(fit) {
  event_rows <- extract_effects(fit, type = "event_time")
  contribution_rows <- extract_effects(fit, type = "event_time_contributions")
  event_vcov <- vcov(fit, type = "event_time")

  expect_true(nrow(event_rows) > 0L)
  expect_true(nrow(contribution_rows) >= nrow(event_rows))
  expect_true(all(c(
    "min_bootstrap_success_rate", "min_bootstrap_reps_valid",
    "max_bootstrap_reps_failed"
  ) %in% names(event_rows)))
  expect_true(all(c(
    "bootstrap_reps_requested", "bootstrap_reps_valid",
    "bootstrap_reps_failed", "bootstrap_success_rate"
  ) %in% names(contribution_rows)))
  expect_true(all(is.finite(event_rows$min_bootstrap_success_rate)))
  expect_true(all(event_rows$min_bootstrap_reps_valid >= 10L))
  expect_equal(
    unique(event_rows$se_aggregation),
    "diagonal_weighted_cohort_se"
  )
  expect_equal(
    unique(event_rows$covariance_assumption),
    "zero_cross_cohort_covariance"
  )
  expect_equal(attr(event_vcov, "se_aggregation"), "diagonal_weighted_cohort_se")
  expect_equal(
    attr(event_vcov, "covariance_assumption"),
    "zero_cross_cohort_covariance"
  )
  expect_equal(unname(diag(event_vcov)), event_rows$se^2, tolerance = 1e-12)
  if (length(event_vcov) > nrow(event_vcov)) {
    expect_true(all(event_vcov[row(event_vcov) != col(event_vcov)] == 0))
  }

  draw_like_names <- c(
    "bootstrap_draws", "bootstrap_att_draws", "boot_atts", "att_boot",
    "att_boots", "bootstrap_replicates", "bootstrap_rep_id"
  )
  expect_false(any(draw_like_names %in% names(fit)))
  expect_false(any(draw_like_names %in% names(fit$att_by_cohort_time)))
  expect_false(any(draw_like_names %in% names(event_rows)))
  expect_false(any(draw_like_names %in% names(contribution_rows)))
  expect_null(fit$vcov_att_event_time)
}

test_that("Castle IPWRA bootstrap keeps analytical row support", {
  analytical <- fit_castle_non_ra("ipwra")
  boot <- fit_castle_non_ra("ipwra", "bootstrap")

  analytical_keys <- analytical$att_by_cohort_time[, c("cohort", "period")]
  boot_keys <- boot$att_by_cohort_time[, c("cohort", "period")]
  boot_diag <- get_diagnostics(boot, "propensity")

  expect_equal(boot_keys, analytical_keys)
  expect_identical(unique(boot$att_by_cohort_time$vce_type), "bootstrap")
  expect_true(all(is.finite(boot$att_by_cohort_time$se)))
  expect_true(all(boot$att_by_cohort_time$bootstrap_reps_requested == 20L))
  expect_true(all(boot$att_by_cohort_time$bootstrap_reps_valid >= 10L))
  expect_equal(
    boot$att_by_cohort_time$bootstrap_success_rate,
    boot$att_by_cohort_time$bootstrap_reps_valid /
      boot$att_by_cohort_time$bootstrap_reps_requested,
    tolerance = 1e-12
  )
  expect_equal(nrow(boot_diag), nrow(boot$att_by_cohort_time))
  expect_true(all(c(
    "bootstrap_reps_requested", "bootstrap_reps_valid",
    "bootstrap_reps_failed", "bootstrap_success_rate"
  ) %in% names(boot_diag)))
  expect_true(max(boot_diag$weights_cv, na.rm = TRUE) > 2)
})

test_that("Castle PSM bootstrap keeps matched row support", {
  default <- fit_castle_non_ra("psm")
  boot <- fit_castle_non_ra("psm", "bootstrap", boot_reps = 50L)

  default_keys <- default$att_by_cohort_time[, c("cohort", "period")]
  boot_keys <- boot$att_by_cohort_time[, c("cohort", "period")]
  boot_diag <- get_diagnostics(boot, "propensity")

  expect_true(
    all(paste(default_keys$cohort, default_keys$period) %in%
          paste(boot_keys$cohort, boot_keys$period))
  )
  expect_gte(nrow(boot_keys), nrow(default_keys))
  expect_identical(unique(boot$att_by_cohort_time$vce_type), "bootstrap")
  expect_true(all(is.finite(boot$att_by_cohort_time$se)))
  expect_true(all(boot$att_by_cohort_time$bootstrap_reps_requested == 50L))
  expect_true(all(boot$att_by_cohort_time$bootstrap_reps_valid >= 10L))
  expect_equal(
    boot$att_by_cohort_time$bootstrap_success_rate,
    boot$att_by_cohort_time$bootstrap_reps_valid /
      boot$att_by_cohort_time$bootstrap_reps_requested,
    tolerance = 1e-12
  )
  expect_equal(nrow(boot_diag), nrow(boot$att_by_cohort_time))
  expect_true(all(c(
    "bootstrap_reps_requested", "bootstrap_reps_valid",
    "bootstrap_reps_failed", "bootstrap_success_rate"
  ) %in% names(boot_diag)))
  expect_true(all(is.finite(boot_diag$match_rate)))
  expect_true(all(boot_diag$match_rate > 0))
})

test_that("Castle IPWRA bootstrap event-time exports replicate diagnostics", {
  fit <- fit_castle_non_ra(
    "ipwra",
    se_method = "bootstrap",
    boot_reps = 20L,
    aggregate = "event_time",
    event_time_range = c(0, 2)
  )

  event_rows <- extract_effects(fit, type = "event_time")

  expect_event_time_bootstrap_boundary(fit)

  summary_file <- tempfile(fileext = ".csv")
  on.exit(unlink(summary_file), add = TRUE)
  summary_csv <- suppressMessages(
    to_csv(fit, file = summary_file, what = "summary")
  )
  dict <- to_dict(fit)

  expect_equal(
    summary_csv$event_time_bootstrap_n_event_time_rows,
    nrow(event_rows)
  )
  expect_equal(
    summary_csv$event_time_bootstrap_min_bootstrap_success_rate,
    min(event_rows$min_bootstrap_success_rate),
    tolerance = 1e-12
  )
  expect_equal(
    summary_csv$event_time_bootstrap_min_bootstrap_reps_valid,
    min(event_rows$min_bootstrap_reps_valid)
  )
  expect_equal(
    summary_csv$event_time_bootstrap_max_bootstrap_reps_failed,
    max(event_rows$max_bootstrap_reps_failed)
  )
  expect_equal(
    dict$event_time_bootstrap_summary$min_bootstrap_success_rate,
    min(event_rows$min_bootstrap_success_rate),
    tolerance = 1e-12
  )
  expect_equal(
    dict$event_time_bootstrap_summary$min_bootstrap_reps_valid,
    min(event_rows$min_bootstrap_reps_valid)
  )
})

test_that("Castle IPW and PSM bootstrap event-time vcov remains diagonal", {
  fits <- list(
    ipw = fit_castle_non_ra(
      "ipw",
      se_method = "bootstrap",
      boot_reps = 20L,
      aggregate = "event_time",
      event_time_range = c(0, 2)
    ),
    psm = fit_castle_non_ra(
      "psm",
      se_method = "bootstrap",
      boot_reps = 50L,
      aggregate = "event_time",
      event_time_range = c(0, 2)
    )
  )

  lapply(fits, expect_event_time_bootstrap_boundary)
})
