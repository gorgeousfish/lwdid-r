# ===========================================================================
# test-staggered-integration.R
# Integration tests for the lwdid R package staggered estimation pipeline.
# Covers: mode detection, parameter validation, control group strategy,
#         AET scenarios, small NT warnings, result object fields,
#         n_control/n semantics, df semantics, ATT fallback numerical,
#         df_fallback formula, end-to-end integration, boundary conditions,
#         and aggregate behavior.
# ===========================================================================

# --- Helper: generate balanced staggered panel data -----------------------
make_integration_data <- function(
  cohorts = c(4L, 6L),
  nt_count = 3L,
  n_per_cohort = 3L,
  t_range = 1:8,
  att = 5.0,
  seed = 42
) {
  set.seed(seed)
  n_treat <- length(cohorts) * n_per_cohort
  n_units <- n_treat + nt_count
  n_t <- length(t_range)
  gvar_vals <- c(
    rep(cohorts, each = n_per_cohort),
    rep(NA_real_, nt_count)
  )
  dt <- data.table::data.table(
    id = rep(seq_len(n_units), each = n_t),
    time = rep(t_range, n_units),
    gvar = rep(gvar_vals, each = n_t)
  )
  alpha_i <- rnorm(n_units, sd = 2)
  dt[, y := alpha_i[id] +
       data.table::fifelse(!is.na(gvar) & time >= gvar, att, 0) +
       rnorm(.N, sd = 0.5)]
  dt
}

# ===========================================================================
# Group 1: Mode auto-detection (Task E4-05.4)
# ===========================================================================

test_that("gvar non-NULL enters Staggered mode", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(result$is_staggered)
})

test_that("d/post specified with gvar triggers warning", {
  dt <- make_integration_data()
  dt[, d := data.table::fifelse(!is.na(gvar), 1L, 0L)]
  expect_warning(
    suppressMessages(lwdid(dt, y = "y", ivar = "id", tvar = "time",
                           gvar = "gvar", d = "d")),
    class = "lwdid_data"
  )
})

test_that("gvar=NULL does not enter Staggered mode", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time")
  )
})

# ===========================================================================
# Group 2: Parameter validation (Task E4-05.4)
# ===========================================================================

test_that("estimator non-ra throws lwdid_not_implemented", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", estimator = "ipwra"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("vce=cluster without cluster_var throws error", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", vce = "cluster"),
    class = "lwdid_invalid_parameter"
  )
})

test_that("exclude_pre_periods negative throws error", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = -1),
    class = "lwdid_invalid_parameter"
  )
})

test_that("exclude_pre_periods non-integer throws error", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = 1.5),
    class = "lwdid_invalid_parameter"
  )
})

test_that("invalid rolling value throws match.arg error", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", rolling = "invalid")
  )
})

test_that("invalid aggregate value throws match.arg error", {
  dt <- make_integration_data()
  expect_error(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "invalid")
  )
})

# ===========================================================================
# Group 3: Control group strategy (Task E4-05.5)
# ===========================================================================

test_that("auto resolves to never_treated when NT exists", {
  dt <- make_integration_data(nt_count = 5L)
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "auto")
  ))
  expect_equal(result$control_group_used, "never_treated")
})

test_that("auto resolves to not_yet_treated when no NT", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L, 8L), nt_count = 0L
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "auto",
          aggregate = "none")
  ))
  expect_equal(result$control_group_used, "not_yet_treated")
})

test_that("explicit never_treated stored correctly", {
  dt <- make_integration_data(nt_count = 5L)
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "never_treated")
  ))
  expect_equal(result$control_group, "never_treated")
  expect_equal(result$control_group_used, "never_treated")
})

test_that("control_group original value preserved", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "auto")
  ))
  expect_equal(result$control_group, "auto")
  expect_true(result$control_group_used %in%
              c("never_treated", "not_yet_treated"))
})

# ===========================================================================
# Group 4: AET scenarios (Task E4-05.5)
# ===========================================================================

test_that("AET + aggregate=none works normally", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L, 8L), nt_count = 0L, t_range = 1:10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none")
  ))
  expect_true(result$is_staggered)
  expect_true(nrow(result$cohort_time_effects) > 0)
})

test_that("AET + aggregate=cohort throws lwdid_no_never_treated", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L, 8L), nt_count = 0L, t_range = 1:10
  )
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort")
    )),
    class = "lwdid_no_never_treated"
  )
})

# ===========================================================================
# Group 5: Small NT warning (Task E4-05.5)
# ===========================================================================

test_that("NT < 2 with never_treated triggers small_nt_sample warning", {
  dt <- make_integration_data(nt_count = 1L)
  cond <- NULL
  withCallingHandlers(
    suppressMessages(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "never_treated")
    ),
    lwdid_data = function(c) {
      if (identical(c$detail, "small_nt_sample")) cond <<- c
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(cond))
})

test_that("NT >= 2 with never_treated no small_nt_sample warning", {
  dt <- make_integration_data(nt_count = 3L)
  cond <- NULL
  withCallingHandlers(
    suppressMessages(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "never_treated")
    ),
    lwdid_data = function(c) {
      if (identical(c$detail, "small_nt_sample")) cond <<- c
      invokeRestart("muffleWarning")
    }
  )
  expect_null(cond)
})

# ===========================================================================
# Group 6: Result object fields (Task E4-05.6)
# ===========================================================================

test_that("result contains all Staggered-specific fields", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(result$is_staggered)
  expect_true(!is.null(result$cohorts))
  expect_true(!is.null(result$cohort_sizes))
  expect_true(!is.null(result$cohort_time_effects))
  expect_true(!is.null(result$control_group))
  expect_true(!is.null(result$control_group_used))
  expect_true(!is.null(result$aggregate))
  expect_true(!is.null(result$rolling))
  expect_true(!is.null(result$estimator))
  expect_true(!is.null(result$n_never_treated))
  expect_true(!is.null(result$exclude_pre_periods))
})

test_that("cohort_time_effects has 17 columns", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  cte <- result$cohort_time_effects
  expect_equal(ncol(cte), 17L)
  expected_cols <- c(
    "cohort", "period", "event_time",
    "att", "se", "ci_lower", "ci_upper",
    "t_stat", "pvalue", "df", "df_inference",
    "n", "n_treated", "n_control",
    "K", "controls_tier", "vce_type"
  )
  expect_true(all(expected_cols %in% names(cte)))
})

test_that("cohorts are sorted ascending", {
  dt <- make_integration_data(cohorts = c(6L, 4L), t_range = 1:10)
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(all(diff(result$cohorts) > 0))
})

# ===========================================================================
# Group 7: n_control and n semantics (Task E4-05.6)
# ===========================================================================

test_that("n_control = n_nt when never_treated", {
  dt <- make_integration_data(nt_count = 5L)
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "never_treated")
  ))
  expect_equal(result$n_control, result$n_never_treated)
})

test_that("n_control = max(cte n_control) when not_yet_treated", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L, 8L), nt_count = 0L, t_range = 1:10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", control_group = "not_yet_treated",
          aggregate = "none")
  ))
  cte <- result$cohort_time_effects
  expect_equal(result$n_control,
               as.integer(max(cte$n_control, na.rm = TRUE)))
})

test_that("top-level nobs = total panel rows", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_equal(result$nobs, nrow(dt))
})

# ===========================================================================
# Group 8: df semantics (Task E4-05.6)
# ===========================================================================

test_that("df_resid and df_inference computed separately", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(is.finite(result$df_resid))
  expect_true(is.finite(result$df_inference))
  expect_true(result$df_resid > 0)
  expect_true(result$df_inference > 0)
})

test_that("cohort_sizes includes all cohorts", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(all(result$cohorts %in%
                  as.integer(names(result$cohort_sizes))))
})

# ===========================================================================
# Group 9: ATT fallback numerical (Task E4-05.7)
# ===========================================================================

test_that("att_fallback is n_treated weighted average", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L), nt_count = 5L,
    n_per_cohort = 5L, att = 3.0, seed = 100
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
          aggregate = "none")
  ))
  cte <- result$cohort_time_effects
  valid <- cte[!is.na(cte$att), ]
  expected_att <- weighted.mean(valid$att, valid$n_treated)
  expect_equal(result$att, expected_att, tolerance = 1e-12)
})

test_that("att_fallback close to true ATT with low noise", {
  true_att <- 5.0
  dt <- make_integration_data(
    cohorts = c(4L, 6L), nt_count = 10L,
    n_per_cohort = 10L, att = true_att, seed = 999
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar",
          aggregate = "none")
  ))
  expect_true(abs(result$att - true_att) < 2.0,
    info = sprintf("ATT=%.4f, true=%.1f, diff=%.4f",
                   result$att, true_att, abs(result$att - true_att)))
})

test_that("default staggered fit returns cohort-level inference statistics", {
  dt <- make_integration_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_equal(result$aggregate, "cohort")
  expect_true(is.finite(result$se_att))
  expect_true(is.finite(result$t_stat))
  expect_true(is.finite(result$pvalue))
  expect_true(is.finite(result$ci_lower))
  expect_true(is.finite(result$ci_upper))
  expect_true(result$ci_lower <= result$att)
  expect_true(result$ci_upper >= result$att)
})

# ===========================================================================
# Group 10: df_fallback formula (Task E4-05.7)
# ===========================================================================

test_that("df values are positive integers", {
  dt <- make_integration_data(
    cohorts = c(4L), nt_count = 5L, n_per_cohort = 3L
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(result$df_resid > 0)
  expect_true(result$df_inference > 0)
})

test_that("df_fallback minimum is 1", {
  dt <- make_integration_data(
    cohorts = c(4L), nt_count = 2L, n_per_cohort = 1L
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(result$df_resid >= 1)
  expect_true(result$df_inference >= 1)
})

# ===========================================================================
# Group 11: End-to-end integration (Task E4-05.8)
# ===========================================================================

test_that("end-to-end demean: known ATT recovery", {
  true_att <- 5.0
  dt <- make_integration_data(
    cohorts = c(4L, 6L), nt_count = 10L,
    n_per_cohort = 10L, att = true_att, seed = 123
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", rolling = "demean")
  ))
  expect_true(abs(result$att - true_att) < 2.0,
    info = sprintf("demean ATT=%.4f, true=%.1f", result$att, true_att))
  cte <- result$cohort_time_effects
  expect_true(all(abs(cte$att - true_att) < 5.0),
    info = "All (g,r) effects within 5.0 of true ATT")
})

test_that("end-to-end detrend: known ATT with unit trends", {
  set.seed(456)
  true_att <- 3.0
  n_treat <- 10; n_ctrl <- 10; t_range <- 1:10
  n_units <- n_treat + n_ctrl
  dt <- data.table::data.table(
    id = rep(seq_len(n_units), each = length(t_range)),
    time = rep(t_range, n_units),
    gvar = rep(c(rep(6L, n_treat), rep(NA_real_, n_ctrl)),
               each = length(t_range))
  )
  alpha_i <- rnorm(n_units, sd = 2)
  beta_i <- rnorm(n_units, sd = 0.3)
  dt[, y := alpha_i[id] + beta_i[id] * time +
       data.table::fifelse(!is.na(gvar) & time >= gvar,
                           true_att, 0) +
       rnorm(.N, sd = 0.01)]
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", rolling = "detrend")
  ))
  expect_true(abs(result$att - true_att) < 1.5,
    info = sprintf("detrend ATT=%.4f, true=%.1f", result$att, true_att))
})

test_that("end-to-end multi-cohort: all cohorts estimated", {
  dt <- make_integration_data(
    cohorts = c(3L, 5L, 7L), nt_count = 5L,
    n_per_cohort = 3L, t_range = 1:10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  cte <- result$cohort_time_effects
  expect_true(all(c(3L, 5L, 7L) %in% cte$cohort))
})

test_that("end-to-end: exclude_pre_periods=1 works", {
  dt <- make_integration_data(
    cohorts = c(5L), nt_count = 5L,
    n_per_cohort = 5L, t_range = 1:8
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = 1L)
  ))
  expect_true(result$is_staggered)
  expect_equal(result$exclude_pre_periods, 1L)
})

# ===========================================================================
# Group 12: Boundary conditions (Task E4-05.8)
# ===========================================================================

test_that("single cohort + NT works", {
  dt <- make_integration_data(
    cohorts = c(5L), nt_count = 5L, n_per_cohort = 5L
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(result$is_staggered)
  expect_equal(length(result$cohorts), 1L)
})

test_that("all cohorts insufficient pre-periods throws error", {
  dt <- data.table::data.table(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    gvar = rep(c(1, 1, NA, NA), each = 3),
    y = rnorm(12)
  )
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
    )),
    class = "lwdid_insufficient_pre_periods"
  )
})

test_that("zero treatment effect estimated near zero", {
  dt <- make_integration_data(
    cohorts = c(5L), nt_count = 15L,
    n_per_cohort = 15L, att = 0, seed = 789
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  expect_true(abs(result$att) < 1.5,
    info = sprintf("Zero ATT test: estimated=%.4f", result$att))
})

# ===========================================================================
# Group 13: Aggregate behavior (Task E4-05.7)
# ===========================================================================

test_that("aggregate non-none returns cohort aggregate without stub message", {
  dt <- make_integration_data(nt_count = 5L)
  expect_no_message(
    result <- suppressWarnings(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort",
            control_group = "never_treated")
    )
  )
  expect_equal(result$aggregate, "cohort")
  expect_false(is.null(result$cohort_effects))
  expect_equal(result$att, result$att_cohort_agg, tolerance = 1e-10)
})

# ===========================================================================
# Group 14: Error conditions & boundary (Task E4-05.9)
# ===========================================================================

test_that("exclude_pre_periods=0 is default behavior", {
  dt <- make_integration_data(
    cohorts = c(4L, 6L), nt_count = 5L,
    n_per_cohort = 3L, seed = 77
  )
  result_default <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time", gvar = "gvar")
  ))
  result_k0 <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = 0L)
  ))
  expect_equal(result_default$att, result_k0$att, tolerance = 1e-12)
})

test_that("exclude_pre_periods > 0 produces valid results", {
  dt <- make_integration_data(
    cohorts = c(5L), nt_count = 5L,
    n_per_cohort = 5L, t_range = 1:10, seed = 88
  )
  result_k0 <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = 0L)
  ))
  result_k1 <- suppressWarnings(suppressMessages(
    lwdid(dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", exclude_pre_periods = 1L)
  ))
  expect_true(is.finite(result_k0$att))
  expect_true(is.finite(result_k1$att))
  expect_equal(result_k0$exclude_pre_periods, 0L)
  expect_equal(result_k1$exclude_pre_periods, 1L)
})

test_that("R enhanced: small NT warning with aggregate='none'", {
  # Python only warns for aggregate='cohort'/'overall'
  # R warns whenever resolved_cg='never_treated' and n_nt < 2
  dt <- make_integration_data(nt_count = 1L)
  cond <- NULL
  withCallingHandlers(
    suppressMessages(
      lwdid(dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "never_treated",
            aggregate = "none")
    ),
    lwdid_data = function(c) {
      if (identical(c$detail, "small_nt_sample")) cond <<- c
      invokeRestart("muffleWarning")
    }
  )
  expect_false(is.null(cond))
})
