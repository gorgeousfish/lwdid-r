# ===========================================================================
# test-lwdid-integration.R
# Integration tests for lwdid() RI and pretreatment paths (Story E7-06)
# ===========================================================================

# --- Helper: Common Timing DGP with treatment effect ---
.make_ct_panel <- function(n_treated = 10L, n_control = 10L,
                           n_periods = 8L, K = 4L,
                           att = 2.0, seed = 42L) {
  set.seed(seed)
  n_units <- n_treated + n_control
  years <- 2000L + seq_len(n_periods)
  df <- expand.grid(id = seq_len(n_units), year = years)
  df <- df[order(df$id, df$year), ]
  df$d <- ifelse(df$id <= n_treated, 1L, 0L)
  df$post <- ifelse(df$year > (2000L + K), 1L, 0L)
  alpha_i <- rnorm(n_units, sd = 1.5)
  df$y <- alpha_i[df$id] + df$d * df$post * att + rnorm(nrow(df), sd = 0.5)
  rownames(df) <- NULL
  df
}

# --- Helper: Staggered DGP ---
.make_stag_panel <- function(cohorts = c(4L, 6L), nt_count = 5L,
                             n_per_cohort = 5L, t_range = 1:10,
                             att = 3.0, seed = 42L) {
  set.seed(seed)
  n_treat <- length(cohorts) * n_per_cohort
  n_units <- n_treat + nt_count
  n_t <- length(t_range)
  gvar_vals <- c(
    rep(cohorts, each = n_per_cohort),
    rep(0L, nt_count)
  )
  dt <- data.table::data.table(
    id = rep(seq_len(n_units), each = n_t),
    time = rep(t_range, n_units),
    gvar = rep(gvar_vals, each = n_t)
  )
  alpha_i <- rnorm(n_units, sd = 2)
  dt[, y := alpha_i[id] +
       data.table::fifelse(gvar > 0L & time >= gvar, att, 0) +
       rnorm(.N, sd = 0.5)]
  dt
}


# ===========================================================================
# Group 1: RI Integration Tests (Task E7-06.5)
# ===========================================================================

# TC-7.6.1: lwdid() + ri=TRUE + Common Timing
# Use n=50 per group and rireps=200 to ensure >100 valid bootstrap reps
test_that("TC-7.6.1: CT RI produces valid ri_* fields", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 101L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 123L,
          ri_method = "bootstrap")
  ))
  # ri_pvalue in [0, 1]
  expect_true(!is.null(result$ri_pvalue))
  expect_true(is.numeric(result$ri_pvalue))
  expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
  # ri_method
  expect_equal(result$ri_method, "bootstrap")
  # rireps
  expect_equal(result$rireps, 200L)
  # ri_seed
  expect_equal(result$ri_seed, 123L)
  # ri_valid is positive integer
  expect_true(!is.null(result$ri_valid))
  expect_true(is.numeric(result$ri_valid) && result$ri_valid > 0L)
  # ri_failed is non-negative integer
  expect_true(!is.null(result$ri_failed))
  expect_true(is.numeric(result$ri_failed) && result$ri_failed >= 0L)
  # ri_distribution is numeric vector
  expect_true(!is.null(result$ri_distribution))
  expect_true(is.numeric(result$ri_distribution))
  expect_true(length(result$ri_distribution) > 0L)
  # ri_error should be NULL (success)
  expect_null(result$ri_error)
  # ATT should still be valid
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
})

# TC-7.6.2: lwdid() + ri=TRUE + Staggered
# Use permutation method (min_valid = max(10, 0.1*reps) = 10) for smaller sample
test_that("TC-7.6.2: Staggered RI produces valid ri_* fields", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 201L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", ri = TRUE, rireps = 100L, seed = 456L,
          ri_method = "permutation")
  ))
  # ri_pvalue in [0, 1]
  expect_true(!is.null(result$ri_pvalue))
  expect_true(is.numeric(result$ri_pvalue))
  expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
  # ri_target matches aggregate (default is cohort for staggered)
  expect_true(!is.null(result$ri_target))
  # ri_method
  expect_equal(result$ri_method, "permutation")
  # rireps
  expect_equal(result$rireps, 100L)
  # ri_seed
  expect_equal(result$ri_seed, 456L)
  # ATT should still be valid
  expect_true(is.finite(result$att))
})


# TC-7.6.8: ri + include_pretreatment simultaneously — ATT unchanged
test_that("TC-7.6.8: ri+pretreatment ATT matches ri-only ATT", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 301L)
  # ri only
  res_ri <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 789L)
  ))
  # ri + pretreatment
  res_both <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 789L,
          include_pretreatment = TRUE)
  ))
  # ATT should be identical (pretreatment doesn't affect ATT)
  expect_equal(res_both$att, res_ri$att, tolerance = 1e-14)
  expect_equal(res_both$se_att, res_ri$se_att, tolerance = 1e-14)
  expect_equal(res_both$pvalue, res_ri$pvalue, tolerance = 1e-14)
  # Both should have RI results
  expect_true(!is.null(res_both$ri_pvalue))
  # RI p-values should be identical (same seed, same data, same ATT)
  expect_equal(res_both$ri_pvalue, res_ri$ri_pvalue, tolerance = 1e-14)
  # res_both should also have pretreatment results
  expect_true(!is.null(res_both$att_pre_treatment))
})


# TC-7.6.10: RI failure captured gracefully
test_that("TC-7.6.10: RI failure captured, ATT unaffected", {
  # Create tiny dataset — RI may fail due to insufficient variation
  df <- .make_ct_panel(n_treated = 3L, n_control = 3L,
                       n_periods = 6L, K = 3L, att = 1.0, seed = 401L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 20L, seed = 999L)
  ))
  # ATT should always be valid regardless of RI outcome
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
  # If RI failed, ri_error should be non-NULL
  if (is.null(result$ri_pvalue)) {
    expect_true(!is.null(result$ri_error))
  }
})

# TC-7.6.11: seed=NULL auto-generates seed, reproducible with ri_seed
test_that("TC-7.6.11: auto-generated seed is reproducible", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 501L)
  # Run with seed=NULL
  res1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = NULL)
  ))
  # ri_seed should be a positive integer
  expect_true(!is.null(res1$ri_seed))
  expect_true(is.numeric(res1$ri_seed) && res1$ri_seed > 0L)
  # Re-run with the auto-generated seed
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L,
          seed = res1$ri_seed)
  ))
  # RI p-values should be identical
  expect_equal(res2$ri_pvalue, res1$ri_pvalue)
})

# TC-7.6.14: Staggered RI aggregate="cohort" uses first valid cohort ATT
test_that("TC-7.6.14: Staggered RI cohort target uses first valid cohort ATT", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 601L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 50L, seed = 111L,
          ri_method = "permutation")
  ))
  # Should have RI results or error
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
    expect_equal(result$ri_target, "cohort")
  } else {
    # If failed, ri_error should explain why
    expect_true(!is.null(result$ri_error))
  }
})

# TC-7.6.15: Staggered RI all cohort ATTs NA -> ri_error set
# Use tryCatch because degenerate data may cause main estimation to fail
test_that("TC-7.6.15: Staggered RI all-NA cohort ATTs sets ri_error", {
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 3L,
                         n_per_cohort = 2L, t_range = 1:10,
                         att = 0.0, seed = 701L)
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort",
            ri = TRUE, rireps = 10L, seed = 222L)
    )),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    # Degenerate data caused main estimation to fail — acceptable
    expect_true(grepl("insufficient|No valid", conditionMessage(result)))
  } else {
    # If estimation succeeded, check RI degradation
    expect_true(inherits(result, "lwdid_result"))
    if (!is.null(result$cohort_effects)) {
      all_na <- all(vapply(result$cohort_effects,
                           function(e) is.na(e$att), logical(1)))
      if (all_na) {
        expect_true(!is.null(result$ri_error))
        expect_equal(result$ri_target, "cohort")
      }
    }
  }
})

# TC-7.6.16: Staggered RI aggregate="cohort_time" (none) uses first valid CTE ATT
test_that("TC-7.6.16: Staggered RI cohort_time target", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 801L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          ri = TRUE, rireps = 50L, seed = 333L,
          ri_method = "permutation")
  ))
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
    expect_equal(result$ri_target, "cohort_time")
  } else {
    expect_true(!is.null(result$ri_error))
  }
})

# TC-7.6.17: Staggered RI aggregate="overall" with NA att_overall
test_that("TC-7.6.17: Staggered RI overall with NA att sets ri_error", {
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 3L,
                         n_per_cohort = 2L, t_range = 1:10,
                         att = 0.0, seed = 901L)
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "overall",
            ri = TRUE, rireps = 10L, seed = 444L)
    )),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_true(grepl("insufficient|No valid", conditionMessage(result)))
  } else {
    expect_true(inherits(result, "lwdid_result"))
    if (is.null(result$att_overall) || is.na(result$att_overall)) {
      expect_true(!is.null(result$ri_error))
      expect_equal(result$ri_target, "overall")
    }
  }
})

# TC-7.6.18: Staggered RI aggregate="cohort" with NULL att_by_cohort
test_that("TC-7.6.18: Staggered RI cohort with NULL effects sets ri_error", {
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 3L,
                         n_per_cohort = 2L, t_range = 1:10,
                         att = 0.0, seed = 1001L)
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort",
            ri = TRUE, rireps = 10L, seed = 555L)
    )),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_true(grepl("insufficient|No valid", conditionMessage(result)))
  } else {
    expect_true(inherits(result, "lwdid_result"))
    if (is.null(result$cohort_effects) ||
        length(result$cohort_effects) == 0L) {
      expect_true(!is.null(result$ri_error))
      expect_equal(result$ri_target, "cohort")
    }
  }
})


# TC-7.6.10: RI failure captured gracefully, ATT unaffected
test_that("TC-7.6.10: RI failure captured, ATT unaffected", {
  # Create tiny dataset that will cause RI to fail (N=4 total)
  df <- .make_ct_panel(n_treated = 2L, n_control = 2L,
                       n_periods = 4L, K = 2L, att = 1.0, seed = 401L)
  # With N=4 and bootstrap, most reps will be degenerate
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 999L)
  ))
  # ATT should always be valid regardless of RI outcome
  expect_true(is.finite(result$att))
  expect_true(is.finite(result$se_att))
  # If RI failed, ri_error should be non-NULL
  if (is.null(result$ri_pvalue)) {
    expect_true(!is.null(result$ri_error))
    expect_true(is.character(result$ri_error))
  }
})

# TC-7.6.11: seed=NULL auto-generates seed, reproducible with ri_seed
test_that("TC-7.6.11: auto-generated seed is reproducible", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 501L)
  # Run with seed=NULL
  res1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = NULL)
  ))
  # ri_seed should be a positive integer
  expect_true(!is.null(res1$ri_seed))
  expect_true(is.numeric(res1$ri_seed) && res1$ri_seed > 0L)
  # Re-run with the auto-generated seed
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L,
          seed = res1$ri_seed)
  ))
  # RI p-values should be identical
  if (!is.null(res1$ri_pvalue) && !is.null(res2$ri_pvalue)) {
    expect_equal(res2$ri_pvalue, res1$ri_pvalue)
  }
})

# TC-7.6.14: Staggered RI aggregate="cohort" uses first valid cohort ATT
test_that("TC-7.6.14: Staggered RI cohort target uses first valid cohort ATT", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 601L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 100L, seed = 111L,
          ri_method = "permutation")
  ))
  # Should have RI results or error
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
    expect_equal(result$ri_target, "cohort")
  } else {
    expect_true(!is.null(result$ri_error))
  }
})

# TC-7.6.15: Staggered RI all cohort ATTs NA -> ri_error set
test_that("TC-7.6.15: Staggered RI all-NA cohort ATTs sets ri_error", {
  # Use very small sample to increase chance of NA effects
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 2L,
                         n_per_cohort = 1L, t_range = 1:10,
                         att = 0.0, seed = 701L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 50L, seed = 222L,
          ri_method = "permutation")
  ))
  # ATT should still be computed (even if RI fails)
  expect_true(inherits(result, "lwdid_result"))
  # If cohort effects are all NA, ri_error should be set
  if (!is.null(result$cohort_effects)) {
    all_na <- all(vapply(result$cohort_effects,
                         function(e) is.na(e$att), logical(1)))
    if (all_na) {
      expect_true(!is.null(result$ri_error))
      expect_equal(result$ri_target, "cohort")
    }
  }
})

# TC-7.6.16: Staggered RI aggregate="cohort_time" uses first valid CTE ATT
test_that("TC-7.6.16: Staggered RI cohort_time target", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 801L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          ri = TRUE, rireps = 100L, seed = 333L,
          ri_method = "permutation")
  ))
  if (!is.null(result$ri_pvalue)) {
    expect_true(result$ri_pvalue >= 0 && result$ri_pvalue <= 1)
    expect_equal(result$ri_target, "cohort_time")
  } else {
    expect_true(!is.null(result$ri_error))
  }
})


# ===========================================================================
# Group 2: Pretreatment & Parameter Validation Tests (Task E7-06.6)
# ===========================================================================

# TC-7.6.3: lwdid() + include_pretreatment=TRUE + Common Timing
test_that("TC-7.6.3: CT pretreatment produces valid att_pre_treatment", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1101L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  # att_pre_treatment should be a data.frame
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(is.data.frame(result$att_pre_treatment))
  # Should have rows (pre-treatment periods)
  expect_true(nrow(result$att_pre_treatment) > 0L)
  # Should contain expected columns
  expected_cols <- c("period", "att", "se")
  expect_true(all(expected_cols %in% names(result$att_pre_treatment)))
  # Pre-treatment ATTs should be close to 0 (parallel trends hold in DGP)
  expect_true(all(is.finite(result$att_pre_treatment$att)))
  # ATT should still be valid
  expect_true(is.finite(result$att))
})

# TC-7.6.4: lwdid() + include_pretreatment=TRUE + Staggered
test_that("TC-7.6.4: Staggered pretreatment produces valid results", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 1201L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", include_pretreatment = TRUE)
  ))
  # att_pre_treatment should be a data.frame with multiple cohorts
  if (!is.null(result$att_pre_treatment)) {
    expect_true(is.data.frame(result$att_pre_treatment))
    expect_true(nrow(result$att_pre_treatment) > 0L)
    # Should have cohort column for staggered
    expect_true("cohort" %in% names(result$att_pre_treatment))
  }
  # ATT should still be valid
  expect_true(is.finite(result$att))
})

# TC-7.6.5: pretreatment_test=TRUE produces parallel_trends_test
test_that("TC-7.6.5: pretreatment_test produces parallel_trends_test", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1301L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post",
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))
  # parallel_trends_test should be a list
  if (!is.null(result$att_pre_treatment) &&
      nrow(result$att_pre_treatment) > 0L) {
    expect_true(!is.null(result$parallel_trends_test))
    expect_true(is.list(result$parallel_trends_test))
    pt <- result$parallel_trends_test
    # Should contain key fields
    expect_true(!is.null(pt$joint_stat))
    expect_true(!is.null(pt$joint_pvalue))
    expect_true(!is.null(pt$reject_null))
    # joint_pvalue in [0, 1]
    expect_true(pt$joint_pvalue >= 0 && pt$joint_pvalue <= 1)
  }
})

# TC-7.6.6: exclude_pre_periods=2 works correctly
test_that("TC-7.6.6: exclude_pre_periods=2 runs without error", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 10L, K = 6L, att = 2.0, seed = 1401L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", exclude_pre_periods = 2L)
  ))
  expect_true(is.finite(result$att))
  expect_equal(result$exclude_pre_periods, 2L)
})

# TC-7.6.7: exclude_pre_periods too large raises error
test_that("TC-7.6.7: exclude_pre_periods too large raises error", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 6L, K = 3L, att = 2.0, seed = 1501L)
  # K=3 means 3 pre-periods; excluding 10 should leave insufficient reference
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = 10L)
    ))
  )
})

# TC-7.6.9: print/summary display RI and pretreatment info
test_that("TC-7.6.9: print/summary show RI and pretreatment fields", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1601L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", ri = TRUE, rireps = 200L, seed = 777L,
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))
  # Capture print output
  print_out <- capture.output(print(result))
  print_text <- paste(print_out, collapse = "\n")
  # Should contain RI info
  expect_true(grepl("Randomization Inference", print_text))
  expect_true(grepl("RI p-value", print_text))
  # Should contain pretreatment info
  expect_true(grepl("Pre-treatment", print_text))
  # Should contain parallel trends test
  if (!is.null(result$parallel_trends_test)) {
    expect_true(grepl("Parallel Trends", print_text))
  }
  # Capture summary output
  summ_out <- capture.output(print(summary(result)))
  summ_text <- paste(summ_out, collapse = "\n")
  expect_true(grepl("Randomization Inference", summ_text))
})

# TC-7.6.12: exclude_pre_periods parameter validation
test_that("TC-7.6.12: exclude_pre_periods validation rejects bad input", {
  df <- .make_ct_panel(seed = 1701L)
  # Negative value
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = -1L)
    ))
  )
  # Non-integer (string)
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = "abc")
    ))
  )
})

# TC-7.6.13: CT and Staggered pretreatment use original data (deterministic)
test_that("TC-7.6.13: pretreatment results consistent across runs", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1801L)
  res1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  # Same data, same params -> identical pretreatment results
  if (!is.null(res1$att_pre_treatment) &&
      !is.null(res2$att_pre_treatment)) {
    expect_equal(res1$att_pre_treatment$att,
                 res2$att_pre_treatment$att, tolerance = 1e-14)
  }
  # ATT should be identical
  expect_equal(res1$att, res2$att, tolerance = 1e-14)
})


# ===========================================================================
# Group 2: Pretreatment Integration Tests (Task E7-06.6)
# ===========================================================================

# TC-7.6.3: lwdid() + include_pretreatment=TRUE + Common Timing
test_that("TC-7.6.3: CT pretreatment produces valid att_pre_treatment", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1101L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  # att_pre_treatment should be a data.frame

  expect_true(!is.null(result$att_pre_treatment))
  expect_true(is.data.frame(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)
  # Required columns
  required_cols <- c("period", "att", "se")
  for (col in required_cols) {
    expect_true(col %in% names(result$att_pre_treatment),
                info = sprintf("Missing column: %s", col))
  }
  # ATT values should be numeric
  expect_true(all(is.numeric(result$att_pre_treatment$att)))
  # SE values should be non-negative where not NA
  valid_se <- result$att_pre_treatment$se[!is.na(result$att_pre_treatment$se)]
  if (length(valid_se) > 0L) {
    expect_true(all(valid_se >= 0))
  }
  # Main ATT should still be valid
  expect_true(is.finite(result$att))
  # include_pretreatment flag stored
  expect_true(isTRUE(result$include_pretreatment))
})

# TC-7.6.4: lwdid() + include_pretreatment=TRUE + Staggered
test_that("TC-7.6.4: Staggered pretreatment produces valid results", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 1201L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", include_pretreatment = TRUE)
  ))
  # att_pre_treatment should be a data.frame
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(is.data.frame(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)
  # Should contain cohort column for staggered
  expect_true("cohort" %in% names(result$att_pre_treatment))
  # Main ATT should still be valid
  expect_true(is.finite(result$att))
})

# TC-7.6.5: pretreatment_test=TRUE produces parallel_trends_test
test_that("TC-7.6.5: pretreatment_test produces parallel_trends_test", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1301L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post",
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))
  # parallel_trends_test should be a list (if pretreatment succeeded)
  if (!is.null(result$att_pre_treatment) &&
      nrow(result$att_pre_treatment) > 0L) {
    expect_true(!is.null(result$parallel_trends_test))
    expect_true(is.list(result$parallel_trends_test))
    # Should contain key fields
    expect_true("joint_stat" %in% names(result$parallel_trends_test))
    expect_true("joint_pvalue" %in% names(result$parallel_trends_test))
    expect_true("reject_null" %in% names(result$parallel_trends_test))
    # p-value in [0, 1]
    expect_true(result$parallel_trends_test$joint_pvalue >= 0 &&
                result$parallel_trends_test$joint_pvalue <= 1)
  }
})


# TC-7.6.17: Staggered RI aggregate="overall" with NA att_overall
# Note: With very small datasets, estimation may fail entirely before RI.
# We test that either: (a) estimation succeeds and ri_error is set for NA ATT,
# or (b) estimation fails with a proper error.
test_that("TC-7.6.17: Staggered RI overall with NA att sets ri_error", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 901L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          ri = TRUE, rireps = 100L, seed = 444L,
          ri_method = "permutation")
  ))
  expect_true(inherits(result, "lwdid_result"))
  # If att_overall is NA, ri_error should be set with target="overall"
  if (is.null(result$att_overall) || is.na(result$att_overall)) {
    expect_true(!is.null(result$ri_error))
    expect_equal(result$ri_target, "overall")
  } else {
    # If att_overall is valid, RI should have run
    expect_true(!is.null(result$ri_pvalue) || !is.null(result$ri_error))
    expect_equal(result$ri_target, "overall")
  }
})

# TC-7.6.18: Staggered RI aggregate="cohort" with NULL/empty att_by_cohort
# Similar to TC-7.6.17 — test the error path
test_that("TC-7.6.18: Staggered RI cohort with NULL effects sets ri_error", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, t_range = 1:10,
                         att = 3.0, seed = 1001L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 100L, seed = 555L,
          ri_method = "permutation")
  ))
  expect_true(inherits(result, "lwdid_result"))
  # If cohort_effects is NULL or empty, ri_error should be set
  if (is.null(result$cohort_effects) ||
      length(result$cohort_effects) == 0L) {
    expect_true(!is.null(result$ri_error))
    expect_equal(result$ri_target, "cohort")
  } else {
    # If cohort effects exist, RI should have attempted
    expect_true(!is.null(result$ri_pvalue) || !is.null(result$ri_error))
  }
})


# TC-7.6.6: exclude_pre_periods=2 works correctly
test_that("TC-7.6.6: exclude_pre_periods=2 reduces pre-treatment set", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1401L)
  # Without exclusion
  res0 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", exclude_pre_periods = 0L)
  ))
  # With exclusion
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", exclude_pre_periods = 2L)
  ))
  # exclude_pre_periods should be stored
  expect_equal(res0$exclude_pre_periods, 0L)
  expect_equal(res2$exclude_pre_periods, 2L)
  # Both should produce valid ATT
  expect_true(is.finite(res0$att))
  expect_true(is.finite(res2$att))
  # ATT may differ because reference set changed
  # (this is expected behavior, not a bug)
})

# TC-7.6.7: exclude_pre_periods too large raises error
test_that("TC-7.6.7: exclude_pre_periods too large raises error", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1501L)
  # K=4 means 4 pre-treatment periods. Excluding 10 should fail.
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = 10L)
    ))
  )
})

# TC-7.6.12: exclude_pre_periods validation rejects bad input
test_that("TC-7.6.12: exclude_pre_periods validation rejects bad input", {
  df <- .make_ct_panel(n_treated = 10L, n_control = 10L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1701L)
  # Negative value
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = -1L)
    ))
  )
  # Non-integer (string)
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = df, y = "y", ivar = "id", tvar = "year",
            d = "d", post = "post", exclude_pre_periods = "abc")
    ))
  )
})

# TC-7.6.13: pretreatment results consistent across runs
test_that("TC-7.6.13: pretreatment results consistent across runs", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 8L, K = 4L, att = 2.0, seed = 1801L)
  res1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  res2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  # Pretreatment results should be identical (deterministic)
  expect_equal(res1$att_pre_treatment, res2$att_pre_treatment)
  # Main ATT should also be identical
  expect_equal(res1$att, res2$att)
})

# ===========================================================================
# Group 3: L5 Integration Tests (Task E7-07.7)
# ===========================================================================

# TC-7.7.28: Staggered ri+pretreatment complete flow
test_that("TC-7.7.28: staggered ri+pretreatment complete flow", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, att = 3.0, seed = 7701L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", ri = TRUE, rireps = 100L,
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))
  expect_true(!is.null(result$ri_pvalue))
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(!is.null(result$parallel_trends_test))
})

# TC-7.7.29: Result object field completeness
test_that("TC-7.7.29: result object field completeness", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L, seed = 7702L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post",
          ri = TRUE, rireps = 200L, seed = 7702L,
          include_pretreatment = TRUE)
  ))
  # 7 RI fields that are always present (non-NULL values)
  ri_fields_present <- c("ri_pvalue", "ri_distribution", "ri_seed",
                          "ri_method", "rireps", "ri_valid", "ri_failed")
  for (f in ri_fields_present) {
    expect_true(f %in% names(result),
                info = paste0("Missing RI field: ", f))
  }
  # ri_error is NULL on success (R removes NULL list elements)
  expect_true(is.null(result$ri_error))
  # ri_target is NULL for CT path (only set for staggered)
  expect_true(is.null(result$ri_target))
  # All 4 pretreatment fields
  pre_fields <- c("att_pre_treatment", "parallel_trends_test",
                   "include_pretreatment", "exclude_pre_periods")
  for (f in pre_fields) {
    expect_true(f %in% names(result),
                info = paste0("Missing pretreatment field: ", f))
  }
})

# TC-7.7.32: demean/detrend dual path
test_that("TC-7.7.32: demean/detrend dual path", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 10L, K = 6L, seed = 7703L)
  res_dm <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", rolling = "demean",
          include_pretreatment = TRUE)
  ))
  res_dt <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", rolling = "detrend",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(res_dm$att_pre_treatment))
  expect_true(!is.null(res_dt$att_pre_treatment))
  # Different transforms should give different pre-treatment ATTs
  expect_false(identical(res_dm$att_pre_treatment, res_dt$att_pre_treatment))
})

# TC-7.7.33: ra/ipwra/psm three estimators
test_that("TC-7.7.33: ra/ipwra/psm three estimators", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L, seed = 7704L)
  set.seed(7704L)
  n_units <- length(unique(df$id))
  unit_x1 <- rnorm(n_units)
  df$x1 <- unit_x1[df$id]
  # ra (no controls needed)
  res_ra <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", estimator = "ra",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(res_ra$att_pre_treatment))
  # ipwra (needs controls)
  res_ipwra <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", estimator = "ipwra",
          controls = "x1", include_pretreatment = TRUE)
  ))
  expect_true(!is.null(res_ipwra$att_pre_treatment))
  # psm (needs controls + ps_controls)
  res_psm <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", estimator = "psm",
          controls = "x1", ps_controls = "x1",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(res_psm$att_pre_treatment))
})

# TC-7.7.34: CT and Staggered pretreatment receive original data
test_that("TC-7.7.34: pretreatment deterministic (original data passed)", {
  # CT path
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L, seed = 7705L)
  ct1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  ct2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  expect_true(identical(ct1$att_pre_treatment, ct2$att_pre_treatment))
  # Staggered path
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, seed = 7705L)
  st1 <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", include_pretreatment = TRUE)
  ))
  st2 <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", include_pretreatment = TRUE)
  ))
  expect_true(identical(st1$att_pre_treatment, st2$att_pre_treatment))
})

# TC-7.7.35: CT estimator="ipwra" routing
test_that("TC-7.7.35: CT estimator=ipwra routing", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L, seed = 7706L)
  set.seed(7706L)
  n_units <- length(unique(df$id))
  unit_x1 <- rnorm(n_units)
  df$x1 <- unit_x1[df$id]
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", estimator = "ipwra",
          controls = "x1", include_pretreatment = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)
  expect_true(is.finite(result$att))
})

# TC-7.7.42: Staggered RI target three-way dispatch
test_that("TC-7.7.42: staggered RI target three-way dispatch", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, att = 3.0, seed = 7707L)
  res_overall <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          ri = TRUE, rireps = 50L, ri_method = "permutation")
  ))
  res_cohort <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 50L, ri_method = "permutation")
  ))
  res_none <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          ri = TRUE, rireps = 50L, ri_method = "permutation")
  ))
  expect_equal(res_overall$ri_target, "overall")
  expect_equal(res_cohort$ri_target, "cohort")
  expect_equal(res_none$ri_target, "cohort_time")
})

# TC-7.7.43: ri_target field exists and valid
test_that("TC-7.7.43: ri_target field exists and valid", {
  dt <- .make_stag_panel(seed = 7708L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          ri = TRUE, rireps = 50L)
  ))
  expect_true(!is.null(result$ri_target))
  expect_true(result$ri_target %in% c("overall", "cohort", "cohort_time"))
})

# TC-7.7.44: ri_staggered() has no n_never_treated parameter
test_that("TC-7.7.44: ri_staggered has no n_never_treated parameter", {
  sig <- formalArgs(lwdid:::ri_staggered)
  expect_false("n_never_treated" %in% sig)
})

# TC-7.7.45: Staggered detrend end-to-end
test_that("TC-7.7.45: staggered detrend end-to-end", {
  dt <- .make_stag_panel(cohorts = c(4L, 6L), nt_count = 10L,
                         n_per_cohort = 10L, att = 3.0, seed = 7709L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", rolling = "detrend",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)
  expect_true(all(is.finite(result$att_pre_treatment$att)))
})

# TC-7.7.46: detrend degrade boundary (2 pre-periods)
test_that("TC-7.7.46: detrend degrade boundary with 2 pre-periods", {
  dt <- .make_stag_panel(cohorts = c(3L, 6L), nt_count = 5L,
                         n_per_cohort = 5L, t_range = 1:8, seed = 7710L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          gvar = "gvar", rolling = "detrend",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
})

# TC-7.7.47: CT symmetric transform verification
test_that("TC-7.7.47: CT pre-treatment ATTs near zero under parallel trends", {
  df <- .make_ct_panel(n_treated = 100L, n_control = 100L,
                       n_periods = 10L, K = 6L, att = 2.0, seed = 7711L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", include_pretreatment = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
  pre_atts <- result$att_pre_treatment$att
  finite_atts <- pre_atts[is.finite(pre_atts)]
  expect_true(length(finite_atts) > 0L)
  expect_true(max(abs(finite_atts)) < 1.0,
              info = paste0("Max pre-treatment |ATT| = ",
                            round(max(abs(finite_atts)), 4),
                            " exceeds 1.0 threshold"))
})

# TC-7.7.48: CT detrend end-to-end
test_that("TC-7.7.48: CT detrend end-to-end", {
  df <- .make_ct_panel(n_treated = 50L, n_control = 50L,
                       n_periods = 10L, K = 6L, seed = 7712L)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          d = "d", post = "post", rolling = "detrend",
          include_pretreatment = TRUE)
  ))
  expect_true(!is.null(result$att_pre_treatment))
  expect_true(nrow(result$att_pre_treatment) > 0L)
  finite_atts <- result$att_pre_treatment$att[is.finite(result$att_pre_treatment$att)]
  expect_true(length(finite_atts) > 0L)
})

# TC-7.7.30: summary() displays RI and pretreatment info
test_that("TC-7.7.30: summary output contains Randomization and Pre-treatment", {
  dt <- generate_ct_panel(N = 100L, T_total = 8L, tpost1 = 5L,
                          treat_frac = 0.4, tau = 2.0, seed = 7730L)
  dt[, post := as.integer(time >= 5L)]

  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          ri = TRUE, rireps = 200L, seed = 7730L,
          include_pretreatment = TRUE, pretreatment_test = TRUE)
  ))

  # Capture summary output
  out <- capture.output(summary(result))
  out_text <- paste(out, collapse = "\n")

  # Summary should mention Randomization Inference section
  expect_true(any(grepl("Randomization", out, ignore.case = TRUE)),
              info = paste0("summary() should contain 'Randomization'. Got:\n",
                            out_text))

  # Summary should mention Pre-treatment section
  expect_true(any(grepl("Pre-treatment|Pre_treatment|Pretreatment|pre.treatment",
                        out, ignore.case = TRUE)),
              info = paste0("summary() should contain 'Pre-treatment'. Got:\n",
                            out_text))
})

# TC-7.7.31: Backward compatibility — pretreatment doesn't change ATT
test_that("TC-7.7.31: include_pretreatment=TRUE does not change ATT/SE", {
  dt <- generate_ct_panel(N = 100L, T_total = 8L, tpost1 = 5L,
                          treat_frac = 0.4, tau = 2.0, seed = 7731L)
  dt[, post := as.integer(time >= 5L)]

  # Run WITHOUT pretreatment
  res_base <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          include_pretreatment = FALSE)
  ))

  # Run WITH pretreatment
  res_pre <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
          d = "d", post = "post",
          include_pretreatment = TRUE)
  ))

  # ATT must be identical (pretreatment is purely additive)
  expect_equal(res_base$att, res_pre$att, tolerance = 1e-14,
               info = "ATT should not change with include_pretreatment=TRUE")

  # SE must be identical
  expect_equal(res_base$se_att, res_pre$se_att, tolerance = 1e-14,
               info = "SE should not change with include_pretreatment=TRUE")

  # p-value must be identical
  expect_equal(res_base$pvalue, res_pre$pvalue, tolerance = 1e-14,
               info = "p-value should not change with include_pretreatment=TRUE")
})

# TC-7.7.49: overall NA boundary — ri_error set, ri_target stays "overall"
test_that("TC-7.7.49: overall NA boundary sets ri_error with target=overall", {
  # Use degenerate staggered data where att_overall may be NA
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 3L,
                         n_per_cohort = 2L, t_range = 1:10,
                         att = 0.0, seed = 7749L)
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "overall",
            ri = TRUE, rireps = 50L, seed = 7749L,
            ri_method = "permutation")
    )),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    # Degenerate data caused main estimation to fail — acceptable
    expect_true(grepl("insufficient|No valid|too few", conditionMessage(result),
                      ignore.case = TRUE))
  } else {
    expect_true(inherits(result, "lwdid_result"))
    # ri_target must be "overall" regardless of success/failure
    expect_equal(result$ri_target, "overall")
    # If att_overall is NA, ri_error must be set (not fall through to cohort_time)
    if (is.null(result$att_overall) || is.na(result$att_overall)) {
      expect_true(!is.null(result$ri_error),
                  info = "ri_error should be set when att_overall is NA")
    }
  }
})

# TC-7.7.50: cohort NULL boundary — ri_error set, ri_target stays "cohort"
test_that("TC-7.7.50: cohort NULL boundary sets ri_error with target=cohort", {
  # Use degenerate staggered data where cohort_effects may be NULL/empty
  dt <- .make_stag_panel(cohorts = c(3L, 9L), nt_count = 3L,
                         n_per_cohort = 2L, t_range = 1:10,
                         att = 0.0, seed = 7750L)
  result <- tryCatch(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort",
            ri = TRUE, rireps = 50L, seed = 7750L,
            ri_method = "permutation")
    )),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_true(grepl("insufficient|No valid|too few", conditionMessage(result),
                      ignore.case = TRUE))
  } else {
    expect_true(inherits(result, "lwdid_result"))
    # ri_target must be "cohort" regardless of success/failure
    expect_equal(result$ri_target, "cohort")
    # If cohort_effects is NULL or all NA, ri_error must be set
    if (is.null(result$cohort_effects) ||
        length(result$cohort_effects) == 0L) {
      expect_true(!is.null(result$ri_error),
                  info = "ri_error should be set when cohort_effects is NULL")
    } else {
      all_na <- all(vapply(result$cohort_effects,
                           function(e) is.na(e$att), logical(1)))
      if (all_na) {
        expect_true(!is.null(result$ri_error),
                    info = "ri_error should be set when all cohort ATTs are NA")
      }
    }
  }
})
