# ============================================================================
# test-e5-05-integration.R — E5-05 Aggregation Integration Tests
# Covers: E5-05.6 (dispatch), E5-05.7 (result+extract), E5-05.8 (print),
#         E5-05.9 (numerical), E5-05.10 (E2E)
# ============================================================================

# --- Helper: create staggered panel with treatment effect ---
make_e5_panel <- function(cohort_sizes = c(g3 = 3, g5 = 2),
                          n_nt = 2, n_periods = 6,
                          y_base = 1.0, treatment_effect = 2.0) {
  set.seed(42L)
  units <- list()
  uid <- 1L
  for (i in seq_along(cohort_sizes)) {
    g_val <- as.integer(gsub("g", "", names(cohort_sizes)[i]))
    for (j in seq_len(cohort_sizes[i])) {
      y_vals <- y_base * uid + seq_len(n_periods) * 0.5 +
        rnorm(n_periods, sd = 0.3)
      post_mask <- seq_len(n_periods) >= g_val
      y_vals[post_mask] <- y_vals[post_mask] + treatment_effect
      units[[uid]] <- data.table::data.table(
        id = uid, time = seq_len(n_periods),
        Y = y_vals, gvar = g_val
      )
      uid <- uid + 1L
    }
  }
  for (j in seq_len(n_nt)) {
    units[[uid]] <- data.table::data.table(
      id = uid, time = seq_len(n_periods),
      Y = y_base * uid + seq_len(n_periods) * 0.3 +
        rnorm(n_periods, sd = 0.3),
      gvar = 0L
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

# ============================================================================
# E5-05.6: Aggregation dispatch tests (via lwdid() API)
# ============================================================================

test_that("lwdid: aggregate=none returns gr effects only", {
  dt <- make_e5_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  expect_true(is.data.frame(result$cohort_time_effects))
  expect_true(nrow(result$cohort_time_effects) > 0L)
  expect_null(result$cohort_effects)
  expect_null(result$overall_effect)
  expect_null(result$event_time_effects)
})

test_that("lwdid: aggregate=cohort populates cohort_effects", {
  dt <- make_e5_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$cohort_effects))
  expect_true(is.list(result$cohort_effects))
  expect_true(length(result$cohort_effects) > 0L)
})

test_that("lwdid: aggregate=overall populates both cohort and overall", {
  dt <- make_e5_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$cohort_effects))
  expect_false(is.null(result$overall_effect))
  expect_true(is.numeric(result$overall_effect$att))
})

test_that("lwdid: aggregate=event_time populates event_time_effects", {
  dt <- make_e5_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "event_time",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$event_time_effects))
  expect_true(is.list(result$event_time_effects))
  expect_true(length(result$event_time_effects) > 0L)
  expect_null(result$overall_effect)
})

test_that("lwdid: cohort aggregate auto-switches NYT to NT", {
  dt <- make_e5_panel()
  expect_warning(
    result <- suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "not_yet_treated",
            aggregate = "cohort")
    ),
    class = "lwdid_control_group_switch"
  )
  expect_equal(result$control_group_used, "never_treated")
})

test_that("lwdid: overall aggregate auto-switches NYT to NT", {
  dt <- make_e5_panel()
  expect_warning(
    result <- suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "not_yet_treated",
            aggregate = "overall")
    ),
    class = "lwdid_control_group_switch"
  )
  expect_equal(result$control_group_used, "never_treated")
})

test_that("lwdid: all_others + cohort auto-switches to NT", {
  dt <- make_e5_panel()
  expect_warning(
    result <- suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "all_others",
            aggregate = "cohort")
    ),
    class = "lwdid_control_group_switch"
  )
  expect_equal(result$control_group_used, "never_treated")
})

test_that("lwdid: event_time does NOT trigger control switch", {
  dt <- make_e5_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "event_time",
          control_group = "not_yet_treated")
  ))
  expect_equal(result$control_group_used, "not_yet_treated")
})

# ============================================================================
# E5-05.7: new_lwdid_result() + extract_effects() tests
# ============================================================================

test_that("new_lwdid_result: returns lwdid_result class", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "none"
  )
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$method, "staggered")
})

test_that("new_lwdid_result: top-level inference computed", {
  result <- new_lwdid_result(
    att = 0.1, se_att = 0.02, df_inference = 50L,
    method = "staggered", is_staggered = TRUE, aggregate = "cohort"
  )
  expect_equal(result$t_stat, 0.1 / 0.02, tolerance = 1e-8)
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
  expect_true(result$ci_lower < 0.1)
  expect_true(result$ci_upper > 0.1)
})

test_that("new_lwdid_result: NA att gives NA inference", {
  result <- new_lwdid_result(
    att = NA_real_, se_att = 0.02, df_inference = 50L,
    method = "staggered", is_staggered = TRUE, aggregate = "none"
  )
  expect_true(is.na(result$t_stat))
  expect_true(is.na(result$pvalue))
})

test_that("new_lwdid_result: cohort_agg fields stored", {
  result <- new_lwdid_result(
    att = 0.15, se_att = 0.03, df_inference = 40L,
    att_cohort_agg = 0.15, se_cohort_agg = 0.03,
    ci_cohort_agg = c(0.09, 0.21),
    method = "staggered", is_staggered = TRUE, aggregate = "cohort"
  )
  expect_equal(result$att_cohort_agg, 0.15)
  expect_equal(result$se_cohort_agg, 0.03)
  expect_equal(result$ci_lower_cohort_agg, 0.09)
  expect_equal(result$ci_upper_cohort_agg, 0.21)
})

test_that("new_lwdid_result: control_group_used defaults", {
  result <- new_lwdid_result(
    control_group = "never_treated",
    method = "staggered", is_staggered = TRUE, aggregate = "none"
  )
  expect_equal(result$control_group_used, "never_treated")
})

test_that("extract_effects: non-lwdid_result raises error", {
  expect_error(
    extract_effects(list(att = 1)),
    class = "lwdid_invalid_input"
  )
})

test_that("extract_effects: type=NULL auto-infers from aggregate", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "none",
    att_by_cohort_time = data.frame(
      cohort = 3L, ref_period = 3L, att = 1.5, se = 0.3,
      n_treated = 3L, n_control = 2L, stringsAsFactors = FALSE
    )
  )
  result$cohort_time_effects <- result$att_by_cohort_time
  df <- extract_effects(result)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0L)
})

test_that("extract_effects: type=gr returns cohort_time_effects", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "none",
    att_by_cohort_time = data.frame(
      cohort = c(3L, 5L), ref_period = c(3L, 5L),
      att = c(1.5, 2.0), se = c(0.3, 0.4),
      n_treated = c(3L, 2L), n_control = c(2L, 2L),
      stringsAsFactors = FALSE
    )
  )
  result$cohort_time_effects <- result$att_by_cohort_time
  df <- extract_effects(result, type = "gr")
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2L)
})

test_that("extract_effects: type=cohort converts list to df", {
  ce <- list(
    list(cohort = 3L, att = 1.5, se = 0.3, ci_lower = 0.9,
         ci_upper = 2.1, t_stat = 5.0, pvalue = 0.001,
         n_periods = 4L, n_units = 3L, n_control = 2L,
         df_resid = 3L, df_inference = 3L),
    list(cohort = 5L, att = 2.0, se = 0.4, ci_lower = 1.2,
         ci_upper = 2.8, t_stat = 5.0, pvalue = 0.001,
         n_periods = 2L, n_units = 2L, n_control = 2L,
         df_resid = 2L, df_inference = 2L)
  )
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "cohort",
    cohort_effects = ce
  )
  df <- extract_effects(result, type = "cohort")
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2L)
  expect_true("cohort" %in% names(df))
  expect_true("att" %in% names(df))
  expect_true("se" %in% names(df))
})

test_that("extract_effects: type=overall returns single-row df", {
  oe <- list(att = 1.8, se = 0.25, ci_lower = 1.3, ci_upper = 2.3,
             t_stat = 7.2, pvalue = 0.0001,
             n_treated = 5L, n_control = 2L,
             df_resid = 5L, df_inference = 5L)
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "overall",
    att_overall = oe$att, se_overall = oe$se
  )
  result$overall_effect <- oe
  df <- extract_effects(result, type = "overall")
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1L)
  expect_equal(df$att, 1.8, tolerance = 1e-8)
})

test_that("extract_effects: type=event_time converts list to df", {
  ete <- list(
    list(event_time = 0L, att = 1.5, se = 0.3, ci_lower = 0.9,
         ci_upper = 2.1, t_stat = 5.0, pvalue = 0.001,
         df_inference = 10L, n_cohorts = 2L, weight_sum = 1.0),
    list(event_time = 1L, att = 2.0, se = 0.4, ci_lower = 1.2,
         ci_upper = 2.8, t_stat = 5.0, pvalue = 0.001,
         df_inference = 10L, n_cohorts = 2L, weight_sum = 1.0)
  )
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "event_time",
    event_time_effects = ete
  )
  df <- extract_effects(result, type = "event_time")
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2L)
  expect_true("event_time" %in% names(df))
})

test_that("extract_effects: missing effects warns + empty df", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "cohort",
    cohort_effects = NULL
  )
  expect_warning(
    df <- extract_effects(result, type = "cohort"),
    class = "lwdid_data"
  )
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 0L)
})

test_that("extract_effects: invalid type raises error", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "none"
  )
  expect_error(
    extract_effects(result, type = "invalid"),
    class = "lwdid_invalid_input"
  )
})

# ============================================================================
# E5-05.8: print.lwdid_result() tests
# ============================================================================

test_that("print: displays header and metadata", {
  result <- new_lwdid_result(
    att = 1.5, se_att = 0.3, df_inference = 20L,
    method = "staggered", is_staggered = TRUE, aggregate = "none",
    rolling = "demean", estimator = "ra", depvar = "Y",
    control_group = "never_treated",
    control_group_used = "never_treated"
  )
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Local Wald DID", out_text))
  expect_true(grepl("staggered", out_text))
  expect_true(grepl("demean", out_text, ignore.case = TRUE))
  expect_true(grepl("none", out_text))
})

test_that("print: displays top-level ATT summary", {
  result <- new_lwdid_result(
    att = 2.345, se_att = 0.5, df_inference = 30L,
    method = "staggered", is_staggered = TRUE, aggregate = "none",
    rolling = "demean", estimator = "ra", depvar = "Y"
  )
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("ATT", out_text))
  expect_true(grepl("2\\.345", out_text))
})

test_that("print: returns invisible(x)", {
  result <- new_lwdid_result(
    method = "staggered", is_staggered = TRUE, aggregate = "none",
    rolling = "demean", estimator = "ra", depvar = "Y"
  )
  vis <- withVisible(capture.output(ret <- print(result)))
  expect_identical(ret, result)
})

test_that("print: shows control group switch info", {
  result <- new_lwdid_result(
    att = 1.0, se_att = 0.2, df_inference = 20L,
    method = "staggered", is_staggered = TRUE, aggregate = "cohort",
    rolling = "demean", estimator = "ra", depvar = "Y",
    control_group = "not_yet_treated",
    control_group_used = "never_treated"
  )
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("auto-switched", out_text, ignore.case = TRUE))
})

test_that("print: cohort aggregate shows cohort table", {
  ce <- list(
    list(cohort = 3L, att = 1.5, se = 0.3, ci_lower = 0.9,
         ci_upper = 2.1, t_stat = 5.0, pvalue = 0.001,
         n_periods = 4L, n_units = 3L, n_control = 2L,
         df_resid = 3L, df_inference = 3L)
  )
  result <- new_lwdid_result(
    att = 1.5, se_att = 0.3, df_inference = 20L,
    method = "staggered", is_staggered = TRUE, aggregate = "cohort",
    rolling = "demean", estimator = "ra", depvar = "Y",
    cohort_effects = ce
  )
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Cohort Effects", out_text))
})

test_that("print: overall aggregate shows overall + hint", {
  oe <- list(att = 1.8, se = 0.25, ci_lower = 1.3, ci_upper = 2.3,
             t_stat = 7.2, pvalue = 0.0001,
             n_treated = 5L, n_control = 2L,
             df_resid = 5L, df_inference = 5L)
  ce <- list(
    list(cohort = 3L, att = 1.5, se = 0.3, ci_lower = 0.9,
         ci_upper = 2.1, t_stat = 5.0, pvalue = 0.001,
         n_periods = 4L, n_units = 3L, n_control = 2L,
         df_resid = 3L, df_inference = 3L)
  )
  result <- new_lwdid_result(
    att = 1.8, se_att = 0.25, df_inference = 5L,
    method = "staggered", is_staggered = TRUE, aggregate = "overall",
    rolling = "demean", estimator = "ra", depvar = "Y",
    att_overall = oe$att, se_overall = oe$se,
    cohort_effects = ce
  )
  result$overall_effect <- oe
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Overall Effect", out_text))
  expect_true(grepl("extract_effects", out_text))
})

test_that("print: event_time shows event time table", {
  ete <- list(
    list(event_time = 0L, att = 1.5, se = 0.3, ci_lower = 0.9,
         ci_upper = 2.1, t_stat = 5.0, pvalue = 0.001,
         df_inference = 10L, n_cohorts = 2L, weight_sum = 1.0)
  )
  result <- new_lwdid_result(
    att = 1.5, se_att = 0.3, df_inference = 10L,
    method = "staggered", is_staggered = TRUE, aggregate = "event_time",
    rolling = "demean", estimator = "ra", depvar = "Y",
    event_time_effects = ete
  )
  out <- capture.output(print(result))
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Event-Time Effects", out_text))
})

# ============================================================================
# E5-05.9: Numerical verification
# ============================================================================

test_that("numerical: cohort_agg weighted average ATT", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  # Extract cohort effects and manually compute weighted average
  ce <- result$cohort_effects
  skip_if(is.null(ce) || length(ce) == 0L,
          "No cohort effects produced")
  atts <- vapply(ce, `[[`, numeric(1), "att")
  n_units <- vapply(ce, `[[`, integer(1), "n_units")
  weights <- n_units / sum(n_units)
  manual_agg <- sum(weights * atts)
  expect_equal(result$att_cohort_agg, manual_agg, tolerance = 1e-8)
})

test_that("numerical: cohort_agg SE delta method", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  ce <- result$cohort_effects
  skip_if(is.null(ce) || length(ce) == 0L,
          "No cohort effects produced")
  ses <- vapply(ce, `[[`, numeric(1), "se")
  n_units <- vapply(ce, `[[`, integer(1), "n_units")
  weights <- n_units / sum(n_units)
  manual_se <- sqrt(sum(weights^2 * ses^2))
  expect_equal(result$se_cohort_agg, manual_se, tolerance = 1e-8)
})

test_that("numerical: cohort_agg t_stat and pvalue consistent", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$att_cohort_agg) || is.na(result$att_cohort_agg),
          "No cohort aggregated ATT")
  skip_if(is.null(result$se_cohort_agg) || is.na(result$se_cohort_agg) ||
          result$se_cohort_agg == 0,
          "No valid cohort aggregated SE")
  # t = att / se
  expected_t <- result$att_cohort_agg / result$se_cohort_agg
  expect_equal(result$t_stat, expected_t, tolerance = 1e-8)
  # pvalue in [0, 1]
  expect_true(result$pvalue >= 0 && result$pvalue <= 1)
  # CI contains att
  expect_true(result$ci_lower <= result$att + 1e-8)
  expect_true(result$ci_upper >= result$att - 1e-8)
})

test_that("numerical: top-level ATT overall uses overall_effect", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$overall_effect),
          "No overall effect produced")
  expect_equal(result$att, result$overall_effect$att, tolerance = 1e-8)
})

test_that("numerical: top-level ATT cohort uses cohort_agg", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$att_cohort_agg) || is.na(result$att_cohort_agg),
          "No cohort aggregated ATT")
  expect_equal(result$att, result$att_cohort_agg, tolerance = 1e-8)
})

test_that("numerical: top-level ATT none uses gr weighted avg", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  skip_if(is.null(result$cohort_time_effects) ||
          !is.data.frame(result$cohort_time_effects) ||
          nrow(result$cohort_time_effects) == 0L,
          "No (g,r) effects produced")
  gr <- result$cohort_time_effects
  valid <- gr[!is.na(gr$att), , drop = FALSE]
  skip_if(nrow(valid) == 0L, "No valid (g,r) effects")
  w <- valid$n_treated
  manual_att <- sum(w * valid$att) / sum(w)
  expect_equal(result$att, manual_att, tolerance = 1e-8)
})

# ============================================================================
# E5-05.10: E2E integration
# ============================================================================

test_that("E2E: aggregate=none full pipeline", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$aggregate, "none")

  # extract_effects works
  df <- extract_effects(result)
  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0L)

  # print works without error
  out <- capture.output(print(result))
  expect_true(length(out) > 0L)
})

test_that("E2E: aggregate=cohort full pipeline", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$aggregate, "cohort")

  # extract cohort effects
  df <- extract_effects(result, type = "cohort")
  expect_true(is.data.frame(df))
  skip_if(nrow(df) == 0L, "No cohort effects to extract")
  expect_true("cohort" %in% names(df))
  expect_true("att" %in% names(df))

  # print works
  out <- capture.output(print(result))
  expect_true(length(out) > 0L)
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Cohort", out_text))
})

test_that("E2E: aggregate=overall full pipeline", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          control_group = "never_treated")
  ))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$aggregate, "overall")

  # extract overall effect
  df <- extract_effects(result, type = "overall")
  expect_true(is.data.frame(df))
  skip_if(nrow(df) == 0L, "No overall effect to extract")
  expect_true("att" %in% names(df))

  # print works
  out <- capture.output(print(result))
  expect_true(length(out) > 0L)
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Overall", out_text))
})

test_that("E2E: aggregate=event_time full pipeline", {
  dt <- make_e5_panel(
    cohort_sizes = c(g4 = 5, g7 = 4), n_nt = 5, n_periods = 10
  )
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "event_time",
          control_group = "never_treated")
  ))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$aggregate, "event_time")

  # extract event-time effects
  df <- suppressWarnings(extract_effects(result, type = "event_time"))
  expect_true(is.data.frame(df))
  skip_if(nrow(df) == 0L, "No event-time effects to extract")
  expect_true("event_time" %in% names(df))
  expect_true("att" %in% names(df))

  # print works
  out <- capture.output(print(result))
  expect_true(length(out) > 0L)
  out_text <- paste(out, collapse = "\n")
  expect_true(grepl("Event-Time", out_text))
})
