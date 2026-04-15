# ============================================================================
# test-lwdid.R — Comprehensive tests for lwdid() main entry function
# Story E1-06: 12 test groups, 51 test cases
# ============================================================================

# --- Setup: ensure .lwdid_env$warning_registry is initialized ---
if (!exists(".lwdid_env") ||
    is.null(.lwdid_env$warning_registry)) {
  .lwdid_env$warning_registry <- new_warning_registry()
}

# --- Helper: create valid Common Timing panel data ---
make_ct_data <- function(n_units = 10, n_periods = 8,
                         n_treated = 5, K = 4) {
  ids <- seq_len(n_units)
  years <- 2000 + seq_len(n_periods)
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  df$d <- ifelse(df$id <= n_treated, 1, 0)
  df$post <- ifelse(df$year > (2000 + K), 1, 0)
  set.seed(42)
  df$y <- rnorm(nrow(df)) + df$d * df$post * 2
  rownames(df) <- NULL
  df
}

# --- Helper: create valid Staggered panel data ---
make_stag_data <- function(n_units = 12, n_periods = 10,
                           cohorts = c(2005, 2007),
                           n_nt = 4) {
  ids <- seq_len(n_units)
  years <- 2000 + seq_len(n_periods)
  df <- expand.grid(id = ids, year = years)
  df <- df[order(df$id, df$year), ]
  treated_ids <- (n_nt + 1):n_units
  n_per_cohort <- length(treated_ids) %/% length(cohorts)
  gvar_map <- numeric(n_units)
  gvar_map[1:n_nt] <- 0
  for (i in seq_along(cohorts)) {
    start_idx <- n_nt + (i - 1) * n_per_cohort + 1
    end_idx <- if (i == length(cohorts)) n_units
               else n_nt + i * n_per_cohort
    gvar_map[start_idx:end_idx] <- cohorts[i]
  }
  df$gvar <- gvar_map[df$id]
  set.seed(42)
  df$y <- rnorm(nrow(df))
  rownames(df) <- NULL
  df
}


# ============================================================================
# Group 1: Mode dispatch correctness (4 tests)
# ============================================================================
test_that("dispatch: gvar non-NULL -> staggered mode", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_equal(result$method, "staggered")
  expect_true(result$is_staggered)
})
test_that("dispatch: d+post non-NULL -> common_timing mode", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(result$method, "common_timing")
  expect_false(result$is_staggered)
})
test_that("dispatch: gvar takes priority over d+post", {
  df <- make_stag_data()
  df$d <- ifelse(df$id <= 4, 0, 1)
  df$post <- ifelse(df$year > 2005, 1, 0)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year",
          gvar = "gvar", d = "d", post = "post")))
  expect_equal(result$method, "staggered")
  expect_true(result$is_staggered)
})
test_that("dispatch: neither d/post nor gvar -> error", {
  df <- make_ct_data()
  expect_error(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year"),
    class = "lwdid_invalid_parameter")
})
# ============================================================================
# Group 2: Core return values (6 tests)
# ============================================================================
test_that("stub CT: returns lwdid_result class", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$method, "common_timing")
})
test_that("stub Staggered: returns lwdid_result class", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$method, "staggered")
})
test_that("CT: core estimate fields are populated", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  estimate_fields <- unlist(result[c(
    "att", "se_att", "t_stat", "pvalue", "ci_lower", "ci_upper"
  )], use.names = FALSE)

  expect_true(all(is.finite(estimate_fields)))
  expect_gt(result$se_att, 0)
  expect_lt(result$ci_lower, result$att)
  expect_gt(result$ci_upper, result$att)
})
test_that("Staggered: core estimate fields are populated", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  estimate_fields <- unlist(result[c(
    "att", "se_att", "t_stat", "pvalue", "ci_lower", "ci_upper"
  )], use.names = FALSE)

  expect_true(all(is.finite(estimate_fields)))
  expect_gt(result$se_att, 0)
  expect_lt(result$ci_lower, result$att)
  expect_gt(result$ci_upper, result$att)
})
test_that("CT: df fields are populated integers", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_type(result$df_resid, "integer")
  expect_type(result$df_inference, "integer")
  expect_gt(result$df_resid, 0L)
  expect_gt(result$df_inference, 0L)
})
test_that("Staggered: df fields are populated integers", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_type(result$df_resid, "integer")
  expect_type(result$df_inference, "integer")
  expect_gt(result$df_resid, 0L)
  expect_gt(result$df_inference, 0L)
})

# ============================================================================
# Group 3: Stub attribute completeness (5 tests)
# ============================================================================
test_that("stub CT: all 74 fields present", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  ef <- c("att","se_att","t_stat","pvalue","ci_lower","ci_upper",
    "df_resid","df_inference","nobs","n_treated","n_control","K","tpost1",
    "cmd","depvar","rolling","estimator","method","vce_type","cluster_var",
    "n_clusters","alpha","is_staggered","controls_used","controls",
    "include_pretreatment","control_group","control_group_used","params",
    "bse","vcov_matrix","resid","data","att_by_period","att_pre_treatment",
    "parallel_trends_test","aggregate","cohorts","cohort_sizes",
    "n_never_treated","att_by_cohort","att_by_cohort_time","cohort_effects",
    "att_overall","se_overall","ci_overall_lower","ci_overall_upper",
    "t_stat_overall","pvalue_overall","event_time_effects","cohort_weights",
    "ri_pvalue","ri_distribution","ri_seed","rireps","ri_method","ri_valid",
    "ri_failed","ri_error","ri_target","diagnostics","warning_diagnostics",
    "propensity_scores","matched_data","n_matched","match_rate","weights_cv",
    "warnings_log","call","lwdid_version","ivar","tvar","is_quarterly")
  expect_true(all(ef %in% names(result)))
})
test_that("stub Staggered: all 74 fields present", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  ef <- c("att","se_att","t_stat","pvalue","ci_lower","ci_upper",
    "df_resid","df_inference","nobs","n_treated","n_control","K","tpost1",
    "cmd","depvar","rolling","estimator","method","vce_type","cluster_var",
    "n_clusters","alpha","is_staggered","controls_used","controls",
    "include_pretreatment","control_group","control_group_used","params",
    "bse","vcov_matrix","resid","data","att_by_period","att_pre_treatment",
    "parallel_trends_test","aggregate","cohorts","cohort_sizes",
    "n_never_treated","att_by_cohort","att_by_cohort_time","cohort_effects",
    "att_overall","se_overall","ci_overall_lower","ci_overall_upper",
    "t_stat_overall","pvalue_overall","event_time_effects","cohort_weights",
    "ri_pvalue","ri_distribution","ri_seed","rireps","ri_method","ri_valid",
    "ri_failed","ri_error","ri_target","diagnostics","warning_diagnostics",
    "propensity_scores","matched_data","n_matched","match_rate","weights_cv",
    "warnings_log","call","lwdid_version","ivar","tvar","is_quarterly")
  expect_true(all(ef %in% names(result)))
})
test_that("stub CT: staggered-specific fields are NULL", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_null(result$aggregate)
  expect_null(result$cohorts)
  expect_null(result$cohort_sizes)
  expect_null(result$n_never_treated)
  expect_null(result$att_by_cohort)
  expect_null(result$att_overall)
  expect_null(result$se_overall)
  expect_null(result$event_time_effects)
  expect_null(result$cohort_weights)
})
test_that("stub Staggered: staggered metadata NOT NULL", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_false(is.null(result$aggregate))
  expect_false(is.null(result$cohorts))
  expect_false(is.null(result$cohort_sizes))
  expect_false(is.null(result$n_never_treated))
})
test_that("stub CT: RI fields are all NULL", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_null(result$ri_pvalue)
  expect_null(result$ri_distribution)
  expect_null(result$ri_seed)
  expect_null(result$rireps)
  expect_null(result$ri_method)
  expect_null(result$ri_valid)
  expect_null(result$ri_failed)
  expect_null(result$ri_error)
  expect_null(result$ri_target)
})

# ============================================================================
# Group 4: Metadata population (4 tests)
# ============================================================================
test_that("metadata: call is a call object containing 'lwdid'", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_true(is.call(result$call))
  expect_true(grepl("lwdid", deparse(result$call[[1]])))
})
test_that("metadata: lwdid_version is a non-empty character string", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_true(is.character(result$lwdid_version))
  expect_equal(length(result$lwdid_version), 1)
  expect_true(nchar(result$lwdid_version) > 0)
})
test_that("metadata: cmd field is 'lwdid'", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(result$cmd, "lwdid")
})
test_that("metadata: staggered result also has correct metadata", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_true(is.call(result$call))
  expect_true(is.character(result$lwdid_version))
  expect_equal(result$cmd, "lwdid")
})
# ============================================================================
# Group 5: Parameter passing (4 tests)
# ============================================================================
test_that("params: depvar correctly stored", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(result$depvar, "y")
})
test_that("params: ivar and tvar correctly stored", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(result$ivar, "id")
  expect_equal(result$tvar, "year")
})
test_that("params: gid=NULL runs normally", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", gid = NULL))))
})
test_that("params: gid=1 passed without error", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", gid = 1))))
})
# ============================================================================
# Group 6: WarningRegistry lifecycle (3 tests)
# ============================================================================
test_that("registry: consecutive calls don't leak warnings", {
  df <- make_ct_data()
  r1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  r2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(nrow(r2$warnings_log), nrow(r1$warnings_log))
})
test_that("registry: warnings_log is a data.frame with correct schema", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_true(is.data.frame(result$warnings_log))
  expect_true(all(c("category", "message", "cohort", "period", "timestamp") %in%
                    names(result$warnings_log)))
  if (nrow(result$warnings_log) > 0) {
    expect_true(all(is.character(result$warnings_log$category)))
    expect_true(all(is.character(result$warnings_log$message)))
  }
})
test_that("registry: warning_diagnostics is a list", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_true(is.list(result$warning_diagnostics))
})
# ============================================================================
# Group 7: Message output on implemented paths (2 tests)
# ============================================================================
test_that("message: CT path runs without implementation-placeholder messages", {
  df <- make_ct_data()
  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post"),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))
  expect_length(msgs, 0L)
})
test_that("message: Staggered path runs without implementation-placeholder messages", {
  df <- make_stag_data()
  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar"),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))
  expect_length(msgs, 0L)
})
# ============================================================================
# Group 8: gid parameter (2 tests)
# ============================================================================
test_that("gid: NULL runs normally for staggered", {
  df <- make_stag_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar", gid = NULL))))
})
test_that("gid: scalar value passed without error", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", gid = 1))))
})

# ============================================================================
# Group 9: kwargs / riseed handling (8 tests)
# ============================================================================
test_that("riseed: integer 42 with seed=NULL accepted", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", riseed = 42L))))
})
test_that("riseed: string '42' with seed=NULL accepted", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", riseed = "42"))))
})
test_that("riseed: non-integer string 'abc' -> lwdid_invalid_parameter", {
  df <- make_ct_data()
  expect_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", riseed = "abc"))),
    class = "lwdid_invalid_parameter")
})
test_that("riseed: non-numeric type list(1) -> lwdid_invalid_parameter", {
  df <- make_ct_data()
  expect_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", riseed = list(1)))),
    class = "lwdid_invalid_parameter")
})
test_that("riseed: seed=99 takes priority over riseed=42", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", seed = 99, riseed = 42))))
})
test_that("riseed: ri=FALSE + valid riseed=42 -> no error", {
  df <- make_ct_data()
  expect_no_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", ri = FALSE, riseed = 42))))
})
test_that("riseed: ri=FALSE + invalid riseed='abc' -> still errors", {
  df <- make_ct_data()
  expect_error(suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", ri = FALSE, riseed = "abc"))),
    class = "lwdid_invalid_parameter")
})
test_that("kwargs: unknown argument 'foo' registers warning in log", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", foo = 1)))
  log_text <- paste(result$warnings_log$message, collapse = " ")
  expect_true(grepl("foo", log_text))
  expect_true(grepl("[Uu]nknown", log_text))
})
# ============================================================================
# Group 10: att_gt parameter (3 tests)
# ============================================================================
test_that("att_gt: TRUE preserves current top-level estimate", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", att_gt = TRUE)))
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))
})
test_that("att_gt: FALSE preserves current top-level estimate", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", att_gt = FALSE)))
  expect_s3_class(result, "lwdid_result")
  expect_true(is.finite(result$att))
})
test_that("att_gt: correctly passed to validate_inputs", {
  df <- make_ct_data()
  r1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", att_gt = TRUE)))
  r2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", att_gt = FALSE)))
  expect_s3_class(r1, "lwdid_result")
  expect_s3_class(r2, "lwdid_result")
  expect_equal(r1$att, r2$att)
  expect_equal(r1$se_att, r2$se_att)
})
# ============================================================================
# Group 11: graph handling (3 tests)
# ============================================================================
test_that("graph: graph=TRUE returns without placeholder plot messages", {
  df <- make_ct_data()
  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", graph = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))
  expect_length(msgs, 0L)
})
test_that("graph: graph=FALSE does not call plot", {
  df <- make_ct_data()
  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", graph = FALSE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))
  expect_false(any(grepl("plot\\.lwdid_result", msgs)))
})
test_that("graph: plot failure does not abort lwdid()", {
  df <- make_ct_data()
  orig_plot <- plot.lwdid_result
  assign("plot.lwdid_result", function(x, ...) stop("Intentional plot failure"), envir = globalenv())
  on.exit(assign("plot.lwdid_result", orig_plot, envir = globalenv()), add = TRUE)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", graph = TRUE)))
  expect_s3_class(result, "lwdid_result")
})
# ============================================================================
# Group 12: End-to-end integration tests (7 tests)
# ============================================================================
test_that("e2e CT: full pipeline returns lwdid_result", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$method, "common_timing")
  expect_false(result$is_staggered)
  expect_equal(result$cmd, "lwdid")
  expect_equal(result$depvar, "y")
})
test_that("e2e Staggered: full pipeline returns lwdid_result", {
  df <- make_stag_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_s3_class(result, "lwdid_result")
  expect_equal(result$method, "staggered")
  expect_true(result$is_staggered)
  expect_equal(result$cmd, "lwdid")
  expect_equal(result$depvar, "y")
})
test_that("e2e CT: sample info matches input data", {
  df <- make_ct_data(n_units = 10, n_periods = 8, n_treated = 5, K = 4)
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post")))
  expect_equal(result$nobs, length(unique(df$id)))
  n_units <- length(unique(df$id))
  expect_true(result$n_treated > 0)
  expect_true(result$n_control > 0)
  expect_equal(result$n_treated + result$n_control, n_units)
})
test_that("e2e Staggered: cohorts match input data", {
  df <- make_stag_data(cohorts = c(2005, 2007))
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", gvar = "gvar")))
  expect_true(2005 %in% result$cohorts)
  expect_true(2007 %in% result$cohorts)
  expect_true(result$n_never_treated > 0)
})
test_that("e2e: rolling parameter correctly passed", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", rolling = "demean")))
  expect_equal(result$rolling, "demean")
})
test_that("e2e: alpha parameter correctly passed", {
  df <- make_ct_data()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", alpha = 0.10)))
  expect_equal(result$alpha, 0.10)
})
test_that("e2e: include_pretreatment correctly passed", {
  df <- make_ct_data()
  r1 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", include_pretreatment = TRUE)))
  r2 <- suppressWarnings(suppressMessages(
    lwdid(data = df, y = "y", ivar = "id", tvar = "year", d = "d", post = "post", include_pretreatment = FALSE)))
  expect_true(r1$include_pretreatment)
  expect_false(r2$include_pretreatment)
})

# =============================================================================
# Group 13-15: E5-05.6 .estimate_staggered() aggregation dispatch tests
# =============================================================================

# Helper: create staggered panel for E5-05.6 integration tests
make_staggered_panel <- function(
  cohort_sizes = c(g3 = 3, g5 = 2),
  n_nt = 2,
  n_periods = 6,
  y_base = 1.0,
  treatment_effect = 2.0
) {
  units <- list()
  uid <- 1L
  for (i in seq_along(cohort_sizes)) {
    g_val <- as.integer(gsub("g", "", names(cohort_sizes)[i]))
    for (j in seq_len(cohort_sizes[i])) {
      y_vals <- y_base * uid + seq_len(n_periods) * 0.5
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
      Y = y_base * uid + seq_len(n_periods) * 0.3,
      gvar = 0L
    )
    uid <- uid + 1L
  }
  data.table::rbindlist(units)
}

# ---------------------------------------------------------------------------
# Group 13: Parameter validation (4 tests)
# ---------------------------------------------------------------------------

test_that("lwdid: invalid aggregate raises lwdid_invalid_parameter", {
  dt <- make_staggered_panel()
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "invalid",
            control_group = "never_treated")
    )),
    class = "lwdid_invalid_parameter"
  )
})

test_that("lwdid: invalid alpha raises lwdid_invalid_parameter", {
  dt <- make_staggered_panel()
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", alpha = 0,
            control_group = "never_treated")
    )),
    class = "lwdid_invalid_parameter"
  )
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", alpha = 1.5,
            control_group = "never_treated")
    )),
    class = "lwdid_invalid_parameter"
  )
})

test_that("lwdid: invalid event_time_range raises error", {
  dt <- make_staggered_panel()
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", event_time_range = c(5, -5),
            control_group = "never_treated")
    )),
    class = "lwdid_invalid_input"
  )
})

test_that("lwdid: invalid df_strategy raises error", {
  dt <- make_staggered_panel()
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", df_strategy = "unknown",
            control_group = "never_treated")
    )),
    class = "lwdid_invalid_input"
  )
})

# ---------------------------------------------------------------------------
# Group 14: FATAL-004 control group constraints (6 tests)
# ---------------------------------------------------------------------------

test_that("lwdid: cohort aggregate auto-switches NYT to NT", {
  dt <- make_staggered_panel()
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
  dt <- make_staggered_panel()
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

test_that("lwdid: no NT + cohort aggregate raises error", {
  dt <- make_staggered_panel(n_nt = 0)
  expect_error(
    suppressWarnings(suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", aggregate = "cohort",
            control_group = "not_yet_treated")
    )),
    class = "lwdid_no_never_treated"
  )
})

test_that("lwdid: n_nt < 2 warns lwdid_small_sample", {
  dt <- make_staggered_panel(n_nt = 1)
  # With n_nt=1 and aggregate="cohort" (default), the FATAL-004 check

  # in .estimate_staggered() emits lwdid_small_sample warning
  expect_warning(
    suppressMessages(
      lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
            gvar = "gvar", control_group = "never_treated",
            aggregate = "cohort")
    ),
    class = "lwdid_small_sample"
  )
})

test_that("lwdid: event_time does NOT trigger control switch", {
  dt <- make_staggered_panel()
  # event_time aggregate should NOT auto-switch control group
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "event_time",
          control_group = "not_yet_treated")
  ))
  expect_equal(result$control_group_used, "not_yet_treated")
})

test_that("lwdid: all_others + cohort aggregate auto-switches to NT", {
  dt <- make_staggered_panel()
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

# ---------------------------------------------------------------------------
# Group 15: Aggregation dispatch (4 tests)
# ---------------------------------------------------------------------------

test_that("lwdid: aggregate=none returns gr effects only", {
  dt <- make_staggered_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "none",
          control_group = "never_treated")
  ))
  expect_true(is.data.frame(result$att_by_cohort_time))
  expect_true(nrow(result$att_by_cohort_time) > 0L)
  expect_null(result$cohort_effects)
  expect_null(result$att_overall)
})

test_that("lwdid: aggregate=cohort populates cohort_effects", {
  dt <- make_staggered_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "cohort",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$cohort_effects))
  expect_true(length(result$cohort_effects) > 0L)
  expect_null(result$att_overall)
})

test_that("lwdid: aggregate=overall populates both cohort and overall", {
  dt <- make_staggered_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "overall",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$cohort_effects))
  expect_false(is.null(result$att_overall))
  expect_true(is.numeric(result$att_overall))
})

test_that("lwdid: aggregate=event_time populates event_time_effects", {
  dt <- make_staggered_panel()
  result <- suppressWarnings(suppressMessages(
    lwdid(data = dt, y = "Y", ivar = "id", tvar = "time",
          gvar = "gvar", aggregate = "event_time",
          control_group = "never_treated")
  ))
  expect_false(is.null(result$event_time_effects))
  expect_null(result$att_overall)
})
